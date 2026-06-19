program fesom_advhordump
    ! M1.1 horizontal-advection operator byte-gate driver. Loads pi 1-rank with the
    ! SAME rotation as the FESOM2 pi run (alpha/beta/gamma=50/15/-90, cyclic 360),
    ! computes geometry, builds helem (linfs full-cell: zbar(nz)-zbar(nz+1)),
    ! PRESCRIBES an analytic velocity + tracer field (identical formula to the
    ! FESOM2 oracle src/fesom_advhor_dump.F90, computed from the byte-identical
    ! coordinates), runs the full horizontal-advection pipeline, and dumps every
    ! field in the FESOM2 oracle's format. tools/advhor_diff.py then compares for
    ! max|delta|=0.
    !
    ! Pipeline (per FESOM2 do_oce_adv_tra horizontal path):
    !   muscl_adv_init -> nboundary_lay, edge_up_dn_tri
    !   tracer_gradient_elements -> tr_xy ; fill_up_dn_grad -> edge_up_dn_grad
    !   adv_tra_hor_upw1/_muscl -> adv_flux_hor (edge) ; scatter -> del_ttf_advhoriz
    !
    !   FESOM3_MESH_DIR    mesh dir   (default: pi)
    !   FESOM3_ADVHOR_OUT  out path   (default: advhor_f3.bin)
    use mpi
    use mod_precision,    only: WP, MP
    use mod_mesh,         only: t_mesh
    use mod_partit,       only: t_partit
    use mod_partitioning, only: par_init, par_ex
    use mod_mesh_read,    only: read_mesh
    use mod_mesh_areas,   only: compute_geometry
    use mod_tracer,       only: t_tracer_work
    use oce_tracer_grad,  only: tracer_gradient_elements
    use oce_muscl_adv,    only: muscl_adv_init, fill_up_dn_grad
    use oce_adv_tra_hor,  only: adv_tra_hor_upw1, adv_tra_hor_muscl
    use oce_adv_tra_flux, only: oce_tra_adv_flux2dtracer
    use mod_advhor_dump,  only: advhor_dump_open, advhor_dump_close, wr_r2, wr_r3, wr_i1, wr_i2
    implicit none

    ! Controlled-input parameters (must equal the FESOM2 oracle's)
    real(kind=WP), parameter :: dt      = 1800.0_WP   ! time step [s]
    real(kind=WP), parameter :: num_ord = 0.75_WP     ! MUSCL 4th-order fraction (exercises both terms)

    character(len=512) :: mesh_dir, out_path
    type(t_partit)      :: partit
    type(t_mesh)        :: mesh
    type(t_tracer_work) :: twork
    integer :: nsw, e, n, nz, nl, u
    real(kind=WP), allocatable :: vel(:,:,:), ttf(:,:), tr_xy(:,:,:), eudg(:,:,:)
    real(kind=WP), allocatable :: aflux_u(:,:), aflux_m(:,:)
    real(kind=WP), allocatable :: dttf_u(:,:), dttf_m(:,:), dttf_v(:,:), flux_v(:,:)
    real(kind=WP) :: lon, lat

    call get_environment_variable('FESOM3_MESH_DIR', mesh_dir)
    if (len_trim(mesh_dir) == 0) &
        mesh_dir = '/home/a/a270088/port2/fesom2/tests/data/MESHES/pi'
    call get_environment_variable('FESOM3_ADVHOR_OUT', out_path)
    if (len_trim(out_path) == 0) out_path = 'advhor_f3.bin'

    call par_init(partit)
    if (partit%npes /= 1) then
        if (partit%mype == 0) write(*,'(a)') 'fesom_advhordump: requires 1 rank'
        call par_ex(partit%MPI_COMM_FESOM, partit%mype); error stop 1
    end if

    call read_mesh(mesh, partit, trim(mesh_dir), 50.0_WP, 15.0_WP, -90.0_WP, 360.0_WP, &
                   force_rotation=.true., n_cw_swaps=nsw)
    call compute_geometry(mesh, partit, cartesian=.false.)
    nl = mesh%nl
    write(*,'(a,i0,a,i0,a,i0,a,i0,a,i0)') 'fesom_advhordump: nod2D=', mesh%nod2D, &
        ' elem2D=', mesh%elem2D, ' edge2D=', mesh%edge2D, ' nl=', nl, ' CW swaps=', nsw

    ! helem (linfs, full cells): layer thickness at element = zbar(nz)-zbar(nz+1)
    ! for nz in [ulevels(elem), nlevels(elem)-1]. (FESOM2 init_thickness_ale linfs,
    ! bottom_elem_thickness = zbar(nle-1)-zbar(nle) for full cells; verified by gate.)
    allocate(mesh%helem(nl-1, mesh%elem2D)); mesh%helem = 0.0_MP
    do e = 1, mesh%elem2D
        do nz = mesh%ulevels(e), mesh%nlevels(e)-1
            mesh%helem(nz, e) = mesh%zbar(nz) - mesh%zbar(nz+1)
        end do
    end do

    ! Prescribe analytic velocity (at elements) and tracer (at nodes). MUST match
    ! FESOM2 src/fesom_advhor_dump.F90 exactly:
    !   ttf(nz,n)   = sin(2*lon_n)*cos(lat_n) + 0.02*nz
    !   vel(1,nz,e) = 0.20*cos(lat_e)*sin(lon_e) + 0.001*nz
    !   vel(2,nz,e) = 0.15*cos(2*lat_e)
    ! lon/lat are the rotated coords; lat_e/lon_e use the element's 1st node.
    allocate(ttf(nl-1, mesh%nod2D), vel(2, nl-1, mesh%elem2D))
    do n = 1, mesh%nod2D
        lon = mesh%coord_nod2D(1, n); lat = mesh%coord_nod2D(2, n)
        do nz = 1, nl-1
            ttf(nz, n) = sin(2.0_WP*lon)*cos(lat) + 0.02_WP*real(nz, WP)
        end do
    end do
    do e = 1, mesh%elem2D
        lon = mesh%coord_nod2D(1, mesh%elem2D_nodes(1, e))
        lat = mesh%coord_nod2D(2, mesh%elem2D_nodes(1, e))
        do nz = 1, nl-1
            vel(1, nz, e) = 0.20_WP*cos(lat)*sin(lon) + 0.001_WP*real(nz, WP)
            vel(2, nz, e) = 0.15_WP*cos(2.0_WP*lat)
        end do
    end do

    ! MUSCL setup: nboundary_lay + edge_up_dn_tri (+ zero edge_up_dn_grad)
    call muscl_adv_init(twork, mesh)

    ! Elemental tracer gradient, then per-edge up/downwind gradient
    allocate(tr_xy(2, nl-1, mesh%elem2D))
    call tracer_gradient_elements(ttf, tr_xy, mesh)
    call fill_up_dn_grad(twork, tr_xy, mesh)
    allocate(eudg(4, nl-1, mesh%edge2D))
    eudg = real(twork%edge_up_dn_grad, WP)

    allocate(aflux_u(nl-1, mesh%edge2D), aflux_m(nl-1, mesh%edge2D))
    allocate(dttf_u(nl-1, mesh%nod2D), dttf_m(nl-1, mesh%nod2D))
    allocate(dttf_v(nl-1, mesh%nod2D), flux_v(nl, mesh%nod2D))
    flux_v = 0.0_WP

    ! --- UPW1 horizontal flux + scatter -> del_ttf_advhoriz ---
    call adv_tra_hor_upw1(vel, ttf, mesh, aflux_u, o_init_zero=.true.)
    dttf_u = 0.0_WP; dttf_v = 0.0_WP
    call oce_tra_adv_flux2dtracer(dt, dttf_u, dttf_v, aflux_u, flux_v, mesh)

    ! --- MUSCL horizontal flux + scatter -> del_ttf_advhoriz ---
    call adv_tra_hor_muscl(vel, ttf, mesh, num_ord, aflux_m, eudg, twork%nboundary_lay, o_init_zero=.true.)
    dttf_m = 0.0_WP; dttf_v = 0.0_WP
    call oce_tra_adv_flux2dtracer(dt, dttf_m, dttf_v, aflux_m, flux_v, mesh)

    ! --- dump (same order/names as the FESOM2 oracle) ---
    call advhor_dump_open(u, trim(out_path), mesh%nod2D, mesh%elem2D, mesh%edge2D, nl)
    call wr_r2(u, 'ttf',                    real(ttf, MP))
    call wr_r3(u, 'vel',                    real(vel, MP))
    call wr_r2(u, 'helem',                  mesh%helem(1:nl-1, 1:mesh%elem2D))
    call wr_i1(u, 'nboundary_lay',          twork%nboundary_lay(1:mesh%nod2D))
    call wr_i2(u, 'edge_up_dn_tri',         twork%edge_up_dn_tri(1:2, 1:mesh%edge2D))
    call wr_r3(u, 'tr_xy',                  real(tr_xy, MP))
    call wr_r3(u, 'edge_up_dn_grad',        twork%edge_up_dn_grad(1:4, 1:nl-1, 1:mesh%edge2D))
    call wr_r2(u, 'adv_flux_hor_upw1',      real(aflux_u, MP))
    call wr_r2(u, 'del_ttf_advhoriz_upw1',  real(dttf_u, MP))
    call wr_r2(u, 'adv_flux_hor_muscl',     real(aflux_m, MP))
    call wr_r2(u, 'del_ttf_advhoriz_muscl', real(dttf_m, MP))
    call advhor_dump_close(u)
    write(*,'(a)') 'fesom_advhordump: wrote '//trim(out_path)

    call par_ex(partit%MPI_COMM_FESOM, partit%mype)
end program fesom_advhordump
