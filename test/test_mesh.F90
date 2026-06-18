program test_mesh
    ! M0.7 gate (self-consistency; the FESOM2 byte-gate is deferred to M1, needs an
    ! instrumented-FESOM2 reference geometry dump). Loads pi (1-rank) and an analytic
    ! mesh, and checks geometric identities that catch transcription bugs:
    !   - clockwise orientation enforced (no element with positive signed area)
    !   - elem_area > 0 and Sum ~ sphere surface (pi) / cell area (analytic)
    !   - gradient_sca annihilates a constant field (Sum of shape-fn gradients = 0)
    !   - on the Cartesian analytic mesh, gradient_sca reproduces a linear field exactly
    use mpi
    use mod_precision,    only: WP, MP
    use mod_constants,    only: r_earth, pi
    use mod_mesh,         only: t_mesh
    use mod_partit,       only: t_partit
    use mod_partitioning, only: par_init, par_ex
    use mod_mesh_read,    only: read_mesh
    use mod_mesh_areas,   only: compute_geometry
    use mod_mesh_analytic, only: generate_analytic_mesh
    use mod_mesh_rotate,  only: trim_cyclic
    implicit none

    character(len=512) :: mesh_dir
    type(t_partit) :: partit
    type(t_mesh)   :: mesh, amesh
    integer :: nfail, nsw, ierr

    call get_environment_variable('FESOM3_MESH_DIR', mesh_dir)
    if (len_trim(mesh_dir) == 0) &
        mesh_dir = '/home/a/a270088/port2/fesom2/tests/data/MESHES/pi'
    nfail = 0
    call par_init(partit)
    if (partit%npes /= 1) then
        if (partit%mype == 0) write(*,'(a)') 'test_mesh: requires 1 rank (M0.7)'
        call par_ex(partit%MPI_COMM_FESOM, partit%mype); error stop 1
    end if

    ! ================= pi (rotated sphere) =================
    call read_mesh(mesh, partit, trim(mesh_dir), 50.0_WP, 15.0_WP, -90.0_WP, 360.0_WP, &
                   force_rotation=.true., n_cw_swaps=nsw)
    call compute_geometry(mesh, partit, cartesian=.false.)
    write(*,'(a,i0,a,i0,a,i0)') 'pi: nod2D=', mesh%nod2D, ' elem2D=', mesh%elem2D, ' nl=', mesh%nl
    write(*,'(a,i0)') 'pi: CW swaps at load = ', nsw
    write(*,'(a,es12.5,a,es12.5)') 'pi: sum(elem_area)=', sum(mesh%elem_area), &
        '  sphere=', 4.0_MP*real(pi,MP)*r_earth*r_earth

    call check(mesh%nod2D == 3140 .and. mesh%elem2D == 5839 .and. mesh%nl == 48, 'pi dims')
    call check(no_positive_orientation(mesh), 'pi: all elements clockwise after enforce')
    call check(all(mesh%elem_area > 0.0_MP), 'pi: elem_area > 0')
    ! pi is OCEAN-only (no land elements); sum(elem_area) ~ ocean fraction of the
    ! sphere (~0.5-0.85 of 4*pi*r^2), not the full sphere.
    call check(ocean_fraction_ok(mesh), 'pi: sum(elem_area) ~ ocean fraction of sphere')
    call check(max_grad_const(mesh) < 1.0e-6_WP, 'pi: gradient_sca annihilates constant')
    call check(all(mesh%nlevels >= 1 .and. mesh%nlevels <= mesh%nl), 'pi: nlevels in range')
    call check(adjacency_consistent(mesh), 'pi: nod_in_elem2D consistent')
    call check(all(mesh%area(1, 1:mesh%nod2D) > 0.0_MP), 'pi: surface control areas > 0')

    ! ================= analytic (Cartesian) =================
    call generate_analytic_mesh(amesh, partit, nx=9, ny=7, nl=5, &
                                Lx=8000.0_WP, Ly=6000.0_WP, max_depth=100.0_WP)
    write(*,'(a,i0,a,i0,a,i0)') 'analytic: nod2D=', amesh%nod2D, ' elem2D=', amesh%elem2D, &
        ' edge2D=', amesh%edge2D
    call check(amesh%elem2D == 2*8*6, 'analytic: elem count')
    call check(amesh%edge2D > 0 .and. amesh%edge2D_in > 0, 'analytic: edges built')
    call check(all(amesh%elem_area > 0.0_MP), 'analytic: elem_area > 0')
    call check(max_grad_const(amesh) < 1.0e-6_WP, 'analytic: gradient annihilates constant')
    call check(linear_gradient_exact(amesh, 3.0_WP, 7.0_WP) < 1.0e-6_WP, &
               'analytic: gradient_sca reproduces linear field (3,7)')

    if (nfail == 0) then
        write(*,'(a)') 'test_mesh: ALL PASS'
    else
        write(*,'(a,i0,a)') 'test_mesh: ', nfail, ' FAILURE(S)'
    end if
    call par_ex(partit%MPI_COMM_FESOM, partit%mype)
    if (nfail /= 0) error stop 1

contains

    subroutine check(cond, name)
        logical, intent(in) :: cond
        character(len=*), intent(in) :: name
        if (.not. cond) then
            nfail = nfail + 1
            write(*,'(a)') '  FAIL: '//name
        end if
    end subroutine

    logical function no_positive_orientation(m)
        type(t_mesh), intent(in) :: m
        integer :: n
        real(kind=WP) :: a1, a2, b1, b2, c1, c2, r
        no_positive_orientation = .true.
        do n = 1, m%elem2D
            a1 = m%coord_nod2D(1, m%elem2D_nodes(1,n)); a2 = m%coord_nod2D(2, m%elem2D_nodes(1,n))
            b1 = m%coord_nod2D(1, m%elem2D_nodes(2,n)) - a1; b2 = m%coord_nod2D(2, m%elem2D_nodes(2,n)) - a2
            c1 = m%coord_nod2D(1, m%elem2D_nodes(3,n)) - a1; c2 = m%coord_nod2D(2, m%elem2D_nodes(3,n)) - a2
            call trim_cyclic(b1); call trim_cyclic(c1)
            r = b1*c2 - b2*c1
            if (r > 1.0e-12_WP) no_positive_orientation = .false.
        end do
    end function

    logical function ocean_fraction_ok(m)
        type(t_mesh), intent(in) :: m
        real(kind=MP) :: frac
        frac = sum(m%elem_area) / (4.0_MP*real(pi,MP)*r_earth*r_earth)
        ocean_fraction_ok = (frac > 0.5_MP) .and. (frac < 0.85_MP)
    end function

    real(kind=WP) function max_grad_const(m)
        ! gradient of f==1 should be ~0 (shape-function gradients sum to zero).
        type(t_mesh), intent(in) :: m
        integer :: e, j
        real(kind=WP) :: gx, gy
        max_grad_const = 0.0_WP
        do e = 1, m%elem2D
            gx = 0.0_WP; gy = 0.0_WP
            do j = 1, 3
                gx = gx + 1.0_WP * real(m%gradient_sca(j,   e), WP)
                gy = gy + 1.0_WP * real(m%gradient_sca(3+j, e), WP)
            end do
            max_grad_const = max(max_grad_const, abs(gx), abs(gy))
        end do
    end function

    real(kind=WP) function linear_gradient_exact(m, a, b)
        ! On a Cartesian mesh, grad of f = a*x + b*y must be (a,b) at every element.
        type(t_mesh), intent(in) :: m
        real(kind=WP), intent(in) :: a, b
        integer :: e, j, nd
        real(kind=WP) :: gx, gy, fnode
        linear_gradient_exact = 0.0_WP
        do e = 1, m%elem2D
            gx = 0.0_WP; gy = 0.0_WP
            do j = 1, 3
                nd = m%elem2D_nodes(j, e)
                ! coords are radians-like; physical position = coord * r_earth
                fnode = a*real(m%coord_nod2D(1, nd), WP)*r_earth &
                      + b*real(m%coord_nod2D(2, nd), WP)*r_earth
                gx = gx + fnode * real(m%gradient_sca(j,   e), WP)
                gy = gy + fnode * real(m%gradient_sca(3+j, e), WP)
            end do
            linear_gradient_exact = max(linear_gradient_exact, abs(gx-a), abs(gy-b))
        end do
    end function

    logical function adjacency_consistent(m)
        type(t_mesh), intent(in) :: m
        integer :: n, j, e, k
        logical :: found
        adjacency_consistent = .true.
        do n = 1, m%nod2D
            do j = 1, m%nod_in_elem2D_num(n)
                e = m%nod_in_elem2D(j, n)
                found = .false.
                do k = 1, m%elem2D_nnodes(e)
                    if (m%elem2D_nodes(k, e) == n) found = .true.
                end do
                if (.not. found) adjacency_consistent = .false.
            end do
        end do
    end function

end program test_mesh
