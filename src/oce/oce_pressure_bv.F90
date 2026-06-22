module oce_pressure_bv
    ! Hydrostatic pressure + Brunt-Vaisala frequency in a single pass, using the
    ! SPLIT form of the Jackett-McDougall equation of state. Transcribed from
    ! FESOM2 v2.7.3 oce_ale_pressure_bv.F90:
    !   pressure_bv          (194-501)  density_m_rho0 + hpressure + bvfreq
    !   densityJM_components (2605-2669) split-form secant bulk modulus + rhopot
    !   smooth_nod3D         (gen_support.F90:99-198) horizontal N^2 smoother
    !
    ! M2.1 SCOPE — the first dynamics kernel; gate target = density_m_rho0 (PGF
    ! density anomaly), hpressure (top-down hydrostatic integration) and bvfreq
    ! (N^2, squared). OMITTED (pure diagnostics that do NOT feed the gated fields):
    ! the dMOC density (ldiag_dMOC), the KPP buoyancy diff dbsfc1/dbsfc, db_max,
    ! and the mixed-layer depths MLD1/2/3. Dropping them leaves density_m_rho0 /
    ! hpressure / bvfreq bit-for-bit unchanged (density_m_rho0 is computed BEFORE
    ! rho_surf/dbsfc1/db_max; the MLD logic only READS bvfreq/rhopot/Z_3d_n).
    !
    ! EOS NOTES (bit-identity critical):
    !  - NEVER linearize alpha/beta — the split form factorizes the secant bulk
    !    modulus (bulk_0 + pz*(bulk_pz + pz*bulk_pz2)) and the bits depend on that
    !    factorization. Transcribe densityJM_components verbatim.
    !  - density anomaly subtracts the density_ref(nz,node) ARRAY (here = density_0
    !    everywhere, since use_density_ref=.false. on pi; init_ref_density is the
    !    deferred use_density_ref=.true. path). N^2 divides by the SCALAR density_0.
    !  - real(state_equation): 1 (full EOS, case 1). Kept a variable so the
    !    0.1*Z*real(state_equation) terms transcribe verbatim (state_equation=0
    !    linear EOS would zero the pressure term — not gated here).
    !
    ! which_ALE='linfs' (pi) routes the hydrostatic pressure integration through
    ! THIS routine (the zstar/zlevel cases compute it in pressure_force_4_zxxxx,
    ! M2.x). use_cavity / partial cells (nzmin>1, the cavity branches below) are
    ! v1-dropped and UNGATED on pi (nzmin==1, full cells) — transcribed faithfully
    ! for the M2.11 cavity re-gate.
    !
    ! 1-rank only: smooth_nod's per-cycle exchange_nod(bvfreq) is a no-op here and
    ! is dropped (lifted with the other M1/M2 halo exchanges at M2.12). The three
    ! outputs are intent(inout) and NOT zeroed inside (FESOM2 pressure_bv does not
    ! initialize the below-bottom entries either) — the CALLER must zero them so the
    ! dump's below-bottom region is a deterministic 0 on both sides of the gate.
    use mod_precision,   only: WP
    use mod_mesh,        only: t_mesh
    use mod_constants,   only: density_0, g
    use mod_config,      only: which_ALE
    use mod_param_phys,  only: state_equation, N2smth_h, N2smth_v, N2smth_hidx
    use mod_partit,      only: t_partit
    use mod_part_bounds, only: owned_bounds, is_multirank
    use mod_halo,        only: exchange_nod
    implicit none
    private
    public :: pressure_bv, densityJM_components, insitu2pot

contains

    !===========================================================================
    subroutine pressure_bv(temp, salt, density_ref, mesh, density_m_rho0, hpressure, bvfreq, partit)
        ! temp/salt/density_ref: (nl-1, nod2D) inputs. density_m_rho0: (nl-1, nod2D)
        ! out. hpressure/bvfreq: (nl, nod2D) out (only 1..nzmax used). The caller
        ! pre-zeros the three outputs (see header).
        ! M2.12c: optional partit -> owned+halo node loop (FESOM2 :237-244 do node=1,
        ! myDim_nod2D+eDim_nod2D — the EOS is per-node so computing the halo here saves
        ! an exchange; the downstream pgf/smoothing then read it locally) + smooth_nod's
        ! per-sweep exchange_nod(bvfreq).
        type(t_mesh),  intent(in)    :: mesh
        real(kind=WP), intent(in)    :: temp(mesh%nl-1, mesh%nod2D)
        real(kind=WP), intent(in)    :: salt(mesh%nl-1, mesh%nod2D)
        real(kind=WP), intent(in)    :: density_ref(mesh%nl-1, mesh%nod2D)
        real(kind=WP), intent(inout) :: density_m_rho0(mesh%nl-1, mesh%nod2D)
        real(kind=WP), intent(inout) :: hpressure(mesh%nl, mesh%nod2D)
        real(kind=WP), intent(inout) :: bvfreq(mesh%nl, mesh%nod2D)
        type(t_partit), intent(in), optional :: partit

        integer       :: node, nz, nzmax, nzmin
        integer       :: nNodO, nNodL, nEdgeO, nElemO
        real(kind=WP) :: zmean, dz_inv, a, rho_up, rho_dn, t, s, smin
        real(kind=WP) :: bulk_up, bulk_dn
        real(kind=WP) :: rhopot(mesh%nl), bulk_0(mesh%nl), bulk_pz(mesh%nl)
        real(kind=WP) :: bulk_pz2(mesh%nl), rho(mesh%nl), bv1(mesh%nl)

        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)

        !_______________________________________________________________________
        ! Screen salinity (diagnostic; s<0 would break sqrt(s) in the EOS).
        smin = 0.0_WP
        do node=1, nNodL
            do nz = mesh%ulevels_nod2D(node), mesh%nlevels_nod2D(node)-1
                smin = min(smin, salt(nz,node))
            end do
        end do
        if (smin < 0.0_WP) write(*,*) ' --> oce_pressure_bv: s<0 happens!', smin

        !_______________________________________________________________________
        do node=1, nNodL
            nzmin = mesh%ulevels_nod2D(node)
            nzmax = mesh%nlevels_nod2D(node)

            rho      = 0.0_WP
            bulk_0   = 0.0_WP
            bulk_pz  = 0.0_WP
            bulk_pz2 = 0.0_WP
            rhopot   = 0.0_WP

            !___________________________________________________________________
            ! apply equation of state
            do nz=nzmin, nzmax-1
                t=temp(nz, node)
                s=salt(nz, node)
                call densityJM_components(t, s, bulk_0(nz), bulk_pz(nz), bulk_pz2(nz), rhopot(nz))
            end do

            !___________________________________________________________________
            ! compute density for PGF
            do nz=nzmin, nzmax-1
                rho(nz) = bulk_0(nz) + mesh%Z_3d_n(nz,node)*(bulk_pz(nz) + mesh%Z_3d_n(nz,node)*bulk_pz2(nz))
                rho(nz) = rho(nz)*rhopot(nz)/(rho(nz)+0.1_WP*mesh%Z_3d_n(nz,node)*real(state_equation,WP))-density_ref(nz,node)
                density_m_rho0(nz,node) = rho(nz)
            end do

            !___________________________________________________________________
            ! fill density levels occupied by the cavity (nzmin>1; v1-dropped, ungated
            ! on pi). Mass that corresponds to T/S at the cavity-ocean interface.
            if (nzmin>1) then
                t=temp(nzmin, node)
                s=salt(nzmin, node)
                do nz=1, nzmin-1
                    call densityJM_components(t, s, bulk_0(nz), bulk_pz(nz), bulk_pz2(nz), rhopot(nz))
                    rho(nz)= bulk_0(nz)   + mesh%Z_3d_n(nz,node)*(bulk_pz(nz)   + mesh%Z_3d_n(nz,node)*bulk_pz2(nz))
                    rho(nz)=rho(nz)*rhopot(nz)/(rho(nz)+0.1_WP*mesh%Z_3d_n(nz,node)*real(state_equation,WP))-density_ref(nz,node)
                    density_m_rho0(nz,node) = rho(nz)
                end do
            end if

            !___________________________________________________________________
            ! calculate pressure (linfs or cavity case)
            if (trim(which_ALE)=='linfs') then
                if (nzmin>1) then ! cavity case (v1-dropped, ungated on pi)
                    hpressure(nzmin, node)=0.5_WP*(mesh%zbar_3d_n(1,node)-mesh%zbar_3d_n(2,node))*rho(1)*g
                    do nz=2,nzmin
                        a=0.5_WP*g*(rho(nz-1)*(mesh%zbar_3d_n(nz-1,node)-mesh%zbar_3d_n(nz,node))+rho(nz)*(mesh%zbar_3d_n(nz,node)-mesh%zbar_3d_n(nz+1,node)))
                        hpressure(nzmin, node)=hpressure(nzmin, node)+a
                    end do
                else
                    hpressure(nzmin, node)=-mesh%Z_3d_n(nzmin,node)*rho(nzmin)*g
                end if

                ! pressure below surface boundary: integrate g*rho*dz over half the
                ! previous + half the actual layer thickness to mid-depth of the layer
                do nz=nzmin+1,nzmax-1
                    a=0.5_WP*g*(rho(nz-1)*mesh%hnode(nz-1,node)+rho(nz)*mesh%hnode(nz,node))
                    hpressure(nz, node)=hpressure(nz-1, node)+a
                end do
            end if

            !___________________________________________________________________
            ! squared Brunt-Vaisala frequency N^2 (defined on full levels except the
            ! first and the last). N^2>0 stable, N^2<0 unstable stratification.
            do nz=nzmin+1,nzmax-1
                zmean   = 0.5_WP*sum(mesh%Z_3d_n(nz-1:nz, node))
                bulk_up = bulk_0(nz-1) + zmean*(bulk_pz(nz-1) + zmean*bulk_pz2(nz-1))
                bulk_dn = bulk_0(nz)   + zmean*(bulk_pz(nz)   + zmean*bulk_pz2(nz))
                rho_up  = bulk_up*rhopot(nz-1) / (bulk_up + 0.1_WP*zmean*real(state_equation,WP))
                rho_dn  = bulk_dn*rhopot(nz)   / (bulk_dn + 0.1_WP*zmean*real(state_equation,WP))
                dz_inv  = 1.0_WP/(mesh%Z_3d_n(nz-1,node)-mesh%Z_3d_n(nz,node))
                bvfreq(nz,node)  = -g*dz_inv*(rho_up-rho_dn)/density_0
            end do

            bvfreq(nzmin,node)=bvfreq(nzmin+1,node)
            bvfreq(nzmax,node)=bvfreq(nzmax-1,node)

            !___________________________________________________________________
            ! optional vertical N^2 smoothing (pinned .false. on the gate)
            if (N2smth_v) then
                do nz=nzmin+1,nzmax-1
                    bv1(nz)=        (mesh%zbar_3d_n(nz-1,node)-mesh%zbar_3d_n(nz,  node))*(bvfreq(nz-1,node)+bvfreq(nz,  node))
                    bv1(nz)=bv1(nz)+(mesh%zbar_3d_n(nz,  node)-mesh%zbar_3d_n(nz+1,node))*(bvfreq(nz,  node)+bvfreq(nz+1,node))
                    bv1(nz)=0.5_WP*bv1(nz)/(mesh%zbar_3d_n(nz-1,node)-mesh%zbar_3d_n(nz+1,  node))
                end do
                do nz=nzmin+1,nzmax-1
                    bvfreq(nz,node)=bv1(nz)
                end do
            end if
        end do

        !_______________________________________________________________________
        ! apply horizontal smoothing of the N^2 buoyancy frequency
        if (N2smth_h) call smooth_nod(bvfreq, N2smth_hidx, mesh, partit)

    end subroutine pressure_bv

    !===========================================================================
    subroutine densityJM_components(t, s, bulk_0, bulk_pz, bulk_pz2, rhopot)
        ! Split form of the Jackett-McDougall (1992, CSIRO) in-situ density EOS.
        ! Verbatim from FESOM2 oce_ale_pressure_bv.F90:2605-2669 (N. Rakowski split).
        ! rho_in-situ(pz) = bulk*rhopot/(bulk+0.1*pz) with bulk = bulk_0 + pz*(bulk_pz
        ! + pz*bulk_pz2). DO NOT pre-multiply/refactor — the bit pattern depends on
        ! the exact Horner factorization below.
        real(kind=WP), intent(in)  :: t, s
        real(kind=WP), intent(out) :: bulk_0, bulk_pz, bulk_pz2, rhopot
        real(kind=WP) :: s_sqrt

        real(kind=WP), parameter   :: a0    = 19092.56,     at   = 209.8925
        real(kind=WP), parameter   :: at2   = -3.041638,    at3  = -1.852732e-3
        real(kind=WP), parameter   :: at4   = -1.361629e-5
        real(kind=WP), parameter   :: as    = 104.4077,     ast  = -6.500517
        real(kind=WP), parameter   :: ast2  = .1553190,     ast3 = 2.326469e-4
        real(kind=WP), parameter   :: ass   = -5.587545,    asst = 0.7390729
        real(kind=WP), parameter   :: asst2 = -1.909078e-2
        real(kind=WP), parameter   :: ap    = -4.721788e-1, apt  = -1.028859e-2
        real(kind=WP), parameter   :: apt2  = 2.512549e-4,  apt3 = 5.939910e-7
        real(kind=WP), parameter   :: aps   = 1.571896e-2,  apst = 2.598241e-4
        real(kind=WP), parameter   :: apst2 = -7.267926e-6, apss = -2.042967e-3
        real(kind=WP), parameter   :: ap2   = 1.045941e-5,  ap2t = -5.782165e-10
        real(kind=WP), parameter   :: ap2t2 = 1.296821e-7
        real(kind=WP), parameter   :: ap2s  = -2.595994e-7,ap2st = -1.248266e-9
        real(kind=WP), parameter   :: ap2st2= -3.508914e-9

        real(kind=WP), parameter   :: b0 = 999.842594,    bt  = 6.793952e-2
        real(kind=WP), parameter   :: bt2 = -9.095290e-3, bt3 = 1.001685e-4
        real(kind=WP), parameter   :: bt4 = -1.120083e-6, bt5 = 6.536332e-9
        real(kind=WP), parameter   :: bs = 0.824493,      bst = -4.08990e-3
        real(kind=WP), parameter   :: bst2 = 7.64380e-5,  bst3 = -8.24670e-7
        real(kind=WP), parameter   :: bst4 = 5.38750e-9
        real(kind=WP), parameter   :: bss = -5.72466e-3,  bsst = 1.02270e-4
        real(kind=WP), parameter   :: bsst2 = -1.65460e-6,bss2 = 4.8314e-4

        s_sqrt = sqrt(s)

        bulk_0 =  a0      + t*(at   + t*(at2  + t*(at3 + t*at4)))      &
                + s* (as  + t*(ast  + t*(ast2 + t*ast3))               &
                     + s_sqrt*(ass  + t*(asst + t*asst2)))

        bulk_pz =  ap  + t*(apt  + t*(apt2 + t*apt3))                  &
                        + s*(aps + t*(apst + t*apst2) + s_sqrt*apss)

        bulk_pz2 = ap2 + t*(ap2t + t*ap2t2)                           &
                      + s *(ap2s + t*(ap2st + t*ap2st2))

        rhopot =  b0 + t*(bt + t*(bt2 + t*(bt3  + t*(bt4  + t*bt5))))  &
                     + s*(bs + t*(bst + t*(bst2 + t*(bst3 + t*bst4)))  &
                        + s_sqrt*(bss + t*(bsst + t*bsst2))            &
                             + s* bss2)
    end subroutine densityJM_components

    !===========================================================================
    subroutine smooth_nod(arr, N_smooth, mesh, partit)
        ! Mass-matrix horizontal smoother, transcribed from FESOM2 g_support
        ! smooth_nod3D (gen_support.F90:99-198): applies the lumped P1 mass matrix
        ! N_smooth times. M2.12c: optional partit -> owned-node loops + the per-sweep
        ! exchange_nod(arr) (the patch accumulation reads arr at an owned node's owned-
        ! element nodes, which can be halo nodes, so arr must be halo-valid each sweep).
        ! Per-level patch areas vary with the bathymetry (a deep node's shallower
        ! neighbour elements drop out level by level). elem_area is geom-proven; the
        ! nod_in_elem2D / elem2D_nodes order is area-gate-proven (L9).
        real(kind=WP), intent(inout)       :: arr(:,:)
        integer,       intent(in)          :: N_smooth
        type(t_mesh),  intent(in), target  :: mesh
        type(t_partit), intent(in), optional :: partit
        integer :: n, q, el, nz, j, nlev, uln, nln, ule, nle
        integer :: nNodO, nNodL, nEdgeO, nElemO
        real(kind=WP), allocatable :: vol(:,:), work_array(:,:)

        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)
        nlev=ubound(arr,1)
        allocate(vol(mesh%nl, nNodL), work_array(nlev, nNodL))

        ! First sweep: precompute the (inverse) patch areas, then smooth.
        do n=1, nNodO
            uln = mesh%ulevels_nod2d(n)
            nln = min(nlev, mesh%nlevels_nod2d(n))
            vol(       1:nln,n) = 0._WP
            work_array(1:nln,n) = 0._WP
            do j=1, mesh%nod_in_elem2D_num(n)
                el  = mesh%nod_in_elem2D(j,n)
                ule = max( uln, mesh%ulevels(el) )
                nle = min( nln, min(nlev, mesh%nlevels(el)) )
                do nz=ule, nle
                    vol(nz,n) = vol(nz,n) + mesh%elem_area(el)
                    work_array(nz,n) = work_array(nz,n) + mesh%elem_area(el) * (arr(nz, mesh%elem2D_nodes(1,el)) &
                                                                              + arr(nz, mesh%elem2D_nodes(2,el)) &
                                                                              + arr(nz, mesh%elem2D_nodes(3,el)))
                end do
            end do
            do nz=uln,nln
                vol(nz,n) = 1._WP / (3._WP * vol(nz,n))  ! inverse, scaled by 1/3
            end do
        end do
        do n=1, nNodO
            uln = mesh%ulevels_nod2d(n)
            nln = min(nlev, mesh%nlevels_nod2d(n))
            do nz=uln,nln
                arr(nz, n) = work_array(nz, n) * vol(nz,n)
            end do
        end do
        if (is_multirank(partit)) call exchange_nod(arr, partit)

        ! Remaining sweeps reuse the precomputed inverse patch areas.
        do q=1,N_smooth-1
            do n=1, nNodO
                uln = mesh%ulevels_nod2d(n)
                nln = min(nlev, mesh%nlevels_nod2d(n))
                work_array(1:nln,n) = 0._WP
                do j=1, mesh%nod_in_elem2D_num(n)
                    el  = mesh%nod_in_elem2D(j,n)
                    ule = max( uln, mesh%ulevels(el) )
                    nle = min( nln, min(nlev, mesh%nlevels(el)) )
                    do nz=ule, nle
                        work_array(nz,n) = work_array(nz,n) + mesh%elem_area(el) * (arr(nz, mesh%elem2D_nodes(1,el)) &
                                                                                  + arr(nz, mesh%elem2D_nodes(2,el)) &
                                                                                  + arr(nz, mesh%elem2D_nodes(3,el)))
                    end do
                end do
            end do
            do n=1, nNodO
                uln = mesh%ulevels_nod2d(n)
                nln = min(nlev, mesh%nlevels_nod2d(n))
                do nz=uln,nln
                    arr(nz, n) = work_array(nz, n) * vol(nz,n)
                end do
            end do
            if (is_multirank(partit)) call exchange_nod(arr, partit)
        end do

        deallocate(vol, work_array)
    end subroutine smooth_nod

    !===========================================================================
    ! Convert in-situ temperature -> potential temperature, IN PLACE (temp), using
    ! salinity (read-only) at reference pressure pr=0. Transcribed from FESOM2
    ! oce_ale_pressure_bv.F90 insitu2pot (3074-3118): per-node loop, the in-situ T
    ! at each level is replaced by ptheta(S, T, |Z(nz)|, 0). FESOM2 uses the 1-D
    ! mid-depth Z(nz) for the pressure proxy pp=abs(Z(nz)) (NOT Z_3d_n — the comment
    ! at :3108 keeps Z for partial-cell stability at init). The IC driver feeds
    ! tracers%data(1)%values (in-situ T) / data(2)%values (S) here (M2.11b do_ic3d,
    ! t_insitu=.true.). M2.12-MVP: optional partit -> the FESOM2 owned+halo node loop
    ! (n=1..myDim+eDim); per-node independent so the halo potential T is computed locally
    ! (= the owner's value, no exchange needed). Absent partit -> n=1..nod2D (1-rank).
    subroutine insitu2pot(temp, salt, mesh, partit)
        type(t_mesh),  intent(in)    :: mesh
        real(kind=WP), intent(inout) :: temp(mesh%nl-1, mesh%nod2D)
        real(kind=WP), intent(in)    :: salt(mesh%nl-1, mesh%nod2D)
        type(t_partit), intent(in), optional :: partit
        integer       :: n, nz, nzmin, nzmax, nNodO, nNodL, nEdgeO, nElemO
        real(kind=WP) :: pp, pr, tt, ss
        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)
        pr = 0.0_WP
        do n = 1, nNodL
            nzmin = mesh%ulevels_nod2D(n)
            nzmax = mesh%nlevels_nod2D(n)
            do nz = nzmin, nzmax-1
                tt = temp(nz, n)
                ss = salt(nz, n)
                pp = abs(mesh%Z(nz))
                temp(nz, n) = ptheta(ss, tt, pp, pr)
            end do
        end do
    end subroutine insitu2pot

    !===========================================================================
    ! Local potential temperature at reference pressure pr, via Bryden-1973
    ! adiabatic lapse rate + 4th-order Runge-Kutta. Verbatim from FESOM2
    ! oce_ale_pressure_bv.F90 ptheta (2674-2714). Args are mutated locally (t,p);
    ! the caller passes scalar temporaries so the mutation is harmless. checkvalue:
    ! theta = 36.89073 C for s=40, t=40, p=10000, pr=0.
    function ptheta(s, t, p, pr) result(theta)
        real(kind=WP) :: theta
        real(kind=WP) :: s, t, p, pr
        real(kind=WP) :: h, xk, q
        h  = pr - p
        xk = h*atg(s, t, p)
        t  = t + 0.5_WP*xk
        q  = xk
        p  = p + 0.5_WP*h
        xk = h*atg(s, t, p)
        t  = t + 0.29289322_WP*(xk-q)
        q  = 0.58578644_WP*xk + 0.121320344_WP*q
        xk = h*atg(s, t, p)
        t  = t + 1.707106781_WP*(xk-q)
        q  = 3.414213562_WP*xk - 4.121320344_WP*q
        p  = p + 0.5_WP*h
        xk = h*atg(s, t, p)
        theta = t + (xk-2.0_WP*q)/6.0_WP
    end function ptheta

    !===========================================================================
    ! Adiabatic temperature gradient deg C / decibar (Bryden 1973). Verbatim from
    ! FESOM2 oce_ale_pressure_bv.F90 atg (2719-2746). checkvalue: atg=3.255976e-4
    ! C/dbar for s=40, t=40, p=10000.
    function atg(s, t, p) result(g_atg)
        real(kind=WP) :: g_atg
        real(kind=WP) :: s, t, p, ds
        ds = s - 35.0_WP
        g_atg = (((-2.1687e-16_WP*t+1.8676e-14_WP)*t-4.6206e-13_WP)*p   &
              +((2.7759e-12_WP*t-1.1351e-10_WP)*ds+((-5.4481e-14_WP*t        &
              +8.733e-12_WP)*t-6.7795e-10_WP)*t+1.8741e-8_WP))*p             &
              +(-4.2393e-8_WP*t+1.8932e-6_WP)*ds                          &
              +((6.6228e-10_WP*t-6.836e-8_WP)*t+8.5258e-6_WP)*t+3.5803e-5_WP
    end function atg

end module oce_pressure_bv
