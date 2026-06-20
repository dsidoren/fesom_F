module oce_mo_conv
    ! Convective adjustment + (deferred) Monin-Obukhov / wind-mixing enhancements.
    ! Transcribed from FESOM2 v2.7.3 oce_mo_conv.F90:4-121 (subroutine mo_convect).
    !
    ! M2.8b SCOPE — the convective (static-instability) adjustment, called AFTER PP
    ! mixing in the step (FESOM2 oce_ale.F90:3729). Where N²<0 the column is statically
    ! unstable, so floor the vertical diffusivity Kv (nodes) and viscosity Av (elements)
    ! to instabmix_kv (=0.1, strong mixing) to mimic convective overturning:
    !
    !   Kv(nz,node) = max(Kv(nz,node), instabmix_kv)   where  bvfreq(nz,node)   < 0   (use_instabmix)
    !   Av(nz,elem) = max(Av(nz,elem), instabmix_kv)   where  any(bvfreq(nz,elnodes)<0) (use_instabmix)
    !
    ! reads:  dyn%work%bvfreq (nl,nod2D, the SMOOTHED N² the step computed in pressure_bv),
    !         mesh level arrays. in/out: dyn%work%Kv (nl,nod2D) + dyn%work%Av (nl,elem2D),
    !         modified IN PLACE (the PP output is the incoming value).
    ! params: use_instabmix, instabmix_kv, use_windmix, windmix_kv, windmix_nl.
    !
    ! DEFERRED to M2.10 (forcing): the use_momix Monin-Obukhov / TB04 block (FESOM2
    ! oce_mo_conv.F90:34-74 + the Av momix term :114) reads forcing/ice fields not ported
    ! yet (water_flux / heat_flux / stress_node_surf / u_ice/v_ice/a_ice / the mo /
    ! mixlength arrays / mo_length / pmlktmo) — OMITTED here (NOT just guarded; the fields
    ! do not exist). pi runs use_momix=.true. in production, but the M2.8b gate FORCES
    ! use_momix=.false. (pre-forcing). The use_windmix block (uses only params + Kv/Av) IS
    ! transcribed, guarded off (use_windmix=.false. — the FESOM2 default; pi does not set it).
    !
    ! 1-rank: node loops 1..nod2D (FESOM2 myDim+eDim, eDim=0); elem loop 1..elem2D (FESOM2
    ! myDim). No halo exchange (purely local). Byte-identical by L9 transitivity (bvfreq
    ! gated M2.1; Kv/Av the gated PP output; max() is order-free). elem2D_nodes MAX_NV=4 in
    ! FESOM3 -> slice (1:3) (the L15 trap). cavity nzmin>1 path transcribed, ungated on pi.
    use mod_precision,  only: WP
    use mod_mesh,       only: t_mesh
    use mod_dyn,        only: t_dyn
    use mod_param_phys, only: use_instabmix, instabmix_kv, use_windmix, windmix_kv, windmix_nl
    implicit none
    private
    public :: mo_convect

contains

    subroutine mo_convect(dyn, mesh)
        type(t_mesh), intent(in)            :: mesh
        type(t_dyn),  intent(inout), target :: dyn
        integer       :: node, elem, nz, nzmin, nzmax, elnodes(3)
        real(kind=WP), dimension(:,:), pointer :: Kv, Av, bvfreq

        Kv     => dyn%work%Kv
        Av     => dyn%work%Av
        bvfreq => dyn%work%bvfreq

        !_______________________________________________________________________
        ! enhance vertical DIFFUSIVITY (nodes): convective adjustment + wind mixing.
        do node = 1, mesh%nod2D
            nzmax = mesh%nlevels_nod2D(node)
            nzmin = mesh%ulevels_nod2D(node)
            do nz = nzmin+1, nzmax-1
                ! force convection where statically unstable (N²<0)
                if (use_instabmix .and. bvfreq(nz,node) < 0._WP) Kv(nz,node) = max(Kv(nz,node), instabmix_kv)
                ! cavity: no wind mixing
                if (nzmin > 1) cycle
                ! enhanced near-surface wind mixing (FESOM1.4 style)
                if (use_windmix .and. nz <= windmix_nl+1) Kv(nz,node) = max(Kv(nz,node), windmix_kv)
            end do
        end do

        !_______________________________________________________________________
        ! enhance vertical VISCOSITY (elements): convective adjustment + wind mixing.
        ! (the use_momix Monin-Obukhov viscosity term is DEFERRED to M2.10 — see header.)
        do elem = 1, mesh%elem2D
            elnodes = mesh%elem2D_nodes(1:3,elem)
            nzmax = mesh%nlevels(elem)
            nzmin = mesh%ulevels(elem)
            do nz = nzmin+1, nzmax-1
                ! force convection where ANY node of the element is statically unstable
                if (use_instabmix .and. any(bvfreq(nz,elnodes) < 0._WP)) Av(nz,elem) = max(Av(nz,elem), instabmix_kv)
                ! cavity: no Monin-Obukhov / wind mixing
                if (nzmin > 1) cycle
                ! enhanced near-surface wind mixing
                if (use_windmix .and. nz <= windmix_nl+1) Av(nz,elem) = max(Av(nz,elem), windmix_kv)
            end do
        end do
    end subroutine mo_convect

end module oce_mo_conv
