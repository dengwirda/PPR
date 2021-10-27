
!   gfortran -pedantic -cpp -O3 -flto ex_6.f90 -o ex_6
!   ./ex_6

!   1d scalar transport in a periodic domain. Fluxes are
!   computed in a flux-form semi-lagrangian sense -- 
!   integrating over the upwind regions "swept" by edges
!   in each time-step. This implementation requires CFL<1,
!   with the upwind regions covering adjacent cells only. 
!   A CFL>=1 variant could be constructed using RMAP1D().
!

#   include "../src/ppr_1d.f90"

    program ex

        use ppr_1d

        implicit none

        integer, parameter :: halo = 4  ! halo cells at boundary
        integer, parameter :: npos = 51 ! no. edge 
        integer, parameter :: nvar = 1  ! no. tracers to trnsprt
        integer, parameter :: ndof = 1  ! no. FV DoF per cell
        integer :: ipos,il,ir,step

    !------------------------------- domain discretisation !
        real(kind=dp) :: xpos(1-halo:npos+halo)
        real(kind=dp) :: xmid,xdel(1),tDEL
        
    !-------------------------------- finite-volume arrays !

    !   INIT: initial cell-wise finite-volume profile
    !   QBAR: dynamic cell-wise finite-volume profile
    !   MASK: cell-wise "land" mask (all TRUE here)
    !   UVEL: edge-wise velocity distribution
    !   FLUX: edge-wise distribution of upwind fluxes
    !   QDIV: cell-wise distribution of divergence

        real(kind=dp) :: init(ndof,nvar,1-halo:npos+halo-1)
        real(kind=dp) :: qbar(ndof,nvar,1-halo:npos+halo-1)
        logical       :: mask(          1-halo:npos+halo-1)
        real(kind=dp) :: uvel(          1-halo:npos+halo)        
        real(kind=dp) :: flux(     nvar,1-halo:npos+halo)        
        real(kind=dp) :: qdiv(     nvar,1-halo:npos+halo-1)
        

    !------------------------------ method data-structures !
        type(rmap_work) :: work
        type(rmap_opts) :: opts
        type(rcon_ends) :: bc_l(nvar)
        type(rcon_ends) :: bc_r(nvar)

    !------------------------------ define a simple domain !

        il =      + 1                   ! 1st real interior cell
        ir = npos - 1                   ! Nth real interior cell

        xpos(il+0) = 0.0d+00
        xpos(ir+1) = 1.0d+00
        
        xdel(1) = (xpos(ir+1)-xpos(il+0))/(npos-1)

        do ipos = il+1, ir-0

            xpos(ipos) = (ipos-1) * xdel(1)         

        end do

    !------------------------------ setup some simple data !

        uvel       = +1.0d+0
        mask       = .true.

        tDEL       = +1.0d-2

        do ipos = +1, npos-1

            xmid = xpos(ipos+0)* 0.5d+00 &
    &            + xpos(ipos+1)* 0.5d+00

            init(1,1,ipos) = &
    &   .8d+0 * exp( -75.0d+0 * (xmid - 0.275d+0) ** 2 ) &
    & + .9d+0 * exp(-100.0d+0 * (xmid - 0.500d+0) ** 2 ) &
    & + 1.d+0 * exp(-125.0d+0 * (xmid - 0.725d+0) ** 2 )
            
        end do

    !------------------------------ specify method options !

        opts%edge_meth = p3e_method     ! 3rd-order edge interp.
        opts%cell_meth = ppm_method     ! PPM method in cells
        opts%cell_lims = weno_limit     ! "non-oscillatory" lim.
        
    !------------------------------ set BC.'s at endpoints !

        bc_l%bcopt = bcon_loose         ! "loose" = extrapolate
        bc_r%bcopt = bcon_loose

    !------------------------------ init. method workspace !

        call work%init(npos+2*halo,nvar,opts)

    !------------------------------ calc. scalar transport !

        qbar = init

        do step = +1, +100        ! 100 steps => full loop

    !------------------------------ periodicity via halo's !

            qbar(1,:,il-halo:il-1) = & 
    &       qbar(1,:,ir-halo+1:ir)
            qbar(1,:,ir+1:ir+halo) = &
    &       qbar(1,:,il:il+halo-1)

    !------------------------------ form lagrangian fluxes !

            call ffsl1d (npos+2*halo,nvar, &
    &            ndof,xdel,tDEL,mask,uvel, &
    &            qbar,flux,bc_l,bc_r,work, &
    &            opts)

  
    !------------------------------ flux divergence eval's !
            qdiv(  1,il:ir) = &
    &              flux(1,il+1:ir+1) &
    &            - flux(1,il+0:ir+0)

    !------------------------------ take a single timestep !

            qbar(1,1,il:ir) = &
    &           qbar(1,1,il:ir) - &
    &               qdiv(1,il:ir) / xdel(1)

        end do

    !------------------------------ clear method workspace !

        call work%free()

    !------------------------------ dump results to stdout !

        print*,"End timestep profile : "

        do ipos = il+0, ir-0

            print *, init(1,:,ipos)  &
    &              , qbar(1,:,ipos)

        end do

        print*,"Conservation defect := " &
    &         , sum(init(1,:,il:ir)) &
    &         - sum(qbar(1,:,il:ir))
        
    end program



