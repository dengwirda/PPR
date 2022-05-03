
!   gfortran -pedantic -cpp -O3 -flto ex_2.f90 -o ex_2
!   ./ex_2

!   Test for the monotone limiter: remap a stairstep profile
!   between randomly perturbed grids. When MONO-LIMIT is
!   selected, all methods should lead to monotone behaviour. 
!

#   include "../src/ppr_1d.f90"

    program ex

        use ppr_1d

        implicit none

        integer, parameter :: npos = 97 ! no. edge (old grid) 
        integer, parameter :: ntmp = 77 ! no. edge (new grid)
        integer, parameter :: nvar = 1  ! no. variables to remap
        integer, parameter :: ndof = 1  ! no. FV DoF per cell
        integer :: ipos

    !------------------------------ position of cell edges !
        real(kind=dp) :: xpos(npos),xtmp(ntmp)
        real(kind=dp) :: xmid
        
    !-------------------------------- finite-volume arrays !

    !   Arrays represent a "block" of finite-volume tracers
    !   to remap. The 1st dim. is the no. of DoF per cell,
    !   NDOF=1 is a standard finite-volume scheme where the
    !   data is specified as cell means. NDOF>1 is reserved
    !   for future use with DG-style schemes. NVAR is the
    !   number of tracers to remap. Processing tracers in a
    !   batch is typically more efficient than one-by-one. 
    !   The last dim. is the no. cells (layers) in the grid.

        real(kind=dp) :: init(ndof,nvar,npos-1)
        real(kind=dp) :: fdat(ndof,nvar,npos-1)        
        real(kind=dp) :: ftmp(ndof,nvar,ntmp-1)
        
    !------------------------------ method data-structures !
        type(rmap_work) :: work
        type(rmap_opts) :: opts
        type(rcon_ends) :: bc_l(nvar)
        type(rcon_ends) :: bc_r(nvar)

    !------------------------------ define a simple domain !

        call linspace(0.d0,1.d0,npos,xpos)
        call rndspace(0.d0,1.d0,ntmp,xtmp)

    !------------------------------ setup some simple data !

        do ipos = +1, npos-1

            xmid = xpos(ipos+0) * 0.5d+0 &
        &        + xpos(ipos+1) * 0.5d+0

            if (xmid .lt. 0.075d0) then
            init(1,1,ipos) = +1.0d0
            else &
            if (xmid .lt. 0.80d0) then
            init(1,1,ipos) = +2.0d0
            else
            init(1,1,ipos) = -0.5d0 * xmid ** 2
            end if

        end do

    !------------------------------ specify method options !

        opts%edge_meth = p5e_method     ! 5th-order edge interp.
        opts%cell_meth = pqm_method     ! PQM method in cells
        opts%cell_lims = mono_limit     ! monotone limiter
       !opts%wall_lims = weno_limit
        
    !------------------------------ set BC.'s at endpoints !

        bc_l%bcopt = bcon_loose         ! "loose" = extrapolate
        bc_r%bcopt = bcon_loose

    !------------------------------ init. method workspace !

        call work%init(npos,nvar,opts)

    !------------------------------ re-map back-and-forth: !

        fdat = init

        do ipos = +1, +1000

    !------------------------------ re-map from dat-to-tmp !

        call rmap1d(npos,ntmp,nvar,ndof, &
        &           xpos,xtmp,fdat,ftmp, &
        &           bc_l,bc_r,work,opts)

    !------------------------------ re-map from tmp-to-dat !

        call rmap1d(ntmp,npos,nvar,ndof, &
        &           xtmp,xpos,ftmp,fdat, &
        &           bc_l,bc_r,work,opts)

        end do

    !------------------------------ clear method workspace !

        call work%free()

    !------------------------------ dump results to stdout !

        print*,"Cell data: [INIT] [RMAP] "

        do ipos = +1, npos-1

            print *, init(1,1,ipos) &
        &          , fdat(1,1,ipos)

        end do

        print*,"Conservation defect := " &
        &     , sum(init) - sum(fdat)

    end program



