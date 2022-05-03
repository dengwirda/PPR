
!   gfortran -pedantic -cpp -O3 -flto ex_5.f90 -o ex_5
!   ./ex_5

!   Assemble high-order interpolants over a uniform domain.
!

#   include "../src/ppr_1d.f90"

    program ex

        use ppr_1d

        implicit none

        integer, parameter :: npos = 37 ! no. edge 
        integer, parameter :: nvar = 1  ! no. variables to build
        integer, parameter :: ndof = 1  ! no. FV DoF per cell
        integer :: ipos,jpos,mdof

    !------------------------------- domain discretisation !
        real(kind=dp) :: xpos(npos),xdel(npos-1)
        real(kind=dp) :: xmid,xhat,xloc,floc
        
    !-------------------------------- finite-volume arrays !

    !   Arrays represent a "block" of finite-volume tracers
    !   to remap. The 1st dim. is the no. of DoF per cell,
    !   NDOF=1 is a standard finite-volume scheme where the
    !   data is specified as cell means. NDOF>1 is reserved
    !   for future use with DG-style schemes. NVAR is the
    !   number of tracers to remap. Processing tracers in a
    !   batch is typically more efficient than one-by-one. 
    !   The last dim. is the no. cells (layers) in the grid.

        real(kind=dp) :: fdat(ndof,nvar,npos-1)
        
    !-------------------------------- reconstruction coeff !

    !   Coeff. for the piecewise polynomial reconstruction.
    !   A polynomial is assembled for each cell w.r.t. a
    !   "local" cell coordinate system: each cell is mapped
    !   onto [-1,+1]. The interpolants can be evaluated by
    !   taking the product FHAT*BVEC, where BVEC is a basis
    !   vector assembled at the interpolation points. Basis
    !   vectors can eb assembled via calls to BFUN1D(). 

        real(kind=dp) :: fhat(   5,nvar,npos-1)
        real(kind=dp) :: bvec(   5)
        real(kind=dp) :: spos(   5)

    !------------------------------ method data-structures !
        type(rcon_work) :: work
        type(rcon_opts) :: opts
        type(rcon_ends) :: bc_l(nvar)
        type(rcon_ends) :: bc_r(nvar)

    !------------------------------ define a simple domain !

        call rndspace(0.d0,1.d0,npos,xpos)
        
        xdel(1) = (xpos(npos)&
    &           -  xpos(   1)) / (npos- 1)

    !------------------------------ setup some simple data !

        do ipos = +1, npos-1

            xdel(ipos) = xpos(ipos+1) - xpos(ipos)
            print*, xdel(ipos)

            xmid = xpos(ipos+0) * 0.5d+0 &
    &            + xpos(ipos+1) * 0.5d+0

            fdat(1,1,ipos) = &
    &   .8d+0 * exp( -75.0d+0 * (xmid - 0.275d+0) ** 2 ) &
    & + .9d+0 * exp(-100.0d+0 * (xmid - 0.500d+0) ** 2 ) &
    & + 1.d+0 * exp(-125.0d+0 * (xmid - 0.725d+0) ** 2 )


            if (xmid .lt. 0.075d0) then
            fdat(1,1,ipos) = +1.0d0
            else &
            if (xmid .lt. 0.80d0) then
            fdat(1,1,ipos) = +2.0d0
            else
            fdat(1,1,ipos) = -0.5d0 * xmid ** 2
            end if
            
        end do

    !------------------------------ specify method options !

        opts%edge_meth = p5e_method     ! 5th-order edge interp.
        opts%cell_meth = ppm_method     ! PPM method in cells
        opts%cell_lims = mono_limit     ! monotone limiter
        
    !------------------------------ set BC.'s at endpoints !

        bc_l%bcopt = bcon_loose         ! "loose" = extrapolate
        bc_r%bcopt = bcon_loose

    !------------------------------ init. method workspace !

        call work%init(npos,nvar,opts)

    !------------------------------ build cell polynomials !

        fhat = 0.d+0; bvec = 0.d+0; spos = 0.d+0

        mdof = ndof1d (opts%cell_meth)

        call rcon1d(npos,nvar,ndof,xdel, &
    &               fdat,bc_l,bc_r,fhat, &
    &               work,opts)

    !------------------------------ clear method workspace !

        call work%free()

    !------------------------------ dump results to stdout !

        print*,"Eval. PPR interpolant: "

        spos(1) = -1.0d+0               ! eval. at local points
        spos(2) = -0.5d+0        
        spos(3) = +0.0d+0
        spos(4) = +0.5d+0
        spos(5) = +1.0d+0

        do ipos = +1, npos-1
        do jpos = +1, +5
    
            xmid = xpos(ipos+1)* 0.5d+0 &
    &            + xpos(ipos+0)* 0.5d+0

            xhat = xpos(ipos+1)* 0.5d+0 &
    &            - xpos(ipos+0)* 0.5d+0

            xloc = xmid + spos(jpos)*xhat
            
            call bfun1d(0,mdof,spos(jpos),bvec)

            floc = dot_product( &
    &       fhat(+1:mdof,1,ipos),bvec(+1:mdof))

            print *, xloc, floc

        end do
        end do

    end program



