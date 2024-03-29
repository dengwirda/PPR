
    !
    ! This program may be freely redistributed under the 
    ! condition that the copyright notices (including this 
    ! entire header) are not removed, and no compensation 
    ! is received through use of the software.  Private, 
    ! research, and institutional use is free.  You may 
    ! distribute modified versions of this code UNDER THE 
    ! CONDITION THAT THIS CODE AND ANY MODIFICATIONS MADE 
    ! TO IT IN THE SAME FILE REMAIN UNDER COPYRIGHT OF THE 
    ! ORIGINAL AUTHOR, BOTH SOURCE AND OBJECT CODE ARE 
    ! MADE FREELY AVAILABLE WITHOUT CHARGE, AND CLEAR 
    ! NOTICE IS GIVEN OF THE MODIFICATIONS.  Distribution 
    ! of this code as part of a commercial system is 
    ! permissible ONLY BY DIRECT ARRANGEMENT WITH THE 
    ! AUTHOR.  (If you are not directly supplying this 
    ! code to a customer, and you are instead telling them 
    ! how they can obtain it for free, then you are not 
    ! required to make any arrangement with me.) 
    !
    ! Disclaimer:  Neither I nor: Columbia University, the 
    ! National Aeronautics and Space Administration, nor 
    ! the Massachusetts Institute of Technology warrant 
    ! or certify this code in any way whatsoever.  This 
    ! code is provided "as-is" to be used at your own risk.
    !
    !

    !    
    ! P1E.f90: set edge estimates via degree-1 polynomials.
    !
    ! Darren Engwirda 
    ! 25-Oct-2021
    ! d [dot] engwirda [at] gmail [dot] com
    !
    !

    subroutine p1e(npos,nvar,ndof,delx, &
        &          fdat,bclo,bchi,edge, &
        &          dfdx,opts,dmin)

    !
    ! NPOS  no. edges over grid.
    ! NVAR  no. state variables.
    ! NDOF  no. degrees-of-freedom per grid-cell.
    ! DELX  grid-cell spacing array. LENGTH(DELX) == +1 if 
    !       spacing is uniform .
    ! FDAT  grid-cell moments array. FDAT is an array with
    !       SIZE = NDOF-by-NVAR-by-NPOS-1 .
    ! BCLO  boundary condition at lower endpoint.
    ! BCHI  boundary condition at upper endpoint.
    ! EDGE  edge-centred interp. for function-value. EDGE
    !       is an array with SIZE = NVAR-by-NPOS .
    ! DFDX  edge-centred interp. for 1st-derivative. DFDX
    !       is an array with SIZE = NVAR-by-NPOS .
    ! OPTS  method parameters. See RCON-OPTS for details .
    ! DMIN  min. grid-cell spacing thresh . 
    !

        implicit none

    !------------------------------------------- arguments !
        integer         , intent(in)    :: npos,nvar,ndof
        real(kind=dp)   , intent(in)    :: delx(:)
        real(kind=dp)   , intent(in)    :: fdat(:,:,:)
        type (rcon_ends), intent(in)    :: bclo(:)
        type (rcon_ends), intent(in)    :: bchi(:)
        real(kind=dp)   , intent(out)   :: edge(:,:)
        real(kind=dp)   , intent(out)   :: dfdx(:,:)
        real(kind=dp)   , intent(in)    :: dmin
        class(rcon_opts), intent(in)    :: opts

    !------------------------------------------- variables !
        integer                         :: ipos,ivar,bopt
        integer                         :: head,tail
        real(kind=dp)                   :: dd10
        real(kind=dp)                   :: delh(-1:+0)

        head = +2; tail = npos-1

        if (npos.lt.2) return
        if (npos.eq.2) then
    !----- default to reduced order if insufficient points !
        do  ivar = 1,nvar
        
            edge(ivar,1) = fdat(1,ivar,1)
            dfdx(ivar,1) = 0.d0
            
            edge(ivar,2) = fdat(1,ivar,1)
            dfdx(ivar,2) = 0.d0
            
        end do
        end if

        if (ndof.le.0) return
        if (npos.le.2) return

        if (opts%edge_meth &
    &  .lt.p1e_method) return
   
    ! Reconstruct edge-centred 2nd-order polynomials. Com- !
    ! pute values/slopes at edges directly. Full-order ex- !
    ! trapolation at endpoints.
    
        if (size(delx).eq.+1) then
        
            do  ipos = head , tail
    
    !--------------- reconstruction: constant grid-spacing !
            
            dd10 = delx(+1) * 2.d0
            
            do  ivar = +1, nvar

                edge(ivar,ipos) = &
        &         + delx(+1) * &
        &       fdat(1,ivar,ipos-1) &
        &         + delx(+1) * &
        &       fdat(1,ivar,ipos+0)

                dfdx(ivar,ipos) = &
        &         - 2.00d+00 * &
        &       fdat(1,ivar,ipos-1) &
        &         + 2.00d+00 * &
        &       fdat(1,ivar,ipos+0)

                edge(ivar,ipos) = &
        &       edge(ivar,ipos) / dd10
                dfdx(ivar,ipos) = &
        &       dfdx(ivar,ipos) / dd10

            end do
            
            end do
            
        else
        
            do  ipos = head , tail
    
    !--------------- reconstruction: variable grid-spacing !
            
            delh(-1) = &
        &       max(delx(ipos-1),dmin)
            delh(+0) = &
        &       max(delx(ipos+0),dmin)

            dd10 = delh(-1)+delh(+0)
            
            do  ivar = +1, nvar

                edge(ivar,ipos) = &
        &         + delh(+0) * &
        &       fdat(1,ivar,ipos-1) &
        &         + delh(-1) * &
        &       fdat(1,ivar,ipos+0)

                dfdx(ivar,ipos) = &
        &         - 2.00d+00 * &
        &       fdat(1,ivar,ipos-1) &
        &         + 2.00d+00 * &
        &       fdat(1,ivar,ipos+0)

                edge(ivar,ipos) = &
        &       edge(ivar,ipos) / dd10
                dfdx(ivar,ipos) = &
        &       dfdx(ivar,ipos) / dd10

            end do
            
            end do  
            
        end if
    
    !------------- 1st-order value/slope BC's at endpoints !

        do  ivar = +1, nvar

    !------------------------------------- at LHS endpoint !

            bopt = bclo(ivar)%bcopt
            dd10 = 1.d0

            if (bopt.eq.bcon_loose) then

                edge(ivar,head-1) = &
        &           fdat(+1,ivar,head-1)

                dfdx(ivar,head-1) = 0.d0

            else &
        &   if (bopt.eq.bcon_value) then

                edge(ivar,head-1) = &
        &           bclo(ivar)%value

                dfdx(ivar,head-1) = 0.d0

            else &
        &   if (bopt.eq.bcon_slope) then

                edge(ivar,head-1) = &
        &           fdat(+1,ivar,head-1)

                dfdx(ivar,head-1) = &
        &           bclo(ivar)%slope * dd10

            end if

    !------------------------------------- at RHS endpoint !
           
            bopt = bchi(ivar)%bcopt
            dd10 = 1.d0

            if (bopt.eq.bcon_loose) then

                edge(ivar,tail+1) = &
        &           fdat(+1,ivar,tail+0)
            
                dfdx(ivar,tail+1) = 0.d0

            else &
        &   if (bopt.eq.bcon_value) then

                edge(ivar,tail+1) = &
        &           bchi(ivar)%value
            
                dfdx(ivar,tail+1) = 0.d0

            else &
        &   if (bopt.eq.bcon_slope) then

                edge(ivar,tail+1) = &
        &           fdat(+1,ivar,tail+0)
            
                dfdx(ivar,tail+1) = &
        &           bchi(ivar)%slope * dd10

            end if

        end do
    
        return
        
    end subroutine
    
    
    
