
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
    ! PQM.f90: a 1d slope-limited, piecewise quartic recon.
    !
    ! Darren Engwirda 
    ! 06-Nov-2021
    ! d [dot] engwirda [at] gmail [dot] com
    !
    !

    ! White, L. and Adcroft, A., A high-order finite volume
    ! remapping scheme for nonuniform grids: The piecewise 
    ! quartic method (PQM), J. Comp. Phys., 227 (15), 2008, 
    ! 7394-7422, https://doi.org/10.1016/j.jcp.2008.04.026.
    !

    pure subroutine pqm(npos,nvar,ndof,delx, &
        &               fdat,fhat,edge,dfdx, &
        &               oscl,ilim,wlim,halo)

    !
    ! NPOS  no. edges over grid.
    ! NVAR  no. state variables.
    ! NDOF  no. degrees-of-freedom per grid-cell.
    ! DELX  grid-cell spacing array. LENGTH(DELX) == +1 if 
    !       spacing is uniform .
    ! FDAT  grid-cell moments array. FDAT is an array with
    !       SIZE = NDOF-by-NVAR-by-NPOS-1 .
    ! FHAT  grid-cell re-con. array. FHAT is an array with
    !       SIZE = MDOF-by-NVAR-by-NPOS-1 .
    ! EDGE  edge-centred interp. for function-value. EDGE
    !       is an array with SIZE = NVAR-by-NPOS .
    ! DFDX  edge-centred interp. for 1st-derivative. DFDX
    !       is an array with SIZE = NVAR-by-NPOS .
    ! OSCL  grid-cell oscil. dof.'s. OSCL is an array with
    !       SIZE = +2  -by-NVAR-by-NPOS-1 .
    ! ILIM  cell slope-limiting selection .
    ! WLIM  wall slope-limiting selection .
    ! HALO  width of re-con. stencil, symmetric about mid. .
    !

        implicit none
        
    !------------------------------------------- arguments !
        integer      , intent(in)  :: npos,nvar,ndof
        integer      , intent(in)  :: ilim,wlim,halo
        real(kind=dp), intent(out) :: fhat(:,:,:)
        real(kind=dp), intent(in)  :: oscl(:,:,:)
        real(kind=dp), intent(in)  :: delx(:)
        real(kind=dp), intent(in)  :: fdat(:,:,:)
        real(kind=dp), intent(in)  :: edge(:,:)
        real(kind=dp), intent(in)  :: dfdx(:,:)

    !------------------------------------------- variables !
        integer       :: ipos,ivar,iill,iirr,head,tail
        real(kind=dp) :: ff00,ffll,ffrr,hh00,hhll,hhrr
        real(kind=dp) :: xhat
        integer       :: mono
        real(kind=dp) :: fell,ferr
        real(kind=dp) :: dell,derr
        real(kind=dp) :: dfds(-1:+1)
        real(kind=dp) :: uhat(+1:+5)
        real(kind=dp) :: lhat(+1:+5)
        real(kind=dp) :: wval(+1:+2)
    
        head = +1; tail = npos - 1

        if (npos.le.2) then
    !----- default to reduced order if insufficient points !
        do  ivar = +1, nvar
            fhat(1,ivar,+1) = fdat(1,ivar,+1)
            fhat(2,ivar,+1) = 0.d0
            fhat(3,ivar,+1) = 0.d0
            fhat(4,ivar,+1) = 0.d0
            fhat(5,ivar,+1) = 0.d0
        end do
        end if

        if (ndof.le.0) return
        if (npos.le.2) return

    !------------------- reconstruct function on each cell !

        do  ipos = +1 , npos-1

        iill = max(head,ipos-1)
        iirr = min(tail,ipos+1)

        do  ivar = +1 , nvar-0

    !----------------------------- cell mean + edge values !
            
            ff00 = fdat(1,ivar,ipos)
            ffll = fdat(1,ivar,iill)
            ffrr = fdat(1,ivar,iirr)
           
            fell = edge(ivar,ipos+0)
            ferr = edge(ivar,ipos+1)
        
    !----------------------------- calc. LL/00/RR gradient !
            
            if (size(delx).gt.+1) then

            hh00 = delx(ipos)
            hhll = delx(iill)
            hhrr = delx(iirr)

            xhat = delx(ipos+0)*.5d+0

            call plsv (dfds,mono_limit, &
    &                  ffll,hhll,ff00 , &
    &                  hh00,ffrr,hhrr)
            else
            
            xhat = delx(    +1)*.5d+0
            
            call plsc (dfds,mono_limit, &
    &                  ffll,ff00,ffrr)

            end if              

            dell = dfdx (ivar,ipos+0)
            dell = dell * xhat
            
            derr = dfdx (ivar,ipos+1)
            derr = derr * xhat

    !----------------------------- calc. cell-wise profile !
            
            select case(ilim)
            case (null_limit)

    !----------------------------- calc. unlimited profile !              
            
            call pqmfn(ff00,ffll,ffrr, &
    &                  fell,ferr,dell, &
    &                  derr,dfds,uhat, &
    &                  lhat,mono)
        
    !----------------------------- pref. unlimited profile !
            
            wval(1) = +1.d+0
            wval(2) = +0.d+0
            
            case (mono_limit)
            
    !----------------------------- calc. monotonic profile !              
            
            call pqmfn(ff00,ffll,ffrr, &
    &                  fell,ferr,dell, &
    &                  derr,dfds,uhat, &
    &                  lhat,mono)
            
    !----------------------------- pref. monotonic profile !
            
            wval(1) = +0.d+0
            wval(2) = +1.d+0
            
            case (weno_limit)

    !----------------------------- calc. monotonic profile !              
            
            call pqmfn(ff00,ffll,ffrr, &
    &                  fell,ferr,dell, &
    &                  derr,dfds,uhat, &
    &                  lhat,mono)
            
            if (mono.gt.+0) then
            
    !----------------------------- calc. WENO-type weights !      
            
            call wenoi(npos,delx,oscl, &
    &                  ipos,ivar,halo, &
    &                  wlim,wval)
            
            else
            
    !----------------------------- pref. unlimited profile !
            
            wval(1) = +1.d+0
            wval(2) = +0.d+0
            
            end if
            
            end select

    !----------------------------- blend "null" and "mono" !
                   
            fhat(1,ivar,ipos) = &
    &       wval(1) * uhat(1) + &
    &       wval(2) * lhat(1)
            fhat(2,ivar,ipos) = &
    &       wval(1) * uhat(2) + &
    &       wval(2) * lhat(2)
            fhat(3,ivar,ipos) = &
    &       wval(1) * uhat(3) + &
    &       wval(2) * lhat(3)
            fhat(4,ivar,ipos) = &
    &       wval(1) * uhat(4) + &
    &       wval(2) * lhat(4)
            fhat(5,ivar,ipos) = &
    &       wval(1) * uhat(5) + &
    &       wval(2) * lhat(5)

        end do
        
        end do
        
        return

    end subroutine
    
    !----------- assemble piecewise quartic reconstruction !
    
    pure subroutine pqmfn(ff00,ffll,ffrr,fell, &
        &                 ferr,dell,derr,dfds, &
        &                 uhat,lhat,mono)

    !
    ! FF00  centred grid-cell mean.
    ! FFLL  left -biased grid-cell mean.
    ! FFRR  right-biased grid-cell mean.
    ! FELL  left -biased edge interp.
    ! FERR  right-biased edge interp.
    ! DELL  left -biased edge df//dx.
    ! DERR  right-biased edge df//dx.
    ! DFDS  piecewise linear gradients in local co-ord.'s.
    !       DFDS(+0) is a centred, slope-limited estimate,
    !       DFDS(-1), DFDS(+1) are left- and right-biased
    !       estimates (unlimited). 
    ! UHAT  unlimited PPM reconstruction coefficients .
    ! LHAT  monotonic PPM reconstruction coefficients .
    ! MONO  slope-limiting indicator, MONO > +0 if some 
    !       limiting has occured .
    !

        implicit none

    !------------------------------------------- arguments !
        real(kind=dp), intent(in)    :: ff00
        real(kind=dp), intent(in)    :: ffll,ffrr
        real(kind=dp), intent(inout) :: fell,ferr
        real(kind=dp), intent(inout) :: dell,derr
        real(kind=dp), intent(in)    :: dfds(-1:+1)
        real(kind=dp), intent(out)   :: uhat(+1:+5)
        real(kind=dp), intent(out)   :: lhat(+1:+5)
        integer      , intent(out)   :: mono
          
    !------------------------------------------- variables !
        integer         :: turn
        real(kind=dp)   :: grad
        real(kind=dp)   :: iflx(+1:+2)
        real(kind=dp)   :: junk(+1:+3)
        logical         :: haveroot
        
    !-------------------------------- "null" slope-limiter !
        
        mono    = 0
        
        uhat(1) = &
    & + (30.d+0 / 16.d+0) * ff00 &
    & - ( 7.d+0 / 16.d+0) *(ferr+fell) &
    & + ( 1.d+0 / 16.d+0) *(derr-dell)
        uhat(2) = &
    & + ( 3.d+0 /  4.d+0) *(ferr-fell) &
    & - ( 1.d+0 /  4.d+0) *(derr+dell)
        uhat(3) = &
    & - (30.d+0 /  8.d+0) * ff00 &
    & + (15.d+0 /  8.d+0) *(ferr+fell) &
    & - ( 3.d+0 /  8.d+0) *(derr-dell)
        uhat(4) = &
    & - ( 1.d+0 /  4.d+0) *(ferr-fell  &
    &                      -derr-dell)
        uhat(5) = &
    & + (30.d+0 / 16.d+0) * ff00 &
    & - (15.d+0 / 16.d+0) *(ferr+fell) &
    & + ( 5.d+0 / 16.d+0) *(derr-dell)
  
    !-------------------------------- "mono" slope-limiter !
    
        if((ffrr - ff00) * & 
    &      (ff00 - ffll) .le. 0.d+0) then

    !----------------------------------- "flatten" extrema !
    
            mono = +1

            lhat(1) = ff00
            lhat(2) = 0.d0
            lhat(3) = 0.d0
            lhat(4) = 0.d0
            lhat(5) = 0.d0
             
            return
              
        end if

    !----------------------------------- limit edge slopes !

        if (dell * dfds(-1) .lt. 0.d+0) then

            lhat(:) =  0.0d+0
            call ppmfn(ff00,ffll,ffrr,fell,&
    &                  ferr,dfds,junk,lhat,&
    &                  mono)
            mono = +1; return

        end if

        if (derr * dfds(+1) .lt. 0.d+0) then

            lhat(:) =  0.0d+0
            call ppmfn(ff00,ffll,ffrr,fell,&
    &                  ferr,dfds,junk,lhat,&
    &                  mono)
            mono = +1; return

        end if

    !----------------------------------- limit edge values !

        if((ffll - fell) * &
    &      (fell - ff00) .le. 0.d+0) then

            mono = +1
            fell = ff00 - dfds(0)
            
        end if

        if((ffrr - ferr) * &
    &      (ferr - ff00) .le. 0.d+0) then

            mono = +1
            ferr = ff00 + dfds(0)
            
        end if
    
    !----------------------------------- limit cell values !
    
        lhat(1) = &
    & + (30.d+0 / 16.d+0) * ff00 &
    & - ( 7.d+0 / 16.d+0) *(ferr+fell) &
    & + ( 1.d+0 / 16.d+0) *(derr-dell)
        lhat(2) = &
    & + ( 3.d+0 /  4.d+0) *(ferr-fell) &
    & - ( 1.d+0 /  4.d+0) *(derr+dell)
        lhat(3) = &
    & - (30.d+0 /  8.d+0) * ff00 &
    & + (15.d+0 /  8.d+0) *(ferr+fell) &
    & - ( 3.d+0 /  8.d+0) *(derr-dell)
        lhat(4) = &
    & - ( 1.d+0 /  4.d+0) *(ferr-fell  &
    &                      -derr-dell)
        lhat(5) = &
    & + (30.d+0 / 16.d+0) * ff00 &
    & - (15.d+0 / 16.d+0) *(ferr+fell) &
    & + ( 5.d+0 / 16.d+0) *(derr-dell)

    !------------------ calc. inflexion via 2nd-derivative !
        
        call roots_2(12.d+0 * lhat(5), &
    &                 6.d+0 * lhat(4), &
    &                 2.d+0 * lhat(3), &
    &                 iflx  , haveroot )

        if (haveroot) then

        turn = +0

        if ( ( iflx(1) .ge. -1.d+0 ) &
    &  .and. ( iflx(1) .le. +1.d+0 ) ) then

    !------------------ check for non-monotonic inflection !

        grad = lhat(2) &
    &+ (iflx(1)**1) * 2.d+0* lhat(3) &
    &+ (iflx(1)**2) * 3.d+0* lhat(4) &
    &+ (iflx(1)**3) * 4.d+0* lhat(5)

        if (grad * dfds(0) .le. 0.d+0) then

            ferr = ffrr
            fell = ffll

            lhat(:) =  0.0d+0
            call ppmfn(ff00,ffll,ffrr,fell,&
    &                  ferr,dfds,junk,lhat,&
    &                  mono)
            mono = +2; return



            if (abs(dfds(-1)) &
    &      .lt. abs(dfds(+1)) ) then

                turn = -1   ! modify L

            else

                turn = +1   ! modify R

            end if

        end if

        end if
            
        if ( ( iflx(2) .ge. -1.d+0 ) &
    &  .and. ( iflx(2) .le. +1.d+0 ) ) then

    !------------------ check for non-monotonic inflection !
                
        grad = lhat(2) &
    &+ (iflx(2)**1) * 2.d+0* lhat(3) &
    &+ (iflx(2)**2) * 3.d+0* lhat(4) &
    &+ (iflx(2)**3) * 4.d+0* lhat(5)

        if (grad * dfds(0) .le. 0.d+0) then

            ferr = ffrr
            fell = ffll

            lhat(:) =  0.0d+0
            call ppmfn(ff00,ffll,ffrr,fell,&
    &                  ferr,dfds,junk,lhat,&
    &                  mono)
            mono = +2; return




            if (abs(dfds(-1)) &
    &      .lt. abs(dfds(+1)) ) then

                turn = -1   ! modify L

            else

                turn = +1   ! modify R

            end if

        end if

        end if
            
    !------------------ pop non-monotone inflexion to edge !
              
        if (turn .eq. -1) then

    !------------------ pop inflection points onto -1 edge !
    
            mono = +2

            ferr = ffrr
            fell = ffll
            derr = &
    &    - ( 5.d+0 /  1.d+0) * ff00 &
    &    + ( 3.d+0 /  1.d+0) * ferr &
    &    + ( 2.d+0 /  1.d+0) * fell
            dell = &
    &    + ( 5.d+0 /  3.d+0) * ff00 &
    &    - ( 1.d+0 /  3.d+0) * ferr &
    &    - ( 4.d+0 /  3.d+0) * fell

            if (dell*dfds(-1) .lt. 0.d+0) then

            lhat(:) =  0.0d+0
            call ppmfn(ff00,ffll,ffrr,fell,&
    &                  ferr,dfds,junk,lhat,&
    &                  mono)
            mono = +2; return

            else &
    &       if (derr*dfds(+1) .lt. 0.d+0) then

            lhat(:) =  0.0d+0
            call ppmfn(ff00,ffll,ffrr,fell,&
    &                  ferr,dfds,junk,lhat,&
    &                  mono)
            mono = +2; return

            end if

            lhat(1) = &
    &    + (30.d+0 / 16.d+0) * ff00 &
    &    - ( 7.d+0 / 16.d+0) *(ferr+fell) &
    &    + ( 1.d+0 / 16.d+0) *(derr-dell)
            lhat(2) = &
    &    + ( 3.d+0 /  4.d+0) *(ferr-fell) &
    &    - ( 1.d+0 /  4.d+0) *(derr+dell)
            lhat(3) = &
    &    - (30.d+0 /  8.d+0) * ff00 &
    &    + (15.d+0 /  8.d+0) *(ferr+fell) &
    &    - ( 3.d+0 /  8.d+0) *(derr-dell)
            lhat(4) = &
    &    - ( 1.d+0 /  4.d+0) *(ferr-fell  &
    &                         -derr-dell)
            lhat(5) = &
    &    + (30.d+0 / 16.d+0) * ff00 &
    &    - (15.d+0 / 16.d+0) *(ferr+fell) &
    &    + ( 5.d+0 / 16.d+0) *(derr-dell)

        end if

        if (turn .eq. +1) then

    !------------------ pop inflection points onto +1 edge !

            mono = +2

            ferr = ffrr
            fell = ffll
            derr = &
    &    - ( 5.d+0 /  3.d+0) * ff00 &
    &    + ( 4.d+0 /  3.d+0) * ferr &
    &    + ( 1.d+0 /  3.d+0) * fell
            dell = &
    &    + ( 5.d+0 /  1.d+0) * ff00 &
    &    - ( 2.d+0 /  1.d+0) * ferr &
    &    - ( 3.d+0 /  1.d+0) * fell

            if (dell*dfds(-1) .lt. 0.d+0) then

            lhat(:) =  0.0d+0
            call ppmfn(ff00,ffll,ffrr,fell,&
    &                  ferr,dfds,junk,lhat,&
    &                  mono)
            mono = +2; return

            else & 
    &       if (derr*dfds(+1) .lt. 0.d+0) then

            lhat(:) =  0.0d+0
            call ppmfn(ff00,ffll,ffrr,fell,&
    &                  ferr,dfds,junk,lhat,&
    &                  mono)
            mono = +2; return

            end if

            lhat(1) = &
    &    + (30.d+0 / 16.d+0) * ff00 &
    &    - ( 7.d+0 / 16.d+0) *(ferr+fell) &
    &    + ( 1.d+0 / 16.d+0) *(derr-dell)
            lhat(2) = &
    &    + ( 3.d+0 /  4.d+0) *(ferr-fell) &
    &    - ( 1.d+0 /  4.d+0) *(derr+dell)
            lhat(3) = &
    &    - (30.d+0 /  8.d+0) * ff00 &
    &    + (15.d+0 /  8.d+0) *(ferr+fell) &
    &    - ( 3.d+0 /  8.d+0) *(derr-dell)
            lhat(4) = &
    &    - ( 1.d+0 /  4.d+0) *(ferr-fell  &
    &                         -derr-dell)
            lhat(5) = &
    &    + (30.d+0 / 16.d+0) * ff00 &
    &    - (15.d+0 / 16.d+0) *(ferr+fell) &
    &    + ( 5.d+0 / 16.d+0) *(derr-dell)

        end if
    
        end if        ! haveroot

        return

    end subroutine
    

    
