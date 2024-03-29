
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
    ! PLM.f90: a 1d, slope-limited piecewise linear method.
    !
    ! Darren Engwirda 
    ! 25-Oct-2021
    ! d [dot] engwirda [at] gmail [dot] com
    !
    !

    pure subroutine plm(npos,nvar,ndof,delx, &
        &               fdat,fhat,ilim)

    !
    ! NPOS  no. edges over grid.
    ! NVAR  no. state variables.
    ! NDOF  no. degrees-of-freedom per grid-cell .
    ! DELX  grid-cell spacing array. LENGTH(DELX) == +1 if 
    !       spacing is uniform .
    ! FDAT  grid-cell moments array. FDAT is an array with
    !       SIZE = NDOF-by-NVAR-by-NPOS-1 .
    ! FHAT  grid-cell re-con. array. FHAT is an array with
    !       SIZE = MDOF-by-NVAR-by-NPOS-1 .
    ! ILIM  cell slope-limiting selection .
    !

        implicit none

    !------------------------------------------- arguments !
        integer      , intent( in) :: npos,nvar
        integer      , intent( in) :: ndof,ilim
        real(kind=dp), intent( in) :: delx(:)
        real(kind=dp), intent(out) :: fhat(:,:,:)
        real(kind=dp), intent( in) :: fdat(:,:,:)

        if (size(delx).gt.+1) then
        
    !------------------------------- variable grid-spacing !
        
            call plmv(npos,nvar,ndof,fdat,fhat,ilim, &
        &             delx)
        
        else
        
    !------------------------------- constant grid-spacing !
        
            call plmc(npos,nvar,ndof,fdat,fhat,ilim)
        
        end if

        return

    end  subroutine
    
    !------------------------- assemble PLM reconstruction !

    pure subroutine plmv(npos,nvar,ndof,fdat,fhat,ilim, &
        &                delx)

    !
    ! *this is the variable grid-spacing variant .
    !
    ! NPOS  no. edges over grid.
    ! NVAR  no. state variables.
    ! NDOF  no. degrees-of-freedom per grid-cell .
    ! DELX  grid-cell spacing array. LENGTH(DELX) == +1 if 
    !       spacing is uniform .
    ! FDAT  grid-cell moments array. FDAT is an array with
    !       SIZE = NDOF-by-NVAR-by-NPOS-1 .
    ! FHAT  grid-cell re-con. array. FHAT is an array with
    !       SIZE = MDOF-by-NVAR-by-NPOS-1 .
    ! ILIM  cell slope-limiting selection .
    !

        implicit none

    !------------------------------------------- arguments !
        integer      , intent( in) :: npos,nvar
        integer      , intent( in) :: ndof,ilim
        real(kind=dp), intent( in) :: delx(:)
        real(kind=dp), intent(out) :: fhat(:,:,:)
        real(kind=dp), intent( in) :: fdat(:,:,:)
        
    !------------------------------------------- variables !
        integer                    :: ipos,ivar
        integer                    :: head,tail
        real(kind=dp)              :: dfds(-1:+1)

        head = +1; tail = npos - 1

        if (npos.eq.2) then
    !----------------------- reduce order if small stencil !
        do  ivar = +1, nvar
            fhat(1,ivar,1) = &
        &   fdat(1,ivar,1)
            fhat(2,ivar,1) = 0.d+0
        end do
        end if

        if (ndof.le.0) return
        if (npos.le.2) return
  
    !-------------------------------------- lower-endpoint !
        
        do  ivar = +1 , nvar-0

            call plsv(dfds,ilim , &
        &   fdat(1,ivar,head+0) , &
        &   delx(head+0), &
        &   fdat(1,ivar,head+0) , &
        &   delx(head+0), &
        &   fdat(1,ivar,head+1) , &
        &   delx(head+1))

            fhat(1,ivar,head) = &
        &   fdat(1,ivar,head)
            fhat(2,ivar,head) = dfds(0)

        end do
        
    !-------------------------------------- upper-endpoint !
        
        do  ivar = +1 , nvar-0

            call plsv(dfds,ilim , &
        &   fdat(1,ivar,tail-1) , &
        &   delx(tail-1), &
        &   fdat(1,ivar,tail+0) , &
        &   delx(tail+0), &
        &   fdat(1,ivar,tail+0) , &
        &   delx(tail+0))

            fhat(1,ivar,tail) = &
        &   fdat(1,ivar,tail)
            fhat(2,ivar,tail) = dfds(0)

        end do

    !-------------------------------------- interior cells !

        do  ipos = +2 , npos-2
        do  ivar = +1 , nvar-0

            call plsv(dfds,ilim , &
        &   fdat(1,ivar,ipos-1) , &
        &   delx(ipos-1), &
        &   fdat(1,ivar,ipos+0) , &
        &   delx(ipos+0), &
        &   fdat(1,ivar,ipos+1) , &
        &   delx(ipos+1))

            fhat(1,ivar,ipos) = &
        &   fdat(1,ivar,ipos)
            fhat(2,ivar,ipos) = dfds(0)

        end do
        end do
        
        return
        
    end  subroutine

    !------------------------- assemble PLM reconstruction !

    pure subroutine plmc(npos,nvar,ndof,fdat,fhat,ilim)

    !
    ! *this is the constant grid-spacing variant .
    !
    ! NPOS  no. edges over grid.
    ! NVAR  no. state variables.
    ! NDOF  no. degrees-of-freedom per grid-cell .
    ! DELX  grid-cell spacing array. LENGTH(DELX) == +1 if 
    !       spacing is uniform .
    ! FDAT  grid-cell moments array. FDAT is an array with
    !       SIZE = NDOF-by-NVAR-by-NPOS-1 .
    ! FHAT  grid-cell re-con. array. FHAT is an array with
    !       SIZE = MDOF-by-NVAR-by-NPOS-1 .
    ! ILIM  cell slope-limiting selection .
    !

        implicit none

    !------------------------------------------- arguments !
        integer      , intent( in) :: npos,nvar
        integer      , intent( in) :: ndof,ilim
        real(kind=dp), intent(out) :: fhat(:,:,:)
        real(kind=dp), intent( in) :: fdat(:,:,:)
        
    !------------------------------------------- variables !
        integer                    :: ipos,ivar
        integer                    :: head,tail
        real(kind=dp)              :: dfds(-1:+1)

        head = +1; tail = npos - 1

        if (npos.eq.2) then
    !----------------------- reduce order if small stencil !
        do  ivar = +1, nvar
            fhat(1,ivar,1) = &
        &   fdat(1,ivar,1)
            fhat(2,ivar,1) = 0.d+0
        end do
        end if

        if (ndof.le.0) return
        if (npos.le.2) return
  
    !-------------------------------------- lower-endpoint !
        
        do  ivar = +1 , nvar-0

            call plsc(dfds,ilim , &
        &   fdat(1,ivar,head+0) , &
        &   fdat(1,ivar,head+0) , &
        &   fdat(1,ivar,head+1))

            fhat(1,ivar,head) = &
        &   fdat(1,ivar,head)
            fhat(2,ivar,head) = dfds(0)

        end do
        
    !-------------------------------------- upper-endpoint !
        
        do  ivar = +1 , nvar-0

            call plsc(dfds,ilim , &
        &   fdat(1,ivar,tail-1) , &
        &   fdat(1,ivar,tail+0) , &
        &   fdat(1,ivar,tail+0))

            fhat(1,ivar,tail) = &
        &   fdat(1,ivar,tail)
            fhat(2,ivar,tail) = dfds(0)

        end do

    !-------------------------------------- interior cells !

        do  ipos = +2 , npos-2
        do  ivar = +1 , nvar-0

            call plsc(dfds,ilim , &
        &   fdat(1,ivar,ipos-1) , &
        &   fdat(1,ivar,ipos+0) , &
        &   fdat(1,ivar,ipos+1))

            fhat(1,ivar,ipos) = &
        &   fdat(1,ivar,ipos)
            fhat(2,ivar,ipos) = dfds(0)

        end do
        end do
        
        return
        
    end  subroutine 
    
    !------------------------------- assemble PLM "slopes" !
    
    pure subroutine plsv(dfds,ilim,ffll,hhll, &
        &                ff00,hh00,ffrr,hhrr)

    !
    ! *this is the variable grid-spacing variant .
    !
    ! DFDS  piecewise linear gradients in local co-ord.'s.
    !       DFDS(+0) is a centred, slope-limited estimate,
    !       DFDS(-1), DFDS(+1) are left- and right-biased
    !       estimates (un-limited).
    ! FFLL  left -biased grid-cell mean.
    ! HHLL  left -biased grid-cell spac.
    ! FF00  centred grid-cell mean.
    ! HH00  centred grid-cell spac.
    ! FFRR  right-biased grid-cell mean.
    ! HHRR  right-biased grid-cell spac.
    !

        implicit none

    !------------------------------------------- arguments !
        integer      , intent( in) :: ilim
        real(kind=dp), intent( in) :: ffll,ff00,ffrr
        real(kind=dp), intent( in) :: hhll,hh00,hhrr
        real(kind=dp), intent(out) :: dfds(-1:+1)

    !------------------------------------------- variables !
        real(kind=dp)  :: fell,ferr,scal

        real(kind=dp), parameter :: ZERO = 1.d-14

    !---------------------------- 2nd-order approximations !

        if (ilim .eq. null_limit) then

    !---------------------------- calc. centred derivative !

            fell = (hh00*ffll+hhll*ff00) & 
        &        / (hhll+hh00)        
            ferr = (hhrr*ff00+hh00*ffrr) &
        &        / (hh00+hhrr)

            dfds(-1) = (ff00 - ffll) &
        &       / (hhll + hh00) * hh00
            dfds(+1) = (ffrr - ff00) &
        &       / (hh00 + hhrr) * hh00

            dfds(+0) = &
        &       0.5d+0 * (ferr - fell)

            return

        end if

    !---------------------------- calc. limited PLM slopes !

        dfds(-1) = ff00-ffll
        dfds(+1) = ffrr-ff00

        if (dfds(-1) * &
        &   dfds(+1) .gt. 0.0d+0) then

    !---------------------------- calc. ll//rr edge values !

            fell = (hh00*ffll+hhll*ff00) & 
        &        / (hhll+hh00)        
            ferr = (hhrr*ff00+hh00*ffrr) &
        &        / (hh00+hhrr)

    !---------------------------- calc. centred derivative !
            
            dfds(+0) = &
        &       0.5d+0 * (ferr - fell)

    !---------------------------- monotonic slope-limiting !
            
            scal = min(abs(dfds(-1)), &
        &              abs(dfds(+1))) &
        &        / max(abs(dfds(+0)), &
                       ZERO)
            scal = min(scal,+1.0d+0)

            dfds(+0) = scal * dfds(+0)

        else

    !---------------------------- flatten if local extrema ! 
      
            dfds(+0) =      +0.0d+0
        
        end if
        
    !---------------------------- scale onto local co-ord. !
        
        dfds(-1) = dfds(-1) &
        &       / (hhll + hh00) * hh00
        dfds(+1) = dfds(+1) &
        &       / (hh00 + hhrr) * hh00

        return

    end  subroutine
    
    !------------------------------- assemble PLM "slopes" !
    
    pure subroutine plsc(dfds,ilim,ffll,ff00,ffrr)

    !
    ! *this is the constant grid-spacing variant .
    !
    ! DFDS  piecewise linear gradients in local co-ord.'s.
    !       DFDS(+0) is a centred, slope-limited estimate,
    !       DFDS(-1), DFDS(+1) are left- and right-biased
    !       estimates (un-limited).
    ! FFLL  left -biased grid-cell mean.
    ! FF00  centred grid-cell mean.
    ! FFRR  right-biased grid-cell mean.
    !

        implicit none

    !------------------------------------------- arguments !
        integer      , intent( in) :: ilim
        real(kind=dp), intent( in) :: ffll,ff00,ffrr
        real(kind=dp), intent(out) :: dfds(-1:+1)

    !------------------------------------------- variables !
        real(kind=dp)  :: fell,ferr,scal

        real(kind=dp), parameter :: ZERO = 1.d-14

    !---------------------------- 2nd-order approximations !

        if (ilim .eq. null_limit) then

    !---------------------------- calc. centred derivative !

            fell = (ffll+ff00) * .5d+0        
            ferr = (ff00+ffrr) * .5d+0

            dfds(-1) = &
        &       0.5d+0 * (ff00 - ffll)
            dfds(+1) = &
        &       0.5d+0 * (ffrr - ff00)

            dfds(+0) = &
        &       0.5d+0 * (ferr - fell)

            return

        end if

    !---------------------------- calc. limited PLM slopes !

        dfds(-1) = ff00-ffll
        dfds(+1) = ffrr-ff00

        if (dfds(-1) * &
        &   dfds(+1) .gt. 0.0d+0) then

    !---------------------------- calc. ll//rr edge values !

            fell = (ffll+ff00) * .5d+0        
            ferr = (ff00+ffrr) * .5d+0

    !---------------------------- calc. centred derivative !
            
            dfds(+0) = &
        &       0.5d+0 * (ferr - fell)

    !---------------------------- monotonic slope-limiting !
            
            scal = min(abs(dfds(-1)), &
        &              abs(dfds(+1))) &
        &        / max(abs(dfds(+0)), &
                       ZERO)
            scal = min(scal,+1.0d+0)

            dfds(+0) = scal * dfds(+0)

        else

    !---------------------------- flatten if local extrema ! 
      
            dfds(+0) =      +0.0d+0
        
        end if
        
    !---------------------------- scale onto local co-ord. !
        
        dfds(-1) = + 0.5d+0 * dfds(-1) 
        dfds(+1) = + 0.5d+0 * dfds(+1)

        return

    end  subroutine   
    
    
    
