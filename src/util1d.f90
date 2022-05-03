
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
    ! UTIL1D.f90: util. func. for 1-dim. grid manipulation.
    !
    ! Darren Engwirda 
    ! 08-Nov-2021
    ! d [dot] engwirda [at] gmail [dot] com
    !
    !

    function sum_kbn(xvec) result(ssum)

    ! 
    ! XVEC  array to sum.
    ! SSUM  sum(XVEC), via Kahan-Babuska-Neumaier sum .
    !

        implicit none

        real(kind=dp), intent(in)   :: xvec(:)

        integer         :: ipos
        real(kind=dp)   :: ssum,serr,stmp

    !------------------------------- compensated summation !

        ssum = 0.d0; serr = 0.d0

        do  ipos = +1, size(xvec) - 0
    
            stmp = ssum + xvec(ipos)

            if (abs(ssum) &
        &       .ge. abs(xvec(ipos))) then

            serr = &
        &   serr + ((ssum-stmp)+xvec(ipos))

            else

            serr = &
        &   serr + ((xvec(ipos)-stmp)+ssum)

            end if

            ssum = stmp
        
        end do

        ssum = ssum + serr

        return

    end function

    subroutine linspace(xxll,xxuu,npos,xpos)

    ! 
    ! XXLL  lower-bound grid position.
    ! NNEW  upper-bound grid position.
    ! NPOS  no. edges in the grid.
    ! XPOS  array of grid edges. XPOS has length NPOS .
    !

        implicit none

        real(kind=dp), intent(in)   :: xxll,xxuu
        integer      , intent(in)   :: npos
        real(kind=dp), intent(out)  :: xpos(:)

        integer                     :: ipos
        real(kind=dp)               :: xdel

        xpos(   1) = xxll
        xpos(npos) = xxuu

        xdel = (xxuu-xxll) / (npos - 1)

        do  ipos = +2, npos-1

            xpos(ipos) = (ipos-1) * xdel         

        end do

        return
    
    end  subroutine

    subroutine rndspace(xxll,xxuu,npos,xpos, &
        &               frac)

    ! 
    ! XXLL  lower-bound grid position.
    ! NNEW  upper-bound grid position.
    ! NPOS  no. edges in the grid.
    ! XPOS  array of grid edges. XPOS has length NPOS .
    ! FRAC  fractional perturbation of cell, OPTIONAL .
    !

        implicit none

        real(kind=dp), intent(in)   :: xxll,xxuu
        integer      , intent(in)   :: npos
        real(kind=dp), intent(out)  :: xpos(:)
        real(kind=dp), intent(in), optional :: frac

        integer                     :: ipos
        real(kind=dp)               :: xdel
        real(kind=dp)               :: rand,move

        if (present(frac)) then
            move = +frac
        else
            move = 3.0d0 / 8.0d0
        end if

        xpos(   1) = xxll
        xpos(npos) = xxuu

        xdel = (xxuu-xxll) / (npos - 1)

        do  ipos = +2, npos-1

            xpos(ipos) = (ipos-1) * xdel         

        end do

        do ipos = +2, npos-1

            call random_number (rand)

            rand = 2.d0 * (rand-.5d0)
            
            xpos(ipos) = &
        &       xpos(ipos) + (move * rand * xdel)         

        end do

        return

    end  subroutine



