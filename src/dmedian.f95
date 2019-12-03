! This fortran code is based on 
! http: /  / www.itl.nist.gov / div898 / software / datapac / median.f      
! revised by Fang Chang and Weiliang Qiu Oct. 16, 2009 
!    (1) changed from fotran77 format to fortan95 format 
!    (2) goto statement has been avoided.
!    (3) replace the constraint:
!     "restrictions--the maximum allowable value of n
!                   for this subroutine is 15000."
!     by 500,000 
!
!  revise 'median' to 'dmedian' by weiliang qiu on jan. 14, 2009


subroutine dmedian(x, n, iwrite, xmed)
!
!---------------------------------------------------------------------
!
!     purpose--this subroutine computes the
!              sample median
!              of the data in the input vector x. 
!              the sample median = that value such that half the
!              data set is below it and half above it.
!     input  arguments--x      = the DOUBLE precision vector of
!                                (unsorted or sorted) observations.
!                     --n      = the integer number of observations
!                                in the vector x. 
!                     --iwrite = an integer flag code which 
!                                (if set to 0) will suppress
!                                the printing of the
!                                sample median
!                                as it is computed;
!                                or (if set to some integer 
!                                value not equal to 0),
!                                like, say, 1) will cause
!                                the printing of the
!                                sample median
!                                at the time it is computed.
!     output arguments--xmed   = the DOUBLE precision value of the
!                                computed sample median.
!     output--the computed double precision value of the
!             sample median.
!     printing--none, unless iwrite has been set to a non-zero
!               integer, or unless an input argument error
!               condition exists.
!     restrictions--the maximum allowable value of n
!                   for this subroutine is 15000. 
!     other datapac   subroutines needed--sort.
!     fortran library subroutines needed--none.
!     mode of internal operations--double precision.
!     language--ansi fortran. 
!     references--kendall and stuart, the advanced theory of
!                 statistics, volume 1, edition 2, 1963, page 326.
!               --kendall and stuart, the advanced theory of
!                 statistics, volume 2, edition 1, 1961, page 49.
!               --david, order statistics, 1970, page 139.
!               --snedecor and cochran, statistical methods,
!                 edition 6, 1967, page 123.
!               --dixon and massey, introduction to statistical
!                 analysis, edition 2, 1957, page 70.
!     written by--james j. filliben
!                 statistical engineering laboratory (205.03)
!                 national bureau of standards
!                 washington, D. C. 20234
!                 phone:  301-921-2315
!     original version--june      1972. 
!     updated         --september 1975. 
!     updated         --november  1975. 
!     updated         --february  1976. 
!---------------------------------------------------------------------
!
    implicit none
    integer, intent(in) :: n, iwrite
    real(8), intent(in) :: x(n)
    real(8), intent(out) :: xmed
    integer :: i, ipr, iupper, iflag, nmid, nmidp1
    logical :: icheck
    real(8) :: y(n), hold 
    character(len = *), parameter :: fm17 = "(a109, i6)"
    character(len = *), parameter :: fm47 = "(a35, i8, a6)"
    character(len = *), parameter :: fm999 = "(1h )"
    character(len = *), parameter :: fm105 = "(1h ,25hthe sample median &
        & of the ,i6,17h observations is e15.8)"
    
    ipr = 6
    iupper = 500000
          
    ! check the input arguments for errors
    if((n < 1) .or. (n > iupper)) then
        !write(ipr, fm17) "***** fatal error--the second input argument to & 
        !    & the median subroutine is outside the allowable ", iupper,  &
        !    & " interval *****"
        !write(ipr, fm47) "***** the value of the argument is ", n, " *****"
        call intpr("**** fatal error--the second input argument to & 
            & the median subroutine is outside the allowable &
            &  interval *****", -1, 0, 0)
        !write(ipr, fm47) "***** the value of the argument is ", n, " *****"
        call intpr("***** the value of the argument is ", -1, n, 6)
        return
    elseif(n > 1) then
        hold = x(1)
        icheck = .FALSE. 
        do i = 2, n
            if(x(i) /= hold) then
                icheck = .TRUE.
                exit
            endif
        enddo
        if (icheck) then
            call dsort(x, n, y)
            iflag = mod(n, 2)
            nmid = n / 2
            nmidp1 = nmid + 1 
            if(iflag == 0) xmed = (y(nmid) + y(nmidp1)) / 2.0
            if(iflag == 1) xmed = y(nmidp1)
        else
            xmed = x(1)
        endif
    else
        xmed = x(1)
    endif
    
    if(iwrite == 0) return

end subroutine dmedian 


