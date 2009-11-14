! This fortran code is based on
! http: /  / www.itl.nist.gov / div898 / software / datapac / sort.f      
! revised by Fang Chang and Weiliang Qiu Oct. 16, 2009 
!    (1) changed from fotran77 format to fortan95 format 
!    (2) goto statement has been avoided.
!
! revise 'sort' to 'dsort' for double data by weiliang qiu on jan 14,
! 2009


subroutine dsort(x, n, y)
!
!--------------------------------------------------------
!
!     purpose--this subroutine sorts (in ascending order)
!              the n elements of the double precision vector x
!              and puts the resulting n sorted values into the
!              double precision vector y.
!     input  arguments--x      = the double precision vector of
!                                observations to be sorted. 
!                     --n      = the integer number of observations
!                                in the vector x. 
!     output arguments--y      = the double precision vector
!                                into which the sorted data values
!                                from x will be placed.
!     output--the double precision vector y
!             containing the sorted
!             (in ascending order) values
!             of the double precision vector x.
!     printing--none unless an input argument error condition exists. 
!     restrictions--the dimensions of the vectors il and iu 
!                   (defined and used internally within
!                   this subroutine) dictate the maximum
!                   allowable value of n for this subroutine.
!                   if il and iu each have dimension k,
!                   then n may not exceed 2**(k+1) - 1.
!                   for this subroutine as written, the dimensions
!                   of il and iu have been set to 36,
!                   thus the maximum allowable value of n is
!                   approximately 137 billion.
!                   since this exceeds the maximum allowable
!                   value for an integer variable in many computers,
!                   and since a sort of 137 billion elements
!                   is presently impractical and unlikely,
!                   then there is no practical restriction
!                   on the maximum value of n for this subroutine.
!                   (in light of the above, no check of the 
!                   upper limit of n has been incorporated
!                   into this subroutine.)
!     other datapac   subroutines needed--none.
!     fortran library subroutines needed--none.
!     mode of internal operations--double precision.
!     language--ansi fortran. 
!     comment--the smallest element of the vector x
!              will be placed in the first position
!              of the vector y,
!              the second smallest element in the vector x
!              will be placed in the second position
!              of the vector y, etc.
!     comment--the input vector x remains unaltered.
!     comment--if the analyst desires a sort 'in place',
!              this is done by having the same
!              output vector as input vector in the calling sequence. 
!              thus, for example, the calling sequence
!              call sort(x,n,x)
!              is allowable and will result in
!              the desired 'in-place' sort.
!     comment--the sorting algorthm used herein
!              is the binary sort.
!              this algorthim is extremely fast as the
!              following time trials indicate.
!              these time trials were carried out on the
!              univac 1108 exec 8 system at nbs
!              in august of 1974.
!              by way of comparison, the time trial values
!              for the easy-to-program but extremely
!              inefficient bubble sort algorithm have
!              also been included--
!              number of random        binary sort       bubble sort
!               numbers sorted
!                n = 10                 .002 sec          .002 sec
!                n = 100                .011 sec          .045 sec
!                n = 1000               .141 sec         4.332 sec
!                n = 3000               .476 sec        37.683 sec
!                n = 10000             1.887 sec      not computed
!     references--cacm march 1969, page 186 (binary sort algorithm
!                 by richard c. singleton).
!               --cacm january 1970, page 54.
!               --cacm october 1970, page 624.
!               --jacm january 1961, page 41.
!     written by--james j. filliben
!                 statistical engineering laboratory (205.03)
!                 national bureau of standards
!                 washington, D. C. 20234
!                 phone--301-921-2315
!     original version--june      1972. 
!     updated         --november  1975. 
!
!---------------------------------------------------------------------
!
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: x(n)
    real(8), intent(out) :: y(n)
    integer :: i, j, l, k, lmi, jmk, m, jmi, mid, ipr, ip1, nm1
    logical :: iflag, iflag310
    real(8) amed, tt, hold
    integer :: iu(36), il(36) 
    character(len = *), parameter :: fm15 = "(a109)"
    character(len = *), parameter :: fm47 = "(a109, i8, a6)"
    
    ipr = 6
    
    ! check the input arguments for errors
          
    if(n < 1) then
        write(ipr, fm15) "***** fatal error--the second input argument to & 
            & the sort   subroutine is non-positive *****"  
        write(ipr, fm47) "***** the value of the argument is", n, "****" 
        return          
    elseif(n == 1) then
        y = x
        return
    else
        hold = x(1)
        iflag = .FALSE.
        do i = 2, n
            if(x(i) /= hold) then
                iflag = .TRUE.
                exit
            endif
        enddo
        if(.not. iflag) then
            y = x
            return
        else
            ! copy the vector x into the vector y
            y = x
           
            ! check to see if the input vector is already sorted
            nm1 = n - 1
            iflag = .FALSE.
            do i = 1, nm1
                ip1 = i + 1
                if(y(i) > y(ip1)) then
                    iflag = .TRUE.
                    exit
                endif
            enddo
            if(.not. iflag) return
            m = 1 
            i = 1 
            j = n 
           
            iflag310 = .FALSE. 
            do 
                if((.not. iflag310) .and. (i >= j)) then
                    m = m - 1
                    if(m == 0) return
                    i = il(m)
                    j = iu(m)
                else
                    iflag310 = .FALSE.
                    k = i
                    mid = (i + j) / 2
                    amed = y(mid)
                    if(y(i) > amed) then
                        y(mid) = y(i)
                        y(i) = amed
                        amed = y(mid)
                    endif
                    l = j 
                    iflag = .FALSE.
                    if(y(j) < amed) then
                        y(mid) = y(j)
                        y(j) = amed
                        amed = y(mid)
                        if(y(i) > amed) then
                            y(mid) = y(i)
                            y(i) = amed
                            amed = y(mid)
                        endif
                    endif
                    do
                        if(iflag) then
                            y(l) = y(k)
                            y(k) = tt
                        endif
                        l = l - 1        
                        if(y(l) > amed) then
                            iflag = .FALSE.
                            cycle
                        endif
                        tt = y(l)
                        do
                            k = k + 1 
                            if(y(k) >= amed) exit
                        enddo 
                        if(k > l) exit
                        iflag = .TRUE.
                    enddo
                    lmi = l - i
                    jmk = j - k
                    if(lmi <= jmk) then
                        il(m) = k
                        iu(m) = j
                        j = l 
                        m = m + 1
                    else
                        il(m) = i
                        iu(m) = l
                        i = k 
                        m = m + 1
                    endif
                endif
                jmi = j - i
                if(jmi >= 11) then
                    iflag310 = .TRUE.
                    cycle
                endif
                if(i == 1) then
                    iflag310 = .FALSE.
                else
                    i = i -  1
                    do
                        i = i + 1
                        if(i == j) then
                            iflag310 = .FALSE.
                            exit
                        endif
                        amed = y(i + 1)
                        if(y(i) <= amed) cycle
                        k = i 
                        do
                            y(k + 1) = y(k)
                            k = k - 1
                            if(amed >= y(k)) exit
                        enddo
                        y(k + 1) = amed
                    enddo
                endif 
            enddo
        endif
    endif
end subroutine dsort

