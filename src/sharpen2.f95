!!!!!!!!!!!!
! revised by Fang Chang and Weiliang Qiu Oct. 16, 2009 
!    (1) changed from fotran77 format to fortan95 format 
!!!!!!!!!!!!


!-------------------------------------------------------------------------
!!!!!!  The subroutine 'DIST' is a revised version of 'DIST'
!!!!!!   at http://astro.u-strasbg.fr/~fmurtagh/mda-sw/knn.f      
!      
!!!!!!  This a brutal-force method. Hence it is not efficient.
!      
!!!!!!! Find the K nearest neighbors of 'obsi'.
!!!!!!! The positions of these neighbors are recorded in 'kList' 
!!!!!!! The corresponding distance between 'obsi' and these neighbors
!!!!!!! are recorded in 'dk'      

!  dat(N,M)    data set, where N is the number of observations (rows),
!                 M is the number of variables (columns)
!  obsi(M)     the i-th observation;    
!  K              number of nearest neighbours to consider;
!  kList(K), dk(K)   are used for storing the K NNs
!                 and their distances to the object under 
!                 consideration.
!
! sometimes, we only interested in the first N2 rows of dat    

subroutine dist(dat, obsi, N, M, N2, K, disMethod2, kList, dk)
    implicit none
    
    integer, intent(in) :: N, M, N2, K, disMethod2
    real(8), intent(in) :: dat(N, M), obsi(M)
    integer, intent(out) :: kList(K)
    real(8), intent(out) :: dk(K)
    
    integer iTrain, iLoc, IiLoc
    real(8) D, sumx, sumx2, sumy, sumy2, sumxy, rho
    real(8) numer, denom1, denom2, tmpx(M)
    
    ! initialize 
    kList = 0
    dk = 1.E+15
    
    sumy = sum(obsi)
    sumy2 = sum(obsi ** 2)
    do iTrain = 1, N2
        tmpx = dat(iTrain, 1:M)
        sumx = sum(tmpx)
        sumx2 = sum(tmpx ** 2)
        sumxy = sum(tmpx * obsi)
     
        ! Euclidean distance
        if(disMethod2 == 1) then
            D = sumx2 + sumy2 - 2 * sumxy
            D = dsqrt(D)
        else 
            !1-correlation or 1-rankcorrelation
            !rho=numer/sqrt(denom1*denom2), where
            !numer=n*sum(x*y) - sum(x)*sum(y) 
            !denom1=n*sum(x^2)-(sum(x))^2
            !denom2=n*sum(y^2)-(sum(y))^2
            
            numer = M * sumxy - sumx * sumy
            denom1 = M * sumx2 - sumx **  2
            denom2 = M * sumy2 - sumy **  2
            rho = numer /  dsqrt(denom1 * denom2)
            ! some times, 1-rho is very close to zero, but be negative
            ! due to computer representation
            ! Hence use abs
            D = abs(1 - rho) 
        endif
     
        do iLoc = 1, K
            if (D < dk(iLoc)) then
                ! Insert at locn. 
                ! ILOC and shift right in the 3 length-k lists we're 
                ! maintaining            
                do IiLoc = K, iLoc + 1,  -1
                    ! Protective measure:      
                    if (IiLoc <= K) then
                        dk(IiLoc) = dk(IiLoc - 1)
                        kList(IiLoc) = kList(IiLoc - 1)
                    endif
                    !Have now freed up space at locn. ILOC        
                end do
                dk(iLoc) = D
                kList(iLoc) = iTrain
                exit
            endif
        end do
    end do
end subroutine dist
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!       nNei2=nNei+1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
 
subroutine sharpen(dat, nObs, nVars, nNei2, nNei, ITMAX, &
    & eps, disMethod2, datnew)
 
    implicit none
    integer, intent(in) :: nObs, nVars, disMethod2, nNei2, nNei, ITMAX
    real(8), intent(in) :: dat(nObs, nVars), eps
    real(8), intent(out) :: datnew(nObs, nVars)
    
    integer :: i, ix, iter
    integer :: kList(nNei2)
    real(8) :: dk(nNei2), varix(nNei), tmp, maxdiff
    real(8) :: datold(nObs, nVars), obsi(nVars)
    
    ! step 1. backup data matrix
    datold = dat
    datnew = dat
    
    ! initialize the iteration number
    iter = 1
    
    ! step 2      
    do while (iter <= ITMAX)
        maxdiff = -1.0
        do i = 1, nObs
            ! find K = nNei nearest neighbors of the data point i         
            ! get i-th observation
            obsi = datold(i, 1:nVars)
         
            ! the K = nNei nearest neighbors of 'obsi' are stored in 'kList'
            ! (not include itsself)
            call dist(datold, obsi, nObs, nVars, nObs, nNei2, &
                & disMethod2, kList, dk)
           
            ! update obsi by using the coordinate-wise median of
            ! the 'nNei' nearest neighbors                
            do ix = 1, nVars
                varix = datold(kList(2:nNei2), ix)
                call dmedian(varix, nNei, 0, tmp)
                datnew(i, ix) = tmp
                tmp = dabs(datnew(i, ix) - datold(i, ix))
                if (tmp > maxdiff ) then
                    maxdiff = tmp
                endif
            end do
        end do
       
        ! Step 3
        ! update datold
        datold = datnew
     
        if (maxdiff < eps) exit
        iter = iter + 1
    enddo
end subroutine sharpen


