!!!!!!!!!!!!
! revised by Fang Chang and Weiliang Qiu Oct. 16, 2009 
!    (1) changed from fotran77 format to fortan95 format 

! created by Weiliang Qiu on Jan. 17, 2009      
!  stwxq@channing.harvard.edu      
!
!!!!!!!!!!!! 
!    
!      
!    input:
!      mem: cluster each data point belonging to
!    output:
!      CH: CH index value

subroutine CHindex(dat, nObs, nVars, mem, nClusters, &
    & clustSize, disMethod2, CH)
 
    implicit none
    integer, intent(in) :: nObs, nVars, nClusters, disMethod2
    integer, intent(in) :: mem(nObs), clustSize(nClusters)
    real(8), intent(in) :: dat(nObs, nVars)
    real(8), intent(out) :: CH
    
    integer :: i, j, k, nk
    real(8) :: mu(nClusters, nVars), bss, wss, muAll(nVars)
    real(8) :: tmp(nClusters), denom
    real(8) :: sumx, sumx2, sumy, sumy2, sumxy, rho
    real(8) :: numer, denom1, denom2, D, colj(nObs), tmpx(nVars)
    real(8) :: tmpy(nVars)
    integer, dimension(:), save, allocatable :: posk
    real(8), dimension(:,:), save, allocatable :: datk
    
    ! initialize mu = 0
    mu = 0
 
    ! for each colum, calculate average
    muAll = sum(dat, dim = 1) / nObs
 
    ! calculate cluster mean 'mu'
    do k = 1, nClusters
        do j = 1, nVars
            ! obtain the j-th colum of 'dat'
            colj = dat(1:nObs, j)
            mu(k, j) = sum(colj, mask = (mem == k)) / clustSize(k)
        end do
    end do
 
    ! calculate between-cluster distance
    bss = 0.0
    sumy = sum(muAll)
    sumy2 = sum(muAll ** 2)
    do k = 1, nClusters
        ! obtain the k-th row of 'mu'
        tmpx = mu(k, 1:nVars)
        sumx = sum(tmpx)
        sumx2 = sum(tmpx ** 2)
        sumxy = sum(tmpx * muAll)
       
!!      Euclidean distance
        if(disMethod2 == 1) then
            D = sumx2 + sumy2 - 2 * sumxy
            D = dsqrt(D)
        else ! 1 - correlation
            numer = nVars * sumxy - sumx * sumy
            denom1 = nVars * sumx2 - sumx ** 2
            denom2 = nVars * sumy2 - sumy ** 2
            rho = numer / dsqrt(denom1 * denom2)
            ! some times, 1-rho is very close to zero, but be negative
            ! due to computer representation
            ! Hence use abs        
            D = abs(1 - rho)
        endif
        bss = bss + clustSize(k) * (D ** 2)
    end do
 
    ! calculate within-cluster distance
    wss = 0.0
    do k = 1, nClusters
        nk = clustSize(k)
        ! get data for the k-th cluster
        j = 0
        allocate(posk(nk))
        do i = 1, nObs
            if(mem(i) == k) then
                j = j + 1
                posk(j) = i
            endif
        end do
        allocate(datk(nk, nVars))
        datk = dat(posk, 1:nVars)
        deallocate(posk)
       
        tmpy = mu(k, 1:nVars)
        sumy = sum(tmpy)
        sumy2 = sum(tmpy ** 2)
       
        tmp(k) = 0.0
       
        do i = 1, nk 
            tmpx = datk(i, 1:nVars)
            sumx = sum(tmpx)
            sumx2 = sum(tmpx ** 2)
            sumxy = sum(tmpx * tmpy)
           
            ! Euclidean distance
            if(disMethod2 == 1) then
                D = sumx2 + sumy2 - 2 * sumxy
                D = dsqrt(D)
            else 
                numer = nVars * sumxy - sumx * sumy
                denom1 = nVars * sumx2 - sumx ** 2
                denom2 = nVars * sumy2 - sumy ** 2
                rho = numer / dsqrt(denom1 * denom2)
                D = 1 - rho 
            endif
            tmp(k) = tmp(k) + D ** 2
        end do
        wss = wss + tmp(k)
        deallocate(datk)
    end do
 
    numer = bss / (nClusters - 1.0)
    denom = wss / (nObs - nClusters)
    CH = numer / denom

end subroutine CHindex

