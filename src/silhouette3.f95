!!!!!!!!!!!!
! revised by Fang Chang and Weiliang Qiu Oct. 16, 2009 
!     (1) changed from fotran77 format to fortan95 format 
!      
! created by Weiliang Qiu on Jan. 15, 2009      
!     stwxq@channing.harvard.edu      
!
!@book{Kaufman:1990,
!    author = "{Kaufman, L.} and {Rousseeuw, P. J.}",
!    title = {Finding Groups in Data: An Introduction to Cluster Analysis},
!    publisher = {Wiley, New York.},
!    year = {1990},
!}
!!!!!!!!!!!! 
! 
!!!!!!!!
! For each observation i, the silhouette width s(i) is defined as
!     follows:       
!     Put a(i) = average dissimilarity between i and all other points of
!     the cluster to which i belongs.  For all other clusters C, put
!     d(i,C) = average dissimilarity of i to all observations of C.  The
!     smallest of these d(i,C) is b(i) := min_C d(i,C), and can be seen
!     as the dissimilarity between i and its ``neighbor'' cluster, i.e.,
!     the nearest one to which it does not belong. Finally,
!    
!             s(i) := ( b(i) - a(i) ) / max( a(i), b(i) ).
!     
!     Observations with a large s(i) (almost 1) are very well clustered,
!     a small s(i) (around 0) means that the observation lies between
!     two clusters, and observations with a negative s(i) are probably
!     placed in the wrong cluster.
! --- the above is copied from R document on the function silhoutte() in
! library(cluster)
!!!!!!!!!
!    
!    input:
!        mem: cluster each data point belonging to
!    output:
!        memNei: nearest neighbor cluster of each data point      
!        sIndex: silhoutte index for each data point
!        avgs: average silhoutte index

subroutine silhouette(dat, nObs, nVars, mem, nClusters, clustSize, &
                    & disMethod2, memNei, sIndex, avgs)
    implicit none
    
    integer, intent(in) :: nObs, nVars, nClusters, disMethod2
    integer, intent(in) :: mem(nObs), clustSize(nClusters)
    real(8), intent(in) :: dat(nObs, nVars)
    integer, intent(out) :: memNei(nObs)
    real(8), intent(out) :: sIndex(nObs), avgs
    
    integer :: i, j, memi, sizemi, memj, j2, sizemj, m
    integer :: i2, ni1, nj
    real(8) :: dik, obsi(nVars), dij, a(nObs), b(nObs)
    real(8) :: sumx, sumx2, sumy, sumy2, sumxy, rho
    real(8) :: numer, denom1, denom2, D
    real(8) :: tmpx(nVars)
    integer, dimension(:), save, allocatable :: posi
    real(8), dimension(:,:), save, allocatable :: dati
    
    m = nVars
    ! initializing memNei
    memNei = 0
    
    avgs = 0.0
    do i = 1, nObs
        ! obtain the i-th observation (row) of dat      
        obsi = dat(i, 1:nVars)
        sumy = sum(obsi)
        sumy2 = sum(obsi ** 2)
     
        memi = mem(i)
        sizemi = clustSize(memi)
        if(sizemi == 1) then
            ! only one data point in cluster 'memi', then a(i) = 0        
            a(i) = 0.0
        else
            ! get the average distance from point i to other data points in the
            ! same cluster as point i
           
            ! get data for the memi-th cluster, except the data point i itself
            ni1 = clustSize(memi) - 1
            j = 0
            allocate(posi(ni1))
            do i2 = 1, nObs
                if(mem(i2) == memi .and. i2 /= i) then
                    j = j + 1
                    posi(j) = i2
                endif
            end do
            allocate(dati(ni1, nVars))
            dati = dat(posi, 1:nVars)
            deallocate(posi)
           
            dik = 0.0
            do j = 1, ni1
                ! get distance between point i and point j                    
                tmpx = dati(j, 1:nVars)
                sumx = sum(tmpx)
                sumx2 = sum(tmpx ** 2)
                sumxy = sum(tmpx * obsi)        
               
                ! Euclidean distance
                if(disMethod2 == 1) then
                    D = sumx2 + sumy2 - 2 * sumxy
                    D = dsqrt(D)
                else 
                    numer = m * sumxy - sumx * sumy
                    denom1 = m * sumx2 - sumx ** 2
                    denom2 = m * sumy2 - sumy ** 2
                    rho = numer / dsqrt(denom1 * denom2)
                    ! some times, 1-rho is very close to zero, but be negative
                    ! due to computer representation
                    ! Hence use abs            
                    D = abs(1 - rho)
                endif
                dij = D
                dik = dik + dij
            end do         
            a(i) = dik / (sizemi - 1.0)
            deallocate(dati)
        endif          
     
        ! get the average distance from point i to data points in the
        ! nearest neighbor cluster besides its own
        b(i) = 1.0E+30
        do j = 1, nClusters
            memj = j
            if (memi /= memj) then
                nj = clustSize(memj)
                ! get data for the memj-th cluster
                j2 = 0
                allocate(posi(nj))
                do i2 = 1, nObs
                    if(mem(i2) == memj) then
                        j2 = j2 + 1
                        posi(j2) = i2
                    endif
                end do
                allocate(dati(nj, nVars))
                dati = dat(posi, 1:nVars)
                deallocate(posi)
               
                dik = 0.0
                do j2 = 1, nj
                    ! distance between point i and point j2 in cluster memj
                    tmpx = dati(j2, 1:nVars)
                    sumx = sum(tmpx)
                    sumx2 = sum(tmpx ** 2)
                    sumxy = sum(tmpx * obsi)        
                    ! Euclidean distance
                    if(disMethod2 == 1) then
                        D = sumx2 + sumy2 - 2 * sumxy
                        D = dsqrt(D)
                    else 
                        numer = m * sumxy - sumx * sumy
                        denom1 = m * sumx2 - sumx ** 2
                        denom2 = m * sumy2 - sumy ** 2
                        rho = numer / dsqrt(denom1 * denom2)
                        D = 1 - rho 
                    endif
                    dij = D
                    dik = dik + dij
                end do            
                deallocate(dati)
                sizemj = clustSize(j)
                dik = dik / (sizemj * 1.0)
                if( b(i) > dik ) then
                    b(i) = dik
                    ! data point i's nearest neighor cluster is cluster 'memj'
                    memNei(i) = memj
                endif
            endif
        end do        
        sIndex(i) = (b(i) - a(i)) / max(b(i), a(i))
        avgs = avgs + sIndex(i)
    end do      
    avgs = avgs / (nObs * 1.0)
end subroutine silhouette

