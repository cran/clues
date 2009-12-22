!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! revised by Fang Chang and Weiliang Qiu Oct. 22, 2009 
!    (1) changed from fotran77 format to fortan95 format 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  nObs1=nObs-1
!
! we can plot db against {1, 2, ..., nObs1}
!   to check the distance between consecutive data points      
!      
!
! outputs:
!    points -- INTEGER. a sequence of labels of data points
!    db -- real*8. distance between consecutive points indicated by 'points'
!    omin -- real*8 if db(k) > omin, then the data point points(k) and
!      points(k+1) will be in two different clusters. That is,
!      a new cluster will start with points(k+1)      
!    nClusters -- INTEGER. number of clusters      
!    mem -- INTEGER. cluster membership      
!    clustSize -- INTEGER. size of clusters (i.e., number of data points
!      in each cluster      
!     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      

subroutine clustering(dat, nObs, nObs1, nVars, disMethod2, points, &
    &  db, omin, nClusters, mem, clustSize)
    implicit none
    
    integer, intent(in) :: nObs, nObs1, nVars, disMethod2
    real(8), intent(in) :: dat(nObs, nVars)
    integer, intent(out) :: points(nObs), nClusters, mem(nObs), clustSize(nObs)
    real(8), intent(out) :: db(nObs1), omin
        
    integer :: nSstar, nNei, nNei2, Sstar(nObs)
    integer :: i, myt, sumS, k
    integer :: m1, m2, nOutliers
    integer :: kList(2), ell, S(nObs)
    integer :: pos1, pos2, L1, L2
    real(8) :: obsi(nVars), upp, low
    real(8) :: dk(2), q1, q3, IQR
    real(8) :: dat2(nObs, nVars), db2(nObs1)
    
    nNei = 1
    nNei2 = 2
    
    ! Step 1. get pair-wise distances dij among shrinked data points
    points = 0
    S = (/ (i, i = 1, nObs, 1) /)
    sumS = sum(S)
    db(1:nObs1) = S(1:nObs1)
    
    ! Step 2. initialization
    myt = 1
    points(myt) = 1
    S(myt) = 0
    
    ! Step 3. get db ( (n-1)x1) vector) recording the distance
    !     between a point and its nearest neighbor point.
    !     The process will start from the first data point (denoted as p0).
    !     Then find its nearest neighbor point (denoted as p1) and
    !     record the distance to db[1]. Then find the nearest neighobr
    !     point (denoted as p2) of p1. record the distance to db[2].
    !     and so on.
    !
    do while(sumS > 0)
        ! find the nearest neighbor x_{\ell} \in S* of 
        ! the data point x_{point_t}  
        nSstar = 0
        do i = 1, nObs
            if (S(i) > 0) then
                nSstar = nSstar + 1
                ! set of remaining points            
                Sstar(nSstar) = S(i)
                dat2(nSstar, 1:nVars) = dat(S(i), 1:nVars)
            endif
        end do
 
        ! find 'nNei2' nearest neighbor of obsi        
        obsi = dat(points(myt), 1:nVars)
        call dist(dat2, obsi, nObs, nVars, nSstar, nNei2, disMethod2, kList, dk)
        ell = Sstar(kList(1))
        db(myt) = dk(1)
        myt = myt + 1
        points(myt) = ell
        S(ell) = 0
     
        sumS = sum(S)
    enddo 
    
    
    ! Step 4
    ! find the outliers of db.  
    ! calculate sample quantiles of db
    ! first sort 'db'
    call dsort(db, nObs1, db2)
    ! 25% percentile
    m1 = mod(nObs1, 4)
    if(m1 > 0) then
        ! round up
        pos1 = nObs1 / 4 + 1
    else
        pos1 = nObs1 / 4
    endif
    q1 = db2(pos1)
    
    ! 75% percentile
    m1 = 3 * nObs1
    m2 = mod(m1, 4) 
    if(m2 > 0) then
        ! round up
        pos2 = m1 / 4 + 1
    else
        pos2 = m1 / 4
    endif
    q3 = db2(pos2)
    
    ! interquantile range
    IQR = q3 - q1
    
    upp = q3 + 1.5 * IQR
    low = max(0.0, q1 - 1.5 * IQR)
    
    ! minimum of outlier distances 
    omin = minval(db, mask = (db > upp) .or. (db < low))
 
    ! outliers are defined as db(k)>q3+1.5*IQR or db(k)<q1-1.5*IQR      
    nOutliers = count((db > upp) .or. (db < low))
    
    ! Step 5. check if there is any cluster structure
    if (nOutliers == 0) then
        mem = 1
        nClusters = 1
    else
        ! Step 6
        ! initialize mem and nClusters
        mem = 0
        nClusters = 1
        myt = 1
        mem(myt) = 1
     
        L1 = 1
        L2 = 2
        do while (L2 > L1) 
            ! Step 7. get number of clusters
            if (db(myt) >= omin) then
                nClusters = nClusters + 1
            endif
            ! Step 8. get cluster membership          
            m1 = myt + 1
            k = points(m1)
            mem(k) = nClusters
            myt = myt + 1
            ! Step 9. all data points have been assigned membership. stop 
            if (myt >= nObs) then
                exit
            endif
        end do
     
        ! Step 10
        ! calculate the cluster size                  
        clustSize = 0
        do k = 1, nClusters
            clustSize(k) = count(mem == k)
        end do   
    endif
end subroutine clustering
   
  
