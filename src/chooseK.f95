!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! revised by Fang Chang and Weiliang Qiu Oct. 22, 2009 
!    (1) changed from fotran77 format to fortan95 format 
!    (2) combining chooseK_CH3.f and chooseK_sil3.f
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Choose an appropriate K of K nearest neighbors based on
! CH index (indexFlag = 0) or silhouette index (indexFlag = 1)
!    input:
!      nClusters0: initial guess of the number of clusters
!      alpha: proportion of data as the K of K nearest neighbors
!      eps: used for the subroutine 'sharpen' 
!      ITMAX: maximum iteration number used for 'sharpen'
!      second: indicating if the second pass of the iteration is required  
!      nNeiVec2: nNeiVec2(1) and nNeiVec2(2) are the lower and upper bounds for 
!        the number of nearest neighbors for the second pass of the iteration.
!      s2: initial value of the lower bound of the silhouette index

!    output:
!      memNei: nearest neighbor cluster of each data point      
!      sIndex: silhoutte index for each data point
!      avgs: average silhoutte index

subroutine chooseK(dat, dat2, nObs, nObs1, nVars, nClusters0, &
    &  alpha, eps, ITMAX, second, nNeiVec2, indInitial, &
    &  disMethod2, indexFlag, quiet, indFinal, sFinal, &
    &  memFinal, nClustersFinal, clustSizeFinal, nNeiFinal, & 
    &  nClustVec, nNeiVec, myt, datold1, datold2, myupdate)
    implicit none
    integer, intent(in) :: nObs, nObs1, nVars, nClusters0, ITMAX 
    integer, intent(in) :: nNeiVec2(2), disMethod2
    real(8), intent(in) :: dat(nObs, nVars), dat2(nObs, nVars), alpha, eps
    real(8), intent(in) :: indInitial
    logical, intent(in) :: second, quiet, indexFlag
 
    integer, intent(out) :: memFinal(nObs), nClustersFinal, clustSizeFinal(nObs)
    integer, intent(out) :: nNeiFinal, nClustVec(nObs), nNeiVec(nObs), myt
    real(8), intent(out) :: indFinal, sFinal(nObs), datold1(nObs, nVars) 
    real(8), intent(out) :: datold2(nObs, nVars)
    logical, intent(inout) :: myupdate
 
    integer :: i, delta, myt2, L1, L2, nNei, nNei2, minsize
    integer :: mem(nObs), clustSize(nObs), memNei(nObs)
    integer :: points(nObs), nClusters
    real(8) :: avgIndex, sIndex(nObs), datnew(nObs, nVars)
    real(8) :: datnew2(nObs, nVars), db(nObs1), omin
        
    ! step 1
    datnew = dat2
    datold1 = dat2
    datold2 = dat2
 
    ! initialize the number K of the nearest neighbors.
    ! and initialize the increment delta of the number of 
    ! the nearest neighbors
    if (.not. second) then
        nNei = dnint(alpha * nObs)
        delta = max(nNei, nClusters0)
    else
        nNei = nNeiVec2(1)
        i = dnint(0.1 * alpha * nObs)
        delta = max(i, 1)
        if(.not. quiet) then
            write(*, *) 'second = TRUE, nNei = ', nNei, ' i = ', i, &
                & ' delta = ', delta
        endif
    endif 
 
    ! initialize the iteration number
    myt = 0
    indFinal = -2.0
    myt2 = 0
 
    L1 = 1
    L2 = 2
    do while (L2 > L1)
        myt = myt + 1
        nNeiVec(myt) = nNei
       
        ! Step 2 shrinking
        nNei2 = nNei + 1
        call sharpen(datnew, nObs, nVars, nNei2, nNei, ITMAX, eps, &
            & disMethod2, datnew2)
        datnew = datnew2
       
        ! step 3 clustering
        call clustering(datnew, nObs, nObs1, nVars, disMethod2, &
           & points, db, omin, nClusters, mem, clustSize)
        nClustVec(myt) = nClusters
       
        if(nClusters > 1) then
            !step 4.1
            !use original data to calculate index
            if(indexFlag) then
                call silhouette(dat, nObs, nVars, mem, nClusters, &
                   & clustSize, disMethod2, memNei, sIndex, avgIndex)
            else
                call CHindex(dat, nObs, nVars, mem, nClusters, &
                   & clustSize, disMethod2, avgIndex)
            endif
           
            minsize = minval(clustSize(1:nClusters))
           
            if(.not. quiet) then
                write(*, *) 'myt = ', myt, ' nNei = ', nNei, & 
                    & ' delta = ', delta, ' minsize = ', minsize, &
                    & ' avgIndex = ', avgIndex, ' nClusters = ', nClusters
            endif
           
            ! step 4.2
            if (avgIndex > indFinal) then 
                indFinal = avgIndex
                memFinal = mem
               
                if(indexFlag) then
                    sFinal = sIndex
                endif
               
                clustSizeFinal = 0
                clustSizeFinal(1:nClusters) = clustSize(1:nClusters)
                nClustersFinal = nClusters
                nNeiFinal = nNei
                myt2 = myt2 + 1
            endif
           
            if (myt2 == 2) then
                myt2 = 1
                datold1 = datold2
            endif
           
            if (avgIndex >= indFinal) then
                datold2 = datnew
            endif
           
            ! step 4.3
            if (nClustersFinal == 2) then
                if(.not. quiet) then
                    write(*, *) 'break 1'
                endif
                exit
            endif
           
            if (nClusters == 2) then
                if(.not. quiet) then
                    write(*, *) 'break 2'
                endif
                exit
            endif
           
            ! step 4.5
            nNei = nNei + delta
            if(.not. second) then
                if(nNei >= nObs-1) then
                    ! stop because nNei >=  nObs
                    ! we require 0<nNei<nObs
                   
                    if(.not. quiet) then
                        write(*, *) 'stop because nNei = ', nNei, & 
                            & '>= nObs = ', nObs          
                    endif
                    exit
                endif
            else
                if(nNei >= nNeiVec2(2)) then
                    ! stop because nNei > nNeiVec2(2)
                    if(.not. quiet) then
                        write(*, *) 'stop because nNei = ', nNei, & 
                            & '>= nNeiVec2(2) = ', nNeiVec2(2)          
                    endif
                    exit
                endif
            endif
        else if (myt == 1 .and. nClusters == 1) then
            nClustersFinal = 1
            indFinal = -999
            memFinal = 1
            if(indexFlag) then
                sFinal = -999
            endif
            clustSizeFinal = 0
            clustSizeFinal(1:nClusters) = nObs
            
            nNeiFinal = nNei
            if(.not. quiet) then
                write(*, *) 'break 3'
            endif
            exit
        else
            if(.not. quiet) then
                write(*, *) 'break 4'
            endif
            exit 
        endif
    end do
        
    myupdate = .FALSE.
 
    if(second .and. (indFinal > indInitial)) then
        !update final partition
        myupdate = .TRUE.
    endif
 
    ! step 5  output the clustering result
    if(.not. quiet) then
        write(*, *) 'final nNei = ', nNeiFinal, &
            & ' final nClusters = ', nClustersFinal
    endif

end subroutine chooseK
  
