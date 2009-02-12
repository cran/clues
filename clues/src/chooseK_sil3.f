C/**********
C created by Weiliang Qiu on Jan. 15, 2009      
C  stwxq@channing.harvard.edu      
C*********/ 

C      program test
C      implicit none
C      integer nObs, nObs1, nVars, i, j, myt
C      integer u, mem(384), nClusters, clustSize(2)
C      integer memNei(384), nClusters0
C      integer ITMAX, memFinal(384), disMethod2
C      integer nClustersFinal, clustSizeFinal(384)
C      integer nNeiFinal, nNeiVec2(2), nClustVec(384), nNeiVec(384)
C      real*8 dat(384, 40), dat2(384, 40), avgs, alpha, eps
C      real*8 sIndex(384), avgsFinal, sFinal(384)
C      real*8 datold1(384,40), datold2(384,40), s2
C      logical myupdate, second, quiet
C
C      nObs=384
C      nObs1=nObs-1
C      nVars=40
C      nClusters=3
CC      nClusters0=5
C      nClusters0=3
C      alpha=0.05
C      eps=1.0E-4
C      ITMAX=200
C      second=.FALSE.
C      nNeiVec2(1)=5
C      nNeiVec2(2)=5
C      s2=-3
C      disMethod2=2
C      quiet=.FALSE.
C
C      u=99
C      open(UNIT=u, STATUS='OLD', FILE='dyeswap.dat')
C      do 10 i=1, nObs
C        read(u,*) (dat(i,j), j=1,nVars)
C10    continue        
C      close(UNIT=u)
C
CC      write(*,*) 'dat>>'
CC      do 30 i=1, nObs
CC        write(*,*) (dat(i,j), j=1,nVars)
CC30    continue        
C      write(*,*) 'nObs=', nObs, ' nVars=', nVars
C
C        do 35 i=1, nObs
C          do 36 j=1, nVars
C            dat2(i,j)=dat(i,j)
C36        continue
C35      continue        
C
CC
C      call chooseKsil(dat, dat2, nObs, nObs1, nVars, nClusters0, alpha, 
C     *  eps,
C     *  ITMAX, second, nNeiVec2, s2, disMethod2, quiet,avgsFinal, 
C     *  sFinal,
C     *  memFinal, nClustersFinal, clustSizeFinal, nNeiFinal, 
C     *  nClustVec, nNeiVec, myt, datold1, datold2, myupdate)
C
C
C      write(*,*) '************* output results *************'
C      write(*,*) 'avgsFinal>>', avgsFinal
C      write(*,*) 'nClustersFinal>>', nClustersFinal
C      write(*,*) 'nNeiFinal>>', nNeiFinal
CC      write(*,*) 'sFinal>>', (sFinal(i), i=1,nObs)
C      write(*,*) 'memFinal>>', (memFinal(i), i=1,nObs)
C      write(*,*) 'clustSizeFinal>>', (clustSizeFinal(i), i=1,
C     *           nClustersFinal)
C      write(*,*) 'myupdate>>', myupdate
C      write(*,*) 'myt>>', myt
C      write(*,*) 'nClustVec>>', (nClustVec(i),i=1,myt)
C      write(*,*) 'nNeiVec>>', (nNeiVec(i),i=1,myt)
C
CC      write(*,*) 'datold1>>'
CC      do 40 i=1, nObs
CC        write(*,*) (datold1(i,j), j=1,nVars)
CC40    continue        
CC
CC      write(*,*) 'datold2>>'
CC      do 50 i=1, nObs
CC        write(*,*) (datold2(i,j), j=1,nVars)
CC50    continue        
CC
CC
C      stop
C      end
C

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C    choose an appropriate K of K nearest neighbors based on
C      silhouette index
C
C    input:
C      nClusters0: initial guess of the number of clusters
C      alpha: proportion of data as the K of K nearest neighbors
C      eps: used for the subroutine 'sharpen' 
C      ITMAX: maximum iteration number used for 'sharpen'
C      second: indicating if the second pass of the iteration is required  
C      nNeiVec2: nNeiVec2(1) and nNeiVec2(2) are the lower and upper bounds for 
C        the number of nearest neighbors for the second pass of the iteration.
C      s2: initial value of the lower bound of the silhouette index

C    output:
C      memNei: nearest neighbor cluster of each data point      
C      sIndex: silhoutte index for each data point
C      avgs: average silhoutte index

      SUBROUTINE chooseKsil(dat, dat2, nObs, nObs1, nVars, nClusters0, 
     *  alpha, eps,
     *  ITMAX, second, nNeiVec2, s2, disMethod2, quiet, 
     *  avgsFinal, sFinal,
     *  memFinal, nClustersFinal, clustSizeFinal, nNeiFinal, 
     *  nClustVec, nNeiVec, myt, datold1, datold2, myupdate)
      implicit none
      integer i, j, k, disMethod2
      integer nClusters0, ITMAX, nNeiVec2(2), delta, myt, myt2
      integer nObs, nObs1, nVars, nClusters
      integer mem(nObs), clustSize(nObs), memNei(nObs)
      integer nNei, nNei2, clustSizeFinal(nObs)
      integer memFinal(nObs), nNeiFinal, L1, L2
      integer nClustersFinal, minsize, nClustVec(nObs)
      integer nNeiVec(nObs), points(nObs)
      real*8 alpha, eps, s2, avgsFinal
      real*8 dat(nObs, nVars), dat2(nObs, nVars), sIndex(nObs), dik
      real*8 obsi(nObs), dij, avgs, a(nObs), b(nObs)
      real*8 datnew(nObs, nVars), datold1(nObs, nVars) 
      real*8 datold2(nObs, nVars), datnew2(nObs, nVars)
      real*8 db(nObs1), omin
      real*8 sFinal(nObs)
      logical myupdate, second, quiet

C     step 1
      do 10 i=1, nObs
        do 15 j=1, nVars
          datnew(i,j)=dat2(i,j)
          datold1(i,j)=dat2(i,j)
          datold2(i,j)=dat2(i,j)
15      continue
10    continue        

C     initialize the number K of the nearest neighbors.
C     and initialize the increment delta of the number of 
C     the nearest neighbors
      if (second .eqv. .FALSE.) then
        nNei=DNINT(alpha*nObs)
        delta=max(nNei, nClusters0)
      else
        nNei=nNeiVec2(1)
        i=DNINT(0.1*alpha*nObs)
        delta=max(i, 1)
        if(quiet .eqv. .FALSE.) then
          write(*,*) 'second=TRUE, nNei=', nNei, ' i=', i,
     *      ' delta=', delta
        endif
      endif 

C     initialize the iteration number
      myt=0
      avgsFinal= -2.0
      myt2=0

      L1=1
      L2=2
20    if (L2 .gt. L1) then
        myt=myt+1
        nNeiVec(myt)=nNei

C       Step 2

C       shrinking
        nNei2=nNei+1
        call sharpen(datnew,nObs,nVars,nNei2,nNei,ITMAX,eps, 
     *    disMethod2, datnew2)
        do 21 i=1, nObs
          do 22 j=1, nVars
            datnew(i,j)=datnew2(i,j)
22        continue
21      continue        


C       clustering
        call clustering(datnew,nObs,nObs1, nVars,disMethod2,
     *    points, db, omin, nClusters,mem, clustSize)
        nClustVec(myt)=nClusters

        if(nClusters .gt. 1) then
C         use original data to calculate sihouette index
          call silhouette(dat,nObs,nVars, mem, nClusters,
     *      clustSize, disMethod2, memNei, sIndex, avgs)

          minsize= nObs+1
          do 35 j=1, nClusters
            if(clustSize(j) .lt. minsize ) then
              minsize=clustSize(j)
            endif
35        continue
          if(quiet .eqv. .FALSE.) then
            write(*,*) 'myt=', myt, ' nNei=', nNei, 
     *        ' delta=', delta, ' minsize=', minsize,
     *        ' avgs=', avgs, ' nClusters=', nClusters
          endif

          if (avgs .gt. avgsFinal) then 
            avgsFinal=avgs
            do 40 i=1, nObs
              sFinal(i)=sIndex(i)
              memFinal(i)=mem(i)
              if (i .le. nClusters) then
                clustSizeFinal(i)=clustSize(i)
              else
                clustSizeFinal(i)=0
              endif
40          continue
            nClustersFinal=nClusters
            nNeiFinal=nNei
            myt2=myt2+1
          endif

          if (myt2 .eq. 2) then
            myt2=1
            do 50 i=1, nObs
              do 60 j=1, nVars
                datold1(i,j)=datold2(i,j)
60            continue
50          continue
          endif

          if (avgs .ge. avgsFinal) then
            do 70 i=1, nObs
              do 80 j=1, nVars
                datold2(i,j)=datnew(i,j)
80            continue
70          continue
          endif

          if (nClustersFinal .eq. 2) then
C           go to step 5
            if(quiet .eqv. .FALSE.) then
              write(*,*) 'break 1'
            endif
            goto 100
          endif

          if (nClusters .eq. 2) then
C           go to step 5
            if(quiet .eqv. .FALSE.) then
              write(*,*) 'break 2'
            endif

            goto 100
          endif

C         step 4
          nNei=nNei+delta
          if(second .eqv. .FALSE.) then
            if(nNei .ge. nObs-1) then
CCC           stop because nNei >= nObs
CCC           we require 0<nNei<nObs

              if(quiet .eqv. .FALSE.) then
                write(*,*) 'stop because nNei=', nNei, 
     *            '>=nObs=', nObs          
              endif
              goto 100
            endif
          else
            if(nNei .ge. nNeiVec2(2)) then
CCC           stop because nNei > nNeiVec2(2)
              if(quiet .eqv. .FALSE.) then
                write(*,*) 'stop because nNei=', nNei, 
     *            '>=nNeiVec2(2)=', nNeiVec2(2)          
              endif

              goto 100
            endif
          endif
        else if (myt .eq. 1 .and. nClusters .eq. 1) then
          nClustersFinal=1
          avgsFinal=-999
          do 90 i=1, nObs
            sFinal(i)=-999
            memFinal(i)=1
            if (i .le. nClusters) then
              clustSizeFinal(i)=nObs
            else
              clustSizeFinal(i)=0
            endif
90        continue
          nNeiFinal=nNei
          if(quiet .eqv. .FALSE.) then
C            write(*,*) 'myt=', myt, ' nNei=', nNei, 
C     *        ' nClusters=', nClusters
            write(*,*) 'break 3'
          endif
          goto 100
        else
          if(quiet .eqv. .FALSE.) then
C            write(*,*) 'myt=', myt, ' nNei=', nNei, 
C     *        ' nClusters=', nClusters
            write(*,*) 'break 4'
          endif

          goto 100
        endif

        goto  20
      endif
      
100   myupdate=.FALSE.

      if(second .eqv. .TRUE. .and. avgsFinal .gt. s2) then
CC      update final partition
        myupdate=.TRUE.
      endif

CC    step 5
CC    output the clustering result
      if(quiet .eqv. .FALSE.) then
        write(*,*) 'final nNei=', nNeiFinal,
     *    ' final nClusters=', nClustersFinal
      endif

      end


