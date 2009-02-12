C      program test
C      implicit none
C      integer nObs, nObs1, nVars, i, j, disMethod2
C      integer u, mem(384), nClusters, clustSize(384), points(384)
C      real*8 dat(384, 40), db(383), omin
C
C      nObs=384
C      nObs1=nObs-1
C      nVars=40
C      disMethod2=2
C
C      u=99
C      open(UNIT=u, STATUS='OLD', FILE='test_sharpened_dye.dat')
C      do 10 i=1, nObs
C        read(u,*) (dat(i,j), j=1,nVars)
C10    continue        
C
CC      write(*,*) 'dat>>'
CC      do 30 i=1, nObs
CC        write(*,*) (dat(i,j), j=1,nVars)
CC30    continue        
C      write(*,*) 'nObs=', nObs, ' nVars=', nVars
C
C      call clustering(dat,nObs,nObs1, nVars,disMethod2,
C     *   points, db, omin, nClusters,mem, clustSize)
C
C      write(*,*) '************* output results *************'
C      write(*,*) 'points>>', (points(i), i=1,nObs)
C      write(*,*) 'db>>', (db(i), i=1,nObs1)
C      write(*,*) 'omin=', omin
C      write(*,*) 'nClusters=', nClusters
C      write(*,*) 'mem>>'
C      write(*,*) (mem(i), i=1,nObs)
C      write(*,*) 'clustSize>>'
C      write(*,*) (clustSize(i), i=1,nClusters)
C
C      stop
C      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCC  nObs1=nObs-1
C
C we can plot db against {1, 2, ..., nObs1}
C   to check the distance between consecutive data points      
C      
C
C outputs:
C    points -- INTEGER. a sequence of labels of data points
C    db -- real*8. distance between consecutive points indicated by 'points'
C    omin -- real*8 if db(k) > omin, then the data point points(k) and
C      points(k+1) will be in two different clusters. That is,
C      a new cluster will start with points(k+1)      
C    nClusters -- INTEGER. number of clusters      
C    mem -- INTEGER. cluster membership      
C    clustSize -- INTEGER. size of clusters (i.e., number of data points
C      in each cluster      
C     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE clustering(dat,nObs,nObs1, nVars, disMethod2,
     *                      points, db, omin, nClusters, mem, clustSize)
      implicit none
      integer nObs, nObs1, nVars, nSstar, nNei, nNei2
      integer Sstar(nObs), disMethod2
      integer i, j, myt, loop, sumS, k
      integer m1, m2, nOutliers, nClusters, mem(nObs)
      integer KLIST(2), ell
      integer S(nObs), points(nObs)
      integer clustSize(nObs), pos1, pos2, L1, L2
      real*8 dat(nObs, nVars), obsi(nVars), dbk, upp, low
      real*8 DK(2), q1, q3, IQR
      real*8 omin, db(nObs1), db2(nObs1)
      real*8 dat2(nObs, nVars)

      nNei=1
      nNei2=nNei+1
CCCC  Step 1. get pair-wise distances dij among shrinked data points
      sumS=0
      do 10 i=1, nObs
        S(i)=i
        sumS=sumS+S(i)
        points(i)=0
        if (i<nObs) then
          db(i)=i
        endif
10    continue

CCCCC Step 2
      myt=1
      points(myt)=1
      S(myt)=0

CCCCC Step 3
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
      loop=0
20    if(sumS .gt. 0) then
        loop=loop+1
C        write(*,*) 'loop=', loop
C        write(*,*) 'S>>', (S(j), j=1, nObs)
C        write(*,*) 'points>>', (points(j), j=1, nObs)

        nSstar=0
        do 30 i=1, nObs
          if (S(i) .gt. 0) then
            nSstar=nSstar+1
C           set of remaining points            
            Sstar(nSstar)=S(i)
          endif
30      continue
C        write(*,*) 'nSstar=', nSstar,' Sstar>>',(Sstar(j), j=1, nSstar)

C    find the nearest neighbor x_{\ell} \in S* of the data point x_{point_t}
        i=points(myt) 
        do 40 j=1, nVars
          obsi(j)=dat(i,j)
40      continue
C        write(*,*) 'obsi>>', (obsi(j), j=1, nVars)
C       data points only in Sstar        
C        write(*,*) 'dat2>>'
        do 42 k=1, nSstar
          m1=Sstar(k)
          do 43 j=1, nVars
            dat2(k,j)=dat(m1, j)
43        continue            
C          write(*,*) (dat2(k,j), j=1, nVars)
42      continue        
        

C       find 'nNei2' nearest neighbor of obsi        
        CALL DIST(dat2,obsi,nObs,nVars,nSstar, nNei2,disMethod2, 
     *    KLIST,DK)
        ell=KLIST(1)
        ell=Sstar(ell)
        db(myt)=DK(1)
C        write(*,*) 'myt=', myt
C        write(*,*) 'KLIST>>', (KLIST(i), i=1, 2)
C        write(*,*) 'DK>>', (DK(i), i=1, 2)
C        write(*, *) 'ell=', ell, ' DK(1)=', DK(1), ' db(myt)=', db(myt)
C
        myt=myt+1
        points(myt)=ell
C        write(*,*) 'loop=', loop
        S(ell)=0

        sumS=0
        do 60 i=1, nObs
          sumS=sumS+S(i)
60      continue
C        if (myt .eq. 23) then
C          stop
C        endif
 
CCCCC goto Step 3
        goto 20
      endif

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
C      write(*,*) 'db>>', (db(i),i=1,nObs1)
C      write(*,*) 'points>>', (points(i),i=1,nObs)
C      write(*,*) 'S>>', (S(i),i=1,nObs)
C
CCCCC Step 4
C     find the outliers of db.
C     calculate sample quantiles
C     first sort 'db'
      CALL DSORT(db, nObs1, db2)
C      write(*,*) 'db2>>', (db2(j), j=1,nObs1)
C     25% percentile
      m1=mod(nObs1, 4)
      if(m1 .gt. 0) then
C       round up
        pos1=nObs1/4+1
      else
        pos1=nObs1/4
      endif
      q1=db2(pos1)

C     75% percentile
      m1=3*nObs1
      m2=mod(m1, 4) 
      if(m2 .gt. 0) then
C       round up
        pos2=m1/4+1
      else
        pos2=m1/4
      endif
      q3=db2(pos2)

C     interquantile range
      IQR=q3-q1

      upp=q3+1.5*IQR
      low=max(0.0, q1-1.5*IQR)
C      write(*,*) 'q1=', q1, ' q3=', q3, ' IQR=', IQR
C      write(*,*) 'low=', low, ' upp=', upp

CC    outliers are defined as db(k)>q3+1.5*IQR or db(k)<q1-1.5*IQR      
      nOutliers=0
C     minimum of outlier distances 
      omin=1.0E+30
      do 70 k=1,nObs1
        dbk=db(k)
        if((dbk .gt. upp) .or. (dbk .lt. low)) then
          nOutliers=nOutliers+1
          if (omin .gt. dbk) then
            omin=dbk
          endif
        endif
70    continue

CCCC  Step 5
      if (nOutliers .eq. 0) then
        do 80 k=1, nObs
          mem(k)=1
80      continue        
        nClusters=1
      else
C       Step 6
C       initialize mem and nClusters
        do 85 k=1, nObs
          mem(k)=0
85      continue        
        nClusters=1
        myt=1
        mem(myt)=1

        L1=1
        L2=2
90      if (L2 .gt. L1) then
C         Step 7
          if (db(myt) .ge. omin) then
            nClusters=nClusters+1
          endif
C         Step 8          
          m1=myt+1
          k=points(m1)
          mem(k)=nClusters
          myt=myt+1
C         Step 9
          if (myt .ge. nObs) then
C           go to Step 10                  
            goto 100
          endif
          goto 90
        endif

C       Step 10
C       calculate the cluster size        
100     do 101 k=1, nClusters
          clustSize(k)=0
101     continue          
        do 102 k=1, nClusters
          m1=0
          do 103 i=1, nObs
            if (mem(i) .eq. k) then
              m1=m1+1
            endif
103       continue            
          clustSize(k)=m1
102     continue        
      endif

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      end

 
