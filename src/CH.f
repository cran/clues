C/**********
C created by Weiliang Qiu on Jan. 17, 2009      
C  stwxq@channing.harvard.edu      
C*********/ 

C      program test
C      implicit none
C      integer nObs, nVars, i, j, disMethod2
C      integer u, mem(384), nClusters, clustSize(56)
C      real*8 dat(384, 40), CH
C
C      nObs=384
C      nVars=40
C      nClusters=56
C      disMethod2=2
C
C      u=99
C      open(UNIT=u, STATUS='OLD', FILE='dyeswap.dat')
C      do 10 i=1, nObs
C        read(u,*) (dat(i,j), j=1,nVars)
C10    continue        
C      close(UNIT=u)
C
C      u=99
C      open(UNIT=u, STATUS='OLD', FILE='dyeswap.mem')
C      do 12 i=1, nObs
C        read(u,*) mem(i)
C12    continue        
C      close(UNIT=u)
C
C      u=99
C      open(UNIT=u, STATUS='OLD', FILE='dyeswap.size')
C      do 15 i=1, nClusters
C        read(u,*) clustSize(i)
C15    continue        
C      close(UNIT=u)
C
CC      write(*,*) 'dat>>'
CC      do 30 i=1, nObs
CC        write(*,*) (dat(i,j), j=1,nVars)
CC30    continue        
C      write(*,*) 'nObs=', nObs, ' nVars=', nVars
C
C      write(*,*) 'mem>>'
C      write(*,*) (mem(i), i=1, nObs)
C
C      write(*,*) 'clustSize>>'
C      write(*,*) (clustSize(i), i=1, nClusters)
C
C
C      call CHindex(dat,nObs,nVars, mem, nClusters,
C     *   clustSize, disMethod2, CH)
C
C      write(*,*) '************* output results *************'
C      write(*,*) 'CH>>', CH
C
C      stop
C      end



C@book{Kaufman:1990,
C  author = "{Kaufman, L.} and {Rousseeuw, P. J.}",
C  title = {Finding Groups in Data: An Introduction to Cluster Analysis},
C  publisher = {Wiley, New York.},
C  year = {1990},
C}
C**********/ 
C 
C/******
C For each observation i, the silhouette width s(i) is defined as
C     follows:       
C     Put a(i) = average dissimilarity between i and all other points of
C     the cluster to which i belongs.  For all other clusters C, put
C     d(i,C) = average dissimilarity of i to all observations of C.  The
C     smallest of these d(i,C) is b(i) := min_C d(i,C), and can be seen
C     as the dissimilarity between i and its ``neighbor'' cluster, i.e.,
C     the nearest one to which it does not belong. Finally,
C
C             s(i) := ( b(i) - a(i) ) / max( a(i), b(i) ).
C  
C     Observations with a large s(i) (almost 1) are very well clustered,
C     a small s(i) (around 0) means that the observation lies between
C     two clusters, and observations with a negative s(i) are probably
C     placed in the wrong cluster.
C --- the above is copied from R document on the function silhoutte() in
C library(cluster)
C*******/
C    
      
C    input:
C      mem: cluster each data point belonging to
C    output:
C      memNei: nearest neighbor cluster of each data point      
C      sIndex: silhoutte index for each data point
C      avgs: average silhoutte index

      SUBROUTINE CHindex(dat, nObs, nVars, mem, nClusters, 
     *  clustSize, disMethod2, CH)
      implicit none
      integer i, j, k, disMethod2
      integer nVars, nClusters, nObs
      integer mem(nObs), clustSize(nClusters)
      real*8 dat(nObs, nVars)
      real*8 mu(nClusters, nVars), bss, wss, muAll(nVars)
      real*8 tmp(nClusters), denom, CH
      real*8 sumx, sumx2, sumy, sumy2, sumxy, rho
      real*8 numer, denom1, denom2, D

      do 2 j=1, nVars
        do 3 k=1, nClusters
          mu(k,j)=0.0
3       continue
2     continue

C     calculate cluster mean 'mu' and overall mean 'muAll'
      do 10 k=1, nClusters
        do 20 j=1, nVars
          do 21 i=1, nObs
            if(mem(i) .eq. k) then
              mu(k,j)=mu(k,j)+dat(i,j)
            endif
21        continue
          mu(k,j)=mu(k,j)/clustSize(k)
20      continue
10    continue

C      write(*,*) 'mu>>'
C      do 25 k=1, nClusters
C        write(*,*) (mu(k,j), j=1, nVars)
C25    continue

      do 26 j=1, nVars
        muAll(j)=0.0
        do 27 i=1, nObs
          muAll(j)=muAll(j)+dat(i,j)
27      continue
        muAll(j)=muAll(j)/nObs 
26    continue
C      write(*,*) 'muAll>>', (muAll(j),j=1,nVars)

C     calculate the between-groups sum of squares and products matrix
C     for CH index


      bss=0.0
      do 70 k=1, nClusters
        sumx = 0.0
        sumx2 = 0.0
        sumy = 0.0
        sumy2 = 0.0
        sumxy = 0.0
        DO 30 j = 1, nVars
          sumx = sumx + mu(k, j)
          sumx2 = sumx2 + mu(k, j)**2
          sumy = sumy + muAll(j)
          sumy2 = sumy2 + muAll(j)**2
          sumxy = sumxy + mu(k, j)*muAll(j)
30      CONTINUE

CCCC    Euclidean distance
        if(disMethod2 .eq. 1) then
          D=sumx2+sumy2-2*sumxy
          D=dsqrt(D)
        else 
CCCC      1-correlation
CC        rho=numer/sqrt(denom1*denom2), where
CC        numer=n*sum(x*y) - sum(x)*sum(y) 
CC        denom1=n*sum(x^2)-(sum(x))^2
CC        denom2=n*sum(y^2)-(sum(y))^2
        
          numer=nVars*sumxy-sumx*sumy
          denom1=nVars*sumx2-sumx**2
          denom2=nVars*sumy2-sumy**2
          rho=numer/sqrt(denom1*denom2)
          D=1-rho 
        endif
        bss=bss+clustSize(k)*(D**2)
70    continue

C     calculate the within-groups sum of squares and products matrices
C     for CH index
      wss=0.0

      do 90 k=1, nClusters
        tmp(k)=0.0
        do 100 i=1, nObs
          if(mem(i) .eq. k) then
            sumx = 0.0
            sumx2 = 0.0
            sumy = 0.0
            sumy2 = 0.0
            sumxy = 0.0
            DO 110 j = 1, nVars
              sumx = sumx + dat(i, j)
              sumx2 = sumx2 + dat(i, j)**2
              sumy = sumy + mu(k,j)
              sumy2 = sumy2 + mu(k,j)**2
              sumxy = sumxy + dat(i, j)*mu(k,j)
110         CONTINUE
         
CCCC        Euclidean distance
            if(disMethod2 .eq. 1) then
              D=sumx2+sumy2-2*sumxy
              D=dsqrt(D)
            else 
CCCC          1-correlation
CC            rho=numer/sqrt(denom1*denom2), where
CC            numer=n*sum(x*y) - sum(x)*sum(y) 
CC            denom1=n*sum(x^2)-(sum(x))^2
CC            denom2=n*sum(y^2)-(sum(y))^2
            
              numer=nVars*sumxy-sumx*sumy
              denom1=nVars*sumx2-sumx**2
              denom2=nVars*sumy2-sumy**2
              rho=numer/sqrt(denom1*denom2)
              D=1-rho 
            endif
            tmp(k)=tmp(k)+D**2

C            do 110 j=1, nVars
C              tmp(k)=tmp(k)+(dat(i,j)-mu(k,j))**2
C110         continue
          endif
100     continue
        wss=wss+tmp(k)
90    continue
C      write(*,*) 'wss=', wss

      numer=bss/(nClusters-1.0)
      denom=wss/(nObs-nClusters)
C      write(*,*) 'numer=', numer, ' denom=', denom

C      write(*,*) 'bss=', bss, ' wss=', wss, 'nClusters=',
C     *  nClusters, ' nObs=', nObs 
      CH=numer/denom

      end

