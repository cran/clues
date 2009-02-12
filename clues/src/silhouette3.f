C/**********
C created by Weiliang Qiu on Jan. 15, 2009      
C  stwxq@channing.harvard.edu      
C*********/ 

C      program test
C      implicit none
C      integer nObs, nObs1, nVars, i, j
C      integer u, mem(384), nClusters, clustSize(56)
C      integer memNei(384), disMethod2
C      real*8 dat(384, 40), avgs
C      real*8 sIndex(384)
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
C      call silhouette(dat,nObs,nVars, mem, nClusters,
C     *   clustSize, disMethod2, memNei, sIndex, avgs)
C
C      write(*,*) '************* output results *************'
C      write(*,*) 'memNei>>', (memNei(i), i=1,nObs)
C      write(*,*) 'sIndex>>', (sIndex(i), i=1,nObs)
C      write(*,*) 'avgs=', avgs
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

      SUBROUTINE SILHOUETTE(dat, nObs, nVars, mem, nClusters, 
     *  clustSize, disMethod2, memNei, sIndex, avgs)
      implicit none
      integer i, j, k, memi, pos, sizemi, memj, j2, sizemj
      integer nObs, nVars, nClusters, disMethod2
      integer mem(nObs), clustSize(nClusters), memNei(nObs)
      integer ICOLS, M
      real*8 dat(nObs, nVars), sIndex(nObs), dik
      real*8 obsi(nObs), dij, avgs, a(nObs), b(nObs)
      real*8 sumx, sumx2, sumy, sumy2, sumxy, rho
      real*8 numer, denom1, denom2, D

      M=nVars
C     initializing memNei
      do 10 i=1, nObs
        memNei(i)=0
10    continue        

      avgs=0.0
      do 20 i=1, nObs
C       obtain the i-th observation (row) of dat      
        do 22 j=1, nVars
          obsi(j)=dat(i,j)
22      continue
        memi=mem(i)
        sizemi=clustSize(memi)
        if(sizemi .eq. 1) then
C         only one data point in cluster 'memi', then a(i)=0        
          a(i)=0.0
        else
C         get the average distance from point i to other data points in the
C         same cluster as point i
          dik=0.0

          do 50 j=1, nObs
            if( (j .ne. i) .and. (mem(j) .eq. memi)) then
C             get distance between point i and point j                    

              sumx = 0.0
              sumx2 = 0.0
              sumy = 0.0
              sumy2 = 0.0
              sumxy = 0.0
              DO 25 ICOLS = 1, nVars
                sumx = sumx + dat(j, ICOLS)
                sumx2 = sumx2 + dat(j, ICOLS)**2
                sumy = sumy + obsi(ICOLS)
                sumy2 = sumy2 + obsi(ICOLS)**2
                sumxy = sumxy + dat(j, ICOLS)*obsi(ICOLS)
25            CONTINUE
             
CCCC          Euclidean distance
              if(disMethod2 .eq. 1) then
                D=sumx2+sumy2-2*sumxy
                D=dsqrt(D)
             
              else 
CCCC            1-correlation
CC              rho=numer/sqrt(denom1*denom2), where
CC              numer=n*sum(x*y) - sum(x)*sum(y) 
CC              denom1=n*sum(x^2)-(sum(x))^2
CC              denom2=n*sum(y^2)-(sum(y))^2
             
                numer=M*sumxy-sumx*sumy
                denom1=M*sumx2-sumx**2
                denom2=M*sumy2-sumy**2
                rho=numer/sqrt(denom1*denom2)
                D=1-rho 
              endif
              dij=D
              dik=dik+dij
            endif
50        continue         
          a(i)=dik/(sizemi-1.0)
        endif          
C        write(*,*) 'i=', i, ' a(i)=', a(i)

C       get the average distance from point i to data points in the
C       nearest neighbor cluster besides its own
        b(i)=1.0E+30
        do 60 j=1, nClusters
          memj=j
          if (memi .ne. memj) then
            dik=0.0
            do 70 j2=1, nObs
              if (mem(j2) .eq. memj) then
C               distance between point i and point j2 in cluster memj

                sumx = 0.0
                sumx2 = 0.0
                sumy = 0.0
                sumy2 = 0.0
                sumxy = 0.0
                DO 73 ICOLS = 1, nVars
                  sumx = sumx + dat(j2, ICOLS)
                  sumx2 = sumx2 + dat(j2, ICOLS)**2
                  sumy = sumy + obsi(ICOLS)
                  sumy2 = sumy2 + obsi(ICOLS)**2
                  sumxy = sumxy + dat(j2, ICOLS)*obsi(ICOLS)
73              CONTINUE
              
CCCC            Euclidean distance
                if(disMethod2 .eq. 1) then
                  D=sumx2+sumy2-2*sumxy
                  D=dsqrt(D)
              
                else 
CCCC              1-correlation
CC                rho=numer/sqrt(denom1*denom2), where
CC                numer=n*sum(x*y) - sum(x)*sum(y) 
CC                denom1=n*sum(x^2)-(sum(x))^2
CC                denom2=n*sum(y^2)-(sum(y))^2
              
                  numer=M*sumxy-sumx*sumy
                  denom1=M*sumx2-sumx**2
                  denom2=M*sumy2-sumy**2
                  rho=numer/sqrt(denom1*denom2)
                  D=1-rho 
                endif
                dij=D
                dik=dik+dij
              endif                      
70          continue            
C            write(*,*) 'j=', j, ' clustSize(j)=', clustSize(j)
            sizemj=clustSize(j)
C           dik=average distance of point i to all points in cluster memj
            dik = dik/(sizemj*1.0)
C            write(*,*) 'dik=', dik, ' b(i)=', b(i)
            if( b(i) .gt. dik ) then
              b(i) = dik
C             data point i's nearest neighor cluster is cluster 'memj'
              memNei(i)=memj
            endif
          endif
60      continue        
C        write(*,*) 'i=', i, ' b(i)=', b(i)
        sIndex(i)=(b(i) - a(i))/max(b(i), a(i))
C        write(*,*) 'i=', i, ' sIndex(i)=', sIndex(i)
        avgs=avgs+sIndex(i)

20    continue      

      avgs=avgs/(nObs*1.0)

      end

