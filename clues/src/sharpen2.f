C      program test
C      implicit none
C      integer nObs, nVars, nNei, nNei2, ITMAX, i, j, k, IX
C      integer KLIST(5), u, disMethod2
C      real*8 dat(384, 40), eps, datnew(384,40), DK(5)
C      real*8 obsi(40)
C
C      nObs=384
C      nVars=40
C      nNei=4
C      nNei2=nNei+1
C      ITMAX=200
C      eps=0.0001
C      disMethod2=2
C
C      u=99
C      open(UNIT=u, STATUS='OLD', FILE='dyeswap.dat')
C      do 10 i=1, nObs
C        read(u,*) (dat(i,j), j=1,nVars)
C10    continue        
C      do 20 j=1, nVars
C        obsi(j)=dat(2,j)
C20    continue        
C
CC      write(*,*) 'dat>>'
CC      do 30 i=1, nObs
CC        write(*,*) (dat(i,j), j=1,nVars)
CC30    continue        
CC      write(*,*) 'obsi>>'
CC      write(*,*) (obsi(j), j=1,nVars)
C      write(*,*) 'nObs=', nObs, ' nVars=', nVars, ' nNei=', nNei
C
C      call DIST(dat, obsi, nObs, nVars, nObs, 
C     *    nNei2, disMethod2, KLIST, DK)
C
C      write(*,*) 'KLIST>>', (KLIST(k),k=1,nNei2)
C      write(*,*) 'DK>>', (DK(k),k=1,nNei2)
C
C
C      call sharpen(dat,nObs,nVars,nNei2,nNei,ITMAX,eps, 
C     *   disMethod2, datnew)
C
C      write(*,*) 'datnew>>'
CC      do 60 i=1, nObs
C      do 60 i=2, 2
C        write(*,*) (datnew(i,j), j=1,nVars)
C60    continue        
C      stop
C      end
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
C      nNei2=nNei+1
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
      SUBROUTINE sharpen(dat,nObs,nVars,nNei2,nNei,ITMAX,
     *  eps, disMethod2, datnew)
      implicit none
      integer nObs, nVars, nNei, iter, ITMAX, nNei2
      integer i, j, IX, posj
      integer KLIST(nNei2), disMethod2
      real*8 DK(nNei2), varIX(nNei), tmp, maxdiff
      real*8 dat(nObs, nVars), datold(nObs, nVars)
      real*8 datnew(nObs, nVars), obsi(nVars)
      real*8 eps 

      nNei=nNei2-1
CCCCCC step 1. backup data matrix
      do 10 i = 1, nObs
        do 20 j = 1, nVars
          datold(i,j)=dat(i,j)
          datnew(i,j)=dat(i,j)
20      continue          
10    continue          

CCCCCC  initialize the iteration number
      iter=1;

CCCCCC step 2      

30    if (iter .le. ITMAX) then 
        maxdiff= -1.0
        do 40 i=1, nObs
CCCCCCCC  find K=nNei nearest neighbors of the data point i         
CCCCCCCCC get i-th observation
          do 60 IX = 1, nVars
            obsi(ix) = datold(i,ix)
60        continue

CCCCCCC   the K=nNei nearest neighbors of 'obsi' are stored in 'KLIST'
CCCCCCC   (not include itsself)
          CALL DIST(datold,obsi,nObs,nVars,nObs, nNei2,
     *      disMethod2, KLIST,DK)
CCCCCCC   update obsi by using the coordinate-wise median of
CCCCCCC   the 'nNei' nearest neighbors          
          do 70 IX = 1, nVars
            do 80 j = 1, nNei
              posj=KLIST(j+1)
              varIX(j)=datold(posj,IX)
80          continue
            CALL DMEDIAN(varIX, nNei, 0, tmp)
            datnew(i,IX) = tmp
            tmp=dabs(datnew(i,IX)-datold(i,IX))
            if (tmp .gt. maxdiff ) then
              maxdiff = tmp
            endif
70        continue
40      continue
CCCCCCC            
CCCCCCC Step 3
C
CCCCCCC update datold
        do 90 i = 1, nObs
          do 100 j = 1, nVars
            datold(i,j)=datnew(i,j)
100       continue            
90      continue            

        if (maxdiff<eps) then
C          write(*,*) 'iter=', iter, ' maxdiff=', maxdiff, 
C     *      ', shrinking coverged'
          goto 200
        endif
        iter=iter+1
      goto 30
      endif

200   end

C-------------------------------------------------------------------------
CCCCCC  The subroutine 'DIST' is a revised version of 'DIST'
CCCCCC   at http://astro.u-strasbg.fr/~fmurtagh/mda-sw/knn.f      
C      
CCCCCC  This a brutal-force method. Hence it is not efficient.
C      
CCCCCCC Find the K nearest neighbors of 'obsi'.
CCCCCCC The positions of these neighbors are recorded in 'KLIST' 
CCCCCCC The corresponding distance between 'obsi' and these neighbors
CCCCCCC are recorded in 'DK'      

C  dat(N,M)    data set, where N is the number of observations (rows),
C                 M is the number of variables (columns)
C  obsi(M)     the i-th observation;    
C  K              number of nearest neighbours to consider;
C  KLIST(K), DK(K)   are used for storing the K NNs
C                 and their distances to the object under 
C                 consideration.
C
C     sometimes, we only interested in the first N2 rows of dat      
      SUBROUTINE DIST(dat,obsi,N,M, N2, K,
     *  disMethod2, KLIST,DK)
      implicit none
      integer N, M, K, ITRAIN, ILOC, ICOLS,  IILOC, IX
      integer KLIST(K), N2, disMethod2
      real*8 D, DK(K), tmpD
      real*8 dat(N,M), obsi(M)
C      real*8 obsj(M)
      real*8 sumx, sumx2, sumy, sumy2, sumxy, rho
      real*8 numer, denom1, denom2
C
C      write(*,*) 'N=', N, ' M=', M, ' K=', K
C     initialize        
      DO 1 IX  = 1, K
        KLIST(IX) = 0
        DK(IX) = 1.E+15
1     CONTINUE

      DO 80 ITRAIN = 1, N2
C       write(*,*) 'ITRAIN=', ITRAIN
C       write(*,*) 'dat(ITRAIN,)>>', (dat(ITRAIN,IX), IX=1, M)
C       write(*,*) 'obsi>>', (obsi(IX), IX=1, M)

C        do 10 ICOLS=1, M
C          obsj(ICOLS)=dat(ITRAIN,ICOLS)
C10      continue
  
        sumx = 0.0
        sumx2 = 0.0
        sumy = 0.0
        sumy2 = 0.0
        sumxy = 0.0
        DO 30 ICOLS = 1, M
          sumx = sumx + dat(ITRAIN, ICOLS)
          sumx2 = sumx2 + dat(ITRAIN, ICOLS)**2
          sumy = sumy + obsi(ICOLS)
          sumy2 = sumy2 + obsi(ICOLS)**2
          sumxy = sumxy + dat(ITRAIN,ICOLS)*obsi(ICOLS)
30      CONTINUE

CCCC    Euclidean distance
        if(disMethod2 .eq. 1) then
          D=sumx2+sumy2-2*sumxy
          D=dsqrt(D)

        else 
CCCC    1-correlation or 1-rankcorrelation
CC      rho=numer/sqrt(denom1*denom2), where
CC      numer=n*sum(x*y) - sum(x)*sum(y) 
CC      denom1=n*sum(x^2)-(sum(x))^2
CC      denom2=n*sum(y^2)-(sum(y))^2

          numer=M*sumxy-sumx*sumy
          denom1=M*sumx2-sumx**2
          denom2=M*sumy2-sumy**2
          rho=numer/sqrt(denom1*denom2)
          D=1-rho 
        endif
C       write(*,*) 'ITRAIN=', ITRAIN, ' D=', D

C
        DO 50 ILOC = 1, K
C         write(*,*) 'ILOC=', ILOC
          IF (D.LT.DK(ILOC)) THEN
C Insert at locn. ILOC and shift right in the 3 length-k lists we''re maintaining
            DO 40 IILOC = K, ILOC+1, -1
C             write(*,*) 'IILOC=', IILOC
C             Protective measure:
              IF (IILOC.LE.K) THEN
                 DK(IILOC)    = DK(IILOC-1)
                 KLIST(IILOC) = KLIST(IILOC-1)
              ENDIF
C             Have now freed up space at locn. ILOC
   40       CONTINUE
            DK(ILOC)    = D
            KLIST(ILOC) = ITRAIN
C            write(*,*) 'DK(ILOC)=', DK(ILOC)
C            write(*,*) 'KLIST(ILOC)=', KLIST(ILOC)
            GOTO 60
          ENDIF
   50   CONTINUE
   60 CONTINUE
   80 CONTINUE
C             
      RETURN
      END
      

