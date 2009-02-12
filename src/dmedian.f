C http://www.itl.nist.gov/div898/software/datapac/MEDIAN.f      
C  revise 'median' to 'dmedian' by weiliang Qiu on Jan. 14, 2009

C      program test
C      implicit none
C
C      integer N, IWRITE, j
C      real*8 X(3), XMED
C
C      N=3
C      IWRITE=1
C
C      X(1)=0.5128151  
C      X(2)=1.3557439  
C      X(3)=0.3444514
CC      X(4)=-0.7638553  
CC      X(5)=0.7372837 
CC      X(6)=-0.3866299
CC
CC      X(1)=-1.1479490  
CC      X(2)=0.5618702  
CC      X(3)=1.2772249 
CC      X(4)=-0.7126383 
CC      X(5)=-0.1055491  
CC      X(6)=1.8723037 
CC      X(7)=-0.8268037
C
C      write(*,*) (X(j),j=1,N)
C
C      call dmedian(X, N, IWRITE, XMED)
C      write(*,*) 'XMED=', XMED
C
C      stop
C      end
C

      SUBROUTINE DMEDIAN(X,N,IWRITE,XMED)
C
C     PURPOSE--THIS SUBROUTINE COMPUTES THE
C              SAMPLE MEDIAN
C              OF THE DATA IN THE INPUT VECTOR X. 
C              THE SAMPLE MEDIAN = THAT VALUE SUCH THAT HALF THE
C              DATA SET IS BELOW IT AND HALF ABOVE IT.
C     INPUT  ARGUMENTS--X      = THE double PRECISION VECTOR OF
C                                (UNSORTED OR SORTED) OBSERVATIONS.
C                     --N      = THE INTEGER NUMBER OF OBSERVATIONS
C                                IN THE VECTOR X. 
C                     --IWRITE = AN INTEGER FLAG CODE WHICH 
C                                (IF SET TO 0) WILL SUPPRESS
C                                THE PRINTING OF THE
C                                SAMPLE MEDIAN
C                                AS IT IS COMPUTED;
C                                OR (IF SET TO SOME INTEGER 
C                                VALUE NOT EQUAL TO 0),
C                                LIKE, SAY, 1) WILL CAUSE
C                                THE PRINTING OF THE
C                                SAMPLE MEDIAN
C                                AT THE TIME IT IS COMPUTED.
C     OUTPUT ARGUMENTS--XMED   = THE double PRECISION VALUE OF THE
C                                COMPUTED SAMPLE MEDIAN.
C     OUTPUT--THE COMPUTED double PRECISION VALUE OF THE
C             SAMPLE MEDIAN.
C     PRINTING--NONE, UNLESS IWRITE HAS BEEN SET TO A NON-ZERO
C               INTEGER, OR UNLESS AN INPUT ARGUMENT ERROR
C               CONDITION EXISTS.
C     RESTRICTIONS--THE MAXIMUM ALLOWABLE VALUE OF N
C                   FOR THIS SUBROUTINE IS 15000. 
C     OTHER DATAPAC   SUBROUTINES NEEDED--SORT.
C     FORTRAN LIBRARY SUBROUTINES NEEDED--NONE.
C     MODE OF INTERNAL OPERATIONS--double PRECISION.
C     LANGUAGE--ANSI FORTRAN. 
C     REFERENCES--KENDALL AND STUART, THE ADVANCED THEORY OF
C                 STATISTICS, VOLUME 1, EDITION 2, 1963, PAGE 326.
C               --KENDALL AND STUART, THE ADVANCED THEORY OF
C                 STATISTICS, VOLUME 2, EDITION 1, 1961, PAGE 49.
C               --DAVID, ORDER STATISTICS, 1970, PAGE 139.
C               --SNEDECOR AND COCHRAN, STATISTICAL METHODS,
C                 EDITION 6, 1967, PAGE 123.
C               --DIXON AND MASSEY, INTRODUCTION TO STATISTICAL
C                 ANALYSIS, EDITION 2, 1957, PAGE 70.
C     WRITTEN BY--JAMES J. FILLIBEN
C                 STATISTICAL ENGINEERING LABORATORY (205.03)
C                 NATIONAL BUREAU OF STANDARDS
C                 WASHINGTON, D. C. 20234
C                 PHONE:  301-921-2315
C     ORIGINAL VERSION--JUNE      1972. 
C     UPDATED         --SEPTEMBER 1975. 
C     UPDATED         --NOVEMBER  1975. 
C     UPDATED         --FEBRUARY  1976. 
C
C---------------------------------------------------------------------
C
       integer N, IWRITE, IPR, IUPPER, I, IFLAG, NMID
       integer NMIDP1
       real*8 X(N), XMED, HOLD
       real*8 Y(15000) 
       COMMON /BLOCK2/ WS(15000)
       EQUIVALENCE (Y(1),WS(1))
C      
       IPR=6
       IUPPER=15000
C      
C      CHECK THE INPUT ARGUMENTS FOR ERRORS
C      
       IF(N .LT. 1 .OR. N .GT. IUPPER) GOTO 50
       IF(N .EQ. 1) GOTO 55
       HOLD=X(1)
       DO 60 I=2,N
       IF(X(I) .NE. HOLD) GOTO 90
60     CONTINUE
C       WRITE(IPR, 9) "***** NON-FATAL DIAGNOSTIC--THE FIRST  INPUT 
C     1 ARGUMENT (A VECTOR) TO THE MEDIAN SUBROUTINE HAS ALL ELEMENTS=",
C     1 HOLD, "*****"
       XMED=X(1)
       GOTO 101
50     WRITE(IPR,17) "***** FATAL ERROR--THE SECOND INPUT ARGUMENT TO 
     1 THE MEDIAN SUBROUTINE IS OUTSIDE THE ALLOWABLE ", IUPPER,  
     1 " INTERVAL *****"
       WRITE(IPR,47) "***** THE VALUE OF THE ARGUMENT IS ", N, " *****"
       RETURN
C55     WRITE(IPR,18) "***** NON-FATAL DIAGNOSTIC--THE SECOND INPUT 
C     1 ARGUMENT TO THE MEDIAN SUBROUTINE HAS THE VALUE 1 *****"
C       XMED=X(1)
55       XMED=X(1)
       GOTO 101
90     CONTINUE
9      FORMAT(A109, E15.8, A6)
17     FORMAT(A109, I6)
18     FORMAT(A109)
47     FORMAT(A35,I8, A6) 
C
C-----START POINT-----------------------------------------------------
C
       CALL DSORT(X,N,Y)
       IFLAG=N-(N/2)*2
       NMID=N/2
       NMIDP1=NMID+1 
       IF(IFLAG .EQ. 0) XMED=(Y(NMID)+Y(NMIDP1))/2.0
       IF(IFLAG .EQ. 1) XMED=Y(NMIDP1)
C
101    IF(IWRITE .EQ. 0) RETURN
       WRITE(IPR,999)
       WRITE(IPR,105) N,XMED
105    FORMAT(1H ,25HTHE SAMPLE MEDIAN OF THE ,I6,17H OBSERVATIONS IS ,
     1 E1
     1 5.8)
999    FORMAT(1H )
       RETURN
       END 

