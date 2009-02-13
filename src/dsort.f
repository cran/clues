C http://www.itl.nist.gov/div898/software/datapac/SORT.f      
C revise 'sort' to 'dsort' for double data by Weiliang Qiu on Jan 14,
C 2009
      SUBROUTINE DSORT(X,N,Y)
C
C     PURPOSE--THIS SUBROUTINE SORTS (IN ASCENDING ORDER)
C              THE N ELEMENTS OF THE double PRECISION VECTOR X
C              AND PUTS THE RESULTING N SORTED VALUES INTO THE
C              double PRECISION VECTOR Y.
C     INPUT  ARGUMENTS--X      = THE double PRECISION VECTOR OF
C                                OBSERVATIONS TO BE SORTED. 
C                     --N      = THE INTEGER NUMBER OF OBSERVATIONS
C                                IN THE VECTOR X. 
C     OUTPUT ARGUMENTS--Y      = THE double PRECISION VECTOR
C                                INTO WHICH THE SORTED DATA VALUES
C                                FROM X WILL BE PLACED.
C     OUTPUT--THE double PRECISION VECTOR Y
C             CONTAINING THE SORTED
C             (IN ASCENDING ORDER) VALUES
C             OF THE double PRECISION VECTOR X.
C     PRINTING--NONE UNLESS AN INPUT ARGUMENT ERROR CONDITION EXISTS. 
C     RESTRICTIONS--THE DIMENSIONS OF THE VECTORS IL AND IU 
C                   (DEFINED AND USED INTERNALLY WITHIN
C                   THIS SUBROUTINE) DICTATE THE MAXIMUM
C                   ALLOWABLE VALUE OF N FOR THIS SUBROUTINE.
C                   IF IL AND IU EACH HAVE DIMENSION K,
C                   THEN N MAY NOT EXCEED 2**(K+1) - 1.
C                   FOR THIS SUBROUTINE AS WRITTEN, THE DIMENSIONS
C                   OF IL AND IU HAVE BEEN SET TO 36,
C                   THUS THE MAXIMUM ALLOWABLE VALUE OF N IS
C                   APPROXIMATELY 137 BILLION.
C                   SINCE THIS EXCEEDS THE MAXIMUM ALLOWABLE
C                   VALUE FOR AN INTEGER VARIABLE IN MANY COMPUTERS,
C                   AND SINCE A SORT OF 137 BILLION ELEMENTS
C                   IS PRESENTLY IMPRACTICAL AND UNLIKELY,
C                   THEN THERE IS NO PRACTICAL RESTRICTION
C                   ON THE MAXIMUM VALUE OF N FOR THIS SUBROUTINE.
C                   (IN LIGHT OF THE ABOVE, NO CHECK OF THE 
C                   UPPER LIMIT OF N HAS BEEN INCORPORATED
C                   INTO THIS SUBROUTINE.)
C     OTHER DATAPAC   SUBROUTINES NEEDED--NONE.
C     FORTRAN LIBRARY SUBROUTINES NEEDED--NONE.
C     MODE OF INTERNAL OPERATIONS--double PRECISION.
C     LANGUAGE--ANSI FORTRAN. 
C     COMMENT--THE SMALLEST ELEMENT OF THE VECTOR X
C              WILL BE PLACED IN THE FIRST POSITION
C              OF THE VECTOR Y,
C              THE SECOND SMALLEST ELEMENT IN THE VECTOR X
C              WILL BE PLACED IN THE SECOND POSITION
C              OF THE VECTOR Y, ETC.
C     COMMENT--THE INPUT VECTOR X REMAINS UNALTERED.
C     COMMENT--IF THE ANALYST DESIRES A SORT 'IN PLACE',
C              THIS IS DONE BY HAVING THE SAME
C              OUTPUT VECTOR AS INPUT VECTOR IN THE CALLING SEQUENCE. 
C              THUS, FOR EXAMPLE, THE CALLING SEQUENCE
C              CALL SORT(X,N,X)
C              IS ALLOWABLE AND WILL RESULT IN
C              THE DESIRED 'IN-PLACE' SORT.
C     COMMENT--THE SORTING ALGORTHM USED HEREIN
C              IS THE BINARY SORT.
C              THIS ALGORTHIM IS EXTREMELY FAST AS THE
C              FOLLOWING TIME TRIALS INDICATE.
C              THESE TIME TRIALS WERE CARRIED OUT ON THE
C              UNIVAC 1108 EXEC 8 SYSTEM AT NBS
C              IN AUGUST OF 1974.
C              BY WAY OF COMPARISON, THE TIME TRIAL VALUES
C              FOR THE EASY-TO-PROGRAM BUT EXTREMELY
C              INEFFICIENT BUBBLE SORT ALGORITHM HAVE
C              ALSO BEEN INCLUDED--
C              NUMBER OF RANDOM        BINARY SORT       BUBBLE SORT
C               NUMBERS SORTED
C                N = 10                 .002 SEC          .002 SEC
C                N = 100                .011 SEC          .045 SEC
C                N = 1000               .141 SEC         4.332 SEC
C                N = 3000               .476 SEC        37.683 SEC
C                N = 10000             1.887 SEC      NOT COMPUTED
C     REFERENCES--CACM MARCH 1969, PAGE 186 (BINARY SORT ALGORITHM
C                 BY RICHARD C. SINGLETON).
C               --CACM JANUARY 1970, PAGE 54.
C               --CACM OCTOBER 1970, PAGE 624.
C               --JACM JANUARY 1961, PAGE 41.
C     WRITTEN BY--JAMES J. FILLIBEN
C                 STATISTICAL ENGINEERING LABORATORY (205.03)
C                 NATIONAL BUREAU OF STANDARDS
C                 WASHINGTON, D. C. 20234
C                 PHONE--301-921-2315
C     ORIGINAL VERSION--JUNE      1972. 
C     UPDATED         --NOVEMBER  1975. 
C
C---------------------------------------------------------------------
C
       implicit none
       integer IPR, N, I, NM1, IP1, J, L, K, LMI, JMK, M, JMI
       integer MID
       real*8 X(1),Y(1), HOLD, AMED, TT
       integer IU(36),IL(36) 
C      
       IPR=6
C      
C      CHECK THE INPUT ARGUMENTS FOR ERRORS
C      
       IF(N .LT. 1) GOTO 50
       IF(N .EQ. 1) GOTO 55
       HOLD=X(1)
       DO 60 I=2,N
         IF(X(I).NE.HOLD) GOTO 90
60     CONTINUE
       DO 61 I=1,N
         Y(I)=X(I)
61     CONTINUE
       RETURN
50     WRITE(IPR,15) "***** FATAL ERROR--THE SECOND INPUT ARGUMENT TO 
     1 THE SORT   SUBROUTINE IS NON-POSITIVE *****"  
       WRITE(IPR,47) "***** THE VALUE OF THE ARGUMENT IS", N, "****" 
       RETURN
55       Y(1)=X(1)
       RETURN
90     CONTINUE
15     FORMAT(A109)
47     FORMAT(A109,I8, A6)
C
C-----START POINT-----------------------------------------------------
C
C     COPY THE VECTOR X INTO THE VECTOR Y
      DO 100 I=1,N
        Y(I)=X(I)
100   CONTINUE
C
C     CHECK TO SEE IF THE INPUT VECTOR IS ALREADY SORTED
C
      NM1=N-1
      DO 200 I=1,NM1
        IP1=I+1
        IF(Y(I) .LE. Y(IP1)) GOTO 200
      GOTO 250
200   CONTINUE
      RETURN
250   M=1 
      I=1 
      J=N 
  305 IF(I.GE.J)GOTO370
  310 K=I 
      MID=(I+J)/2
      AMED=Y(MID)
      IF(Y(I) .LE. AMED) GOTO 320 
      Y(MID)=Y(I)
      Y(I)=AMED
      AMED=Y(MID)
320   L=J 
      IF(Y(J) .GE. AMED) GOTO 340 
      Y(MID)=Y(J)
      Y(J)=AMED
      AMED=Y(MID)
      IF(Y(I) .LE. AMED) GOTO 340 
      Y(MID)=Y(I)
      Y(I)=AMED
      AMED=Y(MID)
      GOTO 340
330   Y(L)=Y(K)
      Y(K)=TT
340   L=L-1
      IF(Y(L) .GT. AMED) GOTO 340 
      TT=Y(L)
350   K=K+1
      IF(Y(K) .LT. AMED) GOTO 350 
      IF(K .LE. L) GOTO 330
      LMI=L-I
      JMK=J-K
      IF(LMI .LE. JMK) GOTO 360
      IL(M)=I
      IU(M)=L
      I=K 
      M=M+1
      GOTO 380
360   IL(M)=K
      IU(M)=J
      J=L 
      M=M+1
      GOTO 380
370   M=M-1
      IF(M .EQ. 0) RETURN
      I=IL(M)
      J=IU(M)
380   JMI=J-I
      IF(JMI .GE. 11) GOTO 310
      IF(I .EQ. 1) GOTO 305
      I=I-1
390   I=I+1
      IF(I .EQ. J) GOTO 370
      AMED=Y(I+1)
      IF(Y(I) .LE. AMED) GOTO 390 
      K=I 
395   Y(K+1)=Y(K)
      K=K-1
      IF(AMED .LT. Y(K)) GOTO 395 
      Y(K+1)=AMED
      GOTO 390
      END 

