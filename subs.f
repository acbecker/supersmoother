      SUBROUTINE DQRSL(X,LDX,N,K,QRAUX,Y,QY,QTY,B,RSD,XB,JOB,INFO)
      INTEGER LDX,N,K,JOB,INFO
      DOUBLE PRECISION X(LDX,1),QRAUX(1),Y(1),QY(1),QTY(1),B(1),RSD(1),
     *                 XB(1)
C
C     DQRSL APPLIES THE OUTPUT OF DQRDC TO COMPUTE COORDINATE
C     TRANSFORMATIONS, PROJECTIONS, AND LEAST SQUARES SOLUTIONS.
C     FOR K .LE. MIN(N,P), LET XK BE THE MATRIX
C
C            XK = (X(JPVT(1)),X(JPVT(2)), ... ,X(JPVT(K)))
C
C     FORMED FROM COLUMNNS JPVT(1), ... ,JPVT(K) OF THE ORIGINAL
C     N X P MATRIX X THAT WAS INPUT TO DQRDC (IF NO PIVOTING WAS
C     DONE, XK CONSISTS OF THE FIRST K COLUMNS OF X IN THEIR
C     ORIGINAL ORDER).  DQRDC PRODUCES A FACTORED ORTHOGONAL MATRIX Q
C     AND AN UPPER TRIANGULAR MATRIX R SUCH THAT
C
C              XK = Q * (R)
C                       (0)
C
C     THIS INFORMATION IS CONTAINED IN CODED FORM IN THE ARRAYS
C     X AND QRAUX.
C
C     ON ENTRY
C
C        X      DOUBLE PRECISION(LDX,P).
C               X CONTAINS THE OUTPUT OF DQRDC.
C
C        LDX    INTEGER.
C               LDX IS THE LEADING DIMENSION OF THE ARRAY X.
C
C        N      INTEGER.
C               N IS THE NUMBER OF ROWS OF THE MATRIX XK.  IT MUST
C               HAVE THE SAME VALUE AS N IN DQRDC.
C
C        K      INTEGER.
C               K IS THE NUMBER OF COLUMNS OF THE MATRIX XK.  K
C               MUST NNOT BE GREATER THAN MIN(N,P), WHERE P IS THE
C               SAME AS IN THE CALLING SEQUENCE TO DQRDC.
C
C        QRAUX  DOUBLE PRECISION(P).
C               QRAUX CONTAINS THE AUXILIARY OUTPUT FROM DQRDC.
C
C        Y      DOUBLE PRECISION(N)
C               Y CONTAINS AN N-VECTOR THAT IS TO BE MANIPULATED
C               BY DQRSL.
C
C        JOB    INTEGER.
C               JOB SPECIFIES WHAT IS TO BE COMPUTED.  JOB HAS
C               THE DECIMAL EXPANSION ABCDE, WITH THE FOLLOWING
C               MEANING.
C
C                    IF A.NE.0, COMPUTE QY.
C                    IF B,C,D, OR E .NE. 0, COMPUTE QTY.
C                    IF C.NE.0, COMPUTE B.
C                    IF D.NE.0, COMPUTE RSD.
C                    IF E.NE.0, COMPUTE XB.
C
C               NOTE THAT A REQUEST TO COMPUTE B, RSD, OR XB
C               AUTOMATICALLY TRIGGERS THE COMPUTATION OF QTY, FOR
C               WHICH AN ARRAY MUST BE PROVIDED IN THE CALLING
C               SEQUENCE.
C
C     ON RETURN
C
C        QY     DOUBLE PRECISION(N).
C               QY CONNTAINS Q*Y, IF ITS COMPUTATION HAS BEEN
C               REQUESTED.
C
C        QTY    DOUBLE PRECISION(N).
C               QTY CONTAINS TRANS(Q)*Y, IF ITS COMPUTATION HAS
C               BEEN REQUESTED.  HERE TRANS(Q) IS THE
C               TRANSPOSE OF THE MATRIX Q.
C
C        B      DOUBLE PRECISION(K)
C               B CONTAINS THE SOLUTION OF THE LEAST SQUARES PROBLEM
C
C                    MINIMIZE NORM2(Y - XK*B),
C
C               IF ITS COMPUTATION HAS BEEN REQUESTED.  (NOTE THAT
C               IF PIVOTING WAS REQUESTED IN DQRDC, THE J-TH
C               COMPONENT OF B WILL BE ASSOCIATED WITH COLUMN JPVT(J)
C               OF THE ORIGINAL MATRIX X THAT WAS INPUT INTO DQRDC.)
C
C        RSD    DOUBLE PRECISION(N).
C               RSD CONTAINS THE LEAST SQUARES RESIDUAL Y - XK*B,
C               IF ITS COMPUTATION HAS BEEN REQUESTED.  RSD IS
C               ALSO THE ORTHOGONAL PROJECTION OF Y ONTO THE
C               ORTHOGONAL COMPLEMENT OF THE COLUMN SPACE OF XK.
C
C        XB     DOUBLE PRECISION(N).
C               XB CONTAINS THE LEAST SQUARES APPROXIMATION XK*B,
C               IF ITS COMPUTATION HAS BEEN REQUESTED.  XB IS ALSO
C               THE ORTHOGONAL PROJECTION OF Y ONTO THE COLUMN SPACE
C               OF X.
C
C        INFO   INTEGER.
C               INFO IS ZERO UNLESS THE COMPUTATION OF B HAS
C               BEEN REQUESTED AND R IS EXACTLY SINGULAR.  IN
C               THIS CASE, INFO IS THE INDEX OF THE FIRST ZERO
C               DIAGONAL ELEMENT OF R AND B IS LEFT UNALTERED.
C
C     THE PARAMETERS QY, QTY, B, RSD, AND XB ARE NOT REFERENCED
C     IF THEIR COMPUTATION IS NOT REQUESTED AND IN THIS CASE
C     CAN BE REPLACED BY DUMMY VARIABLES IN THE CALLING PROGRAM.
C     TO SAVE STORAGE, THE USER MAY IN SOME CASES USE THE SAME
C     ARRAY FOR DIFFERENT PARAMETERS IN THE CALLING SEQUENCE.  A
C     FREQUENTLY OCCURING EXAMPLE IS WHEN ONE WISHES TO COMPUTE
C     ANY OF B, RSD, OR XB AND DOES NOT NEED Y OR QTY.  IN THIS
C     CASE ONE MAY IDENTIFY Y, QTY, AND ONE OF B, RSD, OR XB, WHILE
C     PROVIDING SEPARATE ARRAYS FOR ANYTHING ELSE THAT IS TO BE
C     COMPUTED.  THUS THE CALLING SEQUENCE
C
C          CALL DQRSL(X,LDX,N,K,QRAUX,Y,DUM,Y,B,Y,DUM,110,INFO)
C
C     WILL RESULT IN THE COMPUTATION OF B AND RSD, WITH RSD
C     OVERWRITING Y.  MORE GENERALLY, EACH ITEM IN THE FOLLOWING
C     LIST CONTAINS GROUPS OF PERMISSIBLE IDENTIFICATIONS FOR
C     A SINGLE CALLINNG SEQUENCE.
C
C          1. (Y,QTY,B) (RSD) (XB) (QY)
C
C          2. (Y,QTY,RSD) (B) (XB) (QY)
C
C          3. (Y,QTY,XB) (B) (RSD) (QY)
C
C          4. (Y,QY) (QTY,B) (RSD) (XB)
C
C          5. (Y,QY) (QTY,RSD) (B) (XB)
C
C          6. (Y,QY) (QTY,XB) (B) (RSD)
C
C     IN ANY GROUP THE VALUE RETURNED IN THE ARRAY ALLOCATED TO
C     THE GROUP CORRESPONDS TO THE LAST MEMBER OF THE GROUP.
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     G.W. STEWART, UNIVERSITY OF MARYLAND, ARGONNE NATIONAL LAB.
C
C     DQRSL USES THE FOLLOWING FUNCTIONS AND SUBPROGRAMS.
C
C     BLAS DAXPY,DCOPY,DDOT
C     FORTRAN DABS,MIN0,MOD
C
C     INTERNAL VARIABLES
C
      INTEGER I,J,JJ,JU,KP1
      DOUBLE PRECISION DDOT,T,TEMP
      LOGICAL CB,CQY,CQTY,CR,CXB
C
C
C     SET INFO FLAG.
C
      INFO = 0
C
C     DETERMINE WHAT IS TO BE COMPUTED.
C
      CQY = JOB/10000 .NE. 0
      CQTY = MOD(JOB,10000) .NE. 0
      CB = MOD(JOB,1000)/100 .NE. 0
      CR = MOD(JOB,100)/10 .NE. 0
      CXB = MOD(JOB,10) .NE. 0
      JU = MIN0(K,N-1)
C
C     SPECIAL ACTION WHEN N=1.
C
      IF (JU .NE. 0) GO TO 40
         IF (CQY) QY(1) = Y(1)
         IF (CQTY) QTY(1) = Y(1)
         IF (CXB) XB(1) = Y(1)
         IF (.NOT.CB) GO TO 30
            IF (X(1,1) .NE. 0.0D0) GO TO 10
               INFO = 1
            GO TO 20
   10       CONTINUE
               B(1) = Y(1)/X(1,1)
   20       CONTINUE
   30    CONTINUE
         IF (CR) RSD(1) = 0.0D0
      GO TO 250
   40 CONTINUE
C
C        SET UP TO COMPUTE QY OR QTY.
C
         IF (CQY) CALL DCOPY(N,Y,1,QY,1)
         IF (CQTY) CALL DCOPY(N,Y,1,QTY,1)
         IF (.NOT.CQY) GO TO 70
C
C           COMPUTE QY.
C
            DO 60 JJ = 1, JU
               J = JU - JJ + 1
               IF (QRAUX(J) .EQ. 0.0D0) GO TO 50
                  TEMP = X(J,J)
                  X(J,J) = QRAUX(J)
                  T = -DDOT(N-J+1,X(J,J),1,QY(J),1)/X(J,J)
                  CALL DAXPY(N-J+1,T,X(J,J),1,QY(J),1)
                  X(J,J) = TEMP
   50          CONTINUE
   60       CONTINUE
   70    CONTINUE
         IF (.NOT.CQTY) GO TO 100
C
C           COMPUTE TRANS(Q)*Y.
C
            DO 90 J = 1, JU
               IF (QRAUX(J) .EQ. 0.0D0) GO TO 80
                  TEMP = X(J,J)
                  X(J,J) = QRAUX(J)
                  T = -DDOT(N-J+1,X(J,J),1,QTY(J),1)/X(J,J)
                  CALL DAXPY(N-J+1,T,X(J,J),1,QTY(J),1)
                  X(J,J) = TEMP
   80          CONTINUE
   90       CONTINUE
  100    CONTINUE
C
C        SET UP TO COMPUTE B, RSD, OR XB.
C
         IF (CB) CALL DCOPY(K,QTY,1,B,1)
         KP1 = K + 1
         IF (CXB) CALL DCOPY(K,QTY,1,XB,1)
         IF (CR .AND. K .LT. N) CALL DCOPY(N-K,QTY(KP1),1,RSD(KP1),1)
         IF (.NOT.CXB .OR. KP1 .GT. N) GO TO 120
            DO 110 I = KP1, N
               XB(I) = 0.0D0
  110       CONTINUE
  120    CONTINUE
         IF (.NOT.CR) GO TO 140
            DO 130 I = 1, K
               RSD(I) = 0.0D0
  130       CONTINUE
  140    CONTINUE
         IF (.NOT.CB) GO TO 190
C
C           COMPUTE B.
C
            DO 170 JJ = 1, K
               J = K - JJ + 1
               IF (X(J,J) .NE. 0.0D0) GO TO 150
                  INFO = J
C           ......EXIT
                  GO TO 180
  150          CONTINUE
               B(J) = B(J)/X(J,J)
               IF (J .EQ. 1) GO TO 160
                  T = -B(J)
                  CALL DAXPY(J-1,T,X(1,J),1,B,1)
  160          CONTINUE
  170       CONTINUE
  180       CONTINUE
  190    CONTINUE
         IF (.NOT.CR .AND. .NOT.CXB) GO TO 240
C
C           COMPUTE RSD OR XB AS REQUIRED.
C
            DO 230 JJ = 1, JU
               J = JU - JJ + 1
               IF (QRAUX(J) .EQ. 0.0D0) GO TO 220
                  TEMP = X(J,J)
                  X(J,J) = QRAUX(J)
                  IF (.NOT.CR) GO TO 200
                     T = -DDOT(N-J+1,X(J,J),1,RSD(J),1)/X(J,J)
                     CALL DAXPY(N-J+1,T,X(J,J),1,RSD(J),1)
  200             CONTINUE
                  IF (.NOT.CXB) GO TO 210
                     T = -DDOT(N-J+1,X(J,J),1,XB(J),1)/X(J,J)
                     CALL DAXPY(N-J+1,T,X(J,J),1,XB(J),1)
  210             CONTINUE
                  X(J,J) = TEMP
  220          CONTINUE
  230       CONTINUE
  240    CONTINUE
  250 CONTINUE
      RETURN
      END
      SUBROUTINE DQRDC(X,LDX,N,P,QRAUX,JPVT,WORK,JOB)
      INTEGER LDX,N,P,JOB
      INTEGER JPVT(1)
      DOUBLE PRECISION X(LDX,1),QRAUX(1),WORK(1)
C
C     DQRDC USES HOUSEHOLDER TRANSFORMATIONS TO COMPUTE THE QR
C     FACTORIZATION OF AN N BY P MATRIX X.  COLUMN PIVOTING
C     BASED ON THE 2-NORMS OF THE REDUCED COLUMNS MAY BE
C     PERFORMED AT THE USERS OPTION.
C
C     ON ENTRY
C
C        X       DOUBLE PRECISION(LDX,P), WHERE LDX .GE. N.
C                X CONTAINS THE MATRIX WHOSE DECOMPOSITION IS TO BE
C                COMPUTED.
C
C        LDX     INTEGER.
C                LDX IS THE LEADING DIMENSION OF THE ARRAY X.
C
C        N       INTEGER.
C                N IS THE NUMBER OF ROWS OF THE MATRIX X.
C
C        P       INTEGER.
C                P IS THE NUMBER OF COLUMNS OF THE MATRIX X.
C
C        JPVT    INTEGER(P).
C                JPVT CONTAINS INTEGERS THAT CONTROL THE SELECTION
C                OF THE PIVOT COLUMNS.  THE K-TH COLUMN X(K) OF X
C                IS PLACED IN ONE OF THREE CLASSES ACCORDING TO THE
C                VALUE OF JPVT(K).
C
C                   IF JPVT(K) .GT. 0, THEN X(K) IS AN INITIAL
C                                      COLUMN.
C
C                   IF JPVT(K) .EQ. 0, THEN X(K) IS A FREE COLUMN.
C
C                   IF JPVT(K) .LT. 0, THEN X(K) IS A FINAL COLUMN.
C
C                BEFORE THE DECOMPOSITION IS COMPUTED, INITIAL COLUMNS
C                ARE MOVED TO THE BEGINNING OF THE ARRAY X AND FINAL
C                COLUMNS TO THE END.  BOTH INITIAL AND FINAL COLUMNS
C                ARE FROZEN IN PLACE DURING THE COMPUTATION AND ONLY
C                FREE COLUMNS ARE MOVED.  AT THE K-TH STAGE OF THE
C                REDUCTION, IF X(K) IS OCCUPIED BY A FREE COLUMN
C                IT IS INTERCHANGED WITH THE FREE COLUMN OF LARGEST
C                REDUCED NORM.  JPVT IS NOT REFERENCED IF
C                JOB .EQ. 0.
C
C        WORK    DOUBLE PRECISION(P).
C                WORK IS A WORK ARRAY.  WORK IS NOT REFERENCED IF
C                JOB .EQ. 0.
C
C        JOB     INTEGER.
C                JOB IS AN INTEGER THAT INITIATES COLUMN PIVOTING.
C                IF JOB .EQ. 0, NO PIVOTING IS DONE.
C                IF JOB .NE. 0, PIVOTING IS DONE.
C
C     ON RETURN
C
C        X       X CONTAINS IN ITS UPPER TRIANGLE THE UPPER
C                TRIANGULAR MATRIX R OF THE QR FACTORIZATION.
C                BELOW ITS DIAGONAL X CONTAINS INFORMATION FROM
C                WHICH THE ORTHOGONAL PART OF THE DECOMPOSITION
C                CAN BE RECOVERED.  NOTE THAT IF PIVOTING HAS
C                BEEN REQUESTED, THE DECOMPOSITION IS NOT THAT
C                OF THE ORIGINAL MATRIX X BUT THAT OF X
C                WITH ITS COLUMNS PERMUTED AS DESCRIBED BY JPVT.
C
C        QRAUX   DOUBLE PRECISION(P).
C                QRAUX CONTAINS FURTHER INFORMATION REQUIRED TO RECOVER
C                THE ORTHOGONAL PART OF THE DECOMPOSITION.
C
C        JPVT    JPVT(K) CONTAINS THE INDEX OF THE COLUMN OF THE
C                ORIGINAL MATRIX THAT HAS BEEN INTERCHANGED INTO
C                THE K-TH COLUMN, IF PIVOTING WAS REQUESTED.
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     G.W. STEWART, UNIVERSITY OF MARYLAND, ARGONNE NATIONAL LAB.
C
C     DQRDC USES THE FOLLOWING FUNCTIONS AND SUBPROGRAMS.
C
C     BLAS DAXPY,DDOT,DSCAL,DSWAP,DNRM2
C     FORTRAN DABS,DMAX1,MIN0,DSQRT
C
C     INTERNAL VARIABLES
C
      INTEGER J,JP,L,LP1,LUP,MAXJ,PL,PU
      DOUBLE PRECISION MAXNRM,DNRM2,TT
      DOUBLE PRECISION DDOT,NRMXL,T
      LOGICAL NEGJ,SWAPJ
C
C
      PL = 1
      PU = 0
      IF (JOB .EQ. 0) GO TO 60
C
C        PIVOTING HAS BEEN REQUESTED.  REARRANGE THE COLUMNS
C        ACCORDING TO JPVT.
C
         DO 20 J = 1, P
            SWAPJ = JPVT(J) .GT. 0
            NEGJ = JPVT(J) .LT. 0
            JPVT(J) = J
            IF (NEGJ) JPVT(J) = -J
            IF (.NOT.SWAPJ) GO TO 10
               IF (J .NE. PL) CALL DSWAP(N,X(1,PL),1,X(1,J),1)
               JPVT(J) = JPVT(PL)
               JPVT(PL) = J
               PL = PL + 1
   10       CONTINUE
   20    CONTINUE
         PU = P
         DO 50 JJ = 1, P
            J = P - JJ + 1
            IF (JPVT(J) .GE. 0) GO TO 40
               JPVT(J) = -JPVT(J)
               IF (J .EQ. PU) GO TO 30
                  CALL DSWAP(N,X(1,PU),1,X(1,J),1)
                  JP = JPVT(PU)
                  JPVT(PU) = JPVT(J)
                  JPVT(J) = JP
   30          CONTINUE
               PU = PU - 1
   40       CONTINUE
   50    CONTINUE
   60 CONTINUE
C
C     COMPUTE THE NORMS OF THE FREE COLUMNS.
C
      IF (PU .LT. PL) GO TO 80
      DO 70 J = PL, PU
         QRAUX(J) = DNRM2(N,X(1,J),1)
         WORK(J) = QRAUX(J)
   70 CONTINUE
   80 CONTINUE
C
C     PERFORM THE HOUSEHOLDER REDUCTION OF X.
C
      LUP = MIN0(N,P)
      DO 200 L = 1, LUP
         IF (L .LT. PL .OR. L .GE. PU) GO TO 120
C
C           LOCATE THE COLUMN OF LARGEST NORM AND BRING IT
C           INTO THE PIVOT POSITION.
C
            MAXNRM = 0.0D0
            MAXJ = L
            DO 100 J = L, PU
               IF (QRAUX(J) .LE. MAXNRM) GO TO 90
                  MAXNRM = QRAUX(J)
                  MAXJ = J
   90          CONTINUE
  100       CONTINUE
            IF (MAXJ .EQ. L) GO TO 110
               CALL DSWAP(N,X(1,L),1,X(1,MAXJ),1)
               QRAUX(MAXJ) = QRAUX(L)
               WORK(MAXJ) = WORK(L)
               JP = JPVT(MAXJ)
               JPVT(MAXJ) = JPVT(L)
               JPVT(L) = JP
  110       CONTINUE
  120    CONTINUE
         QRAUX(L) = 0.0D0
         IF (L .EQ. N) GO TO 190
C
C           COMPUTE THE HOUSEHOLDER TRANSFORMATION FOR COLUMN L.
C
            NRMXL = DNRM2(N-L+1,X(L,L),1)
            IF (NRMXL .EQ. 0.0D0) GO TO 180
C              DSIGN requires libF77, so rewrite; jdr July 96
C              IF (X(L,L) .NE. 0.0D0) NRMXL = DSIGN(NRMXL,X(L,L))
               IF (X(L,L) .LT. 0.0D0) NRMXL = -DABS(NRMXL)
               IF (X(L,L) .GT. 0.0D0) NRMXL = DABS(NRMXL)
               CALL DSCAL(N-L+1,1.0D0/NRMXL,X(L,L),1)
               X(L,L) = 1.0D0 + X(L,L)
C
C              APPLY THE TRANSFORMATION TO THE REMAINING COLUMNS,
C              UPDATING THE NORMS.
C
               LP1 = L + 1
               IF (P .LT. LP1) GO TO 170
               DO 160 J = LP1, P
                  T = -DDOT(N-L+1,X(L,L),1,X(L,J),1)/X(L,L)
                  CALL DAXPY(N-L+1,T,X(L,L),1,X(L,J),1)
                  IF (J .LT. PL .OR. J .GT. PU) GO TO 150
                  IF (QRAUX(J) .EQ. 0.0D0) GO TO 150
                     TT = 1.0D0 - (DABS(X(L,J))/QRAUX(J))**2
                     TT = DMAX1(TT,0.0D0)
                     T = TT
                     TT = 1.0D0 + 0.05D0*TT*(QRAUX(J)/WORK(J))**2
                     IF (TT .EQ. 1.0D0) GO TO 130
                        QRAUX(J) = QRAUX(J)*DSQRT(T)
                     GO TO 140
  130                CONTINUE
                        QRAUX(J) = DNRM2(N-L,X(L+1,J),1)
                        WORK(J) = QRAUX(J)
  140                CONTINUE
  150             CONTINUE
  160          CONTINUE
  170          CONTINUE
C
C              SAVE THE TRANSFORMATION.
C
               QRAUX(L) = X(L,L)
               X(L,L) = -NRMXL
  180       CONTINUE
  190    CONTINUE
  200 CONTINUE
      RETURN
      END
      DOUBLE PRECISION FUNCTION DDOT(N,DX,INCX,DY,INCY)
C
C     FORMS THE DOT PRODUCT OF TWO VECTORS.
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE PRECISION DX(1),DY(1),DTEMP
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N
C
      DDOT = 0.0D0
      DTEMP = 0.0D0
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
C
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C          NOT EQUAL TO 1
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DTEMP = DTEMP + DX(IX)*DY(IY)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      DDOT = DTEMP
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DTEMP = DTEMP + DX(I)*DY(I)
   30 CONTINUE
      IF( N .LT. 5 ) GO TO 60
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        DTEMP = DTEMP + DX(I)*DY(I) + DX(I + 1)*DY(I + 1) +
     *   DX(I + 2)*DY(I + 2) + DX(I + 3)*DY(I + 3) + DX(I + 4)*DY(I + 4)
   50 CONTINUE
   60 DDOT = DTEMP
      RETURN
      END
      SUBROUTINE  DCOPY(N,DX,INCX,DY,INCY)
C
C     COPIES A VECTOR, X, TO A VECTOR, Y.
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE PRECISION DX(1),DY(1)
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N
C
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
C
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C          NOT EQUAL TO 1
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 M = MOD(N,7)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DY(I) = DX(I)
   30 CONTINUE
      IF( N .LT. 7 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,7
        DY(I) = DX(I)
        DY(I + 1) = DX(I + 1)
        DY(I + 2) = DX(I + 2)
        DY(I + 3) = DX(I + 3)
        DY(I + 4) = DX(I + 4)
        DY(I + 5) = DX(I + 5)
        DY(I + 6) = DX(I + 6)
   50 CONTINUE
      RETURN
      END
      SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY)
C
C     CONSTANT TIMES A VECTOR PLUS A VECTOR.
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE PRECISION DX(1),DY(1),DA
      INTEGER I,INCX,INCY,M,MP1,N
C
      IF(N.LE.0)RETURN
      IF (DA .EQ. 0.0D0) RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
C
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C          NOT EQUAL TO 1
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = DY(IY) + DA*DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 M = MOD(N,4)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DY(I) = DY(I) + DA*DX(I)
   30 CONTINUE
      IF( N .LT. 4 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
        DY(I) = DY(I) + DA*DX(I)
        DY(I + 1) = DY(I + 1) + DA*DX(I + 1)
        DY(I + 2) = DY(I + 2) + DA*DX(I + 2)
        DY(I + 3) = DY(I + 3) + DA*DX(I + 3)
   50 CONTINUE
      RETURN
      END
      SUBROUTINE  DSWAP (N,DX,INCX,DY,INCY)
C
C     INTERCHANGES TWO VECTORS.
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL ONE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE PRECISION DX(1),DY(1),DTEMP
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N
C
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
C
C       CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS NOT EQUAL
C         TO 1
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DTEMP = DX(IX)
        DX(IX) = DY(IY)
        DY(IY) = DTEMP
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C       CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C       CLEAN-UP LOOP
C
   20 M = MOD(N,3)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DTEMP = DX(I)
        DX(I) = DY(I)
        DY(I) = DTEMP
   30 CONTINUE
      IF( N .LT. 3 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,3
        DTEMP = DX(I)
        DX(I) = DY(I)
        DY(I) = DTEMP
        DTEMP = DX(I + 1)
        DX(I + 1) = DY(I + 1)
        DY(I + 1) = DTEMP
        DTEMP = DX(I + 2)
        DX(I + 2) = DY(I + 2)
        DY(I + 2) = DTEMP
   50 CONTINUE
      RETURN
      END
      SUBROUTINE  DSCAL(N,DA,DX,INCX)
C
C     SCALES A VECTOR BY A CONSTANT.
C     USES UNROLLED LOOPS FOR INCREMENT EQUAL TO ONE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE PRECISION DA,DX(1)
      INTEGER I,INCX,M,MP1,N,NINCX
C
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1)GO TO 20
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
      NINCX = N*INCX
      DO 10 I = 1,NINCX,INCX
        DX(I) = DA*DX(I)
   10 CONTINUE
      RETURN
C
C        CODE FOR INCREMENT EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DX(I) = DA*DX(I)
   30 CONTINUE
      IF( N .LT. 5 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        DX(I) = DA*DX(I)
        DX(I + 1) = DA*DX(I + 1)
        DX(I + 2) = DA*DX(I + 2)
        DX(I + 3) = DA*DX(I + 3)
        DX(I + 4) = DA*DX(I + 4)
   50 CONTINUE
      RETURN
      END
      DOUBLE PRECISION FUNCTION DNRM2 ( N, DX, INCX)
      INTEGER          NEXT
      DOUBLE PRECISION   DX(1), CUTLO, CUTHI, HITEST, SUM, XMAX,ZERO,ONE
      DATA   ZERO, ONE /0.0D0, 1.0D0/
C
C     EUCLIDEAN NORM OF THE N-VECTOR STORED IN DX() WITH STORAGE
C     INCREMENT INCX .
C     IF    N .LE. 0 RETURN WITH RESULT = 0.
C     IF N .GE. 1 THEN INCX MUST BE .GE. 1
C
C           C.L.LAWSON, 1978 JAN 08
C
C     FOUR PHASE METHOD     USING TWO BUILT-IN CONSTANTS THAT ARE
C     HOPEFULLY APPLICABLE TO ALL MACHINES.
C         CUTLO = MAXIMUM OF  DSQRT(U/EPS)  OVER ALL KNOWN MACHINES.
C         CUTHI = MINIMUM OF  DSQRT(V)      OVER ALL KNOWN MACHINES.
C     WHERE
C         EPS = SMALLEST NO. SUCH THAT EPS + 1. .GT. 1.
C         U   = SMALLEST POSITIVE NO.   (UNDERFLOW LIMIT)
C         V   = LARGEST  NO.            (OVERFLOW  LIMIT)
C
C     BRIEF OUTLINE OF ALGORITHM..
C
C     PHASE 1    SCANS ZERO COMPONENTS.
C     MOVE TO PHASE 2 WHEN A COMPONENT IS NONZERO AND .LE. CUTLO
C     MOVE TO PHASE 3 WHEN A COMPONENT IS .GT. CUTLO
C     MOVE TO PHASE 4 WHEN A COMPONENT IS .GE. CUTHI/M
C     WHERE M = N FOR X() REAL AND M = 2*N FOR COMPLEX.
C
C     VALUES FOR CUTLO AND CUTHI..
C     FROM THE ENVIRONMENTAL PARAMETERS LISTED IN THE IMSL CONVERTER
C     DOCUMENT THE LIMITING VALUES ARE AS FOLLOWS..
C     CUTLO, S.P.   U/EPS = 2**(-102) FOR  HONEYWELL.  CLOSE SECONDS ARE
C                   UNIVAC AND DEC AT 2**(-103)
C                   THUS CUTLO = 2**(-51) = 4.44089E-16
C     CUTHI, S.P.   V = 2**127 FOR UNIVAC, HONEYWELL, AND DEC.
C                   THUS CUTHI = 2**(63.5) = 1.30438E19
C     CUTLO, D.P.   U/EPS = 2**(-67) FOR HONEYWELL AND DEC.
C                   THUS CUTLO = 2**(-33.5) = 8.23181D-11
C     CUTHI, D.P.   SAME AS S.P.  CUTHI = 1.30438D19
C     DATA CUTLO, CUTHI / 8.232D-11,  1.304D19 /
C     DATA CUTLO, CUTHI / 4.441E-16,  1.304E19 /
      DATA CUTLO, CUTHI / 8.232D-11,  1.304D19 /
C
      IF(N .GT. 0) GO TO 10
         DNRM2  = ZERO
         GO TO 300
C
   10 ASSIGN 30 TO NEXT
      SUM = ZERO
      NN = N * INCX
C                                                 BEGIN MAIN LOOP
      I = 1
   20    GO TO NEXT,(30, 50, 70, 110)
   30 IF( DABS(DX(I)) .GT. CUTLO) GO TO 85
      ASSIGN 50 TO NEXT
      XMAX = ZERO
C
C                        PHASE 1.  SUM IS ZERO
C
   50 IF( DX(I) .EQ. ZERO) GO TO 200
      IF( DABS(DX(I)) .GT. CUTLO) GO TO 85
C
C                                PREPARE FOR PHASE 2.
      ASSIGN 70 TO NEXT
      GO TO 105
C
C                                PREPARE FOR PHASE 4.
C
  100 I = J
      ASSIGN 110 TO NEXT
      SUM = (SUM / DX(I)) / DX(I)
  105 XMAX = DABS(DX(I))
      GO TO 115
C
C                   PHASE 2.  SUM IS SMALL.
C                             SCALE TO AVOID DESTRUCTIVE UNDERFLOW.
C
   70 IF( DABS(DX(I)) .GT. CUTLO ) GO TO 75
C
C                     COMMON CODE FOR PHASES 2 AND 4.
C                     IN PHASE 4 SUM IS LARGE.  SCALE TO AVOID OVERFLOW.
C
  110 IF( DABS(DX(I)) .LE. XMAX ) GO TO 115
         SUM = ONE + SUM * (XMAX / DX(I))**2
         XMAX = DABS(DX(I))
         GO TO 200
C
  115 SUM = SUM + (DX(I)/XMAX)**2
      GO TO 200
C
C
C                  PREPARE FOR PHASE 3.
C
   75 SUM = (SUM * XMAX) * XMAX
C
C
C     FOR REAL OR D.P. SET HITEST = CUTHI/N
C     FOR COMPLEX      SET HITEST = CUTHI/(2*N)
C
   85 HITEST = CUTHI/FLOAT( N )
C
C                   PHASE 3.  SUM IS MID-RANGE.  NO SCALING.
C
      DO 95 J =I,NN,INCX
      IF(DABS(DX(J)) .GE. HITEST) GO TO 100
   95    SUM = SUM + DX(J)**2
      DNRM2 = DSQRT( SUM )
      GO TO 300
C
  200 CONTINUE
      I = I + INCX
      IF ( I .LE. NN ) GO TO 20
C
C              END OF MAIN LOOP.
C
C              COMPUTE SQUARE ROOT AND ADJUST FOR SCALING.
C
      DNRM2 = XMAX * DSQRT(SUM)
  300 CONTINUE
      RETURN
      END
      subroutine supdbl (n,x,y,w,iper,span,alpha,smo,sc)
c
c------------------------------------------------------------------
c
c super smoother (friedman, 1984).
c
c version 10/10/84
c
c coded  and copywrite <c> 1984 by:
c
c                        jerome h. friedman
c                     department of statistics
c                               and
c                stanford linear accelerator center
c                        stanford university
c
c all rights reserved.
c
c Modified to double precision by: James Reimann, November 1993
c
c input:
c    n : number of observations (x,y - pairs).
c    x(n) : ordered abscissa values.
c    y(n) : corresponding ordinate (response) values.
c    w(n) : weight for each (x,y) observation.
c    iper : periodic variable flag.
c       iper=1 => x is ordered interval variable.
c       iper=2 => x is a periodic variable with values
c                 in the range (0.0,1.0) and peroid 1.0.
c    span : smoother span (fraction of observations in window).
c           span=0.0 => automatic (variable) span selection.
c    alpha : controles high frequency (small span) penality
c            used with automatic span selection (bass tone control).
c            (alpha.le.0.0 or alpha.gt.10.0 => no effect.)
c output:
c   smo(n) : smoothed ordinate (response) values.
c scratch:
c   sc(n,7) : internal working storage.
c
c note:
c    for small samples (n < 40) or if there are substantial serial
c    correlations between obserations close in x - value, then
c    a prespecified fixed span smoother (span > 0) should be
c    used. reasonable span values are 0.3 to 0.5.
c
c------------------------------------------------------------------
c
      implicit double precision (a-h,o-z)
      dimension x(n),y(n),w(n),smo(n),sc(n,7)
      common /spans/ spans(3) /consts/ big,sml,eps
      if (x(n).gt.x(1)) go to 30
      sy=0.0d0
      sw=sy
      do 10 j=1,n
      sy=sy+w(j)*y(j)
      sw=sw+w(j)
 10   continue
      a=0.0d0
      if (sw.gt.0.0d0) a=sy/sw
      do 20 j=1,n
      smo(j)=a
 20   continue
      return
 30   i=n/4
      j=3*i
      scale=x(j)-x(i)
 40   if (scale.gt.0.0d0) go to 50
      if (j.lt.n) j=j+1
      if (i.gt.1) i=i-1
      scale=x(j)-x(i)
      go to 40
 50   vsmlsq=(eps*scale)**2
      jper=iper
      if (iper.eq.2.and.(x(1).lt.0.0d0.or.x(n).gt.1d0)) jper=1
      if (jper.lt.1.or.jper.gt.2) jper=1
      if (span.le.0.0d0) go to 60
      call smooth (n,x,y,w,span,jper,vsmlsq,smo,sc)
      return
 60   do 70 i=1,3
      call smooth (n,x,y,w,spans(i),jper,vsmlsq,sc(1,2*i-1),sc(1,7))
      call smooth (n,x,sc(1,7),w,spans(2),-jper,vsmlsq,sc(1,2*i),h)
 70   continue
      do 90 j=1,n
      resmin=big
      do 80 i=1,3
      if (sc(j,2*i).ge.resmin) go to 80
      resmin=sc(j,2*i)
      sc(j,7)=spans(i)
 80   continue
      if (alpha.gt.0.0d0.and.alpha.le.10.0d0.and.resmin.lt.sc(j,6).
     1and.resmin.gt.0.0d0) sc(j,7)=sc(j,7)+(spans(3)-sc(j,7))*max
     2(sml,resmin/sc(j,6))**(10.0d0-alpha)
 90   continue
      call smooth (n,x,sc(1,7),w,spans(2),-jper,vsmlsq,sc(1,2),h)
      do 110 j=1,n
      if (sc(j,2).le.spans(1)) sc(j,2)=spans(1)
      if (sc(j,2).ge.spans(3)) sc(j,2)=spans(3)
      f=sc(j,2)-spans(2)
      if (f.ge.0.0d0) go to 100
      f=-f/(spans(2)-spans(1))
      sc(j,4)=(1d0-f)*sc(j,3)+f*sc(j,1)
      go to 110
 100  f=f/(spans(3)-spans(2))
      sc(j,4)=(1d0-f)*sc(j,3)+f*sc(j,5)
 110  continue
      call smooth (n,x,sc(1,4),w,spans(1),-jper,vsmlsq,smo,h)
      return
      end
      subroutine smooth (n,x,y,w,span,iper,vsmlsq,smo,acvr)
      implicit double precision (a-h,o-z)
      dimension x(n),y(n),w(n),smo(n),acvr(n)
      integer in,out
      double precision wt,fbo,fbw,xm,ym,tmp,var,cvar,a,h,sy
      xm=0.0d0
      ym=xm
      var=ym
      cvar=var
      fbw=cvar
      jper=iabs(iper)
      ibw=0.5*span*n+0.5
      if (ibw.lt.2) ibw=2
      it=2*ibw+1
      do 20 i=1,it
      j=i
      if (jper.eq.2) j=i-ibw-1
c 666  print *, i, j, ibw
c 667  print *, x(j)
      xti=x(j)
      if (j.ge.1) go to 10
      j=n+j
      xti=x(j)-1d0
 10   wt=w(j)
      fbo=fbw
      fbw=fbw+wt
      if (fbw.gt.0.0d0) xm=(fbo*xm+wt*xti)/fbw
      if (fbw.gt.0.0d0) ym=(fbo*ym+wt*y(j))/fbw
      tmp=0.0d0
      if (fbo.gt.0.0d0) tmp=fbw*wt*(xti-xm)/fbo
      var=var+tmp*(xti-xm)
      cvar=cvar+tmp*(y(j)-ym)
 20   continue
      do 80 j=1,n
      out=j-ibw-1
      in=j+ibw
      if ((jper.ne.2).and.(out.lt.1.or.in.gt.n)) go to 60
      if (out.ge.1) go to 30
      out=n+out
      xto=x(out)-1d0
      xti=x(in)
      go to 50
 30   if (in.le.n) go to 40
      in=in-n
      xti=x(in)+1d0
      xto=x(out)
      go to 50
 40   xto=x(out)
      xti=x(in)
 50   wt=w(out)
      fbo=fbw
      fbw=fbw-wt
      tmp=0.0d0
      if (fbw.gt.0.0d0) tmp=fbo*wt*(xto-xm)/fbw
      var=var-tmp*(xto-xm)
      cvar=cvar-tmp*(y(out)-ym)
      if (fbw.gt.0.0d0) xm=(fbo*xm-wt*xto)/fbw
      if (fbw.gt.0.0d0) ym=(fbo*ym-wt*y(out))/fbw
      wt=w(in)
      fbo=fbw
      fbw=fbw+wt
      if (fbw.gt.0.0d0) xm=(fbo*xm+wt*xti)/fbw
      if (fbw.gt.0.0d0) ym=(fbo*ym+wt*y(in))/fbw
      tmp=0.0d0
      if (fbo.gt.0.0d0) tmp=fbw*wt*(xti-xm)/fbo
      var=var+tmp*(xti-xm)
      cvar=cvar+tmp*(y(in)-ym)
 60   a=0.0d0
      if (var.gt.vsmlsq) a=cvar/var
      smo(j)=a*(x(j)-xm)+ym
      if (iper.le.0) go to 80
      h=0.0d0
      if (fbw.gt.0.0d0) h=1d0/fbw
      if (var.gt.vsmlsq) h=h+(x(j)-xm)**2/var
      acvr(j)=0.0d0
      a=1d0-w(j)*h
      if (a.le.0.0d0) go to 70
      acvr(j)=abs(y(j)-smo(j))/a
      go to 80
 70   if (j.le.1) go to 80
      acvr(j)=acvr(j-1)
 80   continue
      j=1
 90   j0=j
      sy=smo(j)*w(j)
      fbw=w(j)
      if (j.ge.n) go to 110
 100  if (x(j+1).gt.x(j)) go to 110
      j=j+1
      sy=sy+w(j)*smo(j)
      fbw=fbw+w(j)
      if (j.ge.n) go to 110
      go to 100
 110  if (j.le.j0) go to 130
      a=0.0d0
      if (fbw.gt.0.0d0) a=sy/fbw
      do 120 i=j0,j
      smo(i)=a
 120  continue
 130  j=j+1
      if (j.gt.n) go to 140
      go to 90
 140  return
      end
      block data
      implicit double precision (a-h,o-z)
      common /spans/ spans(3) /consts/ big,sml,eps
c
c---------------------------------------------------------------
c
c this sets the compile time (default) values for various
c internal parameters :
c
c spans : span values for the three running linear smoothers.
c spans(1) : tweeter span.
c spans(2) : midrange span.
c spans(3) : woofer span.
c (these span values should be changed only with care.)
c big : a large representable floating point number.
c sml : a small number. should be set so that (sml)**(10.0) does
c       not cause floating point underflow.
c eps : used to numerically stabilize slope calculations for
c       running linear fits.
c
c these parameter values can be changed by declaring the
c relevant labeled common in the main program and resetting
c them with executable statements.
c
c-----------------------------------------------------------------
c
      data spans,big,sml,eps /5.d-2,2.d-1,5.d-1,1.d20,1.d-7,1.d-3/
      end
      subroutine percur(iopt,m,x,y,w,k,s,nest,n,t,c,fp,
     * wrk,lwrk,iwrk,ier)
c  given the set of data points (x(i),y(i)) and the set of positive
c  numbers w(i),i=1,2,...,m-1, subroutine percur determines a smooth
c  periodic spline approximation of degree k with period per=x(m)-x(1).
c  if iopt=-1 percur calculates the weighted least-squares periodic
c  spline according to a given set of knots.
c  if iopt>=0 the number of knots of the spline s(x) and the position
c  t(j),j=1,2,...,n is chosen automatically by the routine. the smooth-
c  ness of s(x) is then achieved by minimalizing the discontinuity
c  jumps of the k-th derivative of s(x) at the knots t(j),j=k+2,k+3,...,
c  n-k-1. the amount of smoothness is determined by the condition that
c  f(p)=sum((w(i)*(y(i)-s(x(i))))**2) be <= s, with s a given non-
c  negative constant, called the smoothing factor.
c  the fit s(x) is given in the b-spline representation (b-spline coef-
c  ficients c(j),j=1,2,...,n-k-1) and can be evaluated by means of
c  subroutine splev.
c
c  calling sequence:
c     call percur(iopt,m,x,y,w,k,s,nest,n,t,c,fp,wrk,
c    * lwrk,iwrk,ier)
c
c  parameters:
c   iopt  : integer flag. on entry iopt must specify whether a weighted
c           least-squares spline (iopt=-1) or a smoothing spline (iopt=
c           0 or 1) must be determined. if iopt=0 the routine will start
c           with an initial set of knots t(i)=x(1)+(x(m)-x(1))*(i-k-1),
c           i=1,2,...,2*k+2. if iopt=1 the routine will continue with
c           the knots found at the last call of the routine.
c           attention: a call with iopt=1 must always be immediately
c           preceded by another call with iopt=1 or iopt=0.
c           unchanged on exit.
c   m     : integer. on entry m must specify the number of data points.
c           m > 1. unchanged on exit.
c   x     : double array of dimension at least (m). before entry, x(i)
c           must be set to the i-th value of the independent variable x,
c           for i=1,2,...,m. these values must be supplied in strictly
c           ascending order. x(m) only indicates the length of the
c           period of the spline, i.e per=x(m)-x(1).
c           unchanged on exit.
c   y     : double array of dimension at least (m). before entry, y(i)
c           must be set to the i-th value of the dependent variable y,
c           for i=1,2,...,m-1. the element y(m) is not used.
c           unchanged on exit.
c   w     : double array of dimension at least (m). before entry, w(i)
c           must be set to the i-th value in the set of weights. the
c           w(i) must be strictly positive. w(m) is not used.
c           see also further comments. unchanged on exit.
c   k     : integer. on entry k must specify the degree of the spline.
c           1<=k<=5. it is recommended to use cubic splines (k=3).
c           the user is strongly dissuaded from choosing k even,together
c           with a small s-value. unchanged on exit.
c   s     : double .on entry (in case iopt>=0) s must specify the smoothing
c           factor. s >=0. unchanged on exit.
c           for advice on the choice of s see further comments.
c   nest  : integer. on entry nest must contain an over-estimate of the
c           total number of knots of the spline returned, to indicate
c           the storage space available to the routine. nest >=2*k+2.
c           in most practical situation nest=m/2 will be sufficient.
c           always large enough is nest=m+2*k,the number of knots needed
c           for interpolation (s=0). unchanged on exit.
c   n     : integer.
c           unless ier = 10 (in case iopt >=0), n will contain the
c           total number of knots of the spline approximation returned.
c           if the computation mode iopt=1 is used this value of n
c           should be left unchanged between subsequent calls.
c           in case iopt=-1, the value of n must be specified on entry.
c   t     : double array of dimension at least (nest).
c           on succesful exit, this array will contain the knots of the
c           spline,i.e. the position of the interior knots t(k+2),t(k+3)
c           ...,t(n-k-1) as well as the position of the additional knots
c           t(1),t(2),...,t(k+1)=x(1) and t(n-k)=x(m),..,t(n) needed for
c           the b-spline representation.
c           if the computation mode iopt=1 is used, the values of t(1),
c           t(2),...,t(n) should be left unchanged between subsequent
c           calls. if the computation mode iopt=-1 is used, the values
c           t(k+2),...,t(n-k-1) must be supplied by the user, before
c           entry. see also the restrictions (ier=10).
c   c     : double array of dimension at least (nest).
c           on succesful exit, this array will contain the coefficients
c           c(1),c(2),..,c(n-k-1) in the b-spline representation of s(x)
c   fp    : double . unless ier = 10, fp contains the weighted sum of
c           squared residuals of the spline approximation returned.
c   wrk   : double array of dimension at least (m*(k+1)+nest*(8+5*k)).
c           used as working space. if the computation mode iopt=1 is
c           used, the values wrk(1),...,wrk(n) should be left unchanged
c           between subsequent calls.
c   lwrk  : integer. on entry,lwrk must specify the actual dimension of
c           the array wrk as declared in the calling (sub)program. lwrk
c           must not be too small (see wrk). unchanged on exit.
c   iwrk  : integer array of dimension at least (nest).
c           used as working space. if the computation mode iopt=1 is
c           used,the values iwrk(1),...,iwrk(n) should be left unchanged
c           between subsequent calls.
c   ier   : integer. unless the routine detects an error, ier contains a
c           non-positive value on exit, i.e.
c    ier=0  : normal return. the spline returned has a residual sum of
c             squares fp such that abs(fp-s)/s <= tol with tol a relat-
c             ive tolerance set to 0.001 by the program.
c    ier=-1 : normal return. the spline returned is an interpolating
c             periodic spline (fp=0).
c    ier=-2 : normal return. the spline returned is the weighted least-
c             squares constant. in this extreme case fp gives the upper
c             bound fp0 for the smoothing factor s.
c    ier=1  : error. the required storage space exceeds the available
c             storage space, as specified by the parameter nest.
c             probably causes : nest too small. if nest is already
c             large (say nest > m/2), it may also indicate that s is
c             too small
c             the approximation returned is the least-squares periodic
c             spline according to the knots t(1),t(2),...,t(n). (n=nest)
c             the parameter fp gives the corresponding weighted sum of
c             squared residuals (fp>s).
c    ier=2  : error. a theoretically impossible result was found during
c             the iteration proces for finding a smoothing spline with
c             fp = s. probably causes : s too small.
c             there is an approximation returned but the corresponding
c             weighted sum of squared residuals does not satisfy the
c             condition abs(fp-s)/s < tol.
c    ier=3  : error. the maximal number of iterations maxit (set to 20
c             by the program) allowed for finding a smoothing spline
c             with fp=s has been reached. probably causes : s too small
c             there is an approximation returned but the corresponding
c             weighted sum of squared residuals does not satisfy the
c             condition abs(fp-s)/s < tol.
c    ier=10 : error. on entry, the input data are controlled on validity
c             the following restrictions must be satisfied.
c             -1<=iopt<=1, 1<=k<=5, m>1, nest>2*k+2, w(i)>0,i=1,...,m-1
c             x(1)<x(2)<...<x(m), lwrk>=(k+1)*m+nest*(8+5*k)
c             if iopt=-1: 2*k+2<=n<=min(nest,m+2*k)
c                         x(1)<t(k+2)<t(k+3)<...<t(n-k-1)<x(m)
c                       the schoenberg-whitney conditions, i.e. there
c                       must be a subset of data points xx(j) with
c                       xx(j) = x(i) or x(i)+(x(m)-x(1)) such that
c                         t(j) < xx(j) < t(j+k+1), j=k+1,...,n-k-1
c             if iopt>=0: s>=0
c                         if s=0 : nest >= m+2*k
c             if one of these conditions is found to be violated,control
c             is immediately repassed to the calling program. in that
c             case there is no approximation returned.
c
c  further comments:
c   by means of the parameter s, the user can control the tradeoff
c   between closeness of fit and smoothness of fit of the approximation.
c   if s is too large, the spline will be too smooth and signal will be
c   lost ; if s is too small the spline will pick up too much noise. in
c   the extreme cases the program will return an interpolating periodic
c   spline if s=0 and the weighted least-squares constant if s is very
c   large. between these extremes, a properly chosen s will result in
c   a good compromise between closeness of fit and smoothness of fit.
c   to decide whether an approximation, corresponding to a certain s is
c   satisfactory the user is highly recommended to inspect the fits
c   graphically.
c   recommended values for s depend on the weights w(i). if these are
c   taken as 1/d(i) with d(i) an estimate of the standard deviation of
c   y(i), a good s-value should be found in the range (m-sqrt(2*m),m+
c   sqrt(2*m)). if nothing is known about the statistical error in y(i)
c   each w(i) can be set equal to one and s determined by trial and
c   error, taking account of the comments above. the best is then to
c   start with a very large value of s ( to determine the least-squares
c   constant and the corresponding upper bound fp0 for s) and then to
c   progressively decrease the value of s ( say by a factor 10 in the
c   beginning, i.e. s=fp0/10, fp0/100,...and more carefully as the
c   approximation shows more detail) to obtain closer fits.
c   to economize the search for a good s-value the program provides with
c   different modes of computation. at the first call of the routine, or
c   whenever he wants to restart with the initial set of knots the user
c   must set iopt=0.
c   if iopt=1 the program will continue with the set of knots found at
c   the last call of the routine. this will save a lot of computation
c   time if percur is called repeatedly for different values of s.
c   the number of knots of the spline returned and their location will
c   depend on the value of s and on the complexity of the shape of the
c   function underlying the data. but, if the computation mode iopt=1
c   is used, the knots returned may also depend on the s-values at
c   previous calls (if these were smaller). therefore, if after a number
c   of trials with different s-values and iopt=1, the user can finally
c   accept a fit as satisfactory, it may be worthwhile for him to call
c   percur once more with the selected value for s but now with iopt=0.
c   indeed, percur may then return an approximation of the same quality
c   of fit but with fewer knots and therefore better if data reduction
c   is also an important objective for the user.
c
c  other subroutines required:
c    fpbacp,fpbspl,fpchep,fpperi,fpdisc,fpgivs,fpknot,fprati,fprota
c
c  references:
c   dierckx p. : algorithms for smoothing data with periodic and
c                parametric splines, computer graphics and image
c                processing 20 (1982) 171-184.
c   dierckx p. : algorithms for smoothing data with periodic and param-
c                etric splines, report tw55, dept. computer science,
c                k.u.leuven, 1981.
c   dierckx p. : curve and surface fitting with splines, monographs on
c                numerical analysis, oxford university press, 1993.
c
c  author:
c    p.dierckx
c    dept. computer science, k.u. leuven
c    celestijnenlaan 200a, b-3001 heverlee, belgium.
c    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
c
c  creation date : may 1979
c  latest update : march 1987
c
c  ..
c  ..scalar arguments..
      double precision s,fp
      integer iopt,m,k,nest,n,lwrk,ier
c  ..array arguments..
      double precision x(m),y(m),w(m),t(nest),c(nest),wrk(lwrk)
      integer iwrk(nest)
c  ..local scalars..
      double precision per,tol
      integer i,ia1,ia2,ib,ifp,ig1,ig2,iq,iz,i1,i2,j1,j2,k1,k2,lwest,
     * maxit,m1,nmin
c  ..subroutine references..
c    perper,pcheck
c  ..
c  we set up the parameters tol and maxit
      maxit = 20
      tol = 0.1d-02
c  before starting computations a data check is made. if the input data
c  are invalid, control is immediately repassed to the calling program.
      ier = 10
      if(k.le.0 .or. k.gt.5) go to 50
      k1 = k+1
      k2 = k1+1
      ier = 11
      if(iopt.lt.(-1) .or. iopt.gt.1) go to 50
      nmin = 2*k1
      ier = 12
      if(m.lt.2 .or. nest.lt.nmin) go to 50
      lwest = m*k1+nest*(8+5*k)
      ier = 13
      if(lwrk.lt.lwest) go to 50
      m1 = m-1
      ier = 14
      do 10 i=1,m1
         if(x(i).gt.x(i+1) .or. w(i).le.0.0d0) go to 50
C        if(x(i).ge.x(i+1) .or. w(i).le.0.0d0) go to 50 (original version)
  10  continue
      ier = 15
      if(iopt.ge.0) go to 30
      if(n.le.nmin .or. n.gt.nest) go to 50
      per = x(m)-x(1)
      j1 = k1
      t(j1) = x(1)
      i1 = n-k
      t(i1) = x(m)
      j2 = j1
      i2 = i1
      do 20 i=1,k
         i1 = i1+1
         i2 = i2-1
         j1 = j1+1
         j2 = j2-1
         t(j2) = t(i2)-per
         t(i1) = t(j1)+per
  20  continue
      call fpchep(x,m,t,n,k,ier)
      if(ier) 50,40,50
  30  ier = 17
      if(s.lt.0.0d0) go to 50
      ier = 18
      if(s.eq.0.0d0 .and. nest.lt.(m+2*k)) go to 50
      ier = 0
c we partition the working space and determine the spline approximation.
  40  ifp = 1
      iz = ifp+nest
      ia1 = iz+nest
      ia2 = ia1+nest*k1
      ib = ia2+nest*k
      ig1 = ib+nest*k2
      ig2 = ig1+nest*k2
      iq = ig2+nest*k1
      call fpperi(iopt,x,y,w,m,k,s,nest,tol,maxit,k1,k2,n,t,c,fp,
     * wrk(ifp),wrk(iz),wrk(ia1),wrk(ia2),wrk(ib),wrk(ig1),wrk(ig2),
     * wrk(iq),iwrk,ier)
  50  return
      end
      subroutine fpperi(iopt,x,y,w,m,k,s,nest,tol,maxit,k1,k2,n,t,c,
     * fp,fpint,z,a1,a2,b,g1,g2,q,nrdata,ier)
c  ..
c  ..scalar arguments..
      double precision s,tol,fp
      integer iopt,m,k,nest,maxit,k1,k2,n,ier
c  ..array arguments..
      double precision x(m),y(m),w(m),t(nest),c(nest),fpint(nest),
     * z(nest),
     * a1(nest,k1),a2(nest,k),b(nest,k2),g1(nest,k2),g2(nest,k1),
     * q(m,k1)
      integer nrdata(nest)
c  ..local scalars..
      double precision acc,cos,c1,d1,fpart,fpms,fpold,fp0,f1,f2,f3,p,
     * per,pinv,
     * piv,p1,p2,p3,sin,store,term,wi,xi,yi,rn,one,con1,con4,con9,half
      integer i,ich1,ich3,ij,ik,it,iter,i1,i2,i3,j,jk,jper,j1,j2,kk,
     * kk1,k3,l,l0,l1,l5,mm,m1,new,nk1,nk2,nmax,nmin,nplus,npl1,
     * nrint,n10,n11,n7,n8
c  ..local arrays..
      double precision h(6),h1(7),h2(6)
c  ..function references..
      double precision abs,fprati
      integer max0,min0
c  ..subroutine references..
c    fpbacp,fpbspl,fpgivs,fpdisc,fpknot,fprota
c  ..
c  set constants
      one = 0.1d+01
      con1 = 0.1d0
      con9 = 0.9d0
      con4 = 0.4d-01
      half = 0.5d0
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  part 1: determination of the number of knots and their position     c
c  **************************************************************      c
c  given a set of knots we compute the least-squares periodic spline   c
c  sinf(x). if the sum f(p=inf) <= s we accept the choice of knots.    c
c  the initial choice of knots depends on the value of s and iopt.     c
c    if s=0 we have spline interpolation; in that case the number of   c
c    knots equals nmax = m+2*k.                                        c
c    if s > 0 and                                                      c
c      iopt=0 we first compute the least-squares polynomial of         c
c      degree k; n = nmin = 2*k+2. since s(x) must be periodic we      c
c      find that s(x) is a constant function.                          c
c      iopt=1 we start with the set of knots found at the last         c
c      call of the routine, except for the case that s > fp0; then     c
c      we compute directly the least-squares periodic polynomial.      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      m1 = m-1
      kk = k
      kk1 = k1
      k3 = 3*k+1
      nmin = 2*k1
c  determine the length of the period of s(x).
      per = x(m)-x(1)
      if(iopt.lt.0) go to 50
c  calculation of acc, the absolute tolerance for the root of f(p)=s.
      acc = tol*s
c  determine nmax, the number of knots for periodic spline interpolation
      nmax = m+2*k
      if(s.gt.0.0d0 .or. nmax.eq.nmin) go to 30
c  if s=0, s(x) is an interpolating spline.
      n = nmax
c  test whether the required storage space exceeds the available one.
      if(n.gt.nest) go to 620
c  find the position of the interior knots in case of interpolation.
   5  if((k/2)*2 .eq. k) go to 20
      do 10 i=2,m1
        j = i+k
        t(j) = x(i)
  10  continue
      if(s.gt.0.0d0) go to 50
      kk = k-1
      kk1 = k
      if(kk.gt.0) go to 50
      t(1) = t(m)-per
      t(2) = x(1)
      t(m+1) = x(m)
      t(m+2) = t(3)+per
      do 15 i=1,m1
        c(i) = y(i)
  15  continue
      c(m) = c(1)
      fp = 0.0d0
      fpint(n) = fp0
      fpint(n-1) = 0.0d0
      nrdata(n) = 0
      go to 630
  20  do 25 i=2,m1
        j = i+k
        t(j) = (x(i)+x(i-1))*half
  25  continue
      go to 50
c  if s > 0 our initial choice depends on the value of iopt.
c  if iopt=0 or iopt=1 and s>=fp0, we start computing the least-squares
c  periodic polynomial. (i.e. a constant function).
c  if iopt=1 and fp0>s we start computing the least-squares periodic
c  spline according the set of knots found at the last call of the
c  routine.
  30  if(iopt.eq.0) go to 35
      if(n.eq.nmin) go to 35
      fp0 = fpint(n)
      fpold = fpint(n-1)
      nplus = nrdata(n)
      if(fp0.gt.s) go to 50
c  the case that s(x) is a constant function is treated separetely.
c  find the least-squares constant c1 and compute fp0 at the same time.
  35  fp0 = 0.0d0
      d1 = 0.0d0
      c1 = 0.0d0
      do 40 it=1,m1
        wi = w(it)
        yi = y(it)*wi
        call fpgivs(wi,d1,cos,sin)
        call fprota(cos,sin,yi,c1)
        fp0 = fp0+yi**2
  40  continue
      c1 = c1/d1
c  test whether that constant function is a solution of our problem.
      fpms = fp0-s
      if(fpms.lt.acc .or. nmax.eq.nmin) go to 640
      fpold = fp0
c  test whether the required storage space exceeds the available one.
      if(nmin.ge.nest) go to 620
c  start computing the least-squares periodic spline with one
c  interior knot.
      nplus = 1
      n = nmin+1
      mm = (m+1)/2
      t(k2) = x(mm)
      nrdata(1) = mm-2
      nrdata(2) = m1-mm
c  main loop for the different sets of knots. m is a save upper
c  bound for the number of trials.
  50  do 340 iter=1,m
c  find nrint, the number of knot intervals.
        nrint = n-nmin+1
c  find the position of the additional knots which are needed for
c  the b-spline representation of s(x). if we take
c      t(k+1) = x(1), t(n-k) = x(m)
c      t(k+1-j) = t(n-k-j) - per, j=1,2,...k
c      t(n-k+j) = t(k+1+j) + per, j=1,2,...k
c  then s(x) is a periodic spline with period per if the b-spline
c  coefficients satisfy the following conditions
c      c(n7+j) = c(j), j=1,...k   (**)   with n7=n-2*k-1.
        t(k1) = x(1)
        nk1 = n-k1
        nk2 = nk1+1
        t(nk2) = x(m)
        do 60 j=1,k
          i1 = nk2+j
          i2 = nk2-j
          j1 = k1+j
          j2 = k1-j
          t(i1) = t(j1)+per
          t(j2) = t(i2)-per
  60    continue
c  compute the b-spline coefficients c(j),j=1,...n7 of the least-squares
c  periodic spline sinf(x). the observation matrix a is built up row
c  by row while taking into account condition (**) and is reduced to
c  triangular form by givens transformations .
c  at the same time fp=f(p=inf) is computed.
c  the n7 x n7 triangularised upper matrix a has the form
c            ! a1 '    !
c        a = !    ' a2 !
c            ! 0  '    !
c  with a2 a n7 x k matrix and a1 a n10 x n10 upper triangular
c  matrix of bandwith k+1 ( n10 = n7-k).
c  initialization.
        do 70 i=1,nk1
          z(i) = 0.0d0
          do 70 j=1,kk1
            a1(i,j) = 0.0d0
  70    continue
        n7 = nk1-k
        n10 = n7-kk
        jper = 0
        fp = 0.0d0
        l = k1
        do 290 it=1,m1
c  fetch the current data point x(it),y(it)
          xi = x(it)
          wi = w(it)
          yi = y(it)*wi
c  search for knot interval t(l) <= xi < t(l+1).
  80      if(xi.lt.t(l+1)) go to 85
          l = l+1
          go to 80
c  evaluate the (k+1) non-zero b-splines at xi and store them in q.
  85      call fpbspl(t,n,k,xi,l,h)
          do 90 i=1,k1
            q(it,i) = h(i)
            h(i) = h(i)*wi
  90      continue
          l5 = l-k1
c  test whether the b-splines nj,k+1(x),j=1+n7,...nk1 are all zero at xi
          if(l5.lt.n10) go to 285
          if(jper.ne.0) go to 160
c  initialize the matrix a2.
          do 95 i=1,n7
          do 95 j=1,kk
              a2(i,j) = 0.0d0
  95      continue
          jk = n10+1
          do 110 i=1,kk
            ik = jk
            do 100 j=1,kk1
              if(ik.le.0) go to 105
              a2(ik,i) = a1(ik,j)
              ik = ik-1
 100        continue
 105        jk = jk+1
 110      continue
          jper = 1
c  if one of the b-splines nj,k+1(x),j=n7+1,...nk1 is not zero at xi
c  we take account of condition (**) for setting up the new row
c  of the observation matrix a. this row is stored in the arrays h1
c  (the part with respect to a1) and h2 (the part with
c  respect to a2).
 160      do 170 i=1,kk
            h1(i) = 0.0d0
            h2(i) = 0.0d0
 170      continue
          h1(kk1) = 0.0d0
          j = l5-n10
          do 210 i=1,kk1
            j = j+1
            l0 = j
 180        l1 = l0-kk
            if(l1.le.0) go to 200
            if(l1.le.n10) go to 190
            l0 = l1-n10
            go to 180
 190        h1(l1) = h(i)
            go to 210
 200        h2(l0) = h2(l0)+h(i)
 210      continue
c  rotate the new row of the observation matrix into triangle
c  by givens transformations.
          if(n10.le.0) go to 250
c  rotation with the rows 1,2,...n10 of matrix a.
          do 240 j=1,n10
            piv = h1(1)
            if(piv.ne.0.0d0) go to 214
            do 212 i=1,kk
              h1(i) = h1(i+1)
 212        continue
            h1(kk1) = 0.0d0
            go to 240
c  calculate the parameters of the givens transformation.
 214        call fpgivs(piv,a1(j,1),cos,sin)
c  transformation to the right hand side.
            call fprota(cos,sin,yi,z(j))
c  transformations to the left hand side with respect to a2.
            do 220 i=1,kk
              call fprota(cos,sin,h2(i),a2(j,i))
 220        continue
            if(j.eq.n10) go to 250
            i2 = min0(n10-j,kk)
c  transformations to the left hand side with respect to a1.
            do 230 i=1,i2
              i1 = i+1
              call fprota(cos,sin,h1(i1),a1(j,i1))
              h1(i) = h1(i1)
 230        continue
            h1(i1) = 0.0d0
 240      continue
c  rotation with the rows n10+1,...n7 of matrix a.
 250      do 270 j=1,kk
            ij = n10+j
            if(ij.le.0) go to 270
            piv = h2(j)
            if(piv.eq.0.0d0) go to 270
c  calculate the parameters of the givens transformation.
            call fpgivs(piv,a2(ij,j),cos,sin)
c  transformations to right hand side.
            call fprota(cos,sin,yi,z(ij))
            if(j.eq.kk) go to 280
            j1 = j+1
c  transformations to left hand side.
            do 260 i=j1,kk
              call fprota(cos,sin,h2(i),a2(ij,i))
 260        continue
 270      continue
c  add contribution of this row to the sum of squares of residual
c  right hand sides.
 280      fp = fp+yi**2
          go to 290
c  rotation of the new row of the observation matrix into
c  triangle in case the b-splines nj,k+1(x),j=n7+1,...n-k-1 are all zero
c  at xi.
 285      j = l5
          do 140 i=1,kk1
            j = j+1
            piv = h(i)
            if(piv.eq.0.0d0) go to 140
c  calculate the parameters of the givens transformation.
            call fpgivs(piv,a1(j,1),cos,sin)
c  transformations to right hand side.
            call fprota(cos,sin,yi,z(j))
            if(i.eq.kk1) go to 150
            i2 = 1
            i3 = i+1
c  transformations to left hand side.
            do 130 i1=i3,kk1
              i2 = i2+1
              call fprota(cos,sin,h(i1),a1(j,i2))
 130        continue
 140      continue
c  add contribution of this row to the sum of squares of residual
c  right hand sides.
 150      fp = fp+yi**2
 290    continue
        fpint(n) = fp0
        fpint(n-1) = fpold
        nrdata(n) = nplus
c  backward substitution to obtain the b-spline coefficients c(j),j=1,.n
        call fpbacp(a1,a2,z,n7,kk,c,kk1,nest)
c  calculate from condition (**) the coefficients c(j+n7),j=1,2,...k.
        do 295 i=1,k
          j = i+n7
          c(j) = c(i)
 295    continue
        if(iopt.lt.0) go to 660
c  test whether the approximation sinf(x) is an acceptable solution.
        fpms = fp-s
        if(abs(fpms).lt.acc) go to 660
c  if f(p=inf) < s accept the choice of knots.
        if(fpms.lt.0.0d0) go to 350
c  if n=nmax, sinf(x) is an interpolating spline.
        if(n.eq.nmax) go to 630
c  increase the number of knots.
c  if n=nest we cannot increase the number of knots because of the
c  storage capacity limitation.
        if(n.eq.nest) go to 620
c  determine the number of knots nplus we are going to add.
        npl1 = nplus*2
        rn = nplus
        if(fpold-fp.gt.acc) npl1 = rn*fpms/(fpold-fp)
        nplus = min0(nplus*2,max0(npl1,nplus/2,1))
        fpold = fp
c  compute the sum(wi*(yi-s(xi))**2) for each knot interval
c  t(j+k) <= xi <= t(j+k+1) and store it in fpint(j),j=1,2,...nrint.
        fpart = 0.0d0
        i = 1
        l = k1
        do 320 it=1,m1
          if(x(it).lt.t(l)) go to 300
          new = 1
          l = l+1
 300      term = 0.0d0
          l0 = l-k2
          do 310 j=1,k1
            l0 = l0+1
            term = term+c(l0)*q(it,j)
 310      continue
          term = (w(it)*(term-y(it)))**2
          fpart = fpart+term
          if(new.eq.0) go to 320
          if(l.gt.k2) go to 315
          fpint(nrint) = term
          new = 0
          go to 320
 315      store = term*half
          fpint(i) = fpart-store
          i = i+1
          fpart = store
          new = 0
 320    continue
        fpint(nrint) = fpint(nrint)+fpart
        do 330 l=1,nplus
c  add a new knot
          call fpknot(x,m,t,n,fpint,nrdata,nrint,nest,1)
c  if n=nmax we locate the knots as for interpolation.
          if(n.eq.nmax) go to 5
c  test whether we cannot further increase the number of knots.
          if(n.eq.nest) go to 340
 330    continue
c  restart the computations with the new set of knots.
 340  continue
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  part 2: determination of the smoothing periodic spline sp(x).       c
c  *************************************************************       c
c  we have determined the number of knots and their position.          c
c  we now compute the b-spline coefficients of the smoothing spline    c
c  sp(x). the observation matrix a is extended by the rows of matrix   c
c  b expressing that the kth derivative discontinuities of sp(x) at    c
c  the interior knots t(k+2),...t(n-k-1) must be zero. the corres-     c
c  ponding weights of these additional rows are set to 1/sqrt(p).      c
c  iteratively we then have to determine the value of p such that      c
c  f(p)=sum(w(i)*(y(i)-sp(x(i)))**2) be = s. we already know that      c
c  the least-squares constant function corresponds to p=0, and that    c
c  the least-squares periodic spline corresponds to p=infinity. the    c
c  iteration process which is proposed here, makes use of rational     c
c  interpolation. since f(p) is a convex and strictly decreasing       c
c  function of p, it can be approximated by a rational function        c
c  r(p) = (u*p+v)/(p+w). three values of p(p1,p2,p3) with correspond-  c
c  ing values of f(p) (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s) are used      c
c  to calculate the new value of p such that r(p)=s. convergence is    c
c  guaranteed by taking f1>0 and f3<0.                                 c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  evaluate the discontinuity jump of the kth derivative of the
c  b-splines at the knots t(l),l=k+2,...n-k-1 and store in b.
 350  call fpdisc(t,n,k2,b,nest)
c  initial value for p.
      p1 = 0.0d0
      f1 = fp0-s
      p3 = -one
      f3 = fpms
      n11 = n10-1
      n8 = n7-1
      p = 0.0d0
      l = n7
      do 352 i=1,k
         j = k+1-i
         p = p+a2(l,j)
         l = l-1
         if(l.eq.0) go to 356
 352  continue
      do 354 i=1,n10
         p = p+a1(i,1)
 354  continue
 356  rn = n7
      p = rn/p
      ich1 = 0
      ich3 = 0
c  iteration process to find the root of f(p) = s.
      do 595 iter=1,maxit
c  form the matrix g  as the matrix a extended by the rows of matrix b.
c  the rows of matrix b with weight 1/p are rotated into
c  the triangularised observation matrix a.
c  after triangularisation our n7 x n7 matrix g takes the form
c            ! g1 '    !
c        g = !    ' g2 !
c            ! 0  '    !
c  with g2 a n7 x (k+1) matrix and g1 a n11 x n11 upper triangular
c  matrix of bandwidth k+2. ( n11 = n7-k-1)
        pinv = one/p
c  store matrix a into g
        do 360 i=1,n7
          c(i) = z(i)
          g1(i,k1) = a1(i,k1)
          g1(i,k2) = 0.0d0
          g2(i,1) = 0.0d0
          do 360 j=1,k
            g1(i,j) = a1(i,j)
            g2(i,j+1) = a2(i,j)
 360    continue
        l = n10
        do 370 j=1,k1
          if(l.le.0) go to 375
          g2(l,1) = a1(l,j)
          l = l-1
 370    continue
 375    do 540 it=1,n8
c  fetch a new row of matrix b and store it in the arrays h1 (the part
c  with respect to g1) and h2 (the part with respect to g2).
          yi = 0.0d0
          do 380 i=1,k1
            h1(i) = 0.0d0
            h2(i) = 0.0d0
 380      continue
          h1(k2) = 0.0d0
          if(it.gt.n11) go to 420
          l = it
          l0 = it
          do 390 j=1,k2
            if(l0.eq.n10) go to 400
            h1(j) = b(it,j)*pinv
            l0 = l0+1
 390      continue
          go to 470
 400      l0 = 1
          do 410 l1=j,k2
            h2(l0) = b(it,l1)*pinv
            l0 = l0+1
 410      continue
          go to 470
 420      l = 1
          i = it-n10
          do 460 j=1,k2
            i = i+1
            l0 = i
 430        l1 = l0-k1
            if(l1.le.0) go to 450
            if(l1.le.n11) go to 440
            l0 = l1-n11
            go to 430
 440        h1(l1) = b(it,j)*pinv
            go to 460
 450        h2(l0) = h2(l0)+b(it,j)*pinv
 460      continue
          if(n11.le.0) go to 510
c  rotate this row into triangle by givens transformations without
c  square roots.
c  rotation with the rows l,l+1,...n11.
 470      do 500 j=l,n11
            piv = h1(1)
c  calculate the parameters of the givens transformation.
            call fpgivs(piv,g1(j,1),cos,sin)
c  transformation to right hand side.
            call fprota(cos,sin,yi,c(j))
c  transformation to the left hand side with respect to g2.
            do 480 i=1,k1
              call fprota(cos,sin,h2(i),g2(j,i))
 480        continue
            if(j.eq.n11) go to 510
            i2 = min0(n11-j,k1)
c  transformation to the left hand side with respect to g1.
            do 490 i=1,i2
              i1 = i+1
              call fprota(cos,sin,h1(i1),g1(j,i1))
              h1(i) = h1(i1)
 490        continue
            h1(i1) = 0.0d0
 500      continue
c  rotation with the rows n11+1,...n7
 510      do 530 j=1,k1
            ij = n11+j
            if(ij.le.0) go to 530
            piv = h2(j)
c  calculate the parameters of the givens transformation
            call fpgivs(piv,g2(ij,j),cos,sin)
c  transformation to the right hand side.
            call fprota(cos,sin,yi,c(ij))
            if(j.eq.k1) go to 540
            j1 = j+1
c  transformation to the left hand side.
            do 520 i=j1,k1
              call fprota(cos,sin,h2(i),g2(ij,i))
 520        continue
 530      continue
 540    continue
c  backward substitution to obtain the b-spline coefficients
c  c(j),j=1,2,...n7 of sp(x).
        call fpbacp(g1,g2,c,n7,k1,c,k2,nest)
c  calculate from condition (**) the b-spline coefficients c(n7+j),j=1,.
        do 545 i=1,k
          j = i+n7
          c(j) = c(i)
 545    continue
c  computation of f(p).
        fp = 0.0d0
        l = k1
        do 570 it=1,m1
          if(x(it).lt.t(l)) go to 550
          l = l+1
 550      l0 = l-k2
          term = 0.0d0
          do 560 j=1,k1
            l0 = l0+1
            term = term+c(l0)*q(it,j)
 560      continue
          fp = fp+(w(it)*(term-y(it)))**2
 570    continue
c  test whether the approximation sp(x) is an acceptable solution.
        fpms = fp-s
        if(abs(fpms).lt.acc) go to 660
c  test whether the maximal number of iterations is reached.
        if(iter.eq.maxit) go to 600
c  carry out one more step of the iteration process.
        p2 = p
        f2 = fpms
        if(ich3.ne.0) go to 580
        if((f2-f3) .gt. acc) go to 575
c  our initial choice of p is too large.
        p3 = p2
        f3 = f2
        p = p*con4
        if(p.le.p1) p = p1*con9 +p2*con1
        go to 595
 575    if(f2.lt.0.0d0) ich3 = 1
 580    if(ich1.ne.0) go to 590
        if((f1-f2) .gt. acc) go to 585
c  our initial choice of p is too small
        p1 = p2
        f1 = f2
        p = p/con4
        if(p3.lt.0.0d0) go to 595
        if(p.ge.p3) p = p2*con1 +p3*con9
        go to 595
 585    if(f2.gt.0.0d0) ich1 = 1
c  test whether the iteration process proceeds as theoretically
c  expected.
 590    if(f2.ge.f1 .or. f2.le.f3) go to 610
c  find the new value for p.
        p = fprati(p1,f1,p2,f2,p3,f3)
 595  continue
c  error codes and messages.
 600  ier = 3
      go to 660
 610  ier = 2
      go to 660
 620  ier = 1
      go to 660
 630  ier = -1
      go to 660
 640  ier = -2
c  the least-squares constant function c1 is a solution of our problem.
c  a constant function is a spline of degree k with all b-spline
c  coefficients equal to that constant c1.
      do 650 i=1,k1
        rn = k1-i
        t(i) = x(1)-rn*per
        c(i) = c1
        j = i+k1
        rn = i-1.0d0
        t(j) = x(m)+rn*per
 650  continue
      n = nmin
      fp = fp0
      fpint(n) = fp0
      fpint(n-1) = 0.0d0
      nrdata(n) = 0
 660  return
      end
      double precision function fprati(p1,f1,p2,f2,p3,f3)
c  given three points (p1,f1),(p2,f2) and (p3,f3), function fprati
c  gives the value of p such that the rational interpolating function
c  of the form r(p) = (u*p+v)/(p+w) equals zero at p.
c  ..
c  ..scalar arguments..
      double precision p1,f1,p2,f2,p3,f3
c  ..local scalars..
      double precision h1,h2,h3,p
c  ..
      if(p3.gt.0.0d0) go to 10
c  value of p in case p3 = infinity.
      p = (p1*(f1-f3)*f2-p2*(f2-f3)*f1)/((f1-f2)*f3)
      go to 20
c  value of p in case p3 ^= infinity.
  10  h1 = f1*(f2-f3)
      h2 = f2*(f3-f1)
      h3 = f3*(f1-f2)
      p = -(p1*p2*h3+p2*p3*h1+p3*p1*h2)/(p1*h1+p2*h2+p3*h3)
c  adjust the value of p1,f1,p3 and f3 such that f1 > 0 and f3 < 0.
  20  if(f2.lt.0.0d0) go to 30
      p1 = p2
      f1 = f2
      go to 40
  30  p3 = p2
      f3 = f2
  40  fprati = p
      return
      end
      subroutine fpdisc(t,n,k2,b,nest)
c  subroutine fpdisc calculates the discontinuity jumps of the kth
c  derivative of the b-splines of degree k at the knots t(k+2)..t(n-k-1)
c  ..scalar arguments..
      integer n,k2,nest
c  ..array arguments..
      double precision t(n),b(nest,k2)
c  ..local scalars..
      double precision an,fac,prod
      integer i,ik,j,jk,k,k1,l,lj,lk,lmk,lp,nk1,nrint
c  ..local array..
      double precision h(12)
c  ..
      k1 = k2-1
      k = k1-1
      nk1 = n-k1
      nrint = nk1-k
      an = nrint
      fac = an/(t(nk1+1)-t(k1))
      do 40 l=k2,nk1
        lmk = l-k1
        do 10 j=1,k1
          ik = j+k1
          lj = l+j
          lk = lj-k2
          h(j) = t(l)-t(lk)
          h(ik) = t(l)-t(lj)
  10    continue
        lp = lmk
        do 30 j=1,k2
          jk = j
          prod = h(j)
          do 20 i=1,k
            jk = jk+1
            prod = prod*h(jk)*fac
  20      continue
          lk = lp+k1
          b(lmk,j) = (t(lk)-t(lp))/prod
          lp = lp+1
  30    continue
  40  continue
      return
      end
      subroutine fpknot(x,m,t,n,fpint,nrdata,nrint,nest,istart)
c  subroutine fpknot locates an additional knot for a spline of degree
c  k and adjusts the corresponding parameters,i.e.
c    t     : the position of the knots.
c    n     : the number of knots.
c    nrint : the number of knotintervals.
c    fpint : the sum of squares of residual right hand sides
c            for each knot interval.
c    nrdata: the number of data points inside each knot interval.
c  istart indicates that the smallest data point at which the new knot
c  may be added is x(istart+1)
c  ..
c  ..scalar arguments..
      integer m,n,nrint,nest,istart
c  ..array arguments..
      double precision x(m),t(nest),fpint(nest)
      integer nrdata(nest)
c  ..local scalars..
      double precision an,am,fpmax
      integer ihalf,j,jbegin,jj,jk,jpoint,k,maxbeg,maxpt,
     * next,nrx,number
c  ..
      k = (n-nrint-1)/2
c  search for knot interval t(number+k) <= x <= t(number+k+1) where
c  fpint(number) is maximal on the condition that nrdata(number)
c  not equals zero.
      fpmax = 0.0d0
      jbegin = istart
      do 20 j=1,nrint
        jpoint = nrdata(j)
        if(fpmax.ge.fpint(j) .or. jpoint.eq.0) go to 10
        fpmax = fpint(j)
        number = j
        maxpt = jpoint
        maxbeg = jbegin
  10    jbegin = jbegin+jpoint+1
  20  continue
c  let coincide the new knot t(number+k+1) with a data point x(nrx)
c  inside the old knot interval t(number+k) <= x <= t(number+k+1).
      ihalf = maxpt/2+1
      nrx = maxbeg+ihalf
      next = number+1
      if(next.gt.nrint) go to 40
c  adjust the different parameters.
      do 30 j=next,nrint
        jj = next+nrint-j
        fpint(jj+1) = fpint(jj)
        nrdata(jj+1) = nrdata(jj)
        jk = jj+k
        t(jk+1) = t(jk)
  30  continue
  40  nrdata(number) = ihalf-1
      nrdata(next) = maxpt-ihalf
      am = maxpt
      an = nrdata(number)
      fpint(number) = fpmax*an/am
      an = nrdata(next)
      fpint(next) = fpmax*an/am
      jk = next+k
      t(jk) = x(nrx)
      n = n+1
      nrint = nrint+1
      return
      end
      subroutine fpbacp(a,b,z,n,k,c,k1,nest)
c  subroutine fpbacp calculates the solution of the system of equations
c  g * c = z  with g  a n x n upper triangular matrix of the form
c            ! a '   !
c        g = !   ' b !
c            ! 0 '   !
c  with b a n x k matrix and a a (n-k) x (n-k) upper triangular
c  matrix of bandwidth k1.
c  ..
c  ..scalar arguments..
      integer n,k,k1,nest
c  ..array arguments..
      double precision a(nest,k1),b(nest,k),z(n),c(n)
c  ..local scalars..
      integer i,i1,j,l,l0,l1,n2
      double precision store
c  ..
      n2 = n-k
      l = n
      do 30 i=1,k
        store = z(l)
        j = k+2-i
        if(i.eq.1) go to 20
        l0 = l
        do 10 l1=j,k
          l0 = l0+1
          store = store-c(l0)*b(l,l1)
  10    continue
  20    c(l) = store/b(l,j-1)
        l = l-1
        if(l.eq.0) go to 80
  30  continue
      do 50 i=1,n2
        store = z(i)
        l = n2
        do 40 j=1,k
          l = l+1
          store = store-c(l)*b(i,j)
  40    continue
        c(i) = store
  50  continue
      i = n2
      c(i) = c(i)/a(i,1)
      if(i.eq.1) go to 80
      do 70 j=2,n2
        i = i-1
        store = c(i)
        i1 = k
        if(j.le.k) i1=j-1
        l = i
        do 60 l0=1,i1
          l = l+1
          store = store-c(l)*a(i,l0+1)
  60    continue
        c(i) = store/a(i,1)
  70  continue
  80  return
      end
      subroutine fpbspl(t,n,k,x,l,h)
c  subroutine fpbspl evaluates the (k+1) non-zero b-splines of
c  degree k at t(l) <= x < t(l+1) using the stable recurrence
c  relation of de boor and cox.
c  ..
c  ..scalar arguments..
      double precision x
      integer n,k,l
c  ..array arguments..
      double precision t(n),h(6)
c  ..local scalars..
      double precision f,one
      integer i,j,li,lj
c  ..local arrays..
      double precision hh(5)
c  ..
      one = 0.1d+01
      h(1) = one
      do 20 j=1,k
        do 10 i=1,j
          hh(i) = h(i)
  10    continue
        h(1) = 0.0d0
        do 20 i=1,j
          li = l+i
          lj = li-j
          f = hh(i)/(t(li)-t(lj))
          h(i) = h(i)+f*(t(li)-x)
          h(i+1) = f*(x-t(lj))
  20  continue
      return
      end
      subroutine fprota(cos,sin,a,b)
c  subroutine fprota applies a givens rotation to a and b.
c  ..
c  ..scalar arguments..
      double precision cos,sin,a,b
c ..local scalars..
      double precision stor1,stor2
c  ..
      stor1 = a
      stor2 = b
      b = cos*stor2+sin*stor1
      a = cos*stor1-sin*stor2
      return
      end
      subroutine fpgivs(piv,ww,cos,sin)
c  subroutine fpgivs calculates the parameters of a givens
c  transformation .
c  ..
c  ..scalar arguments..
      double precision piv,ww,cos,sin
c  ..local scalars..
      double precision dd,one,store
c  ..function references..
      double precision abs,sqrt
c  ..
      one = 0.1d+01
      store = abs(piv)
      if(store.ge.ww) dd = store*sqrt(one+(ww/piv)**2)
      if(store.lt.ww) dd = ww*sqrt(one+(piv/ww)**2)
      cos = ww/dd
      sin = piv/dd
      ww = dd
      return
      end
      subroutine fpchep(x,m,t,n,k,ier)
c  subroutine fpchep verifies the number and the position of the knots
c  t(j),j=1,2,...,n of a periodic spline of degree k, in relation to
c  the number and the position of the data points x(i),i=1,2,...,m.
c  if all of the following conditions are fulfilled, ier is set
c  to zero. if one of the conditions is violated ier is set to ten.
c      1) k+1 <= n-k-1 <= m+k-1
c      2) t(1) <= t(2) <= ... <= t(k+1)
c         t(n-k) <= t(n-k+1) <= ... <= t(n)
c      3) t(k+1) < t(k+2) < ... < t(n-k)
c      4) t(k+1) <= x(i) <= t(n-k)
c      5) the conditions specified by schoenberg and whitney must hold
c         for at least one subset of data points, i.e. there must be a
c         subset of data points y(j) such that
c             t(j) < y(j) < t(j+k+1), j=k+1,...,n-k-1
c  ..
c  ..scalar arguments..
      integer m,n,k,ier
c  ..array arguments..
      double precision x(m),t(n)
c  ..local scalars..
      integer i,i1,i2,j,j1,k1,k2,l,l1,l2,mm,m1,nk1,nk2
      double precision per,tj,tl,xi
c  ..
      k1 = k+1
      k2 = k1+1
      nk1 = n-k1
      nk2 = nk1+1
      m1 = m-1
      ier = 50
c  check condition no 1
      if(nk1.lt.k1 .or. n.gt.m+2*k) go to 130
c  check condition no 2
      j = n
      ier = 51
      do 20 i=1,k
        if(t(i).gt.t(i+1)) go to 130
        if(t(j).lt.t(j-1)) go to 130
        j = j-1
  20  continue
c  check condition no 3
      ier = 52
      do 30 i=k2,nk2
        if(t(i).le.t(i-1)) go to 130
  30  continue
c  check condition no 4
      ier = 53
      if(x(1).lt.t(k1) .or. x(m).gt.t(nk2)) go to 130
c  check condition no 5
      ier = 54
      l1 = k1
      l2 = 1
      do 50 l=1,m
         xi = x(l)
  40     if(xi.lt.t(l1+1) .or. l.eq.nk1) go to 50
         l1 = l1+1
         l2 = l2+1
         if(l2.gt.k1) go to 60
         go to 40
  50  continue
      l = m
  60  per = t(nk2)-t(k1)
      do 120 i1=2,l
         i = i1-1
         mm = i+m1
         do 110 j=k1,nk1
            tj = t(j)
            j1 = j+k1
            tl = t(j1)
  70        i = i+1
            if(i.gt.mm) go to 120
            i2 = i-m1
            if(i2) 80,80,90
  80        xi = x(i)
            go to 100
  90        xi = x(i2)+per
 100        if(xi.le.tj) go to 70
            if(xi.ge.tl) go to 120
 110     continue
         ier = 0
         go to 130
 120  continue
 130  return
      end
      subroutine splev(t,n,c,k,x,y,m,ier)
c  subroutine splev evaluates in a number of points x(i),i=1,2,...,m
c  a spline s(x) of degree k, given in its b-spline representation.
c
c  calling sequence:
c     call splev(t,n,c,k,x,y,m,ier)
c
c  input parameters:
c    t    : array,length n, which contains the position of the knots.
c    n    : integer, giving the total number of knots of s(x).
c    c    : array,length n, which contains the b-spline coefficients.
c    k    : integer, giving the degree of s(x).
c    x    : array,length m, which contains the points where s(x) must
c           be evaluated.
c    m    : integer, giving the number of points where s(x) must be
c           evaluated.
c
c  output parameter:
c    y    : array,length m, giving the value of s(x) at the different
c           points.
c    ier  : error flag
c      ier = 0 : normal return
c      ier =10 : invalid input data (see restrictions)
c
c  restrictions:
c    m >= 1
c    t(k+1) <= x(i) <= x(i+1) <= t(n-k) , i=1,2,...,m-1.
c
c  other subroutines required: fpbspl.
c
c  references :
c    de boor c  : on calculating with b-splines, j. approximation theory
c                 6 (1972) 50-62.
c    cox m.g.   : the numerical evaluation of b-splines, j. inst. maths
c                 applics 10 (1972) 134-149.
c    dierckx p. : curve and surface fitting with splines, monographs on
c                 numerical analysis, oxford university press, 1993.
c
c  author :
c    p.dierckx
c    dept. computer science, k.u.leuven
c    celestijnenlaan 200a, b-3001 heverlee, belgium.
c    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
c
c  latest update : march 1987
c
c  ..scalar arguments..
      integer n,k,m,ier
c  ..array arguments..
      double precision t(n),c(n),x(m),y(m)
c  ..local scalars..
      integer i,j,k1,l,ll,l1,nk1
      double precision arg,sp,tb,te
c  ..local array..
      double precision h(6)
c  ..
c  before starting computations a data check is made. if the input data
c  are invalid control is immediately repassed to the calling program.
      ier = 10
      if(m-1) 100,30,10
  10  do 20 i=2,m
        if(x(i).lt.x(i-1)) go to 100
  20  continue
  30  ier = 0
c  fetch tb and te, the boundaries of the approximation interval.
      k1 = k+1
      nk1 = n-k1
      tb = t(k1)
      te = t(nk1+1)
      l = k1
      l1 = l+1
c  main loop for the different points.
      do 80 i=1,m
c  fetch a new x-value arg.
        arg = x(i)
        if(arg.lt.tb) arg = tb
        if(arg.gt.te) arg = te
c  search for knot interval t(l) <= arg < t(l+1)
  40    if(arg.lt.t(l1) .or. l.eq.nk1) go to 50
        l = l1
        l1 = l+1
        go to 40
c  evaluate the non-zero b-splines at arg.
  50    call fpbspl(t,n,k,arg,l,h)
c  find the value of s(x) at x=arg.
        sp = 0.0d0
        ll = l-k1
        do 60 j=1,k1
          ll = ll+1
          sp = sp+c(ll)*h(j)
  60    continue
        y(i) = sp
  80  continue
 100  return
      end
