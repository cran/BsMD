C>>>>>>>>>>>>>>>>>>>>  bm.f  <<<<<<<<<<<<<<<<<<<<<<<<<C
      SUBROUTINE bm(X,Y,N,COLS,BLKS,MXFAC,MXINT,PI,INDGAM,INDG2,GAM2,
     *                NGAM,GAMMA,NTOP,
     *omdcnt,optop,osigtop,onftop,ojtop,odel,osprob,opgam,oprob,
     *ind)

      INTEGER N,COLS,BLKS,MXFAC,MXINT,INDGAM,INDG2,NGAM,NTOP,ind
      DOUBLE PRECISION Y(N),GAMMA(NGAM),GAM2,X(N,(COLS+BLKS))

      INTEGER omdcnt,onftop(ntop),ojtop(ntop,mxfac)
      DOUBLE PRECISION odel,optop(ntop),osigtop(ntop),p0,osprob(cols+1),
     *opgam(ngam),oprob(cols+1,ngam)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Almost no alteration to the original program but READ statements deleted
C     Most of modifications/additions are typed in small letters.
C
C     ind     Indicator variable, 1 subroutine exit properly
C             otherwise it has the format label number.
C     New defined variables have prefix "o" and the original name,
C     e.g., omdcnt = OMDCNT, osprob = SPROB, etc.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


C        PROGRAM MBCQPI5
C        **********************************************************
C        MBCQPI
C        WRITTEN BY R DANIEL MEYER, THE LUBRIZOL CORPORATION,
C        29400 LAKELAND BLVD, WICKLIFFE, OHIO 44092
C        JUNE 1996
C        **********************************************************
C        WHILE MUCH TIME AND EFFORT HAS BEEN PUT INTO THE PRODUCTION
C        OF RELIABLE SOURCE CODE, IT IS IMPOSSIBLE TO GUARANTEE THAT
C        THE CODE IS ERROR FREE.  COMMENTS ON ERRORS, IRRITATIONS,
C        SUGGESTIONS FOR IMPROVEMENT ARE WELCOMED.
C        COPYRIGHT R DANIEL MEYER 1992
C        **********************************************************
C        THE FOLLOWING PARAMETER STATEMENT SHOULD BE EDITED
C        IF YOU WISH TO INCREASE THE SIZE OF PROBLEM THAT CAN
C        BE TACKLED, CHANGE INPUT/OUTPUT UNIT NUMBERS, OR
C        CHANGE THE NUMBER OF LINES ON AN OUTPUT PAGE.
C        THE PARAMETERS HAVE THE FOLLOWING DEFINITIONS:
C
C        LS      NUMBER OF LINES ON AN OUTPUT PAGE
C        MAXN    MAXIMUM NUMBER OF OBSERVATIONS
C        MAXCOL  MAXIMUM NUMBER OF EXPERIMENTAL FACTORS
C        MAXGAM  MAXIMUM NUMBER OF GAMMA VALUES ALLOWED
C        MAXNMD  MAXIMUM NUMBER OF INDIVIDUAL MODELS TO ID
C        MXMXIN  MAXIMUM ORDER INTERACTION
C        BIG     A LARGE WEIGHT;  WHEN AN INDIVIDUAL MODEL HAS WEIGHT
C                GREATER THAN BIG, ALL PROBABILITIES ARE SCALED BACK
C                TO AVOID OVERFLOW
C        MAXTERM MAXIMUM ORDER OF MODEL MATRIX
C        MAXDIM  MAXIMUM NUMBER OF DISTINCT PROBABILITIES TO COMPUTE
C
      PARAMETER (LS=66,
     *           MAXN=100,
     *           MAXCOL=25,
     *           MAXGAM=20,
     *           MAXNMD=100,
     *           MXMXIN=3,
     *           BIG=1.0E10,
     *           MAXTERM=1+MAXCOL+MAXCOL*(MAXCOL-1)/2,
     *           MAXDIM=MAXGAM*MAXCOL)
C
C
C        DESCRIPTION OF VARIABLES
C        ------------------------
C
C        INUNIT  INPUT UNIT NUMBER
C        OUNIT   OUTPUT UNIT NUMBER
C        N       NUMBER OF OBSERVATIONS
C        COLS    NUMBER OF FACTORS
C        BLKS    NUMBER OF BLOCK VARIABLES
C        MXFAC   MAX NUMBER OF FACTORS CONSIDERED IN ANY MODEL
C        MXINT   MAX ORDER INTERACTION CONSIDERED
C        PI      PRIOR PROBABILITY OF AN ACTIVE FACTOR
C        INDGAM  INDICATOR VARIABLE
C                  1 = MORE THAN ONE GAMMA VALUE
C                  0 = ONLY ONE GAMMA VALUE
C        INDG2   INDICATOR VARIABLE, ONLY APPLIES WHEN INDGAM=0
C                  1 = DIFFERENT GAMMA VALUE FOR INTERACTIONS
C                  0 = SAME GAMMA VALUE FOR ALL EFFECTS
C        NGAM    NUMBER OF GAMMA VALUES
C        WCRIT   USER-SPECIFIED VALUE;  ALL MODELS WITH UNSCALED
C                PROBABILITY GREATER THAN WCRIT ARE IDENTIFIED
C        GAMMA   VECTOR OF GAMMA VALUES
C        PGAM    VECTOR OF VALUES OF POSTERIOR DENSITY OF GAMMA
C        GAM2    SCALAR GAMMA VALUE FOR INTERACTIONS
C        GINVSQ  1/(GAMMA**2)
C        G2INSQ  1/(GAM2**2)
C        X       MATRIX OF FACTORS (N BY COLS)
C        Y       VECTOR OF OBSERVATIONS
C        NORM    VECTOR OF CONSTANTS WHICH FORCE PROBS TO SUM TO 1
C        PROB    MATRIX OF MARGINAL POSTERIOR PROBABILITIES
C        SPROB   VECTOR OF POSTERIOR PROBS, AVERAGED OVER GAMMA
C        ST      '*'  (USED FOR FORMATTING OUTPUT)
C        BL      ' '  (USED FOR FORMATTING OUTPUT)
C        MEAN    MEAN OF Y
C        S       RESIDUAL SUM OF SQUARES FOR NULL MODEL
C        SUMNORM SUM OF NORMS FOR EACH GAMMA VALUE
C        PRATIO  PI/(1-PI)
C        EXPON   -(N-1)/2
C        ROOTN   SQRT(N)
C        DEL     INTERVAL WIDTH OF GAMMA PARTITION
C        NFAC    DO LOOP VARIABLE, 1 TO MXFAC
C        ALL     LOGICAL VARIABLE
C                  .TRUE. = ALL MODELS OF SIZE NFAC HAVE BEEN TRIED
C                 .FALSE. = NOT .TRUE.
C        JFAC    INTEGER VECTOR ID'ING FACTORS IN CURRENT MODEL
C        MULT    INTEGER VECTOR ID'ING FACTORS IN CURRENT INTERACTION
C        NTERM   NUMBER OF TERMS IN CURRENT MODEL
C        A       X-MATRIX AUGMENTED WITH INTERACTION COLUMNS
C        PART    LOGICAL VARIABLE
C                  .TRUE. = ALL INTERACTION COLUMNS HAVE BEEN ADDED
C                  .FALSE. = NOT .TRUE.
C        AA      A'A  + (1/GAMMA**2) ADDED TO ALL DIAGONALS BUT 1,1
C        ATEM    DIAG(A'A)
C        LGAM2   LOG10(GAMMA)
C        PIGAM   (PRATIO/GAMMA)**NFAC
C        B       A'Y;  ALSO SOLUTION TO AA*X=A'Y
C        OCOUNT  COUNT OF OUTPUT LINES (USED FOR FORMATTING)
C        DET1
C        DET2    DET(AA)=DET1*(10**DET2)
C        DT      DT(1)=DET1;  DT(2)=DET2
C        RES     RESIDUALS FROM CURRENT MODEL
C        SR      SUM OF SQUARED BAYESIAN RESIDUALS FROM CURRENT MODEL
C        SR2     SUM OF SQUARED RESIDUALS FROM CURRENT MODEL
C        W       UNSCALED POSTERIOR PROB OF CURRENT MODEL
C        MXNORM  MAX(NORMS FOR EACH GAMMA VALUE)
C        BAR     NUMBER OF *'S IN BAR PLOT
C        OLOOP   NUMBER OF SETS OF 12 GAMMA VALUES (FOR FORMATTING)
C        CC      CARRIAGE CONTROL CHARACTER
C        Z       WORK SPACE VECTOR FOR DPOCO SUBROUTINE
C        INFO    ERROR INDICATOR FROM SUBROUTINE DPOCO
C        NTOP    USER-SPECIFIED NUMBER OF INDIVIDUAL MODELS TO ID
C        PTOP    VECTOR OF PROBABILITIES OF TOP NTOP MODELS
C        JTOP    MATRIX OF FACTOR NUMBERS OF TOP NTOP MODELS
C        SIGTOP  VECTOR OF SIGMA-SQUARED OF TOP NTOP MODELS
C        NLOW    INDEX OF LOWEST PROB MODEL IN TOP NTOP MODELS
C        PSCAL   SCALE FACTOR TO KEEP PROBABILITIES WITHIN MACHINE
C                CONSTRAINTS
C        NFTOP   NUMBER OF FACTORS IN EACH OF NTOP MODELS
C        PINDEX  INTEGER VECTOR OF SORTED INDEXES OF PTOP
C        MDCNT   TOTAL NUMBER OF MODELS EVALUATED
C        DET0    DETERMINANT FOR NULL MODEL
C
      LOGICAL PART,ALL
C
      DOUBLE PRECISION S,PRATIO,LGAM2,P,EXPON,PIGAM,SR,W,PI,DET,
     *DET1,DET2,MEAN,ROOTN,GINVSQ,MXNORM,SUMNORM,DEL,COND,WCRIT,PSCAL,
     *G2INSQ,SR2,PGAMMA,PROBZERO,
     *RES(MAXN),
     *A(MAXN,MAXTERM),AA(MAXTERM,MAXTERM),NORM(MAXGAM),
     *B(MAXTERM),PROB(MAXCOL,MAXGAM),SPROB(MAXGAM),PGAM(MAXGAM),
     *Z(MAXTERM),PTOP(MAXNMD),SIGTOP(MAXNMD),PROB0(MAXGAM),
     *ATEM(MAXTERM,MAXTERM),DT(2)
C
      INTEGER NTERM,NFAC,M,I,J,K,II,BAR,OLOOP,
     *OCOUNT,IGAM,ISTART,INUNIT,OUNIT,NLOW,INFO,
     *MULT(MAXCOL),JFAC(MAXCOL),JTOP(MAXNMD,MAXCOL),NFTOP(MAXNMD),
     *PINDEX(MAXNMD),MDCNT
C
      CHARACTER*1 CC,ST,BL
C
C
C
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c initialization to avoid warning
      G2INSQ=0.D0
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      INUNIT=5
      OUNIT=1

      OPEN(OUNIT,FILE="BsPrint.out")

      WRITE(OUNIT,1000)
      OCOUNT=3
      IF ((N .LT. 1) .OR. (N .GT. MAXN)) THEN
        WRITE(OUNIT,1500)
        WRITE(OUNIT,1501) N,MAXN
	    ind = 1501
        GO TO 700
      ENDIF
C
      IF ((COLS .LT. 1) .OR. (COLS .GT. MAXCOL)) THEN
        WRITE(OUNIT,1500)
        WRITE(OUNIT,1502) COLS,MAXCOL
	    ind = 1502
        GO TO 700
      ENDIF
C
      IF ((BLKS .LT. 0) .OR. (BLKS .GT. MAXCOL)) THEN
        WRITE(OUNIT,1500)
        WRITE(OUNIT,1502) BLKS,MAXCOL
	    ind = 1502
        GO TO 700
      ENDIF
C
      IF ((MXFAC .LT. 1) .OR. (MXFAC .GT. COLS)) THEN
        WRITE(OUNIT,1500)
        WRITE(OUNIT,1503) MXFAC,COLS
	    ind = 1503
        GO TO 700
      ENDIF
C
      IF ((MXINT .LT. 1) .OR. (MXINT .GT. MXMXIN)) THEN
        WRITE(OUNIT,1500)
        WRITE(OUNIT,1504) MXINT,MXMXIN
	    ind = 1504
        GO TO 700
      ENDIF
C
      IF ((MXINT .EQ. 3) .AND.
     & ((MXFAC*(MXFAC-1)*(MXFAC-2)/6) .GT. MAXTERM)) THEN
        WRITE(OUNIT,1500)
        WRITE(OUNIT,1505) MXINT,COLS
	    ind = 1505
        GO TO 700
      ENDIF
C
      IF ((PI .LE. 0) .OR. (PI .GE. 1.0)) THEN
        WRITE(OUNIT,1500)
        WRITE(OUNIT,1506) PI
	    ind = 1506
        GO TO 700
      ENDIF
C
      IF ((INDGAM .NE. 0) .AND. (INDGAM .NE. 1)) THEN
        WRITE(OUNIT,1500)
        WRITE(OUNIT,1507) INDGAM
	    ind = 1507
        GO TO 700
      ENDIF
C
      IF (INDGAM .EQ. 0) THEN
       NGAM=1

       IF (GAMMA(1) .LT. 0.0) THEN
         WRITE(OUNIT,1500)
         WRITE(OUNIT,1508) GAMMA(1)
	     ind = 15081
         GO TO 700
       ENDIF

       IF ((NTOP .LT. 0) .OR. (NTOP .GT. MAXNMD)) THEN
         WRITE(OUNIT,1500)
         WRITE(OUNIT,1513) NTOP,MAXNMD
         NTOP=MAXNMD
       ENDIF
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       IF (INDG2 .EQ. 1) THEN
c this identity replaces a read statement
C      write(*,*) GAM2
         GAM2 = GAM2
       ELSE
         GAM2 = GAMMA(1)
       ENDIF
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       G2INSQ=1.0/(GAM2*GAM2)
      ELSE

       NTOP=0
       IF ((NGAM .LE. 1) .OR. (NGAM .GT. MAXGAM)) THEN
         WRITE(OUNIT,1500)
         WRITE(OUNIT,1509) NGAM,MAXGAM
	     ind = 1509
         GO TO 700
       ENDIF
       IF (GAMMA(1) .LE. 0.0) THEN
         WRITE(OUNIT,1500)
         WRITE(OUNIT,1508) GAMMA(1)
	     ind = 15082
         GO TO 700
       ENDIF
       IF (GAMMA(NGAM) .LE. 0.0) THEN
         WRITE(OUNIT,1500)
         WRITE(OUNIT,1508) GAMMA(1)
	     ind = 15083
         GO TO 700
       ENDIF
       IF (GAMMA(NGAM) .LE. GAMMA(1)) THEN
         WRITE(OUNIT,1500)
         WRITE(OUNIT,1510) GAMMA(1),GAMMA(NGAM)
	     ind = 1510
         GO TO 700
       ENDIF
      ENDIF
C
C
C  **********************************************************
C  * TO SUPPRESS UNDERFLOW WARNINGS ON IBM MAINFRAME RUNNING*
C  * VS FORTRAN, INCLUDE THE FOLLOWING STATEMENT            *
C  **********************************************************
C     CALL ERRSET(208,0,-1,1)
C
C     ONE-TIME INITIALIZATION
C

CC      DATA NORM /MAXGAM*0.0/
CC      DATA PROB /MAXDIM*0.0/
CC      DATA PTOP /MAXNMD*0.0/
CC      DATA SPROB /MAXGAM*0.0/
CC      DATA PROB0 /MAXGAM*0.0/
ccccccccccccccccccccccccccccccccccccccccccccc
      do 50 i = 1,MAXGAM
        NORM(i) = 0.0D0
        SPROB(i) = 0.0D0
        PROB0(i) = 0.0D0
50    continue
      do 51 i = 1, MAXGAM
      do 51 j = 1, MAXCOL
        PROB(i,j) = 0.0D0
51    continue
      do 52 i = 1,MAXNMD
        PTOP(i) = 0.0D0
52    continue
ccccccccccccccccccccccccccccccccccccccccccccc

      ST='*'
      BL=' '
      PSCAL=1.0
      MEAN=0.0
      PROBZERO=0.0
      S=0.0
      SUMNORM=0.0
      MDCNT=0
C     PTOP(1)=(GAMMA(1)**(-3))*DEXP(-1.0/(GAMMA(1)**2))
C     PTOP(1)=1.0
C     NFTOP(1)=0
      PRATIO=(PI/(1-PI))
      EXPON=-FLOAT(N-1)/2.0
      ROOTN=SQRT(FLOAT(N))
      DO 15 I=1,NTOP
 15     PINDEX(I)=I
C
C
C
      DO 100 M=1,N
 100     MEAN=MEAN+Y(M)                                                 
      MEAN=MEAN/FLOAT(N)
      DO 110 M=1,N                                                      
 110   S=S+(Y(M)-MEAN)**2                                               
C     SIGTOP(1)=S/FLOAT(N-1)
C
C
C
      IF (INDGAM .EQ. 1) THEN
        DEL=(GAMMA(NGAM)-GAMMA(1))/FLOAT(NGAM-1)
        DO 120 I=2,NGAM-1
 120      GAMMA(I)=GAMMA(1)+FLOAT(I-1)*DEL
      ENDIF
C
C
      DO 420 NFAC=0,MXFAC
         ALL= .FALSE.                                                   
C                                                                       
      CALL INITIA(JFAC,NFAC,MAXCOL)                                     
C
 200  IF (.NOT. ALL) THEN
C
C     AUGMENT WITH INTERACTION COLUMNS                                  
C                                                                       
      MDCNT=MDCNT+1
      DO 210 I=1,N                                                      
       A(I,1)=1.0                                                       
      DO 205 J=1,BLKS                                                   
 205   A(I,J+1)=X(I,J)
      DO 210 J=1,NFAC
 210   A(I,BLKS+J+1)=X(I,BLKS+JFAC(J))
      NTERM=NFAC+1+BLKS                                                 
      DO 250 M=2,MIN(MXINT,NFAC)
         CALL INITIA(MULT,M,MAXCOL)                                     
         PART=.FALSE.                                                   
 220     IF (.NOT. PART) THEN                                           
           NTERM=NTERM+1                                                
           DO 230 I=1,N                                                 
 230         A(I,NTERM)=A(I,MULT(1)+1+BLKS)*A(I,MULT(2)+1+BLKS)
           DO 240 II=3,M
             DO 240 I=1,N
 240           A(I,NTERM)=A(I,NTERM)*A(I,MULT(II)+1+BLKS)
           CALL INCREM(MULT,PART,M,NFAC,MAXCOL)
           GO TO 220                                                    
         ENDIF                                                          
 250  CONTINUE
C
C      FORM X-PRIME-X MATRIX
C
      DO 270 I=1,NTERM                                                  
      DO 270 J=I,NTERM                                                  
        AA(I,J)=0.0
        DO 260 M=1,N                                                    
 260      AA(I,J)=AA(I,J)+A(M,I)*A(M,J)
        ATEM(I,J)=AA(I,J)
        ATEM(J,I)=AA(I,J)
 270    AA(J,I)=AA(I,J)
C
C   LOOP THROUGH VARIOUS VALUES OF GAMMA
C
      DO 400 IGAM=1,NGAM
        DO 280 I=1,NTERM
          DO 280 J=1,NTERM
 280        AA(I,J)=ATEM(I,J)
      GINVSQ=1.0/(GAMMA(IGAM)**2)
C     PGAMMA=(GAMMA(IGAM)**(-3))*DEXP(-1.0/(GAMMA(IGAM)**2))
      PGAMMA=1.0
      DO 310 I=2,NTERM
 310    AA(I,I)=AA(I,I)+GINVSQ
      IF (INDGAM .EQ. 0) THEN
          DO 311 I=BLKS+NFAC+2,NTERM                                    
 311        AA(I,I)=AA(I,I)-GINVSQ+G2INSQ
      ENDIF
      LGAM2=DLOG10(GAMMA(IGAM))
      IF (INDGAM .EQ. 0) LGAM2=DLOG10(GAM2)
      PIGAM=((PRATIO/GAMMA(IGAM))**NFAC)/(GAMMA(IGAM)**BLKS)
      DO 320 I=1,NTERM
        B(I)=0.0                                                        
        DO 320 M=1,N                                                    
 320      B(I)=B(I)+A(M,I)*Y(M)                                         
      CALL DPOCO(AA,MAXTERM,NTERM,COND,Z,INFO)                          
      IF ((COND+1.0) .EQ. COND) THEN                                    
       WRITE(OUNIT,1511)                                                
       OCOUNT=OCOUNT+1
       GO TO 400                                                        
      ENDIF
      CALL DPOSL(AA,MAXTERM,NTERM,B)                                    
      CALL DPODI(AA,MAXTERM,NTERM,DT,10)                                
      DET1=DT(1)
      DET2=DT(2)
      DET=ROOTN*SQRT(10.0**(-DET2-2.0*LGAM2*(NTERM-NFAC-1-BLKS)))       
      DET=DET/SQRT(DET1)
      SR=0.0                                                            
      SR2=0.0
      DO 340 M=1,N                                                      
        RES(M)=Y(M)
       DO 330 I=1,NTERM                                                 
 330    RES(M)=RES(M)-A(M,I)*B(I)
      SR2=SR2+RES(M)**2
 340  SR=SR+RES(M)*Y(M)
C     DO 350 I=2,NTERM
C350    SR=SR+(B(I)**2)/(GAMMA(IGAM)**2)                                
C
C
C       EXTRA WRITE STATEMENTS FOR BRAD JONES
C     WRITE(OUNIT,*)
C     WRITE(OUNIT,*) 'FACTORS'
C     WRITE(OUNIT,*) (JFAC(I), I=1, NFAC)
C     WRITE(OUNIT,*) 'S= ',SR,' DETERM= ',
C    &DET/(ROOTN*SQRT(10.0**(-2.0*LGAM2*(NTERM-NFAC-1-BLKS)))),
C    &' PI/GAMMA FACTOR= ',PIGAM*(10.0**(-LGAM2*(NTERM-NFAC-1-BLKS)))
C     WRITE(6,*) 'BETA= ',(B(I),I=1,NTERM)
C
            WCRIT=PSCAL*BIG
            CALL IDLOW(PTOP,MAXNMD,NTOP,NLOW,WCRIT)
            W=DEXP(DLOG(PIGAM)+DLOG(DET)+EXPON*DLOG(SR/S)-DLOG(PSCAL))
C     WRITE(6,*) 'PARTIALLY SCALED PROBABILITY= ',W
            IF (((W*PGAMMA) .GT. WCRIT) .AND. (INDGAM .EQ. 0)) THEN
              PTOP(NLOW)=W*PGAMMA
              SIGTOP(NLOW)=SR/FLOAT(N-1)
              NFTOP(NLOW)=NFAC
              DO 355 I=1,NFAC
 355            JTOP(NLOW,I)=JFAC(I)
            ENDIF
            NORM(IGAM)=NORM(IGAM)+W                                     
            IF (NFAC .EQ. 0) THEN
               PROB0(IGAM)=PROB0(IGAM)+W*PGAMMA
               PROBZERO=PROBZERO+W*PGAMMA
            ENDIF
            DO 360 I=1,NFAC
              SPROB(JFAC(I))=SPROB(JFAC(I))+W*PGAMMA
 360          PROB(JFAC(I),IGAM)=PROB(JFAC(I),IGAM)+W
            IF (W .GT. PSCAL*BIG) THEN
              PSCAL=PSCAL*BIG
              DO 370 I=1,COLS
                SPROB(I)=SPROB(I)/BIG
                DO 370 J=1,NGAM
 370             PROB(I,J)=PROB(I,J)/BIG
              DO 380 I=1,NTOP
 380           PTOP(I)=PTOP(I)/BIG
              DO 390 I=1,NGAM
 390           NORM(I)=NORM(I)/BIG
            ENDIF
C
C     END COMPUTE                                                       
C                                                                       
 400     CONTINUE
           CALL INCREM(JFAC,ALL,NFAC,COLS,MAXCOL)
         GO TO 200
      ENDIF
 420  CONTINUE


      WRITE(OUNIT,1101)
      OCOUNT=OCOUNT+1
      OCOUNT=MOD(OCOUNT,LS)
      IF (COLS .LE. 15) CALL OSPACE(LS,OCOUNT,N+3,CC)
      IF (COLS .GT. 15) CALL OSPACE(LS,OCOUNT,2*N+3,CC)
      WRITE(OUNIT,1001) CC
      DO 500 I=1,N
         WRITE(OUNIT,1002) I,(X(I,J), J=1,MIN0(15,COLS+BLKS))
         IF ((COLS+BLKS) .GT. 15)
     &   WRITE(OUNIT,1115) (X(I,J), J=16,COLS+BLKS)
 500  CONTINUE
      WRITE(OUNIT,1101)
      CALL OSPACE(LS,OCOUNT,N+3,CC)
      WRITE(OUNIT,1003) CC
      DO 510 I=1,N
        WRITE(OUNIT,1004) I,Y(I)                                        
 510  CONTINUE
      WRITE(OUNIT,1101)
      CALL OSPACE(LS,OCOUNT,6,CC)
      WRITE(OUNIT,1100) CC                                              
      WRITE(OUNIT,1104) N,COLS,BLKS,PI,MXINT,MXFAC,MDCNT
      IF (INDGAM .EQ. 0) WRITE(OUNIT,1117) GAMMA(1),GAM2,NORM(1)
      IF (INDGAM .EQ. 1) WRITE(OUNIT,1118) GAMMA(1),GAMMA(NGAM),DEL
      WRITE(OUNIT,1101)                                                 
      MXNORM=0.0
C     PROB0=0.0
      DO 530 I=1,NGAM
C       PGAMMA=(GAMMA(I)**(-3))*DEXP(-1.0/(GAMMA(I)**2))
        PGAMMA=1.0
        PGAM(I)=NORM(I)*PGAMMA
C       PROB0=PROB0+PGAMMA
        SUMNORM=SUMNORM+PGAM(I)
        DO 520 J=1,COLS
 520      PROB(J,I)=PROB(J,I)/NORM(I)
 530  MXNORM=DMAX1(MXNORM,PGAM(I))
        DO 540 J=1,COLS
 540      SPROB(J)=SPROB(J)/SUMNORM
      WRITE(OUNIT,1101)                                                 
      IF (INDGAM .EQ. 0) THEN
        DO 545 I=1,NTOP
 545      PTOP(I)=PTOP(I)/SUMNORM
          CALL SSORT(PTOP,PINDEX,NTOP,-2)
        IF (NTOP .GT. MDCNT) NTOP=MDCNT
        CALL OSPACE(LS,OCOUNT,NTOP+4,CC)
        WRITE(OUNIT,1105) CC,NTOP                                       
        DO 547 I=1,NTOP
 547   WRITE(OUNIT,1106) PTOP(I),SIGTOP(PINDEX(I)),NFTOP(PINDEX(I)),
     &(JTOP(PINDEX(I),J),J=1,NFTOP(PINDEX(I)))
       WRITE(OUNIT,1101)                                                
      ENDIF
      IF (INDGAM .EQ. 1) THEN
        CALL OSPACE(LS,OCOUNT,NGAM+2,CC)
        WRITE(OUNIT,1102) CC                                            
        DO 550 I=1,NGAM                                                 
         P=PGAM(I)/MXNORM
         BAR=IDNINT((P-DMOD(P,.05D0))/.05D0)
 550     WRITE(OUNIT,1103)
     &    GAMMA(I),PGAM(I),(ST,J=1,BAR),(BL,J=1,20-BAR)
        WRITE(OUNIT,1101)
      ENDIF
      CALL OSPACE(LS,OCOUNT,COLS+5,CC)
      IF (INDGAM .EQ. 0) WRITE(OUNIT,1107) CC
      IF (INDGAM .EQ. 1) WRITE(OUNIT,1116) CC
      P=DINT((1000.0D0*PROBZERO/SUMNORM)/PSCAL)/1000.0D0
      p0 = P
      BAR=IDNINT((P-DMOD(P,.05D0))/.05D0)
      WRITE(OUNIT,1114) P,(ST,J=1,BAR),(BL,J=1,20-BAR)
      DO 560 I=1,COLS
         P=DINT(1000.0D0*SPROB(I))/1000.0D0
         BAR=IDNINT((P-DMOD(P,.05D0))/.05D0)
 560     WRITE(OUNIT,1109) I,SPROB(I),(ST,J=1,BAR),(BL,J=1,20-BAR)
      WRITE(OUNIT,1101)
      IF (INDGAM .EQ. 1) THEN
        OLOOP=(NGAM-MOD(NGAM,12))/12+1
        DO 610 K=1,OLOOP
          ISTART=1+(K-1)*12
          CALL OSPACE(LS,OCOUNT,COLS+5,CC)
          WRITE(OUNIT,1110) CC
          WRITE(OUNIT,1111) (GAMMA(I),I=ISTART,MIN0(ISTART+11,NGAM))
          WRITE(OUNIT,1112)
     &     ((PROB0(I)/NORM(I))/PSCAL, I=ISTART,MIN0(ISTART+11,NGAM))
          DO 600 I=1,COLS
 600    WRITE(OUNIT,1113) I,(PROB(I,J),J=ISTART,MIN0(ISTART+11,NGAM))
          WRITE(OUNIT,1101)
610    CONTINUE
      ENDIF

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C Whatever is printed is an output subroutine parameter
      omdcnt = MDCNT
      odel = DEL
      do 450 i = 1,NTOP
        optop(i) = PTOP(i)
        osigtop(i) = SIGTOP(PINDEX(i))
        onftop(i) = NFTOP(PINDEX(i))
        do 451 j= 1,MXFAC
451        ojtop(i,j)=0
        do 452 j= 1,NFTOP(PINDEX(i))
452        ojtop(i,j) = JTOP(PINDEX(I),J)
450   continue
      osprob(1) = p0
      do 456 i = 1,COLS
456   osprob(i+1) = SPROB(i)
      do 457 i = 1,NGAM
457   opgam(i) = PGAM(i)
      do 458 j = 1, ngam
458   oprob(1,j) = (PROB0(j)/NORM(j))/PSCAL
      do 459 i = 1, cols
         do 459 j = 1, ngam
459      oprob(i+1,j) = prob(i,j)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      ind = 1

 700  CONTINUE

      CLOSE(OUNIT)

 1000 FORMAT(' MBCQPI5: BAYESIAN ANALYSIS OF CONFOUNDED DATA',/,
     &1X,'WRITTEN BY R. DANIEL MEYER, THE LUBRIZOL CORPORATION',/
     &1X,'ALL RIGHTS RESERVED;    JUNE 1996')
 1001 FORMAT(A1,'X-MATRIX',/,' --------')
 1002 FORMAT(' ',I3,15(1X,F7.3))
 1003 FORMAT(A1,'Y-VECTOR',/,' --------')
 1004 FORMAT(' ',I3,1X,F10.4)
 1005 FORMAT(1X,'WEIGHT= ',E16.6E2,' FACTORS: ',30(1X,I2))
 1100 FORMAT(A1,'NO. OF',3X,'NO. OF',2X,'NO. OF',12X,'MAX ORDER',5X,
     &'MAX NO. OF',5X,'TOTAL NO. OF',/,' RUNS',5X,'FACTORS',2X,'BLOCKS',
     &4X,'PI',4X,'INTERACTION',3X,'ACTIVE FACTORS',5X,'MODELS')
 1101 FORMAT(' ',100('-'))
 1102 FORMAT(A1,5X,'GAMMA',12X,'PGAM')
 1103 FORMAT(1X,F10.3,E16.6E2,2X,'+',20(A1),'+')
 1104 FORMAT(1X,I3,I10,I8,F10.3,I9,I15,I18)
 1105 FORMAT(A1,' BEST ',I3,' MODELS',//,1X,'PROBABILITY   SIGMA-SQ ',
     &' NO OF FACTORS   FACTORS')
 1106 FORMAT(1X,F10.6,F12.4,I6,10X,30(1X,I2))
 1107 FORMAT(A1,5X,'POSTERIOR PROBABILITIES',
     &/,' FACTOR',5X,'POST. PROB.')
 1108 FORMAT(' NONE',F14.3)
 1109 FORMAT(' ',I4,F14.3,2X,'+',20(A1),'+')
 1110 FORMAT(A1,11X,'POSTERIOR PROBABILITIES FOR EACH GAMMA VALUE',/)
 1111 FORMAT(1X,'FACTOR',12(F10.2))
 1112 FORMAT(' NONE',2X,12(F10.3))
 1113 FORMAT(' ',I4,2X,12(F10.3))
 1114 FORMAT(' NONE',F14.3,2X,'+',20(A1),'+')
 1115 FORMAT(' ',3X,15(1X,F7.3))
 1116 FORMAT(A1,'POSTERIOR PROBABILITIES WEIGHT-AVERAGED OVER GAMMA',
     &/,' FACTOR',5X,'POST. PROB.')
 1117 FORMAT(/,1X,'GAMMA= ',F8.3,'   GAMMA2= ',F8.3,'  NORM= ',E16.6E2)
 1118 FORMAT(/,1X,'GAMMA= ',F8.3,' TO ',F8.3,' BY ',F6.4,' INCREMENTS')
 1500 FORMAT(' ***** ERROR *****')
 1501 FORMAT(1X,'N=',I8,' OUT OF RANGE',/,
     &' N MUST BE BETWEEN 1 AND ',I4,/)
 1502 FORMAT(1X,'NO. OF FACTORS = ',I8,' OUT OF RANGE',/,
     &' MUST BE BETWEEN 1 AND ',I4,/)
 1503 FORMAT(1X,'MAX. NO. OF FACTORS = ',I8,' OUT OF RANGE',/,
     &' MUST BE BETWEEN 1 AND TOTAL NO. OF FACTORS= ',I4,/)
 1504 FORMAT(1X,'MAX. ORDER INTERACTION = ',I8,' OUT OF RANGE',/,
     &' MUST BE BETWEEN 1 AND ',I4,/)
 1505 FORMAT(1X,'MAX. ORDER INTERACTION = ',I8,' RESULTS IN TOO',/,
     &' MANY COLUMNS FOR NO. OF FACTORS = ',I4,/)
 1506 FORMAT(1X,'PI = ',F8.4,' OUT OF RANGE',/,
     &' MUST BE BETWEEN 0 AND 1',/)
 1507 FORMAT(1X,'GAMMA INDICATOR = ',I6,' OUT OF RANGE',/,
     &' MUST BE EITHER 0 OR 1',/)
 1508 FORMAT(1X,'GAMMA  = ',I6,' OUT OF RANGE',/,
     &' MUST BE POSITIVE',/)
 1509 FORMAT(1X,'NO. OF GAMMAS FOR SEARCH  = ',I6,' OUT OF RANGE',/,
     &' MUST BE BETWEEN 2 AND ',I5,/)
 1510 FORMAT(1X,'GAMMA(FIRST) = ',F12.4,' GAMMA(LAST)= ',F12.4,/,
     &' GAMMA(FIRST) MUST BE LESS THAN GAMMA(LAST)',/)
 1511 FORMAT(' **** WARNING: SINGULAR MATRIX ENCOUNTERED ****')
 1512 FORMAT(' **** WARNING: MAX NUMBER OF FACTORS TOO LARGE ****',
     &/,' CORRECTIVE ACTION: VALUE REDUCED FROM',I3,' TO',I3,/)
 1513 FORMAT(' **** WARNING: NUMBER OF INDIVIDUAL MODELS TOO BIG ****',
     &/,' CORRECTIVE ACTION: VALUE REDUCED FROM',I4,' TO',I4,/)
 1305 FORMAT(1X,'DET:',E16.6E3,' SR:',E16.6E3,' S:',E16.6E3,
     &' EXPON:',E16.6E3)
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE IDLOW(PTOP,MAXNMD,NTOP,K,P)
      DOUBLE PRECISION PTOP(MAXNMD),P
      INTEGER NTOP,I,K
      DO 10 I=1,NTOP
        IF (PTOP(NTOP-I+1) .LT. P) THEN
          K=NTOP-I+1
          P=PTOP(NTOP-I+1)
        ENDIF
 10   CONTINUE
      RETURN
      END
C
      SUBROUTINE OSPACE(LS,OCOUNT,K,CC)
      INTEGER OCOUNT,K
      CHARACTER*1 CC
      IF (((LS-OCOUNT) .LT. K) .AND. (K .LT. LS)) THEN
        OCOUNT=K
        CC='1'
      ELSE
        OCOUNT=OCOUNT+K+1
        CC='0'
      ENDIF                                                             
      RETURN                                                            
      END                                                               
C
C                                                                       
      SUBROUTINE INCREM(J,ALL,R,N,MAXCOL)                               
C                                                                       
C       INCREMENTS THE INTEGER VECTOR J TO THE NEXT VECTOR IN           
C       LEXICAL ORDER
C                                                                       
      INTEGER M,L,R,N,J(MAXCOL)                                         
      LOGICAL OK,ALL                                                    
      L=R                                                               
      ALL=.FALSE.
      OK=.FALSE.                                                        
 50   IF ((.NOT. OK) .AND. (L .GT. 0)) THEN
         IF (J(L) .LT. N-R+L) THEN                                      
            J(L)=J(L)+1
            DO 101 M=L+1,R                                              
 101          J(M)=J(M-1)+1
            OK=.TRUE.                                                   
         ELSE                                                           
            L=L-1                                                       
         ENDIF                                                          
         GO TO 50                                                       
      ENDIF                                                             
      IF (L .LE. 0) ALL=.TRUE.                                          
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE INITIA(J,R,MAXCOL)                                     
C                                                                       
C       INITIATES THE INTEGER VECTOR J TO THE FIRST VALUE               
C       IN LEXICAL ORDER (1,2,3,...,R)                                  
C                                                                       
        INTEGER J(MAXCOL),R,I                                           
         DO 401 I=1,R                                                   
 401        J(I)=I                                                      
         DO 402 I=R+1,MAXCOL                                            
 402        J(I)=0
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE DPOCO(A,LDA,N,RCOND,Z,INFO)
      INTEGER LDA,N,INFO                                                
      DOUBLE PRECISION A(LDA,1),Z(1)                                    
      DOUBLE PRECISION RCOND                                            
C                                                                       
C     DPOCO FACTORS A DOUBLE PRECISION SYMMETRIC POSITIVE DEFINITE      
C     MATRIX AND ESTIMATES THE CONDITION OF THE MATRIX.                 
C
C     IF  RCOND  IS NOT NEEDED, DPOFA IS SLIGHTLY FASTER.               
C     TO SOLVE  A*X = B , FOLLOW DPOCO BY DPOSL.                        
C     TO COMPUTE  INVERSE(A)*C , FOLLOW DPOCO BY DPOSL.                 
C     TO COMPUTE  DETERMINANT(A) , FOLLOW DPOCO BY DPODI.               
C     TO COMPUTE  INVERSE(A) , FOLLOW DPOCO BY DPODI.                   
C                                                                       
C     ON ENTRY                                                          
C                                                                       
C        A       DOUBLE PRECISION(LDA, N)                               
C                THE SYMMETRIC MATRIX TO BE FACTORED.  ONLY THE
C                DIAGONAL AND UPPER TRIANGLE ARE USED.                  
C
C        LDA     INTEGER                                                
C                THE LEADING DIMENSION OF THE ARRAY  A .                
C
C        N       INTEGER                                                
C                THE ORDER OF THE MATRIX  A .                           
C                                                                       
C     ON RETURN                                                         
C                                                                       
C        A       AN UPPER TRIANGULAR MATRIX  R  SO THAT  A = TRANS(R)*R 
C                WHERE  TRANS(R)  IS THE TRANSPOSE.                     
C                THE STRICT LOWER TRIANGLE IS UNALTERED.                
C                IF  INFO .NE. 0 , THE FACTORIZATION IS NOT COMPLETE.   
C                                                                       
C        RCOND   DOUBLE PRECISION                                       
C                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A .        
C                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS       
C                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE             
C                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND . 
C                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION     
C                           1.0 + RCOND .EQ. 1.0                        
C                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING           
C                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF         
C                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE
C                UNDERFLOWS.  IF INFO .NE. 0 , RCOND IS UNCHANGED.      
C                                                                       
C        Z       DOUBLE PRECISION(N)
C                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.  
C                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS
C                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT           
C                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .                    
C                IF  INFO .NE. 0 , Z  IS UNCHANGED.                     
C                                                                       
C        INFO    INTEGER                                                
C                = 0  FOR NORMAL RETURN.                                
C                = K  SIGNALS AN ERROR CONDITION.  THE LEADING MINOR
C                     OF ORDER  K  IS NOT POSITIVE DEFINITE.            
C                                                                       
C     LINPACK.  THIS VERSION DATED 08/14/78 .                           
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.      
C                                                                       
C     SUBROUTINES AND FUNCTIONS                                         
C                                                                       
C     LINPACK DPOFA                                                     
C     BLAS DAXPY,DDOT,DSCAL,DASUM                                       
C     FORTRAN DABS,DMAX1,DREAL,DSIGN
C                                                                       
C     INTERNAL VARIABLES                                                
C                                                                       
      DOUBLE PRECISION DDOT,EK,T,WK,WKM
      DOUBLE PRECISION ANORM,S,DASUM,SM,YNORM                           
      INTEGER I,J,JM1,K,KB,KP1                                          
C                                                                       
C                                                                       
C     FIND NORM OF A USING ONLY UPPER HALF                              
C                                                                       
      DO 30 J = 1, N                                                    
         Z(J) = DASUM(J,A(1,J),1)                                       
         JM1 = J - 1                                                    
         IF (JM1 .LT. 1) GO TO 20                                       
         DO 10 I = 1, JM1                                               
            Z(I) = Z(I) + DABS(A(I,J))                                  
   10    CONTINUE                                                       
   20    CONTINUE                                                       
   30 CONTINUE                                                          
      ANORM = 0.0D0                                                     
      DO 40 J = 1, N                                                    
         ANORM = DMAX1(ANORM,Z(J))                                      
   40 CONTINUE                                                          
C
C     FACTOR                                                            
C                                                                       
      CALL DPOFA(A,LDA,N,INFO)                                          
      IF (INFO .NE. 0) GO TO 180
C                                                                       
C        RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
C        ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  A*Y = E .      
C        THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL        
C        GROWTH IN THE ELEMENTS OF W  WHERE  TRANS(R)*W = E .           
C        THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW.         
C                                                                       
C        SOLVE TRANS(R)*W = E                                           
C
         EK = 1.0D0                                                     
         DO 50 J = 1, N                                                 
            Z(J) = 0.0D0                                                
   50    CONTINUE                                                       
         DO 110 K = 1, N                                                
            IF (Z(K) .NE. 0.0D0) EK = DSIGN(EK,-Z(K))                   
            IF (DABS(EK-Z(K)) .LE. A(K,K)) GO TO 60                     
               S = A(K,K)/DABS(EK-Z(K))
               CALL DSCAL(N,S,Z,1)                                      
               EK = S*EK
   60       CONTINUE                                                    
            WK = EK - Z(K)                                              
            WKM = -EK - Z(K)
            S = DABS(WK)                                                
            SM = DABS(WKM)                                              
            WK = WK/A(K,K)                                              
            WKM = WKM/A(K,K)                                            
            KP1 = K + 1                                                 
            IF (KP1 .GT. N) GO TO 100                                   
               DO 70 J = KP1, N                                         
                  SM = SM + DABS(Z(J)+WKM*A(K,J))                       
                  Z(J) = Z(J) + WK*A(K,J)                               
                  S = S + DABS(Z(J))                                    
   70          CONTINUE                                                 
               IF (S .GE. SM) GO TO 90                                  
                  T = WKM - WK                                          
                  WK = WKM                                              
                  DO 80 J = KP1, N                                      
                     Z(J) = Z(J) + T*A(K,J)                             
   80             CONTINUE                                              
   90          CONTINUE                                                 
  100       CONTINUE                                                    
            Z(K) = WK
  110    CONTINUE                                                       
         S = 1.0D0/DASUM(N,Z,1)                                         
         CALL DSCAL(N,S,Z,1)                                            
C                                                                       
C        SOLVE R*Y = W                                                  
C                                                                       
         DO 130 KB = 1, N
            K = N + 1 - KB                                              
            IF (DABS(Z(K)) .LE. A(K,K)) GO TO 120                       
               S = A(K,K)/DABS(Z(K))                                    
               CALL DSCAL(N,S,Z,1)                                      
  120       CONTINUE                                                    
            Z(K) = Z(K)/A(K,K)                                          
            T = -Z(K)
            CALL DAXPY(K-1,T,A(1,K),1,Z(1),1)                           
  130    CONTINUE                                                       
         S = 1.0D0/DASUM(N,Z,1)                                         
         CALL DSCAL(N,S,Z,1)                                            
C                                                                       
         YNORM = 1.0D0
C                                                                       
C        SOLVE TRANS(R)*V = Y                                           
C                                                                       
         DO 150 K = 1, N
            Z(K) = Z(K) - DDOT(K-1,A(1,K),1,Z(1),1)                     
            IF (DABS(Z(K)) .LE. A(K,K)) GO TO 140
               S = A(K,K)/DABS(Z(K))                                    
               CALL DSCAL(N,S,Z,1)                                      
               YNORM = S*YNORM                                          
  140       CONTINUE                                                    
            Z(K) = Z(K)/A(K,K)                                          
  150    CONTINUE                                                       
         S = 1.0D0/DASUM(N,Z,1)                                         
         CALL DSCAL(N,S,Z,1)                                            
         YNORM = S*YNORM                                                
C                                                                       
C        SOLVE R*Z = V                                                  
C                                                                       
         DO 170 KB = 1, N                                               
            K = N + 1 - KB                                              
            IF (DABS(Z(K)) .LE. A(K,K)) GO TO 160                       
               S = A(K,K)/DABS(Z(K))                                    
               CALL DSCAL(N,S,Z,1)                                      
               YNORM = S*YNORM                                          
  160       CONTINUE                                                    
            Z(K) = Z(K)/A(K,K)
            T = -Z(K)                                                   
            CALL DAXPY(K-1,T,A(1,K),1,Z(1),1)                           
  170    CONTINUE                                                       
C        MAKE ZNORM = 1.0                                               
         S = 1.0D0/DASUM(N,Z,1)                                         
         CALL DSCAL(N,S,Z,1)                                            
         YNORM = S*YNORM                                                
C
         IF (ANORM .NE. 0.0D0) RCOND = YNORM/ANORM                      
         IF (ANORM .EQ. 0.0D0) RCOND = 0.0D0                            
  180 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE DPODI(A,LDA,N,DET,JOB)                                 
      INTEGER LDA,N,JOB
      DOUBLE PRECISION A(LDA,1)                                         
      DOUBLE PRECISION DET(2)                                           
C                                                                       
C     DPODI COMPUTES THE DETERMINANT AND INVERSE OF A CERTAIN
C     DOUBLE PRECISION SYMMETRIC POSITIVE DEFINITE MATRIX (SEE BELOW)   
C     USING THE FACTORS COMPUTED BY DPOCO, DPOFA OR DQRDC.              
C                                                                       
C     ON ENTRY                                                          
C                                                                       
C        A       DOUBLE PRECISION(LDA, N)
C                THE OUTPUT  A  FROM DPOCO OR DPOFA
C                OR THE OUTPUT  X  FROM DQRDC.                          
C                                                                       
C        LDA     INTEGER                                                
C                THE LEADING DIMENSION OF THE ARRAY  A .                
C                                                                       
C        N       INTEGER                                                
C                THE ORDER OF THE MATRIX  A .                           
C                                                                       
C        JOB     INTEGER                                                
C                = 11   BOTH DETERMINANT AND INVERSE.                   
C                = 01   INVERSE ONLY.                                   
C                = 10   DETERMINANT ONLY.                               
C                                                                       
C     ON RETURN                                                         
C                                                                       
C        A       IF DPOCO OR DPOFA WAS USED TO FACTOR  A  THEN          
C                DPODI PRODUCES THE UPPER HALF OF INVERSE(A) .          
C                IF DQRDC WAS USED TO DECOMPOSE  X  THEN                
C                DPODI PRODUCES THE UPPER HALF OF INVERSE(TRANS(X)*X)   
C                WHERE TRANS(X) IS THE TRANSPOSE.
C                ELEMENTS OF  A  BELOW THE DIAGONAL ARE UNCHANGED.      
C                IF THE UNITS DIGIT OF JOB IS ZERO,  A  IS UNCHANGED.   
C                                                                       
C        DET     DOUBLE PRECISION(2)                                    
C                DETERMINANT OF  A  OR OF  TRANS(X)*X  IF REQUESTED.    
C                OTHERWISE NOT REFERENCED.                              
C                DETERMINANT = DET(1) * 10.0**DET(2)                    
C                WITH  1.0 .LE. DET(1) .LT. 10.0                        
C                OR  DET(1) .EQ. 0.0 .
C                                                                       
C     ERROR CONDITION                                                   
C                                                                       
C        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS     
C        A ZERO ON THE DIAGONAL AND THE INVERSE IS REQUESTED.           
C        IT WILL NOT OCCUR IF THE SUBROUTINES ARE CALLED CORRECTLY      
C        AND IF DPOCO OR DPOFA HAS SET INFO .EQ. 0 .
C                                                                       
C     LINPACK.  THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.      
C                                                                       
C     SUBROUTINES AND FUNCTIONS                                         
C                                                                       
C     BLAS DAXPY,DSCAL                                                  
C     FORTRAN MOD                                                       
C                                                                       
C     INTERNAL VARIABLES
C                                                                       
      DOUBLE PRECISION T                                                
      DOUBLE PRECISION S                                                
      INTEGER I,J,JM1,K,KP1                                             
C                                                                       
C     COMPUTE DETERMINANT                                               
C                                                                       
      IF (JOB/10 .EQ. 0) GO TO 70                                       
         DET(1) = 1.0D0                                                 
         DET(2) = 0.0D0                                                 
         S = 10.0D0                                                     
         DO 50 I = 1, N                                                 
            DET(1) = A(I,I)**2*DET(1)                                   
C        ...EXIT                                                        
            IF (DET(1) .EQ. 0.0D0) GO TO 60                             
   10       IF (DET(1) .GE. 1.0D0) GO TO 20                             
               DET(1) = S*DET(1)                                        
               DET(2) = DET(2) - 1.0D0                                  
            GO TO 10                                                    
   20       CONTINUE
   30       IF (DET(1) .LT. S) GO TO 40                                 
               DET(1) = DET(1)/S                                        
               DET(2) = DET(2) + 1.0D0                                  
            GO TO 30                                                    
   40       CONTINUE                                                    
   50    CONTINUE                                                       
   60    CONTINUE                                                       
   70 CONTINUE                                                          
C                                                                       
C     COMPUTE INVERSE(R)
C                                                                       
      IF (MOD(JOB,10) .EQ. 0) GO TO 140                                 
         DO 100 K = 1, N                                                
            A(K,K) = 1.0D0/A(K,K)                                       
            T = -A(K,K)                                                 
            CALL DSCAL(K-1,T,A(1,K),1)                                  
            KP1 = K + 1
            IF (N .LT. KP1) GO TO 90                                    
            DO 80 J = KP1, N                                            
               T = A(K,J)                                               
               A(K,J) = 0.0D0                                           
               CALL DAXPY(K,T,A(1,K),1,A(1,J),1)                        
   80       CONTINUE                                                    
   90       CONTINUE                                                    
  100    CONTINUE                                                       
C
C        FORM  INVERSE(R) * TRANS(INVERSE(R))
C                                                                       
         DO 130 J = 1, N                                                
            JM1 = J - 1                                                 
            IF (JM1 .LT. 1) GO TO 120                                   
            DO 110 K = 1, JM1                                           
               T = A(K,J)                                               
               CALL DAXPY(K,T,A(1,J),1,A(1,K),1)                        
  110       CONTINUE                                                    
  120       CONTINUE                                                    
            T = A(J,J)                                                  
            CALL DSCAL(J,T,A(1,J),1)                                    
  130    CONTINUE                                                       
  140 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE DPOSL(A,LDA,N,B)                                       
      INTEGER LDA,N                                                     
      DOUBLE PRECISION A(LDA,1),B(1)                                    
C
C     DPOSL SOLVES THE DOUBLE PRECISION SYMMETRIC POSITIVE DEFINITE     
C     SYSTEM A * X = B                                                  
C     USING THE FACTORS COMPUTED BY DPOCO OR DPOFA.                     
C                                                                       
C     ON ENTRY                                                          
C                                                                       
C        A       DOUBLE PRECISION(LDA, N)                               
C                THE OUTPUT FROM DPOCO OR DPOFA.                        
C                                                                       
C        LDA     INTEGER                                                
C                THE LEADING DIMENSION OF THE ARRAY  A .
C                                                                       
C        N       INTEGER                                                
C                THE ORDER OF THE MATRIX  A .                           
C                                                                       
C        B       DOUBLE PRECISION(N)
C                THE RIGHT HAND SIDE VECTOR.                            
C
C     ON RETURN                                                         
C                                                                       
C        B       THE SOLUTION VECTOR  X .                               
C                                                                       
C     ERROR CONDITION                                                   
C                                                                       
C        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS     
C        A ZERO ON THE DIAGONAL.  TECHNICALLY THIS INDICATES
C        SINGULARITY BUT IT IS USUALLY CAUSED BY IMPROPER SUBROUTINE    
C        ARGUMENTS.  IT WILL NOT OCCUR IF THE SUBROUTINES ARE CALLED
C        CORRECTLY AND  INFO .EQ. 0 .                                   
C                                                                       
C     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX                 
C     WITH  P  COLUMNS                                                  
C           CALL DPOCO(A,LDA,N,RCOND,Z,INFO)                            
C           IF (RCOND IS TOO SMALL .OR. INFO .NE. 0) GO TO ...          
C           DO 10 J = 1, P                                              
C              CALL DPOSL(A,LDA,N,C(1,J))                               
C        10 CONTINUE                                                    
C                                                                       
C     LINPACK.  THIS VERSION DATED 08/14/78 .                           
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.      
C                                                                       
C     SUBROUTINES AND FUNCTIONS                                         
C                                                                       
C     BLAS DAXPY,DDOT                                                   
C                                                                       
C     INTERNAL VARIABLES
C                                                                       
      DOUBLE PRECISION DDOT,T                                           
      INTEGER K,KB                                                      
C                                                                       
C     SOLVE TRANS(R)*Y = B                                              
C                                                                       
      DO 10 K = 1, N                                                    
         T = DDOT(K-1,A(1,K),1,B(1),1)                                  
         B(K) = (B(K) - T)/A(K,K)                                       
   10 CONTINUE                                                          
C                                                                       
C     SOLVE R*X = Y
C                                                                       
      DO 20 KB = 1, N                                                   
         K = N + 1 - KB
         B(K) = B(K)/A(K,K)                                             
         T = -B(K)                                                      
         CALL DAXPY(K-1,T,A(1,K),1,B(1),1)                              
   20 CONTINUE
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE DPOFA(A,LDA,N,INFO)                                    
      INTEGER LDA,N,INFO                                                
      DOUBLE PRECISION A(LDA,1)                                         
C
C     DPOFA FACTORS A DOUBLE PRECISION SYMMETRIC POSITIVE DEFINITE      
C     MATRIX.                                                           
C
C     DPOFA IS USUALLY CALLED BY DPOCO, BUT IT CAN BE CALLED            
C     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.          
C     (TIME FOR DPOCO) = (1 + 18/N)*(TIME FOR DPOFA) .                  
C                                                                       
C     ON ENTRY                                                          
C                                                                       
C        A       DOUBLE PRECISION(LDA, N)                               
C                THE SYMMETRIC MATRIX TO BE FACTORED.  ONLY THE         
C                DIAGONAL AND UPPER TRIANGLE ARE USED.                  
C                                                                       
C        LDA     INTEGER                                                
C                THE LEADING DIMENSION OF THE ARRAY  A .                
C                                                                       
C        N       INTEGER                                                
C                THE ORDER OF THE MATRIX  A .                           
C                                                                       
C     ON RETURN
C                                                                       
C        A       AN UPPER TRIANGULAR MATRIX  R  SO THAT  A = TRANS(R)*R 
C                WHERE  TRANS(R)  IS THE TRANSPOSE.                     
C                THE STRICT LOWER TRIANGLE IS UNALTERED.                
C                IF  INFO .NE. 0 , THE FACTORIZATION IS NOT COMPLETE.   
C                                                                       
C        INFO    INTEGER                                                
C                = 0  FOR NORMAL RETURN.                                
C                = K  SIGNALS AN ERROR CONDITION.  THE LEADING MINOR    
C                     OF ORDER  K  IS NOT POSITIVE DEFINITE.            
C                                                                       
C     LINPACK.  THIS VERSION DATED 08/14/78 .                           
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS                                         
C                                                                       
C     BLAS DDOT                                                         
C     FORTRAN DSQRT                                                     
C                                                                       
C     INTERNAL VARIABLES
C                                                                       
      DOUBLE PRECISION DDOT,T                                           
      DOUBLE PRECISION S                                                
      INTEGER J,JM1,K                                                   
C     BEGIN BLOCK WITH ...EXITS TO 40                                   
C
C                                                                       
         DO 30 J = 1, N                                                 
            INFO = J                                                    
            S = 0.0D0
            JM1 = J - 1                                                 
            IF (JM1 .LT. 1) GO TO 20                                    
            DO 10 K = 1, JM1                                            
               T = A(K,J) - DDOT(K-1,A(1,K),1,A(1,J),1)                 
               T = T/A(K,K)                                             
               A(K,J) = T                                               
               S = S + T*T                                              
   10       CONTINUE                                                    
   20       CONTINUE                                                    
            S = A(J,J) - S                                              
C       ....EXIT                                                        
            IF (S .LE. 0.0D0) GO TO 40                                  
            A(J,J) = DSQRT(S)                                           
   30    CONTINUE                                                       
         INFO = 0                                                       
   40 CONTINUE
      RETURN                                                            
      END                                                               
C                                                                       
      DOUBLE PRECISION FUNCTION DASUM(N,DX,INCX)                        
C                                                                       
C     RETURNS SUM OF MAGNITUDES OF DOUBLE PRECISION DX.                 
C     DASUM = SUM FROM 0 TO N-1 OF DABS(DX(1+I*INCX))                   
C                                                                       
      DOUBLE PRECISION DX(1)                                            
      DASUM = 0.D0                                                      
      IF(N.LE.0)RETURN                                                  
      IF(INCX.EQ.1)GOTO 20                                              
C
C        CODE FOR INCREMENTS NOT EQUAL TO 1.
C                                                                       
      NS = N*INCX                                                       
          DO 10 I=1,NS,INCX                                             
          DASUM = DASUM + DABS(DX(I))                                   
   10     CONTINUE                                                      
      RETURN                                                            
C
C        CODE FOR INCREMENTS EQUAL TO 1.                                
C                                                                       
C                                                                       
C        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 6.   
C
   20 M = MOD(N,6)                                                      
      IF( M .EQ. 0 ) GO TO 40                                           
      DO 30 I = 1,M                                                     
         DASUM = DASUM + DABS(DX(I))                                    
   30 CONTINUE
      IF( N .LT. 6 ) RETURN                                             
   40 MP1 = M + 1                                                       
      DO 50 I = MP1,N,6                                                 
         DASUM = DASUM + DABS(DX(I)) + DABS(DX(I+1)) + DABS(DX(I+2))    
     1   + DABS(DX(I+3)) + DABS(DX(I+4)) + DABS(DX(I+5))                
   50 CONTINUE                                                          
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY)                            
C                                                                       
C     OVERWRITE DOUBLE PRECISION DY WITH DOUBLE PRECISION DA*DX + DY.   
C     FOR I = 0 TO N-1, REPLACE  DY(LY+I*INCY) WITH DA*DX(LX+I*INCX) +  
C       DY(LY+I*INCY), WHERE LX = 1 IF INCX .GE. 0, ELSE LX = (-INCX)*N,
C       AND LY IS DEFINED IN A SIMILAR WAY USING INCY.
C                                                                       
      DOUBLE PRECISION DX(1),DY(1),DA                                   
      IF(N.LE.0.OR.DA.EQ.0.D0) RETURN                                   
      IF(INCX.EQ.INCY) IF(INCX-1) 5,20,60                               
    5 CONTINUE                                                          
C                                                                       
C        CODE FOR NONEQUAL OR NONPOSITIVE INCREMENTS.                   
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
C        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 4.   
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
C                                                                       
C        CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS.                  
C                                                                       
   60 CONTINUE                                                          
      NS = N*INCX                                                       
          DO 70 I=1,NS,INCX                                             
          DY(I) = DA*DX(I) + DY(I)                                      
   70     CONTINUE
      RETURN                                                            
      END                                                               
C                                                                       
      DOUBLE PRECISION FUNCTION DDOT(N,DX,INCX,DY,INCY)                 
C                                                                       
C     RETURNS THE DOT PRODUCT OF DOUBLE PRECISION DX AND DY.            
C     DDOT = SUM FOR I = 0 TO N-1 OF  DX(LX+I*INCX) * DY(LY+I*INCY)     
C     WHERE LX = 1 IF INCX .GE. 0, ELSE LX = (-INCX)*N, AND LY IS       
C     DEFINED IN A SIMILAR WAY USING INCY.                              
C                                                                       
      DOUBLE PRECISION DX(1),DY(1)
      DDOT = 0.D0                                                       
      IF(N.LE.0)RETURN                                                  
      IF(INCX.EQ.INCY) IF(INCX-1) 5,20,60                               
    5 CONTINUE                                                          
C
C         CODE FOR UNEQUAL OR NONPOSITIVE INCREMENTS.                   
C                                                                       
      IX = 1                                                            
      IY = 1                                                            
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1                                 
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1                                 
      DO 10 I = 1,N
         DDOT = DDOT + DX(IX)*DY(IY)                                    
        IX = IX + INCX                                                  
        IY = IY + INCY
   10 CONTINUE                                                          
      RETURN                                                            
C                                                                       
C        CODE FOR BOTH INCREMENTS EQUAL TO 1.                           
C                                                                       
C                                                                       
C        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 5.
C                                                                       
   20 M = MOD(N,5)                                                      
      IF( M .EQ. 0 ) GO TO 40                                           
      DO 30 I = 1,M                                                     
         DDOT = DDOT + DX(I)*DY(I)                                      
   30 CONTINUE                                                          
      IF( N .LT. 5 ) RETURN                                             
   40 MP1 = M + 1                                                       
      DO 50 I = MP1,N,5                                                 
         DDOT = DDOT + DX(I)*DY(I) + DX(I+1)*DY(I+1) +                  
     1   DX(I + 2)*DY(I + 2) + DX(I + 3)*DY(I + 3) + DX(I + 4)*DY(I + 4)
   50 CONTINUE                                                          
      RETURN
C                                                                       
C         CODE FOR POSITIVE EQUAL INCREMENTS .NE.1.                     
C                                                                       
   60 CONTINUE                                                          
      NS = N*INCX                                                       
          DO 70 I=1,NS,INCX                                             
          DDOT = DDOT + DX(I)*DY(I)                                     
   70     CONTINUE                                                      
      RETURN                                                            
      END
C                                                                       
      SUBROUTINE DSCAL(N,DA,DX,INCX)                                    
C                                                                       
C     REPLACE DOUBLE PRECISION DX BY DOUBLE PRECISION DA*DX.            
C     FOR I = 0 TO N-1, REPLACE DX(1+I*INCX) WITH  DA * DX(1+I*INCX)    
C                                                                       
      DOUBLE PRECISION DA,DX(1)
      IF(N.LE.0)RETURN                                                  
      IF(INCX.EQ.1)GOTO 20                                              
C                                                                       
C        CODE FOR INCREMENTS NOT EQUAL TO 1.                            
C                                                                       
      NS = N*INCX                                                       
          DO 10 I = 1,NS,INCX
          DX(I) = DA*DX(I)                                              
   10     CONTINUE
      RETURN                                                            
C                                                                       
C        CODE FOR INCREMENTS EQUAL TO 1.                                
C                                                                       
C                                                                       
C        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 5.   
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
C                                                                       
      SUBROUTINE SSORT(X,Y,N,KFLAG)                                     
C***BEGIN PROLOGUE   SSORT                                              
C***REVISION  OCTOBER 1,1980                                            
C***CATEGORY NO.  M1                                                    
C***KEYWORD(S) SORTING,SORT,SINGLETON QUICKSORT,QUICKSORT
C***DATE WRITTEN  NOVEMBER,1976                                         
C***AUTHOR  JONES R.E., WISNIEWSKI J.A. (SLA)                           
C***PURPOSE                                                             
C         SSORT SORTS ARRAY X AND OPTIONALLY MAKES THE SAME             
C         INTERCHANGES IN ARRAY Y.  THE ARRAY X MAY BE SORTED IN        
C         INCREASING ORDER OR DECREASING ORDER.  A SLIGHTLY MODIFIED    
C         QUICKSORT ALGORITHM IS USED.                                  
C***DESCRIPTION                                                         
C     SANDIA MATHEMATICAL PROGRAM LIBRARY
C     APPLIED MATHEMATICS DIVISION 2646                                 
C     SANDIA LABORATORIES                                               
C     ALBUQUERQUE, NEW MEXICO  87185                                    
C     CONTROL DATA 6600/7600  VERSION 8.1  AUGUST 1980                  
C                                                                       
C     WRITTEN BY RONDALL E JONES                                        
C     MODIFIED BY JOHN A. WISNIEWSKI TO USE THE SINGLETON QUICKSORT
C     ALGORITHM. DATE 18 NOVEMBER 1976.
C                                                                       
C     ABSTRACT                                                          
C         SSORT SORTS ARRAY X AND OPTIONALLY MAKES THE SAME             
C         INTERCHANGES IN ARRAY Y.  THE ARRAY X MAY BE SORTED IN        
C         INCREASING ORDER OR DECREASING ORDER.  A SLIGHTLY MODIFIED    
C         QUICKSORT ALGORITHM IS USED.                                  
C                                                                       
C     REFERENCE                                                         
C         SINGLETON,R.C., ALGORITHM 347, AN EFFICIENT ALGORITHM FOR
C         SORTING WITH MINIMAL STORAGE, CACM,12(3),1969,185-7.          
C                                                                       
C     DESCRIPTION OF PARAMETERS                                         
C         X - ARRAY OF VALUES TO BE SORTED   (USUALLY ABSCISSAS)        
C         Y - ARRAY TO BE (OPTIONALLY) CARRIED ALONG                    
C         N - NUMBER OF VALUES IN ARRAY X TO BE SORTED                  
C         KFLAG - CONTROL PARAMETER                                     
C             =2  MEANS SORT X IN INCREASING ORDER AND CARRY Y ALONG.   
C             =1  MEANS SORT X IN INCREASING ORDER (IGNORING Y)         
C             =-1 MEANS SORT X IN DECREASING ORDER (IGNORING Y)         
C             =-2 MEANS SORT X IN DECREASING ORDER AND CARRY Y ALONG.
C                                                                       
C                                                                       
C***REFERENCE(S)                                                        
C         SINGLETON,R.C., ALGORITHM 347, AN EFFICIENT ALGORITHM FOR     
C         SORTING WITH MINIMAL STORAGE, CACM,12(3),1969,185-7.          
C***ROUTINES CALLED  XERROR                                             
C***END PROLOGUE                                                        
C
C     DIMENSION X(N),Y(N),IL(21),IU(21)                                 
      DOUBLE PRECISION X(N)                                             
      INTEGER N,KFLAG,Y(N),IL(21),IU(21)                                
C***FIRST EXECUTABLE STATEMENT    SSORT                                 
      NN = N                                                            
      IF (NN.GE.1) GO TO 10                                             
C     CALL XERROR (58HSSORT- THE NUMBER OF VALUES TO BE SORTED WAS NOT P
C    1OSITIVE.,58,1,1)                                                  
      RETURN                                                            
   10 KK = IABS(KFLAG)                                                  
      IF ((KK.EQ.1).OR.(KK.EQ.2)) GO TO 15
C     CALL XERROR (62HSSORT- THE SORT CONTROL PARAMETER, K, WAS NOT 2, 1
C    1, -1, OR -2.,62,2,1)                                              
      RETURN                                                            
C                                                                       
C ALTER ARRAY X TO GET DECREASING ORDER IF NEEDED                       
C                                                                       
   15 IF (KFLAG.GE.1) GO TO 30
      DO 20 I=1,NN                                                      
   20 X(I) = -X(I)                                                      
   30 GO TO (100,200),KK                                                
C                                                                       
C SORT X ONLY                                                           
C                                                                       
  100 CONTINUE                                                          
      M=1                                                               
      I=1                                                               
      J=NN
      R=.375                                                            
  110 IF (I .EQ. J) GO TO 155                                           
  115 IF (R .GT. .5898437) GO TO 120                                    
      R=R+3.90625E-2                                                    
      GO TO 125                                                         
  120 R=R-.21875                                                        
  125 K=I                                                               
C                                  SELECT A CENTRAL ELEMENT OF THE      
C                                  ARRAY AND SAVE IT IN LOCATION T      
      IJ = I + IFIX (FLOAT (J-I) * R)
      T=X(IJ)                                                           
C                                  IF FIRST ELEMENT OF ARRAY IS GREATER 
C                                  THAN T, INTERCHANGE WITH T           
      IF (X(I) .LE. T) GO TO 130                                        
      X(IJ)=X(I)                                                        
      X(I)=T                                                            
      T=X(IJ)
  130 L=J                                                               
C                                  IF LAST ELEMENT OF ARRAY IS LESS THAN
C                                  T, INTERCHANGE WITH T                
      IF (X(J) .GE. T) GO TO 140                                        
      X(IJ)=X(J)                                                        
      X(J)=T                                                            
      T=X(IJ)                                                           
C                                  IF FIRST ELEMENT OF ARRAY IS GREATER 
C                                  THAN T, INTERCHANGE WITH T           
      IF (X(I) .LE. T) GO TO 140                                        
      X(IJ)=X(I)                                                        
      X(I)=T                                                            
      T=X(IJ)
      GO TO 140                                                         
  135 TT=X(L)                                                           
      X(L)=X(K)                                                         
      X(K)=TT                                                           
C                                  FIND AN ELEMENT IN THE SECOND HALF OF
C                                  THE ARRAY WHICH IS SMALLER THAN T
  140 L=L-1
      IF (X(L) .GT. T) GO TO 140                                        
C                                  FIND AN ELEMENT IN THE FIRST HALF OF 
C                                  THE ARRAY WHICH IS GREATER THAN T    
  145 K=K+1                                                             
      IF (X(K) .LT. T) GO TO 145                                        
C                                  INTERCHANGE THESE ELEMENTS           
      IF (K .LE. L) GO TO 135                                           
C                                  SAVE UPPER AND LOWER SUBSCRIPTS OF   
C                                  THE ARRAY YET TO BE SORTED           
      IF (L-I .LE. J-K) GO TO 150
      IL(M)=I                                                           
      IU(M)=L                                                           
      I=K                                                               
      M=M+1                                                             
      GO TO 160                                                         
  150 IL(M)=K                                                           
      IU(M)=J                                                           
      J=L                                                               
      M=M+1
      GO TO 160                                                         
C                                  BEGIN AGAIN ON ANOTHER PORTION OF    
C                                  THE UNSORTED ARRAY                   
  155 M=M-1                                                             
      IF (M .EQ. 0) GO TO 300                                           
      I=IL(M)
      J=IU(M)                                                           
  160 IF (J-I .GE. 1) GO TO 125                                         
      IF (I .EQ. 1) GO TO 110                                           
      I=I-1                                                             
  165 I=I+1                                                             
      IF (I .EQ. J) GO TO 155                                           
      T=X(I+1)                                                          
      IF (X(I) .LE. T) GO TO 165                                        
      K=I                                                               
  170 X(K+1)=X(K)                                                       
      K=K-1                                                             
      IF (T .LT. X(K)) GO TO 170                                        
      X(K+1)=T                                                          
      GO TO 165                                                         
C
C SORT X AND CARRY Y ALONG                                              
C                                                                       
  200 CONTINUE                                                          
      M=1                                                               
      I=1
      J=NN                                                              
      R=.375
  210 IF (I .EQ. J) GO TO 255                                           
  215 IF (R .GT. .5898437) GO TO 220
      R=R+3.90625E-2                                                    
      GO TO 225                                                         
  220 R=R-.21875                                                        
  225 K=I                                                               
C                                  SELECT A CENTRAL ELEMENT OF THE      
C                                  ARRAY AND SAVE IT IN LOCATION T      
      IJ = I + IFIX (FLOAT (J-I) *R)                                    
      T=X(IJ)
      TY= Y(IJ)                                                         
C                                  IF FIRST ELEMENT OF ARRAY IS GREATER 
C                                  THAN T, INTERCHANGE WITH T           
      IF (X(I) .LE. T) GO TO 230                                        
      X(IJ)=X(I)                                                        
      X(I)=T                                                            
      T=X(IJ)                                                           
       Y(IJ)= Y(I)
       Y(I)=TY                                                          
      TY= Y(IJ)                                                         
  230 L=J                                                               
C                                  IF LAST ELEMENT OF ARRAY IS LESS THAN
C                                  T, INTERCHANGE WITH T
      IF (X(J) .GE. T) GO TO 240                                        
      X(IJ)=X(J)                                                        
      X(J)=T                                                            
      T=X(IJ)                                                           
       Y(IJ)= Y(J)                                                      
       Y(J)=TY                                                          
      TY= Y(IJ)                                                         
C                                  IF FIRST ELEMENT OF ARRAY IS GREATER 
C                                  THAN T, INTERCHANGE WITH T           
      IF (X(I) .LE. T) GO TO 240                                        
      X(IJ)=X(I)                                                        
      X(I)=T                                                            
      T=X(IJ)                                                           
       Y(IJ)= Y(I)                                                      
       Y(I)=TY                                                          
      TY= Y(IJ)                                                         
      GO TO 240
  235 TT=X(L)                                                           
      X(L)=X(K)                                                         
      X(K)=TT                                                           
      TTY= Y(L)
       Y(L)= Y(K)                                                       
       Y(K)=TTY                                                         
C                                  FIND AN ELEMENT IN THE SECOND HALF OF
C                                  THE ARRAY WHICH IS SMALLER THAN T    
  240 L=L-1                                                             
      IF (X(L) .GT. T) GO TO 240                                        
C                                  FIND AN ELEMENT IN THE FIRST HALF OF
C                                  THE ARRAY WHICH IS GREATER THAN T    
  245 K=K+1
      IF (X(K) .LT. T) GO TO 245
C                                  INTERCHANGE THESE ELEMENTS           
      IF (K .LE. L) GO TO 235                                           
C                                  SAVE UPPER AND LOWER SUBSCRIPTS OF
C                                  THE ARRAY YET TO BE SORTED
      IF (L-I .LE. J-K) GO TO 250
      IL(M)=I
      IU(M)=L
      I=K
      M=M+1
      GO TO 260
  250 IL(M)=K
      IU(M)=J
      J=L
      M=M+1
      GO TO 260
C                                  BEGIN AGAIN ON ANOTHER PORTION OF
C                                  THE UNSORTED ARRAY
  255 M=M-1
      IF (M .EQ. 0) GO TO 300
      I=IL(M)
      J=IU(M)
  260 IF (J-I .GE. 1) GO TO 225
      IF (I .EQ. 1) GO TO 210
      I=I-1
  265 I=I+1
      IF (I .EQ. J) GO TO 255
      T=X(I+1)
      TY= Y(I+1)
      IF (X(I) .LE. T) GO TO 265
      K=I
  270 X(K+1)=X(K)
       Y(K+1)= Y(K)
      K=K-1
      IF (T .LT. X(K)) GO TO 270
      X(K+1)=T
       Y(K+1)=TY
      GO TO 265
C
C CLEAN UP
C
  300 IF (KFLAG.GE.1) RETURN
      DO 310 I=1,NN
  310 X(I) = -X(I)
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C>>>>>>>>>>>>>>>>>>>>  mds.f  <<<<<<<<<<<<<<<<<<<<<<<<<C
      SUBROUTINE md(NSTART,NRUNS,ITMAX,INITDES,N0,IND,X0,Y0,
     *GAMMA,GAM2,cBL,cCOLS,N,cX,cNM,cP,cSIGMA2,cNF,MNF,cJFAC,cCUT,
     *MBEST,NTOP,TOPD,TOPDES,EPS,flag)

      INTEGER NSTART,NRUNS,ITMAX,INITDES,N0,IND,MBEST(NSTART,NRUNS)
      INTEGER cBL,cCOLS,N,cNM,cNF(cNM),MNF,cJFAC(cNM,MNF),cCUT
      INTEGER NTOP,TOPDES(NTOP,NRUNS),flag,OUT
      DOUBLE PRECISION GAMMA,GAM2,cP(cNM),cSIGMA2(cNM),EPS
      DOUBLE PRECISION X0(N0,cBL+cCOLS),Y0(N0),cX(N,cBL+cCOLS)
      DOUBLE PRECISION TOPD(NTOP)


      PARAMETER (MAXM=100,MAXCOL=100,MAXN0=40,MAXN1=32,MAXN=256,
     &MXSTRT=100)

      DOUBLE PRECISION BETA(MAXM,MAXCOL),G(MAXM,MAXCOL,MAXCOL)
      DOUBLE PRECISION P(MAXM),X(MAXN,MAXCOL),SIGMA2(MAXM)
      DOUBLE PRECISION DTOP(MXSTRT),XBEST(MAXN1)
      DOUBLE PRECISION D,RCOND,DET(2),DBEST,DSTART
      DOUBLE PRECISION A(MAXN0,MAXCOL),AA(MAXCOL,MAXCOL),Z(MAXCOL)
      DOUBLE PRECISION B(MAXCOL)
      INTEGER IM,NM,I,J,BL,NF(MAXM),II,TK,TOTO,M,CUT
      INTEGER JFAC(MAXM,20),MULT(20),BEST(MAXN1),ROWS(MAXN1),N1,IJ
      INTEGER COLS,INFO,BESTI
      INTEGER TOPROW(MXSTRT,MAXN1),ITOP(MXSTRT),ND,IZ,JJ
      LOGICAL PART,DESIN

      COMMON BETA,G,P,X,SIGMA2,NF,JFAC,BL,CUT,COLS,NM

      OUT = 1
      OPEN(OUT,FILE="MDPrint.out")


      NM = cNM
      BL = cBL
      COLS = cCOLS
      CUT = cCUT

      do 8 i = 1, N
        do 8 j = 1, (BL+COLS)
8       X(i,j) = cX(i,j)
      do 5 i = 1, NM
        NF(i) = cNF(i)
        SIGMA2(i) = cSIGMA2(i)
        P(i) = cP(i)
        do 6 j = 1, MNF
6         JFAC(i,j) = cJFAC(i,j)
5     continue

      WRITE(OUT,800)
      IF (IND .EQ. 1) THEN
         GAM2 = GAM2
      ELSE
         GAM2 = GAMMA
      ENDIF

      WRITE(OUT,1101)
      WRITE(OUT,1000)
      WRITE(OUT,1001) N0,COLS,BL,CUT,GAMMA,GAM2,NM
      WRITE(OUT,1002)
      WRITE(OUT,1003) N,NRUNS,ITMAX,NSTART
      WRITE(OUT,1101)
      WRITE(OUT,1004)
      DO 15 I=1,NM
        WRITE(OUT,1005) I,P(I),SIGMA2(I),NF(I),(JFAC(I,J),J=1,NF(I))
 15   CONTINUE
C  **************************************
C  * SUPPRESS UNDERFLOW WARNINGS ON IBM *
C  **************************************
C     CALL ERRSET(208,0,-1,1)

c     To avoid warnings
      DATA BESTI /0/

C      DATA DTOP /MXSTRT*0.0/
      do 17 i = 1, MXSTRT
17    DTOP(i) = 0.0D0

      WRITE(OUT,1102)
      DO 16 I=1,N
         WRITE(OUT,1103) I,(X(I,J), J=1,COLS+BL)
 16      CONTINUE

C
C     ONE-TIME INITIALIZATION
C
      DO 190 IM=1,NM
C
C     AUGMENT CANDIDATE MATRIX WITH INTERACTION COLUMNS
C
      TOTO=COLS+BL
      DO 920 M=2,CUT
         CALL INITIA2(MULT,M)
         PART=.FALSE.
 925     IF (.NOT. PART) THEN
           TOTO=TOTO+1
           DO 930 I=1,N
 930         X(I,TOTO)=X(I,MULT(1)+BL)*X(I,MULT(2)+BL)
           DO 935 II=3,M
             DO 935 I=1,N
 935           X(I,TOTO)=X(I,TOTO)*X(I,MULT(II)+BL)
           CALL INVREM2(MULT,PART,M,COLS)
           GO TO 925
         ENDIF
 920  CONTINUE
C
      TK=NF(IM)
      DO 110 I=1,N0
       A(I,1)=1.0                                                       
       DO 115 J=1,BL                                                    
 115     A(I,1+J)=X0(I,J)
       DO 110 J=1,TK                                                    
 110     A(I,J+1+BL)=X0(I,JFAC(IM,J)+BL)                                
      TOTO=TK+1+BL
C
C     AUGMENT WITH INTERACTION COLUMNS
C                                                                       
      DO 120 M=2,MIN(CUT,TK)
         CALL INITIA2(MULT,M)
         PART=.FALSE.
 125     IF (.NOT. PART) THEN
           TOTO=TOTO+1
           DO 130 I=1,N0
 130         A(I,TOTO)=A(I,MULT(1)+1+BL)*A(I,MULT(2)+1+BL)
           DO 135 II=3,M
             DO 135 I=1,N0
 135           A(I,TOTO)=A(I,TOTO)*A(I,MULT(II)+1+BL)
           CALL INVREM2(MULT,PART,M,TK)
           GO TO 125
         ENDIF
 120  CONTINUE
C
C      FORM X-PRIME-X MATRIX
C
      DO 140 I=1,TOTO
      DO 140 J=I,TOTO
        AA(I,J)=0.0
        DO 145 M=1,N0
 145      AA(I,J)=AA(I,J)+A(M,I)*A(M,J)
 140   AA(J,I)=AA(I,J)
      DO 150 I=2,TK+1+BL
 150   AA(I,I)=AA(I,I)+1.0/(GAMMA**2)
      DO 155 I=TK+2+BL,TOTO
 155   AA(I,I)=AA(I,I)+1.0/(GAM2**2)
      DO 160 I=1,TOTO
      B(I)=0.0
      DO 160 M=1,N0
 160    B(I)=B(I)+A(M,I)*Y0(M)
C     CALL DLINDS(TOTO,AA,MAXCOL,AA,MAXCOL)
      CALL DPOCO(AA,MAXCOL,TOTO,RCOND,Z,INFO)
      CALL DPODI(AA,MAXCOL,TOTO,DET,1)
      DO 165 I=1,TOTO
        BETA(IM,I)=0.0
        DO 165 J=1,TOTO
          IF (I .GT. J) AA(I,J)=AA(J,I)
          BETA(IM,I)=BETA(IM,I)+AA(I,J)*B(J)
          G(IM,I,J)=AA(I,J)
 165    CONTINUE
 190  CONTINUE
C
C
      IF (ITMAX .EQ. 0) WRITE(OUT,1200)
      NDTOP=1
      IJ=1
      DO 700 ISTART=1,NSTART
      IF (ITMAX .GT. 0) THEN
        WRITE(OUT,1101)
        WRITE(OUT,1199) ISTART
        WRITE(OUT,1200)
      ENDIF
C     CALL RNUND(NRUNS,N,BEST)
      IF (INITDES .GT. 0) THEN
        CALL RANST(NRUNS,N,BEST,0.0D0)
      ELSE
        do 791 I=1,NRUNS
791        BEST(I) = MBEST(ISTART,I)
      ENDIF
C      WRITE(OUT,1211)
      CALL EVAL(NRUNS,BEST,D,NM)
C
C
C
      M=0
      DBEST=D
      IF (DTOP(IJ) .LT. DBEST) THEN
        DESIN=.FALSE.
        DO 207 I=1,NRUNS
 207      XBEST(I)=DFLOAT(BEST(I))
        CALL SSORT(XBEST,BEST,NRUNS,2)
        DO 209 J=1, NDTOP-1
          IF ((DABS((DTOP(J)-DBEST)/DTOP(J)) .LT. 0.00001) .AND.
     &    (DESIN .EQV. .FALSE.)) THEN
            IZ=0
            DO 208 I=1,NRUNS
 208          IZ=IZ+IABS(TOPROW(J,I)-BEST(I))
            IF (IZ .EQ. 0) DESIN=.TRUE.
          ENDIF
 209    CONTINUE
        IF (DESIN .EQV. .FALSE.) THEN
          DTOP(IJ)=DBEST
            DO 210 I=1,NRUNS
              TOPROW(IJ,I)=BEST(I)
 210        CONTINUE
          NDTOP=NDTOP+1
          CALL FINDMIN(NDTOP,IJ,DTOP,MXSTRT)
        ENDIF
      ENDIF
      WRITE(OUT,1201) M,D,(BEST(J), J=1,NRUNS)
      IF (ITMAX .EQ. 0) GO TO 700
C
C    NOW START EXCHANGE ITERATIONS
C
 500  CONTINUE
      DSTART=DBEST
      M=M+1
      N1=NRUNS+1
      DO 410 I=1,NRUNS
 410    ROWS(I)=BEST(I)
C
C     FIRST CYCLE THROUGH THE N POSSIBLE ADDITIONAL POINTS              
C
      DO 450 I=1,N
        ROWS(NRUNS+1)=I
        CALL EVAL(N1,ROWS,D,NM)
        IF (D .GT. DBEST) THEN                                          
          DBEST=D                                                       
          BEST(NRUNS+1)=I
        ENDIF
 450  CONTINUE
      WRITE(OUT,1201) M,DBEST,(BEST(J), J=1,NRUNS+1)
C
C     THEN CYCLE THROUGH THE (NRUNS+1) POSSIBLE DELETED POINTS          
C                                                                       
      N1=NRUNS
      DBEST=DSTART                                                      
      DO 460 I=1,NRUNS+1                                                
        DO 465 J=1,NRUNS
          IF (J .LT. I) ROWS(J)=BEST(J)
          IF (J .GE. I) ROWS(J)=BEST(J+1)                               
 465    CONTINUE                                                        
      CALL EVAL(N1,ROWS,D,NM)
      IF (DTOP(IJ) .LT. D) THEN
        DESIN=.FALSE.
        DO 466 II=1,NRUNS                                                                                         
 466      XBEST(II)=DFLOAT(ROWS(II))                                    
        CALL SSORT(XBEST,ROWS,NRUNS,2)                                  
        DO 468 JJ=1,NDTOP-1
          IF ((DABS((DTOP(JJ)-D)/DTOP(JJ)) .LT. 0.00001) .AND.
     &    (DESIN .EQV. .FALSE.)) THEN
            IZ=0
            DO 467 II=1,NRUNS
 467          IZ=IZ+IABS(TOPROW(JJ,II)-ROWS(II))
            IF (IZ .EQ. 0) DESIN=.TRUE.
          ENDIF
 468    CONTINUE
        IF (DESIN .EQV. .FALSE.) THEN
          DTOP(IJ)=D                                                    
            DO 469 II=1,NRUNS
              TOPROW(IJ,II)=ROWS(II)
 469        CONTINUE
          NDTOP=NDTOP+1
          CALL FINDMIN(NDTOP,IJ,DTOP,MXSTRT)
        ENDIF
      ENDIF
      IF (D .GE. DBEST) THEN
          DBEST=D                                                       
          BESTI=I
      ENDIF                                                             
 460  CONTINUE                                                          
      DELTAD=DBEST-DSTART
      DO 475 J=1,NRUNS
        IF (J .GE. BESTI) BEST(J)=BEST(J+1)
 475  CONTINUE                                                          
      WRITE(OUT,1201) M,DBEST,(BEST(J), J=1,NRUNS)                        
      IF ((DELTAD .GT. EPS).AND. (M .LT. ITMAX)) GO TO 500              
C
C   ITERATIONS ENDED; CONVERGENCE OR MAX ITERATIONS REACHED             
C                                                                       
      IF (DELTAD .LE. EPS) THEN
        WRITE(OUT,1202)
C       CALL SVIGN(NRUNS,BEST,BEST)                                     
C       WRITE(OUT,1204)                                                   
C       DO 690 I=1,NRUNS                                                
C         WRITE(OUT,1205) I,BEST(I),(X(BEST(I),J), J=1,COLS+BL)           
C690    CONTINUE
      ENDIF                                                             
      IF (M .GE. ITMAX) WRITE(OUT,1203)
C     WRITE(OUT,1101)
 700  CONTINUE
      ND=MIN0(MXSTRT,NDTOP-1)                                                 
      DO 701 I=1,ND
 701    ITOP(I)=I                                                       
      CALL SSORT(DTOP,ITOP,ND,-2)                                       
      WRITE(OUT,1209)
      WRITE(OUT,1206) ND
      WRITE(OUT,1209)
      WRITE(OUT,1210)
      DO 710 J=1,ND
        WRITE(OUT,1201) J,DTOP(J),(TOPROW(ITOP(J),K), K=1,NRUNS)
 710  CONTINUE
ccccccccccccc
      NTOP = MIN0(NTOP,ND)
      DO 711 J = 1,NTOP
        TOPD(J) = DTOP(J)
        DO 711 K=1,NRUNS
           TOPDES(J,K) = TOPROW(ITOP(J),K)
711   CONTINUE
ccccccccccccc
      IF (ITMAX .GT. 0) THEN
        DO 720 J=1,ND                                                   
          WRITE(OUT,1207) J,DTOP(J)
          DO 720 I=1,NRUNS
            WRITE(OUT,1205) I,TOPROW(ITOP(J),I),                          
     &           (X(TOPROW(ITOP(J),I),K), K=1,COLS+BL)
 720    CONTINUE                                                        
      ENDIF                                                             
      WRITE(OUT,1006)
 800  FORMAT(7X,' FORTRAN PROGRAM MD: BAYESIAN DESIGN OF EXPERIMENTS',/,                
     &3X,'FOLLOWUP DESIGN / WYNN EXCHANGE / RANDOM START',/,            
     &7X,'WRITTEN BY DAN MEYER',/,7X,'ALL RIGHTS RESERVED',/)
 1000 FORMAT(2X,'          NO OF    NO OF  MAX ORDER',
     &     /,2X,'  N0     FACTORS   BLOCKS INTERACTION  ',
     &'  GAMMA(MAIN)  GAMMA(INT)  NMODELS')                             
 1001 FORMAT(1X,I6,I8,I9,I10,F15.3,F12.3,I12,//)                        
 1002 FORMAT(1X,'NO OF       NO OF   MAX          NO OF RANDOM',/,
     &       1X,'CANDIDATES  RUNS    ITERATIONS   STARTS')
 1003 FORMAT(1X,I5,I9,I9,I12,//)                                        
 1004 FORMAT(2X,'MODEL',8X,'PROB',7X,'SIGSQ',3X,'SIZE',3X,'FACTORS')    
 1005 FORMAT(1X,I6,F12.5,F12.4,I7,3X,12(I4))
 1101 FORMAT(1X,100('-'))
 1102 FORMAT(1X,'CANDIDATE RUNS',/,1X,'--------------')
 1103 FORMAT(1X,I3,2X,12(F5.2,1X))
 1104 FORMAT(1X,I3,2X,I5,F8.5,5X,10(I3,1X))
 1105 FORMAT('1MODEL  SIZE   PROB    FACTORS')
 1006 FORMAT(/,' PROGRAM DONE')
 1199 FORMAT(/,1X,'RANDOM START NUMBER:',I3,/)
 1200 FORMAT(//,5X,'ITERATION    D',6X,'DESIGN(ROWS)',/,
     &5X,9('-'),2X,5('-'),4X,98('-'))
 1203 FORMAT(//,5X,'*** MAX ITERATIONS REACHED ***')
 1204 FORMAT(1X,'RUN  CAND  FACTOR LEVELS',/,1X,35('-'))
 1205 FORMAT(1X,I3,2X,I3,2X,12(F5.2,1X))
 1201 FORMAT(5X,I6,F13.4,2X,24(I4),/,25X,24(I4))
C1202 FORMAT(//,5X,'*** CONVERGENCE ***',//,5X,'DESIGN',/,5X,6('-'),/)
 1202 FORMAT(//,5X,'*** CONVERGENCE ***',//)
 1206 FORMAT(1X,'*  THE ',I3,' BEST DESIGNS  *')
 1207 FORMAT(1X,//,1X,'DESIGN ',I3,' D= ',F13.4,/,
     &       1X,'RUN  CAND  FACTOR LEVELS',/,1X,35('-'))
 1209 FORMAT(1X,26('*'))
 1210 FORMAT(//,8X,'RANK         D',6X,'DESIGN(ROWS)',/,
     &8X,4('-'),7X,5('-'),4X,98('-'))
 1211 FORMAT(1X,'  I  J   P(I)   P(J)  TRACE1  TRACE2 QF(I,J)',
     &' QF(J,I)  N*  TERM(I,J)')

      CLOSE(OUT)

      FLAG = 1
      RETURN
      END
C
      SUBROUTINE FINDMIN(I,J,V,N)
      DOUBLE PRECISION V(N),D
      INTEGER I,J,N,K
      IF (I .LE. N) THEN
        J=I
        RETURN
      ELSE
        D=1.0D20
        DO 100 K=1,N
           IF (V(K) .LT. D) THEN
             J=K
             D=V(K)
           ENDIF
 100      CONTINUE
          RETURN
      ENDIF
      END

C
      SUBROUTINE RANST(N1,N,ROWS,R)
      INTEGER I,N1,N,ROWS(N)
      DOUBLE PRECISION X,R
C
C  THE FUNCTION RAND RETURNS A UNIFORM(0,1) DEVIATE;
C  IF ANOTHER RANDOM NUMBER GENERATOR IS AVAILABLE THAT
C  IS SET UP FOR THE MACHINE THIS IS RUNNING ON, IT CAN BE SUBSTITUTED
C
   	  RR=R
      DO 1 I=1,N1
        X=RANDO(RR)
 1      ROWS(I)=IDINT((N-1)*X)+1
      RETURN
      END
C
C
      SUBROUTINE EVAL(N1,ROWS,D,NM)
      COMMON BETA,G,P,X,SIGMA2,NF,JFAC,BL,CUT,COLS
      DOUBLE PRECISION BETA(100,100),G(100,100,100)
      DOUBLE PRECISION YHAT(100,32),W(32),W1(32)
      DOUBLE PRECISION D,D0,TR,TR1,TR2,DEV,DEV2,DEV1,RCOND,DET(2)
      DOUBLE PRECISION P(100),X(256,100),SIGMA2(100)
      DOUBLE PRECISION A(32,100)
      DOUBLE PRECISION V(32,32),V2(32,32),Z(32)
      DOUBLE PRECISION DV(100,32,32),DV2(100,32,32)
      INTEGER IM,NM,I,J,BL,NF(100),TK,TOTO,M,CUT,I0,I1,I2,CNO
      INTEGER JFAC(100,20),MULT(20),ROWS(32),N1,COLS,INFO
      LOGICAL PART
C
      D=0.0
      DO 210 IM=1,NM
        TK=NF(IM)
C
        DO 215 I=1,N1
          A(I,1)=1.0
          DO 220 J=1,BL
 220        A(I,1+J)=X(ROWS(I),J)
          DO 215 J=1,TK
 215        A(I,J+1+BL)=X(ROWS(I),JFAC(IM,J)+BL)
        TOTO=TK+1+BL
C
C     AUGMENT WITH INTERACTION COLUMNS
C
        DO 225 M=2,MIN(CUT,TK)
          CALL INITIA2(MULT,M)
          PART=.FALSE.
 230      IF (.NOT. PART) THEN
            TOTO=TOTO+1
            I0=JFAC(IM,MULT(1))
            I1=JFAC(IM,MULT(2))
            IF (M .EQ. 2) THEN
              CNO=(I0-1)*COLS-(I0-1)*I0/2+I1-I0+COLS+BL
            ELSE
              I2=JFAC(IM,MULT(3))
              CNO=(COLS-2)*COLS-(COLS-2)*(COLS-1)/2+1+COLS+
     &  ((I0-1)*COLS*COLS-(I0+1)*(I0-1)*COLS+(I0*I0*I0-I0)/3)/2+
     &  (I1-I0-1)*(COLS-I0)-(I1-I0-1)*(I1-I0)/2+I2-I1+BL
            ENDIF
           DO 235 I=1,N1
 235         A(I,TOTO)=X(ROWS(I),CNO)
           CALL INVREM2(MULT,PART,M,TK)
           GO TO 230
         ENDIF
 225  CONTINUE
C
C
C
      DO 240 I=1,TOTO
      DO 240 J=1,N1
         V(I,J)=0.0
       DO 240 M=1,TOTO
 240     V(I,J)=V(I,J)+G(IM,I,M)*A(J,M)
C
      DO 245 I=1,N1
      V2(I,I)=1.0
      DO 245 J=1,N1
         IF (I .NE. J) V2(I,J)=0.0
         DO 245 M=1,TOTO
 245     V2(I,J)=V2(I,J)+A(I,M)*V(M,J)
C
      DO 246 I=1,N1
        DO 246 J=1,N1
 246      V(I,J)=V2(I,J)
C     CALL DLINDS(N1,V2,32,V,32)
      CALL DPOCO(V,32,N1,RCOND,Z,INFO)
      CALL DPODI(V,32,N1,DET,1)
      DO 250 I=1,N1
        DO 250 J=1,N1
          IF (I .GT. J) V(I,J)=V(J,I)
          DV(IM,I,J)=V(I,J)
 250      DV2(IM,I,J)=V2(I,J)
C
      DO 255 I=1,N1
        YHAT(IM,I)=0.0
        DO 255 J=1,TOTO
          YHAT(IM,I)=YHAT(IM,I)+A(I,J)*BETA(IM,J)
 255  continue
 210  CONTINUE
c      return
C
C
      DO 300 IM=1,NM-1
C       WRITE(15,*) 'DES- ',ID,' 1ST MODEL PAIR MEMBER- ',IM
        DO 300 JM=IM+1,NM
C
        TR=0.0
        TR1=0.0
        TR2=0.0
        DO 310 I=1,N1
        DO 310 J=1,N1
         TR1=TR1+0.5*(DV(IM,I,J)*DV2(JM,J,I))
         TR2=TR2+0.5*(DV2(IM,I,J)*DV(JM,J,I))
 310     TR=TR+0.5*(DV(IM,I,J)*DV2(JM,J,I)+DV2(IM,I,J)*DV(JM,J,I))
C
C
      DO 320 I=1,N1
        W(I)=0.0
        W1(I)=0.0
        DO 320 J=1,N1
          W(I)=W(I)+DV(JM,I,J)*(YHAT(IM,J)-YHAT(JM,J))/SIGMA2(IM)
 320      W1(I)=W1(I)+DV(IM,I,J)*(YHAT(IM,J)-YHAT(JM,J))/SIGMA2(JM)
      DEV=0.0
      DEV1=0.0
      DEV2=0.0
      DO 330 I=1,N1
        DEV1=DEV1+(YHAT(IM,I)-YHAT(JM,I))*W1(I)
 330    DEV2=DEV2+(YHAT(IM,I)-YHAT(JM,I))*W(I)
      DEV1=DEV1/2.0
      DEV2=DEV2/2.0
      DEV=(DEV1+DEV2)
C
      D0=P(IM)*P(JM)*(TR+DEV-N1)
      D=D+D0
C      WRITE(OUT,500) IM,JM,P(IM),P(JM),TR1,TR2,DEV1,DEV2,N1,D0
 300  CONTINUE
 500  FORMAT(1X,I3,I3,F7.4,F7.4,F8.2,F8.2,F8.2,F8.2,I3,F8.2)
      RETURN
      END
C
C
C
C
      SUBROUTINE INVREM2(J,ALL,R,N)
      INTEGER M,L,R,N,J(20)
      LOGICAL OK,ALL
      L=R
      ALL=.FALSE.
      OK=.FALSE.
 50   IF ((.NOT. OK) .AND. (L .GT. 0)) THEN
         IF (J(L) .LT. N-R+L) THEN
            J(L)=J(L)+1
            DO 101 M=L+1,R
 101          J(M)=J(M-1)+1
            OK=.TRUE.
         ELSE
            L=L-1
         ENDIF                                                          
         GO TO 50                                                       
      ENDIF                                                             
      IF (L .LE. 0) ALL=.TRUE.                                          
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE INITIA2(J,R)
        INTEGER J(20),R,I                                               
         DO 401 I=1,R                                                   
 401        J(I)=I                                                      
         DO 402 I=R+1,20                                                
 402        J(I)=0                                                      
      RETURN                                                            
      END                                                               
C                                                                       
      FUNCTION RANDO (R)
C APRIL 1977 VERSION.  W. FULLERTON, C3, LOS ALAMOS SCIENTIFIC LAB.
C
C      THIS PSEUDO-RANDOM NUMBER GENERATOR IS PORTABLE AMOUNG A WIDE
C VARIETY OF COMPUTERS.  RAND(R) UNDOUBTEDLY IS NOT AS GOOD AS MANY
C READILY AVAILABLE INSTALLATION DEPENDENT VERSIONS, AND SO THIS
C ROUTINE IS NOT RECOMMENDED FOR WIDESPREAD USAGE.  ITS REDEEMING
C FEATURE IS THAT THE EXACT SAME RANDOM NUMBERS (TO WITHIN FINAL ROUND-
C OFF ERROR) CAN BE GENERATED FROM MACHINE TO MACHINE.  THUS, PROGRAMS
C THAT MAKE USE OF RANDOM NUMBERS CAN BE EASILY TRANSPORTED TO AND      
C CHECKED IN A NEW ENVIRONMENT.                                         
C      THE RANDOM NUMBERS ARE GENERATED BY THE LINEAR CONGRUENTIAL      
C METHOD DESCRIBED, E.G., BY KNUTH IN SEMINUMERICAL METHODS (P.9),      
C ADDISON-WESLEY, 1969.  GIVEN THE I-TH NUMBER OF A PSEUDO-RANDOM
C SEQUENCE, THE I+1 -ST NUMBER IS GENERATED FROM                        
C             X(I+1) = (A*X(I) + C) MOD M,                              
C WHERE HERE M = 2**22 = 4194304, C = 1731 AND SEVERAL SUITABLE VALUES  
C OF THE MULTIPLIER A ARE DISCUSSED BELOW.  BOTH THE MULTIPLIER A AND   
C RANDOM NUMBER X ARE REPRESENTED IN DOUBLE PRECISION AS TWO 11-BIT     
C WORDS.  THE CONSTANTS ARE CHOSEN SO THAT THE PERIOD IS THE MAXIMUM    
C POSSIBLE, 4194304.                                                    
C      IN ORDER THAT THE SAME NUMBERS BE GENERATED FROM MACHINE TO      
C MACHINE, IT IS NECESSARY THAT 23-BIT INTEGERS BE REDUCIBLE MODULO     
C 2**11 EXACTLY, THAT 23-BIT INTEGERS BE ADDED EXACTLY, AND THAT 11-BIT 
C INTEGERS BE MULTIPLIED EXACTLY.  FURTHERMORE, IF THE RESTART OPTION   
C IS USED (WHERE R IS BETWEEN 0 AND 1), THEN THE PRODUCT R*2**22 =      
C R*4194304 MUST BE CORRECT TO THE NEAREST INTEGER.                     
C      THE FIRST FOUR RANDOM NUMBERS SHOULD BE .0004127026,             
C .6750836372, .1614754200, AND .9086198807.  THE TENTH RANDOM NUMBER
C IS .5527787209, AND THE HUNDREDTH IS .3600893021 .  THE THOUSANDTH    
C NUMBER SHOULD BE .2176990509 .                                        
C      IN ORDER TO GENERATE SEVERAL EFFECTIVELY INDEPENDENT SEQUENCES
C WITH THE SAME GENERATOR, IT IS NECESSARY TO KNOW THE RANDOM NUMBER
C FOR SEVERAL WIDELY SPACED CALLS.  THE I-TH RANDOM NUMBER TIMES 2**22, 
C WHERE I=K*P/8 AND P IS THE PERIOD OF THE SEQUENCE (P = 2**22), IS     
C STILL OF THE FORM L*P/8.  IN PARTICULAR WE FIND THE I-TH RANDOM       
C NUMBER MULTIPLIED BY 2**22 IS GIVEN BY                                
C I   =  0  1*P/8  2*P/8  3*P/8  4*P/8  5*P/8  6*P/8  7*P/8  8*P/8      
C RAND=  0  5*P/8  2*P/8  7*P/8  4*P/8  1*P/8  6*P/8  3*P/8  0          
C THUS THE 4*P/8 = 2097152 RANDOM NUMBER IS 2097152/2**22.              
C      SEVERAL MULTIPLIERS HAVE BEEN SUBJECTED TO THE SPECTRAL TEST     
C (SEE KNUTH, P. 82).  FOUR SUITABLE MULTIPLIERS ROUGHLY IN ORDER OF    
C GOODNESS ACCORDING TO THE SPECTRAL TEST ARE
C    3146757 = 1536*2048 + 1029 = 2**21 + 2**20 + 2**10 + 5             
C    2098181 = 1024*2048 + 1029 = 2**21 + 2**10 + 5                     
C    3146245 = 1536*2048 +  517 = 2**21 + 2**20 + 2**9 + 5              
C    2776669 = 1355*2048 + 1629 = 5**9 + 7**7 + 1                       
C                                                                       
C      IN THE TABLE BELOW LOG10(NU(I)) GIVES ROUGHLY THE NUMBER OF      
C RANDOM DECIMAL DIGITS IN THE RANDOM NUMBERS CONSIDERED I AT A TIME.
C C IS THE PRIMARY MEASURE OF GOODNESS.  IN BOTH CASES BIGGER IS BETTER.
C                                                                       
C                   LOG10 NU(I)              C(I)
C       A       I=2  I=3  I=4  I=5    I=2  I=3  I=4  I=5
C
C    3146757    3.3  2.0  1.6  1.3    3.1  1.3  4.6  2.6
C    2098181    3.3  2.0  1.6  1.2    3.2  1.3  4.6  1.7
C    3146245    3.3  2.2  1.5  1.1    3.2  4.2  1.1  0.4
C    2776669    3.3  2.1  1.6  1.3    2.5  2.0  1.9  2.6
C   BEST
C    POSSIBLE   3.3  2.3  1.7  1.4    3.6  5.9  9.7  14.9
C
C             INPUT ARGUMENT --
C R      IF R=0., THE NEXT RANDOM NUMBER OF THE SEQUENCE IS GENERATED.
C        IF R.LT.0., THE LAST GENERATED NUMBER WILL BE RETURNED FOR
C          POSSIBLE USE IN A RESTART PROCEDURE.
C        IF R.GT.0., THE SEQUENCE OF RANDOM NUMBERS WILL START WITH THE
C          SEED R MOD 1.  THIS SEED IS ALSO RETURNED AS THE VALUE OF
C          RAND PROVIDED THE ARITHMETIC IS DONE EXACTLY.
C
C             OUTPUT VALUE --
C RAND   A PSEUDO-RANDOM NUMBER BETWEEN 0. AND 1.
C
C IA1 AND IA0 ARE THE HI AND LO PARTS OF A.  IA1MA0 = IA1 - IA0.
      DATA IA1, IA0, IA1MA0 /1536, 1029, 507/
      DATA IC /1731/
      DATA IX1, IX0 /0, 0/
C
      IF (R.LT.0.) GO TO 10
      IF (R.GT.0.) GO TO 20
C
C           A*X = 2**22*IA1*IX1 + 2**11*(IA1*IX1 + (IA1-IA0)*(IX0-IX1)
C                   + IA0*IX0) + IA0*IX0
C
      IY0 = IA0*IX0
      IY1 = IA1*IX1 + IA1MA0*(IX0-IX1) + IY0
      IY0 = IY0 + IC
      IX0 = MOD (IY0, 2048)
      IY1 = IY1 + (IY0-IX0)/2048
      IX1 = MOD (IY1, 2048)
C
 10   RANDO = IX1*2048 + IX0
      RANDO = RANDO / 4194304.
      RETURN
C
 20   IX1 = AMOD(R,1.)*4194304. + 0.5
      IX0 = MOD (IX1, 2048)
      IX1 = (IX1-IX0)/2048
      GO TO 10
C
      END

