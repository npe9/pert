      PROGRAM GAUSSOR
C  THIS SOLVES A SET OF LINEAR EQUATIONS, Ax = b, WHEN THE MATRIX HAS
C  A STRONG DIAGONAL.   THE SOR (SUCCESSIVE OVER-RELAXATION) METHOD IS
C  USED, WHICH IS THE SAME AS GAUSS-SEIDEL IF THE OVER-RELAXATION
C  FACTOR IS 1.  UP TO 25 EQUATIONS CAN BE HANDLED WITH ARRAY DIMENSIONS
C  PROVIDED.  THE EQUATIONS SHOULD BE READ IN THE ORDER WHICH GIVES
C  THE STRONG DIAGONAL
      REAL A(25,25),B(25),X(25)
      CHARACTER ANS*1
      OPEN(UNIT=9,FILE='LPT1')
      WRITE(6,'(''READ IN NUMBER OF EQUATIONS [<=25] '')')
      READ(5,*)N
      WRITE(6,'(''READ IN TOLERANCE '')')
      READ(5,*)TOL
      WRITE(6,50)N
   50 FORMAT(23H READ IN COEFFICIENTS -,I3,10H PER LINE.) 
      DO 1 I=1,N
      WRITE(6,100)I
  100 FORMAT(13H READ IN ROW ,I2)
      READ(5,*)(A(I,J),J=1,N)
    1 CONTINUE
      WRITE(6,'(''READ IN ELEMENTS OF RHS VECTOR'')')
      READ(5,*)(B(J),J=1,N)
      WRITE(6,'(''READ IN FIRST ESTIMATE OF SOLUTION IN THE FORM'')')
      WRITE(6,'(''OF THE VECTOR ELEMENTS.  MAKING ELEMENTS EQUAL '')')
      WRITE(6,'(''ZERO IS USUALLY SUITABLE.'')')
      READ(5,*)(X(J),J=1,N)
      WRITE(6,'(''READ IN OVER-RELAXATION FACTOR [1.0 TO 2.0]'')')
C  ALL INFORMATION HAS NOW BEEN ENTERED.  THE SOLUTION PROCESS BEGINS.
      READ(5,*)W
   10 WRITE(6,'(''DO YOU WANT PRINTED OUTPUT? [Y/N]'')')
      READ(5,250)ANS
  250 FORMAT(A1)
      IF(ANS.EQ.'Y'.OR.ANS.EQ.'y')GOTO 8
      IF(ANS.EQ.'N'.OR.ANS.EQ.'n')GOTO 9
      GOTO 10
    8 M=9
      GOTO 12
    9 M=6 
   12 ICYCLE=0
    5 ICYCLE=ICYCLE+1
      DIF=0
      DO 7 I=1,N
      SUM=B(I)
      DO 3 J=1,N
      IF(J.EQ.I)GOTO 3
      SUM=SUM-A(I,J)*X(J)
    3 CONTINUE
      EST=SUM/A(I,I)
      Z=W*(EST-X(I))
      X(I)=X(I)+Z
      IF(ABS(Z).GT.DIF)DIF=ABS(Z)
    7 CONTINUE
      IF(DIF.LT.TOL)GOTO 4
C  NOT MORE THAN 100 CYCLES ALLOWED
      IF(ICYCLE.LT.100)GOTO 5
    4 WRITE(M,150)ICYCLE
  150 FORMAT(16H SOLUTION AFTER ,I3,8H CYCLES.)
      WRITE(M,200)(I,X(I),I=1,N)
  200 FORMAT(4(I7,F8.3))
      STOP
      END              