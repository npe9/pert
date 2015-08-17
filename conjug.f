      PROGRAM CONJUG
C  A PROGRAM FOR THE SOLUTION OF LINEAR EQUATIONS WITH A SPARSE 5-BANDED
C  COEFFICIENT MATRIX OF THE FORM (FOR DIMENSION = 10):
C
C            5 -1  0  0 -1  0  0  0  0  0     
C           -1  5 -1  0  0 -1  0  0  0  0 
C            0 -1  5 -1  0  0 -1  0  0  0
C            0  0 -1  5 -1  0  0 -1  0  0
C           -1  0  0 -1  5 -1  0  0 -1  0
C            0 -1  0  0 -1  5 -1  0  0 -1
C            0  0 -1  0  0 -1  5 -1  0  0
C            0  0  0 -1  0  0 -1  5 -1  0
C            0  0  0  0 -1  0  0 -1  5 -1
C            0  0  0  0  0 -1  0  0 -1  5
C
C  AS PROVIDED THE DIMENSION OF THE MATRIX CAN BE UP TO 2000
C      
      DIMENSION NL(10000,2),AL(10000),X(2000),B(2000),Y(2000)
C
C      N IS THE DIMENSION OF THE MATRIX
C      L GIVES THE POSITION OF THE OUTER BAND - i.e. ELEMENT IN (1,L+1)
C      IO GIVES OUTPUT MODE - 6 FOR SCREEN, 9 FOR PRINTER
C      TEST GIVES THE PRECISION REQUIRED IN THE SOLUTION
C      
      DATA N,L,IO,TEST/2000,4,6,0.000001/
      OPEN(UNIT=9,FILE='LPT1')
      DO 11 I=1,N
C  SET DIAGONAL ELEMENTS      
      NL(I,1)=I
      NL(I,2)=I
      AL(I)=5.0
      IF(I.EQ.N)GOTO 11
C  SET ELEMENTS NEIGHBOURING DIAGONAL
      IT=N+2*I-1
      NL(IT,1)=I
      NL(IT,2)=I+1
      NL(IT+1,1)=I+1
      NL(IT+1,2)=I
      AL(IT)=-1.0
      AL(IT+1)=-1.0
      IF(I.GT.N-L)GOTO 11
C  SET OUTER-BAND ELEMENTS
      IT=3*N-3+2*I
      NL(IT,1)=I
      NL(IT,2)=I+L
      NL(IT+1,1)=I+L
      NL(IT+1,2)=I
      AL(IT)=-1.0
      AL(IT+1)=-1.0
   11 CONTINUE
      LIST=5*N-2*L-2
C  SET RIGHT-HAND-SIDE VECTOR
      B(1)=3.0
      B(N)=3.0
      DO 22 I=L+1,N-L
      B(I)=1
   22 CONTINUE
      DO 23 I=2,L
      B(I)=2.0
      B(N-I+1)=2.0
   23 CONTINUE
      DO 60 I=1,N
      X(I)=0
      Y(I)=0
   60 CONTINUE
C  OBTAIN SOLUTION.  MAXIMUM OF 5000 ITERATIONS ALLOWED 
      DO 3 I=1,5000
      CALL CONJUGAT(AL,NL,X,B,N,LIST)
      GREAT=0.0
      DO 50 J=1,N
      DIF=ABS(X(J)-Y(J))
      IF(DIF.GT.GREAT)GREAT=DIF
   50 CONTINUE
      IF(GREAT.LT.TEST)GOTO 30
      GOTO 43
   30 CONTINUE
      WRITE(IO,200)I
      WRITE(IO,100)(X(II),II=1,N)
      GOTO 31
  100 FORMAT(8F8.4)
  200 FORMAT(23H NUMBER OF ITERATIONS =,I6)
   43 DO 53 J=1,N
      Y(J)=X(J)
   53 CONTINUE
    3 CONTINUE
      WRITE(6,'('' TOO MANY ITERATIONS'')')
   31 STOP
      END
      

      SUBROUTINE CONJUGAT(AL,NL,X,B,N,LIST)
C
C  A GENERAL CONJUGATE GRADIENT SOLVER FOR SPARSE-MATRIX LINEAR 
C  EQUATIONS.  THE ELEMENT (NL(K,1),NL(K,2)) HAS VALUE A(K) WHERE K
C  GOES FROM 1 TO NEQ AND NEQ IS THE TOTAL NUMBER OF NON-ZERO ELEMENTS.
C  THE SUBROUTINE CAN HANDLE MATRICES UP TO DIMENSION 'NDIM' AND UP TO
C  'NELEM' NON-ZERO ELEMENTS.
C
      PARAMETER (NDIM=2000,NELEM=10000)
      DIMENSION AL(NELEM),NL(NELEM,2),X(NDIM),B(NDIM),U(NDIM),C(NDIM),
     +D(NDIM),E(NDIM)
      DO 1 I=1,N
      U(I)=0
      C(I)=0
      D(I)=0
      E(I)=0
    1 CONTINUE
      DO 2 I=1,LIST
      N1=NL(I,1)
      N2=NL(I,2)
      C(N1)=C(N1)+AL(I)*X(N2)
    2 CONTINUE
      DO 3 I=1,N
    3 D(I)=C(I)-B(I)
      DO 4 I=1,LIST
      N1=NL(I,1)
      N2=NL(I,2)
      U(N2)=U(N2)+AL(I)*D(N1)
    4 CONTINUE
      DO 5 I=1,LIST      
      N1=NL(I,1)
      N2=NL(I,2)
      E(N1)=E(N1)+AL(I)*U(N2)
    5 CONTINUE
      SUM1=0
      SUM2=0
      DO 6 I=1,N
      SUM1=SUM1+E(I)*D(I)
      SUM2=SUM2+E(I)*E(I)
    6 CONTINUE
      ELAM=SUM1/SUM2
      DO 7 I=1,N
      X(I)=X(I)-ELAM*U(I)
    7 CONTINUE
      RETURN
      END

