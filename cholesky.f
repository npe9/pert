      PROGRAM CHOLESKY
C      
C  SETS UP THE ELEMENTS OF A BANDED SYMMETRICAL MATRIX OF TYPE
C
C         5 -1  0  0 -1  0  0  0  0  0
C        -1  5 -1  0  0 -1  0  0  0  0
C         0 -1  5 -1  0  0 -1  0  0  0
C         0  0 -1  5 -1  0  0 -1  0  0
C        -1  0  0 -1  5 -1  0  0 -1  0
C         0 -1  0  0 -1  5 -1  0  0 -1  
C         0  0 -1  0  0 -1  5 -1  0  0
C         0  0  0 -1  0  0 -1  5 -1  0
C         0  0  0  0 -1  0  0 -1  5 -1
C         0  0  0  0  0 -1  0  0 -1  5
C
C  IN THE FORM REQUIRED FOR THE SUBROUTINE ICCG.  THE NUMBER OF 
C  EQUATIONS (DIMENSION OF THE MATRIX) IS N AND M DESCRIBES THE     
C  POSITION OF THE OUTER BANDS - i.e. THERE IS A NON-ZERO ELEMENT
C  AT (1,M+1).
C
      PARAMETER(N=2000,M=4)
      DIMENSION X(N),Y(N),R(N),S(N),T(N)
      OPEN(UNIT=9,FILE='LPT1')
C  SET UP MATRIX ELEMENTS
      DO 1 K=1,N
      T(K)=-5.0
      IF(K+M.LE.N)R(K+M)=-1.0
      IF(K+1.LE.N)S(K+1)=-1.0
    1 CONTINUE
      DO 2 I=M+1,N-M
      Y(I)=-1
    2 CONTINUE
      Y(1)=-3.0
      Y(N)=-3.0
      DO 3 I=2,M
      Y(I)=-2.0
      Y(N-I+1)=-2.0
    3 CONTINUE
      CALL ICCG(X,Y,R,S,T,M,N,1.0E-6,1.0E-6)
      WRITE(9,100)(X(I),I=1,N)
  100 FORMAT(8F8.4)
      STOP
      END
