      SUBROUTINE ICCG(X,Y,R,S,T,M,N,ACC,ACC1)
C
C  SOLVES THE MATRIX EQUATION:
C  -R(K+M)*X(K+M)-S(K+1)*X(K+1)+T(K)*X(K)-S(K)*X(K-1)-R(K)*X(K-M)=Y(K)
C  USING ICCG
C  The routine is based on the ICCG solution of sparse banded matrices
c  described by Kershaw (J Comp Phys 26, 43-65, 1978)
c  The method works by forming an approximate Cholesky decomposition matrix
c  (R1,S1,T1) of the square banded symmetric matrix (R,S,T) with the same
c  sparsity as the original matrix, the inverse is readily calculated by
c  the usual triangular matrix inversion. By multiplying the input set of
c  equations by the approximate inverse a new set of equations is formed
c  with the same solution but with a matrix form which is close to the
c  identity matrix.
c  The final solution is rapidly formed by a conjugate gradient iteration
c  in a relatively small number of steps.
c  The calculation is terminated at an accuracy:
c  absolute ACC1 and relative ACC for each term.
c  Since the original banded forms are retained throughout the necessary
c  matrix multiplications are rapidly evaluated.
c  The method is easily generalised to more highly structured banded forms.
c  The final result is contained in the X array
C
      PARAMETER (ILT=200,JLT=200)
      DIMENSION X(N),Y(N),R(N),S(N),T(N)
      LOGICAL CHECK
      DIMENSION P(ILT*JLT),Q(ILT*JLT),Q1(ILT*JLT),R1(ILT*JLT),
     + S1(ILT*JLT),T1(ILT*JLT)
C
C  CALCULATE THE CHOLESKY DECOMPOSITION MATRIX (INVERSE MATRIX)
C
      K0=1-M
      K1=1
      T1(1)=T(1)
      DO 1  K=2,N
      K0=K0+1
      S1(K)=-S(K)
      S(K)=-S(K)
      IF (K0) 10,10,11
   10 T1(K)=T(K)-S1(K)**2/T1(K1)
      GO TO 12
   11 R1(K)=-R(K)
      R(K)=-R(K)
      IF (M.EQ.2) S1(K)=S1(K)-R1(K)*S1(K1)/T1(K0)
      T1(K)=T(K)-S1(K)**2/T1(K1)-R1(K)**2/T1(K0)
c  Set a limiting value for the diagonal to prevent later problems
      help=1.0e-20
      help2=abs(t(k))*.001
   12 A=MAX(help,help2)
      IF (ABS(T1(K)).LT.A) T1(K)=A
    1 K1=K
      K0=1-M
      K1=1
      DO 14 K=2,N
      K0=K0+1
      S1(K)=S1(K)/T1(K1)
      IF (K0) 14,14,13
   13 R1(K)=R1(K)/T1(K0)
   14 K1=K
C
C  PERFORM CONJUGATE GRADIENT ITERATION
C
C  X CONTAINS THE UPDATED SOLUTION VECTOR
C  Q CONTAINS THE RESIDUAL (ERROR) VECTOR
C  P CONTAINS THE INTERMEDIATE VECTOR
C
c  Prepare working arrays for entry to conjugate gradient iteration
c  Form the product of the entry solution with the input matrix
      CALL PROD(X,Q,R,S,T,M,N)
c  and generate the residual
      DO 20 K=1,N
   20 Q(K)=Y(K)-Q(K)
c  Form the approximate inverse as the initial intermediate vector
      CALL INVERT(Q,P,R1,S1,T1,M,N)
c  and the scalar product of the residual and its inverse
      B=0.0
      DO 21 K=1,N
   21 B=B+P(K)*Q(K)
c  If value small calculation complete
      IF (ABS(B).LT.1.0E-6) RETURN
c
c  ENTRY TO THE ITERATION
c
c  Number of iterations limited to number of equations
      DO 3 L=1,N
c  Form the product of the intermediate vector
c  to give the gradient direction for the residual
      CALL PROD(P,Q1,R,S,T,M,N)
c  and its scalar product with itself
      A=0.0
      DO 22 K=1,N
   22 A=A+P(K)*Q1(K)
c  to determine the gradient correction term
      A=B/A
      CHECK=.TRUE.
c  Adjust the  solution and residual vectors in the gradient direction
      DO 23 K=1,N
      X(K)=X(K)+A*P(K)
      Q(K)=Q(K)-A*Q1(K)
c  Iteration test on residual
      IF (CHECK) CHECK=ABS(Q(K)).LT.AMAX1(ACC1,ACC*ABS(Y(K)))
   23 CONTINUE
c  If accuracy level not reached continue iteration
      IF (.NOT.CHECK) GO TO 25
c  Check against original matrix equation
      CALL PROD(X,Q1,R,S,T,M,N)
      DO 24 K=1,N
      IF (CHECK) CHECK=ABS(Y(K)-Q1(K)).LT.AMAX1(ACC1,ACC*ABS(Y(K)))
   24 CONTINUE
c   If accuracy level maintained exit from routine
c  PRINCIPAL EXIT
      IF (CHECK) WRITE(6,*)L
      IF (CHECK) RETURN
c  Continue iteration by obtaining approximate inverse of residual
   25 CALL INVERT(Q,Q1,R1,S1,T1,M,N)
      A=B
c  and its scalar product with itself
      B=0.0
      DO 26 K=1,N
   26 B=B+Q(K)*Q1(K)
      IF (ABS(B).LT.1.0E-10) RETURN
c  to generate the next iterate for the intermediate array
      A=B/A
      DO 2 K=1,N
    2 P(K)=Q1(K)+A*P(K)
    3 CONTINUE
c  END OF ITERATIVE LOOP
c  Exit if iteration failed in maximum number of cycles
      RETURN
      END

      SUBROUTINE INVERT(X,Y,R,S,T,M,N)
c  Forms the approximate inverse Y of the right-hand side vector X
c  from the symmetric band triangular forms (R,S,T)
      DIMENSION X(N),Y(N),R(N),S(N),T(N)
c  The first (forward) triangular sweep
      K1=1
      K0=1-M
      Y(1)=X(1)
      DO 1 K=2,N
      K0=K0+1
      IF (K0) 10,10,11
   10 Y(K)=(X(K)-S(K)*Y(K1))
      GO TO 1
   11 Y(K)=(X(K)-S(K)*Y(K1)-R(K)*Y(K0))
    1 K1=K
      DO 2 K=1,N
    2 Y(K)=Y(K)/T(K)
c  The second (backward) triangular sweep
      K=N
      K0=N+M
      DO 3 KK=2,N
      K1=K
      K=K-1
      K0=K0-1
      IF (K0.LE.N) GO TO 30
      Y(K)=(Y(K)-S(K1)*Y(K1))
      GO TO 3
   30 Y(K)=(Y(K)-S(K1)*Y(K1)-R(K0)*Y(K0))
    3 CONTINUE
      RETURN
      END

      SUBROUTINE PROD(X,Y,R,S,T,M,N)
c Forms the product vector Y from the vector X and the symmetric band
c form (R,S,T)
      DIMENSION X(N),Y(N),R(N),S(N),T(N)
      K0=2-M
      K3=1+M
      K1=1
      K=2
      Y(1)=T(1)*X(1)+S(2)*X(2)+R(K3)*X(K3)
      DO 1 K2=3,N
      Y(K)=T(K)*X(K)+S(K)*X(K1)+S(K2)*X(K2)
      IF (K0) 11,11,10
   10 Y(K)=Y(K)+R(K)*X(K0)
   11 IF (K3.GE.N) GO TO 12
      K3=K3+1
      Y(K)=Y(K)+R(K3)*X(K3)
   12 K0=K0+1
      K1=K
    1 K=K2
      Y(N)=T(N)*X(N)+S(N)*X(K1)+R(N)*X(K0)
      RETURN
      END
