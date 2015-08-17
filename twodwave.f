      PROGRAM TWODWAVE
      REAL S(51),Y(51),Z(51),U(51),M,LS
      CHARACTER C*1
      WRITE(6,'('' A circular drumskin of mass m {kg m^[-2]}, of'')')
      WRITE(6,'('' radius S and with tautness T {N m^[-1]} is'')')
      WRITE(6,'('' fixed at its boundary.It is given a circularly'')')
      WRITE(6,'('' symmetric displacement and then released.   The '')')
      WRITE(6,'('' programme finds subsequent circularly symmetric'')')
      WRITE(6,'('' displacements as a function of time.'')')
      WRITE(6,'('' The program shows either on the VDU or printer'')')
      WRITE(6,'('' displacements across a radius each IQ timesteps.'')')
      WRITE(6,'('' If requested it gives six data files, WAVE0.DAT'')')
      WRITE(6,'('' WAVE6.DAT giving values of (s,y) starting from'')')
      WRITE(6,'('' timestep M1 and then every L1 timesteps after'')')
      WRITE(6,'('' which the program terminates.'')')
      WRITE(6,'('' '')')
      WRITE(6,'('' The standard problem is for a drumskin with'')')
      WRITE(6,'('' m = 0.01 kg m^(-2)     S = 0.25m'')')
      WRITE(6,'('' T = 200 Nm^(-1)        r = (c x dt/ds)^2 = 1'')')
      WRITE(6,'('' '')')
      WRITE(6,'('' The problem runs for a simulated 10 ms if not'')')
      WRITE(6,'('' terminated by having given 6 data files.'')')
c Set up standard parameters
      M=0.01
      LS=0.25
      T=200.0
      R=1.0
      WRITE(6,'('' The problem may be run with other parameters'')')
      WRITE(6,'('' Use standard parameters? (Y/N)'' /)')
   20 READ(5,500)C
  500 Format(A1)
      IF(C.EQ.'Y'.OR.C.EQ.'y')GOTO 10
      IF(C.NE.'N'.AND.C.NE.'n')GOTO 20
      WRITE(6,'('' Type in values of m, S, T and r.'' /)')
      READ(5,*)M,LS,T,R
   10 WRITE(6,'('' Input number of divisions of the radius, n.'')')
      WRITE(6,'('' Displacements will be defined at n'')')
      WRITE(6,'('' internal nodes including the centre point.'' /)')
      READ(5,*)N
      NM1 = N - 1
      NP1 = N + 1
      WRITE(6,'('' '')')
c Calculate wave velocity
      CC = SQRT(T/M)
c Now calculate space and time intervals.
      DS = LS/N
      DT = SQRT(R)/CC*DS
c s coordinates of node points including ends
      DO 40 I = 1,NP1
   40 S(I) = (I-1)*DS
      TIM=0
      Y(NP1)=0
      WRITE(6,'('' Input n initial dispacements at the node'')')
      WRITE(6,'('' points, y(0) [centre] to y(n-1)in cm. '' /)')
   32 WRITE(6,'('' If input is via keyboard type K.'')')
      WRITE(6,'('' if by subroutine WAVIN2.for type S.'')')
      READ(5,500)C
      IF(C.EQ.'K'.OR.C.EQ.'k')GOTO 30
      IF(C.EQ.'S'.OR.C.EQ.'s')GOTO 31
      GOTO 32
   31 CALL WAVIN2(N,Y)
      GOTO 33
   30 DO 60 I = 1,N
      J = I - 1
      WRITE(6,520)J
  520 FORMAT(' Input y{',I2,'}')
      READ(5,*)Y(I)
   60 CONTINUE
   33 WRITE(6,'('' Input of displacements complete.'')')
c Decide on output intervals for screen or printer
      DTX=1000.0*DT
      WRITE(6,300)DTX
  300 FORMAT(17H The timestep is ,F10.7,3H ms)
      WRITE(6,'(''How many timestep intervals required between'')')
      WRITE(6,'(''screen or printer output?'')')
      READ(5,*)IQ
c Decide on form of output - screen or printer.
      IOUT=6
   91 WRITE(6,'(''Do you want printed output? [Y/N]'')')
      READ(5,500)C
      IF(C.EQ.'N'.OR.C.EQ.'n')GOTO 92
      IF(C.EQ.'Y'.OR.C.EQ.'y')THEN
      IOUT=9
      OPEN(UNIT=9,FILE='LPT1')
      GOTO 92
      ENDIF
      GOTO 91
   92 IDAT=0
   68 WRITE(6,'(''Do you want data files? [Y/N]'')')
      READ(5,500)C
      IF(C.EQ.'N'.OR.C.EQ.'n')GOTO 67
      IF(C.EQ.'Y'.OR.C.EQ.'y')THEN
      IDAT=1
      OPEN(UNIT=10,FILE='WAVE0.DAT')      
      OPEN(UNIT=11,FILE='WAVE1.DAT')      
      OPEN(UNIT=12,FILE='WAVE2.DAT')      
      OPEN(UNIT=13,FILE='WAVE3.DAT')      
      OPEN(UNIT=14,FILE='WAVE4.DAT')      
      OPEN(UNIT=15,FILE='WAVE5.DAT')      
      WRITE(6,'(''After how many timesteps do you want the '')')
      WRITE(6,'(''first data file?'')')
      READ(5,*)M1
      WRITE(6,'(''How many intervals required between data files?'')')
      READ(5,*)L1
      GOTO 67
      ENDIF
      GOTO 68
   67 WRITE(IOUT,350)TIM
  350 FORMAT(8H TIME = ,F8.5)
      WRITE(IOUT,400)(Y(I),I=1,NP1)
  400 FORMAT(8F8.2)
      IF(M1.EQ.0.AND.IDAT.EQ.1)THEN 
      DO 87 I=1,NP1
      WRITE(10,*)S(I),Y(I)
   87 CONTINUE
      ENDIF
      Z(N+1) = 0
      U(N+1) =0
      KOUNT = 1
c For the first step a different formula is used
c First deal with centre point
      Z(1)=R*Y(2)+(1-R)*Y(1)
      U(1)=Z(1) 
      DO 100 I = 2,N
      Z(I) = R*(Y(I+1)*(1+DS/2/S(I))+Y(I-1)*(1-DS/2/S(I)))/2 
     +  +(1-R)*Y(I)
  100 U(I) = Z(I)
      TIM = TIM+DT
      GOTO 79
c Now calculate the next step if time not exceeded
  200 TIM = TIM + DT
      IF(TIM.GT.0.01)GOTO 1000
      KOUNT = KOUNT + 1
      DO 140 I = 1,N
      Y(I) = Z(I)
  140 Z(I) = U(I)
c First deal with centre point
      U(1)=2*R*Z(2)+2*(1-R)*Z(1)-Y(1)
      DO 160 I = 2,N
  160 U(I) = R*(Z(I+1)*(1+DS/2/S(I))+Z(I-1)*(1-DS/2/S(I))) 
     + + 2*(1-R)*Z(I) - Y(I)
   79 IF(IDAT.EQ.0)GOTO 94
      IF(KOUNT.LT.M1)GOTO 94
      IF(((KOUNT-M1)/L1)*L1.NE.KOUNT-M1)GOTO 94
      NN=(KOUNT-M1)/L1
      DO 95 I=1,NP1
      WRITE(10+NN,*)S(I),U(I)
   95 CONTINUE
      IF(NN.EQ.5)GOTO 1000
   94 IF((KOUNT/IQ)*IQ.EQ.KOUNT)THEN
      WRITE(IOUT,350)TIM
      WRITE(IOUT,400)(U(I),I=1,NP1)   
      ENDIF
      GOTO 200
 1000 STOP
      END
