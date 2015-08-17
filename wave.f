      PROGRAM WAVE
      REAL X(51),Y(51),Z(51),U(51),M,L
      CHARACTER C*1
      WRITE(6,'('' A string of mass m per unit length, of'')')
      WRITE(6,'('' length l and under tension T is fixed at'')')
      WRITE(6,'('' both ends.   It is subjected to an initial'')')
      WRITE(6,'('' displacement and then released.   The '')')
      WRITE(6,'('' programme finds subsequent displacements'')')
      WRITE(6,'('' as a function of time.'')')
      WRITE(6,'('' The program shows either on the VDU or printer'')')
      WRITE(6,'('' the displacements every IQ timesteps.If required'')')
      WRITE(6,'('' it also gives six data files, WAVE0.DAT to '')')
      WRITE(6,'('' WAVE6.DAT giving values of (x,y) starting from'')')
      WRITE(6,'('' timestep M1 and then every L1 timesteps after'')')
      WRITE(6,'('' which the program terminates.'')')
      WRITE(6,'('' '')')
      WRITE(6,'('' The standard problem is for a string with'')')
      WRITE(6,'('' m = 0.001 kg m^(-1)    l = 1.0m'')')
      WRITE(6,'('' T = 200 N              r = (c x dt / dx)^2 = 1'')')
      WRITE(6,'('' '')')
      WRITE(6,'('' The problem runs for a simulated 10 ms if not'')')
      WRITE(6,'('' terminated by giving 6 data files.'')')
c Set up standard parameters
      M=0.001
      L=1.0
      T=200.0
      R=1.0
      WRITE(6,'('' The problem may be run with other parameters'')')
      WRITE(6,'('' Use standard parameters? (Y/N)'' /)')
   20 READ(5,500)C
  500 Format(A1)
      IF(C.EQ.'Y'.OR.C.EQ.'y')GOTO 10
      IF(C.NE.'N'.AND.C.NE.'n')GOTO 20
      WRITE(6,'('' Type in values of m, l, T and r.'' /)')
      READ(5,*)M,L,T,R
   10 WRITE(6,'('' Input number of divisions of the string, n.'')')
      WRITE(6,'('' Displacements will be defined at the n-1'')')
      WRITE(6,'('' internal nodes.'' /)')
      READ(5,*)N
      NM1 = N - 1
      NP1 = N + 1
      WRITE(6,'('' '')')
c Calculate wave velocity
      CC = SQRT(T/M)
c Now calculate space and time intervals.
      DX = L/N
      DT = SQRT(R)/CC*DX
c X coordinates of node points including ends
      DO 40 I = 1,NP1
   40 X(I) = (I-1)*DX
      TIM=0
      Y(1)=0
      Y(NP1)=0
      WRITE(6,'('' Input n-1 initial dispacements at the node'')')
      WRITE(6,'('' points, y(1) to y(n-1) IN CENTIMETRES'' /)')
   32 WRITE(6,'('' If to be done by keyboard type K'')')
      WRITE(6,'('' but if by SUBROUTINE WAVIN type S'')')
      READ(5,500)C
      IF(C.EQ.'K'.OR.C.EQ.'k')GOTO 30
      IF(C.EQ.'S'.OR.C.EQ.'s')GOTO 31
      GOTO 32
   31 CALL WAVIN(N,Y)
      GOTO 33
   30 DO 60 I = 2,N
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
      WRITE(6,'(''How many timestep intervals between screen'')')
      WRITE(6,'(''or printer output?'')')
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
      WRITE(6,'(''After how many timestep do you want the'')') 
      WRITE(6,'(''first data file?'')')
      READ(5,*)M1
      WRITE(6,'(''How many timesteps between data files?'')')
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
      WRITE(10,*)X(I),Y(I)
   87 CONTINUE
      ENDIF
      Z(1) = 0
      Z(N+1) = 0
      U(1) = 0
      U(N+1) =0
      KOUNT = 1
c For the first step a different formula is used
      DO 100 I = 2,N
      Z(I) = R*(Y(I+1)+Y(I-1))/2 +(1-R)*Y(I)
  100 U(I) = Z(I)
      TIM = TIM+DT
      GOTO 79
c Now calculate the next step if time not exceeded
  200 TIM = TIM + DT
      IF(TIM.GT.0.01)GOTO 1000
      KOUNT = KOUNT + 1
      DO 140 I = 2,N
      Y(I) = Z(I)
  140 Z(I) = U(I)
      DO 160 I = 2,N
  160 U(I) = R*(Z(I+1)+Z(I-1)) + 2*(1-R)*Z(I) - Y(I)
   79 IF(IDAT.EQ.0)GOTO 94
      IF(KOUNT.LT.M1)GOTO 94
      IF(((KOUNT-M1)/L1)*L1.NE.KOUNT-M1)GOTO 94
      NN=(KOUNT-M1)/L1
      DO 95 I=1,NP1
      WRITE(10+NN,*)X(I),U(I)
   95 CONTINUE
      IF(NN.EQ.5)GOTO 1000
   94 IF((KOUNT/IQ)*IQ.EQ.KOUNT)THEN
      WRITE(IOUT,350)TIM
      WRITE(IOUT,400)(U(I),I=1,NP1)   
      ENDIF
      GOTO 200
 1000 STOP
      END


