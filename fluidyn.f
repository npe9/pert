      PROGRAM FLUIDYN
C----------------------------------------------------------------------
C
C  USE DOUBLE-LENGTH REAL NUMBERS AND LONG INTEGERS IN THIS PROGRAM
C
C----------------------------------------------------------------------
C  THIS IS A CELL-STRUCTURE FLUID DYNAMICS PROGRAM IN WHICH THE CELL
C  CONTAINS 125 MOLECULES.   THE MOLECULES ARE FIRST PLACED ON A 
C  5 X 5 X 5 REGULAR GRID AND THEN DISPLACED RANDOMLY BY UP TO D/20  
C  IN EACH PRINCIPAL DIRECTION WHERE D IS THE INITIAL SEPARATION.
C  AFTER 50 TIMESTEPS, TO ALLOW THE SYSTEM TO SETTLE DOWN, THE 
C  MOLECULAR CONFIGURATION IS SAMPLED EVERY TIMESTEP TO FIND THE
C  QUANTITY SUM(F(I,J)r(I,J))/3 AND ALSO THE RADIAL DENSITY FUNCTION.
C  THE PROGRAM TRANSFORMS EVERY QUANTITY INTO S.I. UNITS.  SOME FORM OF
C  REDUCED UNITS WOULD KEEP NUMBERS IN A MORE CONVENIENT RANGE BUT MAKE
C  THE PROGRAM MORE DIFFICULT TO FOLLOW. 
      DIMENSION X(125,3),V(125,3),DELX(125,3,0:4),DELV(125,3,0:4),
     +VX(3),SEPAR(125,4),W(4),TAB(0:100),DD(0:100),XT(125,3)
C  DATA FOR THE RANDOM-NUMBER GENERATOR.  IR IS THE SEED
      DATA IR,IX,IY,IM/199,171,11213,53125/
C  WEIGHTS FOR RUNGE-KUTTA
      DATA W/0,0.5,0.5,1.0/
C  DATA FOR THE LENNARD-JONES FORCE.   SIG IN ANGSTROM UNITS AND
C  EPS AS EPS/k WHERE k IS BOLTZMANN'S CONSTANT.  THE FIGURES  
C  ARE FOR ARGON.  
      DATA SIG,EPS/3.45,120/
C  CONVERT TO S.I.
      SIG=SIG*1.0E-10
      EPS=EPS*1.38E-23
C  AM IS THE MASS OF THE MOLECULE IN ATOMIC MASS UNITS AND T THE 
C  TEMPERATURE OF THE LIQUID.  THE MASS GIVEN IS FOR ARGON. 
      DATA AM,T/40,329/
C  THE RANGE OF THE FORCE, DLIM, IS SET AT 2.4 X SIG CORRESPONDING TO
C  WHERE THE FORCE FALLS TO 1% OF ITS MAXIMUM VALUE.   THE SIDE OF THE
C  CUBICAL CELL IS 5 X SIG FOR V* = 1.   THE VALUE OF V* IS READ IN
C  WHICH GIVES THE CELL EDGE, EL.
      SIG2=SIG*SIG
      TFEPS=24.0*EPS
C  MASS OF HYDROGEN ATOM.
      HM=1.667E-27
C  BOLTZMANN CONSTANT
      CAY=1.38E-23
C  NUMBER OF TIMESTEPS TO COMPLETION
      NTOT=150
      PI=4.0*ATAN(1.0)
      DLIM2=(2.4*SIG)**2
      WRITE(6,'('' READ IN THE VALUE OF V* '')')
      READ(5,*)VSTAR
      EL=5.0*SIG*VSTAR**(1.0/3.0)
      GRID=0.2*EL    
C  THE MOLECULES ARE SET ON A UNIFORM GRID AND THEN DISPACED IN THE X,
C  Y AND Z DIRECTIONS BY UP TO EL/20.0 USING A RANDOM-NUMBER GENERATOR.
C  CALCULATE THE MEAN SPEED OF THE MOLECULES = SQRT(3kT/AM)
      VMEAN=SQRT(24835*T/AM)
C  CLEAR TABLE FOR RADIAL DISTRIBUTION
      DO 17 I=1,100
      TAB(I)=0
   17 CONTINUE
      NUM=0
      DO 1 I=1,5
      DO 1 J=1,5
      DO 1 K=1,5
      NUM=NUM+1
      X(NUM,1)=(I-3)*GRID
      X(NUM,2)=(J-3)*GRID
      X(NUM,3)=(K-3)*GRID
C  NOW DISPLACE FROM GRID POSITIONS
      DO 2 L=1,3
      IR=MOD(IR*IX+IY,IM)
      X(NUM,L)=X(NUM,L)+(2*FLOAT(IR)/FLOAT(IM)-1.0)*GRID/20.0
    2 CONTINUE
C  THE SPEEDS OF THE MOLECULES SHOULD HAVE A MAXWELL-BOLTZMANN
C  DISTRIBUTION.  HOWEVER, THE MOLECULES ARE ALL GIVEN THE SAME 
C  AVERAGE SPEED, EQUAL TO SQRT(3kT/AM) BUT IN RANDOM DIRECTIONS.
C  AFTER SOME TIME THE INTERACTIONS SHOULD PRODUCE THE CORRECT 
C  DISTRIBUTION.   THE MODEL IS RUN FOR 50 TIMESTEPS BEFORE ANY
C  INFORMATION IS TAKEN FROM IT TO ALLOW THIS TO HAPPEN.
   19 DO 3 L=1,3
      IR=MOD(IR*IX+IY,IM)
      VX(L)=2*FLOAT(IR)/FLOAT(IM)-1
    3 CONTINUE
      VV=SQRT(VX(1)**2+VX(2)**2+VX(3)**2)
      IF(VV.GT.1.0)GOTO 19
      DO 4 L=1,3
      V(NUM,L)=VMEAN*VX(L)/VV
    4 CONTINUE
    1 CONTINUE
C  A VARIABLE TIMESTEP IS USED IN THIS PROGRAM. THE INITIAL H IS CHOSEN
C  SO THAT THE MAXIMUM DISTANCE TRAVELLED IS 0.05*GRID.
      H=GRID/VMEAN/20.0
      NOSTEP=0
      VIR=0           
  100 DO 5 I=1,125
      DO 5 J=1,3
      DO 5 K=0,4
      DELX(I,J,K)=0
      DELV(I,J,K)=0
    5 CONTINUE
      NOSTEP=NOSTEP+1
      IF((NOSTEP/10)*10.NE.NOSTEP)GOTO 37
      WRITE(6,250)NOSTEP
  250 FORMAT(34H NUMBER OF TIMESTEPS COMPLETED IS ,I5)    
   37 DO 60 IW=1,4
      WT=W(IW)
C  SET ALL VALUES OF DX/DT, DY/DT AND DZ/DT AND TEMPORARY VALUES OF 
C  X, Y AND Z.
      DO 16 I=1,125
      DO 16 K=1,3
      DELX(I,K,IW)=(V(I,K))+WT*H*DELV(I,K,IW-1)
      XT(I,K)=X(I,K)+WT*H*DELX(I,K,IW)
   16 CONTINUE      
      DO 18 I=1,124
C  ALL THE NEAREST NEIGHBOURS ARE NOW IDENTIFIED FOR MOLECULE I.  THIS 
C  INVOLVES FINDING THE MOLECULE J IN THE MAIN CELL OR A GHOST CELL 
C  WHICH IS NEAREST AND THEN CHECKING TO SEE IF IT IS CLOSER THAN DLIM.
      DO 18 J=I+1,125
      CALL DIS(XT,SEPAR,EL,I,J)
      IF(SEPAR(J,4).GT.DLIM2)GOTO 18
C  THE CLOSEST I-J SEPARATION HAS BEEN FOUND TO BE LESS THAN DLIM.
C  CALCULATE CONSTANTS FOR THE LENNARD-JONES FORCE COMPONENTS.
      A=(SIG2/SEPAR(J,4))**3
      B=A/SEPAR(J,4)
      C=TFEPS*B*(1.0-2.0*A)/HM/AM
C  FIND ACCELERATION COMPONENTS.
      DO 31 K=1,3
      DELV(I,K,IW)=DELV(I,K,IW)-C*SEPAR(J,K)
      DELV(J,K,IW)=DELV(J,K,IW)+C*SEPAR(J,K)
   31 CONTINUE
   18 CONTINUE
   60 CONTINUE
C  UPDATE POSITION AND VELOCITY COMPONENTS
      DO 80 I=1,125
      DO 80 K=1,3
      X(I,K)=X(I,K)+H*(DELX(I,K,1)+DELX(I,K,4)+2*(DELX(I,K,2)+
     +DELX(I,K,3)))/6.0
      V(I,K)=V(I,K)+H*(DELV(I,K,1)+DELV(I,K,4)+2*(DELV(I,K,2)+
     +DELV(I,K,3)))/6.0
   80 CONTINUE
C  RESTORE MOLECULE TO BE INSIDE THE BOX
      DO 81 I=1,125
      DO 81,K=1,3
       help=X(I,K)+9.5*EL
      X(I,K)=MOD(help,EL)-0.5*EL
   81 CONTINUE
      IF(NOSTEP.LE.50)GOTO 100
C  FIND NEW TIMESTEP
      VMAX=0
      DO 71 I=1,125
      DO 71 K=1,3
      IF(ABS(V(I,K)).GT.VMAX)VMAX=ABS(V(I,K))
   71 CONTINUE
      H=GRID/VMAX/20.0
C  CALCULATE CONTRIBUTIONS TO THE VIRIAL TERM AND TO THE RADIAL
C  DISTRIBUTION FUNCTION. 
      DO 40 I=1,124       
      DO 41 J=I+1,125
      CALL DIS(X,SEPAR,EL,I,J)
      IF(SEPAR(J,4).GT.DLIM2)GOTO 41
      A=(SIG2/SEPAR(J,4))**3
      C=TFEPS*A*(1.0-2.0*A)
C  ADD CONTRIBUTION TO VIRIAL TERM
      VIR=VIR-C/3.0
C  ADD RADIAL DISTANCE TO TABLE
      L=INT(20.0*SQRT(ABS(SEPAR(J,4)))/SIG)
      TAB(L)=TAB(L)+1
   41 CONTINUE
   40 CONTINUE      
      IF(NOSTEP.LT.NTOT)GOTO 100
C  TAKE THE AVERAGE OF THE VIRIAL TERM OVER 100 TIMESTEPS
      VIR=VIR/100.0
C  NOW CALCULATE THE FACTOR 1 + VIR/NkT            
      FACTOR=1.0+VIR/CAY/T/FLOAT(NUM)
      OPEN(UNIT=9,FILE='LPT1')
C  OUTPUT TEMPERATURE AND VSTAR
      WRITE(9,350)T,VSTAR
  350 FORMAT(15H TEMPERATURE = ,F6.1,9H VSTAR = ,F7.2)    
      WRITE(9,'('' '')')
C  OUTPUT FACTOR
      WRITE(9,200)FACTOR
  200 FORMAT(21H THE FACTOR PV/NkT = ,F7.3)  
C  OUTPUT THE RADIAL DENSITY DISTRIBUTION 4*PI*R**2*RO(R) ON AN 
C  ABSOLUTE SCALE RELATIVE TO THE AVERAGE DENSITY.    
      DO 90 I=0,47
      TAB(I)=0.96774*TAB(I)*VSTAR/PI/((I+1.0)**3-I**3)
      DD(I)=0.05*(I+0.5)
   90 CONTINUE
      WRITE(9,'('' '')')
      WRITE(9,'('' RADIAL DENSITY FUNCTION ON AN ABSOLUTE SCALE '')')
      WRITE(9,300)(DD(I),TAB(I),I=0,47)
  300 FORMAT(2(F10.3,F12.6))
      STOP
      END




      SUBROUTINE DIS(X,S,EL,I,J)
      DIMENSION X(125,3),S(125,4)
      DO 20 K=1,3
      S(J,K)=X(I,K)-X(J,K)
      IF(ABS(S(J,K)+EL).LT.ABS(S(J,K)))S(J,K)=S(J,K)+EL
      IF(ABS(S(J,K)-EL).LT.ABS(S(J,K)))S(J,K)=S(J,K)-EL
   20 CONTINUE
      S(J,4)=S(J,1)**2+S(J,2)**2+S(J,3)**2
      RETURN
      END


