      PROGRAM PIC
c  A simple P.I.C. code for calculating the Langmuir sheath
c  at the edge of a plasma due to electrons only: the ions are
c  stationary.
c  NT 'electron' particles (<=1000) are used to model a one 
c  dimensional electron distribution of density DENS on a mesh 
c  of IT cells and spatial extent XT Debye lengths.
c  The plasma has ion density DENS over X1 Debye lengths, a
c  linear density fall-off over X2 Debye lengths and is zero
c  over the remainder of the mesh.
c  The electrons are assigned a Maxwell-Boltzmann velocity
c  distribution of temperature TEMP in one dimension.
c  The time step is limited by cell size and thermal speed, and
c  by the plasma frequency. For this problem the former is
c  generally dominant.
c  The positions and velocities of the electrons are output every
c  60 timesteps and there are 9 outputs including the initial
c  configuration.   These are in files V10.DAT to V18.DAT and 
c  E10.DAT to E18.DAT.   The V-files give the positions and 
c  velocities of the electrons and the E-files give the position
c  in the cell and the overall field.
c ******************************************************************
c  NOTE - THIS PROGRAM REQUIRES LONG INTEGERS
c ******************************************************************
      COMMON /BLKPOS/ X(1000)
      COMMON /BLKVEL/ V(1000)
      COMMON /BLKFLD/ E(1000)
      COMMON /BLKDEN/ RHO(1000)
      COMMON /BLKION/ RHOP(1000)
      COMMON /BLKCEL/ DX,XT,X1,X2
      COMMON /BLKDAT/ DT,EMT,ECF,DENS,TEMP,VTEMP
      COMMON /BLKTIM/ TIME
      COMMON /BLKINT/ IT,NT
      NOUT=59
      CALL INPUT
      WRITE(6,'(''OUTPUT TIMES'')')
      CALL SET
      CALL FIELD
    1 CALL PUSH
      CALL FIELD
      NOUT=NOUT+1
      IF((NOUT/60)*60.EQ.NOUT)THEN 
      CALL OUTPUT(NOUT)
      WRITE(6,*)TIME
      ENDIF
      TIME=TIME+DT
      IF (NOUT.LT.540) GO TO 1
      STOP
      END



      SUBROUTINE INPUT
      CHARACTER ANS*1
      COMMON /BLKCEL/ DX,XT,X1,X2
      COMMON /BLKDAT/ DT,EMT,ECF,DENS,TEMP,VTEMP
      COMMON /BLKTIM/ TIME
      COMMON /BLKINT/ IT,NT
      XT=10
      IT=100
      NT=1000
      DT=0.2
      X1=2.5
      X2=0
      DENS=1E25
      TEMP=1E4
   13 WRITE(6,'(''The following data are provided'')')
      WRITE(6,101)XT
  101 FORMAT(10H [1] XT = ,F8.3,30H - the extent in Debye lengths)
      WRITE(6,102)IT
  102 FORMAT(10H [2] IT = ,I4,28H - the number of cells in XT)
      WRITE(6,103)NT
  103 FORMAT(10H [3] NT = ,I5,26H - the number of electrons)
      WRITE(6,104)DT
  104 FORMAT(10H [4] DT = ,F5.2,31H - timestep control from plasma)
      WRITE(6,'(''    frequency. See comment in subroutine STEP.'')')
      WRITE(6,105)X1
  105 FORMAT(10H [5] X1 = ,F6.3,35H - the number of Debye lengths with)
      WRITE(6,'(''    density DENS.'')')
      WRITE(6,106)X2
  106 FORMAT(10H [6] X2 = ,F6.3,36H - the number of Debye lengths after)
      WRITE(6,'(''    X1 where density falls linearly to zero.'')')
      WRITE(6,107)DENS
  107 FORMAT(12H [7] DENS = ,E8.3,28H - the density of electrons.)
      WRITE(6,108)TEMP
  108 FORMAT(12H [8] TEMP = ,E8.3,28H - the electron temperature.)
      WRITE(6,'(''Do you want to change any of these? [Y/N]'')')
      READ(5,50)ANS
   50 FORMAT(A1)
      IF(ANS.EQ.'N'.OR.ANS.EQ.'n')GOTO 11
      IF(ANS.EQ.'Y'.OR.ANS.EQ.'y')GOTO 12
      GOTO 13
   12 WRITE(6,'(''Give number of item you wish to change.'')')
      READ(5,*)NITEM
      GOTO(1,2,3,4,5,6,7,8)NITEM
    1 WRITE(6,'(''Read in new value of XT.'')')
      READ(5,*)XT
      GOTO 13
    2 WRITE(6,'(''Read in new value of IT.'')')
      READ(5,*)IT
      GOTO 13
    3 WRITE(6,'(''Read in new value of NT.'')')
      READ(5,*)NT
      GOTO 13
    4 WRITE(6,'(''Read in new value of DT.'')')
      READ(5,*)DT
      GOTO 13
    5 WRITE(6,'(''Read in new value of X1.'')')
      READ(5,*)X1
      GOTO 13
    6 WRITE(6,'(''Read in new value of X2.'')')
      READ(5,*)X2
      GOTO 13
    7 WRITE(6,'(''Read in new value of DENS.'')')
      READ(5,*)DENS
      GOTO 13
    8 WRITE(6,'(''Read in new value of TEMP.'')')
      READ(5,*)TEMP
      GOTO 13
   11 RETURN
      END



      SUBROUTINE SET
c
c  A subroutine to initialise the particle position and velocity
c  and to establish the background ion density. Other calculation
c  parameters are also evaluated.
c  The particle parameters are determined using appropriate
c  random distributions.
c
      COMMON /BLKPOS/ X(1000)
      COMMON /BLKVEL/ V(1000)
      COMMON /BLKION/ RHOP(1000)
      COMMON /BLKCEL/ DX,XT,X1,X2
      COMMON /BLKDAT/ DT,EMT,ECF,DENS,TEMP,VTEMP
      COMMON /BLKTIM/ TIME
      COMMON /BLKINT/ IT,NT
      DATA IDUM/0/,IRAN/0/
c
c  WP is the plasma frequency, VTEMP the thermal velocity,
c  FACTOR the number of electrons per particle, EMT electron
c  EMT electron charge-to-mass ratio and ECF the charge-to-field
c  ratio e/epsilon0.
c
      WP=56.41457936*SQRT(DENS)
      VTEMP=SQRT(1.515623082E7*TEMP)
      DEBYE=VTEMP/WP
      XT=XT*DEBYE
      X1=X1*DEBYE
      X2=X1+X2*DEBYE
      DX=XT/FLOAT(IT)
      FACTOR=DENS*XT/FLOAT(NT)
c  This controls the timestep 
      help=0.1*DX/VTEMP
      DT=MIN(help,(DT/WP))
      EMT=-1.758804786E11*DT
      ECF=1.809527009E-8*FACTOR
c
c  Set up the electron position and velocity distribution
c  The electron density distribution is uniform up to X1 and
c  then falls linearly to X2 and is zero to XT.
c
      ALPHA=(X1+X1)/(X1+X2)
      XX1=0.5*(X1+X2)
      XX2=X2*X2-X1*X1
      DO 1 N=1,NT
      Y=RAN1(IRAN)
      IF (Y.GT.ALPHA) THEN
      X(N)=X2-SQRT(XX2*(1.0-Y))
      ELSE
      X(N)=XX1*Y
      ENDIF
      V(N)=VTEMP*GASDEV(IDUM)
    1 CONTINUE
c  Set up the background ion distribution
      RHOI=FLOAT(NT)*DX/XX1
      RHOT=0.0
      XP=-0.5*DX
      XF=0.0
      DO 2 I=1,IT
      XB=XF
      XP=XP+DX
      XF=XF+DX
      IF (XF.LT.X1) THEN
c  The point lies entirely within the uniform density zone
      RHOP(I)=RHOI
      ELSE IF (XF.LT.X2) THEN
c  The front edge of the cell lies in the falling density
      IF(XB.LT.X1)THEN
c  but the back in the uniform zone
      RHOP(I)=RHOI*((X1-XB)+(XF-X1)*(X2-0.5*(XF+X1))/(X2-X1))/DX
      ELSE
c  The cell is entirely within the falling density region
      RHOP(I)=RHOI*(X2-XP)/(X2-X1)
      ENDIF
      ELSE
c  The front edge is in the zero density region
      IF (XB.LT.X1) THEN
c  whilst the back is still in the uniform zone
      RHOP(I)=RHOI*((X1-XB)+0.5*(X2-X1))/DX
      ELSE IF (XB.LT.X2) THEN
c  or the back is in the falling zone
      RHOP(I)=RHOI*0.5*(X2-XB)*(X2-XB)/(DX*(X2-X1))
      ELSE
c  The cell is entirely in the zero density region
      RHOP(I)=0.0
      ENDIF
      ENDIF
      RHOT=RHOT+RHOP(I)
    2 CONTINUE
      IF (NINT(RHOT).NE.NT) THEN
      I=INT(X1/DX)+1
      RHOP(I)=RHOP(I)+FLOAT(NT)-RHOT
      ENDIF
      RETURN
      END



      SUBROUTINE OUTPUT(NOUT)
      CHARACTER*2 A(9)
      CHARACTER*8 FILENAM1,FILENAM2
      DIMENSION XC(1000),EE(1000)
      COMMON /BLKPOS/ X(1000)
      COMMON /BLKVEL/ V(1000)
      COMMON /BLKFLD/ E(1000)
      COMMON /BLKCEL/ DX,XT,X1,X2
      COMMON /BLKDAT/ DT,EMT,ECF,DENS,TEMP,VTEMP
      COMMON /BLKTIM/ TIME
      COMMON /BLKINT/ IT,NT
      DATA A/'10','11','12','13','14','15','16','17','18'/
      VMAX=0.0
      VMIN=0.0
      DO 1 N=1,NT
      VMAX=AMAX1(V(N),VMAX)
      VMIN=AMIN1(V(N),VMIN)
    1 CONTINUE
      VT=VMAX-VMIN
      EMAX=0.0
      EMIN=0.0
      DO 2 I=1,IT
      EMAX=AMAX1(E(I),EMAX)
      EMIN=AMIN1(E(I),EMIN)
    2 CONTINUE
      ET=EMAX-EMIN
      DO 3 I=1,NT
      X(I)=X(I)/XT
      V(I)=V(I)/VT
    3 CONTINUE
      XC(1)=0.5*DX/XT
      EE(1)=E(1)/ET
      DO 4 I=2,IT
      XC(I)=XC(I-1)+DX/XT
      EE(I)=E(I)/ET
    4 CONTINUE
      FILENAM1='V'//A(NOUT/60)//'.DAT'
      FILENAM2='E'//A(NOUT/60)//'.DAT'
      M1=10+NOUT/60
      M2=11+NOUT/60
      OPEN(UNIT=M1,FILE=FILENAM1)
      OPEN(UNIT=M2,FILE=FILENAM2)
      WRITE(M1,50)(X(I),V(I),I=1,NT)
      WRITE(M2,50)(XC(I),EE(I),I=1,IT)
   50 FORMAT(2E14.5)
      DO 5 I=1,NT
      X(I)=X(I)*XT
      V(I)=V(I)*VT
    5 CONTINUE
      RETURN
      END



      SUBROUTINE PUSH
c
c  Accelerate the particles in the electric field, and
c  move the electron position in response to its velocity
c
      COMMON /BLKPOS/ X(1000)
      COMMON /BLKVEL/ V(1000)
      COMMON /BLKFLD/ E(1000)
      COMMON /BLKDEN/ RHO(1000)
      COMMON /BLKDAT/ DT,EMT,ECF,DENS,TEMP,VTEMP
      COMMON /BLKCEL/ DX,XT,X1,X2
      COMMON /BLKINT/ IT,NT
c  If the electron moves outside the mesh it is returned to
c  maintain overall charge neutrality.
c  If it leaves the dense boundary it is returned with a
c  random velocity from a thermal distribution.
c  If it leaves on the vacuum side it is returned with velocity
c  reversed.
      DO 1 N=1,NT
      I=MIN0((1+INT(X(N)/DX)),IT)
      V(N)=V(N)+EMT*E(I)
      X(N)=X(N)+DT*V(N)
      IF (X(N).GT.XT ) THEN
      V(N)=-V(N)
      X(N)=XT+XT-X(N)
      ELSE IF (X(N).LT.0.0) THEN
      V(N)=VTEMP*ABS(GASDEV(IDUM))
      X(N)=0.0
      ENDIF
    1 CONTINUE
      RETURN
      END



      SUBROUTINE FIELD
c
c  The calculation of the field on the mesh
c
      COMMON /BLKPOS/ X(1000)
      COMMON /BLKFLD/ E(1000)
      COMMON /BLKION/ RHOP(1000)
      COMMON /BLKDEN/ RHO(1000)
      COMMON /BLKCEL/ DX,XT,X1,X2
      COMMON /BLKDAT/ DT,EMT,ECF,DENS,TEMP,VTEMP
      COMMON /BLKTIM/ TIME
      COMMON /BLKINT/ IT,NT
c  Set up the background ion charge
      DO 10 I=1,IT
      RHO(I)=RHOP(I)
   10 CONTINUE
c  Identify the cell containing the particle and assign its 
c  charge to that cell
      DO 1 N=1,NT
      I=MIN0((1+INT(X(N)/DX)),IT)
      RHO(I)=RHO(I)-1.0
    1 CONTINUE
c  Calculate the field at the cell centre with boundary condition
c  zero field at the dense boundary
      E1=0.0
      DO 2 I=1,IT
      E0=E1
      E1=E1+ECF*RHO(I)
      E(I)=0.5*(E0+E1)
    2 CONTINUE
      RETURN
      END
c
      FUNCTION RAN1(IDUM)
c     Returns a uniform deviate between 0.0 and 1.0. Set IDUM
c     to any negative value to initialise or reinintialise
c     the sequence
c
      DIMENSION R(97)
      PARAMETER (M1=259200,IA1=7141,IC1=54773,RM1=1.0/M1)
      PARAMETER (M2=134456,IA2=8121,IC2=28411,RM2=1.0/M2)
      PARAMETER (M3=243000,IA3=4561,IC3=51349)
      SAVE R,IFF,IX1,IX2,IX3
      DATA IFF /0/
c     Initialise on first call even if IDUM is not zero
      IF (IDUM.LT.0.OR.IFF.EQ.0) THEN
      IFF=1
c     Seed the first routine
      IX1=MOD((IC1-IDUM),M1)
      IX1=MOD((IA1*IX1+IC1),M1)
c     and use it to seed the second
      IX2=MOD(IX1,M2)
      IX1=MOD((IA1*IX1+IC1),M1)
c     and the third routines
      IX3=MOD(IX1,M3)
c     Fill the table with sequential uniform deviates generated
c     by the first two routines
      DO 11 J=1,97
      IX1=MOD((IA1*IX1+IC1),M1)
      IX2=MOD((IA2*IX2+IC2),M2)
c     Low and high order pieces combined here
      R(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1
   11 CONTINUE
      IDUM=1
      ENDIF
c     Except when initialising this is where we start.
c     Generate the next number for each sequence
      IX1=MOD((IA1*IX1+IC1),M1)
      IX2=MOD((IA2*IX2+IC2),M2)
      IX3=MOD((IA3*IX3+IC3),M3)
c     Use the third sequence to get an integer between 1 and 97
      J=1+(97*IX3)/M3
      IF (J.GT.97.OR.J.LT.1) WRITE (*,*) ' FAILURE IN J'
c     Return that table entry
      RAN1=R(J)
c     and refill it
      R(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1
      RETURN
      END
c
      FUNCTION GASDEV(IDUM)
c     Returns a normally distributed deviate with zero mean
c     and unit variance, using RAN1(IDUM) as the source of
c     uniform deviates
c
      SAVE ISET,GSET
      DATA ISET/0/
    1 IF (ISET.EQ.0) THEN
c     If no extra deviate is available
c     Pick two uniform numbers in the square extending from
c     -1 to +1 in each direction
      V1=2.0*RAN1(IDUM)-1.0
      V2=2.0*RAN1(IDUM)-1.0
c     Test to check they lie in unit circle
      R=V1*V1+V2*V2
      IF (R.GE.1.0) GO TO 1
c     Make the Box-Muller transformation
      FAC=SQRT(-2.0*ALOG(R)/R)
c     to get two normal deviates
      GASDEV=V1*FAC
c     We have an extra deviate available for the next call
      GSET=V2*FAC
      ISET=1
      ELSE
c     Use the spare deviate left from the previous call
      GASDEV=GSET
c     Unset the flag
      ISET=0
      ENDIF
      RETURN
      END
c
c
      BLOCK DATA
      COMMON /BLKPOS/ X(1000)
      COMMON /BLKVEL/ V(1000)
      COMMON /BLKFLD/ E(1000)
      COMMON /BLKDEN/ RHO(1000)
      COMMON /BLKTIM/ TIME
      DATA X/1000*0.0/
      DATA V/1000*0.0/
      DATA E/1000*0.0/
      DATA RHO/1000*0.0/
      DATA TIME/0.0/
      END

