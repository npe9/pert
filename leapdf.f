      PROGRAM LEAPDF
C *******************************************************************
C  THIS PROGRAM, WHICH IS A DERIVATIVE OF HEATEX, SOLVES THE PROBLEM 
C  EITHER BY THE LEAPFROG OR THE DUFORT-FRANKEL METHOD.   THE ROUTINE
C  CYCLE HAS TWO STATEMENTS LABELLED 70.   THE FIRST OF THESE IS FOR
C  THE LEAPFROG METHOD, THE SECOND FOR DUFORT-FRANKEL.   COMMENT OUT 
C  THE ONE NOT REQUIRED.
C ******************************************************************* 
C  WHERE GRAPHICAL OUTPUT IS REQUESTED UP TO 6 DATA FILES ARE PRODUCED
C  'HEATEXn.DAT' WITH n = 1 TO 6 CONTAINING VALUES OF (X, TEMP).  IF A
C  7th OUTPUT IS REQUESTED THEN IT WILL OVERWRITE THE FIRST AND SO ON.
C  THE TIMES FOR THE OUTPUTS ARE PRINTED.
C  ******************************************************************
C  THIS USES THE LEAPFROG OR DUFORT-FRANKEL METHOD TO SOLVE THE 
C  PROBLEM OF THE TEMPERATURE VARIATION IN A LAGGED BAR WITH BOUNDARY 
C  CONDITIONS GIVING THE TEMPERATURE AT THE ENDS OF THE BAR EITHER
C  FIXED OR AS SOME FUNCTION OF TIME.   THE BOUNDARY CONDITIONS FOR  
C  THE LEFT AND RIGHT-HAND ENDS OF THE BAR ARE GIVEN BY FUNCTION 
C  SUBPROGRAMS BLH AND BRH RESPECTIVELY.  THE FOLLOWING QUANTITIES
C  ARE USED:
C           BLEN    THE LENGTH OF THE BAR
C           CAPPA   THE THERMAL CONDUCTIVITY
C           C       THE SPECIFIC HEAT CAPACITY
C           RO      THE DENSITY
C  FOR THE CALCULATION THE BAR IS DIVIDED INTO N SEGMENTS AND THE
C  THE TEMPERATURES ARE OUTPUT AT M+1 POINTS INCLUDING THE TWO ENDS.
      REAL TEMP(0:50),XX(0:50),TEMPX(2,49)
      CHARACTER ANS*1
      COMMON TEMP,NTIM,TIME,DELT,XX,IT,R,N,TMIN,TMAX,BLEN,TEMPX
C  THE FOLLOWING STATEMENT FUNCTIONS GIVE THE TEMPERATURES AT THE
C  ENDS OF THE BAR AS FUNCTIONS OF TIME
      DATA CAPPA,C,RO/200.,1000.,2700/
      OPEN(UNIT=9,FILE='LPT1')
      BLEN=1.0
      TIME=0
      NG=0
      NW=0
   20 WRITE(6,'('' THE FOLLOWING DATA [IN SI UNITS] ARE BEING USED'')')
      WRITE(6,100)BLEN
  100 FORMAT(19H 1. LENGTH OF BAR  ,F6.2)
      WRITE(6,120)CAPPA
  120 FORMAT(26H 2. THERMAL CONDUCTIVITY  ,F7.1)
      WRITE(6,140)C
  140 FORMAT(28H 3. SPECIFIC HEAT CAPACITY  ,F8.1)
      WRITE(6,160)RO
  160 FORMAT(13H 4. DENSITY  ,F8.1)
    3 WRITE(6,'('' DO YOU WANT TO CHANGE ANY OF THESE? [Y/N]'')')
      READ(5,50)ANS
   50 FORMAT(A1)
      IF(ANS.EQ.'N'.OR.ANS.EQ.'n')GOTO 1
      IF(ANS.EQ.'Y'.OR.ANS.EQ.'y')GOTO 2
      GOTO 3
    2 WRITE(6,'('' INPUT THE NUMBER OF AN ITEM TO BE CHANGED'')')
      READ(5,*)NI
      GOTO(11,12,13,14)NI
   11 WRITE(6,'('' INPUT NEW VALUE FOR LENGTH OF BAR'')')
      READ(5,*)BLEN
      GOTO 20 
   12 WRITE(6,'('' INPUT NEW THERMAL CONDUCTIVITY '')')
      READ(5,*)CAPPA
      GOTO 20
   13 WRITE(6,'('' INPUT NEW SPECIFIC HEAT CAPACITY '')')
      READ(5,*)C
      GOTO 20
   14 WRITE(6,'('' INPUT NEW DENSITY '')')
      READ(5,*)RO
      GOTO 20
    1 WRITE(6,'('' INPUT N, THE NUMBER OF SEGMENTS IN THE BAR'')') 
      WRITE(6,'(''FOR CALCULATION.'')')
      READ(5,*)N
      DELX=BLEN/N
      DO 40 I=0,N
      XX(I)=DELX*I
   40 CONTINUE
    7 WRITE(6,'('' YOU NOW HAVE A CHOICE OF FIXING EITHER THE '')')
      WRITE(6,'('' TIME INTERVAL "T" OR "R"=CAPPA*T/RO*C*DX**2'')')
      WRITE(6,'('' INPUT THE CHOICE [T/R]'')')
      READ(5,50)ANS
      IF(ANS.EQ.'T'.OR.ANS.EQ.'t')GOTO 8
      IF(ANS.EQ.'R'.OR.ANS.EQ.'r')GOTO 9
      GOTO 7
    8 WRITE(6,'('' READ IN VALUE OF "T" '')')
      READ(5,*)DELT
      R=CAPPA*DELT/DELX**2/C/RO
      WRITE(6,'('' '')')
      WRITE(6,600)R
  600 FORMAT(5H R = ,F8.4)
      GOTO 30
    9 WRITE(6,'('' READ IN THE VALUE OF R '')')
      READ(5,*)R
      DELT=R*C*RO*DELX**2/CAPPA
      WRITE(6,'('' '')')
      WRITE(6,620)DELT
  620 FORMAT(8H DELT = ,F9.4)
      WRITE(6,'(''  '')')
   30 WRITE(6,'('' READ IN N-1 INITIAL VALUES OF TEMPERATURE'')')
      WRITE(6,'('' AT INTERNAL POINTS OF THE BAR.'')')   
      WRITE(6,'(''  '')')
      DO 31 I=1,N-1
      WRITE(6,200)I
  200 FORMAT(14H READ IN TEMP[,I2,1H])    
      READ(5,*)TEMP(I)
   31 CONTINUE
      TEMP(0)=BLH(TIME)
      TEMP(N)=BRH(TIME)
C  THE INITIAL STATE OF THE BAR IS NOW FIXED
   62 WRITE(6,'('' GRAPHICAL FILE OR PRINTED OUTPUT? [G/P]'')')
      READ(5,50)ANS
      IF(ANS.EQ.'G'.OR.ANS.EQ.'g')GOTO 60
      IF(ANS.EQ.'P'.OR.ANS.EQ.'p')GOTO 61
      GOTO 62
C  THIS SECTION IS FOR PRINTED OUTPUT
   61 WRITE(6,'(''INPUT M FOR OUTPUT, WHICH MUST BE A FACTOR OF N'')')
      WRITE(6,'(''THE M+1 TEMP VALUES, INCLUDING BAR ENDS,ARE'')')
      WRITE(6,'(''EQUALLY SPACED ALONG THE BAR'')')
      READ(5,*)M
      K=N/M
      IF(TIME.LT.1.0E-15)GOTO 65
   66 WRITE(6,'(''CONTINUE THE CALCULATION? [Y/N]'')')
      READ(5,50)ANS
      IF(ANS.EQ.'Y'.OR.ANS.EQ.'y')GOTO 65
      IF(ANS.EQ.'N'.OR.ANS.EQ.'n')GOTO 500
      GOTO 66
   65 WRITE(6,'('' HOW MANY TIMESTEPS BEFORE OUTPUT?'')')
      READ(5,*)NTIM
      CALL CYCLE
      WRITE(9,300)TIME
  300 FORMAT(8H TIME = ,F7.2)
      IF(NW.NE.0)GOTO 74
      WRITE(9,320)(TEMP(I),I=0,N,K)
  320 FORMAT(6F7.1)
      GOTO 66
C  PRINTED OUTPUT SECTION COMPLETE
C  THIS SECTION IS FOR GRAPHICAL OUTPUT
   60 WRITE(6,'('' HOW MANY TIMESTEPS BEFORE OUTPUT?'')')
      READ(5,*)NTIM
      CALL CYCLE
      OPEN(UNIT=21,FILE='HEATEX1.DAT')
      OPEN(UNIT=22,FILE='HEATEX2.DAT')
      OPEN(UNIT=23,FILE='HEATEX3.DAT')
      OPEN(UNIT=24,FILE='HEATEX4.DAT')
      OPEN(UNIT=25,FILE='HEATEX5.DAT')
      OPEN(UNIT=26,FILE='HEATEX6.DAT')
   74 NG=MOD(NG,6)
      NW=NG+21
      NG=NG+1
      DO 70 I=0,N
      WRITE(NW,*)XX(I),TEMP(I)
   70 CONTINUE
      GOTO 66   
  500 STOP
      END



      REAL FUNCTION BLH(X)
      BLH=300
      RETURN
      END

      REAL FUNCTION BRH(X)
      BRH=400
      RETURN
      END 




      SUBROUTINE CYCLE
C  CONVERTS TO THE LEAPFROG METHOD WITH FIRST STATEMENT 70
C  CONVERTS TO DUFONT-FRANKEL WITH SECOND STATEMENT 70 
      REAL TEMP(0:50),TEMPX(2,49),XX(0:50)
      COMMON TEMP,NTIM,TIME,DELT,XX,IT,R,N,TMIN,TMAX,BLEN,TEMPX
      DO 61 I=1,NTIM
      TIME=TIME+DELT      
      DO 62 J=1,N-1
      IF(TIME.GT.1.1*DELT)GOTO 70
      TEMPX(2,J)=R*(TEMP(J-1)+TEMP(J+1))+(1-2*R)*TEMP(J)
      GOTO 62
C   70 TEMPX(2,J)=TEMPX(1,J)+2*R*(TEMP(J+1)+TEMP(J-1)-2*TEMP(J))
   70 TEMPX(2,J)=((1-R)*TEMPX(1,J)+R*(TEMP(J+1)+TEMP(J-1)))/(1+R)    
   62 CONTINUE
      DO 63 J=1,N-1
      TEMPX(1,J)=TEMP(J)
      TEMP(J)=TEMPX(2,J)
   63 CONTINUE
      TEMP(0)=BLH(TIME)
      TEMP(N)=BRH(TIME)
   61 CONTINUE
      RETURN
      END


