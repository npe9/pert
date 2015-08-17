      PROGRAM MCPLATE
***********************************************************************
C  THIS USES THE MONTE CARLO METHOD TO SOLVE PROBLEMS INVOLVING THE 
C  THERMAL EQUILIBRIUM IN A PLATE WITH INSULATED TOP AND BOTTOM SURFACES
C  AND BOUNDARIES WHICH ARE AT A FIXED TEMPERATURE.
C
C  THE PLATE IS DEFINED ON A SQUARE MESH CONTAINED WITHIN A 
C  RECTANGULAR REGION OF NR ROWS AND NC COLUMNS (I.E.(NR-1) X (NC-1)
C  SQUARE ELEMENTS) WITH BOTH NR AND NC <= 10.  ALL SIDES OF THE 
C  PLATE MUST EITHER BE ALONG THE PRINCIPAL DIRECTIONS OF THE MESH 
C  OR AT 45 DEGREES TO THEM.  EACH SIDE OF THE PLATE MUST CONTAIN 
C  AT LEAST THREE POINTS INCLUDING THE END POINTS.
C
C  WHEN INITIAL CONDITIONS ARE SET UP POINTS OUTSIDE THE PLATE ARE 
C  ENTERED AS X, INTERNAL POINTS, THE TEMPERATURES OF WHICH HAVE TO 
C  BE DETERMINED, ARE ENTERED AS U  AND POINTS AT THE BOUNDARY AS 
C  AS AN INTEGER <= 10,000).  
C
C  THERE IS PROVISION FOR THE PLATE TO BE HEATED (OR COOLED), WHICH
C  REQUIRES THE AMENDMENT OF THE FUNCTION SUBPROGRAM HEAT(X,Y) WHICH
C  GIVES THE RATE OF HEATING PER UNIT VOLUME AS A FUNCTION OF X AND Y.
C  IN ADDITION THE DATA STATEMENT MUST BE CHANGED TO GIVE THE LENGTH
C  OF THE SIDE OF THE SQUARE GRID, H, AND THE CONDUCTIVITY OF THE
C  PLATE, CAPPA.   THE COORDINATES (X,Y) ARE REFERRED TO THE POINT
C  (I,J) = (1,1) AS ORIGIN
C
************************************************************************
      DIMENSION IP(0:12,0:12),Q(0:12,0:12),KK(5),SD(0:12,0:12)
      REAL MCSUM
      CHARACTER EX(60)*1,STAN*13,CN*60,ANS*1
      OPEN(UNIT=9,FILE='LPT1')
      STAN='0123456789 XU'
      DATA H,CAPPA/0.125,400.0/
C DATA FOR THE RANDOM NUMBER GENERATOR.   IR IS THE SEED
      DATA IR,IX,IY,IM/199,171,11213,53125/
      SIG=5.67E-8
      PI=4.0*ATAN(1.0)
      QPF=H*H/4.0/CAPPA
C INITIALLY INDICATE ALL POINTS AS OUTSIDE PLATE AND CLEAR TABLE SD
      DO 55 I=0,12
      DO 55 J=0,12
      IP(I,J)=1
      Q(I,J)=-1
      SD(I,J)=0
   55 CONTINUE
      IOUT=6
   95 WRITE(6,'(''DO YOU WANT PRINTED OUTPUT? [Y/N]'')')
      READ(5,750)ANS
  750 FORMAT(A1)
      IF(ANS.EQ.'N'.OR.ANS.EQ.'n')GOTO 97
      IF(ANS.EQ.'Y'.OR.ANS.EQ.'y')THEN
      IOUT=9
      GOTO 97
      ENDIF
      GOTO 95   
   97 WRITE(6,'(''ARE YOU SURE THAT FUNCTION SUBPROGRAM "HEAT"'')')
      WRITE(6,'(''AND VALUES OF "H" AND "CAPPA" ARE CORRECT '')')
      WRITE(6,'(''FOR THIS "HOTPLATE" APPLICATION? IF SO THEN'')')
      WRITE(6,'(''ENTER NUMBER OF ROWS AND COLUMNS IN MESH'')')
      READ(5,*)NR,NC
C  THE CHARACTERISTICS OF THE PLATE ARE NOW ENTERED
      DO 1 I=1,NR
      WRITE(6,'('' '')')
      WRITE(6,50)I
   50 FORMAT(31H ENTER INFORMATION FOR GRID ROW,I3)
      WRITE(6,'(''X IF POINT IS OUTSIDE PLATE'')')
      WRITE(6,'(''TEMPERATURE ON THE BOUNDARY [INTEGER<=10000]'')')
      WRITE(6,'(''U IF WITHIN THE PLATE AND TO BE DETERMINED'')')
      READ(5,100)CN
  100 FORMAT(A60)
      DO 60 NX=1,60
      EX(NX)=CN(NX:NX)
   60 CONTINUE
      IT=0
      DO 2 J=1,NC
    3 IT=IT+1   
      K=INDEX(STAN,EX(IT))-1
C  TEST OF ILLEGAL CHARACTER
      IF(K.LT.0)GOTO 999
C  TEST FOR BLANK
      IF(K.EQ.10)GOTO 3
C  TEST FOR X - POINT OUTSIDE PLATE
      IF(K.EQ.11)THEN
      IP(I,J)=1
      Q(I,J)=-1
      GOTO 2
      ENDIF
C  TEST FOR U - INTERIOR POINT OF PLATE
      IF(K.EQ.12)THEN
      IP(I,J)=0
      Q(I,J)=0
      GOTO 2
      ENDIF
C  THE CHARACTER MUST BE A DIGIT
      M=0
    4 M=M+1
      KK(M)=K
C  READ NEXT DIGIT OR BLANK
      IT=IT+1
      K=INDEX(STAN,EX(IT))-1
C  TEST FOR ILLEGAL CHARACTER
      IF(K.LT.0)GOTO 999
C  TEST IF X OR U IS ILLEGALLY COMBINED WITH A DIGIT
      IF(K.GT.10)GOTO 999
C  TEST FOR BLANK WHICH IS END OF NUMBER ENTRY      
      IF(K.LE.9)GOTO 4
      IP(I,J)=3
      Q(I,J)=0         
      DO 6 L=1,M               
      Q(I,J)=Q(I,J)+KK(L)*10**(M-L)
    6 CONTINUE
    2 CONTINUE
    1 CONTINUE 
C  NMC MONTE CARLO RANDOM WALKS ARE MADE FOR EACH POINT FOR WHICH THE
C  TEMPERATURE MUST BE DETERMINED.
      WRITE(6,'('' '')')
      WRITE(6,'(''INPUT NUMBER OF MONTE CARLO TRIALS PER POINT'')')
      READ(5,*)NMC
      DO 11 I=1,NR
      DO 11 J=1,NC
      IF(IP(I,J).NE.0)GOTO 11
C  AN INTERIOR POINT HAS BEEN IDENTIFIED. RANDOM WALKS BEGIN.
      TOT=0
      TOT2=0
      DO 12 K=1,NMC
      II=I
      JJ=J
      X=(I-1)*H
      Y=(J-1)*H
      MCSUM=QPF*HEAT(X,Y)
   86 IR=MOD(IR*IX+IY,IM)
      RR=FLOAT(IR)/FLOAT(IM)
      NRR=INT(4.0*RR)+1
      GOTO(80,81,82,83)NRR
   80 II=II+1
      GOTO 84
   81 II=II-1
      GOTO 84
   82 JJ=JJ+1
      GOTO 84
   83 JJ=JJ-1
   84 IF(IP(II,JJ).EQ.3)GOTO 85
C  THE NEXT POINT IS NOT A BOUNDARY POINT
      X=(II-1)*H
      Y=(JJ-1)*H
      MCSUM=MCSUM+QPF*HEAT(X,Y)
      GOTO 86             
C  A BOUNDARY POINT HAS BEEN REACHED
   85 MCSUM=MCSUM+Q(II,JJ)
      TOT=TOT+MCSUM
      TOT2=TOT2+MCSUM*MCSUM
   12 CONTINUE
C  THE RANDOM WALKS ARE FINISHED FOR THE POINT (I, J).
      Q(I,J)=TOT/NMC
C  CALCULATE STANDARD DEVIATION
      SD(I,J)=SQRT(TOT2/NMC-Q(I,J)**2)/SQRT(FLOAT(NMC))
   11 CONTINUE
C  OUTPUT RESULTS
      WRITE(IOUT,'('' '')')
      DO 15 I=1,NR
      WRITE(IOUT,300)(Q(I,J),J=1,NC)   
  300 FORMAT(11F6.0)
  400 FORMAT(11F6.1)
      WRITE(IOUT,'('' '')')
      WRITE(IOUT,'('' '')')
   15 CONTINUE
      DO 16 I=1,NR
      WRITE(IOUT,400)(SD(I,J),J=1,NC)
      WRITE(IOUT,'('' '')')
      WRITE(IOUT,'('' '')')
   16 CONTINUE
      GOTO 998
  999 WRITE(6,'('' ILLEGAL CHARACTER - PROGRAM TERMINATED'')')
  998 STOP
      END
   
      FUNCTION HEAT(X,Y)
      HEAT=0
      END
