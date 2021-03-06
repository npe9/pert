
      PROGRAM RADBAR
      REAL THETA(0:20),THNEW(0:20),X(20),LEN
      CHARACTER Q
      WRITE(6,'('' THIS CALCULATES TEMPERATURE AS A FUNCTION OF'')')
      WRITE(6,'('' POSITION AND TIME FOR A BAR OF UNIFORM CROSS-'')')
      WRITE(6,'('' SECTION EMBEDDED IN A PERFECTLY INSULATING '')')
      WRITE(6,'('' MATERIAL EXCEPT FOR ONE FACE WHICH RADIATES '')')
      WRITE(6,'('' INTO A FIXED-TEMPERATURE ENCLOSURE. INITIALLY '')')
      WRITE(6,'('' THE BAR IS AT A UNIFORM TEMPERATURE THETA(init)'')')
      WRITE(6,'('' AND THE EXTERNAL TEMPERATURE IS THETA(ext).  THE'')')
      WRITE(6,'('' HEAT RADIATED IS A*sig*(THETA**4-THETA(ext)**4) '')')
      WRITE(6,'('' WHERE A IS THE CROSS-SECTION OF THE BAR, sig IS '')')
      WRITE(6,'('' THE STEFAN CONSTANT AND THETA IS THE TEMPERATURE'')')
      WRITE(6,'('' OF THE EXPOSED FACE. '')')
      WRITE(6,'(''  '')')
C  STEFAN CONSTANT
      SIG = 5.67E-8
      WRITE(6,'('' FOR THE STANDARD PROBLEM IN SI UNITS:'')')
      WRITE(6,'('' HEAT CAPACITY, C = 386 '')')
      C = 386
      WRITE(6,'('' THERMAL CONDUCTIVITY, CAP = 401 '')')
      CAP = 401
      WRITE(6,'('' DENSITY, RHO = 8920 '')')
      RHO = 8920
      WRITE(6,'('' EXTERNAL TEMPERATURE, THEX = 300 '')')
      THEX = 300
      WRITE(6,'('' INITIAL TEMPERATURE OF THE BAR, THIN = 500 '')')
      THIN = 500
      WRITE(6,'('' LENGTH OF BAR, LEN = 1.0 '')')
      WRITE(6,'(''      '')')
      WRITE(6,'('' THESE CAN BE CHANGED IN THE SOURCE PROGRAM '')')
      LEN = 1.0
      WRITE(6,'(''   '')')
      WRITE(6,'('' NOW INPUT NUMBER OF DIVISIONS OF THE ROD.''/)')
      READ(5,*)M
      DO 10 I = 0,M
   10 THETA(I) = THIN
      WRITE(6,'('' INPUT RATIO=CAP*DT/(C*RHO*DX*DX) ''/)')
      READ(5,*)RATIO
      DX=LEN/M
      DO 15 I=0,M
   15 X(I) = I * DX
      DT = RATIO*C*RHO*DX*DX/CAP
      WRITE(6,200)DT
  200 FORMAT(' THE TIME INTERVAL IS ',F5.1,' SECONDS')
      WRITE(6,'('' INPUT NUMBER OF INTERVALS BETWEEN OUTPUT'' /)')
      READ(5,*)NOUT
      H = 4*sig*THEX*THEX*THEX
      MM1 = M - 1
      Z1 = 1 - 2*RATIO
      Z3 = 2*RATIO*DX*SIG
      Z2 = THEX**4
      TIME = 0
  110 IF(TIME.LE.1.0E-20)GOTO 20
   30 WRITE(6,'('' DO YOU WISH TO STOP THE CALCULATION? (Y/N)'' /)')
      READ(5,100)Q
  100 FORMAT(A1)
      IF(Q.EQ.'Y'.OR.Q.EQ.'y')GOTO 500
      IF(Q.EQ.'N'.OR.Q.EQ.'n')GOTO 20
      GOTO 30
   20 KOUNT = NOUT
   25 DO 40 I=1,MM1
   40 THNEW(I) = Z1*THETA(I)+RATIO*(THETA(I-1)+THETA(I+1))
      THNEW(0) = Z1*THETA(0)+2*RATIO*THETA(1)-Z3*(THETA(0)**4-Z2)
      THNEW(M) = Z1*THETA(M)+2*RATIO*THETA(M-1)
      DO 50 I=0,M
   50 THETA(I) = THNEW(I)
      TIME=TIME+DT
      KOUNT = KOUNT - 1
      IF(KOUNT.NE.0)GOTO 25
      NO=6
   70 WRITE(6,'('' DO YOU WANT PRINTED OUTPUT? (Y/N)'')')
      WRITE(6,'('' IF NOT THEN OUTPUT IS ON THE SCREEN'' / )')
      READ(5,100)Q
      IF(Q.EQ.'Y'.OR.Q.EQ.'y ')GOTO 60
      IF(Q.EQ.'N'.OR.Q.EQ.'n')GOTO 80
      GOTO 70
   60 NO = 9
      OPEN(UNIT=9,FILE='LPT1')
   80 WRITE(NO,150)TIME,DT,RATIO
  150 FORMAT(' TIME=',F6.1,' with DT=',F6.1,' and RATIO=',F6.3)
      WRITE(NO,151)
  151 FORMAT(1H0)
      WRITE(NO,152)
  152 FORMAT('   X     TEMP  ')
      WRITE(NO,153)(X(I),THETA(I),I=0,M)
  153 FORMAT(F7.3,F10.3)
      GOTO 110
  500 IF(NO.EQ.9)CLOSE(9)
      STOP
      END
