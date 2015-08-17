    program heatcrni
! *******************************************************************
!  this program differs from heatex which can solve the same problem
!  using the explicit method only in the subroutine 'cycle.
! *******************************************************************
!  where graphical output is requested up to 6 data files are produced
!  'heatexn.dat' with n = 1 to 6 containing values of (x, temp).  if a
!  7th output is requested then it will overwrite the first and so on.
!  the times for the outputs are printed.
!  ******************************************************************
!  this uses the crank-nicholson method to solve the problem of the
!  temperature variation in a lagged bar with boundary conditions
!  giving the temperature at the ends of the bar either fixed or
!  as some function of time.   the boundary conditions for the
!  left and right-hand ends of the bar are given by function
!  subprograms blh and brh respectively.  the following quantities
!  are used:
!           blen    the length of the bar
!           cappa   the thermal conductivity
!           c       the specific heat capacity
!           ro      the density
!  for the calculation the bar is divided into n segments and the
!  the temperatures are output at m+1 points including the two ends.
    real :: temp(0:50),xx(0:50)
    character ans*1
    common temp,ntim,time,delt,xx,it,r,n,tmin,tmax,blen
!  the following statement functions give the temperatures at the
!  ends of the bar as functions of time
    data cappa,c,ro/200.,1000.,2700/
    open(unit=9,file='lpt1')
    blen=1.0
    time=0
    ng=0
    nw=0
    20 write(6,'('' the following data [in si units] are being used'')')
    write(6,100)blen
    100 format(19h 1. length of bar  ,f6.2)
    write(6,120)cappa
    120 format(26h 2. thermal conductivity  ,f7.1)
    write(6,140)c
    140 format(28h 3. specific heat capacity  ,f8.1)
    write(6,160)ro
    160 format(13h 4. density  ,f8.1)
    3 write(6,'('' do you want to change any of these? [y/n]'')')
    read(5,50)ans
    50 format(a1)
    if(ans == 'n' .or. ans == 'n')goto 1
    if(ans == 'y' .or. ans == 'y')goto 2
    goto 3
    2 write(6,'('' input the number of an item to be changed'')')
    read(5,*)ni
    goto(11,12,13,14)ni
    11 write(6,'('' input new value for length of bar'')')
    read(5,*)blen
    goto 20
    12 write(6,'('' input new thermal conductivity '')')
    read(5,*)cappa
    goto 20
    13 write(6,'('' input new specific heat capacity '')')
    read(5,*)c
    goto 20
    14 write(6,'('' input new density '')')
    read(5,*)ro
    goto 20
    1 write(6,'('' input n, the number of segments in the bar'')')
    write(6,'(''for calculation.'')')
    read(5,*)n
    delx=blen/n
    do 40 i=0,n
        xx(i)=delx*i
    40 end do
    7 write(6,'('' you now have a choice of fixing either the '')')
    write(6,'('' time interval "t" or "r"=cappa*t/ro*c*dx**2'')')
    write(6,'('' input the choice [t/r]'')')
    read(5,50)ans
    if(ans == 't' .or. ans == 't')goto 8
    if(ans == 'r' .or. ans == 'r')goto 9
    goto 7
    8 write(6,'('' read in value of "t" '')')
    read(5,*)delt
    r=cappa*delt/delx**2/c/ro
    write(6,'('' '')')
    write(6,600)r
    600 format(5h r = ,f8.4)
    goto 30
    9 write(6,'('' read in the value of r '')')
    read(5,*)r
    delt=r*c*ro*delx**2/cappa
    write(6,'('' '')')
    write(6,620)delt
    620 format(8h delt = ,f9.4)
    write(6,'(''  '')')
    30 write(6,'('' read in n-1 initial values of temperature'')')
    write(6,'('' at internal points of the bar.'')')
    write(6,'(''  '')')
    do 31 i=1,n-1
        write(6,200)i
        200 format(14h read in temp[,i2,1h])
        read(5,*)temp(i)
    31 end do
    temp(0)=blh(time)
    temp(n)=brh(time)
!  the initial state of the bar is now fixed
    62 write(6,'('' graphical file or printed output? [g/p]'')')
    read(5,50)ans
    if(ans == 'g' .or. ans == 'g')goto 60
    if(ans == 'p' .or. ans == 'p')goto 61
    goto 62
!  this section is for printed output
    61 write(6,'(''input m for output, which must be a factor of n'')')
    write(6,'(''the m+1 temp values, including bar ends,are'')')
    write(6,'(''equally spaced along the bar'')')
    read(5,*)m
    k=n/m
    if(time < 1.0e-15)goto 65
    66 write(6,'(''continue the calculation? [y/n]'')')
    read(5,50)ans
    if(ans == 'y' .or. ans == 'y')goto 65
    if(ans == 'n' .or. ans == 'n')goto 500
    goto 66
    65 write(6,'('' how many timesteps before output?'')')
    read(5,*)ntim
    call cycle
    write(9,300)time
    300 format(8h time = ,f7.2)
    if(nw /= 0)goto 74
    write(9,320)(temp(i),i=0,n,k)
    320 format(6f7.1)
    goto 66
!  printed output section complete
!  this section is for graphical output
    60 write(6,'('' how many timesteps before output?'')')
    read(5,*)ntim
    call cycle
    open(unit=21,file='heatex1.dat')
    open(unit=22,file='heatex2.dat')
    open(unit=23,file='heatex3.dat')
    open(unit=24,file='heatex4.dat')
    open(unit=25,file='heatex5.dat')
    open(unit=26,file='heatex6.dat')
    74 ng=mod(ng,6)
    nw=ng+21
    ng=ng+1
    do 70 i=0,n
        write(nw,*)xx(i),temp(i)
    70 end do
    goto 66
    500 stop
    end program

    real function blh(x)
    blh=300
    return
    end function blh

    real function brh(x)
    brh=400
    return
    end function brh





    subroutine cycle
    real :: temp(0:50),xx(0:50),a(3,50),c(50),x(50)
    common temp,ntim,time,delt,xx,it,r,n,tmin,tmax,blen
    do 61 j=1,ntim
        time=time+delt
        do 62 i=1,n-1
            c(i)=r*temp(i-1)+2*(1-r)*temp(i)+r*temp(i+1)
            if(i == 1)goto 70
            a(1,i)=-r
            goto 80
            70 c(1)=c(1)+r*blh(time)
            80 a(2,i)=2*(1+r)
            if(i == n-1)goto 71
            a(3,i)=-r
            goto 62
            71 c(n-1)=c(n-1)+r*brh(time)
        62 end do
        n1=n-1
        call tridiag(a,c,x,n1)
        do 74 i=1,n-1
            temp(i)=x(i)
        74 end do
        temp(0)=blh(time)
        temp(n)=brh(time)
    61 end do
    return
    end subroutine cycle
