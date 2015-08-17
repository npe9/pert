
    program radbar
    real :: theta(0:20),thnew(0:20),x(20),len
    character q
    write(6,'('' this calculates temperature as a function of'')')
    write(6,'('' position and time for a bar of uniform cross-'')')
    write(6,'('' section embedded in a perfectly insulating '')')
    write(6,'('' material except for one face which radiates '')')
    write(6,'('' into a fixed-temperature enclosure. initially '')')
    write(6,'('' the bar is at a uniform temperature theta(init)'')')
    write(6,'('' and the external temperature is theta(ext).  the'')')
    write(6,'('' heat radiated is a*sig*(theta**4-theta(ext)**4) '')')
    write(6,'('' where a is the cross-section of the bar, sig is '')')
    write(6,'('' the stefan constant and theta is the temperature'')')
    write(6,'('' of the exposed face. '')')
    write(6,'(''  '')')
!  stefan constant
    sig = 5.67e-8
    write(6,'('' for the standard problem in si units:'')')
    write(6,'('' heat capacity, c = 386 '')')
    c = 386
    write(6,'('' thermal conductivity, cap = 401 '')')
    cap = 401
    write(6,'('' density, rho = 8920 '')')
    rho = 8920
    write(6,'('' external temperature, thex = 300 '')')
    thex = 300
    write(6,'('' initial temperature of the bar, thin = 500 '')')
    thin = 500
    write(6,'('' length of bar, len = 1.0 '')')
    write(6,'(''      '')')
    write(6,'('' these can be changed in the source program '')')
    len = 1.0
    write(6,'(''   '')')
    write(6,'('' now input number of divisions of the rod.''/)')
    read(5,*)m
    do 10 i = 0,m
        theta(i) = thin
    10 end do
    write(6,'('' input ratio=cap*dt/(c*rho*dx*dx) ''/)')
    read(5,*)ratio
    dx=len/m
    do 15 i=0,m
        x(i) = i * dx
    15 end do
    dt = ratio*c*rho*dx*dx/cap
    write(6,200)dt
    200 format(' the time interval is ',f5.1,' seconds')
    write(6,'('' input number of intervals between output'' /)')
    read(5,*)nout
    h = 4*sig*thex*thex*thex
    mm1 = m - 1
    z1 = 1 - 2*ratio
    z3 = 2*ratio*dx*sig
    z2 = thex**4
    time = 0
    110 if(time <= 1.0e-20)goto 20
    30 write(6,'('' do you wish to stop the calculation? (y/n)'' /)')
    read(5,100)q
    100 format(a1)
    if(q == 'y' .or. q == 'y')goto 500
    if(q == 'n' .or. q == 'n')goto 20
    goto 30
    20 kount = nout
    25 do 40 i=1,mm1
        thnew(i) = z1*theta(i)+ratio*(theta(i-1)+theta(i+1))
    40 end do
    thnew(0) = z1*theta(0)+2*ratio*theta(1)-z3*(theta(0)**4-z2)
    thnew(m) = z1*theta(m)+2*ratio*theta(m-1)
    do 50 i=0,m
        theta(i) = thnew(i)
    50 end do
    time=time+dt
    kount = kount - 1
    if(kount /= 0)goto 25
    no=6
    70 write(6,'('' do you want printed output? (y/n)'')')
    write(6,'('' if not then output is on the screen'' / )')
    read(5,100)q
    if(q == 'y' .or. q == 'y ')goto 60
    if(q == 'n' .or. q == 'n')goto 80
    goto 70
    60 no = 9
    open(unit=9,file='lpt1')
    80 write(no,150)time,dt,ratio
    150 format(' time=',f6.1,' with dt=',f6.1,' and ratio=',f6.3)
    write(no,151)
    151 format(1h0)
    write(no,152)
    152 format('   x     temp  ')
    write(no,153)(x(i),theta(i),i=0,m)
    153 format(f7.3,f10.3)
    goto 110
    500 if(no == 9)close(9)
    stop
    end program
