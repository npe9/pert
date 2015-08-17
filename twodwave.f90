    program twodwave
    real :: s(51),y(51),z(51),u(51),m,ls
    character c*1
    write(6,'('' a circular drumskin of mass m {kg m^[-2]}, of'')')
    write(6,'('' radius s and with tautness t {n m^[-1]} is'')')
    write(6,'('' fixed at its boundary.it is given a circularly'')')
    write(6,'('' symmetric displacement and then released.   the '')')
    write(6,'('' programme finds subsequent circularly symmetric'')')
    write(6,'('' displacements as a function of time.'')')
    write(6,'('' the program shows either on the vdu or printer'')')
    write(6,'('' displacements across a radius each iq timesteps.'')')
    write(6,'('' if requested it gives six data files, wave0.dat'')')
    write(6,'('' wave6.dat giving values of (s,y) starting from'')')
    write(6,'('' timestep m1 and then every l1 timesteps after'')')
    write(6,'('' which the program terminates.'')')
    write(6,'('' '')')
    write(6,'('' the standard problem is for a drumskin with'')')
    write(6,'('' m = 0.01 kg m^(-2)     s = 0.25m'')')
    write(6,'('' t = 200 nm^(-1)        r = (c x dt/ds)^2 = 1'')')
    write(6,'('' '')')
    write(6,'('' the problem runs for a simulated 10 ms if not'')')
    write(6,'('' terminated by having given 6 data files.'')')
! set up standard parameters
    m=0.01
    ls=0.25
    t=200.0
    r=1.0
    write(6,'('' the problem may be run with other parameters'')')
    write(6,'('' use standard parameters? (y/n)'' /)')
    20 read(5,500)c
    500 format(a1)
    if(c == 'y' .or. c == 'y')goto 10
    if(c /= 'n' .and. c /= 'n')goto 20
    write(6,'('' type in values of m, s, t and r.'' /)')
    read(5,*)m,ls,t,r
    10 write(6,'('' input number of divisions of the radius, n.'')')
    write(6,'('' displacements will be defined at n'')')
    write(6,'('' internal nodes including the centre point.'' /)')
    read(5,*)n
    nm1 = n - 1
    np1 = n + 1
    write(6,'('' '')')
! calculate wave velocity
    cc = sqrt(t/m)
! now calculate space and time intervals.
    ds = ls/n
    dt = sqrt(r)/cc*ds
! s coordinates of node points including ends
    do 40 i = 1,np1
        s(i) = (i-1)*ds
    40 end do
    tim=0
    y(np1)=0
    write(6,'('' input n initial dispacements at the node'')')
    write(6,'('' points, y(0) [centre] to y(n-1)in cm. '' /)')
    32 write(6,'('' if input is via keyboard type k.'')')
    write(6,'('' if by subroutine wavin2.for type s.'')')
    read(5,500)c
    if(c == 'k' .or. c == 'k')goto 30
    if(c == 's' .or. c == 's')goto 31
    goto 32
    31 call wavin2(n,y)
    goto 33
    30 do 60 i = 1,n
        j = i - 1
        write(6,520)j
        520 format(' input y{',i2,'}')
        read(5,*)y(i)
    60 end do
    33 write(6,'('' input of displacements complete.'')')
! decide on output intervals for screen or printer
    dtx=1000.0*dt
    write(6,300)dtx
    300 format(17h the timestep is ,f10.7,3h ms)
    write(6,'(''how many timestep intervals required between'')')
    write(6,'(''screen or printer output?'')')
    read(5,*)iq
! decide on form of output - screen or printer.
    iout=6
    91 write(6,'(''do you want printed output? [y/n]'')')
    read(5,500)c
    if(c == 'n' .or. c == 'n')goto 92
    if(c == 'y' .or. c == 'y')then
        iout=9
        open(unit=9,file='lpt1')
        goto 92
    endif
    goto 91
    92 idat=0
    68 write(6,'(''do you want data files? [y/n]'')')
    read(5,500)c
    if(c == 'n' .or. c == 'n')goto 67
    if(c == 'y' .or. c == 'y')then
        idat=1
        open(unit=10,file='wave0.dat')
        open(unit=11,file='wave1.dat')
        open(unit=12,file='wave2.dat')
        open(unit=13,file='wave3.dat')
        open(unit=14,file='wave4.dat')
        open(unit=15,file='wave5.dat')
        write(6,'(''after how many timesteps do you want the '')')
        write(6,'(''first data file?'')')
        read(5,*)m1
        write(6,'(''how many intervals required between data files?'')')
        read(5,*)l1
        goto 67
    endif
    goto 68
    67 write(iout,350)tim
    350 format(8h time = ,f8.5)
    write(iout,400)(y(i),i=1,np1)
    400 format(8f8.2)
    if(m1 == 0 .and. idat == 1)then
        do 87 i=1,np1
            write(10,*)s(i),y(i)
        87 end do
    endif
    z(n+1) = 0
    u(n+1) =0
    kount = 1
! for the first step a different formula is used
! first deal with centre point
    z(1)=r*y(2)+(1-r)*y(1)
    u(1)=z(1)
    do 100 i = 2,n
        z(i) = r*(y(i+1)*(1+ds/2/s(i))+y(i-1)*(1-ds/2/s(i)))/2 &
        +(1-r)*y(i)
        u(i) = z(i)
    100 end do
    tim = tim+dt
    goto 79
! now calculate the next step if time not exceeded
    200 tim = tim + dt
    if(tim > 0.01)goto 1000
    kount = kount + 1
    do 140 i = 1,n
        y(i) = z(i)
        z(i) = u(i)
    140 end do
! first deal with centre point
    u(1)=2*r*z(2)+2*(1-r)*z(1)-y(1)
    do 160 i = 2,n
        u(i) = r*(z(i+1)*(1+ds/2/s(i))+z(i-1)*(1-ds/2/s(i))) &
        + 2*(1-r)*z(i) - y(i)
    160 end do
    79 if(idat == 0)goto 94
    if(kount < m1)goto 94
    if(((kount-m1)/l1)*l1 /= kount-m1)goto 94
    nn=(kount-m1)/l1
    do 95 i=1,np1
        write(10+nn,*)s(i),u(i)
    95 end do
    if(nn == 5)goto 1000
    94 if((kount/iq)*iq == kount)then
        write(iout,350)tim
        write(iout,400)(u(i),i=1,np1)
    endif
    goto 200
    1000 stop
    end program
