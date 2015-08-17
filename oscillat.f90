    program oscillate
    real :: m,k,x(0:1000),t(0:1000),wt(4), dx(0:4),du(0:4)
    character ans*1
    data wt/0.0,0.5,0.5,1.0/
!  the following data is for a standard problem.  however, the
!  user can replace any or all of these values.
    data m,f,k,ff,wc,x0,u0/1.0e-3,3.0e-3,100.0,1.0e-2,300.0,0,0/
    pi=4.0*atan(1.0)
    open(unit=9,file='lpt1')
!  open files for graphical data
    open(unit=21,file='osc1.dat')
    open(unit=22,file='osc2.dat')
    open(unit=23,file='osc3.dat')
    open(unit=24,file='osc4.dat')
    open(unit=25,file='osc5.dat')
    open(unit=26,file='osc6.dat')
    8 write(6,'(''the current values of various parameters are:'')')
    write(6,201)m
    201 format(16h 1. mass [m] is ,f8.4,3h kg)
    write(6,202)f
    202 format(28h 2. damping constant [f] is ,f8.4,5h ns/m)
    write(6,203)k
    203 format(36h 3. restoring force constant [k] is ,f8.2,4h n/m)
    write(6,204)ff
    204 format(43h 4. forcing vibration has amplitude [f] is ,f8.4,2h n)
    write(6,205)wc
    205 format(31h 5. and angular frequency [wc] ,f8.4,3h /s)
    write(6,206)x0
    206 format(32h 6. the initial position [x0] = ,f8.4,2h m)
    write(6,207)u0
    207 format(36h 7. and the initial velocity [u0] = ,f8.4,4h m/s)
!  the user may now change these values if required
    3 write(6,'(''do you want to change any of these?  y/n '')')
    read(5,100)ans
    100 format(a1)
    if(ans == 'n' .or. ans == 'n')goto 1
    if(ans == 'y' .or. ans == 'y')goto 2
    goto 3
    2 write(6,'(''indicate by number [1 to 7] the one to change'')')
    read(5,*)no
    goto(21,22,23,24,25,26,27)no
    21 write(6,'(''type in in the new value of m'')')
    read(5,*)m
    goto 8
    22 write(6,'(''type in in the new value of f'')')
    read(5,*)f
    goto 8
    23 write(6,'(''type in in the new value of k'')')
    read(5,*)k
    goto 8
    24 write(6,'(''type in in the new value of f'')')
    read(5,*)ff
    goto 8
    25 write(6,'(''type in in the new value of wc'')')
    read(5,*)wc
    goto 8
    26 write(6,'(''type in in the new value of x0'')')
    read(5,*)x0
    goto 8
    27 write(6,'(''type in in the new value of u0'')')
    read(5,*)u0
    goto 8
!  the value of alpha is now typed in
    1 write(6,'(''type in the value of alpha'')')
    read(5,*)alpha
!  print the parameters for the problem
    write(9,101)m,f,k,ff
    101 format(5h m = ,e10.5,5h f = ,e10.5,5h k = ,e10.5,5h f = ,e10.5)
    write(9,102)wc,x0,u0
    102 format(6h wc = ,e10.5,6h x0 = ,e10.5,6h u0 = ,e10.5)
    write(9,103)alpha
    103 format(9h alpha = ,f6.3)
    write(6,'(''the problem is now defined so integration can '')')
    write(6,'(''begin. the program is run in batches of 10 '')')
    write(6,'(''cycles of the forcing vibration.  after each '')')
    write(6,'(''batch the user may ask for the displacement-time'')')
    write(6,'(''values to be output to a data file oscn.dat '')')
    write(6,'(''where n runs from 1 to 6.  if a seventh output '')')
    write(6,'(''is requested then this will overwrite the first '')')
    write(6,'(''data file, an eighth request will overwrite the '')')
    write(6,'(''second data file and so on.'')')
    nout=-1
    nc1=0
!  the timestep is set at one hundreth of the period of the forcing
!  vibration.  this is more than needed to follow the forcing vibration
!  but allows a fast transient to be followed if it occurs.
    h=pi/50.0/wc
    t0=0
    64 x(0)=x0
    t(0)=t0
    do 4 j=1,1000
        t(j)=t(0)+h*j
        do 5 i=1,4
            tx=t(j-1)+wt(i)*h
            xx=x0+wt(i)*dx(i-1)
            ux=u0+wt(i)*du(i-1)
            dx(i)=h*ux
        !  guard against zero ux since 0**0 is not defined
            if(abs(ux) < 1.0e-16)then
                z=0
            else
                z=abs(ux)**(alpha-1.0)*ux
            endif
            du(i)=h*(ff*cos(wc*tx)-f*z-k*xx)/m
        5 end do
        x0=x0+(dx(1)+dx(4)+2.0*(dx(2)+dx(3)))/6.0
        u0=u0+(du(1)+du(4)+2.0*(du(2)+du(3)))/6.0
        x(j)=x0
    4 end do
    92 nc2=nc1+10
    write(6,300)nc1,nc2
    300 format(35h do you want a data file for cycles,i4,4h to ,i4,1h?)
    read(5,50)ans
    if(ans == 'n' .or. ans == 'n')goto 60
    if(ans == 'y' .or. ans == 'y')goto 91
    goto 92
    91 nout=nout+1
    nn=mod(nout,6)+21
    rewind nn
    do 93 i=0,1000
        write(nn,*)t(i),x(i)
    93 end do
    60 write(6,'('' do you wish to continue? [y/n] '')')
    read(5,50)ans
    if(ans == 'n' .or. ans == 'n')goto 999
    if(ans == 'y' .or. ans == 'y')goto 61
    goto 60
    50 format(a1)
    61 x0=x(1000)
    t0=t(1000)
    nc1=nc1+10
!  provision is made for printed numerical output.  to economize this
!  is for every fifth table entry.
    68 write(6,'('' do you wish to print the displacement and time'')')
    write(6,'('' values for the last cycle?  [y/n]'')')
    read(5,50)ans
    if(ans == 'n' .or. ans == 'n')goto 66
    if(ans == 'y' .or. ans == 'y')goto 67
    goto 68
    67 write(9,200)(t(i),x(i),i=1,1000,5)
    200 format(3(f10.4,e13.4))
    66 goto 64
    999 stop
    end program
