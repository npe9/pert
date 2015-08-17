    program magelec
!  this program computes the trajectory of an electron in the x,y plane
!  subject to a magnetic field b(x,y) in the z direction.   b(x,y) is
!  defined by a statement function.   the user is asked to input the
!  initial position (x,y) and the initial velocity components (vx,vy).
!  the output is in a file 'epath.dat' but the name of the file may be
!  changed by the user.  this file may be printed or saved as input to
!  a graphics program.
    dimension x(0:100),y(0:100),vx(0:100),vy(0:100)
    character ans*1, fname*10
!  set value of e/m for electron
    eom=1.759e11
!  set the name of the output file for (x,y) values
    fname='epath.dat'
!  input initial values of x, y, vx, vy
    write(6,'('' input initial values of x and y '')')
    read(5,*)x(0),y(0)
    write(6,'('' input initial values of vx, and vy '')')
    read(5,*)vx(0),vy(0)
!  input the timestep
    write(6,'('' input the timestep'')')
    read(5,*)dt
!  start calculation
    do 1 i=1,100
    !  predict the position of the electron at the end of the timestep.
    !  this will be used to get an estimate of the average magnetic field
    !  during the timestep.
        xpred=x(i-1)+vx(i-1)*dt
        ypred=y(i-1)+vy(i-1)*dt
    !  estimate of average field in the time interval
        bav=0.5*(b(x(i-1),y(i-1))+b(xpred,ypred))
    !  calculate value of psi [equations (3.18) & (3.19)]
        psi=0.5*eom*bav*dt
    !  calculate coefficient for equations (3.18) & (3.19)
        z1=1-psi*psi
        z2=2-z1
        c1=z1/z2
        c2=psi/z2
    !  advance values of x, y, vx, and vy.
        x(i)=x(i-1)+dt*(vx(i-1)/z2+c2*vy(i-1))
        y(i)=y(i-1)+dt*(vy(i-1)/z2-c2*vx(i-1))
        vx(i)=c1*vx(i-1)+2.0*c2*vy(i-1)
        vy(i)=c1*vy(i-1)-2.0*c2*vx(i-1)
    1 end do
    4 write(6,'(''the values of x & y will be output in a file '')')
    write(6,'(''epath.dat.  do you want to change the name of'')')
    write(6,'(''the file?  [y/n]'')')
    read(5,50)ans
    50 format(a1)
    if(ans == 'n' .or. ans == 'n')goto 2
    if(ans == 'y' .or. ans == 'y')goto 3
    goto 4
    3 write(6,'(''input name of the output file [<= 10 characters]'')')
    read(5,100)fname
    100 format(a10)
    2 open(unit=11,file=fname)
    do 20 i=0,100
        write(11,*)x(i),y(i)
    20 end do
    7 write(6,'(''do you want printed output? [y/n] '')')
    read(5,50)ans
    if(ans == 'n' .or. ans == 'n')goto 5
    if(ans == 'y' .or. ans == 'y')goto 6
    goto 7
    6 open(unit=9,file='lpt1')
    write(6,200)(i,x(i),y(i),i=0,80)
    200 format(3(i6,f8.5,f8.5))
    5 stop
    end program
    function b(x,y)
    data d/0.5/
    b=d*y
    end function b
