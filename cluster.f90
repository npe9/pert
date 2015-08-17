    program cluster
!  this program sets up a cluster of n solar-mass stars in a spherical
!  volume of radius r.   the positions of the stars are chosen randomly
!  and then adjusted to give the centre of mass at the origin.
!  the velocities of the stars are chosed to satisfy the virial theorem
!  and also to give zero momentum. stars are given close to the correct
!  rms speed but in a random direction.
!  the initial step is taken by a predictor-corrector process and then
!  the leapfrog method is started with a timestep equal to
!  rmin/(100*vmax) where rmin is the minimum distance between stars
!  and vmax the maximum speed of any star.  in each cycle rmin and
!  vmax are found.   if the ratio rmin/(vmax*dt) is less than 50 or more
!  than 150 then it is readjusted to rmin/(100*vmax) and a predictor-
!  corrector step taken to restart the leapfrog method.
    dimension x(100),y(100),z(100),vx(100),vy(100),vz(100)
    dimension delvx(100,2),delvy(100,2),delvz(100,2),stat(4)
    character ans*1
    pi=4.0*atan(1.0)
    g=4.0*pi*pi
    write(6,'(''input the number of stars in the cluster [<=100]'')')
    read(5,*)n
    write(6,'(''input the radius of the cluster in a.u.'')')
    read(5,*)r
    write(6,'(''input total time for the simulation in years'')')
    read(5,*)totime
    55 write(6,'(''print initial and final coordinates and '')')
    write(6,'(''velocities [y/n]'')')
    read(5,50)ans
    50 format(a1)
    if(ans == 'n' .or. ans == 'n')then
        key=1
        goto 60
    endif
    if(ans == 'y' .or. ans == 'y')then
        key=0
        goto 60
    endif
    goto 55
    60 call initial(x,y,z,vx,vy,vz,n,r,key)
!     store the initial positions and velocities in file start.dat
    open(unit=10, file='start.dat')
    write(10,*)n,(x(i),y(i),z(i),vx(i),vy(i),vz(i),i=1,n)
!  a choice is made whether or not to determine the geometric
!  moment of inertia, the position and speed of the centre of mass
!  and the energy.
    38 write(6,'(''do you wish to print the initial geometric'')')
    write(6,'(''moment of inertia, position and speed of '')')
    write(6,'(''centre of mass and energy? [y/n]'')')
    read(5,50)ans
    if(ans == 'n' .or. ans == 'n')goto 39
    if(ans == 'y' .or. ans == 'y')then
        call stats(x,y,z,vx,vy,vz,n,g,stat)
        goto 40
    endif
    goto 38
    40 write(9,'('' '')')
    write(9,'('' '')')
    write(9,'(''the initial statistics'')')
    write(9,300)stat
    300 format(6h r2 = ,e10.3,10h compos = ,e10.3,10h comvel = , &
    e10.3,10h toteng = ,f10.4)
!  find minimum distance and greatest speed
    39 rmin2=1.0e12
    vmax2=0
    do 1 i=1,n-1
        vv=vx(i)**2+vy(i)**2+vz(i)**2
        if(vv > vmax2)vmax2=vv
        do 1 j=i+1,n
            rr=(x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2
            if(rr < rmin2)rmin2=rr
    1 end do
    vv=vx(n)**2+vy(n)**2+vz(n)**2
    if(vv > vmax2)vmax2=vv
    time=0
!  now fix either initial timestep or revised timestep
    10 dt=0.01*sqrt(rmin2/vmax2)
!  if the timestep is less than one year the program will run too
!  slowly.  in this case terminate the program with message.
    if(dt < 1)then
        open(unit=9,file='lpt1')
        write(9,200)time
        write(9,'(''try again with a new seed in subroutine initial'')')
        200 format(22h program terminated at,f8.1,7h years.)
        goto 20
    endif
!  now take an initial predictor-corrector step.   this needs to give
!  only the coordinates at time dt.
!  first clear tables in which initial accelerations and estimates of
!  final accelerations will be stored.
    do 2 i=1,n
        do 2 j=1,2
            delvx(i,j)=0
            delvy(i,j)=0
            delvz(i,j)=0
    2 end do
!  calculate accelerations at initial and estimated final positions
    do 3 k=1,2
        do 3 i=1,n-1
            do 3 j=i+1,n
                xd=x(i)-x(j)+(k-1)*dt*(vx(i)-vx(j))
                yd=y(i)-y(j)+(k-1)*dt*(vy(i)-vy(j))
                zd=z(i)-z(j)+(k-1)*dt*(vz(i)-vz(j))
                ax=acc(g,xd,yd,zd)
                delvx(i,k)=delvx(i,k)+ax
                delvx(j,k)=delvx(j,k)-ax
                ax=acc(g,yd,zd,xd)
                delvy(i,k)=delvy(i,k)+ax
                delvy(j,k)=delvy(j,k)-ax
                ax=acc(g,zd,xd,yd)
                delvz(i,k)=delvz(i,k)+ax
                delvz(j,k)=delvz(j,k)-ax
    3 end do
!  now estimate the new positions from average velocities in the timestep
    do 4 i=1,n
        x(i)=x(i)+dt*(vx(i)+0.25*dt*(delvx(i,1)+delvx(i,2)))
        y(i)=y(i)+dt*(vy(i)+0.25*dt*(delvy(i,1)+delvy(i,2)))
        z(i)=z(i)+dt*(vz(i)+0.25*dt*(delvz(i,1)+delvz(i,2)))
    4 end do
!  the leapfrog process can now be started
!  compute new values of velocity
    time=time+dt
    11 do 5 i=1,n
        delvx(i,1)=0
        delvy(i,1)=0
        delvz(i,1)=0
    5 end do
    rmin2=1.0e12
    do 6 i=1,n-1
        do 6 j=i+1,n
            xd=x(i)-x(j)
            yd=y(i)-y(j)
            zd=z(i)-z(j)
            rd2=xd*xd+yd*yd+zd*zd
            if(rd2 < rmin2)rmin2=rd2
            ax=acc(g,xd,yd,zd)
            delvx(i,1)=delvx(i,1)+ax
            delvx(j,1)=delvx(j,1)-ax
            ax=acc(g,yd,zd,xd)
            delvy(i,1)=delvy(i,1)+ax
            delvy(j,1)=delvy(j,1)-ax
            ax=acc(g,zd,xd,yd)
            delvz(i,1)=delvz(i,1)+ax
            delvz(j,1)=delvz(j,1)-ax
    6 end do
    do 7 i=1,n
        vx(i)=vx(i)+2*dt*delvx(i,1)
        vy(i)=vy(i)+2*dt*delvy(i,1)
        vz(i)=vz(i)+2*dt*delvz(i,1)
    7 end do
! compute new values of position
    vmax2=0
    do 8 i=1,n
        x(i)=x(i)+2*dt*vx(i)
        y(i)=y(i)+2*dt*vy(i)
        z(i)=z(i)+2*dt*vz(i)
        vv=vx(i)**2+vy(i)**2+vz(i)**2
        if(vv > vmax2)vmax2=vv
    8 end do
    time=time+2*dt
    if(time > totime)goto 20
!  test for timestep suitability with values of rmin2 and vmax2
!  while they have come from different timesteps it still provides
!  a suitable test.
    test=rmin2/vmax2/dt/dt
    if(test < 2500 .or. test > 12500)goto 10
    goto 11
!  the calculation is complete with total time >= totime.   the final
!  positions and coordinates are stored in file finish.dat.   first
!  subtract one timestep of velocity from the positions to give
!  positions and velocities at the same time.
    20 do 12 i=1,n
        x(i)=x(i)-dt*vx(i)
        y(i)=y(i)-dt*vy(i)
        z(i)=z(i)-dt*vz(i)
    12 end do
    open(unit=11,file='finish.dat')
    write(11,*)n,(x(i),y(i),z(i),vx(i),vy(i),vz(i),i=1,n)
    if(key == 0)write(9,100)(x(i),y(i),z(i),vx(i),vy(i),vz(i),i=1,n)
    100 format(6e11.4)
!  a choice is made whether or not to determine the geometric
!  moment of inertia, the position and speed of the centre of mass
!  and the energy.
    48 write(6,'(''do you wish to print the final geometric'')')
    write(6,'(''moment of inertia, position and speed of '')')
    write(6,'(''centre of mass and energy? [y/n]'')')
    read(5,50)ans
    if(ans == 'n' .or. ans == 'n')goto 49
    if(ans == 'y' .or. ans == 'y')then
        call stats(x,y,z,vx,vy,vz,n,g,stat)
        goto 70
    endif
    goto 48
    70 write(9,'(''  '')')
    write(9,400)time
    400 format(27h the final statistics after,f8.1,6h years)
    write(9,300)stat
    49 stop
    end program
    function acc(g,x,y,z)
    r2=x*x+y*y+z*z
    acc=-g*x/r2**1.5
    end function acc

    subroutine initial(x,y,z,vx,vy,vz,n,r,key)
    real :: one
    dimension x(100),y(100),z(100),vx(100),vy(100),vz(100),t(3)
!  the random number generator is initiated by a seed between 0 and 1.
!  the following data statement can be changed if a new seed is required.
    data seed/0.234/
    open(unit=9,file='lpt1')
    pi=4*atan(1.0)
    g=4*pi*pi
    one=1.0
    ran=mod((seed+pi)**5,one)
    do 1 i=1,n
    !  position each star at random position within spherical volume
        4 do 2 j=1,3
            ran=mod((ran+pi)**5,one)
            t(j)=2*ran-1
        2 end do
        if(t(1)**2+t(2)**2+t(3)**2 > 1.0)goto 4
        x(i)=r*t(1)
        y(i)=r*t(2)
        z(i)=r*t(3)
    1 end do
!  calculate potential energy
    poten=0
    do 3 i=1,n-1
        do 3 j=i+1,n
            rij=sqrt((x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2)
            poten=poten-g/rij
    3 end do
!  calculate mean speed from the virial theorem
    vm=sqrt(-poten/n)
!  give every star the same speed but in a random direction
    do 5 i=1,n
        7 do 6 j=1,3
            ran=mod((ran+pi)**5,one)
            t(j)=2*ran-1
        6 end do
        tt=t(1)**2+t(2)**2+t(3)**2
        if(tt > 1)goto 7
        st=sqrt(tt)
        vx(i)=vm*t(1)/st
        vy(i)=vm*t(2)/st
        vz(i)=vm*t(3)/st
    5 end do
    100 format(6e11.4)
!  now add or subtract constant dx,dy,dz,dvx,dvy,dvz from each star to
!  put the centre of mass at the origin and to give zero momentum
    sumx=0
    sumy=0
    sumz=0
    sumvx=0
    sumvy=0
    sumvz=0
    do 10 i=1,n
        sumx=sumx+x(i)
        sumy=sumy+y(i)
        sumz=sumz+z(i)
        sumvx=sumvx+vx(i)
        sumvy=sumvy+vy(i)
        sumvz=sumvz+vz(i)
    10 end do
    do 12 i=1,n
        x(i)=x(i)-sumx/n
        y(i)=y(i)-sumy/n
        z(i)=z(i)-sumz/n
        vx(i)=vx(i)-sumvx/n
        vy(i)=vy(i)-sumvy/n
        vz(i)=vz(i)-sumvz/n
    12 end do
!  find factor to restore rms speed to vm
    sum=0
    do 13 i=1,n
        sum=sum+vx(i)**2+vy(i)**2+vz(i)**2
    13 end do
    factor=vm/sqrt(sum/n)
    do 14 i=1,n
        vx(i)=factor*vx(i)
        vy(i)=factor*vy(i)
        vz(i)=factor*vz(i)
    14 end do
    if(key == 0)write(9,100)(x(i),y(i),z(i),vx(i),vy(i),vz(i),i=1,n)
    return
    end subroutine initial

    subroutine stats(x,y,z,vx,vy,vz,n,g,stat)
    dimension x(100),y(100),z(100),vx(100),vy(100),vz(100),stat(4)
    sumr2=0
    sumx=0
    sumy=0
    sumz=0
    sumvx=0
    sumvy=0
    sumvz=0
    sumvxyz2=0
    sumpot=0
    do 1 i=1,n
        sumx=sumx+x(i)
        sumy=sumy+y(i)
        sumz=sumz+z(i)
        sumr2=sumr2+x(i)**2+y(i)**2+z(i)**2
        sumvx=sumvx+vx(i)
        sumvy=sumvy+vy(i)
        sumvz=sumvz+vz(i)
        sumvxyz2=sumvxyz2+vx(i)**2+vy(i)**2+vz(i)**2
        if(i == n)goto 1
        do 2 j=i+1,n
            rij=sqrt((x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2)
            sumpot=sumpot-g/rij
        2 end do
    1 end do
    stat(1)=sumr2
    stat(2)=sqrt(sumx**2+sumy**2+sumz**2)/n
    stat(3)=sqrt(sumvx**2+sumvy**2+sumvz**2)/n
    stat(4)=0.5*sumvxyz2+sumpot
    return
    end subroutine stats
