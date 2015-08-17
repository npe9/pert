!  programme trojans
    program nbody
!  this is a general n-body programme where inter-body forces between
!  bodies i and j are of the form cm(i)*cm(j)*f(rij,vij) where f is a
!  function of rij, the distance between the bodies, and vij, the
!  relative velocities of the two bodies.
!  the structure of the programme is;

!  (i)   the main programme "nbody"is which includes the runge-kutta
!        routine with automatic step control.
!  (ii)  subroutine "start" which enables input of the initial boundary
!        conditions.
!  (iii) subroutine "acc" which gives the acceleration of each body
!        due to its interactions with all other bodies.
!  (iv)  subroutine "store" which stores intermediate coordinates and
!        velocity components as the computation progresses.
!  (v)   subroutine "out" which outputs the results to data files.

!  by changing the subroutines different problems may be solved.
!  the cm's can be masses or charges or be made equal to unity while
!  the force law can be inverse-square or anything else -
!  e.g. lennard-jones.  see comment at the beginning of subroutine
!  "acc" for the types of forces operating.

!  the four-step runge-kutta algorithm is used. the results of two
!  steps with timestep h are checked against taking one step with
!  timestep 2*h.  if the difference is within the tolerance then the
!  two steps, each of h, are accepted and the steplength is doubled for
!  the next step.   however, if the tolerance is not satisfied then the
!  step is not accepted and one tries again with a halved steplength.
!  it is advisable, but not essential, to start with a reasonable
!  steplength; the programme quickly finds a suitable value.

!  as provided the programme handles up to 20 bodies but this can be
!  changed from 20 to whatever is required in the dimension statement.

    dimension cm(20),x(20,3),v(20,3),dx(20,3,0:4),dv(20,3,0:4),wt(4), &
    xtemp(2,20,3),vtemp(2,20,3),xt(20,3),vt(20,3),delv(20,3)
    common/a/x,v,tol,h,totime,delv,xt,vt,nb,ist,time,ig,xtemp, &
    vtemp,cm
    data wt/0.0,0.5,0.5,1.0/
    ist=0
    open(unit=9,file='lpt1')

!  setting the initial boundary condition can be done either
!  explicitly as values of (x,y,z) and (u,v,w) for each body or
!  can be computed.   this is controlled by subroutine "start".
!  other parameters are also set in "start" which also indicates
!  the system of units being used.

    call start

    time=0
!  initialize arrays
    do 57 i=1,20
        do 57 j=1,3
            do 57 k=0,4
                dx(i,j,k)=0
                dv(i,j,k)=0
    57 end do

!  we now take two steps with step length h followed by one step
!  with step length 2*h from the same starting point but first we
!  store the orginal space and velocity coordinates and timestep.

    25 do 7 it=1,2
        do 8 j=1,nb
            do 8 k=1,3
                xtemp(it,j,k)=x(j,k)
                vtemp(it,j,k)=v(j,k)
        8 end do
        htemp=h
        do 10 nostep=1,3-it
            do 11 i=1,4
                do 12 j=1,nb
                    do 12 k=1,3
                        xt(j,k)=xtemp(it,j,k)+wt(i)*dx(j,k,i-1)
                        vt(j,k)=vtemp(it,j,k)+wt(i)*dv(j,k,i-1)
                12 end do
            
                call acc
            
                do 13 j=1,nb
                    do 13 k=1,3
                        dv(j,k,i)=it*htemp*delv(j,k)
                        dx(j,k,i)=it*htemp*vt(j,k)
                13 end do
            11 end do
            do 14 j=1,nb
                do 14 k=1,3
                    xtemp(it,j,k)=xtemp(it,j,k)+(dx(j,k,1)+dx(j,k,4)+2* &
                    (dx(j,k,2)+dx(j,k,3)))/6.0
                    vtemp(it,j,k)=vtemp(it,j,k)+(dv(j,k,1)+dv(j,k,4)+2* &
                    (dv(j,k,2)+dv(j,k,3)))/6.0
            14 end do
        10 end do
    7 end do

!  the above has made two steps of h and, from the same starting point,
!  a single step of 2*h.   the results are now compared

    do 20 j=1,nb
        do 20 k=1,3
            if(abs(xtemp(1,j,k)-xtemp(2,j,k)) > tol)then
                h=0.5*h
                goto 25
            endif
    20 end do

!  at this stage the double step with h agrees within tolerance with
!  the single step with 2*h.   the timestep will now be tried with
!  twice the value for the next step.   if it is too big then it will
!  be reduced again.

    h=2*h
    do 80 j=1,nb
        do 80 k=1,3
            x(j,k)=xtemp(1,j,k)
            v(j,k)=vtemp(1,j,k)
    80 end do
    time=time+h

    call store

    if(time >= totime)goto 50
    if(ig > 1000)then
        ig=1000
        goto 50
    endif
    goto 25

    50 call out

    stop
    end program


    subroutine store
    dimension cm(20),x(20,3),v(20,3),xstore(1000,20,3), &
    xtemp(2,20,3),vtemp(2,20,3),delv(20,3),xt(20,3),vt(20,3)
    common/a/x,v,tol,h,totime,delv,xt,vt,nb,ist,time,ig,xtemp, &
    vtemp,cm
    common/b/norig,xstore
    do 21 j=1,nb
        do 21 k=1,3
            x(j,k)=xtemp(1,j,k)
            v(j,k)=vtemp(1,j,k)
    21 end do

!  up to 1000 positions are stored.  these are taken every 50 steps.

    ist=ist+1
    if((ist/50)*50 /= ist)goto 50
    ig=ist/50
    if(ig > 1000)goto 50
    do 22 j=1,nb
        do 22 k=1,3
            xstore(ig,j,k)=x(j,k)
    22 end do
    50 return
    end subroutine store


    subroutine start
    dimension cm(20),x(20,3),v(20,3),xstore(1000,20,3), &
    xtemp(2,20,3),vtemp(2,20,3),delv(20,3),xt(20,3),vt(20,3)
    common/a/x,v,tol,h,totime,delv,xt,vt,nb,ist,time,ig,xtemp, &
    vtemp,cm
    common/b/norig,xstore
    open(unit=21,file='last.dat')
    open(unit=22,file='fast.dat')

!  the programme as provided is for an inverse-square law and
!  uses units for which the unit mass is that of the sun, the
!  the unit of distance is the astronomical unit (mean sun-earth
!  distance),the unit of time is the year and the gravitational
!  constant is 4*pi**2.

    write(6,'('' input the number of bodies'')')
    read(5,*)nb
    write(6,'('' input the values of cm in solar-mass units. '')')
    write(6,'('' the trojan asteroid masses can be put as zero'')')
    do 1 i=1,nb
        write(6,500)i
        500 format(25h read in the value of cm[,i3, 1h])
        read(5,*)cm(i)
    1 end do
    write(6,'('' input the initial timestep [years]'')')
    read(5,*)h
    write(6,'('' input total time for the simulation [years]'')')
    read(5,*)totime

!  the programme asks the user to specify a tolerance, the maximum
!  absolute error that can be tolerated in any positional coordinate
!  (x, y or z).  if this is set too low then the programme can become
!  very slow.  for computations involving planets a tolerance of 1.0e-6
!  (c. 150 km) is usually satisfactory.

    write(6,'('' input the tolerance '')')
    write(6,'('' see comment above this statement in listing'')')
    read(5,*)tol
    write(6,'('' the calculation can be done relative to an '')')
    write(6,'('' arbitrary origin or with respect to one of '')')
    write(6,'('' the bodies as origin.   input zero for an '')')
    write(6,'('' arbitrary origin or the number of the body.'')')
    write(6,'('' if a body is chosen as origin then all its'')')
    write(6,'('' positional and velocity values are set to zero'')')
    read(5,*)norig
    do 31 j=1,nb
        write(6,100)j
        100 format(23h input [x,y,z] for body,i3)
        read(5,*)x(j,1),x(j,2),x(j,3)
        write(6,200)j
        200 format(32h input [xdot,ydot,zdot] for body,i3)
        read(5,*)v(j,1),v(j,2),v(j,3)
    31 end do
    return
    end subroutine start




    subroutine acc
    dimension cm(20),x(20,3),v(20,3),xstore(1000,20,3),r(3), &
    xtemp(2,20,3),vtemp(2,20,3),delv(20,3),xt(20,3),vt(20,3),dd(3)
    common/a/x,v,tol,h,totime,delv,xt,vt,nb,ist,time,ig,xtemp, &
    vtemp,cm
    common/b/norig,xstore

!  the programme as provided is for an inverse-
!  square law and uses units for which the unit mass is that of the sun,
!  the unit of distance is the astronomical unit (mean sun-earth distance),
!  the unit of time is the year and the gravitational constant is 4*pi**2.
!  however, the user may modify the subroutine "acc" to change to any other
!  force law and/or any other system of units.

!  set the value of g in astronomical units
    pi=4.0*atan(1.0)
    g=4*pi*pi
    do 1 j=1,nb
        do 1 k=1,3
            delv(j,k)=0
    1 end do
!  the following pair of do loops finds interactions for all pairs
!  of bodies
    do 2 j=1,nb-1
        do 2 l=j+1,nb
            do 3 k=1,3
                r(k)=xt(j,k)-xt(l,k)
            3 end do
            rrr=(r(1)**2+r(2)**2+r(3)**2)**1.5
            do 4 k=1,3
            !  the next two statements give the contributions to the three
            !  components of acceleration due to body j on body l and due to
            !  body l on body j.
                delv(j,k)=delv(j,k)-g*cm(l)*r(k)/rrr
                delv(l,k)=delv(l,k)+g*cm(j)*r(k)/rrr
            4 end do
    2 end do
!  if one of the bodies is to be the origin then its acceleration is
!  subtracted from that of all other bodies.
    if(norig == 0)goto 10
    dd(1)=delv(norig,1)
    dd(2)=delv(norig,2)
    dd(3)=delv(norig,3)
    do 6 j=1,nb
        do 6 k=1,3
            delv(j,k)=delv(j,k)-dd(k)
    6 end do
    10 return
    end subroutine acc


    subroutine out
    dimension cm(20),x(20,3),v(20,3),xstore(1000,20,3), &
    xtemp(2,20,3),vtemp(2,20,3),delv(20,3),xt(20,3),vt(20,3)
    common/a/x,v,tol,h,totime,delv,xt,vt,nb,ist,time,ig,xtemp, &
    vtemp,cm
    common/b/norig,xstore
!     (x,y) values for the leading asteroid are placed in data file
!     last.dat and for the following asteroid in data file fast.dat.

!  for the trojan asteroid problem the jupiter radius vector is
!  rotated to put it on the y axis.   the positions of the asteroids
!  relative to jupiter are plotted.
    last=min(ig,1000)
    do 60 i=1,last
    !  theta is the angle between the x axis and the jupiter radius vector
        theta=atan2(xstore(i,2,2),xstore(i,2,1))
        xstore(i,2,2)=sqrt(xstore(i,2,2)**2+xstore(i,2,1)**2)
        xstore(i,2,1)=0
    !  now the asteroid radius vectors are rotated by pi/2-theta
        do 61 j=1,2
            aa=xstore(i,j+2,1)*sin(theta)-xstore(i,j+2,2)*cos(theta)
            bb=xstore(i,j+2,1)*cos(theta)+xstore(i,j+2,2)*sin(theta)
            xstore(i,j+2,1)=aa
            xstore(i,j+2,2)=bb
        61 end do
    60 end do
!  the modified positions are now output to data files.
    do 63 j=3,4
        n=18+j
        rewind n
        do 64 i=1,last
            write(n,*)xstore(i,j,1),xstore(i,j,2)
        64 end do
    63 end do
    return
    end subroutine out

