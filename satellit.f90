!  programme satellite

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
!  (v)   subroutine "out" which outputs the results to a data file.

!  by changing the subroutines different problems may be solved.
!  the cm's can be masses or charges or be made equal to unity while
!  the force law can be inverse-square or anything else -
!  e.g. lennard-jones. see comment at beginning of subroutine
!  "acc" for the types of force operating.

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
!  system of units being used.

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

! give output to monitor passage of simulated time
    help=1000.0
    if(abs(mod(time,help)) < h)write(6,66)time
    66 format(15h simulated time, f9.0)
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
    dimension x(20,3),v(20,3),astore(1000),tstore(1000), &
    estore(1000),xtemp(2,20,3),vtemp(2,20,3),cm(20),pos(3), &
    vel(3),xt(20,3),vt(20,3),delv(20,3)
    common/a/x,v,tol,h,totime,delv,xt,vt,nb,ist,time,ig,xtemp, &
    vtemp,cm
    common/b/d,ecc,el0,norig,astore,tstore,estore
    data g/6.667e-11/

!  the values of a and e are calculated every 100 steps
!  and are stored together with the time.

    ist=ist+1
    if((ist/100)*100 /= ist)goto 50
    ig=ist/100
    if(ig > 1000)goto 50

!  first find position and velocity of centre of mass of the
!  satellite.

    do 1 k=1,3
        pos(k)=0
        vel(k)=0
        do 2 j=2,nb
            pos(k)=pos(k)+x(j,k)
            vel(k)=vel(k)+v(j,k)
        2 end do
        pos(k)=pos(k)/(nb-1.0)
        vel(k)=vel(k)/(nb-1.0)
    1 end do

!  calculate orbital distance

    r=sqrt(pos(1)**2+pos(2)**2+pos(3)**2)

!  calculate v**2

    v2=vel(1)**2+vel(2)**2+vel(3)**2

!  calculate intrinsic energy

    totm=cm(1)+cm(2)+cm(3)+cm(4)
    en=-g*totm/r+0.5*v2

!  calculate semi-major axis

    astore(ig)=-0.5*g*totm/en

!  calculate square of intrinsic angular momentum

    d1=pos(2)*vel(3)-pos(3)*vel(2)
    d2=pos(3)*vel(1)-pos(1)*vel(3)
    d3=pos(1)*vel(2)-pos(2)*vel(1)
    h2=d1**2+d2**2+d3**2

!  calculate eccentricity

    e2=1.0-h2/g/totm/astore(ig)
    if(e2 < 0.0)then
        write(6,'('' error - e**2 is negative'')')
        stop
    endif
    estore(ig)=sqrt(e2)
    tstore(ig)=time
    50 return
    end subroutine store


    subroutine start
    dimension x(20,3),v(20,3),astore(1000),tstore(1000), &
    estore(1000),xtemp(2,20,3),vtemp(2,20,3),cm(20),xt(20,3), &
    vt(20,3),delv(20,3)
    common/a/x,v,tol,h,totime,delv,xt,vt,nb,ist,time,ig,xtemp, &
    vtemp,cm
    common/b/d,ecc,el0,norig,astore,tstore,estore
    data g/6.67e-11/
    open(unit=20,file='sat.dat')
    write(6,'('' input the number of bodies'')')
    read(5,*)nb
    write(6,'('' input the values of cm; for programme satellite '')')
    write(6,'('' they are the mass of each body in kg.'')')
    do 1 i=1,nb
        write(6,100)i
        100 format(25h read in the value of cm[,i3, 1h])
        read(5,*)cm(i)
    1 end do
    write(6,'('' input the initial timestep in seconds'')')
    read(5,*)h
    write(6,'('' input the total time for the simulation '')')
    write(6,'('' in seconds.'')')
    read(5,*)totime
    write(6,'('' the calculation can be done relative to an '')')
    write(6,'('' arbitrary origin or with respect to one of '')')
    write(6,'('' the bodies as origin.   input zero for an '')')
    write(6,'('' arbitrary origin or the number of the body.'')')
    write(6,'('' if a body is chosen as origin then all its'')')
    write(6,'('' positional and velocity values are set to zero'')')
    read(5,*)norig

!  the user is required to specify a tolerance, the maximum absolute
!  error that can be tolerated in any positional coordinate (x, y or z).
!  if this is set too low then the programme can become very slow.
!  for this problem a tolerance of 100m should be acceptable.

    write(6,'('' input the tolerance in metres'')')
    write(6,'('' [see comment statement in listing]'')')
    read(5,*)tol
    write(6,'('' input initial distance of satellite in metres'')')
    read(5,*)d
    write(6,'('' input unstretched length of springs in metres'')')
    read(5,*)el0
    write(6,'('' input initial eccentricity'')')
    read(5,*)ecc
    dd=0.5*el0/sqrt(3.0)
    x(1,1)=0
    x(1,2)=0
    x(1,3)=0
    x(2,1)=d-2*dd
    x(2,2)=0
    x(2,3)=0
    x(3,1)=d+dd
    x(3,2)=0.5*el0
    x(3,3)=0
    x(4,1)=d+dd
    x(4,2)=-0.5*el0
    x(4,3)=0
    dis2=x(2,1)
    dis3=sqrt(x(3,1)**2+x(3,2)**2)
    angle=atan2(x(3,2),x(3,1))

!  initial velocities are calculated for the three components of the
!  satellite so that the spin angular velocity of the satellite is
!  approximately equal to the orbital angular velocity.

    vv=sqrt(g*(cm(1)+cm(2)+cm(3)+cm(4))*(1-ecc)/d)
    v(2,2)=vv*dis2/d
    v(3,2)=vv*dis3*cos(angle)/d
    v(4,2)=v(3,2)
    v(1,1)=0
    v(1,3)=0
    v(2,1)=0
    v(2,3)=0
    v(3,1)=-vv*dis3*sin(angle)/d
    v(4,1)=-v(3,1)
    v(3,3)=0
    v(4,3)=0
    return
    end subroutine start




    subroutine acc
    dimension delv(20,3),r(3),xt(20,3),vt(20,3),dif(3), &
    difd(3)
    dimension x(20,3),v(20,3),astore(1000),tstore(1000), &
    estore(1000),xtemp(2,20,3),vtemp(2,20,3),cm(20)
    common/a/x,v,tol,h,totime,delv,xt,vt,nb,ist,time,ig,xtemp, &
    vtemp,cm
    common/b/d,ecc,el0,norig,astore,tstore,estore


!  in this programme normal si units are used.   three kinds of force
!  are operating.   the first is the normal gravitational inverse-square
!  law between all pairs of bodies, the second is due to the elasticity
!  of the material as modelled by the set of three springs and the third
!  operates between the three component bodies of the satellite and
!  depends on the rate at which the springs expand or contract.  the
!  third force provides the dissipation in the system.

!  set g, spring constant, cay,and dissipation constant, cee.

    data g,cay,cee/6.667e-11,3.0e18,1.0e22/
    do 1 j=1,nb
        do 1 k=1,3
            delv(j,k)=0
    1 end do

!  the following pair of do loops finds interactions for all pairs
!  of bodies.   only the gravitational forces are considered in this
!  section

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

    if(norig == 0)goto 60
    do 6 j=1,nb
        do 6 k=1,3
            delv(j,k)=delv(j,k)-delv(norig,k)
    6 end do

!  now we put in the spring and dissipative forces

    do 10 i=2,nb-1
        do 10 j=i+1,nb
            do 11 k=1,3
                dif(k)=xt(i,k)-xt(j,k)
            11 end do
        
        !  calculate length of spring
        
            elp=sqrt(dif(1)**2+dif(2)**2+dif(3)**2)
        
        !  calculate dl/dt
        
            do 12 k=1,3
                difd(k)=vt(i,k)-vt(j,k)
            12 end do
            dlbdt=(dif(1)*difd(1)+dif(2)*difd(2)+dif(3)*difd(3))/elp
            force=cay*(elp-el0)+cee*dlbdt
            do 13 k=1,3
                delv(i,k)=delv(i,k)-force*dif(k)/elp/cm(i)
                delv(j,k)=delv(j,k)+force*dif(k)/elp/cm(j)
            13 end do
    10 end do
    60 return
    end subroutine acc







    subroutine out
    dimension x(20,3),v(20,3),astore(1000),tstore(1000), &
    estore(1000),xtemp(2,20,3),vtemp(2,20,3),cm(20),xt(20,3), &
    vt(20,3),delv(20,3)
    common/a/x,v,tol,h,totime,delv,xt,vt,nb,ist,time,ig,xtemp, &
    vtemp,cm
    common/b/d,ecc,el0,norig,astore,tstore,estore
!  values of (eccentricity, time) are placed in file sat.dat in a form
!  suitable for input to a graphics package.

    rewind 20
    do 60 i=1,ig
        tstore(i)=tstore(i)/86400.
        write(20,*)estore(i),tstore(i)
    60 end do

    return
    end subroutine out

