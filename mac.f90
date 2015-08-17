    program mac
!  a simple two-dimensional program for the application of the marker
!  and cell method.   the enclosure consists of 20 x 20 cells with
!  rigid walls where a full cell contains 16 marker particles initially
!  in a regular grid.   problems should involve no more than 128
!  occupied cells and hence not more than 2048 marker particles.
!  the provided set-up subroutine starts with fluid occupying 10x10
!  cells in the lower left of the container. there is no viscosity
!  included in the simulation. a free-slip boundary condition is
!  applied.
!  the pressure-method algorithm is that due to harlow & welch (1965).
!  up to six data files can be produced - seen.dat for n = 1 to 6 -
!  corresponding to numbers of cycles specified by the user.
    dimension u(-1:21,-1:20),v(-1:20,-1:21),p(-1:20,-1:20),x(2048), &
    y(2048),ivac(-1:20,-1:20),ncell(-1:20,-1:20),ux(0:20,0:19), &
    vx(0:19,0:20),nc(6)
    common x,y
!  fix size of cell, acceleration due to gravity and liquid density.
    data del,g,ro/0.05,9.81,1.0e3/
    open(unit=11,file='see1.dat')
    open(unit=12,file='see2.dat')
    open(unit=13,file='see3.dat')
    open(unit=14,file='see4.dat')
    open(unit=15,file='see5.dat')
    open(unit=16,file='see6.dat')
    npt=1
    write(6,'('' input no of output data files required'')')
    read(5,*)nout
    do 791 i=1,nout
        write(6,400)i
        400 format(44h after how many cycles do you want data file,i2,1h?)
        read(5,*)nc(i)
    791 end do
    time=0
    nocyc=0
!  clear all tables
    do 1 i=-1,20
        do 1 j=-1,21
            u(j,i)=0
            v(i,j)=0
            if(j == 21)goto 1
            p(i,j)=0
    1 end do
!  call routine to set up markers: npart is number of marker particles.
    call setup(del,g,npart)
!  find all vacuum cells
    77 call vacuum(npart,ncell,del,g)
!  find and categorize surface cells
    call surface(ncell,ivac)
!  non-surface cells have ivac=0. surface cells have ivac=1 to 15
!  according to arrangement of neighbouring vacuum cells.

!  calculate timestep.   this is based on the maximum distance moved by
!  any element of fluid being less than one quarter of a cell dimension.
!  gravity is included in the vertical direction and possible pressure
!  gradients in the horizontal direction.
!  find maximum pressure
    pmax=0
    do 57 i=0,19
        do 57 j=0,19
            if(p(i,j) > pmax)pmax=p(i,j)
    57 end do
    umax=0
    vmax=0
    do 2 i=0,19
        do 2 j=0,20
            if(abs(u(j,i)) > umax)umax=abs(u(j,i))
            if(abs(v(i,j)) > vmax)vmax=abs(v(i,j))
    2 end do
    if(pmax < 1.0e-10)pmax=1.0e-10
    tu=ro*del/pmax*(-umax+sqrt(umax**2+pmax/2.0/ro))
    tv=2*(-vmax+sqrt(vmax*vmax+del*g/2))/g
    dt=amin1(tu,tv)
!  calculate pressure to give zero divergence
    call pressure(p,u,v,ncell,ro,g,del,dt)
!  fill in pressure values, ncell and ivac around boundary
    do 79 m=0,19
        p(-1,m)=p(0,m)
        p(20,m)=p(19,m)
        p(m,-1)=p(m,0)
        p(m,20)=p(m,19)
        ncell(-1,m)=ncell(0,m)
        ncell(20,m)=ncell(19,m)
        ncell(m,-1)=ncell(m,0)
        ncell(m,20)=ncell(m,19)
        ivac(-1,m)=ivac(0,m)
        ivac(20,m)=ivac(19,m)
        ivac(m,-1)=ivac(m,0)
        ivac(m,20)=ivac(m,19)
    79 end do
!  advance components of velocity.   there are free-slip boundary walls
!  but components of velocity perpendicular to walls are zero.
!  first deal with components u
    do 4 i=-1,19
        do 4 j=0,19
        !  advance components u(i+1,j) if cell is not a vacuum cell.
            if(ncell(i,j) == 0)goto 4
            a=(u(i+2,j)+u(i,j)+u(i+1,j+1)+u(i+1,j-1))/4
            b=(u(i+2,j)**2-u(i,j)**2)/2/del
            c=(u(i+1,j+1)+u(i+1,j))*(v(i+1,j+1)+v(i,j+1))/4/del
            d=(u(i+1,j-1)+u(i+1,j))*(v(i+1,j)+v(i,j))/4/del
            e=(p(i+1,j)-p(i,j))/del/ro
            ux(i+1,j)=a-(b+c-d+e)*dt
    4 end do
!  now deal with v components
    do 14 i=0,19
        do 14 j=-1,19
        !  advance v(i,j+1) if cell is not a vacuum cell
            if(ncell(i,j) == 0)goto 14
            a=(v(i,j+2)+v(i,j)+v(i+1,j+1)+v(i-1,j+1))/4
            b=(v(i,j+2)**2-v(i,j)**2)/2/del
            c=(v(i+1,j+1)+v(i,j+1))*(u(i+1,j+1)+u(i+1,j))/4/del
            d=(v(i-1,j+1)+v(i,j+1))*(u(i,j+1)+u(i,j))/4/del
            e=(p(i,j+1)-p(i,j))/del/ro
            z=g
            if(j+1 == 0)z=0
            vx(i,j+1)=a-(b+c-d+e+z)*dt
    14 end do
!  now deal with surface cells.  in the above some of their velocity
!  components would have been modified but must now be changed to give
!  correct boundary conditions
    do 8 i=0,19
        do 8 j=0,19
            if(ivac(i,j) == 0)goto 8
        !  surface cell identified
            k=ivac(i,j)
            call boundary(i,j,u,v,ux,vx,k,g,dt)
    8 end do
!  transfer new values
    do 66 i=0,20
        do 66 j=0,20
            if(j == 20)goto 67
            u(i,j)=ux(i,j)
            ux(i,j)=0
            67 if(i == 20)goto 66
            v(i,j)=vx(i,j)
            vx(i,j)=0
    66 end do
!  fill in table edges to give conditions for free-slip boundaries.
    do 78 m=0,20
    !  motions parallel to boundary on opposite sides are equal.
        u(m,-1)=u(m,0)
        u(m,20)=u(m,19)
        v(-1,m)=v(0,m)
        v(20,m)=v(19,m)
    !  motions perpendicular to the boundary on opposite sides are reversed.
        u(-1,m)=-u(1,m)
        u(21,m)=-u(19,m)
        v(m,-1)=-v(m,1)
        v(m,21)=-v(m,19)
    78 end do
!  now move marker particles
    call move(u,v,ncell,npart,del,dt)
    time=time+dt
    nocyc=nocyc+1
    if(nocyc /= nc(npt))goto 77
    nf=10+npt
    do 30 i=1,npart
        write(nf,*)x(i),y(i)
    30 end do
    write(6,*)time
    npt=npt+1
    if(npt <= nout)goto 77
    stop
    end program


    subroutine boundary(i,j,u,v,ux,vx,k,g,dt)
    dimension u(-1:21,-1:20),v(-1:20,-1:21),ux(0:20,0:19), &
    vx(0:19,0:20)
    goto(3,6,3,4,5,6,5,8,9,15,15,8,5,15,15)k
! apply boundary condition according to the value of k
    15 ux(i,j)=u(i,j)
    ux(i+1,j)=u(i+1,j)
    if(k == 11)vx(i,j+1)=vx(i,j)
    if(k == 10 .or. k == 11)goto 100
    if(k == 14)then
        vx(i,j)=vx(i,j+1)
        goto 100
    endif
    5 vx(i,j)=v(i,j)-g*dt
    vx(i,j+1)=v(i,j+1)-g*dt
    if(k == 13)goto 8
    if(k == 15 .or. k == 5)goto 100
    if(k == 7)then
        ux(i,j)=ux(i+1,j)
        goto 100
    endif
    8 ux(i+1,j)=ux(i,j)
    if(k == 12)vx(i,j)=vx(i,j+1)
    goto 100
    9 ux(i+1,j)=ux(i,j)
    3 vx(i,j+1)=vx(i,j)
    if(k == 9 .or. k == 1)goto 100
    ux(i,j+1)=ux(i,j)
    goto 100
    6 ux(i,j)=ux(i+1,j)
    if(k == 2)goto 100
    4 vx(i,j)=vx(i,j+1)
    100 return
    end subroutine boundary


    subroutine vacuum(npart,ncell,del,g)
    dimension ncell(-1:20,-1:20),x(2048),y(2048)
    common x,y
!  clear all cells
    do 1 i=0,19
        do 1 j=0,19
            ncell(i,j)=0
    1 end do
!  find number of particles in each cell.  the particles are constrained
!  to be inside the enclosure.
    do 2 k=1,npart
        i=int(x(k)/del)
        if(i*(i-19) > 0)goto 100
        j=int(y(k)/del)
        if(j*(j-19) > 0)goto 100
        ncell(i,j)=ncell(i,j)+1
    2 end do
!  vacuum cells (i,j) are now recognized by having ncell(i,j)=0.
    goto 50
    100 write(6,'('' error.  particle detected outside enclosure.'')')
    stop
    50 return
    end subroutine vacuum


    subroutine surface(ncell,ivac)
    dimension ivac(-1:20,-1:20),ncell(-1:20,-1:20),x(2048),y(2048)
    common x,y
!  clear array ivac
    do 1 i=0,19
        do 1 j=0,19
            ivac(i,j)=0
    1 end do
!  value of ivac gives the type of non-vacuum cell
!  ivac=0 is a non-surface cell
!  ivac = 1 to 15 describes the type of vacuum cell
    do 2 i=0,19
        do 2 j=0,19
            if(ncell(i,j) == 0)goto 2
            if(j == 19)goto 3
            if(ncell(i,j+1) == 0)ivac(i,j)=ivac(i,j)+1
            3 if(i == 0)goto 4
            if(ncell(i-1,j) == 0)ivac(i,j)=ivac(i,j)+2
            4 if(j == 0)goto 5
            if(ncell(i,j-1) == 0)ivac(i,j)=ivac(i,j)+4
            5 if(i == 19)goto 2
            if(ncell(i+1,j) == 0)ivac(i,j)=ivac(i,j)+8
    2 end do
    return
    end subroutine surface


    subroutine gaussp(a,b,x,n)
!  this solves a set of linear equations, ax = b, when the matrix has
!  a strong diagonal.   the sor (successive over-relaxation) method is
!  used, which is the same as gauss-seidel if the over-relaxation
!  factor is 1.  up to 128 equations can be handled with array dimensions
!  provided.
    real :: a(128,128),b(128),x(128)
!  set tolerance
    tol=1.e-4
!  put initial estimate of solution elements all equal zero
    do 1 i=1,n
        x(i)=0
    1 end do
!  set over-relaxation factor to 1.5
    w=1.5
    icycle=0
    5 icycle=icycle+1
    dif=0
    do 7 i=1,n
        sum=b(i)
        do 3 j=1,n
            if(j == i)goto 3
            sum=sum-a(i,j)*x(j)
        3 end do
        est=sum/a(i,i)
        z=w*(est-x(i))
        x(i)=x(i)+z
        if(abs(z) > dif)dif=abs(z)
    7 end do
    if(dif < tol)goto 4
!  not more than 100 cycles allowed
    if(icycle < 100)goto 5
    4 return
    end subroutine gaussp


    subroutine setup(del,g,no)
    dimension x(2048),y(2048)
    common x,y
!  this puts 16 markers in each of an array of 10x10 cells
    no=0
    do 1 i=0,19
        jz=9
        if(i > 9)goto 1
        do 2 j=0,jz
            do 3 k=0,15
                l=k/4
                m=k-4*l
                no=no+1
                x(no)=(i+0.125+0.25*float(l))*del
                y(no)=(j+0.125+0.25*float(m))*del
            3 end do
        2 end do
    1 end do
    return
    end subroutine setup



    subroutine pressure(p,u,v,ncell,ro,g,del,dt)
    dimension u(-1:21,-1:20),v(-1:20,-1:21),rhs(128),xx(128), &
    ncell(-1:20,-1:20),p(-1:20,-1:20)
    real :: lhs(128,128)
    integer :: pcode(128,2),rp(0:19,0:19)
!  find non-vacuum cells for pressure determination.  allocate sequential
!  code numbers to such cells. array pcode gives the (i,j) for the code
!  number and array rp gives the code number from the (i,j).
!  clear table
    do 1 i=1,128
        do 1 j=1,2
            pcode(i,j)=0
    1 end do
    num=0
    do 2 i=0,19
        do 2 j=0,19
            rp(i,j)=0
            if(ncell(i,j) == 0)goto 2
            num=num+1
            rp(i,j)=num
            pcode(num,1)=i
            pcode(num,2)=j
    2 end do
!  clear lhs which will contain the coefficients of the linear equations
!  for solving for pressures.
    do 10 i=1,128
        do 10 j=1,128
            lhs(i,j)=0
    10 end do
!  find coefficients of equations for solving for pressures
    do 3 k=1,num
        i=pcode(k,1)
        j=pcode(k,2)
        sum=0
        sum=sum+(u(i+2,j)**2-u(i,j)**2)/2/del
        sum=sum-(u(i+1,j)**2-u(i-1,j)**2)/2/del
        sum=sum+(v(i,j+2)**2-v(i,j)**2)/2/del
        sum=sum-(v(i,j+1)**2-v(i,j-1)**2)/2/del
        sum=sum+(u(i+1,j+1)+u(i+1,j))*(v(i+1,j+1)+v(i,j+1))
        sum=sum-(v(i,j)+v(i+1,j))*(u(i+1,j)+u(i+1,j-1))
        sum=sum-(u(i,j+1)+u(i,j))*(v(i,j+1)+v(i-1,j+1))
        sum=sum+(v(i-1,j)+v(i,j))*(u(i,j)+u(i,j-1))
        sump=0
        call dd(u,v,i+1,j,delv)
        sump=sump+delv
        call dd(u,v,i-1,j,delv)
        sump=sump+delv
        call dd(u,v,i,j+1,delv)
        sump=sump+delv
        call dd(u,v,i,j-1,delv)
        sump=sump+delv
        rhs(k)=ro*(del*sump/4/dt-sum/2)
        if(j == 0)rhs(k)=rhs(k)-g*ro*del
        lhs(k,k)=-4
        if(i+1 > 19)then
            lhs(k,k)=lhs(k,k)+1
            goto 4
        endif
        l=rp(i+1,j)
        if(l == 0)goto 4
        lhs(k,l)=1
        4 if(i-1 < 0)then
            lhs(k,k)=lhs(k,k)+1
            goto 5
        endif
        l=rp(i-1,j)
        if(l == 0)goto 5
        lhs(k,l)=1
        5 if(j+1 > 31)then
            lhs(k,k)=lhs(k,k)+1
            goto 6
        endif
        l=rp(i,j+1)
        if(l == 0)goto 6
        lhs(k,l)=1
        6 if(j-1 < 0)then
            lhs(k,k)=lhs(k,k)+1
            goto 3
        endif
        l=rp(i,j-1)
        if(l == 0)goto 3
        lhs(k,l)=1
    3 end do
!  coefficients for all equations found
    call gaussp(lhs,rhs,xx,num)
    do 20 k=1,num
        i=pcode(k,1)
        j=pcode(k,2)
        p(i,j)=xx(k)
    20 end do
    return
    end subroutine pressure


    subroutine move(u,v,ncell,npart,del,dt)
    dimension x(2048),y(2048),u(-1:21,-1:20),v(-1:20,-1:21), &
    ncell(-1:20,-1:20)
    common x,y
    do 1 k=1,npart
    !  find horizontal velocity cell withing which particle resides and
    !  fractional coordinates in that cell
        a=x(k)/del
        b=y(k)/del-0.5
        if(b < 0)b=0
        i=int(a)
        j=int(b)
        alf=a-i
        bet=b-j
        uu=(1-alf)*(1-bet)*u(i,j)+(1-alf)*bet*u(i,j+1)+ &
        alf*(1-bet)*u(i+1,j)+alf*bet*u(i+1,j+1)
    !  test if there is boundary in x direction.   if so just take velocity
    !  of edge away from boundary.
    !    2 if(i.eq.19)then
    !      uu=u(i,j)
    !      goto 3
    !      endif
    !      uu=(1-alf)*u(i,j)+alf*u(i+1,j)
    !  find vertical velocity cell in which particle resides and coordinates
    !  within that cell.
        3 a=x(k)/del-0.5
        if(a < 0)a=0
        b=y(k)/del
        i=int(a)
        j=int(b)
        alf=a-i
        bet=b-j
        vv=(1-alf)*(1-bet)*v(i,j)+alf*(1-bet)*v(i+1,j)+ &
        (1-alf)*bet*v(i,j+1)+alf*bet*v(i+1,j+1)
        goto 10
    !  test if there is a boundary in the y direction.  if so just take
    !  the velocity of the edge away from the boundary.
    !    4 if(j.eq.19)then
    !      vv=v(i,j)
    !      goto 10
    !      endif
    !      vv=(1-bet)*v(i,j)+bet*v(i,j+1)
    !  move particle.   if it goes outside boundary then reflect it
    !  from boundary.
        10 x(k)=x(k)+uu*dt
        if(x(k) >= 20.0*del)x(k)=40.0*del-x(k)-1.0e-3
        if(x(k) <= 0)x(k)=-x(k)+1.0e-3
        y(k)=y(k)+vv*dt
        if(y(k) >= 20.0*del)y(k)=40.0*del-y(k)-1.0e-3
        if(y(k) < 0)y(k)=-y(k)+1.0e-3
    1 end do
    return
    end subroutine move
               
    subroutine dd(u,v,i,j,delv)
    dimension u(-1:21,-1:20),v(-1:20,-1:21)
    delv=v(i,j+1)-v(i,j)+u(i+1,j)-u(i,j)
    return
    end subroutine dd
