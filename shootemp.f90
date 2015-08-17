    program shootemp
!  this provides a solution for a boundary-value problem for a
!  first-order ode with one unknown parameter to be determined.
!  the ode is of the form     dy/dx = f(x,q)  where q is the
!  unknown parameter.   the boundary conditions are y = ya for
!  x = a and y = yb for x = b.
!  the function f(x,q) is provided as a function statement.
    real :: y(0:1000),w(4),dy(4)
    character ans*1
    data w/0.0,0.5,0.5,1.0/
    fun(x,q)=-15.915494*q/(2-x)**2
    write(6,'('' input the first boundary condition as a,ya'')')
    write(6,'('' where y = ya when x = a.'')')
    read(5,*)xi,y(0)
    write(6,'('' input the final boundary condition as b,yb'')')
    write(6,'('' where y = yb when x = b.'')')
    read(5,*)xf,yf
    write(6,'('' input the integration step length h in the '')')
    write(6,'('' form of an integer n where h = [b - a]/n. '')')
    write(6,'('' if you want output of y at intervals of [b-a]/m'')')
    write(6,'('' then make n a multiple of m.'')')
    read(5,*)n
    h=(xf-xi)/n
    write(6,'('' input estimate of the unknown parameter, q.'')')
    read(5,*)q
!  the runge-kutta integration now begins.f
    do 1 i=1,n
        x=(i-1)*h
        do 2 j=1,4
            xx=x+w(j)*h
            dy(j)=h*fun(xx,q)
        2 end do
        y(i)=y(i-1)+(dy(1)+dy(4)+2*(dy(2)+dy(3)))/6.0
    1 end do
    write(6,100)y(n)
    100 format(20h the value of yb is ,f8.2)
    5 write(6,'('' output intermediate values of y? [y/n]'')')
    read(5,50)ans
    50 format(a1)
    if(ans == 'y' .or. ans == 'y')goto 3
    if(ans == 'n' .or. ans == 'n')goto 4
    goto 5
    3 write(6,'('' do you want printed output? [y/n]'')')
    read(5,50)ans
    if(ans == 'y' .or. ans == 'y')then
        iout=9
        open(unit=9, file='lpt1')
        goto 7
    endif
    if(ans == 'n' .or. ans == 'n')then
        iout=6
        goto 7
    endif
    goto 3
    7 write(6,'('' input m where output intervals are h*n/m '')')
    read(5,*)m
    k=n/m
    do 6 i=0,n,k
        xx=i*h
        write(iout,200)xx,y(i)
        200 format(f6.3,2x,f8.2)
    6 end do
    4 stop
    end program
                 
