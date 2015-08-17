    program heatri
!  this gives an explicit solution for temperatures along a bar
!  with truncated conical cross section with the ends at fixed
!  temperatures.
!  n is the number of divisions of the bar for the calculation
!  m is the number of divisions of the bar for output values
!  note: n should be a multiple of m
!  cappa is the thermal conductivity of the bar in si units
!  the length of the bar is xl
!  t0 and tn are the fixed temperatures at the two ends of the bar
    parameter (n=100)
    real :: a(3,n-1),c(n-1),y(n-1)
    data t0,tn,m/500,300,10/
    data cappa,xl,pi/200,1.0,3.141592654/
    k=n/m
!  clear arrays
    do 5 i=1,n-1
        c(i)=0
        do 6 j=1,3
            a(j,i)=0
        6 end do
    5 end do
!  set up the coefficients of the tridiagonal matrix
    do 1 i=1,n-1
        if(i == 1)goto 2
        a(1,i)=(2.0-(i-0.5)/n)**2
        2 a(2,i)=-(2.0-(i+0.5)/n)**2-(2.0-(i-0.5)/n)**2
        if(i == n-1)goto 1
        a(3,i)=(2.0-(i+0.5)/n)**2
    1 end do
!  set up coefficients of right-hand-side vector
    c(1)=-(2-0.5/n)**2*t0
    c(n-1)=-(2.0-(n-0.5)/n)**2*tn
    call tridiag(a,c,y,n-1)
!  now calculate heat flow through bar
    q=cappa*pi*(2-0.5/n)**2*(t0-y(1))*1.0e-4*n/xl
!  output results
    write(6,100)(i,y(i),i=k,n-1,k)
    100 format(3(i4,f8.2))
    write(6,200)q
    200 format(18h the heat flow is ,f10.3)
    stop
    end program
