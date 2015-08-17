    program cholesky

!  sets up the elements of a banded symmetrical matrix of type

!         5 -1  0  0 -1  0  0  0  0  0
!        -1  5 -1  0  0 -1  0  0  0  0
!         0 -1  5 -1  0  0 -1  0  0  0
!         0  0 -1  5 -1  0  0 -1  0  0
!        -1  0  0 -1  5 -1  0  0 -1  0
!         0 -1  0  0 -1  5 -1  0  0 -1
!         0  0 -1  0  0 -1  5 -1  0  0
!         0  0  0 -1  0  0 -1  5 -1  0
!         0  0  0  0 -1  0  0 -1  5 -1
!         0  0  0  0  0 -1  0  0 -1  5

!  in the form required for the subroutine iccg.  the number of
!  equations (dimension of the matrix) is n and m describes the
!  position of the outer bands - i.e. there is a non-zero element
!  at (1,m+1).

    parameter(n=2000,m=4)
    dimension x(n),y(n),r(n),s(n),t(n)
    open(unit=9,file='lpt1')
!  set up matrix elements
    do 1 k=1,n
        t(k)=-5.0
        if(k+m <= n)r(k+m)=-1.0
        if(k+1 <= n)s(k+1)=-1.0
    1 end do
    do 2 i=m+1,n-m
        y(i)=-1
    2 end do
    y(1)=-3.0
    y(n)=-3.0
    do 3 i=2,m
        y(i)=-2.0
        y(n-i+1)=-2.0
    3 end do
    call iccg(x,y,r,s,t,m,n,1.0e-6,1.0e-6)
    write(9,100)(x(i),i=1,n)
    100 format(8f8.4)
    stop
    end program
