    program gaussor
!  this solves a set of linear equations, ax = b, when the matrix has
!  a strong diagonal.   the sor (successive over-relaxation) method is
!  used, which is the same as gauss-seidel if the over-relaxation
!  factor is 1.  up to 25 equations can be handled with array dimensions
!  provided.  the equations should be read in the order which gives
!  the strong diagonal
    real :: a(25,25),b(25),x(25)
    character ans*1
    open(unit=9,file='lpt1')
    write(6,'(''read in number of equations [<=25] '')')
    read(5,*)n
    write(6,'(''read in tolerance '')')
    read(5,*)tol
    write(6,50)n
    50 format(23h read in coefficients -,i3,10h per line.)
    do 1 i=1,n
        write(6,100)i
        100 format(13h read in row ,i2)
        read(5,*)(a(i,j),j=1,n)
    1 end do
    write(6,'(''read in elements of rhs vector'')')
    read(5,*)(b(j),j=1,n)
    write(6,'(''read in first estimate of solution in the form'')')
    write(6,'(''of the vector elements.  making elements equal '')')
    write(6,'(''zero is usually suitable.'')')
    read(5,*)(x(j),j=1,n)
    write(6,'(''read in over-relaxation factor [1.0 to 2.0]'')')
!  all information has now been entered.  the solution process begins.
    read(5,*)w
    10 write(6,'(''do you want printed output? [y/n]'')')
    read(5,250)ans
    250 format(a1)
    if(ans == 'y' .or. ans == 'y')goto 8
    if(ans == 'n' .or. ans == 'n')goto 9
    goto 10
    8 m=9
    goto 12
    9 m=6
    12 icycle=0
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
    4 write(m,150)icycle
    150 format(16h solution after ,i3,8h cycles.)
    write(m,200)(i,x(i),i=1,n)
    200 format(4(i7,f8.3))
    stop
    end program
