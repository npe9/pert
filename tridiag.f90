    subroutine tridiag(a,c,y,n)
!  this programme solves a set of linear equations ay = c where
!  a is a tridiagonal matrix.
!  the parameter max sets the maximum number of equations.
    parameter (max=200)
    real :: work(max),a(3,n),c(n),y(n)
!  forward substitution
    div=a(2,1)
    if(abs(div) < 1.0e-16)goto 50
    y(1)=c(1)/div
    do 1 i=2,n
        work(i)=a(3,i-1)/div
        div=a(2,i)-a(1,i)*work(i)
        if(abs(div) < 1.0e-16)goto 50
        y(i)=(c(i)-a(1,i)*y(i-1))/div
    1 end do
!  back substitution
    do 2 i=n-1,1,-1
        y(i)=y(i)-work(i+1)*y(i+1)
    2 end do
    goto 60
    50 write(6,'('' the tridiagonal algorithm has not worked because'')')
    write(6,'('' it has developed a zero divisor. since realistic'')')
    write(6,'('' physical problems should have solutions this is '')')
    write(6,'('' unexpected.    check the main programme.'')')
    stop
    60 return
    end subroutine tridiag
