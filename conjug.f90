    program conjug
!  a program for the solution of linear equations with a sparse 5-banded
!  coefficient matrix of the form (for dimension = 10):

!            5 -1  0  0 -1  0  0  0  0  0
!           -1  5 -1  0  0 -1  0  0  0  0
!            0 -1  5 -1  0  0 -1  0  0  0
!            0  0 -1  5 -1  0  0 -1  0  0
!           -1  0  0 -1  5 -1  0  0 -1  0
!            0 -1  0  0 -1  5 -1  0  0 -1
!            0  0 -1  0  0 -1  5 -1  0  0
!            0  0  0 -1  0  0 -1  5 -1  0
!            0  0  0  0 -1  0  0 -1  5 -1
!            0  0  0  0  0 -1  0  0 -1  5

!  as provided the dimension of the matrix can be up to 2000

    dimension nl(10000,2),al(10000),x(2000),b(2000),y(2000)

!      n is the dimension of the matrix
!      l gives the position of the outer band - i.e. element in (1,l+1)
!      io gives output mode - 6 for screen, 9 for printer
!      test gives the precision required in the solution

    data n,l,io,test/2000,4,6,0.000001/
    open(unit=9,file='lpt1')
    do 11 i=1,n
    !  set diagonal elements
        nl(i,1)=i
        nl(i,2)=i
        al(i)=5.0
        if(i == n)goto 11
    !  set elements neighbouring diagonal
        it=n+2*i-1
        nl(it,1)=i
        nl(it,2)=i+1
        nl(it+1,1)=i+1
        nl(it+1,2)=i
        al(it)=-1.0
        al(it+1)=-1.0
        if(i > n-l)goto 11
    !  set outer-band elements
        it=3*n-3+2*i
        nl(it,1)=i
        nl(it,2)=i+l
        nl(it+1,1)=i+l
        nl(it+1,2)=i
        al(it)=-1.0
        al(it+1)=-1.0
    11 end do
    list=5*n-2*l-2
!  set right-hand-side vector
    b(1)=3.0
    b(n)=3.0
    do 22 i=l+1,n-l
        b(i)=1
    22 end do
    do 23 i=2,l
        b(i)=2.0
        b(n-i+1)=2.0
    23 end do
    do 60 i=1,n
        x(i)=0
        y(i)=0
    60 end do
!  obtain solution.  maximum of 5000 iterations allowed
    do 3 i=1,5000
        call conjugat(al,nl,x,b,n,list)
        great=0.0
        do 50 j=1,n
            dif=abs(x(j)-y(j))
            if(dif > great)great=dif
        50 end do
        if(great < test)goto 30
        goto 43
        30 continue
        write(io,200)i
        write(io,100)(x(ii),ii=1,n)
        goto 31
        100 format(8f8.4)
        200 format(23h number of iterations =,i6)
        43 do 53 j=1,n
            y(j)=x(j)
        53 end do
    3 end do
    write(6,'('' too many iterations'')')
    31 stop
    end program
          

    subroutine conjugat(al,nl,x,b,n,list)

!  a general conjugate gradient solver for sparse-matrix linear
!  equations.  the element (nl(k,1),nl(k,2)) has value a(k) where k
!  goes from 1 to neq and neq is the total number of non-zero elements.
!  the subroutine can handle matrices up to dimension 'ndim' and up to
!  'nelem' non-zero elements.

    parameter (ndim=2000,nelem=10000)
    dimension al(nelem),nl(nelem,2),x(ndim),b(ndim),u(ndim),c(ndim), &
    d(ndim),e(ndim)
    do 1 i=1,n
        u(i)=0
        c(i)=0
        d(i)=0
        e(i)=0
    1 end do
    do 2 i=1,list
        n1=nl(i,1)
        n2=nl(i,2)
        c(n1)=c(n1)+al(i)*x(n2)
    2 end do
    do 3 i=1,n
        d(i)=c(i)-b(i)
    3 end do
    do 4 i=1,list
        n1=nl(i,1)
        n2=nl(i,2)
        u(n2)=u(n2)+al(i)*d(n1)
    4 end do
    do 5 i=1,list
        n1=nl(i,1)
        n2=nl(i,2)
        e(n1)=e(n1)+al(i)*u(n2)
    5 end do
    sum1=0
    sum2=0
    do 6 i=1,n
        sum1=sum1+e(i)*d(i)
        sum2=sum2+e(i)*e(i)
    6 end do
    elam=sum1/sum2
    do 7 i=1,n
        x(i)=x(i)-elam*u(i)
    7 end do
    return
    end subroutine conjugat

