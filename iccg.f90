    subroutine iccg(x,y,r,s,t,m,n,acc,acc1)

!  solves the matrix equation:
!  -r(k+m)*x(k+m)-s(k+1)*x(k+1)+t(k)*x(k)-s(k)*x(k-1)-r(k)*x(k-m)=y(k)
!  using iccg
!  the routine is based on the iccg solution of sparse banded matrices
!  described by kershaw (j comp phys 26, 43-65, 1978)
!  the method works by forming an approximate cholesky decomposition matrix
!  (r1,s1,t1) of the square banded symmetric matrix (r,s,t) with the same
!  sparsity as the original matrix, the inverse is readily calculated by
!  the usual triangular matrix inversion. by multiplying the input set of
!  equations by the approximate inverse a new set of equations is formed
!  with the same solution but with a matrix form which is close to the
!  identity matrix.
!  the final solution is rapidly formed by a conjugate gradient iteration
!  in a relatively small number of steps.
!  the calculation is terminated at an accuracy:
!  absolute acc1 and relative acc for each term.
!  since the original banded forms are retained throughout the necessary
!  matrix multiplications are rapidly evaluated.
!  the method is easily generalised to more highly structured banded forms.
!  the final result is contained in the x array

    parameter (ilt=200,jlt=200)
    dimension x(n),y(n),r(n),s(n),t(n)
    logical :: check
    dimension p(ilt*jlt),q(ilt*jlt),q1(ilt*jlt),r1(ilt*jlt), &
    s1(ilt*jlt),t1(ilt*jlt)

!  calculate the cholesky decomposition matrix (inverse matrix)

    k0=1-m
    k1=1
    t1(1)=t(1)
    do 1  k=2,n
        k0=k0+1
        s1(k)=-s(k)
        s(k)=-s(k)
        if (k0) 10,10,11
        10 t1(k)=t(k)-s1(k)**2/t1(k1)
        go to 12
        11 r1(k)=-r(k)
        r(k)=-r(k)
        if (m == 2) s1(k)=s1(k)-r1(k)*s1(k1)/t1(k0)
        t1(k)=t(k)-s1(k)**2/t1(k1)-r1(k)**2/t1(k0)
    !  set a limiting value for the diagonal to prevent later problems
        help=1.0e-20
        help2=abs(t(k))*.001
        12 a=max(help,help2)
        if (abs(t1(k)) < a) t1(k)=a
        k1=k
    1 end do
    k0=1-m
    k1=1
    do 14 k=2,n
        k0=k0+1
        s1(k)=s1(k)/t1(k1)
        if (k0) 14,14,13
        13 r1(k)=r1(k)/t1(k0)
        k1=k
    14 end do

!  perform conjugate gradient iteration

!  x contains the updated solution vector
!  q contains the residual (error) vector
!  p contains the intermediate vector

!  prepare working arrays for entry to conjugate gradient iteration
!  form the product of the entry solution with the input matrix
    call prod(x,q,r,s,t,m,n)
!  and generate the residual
    do 20 k=1,n
        q(k)=y(k)-q(k)
    20 end do
!  form the approximate inverse as the initial intermediate vector
    call invert(q,p,r1,s1,t1,m,n)
!  and the scalar product of the residual and its inverse
    b=0.0
    do 21 k=1,n
        b=b+p(k)*q(k)
    21 end do
!  if value small calculation complete
    if (abs(b) < 1.0e-6) return

!  entry to the iteration

!  number of iterations limited to number of equations
    do 3 l=1,n
    !  form the product of the intermediate vector
    !  to give the gradient direction for the residual
        call prod(p,q1,r,s,t,m,n)
    !  and its scalar product with itself
        a=0.0
        do 22 k=1,n
            a=a+p(k)*q1(k)
        22 end do
    !  to determine the gradient correction term
        a=b/a
        check= .true. 
    !  adjust the  solution and residual vectors in the gradient direction
        do 23 k=1,n
            x(k)=x(k)+a*p(k)
            q(k)=q(k)-a*q1(k)
        !  iteration test on residual
            if (check) check=abs(q(k)) < amax1(acc1,acc*abs(y(k)))
        23 end do
    !  if accuracy level not reached continue iteration
        if ( .not. check) go to 25
    !  check against original matrix equation
        call prod(x,q1,r,s,t,m,n)
        do 24 k=1,n
            if (check) check=abs(y(k)-q1(k)) < amax1(acc1,acc*abs(y(k)))
        24 end do
    !   if accuracy level maintained exit from routine
    !  principal exit
        if (check) write(6,*)l
        if (check) return
    !  continue iteration by obtaining approximate inverse of residual
        25 call invert(q,q1,r1,s1,t1,m,n)
        a=b
    !  and its scalar product with itself
        b=0.0
        do 26 k=1,n
            b=b+q(k)*q1(k)
        26 end do
        if (abs(b) < 1.0e-10) return
    !  to generate the next iterate for the intermediate array
        a=b/a
        do 2 k=1,n
            p(k)=q1(k)+a*p(k)
        2 end do
    3 end do
!  end of iterative loop
!  exit if iteration failed in maximum number of cycles
    return
    end subroutine iccg

    subroutine invert(x,y,r,s,t,m,n)
!  forms the approximate inverse y of the right-hand side vector x
!  from the symmetric band triangular forms (r,s,t)
    dimension x(n),y(n),r(n),s(n),t(n)
!  the first (forward) triangular sweep
    k1=1
    k0=1-m
    y(1)=x(1)
    do 1 k=2,n
        k0=k0+1
        if (k0) 10,10,11
        10 y(k)=(x(k)-s(k)*y(k1))
        go to 1
        11 y(k)=(x(k)-s(k)*y(k1)-r(k)*y(k0))
        k1=k
    1 end do
    do 2 k=1,n
        y(k)=y(k)/t(k)
    2 end do
!  the second (backward) triangular sweep
    k=n
    k0=n+m
    do 3 kk=2,n
        k1=k
        k=k-1
        k0=k0-1
        if (k0 <= n) go to 30
        y(k)=(y(k)-s(k1)*y(k1))
        go to 3
        30 y(k)=(y(k)-s(k1)*y(k1)-r(k0)*y(k0))
    3 end do
    return
    end subroutine invert

    subroutine prod(x,y,r,s,t,m,n)
! forms the product vector y from the vector x and the symmetric band
! form (r,s,t)
    dimension x(n),y(n),r(n),s(n),t(n)
    k0=2-m
    k3=1+m
    k1=1
    k=2
    y(1)=t(1)*x(1)+s(2)*x(2)+r(k3)*x(k3)
    do 1 k2=3,n
        y(k)=t(k)*x(k)+s(k)*x(k1)+s(k2)*x(k2)
        if (k0) 11,11,10
        10 y(k)=y(k)+r(k)*x(k0)
        11 if (k3 >= n) go to 12
        k3=k3+1
        y(k)=y(k)+r(k3)*x(k3)
        12 k0=k0+1
        k1=k
        k=k2
    1 end do
    y(n)=t(n)*x(n)+s(n)*x(k1)+r(n)*x(k0)
    return
    end subroutine prod
