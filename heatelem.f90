    program heatelem
!  an application of the finite-element method to a heated uniform
!  rod.   the rod is of length l and its material has thermal
!  conductivity kappa.   the ends of the rod can be at a constant
!  temperature, insulated or exposed to an external temperature.
!  there is provision for heat sources and sinks which can be
!  extended and described by the subroutine source(x,q) or be
!  a point source or sink at a node.
    real :: a(3,199),ll,kappa,c(199),x(7),w(7),y(199), &
    yy(0:200),end(2),xc(0:200)
    dimension kk(5)
    character ans*1,stan*13,cn*6,ex(6)*1
!  data for 7-point gauss integration
    data x/0.025446,0.1292344,0.2970774,0.5,0.7029226,0.8707656, &
    & 0.974554/
    data w/0.0647425,0.1398527,0.190915,0.2089796,0.190915, &
    & 0.1398527,0.0647425/
    stan='0123456789 ie'
    30 write(6,'('' if your program has an extended source or sink'')')
    write(6,'('' then the subroutine source provides this.  if'')')
    write(6,'('' you have forgotten to amend this then abort the'')')
    write(6,'('' program.  if there is an extended source and it'')')
    write(6,'('' is correctly provided then type "y".   if there'')')
    write(6,'('' is no extended source then type "n".'')')
    write(6,'(''  '')')
    write(6,'('' make a note of the form of heating or cooling in'')')
    write(6,'('' subroutine source.  it is not echo-printed'')')
    read(5,50)ans
    if(ans == 'y' .or. ans == 'y')then
        is=1
        goto 31
    endif
    if(ans == 'n' .or. ans == 'n')then
        is=0
        goto 31
    endif
    goto 30
!  define the stefan constant
    31 sig=5.667e-8
!  input the length of the rod and the value of kappa
    ll=1.0
    kappa=400
    write(6,'('' the length of the rod is set at 1m and its '')')
    write(6,'('' thermal conductivity at 400 wm(**-1)k(**-1).'')')
    write(6,'('' do you want to change this? [y/n]'')')
    10 read(5,50)ans
    50 format(a1)
    if(ans == 'n' .or. ans == 'n')goto 11
    if(ans == 'y' .or. ans == 'y')goto 12
    goto 10
    12 write(6,'(''read in required values of the length of '')')
    write(6,'(''the rod and its conductivity.'')')
    read(5,*)ll,kappa
    11 write(6,'(''read in the number of elements in the rod'')')
    read(5,*)n
!  calculate element length
    h=ll/n
!  decide on form of output
    write(6,'(''the final temperatures can be output on the'')')
    write(6,'(''screen and/or on the printer and/or as a data'')')
    write(6,'(''file for input to a graphics package.'')')
    write(6,'(''if there are many elements the screen will not'')')
    write(6,'(''accommodate sufficient data.'')')
    write(6,'('' '')')
    60 write(6,'(''do you want screen output? [y/n]'')')
    read(5,50)ans
    if(ans == 'n' .or. ans == 'n')goto 61
    if(ans == 'y' .or. ans == 'y')goto 62
    goto 60
    61 ios=0
    goto 63
    62 ios=1
    write(6,400)ll,kappa
    400 format(15hlength of rod  ,f5.2,14h  conductivity,f8.1)
    63 write(6,'(''do you want printed output? [y/n]'')')
    read(5,50)ans
    if(ans == 'n' .or. ans == 'n')goto 64
    if(ans == 'y' .or. ans == 'y')goto 65
    goto 63
    64 iop=0
    goto 66
    65 open(unit=9,file='lpt1')
    write(9,400)ll,kappa
    iop=1
    66 write(6,'(''do you want a data file? [y/n]'')')
    read(5,50)ans
    if(ans == 'n' .or. ans == 'n')goto 67
    if(ans == 'y' .or. ans == 'y')goto 68
    goto 66
    67 iod=0
    goto 69
    68 open(unit=10,file='finel.dat')
    iod=1
    write(6,'(''please note - output file is "finel.dat"'')')
!  clear table which will hold the stiffness matrix
    69 do 1 i=1,3
        do 1 j=1,n-1
            a(i,j)=0
    1 end do
!  clear table holding right-hand-side vector
    do 21 j=1,n-1
        c(j)=0
    21 end do
!  add values of -1/h and 2/h to stiffness matrix
    do 2 j=1,n-1
        if(j == 1)goto 3
        a(1,j)=-1/h
        3 a(2,j)=2/h
        if(j == n-1)goto 2
        a(3,j)=-1/h
    2 end do
!  set the boundary conditions
    do 5 i=0,1
        mk=i*(n-2)+1
        mke=i*n
        write(6,'(''note:if there is a radiating boundary then'')')
        write(6,'(''absolute temperatures must be used'')')
        write(6,'(''  '')')
        write(6,300)i+1
        300 format(24hinput boundary condition,i2)
        write(6,'(''if at constant temperature then input temperature'')')
        write(6,'(''if insulated input "i"'')')
        write(6,'(''if exposed to an external temperature tx then '')')
        write(6,'(''input etx - e.g. e1000 [<=5 characters] '')')
        75 format(a5)
        read(5,75)cn
        if(ios == 0)goto 86
        write(6,500)i+1,cn
        write(6,'(''  '')')
        500 format(18hboundary condition,i2,4h is ,a6)
        86 if(iop == 0)goto 87
        write(9,500)i+1,cn
        87 do 6 j=1,6
            ex(j)=cn(j:j)
        6 end do
        it=0
        8 it=it+1
        k=index(stan,ex(it))-1
    !  test for illegal character
        if(k < 0)goto 900
    !  test for blank
        if(k == 10)goto 8
    !  test for insulated end
        if(k /= 11)goto 19
    !  the end is insulated. add term to stiffness matrix
        a(2,mk)=a(2,mk)-1/h
    !  indicate type of end condition
        end(i+1)=1.0e6
        goto 5
    !  test for exposed end
        19 if(k /= 12)goto 37
    !  the end is exposed
        it=it+1
    !  test for non-digit or illegal character)
        k=index(stan,ex(it))-1
        if(k > 9 .or. k < 0)goto 900
    !  the next character is a digit
        m=0
        13 m=m+1
        kk(m)=k
    !  read next digit or blank
        it=it+1
        k=index(stan,ex(it))-1
    !  test for non-digit or illegal character
        if(k > 10 .or. k < 0)goto 900
    !  test for blank which is end of number entry
        if(k < 10)goto 13
    !  blank reached so temperature now read in
        extemp=0
        do 15 l=1,m
            extemp=extemp+kk(l)*10**(m-l)
        15 end do
        cay=4*sig*extemp**3
        rat=kappa/(kappa+h*cay)
    !  add term to stiffness matrix
        a(2,mk)=a(2,mk)-rat/h
    !  add term to right hand side
        c(mk)=c(mk)+(1-rat)/h*extemp
    !  indicate the end condition
        end(i+1)=-1.0e6
        goto 5
    !  the end is at a fixed temperature
        37 m=0
        16 m=m+1
        kk(m)=k
    !  read next digit or blank
        it=it+1
        k=index(stan,ex(it))-1
    !  test for illegal character
        if(k > 10)goto 900
    !  test for blank which is end of number entry
        if(k <= 9)goto 16
    !  number entry complete
        t=0
        do 17 l=1,m
            t=t+kk(l)*10**(m-l)
        17 end do
    !  add term to right-hand-side vector
        c(mk)=c(mk)+t/h
        end(i+1)=t
    5 end do
!  check to see if there is an extended source or sink
    if(is == 0)goto 40
!  there is an extended source or sink.   the contributions to the
!  right-hand-side vector will be found by gauss 7-point quadrature.
    do 41 i=1,n-1
    !  find lower limit for quadrature
        xl=(i-1)*h
    !  carry out gauss quadrature
        sum=0
        do 42 j=1,7
            xx=xl+2*x(j)*h
            call source(xx,q)
            sum=sum+q*w(j)
        42 end do
        c(i)=c(i)+2*sum*h/kappa
    41 end do
!  extended source or sink contribution complete
!  now check for a point source or sink
    40 write(6,'(''is there a point source or sink? [y/n]'')')
    read(5,50)ans
    if(ans == 'n' .or. ans == 'n')goto 51
    if(ans == 'y' .or. ans == 'y')goto 52
    goto 40
    52 write(6,'(''read in source strength "watts per metre cubed".'')')
    write(6,'(''and the node at which it occurs'')')
    read(5,*)q,j
    if(ios == 0)goto 98
    write(6,600)q,j
    600 format(23hpoint source strength  ,e8.2,8h at node,i4)
    98 if(iop == 0)goto 99
    write(9,600)q,j
    99 c(j)=c(j)+q/kappa
!  the equations have now been set up
    51 call tridiag(a,c,y,n-1)
!  build up complete table of temperatures
    do 90 i=1,2
        m=n*(i-1)
        if(end(i) < 9e5 .and. end(i) > -9e5)then
            yy(m)=end(i)
            goto 90
        endif
        mm=(i-1)*(n-2)+1
        if(end(i) > 9e5)then
            yy(m)=y(mm)
            goto 90
        endif
    !  must be an exposed end
        yy(m)=rat*y(mm)+(1-rat)*extemp
    90 end do
    do 91 i=1,n-1
        yy(i)=y(i)
    91 end do
!  output section
    200 format(6(i4,f6.0))
    if(ios == 0)goto 70
    write(6,200)(i,yy(i),i=0,n)
    70 if(iop == 0)goto 71
    write(9,200)(i,yy(i),i=0,n)
    71 if(iod == 0)goto 1000
    do 80 i=0,n
        xc(i)=i*h
    80 end do
    do 81 i=0,n
        write(10,*)xc(i),yy(i)
    81 end do
    close(10)
    goto 1000
    900 write(6,'('' illegal character - program aborted'')')
    1000 stop
    end program

    subroutine source(x,q)
    q=0
    return
    end subroutine source


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
