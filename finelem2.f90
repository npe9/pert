    program finelem
!  this program gives a finite-element solution of the problem of a
!  heated plate insulated on its top and bottom surfaces with all
!  elements triangular.   differential boundary conditions are
!  specified by
!           cappa*(d(theta)/d(n)) = -m*theta + s
!  where n is along the direction of the normal to the surface.
!  values of m and s are input and so are estimated before the
!  program is used.  some nodes, usually at boundaries, can be
!  specified as having fixed temperatures.   there is provision for
!  an extended source of heat which is given by a function statement
!  q(x,y).   if such a source is present then the user must ensure
!  that it gives the required form of heating.   there is also
!  provision for up to 9 point sources which must be situated at nodes.

    dimension nde(100,3),x(100),y(100),area(100),a(100,3),b(100,3), &
    c(100,3),nf(50),tf(50),ietn(50),nte(50),nedge(50,2),cay(50), &
    ss(50),np(50),sps(50),rhs(100),stiff(100,100),side(50), &
    ntemp(50),ttemp(50)
    character ans*1
    ichange=0
    545 write(6,'('' '')')
    write(6,'('' warning warning warning warning warning'')')
    write(6,'('' '')')
    write(6,'(''if there is an extended heat source you must'')')
    write(6,'(''confirm that your extended heat source function'')')
    write(6,'(''statement is correct [y/n]. if not then abort the'')')
    write(6,'(''program, correct routine and re-compile'')')
    write(6,'(''if there is no extended heat source then it does'')')
    write(6,'(''not matter what the function statement contains.'')')
    read(5,50)ans
    if(ans == 'y' .or. ans == 'y')goto 547
    if(ans == 'n' .or. ans == 'n')goto 544
    goto 545
!  input conductivity of the plate
    547 inp=5
    546 write(6,'(''data can be input either from the keyboard or'')')
    write(6,'(''from a data file "finelem.dat" prepared by the'')')
    write(6,'(''program "findata".  do you want input via a data'')')
    write(6,'(''file? [y/n]'')')
    read(5,50)ans
    if(ans == 'n' .or. ans == 'n')goto 543
    if(ans == 'y' .or. ans == 'y')then
        open(unit=11,file='finelem.dat')
        inp=11
        goto 543
    endif
    goto 546
    543 if(inp == 11)goto 444
    write(6,'(''read in conductivity'')')
    444 read(inp,*)cappa
    iout=0
    121 write(6,'(''do you want printed output? [y/n]'')')
    read(5,50)ans
    if(ans == 'n' .or. ans == 'n')goto 120
    if(ans == 'y' .or. ans == 'y')then
        iout=1
        open(unit=9,file='lpt1')
        goto 120
    endif
    goto 121
    120 if(inp == 11)goto 445
    write(6,'(''input number of elements and number of nodes'')')
    445 read(inp,*)ne,nn
    do 1 i=1,ne
        if(inp == 11)goto 461
        write(6,100)i
        100 format(38h enter 3 node numbers defining element,i2)
        461 read(inp,*)(nde(i,j),j=1,3)
    1 end do
    if(inp == 11)goto 446
    write(6,'(''enter coordinates of nodes'')')
    446 do 2 i=1,nn
        if(inp == 11)goto 448
        write(6,200)i,i
        200 format(9h enter x[,i2,8h] and y[,i2,1h])
        448 read(inp,*)x(i),y(i)
    2 end do
!  the area of each element is calculated as are the values of
!  a, b and c [equation(6.66b)] for each node.
    497 do 3 i=1,ne
        do 3 j=1,3
            m1=nde(i,j)
            m2=nde(i,1)
            m3=nde(i,2)
            if(j == 1)then
                m2=nde(i,2)
                m3=nde(i,3)
            endif
            if(j == 2)then
                m2=nde(i,3)
                m3=nde(i,1)
            endif
            a(i,j)=x(m2)*y(m3)-x(m3)*y(m2)
            b(i,j)=y(m2)-y(m3)
            c(i,j)=x(m3)-x(m2)
            area(i)=abs((a(i,j)+x(m1)*b(i,j)+y(m1)*c(i,j))/2.0)
    3 end do
    if(ichange == 1)goto 498
!  clear arrays containing information about fixed-temperature nodes
    do 6 i=1,50
        nf(i)=0
        tf(i)=0
    6 end do
    if(inp == 11)goto 449
    write(6,'(''enter number of nodes with fixed temperatures'')')
    449 read(inp,*)nft
    if(nft == 0)goto 5
    do 4 i=1,nft
        if(inp == 11)goto 460
        write(6,'(''read in node number, fixed temperature'')')
        460 read(inp,*)k,t
        nf(k)=1
        tf(k)=t
    4 end do
!  since nodes with fixed temperatures do not lead to equations node
!  numbers are now linked to equation numbers and vice-versa.
    5 neq=0
    do 7 i=1,nn
        if(nf(i) == 1)goto 7
        neq=neq+1
        ietn(neq)=i
    7 end do
!  neq is now the number of equations to be generated.   the kth
!  equation is linked with variational parameter theta(ietn(k))

!  now clear table nte
    do 8 i=1,nn
        nte(i)=0
    8 end do
!  now generate table nte such that nte(k)=0 if node k is at constant
!  temperature node otherwise nte(k) gives the equation associated
!  with node k.
    do 9 i=1,neq
        k=ietn(i)
        nte(k)=i
    9 end do
    12 write(6,'(''is there an extended heat source [y/n] ?'')')
    iq=0
    read(5,50)ans
    50 format(a1)
    if(ans == 'n' .or. ans == 'n')goto 10
    if(ans == 'y' .or. ans == 'y')goto 11
    goto 12
    11 iq=1
    10 if(inp == 11)goto 451
    write(6,'(''enter the number of point sources'')')
    451 read(inp,*)ip
    if(ip == 0)goto 13
    do 14 i=1,ip
        if(inp == 11)goto 452
        write(6,'(''give node and strength of point source [w/m**3]'')')
        452 read(inp,*)np(i),sps(i)
    14 end do
!  now boundary conditions are to be entered. the edges are identified
!  by the terminating nodes of each element side.   the values of m
!  and s are required for each edge.
    13 if(inp == 11)goto 453
    write(6,'(''how many edges of elements have boundary '')')
    write(6,'(''conditions of the form '')')
    write(6,'('' d(theta)/d(n) = -m*theta + s ?'')')
    453 read(inp,*)nbc
    if(nbc == 0)goto 16
    do 15 i=1,nbc
        if(inp == 11)goto 454
        write(6,'(''read in the edge as the two end nodes and the '')')
        write(6,'(''values of m and s'')')
        454 read(inp,*)nedge(i,1),nedge(i,2),cay(i),ss(i)
    ! calculate the length of the edge
        k1=nedge(i,1)
        k2=nedge(i,2)
        side(i)=sqrt((x(k1)-x(k2))**2+(y(k1)-y(k2))**2)
    15 end do
    write(6,'(''there is now an opportunity to examine most of'')')
    write(6,'(''data which has been entered and to correct it'')')
    write(6,'(''if necessary. this is advisable if the entered'')')
    write(6,'(''data are extensive.'')')
    write(6,'(''  '')')
    107 write(6,'(''do you want to check the node numbers defining'')')
    write(6,'(''elements with a view to correction? [y/n]'')')
    read(5,50)ans
    if(ans == 'n' .or. ans == 'n')goto 117
    if(ans == 'y' .or. ans == 'y')goto 106
    goto 107
    106 write(6,450)(i,(nde(i,j),j=1,3),i=1,ne)
    450 format(3(i7,i6,2i3))
    112 write(6,'(''  '')')
    write(6,'(''do you want to change any of these? [y/n]'')')
    read(5,50)ans
    if(ans == 'n' .or. ans == 'n')goto 117
    if(ans == 'y' .or. ans == 'y')goto 111
    goto 112
    111 ichange=1
    write(6,'(''  '')')
    write(6,'(''enter element number + three node numbers'')')
    read(5,*)nel,(nde(nel,j),j=1,3)
    goto 107
    117 write(6,'(''  '')')
    write(6,'(''do you want to check the coordinates of the'')')
    write(6,'(''nodes with a view to correction? [y/n]'')')
    read(5,50)ans
    if(ans == 'n' .or. ans == 'n')goto 137
    if(ans == 'y' .or. ans == 'y')goto 116
    goto 117
    116 write(6,550)(i,x(i),y(i),i=1,nn)
    550 format(3(i7,2f8.4))
    122 write(6,'(''  '')')
    write(6,'(''do you want to change any of these? [y/n]'')')
    read(5,50)ans
    if(ans == 'n' .or. ans == 'n')goto 137
    if(ans == 'y' .or. ans == 'y')goto 131
    goto 122
    131 ichange=1
    write(6,'(''  '')')
    write(6,'(''enter node number + x and y coordinates'')')
    read(5,*)nod,x(nod),y(nod)
    goto 117
    137 if(ichange == 1)goto 497
    498 if(nft == 0)goto 135
    write(6,'(''  '')')
    write(6,'(''do you want to check the nodes with fixed'')')
    write(6,'(''temperatures with a view to correction? [y/n]'')')
    read(5,50)ans
    if(ans == 'n' .or. ans == 'n')goto 135
    if(ans == 'y' .or. ans == 'y')goto 136
    goto 137
    136 write(6,650)(i,tf(i),i=1,nn)
    650 format(4(i7,f6.1))
    write(6,'(''  '')')
    write(6,'(''if the second column number is zero then'')')
    write(6,'(''the temperature is not fixed'')')
    132 write(6,'(''  '')')
    write(6,'(''do you want to change any of these? [y/n]'')')
    read(5,50)ans
    if(ans == 'n' .or. ans == 'n')goto 135
    if(ans == 'y' .or. ans == 'y')goto 141
    goto 132
    141 write(6,'(''  '')')
    write(6,'(''enter, node number and temperature.  if you '')')
    write(6,'(''wish to cancel a fixed temperature allocation'')')
    write(6,'(''then enter minus node number plus zero.'')')
    read(5,*)i,tt
    if(i < 0)then
        i=-i
        nf(i)=0
        tf(i)=0
        goto 137
    endif
    tf(i)=tt
    nf(i)=1
    goto 137
    135 if(ip == 0)goto 155
    write(6,'(''  '')')
    write(6,'(''do you want to check the nodes which are point'')')
    write(6,'(''sources with a view to correction? [y/n]'')')
    read(5,50)ans
    if(ans == 'n' .or. ans == 'n')goto 155
    if(ans == 'y' .or. ans == 'y')goto 156
    goto 135
    156 do 154 i=1,ip
        write(6,750)np(i),sps(i)
        750 format(i4,e9.4)
    154 end do
    152 write(6,'(''  '')')
    write(6,'(''do you want to change any of these? [y/n]'')')
    read(5,50)ans
    if(ans == 'n' .or. ans == 'n')goto 155
    if(ans == 'y' .or. ans == 'y')goto 161
    goto 152
    161 write(6,'(''  '')')
    write(6,'(''enter list position, node and strength of'')')
    write(6,'(''source.'')')
    read(5,*)i,np(i),sps(i)
    goto 135
    155 if(nbc == 0)goto 16
    write(6,'(''  '')')
    write(6,'(''do you want to check the edges corresponding'')')
    write(6,'(''differential boundary conditions with a view'')')
    write(6,'(''to correction? [y/n]'')')
    read(5,50)ans
    if(ans == 'n' .or. ans == 'n')goto 16
    if(ans == 'y' .or. ans == 'y')goto 176
    goto 155
    176 do 174 i=1,nbc
        write(6,850)nedge(i,1),nedge(i,2),cay(i),ss(i)
        850 format(2i4,2e12.4)
    174 end do
    172 write(6,'(''  '')')
    write(6,'(''do you want to change any of these? [y/n]'')')
    read(5,50)ans
    if(ans == 'n' .or. ans == 'n')goto 16
    if(ans == 'y' .or. ans == 'y')goto 181
    goto 172
    181 write(6,'(''  '')')
    write(6,'(''enter list position, nodes defining edge and'')')
    write(6,'(''values of m and s.'')')
    read(5,*)i,nedge(i,1),nedge(i,2),cay(i),ss(i)
    goto 155
!  the input is now complete and the process of building up the stiffness
!  matrix and the right-hand-side vector can commence.

!  before that, if printing has been specified, the input data will be
!  printed.
    16 if(iout /= 1)goto 17
    write(9,870)cappa
    870 format(17h conductivity is ,f6.1,18h w m**[-1] k**[-1])
    write(9,'(''  '')')
    write(9,'(''node numbers associated with each element'')')
    write(9,960)(i,(nde(i,j),j=1,3),i=1,ne)
    960 format(3(i7,i6,2i4))
    write(9,'(''  '')')
    write(9,'(''coordinates of nodes'')')
    write(9,950)(i,x(i),y(i),i=1,nn)
    950 format(3(i6,f7.4,f7.4))
    if(nft == 0)goto 830
    write(9,'(''  '')')
    write(9,'(''nodes with fixed temperatures'')')
    ny=0
    do 760 i=1,nn
        if(nf(i) == 0)goto 760
        ny=ny+1
        ntemp(ny)=i
        ttemp(ny)=tf(i)
    760 end do
    write(9,970)(ntemp(i),ttemp(i),i=1,ny)
    970 format(4(i6,f7.1))
    write(9,'(''  '')')
    830 if(ip == 0)goto 840
    write(9,'(''positions and strengths of point sources'')')
    do 745 i=1,ip
        write(9,980)np(i),sps(i)
    745 end do
    980 format(i6,e13.4)
    write(9,'(''  '')')
    840 if(nbc == 0)goto 17
    write(9,'(''differential boundary edges and parameters'')')
    write(9,'('' node1 node2     m        s'')')
    write(9,875)(nedge(i,1),nedge(i,2),cay(i),ss(i),i=1,nbc)
    875 format(i5,i6,f8.3,f10.3)
    write(9,'(''  '')')
!  clear arrays
!  divide the values of m and s by cappa
    17 do 115 i=1,nbc
        cay(i)=cay(i)/cappa
        ss(i)=ss(i)/cappa
    115 end do
    do 19 i=1,100
        rhs(i)=0
        do 19 j=1,100
            stiff(i,j)=0
    19 end do
!  first deal with element derivatives.
    do 20 i=1,ne
        do 20 j=1,3
            k=nde(i,j)
            if(nf(k) == 1)goto 20
        !  node k corresponds to a variational parameter. find the equation
        !  it generates.
            ke=nte(k)
            stiff(ke,ke)=stiff(ke,ke)+0.25*(b(i,j)**2+c(i,j)**2)/area(i)
            do 21 l=1,2
                jj=mod(j+l,3)
                if(jj == 0)jj=3
                k=nde(i,jj)
                if(nf(k) == 1)then
                    z=0.25*tf(k)*(b(i,j)*b(i,jj)+c(i,j)*c(i,jj))/area(i)
                    rhs(ke)=rhs(ke)-z
                    goto 21
                endif
                kf=nte(k)
                z=0.25*(b(i,j)*b(i,jj)+c(i,j)*c(i,jj))/area(i)
                stiff(ke,kf)=stiff(ke,kf)+z
            21 end do
    20 end do
!  the extended heating will be added as though it is constant in
!  each element with a value equal to the average at the three nodes
    if(iq == 0)goto 30
    do 23 i=1,ne
        i1=nde(i,1)
        i2=nde(i,2)
        i3=nde(i,3)
        qav=(q(x(i1),y(i1))+q(x(i2),y(i2))+q(x(i3),y(i3)))/3.0
        z=qav*area(i)/3.0/cappa
        if(nf(i1) == 1)goto 24
        k=nte(i1)
        rhs(k)=rhs(k)+z
        24 if(nf(i2) == 1)goto 25
        k=nte(i2)
        rhs(k)=rhs(k)+z
        25 if(nf(i3) == 1)goto 23
        k=nte(i3)
        rhs(k)=rhs(k)+z
    23 end do
! add contributions of differential boundary conditions
    30 if(nbc == 0)goto 60
    do 31 i=1,nbc
        do 31 l=1,2
            if(nf(nedge(i,l)) == 1)goto 31
            ke=nte(nedge(i,l))
            stiff(ke,ke)=stiff(ke,ke)-cay(i)*side(i)/3.0
            if(nf(nedge(i,3-l)) == 0)then
                kf=nte(nedge(i,3-l))
                stiff(ke,kf)=stiff(ke,kf)-cay(i)*side(i)/6.0
                goto 31
            endif
            rhs(ke)=rhs(ke)+cay(i)*side(i)*tf(nedge(i,3-l))/6.0
    31 end do
    do 41 i=1,nbc
        do 41 l=1,2
            if(nf(nedge(i,l)) == 1)goto 41
            ke=nte(nedge(i,l))
            rhs(ke)=rhs(ke)+ss(i)*side(i)/2.0
    41 end do
!  now add point sources
    60 if(ip == 0)goto 61
    do 62 i=1,ip
        ke=nte(np(i))
        rhs(ke)=rhs(ke)+sps(i)/cappa
    62 end do
!  at this point the equations have been set up and must now be solved
    61 call gaussj(stiff,neq,100,rhs,1,1)
!  assign temperatures to all nodes
    ns=0
    do 101 i=1,nn
        if(nf(i) == 1)then
            sps(i)=tf(i)
            goto 101
        endif
        ns=ns+1
        sps(i)=rhs(ns)
    101 end do
    nout=6
    102 write(nout,300)
    300 format(6h node ,6h  temp)
    write(nout,400)(i,sps(i),i=1,nn)
    400 format(4(i4,f8.0))
    if(iout == 1)then
        nout=9
        iout=0
        goto 102
    endif
    544 stop
    end program

    function q(x,y)
    q=1.0e3
    end function q

    subroutine gaussj(a,n,np,b,m,mp)
    parameter (nmax=50)
    dimension a(np,np),b(np,mp),ipiv(nmax),indxr(nmax),indxc(nmax)
    do 11 j=1,n
        ipiv(j)=0
    11 end do
    do 22 i=1,n
        big=0.
        do 13 j=1,n
            if(ipiv(j) /= 1)then
                do 12 k=1,n
                    if (ipiv(k) == 0) then
                        if (abs(a(j,k)) >= big)then
                            big=abs(a(j,k))
                            irow=j
                            icol=k
                        endif
                    else if (ipiv(k) > 1) then
                        pause 'singular matrix'
                    endif
                12 end do
            endif
        13 end do
        ipiv(icol)=ipiv(icol)+1
        if (irow /= icol) then
            do 14 l=1,n
                dum=a(irow,l)
                a(irow,l)=a(icol,l)
                a(icol,l)=dum
            14 end do
            do 15 l=1,m
                dum=b(irow,l)
                b(irow,l)=b(icol,l)
                b(icol,l)=dum
            15 end do
        endif
        indxr(i)=irow
        indxc(i)=icol
        if (abs(a(icol,icol)) < 1.0e-30) pause 'singular matrix.'
        pivinv=1./a(icol,icol)
        a(icol,icol)=1.
        do 16 l=1,n
            a(icol,l)=a(icol,l)*pivinv
        16 end do
        do 17 l=1,m
            b(icol,l)=b(icol,l)*pivinv
        17 end do
        do 21 ll=1,n
            if(ll /= icol)then
                dum=a(ll,icol)
                a(ll,icol)=0.
                do 18 l=1,n
                    a(ll,l)=a(ll,l)-a(icol,l)*dum
                18 end do
                do 19 l=1,m
                    b(ll,l)=b(ll,l)-b(icol,l)*dum
                19 end do
            endif
        21 end do
    22 end do
    do 24 l=n,1,-1
        if(indxr(l) /= indxc(l))then
            do 23 k=1,n
                dum=a(k,indxr(l))
                a(k,indxr(l))=a(k,indxc(l))
                a(k,indxc(l))=dum
            23 end do
        endif
    24 end do
    return
    end subroutine gaussj
