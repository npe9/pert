    program hotplate
!**********************************************************************
!  this solves problems involving the thermal equilibrium in a
!  plate with insulated top and bottom surfaces and boundaries
!  which are at a fixed temperature,insulated or exchanging
!  heat by convection or radiation with the outside.

!  the plate is defined on a square mesh contained within a
!  rectangular region of nr rows and nc columns (i.e.(nr-1) x (nc-1)
!  square elements) with both nr and nc <= 10.  all sides of the
!  plate must either be along the principal directions of the mesh
!  or at 45 degrees to them.  each side of the plate must contain
!  at least three points including the end points.

!  when initial conditions are set up points outside the plate are
!  entered as x, internal points, the temperatures of which have to
!  be determined, are entered as u  and points at the boundary as
!      (1) a fixed temperature (as an integer <= 10,000)
!   or (2) as i, if insulated
!   or (3) as e, if exchanging heat with the outside.  the boundary
!          condition is described by
!            kappa*(d(theta)/dx) = -k*theta + s
!          and values of k and s will be input for each such point.

!  there is provision for the plate to be heated (or cooled), which
!  requires the amendment of the function subprogram heat(x,y) which
!  gives the rate of heating per unit volume as a function of x and y.
!  in addition the data statement must be changed to give the length
!  of the side of the square grid, h, and the conductivity of the
!  plate, cappa.   the coordinates (x,y) are referred to the point
!  (i,j) = (1,1) as origin

!  the method of solution is gauss-seidel with sequential over-relaxation.
!***********************************************************************
    dimension ip(0:12,0:12),q(0:12,0:12),ij(0:12,0:12),kk(5), &
    ji(0:12,0:12),ignore(-1:1,-1:1),extemp(11,11),d(0:12,0:12), &
    bound(11,11,2)
    character ex(60)*1,stan*15,cn*60,ans*1
    open(unit=9,file='lpt1')
    stan='0123456789 xuie'
    data h,cappa/0.125,400.0/
! initially indicate all points as outside plate
    do 55 i=0,12
        do 55 j=0,12
            ip(i,j)=1
            q(i,j)=-1
    55 end do
! initialize array extemp
    iout=6
    95 write(6,'(''do you want printed output? [y/n]'')')
    read(5,750)ans
    750 format(a1)
    if(ans == 'n' .or. ans == 'n')goto 97
    if(ans == 'y' .or. ans == 'y')then
        iout=9
        goto 97
    endif
    goto 95
    97 do 56 i=1,11
        do 56 j=1,11
            extemp(i,j)=1.0e35
    56 end do
    write(6,'(''are you sure that function subprogram "heat"'')')
    write(6,'(''and values of "h" and "cappa" are correct '')')
    write(6,'(''for this "hotplate" application? if so then'')')
    write(6,'(''enter number of rows and columns in mesh'')')
    read(5,*)nr,nc
!  the characteristics of the plate are now entered
    do 1 i=1,nr
        write(6,'('' '')')
        write(6,50)i
        50 format(31h enter information for grid row,i3)
        write(6,'(''x if point is outside plate'')')
        write(6,'(''temperature on the boundary [integer<=10000]'')')
        write(6,'(''i if insulated on boundary'')')
        write(6,'(''e if a boundary point exchanging heat with the '')')
        write(6,'(''outside.'')')
        write(6,'(''u if within the plate and to be determined'')')
        read(5,100)cn
        100 format(a60)
        do 60 nx=1,60
            ex(nx)=cn(nx:nx)
        60 end do
        it=0
        do 2 j=1,nc
            3 it=it+1
            k=index(stan,ex(it))-1
        !  test of illegal character
            if(k < 0)goto 999
        !  test for blank
            if(k == 10)goto 3
        !  test for x - point outside plate
            if(k == 11)then
                ip(i,j)=1
                q(i,j)=-1
                goto 2
            endif
        !  test for u - interior point of plate
            if(k == 12)then
                ip(i,j)=0
                q(i,j)=0
                goto 2
            endif
        !  test for i - insulated boundary point
            if(k == 13)then
                ip(i,j)=2
                q(i,j)=0
                goto 2
            endif
        !  test for point exchanging heat with the outside
            if(k == 14)then
                ip(i,j)=4
                q(i,j)=0
                goto 2
            endif
        !  test for non-digit or illegal character
            if(k > 9)goto 999
            if(k < 0)goto 999
        !  the character must be a digit
            m=0
            4 m=m+1
            kk(m)=k
        !  read next digit or blank
            it=it+1
            k=index(stan,ex(it))-1
        !  test for illegal character
            if(k < 0)goto 999
        !  test if x, u or e is illegally combined with a digit
            if(k > 10)goto 999
        !  test for blank which is end of number entry
            if(k <= 9)goto 4
            ip(i,j)=3
            q(i,j)=0
            do 6 l=1,m
                q(i,j)=q(i,j)+kk(l)*10**(m-l)
            6 end do
        2 end do
    1 end do
!  read in values of k and s for the differential boundary points
    do 64 i=1,nr
        do 64 j=1,nc
            if(ip(i,j) /= 4)goto 64
            write(6,700)i,j
            700 format(38h read in values of k and s for point [,i2,1h ,i2,1h])
            read(5,*)bound(i,j,1),bound(i,j,2)
    64 end do
!  set up equivalent-temperature points for the 'ghost' points
!  required to establish required slope at insulated boundary.
    do 30 i=1,nr
        do 30 j=1,nc
            if((ip(i,j)-2)*(ip(i,j)-4) /= 0)goto 30
        !  an insulated or heat-exchanging boundary point has been identified.
        !  look for neighbouring boundary points.
            do 29 ii=-1,1
                do 29 jj=-1,1
                    ignore(ii,jj)=0
            29 end do
            do 31 ii=-1,1
                do 31,jj=-1,1
                    if(ii == 0 .and. jj == 0)goto 31
                    if(ignore(ii,jj) == 1)goto 31
                    if(ip(i+ii,j+jj) < 2)goto 31
                !  test that points are neighbours on the same side
                    if(ip(i+2*ii,j+2*jj) < 2 .and. ip(i-ii,j-jj) < 2)goto 31
                !  a neighbouring boundary point has been identified
                !  eliminate consideration of point on other side of (i,j) which will
                !  give the same 'ghost' point
                    ignore(-ii,-jj)=0
                !  find the 'ghost' point(s) outside the plate
                    if(ii*jj == 0)then
                    !  the boundary is parallel to a principal direction
                        d(i,j)=2.0
                        ia=iabs(ii)
                        ja=iabs(jj)
                        if(ip(i+ja,j+ia) == 1)goto 33
                        ij(i-ja,j-ia)=i+ja
                        ji(i-ja,j-ia)=j+ia
                        d(i-ja,j-ia)=2.0
                        goto 31
                        33 ij(i+ja,j+ia)=i-ja
                        ji(i+ja,j+ia)=j-ia
                        d(i+ja,j+ia)=2.0
                        goto 31
                    endif
                !  the boundary is at 45 degrees to principal directions
                    d(i,j)=sqrt(2.0)
                    ixj=ii*jj
                    if(ip(i+1,j) /= 1)goto 35
                    ij(i+1,j)=i
                    ji(i+1,j)=j+ixj
                    d(i+1,j)=sqrt(2.0)
                    goto 36
                    35 if(ip(i,j+ixj) /= 1)goto 36
                    ij(i,j+ixj)=i+1
                    ji(i,j+ixj)=j
                    d(i,j+ixj)=sqrt(2.0)
                    36 if(ip(i-1,j) /= 1)goto 37
                    ij(i-1,j)=i
                    ji(i-1,j)=j-ixj
                    d(i-1,j)=sqrt(2.0)
                    goto 31
                    37 if(ip(i,j-ixj) /= 1)goto 31
                    ij(i,j-ixj)=i-1
                    ji(i,j-ixj)=j
                    d(i,j-ixj)=sqrt(2.0)
            31 end do
    30 end do
    write(6,'(''input tolerance - 0.01 will usually suffice'')')
    read(5,*)tol
    write(6,'(''input over-relaxation factor [ 1 to 2 ]'')')
    read(5,*)w
!  the solution iteration now commences. maximum number of cycles = 1000
    icycle=0
    12 icycle=icycle+1
    if(icycle > 1000)goto 20
    dif=0
    do 11 i=1,nr
        do 11 j=1,nc
            if(ip(i,j)*(ip(i,j)-2)*(ip(i,j)-4) /= 0)goto 11
        !  a point whose temperature is to be determined has been identified.
            x=(i-1)*h
            y=(j-1)*h
            est=heat(x,y)*h*h/cappa/4.0
            do 66 in=-1,1
                do 66 jn=-1,1
                    if(in == 0 .and. jn == 0)goto 66
                    if(in*jn /= 0)goto 66
                !  the neighbouring point being taken is along a principal direction.
                    if(ip(i+in,j+jn) == 1)then
                    !  the neighbouring point is outside the plate and the temperature
                    !  is found to give the correct boundary conditions.
                        im=ij(i+in,j+jn)
                        jm=ji(i+in,j+jn)
                        est=est+0.25*q(im,jm)
                        add=0
                        if(extemp(i,j) < 1.0e34)then
                        !  the surface is exchanging heat with the surroundings.  the
                        !  temperature of the false point is found from the gradient
                        !  first find the temperature half way between false and interior
                        !  point.
                            if(mod(im+i+in,2) == 1)then
                            !  the half-way point is in the middle of an element
                                am=float(im+i+in)/2.0
                                bm=float(jm+j+jn)/2.0
                                i1=am+0.5*in*jn+0.1
                                i2=am-0.5*in*jn+0.1
                                j1=bm+0.5*in*jn+0.1
                                j2=bm-0.5*in*jn+0.1
                                temp=0.5*(q(i1,j1)+q(i2,j2))
                                goto 87
                            endif
                            temp=q((im+i+in)/2,(jm+j+jn)/2)
                            87 add=(-bound(i,j,1)*temp+bound(i,j,2))*h*d(i+in,j+jn)/cappa
                        endif
                        est=est+add
                        goto 66
                    endif
                    est=est+0.25*q(i+in,j+jn)
            66 end do
            z=w*(est-q(i,j))
            q(i,j)=q(i,j)+z
            if(abs(z) > dif)dif=abs(z)
    11 end do
    if(dif > tol)goto 12
    20 write(iout,500)tol
    write(iout,550)w
    write(iout,600)icycle
    500 format(17h the tolerance is,f9.5)
    550 format(33h the over-relaxation parameter is,f6.2)
    600 format(24h the number of cycles is,i5)
    write(iout,'('' '')')
    do 15 i=1,nr
        write(iout,300)(q(i,j),j=1,nc)
        300 format(11f6.0)
        write(iout,'('' '')')
        write(iout,'('' '')')
    15 end do
    goto 998
    999 write(6,'('' illegal character - program terminated'')')
    998 stop
    end program
       
    function heat(x,y)
    heat=0
    end function heat

