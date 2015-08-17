    program mcplate
!**********************************************************************
!  this uses the monte carlo method to solve problems involving the
!  thermal equilibrium in a plate with insulated top and bottom surfaces
!  and boundaries which are at a fixed temperature.

!  the plate is defined on a square mesh contained within a
!  rectangular region of nr rows and nc columns (i.e.(nr-1) x (nc-1)
!  square elements) with both nr and nc <= 10.  all sides of the
!  plate must either be along the principal directions of the mesh
!  or at 45 degrees to them.  each side of the plate must contain
!  at least three points including the end points.

!  when initial conditions are set up points outside the plate are
!  entered as x, internal points, the temperatures of which have to
!  be determined, are entered as u  and points at the boundary as
!  as an integer <= 10,000).

!  there is provision for the plate to be heated (or cooled), which
!  requires the amendment of the function subprogram heat(x,y) which
!  gives the rate of heating per unit volume as a function of x and y.
!  in addition the data statement must be changed to give the length
!  of the side of the square grid, h, and the conductivity of the
!  plate, cappa.   the coordinates (x,y) are referred to the point
!  (i,j) = (1,1) as origin

!***********************************************************************
    dimension ip(0:12,0:12),q(0:12,0:12),kk(5),sd(0:12,0:12)
    real :: mcsum
    character ex(60)*1,stan*13,cn*60,ans*1
    open(unit=9,file='lpt1')
    stan='0123456789 xu'
    data h,cappa/0.125,400.0/
! data for the random number generator.   ir is the seed
    data ir,ix,iy,im/199,171,11213,53125/
    sig=5.67e-8
    pi=4.0*atan(1.0)
    qpf=h*h/4.0/cappa
! initially indicate all points as outside plate and clear table sd
    do 55 i=0,12
        do 55 j=0,12
            ip(i,j)=1
            q(i,j)=-1
            sd(i,j)=0
    55 end do
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
    97 write(6,'(''are you sure that function subprogram "heat"'')')
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
        !  the character must be a digit
            m=0
            4 m=m+1
            kk(m)=k
        !  read next digit or blank
            it=it+1
            k=index(stan,ex(it))-1
        !  test for illegal character
            if(k < 0)goto 999
        !  test if x or u is illegally combined with a digit
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
!  nmc monte carlo random walks are made for each point for which the
!  temperature must be determined.
    write(6,'('' '')')
    write(6,'(''input number of monte carlo trials per point'')')
    read(5,*)nmc
    do 11 i=1,nr
        do 11 j=1,nc
            if(ip(i,j) /= 0)goto 11
        !  an interior point has been identified. random walks begin.
            tot=0
            tot2=0
            do 12 k=1,nmc
                ii=i
                jj=j
                x=(i-1)*h
                y=(j-1)*h
                mcsum=qpf*heat(x,y)
                86 ir=mod(ir*ix+iy,im)
                rr=float(ir)/float(im)
                nrr=int(4.0*rr)+1
                goto(80,81,82,83)nrr
                80 ii=ii+1
                goto 84
                81 ii=ii-1
                goto 84
                82 jj=jj+1
                goto 84
                83 jj=jj-1
                84 if(ip(ii,jj) == 3)goto 85
            !  the next point is not a boundary point
                x=(ii-1)*h
                y=(jj-1)*h
                mcsum=mcsum+qpf*heat(x,y)
                goto 86
            !  a boundary point has been reached
                85 mcsum=mcsum+q(ii,jj)
                tot=tot+mcsum
                tot2=tot2+mcsum*mcsum
            12 end do
        !  the random walks are finished for the point (i, j).
            q(i,j)=tot/nmc
        !  calculate standard deviation
            sd(i,j)=sqrt(tot2/nmc-q(i,j)**2)/sqrt(float(nmc))
    11 end do
!  output results
    write(iout,'('' '')')
    do 15 i=1,nr
        write(iout,300)(q(i,j),j=1,nc)
        300 format(11f6.0)
        400 format(11f6.1)
        write(iout,'('' '')')
        write(iout,'('' '')')
    15 end do
    do 16 i=1,nr
        write(iout,400)(sd(i,j),j=1,nc)
        write(iout,'('' '')')
        write(iout,'('' '')')
    16 end do
    goto 998
    999 write(6,'('' illegal character - program terminated'')')
    998 stop
    end program
       
    function heat(x,y)
    heat=0
    end function heat

