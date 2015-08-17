    program pic
!  a simple p.i.c. code for calculating the langmuir sheath
!  at the edge of a plasma due to electrons only: the ions are
!  stationary.
!  nt 'electron' particles (<=1000) are used to model a one
!  dimensional electron distribution of density dens on a mesh
!  of it cells and spatial extent xt debye lengths.
!  the plasma has ion density dens over x1 debye lengths, a
!  linear density fall-off over x2 debye lengths and is zero
!  over the remainder of the mesh.
!  the electrons are assigned a maxwell-boltzmann velocity
!  distribution of temperature temp in one dimension.
!  the time step is limited by cell size and thermal speed, and
!  by the plasma frequency. for this problem the former is
!  generally dominant.
!  the positions and velocities of the electrons are output every
!  60 timesteps and there are 9 outputs including the initial
!  configuration.   these are in files v10.dat to v18.dat and
!  e10.dat to e18.dat.   the v-files give the positions and
!  velocities of the electrons and the e-files give the position
!  in the cell and the overall field.
! ******************************************************************
!  note - this program requires long integers
! ******************************************************************
    common /blkpos/ x(1000)
    common /blkvel/ v(1000)
    common /blkfld/ e(1000)
    common /blkden/ rho(1000)
    common /blkion/ rhop(1000)
    common /blkcel/ dx,xt,x1,x2
    common /blkdat/ dt,emt,ecf,dens,temp,vtemp
    common /blktim/ time
    common /blkint/ it,nt
    nout=59
    call input
    write(6,'(''output times'')')
    call set
    call field
    1 call push
    call field
    nout=nout+1
    if((nout/60)*60 == nout)then
        call output(nout)
        write(6,*)time
    endif
    time=time+dt
    if (nout < 540) go to 1
    stop
    end program



    subroutine input
    character ans*1
    common /blkcel/ dx,xt,x1,x2
    common /blkdat/ dt,emt,ecf,dens,temp,vtemp
    common /blktim/ time
    common /blkint/ it,nt
    xt=10
    it=100
    nt=1000
    dt=0.2
    x1=2.5
    x2=0
    dens=1e25
    temp=1e4
    13 write(6,'(''the following data are provided'')')
    write(6,101)xt
    101 format(10h [1] xt = ,f8.3,30h - the extent in debye lengths)
    write(6,102)it
    102 format(10h [2] it = ,i4,28h - the number of cells in xt)
    write(6,103)nt
    103 format(10h [3] nt = ,i5,26h - the number of electrons)
    write(6,104)dt
    104 format(10h [4] dt = ,f5.2,31h - timestep control from plasma)
    write(6,'(''    frequency. see comment in subroutine step.'')')
    write(6,105)x1
    105 format(10h [5] x1 = ,f6.3,35h - the number of debye lengths with)
    write(6,'(''    density dens.'')')
    write(6,106)x2
    106 format(10h [6] x2 = ,f6.3,36h - the number of debye lengths after)
    write(6,'(''    x1 where density falls linearly to zero.'')')
    write(6,107)dens
    107 format(12h [7] dens = ,e8.3,28h - the density of electrons.)
    write(6,108)temp
    108 format(12h [8] temp = ,e8.3,28h - the electron temperature.)
    write(6,'(''do you want to change any of these? [y/n]'')')
    read(5,50)ans
    50 format(a1)
    if(ans == 'n' .or. ans == 'n')goto 11
    if(ans == 'y' .or. ans == 'y')goto 12
    goto 13
    12 write(6,'(''give number of item you wish to change.'')')
    read(5,*)nitem
    goto(1,2,3,4,5,6,7,8)nitem
    1 write(6,'(''read in new value of xt.'')')
    read(5,*)xt
    goto 13
    2 write(6,'(''read in new value of it.'')')
    read(5,*)it
    goto 13
    3 write(6,'(''read in new value of nt.'')')
    read(5,*)nt
    goto 13
    4 write(6,'(''read in new value of dt.'')')
    read(5,*)dt
    goto 13
    5 write(6,'(''read in new value of x1.'')')
    read(5,*)x1
    goto 13
    6 write(6,'(''read in new value of x2.'')')
    read(5,*)x2
    goto 13
    7 write(6,'(''read in new value of dens.'')')
    read(5,*)dens
    goto 13
    8 write(6,'(''read in new value of temp.'')')
    read(5,*)temp
    goto 13
    11 return
    end subroutine input



    subroutine set

!  a subroutine to initialise the particle position and velocity
!  and to establish the background ion density. other calculation
!  parameters are also evaluated.
!  the particle parameters are determined using appropriate
!  random distributions.

    common /blkpos/ x(1000)
    common /blkvel/ v(1000)
    common /blkion/ rhop(1000)
    common /blkcel/ dx,xt,x1,x2
    common /blkdat/ dt,emt,ecf,dens,temp,vtemp
    common /blktim/ time
    common /blkint/ it,nt
    data idum/0/,iran/0/

!  wp is the plasma frequency, vtemp the thermal velocity,
!  factor the number of electrons per particle, emt electron
!  emt electron charge-to-mass ratio and ecf the charge-to-field
!  ratio e/epsilon0.

    wp=56.41457936*sqrt(dens)
    vtemp=sqrt(1.515623082e7*temp)
    debye=vtemp/wp
    xt=xt*debye
    x1=x1*debye
    x2=x1+x2*debye
    dx=xt/float(it)
    factor=dens*xt/float(nt)
!  this controls the timestep
    help=0.1*dx/vtemp
    dt=min(help,(dt/wp))
    emt=-1.758804786e11*dt
    ecf=1.809527009e-8*factor

!  set up the electron position and velocity distribution
!  the electron density distribution is uniform up to x1 and
!  then falls linearly to x2 and is zero to xt.

    alpha=(x1+x1)/(x1+x2)
    xx1=0.5*(x1+x2)
    xx2=x2*x2-x1*x1
    do 1 n=1,nt
        y=ran1(iran)
        if (y > alpha) then
            x(n)=x2-sqrt(xx2*(1.0-y))
        else
            x(n)=xx1*y
        endif
        v(n)=vtemp*gasdev(idum)
    1 end do
!  set up the background ion distribution
    rhoi=float(nt)*dx/xx1
    rhot=0.0
    xp=-0.5*dx
    xf=0.0
    do 2 i=1,it
        xb=xf
        xp=xp+dx
        xf=xf+dx
        if (xf < x1) then
        !  the point lies entirely within the uniform density zone
            rhop(i)=rhoi
        else if (xf < x2) then
        !  the front edge of the cell lies in the falling density
            if(xb < x1)then
            !  but the back in the uniform zone
                rhop(i)=rhoi*((x1-xb)+(xf-x1)*(x2-0.5*(xf+x1))/(x2-x1))/dx
            else
            !  the cell is entirely within the falling density region
                rhop(i)=rhoi*(x2-xp)/(x2-x1)
            endif
        else
        !  the front edge is in the zero density region
            if (xb < x1) then
            !  whilst the back is still in the uniform zone
                rhop(i)=rhoi*((x1-xb)+0.5*(x2-x1))/dx
            else if (xb < x2) then
            !  or the back is in the falling zone
                rhop(i)=rhoi*0.5*(x2-xb)*(x2-xb)/(dx*(x2-x1))
            else
            !  the cell is entirely in the zero density region
                rhop(i)=0.0
            endif
        endif
        rhot=rhot+rhop(i)
    2 end do
    if (nint(rhot) /= nt) then
        i=int(x1/dx)+1
        rhop(i)=rhop(i)+float(nt)-rhot
    endif
    return
    end subroutine set



    subroutine output(nout)
    character(2) :: a(9)
    character(8) :: filenam1,filenam2
    dimension xc(1000),ee(1000)
    common /blkpos/ x(1000)
    common /blkvel/ v(1000)
    common /blkfld/ e(1000)
    common /blkcel/ dx,xt,x1,x2
    common /blkdat/ dt,emt,ecf,dens,temp,vtemp
    common /blktim/ time
    common /blkint/ it,nt
    data a/'10','11','12','13','14','15','16','17','18'/
    vmax=0.0
    vmin=0.0
    do 1 n=1,nt
        vmax=amax1(v(n),vmax)
        vmin=amin1(v(n),vmin)
    1 end do
    vt=vmax-vmin
    emax=0.0
    emin=0.0
    do 2 i=1,it
        emax=amax1(e(i),emax)
        emin=amin1(e(i),emin)
    2 end do
    et=emax-emin
    do 3 i=1,nt
        x(i)=x(i)/xt
        v(i)=v(i)/vt
    3 end do
    xc(1)=0.5*dx/xt
    ee(1)=e(1)/et
    do 4 i=2,it
        xc(i)=xc(i-1)+dx/xt
        ee(i)=e(i)/et
    4 end do
    filenam1='v'//a(nout/60)//'.dat'
    filenam2='e'//a(nout/60)//'.dat'
    m1=10+nout/60
    m2=11+nout/60
    open(unit=m1,file=filenam1)
    open(unit=m2,file=filenam2)
    write(m1,50)(x(i),v(i),i=1,nt)
    write(m2,50)(xc(i),ee(i),i=1,it)
    50 format(2e14.5)
    do 5 i=1,nt
        x(i)=x(i)*xt
        v(i)=v(i)*vt
    5 end do
    return
    end subroutine output



    subroutine push

!  accelerate the particles in the electric field, and
!  move the electron position in response to its velocity

    common /blkpos/ x(1000)
    common /blkvel/ v(1000)
    common /blkfld/ e(1000)
    common /blkden/ rho(1000)
    common /blkdat/ dt,emt,ecf,dens,temp,vtemp
    common /blkcel/ dx,xt,x1,x2
    common /blkint/ it,nt
!  if the electron moves outside the mesh it is returned to
!  maintain overall charge neutrality.
!  if it leaves the dense boundary it is returned with a
!  random velocity from a thermal distribution.
!  if it leaves on the vacuum side it is returned with velocity
!  reversed.
    do 1 n=1,nt
        i=min0((1+int(x(n)/dx)),it)
        v(n)=v(n)+emt*e(i)
        x(n)=x(n)+dt*v(n)
        if (x(n) > xt ) then
            v(n)=-v(n)
            x(n)=xt+xt-x(n)
        else if (x(n) < 0.0) then
            v(n)=vtemp*abs(gasdev(idum))
            x(n)=0.0
        endif
    1 end do
    return
    end subroutine push



    subroutine field

!  the calculation of the field on the mesh

    common /blkpos/ x(1000)
    common /blkfld/ e(1000)
    common /blkion/ rhop(1000)
    common /blkden/ rho(1000)
    common /blkcel/ dx,xt,x1,x2
    common /blkdat/ dt,emt,ecf,dens,temp,vtemp
    common /blktim/ time
    common /blkint/ it,nt
!  set up the background ion charge
    do 10 i=1,it
        rho(i)=rhop(i)
    10 end do
!  identify the cell containing the particle and assign its
!  charge to that cell
    do 1 n=1,nt
        i=min0((1+int(x(n)/dx)),it)
        rho(i)=rho(i)-1.0
    1 end do
!  calculate the field at the cell centre with boundary condition
!  zero field at the dense boundary
    e1=0.0
    do 2 i=1,it
        e0=e1
        e1=e1+ecf*rho(i)
        e(i)=0.5*(e0+e1)
    2 end do
    return
    end subroutine field

    function ran1(idum)
!     returns a uniform deviate between 0.0 and 1.0. set idum
!     to any negative value to initialise or reinintialise
!     the sequence

    dimension r(97)
    parameter (m1=259200,ia1=7141,ic1=54773,rm1=1.0/m1)
    parameter (m2=134456,ia2=8121,ic2=28411,rm2=1.0/m2)
    parameter (m3=243000,ia3=4561,ic3=51349)
    save r,iff,ix1,ix2,ix3
    data iff /0/
!     initialise on first call even if idum is not zero
    if (idum < 0 .or. iff == 0) then
        iff=1
    !     seed the first routine
        ix1=mod((ic1-idum),m1)
        ix1=mod((ia1*ix1+ic1),m1)
    !     and use it to seed the second
        ix2=mod(ix1,m2)
        ix1=mod((ia1*ix1+ic1),m1)
    !     and the third routines
        ix3=mod(ix1,m3)
    !     fill the table with sequential uniform deviates generated
    !     by the first two routines
        do 11 j=1,97
            ix1=mod((ia1*ix1+ic1),m1)
            ix2=mod((ia2*ix2+ic2),m2)
        !     low and high order pieces combined here
            r(j)=(float(ix1)+float(ix2)*rm2)*rm1
        11 end do
        idum=1
    endif
!     except when initialising this is where we start.
!     generate the next number for each sequence
    ix1=mod((ia1*ix1+ic1),m1)
    ix2=mod((ia2*ix2+ic2),m2)
    ix3=mod((ia3*ix3+ic3),m3)
!     use the third sequence to get an integer between 1 and 97
    j=1+(97*ix3)/m3
    if (j > 97 .or. j < 1) write (*,*) ' failure in j'
!     return that table entry
    ran1=r(j)
!     and refill it
    r(j)=(float(ix1)+float(ix2)*rm2)*rm1
    return
    end function ran1

    function gasdev(idum)
!     returns a normally distributed deviate with zero mean
!     and unit variance, using ran1(idum) as the source of
!     uniform deviates

    save iset,gset
    data iset/0/
    1 if (iset == 0) then
    !     if no extra deviate is available
    !     pick two uniform numbers in the square extending from
    !     -1 to +1 in each direction
        v1=2.0*ran1(idum)-1.0
        v2=2.0*ran1(idum)-1.0
    !     test to check they lie in unit circle
        r=v1*v1+v2*v2
        if (r >= 1.0) go to 1
    !     make the box-muller transformation
        fac=sqrt(-2.0*alog(r)/r)
    !     to get two normal deviates
        gasdev=v1*fac
    !     we have an extra deviate available for the next call
        gset=v2*fac
        iset=1
    else
    !     use the spare deviate left from the previous call
        gasdev=gset
    !     unset the flag
        iset=0
    endif
    return
    end function gasdev


    block data
    common /blkpos/ x(1000)
    common /blkvel/ v(1000)
    common /blkfld/ e(1000)
    common /blkden/ rho(1000)
    common /blktim/ time
    data x/1000*0.0/
    data v/1000*0.0/
    data e/1000*0.0/
    data rho/1000*0.0/
    data time/0.0/
    end block data

