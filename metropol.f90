    program metropol
!----------------------------------------------------------------------

!  use double-length real numbers and long integers in this program

!----------------------------------------------------------------------
!  this is a cell-structure monte carlo program in which the cell
!  contains 125 molecules.   the molecules are first placed on a
!  5 x 5 x 5 regular grid and then displaced randomly by up to d/20
!  in each principal direction where d is the initial separation.
!  after 20000 configurations, to allow the system to settle down,
!  every 50th molecular configuration is sampled to find the quantity
!  sum(f(i,j)r(i,j))/3 and also the radial density function for the
!  next 80000 configurations.
!  the program transforms every quantity into s.i. units.  some form of
!  reduced units would keep numbers in a more convenient range but make
!  the program more difficult to follow.

!  the algorithm is based on that in n. metropolis, a.w.rosenbluth,
!  m.n. rosenbluth and a.h. teller, j. chem. phys. 21, 1087.

    dimension x(125,3),pot(125,125),tempot(125),rr(4), &
    separ(125,4),tab(0:100),dd(0:100)
!  data for the random-number generator.  ir is the seed
    data ir,ix,iy,im/199,171,11213,53125/
!  data for the lennard-jones force.   sig is in nm and
!  eps as eps/k where k is boltzmann's constant.  the figures
!  are for argon.
    data sig,eps/0.345,120/
!  convert to s.i.
    sig=sig*1.0e-9
    eps=eps*1.38e-23
!  t is the temperature of the liquid.
    data t/329/
!  the range of the force, dlim, is set at 2.4 x sig corresponding to
!  where the force falls to 1% of its maximum value.   the side of the
!  cubical cell is 5 x sig for v* = 1.   the value of v* is read in
!  which gives the cell edge, el.
    sig2=sig*sig
    tfeps=24.0*eps
!  boltzmann constant
    cay=1.38e-23
!  number of generated configurations
    ntot=100000
    pi=4.0*atan(1.0)
    dlim2=(2.4*sig)**2
    write(6,'('' read in the value of v* '')')
    read(5,*)vstar
    el=5.0*sig*vstar**(1.0/3.0)
    grid=0.2*el
    vir=0
!  the molecules are initially set on a uniform grid.
!  clear table for radial distribution
    do 17 i=1,100
        tab(i)=0
    17 end do
    num=0
    do 1 i=1,5
        do 1 j=1,5
            do 1 k=1,5
                num=num+1
                x(num,1)=(i-3)*grid
                x(num,2)=(j-3)*grid
                x(num,3)=(k-3)*grid
    1 end do
!  clear potential energy table
    do 15 ii=1,125
        do 15 jj=1,125
            pot(ii,jj)=0
    15 end do
!  calculate initial potential energy table
    do 16 ii=1,124
        do 16 jj=ii+1,125
            call dis(x,separ,el,ii,jj)
            if(separ(jj,4) > dlim2)goto 16
            a=(sig2/separ(jj,4))**3
            pot(ii,jj)=4*eps*a*(a-1.0)
            pot(jj,ii)=pot(ii,jj)
    16 end do
    nostep=0
    100 nostep=nostep+1
    if((nostep/50)*50 == nostep)write(6,125)nostep
    125 format(15h starting step ,i5)
!  choose a molecule
    ir=mod(ir*ix+iy,im)
    mol=int(num*float(ir)/float(im))+1
!  choose a direction
    98 do 19 k=1,3
        ir=mod(ir*ix+iy,im)
        rr(k)=2*float(ir)/float(im)-1.0
    19 end do
    rr(4)=rr(1)**2+rr(2)**2+rr(3)**2
    if(rr(4) > 1.0)goto 98
    rs=sqrt(rr(4))
!  choose a distance for the molecule to be moved
    ir=mod(ir*ix+iy,im)
    del=0.25*grid*float(ir)/float(im)
!  clear temporary potential table
    do 27 j=1,125
        tempot(j)=0
    27 end do
!  calculate new potential contributions
    do 21 j=1,125
        if(j == mol)goto 21
        do 22 k=1,3
            separ(j,k)=x(mol,k)+del*rr(k)/rs-x(j,k)
            if(abs(separ(j,k)+el) < abs(separ(j,k)))separ(j,k)=separ(j,k)+el
            if(abs(separ(j,k)-el) < abs(separ(j,k)))separ(j,k)=separ(j,k)-el
        22 end do
        separ(j,4)=separ(j,1)**2+separ(j,2)**2+separ(j,3)**2
        if(separ(j,4) > dlim2)goto 21
        a=(sig2/separ(j,4))**3
        tempot(j)=4*eps*a*(a-1.0)
    21 end do
!  find deltaphi
    sum1=0
    sum2=0
    do 28 j=1,125
        sum1=sum1+pot(mol,j)
        sum2=sum2+tempot(j)
    28 end do
    delphi=sum2-sum1
!  test delphi
    if(delphi < 0)then
        call change(pot,tempot,x,rr,del,el,mol)
        goto 55
    endif
!  delphi is positive. decide on next configuration
    z=exp(-delphi/cay/t)
    ir=mod(ir*ix+iy,im)
    r=float(ir)/float(im)
    if(z > r)then
        call change(pot,tempot,x,rr,del,el,mol)
    endif
    55 if(nostep <= 50000)goto 100
    if((nostep/10)*10 /= nostep)goto 100
!  calculate contributions to the virial term and to the radial
!  distribution function.
    do 40 i=1,124
        do 41 j=i+1,125
            call dis(x,separ,el,i,j)
            if(separ(j,4) > dlim2)goto 47
            a=(sig2/separ(j,4))**3
            c=tfeps*a*(1.0-2.0*a)
        !  add contribution to virial term
            vir=vir-c/3.0
        !  add radial distance to table
            47 l=int(20.0*sqrt(abs(separ(j,4)))/sig)
            if(l > 48)goto 41
            tab(l)=tab(l)+1
        41 end do
    40 end do
    if(nostep < ntot)goto 100
!  take the average of the virial term over 1600 timesteps
    vir=vir/5000.0
!  now calculate the factor 1 + vir/nkt
    factor=1.0+vir/cay/t/float(num)
    open(unit=9,file='lpt1')
!  output temperature and vstar
    write(9,350)t,vstar
    350 format(15h temperature = ,f6.1,9h vstar = ,f7.2)
    write(9,'('' '')')
!  output factor
    write(9,200)factor
    200 format(21h the factor pv/nkt = ,f7.3)
!  output the radial density distribution 4*pi*r**2*ro(r) on an
!  absolute scale relative to the average density.
    do 90 i=0,47
        tab(i)=6.161e-3*vstar*tab(i)/((i+1.0)**3-i**3)
        dd(i)=0.05*(i+0.5)
    90 end do
    write(9,'('' '')')
    write(9,'('' radial density function on an abolute scale '')')
    write(9,300)(dd(i),tab(i),i=0,47)
    300 format(2(f10.3,f12.6))
    stop
    end program




    subroutine dis(x,s,el,i,j)
    dimension x(125,3),s(125,4)
    do 20 k=1,3
        s(j,k)=x(i,k)-x(j,k)
        if(abs(s(j,k)+el) < abs(s(j,k)))s(j,k)=s(j,k)+el
        if(abs(s(j,k)-el) < abs(s(j,k)))s(j,k)=s(j,k)-el
    20 end do
    s(j,4)=s(j,1)**2+s(j,2)**2+s(j,3)**2
    return
    end subroutine dis



    subroutine change(p,t,x,r,d,e,m)
    dimension p(125,125),t(125),x(125,3),r(4)
    rs=sqrt(r(4))
    do 2 j=1,125
        p(m,j)=t(j)
        p(j,m)=t(j)
    2 end do
!  move molecule m and place in cell
    do 3 k=1,3
        x(m,k)=x(m,k)+d*r(k)/rs
        x(m,k)=amod(x(m,k)+9.5*e,e)-0.5*e
    3 end do
    return
    end subroutine change
