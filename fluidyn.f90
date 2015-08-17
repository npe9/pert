    program fluidyn
!----------------------------------------------------------------------

!  use double-length real numbers and long integers in this program

!----------------------------------------------------------------------
!  this is a cell-structure fluid dynamics program in which the cell
!  contains 125 molecules.   the molecules are first placed on a
!  5 x 5 x 5 regular grid and then displaced randomly by up to d/20
!  in each principal direction where d is the initial separation.
!  after 50 timesteps, to allow the system to settle down, the
!  molecular configuration is sampled every timestep to find the
!  quantity sum(f(i,j)r(i,j))/3 and also the radial density function.
!  the program transforms every quantity into s.i. units.  some form of
!  reduced units would keep numbers in a more convenient range but make
!  the program more difficult to follow.
    dimension x(125,3),v(125,3),delx(125,3,0:4),delv(125,3,0:4), &
    vx(3),separ(125,4),w(4),tab(0:100),dd(0:100),xt(125,3)
!  data for the random-number generator.  ir is the seed
    data ir,ix,iy,im/199,171,11213,53125/
!  weights for runge-kutta
    data w/0,0.5,0.5,1.0/
!  data for the lennard-jones force.   sig in angstrom units and
!  eps as eps/k where k is boltzmann's constant.  the figures
!  are for argon.
    data sig,eps/3.45,120/
!  convert to s.i.
    sig=sig*1.0e-10
    eps=eps*1.38e-23
!  am is the mass of the molecule in atomic mass units and t the
!  temperature of the liquid.  the mass given is for argon.
    data am,t/40,329/
!  the range of the force, dlim, is set at 2.4 x sig corresponding to
!  where the force falls to 1% of its maximum value.   the side of the
!  cubical cell is 5 x sig for v* = 1.   the value of v* is read in
!  which gives the cell edge, el.
    sig2=sig*sig
    tfeps=24.0*eps
!  mass of hydrogen atom.
    hm=1.667e-27
!  boltzmann constant
    cay=1.38e-23
!  number of timesteps to completion
    ntot=150
    pi=4.0*atan(1.0)
    dlim2=(2.4*sig)**2
    write(6,'('' read in the value of v* '')')
    read(5,*)vstar
    el=5.0*sig*vstar**(1.0/3.0)
    grid=0.2*el
!  the molecules are set on a uniform grid and then dispaced in the x,
!  y and z directions by up to el/20.0 using a random-number generator.
!  calculate the mean speed of the molecules = sqrt(3kt/am)
    vmean=sqrt(24835*t/am)
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
            !  now displace from grid positions
                do 2 l=1,3
                    ir=mod(ir*ix+iy,im)
                    x(num,l)=x(num,l)+(2*float(ir)/float(im)-1.0)*grid/20.0
                2 end do
            !  the speeds of the molecules should have a maxwell-boltzmann
            !  distribution.  however, the molecules are all given the same
            !  average speed, equal to sqrt(3kt/am) but in random directions.
            !  after some time the interactions should produce the correct
            !  distribution.   the model is run for 50 timesteps before any
            !  information is taken from it to allow this to happen.
                19 do 3 l=1,3
                    ir=mod(ir*ix+iy,im)
                    vx(l)=2*float(ir)/float(im)-1
                3 end do
                vv=sqrt(vx(1)**2+vx(2)**2+vx(3)**2)
                if(vv > 1.0)goto 19
                do 4 l=1,3
                    v(num,l)=vmean*vx(l)/vv
                4 end do
    1 end do
!  a variable timestep is used in this program. the initial h is chosen
!  so that the maximum distance travelled is 0.05*grid.
    h=grid/vmean/20.0
    nostep=0
    vir=0
    100 do 5 i=1,125
        do 5 j=1,3
            do 5 k=0,4
                delx(i,j,k)=0
                delv(i,j,k)=0
    5 end do
    nostep=nostep+1
    if((nostep/10)*10 /= nostep)goto 37
    write(6,250)nostep
    250 format(34h number of timesteps completed is ,i5)
    37 do 60 iw=1,4
        wt=w(iw)
    !  set all values of dx/dt, dy/dt and dz/dt and temporary values of
    !  x, y and z.
        do 16 i=1,125
            do 16 k=1,3
                delx(i,k,iw)=(v(i,k))+wt*h*delv(i,k,iw-1)
                xt(i,k)=x(i,k)+wt*h*delx(i,k,iw)
        16 end do
        do 18 i=1,124
        !  all the nearest neighbours are now identified for molecule i.  this
        !  involves finding the molecule j in the main cell or a ghost cell
        !  which is nearest and then checking to see if it is closer than dlim.
            do 18 j=i+1,125
                call dis(xt,separ,el,i,j)
                if(separ(j,4) > dlim2)goto 18
            !  the closest i-j separation has been found to be less than dlim.
            !  calculate constants for the lennard-jones force components.
                a=(sig2/separ(j,4))**3
                b=a/separ(j,4)
                c=tfeps*b*(1.0-2.0*a)/hm/am
            !  find acceleration components.
                do 31 k=1,3
                    delv(i,k,iw)=delv(i,k,iw)-c*separ(j,k)
                    delv(j,k,iw)=delv(j,k,iw)+c*separ(j,k)
                31 end do
        18 end do
    60 end do
!  update position and velocity components
    do 80 i=1,125
        do 80 k=1,3
            x(i,k)=x(i,k)+h*(delx(i,k,1)+delx(i,k,4)+2*(delx(i,k,2)+ &
            delx(i,k,3)))/6.0
            v(i,k)=v(i,k)+h*(delv(i,k,1)+delv(i,k,4)+2*(delv(i,k,2)+ &
            delv(i,k,3)))/6.0
    80 end do
!  restore molecule to be inside the box
    do 81 i=1,125
        do 81,k=1,3
            help=x(i,k)+9.5*el
            x(i,k)=mod(help,el)-0.5*el
    81 end do
    if(nostep <= 50)goto 100
!  find new timestep
    vmax=0
    do 71 i=1,125
        do 71 k=1,3
            if(abs(v(i,k)) > vmax)vmax=abs(v(i,k))
    71 end do
    h=grid/vmax/20.0
!  calculate contributions to the virial term and to the radial
!  distribution function.
    do 40 i=1,124
        do 41 j=i+1,125
            call dis(x,separ,el,i,j)
            if(separ(j,4) > dlim2)goto 41
            a=(sig2/separ(j,4))**3
            c=tfeps*a*(1.0-2.0*a)
        !  add contribution to virial term
            vir=vir-c/3.0
        !  add radial distance to table
            l=int(20.0*sqrt(abs(separ(j,4)))/sig)
            tab(l)=tab(l)+1
        41 end do
    40 end do
    if(nostep < ntot)goto 100
!  take the average of the virial term over 100 timesteps
    vir=vir/100.0
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
        tab(i)=0.96774*tab(i)*vstar/pi/((i+1.0)**3-i**3)
        dd(i)=0.05*(i+0.5)
    90 end do
    write(9,'('' '')')
    write(9,'('' radial density function on an absolute scale '')')
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



