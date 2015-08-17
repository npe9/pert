    program reactor
!  programme to calculate the multiplication factor of an
!  finite homogeneous nuclear reactor of radius rr (m.)
!  consisting of a mixture of graphite and uranium in the atomic
!  ratio of a graphite: uranium. the uranium is enriched by a
!  factor e with respect to natural uranium (u235:u238::1:138).
!  the programme uses monte-carlo to set the initial neutron
!  energy via the random number generator disdev of appropriate
!  distribution, and the angle of deflection in collision.
!  the cross sections for scattering, absorption and fission are
!  from experimental data. the neutrons are considered as either
!  fast or thermal. fast neutrons are absorbed by averaged
!  resonance absorption in the uranium, and lose energy by
!  elastic collisions in the moderator and fuel. thermal neutrons
!  give fission, but lose no energy in elastic collisions.
    real :: one=1.d0
    logical :: fast
    character ans*1
    data idum,iran/-50,0/
    write (*,*) 'uranium/graphite atomic fraction'
    read (*,*) a
    write (*,*) 'enrichment factor'
    read (*,*) e
    write (*,*) 'radius of fuel element'
    read (*,*) rr
    write (*,*) 'number of trial neutrons'
    read (*,*) nt
    51 write (*,*) 'do you want printed output (y/n)?'
    read(5,50)ans
    50 format(a1)
    iout=6
    if(ans == 'y' .or. ans == 'y')goto 52
    if(ans == 'n' .or. ans == 'n')goto 53
    goto 51
    52 iout=9
    open(unit=9,file='lpt1')
    53 rr2=rr*rr
    fract35=7.1942e-03*e
    fract38=1.0-fract35
!   calculate the number densities of carbon and uranium atoms
    concc=8.2775
    conc38=concc*a*fract38
    conc35=concc*a*fract35
    concc=concc*(1.0-a)
!   calculate the total elastic scattering probabilty per unit
!   path length
    sigmcs=4.8*concc
    sigmas=sigmcs+10.0*conc35+8.3*conc38
!   calculate the resonance absorption collision probability
!   per unit path length
    sigman=sigmas/conc38
    if (sigman < 50.0) then
        sigmfr=0.02398*((50.0-sigman)*0.601+(sigman-8.3)*1.066)
    else if (sigman < 100.0) then
        sigmfr=0.02*((100.0-sigman)*1.066+(sigman-50.0)*1.506)
    else if (sigman < 300.0 )then
        sigmfr=0.005*((300.0-sigman)*1.506+(sigman-100.0)*2.493)
    else if (sigman < 500.0) then
        sigmfr=0.005*((500.0-sigman)*2.493+(sigman-300.0)*3.096)
    else if (sigman < 800.0) then
        sigmfr=0.003333333*((800.0-sigman)*3.096+ &
        (sigman-500.0)*3.754)
    else if (sigman < 1000.0) then
        sigmfr=0.005*((1000.0-sigman)*3.754+(sigman-800.0)*4.138)
    else if (sigman < 2000.0) then
        sigmfr=0.001*((2000.0-sigman)*4.138+(sigman-1000.0)*5.400)
    else
        sigmfr=15.936-21110.0/sigman
    endif
!   calculate the total interaction probabilities per unit length
!   for fast and slow neutrons, and for fission
    sigmft=sigmas+3.2e-3*concc+sigmfr*conc38
    sigmst=sigmas+3.2e-3*concc+694.0*conc35+2.73*conc38
    sigmsf=582.0*conc35
!  calculate the branching ratios for:
!  carbon/uranium elastic scattering brc
!  fast scattering brf
!  slow scattering brs
!  slow absorption followed by fission brn
    brc=sigmcs/sigmas
    brf=sigmas/sigmft
    brs=sigmas/sigmst
    brn=sigmsf/sigmst
!   calculate the effective radius for diffusion
!            re = rr + 0.71 * mfp
    re=rr+0.71/sigmas
!   initialise the total moderated weight to zero
    wm=0.0
!   initialise the total thermal collision sum to zero
    wn=0.0
    do 1 n=1,nt
        if((n/100)*100 == n)write(*,*)n
    !   set the initial weight to 1 unit
        wt=1.0
        fast= .true. 
    !   calculate the initial energy, direction and radius of the
    !   neutron
        ener=1.0e6*disdev(idum)
        cosp=ran1(iran)
        cosp=-1.0+cosp+cosp
    !   set the current position of the neutron to the last fission
        100 r1=0.318309886*re*dendev(iran)
        if (r1 > rr) go to 100
    !     r1=rr*ran1(iran)
    !     r1=rr*(ran1(iran)**0.3333333333)
        r12=r1*r1
        10 r0=r1
        r02=r12
    !   calculate the minimum path length for escape
        dd=r0*cosp
        cc=rr2-r02
        dd=-(dd+sign(sqrt(dd*dd+cc),dd))
        if (dd < 0.0) dd=-cc/dd
    !   calculate the total cross section taking into account the
    !   neutron energy
        if (ener > 0.025) then
            sigmat=sigmft
        else
            sigmat=sigmst
        endif
        tt=sigmat*dd
    !   calculate the probability of direct escape, and reduce the
    !   neutron weight
        probt=(1.0-exp(-tt))
        wt=wt*probt
    !   calculate the free path within the reactor
    !   guard against zero exponent of alog
        z=1.0-probt*ran1(iran)
        if(z < 1.0e-20)z=1.0e-20
        t=-alog(z)
        s=t/sigmat
        s2=s*s
    !   calculate the neutron collision point
        r12=r02+s2+(r0+r0)*s*cosp
        if (r12 > rr2) r12=rr2
        r1=sqrt(r12)
    !   calculate the direction angle at the scattering point
        cosp=(r12+s2-r02)/((r1+r1)*s)
    !   calculate the angle of scattering in the centre of mass
    !   frame
        cost=ran1(iran)
        cost=-1.0+cost+cost
        phi=3.141592654*ran1(iran)
        cosphi=cos(phi)
        if (ener > 0.025) then
        !   the neutron is fast and will be moderated
        !   the neutron weight is reduced by the branching ratio
            wt=wt*brf
            if (ran1(iran) < brc) then
            !   the neutron has collided with a carbon atom
                ae=0.858088+0.141912*cost
                cosp=(0.07686395*cosp+0.923136049*(cosp*cost+cosphi* &
                sqrt(dim(one,(cost*cost))*dim(one,(cosp*cosp)))))/ &
                sqrt(ae)
                ener=ae*ener
            else
            !   the neutron has collided with a uranium atom
                ae=0.991667854+0.00833214597*cost
                cosp=(0.00420115*cosp+0.995798848*(cosp*cost+cosphi* &
                sqrt(dim(one,(cost*cost))*dim(one,(cosp*cosp)))))/ &
                sqrt(ae)
                ener=ae*ener
            endif
        else
        !   the neutron is thermal and experiences no energy loss
            cosp=cosp*cost+cosphi* &
            sqrt(dim(one,(cosp*cosp))*dim(one,(cost*cost)))
            if (fast) then
            !   identify thermalisation and add weight
                fast= .false. 
                wm=wm+wt
            endif
        !   update the total thermal collision sum
            wn=wn+wt
        !   and decrease the weight to account for absorption
            wt=wt*brs
        endif
    !   if neutron weight is small either terminate path or increase weight
        if (wt > 0.1) go to 10
        if (ran1(iran) < (wt+wt)) then
            wt=0.5
            go to 10
        endif
    1 end do
!   divide the total moderated weight by the number of neutrons
!   to get the thermalisation probability
    wm=wm/float(nt)
!   divide the total thermal collision sum by the number of
!   neutrons and multiply by the fission branching ratio to
!   obtain the total probability of fission
    wn=brn*wn/float(nt)
!   multiply by the average yield (2.47) to obtain the
!   multilplication factor
    wf=2.47*wn
    write (iout,1000) a,e,rr,wm,wn,wf
    1000 format (' atomic fuel ratio',1pe11.3/' enrichment factor', &
    & 1pe11.3/' reactor radius',1pe11.3/ &
    ' thermalisation probability',1pe11.3/ &
    ' fission probability',1pe11.3,4x, &
    ' multiplication constant',1pe11.3)
    stop
    end program


    function disdev(idum)

!   calculates random deviates to the neutron fission product
!   energy probability distribution:
!           n(e)=c*sinh(sqrt(2e))*exp(-e)

!   the deviate is obtained by the rejection method using
!   an exponential deviate as comparison

!  calculate the exponential deviate
    1 y=-alog(ran1(idum))
!  the scaling factor 2.0464 is chosen to be the most
!  efficient value compatable with the ratio function f always
!  being less than 1.
    x=2.0464*y
!  form comparison function
    f=0.76648*sinh(sqrt(x+x))*exp(-0.51134*x)
!  reject if comparison function less than uniform deviate
    if (f < ran1(idum)) go to 1
    disdev=x
    return
    end function disdev


    function dendev(idum)

!   calculates random deviates to the neutron fission product
!   initial density probability distribution:
!           n(r)=c*r*sin(r)

!   the deviate is obtained by the rejection method using
!   a lorentzian deviate as comparison over a restricted range

!  calculate the lorentzian deviate limited in range
    1 y=-0.334750207+0.578242494*ran1(idum)
    y=tan(3.141592654*y)
!  shift mean to peak of density distribution
!  the scaling factor 1.159294016 is chosen to be the smallest
!  value compatable with the ratio function f always being
!  less than 1.
    x=1.159294016*y+2.028757838
!  reject if value negative or value greater than pi
    if ((x < 0.0) .or. (x > 3.141592654)) go to 1
!  form comparison function
    f=0.425183296*x*sin(x)*(1.0+y*y)
!  reject if comparison function less than uniform deviate
    if (f < ran1(idum)) go to 1
    dendev=x
    return
    end function dendev


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
