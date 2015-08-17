      PROGRAM REACTOR
c  Programme to calculate the multiplication factor of an
c  finite homogeneous nuclear reactor of radius rr (m.)
c  consisting of a mixture of graphite and uranium in the atomic
c  ratio of A graphite: uranium. The uranium is enriched by a
c  factor E with respect to natural uranium (U235:U238::1:138).
c  The programme uses Monte-Carlo to set the initial neutron
c  energy via the random number generator DISDEV of appropriate
c  distribution, and the angle of deflection in collision.
c  The cross sections for scattering, absorption and fission are
c  from experimental data. The neutrons are considered as either
c  fast or thermal. Fast neutrons are absorbed by averaged
c  resonance absorption in the uranium, and lose energy by
c  elastic collisions in the moderator and fuel. Thermal neutrons
c  give fission, but lose no energy in elastic collisions.
      real :: one=1.d0
      logical fast
      character ans*1
      data idum,iran/-50,0/
      write (*,*) 'Uranium/graphite atomic fraction'
      read (*,*) a
      write (*,*) 'Enrichment factor'
      read (*,*) e
      write (*,*) 'Radius of fuel element'
      read (*,*) rr
      write (*,*) 'Number of trial neutrons'
      read (*,*) nt
   51 write (*,*) 'Do you want printed output (y/n)?'
      read(5,50)ans
   50 format(a1)
      iout=6
      if(ans.eq.'y'.or.ans.eq.'Y')goto 52
      if(ans.eq.'n'.or.ans.eq.'N')goto 53
      goto 51
   52 iout=9
      open(unit=9,file='lpt1')
   53 rr2=rr*rr
      fract35=7.1942e-03*e
      fract38=1.0-fract35
c   Calculate the number densities of carbon and uranium atoms
      concc=8.2775
      conc38=concc*a*fract38
      conc35=concc*a*fract35
      concc=concc*(1.0-a)
c   Calculate the total elastic scattering probabilty per unit 
c   path length
      sigmcs=4.8*concc
      sigmas=sigmcs+10.0*conc35+8.3*conc38
c   Calculate the resonance absorption collision probability
c   per unit path length
      sigman=sigmas/conc38
      if (sigman.lt.50.0) then
        sigmfr=0.02398*((50.0-sigman)*0.601+(sigman-8.3)*1.066)
      else if (sigman.lt.100.0) then
        sigmfr=0.02*((100.0-sigman)*1.066+(sigman-50.0)*1.506)
      else if (sigman.lt.300.0 )then
        sigmfr=0.005*((300.0-sigman)*1.506+(sigman-100.0)*2.493)
      else if (sigman.lt.500.0) then
        sigmfr=0.005*((500.0-sigman)*2.493+(sigman-300.0)*3.096)
      else if (sigman.lt.800.0) then
        sigmfr=0.003333333*((800.0-sigman)*3.096+
     1   (sigman-500.0)*3.754)
      else if (sigman.lt.1000.0) then
        sigmfr=0.005*((1000.0-sigman)*3.754+(sigman-800.0)*4.138)
      else if (sigman.lt.2000.0) then
        sigmfr=0.001*((2000.0-sigman)*4.138+(sigman-1000.0)*5.400)
      else
        sigmfr=15.936-21110.0/sigman
      endif
c   Calculate the total interaction probabilities per unit length
c   for fast and slow neutrons, and for fission
      sigmft=sigmas+3.2e-3*concc+sigmfr*conc38
      sigmst=sigmas+3.2e-3*concc+694.0*conc35+2.73*conc38
      sigmsf=582.0*conc35
c  Calculate the branching ratios for:
c  carbon/uranium elastic scattering BRC
c  fast scattering BRF
c  slow scattering BRS
c  slow absorption followed by fission BRN
      brc=sigmcs/sigmas
      brf=sigmas/sigmft
      brs=sigmas/sigmst
      brn=sigmsf/sigmst
c   Calculate the effective radius for diffusion
c            re = rr + 0.71 * mfp
      re=rr+0.71/sigmas
c   Initialise the total moderated weight to zero
      wm=0.0
c   Initialise the total thermal collision sum to zero
      wn=0.0
      do 1 n=1,nt
      if((n/100)*100.eq.n)write(*,*)n
c   Set the initial weight to 1 unit
      wt=1.0
      fast=.true.
c   Calculate the initial energy, direction and radius of the
c   neutron
      ener=1.0e6*disdev(idum)
      cosp=ran1(iran)
      cosp=-1.0+cosp+cosp
c   Set the current position of the neutron to the last fission
  100 r1=0.318309886*re*dendev(iran)
      if (r1.gt.rr) go to 100
c     r1=rr*ran1(iran)
c     r1=rr*(ran1(iran)**0.3333333333)
      r12=r1*r1
   10 r0=r1
      r02=r12
c   Calculate the minimum path length for escape
      dd=r0*cosp
      cc=rr2-r02
      dd=-(dd+sign(sqrt(dd*dd+cc),dd))
      if (dd.lt.0.0) dd=-cc/dd
c   Calculate the total cross section taking into account the
c   neutron energy
      if (ener.gt.0.025) then
        sigmat=sigmft
      else
        sigmat=sigmst
      endif
      tt=sigmat*dd
c   Calculate the probability of direct escape, and reduce the
c   neutron weight
      probt=(1.0-exp(-tt))
      wt=wt*probt
c   Calculate the free path within the reactor
c   Guard against zero exponent of alog
      z=1.0-probt*ran1(iran)
      if(z.lt.1.0e-20)z=1.0e-20      
      t=-alog(z)
      s=t/sigmat
      s2=s*s
c   Calculate the neutron collision point
      r12=r02+s2+(r0+r0)*s*cosp
      if (r12.gt.rr2) r12=rr2
      r1=sqrt(r12)
c   Calculate the direction angle at the scattering point
      cosp=(r12+s2-r02)/((r1+r1)*s)
c   Calculate the angle of scattering in the centre of mass
c   frame
      cost=ran1(iran)
      cost=-1.0+cost+cost
      phi=3.141592654*ran1(iran)
      cosphi=cos(phi)
      if (ener.gt.0.025) then
c   The neutron is fast and will be moderated
c   The neutron weight is reduced by the branching ratio
        wt=wt*brf
        if (ran1(iran).lt.brc) then
c   The neutron has collided with a carbon atom
          ae=0.858088+0.141912*cost
          cosp=(0.07686395*cosp+0.923136049*(cosp*cost+cosphi*
     1     sqrt(dim(one,(cost*cost))*dim(one,(cosp*cosp)))))/
     2     sqrt(ae)
          ener=ae*ener
        else
c   The neutron has collided with a uranium atom
          ae=0.991667854+0.00833214597*cost
          cosp=(0.00420115*cosp+0.995798848*(cosp*cost+cosphi*
     1     sqrt(dim(one,(cost*cost))*dim(one,(cosp*cosp)))))/
     2     sqrt(ae)
          ener=ae*ener
        endif
      else
c   The neutron is thermal and experiences no energy loss
        cosp=cosp*cost+cosphi*
     1   sqrt(dim(one,(cosp*cosp))*dim(one,(cost*cost)))
        if (fast) then
c   Identify thermalisation and add weight
          fast=.false.
          wm=wm+wt
        endif
c   Update the total thermal collision sum
        wn=wn+wt
c   and decrease the weight to account for absorption
        wt=wt*brs
      endif
c   If neutron weight is small either terminate path or increase weight
      if (wt.gt.0.1) go to 10
      if (ran1(iran).lt.(wt+wt)) then
        wt=0.5
        go to 10
      endif
    1 continue
c   Divide the total moderated weight by the number of neutrons
c   to get the thermalisation probability
      wm=wm/float(nt)
c   Divide the total thermal collision sum by the number of
c   neutrons and multiply by the fission branching ratio to
c   obtain the total probability of fission
      wn=brn*wn/float(nt)
c   Multiply by the average yield (2.47) to obtain the 
c   multilplication factor
      wf=2.47*wn
      write (iout,1000) a,e,rr,wm,wn,wf
 1000 format (' Atomic fuel ratio',1pe11.3/' Enrichment factor',
     11pe11.3/' Reactor radius',1pe11.3/
     2' Thermalisation probability',1pe11.3/
     3' Fission probability',1pe11.3,4x,
     4' Multiplication constant',1pe11.3)
      stop
      end


      function disdev(idum)
c
c   Calculates random deviates to the neutron fission product
c   energy probability distribution:
c           N(E)=C*SINH(SQRT(2E))*EXP(-E)
c
c   The deviate is obtained by the rejection method using
c   an exponential deviate as comparison
c
c  Calculate the exponential deviate
    1 y=-alog(ran1(idum))
c  The scaling factor 2.0464 is chosen to be the most
c  efficient value compatable with the ratio function f always
c  being less than 1.
      x=2.0464*y
c  Form comparison function
      f=0.76648*sinh(sqrt(x+x))*exp(-0.51134*x)
c  Reject if comparison function less than uniform deviate
      if (f.lt.ran1(idum)) go to 1
      disdev=x
      return
      end


      function dendev(idum)
c
c   Calculates random deviates to the neutron fission product
c   initial density probability distribution:
c           N(R)=C*R*SIN(R)
c
c   The deviate is obtained by the rejection method using
c   a Lorentzian deviate as comparison over a restricted range
c
c  Calculate the Lorentzian deviate limited in range
    1 y=-0.334750207+0.578242494*ran1(idum)
      y=tan(3.141592654*y)
c  Shift mean to peak of density distribution
c  The scaling factor 1.159294016 is chosen to be the smallest
c  value compatable with the ratio function f always being
c  less than 1.
      x=1.159294016*y+2.028757838
c  Reject if value negative or value greater than PI
      if ((x.lt.0.0).or.(x.gt.3.141592654)) go to 1
c  Form comparison function
      f=0.425183296*x*sin(x)*(1.0+y*y)
c  Reject if comparison function less than uniform deviate
      if (f.lt.ran1(idum)) go to 1
      dendev=x
      return
      end


      function ran1(idum)
c     Returns a uniform deviate between 0.0 and 1.0. Set IDUM
c     to any negative value to initialise or reinintialise
c     the sequence
c
      dimension r(97)
      parameter (m1=259200,ia1=7141,ic1=54773,rm1=1.0/m1)
      parameter (m2=134456,ia2=8121,ic2=28411,rm2=1.0/m2)
      parameter (m3=243000,ia3=4561,ic3=51349)
      save r,iff,ix1,ix2,ix3
      data iff /0/
c     Initialise on first call even if IDUM is not zero
      if (idum.lt.0.or.iff.eq.0) then
      iff=1
c     Seed the first routine
      ix1=mod((ic1-idum),m1)
      ix1=mod((ia1*ix1+ic1),m1)
c     and use it to seed the second
      ix2=mod(ix1,m2)
      ix1=mod((ia1*ix1+ic1),m1)
c     and the third routines
      ix3=mod(ix1,m3)
c     Fill the table with sequential uniform deviates generated
c     by the first two routines
      do 11 j=1,97
         ix1=mod((ia1*ix1+ic1),m1)
         ix2=mod((ia2*ix2+ic2),m2)
c     Low and high order pieces combined here
         r(j)=(float(ix1)+float(ix2)*rm2)*rm1
   11 continue
      idum=1
      endif
c     Except when initialising this is where we start.
c     Generate the next number for each sequence
      ix1=mod((ia1*ix1+ic1),m1)
      ix2=mod((ia2*ix2+ic2),m2)
      ix3=mod((ia3*ix3+ic3),m3)
c     Use the third sequence to get an integer between 1 and 97
      j=1+(97*ix3)/m3
      if (j.gt.97.or.j.lt.1) write (*,*) ' failure in j'
c     Return that table entry
      ran1=r(j)
c     and refill it
      r(j)=(float(ix1)+float(ix2)*rm2)*rm1
      return
      end
