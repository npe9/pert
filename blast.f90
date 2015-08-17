
!     a one-dimensional lagrangian code using the von neumann-
!     richtmyer algorithm to calculate the dynamics of a simple
!     blast wave propogating in a cold gas.
!     the calculation is performed to compare with the self-similar
!     model of sedov and taylor.
!     the calculation may be performed for a planar, cylindical or
!     spherical blast wave.

    common /bkit/ it,iout
    common /bkdata/ gama,gama2
    common /artvis/ qp(200)
    common /dens/ rhob(200),rhoa(200),dm(200),dmt(200)
    common /ener/ eb(200),ea(200)
    common /incre/ dt,dta,dtb,dtc,dtd
    common /init/ jm,jmp1,id
    common /radial/ rb(200),ra(200),r2b(200),r2a(200)
    common /tim/ time,ttime,totout,tout
    common /tot/ ekinb,ekina,eintb,einta,etotb,etota,error
    common /value/ rho,ein,rout,dim
    common /vel/ ua(200),fj(200)
    open (21,file='blast.dat')
    open (22,file='blast.out')
    open (30,file='blast0.gra')
    open (31,file='blast1.gra')
    open (32,file='blast2.gra')
    open (33,file='blast3.gra')
    open (34,file='blast4.gra')
    open (35,file='blast5.gra')
    open (36,file='blast6.gra')
    open (37,file='blast7.gra')
    open (38,file='blast8.gra')
    open (39,file='blast9.gra')
    iout=0
    call input
    call set
    1 call output
    call timing
    call accln
    call radius
    call densty
    call visty
    call energy
    call toten
    call update
    it=it+1
    if (time <= ttime) go to 1
    call output
    stop
    end program


    subroutine input

!  reads the input conditions

    common /incre/ dtime,dta,dtb,dtc,dtd
    common /init/ jm,jmp1,id
    common /tim/ time,ttime,totout,tout
    common /value/ rho,ein,rout,dim
    common /bkdata/ gama,gama2
    read (21,1000) jm,id
    1000 format (2i6)
    read (21,1001) rho,ein,rout,gama
    1001 format (4e11.3)
    read (21,1002) dtime,ttime,tout
    1002 format (3e11.3)
    write (22,1101) dtime,ttime,tout,ein,rout,rho,gama,jm,id
    1101 format(1h0,4x,'initial time step',e11.3,4x,'total run time',e11.3, &
    & 4x,'output interval',e11.3/5x,'total energy :',e11.3,4x, &
    'outer radius : ',e11.3/5x,'initial density :',e11.3/5x,'gamma :', &
    e11.3/5x,'no of cells :',i4/5x,'dimension coefficient',i6)
    return
    end subroutine input



    subroutine set

!  sets the constant terms in the flow and initialises the variables

    common /dens/ rhob(200),rhoa(200),dm(200),dmt(200)
    common /ener/ eb(200),ea(200)
    common /init/ jm,jmp1,id
    common /pres/ p(200)
    common /radial/ rb(200),ra(200),r2b(200),r2a(200)
    common /tim/ time,ttime,totout,tout
    common /tot/ ekinb,ekina,eintb,einta,etotb,etota,error
    common /value/ rho,ein,rout,dim
    common /bkdata/ gama,gama2
    common /bkgrap/ alpha,rhosh,rholt
    jmp1=jm+1
    dim=float(iabs(id))
    alpha=sqrt(ein/rho)
    if ((id == 3) .and. (abs(gama-1.4) < 1.0e-3)) &
    alpha=1.083974169*alpha
    rmove=1.2*(alpha*ttime)**(2.0/(3.0+float(id)))
    if (rout < rmove) rout=rmove
    if (id == 2) ein=0.159154943*ein
    if (id == 3) ein=0.079577471*ein
    rholt=rho+rho
    rhosh=(gama+1.0)/(gama-1.0)*rho
    gama=gama-1.0
    gama2=0.5*gama
    totout=0.0
    if ((jm == 1) .or. (id < 0)) then
        beta=2.0*rout/(float(jm)*float(jmp1))
        dr1=0.5*beta
        j=1
        ra(1)=0.0
        do 10 j1=2,jmp1
            dr0=dr1
            dr1=beta*(float(j)-0.5)
            dr=0.5*(dr0+dr1)
            ra(j1)=ra(j)+dr
            rhoa(j)=rho
            rhob(j)=rho
            j=j1
        10 end do
    else
        dr=rout/float(jm)
        ra(1)=0.0
        j=1
        do 11 j1=2,jmp1
            ra(j1)=ra(j)+dr
            rhoa(j)=rho
            rhob(j)=rho
            j=j1
        11 end do
    endif
    id=iabs(id)-1
    do 2 j=2,jmp1
        r2a(j)=ra(j)**id
    2 end do
    do 3 j=1,jmp1
        r2b(j)=r2a(j)
        rb(j)=ra(j)
    3 end do
    rho0=0.0
    do 4 j=1,jm
        rho1=ra(j+1)*r2a(j+1)
        dm(j)=rhoa(j)*(rho1-rho0)/dim
        rho0=rho1
    4 end do
    dm(jmp1)=0.0
    rhoa(jmp1)=0.0
    rhob(jmp1)=0.0
    ea(1)=ein/dm(1)
    eb(1)=ea(1)
    etota=ein
    return
    end subroutine set



    subroutine output

!  output the flow variables

    dimension r0(200),e0(200),v0(200),d0(200)
    common /artvis/ qp(200)
    common /dens/ rhob(200),rhoa(200),dm(200),dmt(200)
    common /ener/ eb(200),ea(200)
    common /incre/ dt,dta,dtb,dtc,dtd
    common /init/ jm,jmp1,id
    common /radial/ rb(200),ra(200),r2b(200),r2a(200)
    common /tim/ time,ttime,totout,tout
    common /tot/ ekinb,ekina,eintb,einta,etotb,etota,error
    common /value/ rho,ein,rout,dim
    common /vel/ ua(200),fj(200)
    common /bkdata/ gama,gama2
    common /bkit/ it,iout
    common /bkgrap/ alpha,rhosh,rholt
    save r0,e0,v0,d0
    save rb0,time0
    data r0/200*0.0/,e0/200*0.0/,v0/200*0.0/,d0/200*0.0/
    data rb0/0.0/,time0/0.0/
!  plot blast wave position
    if ((time > 1.0e-20) .and. (time < ttime)) then
        jj=jm
        do 10 j=2,jm
            jj=jj-1
            if (rhoa(jj) > rholt) then
                xj1=float(jj)-0.5+(rholt-rhoa(jj))/(rhoa(jj+1)-rhoa(jj))
                go to 11
            endif
        10 end do
        xj1=0.5
        11 rb1=xj1/float(jm)*rout
        time0=time
        rb0=rb1
    endif
    if (time < totout) return
!  print out full data report
    iout=iout+1
    write (22,20)it,time,rb1,etota,ekina,einta,error
    20 format(14h1iteration no.,i8,4x,4htime,1pe14.6,'blast wave radius', &
    &  1pe12.4// &
    & 8h energy:,5x,5htotal,1pe11.3,3x,7hkinetic, &
    & 1pe11.3,3x,7hthermal,1pe11.3,5x,5herror,1pe11.3// &
    & 1h ,3x,8hcell no.,5x,6hradius,7x,8hvelocity,5x,7hdensity,5x, &
    & 9hint ener.,4x,9hart.visc.)
    do 2 j=1,jmp1
        write(22,21) j,ra(j),ua(j),rhoa(j),ea(j),qp(j)
    2 end do
    21 format(4x,i4,3x,1pe14.6,1p5e13.3)
    write(22,22) dtd,dtc,dta
    22 format(10h0time step,5x,4hcomp,1pe11.3,5x,3hcfl,1pe11.3, &
    & 5x,5hfinal,1pe11.3)
    totout=totout+tout
    if ((time < 1.0e-20) .or. (time > ttime)) goto 40
    u0=rb1/time
    ee0=(gama+1.0)/gama*((u0/(3.0+float(id)))**2)
    do 4 j=1,jm
        r0(j)=0.5*(ra(j)+ra(j+1))/ra(jmp1)
        d0(j)=rhoa(j)/rhosh
        v0(j)=0.5*(ua(j)+ua(j+1))/u0
        e0(j)=ea(j)/ee0
    4 end do
!  produce files for up to 10 graphical outputs
    40 if(iout > 10)return
    m=29+iout
    do 30 i=1,jmp1
        write(m,*)ra(i),rhoa(i),ua(i),ea(i)
    30 end do
    return
    end subroutine output


    subroutine timing

!  time step adjustment

    common /dens/ rhob(200),rhoa(200),dm(200),dmt(200)
    common /ener/ eb(200),ea(200)
    common /incre/ dt,dta,dtb,dtc,dtd
    common /init/ jm,jmp1,id
    common /radial/ rb(200),ra(200),r2b(200),r2a(200)
    common /pres/ p(200)
    common /tot/ ekinb,ekina,eintb,einta,etotb,etota,error
    common /value/ rho,ein,rout,dim
    common /vel/ ua(200),fj(200)

!  the time step adjustment is performed in four parts

!     the time step may increase by a factor of 2
    dtd=dt+dt
    dtc=dtd*dtd
!  the courant-lewy-friedrichs condition
    do 1 j=1,jm
        speed2=10.0*p(j)/rhoa(j)
        help=1.0e-10
        dt1=(ra(j+1)-ra(j))**2 /max(help,speed2)
        dtc=amin1(dt1,dtc)
    1 end do
    dtc=sqrt(dtc)
!  the density increase is not allowed to exceed 0.1
    j=1
    do 2 j1=2,jm
        dt1=0.1*(ra(j1)-ra(j))/max(help,abs(ua(j1)-ua(j)))
        dtd=amin1(dt1,dtd)
        j=j1
    2 end do
    dtb=amin1(dtc,dtd)
    if (dtb > dta .and. dtb < 1.999*dta) dtb=dta
    dt=0.5*(dta+dtb)
    return
    end subroutine timing


    subroutine accln

!  velocity incremental calculation

    common /artvis/ qp(200)
    common /dens/ rhob(200),rhoa(200),dm(200),dmt(200)
    common /ener/ eb(200),ea(200)
    common /init/ jm,jmp1,id
    common /pres/ p(200)
    common /radial/ rb(200),ra(200),r2b(200),r2a(200)
    common /incre/ dt,dta,dtb,dtc,dtd
    common /vel/ ua(200),fj(200)
    common /bkit/ it,iout
    data t/0.0/

!     the velocity increments are calculated

    if (it /= 0) t=(dta-dtb)/(4.0*dta)
    force=0.0
    j1=jmp1
    j=jmp1
    do 1 jj=1,jm
        j=j-1
        force1=p(j)+qp(j)
        force=2.0*r2a(j1)*(force-force1)/(dm(j)+dm(j1))
        ua(j1)=ua(j1)-((1.0-t)*force+t*fj(j1))*dt
        fj(j1)=force
        force=force1
        j1=j
    1 end do
    ua(1)=0.0
    dt=dtb
    dta=dtb
    return
    end subroutine accln

    subroutine radius

!   calculates the radial increments

    common /init/ jm,jmp1,id
    common /radial/ rb(200),ra(200),r2b(200),r2a(200)
    common /incre/ dt,dta,dtb,dtc,dtd
    common /vel/ ua(200),fj(200)
    do 1 j=2,jmp1
        rb(j)=ra(j)+ua(j)*dt
        r2b(j)=rb(j)**id
    1 end do
    return
    end subroutine radius


    subroutine densty

!  calculates the density

    common /dens/ rhob(200),rhoa(200),dm(200),dmt(200)
    common /init/ jm,jmp1,id
    common /incre/ dt,dta,dtb,dtc,dtd
    common /radial/ rb(200),ra(200),r2b(200),r2a(200)
    common /value/ rho,ein,rout,dim
    rho0=0.0
    j=1
    do 1 j1=2,jmp1
        rho1=rb(j1)*r2b(j1)
        rhob(j)=dim*dm(j)/(rho1-rho0)
        rho0=rho1
        j=j1
    1 end do
    rhob(jmp1)=rhob(jm)
    return
    end subroutine densty


    subroutine visty

!  calculates the von neumann artificial viscosity

    common /artvis/ qp(200)
    common /dens/ rhob(200),rhoa(200),dm(200),dmt(200)
    common /init/ jm,jmp1,id
    common /vel/ ua(200),fj(200)

!     calculation of the artificial viscosity

    do 1 j=1,jm
        j1=j+1
        qp(j)=0.0
        if (rhob(j) > rhoa(j) .and. ua(j1) < ua(j)) &
        qp(j)=0.5*(rhoa(j)+rhob(j))*(ua(j1)-ua(j))**2
    1 end do
    return
    end subroutine visty


    subroutine energy

!  energy transfer without thermal conduction

    common /artvis/ qp(200)
    common /dens/ rhob(200),rhoa(200),dm(200),dmt(200)
    common /ener/ eb(200),ea(200)
    common /incre/ dt,dta,dtb,dtc,dtd
    common /init/ jm,jmp1,id
    common /intern/ a(200),b(200),c(200),d(200)
    common /radial/ rb(200),ra(200),r2b(200),r2a(200)
    common /tim/ time,ttime,totout,tout
    common /bkdata/ gama,gama2
    do 2 j=1,jm
        x=rhob(j)/rhoa(j)
        if (x > 2.0) x=2.0
        if (x < 0.5) x=0.5
        eb(j)=(ea(j)*(1.0-gama2*(1.0-x)/x)-qp(j)/rhob(j)*(1.0-x))/ &
        (1.0+gama2*(1.0-x))
    2 end do
    return
    end subroutine energy


    subroutine toten

!  total energy checks


!     evaluation of the total energy terms

    common /dens/ rhob(200),rhoa(200),dm(200),dmt(200)
    common /ener/ eb(200),ea(200)
    common /incre/ dt,dta,dtb,dtc,dtd
    common /init/ jm,jmp1,id
    common /radial/ rb(200),ra(200),r2b(200),r2a(200)
    common /tim/ time,ttime,totout,tout
    common /tot/ ekinb,ekina,eintb,einta,etotb,etota,error
    common /vel/ ua(200),fj(200)
    ekinb=0.0
    j=1
    do 1 j1=2,jmp1
        ekinb=ekinb+(dm(j)+dm(j1))*ua(j1)**2
        j=j1
    1 end do
    ekinb=0.25*ekinb
    eintb=0.0
    do 2 j=1,jm
        eintb=eintb+dm(j)*eb(j)
    2 end do

!     error is estimated from difference between energy terms before and
!     time step

    etotb=0.5*(eintb+einta)+ekinb
    error=etotb-etota
    if (etotb > 1.0e-10) error=error/etotb
    return
    end subroutine toten


    subroutine update

!  updates the results of run


!     the current values are updated at the end of the run

    common /dens/ rhob(200),rhoa(200),dm(200),dmt(200)
    common /ener/ eb(200),ea(200)
    common /incre/ dt,dta,dtb,dtc,dtd
    common /init/ jm,jmp1,id
    common /pres/ p(200)
    common /radial/ rb(200),ra(200),r2b(200),r2a(200)
    common /tim/ time,ttime,totout,tout
    common /tot/ ekinb,ekina,eintb,einta,etotb,etota,error
    common /bkdata/ gama,gama2
    do 2 j=1,jm
        ea(j)=eb(j)
        rhoa(j)=rhob(j)
        p(j)=rhoa(j)*gama*ea(j)
    2 end do
    do 3 j=1,jmp1
        r2a(j)=r2b(j)
        ra(j)=rb(j)
    3 end do
    einta=eintb
    ekina=ekinb
    time=time+dt
    return
    end subroutine update


    block data
    common /artvis/ qp(200)
    common /dens/ rhob(200),rhoa(200),dm(200),dmt(200)
    common /ener/ eb(200),ea(200)
    common /incre/ dt,dta,dtb,dtc,dtd
    common /init/ jm,jmp1,id
    common /pres/ p(200)
    common /radial/ rb(200),ra(200),r2b(200),r2a(200)
    common /tim/ time,ttime,totout,tout
    common /tot/ ekinb,ekina,eintb,einta,etotb,etota,error
    common /vel/ ua(200),fj(200)
    common /bkit/ it,iout
    data etota,einta,ekina,time,totout,error,ra(1),rb(1)/8*0.0/, &
    ua,fj,p,qp,ea,eb/1200*0.0/,rhoa,rhob,r2a,r2b/800*1.0/, &
    it/0/,dtc,dtd,dta/3*0.0/,dm(1)/0.0/
    end block data
