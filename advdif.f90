    program advdif
!  this programme solves the advection-diffusion equation in one
!  dimension with advection alone or diffusion alone if required.
!  for advection alone the timestep comes from the user-input
!  courant number, velocity and space interval.   for diffusion
!  alone it is given by the dt = r*(dx)**2/d.   with both advection
!  and diffusion the timestep from the courant number is given
!  priority but if it makes r > 0.5 then the timestep is adjusted
!  to give r = 0.5.   periodic boundary conditions are imposed so
!  that a disturbance passing through one barrier reappears at the
!  other.
!  several finite-difference formulae are available and the user
!  must comment out the ones not required.
!  *****************************************************************

    dimension q(0:999),tab(0:999,0:1000),x(0:999)
    character ans*1
    write(6,'(''input velocity, diffusion coefficient and'')')
    write(6,'(''the grid spacing dx'')')
    read(5,*)u,d,dx
    open(unit=11,file='q1.dat')
    open(unit=12,file='q2.dat')
    open(unit=13,file='q3.dat')
    open(unit=14,file='q4.dat')
    open(unit=15,file='q5.dat')
    open(unit=16,file='q6.dat')
!  set values in x array
    do 16 i=0,999
        x(i)=i*dx
    16 end do
!  test if advection is present
    if(u < 1.0e-20)goto 1
    write(6,'(''input the courant number'')')
    read(5,*)c
!  calculate dt
    dt=c*dx/u
!  calculate r
    r=d*dt/dx**2
    if(r > 0.5)dt=0.5*dx**2/d
    goto 2
    1 write(6,'(''input r=d*dt/dx**2'')')
    read(5,*)r
!  calculate dt
    dt=r*dx**2/d
    2 write(6,100)dt
    100 format(6h dt = ,f8.5)
!  clear q and tab tables
    do 11 i=0,999
        q(i)=0
        do 11 j=0,1000
            tab(i,j)=0
    11 end do
!  read in initial distribution of q
    write(6,'(''how many q values to be read in?'')')
    read(5,*)n
    do 12 l=1,n
        write(6,'(''read in i and q[i]'')')
        read(5,*)i,q(i)
    12 end do
!  place initial distribution in tab
    do 15 i=0,999
        tab(i,0)=q(i)
    15 end do
!  set the number of timesteps for the simulation
    write(6,'(''input the number of simulation timesteps'')')
    read(5,*)ntime
    do 5 j=1,ntime
    !  integrate the equations storing the new q values in tab
        do 21 i=0,999
            ip=mod(i+1,1000)
            im=mod(i+999,1000)
        !  pure diffusion time f-d, space c-d (explicit method)
        !      tab(i,j)=(1-2*r)*q(i)+r*(q(ip)+q(im))
        !  equation (7.6) advection only.  time f-d, space c-d.
        !      tab(i,j)=q(i)-c*(q(ip)-q(im))/2
        !  equation (7.7) advection only.  time f-d, space b-d
        !      tab(i,j)=(1-c)*q(i)+c*q(im)
        !  equation (7.8) advection with space b-d, diffusion with space c-d
        !        tab(i,j)=(1-2*r-c)*q(i)+r*q(ip)+(r+c)*q(im)
        !  equation (7.21) both advection and diffusion with space c-d
        !      tab(i,j)=(1-2*r)*q(i)+(r-c/2)*q(ip)+(r+c/2)*q(im)
        !  the lax-wendroff advection equation (7.29)
        !      tab(i,j)=c*(1+c)/2*q(im)+(1-c*c)*q(i)-c*(1-c)/2*q(ip)
        !  quickest method for advection+diffusion
            imm=mod(i+998,1000)
            ipp=mod(i+2,1000)
            c1=0
            c2=r*(1-c)-c*(c-1)*(c-2)/6
            c3=1-r*(2-3*c)+c*(c*c-2*c-1)/2
            c4=r*(1-3*c)-c*(c+1)*(c-2)/2
            c5=r*c+c*(c*c-1)/6
            tab(i,j)=c1*q(ipp)+c2*q(ip)+c3*q(i)+c4*q(im)+c5*q(imm)
        21 end do
    !  update q values
        do 22 i=0,999
            q(i)=tab(i,j)
        22 end do
    5 end do
!  extract up to 6 sets of q values for graphical output
    30 write(6,'(''you can now set up files q1.dat to qk.dat'')')
    write(6,'(''with values of x and q for graphical output'')')
    write(6,'(''for k<=6.  do you want to do this? [y/n]'')')
    read(5,50)ans
    50 format(a1)
    if(ans == 'n' .or. ans == 'n')goto 200
    if(ans == 'y' .or. ans == 'y')goto 34
    goto 30
    34 write(6,'(''input the number of sets of q values you require'')')
    read(5,*)nset
    do 32 i=1,nset
        write(6,'(''read in timestep number of set [0 to 100]'')')
        read(5,*)j
        m=10+i
        do 36 ii=0,999
            write(m,*)x(ii),tab(ii,j)
        36 end do
    32 end do
    200 stop
    end program
