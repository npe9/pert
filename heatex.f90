program heatex
  ! *******************************************************************
  !  this program differs from heatcrni which can solve the same
  !  problem using the implicit crank-nicholson method only in the
  !  subroutine 'cycle.
  ! *******************************************************************
  !  where graphical output is requested up to 6 data files are produced
  !  'heatexn.dat' with n = 1 to 6 containing values of (x, temp).  if a
  !  7th output is requested then it will overwrite the first and so on.
  !  the times for the outputs are printed.
  !  ******************************************************************
  !  this uses the explicit method to solve the problem of the
  !  temperature variation in a lagged bar with boundary conditions
  !  giving the temperature at the ends of the bar either fixed or
  !  as some function of time.   the boundary conditions for the
  !  left and right-hand ends of the bar are given by function
  !  subprograms blh and brh respectively.  the following quantities
  !  are used:
  !           blen    the length of the bar
  !           cappa   the thermal conductivity
  !           c       the specific heat capacity
  !           ro      the density
  !  for the calculation the bar is divided into n segments and the
  !  the temperatures are output at m+1 points including the two ends.
  real :: temp(0:50),xx(0:50)
  character ans*1
  common temp,ntim,time,delt,xx,it,r,n,tmin,tmax,blen
  !  the following statement functions give the temperatures at the
  !  ends of the bar as functions of time
  data cappa,c,ro/200.,1000.,2700/
  open(unit=9,file='lpt1')
  blen=1.0
  time=0
  ng=0
  nw=0
  outer: do
     write(6,'('' the following data [in si units] are being used'')')
     write(6,100) blen
100  format(19h 1. length of bar  ,f6.2)
     write(6,120) cappa
120  format(26h 2. thermal conductivity  ,f7.1)
     write(6,140) c
140  format(28h 3. specific heat capacity  ,f8.1)
     write(6,160) ro
160  format(13h 4. density  ,f8.1)
     inner: do
        write(6,'('' do you want to change any of these? [y/n]'')')
        read(5,50) ans
50      format(a1)
        select case (ans)
        case ('N', 'n')
           exit outer
        case ('Y', 'y')
           exit inner
        case default
           cycle
           end select
        end do inner
        write(6,'('' input the number of an item to be changed'')')
        read(5,*)ni
        goto(11,12,13,14)ni
11      write(6,'('' input new value for length of bar'')')
        read(5,*)blen
        cycle
12      write(6,'('' input new thermal conductivity '')')
        read(5,*)cappa
        cycle
13      write(6,'('' input new specific heat capacity '')')
        read(5,*) c
        cycle
14      write(6,'('' input new density '')')
        read(5,*) ro
     end do  outer
     write(6,'('' input n, the number of segments in the bar'')')
     write(6,'(''for calculation.'')')
     read(5,*)n
     delx=blen/n
     do i=0,n
        xx(i)=delx*i
     end do
     do
        write(6,'('' you now have a choice of fixing either the '')')
        write(6,'('' time interval "t" or "r"=cappa*t/ro*c*dx**2'')')
        write(6,'('' input the choice [t/r]'')')
        read(5,50) ans
        select case (ans)
        case ('T', 't')
           write(6,'('' read in value of "t" '')')
           read(5,*)delt
           r=cappa*delt/delx**2/c/ro
           write(6,'('' '')')
           write(6,600)r
600        format(5h r = ,f8.4)
        case ('R', 'r')
           write(6,'('' read in the value of r '')')
           read(5,*)r
           delt=r*c*ro*delx**2/cappa
           write(6,'('' '')')
           write(6,620)delt
620        format(8h delt = ,f9.4)
           write(6,'(''  '')')

        case default
           cycle
        end select
     end do
30   write(6,'('' read in n-1 initial values of temperature'')')
     write(6,'('' at internal points of the bar.'')')
     write(6,'(''  '')')
     do  i=1,n-1
        write(6,200)i
200     format(14h read in temp[,i2,1h])
        read(5,*)temp(i)
     end do
     temp(0)=blh(time)
     temp(n)=brh(time)
     !  the initial state of the bar is now fixed
     do
        write(6,'('' graphical file or printed output? [g/p]'')')
        read(5,50) ans
        select case(ans)
        case ('G', 'g')
           goto 60
        case ('P', 'p')
           exit
        case default
           cycle
        end select
     end do
     !  this section is for printed output
     write(6,'(''input m for output, which must be a factor of n'')')
     write(6,'(''the m+1 temp values, including bar ends,are'')')
     write(6,'(''equally spaced along the bar'')')
     read(5,*) m
     k=n/m
     if(time < 1.0e-15) goto 65
     do
        write(6,'(''continue the calculation? [y/n]'')')
        read(5,50)ans
        select case (ans)
        case ('Y', 'y')
        case ('N', 'n')
           exit
        case default
           cycle
        end select
65      write(6,'('' how many timesteps before output?'')')
        read(5,*) ntim
        call cycle
        write(9,300) time
300     format(8h time = ,f7.2)
        if(nw == 0) then
           write(9,320)(temp(i),i=0,n,k)
320        format(6f7.1)
           cycle
           !  printed output section complete
           !  this section is for graphical output
60         write(6,'('' how many timesteps before output?'')')
           read(5,*) ntim
           call cycle
           open(unit=21,file='heatex1.dat')
           open(unit=22,file='heatex2.dat')
           open(unit=23,file='heatex3.dat')
           open(unit=24,file='heatex4.dat')
           open(unit=25,file='heatex5.dat')
           open(unit=26,file='heatex6.dat')
        end if
        ng=mod(ng,6)
        nw=ng+21
        ng=ng+1
        do i=0,n
           write(nw,*)xx(i),temp(i)
        end do
     end do
     stop
   end program heatex

   real function blh(x)
     blh=300
     return
   end function blh

   real function brh(x)
     brh=400
     return
   end function brh

   subroutine cycle
     real :: temp(0:50),tempx(49),xx(0:50)
     common temp,ntim,time,delt,xx,it,r,n,tmin,tmax,blen
     do  i=1,ntim
        time=time+delt
        do  j=1,n-1
           tempx(j)=r*(temp(j-1)+temp(j+1))+(1-2*r)*temp(j)
        end do
        do  j=1,n-1
           temp(j)=tempx(j)
        end do
        temp(0)=blh(time)
        temp(n)=brh(time)
     end do
     return
   end subroutine cycle
