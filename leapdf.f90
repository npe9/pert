program leapdf
  ! *******************************************************************
  !  this program, which is a derivative of heatex, solves the problem
  !  either by the leapfrog or the dufort-frankel method.   the routine
  !  cycle has two statements labelled 70.   the first of these is for
  !  the leapfrog method, the second for dufort-frankel.   comment out
  !  the one not required.
  ! *******************************************************************
  !  where graphical output is requested up to 6 data files are produced
  !  'heatexn.dat' with n = 1 to 6 containing values of (x, temp).  if a
  !  7th output is requested then it will overwrite the first and so on.
  !  the times for the outputs are printed.
  !  ******************************************************************
  !  this uses the leapfrog or dufort-frankel method to solve the
  !  problem of the temperature variation in a lagged bar with boundary
  !  conditions giving the temperature at the ends of the bar either
  !  fixed or as some function of time.   the boundary conditions for
  !  the left and right-hand ends of the bar are given by function
  !  subprograms blh and brh respectively.  the following quantities
  !  are used:
  !           blen    the length of the bar
  !           cappa   the thermal conductivity
  !           c       the specific heat capacity
  !           ro      the density
  !  for the calculation the bar is divided into n segments and the
  !  the temperatures are output at m+1 points including the two ends.
  real :: temp(0:50),xx(0:50),tempx(2,49)
  character ans*1
  common temp,ntim,time,delt,xx,it,r,n,tmin,tmax,blen,tempx
  !  the following statement functions give the temperatures at the
  !  ends of the bar as functions of time
  data cappa,c,ro/200.,1000.,2700/
  open(unit=9,file='lpt1')
  blen=1.0
  time=0
  ng=0
  nw=0
20 write(6,'('' the following data [in si units] are being used'')')
  write(6,100)blen
100 format(19h 1. length of bar  ,f6.2)
  write(6,120)cappa
120 format(26h 2. thermal conductivity  ,f7.1)
  write(6,140)c
140 format(28h 3. specific heat capacity  ,f8.1)
  write(6,160)ro
160 format(13h 4. density  ,f8.1)
3 write(6,'('' do you want to change any of these? [y/n]'')')
  read(5,50)ans
50 format(a1)
  select case (ans)
  case ('Y','y')
     write(6,'('' input n, the number of segments in the bar'')')
     write(6,'(''for calculation.'')')
     read(5,*)n
     delx=blen/n
     do i=0,n
        xx(i)=delx*i
     end do
  case ('N','n')
     write(6,'('' input the number of an item to be changed'')')
     read(5,*)ni
     select case (ni)
     case (1)
        write(6,'('' input new value for length of bar'')')
        read(5,*)blen
     case (2)
        write(6,'('' input new thermal conductivity '')')
        read(5,*)cappa
     case (3)
        write(6,'('' input new specific heat capacity '')')
        read(5,*)c
     case (4)
        write(6,'('' input new density '')')
        read(5,*)ro
     end select
  end select
  do
     write(6,'('' you now have a choice of fixing either the '')')
     write(6,'('' time interval "t" or "r"=cappa*t/ro*c*dx**2'')')
     write(6,'('' input the choice [t/r]'')')
     read(5,50)ans
     select case (ans)
     case ('T', 't')
        write(6,'('' read in value of "t" '')')
        read(5,*)delt
        r=cappa*delt/delx**2/c/ro
        write(6,'('' '')')
        write(6,600)r
600     format(5h r = ,f8.4)
     case ('R', 'r')
        write(6,'('' read in the value of r '')')
        read(5,*)r
        delt=r*c*ro*delx**2/cappa
        write(6,'('' '')')
        write(6,620)delt
620     format(8h delt = ,f9.4)
        write(6,'(''  '')')
     end select
  end do
  write(6,'('' read in n-1 initial values of temperature'')')
  write(6,'('' at internal points of the bar.'')')
  write(6,'(''  '')')
  do i=1,n-1
     write(6,200)i
200  format(14h read in temp[,i2,1h])
     read(5,*)temp(i)
  end do
  temp(0)=blh(time)
  temp(n)=brh(time)
  !  the initial state of the bar is now fixed
  do
     write(6,'('' graphical file or printed output? [g/p]'')')
     read(5,50)ans
     select case (ans)
     case ('G','g')
        write(6,'('' how many timesteps before output?'')')
        read(5,*)ntim
        call cycle
        open(unit=21,file='heatex1.dat')
        open(unit=22,file='heatex2.dat')
        open(unit=23,file='heatex3.dat')
        open(unit=24,file='heatex4.dat')
        open(unit=25,file='heatex5.dat')
        open(unit=26,file='heatex6.dat')
74      ng=mod(ng,6)
        nw=ng+21
        ng=ng+1
        do i=0,n
           write(nw,*)xx(i),temp(i)
        end do
     case ('P', 'p')
        write(6,'(''input m for output, which must be a factor of n'')')
        write(6,'(''the m+1 temp values, including bar ends,are'')')
        write(6,'(''equally spaced along the bar'')')
        read(5,*)m
        k=n/m
        if(time < 1.0e-15)goto 65
        !  printed output section complete
        !  this section is for graphical output
     case default
        exit
     end select
     do
66      write(6,'(''continue the calculation? [y/n]'')')
        read(5,50)ans
        select case(ans)
        case ('Y', 'y')
           goto 65
        case ('N', 'n')
           goto 500
        end select
65      write(6,'('' how many timesteps before output?'')')
        read(5,*)ntim
        call cycle
        write(9,300)time
300     format(8h time = ,f7.2)
        if(nw /= 0) goto 74
        write(9,320)(temp(i),i=0,n,k)
320     format(6f7.1)
     end do
  end do
  !  this section is for printed output
500 stop
end program leapdf



real function blh(x)
  blh=300
  return
end function blh

real function brh(x)
  brh=400
  return
end function brh




subroutine cycle
  !  converts to the leapfrog method with first statement 70
  !  converts to dufont-frankel with second statement 70
  real :: temp(0:50),tempx(2,49),xx(0:50)
  common temp,ntim,time,delt,xx,it,r,n,tmin,tmax,blen,tempx
  do 61 i=1,ntim
     time=time+delt
     do 62 j=1,n-1
        if(time > 1.1*delt)goto 70
        tempx(2,j)=r*(temp(j-1)+temp(j+1))+(1-2*r)*temp(j)
        goto 62
        !   70 tempx(2,j)=tempx(1,j)+2*r*(temp(j+1)+temp(j-1)-2*temp(j))
70      tempx(2,j)=((1-r)*tempx(1,j)+r*(temp(j+1)+temp(j-1)))/(1+r)
62   end do
     do 63 j=1,n-1
        tempx(1,j)=temp(j)
        temp(j)=tempx(2,j)
63   end do
     temp(0)=blh(time)
     temp(n)=brh(time)
61 end do
  return
end subroutine cycle
