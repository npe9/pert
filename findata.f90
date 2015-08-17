    program findata
!  this program prepares a data file for the program finelem2 which
!  gives a finite-element solution of the problem of a
!  heated plate insulated on its top and bottom surfaces with all
!  elements triangular.   differential boundary conditions are
!  specified by
!           cappa*(d(theta)/d(n)) = -m*theta + s
!  where n is along the direction of the normal to the surface.
!  values of m and s are input and so are estimated before the
!  program is used.  some nodes, usually at boundaries, can be
!  specified as having fixed temperatures.   there is provision for
!  an extended source of heat which is given by a function statement
!  q(x,y).   if such a source is present then the user must ensure
!  that it gives the required form of heating.   there is also
!  provision for up to 9 point sources which must be situated at nodes.
    dimension nde(100,3),nedge(50,2),cay(50),ss(50),np(50),sps(50)
    open(unit=11,file='finelem.dat')
!  input conductivity of the plate
    write(6,'(''enter conductivity'')')
    read(5,*)q
    write(11,*)q
    iout=0
    write(6,'(''input number of elements and number of nodes'')')
    read(5,*)k,l
    write(11,*)k,l
    write(6,'(''enter node numbers associated with each element'')')
    write(6,'(''  '')')
    do 1 i=1,k
        write(6,100)i
        100 format(38h enter 3 node numbers defining element,i2)
        read(5,*)(nde(i,j),j=1,3)
        write(11,*)(nde(i,j),j=1,3)
    1 end do
    write(6,'(''enter coordinates of nodes'')')
    do 2 i=1,l
        write(6,200)i,i
        200 format(9h enter x[,i2,8h] and y[,i2,1h])
        read(5,*)q,r
        write(11,*)q,r
    2 end do
    write(6,'(''enter number of nodes with fixed temperatures'')')
    read(5,*)nft
    write(11,*)nft
    if(nft == 0)goto 5
    do 4 i=1,nft
        write(6,'(''read in node number, fixed temperature'')')
        read(5,*)k,q
        write(11,*)k,q
    4 end do
    5 write(6,'(''enter the number of point sources'')')
    read(5,*)k
    write(11,*)k
    if(k == 0)goto 13
    do 14 i=1,k
        write(6,'(''give node and strength of point source [w/m**3]'')')
        read(5,*)np(i),sps(i)
        write(11,*)np(i),sps(i)
    14 end do
    13 write(6,'(''how many edges of elements have boundary '')')
    write(6,'(''conditions of the form '')')
    write(6,'('' d(theta)/d(n) = -m*theta + s ?'')')
    read(5,*)nbc
    write(11,*)nbc
    if(nbc == 0)goto 16
    do 15 i=1,nbc
        write(6,'(''read in the edge as the two end nodes and the '')')
        write(6,'(''values of m and s'')')
        read(5,*)nedge(i,1),nedge(i,2),cay(i),ss(i)
        write(11,*)nedge(i,1),nedge(i,2),cay(i),ss(i)
    15 end do
    16 close(11)
    stop
    end program


