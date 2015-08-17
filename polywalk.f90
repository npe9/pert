    program polywalk
!  this programme finds the average distance between the ends of a
!  polymer chain confined to two dimensions under the condition
!  that the polymer is not allowed to cross itself.   the polymer
!  units correspond to the steps of a random walk and a modified
!  random walk is carried out where steps are only in the x or y
!  directions.  a step in a trial that produces a crossed chain is
!  repeated up to four times to find a non-crossing path but if it
!  does not do so then it is abandonned and the rms walk is calculated
!  from those that succeed.   the length of chain should be restricted
!  to 50.

!  for this programme a random number generator is used which
!  is derived from "numerical recipes" by w.h. press, b.p. flannery,
!  s.a. teulkolsky and w.t. vetterling (c.u.p.)

    dimension linktab(-40:40,-40:40),ipos(2)
    data idum/-122/
    write(6,'('' input the length of the chain '')')
    read(5,*)n
!  start the 1000 trials
    nsum=0
    sumd=0
    do 2 i=1,1000
    !  clear linktab
        do 1 k=-40,40
            do 1 l=-40,40
                linktab(k,l)=0
        1 end do
    !  the first step is taken arbitrarily in the +x direction
        linktab(0,0)=1
        linktab(1,0)=1
        ipos(1)=1
        ipos(2)=0
        do 3 j=1,n-1
            itry=0
            13 itry=itry+1
            if(itry > 4)goto 2
            r=ran2(idum)
            if(0.25 >= r)goto 4
            if(0.5 >= r)goto 5
            if(0.75 >= r)goto 6
        !  the step is along +y
            ipos(2)=ipos(2)+1
        !  test to see if a crossing is being made
            if(linktab(ipos(1),ipos(2)) == 1)goto 13
            linktab(ipos(1),ipos(2))=1
            goto 3
        !  the step is along -y
            6 ipos(2)=ipos(2)-1
        !  test to see if a crossing is being made
            if(linktab(ipos(1),ipos(2)) == 1)goto 13
            linktab(ipos(1),ipos(2))=1
            goto 3
        !  the step is along +x
            5 ipos(1)=ipos(1)+1
        !  test to see if a crossing is being made
            if(linktab(ipos(1),ipos(2)) == 1)goto 13
            linktab(ipos(1),ipos(2))=1
            goto 3
        !  the step is along -x
            4 ipos(1)=ipos(1)-1
        !  test to see if a crossing is being made
            if(linktab(ipos(1),ipos(2)) == 1)goto 13
            linktab(ipos(1),ipos(2))=1
        3 end do
        sumd=sumd+float(ipos(1))**2+float(ipos(2))**2
        nsum=nsum+1
    2 end do
    write(6,100)nsum
    100 format(35h the number of successful walks was,i4)
    rms=sqrt(sumd/nsum)
    write(6,200)rms
    200 format(20h the rms distance is,f9.2)
    stop
    end program
                      
    function ran2(idum)
    parameter(m=714025,ia=1366,ic=150889,rm=1./m)
    dimension ir(97)
    data iff/0/
    if(idum < 0 .or. iff == 0)then
        iff=1
        idum=mod(ic-idum,m)
        do 11 j=1,97
            idum=mod(ia*idum+ic,m)
            ir(j)=idum
        11 end do
        idum=mod(ia*idum+ic,m)
        iy=idum
    endif
    j=1+(97*iy)/m
    if(j > 97 .or. j < 1)pause
    iy=ir(j)
    ran2=iy*rm
    idum=mod(ia*idum+ic,m)
    ir(j)=idum
    return
    end function ran2

