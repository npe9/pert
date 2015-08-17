    program drunkard
!  this programme carries out a modified random-walk calculation in
!  which the equal-length steps are only along the x and y directions.
!  to deal with a range of problems the probabilities are defined as
!  for (1) the direction of the last step,
!      (2) the reverse of the direction of the last step
!  and (3) at right angles to the direction of the last step
!          where right and left directions have equal probability.
!  the associated probabilities are p(1), p(2) and p(3).   only p(1)
!  and p(2) need to be defined and they can be modified by the data
!  statement.

!  for this programme a random number generator is used which
!  is derived from "numerical recipes" by w.h.press,b.p. flannery,
!  s.a. teulkolsky and w.t. vetterling (c.u.p.)

!  for each value of n specified the root-mean square distance of
!  travel is estimated from 1000 separate random walks.
    dimension d(1000),p(2),cum(4)
    data p,idum/0.25,0.25,-255/
    cum(1)=p(1)
    cum(2)=p(1)+p(2)
    cum(4)=1.0
    cum(3)=(cum(2)+cum(4))/2.0
    write(6,'('' input the number of steps, n '')')
    read(5,*)n
    do 1 i=1,1000
    !  the first step is taken arbitrarily in the +x direction
        id=0
    !  id=0 if step is in the positive x direction
    !  id=2 if step is in the negative x direction
    !  id=1 if step is in the positive y direction
    !  id=3 if step is in the negative y direction.
        sumx=1
        sumy=0
        do 2 j=1,n-1
            r=ran2(idum)
            iadd=3
            if(cum(3) > r)iadd=1
            if(cum(2) > r)iadd=2
            if(cum(1) > r)iadd=0
            id=mod(id+iadd,4)
            idp=id+1
            goto(3,4,5,6)idp
        ! positive x direction
            3 sumx=sumx+1
            goto 2
        ! positive y direction
            4 sumy=sumy+1
            goto 2
        ! negative x direction
            5 sumx=sumx-1
            goto 2
        ! negative y direction
            6 sumy=sumy-1
        2 end do
        d(i)=sumx**2+sumy**2
    1 end do
    tot=0
    do 10 i=1,1000
        tot=tot+d(i)
    10 end do
    tot=sqrt(tot/1000.0)
    write(6,100)n,tot
    100 format(25h the average distance for,i5,9h steps is,f9.2)
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

          
