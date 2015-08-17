    program markov3
!  this program requires long integers
    dimension ns(3),p(3,3),pp(3)
    idum=-567

!  reading in rows of stochastic matrix
    do 1 i=1,3
        write(6,60)i
        60 format(10h input row,i2,10h of matrix)
        read(5,*)(p(i,j),j=1,3)
        s=p(i,1)+p(i,2)+p(i,3)
        if((s-1.0)**2 < 1.0e-8)goto 1
        write(6,'('' sum of probabilities not equal to 1. '')')
        write(6,'('' program aborted. '')')
        goto 100
    1 end do
!  clear table
    do 4 i=1,3
        ns(i)=0
    4 end do

!  fix an initial variable
    n=1

    do 2 i=1,1000000
        x=ran1(idum)
    !  select next variable if present one is the first
        if(n == 1)then
            if(x < p(1,1)+p(1,2))goto 13
            n1=3
            goto 50
            13 if(x < p(1,1))goto 14
            n1=2
            goto 50
            14 n1=1
            goto 50
        endif
    
    !  select next variable if present one is the second
        if(n == 2)then
            if(x < p(2,1)+p(2,2))goto 15
            n1=3
            goto 50
            15 if(x < p(2,1))goto 16
            n1=2
            goto 50
            16 n1=1
            goto 50
        endif
    
    !  select next variable if present one is the third
        if(n == 3)then
            if(x < p(3,1)+p(3,2))goto 17
            n1=3
            goto 50
            17 if(x < p(3,1))goto 18
            n1=2
            goto 50
            18 n1=1
        endif
    
        50 n=n1
        ns(n)=ns(n)+1
    2 end do
!  output result on screen
    write(6,*)ns
    sumns=float(ns(1)+ns(2)+ns(3))
    do 5 i=1,3
        pp(i)=ns(i)/sumns
    5 end do
    write(6,200)pp
    200 format(22h the probabilities are, 3f8.4)
    100 stop
    end program


    function ran1(idum)
    dimension r(97)
    parameter (m1=259200,ia1=7141,ic1=54773,rm1=3.8580247e-6)
    parameter (m2=134456,ia2=8121,ic2=28411,rm2=7.4373773e-6)
    parameter (m3=243000,ia3=4561,ic3=51349)
    data iff /0/
    if (idum < 0 .or. iff == 0) then
        iff=1
        ix1=mod(ic1-idum,m1)
        ix1=mod(ia1*ix1+ic1,m1)
        ix2=mod(ix1,m2)
        ix1=mod(ia1*ix1+ic1,m1)
        ix3=mod(ix1,m3)
        do 11 j=1,97
            ix1=mod(ia1*ix1+ic1,m1)
            ix2=mod(ia2*ix2+ic2,m2)
            r(j)=(float(ix1)+float(ix2)*rm2)*rm1
        11 end do
        idum=1
    endif
    ix1=mod(ia1*ix1+ic1,m1)
    ix2=mod(ia2*ix2+ic2,m2)
    ix3=mod(ia3*ix3+ic3,m3)
    j=1+(97*ix3)/m3
    if(j > 97 .or. j < 1)pause
    ran1=r(j)
    r(j)=(float(ix1)+float(ix2)*rm2)*rm1
    return
    end function ran1

