program Q5

    implicit none

    integer ::  i , max , iter
    real(8) :: l , d ,r , midpnt ,theta ,hits ,pi_est,P 
    real(8), parameter :: pi_val = 4 * atan(1.0d0)


    open(unit=10, file="buffon.dat", status="replace")
    

    
    print *, "Reminder: Length of needle < Distance"

    print* , "Give the length of needle  : " 
    read(*,*) l

    print* , "Give the dist b/w lines  : " 
    read(*,*) d

    do iter=2,9

        max=10**iter
            
    
        hits=0


        do i=1,max
        call random_number(r)
        midpnt= r*d/2.0d0    !half of the dist between lines to sample mp
        call random_number(r)
        theta=r*pi_val/2.0d0

        if (midpnt - (l/2.0d0)*cos(theta) <0 )then
            hits=hits+1
        end if   



        end do   
        P=hits/max
        pi_est=2*l/(P*d)
        print * , "The value of pi is :" , pi_est

        write(10,*) pi_est


    

    end do    
        
    close(10)


     


end program Q5
