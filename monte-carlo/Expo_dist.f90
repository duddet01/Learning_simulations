program exponential_rng
    implicit none
    integer, parameter :: N = 100000
    integer :: i
    real(8) :: r, x

    call random_seed()
     open(unit=15, file="expo.dat", status="replace")
    
    do i = 1, N
        call random_number(r)
        x = -0.5d0 * log(r)
        write(15,*) x

        
    end do
end program exponential_rng
