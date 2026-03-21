program gaussian_rng
    implicit none
    integer, parameter :: N = 100000
    integer :: i
    real(8) :: r1, r2, z, x
    real(8), parameter :: pi = 4.0d0 * atan(1.0d0)

    call random_seed()
    open(unit=16, file="gauss.dat", status="replace")
    
    do i = 1, N
        call random_number(r1)
        call random_number(r2)

        z = sqrt(-2.0d0 * log(r1)) * cos(2.0d0 * pi * r2)
        x = 2.0d0 * z
        write(16,*) x
        
    end do
end program gaussian_rng
