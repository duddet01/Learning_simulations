program brute_force_MC
    
    implicit none
    integer, parameter :: N = 10000000
    integer :: i
    real(8) :: r1, r2, ct1, ct2, phi1, phi2 ,r12
    real(8) :: f, sum, sum2, val, err
    real(8), parameter :: Rmax = 5.0d0
    real(8) :: volume
    real(8), parameter :: alpha = 2.0d0
    real(8), parameter :: pi = 4.0d0 * atan(1.0d0)
    real(8) :: A=8/pi

    call random_seed()
    sum = 0.0d0
    sum2 = 0.0d0

    volume = (Rmax**2) * (2.0d0)**2 * (2.0d0*pi)**2

    do i = 1, N
        call random_number(r1)
        call random_number(r2)
        call random_number(ct1)
        call random_number(ct2)
        call random_number(phi1)
        call random_number(phi2)

        r1 = Rmax * r1
        r2 = Rmax * r2
        ct1 = 2.0d0*ct1 - 1.0d0
        ct2 = 2.0d0*ct2 - 1.0d0
        phi1 = 2.0d0*pi*phi1
        phi2 = 2.0d0*pi*phi2

        f = A*A*r1*r1 * r2*r2 * exp(-2.0d0*alpha*(r1+r2)) / &
            r12(r1,r2,ct1,ct2,phi1,phi2)

        sum  = sum  + f
        sum2 = sum2 + f*f
    end do

    val = volume * sum / N
    err = volume * sqrt((sum2/N - (sum/N)**2)/N)

    print *, "Brute-force result =", val
    print *, "Brute-force error  =", err
end program brute_force_MC




real(8) function r12(r1, r2, ct1, ct2, phi1, phi2)
    implicit none
    real(8), intent(in) :: r1, r2, ct1, ct2, phi1, phi2
    real(8) :: cosb, st1, st2
    real(8), parameter :: alpha = 2.0d0
    real(8), parameter :: pi = 4.0d0 * atan(1.0d0)

    st1 = sqrt(1.0d0 - ct1*ct1)
    st2 = sqrt(1.0d0 - ct2*ct2)

    cosb = ct1*ct2 + st1*st2*cos(phi1 - phi2)
    r12 = sqrt(r1*r1 + r2*r2 - 2.0d0*r1*r2*cosb)
end function r12
