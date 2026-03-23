program rk4_ring
implicit none

integer, parameter :: nop = 50
integer :: i, kk, kkp1, kkn1, iter

real(8) :: dt, dtby2, t
real(8), parameter :: km = 1.0d0

real(8) :: yim(nop), vim(nop)
real(8) :: y1(nop), y2(nop), y3(nop)
real(8) :: v1(nop), v2(nop), v3(nop)

real(8) :: f0(nop), f1(nop), f2(nop), f3(nop)
real(8) :: f0v(nop), f1v(nop), f2v(nop), f3v(nop)

! ---------------- PARAMETERS ----------------
iter = 2000
dt = 0.02d0
dtby2 = dt / 2.0d0
t = 0.0d0

! ---------------- INITIAL CONDITIONS ----------------
yim = 0.0d0
vim = 0.0d0

yim(1) = 0.8d0
yim(26) = 0.8d0
open(99,file='part1.dat')

! ---------------- TIME LOOP ----------------
do i = 1, iter

    t = t + dt

    ! k1
    do kk = 1, nop
        kkp1 = kk + 1
        kkn1 = kk - 1
        if (kk == 1) kkn1 = nop
        if (kk == nop) kkp1 = 1

        f0(kk) = vim(kk)
        f0v(kk) = km*(yim(kkp1) + yim(kkn1) - 2.0d0*yim(kk))
    end do

    ! k2
    do kk = 1, nop
        y1(kk) = yim(kk) + dtby2*f0(kk)
        v1(kk) = vim(kk) + dtby2*f0v(kk)
    end do

    do kk = 1, nop
        kkp1 = kk + 1  
        kkn1 = kk - 1
        if (kk == 1) kkn1 = nop
        if (kk == nop) kkp1 = 1

        f1(kk) = v1(kk)
        f1v(kk) = km*(y1(kkp1) + y1(kkn1) - 2.0d0*y1(kk))
    end do

    ! k3
    do kk = 1, nop
        y2(kk) = yim(kk) + dtby2*f1(kk)
        v2(kk) = vim(kk) + dtby2*f1v(kk)
    end do

    do kk = 1, nop
        kkp1 = kk + 1
        kkn1 = kk - 1
        if (kk == 1) kkn1 = nop
        if (kk == nop) kkp1 = 1

        f2(kk) = v2(kk)
        f2v(kk) = km*(y2(kkp1) + y2(kkn1) - 2.0d0*y2(kk))
    end do

    ! k4
    do kk = 1, nop
        y3(kk) = yim(kk) + dt*f2(kk)
        v3(kk) = vim(kk) + dt*f2v(kk)
    end do

    do kk = 1, nop
        kkp1 = kk + 1
        kkn1 = kk - 1
        if (kk == 1) kkn1 = nop
        if (kk == nop) kkp1 = 1

        f3(kk) = v3(kk)
        f3v(kk) = km*(y3(kkp1) + y3(kkn1) - 2.0d0*y3(kk))
    end do

    ! update
    do kk = 1, nop
        yim(kk) = yim(kk) + dt*(f0(kk) + 2.0d0*f1(kk) + 2.0d0*f2(kk) + f3(kk))/6.0d0
        vim(kk) = vim(kk) + dt*(f0v(kk) + 2.0d0*f1v(kk) + 2.0d0*f2v(kk) + f3v(kk))/6.0d0
    end do
    
    write(99,*) t , yim(1) ,  yim(26) , yim(40)
end do

! ---------------- FINAL OUTPUT ----------------
write(*,*) 'y1 at t = 40 is:', yim(1) 
close(99)
end program rk4_ring