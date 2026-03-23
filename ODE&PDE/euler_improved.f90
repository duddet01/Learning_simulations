program euler
implicit none

integer, parameter :: dp = kind(1.0d0)

real(dp) :: y, y_new, fxy, fxy1, x, dx
integer :: i, niter

! initialization
x = 0.0d0
y = 0.0d0
dx = 0.001d0
niter = int(1.55d0/dx)

open (70,file="euler_mod.dat")

do i = 1, niter+1

    fxy = 1.0d0 + y*y
    y_new = y + fxy*dx

    fxy1 = 1.0d0 + y_new*y_new
    y = y + dx/(2.0d0)*(fxy + fxy1)

    write (70,*) i, dble(i)*dx, y

end do 

close(70)

end program euler