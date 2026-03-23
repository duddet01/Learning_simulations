program laplace
implicit none

integer, parameter :: lx=34, ly=34

real(8) :: old_temp(lx,ly), temp(lx,ly)
real(8) :: tol
integer :: i, j, ii, jj, counter, test

! ---------------- PARAMETERS ----------------
tol = 1.0d-4

! ---------------- INITIALIZATION ----------------
old_temp = 0.0d0

! ---------------- BOUNDARY CONDITIONS ----------------

! Left boundary (x=1)
do j=1,ly
    old_temp(1,j) = 3.7d0
end do

! Right boundary (x=34)
do j=1,ly
    old_temp(lx,j) = 0.4d0
end do

! Bottom boundary (y=1)
do i=1,lx
    old_temp(i,1) = 3.7d0 + (dble(i-1)/33.0d0)*(0.4d0 - 3.7d0)
end do

! Top boundary (y=34)
do i=1,lx
    old_temp(i,ly) = 3.7d0 + (dble(i-1)/33.0d0)*(0.4d0 - 3.7d0)
end do

temp = old_temp

! ---------------- ITERATION ----------------

counter = 0

do
    counter = counter + 1
    test = 0

    do j=2,ly-1
        do i=2,lx-1
            temp(i,j) = 0.25d0*( old_temp(i-1,j) + old_temp(i+1,j) + old_temp(i,j-1) + old_temp(i,j+1) )
        end do
    end do

    ! convergence check
    do j=2,ly-1
        do i=2,lx-1
            if (abs(temp(i,j) - old_temp(i,j)) > tol) test = 1
        end do
    end do

    if (test == 0) exit

    old_temp = temp

end do

print*, 'Iterations = ', counter
print*, 'Temperature at (20,20) = ', temp(20,20)

! ---------------- SAVE FILE ----------------

open(10,file='laplace.dat')

do i=1,lx
    do j=1,ly
        write(10,*) i, j, temp(i,j)
    end do
end do

close(10)

end program laplace