program laplace_neumann
implicit none

integer, parameter :: lx=34, ly=34

real(8) :: old_temp(lx,ly), temp(lx,ly)
real(8) :: dx, dy, prefactor, shift_temp
real(8) :: A(ly), B(ly), C(lx), D(lx)

integer :: i, j, ii, jj, counter, test

! ---------------- PARAMETERS ----------------

dx = 1.0d0
dy = dx

A = -70.0d0
B = -40.0d0
C =  20.0d0
D = -10.0d0

old_temp = 0.0d0

counter = 0
prefactor = (0.5d0*dx*dx*dy*dy)/(dx*dx + dy*dy)

! ---------------- ITERATION ----------------

do
    counter = counter + 1
    test = 0

    ! ----------- LEFT & RIGHT (Neumann BC) -----------
    do j=2,ly-1
        temp(1,j) = 0.25d0*( 2.0d0*old_temp(2,j) - 2.0d0*dx*A(j) + &
                             old_temp(1,j+1) + old_temp(1,j-1) )

        temp(lx,j) = 0.25d0*( 2.0d0*old_temp(lx-1,j) + 2.0d0*dx*B(j) + &
                              old_temp(lx,j+1) + old_temp(lx,j-1) )
    end do

    ! ----------- BOTTOM & TOP -----------
    do i=2,lx-1
        temp(i,1) = 0.25d0*( old_temp(i+1,1) + old_temp(i-1,1) + &
                             2.0d0*old_temp(i,2) - 2.0d0*dx*C(i) )

        temp(i,ly) = 0.25d0*( old_temp(i+1,ly) + old_temp(i-1,ly) + &
                              2.0d0*old_temp(i,ly-1) + 2.0d0*dx*D(i) )
    end do

    ! ----------- CORNERS -----------
    temp(1,1) = 0.5d0*( old_temp(1,2) - dx*C(1) + old_temp(2,1) - dx*A(1) )

    temp(1,ly) = 0.5d0*( old_temp(1,ly-1) + dx*D(1) + old_temp(2,ly) - dx*A(ly) )

    temp(lx,1) = 0.5d0*( old_temp(lx-1,1) + dx*B(1) + old_temp(lx,2) - dx*C(lx) )

    temp(lx,ly) = 0.5d0*( old_temp(lx-1,ly) + dx*B(ly) + old_temp(lx,ly-1) + dx*D(lx) )

    ! ----------- INTERIOR -----------
    do jj=2,ly-1
        do ii=2,lx-1
            temp(ii,jj) = 0.25d0*( old_temp(ii-1,jj) + old_temp(ii+1,jj) + &
                                   old_temp(ii,jj-1) + old_temp(ii,jj+1) )
        end do
    end do

    ! ----------- CONVERGENCE CHECK -----------
    do jj=1,ly
        do ii=1,lx
            if (abs(temp(ii,jj) - old_temp(ii,jj)) > 1.0d-3) then
                test = 1
            end if
        end do
    end do

    if (test == 0) exit

    ! ----------- SHIFT (FIX GAUGE) -----------
    shift_temp = 2000.0d0 - temp(1,1)
    temp = temp + shift_temp

    old_temp = temp

end do

print*, 'Iterations = ', counter
print*, 'Temperature at (10,10) = ', temp(10,10)

! ---------------- OUTPUT ----------------
open(10,file='neumann.dat')

do i=1,lx
    do j=1,ly
        write(10,*) i, j, temp(i,j)
    end do
end do

close(10)

end program laplace_neumann