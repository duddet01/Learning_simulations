module sim_para
  implicit none

  integer, parameter :: n_part = 2197
  integer, parameter :: nsteps = 20000

  real(8), parameter :: dt = 0.005d0
  real(8), parameter :: mass = 1.0d0

  real(8), parameter :: lx = 20.0d0, ly = 20.0d0, lz = 20.0d0

  real(8), parameter :: sigma = 1.0d0, epsilon = 1.0d0 , kbT=1.0d0
  real(8), parameter :: rc = 2.5d0

  real(8) :: pos(3*n_part), vel(3*n_part), force(3*n_part)
  integer, parameter :: thermo_interval = 500

contains

  ! Periodic Boundary Condition
  real(8) function pbc(x, box)
    real(8), intent(in) :: x, box
    pbc = x - box * dnint(x / box)
  end function

end module

!==================================================================


program md
  use sim_para
  implicit none

  integer :: step

  ! --- open energy file ---
  open(10, file="md_with_thermostat.dat", status="unknown")
  write(10,'(A)') "# Step   KE/N   PE/N   E/N"

  call init_positions()
  call init_velocities()
  call compute_forces()

  ! --- initial configuration ---
  call write_config("pos_init_thermo.dat", pos)
  call write_config("vel_init_thermo.dat", vel)

  ! --- initial energy ---
  call compute_energy(0)

  do step = 1, nsteps

     call velocity_verlet()

     ! --- thermostat ---
     if (mod(step, thermo_interval) == 0) then
        call thermostat_rescale()
     end if

     if (mod(step,100) == 0) then
        call compute_energy(step)
     end if

  end do

  ! --- final configuration ---
  call write_config("pos_final_thermo.dat", pos)
  call write_config("vel_final_thermo.dat", vel)

  close(10)

end program



!==================================================================

subroutine init_positions()
  use sim_para
  implicit none

  integer :: i, nx ,x ,y ,z
  real(8) :: dx 

  nx = int(n_part**(1.0/3.0)) + 1
  dx = lx / nx

  i = 0
  do x = 0, nx-1
     do y = 0, nx-1
        do z = 0, nx-1
           if (i < n_part) then
              pos(3*i+1) = x*dx
              pos(3*i+2) = y*dx
              pos(3*i+3) = z*dx
              i = i + 1
           end if
        end do
     end do
  end do

end subroutine

!==================================================================

subroutine init_velocities()
  use sim_para
  implicit none

  integer :: i
  real(8) :: vx, vy, vz, avgx, avgy, avgz
  real(8) :: vel_const=dsqrt(12*kbT/mass)

  avgx=0; avgy=0; avgz=0

  do i=1,n_part
     call random_number(vx)
     call random_number(vy)
     call random_number(vz)

     vel(3*i-2) = vel_const*(vx - 0.5d0)
     vel(3*i-1) = vel_const*(vy - 0.5d0)
     vel(3*i)   = vel_const*(vz - 0.5d0)

     avgx = avgx + vel(3*i-2)
     avgy = avgy + vel(3*i-1)
     avgz = avgz + vel(3*i)
  end do

  avgx = avgx/n_part
  avgy = avgy/n_part
  avgz = avgz/n_part

  do i=1,n_part
     vel(3*i-2) = vel(3*i-2) - avgx
     vel(3*i-1) = vel(3*i-1) - avgy
     vel(3*i)   = vel(3*i)   - avgz
  end do

end subroutine

!==================================================================

subroutine compute_forces()
  use sim_para
  implicit none

  integer :: i, j
  real(8) :: x, y, z, r
  real(8) :: sigma6, sigma12
  real(8) :: lj_force, fc

  force = 0.0d0

  sigma6  = sigma**6
  sigma12 = sigma**12

  ! Force at cutoff (for shifting)
  fc = 4.0d0*epsilon*(12.0d0*sigma12/(rc**13) - 6.0d0*sigma6/(rc**7))

  do i = 1, n_part-1

     do j = i+1, n_part

        ! Distance with PBC
        x = pos(3*i-2) - pos(3*j-2)
        y = pos(3*i-1) - pos(3*j-1)
        z = pos(3*i)   - pos(3*j)

        x = pbc(x, lx)
        y = pbc(y, ly)
        z = pbc(z, lz)

        r = sqrt(x*x + y*y + z*z)

        if (r < rc) then

           ! Scalar LJ force (reference style)
           lj_force = 4.0d0*epsilon * ( (12.0d0*sigma12/(r**13)) - (6.0d0*sigma6/(r**7)) ) - fc

           ! Convert to vector components
           force(3*i-2) = force(3*i-2) + lj_force * (x/r)
           force(3*i-1) = force(3*i-1) + lj_force * (y/r)
           force(3*i)   = force(3*i)   + lj_force * (z/r)

           force(3*j-2) = force(3*j-2) - lj_force * (x/r)
           force(3*j-1) = force(3*j-1) - lj_force * (y/r)
           force(3*j)   = force(3*j)   - lj_force * (z/r)

        end if

     end do

  end do

end subroutine

!==================================================================

subroutine velocity_verlet()
  use sim_para
  implicit none

  integer :: i

  ! update positions
  do i=1,n_part
     pos(3*i-2) = pos(3*i-2) + vel(3*i-2)*dt + 0.5d0*force(3*i-2)/mass*dt*dt
     pos(3*i-1) = pos(3*i-1) + vel(3*i-1)*dt + 0.5d0*force(3*i-1)/mass*dt*dt
     pos(3*i)   = pos(3*i)   + vel(3*i)*dt   + 0.5d0*force(3*i)/mass*dt*dt

     pos(3*i-2) = modulo(pos(3*i-2), lx)
     pos(3*i-1) = modulo(pos(3*i-1), ly)
     pos(3*i)   = modulo(pos(3*i), lz)
  end do

  ! half velocity update
  do i=1,n_part
     vel(3*i-2) = vel(3*i-2) + 0.5d0*force(3*i-2)/mass*dt
     vel(3*i-1) = vel(3*i-1) + 0.5d0*force(3*i-1)/mass*dt
     vel(3*i)   = vel(3*i)   + 0.5d0*force(3*i)/mass*dt
  end do

  call compute_forces()

  ! second half velocity
  do i=1,n_part
     vel(3*i-2) = vel(3*i-2) + 0.5d0*force(3*i-2)/mass*dt
     vel(3*i-1) = vel(3*i-1) + 0.5d0*force(3*i-1)/mass*dt
     vel(3*i)   = vel(3*i)   + 0.5d0*force(3*i)/mass*dt
  end do

end subroutine

!==================================================================


subroutine compute_energy(step)
  use sim_para
  implicit none

  integer, intent(in) :: step
  integer :: i, j
  real(8) :: ke, pe
  real(8) :: dx, dy, dz, r
  real(8) :: sigma6, sigma12
  real(8) :: v_shift, fc

  ke = 0.0d0
  pe = 0.0d0

  sigma6  = sigma**6
  sigma12 = sigma**12

  ! -----------------------
  ! KINETIC ENERGY
  ! -----------------------
  do i = 1, n_part
     ke = ke + 0.5d0*mass*(vel(3*i-2)**2 + vel(3*i-1)**2 + vel(3*i)**2)
  end do

  ! -----------------------
  ! SHIFT TERMS (same as force)
  ! -----------------------
  fc = 4.0d0*epsilon*(12.0d0*sigma12/(rc**13) - 6.0d0*sigma6/(rc**7))
  v_shift = 4.0d0*epsilon*((sigma12/(rc**12)) - (sigma6/(rc**6)))

  ! -----------------------
  ! POTENTIAL ENERGY
  ! -----------------------
  do i = 1, n_part-1
     do j = i+1, n_part

        dx = pos(3*i-2) - pos(3*j-2)
        dy = pos(3*i-1) - pos(3*j-1)
        dz = pos(3*i)   - pos(3*j)

        dx = pbc(dx, lx)
        dy = pbc(dy, ly)
        dz = pbc(dz, lz)

        r = sqrt(dx*dx + dy*dy + dz*dz)

        if (r < rc) then
           pe = pe + ( 4.0d0*epsilon*((sigma12/(r**12)) - (sigma6/(r**6))) - v_shift +fc*(r-rc))
        end if

     end do
  end do

  ! -----------------------
  ! OUTPUT (per particle)
  ! -----------------------
  write(*,'(A,I6,3F15.6)') "Step:", step, &
       ke/n_part, pe/n_part, (ke+pe)/n_part

  write(10,'(I8,3F15.6)') step, ke/n_part, pe/n_part, (ke+pe)/n_part     

end subroutine

!====================================================
subroutine thermostat_rescale()
  use sim_para
  implicit none

  integer :: i
  real(8) :: ke, ke_target, scale

  ke = 0.0d0

  ! --- compute current KE ---
  do i = 1, n_part
     ke = ke + 0.5d0*mass*(vel(3*i-2)**2 + vel(3*i-1)**2 + vel(3*i)**2)
  end do

  ! --- target KE (k_B T = 1) ---
  ke_target = 1.5d0 * n_part   ! (3/2) N k_B T, with k_B T = 1

  ! --- scaling factor ---
  scale = sqrt(ke_target / ke)

  ! --- rescale velocities ---
  do i = 1, n_part
     vel(3*i-2) = vel(3*i-2) * scale
     vel(3*i-1) = vel(3*i-1) * scale
     vel(3*i)   = vel(3*i)   * scale
  end do

end subroutine

!===========================================================

subroutine write_config(filename, arr)
  use sim_para
  implicit none

  character(len=*), intent(in) :: filename
  real(8), intent(in) :: arr(3*n_part)

  integer :: i, unit

  unit = 20
  open(unit, file=filename, status="unknown")

  do i = 1, n_part
     write(unit,'(3F15.6)') arr(3*i-2), arr(3*i-1), arr(3*i)
  end do

  close(unit)

end subroutine