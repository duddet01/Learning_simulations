
module sim_para

  implicit none

  !=========================================================
  ! Simulation parameters
  !=========================================================

  integer, parameter :: n_part = 100
  integer, parameter :: nsteps = 100000

  integer, parameter :: neigh_update = 20
  integer, parameter :: data_collect = 100

  ! Reduced units
  real(8), parameter :: dt      = 0.001d0
  real(8), parameter :: mass    = 1.0d0
  real(8), parameter :: sigma   = 1.0d0
  real(8), parameter :: epsilon = 1.0d0
  real(8), parameter :: kbT     = 1.0d0

  ! Harmonic spring
  real(8), parameter :: kappa = 100.0d0 * kbT
  real(8), parameter :: r0    = sigma

  ! WCA cutoff
  real(8), parameter :: rc = 2.0d0**(1.0d0/6.0d0) * sigma

  ! Neighbor-list skin
  real(8), parameter :: skin = 0.3d0
  real(8), parameter :: rs   = rc + skin

  ! Simulation box
  real(8), parameter :: lx = 120.0d0
  real(8), parameter :: ly = 20.0d0
  real(8), parameter :: lz = 20.0d0

  ! Maximum number of neighbours
  integer, parameter :: max_neigh = 100

  ! Particle data
  real(8) :: pos(3*n_part)
  real(8) :: pos_unwrap(3*n_part)

  real(8) :: vel(3*n_part)
  real(8) :: force(3*n_part)

  ! Neighbour list
  integer :: no_neighb(n_part)
  integer :: neigh_list(n_part,max_neigh)

end module sim_para


program md

  use sim_para

  implicit none

  integer :: i
  integer :: step

  real(8) :: time

  real(8) :: ke
  real(8) :: uwca
  real(8) :: uspring
  real(8) :: etotal
  real(8) :: temp

  real(8) :: px
  real(8) :: py
  real(8) :: pz
  real(8) :: pmag

  !=========================================================
  ! Open output files
  !=========================================================

  open(unit=10, file='energy_momentum.dat', status='replace')

  open(unit=20, file='polymer_trajectory.dat', status='replace')

  write(10,*) '# Step Time KE U_WCA U_spring E_total Px Py Pz |P|'

  write(20,*) '# Step Time Bead X_wrapped Y_wrapped Z_wrapped', &
              ' X_unwrapped Y_unwrapped Z_unwrapped'

  !=========================================================
  ! Initialization
  !=========================================================

  call init_positions()

  call init_velocities()

  call build_neigh_list()

  call compute_forces()

  !=========================================================
  ! Initial energy
  !=========================================================

  call compute_energy(ke, uwca, uspring, etotal, temp, &
                      px, py, pz, pmag)

  write(10,'(I8,1X,F12.6,1X,8ES18.8)') &
       0, 0.0d0, ke, uwca, uspring, etotal, &
       px, py, pz, pmag

  call write_trajectory(0,0.0d0)

  !=========================================================
  ! Molecular dynamics loop
  !=========================================================

  do step = 1, nsteps

     time = step * dt

     !------------------------------------------------------
     ! Velocity Verlet integration
     !------------------------------------------------------

     call velocity_verlet()

     !------------------------------------------------------
     ! Rebuild neighbour list
     !
     ! The list contains particles within rc + skin.
     ! With this small timestep, rebuilding every 20
     ! steps keeps the list safely updated.
     !------------------------------------------------------

     if (mod(step,neigh_update) == 0) then
        call build_neigh_list()
     endif

     !------------------------------------------------------
     ! Calculate energy, momentum and temperature
     !------------------------------------------------------

     if (mod(step,data_collect) == 0) then

        call compute_energy(ke, uwca, uspring, etotal, temp, &
                            px, py, pz, pmag)

        write(10,'(I8,1X,F12.6,1X,8ES18.8)') &
             step, time, ke, uwca, uspring, etotal, &
             px, py, pz, pmag

        !---------------------------------------------------
        ! Save polymer coordinates for animation
        !---------------------------------------------------

        call write_trajectory(step,time)

     endif

  enddo

  !=========================================================
  ! Close files
  !=========================================================

  close(10)
  close(20)

  write(*,*) 'Simulation complete.'
  write(*,*) 'Energy/momentum data: energy_momentum.dat'
  write(*,*) 'Trajectory data:      polymer_trajectory.dat'

end program md


!===========================================================
! INITIAL POSITIONS
!===========================================================

subroutine init_positions()

  use sim_para

  implicit none

  integer :: i

  !---------------------------------------------------------
  ! Straight polymer initially placed along x direction
  !
  ! Beads are separated by r0 = sigma.
  !---------------------------------------------------------

  do i = 1,n_part

     pos(3*i-2) = 10.0d0 + (i-1)*r0
     pos(3*i-1) = ly/2.0d0
     pos(3*i)   = lz/2.0d0

     ! Initially wrapped and unwrapped positions are identical
     pos_unwrap(3*i-2) = pos(3*i-2)
     pos_unwrap(3*i-1) = pos(3*i-1)
     pos_unwrap(3*i)   = pos(3*i)

  enddo

end subroutine init_positions


!===========================================================
! INITIAL VELOCITIES
!===========================================================

subroutine init_velocities()

  use sim_para

  implicit none

  integer :: i

  real(8) :: ran1
  real(8) :: vx_sum
  real(8) :: vy_sum
  real(8) :: vz_sum

  real(8) :: scale

  !---------------------------------------------------------
  ! Same type of velocity initialization as the original
  ! MD code: uniform random velocities.
  !---------------------------------------------------------

  scale = sqrt(12.0d0*kbT/mass)

  vx_sum = 0.0d0
  vy_sum = 0.0d0
  vz_sum = 0.0d0

  do i = 1,n_part

     call random_number(ran1)
     vel(3*i-2) = (ran1 - 0.5d0)*scale

     call random_number(ran1)
     vel(3*i-1) = (ran1 - 0.5d0)*scale

     call random_number(ran1)
     vel(3*i) = (ran1 - 0.5d0)*scale

     vx_sum = vx_sum + vel(3*i-2)
     vy_sum = vy_sum + vel(3*i-1)
     vz_sum = vz_sum + vel(3*i)

  enddo

  !---------------------------------------------------------
  ! Remove center-of-mass velocity.
  ! This makes total initial momentum approximately zero.
  !---------------------------------------------------------

  vx_sum = vx_sum/n_part
  vy_sum = vy_sum/n_part
  vz_sum = vz_sum/n_part

  do i = 1,n_part

     vel(3*i-2) = vel(3*i-2) - vx_sum
     vel(3*i-1) = vel(3*i-1) - vy_sum
     vel(3*i)   = vel(3*i)   - vz_sum

  enddo

end subroutine init_velocities


!===========================================================
! APPLY PERIODIC BOUNDARY CONDITIONS
!===========================================================

subroutine apply_pbc()

  use sim_para

  implicit none

  integer :: i

  !---------------------------------------------------------
  ! Keep the simulation coordinates inside the box.
  !---------------------------------------------------------

  do i = 1,n_part

     pos(3*i-2) = modulo(pos(3*i-2),lx)
     pos(3*i-1) = modulo(pos(3*i-1),ly)
     pos(3*i)   = modulo(pos(3*i),lz)

  enddo

end subroutine apply_pbc


!===========================================================
! MINIMUM IMAGE
!===========================================================

subroutine minimum_image(dx,dy,dz)

  use sim_para

  implicit none

  real(8), intent(inout) :: dx
  real(8), intent(inout) :: dy
  real(8), intent(inout) :: dz

  !---------------------------------------------------------
  ! Minimum-image convention.
  !
  ! This makes particles interact with the closest periodic
  ! image of another particle.
  !---------------------------------------------------------

  dx = dx - lx*nint(dx/lx)
  dy = dy - ly*nint(dy/ly)
  dz = dz - lz*nint(dz/lz)

end subroutine minimum_image


!===========================================================
! BUILD NEIGHBOUR LIST
!===========================================================

subroutine build_neigh_list()

  use sim_para

  implicit none

  integer :: i
  integer :: j

  real(8) :: dx
  real(8) :: dy
  real(8) :: dz
  real(8) :: r

  !---------------------------------------------------------
  ! Reset neighbour counts
  !---------------------------------------------------------

  do i = 1,n_part
     no_neighb(i) = 0
  enddo

  !---------------------------------------------------------
  ! Construct neighbour list.
  !
  ! Both directions are stored:
  !
  ! i -> j
  ! j -> i
  !
  ! The list uses rs = rc + skin.
  !---------------------------------------------------------

  do i = 1,n_part-1

     do j = i+1,n_part

        dx = pos(3*i-2) - pos(3*j-2)
        dy = pos(3*i-1) - pos(3*j-1)
        dz = pos(3*i)   - pos(3*j)

        call minimum_image(dx,dy,dz)

        r = sqrt(dx*dx + dy*dy + dz*dz)

        if (r < rs) then

           if (no_neighb(i) < max_neigh) then

              no_neighb(i) = no_neighb(i) + 1
              neigh_list(i,no_neighb(i)) = j

           endif

           if (no_neighb(j) < max_neigh) then

              no_neighb(j) = no_neighb(j) + 1
              neigh_list(j,no_neighb(j)) = i

           endif

        endif

     enddo

  enddo

end subroutine build_neigh_list


!===========================================================
! COMPUTE FORCES
!===========================================================

subroutine compute_forces()

  use sim_para

  implicit none

  integer :: i
  integer :: n
  integer :: p

  real(8) :: dx
  real(8) :: dy
  real(8) :: dz

  real(8) :: r
  real(8) :: rinv
  real(8) :: sr2
  real(8) :: sr6
  real(8) :: sr12

  real(8) :: wca_force

  real(8) :: spring_force

  !=========================================================
  ! Reset forces
  !=========================================================

  force = 0.0d0

  !=========================================================
  ! WCA forces
  !=========================================================

  do i = 1,n_part

     do n = 1,no_neighb(i)

        p = neigh_list(i,n)

        !---------------------------------------------------
        ! Avoid calculating the same WCA pair twice.
        !---------------------------------------------------

        if (p <= i) cycle

        dx = pos(3*i-2) - pos(3*p-2)
        dy = pos(3*i-1) - pos(3*p-1)
        dz = pos(3*i)   - pos(3*p)

        call minimum_image(dx,dy,dz)

        r = sqrt(dx*dx + dy*dy + dz*dz)

        if (r < rc) then

           rinv = 1.0d0/r

           sr2  = (sigma*rinv)**2
           sr6  = sr2**3
           sr12 = sr6**2

           !------------------------------------------------
           ! WCA force magnitude:
           !
           ! F = 24 epsilon/r [2(sigma/r)^12
           !                  -(sigma/r)^6]
           !------------------------------------------------

           wca_force = 24.0d0*epsilon*rinv * &
                       (2.0d0*sr12 - sr6)

           !------------------------------------------------
           ! dx = r_i - r_j
           !
           ! Therefore dx/r points from j -> i.
           ! This gives the repulsive force on i.
           !------------------------------------------------

           force(3*i-2) = force(3*i-2) + &
                          wca_force*dx*rinv

           force(3*i-1) = force(3*i-1) + &
                          wca_force*dy*rinv

           force(3*i)   = force(3*i) + &
                          wca_force*dz*rinv

           ! Newton's third law
           ! Force on j is opposite.

           force(3*p-2) = force(3*p-2) - &
                          wca_force*dx*rinv

           force(3*p-1) = force(3*p-1) - &
                          wca_force*dy*rinv

           force(3*p)   = force(3*p) - &
                          wca_force*dz*rinv

        endif

     enddo

  enddo


  !=========================================================
  ! Harmonic spring forces
  !=========================================================

  do i = 1,n_part-1

     p = i+1

     !------------------------------------------------------
     ! Vector from bead i to bead p.
     !------------------------------------------------------

     dx = pos(3*p-2) - pos(3*i-2)
     dy = pos(3*p-1) - pos(3*i-1)
     dz = pos(3*p)   - pos(3*i)

     !------------------------------------------------------
     ! IMPORTANT:
     ! Apply minimum image to the bond as well.
     !
     ! This means a bond crossing a periodic boundary
     ! remains physically short.
     !------------------------------------------------------

     call minimum_image(dx,dy,dz)

     r = sqrt(dx*dx + dy*dy + dz*dz)

     if (r > 1.0d-12) then

        !---------------------------------------------------
        ! Harmonic spring force:
        !
        ! F_i = kappa (r-r0) r_hat
        !
        ! where r_hat points from i -> p.
        !---------------------------------------------------

        spring_force = kappa*(r-r0)

        force(3*i-2) = force(3*i-2) + &
                       spring_force*dx/r

        force(3*i-1) = force(3*i-1) + &
                       spring_force*dy/r

        force(3*i)   = force(3*i) + &
                       spring_force*dz/r

        ! Opposite force on p

        force(3*p-2) = force(3*p-2) - &
                       spring_force*dx/r

        force(3*p-1) = force(3*p-1) - &
                       spring_force*dy/r

        force(3*p)   = force(3*p) - &
                       spring_force*dz/r

     endif

  enddo

end subroutine compute_forces


!===========================================================
! VELOCITY VERLET
!===========================================================

subroutine velocity_verlet()

  use sim_para

  implicit none

  integer :: i

  !=========================================================
  ! 1. Update positions
  !=========================================================

  do i = 1,n_part

     !------------------------------------------------------
     ! Save continuous/unwrapped position.
     !------------------------------------------------------

     pos_unwrap(3*i-2) = pos_unwrap(3*i-2) + &
                         vel(3*i-2)*dt + &
                         0.5d0*force(3*i-2)/mass*dt*dt

     pos_unwrap(3*i-1) = pos_unwrap(3*i-1) + &
                         vel(3*i-1)*dt + &
                         0.5d0*force(3*i-1)/mass*dt*dt

     pos_unwrap(3*i)   = pos_unwrap(3*i) + &
                         vel(3*i)*dt + &
                         0.5d0*force(3*i)/mass*dt*dt

     !------------------------------------------------------
     ! Wrapped position used by the simulation.
     !------------------------------------------------------

     pos(3*i-2) = modulo(pos_unwrap(3*i-2),lx)
     pos(3*i-1) = modulo(pos_unwrap(3*i-1),ly)
     pos(3*i)   = modulo(pos_unwrap(3*i),lz)

  enddo


  !=========================================================
  ! 2. Half velocity update
  !=========================================================

  do i = 1,n_part

     vel(3*i-2) = vel(3*i-2) + &
                  0.5d0*force(3*i-2)/mass*dt

     vel(3*i-1) = vel(3*i-1) + &
                  0.5d0*force(3*i-1)/mass*dt

     vel(3*i)   = vel(3*i) + &
                  0.5d0*force(3*i)/mass*dt

  enddo


  !=========================================================
  ! 3. Calculate new forces
  !=========================================================

  call compute_forces()


  !=========================================================
  ! 4. Second half velocity update
  !=========================================================

  do i = 1,n_part

     vel(3*i-2) = vel(3*i-2) + &
                  0.5d0*force(3*i-2)/mass*dt

     vel(3*i-1) = vel(3*i-1) + &
                  0.5d0*force(3*i-1)/mass*dt

     vel(3*i)   = vel(3*i) + &
                  0.5d0*force(3*i)/mass*dt

  enddo

end subroutine velocity_verlet


!===========================================================
! COMPUTE ENERGY AND MOMENTUM
!===========================================================

subroutine compute_energy(ke,uwca,uspring,etotal,temp, &
                          px,py,pz,pmag)

  use sim_para

  implicit none

  integer :: i
  integer :: j

  real(8), intent(out) :: ke
  real(8), intent(out) :: uwca
  real(8), intent(out) :: uspring
  real(8), intent(out) :: etotal
  real(8), intent(out) :: temp

  real(8), intent(out) :: px
  real(8), intent(out) :: py
  real(8), intent(out) :: pz
  real(8), intent(out) :: pmag

  real(8) :: dx
  real(8) :: dy
  real(8) :: dz

  real(8) :: r
  real(8) :: sr6
  real(8) :: sr12

  !=========================================================
  ! Initialize
  !=========================================================

  ke      = 0.0d0
  uwca    = 0.0d0
  uspring = 0.0d0

  px = 0.0d0
  py = 0.0d0
  pz = 0.0d0

  !=========================================================
  ! Kinetic energy and momentum
  !=========================================================

  do i = 1,n_part

     ke = ke + 0.5d0*mass * &
          (vel(3*i-2)**2 + &
           vel(3*i-1)**2 + &
           vel(3*i)**2)

     px = px + mass*vel(3*i-2)
     py = py + mass*vel(3*i-1)
     pz = pz + mass*vel(3*i)

  enddo

  pmag = sqrt(px*px + py*py + pz*pz)

  !=========================================================
  ! WCA potential energy
  !=========================================================

  do i = 1,n_part-1

     do j = i+1,n_part

        dx = pos(3*i-2) - pos(3*j-2)
        dy = pos(3*i-1) - pos(3*j-1)
        dz = pos(3*i)   - pos(3*j)

        call minimum_image(dx,dy,dz)

        r = sqrt(dx*dx + dy*dy + dz*dz)

        if (r < rc) then

           sr6  = (sigma/r)**6
           sr12 = sr6**2

           uwca = uwca + &
                  4.0d0*epsilon * &
                  (sr12 - sr6 + 0.25d0)

        endif

     enddo

  enddo

  !=========================================================
  ! Harmonic spring potential
  !=========================================================

  do i = 1,n_part-1

     j = i+1

     dx = pos(3*j-2) - pos(3*i-2)
     dy = pos(3*j-1) - pos(3*i-1)
     dz = pos(3*j)   - pos(3*i)

     call minimum_image(dx,dy,dz)

     r = sqrt(dx*dx + dy*dy + dz*dz)

     uspring = uspring + &
               0.5d0*kappa*(r-r0)**2

  enddo

  !=========================================================
  ! Total energy
  !=========================================================

  etotal = ke + uwca + uspring

  !=========================================================
  ! Instantaneous temperature
  !
  ! 3N degrees of freedom minus 3 COM degrees of freedom
  !=========================================================

  temp = 2.0d0*ke / &
         (3.0d0*n_part - 3.0d0)

end subroutine compute_energy


!===========================================================
! WRITE TRAJECTORY
!===========================================================

subroutine write_trajectory(step,time)

  use sim_para

  implicit none

  integer, intent(in) :: step

  real(8), intent(in) :: time

  integer :: i

  !---------------------------------------------------------
  ! Save both wrapped and unwrapped coordinates.
  !
  ! Wrapped:
  !   shows where the bead is inside the periodic box.
  !
  ! Unwrapped:
  !   gives a continuous polymer trajectory and is preferred
  !   for making the time-evolution animation.
  !---------------------------------------------------------

  do i = 1,n_part

     write(20,'(I8,1X,F12.6,1X,I6,1X,6F14.6)') &
          step, time, i, &
          pos(3*i-2), &
          pos(3*i-1), &
          pos(3*i), &
          pos_unwrap(3*i-2), &
          pos_unwrap(3*i-1), &
          pos_unwrap(3*i)

  enddo

end subroutine write_trajectory
