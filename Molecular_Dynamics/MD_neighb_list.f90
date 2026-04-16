module sim_para
  implicit none

  integer, parameter :: n_part     = 3600
  integer, parameter :: nsteps_eq  = 50000     ! equilibration (NVT)
  integer, parameter :: nsteps_run = 250000    ! production   (NVE)
  integer, parameter :: neigh_update = 40      ! rebuild neighbour list every 40 steps
  integer, parameter :: gr_collect   = 100     ! collect g(r) every 100 steps

  real(8), parameter :: dt      = 0.0025d0
  real(8), parameter :: mass    = 1.0d0

  real(8), parameter :: lx = 20.0d0, ly = 20.0d0, lz = 20.0d0

  real(8), parameter :: sigma   = 1.0d0
  real(8), parameter :: epsilon = 1.0d0
  real(8), parameter :: kbT     = 1.0d0
  real(8), parameter :: rc      = 2.5d0          ! force cutoff
  real(8), parameter :: rs      = rc + 2.0d0      ! neighbour list shell = 4.5 sigma

  ! g(r) goes up to L/2 (NOT rc)
  real(8), parameter :: gr_max_r = lx / 2.0d0    ! = 10.0 sigma
  real(8), parameter :: dr_bin   = 0.1d0
  integer, parameter :: n_bins   = int(gr_max_r / dr_bin)  ! = 100 bins

  real(8) :: pos(3*n_part), vel(3*n_part), force(3*n_part)

  ! ---- Neighbour list ----
  integer, parameter :: max_neigh = 400
  integer :: no_neighb(n_part)
  integer :: neigh_list(n_part, max_neigh)

  ! ---- g(r) accumulator ----
  real(8) :: gr_hist(n_bins)
  integer :: gr_count

  ! ---- max-neighbour tracking (for question a) ----
  integer :: max_neigh_seen          ! global max over ALL particles ALL steps
  integer :: max_neigh_step          ! step at which it was seen
  integer :: max_neigh_particle      ! which particle had it

contains

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

  ! --- open output files ---
  open(10, file="md_eq3.dat",       status="unknown")
  open(11, file="md_run3.dat",      status="unknown")
  open(12, file="gr3.dat",          status="unknown")
!   open(13, file="max_neighbours.dat", status="unknown")

  write(10,'(A)') "# Step   KE/N   PE/N   E/N   T_inst"
  write(11,'(A)') "# Step   KE/N   PE/N   E/N   T_inst"
!   write(13,'(A)') "# Step   MaxNeighbours   ParticleIndex"

  gr_hist          = 0.0d0
  gr_count         = 0
!   max_neigh_seen   = 0
!   max_neigh_step   = 0
!   max_neigh_particle = 0

  call init_positions()
  call init_velocities()
  call build_neigh_list()
  call compute_forces()

  ! ============================================================
  ! EQUILIBRATION  (NVT — thermostat every 250 steps)
  ! ============================================================
  write(*,*) "===== EQUILIBRATION (NVT) ====="

  do step = 1, nsteps_eq

     call velocity_verlet()

     if (mod(step, neigh_update) == 0) then
        call build_neigh_list()
        ! Track max neighbours during steps 20000-30000 (question a)
      !   if (step >= 20000 .and. step <= 30000) then
      !      call track_max_neighbours(step)
      !   end if
     end if

     if (mod(step, 250) == 0) call thermostat()

     if (mod(step, 2000) == 0) call compute_energy(step, 10)

  end do

  ! Write summary for question (a)
!   write(13,'(A)')        "# ---- Summary for Question (a) ----"
!   write(13,'(A,I8)')     "# Max neighbours seen (steps 20k-30k) : ", max_neigh_seen
!   write(13,'(A,I8)')     "# At step                             : ", max_neigh_step
!   write(13,'(A,I8)')     "# Particle index                      : ", max_neigh_particle
!   write(13,'(A,F10.4)')  "# Mean expected (rho * 4/3 pi rs^3)   : ", &
!        (real(n_part,8)/(lx*ly*lz)) * (4.0d0/3.0d0) * 3.141592653589793d0 * rs**3

!   write(*,'(A,I6,A,I6,A,I6)') &
!        " Max neighbours in steps 20k-30k = ", max_neigh_seen, &
!        "  (particle ", max_neigh_particle, " at step ", max_neigh_step, ")"

  ! ============================================================
  ! PRODUCTION  (NVE — collect g(r))
  ! ============================================================
  write(*,*) "===== PRODUCTION (NVE + g(r)) ====="

  do step = 1, nsteps_run

     call velocity_verlet()

     if (mod(step, neigh_update) == 0) call build_neigh_list()

     if (mod(step, gr_collect)   == 0) call collect_gr()

     if (mod(step, 5000)         == 0) call compute_energy(nsteps_eq + step, 11)

  end do

  call write_gr()

  close(10); close(11); close(12); close(13)

  write(*,*) "Done.  Results in gr.dat, max_neighbours.dat, md_eq.dat, md_run.dat"

end program

!==================================================================

subroutine init_positions()
  use sim_para
  implicit none

  integer :: i, nx, x, y, z
  real(8) :: dx

  nx = int(n_part**(1.0d0/3.0d0)) + 1
  dx = lx / nx

  i = 0
  do x = 0, nx-1
     do y = 0, nx-1
        do z = 0, nx-1
           if (i < n_part) then
              pos(3*i+1) = x * dx
              pos(3*i+2) = y * dx
              pos(3*i+3) = z * dx
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

  integer  :: i
  real(8)  :: vx, vy, vz, avgx, avgy, avgz
  real(8)  :: vel_const

  vel_const = dsqrt(12.0d0 * kbT / mass)
  avgx = 0.0d0; avgy = 0.0d0; avgz = 0.0d0

  do i = 1, n_part
     call random_number(vx); call random_number(vy); call random_number(vz)
     vel(3*i-2) = vel_const * (vx - 0.5d0)
     vel(3*i-1) = vel_const * (vy - 0.5d0)
     vel(3*i)   = vel_const * (vz - 0.5d0)
     avgx = avgx + vel(3*i-2)
     avgy = avgy + vel(3*i-1)
     avgz = avgz + vel(3*i)
  end do

  avgx = avgx / n_part; avgy = avgy / n_part; avgz = avgz / n_part

  do i = 1, n_part
     vel(3*i-2) = vel(3*i-2) - avgx
     vel(3*i-1) = vel(3*i-1) - avgy
     vel(3*i)   = vel(3*i)   - avgz
  end do

end subroutine

!==================================================================
! Build Verlet neighbour list  (shell radius rs = 4.5 sigma)
! Follows lecture pseudocode exactly: store BOTH i->j and j->i
!==================================================================

subroutine build_neigh_list()
  use sim_para
  implicit none

  integer :: ij, jk
  real(8) :: dx, dy, dz, r

  no_neighb  = 0
  neigh_list = 0

  do ij = 1, n_part - 1
     do jk = ij + 1, n_part

        dx = pbc(pos(3*ij-2) - pos(3*jk-2), lx)
        dy = pbc(pos(3*ij-1) - pos(3*jk-1), ly)
        dz = pbc(pos(3*ij)   - pos(3*jk),   lz)

        r = dsqrt(dx*dx + dy*dy + dz*dz)

        if (r < rs) then
           no_neighb(ij) = no_neighb(ij) + 1
           if (no_neighb(ij) <= max_neigh) &
              neigh_list(ij, no_neighb(ij)) = jk

           no_neighb(jk) = no_neighb(jk) + 1
           if (no_neighb(jk) <= max_neigh) &
              neigh_list(jk, no_neighb(jk)) = ij
        end if

     end do
  end do

end subroutine

!==================================================================
! Track maximum neighbours across all particles  (question a)
!==================================================================

subroutine track_max_neighbours(step)
  use sim_para
  implicit none

  integer, intent(in) :: step
  integer :: i, local_max, local_part

  local_max  = 0
  local_part = 0

  do i = 1, n_part
     if (no_neighb(i) > local_max) then
        local_max  = no_neighb(i)
        local_part = i
     end if
  end do

  ! Write every rebuild step in the window to the dat file
  write(13,'(I10,I8,I10)') step, local_max, local_part

  ! Update global maximum
  if (local_max > max_neigh_seen) then
     max_neigh_seen     = local_max
     max_neigh_step     = step
     max_neigh_particle = local_part
  end if

end subroutine

!==================================================================
! Forces via neighbour list  (lecture Image 2 style)
!==================================================================

subroutine compute_forces()
  use sim_para
  implicit none

  integer :: i, j, p
  real(8) :: dx, dy, dz, r
  real(8) :: sigma6, sigma12, lj_force, fc

  force   = 0.0d0
  sigma6  = sigma**6
  sigma12 = sigma**12
  fc = 4.0d0*epsilon*(12.0d0*sigma12/rc**13 - 6.0d0*sigma6/rc**7)

  do i = 1, n_part
     do j = 1, no_neighb(i)
        p = neigh_list(i, j)
        if (p <= i) cycle          ! avoid double-counting

        dx = pbc(pos(3*i-2) - pos(3*p-2), lx)
        dy = pbc(pos(3*i-1) - pos(3*p-1), ly)
        dz = pbc(pos(3*i)   - pos(3*p),   lz)

        r = dsqrt(dx*dx + dy*dy + dz*dz)

        if (r < rc) then
           lj_force = 4.0d0*epsilon*( 12.0d0*sigma12/r**13 - 6.0d0*sigma6/r**7 ) - fc

           force(3*i-2) = force(3*i-2) + lj_force*(dx/r)
           force(3*i-1) = force(3*i-1) + lj_force*(dy/r)
           force(3*i)   = force(3*i)   + lj_force*(dz/r)

           force(3*p-2) = force(3*p-2) - lj_force*(dx/r)
           force(3*p-1) = force(3*p-1) - lj_force*(dy/r)
           force(3*p)   = force(3*p)   - lj_force*(dz/r)
        end if

     end do
  end do

end subroutine

!==================================================================

subroutine velocity_verlet()
  use sim_para
  implicit none

  integer :: i

  do i = 1, n_part
     pos(3*i-2) = modulo(pos(3*i-2) + vel(3*i-2)*dt + 0.5d0*force(3*i-2)/mass*dt*dt, lx)
     pos(3*i-1) = modulo(pos(3*i-1) + vel(3*i-1)*dt + 0.5d0*force(3*i-1)/mass*dt*dt, ly)
     pos(3*i)   = modulo(pos(3*i)   + vel(3*i)*dt   + 0.5d0*force(3*i)/mass*dt*dt,   lz)
  end do

  do i = 1, n_part
     vel(3*i-2) = vel(3*i-2) + 0.5d0*force(3*i-2)/mass*dt
     vel(3*i-1) = vel(3*i-1) + 0.5d0*force(3*i-1)/mass*dt
     vel(3*i)   = vel(3*i)   + 0.5d0*force(3*i)/mass*dt
  end do

  call compute_forces()

  do i = 1, n_part
     vel(3*i-2) = vel(3*i-2) + 0.5d0*force(3*i-2)/mass*dt
     vel(3*i-1) = vel(3*i-1) + 0.5d0*force(3*i-1)/mass*dt
     vel(3*i)   = vel(3*i)   + 0.5d0*force(3*i)/mass*dt
  end do

end subroutine

!==================================================================

subroutine thermostat()
  use sim_para
  implicit none

  integer :: i
  real(8) :: ke, T_inst, scale

  ke = 0.0d0
  do i = 1, n_part
     ke = ke + 0.5d0*mass*(vel(3*i-2)**2 + vel(3*i-1)**2 + vel(3*i)**2)
  end do

  T_inst = (2.0d0 * ke) / (3.0d0 * n_part)
  scale  = dsqrt(kbT / T_inst)

  do i = 1, n_part
     vel(3*i-2) = vel(3*i-2) * scale
     vel(3*i-1) = vel(3*i-1) * scale
     vel(3*i)   = vel(3*i)   * scale
  end do

end subroutine

!==================================================================

subroutine compute_energy(step, unit)
  use sim_para
  implicit none

  integer, intent(in) :: step, unit
  integer  :: i, j, p
  real(8)  :: ke, pe, dx, dy, dz, r
  real(8)  :: sigma6, sigma12, v_shift, fc, T_inst

  ke = 0.0d0; pe = 0.0d0
  sigma6  = sigma**6;  sigma12 = sigma**12
  fc      = 4.0d0*epsilon*(12.0d0*sigma12/rc**13 - 6.0d0*sigma6/rc**7)
  v_shift = 4.0d0*epsilon*(sigma12/rc**12 - sigma6/rc**6)

  do i = 1, n_part
     ke = ke + 0.5d0*mass*(vel(3*i-2)**2 + vel(3*i-1)**2 + vel(3*i)**2)
  end do

  do i = 1, n_part
     do j = 1, no_neighb(i)
        p = neigh_list(i, j)
        if (p <= i) cycle

        dx = pbc(pos(3*i-2)-pos(3*p-2), lx)
        dy = pbc(pos(3*i-1)-pos(3*p-1), ly)
        dz = pbc(pos(3*i)  -pos(3*p),   lz)
        r  = dsqrt(dx*dx + dy*dy + dz*dz)

        if (r < rc) then
           pe = pe + 4.0d0*epsilon*(sigma12/r**12 - sigma6/r**6) &
                   - v_shift + fc*(r - rc)
        end if
     end do
  end do

  T_inst = (2.0d0 * ke) / (3.0d0 * n_part)

  write(*,'(A,I8,4F14.6)')  "Step:", step, ke/n_part, pe/n_part, (ke+pe)/n_part, T_inst
  write(unit,'(I10,4F14.6)') step, ke/n_part, pe/n_part, (ke+pe)/n_part, T_inst

end subroutine

!==================================================================
! Collect g(r) up to L/2 = 10 sigma  (separate all-pairs loop,
! independent of force cutoff rc)
!==================================================================

subroutine collect_gr()
  use sim_para
  implicit none

  integer :: i, j, ibin
  real(8) :: dx, dy, dz, r

  do i = 1, n_part - 1
     do j = i + 1, n_part

        dx = pbc(pos(3*i-2) - pos(3*j-2), lx)
        dy = pbc(pos(3*i-1) - pos(3*j-1), ly)
        dz = pbc(pos(3*i)   - pos(3*j),   lz)

        r = dsqrt(dx*dx + dy*dy + dz*dz)

        ! Collect up to L/2 — NOT limited to rc
        if (r < gr_max_r) then
           ibin = int(r / dr_bin) + 1
           if (ibin >= 1 .and. ibin <= n_bins) then
              gr_hist(ibin) = gr_hist(ibin) + 2.0d0
           end if
        end if

     end do
  end do

  gr_count = gr_count + 1

end subroutine

!==================================================================
! Normalise and write g(r)
! Also prints first-peak height to screen and gr.dat header
!==================================================================

subroutine write_gr()
  use sim_para
  implicit none

  integer  :: ibin, peak_bin
  real(8)  :: r_lo, r_hi, r_mid, shell_vol, rho, gr_val
  real(8)  :: gr_peak, gr_norm(n_bins)

  rho      = real(n_part, 8) / (lx * ly * lz)
  gr_peak  = 0.0d0
  peak_bin = 1

  write(12,'(A)') "# r(sigma)    g(r)"

  do ibin = 1, n_bins
     r_lo      = (ibin - 1) * dr_bin
     r_hi      =  ibin      * dr_bin
     r_mid     = 0.5d0 * (r_lo + r_hi)
     shell_vol = (4.0d0/3.0d0) * 3.141592653589793d0 * (r_hi**3 - r_lo**3)

     gr_val = gr_hist(ibin) / (rho * shell_vol * real(n_part,8) * real(gr_count,8))
     gr_norm(ibin) = gr_val

     write(12,'(2F14.6)') r_mid, gr_val

     if (gr_val > gr_peak) then
        gr_peak  = gr_val
        peak_bin = ibin
     end if
  end do

  write(*,*)
  write(*,'(A,F8.4,A,F8.4)') &
       " First peak of g(r):  r = ", (peak_bin - 0.5d0)*dr_bin, &
       " sigma,   g(r) = ", gr_peak
  write(*,'(A,F8.4)') &
       " g(r) at large r (should -> 1.0): ", gr_norm(n_bins)

end subroutine

!==================================================================