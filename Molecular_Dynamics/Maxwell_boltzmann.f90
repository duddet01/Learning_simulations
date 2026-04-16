module sim_para
  implicit none

  integer, parameter :: n_part       = 2400
  integer, parameter :: nsteps_eq    = 50000
  integer, parameter :: nsteps_run   = 250000
  integer, parameter :: neigh_update = 40
  integer, parameter :: mb_collect   = 100

  real(8), parameter :: dt      = 0.0025d0
  real(8), parameter :: mass    = 1.0d0

  real(8), parameter :: lx = 20.0d0, ly = 20.0d0, lz = 20.0d0

  real(8), parameter :: sigma   = 1.0d0
  real(8), parameter :: epsilon = 1.0d0
  real(8), parameter :: kbT     = 2.0d0
  real(8), parameter :: rc      = 2.5d0
  real(8), parameter :: rs      = rc + 2.0d0

  real(8) :: pos(3*n_part), vel(3*n_part), force(3*n_part)

  ! Neighbour list
  integer, parameter :: max_neigh = 400
  integer :: no_neighb(n_part)
  integer :: neigh_list(n_part, max_neigh)

  ! Maxwell-Boltzmann speed histogram
  ! Peak of MB is at v = sqrt(2kbT/m) = 1.414 in reduced units
  ! v_max = 10 safely covers the entire tail
  real(8), parameter :: v_max   = 10.0d0
  real(8), parameter :: dv_bin  = 0.05d0
  integer, parameter :: nv_bins = int(v_max / dv_bin)   ! = 200 bins
  real(8) :: mb_hist(nv_bins)
  integer :: mb_count

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

  open(10, file="md_eq2.dat",   status="unknown")
  open(11, file="md_run2.dat",  status="unknown")
  open(14, file="mb_dist2.dat", status="unknown")

  write(10,'(A)') "# Step   KE/N   PE/N   E/N   T_inst"
  write(11,'(A)') "# Step   KE/N   PE/N   E/N   T_inst"
  write(14,'(A)') "# speed(v)   P_sim(v)   P_MB_theory(v)   diff"

  mb_hist  = 0.0d0
  mb_count = 0

  call init_positions()
  call init_velocities()
  call build_neigh_list()
  call compute_forces()

  ! ============================================================
  ! EQUILIBRATION  (NVT)
  ! ============================================================
  write(*,*) "===== EQUILIBRATION (NVT) ====="

  do step = 1, nsteps_eq
     call velocity_verlet()
     if (mod(step, neigh_update) == 0) call build_neigh_list()
     if (mod(step, 250)          == 0) call thermostat()
     if (mod(step, 2000)         == 0) call compute_energy(step, 10)
  end do

  ! ============================================================
  ! PRODUCTION  (NVE — collect MB speed distribution)
  ! ============================================================
  write(*,*) "===== PRODUCTION (NVE) ====="

  do step = 1, nsteps_run
     call velocity_verlet()
     if (mod(step, neigh_update) == 0) call build_neigh_list()
     if (mod(step, mb_collect)   == 0) call collect_mb()
     if (mod(step, 5000)         == 0) call compute_energy(nsteps_eq + step, 11)
  end do

  call write_mb()

  close(10); close(11); close(14)

  write(*,*) "Done. mb_dist2.dat   md_eq2.dat   md_run2.dat"

end program

!==================================================================

subroutine init_positions()
  use sim_para
  implicit none

  integer :: i, nx, x, y, z
  real(8) :: dx

  nx = int(n_part**(1.0d0/3.0d0)) + 1
  dx = lx / nx
  i  = 0

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

  integer :: i
  real(8) :: vx, vy, vz, avgx, avgy, avgz, vel_const

  vel_const = dsqrt(12.0d0*kbT/mass)
  avgx = 0.0d0; avgy = 0.0d0; avgz = 0.0d0

  do i = 1, n_part
     call random_number(vx); call random_number(vy); call random_number(vz)
     vel(3*i-2) = vel_const*(vx - 0.5d0)
     vel(3*i-1) = vel_const*(vy - 0.5d0)
     vel(3*i)   = vel_const*(vz - 0.5d0)
     avgx = avgx + vel(3*i-2)
     avgy = avgy + vel(3*i-1)
     avgz = avgz + vel(3*i)
  end do

  avgx = avgx/n_part; avgy = avgy/n_part; avgz = avgz/n_part

  do i = 1, n_part
     vel(3*i-2) = vel(3*i-2) - avgx
     vel(3*i-1) = vel(3*i-1) - avgy
     vel(3*i)   = vel(3*i)   - avgz
  end do

end subroutine

!==================================================================

subroutine build_neigh_list()
  use sim_para
  implicit none

  integer :: ij, jk
  real(8) :: dx, dy, dz, r

  no_neighb  = 0
  neigh_list = 0

  do ij = 1, n_part-1
     do jk = ij+1, n_part
        dx = pbc(pos(3*ij-2)-pos(3*jk-2), lx)
        dy = pbc(pos(3*ij-1)-pos(3*jk-1), ly)
        dz = pbc(pos(3*ij)  -pos(3*jk),   lz)
        r  = dsqrt(dx*dx + dy*dy + dz*dz)
        if (r < rs) then
           no_neighb(ij) = no_neighb(ij) + 1
           if (no_neighb(ij) <= max_neigh) neigh_list(ij, no_neighb(ij)) = jk
           no_neighb(jk) = no_neighb(jk) + 1
           if (no_neighb(jk) <= max_neigh) neigh_list(jk, no_neighb(jk)) = ij
        end if
     end do
  end do

end subroutine

!==================================================================

subroutine compute_forces()
  use sim_para
  implicit none

  integer :: i, j, p
  real(8) :: dx, dy, dz, r, sigma6, sigma12, lj_force, fc

  force   = 0.0d0
  sigma6  = sigma**6;  sigma12 = sigma**12
  fc = 4.0d0*epsilon*(12.0d0*sigma12/rc**13 - 6.0d0*sigma6/rc**7)

  do i = 1, n_part
     do j = 1, no_neighb(i)
        p = neigh_list(i,j)
        if (p <= i) cycle
        dx = pbc(pos(3*i-2)-pos(3*p-2), lx)
        dy = pbc(pos(3*i-1)-pos(3*p-1), ly)
        dz = pbc(pos(3*i)  -pos(3*p),   lz)
        r  = dsqrt(dx*dx + dy*dy + dz*dz)
        if (r < rc) then
           lj_force = 4.0d0*epsilon*(12.0d0*sigma12/r**13 - 6.0d0*sigma6/r**7) - fc
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
     pos(3*i-2) = modulo(pos(3*i-2)+vel(3*i-2)*dt+0.5d0*force(3*i-2)/mass*dt*dt, lx)
     pos(3*i-1) = modulo(pos(3*i-1)+vel(3*i-1)*dt+0.5d0*force(3*i-1)/mass*dt*dt, ly)
     pos(3*i)   = modulo(pos(3*i)  +vel(3*i)*dt  +0.5d0*force(3*i)/mass*dt*dt,   lz)
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

  T_inst = (2.0d0*ke) / (3.0d0*n_part)
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
  real(8)  :: ke, pe, dx, dy, dz, r, sigma6, sigma12, v_shift, fc, T_inst

  ke = 0.0d0; pe = 0.0d0
  sigma6  = sigma**6;  sigma12 = sigma**12
  fc      = 4.0d0*epsilon*(12.0d0*sigma12/rc**13 - 6.0d0*sigma6/rc**7)
  v_shift = 4.0d0*epsilon*(sigma12/rc**12 - sigma6/rc**6)

  do i = 1, n_part
     ke = ke + 0.5d0*mass*(vel(3*i-2)**2 + vel(3*i-1)**2 + vel(3*i)**2)
  end do

  do i = 1, n_part
     do j = 1, no_neighb(i)
        p = neigh_list(i,j)
        if (p <= i) cycle
        dx = pbc(pos(3*i-2)-pos(3*p-2), lx)
        dy = pbc(pos(3*i-1)-pos(3*p-1), ly)
        dz = pbc(pos(3*i)  -pos(3*p),   lz)
        r  = dsqrt(dx*dx + dy*dy + dz*dz)
        if (r < rc) pe = pe + 4.0d0*epsilon*(sigma12/r**12 - sigma6/r**6) &
                              - v_shift + fc*(r - rc)
     end do
  end do

  T_inst = (2.0d0*ke) / (3.0d0*n_part)
  write(*,'(A,I8,4F14.6)')   "Step:", step, ke/n_part, pe/n_part, (ke+pe)/n_part, T_inst
  write(unit,'(I10,4F14.6)') step, ke/n_part, pe/n_part, (ke+pe)/n_part, T_inst

end subroutine

!==================================================================
! Collect speed |v| = sqrt(vx^2+vy^2+vz^2) into histogram
!==================================================================

subroutine collect_mb()
  use sim_para
  implicit none

  integer :: i, ibin
  real(8) :: spd

  do i = 1, n_part
     spd  = dsqrt(vel(3*i-2)**2 + vel(3*i-1)**2 + vel(3*i)**2)
     ibin = int(spd / dv_bin) + 1
     if (ibin >= 1 .and. ibin <= nv_bins) &
        mb_hist(ibin) = mb_hist(ibin) + 1.0d0
  end do

  mb_count = mb_count + 1

end subroutine

!==================================================================
! Normalise histogram and write alongside exact MB theory
!
!  3D Maxwell-Boltzmann speed PDF:
!
!    P(v) = 4*pi * ( m/(2*pi*kB*T) )^(3/2) * v^2 * exp( -m*v^2/(2*kB*T) )
!
!  Reduced units (m=1, kB=1, T=kbT=1):
!
!    P(v) = 4*pi * (1/(2*pi))^(3/2) * v^2 * exp(-v^2/2)
!
!  Simulation PDF (probability density, normalised to 1):
!
!    P_sim(bin) = count(bin) / ( N_part * n_snapshots * dv_bin )
!==================================================================

subroutine write_mb()
  use sim_para
  implicit none

  integer  :: ibin
  real(8)  :: v_mid, P_sim, P_theory
  real(8)  :: pi, prefactor, norm_factor
  real(8)  :: chi2, integral_sim, integral_theory, max_diff, max_diff_v

  pi          = 3.141592653589793d0
  prefactor   = 4.0d0*pi*(mass/(2.0d0*pi*kbT))**1.5d0
  norm_factor = real(n_part,8) * real(mb_count,8) * dv_bin

  chi2 = 0.0d0; integral_sim = 0.0d0; integral_theory = 0.0d0
  max_diff = 0.0d0; max_diff_v = 0.0d0

  do ibin = 1, nv_bins
     v_mid    = (ibin - 0.5d0) * dv_bin
     P_sim    = mb_hist(ibin) / norm_factor
     P_theory = prefactor * v_mid**2 * dexp(-mass*v_mid**2 / (2.0d0*kbT))

     write(14,'(4F16.8)') v_mid, P_sim, P_theory, P_sim - P_theory

     if (P_theory > 1.0d-10) &
        chi2 = chi2 + (P_sim - P_theory)**2 / P_theory

     if (dabs(P_sim - P_theory) > max_diff) then
        max_diff   = dabs(P_sim - P_theory)
        max_diff_v = v_mid
     end if

     integral_sim    = integral_sim    + P_sim    * dv_bin
     integral_theory = integral_theory + P_theory * dv_bin
  end do

  write(*,*)
  write(*,*) "======= Maxwell-Boltzmann Check (Q8) ======="
  write(*,'(A,F10.6)') " Integral P_sim    (should = 1): ", integral_sim
  write(*,'(A,F10.6)') " Integral P_theory (should = 1): ", integral_theory
  write(*,'(A,F12.6)') " Chi-squared (-> 0 = perfect)  : ", chi2
  write(*,'(A,F10.6,A,F6.3)') &
       " Max |P_sim-P_theory| = ", max_diff, "  at v = ", max_diff_v
  write(*,*) " Plot: gnuplot> plot 'mb_dist.dat' u 1:2 w l t 'sim', '' u 1:3 w l t 'theory'"
  write(*,*) "============================================="

end subroutine

!==================================================================