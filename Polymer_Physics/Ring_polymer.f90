!=================================================================
! RING POLYMER: LANGEVIN DYNAMICS (BAOAB)
!
! Bead-spring ring polymer in free space (no periodic boundaries,
! no walls). Bead i is bonded to i+1, and bead N is bonded back
! to bead 1. All bead pairs interact through the WCA potential.
!
! Outputs, for every N in n_list:
!
!   ring_scaling_results.dat     one line per N: Rg, Rg error,
!                                tau estimate, tau measured,
!                                number of steps, neighbor-list
!                                rebuilds, wall-clock time
!   ring_Rg_trace_N<N>.dat       Rg(t) over the whole run
!   ring_Rouse_AC_N<N>.dat       p=1 ring Rouse-mode autocorrelation
!   ring_MSD_N<N>.dat            g1, g2, g3 MSDs at log-spaced lags
!   ring_track_traj_N<N>.dat     raw positions of the tracked beads
!                                and the center of mass (optional)
!
! PARALLEL VERSION (OpenMP)
!
! Each ring size N is one independent job. Jobs are handed out
! with schedule(dynamic,1): whenever a thread finishes, it takes
! the next unstarted job. The jobs are queued largest N first
! (cost grows roughly as N^3), so the long runs start immediately
! and the short ones fill the gaps at the end.
!
! Every job runs on a single thread, so there is no shared state
! between rings. Each ring is seeded from (seed_value, N) only,
! so the results do not depend on how many threads are used or on
! which thread picked up which job.
!
! Compile:  gfortran -O3 -march=native -fopenmp ring_polymer_omp.f90 -o ring_omp
! Run:      ./ring_omp                     (uses n_threads below)
!           OMP_NUM_THREADS=4 ./ring_omp   (overrides n_threads)
!=================================================================

module sim_para

  implicit none

  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: i8 = selected_int_kind(18)

  real(dp), parameter :: pi = 3.141592653589793d0

  !==============================================================
  ! Ring sizes
  !==============================================================

  integer, parameter :: n_sizes = 16
  integer, parameter :: n_list(n_sizes) = &
       (/100, 110, 120, 130, 140, 150, 160, 170, &
         180, 190, 200, 220, 240, 260 , 280 , 300/)

  ! Number of OpenMP threads (ignored if OMP_NUM_THREADS is set)
  integer, parameter :: n_threads = 8

  !==============================================================
  ! Tracked monomers for the MSD
  !
  ! n_track beads, equally spaced along the ring contour:
  !     i_k = 1 + (k-1)*N/n_track,   k = 1..n_track
  !==============================================================

  integer, parameter :: n_track = 5

  !==============================================================
  ! Run length
  !
  ! Set from the estimated relaxation time of the slowest (p=1)
  ! ring Rouse mode,
  !
  !     tau_ring ~ gam * N^(1+2 nu) / (12 pi^2)
  !
  ! which is 1/4 of the linear-chain value: on a ring the slowest
  ! mode has half the wavelength of the linear p=1 mode.
  !==============================================================

  real(dp), parameter :: nu = 0.588d0

  real(dp), parameter :: n_equil_tau = 10.0d0    ! equilibration, in tau
  real(dp), parameter :: n_prod_tau  = 20.0d0    ! production, in tau
  real(dp), parameter :: sample_tau  = 0.05d0    ! Rg / Rouse sampling, in tau

  ! MSD sampling interval, in absolute time units (not in tau),
  ! so that short-time behaviour is resolved for every N.
  real(dp), parameter :: msd_sample_time = 0.1d0

  ! Points per decade for the log-spaced MSD lag times
  integer, parameter :: n_lag_per_decade = 20

  ! Write the raw tracked-bead trajectories (large files)
  logical, parameter :: write_track_traj = .true.

  real(dp), parameter :: dt   = 0.001d0
  real(dp), parameter :: mass = 1.0d0

  ! Reduced units
  real(dp), parameter :: sigma = 1.0d0
  real(dp), parameter :: eps   = 1.0d0
  real(dp), parameter :: kbT   = 1.0d0

  !==============================================================
  ! Harmonic spring
  !==============================================================

  real(dp), parameter :: kappa = 100.0d0
  real(dp), parameter :: r0    = 1.0d0

  !==============================================================
  ! WCA potential
  !==============================================================

  real(dp), parameter :: rc = 2.0d0**(1.0d0/6.0d0)

  ! Neighbor-list skin
  real(dp), parameter :: skin = 0.3d0
  real(dp), parameter :: rs   = rc + skin

  !==============================================================
  ! Langevin friction
  !==============================================================

  real(dp), parameter :: gam = 1.0d0

  !==============================================================
  ! Random seed (fixed, for reproducibility)
  !==============================================================

  integer, parameter :: seed_value = 100

  !==============================================================
  ! Box-Muller spare deviate: one private copy per thread
  !==============================================================

  real(dp) :: rng_spare = 0.0d0
  logical  :: rng_have_spare = .false.
!$omp threadprivate(rng_spare, rng_have_spare)

end module sim_para


!=================================================================
! MAIN PROGRAM
!=================================================================

program ring_polymer_scaling

  use sim_para
  use omp_lib
  implicit none

  integer :: i, j, itmp
  integer :: ijob, isize, n_part
  integer :: u_res
  integer :: n_used
  integer :: env_len
  integer :: order(n_sizes)

  integer(i8) :: nsteps_done
  integer(i8) :: n_rebuild
  integer(i8) :: res_steps(n_sizes)
  integer(i8) :: res_rebuild(n_sizes)

  real(dp) :: t_start
  real(dp) :: t_job
  real(dp) :: t_total
  real(dp) :: rg_avg
  real(dp) :: rg_err
  real(dp) :: tau_meas
  real(dp) :: tau_est
  real(dp) :: res_rg(n_sizes)
  real(dp) :: res_rg_err(n_sizes)
  real(dp) :: res_tau_est(n_sizes)
  real(dp) :: res_tau_meas(n_sizes)
  real(dp) :: res_time(n_sizes)

  logical :: res_done(n_sizes)

  res_done     = .false.
  res_steps    = 0_i8
  res_rebuild  = 0_i8
  res_rg       = 0.0d0
  res_rg_err   = 0.0d0
  res_tau_est  = 0.0d0
  res_tau_meas = 0.0d0
  res_time     = 0.0d0

  !---------------------------------------------------------------
  ! Job queue: largest N first (insertion sort, descending).
  ! Starting the longest jobs first keeps all threads busy to the
  ! end instead of leaving one thread alone with a big job.
  !---------------------------------------------------------------

  order = [(i, i = 1, n_sizes)]

  do i = 2, n_sizes
     itmp = order(i)
     j = i - 1
     do while (j >= 1)
        if (n_list(order(j)) >= n_list(itmp)) exit
        order(j+1) = order(j)
        j = j - 1
     end do
     order(j+1) = itmp
  end do

  !---------------------------------------------------------------
  ! Threads: n_threads, unless OMP_NUM_THREADS is set. Never more
  ! threads than jobs.
  !---------------------------------------------------------------

  call get_environment_variable('OMP_NUM_THREADS', length=env_len)

  if (env_len == 0) call omp_set_num_threads(min(n_threads, n_sizes))

  n_used = omp_get_max_threads()

  write(*,'(A,I3,A,I3,A)') 'Running ', n_sizes, ' ring sizes on up to ', &
                           n_used, ' threads (dynamic scheduling)'
  write(*,'(A,*(I5))') 'Job queue (N):', (n_list(order(i)), i = 1, n_sizes)
  write(*,*)

  open(newunit=u_res, file='ring_scaling_results.dat', status='replace')

  write(u_res,'(A)') '# N    Rg    Rg_err    tau_est    tau_measured' // &
                     '    n_steps    n_neighlist_rebuilds    wall_time_s'

  t_start = omp_get_wtime()

  !$omp parallel do schedule(dynamic,1) default(shared) &
  !$omp   private(ijob, isize, n_part, t_job, rg_avg, rg_err, &
  !$omp           tau_meas, tau_est, nsteps_done, n_rebuild)
  do ijob = 1, n_sizes

     isize  = order(ijob)
     n_part = n_list(isize)

     !------------------------------------------------------------
     ! A ring needs at least 4 beads for the circle start to be
     ! free of non-bonded overlaps, and at least n_track beads.
     !------------------------------------------------------------

     if (n_part < max(4, n_track)) then
        !$omp critical (io)
        write(*,'(A,I5,A)') 'Skipping N = ', n_part, ' (too small)'
        !$omp end critical (io)
        cycle
     end if

     t_job = omp_get_wtime()

     call run_ring(n_part, rg_avg, rg_err, tau_meas, tau_est, &
                   nsteps_done, n_rebuild)

     t_job = omp_get_wtime() - t_job

     res_rg(isize)       = rg_avg
     res_rg_err(isize)   = rg_err
     res_tau_est(isize)  = tau_est
     res_tau_meas(isize) = tau_meas
     res_steps(isize)    = nsteps_done
     res_rebuild(isize)  = n_rebuild
     res_time(isize)     = t_job
     res_done(isize)     = .true.

     !------------------------------------------------------------
     ! Results are appended as jobs finish (in completion order), so
     ! a partial file exists even if the run is stopped early.
     ! The file is rewritten in order of N at the very end.
     !------------------------------------------------------------

     !$omp critical (io)

     write(u_res,'(I6,4(2X,ES16.8),2(2X,I14),2X,F14.3)') &
          n_part, rg_avg, rg_err, tau_est, tau_meas, &
          nsteps_done, n_rebuild, t_job
     flush(u_res)

     write(*,'(A,I5,A,I2,A)') 'DONE  N = ', n_part, '   (thread ', &
                              omp_get_thread_num(), ')'
     write(*,'(A,F15.8,A,F12.8)') '  Average Rg                 = ', rg_avg, &
                                  ' +/- ', rg_err
     write(*,'(A,F15.4)') '  Measured tau (ring p=1)    = ', tau_meas
     write(*,'(A,I14)')   '  Steps (iterations) run     = ', nsteps_done
     write(*,'(A,I14)')   '  Neighbor-list rebuilds     = ', n_rebuild
     write(*,'(A,F15.3)') '  Wall-clock time (seconds)  = ', t_job
     write(*,*)

     !$omp end critical (io)

  end do
  !$omp end parallel do

  t_total = omp_get_wtime() - t_start

  close(u_res)

  !---------------------------------------------------------------
  ! Rewrite the results file sorted by N
  !---------------------------------------------------------------

  open(newunit=u_res, file='ring_scaling_results.dat', status='replace')

  write(u_res,'(A)') '# N    Rg    Rg_err    tau_est    tau_measured' // &
                     '    n_steps    n_neighlist_rebuilds    wall_time_s'

  do ijob = n_sizes, 1, -1      ! order() is descending in N
     isize = order(ijob)
     if (.not. res_done(isize)) cycle
     write(u_res,'(I6,4(2X,ES16.8),2(2X,I14),2X,F14.3)') &
          n_list(isize), res_rg(isize), res_rg_err(isize), &
          res_tau_est(isize), res_tau_meas(isize), &
          res_steps(isize), res_rebuild(isize), res_time(isize)
  end do

  close(u_res)

  write(*,*) '======================================'
  write(*,*) 'Finished.'
  write(*,'(A,F12.2,A)') ' Total wall-clock time       = ', t_total, ' s'
  write(*,'(A,F12.2,A)') ' Sum of all job times        = ', sum(res_time), ' s'
  write(*,'(A,F12.2)')   ' Speed-up vs one thread      = ', sum(res_time)/t_total
  write(*,*) 'Results written to ring_scaling_results.dat'
  write(*,*) '======================================'

end program ring_polymer_scaling


!=================================================================
! SEED THE RANDOM NUMBER GENERATOR FOR ONE RING
!
! Called at the start of every job. The seed depends only on
! seed_value and N, so a ring gives the same result whichever
! thread runs it and in whatever order. gfortran keeps a separate
! generator state for each OpenMP thread, so seeding here only
! affects the calling thread.
!=================================================================

subroutine seed_run(n_part)

  use sim_para
  implicit none

  integer, intent(in) :: n_part

  integer :: i
  integer :: nseed
  integer, allocatable :: seed(:)

  call random_seed(size=nseed)
  allocate(seed(nseed))

  do i = 1, nseed
     seed(i) = seed_value + 7919*n_part + 104729*i
  end do

  call random_seed(put=seed)

  deallocate(seed)

  ! Throw away any Box-Muller deviate left over from a previous job
  rng_have_spare = .false.

end subroutine seed_run


!=================================================================
! RUN ONE RING
!=================================================================

subroutine run_ring(n_part, rg_avg, rg_err, tau_meas, tau_est, &
                    nsteps, n_rebuild)

  use sim_para
  use omp_lib
  implicit none

  integer, intent(in) :: n_part

  real(dp), intent(out) :: rg_avg
  real(dp), intent(out) :: rg_err
  real(dp), intent(out) :: tau_meas
  real(dp), intent(out) :: tau_est

  integer(i8), intent(out) :: nsteps
  integer(i8), intent(out) :: n_rebuild

  integer(i8) :: step
  integer(i8) :: nequil
  integer(i8) :: sample_every
  integer(i8) :: msd_every

  integer :: k
  integer :: n_samples
  integer :: max_mode_samples
  integer :: n_mode_samples
  integer :: max_msd_samples
  integer :: n_msd_samples
  integer :: max_neigh

  integer :: track_idx(n_track)
  integer :: u_rg
  integer :: u_traj

  integer, allocatable :: no_neighb(:)
  integer, allocatable :: neigh_list(:,:)

  real(dp), allocatable :: pos(:,:)
  real(dp), allocatable :: vel(:,:)
  real(dp), allocatable :: force(:,:)
  real(dp), allocatable :: pos_neigh(:,:)
  real(dp), allocatable :: mode1(:,:)
  real(dp), allocatable :: traj(:,:,:)

  real(dp) :: rg2
  real(dp) :: rg2_sum
  real(dp) :: rg4_sum
  real(dp) :: rg2_avg
  real(dp) :: variance
  real(dp) :: mode_tmp(6)
  real(dp) :: com(3)

  character(len=64) :: fname

  logical :: check_neigh_list

  !---------------------------------------------------------------
  ! Run length from the estimated ring relaxation time
  !---------------------------------------------------------------

  tau_est = gam*real(n_part,dp)**(1.0d0+2.0d0*nu) &
            /(12.0d0*pi*pi)

  nequil = nint(n_equil_tau*tau_est/dt, kind=i8)

  nsteps = nequil + nint(n_prod_tau*tau_est/dt, kind=i8)

  sample_every = max(1_i8, nint(sample_tau*tau_est/dt, kind=i8))

  msd_every = max(1_i8, nint(msd_sample_time/dt, kind=i8))

  !---------------------------------------------------------------
  ! Tracked beads, equally spaced along the ring
  !---------------------------------------------------------------

  do k = 1, n_track
     track_idx(k) = 1 + ((k-1)*n_part)/n_track
  end do

  call seed_run(n_part)

  !$omp critical (io)
  write(*,'(A,I2,A,I5,A,F9.2,A,I12,A,*(I5))') &
       'START (thread ', omp_get_thread_num(), ')  N = ', n_part, &
       '  tau_est = ', tau_est, '  total steps = ', nsteps, &
       '  tracked beads:', track_idx
  !$omp end critical (io)

  !---------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------

  max_neigh = n_part

  allocate(pos(3,n_part))
  allocate(vel(3,n_part))
  allocate(force(3,n_part))

  allocate(no_neighb(n_part))
  allocate(neigh_list(max_neigh,n_part))
  allocate(pos_neigh(3,n_part))

  max_mode_samples = int((nsteps-nequil)/sample_every) + 2
  allocate(mode1(6,max_mode_samples))
  n_mode_samples = 0

  ! traj(:, k, t): k = 1..n_track are the tracked beads,
  !                k = n_track+1 is the center of mass
  max_msd_samples = int((nsteps-nequil)/msd_every) + 2
  allocate(traj(3,n_track+1,max_msd_samples))
  n_msd_samples = 0

  !---------------------------------------------------------------
  ! Initialize ring, neighbor list and forces
  !---------------------------------------------------------------

  call initialize_ring(n_part, pos, vel)

  call build_neigh_list(n_part, pos, no_neighb, neigh_list, &
                        pos_neigh, max_neigh)

  n_rebuild = 0

  call calculate_forces(n_part, pos, force, &
                        no_neighb, neigh_list, max_neigh)

  !---------------------------------------------------------------
  ! Rg averaging: <Rg^2> is averaged, since <Rg^2>^(1/2) ~ N^nu
  !---------------------------------------------------------------

  rg2_sum = 0.0d0
  rg4_sum = 0.0d0
  n_samples = 0

  !---------------------------------------------------------------
  ! Output files written during the run
  !---------------------------------------------------------------

  write(fname,'(A,I0,A)') 'ring_Rg_trace_N', n_part, '.dat'
  open(newunit=u_rg, file=trim(fname), status='replace')
  write(u_rg,'(A)') '# time    Rg'

  if (write_track_traj) then
     write(fname,'(A,I0,A)') 'ring_track_traj_N', n_part, '.dat'
     open(newunit=u_traj, file=trim(fname), status='replace')
     write(u_traj,'(A,*(I6))') '# tracked bead indices:', track_idx
     write(u_traj,'(A)') '# t_prod   x y z of each tracked bead in that order,' // &
                     ' then x y z of the center of mass'
  end if

  !---------------------------------------------------------------
  ! Langevin dynamics, BAOAB
  !
  ! B - half force kick
  ! A - half position update
  ! O - Langevin thermostat
  ! A - half position update
  ! B - half force kick with the new forces
  !---------------------------------------------------------------

  do step = 1, nsteps

     call half_kick(n_part, vel, force)

     call langevin_step(n_part, pos, vel)

     if (check_neigh_list(n_part, pos, pos_neigh)) then

        call build_neigh_list(n_part, pos, no_neighb, neigh_list, &
                              pos_neigh, max_neigh)

        n_rebuild = n_rebuild + 1

     end if

     call calculate_forces(n_part, pos, force, &
                           no_neighb, neigh_list, max_neigh)

     call half_kick(n_part, vel, force)

     !------------------------------------------------------------
     ! Rg and Rouse-mode sampling
     !
     ! The trace is written throughout; averages only after
     ! equilibration.
     !------------------------------------------------------------

     if (mod(step, sample_every) == 0) then

        call calculate_rg(n_part, pos, rg2)

        write(u_rg,'(F15.4,2X,F15.8)') real(step,dp)*dt, sqrt(rg2)

        if (step > nequil) then

           rg2_sum = rg2_sum + rg2
           rg4_sum = rg4_sum + rg2*rg2
           n_samples = n_samples + 1

           call calculate_ring_mode1(n_part, pos, mode_tmp)
           n_mode_samples = n_mode_samples + 1
           mode1(:,n_mode_samples) = mode_tmp

        end if

     end if

     !------------------------------------------------------------
     ! MSD sampling (production only)
     !------------------------------------------------------------

     if (step > nequil) then

        if (mod(step-nequil, msd_every) == 0) then

           n_msd_samples = n_msd_samples + 1

           com = sum(pos, dim=2)/real(n_part,dp)

           do k = 1, n_track
              traj(:,k,n_msd_samples) = pos(:,track_idx(k))
           end do

           traj(:,n_track+1,n_msd_samples) = com

           if (write_track_traj) then
              write(u_traj,'(*(ES16.8,1X))') real(step-nequil,dp)*dt, &
                                         traj(:,:,n_msd_samples)
           end if

        end if

     end if

  end do

  close(u_rg)
  if (write_track_traj) close(u_traj)

  !---------------------------------------------------------------
  ! Final average Rg
  !
  ! The error bar assumes independent samples, so it is a lower
  ! bound. Compare runs with different seeds for the real one.
  !---------------------------------------------------------------

  rg2_avg  = rg2_sum/real(n_samples,dp)
  variance = rg4_sum/real(n_samples,dp) - rg2_avg*rg2_avg

  rg_avg = sqrt(rg2_avg)
  rg_err = 0.5d0*sqrt(max(variance,0.0d0)/real(n_samples,dp))/rg_avg

  !---------------------------------------------------------------
  ! Relaxation time and MSD
  !---------------------------------------------------------------

  call calculate_relaxation_time(n_part, n_mode_samples, &
                                 real(sample_every,dp)*dt, &
                                 mode1(:,1:n_mode_samples), tau_meas)

  call calculate_msd(n_part, n_msd_samples, &
                     real(msd_every,dp)*dt, track_idx, &
                     traj(:,:,1:n_msd_samples))

  !---------------------------------------------------------------
  ! Deallocate
  !---------------------------------------------------------------

  deallocate(pos, vel, force)
  deallocate(no_neighb, neigh_list, pos_neigh)
  deallocate(mode1, traj)

end subroutine run_ring


!=================================================================
! INITIALIZE RING
!
! Beads are placed uniformly on a circle in the xy plane, with
! neighbouring beads exactly r0 apart (chord length):
!
!     2 R sin(pi/N) = r0   =>   R = r0 / (2 sin(pi/N))
!
! The closest non-bonded pair (i, i+2) is then 2 r0 cos(pi/N)
! apart, which is larger than rc for N >= 4, so nothing overlaps.
! The ring is centered at the origin.
!=================================================================

subroutine initialize_ring(n_part, pos, vel)

  use sim_para
  implicit none

  integer, intent(in) :: n_part

  real(dp), intent(out) :: pos(3,n_part)
  real(dp), intent(out) :: vel(3,n_part)

  integer :: i

  real(dp) :: radius
  real(dp) :: theta
  real(dp) :: vcm(3)

  real(dp) :: gaussian_random

  radius = r0/(2.0d0*sin(pi/real(n_part,dp)))

  do i = 1, n_part

     theta = 2.0d0*pi*real(i-1,dp)/real(n_part,dp)

     pos(1,i) = radius*cos(theta)
     pos(2,i) = radius*sin(theta)
     pos(3,i) = 0.0d0

  end do

  !---------------------------------------------------------------
  ! Maxwell-Boltzmann velocities, center-of-mass velocity removed
  !---------------------------------------------------------------

  do i = 1, n_part
     vel(1,i) = gaussian_random()*sqrt(kbT/mass)
     vel(2,i) = gaussian_random()*sqrt(kbT/mass)
     vel(3,i) = gaussian_random()*sqrt(kbT/mass)
  end do

  vcm = sum(vel, dim=2)/real(n_part,dp)

  do i = 1, n_part
     vel(:,i) = vel(:,i) - vcm
  end do

end subroutine initialize_ring


!=================================================================
! HALF FORCE KICK (B step of BAOAB)
!=================================================================

subroutine half_kick(n_part, vel, force)

  use sim_para
  implicit none

  integer, intent(in) :: n_part

  real(dp), intent(inout) :: vel(3,n_part)
  real(dp), intent(in)    :: force(3,n_part)

  vel = vel + 0.5d0*dt*force/mass

end subroutine half_kick


!=================================================================
! LANGEVIN STEP (A O A part of BAOAB)
!
! Free space, so there is no wrapping: pos is always the
! continuous (unwrapped) position, which is what the MSD needs.
!=================================================================

subroutine langevin_step(n_part, pos, vel)

  use sim_para
  implicit none

  integer, intent(in) :: n_part

  real(dp), intent(inout) :: pos(3,n_part)
  real(dp), intent(inout) :: vel(3,n_part)

  integer :: i

  real(dp) :: c1
  real(dp) :: c2

  real(dp) :: gaussian_random

  c1 = exp(-gam*dt/mass)
  c2 = sqrt((kbT/mass)*(1.0d0-c1*c1))

  ! A: first half drift
  pos = pos + 0.5d0*dt*vel

  ! O: Ornstein-Uhlenbeck thermostat
  do i = 1, n_part
     vel(1,i) = c1*vel(1,i) + c2*gaussian_random()
     vel(2,i) = c1*vel(2,i) + c2*gaussian_random()
     vel(3,i) = c1*vel(3,i) + c2*gaussian_random()
  end do

  ! A: second half drift
  pos = pos + 0.5d0*dt*vel

end subroutine langevin_step


!=================================================================
! BUILD NEIGHBOR LIST
!
! Only i -> j (j > i) is stored, since the force loop visits each
! pair once.
!=================================================================

subroutine build_neigh_list(n_part, pos, no_neighb, neigh_list, &
                            pos_neigh, max_neigh)

  use sim_para
  implicit none

  integer, intent(in) :: n_part
  integer, intent(in) :: max_neigh

  real(dp), intent(in) :: pos(3,n_part)

  integer, intent(out) :: no_neighb(n_part)
  integer, intent(out) :: neigh_list(max_neigh,n_part)

  real(dp), intent(out) :: pos_neigh(3,n_part)

  integer :: i
  integer :: j

  real(dp) :: dx, dy, dz
  real(dp) :: r2

  no_neighb = 0

  pos_neigh = pos

  do i = 1, n_part-1
     do j = i+1, n_part

        dx = pos(1,j) - pos(1,i)
        dy = pos(2,j) - pos(2,i)
        dz = pos(3,j) - pos(3,i)

        r2 = dx*dx + dy*dy + dz*dz

        if (r2 < rs*rs) then
           no_neighb(i) = no_neighb(i) + 1
           neigh_list(no_neighb(i),i) = j
        end if

     end do
  end do

end subroutine build_neigh_list


!=================================================================
! CHECK WHETHER NEIGHBOR LIST MUST BE REBUILT
!
! Rebuild when the largest displacement since the last build
! exceeds half the skin.
!=================================================================

logical function check_neigh_list(n_part, pos, pos_neigh)

  use sim_para
  implicit none

  integer, intent(in) :: n_part

  real(dp), intent(in) :: pos(3,n_part)
  real(dp), intent(in) :: pos_neigh(3,n_part)

  integer :: i

  real(dp) :: d(3)
  real(dp) :: max_disp2

  max_disp2 = 0.0d0

  do i = 1, n_part
     d = pos(:,i) - pos_neigh(:,i)
     max_disp2 = max(max_disp2, sum(d*d))
  end do

  check_neigh_list = (max_disp2 > (0.5d0*skin)**2)

end function check_neigh_list


!=================================================================
! FORCE CALCULATION
!
! Total force = WCA excluded volume (all pairs, via neighbor list)
!             + harmonic bonds i -- i+1, including N -- 1
!=================================================================

subroutine calculate_forces(n_part, pos, force, &
                            no_neighb, neigh_list, max_neigh)

  use sim_para
  implicit none

  integer, intent(in) :: n_part
  integer, intent(in) :: max_neigh

  real(dp), intent(in)  :: pos(3,n_part)
  real(dp), intent(out) :: force(3,n_part)

  integer, intent(in) :: no_neighb(n_part)
  integer, intent(in) :: neigh_list(max_neigh,n_part)

  integer :: i
  integer :: j
  integer :: p

  real(dp) :: d(3)
  real(dp) :: r2
  real(dp) :: r
  real(dp) :: rinv
  real(dp) :: sr6
  real(dp) :: fmag
  real(dp) :: fs

  force = 0.0d0

  !===============================================================
  ! WCA forces
  !
  ! F = 24 eps / r [2(sigma/r)^12 - (sigma/r)^6]   for r < rc
  !===============================================================

  do i = 1, n_part
     do p = 1, no_neighb(i)

        j = neigh_list(p,i)

        d  = pos(:,i) - pos(:,j)
        r2 = sum(d*d)

        if (r2 < rc*rc) then

           rinv = 1.0d0/sqrt(r2)
           sr6  = (sigma*rinv)**6

           fmag = 24.0d0*eps*rinv*(2.0d0*sr6*sr6 - sr6)

           force(:,i) = force(:,i) + fmag*d*rinv
           force(:,j) = force(:,j) - fmag*d*rinv

        end if

     end do
  end do

  !===============================================================
  ! Harmonic bonds, closing the ring with N -- 1
  !
  ! F = kappa (r-r0) r_vector/r
  !===============================================================

  do i = 1, n_part

     j = i + 1
     if (j > n_part) j = 1

     d = pos(:,j) - pos(:,i)
     r = sqrt(sum(d*d))

     fs = kappa*(r-r0)/r

     force(:,i) = force(:,i) + fs*d
     force(:,j) = force(:,j) - fs*d

  end do

end subroutine calculate_forces


!=================================================================
! SQUARED RADIUS OF GYRATION
!
! Rg^2 = (1/N) sum_i |r_i - r_CM|^2
!=================================================================

subroutine calculate_rg(n_part, pos, rg2)

  use sim_para
  implicit none

  integer, intent(in) :: n_part
  real(dp), intent(in) :: pos(3,n_part)
  real(dp), intent(out) :: rg2

  integer :: i
  real(dp) :: com(3)
  real(dp) :: d(3)

  com = sum(pos, dim=2)/real(n_part,dp)

  rg2 = 0.0d0
  do i = 1, n_part
     d = pos(:,i) - com
     rg2 = rg2 + sum(d*d)
  end do

  rg2 = rg2/real(n_part,dp)

end subroutine calculate_rg


!=================================================================
! FIRST RING ROUSE MODE
!
! A ring is periodic in the bead index, so its normal modes are
! Fourier modes rather than the cosine modes of a linear chain.
! The p=1 mode is a degenerate pair:
!
!   X_1c = (1/N) sum_i cos(2 pi (i-1)/N) r_i
!   X_1s = (1/N) sum_i sin(2 pi (i-1)/N) r_i
!
! mode(1:3) = X_1c, mode(4:6) = X_1s. Both relax with the same
! time, so all six components are used in the autocorrelation.
!=================================================================

subroutine calculate_ring_mode1(n_part, pos, mode)

  use sim_para
  implicit none

  integer, intent(in) :: n_part
  real(dp), intent(in) :: pos(3,n_part)
  real(dp), intent(out) :: mode(6)

  integer :: i
  real(dp) :: q

  mode = 0.0d0

  do i = 1, n_part
     q = 2.0d0*pi*real(i-1,dp)/real(n_part,dp)
     mode(1:3) = mode(1:3) + cos(q)*pos(:,i)
     mode(4:6) = mode(4:6) + sin(q)*pos(:,i)
  end do

  mode = mode/real(n_part,dp)

end subroutine calculate_ring_mode1


!=================================================================
! MEASURE RELAXATION TIME
!
! C(t) = <dX_1(t0+t) . dX_1(t0)> / <|dX_1|^2>
!
! summed over the cos and sin parts of the p=1 ring mode.
! For C(t) ~ exp(-t/tau), the integrated autocorrelation time
!
!     tau = dt_sample [1/2 + sum_{lag>0} C(lag)]
!
! approaches tau. The sum stops at the first C <= 0.
! C(t) is written to ring_Rouse_AC_N*.dat for fitting in Python.
!=================================================================

subroutine calculate_relaxation_time(n_part, n_samples, dt_sample, &
                                     mode, tau)

  use sim_para
  implicit none

  integer, intent(in) :: n_part
  integer, intent(in) :: n_samples
  real(dp), intent(in) :: dt_sample
  real(dp), intent(in) :: mode(6,n_samples)
  real(dp), intent(out) :: tau

  integer :: i
  integer :: lag
  integer :: count

  real(dp) :: mean_mode(6)
  real(dp) :: dx(6)
  real(dp) :: var0
  real(dp) :: c
  real(dp) :: sum_c

  logical :: summing

  integer :: u_ac

  character(len=64) :: fname

  tau = 0.0d0

  if (n_samples < 10) return

  mean_mode = sum(mode, dim=2)/real(n_samples,dp)

  var0 = 0.0d0
  do i = 1, n_samples
     dx = mode(:,i) - mean_mode
     var0 = var0 + sum(dx*dx)
  end do
  var0 = var0/real(n_samples,dp)

  if (var0 <= 1.0d-30) return

  write(fname,'(A,I0,A)') 'ring_Rouse_AC_N', n_part, '.dat'
  open(newunit=u_ac, file=trim(fname), status='replace')
  write(u_ac,'(A)') '# lag_time    normalized_C_ring_Rouse_p1'

  sum_c   = 0.0d0
  summing = .true.

  ! The full C(t) is written up to half the series, so the
  ! decay can be fitted in Python past the first zero crossing.
  do lag = 0, n_samples/2

     c = 0.0d0
     count = n_samples - lag

     do i = 1, count
        c = c + sum((mode(:,i)-mean_mode)*(mode(:,i+lag)-mean_mode))
     end do

     c = c/real(count,dp)/var0

     write(u_ac,'(F15.6,2X,F15.8)') real(lag,dp)*dt_sample, c

     if (lag > 0 .and. summing) then
        if (c <= 0.0d0) then
           summing = .false.
        else
           sum_c = sum_c + c
        end if
     end if

  end do

  close(u_ac)

  tau = dt_sample*(0.5d0 + sum_c)

end subroutine calculate_relaxation_time


!=================================================================
! MEAN SQUARED DISPLACEMENT
!
! For the tracked beads, at log-spaced lag times, averaged over
! all available time origins:
!
!   g1(t) = < |r_i(t0+t) - r_i(t0)|^2 >                bead, lab frame
!   g2(t) = < |(r_i - r_cm)(t0+t) - (r_i - r_cm)(t0)|^2 >   bead, CM frame
!   g3(t) = < |r_cm(t0+t) - r_cm(t0)|^2 >              center of mass
!
! g1 and g2 are averaged over the tracked beads; g1 is also
! written separately for each tracked bead.
!
! At large lags there are few independent origins, so the last
! part of the curve is noisy.
!=================================================================

subroutine calculate_msd(n_part, n_samples, dt_sample, track_idx, traj)

  use sim_para
  implicit none

  integer, intent(in) :: n_part
  integer, intent(in) :: n_samples
  real(dp), intent(in) :: dt_sample
  integer, intent(in) :: track_idx(n_track)
  real(dp), intent(in) :: traj(3,n_track+1,n_samples)

  integer :: j
  integer :: k
  integer :: lag
  integer :: last_lag
  integer :: t0
  integer :: t1
  integer :: n_orig

  real(dp) :: g1(n_track)
  real(dp) :: g2
  real(dp) :: g3
  real(dp) :: d(3)
  real(dp) :: dcm(3)

  integer :: u_msd

  character(len=64) :: fname

  write(fname,'(A,I0,A)') 'ring_MSD_N', n_part, '.dat'
  open(newunit=u_msd, file=trim(fname), status='replace')

  write(u_msd,'(A,*(I6))') '# tracked bead indices:', track_idx
  write(u_msd,'(A)') '# lag_time    g1_avg    g2_avg    g3_cm' // &
                  '    g1 of each tracked bead in the order above'

  last_lag = 0
  j = 0

  do

     lag = nint(10.0d0**(real(j,dp)/real(n_lag_per_decade,dp)))
     j = j + 1

     if (lag > n_samples-1) exit
     if (lag <= last_lag) cycle

     last_lag = lag

     g1 = 0.0d0
     g2 = 0.0d0
     g3 = 0.0d0

     n_orig = n_samples - lag

     do t0 = 1, n_orig

        t1 = t0 + lag

        dcm = traj(:,n_track+1,t1) - traj(:,n_track+1,t0)
        g3  = g3 + sum(dcm*dcm)

        do k = 1, n_track
           d     = traj(:,k,t1) - traj(:,k,t0)
           g1(k) = g1(k) + sum(d*d)
           d     = d - dcm
           g2    = g2 + sum(d*d)
        end do

     end do

     g1 = g1/real(n_orig,dp)
     g2 = g2/real(n_orig*n_track,dp)
     g3 = g3/real(n_orig,dp)

     write(u_msd,'(*(ES16.8,1X))') real(lag,dp)*dt_sample, &
                                sum(g1)/real(n_track,dp), g2, g3, g1

  end do

  close(u_msd)

end subroutine calculate_msd


!=================================================================
! GAUSSIAN RANDOM NUMBER (Box-Muller)
!
! Mean 0, variance 1. Box-Muller produces two deviates at a time;
! the second is kept in rng_spare and returned on the next call.
! rng_spare / rng_have_spare are threadprivate (see sim_para), so
! every thread has its own.
!=================================================================

function gaussian_random() result(g)

  use sim_para
  implicit none

  real(dp) :: g

  real(dp) :: u1
  real(dp) :: u2
  real(dp) :: fac

  if (rng_have_spare) then
     g = rng_spare
     rng_have_spare = .false.
     return
  end if

  call random_number(u1)
  call random_number(u2)

  if (u1 < 1.0d-12) u1 = 1.0d-12

  fac = sqrt(-2.0d0*log(u1))

  g          = fac*cos(2.0d0*pi*u2)
  rng_spare  = fac*sin(2.0d0*pi*u2)

  rng_have_spare = .true.

end function gaussian_random