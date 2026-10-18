module sim_para

  implicit none

  integer, parameter :: dp = kind(1.0d0)

  real(dp), parameter :: pi = 3.141592653589793d0

  !==============================================================
  ! Polymer sizes
  !==============================================================

  integer, parameter :: n_sizes = 10
  integer, parameter :: n_list(n_sizes) = (/100 ,120, 140, 160, 180, 200, 220, 240, 260, 300 /)

  !==============================================================
  ! Independent runs per chain length
  !
  ! Each run uses a different random seed and is a completely
  ! independent measurement. The final Rg and tau are the mean
  ! over runs, and the error bar is the standard error of those
  ! means.
  !
  ! This is the only reliable error estimate here: within a
  ! single run the Rg samples are strongly correlated, so the
  ! usual sqrt(variance/n) formula underestimates the error by
  ! a factor of about 4.
  !==============================================================

  integer, parameter :: n_runs = 8

  !==============================================================
  ! Simulation parameters
  !
  ! The run length is set from the estimated chain relaxation
  ! time
  !
  !     tau_relax ~ gam * N^(1+2 nu) / (3 pi^2)
  !
  ! so that every chain length is equilibrated and sampled for
  ! the same number of relaxation times.
  !
  ! Equilibration was reduced from 10 tau to 3 tau: the Rg
  ! traces show the chain relaxed within about one tau when
  ! started from a self-avoiding walk, so the extra time was
  ! wasted.
  !==============================================================

  real(dp), parameter :: nu = 0.588d0

  real(dp), parameter :: n_equil_tau  = 3.0d0     ! equilibration
  real(dp), parameter :: n_prod_tau   = 20.0d0    ! production
  real(dp), parameter :: sample_tau   = 0.05d0    ! sampling interval

  real(dp), parameter :: dt = 0.001d0
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
  real(dp), parameter :: rs = rc + skin

  !==============================================================
  ! Langevin friction
  !==============================================================

  real(dp), parameter :: gam = 1.0d0

  !==============================================================
  ! State of the Box-Muller generator
  !
  ! These live here so that init_random_seed can throw away a
  ! leftover deviate. Without this the first number of a new
  ! run would come from the previous run's seed.
  !==============================================================

  real(dp) :: gr_spare = 0.0d0
  logical  :: gr_have_spare = .false.

  !$omp threadprivate(gr_spare, gr_have_spare)

end module sim_para


!=================================================================
! MAIN PROGRAM
!=================================================================

program polymer_scaling

  use sim_para
  use iso_fortran_env, only: int64
  !$ use omp_lib
  implicit none

  integer :: isize
  integer :: irun
  integer :: iseed
  integer :: n_part
  integer(int64) :: c0, c1, clock_rate
  integer(int64) :: wall0, wall1

  real(dp) :: rg_avg
  real(dp) :: rg_err
  real(dp) :: tau_meas
  real(dp) :: tau_est

  !---------------------------------------------------------------
  ! Results for every (chain length, run) pair.
  !
  ! Nothing is written to disk from inside the parallel region
  ! except each run's own trace file, which has a unique name.
  ! The summary files are written afterwards, in order.
  !---------------------------------------------------------------

  real(dp) :: rg_all(n_runs,n_sizes)
  real(dp) :: err_all(n_runs,n_sizes)
  real(dp) :: tau_all(n_runs,n_sizes)
  real(dp) :: time_all(n_runs,n_sizes)
  integer  :: seed_all(n_runs,n_sizes)

  real(dp) :: rg_mean
  real(dp) :: rg_sem
  real(dp) :: tau_mean
  real(dp) :: tau_sem
  real(dp) :: cpu_time_n
  real(dp) :: wall_total

  integer :: nthreads

  nthreads = 1
  !$ nthreads = omp_get_max_threads()

  write(*,'(A,I4)') 'OpenMP threads = ', nthreads
  write(*,'(A,I4,A,I4,A,I5,A)') 'Jobs to run    = ', n_sizes, ' sizes x ', &
        n_runs, ' runs = ', n_sizes*n_runs, ' independent jobs'
  write(*,*)

  call system_clock(count_rate=clock_rate)
  call system_clock(count=wall0)

  !---------------------------------------------------------------
  ! Every (chain length, run) pair is an independent job, so the
  ! two loops are collapsed into one pool of n_sizes*n_runs jobs.
  !
  ! The sizes are visited largest first with a dynamic schedule.
  ! A single N=200 run takes far longer than an N=10 run, so if
  ! the big jobs were handed out last most cores would sit idle
  ! waiting for them at the end.
  !---------------------------------------------------------------

  !$omp parallel do collapse(2) schedule(dynamic,1) default(shared) &
  !$omp   private(isize,irun,n_part,iseed,rg_avg,rg_err, &
  !$omp           tau_meas,tau_est,c0,c1)
  do isize = n_sizes, 1, -1

     do irun = 1, n_runs

        call system_clock(count=c0)

        n_part = n_list(isize)

        !---------------------------------------------------------
        ! A distinct seed for every (N, run) pair.
        !
        ! Each thread has its own generator state, so this gives
        ! exactly the same numbers as the serial version no
        ! matter how the jobs are distributed.
        !---------------------------------------------------------

        iseed = 7919*isize + 104729*irun + 12345

        call init_random_seed(iseed)

        call run_polymer(n_part, irun, rg_avg, rg_err, &
                         tau_meas, tau_est)

        call system_clock(count=c1)

        rg_all(irun,isize)   = rg_avg
        err_all(irun,isize)  = rg_err
        tau_all(irun,isize)  = tau_meas
        seed_all(irun,isize) = iseed
        time_all(irun,isize) = real(c1-c0,dp)/real(clock_rate,dp)

        !$omp critical
        write(*,'(A,I5,A,I3,A,F12.6,A,F14.4,A,F10.1,A)') &
              '  N = ', n_part, '  run ', irun, &
              ' :  Rg = ', rg_avg, '   tau = ', tau_meas, &
              '   (', time_all(irun,isize), ' s)'
        !$omp end critical

     end do

  end do
  !$omp end parallel do

  call system_clock(count=wall1)

  wall_total = real(wall1-wall0,dp)/real(clock_rate,dp)

  !---------------------------------------------------------------
  ! Output files
  !
  ! scaling_results3.dat : one line per N, the averaged result
  ! all_runs3.dat        : every individual run, for checking
  !---------------------------------------------------------------

  open(unit=10, file='scaling_results3.dat', status='replace')

  write(10,'(A)') '# N    Rg    Rg_sem    tau_est    tau_meas    tau_sem    n_runs    cpu_time_s'

  open(unit=11, file='all_runs3.dat', status='replace')

  write(11,'(A)') '# N    run    seed    Rg    Rg_err_naive    tau_meas    time_s'

  do isize = 1, n_sizes

     n_part = n_list(isize)

     tau_est = gam*real(n_part,dp)**(1.0d0+2.0d0*nu) &
               /(3.0d0*pi*pi)

     call mean_and_sem(n_runs, rg_all(:,isize),  rg_mean,  rg_sem)
     call mean_and_sem(n_runs, tau_all(:,isize), tau_mean, tau_sem)

     cpu_time_n = sum(time_all(:,isize))

     write(10,'(I6,5(2X,F15.8),2X,I6,2X,F15.8)') &
           n_part, rg_mean, rg_sem, tau_est, tau_mean, tau_sem, &
           n_runs, cpu_time_n

     do irun = 1, n_runs
        write(11,'(I6,2X,I4,2X,I10,4(2X,F15.8))') &
              n_part, irun, seed_all(irun,isize), &
              rg_all(irun,isize), err_all(irun,isize), &
              tau_all(irun,isize), time_all(irun,isize)
     end do

     write(*,*)
     write(*,'(A,I5)')            'N = ', n_part
     write(*,'(A,F14.6,A,F12.6)') '  Average Rg   = ', rg_mean, &
                                  ' +/- ', rg_sem
     write(*,'(A,F14.4,A,F12.4)') '  Average tau  = ', tau_mean, &
                                  ' +/- ', tau_sem
     write(*,'(A,F14.4)')         '  tau estimate = ', tau_est

  end do

  close(10)
  close(11)

  write(*,*)
  write(*,*) '======================================'
  write(*,'(A,F12.2,A)') ' Total wall-clock time = ', wall_total, ' s'
  write(*,'(A,F12.2,A)') ' Total CPU time        = ', sum(time_all), ' s'
  write(*,'(A,F12.2)')   ' Parallel speed-up     = ', sum(time_all)/wall_total
  write(*,*) 'Averages written to scaling_results3.dat'
  write(*,*) 'Individual runs written to all_runs3.dat'
  write(*,*) '======================================'

end program polymer_scaling


!=================================================================
! MEAN AND STANDARD ERROR OF A SET OF INDEPENDENT VALUES
!
!   mean = (1/M) sum x_i
!
!   sem  = sqrt( sum (x_i-mean)^2 / (M(M-1)) )
!
! Because the runs are independent this needs no assumption
! about correlation times.
!=================================================================

subroutine mean_and_sem(m, x, xmean, xsem)

  use sim_para
  implicit none

  integer, intent(in) :: m
  real(dp), intent(in) :: x(m)
  real(dp), intent(out) :: xmean
  real(dp), intent(out) :: xsem

  integer :: i
  real(dp) :: s2

  xmean = 0.0d0

  do i = 1,m
     xmean = xmean + x(i)
  end do

  xmean = xmean/real(m,dp)

  if (m < 2) then
     xsem = 0.0d0
     return
  end if

  s2 = 0.0d0

  do i = 1,m
     s2 = s2 + (x(i)-xmean)**2
  end do

  xsem = sqrt(s2/(real(m,dp)*real(m-1,dp)))

end subroutine mean_and_sem


!=================================================================
! SEED THE RANDOM NUMBER GENERATOR
!
! The seed is spread over the whole seed array. Filling it with
! one repeated small integer can give a poor sequence.
!=================================================================

subroutine init_random_seed(iseed)

  use sim_para
  implicit none

  integer, intent(in) :: iseed

  integer :: nseed
  integer :: i
  integer, allocatable :: seed(:)

  call random_seed(size=nseed)
  allocate(seed(nseed))

  do i = 1,nseed
     seed(i) = iseed + 37*(i-1) + 1013*i*i
  end do

  call random_seed(put=seed)

  deallocate(seed)

  ! Throw away any leftover Box-Muller deviate from the last run
  gr_have_spare = .false.

end subroutine init_random_seed


!=================================================================
! RUN ONE POLYMER
!=================================================================

subroutine run_polymer(n_part, irun, rg_avg, rg_err, &
                       tau_meas, tau_est)

  use sim_para
  implicit none

  integer, intent(in) :: n_part
  integer, intent(in) :: irun
  real(dp), intent(out) :: rg_avg
  real(dp), intent(out) :: rg_err
  real(dp), intent(out) :: tau_meas
  real(dp), intent(out) :: tau_est

  integer :: step
  integer :: n_samples

  integer :: nsteps
  integer :: nequil
  integer :: sample_every
  integer :: max_mode_samples
  integer :: n_mode_samples

  integer :: max_neigh

  integer, allocatable :: no_neighb(:)
  integer, allocatable :: neigh_list(:,:)

  real(dp), allocatable :: pos(:,:)
  real(dp), allocatable :: pos_unwrap(:,:)
  real(dp), allocatable :: vel(:,:)
  real(dp), allocatable :: force(:,:)
  real(dp), allocatable :: pos_neigh(:,:)
  real(dp), allocatable :: mode1(:,:)

  real(dp) :: lx, ly, lz

  real(dp) :: tau_relax

  real(dp) :: rg2
  real(dp) :: mode_tmp(3)
  real(dp) :: rg2_sum
  real(dp) :: rg4_sum
  real(dp) :: rg2_avg
  real(dp) :: variance

  character(len=64) :: fname

  integer :: iu_trace

  logical :: rebuild
  logical :: check_neigh_list

  !---------------------------------------------------------------
  ! Box dimensions
  !
  ! A cubic box, large enough that the chain never sees its own
  ! periodic images.
  !---------------------------------------------------------------

  lx = 2.0d0*real(n_part,dp) + 20.0d0

  ly = lx
  lz = lx

  !---------------------------------------------------------------
  ! Run length from the estimated relaxation time
  !---------------------------------------------------------------

  tau_relax = gam*real(n_part,dp)**(1.0d0+2.0d0*nu) &
              /(3.0d0*pi*pi)

  nequil = nint(n_equil_tau*tau_relax/dt)

  nsteps = nequil + nint(n_prod_tau*tau_relax/dt)

  sample_every = max(1, nint(sample_tau*tau_relax/dt))

  !---------------------------------------------------------------
  ! Maximum number of neighbors
  !---------------------------------------------------------------

  max_neigh = n_part

  !---------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------

  allocate(pos(3,n_part))
  allocate(pos_unwrap(3,n_part))
  allocate(vel(3,n_part))
  allocate(force(3,n_part))

  allocate(no_neighb(n_part))
  allocate(neigh_list(max_neigh,n_part))
  allocate(pos_neigh(3,n_part))

  max_mode_samples = nsteps/max(1,sample_every) + 2
  allocate(mode1(3,max_mode_samples))
  n_mode_samples = 0

  !---------------------------------------------------------------
  ! Initialize polymer
  !---------------------------------------------------------------

  call initialize_polymer(n_part,lx,ly,lz, &
                          pos,pos_unwrap,vel)

  !---------------------------------------------------------------
  ! Build initial neighbor list
  !---------------------------------------------------------------

  call build_neigh_list(n_part,lx,ly,lz, &
                        pos,pos_unwrap, &
                        no_neighb,neigh_list, &
                        pos_neigh,max_neigh)

  !---------------------------------------------------------------
  ! Initial forces
  !---------------------------------------------------------------

  call calculate_forces(n_part,lx,ly,lz, &
                        pos,force, &
                        no_neighb,neigh_list,max_neigh)

  !---------------------------------------------------------------
  ! Initialize Rg averaging
  !
  ! Rg^2 is averaged, not Rg, because the scaling law is
  ! <Rg^2>^(1/2) ~ N^nu.
  !---------------------------------------------------------------

  rg2_sum = 0.0d0
  rg4_sum = 0.0d0
  n_samples = 0

  !---------------------------------------------------------------
  ! Time series of Rg, written from step 1 so that the
  ! equilibration can be checked by eye.
  !---------------------------------------------------------------

  write(fname,'(A,I0,A,I0,A)') 'Rg_trace_N3', n_part, '_r', irun, '.dat'

  open(newunit=iu_trace, file=trim(fname), status='replace')

  write(iu_trace,'(A)') '# time    Rg'

  !---------------------------------------------------------------
  ! Langevin dynamics
  !
  ! BAOAB:
  !
  ! B - half force kick
  ! A - half position update
  ! O - Langevin thermostat
  ! A - half position update
  ! B - half force kick with the new forces
  !---------------------------------------------------------------

  do step = 1,nsteps

     !------------------------------------------------------------
     ! B: first half kick
     !------------------------------------------------------------

     call half_kick(n_part,vel,force)

     !------------------------------------------------------------
     ! A O A
     !------------------------------------------------------------

     call langevin_step(n_part,lx,ly,lz, &
                        pos,pos_unwrap,vel)

     !------------------------------------------------------------
     ! Check neighbor list
     !------------------------------------------------------------

     rebuild = check_neigh_list(n_part, &
                                pos_unwrap,pos_neigh)

     if (rebuild) then

        call build_neigh_list(n_part,lx,ly,lz, &
                              pos,pos_unwrap, &
                              no_neighb,neigh_list, &
                              pos_neigh,max_neigh)

     end if

     !------------------------------------------------------------
     ! Calculate new forces
     !------------------------------------------------------------

     call calculate_forces(n_part,lx,ly,lz, &
                           pos,force, &
                           no_neighb,neigh_list,max_neigh)

     !------------------------------------------------------------
     ! B: second half kick
     !------------------------------------------------------------

     call half_kick(n_part,vel,force)

     !------------------------------------------------------------
     ! Sampling
     !
     ! The trace is written throughout, the average is
     ! accumulated only after equilibration.
     !------------------------------------------------------------

     if (mod(step,sample_every) == 0) then

        call calculate_rg(n_part,lx,ly,lz,pos,rg2)

        write(iu_trace,'(F15.4,2X,F15.8)') real(step,dp)*dt, sqrt(rg2)

        if (step > nequil) then

           rg2_sum = rg2_sum + rg2

           n_mode_samples = n_mode_samples + 1
           call calculate_rouse_mode1(n_part,pos_unwrap,mode_tmp)
           mode1(:,n_mode_samples) = mode_tmp

           rg4_sum = rg4_sum + rg2*rg2

           n_samples = n_samples + 1

        end if

     end if

  end do

  close(iu_trace)

  !---------------------------------------------------------------
  ! Final average Rg for this run
  !
  ! rg_err here is the naive within-run error. It assumes the
  ! samples are independent, which they are not, so it is a
  ! lower bound and is written to all_runs3.dat for reference
  ! only. The error actually quoted is the spread over runs.
  !---------------------------------------------------------------

  rg2_avg = rg2_sum / real(n_samples,dp)

  variance = rg4_sum/real(n_samples,dp) - rg2_avg*rg2_avg

  rg_avg = sqrt(rg2_avg)

  rg_err = 0.5d0*sqrt(max(variance,0.0d0)/real(n_samples,dp)) &
           /rg_avg

  call calculate_relaxation_time(n_part, irun, n_mode_samples, &
                                 sample_every*dt, mode1, tau_meas)

  tau_est = gam*real(n_part,dp)**(1.0d0+2.0d0*nu) &
            /(3.0d0*pi*pi)

  !---------------------------------------------------------------
  ! Deallocate
  !---------------------------------------------------------------

  deallocate(pos)
  deallocate(pos_unwrap)
  deallocate(vel)
  deallocate(force)

  deallocate(no_neighb)
  deallocate(neigh_list)
  deallocate(pos_neigh)
  deallocate(mode1)

end subroutine run_polymer


!=================================================================
! INITIALIZE POLYMER
!
! A straight chain is the slowest possible starting point, so the
! chain is grown as a self-avoiding random walk instead: each new
! bead is placed at distance r0 from the previous one, in a random
! direction, and the direction is redrawn if it overlaps a bead
! already placed.
!=================================================================

subroutine initialize_polymer(n_part,lx,ly,lz, &
                              pos,pos_unwrap,vel)

  use sim_para
  implicit none

  integer, intent(in) :: n_part

  real(dp), intent(in) :: lx,ly,lz

  real(dp), intent(out) :: pos(3,n_part)
  real(dp), intent(out) :: pos_unwrap(3,n_part)
  real(dp), intent(out) :: vel(3,n_part)

  integer, parameter :: max_try = 1000

  integer :: i
  integer :: j
  integer :: itry
  integer :: nfail

  real(dp) :: ux,uy,uz
  real(dp) :: unorm

  real(dp) :: trial(3)

  real(dp) :: dx,dy,dz
  real(dp) :: r2

  real(dp) :: vxcm
  real(dp) :: vycm
  real(dp) :: vzcm

  real(dp) :: gaussian_random

  logical :: ok

  nfail = 0

  !---------------------------------------------------------------
  ! First bead at the box center
  !---------------------------------------------------------------

  pos_unwrap(1,1) = lx/2.0d0
  pos_unwrap(2,1) = ly/2.0d0
  pos_unwrap(3,1) = lz/2.0d0

  !---------------------------------------------------------------
  ! Grow the rest of the chain
  !---------------------------------------------------------------

  do i = 2,n_part

     ok = .false.

     do itry = 1,max_try

        !---------------------------------------------------------
        ! Random direction on the unit sphere
        !---------------------------------------------------------

        ux = gaussian_random()
        uy = gaussian_random()
        uz = gaussian_random()

        unorm = sqrt(ux*ux + uy*uy + uz*uz)

        trial(1) = pos_unwrap(1,i-1) + r0*ux/unorm
        trial(2) = pos_unwrap(2,i-1) + r0*uy/unorm
        trial(3) = pos_unwrap(3,i-1) + r0*uz/unorm

        !---------------------------------------------------------
        ! Reject if it overlaps an earlier bead
        !---------------------------------------------------------

        ok = .true.

        do j = 1,i-2

           dx = trial(1) - pos_unwrap(1,j)
           dy = trial(2) - pos_unwrap(2,j)
           dz = trial(3) - pos_unwrap(3,j)

           r2 = dx*dx + dy*dy + dz*dz

           if (r2 < rc*rc) then

              ok = .false.

              exit

           end if

        end do

        if (ok) exit

     end do

     !------------------------------------------------------------
     ! If the growth failed the bead is placed anyway, but this
     ! is reported instead of happening silently: an overlapping
     ! bead gives a huge WCA force and can destabilize the start.
     !------------------------------------------------------------

     if (.not. ok) nfail = nfail + 1

     pos_unwrap(:,i) = trial(:)

  end do

  if (nfail > 0) then
     write(*,'(A,I6,A)') '  *** warning: chain growth failed for ', &
                         nfail, ' beads'
  end if

  !---------------------------------------------------------------
  ! Non-periodic simulation: no wrapping.
  !---------------------------------------------------------------

  pos = pos_unwrap

  !---------------------------------------------------------------
  ! Initial Maxwell-Boltzmann velocities
  !---------------------------------------------------------------

  do i = 1,n_part

     vel(1,i) = gaussian_random()*sqrt(kbT/mass)
     vel(2,i) = gaussian_random()*sqrt(kbT/mass)
     vel(3,i) = gaussian_random()*sqrt(kbT/mass)

  end do

  !---------------------------------------------------------------
  ! Remove center-of-mass velocity
  !---------------------------------------------------------------

  vxcm = sum(vel(1,:))/real(n_part,dp)
  vycm = sum(vel(2,:))/real(n_part,dp)
  vzcm = sum(vel(3,:))/real(n_part,dp)

  do i = 1,n_part

     vel(1,i) = vel(1,i) - vxcm
     vel(2,i) = vel(2,i) - vycm
     vel(3,i) = vel(3,i) - vzcm

  end do

end subroutine initialize_polymer


!=================================================================
! HALF FORCE KICK
!
! The B step of BAOAB. It is called twice per step: once before
! the drift with the old forces, once after with the new ones.
!=================================================================

subroutine half_kick(n_part,vel,force)

  use sim_para
  implicit none

  integer, intent(in) :: n_part

  real(dp), intent(inout) :: vel(3,n_part)

  real(dp), intent(in) :: force(3,n_part)

  integer :: i

  do i = 1,n_part

     vel(:,i) = vel(:,i) + &
                0.5d0*dt*force(:,i)/mass

  end do

end subroutine half_kick


!=================================================================
! LANGEVIN STEP
!
! The A O A part of BAOAB:
!
! A - half position update
! O - Langevin thermostat
! A - half position update
!
! The two B half kicks are done outside this routine.
!=================================================================

subroutine langevin_step(n_part,lx,ly,lz, &
                         pos,pos_unwrap,vel)

  use sim_para
  implicit none

  integer, intent(in) :: n_part

  real(dp), intent(in) :: lx,ly,lz

  real(dp), intent(inout) :: pos(3,n_part)
  real(dp), intent(inout) :: pos_unwrap(3,n_part)
  real(dp), intent(inout) :: vel(3,n_part)

  integer :: i

  real(dp) :: c1
  real(dp) :: c2

  real(dp) :: randx
  real(dp) :: randy
  real(dp) :: randz

  real(dp) :: gaussian_random

  !---------------------------------------------------------------
  ! Ornstein-Uhlenbeck coefficients
  !---------------------------------------------------------------

  c1 = exp(-gam*dt/mass)

  c2 = sqrt((kbT/mass)*(1.0d0-c1*c1))

  !---------------------------------------------------------------
  ! A: first half drift
  !---------------------------------------------------------------

  do i = 1,n_part

     pos_unwrap(:,i) = pos_unwrap(:,i) + &
                       0.5d0*dt*vel(:,i)

  end do

  !---------------------------------------------------------------
  ! O: Langevin thermostat
  !---------------------------------------------------------------

  do i = 1,n_part

     randx = gaussian_random()
     randy = gaussian_random()
     randz = gaussian_random()

     vel(1,i) = c1*vel(1,i) + c2*randx
     vel(2,i) = c1*vel(2,i) + c2*randy
     vel(3,i) = c1*vel(3,i) + c2*randz

  end do

  !---------------------------------------------------------------
  ! A: second half drift
  !---------------------------------------------------------------

  do i = 1,n_part

     pos_unwrap(:,i) = pos_unwrap(:,i) + &
                       0.5d0*dt*vel(:,i)

     pos(:,i) = pos_unwrap(:,i)

  end do

end subroutine langevin_step


!=================================================================
! BUILD NEIGHBOR LIST
!=================================================================

subroutine build_neigh_list(n_part,lx,ly,lz, &
                            pos,pos_unwrap, &
                            no_neighb,neigh_list, &
                            pos_neigh,max_neigh)

  use sim_para
  implicit none

  integer, intent(in) :: n_part
  integer, intent(in) :: max_neigh

  real(dp), intent(in) :: lx,ly,lz

  real(dp), intent(in) :: pos(3,n_part)
  real(dp), intent(in) :: pos_unwrap(3,n_part)

  integer, intent(out) :: no_neighb(n_part)

  integer, intent(out) :: neigh_list(max_neigh,n_part)

  real(dp), intent(out) :: pos_neigh(3,n_part)

  integer :: i
  integer :: j

  real(dp) :: dx
  real(dp) :: dy
  real(dp) :: dz
  real(dp) :: r2

  !---------------------------------------------------------------
  ! Reset neighbor counts
  !---------------------------------------------------------------

  no_neighb = 0

  !---------------------------------------------------------------
  ! Store positions at time of neighbor-list construction
  !---------------------------------------------------------------

  pos_neigh = pos_unwrap

  !---------------------------------------------------------------
  ! Search all pairs
  !---------------------------------------------------------------

  do i = 1,n_part-1

     do j = i+1,n_part

        dx = pos(1,j) - pos(1,i)
        dy = pos(2,j) - pos(2,i)
        dz = pos(3,j) - pos(3,i)

        r2 = dx*dx + dy*dy + dz*dz

        !---------------------------------------------------------
        ! Add pair if inside neighbor-list radius
        !
        ! Only i -> j is stored, since the force loop visits
        ! each pair once anyway.
        !---------------------------------------------------------

        if (r2 < rs*rs) then

           no_neighb(i) = no_neighb(i) + 1

           neigh_list(no_neighb(i),i) = j

        end if

     end do

  end do

end subroutine build_neigh_list


!=================================================================
! CHECK WHETHER NEIGHBOR LIST MUST BE REBUILT
!=================================================================

logical function check_neigh_list(n_part, &
                                  pos_unwrap,pos_neigh)

  use sim_para
  implicit none

  integer, intent(in) :: n_part

  real(dp), intent(in) :: pos_unwrap(3,n_part)
  real(dp), intent(in) :: pos_neigh(3,n_part)

  integer :: i

  real(dp) :: dx
  real(dp) :: dy
  real(dp) :: dz

  real(dp) :: max_disp2

  !---------------------------------------------------------------
  ! Maximum squared displacement since last neighbor-list build
  !---------------------------------------------------------------

  max_disp2 = 0.0d0

  do i = 1,n_part

     dx = pos_unwrap(1,i) - pos_neigh(1,i)
     dy = pos_unwrap(2,i) - pos_neigh(2,i)
     dz = pos_unwrap(3,i) - pos_neigh(3,i)

     max_disp2 = max(max_disp2, &
                     dx*dx + dy*dy + dz*dz)

  end do

  !---------------------------------------------------------------
  ! Rebuild when maximum displacement exceeds half the skin
  !---------------------------------------------------------------

  if (max_disp2 > (0.5d0*skin)**2) then

     check_neigh_list = .true.

  else

     check_neigh_list = .false.

  end if

end function check_neigh_list


!=================================================================
! FORCE CALCULATION
!
! Total force =
!
!       WCA excluded-volume force
!       +
!       harmonic bond force
!=================================================================

subroutine calculate_forces(n_part,lx,ly,lz, &
                            pos,force, &
                            no_neighb,neigh_list,max_neigh)

  use sim_para
  implicit none

  integer, intent(in) :: n_part
  integer, intent(in) :: max_neigh

  real(dp), intent(in) :: lx,ly,lz

  real(dp), intent(in) :: pos(3,n_part)

  real(dp), intent(out) :: force(3,n_part)

  integer, intent(in) :: no_neighb(n_part)

  integer, intent(in) :: neigh_list(max_neigh,n_part)

  integer :: i
  integer :: j
  integer :: p

  real(dp) :: dx
  real(dp) :: dy
  real(dp) :: dz

  real(dp) :: r2
  real(dp) :: r
  real(dp) :: rinv

  real(dp) :: fmag

  real(dp) :: fs

  !---------------------------------------------------------------
  ! Set all forces to zero
  !---------------------------------------------------------------

  force = 0.0d0

  !===============================================================
  ! WCA FORCES
  !===============================================================

  do i = 1,n_part

     do p = 1,no_neighb(i)

        j = neigh_list(p,i)

        dx = pos(1,i) - pos(1,j)
        dy = pos(2,i) - pos(2,j)
        dz = pos(3,i) - pos(3,j)

        r2 = dx*dx + dy*dy + dz*dz

        !---------------------------------------------------------
        ! WCA cutoff
        !---------------------------------------------------------

        if (r2 < rc*rc) then

           r = sqrt(r2)

           rinv = 1.0d0/r

           !----------------------------------------------------
           ! WCA force magnitude
           !
           ! F = 24 eps / r
           !     [2(sigma/r)^12 - (sigma/r)^6]
           !----------------------------------------------------

           fmag = 24.0d0*eps*rinv * &
                  (2.0d0*(sigma*rinv)**12 - &
                   (sigma*rinv)**6)

           !----------------------------------------------------
           ! Force on particle i
           !----------------------------------------------------

           force(1,i) = force(1,i) + &
                        fmag*dx*rinv

           force(2,i) = force(2,i) + &
                        fmag*dy*rinv

           force(3,i) = force(3,i) + &
                        fmag*dz*rinv

           !----------------------------------------------------
           ! Newton's third law
           !----------------------------------------------------

           force(1,j) = force(1,j) - &
                        fmag*dx*rinv

           force(2,j) = force(2,j) - &
                        fmag*dy*rinv

           force(3,j) = force(3,j) - &
                        fmag*dz*rinv

        end if

     end do

  end do

  !===============================================================
  ! HARMONIC SPRING FORCES
  !===============================================================

  do i = 1,n_part-1

     j = i+1

     !------------------------------------------------------------
     ! Vector from i to j
     !------------------------------------------------------------

     dx = pos(1,j) - pos(1,i)
     dy = pos(2,j) - pos(2,i)
     dz = pos(3,j) - pos(3,i)

     r2 = dx*dx + dy*dy + dz*dz

     r = sqrt(r2)

     !------------------------------------------------------------
     ! Harmonic force
     !
     ! F = kappa (r-r0) * r_vector/r
     !------------------------------------------------------------

     fs = kappa*(r-r0)/r

     force(1,i) = force(1,i) + fs*dx
     force(2,i) = force(2,i) + fs*dy
     force(3,i) = force(3,i) + fs*dz

     force(1,j) = force(1,j) - fs*dx
     force(2,j) = force(2,j) - fs*dy
     force(3,j) = force(3,j) - fs*dz

  end do

end subroutine calculate_forces


!=================================================================
! CALCULATE SQUARED RADIUS OF GYRATION
!
!              1
! Rg^2 = ----------- sum_i |r_i-r_CM|^2
!             N
!=================================================================

subroutine calculate_rg(n_part,lx,ly,lz,pos,rg2)

  use sim_para
  implicit none

  integer, intent(in) :: n_part
  real(dp), intent(in) :: lx,ly,lz
  real(dp), intent(in) :: pos(3,n_part)
  real(dp), intent(out) :: rg2

  integer :: i
  real(dp) :: com(3)
  real(dp) :: dr(3)

  ! Non-periodic box: use the actual coordinates directly.
  com = 0.0d0
  do i = 1,n_part
     com = com + pos(:,i)
  end do
  com = com/real(n_part,dp)

  rg2 = 0.0d0
  do i = 1,n_part
     dr = pos(:,i) - com
     rg2 = rg2 + sum(dr*dr)
  end do

  rg2 = rg2/real(n_part,dp)

end subroutine calculate_rg


!=================================================================
! FIRST ROUSE MODE
!
! X_1 = (1/N) sum_i cos[pi(i-1/2)/N] r_i
!
! The p=1 mode represents the slowest large-scale internal
! relaxation of the chain.
!=================================================================

subroutine calculate_rouse_mode1(n_part,pos,mode)

  use sim_para
  implicit none

  integer, intent(in) :: n_part
  real(dp), intent(in) :: pos(3,n_part)
  real(dp), intent(out) :: mode(3)

  integer :: i
  real(dp) :: q

  mode = 0.0d0

  do i = 1,n_part
     q = pi*(real(i,dp)-0.5d0)/real(n_part,dp)
     mode = mode + cos(q)*pos(:,i)
  end do

  mode = mode/real(n_part,dp)

end subroutine calculate_rouse_mode1


!=================================================================
! MEASURE RELAXATION TIME
!
! The normalized p=1 Rouse-mode autocorrelation is
!
! C(t) = <delta X_1(t0+t) . delta X_1(t0)> /
!        <|delta X_1|^2>
!
! For C(t) ~ exp(-t/tau), its integrated autocorrelation time
! approaches the exponential relaxation time.
!
! tau = dt_sample [1/2 + sum C(lag)]
!
! Note that a single run of 20 tau gives a noisy tau, biased a
! little low because the sum is cut at the first negative C.
! This is why tau is also averaged over the independent runs.
!
! The autocorrelation is written to Rouse_AC_N3*_r*.dat.
!=================================================================

subroutine calculate_relaxation_time(n_part,irun,n_samples, &
                                     dt_sample,mode,tau)

  use sim_para
  implicit none

  integer, intent(in) :: n_part
  integer, intent(in) :: irun
  integer, intent(in) :: n_samples
  real(dp), intent(in) :: dt_sample
  real(dp), intent(in) :: mode(3,n_samples)
  real(dp), intent(out) :: tau

  integer :: i
  integer :: lag
  integer :: count
  real(dp) :: mean_mode(3)
  real(dp) :: var0
  real(dp) :: c
  real(dp) :: sum_c
  real(dp) :: dx(3)
  integer :: iu_ac
  character(len=64) :: fname

  tau = 0.0d0

  if (n_samples < 10) return

  mean_mode = 0.0d0
  do i = 1,n_samples
     mean_mode = mean_mode + mode(:,i)
  end do
  mean_mode = mean_mode/real(n_samples,dp)

  var0 = 0.0d0
  do i = 1,n_samples
     dx = mode(:,i) - mean_mode
     var0 = var0 + sum(dx*dx)
  end do
  var0 = var0/real(n_samples,dp)

  if (var0 <= 1.0d-30) return

  write(fname,'(A,I0,A,I0,A)') 'Rouse_AC_N3',n_part,'_r',irun,'.dat'
  open(newunit=iu_ac,file=trim(fname),status='replace')
  write(iu_ac,'(A)') '# lag_time    normalized_C_Rouse_p1'

  sum_c = 0.0d0

  do lag = 0,n_samples-1

     c = 0.0d0
     count = n_samples-lag

     do i = 1,count
        dx = mode(:,i) - mean_mode
        c = c + sum(dx*(mode(:,i+lag)-mean_mode))
     end do

     c = c/real(count,dp)/var0

     write(iu_ac,'(F15.6,2X,F15.8)') real(lag,dp)*dt_sample,c

     if (lag > 0) then
        if (c <= 0.0d0) exit
        sum_c = sum_c + c
     end if

  end do

  close(iu_ac)

  tau = dt_sample*(0.5d0 + sum_c)

end subroutine calculate_relaxation_time


!=================================================================
! GAUSSIAN RANDOM NUMBER
!
! Box-Muller transformation
!
! Generates a normally distributed random number with:
!
! mean = 0
! variance = 1
!
! Box-Muller produces two deviates at a time. The second one is
! stored in the module and returned on the next call, which
! halves the cost. init_random_seed discards it so that a new
! run never starts with a number from the previous seed.
!=================================================================

real(8) function gaussian_random()

  use sim_para
  implicit none

  real(8) :: u1
  real(8) :: u2
  real(8) :: fac

  !---------------------------------------------------------------
  ! Return the stored deviate if there is one
  !---------------------------------------------------------------

  if (gr_have_spare) then

     gaussian_random = gr_spare

     gr_have_spare = .false.

     return

  end if

  !---------------------------------------------------------------
  ! Uniform random numbers
  !---------------------------------------------------------------

  call random_number(u1)
  call random_number(u2)

  ! Avoid log(0)

  if (u1 < 1.0d-12) then
     u1 = 1.0d-12
  end if

  !---------------------------------------------------------------
  ! Box-Muller transformation
  !---------------------------------------------------------------

  fac = sqrt(-2.0d0*log(u1))

  gaussian_random = fac*cos(2.0d0*pi*u2)

  gr_spare = fac*sin(2.0d0*pi*u2)

  gr_have_spare = .true.

end function gaussian_random
