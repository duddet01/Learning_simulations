program he_importance
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer :: i, N
  real(dp) :: alpha
  real(dp) :: r1, r2
  real(dp) :: u1, u2, u3
  real(dp) :: cost1, cost2, phi1, phi2
  real(dp) :: beta, r12
  real(dp) :: sum, sum2, avg, err

  N = 10000000
  alpha = 2.0_dp
  sum = 0.0_dp
  sum2 = 0.0_dp

  do i = 1, N

     ! ---- sample r1 from p(r) = 4 a^3 r^2 e^{-2 a r}
     call random_number(u1)
     call random_number(u2)
     call random_number(u3)
     r1 = -(log(u1) + log(u2) + log(u3)) / (2.0_dp * alpha)

     ! ---- sample r2
     call random_number(u1)
     call random_number(u2)
     call random_number(u3)
     r2 = -(log(u1) + log(u2) + log(u3)) / (2.0_dp * alpha)

     ! ---- angles
     call random_number(cost1)
     cost1 = 2.0_dp*cost1 - 1.0_dp

     call random_number(cost2)
     cost2 = 2.0_dp*cost2 - 1.0_dp

     call random_number(phi1)
     phi1 = 2.0_dp * acos(-1.0_dp) * phi1

     call random_number(phi2)
     phi2 = 2.0_dp * acos(-1.0_dp) * phi2

     beta = cost1*cost2 + sqrt(1.0_dp-cost1**2) * &
            sqrt(1.0_dp-cost2**2) * cos(phi1 - phi2)

     r12 = sqrt(r1*r1 + r2*r2 - 2.0_dp*r1*r2*beta)

     sum  = sum  + 1.0_dp / r12
     sum2 = sum2 + (1.0_dp / r12)**2

  end do

  avg = sum / N
  err = sqrt( (sum2/N - avg**2) / N )

  print *, "Importance sampling result = ", avg
  print *, "Importance sampling error  = ", err

end program he_importance
