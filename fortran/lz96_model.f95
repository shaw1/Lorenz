subroutine lz96_model(x, y, b, c, F, h, J, K, n)
  ! This subroutine computes the right-hand-side of the governing
  ! equations for the Lorenz multi-scale model
  !
  ! dX[k] / dt = X[k - 1] * (X[k + 1] - X[k - 2]) - X[k] - (h * c / b)
  !                 * sum_{j = 1}^J Y[j, k] + F                    (1)
  ! dY[j, k] / dt = c * b * Y[j + 1, k] * (Y[j - 1, k] - Y[j + 2, k]) 
  !                 - c * Y[j, k] + (h * c / b) * X[k]             (2)
  !
  ! The state vector x is the collection of all variables X_k and Y_jk
  ! are ordered X_1, X_2, ..., X_K, Y_11, Y_21, ..., Y_J1, Y_12,
  ! ........, Y_JK.
  !
  ! Inputs:
  !     x - Input state vector
  !     b, c, F, h - Model parameters
  !     J - Number of small-scale Y-variables between two X-variables
  !     K - Number of X-variables
  !     n - Length of x
  !
  ! Output:
  !     y - Result of governing equations evaluated at x.
  !
  ! References:
  !     [1] Lorenz, E. N. (1996). Predictability: A problem partly
  !         solved. Proc. Seminar on predictability.
  !
  !     [2] Lorenz, E. N. and Emanuel K. A. (1998). Optimal Sites for
  !         Supplementary Weather Observations: Simulation with a
  !         Small Model. Journal of the Atmospheric Sciences.
  !
  ! Author:             Jeremy Shaw
  ! Institution:        Portland State University
  ! Date Created:       18 November 2015
  ! Last Modified Date: 10 November 2016
  implicit none
  
  integer, intent(in) :: n ! n = (J + 1) * K

  double precision, dimension(n), intent(in) :: x
  double precision, dimension(n), intent(out) :: y
  double precision, intent(in) :: b
  double precision, intent(in) :: c
  double precision, intent(in) :: F
  double precision, intent(in) :: h
  integer, intent(in) :: J
  integer, intent(in) :: K

  integer :: i
  integer :: kk
  integer :: jj
  double precision :: hc_over_b
  double precision :: c_times_b

  hc_over_b = h * c / b
  c_times_b = c * b

  ! Setup of dX_k/dt
  y(1) = x(K) * (x(2) - x(K - 1)) - x(1) - hc_over_b * &
       sum(x(K + 1 : K + J)) + F
  y(2) = x(1) * (x(3) - x(K)) - x(2) - hc_over_b * &
       sum(x(K + J + 1 : 2 * J + K)) + F

  do kk = 3, K - 1
     y(kk) = x(kk - 1) * (x(kk + 1) - x(kk - 2)) - x(kk) - hc_over_b &
          * sum(x((kk - 1) * J + K + 1 : kk * J + K)) + F
  end do
     
  y(K) = x(K - 1) * (x(1) - x(K - 2)) - x(K) - hc_over_b * &
       sum(x((K - 1) * J + K + 1 : (J + 1) * K)) + F

  ! Setup of dY_{j, k}/dt
  y(K + 1) = c_times_b * x(K + 2) * (x((J + 1) * K) - x(K + 3)) - c &
       * x(K + 1) + hc_over_b * x(1)

  do i = K + 2, (J + 1) * K - 2
     y(i) = c_times_b * x(i + 1) * (x(i - 1) - x(i + 2)) - c * x(i) &
          + hc_over_b * x((i - K - 1) / J + 1)
  end do

  jj = (J + 1) * K - 1
  y(jj) = c_times_b * x(jj + 1) * (x(jj - 1) - x(K + 1)) - c * x(jj) &
       + hc_over_b * x(K)
  y(jj + 1) = c_times_b * x(K + 1) * (x(jj) - x(K + 2)) - c * &
       x(jj + 1) + hc_over_b * x(K)
end subroutine lz96_model
