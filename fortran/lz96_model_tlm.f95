subroutine lz96_model_tlm(x, xd, yd, b, c, F, h, J, K, n)
  ! This subroutine computes the tangent linear model of the governing
  ! equations for the Lorenz multi-scale model,
  !
  ! dX[k] / dt = X[k - 1] * (X[k + 1] - X[k - 2]) - X[k] - (h * c / b)
  !                 * sum_{j = 1}^J Y[j, k] + F                    (1)
  ! dY[j, k] / dt = c * b * Y[j + 1, k] * (Y[j - 1, k] - Y[j + 2, k]) 
  !                 - c * Y[j, k] + (h * c / b) * X[k]             (2)
  !
  ! Inputs:
  !     x - Input state vector
  !     xd - Direction vector
  !     b, c, F, h - Model parameters
  !     J - Number of small-scale Y-variables between two X-variables
  !     K - Number of X-variables
  !     n - Length of x, which equals (J + 1) * K
  !
  ! Output:
  !     yd - The tangent linear model evaluated at x in the direction
  !          of xd.
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
  ! Last Modified Date: 4 August 2018
  implicit none
  
  integer, intent(in) :: n
  double precision, dimension(n), intent(in) :: x
  double precision, dimension(n), intent(in) :: xd
  double precision, dimension(n), intent(out) :: yd
  double precision, intent(in) :: b
  double precision, intent(in) :: c
  double precision, intent(in) :: F
  double precision, intent(in) :: h
  integer, intent(in) :: J
  integer, intent(in) :: K
  
  integer :: i
  integer :: kk
  double precision :: hc_over_b
  double precision :: c_times_b

  hc_over_b = h * c / b
  c_times_b = c * b

  ! Tangent linear model for dX[k] / dt
  yd(1) = xd(K) * (x(2) - x(K - 1)) + x(K) * (xd(2) - xd(K - 1)) - &
&         xd(1) - hc_over_b * sum(xd(K + 1 : K + J))

  yd(2) = xd(1) * (x(3) - x(K)) + x(1) * (xd(3) - xd(K)) - xd(2) - &
&         hc_over_b * sum(xd(K + J + 1 : 2 * J + K))

  do kk=3,k-1
     yd(kk) = xd(kk - 1) * (x(kk + 1) - x(kk - 2)) + x(kk - 1) * &
&             (xd(kk + 1) - xd(kk - 2)) - xd(kk) - hc_over_b * &
&             sum(xd((kk - 1) * J + K + 1 : kk * J + K))
  end do
  
  yd(K) = xd(K - 1) * (x(1) - x(K - 2)) + x(K - 1) * (xd(1) - &
&         xd(K - 2)) - xd(K) - hc_over_b * &
&         sum(xd((K - 1) * J + K + 1 : n))

  ! Tangent linear model for dY[j, k] / dt
  yd(K + 1) = c_times_b * (xd(K + 2) * (x(n) - x(K + 3)) + x(K + 2) *&
&            (xd(n) - xd(K + 3))) - c * xd(K + 1) + hc_over_b * xd(1)
  
  do i = K + 2, n - 2
     yd(i) = c_times_b * (xd(i + 1) * (x(i - 1) - x(i + 2)) + x(i + 1)&
&            * (xd(i - 1) - xd(i + 2))) - c * xd(i) + hc_over_b * &
&            xd((i - K - 1) / J + 1)
  end do

  yd(n - 1) = c_times_b * (xd(n) * (x(n - 2) - x(K + 1)) + x(n) * &
&             (xd(n - 2) - xd(K + 1))) - c * xd(n - 1) + hc_over_b &
&             * xd(K)

  yd(n) = c_times_b * (xd(K + 1) * (x(n - 1) - x(K + 2)) + x(K + 1) *&
&         (xd(n - 1) - xd(K + 2))) - c * xd(n) + hc_over_b * xd(k)
end subroutine lz96_model_tlm
