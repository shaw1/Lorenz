subroutine lz96_model_adj(x, yb, xb, b, c, F, h, J, K, n)
  ! This subroutine computes the adjoint model of the governing
  ! equations for the Lorenz multi-scale model.
  !
  ! dX[k] / dt = X[k - 1] * (X[k + 1] - X[k - 2]) - X[k] - (h * c / b)
  !                 * sum_{j = 1}^J Y[j, k] + F                    (1)
  ! dY[j, k] / dt = c * b * Y[j + 1, k] * (Y[j - 1, k] - Y[j + 2, k]) 
  !                 - c * Y[j, k] + (h * c / b) * X[k]             (2)
  !
  ! Inputs:
  !     x - Input state vector
  !     yb - Direction vector
  !     b, c, F, h - Model parameters
  !     J - Number of small-scale Y-variables between two X-variables
  !     K - Number of X-variables
  !     n - Length of x
  !
  ! Output:
  !     xb - The adjoint model evaluated at x in the direction of yb.
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
  double precision, dimension(n), intent(in) :: yb
  double precision, dimension(n), intent(out) :: xb
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

  ! Adjoint for dX[k] / dt
  xb(1) = x(K - 1) * yb(K) - yb(1) + (x(3) - x(K)) * yb(2) - x(2) * & 
&         yb(3) + hc_over_b * sum(yb(K + 1 : J + K))

  xb(2) = x(K) * yb(1) - yb(2) + (x(4) - x(1)) * yb(3) - x(3) * &
&         yb(4) + hc_over_b * sum(yb(J + K + 1 : 2 * J + K))

  do kk = 3, K - 2
     xb(kk) = x(kk - 2) * yb(kk - 1) - yb(kk) + (x(kk + 2) - x(kk - 1))&
&             * yb(kk + 1) - x(kk + 1) * yb(kk + 2) + hc_over_b * &
&             sum(yb((kk - 1) * J + K + 1 : kk * J + K))
  end do
  
  xb(K - 1) = x(K - 3) * yb(K - 2) - yb(K - 1) + (x(1) - x(K - 2)) * &
&             yb(K) - x(K) * yb(1) + hc_over_b * &
&             sum(yb(n - 2 * J + 1 : n - J))
  
  xb(K) = x(K - 2) * yb(K - 1) - yb(K) + (x(2) - x(K - 1)) * &
&         yb(1) - x(1) * yb(2) + hc_over_b * sum(yb(n - J + 1 : n))

  ! Adjoint for dY[j, k] / dt
  xb(K + 1) = c_times_b * ((x(n - 1) - x(K + 2)) * yb(n) - x(n) * &
&             yb(n - 1) + x(K + 3) * yb(K + 2)) - c * yb(K + 1) - &
&             hc_over_b * yb(1)

  xb(K + 2) = c_times_b * ((x(n) - x(K + 3)) * yb(K + 1) - x(K + 1) *&
&             yb(n) + x(K + 4) * yb(K + 3)) - c * yb(K + 2) - &
&             hc_over_b * yb(1)

  do i = K + 3, n - 2
     xb(i) = c_times_b * ((x(i - 2) - x(i + 1)) * yb(i - 1) - x(i - 1)&
&            * yb(i - 2) + x(i + 2) * yb(i + 1)) - c * yb(i) - &
&            hc_over_b * yb((i - K - 1) / J + 1)
  end do

  xb(n - 1) = c_times_b * ((x(n - 3) - x(n)) * yb(n - 2) - x(n - 2) *&
&             yb(n - 3) + x(K + 1) * yb(n)) - c * yb(n - 1) - &
&             hc_over_b * yb(K)
  
  xb(n) = c_times_b * ((x(n - 2) - x(K + 1)) * yb(n - 1) - x(n - 1) *&
&         yb(n - 2) + x(K + 2) * yb(K + 1)) - c * yb(n) - hc_over_b *&
&         yb(K)
end subroutine lz96_model_adj
