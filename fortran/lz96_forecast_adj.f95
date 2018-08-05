subroutine lz96_forecast_adj(x, yb, xb, b, c, F, h, J, K, n, dt)
  ! This subroutine computes the adjoint model of the RK4 integration
  ! of the governing equations for the Lorenz multi-scale model
  !
  ! dX[k] / dt = X[k - 1] * (X[k + 1] - X[k - 2]) - X[k] - (h * c / b)
  !                 * sum_{j = 1}^J Y[j, k] + F                    (1)
  ! dY[j, k] / dt = c * b * Y[j + 1, k] * (Y[j - 1, k] - Y[j + 2, k]) 
  !                 - c * Y[j, k] + (h * c / b) * X[k]             (2)
  !
  ! over one time-step.
  !
  ! Inputs:
  !     x - Input state vector
  !     yb - Direction vector
  !     b, c, F, h - Model parameters
  !     J - Number of small-scale Y-variables between two X-variables
  !     K - Number of X-variables
  !     n - Length of x
  !     dt - Time-step
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
  double precision, intent(in) :: dt
  integer, intent(in) :: J
  integer, intent(in) :: K

  double precision, dimension(n) :: x2
  double precision, dimension(n) :: x3
  double precision, dimension(n) :: k1
  double precision, dimension(n) :: k1hat
  double precision, dimension(n) :: k2
  double precision, dimension(n) :: k2hat
  double precision, dimension(n) :: k3
  double precision, dimension(n) :: k3hat
  double precision, dimension(n) :: k4hat

  call lz96_model(x, k1, b, c, F, h, J, K, n)
  x2 = x + (dt / 2.d0) * k1
  call lz96_model(x2, k2, b, c, F, h, J, K, n)
  x3 = x + (dt / 2.d0) * k2
  call lz96_model(x3, k3, b, c, F, h, J, K, n)

  call lz96_model_adj(x + dt * k3, yb, k4hat, b, c, F, h, J, K, n)
  call lz96_model_adj(x3, yb + (dt / 2.d0) * k4hat, k3hat, b, c, F, &
&                     h, J, K, n)
  call lz96_model_adj(x2, yb + (dt / 2.d0) * k3hat, k2hat, b, c, F, &
&                     h, J, K, n)
  call lz96_model_adj(x, yb + dt * k2hat, k1hat, b, c, F, h, J, K, n)

  xb = yb + (dt/6.d0) * (k1hat + 2.d0 * k2hat + 2.d0 * k3hat + k4hat)
end subroutine lz96_forecast_adj
