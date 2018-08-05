subroutine lz95_forecast_adj(x, yb, xb, n, F, dt)
  ! This subroutine integrates the state vector using the adjoint
  ! model of the RK4 integration of the Lorenz governing equations
  !
  ! dX[k] / dt = X[k] * (X[k + 1] - X[k - 2]) - X[k] + F
  !
  ! over one time-step.
  !
  ! Inputs:
  !     x - Vector of X-variables
  !     xd - Direction vector
  !     n - Length of x
  !     F - Constant forcing term
  !     dt - Time-step
  !
  ! Output:
  !     yd - The tangent model evaluated at x in the direction of yb.
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
  ! Date Created:       10 November 2016
  ! Last Modified Date: 4 August 2018
  implicit none
  
  integer, intent(in) :: n
  double precision, dimension(n), intent(in) :: x
  double precision, dimension(n), intent(in) :: yb
  double precision, dimension(n), intent(out) :: xb
  double precision, intent(in) :: F
  double precision, intent(in) :: dt

  double precision, dimension(n) :: x2
  double precision, dimension(n) :: x3
  double precision, dimension(n) :: k1
  double precision, dimension(n) :: k1hat
  double precision, dimension(n) :: k2
  double precision, dimension(n) :: k2hat
  double precision, dimension(n) :: k3
  double precision, dimension(n) :: k3hat
  double precision, dimension(n) :: k4hat

  call lz95_model(x, k1, n, F)
  x2 = x + (dt / 2.d0) * k1
  call lz95_model(x2, k2, n, F)
  x3 = x + (dt / 2.d0) * k2
  call lz95_model(x3, k3, n, F)

  call lz95_model_adj(x + dt * k3, yb, k4hat, n, F)
  call lz95_model_adj(x3, yb + (dt / 2.d0) * k4hat, k3hat, n, F)
  call lz95_model_adj(x2, yb + (dt / 2.d0) * k3hat, k2hat, n, F)
  call lz95_model_adj(x, yb + dt * k2hat, k1hat, n, F)
  
  xb = yb + (dt/6.d0) * (k1hat + 2.d0 * k2hat + 2.d0 * k3hat + k4hat)
end subroutine lz95_forecast_adj
