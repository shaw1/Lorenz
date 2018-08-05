subroutine lz95_forecast_tlm(x, xd, yd, n, F, dt)
  ! This subroutine integrates the state vector using the the tangent
  ! linear model of the RK4 integration of the Lorenz governing
  ! equations 
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
  ! Date Created:       18 November 2015
  ! Last Modified Date: 4 August 2018
  implicit none
  
  integer, intent(in) :: n
  double precision, dimension(n), intent(in) :: x
  double precision, dimension(n), intent(in) :: xd
  double precision, dimension(n), intent(out) :: yd
  double precision, intent(in) :: F
  double precision, intent(in) :: dt
  
  double precision, dimension(n) :: k1
  double precision, dimension(n) :: k1d
  double precision, dimension(n) :: k2
  double precision, dimension(n) :: k2d
  double precision, dimension(n) :: k3
  double precision, dimension(n) :: k3d
  double precision, dimension(n) :: k4d
  
  call lz95_model_tlm(x, xd, k1d, n, F)
  call lz95_model(x, k1, n, F)
  
  call lz95_model_tlm(x + dt/2.d0*k1, xd + dt*k1d/2.d0, k2d, n, F)
  call lz95_model(x + (dt / 2.d0) * k1, k2, n, F)
  
  call lz95_model_tlm(x + dt/2.d0*k2, xd + dt*k2d/2.d0, k3d, n, F)
  call lz95_model(x + (dt / 2.d0) * k2, k3, n, F)
  
  call lz95_model_tlm(x + dt*k3, xd + dt*k3d, k4d, n, F)
  
  yd = xd + dt * (k1d + 2.d0 * k2d + 2.d0 * k3d + k4d) / 6.d0
end subroutine lz95_forecast_tlm
