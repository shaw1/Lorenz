subroutine lz96_forecast_tlm(x, xd, yd, b, c, F, h, J, K, n, dt)
  ! This subroutine computes the tangent linear model of the RK4
  ! integration of the governing equations for the Lorenz multi-scale
  ! model
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
  !     xd - Direction vector
  !     b, c, F, h - Model parameters
  !     J - Number of small-scale Y-variables between two X-variables
  !     K - Number of X-variables
  !     n - Length of x
  !     dt - Time-step
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
  double precision, intent(in) :: dt
  integer, intent(in) :: J
  integer, intent(in) :: K
  
  double precision, dimension(n) :: k1
  double precision, dimension(n) :: k1d
  double precision, dimension(n) :: k2
  double precision, dimension(n) :: k2d
  double precision, dimension(n) :: k3
  double precision, dimension(n) :: k3d
  double precision, dimension(n) :: k4d

  call lz96_model_tlm(x, xd, k1d, b, c, F, h, J, K, n)
  call lz96_model(x, k1, b, c, F, h, J, K, n)
  
  call lz96_model_tlm(x + dt / 2.d0 * k1, xd + dt * k1d / 2.d0, k2d, &
       b, c, F, h, J, K, n)
  call lz96_model(x + (dt / 2.d0) * k1, k2, b, c, F, h, J, K, n)
  
  call lz96_model_tlm(x + dt / 2.d0 * k2, xd + dt * k2d / 2.d0, k3d, &
&      b, c, F, h, J, K, n)
  call lz96_model(x + (dt / 2.d0) * k2, k3, b, c, F, h, J, K, n)
  
  call lz96_model_tlm(x + dt * k3, xd + dt * k3d, k4d, b, c, F, h, J,&
&      K, n)
  
  yd = xd + dt * (k1d + 2.d0 * k2d + 2.d0 * k3d + k4d) / 6.d0
end subroutine lz96_forecast_tlm
