subroutine lz95_forecast(x, y, n, F, dt)
  ! This subroutine performs the fourth-order Runge-Kutta integration
  ! (RK4) scheme to the Lorenz 40-variable model.
  !
  ! Inputs:
  !     x - Vector of X-variables
  !     n - Length of x
  !     F - Constant forcing term
  !     dt - Time step-size
  !
  ! Output:
  !     y - Result of integrating x using RK4.
  !
  ! Reference:
  !     Lorenz, E. N. and Emanuel K. A. (1998). Optimal Sites for
  !     Supplementary Weather Observations: Simulation with a Small
  !     Model. Journal of the Atmospheric Sciences.
  !
  ! Author:             Jeremy Shaw
  ! Institution:        Portland State University
  ! Date Created:       18 November 2015
  ! Last Modified Date: 10 November 2016
  implicit none

  integer, intent(in) :: n
  double precision, dimension(n), intent(in) :: x
  double precision, dimension(n), intent(out) :: y
  double precision, intent(in) :: F
  double precision, intent(in) :: dt

  double precision, dimension(n) :: k1
  double precision, dimension(n) :: k2
  double precision, dimension(n) :: k3
  double precision, dimension(n) :: k4

  call lz95_model(x, k1, n, F)
  call lz95_model(x + (dt / 2.D0) * k1, k2, n, F)
  call lz95_model(x + (dt / 2.D0) * k2, k3, n, F)
  call lz95_model(x + dt * k3, k4, n, F)

  y = x + (dt / 6.D0) * (k1 + 2.D0 * k2 + 2.D0 * k3 + k4)
end subroutine lz95_forecast
