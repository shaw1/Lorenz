subroutine lz96_forecast(x, y, b, c, F, h, J, K, n, dt)
  ! This subroutine performs the fourth-order Runge-Kutta integration (RK4)
  ! scheme to the Lorenz multi-scale model.
  !
  ! Inputs:
  !     x - Input state vector
  !     b, c, F, h - Model parameters
  !     J - Number of small-scale Y-variables between two X-variables
  !     K - Number of X-variables
  !     n - Length of x
  !     dt - Time step-size
  !
  ! Output:
  !     y - Result of integrating x using RK4.
  !
  ! References:
  !     [1] Lorenz, E. N. (1996). Predictability: A problem partly solved. Proc.
  !         Seminar on predictability.
  !     [2] Lorenz, E. N. and Emanuel K. A. (1998). Optimal Sites for
  !         Supplementary Weather Observations: Simulation with a Small Model.
  !         Journal of the Atmospheric Sciences.
  !
  ! Author:             Jeremy Shaw
  ! Institution:        Portland State University
  ! Date Created:       18 November 2015
  ! Last Modified Date: 20 November 2015
  implicit none

  integer, intent(in) :: n
  double precision, dimension(n), intent(in) :: x
  double precision, dimension(n), intent(out) :: y
  double precision, intent(in) :: b
  double precision, intent(in) :: c
  double precision, intent(in) :: F
  double precision, intent(in) :: h
  double precision, intent(in) :: dt
  integer, intent(in) :: J
  integer, intent(in) :: K

  double precision, dimension(n) :: k1
  double precision, dimension(n) :: k2
  double precision, dimension(n) :: k3
  double precision, dimension(n) :: k4

  call lz96_model(x, k1, b, c, F, h, J, K, n)
  call lz96_model(x + (dt / 2.D0) * k1, k2, b, c, F, h, J, K, n)
  call lz96_model(x + (dt / 2.D0) * k2, k3, b, c, F, h, J, K, n)
  call lz96_model(x + dt * k3, k4, b, c, F, h, J, K, n)

  y = x + (dt / 6.D0) * (k1 + 2.D0 * k2 + 2.D0 * k3 + k4)
end subroutine lz96_forecast
