subroutine lz95_model_tlm(x, xd, yd, n, F)
  ! This subroutine computes the tangent linear model of the Lorenz
  ! governing equations.
  !
  ! dX[k] / dt = X[k] * (X[k + 1] - X[k - 2]) - X[k] + F
  !
  ! Inputs:
  !     x - Vector of X-variables
  !     xd - Direction vector
  !     n - Length of x
  !     F - Constant forcing term
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
  double precision, intent(in) :: F
  
  integer :: i
  
  yd(1) = x(n) * (xd(2) - xd(n - 1)) + xd(n) * (x(2) - x(n - 1)) &
&         - xd(1)
  yd(2) = x(1) * (xd(3) - xd(n)) + xd(1) * (x(3) - x(n)) - xd(2)
  
  do i = 3, n - 1
     yd(i) = x(i - 1) * (xd(i + 1) - xd(i - 2)) + xd(i - 1) * &
&            (x(i + 1) - x(i - 2)) - xd(i)
  end do
 
  yd(n) = x(n - 1) * (xd(1) - xd(n - 2)) + xd(n - 1) * (x(1) - &
&         x(n - 2)) - xd(n)
end subroutine lz95_model_tlm
