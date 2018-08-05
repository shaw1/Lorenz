subroutine lz95_model_adj(x, yb, xb, n, F)
  ! This subroutine computes the adjoint model of the Lorenz governing
  ! equations.
  !
  ! dX[k] / dt = X[k] * (X[k + 1] - X[k - 2]) - X[k] + F
  !
  ! Inputs:
  !     x - Vector of X-variables
  !     yb - Direction vector
  !     n - Length of x
  !     F - Constant forcing term
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
  double precision, intent(in) :: F
  
  integer :: i

  xb(1) = x(n - 1) * yb(n) - yb(1) + (x(3) - x(n)) * yb(2) - x(2) * &
&         yb(3)
  xb(2) = x(n) * yb(1) - yb(2) + (x(4) - x(1)) * yb(3) - x(3) * yb(4)
  
  do i = 3, n - 1
    xb(i) = x(i - 2) * yb(i - 1) - yb(i) + (x(i + 2) - x(i - 1)) * &
&           yb(i + 1) - x(i + 1) * yb(i + 2)
 end do

 xb(n - 1) = x(n - 3) * yb(n - 2) - yb(n - 1) + (x(1) - x(n - 2)) * &
&            yb(n) - x(n) * yb(1)

  xb(n) = x(n - 2) * yb(n - 1) - yb(n) + (x(2) - x(n - 1)) * yb(1) - &
&         x(1) * yb(2)
end subroutine lz95_model_adj
