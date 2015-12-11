subroutine lz95_model(x, y, n, F)
  ! This subroutine computes the right-hand-side of the governing equations
  ! for the Lorenz model
  !
  ! dX[k] / dt = X[k] * (X[k + 1] - X[k - 2]) - X[k] + F.
  !
  ! Inputs:
  !     x - Vector of X-variables
  !     n - Length of x
  !     F - Constant forcing term
  !
  ! Output:
  !     y - Result of governing equations evaluated at x.
  !
  ! Reference:
  !     Lorenz, E. N. and Emanuel K. A. (1998). Optimal Sites for
  !     Supplementary Weather Observations: Simulation with a Small Model.
  !     Journal of the Atmospheric Sciences.
  !
  ! Author:             Jeremy Shaw
  ! Institution:        Portland State University
  ! Date Created:       18 November 2015
  ! Last Modified Date: 20 November 2015
  implicit none

  integer, intent(in) :: n
  double precision, dimension(n), intent(in) :: x
  double precision, dimension(n), intent(out) :: y
  double precision, intent(in) :: F

  integer :: i

  y(1) = x(n) * (x(2) - x(n - 1)) - x(1) + F
  y(2) = x(1) * (x(3) - x(n)) - x(2) + F

  do  i = 3, n - 1
     y(i) = x(i - 1) * (x(i + 1) - x(i - 2)) - x(i) + F
  end do
  
  y(n) = x(n - 1) * (x(1) - x(n - 2)) - x(n) + F
end subroutine lz95_model
