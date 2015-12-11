!        Generated by TAPENADE     (INRIA, Ecuador team)
!  Tapenade 3.10 (r5717) - 30 Jul 2015 16:03
!
!  Differentiation of model in reverse (adjoint) mode:
!   gradient     of useful results: x y
!   with respect to varying inputs: x
SUBROUTINE LZ95_MODEL_ADJ(x, xb, yb, n, f)
  ! This subroutine computes the adjoint model of the Lorenz governing
  ! equations
  !
  ! dX[k] / dt = X[k] * (X[k + 1] - X[k - 2]) - X[k] + F
  !
  ! using Tapenade. The comments above the subroutine were generated by
  ! Tapenade. The auto-generated code was modified to suit the needs of the
  ! user.
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
  !     [1] Lorenz, E. N. and Emanuel K. A. (1998). Optimal Sites for
  !         Supplementary Weather Observations: Simulation with a Small Model.
  !         Journal of the Atmospheric Sciences.
  !
  !     [2] Hascoët, Laurent and Pascual, Valérie. (2013). The Tapenade
  !         Automatic Differentiation Tool: Principles, Model, and
  !         Specification.
  !
  ! Author:             Tapenade 3.10
  ! Edited by:          Jeremy Shaw
  ! Institution:        Portland State University
  ! Date Created:       18 November 2015
  ! Last Modified Date: 20 November 2015
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  DOUBLE PRECISION, DIMENSION(n), INTENT(IN) :: x
  DOUBLE PRECISION, DIMENSION(n), intent(out) :: xb
  DOUBLE PRECISION, DIMENSION(n), intent(in) :: yb
  DOUBLE PRECISION, INTENT(IN) :: f
  INTEGER :: i
  DOUBLE PRECISION :: tempb0
  DOUBLE PRECISION :: tempb
  tempb0 = x(n-1)*yb(n)
  xb(n-1) = xb(n-1) + (x(1)-x(n-2))*yb(n)
  xb(1) = xb(1) + tempb0
  xb(n-2) = xb(n-2) - tempb0
  xb(n) = xb(n) - yb(n)

  DO i=n-1,3,-1
    tempb = x(i-1)*yb(i)
    xb(i-1) = xb(i-1) + (x(i+1)-x(i-2))*yb(i)
    xb(i+1) = xb(i+1) + tempb
    xb(i-2) = xb(i-2) - tempb
    xb(i) = xb(i) - yb(i)

  END DO
  xb(1) = xb(1) + (x(3)-x(n))*yb(2)
  xb(3) = xb(3) + x(1)*yb(2)
  xb(n) = xb(n) - x(1)*yb(2)
  xb(2) = xb(2) - yb(2)

  xb(n) = xb(n) + (x(2)-x(n-1))*yb(1)
  xb(2) = xb(2) + x(n)*yb(1)
  xb(n-1) = xb(n-1) - x(n)*yb(1)
  xb(1) = xb(1) - yb(1)
END SUBROUTINE LZ95_MODEL_ADJ
