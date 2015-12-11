!        Generated by TAPENADE     (INRIA, Ecuador team)
!  Tapenade 3.10 (r5717) - 30 Jul 2015 16:03
!
!  Differentiation of lz96_model in reverse (adjoint) mode:
!   gradient     of useful results: x y
!   with respect to varying inputs: x
SUBROUTINE LZ96_MODEL_ADJ(x, xb, yb, b, c, f, h, j, k, n)
  ! This subroutine computes the adjoint model of the governing equations
  ! for the Lorenz multi-scale model
  !
  !    dX[k] / dt = X[k - 1] * (X[k + 1] - X[k - 2]) - X[k] - (h * c / b)
  !                 * sum_{j = 1}^J Y[j, k] + F                       (1)
  ! dY[j, k] / dt = c * b * Y[j + 1, k] * (Y[j - 1, k] - Y[j + 2, k]) - 
  !                 c * Y[j, k] + (h * c / b) * X[k]                  (2)
  !
  ! using Tapenade. The comments above the subroutine were generated by
  ! Tapenade. The auto-generated code was modified to suit the needs of the
  ! user.
  !
  ! Inputs:
  !     x - Input state vector
  !     yb - Direction vector
  !     b, c, F, h - Model parameters
  !     J - Number of small-scale Y-variables between two X-variables
  !     K - Number of X-variables
  !     n - Length of x
  !
  ! Output:
  !     xb - The adjoint model evaluated at x in the direction of yb.
  !
  ! References:
  !     [1] Lorenz, E. N. (1996). Predictability: A problem partly solved.
  !         Proc. Seminar on predictability.
  !     [2] Lorenz, E. N. and Emanuel K. A. (1998). Optimal Sites for
  !         Supplementary Weather Observations: Simulation with a Small Model.
  !         Journal of the Atmospheric Sciences.
  !     [3] Hascoët, Laurent and Pascual, Valérie. (2013). The Tapenade
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
  DOUBLE PRECISION, INTENT(IN) :: b
  DOUBLE PRECISION, INTENT(IN) :: c
  DOUBLE PRECISION, INTENT(IN) :: f
  DOUBLE PRECISION, INTENT(IN) :: h
  INTEGER, INTENT(IN) :: j
  INTEGER, INTENT(IN) :: k
  INTEGER :: i
  INTEGER :: kk
  INTEGER :: jj
  DOUBLE PRECISION :: hc_over_b
  DOUBLE PRECISION :: c_times_b
  INTRINSIC SUM
  DOUBLE PRECISION :: tempb4
  DOUBLE PRECISION :: tempb3
  DOUBLE PRECISION :: tempb2
  DOUBLE PRECISION :: tempb1
  DOUBLE PRECISION :: tempb0
  DOUBLE PRECISION :: tempb
  hc_over_b = h*c/b
  c_times_b = c*b

  jj = (j+1)*k - 1
  tempb3 = c_times_b*x(k+1)*yb(jj+1)
  xb(k+1) = xb(k+1) + c_times_b*(x(jj)-x(k+2))*yb(jj+1)
  xb(jj) = xb(jj) + tempb3
  xb(k+2) = xb(k+2) - tempb3
  xb(jj+1) = xb(jj+1) - c*yb(jj+1)
  xb(k) = xb(k) + hc_over_b*yb(jj+1)

  tempb4 = c_times_b*x(jj+1)*yb(jj)
  xb(jj+1) = xb(jj+1) + c_times_b*(x(jj-1)-x(k+1))*yb(jj)
  xb(jj-1) = xb(jj-1) + tempb4
  xb(k+1) = xb(k+1) - tempb4
  xb(jj) = xb(jj) - c*yb(jj)
  xb(k) = xb(k) + hc_over_b*yb(jj)

  DO i=(j+1)*k-2,k+2,-1
    tempb2 = c_times_b*x(i+1)*yb(i)
    xb(i+1) = xb(i+1) + c_times_b*(x(i-1)-x(i+2))*yb(i)
    xb(i-1) = xb(i-1) + tempb2
    xb(i+2) = xb(i+2) - tempb2
    xb(i) = xb(i) - c*yb(i)
    xb((i-k-1)/j+1) = xb((i-k-1)/j+1) + hc_over_b*yb(i)

  END DO
  tempb0 = c_times_b*x(k+2)*yb(k+1)
  xb(k+2) = xb(k+2) + c_times_b*(x((j+1)*k)-x(k+3))*yb(k+1)
  xb((j+1)*k) = xb((j+1)*k) + tempb0
  xb(k+3) = xb(k+3) - tempb0
  xb(k+1) = xb(k+1) - c*yb(k+1)
  xb(1) = xb(1) + hc_over_b*yb(k+1)

  tempb1 = x(k-1)*yb(k)
  xb(k-1) = xb(k-1) + (x(1)-x(k-2))*yb(k)
  xb(1) = xb(1) + tempb1
  xb(k-2) = xb(k-2) - tempb1
  xb(k) = xb(k) - yb(k)
  xb((k-1)*j+k+1:(j+1)*k) = xb((k-1)*j+k+1:(j+1)*k) - hc_over_b*yb(k)

  DO kk=k-1,3,-1
    tempb = x(kk-1)*yb(kk)
    xb(kk-1) = xb(kk-1) + (x(kk+1)-x(kk-2))*yb(kk)
    xb(kk+1) = xb(kk+1) + tempb
    xb(kk-2) = xb(kk-2) - tempb
    xb(kk) = xb(kk) - yb(kk)
    xb((kk-1)*j+k+1:kk*j+k) = xb((kk-1)*j+k+1:kk*j+k) - hc_over_b*yb(kk)

  END DO
  xb(1) = xb(1) + (x(3)-x(k))*yb(2)
  xb(3) = xb(3) + x(1)*yb(2)
  xb(k) = xb(k) - x(1)*yb(2)
  xb(2) = xb(2) - yb(2)
  xb(k+j+1:2*j+k) = xb(k+j+1:2*j+k) - hc_over_b*yb(2)

  xb(k) = xb(k) + (x(2)-x(k-1))*yb(1)
  xb(2) = xb(2) + x(k)*yb(1)
  xb(k-1) = xb(k-1) - x(k)*yb(1)
  xb(1) = xb(1) - yb(1)
  xb(k+1:k+j) = xb(k+1:k+j) - hc_over_b*yb(1)
END SUBROUTINE LZ96_MODEL_ADJ
