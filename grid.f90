
!====================================================================
   MODULE grid
!====================================================================
      Use case
!>  Parameter defining the maximum number of points in the grid
      INTEGER, PARAMETER :: npX = 600
!>  Parameter for the set-size and the first non-zero point
      REAL(kind=8), PARAMETER :: H = 1.d0/32, RNT= 1.d0/2**19,  &
                      th=log(1.d0+h),  th3 = th/3.d0, th90=th/90.d0
!>  Values of R, RR, and 1/R  for points on the gridA where
!!  r(1) = RNT/Z, r(i+1) = (1+H)r(i), i= 2, ngr.
      REAL(kind=8), DIMENSION(npX) :: R, RR, RM1, T

   CONTAINS

      SUBROUTINE  radial_grid
      INTEGER :: i

      r(1) =  rnt/z        ; t(1) = log(r(1))
      Do i = 2, npX
        r(i) = (1.d0+h)*r(i-1); t(i) = t(i-1) + th
      END Do
      rr = r*r; rm1=1.d0/r

      END SUBROUTINE radial_grid

!====================================================================
!                                                                      *
      REAL(8) FUNCTION QUAD(ARG, ib, c)
!                                                                      *
!   The result is an approximation  to the integral from R(1) to
!   infinity for the argumemts in the array ARG(:) for R(1) to inifity
!   using the three-point Simpson's rule over the logarithmnic grid.
!   Coef is a array for the intgrand from (0, r(1))  of the form
!       Int(r) = (c(3)*r + c(2))*r + c(1))*r*r
!
!====================================================================
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(8), DIMENSION(:), INTENT(IN) :: arg
      !>  beginning of the logarithmic grid
      Integer, intent(in), optional :: ib
      REAL(8),Dimension(*),intent(in), optional :: c
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(8) :: sum2, sum4, corr
      INTEGER :: i, n
!-----------------------------------------------
!
      Do i = npX, 100, -1
        if (dabs(arg(i)) > 1.d-16) exit
      end do
      n=i+1
      If (n > npX-2 ) then
        Write(0, *) ' Range indequate for the grid'
        STOP
      ENd if
     
      !  sum the even terms
      sum4 = 0.d0
      do i = 4, n, 2
        sum4 = sum4 +arg(i)
      end do
      !  sum the odd terms
      sum2 = 0.d0
      do i = 3, n, 2
        sum2 = sum2 +arg(i)
      end do
      quad = th*(4.d0*sum4 + 2.d0*sum2 - arg(1))/3.d0

      corr = 0.d0
      if (present(ib) .and. present(c))  then
!
!     Contribution from (0,r1) range of integration
!     Integral of  (((c1*r1 + c2)*r1 + c3)*r1)a

      corr = r(ib)**2*(c(1)/2.d0 + r(ib)*(c(2)/3.d0 + r(ib)*c(3)/4.d0))

        Print *, 'logarithmic =', quad
        Print *, 'Correction = ', corr
      End if

      quad = quad + corr

      RETURN
!
      END FUNCTION QUAD

   END MODULE grid
