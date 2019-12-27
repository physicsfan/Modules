   MODULE grid
!====================================================================
!>  Parameter defining the maximum number of points in the grid
      INTEGER, PARAMETER :: npX=400
!>  Parameter for the set-size and the first non-zero point
      REAL(kind=8), PARAMETER :: H = 1.d0/128.d0,  th=log(1+h), &
                                 th3 = th/3.d0
!>  Values of R, RR, and 1/R  for points on the grid  where
!!  s(i) = (1/Z) e^{(i-1)th); r(i) = s(i)-(1/Z); dr = s dt.
      REAL(kind=8), DIMENSION(npX) :: T, S, R, RR, RM1, SM1
      REAL(kind=8) :: Z 

   CONTAINS

      SUBROUTINE  radial_grid
      INTEGER :: i

      t(1) = 0.d0; s(1) = 1.d0/Z ;        r(1) = 0.d0; rm1(1) = 0.d0
      t(2) = th;   s(2) = (1.d0 + h)/Z;   r(2) = h/Z;  rm1(2) = Z/h
      
      Do i = 2, npX
        t(i) = t(i-1) + th
        s(i) = s(i-1)*(1+h)
        r(i) = s(i) - 1.d0/Z
      END Do
 
      rr = r*r
      rm1(1)=0.d0;   rm1(2:npX)= 1.d0/r(2:NPX); sm1 = 1.d0/s
       
      END SUBROUTINE radial_grid

!====================================================================
!                                                                      *
      REAL(kind=8) FUNCTION QUAD(ARG)
!                                                                      *
!   The result is an approximation  to the integral from R(1) to
!   infinity for the argumemts in the array ARG(:) for R(1) to inifity
!   using the three-point Simpson's rule over the logarithmnic grid.
!   Coef is a array for the intgrand from (0, r(1))  of the form
!                                                                     *
!====================================================================
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(kind=8), DIMENSION(:), INTENT(IN) :: ARG
      !>  beginning of the logarithmic grid
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(kind=8) :: sum2, sum4
      INTEGER :: i, n 
!-----------------------------------------------
!
      Do i = npX, 100, -1
        if (dabs(arg(i)) > 1.d-16) exit
      end do
      n=i
      Print *, 'npt = ', n
     
      !  sum the even terms
      sum4 = 0.d0
      do i = 2, n, 2
        sum4 = sum4 +arg(i)
      end do
      !  sum the odd terms
      sum2 = 0.d0
      do i = 1, n, 2
        sum2 = sum2 +arg(i)
      end do
      quad = th3*(4.d0*sum4 + 2.d0*sum2 - arg(1))
      
      Print *, 'Returning from quad'
      RETURN  
!
      END FUNCTION QUAD 
   END MODULE grid
