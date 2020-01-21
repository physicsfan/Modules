
MODULE radial
!***********************************************************************
      Use Grid
      Use orbs_csfs
   
!     Data structure for Radial Functions
      !> P, Q    -- large and small  components of the radial functions
      !> DP, DQ  -- derivatives of large and small component
      REAL(8), Dimension(:,:), allocatable :: P, Q, DP, DQ
      
!     Properties of radial functions
      !> e       -- diagonal energy parameter
      !> scf     -- scf_consistency
      REAL(8), Dimension(:), allocatable :: e,scf
      !> npt      -- number of points
      !> nodes   -- number of amplitudes
      !> sigma   -- orbital screening
      INTEGER, Dimension(:), allocatable :: pz, gamma,  npt, nodes, sigma

!      Nuclear properties
       REAL(kind=8), DIMENSION(:), allocatable :: zz
 
      Real(8), dimension(npX) :: F, yK 

   CONTAINS
      !> allocate memory for radial functions
      SUBROUTINE allocate_radials
         allocate(P(npX,nw),Q(npX,nw),DP(npX,nw), DQ(npX,nw))
         allocate(e(nw),scf(nw), npt(nw), nodes(nw), sigma(nw), zz(npX), pz(nw), gamma(nw))
      END SUBROUTINE allocate_radials


      !> Compute the numerical sDP, sDQ functions using a five-point formula
      !   for i= 3 ... n-2  
      SUBROUTINE  NumericalDP(ia)
         Implicit NONE 
         Integer, intent(in) :: ia
         Integer :: i, n 
         real(8) :: den
    
           
         DP(1:npX,ia) = 0.d0
         DQ(1:npX,ia) = 0.d0

         n = npt(ia)
         i = 2
         den = 1.d0/(2.d0*th)
         dp(i,ia) = (P(i+1,ia) - P(i-1,ia))*den
         dq(i,ia) = (Q(i+1,ia) - Q(i-1,ia))*den
         i = 3
         den = 1.d0/(th*12.d0)
         dp(i,ia) = (8.d0*(P(i+1,ia)-P(i-1,ia)  &
                         - P(i+2,ia)+P(i-2,ia)))*den
         dq(i,ia) = (8.d0*(Q(i+1,ia)-Q(i-1,ia)  &
                         - Q(i+2,ia)+Q(i-2,ia)))*den
         den = 1.d0/(60.d0*th) 
         DO I = 4, NPx-4
           dp(i,ia) = ( 45.d0*(P(i+1,ia)-P(i-1,ia)  &
                        -9.d0*(P(i+2,ia)-P(i-2,ia)) &
                          +   (P(i+3,ia)-P(i-3,ia))))*den
           dq(i,ia) = ( 45.d0*(Q(i+1,ia)-Q(i-1,ia)  &
                        -9.d0*(Q(i+2,ia)-Q(i-2,ia)) &
                          +   (Q(i+3,ia)-Q(i-3,ia))))*den
         end do

       END Subroutine NumericalDP
 
!======================================================================*
!   The value of this  function is the one-electron integral I (ia,ib) *
!   for  orbitals ia,ib. The analytical expression for this quantity   *
!   is given as  eq (9) in  I P Grant, B J McKenzie, P H Norrington,   *
!   D F Mayers, and N C Pyper,  Computer  Phys Commun 21 (1980) 211 .  *
!   The computed result is the sum of four types of integrands.
!                                                                      *
      Subroutine rinti(ia, ib, ans) 
!======================================================================*
      IMPLICIT NONE
      Integer, intent(in) :: ia, ib
      Real(8), intent(out):: ans
   
      Integer :: i,n, kappa 
      Real(8), DIMENSION(npx) :: Int1, Int2, Int3, Int4
      Real(8) :: c

     
      c = c_speed
      If (nak(ia) /= nak(ib)) then
         Write(99, '(A)') 'Error in Rinti -- kappa values not the same'
         stop
      else
        kappa = nak(ia)
      end if
  
      n = min(npt(ia), npt(ib))

      Int1 = 0.d0; Int2 = 0.d0; Int3 = 0.d0; Int4 = 0.d0 
      Int1 = q(:,ia)*dp(:,ib) - p(:,ia)*dq(:,ib)
      Int2 = (q(:,ia)*p(:,ib) + p(:,ia)*q(:,ib))
      Int3 =  q(:,ia)*q(:,ib)*r
      Int4 =  zz(:)*(p(:,ia)*p(:,ib) + q(:,ia)*q(:,ib))
      Print *, ' c = ', c
      Print *, ' quad(Int1) =', quad(int1)
      Print *, ' quad(Int2) =', quad(int2)
      Print *, ' quad(Int3) =', quad(int3)
      Print *, ' quad(Int4) =', quad(int4)
      ans  = c*(quad(Int1) + kappa*quad(Int2)) - 2.d0*c*c*quad(Int3) &
               + quad(Int4)

      End Subroutine rinti 
    
!=======================================================================
!               k
!       Stores Z (ia, ib, ) in the array YK.
!
!
     Subroutine zkf(k, ia, ic)
!-----------------------------------------------------------------------
       Integer, intent(in) :: ia, ic, k
       Integer :: i, m1, m2
       REAL(8) :: C1, C2, A, A2, A3, AI
!
       A  = (1.d0+H)**(-K)
       A2 = A*A
       A3 = A2*A
       AI = 1.d0/A

!       .. integrand for dt
       F =  r*(P(:,ia)*P(:,ic) + Q(:,ia)*Q(:,ic))
       yk(1) = 0.d0
       yk(2) = th*(F(1)+A*F(2))/2.d0

!      Simpson's rule
       i = 3
       yk(i) = th3*(F(i) + 4.d0*A*F(i-1) + A2*F(i-2))

!     Five point rule
      DO  i = 4,npX-1
        yK(i) = YK(i-2)*A2 + th90*(34.d0*(A2*F(i-2) + F(i)) + &
                   114.d0*A*F(i-1) -(A3*f(i-3) + AI*F(i+1))) 
      ENd do

      IF (IABS(ia-ic)  +  IABS(K) .NE. 0) return
!
!  *****  FOR Z0(I,I) SET THE LIMIT TO 1 AND REMOVE OSCILLATIONS
!  *****  INTRODUCED BY THE USE OF SIMPSON'S RULE
!
      M1 = (npX/2)*2 - 3              
      M2 = M1 - 1
      Print *, 'M1, M2 =', m1,m2
      C1 = 1.d0 - yk(M1)
      C2 = 1.d0 - yk(M2)
      
      Print *, 'C1, C2 =', c1, c2
      DO  i = 1,npX,2
        yk(i) = yk(i) + C1
      END DO
      DO  i = 2,npX,2
        yk(i) = yk(i) + C2
      END DO
      yk(npX) = yk(npX-2)
      yk(npX-1) = yk(npX-3)
     
      RETURN
      END SUBROUTINE ZKF
!
!     ------------------------------------------------------------------
!              Y K F
!     ------------------------------------------------------------------
!
!               k
!       Stores Y (i, j; r) in the array zk 
!       cid is a key for the last computed ykf with ykid
!
!
   SUBROUTINE YKF(k,ia,ic)
      Integer, intent(in) :: ia, ic, k
      Integer :: i, c, cid
      Real(8) :: A, A2, A3, AI
      Real(8), dimension(npX) :: error

!    test whether results are already in YK

      cid = ((k*key + ia)*key +ic)
      if ( cid .eq. ykid) return
!
      CALL ZKF(k, ia, ic)

      F=0.d0
      A  = (1.d0+H)**(-(K+1))
      A2 = A*A
      A3 = A2*A
      AI = 1.d0/A
      C = 2*K+1
      DO  i = npX-3,2,-1
        f(i) = th90*( 114.d0*A*YK(i+1) + 34.d0*(YK(i)+A2*YK(i+2)) &
                         - ai*YK(i-1) - a3*YK(i+3))
      End do
      f(npX-2:npX) = f(npX-3)
      ! Simpson's rule
      F(1) = th3*(YK(1) + 4.d0*A*YK(2) + A2*YK(3))
      DO i = npX-3,2,-1
        YK(i) = YK(i+2)*A2 + C*f(i)
      ENd do
      Yk(1) = Yk(3)*A2 + C*f(1)

      !  update the id of the saved result
      ykid = cid
      END subroutine ykf


!     Real(8) Function Slater(k,ia,ib,ic,id)
      Subroutine  Slater(k,ia,ib,ic,id,ans)
        IMPLICIT NONE
        Integer, intent(in) :: k,ia,ib,ic,id
        Real(8), intent(out) :: ans
        Integer :: i
        Real(8), dimension(npX) :: int, error

        call ykf(k, ia, ic)
        
        int = yk*p(:,ib)*p(:,id)
!       slater  = quad(int)
        ans  = quad(int)
      End subroutine  slater
END MODULE radial
