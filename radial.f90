
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
      REAL(8), DIMENSION(:), ALLOCATABLE :: e,scf, pz
      !> npt      -- number of points
      !> nodes   -- number of amplitudes
      !> sigma   -- orbital screening
      INTEGER, Dimension(:), allocatable :: gamma,  npt, nodes, sigma

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
      SUBROUTINE  Slater(k,ia,ib,ic,id,ans)
        IMPLICIT NONE
        Integer, intent(in) :: k,ia,ib,ic,id
        Real(8), intent(out) :: ans
        Integer :: i
        Real(8), dimension(npX) :: int, error

        call ykf(k, ia, ic)
        
        int = yk*p(:,ib)*p(:,id)
!       slater  = quad(int)
       ans  = quad(int)
     END SUBROUTINE  slater

!=======================================================================
!                                                                      *
     SUBROUTINE write_radials  
!                                                                      *
!   Write all subshell radial wavefunctions in binary format to        *
!   "rwfn.out" file.                                                   *
!                                                                      *
!=======================================================================
       IMPLICIT NONE
       INTEGER :: i, j, ierr
  
       OPEN (rwfnout, FILE='rwfn.out', FORM='unformatted', STATUS='replace', iostat=ierr)  

       IF (ierr == 1) THEN 
          WRITE (istde, *) 'Error when opening ''rwfn.out'''
          STOP  
       ENDIF
      
       !   Binary file output
       WRITE (rwfnout) 'G92RWF'
       DO J = 1, NW
          WRITE (rwfnout) INT(NP(J)), INT(NAK(J)), E(J), NPT(j) 
          WRITE (rwfnout) PZ(J), (P(I,J),I=1,NPT(j)), (Q(I,J),I=1,NPT(j))
          WRITE (rwfnout) (R(I),I=1,NPT(J)) 
       END DO
       
       CLOSE(rwfnout) 

    END SUBROUTINE write_radials


!======================================================================
!
    SUBROUTINE load_radials
!                                                                      *
!                                                                      *
!   This subroutine loads radial wavefunctions from the rwfn.inp file. *
!                                                                      *
!   Call(s) to: [Library]: INTRPQ                                      *
!======================================================================  
      IMPLICIT NONE
  
      ! local variables
      INTEGER :: i, ierr, ios, j, k, npty
      INTEGER :: npy, naky
      REAL(kind=8) :: ey, pzy, dnorm, accy
      CHARACTER(len=6) :: g92rwf
      REAL(kind=8), DIMENSION(:), ALLOCATABLE :: py, qy, ry
  
      accy = h**6
  
      OPEN (rwfnin, FILE='rwfn.inp', FORM='unformatted', STATUS='old', iostat=ierr)  
          
      ! Check the file; if not as expected, try again              
      READ (rwfnin) g92rwf 
      IF (G92RWF/='G92RWF') stop  'This is not a Radial WaveFunction File;' 

      ! Read orbital information from Read Orbitals File   
      DO i = 1, nw
         READ (rwfnin) npy, naky, ey, npty
         
         ALLOCATE(py(npty),qy(npty),ry(npty))
         
         READ (rwfnin) pzy, py(1:npty), qy(1:npty)
         READ (rwfnin) ry(1:npty)
         IF (ABS(r(2)-ry(2)).GT.accy) THEN
            e(i) = ey
            pz(i) = pzy
            CALL INTRPQ (py, qy, npty, ry, i, DNORM)
         ELSE
            e(i) = ey
            pz(i) = pzy
            p(:,i) = py(:)
            q(:,i) = qy(:)
         END IF
         !   Determine the effective maximum tabulation point
         k=npx
         DO
            IF (ABS(p(k,i)) .GT. 1.d-16) EXIT
            k = k - 1
         END DO
         npt(i) = k
         DEALLOCATE(py,qy,ry)
         call numericalDP(i)
      END DO
        
      CLOSE(rwfnin)
      
    END SUBROUTINE load_radials



  !***********************************************************************
!                                                                      *
    SUBROUTINE INTRPQ(PA, QA, MA, RA, J, DNORM) 
!                                                                      *
!   This  subprogram  interpolates  the  arrays  PA(1:MA), QA(1:MA),   *
!   tabulated on grid RA(1:MA) into the COMMON arrays PF(1:MF(J),J),   *
!   QF(1:MF(J),J). (Aitken's  algorithm is used. See F B Hildebrand,   *
!   Introduction  to  Numerical  Analysis, 2nd ed., McGraw-Hill, New   *
!   York, NY, 1974.) The orbital is renormalized.                      *
!                                                                      *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 14 Oct 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  15:37:49   2/14/04  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17

      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: MA 
      INTEGER  :: J 
      REAL(kind=8), INTENT(OUT) :: DNORM 
      REAL(kind=8), DIMENSION(ma), INTENT(IN) :: PA, QA, RA
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: MXORD = 13 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, MFJ, NRSTLO, KOUNT, IROW, K, NRSTHI, LLO, LHI, LOCNXT, &
           ILIROK, ILDIAG, ILOTHR, n
      REAL(kind=8) :: RAMA, accy
      INTEGER, DIMENSION(NNNW) :: MF       
      REAL(kind=8), DIMENSION(MXORD) :: X, DX 
      REAL(kind=8), DIMENSION((MXORD*(MXORD + 1))/2) :: POLYP, POLYQ 
      REAL(kind=8) :: XBAR, PESTL, QESTL, DIFF, DIFFT, DXKMN1, DXIROW, &
           FACTOR, PESTT, QESTT, DPBP, DQBQ, DNFAC
      REAL(kind = 8), DIMENSION(npx) :: ta
      LOGICAL :: SET 
      LOGICAL, DIMENSION(npx) :: USED 
!-----------------------------------------------
!
!
!   This function dispenses with the need for a two-dimensional
!   array
!
!
!   Initialization
!
      ACCY = 1.0d-08 
!
!   This is always true in GRASP
!
      P(1,J) = 0.0D00 
      Q(1,J) = 0.0D00 
!
!   Determine end of grid
!
      mfj = npx 
      DO
         IF (R(mfj) <= RA(MA)) EXIT
         mfj = mfj - 1
      END DO
      
!          PRINT *, 'npx:', npx
!          PRINT *, ' ma:', ma
!          PRINT *, 'mfj:', mfj
!          PRINT *, 'RA(MA):', ra(ma)
!          PRINT *, 'RA(MA+1):', ra(ma+1)
!          PRINT *, 'R(mfj):', r(mfj)
!          PRINT *, 'R(mfj+1):', r(mfj+1)
!
!   Overall initialization for interpolation
!
      NRSTLO = 0 
      KOUNT = 0
!
!   Perform interpolation
!
      DO I = 2, MFJ 
!
!   Initialization for interpolation
!
         XBAR = R(I) 
         IROW = 0 
         PESTL = 0.0D00 
         QESTL = 0.0D00 
!
!   Determine the nearest two grid points bounding the present
!   grid point
!
         DO
            k = nrstlo + 1
            IF (RA(k) >= XBAR) EXIT
            nrstlo = k
         END DO
         nrsthi = k
         
!             PRINT *, '   xbar:', xbar
!             PRINT *, 'nrstlo:', nrstlo
!             PRINT *, 'nrsthi:', nrsthi
!
!   Clear relevant piece of use-indicator array
!
         LLO = MAX(NRSTLO - MXORD,1) 
         LHI = MIN(NRSTHI + MXORD,MA) 
         USED(LLO:LHI) = .FALSE. 
!
!   Determine next-nearest grid point
!
4        CONTINUE 
         IROW = IROW + 1 
         LLO = MAX(NRSTLO - IROW + 1,1) 
         LHI = MIN(NRSTHI + IROW - 1,MA) 
         SET = .FALSE. 
         DO K = LLO, LHI 
            IF (USED(K)) CYCLE  
            IF (.NOT.SET) THEN 
               DIFF = RA(K) - XBAR 
               LOCNXT = K 
               SET = .TRUE. 
            ELSE 
               DIFFT = RA(K) - XBAR 
               IF (ABS(DIFFT) < ABS(DIFF)) THEN 
                  DIFF = DIFFT 
                  LOCNXT = K 
               ENDIF
            ENDIF
         END DO
         USED(LOCNXT) = .TRUE. 
         X(IROW) = RA(LOCNXT) 
         DX(IROW) = DIFF 
!
!   Fill table for this row
!
         DO K = 1, IROW 
            ILIROK = (irow*(irow - 1))/2 + k
            IF (K == 1) THEN 
               POLYP(ILIROK) = PA(LOCNXT) 
               POLYQ(ILIROK) = QA(LOCNXT) 
            ELSE 
               ILDIAG = ((k - 1)*(k - 2))/2 + k - 1 
               ILOTHR = (irow*(irow - 1))/2 + k - 1 
               DXKMN1 = DX(K-1) 
               DXIROW = DX(IROW) 
               FACTOR = 1.0D00/(X(IROW)-X(K-1)) 
               POLYP(ILIROK) = (POLYP(ILDIAG)*DXIROW-POLYP(ILOTHR)*DXKMN1)*&
                    FACTOR 
               POLYQ(ILIROK) = (POLYQ(ILDIAG)*DXIROW-POLYQ(ILOTHR)*DXKMN1)*&
                    FACTOR 
            ENDIF
         END DO
!
!   Check for convergence
!
         ILDIAG = (irow*(irow - 1))/2 + irow 
         PESTT = POLYP(ILDIAG) 
         QESTT = POLYQ(ILDIAG)
         IF (PESTT==0.0D00 .OR. QESTT==0.0D00) THEN 
            IF (IROW < MXORD) THEN 
               GO TO 4 
            ELSE 
               P(I,J) = PESTT 
               Q(I,J) = QESTT 
            ENDIF
         ELSE 
            DPBP = ABS((PESTT - PESTL)/PESTT) 
            DQBQ = ABS((QESTT - QESTL)/QESTT) 
            IF (DQBQ<ACCY .AND. DPBP<ACCY) THEN 
               P(I,J) = PESTT 
               Q(I,J) = QESTT 
            ELSE 
               PESTL = PESTT 
               QESTL = QESTT 
               IF (IROW < MXORD) THEN 
                  GO TO 4 
               ELSE 
                  P(I,J) = PESTT 
                  Q(I,J) = QESTT 
                  KOUNT = KOUNT + 1 
               ENDIF
            ENDIF
         ENDIF
!
      END DO
!
!   Ensure that all points of the array are defined by setting the
!   tail to zero
! 
      P(MFJ+1:npx,J) = 0.0D00 
      Q(MFJ+1:npx,J) = 0.0D00
!
!
!   Normalization
!
      ta = r*(p(:,j)*p(:,j) + q(:,j)*q(:,j))           
      DNORM = quad(ta)
      DNFAC = 1.0D00/SQRT(DNORM)
!          PRINT *, 'dnorm', dnorm
!          PRINT *, 'dnfac', dnfac
      P(:MFJ,J) = P(:MFJ,J)*DNFAC 
      Q(:MFJ,J) = Q(:MFJ,J)*DNFAC 
!
      RETURN  
      
    END SUBROUTINE intrpq
    
END MODULE radial
