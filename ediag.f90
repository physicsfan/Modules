!=======================================================================
!  This routine computes contributions to the diagonal energies arising
!  from the filled subshells (stype = 1) or filled-open interactions
!                                                                      *
   Module  diagonal_energies
      Use case
      Use orbs_csfs
      Use factorials
      Use integrals
!     
      Real(8) :: ecore
      Real(8), dimension(:), allocatable :: emt 
      Integer :: ia, ib, k, kmin, nkja, nkjb, qa, qb, &
                 coeff_label, tindex, vindex
 
!=======================================================================
 
   Contains

!>    Allocates memory for diagonal matrix elements and initializes the 
!     array of factorials needed in the calculation of matrix elements.
      Subroutine Initialize
      
      If (.not. allocated(emt)) then
        allocate(emt(ncsf))
        call fact
      end if
      end subroutine initialize 
     
!=======================================================================
!                                                                      *
      REAL(8) FUNCTION fco (K, IR, IA, IB) 
!                                                                      *
!   This routine evaluates a coefficient                               *
!                                                                      *
!                                K                                     *
!                               f   (IA,IB)                            *
!                                IR                                    *
!                                                                      *
!   Here  K  is the multipolarity, IR  is the sequence number of the   *
!   configuration, and  IA  and  IB  are orbital  sequence  numbers.   *
!   ( I P Grant,  B J McKenzie,  P H Norrington, D F Mayers, and N C   *
!   Pyper,  Computer Phys Commun 21 (1980) 207-231, Eqs (6). )         *
!                                                                      *
!   Written by Farid A Parpia, at Oxford  Last revision: 21 Dec 1992   *
!                                                                      *
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER  :: K,IR
      INTEGER , INTENT(IN) :: IA, IB
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: qa, qb
      REAL(8) :: FAC 
!-----------------------------------------------
!
      fco = 0.D0
      IF (IA == IB) THEN 
         qa = IQA (IA,IR)
         IF (K == 0) THEN 
            fco = DBLE((qa*(qa - 1))/2) 
         ELSE 
            qb = NKJ(IA) + 1 
!           .. check if the shell is filled
            IF (qa == qb) THEN 
               FAC = CLRX(nak(ia),K,nak(ia))*DBLE(qa) 
               fco = -0.5D0*FAC*FAC 
            ENDIF 
         ENDIF 
      ELSE 
         IF (K == 0) fco = DBLE (IQA(IA,IR)*IQA(IB,IR))
      ENDIF 
      END FUNCTION fco
      
!=======================================================================
!                                                                      *
      REAL(8) FUNCTION GCO (K, IR, IA, IB) 
!                                                                      *
!   This routine evaluates a coefficient                               *
!                                                                      *
!                                K                                     *
!                               g   (IA,IB)                            *
!                                IR                                    *
!                                                                      *
!                                                                      *
!   Here  K  is the multipolarity, IR  is the sequence number of the   *
!   configuration, and  IA and IB are orbital sequence  numbers. See   *
!   I P Grant,  B J McKenzie,  P H Norrington,  D F  Mayers, and N C   *
!   Pyper, Computer Phys Commun 21 (1980) 207-231, Eq (7).             *
!                                                                      *
!   Written by Farid A Parpia, at Oxford  Last revision: 18 Dec 1992   *
!                                                                      *
!----------------------------------------------------------------------
      INTEGER , INTENT(IN) :: k, ir, ia, ib
      INTEGER :: qa, qb 

      qa = IQA (IA,IR)
      qb = IQA (IB,IR)
      GCO = 0.0d0
!     .. if at least one shell is filled. 
!     .. RKCO needed for two open shells.
      IF (qa==NKJ(IA) + 1 .OR. qb==NKJ(IB)+1) THEN 
         GCO = -(qa*qb)*CLRX(NAK(IA),K,NAK(IB))**2 
      ENDIF 
      END FUNCTION GCO 

!----------------------------------------------------------------------
!>    Computing the common core

    Subroutine  Core_Energy
      Implicit none
 
!     .. energy of the core
      ecore = 0.d0
      DO IA = 1, ncore
        qa = IQA(ia,1) 
        if (qa == 0) cycle
        nkja = nkj(ia)
!         .. I(a, a) 
        coeff_label = ia*key+ia
        tindex = index(-1,coeff_label)
        ecore = ecore + qa*int_value(tindex)
!         .. F0(a,a)
        coeff_label = ((ia*key+ia)*key +ia)*key +ia
        vindex = index(0,coeff_label)
        ecore = ecore + 0.5*qa*qa*int_value(vindex)
!         .. Fk(a,a)
        Do k = 2, nkj(ia)-1, 2
          vindex = index(k,coeff_label)
          ecore = ecore + fco(k,1, ia,ia)* int_value(vindex)
        End do
        DO IB = IA+1, ncore
          qb = iqa(ib,1)
          if (qb == 0) cycle
!         .. F0(a,b)
          coeff_label = ((ia*key+ia)*key +ib)*key +ib
          vindex = index(0,coeff_label)
          ecore = ecore + qa*qb*int_value(vindex)
!          .. GK(a,b)
          nkjb= nkj(ib)
          IF (NAK(IA)*NAK(IB) > 0) THEN 
            KMIN = ABS((NKJA - NKJB)/2) 
          ELSE 
            KMIN = ABS((NKJA - NKJB)/2) + 1 
          ENDIF 
          do k = kmin, (nkja+nkjb)/2, 2
            coeff_label = ((ia*key+ib)*key +ia)*key +ib
            vindex = INDEX(k,coeff_label)
            IF(vindex == -1) cycle
            ecore = ecore - qa*qb*clrx(nak(ia),k,nak(ib))**2 &
                            *int_value(vindex)
          end do
         END DO 
      End do
    End subroutine core_energy

!----------------------------------------------------------------------
!   Diagonal energies for all blocks

    SUBROUTINE  Ediag
      Implicit none
      Integer :: j, sa, sb, qa, qb
      

    emt(1:ncsf) = ecore

      DO IA = 1, nw
        nkja = nkj(ia)
        Do IB = max(ncore+1,IA), nw
           nkjb = nkj(ib)
           DO j = 1,ncsf
             qa = iqa(ia, j);  if (qa == 0) cycle
             sa = stype(ia,j)
             qb = iqa(ib, j);  if (qb == 0) cycle
             sb = stype(ib,j)
!            .. only outside common core
             if ( ia  >  ncore ) then 
!              ..  shell outside the core
!              .. I(a, a) 
               coeff_label = ia*key+ia
               tindex = index(-1,coeff_label)
               emt(j) = emt(j) + qa*int_value(tindex)
!              .. F0(a,a)
               coeff_label = ((ia*key+ia)*key +ia)*key +ia
               vindex = index(0,coeff_label)
               emt(j) = emt(j) + 0.5d0*qa*qa*int_value(vindex)
!              .. Fk(a,a)
               Do k = 2, nkj(ia)-1, 2
                 vindex = index(k, coeff_label)
                 emt(j) = emt(j) + fco(k, j, ia, ia)*int_value(vindex)
               End do
!              .. F0(a,b)
               coeff_label = ((ia*key+ia)*key +ib)*key +ib
               vindex = INDEX(0,coeff_label)
               IF (vindex .NE. -1)
               emt(j) = emt(j) + qa*qb*int_value(vindex)
             end if
!             .. if either ia or ib is full
             IF (sa == 1 .OR. sb ==1) THEN
!              .. GK(a,b)
               IF (NAK(IA)*NAK(IB) > 0) THEN 
                 KMIN = ABS((NKJA - NKJB)/2) 
               ELSE 
                 KMIN = ABS((NKJA - NKJB)/2) + 1 
               ENDIF 
               DO k = kmin, (nkja+nkjb)/2, 2
                  coeff_label = ((ia*key+ib)*key +ia)*key +ib
                  vindex = INDEX(k,coeff_label)
                  IF(vindex == -1) cycle
                 emt(j) = emt(j) - qa*qb*clrx(nak(ia), k, nak(ib))**2 &
                      *int_value(vindex)
               end do
            END IF
         END DO
!         if (ia == ib) PRINT *, ia, ib, emt
         End do     
      End do
    
    END SUBROUTINE ediag

    SUBROUTINE print_diagonals

      k = 1; coeff_label = 2146435
      
      PRINT *, 'Index test:', k, coeff_label, INDEX(k,coeff_label)

    END SUBROUTINE print_diagonals

  End module diagonal_energies
