!
!======================================================================
!                                                                      *
    MODULE factorials
!     
      Implicit none
      INTEGER, PARAMETER :: MFACT = 500 
      REAL(8), DIMENSION(MFACT) :: GAM 

    CONTAINS
!
!----------------------------------------------------------------------
      SUBROUTINE FACT 
!                                                                      *
!   Calculates the logs  of factorials 
!                                                                      *
!   Written by N S Scott                    Last update: 15 Oct 1992   *
!                                                                      *
!-----------------------------------------------
      INTEGER :: I 
      REAL(8) :: X 
!
      GAM(1) = 1.0D00 ; GAM(2) = 1.0D00 ; X = 2.0D00 
!
      DO I = 3, 30 
         GAM(I) = GAM(I-1)*X 
         X = X + 1.0D00 
      END DO 
!
      DO I = 1, 30 
         GAM(I) = LOG(GAM(I)) 
      END DO 
!
      X = 3.0D01 
!
      DO I = 31, MFACT 
         GAM(I) = GAM(I-1) + LOG(X) 
         X = X + 1.0D00 
      END DO 
!
      END SUBROUTINE FACT 
!- --------------------------------------------------------------------
!                                                                      *
      REAL(8) FUNCTION CLRX (KAPPAA, K, KAPPAB) 
!                                                                      *
!   The value of CLRX is the 3-j symbol:                               *
!                                                                      *
!                    ( JA        K        JB  )                        *
!                    ( 1/2       0       -1/2 )                        *
!                                                                      *
!   The  K'S are kappa angular quantum numbers. The formula is taken   *
!   from D M Brink and G R Satchler, <Angular Momentum>, second edi-   *
!   tion (Oxford: Clarendon press, 1968), p 138.   The logarithms of   *
!   the first  MFACT  factorials must be available in  COMMON/FACTS/   *
!   for this program to function correctly. Note that  N!  is stored   *
!   in FACT(N+1)                                                       *
!                                                                      *
!   Written by Farid A Parpia, at Oxford   Last updated: 06 Oct 1992   *
!                                                                      *
!------------------------------------------------------------------------
      INTEGER(1), INTENT(IN) :: KAPPAA,  KAPPAB 
      INTEGER   , intent(in) :: k
      INTEGER :: KA, KB, KAPKB, KABKP, KAMKB, KBMKA 
      REAL(8) :: EXPTRM 
!-----------------------------------------------
!
!   Determine the absolute values of the kappas
      KA = ABS(KAPPAA) 
      KB = ABS(KAPPAB) 
!
!   Perform the triangularity check
      IF (ABS(KA - KB)<=K .AND. KA+KB-1>=K) THEN 
!
!   Triangularity satisfied; compute the 3j coefficient
!
!   Begin with the logarithm of the square of the leading term
         EXPTRM = -LOG(DBLE(KA*KB)) 
!
!   Compute the logarithm of the square root of the leading term
!   and the factorial part that doesn't depend on the parity of
!   KA+KB+K (the delta factor)
!
         KAPKB = KA + KB 
         KABKP = KAPKB + K 
         KAMKB = KA - KB 
         KBMKA = KB - KA 
         EXPTRM = 0.5D00*(EXPTRM + GAM(KAPKB-K)+GAM(KAMKB+K+1)+GAM(KBMKA+K+1)-&
            GAM(KABKP+1)) 
!
!   The remainder depends on the parity of KA+KB+K
         IF (MOD(KABKP,2) == 0) THEN 
!
! Computation for even parity case
!   Include the phase factor: a minus sign if necessary
            IF (MOD(3*KABKP/2,2) == 0) THEN 
               CLRX = 1.0D00 
            ELSE 
               CLRX = -1.0D00 
            ENDIF 
!
!   Include the contribution from the factorials
            EXPTRM = EXPTRM + GAM((KABKP+2)/2) - GAM((KAPKB-K)/2) - GAM((KAMKB+&
               K+2)/2) - GAM((KBMKA+K+2)/2) 
!
         ELSE 
!
! Computation for odd parity case
!   Include the phase factor: a minus sign if necessary
            IF (MOD((3*KABKP - 1)/2,2) == 0) THEN 
               CLRX = 1.0D00 
            ELSE 
               CLRX = -1.0D00 
            ENDIF 
!
!   Include the contribution from the factorials
            EXPTRM = EXPTRM + GAM((KABKP+1)/2) - GAM((KAPKB-K+1)/2) - GAM((&
               KAMKB+K+1)/2) - GAM((KBMKA+K+1)/2) 
!
         ENDIF 
!
! Final assembly
         CLRX = CLRX*EXP(EXPTRM) 
!
      ELSE 
!
! Triangularity violated; set the coefficient to zero
         CLRX = 0.0D00 
!
      ENDIF 
!
      END FUNCTION CLRX 
    END MODULE factorials
