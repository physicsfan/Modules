MODULE racah
  USE orbs_csfs
  use gg
  
  IMPLICIT NONE

  INTEGER, DIMENSION(:), ALLOCATABLE :: NQ1, NQ2 
  INTEGER, DIMENSION(:), ALLOCATABLE :: JJC1, JJC2
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: JJQ1, JJQ2 
  INTEGER, DIMENSION(:), ALLOCATABLE :: JLIST, KLIST 
  INTEGER, DIMENSION(:), ALLOCATABLE :: JLIS, JC1S, JC2S


  !> (from buffer_C)
  INTEGER :: NBDIM, NVCOEF 
  INTEGER, DIMENSION(:,:), ALLOCATABLE  :: label
  REAL(kind=8), DIMENSION(:), ALLOCATABLE :: coeff

    
  
CONTAINS

  SUBROUTINE allocate_racah
    ALLOCATE(NQ1(nw), NQ2(nw),JJC1(nw), JJC2(nw))
    ALLOCATE(JJQ1(3,nw), JJQ2(3,nw)) 
    ALLOCATE(JLIST(nw), KLIST(nw), JLIS(nw), JC1S(nw), JC2S(nw))
    nq1=zero;nq2=zero;jjc1=zero;jjc2=zero;jjq1=zero;jjq2=zero
    jlist=zero;klist=zero;jlis=zero;jc1s=zero;jc2s=zero
  END SUBROUTINE allocate_racah

  SUBROUTINE deallocate_racah
    DEALLOCATE(NQ1,NQ2,JJC1,JJC2,JJQ1,JJQ2) 
    DEALLOCATE(JLIST,KLIST,JLIS,JC1S,JC2S)
  END SUBROUTINE deallocate_racah
  
!*******************************************************************
!                                                                  *
  SUBROUTINE ONESCALAR(JA,JB,IA1,IA2,VSHELL)
!                                                                  *
!   The main program for evaluating the reduced matrix elements of *
!   a one particle operator for configurations in jj-coupling.     *
!                                                                  *
!   The last modification made by G. Gaigalas       October  2017  *
!                                                                  *
!*******************************************************************
!
    IMPLICIT NONE    
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
    INTEGER, INTENT(IN)       :: JA,JB
    INTEGER, INTENT(OUT)      :: IA1,IA2
    REAL(kind=8), INTENT(OUT) :: VSHELL(NW)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
    INTEGER :: KK,IOPAR,IDQ,IJ,JA1,JA2,NDQ,NS,ISH,I,II,I1,IM, &
         JW,KS1,KS2,NX,NPEELM, k
    INTEGER, DIMENSION(2) :: IS,KS
    REAL(kind=8) ::  TCOEFF
    !-----------------------------------------------
    IA1 = 0
    KK = 1
    IOPAR = 1
    IF(ICHKQ1(JA,JB).EQ.0)RETURN

!   Generate the  arrays  defining  the quantum numbers of the 
!   states involved in the  matrix  element  linking  configurations 
!   labelled by JA, JB.                           
    CALL SETQNA (JA,JB)
    
    VSHELL = ZERO

!   Analyse peel shell interactions
    IDQ = 0
    JA1 = 0
    JA2 = 0
    IF (NPEEL .NE. 0) THEN
       DO JW = 1,NPEEL
          IJ = JLIST(JW)
          NDQ = NQ1(IJ)-NQ2(IJ)
          IF (IABS (NDQ) .GT. 1) RETURN
          IF (NDQ .GT. 0) THEN
             JA1 = JW
             IDQ = IDQ + 1
          ELSEIF (NDQ .LT. 0) THEN
             JA2 = JW
             IDQ = IDQ + 1
          ENDIF
       END DO
       
       IF (IDQ .GT. 2) RETURN
       
       !   Evaluate the array VSHELL
       !
       !   Then there are two possibilities IDQ = 0 or IDQ = 2
       !   if IDQ = 0, then loop over all shells by index ISH
       !   if IDQ = 2, then one orbital fixed on each side
       NS = NPEEL
    ENDIF

      
    !   Main computation
    !
    !     JA1, JA2 are the indices of interacting shells in JLIST
    !     IA1, IA2 are the indices of interacting shells in NW

    IF(idq .NE. 0) THEN
       IF(JA1 .GE. 1) IA1 = JLIST(JA1)
       IF(JA2 .GE. 1) IA2 = JLIST(JA2)
       IF(IA1 .GE. 1) KS1 = 2*IABS (INT(NAK(IA1)))
       IF(IA2 .GE. 1) KS2 = 2*IABS (INT(NAK(IA2)))
    END IF
    
    !   Check triangular condition for the active shells
    IF (ITRIG (KS1,KS2,KK).EQ.1) THEN
       
       IF (IDQ .EQ. 0) THEN
          !  Set tables of quantum numbers of non-interacting spectator shells
          !  IDQ = 0 Case         
          IF (NPEEL .NE. 0) THEN
             DO I = 1,NPEEL
                JLIS(I) = JLIST(I)
             END DO
             IF (NPEEL .GT. 1) THEN
                NPEELM = NPEEL-1
                DO I = 1,NPEELM
                   JC1S(I) = JJC1(I)
                   JC2S(I) = JJC2(I)
                END DO
             END IF
          END IF
          
          ISH = 1
          !   Loop over shells when IDQ = 0
          DO WHILE (ish .LE. nw)                  !   If ISH .GT. NW, then loop is over and return
             IF (stype(ISH,JA) .EQ. -1) CYCLE
             IF (stype(ISH,JA) .NE. 0) THEN
                
                !   Case one --- the ISH-th shell is in the core or in the peel and
                !   closed for both sides
                I = 1
                IF (NPEEL.GT.0) THEN
                   DO I = 1,NPEEL
                      IJ = JLIST(I)
                      IF (ISH .LT. IJ) THEN
                         IM = NPEEL-I+1
                         DO II = 1,IM
                            JLIST(NPEEL+2-II) = JLIST(NPEEL+1-II)
                            IF (NPEEL.EQ.II) EXIT
                            JJC1(NPEEL+1-II) = JJC1(NPEEL-II)
                            JJC2(NPEEL+1-II) = JJC2(NPEEL-II)
                         END DO
                         EXIT
                      END IF
                   END DO
                   
                   I = NPEEL+1
                   
                   IF (I .GE. 3) THEN
                      JJC1(I-1) = JJC1(I-2)
                      JJC2(I-1) = JJC2(I-2)
                   ELSE
                      I1 = JLIST(1)
                      JJC1(1) = JJQ1(3,I1)
                      JJC2(1) = JJQ2(3,I1)
                   ENDIF
                ELSE
                   JLIST(I) = ISH
                   JA1 = I
                   JA2 = I
                   NS = NPEEL+1
                END IF
             ELSE
                
                !   Case two --- the ISH-th shell is in the peel and open for either
                !   side
                NS = NPEEL
                DO  JW = 1,NPEEL
                   NX = ISH-JLIST(JW)
                   IF (NX.EQ.0) EXIT
                END DO
                JA1 = JW
                JA2 = JW
             END IF
             
             IF(JA .EQ. JB) THEN
                CALL ONESCALAR1(NS,JA,JB,JA1,JA2,TCOEFF)
             ELSE
                TCOEFF = ZERO
             END IF
             VSHELL(ISH) = TCOEFF
             
             !   Loop over all shells when IDQ = 0
             IF (NPEEL .NE. 0) THEN
                DO I = 1,NPEEL
                   JLIST(I) = JLIS(I)
                END DO
                IF (NPEEL .GT. 1) THEN
                   NPEELM = NPEEL-1
                   DO I = 1,NPEELM
                      JJC1(I) = JC1S(I)
                      JJC2(I) = JC2S(I)
                   END DO
                END IF
             END IF
             ish = ish + 1
          END DO
          
       ELSE
          
          !   IDQ = 2 Case
          CALL ONESCALAR2(JA,JB,JA1,JA2,TCOEFF)
          VSHELL(1) = TCOEFF 
          RETURN
       END IF
       
    ELSE
       IF (IDQ .EQ. 2) RETURN
    ENDIF
    
!    WRITE(*,'(30I3)') (jlist(k), k=1,npeel)      
    
    RETURN

  CONTAINS

      
    SUBROUTINE ONESCALAR1(NS,JJA,JJB,JA,JB,COEFF)
!*******************************************************************
!   --------------  SECTION METWO    SUBPROGRAM 03  -------------  *
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF ONE PARTICLE OPERATOR IN CASE :           N'1 = N1        *
!                                                  N'2 = N2        *
!                                                                  *
!      SUBROUTINE CALLED:                                          *
!                                                                  *
!   Written by  G. Gaigalas                                        *
!   The last modification made by G. Gaigalas       October  2017  *
!   Restructured by A. Senchuk                     September 2019  *     
!                                                                  *
!*******************************************************************
!
    IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
    INTEGER, INTENT(IN)       :: NS,JJA,JJB,JA,JB
    REAL(KIND=8), INTENT(OUT) :: COEFF
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
    INTEGER      :: IAT
    REAL(KIND=8) :: WJ,QM1,QM2,RECOUPL
!-----------------------------------------------
!
!     THE CASE 1111   + + - -
!
    COEFF=ZERO
    IF(JA == JB) THEN
       IF(JJA /= JJB) THEN
          CALL RECOONESCALAR(NS,JA,JA,JA,JA,0,IAT)
          IF(IAT == 0)RETURN
       END IF
       CALL PERKO2(JA,JA,JA,JA,1)
       QM1=HALF
       QM2=-HALF
       CALL WJ1(IK1,BK1,ID1,BD1,0,QM1,QM2,WJ)
       IF(DABS(WJ) > EPS) THEN
          RECOUPL=ONE/DSQRT(DBLE(IK1(6)+1))
          COEFF=WJ*RECOUPL*DSQRT(DBLE(ID1(3)+1))
          COEFF=-COEFF
       END IF
    END IF
    RETURN
  END SUBROUTINE ONESCALAR1
!
!####################################################################
!

  SUBROUTINE ONESCALAR2(JJA,JJB,JA,JB,COEFF)
!*******************************************************************
!   --------------  SECTION METWO    SUBPROGRAM 06  -------------  *
!                                                                  *
!     THIS PACKAGE EVALUATED THE CASES - 2111, 1211 ( + + - - ),   *
!     WHICH APPEARS IN CALCULATION MATRIX ELEMENTS BETWEEN         *
!     CONFIGURATIONS:                               N'1 = N1 - 1   *
!                                                   N'2 = N2 + 1   *
!                                                                  *
!   Written by  G. Gaigalas                                        *
!   The last modification made by G. Gaigalas       October  2017  *
!   Restructured by A. Senchuk                     September 2019  * 
!                                                                  *
!*******************************************************************
!
    IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
    INTEGER, INTENT(IN) :: JJA,JJB,JA,JB
    REAL(KIND=8), INTENT(OUT) :: COEFF
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
    INTEGER      :: IAT,JAA,JBB,NN,IB1,II,INN,IFAZ
    REAL(KIND=8) :: QM1,QM2,REC,WW
!-----------------------------------------------
    COEFF=ZERO
    IF(JA == JB) RETURN
    IF(JA < JB) THEN
       JAA=JA
       JBB=JB
    ELSE
       JAA=JB
       JBB=JA
    END IF
    CALL RECOONESCALAR(-1,JAA,JBB,JBB,JBB,1,IAT)
    IF(IAT == 0)RETURN
    QM1=HALF
    QM2=-HALF
    CALL PERKO2(JA,JB,JA,JA,2)
    IF(ID1(3) /= ID2(3)) RETURN
    CALL RECO2(JAA,JBB,ID2(3),0,IAT,REC)
    IF(IAT == 0)RETURN
    CALL GG12(IK1,IK2,BK1,BK2,ID1,ID2,BD1,BD2,QM1,QM2,WW)
    IF(DABS(WW) > EPS) THEN
       CALL RECO2(JAA,JBB,ID2(3),1,IAT,REC)
       COEFF=WW*REC*DSQRT(DBLE(ID1(3)+1))
       NN=0
       IB1=JBB-1
       DO II=JAA,IB1
          INN=JLIST(II)
          NN=NQ1(INN)+NN
       END DO
       IF((NN/2)*2 == NN) COEFF=-COEFF
!GG         IF(JA.GT.JB) COEFF=-COEFF
         COEFF=-COEFF
!
!     TRANSFORM FANO & RACAH PHASE CONVENTION
!     TO CONDON & SHORTLEY PHASE CONVENTION
!
         IFAZ=IK1(5)*IK1(4)+IK2(5)*IK2(4)-ID1(5)*ID1(4)-ID2(5)*ID2(4)
         IF((IFAZ/4)*4 /= IFAZ)COEFF=-COEFF
      END IF
      RETURN
    END SUBROUTINE ONESCALAR2

      
  END SUBROUTINE ONESCALAR
 

  SUBROUTINE RKCO (JA,JB,INCOR,ICOLBREI)
!***********************************************************************
!   Configurations JA, JB. Analyse the tables of quantum numbers set   *
!   in the COMMON  blocks M0 , M1, M2, M3  to determine all possible   *
!   sets of interacting  orbitals which give a non-vanishing Coulomb   *
!   matrix element,  and  initiates the calculation of coefficients.   *
!   The following conventions are in force: (1) labels 1, 2 refer to   *
!   left, right sides of matrix element respectively;   (2) pointers   *
!   JA1, JB1, JA2, JB2 point to the JLIST array of active  orbitals;   *
!   IA1, IB1, IA2, IB2 point to the complete list of orbitals.         *
!                                                                      *
!                                                                      *
!   Rewrite by G. Gaigalas                                             *
!   The last modification made by G. Gaigalas           October  2017  *
!                                                                      *
!***********************************************************************
  
    IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!      EXTERNAL CORD
    INTEGER, INTENT(IN) :: JA,JB,INCOR,ICOLBREI
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
    INTEGER :: I,IB1,IB2,IDQ,IDQG,II,I1,IJ,IJW,IM,IPCA,IT1,IT2, &
         JA1,JB1,JA2,JB2,J,JW,JT1,JT2,JT3,NDQ, &
         KLAST,KW,KWA,K1,KW1,KW2,NPEELM
!-----------------------------------------------
!
!   The Hamiltonian is an even scalar operator
    IF (ICHKQ2(JA,JB) .EQ. 0) RETURN

    ! initialize
    IDQ = 0; IDQG = 0
    JA1 = 0; JB1 = 0
    JA2 = 0; JB2 = 0

    CALL SETQNA (JA,JB)
!
!   1.0 Analyse peel shell interactions
!
!   1.1 Analyse electron distribution in peel. (The full procedure is
!       needed only if the number of peel orbitals NPEEL .GE. 2)
    IF (NW .LT. 1) THEN
       PRINT *, 'RKCO_GG: No subshells.'
       STOP
    ENDIF
    IF (NPEEL .EQ. 0) RETURN ! incor is always 0
!
!   Find differences in occupations, NDQ, for each peel orbital in
!   turn and use to set up labels of active orbitals maintaining the
!   convention JA1 .LE. JB1, JA2 .LE. JB2.
  IF(npeel .GE. 2) THEN  
     DO JW = 1,NPEEL
        J = JLIST(JW)
        NDQ = NQ1(J) - NQ2(J)
        IF (IABS (NDQ) .GT. 2) RETURN
        IF (NDQ .LT. 0) THEN
           IF (NDQ+1 .GT. 0) THEN
              CYCLE
           ELSE IF (NDQ+1 .EQ. 0) THEN
              IF (JA2 .GT. 0) THEN
                 JB2 = JW
              ELSE
                 JA2 = JW
              END IF
              IDQ = IDQ+1
              IDQG=1+IDQG
           ELSE IF (NDQ+1 .LT. 0) THEN
              JA2 = JW
              IDQ = IDQ+2
              IDQG=20+IDQG
           END IF
        ELSE
           IF (NDQ-1 .GT. 0) THEN
              JA1 = JW
              IDQ = IDQ+2
              IDQG=20+IDQG
              CYCLE
           ELSE IF (NDQ-1 .EQ. 0) THEN
              IF (JA1 .GT. 0) THEN
                 JB1 = JW
              ELSE
                 JA1 = JW
              END IF
              IDQ = IDQ+1
              IDQG=1+IDQG
              CYCLE
           ELSE
              CYCLE
           END IF
        END IF
     END DO
  END IF
    !
!   1.2 Calculate coefficients for all possible sets of active shells.
!
!   There are 4 cases, depending on the value of IDQ, the sum of the
!   absolute differences NDQ:
!
!   1.2.1 IDQ .GT. 4: matrix element null
    IF (IDQ .GT. 4) RETURN

!  1.2.2 IDQ .EQ. 4: matrix element uniquely defined
    IF (IDQ .EQ. 4) THEN
       
       IF (JB1 .EQ. 0) THEN
          JB1 = JA1
       END IF
       IF (JB2 .EQ. 0) THEN
          JB2 = JA2
       END IF
       CONTINUE
       IF(IDQG.NE.40) THEN
          IF(IDQG.NE.22) THEN
             !
!         TARP KONFIGURACIJU
!                        KAI      N'1=N1+-1
!                                 N'2=N2+-1
!                        KAI      N'3=N3+-1
             !                                 N'4=N4+-1
             CALL EL5(JA,JB,JA1,JB1,JA2,JB2,ICOLBREI)
          ELSE
             CALL EL4(JA,JB,JA1,JB1,JA2,JB2,ICOLBREI)
          END IF
       ELSE
          !
!       TARP KONFIGURACIJU
!                        KAI      N'1=N1+1
          !                                 N'2=N2-1
          CALL EL2(JA,JB,JA1,JA2,ICOLBREI)
       END IF
       RETURN
    ENDIF
    

!   1.2.3 IDQ .EQ. 2: One orbital fixed each side include all
!                     possible spectators.
!
!   Also IDQ .EQ. 0 for a matrix element off-diagonal in coupling
!   only. Must sum over all pairs of orbitals excluding core-core
!   terms
    IF (((IDQ.EQ.0 .AND. ja.NE.jb) .OR. IDQ.EQ.2) .AND. (npeel .GE. 2)) THEN
              
       IF (idq .EQ. 0) THEN
          klast = npeel
       ELSE
          klast = 1
       END IF

       DO KWA = 1,KLAST
          IF (IDQ .NE. 2) THEN
             JA1 = KWA
             JA2 = KWA
          END IF
          JT1 = JA1
          JT2 = JA2
          IT1 = JLIST(JA1)
          IT2 = JLIST(JA2)
          DO KW = KWA,NPEEL
             K1 = JLIST(KW)
             IF (NQ1(K1)*NQ2(K1) .EQ. 0) CYCLE
             JB1 = KW
             JB2 = KW
             JA1 = JT1
             JA2 = JT2
             !
             !   Interchange JA1 and JB1 and/or JA2 and JB2 if necessary
             IF (JA1-JB1 .GT. 0) THEN
                JT3 = JB1
                JB1 = JA1
                JA1 = JT3
             ELSE IF (JA1-JB1 .EQ. 0) THEN
                IB1 = JLIST(JB1)
                IF (NQ1(IB1) .LE. 1) CYCLE
             END IF
             IF (JA2-JB2 .GT. 0) THEN
                JT3 = JB2
                JB2 = JA2
                JA2 = JT3
             ELSE IF (JA2-JB2 .EQ. 0) THEN
                IB2 = JLIST(JB2)
                IF (NQ2(IB2) .LE. 1) CYCLE
             END IF
             IF(IDQ.NE.0) THEN
                !
                !     TARP KONFIGURACIJU
                !                        KAI      N'1=N1+1
                !                                 N'2=N2-1
                !                                 N'3=N3
                CALL EL3(JA,JB,JA1,JB1,JA2,JB2,ICOLBREI)
             ELSE
                !
                !     TARP TU PACIU KONFIGURACIJU
                CALL EL1(JA,JB,JA1,JB1,1,ICOLBREI)
             END IF
          END DO
          IF ((IDQ .EQ. 0) .AND. (NCORE .EQ. 0)) CYCLE
          IF ((NCORE .EQ. 0) .OR. (NAK(IT1) .NE. NAK(IT2))) RETURN
          !
          !   This section calculates the terms arising from active electrons
          !   which are in closed shells
          NPEELM = NPEEL-1
          DO I = 1,NPEEL
             JLIS(I) = JLIST(I)
          END DO
          DO I = 1,NPEELM
             JC1S(I) = JJC1(I)
             JC2S(I) = JJC2(I)
          END DO          
          DO KW = 1,NCORE
             IJW = KLIST(KW)
              DO I = 1,NPEEL
                IJ = JLIST(I)
                IF (IJW .LT. IJ) THEN
                   IM = NPEEL-I+1
                   DO II = 1,IM
                      JLIST(NPEEL+2-II) = JLIST(NPEEL+1-II)
                      IF (NPEEL .EQ. II) EXIT
                      JJC1(NPEEL+1-II) = JJC1(NPEEL-II)
                      JJC2(NPEEL+1-II) = JJC2(NPEEL-II)
                   END DO
                   EXIT
                END IF
             END DO
             
             I = NPEEL + 1
             IF (I .LT. 3) THEN
                I1 = JLIST(1)
                JJC1(1) = JJQ1(3,I1)
                JJC2(1) = JJQ2(3,I1)
             ELSE
                JJC1(I-1) = JJC1(I-2)
                JJC2(I-1) = JJC2(I-2)
             END IF
             JLIST(I) = IJW
             JA1 = JT1
             IF (JT1 .GE. I) JA1 = JA1+1
             JB1 = I
             JA2 = JT2
             IF (JT2 .GE. I) JA2 = JA2+1
             JB2 = I
             IF (JA1-JB1 .GT. 0) THEN
                JT3 = JB1
                JB1 = JA1
                JA1 = JT3
             END IF
             IF (JA2-JB2 .GT. 0) THEN
                JT3 = JB2
                JB2 = JA2
                JA2 = JT3
             END IF
             NPEEL = NPEEL+1
             IF(IDQ.NE.0) THEN
                IF(IDQG.NE.40) THEN
                   IF(IDQG.NE.2) THEN
                      WRITE(99,995)
                      RETURN
                   ELSE
                      !
                      !     TARP KONFIGURACIJU
                      !                        KAI      N'1=N1+1
                      !                                 N'2=N2-1
                      !                                 N'3=N3
                      CALL EL3(JA,JB,JA1,JB1,JA2,JB2,ICOLBREI)
                   END IF
                ELSE
                   WRITE(99,994)
                   RETURN
                END IF
             ELSE
                !
                !     TARP TU PACIU KONFIGURACIJU
                CALL EL1(JA,JB,JA1,JB1,1,ICOLBREI)
             END IF
             NPEEL = NPEEL-1
             NPEELM = NPEEL-1
             DO I = 1,NPEEL
                JLIST(I) = JLIS(I)
             END DO
             DO I = 1,NPEELM
                JJC1(I)  = JC1S(I)
                JJC2(I)  = JC2S(I)
             END DO
          END DO
       END DO
       RETURN
       
    ELSE

       IF (IDQ .NE. 2 .AND. IDQ .NE. 0) THEN
          !         3.0 Diagnostic print - NW .LT. 1
          WRITE (*,300)
          STOP
       END IF
       
       !   1.2.4 IDQ .EQ. 0 - diagonal case. Include all pairs with
       !         JA1 = JA2, JB1 = JB2.
       DO KW1 = 1,NPEEL
          K1 = JLIST(KW1)
          JB1 = KW1
          JB2 = KW1
          DO KW2 = 1,KW1
             JA1 = KW2
             IF (JA1 .EQ. JB1) THEN
                IF (NQ1(K1) .LE. 1) CYCLE
             END IF
             JA2 = JA1
             IF(JA.NE.JB) THEN
                IF(IDQG.NE.2) THEN
                   WRITE(99,996)
                   WRITE(99,*)"JA,JB,JA1,JB1",JA,JB,JA1,JB1
                   WRITE(99,*)"IDQG IDQ",IDQG,IDQ
                   RETURN
                ELSE
                   !
                   !                TARP KONFIGURACIJU
                   !                        KAI      N'1=N1+1
                   !                                 N'2=N2-1
                   !                                 N'3=N3
                   CALL EL3(JA,JB,JA1,JB1,JA2,JB2,ICOLBREI)
                END IF
             ELSE
                !
                !             TARP TU PACIU BUSENU
                CALL EL1(JA,JB,JA1,JB1,0,ICOLBREI)
             END IF
          END DO
       END DO
    ENDIF
    
    RETURN
    
300 FORMAT ('RKCO_GG: Error.')
994 FORMAT('   rie zymes 38?? atv N=N-N  !!!!!!')
995 FORMAT('   rie zymes 38 atv N=N-N  !!!!!!')
996 FORMAT('   rie zymes 45 atv N=N-N  !!!!!!')


  CONTAINS
    
  SUBROUTINE EL1(JJA,JJB,JA,JB,IIRE,ICOLBREI)
!*******************************************************************
!   --------------  SECTION METWO    SUBPROGRAM 03  -------------  *
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF TWO PARTICLE OPERATOR IN CASE :           N'1 = N1        *
!                                                  N'2 = N2        *
!                                                                  *
!   Written by  G. Gaigalas                                        *
!   The last modification made by G. Gaigalas       October  2017  *
!                                                                  *
!*******************************************************************
!
    IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
    INTEGER, INTENT(IN) :: JJA,JJB,JA,JB,IIRE,ICOLBREI 
!      DIMENSION CONE(7,20),S(12),IS(4),KAPS(4),KS(4)
!      DIMENSION PMGG(30),RAGG(30),J(2)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
    INTEGER :: II,IA,IB,IAT,IP1,IP2,IP3,IG1,IG2,IG3,IKK,I1,I2,I3,I4,& 
         IFAZ,J12,IBRD,IBRE,KRA,KRA1,L1,L2,MU,N,NU,ND1,ND2,   &
         NE1,NE2,NUP1
    INTEGER, DIMENSION(2) :: J
    INTEGER, DIMENSION(4) :: IS,KAPS,KS
    REAL(KIND=8)          :: QM1,QM2,QM3,QM4,AA,AB,A1,BB,SI,RECC,RAG
    REAL(KIND=8), DIMENSION(12)   :: S
    REAL(KIND=8), DIMENSION(30)   :: PMGG,RAGG
    REAL(KIND=8), DIMENSION(7,20) :: CONE 
!-----------------------------------------------
    IF(JA /= JB)GO TO 9
!
    IF(JJA /= JJB) RETURN
!
!     THE CASE 1111   + + - -
!
    IF(IIRE /= 0) THEN
       CALL RECO(JA,JA,JA,JA,0,IAT)
       IF(IAT /= 0)RETURN
    END IF
    CALL PERKO2(JA,JA,JA,JA,1)
    QM1=HALF
    QM2=HALF
    QM3=-HALF
    QM4=-HALF
    IA=JLIST(JA)
    J(1)=ID1(3)
    IP2=ITREXG(J(1),J(1),J(1),J(1),IKK)+1
    IF(IKK <= 0) RETURN
    IG2=IP2+IKK-1
    L1=(J(1)+1)/2
    IP1=IP2
    IG1=IG2
    IF (ICOLBREI == 2) THEN
       IS(1:4)=IA
       KAPS(1:4)=2*NAK(IS(1:4))
       KS(1:4)  =IABS(KAPS(1:4))
       CALL SNRC(IS,KAPS,KS,ND1,ND2,NE1,NE2,IBRD,IBRE)
       IF(IBRD <= 0)RETURN 
    END IF
    DO I2=IP1,IG1,2
       KRA=(I2-1)/2
       IF (ICOLBREI == 1) THEN
          CALL COULOM(L1,L1,L1,L1,ID1(5),ID1(5),ID1(5),ID1(5),KRA,A1)
          IF(DABS(A1) < EPS)CYCLE
          A1=-A1*HALF
       END IF
       AB=ZERO
       DO I3=IP2,IG2,2
          J12=(I3-1)/2
          IF(IXJTIK(J(1),J(1),KRA*2,J(1),J(1),J12*2) == 0)CYCLE
          CALL WW1(IK1,BK1,ID1,BD1,J12,QM1,QM2,QM3,QM4,AA)
          IF(DABS(AA) < EPS)CYCLE
          CALL SIXJ(J(1),J(1),KRA*2,J(1),J(1),J12*2,0,SI)
          AA=AA*SI*DSQRT(DBLE(I3))
          IFAZ=IK1(3)+J12+KRA
          IF((IFAZ/2)*2 /= IFAZ)AA=-AA
          AB=AB+AA
       END DO
!
!     RECOUPLING COEFFICIENTS
!
       IF (ICOLBREI == 1) THEN
          BB=AB*A1
          BB=BB/DSQRT(DBLE(IK1(6)+1))
          IF(DABS(BB) > EPS)CALL SPEAK(JJA,JJB,IA,IA,IA,IA,KRA,BB)
       ELSE IF (ICOLBREI == 2) THEN
          N=(KRA-ND1)/2+1
          IF(((KRA-ND1)/2)*2 == (KRA-ND1)) THEN
             CALL CXK(S,IS,KAPS,KRA,KRA,3,1)
             IF(DABS(S(1)) > EPS) THEN
                BB =-HALF*S(1)*AB/DSQRT(DBLE(IK1(6)+1))
                IF(DABS(BB) > EPS)CALL TALK(JJA,JJB,KRA,IA,IA,IA,IA,4,BB)
             END IF
          END IF
       END IF
    END DO
    RETURN
!  ............................................................
9   IF(NPEEL <= 1)RETURN
    IF(IIRE /= 0) THEN
       CALL RECO(JA,JB,JB,JB,1,IAT)
       IF(IAT == 0)RETURN
    END IF
    IA=JLIST(JA)
    IB=JLIST(JB)
    QM1=HALF
    QM2=-HALF
    QM3=HALF
    QM4=-HALF
    CALL PERKO2(JA,JB,JA,JA,2)
    J(1)=ID1(3)
    J(2)=ID2(3)
    L1=(J(1)+1)/2
    L2=(J(2)+1)/2
    IP1=ITREXG(J(1),J(1),J(2),J(2),IKK)+1
    IF(IKK <= 0)RETURN
    IG1=IP1+IKK-1
    IP3=IP1
    IG3=IG1
    DO I4=IP1,IG1,2
       KRA=(I4-1)/2
       KRA1=KRA+1
       IF(KRA1 > 30)GO TO 10
       RAGG(KRA1)=ZERO
       PMGG(KRA1)=ZERO
       CALL RECO2(JA,JB,KRA*2,0,IAT,RECC)
       IF(IAT == 0) CYCLE
       CALL GG1122(KRA,KRA,QM1,QM2,QM3,QM4,RAG)
       IF(DABS(RAG) < EPS) CYCLE
       RAGG(KRA1)=RAG
       CALL RECO2(JA,JB,KRA*2,1,IAT,RECC)
       PMGG(KRA1)=RECC
    END DO
! * * *                      * * *                      * * *
!     CASES 1212   + + - -        TRANSFORM TO  1122   + - + -
!           2121                                1122
!
    IF (ICOLBREI == 2) THEN
       IS(1)=IA
       IS(2)=IB
       IS(3)=IA
       IS(4)=IB
       KAPS(1:4)=2*NAK(IS(1:4))
       KS(1:4)=IABS(KAPS(1:4))
       CALL SNRC(IS,KAPS,KS,ND1,ND2,NE1,NE2,IBRD,IBRE)
       CONE = Zero
       IF(IBRD == 0 .AND. IBRE == 0)RETURN 
    END IF
    DO I1=IP1,IG1,2
       KRA=(I1-1)/2
       KRA1=KRA+1
       IF(KRA1 > 30)GO TO 10
       IF (ICOLBREI == 1) THEN
          CALL COULOM(L1,L2,L1,L2,ID1(5),ID2(5),ID1(5),ID2(5),KRA,AA)
          IF(DABS(AA) < EPS) CYCLE
          AA=AA*PMGG(KRA1)
          IF(DABS(AA) < EPS) CYCLE
          AA=AA*RAGG(KRA1)
          IF(DABS(AA) < EPS) CYCLE
          AA=AA/DSQRT(DBLE(I1))
          IF(DABS(AA) > EPS) CALL SPEAK(JJA,JJB,IA,IB,IA,IB,KRA,AA)
       ELSE IF (ICOLBREI == 2) THEN
          N=(KRA-ND1)/2+1
          IF(((KRA-ND1)/2)*2 == (KRA-ND1)) THEN
             CALL CXK(S,IS,KAPS,KRA,KRA,3,1)
             IF(DABS(S(1)) > EPS) THEN
                BB=S(1)*PMGG(KRA1)*RAGG(KRA1)/DSQRT(DBLE(I1))
                IF(DABS(BB) > EPS)CALL TALK(JJA,JJB,KRA,IA,IA,IB,IB,4,BB)
             END IF
          END IF
       END IF
    END DO
! * * *                      * * *                      * * *
!     CASES 1221   + + - -        TRANSFORM TO  1122   + - + -
!           2112                                1122
!
    IP2=ITREXG(J(1),J(2),J(1),J(2),IKK)+1
    IF(IKK <= 0) RETURN
    IG2=IP2+IKK-1
    DO I2=IP2,IG2,2
       KRA=(I2-1)/2
       IF(KRA > 30)GO TO 10
       IF (ICOLBREI == 1) THEN
          CALL COULOM(L1,L2,L2,L1,ID1(5),ID2(5),ID2(5),ID1(5),KRA,A1)
          IF(DABS(A1) < EPS) CYCLE
       END IF
       AB=ZERO
       DO I3=IP3,IG3,2
          J12=(I3-1)/2
          KRA1=J12+1
          IF(KRA1 > 30)GO TO 10
          AA=PMGG(KRA1)
          IF(DABS(AA) < EPS) CYCLE
          AA=AA*RAGG(KRA1)
          IF(DABS(AA) < EPS) CYCLE
          IF(IXJTIK(J(1),J(2),KRA*2,J(2),J(1),J12*2) == 0)CYCLE
          CALL SIXJ(J(1),J(2),KRA*2,J(2),J(1),J12*2,0,SI)
          AA=AA*SI*DSQRT(DBLE(I3))
          AB=AB+AA
       END DO
       IF (ICOLBREI == 1) THEN
          BB=A1*AB
          IF(DABS(BB) > EPS)CALL SPEAK(JJA,JJB,IA,IB,IB,IA,KRA,BB)
       ELSE IF (ICOLBREI == 2) THEN
          NU=KRA 
          IF(((NU-NE1)/2)*2 == (NU-NE1)) THEN
             IF((ITRIG(KS(1),KS(4),NU+NU+1) /= 0) .AND.  &
                  (ITRIG(KS(2),KS(3),NU+NU+1) /= 0)) THEN
                IF(NU > 0) THEN
                   N=(NU-NE1)/2+1
                   CALL CXK(S,IS,KAPS,NU,KRA,4,2)
                   DO MU = 1,3
                      CONE(MU,N)=CONE(MU,N)+AB*S(MU)
                   END DO
                END IF
             END IF
          END IF
          NU=KRA+1
          IF(((NU-NE1)/2)*2 == (NU-NE1)) THEN
             IF((ITRIG(KS(1),KS(4),NU+NU-1) /= 0) .AND.  &
                  (ITRIG(KS(2),KS(3),NU+NU-1) /= 0)) THEN
                IF(NU >= 0) THEN
                   N=(NU-NE1)/2+1
                   IF(N <= NE2) THEN
                      CALL CXK(S,IS,KAPS,NU,KRA,4,2)
                      DO MU = 1,3
                         CONE(MU,N)=CONE(MU,N)+AB*S(MU)
                      END DO
                   END IF
                END IF
             END IF
          END IF
          NU=KRA-1
          IF(((NU-NE1)/2)*2 == (NU-NE1)) THEN
             IF((ITRIG(KS(1),KS(4),NU+NU+3) /= 0) .AND.  &
                  (ITRIG(KS(2),KS(3),NU+NU+3) /= 0)) THEN
                IF(NU >= 0) THEN
                   N=(NU-NE1)/2+1
                   IF(N < NE2) THEN
                      CALL CXK(S,IS,KAPS,NU,KRA,4,2)
                      DO MU = 1,7
                         CONE(MU,N)=CONE(MU,N)+AB*S(MU)
                      END DO
                   END IF
                END IF
             END IF
          END IF
       END IF
    END DO
    IF (ICOLBREI == 2) THEN
       DO N = 1,NE2
          NU=NE1+2*(N-1)
          CALL TALK(JJA,JJB,NU,IB,IA,IB,IA,5,CONE(1,N))
          CALL TALK(JJA,JJB,NU,IA,IB,IB,IA,5,CONE(2,N))
          CALL TALK(JJA,JJB,NU,IA,IB,IA,IB,5,CONE(3,N))
          IF(N == NE2) CYCLE
          NUP1=NU+1
          CALL TALK(JJA,JJB,NUP1,IA,IB,IA,IB,6,CONE(4,N))
          CALL TALK(JJA,JJB,NUP1,IB,IA,IB,IA,6,CONE(5,N))
          CALL TALK(JJA,JJB,NUP1,IA,IB,IB,IA,6,CONE(6,N))
          CALL TALK(JJA,JJB,NUP1,IB,IA,IA,IB,6,CONE(7,N))
       END DO
    END IF
    RETURN
10  WRITE(99,100)
100 FORMAT(5X,'ERRO IN EL1  PMGG RAGG')
    STOP
  END SUBROUTINE EL1
!
!
  SUBROUTINE EL2(JJA,JJB,JA,JB,ICOLBREI)
!*******************************************************************  
!   --------------  SECTION METWO    SUBPROGRAM 04  -------------  *
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF TWO PARTICLE OPERATOR IN CASE :       N'1 = N1 - 2        *
!                                              N'2 = N2 + 2        *
!                                                                  *
!   Written by  G. Gaigalas                                        *
!   The last modification made by G. Gaigalas       October  2017  *
!                                                                  *
!*******************************************************************
!
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: JJA,JJB,JA,JB,ICOLBREI
!      DIMENSION J(2)
!      DIMENSION COND(12,20),S(12),IS(4),KAPS(4),KS(4)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IAT,IA,IB,IBRD,IBRE,IP1,IP2,IG1,IG2,IKK,II,I2,I3, &
                 IFAZ,IFAZ1,IFAZFRCS,J12,JAA,JBB, &
                 KRA,L1,L2,N,ND1,ND2,NE1,NE2,NUP1,NU,MU
      INTEGER, DIMENSION(2) :: J
      INTEGER, DIMENSION(4) :: IS,KAPS,KS
      REAL(KIND=8)          :: AA,A1,AB,BB,QM1,QM2,QM3,QM4,RECC,SI
      REAL(KIND=8), DIMENSION(12)    :: S
      REAL(KIND=8), DIMENSION(12,20) :: COND
!-----------------------------------------------
      IF(NPEEL <= 1)RETURN
      IF(JA > JB) THEN
        JAA=JB
        JBB=JA
      ELSE
        JAA=JA
        JBB=JB
      END IF
      CALL RECO(JAA,JBB,JBB,JBB,1,IAT)
      IF(IAT == 0)RETURN
      IA=JLIST(JA)
      IB=JLIST(JB)
      QM1=HALF
      QM2=HALF
      QM3=-HALF
      QM4=-HALF
      CALL PERKO2(JA,JB,JA,JA,2)
      J(1)=ID1(3)
      J(2)=ID2(3)
      L1=(J(1)+1)/2
      L2=(J(2)+1)/2
      IP1=ITREXG(J(1),J(1),J(2),J(2),IKK)+1
      IF(IKK <= 0)RETURN
      IG1=IP1+IKK-1
! * * *                      * * *                      * * *
!     THE CASE 1122   + + - -
!
      IP2=ITREXG(J(1),J(2),J(1),J(2),IKK)+1
      IF(IKK <= 0) RETURN
      IG2=IP2+IKK-1
      IF (ICOLBREI == 2) THEN
        IS(1:4)=IA
        KAPS(1:4)=2*NAK(IS(1:4))
        KS(1:4)=IABS(KAPS(1:4))
        CALL SNRC(IS,KAPS,KS,ND1,ND2,NE1,NE2,IBRD,IBRE)
        IF(IBRD <= 0)RETURN
        COND = Zero
      END IF
      DO I2=IP2,IG2,2
        KRA=(I2-1)/2
        IF (ICOLBREI == 1) THEN
          CALL COULOM(L1,L1,L2,L2,ID1(5),ID1(5),ID2(5),ID2(5),KRA,A1)
          IF(DABS(A1) < EPS) CYCLE
          A1=-HALF*A1
        END IF
        AB=ZERO
          DO I3=IP1,IG1,2
            J12=(I3-1)/2
            CALL RECO2(JAA,JBB,J12*2,0,IAT,RECC)
            IF(IAT /= 0) THEN
              IF(IXJTIK(J(1),J(2),KRA*2,J(2),J(1),J12*2) /= 0) THEN
                CALL GG1122(J12,J12,QM1,QM2,QM3,QM4,AA)
                IF(DABS(AA) > EPS) THEN
                  CALL RECO2(JAA,JBB,J12*2,1,IAT,RECC)
                  AA=AA*RECC
                  CALL SIXJ(J(1),J(2),KRA*2,J(2),J(1),J12*2,0,SI)
                  AA=AA*SI*DSQRT(DBLE(I3))
                  IFAZ=IK1(3)+IK2(3)+KRA*2+J12*2
                  IF((IFAZ/4)*4 /= IFAZ)AA=-AA
                  AB=AB+AA
                END IF
              END IF
            END IF
          END DO
!
!       TRANSFORM FANO & RACAH PHASE CONVENTION
!       TO CONDON & SHORTLEY PHASE CONVENTION
!
        IFAZFRCS=1
        IFAZ1=IK1(5)*IK1(4)+IK2(5)*IK2(4)-ID1(5)*ID1(4)-ID2(5)*ID2(4)
        IF((IFAZ1/4)*4 /= IFAZ1)IFAZFRCS=-IFAZFRCS
!
        IF (ICOLBREI == 1) THEN
           BB=A1*AB*DBLE(IFAZFRCS)
          IF(DABS(BB) > EPS)CALL SPEAK(JJA,JJB,IA,IA,IB,IB,KRA,BB)
        ELSE IF (ICOLBREI == 2) THEN
          NU=KRA
          IF(((NU-ND1)/2)*2 == (NU-ND1)) THEN
            IF((ITRIG(KS(1),KS(3),NU+NU+1) /= 0) .AND.   &
               (ITRIG(KS(2),KS(4),NU+NU+1) /= 0)) THEN
              N=(NU-ND1)/2+1
              IF(NU > 0) THEN
                CALL CXK(S,IS,KAPS,NU,KRA,1,1)
                DO MU = 1,4
                  COND(MU,N)=COND(MU,N)-HALF*AB*S(MU)
                END DO
              END IF
            END IF
          END IF
          NU=KRA+1
          IF(((NU-ND1)/2)*2 == (NU-ND1)) THEN
            IF((ITRIG(KS(1),KS(3),NU+NU-1) /= 0) .AND.  &
               (ITRIG(KS(2),KS(4),NU+NU-1) /= 0)) THEN
              N=(NU-ND1)/2+1
              IF(N <= ND2) THEN
                CALL CXK(S,IS,KAPS,NU,KRA,1,1)
                DO MU = 1,4
                  COND(MU,N)=COND(MU,N)-HALF*AB*S(MU)
                END DO
              END IF
            END IF
          END IF
          NU=KRA-1
          IF(((NU-ND1)/2)*2 == (NU-ND1)) THEN
            IF((ITRIG(KS(1),KS(3),NU+NU+3) /= 0) .AND.   &
               (ITRIG(KS(2),KS(4),NU+NU+3) /= 0)) THEN
              IF(NU >= 0) THEN
                N=(NU-ND1)/2+1
                IF(N < ND2) THEN
                  CALL CXK(S,IS,KAPS,NU,KRA,1,1)
                  DO MU = 1,12
                    COND(MU,N)=COND(MU,N)-HALF*AB*S(MU)
                  END DO
                END IF
              END IF
            END IF
          END IF
        END IF
      END DO
      IF (ICOLBREI .EQ. 2) THEN
        DO N = 1,ND2
          NU=ND1+2*(N-1)
          CALL TALK(JJA,JJB,NU,IS(1),IS(3),IS(2),IS(4),1,COND(1,N))
          CALL TALK(JJA,JJB,NU,IS(3),IS(1),IS(4),IS(2),1,COND(2,N))
          CALL TALK(JJA,JJB,NU,IS(1),IS(3),IS(4),IS(2),1,COND(3,N))
          CALL TALK(JJA,JJB,NU,IS(3),IS(1),IS(2),IS(4),1,COND(4,N))
          IF(N == ND2) CYCLE
          NUP1=NU+1
          CALL TALK(JJA,JJB,NUP1,IS(1),IS(3),IS(2),IS(4),2,COND(5,N))
          CALL TALK(JJA,JJB,NUP1,IS(2),IS(4),IS(1),IS(3),2,COND(6,N))
          CALL TALK(JJA,JJB,NUP1,IS(3),IS(1),IS(4),IS(2),2,COND(7,N))
          CALL TALK(JJA,JJB,NUP1,IS(4),IS(2),IS(3),IS(1),2,COND(8,N))
          CALL TALK(JJA,JJB,NUP1,IS(1),IS(3),IS(4),IS(2),2,COND(9,N))
          CALL TALK(JJA,JJB,NUP1,IS(4),IS(2),IS(1),IS(3),2,COND(10,N))
          CALL TALK(JJA,JJB,NUP1,IS(3),IS(1),IS(2),IS(4),2,COND(11,N))
          CALL TALK(JJA,JJB,NUP1,IS(2),IS(4),IS(3),IS(1),2,COND(12,N))
        END DO
      END IF
      RETURN
      END SUBROUTINE EL2
      !
      !
!*******************************************************************
!                                                                  *
      SUBROUTINE EL3(JJA,JJB,JA,JB,JC,JD,ICOLBREI)
!                                                                  *
!   --------------  SECTION METWO    SUBPROGRAM 05  -------------  *
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF TWO PARTICLE OPERATOR IN CASE :       N'1 = N1 - 1        *
!                                              N'2 = N2 + 1        *
!                                                                  *
!   Written by  G. Gaigalas                                        *
!   The last modification made by G. Gaigalas       October  2017  *
!                                                                  *
!*******************************************************************
!
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: JJA,JJB,JA,JB,JC,JD,ICOLBREI
!-----------------------------------------------
      IF(NPEEL <= 1)RETURN
      IF(JB == JD) THEN
        IF(JA == JB.OR.JC == JB) THEN
          IF(JA == JC)GO TO 10
          IF(JC /= JB) THEN
            CALL EL32(JJA,JJB,JC,JA,JA,JB,JC,JD,ICOLBREI)
          ELSE
            CALL EL31(JJA,JJB,JC,JA,JA,JB,JC,JD,ICOLBREI)
          END IF
        ELSE
          CALL EL33(JJA,JJB,JC,JA,JB,1,JA,JB,JC,JD,ICOLBREI)
        END IF
        RETURN
      ELSE IF(JA == JC) THEN
        IF(JB == JA.OR.JD == JA) THEN
          IF(JB == JD)GO TO 10
          IF(JD /= JA) THEN
            CALL EL32(JJA,JJB,JD,JB,JA,JB,JC,JD,ICOLBREI)
          ELSE
            CALL EL31(JJA,JJB,JD,JB,JA,JB,JC,JD,ICOLBREI)
          END IF
        ELSE
          CALL EL33(JJA,JJB,JD,JB,JA,1,JA,JB,JC,JD,ICOLBREI)
        END IF
        RETURN
      ELSE IF(JA == JD) THEN
        IF(JB == JA.OR.JC == JA) THEN
          IF(JB == JC)GO TO 10
          IF(JC /= JD) THEN
             CALL EL32(JJA,JJB,JC,JB,JA,JB,JC,JD,ICOLBREI)
          ELSE
             CALL EL31(JJA,JJB,JC,JB,JA,JB,JC,JD,ICOLBREI)
          END IF
        ELSE
          CALL EL33(JJA,JJB,JC,JB,JA,2,JA,JB,JD,JC,ICOLBREI)
        END IF
        RETURN
      ELSE IF(JB == JC) THEN
        IF(JA == JB.OR.JD == JB) THEN
          IF(JA == JD)GO TO 10
          IF(JD /= JB) THEN
            CALL EL32(JJA,JJB,JD,JA,JA,JB,JC,JD,ICOLBREI)
          ELSE
            CALL EL31(JJA,JJB,JD,JA,JA,JB,JC,JD,ICOLBREI)
          END IF
        ELSE
          CALL EL33(JJA,JJB,JD,JA,JB,2,JA,JB,JD,JC,ICOLBREI)
        END IF
        RETURN
      END IF
   10 WRITE(99,100)
  100 FORMAT(5X,'ERRO IN EL3  PMGG RAGG')
      STOP
      END SUBROUTINE EL3
      !
      !
!*******************************************************************
!                                                                  *
      SUBROUTINE EL31(JJJA,JJJB,JA,JB,JJA,JJB,JJC,JJD,ICOLBREI)
!                                                                  *
!   --------------  SECTION METWO    SUBPROGRAM 06  -------------  *
!                                                                  *
!     THIS PACKAGE EVALUATED THE CASES - 2111, 1211 ( + + - - ),   *
!     WHICH APPEARS IN CALCULATION MATRIX ELEMENTS BETWEEN         *
!     CONFIGURATIONS:                               N'1 = N1 - 1   *
!                                                   N'2 = N2 + 1   *
!                                                                  *
!   Written by  G. Gaigalas                                        *
!   The last modification made by G. Gaigalas       October  2017  *
!                                                                  *
!*******************************************************************
!
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: JJJA,JJJB,JA,JB,JJA,JJB,JJC,JJD,ICOLBREI
!      DIMENSION J(2)
!      DIMENSION S(12),IS(4),KAPS(4),KS(4)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER ::  IA,IB,II,I2,I3,IAT,IIA,IIB,IIC,IID,IKK,IP1,IG1, &
                  IBRD,IBRE,IFAZ,IFAZFRCS,INN,JAA,JBB,JB1,J12,L1, &
                  L2,KRA,ND1,ND2,NE1,NE2,N,NN
      INTEGER, DIMENSION(2) :: J
      INTEGER, DIMENSION(4) :: IS,KAPS,KS
      REAL(KIND=8)          :: AA,AB,A1,BB,SI,QM1,QM2,QM3,QM4,RECC
      REAL(KIND=8), DIMENSION(12) :: S
!-----------------------------------------------
      IF(NPEEL <= 1)RETURN
      IIA=JLIST(JJA)
      IIB=JLIST(JJB)
      IIC=JLIST(JJC)
      IID=JLIST(JJD)
      IF(JA > JB) THEN
        JAA=JB
        JBB=JA
      ELSE
        JAA=JA
        JBB=JB
      END IF
      CALL RECO(JAA,JBB,JBB,JBB,1,IAT)
      IF(IAT == 0)RETURN
      IA=JLIST(JA)
      IB=JLIST(JB)
      QM1=HALF
      QM2=HALF
      QM3=-HALF
      QM4=-HALF
      CALL PERKO2(JA,JB,JA,JA,2)
      J(1)=ID1(3)
      J(2)=ID2(3)
      L1=(J(1)+1)/2
      L2=(J(2)+1)/2
      CALL RECO2(JAA,JBB,J(2),0,IAT,RECC)
      IF(IAT == 0)RETURN
      IP1=ITREXG(J(1),J(1),J(1),J(2),IKK)+1
      IF(IKK <= 0)RETURN
      IG1=IP1+IKK-1
      CALL RECO2(JAA,JBB,J(2),1,IAT,RECC)
      IF (ICOLBREI == 2) THEN
        IS(1)=IIA
        IS(2)=IIB
        IS(3)=IIC
        IS(4)=IID
        KAPS(1:4)=2*NAK(IS(1:4))
        KS(1:4)=IABS(KAPS(1:4))
        CALL SNRC(IS,KAPS,KS,ND1,ND2,NE1,NE2,IBRD,IBRE)
        IF(IBRD <= 0)RETURN
      END IF
! * * *                      * * *                      * * *
!     CASES 2111   + + - -        TRANSFORM TO  1112   + - - +
!           1211                                1112
!
      DO I2=IP1,IG1,2
        KRA=(I2-1)/2
        IF (ICOLBREI == 1) THEN
          CALL COULOM(L2,L1,L1,L1,ID2(5),ID1(5),ID1(5),ID1(5),KRA,A1)
          IF(DABS(A1) < EPS) CYCLE
          A1=-A1
        END IF
        AB=ZERO
        DO I3=IP1,IG1,2
          J12=(I3-1)/2
          IFAZ=J(2)-J12+1
          IF((IFAZ/2)*2 /= IFAZ) CYCLE
          IF(IXJTIK(J(2),J(1),KRA*2,J(1),J(1),J12*2) == 0)CYCLE
          CALL GG1222(IK2,IK1,BK2,BK1,ID2,ID1,BD2,BD1,J12,  &
                      QM1,QM2,QM3,QM4,AA)
          IF(DABS(AA) < EPS) CYCLE
          CALL SIXJ(J(2),J(1),KRA*2,J(1),J(1),J12*2,0,SI)
          AA=AA*SI*DSQRT(DBLE(I3))
          IFAZ=2*J(1)+KRA*2+J12*2
          IF((IFAZ/4)*4 /= IFAZ)AA=-AA
          AB=AB+AA
        END DO
        AB=AB*RECC
        IF(DABS(AB) < EPS) CYCLE
!
!       TRANSFORM FANO & RACAH PHASE CONVENTION
!       TO CONDON & SHORTLEY PHASE CONVENTION
!
        IFAZFRCS = 1
        IFAZ=IK1(5)*IK1(4)+IK2(5)*IK2(4)-ID1(5)*ID1(4)-ID2(5)*ID2(4)
        IF((IFAZ/4)*4 /= IFAZ)IFAZFRCS=-IFAZFRCS
!
        NN=0
        JB1=JBB-1
        DO II=JAA,JB1
          INN=JLIST(II)
          NN=NQ1(INN)+NN
        END DO
        IF((NN/2)*2 == NN)AB=-AB
        IF (ICOLBREI == 1) THEN
           BB=A1*AB*DBLE(IFAZFRCS)
           CALL SPEAK(JJJA,JJJB,IIA,IIB,IIC,IID,KRA,BB)
        ELSE IF (ICOLBREI == 2) THEN
          N=(KRA-ND1)/2+1
          IF(((KRA-ND1)/2)*2 == (KRA-ND1)) THEN
            CALL CXK(S,IS,KAPS,KRA,KRA,2,1)
            IF(DABS(S(1)) > EPS) THEN
              BB=-S(1)*AB
              IF(DABS(BB) > EPS)                                &
              CALL TALK(JJJA,JJJB,KRA,IS(1),IS(3),IS(2),IS(4),3,BB)
            END IF
          END IF
        END IF
      END DO
      RETURN
      END SUBROUTINE EL31
      !
      !
!*******************************************************************
!                                                                  *
      SUBROUTINE EL32(JJJA,JJJB,JA,JB,JJA,JJB,JJC,JJD,ICOLBREI)
!                                                                  *
!   --------------  SECTION METWO    SUBPROGRAM 07  -------------  *
!                                                                  *
!     THIS PACKAGE EVALUATED THE CASES - 2221, 2212 ( + + - - ),   *
!     WHICH APPEARS IN CALCULATION MATRIX ELEMENTS BETWEEN         *
!     CONFIGURATIONS:                               N'1 = N1 - 1   *
!                                                   N'2 = N2 + 1   *
!                                                                  *
!   Written by  G. Gaigalas                                        *
!   The last modification made by G. Gaigalas       October  2017  *
!                                                                  *
!*******************************************************************
!
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: JJJA,JJJB,JA,JB,JJA,JJB,JJC,JJD,ICOLBREI
!      DIMENSION J(2)
!      DIMENSION S(12),IS(4),KAPS(4),KS(4)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IIA,IIB,IIC,IID,IATIKK,IP1,IG1,II,I2,I3,IBRD,IBRE,    &
                 IA,IB,IAT,IFAZ,IFAZFRCS,IKK,INN,JB1,JAA,JBB,J12,KRA,  &
                 L1,L2,N,NN,ND1,ND2,NE1,NE2
      INTEGER, DIMENSION(2) :: J
      INTEGER, DIMENSION(4) :: IS,KAPS,KS
      REAL(KIND=8)          :: AA,AB,A1,BB,QM1,QM2,QM3,QM4,SI,RECC
      REAL(KIND=8), DIMENSION(12) :: S
!-----------------------------------------------
      IF(NPEEL <= 1)RETURN
      IIA=JLIST(JJA)
      IIB=JLIST(JJB)
      IIC=JLIST(JJC)
      IID=JLIST(JJD)
      IF(JA > JB) THEN
        JAA=JB
        JBB=JA
      ELSE
        JAA=JA
        JBB=JB
      END IF
      CALL RECO(JAA,JBB,JBB,JBB,1,IAT)
      IF(IAT == 0)RETURN
      IA=JLIST(JA)
      IB=JLIST(JB)
      QM1=HALF
      QM2=HALF
      QM3=-HALF
      QM4=-HALF
      CALL PERKO2(JA,JB,JA,JA,2)
      J(1)=ID1(3)
      J(2)=ID2(3)
      L1=(J(1)+1)/2
      L2=(J(2)+1)/2
      CALL RECO2(JAA,JBB,J(1),0,IAT,RECC)
      IF(IAT == 0)RETURN
      IP1=ITREXG(J(2),J(2),J(1),J(2),IKK)+1
      IF(IKK <= 0)RETURN
      IG1=IP1+IKK-1
      IF (ICOLBREI == 2) THEN
        IS(1)=IIA
        IS(2)=IIB
        IS(3)=IIC
        IS(4)=IID
        KAPS(1:4)=2*NAK(IS(1:4))
        KS(1:4)=IABS(KAPS(1:4))
        CALL SNRC(IS,KAPS,KS,ND1,ND2,NE1,NE2,IBRD,IBRE)
        IF(IBRD <= 0)RETURN
      END IF
      CALL RECO2(JAA,JBB,J(1),1,IAT,RECC)
! * * *                      * * *                      * * *
!     CASES 2221   + + - -        TRANSFORM TO  1222   - + + -
!           2212                                1222
!
      DO I2=IP1,IG1,2
        KRA=(I2-1)/2
        IF (ICOLBREI == 1) THEN
          CALL COULOM(L2,L2,L2,L1,ID2(5),ID2(5),ID2(5),ID1(5),KRA,A1)
          IF(DABS(A1) < EPS) CYCLE
          A1=-A1
        END IF
!
        AB=ZERO
        DO I3=IP1,IG1,2
          J12=(I3-1)/2
          IFAZ=J(1)-J12+1
          IF((IFAZ/2)*2 /= IFAZ) CYCLE
          IF(IXJTIK(J(2),J(2),KRA*2,J(1),J(2),J12*2) == 0) CYCLE
          CALL GG1112(IK2,IK1,BK2,BK1,ID2,ID1,BD2,     &
                      BD1,J12,QM1,QM2,QM3,QM4,AA)
          IF(DABS(AA) < EPS) CYCLE
          CALL SIXJ(J(2),J(2),KRA*2,J(1),J(2),J12*2,0,SI)
          AA=AA*SI*DSQRT(DBLE(I3))
          IFAZ=2*J(2)+KRA*2+J12*2
          IF((IFAZ/4)*4 /= IFAZ)AA=-AA
          AB=AB+AA
        END DO
        AB=AB*RECC
        IF(DABS(AB) < EPS) CYCLE
!
!       TRANSFORM FANO & RACAH PHASE CONVENTION
!       TO CONDON & SHORTLEY PHASE CONVENTION
!
        IFAZFRCS=1
        IFAZ=IK1(5)*IK1(4)+IK2(5)*IK2(4)-ID1(5)*ID1(4)-ID2(5)*ID2(4)
        IF((IFAZ/4)*4 /= IFAZ)IFAZFRCS=-IFAZFRCS
!
        NN=0
        JB1=JBB-1
        DO II=JAA,JB1
          INN=JLIST(II)
          NN=NQ1(INN)+NN
        END DO
        IF((NN/2)*2 == NN)AB=-AB
        IF (ICOLBREI == 1) THEN
          BB=AB*A1*DBLE(IFAZFRCS)
          CALL SPEAK(JJJA,JJJB,IIA,IIB,IIC,IID,KRA,BB)
        ELSE IF (ICOLBREI == 2) THEN
          N=(KRA-ND1)/2+1
          IF(((KRA-ND1)/2)*2 == (KRA-ND1)) THEN
            CALL CXK(S,IS,KAPS,KRA,KRA,IBRD,1)
            IF(DABS(S(1)) > EPS) THEN
              BB=-S(1)*AB
              IF(DABS(BB) > EPS)                              &
              CALL TALK(JJJA,JJJB,KRA,IS(1),IS(3),IS(2),IS(4),3,BB)
            END IF
          END IF
        END IF
      END DO
      RETURN
      END SUBROUTINE EL32
      !
      !
!*******************************************************************
!                                                                  *
      SUBROUTINE EL33(JJJA,JJJB,JA,JB,JC,IREZ,JJA,JJB,JJC,JJD   &
                                                        ,ICOLBREI)
!                                                                  *
!   --------------  SECTION METWO    SUBPROGRAM 08  -------------  *
!                                                                  *
!     THIS PACKAGE EVALUATED THE CASES - 2313, 3231, 3213, 2331    *
!                                                   ( + + - - ),   *
!     WHICH APPEARS IN CALCULATION MATRIX ELEMENTS BETWEEN         *
!     CONFIGURATIONS:                               N'1 = N1 - 1   *
!                                                   N'2 = N2 + 1   *
!                                                                  *
!   Written by  G. Gaigalas                                        *
!   The last modification made by G. Gaigalas       October  2017  *
!                                                                  *
!*******************************************************************
!
        IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: JJJA,JJJB,JA,JB,JC,IREZ,JJA,JJB,JJC,JJD &
                                                        ,ICOLBREI
!      DIMENSION PMGG(30),RAGG(30),J(3)
!      DIMENSION CONE(12,20),S(12),IS(4),KAPS(4),KS(4)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IA,IB,IC,IAT,II,IIA,IIB,IIC,IID,IP1,IP2,IG1,IG2,   &
                 IKK,INN,I1,I2,I3,I4,IBRD,IBRE,IFAZ,IFAZP,IFAZFRCS, &
                 JAA,JBB,JCC,JB1,J12,KRA,KRA1,L1,L2,L3,ND1,ND2,NE1, &
                 NE2,N,NN,NU,NUP1,MU
      INTEGER, DIMENSION(3) :: J
      INTEGER, DIMENSION(4) :: IS,KAPS,KS
      REAL(KIND=8)          :: A1,AA,AB,BB,QM1,QM2,QM3,QM4,RAG,RECC,SI
      REAL(KIND=8), DIMENSION(12) :: S
      REAL(KIND=8), DIMENSION(30) :: PMGG,RAGG
      REAL(KIND=8), DIMENSION(12,20) :: CONE
!-----------------------------------------------
      IF(NPEEL <= 1)RETURN
      CALL EILE(JA,JB,JC,JAA,JBB,JCC)
      CALL RECO(JAA,JCC,JBB,JBB,2,IAT)
      IF(IAT == 0)RETURN
      IA=JLIST(JA)
      IB=JLIST(JB)
      IC=JLIST(JC)
      IIA=JLIST(JJA)
      IIB=JLIST(JJB)
      IIC=JLIST(JJC)
      IID=JLIST(JJD)
      QM1=HALF
      QM2=-HALF
      QM3=HALF
      QM4=-HALF
      CALL PERKO2(JA,JB,JC,JA,3)
      J(1)=ID1(3)
      J(2)=ID2(3)
      J(3)=ID3(3)
      L1=(J(1)+1)/2
      L2=(J(2)+1)/2
      L3=(J(3)+1)/2
      IP1=ITREXG(J(2),J(1),J(3),J(3),IKK)+1
      IF(IKK <= 0)RETURN
      IG1=IP1+IKK-1
      IF (ICOLBREI == 2) THEN
        IS(1)=IIA
        IS(2)=IIB
        IS(3)=IIC
        IS(4)=IID
        KAPS(1:4)=2*NAK(IS(1:4))
        KS(1:4)=IABS(KAPS(1:4))
        CALL SNRC(IS,KAPS,KS,ND1,ND2,NE1,NE2,IBRD,IBRE)
        IF(IBRD <= 0 .AND. IBRE <= 0)RETURN
        CONE = Zero
      END IF
      DO I4=IP1,IG1,2
        KRA=(I4-1)/2
        KRA1=KRA+1
        IF(KRA1 > 30)GO TO 10
        RAGG(KRA1)=ZERO
        PMGG(KRA1)=ZERO
        CALL REC3(JB,JA,JC,J(2),J(1),KRA*2,0,IAT,RECC)
        IF(IAT == 0) CYCLE
        CALL GG1233(IK2,IK1,IK3,BK2,BK1,BK3,ID2,ID1,ID3,BD2,  &
                    BD1,BD3,KRA,QM1,QM2,QM3,QM4,RAG)
        IF(DABS(RAG) < EPS) CYCLE
        RAGG(KRA1)=RAG
        CALL REC3(JB,JA,JC,J(2),J(1),KRA*2,1,IAT,RECC)
        PMGG(KRA1)=RECC
      END DO
      IFAZP=JFAZE(JB,JA,JC,JC)
!
!     TRANSFORM FANO & RACAH PHASE CONVENTION
!     TO CONDON & SHORTLEY PHASE CONVENTION
!
      IFAZFRCS=1
      IFAZ=IK1(5)*IK1(4)+IK2(5)*IK2(4)-ID1(5)*ID1(4)-ID2(5)*ID2(4)
      IF((IFAZ/4)*4 /= IFAZ)IFAZFRCS=-IFAZFRCS
!
      IF(JA > JB) THEN
        JAA=JB
        JBB=JA
      ELSE
        JAA=JA
        JBB=JB
      END IF
      NN=0
      JB1=JBB-1
      DO II=JAA,JB1
        INN=JLIST(II)
        NN=NQ1(INN)+NN
      END DO
      IF((NN/2)*2 == NN)IFAZP=-IFAZP
      IF(IREZ == 2)GO TO 5
! * * *                      * * *                      * * *
!     CASES 2313   + + - -        TRANSFORM TO  2133   + - + -
!           3231                                2133
!
    6 CONTINUE
      DO I1=IP1,IG1,2
        KRA=(I1-1)/2
        KRA1=KRA+1
        IF(KRA1 > 30)GO TO 10
        IF (ICOLBREI == 1) THEN
          CALL COULOM(L2,L3,L1,L3,ID2(5),ID3(5),ID1(5),ID3(5),KRA,A1)
          IF(DABS(A1) < EPS) CYCLE
        END IF
        AA=PMGG(KRA1)
        IF(DABS(AA) < EPS) CYCLE
        AA=AA*RAGG(KRA1)
        IF(DABS(AA) < EPS) CYCLE
        AA=AA/DSQRT(DBLE(I1))
        AA=AA*DBLE(IFAZP)
        IF (ICOLBREI == 1) THEN
          BB=A1*AA*DBLE(IFAZFRCS)
          CALL SPEAK(JJJA,JJJB,IIA,IIB,IIC,IID,KRA,BB)
        ELSE IF (ICOLBREI == 2) THEN
          N=(KRA-ND1)/2+1
          CALL CXK(S,IS,KAPS,KRA,KRA,2,1)
          IF(DABS(S(1)) > EPS) THEN
            BB=S(1)*AA
            IF(DABS(BB) > EPS)                    &
          CALL TALK(JJJA,JJJB,KRA,IS(1),IS(3),IS(2),IS(4),3,BB)
          END IF
        END IF
      END DO
      IF(IREZ == 2) GO TO 7
! * * *                      * * *                      * * *
!     CASES 3213   + + - -        TRANSFORM TO  2133   + - + -
!           2331                                2133
!
    5 IP2=ITREXG(J(3),J(1),J(2),J(3),IKK)+1
      IF(IKK <= 0) RETURN
      IG2=IP2+IKK-1
      DO I2=IP2,IG2,2
        KRA=(I2-1)/2
        IF(KRA > 30)GO TO 10
        IF (ICOLBREI == 1) THEN
          CALL COULOM(L3,L2,L1,L3,ID3(5),ID2(5),ID1(5),ID3(5),KRA,A1)
          IF(DABS(A1) < EPS) CYCLE
        END IF
        AB=ZERO
        DO I3=IP1,IG1,2
          J12=(I3-1)/2
          KRA1=J12+1
          IF(KRA1 > 30)GO TO 10
          AA=PMGG(KRA1)
          IF(DABS(AA) < EPS) CYCLE
          AA=AA*RAGG(KRA1)
          IF(DABS(AA) < EPS) CYCLE
          IF(IXJTIK(J(1),J(3),KRA*2,J(3),J(2),J12*2) == 0) CYCLE
          CALL SIXJ(J(1),J(3),KRA*2,J(3),J(2),J12*2,0,SI)
          AA=AA*SI*DSQRT(DBLE(I3))
          AB=AB+AA
        END DO
        IF(DABS(AB) < EPS) CYCLE
        AB=AB*DBLE(IFAZP)
        IF (ICOLBREI == 1) THEN
          BB=A1*AB*DBLE(IFAZFRCS)
          CALL SPEAK(JJJA,JJJB,IIA,IIB,IID,IIC,KRA,BB)
        ELSE IF (ICOLBREI == 2) THEN
          NU=KRA
          IF(((NU-NE1)/2)*2 == (NU-NE1)) THEN
            IF((ITRIG(KS(1),KS(4),NU+NU+1) /= 0) .AND.  &
               (ITRIG(KS(2),KS(3),NU+NU+1) /= 0)) THEN
              IF(NU > 0) THEN
                N=(NU-NE1)/2+1
                CALL CXK(S,IS,KAPS,NU,KRA,1,2)
                DO MU = 1,4
                  CONE(MU,N)=CONE(MU,N)+AB*S(MU)
                END DO
              END IF
            END IF
          END IF
          NU=KRA+1
          IF(((NU-NE1)/2)*2 == (NU-NE1)) THEN
            IF((ITRIG(KS(1),KS(4),NU+NU-1) /= 0) .AND.  &
               (ITRIG(KS(2),KS(3),NU+NU-1) /= 0)) THEN
              IF(NU >= 0) THEN
                N=(NU-NE1)/2+1
                IF(N <= NE2) THEN
                  CALL CXK(S,IS,KAPS,NU,KRA,1,2)
                  DO MU = 1,4
                    CONE(MU,N)=CONE(MU,N)+AB*S(MU)
                  END DO
                END IF
              END IF
            END IF
          END IF
          NU=KRA-1
          IF(((NU-NE1)/2)*2 == (NU-NE1)) THEN
            IF((ITRIG(KS(1),KS(4),NU+NU+3) /= 0) .AND.  &
               (ITRIG(KS(2),KS(3),NU+NU+3) /= 0)) THEN
              IF(NU >=  0) THEN
                N=(NU-NE1)/2+1
                IF(N < NE2) THEN
                  CALL CXK(S,IS,KAPS,NU,KRA,1,2)
                  DO MU = 1,12
                    CONE(MU,N)=CONE(MU,N)+AB*S(MU)
                  END DO
                END IF
              END IF
            END IF
          END IF
        END IF
      END DO
      IF (ICOLBREI == 2) THEN
        DO N = 1,NE2
          NU=NE1+2*(N-1)
          CALL TALK(JJJA,JJJB,NU,IS(1),IS(4),IS(2),IS(3),1,CONE(1,N))
          CALL TALK(JJJA,JJJB,NU,IS(4),IS(1),IS(3),IS(2),1,CONE(2,N))
          CALL TALK(JJJA,JJJB,NU,IS(1),IS(4),IS(3),IS(2),1,CONE(3,N))
          CALL TALK(JJJA,JJJB,NU,IS(4),IS(1),IS(2),IS(3),1,CONE(4,N))
          IF(N == NE2) CYCLE
          NUP1=NU+1
          CALL TALK(JJJA,JJJB,NUP1,IS(1),IS(4),IS(2),IS(3),2,CONE(5,N))
          CALL TALK(JJJA,JJJB,NUP1,IS(2),IS(3),IS(1),IS(4),2,CONE(6,N))
          CALL TALK(JJJA,JJJB,NUP1,IS(4),IS(1),IS(3),IS(2),2,CONE(7,N))
          CALL TALK(JJJA,JJJB,NUP1,IS(3),IS(2),IS(4),IS(1),2,CONE(8,N))
          CALL TALK(JJJA,JJJB,NUP1,IS(1),IS(4),IS(3),IS(2),2,CONE(9,N))
          CALL TALK(JJJA,JJJB,NUP1,IS(3),IS(2),IS(1),IS(4),2,CONE(10,N))
          CALL TALK(JJJA,JJJB,NUP1,IS(4),IS(1),IS(2),IS(3),2,CONE(11,N))
          CALL TALK(JJJA,JJJB,NUP1,IS(2),IS(3),IS(4),IS(1),2,CONE(12,N))
        END DO
      END IF
      IF(IREZ == 2)GO TO 6
    7 CONTINUE
      RETURN
   10 WRITE(99,100)
  100 FORMAT(5X,'ERRO IN EL33  PMGG RAGG')
      STOP
      END SUBROUTINE EL33
      !
      !

      
!*******************************************************************
!                                                                  *
      SUBROUTINE EL4(JJA,JJB,JA,JB,JC,JD,ICOLBREI)
!                                                                  *
!   --------------  SECTION METWO    SUBPROGRAM 09  -------------  *
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF TWO PARTICLE OPERATOR IN CASE :      N'1 = N1 +- 1        *
!                                             N'2 = N2 +- 1        *
!                                             N'3 = N3 -+ 2        *
!                                                                  *
!   Written by  G. Gaigalas                                        *
!   The last modification made by G. Gaigalas       October  2017  *
!                                                                  *
!*******************************************************************
!
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: JJA,JJB,JA,JB,JC,JD,ICOLBREI
!-----------------------------------------------
      IF(NPEEL <= 2)RETURN
      IF(JA == JB) THEN
        CALL EL41(JJA,JJB,JC,JD,JA,1,JA,JB,JC,JD,ICOLBREI)
      ELSE IF(JC == JD) THEN
        CALL EL41(JJA,JJB,JA,JB,JC,2,JA,JB,JC,JD,ICOLBREI)
      ELSE
        WRITE(99,100)
        STOP
      END IF
      RETURN
  100 FORMAT(5X,'ERRO IN EL4 ')
      END SUBROUTINE EL4
      !
      !
!*******************************************************************
!                                                                  *
      SUBROUTINE EL41(JJJA,JJJB,JA,JB,JC,IREZ,JJA,JJB,JJC,JJD,    &
                                                          ICOLBREI)
!                                                                  *
!   --------------  SECTION METWO    SUBPROGRAM 10  -------------  *
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF TWO PARTICLE OPERATOR IN CASE :       N'1 = N1 + 1        *
!                                              N'2 = N2 + 1        *
!                                              N'3 = N3 - 2,       *
!     WHEN IREZ = 1   . . . . . . . . . . . . . . . . . . .        *
!                                              N'1 = N1 - 1        *
!                                              N'2 = N2 - 1        *
!                                              N'3 = N3 + 2,       *
!     WHEN IREZ = 2   . . . . . . . . . . . . . . . . . . .        *
!                                                                  *
!                                                                  *
!   Written by  G. Gaigalas                                        *
!   The last modification made by G. Gaigalas       October  2017  *
!                                                                  *
!*******************************************************************
!
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: JJJA,JJJB,JA,JB,JC,IREZ,JJA,JJB, &
                             JJC,JJD,ICOLBREI
!      DIMENSION J(3)
!      DIMENSION COND(12,20),S(12),IS(4),KAPS(4),KS(4)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER ::  IA,IB,IC,II,IIA,IIB,IIC,IBRD,IBRE,IP1,IG1,IP2,IG2, &
                  IAT,IID,IKK,IFAZ,IFAZP,IFAZFRCS,INN,I2,I3,JB1,JAA, &
                  JBB,JCC,J12,KRA,L1,L2,L3,ND1,ND2,NE1,NE2,N,NN,NU,  &
                  NUP1,MU
      INTEGER, DIMENSION(3) :: J
      INTEGER, DIMENSION(4) :: IS,KAPS,KS
      REAL(KIND=8)          :: AA,AB,A1,BB,QM1,QM2,QM3,QM4,SI,RECC
      REAL(KIND=8), DIMENSION(12) :: S
      REAL(KIND=8), DIMENSION(12,20) :: COND
!-----------------------------------------------
      CALL EILE(JA,JB,JC,JAA,JBB,JCC)
      IF(NPEEL <= 1)RETURN
      CALL RECO(JAA,JCC,JBB,JBB,2,IAT)
      IF(IAT == 0)RETURN
      IA=JLIST(JA)
      IB=JLIST(JB)
      IC=JLIST(JC)
      IIA=JLIST(JJA)
      IIB=JLIST(JJB)
      IIC=JLIST(JJC)
      IID=JLIST(JJD)
      IF(IREZ == 1) THEN
        QM1=-HALF
        QM2=-HALF
        QM3=HALF
        QM4=HALF
      ELSE
        QM1=HALF
        QM2=HALF
        QM3=-HALF
        QM4=-HALF
      END IF
      CALL PERKO2(JA,JB,JC,JA,3)
      J(1)=ID1(3)
      J(2)=ID2(3)
      J(3)=ID3(3)
      L1=(J(1)+1)/2
      L2=(J(2)+1)/2
      L3=(J(3)+1)/2
      IF (ICOLBREI == 2) THEN
        IS(1)=IIA
        IS(2)=IIB
        IS(3)=IIC
        IS(4)=IID
        KAPS(1:4)=2*NAK(IS(1:4))
        KS(1:4)=IABS(KAPS(1:4))
        CALL SNRC(IS,KAPS,KS,ND1,ND2,NE1,NE2,IBRD,IBRE)
        IF(IBRD <= 0)RETURN
        COND = Zero
      END IF
      IFAZP=JFAZE(JC,JA,JB,JC)
      IFAZFRCS = 1
!
!     TRANSFORM FANO & RACAH PHASE CONVENTION
!     TO CONDON & SHORTLEY PHASE CONVENTION
!
      IFAZ=IK1(5)*IK1(4)+IK2(5)*IK2(4)-ID1(5)*ID1(4)-ID2(5)*ID2(4)+ &
      IK3(5)*IK3(4)-ID3(5)*ID3(4)
      IF((IFAZ/4)*4 /= IFAZ)IFAZFRCS=-IFAZFRCS
!
      IF(JA > JB) THEN
        JAA=JB
        JBB=JA
      ELSE
        JAA=JA
        JBB=JB
      END IF
      NN=0
      JB1=JBB-1
      DO II=JAA,JB1
        INN=JLIST(II)
        NN=NQ1(INN)+NN
      END DO
      IF((NN/2)*2 == NN)IFAZP=-IFAZP
! * * *                      * * *                      * * *
!     CASES 3312   + + - -        TRANSFORM TO  1233   - - + +
!           3321                                1233
!                                                    (IREZ = 1)
!     OR
!     CASES 1233   + + - -        TRANSFORM TO  1233   + + - -
!           2133                                1233
!                                                    (IREZ = 2)
      IP1=ITREXG(J(2),J(1),J(3),J(3),IKK)+1
      IF(IKK <= 0)RETURN
      IG1=IP1+IKK-1
      IP2=ITREXG(J(3),J(1),J(2),J(3),IKK)+1
      IF(IKK <= 0) RETURN
      IG2=IP2+IKK-1
      DO I2=IP2,IG2,2
        KRA=(I2-1)/2
!
        IF (ICOLBREI == 1) THEN
          IF(IREZ == 2) THEN
            CALL COULOM(L1,L2,L3,L3,ID1(5),ID2(5),ID3(5),ID3(5),KRA,A1)
          ELSE
            CALL COULOM(L3,L3,L1,L2,ID3(5),ID3(5),ID1(5),ID2(5),KRA,A1)
          END IF
          IF(DABS(A1) < EPS) CYCLE
        END IF
!
        AB=ZERO
        DO I3=IP1,IG1,2
          J12=(I3-1)/2
          IFAZ=J(2)-J12+1
          IF(IREZ == 2)IFAZ=J(1)-J12+1
          IF((IFAZ/2)*2 /= IFAZ) CYCLE
          CALL REC3(JA,JB,JC,J(1),J(2),J12*2,0,IAT,AA)
          IF(IAT == 0) CYCLE
          IF(IXJTIK(J(3),J(1),KRA*2,J(2),J(3),J12*2) == 0) CYCLE
          CALL GG1233(IK1,IK2,IK3,BK1,BK2,BK3,ID1,ID2,ID3,BD1,    &
                    BD2,BD3,J12,QM1,QM2,QM3,QM4,AA)
          IF(DABS(AA) < EPS) CYCLE
          CALL REC3(JA,JB,JC,J(1),J(2),J12*2,1,IAT,RECC)
          AA=AA*RECC
          CALL SIXJ(J(3),J(1),KRA*2,J(2),J(3),J12*2,0,SI)
          AA=AA*SI*DSQRT(DBLE(I3))
          IFAZ=J(3)+J(1)+2*J12+2*KRA
          IF(IREZ == 2)IFAZ=J(2)+J(3)+2*J12+2*KRA
          IF((IFAZ/4)*4 /= IFAZ)AA=-AA
          AB=AB+AA
        END DO
        IF(DABS(AB) < EPS) CYCLE
        AB=-AB*DBLE(IFAZP)
        IF (ICOLBREI == 1) THEN
          BB=A1*AB*DBLE(IFAZFRCS)
          CALL SPEAK(JJJA,JJJB,IIA,IIB,IIC,IID,KRA,BB)
        ELSE IF (ICOLBREI == 2) THEN
          NU=KRA
          IF(((NU-ND1)/2)*2 == (NU-ND1)) THEN
            IF((ITRIG(KS(1),KS(3),NU+NU+1) /= 0) .AND.  &
               (ITRIG(KS(2),KS(4),NU+NU+1) /= 0)) THEN
              N=(NU-ND1)/2+1
              IF(NU > 0) THEN
                CALL CXK(S,IS,KAPS,NU,KRA,1,1)
                DO MU = 1,4
                  COND(MU,N)=COND(MU,N)+AB*S(MU)
                END DO
              END IF
            END IF
          END IF
          NU=KRA+1
          IF(((NU-ND1)/2)*2 == (NU-ND1)) THEN
            IF((ITRIG(KS(1),KS(3),NU+NU-1) /= 0) .AND.  &
               (ITRIG(KS(2),KS(4),NU+NU-1) /= 0)) THEN
              N=(NU-ND1)/2+1
              IF(N <= ND2) THEN
                CALL CXK(S,IS,KAPS,NU,KRA,1,1)
                DO MU = 1,4
                  COND(MU,N)=COND(MU,N)+AB*S(MU)
                END DO
              END IF
            END IF
          END IF
          NU=KRA-1
          IF(((NU-ND1)/2)*2 == (NU-ND1)) THEN
            IF((ITRIG(KS(1),KS(3),NU+NU+3) /= 0) .AND.  &
               (ITRIG(KS(2),KS(4),NU+NU+3) /= 0)) THEN
              IF(NU >= 0) THEN
                N=(NU-ND1)/2+1
                IF(N < ND2) THEN
                  CALL CXK(S,IS,KAPS,NU,KRA,1,1)
                  DO MU = 1,12
                    COND(MU,N)=COND(MU,N)+AB*S(MU)
                  END DO
                END IF
              END IF
            END IF
          END IF
        END IF
      END DO
      IF (ICOLBREI == 2) THEN
        DO N = 1,ND2
          NU=ND1+2*(N-1)
          CALL TALK(JJJA,JJJB,NU,IS(1),IS(3),IS(2),IS(4),1,COND(1,N))
          CALL TALK(JJJA,JJJB,NU,IS(3),IS(1),IS(4),IS(2),1,COND(2,N))
          CALL TALK(JJJA,JJJB,NU,IS(1),IS(3),IS(4),IS(2),1,COND(3,N))
          CALL TALK(JJJA,JJJB,NU,IS(3),IS(1),IS(2),IS(4),1,COND(4,N))
          IF(N == ND2) CYCLE
          NUP1=NU+1
          CALL TALK(JJJA,JJJB,NUP1,IS(1),IS(3),IS(2),IS(4),2,COND(5,N))
          CALL TALK(JJJA,JJJB,NUP1,IS(2),IS(4),IS(1),IS(3),2,COND(6,N))
          CALL TALK(JJJA,JJJB,NUP1,IS(3),IS(1),IS(4),IS(2),2,COND(7,N))
          CALL TALK(JJJA,JJJB,NUP1,IS(4),IS(2),IS(3),IS(1),2,COND(8,N))
          CALL TALK(JJJA,JJJB,NUP1,IS(1),IS(3),IS(4),IS(2),2,COND(9,N))
          CALL TALK(JJJA,JJJB,NUP1,IS(4),IS(2),IS(1),IS(3),2,COND(10,N))
          CALL TALK(JJJA,JJJB,NUP1,IS(3),IS(1),IS(2),IS(4),2,COND(11,N))
          CALL TALK(JJJA,JJJB,NUP1,IS(2),IS(4),IS(3),IS(1),2,COND(12,N))
        END DO
      END IF
      RETURN
   10 WRITE(99,100)
  100 FORMAT(5X,'ERRO IN EL41  PMGG RAGG')
      STOP
      END SUBROUTINE EL41
      !
      !
!*******************************************************************
!                                                                  *
      SUBROUTINE EL5(JJA,JJB,JA,JB,JC,JD,ICOLBREI)
!                                                                  *
!   --------------  SECTION METWO    SUBPROGRAM 11  -------------  *
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF TWO PARTICLE OPERATOR IN CASE :    N'1 = N1 (+-) 1        *
!                                           N'2 = N2 (+-) 1        *
!                                           N'3 = N3 (+-) 1        *
!                                           N'4 = N4 (+-) 1        *
!                                                                  *
!   Written by  G. Gaigalas                                        *
!   The last modification made by G. Gaigalas       October  2017  *
!                                                                  *
!*******************************************************************
!
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: JJA,JJB,JA,JB,JC,JD,ICOLBREI
!-----------------------------------------------
      IF(NPEEL <= 3)RETURN
      IF(JB < JC) THEN
        CALL EL51(JJA,JJB,JA,JB,JC,JD,1,ICOLBREI)
      ELSE IF(JA > JD.AND.JB > JD) THEN
        CALL EL51(JJA,JJB,JC,JD,JA,JB,2,ICOLBREI)
      ELSE IF(JB > JC.AND.JB < JD.AND.JA < JC) THEN
        CALL EL52(JJA,JJB,JA,JC,JB,JD,1,ICOLBREI)
      ELSE IF(JB > JC.AND.JB > JD.AND.JA > JC) THEN
        CALL EL52(JJA,JJB,JC,JA,JD,JB,2,ICOLBREI)
      ELSE IF(JB > JC.AND.JB > JD.AND.JA < JC) THEN
        CALL EL53(JJA,JJB,JA,JC,JD,JB,1,ICOLBREI)
      ELSE IF(JB > JC.AND.JB < JD.AND.JA > JC) THEN
        CALL EL53(JJA,JJB,JC,JA,JB,JD,2,ICOLBREI)
      ELSE
        WRITE(99,100)
        STOP
      END IF
      RETURN
  100 FORMAT(5X,'ERRO IN EL5 ')
      END SUBROUTINE EL5
      !
      !
!*******************************************************************
!                                                                  *
      SUBROUTINE EL51(JJA,JJB,JA,JB,JC,JD,IREZ,ICOLBREI)
!                                                                  *
!   --------------  SECTION METWO    SUBPROGRAM 12  -------------  *
!                                                                  *
!     THIS PACKAGE EVALUATED THE CASES - 1234, 2134, 1243, 2134    *
!                                                   ( IREZ = 1),   *
!                                                   ( + + - - ),   *
!     WHICH APPEARS IN CALCULATION MATRIX ELEMENTS BETWEEN         *
!     CONFIGURATIONS:                               N'1 = N1 - 1   *
!                                                   N'2 = N2 - 1   *
!                                                   N'3 = N3 + 1   *
!                                                   N'4 = N4 + 1   *
!     AND    3412, 4321, 3421, 4312                 ( IREZ = 2),   *
!                                                   ( + + - - ),   *
!     WHICH APPEARS IN CALCULATION MATRIX ELEMENTS BETWEEN         *
!     CONFIGURATIONS:                               N'1 = N1 + 1   *
!                                                   N'2 = N2 + 1   *
!                                                   N'3 = N3 - 1   *
!                                                   N'4 = N4 - 1   *
!                                                                  *
!   Written by  G. Gaigalas                                        *
!   The last modification made by G. Gaigalas       October  2017  *
!                                                                  *
!*******************************************************************
!
        IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: JJA,JJB,JA,JB,JC,JD,IREZ,ICOLBREI
!      DIMENSION PMGG(30),J(4)
!      DIMENSION COND(12,20),CONE(12,20),S(12),IS(4),KAPS(4),KS(4)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IA,IB,IC,ID,IBRD,IBRE,II,IP1,IP2,IG1,IG2,IKK,I2,I3,  &
                 I4,IFAZ,IFAZP,IFAZFRCS,INN,IAT,KRA,KRA1,L1,L2,L3,L4, &
                 J12,JB1,JD1,ND1,ND2,NE1,NE2,N,NN,NU,NUP1,MU
      INTEGER, DIMENSION(4) :: J
      INTEGER, DIMENSION(4) :: IS,KAPS,KS
      REAL(KIND=8)          :: AA,AB,A1,BB,QM1,QM2,QM3,QM4,RAG,RECC,SI
      REAL(KIND=8), DIMENSION(12) :: S
      REAL(KIND=8), DIMENSION(30) :: PMGG
      REAL(KIND=8), DIMENSION(12,20) :: COND,CONE
!-----------------------------------------------
      IF(NPEEL <= 3)RETURN
      CALL RECO(JA,JD,JC,JB,3,IAT)
      IF(IAT == 0)RETURN
      IA=JLIST(JA)
      IB=JLIST(JB)
      IC=JLIST(JC)
      ID=JLIST(JD)
      IF(IREZ == 2) THEN
        QM1=-HALF
        QM2=-HALF
        QM3=HALF
        QM4=HALF
      ELSE
        QM1=HALF
        QM2=HALF
        QM3=-HALF
        QM4=-HALF
      END IF
      CALL PERKO2(JA,JB,JC,JD,4)
      J(1)=ID1(3)
      J(2)=ID2(3)
      J(3)=ID3(3)
      J(4)=ID4(3)
      L1=(J(1)+1)/2
      L2=(J(2)+1)/2
      L3=(J(3)+1)/2
      L4=(J(4)+1)/2
      IF (ICOLBREI == 2) THEN
        IF(IREZ == 1) THEN
          IS(1)=IA
          IS(2)=IB
          IS(3)=IC
          IS(4)=ID
        ELSE IF(IREZ == 2) THEN
          IS(1)=IC
          IS(2)=ID
          IS(3)=IA
          IS(4)=IB
        END IF
        KAPS(1:4)=2*NAK(IS(1:4))
        KS(1:4)=IABS(KAPS(1:4))
        CALL SNRC(IS,KAPS,KS,ND1,ND2,NE1,NE2,IBRD,IBRE)
        IF(IBRD <= 0 .AND. IBRE <= 0)RETURN
        COND = Zero;  CONE = Zero
      END IF
      CALL GG1234(IK1,IK2,IK3,IK4,BK1,BK2,BK3,BK4,ID1,ID2,   &
                  ID3,ID4,BD1,BD2,BD3,BD4,QM1,QM2,QM3,QM4,RAG)
      IF(DABS(RAG) < EPS) RETURN
      IP1=ITREXG(J(1),J(2),J(3),J(4),IKK)+1
      IF(IKK <= 0)RETURN
      IG1=IP1+IKK-1
      DO I4=IP1,IG1,2
        KRA=(I4-1)/2
        KRA1=KRA+1
        IF(KRA1 > 30)GO TO 10
        PMGG(KRA1)=ZERO
        CALL RECO4(JA,JB,JC,JD,J(1),J(2),J(3),J(4),KRA*2,0,IAT,RECC)
        IF(IAT == 0) CYCLE
        CALL RECO4(JA,JB,JC,JD,J(1),J(2),J(3),J(4),KRA*2,1,IAT,RECC)
        PMGG(KRA1)=RECC
      END DO
!
!     TRANSFORM FANO & RACAH PHASE CONVENTION
!     TO CONDON & SHORTLEY PHASE CONVENTION
!
      IFAZFRCS=1
      IFAZ=IK1(5)*IK1(4)+IK2(5)*IK2(4)-ID1(5)*ID1(4)-ID2(5)*ID2(4)  &
      +IK3(5)*IK3(4)-ID3(5)*ID3(4)+IK4(5)*IK4(4)-ID4(5)*ID4(4)
      IF((IFAZ/4)*4 /= IFAZ)IFAZFRCS=-IFAZFRCS
!
      NN=0
      JB1=JB-1
      IFAZP=1
      DO II=JA,JB1
        INN=JLIST(II)
        NN=NQ1(INN)+NN
      END DO
      IF((NN/2)*2 == NN)IFAZP=-IFAZP
      NN=0
      JD1=JD-1
      DO II=JC,JD1
        INN=JLIST(II)
        NN=NQ1(INN)+NN
      END DO
      IF((NN/2)*2 == NN)IFAZP=-IFAZP
! * * *                      * * *                      * * *
!     CASES 1234   + + - -
!           2134                  TRANSFORM TO  1234   + + - -
!                                                    (IREZ = 1)
!     OR
!     CASES 3412   + + - -        TRANSFORM TO  1234   - - + +
!           3421                                1234
!                                                    (IREZ = 2)
      IP2=ITREXG(J(1),J(3),J(2),J(4),IKK)+1
      IF(IKK <= 0) RETURN
      IG2=IP2+IKK-1
      DO I2=IP2,IG2,2
        KRA=(I2-1)/2
!
        IF (ICOLBREI == 1) THEN
          IF(IREZ == 2) THEN
            CALL COULOM(L3,L4,L1,L2,ID3(5),ID4(5),ID1(5),ID2(5),KRA,A1)
          ELSE
            CALL COULOM(L1,L2,L3,L4,ID1(5),ID2(5),ID3(5),ID4(5),KRA,A1)
          END IF
          IF(DABS(A1) < EPS) CYCLE
        END IF
!
        AB=ZERO
        DO I3=IP1,IG1,2
          J12=(I3-1)/2
          IFAZ=J(1)+J(4)-2*J12
          IF(IREZ == 2)IFAZ=J(2)+J(3)-2*J12
          IF((IFAZ/2)*2 /= IFAZ) CYCLE
          KRA1=J12+1
          IF(KRA1 > 30)GO TO 10
          AA=PMGG(KRA1)
          IF(DABS(AA) < EPS) CYCLE
          IF(IXJTIK(J(1),J(3),KRA*2,J(4),J(2),J12*2) == 0) CYCLE
          AA=AA*RAG
          IF(DABS(AA) < EPS) CYCLE
          CALL SIXJ(J(1),J(3),KRA*2,J(4),J(2),J12*2,0,SI)
          AA=AA*SI*DSQRT(DBLE(I3))
          IFAZ=J(2)+J(3)+2*J12+2*KRA
          IF(IREZ == 2)IFAZ=J(1)+J(4)-2*J12+2*KRA
          IF((IFAZ/4)*4 /= IFAZ)AA=-AA
          AB=AB+AA
        END DO
        IF(DABS(AB) < EPS) CYCLE
        AB=-AB*DBLE(IFAZP)
        IF (ICOLBREI == 1) THEN
          BB=A1*AB*DBLE(IFAZFRCS)
          IF(IREZ == 1)CALL SPEAK(JJA,JJB,IA,IB,IC,ID,KRA,BB)
          IF(IREZ == 2)CALL SPEAK(JJA,JJB,IC,ID,IA,IB,KRA,BB)
        ELSE IF (ICOLBREI == 2) THEN
          NU=KRA
          IF(((NU-ND1)/2)*2 == (NU-ND1)) THEN
            IF((ITRIG(KS(1),KS(3),NU+NU+1) /= 0) .AND.  &
               (ITRIG(KS(2),KS(4),NU+NU+1) /= 0)) THEN
              N=(NU-ND1)/2+1
              IF(NU > 0) THEN
                CALL CXK(S,IS,KAPS,NU,KRA,1,1)
                DO MU = 1,4
                  COND(MU,N)=COND(MU,N)+AB*S(MU)
                END DO
              END IF
            END IF
          END IF
          NU=KRA+1
          IF(((NU-ND1)/2)*2 == (NU-ND1)) THEN
            IF((ITRIG(KS(1),KS(3),NU+NU-1) /= 0) .AND.  &
               (ITRIG(KS(2),KS(4),NU+NU-1) /= 0)) THEN
              N=(NU-ND1)/2+1
              IF(N <= ND2) THEN
                CALL CXK(S,IS,KAPS,NU,KRA,1,1)
                DO MU = 1,4
                  COND(MU,N)=COND(MU,N)+AB*S(MU)
                END DO
              END IF
            END IF
          END IF
          NU=KRA-1
          IF(((NU-ND1)/2)*2 == (NU-ND1)) THEN
            IF((ITRIG(KS(1),KS(3),NU+NU+3) /= 0) .AND.  &
               (ITRIG(KS(2),KS(4),NU+NU+3) /= 0)) THEN
              IF(NU >= 0) THEN
                N=(NU-ND1)/2+1
                IF(N < ND2) THEN
                  CALL CXK(S,IS,KAPS,NU,KRA,1,1)
                  DO MU = 1,12
                    COND(MU,N)=COND(MU,N)+AB*S(MU)
                  END DO
                END IF
              END IF
            END IF
          END IF
        END IF
      END DO
      IF (ICOLBREI == 2) THEN
        DO N = 1,ND2
          NU=ND1+2*(N-1)
          CALL TALK(JJA,JJB,NU,IS(1),IS(3),IS(2),IS(4),1,COND(1,N))
          CALL TALK(JJA,JJB,NU,IS(3),IS(1),IS(4),IS(2),1,COND(2,N))
          CALL TALK(JJA,JJB,NU,IS(1),IS(3),IS(4),IS(2),1,COND(3,N))
          CALL TALK(JJA,JJB,NU,IS(3),IS(1),IS(2),IS(4),1,COND(4,N))
          IF(N == ND2) CYCLE
          NUP1=NU+1
          CALL TALK(JJA,JJB,NUP1,IS(1),IS(3),IS(2),IS(4),2,COND(5,N))
          CALL TALK(JJA,JJB,NUP1,IS(2),IS(4),IS(1),IS(3),2,COND(6,N))
          CALL TALK(JJA,JJB,NUP1,IS(3),IS(1),IS(4),IS(2),2,COND(7,N))
          CALL TALK(JJA,JJB,NUP1,IS(4),IS(2),IS(3),IS(1),2,COND(8,N))
          CALL TALK(JJA,JJB,NUP1,IS(1),IS(3),IS(4),IS(2),2,COND(9,N))
          CALL TALK(JJA,JJB,NUP1,IS(4),IS(2),IS(1),IS(3),2,COND(10,N))
          CALL TALK(JJA,JJB,NUP1,IS(3),IS(1),IS(2),IS(4),2,COND(11,N))
          CALL TALK(JJA,JJB,NUP1,IS(2),IS(4),IS(3),IS(1),2,COND(12,N))
        END DO
      END IF
! * * *                      * * *                      * * *
!     CASES 1243   + + - -        TRANSFORM TO  1234   + + - -
!           2134                                1234
!                                                    (IREZ = 1)
!     OR
!     CASES 3421   + + - -        TRANSFORM TO  1234   - - + +
!           4321                                1234
!                                                    (IREZ = 2)
      IP2=ITREXG(J(1),J(4),J(2),J(3),IKK)+1
      IF(IKK <= 0) RETURN
      IG2=IP2+IKK-1
      DO I2=IP2,IG2,2
        KRA=(I2-1)/2
!
        IF (ICOLBREI == 1) THEN
          IF(IREZ == 2) THEN
            CALL COULOM(L3,L4,L2,L1,ID3(5),ID4(5),ID2(5),ID1(5),KRA,A1)
          ELSE
            CALL COULOM(L1,L2,L4,L3,ID1(5),ID2(5),ID4(5),ID3(5),KRA,A1)
          END IF
          IF(DABS(A1) < EPS) CYCLE
        END IF
!
        AB=ZERO
        DO I3=IP1,IG1,2
          J12=(I3-1)/2
          IFAZ=J(1)-J(4)
          IF(IREZ == 2)IFAZ=J(3)-J(2)
          IF((IFAZ/2)*2 /= IFAZ) CYCLE
          KRA1=J12+1
          IF(KRA1 > 30)GO TO 10
          AA=PMGG(KRA1)
          IF(DABS(AA) < EPS) CYCLE
          AA=AA*RAG
          IF(DABS(AA) < EPS) CYCLE
          IF(IXJTIK(J(1),J(4),KRA*2,J(3),J(2),J12*2) == 0) CYCLE
          CALL SIXJ(J(1),J(4),KRA*2,J(3),J(2),J12*2,0,SI)
          AA=AA*SI*DSQRT(DBLE(I3))
          IFAZ=J(2)+J(3)+2*J(4)+2*KRA
          IF(IREZ == 2)IFAZ=J(1)+2*J(2)+J(4)+4*J12+2*KRA
          IF((IFAZ/4)*4 /= IFAZ)AA=-AA
          AB=AB+AA
        END DO
        IF(DABS(AB) < EPS) CYCLE
        AB=AB*DBLE(IFAZP)
        IF (ICOLBREI == 1) THEN
          BB=A1*AB*DBLE(IFAZFRCS)
          IF(IREZ == 1)CALL SPEAK(JJA,JJB,IA,IB,ID,IC,KRA,BB)
          IF(IREZ == 2)CALL SPEAK(JJA,JJB,IC,ID,IB,IA,KRA,BB)
        ELSE IF (ICOLBREI == 2) THEN
          NU=KRA
          IF(((NU-NE1)/2)*2 == (NU-NE1)) THEN
            IF((ITRIG(KS(1),KS(4),NU+NU+1) /= 0) .AND.   &
               (ITRIG(KS(2),KS(3),NU+NU+1) /= 0)) THEN
              N=(NU-NE1)/2+1
              IF(NU > 0) THEN
                CALL CXK(S,IS,KAPS,NU,KRA,1,2)
                DO MU = 1,4
                  CONE(MU,N)=CONE(MU,N)+AB*S(MU)
                END DO
              END IF
            END IF
          END IF
          NU=KRA+1
          IF(((NU-NE1)/2)*2 == (NU-NE1)) THEN
            IF((ITRIG(KS(1),KS(4),NU+NU-1) /= 0) .AND.  &
               (ITRIG(KS(2),KS(3),NU+NU-1) /= 0)) THEN
              N=(NU-NE1)/2+1
              IF(N <= NE2) THEN
                CALL CXK(S,IS,KAPS,NU,KRA,1,2)
                DO MU = 1,4
                  CONE(MU,N)=CONE(MU,N)+AB*S(MU)
                END DO
              END IF
            END IF
          END IF
          NU=KRA-1
          IF(((NU-NE1)/2)*2 == (NU-NE1)) THEN
            IF((ITRIG(KS(1),KS(4),NU+NU+3) /= 0) .AND.  &
               (ITRIG(KS(2),KS(3),NU+NU+3) /= 0)) THEN
              IF(NU >= 0) THEN
                N=(NU-NE1)/2+1
                IF(N < NE2) THEN
                  CALL CXK(S,IS,KAPS,NU,KRA,1,2)
                  DO MU = 1,12
                    CONE(MU,N)=CONE(MU,N)+AB*S(MU)
                  END DO
                END IF
              END IF
            END IF
          END IF
        END IF
      END DO
      IF (ICOLBREI == 2) THEN
        DO N = 1,NE2
          NU=NE1+2*(N-1)
          CALL TALK(JJA,JJB,NU,IS(1),IS(4),IS(2),IS(3),1,CONE(1,N))
          CALL TALK(JJA,JJB,NU,IS(4),IS(1),IS(3),IS(2),1,CONE(2,N))
          CALL TALK(JJA,JJB,NU,IS(1),IS(4),IS(3),IS(2),1,CONE(3,N))
          CALL TALK(JJA,JJB,NU,IS(4),IS(1),IS(2),IS(3),1,CONE(4,N))
          IF(N == NE2) CYCLE
          NUP1=NU+1
          CALL TALK(JJA,JJB,NUP1,IS(1),IS(4),IS(2),IS(3),2,CONE(5,N))
          CALL TALK(JJA,JJB,NUP1,IS(2),IS(3),IS(1),IS(4),2,CONE(6,N))
          CALL TALK(JJA,JJB,NUP1,IS(4),IS(1),IS(3),IS(2),2,CONE(7,N))
          CALL TALK(JJA,JJB,NUP1,IS(3),IS(2),IS(4),IS(1),2,CONE(8,N))
          CALL TALK(JJA,JJB,NUP1,IS(1),IS(4),IS(3),IS(2),2,CONE(9,N))
          CALL TALK(JJA,JJB,NUP1,IS(3),IS(2),IS(1),IS(4),2,CONE(10,N))
          CALL TALK(JJA,JJB,NUP1,IS(4),IS(1),IS(2),IS(3),2,CONE(11,N))
          CALL TALK(JJA,JJB,NUP1,IS(2),IS(3),IS(4),IS(1),2,CONE(12,N))
        END DO
      END IF
      RETURN
   10 WRITE(99,100)
  100 FORMAT(5X,'ERRO IN EL51  PMGG RAGG')
      STOP
      END SUBROUTINE EL51
      !
      !
!*******************************************************************
!                                                                  *
      SUBROUTINE EL52(JJA,JJB,JA,JB,JC,JD,IREZ,ICOLBREI)
!                                                                  *
!   --------------  SECTION METWO    SUBPROGRAM 13  -------------  *
!                                                                  *
!     THIS PACKAGE EVALUATED THE CASES - 1324, 3142, 1342, 3124    *
!                                                   ( IREZ = 1),   *
!                                                   ( + + - - ),   *
!     WHICH APPEARS IN CALCULATION MATRIX ELEMENTS BETWEEN         *
!     CONFIGURATIONS:                               N'1 = N1 - 1   *
!                                                   N'2 = N2 + 1   *
!                                                   N'3 = N3 - 1   *
!                                                   N'4 = N4 + 1   *
!     AND    2413, 4231, 2431, 4213                 ( IREZ = 2),   *
!                                                   ( + + - - ),   *
!     WHICH APPEARS IN CALCULATION MATRIX ELEMENTS BETWEEN         *
!     CONFIGURATIONS:                               N'1 = N1 + 1   *
!                                                   N'2 = N2 - 1   *
!                                                   N'3 = N3 + 1   *
!                                                   N'4 = N4 - 1   *
!                                                                  *
!                                                                  *
!   Written by  G. Gaigalas                                        *
!   The last modification made by G. Gaigalas       October  2017  *
!                                                                  *
!*******************************************************************
!
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: JJA,JJB,JA,JB,JC,JD,IREZ,ICOLBREI
!      DIMENSION PMGG(30),J(4)
!     DIMENSION COND(12,20),CONE(12,20),S(12),IS(4),KAPS(4),KS(4)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IA,IB,IC,ID,IBRD,IBRE,II,IP1,IP2,IG1,IG2,IKK,I2,I3,  &
                 I4,IFAZ,IFAZP,IFAZFRCS,INN,IAT,KRA,KRA1,L1,L2,L3,L4, &
                 J12,JB1,JD1,ND1,ND2,NE1,NE2,N,NN,NU,NUP1,MU
      INTEGER, DIMENSION(4) :: J
      INTEGER, DIMENSION(4) :: IS,KAPS,KS
      REAL(KIND=8)          :: AA,AB,A1,BB,QM1,QM2,QM3,QM4,RAG,RECC,SI
      REAL(KIND=8), DIMENSION(12) :: S
      REAL(KIND=8), DIMENSION(30) :: PMGG
      REAL(KIND=8), DIMENSION(12,20) :: COND,CONE
!-----------------------------------------------
      IF(NPEEL <= 3)RETURN
      CALL RECO(JA,JD,JC,JB,3,IAT)
      IF(IAT == 0)RETURN
      IA=JLIST(JA)
      IB=JLIST(JB)
      IC=JLIST(JC)
      ID=JLIST(JD)
      IF(IREZ == 2) THEN
        QM1=-HALF
        QM2=HALF
        QM3=-HALF
        QM4=HALF
      ELSE
        QM1=HALF
        QM2=-HALF
        QM3=HALF
        QM4=-HALF
      END IF
      CALL PERKO2(JA,JB,JC,JD,4)
      J(1)=ID1(3)
      J(2)=ID2(3)
      J(3)=ID3(3)
      J(4)=ID4(3)
      L1=(J(1)+1)/2
      L2=(J(2)+1)/2
      L3=(J(3)+1)/2
      L4=(J(4)+1)/2
      IF (ICOLBREI == 2) THEN
        IF(IREZ == 1) THEN
          IS(1)=IA
          IS(2)=IC
          IS(3)=IB
          IS(4)=ID
        ELSE IF(IREZ == 2) THEN
          IS(1)=IB
          IS(2)=ID
          IS(3)=IA
          IS(4)=IC
        END IF
        KAPS(1:4)=2*NAK(IS(1:4))
        KS(1:4)=IABS(KAPS(1:4))
        CALL SNRC(IS,KAPS,KS,ND1,ND2,NE1,NE2,IBRD,IBRE)
        IF(IBRD <= 0 .AND. IBRE <= 0)RETURN
        COND = Zero; CONE = Zero
      END IF
      CALL GG1234(IK1,IK2,IK3,IK4,BK1,BK2,BK3,BK4,ID1,ID2,  &
      ID3,ID4,BD1,BD2,BD3,BD4,QM1,QM2,QM3,QM4,RAG)
      IF(DABS(RAG) < EPS) RETURN
      IP1=ITREXG(J(1),J(2),J(3),J(4),IKK)+1
      IF(IKK <= 0)RETURN
      IG1=IP1+IKK-1
      DO I4=IP1,IG1,2
        KRA=(I4-1)/2
        KRA1=KRA+1
        IF(KRA1 > 30)GO TO 10
        PMGG(KRA1)=ZERO
        CALL RECO4(JA,JB,JC,JD,J(1),J(2),J(3),J(4),KRA*2,0,IAT,RECC)
        IF(IAT == 0) CYCLE
        CALL RECO4(JA,JB,JC,JD,J(1),J(2),J(3),J(4),KRA*2,1,IAT,RECC)
        PMGG(KRA1)=RECC
      END DO
!
!     TRANSFORM FANO & RACAH PHASE CONVENTION
!     TO CONDON & SHORTLEY PHASE CONVENTION
!
      IFAZFRCS=1
      IFAZ=IK1(5)*IK1(4)+IK2(5)*IK2(4)-ID1(5)*ID1(4)-ID2(5)*ID2(4)  &
      +IK3(5)*IK3(4)-ID3(5)*ID3(4)+IK4(5)*IK4(4)-ID4(5)*ID4(4)
      IF((IFAZ/4)*4 /= IFAZ)IFAZFRCS=-IFAZFRCS
!
      NN=0
      JB1=JB-1
      IFAZP=1
      DO II=JA,JB1
        INN=JLIST(II)
        NN=NQ1(INN)+NN
      END DO
      IF((NN/2)*2 == NN)IFAZP=-IFAZP
      NN=0
      JD1=JD-1
      DO II=JC,JD1
        INN=JLIST(II)
        NN=NQ1(INN)+NN
      END DO
      IF((NN/2)*2 == NN)IFAZP=-IFAZP
! * * *                      * * *                      * * *
!     CASES 1324   + + - -        TRANSFORM TO  1234   + - + -
!           1342                                1234
!                                                    (IREZ = 1)
!     OR
!     CASES 2413   + + - -        TRANSFORM TO  1234   - + - +
!           4231                                1234
!                                                    (IREZ = 2)
      DO I3=IP1,IG1,2
        KRA=(I3-1)/2
        IF (ICOLBREI == 1) THEN
          IF(IREZ == 2) THEN
            CALL COULOM(L2,L4,L1,L3,ID2(5),ID4(5),ID1(5),ID3(5),KRA,A1)
          ELSE
            CALL COULOM(L1,L3,L2,L4,ID1(5),ID3(5),ID2(5),ID4(5),KRA,A1)
          END IF
          IF(DABS(A1) < EPS) CYCLE
        END IF
        KRA1=KRA+1
        IF(KRA1 > 30)GO TO 10
        AA=PMGG(KRA1)
        IF(DABS(AA) < EPS) CYCLE
        AA=AA*RAG
        IF(DABS(AA) < EPS) CYCLE
        AA=AA/DSQRT(DBLE(I3))
        AB=AA*DBLE(IFAZP)
        IF(IREZ == 2) THEN
          IFAZ=J(1)+J(2)+J(4)+J(3)+4*KRA
          IF((IFAZ/4)*4 /= IFAZ)AB=-AB
        END IF
        IF (ICOLBREI == 1) THEN
          BB=A1*AB*DBLE(IFAZFRCS)
          IF(IREZ == 1)CALL SPEAK(JJA,JJB,IA,IC,IB,ID,KRA,BB)
          IF(IREZ == 2)CALL SPEAK(JJA,JJB,IB,ID,IA,IC,KRA,BB)
        ELSE IF (ICOLBREI == 2) THEN
          NU=KRA
          IF(((NU-ND1)/2)*2 == (NU-ND1)) THEN
            IF((ITRIG(KS(1),KS(3),NU+NU+1) /= 0) .AND.  &
               (ITRIG(KS(2),KS(4),NU+NU+1) /= 0)) THEN
              N=(NU-ND1)/2+1
              IF(NU >  0) THEN
                CALL CXK(S,IS,KAPS,NU,KRA,1,1)
                DO MU = 1,4
                  COND(MU,N)=COND(MU,N)+AB*S(MU)
                END DO
              END IF
            END IF
          END IF
          NU=KRA+1
          IF(((NU-ND1)/2)*2 == (NU-ND1)) THEN
            IF((ITRIG(KS(1),KS(3),NU+NU-1) /= 0) .AND.  &
               (ITRIG(KS(2),KS(4),NU+NU-1) /= 0)) THEN
              N=(NU-ND1)/2+1
              IF(N <= ND2) THEN
                CALL CXK(S,IS,KAPS,NU,KRA,1,1)
                DO MU = 1,4
                  COND(MU,N)=COND(MU,N)+AB*S(MU)
                END DO
              END IF
            END IF
          END IF
          NU=KRA-1
          IF(((NU-ND1)/2)*2 == (NU-ND1)) THEN
            IF((ITRIG(KS(1),KS(3),NU+NU+3) /= 0) .AND.   &
               (ITRIG(KS(2),KS(4),NU+NU+3) /= 0)) THEN
              IF(NU >= 0) THEN
                N=(NU-ND1)/2+1
                IF(N < ND2) THEN
                  CALL CXK(S,IS,KAPS,NU,KRA,1,1)
                  DO MU = 1,12
                    COND(MU,N)=COND(MU,N)+AB*S(MU)
                  END DO
                END IF
              END IF
            END IF
          END IF
        END IF
      END DO
      IF (ICOLBREI == 2) THEN
        DO N = 1,ND2
          NU=ND1+2*(N-1)
          CALL TALK(JJA,JJB,NU,IS(1),IS(3),IS(2),IS(4),1,COND(1,N))
          CALL TALK(JJA,JJB,NU,IS(3),IS(1),IS(4),IS(2),1,COND(2,N))
          CALL TALK(JJA,JJB,NU,IS(1),IS(3),IS(4),IS(2),1,COND(3,N))
          CALL TALK(JJA,JJB,NU,IS(3),IS(1),IS(2),IS(4),1,COND(4,N))
          IF(N == ND2) CYCLE
          NUP1=NU+1
          CALL TALK(JJA,JJB,NUP1,IS(1),IS(3),IS(2),IS(4),2,COND(5,N))
          CALL TALK(JJA,JJB,NUP1,IS(2),IS(4),IS(1),IS(3),2,COND(6,N))
          CALL TALK(JJA,JJB,NUP1,IS(3),IS(1),IS(4),IS(2),2,COND(7,N))
          CALL TALK(JJA,JJB,NUP1,IS(4),IS(2),IS(3),IS(1),2,COND(8,N))
          CALL TALK(JJA,JJB,NUP1,IS(1),IS(3),IS(4),IS(2),2,COND(9,N))
          CALL TALK(JJA,JJB,NUP1,IS(4),IS(2),IS(1),IS(3),2,COND(10,N))
          CALL TALK(JJA,JJB,NUP1,IS(3),IS(1),IS(2),IS(4),2,COND(11,N))
          CALL TALK(JJA,JJB,NUP1,IS(2),IS(4),IS(3),IS(1),2,COND(12,N))
        END DO
      END IF
! * * *                      * * *                      * * *
!     CASES 1342   + + - -        TRANSFORM TO  1234   + - + -
!           3124                                1234
!                                                    (IREZ = 1)
!     OR
!     CASES 2431   + + - -        TRANSFORM TO  1234   - + - +
!           4213                                1234
!                                                    (IREZ = 2)
      IP2=ITREXG(J(1),J(4),J(2),J(3),IKK)+1
      IF(IKK <= 0) RETURN
      IG2=IP2+IKK-1
      DO I2=IP2,IG2,2
        KRA=(I2-1)/2
!
        IF (ICOLBREI == 1) THEN
          IF(IREZ == 2) THEN
            CALL COULOM(L2,L4,L3,L1,ID2(5),ID4(5),ID3(5),ID1(5),KRA,A1)
          ELSE
            CALL COULOM(L1,L3,L4,L2,ID1(5),ID3(5),ID4(5),ID2(5),KRA,A1)
          END IF
          IF(DABS(A1) < EPS) CYCLE
        END IF
!
        AB=ZERO
        DO I3=IP1,IG1,2
          J12=(I3-1)/2
          IFAZ=J(1)-J(3)+2*J12
          IF(IREZ == 2)IFAZ=J(4)-J(2)-2*J12
          IF((IFAZ/2)*2 /= IFAZ) CYCLE
          KRA1=J12+1
          IF(KRA1 > 30)GO TO 10
          AA=PMGG(KRA1)
          IF(DABS(AA) < EPS) CYCLE
          AA=AA*RAG
          IF(DABS(AA) < EPS) CYCLE
          IF(IXJTIK(J(1),J(4),KRA*2,J(3),J(2),J12*2) == 0) CYCLE
          CALL SIXJ(J(1),J(4),KRA*2,J(3),J(2),J12*2,0,SI)
          AA=AA*SI*DSQRT(DBLE(I3))
          IFAZ=2*J(3)-4*KRA-4*J12
          IF(IREZ == 2)IFAZ=J(1)+J(2)+J(3)-J(4)-4*KRA
          IF((IFAZ/4)*4 /= IFAZ)AA=-AA
          AB=AB+AA
        END DO
        IF(DABS(AB) < EPS) CYCLE
        AB=-AB*DBLE(IFAZP)
        IF (ICOLBREI == 1) THEN
          BB=A1*AB*DBLE(IFAZFRCS)
          IF(IREZ == 1)CALL SPEAK(JJA,JJB,IA,IC,ID,IB,KRA,BB)
          IF(IREZ == 2)CALL SPEAK(JJA,JJB,IB,ID,IC,IA,KRA,BB)
        ELSE IF (ICOLBREI == 2) THEN
          NU=KRA
          IF(((NU-NE1)/2)*2 == (NU-NE1)) THEN
            IF((ITRIG(KS(1),KS(4),NU+NU+1) /= 0) .AND.  &
               (ITRIG(KS(2),KS(3),NU+NU+1) /= 0)) THEN
              N=(NU-NE1)/2+1
              IF(NU > 0) THEN
                CALL CXK(S,IS,KAPS,NU,KRA,1,2)
                DO MU = 1,4
                  CONE(MU,N)=CONE(MU,N)+AB*S(MU)
                END DO
              END IF
            END IF
          END IF
          NU=KRA+1
          IF(((NU-NE1)/2)*2 == (NU-NE1)) THEN
            IF((ITRIG(KS(1),KS(4),NU+NU-1) /= 0) .AND.  &
               (ITRIG(KS(2),KS(3),NU+NU-1) /= 0)) THEN
              N=(NU-NE1)/2+1
              IF(N <= NE2) THEN
                CALL CXK(S,IS,KAPS,NU,KRA,1,2)
                DO MU = 1,4
                  CONE(MU,N)=CONE(MU,N)+AB*S(MU)
                END DO
              END IF
            END IF
          END IF
          NU=KRA-1
          IF(((NU-NE1)/2)*2 == (NU-NE1)) THEN
            IF((ITRIG(KS(1),KS(4),NU+NU+3) /= 0) .AND.  &
               (ITRIG(KS(2),KS(3),NU+NU+3) /= 0)) THEN
              IF(NU >= 0) THEN
                N=(NU-NE1)/2+1
                IF(N < NE2) THEN
                  CALL CXK(S,IS,KAPS,NU,KRA,1,2)
                  DO MU = 1,12
                    CONE(MU,N)=CONE(MU,N)+AB*S(MU)
                  END DO
                END IF
              END IF
            END IF
          END IF
        END IF
      END DO
      IF (ICOLBREI == 2) THEN
        DO N = 1,NE2
          NU=NE1+2*(N-1)
          CALL TALK(JJA,JJB,NU,IS(1),IS(4),IS(2),IS(3),1,CONE(1,N))
          CALL TALK(JJA,JJB,NU,IS(4),IS(1),IS(3),IS(2),1,CONE(2,N))
          CALL TALK(JJA,JJB,NU,IS(1),IS(4),IS(3),IS(2),1,CONE(3,N))
          CALL TALK(JJA,JJB,NU,IS(4),IS(1),IS(2),IS(3),1,CONE(4,N))
          IF(N == NE2) CYCLE
          NUP1=NU+1
          CALL TALK(JJA,JJB,NUP1,IS(1),IS(4),IS(2),IS(3),2,CONE(5,N))
          CALL TALK(JJA,JJB,NUP1,IS(2),IS(3),IS(1),IS(4),2,CONE(6,N))
          CALL TALK(JJA,JJB,NUP1,IS(4),IS(1),IS(3),IS(2),2,CONE(7,N))
          CALL TALK(JJA,JJB,NUP1,IS(3),IS(2),IS(4),IS(1),2,CONE(8,N))
          CALL TALK(JJA,JJB,NUP1,IS(1),IS(4),IS(3),IS(2),2,CONE(9,N))
          CALL TALK(JJA,JJB,NUP1,IS(3),IS(2),IS(1),IS(4),2,CONE(10,N))
          CALL TALK(JJA,JJB,NUP1,IS(4),IS(1),IS(2),IS(3),2,CONE(11,N))
          CALL TALK(JJA,JJB,NUP1,IS(2),IS(3),IS(4),IS(1),2,CONE(12,N))
        END DO
      END IF
      RETURN
   10 WRITE(99,100)
  100 FORMAT(5X,'ERRO IN EL52  PMGG RAGG')
      STOP
      END SUBROUTINE EL52
      !
      !
!*******************************************************************
!                                                                  *
      SUBROUTINE EL53(JJA,JJB,JA,JB,JC,JD,IREZ,ICOLBREI)
!                                                                  *
!   --------------  SECTION METWO    SUBPROGRAM 14  -------------  *
!                                                                  *
!     THIS PACKAGE EVALUATED THE CASES - 1423, 4132, 1432, 4123    *
!                                                   ( IREZ = 1),   *
!                                                   ( + + - - ),   *
!     WHICH APPEARS IN CALCULATION MATRIX ELEMENTS BETWEEN         *
!     CONFIGURATIONS:                               N'1 = N1 - 1   *
!                                                   N'2 = N2 + 1   *
!                                                   N'3 = N3 + 1   *
!                                                   N'4 = N4 - 1   *
!     AND    2314, 3241, 2341, 3214                 ( IREZ = 2),   *
!                                                   ( + + - - ),   *
!     WHICH APPEARS IN CALCULATION MATRIX ELEMENTS BETWEEN         *
!     CONFIGURATIONS:                               N'1 = N1 + 1   *
!                                                   N'2 = N2 - 1   *
!                                                   N'3 = N3 - 1   *
!                                                   N'4 = N4 + 1   *
!                                                                  *
!   Written by  G. Gaigalas                                        *
!   The last modification made by G. Gaigalas       October  2017  *
!                                                                  *
!*******************************************************************
!
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: JJA,JJB,JA,JB,JC,JD,IREZ,ICOLBREI
!      DIMENSION PMGG(30),J(4)
!      DIMENSION COND(12,20),CONE(12,20),S(12),IS(4),KAPS(4),KS(4)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IA,IB,IC,ID,IBRD,IBRE,II,IP1,IP2,IG1,IG2,IKK,I2,I3,  &
                 I4,IFAZ,IFAZP,IFAZFRCS,INN,IAT,KRA,KRA1,L1,L2,L3,L4, &
                 J12,JB1,JD1,ND1,ND2,NE1,NE2,N,NN,NU,NUP1,MU
      INTEGER, DIMENSION(4) :: J
      INTEGER, DIMENSION(4) :: IS,KAPS,KS
      REAL(KIND=8)          :: AA,AB,A1,BB,QM1,QM2,QM3,QM4,RAG,RECC,SI
      REAL(KIND=8), DIMENSION(12) :: S
      REAL(KIND=8), DIMENSION(30) :: PMGG
      REAL(KIND=8), DIMENSION(12,20) :: COND,CONE
!-----------------------------------------------
      IF(NPEEL <= 3)RETURN
      CALL RECO(JA,JD,JC,JB,3,IAT)
      IF(IAT == 0)RETURN
      IA=JLIST(JA)
      IB=JLIST(JB)
      IC=JLIST(JC)
      ID=JLIST(JD)
      IF(IREZ == 2) THEN
        QM1=-HALF
        QM2=HALF
        QM3=HALF
        QM4=-HALF
      ELSE
        QM1=HALF
        QM2=-HALF
        QM3=-HALF
        QM4=HALF
      END IF
      CALL PERKO2(JA,JB,JC,JD,4)
      J(1)=ID1(3)
      J(2)=ID2(3)
      J(3)=ID3(3)
      J(4)=ID4(3)
      L1=(J(1)+1)/2
      L2=(J(2)+1)/2
      L3=(J(3)+1)/2
      L4=(J(4)+1)/2
      IF (ICOLBREI == 2) THEN
        IF(IREZ == 1) THEN
          IS(1)=IA
          IS(2)=ID
          IS(3)=IB
          IS(4)=IC
        ELSE IF(IREZ == 2) THEN
          IS(1)=IB
          IS(2)=IC
          IS(3)=IA
          IS(4)=ID
        END IF
        KAPS(1:4)=2*NAK(IS(1:4))
        KS(1:4)=IABS(KAPS(1:4))
        CALL SNRC(IS,KAPS,KS,ND1,ND2,NE1,NE2,IBRD,IBRE)
        IF(IBRD <= 0 .AND. IBRE <= 0)RETURN
        COND = Zero;  CONE = Zero
      END IF
      CALL GG1234(IK1,IK2,IK3,IK4,BK1,BK2,BK3,BK4,ID1,ID2,  &
      ID3,ID4,BD1,BD2,BD3,BD4,QM1,QM2,QM3,QM4,RAG)
      IF(DABS(RAG) < EPS) RETURN
      IP1=ITREXG(J(1),J(2),J(3),J(4),IKK)+1
      IF(IKK <= 0)RETURN
      IG1=IP1+IKK-1
      DO I4=IP1,IG1,2
        KRA=(I4-1)/2
        KRA1=KRA+1
        IF(KRA1 > 30)GO TO 10
        PMGG(KRA1)=ZERO
        CALL RECO4(JA,JB,JC,JD,J(1),J(2),J(3),J(4),KRA*2,0,IAT,RECC)
        IF(IAT == 0) CYCLE
        CALL RECO4(JA,JB,JC,JD,J(1),J(2),J(3),J(4),KRA*2,1,IAT,RECC)
        PMGG(KRA1)=RECC
      END DO
!
!     TRANSFORM FANO & RACAH PHASE CONVENTION
!     TO CONDON & SHORTLEY PHASE CONVENTION
!
      IFAZFRCS=1
      IFAZ=IK1(5)*IK1(4)+IK2(5)*IK2(4)-ID1(5)*ID1(4)-ID2(5)*ID2(4) &
      +IK3(5)*IK3(4)-ID3(5)*ID3(4)+IK4(5)*IK4(4)-ID4(5)*ID4(4)
      IF((IFAZ/4)*4 /= IFAZ)IFAZFRCS=-IFAZFRCS
!
      NN=0
      JB1=JB-1
      IFAZP=1
      DO II=JA,JB1
        INN=JLIST(II)
        NN=NQ1(INN)+NN
      END DO
      IF((NN/2)*2 == NN)IFAZP=-IFAZP
      NN=0
      JD1=JD-1
      DO II=JC,JD1
        INN=JLIST(II)
        NN=NQ1(INN)+NN
      END DO
      IF((NN/2)*2 == NN)IFAZP=-IFAZP
! * * *                      * * *                      * * *
!     CASES 1423   + + - -        TRANSFORM TO  1234   + - - +
!           4132                                1234
!                                                    (IREZ = 1)
!     OR
!     CASES 2314   + + - -        TRANSFORM TO  1234   - + + -
!           3241                                1234
!                                                    (IREZ = 2)
      DO I3=IP1,IG1,2
        KRA=(I3-1)/2
        IF (ICOLBREI == 1) THEN
          IF(IREZ == 2) THEN
            CALL COULOM(L2,L3,L1,L4,ID2(5),ID3(5),ID1(5),ID4(5),KRA,A1)
          ELSE
            CALL COULOM(L1,L4,L2,L3,ID1(5),ID4(5),ID2(5),ID3(5),KRA,A1)
          END IF
        IF(DABS(A1) < EPS) CYCLE
        END IF
        KRA1=KRA+1
        IF(KRA1 > 30)GO TO 10
        AA=PMGG(KRA1)
        IF(DABS(AA) < EPS) CYCLE
        AA=AA*RAG
        IF(DABS(AA) < EPS) CYCLE
        AA=AA/DSQRT(DBLE(I3))
        IFAZ=J(4)+J(3)-2*KRA+2
        IF(IREZ == 2)IFAZ=J(1)+J(2)-2*KRA+2
        IF((IFAZ/4)*4 /= IFAZ)AA=-AA
        AB=AA*DBLE(IFAZP)
        IF (ICOLBREI == 1) THEN
          BB=A1*AB*DBLE(IFAZFRCS)
          IF(IREZ == 1)CALL SPEAK(JJA,JJB,IA,ID,IB,IC,KRA,BB)
          IF(IREZ == 2)CALL SPEAK(JJA,JJB,IB,IC,IA,ID,KRA,BB)
        ELSE IF (ICOLBREI == 2) THEN
          NU=KRA
          IF(((NU-ND1)/2)*2 == (NU-ND1)) THEN
            IF((ITRIG(KS(1),KS(3),NU+NU+1) /= 0) .AND. &
               (ITRIG(KS(2),KS(4),NU+NU+1) /= 0)) THEN
              N=(NU-ND1)/2+1
              IF(NU > 0) THEN
                CALL CXK(S,IS,KAPS,NU,KRA,1,1)
                DO MU = 1,4
                  COND(MU,N)=COND(MU,N)+AB*S(MU)
                END DO
              END IF
            END IF
          END IF
          NU=KRA+1
          IF(((NU-ND1)/2)*2 == (NU-ND1)) THEN
            IF((ITRIG(KS(1),KS(3),NU+NU-1) /= 0) .AND.  &
               (ITRIG(KS(2),KS(4),NU+NU-1) /= 0)) THEN
              N=(NU-ND1)/2+1
              IF(N <= ND2) THEN
                CALL CXK(S,IS,KAPS,NU,KRA,1,1)
                DO MU = 1,4
                  COND(MU,N)=COND(MU,N)+AB*S(MU)
                END DO
              END IF
            END IF
          END IF
          NU=KRA-1
          IF(((NU-ND1)/2)*2 == (NU-ND1)) THEN
            IF((ITRIG(KS(1),KS(3),NU+NU+3) /= 0) .AND.  &
               (ITRIG(KS(2),KS(4),NU+NU+3) /= 0)) THEN
              IF(NU >= 0) THEN
                N=(NU-ND1)/2+1
                IF(N < ND2) THEN
                  CALL CXK(S,IS,KAPS,NU,KRA,1,1)
                  DO MU = 1,12
                    COND(MU,N)=COND(MU,N)+AB*S(MU)
                  END DO
                END IF
              END IF
            END IF
          END IF
        END IF
      END DO
      IF (ICOLBREI == 2) THEN
        DO N = 1,ND2
          NU=ND1+2*(N-1)
          CALL TALK(JJA,JJB,NU,IS(1),IS(3),IS(2),IS(4),1,COND(1,N))
          CALL TALK(JJA,JJB,NU,IS(3),IS(1),IS(4),IS(2),1,COND(2,N))
          CALL TALK(JJA,JJB,NU,IS(1),IS(3),IS(4),IS(2),1,COND(3,N))
          CALL TALK(JJA,JJB,NU,IS(3),IS(1),IS(2),IS(4),1,COND(4,N))
          IF(N == ND2) CYCLE
          NUP1=NU+1
          CALL TALK(JJA,JJB,NUP1,IS(1),IS(3),IS(2),IS(4),2,COND(5,N))
          CALL TALK(JJA,JJB,NUP1,IS(2),IS(4),IS(1),IS(3),2,COND(6,N))
          CALL TALK(JJA,JJB,NUP1,IS(3),IS(1),IS(4),IS(2),2,COND(7,N))
          CALL TALK(JJA,JJB,NUP1,IS(4),IS(2),IS(3),IS(1),2,COND(8,N))
          CALL TALK(JJA,JJB,NUP1,IS(1),IS(3),IS(4),IS(2),2,COND(9,N))
          CALL TALK(JJA,JJB,NUP1,IS(4),IS(2),IS(1),IS(3),2,COND(10,N))
          CALL TALK(JJA,JJB,NUP1,IS(3),IS(1),IS(2),IS(4),2,COND(11,N))
          CALL TALK(JJA,JJB,NUP1,IS(2),IS(4),IS(3),IS(1),2,COND(12,N))
        END DO
      END IF
! * * *                      * * *                      * * *
!     CASES 1432   + + - -        TRANSFORM TO  1234   + - - +
!           4132                                1234
!                                                    (IREZ = 1)
!     OR
!     CASES 2341   + + - -        TRANSFORM TO  1234   - + + -
!           3214                                1234
!                                                    (IREZ = 2)
      IP2=ITREXG(J(1),J(3),J(4),J(2),IKK)+1
      IF(IKK <= 0) RETURN
      IG2=IP2+IKK-1
      DO I2=IP2,IG2,2
        KRA=(I2-1)/2
        IF (ICOLBREI == 1) THEN
          IF(IREZ == 2) THEN
            CALL COULOM(L2,L3,L4,L1,ID2(5),ID3(5),ID4(5),ID1(5),KRA,A1)
          ELSE
            CALL COULOM(L1,L4,L3,L2,ID1(5),ID4(5),ID3(5),ID2(5),KRA,A1)
          END IF
          IF(DABS(A1) < EPS) CYCLE
        END IF
        AB=ZERO
        DO I3=IP1,IG1,2
          J12=(I3-1)/2
          IFAZ=J(1)+J(4)+2*J12
          IF(IREZ == 2)IFAZ=J(2)+J(3)+2*J12
          IF((IFAZ/2)*2 /= IFAZ) CYCLE
          KRA1=J12+1
          IF(KRA1 > 30)GO TO 10
          AA=PMGG(KRA1)
          IF(DABS(AA) < EPS) CYCLE
          AA=AA*RAG
          IF(DABS(AA) < EPS) CYCLE
          IF(IXJTIK(J(1),J(3),KRA*2,J(4),J(2),J12*2) == 0) CYCLE
          CALL SIXJ(J(1),J(3),KRA*2,J(4),J(2),J12*2,0,SI)
          AA=AA*SI*DSQRT(DBLE(I3))
          IFAZ=J(3)-J(4)-4*KRA+2*J12
          IF(IREZ == 2)IFAZ=J(1)-J(2)-4*KRA-2*J12
          IF((IFAZ/4)*4 /= IFAZ)AA=-AA
          AB=AB+AA
        END DO
        IF(DABS(AB) < EPS) CYCLE
        AB=AB*DBLE(IFAZP)
        IF (ICOLBREI == 1) THEN
          BB=A1*AB*DBLE(IFAZFRCS)
          IF(IREZ == 1)CALL SPEAK(JJA,JJB,IA,ID,IC,IB,KRA,BB)
          IF(IREZ == 2)CALL SPEAK(JJA,JJB,IB,IC,ID,IA,KRA,BB)
        ELSE IF (ICOLBREI == 2) THEN
          NU=KRA
          IF(((NU-NE1)/2)*2 == (NU-NE1)) THEN
            IF((ITRIG(KS(1),KS(4),NU+NU+1) /= 0) .AND.  &
               (ITRIG(KS(2),KS(3),NU+NU+1) /= 0)) THEN
              N=(NU-NE1)/2+1
              IF(NU > 0) THEN
                CALL CXK(S,IS,KAPS,NU,KRA,1,2)
                DO MU = 1,4
                  CONE(MU,N)=CONE(MU,N)+AB*S(MU)
                END DO
              END IF
            END IF
          END IF
          NU=KRA+1
          IF(((NU-NE1)/2)*2 == (NU-NE1)) THEN
            IF((ITRIG(KS(1),KS(4),NU+NU-1) /= 0) .AND.  &
               (ITRIG(KS(2),KS(3),NU+NU-1) /= 0)) THEN
              N=(NU-NE1)/2+1
              IF(N <= NE2) THEN
                CALL CXK(S,IS,KAPS,NU,KRA,1,2)
                DO MU = 1,4
                  CONE(MU,N)=CONE(MU,N)+AB*S(MU)
                END DO
              END IF
            END IF
          END IF
          NU=KRA-1
          IF(((NU-NE1)/2)*2 == (NU-NE1)) THEN
            IF((ITRIG(KS(1),KS(4),NU+NU+3) /= 0) .AND.  &
               (ITRIG(KS(2),KS(3),NU+NU+3) /= 0)) THEN
              IF(NU >= 0) THEN
                N=(NU-NE1)/2+1
                IF(N < NE2) THEN
                  CALL CXK(S,IS,KAPS,NU,KRA,1,2)
                  DO MU = 1,12
                    CONE(MU,N)=CONE(MU,N)+AB*S(MU)
                  END DO
                END IF
              END IF
            END IF
          END IF
        END IF
      END DO
      IF (ICOLBREI == 2) THEN
        DO N = 1,NE2
          NU=NE1+2*(N-1)
          CALL TALK(JJA,JJB,NU,IS(1),IS(4),IS(2),IS(3),1,CONE(1,N))
          CALL TALK(JJA,JJB,NU,IS(4),IS(1),IS(3),IS(2),1,CONE(2,N))
          CALL TALK(JJA,JJB,NU,IS(1),IS(4),IS(3),IS(2),1,CONE(3,N))
          CALL TALK(JJA,JJB,NU,IS(4),IS(1),IS(2),IS(3),1,CONE(4,N))
          IF(N == NE2) CYCLE
          NUP1=NU+1
          CALL TALK(JJA,JJB,NUP1,IS(1),IS(4),IS(2),IS(3),2,CONE(5,N))
          CALL TALK(JJA,JJB,NUP1,IS(2),IS(3),IS(1),IS(4),2,CONE(6,N))
          CALL TALK(JJA,JJB,NUP1,IS(4),IS(1),IS(3),IS(2),2,CONE(7,N))
          CALL TALK(JJA,JJB,NUP1,IS(3),IS(2),IS(4),IS(1),2,CONE(8,N))
          CALL TALK(JJA,JJB,NUP1,IS(1),IS(4),IS(3),IS(2),2,CONE(9,N))
          CALL TALK(JJA,JJB,NUP1,IS(3),IS(2),IS(1),IS(4),2,CONE(10,N))
          CALL TALK(JJA,JJB,NUP1,IS(4),IS(1),IS(2),IS(3),2,CONE(11,N))
          CALL TALK(JJA,JJB,NUP1,IS(2),IS(3),IS(4),IS(1),2,CONE(12,N))
        END DO
      END IF
      RETURN
   10 WRITE(99,100)KRA1
  100 FORMAT(5X,'ERRO IN EL53  PMGG RAGG KRA1=',I100)
      STOP
      END SUBROUTINE EL53
    
  END SUBROUTINE RKCO

  


!*******************************************************************
!                                                                  *
      SUBROUTINE COULOM(J1,J2,J3,J4,L1,L2,L3,L4,K,AA)
!                                                                  *
!   --------------  SECTION METWO    SUBPROGRAM 01  ------------   *
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF COULOMB INTERACTIONS BETWEEN THE ELECTRONS                *
!                                                                  *
!                          k   k+1  (k) (k)                        *
!     (n l j T  n l j T ::r  / r  ( C   C )::n l j T  n l j T )    *
!       1 1 1 1  2 2 2 2   <    >             3 3 3 3  4 4 4 4     *
!                                                                  *
!                                                                  *
!   Written by G. Gaigalas,                                        *
!   The last modification made by G. Gaigalas       October  2017  *
!                                                                  *
!*******************************************************************
!
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)        :: J1,J2,J3,J4,L1,L2,L3,L4,K
      REAL(kind=8), INTENT(OUT)  :: AA
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, IFAZ, ITTK
!-----------------------------------------------
      AA=ZERO
      IF(ITTK(L1,L3,K) == 0)RETURN
      IF(ITTK(L2,L4,K) == 0)RETURN
      I=(2*K+1)/2
      AA=CRE (J1,I,J3)
      IF(DABS(AA) < EPS)RETURN
      AA=AA*CRE (J2,I,J4)
      IF(DABS(AA) < EPS)RETURN
      IFAZ=L3-2*K-L1+L4-L2
      IF((IFAZ/4)*4 /= IFAZ)AA=-AA
      RETURN

    CONTAINS

!***********************************************************************
!                                                                      *
      REAL(KIND(0.0D0)) FUNCTION CRE (KAP1, K, KAP2) 
!-----------------------------------------------
!                                                                      *
!   Computes the relativistic reduced matrix element                   *
!                                                                      *
!                         (j1 || C(K) || j2),                          *
!                                                                      *
!   Eq. (5.15) of I P Grant, Advances in Physics 19 (1970) 762. KAP1,  *
!   KAP2 are the kappa values corresponding to j1, j2.  The triangle   *
!   conditions are tested by the routine CLRX.                         *
!                                                                      *
!                                           Last update: 09 Oct 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:47:10   2/14/04  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: KAP1 
      INTEGER  :: K 
      INTEGER  :: KAP2 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K1 
      REAL(kind=8) :: DK1K2 
!-----------------------------------------------
!
      K1 = ABS(KAP1) 
      DK1K2 = DBLE(4*K1*IABS(KAP2)) 
      CRE = SQRT(DK1K2)*CLRX(INT(KAP1,kind=1),K,INT(KAP2,kind=1)) 
      IF (MOD(K1,2) == 1) CRE = -CRE 
!
      RETURN  
      END FUNCTION CRE 


    END SUBROUTINE COULOM


!***********************************************************************
!                                                                      *
      SUBROUTINE CXK (S,IS,KAPS,NU,K,IBR,IEX)
!                                                                      *
!***********************************************************************
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
      IMPLICIT NONE
!      DIMENSION IS(4),KAPS(4),S(12)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: NU
      INTEGER :: K
      INTEGER, INTENT(IN) :: IBR
      INTEGER, INTENT(IN) :: IEX
      INTEGER, INTENT(IN) :: IS(4)
      INTEGER, INTENT(IN) :: KAPS(4)
      REAL(kind=8) , INTENT(INOUT) :: S(12)
!-----------------------------------------------
!
      STOP 'CXK: Error '
    END SUBROUTINE CXK

  SUBROUTINE DIAGA1(JA1,KA,IRE,IAT,RECC)
!*******************************************************************
!                                                                  *
!     SUBROUTINE CALLED:  IXJTIK, SIXJ                             *
!                                                                  *
!                                                                  *
!*******************************************************************
    IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
    INTEGER, INTENT(IN)       :: JA1, KA, IRE
    INTEGER, INTENT(OUT)      :: IAT
    REAL(KIND=8), INTENT(OUT) :: RECC
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
    INTEGER :: IJ1,IA1,IB1,K1,J1,IT1,IT1S,IFAZ,LL1
    REAL(KIND=8) :: A1
!-----------------------------------------------
    RECC = ZERO
    IJ1=JLIST(JA1)
    IA1=JJQ1(3,IJ1)-1
    IB1=JJQ2(3,IJ1)-1
    IF(JA1 <= 2) THEN
       LL1=JLIST(1)
       IF(JA1 == 1)LL1=JLIST(2)
       J1=JJQ1(3,LL1)-1
       IT1=JJC1(1)-1
       IT1S=JJC2(1)-1
    ELSE
       K1=JA1-2
       J1=JJC1(K1)-1
       K1=JA1-1
       IT1=JJC1(K1)-1
       IT1S=JJC2(K1)-1
    END IF
    IF(IRE /= 0) THEN
       CALL SIXJ(KA,IB1,IA1,J1,IT1,IT1S,0,A1)
       A1=A1*DSQRT(DBLE((IA1+1)*(IT1S+1)))
       IFAZ=J1+IT1+IB1+KA
       IF((IFAZ/4)*4 /= IFAZ)A1=-A1
       RECC=A1
       IAT=1
       IF(JA1 /= 1)RETURN
       IFAZ=IA1+IB1+2*J1-IT1-IT1S
       IF((IFAZ/4)*4 /= IFAZ)RECC=-RECC
    ELSE
       IF(IXJTIK(KA,IB1,IA1,J1,IT1,IT1S) == 0)RETURN
       IAT=1
    END IF
    RETURN
  END SUBROUTINE DIAGA1
!
!###################################################################
!
  SUBROUTINE DIAGA2(JA1,JA2,KA,IRE,IAT,RECC)
!*******************************************************************
!                                                                  *
!     SUBROUTINE CALLED:  IXJTIK, SIXJ                             *
!                                                                  *
!                                                                  *
!*******************************************************************
    IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
    INTEGER, INTENT(IN)       :: JA1, JA2, KA, IRE
    INTEGER, INTENT(OUT)      :: IAT
    REAL(KIND=8), INTENT(OUT) :: RECC
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
    INTEGER :: IJ1,IJ2,IA1,IA2,IB1,IB2,N1,N2,J2,IT2,IT2S,IFAZ
    REAL(KIND=8) :: A2
!-----------------------------------------------
    RECC = ZERO
    IJ1=JLIST(JA1)
    IJ2=JLIST(JA2)
    IA1=JJQ1(3,IJ1)-1
    IA2=JJQ1(3,IJ2)-1
    IB1=JJQ2(3,IJ1)-1
    IB2=JJQ2(3,IJ2)-1
    IF(JA1 == 1.AND.JA2 == 2) THEN
       IT2=IA1
       IT2S=IB1
       J2=JJC1(1)-1
    ELSE
       N1=JA2-1
       J2=JJC1(N1)-1
       N2=JA2-2
       IT2=JJC1(N2)-1
       IT2S=JJC2(N2)-1
    END IF
    IF(IRE /= 0) THEN
       CALL SIXJ(KA,IB2,IA2,J2,IT2,IT2S,0,A2)
       RECC=A2*DSQRT(DBLE((IB2+1)*(IT2+1)))
       IFAZ=J2+IT2S+IA2+KA
       IF((IFAZ/4)*4 /= IFAZ)RECC=-RECC
       IAT=1
    ELSE
       IF(IXJTIK(KA,IB2,IA2,J2,IT2,IT2S) == 0)RETURN
       IAT=1
    END IF
    RETURN
  END SUBROUTINE DIAGA2
!
!###################################################################
!
  SUBROUTINE DIAGA3(JA1,JA2,KA,IRE,IAT,REC)
!*******************************************************************
!                                                                  *
!     SUBROUTINE CALLED:  IXJTIK, SIXJ                             *
!                                                                  *
!                                                                  *
!*******************************************************************
    IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
    INTEGER, INTENT(IN)  :: JA1, JA2, KA, IRE
    INTEGER, INTENT(OUT) :: IAT
    REAL(KIND=8), INTENT(OUT) :: REC
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
    INTEGER :: I,LL1,JI,KK1,KK2,ITI,ITI1,ITIS,ITI1S,IFAZ
    REAL(KIND=8) :: AA, A3
!-----------------------------------------------
    REC = ZERO
    AA=ONE
    I=JA1+1
    IF(JA1 == 1)I=I+1
    IF(I >= JA2) THEN
       REC=AA
       IAT=1
    ELSE
1      LL1=JLIST(I)
       JI=JJQ1(3,LL1)-1
       KK2=I-2
       ITI=JJC1(KK2)-1
       ITIS=JJC2(KK2)-1
       KK1=I-1
       ITI1=JJC1(KK1)-1
       ITI1S=JJC2(KK1)-1
       IF(IRE /= 0) THEN
          CALL SIXJ(KA,ITIS,ITI,JI,ITI1,ITI1S,0,A3)
          A3=A3*SQRT(DBLE((ITI+1)*(ITI1S+1)))
          IFAZ=KA+JI+ITI+ITI1S
          IF((IFAZ/4)*4 /= IFAZ)A3=-A3
          AA=AA*A3
       ELSE
          IF(IXJTIK(KA,ITIS,ITI,JI,ITI1,ITI1S) == 0)RETURN
       END IF
       I=I+1
       IF(I == JA2) THEN
          REC=AA
          IAT=1
       ELSE
          GOTO 1
       END IF
    END IF
    RETURN
  END SUBROUTINE DIAGA3
! 
!###################################################################
!
  SUBROUTINE DIAGA4(JA1,JA2,K1,K2,KA,IRE,IAT,RECC)
!*******************************************************************
!                                                                  *
!     SUBROUTINE CALLED:  NINEJ                                     *
!                                                                  *
!*******************************************************************
    IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
    INTEGER, INTENT(IN)       :: JA1,JA2,K1,K2,KA,IRE
    INTEGER, INTENT(OUT)      :: IAT
    REAL(KIND=8), INTENT(OUT) :: RECC
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
    INTEGER :: IJ1,IJ2,IA1,IA2,IB1,IB2,N1,N2,J2,J2S,IT2,IT2S
    REAL(KIND=8) :: A2
!-----------------------------------------------
    RECC = ZERO
    IJ1=JLIST(JA1)
    IJ2=JLIST(JA2)
    IA1=JJQ1(3,IJ1)-1
    IA2=JJQ1(3,IJ2)-1
    IB1=JJQ2(3,IJ1)-1
    IB2=JJQ2(3,IJ2)-1
    IF(JA1 == 1.AND.JA2 == 2) THEN
       IT2=IA1
       IT2S=IB1
       J2=JJC1(1)-1
       J2S=JJC2(1)-1
    ELSE
       N1=JA2-1
       J2=JJC1(N1)-1
       J2S=JJC2(N1)-1
       N2=JA2-2
       IT2=JJC1(N2)-1
       IT2S=JJC2(N2)-1
    END IF
    IF(IRE /= 0) THEN
       CALL NINEJ(IT2S,K1,IT2,IB2,K2,IA2,J2S,KA,J2,0,IAT,A2)
       RECC=A2*DSQRT(DBLE((IT2+1)*(KA+1)*(IA2+1)*(J2S+1)))
    ELSE
       CALL NINEJ(IT2S,K1,IT2,IB2,K2,IA2,J2S,KA,J2,1,IAT,A2)
    END IF
    RETURN
  END SUBROUTINE DIAGA4
!
!###################################################################
!
  SUBROUTINE DIAGA5(NPEELGG,JA1,KA,IRE,IAT,RECC)
!*******************************************************************
!                                                                  *
!     SUBROUTINE CALLED:  IXJTIK, SIXJ                             *
!                                                                  *
!                                                                  *
!*******************************************************************
    IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
    INTEGER, INTENT(IN)       :: NPEELGG,JA1,KA,IRE
    INTEGER, INTENT(OUT)      :: IAT
    REAL(KIND=8), INTENT(OUT) :: RECC
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
    INTEGER      :: ITI1,ITI1S,IJ1,ITI,ITIS,JI
    REAL(KIND=8) :: A3
!-----------------------------------------------
    RECC = ZERO
    ITI1=JJC1(NPEELGG-1)-1
    ITI1S=JJC2(NPEELGG-1)-1
    IJ1=JLIST(NPEELGG)
    IF(JA1 == NPEELGG) THEN
       ITI=JJQ1(3,IJ1)-1
       ITIS=JJQ2(3,IJ1)-1
       JI=JJC1(NPEELGG-2)-1
    ELSE
       JI=JJQ1(3,IJ1)-1
       ITI=JJC1(NPEELGG-2)-1
       ITIS=JJC2(NPEELGG-2)-1
    END IF
    IF(IRE == 0) THEN
       IF(IXJTIK(KA,ITIS,ITI,JI,ITI1,ITI1S) /= 0) IAT=1
    ELSE
       CALL SIXJ(KA,ITIS,ITI,JI,ITI1,ITI1S,0,A3)
       RECC=A3*DSQRT(DBLE((ITI+1)*(ITI1S+1)))
       IF(MOD(KA+JI+ITIS+ITI1,4) /= 0)RECC=-RECC
       IAT=1
       IF(JA1 == NPEELGG)RETURN
       IF(MOD(ITI+ITIS-ITI1S-ITI1+2*JI,4) /= 0)RECC=-RECC
    END IF
    RETURN
  END SUBROUTINE DIAGA5
  

!*******************************************************************
!                                                                  *
      SUBROUTINE EILE(JA,JB,JC,JAA,JBB,JCC)
!                                                                  *
!     ------------  SECTION METWO    SUBPROGRAM 02  ------------   *
!                                                                  *
!                                                                  *
!   Written by G. Gaigalas,                                        *
!   The last modification made by G. Gaigalas       October  2017  *
!                                                                  *
!*******************************************************************
!
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)  :: JA, JB, JC 
      INTEGER, INTENT(OUT) :: JAA, JBB, JCC
!-----------------------------------------------
      JAA=JA
      JCC=JA
      IF(JAA > JB)JAA=JB
      IF(JCC < JB)JCC=JB
      IF(JAA > JC)JAA=JC
      IF(JCC < JC)JCC=JC
      IF((JA > JAA).AND.(JA < JCC))JBB=JA
      IF((JB > JAA).AND.(JB < JCC))JBB=JB
      IF((JC > JAA).AND.(JC < JCC))JBB=JC
      RETURN
      END SUBROUTINE EILE

      

!***********************************************************************
!                                                                      *
      INTEGER FUNCTION ICHKQ1 (JA, JB) 
!                                                                      *
!   This routine is to check the occupation condition for one electron *
!   operator.                                                          *
!                                                                      *
!   Yu Zou                                Last revision: 8/16/00       *
!                                                                      *
!***********************************************************************
        IMPLICIT NONE        
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER :: JA 
      INTEGER :: JB 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K, I, I_IQA, I_IQB 
!-----------------------------------------------
!
!
      ICHKQ1 = 0 
      K = 0 
      DO I = 1, NW 
         I_IQA = IQA(I,JA) 
         I_IQB = IQA(I,JB) 
         IF (I_IQA == I_IQB) CYCLE  
         K = K + 1 
         IF (K > 2) RETURN  
         IF (IABS(I_IQA - I_IQB) <= 1) CYCLE  
         RETURN  
      END DO 
      IF (K==2 .OR. K==0) ICHKQ1 = 1 
      RETURN  
      END FUNCTION ICHKQ1 


!***********************************************************************
!                                                                      *
      INTEGER FUNCTION ICHKQ2 (JA, JB) 
!                                                                      *
!   This routine is to check the occupation condition for two electron *
!   operator.                                                          *
!                                                                      *
!   Yu Zou                                Last revision: 8/21/00       *
!                                                                      *
!***********************************************************************
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER :: JA 
      INTEGER :: JB 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K, I, I_IQA, I_IQB 
!-----------------------------------------------
!
!
      ICHKQ2 = 0 
      K = 0 
      DO I = 1, NW 
         I_IQA = IQA(I,JA) 
         I_IQB = IQA(I,JB) 
         IF (I_IQA == I_IQB) CYCLE  
         K = K + 1 
         IF (K > 4) RETURN  
         IF (IABS(I_IQA - I_IQB) <= 2) CYCLE  
         RETURN  
      END DO 
      ICHKQ2 = 1 
      RETURN  
      END FUNCTION ICHKQ2 



!*******************************************************************
!                                                                  *
      INTEGER FUNCTION ITREXG(I1,I2,I3,I4,K)
!
!   Written by G. Gaigalas,                                        * 
!   Vanderbilt University,  Nashville               October  1996  *
!   The last modification made by G. Gaigalas       October  2017  *
!                                                                  *
!*******************************************************************
!
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)  :: I1, I2, I3, I4
      INTEGER, INTENT(OUT) :: K
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J
!-----------------------------------------------
      J=MAX0(IABS(I1-I2),IABS(I3-I4))
      K=MIN0(IABS(I1+I2),IABS(I3+I4))-J+1
      ITREXG=J
      RETURN
      END FUNCTION ITREXG

      
!***********************************************************************
!                                                                      *
      INTEGER FUNCTION ITRIG (I1, I2, I3) 
!                                                                      *
!   The  triangular delta. Input: Values of 2*J+1; Output: 1, IF J'S   *
!   form a triangle; 0, otherwise.                                     *
!                                           Last update: 09 Oct 1992   *
!                                                                      *
!***********************************************************************
        IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: I1 
      INTEGER, INTENT(IN) :: I2 
      INTEGER, INTENT(IN) :: I3 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I4 
!-----------------------------------------------
!
      I4 = I2 - I3 
      IF (I1>=ABS(I4) + 1 .AND. I1<=I2+I3-1) THEN 
         ITRIG = 1 
      ELSE 
         ITRIG = 0 
      ENDIF 
!
      RETURN  
      END FUNCTION ITRIG 



      


!*******************************************************************
!                                                                  *
      INTEGER FUNCTION JFAZE(I1,I2,I3,I4)
!                                                                  *
!     ------------  SECTION METWO    SUBPROGRAM 20  -------------  *
!                                                                  *
!     DETERMINATE THE PHASE FACTOR WHICH APPEAR FROM PERMUTATION   *
!     OF OPERATORS OF SECOND QUANTIZATION                          *
!                                                                  *
!                                                                  *
!   Written by  G. Gaigalas                                        *
!   The last modification made by G. Gaigalas       October  2017  *
!                                                                  *
!*******************************************************************
!
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: I1, I2, I3, I4
!-----------------------------------------------
      JFAZE=1
      IF(I1 > I2)JFAZE=-JFAZE
      IF(I1 > I3)JFAZE=-JFAZE
      IF(I1 > I4)JFAZE=-JFAZE
      IF(I2 > I3)JFAZE=-JFAZE
      IF(I2 > I4)JFAZE=-JFAZE
      IF(I3 > I4)JFAZE=-JFAZE
      RETURN
      END FUNCTION JFAZE

      

!*******************************************************************
!                                                                  *
      SUBROUTINE PERKO2(JA1,JA2,JA3,JA4,I)
!                                                                  *
!   --------------  SECTION METWO    SUBPROGRAM 23  -------------  *
!                                                                  *
!     INTERFACE BETWEEN "GRASP" AND BOLCK "SQ"                     *
!                                               (GENERAL CASE)     *
!                                                                  *
!   Written by  G. Gaigalas                                        *
!   The last modification made by G. Gaigalas       October  2017  *
!                                                                  *
!*******************************************************************

        IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: JA1, JA2, JA3, JA4, I
!-----------------------------------------------
      CALL PERKO1(JA1,BK1,IK1,BD1,ID1)
      IF(I == 1)RETURN
      CALL PERKO1(JA2,BK2,IK2,BD2,ID2)
      IF(I == 2)RETURN
      CALL PERKO1(JA3,BK3,IK3,BD3,ID3)
      IF(I == 3)RETURN
      CALL PERKO1(JA4,BK4,IK4,BD4,ID4)
      RETURN

    CONTAINS

      !*******************************************************************
!                                                                  *
      SUBROUTINE PERKO1(JA,BK,IK,BD,ID)
!                                                                  *
!     ------------  SECTION METWO    SUBPROGRAM 22  -------------  *
!                                                                  *
!     INTERFACE BETWEEN "GRASP" AND BOLCK "SQ"                     *
!                                               (FOR ONE SHELL)    *
!                                                                  *
!                                                                  *
!   Written by  G. Gaigalas                                        *
!   The last modification made by G. Gaigalas       October  2017  *
!                                                                  *
!*******************************************************************

        IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER,      INTENT(IN)                :: JA
      INTEGER,      INTENT(OUT), DIMENSION(7) :: IK, ID
      REAL(kind=8), INTENT(OUT), DIMENSION(3) :: BK, BD
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IJ
!-----------------------------------------------
      IJ=JLIST(JA)
      IK(2)=NP(IJ)
      ID(2)=IK(2)
      IK(3)=(IABS(int(NAK(IJ)))*2)-1
      ID(3)=IK(3)
      IK(4)=NQ1(IJ)
      ID(4)=NQ2(IJ)
      IK(5)=(IK(3)+int(NAK(IJ))/IABS(int(NAK(IJ))))/2
      ID(5)=IK(5)
      IK(6)=JJQ1(3,IJ)-1
      ID(6)=JJQ2(3,IJ)-1
      IK(7)=IABS(int(NAK(IJ)))-JJQ1(1,IJ)
      ID(7)=IABS(int(NAK(IJ)))-JJQ2(1,IJ)
      BK(1)=HALF*DBLE(IK(7))
      BD(1)=HALF*DBLE(ID(7))
      BK(2)=HALF*DBLE(IK(6))
      BD(2)=HALF*DBLE(ID(6))
      BK(3)=-HALF*DBLE(IABS(int(NAK(IJ)))-IK(4))
      BD(3)=-HALF*DBLE(IABS(int(NAK(IJ)))-ID(4))
      IK(1)=NMTEJJ(IK(7),IK(6),IK(3),ID(4),IK(4))
      ID(1)=NMTEJJ(ID(7),ID(6),ID(3),ID(4),IK(4))
      RETURN
      END SUBROUTINE PERKO1

    END SUBROUTINE PERKO2


    !*******************************************************************
!                                                                  *
      SUBROUTINE RECO(JA1,JA2,JA3,JA4,KA,IAT)
!                                                                  *
!     -------------  SECTION REC    SUBPROGRAM 05  --------------  *
!                                                                  *
!                                                                  *
!   Written by  G. Gaigalas                                        *
!   The last modification made by G. Gaigalas       October  2017  *
!                                                                  *
!*******************************************************************
!
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)  :: JA1, JA2, JA3, JA4, KA
      INTEGER, INTENT(OUT) :: IAT
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, IJ1, IJ2, IJ, J, IA1, IA2
!-----------------------------------------------
      IAT=1
      IF(NPEEL == 1)RETURN
      IJ1=JLIST(JA1)
      IJ2=JLIST(JA2)
      IF(JA1 == 1.AND.JA2 == 2) GO TO 1
      IF(KA /= 0)GO TO 5
!
!  CASES WHEN :          KA = 0
!                  OR    JA1 = JA2
!                  OR    JA1 = 1    JA2 = 2
!
    1 DO I=1,NPEEL
        IJ=JLIST(I)
        IF(I < NPEEL-1) THEN
          IF(JJC1(I) /= JJC2(I))IAT=0
        END IF
        IF(KA /= 0) THEN
          IF(I == JA1) CYCLE
          IF(I == JA2) CYCLE
        END IF
        DO J=1,3
          IF(JJQ1(J,IJ) /= JJQ2(J,IJ))IAT=0
        END DO
      END DO
      RETURN
!
!  OTHER CASES
!
    5 CONTINUE
      DO I=1,NPEEL
        IJ=JLIST(I)
        IF(I < NPEEL-1) THEN
          IA1=JA1-1
          IA2=JA2-1
          IF(JA1 == 1)IA1=JA1
          IF(I >= IA1.AND.I < IA2)GO TO 7
          IF(JJC1(I) /= JJC2(I))IAT=0
        END IF
    7   IF(I == JA1) CYCLE
        IF(I == JA2) CYCLE
        IF((KA == 2).AND.(I == JA3)) CYCLE
        IF((KA == 3).AND.(I == JA3)) CYCLE
        IF((KA == 3).AND.(I == JA4)) CYCLE
        DO J=1,3
          IF(JJQ1(J,IJ) /= JJQ2(J,IJ))IAT=0
        END DO
      END DO
      RETURN
      END SUBROUTINE RECO

      
  SUBROUTINE RECO2(JA1,JA2,KA,IRE,IAT,RECC)
!*******************************************************************
!                                                                  *
!     SUBROUTINE CALLED:  DIAGA1,DIAGA2,DIAGA3                     *
!                                                                  *
!   Written by  G. Gaigalas                                        *
!   The last modification made by G. Gaigalas       October  2017  *
!                                                                  *
!*******************************************************************
!
    IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
    INTEGER, INTENT(IN)       :: JA1,JA2,KA,IRE
    INTEGER, INTENT(OUT)      :: IAT
    REAL(KIND=8), INTENT(OUT) :: RECC
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
    INTEGER      :: IA1,IA2,IB1,IB2,IJ1,IJ2,ISKR
    REAL(KIND=8) :: S, SS, RE
!-----------------------------------------------
    IAT=0
    IJ1=JLIST(JA1)
    IJ2=JLIST(JA2)
    S=DBLE(JJQ1(3,IJ1))
    SS=DBLE(JJQ1(3,IJ2))
    SS=S*SS
    RECC=ONE/DSQRT(SS)
    IF(IRE == 0) THEN
       IAT=0
    ELSE IF(KA /= 0) THEN
       IAT=0
    ELSE
       IAT=1
       RETURN
    END IF
    IA1=JJQ1(3,IJ1)-1
    IA2=JJQ1(3,IJ2)-1
    IB1=JJQ2(3,IJ1)-1
    IB2=JJQ2(3,IJ2)-1
!
    CALL DIAGA2(JA1,JA2,KA,IRE,IAT,RE)
    IF(IAT == 0)RETURN
    RECC=RE*RECC*DSQRT(DBLE(IA2+1))/DSQRT(DBLE((KA+1)*(IB2+1)))
    IF(JA1 == 1.AND.JA2 == 2)RETURN
!
    IAT=0
    CALL DIAGA1(JA1,KA,IRE,IAT,RE)
    IF(IAT == 0)RETURN
    RECC=RE*RECC
    ISKR=JA2-JA1
    IF(JA1 == 1)ISKR=JA2-1-JA1
    IF(ISKR <= 1)RETURN
!
    IAT=0
    CALL DIAGA3(JA1,JA2,KA,IRE,IAT,RE)
    RECC=RE*RECC
    RETURN
  END SUBROUTINE RECO2
  !
  !##############################################################
  !
  SUBROUTINE REC3(JA1,JA2,JA3,K1,K2,K3,IRE,IAT,RECC)
!*******************************************************************
!                                                                  *
!                                                                  *
!   Written by  G. Gaigalas                                        *
!   The last modification made by G. Gaigalas       October  2017  *
!   Restructured by A. Senchuk                     September 2019  *  
!                                                                  *
!*******************************************************************
!
    IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
    INTEGER, INTENT(IN)       :: JA1,JA2,JA3,K1,K2,K3,IRE
    INTEGER, INTENT(OUT)      :: IAT
    REAL(KIND=8), INTENT(OUT) :: RECC
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
    INTEGER :: IFAZ
!-----------------------------------------------
    IF((JA3 > JA1).AND.(JA3 > JA2)) THEN
       IF(JA1-JA2 < 0) THEN
          CALL RECO3(JA1,JA2,JA3,K1,K2,K3,IRE,IAT,RECC)
       ELSE IF(JA1-JA2 > 0) THEN
          CALL RECO3(JA2,JA1,JA3,K2,K1,K3,IRE,IAT,RECC)
          IFAZ=K1+K2-K3
          IF((IFAZ/4)*4 /= IFAZ)RECC=-RECC
       ELSE
          GO TO 10
       END IF
    ELSE IF((JA3 < JA1).AND.(JA3 < JA2)) THEN
       IF(JA1-JA2 < 0) THEN
          CALL RECO3(JA3,JA1,JA2,K3,K1,K2,IRE,IAT,RECC)
          IF((K3/2)*2 /= K3)RECC=-RECC
       ELSE IF(JA1-JA2 > 0) THEN
          CALL RECO3(JA3,JA2,JA1,K3,K2,K1,IRE,IAT,RECC)
          IFAZ=K1+K2+K3
          IF((IFAZ/4)*4 /= IFAZ)RECC=-RECC
       ELSE
          GO TO 10
       END IF
    ELSE
       IF(JA1-JA2 < 0)THEN
          CALL RECO3(JA1,JA3,JA2,K1,K3,K2,IRE,IAT,RECC)
          IFAZ=K1-K2-K3
          IF((IFAZ/4)*4 /= IFAZ)RECC=-RECC
       ELSE IF(JA1-JA2 > 0) THEN
          CALL RECO3(JA2,JA3,JA1,K2,K3,K1,IRE,IAT,RECC)
          IF((K1/2)*2 /= K1)RECC=-RECC
       ELSE
          GO TO 10
       END IF
    END IF
    RETURN
10  WRITE(99,100)
100 FORMAT(5X,'ERRO IN REC')
    STOP
    
  CONTAINS

    SUBROUTINE RECO3(JA1,JA2,JA3,K1,K2,KA,IRE,IAT,RECC)
!*******************************************************************
!                                                                  *
!     SUBROUTINE CALLED:  DIAGA1,DIAGA2,DIAGA3,DIAGA4              *
!                                                                  *
!*******************************************************************
!
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)       :: JA1,JA2,JA3,K1,K2,KA,IRE
      INTEGER, INTENT(OUT)      :: IAT
      REAL(KIND=8), INTENT(OUT) :: RECC
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER      :: IA3, IB3, IJ1, IJ2, IJ3, ISKR
      REAL(KIND=8) :: S, S1, S2, S3, RE
!-----------------------------------------------
      IJ1=JLIST(JA1)
      IJ2=JLIST(JA2)
      IJ3=JLIST(JA3)
      S1=JJQ1(3,IJ1)
      S2=JJQ1(3,IJ2)
      S3=JJQ1(3,IJ3)
      S=S1*S2*S3
      RECC=ONE/DSQRT(S)
      IA3=JJQ1(3,IJ3)-1
      IB3=JJQ2(3,IJ3)-1
      RECC=RECC*DSQRT(DBLE(IA3+1))/DSQRT(DBLE((KA+1)*(IB3+1)))
!
      IAT=0
      ISKR=JA3-JA2
      IF(ISKR > 1) THEN
         CALL DIAGA3(JA2,JA3,KA,IRE,IAT,RE)
         IF(IAT == 0)RETURN
         RECC=RE*RECC
      END IF
!
      IAT=0
      CALL DIAGA2(JA1,JA3,KA,IRE,IAT,RE)
      IF(IAT == 0)RETURN
      RECC=RE*RECC
!
      IAT=0
      CALL DIAGA4(JA1,JA2,K1,K2,KA,IRE,IAT,RE)
      IF(IAT == 0)RETURN
      RECC=RE*RECC
      IF(JA1 == 1.AND.JA2 == 2)RETURN
!
      IAT=0
      CALL DIAGA1(JA1,K1,IRE,IAT,RE)
      IF(IAT == 0)RETURN
      RECC=RE*RECC
!
      ISKR=JA2-JA1
      IF(JA1 == 1)ISKR=JA2-1-JA1
      IF(ISKR <= 1)RETURN
      IAT=0
      CALL DIAGA3(JA1,JA2,K1,IRE,IAT,RE)
      RECC=RE*RECC
      RETURN
    END SUBROUTINE RECO3
      
  END SUBROUTINE REC3
  !
  !######################################################################
  !
  SUBROUTINE RECO4(JA1,JA2,JA3,JA4,K1,K2,K3,K4,KA,IRE,IAT,RECC)
!                                                                  *
!   ---------------  SECTION REC    SUBPROGRAM 09  --------------  *
!                                                                  *
!                                                                  *
!   Written by  G. Gaigalas                                        *
!   The last modification made by G. Gaigalas       October  2017  *
!                                                                  *
 !*******************************************************************
!
    IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
    INTEGER, INTENT(IN)  :: JA1,JA2,JA3,JA4,K1,K2,K3,K4,KA,IRE
    INTEGER, INTENT(OUT) :: IAT
    REAL(KIND=8), INTENT(OUT) :: RECC
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
    INTEGER      :: IA4, IB4, IJ1, IJ2, IJ3, IJ4, ISKR
    REAL(KIND=8) :: S, S1, S2, S3, S4, RE
!-----------------------------------------------
    IJ1=JLIST(JA1)
    IJ2=JLIST(JA2)
    IJ3=JLIST(JA3)
    IJ4=JLIST(JA4)
    S1=JJQ1(3,IJ1)
    S2=JJQ1(3,IJ2)
    S3=JJQ1(3,IJ3)
    S4=JJQ1(3,IJ4)
    S=S1*S2*S3*S4
    RECC=ONE/DSQRT(S)
    IA4=JJQ1(3,IJ4)-1
    IB4=JJQ2(3,IJ4)-1
    RECC=RECC*DSQRT(DBLE(IA4+1))/DSQRT(DBLE((K4+1)*(IB4+1)))
!
    ISKR=JA3-JA2
    IF(ISKR > 1) THEN
       IAT=0
       CALL DIAGA3(JA2,JA3,KA,IRE,IAT,RE)
       IF(IAT == 0)RETURN
       RECC=RE*RECC
    END IF
!
    ISKR=JA4-JA3
    IF(ISKR > 1) THEN
       IAT=0
       CALL DIAGA3(JA3,JA4,K4,IRE,IAT,RE)
       IF(IAT == 0)RETURN
       RECC=RE*RECC
    END IF
!
    IAT=0
    CALL DIAGA2(JA1,JA4,K4,IRE,IAT,RE)
    IF(IAT == 0)RETURN
    RECC=RE*RECC
!
    IAT=0
    CALL DIAGA4(JA1,JA2,K1,K2,KA,IRE,IAT,RE)
    IF(IAT == 0)RETURN
    RECC=RE*RECC
!
    IAT=0
    CALL DIAGA4(JA2,JA3,KA,K3,K4,IRE,IAT,RE)
    IF(IAT == 0)RETURN
    RECC=RE*RECC
    IF(JA1 == 1.AND.JA2 == 2)RETURN
!
    IAT=0
    CALL DIAGA1(JA1,K1,IRE,IAT,RE)
    IF(IAT == 0)RETURN
    RECC=RE*RECC
!
    ISKR=JA2-JA1
    IF(JA1 == 1)ISKR=JA2-1-JA1
    IF(ISKR <= 1)RETURN
    IAT=0
    CALL DIAGA3(JA1,JA2,K1,IRE,IAT,RE)
    RECC=RE*RECC
    RETURN
  END SUBROUTINE RECO4

  
!*******************************************************************
!                                                                  *
      SUBROUTINE RECOONESCALAR(NS,JA1,JA2,JA3,JA4,KA,IAT)
!                                                                  *
!     -------------  SECTION REC    SUBPROGRAM 05  --------------  *
!                                                                  *
!                                                                  *
!   Written by  G. Gaigalas                                        *
!   The last modification made by G. Gaigalas       October  2017  *
!                                                                  *
!*******************************************************************
        !
        IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)  :: NS, JA1, JA2, JA3, JA4, KA
      INTEGER, INTENT(OUT) :: IAT
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, IA1, IA2, IJ, IJ1, IJ2, J, NPEELGG
!-----------------------------------------------
      IAT=1
      IF(NPEEL == 1 .AND. NS == -1)RETURN
      IF(NS == -1) THEN
         NPEELGG = NPEEL
      ELSE
         NPEELGG = NS
      END IF
      IJ1=JLIST(JA1)
      IJ2=JLIST(JA2)
      IF(JA1 == 1.AND.JA2 == 2) GO TO 1
      IF(KA /= 0)GO TO 5
!
!  CASES WHEN :          KA = 0
!                  OR    JA1 = JA2
!                  OR    JA1 = 1    JA2 = 2
!
    1 DO I=1,NPEELGG
        IJ=JLIST(I)
        IF(I < NPEELGG-1) THEN
         IF(JJC1(I) /= JJC2(I))IAT=0
        END IF
        IF(KA /= 0) THEN
          IF(I == JA1) CYCLE
          IF(I == JA2) CYCLE
        END IF
        DO J=1,3
          IF(JJQ1(J,IJ) /= JJQ2(J,IJ))IAT=0
        END DO
      END DO
      RETURN
!
!  OTHER CASES
!
    5 CONTINUE
      DO I=1,NPEELGG
        IJ=JLIST(I)
        IF(I < NPEELGG-1) THEN
          IA1=JA1-1
          IA2=JA2-1
          IF(JA1 == 1)IA1=JA1
          IF(I >= IA1.AND.I < IA2)GO TO 7
          IF(JJC1(I) /= JJC2(I))IAT=0
        END IF
    7   IF(I == JA1) CYCLE
        IF(I == JA2) CYCLE
        IF((KA == 2).AND.(I == JA3)) CYCLE
        IF((KA == 3).AND.(I == JA3)) CYCLE
        IF((KA == 3).AND.(I == JA4)) CYCLE
        DO J=1,3
          IF(JJQ1(J,IJ) /= JJQ2(J,IJ))IAT=0
        END DO
      END DO
      RETURN
      END SUBROUTINE RECOONESCALAR




!***********************************************************************
!                                                                      *
      SUBROUTINE SETQNA(JA, JB) 
!                                                                      *
!   This generates the  arrays  defining  the quantum numbers of the   *
!   states involved in the  matrix  element  linking  configurations   *
!   labelled by JA, JB.                                                *
!                                                                      *
!                                           Last update: 30 Oct 1987   *
!                                                                      *
!***********************************************************************

      IMPLICIT NONE        
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: JA 
      INTEGER  :: JB 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J, K, JCNT, JCNTOP, JW1, JW2, JW 
!-----------------------------------------------
!
!
!   List parameters defining all shells in both configurations, whether
!   participating or not
!
      DO J = 1, NW 
         NQ1(J) = IQA(J,JA) 
         NQ2(J) = IQA(J,JB)
         DO K = 1, 3 
            JJQ1(K,J) = JQSA(J,K,JA) 
            JJQ2(K,J) = JQSA(J,K,JB) 
         END DO 
      END DO 
!
!   Define coupling schemes: set JLIST array to define those shells
!   which are open in either configuration, and KLIST array to locate
!   the rest. Exclude shells which are empty in both configurations
!
      NPEEL = 0 
      NCORE = 0 
      DO J = 1, NW
!         PRINT *, stype(j,ja), stype(j,jb)
         IF (stype(J,JA)==(-1) .AND. stype(J,JB)==(-1)) CYCLE  
         IF (stype(J,JA)==1 .AND. stype(J,JB)==1) THEN 
            NCORE = NCORE + 1 
            KLIST(NCORE) = J
         ELSE 
            NPEEL = NPEEL + 1 
            JLIST(NPEEL) = J
         ENDIF 
      END DO
!      PRINT *, 'setqna:', npeel, (jlist(j), j=1,nw)
!
!   Return if not more than one shell is open
!
      IF (NPEEL <= 1) RETURN  
!
!   Set arrays of coupling angular momenta interpolating closed
!   shells where necessary. Left hand side first ...
!
      JCNT = 1 
      JCNTOP = 0 
      JW1 = JLIST(1) 
      JW2 = JLIST(2) 
      IF (stype(JW1,JA) /= 0) THEN 
         JJC1(1) = JQSA(JW2,3,JA) 
         IF (stype(JW2,JA) == 0) JCNTOP = 1 
      ELSE 
         JCNTOP = 1 
         IF (stype(JW2,JA) == 0) THEN 
            JJC1(1) = JCUPA(JCNT,JA) 
            JCNT = JCNT + 1 
         ELSE 
            JJC1(1) = JQSA(JW1,3,JA) 
         ENDIF 
      ENDIF 
!
      DO J = 3, NPEEL 
         JW = JLIST(J) 
         IF (stype(JW,JA) /= 0) THEN 
            JJC1(J-1) = JJC1(J-2) 
         ELSE 
            IF (JCNTOP /= 0) THEN 
               JJC1(J-1) = JCUPA(JCNT,JA) 
               JCNT = JCNT + 1 
            ELSE 
               JJC1(J-1) = JQSA(JW,3,JA) 
            ENDIF 
            JCNTOP = JCNTOP + 1 
         ENDIF 
      END DO 
!
!   ... and repeat for right hand side
!
      JCNT = 1 
      JCNTOP = 0 
      JW1 = JLIST(1) 
      JW2 = JLIST(2) 
      IF (stype(JW1,JB) /= 0) THEN 
         JJC2(1) = JQSA(JW2,3,JB) 
         IF (stype(JW2,JB) == 0) JCNTOP = 1 
      ELSE 
         JCNTOP = 1 
         IF (stype(JW2,JB) == 0) THEN 
            JJC2(1) = JCUPA(JCNT,JB) 
            JCNT = JCNT + 1 
         ELSE 
            JJC2(1) = JQSA(JW1,3,JB) 
         ENDIF 
      ENDIF 
!
      DO J = 3, NPEEL 
         JW = JLIST(J) 
         IF (stype(JW,JB) /= 0) THEN 
            JJC2(J-1) = JJC2(J-2) 
         ELSE 
            IF (JCNTOP /= 0) THEN 
               JJC2(J-1) = JCUPA(JCNT,JB) 
               JCNT = JCNT + 1 
            ELSE 
               JJC2(J-1) = JQSA(JW,3,JB) 
            ENDIF 
            JCNTOP = JCNTOP + 1 
         ENDIF 
      END DO 
!
      RETURN  
      END SUBROUTINE SETQNA 

      


!***********************************************************************
!                                                                      *
      SUBROUTINE SNRC(IS, KAPS, KS, ND1, ND2, NE1, NE2, IBRD, IBRE) 
!                                                                      *
!   Determines the range of tensor rank NU for direct/exchange terms,  *
!   and classifies the types of radial integral.                       *
!                                                                      *
!   Input variables:                                                   *
!                                                                      *
!      IS      : Orbital labels                                        *
!      KAPS    : Values of 2*kappa                                     *
!      KS      : Values of 2*J+1                                       *
!                                                                      *
!   Outputs:                                                           *
!                                                                      *
!      ND1/NE1 : Lowest NU value for direct/exchange types             *
!      ND2/NE2 : Corresponding number of contributing NU values: NU    *
!                = ND1, ND1+2 ,..., ND1+2*(ND2-1) etc                  *
!      IBRD    : Classify types of  radial  integrals  contributing;   *
!      IBRE      negative value implies null contribution              *
!                                                                      *
!                                           Last update: 09 Oct 1992   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(OUT) :: ND1 
      INTEGER, INTENT(OUT) :: ND2 
      INTEGER, INTENT(OUT) :: NE1 
      INTEGER, INTENT(OUT) :: NE2 
      INTEGER, INTENT(OUT) :: IBRD 
      INTEGER, INTENT(OUT) :: IBRE 
      INTEGER, INTENT(IN) :: IS(4) 
      INTEGER, INTENT(IN) :: KAPS(4) 
      INTEGER, INTENT(IN) :: KS(4) 
!----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IAC, IAD, ND1A, ND2A, NE1A, NE2A 
!-----------------------------------------------
!
      ND2 = 0 
      NE2 = 0 
!
!   2.0  Form limits for direct terms
!
      IAC = 1 
      IF (KAPS(1)*KAPS(3) < 0) IAC = -1 
      IAD = 1 
      IF (KAPS(2)*KAPS(4) < 0) IAD = -1 
      ND1 = ABS(KS(1)-KS(3))/2 - 1 
      IF (IAC == (-1)) ND1 = ND1 + 1 
      IF (ND1 == (-1)) ND1 = 1 
      ND1A = ABS(KS(2)-KS(4))/2 - 1 
      IF (IAD == (-1)) ND1A = ND1A + 1 
      IF (ND1A == (-1)) ND1A = 1 
      IF (MOD(ND1 - ND1A,2) /= 0) THEN 
         IBRD = -1 
      ELSE 
         ND2 = ABS(KS(1)+KS(3))/2 
         IF (IAC == (-1)) ND2 = ND2 + 1 
         ND2A = ABS(KS(2)+KS(4))/2 
         IF (IAD == (-1)) ND2A = ND2A + 1 
         ND1 = MAX(ND1,ND1A) 
         ND2 = MIN(ND2,ND2A) 
         ND2 = (ND2 - ND1)/2 + 1 
!
!   2.1  Identify type of radial integrals
!
         IBRD = 1 
         IF (IS(1)==IS(3) .AND. IS(2)/=IS(4) .OR. IS(1)/=IS(3) .AND. IS(2)==IS(&
            4)) IBRD = 2 
         IF (IS(1)==IS(3) .AND. IS(2)==IS(4)) IBRD = 3 
      ENDIF 
!
!   3.0  Form limits for exchange terms
!
      IF (IS(1)==IS(2) .OR. IS(3)==IS(4)) THEN 
         IBRE = -1 
         RETURN  
      ENDIF 
      IAC = 1 
      IF (KAPS(1)*KAPS(4) < 0) IAC = -1 
      IAD = 1 
      IF (KAPS(2)*KAPS(3) < 0) IAD = -1 
      NE1 = IABS(KS(1)-KS(4))/2 - 1 
      IF (IAC == (-1)) NE1 = NE1 + 1 
      IF (NE1 == (-1)) NE1 = 1 
      NE1A = ABS(KS(2)-KS(3))/2 - 1 
      IF (IAD == (-1)) NE1A = NE1A + 1 
      IF (NE1A == (-1)) NE1A = 1 
      IF (MOD(NE1 - NE1A,2) /= 0) THEN 
         IBRE = -1 
         RETURN  
      ENDIF 
!
      NE2 = ABS(KS(1)+KS(4))/2 
      IF (IAC == (-1)) NE2 = NE2 + 1 
      NE2A = ABS(KS(2)+KS(3))/2 
      IF (IAD == (-1)) NE2A = NE2A + 1 
      NE1 = MAX(NE1,NE1A) 
      NE2 = MIN(NE2,NE2A) 
      NE2 = (NE2 - NE1)/2 + 1 
!
!   3.1  Identify type of radial integrals
!
      IBRE = 1 
      IF (IS(1)==IS(4) .AND. IS(2)/=IS(3) .OR. IS(1)/=IS(4) .AND. IS(2)==IS(3)&
         ) IBRE = 2 
      IF (IS(1)==IS(3) .AND. IS(2)==IS(4)) IBRE = 4 
      RETURN  
!
      END SUBROUTINE SNRC 

      
!***********************************************************************
!                                                                      *
      SUBROUTINE SPEAK(JA, JB, IA1, IB1, IA2, IB2, K, X) 
!                                                                      *
!   Output MCP  coefficients and integral parameters to COMMON block   *
!   /BUFFER/. Also print these if  IBUG1 = 1 .                         *
!                                                                      *
!                                           Last Update: 14 Oct 1992   *
!                                                                      *
!***********************************************************************
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: JA 
      INTEGER, INTENT(IN) :: JB 
      INTEGER, INTENT(IN) :: IA1 
      INTEGER, INTENT(IN) :: IB1 
      INTEGER, INTENT(IN) :: IA2 
      INTEGER, INTENT(IN) :: IB2 
      INTEGER, INTENT(IN) :: K 
      REAL(KIND=8), INTENT(IN) :: X 
!-----------------------------------------------
!
!
!      IF (IBUG1 /= 0) WRITE (99, 300) JA, JB, NP(IA1), NH(IA1), NP(IB1), NH(IB1&
!         ), NP(IA2), NH(IA2), NP(IB2), NH(IB2), K, X 
!
!   Increment counter
!
      NVCOEF = NVCOEF + 1 
!
!   Store integral indices and coefficient in COMMON/BUFFER/
!
      LABEL(1,NVCOEF) = IA1 
      LABEL(2,NVCOEF) = IB1 
      LABEL(3,NVCOEF) = IA2 
      LABEL(4,NVCOEF) = IB2 
      LABEL(5,NVCOEF) = K 
      COEFF(NVCOEF) = X 
!
      RETURN  
!
!  300 FORMAT(2(1X,1I2),4(1X,I2,A2),1X,I2,1X,1P,D19.12) 
!      RETURN  
!
      END SUBROUTINE SPEAK 


  
!***********************************************************************
!                                                                      *
      SUBROUTINE TALK(JA, JB, NU, IA, IB, IC, ID, ITYPE, COEF) 
!                                                                      *
!   Print  coefficients  and  integral  parameters  if IBUG1 > 0 and   *
!   write to disk.                                                     *
!                                                                      *
!                                           Last update: 14 Oct 1992   *
!                                                                      *
!***********************************************************************
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: JA, JB, NU, IA, IB, IC, ID, ITYPE 
      REAL(KIND=8), INTENT(IN) :: COEF 

!   Increment counter
!
      NVCOEF = NVCOEF + 1 

!   Store integral indices and coefficient in COMMON/BUFFER/
!
      LABEL(1,NVCOEF) = IA 
      LABEL(2,NVCOEF) = IB 
      LABEL(3,NVCOEF) = IC 
      LABEL(4,NVCOEF) = ID 
      LABEL(5,NVCOEF) = NU 
      LABEL(6,NVCOEF) = ITYPE 
      COEFF(NVCOEF) = COEF 
!
      RETURN  

    END SUBROUTINE TALK

END MODULE racah
