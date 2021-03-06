
!***********************************************************************
!>   Funftion that converts a binary integer to a character string
      FUNCTION Convert (intnum) 
!                                                                      *
!>   Written by  C. Froese Fischer                          June, 2019  
!                                                                      *
!***********************************************************************
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!>    binary integer to be converted
      INTEGER,   INTENT(IN)    :: intnum
!>    character string for the decimal number
      Character(Len=12)        :: convert
!-----------------------------------------------

      write(convert,*) intnum
      convert = adjustl(convert)

      Return
    END FUNCTION Convert


    
SUBROUTINE  getrlist(nlist, list, str, nrlist, rlist)


!>   Determine the list of orbitals the user selected.
!>   List   :: orbitls to be determined
!>   nlist  :: number of orbitals
!>   str    :: user's reponse to orbitals for an option
!!     The basic input format is a list of orbitals, separated by 
!!     space or comma, with a  possible wild card.  * is equivalent to 
!!     all orbitals, n* to all orbitals with a given n, or orbitals 
!!     in non-relativistic notation
!>   nrlist :: number of orbitals selected
!>   rlist  :: list of orbitals selected

  USE  orbs_csfs
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nlist
  INTEGER, DIMENSION(nw), INTENT(in) :: list
  CHARACTER(len=4*NWX), INTENT(inout)      :: str
  INTEGER, INTENT(out) ::  nrlist
  INTEGER, DIMENSION(nw), INTENT(out) :: rlist
  
  !      local variables
  INTEGER :: i, ipos, iorb,  len, lens, n
  
  
  nrlist = 0
  rlist = 0
  DO
     str = ADJUSTL(str)   
     len = len_TRIM(str)
     PRINT *, 'len =',  len
     IF (len .EQ. 0) EXIT
     ! scan for blank or comma
     ipos = SCAN(str(1:len), ' ,')   
     PRINT *, 'S1: ipos =', str(1:len), ipos
     IF (ipos .GT. 0) THEN
        lens = ipos-1
     ELSE
        str(len+1:len+1) = ' '
        lens = len + 1
     END IF
     ! when all orbitals are slected
     IF (str(1:lens) .EQ. '*') THEN
        nrlist = nlist
        rlist = list
        ! all orbitals selected
        EXIT   
     ELSE
        ipos = INDEX(str(1:lens),'*')       
        !  we have n* format
        IF (ipos .GT. 0) THEN     
           PRINT *, ipos
           READ(str(1:ipos-1), *) n 
           ! find all orbits with np(i) = n 
           DO i = 1,nw
              IF (np(i) .EQ. n) THEN
                 nrlist = nrlist +1
                 rlist(nrlist) = i
              END IF
           END DO
        ELSE
           ! we have a separate orbital  (say 2p which includes ( 2p-, 2p )
           iorb = INDEX(orblist, str(1:lens)) 
           IF (iorb .GT. 0) THEN
              nrlist = nrlist+1
              rlist(nrlist) = iorb/4 +1
              IF (lens .LT. 3) THEN
                 iorb = INDEX(orblist(ipos+1:nw), str(1:lens))
                 IF (iorb .GT. 0) THEN
                    nrlist = nrlist+1
                    rlist(nrlist) = iorb/4 + 1
                 END IF
              END IF
           ELSE 
              PRINT *, 'Orbital ', str(1:lens), ' not found in list'
              STOP
           END IF
        END IF
     END IF
     ! blank out the first lens characters of str
     
     DO i = 1, MIN(lens+1,len) 
        str(i:i) = ' '
     END DO
  END DO
  
END SUBROUTINE getrlist


!> function to find the location of an integral in the list
INTEGER FUNCTION IINDEX ( k, label )
  USE integrals,  ONLY: int_label, int_end  
  !> k     -- integral type, k<0 is I integral
  !> iindex -  iindex of an integral in canonical order
  !           (ia <= ib; if k <0)
  !           (ia <= ic; ib <= id; ia <= ib  if k>=0)
  INTEGER, INTENT(in)  :: k, label
  
  !>  il -- lower bound
  !>  iu -- upper bound
  !>  im -- mid-point
  INTEGER  :: il, im, iu
  
  IF (k < 0) THEN
     il =1; iu = Int_end(1)
  ELSE
     il = Int_end(k+1)+1 ; iu = Int_end(k+2)
  END IF
  
  
  ! perform binary search
  DO
     IF (iu-il >  1) THEN
        im = (iu + il)/2
        IF ( label > int_label(im)) THEN
           il = im
        ELSE IF ( label < int_label(im)) THEN
           iu = im
        ELSE IF (label == int_label(im)) THEN
           iindex = im    
           EXIT
        END IF
     ELSE IF (iu-il == 1) THEN
        IF ( label == int_label(iu)) THEN
           iindex = iu
           EXIT
        ELSE IF ( label ==  int_label(il)) THEN
           iindex = il 
           EXIT
        ELSE
           iindex = -1
           PRINT *, iindex, ' not found'
           EXIT
        END IF
     END IF
  END DO
  
END FUNCTION iindex


  
!*******************************************************************
!                                                                  *
      INTEGER FUNCTION ITTK (I, J, K) 
!                                                                  *
!     CHECKED TRIANGULAR CONDITIONS FOR   I/2, J/2, K/2.           *
!     I+J>=K, I+K>=J, J+K>=I,                                      *
!     I/2+J/2+K/2 - WHOLE NUMBER                                   *
!     ITTK=1 -   IF NOT SATISFY                                    *
!     ITTK=0 -   IN OVER CASES                                     *
!                                                                  *
!   Written by G. Gaigalas,                                        *
!   Vilnius,  Lithuania                             December 1993  *
!   The last modification made by G. Gaigalas       October  2017  *
!                                                                  *
!*******************************************************************

        IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: I, J, K
!-----------------------------------------------
      ITTK = 0 
      IF (IABS(I - J) > K) RETURN  
      IF (I + J < K) RETURN  
      IF (MOD(I + J + K,2) /= 0) RETURN  
      ITTK = 1 
      RETURN  
      END FUNCTION ITTK 

      
!***********************************************************************      
!                                                                      *
    LOGICAL FUNCTION TRIANGRK (LA, K, LB) 
!                                                                      *
!   Written by Per Jonsson                                             *
!                                                                      *
!***********************************************************************
    IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
    INTEGER(kind=1), INTENT(IN) :: LA, LB
    INTEGER, INTENT(IN) :: K 
!
!   Perform the triangularity check
!
    IF (MOD(K + LA + LB,2) /= 0) THEN 
       TRIANGRK = .FALSE. 
    ELSE 
       IF (ABS(LA - LB) > K) THEN 
          TRIANGRK = .FALSE. 
       ELSE IF (LA + LB < K) THEN 
          TRIANGRK = .FALSE. 
       ELSE 
          TRIANGRK = .TRUE. 
       ENDIF
    ENDIF
    RETURN  
  END FUNCTION TRIANGRK



  REAL(kind=8) FUNCTION wt(j, blck)
    
    IMPLICIT NONE

    INTEGER, INTENT(in) :: j, blck


  END FUNCTION wt
