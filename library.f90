
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
