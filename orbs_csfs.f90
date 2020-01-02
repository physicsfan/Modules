MODULE orbs_csfs
!***********************************************************************
!      Use constants
      Use case

!     Data structures for orbitals
      !> orblist -- list of orbitals as one string
      CHARACTER(LEN=4*NNNW) :: orblist
      !> nh      -- the symmetry characters for relativistic orbitals
      CHARACTER(LEN=2),DIMENSION(:), allocatable :: nh
      !> np      -- n principal quantum number
      !> nkl     -- l quantum number
      !> nkj     -- 2*j quantum number
      !> nak     -- kappa quantum number
      !> Oparity -- orbital parity
      INTEGER(kind=1), DIMENSION(:), allocatable :: NP, NKL, NKJ, NAK, Oparity 
      !> el      -- character symbol for the orbital
      CHARACTER(LEN=4),DImension(:), allocatable  :: el  

      !> Data structures for CSFs
      !> iqa     -- occupation numbers 
      !> stype   -- shell type -- unocupied, open, full
      !> sparity -- shell parity
      !> jcupa   -- coupling of open shells (others are zero)
      INTEGER(kind=1), DIMENSION(:,:), allocatable :: iqa, stype, sparity, jcupa
      !> jqsa    -- shell quantum numbers of three types: 
      !              1) seniority or quasi-spin,
      !              2) q_M (needed for some g^N CSFs)
      !              3) J quantum number (in 2J+1) representation 
      INTEGER(kind=1), DIMENSION(:,:,:), allocatable :: jqsa
!
!-----------------------------------------------

CONTAINS
 
   SUBROUTINE load_orbs(nfile)           
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nfile
      INTEGER :: i, ii, iib2, npeel
      CHARACTER(LEN=150) :: str1    !> temporary space for core orbs
      CHARACTER(LEN=500) :: str2    !> temporary space for peel shells

!
!   Read the records
      READ (nfile,'(A)') str1 
      Print *, str1
      If (str1(1:5) .ne. 'Core') then
         Write(6, *) 'Error in rcsf.inp file'
         stop
      End IF

!    .. Determine the core orbitals     
      Read (nfile,'(A)') str1
      ncore = ( len_trim(str1)+ 1)/5
!
!    .. Determine the core orbitals (skip line of input)
      Read (nfile,'(/(A))')  str2
      npeel = ( len_trim(str2)+ 1)/5
      nw = ncore + npeel
      if ( nw > NNNW) then
         write(6, '(A, i3, A)')                                             &
                     'Error:', nw, ' orbitals exceeds the maximum of 127'
         stop
      Else
!        We know the number of orbitals, allocate memory
         allocate(nh(nw), np(nw), nkl(nw), nkj(nw), nak(nw), oparity(nw))
         np=0; nkl=0; nkj=0; nak=0; oparity=0
         allocate(el(nw))
      End if

      if (ncore > 0) read(str1,'(100(1X,A4))') (el(i),i=1,ncore)
      read(str2,'(100(1X,A4))') (el(ncore+i),i=1,npeel)
      
!     ...  generate orblist
      write(orblist, '(100(A4))') (el(i),i=1,nw)

!     ...  Load the CSFs (we are assuming only one block)
        
      READ (orblist, '(100(I2,A2))') ( NP(i), NH(i), i=1,nw)
      !Print *, 'Orblist =', trim(orblist)
      !Print *, 'NP(i)=', np(1:nw)
      !Print *, 'NH(i)=', nh(1:nw)
      Do  i=1,nw
          ii = (index(symlist,nh(i)) +1)/2
          IIB2 = II/2
          NKL(i) = IIB2
          IF (MOD(II,2) == 1) THEN
            NAK(i) = (-IIB2) - 1
            NKJ(i) = II
          ELSE
            NAK(i) = IIB2
            NKJ(i) = II - 1
          ENDIF
          oparity(i) = 1
          If ( mod(nkl(i),2) .gt. 0 ) oparity(i) = -1
      End do
      !print *, 'nkl = ', nkl(1:nw)
      !print *, 'nak = ', nak(1:nw)
      !print *, 'nkj = ', nkj(1:nw)
      !print *, 'oparity =', oparity(1:nw)
    
    RETURN
  END SUBROUTINE load_orbs

  SUBROUTINE Load_csfs(nfile)
!                                                                      *
    IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
    INTEGER, INTENT(in) :: nfile
    INTEGER :: I, j, II, nsub, ios, loc, ishell, iopen
    INTEGER, DIMENSION(MX) :: jcup, jqs1, jqs3
    CHARACTER(LEN=9)      :: strsub 
    CHARACTER(LEN=9), DIMENSION(MX) :: jstr
    CHARACTER(LEN=9*MX)   :: str
    CHARACTER(LEN=4),DIMENSION(:), ALLOCATABLE  :: elsym
    INTEGER(kind=1), DIMENSION(:), ALLOCATABLE  ::  iocc
    CHARACTER(Len=12)        :: convert
!
!-----------------------------------------------

!   Load orbitals first....
    CALL load_orbs(nfile)
      
!     .. initialize CSF data
    ALLOCATE(iqa(nw,ncsf), stype(nw,ncsf), sparity(nw,ncsf), jcupa(nw, ncsf))
    ALLOCATE(jqsa(nw, 3, ncsf))
    IQA = 0; JQSA = 0;  JCUPA =0; stype = -1; sparity = 1
!
!   Set the occupation and subshell quantum number array elements
!   in IQ, JQS for the core subshells for all CSFs
!       print *, 'ncore', ncore
! 
    DO I = 1, NCORE
       IQA   (i, 1:NCSF) =  NKJ(I)+1
       JQSA (i,3,1:NCSF) = 1
       stype (i, 1:NCSF) = 1
       sparity(i,1:NCSF) = 1
    END DO

    READ (nfile, '(8X,I3)') nblock
    READ (nfile, *)

    ALLOCATE(jend(nblock), jblk(nblock), pblk(nblock))

    ! allocate local arrays and initialize
    ALLOCATE(elsym(nw), iocc(nw))
    jqs1=0; jqs3=1; jcup=0;  iocc=0

    ncsf =0
    DO j = 1,nblock

!     This READ makes the routine load the entire file (all blocks)
       DO
          READ (NFILE, '(A)', IOSTAT=IOS) STR

!     .. deal with case where a line2 is empty
          IF (len_TRIM(str) .EQ. 0) THEN
             jqs3(1:nsub) = 1
             CYCLE
          END IF

!     .. the end of the block has been reached
          IF (IOS .NE. 0 .OR. str(1:2) .EQ. ' *' )  THEN
             ! Print '(A, I3 //)', 'End of CSF', ncsf
             jend(j) = NCSF
             jblk(j) = jcup(nsub)
             pblk(j) = PRODUCT(sparity(1:nw, ncsf), MASK=stype(1:nw,ncsf)==0)
             ! print *, 'Block parity =', pblk(j)
             PRINT *, 'Properties for Block ', len_TRIM(convert(j))
             PRINT *, '  Ending CSF         ', jend(j)
             PRINT *, '  J and parity       ', jstr(nsub)(7:9)
             jqs1=0; jqs3=1; jcup=0;  iocc=0
             IF (str(1:2) .EQ. ' *')  CYCLE
             IF (IOS .NE. 0) EXIT
          END IF
!
!     Determine what kind of line we have;
          IF (IOS==0 .AND. STR(1:2)/=' *') THEN

!        The string specifies the subshells and their occupation: line1 
             IF(INDEX(STR,'(')/=0) THEN
                PRINT *, 'Line1:',  str
                nsub = (len_TRIM(str))/9 
                !          .. shells line
                READ (str, '(127(1X,A4,1X,I2,1X))') (elsym(i), iocc(i),i=1,nsub)
                CYCLE
             ENDIF

!        Check for the left to right coupling, total J and parity string: line3
             IF (INDEX(STR,'+')/=0 .OR. INDEX(STR,'-')/=0) THEN
                ! Read( str, '(3x,127(A9))' ) (Jstr(i), i=1,nsub)
                PRINT *,  'Line3: ', str
                DO i=1, nsub
                   strsub = jstr(i)
                   IF (len_TRIM(strsub) > 0) THEN
                      loc = INDEX(strsub,'/') 
                      IF ( loc  > 0 ) THEN
                         READ(strsub(5:7), '(I3)') jcup(i)
                         jcup(i) = jcup(i) + 1
                      ELSE 
                         IF (i == nsub) THEN
                            READ(strsub(5:7), '(I3)') jcup(i)
                         ELSE 
                            READ(strsub(7:9), '(I3)') jcup(i)
                         END IF
                         jcup(i) = 2*jcup(i) +1
                      END IF
                   ELSE
                      jcup(i) = 1              
                   END IF
                   ! Print *, 'Line3:i', i, iocc(i),jqs1(i),jqs3(i),jcup(i)
                END DO
!
!          we have a CSF -- store the values
                ncsf = ncsf+1
                iopen = 0
                DO i = 1,nsub
                   ishell = INDEX(orblist, elsym(i))/4 +1
                   iqa(ishell,ncsf) = iocc(i)               
                   jqsa(ishell,1,ncsf) = jqs1(i)
                   jqsa(ishell,3,ncsf) = jqs3(i)
                   sparity(ishell, ncsf) = oparity(ishell)**iocc(i)
                   stype(ishell, ncsf) = 1  
                   IF (iocc(i) .LT. nkj(ishell) +1) THEN
                      stype(ishell, ncsf) = 0
                      IF (iopen > 0) jcupa(iopen, ncsf) = jcup(i)
                      iopen = iopen +1
                   END IF
                   ! Print *, 'i, iocc, jqs(2), jcup',i, iocc(i), jqs1(i), jqs3(i), jcup(i)
                END DO

                CYCLE
             END IF

!   The string must contain the subshell quantum numbers: line2
             PRINT *, 'Line2: str=',str
             READ (str, '(127(A9))') (jstr(i), i=1, nsub)
             DO i =1, nsub
                IF (len_TRIM(jstr(i)) > 0) THEN
                   strsub = jstr(i)
                   loc = INDEX(strsub,'/') 
                   IF ( loc  > 0 ) THEN
                      ! print *, 'strsub(5:7)', strsub(5:7)
                      READ(strsub(5:7), '(I3)') jqs3(i)
                      ! print *, 'jqs3(i) =', jqs3(i), 'i=',i
                      jqs3(i) = jqs3(i)+1
                      ! print *, 'jqs3(i) =', jqs3(i), 'i=',i
                   ELSE 
                      READ(strsub(7:9), '(I3)') jqs3(i)
                      jqs3(i) = 2*jqs3(i)+1
                      loc = INDEX(strsub, ';')
                      IF (loc >0) READ(strsub(loc-3:loc-1), '(I3)') jqs1(i)
                   END IF
                ELSE
                   jqs3(i) = 1
                END IF
             END DO
             CYCLE
          END IF
       END DO
    END DO
    !  Clean up
    DEALLOCATE(elsym, iocc)
    RETURN
  END SUBROUTINE Load_CSFs
END MODULE orbs_csfs
