module constants 

  ! Constants or parameters that define the implementation

  IMPLICIT NONE

  
   ! Basic Data 
   !> NNNW  -maximum number of orbitals 
   INTEGER, PARAMETER :: NNNW=127 
   !>  key  -- key for encoding data
   INTEGER, PARAMETER :: key = 128 
   !> symlist :: list of symmetries supported in the current version 
   CHARACTER(LEN=38)  :: SYMlist='s p-p d-d f-f g-g h-h i-i k-k l-l m-m '  
   !> kaplist :: list establishing the kappa order 
   INTEGER, DIMENSION(19) :: kaplist 
   DATA kaplist / -1, 1, -2, 2, -3, 3, -4, 4, -5, 5, -6, 6, -7, 7, -8, & 
        8, -9, 9, -10 / 
   !Note: index(kaplist, kappa) defines the kappa order in the list. 

   ! Parameters for matrix elements 
   REAL(kind=8) :: cutoff = 1.0d-10 

  ! Define UNIT numbers for I/O in grasp: 
   INTEGER ::  istderr=0, istdin=5,  istdout=6
   
  ! Define UNIT numbers for RANGULAR I/O: 
  INTEGER, PARAMETER :: nfile = 10 
  INTEGER, PARAMETER :: outfile=13 
  INTEGER, PARAMETER :: isodatafile = 22
   
  ! ... some useful physical constants
  ! ... (based on NIST compilation 2001) from DBSR_HF   (variables in GRASP)
  
  REAL(kind=8), PARAMETER ::                           &
       bohr_radius_in_cm      =    0.5291772083d-8,  & ! AINFCM
       ao_cm                  =    0.5291772083d-8,  &
       one_over_alpha         =    137.03599976d+0,  &  !ALFAI
       alpha                  =    7.297352533d-3,   &
       c_au                   =    137.03599976d+0,  & !CCPMS (cm/s)
       electron_charge_in_C   =    1.602176462d-19,  &
       electron_charge_in_esu =    4.803204673d-10   &  !EESU
       electron_mass_in_g     =    9.10938188d-28,   &  !EMEG
       electron_mass_in_amu   =    5.485799110d-4,   &  !EMEAMU
       proton_mass_in_amu     =    1.00727646688d+0, &  !EMPAMU
       proton_electron        =    1836.1526675d+0,  &
       electron_proton        =    5.446170232d-4,   &
       hbar_in_J_s            =    1.054571596d-27,  &  !HBARES
       hbar_in_eV_s           =    6.58211889d-16,   &
       au_eV_inf              =    27.2113834d+0,    &
       au_cm_inf              =    219471.62d+0,     &
       RinfeV                 =    13.605693009d+0,  &  !Ryberg in eV
       Rinfk                  =    109737.31568508d+0 & !Rydberg in Kaysers
       
end module constants 

MODULE case
  IMPLICIT NONE

  ! About the Atom    

  !> Z       :: atomic number
  !> A       :: nuclear mass
  !> APRAM, CPRAM ::Paremeters defining the nuclear skin
  !> RRMS    :: root means square radius
  INTEGER    :: z, a
  REAL(kind=8)   :: apram, cpram, rrms
  
  ! Type  of calculation
  ! c_speed    :: Speed of light used in the calculation
  !               (infinity for non-relativistic calculation)
  REAL (kind=8)  :: c_speed =  c_au

  ! Basic Size dimensions of a calculation
  !> ncore   -- number of core orbitals
  !> nw      -- number of wavefunctions or orbitals
  !> nblock  -- number of Jblocks of the same total J and parity
  !> lmax    -- maximum value of l quantum number 
  !> Int_num -- number of integrals
  !> kmax    -- maximum value of k = 2*lmax
  !> mx      -- maximun number of subshells in a csf (<= NNNW)
  !> ncsf    -- number of CSFs
  INTEGER :: ncore, nw, nblock, lmax, int_num, kmax, ncsf, mx

  ! JBlock data structures
  !> jend    -- array indicating the last CSF in the i'th JBLOCK
  !> jblk    -- array for the J-value of each block
  !> pblk    -- parity of each block
  INTEGER, DIMENSION(:), ALLOCATABLE :: jend, jblk, pblk
      

END MODULE case
MODULE load_mod
!***********************************************************************
      Use constants
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
END MODULE Load_mod
   
!***********************************************************************
!>   Funftion that converts a binary integer to a character string
      Function Convert (intnum) 
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
      End Function Convert


module GenInt_module

  Use load_mod

  implicit none

! Data structures for integral lists
  !> Int_label --  I(a,b) label of ia, ib encoded for category 1
  !>               Rk(abcd) label ia, ib, ic, id for category k+1
  !>Int_end    --  integer index array indicating the last position
  !                 for integrals in category k+2, .. j2max+2
  Integer, dimension(:), allocatable      ::  Int_label, Int_end

  ! Needed functions
  !> triangrk  -- computes the triangle relationship
  Logical, external :: triangrk

contains

  subroutine genint

!>  Generate a list of integral with the following types:
!> 1. I (or T) integrals  I(ab)
!> 2. Rk integrals for k= 0,2*kmax of Symmetry Rk(ab,cd)

!> Strategy: Each type of integral has a last value. The 
!> first category starts at i=1 whereas all others start
!> at the last value of the previous category +1

    !   L o c a l  V a r i a b l e s
    LOGICAL :: gen
    integer :: i, ia, ib, ic, id, k, n, m
    
    kmax = 2*maxval(nkl)
    
    allocate(Int_end(kmax+2))
    
! Generate the list of I Integrals that arise from a set of orbitals.
! The subroutine has two phases, either return the list with unevaluated
! integrals or a list with evaluated integrals.
!
! This is controlled by the argument: "evaluated"

      gen = .false.
      n=0
      do i = 1, 2
         do ia = 1, nw
            do ib = ia, nw
               if(nak(ia) == nak(ib)) then
                  n = n + 1
                  if (gen) then
                     Int_Label(n)  = ia*key + ib
                  endif
               end if
            end do
         end do
         Int_end(1)  = n
      !
      ! Generate the list of Rk Integrals that arise from a set of orbitals.
      ! The subroutine has two phases, either return the list with unevaluated
      ! integrals or a list with evaluated integrals.
      !
         do k = 0, kmax
            do ia = 1, nw
               do ic = ia, nw
                  if (triangrk(nkl(ia), k, nkl(ic))) then
                     ! if (evaluate)  !Compute Yk(ac)
                     do ib = ia, nw           
                        do id = ib, nw
                           if (triangrk(nkl(ib), k, nkl(id))) then
                              n = n + 1
                              if (gen) then
                                 int_Label(n) = ((ia*key+ic)*key+ib)*key+id
                              end if
                           end if
                        end do
                     end do
                  end if
               end do
            end do
            if (gen) Int_end(k+2) = n 
         end do

         
         
         ! Allocate memory for integral book keeping
         if (.not. gen) then
            ! Number of integrals
            Int_num = n
            allocate(int_Label(1:n))
            gen = .true.
            n=0
         end if
      end do
      
    end subroutine genint
          
    subroutine write_labels
     
    Integer :: i, lab, k, ia, ib, ic, id, istart, n

    Print *, int_end(1:kmax+2)
    n=1
    do i =1, Int_end(1)
        ib =  MOD(int_Label(i), key)
        ia = int_label(i)/key
        Write (6, FMT="(I2, 2X, 'I(', A4, ',', A4,')')" ) n,  el(ia), el(ib)
        n = n+1
    end do  
    istart = int_end(1)+1
    do k = 0, kmax
       do i = istart, Int_end(k+2)
           lab = int_label(i)
           id = mod(lab, key)
           lab = lab/key
           ib = mod(lab,key)
           lab = lab/key
           ic = mod(lab,key)
           ia = lab/key
           write (6, '(I2,2X,A,I2,A,4(A4,1X),A)' ) n, 'R', k,'(', el(ia), el(ib), el(ic), el(id), ')'
           n = n+1
        ENd do
        istart = Int_end(k+2)+1
    End Do
    end subroutine write_labels
           
EnD Module Genint_module
       
       
  LOGICAL FUNCTION TRIANGRK (LA, K, LB) 
!                                                                      *
!   Written by Per Jonsson                                             *
!                                                                      *
!***********************************************************************
!   M o d u l e s 
!-----------------------------------------------
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

MODULE matrix_elements
  USE genint_module

! Parameters for matrix elements
  INTEGER :: keysq = 128**2
  
! Needed subroutines
  !> onescalar  -- computes the Tcoeffs and returns corresponding ia, ib
  !> rkco    -- computes the Vcoeffs and returns ia, ib, ic, id

! Development
  INTEGER :: tdatfile, vdatfile
  
  
CONTAINS
  
  SUBROUTINE genmat(outfile)

!   Generate a list of coeffisients with the following types:
!   1. Tcoeffs(ab) associated with I integrals
!   2. Vkcoeffs(ab,cd) associated with Rk integrals

    !  L o c a l  V a r i a b l e s
    INTEGER      :: blk, coeff_label, i, ia, ib, ic, id, iswap, outfile
    INTEGER      :: j, ja, jb, k, n, nz, nswap, nt, nv, t_total, v_total
    REAL(kind=8), DIMENSION(:), ALLOCATABLE :: tshell, vshell
    INTEGER, DIMENSION(:), ALLOCATABLE :: tindex, vindex
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: label

    
    ! initialize
    ALLOCATE(tshell(nnnw), tindex(nnnw))
    ALLOCATE(vshell(nnnw), vindex(nnnw))
    ALLOCATE(label(5,nnnw))
    t_total = 0; v_total = 0


    WRITE(6, '(a,i3)' ) 'NW', nw 
    WRITE(6, '(a,i3)') "Nblocks", nblock
    WRITE(6, '(a,i3)') "Ncsfs", ncsf
 
    WRITE(6, '(//,a)') "Integral index output"
    WRITE(6, '(/,a)') " index   k   ia   ib   ic   id"

    
    DO blk = 1, nblock

!   JA and JB respectively refer to the initial and final states
!   in the list of NCF configurations
       !
outer:   DO jb = 1, ncsf
inner:     DO ja = 1, jb           ! ja and jb loop over all upper  matrix elements
             
            label = 0
!
!  Generate t coefficients.
!
             nt = 0
             tshell = 0.0; tindex = 0
             ! onescalar is a fake function that reads in data from old rangular
             CALL onescalar(ja,jb,ia,ib,tshell)      
             DO i = 1, SIZE(tshell)
                IF ((ia .NE. 0) .AND. (ia .NE. ib) .AND. (ABS(tshell(i)) .GT. cutoff)) THEN             
                   nt = nt + 1
                   IF (ia .GT. ib) THEN   !swap if necessary
                      iswap = ib
                      ib = ia
                      ia = iswap
                   ENDIF
                   coeff_label= ia*key + ib
                   tindex(i) = find_index(-1, coeff_label)     ! found in genint
                ENDIF
             END DO
!
!  Generate V coefficients; ac and bd are the density pairs.
!      
             nv = 0
             vshell = 0.0; vindex = 0
             ! rkco ia a fake function that reads in data from old rangular
             CALL rkco(ja,jb,label,vshell)
             IF(label(1,1) == -1) EXIT outer  ! end of a block has been reached
             DO j = 1, SIZE(vshell)
                IF(ABS(vshell(j)) .GT. cutoff) THEN
                   nv = nv + 1
                   k  = label(1,j)
                   ia = label(2,j)
                   ib = label(3,j)
                   ic = label(4,j)
                   id = label(5,j)
                   IF (.NOT.((k.EQ.0).AND.(ia.EQ.ic).AND.(ib.EQ.id))) THEN
                      ! Swap index to make sure IA <= IC, IB <= ID, IA <= IB
                      IF (ia .GT. ic) THEN
                         iswap = ic
                         ic = ia
                         ia = iswap
                      ENDIF
                      IF (ib .GT. id) THEN
                         iswap = id
                         id = ib
                         ib = iswap
                      ENDIF
                      IF (ia .GT. ib) THEN
                         ! need to swap both pairs of indicies
                         iswap = ib
                         ib = ia
                         ia = iswap
                         iswap = id
                         id = ic
                         ic = iswap
                      END IF
                   ENDIF
                   coeff_label = ((ia*key+ic)*key+ib)*key+id                        
                   vindex(j) = find_index(k,coeff_label)
                ENDIF
             END DO

             ! output matrix element data to file and screen
             ! row, col, num_int(t+V), coeff(), int_index()
             WRITE(outfile,'(i3,i3,i3,i3)') ja, jb, nt, nv
             IF(nt > 0) THEN
                DO i = 1, nt
                   WRITE(outfile, '(f9.4,i8)') tshell(i), tindex(i)
                END DO
             ENDIF
             IF(nv > 0) THEN
                DO i = 1, nv
                   WRITE(outfile, '(f9.4,i8)') vshell(i), vindex(i)
                   WRITE(6, '(i5,5i5)') vindex(i), (label(j,i), j=1,5)  
                END DO
             END IF
             WRITE(outfile,*)
             WRITE(6,*)
             t_total = t_total + nt
          END DO inner
          v_total = v_total + nv
       END DO outer
       WRITE(outfile,'(a,/)') "*"
       WRITE(6, '(a)') "*"
    END DO

    WRITE(6,'(//,a,i3)') 'Total number of t coefficients generated:', t_total
    WRITE(6,'(a,i3)') 'Total number of v coefficients generated:', v_total
    
    deallocate(tshell, tindex, vshell, vindex, label)
 
  END SUBROUTINE genmat

    !> function to find the location of an integral in the list
    INTEGER FUNCTION find_index ( k, label )

    !> k     -- integral type, k<0 is I integral
    !> index -  index of an integral in canonical order
    !           (ia <= ib; if k <0)
    !           (ia <= ic; ib <= id; ia <= ib  if k>=0)
    Integer, intent(in)  :: k, label

    !>  il -- lower bound
    !>  iu -- upper bound
    !>  im -- mid-point
    Integer  :: il, im, iu
    
    if (k < 0) then
       il =1; iu = Int_end(1)
    else
       il = Int_end(k+1)+1 ; iu = Int_end(k+2)
    end if

    do
      if (iu-il >  1) then
         im = (iu + il)/2
         if ( label > int_label(im)) then
           il = im
         else if ( label < int_label(im)) then
           iu = im
         else if (label == int_label(im)) then
           index = im    
           exit
         end if
       else if (iu-il == 1) then
         if ( label == int_label(iu)) then
           index = iu
           exit
         else if ( label ==  int_label(il)) then
           index = il 
           exit
         else
            index = -1
            Print *, index, ' not found'
            exit
         end if
      end if
    end do
  END FUNCTION find_index
  
End MODULE matrix_elements

!====================================================================
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


MODULE radial
!***********************************************************************
      Use Grid
   
      IMPLICIT NONE
   
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
      INTEGER, Dimension(:), allocatable :: npt, nodes

 
   CONTAINS
      !> allocate memory for radial functions
      SUBROUTINE allocate_radials
         allocate(P(npX,nw),Q(npX,nw),DP(npX,nw), DQ(npX,ns)
         allocate(e(nw),scf(nw),npt(nw), nodes(nw))
      END SUBROUTINE allocate_radials


END MODULE radial
    
