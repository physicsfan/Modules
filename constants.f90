  module constants 

  ! Constants or parameters that define the implementation

  IMPLICIT NONE
  
   ! Basic Data 
   !> NNNW  -maximum number of orbitals 
   INTEGER, PARAMETER :: NWX=127 
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
   REAL(kind=8) :: eps    = 1.0d-8

   
  ! Define UNIT numbers for I/O: 
  INTEGER, PARAMETER :: nfile   = 10 
  INTEGER, PARAMETER :: outfile = 13
  INTEGER, PARAMETER :: infile  = 14 
  INTEGER, PARAMETER :: isofile = 22
  INTEGER, PARAMETER :: rwfndat = 23
  INTEGER, PARAMETER :: rwfnout = 24
  INTEGER, PARAMETER :: rwfnin  = 25
  INTEGER, PARAMETER :: logfile = 26
  INTEGER, PARAMETER :: sumfile = 27
  
  ! Define UNIT numbers for I/O in grasp:
  INTEGER, PARAMETER ::  istde=0, istdi=5,  istdo=6

  ! Define kind-constants
!!$  INTEGER, PARAMETER :: int4 = selected_int_KIND(9)
!!$  INTEGER, PARAMETER :: int2 = selected_int_KIND(4)
!!$  INTEGER, PARAMETER :: int1 = selected_int_KIND(2)
!!$  
!!$  INTEGER, PARAMETER :: sp = KIND(1.0)
!!$  INTEGER, PARAMETER :: dp = selected_real_KIND(2*PRECISION(1.0_sp))
!!$  INTEGER, PARAMETER :: qp = selected_real_KIND(2*PRECISION(1.0_dp))
  
  ! Define simple constants
  REAL(kind=8), PARAMETER :: zero  = 0.0,    &
       one   = 1.0d0,    &
       two   = 2.0d0,    &
       three = 3.0d0,    &
       four  = 4.0d0,    &
       five  = 5.0d0,    &
       six   = 6.0d0,    &
       seven = 7.0d0,    &
       eight = 8.0d0,    &
       nine  = 9.0d0,    &
       ten   = 10.0d0,   &
       eleven= 11.0d0,   &
       tenth = one/ten,  &
       half  = one/two,   &
       third = one/three

  real(8) :: PI = 4.0D0*ATAN(1.0D0)

! ... some useful physical constants
! ... (based on NIST compilation 2001) from DBSR_HF   (variables in GRASP)

      real(kind=8), parameter ::                        &
      bohr_radius_in_cm      =    0.5291772083d-8,  &  !AINFCM
      ao_cm                  =    0.5291772083d-8,  &
      one_over_alpha         =    137.035999084d+0,  &  !ALFAI
      alpha                  =    7.2973525693d-3,   &
      c_au                   =    137.035999084d0,   &  !CCPMS (cm/s)
      electron_charge_in_C   =    1.602176462d-19,  &
      electron_charge_in_esu =    4.803204673d-10,   &  !EESU
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
      Rinfk                  =    109737.31568508d+0  !Rydberg in Kaysers

end module constants 
