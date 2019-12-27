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
