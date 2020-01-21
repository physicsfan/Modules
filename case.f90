MODULE case
  Use constants
  IMPLICIT NONE

  ! Basic Size dimensions of a calculation
  !> ncore   -- number of core orbitals
  !> nw      -- number of wavefunctions or orbitals
  !> nblock  -- number of Jblocks of the same total J and parity
  !> lmax    -- maximum value of l quantum number 
  !> Int_num -- number of integrals
  !> kmax    -- maximum value of k = 2*lmax
  !> mx      -- maximun number of subshells in a csf (<= NNNW)
  !> ncsf    -- number of CSFs
  !> m_el    -- number of matrix elements
  INTEGER :: ncore, nw, nblock, lmax, int_num, kmax, ncsf, mx, m_el

  ! JBlock data structures
  !> jend    -- array indicating the last CSF in the i'th JBLOCK
  !> jblk    -- array for the J-value of each block
  !> pblk    -- parity of each block
  !> m_end   -- last matrix element for each block
  INTEGER, DIMENSION(:), ALLOCATABLE :: jend, jblk, pblk, m_end

! About the Atom


  !> Z       :: atomic number
  !> A_mass  :: nuclear mass number
   Character(20) :: nuclear = 'Fermi'
  !> nuclear = point   - point nuclear
  !> nuclear = uniform - uniform distribution
  !> nuclear = fermi   - Fermi distibution
  !> r_uniform :: parameter for uniform distribution
  !> ro_uniform :: parameter for uniform distriubtion
  !> A_fermi, C_fermi ::Paremeters defining the nuclear skin
  !> RRMS    :: root means square radius
  !> I_nuc   :: nuclear spin (I) (in units of h/2*pi)
  !> D_nuc   :: nuclear dipole moment (in nuclear magnetoms)
  !> Q_nuc   :: nuclear quadrupole moment (in barns)
  Integer  :: z, a_mass
  Real(8)  :: r_uniform,  a_fermi, c_fermi, rrms,  &
                  I_nuc,  D_nuc, Q_nuc
  Real(8)  :: ro_uniform  

  ! Type  of calculation
  ! c_speed    :: Speed of light used in the calculation
  !               (infinity for non-relativistic calculation)
   Real (kind=8)  :: c_speed =  c_au

END MODULE case
