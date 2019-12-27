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
