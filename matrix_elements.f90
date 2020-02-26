MODULE matrix_elements
  USE racah
  
  IMPLICIT NONE
  
  
  !> vmax -- maximum number of Vcoeffs in a matrix-element
  INTEGER :: vmax

! Needed subroutines
  !> onescalar  -- computes the Tcoeffs and returns corresponding ia, ib
  !> rkco    -- computes the Vcoeffs and returns ia, ib, ic, id
  REAL(kind=8), DIMENSION(:), ALLOCATABLE :: tshell, tcoeff, vcoeff
  INTEGER, DIMENSION(:), ALLOCATABLE :: tindex, vindex
  INTEGER, DIMENSION(:), ALLOCATABLE :: m_ja, m_jb, m_nt, m_nv
  REAL(kind=8), DIMENSION(:,:), ALLOCATABLE :: m_tcoeff, m_vcoeff
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: m_tindex, m_vindex
  INTEGER :: t_total, v_total  


  ! Needed functions
  !> index  -- finds the location of an integral in the integral list
  integer, external :: iindex 

  
CONTAINS
  
  SUBROUTINE genmat(outfile)

!   Generate a list of coeffisients with the following types:
!   1. Tcoeffs(ab) associated with I integrals
!   2. Vkcoeffs(ab,cd) associated with Rk integrals

    !  L o c a l  V a r i a b l e s
    INTEGER      :: blk, coeff_label, dat, i, ia, ib, ic, id, iswap, outfile, tmp
    INTEGER      :: dattmp, j, ja, jb, jprev, k, m, n, nz, nswap, nt, nv
    CHARACTER(Len=12):: convert
    
    
    ! initialize
    vmax = 4*(2*kmax+1)
    CALL allocate_racah
    CALL fact
    ALLOCATE(tshell(nw), tindex(nw), tcoeff(nw))
    ALLOCATE(coeff(vmax), label(5,vmax), vcoeff(vmax), vindex(vmax))
    ALLOCATE(m_end(nblock))
    t_total = 0; v_total = 0
    
   
    ! Define and open temporary files
    tmp = 999
    OPEN(tmp,  FORM='UNFORMATTED', STATUS='unknown')

    ! Define files for .dmel files
    dattmp = 998; dat = 500
    OPEN(dattmp,  form='formatted',status='unknown')
    OPEN(dat, FILE='matelems.dmel', form='formatted',status='unknown')

            
    jprev = 0; m = 0
    DO blk = 1, nblock
       
       !   JA and JB respectively refer to the initial and final states
       !   in the list of NCF configurations
       !
       outer:   DO jb = 1+jprev, jend(blk)
          inner:     DO ja = 1+jprev, jb           ! ja and jb loop over all matrix elements
             
             !
             !  Generate t coefficients.
             !
             nt = 0
             tshell = 0.0; tindex = 0; tcoeff=0.0
             CALL onescalar(ja,jb,ia,ib,tshell)
             IF ((ia .NE. 0) .AND. (ia .NE. ib) .AND. (ABS(tshell(1)) .GT. cutoff)) THEN
                nt = nt + 1
                IF (ia .GT. ib) THEN   !swap if necessary
                   iswap = ib
                   ib = ia
                   ia = iswap
                ENDIF
                coeff_label= ia*key + ib
                tindex(nt) = iINDEX(-1, coeff_label)     ! found in genint      
                tcoeff(nt) = tshell(1)                                
             ENDIF
             !
             !  Generate V coefficients; ac and bd are the density pairs.
             !      
             nv = 0
             vcoeff = 0.0; vindex = 0
             label = 0;coeff = 0.d0;nvcoef = 0
             CALL rkco(ja,jb,0,1)
             ! IF(label(5,1) == -1) EXIT outer  ! end of a block has been reached
!             PRINT *, 'kmax', kmax, 'nvcoef:', nvcoef
             DO j = 1, nvcoef
                IF(ABS(coeff(j)) .GT. cutoff) THEN
                   ia = label(1,j)
                   ib = label(2,j)
                   ic = label(3,j)
                   id = label(4,j)
                    k = label(5,j)
                    IF (.NOT.((k.EQ.0).AND.(ia.EQ.ic).AND.(ib.EQ.id))) THEN
                       nv = nv + 1
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
                          ! swap ia, ib 
                          iswap = ib
                          ib = ia
                          ia = iswap
                          ! swap ic, id
                          iswap = id
                          id = ic
                          ic = iswap
                       END IF
                       coeff_label = ((ia*key+ic)*key+ib)*key+id                        
                       vindex(nv) = iINDEX(k,coeff_label)
                       vcoeff(nv) = coeff(j)
                    ENDIF
                 ENDIF
              END DO
              
             ! output matrix element data to a temporary file and screen
             ! row, col, num_int(t+V), coeff(), int_index()
             
             WRITE(tmp) ja, jb, nt, nv
             IF(nt > 0) THEN
                WRITE(tmp) tcoeff(1:nt), tindex(1:nt)
             ENDIF
             IF(nv > 0) THEN
                WRITE(tmp) vcoeff(1:nv), vindex(1:nv)
             END IF

             ! output matrix element data to a temporary file and screen
             ! row, col, num_int(t+V), coeff(), int_index()
             
             WRITE(dattmp,*) ja, jb, nt, nv
             IF(nt > 0) THEN
                WRITE(dattmp,*) tcoeff(1:nt), tindex(1:nt)
             ENDIF
             IF(nv > 0) THEN
                WRITE(dattmp,*) vcoeff(1:nv), vindex(1:nv)
             END IF

             
             m = m + 1
             t_total = t_total + nt
             v_total = v_total + nv
           END DO inner
        END DO outer
        m_end(blk) = m
        jprev = jend(blk)
     END DO
     
     REWIND(tmp)
     REWIND(dattmp)
    
    !   output in binary mode  to outfile
    
    WRITE(outfile)  'MCP    '
    WRITE(outfile)  'Nblocks', nblock, 'Ncsfs  ', ncsf, 'NW     ', nw, 'ME     ', m
    WRITE(outfile)   m_end(1:nblock), t_total, v_total
    DO i = 1, m
       READ(tmp) ja, jb, nt, nv
       WRITE(outfile) ja, jb, nt, nv
       IF (nt > 0) THEN
          READ(tmp) tcoeff(1:nt), tindex(1:nt)
          WRITE(outfile) tcoeff(1:nt), tindex(1:nt)
       END IF
       IF (nv > 0) THEN
          READ(tmp) vcoeff(1:nv), vindex(1:nv)
          WRITE(outfile) vcoeff(1:nv), vindex(1:nv)
       END IF
    END DO
    CLOSE(tmp)

        !   output in decimal mode to .dat file
    
    WRITE(dat,*)  'MCP    '
    WRITE(dat,*)  'Nblocks', nblock, 'Ncsfs  ', ncsf, 'NW     ', nw, 'ME     ', m
    WRITE(dat,*)   m_end(1:nblock), t_total, v_total
    DO i = 1, m
       READ(dattmp,*) ja, jb, nt, nv
       WRITE(dat,*) ja, jb, nt, nv
       IF (nt > 0) THEN
          READ(dattmp,*) tcoeff(1:nt), tindex(1:nt)
          WRITE(dat,*) tcoeff(1:nt), tindex(1:nt)
       END IF
       IF (nv > 0) THEN
          READ(dattmp,*) vcoeff(1:nv), vindex(1:nv)
          WRITE(dat,*) vcoeff(1:nv), vindex(1:nv)
       END IF
    END DO
    CLOSE(dattmp)
    CLOSE(dat)

    WRITE(6, *) 'There are ', TRIM(convert(nblock)), ' blocks'
    WRITE(6, *) '          ', TRIM(convert(ncsf)), ' CSFs'
    WRITE(6, *) '          ', TRIM(convert(m)),' matrix elements'

    
    DEALLOCATE(coeff, label, tshell, tindex, tcoeff, vcoeff, vindex, m_end)
    CALL deallocate_racah
    
  END SUBROUTINE genmat


  SUBROUTINE readmat(infile)

     IMPLICIT NONE
     Integer, intent(in) :: infile

!   Read a list of coefficients with the following types:
!   1. Tcoeffs(ab) associated with I integrals
!   2. Vkcoeffs(ab,cd) associated with Rk integrals

    !  L o c a l  V a r i a b l e s
    INTEGER      :: i, m, nt, nv, t_pos, v_pos, t_total, v_total
    Character(Len=7) :: arg1, arg2, arg3, arg4
    CHARACTER(Len=12):: convert

    REWIND(infile)
    
    READ(infile)  arg1
    if (arg1(1:3) .ne. 'MCP')  Stop 'Not an MCP file'
    READ(infile)  arg1, nblock, arg2, ncsf, arg3, nw, arg4, m
    ALLOCATE(m_end(nblock))
    READ(infile)  m_end(1:nblock), t_total, v_total
    m = m_end(nblock)
    
    ALLOCATE(m_ja(m), m_jb(m), m_nt(m), m_nv(m))
    ALLOCATE(tcoeff(t_total), tindex(t_total), vcoeff(v_total), vindex(v_total))

    ! position markers for the t and v arrays
    t_pos = 0; v_pos= 0

    DO i = 1, m
       READ(infile) m_ja(i), m_jb(i), m_nt(i), m_nv(i)
       IF (m_nt(i) > 0) THEN
          READ(infile) tcoeff(t_pos+1:t_pos+m_nt(i)), tindex(t_pos+1:t_pos+m_nt(i))
          t_pos = t_pos + m_nt(i)
       END IF
       IF (m_nv(i) > 0) THEN
          READ(infile) vcoeff(v_pos+1:v_pos+m_nv(i)), vindex(v_pos+1:v_pos+m_nv(i))
          v_pos = v_pos + m_nv(i)
       END IF
    END DO

    WRITE(6, *) 'There are ', TRIM(convert(nblock)), ' blocks'
    WRITE(6, *) '          ', TRIM(convert(ncsf)), ' CSFs'
    WRITE(6, *) '          ', TRIM(convert(m)),' matrix elements'
    
  END SUBROUTINE readmat
  
END MODULE matrix_elements
