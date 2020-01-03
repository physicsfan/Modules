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
 
End MODULE matrix_elements
