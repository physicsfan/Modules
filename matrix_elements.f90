MODULE matrix_elements
  USE integrals

! Needed subroutines
  !> onescalar  -- computes the Tcoeffs and returns corresponding ia, ib
  !> rkco    -- computes the Vcoeffs and returns ia, ib, ic, id
    Integer, dimension(:), allocatable :: m_ja, m_jb, m_nt, m_nv
    REAL(8), DIMENSION(:), ALLOCATABLE :: coeff, tcoeff, vcoeff
    INTEGER, DIMENSION(:), ALLOCATABLE :: tindex, vindex
    INTEGER :: t_total, v_total
  
  
CONTAINS
 
!>  Generate angular data and write to outfile
  SUBROUTINE genmat(outfile)
    implicit none

!   Generate a list of coefficients with the following types:
!   1. Tcoeffs(ab) associated with I integrals
!   2. Vkcoeffs(ab,cd) associated with Rk integrals

    !  L o c a l  V a r i a b l e s
    INTEGER   :: blk, coeff_label, i, ia, ib, ic, id, iswap, outfile, tmp
    INTEGER   :: j, ja, jb, k, m, n, nz, nswap, nt, nv
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: label

    ! initialize
    ALLOCATE(coeff(nw), tcoeff(nw), tindex(nw), vcoeff(nw), vindex(nw))
    t_total = 0; v_total = 0
    Allocate(m_end(nblock))
    ALLOCATE(label(5,nw))

    tmp = 999; !tdatfile = 100; vdatfile = 101
    OPEN(tmp,  FORM='UNFORMATTED', STATUS='unknown')


    WRITE(6, '(a,i3)' ) 'NW', nw 
    WRITE(6, '(a,i3)') "Nblocks", nblock
    WRITE(6, '(a,i3)') "Ncsfs", ncsf
 
    WRITE(6, '(//,a)') "Integral index output"
    WRITE(6, '(/,a)') " index   k   ia   ib   ic   id"

    DO blk = 1, nblock

!   JA and JB respectively refer to the initial and final states
!   in the list of NCF configurations
!
         m = 0

outer:   DO jb = 1, ncsf
inner:     DO ja = 1, jb           ! ja and jb loop over all upper  matrix elements
             
            label = 0
!
!  Generate t coefficients.
!
            nt = 0
             tcoeff = 0.0; tindex = 0
             ! onescalar is a fake function that reads in data from old rangular
             CALL onescalar(ja,jb,ia,ib,tcoeff)
             !PRINT *, tcoeff
             DO i = 1, SIZE(tcoeff)
                IF ((ia .NE. 0) .AND. (ia .NE. ib) .AND. (ABS(tcoeff(i)) .GT. cutoff)) THEN             
                   nt = nt + 1
                   IF (ia .GT. ib) THEN   !swap if necessary
                      iswap = ib
                      ib = ia
                      ia = iswap
                   ENDIF
                   coeff_label= ia*key + ib
                   tindex(i) = INDEX(-1, coeff_label)     ! found in genint
                ENDIF
             END DO
!
!  Generate V coefficients; ac and bd are the density pairs.
!      
             nv = 0
             coeff = 0.0; vcoeff = 0.0; vindex = 0
             ! rkco ia a fake function that reads in data from old rangular
             CALL rkco(ja,jb,label,coeff)
             IF(label(1,1) == -1) EXIT outer  ! end of a block has been reached
             DO j = 1, SIZE(coeff)
                IF(ABS(coeff(j)) .GT. cutoff) THEN
                   k  = label(1,j)
                   ia = label(2,j)
                   ib = label(3,j)
                   ic = label(4,j)
                   id = label(5,j)
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
                         iswap = ib
                         ib = ia
                         ia = iswap
                         iswap = id
                         id = ic
                         ic = iswap
                      END IF
                      coeff_label = ((ia*key+ic)*key+ib)*key+id                        
                      vindex(j) = INDEX(k,coeff_label)
                      vcoeff(j) = coeff(j)
                   ENDIF
                ENDIF
             END DO
             ! output or save matrix element data 
             ! row, col, num_int(t+V), coeff(), int_index()
             WRITE(tmp) ja, jb, nt, nv
             IF(nt > 0) WRITE(tmp) tcoeff(1:nt), tindex(1:nt)
             IF(nv > 0) THEN
                WRITE(tmp) vcoeff(1:nv), vindex(1:nv)
                DO i = 1, nv
                   WRITE(6, '(i5,5i5)') vindex(i), (label(j,i), j=1,5)
                END DO
             END IF
             m = m+1
             t_total = t_total + nt
             v_total = v_total + nv
          END DO inner
       END DO outer
       m_end(blk) = m
    END DO
    rewind(tmp)

!   output in binary mode  to outfile

    write(outfile)  'MCP'
    write(outfile)  'Nblocks', nblock, 'Ncsfs  ', ncsf, 'NW     ', nw, 'ME     ', m
    write(outfile)   m_end(1:nblock), t_total, v_total
    Do i = 1, m
      read(tmp) ja, jb, nt, nv
      write(outfile) ja, jb, nt, nv
      if (nt > 0) then
         read(tmp) tcoeff(1:nt), tindex(1:nt)
         write(outfile) tcoeff(1:nt), tindex(1:nt)
      end if
      if (nv > 0) then
         read(tmp) vcoeff(1:nt), vindex(1:nv)
         write(outfile) vcoeff(1:nt), vindex(1:nv)
      end if
    ENd do
    close(tmp)

    WRITE(6,'(//,a,i3)') 'Total number of t coefficients generated:', t_total
    WRITE(6,'(a,i3)') 'Total number of v coefficients generated:', v_total

    deallocate(tcoeff, tindex, vcoeff,  vindex, label)

    
  END SUBROUTINE genmat
 
  SUBROUTINE readmat(infile)

     IMPLICIT NONE
     Integer, intent(in) :: infile

!   Read a list of coefficients with the following types:
!   1. Tcoeffs(ab) associated with I integrals
!   2. Vkcoeffs(ab,cd) associated with Rk integrals

    !  L o c a l  V a r i a b l e s
    INTEGER      :: i, m, nt, nv, t_total, v_total
    Character(Len=7) :: arg1, arg2, arg3, arg4
    Character(Len=12):: convert

    
    read(infile)  arg1
    if (arg1(1:3) .ne. 'MCP')  Stop 'Not an MCP file'
    read(infile)  arg1, nblock, arg2, ncsf, arg3, nw, arg4, m
    read(infile)  m_end(1:nblock), t_total, v_total
    m = m_end(nblock)

    allocate(tcoeff(t_total), tindex(t_total), vcoeff(v_total), vindex(v_total))
    allocate(m_ja(m), m_jb(m), m_nt(m), m_nv(m))
    allocate(m_end(nblock))
    
    do i = 1, m
      read(infile) m_ja(m), m_jb(m), m_nt(m), m_nv(m)
      if (nt > 0) read(infile) tcoeff(1:nt), tindex(1:nt)
      if (nv > 0) read(infile) vcoeff(1:nt), vindex(1:nv)
    end do
 
    Write(istde, *) 'There are ', trim(convert(nblock)), ' blocks'
    write(istde, *) '          ', trim(convert(ncsf)), ' CSFs'
    write(istde, *) '          ', trim(convert(m)),' matrix elements'

   End subroutine readmat


   
   Subroutine update_wt
     implicit none
     integer :: i, j, ja, jb, m, nt, nv, it, iv, iit, iiv, blk
     real(8) :: c
     real(8), external :: wt

     ! start with the first block
     blk = 1
     ! we have processed zero coefficients
     it =0; iv  =0
     ! we have m matrix elements
     do i = 1,m
         ja = m_ja(i); jb = m_jb(i); c= wt(ja,blk)*wt(jb,blk)
         nt = m_nt(i); nv = m_nv(i)
         if ( ja .ne. jb) c = c+c
         if (nt > 0) then
            do j = 1,nt
              iit = tindex(it+j)
              int_wt(iit) = int_wt(iit) + c*tcoeff(it+j) 
            end do
            it = it +nt
          end if
          if (nv>0) then
            do j = 1,nv
              iiv = vindex(iv+j)
              int_wt(iiv) = int_wt(iiv) + c*vcoeff(iv+j)
            end do
            iv = iv + nv
          end if
          if (i == m_end(blk)) blk = blk+1
       end do
    End subroutine update_wt


 END module matrix_elements
