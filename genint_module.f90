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
