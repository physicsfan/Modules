module Integrals

  Use radial
  implicit none

! Data structures for integral lists
  !> Int_label --  I(a,b) label of ia, ib encoded for category 1
  !>               Rk(abcd) label ia, ib, ic, id for category k+1
  !>Int_end    --  integer index array indicating the last position
  !                 for integrals in category k+2, .. j2max+2
  Integer, dimension(:), allocatable   ::  Int_label, Int_end 
  Real(8), dimension(:), allocatable   ::  Int_value,  Int_wt
  Logical, dimension(:), allocatable   ::  lval, int_varied

  ! Needed functions
  !> triangrk  -- computes the triangle relationship
  Logical, external :: triangrk 
 
  !> ykid :: key for the last computed yk(i, i, j)
  Integer :: ykid = 0

contains


  subroutine genint (mode)

!>  Generate a list of integral with the following types:
!> 1. I (or T) integrals  I(ab)
!> 2. Rk integrals for k= 0,2*kmax of Symmetry Rk(ab,cd)

!> Strategy: Each type of integral has a last value. The 
!> first category starts at i=1 whereas all others start
!> at the last value of the previous category +1

    implicit none
    integer, intent(in) :: mode
    !   L o c a l  V a r i a b l e s
    LOGICAL :: gen
    integer :: i, ia, ib, ic, id, k, n, m
    
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
                     if (mode==1)  call rinti(ia,ib,int_value(n))
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
                                 if (mode==1) call slater(k,ia,ib,ic,id, int_value(n))
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
            if (mode == 1) then
               allocate(int_Value(1:n), int_varied(1:n), int_wt(1:n))
               int_Value =0.d0; int_varied= .true.
!              when lval=.true. the Int_val is updated
            end if
            gen = .true.
            n=0
         end if
      end do

    end subroutine genint

!>  Update Int_val data
    subroutine update_int
       implicit none 
       integer :: i, ia, ib, ic, id, istart, k, lab 
     
    do i =1, Int_end(1)
       if (int_varied(i) ) then
           ib =  MOD(int_Label(i), key)
           ia = int_label(i)/key
           call rinti(ia,ib,int_value(i))
        end if
    end do  
    istart = int_end(1)+1
    do k = 0, kmax
       do i = istart, Int_end(k+2)
         if (int_varied(i) ) then
           lab = int_label(i)
           id = mod(lab, key)
           lab = lab/key
           ib = mod(lab,key)
           lab = lab/key
           ic = mod(lab,key)
           ia = lab/key
           call slater(k,ia,ib,ic,id, int_value(i))
         end if
       ENd do
       istart = Int_end(k+2)+1
     end do
    End subroutine update_int

    !> function to find the location of an integral in the list
    Integer Function index ( k, label )

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
    end function index
        
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

           
EnD Module integrals
