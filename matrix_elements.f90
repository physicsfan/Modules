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
