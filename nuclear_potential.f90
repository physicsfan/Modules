  Module nuclear_Potential
   Use radial
   Use case,  b => r_uniform, a => a_fermi, c => c_fermi

!  Real(8), dimension(npX) :: ZZ__nuc

  CONTAINS 
      Real(8) Function Z_nucleus(r)
!======================================================================
!     Nucleus potential at point "r" for different charge distributions
!     F. A. Parpia and A. K. Mohanty
!     "Relativistic basis set calculations for atoms with Fermi nuclei"
!     Phys Rev A46,3735 (1992) 
!     ZZ_nuc(r) = -r*V(r)
!----------------------------------------------------------------------
      Implicit none
      Real(8), dimension(*), intent(in) :: r
      Real(8) :: V,  N, ri, ra,rc, ac,ac2,ac3, s,s2,s3, p
!     Real(8), external :: SKFUN
      Integer :: i

      
    if (nuclear.eq.'point') then       ! point nucleus

       ZZ_nuc = Z 

    else
      do i = 1, npx
        ri = r(i)
        if(nuclear.eq.'uniform') then     ! uniform distribution
          if(ri.gt.b) then
            ZZ_nuc(i)  = Z 
          else
            ZZ_nuc(i) = -ri*(three*Z)/(two*b)*(one - ri*ri/(three*b*b))
          end if

        else if(nuclear.eq.'Fermi') then       ! Fermi distribution

          ac = a/c; ac2=ac*ac; ac3=ac2*ac
          p = pi*pi*ac2; s = six*ac3*SKFUN(3,-c/a)
          N = one + p - s
          rc = ri/c

          if(ri .le.c) then
            ra = (ri-c)/a; s2 = ac2*SKFUN(2,ra); s3 = ac3*SKFUN(3,ra)
             V = (three - rc*rc + p)/two - three*s2 - s/rc + six/rc*s3
             ZZ_nuc(i) = ri*(Z*V)/(N*c)
          else
            ra = (c-ri)/a; s2 = ac2*SKFUN(2,ra); s3 = ac3*SKFUN(3,ra)
             V = N + six*s3 + three*rc*s2
             ZZ_nuc(i) = Z * (V/N)
          end if
        else 
          Stop 'Z_nucleus: unknown charge distribution '
        end if
      end do
    end if

    End Function Z_nucleus


!=======================================================================
      Real(8) Function SKFUN (k,x)
!=======================================================================
!              infinity   (-1)^n e^nx
!      S (x) =   Sum      -----------
!       k        n=1          n^k
!======================================================================
      Implicit none
      Integer, intent(in) :: k
      Real(8), intent(in) :: x
      Integer :: n
      Real(8) :: eps, base, dnum, dk, delta

      SKFUN = 0.d0
      if(x.lt.-500d0) Return

      eps = 10.d0*epsilon(1.d0)
      BASE = -exp(x)
      DNUM = BASE
      SKFUN = BASE
      n = 1
      Do
       n = n + 1
       DNUM = DNUM*BASE
       dk = n; dk = dk ** k
       DELTA = DNUM/dk
       SKFUN = SKFUN+DELTA
       if(abs(DELTA/SKFUN).le.eps) Exit
      End do

      End Function SKFUN

   END module nuclear_potential
