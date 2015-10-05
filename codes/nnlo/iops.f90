! This source contains the routines to calculate the integrated
! subtraction terms: 
module IopInvariants
use momenta
implicit none
!
  logical :: init = .true.
!
  real(kind(1d0)) :: Q2
  real(kind(1d0)) , dimension(:) , allocatable :: yiQ
  real(kind(1d0)) , dimension(:,:) , allocatable :: YikQ
  type(mom) :: Q
!
contains
!
subroutine CalcIopInvariants(p)
use process
use particles
use QCDparams
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
!
  integer :: istat
  integer :: npart,ipart,kpart
  real(kind(1d0)) :: yik
!
!
! We have to allocate the arrays if no allocation was takes place
! already:
  if (init) then
    allocate(yiQ(nleg_born+1), &
             YikQ(nleg_born+1,nleg_born+1), &
             stat=istat)
    if (istat.ne.0) then
      print *,"Problem occured during allocation in IopInvariants..."
      stop
    end if
    init = .false.
  end if
!
  yiQ = -1234d0
  YikQ = -1234d0
!
  Q = p(1)%p + p(2)%p
  Q2 = Q*Q
  npart = size(p)
  do ipart=1,npart
    if (abs(p(ipart)%flv).gt.qcd_nf) cycle
    yiQ(ipart) = 2d0*p(ipart)%p*Q/Q2
  end do
!
  do ipart=1,npart
    if (abs(p(ipart)%flv).gt.qcd_nf) cycle
    do kpart=ipart+1,npart
      if (abs(p(kpart)%flv).gt.qcd_nf) cycle
      yik = 2d0*p(ipart)%p*p(kpart)%p/Q2
      YikQ(ipart,kpart) = yik/(yiQ(ipart)*yiQ(kpart))
      YikQ(kpart,ipart) = YikQ(ipart,kpart)
    end do
  end do
!
end subroutine CalcIopInvariants
!
end module IopInvariants
!
subroutine SetupQ2(Q2in)
use IopInvariants
implicit none
!
  real(kind(1d0)) , intent(in) :: Q2in
!
!
!
  Q2 = Q2in
!
end subroutine SetupQ2
!
! This is the main routine it calculates the I1 operator times the
! corresponding underlying Born SME:
subroutine CalcI1(p,Bij,CalcSMEB,CalcSMEBij,I1term,I1Laurent)
use IopInvariants
use particles
use QCDparams
use math
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
  real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
  real(kind(1d0)) , intent(out) :: I1term
  real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: I1Laurent
!
  integer :: ipart,kpart
  real(kind(1d0)) :: Ccont,Scont
  real(kind(1d0)) :: B
  real(kind(1d0)) , dimension(-4:2) :: CLaurent,SLaurent
!
  interface
    subroutine CalcSMEB(p,smeB)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: smeB
!
    end subroutine CalcSMEB
!
    subroutine CalcSMEBij(p,Bij)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Bij
!
    end subroutine CalcSMEBij
!
    subroutine CalcC1i0F(x,fi,Ccont,CLaurent)
    use particles
    implicit none
!
      real(kind(1d0)) , intent(in) :: x
      character , intent(in) :: fi
      real(kind(1d0)) , intent(out) :: Ccont
      real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: CLaurent
!
    end subroutine CalcC1i0F
!
    subroutine CalcS10ikFF(Y,Scont,SLaurent)
    use particles
    implicit none
!
      real(kind(1d0)) , intent(in) :: Y
      real(kind(1d0)) , intent(out) :: Scont
      real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: &
        SLaurent
!
    end subroutine CalcS10ikFF
  end interface
!
  I1term = 0d0
  if (present(I1Laurent)) I1Laurent = 0d0
! We calculate the invariants:
  call calcIopInvariants(p)
! We calculate the Born matrix element since it will be used
! everywhere:
  call CalcSMEB(p,B)
! We also calculate the color-correlated Born:
  call CalcSMEBij(p,Bij)
! We loop over all particles:
  do ipart=1,size(p)
! But only considering massless partons:
    if (abs(p(ipart)%flv).gt.qcd_nf) cycle
! We differentiate between two scenarios: Initial and Final
    if (ipart.le.2) then
! Initial State:
      Print *,"Initial state is not implemented yet in CalcI1..."
      stop
    else
! Final State:
! The full Laurent series is needed:
      if (present(I1Laurent)) then
! (anti)Quark:
        if (p(ipart)%flv.ne.0) then
          call CalcC1i0F(yiQ(ipart),'q',Ccont,CLaurent)
          I1term    = I1term    + 1d0/2d0/pi*Ccont*qcd_cf*B
          I1Laurent = I1Laurent + 1d0/2d0/pi*CLaurent*qcd_cf*B
! Gluon:
        else
          call CalcC1i0F(yiQ(ipart),'g',Ccont,CLaurent)
          I1term    = I1term    + 1d0/2d0/pi*Ccont*qcd_ca*B
          I1Laurent = I1Laurent + 1d0/2d0/pi*CLaurent*qcd_ca*B
        end if
! Only the finite part is needed:
      else
! (anti)Quark:
        if (p(ipart)%flv.ne.0) then
          call CalcC1i0F(yiQ(ipart),'q',Ccont)
          I1term = I1term + 1d0/2d0/pi*Ccont*qcd_cf*B
! Gluon:
        else
          call CalcC1i0F(yiQ(ipart),'g',Ccont)
          I1term = I1term + 1d0/2d0/pi*Ccont*qcd_ca*B
        end if
      end if
    end if
  end do
! We consider the contributions coming from the soft integrated
! counterterms:
  do ipart=1,size(p)
    if (abs(p(ipart)%flv).gt.qcd_nf) cycle
    do kpart=ipart+1,size(p)
      if (abs(p(kpart)%flv).gt.qcd_nf) cycle
      if (present(I1Laurent)) then
        call CalcS10ikFF(YikQ(ipart,kpart),Scont,SLaurent)
        I1term    = I1term    + 2d0/2d0/pi*Scont &
                  * Bij(ipart,kpart)
        I1Laurent = I1Laurent + 2d0/2d0/pi*SLaurent &
                  * Bij(ipart,kpart)
      else
        call CalcS10ikFF(YikQ(ipart,kpart),Scont)
        I1term = I1term + 2d0/2d0/pi*Scont &
               * Bij(ipart,kpart)
      end if
    end do
  end do
!
end subroutine CalcI1
!
! This is the same as the previous one, but this version
! takes the Born SME and the color-correlated Born SME as
! input to gain some speed:
subroutine CalcI1nlo(p,smeB,Bij,I1term,I1Laurent)
use IopInvariants
use particles
use QCDparams
use math
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
  real(kind(1d0)) , intent(in) :: smeB
  real(kind(1d0)) , dimension(:,:) , intent(in) :: Bij
  real(kind(1d0)) , intent(out) :: I1term
  real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: I1Laurent
!
  integer :: ipart,kpart,npart
  real(kind(1d0)) :: Ccont,Scont
  real(kind(1d0)) , dimension(-4:2) :: CLaurent,SLaurent
!
  interface
    subroutine CalcC1i0F(x,fi,Ccont,CLaurent)
    use particles
    implicit none
!
      real(kind(1d0)) , intent(in) :: x
      character , intent(in) :: fi
      real(kind(1d0)) , intent(out) :: Ccont
      real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: CLaurent
!
    end subroutine CalcC1i0F
!
    subroutine CalcS10ikFF(Y,Scont,SLaurent)
    use particles
    implicit none
!
      real(kind(1d0)) , intent(in) :: Y
      real(kind(1d0)) , intent(out) :: Scont
      real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: &
        SLaurent
!
    end subroutine CalcS10ikFF
  end interface
!
  npart = size(p)
!
  I1term = 0d0
  if (present(I1Laurent)) I1Laurent = 0d0
! We calculate the invariants:
  call calcIopInvariants(p)
! We loop over all particles:
  do ipart=1,npart
! But only considering massless partons:
    if (abs(p(ipart)%flv).gt.qcd_nf) cycle
! We differentiate between two scenarios: Initial and Final
    if (ipart.le.2) then
! Initial State:
      Print *,"Initial state is not implemented yet in CalcI1..."
      stop
    else
! Final State:
! The full Laurent series is needed:
      if (present(I1Laurent)) then
! (anti)Quark:
        if (p(ipart)%flv.ne.0) then
          call CalcC1i0F(yiQ(ipart),'q',Ccont,CLaurent)
          I1term    = I1term    + 1d0/2d0/pi*Ccont*qcd_cf*smeB
          I1Laurent = I1Laurent + 1d0/2d0/pi*CLaurent*qcd_cf*smeB
! Gluon:
        else
          call CalcC1i0F(yiQ(ipart),'g',Ccont,CLaurent)
          I1term    = I1term    + 1d0/2d0/pi*Ccont*qcd_ca*smeB
          I1Laurent = I1Laurent + 1d0/2d0/pi*CLaurent*qcd_ca*smeB
        end if
! Only the finite part is needed:
      else
! (anti)Quark:
        if (p(ipart)%flv.ne.0) then
          call CalcC1i0F(yiQ(ipart),'q',Ccont)
          I1term = I1term + 1d0/2d0/pi*Ccont*qcd_cf*smeB
! Gluon:
        else
          call CalcC1i0F(yiQ(ipart),'g',Ccont)
          I1term = I1term + 1d0/2d0/pi*Ccont*qcd_ca*smeB
        end if
      end if
    end if
  end do
! We consider the contributions coming from the soft integrated
! counterterms:
  do ipart=1,npart
    if (abs(p(ipart)%flv).gt.qcd_nf) cycle
    do kpart=ipart+1,npart
      if (abs(p(kpart)%flv).gt.qcd_nf) cycle
      if (present(I1Laurent)) then
        call CalcS10ikFF(YikQ(ipart,kpart),Scont,SLaurent)
        I1term    = I1term    + 2d0/2d0/pi*Scont &
                  * Bij(ipart,kpart)
        I1Laurent = I1Laurent + 2d0/2d0/pi*SLaurent &
                  * Bij(ipart,kpart)
      else
        call CalcS10ikFF(YikQ(ipart,kpart),Scont)
        I1term = I1term + 2d0/2d0/pi*Scont &
               * Bij(ipart,kpart)
      end if
    end do
  end do
!
end subroutine CalcI1nlo
!
! This corresponds to the LHS of (3.3): 
subroutine CalcC1i0F(x,fi,Ccont,CLaurent)
implicit none
!
  real(kind(1d0)) , intent(in) :: x
  character , intent(in) :: fi
  real(kind(1d0)) , intent(out) :: Ccont
  real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: CLaurent
!
  real(kind(1d0)) :: CScont
  real(kind(1d0)) , dimension(-4:2) :: CSLaurent
!
  interface
    subroutine CalcCir0qF(x,Ccont,CLaurent)
    implicit none
!
      real(kind(1d0)) , intent(in) :: x
      real(kind(1d0)) , intent(out) :: Ccont
      real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: CLaurent
!
    end subroutine CalcCir0qF
!
    subroutine CalcCir0gF(x,Ccont,CLaurent)
    implicit none
!
      real(kind(1d0)) , intent(in) :: x
      real(kind(1d0)) , intent(out) :: Ccont
      real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: CLaurent
!
    end subroutine CalcCir0gF
!
    subroutine CalcCirSr0F(CScont,CSLaurent)
    implicit none
!
      real(kind(1d0)) , intent(out) :: CScont
      real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: CSLaurent
!
    end subroutine CalcCirSr0F
  end interface
!
  if (present(CLaurent)) then
    if (fi.eq.'q') then
      call CalcCir0qF(x,Ccont,CLaurent)
    else if (fi.eq.'g') then
      call CalcCir0gF(x,Ccont,CLaurent)
    else
      print *,"Invalid flavor type in CalcC1i0F..."
      print *,"fi: ",fi
      stop
    end if
  else
    if (fi.eq.'q') then
      call CalcCir0qF(x,Ccont)
    else if (fi.eq.'g') then
      call CalcCir0gF(x,Ccont)
    else
      print *,"Invalid flavor type in CalcC1i0F..."
      print *,"fi: ",fi
      stop
    end if
  end if
!
! The CS term does not depend upon the parent flavor:
  if (present(CLaurent)) then
    call CalcCirSr0F(CScont,CSLaurent)
    Ccont    = Ccont    - CScont
    CLaurent = CLaurent - CSLaurent
  else
    call CalcCirSr0F(CScont)
    Ccont = Ccont - CScont
  end if
!
end subroutine CalcC1i0F
!
! This corresponds to the LHS of (3.4):
subroutine CalcS10ikFF(Y,Scont,SLaurent)
implicit none
!
  real(kind(1d0)) , intent(in) :: Y
  real(kind(1d0)) , intent(out) :: Scont
  real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: SLaurent
!
!
  interface
    subroutine CalcSr0ikFF(Y,Scont,SLaurent)
    implicit none
!
      real(kind(1d0)) , intent(in) :: Y
      real(kind(1d0)) , intent(out) :: Scont
      real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: &
        SLaurent
    end subroutine CalcSr0ikFF
  end interface
!
  if (present(SLaurent)) then
    call CalcSr0ikFF(Y,Scont,SLaurent)
  else
    call CalcSr0ikFF(Y,Scont)
  end if
!
end subroutine CalcS10ikFF
!
! This corresponds to the first term on the RHS of (3.3) if f_i = q
subroutine CalcCir0qF(x,Ccont,CLaurent)
use alphamax
use QCDparams
use IopInvariants
use scales
use math
implicit none
!
  real(kind(1d0)) , intent(in) :: x
  real(kind(1d0)) , intent(out) :: Ccont
  real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: CLaurent
!
  real(kind(1d0)) :: al0,mu2
!
  real(kind(1d0)) , external :: ddilog
!
  if (x.eq.-1234d0) then
    print *,"x is not initialized properly for Cir0qF..."
    print *,"x = ",x
    stop
  end if
!
  Ccont = 0d0
  if (present(CLaurent)) CLaurent = 0d0
!
  al0 = alpha0
  mu2 = mur*mur
!
! The finite piece should be always calculated:
  if (x.lt.0.9d0) then
    Ccont = &
      -pi**2/2.d0+(84+(-336-72*al0+116*al0**2-52*al0**3+ &
      9*al0**4)*x+(504+224*al0-310*al0**2+148*al0**3- &
      27*al0**4)*x**2+(-336-188*al0+276*al0**2-140*al0**3+ &
      27*al0**4)*x**3+(84+52*al0-82*al0**2+44*al0**3- &
      9*al0**4)*x**4)/(2.4d1*(-1+x)**4)+(2-2/(-1+ &
      x)**5)*ddilog((al0-al0*x)/(al0+x-al0*x))-log(al0)**2+ &
      log(mu2/Q2)**2/2.d0+(log(mu2/Q2)*(3-4*log(x)))/2.d0+ &
      2*log(x)**2+((-9+3*al0**4*(-2+x)*(-1+x)**4+96*x-192*x**2+ &
      212*x**3-116*x**4+25*x**5-4*al0**3*(-1+x)**3*(8-11*x+4*x**2) &
      +6*al0**2*(-1+x)**2*(-12+26*x-21*x**2+6*x**3)-12*al0*(-8+ &
      30*x-50*x**2+45*x**3-21*x**4+4*x**5))*log(al0+x- &
      al0*x))/(6.d0*(-1+x)**5)+(1-2/(-1+x)**5)*log(al0+x-al0*x)**2 &
      +log(x)*((18-141*x+282*x**2-302*x**3+161*x**4- &
      34*x**5)/(6.d0*(-1+x)**5)+2*(-1+(-1+x)**(-5))*log(al0+x- &
      al0*x))+log(al0)*((-9+(al0*(-3*al0**3*(-2+x)*(-1+x)**3+ &
      4*al0**2*(-1+x)**2*(8-11*x+4*x**2)+12*(8-22*x+28*x**2- &
      17*x**3+4*x**4)-6*al0*(12-38*x+47*x**2-27*x**3+6*x**4)))/(-1 &
      +x)**4)/6.d0+(2-2/(-1+x)**5)*log(x)+(2*log(al0+x-al0*x))/(-1+x)**5)
  else
    Ccont = &
      (-70*(-1+x)**4*(-58590-151200*al0+15120*al0**2+423360*al0**3 &
      -979020*al0**4+1185408*al0**5-882000*al0**6+406080*al0**7- &
      106785*al0**8+12320*al0**9+2520*(9+36*al0-144*al0**2+ &
      336*al0**3-504*al0**4+504*al0**5-336*al0**6+144*al0**7- &
      36*al0**8+4*al0**9)*log(al0)-22680*log(mu2/Q2))-3240*(-1+ &
      x)**2*(-1715-6370*al0+5880*al0**2-2450*al0**3-1225*al0**4+ &
      2058*al0**5-980*al0**6+170*al0**7+140*(7+14*al0-42*al0**2+ &
      70*al0**3-70*al0**4+42*al0**5-14*al0**6+2*al0**7)*log(al0)- &
      980*log(mu2/Q2))-17640*(-1+x)*(270+2070*al0-2025*al0**2+ &
      1300*al0**3-450*al0**4+54*al0**5+5*al0**6+60*(-6-6*al0+ &
      15*al0**2-20*al0**3+15*al0**4-6*al0**5+al0**6)*log(al0)+ &
      360*log(mu2/Q2))-756*(-1+x)**5*(4760+10640*al0+8820*al0**2- &
      73920*al0**3+173460*al0**4-239904*al0**5+217560*al0**6- &
      131520*al0**7+51345*al0**8-11760*al0**9+1204*al0**10+840*(-2 &
      -10*al0+45*al0**2-120*al0**3+210*al0**4-252*al0**5+ &
      210*al0**6-120*al0**7+45*al0**8-10*al0**9+al0**10)*log(al0)+ &
      1680*log(mu2/Q2))-315*(-1+x)**3*(15120+45360*al0- &
      29400*al0**2-27440*al0**3+88200*al0**4-96432*al0**5+ &
      56840*al0**6-18000*al0**7+2415*al0**8+840*(-8-24*al0+ &
      84*al0**2-168*al0**3+210*al0**4-168*al0**5+84*al0**6- &
      24*al0**7+3*al0**8)*log(al0)+6720*log(mu2/Q2))+3528*(3150+ &
      1350*al0-5850*al0**2+4900*al0**3-2025*al0**4+342*al0**5- &
      450*pi**2-30*(45-180*al0+60*al0**2+40*al0**3-45*al0**4+ &
      12*al0**5)*log(al0)-900*log(al0)**2+1350*log(mu2/Q2)+ &
      450*log(mu2/Q2)**2))/3.1752d6
  end if
!
  if (present(CLaurent)) then
    CLaurent(-2) = 1d0
    CLaurent(-1) = log(mu2/Q2)+(3-4*log(x))/2.d0
    CLaurent( 0) = Ccont
  end if
!
end subroutine CalcCir0qF
!
! This corresponds to the first term on the RHS of (3.3) if f_i = g
subroutine CalcCir0gF(x,Ccont,CLaurent)
use alphamax
use QCDparams
use IopInvariants
use scales
use math
implicit none
!
  real(kind(1d0)) , intent(in) :: x
  real(kind(1d0)) , intent(out) :: Ccont
  real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: CLaurent
!
  real(kind(1d0)) :: al0,mu2
  real(kind(1d0)) :: CF,CA,Nf,TR
!
  real(kind(1d0)) , external :: ddilog
!
  if (x.eq.-1234d0) then
    print *,"x is not initialized properly for Cir0gF..."
    print *,"x = ",x
    stop
  end if
!
  Ccont = 0d0
  if (present(CLaurent)) CLaurent = 0d0
!
  al0 = alpha0
  mu2 = mur*mur
!
  CF  = qcd_cf
  CA  = qcd_ca
  TR  = qcd_tr
  Nf  = qcd_nf
!
! The finite piece should be always calculated:
  if (x.lt.0.9d0) then
!    goto 100
    Ccont = &
      log(mu2/Q2)**2/2.d0+log(mu2/Q2)*((-2*Nf*TR)/(3.d0*CA)+ &
      (11d0/3d0-4*log(x))/2.d0)+(6-pi**2+(26+48*al0- &
      36*al0**2+16*al0**3-3*al0**4+(12*(-1+al0)**5*al0)/(al0*(-2+ &
      x)-x)-(960*al0)/(-2+x)**5-(240*al0**2)/(-2+x)**4- &
      (40*al0**2*(-3+2*al0))/(-2+x)**3-(10*al0**2*(6-8*al0+ &
      3*al0**2))/(-2+x)**2-(3*al0**2*(-10+20*al0-15*al0**2+ &
      4*al0**3))/(-2+x)-(12*al0)/(-1+x)**4-(6*(-2+al0)*al0)/(-1+ &
      x)**3-(4*al0*(3-3*al0+al0**2))/(-1+x)**2-(3*al0*(-4+6*al0- &
      4*al0**2+al0**3))/(-1+x))/1.8d1+ &
      (3*al0**4*x)/(4-4*x)+(10*al0*x)/(-1+x)**4+(5*al0**2*x)/(-1+ &
      x)**3+(5*al0**3*x)/(3.d0*(-1+x)**2)-(20*al0*x**2)/(-1+x)**4- &
      (15*al0**2*x**2)/(2.d0*(-1+x)**3)-(4*al0**3*x**2)/(3.d0*(-1+ &
      x)**2)+(15*al0*x**3)/(-1+x)**4+(3*al0**2*x**3)/(-1+x)**3- &
      (4*al0*x**4)/(-1+x)**4+(4*al0**3*x*(-5+4*x))/(3.d0*(-1+ &
      x)**2)-(2*al0**2*x*(10-15*x+6*x**2))/(-1+x)**3+(4*al0*x*(-10 &
      +20*x-15*x**2+4*x**3))/(-1+x)**4+(al0*x*(144-2*al0**2*(-2+ &
      x)*(-1+x)**2-248*x+176*x**2-46*x**3+al0*(-32+70*x-51*x**2+ &
      13*x**3)))/(6.d0*(-1+x)**4)+(4-4/(-1+x)**5)*ddilog((al0- &
      al0*x)/(al0+x-al0*x))-2*log(al0)**2-(11*log(x))/3.d0+ &
      (320*log(x))/(3.d0*(-2+x)**6)+(160*log(x))/(3.d0*(-2+x)**5)+ &
      (11*log(x))/(3.d0*(-1+x)**5)-(32*x*log(x))/(-1+x)**5+ &
      (64*x**2*log(x))/(-1+x)**5-(212*x**3*log(x))/(3.d0*(-1+ &
      x)**5)+(116*x**4*log(x))/(3.d0*(-1+x)**5)- &
      (25*x**5*log(x))/(3.d0*(-1+x)**5)+4*log(x)**2+((-11+ &
      3*al0**4*(-2+x)*(-1+x)**4+96*x-192*x**2+212*x**3-116*x**4+ &
      25*x**5-4*al0**3*(-1+x)**3*(8-11*x+4*x**2)+6*al0**2*(-1+ &
      x)**2*(-12+26*x-21*x**2+6*x**3)-12*al0*(-8+30*x-50*x**2+ &
      45*x**3-21*x**4+4*x**5))*log(al0+x-al0*x))/(3.d0*(-1+x)**5)- &
      4*log(x)*log(al0+x-al0*x)+(4*log(x)*log(al0+x-al0*x))/(-1+ &
      x)**5+(2-4/(-1+x)**5)*log(al0+x-al0*x)**2+log(al0)*((-11+ &
      (al0*(-3*al0**3*(-2+x)*(-1+x)**3+4*al0**2*(-1+x)**2*(8-11*x+ &
      4*x**2)+12*(8-22*x+28*x**2-17*x**3+4*x**4)-6*al0*(12-38*x+ &
      47*x**2-27*x**3+6*x**4)))/(-1+x)**4)/3.d0+(4-4/(-1+ &
      x)**5)*log(x)+(4*log(al0+x-al0*x))/(-1+x)**5)- &
      (160*x*log(2*al0+x-al0*x))/(3.d0*(-2+x)**6))/2.d0+ &
      Nf*((TR*(-60+2*(-26-48*al0+36*al0**2-16*al0**3+3*al0**4- &
      (12*(-1+al0)**5*al0)/(al0*(-2+x)-x)+(960*al0)/(-2+x)**5+ &
      (240*al0**2)/(-2+x)**4+(40*al0**2*(-3+2*al0))/(-2+x)**3+ &
      (10*al0**2*(6-8*al0+3*al0**2))/(-2+x)**2+(3*al0**2*(-10+ &
      20*al0-15*al0**2+4*al0**3))/(-2+x)+(12*al0)/(-1+x)**4+(6*(-2 &
      +al0)*al0)/(-1+x)**3+(4*al0*(3-3*al0+al0**2))/(-1+x)**2+ &
      (3*al0*(-4+6*al0-4*al0**2+al0**3))/(-1+x))+(9*al0**4*x)/(-1+ &
      x)-(12*al0**3*x*(-5+4*x))/(-1+x)**2+(18*al0**2*x*(10-15*x+ &
      6*x**2))/(-1+x)**3-(36*al0*x*(-10+20*x-15*x**2+4*x**3))/(-1+ &
      x)**4-3*(-24+(3*al0**4*x)/(-1+x)-(4*al0**3*x*(-5+4*x))/(-1+ &
      x)**2+(6*al0**2*x*(10-15*x+6*x**2))/(-1+x)**3-(12*al0*x*(-10 &
      +20*x-15*x**2+4*x**3))/(-1+ &
      x)**4)))/(3.6d1*CA)+(2*TR*log(al0))/(3.d0*CA) &
      +(TR*(2-320/(-2+x)**6-160/(-2+x)**5-2/(-1+ &
      x)**5)*log(x))/(3.d0*CA)+(2*TR*log(al0+x- &
      al0*x))/(3.d0*CA*(-1+x)**5)+(160*TR*x*log(2*al0+x- &
      al0*x))/(3.d0*CA*(-2+x)**6))
      goto 999
100 continue
! Only qq:
    Ccont = &
      (TR*(-60+2*(-26-48*al0+36*al0**2-16*al0**3+3*al0**4-(12*(-1+ &
      al0)**5*al0)/(al0*(-2+x)-x)+(960*al0)/(-2+x)**5+ &
      (240*al0**2)/(-2+x)**4+(40*al0**2*(-3+2*al0))/(-2+x)**3+ &
      (10*al0**2*(6-8*al0+3*al0**2))/(-2+x)**2+(3*al0**2*(-10+ &
      20*al0-15*al0**2+4*al0**3))/(-2+x)+(12*al0)/(-1+x)**4+(6*(-2 &
      +al0)*al0)/(-1+x)**3+(4*al0*(3-3*al0+al0**2))/(-1+x)**2+ &
      (3*al0*(-4+6*al0-4*al0**2+al0**3))/(-1+x))+(9*al0**4*x)/(-1+ &
      x)-(12*al0**3*x*(-5+4*x))/(-1+x)**2+(18*al0**2*x*(10-15*x+ &
      6*x**2))/(-1+x)**3-(36*al0*x*(-10+20*x-15*x**2+4*x**3))/(-1+ &
      x)**4-3*(-24+(3*al0**4*x)/(-1+x)-(4*al0**3*x*(-5+4*x))/(-1+ &
      x)**2+(6*al0**2*x*(10-15*x+6*x**2))/(-1+x)**3-(12*al0*x*(-10 &
      +20*x-15*x**2+4*x**3))/(-1+ &
      x)**4)))/(3.6d1*CA)+(2*TR*log(al0))/(3.d0*CA) &
      -(2*TR*log(mu2/Q2))/(3.d0*CA)+(TR*(2-320/(-2+x)**6-160/(-2+ &
      x)**5-2/(-1+x)**5)*log(x))/(3.d0*CA)+(2*TR*log(al0+x- &
      al0*x))/(3.d0*CA*(-1+x)**5)+(160*TR*x*log(2*al0+x- &
      al0*x))/(3.d0*CA*(-2+x)**6)
      goto 999
110 continue
! Only gg:
    Ccont = &
      6-pi**2+(26+48*al0-36*al0**2+16*al0**3-3*al0**4+(12*(-1+ &
      al0)**5*al0)/(al0*(-2+x)-x)-(960*al0)/(-2+x)**5- &
      (240*al0**2)/(-2+x)**4-(40*al0**2*(-3+2*al0))/(-2+x)**3- &
      (10*al0**2*(6-8*al0+3*al0**2))/(-2+x)**2-(3*al0**2*(-10+ &
      20*al0-15*al0**2+4*al0**3))/(-2+x)-(12*al0)/(-1+x)**4-(6*(-2 &
      +al0)*al0)/(-1+x)**3-(4*al0*(3-3*al0+al0**2))/(-1+x)**2- &
      (3*al0*(-4+6*al0-4*al0**2+al0**3))/(-1+ &
      x))/1.8d1+(3*al0**4*x)/(4-4*x)+(10*al0*x)/(-1 &
      +x)**4+(5*al0**2*x)/(-1+x)**3+(5*al0**3*x)/(3.d0*(-1+x)**2)- &
      (20*al0*x**2)/(-1+x)**4-(15*al0**2*x**2)/(2.d0*(-1+x)**3)- &
      (4*al0**3*x**2)/(3.d0*(-1+x)**2)+(15*al0*x**3)/(-1+x)**4+ &
      (3*al0**2*x**3)/(-1+x)**3-(4*al0*x**4)/(-1+x)**4+ &
      (4*al0**3*x*(-5+4*x))/(3.d0*(-1+x)**2)-(2*al0**2*x*(10-15*x+ &
      6*x**2))/(-1+x)**3+(4*al0*x*(-10+20*x-15*x**2+4*x**3))/(-1+ &
      x)**4+(al0*x*(144-2*al0**2*(-2+x)*(-1+x)**2-248*x+176*x**2- &
      46*x**3+al0*(-32+70*x-51*x**2+13*x**3)))/(6.d0*(-1+x)**4)+(4 &
      -4/(-1+x)**5)*ddilog((al0-al0*x)/(al0+x-al0*x))- &
      2*log(al0)**2+log(mu2/Q2)**2+ &
      log(mu2/Q2)*(11d0/3d0-4*log(x))-(11*log(x))/3.d0 &
      +(320*log(x))/(3.d0*(-2+x)**6)+(160*log(x))/(3.d0*(-2+x)**5) &
      +(11*log(x))/(3.d0*(-1+x)**5)-(32*x*log(x))/(-1+x)**5+ &
      (64*x**2*log(x))/(-1+x)**5-(212*x**3*log(x))/(3.d0*(-1+ &
      x)**5)+(116*x**4*log(x))/(3.d0*(-1+x)**5)- &
      (25*x**5*log(x))/(3.d0*(-1+x)**5)+4*log(x)**2+((-11+ &
      3*al0**4*(-2+x)*(-1+x)**4+96*x-192*x**2+212*x**3-116*x**4+ &
      25*x**5-4*al0**3*(-1+x)**3*(8-11*x+4*x**2)+6*al0**2*(-1+ &
      x)**2*(-12+26*x-21*x**2+6*x**3)-12*al0*(-8+30*x-50*x**2+ &
      45*x**3-21*x**4+4*x**5))*log(al0+x-al0*x))/(3.d0*(-1+x)**5)- &
      4*log(x)*log(al0+x-al0*x)+(4*log(x)*log(al0+x-al0*x))/(-1+ &
      x)**5+(2-4/(-1+x)**5)*log(al0+x-al0*x)**2+log(al0)*((-11+ &
      (al0*(-3*al0**3*(-2+x)*(-1+x)**3+4*al0**2*(-1+x)**2*(8-11*x+ &
      4*x**2)+12*(8-22*x+28*x**2-17*x**3+4*x**4)-6*al0*(12-38*x+ &
      47*x**2-27*x**3+6*x**4)))/(-1+x)**4)/3.d0+(4-4/(-1+ &
      x)**5)*log(x)+(4*log(al0+x-al0*x))/(-1+x)**5)- &
      (160*x*log(2*al0+x-al0*x))/(3.d0*(-2+x)**6)
      Ccont = Ccont / 2d0
      goto 999
  else
!    goto 101
    Ccont = &
      -(Nf*TR*(50+2630*al0+1125*al0**2-325*al0**3+125*al0**4- &
      39*al0**5+6*al0**6-30*(1+al0)*log(al0)-2400*(1+al0)*log(1+ &
      al0)))/(4.5d1*(1+al0)*CA)-(-3350-30500*al0-6750*al0**2+ &
      4200*al0**3-4125*al0**4+2073*al0**5-402*al0**6+450*pi**2+ &
      450*al0*pi**2+30*(55-125*al0-120*al0**2+100*al0**3-5*al0**4- &
      33*al0**5+12*al0**6)*log(al0)+900*(1+al0)*log(al0)**2+ &
      24000*(1+al0)*log(1+al0))/(9.d2*(1+al0))+((11*CA- &
      4*Nf*TR)*log(mu2/Q2))/(6.d0*CA)+log(mu2/Q2)**2/2.d0-((-1+ &
      x)*(330*CA-30750*al0*CA-48195*al0**2*CA-11420*al0**3*CA+ &
      2300*al0**4*CA+4*al0**5*CA-332*al0**6*CA+104*al0**7*CA- &
      5*al0**8*CA-120*Nf*TR+66720*al0*Nf*TR+101160*al0**2*Nf*TR+ &
      21480*al0**3*Nf*TR-4350*al0**4*Nf*TR+900*al0**5*Nf*TR- &
      10*al0**6*Nf*TR-80*al0**7*Nf*TR+20*al0**8*Nf*TR+60*(1+ &
      al0)**2*(-6-6*al0+15*al0**2-20*al0**3+15*al0**4-6*al0**5+ &
      al0**6)*CA*log(al0)+33600*(1+al0)**2*(CA-2*Nf*TR)*log(1+al0) &
      +360*CA*log(mu2/Q2)+720*al0*CA*log(mu2/Q2)+ &
      360*al0**2*CA*log(mu2/Q2)))/(1.8d2*(1+ &
      al0)**2*CA)-((-1+x)**2*(-5635*CA-2154775*al0*CA- &
      5344185*al0**2*CA-3910445*al0**3*CA-499800*al0**4*CA+ &
      69384*al0**5*CA-13818*al0**6*CA+6642*al0**7*CA- &
      1536*al0**8*CA-850*al0**9*CA+370*al0**10*CA+980*Nf*TR+ &
      4240460*al0*Nf*TR+10578120*al0**2*Nf*TR+7787080*al0**3*Nf*TR &
      +1015770*al0**4*Nf*TR-157290*al0**5*Nf*TR+22050*al0**6*Nf*TR &
      -210*al0**7*Nf*TR+840*al0**8*Nf*TR-1120*al0**9*Nf*TR+ &
      280*al0**10*Nf*TR+420*(1+al0)**3*(7+14*al0-42*al0**2+ &
      70*al0**3-70*al0**4+42*al0**5-14*al0**6+ &
      2*al0**7)*CA*log(al0)+2116800*(1+al0)**3*(CA-2*Nf*TR)*log(1+ &
      al0)-2940*CA*log(mu2/Q2)-8820*al0*CA*log(mu2/Q2)- &
      8820*al0**2*CA*log(mu2/Q2)- &
      2940*al0**3*CA*log(mu2/Q2)))/(2.94d3*(1+al0)**3*CA)-((-1+ &
      x)**3*(16240*CA-20580560*al0*CA-72208920*al0**2*CA- &
      89434800*al0**3*CA-43249080*al0**4*CA-3928512*al0**5*CA+ &
      558992*al0**6*CA-112192*al0**7*CA-13053*al0**8*CA+ &
      27548*al0**9*CA-1510*al0**10*CA-6660*al0**11*CA+ &
      1995*al0**12*CA-2240*Nf*TR+41372800*al0*Nf*TR+ &
      144903360*al0**2*Nf*TR+179244800*al0**3*Nf*TR+ &
      86495360*al0**4*Nf*TR+7896000*al0**5*Nf*TR- &
      995680*al0**6*Nf*TR+136640*al0**7*Nf*TR-26040*al0**8*Nf*TR+ &
      10080*al0**9*Nf*TR+1680*al0**10*Nf*TR-3360*al0**11*Nf*TR+ &
      840*al0**12*Nf*TR+840*(1+al0)**4*(-8-24*al0+84*al0**2- &
      168*al0**3+210*al0**4-168*al0**5+84*al0**6-24*al0**7+ &
      3*al0**8)*CA*log(al0)+20697600*(1+al0)**4*(CA-2*Nf*TR)*log(1 &
      +al0)+6720*CA*log(mu2/Q2)+26880*al0*CA*log(mu2/Q2)+ &
      40320*al0**2*CA*log(mu2/Q2)+26880*al0**3*CA*log(mu2/Q2)+ &
      6720*al0**4*CA*log(mu2/Q2)))/(1.008d4*(1+al0)**4*CA)-((-1+ &
      x)**4*(-62370*CA-220640490*al0*CA-991943820*al0**2*CA- &
      1726518780*al0**3*CA-1412225010*al0**4*CA- &
      504391482*al0**5*CA-35474040*al0**6*CA+4037400*al0**7*CA- &
      59625*al0**8*CA-265165*al0**9*CA-138122*al0**10*CA+ &
      176870*al0**11*CA-6325*al0**12*CA-38465*al0**13*CA+ &
      10640*al0**14*CA+7560*Nf*TR+440392680*al0*Nf*TR+ &
      1981234080d0*al0**2*Nf*TR+3449839680d0*al0**3*Nf*TR+ &
      2823418080d0*al0**4*Nf*TR+1008504000d0*al0**5*Nf*TR+ &
      69773760*al0**6*Nf*TR-7691040*al0**7*Nf*TR+ &
      1091160*al0**8*Nf*TR-135240*al0**9*Nf*TR-63840*al0**10*Nf*TR &
      +53760*al0**11*Nf*TR+3360*al0**12*Nf*TR-13440*al0**13*Nf*TR+ &
      3360*al0**14*Nf*TR+2520*(1+al0)**5*(9+36*al0-144*al0**2+ &
      336*al0**3-504*al0**4+504*al0**5-336*al0**6+144*al0**7- &
      36*al0**8+4*al0**9)*CA*log(al0)+220147200*(1+al0)**5*(CA- &
      2*Nf*TR)*log(1+al0)-22680*CA*log(mu2/Q2)- &
      113400*al0*CA*log(mu2/Q2)-226800*al0**2*CA*log(mu2/Q2)- &
      226800*al0**3*CA*log(mu2/Q2)-113400*al0**4*CA*log(mu2/Q2)- &
      22680*al0**5*CA*log(mu2/Q2)))/(4.536d4*(1+al0)**5*CA)-((-1+ &
      x)**5*(5040*CA-42292320*al0*CA-232708140*al0**2*CA- &
      521855600*al0**3*CA-603263710*al0**4*CA-368017524*al0**5*CA- &
      103941754*al0**6*CA-5757560*al0**7*CA+458595*al0**8*CA- &
      19310*al0**9*CA+62745*al0**10*CA-44580*al0**11*CA- &
      15725*al0**12*CA+21430*al0**13*CA-1155*al0**14*CA- &
      3976*al0**15*CA+1064*al0**16*CA-560*Nf*TR+84663040*al0*Nf*TR &
      +465704400*al0**2*Nf*TR+1044178800*al0**3*Nf*TR+ &
      1206820300*al0**4*Nf*TR+736148280*al0**5*Nf*TR+ &
      208088580*al0**6*Nf*TR+11513600*al0**7*Nf*TR- &
      1127700*al0**8*Nf*TR+124600*al0**9*Nf*TR+10220*al0**10*Nf*TR &
      -10080*al0**11*Nf*TR-5600*al0**12*Nf*TR+5600*al0**13*Nf*TR- &
      1120*al0**15*Nf*TR+280*al0**16*Nf*TR+840*(1+al0)**6*(-2- &
      10*al0+45*al0**2-120*al0**3+210*al0**4-252*al0**5+210*al0**6 &
      -120*al0**7+45*al0**8-10*al0**9+al0**10)*CA*log(al0)+ &
      42336000*(1+al0)**6*(CA-2*Nf*TR)*log(1+al0)+ &
      1680*CA*log(mu2/Q2)+10080*al0*CA*log(mu2/Q2)+ &
      25200*al0**2*CA*log(mu2/Q2)+33600*al0**3*CA*log(mu2/Q2)+ &
      25200*al0**4*CA*log(mu2/Q2)+10080*al0**5*CA*log(mu2/Q2)+ &
      1680*al0**6*CA*log(mu2/Q2)))/(4.2d3*(1+al0)**6*CA)
      goto 999
101 continue 
! only qq:
    Ccont = &
      (TR*(-210*(1+al0)**4*(-1+x)*(-12+6672*al0+10116*al0**2+ &
      2148*al0**3-435*al0**4+90*al0**5-al0**6-8*al0**7+2*al0**8- &
      6720*(1+al0)**2*log(1+al0))-90*(1+al0)**3*(-1+x)**2*(14+ &
      60578*al0+151116*al0**2+111244*al0**3+14511*al0**4- &
      2247*al0**5+315*al0**6-3*al0**7+12*al0**8-16*al0**9+ &
      4*al0**10-60480*(1+al0)**3*log(1+al0))-105*(1+al0)**2*(-1+ &
      x)**3*(-8+147760*al0+517512*al0**2+640160*al0**3+ &
      308912*al0**4+28200*al0**5-3556*al0**6+488*al0**7-93*al0**8+ &
      36*al0**9+6*al0**10-12*al0**11+3*al0**12-147840*(1+ &
      al0)**4*log(1+al0))-70*(1+al0)*(-1+x)**4*(9+524277*al0+ &
      2358612*al0**2+4106952*al0**3+3361212*al0**4+1200600*al0**5+ &
      83064*al0**6-9156*al0**7+1299*al0**8-161*al0**9-76*al0**10+ &
      64*al0**11+4*al0**12-16*al0**13+4*al0**14-524160*(1+ &
      al0)**5*log(1+al0))-126*(-1+x)**5*(-4+604736*al0+ &
      3326460*al0**2+7458420*al0**3+8620145*al0**4+5258202*al0**5+ &
      1486347*al0**6+82240*al0**7-8055*al0**8+890*al0**9+ &
      73*al0**10-72*al0**11-40*al0**12+40*al0**13-8*al0**15+ &
      2*al0**16-604800*(1+al0)**6*log(1+al0))-84*(1+al0)**5*(50+ &
      2630*al0+1125*al0**2-325*al0**3+125*al0**4-39*al0**5+ &
      6*al0**6-30*(1+al0)*log(al0)-2400*(1+al0)*log(1+al0)+ &
      30*log(mu2/Q2)+ &
      30*al0*log(mu2/Q2))))/(3.78d3*(1+al0)**6*CA)
      goto 999
102 continue
! Only gg:
    Ccont = &
      -(17640*(1+al0)**4*(-1+x)*(330-30750*al0-48195*al0**2- &
      11420*al0**3+2300*al0**4+4*al0**5-332*al0**6+104*al0**7- &
      5*al0**8+60*(1+al0)**2*(-6-6*al0+15*al0**2-20*al0**3+ &
      15*al0**4-6*al0**5+al0**6)*log(al0)+33600*(1+al0)**2*log(1+ &
      al0)+360*log(mu2/Q2)+720*al0*log(mu2/Q2)+ &
      360*al0**2*log(mu2/Q2))+1080*(1+al0)**3*(-1+x)**2*(-5635- &
      2154775*al0-5344185*al0**2-3910445*al0**3-499800*al0**4+ &
      69384*al0**5-13818*al0**6+6642*al0**7-1536*al0**8-850*al0**9 &
      +370*al0**10+420*(1+al0)**3*(7+14*al0-42*al0**2+70*al0**3- &
      70*al0**4+42*al0**5-14*al0**6+2*al0**7)*log(al0)+2116800*(1+ &
      al0)**3*log(1+al0)-2940*log(mu2/Q2)-8820*al0*log(mu2/Q2)- &
      8820*al0**2*log(mu2/Q2)-2940*al0**3*log(mu2/Q2))+315*(1+ &
      al0)**2*(-1+x)**3*(16240-20580560*al0-72208920*al0**2- &
      89434800*al0**3-43249080*al0**4-3928512*al0**5+558992*al0**6 &
      -112192*al0**7-13053*al0**8+27548*al0**9-1510*al0**10- &
      6660*al0**11+1995*al0**12+840*(1+al0)**4*(-8-24*al0+ &
      84*al0**2-168*al0**3+210*al0**4-168*al0**5+84*al0**6- &
      24*al0**7+3*al0**8)*log(al0)+20697600*(1+al0)**4*log(1+al0)+ &
      6720*log(mu2/Q2)+26880*al0*log(mu2/Q2)+ &
      40320*al0**2*log(mu2/Q2)+26880*al0**3*log(mu2/Q2)+ &
      6720*al0**4*log(mu2/Q2))-70*(1+al0)*(-1+x)**4*(62370+ &
      220640490*al0+991943820*al0**2+1726518780*al0**3+ &
      1412225010*al0**4+504391482*al0**5+35474040*al0**6- &
      4037400*al0**7+59625*al0**8+265165*al0**9+138122*al0**10- &
      176870*al0**11+6325*al0**12+38465*al0**13-10640*al0**14- &
      2520*(1+al0)**5*(9+36*al0-144*al0**2+336*al0**3-504*al0**4+ &
      504*al0**5-336*al0**6+144*al0**7-36*al0**8+ &
      4*al0**9)*log(al0)-220147200*(1+al0)**5*log(1+al0)+ &
      22680*log(mu2/Q2)+113400*al0*log(mu2/Q2)+ &
      226800*al0**2*log(mu2/Q2)+226800*al0**3*log(mu2/Q2)+ &
      113400*al0**4*log(mu2/Q2)+22680*al0**5*log(mu2/Q2))+756*(-1+ &
      x)**5*(5040-42292320*al0-232708140*al0**2-521855600*al0**3- &
      603263710*al0**4-368017524*al0**5-103941754*al0**6- &
      5757560*al0**7+458595*al0**8-19310*al0**9+62745*al0**10- &
      44580*al0**11-15725*al0**12+21430*al0**13-1155*al0**14- &
      3976*al0**15+1064*al0**16+840*(1+al0)**6*(-2-10*al0+ &
      45*al0**2-120*al0**3+210*al0**4-252*al0**5+210*al0**6- &
      120*al0**7+45*al0**8-10*al0**9+al0**10)*log(al0)+42336000*(1 &
      +al0)**6*log(1+al0)+1680*log(mu2/Q2)+10080*al0*log(mu2/Q2)+ &
      25200*al0**2*log(mu2/Q2)+33600*al0**3*log(mu2/Q2)+ &
      25200*al0**4*log(mu2/Q2)+10080*al0**5*log(mu2/Q2)+ &
      1680*al0**6*log(mu2/Q2))-3528*(1+al0)**5*(3350+30500*al0+ &
      6750*al0**2-4200*al0**3+4125*al0**4-2073*al0**5+402*al0**6- &
      450*pi**2-450*al0*pi**2-30*(55-125*al0-120*al0**2+100*al0**3 &
      -5*al0**4-33*al0**5+12*al0**6)*log(al0)-900*(1+ &
      al0)*log(al0)**2-24000*(1+al0)*log(1+al0)+1650*log(mu2/Q2)+ &
      1650*al0*log(mu2/Q2)+450*log(mu2/Q2)**2+ &
      450*al0*log(mu2/Q2)**2))/(1.5876d6*(1+al0)**6)
      Ccont = Ccont / 2d0
      goto 999
  end if
!
999 continue
!
  if (present(CLaurent)) then
    CLaurent(-2) = 1d0
    CLaurent(-1) =  &
      (-2*qcd_nf*qcd_tr)/(3.d0*qcd_ca)+log(mu2/Q2)+(11d0/3.d0 - &
      4*log(x))/2.d0
    CLaurent( 0) = Ccont
  end if
!
end subroutine CalcCir0gF
!
! This corresponds to the second term on the RHS of (3.3) if f_i = q
! also note that it is independent of the flavor hence that index
! is omitted and this routine is used for both of the contributions:
subroutine CalcCirSr0F(CScont,CSLaurent)
use alphamax
use QCDparams
use IopInvariants
use scales
use math
implicit none
!
  real(kind(1d0)) , intent(out) :: CScont
  real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: CSLaurent
!
  real(kind(1d0)) :: al0,mu2
!
  real(kind(1d0)) , external :: ddilog
!
  CScont = 0d0
  if (present(CSLaurent)) CSLaurent = 0d0
!
  al0 = alpha0
  mu2 = mur*mur
!
! The finite piece should be always calculated:
  if (y0.ne.1d0) then
    CScont = &
      log(mu2/Q2)**2/2.d0+log(mu2/Q2)*(6*y0-3*y0**2+(2*y0**3)/3.d0 &
      -2*log(y0))+(-3*pi**2+450*y0-153*y0**2+32*y0**3- &
      108*ddilog(y0)+198*log(1-y0)-324*y0*log(1-y0)+ &
      162*y0**2*log(1-y0)-36*y0**3*log(1-y0)-216*y0*log(y0)+ &
      108*y0**2*log(y0)-24*y0**3*log(y0)+ &
      36*log(y0)**2)/1.8d1
  else
    CScont = &
      (329-21*pi**2+66*log(mu2/Q2)+  &
      9*log(mu2/Q2)**2)/1.8d1
  end if
!
  if (present(CSLaurent)) then
    CSLaurent(-2) = 1d0
    CSLaurent(-1) = 6*y0-3*y0**2+(2*y0**3)/3.d0+log(mu2/Q2)-2*log(y0)
    CSLaurent( 0) = CScont
  end if
!
end subroutine CalcCirSr0F
!
! This contribution corresponds to (3.4):
subroutine CalcSr0ikFF(Y,Scont,SLaurent)
use alphamax
use QCDparams
use IopInvariants
use scales
use math
implicit none
!
  real(kind(1d0)) , intent(in) :: Y
  real(kind(1d0)) , intent(out) :: Scont
  real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: &
    SLaurent
!
  real(kind(1d0)) :: mu2
!
  real(kind(1d0)) , external :: ddilog
!
  Scont = 0d0
  if (present(SLaurent)) SLaurent = 0d0
!
  mu2 = mur*mur
!
  if (Y.eq.-1234d0) then
    print *,"Y is not initialized properly for Sr0ik..."
    print *,"Y = ",Y
    stop
  end if
!
  if (y0.ne.1d0) then
    Scont = &
      -log(mu2/Q2)**2/2.d0+log(mu2/Q2)*(-6*y0+3*y0**2- &
      (2*y0**3)/3.d0+log(Y)+2*log(y0))+(3*pi**2-414*y0+117*y0**2- &
      20*y0**3-18*ddilog(1-Y)+108*ddilog(y0)+108*y0*log(Y)- &
      54*y0**2*log(Y)+12*y0**3*log(Y)-9*log(Y)**2-198*log(1-y0)+ &
      324*y0*log(1-y0)-162*y0**2*log(1-y0)+36*y0**3*log(1-y0)+ &
      216*y0*log(y0)-108*y0**2*log(y0)+24*y0**3*log(y0)- &
      36*log(Y)*log(y0)-36*log(y0)**2)/1.8d1
  else
    Scont = &
      (-317+21*pi**2-18*ddilog(1-Y)-66*log(mu2/Q2)-  &
      9*log(mu2/Q2)**2+66*log(Y)+18*log(mu2/Q2)*log(Y)-  &
      9*log(Y)**2)/1.8d1
  end if
!
  if (present(SLaurent)) then
    SLaurent(-2) = -1d0
    SLaurent(-1) = -6*y0+3*y0**2-(2*y0**3)/3.d0-log(mu2/Q2)+log(Y)+2*log(y0)
    SLaurent( 0) = Scont
  end if
!
end subroutine CalcSr0ikFF
!
! New versions for C and S routines containing terms up to 
! O(ep^2):

!*****C1q(x) routine*****
subroutine CalcC1qF(x,laurent)
  use math  
  implicit none

!parameterlist
  real(kind(1d0)), intent(in) :: x
  real(kind(1d0)), dimension(-4:2), intent(out) :: laurent

!variables
  real(kind(1d0)) :: logX, log1mX
  real(kind(1d0)), dimension(2:4) :: LiN
  real(kind(1d0)) :: Li2arg1, Li3arg1, Li4arg1
  real(kind(1d0)) :: Li2arg2, Li3arg2, Li4arg2
  real(kind(1d0)) :: Li4arg3

!Calculate logs and Li's once for all
  logX = log(x)
  log1mX = log(1-x)

!arg1: x
  call PolyLogs(x,LiN)

  Li2arg1 = LiN(2)
  Li3arg1 = LiN(3)
  Li4arg1 = LiN(4)

!arg2: 1 - x
  call PolyLogs(1d0 - x,LiN)

  Li2arg2 = LiN(2)
  Li3arg2 = LiN(3)
  Li4arg2 = LiN(4)

!arg3: x/(x - 1)
  call PolyLogs(x/(x - 1d0),LiN)

  Li4arg3 = LiN(4)

  laurent = 0

!Calculate series
  laurent(-2) = 0

  laurent(-1) = -13d0/6d0 - 2*logX

!In case x<=0.9 we can use the original formula
  if (x.lt.0.9d0) then
   laurent(0) = (1064 - 96*pi**2 + 144*logX**2*(-1 + x)**5 - 5323*x + &
    360*pi**2*x + 10538*x**2 - 720*pi**2*x**2 - 10460*x**3 + &
    720*pi**2*x**3 + 5230*x**4 - 360*pi**2*x**4 - &
    1049*x**5 + 72*pi**2*x**5 - &
    12*logX*(-18 + 141*x - 282*x**2 + 302*x**3 - 161*x**4 + &
       34*x**5 + 12*log1mX*&
        (-2 + 5*x - 10*x**2 + 10*x**3 - 5*x**4 + x**5)) - &
    144*(-2 + 5*x - 10*x**2 + 10*x**3 - 5*x**4 + x**5)*&
     Li2arg1)/(72.*(-1 + x)**5) 

   laurent(1) = (36092 + 3024*logX - 9216*log1mX*logX - 1296*logX**2 - &
    1728*log1mX*logX**2 + 576*logX**3 - 1956*pi**2 + &
    144*logX*pi**2 - 186952*x - 49608*logX*x + &
    12024*log1mX*logX*x + 10152*logX**2*x + &
    4320*log1mX*logX**2*x - 2880*logX**3*x + 11760*pi**2*x + &
    720*logX*pi**2*x + 373598*x**2 + 100512*logX*x**2 - &
    6768*log1mX*logX*x**2 - 20304*logX**2*x**2 - &
    8640*log1mX*logX**2*x**2 + 5760*logX**3*x**2 - &
    23520*pi**2*x**2 - 1440*logX*pi**2*x**2 - 372755*x**3 - &
    114648*logX*x**3 - 6192*log1mX*logX*x**3 + &
    21744*logX**2*x**3 + 8640*log1mX*logX**2*x**3 - &
    5760*logX**3*x**3 + 24720*pi**2*x**3 + &
    1440*logX*pi**2*x**3 + 188074*x**4 + 63924*logX*x**4 + &
    8136*log1mX*logX*x**4 - 11592*logX**2*x**4 - &
    4320*log1mX*logX**2*x**4 + 2880*logX**3*x**4 - &
    12960*pi**2*x**4 - 720*logX*pi**2*x**4 - 38057*x**5 - &
    13962*logX*x**5 - 2448*log1mX*logX*x**5 + &
    2448*logX**2*x**5 + 864*log1mX*logX**2*x**5 - &
    576*logX**3*x**5 + 2700*pi**2*x**5 + &
    144*logX*pi**2*x**5 - &
    72*(128 - 167*x + 94*x**2 + 86*x**3 - 113*x**4 + &
       34*x**5 + 48*logX*&
        (-2 + 5*x - 10*x**2 + 10*x**3 - 5*x**4 + x**5))*&
     Li2arg1 + 864*&
     (-18 + 25*x - 50*x**2 + 50*x**3 - 25*x**4 + 5*x**5)*&
     Li3arg2 - 17280*Li3arg1 + &
    43200*x*Li3arg1 - 86400*x**2*Li3arg1 + &
    86400*x**3*Li3arg1 - 43200*x**4*Li3arg1 + &
    8640*x**5*Li3arg1 + 3456*Zeta3 + &
    25920*x*Zeta3 - 51840*x**2*Zeta3 + &
    51840*x**3*Zeta3 - 25920*x**4*Zeta3 + &
    5184*x**5*Zeta3)/(432.*(-1 + x)**5) 

    laurent(2) = -(-2819570 + 23760*log1mX**4 - 90720*logX + &
     1604160*log1mX*logX - 95040*log1mX**3*logX + &
     45360*logX**2 - 138240*log1mX*logX**2 - &
     90720*log1mX**2*logX**2 - 12960*logX**3 - &
     17280*log1mX*logX**3 + 4320*logX**4 + 141960*pi**2 + &
     47520*log1mX**2*pi**2 - 3240*logX*pi**2 + &
     47520*log1mX*logX*pi**2 + 2160*logX**2*pi**2 - &
     9792*pi**4 + 15504760*x + 48600*log1mX**4*x + &
     4535280*logX*x - 2456820*log1mX*logX*x - &
     194400*log1mX**3*logX*x - 744120*logX**2*x + &
     180360*log1mX*logX**2*x + 1393200*log1mX**2*logX**2*x + &
     101520*logX**3*x + 43200*log1mX*logX**3*x - &
     21600*logX**4*x - 1103685*pi**2*x + &
     97200*log1mX**2*pi**2*x - 43740*logX*pi**2*x - &
     550800*log1mX*logX*pi**2*x + 10800*logX**2*pi**2*x + &
     10440*pi**4*x - 31504880*x**2 - 97200*log1mX**4*x**2 - &
     9079740*logX*x**2 + 2382120*log1mX*logX*x**2 + &
     388800*log1mX**3*logX*x**2 + 1507680*logX**2*x**2 - &
     101520*log1mX*logX**2*x**2 - &
     2786400*log1mX**2*logX**2*x**2 - 203040*logX**3*x**2 - &
     86400*log1mX*logX**3*x**2 + 43200*logX**4*x**2 + &
     2208990*pi**2*x**2 - 194400*log1mX**2*pi**2*x**2 + &
     87480*logX*pi**2*x**2 + &
     1101600*log1mX*logX*pi**2*x**2 - &
     21600*logX**2*pi**2*x**2 - 20880*pi**4*x**2 + &
     31697960*x**3 + 97200*log1mX**4*x**3 + &
     10807530*logX*x**3 - 522360*log1mX*logX*x**3 - &
     388800*log1mX**3*logX*x**3 - 1719720*logX**2*x**3 - &
     92880*log1mX*logX**2*x**3 + &
     2786400*log1mX**2*logX**2*x**3 + 217440*logX**3*x**3 + &
     86400*log1mX*logX**3*x**3 - 43200*logX**4*x**3 - &
     2373720*pi**2*x**3 + 194400*log1mX**2*pi**2*x**3 - &
     98280*logX*pi**2*x**3 - &
     1101600*log1mX*logX*pi**2*x**3 + &
     21600*logX**2*pi**2*x**3 + 20880*pi**4*x**3 - &
     16191565*x**4 - 48600*log1mX**4*x**4 - &
     6204120*logX*x**4 - 457740*log1mX*logX*x**4 + &
     194400*log1mX**3*logX*x**4 + 958860*logX**2*x**4 + &
     122040*log1mX*logX**2*x**4 - &
     1393200*log1mX**2*logX**2*x**4 - 115920*logX**3*x**4 - &
     43200*log1mX*logX**3*x**4 + 21600*logX**4*x**4 + &
     1268640*pi**2*x**4 - 97200*log1mX**2*pi**2*x**4 + &
     54540*logX*pi**2*x**4 + 550800*log1mX*logX*pi**2*x**4 - &
     10800*logX**2*pi**2*x**4 - 10440*pi**4*x**4 + &
     3313295*x**5 + 9720*log1mX**4*x**5 + &
     1385085*logX*x**5 + 209430*log1mX*logX*x**5 - &
     38880*log1mX**3*logX*x**5 - 209430*logX**2*x**5 - &
     36720*log1mX*logX**2*x**5 + &
     278640*log1mX**2*logX**2*x**5 + 24480*logX**3*x**5 + &
     8640*log1mX*logX**3*x**5 - 4320*logX**4*x**5 - &
     268650*pi**2*x**5 + 19440*log1mX**2*pi**2*x**5 - &
     11880*logX*pi**2*x**5 - 110160*log1mX*logX*pi**2*x**5 + &
     2160*logX**2*pi**2*x**5 + 2088*pi**4*x**5 - &
     90*(-17824 + 27298*x - 26468*x**2 + 5804*x**3 + &
        5086*x**4 - 2327*x**5 + &
        576*logX**2*(-2 + 5*x - 10*x**2 + 10*x**3 - 5*x**4 + &
           x**5) + 792*pi**2*&
         (-2 + 5*x - 10*x**2 + 10*x**3 - 5*x**4 + x**5) - &
        24*logX*(-400 + 847*x - 1214*x**2 + 914*x**3 - &
           347*x**4 + 52*x**5 + &
           12*log1mX*&
            (-18 + 85*x - 170*x**2 + 170*x**3 - 85*x**4 + &
              17*x**5)))*Li2arg1 + &
     12960*(-18 + 85*x - 170*x**2 + 170*x**3 - 85*x**4 + &
        17*x**5)*Li2arg1**2 + &
     1080*(-496 + 937*x - 674*x**2 - 186*x**3 + 423*x**4 - &
        138*x**5 + 192*logX*&
         (-2 + 15*x - 30*x**2 + 30*x**3 - 15*x**4 + 3*x**5))*&
      Li3arg2 + 1451520*Li3arg1 - &
     414720*logX*Li3arg1 - 3298320*x*Li3arg1 + &
     1036800*logX*x*Li3arg1 + &
     5041440*x**2*Li3arg1 - &
     2073600*logX*x**2*Li3arg1 - &
     4134240*x**3*Li3arg1 + &
     2073600*logX*x**3*Li3arg1 + &
     1743120*x**4*Li3arg1 - &
     1036800*logX*x**4*Li3arg1 - &
     298080*x**5*Li3arg1 + &
     207360*logX*x**5*Li3arg1 + &
     1373760*Li4arg2 - 2656800*x*Li4arg2 + &
     5313600*x**2*Li4arg2 - &
     5313600*x**3*Li4arg2 + &
     2656800*x**4*Li4arg2 - &
     531360*x**5*Li4arg2 + 1088640*Li4arg1 - &
     129600*x*Li4arg1 + 259200*x**2*Li4arg1 - &
     259200*x**3*Li4arg1 + 129600*x**4*Li4arg1 - &
     25920*x**5*Li4arg1 + &
     570240*Li4arg3 + &
     1166400*x*Li4arg3 - &
     2332800*x**2*Li4arg3 + &
     2332800*x**3*Li4arg3 - &
     1166400*x**4*Li4arg3 + &
     233280*x**5*Li4arg3 + 1257120*Zeta3 + &
     466560*logX*Zeta3 - 6585840*x*Zeta3 - &
     3758400*logX*x*Zeta3 + 11875680*x**2*Zeta3 + &
     7516800*logX*x**2*Zeta3 - 11832480*x**3*Zeta3 - &
     7516800*logX*x**3*Zeta3 + 6002640*x**4*Zeta3 + &
     3758400*logX*x**4*Zeta3 - 1222560*x**5*Zeta3 - &
     751680*logX*x**5*Zeta3)/(6480.*(-1 + x)**5) 
!Else we must use a series around x=1
   else
    laurent(0) = (-4861)/300. + (34*(1 - x))/5. + &
     (1158*(1 - x)**2)/245. + (12221*(1 - x)**3)/3360. + &
     (135307*(1 - x)**4)/45360. + (2137*(1 - x)**5)/840. + &
     (2037547*(1 - x)**6)/914760. + &
     (2313539*(1 - x)**7)/1.16424e6 + &
     (11230871*(1 - x)**8)/6.24624e6 + &
     (4667798*(1 - x)**9)/2.837835e6 + &
     (1367573*(1 - x)**10)/900900. + &
     (12785009*(1 - x)**11)/9.06048e6 + &
     (1649107517*(1 - x)**12)/1.24972848e9 + (2*pi**2)/3.

    laurent(1) = ((1 - x)*(33269 - 2300*pi**2))/900. + &
     (1 - x)**2*(24.292933268545514d0 - (19*pi**2)/14.) + &
     (1 - x)**3*(19.045174201625095d0 - (17*pi**2)/18.) + &
     (1 - x)**4*(15.982499378866773d0 - (79*pi**2)/108.) + &
     (1 - x)**5*(13.920304201310154d0 - (3*pi**2)/5.) + &
     (1 - x)**6*(12.414757321241563d0 - (101*pi**2)/198.) + &
     (1 - x)**7*(11.25635033451521d0 - (4*pi**2)/9.) + &
     (1 - x)**8*(10.331444983820784d0 - (41*pi**2)/104.) + &
     (1 - x)**9*(9.572333309001118d0 - (67*pi**2)/189.) + &
     (1 - x)**10*(8.935831202852597d0 - (29*pi**2)/90.) + &
     (1 - x)**11*(8.392932669798085d0 - (13*pi**2)/44.) + &
     (1 - x)**12*(7.923344499259915d0 - (167*pi**2)/612.) + &
     (-112.30625925925926d0 + (907*pi**2)/180. + 32*Zeta3)

    laurent(2) = (1 - x)*(222.41511728395062d0 - (1154*pi**2)/135. - &
        (308*Zeta3)/3.) + &
     (1 - x)**2*(144.41899514476992d0 - (4397*pi**2)/735. - &
     (396*Zeta3)/7.) + &
     (1 - x)**3*(111.78366853395625d0 - (278809*pi**2)/60480. - &
     (121*Zeta3)/3.) + &
     (1 - x)**4*(93.01360223061745d0 - &
     (1026901*pi**2)/272160. - (286*Zeta3)/9.) + &
     (1 - x)**5*(80.52128610449985d0 - (242369*pi**2)/75600. - &
     (132*Zeta3)/5.) + &
     (1 - x)**6*(71.48236216870302d0 - &
     (15334901*pi**2)/5.48856d6 - (68*Zeta3)/3.) + &
     (1 - x)**7*(64.57750805992605d0 - &
     (17324597*pi**2)/6.98544d6 - (418*Zeta3)/21.) + &
     (1 - x)**8*(59.09762779754413d0 - &
     (250996139*pi**2)/1.1243232d8 - (231*Zeta3)/13.) + &
     (1 - x)**9*(54.62338374526518d0 - &
     (17296457*pi**2)/8.513505d6 - (1012*Zeta3)/63.) + &
     (1 - x)**10*(50.88897349792131d0 - &
     (916661*pi**2)/491400. - (44*Zeta3)/3.) + &
     (1 - x)**11*(47.7167928746431d0 - &
     (218855843*pi**2)/1.2684672d8 - (27*Zeta3)/2.) + &
     (1 - x)**12*(44.98314626600186d0 - &
     (4013303857d0*pi**2)/2.49945696d9 - (638*Zeta3)/51.)&
     + (-205476287 + 8842650*pi**2 + 354000*pi**4 + 58608000*Zeta3)/270000. 
   end if
end subroutine CalcC1qF

!*****C1g routine*****
subroutine CalcC1gF(x,laurent)
  use math  
  implicit none

!parameterlist
  real(kind(1d0)), intent(in) :: x
  real(kind(1d0)), dimension(-4:2), intent(out) :: laurent

!variables
  real(kind(1d0)) :: logX, log1mX, log2

  real(kind(1d0)), dimension(2:4) :: LiN
  real(kind(1d0)) :: Li2arg1, Li3arg1, Li4arg1
  real(kind(1d0)) :: Li2arg2, Li3arg2, Li4arg2
  real(kind(1d0)) :: Li2arg3, Li3arg3, Li4arg3
  real(kind(1d0)) :: Li2arg4, Li3arg4
  real(kind(1d0)) :: Li3arg5
  real(kind(1d0)) :: Li2arg6, Li3arg6
  real(kind(1d0)) :: Li3arg7

!Calculate logs and Li's once for all
  logX = log(x)
  log1mX = log(1-x)
  log2 = log(2d0)

!arg1: x
  call PolyLogs(x,LiN)

  Li2arg1 = LiN(2)
  Li3arg1 = LiN(3)
  Li4arg1 = LiN(4)

!arg2: 1 - x
  call PolyLogs(1d0 - x,LiN)

  Li2arg2 = LiN(2)
  Li3arg2 = LiN(3)
  Li4arg2 = LiN(4)

!arg3: x/(x - 1)
  call PolyLogs(x/(x - 1d0),LiN)

  Li2arg3 = LiN(2)
  Li3arg3 = LiN(3)
  Li4arg3 = LiN(4)

!arg4: 1 - x/2
  call PolyLogs(1d0 - x/2d0,LiN)

  Li2arg4 = LiN(2)
  Li3arg4 = LiN(3)

!arg5: x/(x - 2)
  call PolyLogs(-x/(x - 2d0),LiN)

  Li3arg5 = LiN(3)

!arg6: 1/(2 - x)
  call PolyLogs(1d0/(2d0 - x),LiN)

  Li2arg6 = LiN(2)
  Li3arg6 = LiN(3)

!arg5: x/2
  call PolyLogs(x/2d0,LiN)

  Li3arg7 = LiN(3)

  laurent = 0

!Calculate series
  laurent(-2) = 0
  
  laurent(-1) = -43d0/18d0 - 2*logX

!In case x<=0.9 we can use the original formula
  if (x.lt.0.9d0) then

   laurent(0) = (215552 - 18432*pi**2 + 432*logX**2*(-2 + x)**6*(-1 + x)**5 - &
    1735936*x - 3840*log2*x + 124416*pi**2*x + 6239744*x**2 + &
    19200*log2*x**2 - 414720*pi**2*x**2 - 13269104*x**3 - &
    38400*log2*x**3 + 858240*pi**2*x**3 + 18587424*x**4 + &
    38400*log2*x**4 - 1192320*pi**2*x**4 - 18022944*x**5 - &
    19200*log2*x**5 + 1153440*pi**2*x**5 + 12344056*x**6 + &
    3840*log2*x**6 - 789120*pi**2*x**6 - 5970557*x**7 + &
    381240*pi**2*x**7 + 1998262*x**8 - 127440*pi**2*x**8 - &
    440812*x**9 + 28080*pi**2*x**9 + 57714*x**10 - &
    3672*pi**2*x**10 - 3399*x**11 + 216*pi**2*x**11 - &
    12*logX*(-2944 + 34304*x - 138400*x**2 + 311056*x**3 - &
    453384*x**4 + 455100*x**5 - 321418*x**6 + 159635*x**7 - &
    54658*x**8 + 12302*x**9 - 1639*x**10 + 98*x**11 + &
    36*log1mX*(-2 + x)**7*(1 - 2*x + 4*x**2 - 3*x**3 + x**4)) - &
    432*(-2 + x)**7*(1 - 2*x + 4*x**2 - 3*x**3 + x**4)*&
    Li2arg1)/(216.*(-2 + x)**6*(-1 + x)**5)
 
   laurent(1) = (7080704 + 445440*logX - 1769472*log1mX*logX - 211968*logX**2 - &
    331776*log1mX*logX**2 + 110592*logX**3 - 397056*pi**2 + &
    27648*logX*pi**2 - 58633408*x - 82176*log2*x - &
    23040*log2**2*x - 10440960*logX*x + 7617024*log1mX*logX*x + &
    46080*log2*logX*x + 2469888*logX**2*x + &
    1824768*log1mX*logX**2*x - 884736*logX**3*x + &
    3558528*pi**2*x + 55296*logX*pi**2*x + 213848960*x**2 + &
    410880*log2*x**2 + 115200*log2**2*x**2 + 47373312*logX*x**2 - &
    14860800*log1mX*logX*x**2 - 230400*log2*logX*x**2 - &
    9964800*logX**2*x**2 - 5391360*log1mX*logX**2*x**2 + &
    3179520*logX**3*x**2 - 13325760*pi**2*x**2 - &
    587520*logX*pi**2*x**2 - 458671680*x**3 - 821760*log2*x**3 - &
    230400*log2**2*x**3 - 111635328*logX*x**3 + &
    15790464*log1mX*logX*x**3 + 460800*log2*logX*x**3 + &
    22396032*logX**2*x**3 + 10575360*log1mX*logX**2*x**3 - &
    6773760*logX**3*x**3 + 29036160*pi**2*x**3 + &
    1555200*logX*pi**2*x**3 + 646468912*x**4 + 821760*log2*x**4 + &
    230400*log2**2*x**4 + 167651136*logX*x**4 - &
    7174656*log1mX*logX*x**4 - 460800*log2*logX*x**4 - &
    32643648*logX**2*x**4 - 14411520*log1mX*logX**2*x**4 + &
    9504000*logX**3*x**4 - 41525520*pi**2*x**4 - &
    2324160*logX*pi**2*x**4 - 630149156*x**5 - 410880*log2*x**5 - &
    115200*log2**2*x**5 - 172305504*logX*x**5 - &
    3869856*log1mX*logX*x**5 + 230400*log2*logX*x**5 + &
    32767200*logX**2*x**5 + 13862016*log1mX*logX**2*x**5 - &
    9220608*logX**3*x**5 + 41061072*pi**2*x**5 + &
    2294784*logX*pi**2*x**5 + 433685608*x**6 + 82176*log2*x**6 + &
    23040*log2**2*x**6 + 124149840*logX*x**6 + &
    8561376*log1mX*logX*x**6 - 46080*log2*logX*x**6 - &
    23142096*logX**2*x**6 - 9471168*log1mX*logX**2*x**6 + &
    6312384*logX**3*x**6 - 28644156*pi**2*x**6 - &
    1577232*logX*pi**2*x**6 - 210674204*x**7 - &
    62732376*logX*x**7 - 6502680*log1mX*logX*x**7 + &
    11493720*logX**2*x**7 + 4574880*log1mX*logX**2*x**7 - &
    3049920*logX**3*x**7 + 14089440*pi**2*x**7 + &
    762480*logX*pi**2*x**7 + 70769674*x**8 + 21801024*logX*x**8 + &
    2842128*log1mX*logX*x**8 - 3935376*logX**2*x**8 - &
    1529280*log1mX*logX**2*x**8 + 1019520*logX**3*x**8 - &
    4787520*pi**2*x**8 - 254880*logX*pi**2*x**8 - 15660033*x**9 - &
    4969128*logX*x**9 - 752112*log1mX*logX*x**9 + &
    885744*logX**2*x**9 + 336960*log1mX*logX**2*x**9 - &
    224640*logX**3*x**9 + 1070400*pi**2*x**9 + &
    56160*logX*pi**2*x**9 + 2056010*x**10 + 669132*logX*x**10 + &
    112536*log1mX*logX*x**10 - 118008*logX**2*x**10 - &
    44064*log1mX*logX**2*x**10 + 29376*logX**3*x**10 - &
    141792*pi**2*x**10 - 7344*logX*pi**2*x**10 - 121387*x**11 - &
    40374*logX*x**11 - 7344*log1mX*logX*x**11 + &
    7056*logX**2*x**11 + 2592*log1mX*logX**2*x**11 - &
    1728*logX**3*x**11 + 8436*pi**2*x**11 + &
    432*logX*pi**2*x**11 - &
    288*(-512 + 1776*x - 3120*x**2 + 4240*x**3 - 5120*x**4 + &
       5020*x**5 - 3580*x**6 + 1765*x**7 - 590*x**8 + 130*x**9 - &
       17*x**10 + x**11)*Li2arg2 + &
    184320*(-1 + x)**5*x*Li2arg4 - &
    1769472*Li2arg1 + 1327104*logX*Li2arg1 + &
    7617024*x*Li2arg1 - 7299072*logX*x*Li2arg1 - &
    14860800*x**2*Li2arg1 + &
    21565440*logX*x**2*Li2arg1 + &
    15790464*x**3*Li2arg1 - &
    42301440*logX*x**3*Li2arg1 - 7174656*x**4*Li2arg1 + &
    57646080*logX*x**4*Li2arg1 - 3869856*x**5*Li2arg1 - &
    55448064*logX*x**5*Li2arg1 + 8561376*x**6*Li2arg1 + &
    37884672*logX*x**6*Li2arg1 - 6502680*x**7*Li2arg1 - &
    18299520*logX*x**7*Li2arg1 + 2842128*x**8*Li2arg1 + &
    6117120*logX*x**8*Li2arg1 - 752112*x**9*Li2arg1 - &
    1347840*logX*x**9*Li2arg1 + 112536*x**10*Li2arg1 + &
    176256*logX*x**10*Li2arg1 - 7344*x**11*Li2arg1 - &
    10368*logX*x**11*Li2arg1 - 2985984*Li3arg2 + &
    13105152*x*Li3arg2 - 31933440*x**2*Li3arg2 + &
    56194560*x**3*Li3arg2 - 73301760*x**4*Li3arg2 + &
    69558912*x**5*Li3arg2 - 47376576*x**6*Li3arg2 + &
    22874400*x**7*Li3arg2 - 7646400*x**8*Li3arg2 + &
    1684800*x**9*Li3arg2 - 220320*x**10*Li3arg2 + &
    12960*x**11*Li3arg2 - 3317760*Li3arg1 + &
    18247680*x*Li3arg1 - 53913600*x**2*Li3arg1 + &
    105753600*x**3*Li3arg1 - 144115200*x**4*Li3arg1 + &
    138620160*x**5*Li3arg1 - 94711680*x**6*Li3arg1 + &
    45748800*x**7*Li3arg1 - 15292800*x**8*Li3arg1 + &
    3369600*x**9*Li3arg1 - 440640*x**10*Li3arg1 + &
    25920*x**11*Li3arg1 + 663552*Zeta3 + &
    2985984*x*Zeta3 - 22394880*x**2*Zeta3 + &
    56816640*x**3*Zeta3 - 83980800*x**4*Zeta3 + &
    82674432*x**5*Zeta3 - 56785536*x**6*Zeta3 + &
    27449280*x**7*Zeta3 - 9175680*x**8*Zeta3 + &
    2021760*x**9*Zeta3 - 264384*x**10*Zeta3 + &
    15552*x**11*Zeta3)/(1296.*(-2 + x)**6*(-1 + x)**5)

    laurent(2) = -545.3567386831276d0 - (3*log1mX**4)/2. - (269801*logX)/1296. - &
      (2351*log1mX*logX)/72. + 6*log1mX**3*logX + (2243*logX**2)/72. + &
      (37*log1mX*logX**2)/9. - 43*log1mX**2*logX**2 - (98*logX**3)/27. - &
      (4*log1mX*logX**3)/3. + (2*logX**4)/3. + (343*pi**2)/8. - &
      3*log1mX**2*pi**2 + (113*logX*pi**2)/54. + 17*log1mX*logX*pi**2 - &
      (logX**2*pi**2)/3. - (29*pi**4)/90. + (142*logX)/(6 - 3*x) + &
      (19232*log2)/(81.*(-2 + x)**6) + (3424*log2**2)/(27.*(-2 + x)**6) - &
      (4480*log2**3)/(27.*(-2 + x)**6) - (19232*logX)/(81.*(-2 + x)**6) - &
      (6848*log2*logX)/(27.*(-2 + x)**6) + (3200*log2**2*logX)/(9.*(-2 + x)**6) + &
      (3424*logX**2)/(27.*(-2 + x)**6) - (2240*log1mX*logX**2)/(9.*(-2 + x)**6) - &
      (640*log2*logX**2)/(3.*(-2 + x)**6) - (640*logX**3)/(27.*(-2 + x)**6) - &
      (856*pi**2)/(81.*(-2 + x)**6) + (640*log2*pi**2)/(27.*(-2 + x)**6) + &
      (320*logX*pi**2)/(27.*(-2 + x)**6) + 142864/(81.*(-2 + x)**5) + &
      (9616*log2)/(81.*(-2 + x)**5) + (1712*log2**2)/(27.*(-2 + x)**5) - &
      (2240*log2**3)/(27.*(-2 + x)**5) - (42928*logX)/(81.*(-2 + x)**5) - &
      (3424*log2*logX)/(27.*(-2 + x)**5) + (1600*log2**2*logX)/(9.*(-2 + x)**5) + &
      (2672*logX**2)/(27.*(-2 + x)**5) - (1120*log1mX*logX**2)/(9.*(-2 + x)**5) - &
      (320*log2*logX**2)/(3.*(-2 + x)**5) - (320*logX**3)/(27.*(-2 + x)**5) - &
      (2108*pi**2)/(81.*(-2 + x)**5) + (320*log2*pi**2)/(27.*(-2 + x)**5) + &
      (160*logX*pi**2)/(27.*(-2 + x)**5) + 58132/(81.*(-2 + x)**4) - &
      (140*pi**2)/(27.*(-2 + x)**4) - 2014/(27.*(-2 + x)**3) - &
      (4136*logX)/(81.*(-2 + x)**3) + (80*logX**2)/(27.*(-2 + x)**3) + &
      (70*pi**2)/(81.*(-2 + x)**3) + 1874/(81.*(-2 + x)**2) + &
      (4136*logX)/(81.*(-2 + x)**2) - (80*logX**2)/(27.*(-2 + x)**2) - &
      (35*pi**2)/(162.*(-2 + x)**2) - 1975/(162.*(-2 + x)) + &
      (8*logX**2)/(3.*(-2 + x)) + (7*pi**2)/(108.*(-2 + x)) - &
      (31*log1mX**4)/(6.*(-1 + x)**5) - (277879*logX)/(1296.*(-1 + x)**5) - &
      (8263*log1mX*logX)/(72.*(-1 + x)**5) + (62*log1mX**3*logX)/(3.*(-1 + x)**5) + &
      (1877*logX**2)/(72.*(-1 + x)**5) + (95*log1mX*logX**2)/(9.*(-1 + x)**5) - &
      (29*log1mX**2*logX**2)/(-1 + x)**5 - (52*logX**3)/(27.*(-1 + x)**5) + &
      (4*log1mX*logX**3)/(3.*(-1 + x)**5) + (8263*pi**2)/(432.*(-1 + x)**5) - &
      (31*log1mX**2*pi**2)/(3.*(-1 + x)**5) + (68*logX*pi**2)/(27.*(-1 + x)**5) + &
      (29*log1mX*logX*pi**2)/(3.*(-1 + x)**5) - &
      (2*logX**2*pi**2)/(3.*(-1 + x)**5) + (107*pi**4)/(90.*(-1 + x)**5) + &
      154093/(1296.*(-1 + x)**4) - (20659*logX)/(144.*(-1 + x)**4) + &
      (161*log1mX*logX)/(24.*(-1 + x)**4) + (377*logX**2)/(24.*(-1 + x)**4) - &
      (log1mX*logX**2)/(2.*(-1 + x)**4) - logX**3/(-1 + x)**4 + &
      (4925*pi**2)/(432.*(-1 + x)**4) + (3*logX*pi**2)/(4.*(-1 + x)**4) + &
      39779/(648.*(-1 + x)**3) + (2395*logX)/(324.*(-1 + x)**3) - &
      (905*log1mX*logX)/(36.*(-1 + x)**3) - (157*logX**2)/(108.*(-1 + x)**3) + &
      (7*log1mX*logX**2)/(3.*(-1 + x)**3) + (2*logX**3)/(9.*(-1 + x)**3) - &
      (137*pi**2)/(72.*(-1 + x)**3) - (logX*pi**2)/(6.*(-1 + x)**3) - &
      6923/(648.*(-1 + x)**2) + (1951*logX)/(108.*(-1 + x)**2) + &
      (1439*log1mX*logX)/(36.*(-1 + x)**2) - (97*logX**2)/(36.*(-1 + x)**2) - &
      (13*log1mX*logX**2)/(3.*(-1 + x)**2) + (2*logX**3)/(9.*(-1 + x)**2) - &
      (745*pi**2)/(324.*(-1 + x)**2) - (logX*pi**2)/(6.*(-1 + x)**2) - &
      10459/(216.*(-1 + x)) - (8717*logX)/(144.*(-1 + x)) - &
      (2183*log1mX*logX)/(24.*(-1 + x)) + (85*logX**2)/(8.*(-1 + x)) + &
      (19*log1mX*logX**2)/(2.*(-1 + x)) - logX**3/(-1 + x) + &
      (2501*pi**2)/(216.*(-1 + x)) + (3*logX*pi**2)/(4.*(-1 + x)) + &
      (1280*log2**2*Log(2 - x))/(9.*(-2 + x)**6) - &
      (2560*log2*logX*Log(2 - x))/(9.*(-2 + x)**6) + &
      (2560*logX**2*Log(2 - x))/(9.*(-2 + x)**6) + &
      (640*log2**2*Log(2 - x))/(9.*(-2 + x)**5) - &
      (1280*log2*logX*Log(2 - x))/(9.*(-2 + x)**5) + &
      (1280*logX**2*Log(2 - x))/(9.*(-2 + x)**5) + &
      ((72*logX*(-1 - 160/(-2 + x)**6 - 80/(-2 + x)**5 - (-1 + x)**(-5)) + &
      (31744 - 95232*x + 130560*x**2 - 132416*x**3 + 149824*x**4 - &
      175106*x**5 + 156246*x**6 - 93900*x**7 + 36678*x**8 - &
      9037*x**9 + 1287*x**10 - 81*x**11)/((-2 + x)**6*(-1 + x)**5))*&
      Li2arg2)/54. + &
      (1280*logX*x*Li2arg6)/(9.*(-2 + x)**6) + &
      (27392*Li2arg4)/(27.*(-2 + x)**6) + &
      (13696*Li2arg4)/(27.*(-2 + x)**5) - &
      (2351*Li2arg1)/72. - (52*logX*Li2arg1)/3. - &
      68*log1mX*logX*Li2arg1 + 8*logX**2*Li2arg1 + &
      11*pi**2*Li2arg1 - (8263*Li2arg1)/(72.*(-1 + x)**5) + &
      (148*logX*Li2arg1)/(3.*(-1 + x)**5) + (4*log1mX*logX*Li2arg1)/(-1 + x)**5 - &
      (8*logX**2*Li2arg1)/(-1 + x)**5 - (11*pi**2*Li2arg1)/(-1 + x)**5 + &
      (161*Li2arg1)/(24.*(-1 + x)**4) - (11*logX*Li2arg1)/(-1 + x)**4 - &
      (905*Li2arg1)/(36.*(-1 + x)**3) + (34*logX*Li2arg1)/(3.*(-1 + x)**3) + &
      (1439*Li2arg1)/(36.*(-1 + x)**2) - (46*logX*Li2arg1)/(3.*(-1 + x)**2) - &
      (2183*Li2arg1)/(24.*(-1 + x)) + (29*logX*Li2arg1)/(-1 + x) - &
      34*Li2arg1**2 + (2*Li2arg1**2)/(-1 + x)**5 + &
      (193*Li3arg2)/9. - 96*logX*Li3arg2 - &
      (2240*Li3arg2)/(9.*(-2 + x)**6) - (1120*Li3arg2)/(9.*(-2 + x)**5) + &
      (299*Li3arg2)/(9.*(-1 + x)**5) - (32*logX*Li3arg2)/(-1 + x)**5 - &
      (11*Li3arg2)/(2.*(-1 + x)**4) + (37*Li3arg2)/(3.*(-1 + x)**3) - &
      (21*Li3arg2)/(-1 + x)**2 + (89*Li3arg2)/(2.*(-1 + x)) + &
      (2560*Li3arg6)/(9.*(-2 + x)**6) + (1280*Li3arg6)/(9.*(-2 + x)**5) + &
      (17920*Li3arg4)/(9.*(-2 + x)**6) + (8960*Li3arg4)/(9.*(-2 + x)**5) + &
      (2560*Li3arg7)/(9.*(-2 + x)**6) + &
      (1280*Li3arg7)/(9.*(-2 + x)**5) + (386*Li3arg1)/9. - &
      32*logX*Li3arg1 - (640*Li3arg1)/(3.*(-2 + x)**6) - &
      (320*Li3arg1)/(3.*(-2 + x)**5) - (698*Li3arg1)/(9.*(-1 + x)**5) + &
      (32*logX*Li3arg1)/(-1 + x)**5 + (21*Li3arg1)/(-1 + x)**4 - &
      (18*Li3arg1)/(-1 + x)**3 + (22*Li3arg1)/(-1 + x)**2 - &
      (39*Li3arg1)/(-1 + x) - &
      (2560*Li3arg5)/(9.*(-2 + x)**6) - &
      (1280*Li3arg5)/(9.*(-2 + x)**5) + &
      82*Li4arg2 - (130*Li4arg2)/(-1 + x)**5 + &
      4*Li4arg1 - (164*Li4arg1)/(-1 + x)**5 - &
      36*Li4arg3 - (124*Li4arg3)/(-1 + x)**5 + &
      (1814*Zeta3)/9. + 116*logX*Zeta3 - &
      (2560*Zeta3)/(9.*(-2 + x)**6) - (1280*Zeta3)/(9.*(-2 + x)**5) + &
      (698*Zeta3)/(9.*(-1 + x)**5) + (44*logX*Zeta3)/(-1 + x)**5 + &
      (67*Zeta3)/(-1 + x)**4 - (26*Zeta3)/(-1 + x)**3 + &
      (22*Zeta3)/(3.*(-1 + x)**2) + (17*Zeta3)/(-1 + x)
!Else we must use a series around x=1
    else
     laurent(0) = (1 - x)**11*(90831.38203373332d0 - 131040*log2) + &
      (1 - x)**9*(40530.68805462615d0 - (526240*log2)/9.) + &
      (1 - x)**7*(15454.495103243318d0 - (66880*log2)/3.) + &
      (1 - x)**5*(4660.421825396826d0 - 6720*log2) + &
      (1 - x)**3*(952.3686838624338d0 - (12320*log2)/9.) + &
      (1 - x)*(92.78148148148148d0 - (1120*log2)/9.) + &
      (-29.766296296296296d0 + (2*pi**2)/3. + &
         (160*log2)/9.) + &
      (1 - x)**2*(-328.1385487528345d0 + 480*log2) + &
      (1 - x)**4*(-2239.8195105820105d0 + (29120*log2)/9.) + &
      (1 - x)**6*(-8796.181679347588d0 + (38080*log2)/3.) + &
      (1 - x)**8*(-25616.96907296763d0 + 36960*log2) + &
      (1 - x)**10*(-61673.21656035323d0 + (800800*log2)/9.) + &
      (1 - x)**12*(-130076.2535441429d0 + (1688960*log2)/9.)

     laurent(1) = (1 - x)**12*(42358.44112314064d0 + &
      (200985739d0*pi**2)/1836. - (1095561232d0*log2)/891. - &
      (1688960d0*log2**2)/3.) + &
      (1 - x)**10*(-43864.680117656244d0 + &
      (14013913d0*pi**2)/270. - (39748984d0*log2)/81. - &
      (800800*log2**2)/3.) + &
      (1 - x)**8*(-48705.50578971287d0 + (2242199*pi**2)/104. - &
      (1438664*log2)/9. - 110880*log2**2) + &
      (1 - x)**6*(-28997.850462760773d0 + (1465979*pi**2)/198. - &
      (1003952*log2)/27. - 38080*log2**2) + &
      (1 - x)**4*(-11155.573786139455d0 + (203761*pi**2)/108. - &
      (36272*log2)/9. - (29120*log2**2)/3.) + &
      (1 - x)**2*(-2335.586034985423d0 + (3901*pi**2)/14. + &
      (3728*log2)/9. - 1440*log2**2) + &
      (1 - x)*(756.8121604938272d0 - (2029*pi**2)/27. - &
      (7184*log2)/27. + (1120*log2**2)/3.) + &
      (1 - x)**3*(5683.416261062295d0 - (43171*pi**2)/54. + &
      (3152*log2)/9. + (12320*log2**2)/3.) + &
      (1 - x)**5*(19092.80789238473d0 - (19603*pi**2)/5. + &
      (128848*log2)/9. + 20160*log2**2) + &
      (1 - x)**7*(39715.32381714167d0 - (117044*pi**2)/9. + &
      (15408944d0*log2)/189. + 66880*log2**2) + &
      (1 - x)**9*(52071.36842738248d0 - (6446507*pi**2)/189. + &
      (163841992d0*log2)/567. + (526240*log2**2)/3.) + &
      (1 - x)**11*(16013.680067713136d0 - (3363373*pi**2)/44. + &
      (78490648*log2)/99. + 393120*log2**2) + &
      (-238.78391358024692d0 + (8461*pi**2)/540. + &
      (1712*log2)/27. - (160*log2**2)/3. + 32*Zeta3)

     laurent(2) = -1655.9840279835391d0 + (59*pi**4)/45. + &
      pi**2*(85.56166666666667d0 - (2240*log2)/27.) + &
      (9616*log2)/81. - (1712*log2**2)/9. + &
      (320*log2**3)/3. + &
      (1 - x)**11*(1.6278122525122717d6 - &
      (58919199568d0*log2)/31185. - &
      (78490648*log2**2)/33. - 786240*log2**3 + &
      pi**2*(356516.95835538354d0 + 611520*log2) - &
      (11007387*Zeta3)/2.) + &
      (1 - x)**9*(769173.3446367856d0 - (852653800d0*log2)/1701. - &
      (163841992d0*log2**2)/189. - &
      (1052480*log2**3)/3. + &
      pi**2*(121275.80477039127d0 + (7367360*log2)/27.) - &
      (154715572*Zeta3)/63.) + &
      (1 - x)**7*(336761.99443528714d0 - (42784552*log2)/567. - &
      (15408944d0*log2**2)/63. - 133760*log2**3 + &
      pi**2*(29528.06310883781d0 + (936320*log2)/9.) - &
      (19663138d0*Zeta3)/21.) + &
      (1 - x)**5*(127534.10097568661d0 + (99256*log2)/27. - &
      (128848*log2**2)/3. - 40320*log2**3 + &
      pi**2*(2913.8625793650795d0 + 31360*log2) - &
      (1411332*Zeta3)/5.) + &
      (1 - x)**3*(35747.2860591314d0 + (24896*log2)/9. - &
      (3152*log2**2)/3. - (24640*log2**3)/3. + &
      pi**2*(-907.1670359347443d0 + (172480*log2)/27.) - &
      (172601*Zeta3)/3.) + &
      (1 - x)*(5050.9001327160495d0 - (15952*log2)/81. + &
      (7184*log2**2)/9. - (2240*log2**3)/3. + &
      (pi**2*(-427793 + 940800*log2))/1620. - &
      (15988*Zeta3)/3.) + (43808*Zeta3)/45. + &
      (1 - x)**2*(-14978.130310410786d0 + &
      pi**2*(623.9899092970521d0 - 2240*log2) - &
      (14432*log2)/27. - (3728*log2**2)/3. + &
      2880*log2**3 + (140724*Zeta3)/7.) + &
      (1 - x)**4*(-71034.39334032223d0 + &
      pi**2*(261.8667658730159d0 - (407680*log2)/27.) - &
      (149096*log2)/27. + (36272*log2**2)/3. + &
      (58240*log2**3)/3. + (1222754*Zeta3)/9.) + &
      (1 - x)**6*(-212346.74149328744d0 + &
      pi**2*(-11428.304356880493d0 - (533120*log2)/9.) + &
      (1240840*log2)/81. + (1003952*log2**2)/9. + &
      76160*log2**3 + (1599292*Zeta3)/3.) + &
      (1 - x)**8*(-515239.0574650315d0 + &
      pi**2*(-63360.411550869d0 - 172480*log2) + &
      (13663928*log2)/63. + (1438664*log2**2)/3. + &
      221760*log2**3 + (20179929*Zeta3)/13.) + &
      (1 - x)**10*(-1.126928974061584d6 + &
      pi**2*(-214306.12363745307d0 - &
      (11211200*log2)/27.) + &
      (8655660032d0*log2)/8505. + &
      (39748984*log2**2)/27. + (1601600*log2**3)/3. + &
      3737052*Zeta3) + &
      (1 - x)**12*(-2.3229349522564285d6 + &
      pi**2*(-565503.006431747d0 - (23645440d0*log2)/27.) + &
      (102296782748d0*log2)/31185. + &
      (1095561232d0*log2**2)/297. + &
      (3377920*log2**3)/3. + (133990614d0*Zeta3)/17.)
    end if
end subroutine CalcC1gF

!*****S1ik routine*****
subroutine CalcS1ikFF(Y,laurent)
  use math  
  implicit none
!
!parameterlist
  real(kind(1d0)), intent(in) :: Y
  real(kind(1d0)), dimension(-4:2), intent(out) :: laurent
!
!variables
  real(kind(1d0)) :: logY, log1mY
!
  real(kind(1d0)), dimension(2:4) :: LiN
  real(kind(1d0)) :: Li3arg1
  real(kind(1d0)) :: Li2arg2, Li3arg2
  real(kind(1d0)) :: Li4arg3
!
!Calculate logs and Li's once for all
  logY = log(Y)
  log1mY = log(1-Y)

!arg1: Y
  call PolyLogs(Y,LiN)

  Li3arg1 = LiN(3)

!arg2: 1 - Y
  call PolyLogs(1d0 - Y,LiN)

  Li2arg2 = LiN(2)
  Li3arg2 = LiN(3)

!arg3: (Y - 1)/Y
  call PolyLogs((Y - 1d0)/Y,LiN)

  Li4arg3 = LiN(4)

  laurent = 0

!Calculate series
  laurent(-2) = -1

  laurent(-1) = -11d0/3d0 + logY

  laurent(0) = (-317 + 66*logY - 9*logY**2 + 21*pi**2 - 18*Li2arg2)/18.

!In case x<=0.9 we can use the original formula
  if (Y.lt.0.9d0) then

   laurent(1) = - 9299d0/108d0 + (317*logY)/18. - (11*logY**2)/6. - (log1mY*logY**2)/2. + &
    logY**3/6. + (77*pi**2)/18. - logY*pi**2 - (11*Li2arg2)/3. - &
    Li3arg2 - Li3arg1 + 33*Zeta3

   laurent(2) = -275645d0/648d0 + (9299*logY)/108. - (317*logY**2)/36. - &
    (11*log1mY*logY**2)/6. + (11*logY**3)/18. + (2219*pi**2)/108. - &
    (11*logY*pi**2)/3. + (7*logY**2*pi**2)/12. + (33*pi**4)/40. + &
    ((-317 + 21*pi**2)*Li2arg2)/18. - (11*Li3arg2)/3. - &
    (11*Li3arg1)/3. + Li4arg3 + (121 - 32*logY)*Zeta3

  else
   laurent(1) =  -86.10185185185185d0 + (77*pi**2)/18. + &
    (1 - Y)**12*(-2.882249736703425d0 + (7*pi**2)/72.) + &
    (1 - Y)**11*(-3.0929854904576963d0 + (7*pi**2)/66.) + &
    (1 - Y)**10*(-3.3414972757621566d0 + (7*pi**2)/60.) + &
    (1 - Y)**9*(-3.6394899183953417d0 + (7*pi**2)/54.) + &
    (1 - Y)**8*(-4.004209006519274d0 + (7*pi**2)/48.) + &
    (1 - Y)**7*(-4.46222951085196d0 + pi**2/6.) + &
    (1 - Y)**6*(-5.056898148148148d0 + (7*pi**2)/36.) + &
    (1 - Y)**5*(-5.864388888888889d0 + (7*pi**2)/30.) + &
    (1 - Y)**4*(-7.032986111111111d0 + (7*pi**2)/24.) + &
    (1 - Y)**3*(-8.898148148148149d0 + (7*pi**2)/18.) + &
    (1 - Y)**2*(-12.430555555555555d0 + (7*pi**2)/12.) + &
    (1 - Y)*(-22.27777777777778d0 + (7*pi**2)/6.) + 32*Zeta3

   laurent(2) =   -425.37808641975306d0 + (2219*pi**2)/108. + &
    (33*pi**4)/40. + (352*Zeta3)/3. + &
    (1 - Y)**12*(-14.090438951680945d0 + (187661*pi**2)/285120. + &
    (8*Zeta3)/3.) + &
    (1 - Y)**11*(-15.109365243951629d0 + (185351*pi**2)/261360. + &
    (32*Zeta3)/11.) + &
    (1 - Y)**10*(-16.31100321930102d0 + (16621*pi**2)/21600. + &
    (16*Zeta3)/5.) + &
    (1 - Y)**9*(-17.752059435249784d0 + (16369*pi**2)/19440. + &
    (32*Zeta3)/9.) + &
    (1 - Y)**8*(-19.51613062485659d0 + (5363*pi**2)/5760. + &
    4*Zeta3) + (1 - Y)**7*&
    (-21.732119427476206d0 + (2629*pi**2)/2520. + &
    (32*Zeta3)/7.) + &
    (1 - Y)**6*(-24.610434413580247d0 + (2569*pi**2)/2160. + &
    (16*Zeta3)/3.) + &
    (1 - Y)**5*(-28.521141666666665d0 + (833*pi**2)/600. + &
    (32*Zeta3)/5.) + &
    (1 - Y)**4*(-34.185329861111114d0 + (161*pi**2)/96. + &
    8*Zeta3) + (1 - Y)**2*&
    (-60.405092592592595d0 + (217*pi**2)/72. + 16*Zeta3) + &
    (1 - Y)*(-108.37962962962963d0 + (49*pi**2)/9. + 32*Zeta3) + &
    ((1 - Y)**3*(-3113 + 154*pi**2 + 768*Zeta3))/72.
  end if
 end subroutine CalcS1ikFF
!
subroutine CalcInnlo(p,                  &
                     smeB,BornLaurent,   &
                     Bij,BijLaurent,     &
                     Bijkl,BijklLaurent, &
                     smeV,VijLaurent,    &
                     Inloterm,Innloterm,InloLaurent,InnloLaurent)
use IopInvariants
use particles
use QCDparams
use scales
use math
use misc
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
  real(kind(1d0)) , intent(in) :: smeB,smeV
  real(kind(1d0)) , dimension(-4:2) :: BornLaurent
  real(kind(1d0)) , dimension(:,:) , intent(in) :: Bij
  real(kind(1d0)) , dimension(:,:,-4:) :: BijLaurent
  real(kind(1d0)) , dimension(:,:,:,:) , intent(in) :: Bijkl
  real(kind(1d0)) , dimension(:,:,:,:,-4:) :: BijklLaurent
  real(kind(1d0)) , dimension(:,:,-4:) , intent(in) :: VijLaurent
  real(kind(1d0)) , intent(out) :: Inloterm,Innloterm
  real(kind(1d0)) , dimension(-4:2) , intent(out) :: InloLaurent
  real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: &
    InnloLaurent
!
  integer :: ipart,kpart,jpart,lpart,i,npart
  real(kind(1d0)) :: Innlofit
  real(kind(1d0)) :: L
  real(kind(1d0)) , dimension(-4:2) :: InnlofitLaurent
  real(kind(1d0)) , dimension(-4:2) :: smeVLaurent,I1VLaurent
  real(kind(1d0)) , dimension(-4:2) :: Laurent1,Laurent2
  real(kind(1d0)) , dimension(-4:2) :: CLaurent,SLaurent,I1I1Laurent
!
! DEBUG!!!!
  real(kind(1d0)) :: y12,y13,y23
!
  interface
    subroutine InnloUser(p,yiQ,YikQ,         &
                         smeB,BornLaurent,   &
                         Bij,BijLaurent,     &
                         Bijkl,bijklLaurent, &
                         Ifit,IfitLaurent)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:) , intent(in) :: yiQ
      real(kind(1d0)) , dimension(:,:) , intent(in) :: YikQ
      real(kind(1d0)) , intent(in) :: smeB
      real(kind(1d0)) , dimension(-4:2) :: BornLaurent
      real(kind(1d0)) , dimension(:,:) , intent(in) :: Bij
      real(kind(1d0)) , dimension(:,:,-4:) :: BijLaurent
      real(kind(1d0)) , dimension(:,:,:,:) , intent(in) :: Bijkl
      real(kind(1d0)) , dimension(:,:,:,:,-4:) :: BijklLaurent
      real(kind(1d0)) , intent(out) :: Ifit
      real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: IfitLaurent
!
    end subroutine InnloUser
!
    subroutine CalcC1qF(x,laurent)
    implicit none
!
      real(kind(1d0)) , intent(in) :: x
      real(kind(1d0)) , dimension(-4:2) , intent(out) :: laurent
!
    end subroutine CalcC1qF
!
    subroutine CalcC1gF(x,laurent)
    implicit none
!
      real(kind(1d0)) , intent(in) :: x
      real(kind(1d0)) , dimension(-4:2) , intent(out) :: laurent
!
    end subroutine CalcC1gF
!
    subroutine CalcS1ikFF(Y,laurent)
    implicit none
!
      real(kind(1d0)) , intent(in) :: Y
      real(kind(1d0)) , dimension(-4:2) , intent(out) :: laurent
!
    end subroutine CalcS1ikFF
  end interface
!
  npart = size(p)
!
  Inloterm    = 0
  InloLaurent = 0
  Innloterm   = 0
  I1VLaurent  = 0
  if (present(InnloLaurent)) InnloLaurent = 0
! We calculate the invariants:
  call CalcIopInvariants(p)
!
! log factor for scale dependence:
  L = log(mur**2/Q2)
!
! DEBUG!!!!
  y12 = YikQ(3,4)*yiQ(3)*yiQ(4)
  y13 = YikQ(3,5)*yiQ(3)*yiQ(5)
  y23 = YikQ(4,5)*yiQ(4)*yiQ(5)
!  print *,"y12,y13,y23: ",y12,y13,y23
!  print *,"YikQ: ",YikQ(3,5),YikQ(4,5)
!  print *,"yiQ: ",yiQ(3:5)
!
!***********************************************************************
! Calculation of I_1^{(0)} \times \ud\sigma^B with its full Laurent 
! series expansion up to O(ep^2) terms needed not only for pole 
! cancellation checks but for renormalization of the I_1^1.
!***********************************************************************
! We loop over all particles:
  do ipart=1,npart
! But only considering massless partons:
    if (abs(p(ipart)%flv).gt.qcd_nf) cycle
! We differentiate between two scenarios: Initial and Final
    if (ipart.le.2) then
! Initial State:
      Print *,"Initial state is not implemented yet..."
      stop "CalcInnlo"
    else
! Final State:
! (anti)Quark:
      if (p(ipart)%flv.ne.0) then
        call CalcC1qF(yiQ(ipart),CLaurent)
!        InloLaurent = InloLaurent + CLaurent*qcd_cf*smeB
        InloLaurent = InloLaurent + qcd_cf*SeriesProd(CLaurent,BornLaurent)
! Gluon:
      else
        call CalcC1gF(yiQ(ipart),CLaurent)
!        InloLaurent = InloLaurent + CLaurent*qcd_ca*smeB
        InloLaurent = InloLaurent + qcd_ca*SeriesProd(CLaurent,BornLaurent)
      end if
! The Ti2 factor is obtained through the color-correlated SME:
      I1VLaurent = I1VLaurent &
                 + SeriesProd(CLaurent,VijLaurent(ipart,ipart,-4:))
    end if
  end do
! We consider the contributions coming from the soft integrated
! counterterms:
  do ipart=1,npart
    if (abs(p(ipart)%flv).gt.qcd_nf) cycle
    do kpart=ipart+1,npart
      if (abs(p(kpart)%flv).gt.qcd_nf) cycle
      call CalcS1ikFF(YikQ(ipart,kpart),SLaurent)
!      InloLaurent = InloLaurent + 2*SLaurent &
!                  * Bij(ipart,kpart)
      InloLaurent = InloLaurent + 2*SeriesProd(SLaurent, &
                    BijLaurent(ipart,kpart,:))
      I1VLaurent = I1VLaurent &
                 + 2*SeriesProd(SLaurent,VijLaurent(ipart,kpart,-4:))
    end do
  end do
!***********************************************************************
! Calculation of I_1^{(0)} * I_1^{(0)} \times \ud\sigma^B
!***********************************************************************
  I1I1Laurent = 0
! Calculation of the C_{1,i}^{(0)}C_{1,k}^{(0)}Ti2*Tk2 contribution:
  do ipart=1,npart
    if (abs(p(ipart)%flv).gt.qcd_nf) cycle
! (anti)Quark:
    if (p(ipart)%flv.ne.0) then
      call CalcC1qF(yiQ(ipart),Laurent1)
      Laurent1 = Laurent1 * qcd_cf
! Gluon:
    else
      call CalcC1gF(yiQ(ipart),Laurent1)
      Laurent1 = Laurent1 * qcd_ca
    end if
    do kpart=ipart,npart
      if (abs(p(kpart)%flv).gt.qcd_nf) cycle
! (anti)Quark:
      if (p(kpart)%flv.ne.0) then
        call CalcC1qF(yiQ(kpart),Laurent2)
        Laurent2 = Laurent2 * qcd_cf
! Gluon:
      else
        call CalcC1gF(yiQ(kpart),Laurent2)
        Laurent2 = Laurent2 * qcd_ca
      end if
! Note the factor of 2 due to the ordered summation:
      if (ipart.ne.kpart) Laurent2  = 2*Laurent2
!      I1I1Laurent = I1I1Laurent + SeriesProd(Laurent1,Laurent2)*smeB
      I1I1Laurent = I1I1Laurent + SeriesProd(SeriesProd(Laurent1,Laurent2),BornLaurent)
    end do
  end do
! Calculation of the S_1^{(0)(i,k)}S_1^{(0)(j,l)} (TiTk)(TjTl):
  do ipart=1,npart-1
    if (abs(p(ipart)%flv).gt.qcd_nf) cycle
    do kpart=ipart+1,npart
      if (abs(p(kpart)%flv).gt.qcd_nf) cycle
      call CalcS1ikFF(YikQ(ipart,kpart),Laurent1)
      do jpart=1,npart-1
        if (abs(p(jpart)%flv).gt.qcd_nf) cycle
        do lpart=jpart+1,npart
          if (abs(p(lpart)%flv).gt.qcd_nf) cycle
          call CalcS1ikFF(YikQ(jpart,lpart),Laurent2)
! Note the factor of due to ordered sums:
!          I1I1Laurent = I1I1Laurent + 2*SeriesProd(Laurent1,Laurent2) &
!                      * Bijkl(ipart,kpart,jpart,lpart)
          I1I1Laurent = I1I1Laurent + 2*SeriesProd(SeriesProd(Laurent1,Laurent2), &
                        BijklLaurent(ipart,kpart,jpart,lpart,:))
        end do
      end do
    end do
  end do
! Calculation of the C_{1,i}^{(0)}S_1^{(0)(j,l)} Ti2*(TjTl):
  do ipart=1,npart
    if (abs(p(ipart)%flv).gt.qcd_nf) cycle
! (anti)Quark:
    if (p(ipart)%flv.ne.0) then
      call CalcC1qF(yiQ(ipart),Laurent1)
      Laurent1 = Laurent1 * qcd_cf
! Gluon:
    else
      call CalcC1gF(yiQ(ipart),Laurent1)
      Laurent1 = Laurent1 * qcd_ca
    end if
    do jpart=1,npart-1
      if (abs(p(jpart)%flv).gt.qcd_nf) cycle
      do lpart=jpart+1,npart
        if (abs(p(lpart)%flv).gt.qcd_nf) cycle
        call CalcS1ikFF(YikQ(jpart,lpart),Laurent2)
! Factor of 4 due to the ordered soft sum and because of squaring
! the I_1^{(0)}:
!        I1I1Laurent = I1I1Laurent + 4*SeriesProd(Laurent1,Laurent2) &
!                    * Bij(jpart,lpart)
        I1I1Laurent = I1I1Laurent + 4*SeriesProd(SeriesProd(Laurent1,Laurent2), &
                      BijLaurent(jpart,lpart,:))
      end do
    end do
  end do
!
  call InnloUser(p,yiQ,YikQ,         &
                 smeB,BornLaurent,   &
                 Bij,BijLaurent,     &
                 Bijkl,BijklLaurent, &
                 Innlofit,InnlofitLaurent)
! Giving back a factor of 1/(2*pi)^2:
  InnlofitLaurent = InnlofitLaurent/(2*pi)**2
!
!*********************************
! Reinstating scale dependence and r(eps):
!*********************************
!
! A factor of 1/(2pi) was factored out, giving it back...
  Laurent1 = InloLaurent/(2*pi)
!
  InloLaurent(-2) = Laurent1(-2)
  InloLaurent(-1) = Laurent1(-1) &
                  + Laurent1(-2)*L
  InloLaurent( 0) = Laurent1( 0) &
                  + Laurent1(-1)*L &
                  + Laurent1(-2)*L**2/2d0 &
                  - Laurent1(-2)*pi**2/12d0
  InloLaurent( 1) = Laurent1( 1) &
                  + Laurent1( 0)*L &
                  + Laurent1(-1)*L**2/2d0 &
                  - Laurent1(-1)*pi**2/12d0 &
                  + Laurent1(-2)*L**3/6d0 &
                  - Laurent1(-2)*L*pi**2/12d0 &
                  - Laurent1(-2)*Zeta3/3d0
  InloLaurent( 2) = Laurent1( 2) &
                  + Laurent1( 1)*L &
                  + Laurent1( 0)*L**2/2d0 &
                  - Laurent1( 0)*pi**2/12d0 &
                  + Laurent1(-1)*L**3/6d0 &
                  - Laurent1(-1)*L*pi**2/12d0 &
                  - Laurent1(-1)*Zeta3/3d0 &
                  + Laurent1(-2)*L**4/24d0 &
                  - Laurent1(-2)*L**2*pi**2/24d0 &
                  + Laurent1(-2)*pi**4/1440d0 &
                  - Laurent1(-2)*L*Zeta3/3d0
!
! We give back the finite term separately:
  Inloterm = InloLaurent(0)
!
! A factor of 1/(2*pi)^2 was factored out, giving it back...
  Laurent1 = I1I1Laurent/(2*pi)**2
! I_1^{(0)} I_1^{(0)}\times Born:
  I1I1Laurent(-4) = Laurent1(-4)
  I1I1Laurent(-3) = Laurent1(-3) &
                  + 2*Laurent1(-4)*L
  I1I1Laurent(-2) = Laurent1(-2) &
                  + 2*Laurent1(-3)*L &
                  + 2*Laurent1(-4)*L**2 &
                  - Laurent1(-4)*pi**2/6d0
  I1I1Laurent(-1) = Laurent1(-1) &
                  + 2*Laurent1(-2)*L &
                  + 2*Laurent1(-3)*L**2 &
                  - Laurent1(-3)*pi**2/6d0 &
                  + 4*Laurent1(-4)*L**3/3d0 &
                  - Laurent1(-4)*L*pi**2/3d0 &
                  - 2*Laurent1(-4)*Zeta3/3d0
  I1I1Laurent( 0) = Laurent1( 0) &
                  + 2*Laurent1(-1)*L &
                  + 2*Laurent1(-2)*L**2 &
                  - Laurent1(-2)*pi**2/6d0 &
                  + 4*Laurent1(-3)*L**3/3d0 &
                  - Laurent1(-3)*L*pi**2/3d0 &
                  - 2*Laurent1(-3)*Zeta3/3d0 &
                  + 2*Laurent1(-4)*L**4/3d0 &
                  - Laurent1(-4)*L**2*pi**2/3d0 &
                  + Laurent1(-4)*pi**4/120d0 &
                  - 4*Laurent1(-4)*L*Zeta3/3d0
!
! I_1^{(0)}\times V:
! To get correct scale dependence we have to multiply by
! S_\eps/S_\eps^{\bar{\rm MS}} (\mu^2/Q^2)^\eps, where
! S_\eps^{\bar{\rm MS}} = (4\pi)^\eps \ue^{-\eps\gamma_E}
! S_\eps = (4\pi)^\eps/Gamma(1 - \eps):
! Giving back a 1/(2*pi) factor:
  I1VLaurent = I1VLaurent/(2*pi)
  Laurent1 = I1VLaurent
  I1VLaurent(-4) = Laurent1(-4)
  I1VLaurent(-3) = Laurent1(-3) &
                 + Laurent1(-4)*L
  I1VLaurent(-2) = Laurent1(-2) &
                 + Laurent1(-3)*L &
                 + Laurent1(-4)*L**2/2d0 &
                 - Laurent1(-4)*pi**2/12d0
  I1VLaurent(-1) = Laurent1(-1) &
                 + Laurent1(-2)*L &
                 + Laurent1(-3)*L**2/2d0 &
                 - Laurent1(-3)*pi**2/12d0 &
                 + Laurent1(-4)*L**3/6d0 &
                 - Laurent1(-4)*L*pi**2/12d0 &
                 - Laurent1(-4)*Zeta3/3d0
  I1VLaurent( 0) = Laurent1( 0) &
                 + Laurent1(-1)*L &
                 + Laurent1(-2)*L**2/2d0 &
                 - Laurent1(-2)*pi**2/12d0 &
                 + Laurent1(-3)*L**3/6d0 &
                 - Laurent1(-3)*L*pi**2/12d0 &
                 - Laurent1(-3)*Zeta3/3d0 &
                 + Laurent1(-4)*L**4/24d0 &
                 - Laurent1(-4)*L**2*pi**2/24d0 &
                 + Laurent1(-4)*pi**4/1440d0 &
                 - Laurent1(-4)*L*Zeta3/3d0
!
  Laurent1 = InnlofitLaurent
!
  InnlofitLaurent(-4) = Laurent1(-4)
  InnlofitLaurent(-3) = Laurent1(-3) &
                      + 2*Laurent1(-4)*L
  InnlofitLaurent(-2) = Laurent1(-2) &
                      + 2*Laurent1(-3)*L &
                      + 2*Laurent1(-4)*L**2 &
                      - Laurent1(-4)*pi**2/6d0
  InnlofitLaurent(-1) = Laurent1(-1) &
                      + 2*Laurent1(-2)*L &
                      + 2*Laurent1(-3)*L**2 &
                      - Laurent1(-3)*pi**2/6d0 &
                      + 4*Laurent1(-4)*L**3/3d0 &
                      - Laurent1(-4)*L*pi**2/3d0 &
                      - 2*Laurent1(-4)*Zeta3/3d0
  InnlofitLaurent( 0) = Laurent1( 0) &
                      + 2*Laurent1(-1)*L &
                      + 2*Laurent1(-2)*L**2 &
                      - Laurent1(-2)*pi**2/6d0 &
                      + 4*Laurent1(-3)*L**3/3d0 &
                      - Laurent1(-3)*L*pi**2/3d0 &
                      - 2*Laurent1(-3)*Zeta3/3d0 &
                      + 2*Laurent1(-4)*L**4/3d0 &
                      - Laurent1(-4)*L**2*pi**2/3d0 &
                      + Laurent1(-4)*pi**4/120d0 &
                      - 4*Laurent1(-4)*L*Zeta3/3d0
!
  Innlofit = InnlofitLaurent(0)
!
  Innloterm = InnlofitLaurent(0) + I1I1Laurent(0) &
            + I1VLaurent(0) - qcd_beta0/2d0*InloLaurent(1)/(2*pi)
  if (present(InnloLaurent)) then
    do i=-4,2
      InnloLaurent(i) = InnlofitLaurent(i) + I1I1Laurent(i) &
                      + I1VLaurent(i)
      if (i.lt.2) then
        InnloLaurent(i) = InnloLaurent(i) &
                        - qcd_beta0/2d0*InloLaurent(i+1)/(2*pi)
      end if
    end do
  end if
!
end subroutine CalcInnlo
!
! This routine is the d-dimensional equivalent of CalcI1nlo,
! it uses d-dimensional Born to construct the 
! I_1^{(0)}\times\ud\sigma^{B} term.
subroutine CalcI1nloddim(p,BLaurent,BijLaurent,I1term,I1Laurent)
use IopInvariants
use particles
use QCDparams
use math
use misc
!
use scales
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
  real(kind(1d0)) , dimension(-4:2) , intent(in) :: BLaurent
  real(kind(1d0)) , dimension(:,:,-4:) , intent(in) :: BijLaurent
  real(kind(1d0)) , intent(out) :: I1term
  real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: I1Laurent
!
  integer :: ipart,kpart,npart
  real(kind(1d0)) :: Ccont,Scont
  real(kind(1d0)) , dimension(-4:2) :: CLaurent,SLaurent,Laurent
!
! DEBUG:
  real(kind(1d0)) :: mur_tmp
!
  interface
    subroutine CalcC1i0F(x,fi,Ccont,CLaurent)
    use particles
    implicit none
!
      real(kind(1d0)) , intent(in) :: x
      character , intent(in) :: fi
      real(kind(1d0)) , intent(out) :: Ccont
      real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: CLaurent
!
    end subroutine CalcC1i0F
!
    subroutine CalcS10ikFF(Y,Scont,SLaurent)
    use particles
    implicit none
!
      real(kind(1d0)) , intent(in) :: Y
      real(kind(1d0)) , intent(out) :: Scont
      real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: &
        SLaurent
!
    end subroutine CalcS10ikFF
  end interface
!
  npart = size(p)
!
  I1term = 0
  if (present(I1Laurent)) I1Laurent = 0
! When d-dimensional Born is present we always have to calculate the
! Laurent series of the I1xB, hence a local array is introduced to
! hold the Laurent series of the temporary result:
  Laurent = 0
! We calculate the invariants:
  call calcIopInvariants(p)
! We loop over all particles:
  do ipart=1,npart
! But only considering massless partons:
    if (abs(p(ipart)%flv).gt.qcd_nf) cycle
! We differentiate between two scenarios: Initial and Final
    if (ipart.le.2) then
! Initial State:
      Print *,"Initial state is not implemented yet in CalcI1ddim..."
      stop
    else
! Final State:
! (anti)Quark:
      if (p(ipart)%flv.ne.0) then
        call CalcC1i0F(yiQ(ipart),'q',Ccont,CLaurent)
        Laurent = Laurent &
                + 1d0/(2*pi)*qcd_cf*SeriesProd(CLaurent,BLaurent)
! Gluon:
      else
        call CalcC1i0F(yiQ(ipart),'g',Ccont,CLaurent)
        Laurent = Laurent &
                + 1d0/(2*pi)*qcd_ca*SeriesProd(CLaurent,BLaurent)
      end if
    end if
  end do
! We consider the contributions coming from the soft integrated
! counterterms:
  do ipart=1,npart
    if (abs(p(ipart)%flv).gt.qcd_nf) cycle
    do kpart=ipart+1,npart
      if (abs(p(kpart)%flv).gt.qcd_nf) cycle
      call CalcS10ikFF(YikQ(ipart,kpart),Scont,SLaurent)
      Laurent = Laurent + 2d0/(2*pi) &
              * SeriesProd(SLaurent,BijLaurent(ipart,kpart,:))
    end do
  end do
!
  I1term = Laurent(0)
!
! Additional term coming from the eps expansion of Seps/SepsMSbar:
  I1term = I1term - pi**2/12d0*Laurent(-2)
!
  if (present(I1Laurent)) I1Laurent = Laurent
!
end subroutine CalcI1nloddim
