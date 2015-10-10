! This source contains the momentum mappings and the explicit
! calculations for all the subtraction terms:
!
! This module holds all the subtractions needed in an NLO
! calculation:
module nlo_subtractions
implicit none
!
contains
! Final-Final Collinear subtraction term:
subroutine CalcCirFF(p,ptilde,smeB,Bmunu,emit,radi,radr, &
                     CirTerm)
use particles
use regions
use math
use QCDparams
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
  type(particle) , dimension(:) , intent(in) :: ptilde
  real(kind(1d0)) , intent(in) :: smeB
  real(kind(1d0)) , dimension(0:3,0:3) , intent(in) :: Bmunu
  integer , intent(in) :: emit,radi,radr
  real(kind(1d0)) , intent(out) :: CirTerm
!
  integer :: d0,m
!
  real(kind(1d0)) , external :: Pqg0,Pgg0,Pqq0
!
!
!
  CirTerm = 0d0
!
  d0 = sub_d0
  m  = size(ptilde) - 3
!
! q(~) -> q(~) + g
  if ((p(radi)%flv.ne.0).and.(p(radr)%flv.eq.0)) then
!    print *,"q(~) -> q(~) + g"
    CirTerm = 8d0*pi/sir*Pqg0(ztildei,ztilder)*smeB
! g -> g + g
  elseif ((p(radi)%flv.eq.0).and.(p(radr)%flv.eq.0)) then 
!    print *,"g -> g + g"
    CirTerm = 8d0*pi/sir*Pgg0(ztildei,ztilder,kttildei,Bmunu)
! g -> q + q~ || g -> q~ + q
  elseif ((p(radi)%flv.ne.0).and.(p(radr)%flv.ne.0)) then 
!    print *,"g -> q + q~"
    CirTerm = 8d0*pi/sir*Pqq0(ztildei,ztilder,kttildei,Bmunu)
  else
  end if
  CirTerm = CirTerm * (1d0 - alphair)**(2*d0 - 2*(m - 0*1))
!
!
end subroutine CalcCirFF
!
subroutine CalcSr(p,ptilde,Bij,radr,SrTerm)
use particles
use regions
use math
use QCDparams
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
  type(particle) , dimension(:) , intent(in) :: ptilde
  real(kind(1d0)) , dimension(:,:) , intent(in) :: Bij
  integer , intent(in) :: radr
  real(kind(1d0)) , intent(out) :: SrTerm
!
  integer :: i,k,r
  integer :: ipart,kpart
  integer :: d0pr,m
!
  real(kind(1d0)) :: Sik
!
  real(kind(1d0)) :: smeR
!
!
  Sik(i,k,r) = 2d0*sir_arr(i,k)/sir_arr(i,r)/sir_arr(k,r)
!
  SrTerm = 0d0
!
  d0pr = sub_d0pr
  m  = size(ptilde) - 3
!
  r = radr
!
! We have to loop over all partons:
! The 1/2 from the front is missing since the summation is ordered:
! The summation runs over the momenta defining the real-emission PS,
! omitting the one corresponding to the momentum going soft (r).
! In order to obtain the corresponding color-ordered SME we introduced
! two separated, related running indices: ipart and kpart which are
! obtained by mapping the real momenta to the underlying Born ones:
  ipart = 0
  do i=1,size(p)-1
    if (i.eq.r) cycle
    ipart = ipart + 1
    if (abs(p(i)%flv).gt.6) cycle
    kpart = ipart
    do k=i+1,size(p)
      if (k.eq.r) cycle
      kpart = kpart + 1
      if (abs(p(k)%flv).gt.6) cycle
      SrTerm = SrTerm + Sik(i,k,r) * Bij(ipart,kpart)
    end do
  end do
!
  SrTerm = -8d0*pi * SrTerm
!
  SrTerm = SrTerm * (1d0 - yrQ)**(d0pr - (m - 0*1))
!
end subroutine CalcSr
!
subroutine CalcCirFFSr(p,ptilde,smeB,emit,radi,radr,CSirTerm)
use particles
use regions
use math
use QCDparams
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
  type(particle) , dimension(:) , intent(in) :: ptilde
  real(kind(1d0)) , intent(in) :: smeB
  integer , intent(in) :: emit,radi,radr
  real(kind(1d0)) , intent(out) :: CSirTerm
!
  integer :: d0pr,m
  real(kind(1d0)) :: T2smeB
!
!
  CSirTerm = 0d0
!
  d0pr = sub_d0pr
  m  = size(ptilde) - 3
!
! We have to calculate the momentum fractions:
  yirQ = yiQ_arr(radi) + yiQ_arr(radr)
  ztildei = yiQ_arr(radi)/yirQ
  ztilder = 1d0 - ztildei
  sir     = sir_arr(radi,radr)
!
  if (abs(p(emit)%flv).ne.0) T2smeB = qcd_cf*smeB
  if (abs(p(emit)%flv).eq.0) T2smeB = qcd_ca*smeB
!
  CSirTerm = 8d0*pi / sir*2d0*ztildei/ztilder * T2smeB
!
  CSirTerm = CSirTerm * (1d0 - yrQ)**(d0pr - (m - 0*1))
!
end subroutine CalcCirFFSr
!
end module nlo_subtractions
!
module nnlo_subtractions
implicit none
!
!
contains
subroutine CalcCir01FF(p,ptilde,emit,radi,radr,CirTerm,CirLaurent)
use particles
use regions
use math
use QCDparams
use misc
use process
use scales
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
  type(particle) , dimension(:) , intent(in) :: ptilde
  integer , intent(in) :: emit,radi,radr
  real(kind(1d0)) , intent(out) :: CirTerm
  real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: CirLaurent
!
  integer :: d0,m
  real(kind(1d0)) :: smeV
  real(kind(1d0)) , dimension(-4:2) :: smeVlaurent
  real(kind(1d0)) , dimension(0:3,0:3) :: Vmunu
  real(kind(1d0)) , dimension(:,:,:) , allocatable :: VmunuLaurent
!
  real(kind(1d0)) , external :: Pqg0,Pgg0,Pqq0
!
  interface 
    subroutine CalcV(p,smeV,smeVLaurent)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: smeV
      real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: smeVLaurent
!
    end subroutine CalcV
!
    subroutine CalcVmunu(ileg,p,Vmunu,VmunuLaurent)
    use particles
    implicit none
!
      integer , intent(in) :: ileg
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(0:3,0:3) , intent(out) :: Vmunu
      real(kind(1d0)) , optional , dimension(0:3,0:3,-4:2) , intent(out) :: VmunuLaurent
!
    end subroutine CalcVmunu
  end interface
!
  CirTerm = 0d0
!
  d0 = sub_d0
  m  = size(ptilde) - 2
!
! We only allocate the array for the Laurent expansion if it is
! needed:
  if (present(CirLaurent)) allocate(VmunuLaurent(0:3,0:3,-4:2))
!
! q(~) -> q(~) + g
  if ((p(radi)%flv.ne.0).and.(p(radr)%flv.eq.0)) then
!    print *,"q(~) -> q(~) + g"
    if (.not.present(CirLaurent)) then
      call CalcV(ptilde,smeV)
    else
      call CalcV(ptilde,smeV,smeVlaurent)
    end if
    CirTerm = 8d0*pi/sir*Pqg0(ztildei,ztilder)*smeV
    if (present(CirLaurent)) CirLaurent = &
      8d0*pi/sir*Pqg0(ztildei,ztilder)*smeVlaurent
! g -> g + g
  elseif ((p(radi)%flv.eq.0).and.(p(radr)%flv.eq.0)) then 
!    print *,"g -> g + g"
    if (.not.present(CirLaurent)) then
      call CalcVmunu(emit,ptilde,Vmunu)
    else
      call CalcVmunu(emit,ptilde,Vmunu,VmunuLaurent)
    end if
    CirTerm = 8d0*pi/sir*Pgg0(ztildei,ztilder,kttildei,Vmunu)
!
    if (present(CirLaurent)) then
      CirLaurent = 0d0
      CirLaurent(-2) = 8d0*pi/sir*Pgg0(ztildei,ztilder,kttildei, &
                                       VmunuLaurent(:,:,-2))
      CirLaurent(-1) = 8d0*pi/sir*Pgg0(ztildei,ztilder,kttildei, &
                                       VmunuLaurent(:,:,-1))
    end if
! g -> q + q~ || g -> q~ + q
  elseif ((p(radi)%flv.ne.0).and.(p(radr)%flv.ne.0).and. &
          (p(radi)%flv.eq.-p(radr)%flv)) then 
!    print *,"g -> q + q~"
    if (.not.present(CirLaurent)) then
      call CalcVmunu(emit,ptilde,Vmunu)
    else
      call CalcVmunu(emit,ptilde,Vmunu,VmunuLaurent)
    end if
!
    CirTerm = 8d0*pi/sir*Pqq0(ztildei,ztilder,kttildei,Vmunu)
!
    if (present(CirLaurent)) then
      CirLaurent = 0d0
      CirLaurent(-2) = 8d0*pi/sir*Pqq0(ztildei,ztilder,kttildei, &
                                       VmunuLaurent(:,:,-2))
      CirLaurent(-1) = 8d0*pi/sir*Pqq0(ztildei,ztilder,kttildei, &
                                       VmunuLaurent(:,:,-1))
    end if
  else
  end if
  CirTerm = CirTerm * (1d0 - alphair)**(2*d0 - 2*(m - 1))
!
  if (present(CirLaurent)) deallocate(VmunuLaurent)
!
end subroutine CalcCir01FF
!
subroutine CalcCir10FF(p,ptilde,smeB,Bmunu,emit,radi,radr,CirTerm,CirLaurent)
use particles
use regions
use math
use QCDparams
use scales
use misc
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
  type(particle) , dimension(:) , intent(in) :: ptilde
  real(kind(1d0)) , intent(in) :: smeB
  real(kind(1d0)) , dimension(0:3,0:3) , intent(in) :: Bmunu
  integer , intent(in) :: emit,radi,radr
  real(kind(1d0)) , intent(out) :: CirTerm
  real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: CirLaurent
!
  integer :: d0,m
!
!
  CirTerm = 0d0
!
  d0 = sub_d0
  m  = size(ptilde) - 2
!
! q(~) -> q(~) + g
  if ((p(radi)%flv.ne.0).and.(p(radr)%flv.eq.0)) then
!    print *,"q(~) -> q(~) + g"
    if (.not.present(CirLaurent)) then
      call CalcCqg10FF(ztildei,ztilder,sir,mur**2,CirTerm)
    else
      call CalcCqg10FF(ztildei,ztilder,sir,mur**2,CirTerm,CirLaurent)
      CirLaurent = CirLaurent * smeB
    end if
    CirTerm = CirTerm * smeB
! g -> g + g
  elseif ((p(radi)%flv.eq.0).and.(p(radr)%flv.eq.0)) then 
!    print *,"g -> g + g"
    if (.not.present(CirLaurent)) then
      call CalcCgg10FF(ztildei,ztilder,kttildei,Bmunu,sir,mur**2, &
                       CirTerm)
    else
      call CalcCgg10FF(ztildei,ztilder,kttildei,Bmunu,sir,mur**2, &
                       CirTerm,CirLaurent)
    end if
    CirTerm = CirTerm
! g -> q + q~ || g -> q~ + q
  elseif ((p(radi)%flv.ne.0).and.(p(radr)%flv.ne.0)) then 
!    print *,"g -> q + q~"
    if (.not.present(CirLaurent)) then
      call CalcCqq10FF(ztildei,ztilder,kttildei,Bmunu,sir,mur**2, &
                       CirTerm)
    else
      call CalcCqq10FF(ztildei,ztilder,kttildei,Bmunu,sir,mur**2, &
                       CirTerm,CirLaurent)
    end if
    CirTerm = CirTerm
  else
  end if
  CirTerm = CirTerm * (1d0 - alphair)**(2*d0 - 2*(m - 1))
!
end subroutine CalcCir10FF
!
subroutine CalcCqg10FF(zi,zr,sir,mu2,Cir10Term,Cir10Laurent)
use QCDparams
use math
implicit none
!
  real(kind(1d0)) , intent(in) :: zi,zr
  real(kind(1d0)) , intent(in) :: sir,mu2
  real(kind(1d0)) , intent(out) :: Cir10Term
  real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: Cir10Laurent
!
  real(kind(1d0)) , external :: Pqg0,ddilog
!
  Cir10Term = &
    cGamma*(8d0*pi)**2*(qcd_ca*qcd_cf-qcd_cf**2+  &
    Pqg0(zi,zr)*((qcd_ca*(2*pi**2-12*ddilog(-(zr/zi))-3*log(mu2/sir)**2  &
    -6*log(mu2/sir)*log(zi/zr)-3*log(zi/zr)**2))/6.d0+  &
    (qcd_cf*(12*ddilog(-(zr/zi))-12*log(mu2/sir)*log(1+ &
    zr/zi)))/6.d0))
! The eps poles are only calculated if needed:
  if (present(Cir10Laurent)) then
    Cir10Laurent = 0d0
    Cir10Laurent( 0) = Cir10Term
    Cir10Laurent(-1) = &
      -(qcd_beta0*(8d0*pi)**2*Pqg0(zi,zr))/(3.2d1*pi**2)+  &
      cGamma*(8d0*pi)**2*Pqg0(zi,zr)*(qcd_ca*(-log(mu2/sir)-  &
      log(zi/zr))-2*qcd_cf*log(1+zr/zi))
    Cir10Laurent(-2) = -qcd_ca*cGamma*(8d0*pi)**2*Pqg0(zi,zr)
  end if
  Cir10Term = Cir10Term / sir
  if (present(Cir10Laurent)) Cir10Laurent = Cir10Laurent / sir
!
end subroutine CalcCqg10FF
!
subroutine CalcCgg10FF(zi,zr,kt,Bmunu,sir,mu2,Cir10Term,Cir10Laurent)
use particles
use QCDparams
use math
implicit none
!
  real(kind(1d0)) , intent(in) :: zi,zr
  type(particle) , intent(in) :: kt
  real(kind(1d0)) , intent(in) :: sir,mu2
  real(kind(1d0)) , dimension(0:3,0:3) , intent(in) :: Bmunu
  real(kind(1d0)) , intent(out) :: Cir10Term
  real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: Cir10Laurent
!
  real(kind(1d0)) :: kt2,P0gg
!
  real(kind(1d0)) , external :: Pgg0
!
  kt2 = kt%p*kt%p
  P0gg = Pgg0(zi,zr,kt,Bmunu)
!
  Cir10Term = &
    qcd_ca*cGamma*(8d0*pi)**2*((-2*(kt%p*Bmunu*kt%p)*(qcd_ca-  &
    2*qcd_nf*qcd_tr))/(3.d0*kt2)+(P0gg*(2*pi**2-  &
    3*log(mu2/sir)**2-6*log(mu2/sir)*log(zi/zr)-3*log(zi/zr)**2+ &
    12*log(mu2/sir)*log(zi)))/6.d0)
! The eps poles are only calculated if needed:
  if (present(Cir10Laurent)) then
    Cir10Laurent = 0d0
    Cir10Laurent( 0) = Cir10Term
    Cir10Laurent(-1) = &
      -(qcd_beta0*(8d0*pi)**2*P0gg*Seps)/(3.2d1*pi**2)+  &
      qcd_ca*cGamma*(8d0*pi)**2*P0gg*(-log(mu2/sir)-log(zi/zr)+  &
      2*log(zi))
    Cir10Laurent(-2) = -qcd_ca*cGamma*(8d0*pi)**2*P0gg
  end if
  Cir10Term = Cir10Term / sir
  if (present(Cir10Laurent)) Cir10Laurent = Cir10Laurent / sir
!
end subroutine CalcCgg10FF
!
subroutine CalcCqq10FF(zi,zr,kt,Bmunu,sir,mu2,Cir10Term,Cir10Laurent)
use particles
use QCDparams
use math
implicit none
!
  real(kind(1d0)) , intent(in) :: zi,zr
  type(particle) , intent(in) :: kt
  real(kind(1d0)) , intent(in) :: sir,mu2
  real(kind(1d0)) , dimension(0:3,0:3) , intent(in) :: Bmunu
  real(kind(1d0)) , intent(out) :: Cir10Term
  real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: Cir10Laurent
!
  real(kind(1d0)) :: P0qq
!
  real(kind(1d0)) , external :: Pqq0,ddilog
!
  P0qq = Pqq0(zi,zr,kt,Bmunu)
!
  Cir10Term = &
    qcd_ca*cGamma*(8d0*pi)**2*P0qq*(-(qcd_nf*(20+  &
    12*log(mu2/sir)))/(1.8d1*qcd_nc)-(-72+9*pi**2-  &
    27*log(mu2/sir)-  &
    9*log(mu2/sir)**2)/(1.8d1*qcd_nc**2)+(80-3*pi**2+ &
    39*log(mu2/sir)+36*log(mu2/sir)*log(zr)- &
    18*log(mu2/sir)*log(zr/zi)- &
    9*log(zr/zi)**2)/1.8d1)
! The eps poles are only calculated if needed:
  if (present(Cir10Laurent)) then
    Cir10Laurent = 0d0
    Cir10Laurent( 0) = Cir10Term
    Cir10Laurent(-1) = &
      -(qcd_beta0*(8d0*pi)**2*P0qq*Seps)/(3.2d1*pi**2)+  &
      cGamma*(8d0*pi)**2*P0qq*((-4*qcd_nf*qcd_tr)/3.d0+qcd_cf*(-3-  &
      2*log(mu2/sir))+qcd_ca*(11d0/3.d0+log((mu2*zi)/sir)+ &
      log(zr)))
    Cir10Laurent(-2) = (qcd_ca-2*qcd_cf)*cGamma*(8d0*pi)**2*P0qq
  end if
  Cir10Term = Cir10Term / sir
  if (present(Cir10Laurent)) Cir10Laurent = Cir10Laurent / sir
!
end subroutine CalcCqq10FF
!
subroutine CalcSr01(p,ptilde,Vij,radr,SrTerm,SrLaurent)
use particles
use regions
use math
use QCDparams
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
  type(particle) , dimension(:) , intent(in) :: ptilde
  real(kind(1d0)) , dimension(:,:) , intent(inout) :: Vij
  integer , intent(in) :: radr
  real(kind(1d0)) , intent(out) :: SrTerm
  real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: SrLaurent
!
  integer :: i,k,r,n
  integer :: ipart,kpart
  integer :: d0pr,m
  real(kind(1d0)) , dimension(:,:,:) , allocatable :: VijLaurent
!
  real(kind(1d0)) :: Sik
!
  interface
    subroutine CalcVij(p,Vij,VijLaurent)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Vij
      real(kind(1d0)) , optional , dimension(:,:,-4:) , intent(out) :: VijLaurent
!
    end subroutine CalcVij
  end interface
!
  Sik(i,k,r) = 2d0*sir_arr(i,k)/sir_arr(i,r)/sir_arr(k,r)
!
  SrTerm = 0d0
!
  d0pr = sub_d0pr
  m  = size(ptilde) - 2
!
  n = size(p)
!
  r = radr
!
! For details on the next part take a look at the comments in 
! CalcSr:
  ipart = 0
  do i=1,n
    if (i.eq.r) cycle
    ipart = ipart + 1
    if (abs(p(i)%flv).gt.6) cycle
    kpart = ipart
    do k=i+1,n
      if (k.eq.r) cycle
      kpart = kpart + 1
      if (abs(p(k)%flv).gt.6) cycle
      SrTerm = SrTerm + Sik(i,k,r) * Vij(ipart,kpart)
    end do
  end do
!
  SrTerm = -8d0*pi * SrTerm * (1d0 - yrQ)**(d0pr - (m - 1))
!
  if (present(SrLaurent)) then
! Only allocate if needed:
    allocate(VijLaurent(n-1,n-1,-4:2))
    call CalcVij(ptilde,Vij,VijLaurent)
    SrLaurent = 0d0
    ipart = 0
    do i=1,n
      if (i.eq.r) cycle
      ipart = ipart + 1
      if (abs(p(i)%flv).gt.6) cycle
      kpart = ipart
      do k=i+1,n
        if (k.eq.r) cycle
        kpart = kpart + 1
        if (abs(p(k)%flv).gt.6) cycle
        SrLaurent(-2:-1) = SrLaurent(-2:-1) &
                  + Sik(i,k,r)*VijLaurent(ipart,kpart,-2:-1)
      end do
    end do
    SrLaurent = -8d0*pi*SrLaurent
    deallocate(VijLaurent)
  end if
!
end subroutine CalcSr01
!
subroutine CalcSr10(p,ptilde,Bij,Bijk,radr,SrTerm,SrLaurent)
use flags
use particles
use regions
use math
use QCDparams
use scales
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
  type(particle) , dimension(:) , intent(in) :: ptilde
  real(kind(1d0)) , dimension(:,:) , intent(in) :: Bij
  real(kind(1d0)) , dimension(:,:,:) , intent(out) :: Bijk
  integer , intent(in) :: radr
  real(kind(1d0)) , intent(out) :: SrTerm
  real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: SrLaurent
!
  integer :: i,k,l,r
  integer :: ipart,kpart,lpart
  integer :: numcol
  integer :: d0pr,m
  integer :: n
  real(kind(1d0)) :: mu2
!
  real(kind(1d0)) :: Sik
!
  interface
    subroutine CalcBijk(p,Bijk)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:,:,:) , intent(out) :: Bijk
!
    end subroutine CalcBijk
  end interface
!
  Sik(i,k,r) = 2d0*sir_arr(i,k)/sir_arr(i,r)/sir_arr(k,r)
!
! When the routine is called the first time a small initialization
! is allowed:
  if (init_CalcSr10) then
! Counting down colored particles at the Born level:
    numcol = 0
    do i=1,size(ptilde)
      if (abs(ptilde(i)%flv).le.6) numcol = numcol + 1
    end do
! If the number of colored particles at the Born level is
! less than 4 the contribution proportional to Bijk is not
! needed provided it is zero:
    flg_Bijk = .true.
    if (numcol.lt.4) flg_Bijk = .false.
!
    init_CalcSr10 = .false.
  end if
!
  SrTerm = 0d0
  if (present(SrLaurent)) SrLaurent = 0d0
!
  d0pr = sub_d0pr
  m  = size(ptilde) - 2
!
  r = radr
!
  n = size(p)
!
  mu2 = mur*mur
!
  if (flg_Bijk) call CalcBijk(ptilde,Bijk)
!
! For details on the next part take a look at the comments in 
! CalcSr:
  ipart = 0
! i runs through momenta in the real-emission-like configuration:
  do i=1,n
    if (i.eq.r) cycle
! ipart is the corresponding momenutm in the ptilde array:
    ipart = ipart + 1
! Only partons matter:
    if (abs(p(i)%flv).gt.6) cycle
    kpart = ipart
! Ordered pairs are created in (i,k):
    do k=i+1,n
      if (k.eq.r) cycle
! kpart is the correspodent for k used with ptilde:
      kpart = kpart + 1
! Only partons!
      if (abs(p(k)%flv).gt.6) cycle
      SrTerm = SrTerm + &
        qcd_ca*Sik(i,k,r)*(-2*pi**2+  &
        3*(-log(mu2/sir_arr(i,k))+log(mu2/sir_arr(i,r))+  &
        log(mu2/sir_arr(k,r)))**2)/6d0*Bij(ipart,kpart)
      if (present(SrLaurent)) then
        SrLaurent(-2) = SrLaurent(-2) &
                      + cGamma*qcd_ca*Sik(i,k,r)*Bij(ipart,kpart)
        SrLaurent(-1) = SrLaurent(-1) &
                      + Sik(i,k,r)*((qcd_beta0*Seps)/(3.2d1*pi**2) &
                      + qcd_ca*cGamma*(-log(mu2/sir_arr(i,k)) &
                      + log(mu2/sir_arr(i,r)) &
                      + log(mu2/sir_arr(k,r))))*Bij(ipart,kpart)
      end if
! This part is only evaluated if Bijk is not identically zero:
      if (flg_Bijk) then
        lpart = 0
        do l=1,n-1
          if (l.eq.r) cycle
          lpart = lpart + 1
          if ((l.eq.i).or.(l.eq.k)) cycle
          if (abs(p(l)%flv).gt.6) cycle
          SrTerm = SrTerm + &
            2*pi*Sik(i,k,r)*(log(mu2/sir_arr(k,l))-  &
            log(mu2/sir_arr(k,r))-  &
            log(mu2/sir_arr(l,r)))*Bijk(ipart,kpart,lpart)
          if (present(SrLaurent)) then
            SrLaurent(-1) = SrLaurent(-1) &
                          - 2*cGamma*pi*Sik(i,k,r) &
                          * Bijk(ipart,kpart,lpart)
          end if
        end do
      end if
    end do
  end do
  SrTerm = cGamma*(8d0*pi)**2*SrTerm * (1d0 - yrQ)**(d0pr - (m - 1))
!
! Note, that cGamma cannot be factored out due to the renormalization
! related term:
  if (present(SrLaurent)) SrLaurent = (8d0*pi)**2*SrLaurent
!
end subroutine CalcSr10
!
subroutine CalcCirFFSr01(p,ptilde,smeV,emit,radi,radr,CirSrTerm,CirSrLaurent)
use particles
use regions
use math
use QCDparams
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
  type(particle) , dimension(:) , intent(in) :: ptilde
  real(kind(1d0))  , intent(inout) :: smeV
  integer , intent(in) :: emit,radi,radr
  real(kind(1d0)) , intent(out) :: CirSrTerm
  real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: CirSrLaurent
!
  integer :: d0pr,m
  real(kind(1d0)) :: zi,zr,Ti2
  real(kind(1d0)) , dimension(-4:2) :: VirtLaurent
!
  interface
    subroutine CalcV(p,smeV,smeVLaurent)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: smeV
      real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: smeVLaurent
!
    end subroutine CalcV
  end interface
!
  CirSrTerm = 0d0
!
  d0pr = sub_d0pr
  m  = size(ptilde) - 2
!
  if (present(CirSrLaurent)) then
    call CalcV(ptilde,smeV,VirtLaurent)
  end if
! We have to calculate the momentum fractions:
  yirQ = yiQ_arr(radi) + yiQ_arr(radr)
  zi = yiQ_arr(radi)/yirQ
  zr = 1d0 - zi
  sir     = sir_arr(radi,radr)
!
  if (abs(p(emit)%flv).ne.0) Ti2 = qcd_cf
  if (abs(p(emit)%flv).eq.0) Ti2 = qcd_ca
!
  CirSrTerm = 8d0*pi / sir*2d0*zi/zr * Ti2*smeV &
            * (1d0 - yrQ)**(d0pr - (m - 1))
!
  if (present(CirSrLaurent)) then
    CirSrLaurent = 8d0*pi/sir &
                 * 2d0*zi/zr*Ti2*VirtLaurent
  end if
!
end subroutine CalcCirFFSr01
!
subroutine CalcCirFFSr10(p,ptilde,smeB,emit,radi,radr,CirSrTerm,CirSrLaurent)
use particles
use regions
use math
use QCDparams
use scales
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
  type(particle) , dimension(:) , intent(in) :: ptilde
  real(kind(1d0)) , intent(in) :: smeB
  integer , intent(in) :: emit,radi,radr
  real(kind(1d0)) , intent(out) :: CirSrTerm
  real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: CirSrLaurent
!
  integer :: d0pr,m
!
  real(kind(1d0)) :: mu2,zi,zr,Ti2
!
  CirSrTerm = 0d0
!
  d0pr = sub_d0pr
  m  = size(ptilde) - 2
!
  mu2 = mur*mur
!
! We have to calculate the momentum fractions:
  yirQ = yiQ_arr(radi) + yiQ_arr(radr)
  zi = yiQ_arr(radi)/yirQ
  zr = 1d0 - zi
  sir     = sir_arr(radi,radr)
!
  if (abs(p(emit)%flv).ne.0) Ti2 = qcd_cf
  if (abs(p(emit)%flv).eq.0) Ti2 = qcd_ca
!
  CirSrTerm = &
    (qcd_ca*cGamma*(8d0*pi)**2*smeB*Ti2*zi*(2*pi**2-  &
    3*(log(mu2/sir)+log(zi/zr))**2))/(3.d0*sir*zr) &
    *(1d0 - yrQ)**(d0pr - (m - 1))
!
  if (present(CirSrLaurent)) then
    CirSrLaurent = 0d0
    CirSrLaurent(-2) = &
      -2*qcd_ca*cGamma*(8d0*pi)**2*Ti2*smeB*zi/(sir*zr)
    CirSrLaurent(-1) = &
      (8d0*pi)**2*smeB*Ti2*zi*(-(qcd_beta0*Seps)/(1.6d1*pi**2)-  &
      2*qcd_ca*cGamma*(log(mu2/sir)+log(zi/zr)))/(sir*zr)
  end if
!
end subroutine CalcCirFFSr10
!
subroutine CalcRirF(flv_i,flv_r,flv_ir,zi,zr,yir,yirQ,Rir,RirLaurent)
use QCDparams
use math
implicit none
!
  integer , intent(in) :: flv_i,flv_r,flv_ir
  real(kind(1d0)) , intent(in) :: zi,zr
  real(kind(1d0)) , intent(in) :: yir,yirQ
  real(kind(1d0)) , intent(out) :: Rir
  real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: RirLaurent
!
  real(kind(1d0)) :: Ti2,Tr2,Tir2
  real(kind(1d0)) :: Ci,Cr,Cir,Sir
  real(kind(1d0)) :: Y
  real(kind(1d0)) , dimension(-4:2) :: CiLaurent,CrLaurent,CirLaurent,SirLaurent
!
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
  Rir = 0d0
!
!  print *,"yirQ: ",yirQ
!  print *,"zi*yirQ: ",zi*yirQ
!  print *,"zr*yirQ: ",zr*yirQ
!  print *,"Y: ",yir/(zi*zr*yirQ**2)
!
  if (.not.present(RirLaurent)) then
    if (flv_i.ne.0) then
      Ti2 = qcd_cf
      call CalcC1i0F(zi*yirQ,'q',Ci)
    else
      Ti2 = qcd_ca
      call CalcC1i0F(zi*yirQ,'g',Ci)
    end if
    if (flv_r.ne.0) then
      Tr2 = qcd_cf
      call CalcC1i0F(zr*yirQ,'q',Cr)
    else
      Tr2 = qcd_ca
      call CalcC1i0F(zr*yirQ,'g',Cr)
    end if
    if (flv_ir.ne.0) then
      Tir2 = qcd_cf
      call CalcC1i0F(yirQ,'q',Cir)
    else
      Tir2 = qcd_ca
      call CalcC1i0F(yirQ,'g',Cir)
    end if
  else
    if (flv_i.ne.0) then
      Ti2 = qcd_cf
      call CalcC1i0F(zi*yirQ,'q',Ci,CiLaurent)
    else
      Ti2 = qcd_ca
      call CalcC1i0F(zi*yirQ,'g',Ci,CiLaurent)
    end if
    if (flv_r.ne.0) then
      Tr2 = qcd_cf
      call CalcC1i0F(zr*yirQ,'q',Cr,CrLaurent)
    else
      Tr2 = qcd_ca
      call CalcC1i0F(zr*yirQ,'g',Cr,CrLaurent)
    end if
    if (flv_ir.ne.0) then
      Tir2 = qcd_cf
      call CalcC1i0F(yirQ,'q',Cir,CirLaurent)
    else
      Tir2 = qcd_ca
      call CalcC1i0F(yirQ,'g',Cir,CirLaurent)
    end if
  end if
!
  Y = yir/(zi*yirQ*zr*yirQ)
  if (.not.present(RirLaurent)) then
    call CalcS10ikFF(Y,Sir)
  else
    call CalcS10ikFF(Y,Sir,SirLaurent)
  end if
!
  Rir = 1d0/(2d0*pi)*(Ci*Ti2 + Cr*Tr2 - Cir*Tir2 &
      + (Tir2 - Ti2 - Tr2)*Sir)
  if (present(RirLaurent)) then
    RirLaurent = 1d0/(2d0*pi)*(CiLaurent*Ti2 + CrLaurent*Tr2  &
               - CirLaurent*Tir2 + (Tir2 - Ti2 - Tr2)*SirLaurent)
  end if
!
end subroutine CalcRirF
!
subroutine CalcRikrF(yik,yir,ykr,yiQ,ykQ,yrQ,Rikr,RikrLaurent)
use QCDparams
use math
implicit none
!
  real(kind(1d0)) , intent(in) :: yik,yir,ykr
  real(kind(1d0)) , intent(in) :: yiQ,ykQ,yrQ
  real(kind(1d0)) , intent(out) :: Rikr
  real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: RikrLaurent
!
  real(kind(1d0)) :: Y
  real(kind(1d0)) :: Cg,Sik,Sir,Skr
!
  real(kind(1d0)) , dimension(:) , allocatable :: CgLaurent, &
                                                  SikLaurent, &
                                                  SirLaurent, &
                                                  SkrLaurent
!
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
  Rikr = 0d0
!
  if (.not.present(RikrLaurent)) then
    call CalcC1i0F(yrQ,'g',Cg)
!
    Y = yik/(yiQ*ykQ)
    call CalcS10ikFF(Y,Sik)
    Y = yir/(yiQ*yrQ)
    call CalcS10ikFF(Y,Sir)
    Y = ykr/(ykQ*yrQ)
    call CalcS10ikFF(Y,Skr)
  else
    allocate(CgLaurent(-4:2), &
             SikLaurent(-4:2), &
             SirLaurent(-4:2), &
             SkrLaurent(-4:2))
    call CalcC1i0F(yrQ,'g',Cg,CgLaurent)
!
    Y = yik/(yiQ*ykQ)
    call CalcS10ikFF(Y,Sik,SikLaurent)
    Y = yir/(yiQ*yrQ)
    call CalcS10ikFF(Y,Sir,SirLaurent)
    Y = ykr/(ykQ*yrQ)
    call CalcS10ikFF(Y,Skr,SkrLaurent)
  end if
!
  Rikr = qcd_ca/(2d0*pi)*(Cg + Sik - Sir - Skr)
!
  if (present(RikrLaurent)) then
    RikrLaurent = qcd_ca/(2d0*pi)*(CgLaurent &
                + SikLaurent - SirLaurent - SkrLaurent)
    deallocate(CgLaurent,SikLaurent,SirLaurent,SkrLaurent)
  end if
!
end subroutine CalcRikrF
!
subroutine CalcCir00IFF(p,ptilde,smeB,Bij,Bmunu,Bmunuij,emit,radi,radr, &
                        CirTerm,CirLaurent)
use particles
use regions
use math
use QCDparams
use misc
use process
use scales
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
  type(particle) , dimension(:) , intent(in) :: ptilde
  real(kind(1d0)) , intent(in) :: smeB
  real(kind(1d0)) , dimension(:,:) , intent(in) :: Bij
  real(kind(1d0)) , dimension(0:3,0:3) , intent(in) :: Bmunu
  real(kind(1d0)) , dimension(0:,0:,:,:) , intent(in) :: Bmunuij
  integer , intent(in) :: emit,radi,radr
  real(kind(1d0)) , intent(out) :: CirTerm
  real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: CirLaurent
!
  integer :: d0,m
  integer :: jpart,kpart,ntilde
  real(kind(1d0)) :: Sjk,CS,Cj
  real(kind(1d0)) :: Yjk
  real(kind(1d0)) :: Tj2
  real(kind(1d0)) :: mu2
  real(kind(1d0)) :: P0qg,P0gg,P0qq
  real(kind(1d0)) , dimension(0:3,0:3) :: Bmunu_tmp
!
  real(kind(1d0)) , external :: Pqg0,Pgg0,Pqq0
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
    subroutine CalcCirSr0F(CScont,CSLaurent)
    implicit none
!
      real(kind(1d0)) , intent(out) :: CScont
      real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: CSLaurent
!
    end subroutine CalcCirSr0F
  end interface
!
  CirTerm = 0d0
!
  d0 = sub_d0
  m  = size(ptilde) - 2
!
  ntilde = size(ptilde)
!
! Calculation of the CS term which does not depend upon the
! kinematics:
  call CalcCirSr0F(CS)
!
! q(~) -> q(~) + g
  if ((p(radi)%flv.ne.0).and.(p(radr)%flv.eq.0)) then
!    print *,"q(~) -> q(~) + g"
    do jpart=1,ntilde
      if (abs(ptilde(jpart)%flv).gt.6) cycle
      if (ptilde(jpart)%flv.ne.0) then
        call CalcC1i0F(yitildeQ_arr(jpart),'q',Cj)
        Tj2 = qcd_cf
      else
        call CalcC1i0F(yitildeQ_arr(jpart),'g',Cj)
        Tj2 = qcd_ca
      end if
      CirTerm = CirTerm + (Cj + CS)*Tj2*smeB
      do kpart=jpart+1,ntilde
        if (abs(ptilde(kpart)%flv).gt.6) cycle
        Yjk = yirtilde_arr(jpart,kpart) &
            / (yitildeQ_arr(jpart)*yitildeQ_arr(kpart))
        call CalcS10ikFF(Yjk,Sjk)
        CirTerm = CirTerm + 2d0*(Sjk + CS)*Bij(jpart,kpart)
      end do
    end do
    CirTerm = Pqg0(ztildei,ztilder)*CirTerm
!
! If the Laurent expansion is also needed the poles of the
! I operator is calculated:
    if (present(CirLaurent)) then
      CirLaurent = 0d0
      mu2 = mur*mur
      P0qg = Pqg0(ztildei,ztilder)
      do jpart=1,ntilde
        if (abs(ptilde(jpart)%flv).gt.6) cycle
        if (ptilde(jpart)%flv.eq.0) then
          CirLaurent(-2) = CirLaurent(-2) &
                         + qcd_ca*P0qg*smeB
          CirLaurent(-1) = CirLaurent(-1) &
                         + qcd_ca*log(mu2/Q2)*P0qg*smeB &
                         + 0.5d0*qcd_beta0*P0qg*smeB
        else 
          CirLaurent(-2) = CirLaurent(-2) &
                         + qcd_cf*P0qg*smeB
          CirLaurent(-1) = CirLaurent(-1) &
                         + qcd_cf*log(mu2/Q2)*P0qg*smeB &
                         + 1.5d0*qcd_cf*P0qg*smeB
        end if
        do kpart=jpart+1,ntilde
          if (abs(ptilde(kpart)%flv).gt.6) cycle
          CirLaurent(-1) = CirLaurent(-1) &
                         + 2*log(yirtilde_arr(jpart,kpart)) &
                         * P0qg*Bij(jpart,kpart)
        end do
      end do
    end if
!
! g -> g + g
  elseif ((p(radi)%flv.eq.0).and.(p(radr)%flv.eq.0)) then 
!    print *,"g -> g + g"
    Bmunu_tmp = 0d0
    do jpart=1,ntilde
      if (abs(ptilde(jpart)%flv).gt.6) cycle
      if (ptilde(jpart)%flv.ne.0) then
        call CalcC1i0F(yitildeQ_arr(jpart),'q',Cj)
        Tj2 = qcd_cf
      else
        call CalcC1i0F(yitildeQ_arr(jpart),'g',Cj)
        Tj2 = qcd_ca
      end if
      Bmunu_tmp = Bmunu_tmp + (Cj + CS)*Tj2*Bmunu
      do kpart=jpart+1,ntilde
        if (abs(ptilde(kpart)%flv).gt.6) cycle
        Yjk = yirtilde_arr(jpart,kpart) &
            / (yitildeQ_arr(jpart)*yitildeQ_arr(kpart))
        call CalcS10ikFF(Yjk,Sjk)
        Bmunu_tmp(:,:) = Bmunu_tmp(:,:) &
                       + 2d0*(Sjk + CS)*Bmunuij(:,:,jpart,kpart)
      end do
    end do
    CirTerm = Pgg0(ztildei,ztilder,kttildei,Bmunu_tmp)
!
! If the Laurent expansion is also needed the poles of the
! I operator is calculated:
    if (present(CirLaurent)) then
      CirLaurent = 0d0
      mu2 = mur*mur
      P0gg = Pgg0(ztildei,ztilder,kttildei,Bmunu)
      Bmunu_tmp = 0d0
      do jpart=1,ntilde
        if (abs(ptilde(jpart)%flv).gt.6) cycle
        if (ptilde(jpart)%flv.eq.0) then
          CirLaurent(-2) = CirLaurent(-2) &
                         + qcd_ca*P0gg
          CirLaurent(-1) = CirLaurent(-1) &
                         + qcd_ca*log(mu2/Q2)*P0gg &
                         + 0.5d0*qcd_beta0*P0gg
        else 
          CirLaurent(-2) = CirLaurent(-2) &
                         + qcd_cf*P0gg
          CirLaurent(-1) = CirLaurent(-1) &
                         + qcd_cf*log(mu2/Q2)*P0gg &
                         + 1.5d0*qcd_cf*P0gg
        end if
        do kpart=jpart+1,ntilde
          if (abs(ptilde(kpart)%flv).gt.6) cycle
          Bmunu_tmp = Bmunu_tmp + 2*log(yirtilde_arr(jpart,kpart)) &
                    * Bmunuij(:,:,jpart,kpart)
        end do
      end do
      CirLaurent(-1) = CirLaurent(-1) &
                     + Pgg0(ztildei,ztilder,kttildei,Bmunu_tmp)
    end if
!
! g -> q + q~ || g -> q~ + q
  elseif ((p(radi)%flv.ne.0).and.(p(radr)%flv.ne.0).and. &
          (p(radi)%flv.eq.-p(radr)%flv)) then 
!    print *,"g -> q + q~"
    Bmunu_tmp = 0d0
    do jpart=1,ntilde
      if (abs(ptilde(jpart)%flv).gt.6) cycle
      if (ptilde(jpart)%flv.ne.0) then
        call CalcC1i0F(yitildeQ_arr(jpart),'q',Cj)
        Tj2 = qcd_cf
      else
        call CalcC1i0F(yitildeQ_arr(jpart),'g',Cj)
        Tj2 = qcd_ca
      end if
      Bmunu_tmp = Bmunu_tmp + (Cj + CS)*Tj2*Bmunu
      do kpart=jpart+1,ntilde
        if (abs(ptilde(kpart)%flv).gt.6) cycle
        Yjk = yirtilde_arr(jpart,kpart) &
            / (yitildeQ_arr(jpart)*yitildeQ_arr(kpart))
        call CalcS10ikFF(Yjk,Sjk)
        Bmunu_tmp(:,:) = Bmunu_tmp(:,:) &
                       + 2d0*(Sjk + CS)*Bmunuij(:,:,jpart,kpart)
      end do
    end do
    CirTerm = Pqq0(ztildei,ztilder,kttildei,Bmunu_tmp)
!
! If the Laurent expansion is also needed the poles of the
! I operator is calculated:
    if (present(CirLaurent)) then
      CirLaurent = 0d0
      mu2 = mur*mur
      P0qq = Pqq0(ztildei,ztilder,kttildei,Bmunu)
      Bmunu_tmp = 0d0
      do jpart=1,ntilde
        if (abs(ptilde(jpart)%flv).gt.6) cycle
        if (ptilde(jpart)%flv.eq.0) then
          CirLaurent(-2) = CirLaurent(-2) &
                         + qcd_ca*P0qq
          CirLaurent(-1) = CirLaurent(-1) &
                         + qcd_ca*log(mu2/Q2)*P0qq &
                         + 0.5d0*qcd_beta0*P0qq
        else 
          CirLaurent(-2) = CirLaurent(-2) &
                         + qcd_cf*P0qq
          CirLaurent(-1) = CirLaurent(-1) &
                         + qcd_cf*log(mu2/Q2)*P0qq &
                         + 1.5d0*qcd_cf*P0qq
        end if
        do kpart=jpart+1,ntilde
          if (abs(ptilde(kpart)%flv).gt.6) cycle
          Bmunu_tmp = Bmunu_tmp + 2*log(yirtilde_arr(jpart,kpart)) &
                    * Bmunuij(:,:,jpart,kpart)
        end do
      end do
      CirLaurent(-1) = CirLaurent(-1) &
                     + Pqq0(ztildei,ztilder,kttildei,Bmunu_tmp)
    end if
  else
    print *,"Erroneous config..."
    print *,"emit: ",emit
    print *,"radi: ",radi
    print *,"radr: ",radr
    call PrintParts(p)
    call PrintParts(ptilde)
    stop "CalcCir00IFF"
  end if
  CirTerm = 4d0/sir * CirTerm * (1d0 - alphair)**(2*d0 - 2*(m - 1))
!
  if (present(CirLaurent)) CirLaurent = CirLaurent * 4d0/sir
!
end subroutine CalcCir00IFF
!
subroutine CalcCirR00FF(p,ptilde,smeB,Bmunu,emit,radi,radr,CirTerm,CirLaurent)
use particles
use regions
use math
use QCDparams
use misc
use process
use scales
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
  type(particle) , dimension(:) , intent(in) :: ptilde
  real(kind(1d0)) , intent(in) :: smeB
  real(kind(1d0)) , dimension(0:3,0:3) , intent(in) :: Bmunu
  integer , intent(in) :: emit,radi,radr
  real(kind(1d0)) , intent(out) :: CirTerm
  real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: CirLaurent
!
  integer :: d0,m
  integer :: flv_i,flv_r,flv_ir
  real(kind(1d0)) :: Rir
!
  real(kind(1d0)) , dimension(-4:2) :: RirLaurent
!
  real(kind(1d0)) , external :: Pqg0,Pgg0,Pqq0
!
  CirTerm = 0d0
!
  d0 = sub_d0
  m  = size(ptilde) - 2
!
  flv_i  = p(radi)%flv
  flv_r  = p(radr)%flv
! The position of the emitter in the underlying Born always
! corresponds to the radiated leg lying in a lower position:
  flv_ir = ptilde(min(radi,radr))%flv
!
  if (.not.present(CirLaurent)) then
    call CalcRirF(flv_i,flv_r,flv_ir,ztildei,ztilder,yir,yirtildeQ, &
                  Rir)
  else
    call CalcRirF(flv_i,flv_r,flv_ir,ztildei,ztilder,yir,yirtildeQ, &
                  Rir,RirLaurent)
  end if
!
! q(~) -> q(~) + g
  if ((p(radi)%flv.ne.0).and.(p(radr)%flv.eq.0)) then
!    print *,"q(~) -> q(~) + g"
    CirTerm = 8d0*pi/sir*Rir*Pqg0(ztildei,ztilder)*smeB
!
    if (present(CirLaurent)) then
      CirLaurent = 8d0*pi/sir*RirLaurent*Pqg0(ztildei,ztilder)*smeB
    end if
!
! g -> g + g
  elseif ((p(radi)%flv.eq.0).and.(p(radr)%flv.eq.0)) then 
!    print *,"g -> g + g"
    CirTerm = 8d0*pi/sir*Rir*Pgg0(ztildei,ztilder,kttildei,Bmunu)
!
    if (present(CirLaurent)) then
      CirLaurent = 8d0*pi/sir*RirLaurent &
                 * Pgg0(ztildei,ztilder,kttildei,Bmunu)
    end if
!
! g -> q + q~ || g -> q~ + q
  elseif ((p(radi)%flv.ne.0).and.(p(radr)%flv.ne.0).and. &
          (p(radi)%flv.eq.-p(radr)%flv)) then 
!    print *,"g -> q + q~"
    CirTerm = 8d0*pi/sir*Rir*Pqq0(ztildei,ztilder,kttildei,Bmunu)
!
    if (present(CirLaurent)) then
      CirLaurent = 8d0*pi/sir*RirLaurent &
                 * Pqq0(ztildei,ztilder,kttildei,Bmunu)
    end if
!
  else
    print *,"Unknown branching in CalcCirR00FF..."
    stop
  end if
  CirTerm = CirTerm * (1d0 - alphair)**(2*d0 - 2*(m - 1))
!
end subroutine CalcCirR00FF
!
subroutine CalcSr00I(p,ptilde,Bij,Bijkl,radr,SrTerm,SrLaurent)
use particles
use regions
use math
use scales
use QCDparams
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
  type(particle) , dimension(:) , intent(in) :: ptilde
  real(kind(1d0)) , dimension(:,:) , intent(in) :: Bij
  real(kind(1d0)) , dimension(:,:,:,:) , intent(in) :: Bijkl
  integer , intent(in) :: radr
  real(kind(1d0)) , intent(out) :: SrTerm
  real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: SrLaurent
!
  integer :: i,k,r
  integer :: ipart,kpart,jpart,lpart,n,ntilde
  integer :: d0pr,m
  real(kind(1d0)) :: Sjl,CS,Cj
  real(kind(1d0)) :: Yjl
  real(kind(1d0)) :: Tj2,gammai
  real(kind(1d0)) :: mu2
!
  real(kind(1d0)) :: Sik
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
    subroutine CalcCirSr0F(CScont,CSLaurent)
    implicit none
!
      real(kind(1d0)) , intent(out) :: CScont
      real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: CSLaurent
!
    end subroutine CalcCirSr0F
  end interface
!
  Sik(i,k,r) = 2d0*sir_arr(i,k)/sir_arr(i,r)/sir_arr(k,r)
!
  SrTerm = 0d0
!
  d0pr = sub_d0pr
  m  = size(ptilde) - 2
!
  r = radr
!
! Calculation of the CS term which does not depend upon the
! kinematics:
  call CalcCirSr0F(CS)
!
! size of the original momentum array:
  n = size(p)
! size of the tilde momentum array:
  ntilde = size(ptilde)
!
  ipart = 0
! Sum over all momenta in the real-emission-like configuration:
  do i=1,n
! Omitting the radiated gluon:
    if (i.eq.r) cycle
! ipart is for book keeping the corresponding momentum
! in the underlying Born kinematics:
    ipart = ipart + 1
! Skipping non-QCD particles:
    if (abs(p(i)%flv).gt.6) cycle
! kpart is the same as ipart:
    kpart = ipart
! Loop over all possible ordered-pairs in the real-emission-like case:
    do k=i+1,n
! k cannot coincide with the radiated gluon:
      if (k.eq.r) cycle
      kpart = kpart + 1
! Only partons count:
      if (abs(p(k)%flv).gt.6) cycle
! loop over underlying Born momenta:
      do jpart=1,ntilde
! Only partons:
        if (abs(ptilde(jpart)%flv).gt.6) cycle
! (anti)Quarks:
        if (ptilde(jpart)%flv.ne.0) then
          call CalcC1i0F(yitildeQ_arr(jpart),'q',Cj)
          Tj2 = qcd_cf
! Gluons:
        else
          call CalcC1i0F(yitildeQ_arr(jpart),'g',Cj)
          Tj2 = qcd_ca
        end if
        SrTerm = SrTerm + Sik(i,k,r)*(Cj + CS)*Tj2*Bij(ipart,kpart)
! Ordered-pairs in the underlying Born kinematics:
        do lpart=jpart+1,ntilde
          if (abs(ptilde(lpart)%flv).gt.6) cycle
          Yjl = yirtilde_arr(jpart,lpart) &
              / (yitildeQ_arr(jpart)*yitildeQ_arr(lpart))
          call CalcS10ikFF(Yjl,Sjl)
          SrTerm = SrTerm + Sik(i,k,r)*(Sjl + CS) &
                 * Bijkl(ipart,kpart,jpart,lpart)
        end do
      end do
    end do
  end do
!
  SrTerm = -4d0 * SrTerm * (1d0 - yrQ)**(d0pr - (m - 1))
!
! For testing purposes it can be desirable to calculate the
! pole pieces too, this involves the pole part of the I operator,
! hence it has to be constructed:
  if (present(SrLaurent)) then
    SrLaurent = 0d0
    mu2 = mur*mur
    ipart = 0
    do i=1,n
      if (i.eq.r) cycle
      ipart = ipart + 1
      if (abs(p(i)%flv).gt.6) cycle
      kpart = ipart
      do k=i+1,n
        if (k.eq.r) cycle
        kpart = kpart + 1
        if (abs(p(k)%flv).gt.6) cycle
        do jpart=1,ntilde
          if (abs(ptilde(jpart)%flv).gt.6) cycle
          if (ptilde(jpart)%flv.ne.0) then
            Tj2 = qcd_cf
            gammai = 1.5d0*qcd_cf
          else
            Tj2 = qcd_ca
            gammai = 0.5d0*qcd_beta0
          end if
          SrLaurent(-2) = SrLaurent(-2) &
                        + Sik(i,k,r)*Tj2*Bij(ipart,kpart)
          SrLaurent(-1) = SrLaurent(-1) &
                        + Sik(i,k,r)*(Tj2*log(mu2/Q2) + gammai) &
                        * Bij(ipart,kpart)
          do lpart=jpart+1,ntilde
            if (abs(ptilde(lpart)%flv).gt.6) cycle
            SrLaurent(-1) = SrLaurent(-1) &
                          + Sik(i,k,r)*log(yirtilde_arr(jpart,lpart)) &
                          * Bijkl(ipart,kpart,jpart,lpart)
          end do
        end do
      end do
    end do
    SrLaurent = -4d0*SrLaurent
  end if
!
end subroutine CalcSr00I
!
subroutine CalcSrR00(p,ptilde,Bij,radr,SrTerm,SrLaurent)
use particles
use regions
use math
use QCDparams
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
  type(particle) , dimension(:) , intent(in) :: ptilde
  real(kind(1d0)) , dimension(:,:) , intent(in) :: Bij
  integer , intent(in) :: radr
  real(kind(1d0)) , intent(out) :: SrTerm
  real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: SrLaurent
!
  integer :: i,k,r
  integer :: ipart,kpart
  integer :: d0pr,m
  integer :: n
  real(kind(1d0)) :: Rikr
  real(kind(1d0)) , dimension(-4:2) :: RikrLaurent
!
  real(kind(1d0)) :: Sik
!
  Sik(i,k,r) = 2d0*sir_arr(i,k)/sir_arr(i,r)/sir_arr(k,r)
!
  SrTerm = 0d0
  if (present(SrLaurent)) SrLaurent = 0d0
!
  d0pr = sub_d0pr
  m  = size(ptilde) - 2
!
  n = size(p)
!
  r = radr
!
  ipart = 0
  do i=1,n-1
    if (i.eq.r) cycle
    ipart = ipart + 1
    if (abs(p(i)%flv).gt.6) cycle
    kpart = ipart
    do k=i+1,n
      if (k.eq.r) cycle
      kpart = kpart + 1
      if (abs(p(k)%flv).gt.6) cycle
      if (.not.present(SrLaurent)) then
        call CalcRikrF(yir_arr(i,k),yir_arr(i,r),yir_arr(k,r), &
                       yiQ_arr(i),yiQ_arr(k),yiQ_arr(r), &
                       Rikr)
        SrTerm = SrTerm &
               + Sik(i,k,r) * Rikr * Bij(ipart,kpart)
      else
        call CalcRikrF(yir_arr(i,k),yir_arr(i,r),yir_arr(k,r), &
                       yiQ_arr(i),yiQ_arr(k),yiQ_arr(r), &
                       Rikr,RikrLaurent)
        SrTerm = SrTerm &
               + Sik(i,k,r) * Rikr * Bij(ipart,kpart)
        SrLaurent = SrLaurent &
                  + Sik(i,k,r) * RikrLaurent * Bij(ipart,kpart)
      end if
    end do
  end do
!
  SrTerm = -8d0*pi * SrTerm * (1d0 - yrQ)**(d0pr - (m - 1))
!
  if (present(SrLaurent)) SrLaurent = -8d0*pi*SrLaurent
!
end subroutine CalcSrR00
!
subroutine CalcCirFFSr00I(p,ptilde,smeB,Bij,emit,radi,radr, &
                          CirSrTerm,CirSrLaurent)
use particles
use regions
use math
use QCDparams
use scales
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
  type(particle) , dimension(:) , intent(in) :: ptilde
  real(kind(1d0)) , intent(in) :: smeB
  real(kind(1d0)) , dimension(:,:) , intent(in) :: Bij
  integer , intent(in) :: emit,radi,radr
  real(kind(1d0)) , intent(out) :: CirSrTerm
  real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: CirSrLaurent
!
  integer :: d0pr,m
  integer :: jpart,lpart,ntilde
  real(kind(1d0)) :: zi,zr,Ti2,Tj2,Yjl
  real(kind(1d0)) :: Sjl,CS,Cj
  real(kind(1d0)) :: mu2,gammai
!
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
  CirSrTerm = 0d0
!
  d0pr = sub_d0pr
  m  = size(ptilde) - 2
!
! We have to calculate the momentum fractions:
  yirQ = yiQ_arr(radi) + yiQ_arr(radr)
  zi = yiQ_arr(radi)/yirQ
  zr = 1d0 - zi
  sir     = sir_arr(radi,radr)
!
  if (abs(p(emit)%flv).ne.0) Ti2 = qcd_cf
  if (abs(p(emit)%flv).eq.0) Ti2 = qcd_ca
!
! Calculation of the CS term which does not depend upon the
! kinematics:
  call CalcCirSr0F(CS)
!
! To construct the finite part of the subtraction term the tilde
! momenta have to be used:
  ntilde = size(ptilde)
  do jpart=1,ntilde
    if (abs(ptilde(jpart)%flv).gt.qcd_nf) cycle
! Calculation of the Cj term:
    if (ptilde(jpart)%flv.ne.0) then
      call CalcC1i0F(yitildeQ_arr(jpart),'q',Cj)
      Tj2 = qcd_cf
    else
      call CalcC1i0F(yitildeQ_arr(jpart),'g',Cj)
      Tj2 = qcd_ca
    end if
    CirSrTerm = CirSrTerm + (Cj + CS)*Tj2*smeB
    do lpart=jpart+1,ntilde
      if (abs(ptilde(lpart)%flv).gt.qcd_nf) cycle
      Yjl = yirtilde_arr(jpart,lpart)/(yitildeQ_arr(jpart)*yitildeQ_arr(lpart))
      call CalcS10ikFF(Yjl,Sjl)
! Factor of two is introduced since only half of the sum is
! evaluated:
      CirSrTerm = CirSrTerm + 2d0*(Sjl + CS)*Bij(jpart,lpart)
    end do
  end do
!
  CirSrTerm = 4d0 / sir*2d0*zi/zr * Ti2 * CirSrTerm &
            * (1d0 - yrQ)**(d0pr - (m - 1))
!
  if (present(CirSrLaurent)) then
    CirSrLaurent = 0d0 
    mu2 = mur*mur
    do jpart=1,ntilde
      if (abs(ptilde(jpart)%flv).gt.6) cycle
      if (ptilde(jpart)%flv.ne.0) then
        Tj2 = qcd_cf
        gammai = 1.5d0*qcd_cf
      else
        Tj2 = qcd_ca
        gammai = 0.5d0*qcd_beta0
      end if
      CirSrLaurent(-2) = CirSrLaurent(-2) + Tj2*smeB
      CirSrLaurent(-1) = CirSrLaurent(-1) &
                       + (Tj2*log(mu2/Q2) + gammai)*smeB
      do lpart=jpart+1,ntilde
        if (abs(ptilde(lpart)%flv).gt.6) cycle
        CirSrLaurent(-1) = CirSrLaurent(-1) &
                         + 2d0*log(yirtilde_arr(jpart,lpart)) &
                         * Bij(jpart,lpart)
      end do
    end do
    CirSrLaurent = 8d0/sir*zi/zr*Ti2*cirSrLaurent
  end if
!
end subroutine CalcCirFFSr00I
!
subroutine CalcCirFFSrR00(p,ptilde,smeB,emit,radi,radr,CirSrTerm,CirSrLaurent)
use particles
use regions
use math
use QCDparams
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
  type(particle) , dimension(:) , intent(in) :: ptilde
  real(kind(1d0)) , intent(in) :: smeB
  integer , intent(in) :: emit,radi,radr
  real(kind(1d0)) , intent(out) :: CirSrTerm 
  real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: CirSrLaurent 
!
  integer :: d0pr,m
  real(kind(1d0)) :: zi,zr,T2smeB
  real(kind(1d0)) :: Y,Cg,SSir
  real(kind(1d0)) , dimension(:) , allocatable :: CgLaurent,SirLaurent
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
  CirSrTerm = 0d0
!
  d0pr = sub_d0pr
  m  = size(ptilde) - 2
!
! We have to calculate the momentum fractions:
  yirQ = yiQ_arr(radi) + yiQ_arr(radr)
  zi = yiQ_arr(radi)/yirQ
  zr = 1d0 - zi
  sir = sir_arr(radi,radr)
!
  if (abs(p(emit)%flv).ne.0) T2smeB = qcd_cf*smeB
  if (abs(p(emit)%flv).eq.0) T2smeB = qcd_ca*smeB
!
  Y = sir/(Q2*zi*yirQ*zr*yirQ)
!
  if (.not.present(CirSrLaurent)) then
    call CalcC1i0F(zr*yirQ,'g',Cg)
    call CalcS10ikFF(Y,SSir)
  else
    allocate(CgLaurent(-4:2),SirLaurent(-4:2))
    call CalcC1i0F(zr*yirQ,'g',Cg,CgLaurent)
    call CalcS10ikFF(Y,SSir,SirLaurent)
    CirSrLaurent = 8d0*pi / sir*2d0*zi/zr * T2smeB &
                 * qcd_ca/(2d0*pi)*(CgLaurent - SirLaurent)
    deallocate(CgLaurent,SirLaurent)
  end if
!
  CirSrTerm = 8d0*pi / sir*2d0*zi/zr * T2smeB &
            * qcd_ca/(2d0*pi) * (Cg - SSir) &
            * (1d0 - yrQ)**(d0pr - (m - 1))
!
!
end subroutine CalcCirFFSrR00
!
subroutine CalcCirsFF(p,ptilde,smeB,Bmunu,emit,radi,radr,rads,CirsTerm)
use particles
use regions
use math
use utils
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p,ptilde
  real(kind(1d0)) , intent(in) :: smeB
  real(kind(1d0)) , dimension(0:3,0:3) , intent(in) :: Bmunu
  integer , intent(in) :: emit,radi,radr,rads
  real(kind(1d0)) , intent(out) :: CirsTerm
!
  real(kind(1d0)) :: s_ir,s_is,s_rs
!
  real(kind(1d0)) , external :: Prbrq0,Pqbqq0,Pggq0,Pgqqb0,Pggg0
!
!
  CirsTerm = 0d0
!
  s_ir = sir_arr(radi,radr)
  s_is = sir_arr(radi,rads)
  s_rs = sir_arr(radr,rads)
!
! The orderings are always:
! q(~) -> q(~) + r + r~
! q(~) -> q(~) + q + q~
! q(~) -> q(~) + g + g
!
! Preselection according to the flavor of the emitter:
! (anti)quark:
  if (ptilde(emit)%flv.ne.0) then
! q(~) -> q(~) + r + r~
    if ((p(radr)%flv.ne.0).and. &
        (abs(p(radi)%flv).ne.abs(p(radr)%flv)).and. &
        (p(radr)%flv.eq.-p(rads)%flv)) then
!      print *,"q(~) -> q(~) + r + r~"
      CirsTerm = Prbrq0(sirs,s_rs,s_is,s_ir,ztildes,ztilder,ztildei)*smeB
! q -> q + q + q~
    elseif ((p(radr)%flv.ne.0).and. &
        (p(radi)%flv.eq.p(radr)%flv).and. &
        (p(radr)%flv.eq.-p(rads)%flv)) then
!      print *,"q -> q + q + q~"
      CirsTerm = Pqbqq0(sirs,s_rs,s_is,s_ir,ztildes,ztilder,ztildei)*smeB
! q~ -> q~ + q + q~
    elseif ((p(radr)%flv.ne.0).and. &
        (p(radi)%flv.eq.p(rads)%flv).and. &
        (p(radr)%flv.eq.-p(rads)%flv)) then
!      print *,"q~ -> q~ + q + q~"
! The splitting kernel for q~ -> q~ + q + q~ is obtained from the one for
! q -> q + q + q~ with a r <-> s change:
      CirsTerm = Pqbqq0(sirs,s_rs,s_ir,s_is,ztilder,ztildes,ztildei)*smeB
! q(~) -> q(~) + g + g
    elseif ((p(radr)%flv.eq.0).and. &
        (p(radr)%flv.eq.p(rads)%flv)) then
!      print *,"q(~) -> q(~) + g + g"
      CirsTerm = Pggq0(sirs,s_rs,s_is,s_ir,ztildes,ztilder,ztildei)*smeB
    else
      print *,"Error in particle ordering in CalcCirsFF..."
      print *,ConvertFromPDG(ptilde(emit)%flv)," -> ", &
              ConvertFromPDG(p(radi)%flv)," , ", &
              ConvertFromPDG(p(radr)%flv)," , ", &
              ConvertFromPDG(p(rads)%flv)
      stop "CalcCirsFF"
    end if
! gluon:
  else
! g -> g + q + q~
    if ((p(radr)%flv.ne.0).and. &
        (p(radr)%flv.eq.-p(rads)%flv)) then
!      print *,"g -> g + q + q~"
      CirsTerm = Pgqqb0(sirs,s_ir,s_is,s_rs,ztildei,ztilder,ztildes, &
                        kttildei,kttilder,kttildes,smeB,Bmunu)
! g -> g + g + g
    elseif ((p(radr)%flv.eq.0).and. &
        (p(radr)%flv.eq.p(rads)%flv)) then
!      print *,"g -> g + g + g"
!
      CirsTerm = Pggg0(sirs,s_ir,s_is,s_rs,ztildei,ztilder,ztildes, &
                       kttildei,kttilder,kttildes,smeB,Bmunu)
    else
      print *,"Error in particle ordering in CalcCirsFF..."
      print *,ConvertFromPDG(ptilde(emit)%flv)," -> ", &
              ConvertFromPDG(p(radi)%flv)," , ", &
              ConvertFromPDG(p(radr)%flv)," , ", &
              ConvertFromPDG(p(rads)%flv)
      stop "CalcCirsFF"
    end if
  end if
!
  CirsTerm = (8d0*pi)**2/sirs**2*CirsTerm
!
end subroutine CalcCirsFF
!
subroutine CalcCirjsFF(p,ptilde,smeB,Bmunu1,Bmunu2,Balbemunu, &
                       emiti,radi,radr,emitj,radj,rads, &
                       CirjsTerm)
use particles
use regions
use math
use utils
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p,ptilde
  real(kind(1d0)) , intent(in) :: smeB
  real(kind(1d0)) , dimension(0:3,0:3) , intent(in) :: Bmunu1,Bmunu2
  real(kind(1d0)) , dimension(0:3,0:3,0:3,0:3) , intent(in) :: Balbemunu
  integer , intent(in) :: emiti,emitj,radi,radr,radj,rads
  real(kind(1d0)) , intent(out) :: CirjsTerm
!
!
  real(kind(1d0)) , external :: Pqg0,Pgg0,Pqq0
!
  CirjsTerm = 0d0
!
! Both of the emitters are quarks:
  if ((ptilde(emiti)%flv.ne.0).and. &
      (ptilde(emitj)%flv.ne.0)) then
!    print *,"q(~) -> q(~) g ; r(~) -> r(~) g"
    CirjsTerm = Pqg0(ztildei,ztilder) &
              * Pqg0(ztildej,ztildes)*smeB
! The first emitter is a quark, the second is a gluon splitting 
! into a quark-pair:
  elseif ((ptilde(emiti)%flv.ne.0).and. &
          (ptilde(emitj)%flv.eq.0).and. &
          (p(radj)%flv.ne.0)) then
!    print *,"q(~) -> q(~) g ; g -> r r(~)"
    CirjsTerm = Pqg0(ztildei,ztilder) &
              * Pqq0(ztildej,ztildes,kttildej,Bmunu2)
! The first emitter is a quark, the second is a gluon splitting 
! into a gluon-pair:
  elseif ((ptilde(emiti)%flv.ne.0).and. &
          (ptilde(emitj)%flv.eq.0).and. &
          (p(radj)%flv.eq.0)) then
!    print *,"q(~) -> q(~) g ; g -> g g"
    CirjsTerm = Pqg0(ztildei,ztilder) &
              * Pgg0(ztildej,ztildes,kttildej,Bmunu2)
! The first emitter is a gluon splitting into a quark-pair,
! the second is a quark:
  elseif ((ptilde(emiti)%flv.eq.0).and. &
          (ptilde(emitj)%flv.ne.0).and. &
          (p(radi)%flv.ne.0)) then
!    print *,"g -> r r(~) ; q(~) -> q(~) g"
    CirjsTerm = Pqq0(ztildei,ztilder,kttildei,Bmunu1) &
              * Pqg0(ztildej,ztildes)
! The first emitter is a gluon splitting into a gluon-pair, 
! the second is a quark:
  elseif ((ptilde(emiti)%flv.eq.0).and. &
          (ptilde(emitj)%flv.ne.0).and. &
          (p(radi)%flv.eq.0)) then
!    print *,"g -> g g ; q(~) -> q(~) g"
    CirjsTerm = Pgg0(ztildei,ztilder,kttildei,Bmunu1) &
              * Pqg0(ztildej,ztildes)
  else
    print *,"Error in particle ordering in CalcCirjsFF..."
    print *,ConvertFromPDG(ptilde(emiti)%flv)," -> ", &
            ConvertFromPDG(p(radi)%flv)," , ", &
            ConvertFromPDG(p(radr)%flv)
    print *,ConvertFromPDG(ptilde(emitj)%flv)," -> ", &
            ConvertFromPDG(p(radj)%flv)," , ", &
            ConvertFromPDG(p(rads)%flv)
    stop "CalcCirjsFF"
  end if
!
  CirjsTerm = (8d0*pi)**2/(sir*sjs)*CirjsTerm
!
end subroutine CalcCirjsFF
!
subroutine CalcCSirsFF(p,ptilde,smeB,Bij,Bmunuij, &
                       emitir,radi,radr,rads,CSirsTerm)
use particles
use regions
use math
use utils
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p,ptilde
  real(kind(1d0)) , intent(in) :: smeB
  real(kind(1d0)) , dimension(:,:) , intent(in) :: Bij
  real(kind(1d0)) , dimension(0:,0:,:,:) , intent(in) :: Bmunuij
  integer , intent(in) :: emitir,radi,radr,rads
  real(kind(1d0)) , intent(out) :: CSirsTerm
!
  integer :: j,k,s,n,emit,rad
  integer :: jpart,kpart
  real(kind(1d0)) :: Pqg,Pgg,Pqq
  real(kind(1d0)) , dimension(0:3,0:3) :: Bmunu
!
  real(kind(1d0)) :: Sik
  real(kind(1d0)) , external :: Pqg0,Pgg0,Pqq0
!
  Sik(j,k,s) = 2d0*sir_arr(j,k)/(sir_arr(j,s)*sir_arr(k,s))
!
  CSirsTerm = 0d0
!
  s = rads
  n = size(p)
! The position of the radiated parton originating from the 
! collinear pair:
  rad  = max(radi,radr)
  emit = min(radi,radr)
!
! Determining the emitter for the collinear-pair:
! q(~) -> q(~) + g
  if (ptilde(emitir)%flv.ne.0) then
!    print *,"q(~) -> q(~) + g"
    Pqg = Pqg0(ztildei,ztilder)
    jpart = 0
    do j=1,n
      if ((j.eq.s).or.(j.eq.rad)) cycle
      jpart = jpart + 1
      if (j.eq.emit) cycle
      if (abs(p(j)%flv).gt.6) cycle
      kpart = jpart
      do k=j+1,n
        if ((k.eq.s).or.(k.eq.rad)) cycle
        kpart = kpart + 1
        if (k.eq.emit) cycle
        if (abs(p(k)%flv).gt.6) cycle
        CSirsTerm = CSirsTerm + Sik(j,k,s)*Pqg*Bij(jpart,kpart)
      end do
    end do
! Contribution is also coming from the (ir) leg, but for this
! particular one the form of \mathcal{S} changes,
! for the precise definition take a look at Eq. (6.26) of 
! arXiv:hep-ph/0609042.
    kpart = 0
    do k=1,n
      if ((k.eq.s).or.(k.eq.rad)) cycle
      kpart = kpart + 1
      if (k.eq.emit) cycle
      if (abs(p(k)%flv).gt.6) cycle
      CSirsTerm = CSirsTerm &
                + 2*(sir_arr(radi,k) + sir_arr(radr,k)) &
                /((sir_arr(radi,rads) + sir_arr(radr,rads))*sir_arr(k,rads)) &
                * Pqg*Bij(emitir,kpart)
    end do
!
! g -> g + g || g -> q + q~
  else
!    print *,"g -> g + g || g -> q + q~"
    Bmunu = 0d0
    jpart = 0
    do j=1,n
      if ((j.eq.s).or.(j.eq.rad)) cycle
      jpart = jpart + 1
      if (j.eq.emit) cycle
      if (abs(p(j)%flv).gt.6) cycle
      kpart = jpart
      do k=j+1,n
        if ((k.eq.s).or.(k.eq.rad)) cycle
        kpart = kpart + 1
        if (k.eq.emit) cycle
        if (abs(p(k)%flv).gt.6) cycle
        Bmunu = Bmunu + Sik(j,k,s)*Bmunuij(:,:,jpart,kpart)
      end do
    end do
    kpart = 0
    do k=1,n
      if ((k.eq.s).or.(k.eq.rad)) cycle
      kpart = kpart + 1
      if (k.eq.emit) cycle
      if (abs(p(k)%flv).gt.6) cycle
      Bmunu = Bmunu &
            + 2*(sir_arr(radi,k) + sir_arr(radr,k)) &
            /((sir_arr(radi,rads) + sir_arr(radr,rads))*sir_arr(k,rads)) &
            * Bmunuij(:,:,emitir,kpart)
    end do
! g -> q + q~
    if (p(radi)%flv.ne.0) then
      Pqq = Pqq0(ztildei,ztilder,kttildei,Bmunu)
      CSirsTerm = Pqq
!
! g -> g + g
    else
      Pgg = Pgg0(ztildei,ztilder,kttildei,Bmunu)
      CSirsTerm = Pgg
    end if
!
  end if
!
  CSirsTerm = -(8d0*pi)**2/sir*CSirsTerm
!
end subroutine CalcCSirsFF
!
subroutine CalcCirsCSirsFF(p,ptilde,smeB,Bmunu, &
                           emit,radi,radr,rads,CirsCSirsTerm)
use QCDparams
use particles
use regions
use math
use utils
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p,ptilde
  real(kind(1d0)) , intent(in) :: smeB
  real(kind(1d0)) , dimension(0:3,0:3) , intent(in) :: Bmunu
  integer , intent(in) :: emit,radi,radr,rads
  real(kind(1d0)) , intent(out) :: CirsCSirsTerm
!
  real(kind(1d0)) :: s_ir_s
!
  real(kind(1d0)) , external :: Pqg0,Pgg0,Pqq0
!
  CirsCSirsTerm = 0d0
!
! Additional invariants have to be calculated:
  s_ir_s = sir_arr(radi,rads) + sir_arr(radr,rads)
  yiQ  = yiQ_arr(radi)
  yrQ  = yiQ_arr(radr)
  ysQ  = yiQ_arr(rads)
  yirsQ = yiQ + yrQ + ysQ
  ztildes = ysQ/yirsQ
!
! The Collinear-pair can be of several different types:
! q(~) -> q(~) + g
  if (ptilde(emit)%flv.ne.0) then
!    print *,"q(~) -> q(~) + g"
    CirsCSirsTerm = 2d0/(s_ir_s*ztildes)*(1 - ztildes) * qcd_cf &
                  * Pqg0(ztildei,ztilder)*smeB
!
! g -> g -> q + q~
  elseif (p(radi)%flv.ne.0) then
!    print *,"g -> g -> q + q~"
    CirsCSirsTerm = 2d0/(s_ir_s*ztildes)*(1 - ztildes) * qcd_ca &
                  * Pqq0(ztildei,ztilder,kttildei,Bmunu)
!
! g -> g + g
  elseif (p(radi)%flv.eq.0) then
!    print *,"g -> g + g"
    CirsCSirsTerm = 2d0/(s_ir_s*ztildes)*(1 - ztildes) * qcd_ca &
                  * Pgg0(ztildei,ztilder,kttildei,Bmunu)
  end if
!
  CirsCSirsTerm = (8d0*pi)**2/sir*CirsCSirsTerm
!
end subroutine CalcCirsCSirsFF
!
subroutine CalcCirjsCSirsFF(p,ptilde,smeB,Bmunu, &
                            emitir,radi,radr,    &
                            emitjs,radj,rads,    &
                            CirjsCSirsTerm)
use QCDparams
use particles
use regions
use math
use utils
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p,ptilde
  real(kind(1d0)) , intent(in) :: smeB
  real(kind(1d0)) , dimension(0:3,0:3) , intent(in) :: Bmunu
  integer , intent(in) :: emitir,emitjs,radi,radr,radj,rads
  real(kind(1d0)) , intent(out) :: CirjsCSirsTerm
!
  real(kind(1d0)) :: Tj2
!
  real(kind(1d0)) , external :: Pqg0,Pgg0,Pqq0
!
  CirjsCSirsTerm = 0d0
!
! Additional invariants have to be calculated:
  sjs  = sir_arr(radj,rads)
  yjQ  = yiQ_arr(radj)
  ysQ  = yiQ_arr(rads)
  yjsQ  = yjQ + ysQ
  ztildej = yjQ/yjsQ
  ztildes = 1 - ztildej
!
! Color charge for the parton j:
  if (p(radj)%flv.eq.0) then
    Tj2 = qcd_ca
  else
    Tj2 = qcd_cf
  end if
!
! The (ir) -> i + r Collinear-pair can be of several different types:
! q(~) -> q(~) + g
  if (ptilde(emitir)%flv.ne.0) then
!    print *,"q(~) -> q(~) + g"
    CirjsCSirsTerm = 2d0/(sjs*ztildes)*(1 - ztildes) * Tj2 &
                   * Pqg0(ztildei,ztilder)*smeB
!
! g -> g -> q + q~
  elseif (p(radi)%flv.ne.0) then
!    print *,"g -> g -> q + q~"
    CirjsCSirsTerm = 2d0/(sjs*ztildes)*(1 - ztildes) * Tj2 &
                   * Pqq0(ztildei,ztilder,kttildei,Bmunu)
!
! g -> g + g
  elseif (p(radi)%flv.eq.0) then
!    print *,"g -> g + g"
    CirjsCSirsTerm = 2d0/(sjs*ztildes)*(1 - ztildes) * Tj2 &
                   * Pgg0(ztildei,ztilder,kttildei,Bmunu)
  end if
!
  CirjsCSirsTerm = (8d0*pi)**2/sir*CirjsCSirsTerm
!
end subroutine CalcCirjsCSirsFF
!
subroutine CalcSrs(p,ptilde,Bij,Bijkl,radr,rads,SrsTerm)
use particles
use regions
use math
use QCDparams
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
  type(particle) , dimension(:) , intent(in) :: ptilde
  real(kind(1d0)) , dimension(:,:) , intent(in) :: Bij
  real(kind(1d0)) , dimension(:,:,:,:) , intent(in) :: Bijkl
  integer , intent(in) :: radr,rads
  real(kind(1d0)) , intent(out) :: SrsTerm
!
  integer :: i,j,k,l,r,s,n
  integer :: ipart,jpart,kpart,lpart
  real(kind(1d0)) :: srs
!
  real(kind(1d0)) :: Sikr,SikrsSO,Sik_rs,Sikrs,Sqq
!
  Sikr(i,k,r) = 2*sir_arr(i,k)/(sir_arr(i,r)*sir_arr(k,r))
! As defined in Eq. (6.39) of arXiv:hep-ph/0609042:
  SikrsSO(i,k,r,s) = Sikr(i,k,s) &
                   * (Sikr(i,s,r) + Sikr(k,s,r) - Sikr(i,k,r))
! As defined in Eq. (6.40) of arXiv:hep-ph/0609042:
  Sik_rs(i,k,r,s) = 2*sir_arr(i,k) &
                  / ((sir_arr(i,r) + sir_arr(i,s)) &
                  *  (sir_arr(k,r) + sir_arr(k,s)))
! As defined in Eq. (6.38) of arXiv:hep-ph/0609042:
  Sikrs(i,k,r,s) = SikrsSO(i,k,r,s) &
                 + 4*(sir_arr(i,r)*sir_arr(k,s) &
                 +    sir_arr(k,r)*sir_arr(i,s)) &
                 / ((sir_arr(i,r) + sir_arr(i,s)) &
                 *  (sir_arr(k,r) + sir_arr(k,s))) &
                 * (1d0/srs**2 - 0.125d0*SikrsSO(i,k,r,s)) &
                 - 4*Sik_rs(i,k,r,s)/srs
! This statement function corresponds to the terms in the parentheses
! of Eq. (6.41) of arXiv:hep-ph/0609042
  Sqq(i,k,r,s) = (sir_arr(i,r)*sir_arr(k,s) &
               +  sir_arr(k,r)*sir_arr(i,s) &
               -  sir_arr(i,k)*sir_arr(r,s)) &
               / ((sir_arr(i,r) + sir_arr(i,s)) &
               *  (sir_arr(k,r) + sir_arr(k,s))) &
               - 2*sir_arr(i,r)*sir_arr(i,s) &
               / (sir_arr(i,r) + sir_arr(i,s))**2
!
  SrsTerm = 0d0
!
  srs = yrs*Q2
!
  n = size(p)
!
  r = radr
  s = rads
!
  if ((p(radr)%flv.ne.0).and.(p(rads)%flv.ne.0)) then
    ipart = 0
    do i=1,n-1
      if ((i.eq.r).or.(i.eq.s)) cycle
      ipart = ipart + 1
      if (abs(p(i)%flv).gt.6) cycle
      kpart = ipart
      do k=i+1,n
        if ((k.eq.r).or.(k.eq.s)) cycle
        kpart = kpart + 1
        if (abs(p(k)%flv).gt.6) cycle
        SrsTerm = SrsTerm &
                + Sqq(i,k,r,s)*Bij(ipart,kpart) &
                + Sqq(k,i,r,s)*Bij(kpart,ipart)
      end do
    end do
    SrsTerm = qcd_tr/srs**2*SrsTerm
  elseif ((p(radr)%flv.eq.0).and.(p(rads)%flv.eq.0)) then
    ipart = 0
    do i=1,n
      if ((i.eq.r).or.(i.eq.s)) cycle
      ipart = ipart + 1
      if (abs(p(i)%flv).gt.6) cycle
! Sikrs has a part which does not vanish if i = k, hence the following
! term gives a non-vanishing contribution:
      SrsTerm = SrsTerm &
              - 2*qcd_ca*sir_arr(i,r)*sir_arr(i,s) &
              / (sir_arr(i,r) + sir_arr(i,s))**2 &
              / srs**2 &
              * Bij(ipart,ipart)
      kpart = ipart
      do k=i+1,n
        if ((k.eq.r).or.(k.eq.s)) cycle
        kpart = kpart + 1
        if (abs(p(k)%flv).gt.6) cycle
        SrsTerm = SrsTerm - 0.5d0*qcd_ca*Sikrs(i,k,r,s)*Bij(ipart,kpart)
! Contribution coming from the case when i = j and k = l:
        SrsTerm = SrsTerm &
                + 0.5d0*Sikr(i,k,r)*Sikr(i,k,s) &
                * Bijkl(ipart,kpart,ipart,kpart)
        jpart = ipart - 1
        do j=i,n-1
          if ((j.eq.r).or.(j.eq.s)) cycle
          jpart = jpart + 1
          if (abs(p(j)%flv).gt.6) cycle
          lpart = kpart - 1
          do l=k,n
            if ((l.eq.r).or.(l.eq.s)) cycle
            lpart = lpart + 1
            if (j.ge.l) cycle
            if ((i.eq.j).and.(k.eq.l)) cycle
            if (abs(p(l)%flv).gt.6) cycle
            SrsTerm = SrsTerm &
                    + 0.5d0*(Sikr(i,k,r)*Sikr(j,l,s) &
                    + Sikr(i,k,s)*Sikr(j,l,r)) &
                    * Bijkl(ipart,kpart,jpart,lpart)
          end do
        end do
      end do
    end do
  else
    print *,"Error in CalcSrs..."
    print *,"radr,rads: ",radr,rads
    stop "CalcSrs"
  end if
!
  SrsTerm = (8d0*pi)**2 * SrsTerm
!
end subroutine CalcSrs
!
subroutine CalcCirsSrs(p,ptilde,smeB,emitirs,radi,radr,rads,CirsSrsTerm)
use particles
use regions
use math
use QCDparams
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
  type(particle) , dimension(:) , intent(in) :: ptilde
  real(kind(1d0)) , intent(in) :: smeB
  integer , intent(in) :: emitirs,radi,radr,rads
  real(kind(1d0)) , intent(out) :: CirsSrsTerm
!
  integer :: i,r,s
  real(kind(1d0)) :: sis,srs,s_i_rs
  real(kind(1d0)) :: Ti2
!
!
  CirsSrsTerm = 0d0
!
  i = radi
  r = radr
  s = rads
!
  srs = yrs*Q2
  sir = sir_arr(i,r)
  sis = sir_arr(i,s)
  s_i_rs = sir + sis
  yirsQ = yiQ_arr(radi) + yiQ_arr(radr) + yiQ_arr(rads)
  ztildei = yiQ_arr(radi)/yirsQ
  ztilder = yiQ_arr(radr)/yirsQ
  ztildes = yiQ_arr(rads)/yirsQ
!
  if (p(radi)%flv.ne.0) then
    Ti2 = qcd_cf
  else
    Ti2 = qcd_ca
  end if
!
  if ((p(radr)%flv.ne.0).and.(p(rads)%flv.ne.0)) then
    CirsSrsTerm = 2*Ti2*qcd_tr/(s_i_rs*srs) &
                * (ztildei/(ztilder + ztildes) &
                -  (sir*ztildes - sis*ztilder)**2 &
                /  (s_i_rs*srs*(ztilder + ztildes)**2)) &
                * smeB
  elseif ((p(radr)%flv.eq.0).and.(p(rads)%flv.eq.0)) then
    CirsSrsTerm = 4*Ti2*ztildei**2/(sir*sis*ztilder*ztildes) &
                + qcd_ca &
                * ( &
                   (sir*ztildes - sis*ztilder)**2 &
                /  ((s_i_rs*srs)**2*(ztilder + ztildes)**2) &
                -  ztildei/(s_i_rs*srs)*(4d0/(ztilder + ztildes) &
                -   1d0/ztilder) &
                -  2*ztildei**2/(s_i_rs*sir*ztilder*(ztilder + ztildes)) &
                -  ztildei**2/(s_i_rs*sis*ztilder*(ztilder + ztildes)) &
                +  ztildei/(sir*srs) &
                *   (1d0/ztildes + 1d0/(ztilder + ztildes)) &
! r <-> s:
                +  (sis*ztilder - sir*ztildes)**2 &
                /  ((s_i_rs*srs)**2*(ztildes + ztilder)**2) &
                -  ztildei/(s_i_rs*srs)*(4d0/(ztildes + ztilder) &
                -   1d0/ztildes) &
                -  2*ztildei**2/(s_i_rs*sis*ztildes*(ztildes + ztilder)) &
                -  ztildei**2/(s_i_rs*sir*ztildes*(ztildes + ztilder)) &
                +  ztildei/(sis*srs) &
                *   (1d0/ztilder + 1d0/(ztildes + ztilder)) &
                  )
    CirsSrsTerm = CirsSrsTerm * Ti2*smeB
  else
    print *,"Error in CalcCirsSrs..."
    print *,"emitirs,radi,radr,rads: ",emitirs,radi,radr,rads
    stop "CalcCirsSrs"
  end if
!
  CirsSrsTerm = (8d0*pi)**2 * CirsSrsTerm
!
end subroutine CalcCirsSrs
!
subroutine CalcCSirsSrs(p,ptilde,Bij,emitir,radi,radr,rads,CSirsSrsTerm)
use particles
use regions
use math
use QCDparams
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
  type(particle) , dimension(:) , intent(in) :: ptilde
  real(kind(1d0)) , dimension(:,:) , intent(in) :: Bij
  integer , intent(in) :: emitir,radi,radr,rads
  real(kind(1d0)) , intent(out) :: CSirsSrsTerm
!
  integer :: i,k,r,j,l,n
  integer :: jpart,lpart
  real(kind(1d0)) :: Ti2
!
  real(kind(1d0)) :: Sikr,S_ij_k_r
!
  Sikr(i,k,r) = 2*sir_arr(i,k)/(sir_arr(i,r)*sir_arr(k,r))
! Modified version of eikonal factor as of Eq. (6.26) in 
! arXiv:hep-ph/0609042:
  S_ij_k_r(i,j,k,r) = 2*(sir_arr(i,k) + sir_arr(j,k)) &
                    / ((sir_arr(i,r) + sir_arr(j,r))*sir_arr(k,r))
!
  CSirsSrsTerm = 0d0
!
  n = size(p)
!
  yirQ    = yiQ_arr(radi) + yiQ_arr(radr)
  ztildei = yiQ_arr(radi)/yirQ
  ztilder = yiQ_arr(radr)/yirQ
  sir     = sir_arr(radi,radr)
!
  if (p(radi)%flv.ne.0) then
    Ti2 = qcd_cf
  else
    Ti2 = qcd_ca
  end if
!
  jpart = 0
  do j=1,n-1
    if ((j.eq.radr).or.(j.eq.rads)) cycle
    jpart = jpart + 1
    if (abs(p(j)%flv).gt.6) cycle
    lpart = jpart
    do l=j+1,n
      if ((l.eq.radr).or.(l.eq.rads)) cycle
      lpart = lpart + 1
      if (abs(p(l)%flv).gt.6) cycle
! When summation hits the ir leg a modified version of the
! Eikonal factor has to be used:
      if ((j.ne.radi).and.(l.ne.radi)) then
        CSirsSrsTerm = CSirsSrsTerm &
                     + Sikr(j,l,rads)*Bij(jpart,lpart)
      elseif (j.eq.radi) then
        CSirsSrsTerm = CSirsSrsTerm &
                     + S_ij_k_r(radi,radr,l,rads)*Bij(jpart,lpart)
      elseif (l.eq.radi) then
        CSirsSrsTerm = CSirsSrsTerm &
                     + S_ij_k_r(radi,radr,j,rads)*Bij(jpart,lpart)
      end if
    end do
  end do
!
  CSirsSrsTerm = -(8d0*pi)**2 * 2*ztildei/(sir*ztilder)*Ti2 &
               * CSirsSrsTerm
!
end subroutine CalcCSirsSrs
!
subroutine CalcCirjsSrs(p,ptilde,smeB,    &
                        emitir,radi,radr, &
                        emitjs,radj,rads, &
                        CirjsSrsTerm)
use particles
use regions
use math
use QCDparams
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
  type(particle) , dimension(:) , intent(in) :: ptilde
  real(kind(1d0)) , intent(in) :: smeB
  integer , intent(in) :: emitir,emitjs,radi,radr,radj,rads
  real(kind(1d0)) , intent(out) :: CirjsSrsTerm
!
  real(kind(1d0)) :: Ti2,Tj2
!
!
  CirjsSrsTerm = 0d0
!
  sir = sir_arr(radi,radr)
  sjs = sir_arr(radj,rads)
  yirQ = yiQ_arr(radi) + yiQ_arr(radr)
  yjsQ = yiQ_arr(radj) + yiQ_arr(rads)
  ztildei = yiQ_arr(radi)/yirQ
  ztilder = yiQ_arr(radr)/yirQ
  ztildej = yiQ_arr(radj)/yjsQ
  ztildes = yiQ_arr(rads)/yjsQ
!
  if (p(radi)%flv.ne.0) then
    Ti2 = qcd_cf
  else
    Ti2 = qcd_ca
  end if
!
  if (p(radj)%flv.ne.0) then
    Tj2 = qcd_cf
  else
    Tj2 = qcd_ca
  end if
!
  CirjsSrsTerm = (8*pi)**2*4*ztildei*ztildej &
               / (sir*sjs*ztilder*ztildes)   &
               * Ti2*Tj2*smeB
!
end subroutine CalcCirjsSrs
!
subroutine CalcCirsCSirsSrs(p,ptilde,smeB,emitirs,radi,radr,rads, &
                            CirsCSirsSrsTerm)
use particles
use regions
use math
use QCDparams
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
  type(particle) , dimension(:) , intent(in) :: ptilde
  real(kind(1d0)) , intent(in) :: smeB
  integer , intent(in) :: emitirs,radi,radr,rads
  real(kind(1d0)) , intent(out) :: CirsCSirsSrsTerm
!
  real(kind(1d0)) :: sis,srs
  real(kind(1d0)) :: Ti2
!
!
  CirsCSirsSrsTerm = 0d0
!
  srs = yrs*Q2
  sir = sir_arr(radi,radr)
  sis = sir_arr(radi,rads)
  yirsQ = yiQ_arr(radi) + yiQ_arr(radr) + yiQ_arr(rads)
  ztildei = yiQ_arr(radi)/yirsQ
  ztilder = yiQ_arr(radr)/yirsQ
  ztildes = yiQ_arr(rads)/yirsQ
!
  if (p(radi)%flv.ne.0) then
    Ti2 = qcd_cf
  else
    Ti2 = qcd_ca
  end if
!
  CirsCSirsSrsTerm = (8*pi)**2*4*ztildei*(ztildei + ztilder) &
                   / (sir*(sis + srs)*ztilder*ztildes) &
                   * Ti2**2*smeB
!
end subroutine CalcCirsCSirsSrs
!
subroutine CalcCktCktrFF(p,phat,ptilde,smeB,Bmunu, &
                         emitktr,radk,radt,radr,rad_kt,rad_r, &
                         CktCktrTerm)
use particles
use regions
use math
use utils
use misc
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p,phat,ptilde
  real(kind(1d0)) , intent(in) :: smeB
  real(kind(1d0)) , dimension(0:3,0:3) , intent(in) :: Bmunu
  integer , intent(in) :: emitktr,radk,radt,radr,rad_kt,rad_r
  real(kind(1d0)) , intent(out) :: CktCktrTerm
!
  integer :: i
  real(kind(1d0)) :: srkt,kt12,kt22,kt_tilde_k2
  type(particle) :: kt_tilde_k
  real(kind(1d0)) , dimension(0:3,0:3) :: Bmunu_tmp,dmunu
!
  real(kind(1d0)) , external :: Pqg0,Pqq0, &
                                PqqbqSO0,PqggSO0,PgqqbSO0,PgggSO0,Pggg0
!
  real(kind(1d0)) :: skt,skr,str
!
  interface
    subroutine CalcCggCgggSG(p,phat,ptilde,radi,radr,rads,CrsCirs)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p,phat,ptilde
      integer , intent(in) :: radi,radr,rads
      real(kind(1d0)) , intent(out) :: CrsCirs
!
    end subroutine CalcCggCgggSG
  end interface
!
  CktCktrTerm = 0d0
!
! The ordering among radiated parton is: k,t,r
! q(~) -> r r~ q(~):
  if ((p(radk)%flv.ne.0).and. &
      (p(radt)%flv.ne.0).and. &
      (p(radr)%flv.ne.0)) then
!    print *,"q(~) -> r r~ q(~)"
!
! The original strongly ordered splitting kernel is formulated
! for q_r \bar{q}_k' q_t'
! s_{\hat{r} k_{\bot,k,t}}:
    srkt = 2*phat(rad_r)%p*kttildek%p
! k_{\bot,k,t}^2:
    kt12 = kttildek%p**2
!
    CktCktrTerm = PqqbqSO0(ztildek,ztildet,zhatkt,zhatr,s_kt_r,srkt,kt12)*smeB
!
! q(~) -> q(~) g g:
  elseif ((p(radk)%flv.ne.0).and. &
          (p(radt)%flv.eq.0).and. &
          (p(radr)%flv.eq.0)) then
!    print *,"q(~) -> q(~) g g"
!
    CktCktrTerm = Pqg0(ztildek,ztildet)*Pqg0(zhatkt,zhatr)*smeB
!
! q(~) -> g g q(~):
  elseif ((p(radk)%flv.eq.0).and. &
          (p(radt)%flv.eq.0).and. &
          (p(radr)%flv.ne.0)) then
!    print *,"q(~) -> g g q(~)"
!
! s_{\hat{r} k_{\bot,k,t}}:
    srkt = 2*phat(rad_r)%p*kttildek%p
! k_{\bot,k,t}^2:
    kt12 = kttildek%p**2
!
    CktCktrTerm = PqggSO0(ztildek,ztildet,zhatkt,zhatr,s_kt_r,srkt,kt12)*smeB
!
! g -> q q~ g:
  elseif ((p(radk)%flv.ne.0).and. &
          (p(radt)%flv.ne.0).and. &
          (p(radr)%flv.eq.0)) then
!    print *,"g -> q q~ g"
!
    kt_tilde_k = kttildek
    call ChangeKt(kt_tilde_k%p,ptilde(emitktr)%p,Q)
!
! s_{\hat{r} k_{\bot,k,t}}:
    srkt = 2*phat(rad_r)%p*kttildek%p
! k_{\bot,k,t}^2:
    kt12 = kttildek%p**2
! k_{\bot,\hat{r},\hat{kt}}^2:
    kt22 = kthatr%p**2
! \tilde{k}_{\bot,k,t}^2:
    kt_tilde_k2 = kt_tilde_k%p**2
!
    CktCktrTerm = PgqqbSO0(ztildek,ztildet,zhatkt,zhatr, &
                           kt_tilde_k,kt_tilde_k2,kthatr, &
                           s_kt_r,srkt,kt12,kt22,Bmunu,Bmunu)
!
! g -> q g q~ or g -> q~ g q 
  elseif ((p(radk)%flv.ne.0).and. &
          (p(radt)%flv.eq.0).and. &
          (p(radr)%flv.ne.0)) then
!    print *,"g -> q g q~"
!
    CktCktrTerm = Pqg0(ztildek,ztildet) &
                * Pqq0(zhatkt,zhatr,kthatkt,Bmunu)
!
! g -> g g g:
  elseif ((p(radk)%flv.eq.0).and. &
          (p(radt)%flv.eq.0).and. &
          (p(radr)%flv.eq.0)) then
!    print *,"g -> g g g"
!
    kt_tilde_k = kttildek
    call ChangeKt(kt_tilde_k%p,ptilde(emitktr)%p,Q)
!
! s_{\hat{r} k_{\bot,k,t}}:
    srkt = 2*phat(rad_r)%p*kttildek%p
! k_{\bot,k,t}^2:
    kt12 = kttildek%p**2
! k_{\bot,\hat{r},\hat{kt}}^2:
    kt22 = kthatr%p**2
! \tilde{k}_{\bot,k,t}^2:
    kt_tilde_k2 = kt_tilde_k%p**2
!
    CktCktrTerm = PgggSO0(ztildek,ztildet,zhatkt,zhatr,kt_tilde_k,kt_tilde_k2, &
                          kthatr,s_kt_r,srkt,kt12,kt22,Bmunu)
!
!    print *,"*****"
!    print *,"CktCktrTerm: ",CktCktrTerm
!    call CalcCggCgggSG(p,phat,ptilde,radr,radk,radt,CktCktrTerm)
!    print *,"CktCktrTerm: ",CktCktrTerm
!
  else
    print *,"Error in CalcCktCktrFF..."
    print *,"radk,radt,radr: ",radk,radt,radr
    print *,"IDs: ",p(radk)%flv,p(radt)%flv,p(radr)%flv
    stop "CalcCktCktrFF"
  end if
!
  CktCktrTerm = (8d0*pi)**2/(s_k_t*s_kt_r)*CktCktrTerm
!
end subroutine CalcCktCktrFF
!
subroutine CalcCktCirktFF(p,phat,ptilde,smeB,Bmunu1,Bmunu2,Balbemunu, &
                          radi,radr,radk,radt, &
                          CktCirktTerm)
use particles
use regions
use math
use utils
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p,phat,ptilde
  real(kind(1d0)) , intent(in) :: smeB
  real(kind(1d0)) , dimension(0:3,0:3) , intent(in) :: Bmunu1,Bmunu2
  real(kind(1d0)) , dimension(0:3,0:3,0:3,0:3) , intent(in) :: Balbemunu
  integer , intent(in) :: radi,radr,radk,radt
  real(kind(1d0)) , intent(out) :: CktCirktTerm
!
!
  real(kind(1d0)) , external :: Pqg0,Pqq0,Pgg0
!
  CktCirktTerm = 0d0
!
! q(~) -> q(~) + g && r(~) -> r(~) + g
  if (((p(radi)%flv+p(radr)%flv).ne.0).and. &
      ((p(radk)%flv+p(radt)%flv).ne.0)) then
!
    CktCirktTerm = Pqg0(zhati,zhatr) &
                 * Pqg0(ztildek,ztildet) * smeB
!
! q(~) -> q(~) + g && g -> r + r~
  elseif (((p(radi)%flv+p(radr)%flv).ne.0).and. &
          ((p(radk)%flv+p(radt)%flv).eq.0).and. &
          (p(radk)%flv.ne.0)) then
!
    CktCirktTerm = Pqg0(zhati,zhatr) &
                 * Pqq0(ztildek,ztildet,kttildek,Bmunu2)
!
! q(~) -> q(~) + g && g -> g + g
  elseif (((p(radi)%flv+p(radr)%flv).ne.0).and. &
          ((p(radk)%flv+p(radt)%flv).eq.0).and. &
          (p(radk)%flv.eq.0)) then
!
    CktCirktTerm = Pqg0(zhati,zhatr) &
                 * Pgg0(ztildek,ztildet,kttildek,Bmunu2)
!
! g -> q + q~ && r(~) -> r(~) + g
  elseif (((p(radi)%flv+p(radr)%flv).eq.0).and. &
          ((p(radk)%flv+p(radt)%flv).ne.0).and. &
          (p(radi)%flv.ne.0)) then
!
    CktCirktTerm = Pqq0(zhati,zhatr,kthati,Bmunu1) &
                 * Pqg0(ztildek,ztildet)
!
! g -> g + g && r(~) -> r(~) + g
  elseif (((p(radi)%flv+p(radr)%flv).eq.0).and. &
          ((p(radk)%flv+p(radt)%flv).ne.0).and. &
          (p(radi)%flv.eq.0)) then
!
    CktCirktTerm = Pgg0(zhati,zhatr,kthati,Bmunu1) &
                 * Pqg0(ztildek,ztildet)
!
! g -> q + q~ && g -> r + r~
  elseif (((p(radi)%flv+p(radr)%flv).eq.0).and. &
          ((p(radk)%flv+p(radt)%flv).eq.0).and. &
          (p(radi)%flv.ne.0).and.(p(radk)%flv.ne.0)) then
!
    print *,"Not implemented yet..."
    stop "CalcCktCirkt..."
!
! g -> q + q~ && g -> g + g
  elseif (((p(radi)%flv+p(radr)%flv).eq.0).and. &
          ((p(radk)%flv+p(radt)%flv).eq.0).and. &
          (p(radi)%flv.ne.0).and.(p(radk)%flv.eq.0)) then
!
    print *,"Not implemented yet..."
    stop "CalcCktCirkt..."
!
! g -> g + g && g -> q + q~
  elseif (((p(radi)%flv+p(radr)%flv).eq.0).and. &
          ((p(radk)%flv+p(radt)%flv).eq.0).and. &
          (p(radi)%flv.eq.0).and.(p(radk)%flv.ne.0)) then
!
    print *,"Not implemented yet..."
    stop "CalcCktCirkt..."
!
! g -> g + g && g -> g + g~
  elseif (((p(radi)%flv+p(radr)%flv).eq.0).and. &
          ((p(radk)%flv+p(radt)%flv).eq.0).and. &
          (p(radi)%flv.eq.0).and.(p(radk)%flv.eq.0)) then
!
    print *,"Not implemented yet..."
    stop "CalcCktCirkt..."
!
  else
    print *,"Error in CalcCktCirktFF..."
    print *,"radi,radr,radk,radt: ",radi,radr,radk,radt
    print *,"IDs: ",p(radi)%flv,p(radr)%flv,p(radk)%flv,p(radt)%flv
    stop "CalcCktCirktFF"
  end if
!
  CktCirktTerm = (8d0*pi)**2/(s_k_t*sir)*CktCirktTerm
!
end subroutine CalcCktCirktFF
!
subroutine CalcCktCSktrFF(p,phat,ptilde,Bij,Bmunuij, &
                          radk,radt,rad_r,CktCSktrTerm)
use particles
use regions
use math
use utils
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p,phat,ptilde
  real(kind(1d0)) , dimension(:,:) , intent(in) :: Bij
  real(kind(1d0)) , dimension(0:,0:,:,:) , intent(in) :: Bmunuij
  integer , intent(in) :: radk,radt,rad_r
  real(kind(1d0)) , intent(out) :: CktCSktrTerm
!
  integer :: j,l,r,n,jpart,lpart,kt
  real(kind(1d0)) :: Pqg,Pqq,Pgg
  real(kind(1d0)) , dimension(0:3,0:3) :: Bmunu
!
  real(kind(1d0)) :: Sjl
  real(kind(1d0)) , external :: Pqg0,Pqq0,Pgg0
!
! This eikonal factor is defined on the hatted momenta:
  Sjl(j,l,r) = 2d0*sirhat_arr(j,l)/(sirhat_arr(j,r)*sirhat_arr(l,r))
!
  CktCSktrTerm = 0d0
!
  n = size(p)
  r = rad_r
  kt = min(radk,radt)
!
! Since parton r can only be a gluon it is sufficient to 
! test only radk and radt:
! q(~) -> q(~) + g
  if ((p(radk)%flv.ne.0).and.(p(radt)%flv.eq.0)) then
!
    Pqg = Pqg0(ztildek,ztildet)
!
    jpart = 0
    do j=1,n-1
      if (j.eq.r) cycle
      jpart = jpart + 1
      if (abs(phat(j)%flv).gt.6) cycle
      lpart = jpart
      do l=j+1,n-1
        if (l.eq.r) cycle
        lpart = lpart + 1
        if (abs(phat(l)%flv).gt.6) cycle
        CktCSktrTerm = CktCSktrTerm + Sjl(j,l,r)*Pqg*Bij(jpart,lpart)
      end do
    end do
!
! g -> g + g || g -> q + q~
  elseif (((p(radk)%flv.eq.0).and.(p(radt)%flv.eq.0)).or. &
          ((p(radk)%flv.ne.0).and.(p(radt)%flv.ne.0))) then
!
    Bmunu = 0
!
    jpart = 0
    do j=1,n-1
      if (j.eq.r) cycle
      jpart = jpart + 1
      if (abs(phat(j)%flv).gt.6) cycle
      lpart = jpart
      do l=j+1,n-1
        if (l.eq.r) cycle
        lpart = lpart + 1
        if (abs(phat(l)%flv).gt.6) cycle
        Bmunu = Bmunu + Sjl(j,l,r)*Bmunuij(:,:,jpart,lpart)
      end do
    end do
!
! g -> q + q~
    if (p(radk)%flv.ne.0) then
      Pqq = Pqq0(ztildek,ztildet,kttildek,Bmunu)
      CktCSktrTerm = Pqq
!
! g -> g + g
    else
      Pgg = Pgg0(ztildek,ztildet,kttildek,Bmunu)
      CktCSktrTerm = Pgg
    end if
!
  else
    print *,"Error in CalcCktCSktrFF..."
    print *,"radk,radt,rad_r: ",radk,radt,rad_r
    print *,"IDs: ",p(radk)%flv,p(radt)%flv,phat(rad_r)%flv
    stop "CalcCktCSktrFF"
  end if
!
  CktCSktrTerm = -(8d0*pi)**2/s_k_t*CktCSktrTerm
!
end subroutine CalcCktCSktrFF
!
! In rad_i and rad_r the underscore reminds us that those
! indices are defined on hatted momenta (phat).
subroutine CalcCktCirktCSktrFF(p,phat,ptilde,smeB,Bmunu, &
                               rad_i,rad_r,radk,radt, &
                               CktCirktCSktrTerm)
use particles
use regions
use math
use utils
use QCDparams
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p,phat,ptilde
  real(kind(1d0)) , intent(in) :: smeB
  real(kind(1d0)) , dimension(0:3,0:3) , intent(in) :: Bmunu
  integer , intent(in) :: rad_i,rad_r,radk,radt
  real(kind(1d0)) , intent(out) :: CktCirktCSktrTerm
!
  real(kind(1d0)) :: Ti2
  real(kind(1d0)) :: zi,s_i_r
!
  real(kind(1d0)) , external :: Pqg0,Pqq0,Pgg0
!
  CktCirktCSktrTerm = 0d0
!
  zi = yihatQ_arr(rad_i)/(yihatQ_arr(rad_i) + yihatQ_arr(rad_r))
  s_i_r = sirhat_arr(rad_i,rad_r)
!
  if (phat(rad_i)%flv.ne.0) Ti2 = qcd_cf
  if (phat(rad_i)%flv.eq.0) Ti2 = qcd_ca
!
! q(~) -> q(~) + g
  if ((p(radk)%flv.ne.0).and. &
      (p(radt)%flv.eq.0)) then
!
    CktCirktCSktrTerm = 2*zi/(1 - zi)*Ti2 &
                      * Pqg0(ztildek,ztildet)*smeB
!
! g -> q + q~
  elseif  ((p(radk)%flv.ne.0).and. &
           (p(radk)%flv.eq.-p(radt)%flv)) then
!
    CktCirktCSktrTerm = 2*zi/(1 - zi)*Ti2 &
                      * Pqq0(ztildek,ztildet,kttildek,Bmunu)
!
! g -> g + g
  elseif ((p(radk)%flv.eq.0).and. &
          (p(radt)%flv.eq.0)) then
!
    CktCirktCSktrTerm = 2*zi/(1 - zi)*Ti2 &
                      * Pgg0(ztildek,ztildet,kttildek,Bmunu)
!
  else
    print *,"Error in CalcCktCirktCSktr..."
    print *,"rad_i,rad_r,radk,radt: ",rad_i,rad_r,radk,radt
    call PrintSubProc(p)
    stop "CalcCktCirktCSktr"
  end if
!
  CktCirktCSktrTerm = (8d0*pi)**2/(s_k_t*s_i_r)*CktCirktCSktrTerm
!
end subroutine CalcCktCirktCSktrFF
!
! In rad_kt and rad_r the underscore reminds us that those
! indices are defined on hatted momenta (phat).
subroutine CalcCktCktrCSktrFF(p,phat,ptilde,smeB,Bmunu, &
                              rad_kt,rad_r,radk,radt, &
                              CktCktrCSktrTerm)
use particles
use regions
use math
use utils
use QCDparams
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p,phat,ptilde
  real(kind(1d0)) , intent(in) :: smeB
  real(kind(1d0)) , dimension(0:3,0:3) , intent(in) :: Bmunu
  integer , intent(in) :: rad_kt,rad_r,radk,radt
  real(kind(1d0)) , intent(out) :: CktCktrCSktrTerm
!
  real(kind(1d0)) :: Tkt2
  real(kind(1d0)) :: zkt
!
  real(kind(1d0)) , external :: Pqg0,Pqq0,Pgg0
!
  CktCktrCSktrTerm = 0d0
!
  zkt = yihatQ_arr(rad_kt)/(yihatQ_arr(rad_kt) + yihatQ_arr(rad_r))
  s_kt_r = sirhat_arr(rad_kt,rad_r)
!
  if (phat(rad_kt)%flv.ne.0) Tkt2 = qcd_cf
  if (phat(rad_kt)%flv.eq.0) Tkt2 = qcd_ca
!
! q(~) -> q(~) + g
  if ((p(radk)%flv.ne.0).and. &
      (p(radt)%flv.eq.0)) then
!
    CktCktrCSktrTerm = 2*zkt/(1 - zkt)*Tkt2 &
                     * Pqg0(ztildek,ztildet)*smeB
!
! g -> q + q~
  elseif  ((p(radk)%flv.ne.0).and. &
           (p(radk)%flv.eq.-p(radt)%flv)) then
!
    CktCktrCSktrTerm = 2*zkt/(1 - zkt)*Tkt2 &
                     * Pqq0(ztildek,ztildet,kttildek,Bmunu)
!
! g -> g + g
  elseif ((p(radk)%flv.eq.0).and. &
          (p(radt)%flv.eq.0)) then
!
    CktCktrCSktrTerm = 2*zkt/(1 - zkt)*Tkt2 &
                     * Pgg0(ztildek,ztildet,kttildek,Bmunu)
!
  else
    print *,"Error in CalcCktCktrCSktr..."
    print *,"rad_kt,rad_r,radk,radt: ",rad_kt,rad_r,radk,radt
    call PrintSubProc(p)
    stop "CalcCktCktrCSktr"
  end if
!
  CktCktrCSktrTerm = (8d0*pi)**2/(s_k_t*s_kt_r)*CktCktrCSktrTerm
!
end subroutine CalcCktCktrCSktrFF
!
subroutine CalcCktSktFF(p,phat,ptilde,Bij, &
                        rad_kt,radk,radt,  &
                        CktSktTerm)
use particles
use regions
use math
use utils
use QCDparams
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p,phat,ptilde
  real(kind(1d0)) , dimension(:,:) , intent(in) :: Bij
  integer , intent(in) :: rad_kt,radk,radt
  real(kind(1d0)) , intent(out) :: CktSktTerm
!
  integer :: n,j,l,jpart,lpart
  real(kind(1d0)) , dimension(0:3,0:3) :: Smunu,Smunu_tmp
!
  real(kind(1d0)) , external :: Pqq0,Pgg0
!
  CktSktTerm = 0d0
!
! Note that n is the number of particles for hatted momenta!!!!
  n = size(phat)
!
! Calculating \mathcal{S}^{\mu\nu}_{\hat{j}\hat{l}} which definition
! can be found in Eq. (7.27) in arXiv:hep-ph/0609042:
  Smunu = 0
  jpart = 0
  do j=1,n
    if (j.eq.rad_kt) cycle
    jpart = jpart + 1
    if (abs(phat(j)%flv).gt.6) cycle
    Smunu_tmp = 0
    Smunu_tmp(0,0) = phat(j)%p%E*phat(j)%p%E
    Smunu_tmp(0,1) = phat(j)%p%E*phat(j)%p%px
    Smunu_tmp(0,2) = phat(j)%p%E*phat(j)%p%py
    Smunu_tmp(0,3) = phat(j)%p%E*phat(j)%p%pz
    Smunu_tmp(1,1) = phat(j)%p%px*phat(j)%p%px
    Smunu_tmp(1,2) = phat(j)%p%px*phat(j)%p%py
    Smunu_tmp(1,3) = phat(j)%p%px*phat(j)%p%pz
    Smunu_tmp(2,2) = phat(j)%p%py*phat(j)%p%py
    Smunu_tmp(2,3) = phat(j)%p%py*phat(j)%p%pz
    Smunu_tmp(3,3) = phat(j)%p%pz*phat(j)%p%pz
    Smunu = Smunu + 2*Smunu_tmp &
          / (sirhat_arr(rad_kt,j)*sirhat_arr(rad_kt,j)) &
          * Bij(jpart,jpart)
    lpart = jpart
    do l=j+1,n
      if (l.eq.rad_kt) cycle
      lpart = lpart + 1
      if (abs(phat(l)%flv).gt.6) cycle
      Smunu_tmp(0,0) = phat(j)%p%E*phat(l)%p%E &
                     + phat(l)%p%E*phat(j)%p%E
      Smunu_tmp(0,1) = phat(j)%p%E*phat(l)%p%px &
                     + phat(l)%p%E*phat(j)%p%px
      Smunu_tmp(0,2) = phat(j)%p%E*phat(l)%p%py &
                     + phat(l)%p%E*phat(j)%p%py
      Smunu_tmp(0,3) = phat(j)%p%E*phat(l)%p%pz &
                     + phat(l)%p%E*phat(j)%p%pz
      Smunu_tmp(1,1) = phat(j)%p%px*phat(l)%p%px &
                     + phat(l)%p%px*phat(j)%p%px
      Smunu_tmp(1,2) = phat(j)%p%px*phat(l)%p%py &
                     + phat(l)%p%px*phat(j)%p%py
      Smunu_tmp(1,3) = phat(j)%p%px*phat(l)%p%pz &
                     + phat(l)%p%px*phat(j)%p%pz
      Smunu_tmp(2,2) = phat(j)%p%py*phat(l)%p%py &
                     + phat(l)%p%py*phat(j)%p%py
      Smunu_tmp(2,3) = phat(j)%p%py*phat(l)%p%pz &
                     + phat(l)%p%py*phat(j)%p%pz
      Smunu_tmp(3,3) = phat(j)%p%pz*phat(l)%p%pz &
                     + phat(l)%p%pz*phat(j)%p%pz
      Smunu = Smunu + 2*Smunu_tmp &
            / (sirhat_arr(rad_kt,j)*sirhat_arr(rad_kt,l)) &
            * Bij(jpart,lpart)
    end do
  end do
  Smunu(1,0) = Smunu(0,1)
  Smunu(2,0) = Smunu(0,2)
  Smunu(3,0) = Smunu(0,3)
  Smunu(2,1) = Smunu(1,2)
  Smunu(3,1) = Smunu(1,3)
  Smunu(3,2) = Smunu(2,3)
!
! g -> q + q~
  if  ((p(radk)%flv.ne.0).and. &
       (p(radk)%flv.eq.-p(radt)%flv)) then
!
    CktSktTerm = Pqq0(ztildek,ztildet,kttildek,Smunu)
!
! g -> g + g
  elseif ((p(radk)%flv.eq.0).and. &
          (p(radt)%flv.eq.0)) then
!
    CktSktTerm = Pgg0(ztildek,ztildet,kttildek,Smunu)
!
  else
    print *,"Error in CalcCktSkt..."
    print *,"rad_kt,radk,radt: ",rad_kt,radk,radt
    call PrintSubProc(p)
    stop "CalcCktSkt"
  end if
!
  CktSktTerm = (8d0*pi)**2/s_k_t*CktSktTerm
!
end subroutine CalcCktSktFF
!
subroutine CalcCktCrktSktFF(p,phat,ptilde,smeB,          &
                            rad_r,rad_kt,radr,radk,radt, &
                            CktCrktSktTerm)
use particles
use regions
use math
use utils
use QCDparams
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p,phat,ptilde
  integer , intent(in) :: rad_r,rad_kt,radr,radk,radt
  real(kind(1d0)) , intent(in) :: smeB
  real(kind(1d0)) , intent(out) :: CktCrktSktTerm
!
  real(kind(1d0)) :: Tr2
  real(kind(1d0)) :: s_r_kt_k_t,kt_k_t2
!
  if (p(radr)%flv.ne.0) then
    Tr2 = qcd_cf
  else
    Tr2 = qcd_ca
  end if
!
! s_kt_r : s_{\hat{kt},\hat{r}}
! zhatkt: z_{\hat{kt}}
! zhatr: z_{\hat{r}}
! s_r_kt_k_t : s_{\hat{r},k_{\bot,k,t}}
! kt_k_t2 : k_{\bot,k,t}^2
!  s_kt_r = 2*phat(rad_kt)%p*phat(rad_r)%p
  s_kt_r = sirhat_arr(rad_kt,rad_r)
  zhatkt = yihatQ_arr(rad_kt)/(yihatQ_arr(rad_kt) + yihatQ_arr(rad_r))
  zhatr  = 1 - zhatkt
  s_r_kt_k_t = 2*phat(rad_r)%p*kttildek%p
  kt_k_t2 = kttildek%p**2
!
  CktCrktSktTerm = 0d0
!
! g -> q + q~
  if  ((p(radk)%flv.ne.0).and. &
       (p(radk)%flv.eq.-p(radt)%flv)) then
!
    CktCrktSktTerm = 2*Tr2*qcd_tr &
                   * (zhatr/zhatkt + ztildek*ztildet &
                   *  s_r_kt_k_t**2/(kt_k_t2*s_kt_r)) &
                   * smeB
!
! g -> g + g
  elseif ((p(radk)%flv.eq.0).and. &
          (p(radt)%flv.eq.0)) then
!
    CktCrktSktTerm = 2*Tr2*qcd_ca &
                   * (2*zhatr/zhatkt &
                   *  (ztildek/ztildet + ztildet/ztildek) &
                   -  ztildek*ztildet &
                   *  s_r_kt_k_t**2/(kt_k_t2*s_kt_r)) &
                   * smeB
!
  else
    print *,"Error in CalcCktCrktSkt..."
    print *,"rad_r,rad_kt,radk,radt: ",rad_r,rad_kt,radk,radt
    call PrintSubProc(p)
    stop "CalcCktCrktSkt"
  end if
!
  CktCrktSktTerm = (8d0*pi)**2/(s_k_t*s_kt_r)*CktCrktSktTerm
!
end subroutine CalcCktCrktSktFF
!
subroutine CalcStCirtFF(p,phat,ptilde,smeB,Bmunu, &
                        radi,radr,radt,       &
                        StCirtTerm)
use particles
use regions
use math
use utils
use QCDparams
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p,phat,ptilde
  integer , intent(in) :: radi,radr,radt
  real(kind(1d0)) , intent(in) :: smeB
  real(kind(1d0)) , dimension(0:3,0:3) , intent(in) :: Bmunu
  real(kind(1d0)) , intent(out) :: StCirtTerm
!
  real(kind(1d0)) :: yirtQ
  real(kind(1d0)) :: zi,zr,zt,s_i_r,s_i_t,s_r_t
!
  real(kind(1d0)) , external :: Pqg0,Pqq0,Pgg0,PSqgg,PSgqq,PSggg
!
  StCirtTerm = 0d0
!
  yirtQ = yiQ_arr(radi) + yiQ_arr(radr) + yiQ_arr(radt)
!
  zi = yiQ_arr(radi)/yirtQ
  zr = yiQ_arr(radr)/yirtQ
  zt = 1 - zi - zr
!
!  s_i_r = 2*p(radi)%p*p(radr)%p
!  s_i_t = 2*p(radi)%p*p(radt)%p
!  s_r_t = 2*p(radr)%p*p(radt)%p
  s_i_r = sir_arr(radi,radr)
  s_i_t = sir_arr(radi,radt)
  s_r_t = sir_arr(radr,radt)
!
! Note that the ordering in the argument of the various P^{(S)}
! functions is different in arXiv:hep-ph/0609042 compared to
! arXiv:hep-ph/0502226:
!
! Note that no criterion is put on parton t which is always
! a gluon:
! q(~) -> q(~) + g + g
  if ((p(radi)%flv.ne.0).and. &
      (p(radr)%flv.eq.0)) then
!
! The argument is ordered differently, to get agreement:
! r -> t ; s -> r
    StCirtTerm = PSqgg(zi,zt,zr,s_i_t,s_i_r,s_r_t) &
               * Pqg0(zhati,zhatr)*smeB
!
! g -> q + q~ + g
  elseif ((p(radi)%flv.ne.0).and. &
          (p(radr)%flv.ne.0).and. &
          (p(radi)%flv.eq.-p(radr)%flv)) then
!
! The argument is ordered differently, to get agreement:
! s -> i ; i -> r ; r -> t
    StCirtTerm = PSgqq(zr,zt,zi,s_r_t,s_i_r,s_i_t) &
               * Pqq0(zhati,zhatr,kthati,Bmunu)
!
! g -> g + g + g
  elseif ((p(radi)%flv.eq.0).and. &
          (p(radr)%flv.eq.0)) then
!
! The argument is ordered differently, to get agreement:
! r -> t s -> r
    StCirtTerm = PSggg(zi,zt,zr,s_i_t,s_i_r,s_r_t) &
               * Pgg0(zhati,zhatr,kthati,Bmunu)
!
  else
    print *,"Error in CalcStCirt..."
    print *,"radi,radr,radt: ",radi,radr,radt
    call PrintSubProc(p)
    stop "CalcStCirt"
  end if
!
  StCirtTerm = (8d0*pi)**2/sir*StCirtTerm
!
end subroutine CalcStCirtFF
!
subroutine CalcStCSirtFF(p,phat,ptilde,Bij,Bmunuij, &
                         emit_ir,radi,radr,radt,    &
                         StCSirtTerm)
use particles
use regions
use math
use utils
use QCDparams
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p,phat,ptilde
  integer , intent(in) :: emit_ir,radi,radr,radt
  real(kind(1d0)) , dimension(:,:) , intent(in) :: Bij
  real(kind(1d0)) , dimension(0:,0:,:,:) , intent(in) :: Bmunuij
  real(kind(1d0)) , intent(out) :: StCSirtTerm
!
  integer :: j,l,s,n,jpart,lpart
  integer :: emit,rad
  real(kind(1d0)) :: Pqg,Pqq,Pgg
  real(kind(1d0)) , dimension(0:3,0:3) :: Bmunu
!
  real(kind(1d0)) :: Sjl
  real(kind(1d0)) , external :: Pqg0,Pqq0,Pgg0
!
  Sjl(j,l,s) = 2d0*sir_arr(j,l)/(sir_arr(j,s)*sir_arr(l,s))
!
  StCSirtTerm = 0d0
!
  n = size(p)
! The position of the radiated parton originating from the 
! collinear pair:
  rad  = max(radi,radr)
  emit = min(radi,radr)
!
! q(~) -> q(~) + g + g
  if ((p(radi)%flv.ne.0).and. &
      (p(radr)%flv.eq.0)) then
!
    Pqg = Pqg0(zhati,zhatr)
    jpart = 0
    do j=1,n
      if ((j.eq.radt).or.(j.eq.rad)) cycle
      jpart = jpart + 1
      if (j.eq.emit) cycle
      if (abs(p(j)%flv).gt.6) cycle
      lpart = jpart
      do l=j+1,n
        if ((l.eq.radt).or.(l.eq.rad)) cycle
        lpart = lpart + 1
        if (l.eq.emit) cycle
        if (abs(p(l)%flv).gt.6) cycle
        StCSirtTerm = StCSirtTerm &
                    + Sjl(j,l,radt)*Pqg*Bij(jpart,lpart)
      end do
    end do
! Contribution is also coming from the (ir) leg, but for this
! particular one the form of \mathcal{S} changes,
! for the precise definition take a look at Eq. (6.26) of 
! arXiv:hep-ph/0609042.
    lpart = 0
    do l=1,n
      if ((l.eq.radt).or.(l.eq.rad)) cycle
      lpart = lpart + 1
      if (l.eq.emit) cycle
      if (abs(p(l)%flv).gt.6) cycle
      StCSirtTerm = StCSirtTerm &
                  + 2*(sir_arr(radi,l) + sir_arr(radr,l)) &
                  / ((sir_arr(radi,radt) + sir_arr(radr,radt))*sir_arr(l,radt)) &
                  * Pqg*Bij(emit_ir,lpart)
    end do
!
! g -> q + q~ + g || g -> g + g + g
  elseif (p(radi)%flv.eq.-p(radr)%flv) then
!
    Bmunu = 0d0
    jpart = 0
    do j=1,n
      if ((j.eq.radt).or.(j.eq.rad)) cycle
      jpart = jpart + 1
      if (j.eq.emit) cycle
      if (abs(p(j)%flv).gt.6) cycle
      lpart = jpart
      do l=j+1,n
        if ((l.eq.radt).or.(l.eq.rad)) cycle
        lpart = lpart + 1
        if (l.eq.emit) cycle
        if (abs(p(l)%flv).gt.6) cycle
        Bmunu = Bmunu + Sjl(j,l,radt)*Bmunuij(:,:,jpart,lpart)
      end do
    end do
    lpart = 0
    do l=1,n
      if ((l.eq.radt).or.(l.eq.rad)) cycle
      lpart = lpart + 1
      if (l.eq.emit) cycle
      if (abs(p(l)%flv).gt.6) cycle
      Bmunu = Bmunu &
            + 2*(sir_arr(radi,l) + sir_arr(radr,l)) &
            /((sir_arr(radi,radt) + sir_arr(radr,radt))*sir_arr(l,radt)) &
            * Bmunuij(:,:,emit_ir,lpart)
    end do
! g -> q + q~
    if (p(radi)%flv.ne.0) then
      Pqq = Pqq0(zhati,zhatr,kthati,Bmunu)
      StCSirtTerm = Pqq
!
! g -> g + g
    else
      Pgg = Pgg0(zhati,zhatr,kthati,Bmunu)
      StCSirtTerm = Pgg
    end if
! 
  else
    print *,"Error in CalcStCSirt..."
    print *,"radi,radr,radt: ",radi,radr,radt
    call PrintSubProc(p)
    stop "CalcStCSirt"
  end if
!
  StCSirtTerm = -(8d0*pi)**2/sir*StCSirtTerm
!
end subroutine CalcStCSirtFF
!
subroutine CalcStCirtCSirtFF(p,phat,ptilde,smeB,Bmunu, &
                             emit_ir,radi,radr,radt,   &
                             StCirtCSirtTerm)
use particles
use regions
use math
use utils
use QCDparams
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p,phat,ptilde
  real(kind(1d0)) , intent(in) :: smeB
  real(kind(1d0)) , dimension(0:3,0:3) , intent(in) :: Bmunu
  integer , intent(in) :: emit_ir,radi,radr,radt
  real(kind(1d0)) , intent(out) :: StCirtCSirtTerm
!
  real(kind(1d0)) :: Tir2
  real(kind(1d0)) :: yirtQ
  real(kind(1d0)) :: zt,s_ir_t
!
  real(kind(1d0)) , external :: Pqg0,Pqq0,Pgg0
!
  StCirtCSirtTerm = 0d0
!
  yirtQ = yiQ_arr(radi) + yiQ_arr(radr) + yiQ_arr(radt)
!
  zt = yiQ_arr(radt)/yirtQ
!
!  s_ir_t = 2*(p(radi)%p + p(radr)%p)*p(radt)%p
  s_ir_t = sir_arr(radi,radt) + sir_arr(radr,radt)
!
  if (ptilde(emit_ir)%flv.ne.0) then
    Tir2 = qcd_cf
  else
    Tir2 = qcd_ca
  end if
!
! q(~) -> q(~) + g + g
  if ((p(radi)%flv.ne.0).and. &
      (p(radr)%flv.eq.0)) then
!
    StCirtCSirtTerm = 2*(1 - zt)/(s_ir_t*zt)*Tir2 &
                    * Pqg0(zhati,zhatr)*smeB
!
! g -> q + q~ + g
  elseif ((p(radi)%flv.ne.0).and. &
          (p(radr)%flv.ne.0).and. &
          (p(radi)%flv.eq.-p(radr)%flv)) then
!
    StCirtCSirtTerm = 2*(1 - zt)/(s_ir_t*zt)*Tir2 &
                    * Pqq0(zhati,zhatr,kthati,Bmunu)
!
! g -> g + g + g
  elseif ((p(radi)%flv.eq.0).and. &
          (p(radr)%flv.eq.0)) then
!
    StCirtCSirtTerm = 2*(1 - zt)/(s_ir_t*zt)*Tir2 &
                    * Pgg0(zhati,zhatr,kthati,Bmunu)
!
  else
    print *,"Error in CalcStCirtCSirt..."
    print *,"radi,radr,radt: ",radi,radr,radt
    call PrintSubProc(p)
    stop "CalcStCirtCSirt"
  end if
!
  StCirtCSirtTerm = (8d0*pi)**2/sir*StCirtCSirtTerm
!
end subroutine CalcStCirtCSirtFF
!
subroutine CalcStCirtSrtFF(p,phat,ptilde,smeB,         &
                           rad_i,rad_r,radi,radr,radt, &
                           StCirtSrtTerm)
use particles
use regions
use math
use utils
use QCDparams
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p,phat,ptilde
  real(kind(1d0)) , intent(in) :: smeB
  integer , intent(in) :: rad_i,rad_r,radi,radr,radt
  real(kind(1d0)) , intent(out) :: StCirtSrtTerm
!
  real(kind(1d0)) :: Ti2
  real(kind(1d0)) :: yirtQ,yirhatQ
  real(kind(1d0)) :: zi,zr,zt,s_i_r,s_i_t,s_r_t
!
!
  StCirtSrtTerm = 0d0
!
  yirtQ = yiQ_arr(radi) + yiQ_arr(radr) + yiQ_arr(radt)
!
! Quantities defined with the help of the original momenta:
  zi = yiQ_arr(radi)/yirtQ
  zr = yiQ_arr(radr)/yirtQ
  zt = 1 - zi - zr
!
!  s_i_r = 2*p(radi)%p*p(radr)%p
!  s_i_t = 2*p(radi)%p*p(radt)%p
!  s_r_t = 2*p(radr)%p*p(radt)%p
  s_i_r = sir_arr(radi,radr)
  s_i_t = sir_arr(radi,radt)
  s_r_t = sir_arr(radr,radt)
!
! Quantities defined with the help of the hatted momenta:
  yirhatQ = yihatQ_arr(rad_i) + yihatQ_arr(rad_r)
!
  zhati = yihatQ_arr(rad_i)/yirhatQ
  zhatr = 1 - zhati
!
!  sir = 2*phat(rad_i)%p*phat(rad_r)%p
  sir = sirhat_arr(rad_i,rad_r)
!
  if (p(radi)%flv.ne.0) then
    Ti2 = qcd_cf
  else
    Ti2 = qcd_ca
  end if
!
  StCirtSrtTerm = &
                + qcd_ca*2*zhati/(sir*zhatr)*( &
                +  s_i_r/(s_i_t*s_r_t) &
                +  zr/(s_r_t*zt) &
                -  zi/(s_i_t*zt)) &
                + Ti2*4*zi*zhati/(s_i_t*sir*zt*zhatr)
  StCirtSrtTerm = StCirtSrtTerm * Ti2*smeB
!
  StCirtSrtTerm = (8d0*pi)**2*StCirtSrtTerm
!
end subroutine CalcStCirtSrtFF
!
subroutine CalcStCSirtSrtFF(p,phat,ptilde,Bij,                  &
                            emit_ir,rad_i,rad_r,radi,radr,radt, &
                            StCSirtSrtTerm)
use particles
use regions
use math
use utils
use QCDparams
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p,phat,ptilde
  real(kind(1d0)) , dimension(:,:) , intent(in) :: Bij
  integer , intent(in) :: emit_ir,rad_i,rad_r,radi,radr,radt
  real(kind(1d0)) , intent(out) :: StCSirtSrtTerm
!
  integer :: n,j,l,s,emit,rad,jpart,lpart
  real(kind(1d0)) :: Ti2
  real(kind(1d0)) :: yirhatQ
!
  real(kind(1d0)) :: Sjl
!
  Sjl(j,l,s) = 2d0*sir_arr(j,l)/(sir_arr(j,s)*sir_arr(l,s))
!
  StCSirtSrtTerm = 0d0
!
  n = size(p)
! The position of the radiated parton originating from the 
! collinear pair:
  rad  = max(radi,radr)
  emit = min(radi,radr)
!
! Quantities defined with the help of the hatted momenta:
  yirhatQ = yihatQ_arr(rad_i) + yihatQ_arr(rad_r)
!
  zhati = yihatQ_arr(rad_i)/yirhatQ
  zhatr = 1 - zhati
!
!  sir = 2*phat(rad_i)%p*phat(rad_r)%p
  sir = sirhat_arr(rad_i,rad_r)
!
  if (p(radi)%flv.ne.0) then
    Ti2 = qcd_cf
  else
    Ti2 = qcd_ca
  end if
!
  jpart = 0
  do j=1,n
    if ((j.eq.radt).or.(j.eq.rad)) cycle
    jpart = jpart + 1
    if (j.eq.emit) cycle
    if (abs(p(j)%flv).gt.6) cycle
    lpart = jpart
    do l=j+1,n
      if ((l.eq.radt).or.(l.eq.rad)) cycle
      lpart = lpart + 1
      if (l.eq.emit) cycle
      if (abs(p(l)%flv).gt.6) cycle
      StCSirtSrtTerm = StCSirtSrtTerm &
                     + Sjl(j,l,radt)*Bij(jpart,lpart)
    end do
  end do
! Contribution is also coming from the (ir) leg, but for this
! particular one the form of \mathcal{S} changes,
! for the precise definition take a look at Eq. (6.26) of 
! arXiv:hep-ph/0609042.
  lpart = 0
  do l=1,n
    if ((l.eq.radt).or.(l.eq.rad)) cycle
    lpart = lpart + 1
    if (l.eq.emit) cycle
    if (abs(p(l)%flv).gt.6) cycle
    StCSirtSrtTerm = StCSirtSrtTerm &
                   + 2*(sir_arr(radi,l) + sir_arr(radr,l)) &
                   / ((sir_arr(radi,radt) + sir_arr(radr,radt))*sir_arr(l,radt)) &
                   * Bij(emit_ir,lpart)
  end do
!
  StCSirtSrtTerm = -(8d0*pi)**2*2*zhati/(sir*zhatr)*Ti2*StCSirtSrtTerm
!
end subroutine CalcStCSirtSrtFF
!
subroutine CalcStCirtCSirtSrtFF(p,phat,ptilde,smeB,         &
                                rad_i,rad_r,radi,radr,radt, &
                                StCirtCSirtSrtTerm)
use particles
use regions
use math
use utils
use QCDparams
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p,phat,ptilde
  real(kind(1d0)) , intent(in) :: smeB
  integer , intent(in) :: rad_i,rad_r,radi,radr,radt
  real(kind(1d0)) , intent(out) :: StCirtCSirtSrtTerm
!
  real(kind(1d0)) :: Ti2
  real(kind(1d0)) :: yirtQ,yirhatQ
  real(kind(1d0)) :: zt,s_ir_t
!
!
  StCirtCSirtSrtTerm = 0d0
!
! Quantities defined with the help of the original momenta:
  yirtQ = yiQ_arr(radi) + yiQ_arr(radr) + yiQ_arr(radt)
!
  zt = yiQ_arr(radt)/yirtQ
!
  s_ir_t = 2*(p(radi)%p + p(radr)%p)*p(radt)%p
!
! Quantities defined with the help of the hatted momenta:
  yirhatQ = yihatQ_arr(rad_i) + yihatQ_arr(rad_r)
!
  zhati = yihatQ_arr(rad_i)/yirhatQ
  zhatr = 1 - zhati
!
  sir = 2*phat(rad_i)%p*phat(rad_r)%p
!
  if (p(radi)%flv.ne.0) then
    Ti2 = qcd_cf
  else
    Ti2 = qcd_ca
  end if
!
  StCirtCSirtSrtTerm = 4*Ti2**2*zhati*(1 - zt)/(sir*s_ir_t*zhatr*zt) &
                     * smeB
!
  StCirtCSirtSrtTerm = (8d0*pi)**2*StCirtCSirtSrtTerm
!
end subroutine CalcStCirtCSirtSrtFF
!
subroutine CalcStSrtFF(p,phat,ptilde,Bij,Bijkl, &
                       rad_r,radr,radt,         &
                       StSrtTerm)
use particles
use regions
use math
use utils
use QCDparams
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p,phat,ptilde
  real(kind(1d0)) , dimension(:,:) , intent(in) :: Bij
  real(kind(1d0)) , dimension(:,:,:,:) , intent(in) :: Bijkl
  integer , intent(in) :: rad_r,radr,radt
  real(kind(1d0)) , intent(out) :: StSrtTerm
!
  integer :: i,j,k,l,r,t
  integer :: n
  integer :: ihat,khat
  integer :: itilde,jtilde,ktilde,ltilde
!
  real(kind(1d0)) :: Sjl,Sik
!
! This eikonal is defined with the original indices:
  Sjl(j,l,t) = 2*sir_arr(j,l)/(sir_arr(j,t)*sir_arr(l,t))
!
! This eikonal factor is defined on the hatted momenta:
  Sik(i,k,r) = 2*sirhat_arr(i,k)/(sirhat_arr(i,r)*sirhat_arr(k,r))
!
  StSrtTerm = 0d0
!
  n = size(p)
!
  itilde = 0
! i is defined on hat indices
  do i=1,n-1
! skipping parton \hat{r}
    if (i.eq.rad_r) cycle
    itilde = itilde + 1
! omitting non-QCD particles:
    if (abs(phat(i)%flv).gt.6) cycle
    ktilde = itilde
! k is defined on hat indices:
    do k=i+1,n-1
      if (k.eq.rad_r) cycle
      ktilde = ktilde + 1
      if (abs(phat(k)%flv).gt.6) cycle
      jtilde = 0
! j is defined on the original indices:
      do j=1,n
        if (j.eq.radr) cycle
        if (j.eq.radt) cycle
        jtilde = jtilde + 1
        if (abs(p(j)%flv).gt.6) cycle
        ltilde = jtilde
! l is defined on the original indices:
        do l=j+1,n
          if (l.eq.radr) cycle
          if (l.eq.radt) cycle
          ltilde = ltilde + 1
          if (abs(p(l)%flv).gt.6) cycle
          StSrtTerm = StSrtTerm &
                    + Sik(i,k,rad_r) &
                    * Sjl(j,l,radt) &
                    * Bijkl(itilde,ktilde,jtilde,ltilde)
        end do
      end do
    end do
  end do
!
! ihat is defined on indices with hat on them:
  ihat   = 0
! itilde is defined on indices with tilde on them:
  itilde = 0
! i is defined on the original indices:
  do i=1,n
    if (i.eq.radt) cycle
! By dropping parton t \hat{i} can be obtained:
    ihat = ihat + 1
! By dropping parton r \tilde{i} can be obtained:
    if (i.eq.radr) cycle
    itilde = itilde + 1
    if (abs(p(i)%flv).gt.6) cycle
    khat   = ihat
    ktilde = itilde
    do k=i+1,n
      if (k.eq.radt) cycle
      khat = khat + 1
      if (k.eq.radr) cycle
      ktilde = ktilde + 1
      if (abs(p(k)%flv).gt.6) cycle
      StSrtTerm = StSrtTerm &
                - qcd_ca*Sik(ihat,khat,rad_r) &
                * (Sjl(i,radr,radt) + Sjl(k,radr,radt) - Sjl(i,k,radt)) &
                * Bij(itilde,ktilde)
    end do
  end do
!
! The extra 1/2 is the remnant coming from the 1/8 and 1/4
! sitting at the front of the eikonals a 1/4 and a 1/2
! is missing provided the summation is ordered.
  StSrtTerm = 0.5d0*(8d0*pi)**2*StSrtTerm
!
end subroutine CalcStSrtFF
!
subroutine CalcCitStCirtFF(p,phat,ptilde,smeB,Bmunu,   &
                           rad_i,rad_r,radi,radr,radt, &
                           CitStCirtTerm)
use particles
use regions
use math
use utils
use QCDparams
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p,phat,ptilde
  real(kind(1d0)) , intent(in) :: smeB
  real(kind(1d0)) , dimension(0:3,0:3) , intent(in) :: Bmunu
  integer , intent(in) :: rad_i,rad_r,radi,radr,radt
  real(kind(1d0)) , intent(out) :: CitStCirtTerm
!
  real(kind(1d0)) :: Tit2
  real(kind(1d0)) :: zi,s_i_t
!
  real(kind(1d0)) , external :: Pqg0,Pqq0,Pgg0
!
  CitStCirtTerm = 0d0
!
  zi = yiQ_arr(radi)/(yiQ_arr(radi) + yiQ_arr(radt))
  s_i_t = sir_arr(radi,radt)
!
  if ((p(radi)%flv+p(radt)%flv).ne.0) then
    Tit2 = qcd_cf
  else
    Tit2 = qcd_ca
  end if
!
! q(~) -> q(~) + g
  if ((p(radi)%flv.ne.0).and. &
      (p(radr)%flv.eq.0)) then
!
    CitStCirtTerm = 2*zi/(1 - zi)*Tit2 &
                  * Pqg0(zhati,zhatr)*smeB
!
! q(~) -> g + q(~)
  elseif ((p(radi)%flv.eq.0).and. &
          (p(radr)%flv.ne.0)) then
!
    CitStCirtTerm = 2*zi/(1 - zi)*Tit2 &
                  * Pqg0(zhatr,zhati)*smeB
!
! g -> q + q~
  elseif  ((p(radi)%flv.ne.0).and. &
           (p(radi)%flv.eq.-p(radr)%flv)) then
!
    CitStCirtTerm = 2*zi/(1 - zi)*Tit2 &
                  * Pqq0(zhati,zhatr,kthati,Bmunu)
!
! g -> g + g
  elseif ((p(radi)%flv.eq.0).and. &
          (p(radr)%flv.eq.0)) then
!
    CitStCirtTerm = 2*zi/(1 - zi)*Tit2 &
                  * Pgg0(zhati,zhatr,kthati,Bmunu)
!
  else
    print *,"Error in CalcCitStCirt..."
    print *,"rad_i,rad_r,radi,radr,radt: ",rad_i,rad_r,radi,radr,radt
    call PrintSubProc(p)
    stop "CalcCitStCirt"
  end if
!
  CitStCirtTerm = (8d0*pi)**2/(s_i_t*sir)*CitStCirtTerm
!
end subroutine CalcCitStCirtFF
!
! This is just the previous routine but with an interchange of
! i and r:
subroutine CalcCrtStCirtFF(p,phat,ptilde,smeB,Bmunu,   &
                           rad_i,rad_r,radi,radr,radt, &
                           CrtStCirtTerm)
use particles
use regions
use math
use utils
use QCDparams
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p,phat,ptilde
  real(kind(1d0)) , intent(in) :: smeB
  real(kind(1d0)) , dimension(0:3,0:3) , intent(in) :: Bmunu
  integer , intent(in) :: rad_i,rad_r,radi,radr,radt
  real(kind(1d0)) , intent(out) :: CrtStCirtTerm
!
  real(kind(1d0)) :: Trt2
  real(kind(1d0)) :: zr,s_r_t
!
  real(kind(1d0)) , external :: Pqg0,Pqq0,Pgg0
!
  CrtStCirtTerm = 0d0
!
  zr = yiQ_arr(radr)/(yiQ_arr(radr) + yiQ_arr(radt))
  s_r_t = sir_arr(radr,radt)
!
  if ((p(radr)%flv+p(radt)%flv).ne.0) then
    Trt2 = qcd_cf
  else
    Trt2 = qcd_ca
  end if
!
! q(~) -> q(~) + g
  if ((p(radr)%flv.ne.0).and. &
      (p(radi)%flv.eq.0)) then
!
    CrtStCirtTerm = 2*zr/(1 - zr)*Trt2 &
                  * Pqg0(zhatr,zhati)*smeB
!
! q(~) -> g + q(~)
  elseif ((p(radr)%flv.eq.0).and. &
          (p(radi)%flv.ne.0)) then
!
    CrtStCirtTerm = 2*zr/(1 - zr)*Trt2 &
                  * Pqg0(zhati,zhatr)*smeB
!
! g -> q + q~
  elseif  ((p(radr)%flv.ne.0).and. &
           (p(radr)%flv.eq.-p(radi)%flv)) then
!
    CrtStCirtTerm = 2*zr/(1 - zr)*Trt2 &
                  * Pqq0(zhatr,zhati,kthatr,Bmunu)
!
! g -> g + g
  elseif ((p(radr)%flv.eq.0).and. &
          (p(radi)%flv.eq.0)) then
!
    CrtStCirtTerm = 2*zr/(1 - zr)*Trt2 &
                  * Pgg0(zhatr,zhati,kthatr,Bmunu)
!
  else
    print *,"Error in CalcCrtStCirt..."
    print *,"rad_i,rad_r,radi,radr,radt: ",rad_i,rad_r,radi,radr,radt
    call PrintSubProc(p)
    stop "CalcCrtStCirt"
  end if
!
  CrtStCirtTerm = (8d0*pi)**2/(s_r_t*sir)*CrtStCirtTerm
!
end subroutine CalcCrtStCirtFF
!
subroutine CalcCktStCSirtFF(p,phat,ptilde,smeB,Bmunu, &
                            radi,radr,radk,radt,      &
                            CktStCSirtTerm)
use particles
use regions
use math
use utils
use QCDparams
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p,phat,ptilde
  real(kind(1d0)) , intent(in) :: smeB
  real(kind(1d0)) , dimension(0:3,0:3) , intent(in) :: Bmunu
  integer , intent(in) :: radi,radr,radk,radt
  real(kind(1d0)) , intent(out) :: CktStCSirtTerm
!
  real(kind(1d0)) :: Tkt2
  real(kind(1d0)) :: zk
!
  real(kind(1d0)) , external :: Pqg0,Pqq0,Pgg0
!
  CktStCSirtTerm = 0d0
!
  zk = yiQ_arr(radk)/(yiQ_arr(radk) + yiQ_arr(radt))
  s_k_t = sir_arr(radk,radt)
!
  if ((p(radk)%flv+p(radt)%flv).ne.0) then
    Tkt2 = qcd_cf
  else
    Tkt2 = qcd_ca
  end if
!
! q(~) -> q(~) + g
  if ((p(radi)%flv.ne.0).and. &
      (p(radr)%flv.eq.0)) then
!
    CktStCSirtTerm = 2*zk/(1 - zk)*Tkt2 &
                   * Pqg0(zhati,zhatr)*smeB
!
! g -> q + q~
  elseif  ((p(radi)%flv.ne.0).and. &
           (p(radi)%flv.eq.-p(radr)%flv)) then
!
    CktStCSirtTerm = 2*zk/(1 - zk)*Tkt2 &
                   * Pqq0(zhati,zhatr,kthati,Bmunu)
!
! g -> g + g
  elseif ((p(radi)%flv.eq.0).and. &
          (p(radr)%flv.eq.0)) then
!
    CktStCSirtTerm = 2*zk/(1 - zk)*Tkt2 &
                   * Pgg0(zhati,zhatr,kthati,Bmunu)
!
  else
    print *,"Error in CalcCktStCSirt..."
    print *,"radi,radr,radk,radt: ",radi,radr,radk,radt
    call PrintSubProc(p)
    stop "CalcCktStCSirt"
  end if
!
  CktStCSirtTerm = (8d0*pi)**2/(s_k_t*sir)*CktStCSirtTerm
!
end subroutine CalcCktStCSirtFF
!
subroutine CalcCktStCSirtSrtFF(p,phat,ptilde,smeB,  &
                               rad_i,rad_r,         &
                               radi,radr,radk,radt, &
                               CktStCSirtSrtTerm)
use particles
use regions
use math
use utils
use QCDparams
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p,phat,ptilde
  real(kind(1d0)) , intent(in) :: smeB
  integer , intent(in) :: rad_i,rad_r,radi,radr,radk,radt
  real(kind(1d0)) , intent(out) :: CktStCSirtSrtTerm
!
  real(kind(1d0)) :: Tir2,Tk2
  real(kind(1d0)) :: zi,zk
!
!
  zk = yiQ_arr(radk)/(yiQ_arr(radk) + yiQ_arr(radt))
  s_k_t = sir_arr(radk,radt)
!
  zi = yihatQ_arr(rad_i)/(yihatQ_arr(rad_i) + yihatQ_arr(rad_r))
  sir = sirhat_arr(rad_i,rad_r)
!
  if ((p(radi)%flv+p(radr)%flv).ne.0) then
    Tir2 = qcd_cf
  else
    Tir2 = qcd_ca
  end if
  if (p(radk)%flv.ne.0) then
    Tk2 = qcd_cf
  else
    Tk2 = qcd_ca
  end if
!
  CktStCSirtSrtTerm = (8d0*pi)**2 &
                    * 4*zi*zk/(s_k_t*sir*(1 - zi)*(1 - zk)) &
                    * Tir2*Tk2*smeB
!
end subroutine CalcCktStCSirtSrtFF
!
subroutine CalcCktStCkrtSrtFF(p,phat,ptilde,smeB, &
                               rad_k,rad_r,       &
                               radk,radr,radt,    &
                               CktStCkrtSrtTerm)
use particles
use regions
use math
use utils
use QCDparams
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p,phat,ptilde
  real(kind(1d0)) , intent(in) :: smeB
  integer , intent(in) :: rad_k,rad_r,radk,radr,radt
  real(kind(1d0)) , intent(out) :: CktStCkrtSrtTerm
!
  real(kind(1d0)) :: Tkr2,Tk2
  real(kind(1d0)) :: skr,z_k_r,z_k_t
!
!
  z_k_t = yiQ_arr(radk)/(yiQ_arr(radk) + yiQ_arr(radt))
  s_k_t = sir_arr(radk,radt)
!
  z_k_r = yihatQ_arr(rad_k)/(yihatQ_arr(rad_k) + yihatQ_arr(rad_r))
  skr = sirhat_arr(rad_k,rad_r)
!
  if ((p(radk)%flv+p(radr)%flv).ne.0) then
    Tkr2 = qcd_cf
  else
    Tkr2 = qcd_ca
  end if
  if (p(radk)%flv.ne.0) then
    Tk2 = qcd_cf
  else
    Tk2 = qcd_ca
  end if
!
  CktStCkrtSrtTerm = (8d0*pi)**2 &
                   * 4*z_k_r*z_k_t/(s_k_t*skr*(1 - z_k_r)*(1 - z_k_t)) &
                   * Tkr2*Tk2*smeB
!
end subroutine CalcCktStCkrtSrtFF
!
subroutine CalcCrtStCkrtSrtFF(p,phat,ptilde,smeB, &
                               rad_k,rad_r,       &
                               radk,radr,radt,    &
                               CrtStCkrtSrtTerm)
use particles
use regions
use math
use utils
use QCDparams
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p,phat,ptilde
  real(kind(1d0)) , intent(in) :: smeB
  integer , intent(in) :: rad_k,rad_r,radk,radr,radt
  real(kind(1d0)) , intent(out) :: CrtStCkrtSrtTerm
!
  real(kind(1d0)) :: Tkr2
  real(kind(1d0)) :: skr,srt,z_k_r,z_r_t
!
!
  z_r_t = yiQ_arr(radr)/(yiQ_arr(radr) + yiQ_arr(radt))
  srt = sir_arr(radr,radt)
!
  z_k_r = yihatQ_arr(rad_k)/(yihatQ_arr(rad_k) + yihatQ_arr(rad_r))
  skr = sirhat_arr(rad_k,rad_r)
!
  if ((p(radk)%flv+p(radr)%flv).ne.0) then
    Tkr2 = qcd_cf
  else
    Tkr2 = qcd_ca
  end if
!
  CrtStCkrtSrtTerm = (8d0*pi)**2 &
                   * 4*z_k_r*z_r_t/(srt*skr*(1 - z_k_r)*(1 - z_r_t)) &
                   * Tkr2*qcd_ca*smeB
!
end subroutine CalcCrtStCkrtSrtFF
!
subroutine CalcCktStSrtFF(p,phat,ptilde,Bij,    &
                          rad_r,radk,radr,radt, &
                          CktStSrtTerm)
use particles
use regions
use math
use utils
use QCDparams
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p,phat,ptilde
  real(kind(1d0)) , dimension(:,:) , intent(in) :: Bij
  integer , intent(in) :: rad_r,radk,radr,radt
  real(kind(1d0)) , intent(out) :: CktStSrtTerm
!
  integer :: j,l,r,jpart,lpart,n
  real(kind(1d0)) :: Tk2
  real(kind(1d0)) :: z_k_t
!
  real(kind(1d0)) :: Sjl
!
! This eikonal factor is defined on the hatted momenta:
  Sjl(j,l,r) = 2d0*sirhat_arr(j,l)/(sirhat_arr(j,r)*sirhat_arr(l,r))
!
  z_k_t = yiQ_arr(radk)/(yiQ_arr(radk) + yiQ_arr(radt))
  s_k_t = sir_arr(radk,radt)
!
  if (p(radk)%flv.ne.0) then
    Tk2 = qcd_cf
  else
    Tk2 = qcd_ca
  end if
!
  CktStSrtTerm = 0d0
!
! parton index j,l and r are defined on momenta with hat:
  n = size(p)
  r = rad_r
! jpart and lpart are defined on underlying Born momenta:
  jpart = 0
  do j=1,n-1
    if (j.eq.r) cycle
    jpart = jpart + 1
    if (abs(phat(j)%flv).gt.6) cycle
    lpart = jpart
    do l=j+1,n-1
      if (l.eq.r) cycle
      lpart = lpart + 1
      if (abs(phat(l)%flv).gt.6) cycle
      CktStSrtTerm = CktStSrtTerm &
                   + Sjl(j,l,r)*Bij(jpart,lpart)   
    end do
  end do
!
  CktStSrtTerm = -(8d0*pi)**2 &
               * 2*z_k_t/(s_k_t*(1 - z_k_t)) &
               * Tk2 &
               * CktStSrtTerm
!
end subroutine CalcCktStSrtFF
!
subroutine CalcCrtStSrtFF(p,phat,ptilde,Bij, &
                          rad_r,radr,radt,   &
                          CrtStSrtTerm)
use particles
use regions
use math
use utils
use QCDparams
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p,phat,ptilde
  real(kind(1d0)) , dimension(:,:) , intent(in) :: Bij
  integer , intent(in) :: rad_r,radr,radt
  real(kind(1d0)) , intent(out) :: CrtStSrtTerm
!
  integer :: j,l,r,jpart,lpart,n
  real(kind(1d0)) :: s_r_t,z_r_t
!
  real(kind(1d0)) :: Sjl
!
! This eikonal factor is defined on the hatted momenta:
  Sjl(j,l,r) = 2d0*sirhat_arr(j,l)/(sirhat_arr(j,r)*sirhat_arr(l,r))
!
  z_r_t = yiQ_arr(radr)/(yiQ_arr(radr) + yiQ_arr(radt))
  s_r_t = sir_arr(radr,radt)
!
  CrtStSrtTerm = 0d0
!
! parton index j,l and r are defined on momenta with hat:
  n = size(p)
  r = rad_r
! jpart and lpart are defined on underlying Born momenta:
  jpart = 0
  do j=1,n-1
    if (j.eq.r) cycle
    jpart = jpart + 1
    if (abs(phat(j)%flv).gt.6) cycle
    lpart = jpart
    do l=j+1,n-1
      if (l.eq.r) cycle
      lpart = lpart + 1
      if (abs(phat(l)%flv).gt.6) cycle
      CrtStSrtTerm = CrtStSrtTerm &
                   + Sjl(j,l,r)*Bij(jpart,lpart)   
    end do
  end do
!
  CrtStSrtTerm = -(8d0*pi)**2 &
               * 2*z_r_t/(s_r_t*(1 - z_r_t)) &
               * qcd_ca &
               * CrtStSrtTerm
!
end subroutine CalcCrtStSrtFF
!
end module nnlo_subtractions
!
module nlo_mom_maps
implicit none
!
contains
!
! This routine remaps momenta, taking them from p and placing
! them in ptilde multiplying with a prefact factor the emit-th
! momentum is taken from pp, no prefactor is applied for it
! and momentum at position rad1 is omitted:
subroutine RemapMomentaO1(p,ptilde,prefact,emit,pp,rad1)
use particles
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
  type(particle) , dimension(:) , intent(out) :: ptilde
  real(kind(1d0)) , intent(in) :: prefact
  integer , intent(in) :: emit
  type(particle) , intent(in) :: pp
  integer , intent(in) :: rad1
!
  integer :: ipart,jpart
  integer :: npart
!
!
! The number of particles is taken from the original set
! of momenta:
  npart = size(p)
! Position in the original momentum array:
  ipart = 1
! Position in the new momentum array:
  jpart = 1
  do while (.true.)
! We have npart - 1 momenta in the ptilde array, do loop
! should be quit if we pass this number:
    if (jpart.gt.npart-1) exit
! We have to omit the momentum which corresponds to the
! emitted particle:
    if (ipart.eq.rad1) then
      ipart = ipart + 1
      cycle
    end if
! If dealing with a momentum different from the emitter
! we include a prefactor if situating in the final state:
    if (ipart.ne.emit) then
      if (ipart.gt.2) ptilde(jpart) = prefact*p(ipart)
      if (ipart.lt.3) ptilde(jpart) = p(ipart)
! For the emitter we substitute pp:
    else
      ptilde(jpart) = pp
    end if
    ipart = ipart + 1
    jpart = jpart + 1
  end do
!
end subroutine RemapMomentaO1
!
subroutine MapMomCirFF(p,ptilde,emitir,emitID,radi,radiID,radr,radrID)
use particles
use regions
use math
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
  type(particle) , dimension(:) , intent(out) :: ptilde
  integer , intent(in) :: emitir,radi,radr
  integer , intent(in) :: emitID,radiID,radrID
!
  integer :: ipart,emit,rad
!
!
!
!
! Obtaining y variables:
  yiQ  = yiQ_arr(radi)
  yrQ  = yiQ_arr(radr)
  yirQ = yiQ + yrQ
  sir  = sir_arr(radi,radr)
  yir  = sir/Q2
  ztildei = yiQ/yirQ
  ztilder = 1d0 - ztildei
  alphair = 0.5d0*(yirQ - sqrt(abs(yirQ**2 - 4d0*yir)))
! The new momentum:
  ptildeir%p = (p(radi)%p + p(radr)%p - alphair*Q)/(1d0 - alphair)
!
! We calculate the transverse momentum, too:
  zetai_r    = ztildei - yir/(alphair*yirQ)
  zetar_i    = ztilder - yir/(alphair*yirQ)
  yirtildeQ  = 2d0*ptildeir%p*Q/Q2
  zetair     = yir/(alphair*yirtildeQ)*(ztilder - ztildei)
  kttildei%p = zetai_r*p(radr)%p - zetar_i*p(radi)%p &
             + 0*zetair*ptildeir%p
  kttildei%flv = emitID
  if (ptildeir%p%E.lt.0) then
    print *,"ptildeir: "
    call PrintMom(ptildeir%p)
    print *,"sir: ",sir
    print *,"yir: ",yir
    print *,"yiQ,yrQ: ",yiQ,yrQ
    print *,"yirQ: ",yirQ
    print *,"ztildei,ztilder: ",ztildei,ztilder
    print *,"alphair: ",alphair
    print *,"zetai_r,zetar_r,zetair: ",zetai_r,zetar_i,zetair
    call PrintParts(kttildei)
  end if
! Safety check upon branching:
  if ((radiID+radrID).ne.emitID) then
    print *,"Problem in MapMomCirFF..."
    write(*,'(a,3(1x,I0))') "emitir,radi,radr: ",emitir,radi,radr
    write(*,'(a,3(1x,I0))') "emitID,radiID,radrID: ",emitID,radiID,radrID
    stop "MapMomCirFF"
  end if
  ptildeir%flv = p(radi)%flv + p(radr)%flv
! it is possible that radi > radr, but we want to remove the parton
! with higher position always, hence
  emit = min(radi,radr)
  rad  = max(radi,radr)
  call RemapMomentaO1(p,ptilde,1d0/(1d0 - alphair),emit,ptildeir,rad)
!
!
end subroutine MapMomCirFF
!
subroutine MapMomSr(p,ptilde,radr,radrID)
use particles
use regions
use misc
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
  type(particle) , dimension(:) , intent(out) :: ptilde
  integer , intent(in) :: radr
  integer , intent(in) :: radrID
!
!
  integer :: ipart,jpart
!
  real(kind(1d0)) :: lambdar
!
! Obtaining y variables:
  yrQ  = yiQ_arr(radr)
  lambdar = sqrt(1d0 - yrQ)
! Safety check, the radiated parton should be a gluon:
  if (radrID.ne.0) then
    print *,"Problem in MapmomSr..."
    write(*,'(a,3(1x,I0))') "radr: ",radr
    write(*,'(a,3(1x,I0))') "radrID: ",radr
  end if
! The soft-type momentum mapping is fairly simple to omit the
! declaration of a separate routine just to deal with the mapping
! which only involves a unique Lorentz transformation and a rescaling:
  jpart = 1
! We run through all the momenta:
  do ipart=1,size(p)
! We skip the gluon:
    if (ipart.eq.radr) cycle
! For the initial state we simply copy:
    if (ipart.lt.3) ptilde(jpart) = p(ipart)
! Otherwise perform a Lorentz boost:
    if (ipart.gt.2) call LorentzLambda(Q,(Q - p(radr)%p)/lambdar, &
                                       p(ipart)/lambdar,ptilde(jpart))
    jpart = jpart + 1
  end do
!
!  print *,"Soft momentum mapping: "
!  call PrintParts(p)
!  call PrintParts(ptilde)
!
!
end subroutine MapMomSr
!
end module nlo_mom_maps
!
module nnlo_mom_maps
implicit none
!
interface RemapMomentaO2
  module procedureRemapMomentaO2_dou
  module procedureRemapMomentaO2_tri
end interface RemapMomentaO2
!
contains
! This routine is the same as RemapMomentaO1 but it
! omits two momenta (radi and radj) from the ptilde array:
! For the curious one O2 means: omit 2...
subroutine RemapMomentaO2_tri(p,ptilde,prefact,emit,pp,radi,radj)
use particles
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
  type(particle) , dimension(:) , intent(out) :: ptilde
  real(kind(1d0)) , intent(in) :: prefact
  integer , intent(in) :: emit
  type(particle) , intent(in) :: pp
  integer , intent(in) :: radi,radj
!
  integer :: ipart,jpart
  integer :: npart
!
!
! The number of particles is taken from the original set
! of momenta:
  npart = size(p)
! Position in the original momentum array:
  ipart = 1
! Position in the new momentum array:
  jpart = 1
  do while (.true.)
! We have npart - 2 momenta in the ptilde array, do loop
! should quit if we pass this number:
    if (jpart.gt.npart-2) exit
! We have to omit the momentum which corresponds to the
! emitted particles:
    if ((ipart.eq.radi).or.(ipart.eq.radj)) then
      ipart = ipart + 1
      cycle
    end if
! If dealing with a momentum different from the emitter
! we include a prefactor if situating in the final state:
    if (ipart.ne.emit) then
      if (ipart.gt.2) ptilde(jpart) = prefact*p(ipart)
      if (ipart.lt.3) ptilde(jpart) = p(ipart)
! For the emitter we substitute pp:
    else
      ptilde(jpart) = pp
    end if
    ipart = ipart + 1
    jpart = jpart + 1
  end do
!
end subroutine RemapMomentaO2_tri
!
! This routine is the same as RemapMomentaO1 but it
! omits two momenta (radi and radj) from the ptilde array,
! while replacing momenta at positions emit1 and emit2 with
! pp1 and pp2:
subroutine RemapMomentaO2_dou(p,ptilde,prefact, &
                              emit1,pp1,emit2,pp2,radi,radj)
use particles
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
  type(particle) , dimension(:) , intent(out) :: ptilde
  real(kind(1d0)) , intent(in) :: prefact
  integer , intent(in) :: emit1,emit2
  type(particle) , intent(in) :: pp1,pp2
  integer , intent(in) :: radi,radj
!
  integer :: ipart,jpart
  integer :: npart
!
!
! The number of particles is taken from the original set
! of momenta:
  npart = size(p)
! Position in the original momentum array:
  ipart = 1
! Position in the new momentum array:
  jpart = 1
  do while (.true.)
! We have npart - 2 momenta in the ptilde array, do loop
! should quit if we pass this number:
    if (jpart.gt.npart-2) exit
! We have to omit the momentum which corresponds to the
! emitted particles:
    if ((ipart.eq.radi).or.(ipart.eq.radj)) then
      ipart = ipart + 1
      cycle
    end if
! If the position coincides with the position of the first
! emitter substitute momentum pp1, instead:
    if (ipart.eq.emit1) then
      ptilde(jpart) = pp1
! If the position coincides with the position of the second
! emitter substitute momentum pp2, instead:
    elseif (ipart.eq.emit2) then
      ptilde(jpart) = pp2
! Otherwise substitute the momentum can be found in the original
! array multiplied by a prefactor if it is a final state particle:
    else
      if (ipart.gt.2) ptilde(jpart) = prefact*p(ipart)
      if (ipart.lt.3) ptilde(jpart) = p(ipart)
    end if
    ipart = ipart + 1
    jpart = jpart + 1
  end do
!
end subroutine RemapMomentaO2_dou
!
! This routine can be used to bring the gluon legs in order:
subroutine ReOrderRads(flv,radi,radr,rads,rad1,rad2,rad3)
use utils
implicit none
!
  integer , dimension(:) , intent(in) :: flv
  integer , intent(in) :: radi,radr,rads
  integer , intent(out) :: rad1,rad2,rad3
!
!
!
!  print *,"ReOrderRads..."
!  print *,"radi,radr,rads: ",radi,radr,rads
!  print *,ConvertFromPDG(flv(radi))," , ", &
!          ConvertFromPDG(flv(radr))," , ", &
!          ConvertFromPDG(flv(rads))
!
  rad1 = radi
  rad2 = radr
  rad3 = rads
!
! q(~) g g:
  if ((flv(radi).ne.0).and.(flv(radr).eq.0)) then
!    print *,"q(~) g g"
    if (radr.gt.rads) then
      rad1 = radi
      rad2 = rads
      rad3 = radr
    end if
!
! q q~ g:
  elseif ((flv(radi).ne.0).and.(flv(radr).ne.0)) then
!    print *,"q q~ g"
    rad1 = rads
    rad2 = radi
    rad3 = radr
! 
! g g g:
  elseif ((flv(radi).eq.0).and.(flv(radr).eq.0)) then
!    print *,"g g g"
! 2 3 1
    if ((radr.lt.rads).and.(rads.lt.radi)) then
      rad1 = radr
      rad2 = rads
      rad3 = radi
! 3 1 2
    elseif ((rads.lt.radi).and.(radi.lt.radr)) then
      rad1 = rads
      rad2 = radi
      rad3 = radr
! 2 1 3
    elseif ((radr.lt.radi).and.(radi.lt.rads)) then
      rad1 = radr
      rad2 = radi
      rad3 = rads
! 1 3 2
    elseif ((radi.lt.rads).and.(rads.lt.radr)) then
      rad1 = radi
      rad2 = rads
      rad3 = radr
! 3 2 1
    elseif ((rads.lt.radr).and.(radr.lt.radi)) then
      rad1 = rads
      rad2 = radr
      rad3 = radi
    end if
  end if
!
!  print *,"rad1,rad2,rad3: ",rad1,rad2,rad3
!
!  read(*,*)
!
end subroutine ReOrderRads
! 
! This routine is also a carbon-copy of the CalcSubInvariants
! routine, this is mainly used in the case of iterated mappings
! , that is in the construction of A12 counterterms:
subroutine CalcSubhatInvariants(p)
use regions
use particles
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
!
  integer :: ipart,jpart
  integer :: npart
!
!
  npart = size(p)
!
  do ipart=1,npart
    yihatQ_arr(ipart) = 2*p(ipart)%p*Q/Q2
    do jpart=ipart+1,npart
      sirhat_arr(ipart,jpart) = 2*p(ipart)%p*p(jpart)%p
      sirhat_arr(jpart,ipart) = sirhat_arr(ipart,jpart)
    end do
  end do
!
end subroutine CalcSubhatInvariants
!
! This is a carbon-copy of the MapmomSr routine can be found
! in the previous module with one exception: this can be
! used in the NNLO calculation where iterative mappings are
! used while radr denotes the position of the soft particle
! in the original PS point rad_r denotes the position in the
! hated PS point, that is where a momentum mapping is already
! used:
subroutine MapMomSrNNLO(p,ptilde,radr,rad_r,radrID)
use particles
use regions
use misc
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
  type(particle) , dimension(:) , intent(out) :: ptilde
  integer , intent(in) :: radr,rad_r
  integer , intent(in) :: radrID
!
!
  integer :: ipart,jpart
!
  real(kind(1d0)) :: lambdar
!
! Obtaining y variables:
  yrQ  = 2*p(rad_r)%p*Q/Q2
  lambdar = sqrt(1d0 - yrQ)
! Safety check, the radiated parton should be a gluon:
  if (radrID.ne.0) then
    print *,"Problem in MapmomSrNNLO..."
    write(*,'(a,3(1x,I0))') "radr,rad_r: ",radr,rad_r
    write(*,'(a,3(1x,I0))') "radrID: ",radrID
  end if
! The soft-type momentum mapping is fairly simple to omit the
! declaration of a separate routine just to deal with the mapping
! which only involves a unique Lorentz transformation and a rescaling:
  jpart = 1
! We run through all the momenta:
  do ipart=1,size(p)
! We skip the gluon:
    if (ipart.eq.rad_r) cycle
! For the initial state we simply copy:
    if (ipart.lt.3) ptilde(jpart) = p(ipart)
! Otherwise perform a Lorentz boost:
    if (ipart.gt.2) call LorentzLambda(Q,(Q - p(rad_r)%p)/lambdar, &
                                       p(ipart)/lambdar,ptilde(jpart))
    jpart = jpart + 1
  end do
!
! Additionally the transverse momentum created in the previous
! mapping has to be boosted, too:
  call LorentzLambda(Q,(Q - p(rad_r)%p)/lambdar, &
                     kttildei/lambdar,kttildei)
!
end subroutine MapMomSrNNLO
!
! This routine is the same as the one used in the NLO case,
! but it takes all variables from the argument list hence 
! can be conveniently used for iterated mappings such as in
! the case of the A12 terms:
subroutine MapMomCirFF_A12(p,ptilde,radi,radr,Q,Q2,yiQ_arr,sir_arr, &
                           sir,zi,zr,kti,ktr,alphair)
use particles
use nlo_mom_maps
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
  type(particle) , dimension(:) , intent(out) :: ptilde
  integer , intent(in) :: radi,radr
  type(mom) , intent(in) :: Q
  real(kind(1d0)) , intent(in) :: Q2
  real(kind(1d0)) , dimension(:) , intent(in) :: yiQ_arr
  real(kind(1d0)) , dimension(:,:) , intent(in) :: sir_arr
  real(kind(1d0)) , intent(out) :: sir
  real(kind(1d0)) , intent(out) :: zi,zr
  type(particle) , intent(out) :: kti,ktr
  real(kind(1d0)) , intent(out) :: alphair
!
  integer :: emit,rad
  real(kind(1d0)) :: yiQ,yrQ,yirQ,yir,yirtildeQ
  real(kind(1d0)) :: zetai_r,zetar_i,zetair,zetari
  type(particle) :: pir
!
!
! Obtaining y variables:
  yiQ  = yiQ_arr(radi)
  yrQ  = yiQ_arr(radr)
  yirQ = yiQ + yrQ
  sir  = sir_arr(radi,radr)
  yir  = sir/Q2
! fractions:
  zi = yiQ/yirQ
  zr = 1 - zi
!
  alphair = 0.5d0*(yirQ - sqrt(abs(yirQ**2 - 4*yir)))
! The new momentum:
  pir%p = (p(radi)%p + p(radr)%p - alphair*Q)/(1d0 - alphair)
!
! Transverse momentum:
  zetai_r    = zi - yir/(alphair*yirQ)
  zetar_i    = zr - yir/(alphair*yirQ)
  yirtildeQ  = 2*pir%p*Q/Q2
  zetair     = yir/(alphair*yirtildeQ)*(zr - zi)
  zetari     = yir/(alphair*yirtildeQ)*(zi - zr)
  kti%p = zetai_r*p(radr)%p - zetar_i*p(radi)%p &
             + zetair*pir%p
  ktr%p = zetar_i*p(radi)%p - zetai_r*p(radr)%p &
             + zetari*pir%p
  pir%flv = p(radi)%flv + p(radr)%flv
  kti%flv = pir%flv
  ktr%flv = pir%flv
!
! Mapping n momenta to n-1:
! Note that in the present routine there is no implicit ordering
! defined for radi and radr, hence the emitter and the radiated
! parton is defined according to their position only:
  emit = min(radi,radr)
  rad  = max(radi,radr)
  call RemapMomentaO1(p,ptilde,1d0/(1d0 - alphair),emit,pir,rad)
!
end subroutine MapMomCirFF_A12
!
! This routine is the same as the one used in the NLO case,
! but it takes all variables from the argument list hence 
! can be conveniently used for iterated mappings such as in
! the case of the A12 terms:
subroutine MapMomSr_A12(p,ptilde,radr,Q,Q2,yiQ_arr,yrhatQ,lambdar)
use particles
use nlo_mom_maps
use misc
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
  type(particle) , dimension(:) , intent(out) :: ptilde
  integer , intent(in) :: radr
  type(mom) , intent(in) :: Q
  real(kind(1d0)) , intent(in) :: Q2
  real(kind(1d0)) , dimension(:) , intent(in) :: yiQ_arr
  real(kind(1d0)) , intent(out) :: yrhatQ,lambdar
!
  integer :: ipart,jpart
  type(mom) :: Qmpr
!
!
  yrhatQ  = yiQ_arr(radr)
  lambdar = sqrt(abs(1 - yrhatQ))
! The contribution is trivially zero:
  if (lambdar.eq.0) return
  Qmpr = Q - p(radr)%p
! Map momenta:
  jpart = 1
! We run through all the momenta:
  do ipart=1,size(p)
! We skip the gluon:
    if (ipart.eq.radr) cycle
! For the initial state we simply copy:
    if (ipart.lt.3) ptilde(jpart) = p(ipart)
! Otherwise perform a Lorentz boost:
    if (ipart.gt.2) call LorentzLambda(Q,Qmpr/lambdar, &
                                       p(ipart)/lambdar,ptilde(jpart))
    jpart = jpart + 1
  end do
!
end subroutine MapMomSr_A12
!
! Triply-collinear momentum mapping as in arXiv:0609042:
subroutine MapMomCirsFF(p,ptilde,emit,emitID,radi,radiID, &
                        radr,radrID,rads,radsID)
use regions
use particles
use misc
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
  type(particle) , dimension(:) , intent(out) :: ptilde
  integer , intent(in) :: emit,radi,radr,rads
  integer , intent(in) :: emitID,radiID,radrID,radsID
!
! The declarations in the next line correspond to y_{(ij)k}:
  real(kind(1d0)) :: y_rs_i,y_is_r,y_ir_s
  real(kind(1d0)) :: yis
  real(kind(1d0)) :: yirstildeQ
  real(kind(1d0)) :: zetai_rs,zetar_is,zetas_ir
  real(kind(1d0)) :: zetaris,zetairs,zetasir
!
!
! Obtaining y variables:
  yiQ  = yiQ_arr(radi)
  yrQ  = yiQ_arr(radr)
  ysQ  = yiQ_arr(rads)
  yirsQ = yiQ + yrQ + ysQ
  sirs  = sir_arr(radi,radr) + sir_arr(radi,rads) + sir_arr(radr,rads)
  yirs  = sirs/Q2
!
  yir = sir_arr(radi,radr)/Q2
  yis = sir_arr(radi,rads)/Q2
  yrs = sir_arr(radr,rads)/Q2
!
! y_ij_k := y_{(ij)k}:
  y_rs_i = yir + yis
  y_is_r = yir + yrs
  y_ir_s = yis + yrs
!
! momentum fractions (ztildei := z_{i,rs}):
  ztildei = yiQ/yirsQ
  ztilder = yrQ/yirsQ
  ztildes = ysQ/yirsQ
!
  alphairs = 0.5d0*(yirsQ - sqrt(abs(yirsQ**2 - 4d0*yirs)))
!
  ptildeirs%p = (p(radi)%p + p(radr)%p + p(rads)%p - alphairs*Q) &
              / (1d0 - alphairs)
  ptildeirs%flv = emitID
!
  yirstildeQ = ptildeirs%p*Q/Q2
!
  zetai_rs = ztildei - y_rs_i/(alphairs*yirsQ)
  zetar_is = ztilder - y_is_r/(alphairs*yirsQ)
  zetas_ir = ztildes - y_ir_s/(alphairs*yirsQ)
!
  zetaris = (yir + yrs - 2d0*ztilder*yirs) / (alphairs*yirstildeQ)
  zetairs = (yir + yis - 2d0*ztildei*yirs) / (alphairs*yirstildeQ)
  zetasir = (yis + yrs - 2d0*ztildes*yirs) / (alphairs*yirstildeQ)
!
! Note that the contributions proportional to zetaxxx are set to
! zero since in an NNLO calculation they are not needed.
  kttilder%p = zetar_is*p(radi)%p - zetai_rs*p(radr)%p &
             + zetar_is*p(rads)%p - zetas_ir*p(radr)%p &
             + 0*zetaris*ptildeirs%p
!
  kttildes%p = zetas_ir*p(radi)%p - zetai_rs*p(rads)%p &
             + zetas_ir*p(radr)%p - zetar_is*p(rads)%p &
             + 0*zetasir*ptildeirs%p
!
  kttildei%p = zetai_rs*p(radr)%p - zetar_is*p(radi)%p &
             + zetai_rs*p(rads)%p - zetas_ir*p(radi)%p &
             + 0*zetairs*ptildeirs%p
!
! Construction of the tilde momenta such a way that the
! momenta corresponding to the emitted partons are omitted:
  if (emit.eq.radi) then
    call RemapMomentaO2(p,ptilde,1d0/(1d0 - alphairs), &
                        emit,ptildeirs,radr,rads)
  elseif (emit.eq.radr) then
    call RemapMomentaO2(p,ptilde,1d0/(1d0 - alphairs), &
                        emit,ptildeirs,radi,rads)
  elseif (emit.eq.rads) then
    call RemapMomentaO2(p,ptilde,1d0/(1d0 - alphairs), &
                        emit,ptildeirs,radi,radr)
  end if
!
end subroutine MapMomCirsFF
!
! Double-collinear momentum mapping as in arXiv:0609042:
subroutine MapMomCirjsFF(p,ptilde, &
                         emiti,emitiID,radi,radiID,radr,radrID, &
                         emitj,emitjID,radj,radjID,rads,radsID)
use regions
use particles
use misc
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
  type(particle) , dimension(:) , intent(out) :: ptilde
  integer , intent(in) :: emiti,emitj,radi,radr,radj,rads
  integer , intent(in) :: emitiID,emitJID,radiID,radrID,radjID,radsID
!
!
!
! Obtaining y variables:
  yiQ  = yiQ_arr(radi)
  yrQ  = yiQ_arr(radr)
  yjQ  = yiQ_arr(radj)
  ysQ  = yiQ_arr(rads)
  yirQ = yiQ + yrQ
  yjsQ = yjQ + ysQ
  sir  = sir_arr(radi,radr)
  sjs  = sir_arr(radj,rads)
  yir  = sir/Q2
  yjs  = sjs/Q2
  ztildei = yiQ/yirQ
  ztilder = 1d0 - ztildei
  ztildej = yjQ/yjsQ
  ztildes = 1d0 - ztildej
  alphair = 0.5d0*(yirQ - sqrt(abs(yirQ**2 - 4d0*yir)))
  alphajs = 0.5d0*(yjsQ - sqrt(abs(yjsQ**2 - 4d0*yjs)))
!
! The new momenta:
  ptildeir%p = (p(radi)%p + p(radr)%p - alphair*Q)/(1d0 - alphair - alphajs)
  ptildejs%p = (p(radj)%p + p(rads)%p - alphajs*Q)/(1d0 - alphair - alphajs)
!
! We calculate the transverse momenta, too:
! i-r pair:
  zetai_r    = ztildei - yir/(alphair*yirQ)
  zetar_i    = ztilder - yir/(alphair*yirQ)
!  zetair     = yir/(alphair*yirtildeQ)*(ztilder - ztildei)
! \zeta_{ir} can be chosen to be zero in NNLO:
  zetair     = 0
  kttildei%p = zetai_r*p(radr)%p - zetar_i*p(radi)%p &
             + zetair*ptildeir%p
  kttildei%flv = emitiID
! j-s pair:
  zetaj_s    = ztildej - yjs/(alphajs*yjsQ)
  zetas_j    = ztildes - yjs/(alphajs*yjsQ)
!  zetajs     = yjs/(alphajs*yjstildeQ)*(ztildes - ztildej)
! \zeta_{js} can be chosen to be zero in NNLO:
  zetajs     = 0
  kttildej%p = zetaj_s*p(rads)%p - zetas_j*p(radj)%p &
             + zetajs*ptildejs%p
  kttildej%flv = emitjID
! Safety check upon branching:
  if ((radiID+radrID).ne.emitiID) then
    print *,"Problem in MapMomCirjsFF..."
    write(*,'(a,3(1x,I0))') "emiti,radi,radr: ",emiti,radi,radr
    write(*,'(a,3(1x,I0))') "emitiID,radiID,radrID: ",emiti,radi,radr
    stop "MapMomCirjsFF"
  end if
  if ((radjID+radsID).ne.emitjID) then
    print *,"Problem in MapMomCirjsFF..."
    write(*,'(a,3(1x,I0))') "emitj,radj,rads: ",emitj,radj,rads
    write(*,'(a,3(1x,I0))') "emitjID,radjID,radsID: ",emitj,radj,rads
    stop "MapMomCirjsFF"
  end if
  ptildeir%flv = p(radi)%flv + p(radr)%flv
  ptildejs%flv = p(radj)%flv + p(rads)%flv
!
  call RemapMomentaO2(p,ptilde,1d0/(1d0 - alphair - alphajs), &
                      min(radi,radr),ptildeir, &
                      min(radj,rads),ptildejs, &
                      max(radi,radr),max(radj,rads))
!
end subroutine MapMomCirjsFF
!
subroutine MapMomSrs(p,ptilde,radr,rads)
use particles
use regions
use misc
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
  type(particle) , dimension(:) , intent(out) :: ptilde
  integer , intent(in) :: radr,rads
!
  integer :: ipart,jpart
  type(mom) :: Qmprs
!
!
! Obtaining y variables:
  yrQ = yiQ_arr(radr)
  ysQ = yiQ_arr(rads)
  yrsQ = yrQ + ysQ
  yrs  = sir_arr(radr,rads)/Q2
  lambdars = sqrt(abs(yrs + (1 - yrsQ)))
! If the lambdars is zero we simply quit because the subtraction term
! will acquire a zero prefactor hence not giving any contribution at all:
  if (lambdars.eq.0) return
  Qmprs = Q - p(radr)%p - p(rads)%p
! For the double soft mapping two momenta have to be omitted in the
! tilde PS:
  jpart = 1
! Running through all the momenta:
  do ipart=1,size(p)
! Skipping the soft momenta:
    if ((ipart.eq.radr).or.(ipart.eq.rads)) cycle
! For the initial state we simply copy:
    if (ipart.lt.3) ptilde(jpart) = p(ipart)
! Otherwise perform a Lorentz boost:
    if (ipart.gt.2) call LorentzLambda(Q,Qmprs/lambdars, &
                                       p(ipart)/lambdars,ptilde(jpart))
    jpart = jpart + 1
  end do
!
end subroutine MapMomSrs
!
end module nnlo_mom_maps
!
function ffunc(z0,z,p) result(f)
implicit none
!
  real(kind(1d0)) :: f
!
  real(kind(1d0)) , intent(in) :: z0,z,p
!
!
!
  f = 0
!
  if (z.gt.z0) return
  f = (1 - z)**(-p)
!
end function ffunc
!
function dfunc(m) result(d)
use regions
implicit none
!
  real(kind(1d0)) :: d
!
  integer , intent(in) :: m
!
!
  d = 2*m - 2*sub_d0
!
end function dfunc
!
function dprfunc(m) result(dpr)
use regions
implicit none
!
  real(kind(1d0)) :: dpr
!
  integer , intent(in) :: m
!
!
  dpr = m - sub_d0pr
!
end function dprfunc
!
! This routine calculates all the invariants needed to 
! construct the various subtraction terms and are general
! enough to be used at various different subtraction terms
! hence making recycling worthy.
subroutine CalcSubInvariants(p)
use regions
use particles
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
!
  integer :: ipart,jpart
  integer :: npart
!
!
! We take the number of particles from the array directly:
  npart = size(p)
!  call PrintParts(p)
! We fill up the array with a strange numerical value:
  yiQ_arr  = -1234d0
  sir_arr  = -1234d0
  yir_arr  = -1234d0
! Calculating the CM energy squared:
  Q = p(1)%p + p(2)%p
  Q2 = Q*Q
  do ipart=1,npart
    yiQ_arr(ipart) = 2d0*p(ipart)%p*Q/Q2
    do jpart=ipart+1,npart
      sir_arr(ipart,jpart) = 2d0*p(ipart)%p*p(jpart)%p
      yir_arr(ipart,jpart) = sir_arr(ipart,jpart)/Q2
      sir_arr(jpart,ipart) = sir_arr(ipart,jpart)
      yir_arr(jpart,ipart) = yir_arr(ipart,jpart)
    end do
  end do
!
end subroutine CalcSubInvariants
!
! This is a carbon-copy of the above routine but it sets up
! ivariants with tilde momenta:
subroutine CalcSubtildeInvariants(p)
use regions
use particles
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
!
  integer :: ipart,jpart
  integer :: npart
!
!
! We take the number of particles from the array directly:
  npart = size(p)
!  call PrintParts(p)
! We fill up the array with a strange numerical value:
  yitildeQ_arr  = -1234d0
  sirtilde_arr  = -1234d0
  yirtilde_arr  = -1234d0
! Note that Q2 is already available.
  do ipart=1,npart
    yitildeQ_arr(ipart) = 2d0*p(ipart)%p*Q/Q2
    do jpart=ipart+1,npart
      sirtilde_arr(ipart,jpart) = 2d0*p(ipart)%p*p(jpart)%p
      yirtilde_arr(ipart,jpart) = sirtilde_arr(ipart,jpart)/Q2
      sirtilde_arr(jpart,ipart) = sirtilde_arr(ipart,jpart)
      yirtilde_arr(jpart,ipart) = yirtilde_arr(ipart,jpart)
    end do
  end do
!
end subroutine CalcSubtildeInvariants
!
subroutine CalcPterm(yrQ,Pterm)
use alphamax
implicit none
!
  real(kind(1d0)) , intent(in) :: yrQ
  real(kind(1d0)) , intent(out) :: Pterm
!
  real(kind(1d0)) :: a1,a2
!
!
! If y0 is one we use the simplified form:
  if (y0.eq.1d0) then
    a1 = -3448d0/5d0 + 1008*log(2d0)
    a2 = 1734 - 2520*log(2d0)
  else
!    print *,"This option is not implemented yet..."
!    stop
  end if
!
!  Pterm = 1d0 - a1*yrQ - a2*yrQ**2
  Pterm = 1d0
!
end subroutine CalcPterm
!
subroutine CalcA1subtraction(p,ptilde,Bij,weightPS,Cir,Sr,CSir, &
                             CalcSMEB,CalcSMEBmunu,CalcSMEBij,A1term)
use process
use flags
use regions
use particles
use scales
use my_scale
use coupling
use misc
use phasespace
use alphamax
use nlo_mom_maps
use nlo_subtractions
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
  type(particle) , dimension(:) , intent(out) :: ptilde
  real(kind(1d0)) , dimension(:,:) , intent(out) :: Bij
  real(kind(1d0)) , intent(in) :: weightPS
  type(subterms) , intent(in) :: Cir,Sr,CSir
  real(kind(1d0)) , dimension(:) , intent(out) :: A1term
!
  integer :: iterm,jterm,iscale,icut
  integer :: nleg
  real(kind(1d0)) :: cfunc,weight
  real(kind(1d0)) :: fluxfact
  real(kind(1d0)) :: smeB
  real(kind(1d0)) , dimension(0:3,0:3) :: Bmunu
  real(kind(1d0)) :: CirTerm,SrTerm,CSirTerm
!
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
    subroutine CalcSMEBmunu(ileg,p,Bmunu)
    use particles
    implicit none
!
      integer , intent(in) :: ileg
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(0:3,0:3) , intent(out) :: Bmunu
!
    end subroutine CalcSMEBmunu
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
    subroutine apply_cuts(p,icut)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      integer , intent(out) :: icut
!
    end subroutine apply_cuts
!
    subroutine cutfunc(p,cfunc)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: cfunc
!
    end subroutine cutfunc
!
    subroutine analysis(p,weights)
    use particles
    implicit none
    !
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:) , intent(in) :: weights
    !
    end subroutine analysis
  end interface
!
  A1term = 0d0
!
! If the weight is zero we simply leave the routine:
  if (weightPS.eq.0d0) return
!
! This is the number of legs for the real process:
  nleg = size(p)
!
! We go through all the subtraction terms:
! C-type terms:
  do iterm=1,Cir%numterm
    weight = weightPS
! FF:
    if (Cir%term(iterm)%emit(1)%i.gt.2) then
      call MapMomCirFF(p,ptilde, &
                       Cir%term(iterm)%emit(1)%i,  &
                       Cir%term(iterm)%emit(1)%ID, &
                       Cir%term(iterm)%rad(1)%i,   &
                       Cir%term(iterm)%rad(1)%ID,  &
                       Cir%term(iterm)%rad(2)%i,   &
                       Cir%term(iterm)%rad(2)%ID)
      if (alphair.gt.alpha0) cycle
! We apply our cuts if any:
      if (flg_cuts) call apply_cuts(ptilde,icut)
! Or a cut function:
      if (flg_cutfunc) call cutfunc(ptilde,cfunc)
! If the PS point does not pass the cuts the term is dropped:
      if (flg_cuts.and.(icut.eq.0)) cycle
! We obtain the flux factor:
      fluxfact = CalcFluxFact(ptilde)
! Obtaining the underlying Born SME:
! If the emitter is a gluon we always need spin-correlated SME:
      if (Cir%term(iterm)%emit(1)%ID.eq.0) then
        call CalcSMEBmunu(Cir%term(iterm)%emit(1)%i,ptilde,Bmunu)
! Otherwise just the SME is requested:
      else
        call CalcSMEB(ptilde,smeB)
      end if
! We can have several different scale choices:
      SubTerm_multiscale = 0d0
      call StoreXis_scales
      do iscale=1,nscales
! We change the scales and recalculate the couplings, in our
! convention the subtraction terms use a scale derived from the
! corresponding real emission hence no call is needed to 
! calcmyscales:
        call ChangeXis_scales(iscale)
! couplings:
        call calcmyscales(p)
        call calc_couplings
!
        call CalcCirFF(p,ptilde,smeB,Bmunu,       &
                       Cir%term(iterm)%emit(1)%i, &
                       Cir%term(iterm)%rad(1)%i,  &
                       Cir%term(iterm)%rad(2)%i,  &
                       CirTerm)
! Including a tower of alphas:
! The order in alphaS is not explicitly given since this 
! routine is also used to construct NLO-type subtractions
! for the RR contribution:
        CirTerm = CirTerm & 
                * weight &
                * fluxfact &
                * alphas**((nleg - nleg_born)+border_as) &
                * alphaEM**border_aEM &
                * Cir%term(iterm)%symfact
        SubTerm_multiscale(iscale) = CirTerm
      end do
!
! If there is a request for an analysis we call it with the
! contribution just calculated:
      if (.not.flg_optim.and.flg_analysis.and.(weightPS.ne.0)) then
        call analysis(ptilde,-SubTerm_multiscale)
      end if
! If needed we multiply by the cutfunction:
      if (flg_cutfunc) SubTerm_multiscale  = SubTerm_multiscale * cfunc
      A1term = A1term + SubTerm_multiscale
! Restoring the old xi values:
      call RestoreXis_scales
    end if
  end do
  do iterm=1,Sr%numterm
    weight = weightPS
    call MapMomSr(p,ptilde, &
                  Sr%term(iterm)%rad(1)%i, &
                  Sr%term(iterm)%rad(1)%ID)
    if (yrQ.gt.y0) cycle
! We apply our cuts if any:
    if (flg_cuts) call apply_cuts(ptilde,icut)
! Or a cut function:
    if (flg_cutfunc) call cutfunc(ptilde,cfunc)
! If the PS point does not pass the cuts, the cut function is
! zero or the weight is zero we simply skip the term:
    if ((flg_cuts.and.(icut.eq.0)).or. &
        (flg_cutfunc.and.(cfunc.eq.0d0)).or. &
        (weight.eq.0d0)) cycle
! We obtain the flux factor:
    fluxfact = CalcFluxFact(ptilde)
! Only calculating the color-correlated underlying Born SME,
! the Born SME is calculated from it:
    call CalcSMEBij(ptilde,Bij)
    call CastSMEijToSME(ptilde,Bij,smeB)
!
    SubTerm_multiscale = 0d0
    call StoreXis_scales
    do iscale=1,nscales
      call ChangeXis_scales(iscale)
      call calcmyscales(p)
      call calc_couplings
!
      call CalcSr(p,ptilde,Bij, &
                  Sr%term(iterm)%rad(1)%i, &
                  SrTerm)
!
      SrTerm = SrTerm &
             * weight &
             * fluxfact &
             * alphas**((nleg - nleg_born)+border_as) &
             * alphaEM**border_aEM &
             * Sr%term(iterm)%symfact
      SubTerm_multiscale(iscale) = SrTerm
! CS-type terms:
! If the soft gluon is the same the momentum mapping is essentially the
! same for both the Sr and CirSr term hence it is not needed to calculate
! the mapping igain:
      do jterm=1,CSir%numterm
        if ((Sr%term(iterm)%rad(1)%i.eq.CSir%term(jterm)%rad(2)%i).and. &
            (Sr%term(iterm)%rad(1)%ID.eq.CSir%term(jterm)%rad(2)%ID)) then
          if (CSir%term(jterm)%emit(1)%i.gt.2) then
            call CalcCirFFSr(p,ptilde,smeB, &
                             CSir%term(jterm)%emit(1)%i, &
                             CSir%term(jterm)%rad(1)%i,  &
                             CSir%term(jterm)%rad(2)%i,  &
                             CSirTerm)
          else
            print *,"Error in CalcA1subtraction CSirIF..."
            stop
          end if
          CSirTerm = CSirTerm &
                   * weight &
                   * fluxfact &
                   * alphas**((nleg - nleg_born)+border_as) &
                   * alphaEM**border_aEM &
                   * CSir%term(jterm)%symfact
          SubTerm_multiscale(iscale) = SubTerm_multiscale(iscale) &
                                     - CSirTerm
        end if
      end do
    end do
    if (.not.flg_optim.and.flg_analysis.and.(weightPS.ne.0)) then
      call analysis(ptilde,-SubTerm_multiscale)
    end if
! Multiplying with the cut function:
    if (flg_cutfunc) SubTerm_multiscale = SubTerm_multiscale * cfunc
! Adding the contributions just calculated:
    A1term = A1term + SubTerm_multiscale
! Restoring the old xi values:
    call RestoreXis_scales
  end do
!
end subroutine CalcA1subtraction
!
! This routine calculates all the subtrations which are
! defined on m+1 parton final states, that is for R,RV and
! RR,A1. The behavior is governed by flags accordingly:
subroutine Calc_mp1_subs(p,ptilde,Bij,Vij,Bijk,Bijkl,Bmunuij,weightPS, &
                         Cir,Sr,CirSr, &
                         subsR,subsRV,subsRRA1)
use regions
use flags
use particles
use scales
use my_scale
use coupling
use phasespace
use alphamax
use misc
use histo
use nlo_subtractions
use nnlo_subtractions
use nlo_mom_maps
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
  type(particle) , dimension(:) , intent(out) :: ptilde
  real(kind(1d0)) , dimension(:,:) , intent(out) :: Bij
  real(kind(1d0)) , dimension(:,:) , intent(out) :: Vij
  real(kind(1d0)) , dimension(:,:,:) , intent(out) :: Bijk
  real(kind(1d0)) , dimension(:,:,:,:) , intent(out) :: Bijkl
  real(kind(1d0)) , dimension(0:,0:,:,:) , intent(out) :: Bmunuij
  real(kind(1d0)) , intent(inout) :: weightPS
  type(subterms) , intent(in) :: Cir,Sr,CirSr
  real(kind(1d0)) , dimension(:) , intent(out) :: subsR
  real(kind(1d0)) , dimension(:) , intent(out) :: subsRV
  real(kind(1d0)) , dimension(:) , intent(out) :: subsRRA1
!
  integer :: iterm,jterm,iscale,icut,inan
  real(kind(1d0)) :: fluxfact
  real(kind(1d0)) :: cfunc
  real(kind(1d0)) :: smeB,smeV
  real(kind(1d0)) , dimension(0:3,0:3) :: Bmunu
  real(kind(1d0)) :: SubTerm,SubCont
!
!
  interface
    subroutine CalcSubInvariants(p)
    use regions
    use particles
!
      type(particle) , dimension(:) , intent(in) :: p
!
    end subroutine CalcSubInvariants
!
    subroutine CalcSubtildeInvariants(p)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
!
    end subroutine CalcSubtildeInvariants
!
    subroutine CalcB(p,smeB)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: smeB
!
    end subroutine CalcB
!
    subroutine CalcBmunu(ileg,p,Bmunu)
    use particles
    implicit none
!
      integer , intent(in) :: ileg
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(0:3,0:3) , intent(out) :: Bmunu
!
    end subroutine CalcBmunu
!
    subroutine CalcBij(p,Bij)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Bij
!
    end subroutine CalcBij
!
    subroutine CalcVij(p,Vij,VijLaurent)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Vij
      real(kind(1d0)) , optional , dimension(:,:,-4:) , intent(out) :: VijLaurent
!
    end subroutine CalcVij
!
    subroutine CalcBmunuij(ileg,p,Bmunuij)
    use particles
    implicit none
!
      integer , intent(in) :: ileg
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(0:,0:,:,:) , intent(out) :: Bmunuij
!
    end subroutine CalcBmunuij
!
    subroutine CalcBijkl(p,Bijkl)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:,:,:,:) , intent(out) :: Bijkl
!
    end subroutine CalcBijkl
!
    subroutine Cuts(p,cfunc)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: cfunc
!
    end subroutine Cuts
!
    subroutine yminCut(p,icut)
    use momenta
    implicit none
!
      type(mom) , dimension(:) , intent(in) :: p
      integer , intent(out) :: icut
!
    end subroutine yminCut
!
    subroutine analysis(p,weights)
    use particles
    implicit none
    !
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:) , intent(in) :: weights
    !
    end subroutine analysis
  end interface
!
  if (flg_NLO_R_A1)        subsR = 0d0
  if (flg_NNLO_RV_A1)     subsRV = 0d0
  if (flg_NNLO_RRA1_A1) subsRRA1 = 0d0
! If no subtractions are needed for the m+1 parton final state
! we simply quit the routine:
  if (((.not.flg_NLO_R_A1).and. &
       (.not.flg_NNLO_RV_A1).and. &
       (.not.flg_NNLO_RRA1_A1)).or.(weightPS.eq.0d0)) return
!
! Calculating invariants:
  call CalcSubInvariants(p)
  call calcmyscales(p)
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@ Cir @@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! Considering Cir-type subtractions:
  do iterm=1,Cir%numterm
! FF:
    if (Cir%term(iterm)%emit(1)%i.gt.2) then
      call MapMomCirFF(p,ptilde, &
                       Cir%term(iterm)%emit(1)%i,  &
                       Cir%term(iterm)%emit(1)%ID, &
                       Cir%term(iterm)%rad(1)%i,   &
                       Cir%term(iterm)%rad(1)%ID,  &
                       Cir%term(iterm)%rad(2)%i,   &
                       Cir%term(iterm)%rad(2)%ID)
! ymin and NaN tests:
      call yminCut(ptilde(:)%p,icut)
      call IsNaNmom(ptilde,inan)
      if ((icut.eq.0).or.(inan.eq.0)) then
        call drop_hist
        weightPS = 0
        return
      end if
! alphamax cut:
      if (alphair.gt.alpha0) cycle
! Applying cuts, if any and cycle if not passing:
      call Cuts(ptilde,cfunc)
      if (cfunc.eq.0d0) cycle
! Obtaining a flux factor:
      fluxfact = CalcFluxFact(ptilde)
! Since the same underlying Born quantity is used at various
! places it is better to calculate it once and reuse it if needed:
! If the emitter was a gluon spin-correlated SMEs are naturally needed:
      if (Cir%term(iterm)%emit(1)%ID.eq.0) then
! If no subtraction is issued for the I_1\otimes R term only the
! spin-correlated SME is needed:
        if (.not.flg_NNLO_RRA1_A1) then
          call CalcBmunu(Cir%term(iterm)%emit(1)%i,ptilde,Bmunu)
! If subtractions are needed for the I1_R term also the simultaneously
! spin- and color-correlated SME is needed, 
! all the others can be obtained from this one:
        else
          call CalcBmunuij(Cir%term(iterm)%emit(1)%i,ptilde,Bmunuij)
          call CastSMEmunuijToSMEmunu(.false.,ptilde,Bmunuij,Bmunu)
          call CastSMEmunuijToSMEij(ptilde,Bmunuij,Bij)
          call CastSMEijToSME(ptilde,Bij,smeB)
        end if
! If the emitter was a quark the only possible splitting is q -> q g,
! hence the AP kernels only involve the underlying Born SME and not the
! spin-correlated one:
      else
! Without subtractions for R_I1 the Born SME is calculated:
        if (.not.flg_NNLO_RRA1_A1) then
          call CalcB(ptilde,smeB)
! If subtractions are requested for the I1_R term the color-correlated
! underlying Born SME has to be calculated:
        else
          call CalcBij(ptilde,Bij)
          call CastSMEijToSME(ptilde,Bij,smeB)
        end if
      end if
      SubTerm_multiscale = 0d0
      call StoreXis_scales
      do iscale=1,nscales
        call ChangeXis_scales(iscale)
        call calcmyscales(p)
        call calc_couplings
! Calculation of Cir for R, if needed:
!###########################################
!############### Cir for R #################
!###########################################
        if (flg_NLO_R_A1) then
          call CalcCirFF(p,ptilde,smeB,Bmunu,       &
                         Cir%term(iterm)%emit(1)%i, &
                         Cir%term(iterm)%rad(1)%i,  &
                         Cir%term(iterm)%rad(2)%i,  &
                         SubTerm)
          SubTerm = SubTerm               & 
                  * weightPS              &
                  * fluxfact              &
                  * alphas**(border_as+1) &
                  * alphaEM**border_aEM   &
                  * Cir%term(iterm)%symfact
!
          subsR(iscale) = subsR(iscale) &
                        + SubTerm * cfunc
!
          SubTerm_multiscale(iscale) = SubTerm
        end if
! Calculation of Cir for RV, if needed:
!###########################################
!############## Cir for RV #################
!###########################################
        if (flg_NNLO_RV_A1) then
          call CalcCir01FF(p,ptilde, &
                           Cir%term(iterm)%emit(1)%i, &
                           Cir%term(iterm)%rad(1)%i,  &
                           Cir%term(iterm)%rad(2)%i,  &
                           SubTerm)
          call CalcCir10FF(p,ptilde,smeB,Bmunu, &
                           Cir%term(iterm)%emit(1)%i, &
                           Cir%term(iterm)%rad(1)%i,  &
                           Cir%term(iterm)%rad(2)%i,  &
                           SubCont)
          SubTerm = SubTerm + SubCont
          SubTerm = SubTerm               & 
                  * weightPS              &
                  * fluxfact              &
                  * alphas**(border_as+2) &
                  * alphaEM**border_aEM   &
                  * Cir%term(iterm)%symfact
!
          subsRV(iscale) = subsRV(iscale) &
                         + SubTerm * cfunc
!
          SubTerm_multiscale(iscale) = SubTerm_multiscale(iscale) &
                                     + SubTerm
        end if
! Calculation of Cir for RRA1, if needed:
!###########################################
!############## Cir for RRA1 ###############
!###########################################
        if (flg_NNLO_RRA1_A1) then
          call CalcSubtildeInvariants(ptilde)
          call CalcCir00IFF(p,ptilde,smeB,Bij,Bmunu,Bmunuij, &
                            Cir%term(iterm)%emit(1)%i,       &
                            Cir%term(iterm)%rad(1)%i,        &
                            Cir%term(iterm)%rad(2)%i,        &
                            SubTerm)
          call CalcCirR00FF(p,ptilde,smeB,Bmunu,             &
                            Cir%term(iterm)%emit(1)%i, &
                            Cir%term(iterm)%rad(1)%i,  &
                            Cir%term(iterm)%rad(2)%i,  &
                            SubCont)
          SubTerm = SubTerm + SubCont
          SubTerm = SubTerm               & 
                  * weightPS              &
                  * fluxfact              &
                  * alphas**(border_as+2) &
                  * alphaEM**border_aEM   &
                  * Cir%term(iterm)%symfact
!
          subsRRA1(iscale) = subsRRA1(iscale) &
                           + SubTerm * cfunc
!
          SubTerm_multiscale(iscale) = SubTerm_multiscale(iscale) &
                                     + SubTerm
        end if
      end do
! If an analysis is enabled and not being in optimization mode
! the analysis routine is called:
      if (.not.flg_optim.and.flg_analysis.and.(weightPS.ne.0)) then
        call analysis(ptilde,-SubTerm_multiscale)
      end if
! Restoring the old xi values:
      call RestoreXis_scales
! 
    else
      print *,"Initial state radiation is not treated yet..."
      stop "Calc_mp1_subs"
    end if
  end do
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@ Sr @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! Considering Sr-type subtractions:
  do iterm=1,Sr%numterm
    call MapMomSr(p,ptilde, &
                  Sr%term(iterm)%rad(1)%i, &
                  Sr%term(iterm)%rad(1)%ID)
! ymin and NaN tests:
    call yminCut(ptilde(:)%p,icut)
    call IsNaNmom(ptilde,inan)
    if ((icut.eq.0).or.(inan.eq.0)) then
      call drop_hist
      weightPS = 0
      return
    end if
! alphamax cut:
    if (yrQ.gt.y0) cycle
! Applying cuts, if any and cycle if not passing:
    call Cuts(ptilde,cfunc)
    if (cfunc.eq.0d0) cycle
! Obtaining a flux factor:
    fluxfact = CalcFluxFact(ptilde)
! For better performance the SMEs are calculated only onces:
    if (flg_NNLO_RV_A1) then
      call CalcVij(ptilde,Vij)
      call CastSMEijToSME(ptilde,Vij,smeV)
    end if
    if (.not.flg_NNLO_RRA1_A1) then
      call CalcBij(ptilde,Bij)
      call CastSMEijToSME(ptilde,Bij,smeB)
    else
      call CalcBijkl(ptilde,Bijkl)
      call CastSMEijklToSMEij(.false.,ptilde,Bijkl,Bij)
      call CastSMEijToSME(ptilde,Bij,smeB)
    end if
! Multiple scale choices are allowed, nullify collector:
    SubTerm_multiscale = 0d0
    call StoreXis_scales
    do iscale=1,nscales
      call ChangeXis_scales(iscale)
      call calcmyscales(p)
      call calc_couplings
! Calculation of Sr for R, if needed:
!###########################################
!############### Sr for R ##################
!###########################################
      if (flg_NLO_R_A1) then
        call CalcSr(p,ptilde,Bij, &
                    Sr%term(iterm)%rad(1)%i, &
                    SubTerm)
        SubTerm = SubTerm               &
                * weightPS              &
                * fluxfact              &
                * alphas**(border_as+1) &
                * alphaEM**border_aEM   &
                * Sr%term(iterm)%symfact
!
        subsR(iscale) = subsR(iscale) &
                      + SubTerm * cfunc
!
        SubTerm_multiscale(iscale) = SubTerm
      end if
! Calculation of Sr for RV, if needed:
!###########################################
!############### Sr for RV #################
!###########################################
      if (flg_NNLO_RV_A1) then
        call CalcSr01(p,ptilde,Vij,            &
                      Sr%term(iterm)%rad(1)%i, &
                      SubTerm)
        call CalcSr10(p,ptilde,Bij,Bijk,       &
                      Sr%term(iterm)%rad(1)%i, &
                      SubCont)
        SubTerm = SubTerm + SubCont
        SubTerm = SubTerm               &
                * weightPS              &
                * fluxfact              &
                * alphas**(border_as+2) &
                * alphaEM**border_aEM   &
                * Sr%term(iterm)%symfact
!
        subsRV(iscale) = subsRV(iscale) &
                       + SubTerm * cfunc
!
        SubTerm_multiscale(iscale) = SubTerm_multiscale(iscale) &
                                   + SubTerm
      end if
! Calculation of Sr for RRA1, if needed:
!###########################################
!############## Sr for RRA1 ################
!###########################################
      if (flg_NNLO_RRA1_A1) then
! Additional invariants are needed:
        call CalcSubtildeInvariants(ptilde)
        call CalcSr00I(p,ptilde,Bij,Bijkl,      &
                       Sr%term(iterm)%rad(1)%i, &
                       SubTerm)
        call CalcSrR00(p,ptilde,Bij,            &
                       Sr%term(iterm)%rad(1)%i, &
                       SubCont)
        SubTerm = SubTerm + SubCont
        SubTerm = SubTerm               &
                * weightPS              &
                * fluxfact              &
                * alphas**(border_as+2) &
                * alphaEM**border_aEM   &
                * Sr%term(iterm)%symfact
!
        subsRRA1(iscale) = subsRRA1(iscale) &
                         + SubTerm * cfunc
!
        SubTerm_multiscale(iscale) = SubTerm_multiscale(iscale) &
                                   + SubTerm
      end if
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@ CirSr @@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! Considering CirSr-type subtractions:
      do jterm=1,CirSr%numterm
        if ((Sr%term(iterm)%rad(1)%i.eq.CirSr%term(jterm)%rad(2)%i).and. &
            (Sr%term(iterm)%rad(1)%ID.eq.CirSr%term(jterm)%rad(2)%ID)) then
          if (CirSr%term(jterm)%emit(1)%i.gt.2) then
! Calculation of CirSr for R, if needed:
!###########################################
!############### CirSr for R ###############
!###########################################
            if (flg_NLO_R_A1) then
              call CalcCirFFSr(p,ptilde,smeB, &
                               CirSr%term(jterm)%emit(1)%i, &
                               CirSr%term(jterm)%rad(1)%i,  &
                               CirSr%term(jterm)%rad(2)%i,  &
                               SubTerm)
              SubTerm = SubTerm               &
                      * weightPS              &
                      * fluxfact              &
                      * alphas**(border_as+1) &
                      * alphaEM**border_aEM   &
                      * CirSr%term(jterm)%symfact
!
              subsR(iscale) = subsR(iscale) &
                            - SubTerm * cfunc
!
              SubTerm_multiscale(iscale) = SubTerm_multiscale(iscale) &
                                         - SubTerm
            end if
! Calculation of CirSr for RV, if needed:
!###########################################
!############## CirSr for RV ###############
!###########################################
            if (flg_NNLO_RV_A1) then
              call CalcCirFFSr01(p,ptilde,smeV, &
                                 CirSr%term(jterm)%emit(1)%i, &
                                 CirSr%term(jterm)%rad(1)%i,  &
                                 CirSr%term(jterm)%rad(2)%i,  &
                                 SubTerm)
              call CalcCirFFSr10(p,ptilde,smeB, &
                                 CirSr%term(jterm)%emit(1)%i, &
                                 CirSr%term(jterm)%rad(1)%i,  &
                                 CirSr%term(jterm)%rad(2)%i,  &
                                 SubCont)
              SubTerm = SubTerm + SubCont
              SubTerm = SubTerm               &
                      * weightPS              &
                      * fluxfact              &
                      * alphas**(border_as+2) &
                      * alphaEM**border_aEM   &
                      * CirSr%term(jterm)%symfact
!
              subsRV(iscale) = subsRV(iscale) &
                             - SubTerm * cfunc
!
              SubTerm_multiscale(iscale) = SubTerm_multiscale(iscale) &
                                         - SubTerm
            end if
! Calculation of CirSr for RRA1, if needed:
!###########################################
!############# CirSr for RRA1 ##############
!###########################################
            if (flg_NNLO_RRA1_A1) then
              call CalcCirFFSr00I(p,ptilde,smeB,Bij,           &
                                  CirSr%term(jterm)%emit(1)%i, &
                                  CirSr%term(jterm)%rad(1)%i,  &
                                  CirSr%term(jterm)%rad(2)%i,  &
                                  SubTerm)
              call CalcCirFFSrR00(p,ptilde,smeB,               &
                                  CirSr%term(jterm)%emit(1)%i, &
                                  CirSr%term(jterm)%rad(1)%i,  &
                                  CirSr%term(jterm)%rad(2)%i,  &
                                  SubCont)
              SubTerm = SubTerm + SubCont
              SubTerm = SubTerm               &
                      * weightPS              &
                      * fluxfact              &
                      * alphas**(border_as+2) &
                      * alphaEM**border_aEM   &
                      * CirSr%term(jterm)%symfact
!
              subsRRA1(iscale) = subsRRA1(iscale) &
                               - SubTerm * cfunc
!
              SubTerm_multiscale(iscale) = SubTerm_multiscale(iscale) &
                                         - SubTerm
            end if
          else
            print *,"Initial state ratiadiot is not yet treated..."
            stop "Calc_mp1_subs"
          end if
        end if
      end do
    end do
! If an analysis is enabled and not being in optimization mode
! the analysis routine is called:
    if (.not.flg_optim.and.flg_analysis.and.(weightPS.ne.0)) then
      call analysis(ptilde,-SubTerm_multiscale)
    end if
! Restoring the old xi values:
    call RestoreXis_scales
  end do
!
!
end subroutine Calc_mp1_subs
!
subroutine Calc_mp2_subs(p,phat,ptilde,         &
                         Bij,Bijkl,Bmunuij,Rij, &
                         weightPS,              &
                         Cir,Sr,CirSr,          &
                         Cirs,CSirs,Cirjs,Srs,  &
                         subsRR)
use regions
use flags
use particles
use scales
use my_scale
use coupling
use phasespace
use alphamax
use misc
use QCDparams
use nlo_subtractions
use nnlo_subtractions
use nlo_mom_maps
use nnlo_mom_maps
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
  type(particle) , dimension(:) , intent(out) :: phat,ptilde
  real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
  real(kind(1d0)) , dimension(:,:,:,:) , intent(inout) :: Bijkl
  real(kind(1d0)) , dimension(0:,0:,:,:) , intent(inout) :: Bmunuij
  real(kind(1d0)) , dimension(:,:) , intent(out) :: Rij
  real(kind(1d0)) , intent(in) :: weightPS
  type(subterms) , intent(in) :: Cir,Sr,CirSr
  type(subterms) , intent(in) :: Cirs,CSirs,Cirjs,Srs
  real(kind(1d0)) , dimension(:) , intent(out) :: subsRR
!
  integer :: iterm,jterm,iscale
  integer :: radi,radj,radr,rads,radk,radt,rad1
  integer :: rad_kt,rad_i,rad_r,rad_k,rad_ir
  integer :: emitir,emitkt,emitirs,emitjs,emitirt
  real(kind(1d0)) :: cfunc,fluxfact,ff
  real(kind(1d0)) :: SME,SubTerm,SubCont
  real(kind(1d0)) :: ytQ,yrhatQ
  real(kind(1d0)) , dimension(0:3,0:3) :: SMEmunu,Bmunu
  real(kind(1d0)) , dimension(0:3,0:3,0:3,0:3) :: Balbemunu
!
  logical , parameter :: show_terms = .false.
!
  real(kind(1d0)) , external :: ffunc,dfunc,dprfunc
!
  interface
    subroutine CalcSubInvariants(p)
    use regions
    use particles
!
      type(particle) , dimension(:) , intent(in) :: p
!
    end subroutine CalcSubInvariants
!
    subroutine CalcB(parts,smeB)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: parts
      real(kind(1d0)) , intent(out) :: smeB
!
    end subroutine CalcB
!
    subroutine CalcBij(p,Bij)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Bij
!
    end subroutine CalcBij
!
    subroutine CalcBmunu(ileg,p,Bmunu)
    use particles
    implicit none
!
      integer , intent(in) :: ileg
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(0:3,0:3) , intent(out) :: Bmunu
!
    end subroutine CalcBmunu
!
    subroutine CalcBmunuij(ileg,p,Bmunuij)
    use particles
    implicit none
!
      integer , intent(in) :: ileg
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(0:,0:,:,:) , intent(out) :: Bmunuij
!
    end subroutine CalcBmunuij
!
    subroutine CalcBijkl(p,Bijkl)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:,:,:,:) , intent(out) :: Bijkl
!
    end subroutine CalcBijkl
!
    subroutine CalcBalbemunu(ileg,jleg,p,Balbemunu)
    use particles
    implicit none
!
      integer , intent(in) :: ileg,jleg
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(0:3,0:3,0:3,0:3) , intent(out) :: Balbemunu
!
    end subroutine CalcBalbemunu
!
    subroutine CalcR(parts,smeR)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: parts
      real(kind(1d0)) , intent(out) :: smeR
!
    end subroutine CalcR
!
    subroutine CalcRmunu(ileg,p,Rmunu)
    use particles
    implicit none
!
      integer , intent(in) :: ileg
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(0:3,0:3) , intent(out) :: Rmunu
!
    end subroutine CalcRmunu
!
    subroutine CalcRij(p,Rij)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Rij
!
    end subroutine CalcRij
!
    subroutine Cuts(p,cfunc)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: cfunc
!
    end subroutine Cuts
!
    subroutine analysis(p,weights)
    use particles
    implicit none
    !
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:) , intent(in) :: weights
    !
    end subroutine analysis
  end interface
!
  subsRR = 0
!
! If no NNLO subtraction is needed for the RR we return:
  if (.not.(flg_NNLO_RR_A1.or. &
         flg_NNLO_RR_A12.or. &
         flg_NNLO_RR_A2) &
      .or.(weightPS.eq.0)) return
!
! Calculating invariants:
  call CalcSubInvariants(p)
  call calcmyscales(p)
!
! *************************************************************
! ************************* A1-type ***************************
! *************************************************************
!
  if (flg_NNLO_RR_A1) then
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@ Cir @@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    do iterm=1,Cir%numterm
! FF:
      if (Cir%term(iterm)%emit(1)%i.gt.2) then
        call MapMomCirFF(p,phat, &
                         Cir%term(iterm)%emit(1)%i,  &
                         Cir%term(iterm)%emit(1)%ID, &
                         Cir%term(iterm)%rad(1)%i,   &
                         Cir%term(iterm)%rad(1)%ID,  &
                         Cir%term(iterm)%rad(2)%i,   &
                         Cir%term(iterm)%rad(2)%ID)
! alphamax cut:
        if (alphair.gt.alpha0) cycle
! Applying cuts, if any and cycle if not passing:
        call Cuts(phat,cfunc)
        if (cfunc.eq.0d0) cycle
! Obtaining a flux factor:
        fluxfact = CalcFluxFact(phat)
! If the emitter was a gluon spin-correlated SMEs are naturally needed:
        if (Cir%term(iterm)%emit(1)%ID.eq.0) then
          call CalcRmunu(Cir%term(iterm)%emit(1)%i,phat,SMEmunu)
! If the emitter was a quark the only possible splitting is q -> q g,
! hence the AP kernels only involve the underlying Born SME and not the
! spin-correlated one:
        else
          call CalcR(phat,SME)
        end if
        SubTerm_multiscale = 0d0
        call StoreXis_scales
        do iscale=1,nscales
          call ChangeXis_scales(iscale)
          call calcmyscales(p)
          call calc_couplings
          call CalcCirFF(p,phat,SME,SMEmunu,        &
                         Cir%term(iterm)%emit(1)%i, &
                         Cir%term(iterm)%rad(1)%i,  &
                         Cir%term(iterm)%rad(2)%i,  &
                         SubTerm)
          if (show_terms) then
            write(*,"(A,2(I0,1x),G0)")                   &
              "Cir,i,r,Term: ",Cir%term(iterm)%rad(1)%i, &
                               Cir%term(iterm)%rad(2)%i, &
                               SubTerm*Cir%term(iterm)%symfact
!            write(*,"(G0)") SubTerm               & 
!                          * weightPS              &
!                          * fluxfact              &
!                          * alphas**(border_as+2) &
!                          * alphaEM**border_aEM   &
!                          * Cir%term(iterm)%symfact
          end if
          SubTerm = SubTerm               & 
                  * weightPS              &
                  * fluxfact              &
                  * alphas**(border_as+2) &
                  * alphaEM**border_aEM   &
                  * Cir%term(iterm)%symfact
!
          subsRR(iscale) = subsRR(iscale) &
                         + SubTerm * cfunc
!
          SubTerm_multiscale(iscale) = SubTerm
        end do
! If an analysis is enabled and not being in optimization mode
! the analysis routine is called:
        if (.not.flg_optim.and.flg_analysis.and.(weightPS.ne.0)) then
          call analysis(phat,-SubTerm_multiscale)
        end if
! Restoring the old xi values:
        call RestoreXis_scales
      end if
    end do
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@ Sr @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    do iterm=1,Sr%numterm
      call MapMomSr(p,phat, &
                    Sr%term(iterm)%rad(1)%i, &
                    Sr%term(iterm)%rad(1)%ID)
! alphamax cut:
      if (yrQ.gt.y0) cycle
! Applying cuts, if any and cycle if not passing:
      call Cuts(phat,cfunc)
      if (cfunc.eq.0d0) cycle
! Obtaining a flux factor:
      fluxfact = CalcFluxFact(phat)
      call CalcRij(phat,Rij)
      call CastSMEijToSME(phat,Rij,SME)
! Multiple scale choices are allowed, nullify collector:
      SubTerm_multiscale = 0d0
      call StoreXis_scales
      do iscale=1,nscales
        call ChangeXis_scales(iscale)
        call calcmyscales(p)
        call calc_couplings
        call CalcSr(p,phat,Rij, &
                    Sr%term(iterm)%rad(1)%i, &
                    SubTerm)
        if (show_terms) then
          write(*,"(A,I0,1x,G0)")                  &
            "Sr,r,Term: ",Sr%term(iterm)%rad(1)%i, &
                          SubTerm*Sr%term(iterm)%symfact
!          write(*,"(G0)") SubTerm               & 
!                        * weightPS              &
!                        * fluxfact              &
!                        * alphas**(border_as+2) &
!                        * alphaEM**border_aEM   &
!                        * Sr%term(iterm)%symfact
        end if
        SubTerm = SubTerm               &
                * weightPS              &
                * fluxfact              &
                * alphas**(border_as+2) &
                * alphaEM**border_aEM   &
                * Sr%term(iterm)%symfact
!
        subsRR(iscale) = subsRR(iscale) &
                       + SubTerm * cfunc
!
        SubTerm_multiscale(iscale) = SubTerm
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@ CirSr @@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        do jterm=1,CirSr%numterm
          if ((Sr%term(iterm)%rad(1)%i.eq.CirSr%term(jterm)%rad(2)%i).and. &
              (Sr%term(iterm)%rad(1)%ID.eq.CirSr%term(jterm)%rad(2)%ID)) then
            if (CirSr%term(jterm)%emit(1)%i.gt.2) then
              call CalcCirFFSr(p,phat,SME,                  &
                               CirSr%term(jterm)%emit(1)%i, &
                               CirSr%term(jterm)%rad(1)%i,  &
                               CirSr%term(jterm)%rad(2)%i,  &
                               SubTerm)
              if (show_terms) then
                write(*,"(A,2(I0,1x),G0)")                       &
                  "CirSr,i,r,Term: ",CirSr%term(jterm)%rad(1)%i, &
                                     CirSr%term(jterm)%rad(2)%i, &
                                     -SubTerm * CirSr%term(jterm)%symfact
!                write(*,"(G0)") -SubTerm              & 
!                              * weightPS              &
!                              * fluxfact              &
!                              * alphas**(border_as+2) &
!                              * alphaEM**border_aEM   &
!                              * CirSr%term(jterm)%symfact
              end if
              SubTerm = SubTerm               &
                      * weightPS              &
                      * fluxfact              &
                      * alphas**(border_as+2) &
                      * alphaEM**border_aEM   &
                      * CirSr%term(jterm)%symfact
!
              subsRR(iscale) = subsRR(iscale) &
                             - SubTerm * cfunc
!
              SubTerm_multiscale(iscale) = SubTerm_multiscale(iscale) &
                                         - SubTerm
            end if
          end if
        end do
      end do
! If an analysis is enabled and not being in optimization mode
! the analysis routine is called:
      if (.not.flg_optim.and.flg_analysis.and.(weightPS.ne.0)) then
        call analysis(phat,-SubTerm_multiscale)
      end if
! Restoring the old xi values:
      call RestoreXis_scales
    end do
  end if
!
! *************************************************************
! ************************* A2-type ***************************
! *************************************************************
!
  if (flg_NNLO_RR_A2) then
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@ Cirs @@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    do iterm=1,Cirs%numterm
      if (Cirs%term(iterm)%emit(1)%i.lt.3) then
        print *,"Initial state is not treated yet..."
        stop "Calc_mp2_subs"
      else
        call MapMomCirsFF(p,ptilde,                    &
                          Cirs%term(iterm)%emit(1)%i,  &
                          Cirs%term(iterm)%emit(1)%ID, &
                          Cirs%term(iterm)%rad(1)%i,   &
                          Cirs%term(iterm)%rad(1)%ID,  &
                          Cirs%term(iterm)%rad(2)%i,   &
                          Cirs%term(iterm)%rad(2)%ID,  &
                          Cirs%term(iterm)%rad(3)%i,   &
                          Cirs%term(iterm)%rad(3)%ID)
! alphamax cut:
        ff = ffunc(alpha0,alphairs,dfunc(size(p)-4))
        if (ff.eq.0) cycle
! Applying cuts, if any and cycle if not passing:
        call Cuts(ptilde,cfunc)
        if (cfunc.eq.0d0) cycle
! Obtaining a flux factor:
        fluxfact = CalcFluxFact(ptilde)
! Obtaining the underlying Born SME:
! If the emitter is a gluon we always need spin-correlated SME:
        if (Cirs%term(iterm)%emit(1)%ID.eq.0) then
          call CalcBmunu(Cirs%term(iterm)%emit(1)%i,ptilde,SMEmunu)
          call CastSMEmunuToSME(ptilde,SMEmunu,SME)
! Otherwise just the SME is requested:
        else
          call CalcB(ptilde,SME)
        end if
! The Counterterm is computationally heavy it is better to
! calculate it only once:
        call CalcCirsFF(p,ptilde,SME,SMEmunu,       &
                        Cirs%term(iterm)%emit(1)%i, &
                        Cirs%term(iterm)%rad(1)%i,  &
                        Cirs%term(iterm)%rad(2)%i,  &
                        Cirs%term(iterm)%rad(3)%i,  &
                        SubTerm)
        if (show_terms) then
          write(*,"(A,3(I0,1x),G0)")                   &
            "Cirs,i,r,s,Term: ",Cirs%term(iterm)%rad(1)%i, & 
                                Cirs%term(iterm)%rad(2)%i, & 
                                Cirs%term(iterm)%rad(3)%i, & 
                                SubTerm*Cirs%term(iterm)%symfact
        end if
      end if
! Including prefactors not changing with scales:
      SubTerm = SubTerm   &
              * weightPS  &
              * fluxfact  &
              * ff        &
              * Cirs%term(iterm)%symfact
      SubTerm_multiscale = 0d0
      call StoreXis_scales
      do iscale=1,nscales
        call ChangeXis_scales(iscale)
        call calcmyscales(p)
        call calc_couplings
!
        subsRR(iscale) = subsRR(iscale)        &
                       + SubTerm               &
                       * cfunc                 &
                       * alphas**(border_as+2) &
                       * alphaEM**border_aEM
!
        SubTerm_multiscale(iscale) = SubTerm               & 
                                   * alphas**(border_as+2) &
                                   * alphaEM**border_aEM
      end do
! If an analysis is enabled and not being in optimization mode
! the analysis routine is called:
      if (.not.flg_optim.and.flg_analysis.and.(weightPS.ne.0)) then
        call analysis(ptilde,-SubTerm_multiscale)
      end if
! Restoring the old xi values:
      call RestoreXis_scales
    end do
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@ Cirjs @@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    do iterm=1,Cirjs%numterm
      if (Cirjs%term(iterm)%emit(1)%i.lt.3) then
        print *,"Initial state is not treated yet..."
        stop "Calc_mp2_subs"
      else
        call MapMomCirjsFF(p,ptilde, &
                           Cirjs%term(iterm)%emit(1)%i,  &
                           Cirjs%term(iterm)%emit(1)%ID, &
                           Cirjs%term(iterm)%rad(1)%i,   &
                           Cirjs%term(iterm)%rad(1)%ID,  &
                           Cirjs%term(iterm)%rad(2)%i,   &
                           Cirjs%term(iterm)%rad(2)%ID,  &
                           Cirjs%term(iterm)%emit(2)%i,  &
                           Cirjs%term(iterm)%emit(2)%ID, &
                           Cirjs%term(iterm)%rad(3)%i,   &
                           Cirjs%term(iterm)%rad(3)%ID,  &
                           Cirjs%term(iterm)%rad(4)%i,   &
                           Cirjs%term(iterm)%rad(4)%ID)
! alphamax cut:
        ff = ffunc(alpha0,alphair+alphajs,dfunc(size(p)-4))
        if (ff.eq.0) cycle
! Applying cuts, if any and cycle if not passing:
        call Cuts(ptilde,cfunc)
        if (cfunc.eq.0d0) cycle
! Obtaining a flux factor:
        fluxfact = CalcFluxFact(ptilde)
! Obtaining the underlying Born SME:
! If both emitters are quarks only the Born SME has to be
! calculated:
        if ((Cirjs%term(iterm)%emit(1)%ID.ne.0).and. &
            (Cirjs%term(iterm)%emit(2)%ID.ne.0)) then
          call CalcB(ptilde,SME)
! The first emitter is a quark, the second is a gluon:
        elseif ((Cirjs%term(iterm)%emit(1)%ID.ne.0).and. &
                (Cirjs%term(iterm)%emit(2)%ID.eq.0)) then
          call CalcBmunu(Cirjs%term(iterm)%emit(2)%i,ptilde,Bmunu)
          call CastSMEmunuToSME(ptilde,Bmunu,SME)
! The first emitter is a gluon, the second is a quark:
        elseif ((Cirjs%term(iterm)%emit(1)%ID.eq.0).and. &
                (Cirjs%term(iterm)%emit(2)%ID.ne.0)) then
          call CalcBmunu(Cirjs%term(iterm)%emit(1)%i,ptilde,SMEmunu)
          call CastSMEmunuToSME(ptilde,SMEmunu,SME)
! Both of the emitters are gluons:
        elseif ((Cirjs%term(iterm)%emit(1)%ID.eq.0).and. &
                (Cirjs%term(iterm)%emit(2)%ID.eq.0)) then
          call CalcBalbemunu(Cirjs%term(iterm)%emit(1)%i, &
                             Cirjs%term(iterm)%emit(2)%i, &
                             ptilde,Balbemunu)
          call CastSMEalbemunuToSMEmunu(ptilde,Balbemunu,SMEmunu,Bmunu)
          call CastSMEmunuToSME(ptilde,SMEmunu,SME)
        end if
! The Counterterm is computationally heavy it is better to
! calculate it only once:
        call CalcCirjsFF(p,ptilde,                    &
                         sme,SMEmunu,Bmunu,Balbemunu, &
                         Cirjs%term(iterm)%emit(1)%i, &
                         Cirjs%term(iterm)%rad(1)%i,  &
                         Cirjs%term(iterm)%rad(2)%i,  &
                         Cirjs%term(iterm)%emit(2)%i, &
                         Cirjs%term(iterm)%rad(3)%i,  &
                         Cirjs%term(iterm)%rad(4)%i,  &
                         SubTerm)
        if (show_terms) then
          write(*,"(A,4(I0,1x),G0)")                   &
            "Cirjs,i,r,j,s,Term: ",Cirjs%term(iterm)%rad(1)%i, & 
                                   Cirjs%term(iterm)%rad(2)%i, & 
                                   Cirjs%term(iterm)%rad(3)%i, & 
                                   Cirjs%term(iterm)%rad(4)%i, & 
                                   SubTerm*Cirjs%term(iterm)%symfact
        end if
      end if
! Including prefactors not changing with scales:
      SubTerm = SubTerm   &
              * weightPS  &
              * fluxfact  &
              * ff        &
              * Cirjs%term(iterm)%symfact
      SubTerm_multiscale = 0d0
      call StoreXis_scales
      do iscale=1,nscales
        call ChangeXis_scales(iscale)
        call calcmyscales(p)
        call calc_couplings
!
        subsRR(iscale) = subsRR(iscale)        &
                       + SubTerm               &
                       * cfunc                 &
                       * alphas**(border_as+2) &
                       * alphaEM**border_aEM
!
        SubTerm_multiscale(iscale) = SubTerm               & 
                                   * alphas**(border_as+2) &
                                   * alphaEM**border_aEM
      end do
! If an analysis is enabled and not being in optimization mode
! the analysis routine is called:
      if (.not.flg_optim.and.flg_analysis.and.(weightPS.ne.0)) then
        call analysis(ptilde,-SubTerm_multiscale)
      end if
! Restoring the old xi values:
      call RestoreXis_scales
    end do
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@ CSirs,CirsCSirs,CirjsCSirs @@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    do iterm=1,CSirs%numterm
      if (CSirs%term(iterm)%emit(1)%i.lt.3) then
        print *,"Initial state is not treated yet..."
        stop "Calc_mp2_subs"
      else
        call MapMomCirFF(p,phat,                       &
                         CSirs%term(iterm)%emit(1)%i,  &
                         CSirs%term(iterm)%emit(1)%ID, &
                         CSirs%term(iterm)%rad(1)%i,   &
                         CSirs%term(iterm)%rad(1)%ID,  &
                         CSirs%term(iterm)%rad(2)%i,   &
                         CSirs%term(iterm)%rad(2)%ID)
! As for the soft mapping two possibilities can arise:
! The soft leg has the highest position, hence the intermediate
! position is shifted by one:
        if (max(CSirs%term(iterm)%rad(1)%i,     &
                CSirs%term(iterm)%rad(2)%i).lt. &
                CSirs%term(iterm)%rad(3)%i) then
          call MapMomSrNNLO(phat,ptilde,                    &
                            CSirs%term(iterm)%rad(3)%i,     &
                            CSirs%term(iterm)%rad(3)%i - 1, &
                            CSirs%term(iterm)%rad(3)%ID)
        else
          call MapMomSrNNLO(phat,ptilde,                &
                            CSirs%term(iterm)%rad(3)%i, &
                            CSirs%term(iterm)%rad(3)%i, &
                            CSirs%term(iterm)%rad(3)%ID)
        end if
! alphamax cut:
! Note that in arXiv:1301.3504 y_{\hat{s}Q} was used in the 
! second instance of the y function, but due to the 
! organization of the calculation the value of ysQ is stored 
! in yrQ:
        ff = ffunc(alpha0,alphair,dfunc(size(p)-4)) &
           * ffunc(y0,yrQ,dprfunc(size(p)-4))
        if (ff.eq.0) cycle
! Applying cuts, if any and cycle if not passing:
        call Cuts(ptilde,cfunc)
        if (cfunc.eq.0d0) cycle
! Obtaining a flux factor:
        fluxfact = CalcFluxFact(ptilde)
! Obtaining the underlying Born SME:
! If the emitter is a gluon we always need color- and 
! spin-correlated SME:
        if (CSirs%term(iterm)%emit(1)%ID.eq.0) then
          call CalcBmunuij(CSirs%term(iterm)%emit(1)%i,ptilde,Bmunuij)
          call CastSMEmunuijToSMEij(ptilde,Bmunuij,Bij)
          call CastSMEmunuijToSMEmunu(.false.,ptilde,Bmunuij,SMEmunu)
          call CastSMEijToSME(ptilde,Bij,SME)
! Otherwise just the color-correlated SME is requested:
        else
          call CalcBij(ptilde,Bij)
          call CastSMEijToSME(ptilde,Bij,SME)
        end if
! The Counterterm is computationally heavy it is better to
! calculate it only once:
!##############################################################
!########################## CSirs #############################
!##############################################################
        call CalcCSirsFF(p,ptilde,SME,Bij,Bmunuij,    &
                         CSirs%term(iterm)%emit(1)%i, &
                         CSirs%term(iterm)%rad(1)%i,  &
                         CSirs%term(iterm)%rad(2)%i,  &
                         CSirs%term(iterm)%rad(3)%i,  &
                         SubTerm)
        if (show_terms) then
          write(*,"(A,3(I0,1x),G0)")                   &
            "CSirs,i,r,s,Term: ",CSirs%term(iterm)%rad(1)%i, & 
                                 CSirs%term(iterm)%rad(2)%i, & 
                                 CSirs%term(iterm)%rad(3)%i, & 
                                 SubTerm*CSirs%term(iterm)%symfact
        end if
!##############################################################
!####################### CirsCSirs ############################
!##############################################################
        call CalcCirsCSirsFF(p,ptilde,SME,SMEmunu,        &
                             CSirs%term(iterm)%emit(1)%i, &
                             CSirs%term(iterm)%rad(1)%i,  &
                             CSirs%term(iterm)%rad(2)%i,  &
                             CSirs%term(iterm)%rad(3)%i,  &
                             SubCont)
        if (show_terms) then
          write(*,"(A,3(I0,1x),G0)")                   &
            "CirsCSirs,i,r,s,Term: ",CSirs%term(iterm)%rad(1)%i, & 
                                     CSirs%term(iterm)%rad(2)%i, & 
                                     CSirs%term(iterm)%rad(3)%i, & 
                                     -SubCont*CSirs%term(iterm)%symfact
        end if
        SubTerm = SubTerm - SubCont
!##############################################################
!####################### CirjsCSirs ###########################
!##############################################################
! For the CirjsCSirs counterterm an additional leg has to be
! selected:
        do radj=1,size(p)
! Only partons count:
          if (abs(p(radj)%flv).gt.qcd_nf) cycle
! Omit already selected legs:
          if ((radj.eq.CSirs%term(iterm)%rad(1)%i).or. &
              (radj.eq.CSirs%term(iterm)%rad(2)%i).or. & 
              (radj.eq.CSirs%term(iterm)%rad(3)%i)) cycle
!          print *,"radj: ",radj
! Determining the emitter for the radj rads parton-pair:
          emitjs = min(radj,CSirs%term(iterm)%rad(3)%i)
! The position of the emitter for the j-s pair can change
! at the underlying Born level if the first pair comes 
! before:
          if (max(CSirs%term(iterm)%rad(1)%i, &
                  CSirs%term(iterm)%rad(2)%i).lt.emitjs) then
            emitjs = emitjs - 1 
          end if
!          print *,"emitjs: ",emitjs
          call CalcCirjsCSirsFF(p,ptilde,SME,SMEmunu,        &
                                CSirs%term(iterm)%emit(1)%i, &
                                CSirs%term(iterm)%rad(1)%i,  &
                                CSirs%term(iterm)%rad(2)%i,  &
                                emitjs,radj,                 &
                                CSirs%term(iterm)%rad(3)%i,  &
                                SubCont)
          if (show_terms) then
            write(*,"(A,4(I0,1x),G0)")                                &
              "CirjsCSirs,i,r,j,s,Term: ",CSirs%term(iterm)%rad(1)%i, & 
                                          CSirs%term(iterm)%rad(2)%i, & 
                                          radj,                       &
                                          CSirs%term(iterm)%rad(3)%i, & 
                                          -SubCont*CSirs%term(iterm)%symfact
          end if
          SubTerm = SubTerm - SubCont
        end do
      end if
! Including prefactors not changing with scales:
      SubTerm = SubTerm   &
              * weightPS  &
              * fluxfact  &
              * ff        &
              * CSirs%term(iterm)%symfact
      SubTerm_multiscale = 0d0
      call StoreXis_scales
      do iscale=1,nscales
        call ChangeXis_scales(iscale)
        call calcmyscales(p)
        call calc_couplings
!
        subsRR(iscale) = subsRR(iscale)        &
                       + SubTerm               &
                       * cfunc                 &
                       * alphas**(border_as+2) &
                       * alphaEM**border_aEM
!
        SubTerm_multiscale(iscale) = SubTerm               & 
                                   * alphas**(border_as+2) &
                                   * alphaEM**border_aEM
      end do
! If an analysis is enabled and not being in optimization mode
! the analysis routine is called:
      if (.not.flg_optim.and.flg_analysis.and.(weightPS.ne.0)) then
        call analysis(ptilde,-SubTerm_multiscale)
      end if
! Restoring the old xi values:
      call RestoreXis_scales
    end do
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@ Srs,CirsSrs,CSirsSrs,CirsCSirsSrs,CirjsSrs @@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    do iterm=1,Srs%numterm
      call MapMomSrs(p,ptilde, &
                     Srs%term(iterm)%rad(1)%i,   &
                     Srs%term(iterm)%rad(2)%i)
! If lambdars is zero the subtraction term does not give
! any contribution at all due to zero prefactor:
      if (lambdars.eq.0) cycle
! alphamax cut:
      ff = ffunc(y0,yrQ + ysQ - yrs,dprfunc(size(p)-4))
      if (ff.eq.0) cycle
! Applying cuts, if any and cycle if not passing:
      call Cuts(ptilde,cfunc)
      if (cfunc.eq.0d0) cycle
! Obtaining a flux factor:
      fluxfact = CalcFluxFact(ptilde)
! Obtaining the underlying Born SME:
      call CalcBijkl(ptilde,Bijkl)
      call CastSMEijklToSMEij(.false.,ptilde,Bijkl,Bij)
      call CastSMEijToSME(ptilde,Bij,SME)
! The Counterterm is computationally heavy it is better to
! calculate it only once:
!##############################################################
!############################ Srs #############################
!##############################################################
      call CalcSrs(p,ptilde,Bij,Bijkl,       &
                   Srs%term(iterm)%rad(1)%i, &
                   Srs%term(iterm)%rad(2)%i, &
                   SubTerm)
      if (show_terms) then
        write(*,"(A,2(I0,1x),G0)")                   &
          "Srs,r,s,Term: ",Srs%term(iterm)%rad(1)%i, & 
                           Srs%term(iterm)%rad(2)%i, & 
                           SubTerm*Srs%term(iterm)%symfact
      end if
! For the upcoming counterterms an additional leg has to be
! selected:
      do radi=1,size(p)
! Only partons count:
        if (abs(p(radi)%flv).gt.qcd_nf) cycle
! Omit already selected legs:
        if ((radi.eq.Srs%term(iterm)%rad(1)%i).or. &
            (radi.eq.Srs%term(iterm)%rad(2)%i)) cycle
!        print *,"radi: ",radi
        emitirs = min(radi,Srs%term(iterm)%rad(1)%i, &
                           Srs%term(iterm)%rad(2)%i)
!##############################################################
!########################## CirsSrs ###########################
!##############################################################
        call CalcCirsSrs(p,ptilde,SME,             &
                         emitirs,radi,             &
                         Srs%term(iterm)%rad(1)%i, &
                         Srs%term(iterm)%rad(2)%i, &
                         SubCont)
        if (show_terms) then
          write(*,"(A,3(I0,1x),G0)")                              &
            "CirsSrs,i,r,s,Term: ",radi,Srs%term(iterm)%rad(1)%i, & 
                                   Srs%term(iterm)%rad(2)%i,      &
                                   -SubCont*Srs%term(iterm)%symfact
        end if
        SubTerm = SubTerm - SubCont
! All the other counterterms are defined only if the soft legs
! belong to gluons:
        if (Srs%term(iterm)%rad(1)%ID.eq.0) then
! When CSirsSrs terms are considered as a distinction has to be
! made between legs r and s, since the i-r pair becomes a 
! collinear pair:
          do rad1=1,2
            if (rad1.eq.1) then
              radr = Srs%term(iterm)%rad(1)%i
              rads = Srs%term(iterm)%rad(2)%i
            else
              radr = Srs%term(iterm)%rad(2)%i
              rads = Srs%term(iterm)%rad(1)%i
            end if
            emitir = min(radi,radr)
            if (emitir.gt.rads) emitir = emitir - 1
!            print *,"radr,rads: ",radr,rads
!##############################################################
!######################### CSirsSrs ###########################
!##############################################################
            call CalcCSirsSrs(p,ptilde,Bij,          &
                              emitir,radi,radr,rads, &
                              SubCont)
            if (show_terms) then
              write(*,"(A,3(I0,1x),G0)")                &
                "CSirsSrs,i,r,s,Term: ",radi,radr,rads, & 
                                        -SubCont*Srs%term(iterm)%symfact
            end if
            SubTerm = SubTerm - SubCont
!##############################################################
!######################## CirsCSirsSrs ########################
!##############################################################
            call CalcCirsCSirsSrs(p,ptilde,SME,           &
                                  emitirs,radi,radr,rads, &
                                  SubCont)
            if (show_terms) then
              write(*,"(A,3(I0,1x),G0)")                    &
                "CirsCSirsSrs,i,r,s,Term: ",radi,radr,rads, & 
                                            SubCont*Srs%term(iterm)%symfact
            end if
            SubTerm = SubTerm + SubCont
! For the CirjsSrs terms one additional leg appears:
            do radj=radi+1,size(p)
              if (abs(p(radj)%flv).gt.qcd_nf) cycle
              if ((radj.eq.radi).or. &
                  (radj.eq.radr).or. &
                  (radj.eq.rads)) cycle
              emitir = min(radi,radr)
              emitjs = min(radj,rads)
              if (emitir.gt.max(radj,rads)) emitir = emitir - 1
              if (emitjs.gt.max(radi,radr)) emitjs = emitjs - 1
!##############################################################
!######################### CirjsSrs ###########################
!##############################################################
              call CalcCirjsSrs(p,ptilde,SME,     &
                                emitir,radi,radr, &
                                emitjs,radj,rads, &
                                SubCont)
              if (show_terms) then
                write(*,"(A,4(I0,1x),G0)")                       &
                  "CirjsSrs,i,r,j,s,Term: ",radi,radr,radj,rads, & 
                                            SubCont*Srs%term(iterm)%symfact  
              end if
              SubTerm = SubTerm + SubCont
            end do
          end do
        end if
      end do
! Including prefactors not changing with scales:
      SubTerm = SubTerm  &
              * weightPS &
              * fluxfact &
              * ff       &
              * Srs%term(iterm)%symfact
      SubTerm_multiscale = 0d0
      call StoreXis_scales
      do iscale=1,nscales
        call ChangeXis_scales(iscale)
        call calcmyscales(p)
        call calc_couplings
!
        subsRR(iscale) = subsRR(iscale)        &
                       + SubTerm               &
                       * cfunc                 &
                       * alphas**(border_as+2) &
                       * alphaEM**border_aEM
!
        SubTerm_multiscale(iscale) = SubTerm               & 
                                   * alphas**(border_as+2) &
                                   * alphaEM**border_aEM
      end do
! If an analysis is enabled and not being in optimization mode
! the analysis routine is called:
      if (.not.flg_optim.and.flg_analysis.and.(weightPS.ne.0)) then
        call analysis(ptilde,-SubTerm_multiscale)
      end if
! Restoring the old xi values:
      call RestoreXis_scales
    end do
  end if
!
! *************************************************************
! ************************* A12-type **************************
! *************************************************************
!
  if (flg_NNLO_RR_A12) then
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@ CktCktr @@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    do iterm=1,Cirs%numterm
! The CktCktr counterterm is derived from the Cirs one,
! two legs are dedicated (k-t) these two have to be picked:
! Have to check whether the k-t pair can be produced in a 
! splitting.
! If the splitting is possible a possible reordering takes 
! place.
      do rad1=1,3
        if (.not.(((Cirs%term(iterm)%rad(mod(rad1+2,3)+1)%ID            &
                  + Cirs%term(iterm)%rad(mod(rad1+3,3)+1)%ID).eq.0).or. &
                   (Cirs%term(iterm)%rad(mod(rad1+2,3)+1)%ID.eq.0).or.  &
                   (Cirs%term(iterm)%rad(mod(rad1+3,3)+1)%ID.eq.0))) cycle
        call ReorderRadir(Cirs%term(iterm)%rad(mod(rad1+2,3)+1)%i,  &
                          Cirs%term(iterm)%rad(mod(rad1+2,3)+1)%ID, &
                          Cirs%term(iterm)%rad(mod(rad1+3,3)+1)%i,  &
                          Cirs%term(iterm)%rad(mod(rad1+3,3)+1)%ID, &
                          radk,radt)
        radr = Cirs%term(iterm)%rad(mod(rad1+1,3)+1)%i
! Leg positions in the intermediate step:
        rad_kt = min(radk,radt)
        rad_r = radr
        if (radr.gt.max(radk,radt)) rad_r = radr - 1
        if (Cirs%term(iterm)%emit(1)%i.lt.3) then
          print *,"Initial state is not treated yet..."
          stop "Calc_mp2_subs"
        else
          call MapMomCirFF_A12(p,phat,                &
                               radk,radt,             &
                               Q,Q2,yiQ_arr,sir_arr,  &
                               s_k_t,ztildek,ztildet, &
                               kttildek,kttildet,     &
                               alphakt)
          call CalcSubhatInvariants(phat)
          call MapMomCirFF_A12(phat,ptilde,                &
                               rad_kt,rad_r,               &
                               Q,Q2,yihatQ_arr,sirhat_arr, &
                               s_kt_r,zhatkt,zhatr,        &
                               kthatkt,kthatr,             &
                               alphaktr)
        end if
! alphamax cut:
        ff = ffunc(alpha0,alphakt,dfunc(size(p)-4)) &
           * ffunc(alpha0,alphaktr,dfunc(size(p)-4))
        if (ff.eq.0) cycle
! Applying cuts, if any and cycle if not passing:
        call Cuts(ptilde,cfunc)
        if (cfunc.eq.0d0) cycle
! Obtaining a flux factor:
        fluxfact = CalcFluxFact(ptilde)
! If the emitter is a gluon we always need spin-correlated SME:
        if (Cirs%term(iterm)%emit(1)%ID.eq.0) then
          call CalcBmunu(Cirs%term(iterm)%emit(1)%i,ptilde,SMEmunu)
          call CastSMEmunuToSME(ptilde,SMEmunu,SME)
! Otherwise just the SME is requested:
        else
          call CalcB(ptilde,SME)
        end if
        if (Cirs%term(iterm)%emit(1)%i.lt.3) then
          print *,"Initial state is not treated yet..."
          stop "Calc_mp2_subs"
        else
!##############################################################
!######################### CktCktr ############################
!##############################################################
          call CalcCktCktrFF(p,phat,ptilde,SME,SMEmunu,   &
                             Cirs%term(iterm)%emit(1)%i,  &
                             radk,radt,radr,rad_kt,rad_r, &
                             SubTerm)
          if (show_terms) then
            write(*,"(A,3(I0,1x),G0)")               &
              "CktCktr,k,t,r,Term: ",radk,radt,radr, &
                                     -SubTerm*Cirs%term(iterm)%symfact
!            write(*,"(G0)") -SubTerm              & 
!                          * weightPS              &
!                          * fluxfact              &
!                          * alphas**(border_as+2) &
!                          * alphaEM**border_aEM   &
!                          * Cirs%term(iterm)%symfact
          end if
        end if
! Including prefactors not changing with scales:
        SubTerm = SubTerm  &
                * weightPS &
                * fluxfact &
                * ff       &
                * Cirs%term(iterm)%symfact
        SubTerm_multiscale = 0d0
        call StoreXis_scales
        do iscale=1,nscales
          call ChangeXis_scales(iscale)
          call calcmyscales(p)
          call calc_couplings
!
          subsRR(iscale) = subsRR(iscale)        &
                         - SubTerm               &
                         * cfunc                 &
                         * alphas**(border_as+2) &
                         * alphaEM**border_aEM
!
          SubTerm_multiscale(iscale) = SubTerm               & 
                                     * alphas**(border_as+2) &
                                     * alphaEM**border_aEM
        end do
! If an analysis is enabled and not being in optimization mode
! the analysis routine is called:
        if (.not.flg_optim.and.flg_analysis.and.(weightPS.ne.0)) then
          call analysis(ptilde,+SubTerm_multiscale)
        end if
! Restoring the old xi values:
        call RestoreXis_scales
      end do
    end do
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@ CktCirkt @@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    do iterm=1,Cirjs%numterm
! The CktCirkt counterterm is derived from the Cirjs one,
! two legs are dedicated (k-t) these two have to be picked:
! Note that in this case no reordering nor check upon 
! splitting is needed:
      do rad1=1,2
        radi = Cirjs%term(iterm)%rad(2*rad1-1)%i
        radr = Cirjs%term(iterm)%rad(2*rad1)%i
        radk = Cirjs%term(iterm)%rad(mod(2*rad1+1,4))%i
        radt = Cirjs%term(iterm)%rad(mod(2*rad1+1,4)+1)%i
        emitir = Cirjs%term(iterm)%emit(rad1)%i
        emitkt = Cirjs%term(iterm)%emit(mod(rad1,2)+1)%i
! The first mapping combines parton k and t, hence the 
! position of i and r can change:
        rad_i = radi
        rad_r = radr
        if (radi.gt.max(radk,radt)) rad_i = radi - 1
        if (radr.gt.max(radk,radt)) rad_r = radr - 1
! Both mappings are of collniear-type, hence they can be in
! the intial or final state:
        if (emitkt.lt.3) then
          print *,"ISR is not treated in Calc_mp2_subs..."
          stop "Calc_mp2_subs"
        else
          call MapMomCirFF_A12(p,phat,                &
                               radk,radt,             &
                               Q,Q2,yiQ_arr,sir_arr,  &
                               s_k_t,ztildek,ztildet, &
                               kttildek,kttildet,     &
                               alphakt)
        end if
        call CalcSubhatInvariants(phat)
        if (emitir.lt.3) then
          print *,"ISR is not treated in Calc_mp2_subs..."
          stop "Calc_mp2_subs"
        else
          call MapMomCirFF_A12(phat,ptilde,                &
                               rad_i,rad_r,                &  
                               Q,Q2,yihatQ_arr,sirhat_arr, &
                               sir,zhati,zhatr,            &
                               kthati,kthatr,              &
                               alphair)
        end if
! alphamax cut:
        ff = ffunc(alpha0,alphakt,dfunc(size(p)-4)) &
           * ffunc(alpha0,alphair,dfunc(size(p)-4))
        if (ff.eq.0) cycle
! Applying cuts, if any and cycle if not passing:
        call Cuts(ptilde,cfunc)
        if (cfunc.eq.0d0) cycle
! Obtaining a flux factor:
        fluxfact = CalcFluxFact(ptilde)
! Obtaining the underlying Born SME:
! If both emitters are quarks only the Born SME has to be
! calculated:
        if ((ptilde(emitir)%flv.ne.0).and. &
            (ptilde(emitkt)%flv.ne.0)) then
!          print *,"q-q"
          call CalcB(ptilde,SME)
! The first emitter is a quark, the second is a gluon:
        elseif ((ptilde(emitir)%flv.ne.0).and. &
                (ptilde(emitkt)%flv.eq.0)) then
!          print *,"q-g"
          call CalcBmunu(emitkt,ptilde,Bmunu)
          call CastSMEmunuToSME(ptilde,Bmunu,SME)
! The first emitter is a gluon, the second is a quark:
        elseif ((ptilde(emitir)%flv.eq.0).and. &
                (ptilde(emitkt)%flv.ne.0)) then
!          print *,"g-q"
          call CalcBmunu(emitir,ptilde,SMEmunu)
          call CastSMEmunuToSME(ptilde,SMEmunu,SME)
! Both of the emitters are gluons:
        elseif ((ptilde(emitir)%flv.eq.0).and. &
                (ptilde(emitkt)%flv.eq.0)) then
!          print *,"g-g"
          call CalcBalbemunu(emitir,emitkt, &
                             ptilde,Balbemunu)
          call CastSMEalbemunuToSMEmunu(ptilde,Balbemunu,SMEmunu,Bmunu)
          call CastSMEmunuToSME(ptilde,SMEmunu,SME)
        end if
!##############################################################
!######################### CktCirkt ###########################
!##############################################################
        if ((emitkt.gt.2).and.(emitir.gt.2)) then
          call CalcCktCirktFF(p,phat,ptilde,SME,SMEmunu,Bmunu,Balbemunu, &
                              radi,radr,radk,radt, &
                              SubTerm)
          if (show_terms) then
            write(*,"(A,4(I0,1x),G0)")                       &
              "CktCirkt,i,r,k,t,Term: ",radi,radr,radk,radt, &
                                        -SubTerm*Cirjs%term(iterm)%symfact
!            write(*,"(G0)") -SubTerm              & 
!                          * weightPS              &
!                          * fluxfact              &
!                          * alphas**(border_as+2) &
!                          * alphaEM**border_aEM   &
!                          * Cirjs%term(iterm)%symfact
          end if
        else
          print *,"ISR is not treated in Calc_mp2_subs..."
          stop "Calc_mp2_subs"
        end if
! Including prefactors not changing with scales:
        SubTerm = SubTerm  &
                * weightPS &
                * fluxfact &
                * ff       &
                * Cirjs%term(iterm)%symfact
        SubTerm_multiscale = 0d0
        call StoreXis_scales
        do iscale=1,nscales
          call ChangeXis_scales(iscale)
          call calcmyscales(p)
          call calc_couplings
!
          subsRR(iscale) = subsRR(iscale)        &
                         - SubTerm               &
                         * cfunc                 &
                         * alphas**(border_as+2) &
                         * alphaEM**border_aEM
!
          SubTerm_multiscale(iscale) = SubTerm               & 
                                     * alphas**(border_as+2) &
                                     * alphaEM**border_aEM
        end do
! If an analysis is enabled and not being in optimization mode
! the analysis routine is called:
        if (.not.flg_optim.and.flg_analysis.and.(weightPS.ne.0)) then
          call analysis(ptilde,+SubTerm_multiscale)
        end if
! Restoring the old xi values:
        call RestoreXis_scales
      end do
    end do
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@ CktCSktr,CktCktrCSktr,CktCirktCSktr @@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    do iterm=1,CSirs%numterm
! Position of leg kt in the hat-ed configuration:
      rad_kt = min(CSirs%term(iterm)%rad(1)%i, &
                   CSirs%term(iterm)%rad(2)%i)
      rad_r = CSirs%term(iterm)%rad(3)%i
      if (rad_r.gt.max(CSirs%term(iterm)%rad(1)%i, &
                       CSirs%term(iterm)%rad(2)%i)) &
        rad_r = rad_r - 1
      if (CSirs%term(iterm)%emit(1)%i.lt.3) then
        print *,"ISR is not treated in Calc_mp2_subs..."
        stop "Calc_mp2_subs"
      else
        call MapMomCirFF_A12(p,phat,                     &
                             CSirs%term(iterm)%rad(1)%i, &
                             CSirs%term(iterm)%rad(2)%i, &
                             Q,Q2,yiQ_arr,sir_arr,       &
                             s_k_t,ztildek,ztildet,      &
                             kttildek,kttildet,          &
                             alphakt)
      end if
      call CalcSubhatInvariants(phat)
      call MapMomSr_A12(phat,ptilde,     &
                        rad_r,           &
                        Q,Q2,yihatQ_arr, &
                        yrhatQ,lambda_r)
      if (lambda_r.eq.0) cycle
! alphamax cut:
      ff = ffunc(alpha0,alphakt,dfunc(size(p)-4)) &
         * ffunc(y0,yrhatQ,dprfunc(size(p)-4))
      if (ff.eq.0) cycle
! Applying cuts, if any and cycle if not passing:
      call Cuts(ptilde,cfunc)
      if (cfunc.eq.0d0) cycle
! The soft mapping involves a boost this has to be applied
! even on the transverse momenta:
      call LorentzLambda(Q,(Q - phat(rad_r)%p)/lambda_r, &
                         kttildek/lambda_r,kttildek)
      call LorentzLambda(Q,(Q - phat(rad_r)%p)/lambda_r, &
                         kttildet/lambda_r,kttildet)
! Obtaining a flux factor:
      fluxfact = CalcFluxFact(ptilde)
! If the emitter is a gluon we always need spin- and 
! color-correlated SME:
      if (CSirs%term(iterm)%emit(1)%ID.eq.0) then
        call CalcBmunuij(CSirs%term(iterm)%emit(1)%i,ptilde,Bmunuij)
        call CastSMEmunuijToSMEmunu(.false.,ptilde,Bmunuij,SMEmunu)
        call CastSMEmunuToSME(ptilde,SMEmunu,SME)
! Otherwise just the color-correlated is requested:
      else
        call CalcBij(ptilde,Bij)
        call CastSMEijToSME(ptilde,Bij,SME)
      end if
      if (CSirs%term(iterm)%emit(1)%i.lt.3) then
        print *,"ISR is not treated in Calc_mp2_subs..."
        stop "Calc_mp2_subs"
      else
!##############################################################
!######################### CktCSktr ###########################
!##############################################################
        call CalcCktCSktrFF(p,phat,ptilde,Bij,Bmunuij,  &
                            CSirs%term(iterm)%rad(1)%i, &
                            CSirs%term(iterm)%rad(2)%i, &
                            rad_r, &
                            SubTerm)
        if (show_terms) then
          write(*,"(A,3(I0,1x),G0)")                            &
            "CktCSktr,k,t,r,Term: ",CSirs%term(iterm)%rad(1)%i, &
                                    CSirs%term(iterm)%rad(2)%i, &
                                    CSirs%term(iterm)%rad(3)%i, &
                                    -SubTerm*CSirs%term(iterm)%symfact
!          write(*,"(G0)") -SubTerm              & 
!                        * weightPS              &
!                        * fluxfact              &
!                        * alphas**(border_as+2) &
!                        * alphaEM**border_aEM   &
!                        * CSirs%term(iterm)%symfact
        end if
!##############################################################
!####################### CktCktrCSktr #########################
!##############################################################
        call CalcCktCktrCSktrFF(p,phat,ptilde,SME,SMEmunu, &
                                rad_kt,rad_r, &
                                CSirs%term(iterm)%rad(1)%i, &
                                CSirs%term(iterm)%rad(2)%i, &
                                SubCont)
        if (show_terms) then
          write(*,"(A,3(I0,1x),G0)")                                &
            "CktCktrCSktr,k,t,r,Term: ",CSirs%term(iterm)%rad(1)%i, &
                                        CSirs%term(iterm)%rad(2)%i, &
                                        CSirs%term(iterm)%rad(3)%i, &
                                        SubCont*CSirs%term(iterm)%symfact
!          write(*,"(G0)") SubCont              & 
!                        * weightPS              &
!                        * fluxfact              &
!                        * alphas**(border_as+2) &
!                        * alphaEM**border_aEM   &
!                        * CSirs%term(iterm)%symfact
        end if
        SubTerm = SubTerm - SubCont
        do radi=1,size(p)
! Only partons count:
          if (abs(p(radi)%flv).gt.qcd_nf) cycle
! Omit already selected legs:
          if ((radi.eq.CSirs%term(iterm)%rad(1)%i).or. &
              (radi.eq.CSirs%term(iterm)%rad(2)%i).or. &
              (radi.eq.CSirs%term(iterm)%rad(3)%i)) cycle
          rad_i = radi
          if (rad_i.gt.max(CSirs%term(iterm)%rad(1)%i, &
                           CSirs%term(iterm)%rad(2)%i)) &
            rad_i = rad_i - 1
!##############################################################
!####################### CktCirktCSktr ########################
!##############################################################
          call CalcCktCirktCSktrFF(p,phat,ptilde,SME,SMEmunu,  &
                                   rad_i,rad_r,                &
                                   CSirs%term(iterm)%rad(1)%i, &
                                   CSirs%term(iterm)%rad(2)%i, &
                                   SubCont)
          if (show_terms) then
            write(*,"(A,4(I0,1x),G0)")                                        &
              "CktCirktCSktr,i,r,k,t,Term: ",radi,CSirs%term(iterm)%rad(3)%i, &
                                             CSirs%term(iterm)%rad(1)%i,      &
                                             CSirs%term(iterm)%rad(2)%i,      &
                                             SubCont*CSirs%term(iterm)%symfact
!            write(*,"(G0)") SubCont              & 
!                          * weightPS              &
!                          * fluxfact              &
!                          * alphas**(border_as+2) &
!                          * alphaEM**border_aEM   &
!                          * CSirs%term(iterm)%symfact
          end if
          SubTerm = SubTerm - SubCont
        end do
      end if
! Including prefactors not changing with scales:
      SubTerm = SubTerm  &
              * weightPS &
              * fluxfact &
              * ff       &
              * CSirs%term(iterm)%symfact
      SubTerm_multiscale = 0d0
      call StoreXis_scales
      do iscale=1,nscales
        call ChangeXis_scales(iscale)
        call calcmyscales(p)
        call calc_couplings
!
        subsRR(iscale) = subsRR(iscale)        &
                       - SubTerm               &
                       * cfunc                 &
                       * alphas**(border_as+2) &
                       * alphaEM**border_aEM
!
        SubTerm_multiscale(iscale) = SubTerm               & 
                                   * alphas**(border_as+2) &
                                   * alphaEM**border_aEM
      end do
! If an analysis is enabled and not being in optimization mode
! the analysis routine is called:
      if (.not.flg_optim.and.flg_analysis.and.(weightPS.ne.0)) then
        call analysis(ptilde,+SubTerm_multiscale)
      end if
! Restoring the old xi values:
      call RestoreXis_scales
    end do
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@ CktSkt,CktCrktSkt @@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    do iterm=1,Srs%numterm
      rad_kt = min(Srs%term(iterm)%rad(1)%i, &
                   Srs%term(iterm)%rad(2)%i)
      call MapMomCirFF_A12(p,phat,                   &
                           Srs%term(iterm)%rad(1)%i, &
                           Srs%term(iterm)%rad(2)%i, &
                           Q,Q2,yiQ_arr,sir_arr,     &
                           s_k_t,ztildek,ztildet,    &
                           kttildek,kttildet,        &
                           alphakt)
      call CalcSubhatInvariants(phat)
! Note that the notation can be misleading instead of 
! yrhatQ and lambda_r we should have used ykthatQ and
! lambda_kt:
      call MapMomSr_A12(phat,ptilde,     &
                        rad_kt,          &
                        Q,Q2,yihatQ_arr, &
                        yrhatQ,lambda_r)
      if (lambda_r.eq.0) cycle
! alphamax cut:
      ff = ffunc(alpha0,alphakt,dfunc(size(p)-4)) &
         * ffunc(y0,yrhatQ,dprfunc(size(p)-4))
      if (ff.eq.0) cycle
! Applying cuts, if any and cycle if not passing:
      call Cuts(ptilde,cfunc)
      if (cfunc.eq.0d0) cycle
! Obtaining a flux factor:
      fluxfact = CalcFluxFact(ptilde)
! Obtaining the underlying Born SME:
      call CalcBij(ptilde,Bij)
      call CastSMEijToSME(ptilde,Bij,SME)
!##############################################################
!######################### CktSkt #############################
!##############################################################
      call CalcCktSktFF(p,phat,ptilde,Bij, &
                        rad_kt, &
                        Srs%term(iterm)%rad(1)%i, &
                        Srs%term(iterm)%rad(2)%i, &
                        SubTerm)
      if (show_terms) then
        write(*,"(A,2(I0,1x),G0)")                      &
          "CktSkt,k,t,Term: ",Srs%term(iterm)%rad(1)%i, &
                              Srs%term(iterm)%rad(2)%i, & 
                              -SubTerm*Srs%term(iterm)%symfact
!        write(*,"(G0)") -SubTerm              & 
!                      * weightPS              &
!                      * fluxfact              &
!                      * alphas**(border_as+2) &
!                      * alphaEM**border_aEM   &
!                      * Srs%term(iterm)%symfact
      end if
      do radr=1,size(p)
! Only partons count:
        if (abs(p(radr)%flv).gt.qcd_nf) cycle
! Omit already selected legs:
        if ((radr.eq.Srs%term(iterm)%rad(1)%i).or. &
            (radr.eq.Srs%term(iterm)%rad(2)%i)) cycle
        rad_r = radr
        if (rad_r.gt.max(Srs%term(iterm)%rad(1)%i, &
                         Srs%term(iterm)%rad(2)%i)) &
          rad_r = rad_r - 1
!##############################################################
!######################## CktCrktSkt ##########################
!##############################################################
        call CalcCktCrktSktFF(p,phat,ptilde,SME,        &
                              rad_r,rad_kt,radr,        &
                              Srs%term(iterm)%rad(1)%i, &
                              Srs%term(iterm)%rad(2)%i, &
                              SubCont)
        if (show_terms) then
          write(*,"(A,3(I0,1x),G0)")                            &
            "CktCrktSkt,r,k,t,Term: ",radr,                     &
                                      Srs%term(iterm)%rad(1)%i, &
                                      Srs%term(iterm)%rad(2)%i, & 
                                      SubCont*Srs%term(iterm)%symfact
!          write(*,"(G0)") SubCont              & 
!                        * weightPS              &
!                        * fluxfact              &
!                        * alphas**(border_as+2) &
!                        * alphaEM**border_aEM   &
!                        * Srs%term(iterm)%symfact
        end if
        SubTerm = SubTerm - SubCont
      end do
! Including prefactors not changing with scales:
      SubTerm = SubTerm  &
              * weightPS &
              * fluxfact &
              * ff       &
              * Srs%term(iterm)%symfact
      SubTerm_multiscale = 0d0
      call StoreXis_scales
      do iscale=1,nscales
        call ChangeXis_scales(iscale)
        call calcmyscales(p)
        call calc_couplings
!
        subsRR(iscale) = subsRR(iscale)        &
                       - SubTerm               &
                       * cfunc                 &
                       * alphas**(border_as+2) &
                       * alphaEM**border_aEM
!
        SubTerm_multiscale(iscale) = SubTerm               & 
                                   * alphas**(border_as+2) &
                                   * alphaEM**border_aEM
      end do
! If an analysis is enabled and not being in optimization mode
! the analysis routine is called:
      if (.not.flg_optim.and.flg_analysis.and.(weightPS.ne.0)) then
        call analysis(ptilde,+SubTerm_multiscale)
      end if
! Restoring the old xi values:
      call RestoreXis_scales
    end do
!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!&&&&&&&&&&&&&&&&&&&&&&&&&& St-type &&&&&&&&&&&&&&&&&&&&&&&&&&&
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@ StCirt,StCSirt,StCirtCSirt @@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    do iterm=1,Cirs%numterm
! Loop over all legs of the subtractian term to determine
! candidates for t:
      do rad1=1,3
        radt = Cirs%term(iterm)%rad(mod(rad1+1,3)+1)%i
! Only gluons contribute:
        if (p(radt)%flv.ne.0) cycle
        call ReorderRadir(Cirs%term(iterm)%rad(mod(rad1+2,3)+1)%i,  &
                          Cirs%term(iterm)%rad(mod(rad1+2,3)+1)%ID, &
                          Cirs%term(iterm)%rad(mod(rad1+3,3)+1)%i,  &
                          Cirs%term(iterm)%rad(mod(rad1+3,3)+1)%ID, &
                          radi,radr)
        rad_i = radi
        rad_r = radr
        if (rad_i.gt.radt) rad_i = rad_i - 1
        if (rad_r.gt.radt) rad_r = rad_r - 1
        rad_ir = min(rad_i,rad_r)
        call MapMomSr_A12(p,phat,       &
                          radt,         &
                          Q,Q2,yiQ_arr, &
                          ytQ,lambda_t)
        call CalcSubhatInvariants(phat)
        if (Cirs%term(iterm)%emit(1)%i.lt.3) then
          print *,"ISR is not treated in Calc_mp2_subs..."
          stop "Calc_mp2_subs"
        else
          call MapMomCirFF_A12(phat,ptilde,                &
                               rad_i,rad_r,                &
                               Q,Q2,yihatQ_arr,sirhat_arr, &
                               sir,zhati,zhatr,            &
                               kthati,kthatr,              &
                               alphair)
        end if
! alphamax cut:
        ff = ffunc(y0,ytQ,dprfunc(size(p)-4)) &
           * ffunc(alpha0,alphair,dfunc(size(p)-4))
        if (ff.eq.0) cycle
! Applying cuts, if any and cycle if not passing:
        call Cuts(ptilde,cfunc)
        if (cfunc.eq.0d0) cycle
! Obtaining a flux factor:
        fluxfact = CalcFluxFact(ptilde)
! It cannot be avoided to call the SME two times (at most),
! the reason is just that we would like to treat the StCirt
! and StCSirt terms on the same footing, but the emitters
! can be different:
!        if (Cirs%term(iterm)%emit(1)%ID.eq.0) then
!          call CalcBmunu(Cirs%term(iterm)%emit(1)%i,ptilde,SMEmunu)
!          call CastSMEmunuToSME(ptilde,SMEmunu,SME)
!        end if
! If the ir parton corresponds to a gluon we have to 
! calculate the simultaneously color- and spin-correlated SME,
! note that the position of the ir parton does not coincide
! with the position of the emitter coming from Cirs subtraction
! term, hence it has to be calculated:
        if (ptilde(rad_ir)%flv.eq.0) then
          call CalcBmunuij(rad_ir,ptilde,Bmunuij)
          call CastSMEmunuijToSMEmunu(.false.,ptilde,Bmunuij,SMEmunu)
          call CastSMEmunuToSME(ptilde,SMEmunu,SME)
        else 
          call CalcBij(ptilde,Bij)
          call CastSMEijToSME(ptilde,Bij,SME)
        end if
        if (rad_ir.lt.3) then
          print *,"Initial state is not treated yet..."
          stop "Calc_mp2_subs"
        else
!##############################################################
!########################## StCirt ############################
!##############################################################
          call CalcStCirtFF(p,phat,ptilde,SME,SMEmunu, &
                            radi,radr,radt,            &
                            SubTerm)
          if (show_terms) then
            write(*,"(A,3(I0,1x),G0)")              &
              "StCirt,i,r,t,Term: ",radi,radr,radt, &
                                    -SubTerm*Cirs%term(iterm)%symfact
!            write(*,"(G0)") -SubTerm              & 
!                          * weightPS              &
!                          * fluxfact              &
!                          * alphas**(border_as+2) &
!                          * alphaEM**border_aEM   &
!                          * Cirs%term(iterm)%symfact
          end if
!##############################################################
!########################## StCSirt ###########################
!##############################################################
          call CalcStCSirtFF(p,phat,ptilde,Bij,Bmunuij, &
                             rad_ir,radi,radr,radt,     &
                             SubCont)
          SubTerm = SubTerm + SubCont
          if (show_terms) then
            write(*,"(A,3(I0,1x),G0)")               &
              "StCSirt,i,r,t,Term: ",radi,radr,radt, &
                                     -SubCont*Cirs%term(iterm)%symfact
!            write(*,"(G0)") -SubCont              & 
!                          * weightPS              &
!                          * fluxfact              &
!                          * alphas**(border_as+2) &
!                          * alphaEM**border_aEM   &
!                          * Cirs%term(iterm)%symfact
          end if
!##############################################################
!######################### StCirtCSirt ########################
!##############################################################
          call CalcStCirtCSirtFF(p,phat,ptilde,SME,SMEmunu, &
                                 rad_ir,radi,radr,radt,     &
                                 SubCont)
          SubTerm = SubTerm - SubCont
          if (show_terms) then
            write(*,"(A,3(I0,1x),G0)")                   &
              "StCirtCSirt,i,r,t,Term: ",radi,radr,radt, &
                                         SubCont*Cirs%term(iterm)%symfact
!            write(*,"(G0)") SubCont              & 
!                          * weightPS              &
!                          * fluxfact              &
!                          * alphas**(border_as+2) &
!                          * alphaEM**border_aEM   &
!                          * Cirs%term(iterm)%symfact
          end if
        end if
! Including prefactors not changing with scales:
        SubTerm = SubTerm  &
                * weightPS &
                * fluxfact &
                * ff       &
                * Cirs%term(iterm)%symfact
        SubTerm_multiscale = 0d0
        call StoreXis_scales
        do iscale=1,nscales
          call ChangeXis_scales(iscale)
          call calcmyscales(p)
          call calc_couplings
!
          subsRR(iscale) = subsRR(iscale)        &
                         - SubTerm               &
                         * cfunc                 &
                         * alphas**(border_as+2) &
                         * alphaEM**border_aEM
!
          SubTerm_multiscale(iscale) = SubTerm               & 
                                     * alphas**(border_as+2) &
                                     * alphaEM**border_aEM
        end do
! If an analysis is enabled and not being in optimization mode
! the analysis routine is called:
        if (.not.flg_optim.and.flg_analysis.and.(weightPS.ne.0)) then
          call analysis(ptilde,+SubTerm_multiscale)
        end if
! Restoring the old xi values:
        call RestoreXis_scales
      end do
    end do
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@ StSrt,StCirtSrt,StCSirtSrt,StCirtCSirtSrt @@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    do iterm=1,Srs%numterm
! Only terms defined for gluons can contribute:
      if (Srs%term(iterm)%rad(1)%ID.ne.0) cycle
! Only those Srs terms are stored where r < s, but for
! the St-type we have to have both cases:
      do rad1=1,2
        radr = Srs%term(iterm)%rad(rad1)%i
        radt = Srs%term(iterm)%rad(mod(rad1,2)+1)%i
        rad_r = radr
        if (radt.lt.radr) rad_r = rad_r - 1
        call MapMomSr_A12(p,phat,       &
                          radt,         &
                          Q,Q2,yiQ_arr, &
                          ytQ,lambda_t)
        call CalcSubhatInvariants(phat)
        call MapMomSr_A12(phat,ptilde,     &
                          rad_r,           &
                          Q,Q2,yihatQ_arr, &
                          yrhatQ,lambda_r)
! alphamax cut:
        ff = ffunc(y0,ytQ,dprfunc(size(p)-4))    &
           * ffunc(y0,yrhatQ,dprfunc(size(p)-4))
        if (ff.eq.0) cycle
! Applying cuts, if any and cycle if not passing:
        call Cuts(ptilde,cfunc)
        if (cfunc.eq.0d0) cycle
! Obtaining a flux factor:
        fluxfact = CalcFluxFact(ptilde)
! Obtaining the underlying Born doubly color-correlated SME:
        call CalcBijkl(ptilde,Bijkl)
! And deriving all the others from this:
        call CastSMEijklToSMEij(.false.,p,Bijkl,Bij)
        call CastSMEijToSME(p,Bij,SME)
!##############################################################
!########################### StSrt ############################
!##############################################################
        call CalcStSrtFF(p,phat,ptilde,Bij,Bijkl, &
                         rad_r,radr,radt,         &
                         SubTerm)
        if (show_terms) then
          write(*,"(A,2(I0,1x),G0)")      &
            "StSrt,r,t,Term: ",radr,radt, &
                               -SubTerm*Srs%term(iterm)%symfact
!          write(*,"(G0)") -SubTerm              & 
!                        * weightPS              &
!                        * fluxfact              &
!                        * alphas**(border_as+2) &
!                        * alphaEM**border_aEM   &
!                        * Srs%term(iterm)%symfact
        end if
! We have to pick up one additional parton, note that r and t
! are both gluons hence i can be anything (partonic):
        do radi=1,size(p)
          if (abs(p(radi)%flv).gt.qcd_nf) cycle
          if ((radi.eq.radr).or.(radi.eq.radt)) cycle
          rad_i = radi
          if (radt.lt.radi) rad_i = rad_i - 1
          rad_ir = min(rad_i,rad_r)
!##############################################################
!########################## StCirtSrt #########################
!##############################################################
          call CalcStCirtSrtFF(p,phat,ptilde,SME,          &
                               rad_i,rad_r,radi,radr,radt, &
                               SubCont)
          SubTerm = SubTerm - SubCont
          if (show_terms) then
            write(*,"(A,3(I0,1x),G0)")                 &
              "StCirtSrt,i,r,t,Term: ",radi,radr,radt, &
                                       SubCont*Srs%term(iterm)%symfact
!            write(*,"(G0)") SubCont              & 
!                          * weightPS              &
!                          * fluxfact              &
!                          * alphas**(border_as+2) &
!                          * alphaEM**border_aEM   &
!                          * Srs%term(iterm)%symfact
          end if
!##############################################################
!######################### StCSirtSrt #########################
!##############################################################
          call CalcStCSirtSrtFF(p,phat,ptilde,Bij,                 &
                                rad_ir,rad_i,rad_r,radi,radr,radt, &
                                SubCont)
          SubTerm = SubTerm - SubCont
          if (show_terms) then
            write(*,"(A,3(I0,1x),G0)")                  &
              "StCSirtSrt,i,r,t,Term: ",radi,radr,radt, &
                                        SubCont*Srs%term(iterm)%symfact
!            write(*,"(G0)") SubCont              & 
!                          * weightPS              &
!                          * fluxfact              &
!                          * alphas**(border_as+2) &
!                          * alphaEM**border_aEM   &
!                          * Srs%term(iterm)%symfact
          end if
!##############################################################
!####################### StCirtCSirtSrt #######################
!##############################################################
          call CalcStCirtCSirtSrtFF(p,phat,ptilde,SME,          &
                                    rad_i,rad_r,radi,radr,radt, &
                                    SubCont)
          SubTerm = SubTerm + SubCont
          if (show_terms) then
            write(*,"(A,3(I0,1x),G0)")                      &
              "StCirtCSirtSrt,i,r,t,Term: ",radi,radr,radt, &
                                            -SubCont*Srs%term(iterm)%symfact
!            write(*,"(G0)") -SubCont              & 
!                          * weightPS              &
!                          * fluxfact              &
!                          * alphas**(border_as+2) &
!                          * alphaEM**border_aEM   &
!                          * Srs%term(iterm)%symfact
          end if
        end do
! Including prefactors not changing with scales:
        SubTerm = SubTerm   &
                * weightPS  &
                * fluxfact  &
                * ff        &
                * Srs%term(iterm)%symfact
        SubTerm_multiscale = 0d0
        call StoreXis_scales
        do iscale=1,nscales
          call ChangeXis_scales(iscale)
          call calcmyscales(p)
          call calc_couplings
!
          subsRR(iscale) = subsRR(iscale)        &
                         - SubTerm               &
                         * cfunc                 &
                         * alphas**(border_as+2) &
                         * alphaEM**border_aEM
!
          SubTerm_multiscale(iscale) = SubTerm               & 
                                     * alphas**(border_as+2) &
                                     * alphaEM**border_aEM
        end do
! If an analysis is enabled and not being in optimization mode
! the analysis routine is called:
        if (.not.flg_optim.and.flg_analysis.and.(weightPS.ne.0)) then
          call analysis(ptilde,+SubTerm_multiscale)
        end if
! Restoring the old xi values:
        call RestoreXis_scales
      end do
    end do
!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!&&&&&&&&&&&&&&&&&&&&&&&& CktSt-type &&&&&&&&&&&&&&&&&&&&&&&&&&
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@ CitStCirt,CktStCSirt @@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    do iterm=1,Cirs%numterm
! Picking a candidate for t:
      do rad1=1,3
        radt = Cirs%term(iterm)%rad(mod(rad1+1,3)+1)%i
! Only gluons contribute:
        if (p(radt)%flv.ne.0) cycle
! Obtaining the other legs as well:
        call ReorderRadir(Cirs%term(iterm)%rad(mod(rad1+2,3)+1)%i,  &
                          Cirs%term(iterm)%rad(mod(rad1+2,3)+1)%ID, &
                          Cirs%term(iterm)%rad(mod(rad1+3,3)+1)%i,  &
                          Cirs%term(iterm)%rad(mod(rad1+3,3)+1)%ID, &
                          radi,radr)
        rad_i = radi
        rad_r = radr
        if (rad_i.gt.radt) rad_i = rad_i - 1
        if (rad_r.gt.radt) rad_r = rad_r - 1
        rad_ir = min(rad_i,rad_r)
        call MapMomSr_A12(p,phat,       &
                          radt,         &
                          Q,Q2,yiQ_arr, &
                          ytQ,lambda_t)
        call CalcSubhatInvariants(phat)
        if (rad_ir.lt.3) then
          print *,"ISR is not treated in Calc_mp2_subs..."
          stop "Calc_mp2_subs"
        else
          call MapMomCirFF_A12(phat,ptilde,                &
                               rad_i,rad_r,                &
                               Q,Q2,yihatQ_arr,sirhat_arr, &
                               sir,zhati,zhatr,            &
                               kthati,kthatr,              &
                               alphairt)
        end if
! alphamax cut:
        ff = ffunc(y0,ytQ,dprfunc(size(p)-4)) &
           * ffunc(alpha0,alphairt,dfunc(size(p)-4))
        if (ff.eq.0) cycle
! Applying cuts, if any and cycle if not passing:
        call Cuts(ptilde,cfunc)
        if (cfunc.eq.0d0) cycle
! Obtaining a flux factor:
        fluxfact = CalcFluxFact(ptilde)
! Note that because of the iterated PS mapping if a spin-
! correlated SME is needed it is calculated w.r.t the leg
! rad_ir:
        if (ptilde(rad_ir)%flv.eq.0) then
          call CalcBmunu(rad_ir,ptilde,SMEmunu)
          call CastSMEmunutoSME(ptilde,SMEmunu,SME)
        else
          call CalcB(ptilde,SME)
        end if
!##############################################################
!######################### CitStCirt ##########################
!##############################################################
        call CalcCitStCirtFF(p,phat,ptilde,SME,SMEmunu,  &
                             rad_i,rad_r,radi,radr,radt, &
                             SubTerm)
        if (show_terms) then
          write(*,"(A,3(I0,1x),G0)")                 &
            "CitStCirt,i,r,t,Term: ",radi,radr,radt, &
                                     SubTerm*Cirs%term(iterm)%symfact
!          write(*,"(G0)") SubTerm              & 
!                        * weightPS              &
!                        * fluxfact              &
!                        * alphas**(border_as+2) &
!                        * alphaEM**border_aEM   &
!                        * Cirs%term(iterm)%symfact
        end if
!##############################################################
!######################### CrtStCirt ##########################
!##############################################################
        call CalcCrtStCirtFF(p,phat,ptilde,SME,SMEmunu,  &
                             rad_i,rad_r,radi,radr,radt, &
                             SubCont)
        SubTerm = SubTerm + SubCont
        if (show_terms) then
          write(*,"(A,3(I0,1x),G0)")                 &
            "CrtStCirt,i,r,t,Term: ",radi,radr,radt, &
                                     SubCont*Cirs%term(iterm)%symfact
!          write(*,"(G0)") SubCont              & 
!                        * weightPS              &
!                        * fluxfact              &
!                        * alphas**(border_as+2) &
!                        * alphaEM**border_aEM   &
!                        * Cirs%term(iterm)%symfact
        end if
! To include the contribution from CktStCSirt terms a new
! leg has to be introduced:
        do radk=1,size(p)
! Only QCD partons:
          if (abs(p(radk)%flv).gt.qcd_nf) cycle
! Cannot coincide with the other selected legs:
          if ((radk.eq.radi).or.(radk.eq.radr).or.(radk.eq.radt)) cycle
!##############################################################
!######################### CktStCSirt #########################
!##############################################################
          call CalcCktStCSirtFF(p,phat,ptilde,SME,SMEmunu, &
                                radi,radr,radk,radt,       &
                                SubCont)
          if (show_terms) then
            write(*,"(A,4(I0,1x),G0)")                         &
              "CktStCSirt,k,i,r,t,Term: ",radk,radi,radr,radt, &
                                          SubCont*Cirs%term(iterm)%symfact
!            write(*,"(G0)") SubCont              & 
!                          * weightPS              &
!                          * fluxfact              &
!                          * alphas**(border_as+2) &
!                          * alphaEM**border_aEM   &
!                          * Cirs%term(iterm)%symfact
          end if
          SubTerm = SubTerm + SubCont
        end do
! Including prefactors not changing with scales:
        SubTerm = SubTerm  &
                * weightPS &
                * fluxfact &
                * ff       &
                * Cirs%term(iterm)%symfact
        SubTerm_multiscale = 0d0
        call StoreXis_scales
        do iscale=1,nscales
          call ChangeXis_scales(iscale)
          call calcmyscales(p)
          call calc_couplings
!
          subsRR(iscale) = subsRR(iscale)        &
                         + SubTerm               &
                         * cfunc                 &
                         * alphas**(border_as+2) &
                         * alphaEM**border_aEM
!
          SubTerm_multiscale(iscale) = SubTerm               & 
                                     * alphas**(border_as+2) &
                                     * alphaEM**border_aEM
        end do
! If an analysis is enabled and not being in optimization mode
! the analysis routine is called:
        if (.not.flg_optim.and.flg_analysis.and.(weightPS.ne.0)) then
          call analysis(ptilde,-SubTerm_multiscale)
        end if
! Restoring the old xi values:
        call RestoreXis_scales
      end do
    end do
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@ CrtStSrt @@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    do iterm=1,Srs%numterm
! Only terms defined for gluons can contribute:
      if (Srs%term(iterm)%rad(1)%ID.ne.0) cycle
! Only those Srs terms are stored where r < s, but for
! the St-type we have to have both cases:
      do rad1=1,2
        radr = Srs%term(iterm)%rad(rad1)%i
        radt = Srs%term(iterm)%rad(mod(rad1,2)+1)%i
        rad_r = radr
        if (rad_r.gt.radt) rad_r = rad_r - 1
        call MapMomSr_A12(p,phat,       &
                          radt,         &
                          Q,Q2,yiQ_arr, &
                          ytQ,lambda_t)
        call CalcSubhatInvariants(phat)
        call MapMomSr_A12(phat,ptilde,     &
                          rad_r,           &
                          Q,Q2,yihatQ_arr, &
                          yrhatQ,lambda_r)
! alphamax cut:
        ff = ffunc(y0,ytQ,dprfunc(size(p)-4)) &
           * ffunc(y0,yrhatQ,dprfunc(size(p)-4))
        if (ff.eq.0) cycle
! Applying cuts, if any and cycle if not passing:
        call Cuts(ptilde,cfunc)
        if (cfunc.eq.0d0) cycle
! Obtaining a flux factor:
        fluxfact = CalcFluxFact(ptilde)
        call CalcBij(ptilde,Bij)
        call CastSMEijToSME(ptilde,Bij,SME)
!##############################################################
!########################## CrtStSrt ##########################
!##############################################################
        call CalcCrtStSrtFF(p,phat,ptilde,Bij, &
                            rad_r,radr,radt,   &
                            SubTerm)
        if (show_terms) then
          write(*,"(A,2(I0,1x),G0)")         &
            "CrtStSrt,r,t,Term: ",radr,radt, &
                                  SubTerm*Srs%term(iterm)%symfact
!          write(*,"(G0)") SubTerm              & 
!                        * weightPS              &
!                        * fluxfact              &
!                        * alphas**(border_as+2) &
!                        * alphaEM**border_aEM   &
!                        * Srs%term(iterm)%symfact
        end if
        do radk=1,size(p)
! Only QCD partons:
          if (abs(p(radk)%flv).gt.qcd_nf) cycle
! Cannot coincide with the other selected legs:
          if ((radk.eq.radr).or.(radk.eq.radt)) cycle
          rad_k = radk
          if (rad_k.gt.radt) rad_k = rad_k - 1
!##############################################################
!######################## CktStCkrtSrt ########################
!##############################################################
          call CalcCktStCkrtSrtFF(p,phat,ptilde,SME, &
                                  rad_k,rad_r,       &
                                  radk,radr,radt,    &
                                  SubCont)
          if (show_terms) then
            write(*,"(A,3(I0,1x),G0)")                    &
              "CktStCkrtSrt,k,r,t,Term: ",radk,radr,radt, &
                                          -SubCont*Srs%term(iterm)%symfact
!            write(*,"(G0)") -SubCont              & 
!                          * weightPS              &
!                          * fluxfact              &
!                          * alphas**(border_as+2) &
!                          * alphaEM**border_aEM   &
!                          * Srs%term(iterm)%symfact
          end if
          SubTerm = SubTerm - SubCont
!##############################################################
!######################## CrtStCkrtSrt ########################
!##############################################################
          call CalcCrtStCkrtSrtFF(p,phat,ptilde,SME, &
                                  rad_k,rad_r,       &
                                  radk,radr,radt,    &
                                  SubCont)
          if (show_terms) then
            write(*,"(A,3(I0,1x),G0)")                    &
              "CrtStCkrtSrt,k,r,t,Term: ",radk,radr,radt, &
                                          -SubCont*Srs%term(iterm)%symfact
!            write(*,"(G0)") -SubCont              & 
!                          * weightPS              &
!                          * fluxfact              &
!                          * alphas**(border_as+2) &
!                          * alphaEM**border_aEM   &
!                          * Srs%term(iterm)%symfact
          end if
          SubTerm = SubTerm - SubCont
!##############################################################
!########################## CktStSrt ##########################
!##############################################################
          call CalcCktStSrtFF(p,phat,ptilde,Bij,    &
                              rad_r,radk,radr,radt, &
                              SubCont)
          if (show_terms) then
            write(*,"(A,3(I0,1x),G0)")                &
              "CktStSrt,k,r,t,Term: ",radk,radr,radt, &
                                      SubCont*Srs%term(iterm)%symfact
!            write(*,"(G0)") SubCont              & 
!                          * weightPS              &
!                          * fluxfact              &
!                          * alphas**(border_as+2) &
!                          * alphaEM**border_aEM   &
!                          * Srs%term(iterm)%symfact
          end if
          SubTerm = SubTerm + SubCont
          do radi=1,size(p)
! Only QCD partons:
            if (abs(p(radi)%flv).gt.qcd_nf) cycle
! Cannot coincide with the other selected legs:
            if ((radi.eq.radk).or.(radi.eq.radr).or.(radi.eq.radt)) cycle
            rad_i = radi
            if (rad_i.gt.radt) rad_i = rad_i - 1
!##############################################################
!####################### CktStCSirtSrt ########################
!##############################################################
            call CalcCktStCSirtSrtFF(p,phat,ptilde,SME,   &
                                     rad_i,rad_r,         &
                                     radi,radr,radk,radt, &
                                     SubCont)
            if (show_terms) then
              write(*,"(A,4(I0,1x),G0)")                            &
                "CktStCSirtSrt,k,i,r,t,Term: ",radk,radi,radr,radt, &
                                               -SubCont*Srs%term(iterm)%symfact
!              write(*,"(G0)") -SubCont              & 
!                            * weightPS              &
!                            * fluxfact              &
!                            * alphas**(border_as+2) &
!                            * alphaEM**border_aEM   &
!                            * Srs%term(iterm)%symfact
            end if
            SubTerm = SubTerm - SubCont
          end do
        end do
! Including prefactors not changing with scales:
        SubTerm = SubTerm  &
                * weightPS &
                * fluxfact &
                * ff       &
                * Srs%term(iterm)%symfact
        SubTerm_multiscale = 0d0
        call StoreXis_scales
        do iscale=1,nscales
          call ChangeXis_scales(iscale)
          call calcmyscales(p)
          call calc_couplings
!
          subsRR(iscale) = subsRR(iscale)        &
                         + SubTerm               &
                         * cfunc                 &
                         * alphas**(border_as+2) &
                         * alphaEM**border_aEM
!
          SubTerm_multiscale(iscale) = SubTerm               & 
                                     * alphas**(border_as+2) &
                                     * alphaEM**border_aEM
        end do
! If an analysis is enabled and not being in optimization mode
! the analysis routine is called:
        if (.not.flg_optim.and.flg_analysis.and.(weightPS.ne.0)) then
          call analysis(ptilde,-SubTerm_multiscale)
        end if
! Restoring the old xi values:
        call RestoreXis_scales
      end do
    end do
  end if
!
end subroutine Calc_mp2_subs
!
! This momentum mapping is essentially the same as of MapMomCirFF
! with the exception of constructing the transverse momentum.
! In this definition of mapping Zetair is chosen to be zero,
! hence corresponding to the mapping defined in arXiv:hep-ph/0609043:
subroutine MapMomCir0FF(p,ptilde,emitir,emitID,radi,radiID,radr,radrID)
use particles
use regions
use math
use nlo_mom_maps
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
  type(particle) , dimension(:) , intent(out) :: ptilde
  integer , intent(in) :: emitir,radi,radr
  integer , intent(in) :: emitID,radiID,radrID
!
  integer :: ipart,emit,rad
!
!
!
!
! Obtaining y variables:
  yiQ  = yiQ_arr(radi)
  yrQ  = yiQ_arr(radr)
  yirQ = yiQ + yrQ
  sir  = sir_arr(radi,radr)
  yir  = sir/Q2
  ztildei = yiQ/yirQ
  ztilder = 1d0 - ztildei
  alphair = 0.5d0*(yirQ - sqrt(yirQ**2 - 4d0*yir))
! The new momentum:
  ptildeir%p = (p(radi)%p + p(radr)%p - alphair*Q)/(1d0 - alphair)
!
! We calculate the transverse momentum, too:
  zetai_r    = ztildei - yir/(alphair*yirQ)
  zetar_i    = ztilder - yir/(alphair*yirQ)
  yirtildeQ  = 2d0*ptildeir%p*Q/Q2
  zetair     = 0d0
  kttildei%p = zetai_r*p(radr)%p - zetar_i*p(radi)%p &
             + zetair*ptildeir%p
  kttildei%flv = emitID
! Safety check upon branching:
  if ((radiID+radrID).ne.emitID) then
    print *,"Problem in MapMomCir0FF..."
    write(*,'(a,3(1x,I0))') "emitir,radi,radr: ",emitir,radi,radr
    write(*,'(a,3(1x,I0))') "emitID,radiID,radrID: ",emitID,radiID,radrID
  end if
  ptildeir%flv = p(radi)%flv + p(radr)%flv
! it is possible that radi > radr, the radiated parton should always
! have a higher position, hence:
  emit = min(radi,radr)
  rad  = max(radi,radr)
  call RemapMomentaO1(p,ptilde,1d0/(1d0 - alphair),emit,ptildeir,rad)
!
!
end subroutine MapMomCir0FF
!
subroutine PickCir(p,ptilde,emit,rad,Cir,CalcSMEB,CalcSMEBmunu, &
                   Cirterm)
use regions
use particles
use nlo_subtractions
use nlo_mom_maps
implicit none
!
  type(particle), dimension(:), intent(in) :: p
  type(particle), dimension(:), intent(out) :: ptilde
  integer , intent(in) :: emit,rad
  type(subterms) , intent(in) :: Cir
  real(kind(1d0)) , intent(out) :: Cirterm
!
  integer :: iterm
  real(kind(1d0)) :: smeB
  real(kind(1d0)) , dimension(0:3,0:3) :: Bmunu
  real(kind(1d0)) :: CirCont
!
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
    subroutine CalcSMEBmunu(ileg,p,Bmunu)
    use particles
    implicit none
!
      integer , intent(in) :: ileg
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(0:3,0:3) , intent(out) :: Bmunu
!
    end subroutine CalcSMEBmunu
  end interface
!
  Cirterm = 0d0
! We go through all the subtraction terms and pick those which
! are contributing:
  do iterm=1,Cir%numterm
    if ((Cir%term(iterm)%rad(1)%i.eq.emit).and.  &
        (Cir%term(iterm)%rad(2)%i.eq.rad)) then
      call MapMomCirFF(p,ptilde, &
                       Cir%term(iterm)%emit(1)%i,  &
                       Cir%term(iterm)%emit(1)%ID, &
                       Cir%term(iterm)%rad(1)%i,   &
                       Cir%term(iterm)%rad(1)%ID,  &
                       Cir%term(iterm)%rad(2)%i,   &
                       Cir%term(iterm)%rad(2)%ID)
! Obtaining the underlying Born SME:
! If the emitter is a gluon we always need spin-correlated SME:
      if (Cir%term(iterm)%emit(1)%ID.eq.0) then
        call CalcSMEBmunu(Cir%term(iterm)%emit(1)%i,ptilde,Bmunu)
! Otherwise just the SME is requested:
      else
        call CalcSMEB(ptilde,smeB)
      end if
      call CalcCirFF(p,ptilde,smeB,Bmunu,       &
                     Cir%term(iterm)%emit(1)%i, &
                     Cir%term(iterm)%rad(1)%i,  &
                     Cir%term(iterm)%rad(2)%i,  &
                     CirCont)
      Cirterm = Cirterm + Cir%term(iterm)%symfact*CirCont
    end if
  end do
!
end subroutine PickCir
!
subroutine PickSr(p,ptilde,Bij,radr,Sr,CalcSMEBij,Srterm)
use regions
use particles
use nlo_subtractions
use nlo_mom_maps
implicit none
!
  type(particle), dimension(:), intent(in) :: p
  type(particle), dimension(:), intent(out) :: ptilde
  real(kind(1d0)) , dimension(:,:) , intent(out) :: Bij
  integer , intent(in) :: radr
  type(subterms) , intent(in) :: Sr
  real(kind(1d0)) , intent(out) :: Srterm
!
  integer :: iterm
  real(kind(1d0)) :: SrCont
!
!
  interface 
    subroutine CalcSMEBij(p,Bij)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Bij
!
    end subroutine CalcSMEBij
  end interface
!
  Srterm = 0d0
! We go through all the subtraction terms and pick those which
! are contributing:
  do iterm=1,Sr%numterm
    if (Sr%term(iterm)%rad(1)%i.eq.radr) then
      call MapMomSr(p,ptilde, &
                    Sr%term(iterm)%rad(1)%i, &
                    Sr%term(iterm)%rad(1)%ID)
      call CalcSMEBij(ptilde,Bij)
      call CalcSr(p,ptilde,Bij, &
                  Sr%term(iterm)%rad(1)%i, &
                  SrCont)
      Srterm = Srterm + Sr%term(iterm)%symfact*SrCont
    end if
  end do
!
end subroutine PickSr
!
subroutine PickCirSr(p,ptilde,emit,rad,CSir,CalcSMEB,CSirterm)
use regions
use particles
use nlo_mom_maps
use nlo_subtractions
implicit none
!
  type(particle), dimension(:), intent(in) :: p
  type(particle), dimension(:), intent(out) :: ptilde
  integer , intent(in) :: emit,rad
  type(subterms) , intent(in) :: CSir
  real(kind(1d0)) , intent(out) :: CSirterm
!
  integer :: iterm
  real(kind(1d0)) :: smeB
  real(kind(1d0)) :: CSirCont
!
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
  end interface
!
  CSirterm = 0d0
! We go through all the subtraction terms and pick those which
! are contributing:
  do iterm=1,CSir%numterm
    if ((CSir%term(iterm)%rad(1)%i.eq.emit).and. &
        (CSir%term(iterm)%rad(2)%i.eq.rad)) then
      call MapMomSr(p,ptilde, &
                    CSir%term(iterm)%rad(2)%i, &
                    CSir%term(iterm)%rad(2)%ID)
      call CalcSMEB(ptilde,smeB)
      if (CSir%term(iterm)%emit(1)%i.gt.2) then
        call CalcCirFFSr(p,ptilde,smeB, &
                         CSir%term(iterm)%emit(1)%i, &
                         CSir%term(iterm)%rad(1)%i,  &
                         CSir%term(iterm)%rad(2)%i,  &
                         CSirCont)
      end if
      CSirterm = CSirterm + CSir%term(iterm)%symfact*CSirCont
    end if
  end do
!
end subroutine PickCirSr
!
subroutine PickCirRV(p,ptilde,emit,rad,Cir,Cirterm)
use regions
use particles
use nnlo_subtractions
implicit none
!
  type(particle), dimension(:), intent(in) :: p
  type(particle), dimension(:), intent(out) :: ptilde
  integer , intent(in) :: emit,rad
  type(subterms) , intent(in) :: Cir
  real(kind(1d0)) , intent(out) :: Cirterm
!
  integer :: iterm
  real(kind(1d0)) :: smeB
  real(kind(1d0)) , dimension(0:3,0:3) :: Bmunu
  real(kind(1d0)) :: CirCont
!
!
  interface 
    subroutine MapMomCir0FF(p,ptilde,emit,emitID,radi,radiID, &
                            radr,radrID)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      type(particle) , dimension(:) , intent(out) :: ptilde
      integer , intent(in) :: emit,radi,radr
      integer , intent(in) :: emitID,radiID,radrID
!
    end subroutine MapMomCir0FF
!
    subroutine CalcB(p,smeB)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: smeB
!
    end subroutine CalcB
!
    subroutine CalcBmunu(ileg,p,Bmunu)
    use particles
    implicit none
!
      integer , intent(in) :: ileg
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(0:3,0:3) , intent(out) :: Bmunu
!
    end subroutine CalcBmunu
  end interface
!
  Cirterm = 0d0
! We go through all the subtraction terms and pick those which
! are contributing:
  do iterm=1,Cir%numterm
    if ((Cir%term(iterm)%rad(1)%i.eq.emit).and.  &
        (Cir%term(iterm)%rad(2)%i.eq.rad)) then
      call MapMomCir0FF(p,ptilde, &
                        Cir%term(iterm)%emit(1)%i,  &
                        Cir%term(iterm)%emit(1)%ID, &
                        Cir%term(iterm)%rad(1)%i,   &
                        Cir%term(iterm)%rad(1)%ID,  &
                        Cir%term(iterm)%rad(2)%i,   &
                        Cir%term(iterm)%rad(2)%ID)
      if (Cir%term(iterm)%emit(1)%ID.eq.0) then
        call CalcBmunu(Cir%term(iterm)%emit(1)%i,ptilde,Bmunu)
      else
        call CalcB(ptilde,smeB)
      end if
      call CalcCir01FF(p,ptilde, &
                       Cir%term(iterm)%emit(1)%i,  &
                       Cir%term(iterm)%rad(1)%i,   &
                       Cir%term(iterm)%rad(2)%i,   &
                       CirCont)
      Cirterm = Cirterm + Cir%term(iterm)%symfact*CirCont
      call CalcCir10FF(p,ptilde,smeB,Bmunu, &
                       Cir%term(iterm)%emit(1)%i,  &
                       Cir%term(iterm)%rad(1)%i,   &
                       Cir%term(iterm)%rad(2)%i,   &
                       CirCont)
      Cirterm = Cirterm + Cir%term(iterm)%symfact*CirCont
    end if
  end do
!
end subroutine PickCirRV
!
subroutine PickSrRV(p,ptilde,Bij,Vij,Bijk,radr,Sr,Srterm)
use regions
use particles
use nlo_mom_maps
use nnlo_subtractions
implicit none
!
  type(particle), dimension(:), intent(in) :: p
  type(particle), dimension(:), intent(out) :: ptilde
  real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
  real(kind(1d0)) , dimension(:,:) , intent(inout) :: Vij
  real(kind(1d0)) , dimension(:,:,:) , intent(inout) :: Bijk
  integer , intent(in) :: radr
  type(subterms) , intent(in) :: Sr
  real(kind(1d0)) , intent(out) :: Srterm
!
  integer :: iterm
  real(kind(1d0)) :: SrCont
!
!
  interface 
    subroutine CalcBij(p,Bij)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Bij
!
    end subroutine CalcBij
!
    subroutine CalcVij(p,Vij,VijLaurent)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Vij
      real(kind(1d0)) , optional , dimension(:,:,-4:) , intent(out) :: VijLaurent
!
    end subroutine CalcVij
  end interface
!
  Srterm = 0d0
! We go through all the subtraction terms and pick those which
! are contributing:
  do iterm=1,Sr%numterm
    if (Sr%term(iterm)%rad(1)%i.eq.radr) then
      call MapMomSr(p,ptilde, &
                    Sr%term(iterm)%rad(1)%i, &
                    Sr%term(iterm)%rad(1)%ID)
!
      call CalcBij(ptilde,Bij)
      call CalcVij(ptilde,Vij)
!
      call CalcSr01(p,ptilde,Vij, &
                    Sr%term(iterm)%rad(1)%i,   &
                    SrCont)
      Srterm = Srterm + Sr%term(iterm)%symfact*SrCont
      call CalcSr10(p,ptilde,Bij,Bijk, &
                    Sr%term(iterm)%rad(1)%i,   &
                    SrCont)
      Srterm = Srterm + Sr%term(iterm)%symfact*SrCont
    end if
  end do
!
end subroutine PickSrRV
!
subroutine PickCirSrRV(p,ptilde,emit,rad,CirSr,CirSrterm)
use regions
use particles
use nlo_mom_maps
use nnlo_subtractions
implicit none
!
  type(particle), dimension(:), intent(in) :: p
  type(particle), dimension(:), intent(out) :: ptilde
  integer , intent(in) :: emit,rad
  type(subterms) , intent(in) :: CirSr
  real(kind(1d0)) , intent(out) :: CirSrterm
!
  integer :: iterm
  real(kind(1d0)) :: smeB,smeV
  real(kind(1d0)) :: CirSrCont
!
!
  interface 
    subroutine CalcB(p,smeB)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: smeB
!
    end subroutine CalcB
!
    subroutine CalcV(p,smeV,smeVLaurent)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: smeV
      real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: smeVLaurent
!
    end subroutine CalcV
  end interface
!
  CirSrterm = 0d0
! We go through all the subtraction terms and pick those which
! are contributing:
  do iterm=1,CirSr%numterm
    if ((CirSr%term(iterm)%rad(1)%i.eq.emit).and. &
        (CirSr%term(iterm)%rad(2)%i.eq.rad)) then
      call MapMomSr(p,ptilde, &
                    CirSr%term(iterm)%rad(2)%i, &
                    CirSr%term(iterm)%rad(2)%ID)
!
      call CalcB(ptilde,smeB)
      call CalcV(ptilde,smeV)
!
      if (CirSr%term(iterm)%emit(1)%i.gt.2) then
        call CalcCirFFSr01(p,ptilde,smeV, &
                           CirSr%term(iterm)%emit(1)%i, &
                           CirSr%term(iterm)%rad(1)%i,  &
                           CirSr%term(iterm)%rad(2)%i,  &
                           CirSrCont)
        CirSrterm = CirSrterm + CirSr%term(iterm)%symfact*CirSrCont
        call CalcCirFFSr10(p,ptilde,smeB, &
                           CirSr%term(iterm)%emit(1)%i, &
                           CirSr%term(iterm)%rad(1)%i,  &
                           CirSr%term(iterm)%rad(2)%i,  &
                           CirSrCont)
        CirSrterm = CirSrterm + CirSr%term(iterm)%symfact*CirSrCont
      else
        print *,"Error occured in PickCirSrRV..."
        print *,"Initial state cannot be treated yet..."
        stop
      end if
    end if
  end do
!
end subroutine PickCirSrRV
!
subroutine PickCirRRA1(p,ptilde,Bij,Bmunuij,emit,rad,Cir,Cirterm)
use regions
use particles
use nnlo_subtractions
implicit none
!
  type(particle), dimension(:), intent(in) :: p
  type(particle), dimension(:), intent(out) :: ptilde
  real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
  real(kind(1d0)) , dimension(0:,0:,:,:) , intent(inout) :: Bmunuij
  integer , intent(in) :: emit,rad
  type(subterms) , intent(in) :: Cir
  real(kind(1d0)) , intent(out) :: Cirterm
!
  integer :: iterm
  real(kind(1d0)) :: smeB
  real(kind(1d0)) , dimension(0:3,0:3) :: Bmunu
  real(kind(1d0)) :: CirCont
!
!
  interface 
    subroutine MapMomCir0FF(p,ptilde,emit,emitID,radi,radiID, &
                            radr,radrID)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      type(particle) , dimension(:) , intent(out) :: ptilde
      integer , intent(in) :: emit,radi,radr
      integer , intent(in) :: emitID,radiID,radrID
!
    end subroutine MapMomCir0FF
!
    subroutine CalcSubtildeInvariants(p)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
!
    end subroutine CalcSubtildeInvariants
!
    subroutine CalcB(p,smeB)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: smeB
!
    end subroutine CalcB
!
    subroutine CalcBmunu(ileg,p,Bmunu)
    use particles
    implicit none
!
      integer , intent(in) :: ileg
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(0:3,0:3) , intent(out) :: Bmunu
!
    end subroutine CalcBmunu
!
    subroutine CalcBij(p,Bij)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Bij
!
    end subroutine CalcBij
!
    subroutine CalcBmunuij(ileg,p,Bmunuij)
    use particles
    implicit none
!
      integer , intent(in) :: ileg
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(0:,0:,:,:) , intent(out) :: Bmunuij
!
    end subroutine CalcBmunuij
  end interface
!
  Cirterm = 0d0
! We go through all the subtraction terms and pick those which
! are contributing:
  do iterm=1,Cir%numterm
    if ((Cir%term(iterm)%rad(1)%i.eq.emit).and.  &
        (Cir%term(iterm)%rad(2)%i.eq.rad)) then
      call MapMomCir0FF(p,ptilde, &
                        Cir%term(iterm)%emit(1)%i,  &
                        Cir%term(iterm)%emit(1)%ID, &
                        Cir%term(iterm)%rad(1)%i,   &
                        Cir%term(iterm)%rad(1)%ID,  &
                        Cir%term(iterm)%rad(2)%i,   &
                        Cir%term(iterm)%rad(2)%ID)
! Additional invariants are needed:
      call CalcSubtildeInvariants(ptilde)
!
      call CalcB(ptilde,smeB)
      call CalcBij(ptilde,Bij)
      if (Cir%term(iterm)%emit(1)%ID.eq.0) then
        call CalcBmunu(Cir%term(iterm)%emit(1)%i,ptilde,Bmunu)
        call CalcBmunuij(Cir%term(iterm)%emit(1)%i,ptilde,Bmunuij)
      end if
!
      call CalcCir00IFF(p,ptilde,smeB,Bij,Bmunu,Bmunuij, &
                       Cir%term(iterm)%emit(1)%i,        &
                       Cir%term(iterm)%rad(1)%i,         &
                       Cir%term(iterm)%rad(2)%i,         &
                       CirCont)
      Cirterm = Cirterm + Cir%term(iterm)%symfact*CirCont
      call CalcCirR00FF(p,ptilde,smeB,Bmunu, &
                       Cir%term(iterm)%emit(1)%i,  &
                       Cir%term(iterm)%rad(1)%i,   &
                       Cir%term(iterm)%rad(2)%i,   &
                       CirCont)
      Cirterm = Cirterm + Cir%term(iterm)%symfact*CirCont
    end if
  end do
!
end subroutine PickCirRRA1
!
subroutine PickSrRRA1(p,ptilde,Bij,Bijkl,radr,Sr,Srterm)
use regions
use particles
use nlo_mom_maps
use nnlo_subtractions
implicit none
!
  type(particle), dimension(:), intent(in) :: p
  type(particle), dimension(:), intent(out) :: ptilde
  real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
  real(kind(1d0)) , dimension(:,:,:,:) , intent(inout) :: Bijkl
  integer , intent(in) :: radr
  type(subterms) , intent(in) :: Sr
  real(kind(1d0)) , intent(out) :: Srterm
!
  integer :: iterm
  real(kind(1d0)) :: SrCont
!
!
  interface 
    subroutine CalcSubtildeInvariants(p)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
!
    end subroutine CalcSubtildeInvariants
!
    subroutine CalcBij(p,Bij)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Bij
!
    end subroutine CalcBij
!
    subroutine CalcBijkl(p,Bijkl)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:,:,:,:) , intent(out) :: Bijkl
!
    end subroutine CalcBijkl
  end interface
!
  Srterm = 0d0
! We go through all the subtraction terms and pick those which
! are contributing:
  do iterm=1,Sr%numterm
    if (Sr%term(iterm)%rad(1)%i.eq.radr) then
      call MapMomSr(p,ptilde, &
                    Sr%term(iterm)%rad(1)%i,   &
                    Sr%term(iterm)%rad(1)%ID)
! Additional invariants are needed:
      call CalcSubtildeInvariants(ptilde)
!
      call CalcBij(ptilde,Bij)
      call CalcBijkl(ptilde,Bijkl)
!
      call CalcSr00I(p,ptilde,Bij,Bijkl, &
                     Sr%term(iterm)%rad(1)%i,   &
                     SrCont)
      Srterm = Srterm + Sr%term(iterm)%symfact*SrCont
      call CalcSrR00(p,ptilde,Bij, &
                     Sr%term(iterm)%rad(1)%i,   &
                     SrCont)
      Srterm = Srterm + Sr%term(iterm)%symfact*SrCont
    end if
  end do
!
end subroutine PickSrRRA1
!
subroutine PickCirSrRRA1(p,ptilde,Bij,emit,rad,CirSr,CirSrterm)
use misc
use regions
use particles
use nlo_mom_maps
use nnlo_subtractions
implicit none
!
  type(particle), dimension(:), intent(in) :: p
  type(particle), dimension(:), intent(out) :: ptilde
  real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
  integer , intent(in) :: emit,rad
  type(subterms) , intent(in) :: CirSr
  real(kind(1d0)) , intent(out) :: CirSrterm
!
  integer :: iterm
  real(kind(1d0)) :: CirSrCont
  real(kind(1d0)) :: smeB
!
!
  interface 
    subroutine CalcBij(p,Bij)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Bij
!
    end subroutine CalcBij
!
    subroutine CalcSubtildeInvariants(p)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
!
    end subroutine CalcSubtildeInvariants
  end interface
!
  CirSrterm = 0d0
! We go through all the subtraction terms and pick those which
! are contributing:
  do iterm=1,CirSr%numterm
    if ((CirSr%term(iterm)%rad(1)%i.eq.emit).and.  &
        (CirSr%term(iterm)%rad(2)%i.eq.rad)) then
      call MapMomSr(p,ptilde, &
                    CirSr%term(iterm)%rad(2)%i,   &
                    CirSr%term(iterm)%rad(2)%ID)
! Additional invariants are needed:
      call CalcSubtildeInvariants(ptilde)
! The Born SME is calculated for once and for all:
      call CalcBij(ptilde,Bij)
      call CastSMEijToSME(ptilde,Bij,smeB)
      call CalcCirFFSr00I(p,ptilde,smeB,Bij, &
                          CirSr%term(iterm)%emit(1)%i,  &
                          CirSr%term(iterm)%rad(1)%i,   &
                          CirSr%term(iterm)%rad(2)%i,   &
                          CirSrCont)
      CirSrterm = CirSrterm + CirSr%term(iterm)%symfact*CirSrCont
      call CalcCirFFSrR00(p,ptilde,smeB, &
                          CirSr%term(iterm)%emit(1)%i,  &
                          CirSr%term(iterm)%rad(1)%i,   &
                          CirSr%term(iterm)%rad(2)%i,   &
                          CirSrCont)
      CirSrterm = CirSrterm + CirSr%term(iterm)%symfact*CirSrCont
    end if
  end do
!
end subroutine PickCirSrRRA1
!
subroutine PickCirs(p,ptilde,emit,radi,radr,rads,Cirs, &
                    Cirsterm)
use regions
use particles
use misc
use nnlo_subtractions
use nnlo_mom_maps
implicit none
!
  type(particle), dimension(:), intent(in) :: p
  type(particle), dimension(:), intent(out) :: ptilde
  integer , intent(in) :: emit,radi,radr,rads
  type(subterms) , intent(in) :: Cirs
  real(kind(1d0)) , intent(out) :: Cirsterm
!
  integer :: iterm
  real(kind(1d0)) :: smeB
  real(kind(1d0)) , dimension(0:3,0:3) :: Bmunu
  real(kind(1d0)) :: CirsCont
!
!
  interface 
    subroutine CalcB(p,smeB)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: smeB
!
    end subroutine CalcB
!
    subroutine CalcBmunu(ileg,p,Bmunu)
    use particles
    implicit none
!
      integer , intent(in) :: ileg
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(0:3,0:3) , intent(out) :: Bmunu
!
    end subroutine CalcBmunu
!
    subroutine CalcBmunu_q1q2(ileg,q1,q2,p,Bmunu)
    use particles
    implicit none
!
      integer , intent(in) :: ileg
      type(mom) , intent(in) :: q1,q2
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(0:3,0:3) , intent(out) :: Bmunu
!
    end subroutine CalcBmunu_q1q2
  end interface
!
  Cirsterm = 0d0
! We go through all the subtraction terms and pick those which
! are contributing:
  do iterm=1,Cirs%numterm
    if ((Cirs%term(iterm)%emit(1)%i.eq.emit).and. &
        (Cirs%term(iterm)%rad(1)%i.eq.radi).and.  &
        (Cirs%term(iterm)%rad(2)%i.eq.radr).and. &
        (Cirs%term(iterm)%rad(3)%i.eq.rads)) then
      call MapMomCirsFF(p,ptilde, &
                        Cirs%term(iterm)%emit(1)%i,  &
                        Cirs%term(iterm)%emit(1)%ID, &
                        Cirs%term(iterm)%rad(1)%i,   &
                        Cirs%term(iterm)%rad(1)%ID,  &
                        Cirs%term(iterm)%rad(2)%i,   &
                        Cirs%term(iterm)%rad(2)%ID,  &
                        Cirs%term(iterm)%rad(3)%i,   &
                        Cirs%term(iterm)%rad(3)%ID)
! Obtaining the underlying Born SME:
! If the emitter is a gluon we always need spin-correlated SME:
      if (Cirs%term(iterm)%emit(1)%ID.eq.0) then
        call CalcBmunu(Cirs%term(iterm)%emit(1)%i,ptilde,Bmunu)
        call CastSMEmunuToSME(ptilde,Bmunu,smeB)
! Otherwise just the SME is requested:
      else
        call CalcB(ptilde,smeB)
      end if
      call CalcCirsFF(p,ptilde,smeB,Bmunu,       &
                      Cirs%term(iterm)%emit(1)%i, &
                      Cirs%term(iterm)%rad(1)%i,  &
                      Cirs%term(iterm)%rad(2)%i,  &
                      Cirs%term(iterm)%rad(3)%i,  &
                      CirsCont)
      Cirsterm = Cirsterm + Cirs%term(iterm)%symfact*CirsCont
    end if
  end do
!
end subroutine PickCirs
!
subroutine PickCirjs(p,ptilde,emiti,radi,radr,emitj,radj,rads, &
                     Cirjs,Cirjsterm)
use regions
use particles
use misc
use nnlo_subtractions
use nnlo_mom_maps
implicit none
!
  type(particle), dimension(:), intent(in) :: p
  type(particle), dimension(:), intent(out) :: ptilde
  integer , intent(in) :: emiti,emitj,radi,radr,radj,rads
  type(subterms) , intent(in) :: Cirjs
  real(kind(1d0)) , intent(out) :: Cirjsterm
!
  integer :: iterm
  real(kind(1d0)) :: smeB
  real(kind(1d0)) , dimension(0:3,0:3) :: Bmunu1,Bmunu2
  real(kind(1d0)) , dimension(0:3,0:3,0:3,0:3) :: Balbemunu
  real(kind(1d0)) :: CirjsCont
!
!
  interface 
    subroutine CalcB(p,smeB)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: smeB
!
    end subroutine CalcB
!
    subroutine CalcBmunu(ileg,p,Bmunu)
    use particles
    implicit none
!
      integer , intent(in) :: ileg
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(0:3,0:3) , intent(out) :: Bmunu
!
    end subroutine CalcBmunu
!
    subroutine CalcBalbemunu(ileg,jleg,p,Balbemunu)
    use particles
    implicit none
!
      integer , intent(in) :: ileg,jleg
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(0:3,0:3,0:3,0:3) , intent(out) :: Balbemunu
!
    end subroutine CalcBalbemunu
  end interface
!
  Cirjsterm = 0d0
! We go through all the subtraction terms and pick those which
! are contributing:
  do iterm=1,Cirjs%numterm
    if ((Cirjs%term(iterm)%emit(1)%i.eq.emiti).and. &
        (Cirjs%term(iterm)%emit(2)%i.eq.emitj).and. &
        (Cirjs%term(iterm)%rad(1)%i.eq.radi).and.  &
        (Cirjs%term(iterm)%rad(2)%i.eq.radr).and. &
        (Cirjs%term(iterm)%rad(3)%i.eq.radj).and. &
        (Cirjs%term(iterm)%rad(4)%i.eq.rads)) then
      call MapMomCirjsFF(p,ptilde, &
                         Cirjs%term(iterm)%emit(1)%i,  &
                         Cirjs%term(iterm)%emit(1)%ID, &
                         Cirjs%term(iterm)%rad(1)%i,   &
                         Cirjs%term(iterm)%rad(1)%ID,  &
                         Cirjs%term(iterm)%rad(2)%i,   &
                         Cirjs%term(iterm)%rad(2)%ID,  &
                         Cirjs%term(iterm)%emit(2)%i,  &
                         Cirjs%term(iterm)%emit(2)%ID, &
                         Cirjs%term(iterm)%rad(3)%i,   &
                         Cirjs%term(iterm)%rad(3)%ID,  &
                         Cirjs%term(iterm)%rad(4)%i,   &
                         Cirjs%term(iterm)%rad(4)%ID)
! Obtaining the underlying Born SME:
! If both emitters are quarks only the Born SME has to be
! calculated:
      if ((Cirjs%term(iterm)%emit(1)%ID.ne.0).and. &
          (Cirjs%term(iterm)%emit(2)%ID.ne.0)) then
!        print *,"q-q"
        call CalcB(ptilde,smeB)
! The first emitter is a quark, the second is a gluon:
      elseif ((Cirjs%term(iterm)%emit(1)%ID.ne.0).and. &
              (Cirjs%term(iterm)%emit(2)%ID.eq.0)) then
!        print *,"q-g"
        call CalcBmunu(Cirjs%term(iterm)%emit(2)%i,ptilde,Bmunu2)
        call CastSMEmunuToSME(ptilde,Bmunu2,smeB)
! The first emitter is a gluon, the second is a quark:
      elseif ((Cirjs%term(iterm)%emit(1)%ID.eq.0).and. &
              (Cirjs%term(iterm)%emit(2)%ID.ne.0)) then
!        print *,"g-q"
        call CalcBmunu(Cirjs%term(iterm)%emit(1)%i,ptilde,Bmunu1)
        call CastSMEmunuToSME(ptilde,Bmunu1,smeB)
! Both of the emitters are gluons:
      elseif ((Cirjs%term(iterm)%emit(1)%ID.eq.0).and. &
              (Cirjs%term(iterm)%emit(2)%ID.eq.0)) then
!        print *,"g-g"
        call CalcBalbemunu(Cirjs%term(iterm)%emit(1)%i, &
                           Cirjs%term(iterm)%emit(2)%i, &
                           ptilde,Balbemunu)
        call CastSMEalbemunuToSMEmunu(ptilde,Balbemunu,Bmunu1,Bmunu2)
        call CastSMEmunuToSME(ptilde,Bmunu1,smeB)
      end if
      call CalcCirjsFF(p,ptilde,                      &
                       smeB,Bmunu1,Bmunu2,Balbemunu,  &
                       Cirjs%term(iterm)%emit(1)%i,   &
                       Cirjs%term(iterm)%rad(1)%i,    &
                       Cirjs%term(iterm)%rad(2)%i,    &
                       Cirjs%term(iterm)%emit(2)%i,   &
                       Cirjs%term(iterm)%rad(3)%i,    &
                       Cirjs%term(iterm)%rad(4)%i,    &
                       CirjsCont)
      Cirjsterm = Cirjsterm + Cirjs%term(iterm)%symfact*CirjsCont
    end if
  end do
!
end subroutine PickCirjs
!
subroutine PickCSirs(p,phat,ptilde,Bij,Bmunuij,emit,radi,radr,rads, &
                     CSirs,CSirsterm)
use regions
use particles
use misc
use nnlo_subtractions
use nlo_mom_maps
use nnlo_mom_maps
implicit none
!
  type(particle), dimension(:), intent(in) :: p
  type(particle), dimension(:), intent(out) :: phat,ptilde
  real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
  real(kind(1d0)) , dimension(0:,0:,:,:) , intent(inout) :: Bmunuij
  integer , intent(in) :: emit,radi,radr,rads
  type(subterms) , intent(in) :: CSirs
  real(kind(1d0)) , intent(out) :: CSirsterm
!
  integer :: iterm,rad_s
  real(kind(1d0)) :: smeB
  real(kind(1d0)) :: CSirsCont
!
!
  interface 
    subroutine CalcBij(p,Bij)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Bij
!
    end subroutine CalcBij
!
    subroutine CalcBmunuij(ileg,p,Bmunuij)
    use particles
    implicit none
!
      integer , intent(in) :: ileg
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(0:,0:,:,:) , intent(out) :: Bmunuij
!
    end subroutine CalcBmunuij
  end interface
!
  CSirsterm = 0d0
! We go through all the subtraction terms and pick those which
! are contributing:
  do iterm=1,CSirs%numterm
    if ((CSirs%term(iterm)%emit(1)%i.eq.emit).and. &
        (CSirs%term(iterm)%rad(1)%i.eq.radi).and.  &
        (CSirs%term(iterm)%rad(2)%i.eq.radr).and. &
        (CSirs%term(iterm)%rad(3)%i.eq.rads)) then
! If the collinear emission is before the soft one in position
! the position of the soft parton has to be decreased by one:
      if (max(radi,radr).lt.rads) then
        rad_s = rads - 1
      else
        rad_s = rads
      end if
      call MapMomCirFF(p,phat,   &
                       CSirs%term(iterm)%emit(1)%i,  &
                       CSirs%term(iterm)%emit(1)%ID, &
                       CSirs%term(iterm)%rad(1)%i,   &
                       CSirs%term(iterm)%rad(1)%ID,  &
                       CSirs%term(iterm)%rad(2)%i,   &
                       CSirs%term(iterm)%rad(2)%ID)
      call MapMomSrNNLO(phat,ptilde, &
                        rads,rad_s, &
                        CSirs%term(iterm)%rad(3)%ID)
! Obtaining the underlying Born SME:
! If the emitter is a gluon we always need color- and 
! spin-correlated SME:
      if (CSirs%term(iterm)%emit(1)%ID.eq.0) then
        call CalcBmunuij(CSirs%term(iterm)%emit(1)%i,ptilde,Bmunuij)
        call CastSMEmunuijToSMEij(ptilde,Bmunuij,Bij)
! Otherwise just the color-correlated SME is requested:
      else
        call CalcBij(ptilde,Bij)
        call CastSMEijToSME(ptilde,Bij,smeB)
      end if
      call CalcCSirsFF(p,ptilde,smeB,Bij,Bmunuij,   &
                       CSirs%term(iterm)%emit(1)%i, &
                       CSirs%term(iterm)%rad(1)%i,  &
                       CSirs%term(iterm)%rad(2)%i,  &
                       CSirs%term(iterm)%rad(3)%i,  &
                       CSirsCont)
      CSirsterm = CSirsterm + CSirs%term(iterm)%symfact*CSirsCont
    end if
  end do
!
end subroutine PickCSirs
!
subroutine PickCirsCSirs(p,phat,ptilde,emit,radi,radr,rads, &
                         CSirs,CirsCSirsterm)
use regions
use particles
use misc
use nnlo_subtractions
use nlo_mom_maps
use nnlo_mom_maps
implicit none
!
  type(particle), dimension(:), intent(in) :: p
  type(particle), dimension(:), intent(out) :: phat,ptilde
  integer , intent(in) :: emit,radi,radr,rads
  type(subterms) , intent(in) :: CSirs
  real(kind(1d0)) , intent(out) :: CirsCSirsterm
!
  integer :: iterm,rad_s
  real(kind(1d0)) :: smeB
  real(kind(1d0)) , dimension(0:3,0:3) :: Bmunu
  real(kind(1d0)) :: CirsCSirsCont
!
!
  interface 
    subroutine CalcB(p,B)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: B
!
    end subroutine CalcB
!
    subroutine CalcBmunu(ileg,p,Bmunu)
    use particles
    implicit none
!
      integer , intent(in) :: ileg
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(0:3,0:3) , intent(out) :: Bmunu
!
    end subroutine CalcBmunu
  end interface
!
  CirsCSirsterm = 0d0
! We go through all the subtraction terms and pick those which
! are contributing:
  do iterm=1,CSirs%numterm
    if ((CSirs%term(iterm)%emit(1)%i.eq.emit).and. &
        (CSirs%term(iterm)%rad(1)%i.eq.radi).and.  &
        (CSirs%term(iterm)%rad(2)%i.eq.radr).and. &
        (CSirs%term(iterm)%rad(3)%i.eq.rads)) then
! If the collinear emission is before the soft one in position
! the position of the soft parton has to be decreased by one:
      if (max(radi,radr).lt.rads) then
        rad_s = rads - 1
      else
        rad_s = rads
      end if
      call MapMomCirFF(p,phat,   &
                       CSirs%term(iterm)%emit(1)%i,  &
                       CSirs%term(iterm)%emit(1)%ID, &
                       CSirs%term(iterm)%rad(1)%i,   &
                       CSirs%term(iterm)%rad(1)%ID,  &
                       CSirs%term(iterm)%rad(2)%i,   &
                       CSirs%term(iterm)%rad(2)%ID)
      call MapMomSrNNLO(phat,ptilde, &
                        rads,rad_s, &
                        CSirs%term(iterm)%rad(3)%ID)
! Obtaining the underlying Born SME:
! If the emitter is a gluon we always need spin-correlated SME:
      if (CSirs%term(iterm)%emit(1)%ID.eq.0) then
        call CalcBmunu(CSirs%term(iterm)%emit(1)%i,ptilde,Bmunu)
        call CastSMEmunuToSME(ptilde,Bmunu,smeB)
! Otherwise just the SME is requested:
      else
        call CalcB(ptilde,smeB)
      end if
      call CalcCirsCSirsFF(p,ptilde,smeB,Bmunu,         &
                           CSirs%term(iterm)%emit(1)%i, &
                           CSirs%term(iterm)%rad(1)%i,  &
                           CSirs%term(iterm)%rad(2)%i,  &
                           CSirs%term(iterm)%rad(3)%i,  &
                           CirsCSirsCont)
      CirsCSirsterm = CirsCSirsterm &
                    + CSirs%term(iterm)%symfact*CirsCSirsCont
    end if
  end do
!
end subroutine PickCirsCSirs
!
subroutine PickCirjsCSirs(p,phat,ptilde, &
                          emitir,radi,radr,emitjs,radj,rads, &
                          CSirs,CirjsCSirsterm)
use regions
use particles
use misc
use nnlo_subtractions
use nlo_mom_maps
use nnlo_mom_maps
implicit none
!
  type(particle), dimension(:), intent(in) :: p
  type(particle), dimension(:), intent(out) :: phat,ptilde
  integer , intent(in) :: emitir,emitjs,radi,radr,radj,rads
  type(subterms) , intent(in) :: CSirs
  real(kind(1d0)) , intent(out) :: CirjsCSirsterm
!
  integer :: iterm,rad_s
  real(kind(1d0)) :: smeB
  real(kind(1d0)) , dimension(0:3,0:3) :: Bmunu
  real(kind(1d0)) :: CirjsCSirsCont
!
!
  interface 
    subroutine CalcB(p,B)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: B
!
    end subroutine CalcB
!
    subroutine CalcBmunu(ileg,p,Bmunu)
    use particles
    implicit none
!
      integer , intent(in) :: ileg
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(0:3,0:3) , intent(out) :: Bmunu
!
    end subroutine CalcBmunu
  end interface
!
  CirjsCSirsterm = 0d0
! We go through all the subtraction terms and pick those which
! are contributing:
  do iterm=1,CSirs%numterm
    if ((CSirs%term(iterm)%emit(1)%i.eq.emitir).and. &
        (CSirs%term(iterm)%rad(1)%i.eq.radi).and.  &
        (CSirs%term(iterm)%rad(2)%i.eq.radr).and. &
        (CSirs%term(iterm)%rad(3)%i.eq.rads)) then
! If the collinear emission is before the soft one in position
! the position of the soft parton has to be decreased by one:
      if (max(radi,radr).lt.rads) then
        rad_s = rads - 1
      else
        rad_s = rads
      end if
      call MapMomCirFF(p,phat,   &
                       CSirs%term(iterm)%emit(1)%i,  &
                       CSirs%term(iterm)%emit(1)%ID, &
                       CSirs%term(iterm)%rad(1)%i,   &
                       CSirs%term(iterm)%rad(1)%ID,  &
                       CSirs%term(iterm)%rad(2)%i,   &
                       CSirs%term(iterm)%rad(2)%ID)
      call MapMomSrNNLO(phat,ptilde, &
                        rads,rad_s, &
                        CSirs%term(iterm)%rad(3)%ID)
! Obtaining the underlying Born SME:
! If the emitter is a gluon we always need spin-correlated SME:
      if (CSirs%term(iterm)%emit(1)%ID.eq.0) then
        call CalcBmunu(CSirs%term(iterm)%emit(1)%i,ptilde,Bmunu)
        call CastSMEmunuToSME(ptilde,Bmunu,smeB)
! Otherwise just the SME is requested:
      else
        call CalcB(ptilde,smeB)
      end if
      call CalcCirjsCSirsFF(p,ptilde,smeB,Bmunu, &
                            emitir,radi,radr,    &
                            emitjs,radj,rads,    &
                            CirjsCSirsCont)
      CirjsCSirsterm = CirjsCSirsterm &
                    + CSirs%term(iterm)%symfact*CirjsCSirsCont
    end if
  end do
!
end subroutine PickCirjsCSirs
!
subroutine PickSrs(p,ptilde,Bij,Bijkl,radr,rads,Srs,Srsterm)
use regions
use particles
use misc
use nnlo_subtractions
use nnlo_mom_maps
implicit none
!
  type(particle), dimension(:), intent(in) :: p
  type(particle), dimension(:), intent(out) :: ptilde
  real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
  real(kind(1d0)) , dimension(:,:,:,:) , intent(inout) :: Bijkl
  integer , intent(in) :: radr,rads
  type(subterms) , intent(in) :: Srs
  real(kind(1d0)) , intent(out) :: Srsterm
!
  integer :: iterm
  real(kind(1d0)) :: SrsCont
!
!
  interface 
    subroutine CalcBij(p,Bij)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Bij
!
    end subroutine CalcBij
!
    subroutine CalcBijkl(p,Bijkl)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:,:,:,:) , intent(out) :: Bijkl
!
    end subroutine CalcBijkl
  end interface
!
  Srsterm = 0d0
! We go through all the subtraction terms and pick those which
! are contributing:
  do iterm=1,Srs%numterm
    if (((Srs%term(iterm)%rad(1)%i.eq.radr).and.  &
         (Srs%term(iterm)%rad(2)%i.eq.rads)).or.  &
        ((Srs%term(iterm)%rad(1)%i.eq.rads).and.  &
         (Srs%term(iterm)%rad(2)%i.eq.radr))) then
      call MapMomSrs(p,ptilde, &
                     Srs%term(iterm)%rad(1)%i,   &
                     Srs%term(iterm)%rad(2)%i)
!      call PrintParts(p)
!      call PrintParts(ptilde)
!      call CheckParts(ptilde)
      call CalcBijkl(ptilde,Bijkl)
      call CastSMEijklToSMEij(.false.,ptilde,Bijkl,Bij)
      call CalcSrs(p,ptilde,Bij,Bijkl,       &
                   Srs%term(iterm)%rad(1)%i, &
                   Srs%term(iterm)%rad(2)%i, &
                   SrsCont)
      Srsterm = Srsterm + Srs%term(iterm)%symfact*SrsCont
    end if
  end do
!
end subroutine PickSrs
!
subroutine PickCirsSrs(p,ptilde,emitirs,radi,radr,rads,Srs,CirsSrsterm)
use regions
use particles
use misc
use nnlo_subtractions
use nnlo_mom_maps
implicit none
!
  type(particle), dimension(:), intent(in) :: p
  type(particle), dimension(:), intent(out) :: ptilde
  integer , intent(in) :: emitirs,radi,radr,rads
  type(subterms) , intent(in) :: Srs
  real(kind(1d0)) , intent(out) :: CirsSrsterm
!
  integer :: iterm
  real(kind(1d0)) :: smeB,CirsSrsCont
!
!
  interface 
    subroutine CalcB(p,smeB)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: smeB
!
    end subroutine CalcB
  end interface
!
  CirsSrsterm = 0d0
! We go through all the subtraction terms and pick those which
! are contributing:
  do iterm=1,Srs%numterm
    if (((Srs%term(iterm)%rad(1)%i.eq.radr).and.  &
         (Srs%term(iterm)%rad(2)%i.eq.rads)).or. &
        ((Srs%term(iterm)%rad(1)%i.eq.rads).and.  &
         (Srs%term(iterm)%rad(2)%i.eq.radr))) then
      call MapMomSrs(p,ptilde, &
                     Srs%term(iterm)%rad(1)%i,   &
                     Srs%term(iterm)%rad(2)%i)
!      call PrintParts(p)
!      call PrintParts(ptilde)
!      call CheckParts(ptilde)
      call CalcB(ptilde,smeB)
      call CalcCirsSrs(p,ptilde,smeB,          &
                       emitirs,radi,radr,rads, &
                       CirsSrsCont)
      CirsSrsterm = CirsSrsterm &
                  + Srs%term(iterm)%symfact*CirsSrsCont
    end if
  end do
!
end subroutine PickCirsSrs
!
subroutine PickCSirsSrs(p,ptilde,Bij,emitir,radi,radr,rads,Srs, &
                        CSirsSrsterm)
use regions
use particles
use misc
use nnlo_subtractions
use nnlo_mom_maps
implicit none
!
  type(particle), dimension(:), intent(in) :: p
  type(particle), dimension(:), intent(out) :: ptilde
  real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
  integer , intent(in) :: emitir,radi,radr,rads
  type(subterms) , intent(in) :: Srs
  real(kind(1d0)) , intent(out) :: CSirsSrsterm
!
  integer :: iterm
  real(kind(1d0)) :: CSirsSrsCont
!
!
  interface 
    subroutine CalcBij(p,Bij)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Bij
!
    end subroutine CalcBij
  end interface
!
  CSirsSrsterm = 0d0
! We go through all the subtraction terms and pick those which
! are contributing:
  do iterm=1,Srs%numterm
    if (((Srs%term(iterm)%rad(1)%i.eq.radr).and.  &
         (Srs%term(iterm)%rad(2)%i.eq.rads)).or.  &
        ((Srs%term(iterm)%rad(1)%i.eq.rads).and.  &
         (Srs%term(iterm)%rad(2)%i.eq.radr))) then
      call MapMomSrs(p,ptilde,radr,rads)
!      call PrintParts(p)
!      call PrintParts(ptilde)
!      call CheckParts(ptilde)
      call CalcBij(ptilde,Bij)
      call CalcCSirsSrs(p,ptilde,Bij,          &
                        emitir,radi,radr,rads, &
                        CSirsSrsCont)
      CSirsSrsterm = CSirsSrsterm &
                   + Srs%term(iterm)%symfact*CSirsSrsCont
    end if
  end do
!
end subroutine PickCSirsSrs
!
subroutine PickCirjsSrs(p,ptilde, &
                        emitir,radi,radr, &
                        emitjs,radj,rads, &
                        Srs,CirjsSrsterm)
use regions
use particles
use misc
use nnlo_subtractions
use nnlo_mom_maps
implicit none
!
  type(particle), dimension(:), intent(in) :: p
  type(particle), dimension(:), intent(out) :: ptilde
  integer , intent(in) :: emitir,emitjs,radi,radr,radj,rads
  type(subterms) , intent(in) :: Srs
  real(kind(1d0)) , intent(out) :: CirjsSrsterm
!
  integer :: iterm
  real(kind(1d0)) :: smeB,CirjsSrsCont
!
!
  interface 
    subroutine CalcB(p,smeB)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: smeB
!
    end subroutine CalcB
  end interface
!
  CirjsSrsterm = 0d0
! We go through all the subtraction terms and pick those which
! are contributing:
  do iterm=1,Srs%numterm
    if ((Srs%term(iterm)%rad(1)%i.eq.radr).and.  &
        (Srs%term(iterm)%rad(2)%i.eq.rads)) then
      call MapMomSrs(p,ptilde, &
                     Srs%term(iterm)%rad(1)%i,   &
                     Srs%term(iterm)%rad(2)%i)
!      call PrintParts(p)
!      call PrintParts(ptilde)
!      call CheckParts(ptilde)
      call CalcB(ptilde,smeB)
      call CalcCirjsSrs(p,ptilde,smeB,    &
                        emitir,radi,radr, &
                        emitjs,radj,rads, &
                        CirjsSrsCont)
      CirjsSrsterm = CirjsSrsterm &
                   + Srs%term(iterm)%symfact*CirjsSrsCont
    end if
  end do
!
end subroutine PickCirjsSrs
!
subroutine PickCirsCSirsSrs(p,ptilde,emitirs,radi,radr,rads,Srs, &
                            CirsCSirsSrsterm)
use regions
use particles
use misc
use nnlo_subtractions
use nnlo_mom_maps
implicit none
!
  type(particle), dimension(:), intent(in) :: p
  type(particle), dimension(:), intent(out) :: ptilde
  integer , intent(in) :: emitirs,radi,radr,rads
  type(subterms) , intent(in) :: Srs
  real(kind(1d0)) , intent(out) :: CirsCSirsSrsterm
!
  integer :: iterm
  real(kind(1d0)) :: smeB,CirsCSirsSrsCont
!
!
  interface 
    subroutine CalcB(p,smeB)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: smeB
!
    end subroutine CalcB
  end interface
!
  CirsCSirsSrsterm = 0d0
! We go through all the subtraction terms and pick those which
! are contributing:
  do iterm=1,Srs%numterm
    if (((Srs%term(iterm)%rad(1)%i.eq.radr).and.  &
         (Srs%term(iterm)%rad(2)%i.eq.rads)).or.  &
        ((Srs%term(iterm)%rad(1)%i.eq.rads).and.  &
         (Srs%term(iterm)%rad(2)%i.eq.radr))) then
      call MapMomSrs(p,ptilde, &
                     Srs%term(iterm)%rad(1)%i,   &
                     Srs%term(iterm)%rad(2)%i)
!      call PrintParts(p)
!      call PrintParts(ptilde)
!      call CheckParts(ptilde)
      call CalcB(ptilde,smeB)
      call CalcCirsCSirsSrs(p,ptilde,smeB,          &
                            emitirs,radi,radr,rads, &
                            CirsCSirsSrsCont)
      CirsCSirsSrsterm = CirsCSirsSrsterm &
                       + Srs%term(iterm)%symfact*CirsCSirsSrsCont
    end if
  end do
!
end subroutine PickCirsCSirsSrs
!
subroutine PickCktCktr(p,phat,ptilde, &
                       emitktr,radk,radt,radr,rad1,rad2,rad3, &
                       Cktr,CktCktrterm)
use regions
use particles
use misc
use nnlo_subtractions
use nnlo_mom_maps
implicit none
!
  type(particle), dimension(:), intent(in) :: p
  type(particle), dimension(:), intent(out) :: phat,ptilde
  integer , intent(in) :: emitktr,radk,radt,radr,rad1,rad2,rad3
  type(subterms) , intent(in) :: Cktr
  real(kind(1d0)) , intent(out) :: CktCktrterm
!
  integer :: iterm
  integer :: rad_kt,rad_r
  real(kind(1d0)) :: smeB
  real(kind(1d0)) , dimension(0:3,0:3) :: Bmunu
  real(kind(1d0)) :: CktCktrCont
!
!
  interface 
    subroutine CalcB(p,smeB)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: smeB
!
    end subroutine CalcB
!
    subroutine CalcBmunu(ileg,p,Bmunu)
    use particles
    implicit none
!
      integer , intent(in) :: ileg
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(0:3,0:3) , intent(out) :: Bmunu
!
    end subroutine CalcBmunu
  end interface
!
  CktCktrterm = 0d0
! Due to the iterated mappings the first mapping can change the
! position of the radiated partons, hence if needed the position
! is shifted:
  rad_kt = min(radk,radt)
  rad_r  = radr
  if (radr.gt.max(radk,radt)) rad_r = radr - 1
! Going through the mother subtraction terms to identify the
! one which corresponds to Cktr in order to get the correct
! prefactor, to find the correct Cktr term the radiated legs
! have to be in a certain order this is the reason why two sets
! of rad legs are supplied to the routine:
  do iterm=1,Cktr%numterm
    if ((Cktr%term(iterm)%rad(1)%i.eq.rad1).and.  &
        (Cktr%term(iterm)%rad(2)%i.eq.rad2).and. &
        (Cktr%term(iterm)%rad(3)%i.eq.rad3)) then
      call MapMomCirFF_A12(p,phat,                &
                           radk,radt,             &
                           Q,Q2,yiQ_arr,sir_arr,  &
                           s_k_t,ztildek,ztildet, &
                           kttildek,kttildet,     &
                           alphakt)
      call CalcSubhatInvariants(phat)
      call MapMomCirFF_A12(phat,ptilde,                &
                           rad_kt,rad_r,               &
                           Q,Q2,yihatQ_arr,sirhat_arr, &
                           s_kt_r,zhatkt,zhatr,        &
                           kthatkt,kthatr,             &
                           alphaktr)
! Obtaining the underlying Born SME:
! If the emitter is a gluon we always need spin-correlated SME:
      if (Cktr%term(iterm)%emit(1)%ID.eq.0) then
        call CalcBmunu(Cktr%term(iterm)%emit(1)%i,ptilde,Bmunu)
        call CastSMEmunuToSME(ptilde,Bmunu,smeB)
! Otherwise just the SME is requested:
      else
        call CalcB(ptilde,smeB)
      end if
      call CalcCktCktrFF(p,phat,ptilde,smeB,Bmunu,            &
                         emitktr,radk,radt,radr,rad_kt,rad_r, &
                         CktCktrCont)
      CktCktrterm = CktCktrterm + Cktr%term(iterm)%symfact*CktCktrCont
    end if
  end do
!
end subroutine PickCktCktr
!
subroutine PickCktCirkt(p,phat,ptilde, &
                        emitir,radi,radr,emitkt,radk,radt, &
                        rad1,rad2,rad3,rad4, &
                        Cirkt,CktCirktterm)
use regions
use particles
use misc
use nnlo_subtractions
use nnlo_mom_maps
implicit none
!
  type(particle), dimension(:), intent(in) :: p
  type(particle), dimension(:), intent(out) :: phat,ptilde
  integer , intent(in) :: emitir,emitkt,radi,radr,radk,radt,rad1,rad2,rad3,rad4
  type(subterms) , intent(in) :: Cirkt
  real(kind(1d0)) , intent(out) :: CktCirktterm
!
  integer :: iterm
  integer :: rad_i,rad_r
  real(kind(1d0)) :: smeB
  real(kind(1d0)) , dimension(0:3,0:3) :: Bmunu1,Bmunu2
  real(kind(1d0)) , dimension(0:3,0:3,0:3,0:3) :: Balbemunu
  real(kind(1d0)) :: CktCirktCont
!
!
  interface 
    subroutine CalcB(p,smeB)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: smeB
!
    end subroutine CalcB
!
    subroutine CalcBmunu(ileg,p,Bmunu)
    use particles
    implicit none
!
      integer , intent(in) :: ileg
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(0:3,0:3) , intent(out) :: Bmunu
!
    end subroutine CalcBmunu
!
    subroutine CalcBalbemunu(ileg,jleg,p,Balbemunu)
    use particles
    implicit none
!
      integer , intent(in) :: ileg,jleg
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(0:3,0:3,0:3,0:3) , intent(out) :: Balbemunu
!
    end subroutine CalcBalbemunu
  end interface
!
  CktCirktterm = 0d0
! Due to the iterated mappings the first mapping can change the
! position of the radiated partons, hence if needed the position
! is shifted:
  rad_i = radi
  rad_r = radr
  if (radi.gt.max(radk,radt)) rad_i = radi - 1
  if (radr.gt.max(radk,radt)) rad_r = radr - 1
! Going through the mother subtraction terms to identify the
! one which corresponds to Cirkt in order to get the correct
! prefactor, to find the correct Cirkt term the radiated legs
! have to be in a certain order this is the reason why two sets
! of rad legs are supplied to the routine:
  do iterm=1,Cirkt%numterm
    if ((Cirkt%term(iterm)%rad(1)%i.eq.rad1).and.  &
        (Cirkt%term(iterm)%rad(2)%i.eq.rad2).and. &
        (Cirkt%term(iterm)%rad(3)%i.eq.rad3).and. &
        (Cirkt%term(iterm)%rad(4)%i.eq.rad4)) then
      call MapMomCirFF_A12(p,phat,                &
                           radk,radt,             &
                           Q,Q2,yiQ_arr,sir_arr,  &
                           s_k_t,ztildek,ztildet, &
                           kttildek,kttildet,     &
                           alphakt)
      call CalcSubhatInvariants(phat)
      call MapMomCirFF_A12(phat,ptilde,                &
                           rad_i,rad_r,                &  
                           Q,Q2,yihatQ_arr,sirhat_arr, &
                           sir,zhati,zhatr,            &
                           kthati,kthatr,              &
                           alphair)
! Obtaining the underlying Born SME:
! If both emitters are quarks only the Born SME has to be
! calculated:
      if ((ptilde(emitir)%flv.ne.0).and. &
          (ptilde(emitkt)%flv.ne.0)) then
!        print *,"q-q"
        call CalcB(ptilde,smeB)
! The first emitter is a quark, the second is a gluon:
      elseif ((ptilde(emitir)%flv.ne.0).and. &
              (ptilde(emitkt)%flv.eq.0)) then
!        print *,"q-g"
        call CalcBmunu(emitkt,ptilde,Bmunu2)
        call CastSMEmunuToSME(ptilde,Bmunu2,smeB)
! The first emitter is a gluon, the second is a quark:
      elseif ((ptilde(emitir)%flv.eq.0).and. &
              (ptilde(emitkt)%flv.ne.0)) then
!        print *,"g-q"
        call CalcBmunu(emitir,ptilde,Bmunu1)
        call CastSMEmunuToSME(ptilde,Bmunu1,smeB)
! Both of the emitters are gluons:
      elseif ((ptilde(emitir)%flv.eq.0).and. &
              (ptilde(emitkt)%flv.eq.0)) then
!        print *,"g-g"
        call CalcBalbemunu(emitir,emitkt, &
                           ptilde,Balbemunu)
        call CastSMEalbemunuToSMEmunu(ptilde,Balbemunu,Bmunu1,Bmunu2)
        call CastSMEmunuToSME(ptilde,Bmunu1,smeB)
      end if
      call CalcCktCirktFF(p,phat,ptilde,smeB,Bmunu1,Bmunu2,Balbemunu, &
                          radi,radr,radk,radt, &
                          CktCirktCont)
      CktCirktterm = CktCirktterm + Cirkt%term(iterm)%symfact*CktCirktCont
    end if
  end do
!
end subroutine PickCktCirkt
!
subroutine PickCktCSktr(p,phat,ptilde,Bij,Bmunuij, &
                        emitkt,radk,radt,radr, &
                        CSktr,CktCSktrterm)
use regions
use particles
use misc
use nnlo_subtractions
use nnlo_mom_maps
implicit none
!
  type(particle), dimension(:), intent(in) :: p
  type(particle), dimension(:), intent(out) :: phat,ptilde
  real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
  real(kind(1d0)) , dimension(0:,0:,:,:) , intent(inout) :: Bmunuij
  integer , intent(in) :: emitkt,radk,radt,radr
  type(subterms) , intent(in) :: CSktr
  real(kind(1d0)) , intent(out) :: CktCSktrterm
!
  integer :: iterm
  integer :: rad_kt,rad_r
  real(kind(1d0)) :: yrhatQ
  real(kind(1d0)) :: CktCSktrCont
!
!
  interface 
    subroutine CalcBij(p,Bij)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Bij
!
    end subroutine CalcBij
!
    subroutine CalcBmunuij(ileg,p,Bmunuij)
    use particles
    implicit none
!
      integer , intent(in) :: ileg
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(0:,0:,:,:) , intent(out) :: Bmunuij
!
    end subroutine CalcBmunuij
  end interface
!
  CktCSktrterm = 0d0
! Due to the iterated mappings the first mapping can change the
! position of the radiated partons, hence if needed the position
! is shifted:
  rad_kt = min(radk,radt)
  rad_r  = radr
  if (radr.gt.max(radk,radt)) rad_r = radr - 1
!  print *,"rad_kt,rad_r: ",rad_kt,rad_r
!  print *,"radk,radt: ",radk,radt
!  print *,"rad1,rad2,rad3: ",rad1,rad2,rad3
! Going through the mother subtraction terms to identify the
! one which corresponds to CSktr in order to get the correct
! prefactor, to find the correct Cktr term the radiated legs
! have to be in a certain order this is the reason why two sets
! of rad legs are supplied to the routine:
  do iterm=1,CSktr%numterm
    if ((CSktr%term(iterm)%rad(1)%i.eq.radk).and.  &
        (CSktr%term(iterm)%rad(2)%i.eq.radt).and. &
        (CSktr%term(iterm)%rad(3)%i.eq.radr)) then
      call MapMomCirFF_A12(p,phat,                &
                           radk,radt,             &
                           Q,Q2,yiQ_arr,sir_arr,  &
                           s_k_t,ztildek,ztildet, &
                           kttildek,kttildet,     &
                           alphakt)
      call CalcSubhatInvariants(phat)
      call MapMomSr_A12(phat,ptilde,     &
                        rad_r,           &
                        Q,Q2,yihatQ_arr, &
                        yrhatQ,lambda_r)
! The soft mapping involves a boost this has to be applied
! even on the transverse momenta:
      call LorentzLambda(Q,(Q - phat(rad_r)%p)/lambda_r, &
                         kttildek/lambda_r,kttildek)
      call LorentzLambda(Q,(Q - phat(rad_r)%p)/lambda_r, &
                         kttildet/lambda_r,kttildet)
! Obtaining the underlying Born SME:
! If the emitter is a gluon we always need spin- and 
! color-correlated SME:
      if (CSktr%term(iterm)%emit(1)%ID.eq.0) then
        call CalcBmunuij(CSktr%term(iterm)%emit(1)%i,ptilde,Bmunuij)
! Otherwise just the color-correlated is requested:
      else
        call CalcBij(ptilde,Bij)
      end if
      call CalcCktCSktrFF(p,phat,ptilde,Bij,Bmunuij,           &
                          radk,radt,rad_r, &
                          CktCSktrCont)
      CktCSktrterm = CktCSktrterm + CSktr%term(iterm)%symfact*CktCSktrCont
    end if
  end do
!
end subroutine PickCktCSktr
!
subroutine PickCktCirktCSktr(p,phat,ptilde, &
                             emitir,radi,radr,emitkt,radk,radt, &
                             CSktr,CktCirktCSktrterm)
use regions
use particles
use misc
use nnlo_subtractions
use nnlo_mom_maps
implicit none
!
  type(particle), dimension(:), intent(in) :: p
  type(particle), dimension(:), intent(out) :: phat,ptilde
  integer , intent(in) :: emitir,emitkt,radi,radr,radk,radt
  type(subterms) , intent(in) :: CSktr
  real(kind(1d0)) , intent(out) :: CktCirktCSktrterm
!
  integer :: iterm
  integer :: rad_i,rad_r
  real(kind(1d0)) :: smeB
  real(kind(1d0)) :: yrhatQ
  real(kind(1d0)) , dimension(0:3,0:3) :: Bmunu
  real(kind(1d0)) :: CktCirktCSktrCont
!
!
  interface 
    subroutine CalcB(p,smeB)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: smeB
!
    end subroutine CalcB
!
    subroutine CalcBmunu(ileg,p,Bmunu)
    use particles
    implicit none
!
      integer , intent(in) :: ileg
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(0:3,0:3) , intent(out) :: Bmunu
!
    end subroutine CalcBmunu
  end interface
!
  CktCirktCSktrterm = 0d0
! Due to the iterated mappings the first mapping can change the
! position of the radiated partons, hence if needed the position
! is shifted:
  rad_i = radi
  rad_r = radr
  if (radi.gt.max(radk,radt)) rad_i = radi - 1
  if (radr.gt.max(radk,radt)) rad_r = radr - 1
! To identify the correct mother counterterm there is no
! need to alter the ordering of the legs provided we started
! with the CSktr counterterm and added an extra leg only.
  do iterm=1,CSktr%numterm
    if ((CSktr%term(iterm)%rad(1)%i.eq.radk).and.  &
        (CSktr%term(iterm)%rad(2)%i.eq.radt).and. &
        (CSktr%term(iterm)%rad(3)%i.eq.radr)) then
      call MapMomCirFF_A12(p,phat,                &
                           radk,radt,             &
                           Q,Q2,yiQ_arr,sir_arr,  &
                           s_k_t,ztildek,ztildet, &
                           kttildek,kttildet,     &
                           alphakt)
      call CalcSubhatInvariants(phat)
      call MapMomSr_A12(phat,ptilde,     &
                        rad_r,           &
                        Q,Q2,yihatQ_arr, &
                        yrhatQ,lambda_r)
! The soft mapping involves a boost this has to be applied
! even on the transverse momenta:
      call LorentzLambda(Q,(Q - phat(rad_r)%p)/lambda_r, &
                         kttildek/lambda_r,kttildek)
      call LorentzLambda(Q,(Q - phat(rad_r)%p)/lambda_r, &
                         kttildet/lambda_r,kttildet)
! Obtaining the underlying Born SME:
! If the emitter is a gluon we always need spin-correlated SME:
      if (CSktr%term(iterm)%emit(1)%ID.eq.0) then
        call CalcBmunu(CSktr%term(iterm)%emit(1)%i,ptilde,Bmunu)
! Otherwise just the color-correlated is requested:
      else
        call CalcB(ptilde,smeB)
      end if
      call CalcCktCirktCSktrFF(p,phat,ptilde,smeB,Bmunu, &
                               rad_i,rad_r,radk,radt, &
                               CktCirktCSktrCont)
      CktCirktCSktrterm = CktCirktCSktrterm &
                        + CSktr%term(iterm)%symfact*CktCirktCSktrCont
    end if
  end do
!
end subroutine PickCktCirktCSktr
!
subroutine PickCktCktrCSktr(p,phat,ptilde, &
                            emitktr,radk,radt,radr, &
                            CSktr,CktCktrCSktrterm)
use regions
use particles
use misc
use nnlo_subtractions
use nnlo_mom_maps
implicit none
!
  type(particle), dimension(:), intent(in) :: p
  type(particle), dimension(:), intent(out) :: phat,ptilde
  integer , intent(in) :: emitktr,radk,radt,radr
  type(subterms) , intent(in) :: CSktr
  real(kind(1d0)) , intent(out) :: CktCktrCSktrterm
!
  integer :: iterm
  integer :: rad_kt,rad_r
  real(kind(1d0)) :: smeB
  real(kind(1d0)) :: yrhatQ
  real(kind(1d0)) , dimension(0:3,0:3) :: Bmunu
  real(kind(1d0)) :: CktCktrCSktrCont
!
!
  interface 
    subroutine CalcB(p,smeB)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: smeB
!
    end subroutine CalcB
!
    subroutine CalcBmunu(ileg,p,Bmunu)
    use particles
    implicit none
!
      integer , intent(in) :: ileg
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(0:3,0:3) , intent(out) :: Bmunu
!
    end subroutine CalcBmunu
  end interface
!
  CktCktrCSktrterm = 0d0
! Due to the iterated mappings the first mapping can change the
! position of the radiated partons, hence if needed the position
! is shifted:
  rad_kt = min(radk,radt)
  rad_r  = radr
  if (radr.gt.max(radk,radt)) rad_r = radr - 1
! To identify the correct mother counterterm there is no
! need to alter the ordering of the legs provided we started
! with the CSktr counterterm and added an extra leg only.
  do iterm=1,CSktr%numterm
    if ((CSktr%term(iterm)%rad(1)%i.eq.radk).and.  &
        (CSktr%term(iterm)%rad(2)%i.eq.radt).and. &
        (CSktr%term(iterm)%rad(3)%i.eq.radr)) then
      call MapMomCirFF_A12(p,phat,                &
                           radk,radt,             &
                           Q,Q2,yiQ_arr,sir_arr,  &
                           s_k_t,ztildek,ztildet, &
                           kttildek,kttildet,     &
                           alphakt)
      call CalcSubhatInvariants(phat)
      call MapMomSr_A12(phat,ptilde,     &
                        rad_r,           &
                        Q,Q2,yihatQ_arr, &
                        yrhatQ,lambda_r)
! The soft mapping involves a boost this has to be applied
! even on the transverse momenta:
      call LorentzLambda(Q,(Q - phat(rad_r)%p)/lambda_r, &
                         kttildek/lambda_r,kttildek)
      call LorentzLambda(Q,(Q - phat(rad_r)%p)/lambda_r, &
                         kttildet/lambda_r,kttildet)
! Obtaining the underlying Born SME:
! If the emitter is a gluon we always need spin-correlated SME:
      if (CSktr%term(iterm)%emit(1)%ID.eq.0) then
        call CalcBmunu(CSktr%term(iterm)%emit(1)%i,ptilde,Bmunu)
! Otherwise just the color-correlated is requested:
      else
        call CalcB(ptilde,smeB)
      end if
      call CalcCktCktrCSktrFF(p,phat,ptilde,smeB,Bmunu, &
                              rad_kt,rad_r,radk,radt, &
                              CktCktrCSktrCont)
      CktCktrCSktrterm = CktCktrCSktrterm &
                       + CSktr%term(iterm)%symfact*CktCktrCSktrCont
    end if
  end do
!
end subroutine PickCktCktrCSktr
!
subroutine PickCktSkt(p,phat,ptilde,Bij, &
                      emitkt,radk,radt,  &
                      Skt,CktSktterm)
use regions
use particles
use misc
use nnlo_subtractions
use nnlo_mom_maps
implicit none
!
  type(particle), dimension(:), intent(in) :: p
  type(particle), dimension(:), intent(out) :: phat,ptilde
  real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
  integer , intent(in) :: emitkt,radk,radt
  type(subterms) , intent(in) :: Skt
  real(kind(1d0)) , intent(out) :: CktSktterm
!
  integer :: iterm
  integer :: rad_kt
  real(kind(1d0)) :: yrhatQ
  real(kind(1d0)) :: CktSktCont
!
!
  interface 
    subroutine CalcBij(p,Bij)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Bij
!
    end subroutine CalcBij
  end interface
!
  CktSktterm = 0d0
  rad_kt = min(radk,radt)
  do iterm=1,Skt%numterm
    if ((Skt%term(iterm)%rad(1)%i.eq.radk).and.  &
        (Skt%term(iterm)%rad(2)%i.eq.radt)) then
      call MapMomCirFF_A12(p,phat,                &
                           radk,radt,             &
                           Q,Q2,yiQ_arr,sir_arr,  &
                           s_k_t,ztildek,ztildet, &
                           kttildek,kttildet,     &
                           alphakt)
      call CalcSubhatInvariants(phat)
      call MapMomSr_A12(phat,ptilde,     &
                        rad_kt,          &
                        Q,Q2,yihatQ_arr, &
                        yrhatQ,lambda_r)
! Obtaining the underlying Born SME:
      call CalcBij(ptilde,Bij)
      call CalcCktSktFF(p,phat,ptilde,Bij, &
                        rad_kt,radk,radt,  &
                        CktSktCont)
      CktSktterm = CktSktterm &
                 + Skt%term(iterm)%symfact*CktSktCont
    end if
  end do
!
end subroutine PickCktSkt
!
subroutine PickCktCrktSkt(p,phat,ptilde,      &
                          radr,radk,radt,     &
                          Skt,CktCrktSktterm)
use regions
use particles
use misc
use nnlo_subtractions
use nnlo_mom_maps
implicit none
!
  type(particle), dimension(:), intent(in) :: p
  type(particle), dimension(:), intent(out) :: phat,ptilde
  integer , intent(in) :: radr,radk,radt
  type(subterms) , intent(in) :: Skt
  real(kind(1d0)) , intent(out) :: CktCrktSktterm
!
  integer :: iterm
  integer :: rad_kt,rad_r
  real(kind(1d0)) :: smeB
  real(kind(1d0)) :: yrhatQ
  real(kind(1d0)) :: CktCrktSktCont
!
!
  interface
    subroutine CalcB(p,smeB)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: smeB
!
    end subroutine CalcB
  end interface
!
  CktCrktSktterm = 0d0
! Due to the iterated mappings the first mapping can change the
! position of the radiated partons, hence if needed the position
! is shifted:
  rad_kt = min(radk,radt)
  rad_r  = radr
  if (rad_r.gt.max(radk,radt)) rad_r = rad_r - 1
! To identify the correct mother counterterm there is no
! need to alter the ordering of the legs provided we started
! with the Skt counterterm and added an extra leg only.
  do iterm=1,Skt%numterm
    if ((Skt%term(iterm)%rad(1)%i.eq.radk).and.  &
        (Skt%term(iterm)%rad(2)%i.eq.radt)) then
      call MapMomCirFF_A12(p,phat,                &
                           radk,radt,             &
                           Q,Q2,yiQ_arr,sir_arr,  &
                           s_k_t,ztildek,ztildet, &
                           kttildek,kttildet,     &
                           alphakt)
      call CalcSubhatInvariants(phat)
      call MapMomSr_A12(phat,ptilde,     &
                        rad_kt,          &
                        Q,Q2,yihatQ_arr, &
                        yrhatQ,lambda_r)
! Obtaining the underlying Born SME:
      call CalcB(ptilde,smeB)
      call CalcCktCrktSktFF(p,phat,ptilde,smeB,           &
                            rad_r,rad_kt,radr,radk,radt,  &
                            CktCrktSktCont)
      CktCrktSktterm = CktCrktSktterm &
                 + Skt%term(iterm)%symfact*CktCrktSktCont
    end if
  end do
!
end subroutine PickCktCrktSkt
!
subroutine PickStCirt(p,phat,ptilde, &
                      emitirt,radi,radr,radt,rad1,rad2,rad3, &
                      Cirs,StCirtterm)
use regions
use particles
use misc
use nnlo_subtractions
use nnlo_mom_maps
implicit none
!
  type(particle), dimension(:), intent(in) :: p
  type(particle), dimension(:), intent(out) :: phat,ptilde
  integer , intent(in) :: emitirt,radi,radr,radt,rad1,rad2,rad3
  type(subterms) , intent(in) :: Cirs
  real(kind(1d0)) , intent(out) :: StCirtterm
!
  integer :: iterm
  integer :: rad_i,rad_r,rad_ir
  real(kind(1d0)) :: smeB
  real(kind(1d0)) :: ytQ
  real(kind(1d0)) , dimension(0:3,0:3) :: Bmunu
  real(kind(1d0)) :: StCirtCont
!
!
  interface 
    subroutine CalcB(p,smeB)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: smeB
!
    end subroutine CalcB
!
    subroutine CalcBmunu(ileg,p,Bmunu)
    use particles
    implicit none
!
      integer , intent(in) :: ileg
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(0:3,0:3) , intent(out) :: Bmunu
!
    end subroutine CalcBmunu
  end interface
!
  StCirtterm = 0d0
! Position of parton i and r among hat momenta: 
  rad_i = radi
  rad_r = radr
  if (rad_i.gt.radt) rad_i = rad_i - 1
  if (rad_r.gt.radt) rad_r = rad_r - 1
! Position for the \hat{ir} parton 
  rad_ir = min(rad_i,rad_r)
  do iterm=1,Cirs%numterm
    if ((Cirs%term(iterm)%rad(1)%i.eq.rad1).and.  &
        (Cirs%term(iterm)%rad(2)%i.eq.rad2).and. &
        (Cirs%term(iterm)%rad(3)%i.eq.rad3)) then
      call MapMomSr_A12(p,phat,       &
                        radt,         &
                        Q,Q2,yiQ_arr, &
                        ytQ,lambda_t)
      call CalcSubhatInvariants(phat)
      call MapMomCirFF_A12(phat,ptilde,                &
                           rad_i,rad_r,                &
                           Q,Q2,yihatQ_arr,sirhat_arr, &
                           sir,zhati,zhatr,            &
                           kthati,kthatr,              &
                           alphairt)
! Obtaining the underlying Born SME:
! If the emitter is a gluon we always need spin-correlated SME:
      if (Cirs%term(iterm)%emit(1)%ID.eq.0) then
        call CalcBmunu(Cirs%term(iterm)%emit(1)%i,ptilde,Bmunu)
        call CastSMEmunuToSME(ptilde,Bmunu,smeB)
! Otherwise just the SME is requested:
      else
        call CalcB(ptilde,smeB)
      end if
      call CalcStCirtFF(p,phat,ptilde,smeB,Bmunu, &
                        radi,radr,radt,           &
                        StCirtCont)
      StCirtterm = StCirtterm + Cirs%term(iterm)%symfact*StCirtCont
    end if
  end do
!
end subroutine PickStCirt
!
subroutine PickStCSirt(p,phat,ptilde,Bij,Bmunuij, &
                       emitir,radi,radr,radt,     &
                       CSirt,StCSirtterm)
use regions
use particles
use misc
use nnlo_subtractions
use nnlo_mom_maps
implicit none
!
  type(particle), dimension(:), intent(in) :: p
  type(particle), dimension(:), intent(out) :: phat,ptilde
  real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
  real(kind(1d0)) , dimension(0:,0:,:,:) , intent(inout) :: Bmunuij
  integer , intent(in) :: emitir,radi,radr,radt
  type(subterms) , intent(in) :: CSirt
  real(kind(1d0)) , intent(out) :: StCSirtterm
!
  integer :: iterm
  integer :: rad_i,rad_r,rad_ir
  real(kind(1d0)) :: ytQ
  real(kind(1d0)) :: StCSirtCont
!
!
  interface 
    subroutine CalcBij(p,Bij)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Bij
!
    end subroutine CalcBij
!
    subroutine CalcBmunuij(ileg,p,Bmunuij)
    use particles
    implicit none
!
      integer , intent(in) :: ileg
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(0:,0:,:,:) , intent(out) :: Bmunuij
!
    end subroutine CalcBmunuij
  end interface
!
  StCSirtterm = 0d0
! Position of parton i and r among hat momenta: 
  rad_i = radi
  rad_r = radr
  if (rad_i.gt.radt) rad_i = rad_i - 1
  if (rad_r.gt.radt) rad_r = rad_r - 1
! Position for the \hat{ir} parton 
  rad_ir = min(rad_i,rad_r)
  do iterm=1,CSirt%numterm
    if ((CSirt%term(iterm)%rad(1)%i.eq.radi).and.  &
        (CSirt%term(iterm)%rad(2)%i.eq.radr).and. &
        (CSirt%term(iterm)%rad(3)%i.eq.radt)) then
      call MapMomSr_A12(p,phat,       &
                        radt,         &
                        Q,Q2,yiQ_arr, &
                        ytQ,lambda_t)
      call CalcSubhatInvariants(phat)
      call MapMomCirFF_A12(phat,ptilde,                &
                           rad_i,rad_r,                &
                           Q,Q2,yihatQ_arr,sirhat_arr, &
                           sir,zhati,zhatr,            &
                           kthati,kthatr,              &
                           alphairt)
! Obtaining the underlying Born SME:
! If the emitter is a gluon we always need spin- and color-correlated SME:
      if (CSirt%term(iterm)%emit(1)%ID.eq.0) then
        call CalcBmunuij(CSirt%term(iterm)%emit(1)%i,ptilde,Bmunuij)
! Otherwise only the color-correlated SME is requested:
      else
        call CalcBij(ptilde,Bij)
      end if
      call CalcStCSirtFF(p,phat,ptilde,Bij,Bmunuij, &
                         rad_ir,radi,radr,radt,     &
                         StCSirtCont)
      StCSirtterm = StCSirtterm + CSirt%term(iterm)%symfact*StCSirtCont
    end if
  end do
!
end subroutine PickStCSirt
!
subroutine PickStCirtCSirt(p,phat,ptilde,         &
                           emitir,radi,radr,radt, &
                           CSirt,StCirtCSirtterm)
use regions
use particles
use misc
use nnlo_subtractions
use nnlo_mom_maps
implicit none
!
  type(particle), dimension(:), intent(in) :: p
  type(particle), dimension(:), intent(out) :: phat,ptilde
  integer , intent(in) :: emitir,radi,radr,radt
  type(subterms) , intent(in) :: CSirt
  real(kind(1d0)) , intent(out) :: StCirtCSirtterm
!
  integer :: iterm
  integer :: rad_i,rad_r,rad_ir
  real(kind(1d0)) :: StCirtCSirtCont
  real(kind(1d0)) :: smeB
  real(kind(1d0)) :: ytQ
  real(kind(1d0)) , dimension(0:3,0:3) :: Bmunu
!
!
  interface 
    subroutine CalcB(p,smeB)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: smeB
!
    end subroutine CalcB
!
    subroutine CalcBmunu(ileg,p,Bmunu)
    use particles
    implicit none
!
      integer , intent(in) :: ileg
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(0:3,0:3) , intent(out) :: Bmunu
!
    end subroutine CalcBmunu
  end interface
!
  StCirtCSirtterm = 0d0
! Position of parton i and r among hat momenta: 
  rad_i = radi
  rad_r = radr
  if (rad_i.gt.radt) rad_i = rad_i - 1
  if (rad_r.gt.radt) rad_r = rad_r - 1
! Position for the \hat{ir} parton 
  rad_ir = min(rad_i,rad_r)
  do iterm=1,CSirt%numterm
    if ((CSirt%term(iterm)%rad(1)%i.eq.radi).and.  &
        (CSirt%term(iterm)%rad(2)%i.eq.radr).and. &
        (CSirt%term(iterm)%rad(3)%i.eq.radt)) then
      call MapMomSr_A12(p,phat,       &
                        radt,         &
                        Q,Q2,yiQ_arr, &
                        ytQ,lambda_t)
      call CalcSubhatInvariants(phat)
      call MapMomCirFF_A12(phat,ptilde,                &
                           rad_i,rad_r,                &
                           Q,Q2,yihatQ_arr,sirhat_arr, &
                           sir,zhati,zhatr,            &
                           kthati,kthatr,              &
                           alphairt)
! Obtaining the underlying Born SME:
! If the emitter is a gluon we always need spin-correlated SME:
      if (CSirt%term(iterm)%emit(1)%ID.eq.0) then
        call CalcBmunu(CSirt%term(iterm)%emit(1)%i,ptilde,Bmunu)
! Otherwise only the SME is requested:
      else
        call CalcB(ptilde,smeB)
      end if
      call CalcStCirtCSirtFF(p,phat,ptilde,smeB,Bmunu, &
                             rad_ir,radi,radr,radt,    &
                             StCirtCSirtCont)
      StCirtCSirtterm = StCirtCSirtterm + CSirt%term(iterm)%symfact*StCirtCSirtCont
    end if
  end do
!
end subroutine PickStCirtCSirt
!
subroutine PickStCirtSrt(p,phat,ptilde,          &
                         emitirt,radi,radr,radt, &
                         Srt,StCirtSrtterm)
use regions
use particles
use misc
use nnlo_subtractions
use nnlo_mom_maps
implicit none
!
  type(particle), dimension(:), intent(in) :: p
  type(particle), dimension(:), intent(out) :: phat,ptilde
  integer , intent(in) :: emitirt,radi,radr,radt
  type(subterms) , intent(in) :: Srt
  real(kind(1d0)) , intent(out) :: StCirtSrtterm
!
  integer :: iterm
  integer :: rad_i,rad_r,rad1,rad2
  real(kind(1d0)) :: smeB
  real(kind(1d0)) :: yrhatQ
  real(kind(1d0)) :: ytQ
  real(kind(1d0)) :: StCirtSrtCont
!
!
  interface 
    subroutine CalcB(p,smeB)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: smeB
!
    end subroutine CalcB
  end interface
!
  StCirtSrtterm = 0d0
! The position of parton i and r can change as we go
! to hat momenta if parton t precedes at least one of them:
  rad_i = radi
  rad_r = radr
! This is needed since the indices are not ordered due to the
! presence of the St operator:
  rad1 = min(radr,radt)
  rad2 = max(radr,radt)
  if (radt.lt.radi) rad_i = rad_i - 1
  if (radt.lt.radr) rad_r = rad_r - 1
  do iterm=1,Srt%numterm
    if ((Srt%term(iterm)%rad(1)%i.eq.rad1).and.  &
        (Srt%term(iterm)%rad(2)%i.eq.rad2)) then
      call MapMomSr_A12(p,phat,       &
                        radt,         &
                        Q,Q2,yiQ_arr, &
                        ytQ,lambda_t)
      call CalcSubhatInvariants(phat)
      call MapMomSr_A12(phat,ptilde,     &
                        rad_r,           &
                        Q,Q2,yihatQ_arr, &
                        yrhatQ,lambda_r)
! Obtaining the underlying Born SME:
      call CalcB(ptilde,smeB)
      call CalcStCirtSrtFF(p,phat,ptilde,smeB,         &
                           rad_i,rad_r,radi,radr,radt, &
                           StCirtSrtCont)
      StCirtSrtterm = StCirtSrtterm + Srt%term(iterm)%symfact*StCirtSrtCont
    end if
  end do
!
end subroutine PickStCirtSrt
!
subroutine PickStCSirtSrt(p,phat,ptilde,Bij,     &
                          emitir,radi,radr,radt, &
                          Srt,StCSirtSrtterm)
use regions
use particles
use misc
use nnlo_subtractions
use nnlo_mom_maps
implicit none
!
  type(particle), dimension(:), intent(in) :: p
  type(particle), dimension(:), intent(out) :: phat,ptilde
  real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
  integer , intent(in) :: emitir,radi,radr,radt
  type(subterms) , intent(in) :: Srt
  real(kind(1d0)) , intent(out) :: StCSirtSrtterm
!
  integer :: iterm
  integer :: rad_ir,rad_i,rad_r,rad1,rad2
  real(kind(1d0)) :: ytQ
  real(kind(1d0)) :: yrhatQ
  real(kind(1d0)) :: StCSirtSrtCont
!
!
  interface 
    subroutine CalcBij(p,Bij)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Bij
!
    end subroutine CalcBij
  end interface
!
  StCSirtSrtterm = 0d0
! The position of parton i and r can change as we go
! to hat momenta if parton t precedes at least one of them:
  rad_i = radi
  rad_r = radr
! This is needed since the indices are not ordered due to the
! presence of the St operator:
  rad1 = min(radr,radt)
  rad2 = max(radr,radt)
  if (radt.lt.radi) rad_i = rad_i - 1
  if (radt.lt.radr) rad_r = rad_r - 1
! Position for the \hat{ir} parton 
  rad_ir = min(rad_i,rad_r)
  do iterm=1,Srt%numterm
    if ((Srt%term(iterm)%rad(1)%i.eq.rad1).and.  &
        (Srt%term(iterm)%rad(2)%i.eq.rad2)) then
      call MapMomSr_A12(p,phat,       &
                        radt,         &
                        Q,Q2,yiQ_arr, &
                        ytQ,lambda_t)
      call CalcSubhatInvariants(phat)
      call MapMomSr_A12(phat,ptilde,     &
                        rad_r,           &
                        Q,Q2,yihatQ_arr, &
                        yrhatQ,lambda_r)
! Obtaining the underlying Born color-correlated SME:
      call CalcBij(ptilde,Bij)
      call CalcStCSirtSrtFF(p,phat,ptilde,Bij,                 &
                            rad_ir,rad_i,rad_r,radi,radr,radt, &
                            StCSirtSrtCont)
      StCSirtSrtterm = StCSirtSrtterm + Srt%term(iterm)%symfact*StCSirtSrtCont
    end if
  end do
!
end subroutine PickStCSirtSrt
!
subroutine PickStCirtCSirtSrt(p,phat,ptilde,          &
                              emitirt,radi,radr,radt, &
                              Srt,StCirtCSirtSrtterm)
use regions
use particles
use misc
use nnlo_subtractions
use nnlo_mom_maps
implicit none
!
  type(particle), dimension(:), intent(in) :: p
  type(particle), dimension(:), intent(out) :: phat,ptilde
  integer , intent(in) :: emitirt,radi,radr,radt
  type(subterms) , intent(in) :: Srt
  real(kind(1d0)) , intent(out) :: StCirtCSirtSrtterm
!
  integer :: iterm
  integer :: rad_ir,rad_i,rad_r,rad1,rad2
  real(kind(1d0)) :: smeB
  real(kind(1d0)) :: ytQ
  real(kind(1d0)) :: yrhatQ
  real(kind(1d0)) :: StCirtCSirtSrtCont
!
!
  interface 
    subroutine CalcB(p,smeB)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: smeB
!
    end subroutine CalcB
  end interface
!
  StCirtCSirtSrtterm = 0d0
! The position of parton i and r can change as we go
! to hat momenta if parton t precedes at least one of them:
  rad_i = radi
  rad_r = radr
! This is needed since the indices are not ordered due to the
! presence of the St operator:
  rad1 = min(radr,radt)
  rad2 = max(radr,radt)
  if (radt.lt.radi) rad_i = rad_i - 1
  if (radt.lt.radr) rad_r = rad_r - 1
! Position for the \hat{ir} parton 
  rad_ir = min(rad_i,rad_r)
  do iterm=1,Srt%numterm
    if ((Srt%term(iterm)%rad(1)%i.eq.rad1).and.  &
        (Srt%term(iterm)%rad(2)%i.eq.rad2)) then
      call MapMomSr_A12(p,phat,       &
                        radt,         &
                        Q,Q2,yiQ_arr, &
                        ytQ,lambda_t)
      call CalcSubhatInvariants(phat)
      call MapMomSr_A12(phat,ptilde,     &
                        rad_r,           &
                        Q,Q2,yihatQ_arr, &
                        yrhatQ,lambda_r)
! Obtaining the underlying Born color-correlated SME:
      call CalcB(ptilde,smeB)
      call CalcStCirtCSirtSrtFF(p,phat,ptilde,smeB,         &
                                rad_i,rad_r,radi,radr,radt, &
                                StCirtCSirtSrtCont)
      StCirtCSirtSrtterm = StCirtCSirtSrtterm + Srt%term(iterm)%symfact*StCirtCSirtSrtCont
    end if
  end do
!
end subroutine PickStCirtCSirtSrt
!
subroutine PickStSrt(p,phat,ptilde,Bij,Bijkl, &
                     radr,radt,Srt,StSrtterm)
use regions
use particles
use misc
use nnlo_subtractions
use nnlo_mom_maps
implicit none
!
  type(particle), dimension(:), intent(in) :: p
  type(particle), dimension(:), intent(out) :: phat,ptilde
  real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
  real(kind(1d0)) , dimension(:,:,:,:) , intent(inout) :: Bijkl
  integer , intent(in) :: radr,radt
  type(subterms) , intent(in) :: Srt
  real(kind(1d0)) , intent(out) :: StSrtterm
!
  integer :: iterm
  integer :: rad_r,rad1,rad2
  real(kind(1d0)) :: ytQ
  real(kind(1d0)) :: yrhatQ
  real(kind(1d0)) :: StSrtCont
!
!
  interface 
    subroutine CalcBijkl(p,Bijkl)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:,:,:,:) , intent(out) :: Bijkl
!
    end subroutine CalcBijkl
  end interface
!
  StSrtterm = 0d0
! This is needed since the indices are not ordered due to the
! presence of the St operator:
  rad1 = min(radr,radt)
  rad2 = max(radr,radt)
! The position of parton r can change in the hatted configuration:
  rad_r = radr
  if (radt.lt.radr) rad_r = rad_r - 1
  do iterm=1,Srt%numterm
    if ((Srt%term(iterm)%rad(1)%i.eq.rad1).and.  &
        (Srt%term(iterm)%rad(2)%i.eq.rad2)) then
      call MapMomSr_A12(p,phat,       &
                        radt,         &
                        Q,Q2,yiQ_arr, &
                        ytQ,lambda_t)
      call CalcSubhatInvariants(phat)
      call MapMomSr_A12(phat,ptilde,     &
                        rad_r,           &
                        Q,Q2,yihatQ_arr, &
                        yrhatQ,lambda_r)
! Obtaining the underlying Born doubly color-correlated SME:
      call CalcBijkl(ptilde,Bijkl)
! The color-correlated SME can be conveniently obtained from the
! doubly color-correlated one:
      call CastSMEijklToSMEij(.false.,p,Bijkl,Bij)
      call CalcStSrtFF(p,phat,ptilde,Bij,Bijkl, &
                       rad_r,radr,radt,         &
                       StSrtCont)
      StSrtterm = StSrtterm + Srt%term(iterm)%symfact*StSrtCont
    end if
  end do
!
end subroutine PickStSrt
!
subroutine PickCitStCirt(p,phat,ptilde,                 &
                         radi,radr,radt,rad1,rad2,rad3, &
                         Cirt,CitStCirtterm)
use regions
use particles
use misc
use nnlo_subtractions
use nnlo_mom_maps
implicit none
!
  type(particle), dimension(:), intent(in) :: p
  type(particle), dimension(:), intent(out) :: phat,ptilde
  integer , intent(in) :: radi,radr,radt,rad1,rad2,rad3
  type(subterms) , intent(in) :: Cirt
  real(kind(1d0)) , intent(out) :: CitStCirtterm
!
  integer :: iterm
  integer :: rad_ir,rad_i,rad_r
  real(kind(1d0)) :: smeB
  real(kind(1d0)) :: ytQ
  real(kind(1d0)) , dimension(0:3,0:3) :: Bmunu
  real(kind(1d0)) :: CitStCirtCont
!
!
  interface 
    subroutine CalcB(p,smeB)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: smeB
!
    end subroutine CalcB
!
    subroutine CalcBmunu(ileg,p,Bmunu)
    use particles
    implicit none
!
      integer , intent(in) :: ileg
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(0:3,0:3) , intent(out) :: Bmunu
!
    end subroutine CalcBmunu
  end interface
!
  CitStCirtterm = 0d0
! Position of parton i and r among hat momenta: 
  rad_i = radi
  rad_r = radr
  if (rad_i.gt.radt) rad_i = rad_i - 1
  if (rad_r.gt.radt) rad_r = rad_r - 1
! Position for the \hat{ir} parton 
  rad_ir = min(rad_i,rad_r)
  do iterm=1,Cirt%numterm
    if ((Cirt%term(iterm)%rad(1)%i.eq.rad1).and.  &
        (Cirt%term(iterm)%rad(2)%i.eq.rad2).and. &
        (Cirt%term(iterm)%rad(3)%i.eq.rad3)) then
      call MapMomSr_A12(p,phat,       &
                        radt,         &
                        Q,Q2,yiQ_arr, &
                        ytQ,lambda_t)
      call CalcSubhatInvariants(phat)
      call MapMomCirFF_A12(phat,ptilde,                &
                           rad_i,rad_r,                &
                           Q,Q2,yihatQ_arr,sirhat_arr, &
                           sir,zhati,zhatr,            &
                           kthati,kthatr,              &
                           alphairt)
! Obtaining the underlying Born SME:
! If the emitter is a gluon we always need spin-correlated SME:
!      if (Cirt%term(iterm)%emit(1)%ID.eq.0) then
! This form should be used because of the iterated mapping:
      if (ptilde(rad_ir)%flv.eq.0) then
        call CalcBmunu(Cirt%term(iterm)%emit(1)%i,ptilde,Bmunu)
        call CastSMEmunuToSME(ptilde,Bmunu,smeB)
! Otherwise just the SME is requested:
      else
        call CalcB(ptilde,smeB)
      end if
      call CalcCitStCirtFF(p,phat,ptilde,smeB,Bmunu,   &
                           rad_i,rad_r,radi,radr,radt, &
                           CitStCirtCont)
      CitStCirtterm = CitStCirtterm + Cirt%term(iterm)%symfact*CitStCirtCont
    end if
  end do
!
end subroutine PickCitStCirt
!
subroutine PickCktStCSirt(p,phat,ptilde,       &
                          radi,radr,radk,radt, &
                          CSirt,CktStCSirtterm)
use regions
use particles
use misc
use nnlo_subtractions
use nnlo_mom_maps
implicit none
!
  type(particle), dimension(:), intent(in) :: p
  type(particle), dimension(:), intent(out) :: phat,ptilde
  integer , intent(in) :: radi,radr,radk,radt
  type(subterms) , intent(in) :: CSirt
  real(kind(1d0)) , intent(out) :: CktStCSirtterm
!
  integer :: iterm
  integer :: rad_ir,rad_i,rad_r
  real(kind(1d0)) :: smeB
  real(kind(1d0)) :: ytQ
  real(kind(1d0)) , dimension(0:3,0:3) :: Bmunu
  real(kind(1d0)) :: CktStCSirtCont
!
!
  interface 
    subroutine CalcB(p,smeB)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: smeB
!
    end subroutine CalcB
!
    subroutine CalcBmunu(ileg,p,Bmunu)
    use particles
    implicit none
!
      integer , intent(in) :: ileg
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(0:3,0:3) , intent(out) :: Bmunu
!
    end subroutine CalcBmunu
  end interface
!
  CktStCSirtterm = 0d0
! Position of parton i and r among hat momenta: 
  rad_i = radi
  rad_r = radr
  if (rad_i.gt.radt) rad_i = rad_i - 1
  if (rad_r.gt.radt) rad_r = rad_r - 1
! Position for the \hat{ir} parton 
  rad_ir = min(rad_i,rad_r)
  do iterm=1,CSirt%numterm
    if ((CSirt%term(iterm)%rad(1)%i.eq.radi).and.  &
        (CSirt%term(iterm)%rad(2)%i.eq.radr).and. &
        (CSirt%term(iterm)%rad(3)%i.eq.radt)) then
      call MapMomSr_A12(p,phat,       &
                        radt,         &
                        Q,Q2,yiQ_arr, &
                        ytQ,lambda_t)
      call CalcSubhatInvariants(phat)
      call MapMomCirFF_A12(phat,ptilde,                &
                           rad_i,rad_r,                &
                           Q,Q2,yihatQ_arr,sirhat_arr, &
                           sir,zhati,zhatr,            &
                           kthati,kthatr,              &
                           alphairt)
! Obtaining the underlying Born SME:
! If the emitter is a gluon we always need spin-correlated SME:
      if (CSirt%term(iterm)%emit(1)%ID.eq.0) then
        call CalcBmunu(CSirt%term(iterm)%emit(1)%i,ptilde,Bmunu)
        call CastSMEmunuToSME(ptilde,Bmunu,smeB)
! Otherwise just the SME is requested:
      else
        call CalcB(ptilde,smeB)
      end if
      call CalcCktStCSirtFF(p,phat,ptilde,smeB,Bmunu, &
                            radi,radr,radk,radt,      &
                            CktStCSirtCont)
      CktStCSirtterm = CktStCSirtterm + CSirt%term(iterm)%symfact*CktStCSirtCont
    end if
  end do
!
end subroutine PickCktStCSirt
!
subroutine PickCktStCSirtSrt(p,phat,ptilde,       &
                             radi,radr,radk,radt, &
                             Srt,CktStCSirtSrtterm)
use regions
use particles
use misc
use nnlo_subtractions
use nnlo_mom_maps
implicit none
!
  type(particle), dimension(:), intent(in) :: p
  type(particle), dimension(:), intent(out) :: phat,ptilde
  integer , intent(in) :: radi,radr,radk,radt
  type(subterms) , intent(in) :: Srt
  real(kind(1d0)) , intent(out) :: CktStCSirtSrtterm
!
  integer :: iterm
  integer :: rad_i,rad_r
  real(kind(1d0)) :: smeB
  real(kind(1d0)) :: ytQ
  real(kind(1d0)) :: yrhatQ
  real(kind(1d0)) :: CktStCSirtSrtCont
!
!
  interface 
    subroutine CalcB(p,smeB)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: smeB
!
    end subroutine CalcB
  end interface
!
  CktStCSirtSrtterm = 0d0
! Position of parton i and r among hat momenta: 
  rad_i = radi
  rad_r = radr
  if (rad_i.gt.radt) rad_i = rad_i - 1
  if (rad_r.gt.radt) rad_r = rad_r - 1
  do iterm=1,Srt%numterm
    if (((Srt%term(iterm)%rad(1)%i.eq.radr).and. &
         (Srt%term(iterm)%rad(2)%i.eq.radt)).or. &
        ((Srt%term(iterm)%rad(1)%i.eq.radt).and. &
         (Srt%term(iterm)%rad(2)%i.eq.radr))) then
      call MapMomSr_A12(p,phat,       &
                        radt,         &
                        Q,Q2,yiQ_arr, &
                        ytQ,lambda_t)
      call CalcSubhatInvariants(phat)
      call MapMomSr_A12(phat,ptilde,     &
                        rad_r,           &
                        Q,Q2,yihatQ_arr, &
                        yrhatQ,lambda_r)
! Obtaining the underlying Born SME:
      call CalcB(ptilde,smeB)
      call CalcCktStCSirtSrtFF(p,phat,ptilde,smeB,  &
                               rad_i,rad_r,         &
                               radi,radr,radk,radt, &
                               CktStCSirtSrtCont)
      CktStCSirtSrtterm = CktStCSirtSrtterm + Srt%term(iterm)%symfact*CktStCSirtSrtCont
    end if
  end do
!
end subroutine PickCktStCSirtSrt
!
subroutine PickCktStCkrtSrt(p,phat,ptilde,  &
                            radk,radr,radt, &
                            Srt,CktStCkrtSrtterm)
use regions
use particles
use misc
use nnlo_subtractions
use nnlo_mom_maps
implicit none
!
  type(particle), dimension(:), intent(in) :: p
  type(particle), dimension(:), intent(out) :: phat,ptilde
  integer , intent(in) :: radk,radr,radt
  type(subterms) , intent(in) :: Srt
  real(kind(1d0)) , intent(out) :: CktStCkrtSrtterm
!
  integer :: iterm
  integer :: rad_k,rad_r
  real(kind(1d0)) :: smeB
  real(kind(1d0)) :: ytQ
  real(kind(1d0)) :: yrhatQ
  real(kind(1d0)) :: CktStCkrtSrtCont
!
!
  interface 
    subroutine CalcB(p,smeB)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: smeB
!
    end subroutine CalcB
  end interface
!
  CktStCkrtSrtterm = 0d0
! Position of parton k and r among hat momenta: 
  rad_k = radk
  rad_r = radr
  if (rad_k.gt.radt) rad_k = rad_k - 1
  if (rad_r.gt.radt) rad_r = rad_r - 1
  do iterm=1,Srt%numterm
    if (((Srt%term(iterm)%rad(1)%i.eq.radr).and. &
         (Srt%term(iterm)%rad(2)%i.eq.radt)).or. &
        ((Srt%term(iterm)%rad(1)%i.eq.radt).and. &
         (Srt%term(iterm)%rad(2)%i.eq.radr))) then
      call MapMomSr_A12(p,phat,       &
                        radt,         &
                        Q,Q2,yiQ_arr, &
                        ytQ,lambda_t)
      call CalcSubhatInvariants(phat)
      call MapMomSr_A12(phat,ptilde,     &
                        rad_r,           &
                        Q,Q2,yihatQ_arr, &
                        yrhatQ,lambda_r)
! Obtaining the underlying Born SME:
      call CalcB(ptilde,smeB)
      call CalcCktStCkrtSrtFF(p,phat,ptilde,smeB, &
                              rad_k,rad_r,        &
                              radk,radr,radt,     &
                              CktStCkrtSrtCont)
      CktStCkrtSrtterm = CktStCkrtSrtterm + Srt%term(iterm)%symfact*CktStCkrtSrtCont
    end if
  end do
!
end subroutine PickCktStCkrtSrt
!
subroutine PickCrtStCkrtSrt(p,phat,ptilde,  &
                            radk,radr,radt, &
                            Srt,CrtStCkrtSrtterm)
use regions
use particles
use misc
use nnlo_subtractions
use nnlo_mom_maps
implicit none
!
  type(particle), dimension(:), intent(in) :: p
  type(particle), dimension(:), intent(out) :: phat,ptilde
  integer , intent(in) :: radk,radr,radt
  type(subterms) , intent(in) :: Srt
  real(kind(1d0)) , intent(out) :: CrtStCkrtSrtterm
!
  integer :: iterm
  integer :: rad_k,rad_r
  real(kind(1d0)) :: smeB
  real(kind(1d0)) :: ytQ
  real(kind(1d0)) :: yrhatQ
  real(kind(1d0)) :: CrtStCkrtSrtCont
!
!
  interface 
    subroutine CalcB(p,smeB)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: smeB
!
    end subroutine CalcB
  end interface
!
  CrtStCkrtSrtterm = 0d0
! Position of parton k and r among hat momenta: 
  rad_k = radk
  rad_r = radr
  if (rad_k.gt.radt) rad_k = rad_k - 1
  if (rad_r.gt.radt) rad_r = rad_r - 1
  do iterm=1,Srt%numterm
    if (((Srt%term(iterm)%rad(1)%i.eq.radr).and. &
         (Srt%term(iterm)%rad(2)%i.eq.radt)).or. &
        ((Srt%term(iterm)%rad(1)%i.eq.radt).and. &
         (Srt%term(iterm)%rad(2)%i.eq.radr))) then
      call MapMomSr_A12(p,phat,       &
                        radt,         &
                        Q,Q2,yiQ_arr, &
                        ytQ,lambda_t)
      call CalcSubhatInvariants(phat)
      call MapMomSr_A12(phat,ptilde,     &
                        rad_r,           &
                        Q,Q2,yihatQ_arr, &
                        yrhatQ,lambda_r)
! Obtaining the underlying Born SME:
      call CalcB(ptilde,smeB)
      call CalcCrtStCkrtSrtFF(p,phat,ptilde,smeB, &
                              rad_k,rad_r,        &
                              radk,radr,radt,     &
                              CrtStCkrtSrtCont)
      CrtStCkrtSrtterm = CrtStCkrtSrtterm + Srt%term(iterm)%symfact*CrtStCkrtSrtCont
    end if
  end do
!
end subroutine PickCrtStCkrtSrt
!
subroutine PickCktStSrt(p,phat,ptilde,Bij, &
                        radk,radr,radt,    &
                        Srt,CktStSrtterm)
use regions
use particles
use misc
use nnlo_subtractions
use nnlo_mom_maps
implicit none
!
  type(particle), dimension(:), intent(in) :: p
  type(particle), dimension(:), intent(out) :: phat,ptilde
  real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
  integer , intent(in) :: radk,radr,radt
  type(subterms) , intent(in) :: Srt
  real(kind(1d0)) , intent(out) :: CktStSrtterm
!
  integer :: iterm
  integer :: rad_r
  real(kind(1d0)) :: ytQ
  real(kind(1d0)) :: yrhatQ
  real(kind(1d0)) :: CktStSrtCont
!
!
  interface 
    subroutine CalcBij(p,Bij)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Bij
!
    end subroutine CalcBij
  end interface
!
  CktStSrtterm = 0d0
! Position of r among hat momenta: 
  rad_r = radr
  if (rad_r.gt.radt) rad_r = rad_r - 1
  do iterm=1,Srt%numterm
    if (((Srt%term(iterm)%rad(1)%i.eq.radr).and. &
         (Srt%term(iterm)%rad(2)%i.eq.radt)).or. &
        ((Srt%term(iterm)%rad(1)%i.eq.radt).and. &
         (Srt%term(iterm)%rad(2)%i.eq.radr))) then
      call MapMomSr_A12(p,phat,       &
                        radt,         &
                        Q,Q2,yiQ_arr, &
                        ytQ,lambda_t)
      call CalcSubhatInvariants(phat)
      call MapMomSr_A12(phat,ptilde,     &
                        rad_r,           &
                        Q,Q2,yihatQ_arr, &
                        yrhatQ,lambda_r)
! Obtaining the underlying color-correlated Born SME:
      call CalcBij(ptilde,Bij)
      call CalcCktStSrtFF(p,phat,ptilde,Bij,    &
                          rad_r,radk,radr,radt, &
                          CktStSrtCont)
      CktStSrtterm = CktStSrtterm + Srt%term(iterm)%symfact*CktStSrtCont
    end if
  end do
!
end subroutine PickCktStSrt
!
subroutine PickCrtStSrt(p,phat,ptilde,Bij, &
                        radr,radt,         &
                        Srt,CrtStSrtterm)
use regions
use particles
use misc
use nnlo_subtractions
use nnlo_mom_maps
implicit none
!
  type(particle), dimension(:), intent(in) :: p
  type(particle), dimension(:), intent(out) :: phat,ptilde
  real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
  integer , intent(in) :: radr,radt
  type(subterms) , intent(in) :: Srt
  real(kind(1d0)) , intent(out) :: CrtStSrtterm
!
  integer :: iterm
  integer :: rad_r
  real(kind(1d0)) :: ytQ
  real(kind(1d0)) :: yrhatQ
  real(kind(1d0)) :: CrtStSrtCont
!
!
  interface 
    subroutine CalcBij(p,Bij)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Bij
!
    end subroutine CalcBij
  end interface
!
  CrtStSrtterm = 0d0
! Position of r among hat momenta: 
  rad_r = radr
  if (rad_r.gt.radt) rad_r = rad_r - 1
  do iterm=1,Srt%numterm
    if (((Srt%term(iterm)%rad(1)%i.eq.radr).and. &
         (Srt%term(iterm)%rad(2)%i.eq.radt)).or. &
        ((Srt%term(iterm)%rad(1)%i.eq.radt).and. &
         (Srt%term(iterm)%rad(2)%i.eq.radr))) then
      call MapMomSr_A12(p,phat,       &
                        radt,         &
                        Q,Q2,yiQ_arr, &
                        ytQ,lambda_t)
      call CalcSubhatInvariants(phat)
      call MapMomSr_A12(phat,ptilde,     &
                        rad_r,           &
                        Q,Q2,yihatQ_arr, &
                        yrhatQ,lambda_r)
! Obtaining the underlying color-correlated Born SME:
      call CalcBij(ptilde,Bij)
      call CalcCrtStSrtFF(p,phat,ptilde,Bij, &
                          rad_r,radr,radt,   &
                          CrtStSrtCont)
      CrtStSrtterm = CrtStSrtterm + Srt%term(iterm)%symfact*CrtStSrtCont
    end if
  end do
!
end subroutine PickCrtStSrt
!
subroutine Calcdmunu(k,n,dmunu)
use momenta
implicit none
!
  type(mom) , intent(in) :: k,n
  real(kind(1d0)) , dimension(0:3,0:3) , intent(out) :: dmunu
!
  real(kind(1d0)) :: kn
!
!
  kn = k*n
!
  dmunu(0,0) = -1 + (k%E*n%E + k%E*n%E)/kn
  dmunu(1,1) = +1 + (k%px*n%px + k%px*n%px)/kn
  dmunu(2,2) = +1 + (k%py*n%py + k%py*n%py)/kn
  dmunu(3,3) = +1 + (k%pz*n%pz + k%pz*n%pz)/kn
!
  dmunu(0,1) = (k%E*n%px + k%px*n%E)/kn
  dmunu(1,0) = dmunu(0,1)
  dmunu(0,2) = (k%E*n%py + k%py*n%E)/kn
  dmunu(2,0) = dmunu(0,2)
  dmunu(0,3) = (k%E*n%pz + k%pz*n%E)/kn
  dmunu(3,0) = dmunu(0,3)
  dmunu(1,2) = (k%px*n%py + k%py*n%px)/kn
  dmunu(2,1) = dmunu(1,2)
  dmunu(1,3) = (k%px*n%pz + k%pz*n%px)/kn
  dmunu(3,1) = dmunu(1,3)
  dmunu(2,3) = (k%py*n%pz + k%pz*n%py)/kn
  dmunu(3,2) = dmunu(2,3)
!
end subroutine Calcdmunu
!
! This routine takes a transverse momentum and changes it to
! be exactly orthogonal to p using a 4-vector Q:
subroutine ChangeKt(kt,p,Q)
use momenta
implicit none
!
  type(mom) , intent(inout) :: kt
  type(mom) , intent(in) :: p,Q
!
  real(kind(1d0)) :: ktp,pQ
!
!
  ktp = kt*p
  pQ  = p*Q
  kt  = kt - ktp/pQ*Q
!
end subroutine ChangeKt
!
subroutine MakeMassless(kin,k1,k2)
use momenta
implicit none
!
  type(mom) , intent(in) :: kin
  type(mom) , intent(out) :: k1,k2
!
  real(kind(1d0)) :: K
  type(mom) :: n
!
!
  call PrintMom(kin)
!
  K = sqrt(-kin*kin)
  k1 = kin/2d0
  k2 = kin/2d0
  k1%E = K/2d0
  k2%E = -K/2d0
!
  print *,"Original momentum: "
  call PrintMom(kin)
  print *,"k1: "
  call PrintMom(k1)
  print *,"k2: "
  call PrintMom(k2)
  call PrintMom(k1 + k2)
  print *,"k1.k1: ",k1*k1
  print *,"k2.k2: ",k2*k2
!
end subroutine MakeMassless
!
subroutine CalcCgggSG(p,ptilde,radi,radr,rads,Cirs)
use particles
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p,ptilde
  integer , intent(in) :: radi,radr,rads
  real(kind(1d0)) , intent(out) :: Cirs
!
  integer :: emit
  real(kind(1d0)) :: CirsCont
!
  interface
    subroutine CalcCgggSGORD(emit,p_i,p_r,p_s,ptilde,CirsCont)
    use particles
    implicit none
!
      integer , intent(in) :: emit
      type(mom) , intent(in) :: p_i,p_r,p_s
      type(particle) , dimension(:) , intent(in) :: ptilde
      real(kind(1d0)) , intent(out) :: CirsCont
!
    end subroutine CalcCgggSGORD
  end interface
!
  Cirs = 0
!
  emit = min(radi,radr,rads)
!
! irs:
!  print *,"123: "
  call CalcCgggSGORD(emit,p(radi)%p,p(radr)%p,p(rads)%p,ptilde,CirsCont)
  Cirs = Cirs + CirsCont
! rsi:
!  print *,"231: "
  call CalcCgggSGORD(emit,p(radr)%p,p(rads)%p,p(radi)%p,ptilde,CirsCont)
  Cirs = Cirs + CirsCont
! sir:
!  print *,"312: "
  call CalcCgggSGORD(emit,p(rads)%p,p(radi)%p,p(radr)%p,ptilde,CirsCont)
  Cirs = Cirs + CirsCont
! sri:
!  print *,"321: "
  call CalcCgggSGORD(emit,p(rads)%p,p(radr)%p,p(radi)%p,ptilde,CirsCont)
  Cirs = Cirs + CirsCont
! isr:
!  print *,"132: "
  call CalcCgggSGORD(emit,p(radi)%p,p(rads)%p,p(radr)%p,ptilde,CirsCont)
  Cirs = Cirs + CirsCont
! ris:
!  print *,"213: "
  call CalcCgggSGORD(emit,p(radr)%p,p(radi)%p,p(rads)%p,ptilde,CirsCont)
  Cirs = Cirs + CirsCont
!
  Print *,"CirsSG:   ",Cirs
!
end subroutine CalcCgggSG
!
subroutine CalcCgggSGORD(emit,p_i,p_r,p_s,ptilde,CirsCont)
use particles
use QCDparams
use math
implicit none
!
  integer , intent(in) :: emit
  type(mom) , intent(in) :: p_i,p_r,p_s
  type(particle) , dimension(:) , intent(in) :: ptilde
  real(kind(1d0)) , intent(out) :: CirsCont
!
  real(kind(1d0)) :: sir,sis,srs,sirs
  real(kind(1d0)) :: zirs,zris,zsir
  real(kind(1d0)) :: tirs
  real(kind(1d0)) :: Q2
  real(kind(1d0)) :: pirs2,pirsQ
  real(kind(1d0)) :: siQ,srQ,ssQ
  real(kind(1d0)) :: alpirs
  real(kind(1d0)) :: zirsv,zrisv,zsirv
  real(kind(1d0)) :: kir2,ki2,kr2,ks2,krks2,krs2,kikr2
  type(mom) :: Q,p_irs,kti,ktr,kts,kir,krs
  real(kind(1d0)) , dimension(4) :: pirs_arr, &
                                    kir_arr,krs_arr,ktr_arr,kts_arr,ki_arr
  complex(kind(1d0)) :: kirPP,kirPM,kirMP,kirMM, &
                        krsPP,krsPM,krsMP,krsMM, &
                        krPP,krPM,krMP,krMM, &
                        ksPP,ksPM,ksMP,ksMM, &
                        krksPP,krksPM,krksMP,krksMM, &
                        kiPP,kiMM,kiPM,kiMP, &
                        kikrPP,kikrPM,kikrMP,kikrMM, &
                        VPP,VMM,VPM,VMP
  complex(kind(1d0)) , dimension(2,2) :: M2IJ
!
!
  interface
    subroutine CalcBmunu_hel(ileg,parts,Balbe)
    use particles
    implicit none
!
      integer , intent(in) :: ileg
      type(particle) , dimension(:) , intent(in) :: parts
      complex(kind(1d0)) , dimension(2,2) , intent(out) :: Balbe
!
    end subroutine CalcBmunu_hel
  end interface
!
  CirsCont = 0
!
  sir = 2*p_i*p_r
  sis = 2*p_i*p_s
  srs = 2*p_r*p_s
  sirs = sir + sis + srs
!  print *,"sir:  ",sir
!  print *,"sis:  ",sis
!  print *,"srs:  ",srs
!  print *,"sirs: ",sirs
!
  Q = ptilde(1)%p + ptilde(2)%p
  Q2 = Q*Q
  zirs = 2*p_i*Q &
      /(2*p_i*Q + 2*p_r*Q + 2*p_s*Q)
  zris = 2*p_r*Q &
     /(2*p_i*Q + 2*p_r*Q + 2*p_s*Q)
  zsir = 2*p_s*Q &
     /(2*p_i*Q + 2*p_r*Q + 2*p_s*Q)      
!  print *,"zirs: ",zirs
!  print *,"zris: ",zris
!  print *,"zsir: ",zsir
!
  tirs = 2*(zirs * srs - zris * sis)/(zirs + zris) &
       +(zirs - zris)/(zirs + zris) * sir
!  print *,"tirs: ",tirs
!
  siQ = 2*p_i*Q
  srQ = 2*p_r*Q
  ssQ = 2*p_s*Q
!
  pirs2 = (p_i + p_r + p_s)**2
  pirsQ = (p_i + p_r + p_s)*Q
  alpirs = (sqrt(pirsQ**2 - pirs2*Q2) - pirsQ)/Q2
  p_irs = (p_i + p_r + p_s + alpirs*Q)/(1 + alpirs)
!  print *,"p_irs: "
!  call PrintMom(p_irs)
!
  zirsv = (zirs + (sir + sis)/(alpirs*(siQ+srQ+ssQ)))
  zrisv = (zris + (sir + srs)/(alpirs*(siQ+srQ+ssQ)))
  zsirv = (zsir + (sis + srs)/(alpirs*(siQ+srQ+ssQ)))
!
  kti = zrisv * p_i - zirsv * p_r &
      + zsirv * p_i - zirsv * p_s &
      +0*(2*zirs*sirs-sir-sis)/sirs * p_irs
  ktr = zirsv * p_r - zrisv * p_i &
      + zsirv * p_r - zrisv * p_s &
      +0*(2*zris*sirs-sir-srs)/sirs * p_irs
  kts = zirsv * p_s - zsirv * p_i &
      + zrisv * p_s - zsirv * p_r &
      +0*(2*zsir*sirs-sis-srs)/sirs * p_irs
!
!  print *,"kt_i: "
!  call PrintMom(kti)
!  print *,"kt_r: "
!  call PrintMom(ktr)
!  print *,"kt_s: "
!  call PrintMom(kts)
!
  kir = kti/zirs - ktr/zris
  krs = ktr/zris - kts/zsir
!
  kir2 = -sir/(zirs * zris)
  ki2  = -zirs*((1-zirs)*sirs - srs)
  kr2  = -zris*((1-zris)*sirs - sis)
  ks2  = -zsir*((1-zsir)*sirs - sir)
  krks2 = -(zirs*(1-zirs) - zris*(1-zris) - zsir*(1-zsir)) &
        * sirs + zirs * srs - zris * sis - zsir * sir
  kikr2 = -(zsir*(1-zsir) - zirs*(1-zirs) - zris*(1-zris)) &
        * sirs + zsir * sir - zirs * srs - zris * sis
! Equivalent definition:
!  krks2 = srs + 2*sirs*zris*zsir - zsir*(sir+srs) - zris*(sis+srs)
  krs2 = kr2/zris**2 + ks2/zsir**2 - krks2/(zris*zsir)
!
!  print *,"kr2: ",kr2
!  print *,"ks2: ",ks2
!  print *,"krks2: ",krks2/2d0
!
  ki_arr  = kti
  kir_arr = kir
  krs_arr = krs
  ktr_arr = ktr
  kts_arr = kts
  pirs_arr = p_irs
!
  print *,"sirs,sir,sis,srs,zirs,zris,zsir: ",sirs,sir,sis,srs,zirs,zris,zsir
!
  call LtoH(kir_arr,pirs_arr,kirPP,kirPM,kirMP,kirMM)
  call LtoH(krs_arr,pirs_arr,krsPP,krsPM,krsMP,krsMM)
  call LtoH(ktr_arr,pirs_arr,krPP,krPM,krMP,krMM)
  call LtoH(kts_arr,pirs_arr,ksPP,ksPM,ksMP,ksMM)
  call LtoH(ki_arr,pirs_arr,kiPP,kiPM,kiMP,kiMM)
!
!  krksPP = zris * zsir/(krks2) &
!         * ( kr2/zris**2 * (krPP - krsPP) &
!         +ks2/zsir**2 * (ksPP - krsPP) &
!         +krks2/(zris*zsir) * krsPP)
!
  krksPP = zris * zsir/(2*ktr*kts) &
         * ( ktr**2/zris**2 * (krPP - krsPP) &
         +kts**2/zsir**2 * (ksPP - krsPP) &
         +2*ktr*kts/(zris*zsir) * krsPP)
  krksMM = krksPP
!  krksPM = zris * zsir/(krks2) &
!         * ( kr2/zris**2 * (krPM - krsPM) & 
!         +ks2/zsir**2 * (ksPM - krsPM) &
!         +krks2/(zris*zsir) * krsPM)
  krksPM = zris * zsir/(2*ktr*kts) &
         * ( ktr**2/zris**2 * (krPM - krsPM) & 
         +kts**2/zsir**2 * (ksPM - krsPM) &
         +2*ktr*kts/(zris*zsir) * krsPM)
  krksMP = conjg(krksPM)
!
  kikrPP = zirs * zris/(2*kti*ktr) &
         * ( kti**2/zirs**2 * (kiPP - kirPP) &
         +ktr**2/zris**2 * (krPP - kirPP) &
         +2*kti*ktr/(zirs*zris) * kirPP)
  kikrMM = kikrPP
!
  kikrPM = zirs * zris/(2*kti*ktr) &
         * ( kti**2/zirs**2 * (kiPM - kirPM) &
         +ktr**2/zris**2 * (krPM - kirPM) &
         +2*kti*ktr/(zirs*zris) * kirPM)
  kikrMP = conjg(kikrPM)
!
! Just to test the diagonal part:
!  kirPP = 0 ; kirPM = 0 ; kirMP = 0 ; kirMM = 0
!  krsPP = 0 ; krsPM = 0 ; krsMP = 0 ; krsMM = 0
!  krPP = 0 ; krPM = 0 ; krMP = 0 ; krMM = 0
!  ksPP = i ; ksPM = 0 ; ksMP = 0 ; ksMM = 0
!  krksPP = 0 ; krksPM = 0 ; krksMP = 0 ; krksMM = 0
!
!  kirPP = kikrPP ; kirPM = kikrPM ; kirMP = kikrMP ; kirMM = kikrMM
!  kir2 = -kikr2/zirs/zris
!
! New stuff:
!  kirPP = -kikrPP*kikr2/zirs/zris &
!        +  kiPP*ki2/zirs**2 &
!        +  krPP*kr2/zris**2
!  kirPM = -kikrPM*kikr2/zirs/zris &
!        +  kiPM*ki2/zirs**2 & 
!        +  krPM*kr2/zris**2
!  kirMP = -kikrMP*kikr2/zirs/zris &
!        +  kiMP*ki2/zirs**2 &
!        +  krMP*kr2/zris**2
!  kirMM = -kikrMM*kikr2/zirs/zris &
!        +  kiMM*ki2/zirs**2 & 
!        +  krMM*kr2/zris**2
!  kir2 = 1d0
!
!  krksPP = -1d0 ; krksPM = 0 ; krksMP = 0 ; krksMM = -1d0
!  kirPP = -1d0 ; kirPM = 0 ; kirMP = 0 ; kirMM = -1d0
!
  VPP = 1d0/(4d0*sir**2) &
      *(tirs**2 + 16d0*sirs*(zirs**2 * zris**2)/(zsir*(1d0-zsir)) &
      *kirPP * kir2) &
      +3d0/4d0 &
      -sirs/sir*1d0/zsir*((2d0*(1d0-zsir)+4d0*zsir**2)/(1d0-zsir) &
      -(1d0-2d0*zsir*(1d0-zsir))/(zirs*(1d0-zirs))) &
      +sirs/(sir*sis)*(2d0*zirs*(krPP * kr2 * (1d0-2d0*zsir) &
      /(zsir*(1d0-zsir)) &
      +ksPP * ks2 * (1d0-2d0*zris) &
      /(zris*(1d0-zris))) &
      -sirs/2d0*((4d0*zris*zsir+2d0*zirs*(1d0-zirs)-1d0) &
      /((1d0-zris)*(1d0-zsir)) &
      -(1d0-2d0*zirs*(1d0-zirs))/(zris*zsir)) &
      +krksPP * krks2 * (2d0*zris*(1d0-zris) &
      /(zsir*(1d0-zsir)) - 3d0))
!
!  VPP = VPP - (&
!       1d0/(4d0*sir**2) &
!      *(tirs**2 + 16d0*sirs*(zirs**2 * zris**2)/(zsir*(1d0-zsir)) &
!      *0 * kir2) &
!      +3d0/4d0 &
!      -sirs/sir*1d0/zsir*((2d0*(1d0-zsir)+4d0*zsir**2)/(1d0-zsir) &
!      -(1d0-2d0*zsir*(1d0-zsir))/(zirs*(1d0-zirs))) &
!      +sirs/(sir*sis)*(2d0*zirs*(0 * kr2 * (1d0-2d0*zsir) &
!      /(zsir*(1d0-zsir)) &
!      +0 * ks2 * (1d0-2d0*zris) &
!      /(zris*(1d0-zris))) &
!      -sirs/2d0*((4d0*zris*zsir+2d0*zirs*(1d0-zirs)-1d0) &
!      /((1d0-zris)*(1d0-zsir)) &
!      -(1d0-2d0*zirs*(1d0-zirs))/(zris*zsir)) &
!      +0 * krks2 * (2d0*zris*(1d0-zris) &
!      /(zsir*(1d0-zsir)) - 3d0)) &
!      )
!
  VMM = VPP
!
  VPM = 1d0/(4d0*sir**2)  &
      *(16d0*sirs*(zirs**2 * zris**2)/(zsir*(1d0-zsir)) &
      *kirPM * kir2) &
      +sirs/(sir*sis)*(2d0*zirs*(krPM * kr2 * (1d0-2d0*zsir) &
      /(zsir*(1d0-zsir)) &
      +ksPM * ks2 * (1d0-2d0*zris) &
      /(zris*(1d0-zris))) &
      +krksPM * krks2 * (2d0*zris*(1d0-zris) &
      /(zsir*(1d0-zsir)) - 3d0))
!
!  VPM = VPM - ( &
!      1d0/(4d0*sir**2)  &
!      *(16d0*sirs*(zirs**2 * zris**2)/(zsir*(1d0-zsir)) &
!      *0 * kir2) &
!      +sirs/(sir*sis)*(2d0*zirs*(0 * kr2 * (1d0-2d0*zsir) &
!      /(zsir*(1d0-zsir)) &
!      +0 * ks2 * (1d0-2d0*zris) &
!      /(zris*(1d0-zris))) &
!      +0 * krks2 * (2d0*zris*(1d0-zris) &
!      /(zsir*(1d0-zsir)) - 3d0)) &
!      )
!
  VMP = conjg(VPM)
!
!  VPM = -VPM
!
  VMP = conjg(VPM)
!
!  VPP = kirPP ; VMM = kirMM ; VPM = kirPM ; VMP = kirMP
!  VPP = krsPP ; VMM = krsMM ; VPM = krsPM ; VMP = krsMP
!
!  print *,"VPP: ",VPP
!  print *,"VMM: ",VMM
!  print *,"VPM: ",VPM
!  print *,"VMP: ",VMP
!
  call CalcBmunu_hel(emit,ptilde,M2IJ)
!
  CirsCont = M2IJ(1,1)*VPP + M2IJ(2,2)*VMM &
           + M2IJ(1,2)*VMP + M2IJ(2,1)*VPM
!  print *,"CirsCont: ",CirsCont
!
  CirsCont = (8d0*pi)**2/sirs**2*CirsCont*qcd_ca**2/6d0
!
  print *,"smeB from hel: ",M2IJ(1,1) + M2IJ(2,2)
!
end subroutine CalcCgggSGORD
!
! CrsCirs(ggg)
subroutine CalcCggCgggSG(p,phat,ptilde,radi,radr,rads,CrsCirs)
use particles
use QCDparams
use math
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p,phat,ptilde
  integer , intent(in) :: radi,radr,rads
  real(kind(1d0)) , intent(out) :: CrsCirs
!
  integer :: emit
  integer :: rad_rs,rad_i
  real(kind(1d0)) :: srs,zrs,zsr,Q2
  real(kind(1d0)) :: srQ,ssQ,prsQ,pitrstQ,alprs,alpitrst
  real(kind(1d0)) :: zrsv,zsrv
  real(kind(1d0)) :: sitQ,srstQ,sitrst
  real(kind(1d0)) :: zitrst,zrstit,zitrstv,zrstitv
  real(kind(1d0)) :: Ktrs2,Ktitrst2,skk,sitirst,srstirst,skp
  real(kind(1d0)) :: na,nb,nc,skn,spn
  real(kind(1d0)) , dimension(4) :: p_irs_arr,Ktrs_arr,Ktitrst_arr
  real(kind(1d0)) , dimension(0:3,0:3) :: Bmunu
  complex(kind(1d0)) :: gPP,gMM,gPM,gMP
  complex(kind(1d0)) :: kPP,kMM,kPM,kMP, &
                        VPP,VMM,VPM,VMP, &
                        khPP,khMM,khPM,khMP
  complex(kind(1d0)) , dimension(2,2) :: M2IJ
  type(mom) :: Q,p_i,p_r,p_s,phat_rs,phat_i,p_irs,Ktrs,Ktitrst
  type(mom) :: nirs
!
!
  interface
    subroutine CalcBmunu_hel(ileg,parts,Balbe)
    use particles
    implicit none
!
      integer , intent(in) :: ileg
      type(particle) , dimension(:) , intent(in) :: parts
      complex(kind(1d0)) , dimension(2,2) , intent(out) :: Balbe
!
    end subroutine CalcBmunu_hel
!
    subroutine CalcBmunu(ileg,p,Bmunu)
    use particles
    implicit none
!
      integer , intent(in) :: ileg
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(0:3,0:3) , intent(out) :: Bmunu
!
    end subroutine CalcBmunu
  end interface
!
  CrsCirs = 0
!
!  print *,"radi,radr,rads: ",radi,radr,rads
!
  emit = min(radi,radr,rads)
!
  rad_rs = min(radr,rads)
  rad_i = radi
  if (rad_i.gt.max(radr,rads)) rad_i = rad_i - 1
!
  p_i = p(radi)%p
  p_r = p(radr)%p
  p_s = p(rads)%p
!
  phat_rs = phat(rad_rs)%p
  phat_i  = phat(rad_i)%p
!
  p_irs = ptilde(emit)%p
!
  Q = p(1)%p + p(2)%p
  Q2 = Q**2
!
  srs = 2*p_r*p_s
  zrs = p_r*Q/((p_r + p_s)*Q)
  zsr = p_s*Q/((p_r + p_s)*Q)
!
!  print *,"srs: ",srs
!  print *,"zrs,zsr: ",zrs,zsr
!
  sitrst = 2*phat_i*phat_rs 
  zitrst = phat_i*Q/((phat_i + phat_rs)*Q)
  zrstit = phat_rs*Q/((phat_i + phat_rs)*Q)
!
!  print *,"sitrst: ",sitrst
!  print *,"zitrst,zrstit: ",zitrst,zrstit
!
  srQ = 2*p_r*Q
  ssQ = 2*p_s*Q
  prsQ = (p_r + p_s)*Q
  alprs = (sqrt(prsQ**2 - srs*Q2) - prsQ)/Q2
!
  alprs = (sqrt(prsQ**2 - srs*Q2) - prsQ)/Q2
  zrsv = (zrs + srs/(alprs*(srQ+ssQ)))
  zsrv = (zsr + srs/(alprs*(srQ+ssQ)))
!
  Ktrs = zrsv*p_s - zsrv*p_r &
       + (zsr - zrs)*srs/(-alprs*(2*phat_rs*Q))*phat_rs
!
!  print *,"Ktrs: "
!  call PrintMom(Ktrs)
!
  sitQ = 2*phat_i*Q
  srstQ = 2*phat_rs*Q
!
  pitrstQ = phat_i*Q + phat_rs*Q
!
  alpitrst = (sqrt(pitrstQ**2 - sitrst*Q2) - pitrstQ)/Q2
!
  zitrstv = (zitrst + sitrst/(alpitrst*(sitQ+srstQ)))
  zrstitv = (zrstit + sitrst/(alpitrst*(sitQ+srstQ)))
!
  Ktitrst = zitrstv*phat_rs - zrstitv*phat_i &
          + (zrstit - zitrst)*sitrst/(-alpitrst*(2*p_irs*Q))*p_irs
!
  Ktrs2 = Ktrs*Ktrs
  Ktitrst2 = Ktitrst*Ktitrst
  skk = 2*Ktrs*Ktitrst
  sitirst = 2*phat_i*p_irs
  srstirst = 2*phat_rs*p_irs
  skp = 2*Ktrs*p_irs
!
!  print *,"Ktrs2: ",Ktrs2
!  print *,"Ktitrst2: ",Ktitrst2
!
  na = (zitrstv*sitrst + (zitrst - zrstit)*sitirst) &
     / (zitrstv*sitrst*sitirst)
  nb = (zitrstv*sitrst - (zitrst - zrstit)*sitirst) &
     / (zitrstv*sitrst*srstirst)
  nc = - na*nb*sitrst/2d0
!
  nirs = na*phat_i + nb*phat_rs + nc*p_irs
!
  skn = 2*Ktrs*nirs
  spn = 2*p_irs*nirs
!
  Ktrs_arr = Ktrs
  Ktitrst_arr = Ktitrst
  p_irs_arr = p_irs
!
  call LtoH(Ktrs_arr,p_irs_arr,kPP,kPM,kMP,kMM)
  call LtoH(Ktitrst_arr,p_irs_arr,khPP,khPM,khMP,khMM)
!
!  print *,"zrs: ",zrs
!  print *,"zsr: ",zsr
!  print *,"2pi.ktrs: ",2*phat_i*Ktrs
!  print *,"ktrs2: ",ktrs2
!  print *,"sitrst: ",sitrst
!
!  print *,"ktitrst.p_irs: ",ktitrst*p_irs
!
  gPP = (-1d0,0d0)
  gMM = gPP
  gPM = (0d0,0d0)
  gMP = conjg(gPM)
!
!  print *,"Ktitrst: "
!  call PrintMom(Ktitrst)
!
!  kPP = 0 ; kMM = 0 ; kPM = 0 ; kMP = 0
!  khPP = 0 ; khMM = 0 ; khPM = 0 ; khMP = 0
!  gPP = 0 ; gMM = 0 ; gPM = 0 ; gMP = 0
!
  VPP = 4d0*qcd_ca**2 &
      * (-gPP*(zrs/zsr+zsr/zrs)*(zitrst/zrstit+zrstit/zitrst) &
      -gPP*(2*phat_i*Ktrs)**2/(2d0*Ktrs2*sitrst)*(-zrs*zsr) &
      -2d0*zrs*zsr*zrstit/zitrst*kPP &
      -2d0*zitrst*zrstit*(zrs/zsr+zsr/zrs+zrs*zsr)*khPP)
!
  VMM = VPP
!
  VPM = 4d0*qcd_ca**2 &
      * (-gPM*(zrs/zsr+zsr/zrs)*(zitrst/zrstit+zrstit/zitrst) &
      -gPM*(2*phat_i*Ktrs)**2/(2d0*Ktrs2*sitrst)*(-zrs*zsr) &
      -2d0*zrs*zsr*zrstit/zitrst*kPM &
      -2d0*zitrst*zrstit*(zrs/zsr+zsr/zrs+zrs*zsr)*khPM)
!
  VMP = conjg(VPM)
!
!  VPP = khPP ; VMM = khMM ; VPM = khPM ; VMP = khMP
!  VPP = kPP ; VMM = kMM ; VPM = kPM ; VMP = kMP
!
!  print *,"emit: ",emit
!  print *,"ptilde: "
!  call PrintParts(ptilde)
!
  call CalcBmunu_hel(emit,ptilde,M2IJ)
!
  CrsCirs = M2IJ(1,1)*VPP + M2IJ(2,2)*VMM &
          + M2IJ(1,2)*VMP + M2IJ(2,1)*VPM
!
!  CrsCirs = CrsCirs * (8*pi)**2/(srs*sitrst)/6d0
!
!  print *,"CrsCirsSG:   ",CrsCirs
!
end subroutine CalcCggCgggSG
