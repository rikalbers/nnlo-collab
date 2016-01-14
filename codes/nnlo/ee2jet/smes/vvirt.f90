    subroutine VVirtSME(iptrn,p,mur,VVLaurent,VLaurent,BLaurent)
use process
use particles
use QCDparams
use my_model
use misc
implicit none
!
  integer , intent(in) :: iptrn
  type(particle) , dimension(:) , intent(in) :: p
  real(kind(1d0)) , intent(in) :: mur
  real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: &
    VVLaurent,VLaurent,BLaurent
!
  real(kind(1d0)) ::  NC,nf,TR,nfg,pi,z3,eg,K,pg21,c1
  common/piez/ pi,z3
  common/const/ NC,nf,TR,eg
!
  integer :: i,j
  real(kind(1d0)) :: Q2,y12,y15,y25
  real(kind(1d0)) :: L
  real(kind(1d0)) , dimension(-4:2) :: LaurentMuComp,Tqq61x1,Tqq62x0,Tqqsum
  real(kind(1d0)) , dimension(-4:4) :: rf
  complex(kind(1d0)) , dimension(-4:2) :: f1yz,f1zy,f2yz,f2zy
  type(mom) :: Q
!
  pi = 3.14159265358979324d0
  z3 = 1.20205690315959429d0 !zeta(3)
  eg = 0.57721566490153286d0 !eulergamma
  K = (67d0/18d0 - pi**2/6d0)*qcd_ca - 10d0/9d0*qcd_tr*qcd_nf
  pg21 = -2d0*z3 !polygamma(2,1)
  c1 = 43d0/4d0 - pi**2 / 3d0
!
! Calculation of invariants:
  Q = p(1)%p + p(2)%p
  Q2 = Q**2
!
  y12 = 2*p(3)%p*p(4)%p/Q2
  y15 = 2*p(3)%p*p(2)%p/Q2
  y25 = 2*p(4)%p*p(2)%p/Q2
!
  BLaurent = 0
  BLaurent(0) = 2d0*qcd_nc*(y15**2 + y25**2)/y12**2
  BLaurent(1) = -2d0*qcd_nc
  BLaurent = BLaurent/ 9d0
!
  VLaurent = 0
  VLaurent(-2) = -1d0
  VLaurent(-1) = -3d0/2d0
  VLaurent(0)  = -4d0+ 7d0*pi**2 / 12d0
  VLaurent(1)  = -8d0+ 7d0*pi**2 / 8d0 + 7d0*z3/3d0
  VLaurent(2)  = -16d0 +7d0*pi**2 / 3d0 +7d0*z3/2d0 - 73d0*pi**4 / 1440d0
  VLaurent = VLaurent * (qcd_nc- 1d0/qcd_nc)
  VLaurent = SeriesProd(VLaurent,BLaurent)
!
! Factor to change the normalization (Exp[-ep gamma] to Gamma[1-ep])
!  rf = 0d0
!  rf(0) = 1d0
!  rf(2) = pi**2 /6d0
!  rf(3) = -pg21/3d0
!  rf(4) = 7d0*pi**4 /360d0
! Including the Born expansion in the overal factor
!  do i=0,4
!    rf(4-i) = rf(4-i)*BLaurent(0) + rf(3-i)*BLaurent(1)
!  end do
! Turn off the r[e]^2 term
  rf = 0d0
  rf(0) = BLaurent(0)
  rf(1) = BLaurent(1)
!
  Tqq61x1 = 0d0
  Tqq61x1(-4) = 1d0/4d0
  Tqq61x1(-3) = 3d0/4d0
  Tqq61x1(-2) = (123d0- 2d0*pi**2)/48d0
  Tqq61x1(-1) = (168d0- 3d0*pi**2- 28d0*z3)/24d0
  Tqq61x1(0)  = (8640d0- 205d0*pi**2- 7d0*pi**4- 1680d0*z3)/480d0
  Tqq61x1 = Tqq61x1*(qcd_nc**2- 1d0)**2 /(qcd_nc**2)
!  Tqq61x1 = SeriesProd(BLaurent,Tqq61x1)
!
  Tqq62x0 = 0d0
  Tqq62x0(-4) = (qcd_nc**2- 1d0)/4d0
  Tqq62x0(-3) = (17d0*qcd_nc**2- 2d0*qcd_nc*qcd_nf- 6d0)/8d0
  Tqq62x0(-2) = (-369d0+ 433d0*qcd_nc**2- 16d0*qcd_nc*qcd_nf+ 78d0*pi**2      &
                -72d0*qcd_nc**2 * pi**2)/144d0
  Tqq62x0(-1) = (-5967d0+ 4045d0*qcd_nc**2+ 260d0*qcd_nc*qcd_nf+ 1296d0*pi**2 &
                -1494d0*qcd_nc**2 *pi**2+ 36d0*qcd_nc*qcd_nf*pi**2+ 2304d0*z3 &
                +504d0*qcd_nc**2 *z3)/864d0
  Tqq62x0(0)  = (-466155d0- 45415d0*qcd_nc**2+ 81700d0*qcd_nc*qcd_nf          &
                +128250d0*pi**2- 64590d0*qcd_nc**2 *pi**2                     &
                - 10920d0*qcd_nc*qcd_nf*pi**2- 5310d0*pi**4                   &
                + 4734d0*qcd_nc**2 *pi**4+ 187920d0*z3+ 37440d0*qcd_nc**2 *z3 &
                + 1440d0*qcd_nc*qcd_nf*z3)/25920d0
  Tqq62x0 = Tqq62x0*(qcd_nc**2- 1d0)/(qcd_nc**2)
!  Tqq62x0 = SeriesProd(BLaurent,Tqq62x0)
!
! Compute VVLaurent and fix the normalization
  Tqqsum = Tqq62x0 + Tqq61x1
!
  VVLaurent = 0
  VVLaurent(-4) = Tqqsum(-4)*rf(0)
  VVLaurent(-3) = Tqqsum(-4)*rf(1)+ Tqqsum(-3)*rf(0)
  VVLaurent(-2) = Tqqsum(-4)*rf(2)+ Tqqsum(-3)*rf(1)      &
                + Tqqsum(-2)*rf(0)
  VVLaurent(-1) = Tqqsum(-4)*rf(3)+ Tqqsum(-3)*rf(2)      &
                + Tqqsum(-2)*rf(1)+ Tqqsum(-1)*rf(0)
  VVLaurent(0)  = Tqqsum(-4)*rf(4)+Tqqsum(-3)*rf(3)       &
                +Tqqsum(-2)*rf(2)                         &
                +Tqqsum(-1)*rf(1)+Tqqsum(0)*rf(0)
!
!  if (.true.) then
!    do i=-4,2
!    print *,i, ": ",VVLaurent(i)
!   end do
!  print *,""
!  end if
!  stop
!
! Pattern 1 is for d d~ g and pattern 2 is for u u~ g:
  if (iptrn.eq.2) then
    VVLaurent = 4*VVLaurent
    VLaurent  = 4*VLaurent
    BLaurent  = 4*BLaurent
  end if
!
! Normalization of 1/(2pi) with respect to the born matrix element
  VLaurent = VLaurent/(2d0*pi)
  VVLaurent = VVLaurent/(2d0*pi)**2
!
! Factoring out gs instead of alpha_s to get agreement with
! the normalization in contvv.f90
  VLaurent = VLaurent/(4d0*pi)
  VVLaurent = VVLaurent/(4d0*pi)**2
!
end subroutine VVirtSME
