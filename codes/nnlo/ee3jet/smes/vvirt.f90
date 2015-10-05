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
  real(kind(1d0)) ::  NC,nf,TR,nfg,pi,z3,eg,K
  common/piez/ pi,z3
  common/const/ NC,nf,TR,eg
!
  integer :: i,j
  real(kind(1d0)) :: Q2,y12,y13,y23
  real(kind(1d0)) :: H2qqgm1,polrem
  real(kind(1d0)) :: L
  real(kind(1d0)) , dimension(-4:2) :: I1qqg,I1qqg2,I1qqgsq,pol1,pol2
  real(kind(1d0)) , dimension(-4:2) :: LaurentMuComp
  complex(kind(1d0)) , dimension(-4:2) :: f1yz,f1zy,f2yz,f2zy
  type(mom) :: Q
!
  real(kind(1d0)) , external :: lof
  real(kind(1d0)) , external :: fin11,fin20
!
  interface
    subroutine VirtSMEddim_new(iptrn,pvirt,mur,VirtLaurent,BornLaurent)
    use particles
    implicit none
!
      integer , intent(in) :: iptrn
      type(particle) , dimension(:) , intent(in) :: pvirt
      real(kind(1d0)) , intent(in) :: mur
      real(kind(1d0)) , dimension(-4:2) , intent(out) :: VirtLaurent,BornLaurent
!
    end subroutine VirtSMEddim_new
  end interface
!
  pi = 3.14159265358979324d0
  z3 = 1.20205690315959429d0 !zeta(3)
  eg = 0.57721566490153286d0 !eulergamma
  K = (67d0/18d0 - pi**2/6d0)*qcd_ca - 10d0/9d0*qcd_tr*qcd_nf
!
  if (qcd_nf.eq.6.0d0) then
    nfg = 9/15d0
  else if (qcd_nf.eq.5.0d0) then
    nfg = 1/11d0
  else if (qcd_nf.eq.4.0d0) then
    nfg = 2/5d0
  else if (qcd_nf.eq.3.0d0) then
    nfg = 0d0
  else if (qcd_nf.eq.2.0d0) then
    nfg = 1/5d0
  else if (qcd_nf.eq.1.0d0) then
    nfg = 1d0
  end if
!
! Calculation of invariants:
!
  Q = p(1)%p + p(2)%p
  Q2 = Q**2
!
  y12 = 2*p(3)%p*p(4)%p/Q2
  y13 = 2*p(3)%p*p(5)%p/Q2
  y23 = 2*p(4)%p*p(5)%p/Q2
!
  BLaurent = 0
  BLaurent(0) = 8*qcd_nc*qcd_cf*(y13/y23+y23/y13+2*y12/(y13*y23))
  BLaurent(1) = -16*qcd_nc*qcd_cf*(y13/y23+y23/y13+y12/(y13*y23)+1)
  BLaurent(2) = 8*qcd_nc*qcd_cf*(y13/y23+y23/y13+2)
!
  call f1(f1yz,y13,y23)
  call f1(f1zy,y23,y13)
  call f2(f2yz,y13,y23)
  call f2(f2zy,y23,y13)
  do i=-4,2
    VLaurent(i) = 4*qcd_nc*qcd_cf*real(qcd_nc*(f1yz(i)+f1zy(i)) &
                + 1d0/qcd_nc*(f2yz(i)+f2zy(i)))
  end do
  do i=-4,1
    VLaurent(i) = VLaurent(i) - qcd_beta0/2d0*BLaurent(i+1)
  end do
!
  I1qqg(-4) = 0
  I1qqg(-3) = 0
  I1qqg(-2) = -qcd_nc + 1/(2*qcd_nc)

  I1qqg(-1) = qcd_nc/6d0*(-10 + 3*lof(y13,1) + 3*lof(y23,1))   &
            - 1d0/qcd_nc*(lof(y12,1)/2d0 - 3d0/4d0) + qcd_nf/6d0

  I1qqg(0) = qcd_nc/12d0*(pi**2 + 10*(lof(y13,1) + lof(y23,1))     &
             - 3*(lof(y13,2) + lof(y23,2))) - 1/qcd_nc*(pi**2/24d0 &
             + 3*lof(y12,1)/4d0 - lof(y12,2)/4d0)                  &
           - qcd_nf/12d0*(lof(y13,1) + lof(y23,1))

  I1qqg(1) = qcd_nc/72d0*(10*pi**2 + 24*z3 - 3*pi**2*(lof(y13,1) + lof(y23,1))  &
           - 30*(lof(y13,2) + lof(y23,2)) + 6*(lof(y13,3) + lof(y23,3)))        &
           - 1/(48*qcd_nc)*(3*pi**2 + 8*z3 - 2*pi**2*lof(y12,1) - 18*lof(y12,2) &
           + 4*lof(y12,3)) + qcd_nf/72d0*(3*(lof(y13,2) + lof(y23,2)) - pi**2)

  I1qqg(2) = qcd_nc/1440d0*(-pi**4 + 800*z3 + 10*(                              &
           - (10*pi**2 + 24*z3)*(lof(y13,1) + lof(y23,1))                       &
           + 3*pi**2*(lof(y13,2) + lof(y23,2)) + 20*(lof(y13,3) + lof(y23,3))   &
           - 3*(lof(y13,4) + lof(y23,4))))                                      &
           - 1/(2880*qcd_nc)*(-pi**4 + 720*z3 - 60*((3*pi**2 + 8*z3)*lof(y12,1) &
           - pi**2*lof(y12,2) - 6*lof(y12,3) + lof(y12,4)))                     &
           + qcd_nf/144d0*(-8*z3 + pi**2*(lof(y13,1) + lof(y23,1))              &
           - 2*(lof(y13,3) + lof(y23,3)))

  do i=-4,2
    I1qqg2(i) = I1qqg(i)*(2d0**i)
  end do
!
  I1qqgsq = SeriesProd(I1qqg,I1qqg)
!
  H2qqgm1 = ((4*z3 + 589d0/432d0 - 11*pi**2/72d0)*qcd_nc**2                            &       
          - z3/2d0 - 41d0/54d0 - pi**2/48d0 + (pi**2/4d0 - 3*z3 - 3/16d0)/qcd_nc**2    &
          + (pi**2/36d0 - 19/18d0)*qcd_nc*qcd_nf - (1/54d0 + pi**2/24d0)*qcd_nf/qcd_nc &
          + 5d0/27d0*qcd_nf**2 )/4d0
!
! poles for 2Re<M0|M2> + |M1|^2
  pol1 = 0
  pol2 = 0
! 2x0 part
  do i = -4,0
    pol2(i) = 2*(-I1qqgsq(i) + qcd_beta0/2d0*(I1qqg2(i+1) - I1qqg(i+1)) &
            + K*I1qqg2(i))
    if (i.gt.-4) then
      pol2(i) = pol2(i) + 2*(qcd_beta0/2d0*pi**2/4d0*I1qqg2(i-1))
      if (i.gt.-3) then
        pol2(i) = pol2(i) + 2*(K*pi**2/4d0 + 7*qcd_beta0/2d0/3d0*z3)*I1qqg2(i-2)
      end if
    end if
  end do
  pol2(-1) = pol2(-1) + 2*H2qqgm1
!
  do i=0,2
    pol2(-i) = pol2(-i)*BLaurent(0) + pol2(-i-1)*BLaurent(1) + &
    pol2(-i-2)*BLaurent(2)
  end do
  pol2(-3) = pol2(-3)*BLaurent(0) + pol2(-4)*BLaurent(1)
  pol2(-4) = pol2(-4)*BLaurent(0)

! 1x1 part
  do i=-4,0
    do j=0,2
      pol1(i) = pol1(i) + 2*I1qqg(-j)*VLaurent(i+j)
    end do
    if (i.gt.-4) then
      pol1(i) = pol1(i) + 2*I1qqg(1)*VLaurent(i-1)
      if (i.gt.-3) then
        pol1(i) = pol1(i) + 2*I1qqg(2)*VLaurent(i-2)
      end if
    end if
  end do
!
  do i=-4,1
    VVLaurent(i) = pol2(i) + pol1(i)
  end do
!
  polrem = pol2(0) + pol1(0)
  VVLaurent(0) = polrem + (qcd_nc**2 - 1)*(                          &
      qcd_nc**2*(fin20(1,y13,y23) + fin20(1,y23,y13))                &
    + fin20(2,y13,y23) + fin20(2,y23,y13)                            &
    + (fin20(3,y13,y23) + fin20(3,y23,y13))/(qcd_nc**2)              &
    + qcd_nc*qcd_nf*(fin20(4,y13,y23) + fin20(4,y23,y13))            &
    + qcd_nf/qcd_nc*(fin20(5,y13,y23) + fin20(5,y23,y13))            &
    + qcd_nf**2*(fin20(6,y13,y23) + fin20(6,y23,y13))                &
    + nfg*(4/qcd_nc - qcd_nc)*(fin20(7,y13,y23) + fin20(7,y23,y13))) &
    + (qcd_nc**2 - 1)*(                                              &
      qcd_nc**2*(fin11(1,y13,y23) + fin11(1,y23,y13))                &
    + fin11(2,y13,y23) + fin11(2,y23,y13)                            &
    + (fin11(3,y13,y23) + fin11(3,y23,y13))/(qcd_nc**2)              &
    + qcd_nc*qcd_nf*(fin11(4,y13,y23) + fin11(4,y23,y13))            &
    + qcd_nf/qcd_nc*(fin11(5,y13,y23) + fin11(5,y23,y13))            &
    + qcd_nf**2*(fin11(6,y13,y23) + fin11(6,y23,y13)))
!
! Calculation of the mu dependent part:
!
! xi = mur**2/Q**2:
  L = log(mur**2/Q2)
!
!
  LaurentMuComp = 0
  LaurentMuComp(-4) = 0
  LaurentMuComp(-3) = 2*VVLaurent(-4)*L
  LaurentMuComp(-2) = 2*VVLaurent(-4)*L**2       &
                    + 2*VVLaurent(-3)*L          &
                    + 1*qcd_beta0*VLaurent(-2)*L
  LaurentMuComp(-1) = 4d0/3d0*VVLaurent(-4)*L**3    &
                    + 2*VVLaurent(-3)*L**2          &
                    + 3d0/2d0*qcd_beta0*VLaurent(-2)*L**2 &
                    + 2*VVLaurent(-2)*L             &
                    + 1*qcd_beta0*VLaurent(-1)*L
  LaurentMuComp( 0) = 2d0/3d0*VVLaurent(-4)*L**4          &
                    + 4d0/3d0*VVLaurent(-3)*L**3          & 
                    + 7d0/6d0*qcd_beta0*VLaurent(-2)*L**3 &
                    + 2*VVLaurent(-2)*L**2                &
                    + 3d0/2d0*qcd_beta0*VLaurent(-1)*L**2       &
                    + 0.25d0*qcd_beta0**2*BLaurent( 0)*L**2      &
                    + 2*VVLaurent(-1)*L                   &
                    + 1*qcd_beta0*VLaurent( 0)*L          &
                    + qcd_beta1/2d0*BLaurent( 0)*L 
!
  VVLaurent = LaurentMuComp + VVLaurent
!
  LaurentMuComp = 0
!
  LaurentMuComp(-2) = 0
  LaurentMuComp(-1) = VLaurent(-2)*L
  LaurentMuComp( 0) = 0.5d0*VLaurent(-2)*L**2 &
                    + (VLaurent(-1) &
                    + 0.5d0*qcd_beta0*BLaurent(0))*L
  LaurentMuComp( 1) = 1d0/6d0*VLaurent(-2)*L**3 &
                    + 0.5d0*(VLaurent(-1) &
                    + 0.5d0*qcd_beta0*BLaurent(0))*L**2 &
                    + (VLaurent( 0) &
                    + 0.5d0*qcd_beta0*BLaurent(1))*L
  LaurentMuComp( 2) = 1d0/24d0*VLaurent(-2)*L**4 &
                    + 1d0/6d0*(VLaurent(-1) &
                    + 0.5d0*qcd_beta0*BLaurent(0))*L**3 &
                    + 0.5d0*(VLaurent( 0) &
                    + 0.5d0*qcd_beta0*BLaurent(1))*L**2 &
                    + (VLaurent( 1) &
                    + 0.5d0*qcd_beta0*BLaurent(2))*L
!
  VLaurent = VLaurent + LaurentMuComp
!
! Pattern 1 is for d d~ g and pattern 2 is for u u~ g:
  if (iptrn.eq.2) then
    VVLaurent = 4*VVLaurent
    VLaurent  = 4*VLaurent
    BLaurent  = 4*BLaurent
  end if
  VVLaurent = VVLaurent/Q2/27d0/(2*pi)**2
  VLaurent  = VLaurent/Q2/27d0/(2*pi)**2
  BLaurent  = BLaurent/Q2/27d0/(2*pi)**2
!
end subroutine VVirtSME
