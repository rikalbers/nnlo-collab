! This source contains routines to calculate the strong coupling
! at various loop orders:
function calc_alphas(order,q2,lam,nf) result(as)
use math
implicit none
!
  real(kind(1d0)) :: as
!
  integer , intent(in) :: order,nf
  real(kind(1d0)) , intent(in) :: q2,lam
!
  real(kind(1d0)) , external :: alphas1loop, &
                                alphas2loop, &
                                alphas3loop
!
  if (order.eq.1) then
    as = alphas1loop(q2,lam,nf)
  elseif (order.eq.2) then
    as = alphas2loop(q2,lam,nf)
  elseif (order.eq.3) then
    as = alphas3loop(q2,lam,nf)
  else
    print *,"Problem in alphas..."
    print *,"order is not implemented yet..."
    print *,"order: ",order
    stop
  end if
!
end function calc_alphas
!
! This function calculates alphaS with one-loop accuracy:
function alphas1loop(q2,lam,nf) result(as)
use math
implicit none
!
  real(kind(1d0)) :: as
!
  real(kind(1d0)) , intent(in) :: q2,lam
  integer , intent(in) :: nf
!
  real(kind(1d0)) :: b0
!
!
  b0 = (33d0 - 2d0*nf)/(12d0*pi)
  as = 1d0 / (b0 * log(q2/lam**2))
  return
!
end function alphas1loop
!
function alphas2loop(q2,lam,nf) result(as)
use math
implicit none
!
  real(kind(1d0)) :: as
!
  real(kind(1d0)) , intent(in) :: q2,lam
  integer , intent(in) :: nf
!
  real(kind(1d0)) , save :: olam = 0d0
  real(kind(1d0)) , parameter :: xmc = 1.42d0
  real(kind(1d0)) , parameter :: xmb = 4.7d0
  real(kind(1d0)) , save :: b5,bp5,b4,bp4,b3,bp3,xlc,xlb,xllc,xllb,c45,c35
  real(kind(1d0)) :: q
  real(kind(1d0)) :: xlq,xllq
!
!
! If there was a change in Lambda_{QCD} we recalculate: 
  if (olam.ne.lam) then
    olam = lam
    b5  = (33-2*5)/pi/12
    bp5 = (153 - 19*5) / pi / 2 / (33 - 2*5)
    b4  = (33-2*4)/pi/12
    bp4 = (153 - 19*4) / pi / 2 / (33 - 2*4)
    b3  = (33-2*3)/pi/12
    bp3 = (153 - 19*3) / pi / 2 / (33 - 2*3)
    xlc = 2 * log(xmc/lam)
    xlb = 2 * log(xmb/lam)
    xllc = log(xlc)
    xllb = log(xlb)
    c45  =  1/( 1/(b5 * xlb) - xllb*bp5/(b5 * xlb)**2 ) &
         - 1/( 1/(b4 * xlb) - xllb*bp4/(b4 * xlb)**2 )
    c35  =  1/( 1/(b4 * xlc) - xllc*bp4/(b4 * xlc)**2 ) &
         - 1/( 1/(b3 * xlc) - xllc*bp3/(b3 * xlc)**2 ) + c45
  end if
  q = sqrt(q2)
  xlq = 2 * log(q/lam)
  xllq = log(xlq)
  if (nf.eq.5) then
    as = 1/(b5 * xlq) -  bp5/(b5 * xlq)**2 * xllq
  elseif(nf.eq.4) then
    as = 1/( 1/(1/(b4 * xlq) - bp4/(b4 * xlq)**2 * xllq) + c45 )
  elseif(nf.eq.3) then
    as = 1/( 1/(1/(b3 * xlq) - bp3/(b3 * xlq)**2 * xllq) + c35 )
  else
    print *,"2-loop alphaS cannot be calculated for nf = ",nf
    stop
  endif
!
end function alphas2loop
!
function alphas3loop(q2,lam,nf) result(as)
use qcdparams
use math
implicit none
!
  real(kind(1d0)) :: as
!
  real(kind(1d0)) , intent(in) :: q2,lam
  integer , intent(in) :: nf
!
  real(kind(1d0)) :: q,xlq,xllq
  real(kind(1d0)) :: b0,b1,b2
!
!
! As it was taken from arXiv:1402.4140:
  b0 = (11*qcd_ca - 4*qcd_tr*nf)/6d0
  b1 = (17*qcd_ca**2 - 10*qcd_ca*qcd_tr*nf &
     - 6*qcd_cf*qcd_tr*nf)/6d0
  b2 = (2857*qcd_ca**3 + 108*qcd_cf**2*qcd_tr*nf &
     - 1230*qcd_cf*qcd_ca*qcd_tr*nf &
     - 2830*qcd_ca**2*qcd_tr*nf &
     + 264*qcd_cf*qcd_tr**2*nf**2 &
     + 316*qcd_ca*qcd_tr**2*nf**2)/432d0
!
  q = sqrt(q2)
  xlq = 2 * log(q/lam)
  xllq = log(xlq)
!
  as = 2*pi/(b0*xlq)*( &
       1 - b1/b0**2*xllq/xlq + 1/(b0*xlq)**2*( &
         b1**2/b0**2*(xllq**2 - xllq - 1) + b2/b0))
!
end function alphas3loop
!
! This function is taken from POWHEG-BOX:
function genericxlambdL(as,q,nf)
use math
implicit none
!
  real(kind(1d0)) :: genericxlambdL
!  
  real(kind(1d0)) , intent(in) :: as,q
  integer , intent(in) :: nf
!  
  real(kind(1d0)) ::  b,t,xlt,ot,as0,as1
!
!
  b  = (33-2*nf)/pi/12
  t  = 1/b/as
  do while (.true.)
    xlt = log(t)
    ot = t
!-----------------------------------------------------------
! Value and Derivative of alfa with respect to t
    as0  = 1/b/t
    as1  = - 1/b/t**2
    t  = (as-as0)/as1 + t
    if(abs(ot-t)/ot.lt.0.00000001d0) exit
  end do
  genericxlambdL = q/exp(t/2)
  return
end function genericxlambdL
!
function genericxlambdNL(as,q,nf)
use math
implicit none
!
  real(kind(1d0)) :: genericxlambdNL
!
  real(kind(1d0)) , intent(in) :: as,q
  integer , intent(in) :: nf
!
  real(kind(1d0)) :: b,bp,t,xlt,ot,as0,as1
!
!
  b  = (33-2*nf)/pi/12
  bp = (153 - 19*nf) / pi / 2 / (33 - 2*nf)
  t  = 1/b/as
  do while (.true.)
    xlt = log(t)
    ot = t
!-----------------------------------------------------------
! Value and Derivative of alfa with respect to t
    as0  = 1/b/t - bp*xlt/(b*t)**2
    as1  = - 1/b/t**2 -bp/b**2*(1-2*xlt)/t**3
    t  = (as-as0)/as1 + t
    if(abs(ot-t)/ot.lt.0.00000001d0) exit
  end do
  genericxlambdNL = q/exp(t/2)
  return
end function genericxlambdNL
!
function genericxlambdNNL(as,q,nf)
use math
implicit none
!
  real(kind(1d0)) :: genericxlambdNNL
!
  real(kind(1d0)) , intent(in) :: as,q
  integer , intent(in) :: nf
!
  real(kind(1d0)) :: b0,b1,b2,t,xlt,ot,as0,as1
  integer :: icount
!
!
  b0 = (33.d0-2.d0*nf)/(12.d0*pi)
  b1 = (153.d0 - 19.d0*nf) / (24.d0*pi**2)
  b2 = (2857.d0/2.d0-5033.d0/18.d0*nf+325.d0/54.d0*nf**2) &
     / (4.d0*pi)**3
  t = 1/b0/as
  icount = 0
  do while (.true.)
    xlt = log(t)
    if (icount.gt.10000) then
      write(*,*) ' xlambd: cannot converge '
      stop
    end if
    icount = icount + 1
    ot = t
!-----------------------------------------------------------
! Value and Derivative of alfa with respect to t
    as0 = 1/(t*b0)*(1-b1/b0**2*log(t)/t        &
        + (b1/b0**2*log(t)/t)**2               &
        - (b1**2*(log(t)+1)-b0*b2)/b0**4/t**2)
    as1 = (-2*b1**2*log(t)**2/(b0**4*t**3)+2*(b1**2*(log(t)+1)-b0*b2)     &
        / (b0**4*t**3)+b1*log(t)/(b0**2*t**2)+2*b1**2*log(t)/(b0**4*t**3) &
        - b1/(b0**2*t**2)-b1**2/(b0**4*t**3))/(b0*t)-(b1**2*log(t)**2     &
        / (b0**4*t**2)-(b1**2*(log(t)+1)-b0*b2)/(b0**4*t**2)-b1*log(t)    &
        / (b0**2*t)+1)/(b0*t**2)
    t = (as-as0)/as1 + t
    if (abs(ot-t)/ot.lt.0.00000001d0) exit
  end do
  genericxlambdNNL = q/exp(t/2)
  return
end function genericxlambdNNL
