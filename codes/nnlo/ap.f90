! This source contains the Altarelli-Parisi splitting
! kernels at various orders:
function Pqg0(zi,zr) result(ans)
use QCDparams
implicit none
!
  real(kind(1d0)) :: ans
!
  real(kind(1d0)) , intent(in) :: zi,zr
!
!
  ans = qcd_cf * (1 + zi**2)/(1 - zi)
  return
!
end function Pqg0
!
function Pgg0(zi,zr,kt,M2munu) result(ans)
use QCDparams
use momenta
use particles
use misc
implicit none
!
  real(kind(1d0)) :: ans
!
  real(kind(1d0)) , intent(in) :: zi,zr
  type(particle) , intent(in) :: kt
  real(kind(1d0)) , dimension(0:3,0:3) , intent(in) :: M2munu
!
  real(kind(1d0)) :: kt2
!
!
  kt2 = kt%p*kt%p
!
  ans = 2*qcd_ca*((zi/zr + zr/zi) &
      * (-M2munu(0,0) + M2munu(1,1) + M2munu(2,2) + M2munu(3,3)) &
      - 2*zi*zr*(kt%p*M2munu*kt%p)/kt2)
!
end function Pgg0
!
function Pqq0(zi,zr,kt,M2munu) result(ans)
use QCDparams
use momenta
use particles
use misc
implicit none
!
  real(kind(1d0)) :: ans
!
  real(kind(1d0)) , intent(in) :: zi,zr
  type(particle) , intent(in) :: kt
  real(kind(1d0)) , dimension(0:3,0:3) , intent(in) :: M2munu
!
  real(kind(1d0)) :: kt2
!
!
  kt2 = kt%p*kt%p
!
  ans = qcd_tr*( &
        (-M2munu(0,0) + M2munu(1,1) + M2munu(2,2) + M2munu(3,3)) &
      + 4*zi*zr*(kt%p*M2munu*kt%p)/kt2)
!
end function Pqq0
!
! This is the spin-averaged version of the q-pair splitting
! function:
function Pqq0sa(zi,zr) result(ans)
use QCDparams
implicit none
!
  real(kind(1d0)) :: ans
!
  real(kind(1d0)) , intent(in) :: zi,zr
!
!
  ans = qcd_tr*(1 - 2*zi*zr)
!
end function Pqq0sa
!
! This is the spin-averaged version of the g-g splitting
! function:
function Pgg0sa(zi,zr) result(ans)
use QCDparams
implicit none
!
  real(kind(1d0)) :: ans
!
  real(kind(1d0)) , intent(in) :: zi,zr
!
!
  ans = 2*qcd_ca*(zi/zr + zr/zi + zi*zr)
!
end function Pgg0sa
! 
! This is <\hat{P}_{\bar{q}_1'q_2'q_3}> as presented in Eq. (57) in
! arXiv:hep-ph/9908523:
function Prbrq0(s123,s12,s13,s23,z1,z2,z3) result(ans)
use QCDparams
implicit none
!
  real(kind(1d0)) :: ans
!
  real(kind(1d0)) , intent(in) :: s123
  real(kind(1d0)) , intent(in) :: s12,s13,s23
  real(kind(1d0)) , intent(in) :: z1,z2,z3
!
  real(kind(1d0)) :: t123
!
!
  ans = 0d0
!
! Definition is in Eq. (22) of arXiv:hep-ph/9908523:
  t123 = 2*(z1*s23 - z2*s13)/(z1 + z2) + (z1 - z2)/(z1 + z2)*s12
!
  ans = -t123**2/(s12*s123) + (4*z3 + (z1 - z2)**2)/(z1 + z2) &
      + z1 + z2 - s12/s123
  ans = 0.5d0*qcd_cf*qcd_tr*s123/s12*ans
!
end function Prbrq0
!
! This is <\hat{P}_{\bar{q}_1 q_2 q_3}> as presented in Eq. (58) in
! arXiv:hep-ph/9908523:
function Pqbqq0(s123,s12,s13,s23,z1,z2,z3) result(ans)
use QCDparams
implicit none
!
  real(kind(1d0)) :: ans
!
  real(kind(1d0)) , intent(in) :: s123
  real(kind(1d0)) , intent(in) :: s12,s13,s23
  real(kind(1d0)) , intent(in) :: z1,z2,z3
!
!
  real(kind(1d0)) , external :: Prbrq0
!
  ans = 0d0
!
  ans = Prbrq0(s123,s12,s13,s23,z1,z2,z3) &
      + Prbrq0(s123,s13,s12,s23,z1,z3,z2)
  ans = ans &
      + qcd_cf*(qcd_cf - 0.5d0*qcd_ca) &
      * ( &
          2*s23/s12 &
      +   s123/s12*((1 + z1**2)/(1 - z2) - 2*z2/(1 - z3)) &
      -   s123**2*z1*(1 + z1**2)/(2*s12*s13*(1 - z2)*(1 - z3)) &
! The same as the previous 3 lines but with 2 <-> 3:
      +   2*s23/s13 &
      +   s123/s13*((1 + z1**2)/(1 - z3) - 2*z3/(1 - z2)) &
      -   s123**2*z1*(1 + z1**2)/(2*s13*s12*(1 - z3)*(1 - z2)) &
        )
!
end function Pqbqq0
!
! This is <\hat{P}^{(ab)}_{g_1g_2q_3}> as presented in Eq. (61) in
! arXiv:hep-ph/9908523:
function Pabggq0(s123,s12,s13,s23,z1,z2,z3) result(ans)
implicit none
!
  real(kind(1d0)) :: ans
!
  real(kind(1d0)) , intent(in) :: s123
  real(kind(1d0)) , intent(in) :: s12,s13,s23
  real(kind(1d0)) , intent(in) :: z1,z2,z3
!
!
!
  ans = 0d0
!
  ans = ans &
      + s123**2*z3*(1 + z3**2)/(2*s13*s23*z1*z2) &
      + s123*(z3*(1 - z1) + (1 - z2)**3)/(s13*z1*z2) &
      - s23/s13 &
! The same as the previous three lines but weth the 1 <-> 2 interchange:
      + s123**2*z3*(1 + z3**2)/(2*s23*s13*z2*z1) &
      + s123*(z3*(1 - z2) + (1 - z1)**3)/(s23*z2*z1) &
      - s13/s23
!
end function Pabggq0
!
! This is <\hat{P}^{(nab)}_{g_1g_2q_3}> as presented in Eq. (61) in
! arXiv:hep-ph/9908523:
function Pnabggq0(s123,s12,s13,s23,z1,z2,z3) result(ans)
implicit none
!
  real(kind(1d0)) :: ans
!
  real(kind(1d0)) , intent(in) :: s123
  real(kind(1d0)) , intent(in) :: s12,s13,s23
  real(kind(1d0)) , intent(in) :: z1,z2,z3
!
  real(kind(1d0)) :: t123
!
!
  ans = 0d0
!
! Definition is in Eq. (22) of arXiv:hep-ph/9908523:
  t123 = 2*(z1*s23 - z2*s13)/(z1 + z2) + (z1 - z2)/(z1 + z2)*s12
!
  ans = ans &
      + t123**2/(4*s12**2) + 0.25d0 + s123**2/(2*s12*s13) &
      * (((1 - z3)**2 + 2*z3)/z2 + (z2**2 + 2*(1 - z2))/(1 - z3)) &
      - s123**2*z3*((1 - z3)**2 + 2*z3)/(4*s13*s23*z1*z2) &
      + s123*(z1*(2 - 2*z1 + z1**2) - z2*(6 - 6*z2 + z2**2)) &
      /   (2*s12*z2*(1 - z3)) &
      + s123/(2*s13)*(((1 - z2)**3 + z3**2 - z2)/(z2*(1 - z3)) &
      -   (z3*(1 - z1) + (1 - z2)**3)/(z1*z2)) &
! The same as the previous three lines but weth the 1 <-> 2 interchange:
      + t123**2/(4*s12**2) + 0.25d0 + s123**2/(2*s12*s23) &
      * (((1 - z3)**2 + 2*z3)/z1 + (z1**2 + 2*(1 - z1))/(1 - z3)) &
      - s123**2*z3*((1 - z3)**2 + 2*z3)/(4*s23*s13*z2*z1) &
      + s123*(z2*(2 - 2*z2 + z2**2) - z1*(6 - 6*z1 + z1**2)) &
      /   (2*s12*z1*(1 - z3)) &
      + s123/(2*s23)*(((1 - z1)**3 + z3**2 - z1)/(z1*(1 - z3)) &
      -   (z3*(1 - z2) + (1 - z1)**3)/(z2*z1))
!
end function Pnabggq0
!
! This is <\hat{P}_{g_1g_2q_3}> as presented in Eq. (60) in
! arXiv:hep-ph/9908523:
function Pggq0(s123,s12,s13,s23,z1,z2,z3) result(ans)
use QCDparams
implicit none
!
  real(kind(1d0)) :: ans
!
  real(kind(1d0)) , intent(in) :: s123
  real(kind(1d0)) , intent(in) :: s12,s13,s23
  real(kind(1d0)) , intent(in) :: z1,z2,z3
!
!
  real(kind(1d0)) , external :: Pabggq0,Pnabggq0
!
  ans = qcd_cf**2*Pabggq0(s123,s12,s13,s23,z1,z2,z3) &
      + qcd_cf*qcd_ca*Pnabggq0(s123,s12,s13,s23,z1,z2,z3)
!
end function Pggq0
! 
! This is \hat{P}^{\mu\nu}_{g_1q_2\bar{q}_3} as presented in
! Eq. (63) of arXiv:hep-ph/9908523: 
! Since this kernel contains spin-correlated MEs the abelian and
! non-abelian part is considered in the same function to reduce
! the overhead produced by re-evaluating the various products:
function Pgqqb0(s123,s12,s13,s23,z1,z2,z3,kt1,kt2,kt3,SME,SMEmunu) result(ans)
use QCDparams
use particles
implicit none
!
  real(kind(1d0)) :: ans
!
  real(kind(1d0)) , intent(in) :: s123
  real(kind(1d0)) , intent(in) :: s12,s13,s23
  real(kind(1d0)) , intent(in) :: z1,z2,z3
  type(particle) , intent(in) :: kt1,kt2,kt3
  real(kind(1d0)) , intent(in) :: SME
  real(kind(1d0)) , dimension(0:3,0:3) , intent(in) :: SMEmunu
!
  real(kind(1d0)) Pabgqqb0,Pnabgqqb0
  real(kind(1d0)) :: t231
  real(kind(1d0)) :: kt1kt1,kt2kt2,kt3kt3,kt23kt23
  real(kind(1d0)) :: kt1SMEkt1,kt2SMEkt2,kt3SMEkt3,kt23SMEkt23
  type(mom) :: kt23
!
!
! Definition is in Eq. (22) of arXiv:hep-ph/9908523:
  t231 = 2*(z2*s13 - z3*s12)/(z2 + z3) + (z2 - z3)/(z2 + z3)*s23
!
! Constructing k_{\bot,rs} = k_{\bot,r}/z_r - k_{\bot,s}/z_s:
  kt23 = kt2%p/z2 - kt3%p/z3
!
! dot products of the various transverse momenta:
  kt1kt1 = kt1%p*kt1%p
  kt2kt2 = kt2%p*kt2%p
  kt3kt3 = kt3%p*kt3%p
  kt23kt23 = kt23*kt23
! contractions with the spin correlated ME:
  kt1SMEkt1 = kt1%p*SMEmunu*kt1%p
  kt2SMEkt2 = kt2%p*SMEmunu*kt2%p
  kt3SMEkt3 = kt3%p*SMEmunu*kt3%p
!
  kt23SMEkt23 = kt23*SMEmunu*kt23
!
  Pabgqqb0 = &
    (-2+((s123-s23)**2+2*s123*s23)/(s12*s13))*SME+  &
    (4*s123*(-((kt1SMEkt1*(s23*z1-s123*(1-z1)*z1))/kt1kt1)+  &
    z2*z3*((kt2SMEkt2*(s13*z2-s123*(1-z2)*z2))/(kt2kt2*z2**2)+ &
    (kt23SMEkt23*s23)/(kt23kt23*z2*z3)+(kt3SMEkt3*(s12*z3- &
    s123*(1-z3)*z3))/(kt3kt3*z3**2))))/(s12*s13)
!
  Pnabgqqb0 = &
    (SME*(-1+(2*s123*(1-z1+2*z1**2))/(s23*(1-z1)*z1)+(2*s123*(1-  &
    z3))/(s12*(1-z1)*z1))+(s123*(-((SME*t231**2)/s123)+  &
    (16*kt23SMEkt23*s23*z2*z3)/(kt23kt23*(1-z1)*z1)))/s23**2+ &
    (s123*((2*s123*SME*(1-2*z1)*z2)/((1-z1)*z1)+  &
    (8*kt2SMEkt2*(s13*z2-s123*(1-z2)*z2))/kt2kt2- &
    (16*kt3SMEkt3*z2**2*(s12*z3-s123*(1-z3)*z3))/(kt3kt3*(1- &
    z1)*z1)+4*(z2+(2*z2**2*(-z1+z3))/((1- &
    z1)*z1))*(((kt2SMEkt2*(s13*z2-s123*(1- &
    z2)*z2))/(kt2kt2*z2**2)+ &
    (kt23SMEkt23*s23)/(kt23kt23*z2*z3))*z3+(kt3SMEkt3*(s12*z3- &
    s123*(1-z3)*z3))/(kt3kt3*z3))))/(s12*s23)+(s123*(-2*s123*SME &
    -4*(-((kt1SMEkt1*(s23*z1-s123*(1-z1)*z1))/kt1kt1)+ &
    z3*((kt2SMEkt2*(s13*z2-s123*(1-z2)*z2))/(kt2kt2*z2)+ &
    z2*((kt23SMEkt23*s23)/(kt23kt23*z2*z3)+(kt3SMEkt3*(s12*z3- &
    s123*(1-z3)*z3))/(kt3kt3*z3**2))))))/(s12*s13))/4.d0 + &
! replace 2 with 3 and vice versa:
    (SME*(-1+(2*s123*(1-z1+2*z1**2))/(s23*(1-z1)*z1)+(2*s123*(1-  &
    z2))/(s13*(1-z1)*z1))+(s123*(-((SME*t231**2)/s123)+  &
    (16*kt23SMEkt23*s23*z3*z2)/(kt23kt23*(1-z1)*z1)))/s23**2+ &
    (s123*((2*s123*SME*(1-2*z1)*z3)/((1-z1)*z1)+  &
    (8*kt3SMEkt3*(s12*z3-s123*(1-z3)*z3))/kt3kt3- &
    (16*kt2SMEkt2*z3**2*(s13*z2-s123*(1-z2)*z2))/(kt2kt2*(1- &
    z1)*z1)+4*(z3+(2*z3**2*(-z1+z2))/((1- &
    z1)*z1))*(((kt3SMEkt3*(s12*z3-s123*(1- &
    z3)*z3))/(kt3kt3*z3**2)+ &
    (kt23SMEkt23*s23)/(kt23kt23*z3*z2))*z2+(kt2SMEkt2*(s13*z2- &
    s123*(1-z2)*z2))/(kt2kt2*z2))))/(s13*s23)+(s123*(-2*s123*SME &
    -4*(-((kt1SMEkt1*(s23*z1-s123*(1-z1)*z1))/kt1kt1)+ &
    z2*((kt3SMEkt3*(s12*z3-s123*(1-z3)*z3))/(kt3kt3*z3)+ &
    z3*((kt23SMEkt23*s23)/(kt23kt23*z3*z2)+(kt2SMEkt2*(s13*z2- &
    s123*(1-z2)*z2))/(kt2kt2*z2**2))))))/(s13*s12))/4.d0
!
  ans = qcd_cf*qcd_tr*Pabgqqb0 + qcd_ca*qcd_tr*Pnabgqqb0
!
end function Pgqqb0
! 
! This is \hat{P}^{\mu\nu}_{g_1g_2g_3} as presented in
! Eq. (66) of arXiv:hep-ph/9908523: 
function Pggg0(s123,s12,s13,s23,z1,z2,z3,kt1,kt2,kt3,SME,SMEmunu) result(ans)
use QCDparams
use particles
implicit none
!
  real(kind(1d0)) :: ans
!
  real(kind(1d0)) , intent(in) :: s123
  real(kind(1d0)) , intent(in) :: s12,s13,s23
  real(kind(1d0)) , intent(in) :: z1,z2,z3
  type(particle) , intent(in) :: kt1,kt2,kt3
  real(kind(1d0)) , intent(inout) :: SME
  real(kind(1d0)) , dimension(0:3,0:3) , intent(in) :: SMEmunu
!
  real(kind(1d0)) :: t123,t231,t312,t213,t321,t132
  real(kind(1d0)) :: kt1kt1,kt2kt2,kt3kt3, &
                     kt12kt12,kt13kt13,kt23kt23
  real(kind(1d0)) :: kt1SMEkt1,kt2SMEkt2,kt3SMEkt3, &
                     kt12SMEkt12,kt13SMEkt13,kt23SMEkt23
  type(mom) :: kt12,kt13,kt23
!
!
  ans = 0d0
!
! Definition is in Eq. (22) of arXiv:hep-ph/9908523:
  t123 = 2*(z1*s23 - z2*s13)/(z1 + z2) + (z1 - z2)/(z1 + z2)*s12
  t213 = -t123
  t231 = 2*(z2*s13 - z3*s12)/(z2 + z3) + (z2 - z3)/(z2 + z3)*s23
  t321 = -t231
  t312 = 2*(z3*s12 - z1*s23)/(z3 + z1) + (z3 - z1)/(z3 + z1)*s13
  t132 = -t312
!
! Constructing k_{\bot,rs} = k_{\bot,r}/z_r - k_{\bot,s}/z_s:
  kt12 = kt1%p/z1 - kt2%p/z2
  kt13 = kt1%p/z1 - kt3%p/z3
  kt23 = kt2%p/z2 - kt3%p/z3
!
! dot products of the various transverse momenta:
  kt1kt1 = kt1%p*kt1%p
  kt2kt2 = kt2%p*kt2%p
  kt3kt3 = kt3%p*kt3%p
  kt12kt12  = kt12*kt12
  kt13kt13  = kt13*kt13
  kt23kt23  = kt23*kt23
! contractions with the spin correlated ME:
  kt1SMEkt1 = kt1%p*SMEmunu*kt1%p
  kt2SMEkt2 = kt2%p*SMEmunu*kt2%p
  kt3SMEkt3 = kt3%p*SMEmunu*kt3%p
!
  kt12SMEkt12 = kt12*SMEmunu*kt12
  kt13SMEkt13 = kt13*SMEmunu*kt13
  kt23SMEkt23 = kt23*SMEmunu*kt23
!
  ans = &
! 123
    qcd_ca**2*((3*SME)/4.d0+(SME*t123**2-  &
    (16*kt12SMEkt12*s12*s123*z1*z2)/(kt12kt12*(1-  &
    z3)*z3))/(4.d0*s12**2)+(s123*(-(s123*SME*(-((1-2*(z1-  &
    z1**2))/(z2*z3))+(-1+2*(z1-z1**2)+4*z2*z3)/((1-z2)*(1- &
    z3))))/2.d0+2*z1*((kt2SMEkt2*(s13*z2-s123*(1-z2)*z2)*(-2+ &
    1/z3))/(kt2kt2*(1-z3))+(kt3SMEkt3*(-2+1/z2)*(s12*z3-s123*(1- &
    z3)*z3))/(kt3kt3*(1-z2)))+((2*(z2**2-z2**3))/(1-z3)- &
    3*z2*z3)*((kt2SMEkt2*(s13*z2-s123*(1-z2)*z2))/(kt2kt2*z2**2) &
    +(kt23SMEkt23*s23)/(kt23kt23*z2*z3)+(kt3SMEkt3*(s12*z3- &
    s123*(1-z3)*z3))/(kt3kt3*z3**2))))/(s12*s13)- &
    (s123*SME*((2*(1-z3)+4*z3**2)/(1-z3)-(1-2*(z3-z3**2))/((1- &
    z1)*z1)))/(s12*z3)) + &
! 231
    qcd_ca**2*((3*SME)/4.d0-(s123*SME*((2*(1-z1)+4*z1**2)/(1-z1)  &
    -(1-2*(z1-z1**2))/((1-z2)*z2)))/(s23*z1)+(SME*t231**2-  &
    (16*kt23SMEkt23*s123*s23*z2*z3)/(kt23kt23*(1-  &
    z1)*z1))/(4.d0*s23**2)+(s123*(-(s123*SME*(-((1-2*(z2- &
    z2**2))/(z1*z3))+(-1+2*(z2-z2**2)+4*z1*z3)/((1-z1)*(1- &
    z3))))/2.d0+2*z2*((kt1SMEkt1*(s23*z1-s123*(1-z1)*z1)*(-2+ &
    1/z3))/(kt1kt1*(1-z3))+(kt3SMEkt3*(-2+1/z1)*(s12*z3-s123*(1- &
    z3)*z3))/(kt3kt3*(1-z1)))+((kt1SMEkt1*(s23*z1-s123*(1- &
    z1)*z1))/(kt1kt1*z1**2)+(kt13SMEkt13*s13)/(kt13kt13*z1*z3)+ &
    (kt3SMEkt3*(s12*z3-s123*(1- &
    z3)*z3))/(kt3kt3*z3**2))*(-3*z1*z3+(2*(z3**2-z3**3))/(1- &
    z1))))/(s12*s23)) + &
! 312
    qcd_ca**2*((3*SME)/4.d0-(s123*SME*((2*(1-z2)+4*z2**2)/(1-z2)  &
    -(1-2*(z2-z2**2))/((1-z3)*z3)))/(s13*z2)+(SME*t312**2-  &
    (16*kt13SMEkt13*s123*s13*z1*z3)/(kt13kt13*(1-  &
    z2)*z2))/(4.d0*s13**2)+(s123*(((2*(z1**2-z1**3))/(1-z2)- &
    3*z1*z2)*((kt1SMEkt1*(s23*z1-s123*(1-z1)*z1))/(kt1kt1*z1**2) &
    +(kt12SMEkt12*s12)/(kt12kt12*z1*z2)+(kt2SMEkt2*(s13*z2- &
    s123*(1-z2)*z2))/(kt2kt2*z2**2))+2*((kt1SMEkt1*(s23*z1- &
    s123*(1-z1)*z1)*(-2+1/z2))/(kt1kt1*(1-z2))+(kt2SMEkt2*(-2+ &
    1/z1)*(s13*z2-s123*(1-z2)*z2))/(kt2kt2*(1-z1)))*z3- &
    (s123*SME*(-((1-2*(z3-z3**2))/(z1*z2))+(-1+4*z1*z2+2*(z3- &
    z3**2))/((1-z1)*(1-z2))))/2.d0))/(s13*s23)) + &
! 321
    qcd_ca**2*((3*SME)/4.d0-(s123*SME*((2*(1-z1)+4*z1**2)/(1-z1)  &
    -(1-2*(z1-z1**2))/((1-z3)*z3)))/(s23*z1)+(SME*t321**2-  &
    (16*kt23SMEkt23*s123*s23*z2*z3)/(kt23kt23*(1-  &
    z1)*z1))/(4.d0*s23**2)+(s123*(((kt1SMEkt1*(s23*z1-s123*(1- &
    z1)*z1))/(kt1kt1*z1**2)+(kt12SMEkt12*s12)/(kt12kt12*z1*z2)+ &
    (kt2SMEkt2*(s13*z2-s123*(1-  &
    z2)*z2))/(kt2kt2*z2**2))*(-3*z1*z2+(2*(z2**2-z2**3))/(1-z1)) &
    +2*((kt1SMEkt1*(s23*z1-s123*(1-z1)*z1)*(-2+1/z2))/(kt1kt1*(1 &
    -z2))+(kt2SMEkt2*(-2+1/z1)*(s13*z2-s123*(1- &
    z2)*z2))/(kt2kt2*(1-z1)))*z3-(s123*SME*(-((1-2*(z3- &
    z3**2))/(z1*z2))+(-1+4*z1*z2+2*(z3-z3**2))/((1-z1)*(1- &
    z2))))/2.d0))/(s13*s23)) + &
! 132
    qcd_ca**2*((3*SME)/4.d0-(s123*SME*((2*(1-z2)+4*z2**2)/(1-z2)  &
    -(1-2*(z2-z2**2))/((1-z1)*z1)))/(s13*z2)+(SME*t132**2-  &
    (16*kt13SMEkt13*s123*s13*z1*z3)/(kt13kt13*(1-  &
    z2)*z2))/(4.d0*s13**2)+(s123*(-(s123*SME*(-((1-2*(z1- &
    z1**2))/(z2*z3))+(-1+2*(z1-z1**2)+4*z2*z3)/((1-z2)*(1- &
    z3))))/2.d0+2*z1*((kt2SMEkt2*(s13*z2-s123*(1-z2)*z2)*(-2+ &
    1/z3))/(kt2kt2*(1-z3))+(kt3SMEkt3*(-2+1/z2)*(s12*z3-s123*(1- &
    z3)*z3))/(kt3kt3*(1-z2)))+((kt2SMEkt2*(s13*z2-s123*(1- &
    z2)*z2))/(kt2kt2*z2**2)+(kt23SMEkt23*s23)/(kt23kt23*z2*z3)+ &
    (kt3SMEkt3*(s12*z3-s123*(1- &
    z3)*z3))/(kt3kt3*z3**2))*(-3*z2*z3+(2*(z3**2-z3**3))/(1- &
    z2))))/(s12*s13)) + &
! 213
    qcd_ca**2*((3*SME)/4.d0+(SME*t213**2-  &
    (16*kt12SMEkt12*s12*s123*z1*z2)/(kt12kt12*(1-  &
    z3)*z3))/(4.d0*s12**2)+(s123*(-(s123*SME*(-((1-2*(z2-  &
    z2**2))/(z1*z3))+(-1+2*(z2-z2**2)+4*z1*z3)/((1-z1)*(1- &
    z3))))/2.d0+2*z2*((kt1SMEkt1*(s23*z1-s123*(1-z1)*z1)*(-2+ &
    1/z3))/(kt1kt1*(1-z3))+(kt3SMEkt3*(-2+1/z1)*(s12*z3-s123*(1- &
    z3)*z3))/(kt3kt3*(1-z1)))+((2*(z1**2-z1**3))/(1-z3)- &
    3*z1*z3)*((kt1SMEkt1*(s23*z1-s123*(1-z1)*z1))/(kt1kt1*z1**2) &
    +(kt13SMEkt13*s13)/(kt13kt13*z1*z3)+(kt3SMEkt3*(s12*z3- &
    s123*(1-z3)*z3))/(kt3kt3*z3**2))))/(s12*s23)- &
    (s123*SME*((2*(1-z3)+4*z3**2)/(1-z3)-(1-2*(z3-z3**2))/((1- &
    z2)*z2)))/(s12*z3))
!
end function Pggg0
!
! This is the strongly ordered quark splitting for:
! q_r \bar{r}_k r_t
! The notation is quite cryptic here, hence to be more
! understanable:
! z_k_t : z_{k,t}
! z_t_k : z_{t,k}
! z_kt_r : z_{\hat{kt},\hat{r}}
! z_r_kt : z_{\hat{r},\hat{kt}}
! s_kt_r : s_{\hat{kt},\hat{r}}
! s_r_kt_k_t : s_{\hat{r},k_{\bot,k,t}}
! kt_k_t2 : k_{\bot,k,t}^2
function PqqbqSO0(z_k_t,z_t_k,z_kt_r,z_r_kt, &
                  s_kt_r,s_r_kt_k_t,kt_k_t2) result(ans)
use QCDparams
implicit none
!
  real(kind(1d0)) :: ans
!
  real(kind(1d0)) , intent(in) :: z_k_t,z_t_k
  real(kind(1d0)) , intent(in) :: z_kt_r,z_r_kt
  real(kind(1d0)) , intent(in) :: s_r_kt_k_t,s_kt_r
  real(kind(1d0)) , intent(in) :: kt_k_t2
!
!
  real(kind(1d0)) , external :: Pqg0
!
  ans = qcd_tr*(Pqg0(z_r_kt,z_kt_r) &
      - 2*qcd_cf*z_k_t*z_t_k*(z_kt_r - s_r_kt_k_t**2 &
      / (kt_k_t2*s_kt_r)))
!
end function PqqbqSO0
!
! This is the strongly ordered quark splitting for:
! q_r g_k g_t
! z_k_t : z_{k,t}
! z_t_k : z_{t,k}
! z_kt_r : z_{\hat{kt},\hat{r}}
! z_r_kt : z_{\hat{r},\hat{kt}}
! s_kt_r : s_{\hat{kt},\hat{r}}
! s_r_kt_k_t : s_{\hat{r},k_{\bot,k,t}}
! kt_k_t2 : k_{\bot,k,t}^2
function PqggSO0(z_k_t,z_t_k,z_kt_r,z_r_kt, &
                 s_kt_r,s_r_kt_k_t,kt_k_t2) result(ans)
use QCDparams
implicit none
!
  real(kind(1d0)) :: ans
!
  real(kind(1d0)) , intent(in) :: z_k_t,z_t_k
  real(kind(1d0)) , intent(in) :: z_kt_r,z_r_kt
  real(kind(1d0)) , intent(in) :: s_r_kt_k_t,s_kt_r
  real(kind(1d0)) , intent(in) :: kt_k_t2
!
!
  real(kind(1d0)) , external :: Pqg0
!
  ans = 2*qcd_ca*(Pqg0(z_r_kt,z_kt_r) &
      * (z_k_t/z_t_k + z_t_k/z_k_t)  &
      + qcd_cf*z_k_t*z_t_k*(z_kt_r - s_r_kt_k_t**2/(kt_k_t2*s_kt_r)))
!
end function PqggSO0
!
! This is the strongly ordered gluon splitting for:
! g_r q_k \bar{q}_t
! Notation legend:
! z_k_t : z_{k,t}
! z_t_k : z_{t,k}
! z_kt_r : z_{\hat{kt},\hat{r}}
! z_r_kt : z_{\hat{r},\hat{kt}}
! s_kt_r : s_{\hat{kt},\hat{r}}
! s_r_kt_k_t : s_{\hat{r},k_{\bot,k,t}}
! kt_k_t2 : k_{\bot,k,t}^2
function PgqqbSO0(z_k_t,z_t_k,z_kt_r,z_r_kt,kt_tilde_k_t,kt_tilde_k_t2,kt_r_kt, &
                  s_kt_r,s_r_kt_k_t,kt_k_t2,kt_r_kt2,M2munu) result(ans)
use particles
use QCDparams
!
  real(kind(1d0)) :: ans
!
  real(kind(1d0)) , intent(in) :: z_k_t,z_t_k
  real(kind(1d0)) , intent(in) :: z_kt_r,z_r_kt
  real(kind(1d0)) , intent(in) :: s_r_kt_k_t,s_kt_r
  real(kind(1d0)) , intent(in) :: kt_k_t2,kt_r_kt2,kt_tilde_k_t2
  type(particle) , intent(in) :: kt_tilde_k_t,kt_r_kt
  real(kind(1d0)) , dimension(0:3,0:3) , intent(in) :: M2munu
!
!
  real(kind(1d0)) , external :: Pqq0sa
!
  ans = 2*qcd_ca*qcd_tr*( &
        (-M2munu(0,0) + M2munu(1,1) + M2munu(2,2) + M2munu(3,3)) &
      * (z_r_kt/z_kt_r + z_kt_r/z_r_kt + z_k_t*z_t_k*s_r_kt_k_t**2/(kt_k_t2*s_kt_r)) &
      + 4*z_k_t*z_t_k*z_kt_r/z_r_kt*(kt_tilde_k_t%p*M2munu*kt_tilde_k_t%p)/kt_tilde_k_t2) &
      - 4*qcd_ca*z_r_kt*z_kt_r*Pqq0sa(z_k_t,z_t_k)*(kt_r_kt%p*M2munu*kt_r_kt%p)/kt_r_kt2
!
end function PgqqbSO0
!
! This is the strongly ordered gluon splitting into:
! g_k g_t g_r:
function PgggSO0(z_k_t,z_t_k,z_kt_r,z_r_kt,kt_k_t,kt_tilde_k_t2,kt_r_kt, &
                 s_kt_r,s_r_kt_k_t,kt_k_t2,kt_r_kt2,M2munu) result(ans)
use particles
use QCDparams
use misc
use math
!
  real(kind(1d0)) :: ans
!
  real(kind(1d0)) , intent(in) :: z_k_t,z_t_k
  real(kind(1d0)) , intent(in) :: z_kt_r,z_r_kt
  real(kind(1d0)) , intent(in) :: s_r_kt_k_t,s_kt_r
  real(kind(1d0)) , intent(in) :: kt_k_t2,kt_r_kt2,kt_tilde_k_t2
  type(particle) , intent(in) :: kt_k_t,kt_r_kt
  real(kind(1d0)) , dimension(0:3,0:3) , intent(in) :: M2munu
!
  real(kind(1d0)) , external :: Pgg0sa
!
!  print *,"z_k_t: ",z_k_t
!  print *,"z_t_k: ",z_t_k
!  print *,"s_r_kt_k_t: ",s_r_kt_k_t
!  print *,"kt_k_t2: ",kt_k_t2
!  print *,"s_kt_r: ",s_kt_r
!
!  print *,"kt_r_kt: "
!  call PrintMom(kt_r_kt%p)
!
  ans = 4*qcd_ca**2*( &
        (-M2munu(0,0) + M2munu(1,1) + M2munu(2,2) + M2munu(3,3)) &
      * ((z_r_kt/z_kt_r + z_kt_r/z_r_kt)*(z_k_t/z_t_k + z_t_k/z_k_t) &
      - z_k_t*z_t_k*s_r_kt_k_t**2/(2*kt_k_t2*s_kt_r)) &
      - 2*z_k_t*z_t_k*z_kt_r/z_r_kt*(kt_k_t%p*M2munu*kt_k_t%p)/kt_tilde_k_t2) &
      - 4*qcd_ca*z_r_kt*z_kt_r*Pgg0sa(z_k_t,z_t_k)*(kt_r_kt%p*M2munu*kt_r_kt%p)/kt_r_kt2
!
!  ans = (kt_r_kt%p*M2munu*kt_r_kt%p)/kt_r_kt2
!  ans = (kt_k_t%p*M2munu*kt_k_t%p)/kt_k_t2
!
end function PgggSO0
!
! The soft gluon corresponds to leg r, ordering: q_i,g_r,g_s
! As defined in arXiv:hep-ph/0502226, Eq.(7.11).
function PSqgg(zi,zr,zs,sir,sis,srs) result(ans)
use QCDparams
implicit none
!
  real(kind(1d0)) :: ans
! 
  real(kind(1d0)) , intent(in) :: zi,zr,zs,sir,sis,srs
!
!
!
  ans = 2*qcd_cf*zi/(zr*sir) &
      + qcd_ca*(sis/(sir*srs) + zs/(srs*zr) - zi/(sir*zr))
!
end function PSqgg
!
! The soft gluon corresponds to leg r, ordering: g_r,q_s,\bar{q}_i
! As defined in arXiv:hep-ph/0502226, Eq.(7.12).
function PSgqq(zi,zr,zs,sir,sis,srs) result(ans)
use QCDparams
implicit none
!
  real(kind(1d0)) :: ans
! 
  real(kind(1d0)) , intent(in) :: zi,zr,zs,sir,sis,srs
!
!
!
  ans = 2*qcd_cf*sis/(sir*srs) &
      + qcd_ca*(zs/(srs*zr) + zi/(sir*zr) - sis/(sir*srs))
!
end function PSgqq
!
! The soft gluon corresponds to leg r, ordering: g_i,g_r,g_s
! As defined in arXiv:hep-ph/0502226, Eq.(7.13).
function PSggg(zi,zr,zs,sir,sis,srs) result(ans)
use QCDparams
implicit none
!
  real(kind(1d0)) :: ans
! 
  real(kind(1d0)) , intent(in) :: zi,zr,zs,sir,sis,srs
!
!
!
  ans = qcd_ca*(sis/(sir*srs) + zi/(sir*zr) + zs/(srs*zr))
!
end function PSggg
