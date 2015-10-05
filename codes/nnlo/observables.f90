! This file is a mere collection of various physical observables which
! can be useful in various physics analyses:
!
module observables
implicit none
!
!
interface get_pt
  module procedure get_pt_mom
  module procedure get_pt_part
end interface get_pt
!
interface get_pabs
  module procedure get_pabs_mom
  module procedure get_pabs_part
end interface get_pabs
!
interface get_invm
  module procedure get_invm_mom
  module procedure get_invm_part
end interface get_invm
!
interface get_rapidity
  module procedure get_rapidity_mom
  module procedure get_rapidity_part
end interface get_rapidity
!
interface get_azimuth
  module procedure get_azimuth_mom
  module procedure get_azimuth_part
end interface get_azimuth
!
interface get_dR
  module procedure get_dR_mm
  module procedure get_dR_pp
  module procedure get_dR_mp
  module procedure get_dR_pm
end interface get_dR
!
interface get_costheta
  module procedure get_costheta_mm
  module procedure get_costheta_pp
  module procedure get_costheta_mp
  module procedure get_costheta_pm
end interface get_costheta
!
contains
!
function get_pt_mom(p) result(pt)
use momenta
implicit none
!
  real(kind(1d0)) :: pt
!
  type(mom) , intent(in) :: p
!
  real(kind(1d0)) :: pt2
!
  pt2 = p%px**2 + p%py**2
!
  if (pt2.eq.0) then
    pt = 0d0
    return
  else
    pt = sqrt(abs(pt2))
  end if
  return
!
end function get_pt_mom
!
function get_pt_part(p) result(pt)
use particles
implicit none
!
  real(kind(1d0)) :: pt
!
  type(particle) , intent(in) :: p
!
!
  pt = get_pt_mom(p%p)
  return
!
end function get_pt_part
!
function get_pabs_mom(p) result(pabs)
use particles
implicit none
!
  real(kind(1d0)) :: pabs
!
  type(mom) , intent(in) :: p
!
!
  pabs = sqrt(p%px**2 + p%py**2 + p%pz**2)
!
end function get_pabs_mom
!
function get_pabs_part(p) result(pabs)
use particles
implicit none
!
  real(kind(1d0)) :: pabs
!
  type(particle) , intent(in) :: p
!
!
  pabs = get_pabs_mom(p%p)
!
end function get_pabs_part
!
function get_invm_mom(p) result(m)
use momenta
implicit none
!
  real(kind(1d0)) :: m
!
  type(mom) , intent(in) :: p
!
!
  m = sqrt(abs(p%E**2 - p%px**2 - p%py**2 - p%pz**2))
  return
!
end function get_invm_mom
!
function get_invm_part(p) result(m)
use particles
implicit none
!
  real(kind(1d0)) :: m
!
  type(particle) , intent(in) :: p
!
!
  m = get_invm_mom(p%p)
  return
!
end function get_invm_part
!
function get_rapidity_mom(p) result(y)
use momenta
implicit none
!
  real(kind(1d0)) :: y
!
  type(mom) , intent(in) :: p
!
  real(kind(1d0)) :: xplus,xminus
  real(kind(1d0)) , parameter :: tiny = 1d-10
!
  xplus  = p%E + p%pz
  xminus = p%E - p%pz
!
  if ((xplus.gt.tiny).and.(xminus.gt.tiny)) then
    if (xplus/xminus.gt.tiny) then
      y = 0.5d0*log(xplus/xminus)
    else
      y = sign(1d8,p%pz)
    end if
  else
    y = sign(1d8,p%pz)
  end if
  return
!
end function get_rapidity_mom
!
function get_rapidity_part(p) result(y)
use particles
implicit none
!
  real(kind(1d0)) :: y
!
  type(particle) , intent(in) :: p
!
  y = get_rapidity_mom(p%p)
  return
!
end function get_rapidity_part
!
function get_azimuth_mom(p) result(azi)
use math
use momenta
implicit none
!
  real(kind(1d0)) :: azi
!
  type(mom) , intent(in) :: p
!
  azi = atan2(p%py,p%px)
!
  return
!
end function get_azimuth_mom
!
function get_azimuth_part(p) result(azi)
use particles
implicit none
!
  real(kind(1d0)) :: azi
!
  type(particle) , intent(in) :: p
!
  azi = get_azimuth_mom(p%p)
  return
!
end function get_azimuth_part
!
function get_dR_mm(p1,p2) result(dR)
use math
use momenta
implicit none
!
  real(kind(1d0)) :: dR
!
  type(mom) , intent(in) :: p1,p2
!
  real(kind(1d0)) :: y1,y2
  real(kind(1d0)) :: azi1,azi2
  real(kind(1d0)) :: dphi
!
  y1 = get_rapidity(p1)
  y2 = get_rapidity(p2)
!
  azi1 = get_azimuth(p1)
  azi2 = get_azimuth(p2)
!
  dphi = abs(azi1 - azi2)
  if (dphi.gt.pi) dphi = 2d0*pi - dphi
!
  if ((dphi.lt.0d0).or.(dphi.gt.pi)) then
    print *,"Problem with dphi in  get_dR..."
    print *,"dphi = ",dphi
    stop
  end if
!
  dR = sqrt((y1 - y2)**2 + dphi**2)
  return
!
end function get_dR_mm
!
function get_dR_pp(p1,p2) result(dR)
use math
use particles
implicit none
!
  real(kind(1d0)) :: dR
!
  type(particle) , intent(in) :: p1,p2
!
!
  dR = get_dR_mm(p1%p,p2%p)
  return
!
end function get_dR_pp
!
function get_dR_mp(p1,p2) result(dR)
use math
use momenta
use particles
implicit none
!
  real(kind(1d0)) :: dR
!
  type(mom)      , intent(in) :: p1
  type(particle) , intent(in) :: p2
!
!
  dR = get_dR_mm(p1,p2%p)
  return
!
end function get_dR_mp
!
function get_dR_pm(p1,p2) result(dR)
use math
use momenta
use particles
implicit none
!
  real(kind(1d0)) :: dR
!
  type(particle) , intent(in) :: p1
  type(mom)      , intent(in) :: p2
!
!
  dR = get_dR_mp(p2,p1)
  return
!
end function get_dR_pm
!
function get_costheta_mm(p1,p2) result(costheta)
use particles
implicit none
!
  real(kind(1d0)) :: costheta
!
  type(mom) , intent(in) :: p1,p2
!
  real(kind(1d0)) :: p1abs,p2abs
!
!
  p1abs = sqrt(p1%px**2 + p1%py**2 + p1%pz**2)
  p2abs = sqrt(p2%px**2 + p2%py**2 + p2%pz**2)
!
  if ((p1abs.eq.0).or.(p2abs.eq.0)) then
    costheta = 1
    return
  end if
!
  costheta = (p1%E*p2%E - p1*p2)/(p1abs*p2abs)
!
end function get_costheta_mm
!
function get_costheta_pp(p1,p2) result(costheta)
use particles
implicit none
!
  real(kind(1d0)) :: costheta
!
  type(particle) , intent(in) :: p1,p2
!
!
  costheta = get_costheta_mm(p1%p,p2%p)
!
end function get_costheta_pp
!
function get_costheta_pm(p1,p2) result(costheta)
use particles
implicit none
!
  real(kind(1d0)) :: costheta
!
  type(particle) , intent(in) :: p1
  type(mom) , intent(in) :: p2
!
!
  costheta = get_costheta_mm(p1%p,p2)
!
end function get_costheta_pm
!
function get_costheta_mp(p1,p2) result(costheta)
use particles
implicit none
!
  real(kind(1d0)) :: costheta
!
  type(mom) , intent(in) :: p1
  type(particle) , intent(in) :: p2
!
!
  costheta = get_costheta_pm(p2,p1)
!
end function get_costheta_mp
!
! This routine calculates the C and D parameters from particles
! present in p only cosidering those in the final state for which
! |flv| <= nf, that is light partons, gluon is 0.
!
! The original version of this routine was taken from 
! Guenther Dissertori
subroutine CDpars(p,nf,Cpar,Dpar)
use particles
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
  integer , intent(in) :: nf
  real(kind(1d0)) , intent(out) :: Cpar,Dpar
!
  integer :: npart
  integer :: ipart
  real(kind(1d0)) :: pi
  real(kind(1d0)) :: txx,tyy,tzz,txy,txz,tyz
  real(kind(1d0)) :: sum
!
  integer :: jpart
  real(kind(1d0)) , save :: minCpar = 1d99
  real(kind(1d0)) :: Q2
!
  sum = 0d0
  txx = 0d0
  tyy = 0d0
  tzz = 0d0
  txy = 0d0
  txz = 0d0
  tyz = 0d0
!
! Size is taken directly from the array
  npart = size(p)
! Only considering final state particles
  do ipart=3,npart
! Skip, if not light parton
    if (abs(p(ipart)%flv).gt.nf) cycle
    pi  = sqrt(abs(p(ipart)%p%px**2 + p(ipart)%p%py**2 + p(ipart)%p%pz**2))
    sum = sum + pi
    if (pi.ne.0d0) then
      txx = txx + p(ipart)%p%px*p(ipart)%p%px/pi
      tyy = tyy + p(ipart)%p%py*p(ipart)%p%py/pi
      tzz = tzz + p(ipart)%p%pz*p(ipart)%p%pz/pi
      txy = txy + p(ipart)%p%px*p(ipart)%p%py/pi
      txz = txz + p(ipart)%p%px*p(ipart)%p%pz/pi
      tyz = tyz + p(ipart)%p%py*p(ipart)%p%pz/pi
    endif
  end do
!
! It is possible that all the components are zero, this can tipically happen
! in a discarded PS points where the momenta are is identically zero, to
! avoid NaNs an extra check is built in:
  if (max(abs(txx),abs(tyy),abs(tzz), &
          abs(txy),abs(txz),abs(tyz)).eq.0d0) then
    Cpar = 0d0
    return
  end if
!
  Cpar = txx*tyy + tyy*tzz + tzz*txx - txy*txy - txz*txz - tyz*tyz
  Cpar = 3d0*Cpar/sum/sum
!
!  if (Cpar.lt.minCpar) then
!    minCpar = Cpar
!    print *,"minCpar is adjusted to: ",minCpar
!    call PrintParts(p)
!    print *,"ts: ",txx,tyy,tzz,txy,txz,tyz
!    print *,"sum: ",sum
!    Q2 = 2d0*p(1)%p*p(2)%p
!    do ipart=1,npart-1
!      do jpart=ipart+1,npart
!        print *,"i,j,yij: ",ipart,jpart,2d0*p(ipart)%p*p(jpart)%p/Q2
!      end do
!    end do
!  end if
! 
  Dpar = txx*tyy*tzz + 2*txy*txz*tyz - tyy*(txz**2) - txx*(tyz**2) &
       - tzz*(txy**2)
  Dpar = 27d0*Dpar/sum/sum/sum
!
end subroutine CDpars
!
! This routine calculates the sphericity which only differs from the
! C parameter how it is normalized.
! It is worth mentioning that this observable is *NOT* IR-safe.
subroutine Sphericity(p,nf,S)
use particles
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
  integer , intent(in) :: nf
  real(kind(1d0)) , intent(out) :: S
!
  integer :: npart
  integer :: ipart
  real(kind(1d0)) :: pi
  real(kind(1d0)) :: txx,tyy,tzz,txy,txz,tyz
!
!
  txx = 0d0
  tyy = 0d0
  tzz = 0d0
  txy = 0d0
  txz = 0d0
  tyz = 0d0
!
! Size is taken directly from the array
  npart = size(p)
! Only considering final state particles
  do ipart=3,npart
! Skip, if not light parton
    if (abs(p(ipart)%flv).gt.nf) cycle
    pi  = sqrt(p(ipart)%p%px**2 + p(ipart)%p%py**2 + p(ipart)%p%pz**2)
    if (pi.ne.0d0) then
      txx = txx + p(ipart)%p%px*p(ipart)%p%px/(pi*pi)
      tyy = tyy + p(ipart)%p%py*p(ipart)%p%py/(pi*pi)
      tzz = tzz + p(ipart)%p%pz*p(ipart)%p%pz/(pi*pi)
      txy = txy + p(ipart)%p%px*p(ipart)%p%py/(pi*pi)
      txz = txz + p(ipart)%p%px*p(ipart)%p%pz/(pi*pi)
      tyz = tyz + p(ipart)%p%py*p(ipart)%p%pz/(pi*pi)
    endif
  end do
!
  S = txx*tyy + tyy*tzz + tzz*txx - txy*txy - txz*txz - tyz*tyz
  S = 3d0*S
!
end subroutine Sphericity
!
! If i or j is zero the total EEC is returned,
! Otherwise just for the i-j pair, the factor of 2 due to
! symmetry is also included.
function CalcEEC(i,j,nf,p) result(EEC)
use particles
implicit none
!
  real(kind(1d0)) :: EEC
!
  integer , intent(in) :: i,j
  integer , intent(in) :: nf
  type(particle) , dimension(:) , intent(in) :: p
!
  integer :: ipart,jpart,npart
  real(kind(1d0)) :: Q2,coschi
!
!
  EEC = 0
!
  Q2 = (p(1)%p + p(2)%p)**2
!
  if (Q2.eq.0) return
!
  if ((i.eq.0).or.(j.eq.0)) then
    npart = size(p)
    do ipart=3,npart-1
      if (abs(p(ipart)%flv).gt.nf) cycle
      do jpart=ipart+1,npart
        if (abs(p(jpart)%flv).gt.nf) cycle
        coschi = get_costheta(p(ipart),p(jpart))
! Factor of 2 is because only half of the sum is calculated:
        EEC = EEC &
            + 2*p(ipart)%p%E*p(jpart)%p%E/Q2*(1 - coschi**2)
      end do
    end do
    return
  end if
!
! Otherwise only for the specific pair:
  coschi = get_costheta(p(i),p(j))
  EEC = 2*p(i)%p%E*p(j)%p%E/Q2*(1 - coschi**2)
!
end function CalcEEC
!
function CalcThrust(p,nf) result(thrust)
use particles
implicit none
!
  real(kind(1d0)) :: thrust
!
  type(particle) , dimension(:) , intent(in) :: p
  integer , intent(in) :: nf
!
  integer :: n,i,j,k
  integer :: si,sj,sk
  real(kind(1d0)) :: maxabssp,pkpixpj,sumabsmom,abssp
  real(kind(1d0)) , dimension(3) :: nvec,pixpj
  type(mom) :: sumskpk,sumsp,maxsp
!
  thrust = 0
!
  maxabssp = 0
  nvec = 0
  sumabsmom = 0
!
  n = size(p)
! Only considering final state particles
  do i=3,n
! And only those which belong to massless quark flavors:
    if (abs(p(i)%flv).gt.nf) cycle
    sumabsmom = sumabsmom + get_pabs(p(i))
    do j=i+1,n
      if (abs(p(j)%flv).gt.nf) cycle
! pi x pj:
      pixpj = cross(p(i),p(j))
      sumskpk = 0d0
      do k=3,n
        if ((k.eq.i).or.(k.eq.j)) cycle
        if (abs(p(k)%flv).gt.nf) cycle
        pkpixpj = p(k)%p%px*pixpj(1) &
                + p(k)%p%py*pixpj(2) &
                + p(k)%p%pz*pixpj(3)
        sk = sign(1d0,pkpixpj)
        sumskpk = sumskpk + sk*p(k)%p
      end do
      do si=-1,1,2
        do sj=-1,1,2
          sumsp = sumskpk + si*p(i)%p + sj*p(j)%p
          abssp = get_pabs(sumsp)
          if (maxabssp.lt.abssp) then
            maxabssp = abssp
            nvec(1) = sumsp%px/abssp
            nvec(2) = sumsp%py/abssp
            nvec(3) = sumsp%pz/abssp
          end if
        end do
      end do
    end do
  end do
!
  if (sumabsmom.eq.0) return
!
  thrust = maxabssp/sumabsmom
!
end function CalcThrust
!
end module observables
