! mode is a character variable:
! mode = 'f' the full virtual part is needed,
! mode = 'i' the mu independent part is needed,
! mode = 'd' the mu dependent part is needed,
! mode = 'b' both the mu dependent and indepent part is calculated.
subroutine VirtSME(iptrn,pin,mur,mode,Virt,Virt_dep,Virt_indep, &
                   VirtLaurent)
use process
use particles
use QCDparams
use my_model
use math
implicit none
!
  integer , intent(in) :: iptrn
  type(particle) , dimension(:) , intent(in) :: pin
  real(kind(1d0)) , intent(in) :: mur
  character , intent(in) :: mode
  real(kind(1d0)) , intent(out) :: Virt
  real(kind(1d0)) , intent(out) :: Virt_dep,Virt_indep
  real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: VirtLaurent
!
  integer :: ipart
  logical , save :: init = .true.
  real(kind(1d0)) :: Born
  real(kind(1d0)) :: mu2,mH2,s12
!
!
  if (init) then
    init = .false.
  end if
!
  if (present(VirtLaurent)) VirtLaurent = 0d0
  Virt = 0d0
  Virt_dep = 0d0
  Virt_indep = 0d0
!
  mu2 = mur*mur
  mH2 = phys_Hmass**2
  s12 = 2d0*pin(1)%p*pin(2)%p
!
  Born = 2d0*phys_yb**2*mH2*qcd_nc
!
  if (present(VirtLaurent)) then
    VirtLaurent(-2) = Born/pi*qcd_cf
    VirtLaurent(-1) = Born/pi*qcd_cf*(log(mu2/s12) + 3d0/2d0)
    VirtLaurent( 0) = Born/pi*qcd_cf*(1d0 - pi**2/2d0 + 0.5d0*log(mu2/s12)**2)
    Virt = VirtLaurent( 0)
    Virt_dep = Born/pi*qcd_cf*0.5d0*log(mu2/s12)**2
    Virt_indep = Born/pi*qcd_cf*(1d0 - pi**2/2d0)
    return
  end if
!
  if ((mode.eq.'f').or.(mode.eq.'b')) then
    Virt_dep = Born/pi*qcd_cf*0.5d0*log(mu2/s12)**2
    Virt_indep = Born/pi*qcd_cf*(1d0 - pi**2/2d0)
    Virt = Virt_dep + Virt_indep
  elseif (mode.eq.'d') then
    Virt_dep = Born/pi*qcd_cf*0.5d0*log(mu2/s12)**2
  elseif (mode.eq.'i') then
    Virt_indep = Born/pi*qcd_cf*(1d0 - pi**2/2d0)
  end if
!
end subroutine VirtSME
