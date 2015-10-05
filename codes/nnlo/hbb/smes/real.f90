subroutine RealSME(iptrn,pin,smeR)
use process
use particles
use QCDparams
use my_model
use math
implicit none
!
  integer , intent(in) :: iptrn
  type(particle) , dimension(:) , intent(in) :: pin
  real(kind(1d0)) , intent(out) :: smeR
!
  integer :: ipart
  logical , save :: init = .true.
!
  real(kind(1d0)) :: s12,s13,s23,z,Q2
!
  real(kind(1d0)) , external :: Pqg0
!
  if (init) then
    init = .false.
  end if
!
  if (iptrn.ne.1) then
    print *,"We can only have one pattern..."
    print *,"iptrn: ",iptrn
    stop
  end if
!
  s12 = 2d0*pin(3)%p*pin(4)%p
  s13 = 2d0*pin(3)%p*pin(5)%p
  s23 = 2d0*pin(4)%p*pin(5)%p
  Q2  = 2d0*pin(1)%p*pin(2)%p
  z   = s12/Q2
!
  smeR = 1d0/pi*phys_yb**2*phys_Hmass**2*qcd_nc*(4d0*pi)**2 &
       * Pqg0(z,1d0-z)*(1d0/s13 + 1d0/s23)
!
end subroutine RealSME
