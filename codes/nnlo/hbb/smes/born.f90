subroutine BornSME(iptrn,pin,Born)
use process
use particles
use QCDparams
use my_model
implicit none
!
  integer , intent(in) :: iptrn
  type(particle) , dimension(:) , intent(in) :: pin
  real(kind(1d0)) , intent(out) :: Born
!
  integer :: ipart
  logical , save :: init = .true.
!
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
  Born = 2d0*phys_yb**2*phys_Hmass**2*qcd_nc
!
end subroutine BornSME
!
subroutine BmunuSME(iptrn,ileg,pin,Bmunu)
use process
use particles
use QCDparams
use my_model
use misc
implicit none
!
  integer , intent(in) :: iptrn,ileg
  type(particle) , dimension(:) , intent(in) :: pin
  real(kind(1d0)) , dimension(0:3,0:3) , intent(out) :: Bmunu
!
  integer :: ipart
  logical , save :: init = .true.
  real(kind(1d0)) :: Born
  complex(kind(1d0)) , dimension(2,2) :: m2ij
!
!
  if (init) then
    init = .false.
  end if
!
!
end subroutine BmunuSME
!
subroutine BijSME(iptrn,pin,Bij)
use process
use particles
use QCDparams
use my_model
implicit none
!
  integer , intent(in) :: iptrn
  type(particle) , dimension(:) , intent(in) :: pin
  real(kind(1d0)) , dimension(:,:) , intent(out) :: Bij
!
  integer :: ipart
  logical , save :: init = .true.
  real(kind(1d0)) :: Born
!
!
  interface
    subroutine BornSME(iptrn,pborn,Born)
    use particles
    implicit none
!
      integer , intent(in) :: iptrn
      type(particle) , dimension(:) , intent(in) :: pborn
      real(kind(1d0)) , intent(out) :: Born
!
    end subroutine BornSME
  end interface
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
  call BornSME(iptrn,pin,Born)
!
  Bij = 0d0
  Bij(3,3) = qcd_cf*Born
  Bij(4,4) = Bij(3,3)
  Bij(3,4) = -Bij(3,3)
  Bij(4,3) = -Bij(3,3)
!
end subroutine BijSME
