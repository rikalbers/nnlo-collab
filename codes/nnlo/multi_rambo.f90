module multi_rambo
  use rambo_module
  implicit none
!
  type(rambo_type) , dimension(:) , allocatable :: obj
!
end module multi_rambo
!
subroutine InitRamboArrays(numproc)
use multi_rambo
implicit none
!
  integer , intent(in) :: numproc
!
  integer :: istat
!
!
  allocate(obj(numproc),stat=istat)
  if (istat.ne.0) then
    print *,"Error in allocation of obj in multi_rambo..."
    stop
  end if
!
end subroutine InitRamboArrays
!
subroutine multi_Rambo_init(iproc,process,nfinst,ecm,option,masses)
use multi_rambo
use rambo_module
use kaleu_particles
implicit none
!
  integer , intent(in) :: iproc
  integer , intent(in) , dimension(-2:17) :: process
  integer , intent(in) :: nfinst
  real(kind(1d0)) , intent(in) :: ecm
  integer , intent(in) :: option
  real(kind(1d0)) , dimension(nfirst:nlast) , intent(in) :: masses
!
  integer , dimension(-2:17) :: prcss
!
!
! No need to store the sign of particles:
  prcss(-2:17) = abs(process(-2:17))
! Booking the process:
  call rambo_init(obj(iproc),prcss,nfinst,ecm,option,masses)
!
end subroutine multi_Rambo_init
!
subroutine multi_Rambo_gnrt(iproc,discard,x1rambo,x2rambo,prambo)
use multi_rambo
use rambo_module
implicit none
!
  integer , intent(in) :: iproc
  logical , intent(out) :: discard
  real(kind(1d0)) , intent(out) :: x1rambo,x2rambo
  real(kind(1d0)) , intent(out) :: prambo(0:3,-2:17)
!
!
!
! Right now there is no possibility to treat Bjorken xs:
  x1rambo = 1d0
  x2rambo = 1d0
!
  call rambo_gnrt(obj(iproc),prambo)
! For Rambo there is no discard, it is always false, it is only
! kept for consistency:
  discard = .false.
!
end subroutine multi_Rambo_gnrt
!
subroutine multi_Rambo_wght(iproc,weight)
use multi_rambo
use rambo_module
implicit none
!
  integer , intent(in) :: iproc
  real(kind(1d0)) , intent(out) :: weight
!
!
!
  call rambo_wght(obj(iproc),weight)
!
end subroutine multi_Rambo_wght
