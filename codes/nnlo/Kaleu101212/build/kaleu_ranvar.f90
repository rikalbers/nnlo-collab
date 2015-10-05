module avh_kaleu_ranvar
  use avh_kaleu_model
!
! The unit to which messages are send
  integer ,private ,parameter :: nunit = 6
!
!
  type :: ranvar_type
    private
    integer :: fl
  end type
!
!
contains
!
!
  subroutine ranvar_init( obj ,fl )
!********************************************************************
!* Initialization
!********************************************************************
  implicit none
  type(ranvar_type) ,intent(inout) :: obj
  integer           ,intent(in)    :: fl
  obj%fl = fl
  end subroutine
!
!
  subroutine ranvar_close( obj )
!********************************************************************
!********************************************************************
  implicit none
  type(ranvar_type) ,intent(inout) :: obj
  obj%fl = 0
  end subroutine
!
!
  subroutine ranvar_gnrt( model ,obj ,discard,ss ,slow,supp )
!********************************************************************
!* Generation
!********************************************************************
  implicit none
  type(model_type)  ,intent(in)    :: model
  type(ranvar_type) ,intent(inout) :: obj
  logical           ,intent(inout) :: discard
  real(kind(1d0))   ,intent(out)   :: ss
  real(kind(1d0))   ,intent(in)    :: slow,supp
  real(kind(1d0)) :: xx
!
  call avh_random( xx )
  call model_gnrt_above( model,obj%fl ,discard,ss ,xx,slow,supp )
  end subroutine
!
!
  subroutine ranvar_wght( model ,obj ,ww,ss ,slow,supp )
!********************************************************************
!* Weight calculation
!********************************************************************
  implicit none
  type(model_type)  ,intent(in)    :: model
  type(ranvar_type) ,intent(inout) :: obj
  real(kind(1d0))   ,intent(out)   :: ww
  real(kind(1d0))   ,intent(in)    :: ss,slow,supp
!
  call model_wght_above( model,obj%fl ,ww ,ss,slow,supp )
  end subroutine
!
!
  subroutine ranvar_collect( obj ,weight )
!********************************************************************
!********************************************************************
  implicit none
  type(ranvar_type) ,intent(inout) :: obj
  real(kind(1d0))   ,intent(in)    :: weight
  return
  end subroutine
!
!
  subroutine ranvar_plot( obj ,iunit )
!********************************************************************
!********************************************************************
  implicit none
  type(ranvar_type) ,intent(in) :: obj
  integer           ,intent(in) :: iunit
  return
  end subroutine
!
!
end module
