module avh_kaleu_ranvar
  use avh_kaleu_model
!
! The unit to which messages are send
  integer ,private ,parameter :: nunit = 6
!
! The following variables can be changed with the routine ranvar_global
!
! Batch size for optimization choice above/below
  integer ,private :: nbatch_cut=1000
!
!
  type :: ranvar_type
    private
    real(kind(1d0)) :: prob=0.5d0,w0=0d0,w1=0d0 
    integer         :: idat=0
    logical         :: abov=.true.
    integer         :: fl=0
  end type
!
!
contains
!
!
  subroutine ranvar_global( nbatch_c )
!********************************************************************
!********************************************************************
  implicit none
  integer ,intent(in) :: nbatch_c
  nbatch_cut  = nbatch_c
  end subroutine
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
  return
  end subroutine
!
!
  subroutine ranvar_gnrt( model ,obj ,discard,ss ,slow,smid,supp )
!********************************************************************
!* Generation
!********************************************************************
  implicit none
  type(model_type)  ,intent(in)    :: model
  type(ranvar_type) ,intent(inout) :: obj
  logical           ,intent(inout) :: discard
  real(kind(1d0))   ,intent(out)   :: ss
  real(kind(1d0))   ,intent(in)    :: slow,smid,supp
  real(kind(1d0)) :: xx ,s0a,s1a,s0b,s1b
!
  if (dabs(slow).lt.supp) then
    s0b = slow
    s1b = smid
    s0a = smid
    s1a = supp
    if (supp.le.smid) then
      s0b = slow
      s1b = smid-slow ! assumes  slow  is very small
      s0a = smid-slow
      s1a = smid
    endif
  else!if (slow.lt.-dabs(supp)) then
    s0a = slow
    s1a = smid
    s0b = smid
    s1b = supp
    if (smid.le.slow) then
      s0a = smid
      s1a = smid+supp ! assumes |supp|=-supp is very small
      s0b = smid+supp
      s1b = supp
    endif
  endif
!
  call avh_random( xx )
  if (xx.ge.obj%prob) then
    obj%abov = .true.
    call avh_random( xx )
    call model_gnrt_above( model,obj%fl ,discard,ss ,xx,s0a,s1a )
  else
    obj%abov = .false.
    call avh_random( xx )
    call model_gnrt_below( model,obj%fl ,discard,ss ,xx,s0b,s1b )
  endif
  end subroutine
!
!
  subroutine ranvar_wght( model ,obj ,ww,ss ,slow,smid,supp )
!********************************************************************
!* Weight calculation
!********************************************************************
  implicit none
  type(model_type)  ,intent(in)    :: model
  type(ranvar_type) ,intent(inout) :: obj
  real(kind(1d0))   ,intent(out)   :: ww
  real(kind(1d0))   ,intent(in)    :: ss,slow,smid,supp
  real(kind(1d0)) :: xx ,s0a,s1a,s0b,s1b
!
  if (dabs(slow).lt.supp) then
    s0b = slow
    s1b = smid
    s0a = smid
    s1a = supp
    if (supp.le.smid) then
      s0b = slow
      s1b = smid-slow ! assumes  slow  is very small
      s0a = smid-slow
      s1a = smid
    endif
    obj%abov = (smid.lt.ss)
  else!if (slow.lt.-dabs(supp)) then
    s0a = slow
    s1a = smid
    s0b = smid
    s1b = supp
    if (smid.le.slow) then
      s0a = smid
      s1a = smid+supp ! assumes |supp|=-supp is very small
      s0b = smid+supp
      s1b = supp
    endif
    obj%abov = (smid.gt.ss)
  endif
!
  if (obj%abov) then
    call model_wght_above( model,obj%fl ,ww ,ss,s0a,s1a )
    ww = ww/(1d0-obj%prob)
  else
    call model_wght_below( model,obj%fl ,ww ,ss,s0b,s1b )
    ww = ww/obj%prob
  endif
  end subroutine
!
!
  subroutine ranvar_collect( obj ,weight )
!********************************************************************
!********************************************************************
  implicit none
  type(ranvar_type) ,intent(inout) :: obj
  real(kind(1d0))   ,intent(in)    :: weight
  real(kind(1d0)) :: hh
  if (obj%abov) then
    obj%w0 = obj%w0 + weight
  else
    obj%w1 = obj%w1 + weight
  endif
  obj%idat = obj%idat+1
  if (mod(obj%idat,nbatch_cut).eq.0) then
    hh = obj%w0 + obj%w1
    if (hh.ne.0d0) obj%prob = obj%w1/hh
!Security not to eliminate the channel too soon:
    obj%prob = max( obj%prob ,1d0/dble(obj%idat) )
  endif
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
