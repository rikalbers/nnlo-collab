module avh_kaleu_strf
!  use avh_parni
  use avh_kaleu_grid
!
  private ! Everything is private, except the following list
  public :: strf_type,strf_init,strf_updt_xmin,strf_set_xmin &
           ,strf_gnrt,strf_wght,strf_collect,strf_close,strf_plot
!
!! The array dimension of random numbers from parni
!  integer ,parameter ,private :: size_dim = avh_parni_dim
! Unit to which messages are send
  integer ,parameter ,private :: nunit = 6
! Batch size and maximum number of channels for grids
  integer ,parameter ,private :: nbatch=1000 ,nchmax=200
! Fraction of channels fused in grid during optimization
  real(kind(1d0)) ,parameter ,private :: fusefrac=2d0/(2*nchmax-1)
! Minimum value for xmin=sqrt(tau0)
  real(kind(1d0)) ,parameter ,private :: cutoff=1d-6
!
  type :: strf_type
    private
    real(kind(1d0)) :: xmin,lntau0,pl,wl,wr,w
!    type(avh_parni_type) :: orr,oxx             ! Instances of parni
    type(grid_type) :: orr,oxx             ! Instances of 1d grids
  end type
!
!
contains
!
!
  subroutine strf_init( obj ,xmin )
!*********************************************************************
!* xmin is the value for sqrt(tau0)
!*
!* x1x2=exp(xx) is generated following (the normalized version of)
!*
!*             /  1/tau0  for  0 < x1x2 < tau0
!*  f(x1x2) = <
!*             \  1/x1x2  for  tau0 < x1x2 < 1
!*
!* with an adaptive grid underneath
!*********************************************************************
  implicit none
  type(strf_type) ,intent(inout) :: obj
  real(kind(1d0)) ,intent(in)    :: xmin
!
  if (xmin.lt.0d0.or.xmin.ge.1d0) then
    write(*,*) 'ERROR in kaleu_strf: xmin =',xmin ,&
               ', but should be between 0 and 1'
    stop
  elseif (xmin.lt.cutoff) then
!    if (nunit.gt.0) write(nunit,*) &
!      'WARNING from kaleu_strf: ' ,&
!      'overruling xmin =',xmin,' to minimum =',cutoff
    call strf_set_xmin( obj ,cutoff )
  else
    call strf_set_xmin( obj ,xmin   )
  endif
!
!
!  if (nunit.gt.0) then
  if (.false.) then
    write(nunit,*) 'MESSAGE from kaleu_strf: ' ,&
       'dont forget to put the flux factor 1/(2*ecm^2) yourself'
    write(nunit,*) 'MESSAGE from kaleu_strf: ' ,&
       'and dont put  1/(2*ehat^2)   but   1/(2*emc^2)'
  endif
!
!  call avh_parni_init( obj%orr ,11 ,1 ,nbatch,nchmax )
!  call avh_parni_init( obj%oxx ,11 ,1 ,nbatch,nchmax )
  call grid_init( obj%orr ,nbatch,nchmax,fusefrac ,'STRFrvar' )
  call grid_init( obj%oxx ,nbatch,nchmax,fusefrac ,'STRFxvar' )
  end subroutine
!
!
  subroutine strf_updt_xmin( obj ,xmin )
!*********************************************************************
!*********************************************************************
  implicit none
  type(strf_type) ,intent(inout) :: obj
  real(kind(1d0)) ,intent(in)    :: xmin
  if (xmin.gt.obj%xmin) call strf_set_xmin( obj ,xmin )
  end subroutine
!
!
  subroutine strf_set_xmin( obj ,xmin )
!*********************************************************************
!*********************************************************************
  implicit none
  type(strf_type) ,intent(inout) :: obj
  real(kind(1d0)) ,intent(in)    :: xmin
  real(kind(1d0)) :: tau0
  obj%xmin = xmin
  tau0 = obj%xmin**2
  obj%lntau0 = dlog( tau0 )
  obj%pl = 1d0/(1d0-obj%lntau0)
  obj%wl = tau0/obj%pl
  obj%wr = obj%lntau0/(obj%pl-1d0)
  if (nunit.gt.0) write(nunit,*) &
    'MESSAGE from kaleu_strf: xmin set to',obj%xmin
  end subroutine
!
!
  subroutine strf_gnrt( obj ,x1,x2 )
!*********************************************************************
!*********************************************************************
  implicit none
  type(strf_type) ,intent(inout) :: obj
  real(kind(1d0)) ,intent(out)   :: x1,x2
  real(kind(1d0)) :: rr,xx
!  real(kind(1d0)) :: xparni(size_dim)
!  call avh_parni_generate( obj%orr ,xparni )
!  rr = xparni(1)
!  call avh_parni_generate( obj%oxx ,xparni )
!  xx = xparni(1)
  call grid_gnrt( obj%orr ,rr )
  call grid_gnrt( obj%oxx ,xx )
!
  if (xx.le.0d0) then
    x1 = 1d0
    x2 = 1d0
    obj%w = 0d0
    return
  endif
!
  if (xx.le.obj%pl) then
    obj%w = xx*obj%wl
    xx = dlog( obj%w )
    obj%w = obj%wl/obj%w
  else
    xx = (xx-1d0)*obj%wr
    obj%w = obj%wr
  endif
  x1 = dexp( xx*(1d0-rr) )
  x2 = dexp( xx*(    rr) )
  obj%w = -obj%w*xx ! jacobian
  end subroutine
!
!
  subroutine strf_wght( obj ,weight )
!*********************************************************************
!*********************************************************************
  implicit none
  type(strf_type) ,intent(inout) :: obj
  real(kind(1d0)) ,intent(out)   :: weight
  real(kind(1d0)) :: ww
!  real(kind(1d0)) :: xparni(size_dim)
!
  weight = obj%w
!
!  call avh_parni_weight( obj%orr ,ww ,xparni )
!  weight = weight*ww
!  call avh_parni_weight( obj%oxx ,ww ,xparni )
!  weight = weight*ww
!
  call grid_wght( obj%orr ,ww ,-1d0 )
  weight = weight*ww
  call grid_wght( obj%oxx ,ww ,-1d0 )
  weight = weight*ww
  end subroutine
!
!
  subroutine strf_collect( obj ,weight )
!*********************************************************************
!*********************************************************************
  implicit none
  type(strf_type) ,intent(inout) :: obj
  real(kind(1d0)) ,intent(in)    :: weight
!
!  call avh_parni_adaptone( obj%orr ,weight )
!  call avh_parni_adaptone( obj%oxx ,weight )
  call grid_collect( obj%orr ,weight ,-1d0 )
  call grid_collect( obj%oxx ,weight ,-1d0 )
  end subroutine
!
!
  subroutine strf_close( obj )
!*********************************************************************
!*********************************************************************
  implicit none
  type(strf_type) ,intent(inout) :: obj
!  call avh_parni_delete( obj%orr )
!  call avh_parni_delete( obj%oxx )
  call grid_close( obj%orr )
  call grid_close( obj%oxx )
  end subroutine
!
!
  subroutine strf_plot( obj ,iunit )
!*********************************************************************
!*********************************************************************
  implicit none
  type(strf_type) ,intent(in) :: obj
  integer         ,intent(in) :: iunit
!  call avh_parni_plot( obj%orr ,iunit ,1 )
!  call avh_parni_plot( obj%oxx ,iunit ,2 )
  call grid_plot( obj%orr ,iunit )
  call grid_plot( obj%oxx ,iunit )
  end subroutine
!
!
end module

!*********************************************************************
!*********************************************************************

! program test
!   use avh_kaleu_strf
!   implicit none
!   type(strf_type) :: obj
!   integer ,parameter :: nev = 1000000
!   real(kind(1d0)) :: s0=0d0,s1=0d0,s2=0d0
!   real(kind(1d0)) :: ww,x1,x2,xmin
!   integer :: iev
! !
!   xmin = 0.3d0
!   call strf_init( obj ,xmin )
! !
!   do iev=1,nev
!     call strf_gnrt( obj ,x1,x2 )
!     call strf_wght( obj ,ww )
!     ww = ww*x1*x2
!     ww = ww*( 2*x1 + 3*x2*x2 )/2
! !    ww = ww/max( xmin*xmin ,x1*x2 )
!     call strf_collect( obj ,ww )
! !
!     s0 = s0 + 1d0
!     s1 = s1 + ww
!     s2 = s2 + ww*ww
!   enddo
!   call strf_plot( obj ,21 )
!   s1 = s1/s0
!   s2 = dsqrt( ( s2/s0-s1*s1 )/( s0-1d0 ) )
!   write(6,*) 'result:',s1,s2,s2/s1
! !
! end program
