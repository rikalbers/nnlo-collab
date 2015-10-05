module avh_toyamp_kinem
  use avh_bint
!
! The unit to which messages are send
  integer ,private ,parameter :: nunit = 6
! Overall initialization parameter
  integer ,private :: nfinst_old=0
!
!
! Work arrays
  real(kind(1d0)) ,allocatable ,save :: wrkp(:,:),wrkm(:),wrks(:)
!
!
  type :: kinem_type
    real(kind(1d0)) ,dimension(0:3,-2:avh_bint_maxn) :: p
    real(kind(1d0)) ,dimension(    -2:avh_bint_maxn) :: m,s
    integer         ,dimension(    -2:avh_bint_maxn) :: b
    integer :: n1 = 0
  end type
!
!
contains
!
!
  subroutine kinem_init( nfinst )
!********************************************************************
!********************************************************************
  implicit none
  integer ,intent(in) :: nfinst
!
  if (nfinst_old.ge.nfinst) return
  nfinst_old = nfinst
!
  call avh_bint_base_init( nfinst )
!
  if (allocated( wrkp )) deallocate( wrkp )
  if (allocated( wrks )) deallocate( wrks )
  if (allocated( wrkm )) deallocate( wrkm )
  allocate( wrkp(0:3 ,avh_bint_b(nfinst+1)) )
  allocate( wrks(     avh_bint_b(nfinst+1)) )
  allocate( wrkm(     avh_bint_b(nfinst+1)) )
  end subroutine
!
!
  subroutine kinem_put_nfinst( kin ,nfinst )
!********************************************************************
!********************************************************************
  implicit none
  type(kinem_type) ,intent(inout) :: kin
  integer          ,intent(in)    :: nfinst
  integer :: ii
  kin%n1 = nfinst+1
  kin%b(-2) = avh_bint_b( kin%n1 )
  kin%b(-1) = avh_bint_b( 0 )
  kin%b( 0) = avh_bint_b( 0 )
  kin%b(1:kin%n1) = avh_bint_b(1:kin%n1)
  end subroutine
!
!
  subroutine kinem_closeall
!********************************************************************
!********************************************************************
  implicit none
  nfinst_old = 0
  if (allocated( wrkp )) deallocate( wrkp )
  if (allocated( wrks )) deallocate( wrks )
  if (allocated( wrkm )) deallocate( wrkm )
  end subroutine
!
!
  subroutine kinem_put_mass( kin ,mass,ii )
!********************************************************************
!********************************************************************
  implicit none
  type(kinem_type)  ,intent(inout) :: kin
  real(kind(1d0))   ,intent(in)    :: mass
  integer           ,intent(in)    :: ii
  if (kin%n1.le.0) then
    write(*,*) 'ERROR in toyamp_kinem_put_mass: ' ,&
               'you must call  put_nfinst  before  put_mass'
    stop
  endif
  if (ii.lt.-2.or.ii.eq.0.or.ii.gt.kin%n1) then
    write(*,*) 'ERROR in toyamp_kinem_put_mass: ii =',ii ,&
               ', but should be -2,-1, or between 0 and',kin%n1
    stop
  endif
  if (mass.lt.0d0) then
    write(*,*) 'ERROR in toyamp_kinem_put_mass: mass =',mass ,&
               ', but should not be negative'
    stop
  endif
  kin%m(ii) = mass
  kin%s(ii) = kin%m(ii)**2
  end subroutine 
!
!
  subroutine kinem_result( kin )
!*********************************************************************
!* Put the external momenta from "wrk" to "kin"
!*********************************************************************
  implicit none
  type(kinem_type) ,intent(inout) :: kin
  integer :: ii
  do ii=-2,kin%n1-1
    kin%p(0:3,ii) = wrkp(0:3,kin%b(ii))
  enddo
  end subroutine
!
!
  subroutine kinem_calcall( kin ,discard )
!********************************************************************
! Calculate all, internal, momenta and their squares
!********************************************************************
  implicit none
  type(kinem_type) ,intent(inout) :: kin
  logical          ,intent(out)   :: discard
  real(kind(1d0)) :: stot,small
  integer :: ii,i0,i1,i2,nout
!
  discard = .false.
!
  nout = kin%n1-1
  small = 1d6*epsilon(1d0)*sum( kin%p(0,1:nout) )**2
!
  do ii=-2,kin%n1-1
    if (ii.eq.0) cycle
    i0 = kin%b(ii)
    wrkp(0:3,i0) = kin%p(0:3,ii)
    wrks(i0) = wrkp(0,i0)**2 &
             - wrkp(1,i0)**2 - wrkp(2,i0)**2 - wrkp(3,i0)**2
    if (dabs(wrks(i0)-kin%s(ii)).gt.small) then
       if (nunit.gt.0) then
         write(nunit,*) &
         'WARNING from toyamp_kinem_calcall: momentum',ii,'on diffenent' ,&
         ' mass-shell than initially defined.'
       endif
!       wrkm(i0) = dsqrt(dabs( wrks(i0) ))
      wrks(i0) = kin%s(ii)
      wrkm(i0) = kin%m(ii)
    else
      wrks(i0) = kin%s(ii)
      wrkm(i0) = kin%m(ii)
    endif
  enddo
  i0 = kin%b(-2)
  wrkp(0:3,i0-1) = -wrkp(0:3,i0)
  wrks(    i0-1) =  wrks(    i0)
  wrkm(    i0-1) =  wrkm(    i0)
!
  stot = ( wrkp(0,1) + wrkp(0,kin%b(nout+1)) )**2
  small = epsilon(1d0)*stot
  do ii=0,nout
    i1 = kin%b(ii)
    do i2=1,i1-1
      i0 = i1+i2
      wrkp(0:3,i0) = wrkp(0:3,i1) + wrkp(0:3,i2)
      wrks(i0) = wrkp(0,i0)**2 &
               - wrkp(1,i0)**2 - wrkp(2,i0)**2 - wrkp(3,i0)**2
!      if (i0.eq.i0/2*2.and.wrks(i0).lt.small) then 
!        if (nunit.gt.0) write(nunit,*) &
!          'WARNING from toyamp_kinem_calcall: s(',i0,')/ehat^2 =' ,&
!          wrks(i0)/stot,', returning weight 0' 
!        discard = .true.
!        return
!      endif
      wrkm(i0) = dsqrt(dabs( wrks(i0) ))
    enddo
  enddo
  end subroutine
!
!
end module
