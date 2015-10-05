module avh_kaleu_kinem
  use avh_bint
!
! The unit to which messages are send
  integer ,private ,parameter :: nunit = 6
! Overall initialization parameter
  integer ,private :: nfinst_old=0
! Hard cutoff, in units of ecm^2
  real(kind(1d0)) ,private ,parameter :: cutoff = 0d0 !1d-12
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
    real(kind(1d0)) ,allocatable :: smin(:)
    real(kind(1d0)) :: ecm=0d0 ,techcut=0d0
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
  subroutine kinem_put_nfinst( kin ,nfinst ,ecm )
!********************************************************************
!********************************************************************
  implicit none
  type(kinem_type) ,intent(inout) :: kin
  integer          ,intent(in)    :: nfinst
  real(kind(1d0))  ,intent(in)    :: ecm
  integer :: ii
  kin%n1 = nfinst+1
  kin%ecm = ecm
  kin%techcut = cutoff*ecm**2
  kin%b(-2) = avh_bint_b( kin%n1 )
  kin%b(-1) = avh_bint_b( 0 )
  kin%b( 0) = avh_bint_b( 0 )
  kin%b(1:kin%n1) = avh_bint_b(1:kin%n1)
  allocate( kin%smin(1:kin%b(kin%n1)) )
! Without initialization results in SIGFPE on Mac OSX, fixed by AK
  kin%smin(:) = cutoff
  end subroutine
!
!
  subroutine kinem_close( kin )
!********************************************************************
!********************************************************************
  implicit none
  type(kinem_type) ,intent(inout) :: kin
  if (allocated( kin%smin )) deallocate( kin%smin )
  end subroutine
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
    write(*,*) 'ERROR in kaleu_kinem_put_mass: ' ,&
               'you must call  put_nfinst  before  put_mass'
    stop
  endif
  if (ii.lt.-2.or.ii.eq.0.or.ii.gt.kin%n1) then
    write(*,*) 'ERROR in kaleu_kinem_put_mass: ii =',ii ,&
               ', but should be -2,-1, or between 0 and',kin%n1
    stop
  endif
  if (mass.lt.0d0) then
    write(*,*) 'ERROR in kaleu_kinem_put_mass: mass =',mass ,&
               ', but should not be negative'
    stop
  endif
  kin%m(ii) = mass
  kin%s(ii) = kin%m(ii)**2
  end subroutine 
!
!
  subroutine kinem_init_smin( kin )
!********************************************************************
! Put lower kinematical limits, necessary for positive invariants
!********************************************************************
  implicit none
  type(kinem_type) ,intent(inout) :: kin
  integer :: ii,i0,i1,i2,nout
  real(kind(1d0)) :: mass( kin%b(kin%n1) ),sminmin
!
  nout = kin%n1-1
  do ii=-1,nout
    if (ii.eq.0) cycle
    i1 = kin%b(ii)
    mass(i1) = kin%m(ii)
    kin%smin(i1) = mass(i1)**2
    do i2=1,i1-1
      i0 = i1+i2
      mass(i0) = mass(i1) + mass(i2)
      kin%smin(i0) = max( mass(i0)**2 ,2*kin%techcut )
    enddo
  enddo
!
  do ii=1,kin%b(kin%n1),2
    kin%smin(ii) = -2*kin%techcut
  enddo
!
  end subroutine
!
  subroutine kinem_updt_smin( kin ,smin )
!********************************************************************
!* Update lower kinematical limits on positive invariants, given a
!* matrix of lower limits on 2-particle invariants
!********************************************************************
  implicit none
  type(kinem_type) ,intent(inout) :: kin
  real(kind(1d0))  ,intent(in)    :: smin(-2:17,-2:17) ! (p_i+p_j)^2
  real(kind(1d0)) :: mass( kin%b(kin%n1) )
  integer :: ii,i0,nn
!
  nn = kin%n1
!  print *,"nn: ",kin%n1
!  print *,"dim: ",kin%b(kin%n1)
!  print *,"size: ",size(mass)
!
  call kinem_mass_smin( kin ,mass ,smin )
  do i0=2,kin%b(nn),2
!    print *,"i0: ",i0
!    print *,"mass: ",mass(i0)
!    print *,"smin: ",kin%smin(i0)
    kin%smin(i0) = max( kin%smin(i0) ,mass(i0)**2 )
  enddo
!
  do ii=1,nn-1
    i0 = kin%b(ii)+1
    kin%smin(i0) = min( kin%smin(i0) ,smin(-1,ii) ,smin(ii,-1) )
    i0 = kin%b(nn)-i0
    kin%smin(i0) = min( kin%smin(i0) ,smin(-2,ii) ,smin(ii,-2) )
  enddo
!
  end subroutine
!
  subroutine kinem_mass_smin( kin ,mass ,smin )
!********************************************************************
!* This routine does not alter "kin"
!* Get "mass", ie the minima on the positive invariants, given the
!* final-state masses in "kin" and given "smin"
!********************************************************************
  implicit none
  type(kinem_type) ,intent(in)  :: kin
  real(kind(1d0))  ,intent(in)  :: smin(-2:17,-2:17)
  real(kind(1d0))  ,intent(out) :: mass( 2**(kin%n1) )
  integer :: nn,nkinv,ii,i0,i1,i2,n1,n2,i0prev
  integer :: kinv(0:2,avh_bint_maxv)
  integer :: l1(avh_bint_maxn+1),l2(avh_bint_maxn)
!
  mass(1:2**(kin%n1)) = 0d0
!
  nn = kin%n1
  call avh_bint_kinv( kinv,nkinv ,nn )
!
  i0prev = 0
  do ii=1,nkinv
    i0 = kinv(0,ii)
    i1 = kinv(1,ii)
    i2 = kinv(2,ii)
    if (avh_bint_l(i0).eq.1.or.i0.ne.i0/2*2) then
      mass(i0) = 0d0
    elseif (avh_bint_l(i0).eq.2) then
      call avh_bint_esab( l1,n1,l2,n2 ,i0 ,nn )
      mass(i0) = max( kin%m(l1(1)) + kin%m(l1(2)) &
                     ,dsqrt(smin(l1(1),l1(2))) &
                     ,dsqrt(smin(l1(2),l1(1))) )
    else
      if (i0.ne.i0prev) then
        i0prev = i0
        mass(i0) = mass(i1) + mass(i2)
      else
        mass(i0) = max( mass(i0) ,mass(i1) + mass(i2) )
      endif
    endif
  enddo
  end subroutine
!
!
  subroutine kinem_inst( kin ,discard )
!********************************************************************
!* Put initial-state momenta
!********************************************************************
  implicit none
  type(kinem_type) ,intent(inout) :: kin
  logical          ,intent(out)   :: discard
  integer         :: i1,i2,i3,i4,ii,n1
  real(kind(1d0)) :: small,h0,h1,h2
  logical         :: init = .true.
  save small
  if (init) then
    init = .false.
    small = epsilon(1d0)**(14d0/16d0)
  endif
!
  discard = .false.
  n1 = kin%n1
!
! Put incoming momenta
  i1 = kin%b(-1)
  i2 = kin%b(-2)
  i3 = kin%b(n1)-2
  i4 = kin%b(n1)-1
  wrkp(0:3,i1) = kin%p(0:3,-1)
  wrkp(0:3,i2) = kin%p(0:3,-2)
  wrkp(0:3,i3) =-kin%p(0:3,-1)-kin%p(0:3,-2)
  wrkp(0:3,i4) =-kin%p(0:3,-2)
  do ii=-2,-1
    h2 = kin%p(1,ii)**2 + kin%p(2,ii)**2 + kin%p(3,ii)**2     
    h0 = kin%p(0,ii)**2
    h1 = h0 - h2
    h2 = h0 + h2
    if (h1.lt.-small*h2) then
!      if (nunit.gt.0) write(nunit,*) &
!        'ERROR in kaleu_kinem_inst: incoming momentum' ,&
!        ' is space-like: p(',ii,')^2 =',h1,', discard event.'
!      discard = .true.
!      return
    elseif (h1.gt.small*h2) then
!          if (nunit.gt.0) write(nunit,*) &
!             'WARNING from kaleu_kinem_inst: incoming momentum' ,&
!            ,' is not light-like: p(',ii,')^2 =',h1
    else
      h1 = 0d0
    endif
    if     (ii.eq.-1) then
      i1 = kin%b(ii)
      wrks(i1) = h1
      wrkm(i1) = dsqrt(h1)
    elseif (ii.eq.-2) then
      i1 = kin%b(ii)
      wrks(i1) = h1
      wrkm(i1) = dsqrt(h1)
      wrks(i1-1) = wrks(i1)
      wrkm(i1-1) = wrkm(i1)
      i1 = i1-2
      wrks(i1) = wrkp(0,i1)**2-wrkp(1,i1)**2-wrkp(2,i1)**2-wrkp(3,i1)**2
      wrkm(i1) = dsqrt(dabs(wrks(i1)))
    endif
  enddo
  do ii=1,kin%n1-1
    wrks(kin%b(ii)) = kin%s(ii)
    wrkm(kin%b(ii)) = kin%m(ii)
  enddo
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
  use avh_print
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
           'WARNING from kaleu_kinem_calcall: momentum ',trim(printint(ii)) &
          ,' on diffenent mass-shell than initially defined.'
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
      if (i0.eq.i0/2*2.and.wrks(i0).lt.small) then 
        if (nunit.gt.0) write(nunit,*) &
          'WARNING from kaleu_kinem_calcall: s(',trim(printint(i0)),')/ehat^2 =' ,&
          printsdbl(0,2,wrks(i0)/stot),', returning weight 0' 
        discard = .true.
        return
      endif
      wrkm(i0) = dsqrt(dabs( wrks(i0) ))
    enddo
  enddo
  end subroutine
!
!
end module
