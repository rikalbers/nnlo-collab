module avh_toyamp
  use avh_toyamp_model
  use avh_toyamp_tree
  use avh_toyamp_kinem
  use avh_toyamp_vertex
!
  private ! everything is private, except the following list
  public :: toyamp_type &
           ,toyamp_put_process &
           ,toyamp_put_mom &
           ,toyamp_calc &
           ,toyamp_close &
           ,toyamp_closeall &
           ,toyamp_get_nfinst &
           ,toyamp_get_m &
           ,toyamp_get_s
!
! The unit to which messages are send
  integer ,private ,parameter :: nunit = 6
!
  type :: toyamp_type
    private
    type(kinem_type) :: kin
    type(vertex_type) ,allocatable :: vtx(:)
    integer :: nv,n1,no ! #vertices ,#external-1 ,#off-shell currents
  end type
!
contains
!
!
  subroutine toyamp_hello
!**********************************************************************
!**********************************************************************
  implicit none
  logical :: init = .true.
  if (init) then
    init = .false.
    write(*,'(a72)') '########################################################################'
    write(*,'(a72)') '#                                                                      #'
    write(*,'(a72)') '#                         You are using ToyAmp                         #'
    write(*,'(a72)') '#                                                                      #'
    write(*,'(a72)') '# author: Andreas van Hameren <hamerenREMOVETHIS@ifj.edu.pl>           #'
    write(*,'(a72)') '#   date: 30-07-2010                                                   #'
    write(*,'(a72)') '#                                                                      #'
    write(*,'(a72)') '########################################################################'
  endif
  end subroutine
!
!
  subroutine toyamp_put_process( mdl,mlv,obj ,process,nfinst ,i1,i2 ,cancel )
!**********************************************************************
!**********************************************************************
  use avh_bint
  implicit none
  type(model_type) ,intent(inout) :: mdl
  type(vertx_type) ,intent(inout) :: mlv ! this is not vertex_type !
  type(toyamp_type) ,intent(inout) :: obj
  integer          ,intent(in)    :: process(-2:17),nfinst
  integer ,optional,intent(in)    :: i1,i2
  logical ,optional,intent(out)   :: cancel
  real(kind(1d0)) :: ptoyamp(0:3,-2:17),stot,m1,m2,s1,s2
  type(tree_type) :: tree
  integer :: ii,jj,prcss(-2:17)
  logical :: discard
!
  call toyamp_hello
!
! Initialize some sizes
  call kinem_init( nfinst )
!
! No anti-particles defined in ToyAmp. Just check charge conservation.
  jj = 0
  do ii=-2,nfinst
    if     (ii.lt.0) then
      if (process(ii).le.13) jj = jj+sign(1,-process(ii))
    elseif (ii.gt.0) then
      if (process(ii).le.13) jj = jj+sign(1, process(ii))
    endif
  enddo
  if (jj.ne.0) then
    cancel = .true.
    return
  endif
  prcss(-2:nfinst) = abs(process(-2:nfinst))
!
! Construct tree
  call tree_checkinput( mdl ,prcss, nfinst )
  call tree_dress( mdl,mlv ,tree ,prcss,nfinst )
  call tree_checkoutput( mdl ,tree ,prcss, nfinst ,cancel )
  if (cancel) return
  if (present(i1).and.present(i2)) then
    call tree_ij( tree ,cancel ,i1,i2 )
    call tree_checkoutput( mdl ,tree ,prcss, nfinst ,cancel )
    if (cancel) return
  endif
!
! Allocate and intialize vertices
  obj%n1 = tree%n1
  obj%no = tree%no
  obj%nv = tree%nv
  allocate( obj%vtx(1:obj%nv) )
  do ii=1,obj%nv
    call vertex_init( tree,ii ,obj%vtx(ii) )
  enddo
!
! Put on-shell invariants
  call kinem_put_nfinst( obj%kin ,nfinst )
  do ii=-2,nfinst
    if (ii.eq.0) cycle
    call kinem_put_mass( obj%kin ,mparticle(mdl,prcss(ii)) ,ii )
  enddo
!
  end subroutine
!
!
  subroutine toyamp_put_mom( obj ,pkaleu )
!*********************************************************************
!* Put final state momenta
!*********************************************************************
  implicit none
  type(toyamp_type) ,intent(inout) :: obj
  real(kind(1d0))   ,intent(in)    :: pkaleu(0:3,-2:17)
  integer :: nfinst
  nfinst = obj%kin%n1-1
  obj%kin%p(0:3,-2:-1    ) = pkaleu(0:3,-2:-1)
  obj%kin%p(0:3, 1:nfinst) = pkaleu(0:3, 1:nfinst)
  end subroutine
!
!
  subroutine toyamp_calc( mdl,obj ,amp )
!*********************************************************************
!*********************************************************************
  implicit none
  type(model_type) ,intent(in)     :: mdl
  type(toyamp_type) ,intent(inout) :: obj
  complex(kind(1d0)) ,intent(out)  :: amp
  complex(kind(1d0)) :: osc(obj%no)
  complex(kind(1d0)) ,parameter :: one=(1d0,0d0) ,zero=(0d0,0d0)
  integer :: ii
  logical :: discard
  amp = zero
  call kinem_calcall( obj%kin ,discard )
  if (discard) return
  osc(        1:obj%n1 ) = one
  osc( obj%n1+1:obj%no ) = zero
  do ii=1,obj%nv
    call vertex_wght( mdl ,obj%kin ,obj%vtx(ii) ,osc,obj%no )
  enddo
  amp = osc( obj%no )
!  write(6,*) 'amp',obj%no,amp !DEBUG
  end subroutine
!
!
  subroutine toyamp_close( obj )
!*********************************************************************
!*********************************************************************
  implicit none
  type(toyamp_type) ,intent(inout) :: obj
  integer :: ii
  deallocate( obj%vtx )
  end subroutine
!
  subroutine toyamp_closeall
!*********************************************************************
!*********************************************************************
  implicit none
  call kinem_closeall
  end subroutine
!
!
  function toyamp_get_nfinst( obj ) result(value)
!*********************************************************************
!*********************************************************************
  type(toyamp_type) ,intent(in) :: obj
  integer :: value
  value = obj%kin%n1-1
  end function
!
  function toyamp_get_m( obj, ii ) result(value)
!*********************************************************************
!*********************************************************************
  type(toyamp_type) ,intent(in) :: obj
  integer          ,intent(in) :: ii
  real(kind(1d0)) :: value
  value = obj%kin%m( ii )
  end function
!
  function toyamp_get_s( obj, ii ) result(value)
!*********************************************************************
!*********************************************************************
  type(toyamp_type) ,intent(in) :: obj
  integer          ,intent(in) :: ii
  real(kind(1d0)) :: value
  value = obj%kin%s( ii )
  end function
!
!
end module
