module avh_toyamp_vertex
  use avh_toyamp_model
  use avh_toyamp_tree ,only : tree_type
  use avh_toyamp_kinem
!
  private
  public :: vertex_type,vertex_init,vertex_wght
!
! Unit messages are send to
  integer ,parameter :: unitw=6
!
  real(kind(1d0)) ,parameter :: eps = 1d-32
!
!
  type :: vertex_type
    private
    integer :: p0,p1,p2
    integer :: f0,f1,f2
    integer :: o0,o1,o2
    logical :: putprop=.false.
  end type
!
!
contains
!
!
  subroutine vertex_init( tree,iv ,obj )
!********************************************************************
!********************************************************************
  use avh_bint
  implicit none
  type(tree_type)   ,intent(in)    :: tree
  type(vertex_type) ,intent(inout) :: obj
  integer           ,intent(in)    :: iv
  obj%p0 = tree%po( tree%ov(0,iv) )
  obj%p1 = tree%po( tree%ov(1,iv) )
  obj%p2 = tree%po( tree%ov(2,iv) )
  obj%f0 = tree%fo( tree%ov(0,iv) )
  obj%f1 = tree%fo( tree%ov(1,iv) )
  obj%f2 = tree%fo( tree%ov(2,iv) )
  obj%o0 =        ( tree%ov(0,iv) )
  obj%o1 =        ( tree%ov(1,iv) )
  obj%o2 =        ( tree%ov(2,iv) )
  obj%putprop = ( avh_bint_l(obj%p0).gt.1 .and. &
                  avh_bint_l(obj%p0).lt.tree%n1 )
  end subroutine
!
!
  subroutine vertex_wght( mdl ,kin ,obj ,osc,nsize )
!********************************************************************
!********************************************************************
  implicit none
  type(model_type)   ,intent(in)    :: mdl
  type(vertex_type)  ,intent(inout) :: obj
  type(kinem_type)   ,intent(inout) :: kin
  integer            ,intent(in)    :: nsize
  complex(kind(1d0)) ,intent(inout) :: osc(nsize)
  complex(kind(1d0)) :: prop
  real(kind(1d0))    :: mass,width
!
  prop = cmplx( 1d0 )
  if (obj%putprop) then
    mass  = mparticle( mdl ,obj%f0 )
    width = wparticle( mdl ,obj%f0 )
    if (width.gt.0d0) then
      prop = prop/cmplx( wrks(obj%p0)-mass*mass ,mass*width )
    else
      prop = prop/sqrt(cmplx( wrks(obj%p0)-mass*mass ,eps ))
    endif
  endif
!
  osc(obj%o0) = osc(obj%o0) + prop*osc(obj%o1)*osc(obj%o2)
!  write(6,*) obj%o0,prop,obj%o1,obj%o1,osc(obj%o0) !DEBUG
!
  end subroutine
!
!
end module
