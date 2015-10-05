module cuts_module
  use particles
  private
  public :: cuts_type,cuts_init,cuts_pass
!
  real(kind(1d0)) ,parameter :: cutoff=1d-3
  real(kind(1d0)) ,parameter :: twopi=6.2831853071795864769252867665590d0
!
  type :: cuts_type
    private
    integer :: nfinst=0
    real(kind(1d0)) ,allocatable :: ec(:),c1(:),c2(:),cc(:,:),gmas(:,:)
    real(kind(1d0)) ,allocatable :: pt(:),eta(:),dr(:,:)
  end type
!
contains
!
!
  include 'cuts.h90'
!
!
  subroutine cuts_alloc( cuts ,nfinst )
  type(cuts_type) ,intent(inout) :: cuts
  integer         ,intent(in)    :: nfinst
  call cuts_close( cuts )
  cuts%nfinst = nfinst
  allocate( cuts%ec(1:nfinst) )
  allocate( cuts%c1(1:nfinst) )
  allocate( cuts%c2(1:nfinst) )
  allocate( cuts%cc(1:nfinst,1:nfinst) )
  allocate( cuts%gmas(1:nfinst,1:nfinst) )
  allocate( cuts%pt(1:nfinst) )
  allocate( cuts%eta(1:nfinst) )
  allocate( cuts%dr(1:nfinst,1:nfinst) )
  cuts%ec(1:nfinst) = 0d0
  cuts%c1(1:nfinst) = 1d0
  cuts%c2(1:nfinst) = 1d0
  cuts%cc(1:nfinst,1:nfinst) = 1d0
  cuts%gmas(1:nfinst,1:nfinst) = cutoff
  cuts%pt(1:nfinst)  = 0d0
  cuts%eta(1:nfinst) = 9d99
  cuts%dr(1:nfinst,1:nfinst) = 0d0
  end subroutine
!
!
  subroutine cuts_close( cuts )
  type(cuts_type) ,intent(inout) :: cuts
  if (allocated(cuts%ec  )) deallocate( cuts%ec   )
  if (allocated(cuts%c1  )) deallocate( cuts%c1   )
  if (allocated(cuts%c2  )) deallocate( cuts%c2   )
  if (allocated(cuts%cc  )) deallocate( cuts%cc   )
  if (allocated(cuts%gmas)) deallocate( cuts%gmas )
  if (allocated(cuts%pt  )) deallocate( cuts%pt   )
  if (allocated(cuts%eta )) deallocate( cuts%eta  )
  if (allocated(cuts%dr  )) deallocate( cuts%dr   )
  cuts%nfinst = 0
  end subroutine
!
!
  subroutine putsmin( cuts ,process,nfinst ,masses ,smin )
  implicit none
  type(cuts_type) ,intent(in) :: cuts
  integer         ,intent(in) :: process(-2:17),nfinst
  real(kind(1d0)) ,intent(in) :: masses(nfirst:nlast)
  real(kind(1d0)) ,intent(out) :: smin(-2:17,-2:17)
  real(kind(1d0)) :: e1,e2,p1,p2,rm1,rm2,xx
  integer :: i1,j1,j2

! Determine cuts on 2-particle invariants
  do j1=   1,nfinst-1
  do j2=j1+1,nfinst
    rm1 = masses(abs(process(j1))); rm2 = masses(abs(process(j2)))
    if(rm1.eq.0.and.rm2.eq.0)then
      e1 = cuts%ec(j1); e2 = cuts%ec(j2)
      p1 = 0;           p2 = 0
      if (e1.gt.0) p1 = e1*dsqrt( 1-(rm1/e1)**2 )
      if (e2.gt.0) p2 = e2*dsqrt( 1-(rm2/e2)**2 )
      xx = rm1**2 + rm2**2 + 2*e1*e2 - 2*p1*p2*cuts%cc(j1,j2)
      smin(j1,j2) = max( xx ,cuts%gmas(j1,j2)**2 ,(rm1+rm2)**2 )
      smin(j2,j1) = smin(j1,j2)
    else
      smin(j1,j2) = max( cuts%gmas(j1,j2)**2 ,(rm1+rm2)**2 )
      smin(j2,j1) = smin(j1,j2)
    endif
  enddo
  enddo
  do j1=1,nfinst
    rm1 = masses(abs(process(j1)))
    j2 = -1
      rm2 = masses(abs(process(j2)))
      smin(j1,j2) = rm1**2 + rm2**2 - cuts%pt(j1)**2/2
      smin(j2,j1) = smin(j1,j2)
    j2 = -2
      rm2 = masses(abs(process(j2)))
      smin(j1,j2) = rm1**2 + rm2**2 - cuts%pt(j1)**2/2
      smin(j2,j1) = smin(j1,j2)
  enddo
!
  end subroutine
!
!
end module
