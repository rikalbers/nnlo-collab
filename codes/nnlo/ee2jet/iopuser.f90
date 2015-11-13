! This file contains the process dependent part of the I
! operator provided by the user...
subroutine InnloUser(p,yiQ,YikQ,         &
                     smeB,BornLaurent,   &
                     Bij,BijLaurent,     &
                     Bijkl,BijklLaurent, &
                     Ifit,IfitLaurent)
use particles
use Ioperators
use misc
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
  real(kind(1d0)) , dimension(:) , intent(in) :: yiQ
  real(kind(1d0)) , dimension(:,:) , intent(in) :: YikQ
  real(kind(1d0)) , intent(in) :: smeB
  real(kind(1d0)) , dimension(-4:2) :: BornLaurent
  real(kind(1d0)) , dimension(:,:) , intent(in) :: Bij
  real(kind(1d0)) , dimension(:,:,-4:) :: BijLaurent
  real(kind(1d0)) , dimension(:,:,:,:) , intent(in) :: Bijkl
  real(kind(1d0)) , dimension(:,:,:,:,-4:) :: BijklLaurent
  real(kind(1d0)) , intent(out) :: Ifit
  real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: IfitLaurent
!
  real(kind(1d0)) :: y13,y23
  real(kind(1d0)) , dimension(-4:2) :: Laurent
!
!
  y13 = YikQ(3,5)*yiQ(3)*yiQ(5)
  y23 = YikQ(4,5)*yiQ(4)*yiQ(5)
!
!  print *,"y13,y23: ",y13,y23
!  print *,"YikQ: ",YikQ(3,5),YikQ(4,5)
!  print *,"yiQ: ",yiQ(3:5)
!
  call Ifit3jet(y13,y23,Ifit,Laurent)
  if (present(IfitLaurent)) then
    IfitLaurent = SeriesProd(Laurent,BornLaurent)
    Ifit        = IfitLaurent(0)
  else
    Laurent = SeriesProd(Laurent,BornLaurent)
    Ifit    = IfitLaurent(0)
  end if
!
end subroutine InnloUser
