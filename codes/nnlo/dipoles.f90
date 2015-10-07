! This source mainly contains the dipole generating routines used
! in the phasespace integrator Kaleu:
subroutine GenerateDipoles(flv,nflv_UB,flv_UB_arr,ndip,dips)
use subprocesses
use particles
use QCDparams
use utils
use regions
implicit none
!
  integer , dimension(:) , intent(in) :: flv
  integer , intent(in) :: nflv_UB
  integer , dimension(:,:) , intent(in) :: flv_UB_arr
  integer , intent(out) :: ndip
  integer , dimension(:,:) , allocatable , intent(out) :: dips
!
  integer :: istat
  integer :: i,j,k,n,nparton,ipart,partij
  integer :: ndip_tmp
  integer , dimension(size(flv)-1) :: flv_UB
  integer , dimension(:,:) , allocatable :: dips_tmp
!
!
!  call PrintSubproc(flv)
!
  ndip = 0
!
  n = size(flv)
!
! In principle dips should not be allocate at this moment but it
! is worth checking:
  if (allocated(dips)) deallocate(dips)
!
! The number of possible dipoles can be calculated a priori:
! The number of light partons:
  nparton = 0
  do ipart=1,n
    if (abs(flv(ipart)).gt.qcd_nf) cycle
    nparton = nparton + 1
  end do
!  print *,"Number of light partons: ",nparton
!
  ndip_tmp = nparton*(nparton-1)*(nparton-2)/2
!  print *,"Number of possible dipoles are: ",ndip_tmp
!
! Allocating a temporary array to hold dipoles:
  allocate(dips_tmp(3,ndip_tmp),stat=istat)
  if (istat.ne.0) then
    print *,"Allocation of dips_tmp went wrong..."
    stop "GenerateDipoles"
  end if
!
! Counting down dipoles , D_{ij,k}:
  do i=1,n-1
! Only QCD particles are considered:
    if (abs(flv(i)).gt.qcd_nf) cycle
    do j=i+1,n
! Only one parton can be in the initial state, chosen to be the
! first one:
      if (j.eq.2) cycle
      if (abs(flv(j)).gt.qcd_nf) cycle
! Looking for possible branchings:
! ISR:
      if (i.lt.3) then
! To have an allowed branching one of the partons has to be
! a gluon or the difference should be zero:
        partij = flv(i) - flv(j)
        if (.not.((flv(i).eq.0).or.(flv(j).eq.0).or. &
                  (partij.eq.0))) cycle
! FSR:
      else
! For the final state it still holds that one of the partons
! has to be a gluon or the sum of the partons should be zero:
        partij = flv(i) + flv(j)
        if (.not.((flv(i).eq.0).or.(flv(j).eq.0).or. &
                  (partij.eq.0))) cycle
      end if
!      print *,"pair found: ",i,j
!      print *,"parti,partj: ",ConvertFromPDG(flv(i)), &
!                              ConvertFromPDG(flv(j))
! Obtaining the underlying Born flavor configuration:
      call UndoBranchingMod(flv,flv_UB,i,j)
!      call PrintSubproc(flv_UB)
! We have to check whether this underlying Born is present among our
! subprocesses:
      if (.not.CheckEquivProc(n-1,nflv_UB,flv_UB_arr,flv_UB)) cycle
      do k=1,n
        if ((i.eq.k).or.(j.eq.k)) cycle
        if (abs(flv(k)).gt.qcd_nf) cycle
! If this point is reached a dipole is found, this has to be stored:
        ndip = ndip + 1
! When saving the dipoles we have to follow the convention of Kaleu, 
! -2,-1 are for initial state and 1,2,... for final state:
        if (i.lt.3) then
          dips_tmp(1,ndip) = i - 3
        else
          dips_tmp(1,ndip) = i - 2
        end if
        if (j.lt.3) then
          dips_tmp(2,ndip) = j - 3
        else
          dips_tmp(2,ndip) = j - 2
        end if
        if (k.lt.3) then
          dips_tmp(3,ndip) = k - 3
        else
          dips_tmp(3,ndip) = k - 2
        end if
      end do
    end do
  end do
!
! Right now we know the number of contributing dipoles the real array
! for dipoles can be allocated:
  allocate(dips(3,ndip),stat=istat)
  if (istat.ne.0) then
    print *,"Error during allocating dips..."
    stop "GenerateDipoles"
  end if
!
! Copy what is relevant:
  dips(1:3,1:ndip) = dips_tmp(1:3,1:ndip)
!
! Deallocating the temporary one:
  deallocate(dips_tmp)
!
end subroutine GenerateDipoles
