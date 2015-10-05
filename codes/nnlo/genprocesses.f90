! This source contains routines to generate subprocesses:
subroutine gen_processes(nleg,proc_string,nproc,maxproc,proc_arr)
use process
use utils
use QCDparams
implicit none
!
  integer , intent(in) :: nleg
  character*2, dimension(:) , intent(in) :: proc_string
  integer , intent(out) :: nproc
  integer , intent(in) :: maxproc
  integer , dimension(:,:) , intent(out) :: proc_arr
!
  integer :: ileg,iproc
  integer :: istat
  integer :: numsubproc
  integer , dimension(:) , allocatable :: subproc_template
  integer , dimension(:) , allocatable :: subproc
  character*2 particle
!
  integer , dimension(:,:) , allocatable :: lst_subprocesses
!
  interface
    function UserCheckConfig(nleg,subproc) result(ans)
    implicit none
    !
      integer , intent(in) :: nleg
      integer , dimension(:) , intent(in) :: subproc
    !
      logical :: ans
    !
    end function UserCheckConfig
!
    subroutine PutNormalOrder(nleg,subproc,subproc_out)
    implicit none
    !
      integer , intent(in) :: nleg
      integer , dimension(:) , intent(in) :: subproc
      integer , dimension(:) , intent(out) :: subproc_out
    !
    end subroutine PutNormalOrder
  end interface
!
  nproc = 0
!
! A template has to be created for the subprocesses:
  allocate(subproc_template(nleg),subproc(nleg),stat=istat)
  if (istat.ne.0) then
    print *,"Allocation went wrong for subproc_template or subproc"
    stop
  end if
! We allocate an array for the subprocesses:
  allocate(lst_subprocesses(nleg,maxproc),stat=istat)
  if (istat.ne.0) then
    print *,"Allocation went wrong for lst_subprocesses"
    stop
  end if
!
! We fill up the template with -10000 if the term is a 
! (anti)proton or jet otherwise the PDG code is used:
  do ileg=1,nleg
    if ((proc_string(ileg).eq.'p ').or. &
        (proc_string(ileg).eq.'p~').or. &
        (proc_string(ileg).eq.'j ').or. &
        (proc_string(ileg).eq.'J ')) then
      subproc_template(ileg) = -10000
    else
      subproc_template(ileg) = ConvertToPDG(proc_string(ileg))
    end if
  end do
!
  subproc = subproc_template
!
  where (subproc.eq.-10000) 
    subproc = qcd_nf
  end where
!
  numsubproc = 0
!
! We add all the other possibilities
  ileg = nleg
  do while (.true.)
!    print *,subproc
!    read(*,*)
    if (CheckConfig(nleg,subproc).and. &
        UserCheckConfig(nleg,subproc)) then
      if (.not.CheckEquivProc(nleg,numsubproc,lst_subprocesses, &
          subproc)) then
        numsubproc = numsubproc + 1
        if (numsubproc.gt.maxproc) then
          print *,"The number of subprocesses exceeds maxproc..."
          stop
        end if
! We store the subprocess in normal order:
        call PutNormalOrder(nleg,subproc, &
                            lst_subprocesses(:,numsubproc))
      end if
    end if
! When we pass the first light parton we can safely jump into the 
! initial state: 
!    if ((ileg.gt.2).and.(ileg.lt.proc_firstlight)) ileg=proc_firstlight
    if ((ileg.gt.2).and.(ileg.lt.proc_firstlight)) ileg = 2
! This is included for non-hadronic initial states:
    if (subproc_template(ileg).ne.-10000) then
      ileg = ileg - 1
      if (ileg.eq.0) exit
      cycle
    end if
    if (subproc(ileg).gt.-qcd_nf) then
      subproc(ileg) = subproc(ileg) - 1
      ileg = nleg
      cycle
    else
      if (ileg.eq.1) exit
      subproc(ileg) = qcd_nf
      ileg = ileg - 1
      cycle
    end if
  end do
!
  nproc = numsubproc
  proc_arr(1:nleg,1:nproc) = lst_subprocesses(1:nleg,1:nproc)
!
  deallocate(subproc_template,subproc,lst_subprocesses)
!
end subroutine gen_processes
!
! This routine takes a flavor config at puts the final state particles 
! in the order of non-QCD, massive partons, up-type pairs, down-type
! pairs and gluons:
subroutine PutNormalOrder(nleg,flv,flv_out)
use process
implicit none
!
  integer , intent(in) :: nleg
  integer , dimension(:) , intent(in) :: flv
  integer , dimension(:) , intent(out) :: flv_out
!
  integer :: ipart,jpart
  integer :: flv_tmp
!
  flv_out = flv
!
!  print *,"Original ordering: ",flv
! Our only concern is the light QCD partons:
  do ipart=proc_firstlight,nleg-1
! We simply put up type (anti)quarks first:
    if ((mod(flv_out(ipart),2).eq.0).or. &
        (flv_out(ipart).eq.0)) then
      do jpart=ipart+1,nleg
        if (mod(flv_out(jpart),2).ne.0) then
          call swap(ipart,jpart,flv_out)
          exit
        end if
      end do
    end if
  end do
! By putting the gluons at the very end we arrive at the
! the desired ordering:
  do ipart=nleg,proc_firstlight+1,-1
    if (flv_out(ipart).ne.0) then
      do jpart=proc_firstlight,ipart-1
        if (flv_out(jpart).eq.0) then
          call swap(ipart,jpart,flv_out)
        end if
      end do
    end if
  end do
!
! We put lower lying quarks first:
  do ipart=proc_firstlight,nleg-1
    if (flv_out(ipart).eq.0) cycle
    do jpart=ipart+1,nleg
      if (flv_out(jpart).eq.0) cycle
      if ((mod(abs(flv_out(ipart)),2).eq. &
             mod(abs(flv_out(jpart)),2)).and. &
          (abs(flv_out(ipart)).gt.abs(flv_out(jpart)))) then
        call swap(ipart,jpart,flv_out)
      end if
    end do
  end do
!
! We put everything into q-qbar pairs:
  do ipart=proc_firstlight,nleg-1
    if (flv_out(ipart).eq.0) exit
! Always try to swap if the first is an antiquark or the
! first in one flavor is antiquark:
    if (((ipart.eq.proc_firstlight).or. &
        (abs(flv_out(ipart-1)).ne.abs(flv_out(ipart)))).and. &
        flv_out(ipart).lt.0) then
      do jpart=ipart+1,nleg
        if (flv_out(ipart).eq.0) exit
        if (flv_out(ipart).eq.-flv_out(jpart)) then
          call swap(ipart,jpart,flv_out)
        end if
      end do
    end if
! Otherwise we only try to swap if two quarks of the same 
! flavor are adjacent:
    if ((flv_out(ipart).gt.0).and. &
        (flv_out(ipart).eq.flv_out(ipart-1))) then
      do jpart=ipart+1,nleg
        if (flv_out(ipart).eq.-flv_out(jpart)) then
          call swap(ipart,jpart,flv_out)
        end if
      end do
    end if
  end do
!
!  print *,flv_out
!  read(*,*)
!
  contains
  subroutine swap(i,j,flv)
  implicit none
!
  integer , intent(in) :: i,j
  integer , dimension(:) , intent(inout) :: flv
!
  integer :: flv_tmp
!
  flv_tmp = flv(i)
  flv(i)  = flv(j)
  flv(j)  = flv_tmp
!
  end subroutine swap
end subroutine PutNormalOrder
