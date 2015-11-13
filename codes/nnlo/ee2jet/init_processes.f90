subroutine init_processes(modselect)
use process
use subprocesses
use coupling
implicit none
!
  character*4 , intent(in) :: modselect
!
  integer :: istat
  integer :: ileg,iproc
  integer :: num_tmp
  integer , dimension(:,:) , allocatable , save :: flv_B_tmp, &
                                                   flv_R_tmp, &
                                                   flv_RR_tmp
!
  character*2 , dimension(nleg_born)   :: proc_LO
  character*2 , dimension(nleg_born+1) :: proc_NLO
  character*2 , dimension(nleg_born+2) :: proc_NNLO
!
  integer , parameter :: numproc_max = 2000
!
  interface
    subroutine gen_processes(nleg,proc_string,nproc,maxproc,proc_arr)
    implicit none
!
      integer , intent(in) :: nleg
      character*2 , dimension(:) , intent(in) :: proc_string
      integer , intent(out) :: nproc
      integer , intent(in) :: maxproc
      integer , dimension(:,:) , intent(out) :: proc_arr
!
    end subroutine gen_processes
  end interface
!
  border_as  = 0
  border_aEM = 2
!
! If modsecelct equals 'gen ' we generate all the subprocesses:
  if (modselect.eq.'gen ') then
    print *,"Generation of the subprocesses..."
! We define the process at the Born level:
    proc_LO(1) = 'e+'
    proc_LO(2) = 'e-'
    proc_LO(3) = 'j '
    proc_LO(4) = 'j '
! The first light parton resides at position 3:
    proc_firstlight = 3
!
! We allocate an array for the masses:
    allocate(fs_masses(nleg_born+2))
! For this process it is trivially zero:
    fs_masses = 0d0
!
    print *,"The process at the Born level: "
    do ileg=1,nleg_born
      write(*,'(A,1x)',advance='no') proc_LO(ileg)
      if (ileg.eq.2) write(*,'(A)',advance='no') "-> " 
      if (ileg.eq.nleg_born) write(*,*)
    end do
! Allocate an array to hold the Born subprocesses:
    allocate(flv_B_tmp(nleg_born,numproc_max),stat=istat)
    if (istat.ne.0) then
      print *,"Problem occured during allocating the Born temp aray..."
      stop
    end if
    call gen_processes(nleg_born,proc_LO,num_flv_LO, &
                       numproc_max,flv_B_tmp)
! For testing:
!    num_flv_LO = 2
!    flv_B_tmp(1,2) = -11
!    flv_B_tmp(2,2) =  11
!    flv_B_tmp(3,2) =   1
!    flv_B_tmp(4,2) =  -1
!    flv_B_tmp(5,2) =   0
! In this case the Born subprocesses coincide with the virtual ones:
    num_flv_NLO_V = num_flv_LO
! And the number of VV subprocesses coincide with the Born ones:
    num_flv_NNLO_VV = num_flv_LO
! We define the process at the NLO level:
    proc_NLO(1:nleg_born) = proc_LO(1:nleg_born)
    proc_NLO(5) = 'j '
! For testing:
!    proc_NLO(3) = 'u '
!    proc_NLO(4) = 'u~'
!    proc_NLO(5) = 'g '
!    proc_NLO(6) = 'g '
!
    print *,"The process at the NLO-R level: "
    do ileg=1,nleg_born+1
      write(*,'(A,1x)',advance='no') proc_NLO(ileg)
      if (ileg.eq.2) write(*,'(A)',advance='no') "-> " 
      if (ileg.eq.nleg_born+1) write(*,*)
    end do
! Allocate an array to hold the Real subprocesses:
    allocate(flv_R_tmp(nleg_born+1,numproc_max),stat=istat)
    if (istat.ne.0) then
      print *,"Problem occured during allocating the Real temp aray..."
      stop
    end if
    call gen_processes(nleg_born+1,proc_NLO,num_flv_NLO_R, &
                       numproc_max,flv_R_tmp)
! For testing:
!    num_flv_NLO_R = 2
!    flv_R_tmp(1,2) = -11
!    flv_R_tmp(2,2) =  11
!    flv_R_tmp(3,2) =   1
!    flv_R_tmp(4,2) =  -1
!    flv_R_tmp(5,2) =   1
!    flv_R_tmp(6,2) =  -1
! The number of RV subprocesses coincide with the NLO-R ones:
    num_flv_NNLO_RV = num_flv_NLO_R
! We define the process at the NNLO level:
    proc_NNLO(1:nleg_born+1) = proc_NLO(1:nleg_born+1)
    proc_NNLO(6) = 'j '
!    proc_NNLO(7) = 'g '
!
!    proc_NNLO(3) = 'u '
!    proc_NNLO(4) = 'u~'
!    proc_NNLO(5) = 'd '
!    proc_NNLO(6) = 'd~'
!
    print *,"The process at the NNLO-RR level: "
    do ileg=1,nleg_born+2
      write(*,'(A,1x)',advance='no') proc_NNLO(ileg)
      if (ileg.eq.2) write(*,'(A)',advance='no') "-> " 
      if (ileg.eq.nleg_born+2) write(*,*)
    end do
! Allocate an array to hold the RR subprocesses:
    allocate(flv_RR_tmp(nleg_born+2,numproc_max),stat=istat)
    if (istat.ne.0) then
      print *,"Problem occured during allocating the Real temp aray..."
      stop
    end if
    call gen_processes(nleg_born+2,proc_NNLO,num_flv_NNLO_RR, &
                       numproc_max,flv_RR_tmp)
  end if
!
  if (modselect.eq.'load') then
! We load the Born subprocesses:
    flv_LO(1:nleg_born,1:num_flv_LO) = flv_B_tmp(1:nleg_born, &
                                       1:num_flv_LO)
! We load the virtual subprocesses:
    flv_NLO_V(1:nleg_born,1:num_flv_LO) = flv_B_tmp(1:nleg_born, &
                                          1:num_flv_LO)
! We load the real subprocesses:
    flv_NLO_R(1:nleg_born+1, &
              1:num_flv_NLO_R) = flv_R_tmp(1:nleg_born+1, &
                                           1:num_flv_NLO_R)
! We load the VV subprocesses:
    flv_NNLO_VV(1:nleg_born,1:num_flv_LO) = flv_B_tmp(1:nleg_born, &
                                            1:num_flv_LO)
! We load the RV subprocesses:
    flv_NNLO_RV(1:nleg_born+1, &
                1:num_flv_NLO_R) = flv_R_tmp(1:nleg_born+1, &
                                             1:num_flv_NLO_R)
! We load the RR subprocesses:
    flv_NNLO_RR(1:nleg_born+2, &
                1:num_flv_NNLO_RR) = flv_RR_tmp(1:nleg_born+2, &
                                                1:num_flv_NNLO_RR)
! We loaded all processes, the arrays can be deallocated:
    deallocate(flv_B_tmp,flv_R_tmp,flv_RR_tmp)
  end if
!
end subroutine init_processes
!
function UserCheckConfig(nleg,subproc) result(ans)
implicit none
!
  logical :: ans
!
  integer , intent(in) :: nleg
  integer , dimension(:) , intent(in) :: subproc
!
  ans = .false.
!
! The initial state for this process is not partonic, hence
! a quark pair is required in the final state:
  if (sum(abs(subproc(3:nleg))).gt.0) ans = .true.
!
end function UserCheckConfig
!
! This routine has to be called in order to get initialized all
! the couplings and parameters for each and every contribution:
subroutine init_contribs()
use flags
use Born_data
use Real_data
use Virt_data
use VVirt_data
use RVirt_data
use RReal_data
use process
implicit none
!
!
  nleg_born = 4
!
  if (flg_LO) then
    call init_Born()
  end if
!
  if (flg_NLO_R) then
    call init_Real()
! We have to initialize the Born even when we only calculate 
! the real contribution since the Born and related correlated
! SMEs are needed to construct the counterterms:
    call init_Born()
  end if
!
  if (flg_NLO_V) then
    call init_Virt()
! The Born has to be initialized too because it is used to 
! convert between renormalization schemes:
    call init_Born()
  end if
!
  if (flg_NNLO_VV) then
    call init_VVirt()
! The Born and the virtual have to be initialized too because they are 
! used to for renormalization and for the construction of the 
! integrated subtraction terms.
    call init_Born()
    call init_Virt()
  end if
  if (flg_NNLO_RR) then
    call init_RReal()
! The real is only needed if subtraction terms are needed too:
    if (flg_NNLO_RR_A1.or.flg_NNLO_RR_A12) then
      call init_Real()
    end if
    if (flg_NNLO_RR_A2.or.flg_NNLO_RR_A12) then
      call init_Born()
    end if
  end if
  if (flg_NNLO_RV) then
    call init_RVirt()
! The real contribution has to be initialized since scheme
! conversion needs it:
    call init_Real()
! The virtual is only needed if subtraction terms are needed too:
    if (flg_NNLO_RV_A1) then
      call init_Born()
      call init_Virt()
    end if
  end if
!
end subroutine init_contribs
!
! This routine is called from init_Regions in order to change
! d0 and d0prime which can depend upon the process but we do not
! want to control their value through the input card:
subroutine init_procsubtract()
use regions
!
!
!
!
  sub_d0 = 3 
  sub_d0pr = 3
!
end subroutine init_procsubtract
