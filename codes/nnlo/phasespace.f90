! This source contains routines related to the phase space
! generation
module KaleuPS
implicit none
!
  integer :: PDF_option
  real(kind(1d0)) :: ymin
  real(kind(1d0)) :: smin(-2:17,-2:17)
  real(kind(1d0)) , dimension(:) , allocatable :: masses,widths
  character(20)   , dimension(:) , allocatable :: labels
  logical :: onlyqcd,withqcd,withhiggs
!
  integer :: nbatchBkin,nbatchRkin,nbatchRRkin
  integer :: nstepBkin,nstepRkin,nstepRRkin
  integer :: niterBkin,niterRkin,niterRRkin
  integer :: neventsBkin,neventsRkin,neventsRRkin
  integer :: noptimBkin,noptimRkin,noptimRRkin
  real(kind(1d0)) :: thrsBkin,thrsRkin,thrsRRkin
  logical :: optim = .true.
!
contains
!
function PDGtoKaleu(PDGid) result(Kid)
implicit none
!
  integer , intent(in) :: PDGid
!
  integer :: Kid
!
! d quark:
  if (abs(PDGid).eq.1) then
    Kid = sign(4,PDGid)
! u quark:
  elseif (abs(PDGid).eq.2) then
    Kid = sign(3,PDGid)
! s quark:
  elseif (abs(PDGid).eq.3) then
    Kid = sign(8,PDGid)
! c quark:
  elseif (abs(PDGid).eq.4) then
    Kid = sign(7,PDGid)
! b quark:
  elseif (abs(PDGid).eq.5) then
    Kid = sign(12,PDGid)
! t quark:
  elseif (abs(PDGid).eq.6) then
    Kid = sign(11,PDGid)
! gluon:
  elseif (abs(PDGid).eq.0) then
    Kid = 17
! electron:
  elseif (abs(PDGid).eq.11) then
    Kid = sign(2,PDGid)
  else
    write(*,'(A,I3)') "No Kaleu alias can be given for: ",PDGid
    stop
  end if
!
end function PDGtoKaleu
!
subroutine init_PSgen()
use process
use collider
use subprocesses
use flags
use input
implicit none
!
!
  integer :: iproc,jproc,ipart
  integer :: numproc
  integer , dimension(-2:17) :: subproc
!
  integer :: ndip,ndipsqr
  integer , dimension(:,:) , allocatable :: dip_tmp,dipsqr_tmp
!
  real(kind(1d0)) , parameter :: alphamax_kaleu = 0.999d0
!
  interface
    subroutine GenerateDipoles(flv,nflv_UB,flv_UB_arr,ndip,dips)
    implicit none
!
      integer , dimension(:) , intent(in) :: flv
      integer , intent(in) :: nflv_UB
      integer , dimension(:,:) , intent(in) :: flv_UB_arr
      integer , intent(out) :: ndip
      integer , dimension(:,:) , allocatable , intent(out) :: dips
!
    end subroutine GenerateDipoles
!
    subroutine GenerateDipolesSquared(flv,nflv_UB,flv_UB_arr, &
                                      ndipsqr,dipsqr)
    implicit none
!
      integer , dimension(:) , intent(in) :: flv
      integer , intent(in) :: nflv_UB
      integer , dimension(:,:) , intent(in) :: flv_UB_arr
      integer , intent(out) :: ndipsqr
      integer , dimension(:,:) , allocatable , intent(out) :: dipsqr
!
    end subroutine GenerateDipolesSquared
  end interface
!
  print *,"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
  print *,"@                                                          @"
  print *,"@       We are going to use KALEU as PS generator...       @"
  print *,"@                                                          @"
  print *,"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
!
! We have to calculate the total number of subprocesses independent
! of the multiplicity to set up the needed arrays for Kaleu:
  numproc = num_flv_ch_Bkin + num_flv_ch_Rkin + num_flv_ch_RRkin
!
! We have to allocate those arrays which will hold information for the
! subprocesses in Kaleu:
  call InitKaleuArrays(numproc)
!
! We perform process-dependent initialization too:
  call Kaleu_initmyprocess
!
! We read in generator dependent parameters:
  nbatchBkin  = nnloinput("#nbatch")
  nbatchRkin  = nbatchBkin
  nbatchRRkin = nbatchBkin
!
  nstepBkin  = nnloinput("#nstep")
  nstepRkin  = nstepBkin
  nstepRRkin = nstepBkin
!
  thrsBkin   = nnloinput("#thrs")
  thrsRkin   = thrsBkin
  thrsRRkin  = thrsBkin
!
  niterBkin  = nnloinput("#niter")
  niterRkin  = niterBkin
  niterRRkin = niterBkin
!
!
  if ((nbatchBkin.lt.0).or.(nnloinput("#nbatchBkin").gt.0)) then
    nbatchBkin = nnloinput("nbatchBkin")
  end if
  if ((nbatchRkin.lt.0).or.(nnloinput("#nbatchRkin").gt.0)) then
    nbatchRkin = nnloinput("nbatchRkin")
  end if
  if ((nbatchRRkin.lt.0).or.(nnloinput("#nbatchRRkin").gt.0)) then
    nbatchRRkin = nnloinput("nbatchRRkin")
  end if
!
  if ((nstepBkin.lt.0).or.(nnloinput("#nstepBkin").gt.0)) then
    nstepBkin = nnloinput("nstepBkin")
  end if
  if ((nstepRkin.lt.0).or.(nnloinput("#nstepRkin").gt.0)) then
    nstepRkin = nnloinput("nstepRkin")
  end if
  if ((nstepRRkin.lt.0).or.(nnloinput("#nstepRRkin").gt.0)) then
    nstepRRkin = nnloinput("nstepRRkin")
  end if
!
  noptimBkin  = nbatchBkin*nstepBkin
  noptimRkin  = nbatchRkin*nstepRkin
  noptimRRkin = nbatchRRkin*nstepRRkin
!
  if ((thrsBkin.lt.0).or.(nnloinput("#thrsBkin").gt.0)) then
    thrsBkin = nnloinput("thrsBkin")
  end if
  if ((thrsRkin.lt.0).or.(nnloinput("#thrsRkin").gt.0)) then
    thrsRkin = nnloinput("thrsRkin")
  end if
  if ((thrsRRkin.lt.0).or.(nnloinput("#thrsRRkin").gt.0)) then
    thrsRRkin = nnloinput("thrsRRkin")
  end if
!
  if ((niterBkin.lt.0).or.(nnloinput("#niterBkin").gt.0)) then
    niterBkin = nnloinput("niterBkin")
  end if
  if ((niterRkin.lt.0).or.(nnloinput("#niterRkin").gt.0)) then
    niterRkin = nnloinput("niterRkin")
  end if
  if ((niterRRkin.lt.0).or.(nnloinput("#niterRRkin").gt.0)) then
    niterRRkin = nnloinput("niterRRkin")
  end if
!
! We need the number of events needed to be generated:
! If only nevents supplied all the contributions are
! calculated with the same number of events:
  neventsBkin  = nnloinput("#nevents")
!
  neventsRkin  = neventsBkin
  neventsRRkin = neventsBkin
!
  if ((neventsBkin.lt.0).or.(nnloinput("#neventsBkin").gt.0)) then
    neventsBkin = nnloinput("neventsBkin")
  end if
  if ((neventsRkin.lt.0).or.(nnloinput("#neventsRkin").gt.0)) then
    neventsRkin = nnloinput("neventsRkin")
  end if
  if ((neventsRRkin.lt.0).or.(nnloinput("#neventsRRkin").gt.0)) then
    neventsRRkin = nnloinput("neventsRRkin")
  end if
!
  subproc = 0
!
! We calculate the minimal kinematic invariant from the 
! dimensionless one:
  smin = ymin * stot
!
! We set up all the subprocesses:
!=======================================================================
! Born kinematics:
!=======================================================================
  if (flg_LO.or.flg_NLO_V.or.flg_NNLO_VV) then
    do iproc=1,num_flv_ch_Bkin
! We have to convert the flavor list into a form which is
! acceptible for Kaleu:
      do ipart=1,nleg_born
! This extra statement is needed since in Kaleu the ordering is
! different: 
! initial state: -2:-1
! final state:    1:nfinst
        if (ipart.le.2) then
          subproc(ipart-3) = PDGtoKaleu(flv_ch_Bkin(ipart,iproc))
        else
          subproc(ipart-2) = PDGtoKaleu(flv_ch_Bkin(ipart,iproc))
        end if
      end do
!
!      print *,"iproc: ",iproc
!      print *,"flv: ",subproc
!
      jproc = iproc
!
! For Born-like kinematics no Dipole has to be defined:
      ndip = 0
!
! Note that the alphamax here is set to one and ymin is set
! at a different point, if the PS point is deeper than ymin it is
! returned to the generator was a point with zero weight, hence to
! optimize it away:
      call multi_kaleu_init(iproc,subproc,nleg_born-2,          &
                            rstot,PDF_option,                   &
                            masses,widths,labels,               &
                            onlyqcd,withqcd,withhiggs,          &
                            nbatchBkin,nstepBkin,thrsBkin,smin, &
                            alphamax_kaleu,alphamax_kaleu,      &
                            alphamax_kaleu,alphamax_kaleu,      &
                            ndip,dip_tmp)
    end do
  end if
!=======================================================================
! Real kinematics:
!=======================================================================
  if (flg_NLO_R.or.flg_NNLO_RV) then
    do iproc=1,num_flv_ch_Rkin
! We have to convert the flavor list into a form which is
! acceptible for Kaleu:
      do ipart=1,nleg_born+1
! This extra statement is needed since in Kaleu the ordering is
! different: 
! initial state: -2:-1
! final state:    1:nfinst
        if (ipart.le.2) then
          subproc(ipart-3) = PDGtoKaleu(flv_ch_Rkin(ipart,iproc))
        else
          subproc(ipart-2) = PDGtoKaleu(flv_ch_Rkin(ipart,iproc))
        end if
      end do
!
!      print *,"iproc: ",iproc
!      print *,"flv: ",subproc
!
      jproc = iproc + num_flv_ch_Bkin
!
! This version of Kaleu allows for multichannel optimization including
! dipole channels too:
      ndip = 0
      if (flg_dipolechan) then
        call GenerateDipoles(flv_ch_Rkin(:,iproc), &
                             num_flv_LO,flv_LO,    &
                             ndip,dip_tmp)
      end if
!
      call multi_kaleu_init(jproc,subproc,nleg_born-1,          &
                            rstot,PDF_option,                   &
                            masses,widths,labels,               &
                            onlyqcd,withqcd,withhiggs,          &
                            nbatchRkin,nstepRkin,thrsRkin,smin, &
                            alphamax_kaleu,alphamax_kaleu,      &
                            alphamax_kaleu,alphamax_kaleu,      &
                            ndip,dip_tmp)
    end do
  end if
!=======================================================================
! Double-real kinematics:
!=======================================================================
  if (flg_NNLO_RR) then
    do iproc=1,num_flv_ch_RRkin
! We have to convert the flavor list into a form which is
! acceptible for Kaleu:
      do ipart=1,nleg_born+2
! This extra statement is needed since in Kaleu the ordering is
! different: 
! initial state: -2:-1
! final state:    1:nfinst
        if (ipart.le.2) then
          subproc(ipart-3) = PDGtoKaleu(flv_ch_RRkin(ipart,iproc))
        else
          subproc(ipart-2) = PDGtoKaleu(flv_ch_RRkin(ipart,iproc))
        end if
      end do
!
!      print *,"iproc: ",iproc
!      print *,"flv: ",subproc
!
      jproc = iproc + num_flv_ch_Bkin + num_flv_ch_Rkin
!
! For the RR part we need dipole and squared dipole channels as well:
      ndip = 0
      if (flg_dipolechan) then
        call GenerateDipoles(flv_ch_RRkin(:,iproc),   &
                             num_flv_NLO_R,flv_NLO_R, &
                             ndip,dip_tmp)
      end if
      ndipsqr = 0
      if (flg_dipsqrchan) then
        call GenerateDipolesSquared(flv_ch_RRkin(:,iproc), &
                                    num_flv_LO,flv_LO,     &
                                    ndipsqr,dipsqr_tmp)
      end if
!
      call multi_kaleu_init(jproc,subproc,nleg_born,               &
                            rstot,PDF_option,                      &
                            masses,widths,labels,                  &
                            onlyqcd,withqcd,withhiggs,             &
                            nbatchRRkin,nstepRRkin,thrsRRkin,smin, &
                            alphamax_kaleu,alphamax_kaleu,         &
                            alphamax_kaleu,alphamax_kaleu,         &
                            ndip,dip_tmp)
    end do
  end if
!
end subroutine init_PSgen
!
subroutine integrate_PSgen
use flags
use statistics
use histo
implicit none
!
!
  integer :: ievent
  character (len=5) :: cont
  integer :: hrs,mins,secs
  real :: start_time,finish_time
!
! For our convenience we log the elapsed time:
  call cpu_time(start_time)
! When the Kaleu generator is used we only specify the number
! of events hence we only have to loop over them:
! In NNLO we can have three different multiplicities:
!=======================================================================
! Born kinematics:
!=======================================================================
  if (flg_LO.or.flg_NLO_V.or.flg_NNLO_VV) then
    cont = 'Bkin '
!    print *,"cont: ",cont
! Optimization mode is active:
    flg_optim = .true.
    ievent = 0
    do while (ievent.lt.neventsBkin)
      ievent = ievent + 1
!      print *,"ievent: ",ievent
      call CalcSigma(cont)
      if (mod(ievent,neventsBkin/niterBkin).eq.0) then
        call finalize_stat(6)
        call cpu_time(finish_time)
        hrs  = int(finish_time - start_time)/3600
        mins = mod(int(finish_time - start_time),3600)/60
        secs = int(finish_time - start_time) - 3600*hrs - 60*mins
        write(*,fmt="(a,I0,a,I0,a,I0,a)") &
          "****************** Elapsed time: ",hrs," hrs, ", &
          mins," mins, ",secs," secs ********************"
      end if
! When we reached the end of optimization phase we nullify the
! statistics:
      if (flg_optim.and.(noptimBkin.eq.ievent)) then
        call null_stati(cont)
! Optimization mode is over:
        flg_optim = .false.
! Nullifying the number of events too:
        ievent = 0
      end if
    end do
! We add the partial result to the output histos:
    if (flg_analysis) call store_hist
  end if
!=======================================================================
! Real kinematics:
!=======================================================================
  if (flg_NLO_R.or.flg_NNLO_RV) then
    cont = 'Rkin '
!    print *,"cont: ",cont
! Optimization mode is active:
    flg_optim = .true.
    ievent = 0
    do while (ievent.lt.neventsRkin)
      ievent = ievent + 1
!      print *,"ievent: ",ievent
      call CalcSigma(cont)
      if (mod(ievent,neventsRkin/niterRkin).eq.0) then
        call finalize_stat(6)
        call cpu_time(finish_time)
        hrs  = int(finish_time - start_time)/3600
        mins = mod(int(finish_time - start_time),3600)/60
        secs = int(finish_time - start_time) - 3600*hrs - 60*mins
        write(*,fmt="(a,I0,a,I0,a,I0,a)") &
          "****************** Elapsed time: ",hrs," hrs, ", &
          mins," mins, ",secs," secs ********************"
      end if
! When we reached the end of optimization phase we nullify the
! statistics:
      if (flg_optim.and.(noptimRkin.eq.ievent)) then
        call null_stati(cont)
! Optimization mode is over:
        flg_optim = .false.
! Nullifying the number of events too:
        ievent = 0
      end if
    end do
! We add the partial result to the output histos:
    if (flg_analysis) call store_hist
  end if
!=======================================================================
! Double-real kinematics:
!=======================================================================
  if (flg_NNLO_RR) then
    cont = 'RRkin'
!    print *,"cont: ",cont
! Optimization mode is active:
    flg_optim = .true.
    ievent = 0
    do while (ievent.lt.neventsRRkin)
      ievent = ievent + 1
!      print *,"ievent: ",ievent
      call CalcSigma(cont)
      if (mod(ievent,neventsRRkin/niterRRkin).eq.0) then
        call finalize_stat(6)
        call cpu_time(finish_time)
        hrs  = int(finish_time - start_time)/3600
        mins = mod(int(finish_time - start_time),3600)/60
        secs = int(finish_time - start_time) - 3600*hrs - 60*mins
        write(*,fmt="(a,I0,a,I0,a,I0,a)") &
          "****************** Elapsed time: ",hrs," hrs, ", &
          mins," mins, ",secs," secs ********************"
      end if
! When we reached the end of optimization phase we nullify the
! statistics:
      if (flg_optim.and.(noptimRRkin.eq.ievent)) then
        call null_stati(cont)
! Optimization mode is over:
        flg_optim = .false.
! Nullifying the number of events too:
        ievent = 0
      end if
    end do
! We add the partial result to the output histos:
    if (flg_analysis) call store_hist
  end if
!
end subroutine integrate_PSgen
!
subroutine gen_PSgen(iproc,cont,discard,p,weight)
use collider
use math
use momenta
use subprocesses
implicit none
!
  integer , intent(in) :: iproc
  character (len=5) , intent(in) :: cont
  logical , intent(out) :: discard
  type(mom) , dimension(:) , intent(out) :: p
  real(kind(1d0)) , intent(out) :: weight
!
  integer :: jproc,ipart,jpart,npart
  integer :: nfinst
  real(kind(1d0)) :: x1,x2
  real(kind(1d0)) :: yij,y1,y2
  real(kind(1d0)) , dimension(4) :: p_tmp
  real(kind(1d0)) , dimension(0:3,-2:17) :: p_kaleu
!
  if (cont.eq.'Bkin ') jproc = iproc
  if (cont.eq.'Rkin ') jproc = num_flv_ch_Bkin + iproc
  if (cont.eq.'RRkin') jproc = num_flv_ch_Bkin + num_flv_ch_Rkin + iproc
!
  call multi_kaleu_gnrt(jproc,discard,x1,x2,p_kaleu)
!
  npart = size(p)
!
! We fill up the momentum array:
  do ipart=1,npart
    if (ipart.le.2) then
! The extra minus sign is needed because in Kaleu momenta are
! considered as out-going:
      p_tmp(4)   = -p_kaleu(0,ipart-3)
      p_tmp(1:3) = -p_kaleu(1:3,ipart-3)
    else
      p_tmp(4)   = p_kaleu(0,ipart-2)
      p_tmp(1:3) = p_kaleu(1:3,ipart-2)
    end if
    p(ipart) = p_tmp
  end do
!
! If discard is true the weight is simply zero:
  if (discard) then
! When discard is true we fill up the momentum array with zeros:
    do ipart=1,npart
      p_tmp = 0d0
      p(ipart) = p_tmp
    end do
! weight is zero:
    weight = 0d0
    return
  end if
!
! Extra technical cut is introduced:
! Beside the NLO-style technical cut an NNLO-style one is
! needed, too: that is the product of any yij pair cannot
! be smaller than ymin:
!  y1 = 1d99
!  y2 = 1d99
  do ipart=3,npart-1
    do jpart=ipart+1,npart
      yij = (p(ipart) + p(jpart))*(p(ipart) + p(jpart))/stot
      yij = abs(yij)
!      if (yij.lt.y1) then
!        y2 = y1
!        y1 = yij
! When the product of the two smallest invariants is
! less then ymin the PS point is rejected:
!        if (y1*y2.lt.ymin) then
!          weight = 0
!          return
!        end if
!      elseif ((yij.gt.y1).and.(yij.lt.y2)) then
!        y2 = yij
!        if (y1*y2.lt.ymin) then
!          weight = 0
!          return
!        end if
!      end if
      if (yij.lt.ymin) then
        weight = 0
        return
      end if
    end do
  end do
!
! Otherwise we obtain the weight too and give it back:
  call multi_kaleu_wght(jproc,weight)
!
! In Kaleu a tower of pis are not included, including them:
  weight = weight * (2d0*pi)**(4 - 3*(npart-2))
!
!  print *,"The momenta are: "
!  call PrintMom(p)
!
end subroutine gen_PSgen
!
subroutine PutWeight_PSgen(iproc,cont,weight)
use subprocesses
implicit none
!
  integer , intent(in) :: iproc
  character (len=5) , intent(in) :: cont
  real(kind(1d0)) , intent(in) :: weight
!
  integer :: jproc
!
  if (cont.eq.'Bkin ') jproc = iproc
  if (cont.eq.'Rkin ') jproc = num_flv_ch_Bkin + iproc
  if (cont.eq.'RRkin') jproc = num_flv_ch_Bkin + num_flv_ch_Rkin + iproc
!
  call multi_kaleu_collect(jproc,abs(weight))
!
end subroutine PutWeight_PSgen
!
end module KaleuPS
!
module RamboPS
implicit none
!
  integer :: PDF_option
  real(kind(1d0)) :: ymin
  real(kind(1d0)) :: smin(-2:17,-2:17)
  real(kind(1d0)) , dimension(:) , allocatable :: masses
  character(20)   , dimension(:) , allocatable :: labels
!
  integer :: niter
  integer :: neventsBkin,neventsRkin,neventsRRkin
!
contains
!
function PDGtoKaleu(PDGid) result(Kid)
implicit none
!
  integer , intent(in) :: PDGid
!
  integer :: Kid
!
! d quark:
  if (abs(PDGid).eq.1) then
    Kid = sign(4,PDGid)
! u quark:
  elseif (abs(PDGid).eq.2) then
    Kid = sign(3,PDGid)
! s quark:
  elseif (abs(PDGid).eq.3) then
    Kid = sign(8,PDGid)
! c quark:
  elseif (abs(PDGid).eq.4) then
    Kid = sign(7,PDGid)
! b quark:
  elseif (abs(PDGid).eq.5) then
    Kid = sign(12,PDGid)
! t quark:
  elseif (abs(PDGid).eq.6) then
    Kid = sign(11,PDGid)
! gluon:
  elseif (abs(PDGid).eq.0) then
    Kid = 17
! electron:
  elseif (abs(PDGid).eq.11) then
    Kid = sign(2,PDGid)
  else
    write(*,'(A,I3)') "No Kaleu alias can be given for: ",PDGid
    stop
  end if
!
end function PDGtoKaleu
!
subroutine init_PSgen()
use process
use collider
use subprocesses
use flags
use input
implicit none
!
!
  integer :: iproc,jproc,ipart
  integer :: numproc
  integer , dimension(-2:17) :: subproc
!
!
  print *,"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
  print *,"@                                                          @"
  print *,"@       We are going to use RAMBO as PS generator...       @"
  print *,"@                                                          @"
  print *,"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
!
! We have to calculate the total number of subprocesses independent
! of the multiplicity to set up the needed arrays for Rambo:
  numproc = num_flv_ch_Bkin + num_flv_ch_Rkin + num_flv_ch_RRkin
!
! We have to allocate those arrays which will hold information for the
! subprocesses in Kaleu:
  call InitRamboArrays(numproc)
!
! We perform process-dependent initialization too:
  call Rambo_initmyprocess
!
  niter = nnloinput("niter")
!
! We need the number of events needed to be generated:
! If only nevents supplied all the contributions are
! calculated with the same number of events:
  neventsBkin   = nnloinput("nevents")
  neventsRkin   = neventsBkin
  neventsRRkin  = neventsBkin
! When neventsBkin, neventsRkin or neventsRRkin turns up in the
! input card we modify the corresponding number of events
! accordingly:
  if (nnloinput("#neventsBkin").gt.0) then
    neventsBkin = nnloinput("#neventsBkin")
  end if
  if (nnloinput("#neventsRkin").gt.0) then
  neventsRkin  = nnloinput("#neventsRkin")
  end if
  if (nnloinput("#neventsRRkin").gt.0) then
  neventsRRkin = nnloinput("#neventsRRkin")
  end if
!
  subproc = 0
!
! We calculate the minimal kinematic invariant from the 
! dimensionless one:
  smin = ymin * stot
!
! We set up all the subprocesses:
!=======================================================================
! Born kinematics:
!=======================================================================
  if (flg_LO.or.flg_NLO_V.or.flg_NNLO_VV) then
    do iproc=1,num_flv_ch_Bkin
! We have to convert the flavor list into a form which is
! acceptible for Kaleu:
      do ipart=1,nleg_born
! This extra statement is needed since in Kaleu the ordering is
! different: 
! initial state: -2:-1
! final state:    1:nfinst
        if (ipart.le.2) then
          subproc(ipart-3) = PDGtoKaleu(flv_ch_Bkin(ipart,iproc))
        else
          subproc(ipart-2) = PDGtoKaleu(flv_ch_Bkin(ipart,iproc))
        end if
      end do
!
!      print *,"iproc: ",iproc
!      print *,"flv: ",subproc
!
      jproc = iproc
!
      call multi_Rambo_init(jproc,subproc,nleg_born-2, &
                            rstot,PDF_option, &
                            masses)
    end do
  end if
!=======================================================================
! Real kinematics:
!=======================================================================
  if (flg_NLO_R.or.flg_NNLO_RV) then
    do iproc=1,num_flv_ch_Rkin
! We have to convert the flavor list into a form which is
! acceptible for Kaleu:
      do ipart=1,nleg_born+1
! This extra statement is needed since in Kaleu the ordering is
! different: 
! initial state: -2:-1
! final state:    1:nfinst
        if (ipart.le.2) then
          subproc(ipart-3) = PDGtoKaleu(flv_ch_Rkin(ipart,iproc))
        else
          subproc(ipart-2) = PDGtoKaleu(flv_ch_Rkin(ipart,iproc))
        end if
      end do
!
!      print *,"iproc: ",iproc
!      print *,"flv: ",subproc
!
      jproc = iproc + num_flv_ch_Bkin
!
      call multi_Rambo_init(jproc,subproc,nleg_born-1, &
                            rstot,PDF_option, &
                            masses)
    end do
  end if
!=======================================================================
! Double-real kinematics:
!=======================================================================
  if (flg_NNLO_RR) then
    do iproc=1,num_flv_ch_RRkin
! We have to convert the flavor list into a form which is
! acceptible for Kaleu:
      do ipart=1,nleg_born+2
! This extra statement is needed since in Kaleu the ordering is
! different: 
! initial state: -2:-1
! final state:    1:nfinst
        if (ipart.le.2) then
          subproc(ipart-3) = PDGtoKaleu(flv_ch_RRkin(ipart,iproc))
        else
          subproc(ipart-2) = PDGtoKaleu(flv_ch_RRkin(ipart,iproc))
        end if
      end do
!
!      print *,"iproc: ",iproc
!      print *,"flv: ",subproc
!
      jproc = iproc + num_flv_ch_Bkin + num_flv_ch_Rkin
!
      call multi_Rambo_init(jproc,subproc,nleg_born, &
                            rstot,PDF_option, &
                            masses)
    end do
  end if
!
end subroutine init_PSgen
!
subroutine integrate_PSgen
use flags
use statistics
use histo
implicit none
!
!
  integer :: ievent
  character (len=5) :: cont
!
! When the Rambo generator is used we only specify the number
! of events hence we only have to loop over them:
! In NNLO we can have three different multiplicities:
!=======================================================================
! Born kinematics:
!=======================================================================
  if (flg_LO.or.flg_NLO_V.or.flg_NNLO_VV) then
    cont = 'Bkin '
! When Rambo is used there is no need to optimize:
    flg_optim = .false.
    do ievent=1,neventsBkin
!      print *,"ievent: ",ievent
      call CalcSigma(cont)
      if (mod(ievent,neventsBkin/niter).eq.0) then
        call finalize_stat(6)
      end if
    end do
! We add the partial result to the output histos:
    if (flg_analysis) call store_hist
  end if
!=======================================================================
! Real kinematics:
!=======================================================================
  if (flg_NLO_R.or.flg_NNLO_RV) then
    cont = 'Rkin '
! When Rambo is used there is no need to optimize:
    flg_optim = .false.
    do ievent=1,neventsRkin
!      print *,"ievent: ",ievent
      call CalcSigma(cont)
      if (mod(ievent,neventsRkin/niter).eq.0) then
        call finalize_stat(6)
      end if
    end do
! We add the partial result to the output histos:
    if (flg_analysis) call store_hist
  end if
!=======================================================================
! Double-real kinematics:
!=======================================================================
  if (flg_NNLO_RR) then
    cont = 'RRkin'
! When Rambo is used there is no need to optimize:
    flg_optim = .false.
    do ievent=1,neventsRRkin
!      print *,"ievent: ",ievent
      call CalcSigma(cont)
      if (mod(ievent,neventsRRkin/niter).eq.0) then
        call finalize_stat(6)
      end if
    end do
! We add the partial result to the output histos:
    if (flg_analysis) call store_hist
  end if
!
end subroutine integrate_PSgen
!
subroutine gen_PSgen(iproc,cont,discard,p,weight)
use math
use momenta
use subprocesses
implicit none
!
  integer , intent(in) :: iproc
  character (len=5) , intent(in) :: cont
  logical , intent(out) :: discard
  type(mom) , dimension(:) , intent(out) :: p
  real(kind(1d0)) , intent(out) :: weight
!
  integer :: jproc,ipart
  integer :: nfinst
  real(kind(1d0)) :: x1,x2
  real(kind(1d0)) , dimension(4) :: p_tmp
  real(kind(1d0)) , dimension(0:3,-2:17) :: p_rambo
!
  if (cont.eq.'Bkin ') jproc = iproc
  if (cont.eq.'Rkin ') jproc = num_flv_ch_Bkin + iproc
  if (cont.eq.'RRkin') jproc = num_flv_ch_Bkin + num_flv_ch_Rkin + iproc
!
  call multi_Rambo_gnrt(jproc,discard,x1,x2,p_rambo)
!
! We fill up the momentum array:
  do ipart=1,size(p)
    if (ipart.le.2) then
! The extra minus sign is needed because in Rambo momenta are
! considered as out-going:
      p_tmp(4)   = -p_rambo(0,ipart-3)
      p_tmp(1:3) = -p_rambo(1:3,ipart-3)
    else
      p_tmp(4)   = p_rambo(0,ipart-2)
      p_tmp(1:3) = p_rambo(1:3,ipart-2)
    end if
    p(ipart) = p_tmp
  end do
!
! If discard is true the weight is simply zero:
  if (discard) then
! When discard is true we fill up the momentum array with zeros:
    do ipart=1,size(p)
      p_tmp = 0d0
      p(ipart) = p_tmp
    end do
! weight is zero:
    weight = 0d0
    return
  end if
!
! Otherwise we obtain the weight too and give it back:
  call multi_Rambo_wght(jproc,weight)
!
! In Rambo a tower of pis are not included, including them:
  weight = weight * (2d0*pi)**(4 - 3*(size(p)-2))
!
!  print *,"The momenta are: "
!  call PrintMom(p)
!
end subroutine gen_PSgen
!
subroutine PutWeight_PSgen(iproc,cont,weight)
use subprocesses
implicit none
!
  integer , intent(in) :: iproc
  character (len=5) , intent(in) :: cont
  real(kind(1d0)) , intent(in) :: weight
!
  integer :: jproc
!
! This routine is empty for Rambo we do not have to give back any
! weight
!
end subroutine PutWeight_PSgen
!
end module RamboPS
!
module phasespace
use KaleuPS
!use RamboPS
use momenta
use particles
implicit none
!
  type(mom) , dimension(:) , allocatable :: ps_B
  type(mom) , dimension(:) , allocatable :: ps_R
  type(mom) , dimension(:) , allocatable :: ps_RR
!
  type(particle) , dimension(:) , allocatable :: parts_B
  type(particle) , dimension(:) , allocatable :: parts_R
  type(particle) , dimension(:) , allocatable :: parts_RR
!
  real(kind(1d0)) , parameter :: gevm22pb = 3.8937966d8
!
  real(kind(1d0)) , dimension(:) , allocatable :: SubTerm_multiscale
!
  real(kind(1d0)) , dimension(:) , allocatable :: dsigmaB,dsigmaBi
  real(kind(1d0)) , dimension(:) , allocatable :: dsigmaV,dsigmaVi
  real(kind(1d0)) , dimension(:) , allocatable :: dsigmaVV,dsigmaVVi
  real(kind(1d0)) , dimension(:) , allocatable :: dsigmaR,dsigmaRi
  real(kind(1d0)) , dimension(:) , allocatable :: dsigmaRV,dsigmaRVi
  real(kind(1d0)) , dimension(:) , allocatable :: dsigmaRR,dsigmaRRi
  real(kind(1d0)) , dimension(:) , allocatable :: dsigmaR_A1
  real(kind(1d0)) , dimension(:) , allocatable :: dsigmaRV_A1
  real(kind(1d0)) , dimension(:) , allocatable :: dsigmaRR_A1
  real(kind(1d0)) , dimension(:) , allocatable :: dsigmaB_I1
  real(kind(1d0)) , dimension(:) , allocatable :: dsigma_Innlo
  real(kind(1d0)) , dimension(:) , allocatable :: dsigmaR_I1
!
  real(kind(1d0)) , dimension(:) , allocatable :: subs_R
  real(kind(1d0)) , dimension(:) , allocatable :: subs_RV
!FIXME this variable will not be used after finishing the 
! Calc_mp2_subs routine...
  real(kind(1d0)) , dimension(:) , allocatable :: subs_RRA1
  real(kind(1d0)) , dimension(:) , allocatable :: subs_RR
!
  real(kind(1d0)) , dimension(:,:) , allocatable :: Bij_arr
  real(kind(1d0)) , dimension(:,:,:) , allocatable :: BijLaurent_arr
  real(kind(1d0)) , dimension(:,:,:) , allocatable :: Bijk_arr
  real(kind(1d0)) , dimension(:,:,:,:) , allocatable :: Bijkl_arr
  real(kind(1d0)) , dimension(:,:,:,:,:) , allocatable :: BijklLaurent_arr
  real(kind(1d0)) , dimension(:,:,:,:) , allocatable :: Bmunuij_arr
  real(kind(1d0)) , dimension(:,:) , allocatable :: Vij_arr
  real(kind(1d0)) , dimension(:,:,:) , allocatable :: VijLaurent_arr
  real(kind(1d0)) , dimension(:,:) , allocatable :: Rij_arr
!
contains
!
subroutine init_PS()
use process
use momenta
use scales
implicit none
!
!
  integer :: istat
!
! We allocate the arrays which will hold the momenta:
! Born-like kinematics:
  allocate(ps_B(nleg_born),stat=istat)
  if (istat.ne.0) then
    print *,"Problem with the allocation of ps_B..."
    stop
  end if
! Real-like kinematics:
  allocate(ps_R(nleg_born+1),stat=istat)
  if (istat.ne.0) then
    print *,"Problem with the allocation of ps_R..."
    stop
  end if
! Double-real-like kinematics:
  allocate(ps_RR(nleg_born+2),stat=istat)
  if (istat.ne.0) then
    print *,"Problem with the allocation of ps_RR..."
    stop
  end if
!
! We also allocate arrays which will hold momenta and flavours:
! Born-like kinematics:
  allocate(parts_B(nleg_born),stat=istat)
  if (istat.ne.0) then
    print *,"Problem with the allocation of parts_B..."
    stop
  end if
! Real-like kinematics:
  allocate(parts_R(nleg_born+1),stat=istat)
  if (istat.ne.0) then
    print *,"Problem with the allocation of parts_R..."
    stop
  end if
! Double-real-like kinematics:
  allocate(parts_RR(nleg_born+2),stat=istat)
  if (istat.ne.0) then
    print *,"Problem with the allocation of parts_RR..."
    stop
  end if
!
! We have to initialize the generator:
  call init_PSgen()
!
! We also allocate an array which will hold a subtraction term
! at different scale choices:
  allocate(SubTerm_multiscale(nscales),stat=istat)
  if (istat.ne.0) then
    print *,"Problem occured during allocating SubTerm_multiscale..."
    stop
  end if
! Allocation of array holding combined subtraction terms:
  allocate(subs_R(nscales),    &
           subs_RV(nscales),   &
           subs_RR(nscales),   &
           subs_RRA1(nscales), &
           stat=istat)
  if (istat.ne.0) then
    print *,"Problem with allocation of Subtractions..."
    stop "init_PS"
  end if
  subs_R    = 0d0
  subs_RV   = 0d0
  subs_RR   = 0d0
  subs_RRA1 = 0d0
! Allocation of the differential Born-like contributions:
  allocate(dsigmaB(nscales),    &
           dsigmaBi(nscales),   &
           dsigmaV(nscales),    &
           dsigmaVi(nscales),   &
           dsigmaVV(nscales),   &
           dsigmaVVi(nscales),  &
           stat=istat)
  if (istat.ne.0) then
    print *,"Problem occured during allocating dsigmaB..."
    stop
  end if
  dsigmaB   = 0d0
  dsigmaV   = 0d0
  dsigmaVV  = 0d0
  dsigmaBi  = 0d0
  dsigmaVi  = 0d0
  dsigmaVVi = 0d0
! Allocation of the differential Real-like contributions:
  allocate(dsigmaR(nscales), &
           dsigmaRi(nscales), &
           dsigmaRV(nscales), &
           dsigmaRVi(nscales), &
           stat=istat)
  if (istat.ne.0) then
    print *,"Problem occured during allocating dsigmaR..."
    stop
  end if
  dsigmaR  = 0d0
  dsigmaRV = 0d0
! Allocation of the I operators:
  allocate(dsigmaB_I1(nscales),   &
           dsigma_Innlo(nscales), &
           dsigmaR_I1(nscales),   &
           stat=istat)
  if (istat.ne.0) then
    print *,"Problem occured during allocating I op arrays..."
    stop
  end if
  dsigmaB_I1    = 0d0
  dsigma_Innlo  = 0d0
  dsigmaR_I1    = 0d0
! Allocation of the Bij arrays:
  allocate(Bij_arr(nleg_born,nleg_born),                                   &
           BijLaurent_arr(nleg_born,nleg_born,-4:2),                       &
           Bijk_arr(nleg_born,nleg_born,nleg_born),                        &
           Bijkl_arr(nleg_born,nleg_born,nleg_born,nleg_born),             &
           BijklLaurent_arr(nleg_born,nleg_born,nleg_born,nleg_born,-4:2), &
           Bmunuij_arr(0:3,0:3,nleg_born,nleg_born),                       &
           Vij_arr(nleg_born,nleg_born),                                   &
           VijLaurent_arr(nleg_born,nleg_born,-4:2),                       &
           Rij_arr(nleg_born+1,nleg_born+1),                               &
           stat=istat)
  if (istat.ne.0) then
    print *,"Problem occured during allocating color-cor. arrays..."
    stop
  end if
  Bij_arr          = 0d0
  BijLaurent_arr   = 0d0
  Bijk_arr         = 0d0
  Bijkl_arr        = 0d0
  BijklLaurent_arr = 0d0
  Bmunuij_arr      = 0d0
  Vij_arr          = 0d0
  VijLaurent_arr   = 0d0
  Rij_arr          = 0d0
! Allocation of the A1 subtractions for the Real, Real-Virtual and
! double-Real contributions:
  allocate(dsigmaR_A1(nscales), &
           dsigmaRV_A1(nscales), &
           dsigmaRR_A1(nscales), &
           stat=istat)
  if (istat.ne.0) then
    print *,"Problem occured during allocating A1 terms..."
    stop
  end if
  dsigmaR_A1  = 0d0
  dsigmaRV_A1 = 0d0
  dsigmaRR_A1 = 0d0
! Allocation of the differential Double-Real contributions:
  allocate(dsigmaRR(nscales), &
           dsigmaRRi(nscales), &
           stat=istat)
  if (istat.ne.0) then
    print *,"Problem occured during allocating dsigmaRR..."
    stop
  end if
  dsigmaRR  = 0d0
  dsigmaRRi = 0d0
!
end subroutine init_PS
!
subroutine integrate
implicit none
!
!
! We call the PS generator and it will decide how to integrate:
  call integrate_PSgen()
!
end subroutine integrate
!
subroutine GetPSpoint(iproc,cont,discard,weight)
use momenta
use collider
implicit none
!
  integer , intent(in) :: iproc
  character (len=5) , intent(in) :: cont
  logical , intent(out) :: discard
  real(kind(1d0)) , intent(out) :: weight
!
!
  if (cont.eq.'Bkin ') then
    call gen_PSgen(iproc,cont,discard,ps_B,weight)
  elseif (cont.eq.'Rkin ') then
    call gen_PSgen(iproc,cont,discard,ps_R,weight)
  elseif (cont.eq.'RRkin') then
    call gen_PSgen(iproc,cont,discard,ps_RR,weight)
  end if
!
! We include the conversion factor from GeV^{-2} to pbs:
  weight = weight * gevm22pb
!
end subroutine GetPSpoint
!
subroutine PutWeight(iproc,cont,weight)
implicit none
!
  integer , intent(in) :: iproc
  character (len=5) , intent(in) :: cont
  real(kind(1d0)) , intent(in) :: weight
!
!
  call PutWeight_PSgen(iproc,cont,weight)
!
end subroutine PutWeight
!
end module phasespace
