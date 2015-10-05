module regions
use particles
implicit none
!
  type parton
    integer :: i
    integer :: ID
  end type parton
!
  type subtractterm
    integer , dimension(:) , allocatable :: uborn
    type(parton) , dimension(:) , allocatable :: emit,rad
    real(kind(1d0)) :: symfact
  end type subtractterm
!
  type subterms
    integer :: numterm
    type(subtractterm) , dimension(:) , allocatable :: term
  end type subterms
!
  type(particle) , dimension(:) , allocatable :: psubB
  type(particle) , dimension(:) , allocatable :: psubR
  type(particle) , dimension(:) , allocatable :: psubRR
!
  real(kind(1d0)) , dimension(:,:) , allocatable :: subBij
  real(kind(1d0)) , dimension(:,:) , allocatable :: subRij
!
  type(subterms) , dimension(:) , allocatable :: subterms_Cir_R
  type(subterms) , dimension(:) , allocatable :: subterms_Sr_R
  type(subterms) , dimension(:) , allocatable :: subterms_CSir_R
!
  type(subterms) , dimension(:) , allocatable :: subterms_Cir_RR
  type(subterms) , dimension(:) , allocatable :: subterms_Sr_RR
  type(subterms) , dimension(:) , allocatable :: subterms_CSir_RR
!
  type(subterms) , dimension(:) , allocatable :: subterms_Cirs_RR
  type(subterms) , dimension(:) , allocatable :: subterms_Cirjs_RR
  type(subterms) , dimension(:) , allocatable :: subterms_CSirs_RR
  type(subterms) , dimension(:) , allocatable :: subterms_Srs_RR
!
! Declaration of kinematics related quantities:
  integer :: sub_d0,sub_d0pr
!
  type(mom) :: Q
  real(kind(1d0)) :: Q2
  real(kind(1d0)) , dimension(:) , allocatable :: yiQ_arr
  real(kind(1d0)) , dimension(:,:) , allocatable :: sir_arr
  real(kind(1d0)) , dimension(:,:) , allocatable :: yir_arr
!
  real(kind(1d0)) , dimension(:) , allocatable :: yihatQ_arr
  real(kind(1d0)) , dimension(:,:) , allocatable :: sirhat_arr
!
  real(kind(1d0)) , dimension(:) , allocatable :: yitildeQ_arr
  real(kind(1d0)) , dimension(:,:) , allocatable :: sirtilde_arr
  real(kind(1d0)) , dimension(:,:) , allocatable :: yirtilde_arr
!
  real(kind(1d0)) :: yiQ,yrQ,yjQ,ysQ,yirQ,yjsQ,yrsQ,yirsQ,yirtildeQ
  real(kind(1d0)) :: sir,yir,sjs,yjs,yrs,sirs,yirs
  real(kind(1d0)) :: s_k_t,s_kt_r
  real(kind(1d0)) :: ztildei,ztilder,ztildej,ztildes, &
                     ztildek,ztildet,zhati,zhatr,zhatt,zhatkt
  real(kind(1d0)) :: alphair,alphajs,alphairs,alphairt,alphakt,alphaktr
  real(kind(1d0)) :: lambda_r,lambda_t,lambdars
  real(kind(1d0)) :: zetai_r,zetar_i,zetair
  real(kind(1d0)) :: zetaj_s,zetas_j,zetajs

  type(particle) :: kttildei,kttilder,kttildej,kttildes, &
                    kttildek,kttildet, &
                    kthati,kthatkt,kthatr, &
                    kthatir,kthatt
  type(particle) :: ptildeir,ptildejs
  type(particle) :: ptildeirs
!
! logicals for initializations:
  logical :: init_CalcSr10 = .true.
!
contains
!
subroutine PrintSubTerms
use flags
use subprocesses
use particles
use utils
implicit none
!
!
  integer , parameter :: iun = 99
  integer :: iproc,nproc,iterm
  character (len=128) :: fname
!
!
! We open a file and we dump all the subtraction terms there:
  fname = 'subtractionterms.txt'
  open(file=fname,status='unknown',unit=iun)
  write(iun,'(a)') &
  "********************************************************************"
  write(iun,'(a)') &
  "************************ Subtraction Terms *************************"
  write(iun,'(a)') &
  "********************************************************************"
  write(iun,*) 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! NLO R - A1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  if (flg_NLO_R_A1) then
    write(iun,'(a)') &
  "####################################################################"
    write(iun,'(a)') &
  "############################# NLO R ################################"
    write(iun,'(a)') &
  "####################################################################"
    do iproc=1,num_flv_irr_NLO_R
      write(iun,fmt='(a,I0,a,F0.3,a)',advance='no') "iproc: ",iproc, &
      " , weight: ",weights_Rkin(iproc)," , "
      call PrintSubproc(iun,flv_ch_Rkin(:,iproc))
      if (subterms_Cir_R(iproc)%numterm.ne.0) then
        write(iun,fmt='(a)') "~~~~ Cir ~~~~"
      end if
      do iterm=1,subterms_Cir_R(iproc)%numterm
        write(iun,fmt='(a,I4,3(3(a),I0),a)') &
        "iterm: ",iterm," , ", &
        ConvertFromPDG(subterms_Cir_R(iproc)%term(iterm)%emit(1)%ID), &
        "(",subterms_Cir_R(iproc)%term(iterm)%emit(1)%i,") -> ", &
        ConvertFromPDG(subterms_Cir_R(iproc)%term(iterm)%rad(1)%ID), &
        "(",subterms_Cir_R(iproc)%term(iterm)%rad(1)%i,") || ", &
        ConvertFromPDG(subterms_Cir_R(iproc)%term(iterm)%rad(2)%ID), &
        "(",subterms_Cir_R(iproc)%term(iterm)%rad(2)%i,")" 
        call PrintBranching(iun,subterms_Cir_R(iproc)%term(iterm))
      end do
      if (subterms_Sr_R(iproc)%numterm.ne.0) then
        write(iun,fmt='(a)') "~~~~ Sr ~~~~"
      end if
      do iterm=1,subterms_Sr_R(iproc)%numterm
        write(iun,fmt='(a,I4,3(a),I0,a)') &
        "iterm: ",iterm," , ", &
        ConvertFromPDG(subterms_Sr_R(iproc)%term(iterm)%rad(1)%ID), &
        "(",subterms_Sr_R(iproc)%term(iterm)%rad(1)%i,") -> 0" 
! We print the Underlying Born:
        write(iun,fmt='(a)',advance='no') "UBorn: "
        call PrintSubproc(iun, &
          subterms_Sr_R(iproc)%term(iterm)%uborn(:))
      end do
      if (subterms_CSir_R(iproc)%numterm.ne.0) then
        write(iun,fmt='(a)') "~~~~ CSir ~~~~"
      end if
      do iterm=1,subterms_CSir_R(iproc)%numterm
        write(iun,fmt='(a,I4,4(3(a),I0),a)') &
        "iterm: ",iterm," , ", &
        ConvertFromPDG(subterms_CSir_R(iproc)%term(iterm)%emit(1)%ID), &
        "(",subterms_CSir_R(iproc)%term(iterm)%emit(1)%i,") -> ", &
        ConvertFromPDG(subterms_CSir_R(iproc)%term(iterm)%rad(1)%ID), &
        "(",subterms_CSir_R(iproc)%term(iterm)%rad(1)%i,") || ", &
        ConvertFromPDG(subterms_CSir_R(iproc)%term(iterm)%rad(2)%ID), &
        "(",subterms_CSir_R(iproc)%term(iterm)%rad(2)%i,") , ", & 
        ConvertFromPDG(subterms_CSir_R(iproc)%term(iterm)%rad(2)%ID), &
        "(",subterms_CSir_R(iproc)%term(iterm)%rad(2)%i,") -> 0" 
        call PrintBranching(iun,subterms_CSir_R(iproc)%term(iterm))
      end do
    end do
  end if
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! NNLO RR - A1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  if (flg_NNLO_RR_A1) then
    write(iun,'(a)') &
  "####################################################################"
    write(iun,'(a)') &
  "############################ NNLO RR - A1 ##########################"
    write(iun,'(a)') &
  "####################################################################"
    do iproc=1,num_flv_irr_NNLO_RR
      write(iun,fmt='(a,I0,a,F0.3,a)',advance='no') "iproc: ",iproc, &
      " , weight: ",weights_RRkin(iproc)," , "
      call PrintSubproc(iun,flv_ch_RRkin(:,iproc))
      if (subterms_Cir_RR(iproc)%numterm.ne.0) then
        write(iun,fmt='(a)') "~~~~ Cir ~~~~"
      end if
      do iterm=1,subterms_Cir_RR(iproc)%numterm
        write(iun,fmt='(a,I4,3(3(a),I0),a)') &
        "iterm: ",iterm," , ", &
        ConvertFromPDG(subterms_Cir_RR(iproc)%term(iterm)%emit(1)%ID), &
        "(",subterms_Cir_RR(iproc)%term(iterm)%emit(1)%i,") -> ", &
        ConvertFromPDG(subterms_Cir_RR(iproc)%term(iterm)%rad(1)%ID), &
        "(",subterms_Cir_RR(iproc)%term(iterm)%rad(1)%i,") || ", &
        ConvertFromPDG(subterms_Cir_RR(iproc)%term(iterm)%rad(2)%ID), &
        "(",subterms_Cir_RR(iproc)%term(iterm)%rad(2)%i,")" 
        call PrintBranching(iun,subterms_Cir_RR(iproc)%term(iterm))
      end do
      if (subterms_Sr_RR(iproc)%numterm.ne.0) then
        write(iun,fmt='(a)') "~~~~ Sr ~~~~"
      end if
      do iterm=1,subterms_Sr_RR(iproc)%numterm
        write(iun,fmt='(a,I4,3(a),I0,a)') &
        "iterm: ",iterm," , ", &
        ConvertFromPDG(subterms_Sr_RR(iproc)%term(iterm)%rad(1)%ID), &
        "(",subterms_Sr_RR(iproc)%term(iterm)%rad(1)%i,") -> 0" 
! We print the Underlying Born:
        write(iun,fmt='(a)',advance='no') "UBorn: "
        call PrintSubproc(iun, &
          subterms_Sr_RR(iproc)%term(iterm)%uborn(:))
      end do
      if (subterms_CSir_RR(iproc)%numterm.ne.0) then
        write(iun,fmt='(a)') "~~~~ CSir ~~~~"
      end if
      do iterm=1,subterms_CSir_RR(iproc)%numterm
        write(iun,fmt='(a,I4,4(3(a),I0),a)') &
        "iterm: ",iterm," , ", &
        ConvertFromPDG(subterms_CSir_RR(iproc)%term(iterm)%emit(1)%ID), &
        "(",subterms_CSir_RR(iproc)%term(iterm)%emit(1)%i,") -> ", &
        ConvertFromPDG(subterms_CSir_RR(iproc)%term(iterm)%rad(1)%ID), &
        "(",subterms_CSir_RR(iproc)%term(iterm)%rad(1)%i,") || ", &
        ConvertFromPDG(subterms_CSir_RR(iproc)%term(iterm)%rad(2)%ID), &
        "(",subterms_CSir_RR(iproc)%term(iterm)%rad(2)%i,") , ", & 
        ConvertFromPDG(subterms_CSir_RR(iproc)%term(iterm)%rad(2)%ID), &
        "(",subterms_CSir_RR(iproc)%term(iterm)%rad(2)%i,") -> 0" 
        call PrintBranching(iun,subterms_CSir_RR(iproc)%term(iterm))
      end do
    end do
  end if
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! NNLO RV - A1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  if (flg_NNLO_RV_A1) then
    nproc = num_flv_irr_NLO_R
    if (.not.flg_NLO_R) nproc = num_flv_irr_NNLO_RV
    write(iun,'(a)') &
  "####################################################################"
    write(iun,'(a)') &
  "############################ NNLO RV - A1 ##########################"
    write(iun,'(a)') &
  "####################################################################"
    do iproc=1,nproc
      write(iun,fmt='(a,I0,a,F0.3,a)',advance='no') "iproc: ",iproc, &
      " , weight: ",weights_Rkin(iproc)," , "
      call PrintSubproc(iun,flv_ch_Rkin(:,iproc))
      if (subterms_Cir_R(iproc)%numterm.ne.0) then
        write(iun,fmt='(a)') "~~~~ Cir ~~~~"
      end if
      do iterm=1,subterms_Cir_R(iproc)%numterm
        write(iun,fmt='(a,I4,3(3(a),I0),a)') &
        "iterm: ",iterm," , ", &
        ConvertFromPDG(subterms_Cir_R(iproc)%term(iterm)%emit(1)%ID), &
        "(",subterms_Cir_R(iproc)%term(iterm)%emit(1)%i,") -> ", &
        ConvertFromPDG(subterms_Cir_R(iproc)%term(iterm)%rad(1)%ID), &
        "(",subterms_Cir_R(iproc)%term(iterm)%rad(1)%i,") || ", &
        ConvertFromPDG(subterms_Cir_R(iproc)%term(iterm)%rad(2)%ID), &
        "(",subterms_Cir_R(iproc)%term(iterm)%rad(2)%i,")" 
        call PrintBranching(iun,subterms_Cir_R(iproc)%term(iterm))
      end do
      if (subterms_Sr_R(iproc)%numterm.ne.0) then
        write(iun,fmt='(a)') "~~~~ Sr ~~~~"
      end if
      do iterm=1,subterms_Sr_R(iproc)%numterm
        write(iun,fmt='(a,I4,3(a),I0,a)') &
        "iterm: ",iterm," , ", &
        ConvertFromPDG(subterms_Sr_R(iproc)%term(iterm)%rad(1)%ID), &
        "(",subterms_Sr_R(iproc)%term(iterm)%rad(1)%i,") -> 0" 
! We print the Underlying Born:
        write(iun,fmt='(a)',advance='no') "UBorn: "
        call PrintSubproc(iun, &
          subterms_Sr_R(iproc)%term(iterm)%uborn(:))
      end do
      if (subterms_CSir_R(iproc)%numterm.ne.0) then
        write(iun,fmt='(a)') "~~~~ CSir ~~~~"
      end if
      do iterm=1,subterms_CSir_R(iproc)%numterm
        write(iun,fmt='(a,I4,4(3(a),I0),a)') &
        "iterm: ",iterm," , ", &
        ConvertFromPDG(subterms_CSir_R(iproc)%term(iterm)%emit(1)%ID), &
        "(",subterms_CSir_R(iproc)%term(iterm)%emit(1)%i,") -> ", &
        ConvertFromPDG(subterms_CSir_R(iproc)%term(iterm)%rad(1)%ID), &
        "(",subterms_CSir_R(iproc)%term(iterm)%rad(1)%i,") || ", &
        ConvertFromPDG(subterms_CSir_R(iproc)%term(iterm)%rad(2)%ID), &
        "(",subterms_CSir_R(iproc)%term(iterm)%rad(2)%i,") , ", & 
        ConvertFromPDG(subterms_CSir_R(iproc)%term(iterm)%rad(2)%ID), &
        "(",subterms_CSir_R(iproc)%term(iterm)%rad(2)%i,") -> 0" 
        call PrintBranching(iun,subterms_CSir_R(iproc)%term(iterm))
      end do
    end do
  end if
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! NNLO RR - A2
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  if (flg_NNLO_RR_A2) then
    write(iun,'(a)') &
  "####################################################################"
    write(iun,'(a)') &
  "############################ NNLO RR - A2 ##########################"
    write(iun,'(a)') &
  "####################################################################"
    do iproc=1,num_flv_irr_NNLO_RR
      write(iun,fmt='(a,I0,a,F0.3,a)',advance='no') "iproc: ",iproc, &
      " , weight: ",weights_RRkin(iproc)," , "
      call PrintSubproc(iun,flv_ch_RRkin(:,iproc))
      if (subterms_Cirs_RR(iproc)%numterm.ne.0) then
        write(iun,fmt='(a)') "~~~~ Cirs ~~~~"
      end if
      do iterm=1,subterms_Cirs_RR(iproc)%numterm
        write(iun,fmt='(a,I4,4(3(a),I0),a)') &
        "iterm: ",iterm," , ", &
        ConvertFromPDG(subterms_Cirs_RR(iproc)%term(iterm)%emit(1)%ID), &
        "(",subterms_Cirs_RR(iproc)%term(iterm)%emit(1)%i,") -> ", &
        ConvertFromPDG(subterms_Cirs_RR(iproc)%term(iterm)%rad(1)%ID), &
        "(",subterms_Cirs_RR(iproc)%term(iterm)%rad(1)%i,") || ", &
        ConvertFromPDG(subterms_Cirs_RR(iproc)%term(iterm)%rad(2)%ID), &
        "(",subterms_Cirs_RR(iproc)%term(iterm)%rad(2)%i,") || ", & 
        ConvertFromPDG(subterms_Cirs_RR(iproc)%term(iterm)%rad(3)%ID), &
        "(",subterms_Cirs_RR(iproc)%term(iterm)%rad(3)%i,")" 
        call PrintBranching(iun,subterms_Cirs_RR(iproc)%term(iterm))
      end do
      if (subterms_Cirjs_RR(iproc)%numterm.ne.0) then
        write(iun,fmt='(a)') "~~~~ Cirjs ~~~~"
      end if
      do iterm=1,subterms_Cirjs_RR(iproc)%numterm
        write(iun,fmt='(a,I4,6(3(a),I0),a)') &
        "iterm: ",iterm," , ", &
        ConvertFromPDG(subterms_Cirjs_RR(iproc)%term(iterm)%emit(1)%ID), &
        "(",subterms_Cirjs_RR(iproc)%term(iterm)%emit(1)%i,") -> ", &
        ConvertFromPDG(subterms_Cirjs_RR(iproc)%term(iterm)%rad(1)%ID), &
        "(",subterms_Cirjs_RR(iproc)%term(iterm)%rad(1)%i,") || ", &
        ConvertFromPDG(subterms_Cirjs_RR(iproc)%term(iterm)%rad(2)%ID), &
        "(",subterms_Cirjs_RR(iproc)%term(iterm)%rad(2)%i,")  ,  ", & 
        ConvertFromPDG(subterms_Cirjs_RR(iproc)%term(iterm)%emit(2)%ID), &
        "(",subterms_Cirjs_RR(iproc)%term(iterm)%emit(2)%i,") -> ", &
        ConvertFromPDG(subterms_Cirjs_RR(iproc)%term(iterm)%rad(3)%ID), &
        "(",subterms_Cirjs_RR(iproc)%term(iterm)%rad(3)%i,") || ", & 
        ConvertFromPDG(subterms_Cirjs_RR(iproc)%term(iterm)%rad(4)%ID), &
        "(",subterms_Cirjs_RR(iproc)%term(iterm)%rad(4)%i,")" 
        call PrintBranching(iun,subterms_Cirjs_RR(iproc)%term(iterm))
      end do
!
      if (subterms_CSirs_RR(iproc)%numterm.ne.0) then
        write(iun,fmt='(a)') "~~~~ CSirs ~~~~"
      end if
      do iterm=1,subterms_CSirs_RR(iproc)%numterm
        write(iun,fmt='(a,I4,6(3(a),I0),a)') &
        "iterm: ",iterm," , ", &
        ConvertFromPDG(subterms_CSirs_RR(iproc)%term(iterm)%emit(1)%ID), &
        "(",subterms_CSirs_RR(iproc)%term(iterm)%emit(1)%i,") -> ", &
        ConvertFromPDG(subterms_CSirs_RR(iproc)%term(iterm)%rad(1)%ID), &
        "(",subterms_CSirs_RR(iproc)%term(iterm)%rad(1)%i,") || ", &
        ConvertFromPDG(subterms_CSirs_RR(iproc)%term(iterm)%rad(2)%ID), &
        "(",subterms_CSirs_RR(iproc)%term(iterm)%rad(2)%i,")  ,  ", & 
        ConvertFromPDG(subterms_CSirs_RR(iproc)%term(iterm)%rad(3)%ID), &
        "(",subterms_CSirs_RR(iproc)%term(iterm)%rad(3)%i,") -> 0"
        call PrintBranching(iun,subterms_CSirs_RR(iproc)%term(iterm))
      end do
!
      if (subterms_Srs_RR(iproc)%numterm.ne.0) then
        write(iun,fmt='(a)') "~~~~ Srs ~~~~"
      end if
      do iterm=1,subterms_Srs_RR(iproc)%numterm
        write(iun,fmt='(a,I4,4(3(a),I0),a)') &
        "iterm: ",iterm," , ", &
        ConvertFromPDG(subterms_Srs_RR(iproc)%term(iterm)%rad(1)%ID), &
        "(",subterms_Srs_RR(iproc)%term(iterm)%rad(1)%i,") -> 0  ,  ", &
        ConvertFromPDG(subterms_Srs_RR(iproc)%term(iterm)%rad(2)%ID), &
        "(",subterms_Srs_RR(iproc)%term(iterm)%rad(2)%i,") -> 0"
        call PrintBranching(iun,subterms_Srs_RR(iproc)%term(iterm))
      end do
    end do
  end if
!
  close(iun)
!
end subroutine PrintSubTerms
!
subroutine PrintBranching(iun,subterm)
use utils
use particles
implicit none
!
  integer , intent(in) :: iun
  type(subtractterm) , intent(in) :: subterm
!
  integer :: i,ibranch,jbranch,kbranch,ich,iemit,jemit
!
! We print the Underlying Born:
  write(iun,fmt='(a)',advance='no') "UBorn: "
  call PrintSubproc(iun,subterm%uborn(:))
!
! TODO the initial state is not treated properly!
! We need to write as many lines as many branchings occur:
  do ibranch=size(subterm%emit(:)),1,-1
    iemit = subterm%emit(ibranch)%i
! When multiple branchings take place the position of the emitter has to be
! decreased by one unit if the position for the previous radiated parton comes
! before the current emitter:
! Became obsolete:
!    do jbranch=1,ibranch-1
!      if (max(subterm%rad(2*(jbranch-1)+1)%i, &
!              subterm%rad(2*(jbranch-1)+2)%i).lt.iemit) then
!        iemit = iemit - 1
!      end if
!    end do
! Initial State:
    if (iemit.lt.3) then
! Final State:
    else
! From the number of branchings and the position of the emitter the
! number of spaces mandatory included can be calculated:
      do ich=1,16 + 3*(iemit - 3) - (ibranch-1)
        write(iun,fmt='(a)',advance='no') ' '
! For all those branchings which have a smaller position a vertical
! bar has to be included, but care should be taken that the position
! of the emitter can change due to other branchings:
        do jbranch=1,ibranch-1
          jemit = subterm%emit(jbranch)%i
! Became obsolete:
!          do kbranch=1,jbranch-1
!            if (max(subterm%rad(2*(kbranch-1)+1)%i, &
!                    subterm%rad(2*(kbranch-1)+2)%i).lt.jemit) then
!              jemit = jemit - 1
!            end if
!          end do
          if (ich.eq.16 + 3*(jemit - 3)) then
            write(iun,fmt='(a)',advance='no') '|'
          end if
        end do
      end do
      write(iun,fmt='(a)',advance='no') '\-> '
! If only branching is possible just dump the radiated partons:
      if (size(subterm%emit(:)).eq.1) then
        write(iun,fmt='(3(1x,a))') &
          (ConvertFromPDG(subterm%rad(i)%ID),i=1, &
          size(subterm%rad(:)))
      else
        write(iun,fmt='(2(1x,a))') &
          (ConvertFromPDG(subterm%rad(i)%ID), i=2*ibranch-1,2*ibranch)
      end if
    end if
  end do
!
end subroutine PrintBranching
!
subroutine PrintBranching_old(iun,subterm)
use utils
use particles
implicit none
!
  integer , intent(in) :: iun
  type(subtractterm) , intent(in) :: subterm
!
  integer :: i,ibranch,jbranch,kbranch,ich,iemit,irad
!
! We print the Underlying Born:
  write(iun,fmt='(a)',advance='no') "UBorn: "
  call PrintSubproc(iun,subterm%uborn(:))
!
! TODO the initial state is not treated properly!
  do ibranch=size(subterm%emit(:)),1,-1
    ich = 1
    do while (.true.)
      write(iun,fmt='(a)',advance='no') ' '
      ich = ich + 1
      do jbranch=1,1!size(subterm%emit(:))
        iemit = subterm%emit(jbranch)%i
        irad  = subterm%rad(2*(jbranch-1)+1)%i
! When multiple emissions are present and the radiated parton has
! a position less than the upcoming emitters the position of the
! up-coming emitters have to be decreased by one:
        do kbranch=1,jbranch-1
          if (subterm%rad(2*(kbranch-1)+1)%i.lt.iemit) iemit = iemit - 1
        end do
        if (min(iemit,irad).lt.3) then
          if (ich.eq.(7 + 3*(min(iemit,irad) - 1))) then
            if (ibranch.ne.jbranch) then
              write(iun,fmt='(a)',advance='no') '|'
            else
              write(iun,fmt='(a)',advance='no') char(192)
            end if
            if (jbranch.eq.size(subterm%emit(:))) then
              goto 100
            end if
            ich = ich + 1
          end if
        else
          if (ich.eq.(16 + 3*(min(iemit,irad) - 3) + 1)) then
            if (ibranch.ne.jbranch) then
              write(iun,fmt='(a)',advance='no') '|'
            else
              write(iun,fmt='(a)',advance='no') '\-> '
              if (size(subterm%emit(:)).eq.1) then
                write(iun,fmt='(3(1x,a))') &
                  (ConvertFromPDG(subterm%rad(i)%ID),i=1, &
                  size(subterm%rad(:)))
              else
                write(iun,fmt='(1x,a)') &
                  ConvertFromPDG(subterm%rad(1)%ID)
              end if
            end if
            if (jbranch.eq.size(subterm%emit(:))) then
              write(iun,*)
              goto 100
            end if
            ich = ich + 1
          end if
        end if
      end do
    end do
100 continue
  end do
!
end subroutine PrintBranching_old
!
! ileg and jleg should be ordered: ileg < jleg
subroutine UndoBranching(ileg,jleg,flv_in,flv_out)
implicit none
!
  integer , intent(in) :: ileg,jleg
  integer , dimension(:) , intent(in) :: flv_in
  integer , dimension(:) , intent(out) :: flv_out
!
  integer :: ipart,jpart
  integer :: nleg,merged
!
!
! We obtain the parton type before branching:
! ISR:
  if (ileg.lt.3) then
    merged = flv_in(ileg) - flv_in(jleg)
! FSR:
  else
    merged = flv_in(ileg) + flv_in(jleg)
  end if
!
  nleg = size(flv_in)
  jpart = 0
!
! If ileg equals jleg we simply omit that particular leg:
  if (ileg.eq.jleg) then
    do ipart=1,nleg
      if (ipart.ne.ileg) then
        jpart = jpart + 1
        flv_out(jpart) = flv_in(ipart)
      end if
    end do
    return
  end if
!
  do ipart=1,nleg
! Non-affected leg:
    if ((ipart.ne.ileg).and.(ipart.ne.jleg)) then
      jpart = jpart + 1
      flv_out(jpart) = flv_in(ipart)
! Otherwise we substitute the merged partonID:
    elseif (ipart.eq.ileg) then
      jpart = jpart + 1
      flv_out(jpart) = merged
! if ipart coincides with jleg we do not do anything, it is neglected
    end if
  end do
!
end subroutine UndoBranching
!
! The legs are strongly ordered: ileg < jleg < kleg < lleg
subroutine UndoBranchingMod(flv_in,flv_out,ileg,jleg,kleg,lleg)
implicit none
!
  integer , dimension(:) , intent(in) :: flv_in
  integer , dimension(:) , intent(out) :: flv_out
  integer , intent(in) :: ileg,jleg
  integer , optional , intent(in) :: kleg,lleg
!
  integer :: ipart,jpart,npart
  integer :: mergedij,mergedijk,mergedkl
!
!
  npart = size(flv_in)
!
! One splitting:
  if (.not.present(kleg).and..not.present(lleg)) then
! ISR:
    if (ileg.lt.3) then
      mergedij = flv_in(ileg) - flv_in(jleg)
! FSR:
    else
      mergedij = flv_in(ileg) + flv_in(jleg)
    end if
! Constructing the underlying configuration:
    jpart = 0
    do ipart=1,npart
! If the position does not correspond neither to the
! emitter nor the radiated one simply copy the particle:
      if ((ipart.ne.ileg).and.(ipart.ne.jleg)) then
        jpart = jpart + 1
        flv_out(jpart) = flv_in(ipart)
      elseif (ipart.eq.ileg) then
        jpart = jpart + 1
        flv_out(jpart) = mergedij
      end if
    end do
    return
! One splitting, triple collinear:
  elseif (present(kleg).and..not.present(lleg)) then
! ISR:
    if (ileg.lt.3) then
      mergedijk = flv_in(ileg) - flv_in(jleg) - flv_in(kleg)
! FSR:
    else
      mergedijk = flv_in(ileg) + flv_in(jleg) + flv_in(kleg)
    end if
    jpart=0
    do ipart=1,npart
      if ((ipart.ne.ileg).and. &
          (ipart.ne.jleg).and. &
          (ipart.ne.kleg)) then
        jpart = jpart + 1
        flv_out(jpart) = flv_in(ipart)
      elseif (ipart.eq.ileg) then
        jpart = jpart + 1
        flv_out(jpart) = mergedijk
      end if
    end do
    return
! Two splittings:
  elseif (present(kleg).and.present(lleg)) then
! ISR - ISR:
    if ((ileg.lt.3).and.(kleg.lt.3)) then
      mergedij = flv_in(ileg) - flv_in(jleg)
      mergedkl = flv_in(kleg) - flv_in(lleg)
! ISR - FSR:
    elseif ((ileg.lt.3).and.(kleg.ge.3)) then
      mergedij = flv_in(ileg) - flv_in(jleg)
      mergedkl = flv_in(kleg) + flv_in(lleg)
! FSR - ISR:
    elseif ((ileg.ge.3).and.(kleg.lt.3)) then
      mergedij = flv_in(ileg) + flv_in(jleg)
      mergedkl = flv_in(kleg) - flv_in(lleg)
! FSR - FSR:
    elseif ((ileg.ge.3).and.(kleg.ge.3)) then
      mergedij = flv_in(ileg) + flv_in(jleg)
      mergedkl = flv_in(kleg) + flv_in(lleg)
    end if
    jpart = 0
    do ipart=1,npart
      if ((ipart.ne.ileg).and. &
          (ipart.ne.jleg).and. &
          (ipart.ne.kleg).and. &
          (ipart.ne.lleg)) then
        jpart = jpart + 1
        flv_out(jpart) = flv_in(ipart)
      elseif (ipart.eq.ileg) then
        jpart = jpart + 1
        flv_out(jpart) = mergedij
      elseif (ipart.eq.kleg) then
        jpart = jpart + 1
        flv_out(jpart) = mergedkl
      end if
    end do
    return
  end if
!
end subroutine UndoBranchingMod
!
subroutine FindSubTermsNLO(nBorn,flv_arr_Born,flv_Real,Cir,Sr,CSir)
use QCDparams
use utils
use particles
implicit none
!
  integer , intent(in) :: nBorn
  integer , dimension(:,:) , intent(in) :: flv_arr_Born
  integer , dimension(:) , intent(in) :: flv_Real
  type(subterms) , intent(out) :: Cir,Sr,CSir
!
  integer :: ileg,jleg
  integer :: ipart,jpart
  integer :: iterm,irad
  integer :: istat
  integer :: nleg
  integer :: Smp1,Sm
  integer , dimension(:) , allocatable :: uborn_tmp
  type(subterms) :: Cir_tmp,Sr_tmp,CSir_tmp
!
!
! We determine the number of legs for the Real-type contribution:
  nleg = size(flv_Real)
!  call PrintSubProc(flv_Real)
  allocate(uborn_tmp(nleg-1),stat=istat)
  if (istat.ne.0) then
    print *,"Allocation of uborn_tmp in FindSubTermsNLO went wrong..."
    stop
  end if
! ----------------------------------------------------------------------
! We start with C type subtractions:
! ----------------------------------------------------------------------
! We allocate a temporary array to hold the Cir terms:
  allocate(Cir_tmp%term(1:nleg*(nleg-1)),stat=istat)
  if (istat.ne.0) then
    print *,"Problem occured during allocating Cir_tmp..."
    stop
  end if
!
  Cir_tmp%numterm = 0
!
! TODO thorough checks are needed for ISR:
!***********
!*** ISR ***
!***********
  do ileg=1,2
! Collinear emission can only come from massless partons:
    if (abs(flv_Real(ileg)).gt.qcd_nf) cycle
    do jleg=3,nleg
      if (abs(flv_Real(jleg)).gt.qcd_nf) cycle
! g -> g  g
! g -> q  q~ or q~ q
      if (((flv_Real(ileg).eq.0).and.(flv_Real(jleg).eq.0)).or. &
          (flv_Real(ileg).eq.flv_Real(jleg))) then
! We construct the possible underlying Born:
        call UndoBranching(ileg,jleg,flv_Real,uborn_tmp)
! We have to check whether the underlying Born really exists or not:
        if (.not.CheckEquivProc(nleg-1,nBorn,flv_arr_Born, &
                                uborn_tmp)) then
          cycle
        end if
        Cir_tmp%numterm = Cir_tmp%numterm + 1
        allocate(Cir_tmp%term(Cir_tmp%numterm)%emit(1:1))
        allocate(Cir_tmp%term(Cir_tmp%numterm)%rad(1:2))
        allocate(Cir_tmp%term(Cir_tmp%numterm)%uborn(1:nleg-1))
! In these cases the radiator is a gluon:
        Cir_tmp%term(Cir_tmp%numterm)%emit(1)%ID = 0
! If i is zero it indicates that the radiation is ISR:
        Cir_tmp%term(Cir_tmp%numterm)%emit(1)%i  = 0
        Cir_tmp%term(Cir_tmp%numterm)%rad(1)%ID = flv_Real(ileg)
        Cir_tmp%term(Cir_tmp%numterm)%rad(2)%ID = flv_Real(jleg)
        Cir_tmp%term(Cir_tmp%numterm)%rad(1)%i  = ileg
        Cir_tmp%term(Cir_tmp%numterm)%rad(2)%i  = jleg
! We copy even the underlying Born:
        Cir_tmp%term(Cir_tmp%numterm)%uborn = uborn_tmp
      end if
! q  -> q  g
! q~ -> q~ g
      if (((flv_Real(ileg).ne.0).and.(flv_Real(jleg).eq.0)).or. &
          ((flv_Real(ileg).eq.0).and.(flv_Real(jleg).ne.0))) then
! We construct the possible underlying Born:
        call UndoBranching(ileg,jleg,flv_Real,uborn_tmp)
! We have to check whether the underlying Born really exists or not:
        if (.not.CheckEquivProc(nleg-1,nBorn,flv_arr_Born, &
                                uborn_tmp)) then
          cycle
        end if
        Cir_tmp%numterm = Cir_tmp%numterm + 1
        allocate(Cir_tmp%term(Cir_tmp%numterm)%emit(1:1))
        allocate(Cir_tmp%term(Cir_tmp%numterm)%rad(1:2))
        allocate(Cir_tmp%term(Cir_tmp%numterm)%uborn(1:nleg-1))
! In these cases the radiator is a quark or antiquark:
! We have to put the (anti)quark as the emitter:
        if (flv_Real(ileg).ne.0) then
          Cir_tmp%term(Cir_tmp%numterm)%emit(1)%ID = flv_Real(ileg)
        else
          Cir_tmp%term(Cir_tmp%numterm)%emit(1)%ID = flv_Real(jleg)
        end if
! If i is zero it indicates that the radiation is ISR:
        Cir_tmp%term(Cir_tmp%numterm)%emit(1)%i  = 0
        Cir_tmp%term(Cir_tmp%numterm)%rad(1)%ID = flv_Real(ileg)
        Cir_tmp%term(Cir_tmp%numterm)%rad(2)%ID = flv_Real(jleg)
        Cir_tmp%term(Cir_tmp%numterm)%rad(1)%i  = ileg
        Cir_tmp%term(Cir_tmp%numterm)%rad(2)%i  = jleg
! We copy even the underlying Born:
        Cir_tmp%term(Cir_tmp%numterm)%uborn = uborn_tmp
      end if
    end do
  end do
!***********
!*** FSR ***
!***********
  do ileg=3,nleg-1
! Collinear emission can only come from massless partons:
    if (abs(flv_Real(ileg)).gt.qcd_nf) cycle
    do jleg=ileg+1,nleg
      if (abs(flv_Real(jleg)).gt.qcd_nf) cycle
! g -> g  g
! g -> q  q~ or q~ q
      if (((flv_Real(ileg).eq.0).and.(flv_Real(jleg).eq.0)).or. &
          (flv_Real(ileg).eq.-flv_Real(jleg))) then
! We construct the possible underlying Born:
        call UndoBranching(ileg,jleg,flv_Real,uborn_tmp)
!
!        print *,flv_Real(ileg),flv_Real(jleg)
!        print *,ileg,jleg
!        print *,"underlying Born: "
!        call PrintSubProc(uborn_tmp)
!
! We have to check whether the underlying Born really exists or not:
        if (.not.CheckEquivProc(nleg-1,nBorn,flv_arr_Born, &
                                uborn_tmp)) then
          cycle
        end if
!
        Cir_tmp%numterm = Cir_tmp%numterm + 1
        allocate(Cir_tmp%term(Cir_tmp%numterm)%emit(1:1))
        allocate(Cir_tmp%term(Cir_tmp%numterm)%rad(1:2))
        allocate(Cir_tmp%term(Cir_tmp%numterm)%uborn(1:nleg-1))
! In these cases the radiator is a gluon:
        Cir_tmp%term(Cir_tmp%numterm)%emit(1)%ID = 0
! ileg and jleg are ordered: ileg < jleg, the emitter is always
! ileg:
        Cir_tmp%term(Cir_tmp%numterm)%emit(1)%i  = ileg 
! The ordering among the emitted partons is to make the quark the
! first radiated parton, that is: q + q~
        if (flv_Real(ileg).ge.flv_Real(jleg)) then
          Cir_tmp%term(Cir_tmp%numterm)%rad(1)%ID  = flv_Real(ileg)
          Cir_tmp%term(Cir_tmp%numterm)%rad(2)%ID  = flv_Real(jleg)
          Cir_tmp%term(Cir_tmp%numterm)%rad(1)%i   = ileg
          Cir_tmp%term(Cir_tmp%numterm)%rad(2)%i   = jleg
        else
          Cir_tmp%term(Cir_tmp%numterm)%rad(1)%ID  = flv_Real(jleg)
          Cir_tmp%term(Cir_tmp%numterm)%rad(2)%ID  = flv_Real(ileg)
          Cir_tmp%term(Cir_tmp%numterm)%rad(1)%i   = jleg
          Cir_tmp%term(Cir_tmp%numterm)%rad(2)%i   = ileg
        end if
! We copy even the underlying Born:
        Cir_tmp%term(Cir_tmp%numterm)%uborn = uborn_tmp
! We calculate the symmetry factors:
        Smp1 = CalcSymFact(flv_Real)
        Sm   = CalcSymFact(uborn_tmp)
!
!        print *,"sym facts: ",Smp1,Sm
!
! We store the ratio:
        Cir_tmp%term(Cir_tmp%numterm)%symfact = real(Sm)/real(Smp1)
! q -> q  g
! q~ -> q~ g
      elseif (((flv_Real(ileg).eq.0).and.(flv_Real(jleg).ne.0)).or. &
              ((flv_Real(ileg).ne.0).and.(flv_Real(jleg).eq.0))) then
! We construct the possible underlying Born:
        call UndoBranching(ileg,jleg,flv_Real,uborn_tmp)
!
!        print *,flv_Real(ileg),flv_Real(jleg)
!        print *,ileg,jleg
!        print *,"underlying Born: "
!        call PrintSubProc(uborn_tmp)
!
! We have to check whether the underlying Born really exists or not:
        if (.not.CheckEquivProc(nleg-1,nBorn,flv_arr_Born, &
                                uborn_tmp)) then
          cycle
        end if
!
        Cir_tmp%numterm = Cir_tmp%numterm + 1
        allocate(Cir_tmp%term(Cir_tmp%numterm)%emit(1:1))
        allocate(Cir_tmp%term(Cir_tmp%numterm)%rad(1:2))
        allocate(Cir_tmp%term(Cir_tmp%numterm)%uborn(1:nleg-1))
! In this case the emitter is the (anti)quark:
        if (flv_Real(ileg).ne.0) then
          Cir_tmp%term(Cir_tmp%numterm)%emit(1)%ID = flv_Real(ileg)
          Cir_tmp%term(Cir_tmp%numterm)%rad(1)%ID  = flv_Real(ileg)
          Cir_tmp%term(Cir_tmp%numterm)%rad(2)%ID  = flv_Real(jleg)
          Cir_tmp%term(Cir_tmp%numterm)%emit(1)%i  = ileg 
          Cir_tmp%term(Cir_tmp%numterm)%rad(1)%i   = ileg
          Cir_tmp%term(Cir_tmp%numterm)%rad(2)%i   = jleg
        else
          Cir_tmp%term(Cir_tmp%numterm)%emit(1)%ID = flv_Real(jleg)
          Cir_tmp%term(Cir_tmp%numterm)%rad(1)%ID  = flv_Real(jleg)
          Cir_tmp%term(Cir_tmp%numterm)%rad(2)%ID  = flv_Real(ileg)
          Cir_tmp%term(Cir_tmp%numterm)%emit(1)%i  = jleg 
          Cir_tmp%term(Cir_tmp%numterm)%rad(1)%i   = jleg
          Cir_tmp%term(Cir_tmp%numterm)%rad(2)%i   = ileg
        end if
! We copy even the underlying Born:
        Cir_tmp%term(Cir_tmp%numterm)%uborn = uborn_tmp
! We calculate the symmetry factors:
        Smp1 = CalcSymFact(flv_Real)
        Sm   = CalcSymFact(uborn_tmp)
! We store the ratio:
        Cir_tmp%term(Cir_tmp%numterm)%symfact = real(Sm)/real(Smp1)
      end if
    end do
  end do
!
! Right now we know the dimensionality of the array holding the
! subtraction terms hence we can allocate it and store the terms:
  allocate(Cir%term(1:Cir_tmp%numterm),stat=istat)
  if (istat.ne.0) then
    print *,"Problem occured during allocating Cir..."
    stop
  end if
!
! We copy everything there:
  Cir%numterm = Cir_tmp%numterm
  do iterm=1,Cir%numterm
    allocate(Cir%term(iterm)%uborn(1:nleg-1))
    allocate(Cir%term(iterm)%emit(1:1))
    allocate(Cir%term(iterm)%rad(1:2))
    Cir%term(iterm) = Cir_tmp%term(iterm)
  end do
!
! ----------------------------------------------------------------------
! S type subtractions:
! ----------------------------------------------------------------------
! We allocate a temporary array to hold the Sr terms:
  allocate(Sr_tmp%term(1:nleg),stat=istat)
  if (istat.ne.0) then
    print *,"Problem occured during allocating Sr_tmp..."
    stop
  end if
!
  Sr_tmp%numterm = 0
!
  do ileg=1,nleg
! Only gluons can go soft:
    if (abs(flv_Real(ileg)).ne.0) cycle
!    print *,"S-type: "
!    write(*,'(a,a,I0,a)') ConvertFromPDG(flv_Real(ileg)),"(",ileg,")"
! We construct the possible underlying Born:
    call UndoBranching(ileg,ileg,flv_Real,uborn_tmp)
!
!    print *,"underlying Born: "
!    call PrintSubProc(uborn_tmp)
!
! We have to check whether the underlying Born really exists or not:
    if (.not.CheckEquivProc(nleg-1,nBorn,flv_arr_Born, &
                            uborn_tmp)) then
      cycle
    end if
    Sr_tmp%numterm = Sr_tmp%numterm + 1
    allocate(Sr_tmp%term(Sr_tmp%numterm)%emit(1:1))
    allocate(Sr_tmp%term(Sr_tmp%numterm)%rad(1:1))
    allocate(Sr_tmp%term(Sr_tmp%numterm)%uborn(1:nleg-1))
! This time the -1 for the emitter leg indicates that the subtraction term
! is for a soft-type of radiation:
    Sr_tmp%term(Sr_tmp%numterm)%emit(1)%ID = 0
    Sr_tmp%term(Sr_tmp%numterm)%emit(1)%i  = -1
    Sr_tmp%term(Sr_tmp%numterm)%rad(1)%ID = flv_Real(ileg)
    Sr_tmp%term(Sr_tmp%numterm)%rad(1)%i  = ileg
! We copy even the underlying Born:
    Sr_tmp%term(Sr_tmp%numterm)%uborn = uborn_tmp
! We calculate the symmetry factors:
    Smp1 = CalcSymFact(flv_Real)
    Sm   = CalcSymFact(uborn_tmp)
! We store the ratio:
    Sr_tmp%term(Sr_tmp%numterm)%symfact = real(Sm)/real(Smp1)
  end do
!
! Right now we know the dimensionality of the array holding the
! subtraction terms hence we can allocate it and store the terms:
  allocate(Sr%term(1:Sr_tmp%numterm),stat=istat)
  if (istat.ne.0) then
    print *,"Problem occured during allocating Sr..."
    stop
  end if
! We copy everything there:
  Sr%numterm = Sr_tmp%numterm
  do iterm=1,Sr%numterm
    allocate(Sr%term(iterm)%uborn(1:nleg-1))
    allocate(Sr%term(iterm)%emit(1:1))
    allocate(Sr%term(iterm)%rad(1:1))
    Sr%term(iterm) = Sr_tmp%term(iterm)
  end do
!
! ----------------------------------------------------------------------
! CS type subtractions:
! ----------------------------------------------------------------------
! We allocate a temporary array to hold the CSir terms:
  allocate(CSir_tmp%term(1:nleg*(nleg-1)),stat=istat)
  if (istat.ne.0) then
    print *,"Problem occured during allocating CSir_tmp..."
    stop
  end if
!
  CSir_tmp%numterm = 0
!
! CS is a subselection of the C-type subtraction terms, hence we 
! only have to go through the C-type terms and isolate the 
! contributing ones:
!  print *,"Number of Cir terms: ",Cir%numterm
  do iterm=1,Cir%numterm
!    print *,"iterm: ",iterm
! To have a nonzero contribution one of the radiated partons has to be
! a gluon:
    if ((Cir%term(iterm)%rad(1)%ID.ne.0).and. &
        (Cir%term(iterm)%rad(2)%ID.ne.0)) cycle
    do irad=1,2
! We found one CS term:
      if (Cir%term(iterm)%rad(irad)%ID.eq.0) then
!        print *,"CS-type: "
!        print *,"Cir: ",Cir%term(iterm)%rad
!        print *,"uborn: ",Cir%term(iterm)%uborn
!        print *,"irad: ",irad
        CSir_tmp%numterm = CSir_tmp%numterm + 1
        allocate(CSir_tmp%term(CSir_tmp%numterm)%emit(1:1))
        allocate(CSir_tmp%term(CSir_tmp%numterm)%rad(1:2))
        allocate(CSir_tmp%term(CSir_tmp%numterm)%uborn(1:nleg-1))
        CSir_tmp%term(CSir_tmp%numterm)%emit = &
          Cir_tmp%term(iterm)%emit
! The ordering is such that the second index corresponds to r:
        CSir_tmp%term(CSir_tmp%numterm)%rad(1) = &
          Cir_tmp%term(iterm)%rad(mod(irad,2) + 1)
        CSir_tmp%term(CSir_tmp%numterm)%rad(2) = &
          Cir_tmp%term(iterm)%rad(irad)
        CSir_tmp%term(CSir_tmp%numterm)%uborn = &
          Cir_tmp%term(iterm)%uborn
! We calculate the symmetry factors:
        Smp1 = CalcSymFact(flv_Real)
        Sm   = CalcSymFact(uborn_tmp)
! We store the ratio:
        CSir_tmp%term(CSir_tmp%numterm)%symfact = real(Sm)/real(Smp1)
      end if
    end do
  end do
!
! Right now we know the dimensionality of the array holding the
! subtraction terms hence we can allocate it and store the terms:
  allocate(CSir%term(1:CSir_tmp%numterm),stat=istat)
  if (istat.ne.0) then
    print *,"Problem occured during allocating CSir..."
    stop
  end if
!
! We copy everything there:
  CSir%numterm = CSir_tmp%numterm
  do iterm=1,CSir%numterm
    allocate(CSir%term(iterm)%uborn(1:nleg-1))
    allocate(CSir%term(iterm)%emit(1:1))
    allocate(CSir%term(iterm)%rad(1:2))
    CSir%term(iterm) = CSir_tmp%term(iterm)
  end do
!
  deallocate(Cir_tmp%term,  &
             Sr_tmp%term,   &
             CSir_tmp%term, &
             uborn_tmp)
!
end subroutine FindSubTermsNLO
!
subroutine FindSubTermsNNLO(nBorn,flv_arr_B,flv_RR, &
                            Cirs,Cirjs,CSirs,Srs)
use QCDparams
use utils
use particles
implicit none
!
  integer , intent(in) :: nBorn
  integer , dimension(:,:) , intent(in) :: flv_arr_B
  integer , dimension(:) , intent(in) :: flv_RR
  type(subterms) , intent(out) :: Cirs,Cirjs, &
                                  CSirs,Srs
!
  integer :: ileg,rleg,jleg,sleg
  integer :: ileg_uborn,jleg_uborn,sleg_uborn
  integer :: iterm
  integer :: istat
  integer :: nleg
  integer :: Smp2,Sm
  integer , dimension(:) , allocatable :: uborn_tmp
  type(subterms) :: Cirs_tmp,Cirjs_tmp,CSirs_tmp,Srs_tmp
!
!
! Obtaining the number of legs for the double-real configuration:
  nleg = size(flv_RR)
!  call PrintSubproc(flv_RR)
!  print *,"number of legs: ",nleg
! Allocating a temporary array to hold the underlying Born
! configuration:
  allocate(uborn_tmp(nleg-2),stat=istat)
  if (istat.ne.0) then
    print *,"Allocation of uborn_tmp in FindSubTermsNNLO went wrong..."
    stop
  end if
! ----------------------------------------------------------------------
! ------------------------------- Cirs: --------------------------------
! ----------------------------------------------------------------------
! Allocating a temporary array to hold the Cirs subtraction terms,
! the size is only guessed:
  allocate(Cirs_tmp%term(1:nleg*(nleg-1)*(nleg-2)),stat=istat)
  if (istat.ne.0) then
    print *,"Allocation went wrong for Cirs_tmp..."
    stop
  end if
!
  Cirs_tmp%numterm = 0
!
! TODO fix the initial state too:
!***********
!*** ISR ***
!***********
!***********
!*** FSR ***
!***********
  do ileg=3,nleg-2
    if (abs(flv_RR(ileg)).gt.qcd_nf) cycle
    do rleg=ileg+1,nleg-1
      if (abs(flv_RR(rleg)).gt.qcd_nf) cycle
      do sleg=rleg+1,nleg
        if (abs(flv_RR(sleg)).gt.qcd_nf) cycle
! All possible branchings can be tested by examining all
! possible pairs and looking for equal flavors but opposite
! in sign:
        if (.not.((flv_RR(ileg).eq.-flv_RR(rleg)).or. &
                  (flv_RR(ileg).eq.-flv_RR(sleg)).or. &
                  (flv_RR(rleg).eq.-flv_RR(sleg)))) cycle
! Constructing the underlying Born:
        call UndoBranchingMod(flv_RR,uborn_tmp,ileg,rleg,sleg)
        if (.not.CheckEquivProc(nleg-2,nBorn,flv_arr_B, &
                                uborn_tmp)) then
          cycle
        end if
!
!        print *,"ileg,rleg,sleg: ",ileg,rleg,sleg
!        print *," ileg: ",ConvertFromPDG(flv_RR(ileg)), &
!                " rleg: ",ConvertFromPDG(flv_RR(rleg)), &
!                " sleg: ",ConvertFromPDG(flv_RR(sleg))
!        call PrintSubProc(uborn_tmp)
!
        Cirs_tmp%numterm = Cirs_tmp%numterm + 1
        allocate(Cirs_tmp%term(Cirs_tmp%numterm)%emit(1:1))
        allocate(Cirs_tmp%term(Cirs_tmp%numterm)%rad(1:3))
        allocate(Cirs_tmp%term(Cirs_tmp%numterm)%uborn(1:nleg-2))
! The particle ID of the emitter can be easily obtained,
! since it is at position ileg in uborn_tmp:
        Cirs_tmp%term(Cirs_tmp%numterm)%emit(1)%ID = uborn_tmp(ileg)
        Cirs_tmp%term(Cirs_tmp%numterm)%emit(1)%i  = ileg 
! The ordering between radiated partons depend upon the emitter itself:
! (anti)quark:
! Ordering: 
! q(~) -> q(~) g g or 
! q(~) -> q(~) r r~ or 
! g    -> g q q~ or 
! g    -> g g g
        if (uborn_tmp(ileg).eq.flv_RR(ileg)) then
          Cirs_tmp%term(Cirs_tmp%numterm)%rad(1)%ID  = flv_RR(ileg)
          Cirs_tmp%term(Cirs_tmp%numterm)%rad(1)%i   = ileg
          if (max(flv_RR(rleg),flv_RR(sleg)).eq.flv_RR(rleg)) then
            Cirs_tmp%term(Cirs_tmp%numterm)%rad(2)%ID  = flv_RR(rleg)
            Cirs_tmp%term(Cirs_tmp%numterm)%rad(2)%i   = rleg
            Cirs_tmp%term(Cirs_tmp%numterm)%rad(3)%ID  = flv_RR(sleg)
            Cirs_tmp%term(Cirs_tmp%numterm)%rad(3)%i   = sleg
          else
            Cirs_tmp%term(Cirs_tmp%numterm)%rad(3)%ID  = flv_RR(rleg)
            Cirs_tmp%term(Cirs_tmp%numterm)%rad(3)%i   = rleg
            Cirs_tmp%term(Cirs_tmp%numterm)%rad(2)%ID  = flv_RR(sleg)
            Cirs_tmp%term(Cirs_tmp%numterm)%rad(2)%i   = sleg
          end if
        elseif (uborn_tmp(ileg).eq.flv_RR(rleg)) then
          Cirs_tmp%term(Cirs_tmp%numterm)%rad(1)%ID  = flv_RR(rleg)
          Cirs_tmp%term(Cirs_tmp%numterm)%rad(1)%i   = rleg
          if (max(flv_RR(ileg),flv_RR(sleg)).eq.flv_RR(ileg)) then
            Cirs_tmp%term(Cirs_tmp%numterm)%rad(2)%ID  = flv_RR(ileg)
            Cirs_tmp%term(Cirs_tmp%numterm)%rad(2)%i   = ileg
            Cirs_tmp%term(Cirs_tmp%numterm)%rad(3)%ID  = flv_RR(sleg)
            Cirs_tmp%term(Cirs_tmp%numterm)%rad(3)%i   = sleg
          else
            Cirs_tmp%term(Cirs_tmp%numterm)%rad(3)%ID  = flv_RR(ileg)
            Cirs_tmp%term(Cirs_tmp%numterm)%rad(3)%i   = ileg
            Cirs_tmp%term(Cirs_tmp%numterm)%rad(2)%ID  = flv_RR(sleg)
            Cirs_tmp%term(Cirs_tmp%numterm)%rad(2)%i   = sleg
          end if
        elseif (uborn_tmp(ileg).eq.flv_RR(sleg)) then
          Cirs_tmp%term(Cirs_tmp%numterm)%rad(1)%ID  = flv_RR(sleg)
          Cirs_tmp%term(Cirs_tmp%numterm)%rad(1)%i   = sleg
          if (max(flv_RR(rleg),flv_RR(ileg)).eq.flv_RR(rleg)) then
            Cirs_tmp%term(Cirs_tmp%numterm)%rad(2)%ID  = flv_RR(rleg)
            Cirs_tmp%term(Cirs_tmp%numterm)%rad(2)%i   = rleg
            Cirs_tmp%term(Cirs_tmp%numterm)%rad(3)%ID  = flv_RR(ileg)
            Cirs_tmp%term(Cirs_tmp%numterm)%rad(3)%i   = ileg
          else
            Cirs_tmp%term(Cirs_tmp%numterm)%rad(3)%ID  = flv_RR(rleg)
            Cirs_tmp%term(Cirs_tmp%numterm)%rad(3)%i   = rleg
            Cirs_tmp%term(Cirs_tmp%numterm)%rad(2)%ID  = flv_RR(ileg)
            Cirs_tmp%term(Cirs_tmp%numterm)%rad(2)%i   = ileg
          end if
        end if
! We copy even the underlying Born:
        Cirs_tmp%term(Cirs_tmp%numterm)%uborn = uborn_tmp
! We calculate the symmetry factors:
        Smp2 = CalcSymFact(flv_RR)
        Sm   = CalcSymFact(uborn_tmp)
!
!        print *,"sym facts: ",Smp2,Sm
!        print *,"emit: ",Cirs_tmp%term(Cirs_tmp%numterm)%emit(1)%i, &
!                         ConvertFromPDG(Cirs_tmp%term(Cirs_tmp%numterm)%emit(1)%ID)
!        print *,"rad1: ",Cirs_tmp%term(Cirs_tmp%numterm)%rad(1)%i, &
!                         ConvertFromPDG(Cirs_tmp%term(Cirs_tmp%numterm)%rad(1)%ID)
!        print *,"rad2: ",Cirs_tmp%term(Cirs_tmp%numterm)%rad(2)%i, &
!                         ConvertFromPDG(Cirs_tmp%term(Cirs_tmp%numterm)%rad(2)%ID)
!        print *,"rad3: ",Cirs_tmp%term(Cirs_tmp%numterm)%rad(3)%i, &
!                         ConvertFromPDG(Cirs_tmp%term(Cirs_tmp%numterm)%rad(3)%ID)
!
! We store the ratio:
        Cirs_tmp%term(Cirs_tmp%numterm)%symfact = real(Sm)/real(Smp2)
      end do
    end do
  end do
!
! Right now we know the dimensionality of the array holding the
! subtraction terms hence we can allocate it and store the terms:
  allocate(Cirs%term(1:Cirs_tmp%numterm),stat=istat)
  if (istat.ne.0) then
    print *,"Problem occured during allocating Cirs..."
    stop
  end if
!
! We copy everything there:
  Cirs%numterm = Cirs_tmp%numterm
  do iterm=1,Cirs%numterm
    allocate(Cirs%term(iterm)%uborn(1:nleg-2))
    allocate(Cirs%term(iterm)%emit(1:1))
    allocate(Cirs%term(iterm)%rad(1:3))
    Cirs%term(iterm) = Cirs_tmp%term(iterm)
  end do
!
! ----------------------------------------------------------------------
! ------------------------------ Cirjs: --------------------------------
! ----------------------------------------------------------------------
! Allocating a temporary array to hold the Cirjs subtraction terms,
! the size is only guessed:
  allocate(Cirjs_tmp%term(1:nleg*(nleg-1)*(nleg-2)*(nleg-3)),stat=istat)
  if (istat.ne.0) then
    print *,"Allocation went wrong for Cirjs_tmp..."
    stop
  end if
!
  Cirjs_tmp%numterm = 0
!
! TODO fix the initial state too:
!***********
!*** ISR ***
!***********
!***********
!*** FSR ***
!***********
  do ileg=3,nleg-2
    if (abs(flv_RR(ileg)).gt.qcd_nf) cycle
    do rleg=ileg+1,nleg
      if (abs(flv_RR(rleg)).gt.qcd_nf) cycle
! Test upon the first pairing:
      if (.not.(((flv_RR(ileg)+flv_RR(rleg)).eq.0).or. &
                ((flv_RR(ileg)+flv_RR(rleg)).eq.flv_RR(ileg)).or. &
                ((flv_RR(ileg)+flv_RR(rleg)).eq.flv_RR(rleg)))) cycle
      do jleg=ileg+1,nleg-1
        if (abs(flv_RR(jleg)).gt.qcd_nf) cycle
! The second emitter cannot coincide with the radiated one from
! the first branching:
        if (jleg.eq.rleg) cycle
        do sleg=jleg+1,nleg
          if (abs(flv_RR(sleg)).gt.qcd_nf) cycle
! The second radiated parton cannot coincide with the one coming
! from the first pairing:
          if (rleg.eq.sleg) cycle
! Test upon the second pairing:
          if (.not.(((flv_RR(jleg)+flv_RR(sleg)).eq.0).or. &
                    ((flv_RR(jleg)+flv_RR(sleg)).eq.flv_RR(jleg)).or. &
                    ((flv_RR(jleg)+flv_RR(sleg)).eq.flv_RR(sleg)))) cycle
!          print *,"ileg,rleg,jleg,sleg: ",ileg,rleg,jleg,sleg
! Constructing the underlying Born:
          call UndoBranchingMod(flv_RR,uborn_tmp,ileg,rleg,jleg,sleg)
          if (.not.CheckEquivProc(nleg-2,nBorn,flv_arr_B, &
                                  uborn_tmp)) then
            cycle
          end if
!
!          print *,"ileg,rleg,jleg,sleg: ",ileg,rleg,jleg,sleg
!          print *," ileg: ",ConvertFromPDG(flv_RR(ileg)), &
!                  " rleg: ",ConvertFromPDG(flv_RR(rleg)), &
!                  " jleg: ",ConvertFromPDG(flv_RR(jleg)), &
!                  " sleg: ",ConvertFromPDG(flv_RR(sleg))
!          call PrintSubProc(uborn_tmp)
!
! Storing the information on the singular region:
! Increase the counter by one unit:
          Cirjs_tmp%numterm = Cirjs_tmp%numterm + 1
          allocate(Cirjs_tmp%term(Cirjs_tmp%numterm)%emit(1:2))
          allocate(Cirjs_tmp%term(Cirjs_tmp%numterm)%rad(1:4))
          allocate(Cirjs_tmp%term(Cirjs_tmp%numterm)%uborn(1:nleg-2))
! The legs were ordered such ileg < rleg , jleg < sleg and ileg < jleg,
! hence the corresponding emitter legs are ileg and jleg:
          Cirjs_tmp%term(Cirjs_tmp%numterm)%emit(1)%ID = uborn_tmp(ileg)
          Cirjs_tmp%term(Cirjs_tmp%numterm)%emit(1)%i  = ileg 
! It is possible that rleg is less than jleg in that case jleg has to be
! decreased by one:
          if (rleg.lt.jleg) then
            Cirjs_tmp%term(Cirjs_tmp%numterm)%emit(2)%ID = uborn_tmp(jleg-1)
            Cirjs_tmp%term(Cirjs_tmp%numterm)%emit(2)%i  = jleg - 1
            jleg_uborn = jleg - 1
          else
            Cirjs_tmp%term(Cirjs_tmp%numterm)%emit(2)%ID = uborn_tmp(jleg)
            Cirjs_tmp%term(Cirjs_tmp%numterm)%emit(2)%i  = jleg 
            jleg_uborn = jleg
          end if
! The ordering among the radiated partons is:
! * q(~) -> q(~) g
! * g    -> g    g
! * g    -> q    q~
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!! Considering the first radiated pair: !!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The emitter was a(n) (anti)quark:
          if (uborn_tmp(ileg).ne.0) then
! The (anti)quark situates at position ileg:
! The other parton is a gluon hence its ID is zero...
            if (flv_RR(ileg).ne.0) then
              Cirjs_tmp%term(Cirjs_tmp%numterm)%rad(1)%ID = flv_RR(ileg)
              Cirjs_tmp%term(Cirjs_tmp%numterm)%rad(1)%i  = ileg 
              Cirjs_tmp%term(Cirjs_tmp%numterm)%rad(2)%ID = flv_RR(rleg)
              Cirjs_tmp%term(Cirjs_tmp%numterm)%rad(2)%i  = rleg 
            else
              Cirjs_tmp%term(Cirjs_tmp%numterm)%rad(1)%ID = flv_RR(rleg)
              Cirjs_tmp%term(Cirjs_tmp%numterm)%rad(1)%i  = rleg 
              Cirjs_tmp%term(Cirjs_tmp%numterm)%rad(2)%ID = flv_RR(ileg)
              Cirjs_tmp%term(Cirjs_tmp%numterm)%rad(2)%i  = ileg 
            end if
! Having a gluon as emitter:
          else
! If a pair of gluons is radiated simply store the original order:
            if (flv_RR(ileg).eq.0) then
              Cirjs_tmp%term(Cirjs_tmp%numterm)%rad(1)%ID = flv_RR(ileg)
              Cirjs_tmp%term(Cirjs_tmp%numterm)%rad(1)%i  = ileg 
              Cirjs_tmp%term(Cirjs_tmp%numterm)%rad(2)%ID = flv_RR(rleg)
              Cirjs_tmp%term(Cirjs_tmp%numterm)%rad(2)%i  = rleg 
! If a quark-antiquark pair is emitted the quark comes first:
            else
              if (flv_RR(ileg).gt.0) then
                Cirjs_tmp%term(Cirjs_tmp%numterm)%rad(1)%ID = flv_RR(ileg)
                Cirjs_tmp%term(Cirjs_tmp%numterm)%rad(1)%i  = ileg 
                Cirjs_tmp%term(Cirjs_tmp%numterm)%rad(2)%ID = flv_RR(rleg)
                Cirjs_tmp%term(Cirjs_tmp%numterm)%rad(2)%i  = rleg 
              else
                Cirjs_tmp%term(Cirjs_tmp%numterm)%rad(1)%ID = flv_RR(rleg)
                Cirjs_tmp%term(Cirjs_tmp%numterm)%rad(1)%i  = rleg 
                Cirjs_tmp%term(Cirjs_tmp%numterm)%rad(2)%ID = flv_RR(ileg)
                Cirjs_tmp%term(Cirjs_tmp%numterm)%rad(2)%i  = ileg 
              end if
            end if
          end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!! Considering the second radiated pair: !!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The emitter was a(n) (anti)quark:
          if (uborn_tmp(jleg_uborn).ne.0) then
! The (anti)quark situates at position jleg:
! The other parton is a gluon hence its ID is zero...
            if (flv_RR(jleg).ne.0) then
              Cirjs_tmp%term(Cirjs_tmp%numterm)%rad(3)%ID = flv_RR(jleg)
              Cirjs_tmp%term(Cirjs_tmp%numterm)%rad(3)%i  = jleg 
              Cirjs_tmp%term(Cirjs_tmp%numterm)%rad(4)%ID = flv_RR(sleg)
              Cirjs_tmp%term(Cirjs_tmp%numterm)%rad(4)%i  = sleg 
            else
              Cirjs_tmp%term(Cirjs_tmp%numterm)%rad(3)%ID = flv_RR(sleg)
              Cirjs_tmp%term(Cirjs_tmp%numterm)%rad(3)%i  = sleg 
              Cirjs_tmp%term(Cirjs_tmp%numterm)%rad(4)%ID = flv_RR(jleg)
              Cirjs_tmp%term(Cirjs_tmp%numterm)%rad(4)%i  = jleg 
            end if
! Having a gluon as emitter:
          else
! If a pair of gluons is radiated simply store the original order:
            if (flv_RR(jleg).eq.0) then
              Cirjs_tmp%term(Cirjs_tmp%numterm)%rad(3)%ID = flv_RR(jleg)
              Cirjs_tmp%term(Cirjs_tmp%numterm)%rad(3)%i  = jleg 
              Cirjs_tmp%term(Cirjs_tmp%numterm)%rad(4)%ID = flv_RR(sleg)
              Cirjs_tmp%term(Cirjs_tmp%numterm)%rad(4)%i  = sleg 
! If a quark-antiquark pair is emitted the quark comes first:
            else
              if (flv_RR(jleg).gt.0) then
                Cirjs_tmp%term(Cirjs_tmp%numterm)%rad(3)%ID = flv_RR(jleg)
                Cirjs_tmp%term(Cirjs_tmp%numterm)%rad(3)%i  = jleg 
                Cirjs_tmp%term(Cirjs_tmp%numterm)%rad(4)%ID = flv_RR(sleg)
                Cirjs_tmp%term(Cirjs_tmp%numterm)%rad(4)%i  = sleg 
              else
                Cirjs_tmp%term(Cirjs_tmp%numterm)%rad(3)%ID = flv_RR(sleg)
                Cirjs_tmp%term(Cirjs_tmp%numterm)%rad(3)%i  = sleg 
                Cirjs_tmp%term(Cirjs_tmp%numterm)%rad(4)%ID = flv_RR(jleg)
                Cirjs_tmp%term(Cirjs_tmp%numterm)%rad(4)%i  = jleg 
              end if
            end if
          end if
! We copy even the underlying Born:
          Cirjs_tmp%term(Cirjs_tmp%numterm)%uborn = uborn_tmp
! We calculate the symmetry factors:
          Smp2 = CalcSymFact(flv_RR)
          Sm   = CalcSymFact(uborn_tmp)
! We store the ratio:
          Cirjs_tmp%term(Cirjs_tmp%numterm)%symfact = real(Sm)/real(Smp2)
        end do
      end do
    end do
  end do
!
! Right now we know the dimensionality of the array holding the
! subtraction terms hence we can allocate it and store the terms:
  allocate(Cirjs%term(1:Cirjs_tmp%numterm),stat=istat)
  if (istat.ne.0) then
    print *,"Problem occured during allocating Cirjs..."
    stop
  end if
!
! We copy everything there:
  Cirjs%numterm = Cirjs_tmp%numterm
  do iterm=1,Cirjs%numterm
    allocate(Cirjs%term(iterm)%uborn(1:nleg-2))
    allocate(Cirjs%term(iterm)%emit(1:2))
    allocate(Cirjs%term(iterm)%rad(1:4))
    Cirjs%term(iterm) = Cirjs_tmp%term(iterm)
  end do
!
! ----------------------------------------------------------------------
! ------------------------------ CSirs: --------------------------------
! ----------------------------------------------------------------------
! Allocating a temporary array to hold the CSirs subtraction terms,
! the size is only guessed:
  allocate(CSirs_tmp%term(1:nleg*(nleg-1)*(nleg-2)),stat=istat)
  if (istat.ne.0) then
    print *,"Allocation went wrong for CSirs_tmp..."
    stop
  end if
!
  CSirs_tmp%numterm = 0
!
! TODO fix the initial state too:
!***********
!*** ISR ***
!***********
!***********
!*** FSR ***
!***********
  do ileg=3,nleg-1
    if (abs(flv_RR(ileg)).gt.qcd_nf) cycle
    do rleg=ileg+1,nleg
      if (abs(flv_RR(rleg)).gt.qcd_nf) cycle
! Test upon pairing:
      if (.not.(((flv_RR(ileg)+flv_RR(rleg)).eq.0).or. &
                ((flv_RR(ileg)+flv_RR(rleg)).eq.flv_RR(ileg)).or. &
                ((flv_RR(ileg)+flv_RR(rleg)).eq.flv_RR(rleg)))) cycle
! An additional loop is needed to find a the soft leg:
      do sleg=3,nleg
! The soft leg is different from the collinear pair:
        if ((ileg.eq.sleg).or.(rleg.eq.sleg)) cycle
! It can only be a gluon:
        if (flv_RR(sleg).ne.0) cycle
! Since the soft parton can only be a gluon the same routine
! can be used to determine the underlying Born as for the triple-Collinear:
        call UndoBranchingMod(flv_RR,uborn_tmp,ileg,rleg,sleg)
        if (.not.CheckEquivProc(nleg-2,nBorn,flv_arr_B, &
                                uborn_tmp)) then
          cycle
        end if
!
!        print *,"ileg,rleg,sleg: ",ileg,rleg,sleg
!        print *," ileg: ",ConvertFromPDG(flv_RR(ileg)), &
!                " rleg: ",ConvertFromPDG(flv_RR(rleg)), &
!                " sleg: ",ConvertFromPDG(flv_RR(sleg))
!        call PrintSubProc(uborn_tmp)
! 
! At this point the subtraction term can be stored:
        CSirs_tmp%numterm = CSirs_tmp%numterm + 1
        allocate(CSirs_tmp%term(CSirs_tmp%numterm)%emit(1:2))
        allocate(CSirs_tmp%term(CSirs_tmp%numterm)%rad(1:3))
        allocate(CSirs_tmp%term(CSirs_tmp%numterm)%uborn(1:nleg-2))
! The ileg and rleg legs are ordered: ileg < rleg, but sleg
! can be everywhere, hence the positions can change:
        ileg_uborn = ileg
! If the soft leg comes before the emitter for the collinear pair 
! it has to be rescaled:
        if (sleg.lt.ileg) ileg_uborn = ileg_uborn - 1
!
        CSirs_tmp%term(CSirs_tmp%numterm)%emit(1)%ID = uborn_tmp(ileg_uborn)
        CSirs_tmp%term(CSirs_tmp%numterm)%emit(1)%i  = ileg_uborn 
        CSirs_tmp%term(CSirs_tmp%numterm)%emit(2)%ID = 0
        CSirs_tmp%term(CSirs_tmp%numterm)%emit(2)%i  = -1
! For the collinear pair the usual convention is followed:
! * q(~) -> q(~) g
! * g    -> g    g
! * g    -> q    q~
! The emitter was a(n) (anti)quark:
        if (uborn_tmp(ileg_uborn).ne.0) then
! The (anti)quark situates at position ileg:
! The other parton is a gluon hence its ID is zero...
          if (flv_RR(ileg).ne.0) then
            CSirs_tmp%term(CSirs_tmp%numterm)%rad(1)%ID = flv_RR(ileg)
            CSirs_tmp%term(CSirs_tmp%numterm)%rad(1)%i  = ileg 
            CSirs_tmp%term(CSirs_tmp%numterm)%rad(2)%ID = flv_RR(rleg)
            CSirs_tmp%term(CSirs_tmp%numterm)%rad(2)%i  = rleg 
          else
            CSirs_tmp%term(CSirs_tmp%numterm)%rad(1)%ID = flv_RR(rleg)
            CSirs_tmp%term(CSirs_tmp%numterm)%rad(1)%i  = rleg 
            CSirs_tmp%term(CSirs_tmp%numterm)%rad(2)%ID = flv_RR(ileg)
            CSirs_tmp%term(CSirs_tmp%numterm)%rad(2)%i  = ileg 
          end if
! Having a gluon as emitter:
        else
! If a pair of gluons is radiated simply store the original order:
          if (flv_RR(ileg).eq.0) then
            CSirs_tmp%term(CSirs_tmp%numterm)%rad(1)%ID = flv_RR(ileg)
            CSirs_tmp%term(CSirs_tmp%numterm)%rad(1)%i  = ileg 
            CSirs_tmp%term(CSirs_tmp%numterm)%rad(2)%ID = flv_RR(rleg)
            CSirs_tmp%term(CSirs_tmp%numterm)%rad(2)%i  = rleg 
! If a quark-antiquark pair is emitted the quark comes first:
          else
            if (flv_RR(ileg).gt.0) then
              CSirs_tmp%term(CSirs_tmp%numterm)%rad(1)%ID = flv_RR(ileg)
              CSirs_tmp%term(CSirs_tmp%numterm)%rad(1)%i  = ileg 
              CSirs_tmp%term(CSirs_tmp%numterm)%rad(2)%ID = flv_RR(rleg)
              CSirs_tmp%term(CSirs_tmp%numterm)%rad(2)%i  = rleg 
            else
              CSirs_tmp%term(CSirs_tmp%numterm)%rad(1)%ID = flv_RR(rleg)
              CSirs_tmp%term(CSirs_tmp%numterm)%rad(1)%i  = rleg 
              CSirs_tmp%term(CSirs_tmp%numterm)%rad(2)%ID = flv_RR(ileg)
              CSirs_tmp%term(CSirs_tmp%numterm)%rad(2)%i  = ileg 
            end if
          end if
        end if
! For the soft leg the situation is easy:
        CSirs_tmp%term(CSirs_tmp%numterm)%rad(3)%ID = flv_RR(sleg)
        CSirs_tmp%term(CSirs_tmp%numterm)%rad(3)%i  = sleg 
! We copy even the underlying Born:
        CSirs_tmp%term(CSirs_tmp%numterm)%uborn = uborn_tmp
! We calculate the symmetry factors:
        Smp2 = CalcSymFact(flv_RR)
        Sm   = CalcSymFact(uborn_tmp)
! We store the ratio:
        CSirs_tmp%term(CSirs_tmp%numterm)%symfact = real(Sm)/real(Smp2)
      end do
    end do
  end do
!
  allocate(CSirs%term(1:CSirs_tmp%numterm),stat=istat)
  if (istat.ne.0) then
    print *,"Problem occured during allocating CSirs..."
    stop
  end if
!
! We copy everything there:
  CSirs%numterm = CSirs_tmp%numterm
  do iterm=1,CSirs%numterm
    allocate(CSirs%term(iterm)%uborn(1:nleg-2))
    allocate(CSirs%term(iterm)%emit(1:2))
    allocate(CSirs%term(iterm)%rad(1:4))
    CSirs%term(iterm) = CSirs_tmp%term(iterm)
  end do
!
! ----------------------------------------------------------------------
! ------------------------------- Srs: ---------------------------------
! ----------------------------------------------------------------------
! Allocating a temporary array to hold the Srs subtraction terms,
! the size is only guessed:
  allocate(Srs_tmp%term(1:nleg*(nleg-1)),stat=istat)
  if (istat.ne.0) then
    print *,"Allocation went wrong for Srs_tmp..."
    stop
  end if
!
  Srs_tmp%numterm = 0
!
! Dobule loop over final state particles to find gluon- or
! quark pairs for suitable Srs counterterms:
  do rleg=3,nleg-1
    if (abs(flv_RR(rleg)).gt.qcd_nf) cycle
    do sleg=rleg+1,nleg
      if (abs(flv_RR(sleg)).gt.qcd_nf) cycle
! Test upon pairing:
! Only looking for g -> g + g or g -> q + q~ splittings:
      if (.not.(((flv_RR(rleg)+flv_RR(sleg)).eq.0))) cycle
!
! The construction of the underlying Born is easy provided 
! only partons rleg and sleg have to be omitted:
      jleg = 0
      do ileg=1,nleg
        if ((ileg.eq.rleg).or.(ileg.eq.sleg)) cycle
        jleg = jleg + 1
        uborn_tmp(jleg) = flv_RR(ileg)
      end do
      if (.not.CheckEquivProc(nleg-2,nBorn,flv_arr_B, &
                              uborn_tmp)) then
        cycle
      end if
!
!      print *,"rleg,sleg: ",rleg,sleg
!      print *," rleg: ",ConvertFromPDG(flv_RR(rleg)), &
!              " sleg: ",ConvertFromPDG(flv_RR(sleg))
!      call PrintSubProc(uborn_tmp)
!
      Srs_tmp%numterm = Srs_tmp%numterm + 1
      allocate(Srs_tmp%term(Srs_tmp%numterm)%emit(1:2))
      allocate(Srs_tmp%term(Srs_tmp%numterm)%rad(1:2))
      allocate(Srs_tmp%term(Srs_tmp%numterm)%uborn(1:nleg-2))
! Both radiated partons are soft hence no emitter is specified:
      Srs_tmp%term(Srs_tmp%numterm)%emit(1:2)%ID = 0
      Srs_tmp%term(Srs_tmp%numterm)%emit(1:2)%i  = -1
! The radiated partons are stored in the order of:
! * g -> g + g
! * g -> q + q~
! hence reordering is needed if a quark-pair is emitted but in
! the wrong order, otherwise simply dump the positions in:
      if (flv_RR(rleg).ge.flv_RR(sleg)) then
        Srs_tmp%term(Srs_tmp%numterm)%rad(1)%ID = flv_RR(rleg)
        Srs_tmp%term(Srs_tmp%numterm)%rad(2)%ID = flv_RR(sleg)
        Srs_tmp%term(Srs_tmp%numterm)%rad(1)%i  = rleg
        Srs_tmp%term(Srs_tmp%numterm)%rad(2)%i  = sleg
      else
        Srs_tmp%term(Srs_tmp%numterm)%rad(1)%ID = flv_RR(sleg)
        Srs_tmp%term(Srs_tmp%numterm)%rad(2)%ID = flv_RR(rleg)
        Srs_tmp%term(Srs_tmp%numterm)%rad(1)%i  = sleg
        Srs_tmp%term(Srs_tmp%numterm)%rad(2)%i  = rleg
      end if
! Copying the underlying Born:
      Srs_tmp%term(Srs_tmp%numterm)%uborn = uborn_tmp
! Determining symmetry factors:
      Smp2 = CalcSymFact(flv_RR)
      Sm   = CalcSymFact(uborn_tmp)
! Storing the ratio of it:
      Srs_tmp%term(Srs_tmp%numterm)%symfact = real(Sm)/real(Smp2)
    end do
  end do
!
  allocate(Srs%term(1:Srs_tmp%numterm),stat=istat)
  if (istat.ne.0) then
    print *,"Problem occured during allocating Srs..."
    stop
  end if
!
! We copy everything there:
  Srs%numterm = Srs_tmp%numterm
  do iterm=1,Srs%numterm
    allocate(Srs%term(iterm)%uborn(1:nleg-2))
    allocate(Srs%term(iterm)%emit(1:2))
    allocate(Srs%term(iterm)%rad(1:2))
    Srs%term(iterm) = Srs_tmp%term(iterm)
  end do
!
  deallocate(uborn_tmp, &
             Cirs_tmp%term, &
             Cirjs_tmp%term, &
             CSirs_tmp%term, &
             Srs_tmp%term)
!
!  stop "FindSubTermsNNLO"
!
end subroutine FindSubTermsNNLO
!
subroutine PrintCir(iun,iterm,subterm)
use utils
implicit none
!
  integer , intent(in) :: iun,iterm
  type(subtractterm) , intent(in) :: subterm
!
!
!
  write(iun,fmt='(a,I4,3(3(a),I0),a)') &
  "iterm: ",iterm," , ", &
  ConvertFromPDG(subterm%emit(1)%ID),"(",subterm%emit(1)%i,") -> ", &
  ConvertFromPDG(subterm%rad(1)%ID),"(",subterm%rad(1)%i,") || ", &
  ConvertFromPDG(subterm%rad(2)%ID),"(",subterm%rad(2)%i,")" 
  call PrintBranching(iun,subterm)
!
end subroutine PrintCir
!
subroutine PrintSr(iun,iterm,subterm)
use utils
implicit none
!
  integer , intent(in) :: iun,iterm
  type(subtractterm) , intent(in) :: subterm
!
!
!
  write(iun,fmt='(a,I4,3(a),I0,a)') "iterm: ",iterm," , ", &
    ConvertFromPDG(subterm%rad(1)%ID),"(",subterm%rad(1)%i,") -> 0" 
! We print the Underlying Born:
  write(iun,fmt='(a)',advance='no') "UBorn: "
  call PrintSubproc(iun,subterm%uborn(:))
!
end subroutine PrintSr
!
subroutine PrintCirSr(iun,iterm,subterm)
use utils
implicit none
!
  integer , intent(in) :: iun,iterm
  type(subtractterm) , intent(in) :: subterm
!
!
!
  write(iun,fmt='(a,I4,a,a,a,I0,a,a,a,I0,a,a,a,I0,a,a,a,I0,a)') &
    "iterm: ",iterm," , ",ConvertFromPDG(subterm%emit(1)%ID), &
    "(",subterm%emit(1)%i,") -> ",ConvertFromPDG(subterm%rad(1)%ID), &
    "(",subterm%rad(1)%i,") || ",ConvertFromPDG(subterm%rad(2)%ID), &
    "(",subterm%rad(2)%i,") , ",ConvertFromPDG(subterm%rad(2)%ID), &
    "(",subterm%rad(2)%i,") -> 0" 
  call PrintBranching(iun,subterm)
!
end subroutine PrintCirSr
!
subroutine PrintCirs(iun,iterm,subterm)
use utils
implicit none
!
  integer , intent(in) :: iun,iterm
  type(subtractterm) , intent(in) :: subterm
!
!
!
  write(iun,fmt='(a,I4,4(3(a),I0),a)') &
  "iterm: ",iterm," , ", &
  ConvertFromPDG(subterm%emit(1)%ID),"(",subterm%emit(1)%i,") -> ", &
  ConvertFromPDG(subterm%rad(1)%ID), "(",subterm%rad(1)%i,") || ", &
  ConvertFromPDG(subterm%rad(2)%ID), "(",subterm%rad(2)%i,") || ", & 
  ConvertFromPDG(subterm%rad(3)%ID), "(",subterm%rad(3)%i,")" 
  call PrintBranching(iun,subterm)
!
end subroutine PrintCirs
!
subroutine PrintCSirs(iun,iterm,subterm)
use utils
implicit none
!
  integer , intent(in) :: iun,iterm
  type(subtractterm) , intent(in) :: subterm
!
!
!
  write(iun,fmt='(a,I4,6(3(a),I0),a)') &
  "iterm: ",iterm," , ", &
  ConvertFromPDG(subterm%emit(1)%ID),"(",subterm%emit(1)%i,") -> ", &
  ConvertFromPDG(subterm%rad(1)%ID),"(",subterm%rad(1)%i,") || ", &
  ConvertFromPDG(subterm%rad(2)%ID),"(",subterm%rad(2)%i,")  ,  ", & 
  ConvertFromPDG(subterm%rad(3)%ID),"(",subterm%rad(3)%i,") -> 0"
  call PrintBranching(iun,subterm)
!
end subroutine PrintCSirs
!
subroutine PrintCirsCSirs(iun,iterm,subterm)
use utils
implicit none
!
  integer , intent(in) :: iun,iterm
  type(subtractterm) , intent(in) :: subterm
!
!
!
  write(iun,fmt='(a,I4,6(3(a),I0),a)') &
  "iterm: ",iterm," , ", &
  ConvertFromPDG(subterm%emit(1)%ID),"(",subterm%emit(1)%i,") -> ", &
  ConvertFromPDG(subterm%rad(1)%ID), "(",subterm%rad(1)%i,") || ",  &
  ConvertFromPDG(subterm%rad(2)%ID), "(",subterm%rad(2)%i,") || ",  & 
  ConvertFromPDG(subterm%rad(3)%ID), "(",subterm%rad(3)%i,")  ,  ", & 
  ConvertFromPDG(subterm%rad(3)%ID),"(",subterm%rad(3)%i,") -> 0"
  call PrintBranching(iun,subterm)
!
end subroutine PrintCirsCSirs
!
subroutine PrintCirjs(iun,iterm,subterm)
use utils
implicit none
!
  integer , intent(in) :: iun,iterm
  type(subtractterm) , intent(in) :: subterm
!
!
!
  write(iun,fmt='(a,I4,6(3(a),I0),a)') &
  "iterm: ",iterm," , ", &
  ConvertFromPDG(subterm%emit(1)%ID),"(",subterm%emit(1)%i,") -> ", &
  ConvertFromPDG(subterm%rad(1)%ID),"(",subterm%rad(1)%i,") || ", &
  ConvertFromPDG(subterm%rad(2)%ID),"(",subterm%rad(2)%i,")  ,  ", & 
  ConvertFromPDG(subterm%emit(2)%ID),"(",subterm%emit(2)%i,") -> ", &
  ConvertFromPDG(subterm%rad(3)%ID),"(",subterm%rad(3)%i,") || ", & 
  ConvertFromPDG(subterm%rad(4)%ID),"(",subterm%rad(4)%i,")" 
  call PrintBranching(iun,subterm)
!
end subroutine PrintCirjs
!
subroutine PrintCirjsCSirs(iun,iterm,emitjs,radj,radjID,subterm)
use utils
implicit none
!
  integer , intent(in) :: iun,iterm
  integer , intent(in) :: emitjs,radj,radjID
  type(subtractterm) , intent(in) :: subterm
!
  integer :: emitir,radi,radr,rads
  integer :: emitirID,emitjsID,radiID,radrID,radsID
  type(subtractterm) :: subterm_tmp
!
!
  allocate(subterm_tmp%UBorn(size(subterm%UBorn(:))), &
           subterm_tmp%emit(2), &
           subterm_tmp%rad(4))
!
  subterm_tmp%UBorn = subterm%UBorn
!
  emitir   = subterm%emit(1)%i
  emitirID = subterm%emit(1)%ID
  emitjsID = subterm%UBorn(emitjs)
  radi     = subterm%rad(1)%i
  radr     = subterm%rad(2)%i
  rads     = subterm%rad(3)%i
  radiID   = subterm%rad(1)%ID
  radrID   = subterm%rad(2)%ID
  radsID   = subterm%rad(3)%ID
!
  if (emitir.lt.emitjs) then
    subterm_tmp%emit(1)    = subterm%emit(1)
    subterm_tmp%emit(2)%i  = emitjs
    subterm_tmp%emit(2)%ID = subterm%UBorn(emitjs)
    subterm_tmp%rad(1)     = subterm%rad(1)
    subterm_tmp%rad(2)     = subterm%rad(2)
    subterm_tmp%rad(4)     = subterm%rad(3)
    subterm_tmp%rad(3)%i   = radj
    subterm_tmp%rad(3)%ID  = radjID
  else
    subterm_tmp%emit(2)    = subterm%emit(1)
    subterm_tmp%emit(1)%i  = emitjs
    subterm_tmp%emit(1)%ID = subterm%UBorn(emitjs)
    subterm_tmp%rad(3)     = subterm%rad(1)
    subterm_tmp%rad(4)     = subterm%rad(2)
    subterm_tmp%rad(2)     = subterm%rad(3)
    subterm_tmp%rad(1)%i   = radj
    subterm_tmp%rad(1)%ID  = radjID
  end if
!
  write(iun,fmt='(a,I4,7(3(a),I0),a)') &
  "iterm: ",iterm," , ", &
  ConvertFromPDG(emitirID),"(",emitir,") -> ", &
  ConvertFromPDG(radiID),"(",radi,") || ", &
  ConvertFromPDG(radrID),"(",radr,")  ,  ", & 
  ConvertFromPDG(emitjsID),"(",emitjs,") -> ", &
  ConvertFromPDG(radjID),"(",radj,") || ", & 
  ConvertFromPDG(radsID),"(",rads,")  ,  ", &
  ConvertFromPDG(radsID),"(",rads,") -> 0"
  call PrintBranching(iun,subterm_tmp)
!
end subroutine PrintCirjsCSirs
!
subroutine PrintSrs(iun,iterm,subterm)
use utils
implicit none
!
  integer , intent(in) :: iun,iterm
  type(subtractterm) , intent(in) :: subterm
!
!
!
  write(iun,fmt='(a,I4,2(3(a),I0),a)') "iterm: ",iterm," , ", &
    ConvertFromPDG(subterm%rad(1)%ID),"(",subterm%rad(1)%i,") -> 0  ,  ", &
    ConvertFromPDG(subterm%rad(2)%ID),"(",subterm%rad(2)%i,") -> 0"
! We print the Underlying Born:
  write(iun,fmt='(a)',advance='no') "UBorn: "
  call PrintSubproc(iun,subterm%uborn(:))
!
end subroutine PrintSrs
!
subroutine PrintCirsSrs(iun,iterm,emitirs,radi,radiID,subterm)
use utils
implicit none
!
  integer , intent(in) :: iun,iterm
  integer , intent(in) :: emitirs,radi,radiID
  type(subtractterm) , intent(in) :: subterm
!
  integer :: radr,rads
  integer :: emitirsID,radrID,radsID
  type(subtractterm) :: subterm_tmp
!
!
  allocate(subterm_tmp%UBorn(size(subterm%UBorn(:))), &
           subterm_tmp%emit(1), &
           subterm_tmp%rad(3))
!
  subterm_tmp%UBorn = subterm%UBorn
!
  emitirsID = subterm%UBorn(emitirs)
  radr     = subterm%rad(1)%i
  rads     = subterm%rad(2)%i
  radrID   = subterm%rad(1)%ID
  radsID   = subterm%rad(2)%ID
!
  subterm_tmp%emit(1)%i  = emitirs
  subterm_tmp%emit(1)%ID = emitirsID
!
! In the new subtraction term the radiated partons have to 
! be ordered, radr and rads are ordered:
  if (radi.eq.min(radi,radr,rads)) then
    subterm_tmp%rad(1)%i   = radi
    subterm_tmp%rad(1)%ID  = radiID
    subterm_tmp%rad(2)     = subterm%rad(1)
    subterm_tmp%rad(3)     = subterm%rad(2)
  elseif (radi.eq.max(radi,radr,rads)) then
    subterm_tmp%rad(1)     = subterm%rad(1)
    subterm_tmp%rad(2)     = subterm%rad(2)
    subterm_tmp%rad(3)%i   = radi
    subterm_tmp%rad(3)%ID  = radiID
  else
    subterm_tmp%rad(1)     = subterm%rad(1)
    subterm_tmp%rad(2)%i   = radi
    subterm_tmp%rad(2)%ID  = radiID
    subterm_tmp%rad(3)     = subterm%rad(2)
  end if
!
  write(iun,fmt='(a,I4,7(3(a),I0),a)') &
  "iterm: ",iterm," , ", &
  ConvertFromPDG(subterm_tmp%emit(1)%ID),"(",subterm_tmp%emit(1)%i,") -> ", &
  ConvertFromPDG(subterm_tmp%rad(1)%ID), "(",subterm_tmp%rad(1)%i,") || ",  &
  ConvertFromPDG(subterm_tmp%rad(2)%ID), "(",subterm_tmp%rad(2)%i,") || ",  & 
  ConvertFromPDG(subterm_tmp%rad(3)%ID), "(",subterm_tmp%rad(3)%i,")  ,  ", & 
  ConvertFromPDG(radrID),"(",radr,") -> 0  ,  ", &
  ConvertFromPDG(radsID),"(",rads,") -> 0"
  call PrintBranching(iun,subterm_tmp)
!
end subroutine PrintCirsSrs
!
subroutine PrintCSirsSrs(iun,iterm,emitir,radi,radiID,radr,rads,subterm)
use utils
implicit none
!
  integer , intent(in) :: iun,iterm
  integer , intent(in) :: emitir,radi,radiID,radr,rads
  type(subtractterm) , intent(in) :: subterm
!
  integer :: emitirID,radrID,radsID
  type(subtractterm) :: subterm_tmp
!
!
  allocate(subterm_tmp%UBorn(size(subterm%UBorn(:))), &
           subterm_tmp%emit(2), &
           subterm_tmp%rad(3))
!
  subterm_tmp%UBorn = subterm%UBorn
!
  emitirID = subterm%UBorn(emitir)
! In Srs the radiated partons are ordered:
  if (radr.lt.rads) then
    radrID   = subterm%rad(1)%ID
    radsID   = subterm%rad(2)%ID
  else
    radsID   = subterm%rad(1)%ID
    radrID   = subterm%rad(2)%ID
  end if
!
  subterm_tmp%emit(1)%i  = emitir
  subterm_tmp%emit(1)%ID = emitirID
  subterm_tmp%emit(2)%i  = -1
  subterm_tmp%emit(2)%ID = 0
  subterm_tmp%rad(3)%i   = rads
  subterm_tmp%rad(3)%ID  = radsID
!
! In the new subtraction term the radiated partons have to 
! be ordered, radr and rads are ordered:
  if (radi.eq.min(radi,radr)) then
    subterm_tmp%rad(1)%i   = radi
    subterm_tmp%rad(1)%ID  = radiID
    subterm_tmp%rad(2)%i   = radr
    subterm_tmp%rad(2)%ID  = radrID
  else
    subterm_tmp%rad(2)%i   = radi
    subterm_tmp%rad(2)%ID  = radiID
    subterm_tmp%rad(1)%i   = radr
    subterm_tmp%rad(1)%ID  = radrID
  end if
!
  write(iun,fmt='(a,I4,7(3(a),I0),a)') &
  "iterm: ",iterm," , ", &
  ConvertFromPDG(subterm_tmp%emit(1)%ID),"(",subterm_tmp%emit(1)%i,") -> ", &
  ConvertFromPDG(subterm_tmp%rad(1)%ID),"(",subterm_tmp%rad(1)%i,") || ", &
  ConvertFromPDG(subterm_tmp%rad(2)%ID),"(",subterm_tmp%rad(2)%i,")  ,  ", & 
  ConvertFromPDG(radrID),"(",radr,") -> 0  ,  ",&
  ConvertFromPDG(radsID),"(",rads,") -> 0"
  call PrintBranching(iun,subterm_tmp)
!
end subroutine PrintCSirsSrs
!
subroutine PrintCirjsSrs(iun,iterm, &
                         emitir,radi,radiID, &
                         emitjs,radj,radjID,subterm)
use utils
implicit none
!
  integer , intent(in) :: iun,iterm
  integer , intent(in) :: emitir,emitjs,radi,radj,radiID,radjID
  type(subtractterm) , intent(in) :: subterm
!
  integer :: radr,rads
  integer :: emitirID,emitjsID,radrID,radsID
  type(subtractterm) :: subterm_tmp
!
!
  allocate(subterm_tmp%UBorn(size(subterm%UBorn(:))), &
           subterm_tmp%emit(2), &
           subterm_tmp%rad(4))
!
  subterm_tmp%UBorn = subterm%UBorn
!
  emitirID = subterm%UBorn(emitir)
  emitjsID = subterm%UBorn(emitjs)
  radr     = subterm%rad(1)%i
  rads     = subterm%rad(2)%i
  radrID   = subterm%rad(1)%ID
  radsID   = subterm%rad(2)%ID
!
  if (emitir.lt.emitjs) then
    subterm_tmp%emit(1)%i  = emitir
    subterm_tmp%emit(1)%ID = emitirID
    subterm_tmp%emit(2)%i  = emitjs
    subterm_tmp%emit(2)%ID = emitjsID
!
    subterm_tmp%rad(1)%i   = radi
    subterm_tmp%rad(1)%ID  = radiID
    subterm_tmp%rad(2)     = subterm%rad(1)
    subterm_tmp%rad(3)%i   = radj
    subterm_tmp%rad(3)%ID  = radjID
    subterm_tmp%rad(4)     = subterm%rad(2)
  else
    subterm_tmp%emit(2)%i  = emitir
    subterm_tmp%emit(2)%ID = emitirID
    subterm_tmp%emit(1)%i  = emitjs
    subterm_tmp%emit(1)%ID = emitjsID
!
    subterm_tmp%rad(3)%i   = radi
    subterm_tmp%rad(3)%ID  = radiID
    subterm_tmp%rad(4)     = subterm%rad(1)
    subterm_tmp%rad(1)%i   = radj
    subterm_tmp%rad(1)%ID  = radjID
    subterm_tmp%rad(2)     = subterm%rad(2)
  end if
!
  write(iun,fmt='(a,I4,8(3(a),I0),a)') &
  "iterm: ",iterm," , ", &
  ConvertFromPDG(emitirID),"(",emitir,") -> ", &
  ConvertFromPDG(radiID),"(",radi,") || ", &
  ConvertFromPDG(radrID),"(",radr,")  ,  ", & 
  ConvertFromPDG(emitjsID),"(",emitjs,") -> ", &
  ConvertFromPDG(radjID),"(",radj,") || ", & 
  ConvertFromPDG(radsID),"(",rads,")  ,  ", &
  ConvertFromPDG(radrID),"(",radr,") -> 0  ,  ", &
  ConvertFromPDG(radsID),"(",rads,") -> 0"
  call PrintBranching(iun,subterm_tmp)
!
end subroutine PrintCirjsSrs
!
subroutine PrintCirsCSirsSrs(iun,iterm,emitirs,radi,radiID, &
                             radr,rads,subterm)
use utils
implicit none
!
  integer , intent(in) :: iun,iterm
  integer , intent(in) :: emitirs,radi,radiID,radr,rads
  type(subtractterm) , intent(in) :: subterm
!
  integer :: emitirsID,radrID,radsID
  type(subtractterm) :: subterm_tmp
!
!
  allocate(subterm_tmp%UBorn(size(subterm%UBorn(:))), &
           subterm_tmp%emit(1), &
           subterm_tmp%rad(3))
!
  subterm_tmp%UBorn = subterm%UBorn
!
  emitirsID = subterm%UBorn(emitirs)
! In Srs the radiated partons are ordered:
  if (radr.lt.rads) then
    radrID   = subterm%rad(1)%ID
    radsID   = subterm%rad(2)%ID
  else
    radsID   = subterm%rad(1)%ID
    radrID   = subterm%rad(2)%ID
  end if
!
  subterm_tmp%emit(1)%i  = emitirs
  subterm_tmp%emit(1)%ID = emitirsID
!
! In the new subtraction term the radiated partons have to 
! be ordered, radr and rads are ordered:
  if (radi.eq.min(radi,radr,rads)) then
    subterm_tmp%rad(1)%i   = radi
    subterm_tmp%rad(1)%ID  = radiID
    subterm_tmp%rad(2)     = subterm%rad(1)
    subterm_tmp%rad(3)     = subterm%rad(2)
  elseif (radi.eq.max(radi,radr,rads)) then
    subterm_tmp%rad(1)     = subterm%rad(1)
    subterm_tmp%rad(2)     = subterm%rad(2)
    subterm_tmp%rad(3)%i   = radi
    subterm_tmp%rad(3)%ID  = radiID
  else
    subterm_tmp%rad(1)     = subterm%rad(1)
    subterm_tmp%rad(2)%i   = radi
    subterm_tmp%rad(2)%ID  = radiID
    subterm_tmp%rad(3)     = subterm%rad(2)
  end if
!
  write(iun,fmt='(a,I4,7(3(a),I0),a)') &
  "iterm: ",iterm," , ", &
  ConvertFromPDG(subterm_tmp%emit(1)%ID),"(",subterm_tmp%emit(1)%i,") -> ", &
  ConvertFromPDG(subterm_tmp%rad(1)%ID), "(",subterm_tmp%rad(1)%i,") || ",  &
  ConvertFromPDG(subterm_tmp%rad(2)%ID), "(",subterm_tmp%rad(2)%i,") || ",  & 
  ConvertFromPDG(subterm_tmp%rad(3)%ID), "(",subterm_tmp%rad(3)%i,")  ,  ", & 
  ConvertFromPDG(radrID),"(",radr,") -> 0  ,  ", &
  ConvertFromPDG(radsID),"(",rads,") -> 0"
  call PrintBranching(iun,subterm_tmp)
!
end subroutine PrintCirsCSirsSrs
!
subroutine PrintCktCktr(iun,iterm,subterm, &
                        radk,radkID,radt,radtID,radr,radrID)
use utils
implicit none
!
  integer , intent(in) :: iun,iterm
  type(subtractterm) , intent(in) :: subterm
  integer , intent(in) :: radk,radt,radr,radkID,radtID,radrID
!
!
!
  write(iun,fmt='(a,I4,4(3(a),I0),a)') &
  "iterm: ",iterm," , ", &
  ConvertFromPDG(subterm%emit(1)%ID),"(",subterm%emit(1)%i,") -> [[", &
  ConvertFromPDG(radkID), "(",radk,") || ", &
  ConvertFromPDG(radtID), "(",radt,")]] || ", & 
  ConvertFromPDG(radrID), "(",radr,")" 
  call PrintBranching(iun,subterm)
!
end subroutine PrintCktCktr
!
subroutine PrintCktCirkt(iun,iterm,subterm, &
                         emitir,emitirID,radi,radiID,radr,radrID, &
                         emitkt,emitktID,radk,radkID,radt,radtID)
use utils
implicit none
!
  integer , intent(in) :: iun,iterm
  type(subtractterm) , intent(in) :: subterm
  integer , intent(in) :: emitir,emitkt,emitirID,emitktID
  integer , intent(in) :: radi,radr,radk,radt, &
                          radiID,radrID,radkID,radtID
!
!
!
  write(iun,fmt='(a,I4,6(3(a),I0),a)') &
  "iterm: ",iterm," , ", &
  ConvertFromPDG(emitktID),"(",emitkt,") -> [[", &
  ConvertFromPDG(radkID), "(",radk,") || ", &
  ConvertFromPDG(radtID), "(",radt,")]]  ,  ", & 
  ConvertFromPDG(emitirID),"(",emitir,") -> ", &
  ConvertFromPDG(radiID),"(",radi,") || ", & 
  ConvertFromPDG(radrID),"(",radr,")" 
  call PrintBranching(iun,subterm)
!
end subroutine PrintCktCirkt
!
subroutine PrintCktCSktr(iun,iterm,subterm)
use utils
implicit none
!
  integer , intent(in) :: iun,iterm
  type(subtractterm) , intent(in) :: subterm
!
!
!
  write(iun,fmt='(a,I4,6(3(a),I0),a)') &
  "iterm: ",iterm," , ", &
  ConvertFromPDG(subterm%emit(1)%ID),"(",subterm%emit(1)%i,") -> [[", &
  ConvertFromPDG(subterm%rad(1)%ID),"(",subterm%rad(1)%i,") || ", &
  ConvertFromPDG(subterm%rad(2)%ID),"(",subterm%rad(2)%i,")]]  ,  ", & 
  ConvertFromPDG(subterm%rad(3)%ID),"(",subterm%rad(3)%i,") -> 0"
  call PrintBranching(iun,subterm)
!
end subroutine PrintCktCSktr
!
subroutine PrintCktCirktCSktr(iun,iterm,emitir,radi,radiID,subterm)
use utils
implicit none
!
  integer , intent(in) :: iun,iterm
  integer , intent(in) :: emitir,radi,radiID
  type(subtractterm) , intent(in) :: subterm
!
  integer :: emitkt,radk,radt,radr
  integer :: emitirID,emitktID,radkID,radtID,radrID
  type(subtractterm) :: subterm_tmp
!
!
  allocate(subterm_tmp%UBorn(size(subterm%UBorn(:))), &
           subterm_tmp%emit(2), &
           subterm_tmp%rad(4))
!
  subterm_tmp%UBorn = subterm%UBorn
!
  emitkt   = subterm%emit(1)%i
  emitktID = subterm%emit(1)%ID
  emitirID = subterm%UBorn(emitir)
  radk     = subterm%rad(1)%i
  radt     = subterm%rad(2)%i
  radr     = subterm%rad(3)%i
  radkID   = subterm%rad(1)%ID
  radtID   = subterm%rad(2)%ID
  radrID   = subterm%rad(3)%ID
!
  if (emitkt.lt.emitir) then
    subterm_tmp%emit(1)    = subterm%emit(1)
    subterm_tmp%emit(2)%i  = emitir
    subterm_tmp%emit(2)%ID = subterm%UBorn(emitir)
    subterm_tmp%rad(1)     = subterm%rad(1)
    subterm_tmp%rad(2)     = subterm%rad(2)
    subterm_tmp%rad(4)     = subterm%rad(3)
    subterm_tmp%rad(3)%i   = radi
    subterm_tmp%rad(3)%ID  = radiID
  else
    subterm_tmp%emit(2)    = subterm%emit(1)
    subterm_tmp%emit(1)%i  = emitir
    subterm_tmp%emit(1)%ID = subterm%UBorn(emitir)
    subterm_tmp%rad(3)     = subterm%rad(1)
    subterm_tmp%rad(4)     = subterm%rad(2)
    subterm_tmp%rad(2)     = subterm%rad(3)
    subterm_tmp%rad(1)%i   = radi
    subterm_tmp%rad(1)%ID  = radiID
  end if
!
  write(iun,fmt='(a,I4,7(3(a),I0),a)') &
  "iterm: ",iterm," , ", &
  ConvertFromPDG(emitirID),"(",emitir,") -> ", &
  ConvertFromPDG(radiID),"(",radi,") || ", &
  ConvertFromPDG(radrID),"(",radr,")  ,  ", & 
  ConvertFromPDG(emitktID),"(",emitkt,") -> [[", &
  ConvertFromPDG(radkID),"(",radk,") || ", & 
  ConvertFromPDG(radtID),"(",radt,")]]  ,  ", &
  ConvertFromPDG(radrID),"(",radr,") -> 0"
  call PrintBranching(iun,subterm_tmp)
!
end subroutine PrintCktCirktCSktr
!
subroutine PrintCktCktrCSktr(iun,iterm,subterm)
use utils
implicit none
!
  integer , intent(in) :: iun,iterm
  type(subtractterm) , intent(in) :: subterm
!
!
!
  write(iun,fmt='(a,I4,6(3(a),I0),a)') &
  "iterm: ",iterm," , ", &
  ConvertFromPDG(subterm%emit(1)%ID),"(",subterm%emit(1)%i,") -> [[", &
  ConvertFromPDG(subterm%rad(1)%ID), "(",subterm%rad(1)%i,") || ",  &
  ConvertFromPDG(subterm%rad(2)%ID), "(",subterm%rad(2)%i,")]] || ",  & 
  ConvertFromPDG(subterm%rad(3)%ID), "(",subterm%rad(3)%i,")  ,  ", & 
  ConvertFromPDG(subterm%rad(3)%ID),"(",subterm%rad(3)%i,") -> 0"
  call PrintBranching(iun,subterm)
!
end subroutine PrintCktCktrCSktr
!
subroutine PrintCktSkt(iun,iterm,subterm)
use utils
implicit none
!
  integer , intent(in) :: iun,iterm
  type(subtractterm) , intent(in) :: subterm
!
!
!
  write(iun,fmt='(a,I4,4(3(a),I0),a)') "iterm: ",iterm," , [[", &
    ConvertFromPDG(subterm%rad(1)%ID), "(",subterm%rad(1)%i,") || ",  &
    ConvertFromPDG(subterm%rad(2)%ID), "(",subterm%rad(2)%i,")]] , ",  & 
    ConvertFromPDG(subterm%rad(1)%ID),"(",subterm%rad(1)%i,") -> 0  ,  ", &
    ConvertFromPDG(subterm%rad(2)%ID),"(",subterm%rad(2)%i,") -> 0"
! We print the Underlying Born:
  write(iun,fmt='(a)',advance='no') "UBorn: "
  call PrintSubproc(iun,subterm%uborn(:))
!
end subroutine PrintCktSkt
!
subroutine PrintCktCrktSkt(iun,iterm,subterm,emitrkt,radr,radrID)
use utils
implicit none
!
  integer , intent(in) :: iun,iterm
  integer , intent(in) :: emitrkt,radr,radrID
  type(subtractterm) , intent(in) :: subterm
!
  integer :: radk,radt
  integer :: emitrktID,radkID,radtID
  type(subtractterm) :: subterm_tmp
!
!
  allocate(subterm_tmp%UBorn(size(subterm%UBorn(:))), &
           subterm_tmp%emit(1), &
           subterm_tmp%rad(3))
!
  subterm_tmp%UBorn = subterm%UBorn
!
  emitrktID = subterm%UBorn(emitrkt)
  radk     = subterm%rad(1)%i
  radt     = subterm%rad(2)%i
  radkID   = subterm%rad(1)%ID
  radtID   = subterm%rad(2)%ID
!
  subterm_tmp%emit(1)%i  = emitrkt
  subterm_tmp%emit(1)%ID = emitrktID
!
! In the new subtraction term the radiated partons have to 
! be ordered, radr and rads are ordered:
  if (radr.eq.min(radr,radk,radt)) then
    subterm_tmp%rad(1)%i   = radr
    subterm_tmp%rad(1)%ID  = radrID
    subterm_tmp%rad(2)     = subterm%rad(1)
    subterm_tmp%rad(3)     = subterm%rad(2)
  elseif (radr.eq.max(radr,radk,radt)) then
    subterm_tmp%rad(1)     = subterm%rad(1)
    subterm_tmp%rad(2)     = subterm%rad(2)
    subterm_tmp%rad(3)%i   = radr
    subterm_tmp%rad(3)%ID  = radrID
  else
    subterm_tmp%rad(1)     = subterm%rad(1)
    subterm_tmp%rad(2)%i   = radr
    subterm_tmp%rad(2)%ID  = radrID
    subterm_tmp%rad(3)     = subterm%rad(2)
  end if
!
  write(iun,fmt='(a,I4,7(3(a),I0),a)') &
  "iterm: ",iterm," , ", &
  ConvertFromPDG(emitrktID),"(",emitrkt,") -> ", &
  ConvertFromPDG(radrID), "(",radr,") || [[",  &
  ConvertFromPDG(radkID), "(",radk,") || ",  & 
  ConvertFromPDG(radtID), "(",radt,")]]  ,  ", & 
  ConvertFromPDG(radkID),"(",radk,") -> 0  ,  ", &
  ConvertFromPDG(radtID),"(",radt,") -> 0"
  call PrintBranching(iun,subterm_tmp)
!
end subroutine PrintCktCrktSkt
!
subroutine PrintStCirt(iun,iterm,subterm, &
                       radi,radiID,radr,radrID,radt,radtID)
use utils
implicit none
!
  integer , intent(in) :: iun,iterm
  type(subtractterm) , intent(in) :: subterm
  integer , intent(in) :: radi,radr,radt,radiID,radrID,radtID
!
!
!
  write(iun,fmt='(a,I4,4(3(a),I0),a)') &
  "iterm: ",iterm," , ", &
  ConvertFromPDG(subterm%emit(1)%ID),"(",subterm%emit(1)%i,") -> ", &
  ConvertFromPDG(radiID), "(",radi,") || ", &
  ConvertFromPDG(radrID), "(",radr,") || [[", & 
  ConvertFromPDG(radtID), "(",radt,")]]" 
  call PrintBranching(iun,subterm)
!
end subroutine PrintStCirt
!
subroutine PrintStCSirt(iun,iterm,subterm)
use utils
implicit none
!
  integer , intent(in) :: iun,iterm
  type(subtractterm) , intent(in) :: subterm
!
!
!
  write(iun,fmt='(a,I4,6(3(a),I0),a)') &
  "iterm: ",iterm," , ", &
  ConvertFromPDG(subterm%emit(1)%ID),"(",subterm%emit(1)%i,") -> ", &
  ConvertFromPDG(subterm%rad(1)%ID),"(",subterm%rad(1)%i,") || ", &
  ConvertFromPDG(subterm%rad(2)%ID),"(",subterm%rad(2)%i,")  ,  [[", & 
  ConvertFromPDG(subterm%rad(3)%ID),"(",subterm%rad(3)%i,") -> 0]]"
  call PrintBranching(iun,subterm)
!
end subroutine PrintStCSirt
!
subroutine PrintStCirtCSirt(iun,iterm,subterm)
use utils
implicit none
!
  integer , intent(in) :: iun,iterm
  type(subtractterm) , intent(in) :: subterm
!
!
!
  write(iun,fmt='(a,I4,6(3(a),I0),a)') &
  "iterm: ",iterm," , ", &
  ConvertFromPDG(subterm%emit(1)%ID),"(",subterm%emit(1)%i,") -> ", &
  ConvertFromPDG(subterm%rad(1)%ID), "(",subterm%rad(1)%i,") || ",  &
  ConvertFromPDG(subterm%rad(2)%ID), "(",subterm%rad(2)%i,") || ",  & 
  ConvertFromPDG(subterm%rad(3)%ID), "(",subterm%rad(3)%i,")  ,  [[", & 
  ConvertFromPDG(subterm%rad(3)%ID),"(",subterm%rad(3)%i,") -> 0]]"
  call PrintBranching(iun,subterm)
!
end subroutine PrintStCirtCSirt
!
subroutine PrintStCirtSrt(iun,iterm,subterm,emitirt, &
                          radi,radiID,radr,radrID,radt,radtID)
use utils
implicit none
!
  integer , intent(in) :: iun,iterm
  integer , intent(in) :: emitirt,radi,radiID,radr,radrID,radt,radtID
  type(subtractterm) , intent(in) :: subterm
!
  integer :: emitirtID
  type(subtractterm) :: subterm_tmp
!
!
  allocate(subterm_tmp%UBorn(size(subterm%UBorn(:))), &
           subterm_tmp%emit(1), &
           subterm_tmp%rad(3))
!
  subterm_tmp%UBorn = subterm%UBorn
!
  emitirtID = subterm%UBorn(emitirt)
!
  subterm_tmp%emit(1)%i  = emitirt
  subterm_tmp%emit(1)%ID = emitirtID
!
! In the new subtraction term the radiated partons have to 
! be ordered, radr and rads are ordered:
  if (radi.eq.min(radi,radr,radt)) then
    subterm_tmp%rad(1)%i   = radi
    subterm_tmp%rad(1)%ID  = radiID
    subterm_tmp%rad(2)     = subterm%rad(1)
    subterm_tmp%rad(3)     = subterm%rad(2)
  elseif (radi.eq.max(radi,radr,radt)) then
    subterm_tmp%rad(1)     = subterm%rad(1)
    subterm_tmp%rad(2)     = subterm%rad(2)
    subterm_tmp%rad(3)%i   = radi
    subterm_tmp%rad(3)%ID  = radiID
  else
    subterm_tmp%rad(1)     = subterm%rad(1)
    subterm_tmp%rad(2)%i   = radi
    subterm_tmp%rad(2)%ID  = radiID
    subterm_tmp%rad(3)     = subterm%rad(2)
  end if
!
  write(iun,fmt='(a,I4,7(3(a),I0),a)') &
  "iterm: ",iterm," , ", &
  ConvertFromPDG(emitirtID),"(",emitirt,") -> ", &
  ConvertFromPDG(radiID), "(",radi,") || ",  &
  ConvertFromPDG(radrID), "(",radr,") || ",  & 
  ConvertFromPDG(radtID), "(",radt,")  ,  ", & 
  ConvertFromPDG(radrID),"(",radr,") -> 0  ,  [[", &
  ConvertFromPDG(radtID),"(",radt,") -> 0]]"
  call PrintBranching(iun,subterm_tmp)
!
end subroutine PrintStCirtSrt
!
subroutine PrintStCSirtSrt(iun,iterm,subterm,emitir,radi,radiID, &
                           radr,radrID,radt,radtID)
use utils
implicit none
!
  integer , intent(in) :: iun,iterm
  integer , intent(in) :: emitir,radi,radiID,radr,radrID,radt,radtID
  type(subtractterm) , intent(in) :: subterm
!
  integer :: emitirID
  type(subtractterm) :: subterm_tmp
!
!
  allocate(subterm_tmp%UBorn(size(subterm%UBorn(:))), &
           subterm_tmp%emit(2), &
           subterm_tmp%rad(3))
!
  subterm_tmp%UBorn = subterm%UBorn
!
  emitirID = subterm%UBorn(emitir)
!
  subterm_tmp%emit(1)%i  = emitir
  subterm_tmp%emit(1)%ID = emitirID
  subterm_tmp%emit(2)%i  = -1
  subterm_tmp%emit(2)%ID = 0
  subterm_tmp%rad(3)%i   = radt
  subterm_tmp%rad(3)%ID  = radtID
!
! In the new subtraction term the radiated partons have to 
! be ordered, radi and radr are ordered:
  if (radi.eq.min(radi,radr)) then
    subterm_tmp%rad(1)%i   = radi
    subterm_tmp%rad(1)%ID  = radiID
    subterm_tmp%rad(2)%i   = radr
    subterm_tmp%rad(2)%ID  = radrID
  else
    subterm_tmp%rad(2)%i   = radi
    subterm_tmp%rad(2)%ID  = radiID
    subterm_tmp%rad(1)%i   = radr
    subterm_tmp%rad(1)%ID  = radrID
  end if
!
  write(iun,fmt='(a,I4,7(3(a),I0),a)') &
  "iterm: ",iterm," , ", &
  ConvertFromPDG(subterm_tmp%emit(1)%ID),"(",subterm_tmp%emit(1)%i,") -> ", &
  ConvertFromPDG(subterm_tmp%rad(1)%ID),"(",subterm_tmp%rad(1)%i,") || ", &
  ConvertFromPDG(subterm_tmp%rad(2)%ID),"(",subterm_tmp%rad(2)%i,")  ,  ", & 
  ConvertFromPDG(radrID),"(",radr,") -> 0  ,  [[",&
  ConvertFromPDG(radtID),"(",radt,") -> 0]]"
  call PrintBranching(iun,subterm_tmp)
!
end subroutine PrintStCSirtSrt
!
subroutine PrintStCirtCSirtSrt(iun,iterm,emitirt,radi,radiID, &
                               radr,radt,subterm)
use utils
implicit none
!
  integer , intent(in) :: iun,iterm
  integer , intent(in) :: emitirt,radi,radiID,radr,radt
  type(subtractterm) , intent(in) :: subterm
!
  integer :: emitirtID,radrID,radtID
  type(subtractterm) :: subterm_tmp
!
!
  allocate(subterm_tmp%UBorn(size(subterm%UBorn(:))), &
           subterm_tmp%emit(1), &
           subterm_tmp%rad(3))
!
  subterm_tmp%UBorn = subterm%UBorn
!
  emitirtID = subterm%UBorn(emitirt)
! In Srs the radiated partons are ordered:
  if (radr.lt.radt) then
    radrID   = subterm%rad(1)%ID
    radtID   = subterm%rad(2)%ID
  else
    radtID   = subterm%rad(1)%ID
    radrID   = subterm%rad(2)%ID
  end if
!
  subterm_tmp%emit(1)%i  = emitirt
  subterm_tmp%emit(1)%ID = emitirtID
!
! In the new subtraction term the radiated partons have to 
! be ordered, radr and rads are ordered:
  if (radi.eq.min(radi,radr,radt)) then
    subterm_tmp%rad(1)%i   = radi
    subterm_tmp%rad(1)%ID  = radiID
    subterm_tmp%rad(2)     = subterm%rad(1)
    subterm_tmp%rad(3)     = subterm%rad(2)
  elseif (radi.eq.max(radi,radr,radt)) then
    subterm_tmp%rad(1)     = subterm%rad(1)
    subterm_tmp%rad(2)     = subterm%rad(2)
    subterm_tmp%rad(3)%i   = radi
    subterm_tmp%rad(3)%ID  = radiID
  else
    subterm_tmp%rad(1)     = subterm%rad(1)
    subterm_tmp%rad(2)%i   = radi
    subterm_tmp%rad(2)%ID  = radiID
    subterm_tmp%rad(3)     = subterm%rad(2)
  end if
!
  write(iun,fmt='(a,I4,7(3(a),I0),a)') &
  "iterm: ",iterm," , ", &
  ConvertFromPDG(subterm_tmp%emit(1)%ID),"(",subterm_tmp%emit(1)%i,") -> ", &
  ConvertFromPDG(subterm_tmp%rad(1)%ID), "(",subterm_tmp%rad(1)%i,") || ",  &
  ConvertFromPDG(subterm_tmp%rad(2)%ID), "(",subterm_tmp%rad(2)%i,") || ",  & 
  ConvertFromPDG(subterm_tmp%rad(3)%ID), "(",subterm_tmp%rad(3)%i,")  ,  ", & 
  ConvertFromPDG(radrID),"(",radr,") -> 0  ,  [[", &
  ConvertFromPDG(radtID),"(",radt,") -> 0]]"
  call PrintBranching(iun,subterm_tmp)
!
end subroutine PrintStCirtCSirtSrt
!
subroutine PrintStSrt(iun,iterm,radr,radt,subterm)
use utils
implicit none
!
  integer , intent(in) :: iun,iterm
  integer , intent(in) :: radr,radt
  type(subtractterm) , intent(in) :: subterm
!
  integer :: radrID,radtID
!
!
  if (radr.eq.subterm%rad(1)%i) then
    radrID = subterm%rad(1)%ID
    radtID = subterm%rad(2)%ID
  else
    radrID = subterm%rad(2)%ID
    radtID = subterm%rad(1)%ID
  end if
!
  write(iun,fmt='(a,I4,2(3(a),I0),a)') "iterm: ",iterm," , ", &
    ConvertFromPDG(radrID),"(",radr,") -> 0  ,  [[", &
    ConvertFromPDG(radtID),"(",radt,") -> 0]]"
! We print the Underlying Born:
  write(iun,fmt='(a)',advance='no') "UBorn: "
  call PrintSubproc(iun,subterm%uborn(:))
!
end subroutine PrintStSrt
!
subroutine PrintCitStCirt(iun,iterm,subterm, &
                          radi,radiID,radr,radrID,radt,radtID)
use utils
implicit none
!
  integer , intent(in) :: iun,iterm
  type(subtractterm) , intent(in) :: subterm
  integer , intent(in) :: radi,radr,radt,radiID,radrID,radtID
!
!
!
  write(iun,fmt='(a,I4,7(3(a),I0),a)') &
  "iterm: ",iterm," , ", &
  ConvertFromPDG(subterm%emit(1)%ID),"(",subterm%emit(1)%i,") -> ", &
  ConvertFromPDG(subterm%rad(1)%ID), "(",subterm%rad(1)%i,") || ", &
  ConvertFromPDG(subterm%rad(2)%ID), "(",subterm%rad(2)%i,") || ", & 
  ConvertFromPDG(subterm%rad(3)%ID), "(",subterm%rad(3)%i,") , [[", &
  ConvertFromPDG(radiID), "(",radi,") || ", &
  ConvertFromPDG(radtID), "(",radt,") , ", & 
  ConvertFromPDG(radtID), "(",radt,") -> 0]]" 
  call PrintBranching(iun,subterm)
!
end subroutine PrintCitStCirt
!
subroutine PrintCktStCSirt(iun,iterm,subterm, &
                           radk,radkID,radt,radtID)
use utils
implicit none
!
  integer , intent(in) :: iun,iterm
  type(subtractterm) , intent(in) :: subterm
  integer , intent(in) :: radk,radt,radkID,radtID
!
!
!
  write(iun,fmt='(a,I4,7(3(a),I0),a)') &
  "iterm: ",iterm," , ", &
  ConvertFromPDG(subterm%emit(1)%ID),"(",subterm%emit(1)%i,") -> ", &
  ConvertFromPDG(subterm%rad(1)%ID), "(",subterm%rad(1)%i,") || ", &
  ConvertFromPDG(subterm%rad(2)%ID), "(",subterm%rad(2)%i,") , ", & 
  ConvertFromPDG(subterm%rad(3)%ID), "(",subterm%rad(3)%i,") -> 0 , [[", &
  ConvertFromPDG(radkID), "(",radk,") || ", &
  ConvertFromPDG(radtID), "(",radt,") , ", & 
  ConvertFromPDG(radtID), "(",radt,") -> 0]]" 
  call PrintBranching(iun,subterm)
!
end subroutine PrintCktStCSirt
!
subroutine PrintCktStCSirtSrt(iun,iterm,subterm,emitir,radi,radiID, &
                              radr,radrID,radk,radkID,radt,radtID)
use utils
implicit none
!
  integer , intent(in) :: iun,iterm
  integer , intent(in) :: emitir
  integer , intent(in) :: radi,radiID,radr,radrID,radk,radkID,radt,radtID
  type(subtractterm) , intent(in) :: subterm
!
  integer :: emitirID
  type(subtractterm) :: subterm_tmp
!
!
  allocate(subterm_tmp%UBorn(size(subterm%UBorn(:))), &
           subterm_tmp%emit(2), &
           subterm_tmp%rad(3))
!
  subterm_tmp%UBorn = subterm%UBorn
!
  emitirID = subterm%UBorn(emitir)
!
  subterm_tmp%emit(1)%i  = emitir
  subterm_tmp%emit(1)%ID = emitirID
  subterm_tmp%emit(2)%i  = -1
  subterm_tmp%emit(2)%ID = 0
  subterm_tmp%rad(3)%i   = radt
  subterm_tmp%rad(3)%ID  = radtID
!
! In the new subtraction term the radiated partons have to 
! be ordered, radi and radr are ordered:
  if (radi.eq.min(radi,radr)) then
    subterm_tmp%rad(1)%i   = radi
    subterm_tmp%rad(1)%ID  = radiID
    subterm_tmp%rad(2)%i   = radr
    subterm_tmp%rad(2)%ID  = radrID
  else
    subterm_tmp%rad(2)%i   = radi
    subterm_tmp%rad(2)%ID  = radiID
    subterm_tmp%rad(1)%i   = radr
    subterm_tmp%rad(1)%ID  = radrID
  end if
!
  write(iun,fmt='(a,I4,10(3(a),I0),a)') &
  "iterm: ",iterm," , ", &
  ConvertFromPDG(subterm_tmp%emit(1)%ID),"(",subterm_tmp%emit(1)%i,") -> ", &
  ConvertFromPDG(subterm_tmp%rad(1)%ID),"(",subterm_tmp%rad(1)%i,") || ", &
  ConvertFromPDG(subterm_tmp%rad(2)%ID),"(",subterm_tmp%rad(2)%i,")  ,  ", & 
  ConvertFromPDG(subterm%rad(1)%ID),"(",subterm%rad(1)%i,") -> 0  ,  ", &
  ConvertFromPDG(subterm%rad(2)%ID),"(",subterm%rad(2)%i,") -> 0 , [[", &
  ConvertFromPDG(radkID), "(",radk,") || ", &
  ConvertFromPDG(radtID), "(",radt,") , ", & 
  ConvertFromPDG(radtID), "(",radt,") -> 0]]" 
  call PrintBranching(iun,subterm_tmp)
!
end subroutine PrintCktStCSirtSrt
!
subroutine PrintCktStCkrtSrt(iun,iterm,subterm,emitkrt, &
                             radk,radkID,radr,radrID,radt,radtID)
use utils
implicit none
!
  integer , intent(in) :: iun,iterm
  integer , intent(in) :: emitkrt,radk,radkID,radr,radrID,radt,radtID
  type(subtractterm) , intent(in) :: subterm
!
  integer :: emitkrtID
  type(subtractterm) :: subterm_tmp
!
!
  allocate(subterm_tmp%UBorn(size(subterm%UBorn(:))), &
           subterm_tmp%emit(1), &
           subterm_tmp%rad(3))
!
  subterm_tmp%UBorn = subterm%UBorn
!
  emitkrtID = subterm%UBorn(emitkrt)
!
  subterm_tmp%emit(1)%i  = emitkrt
  subterm_tmp%emit(1)%ID = emitkrtID
!
! In the new subtraction term the radiated partons have to 
! be ordered, radr and radt are ordered:
  if (radk.eq.min(radk,radr,radt)) then
    subterm_tmp%rad(1)%i   = radk
    subterm_tmp%rad(1)%ID  = radkID
    subterm_tmp%rad(2)     = subterm%rad(1)
    subterm_tmp%rad(3)     = subterm%rad(2)
  elseif (radk.eq.max(radk,radr,radt)) then
    subterm_tmp%rad(1)     = subterm%rad(1)
    subterm_tmp%rad(2)     = subterm%rad(2)
    subterm_tmp%rad(3)%i   = radk
    subterm_tmp%rad(3)%ID  = radkID
  else
    subterm_tmp%rad(1)     = subterm%rad(1)
    subterm_tmp%rad(2)%i   = radk
    subterm_tmp%rad(2)%ID  = radkID
    subterm_tmp%rad(3)     = subterm%rad(2)
  end if
!
  write(iun,fmt='(a,I4,9(3(a),I0),a)') &
  "iterm: ",iterm," , ", &
  ConvertFromPDG(emitkrtID),"(",emitkrt,") -> ", &
  ConvertFromPDG(radkID), "(",radk,") || ",  &
  ConvertFromPDG(radrID), "(",radr,") || ",  & 
  ConvertFromPDG(radtID), "(",radt,")  ,  ", & 
  ConvertFromPDG(subterm%rad(1)%ID),"(",subterm%rad(1)%i,") -> 0  ,  ", &
  ConvertFromPDG(subterm%rad(2)%ID),"(",subterm%rad(2)%i,") -> 0 , [[", &
  ConvertFromPDG(radkID), "(",radk,") || ", &
  ConvertFromPDG(radtID), "(",radt,") , ", & 
  ConvertFromPDG(radtID), "(",radt,") -> 0]]" 
  call PrintBranching(iun,subterm_tmp)
!
end subroutine PrintCktStCkrtSrt
!
subroutine PrintCrtStCkrtSrt(iun,iterm,subterm,emitkrt, &
                             radk,radkID,radr,radrID,radt,radtID)
use utils
implicit none
!
  integer , intent(in) :: iun,iterm
  integer , intent(in) :: emitkrt,radk,radkID,radr,radrID,radt,radtID
  type(subtractterm) , intent(in) :: subterm
!
  integer :: emitkrtID
  type(subtractterm) :: subterm_tmp
!
!
  allocate(subterm_tmp%UBorn(size(subterm%UBorn(:))), &
           subterm_tmp%emit(1), &
           subterm_tmp%rad(3))
!
  subterm_tmp%UBorn = subterm%UBorn
!
  emitkrtID = subterm%UBorn(emitkrt)
!
  subterm_tmp%emit(1)%i  = emitkrt
  subterm_tmp%emit(1)%ID = emitkrtID
!
! In the new subtraction term the radiated partons have to 
! be ordered, radr and radt are ordered:
  if (radk.eq.min(radk,radr,radt)) then
    subterm_tmp%rad(1)%i   = radk
    subterm_tmp%rad(1)%ID  = radkID
    subterm_tmp%rad(2)     = subterm%rad(1)
    subterm_tmp%rad(3)     = subterm%rad(2)
  elseif (radk.eq.max(radk,radr,radt)) then
    subterm_tmp%rad(1)     = subterm%rad(1)
    subterm_tmp%rad(2)     = subterm%rad(2)
    subterm_tmp%rad(3)%i   = radk
    subterm_tmp%rad(3)%ID  = radkID
  else
    subterm_tmp%rad(1)     = subterm%rad(1)
    subterm_tmp%rad(2)%i   = radk
    subterm_tmp%rad(2)%ID  = radkID
    subterm_tmp%rad(3)     = subterm%rad(2)
  end if
!
  write(iun,fmt='(a,I4,9(3(a),I0),a)') &
  "iterm: ",iterm," , ", &
  ConvertFromPDG(emitkrtID),"(",emitkrt,") -> ", &
  ConvertFromPDG(radkID), "(",radk,") || ",  &
  ConvertFromPDG(radrID), "(",radr,") || ",  & 
  ConvertFromPDG(radtID), "(",radt,")  ,  ", & 
  ConvertFromPDG(subterm%rad(1)%ID),"(",subterm%rad(1)%i,") -> 0  ,  ", &
  ConvertFromPDG(subterm%rad(2)%ID),"(",subterm%rad(2)%i,") -> 0 , [[", &
  ConvertFromPDG(radrID), "(",radr,") || ", &
  ConvertFromPDG(radtID), "(",radt,") , ", & 
  ConvertFromPDG(radtID), "(",radt,") -> 0]]" 
  call PrintBranching(iun,subterm_tmp)
!
end subroutine PrintCrtStCkrtSrt
!
subroutine PrintCktStSrt(iun,iterm,subterm, &
                         radk,radkID,radr,radrID,radt,radtID)
use utils
implicit none
!
  integer , intent(in) :: iun,iterm
  integer , intent(in) :: radk,radkID,radr,radrID,radt,radtID
  type(subtractterm) , intent(in) :: subterm
!
!
!
!
  write(iun,fmt='(a,I4,5(3(a),I0),a)') &
  "iterm: ",iterm," , ", &
  ConvertFromPDG(subterm%rad(1)%ID),"(",subterm%rad(1)%i,") -> 0  ,  ", &
  ConvertFromPDG(subterm%rad(2)%ID),"(",subterm%rad(2)%i,") -> 0 , [[", &
  ConvertFromPDG(radkID), "(",radk,") || ", &
  ConvertFromPDG(radtID), "(",radt,") , ", & 
  ConvertFromPDG(radtID), "(",radt,") -> 0]]" 
  call PrintBranching(iun,subterm)
!
end subroutine PrintCktStSrt
!
subroutine PrintCrtStSrt(iun,iterm,subterm, &
                         radr,radrID,radt,radtID)
use utils
implicit none
!
  integer , intent(in) :: iun,iterm
  integer , intent(in) :: radr,radrID,radt,radtID
  type(subtractterm) , intent(in) :: subterm
!
!
!
!
  write(iun,fmt='(a,I4,5(3(a),I0),a)') &
  "iterm: ",iterm," , ", &
  ConvertFromPDG(subterm%rad(1)%ID),"(",subterm%rad(1)%i,") -> 0  ,  ", &
  ConvertFromPDG(subterm%rad(2)%ID),"(",subterm%rad(2)%i,") -> 0 , [[", &
  ConvertFromPDG(radrID), "(",radr,") || ", &
  ConvertFromPDG(radtID), "(",radt,") , ", & 
  ConvertFromPDG(radtID), "(",radt,") -> 0]]" 
  call PrintBranching(iun,subterm)
!
end subroutine PrintCrtStSrt
!
end module 
!
! This routine goes through all the real real-virtual and
! double real contributions and detects the singular
! regions:
subroutine init_Regions
use process
use flags
use regions
use subprocesses
use input
implicit none
!
!
  integer :: iproc,nproc
  integer :: istat
!
! We set up the parameters:
  sub_d0 = 3
  sub_d0pr = 3
! We also call init_procsubtract to have the possibility to 
! alter the predefined values:
  call init_procsubtract()
!
  write(*,*) "*********************************************************"
  write(*,*) "*                                                       *"
  write(*,"(1x,a,I0,a)") "*                 d0 is:     ",sub_d0, &
          "                          *"
  write(*,"(1x,a,I0,a)") "*                 d0pr is:   ",sub_d0, &
          "                          *"
  write(*,*) "*                                                       *"
  write(*,*) "*********************************************************"
!
! We allocate arrays to hold momenta:
  allocate(psubB(1:nleg_born),    &
           psubR(1:nleg_born+1),  &
           psubRR(1:nleg_born+2), &
           stat=istat)
  if (istat.ne.0) then
    print *,"Allocational problem for psub arrays..."
    stop
  end if
! We allocate the arrays holding color-correlated MEs:
  allocate(subBij(nleg_born,nleg_born), &
           subRij(nleg_born+1,nleg_born+1), &
           stat=istat)
  if (istat.ne.0) then
    print *,"Allocation not working for sub[B,R]ij..."
    stop
  end if
! We also allocate the invariants:
  allocate(yiQ_arr(nleg_born+2),  &
           sir_arr(nleg_born+2,nleg_born+2),  &
           yir_arr(nleg_born+2,nleg_born+2),  &
           stat=istat)
  if (istat.ne.0) then
    print *,"Allocational problem for subtraction related invariants..."
    stop
  end if
! Allocate tilde invariants only if they are really needed:
  if (flg_NNLO_RRA1_A1) then
    allocate(yitildeQ_arr(nleg_born), &
             sirtilde_arr(nleg_born,nleg_born),  &
             yirtilde_arr(nleg_born,nleg_born),  &
             stat=istat)
  end if
  if (istat.ne.0) then 
    print *,"Problem with allocating the tilde invariants..."
    stop
  end if
! Allocate hat invariants only if they are really needed:
  if (flg_NNLO_RR_A12) then
    allocate(yihatQ_arr(nleg_born+1), &
             sirhat_arr(nleg_born+1,nleg_born+1),  &
             stat=istat)
  end if
  if (istat.ne.0) then 
    print *,"Problem with allocating the hat invariants..."
    stop
  end if
!
! if debugscsme is present in the input card we check the 
! spin-correlated squared matrix elements.
  if (nnloinput("#debugscsme").gt.0) call DebugSCSME
! if debugccsme is present we check the color-correlated MEs:
  if (nnloinput("#debugccsme").gt.0) call DebugCCSME
! if debugsccsme is present we check the simultaneously spin- and 
! color-correlated MEs:
  if (nnloinput("#debugsccsme").gt.0) call DebugSCCSME
! if debugdccsme is present we check the doubly color-correlated MEs:
  if (nnloinput("#debugdccsme").gt.0) call DebugDCCSME
!
! In the following we obtain all the singular regions relevant for 
! NLO-type subtractions, for this type of regions subtractions are
! assigned in three different contributions: Real(R), Real-Virtual(RV)
! and Double-Real(RR):
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!///////////////////////////////////////////////////////////////////////
! NLO-type subtractions:
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!***********************************************************************
! Real contribution:
!***********************************************************************
!
! We not only need the subtraction terms for the Real contribution if
! it is required to be calculated but when subtraction terms are needed
! for the integrated A_1-type subtraction terms associated to the 
! double-real contribution:
  if (flg_NLO_R.or.flg_NNLO_RV_A1.or.flg_NNLO_RRA1_A1) then
! We allocate the arrays which will hold the subtraction terms:
    nproc = num_flv_irr_NLO_R
    if (.not.flg_NLO_R.and.flg_NNLO_RV) nproc = num_flv_irr_NNLO_RV
    allocate(subterms_Cir_R(nproc),  &
             subterms_Sr_R(nproc),   &
             subterms_CSir_R(nproc), &
             stat=istat)
    if (istat.ne.0) then
      print *,"Problem occured during allocating subterms_*_R..."
      stop
    end if
! We go through all the subprocesses at the real emission level
! and try to isolate NLO-type subtraction terms (C,S,CS):
    do iproc=1,nproc
      call FindSubTermsNLO(num_flv_LO,flv_LO,flv_ch_Rkin(:,iproc), &
                           subterms_Cir_R(iproc),subterms_Sr_R(iproc), &
                           subterms_CSir_R(iproc))
!
    end do
  end if
!
!***********************************************************************
! Double-real contribution:
!***********************************************************************
!
  if (flg_NNLO_RR_A1.or.flg_NNLO_RR_A12) then
! We allocate the arrays which will hold the subtraction terms:
    allocate(subterms_Cir_RR(num_flv_irr_NNLO_RR),  &
             subterms_Sr_RR(num_flv_irr_NNLO_RR),   &
             subterms_CSir_RR(num_flv_irr_NNLO_RR), &
             stat=istat)
    if (istat.ne.0) then
      print *,"Problem occured during allocating subterms_*_RR..."
      stop
    end if
! We go through all the subprocesses at the double-real emission level
! and try to isolate NLO-type subtraction terms (C,S,CS):
    do iproc=1,num_flv_irr_NNLO_RR
      call FindSubTermsNLO(num_flv_NLO_R,flv_NLO_R,flv_ch_RRkin(:,iproc), &
                           subterms_Cir_RR(iproc),subterms_Sr_RR(iproc),  &
                           subterms_CSir_RR(iproc))
!
    end do
  end if
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!///////////////////////////////////////////////////////////////////////
! NNLO-type subtractions:
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!***********************************************************************
! Real-virtual contribution:
!***********************************************************************
! Note that the kinematically degenerate regions of RV correspond to
! those of the R, hence no new array is needed.
!
!***********************************************************************
! Double-Real contribution:
!***********************************************************************
!
! A2 subtraction terms are needed not only when A2 are requested 
! but when A12 terms are needed to be calculated too. Since
! for the A12 terms no separate subtraction terms are defined
! and stored in memory but determined on the fly.
  if (flg_NNLO_RR_A2.or.flg_NNLO_RR_A12) then
! We allocate the arrays which will hold the subtraction terms:
    allocate(subterms_Cirs_RR(num_flv_irr_NNLO_RR),  &
             subterms_Cirjs_RR(num_flv_irr_NNLO_RR), &
             subterms_CSirs_RR(num_flv_irr_NNLO_RR), &
             subterms_Srs_RR(num_flv_irr_NNLO_RR),   &
             stat=istat)
    if (istat.ne.0) then
      print *,"Problem occured during allocating subterms_*_RR..."
      stop
    end if
! We go through all the subprocesses at the double-real emission level
! and try to isolate NNLO-type subtraction terms:
    do iproc=1,num_flv_irr_NNLO_RR
      call FindSubTermsNNLO(num_flv_LO,flv_LO,flv_ch_RRkin(:,iproc), &
                            subterms_Cirs_RR(iproc),       &
                            subterms_Cirjs_RR(iproc),      &
                            subterms_CSirs_RR(iproc),      &
                            subterms_Srs_RR(iproc)         &
                           )
!
    end do
  end if
!
! We print out all the Subtraction terms:
  if (flg_NLO_R.or. &
      flg_NNLO_RR.or. &
      flg_NNLO_RV) then
    call PrintSubTerms
  end if
!  stop "init_Regions"
!
! We test the limits, too:
  if (flg_NLO_R.or. &
      flg_NNLO_RR.or. &
      flg_NNLO_RV) then
! We offer the possibility to make scatter plots: 
!    if (nnloinput("#scatterplot").gt.0) call MakeScatterPlot
    if (nnloinput("#scatterplot").gt.0) call MakeRatioPlot
! Checking the various limits and writting the various ratios
! to a file:
    call TestLimits
! We can also check the regions in an entire line:
    if (nnloinput("#regioncheck").gt.0) call RegionCheck
  end if
!
! It is possible to check the cancellation of poles among
! subtractions for the RV and RR,A1:
  if (flg_NNLO_RV.and.(nnloinput("#debugpoles").gt.0)) then
    call CheckSubtractPoles
    stop
  end if
!
end subroutine init_Regions
!
subroutine OrderRadir(radi,radiID,radr,radrID,rad1,rad1ID,rad2,rad2ID)
implicit none
!
  integer , intent(in) :: radi,radr
  integer , intent(in) :: radiID,radrID
  integer , intent(out) :: rad1,rad2
  integer , intent(out) :: rad1ID,rad2ID
!
!
!
! If having gluons the ordering goes with the position:
  if (min(radiID,radrID).eq.max(radiID,radrID)) then
    if (min(radi,radr).eq.radi) then
      rad1 = radi
      rad2 = radr
      rad1ID = radiID
      rad2ID = radrID
    else
      rad1 = radr
      rad2 = radi
      rad1ID = radrID
      rad2ID = radiID
    end if
! Otherwise the quark always comes first:
! q(~) + g splitting:
  elseif (min(abs(radiID),abs(radrID)).eq.0) then
    if (abs(radiID).gt.0) then
      rad1 = radi
      rad2 = radr
      rad1ID = radiID
      rad2ID = radrID
    else
      rad1 = radr
      rad2 = radi
      rad1ID = radrID
      rad2ID = radiID
    end if
! q + q(~) splitting:
  elseif (radiID.eq.-radrID) then
    if (radiID.gt.0) then
      rad1 = radi
      rad2 = radr
      rad1ID = radiID
      rad2ID = radrID
    else
      rad1 = radr
      rad2 = radi
      rad1ID = radrID
      rad2ID = radiID
    end if
  end if
!
end subroutine OrderRadir
!
subroutine ReorderRadir(radi,radiID,radr,radrID,rad1,rad2)
implicit none
!
  integer , intent(in) :: radi,radr
  integer , intent(in) :: radiID,radrID
  integer , intent(out) :: rad1,rad2
!
!
!
! If having gluons the ordering goes with the position:
  if (min(radiID,radrID).eq.max(radiID,radrID)) then
    if (min(radi,radr).eq.radi) then
      rad1 = radi
      rad2 = radr
    else
      rad1 = radr
      rad2 = radi
    end if
! Otherwise the quark always comes first:
! q(~) + g splitting:
  elseif (min(abs(radiID),abs(radrID)).eq.0) then
    if (abs(radiID).gt.0) then
      rad1 = radi
      rad2 = radr
    else
      rad1 = radr
      rad2 = radi
    end if
! q + q(~) splitting:
  elseif (radiID.eq.-radrID) then
    if (radiID.gt.0) then
      rad1 = radi
      rad2 = radr
    else
      rad1 = radr
      rad2 = radi
    end if
  end if
!
end subroutine ReorderRadir
!
! This routine takes three radiated partons: rad1,rad2,rad3
! and reorder them to follow the convention used to set up
! the Cirs counterterms: 
! q(~) -> q(~) g g or 
! q(~) -> q(~) r r~ or 
! g    -> g q q~ or 
! g    -> g g g
subroutine ReorderRadirs(rad1,rad1ID,rad2,rad2ID,rad3,rad3ID, &
                         radi,radr,rads)
implicit none
!
  integer , intent(in) :: rad1,rad2,rad3
  integer , intent(in) :: rad1ID,rad2ID,rad3ID
  integer , intent(out) :: radi,radr,rads
!
!
!
! In the case of g -> g g g the positions have to be put in
! increasing order:
  if ((abs(rad1ID) + abs(rad2ID) + abs(rad3ID)).eq.0) then
!    print *,"g -> g g g"
    radi = min(rad1,rad2,rad3)
    rads = max(rad1,rad2,rad3)
    if (((radi.eq.rad1).and.(rads.eq.rad2)).or. &
        ((radi.eq.rad2).and.(rads.eq.rad1))) radr = rad3
    if (((radi.eq.rad1).and.(rads.eq.rad3)).or. &
        ((radi.eq.rad3).and.(rads.eq.rad1))) radr = rad2
    if (((radi.eq.rad2).and.(rads.eq.rad3)).or. &
        ((radi.eq.rad3).and.(rads.eq.rad2))) radr = rad1
    return
! g -> g q q~ 
  elseif ((max(rad1ID,rad2ID,rad3ID).eq. &
           -min(rad1ID,rad2ID,rad3ID)).and. &
          (min(abs(rad1ID),abs(rad2ID),abs(rad3ID)).eq.0)) then
!    print *,"g -> g q q~"
    if (rad1ID.eq.0) then
      radi = rad1
! rad2 corresponds to the q quark:
      if (rad2ID.gt.0) then
        radr = rad2
        rads = rad3
      else
        radr = rad3
        rads = rad2
      end if
    elseif (rad2ID.eq.0) then
      radi = rad2
! rad1 corresponds to the q quark:
      if (rad1ID.gt.0) then
        radr = rad1
        rads = rad3
      else
        radr = rad3
        rads = rad1
      end if
    elseif (rad3ID.eq.0) then
      radi = rad3
! rad1 corresponds to the q quark:
      if (rad1ID.gt.0) then
        radr = rad1
        rads = rad2
      else
        radr = rad2
        rads = rad1
      end if
    end if
    return
! q(~) -> q(~) g g or q(~) -> q(~) q q~ or q(~) -> q(~) r r~
  elseif ((rad1ID + rad2ID + rad3ID).ne.0) then
! q(~) -> q(~) g g
    if (min(abs(rad1ID),abs(rad2ID),abs(rad3ID)).eq.0) then
!      print *,"q(~) -> q(~) g g"
! When q(~) -> q(~) g g is considered only the quark has to be
! identified the postions of the gluons are put in increasing
! order:
      if (rad1ID.ne.0) then
        radi = rad1
        radr = min(rad2,rad3)
        rads = max(rad2,rad3)
      elseif (rad2ID.ne.0) then
        radi = rad2
        radr = min(rad1,rad3)
        rads = max(rad1,rad3)
      elseif (rad3ID.ne.0) then
        radi = rad3
        radr = min(rad1,rad2)
        rads = max(rad1,rad2)
      end if
! If a q(~) -> q(~) q q~ is present
    elseif ((abs(rad1ID).eq.abs(rad2ID)).and. &
            (abs(rad1ID).eq.abs(rad3ID)).and. &
            (abs(rad2ID).eq.abs(rad3ID))) then
!      print *,"q(~) -> q(~) q q~"
! The trick is to find those partons which are the same kind,
! whichever out of these two has the lower position it will
! become radi the remaining two are ordered as q q~:
      if (rad1ID.eq.rad2ID) then
        radi = min(rad1,rad2)
! Having q -> q q~
        if (rad1ID.gt.0) then
          radr = max(rad1,rad2)
          rads = rad3
! Having q~ -> q q~
        else
          rads = max(rad1,rad2)
          radr = rad3
        end if
      elseif (rad1ID.eq.rad3ID) then
        radi = min(rad1,rad3)
! Having q -> q q~
        if (rad1ID.gt.0) then
          radr = max(rad1,rad3)
          rads = rad2
! Having q~ -> q q~
        else
          rads = max(rad1,rad3)
          radr = rad2
        end if
      elseif (rad2ID.eq.rad3ID) then
        radi = min(rad2,rad3)
! Having q -> q q~
        if (rad2ID.gt.0) then
          radr = max(rad2,rad3)
          rads = rad1
! Having q~ -> q q~
        else
          rads = max(rad2,rad3)
          radr = rad1
        end if
      end if
      return
! The remaining is the case of q(~) -> q(~) r r~, q != r
    else
!      print *,"q(~) -> q(~) r r~"
! Only have to identify the r-r~ pair:
      if (rad1ID.eq.-rad2ID) then
        radi = rad3
! rad1 corresponds to the r quark:
        if (rad1ID.gt.0) then
        radr = rad1
        rads = rad2
      else
        radr = rad2
        rads = rad1
      end if
      elseif (rad1ID.eq.-rad3ID) then
        radi = rad2
! rad1 corresponds to the r quark:
        if (rad1ID.gt.0) then
          radr = rad1
          rads = rad3
        else
          radr = rad3
          rads = rad1
        end if
      elseif (rad2ID.eq.-rad3ID) then
        radi = rad1
! rad1 corresponds to the r quark:
        if (rad2ID.gt.0) then
          radr = rad2
          rads = rad3
        else
          radr = rad3
          rads = rad2
        end if
      end if
    end if
  end if
!
end subroutine ReorderRadirs
!
! This routine takes two pairs of radiated partons: rad1-rad2
! and rad3-rad4 and orders them such the resulting ordering
! is the same used to obtain the Cirjs terms.
! The legs within the pairs are not needed to be reordered,
! provided they are already in order, that is:
! g -> g g || g -> q q~ || q(~) -> q(~) g
subroutine ReorderRadirjs(emit12,rad1,rad2,emit34,rad3,rad4, &
                          emitir,radi,radr,emitjs,radj,rads)
implicit none
!
  integer , intent(in) :: emit12,emit34
  integer , intent(in) :: rad1,rad2,rad3,rad4
  integer , intent(out) :: emitir,emitjs
  integer , intent(out) :: radi,radr,radj,rads
!
!
! rad1-rad2 comes before rad3-rad4:
  if (min(rad1,rad2).lt.min(rad3,rad4)) then
    emitir = emit12
    emitjs = emit34
    radi = rad1
    radr = rad2
    radj = rad3
    rads = rad4
  else
    emitir = emit34
    emitjs = emit12
    radi = rad3
    radr = rad4
    radj = rad1
    rads = rad2
  end if
!
end subroutine ReorderRadirjs
