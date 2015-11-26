module math
implicit none
!
  real(kind(1d0)) , parameter :: pi = 3.14159265358979323846264338327d0
  real(kind(1d0)) , parameter :: EulerGamma = 0.57721566490153286060651209008240d0
  real(kind(1d0)) , parameter :: polygamma1 = -2.40411380631918857079947632302d0
  real(kind(1d0)) , parameter :: Zeta3 = 1.20205690315959428539973816151d0
  real(kind(1d0)) , parameter :: cGamma = 1d0/(4d0*pi)**2
  real(kind(1d0)) , parameter :: Seps = 1d0
!
contains
subroutine PolyLogs(x,Li234)
  implicit none

!parameters
  real(kind(1d0)) , intent(in) :: x
  real(kind(1d0)) , dimension(2:4), intent(out) :: Li234

!variables
  integer , parameter :: n1 = 0
  integer , parameter :: n2 = 1
  integer , parameter :: weight = 4

  complex(kind(1d0)) , dimension(n1:n2) :: Hc1
  complex(kind(1d0)) , dimension(n1:n2,n1:n2) :: Hc2
  complex(kind(1d0)) , dimension(n1:n2,n1:n2,n1:n2) :: Hc3
  complex(kind(1d0)) , dimension(n1:n2,n1:n2,n1:n2,n1:n2) :: Hc4

  real(kind(1d0)) , dimension(n1:n2) :: Hr1
  real(kind(1d0)) , dimension(n1:n2,n1:n2) :: Hr2
  real(kind(1d0)) , dimension(n1:n2,n1:n2,n1:n2) :: Hr3
  real(kind(1d0)) , dimension(n1:n2,n1:n2,n1:n2,n1:n2) :: Hr4

  real(kind(1d0)) , dimension(n1:n2) :: Hi1
  real(kind(1d0)) , dimension(n1:n2,n1:n2) :: Hi2
  real(kind(1d0)) , dimension(n1:n2,n1:n2,n1:n2) :: Hi3
  real(kind(1d0)) , dimension(n1:n2,n1:n2,n1:n2,n1:n2) :: Hi4

!call HPL routine
  call hplog(x,weight,Hc1,Hc2,Hc3,Hc4,Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2)

!Take the proper HPLs to the Li's
  Li234(2) = Hr2(0,1)
  Li234(3) = Hr3(0,0,1)
  Li234(4) = Hr4(0,0,0,1)
   
!Happy!

end subroutine PolyLogs
end module math
!
module process
implicit none
!
  integer :: nleg_born
! FIXME maybe this should be moved elsewhere...
  integer :: proc_firstlight
  character(72) :: process_name
! FIXME this is not really nice here too...
  real(kind(1d0)) , dimension(:) , allocatable :: fs_masses
!
contains
!
subroutine init_process()
implicit none
!
!
! Read in process name:
  write(*,'(A,1x)',advance='no') "Enter the name of the process: "
  read(*,*) process_name
  print *,"The process name is: ",process_name
!
end subroutine init_process
!
end module process
!
module input
!
contains
!
function nnloinput(stringa)
use process
implicit none
!
  real(kind(1d0)) :: nnloinput
!
  character (len=*) :: stringa
!
  integer , parameter :: maxnum = 500
  character (len = 100) , save :: fname
  character (len = 100) :: line,line0
  character (len = 20) :: string
  character (len = 20) , dimension(maxnum) , save :: keywords
  real(kind(1d0)) , dimension(maxnum) , save :: values
  integer :: ios,j,k,l,imode
  integer , save :: ini = 0
  integer , save :: numvalues
!
  string = stringa
  if (ini.eq.0) then
    fname = process_name(1:len_trim(process_name))//'-nnlo.input'
    open(unit=33,file=fname,status='old',iostat=ios)
    if (ios.ne.0) then
      print *,"Problem occured with input card..."
      print *,fname(1:len_trim(fname))," is not present..."
      stop
    end if
    numvalues = 0
    do l=1,1000000
      line0 = ' '
      read(unit=33,fmt='(a)',iostat=ios) line0
      if (ios.ne.0.and.line0.eq.' ') goto 10
      line = line0
      do k=1,100
        if (line(k:k).eq.'#'.or.line(k:k).eq.'!') then
          line(k:) = ' '
        endif
      enddo
      if (line.ne.' ') then
        if (numvalues.eq.maxnum) then
          write(*,*) ' too many entries in nnloinput.dat'
          stop
        endif
        numvalues = numvalues+1
! skip blanks
 12     if (line(1:1).eq.' ') then
          line=line(2:)
          goto 12
        endif
        k = index(line,' ')
        keywords(numvalues) = line(1:k-1)
        line = line(k+1:)
        read(unit=line,fmt=*,iostat=ios) values(numvalues)
        if (ios.ne.0) then
          write(*,*) ' nnloinput error: cannot parse '
          write(*,'(a)') line0
          stop
        endif
      endif
    enddo
10  continue
    close(33)
    ini=1
  endif
  if (string(1:1).eq.'#') then
    string=string(2:)
    imode=0
  else
    imode=1
  endif
  do j=1,numvalues
    if (string.eq.keywords(j)) then
      nnloinput=values(j)
      return
    endif
  enddo
  if (imode.eq.1) then
    write(*,*) ' nnloinput: keyword ',string,' not found'
    stop
  endif
  nnloinput=-1d6
!
end function nnloinput
!
end module input
!
module collider
implicit none
!
  integer :: beam1,beam2
  real(kind(1d0)) :: ebeam1,ebeam2
  real(kind(1d0)) :: stot,rstot
!
contains
!
subroutine init_collider()
use input
implicit none
!
!
!
! defining beam type:
  beam1 = nnloinput("beam1")
  beam2 = nnloinput("beam2")
! beam energies:
  ebeam1 = nnloinput("ebeam1")
  ebeam2 = nnloinput("ebeam2")
! We calculate CM energies too:
  stot = 2d0*(ebeam1**2 + ebeam2**2)
  rstot = sqrt(stot)
! FIXME PDF...
!
end subroutine init_collider
!
end module collider
!
module flags
implicit none
!
  logical :: flg_LO
  logical :: flg_NLO_V,flg_NLO_R
  logical :: flg_NNLO_VV,flg_NNLO_RV,flg_NNLO_RR
!
  logical :: flg_NLO,flg_NNLO
!
! These flags are meant to control the various subtraction terms:
  logical flg_NLO_R_A1
  logical flg_NNLO_RR_A1
  logical flg_NNLO_RR_A2
  logical flg_NNLO_RR_A12
  logical flg_NNLO_RRA1_A1
  logical flg_NNLO_RV_A1
!
  logical flg_NNLO_R_I1
!
  logical :: flg_cuts,flg_cutfunc
!
  logical :: flg_scalestudy
!
  logical :: flg_analysis
!
  logical :: flg_optim
!
  logical :: flg_scaledep
!
  logical :: flg_manyseeds
!
  logical :: flg_Bijk
!
  logical :: flg_dipolechan,flg_dipsqrchan
!
contains
!
subroutine init_flags()
use input
implicit none
!
!
!
  flg_LO = .true.
  flg_NLO_V = .true.
  flg_NLO_R = .true.
  flg_NNLO_VV = .true.
  flg_NNLO_RV = .true.
  flg_NNLO_RR = .true.
!
  if (nnloinput("#LO").eq.0) flg_LO = .false.
  if (nnloinput("#NLO_V").eq.0) flg_NLO_V = .false.
  if (nnloinput("#NLO_R").eq.0) flg_NLO_R = .false.
  if (nnloinput("#NNLO_VV").eq.0) flg_NNLO_VV = .false.
  if (nnloinput("#NNLO_RV").eq.0) flg_NNLO_RV = .false.
  if (nnloinput("#NNLO_RR").eq.0) flg_NNLO_RR = .false.
  if (nnloinput("#NLO").eq.0) then
    flg_NLO_V = .false.
    flg_NLO_R = .false.
  end if
  if (nnloinput("#NNLO").eq.0) then
    flg_NNLO_VV = .false.
    flg_NNLO_RV = .false.
    flg_NNLO_RR = .false.
  end if
!
  print *,"Calculating contributions: "
  print *,"LO:      ",flg_LO
  print *,"NLO V:   ",flg_NLO_V
  print *,"NLO R:   ",flg_NLO_R
  print *,"NLO:     ",flg_NLO_V.and.flg_NLO_R
  print *,"NNLO VV: ",flg_NNLO_VV
  print *,"NNLO RV: ",flg_NNLO_RV
  print *,"NNLO RR: ",flg_NNLO_RR
  print *,"NNLO:    ",flg_NNLO_VV.and.flg_NNLO_RV.and.flg_NNLO_RR
!
! Flags related to cutting the phase space:
  flg_cuts    = .false.
  flg_cutfunc = .false.
  if (nnloinput("#cuts").gt.0) flg_cuts = .true.
  if (nnloinput("#cutfunc").gt.0) flg_cutfunc = .true.
  if (flg_cuts.and.flg_cutfunc) then
    print *,"Both cuts and cutfunc are active, only one is allowed..."
    stop
  end if
! Flag related to the scale study:
  flg_scalestudy = .false.
  if (nnloinput("#scalestudy").gt.0) flg_scalestudy = .true.
!
! We can turn off histogramming:
  flg_analysis = .true.
  if (nnloinput("#noanalysis").gt.0) flg_analysis = .false.
!
! This flag is used to indicate that only scale dependent
! contributions are needed to be calculated provided by the
! fact that the PS point is the same and the only change
! happened in the unphysical scales:
  flg_scaledep = .false.
!
! By default we always calculate the subtraction terms, but
! for development we can fine-tune the behavior by setting the
! following flags false:
  flg_NLO_R_A1 = flg_NLO_R
  if (nnloinput("#noNLO_R_A1").eq.1) flg_NLO_R_A1 = .false.
  flg_NNLO_RR_A1 = flg_NNLO_RR
  if (nnloinput("#noNNLO_RR_A1").eq.1) flg_NNLO_RR_A1 = .false.
  flg_NNLO_RR_A2 = flg_NNLO_RR
  if (nnloinput("#noNNLO_RR_A2").eq.1) flg_NNLO_RR_A2 = .false.
  flg_NNLO_RR_A12 = flg_NNLO_RR
  if (nnloinput("#noNNLO_RR_A12").eq.1) flg_NNLO_RR_A12 = .false.
! This flag controls the presence of the term I_1^{(0)}\otimes R:
  flg_NNLO_R_I1 = flg_NNLO_RV
  if (nnloinput("#noNNLO_R_I1").eq.1) flg_NNLO_R_I1 = .false.
! the A_1-type subtraction for RR,A1 is associated to the second
! line of the NNLO cross section formula (m+1 parton final state),
! and it regularizes the kinematical singularities of the R_I1 term
! hence only present if R_I1 is also calculated:
  flg_NNLO_RRA1_A1 = flg_NNLO_R_I1
  if (nnloinput("#noNNLO_RRA1_A1").eq.1) flg_NNLO_RRA1_A1 = .false.
  flg_NNLO_RV_A1 = flg_NNLO_RV
  if (nnloinput("#noNNLO_RV_A1").eq.1) flg_NNLO_RV_A1 = .false.
!
  flg_NLO = flg_NLO_V.and.flg_NLO_R.and.flg_NLO_R_A1
!
! FIXME this is incomplete we need the subtractions too:
  flg_NNLO = flg_NNLO_VV.and.flg_NNLO_RV.and.flg_NNLO_RR.and. &
             flg_NNLO_RR_A1.and.flg_NNLO_RR_A2
!
! This flag is implemented to run the code in several instances
! each with a different random seed:
  flg_manyseeds = .false.
  if (nnloinput("#manyseeds").eq.1) flg_manyseeds = .true.
!
! We offer the possibility to optimize also on dipole channels:
! This is turned on by default:
  flg_dipolechan = .true.
  if (nnloinput("#optdipchan").eq.0) flg_dipolechan = .false.
! For the RR contribution not only dipole but squared dipole
! channels are needed to achieve the most out of Kaleu:
  flg_dipsqrchan = .true.
  if (nnloinput("#optdipsqrchan").eq.0) flg_dipsqrchan = .false.
!
end subroutine init_flags
!
end module flags
!
module scales
implicit none
!
  real(kind(1d0)) :: muscale
  real(kind(1d0)) :: mur,muf
  real(kind(1d0)) :: xir,xif
!
  real(kind(1d0)) :: xir_tmp,xif_tmp
!
  logical :: scaleoffdiag
  integer :: nscales
  integer :: centscale
  real(kind(1d0)) :: scalevarpar
  real(kind(1d0)) , dimension(:,:) , allocatable :: xirf,xirf_tmp
!
contains
!
subroutine init_scales()
use collider
use flags
use input
implicit none
!
!
  integer :: iscale,expnR,expnF
  real(kind(1d0)) :: baseR,baseF
!
!
  muscale = nnloinput("muscale")
!
  if (muscale.gt.0) then
    print *,"*********************************************************"
    print *,"*                                                       *"
    print *,"*                  Using static scale...                *"
    print *,"*       mur = ",muscale
    print *,"*                                                       *"
    print *,"*********************************************************"
  else
    print *,"*********************************************************"
    print *,"*                                                       *"
    print *,"*                  Using dynamic scale...               *"
    print *,"*                                                       *"
    print *,"*********************************************************"
  end if
! If xiR or xiF is present we read it in:
  if (nnloinput("#renscalefact").gt.0) then
    xir = nnloinput("#renscalefact")
  else
    xir = 1d0
  end if
  if (nnloinput("#facscalefact").gt.0) then
    xif = nnloinput("#facscalefact")
  else
    xif = 1d0
  end if
! baseF has to be initialized here to avoid an annoying SIGFPE with
! certain versions of gfortran:
  baseF = -1
! We not only initialize the scale but if needed set up everything 
! for scale studies:
  if (flg_scalestudy) then
! If scale variation is demanded we have to know the scale
! variation parameter:
    scalevarpar = nnloinput("scalevarpar")
    write(*,'(a)') "Scale variation will be performed"
! By default there is no off diagonal:
    scaleoffdiag = .false.
    if (nnloinput("scaleoffdiag").gt.0) scaleoffdiag = .true.
! We count down the various, possible scale combinations:
    nscales = 0
    expnR = -1
    expnF = -1
    baseR = scalevarpar
! We only vary the factorization scale if it is necessary.
    if ((abs(beam1).eq.2212).or.(abs(beam2).eq.2212)) then
      baseF = scalevarpar
    else
      baseF = -1
    end if
! We start with the ansatz that we have only one scale:
    allocate(xirf(2,1))
    do while (.true.)
      if ((baseR*expnR.le.scalevarpar).and. &
          (baseF*expnF.le.scalevarpar)) then
        nscales = nscales + 1
! If we have more than one scale we have to reallocate xirf:
        allocate(xirf_tmp(2,nscales-1))
! We dump all the contents of xirf into a temporary array:
        xirf_tmp = xirf
! Than simply deallocate and reallocate:
        deallocate(xirf)
        allocate(xirf(2,nscales))
! Than we move back what we had:
        xirf(1:2,1:nscales-1) = xirf_tmp(1:2,1:nscales-1)
! Deallocating the temp array:
        deallocate(xirf_tmp)
!        print *,"nscales: ",nscales
!        print *,"xiR,xiF: ",baseR**expnR,baseF**expnF
!        print *,"baseR,baseF: ",baseR,baseF
!        print *,"expnR,expnF: ",expnR,expnF
! We store the xi values:
        xirf(1,nscales) = baseR**expnR
        if (baseF.gt.0) then
          xirf(2,nscales) = baseF**expnF
        else
          xirf(2,nscales) = 1d0
        end if
! We always start with incrementing the factorization scale:
        if ((baseF.lt.0).or.(baseF*expnF.eq.scalevarpar).and. &
            scaleoffdiag) then
! If we reached the up-most factorization scale parameter we
! increment the renormalization scale parameter:
          expnF = -1
          baseR = baseR + expnR
! When baseR reaches 1 the exponent should change sign:
          if (baseR.eq.1) expnR = 1
        else
! Increment the factorization scale parameter and flip sign
! for the exponent if needed:
          baseF = baseF + expnF
          if (baseF.eq.1) expnF = 1
! if no off-diagonal terms are considered we simply increase
! the renormalization scale parameter too:
          if (.not.scaleoffdiag) then
            baseR = baseR + expnR
            if (baseR.eq.1) expnR = 1
          end if
        end if
      else
        exit
      end if
    end do
  else
! If no scale variation is performed we allocate only one entry for
! the \xis:
    nscales = 1
    allocate(xirf(2,nscales))
    xirf(1,1) = xir
    xirf(2,1) = xif
  end if
!
! We look for the scale choice which corresponds to the central one:
  do iscale=1,nscales
    if ((xirf(1,iscale).eq.1d0).and. &
        (xirf(2,iscale).eq.1d0)) then
      centscale = iscale
      exit
    end if
  end do
!
  write(*,'(a,1x,I3,1x,a)') "Histograms will be produced for ",nscales,&
                            "scale settings"
  print *,"The following xiR,xiF values will be used:"
  do iscale=1,nscales
    if (baseF.gt.0) then
      write(*,'(2(a,1x,F5.3,1x))') "xiR:",xirf(1,iscale), &
                                   "xiF:",xirf(2,iscale)
    else
      write(*,'(a,1x,F5.3,1x)') "xiR:",xirf(1,iscale)
    end if
  end do
!
end subroutine init_scales
!
subroutine StoreXis_scales
implicit none
!
!
!
! We stash away the original xi values:
  xir_tmp = xir
  xif_tmp = xif
!
end subroutine StoreXis_scales
!
subroutine RestoreXis_scales
use flags
implicit none
!
!
!
! We copy back the old values:
  xir = xir_tmp
  xif = xif_tmp
!
! We also have to turn flg_scaledep back to false:
  flg_scaledep = .false.
!
end subroutine RestoreXis_scales
!
subroutine ChangeXis_scales(iscale)
use flags
implicit none
!
  integer , intent(in) :: iscale
!
!
!
  xir = xirf(1,iscale)
  xif = xirf(2,iscale)
!
! We also change flg_scaledep to true if iscale is greater
! than 1:
  if (iscale.gt.1) then
    flg_scaledep = .true.
  end if
!
end subroutine ChangeXis_scales
!
end module scales
!
module QCDparams
implicit none
!
  integer :: qcd_nf,qcd_nu,qcd_nd
  real(kind(1d0)) :: qcd_nc,qcd_tr,qcd_ca,qcd_cf
  real(kind(1d0)) :: qcd_beta0,qcd_beta1
  real(kind(1d0)) :: qcd_lambda
!
contains
!
subroutine init_QCDparams()
use input
implicit none
!
!
!
  qcd_nu = 2
  qcd_nd = 3
  qcd_nf = qcd_nu + qcd_nd
!
  qcd_nc = 3d0
  qcd_tr = 0.5d0
  qcd_cf = qcd_tr*(qcd_nc**2 - 1d0)/qcd_nc
  qcd_ca = 2d0*qcd_tr*qcd_nc
!
  qcd_beta0 = (11d0*qcd_ca - 4d0*qcd_tr*qcd_nf)/3d0
!
  qcd_beta1 = (17*qcd_ca**2 - 10*qcd_ca*qcd_tr*qcd_nf &
            - 6*qcd_cf*qcd_tr*qcd_nf)/3d0
!
! LambdaQCD can be taken from the input card:
  if (nnloinput("#lambdaQCD").gt.0d0) then
    qcd_lambda = nnloinput("#lambdaQCD")
  end if
!
end subroutine init_QCDparams
!
end module QCDparams
!
module utils
!
contains
!
integer function ConvertToPDG(ch) Result(PDG)
implicit none
!
  character*2 , intent(in) :: ch
!
!
!
  if (ch.eq.'u ') then
    PDG = 2
  elseif (ch.eq.'u~') then
    PDG = -2
  elseif (ch.eq.'d ') then
    PDG =  1
  elseif (ch.eq.'d~') then
    PDG = -1
  elseif (ch.eq.'c ') then
    PDG =  4
  elseif (ch.eq.'c~') then
    PDG = -4
  elseif (ch.eq.'s ') then
    PDG =  3
  elseif (ch.eq.'s~') then
    PDG = -3
  elseif (ch.eq.'b ') then
    PDG =  5
  elseif (ch.eq.'b~') then
    PDG = -5
  elseif (ch.eq.'t ') then
    PDG =  6
  elseif (ch.eq.'t~') then
    PDG = -6
  elseif (ch.eq.'e+') then
    PDG = -11
  elseif (ch.eq.'e-') then
    PDG = 11
  elseif (ch.eq.'A ') then
    PDG = 22
  elseif (ch.eq.'g ') then
    PDG =  0
  else
    print *,"Unknown particle is supplied to ConvertToPDG..."
    print *,"ch = ",ch
    stop
  end if
!
  return
!
end function ConvertToPDG
!
function ConvertFromPDG(PDG) result(ans)
implicit none
!
  integer , intent(in) :: PDG
!
  character*2 :: ans
!
  ans = 'NA'
!
  if (PDG.eq.2) then
    ans = 'u '
  elseif (PDG.eq.-2) then
    ans = 'u~'
  elseif (PDG.eq.1) then
    ans = 'd '
  elseif (PDG.eq.-1) then
    ans = 'd~'
  elseif (PDG.eq.3) then
    ans = 's '
  elseif (PDG.eq.-3) then
    ans = 's~'
  elseif (PDG.eq.4) then
    ans = 'c '
  elseif (PDG.eq.-4) then
    ans = 'c~'
  elseif (PDG.eq.5) then
    ans = 'b '
  elseif (PDG.eq.-5) then
    ans = 'b~'
  elseif (PDG.eq.6) then
    ans = 't '
  elseif (PDG.eq.-6) then
    ans = 't~'
  elseif (PDG.eq.11) then
    ans = 'e-'
  elseif (PDG.eq.-11) then
    ans = 'e+'
  elseif (PDG.eq.22) then
    ans = 'A '
  elseif (PDG.eq.0) then
    ans = 'g '
  else
    print *,"Unknown PDG number is supplied..."
    print *,"PDG = ",PDG
    stop "ConvertFromPDG"
  end if
!
end function ConvertFromPDG
!
! Note that in the first version of this routine the 
! check upon final state started ot the position of the 
! first light parton position it is modified to start the
! check at the very first final state parton instead.
function CheckEquivProc(nleg,numproc,lst_proc,flv) result(answer)
use process
implicit none
!
  logical :: answer
  integer , intent(in) :: nleg,numproc
  integer , dimension(:,:) , intent(in) :: lst_proc
  integer , dimension(nleg) , intent(in) :: flv
!
  integer :: ileg,jleg,iproc
  integer , dimension(nleg) :: flv_tmp
!
!
  answer = .true.
!
!
!  print *,"************************************************"
!  write(*,'(a,20(I0,1x))') "flv:    ",flv
! We loop over all subprocesses we have so far:
  iproc = 0
  do while (.true.)
100 continue
    iproc = iproc + 1
    if (iproc.gt.numproc) exit
! We copy the light part into flv_tmp:
    flv_tmp = lst_proc(:,iproc)
!    write(*,'(a,20(I2,1x))') "flv_tmp:",flv_tmp
! The initial state should be the same:
    if ((flv_tmp(1).ne.flv(1)).or.(flv_tmp(2).ne.flv(2))) goto 100
    do ileg=3,nleg
      do jleg=3,nleg
        if (flv_tmp(jleg).eq.-10000) cycle
        if (flv_tmp(jleg).eq.flv(ileg)) then
          flv_tmp(jleg) = -10000
          goto 110
        end if
      end do
110   continue
    end do
! Checking for equivalence:
!    write(*,'(a,20(I6,1x))') "flv_tmp after: ",flv_tmp
!    read(*,*)
    do ileg=3,nleg
      if (flv_tmp(ileg).ne.-10000) goto 100
    end do
! If we are at this point the configuration is already considered.
    return
  end do
!
  answer = .false.
!
  return
!
end function CheckEquivProc
!
! This routine checks the flavour configuration for physicalness:
function CheckConfig(nleg,flv) result(answer)
use process
use QCDparams
implicit none
!
  logical :: answer
  integer , intent(in) :: nleg
  integer , dimension(nleg) , intent(in) :: flv
!
  integer :: ileg,iflv
  integer , parameter :: numflv = 6
  integer , dimension(numflv) :: quark_arr
!
!
  answer = .false.
!
! We go through all the legs which are connected to massless
! partons and and increase or decrease the content of quark_arr.
! We work in an all-outgoing scheme:
  quark_arr = 0
  do ileg=1,nleg
    if (flv(ileg).eq.0) cycle
    if ((ileg.gt.2).and.(ileg.lt.proc_firstlight)) cycle
    if (abs(flv(ileg)).gt.qcd_nf) cycle
    if (ileg.lt.3) then
      quark_arr(abs(flv(ileg))) = quark_arr(abs(flv(ileg))) &
                                - sign(1,flv(ileg))
    else
      quark_arr(abs(flv(ileg))) = quark_arr(abs(flv(ileg))) &
                                + sign(1,flv(ileg))
    end if
  end do
!
! If the quark array is not zero we have a bad configuration:
  do iflv=1,numflv
    if (quark_arr(iflv).ne.0) return
  end do
!
!  write(*,'(a,20(I3,1x))') "flv: ",flv
!  write(*,'(a,20(I3,1x))') "quark_arr: ",quark_arr
!
  answer = .true.
!
end function CheckConfig
!
! Calculates the symmetry factor for a subprocess:
function CalcSymFact(flv) result(Sm)
implicit none
!
  integer :: Sm
!
  integer , dimension(:) , intent(in) :: flv
!
  integer :: ileg,iflv,nleg
  integer , dimension(-6:6) :: flv_arr
!
  integer , external :: ifactorial
!
! We take the number of legs from the flavor array itself:
  nleg = size(flv)
! We nullify the number of each flavor:
  flv_arr = 0
! Loop over final state particles skip non-QCD ones:
  do ileg=3,nleg
    if (abs(flv(ileg)).gt.6) cycle
! Increment:
    flv_arr(flv(ileg)) = flv_arr(flv(ileg)) + 1
  end do
  Sm = 1
! The symmetry factor is the product of the factorials taken
! for each item in flv_arr:
  do iflv=-6,6
    Sm = Sm * ifactorial(flv_arr(iflv))
  end do
!
end function CalcSymFact
!
end module utils
!
module random_RN
implicit none
!
  integer , private :: seed
!
contains
!
subroutine put_seed(iseed)
implicit none
!
  integer , intent(in) :: iseed
!
!
!
  seed = iseed
!
end subroutine put_seed
!
function get_random() result(randnum)
implicit none
!
  real(kind(1d0)) :: randnum
!
  integer :: hi,lo
!
  integer , parameter :: m = 2147483647
  integer , parameter :: a = 16807
  integer , parameter :: q = 127773
  integer , parameter :: r = 2836
  real(kind(1d0)) , parameter :: minv = 0.46566128752458d-09
!
  hi = seed/q
  lo = mod(seed,q)
  seed = a*lo - r*hi
  if (seed.le.0) seed = seed + m
  randnum = seed*minv
!
end function get_random
!
end module random_RN
!
module random_RANLUX
implicit none
!
!
contains
!
subroutine put_seed(iseed)
implicit none
!
  integer , intent(in) :: iseed
!
!
!
  call RLUXGO(4,iseed,0,0)
!
end subroutine put_seed
!
function get_random() result(randnum)
implicit none
!
  real(kind(1d0)) :: randnum
!
  real(kind(1d0)) , dimension(1) :: ran
!
!
!
  call RANLUX(ran,1) 
  randnum = ran(1)
!
end function get_random
!
end module random_RANLUX
!
module random
use random_RN
!use random_RANLUX
use flags
implicit none
!
  integer :: seed
  integer :: iseed_manyseeds
!
contains
!
subroutine init_seed()
use process
use flags
use input
implicit none
!
!
  integer , parameter :: iun = 99
  character (len=256) :: fname
  integer :: istat,iline,iseed,seed_tmp
!
  seed = 123456
!
  if (nnloinput("#seed").gt.0) then
    seed = nnloinput("#seed")
  end if
!
! If the manyseeds option is activated the seed is taken from a file:
  if (flg_manyseeds) then
    fname = trim(process_name)//"-seeds.dat"
    open(unit=iun,file=fname,status='old',iostat=istat)
    if (istat.ne.0) then
      print *,"Problem occured during opening ",fname
      stop
    end if
    write(*,fmt='(a)',advance='no') "Which seed do you want to use: "
    read(*,*) iseed
    iseed_manyseeds = iseed
    do iline=1,iseed
      read(iun,*) seed_tmp
      if (iline.eq.iseed) seed = seed_tmp
    end do
    close(iun)
  end if
!
  print *,"rrrrrrrrrrrrrrrrrrrr*********************rrrrrrrrrrrrrrrrrr"
  print *,"The random number generator is initialized with: ",seed
  print *,"rrrrrrrrrrrrrrrrrrrr*********************rrrrrrrrrrrrrrrrrr"
!
  call put_seed(seed)
!
end subroutine init_seed
!
! This is a crude f90 compliant realization of the RN algorithm:
function gen_rand() result(randnum)
implicit none
!
  real(kind(1d0)) :: randnum
!
!
!
!
  randnum = get_random()
!
end function gen_rand
!
end module random
!
module momenta
implicit none
!
! We define a new type called mom to hold momenta in a consise way:
  type mom
    real(kind(1d0)) :: E,px,py,pz
  end type mom
! We also define a new type called cmom which is the complex equivalent
! of the previous type definition:
  type cmom
    complex(kind(1d0)) :: E,px,py,pz
  end type cmom
!
  interface operator (+)
    module procedure add_momenta
  end interface operator (+)
  interface operator (*)
    module procedure mult_mom_int
    module procedure mult_mom_arr_int
    module procedure mult_mom_real
    module procedure mult_mom_arr_real
    module procedure mult_cmom
    module procedure mult_momenta
    module procedure mult_cmom_cmom
    module procedure mult_mom_cmom
    module procedure mult_cmom_mom
    module procedure mult_Amunu_mom
    module procedure mult_mom_Amunu
  end interface operator (*)
  interface operator (**)
    module procedure exp_mom
  end interface operator (**)
  interface operator (/)
    module procedure div_mom
    module procedure div_cmom_c
    module procedure div_cmom_r
  end interface operator (/)
  interface operator (-)
    module procedure subtr_momenta
    module procedure neg_mom
    module procedure neg_cmom
  end interface operator (-)
  interface assignment (=)
    module procedure init_zero_mom
    module procedure init_zero_cmom
    module procedure convertto_arr14
    module procedure convertfrom_arr14
  end interface assignment (=)
!
  interface PrintMom
    module procedure PrintMom_sngle
    module procedure PrintMom_arr
  end interface
!
contains
!
function add_momenta(p1,p2) result(p12)
implicit none
!
  type(mom) :: p12
!
  type(mom) , intent(in) :: p1,p2
!
  p12%E  = p1%E + p2%E
  p12%px = p1%px + p2%px
  p12%py = p1%py + p2%py
  p12%pz = p1%pz + p2%pz
!
end function add_momenta
!
function subtr_momenta(p1,p2) result(p12)
implicit none
!
  type(mom) :: p12
!
  type(mom) , intent(in) :: p1,p2
!
  p12%E  = p1%E - p2%E
  p12%px = p1%px - p2%px
  p12%py = p1%py - p2%py
  p12%pz = p1%pz - p2%pz
!
end function subtr_momenta
!
function neg_mom(p) result(q)
implicit none
!
  type(mom) :: q
!
  type(mom) , intent(in) :: p
!
  q%E  = - p%E
  q%px = - p%px
  q%py = - p%py
  q%pz = - p%pz
!
end function neg_mom
!
function neg_cmom(p) result(q)
implicit none
!
  type(cmom) :: q
!
  type(cmom) , intent(in) :: p
!
  q%E  = - p%E
  q%px = - p%px
  q%py = - p%py
  q%pz = - p%pz
!
end function neg_cmom
!
subroutine init_zero_mom(p,x)
implicit none
!
  type(mom) , intent(out) :: p
  real(kind(1d0)) , intent(in) :: x
!
  p%E  = 0d0
  p%px = 0d0
  p%py = 0d0
  p%pz = 0d0
!
end subroutine init_zero_mom
!
subroutine init_zero_cmom(p,x)
implicit none
!
  type(cmom) , intent(out) :: p
  real(kind(1d0)) , intent(in) :: x
!
  p%E  = cmplx(0d0,0d0)
  p%px = cmplx(0d0,0d0)
  p%py = cmplx(0d0,0d0)
  p%pz = cmplx(0d0,0d0)
!
end subroutine init_zero_cmom
!
function mult_mom_int(i,p) result(q)
implicit none
!
  type(mom) :: q
!
  integer , intent(in) :: i
  type(mom) , intent(in) :: p

!
  q%E  = i*p%E
  q%px = i*p%px
  q%py = i*p%py
  q%pz = i*p%pz
!
end function mult_mom_int
!
function mult_mom_real(x,p) result(q)
implicit none
!
  type(mom) :: q
!
  real(kind(1d0)) , intent(in) :: x
  type(mom) , intent(in) :: p

!
  q%E  = x*p%E
  q%px = x*p%px
  q%py = x*p%py
  q%pz = x*p%pz
!
end function mult_mom_real
!
! This function realizes (p.p)^(n/2)
function exp_mom(p,n) result(res)
implicit none
!
  real(kind(1d0)) :: res
!
  type(mom) , intent(in) :: p
  integer , intent(in) :: n

!
! Sanity check:
  if (mod(n,2).ne.0) then
    print *,"Error in exp_mom"
    print *,"momentum can only be raised to an even power..."
    stop "exp_mom"
    stop
  end if
  res = p*p
  res = res**(n/2)
!
end function exp_mom
!
function mult_mom_arr_int(i,p_arr) result(q_arr)
implicit none
!
  integer , intent(in) :: i
  type(mom) , dimension(:) , intent(in) :: p_arr
  type(mom) , dimension(size(p_arr)) :: q_arr
!
  integer :: ipart,npart
!
  npart = size(p_arr)
  do ipart=1,npart
    q_arr(ipart) = i*p_arr(ipart)
  end do
!
end function mult_mom_arr_int
!
function mult_mom_arr_real(x,p_arr) result(q_arr)
implicit none
!
  real(kind(1d0)) , intent(in) :: x
  type(mom) , dimension(:) , intent(in) :: p_arr
  type(mom) , dimension(size(p_arr)) :: q_arr
!
  integer :: ipart,npart
!
  npart = size(p_arr)
  do ipart=1,npart
    q_arr(ipart) = x*p_arr(ipart)
  end do
!
end function mult_mom_arr_real
!
function mult_cmom(x,p) result(q)
implicit none
!
  type(cmom) :: q
!
  real(kind(1d0)) , intent(in) :: x
  type(cmom) , intent(in) :: p

!
  q%E  = x*p%E
  q%px = x*p%px
  q%py = x*p%py
  q%pz = x*p%pz
!
end function mult_cmom
!
function mult_momenta(p1,p2) result(p1p2)
implicit none
!
  real(kind(1d0)) :: p1p2
!
  type(mom) , intent(in) :: p1,p2
!
  p1p2 = p1%E*p2%E - p1%px*p2%px - p1%py*p2%py - p1%pz*p2%pz
!
end function mult_momenta
!
function mult_cmom_cmom(p1,p2) result(p1p2)
implicit none
!
  complex(kind(1d0)) :: p1p2
!
  type(cmom) , intent(in) :: p1,p2
!
  p1p2 = p1%E*p2%E - p1%px*p2%px - p1%py*p2%py - p1%pz*p2%pz
!
end function mult_cmom_cmom
!
function mult_mom_cmom(p1,p2) result(p1p2)
implicit none
!
  complex(kind(1d0)) :: p1p2
!
  type(mom) , intent(in) :: p1
  type(cmom) , intent(in) :: p2
!
  p1p2 = p1%E*p2%E - p1%px*p2%px - p1%py*p2%py - p1%pz*p2%pz
!
end function mult_mom_cmom
!
function mult_cmom_mom(p1,p2) result(p1p2)
implicit none
!
  complex(kind(1d0)) :: p1p2
!
  type(cmom) , intent(in) :: p1
  type(mom) , intent(in) :: p2
!
  p1p2 = p1%E*p2%E - p1%px*p2%px - p1%py*p2%py - p1%pz*p2%pz
!
end function mult_cmom_mom
!
function div_mom(p,x) result(q)
implicit none
!
  type(mom) :: q
!
  type(mom) , intent(in) :: p
  real(kind(1d0)) , intent(in) :: x

!
  q%E  = p%E/x
  q%px = p%px/x
  q%py = p%py/x
  q%pz = p%pz/x
!
end function div_mom
!
function div_cmom_c(p,x) result(q)
implicit none
!
  type(cmom) :: q
!
  type(cmom) , intent(in) :: p
  complex(kind(1d0)) , intent(in) :: x

!
  q%E  = p%E/x
  q%px = p%px/x
  q%py = p%py/x
  q%pz = p%pz/x
!
end function div_cmom_c
!
function div_cmom_r(p,x) result(q)
implicit none
!
  type(cmom) :: q
!
  type(cmom) , intent(in) :: p
  real(kind(1d0)) , intent(in) :: x

!
  q%E  = p%E/x
  q%px = p%px/x
  q%py = p%py/x
  q%pz = p%pz/x
!
end function div_cmom_r
!
subroutine convertto_arr14(pout,pin)
implicit none
!
  type(mom) , intent(in) :: pin
  real(kind(1d0)) , dimension(4) , intent(out) :: pout
!
  pout(4) = pin%E
  pout(1) = pin%px
  pout(2) = pin%py
  pout(3) = pin%pz
!
end subroutine convertto_arr14
!
subroutine convertfrom_arr14(pout,pin)
implicit none
!
  real(kind(1d0)) , dimension(4) , intent(in) :: pin
  type(mom) , intent(out) :: pout
!
  pout%E  = pin(4)
  pout%px = pin(1)
  pout%py = pin(2)
  pout%pz = pin(3)
!
end subroutine convertfrom_arr14
!
! The input is considered as A^{\mu\nu} and p^\mu
function mult_Amunu_mom(A,p) result(q)
implicit none
!
  type(mom) :: q
!
  real(kind(1d0)) , dimension(0:3,0:3) , intent(in) :: A
  type(mom) , intent(in) :: p
!
!
!
  q%E  = A(0,0)*p%E - A(0,1)*p%px - A(0,2)*p%py - A(0,3)*p%pz
  q%px = A(1,0)*p%E - A(1,1)*p%px - A(1,2)*p%py - A(1,3)*p%pz
  q%py = A(2,0)*p%E - A(2,1)*p%px - A(2,2)*p%py - A(2,3)*p%pz
  q%pz = A(3,0)*p%E - A(3,1)*p%px - A(3,2)*p%py - A(3,3)*p%pz
!
end function mult_Amunu_mom
!
! The input is considered as A^{\mu\nu} and p^\mu
function mult_mom_Amunu(p,A) result(q)
implicit none
!
  type(mom) :: q
!
  type(mom) , intent(in) :: p
  real(kind(1d0)) , dimension(0:3,0:3) , intent(in) :: A
!
!
!
  q%E  = A(0,0)*p%E - A(1,0)*p%px - A(2,0)*p%py - A(3,0)*p%pz
  q%px = A(0,1)*p%E - A(1,1)*p%px - A(2,1)*p%py - A(3,1)*p%pz
  q%py = A(0,2)*p%E - A(1,2)*p%px - A(2,2)*p%py - A(3,2)*p%pz
  q%pz = A(0,3)*p%E - A(1,3)*p%px - A(2,3)*p%py - A(3,3)*p%pz
!
end function mult_mom_Amunu
!
! FIXME this routine can be made more fancy with pt, rapidity,azimuth...
subroutine PrintMom_sngle(p)
implicit none
!
  type(mom) , intent(in) :: p
!
  write(*,'(A,D22.12,1x,A,D22.12,1x,A,D22.12,1x,A,D22.12)') &
    "E= ",p%E,"px= ",p%px,"py= ",p%py,"pz= ",p%pz
!
end subroutine PrintMom_sngle
!
subroutine PrintMom_arr(p)
implicit none
!
  type(mom) , dimension(:) , intent(in) :: p
!
  integer :: ipart
!
  write(*,'(A,1x,I2,1x,A)') "The array contains",size(p),"momenta: "
  do ipart=1,size(p)
    call PrintMom_sngle(p(ipart))
  end do
!
end subroutine PrintMom_arr
!
function createmom(qE,qx,qy,qz) result(p)
implicit none
!
  type(mom) :: p
!
  real(kind(1d0)) , intent(in) :: qE,qx,qy,qz
!
  p%E  = qE
  p%px = qx
  p%py = qy
  p%pz = qz
!
end function createmom
!
subroutine CalcSij(p)
implicit none
!
  type(mom) , dimension(:) :: p
!
  integer :: npart
  integer :: ipart,jpart
!
!
  npart = size(p)
  write(*,'(a,I0,a)') "The array holds ",npart," momenta"
  do ipart=1,npart
    do jpart=1,npart
      write(*,fmt='(E14.6)',advance='no') (p(ipart) + p(jpart)) * &
                                          (p(ipart) + p(jpart))
    end do
    write(*,*)
  end do
!
end subroutine CalcSij
!
subroutine gen_ranmom(nleg,parr)
use process
use collider
use random
implicit none
!
  integer , intent(in) :: nleg
  type(mom) , dimension(:) , intent(out) :: parr
!
  integer :: ileg
  real(kind(1d0)) :: avgEn
  real(kind(1d0)) :: xi1,xi2,xi3
  real(kind(1d0)) :: px,py,pz,pabs
  real(kind(1d0)) :: beta
  real(kind(1d0)) , dimension(3) :: vec
  type(mom) :: p_tmp
!
! To generate a set of random momenta we need a scale,
! we use the hadronic CM energy
  avgEn = rstot / (nleg - 2)
!
  p_tmp = 0d0
  do ileg=3,nleg-1
    xi1 = gen_rand()
    xi2 = gen_rand()
    xi3 = gen_rand()
    px  = avgEn*(1d0 - 2d0*xi1)
    py  = avgEn*(1d0 - 2d0*xi2)
    pz  = avgEn*(1d0 - 2d0*xi3)
    pabs = sqrt(px**2 + py**2 + pz**2)
    parr(ileg) = createmom(sqrt(fs_masses(ileg)**2 + pabs**2),px,py,pz)
    p_tmp = p_tmp + parr(ileg)
  end do
! The transverse components of the last momentum are given 
! by conservation:
  parr(nleg)%px = -p_tmp%px
  parr(nleg)%py = -p_tmp%py
  parr(nleg)%pz = avgEn*(1d0 - 2d0*gen_rand())
  pabs = sqrt(parr(nleg)%px**2 + parr(nleg)%py**2 + parr(nleg)%pz**2)
  parr(nleg)%E = sqrt(fs_masses(nleg)**2 + pabs**2)
  p_tmp = 0d0
! In order to boost to the CM frame we calculate the total momentum:
  do ileg=3,nleg
    p_tmp = p_tmp + parr(ileg)
  end do
  beta = p_tmp%pz/p_tmp%E
  vec = 0d0
  vec(3) = -1d0
  do ileg=3,nleg
    call boost4vec(beta,vec,parr(ileg),parr(ileg))
  end do
  p_tmp = 0d0
  do ileg=3,nleg
    p_tmp = p_tmp + parr(ileg)
  end do
  parr(1) = createmom(p_tmp%E/2d0,0d0,0d0,p_tmp%E/2d0)
  parr(2) = createmom(p_tmp%E/2d0,0d0,0d0,-p_tmp%E/2d0)
!
end subroutine gen_ranmom
!
function CheckInvariants(p,smin) result(ans)
implicit none
!
  logical :: ans
!
  type(mom) , dimension(:) , intent(in) :: p
  real(kind(1d0)) , intent(in) :: smin
!
  integer :: ipart,jpart,npart
  real(kind(1d0)) :: sij
!
  ans = .false.
!
  npart = size(p)
  do ipart=1,npart-1
    do jpart=ipart+1,npart
      sij = (p(ipart) + p(jpart))*(p(ipart) + p(jpart))
      if (sij.lt.smin) return
    end do
  end do
!
  ans = .true.
!
end function CheckInvariants
!
end module momenta
!
module particles
use momenta
implicit none
!
type particle 
  type(mom) :: p
  integer   :: flv
end type particle
!
interface operator (*)
  module procedure mult_part
end interface operator (*)
  interface operator (/)
    module procedure div_part
end interface operator (/)
interface cross
  module procedure cross_pi_pj_mom_mom
  module procedure cross_pi_pj_part_part
end interface cross
interface CreateParts
  module procedure CreateParts_mom_arr
  module procedure CreateParts_mom_sngle
  module procedure CreateParts_sngle
  module procedure CreateParts_arr
end interface CreateParts
interface PrintParts
  module procedure PrintParts_arr
  module procedure PrintParts_sngle
end interface PrintParts
!
interface PrintSubproc
  module procedure PrintSubproc_flv
  module procedure PrintSubproc_file_flv
  module procedure PrintSubproc_parts
end interface PrintSubproc
!
interface gen_ranparts
  module procedure gen_ranparts_mom
  module procedure gen_ranparts_wflv
end interface gen_ranparts
!
interface PrintInvariants
  module procedure PrintInvariants_parts
end interface PrintInvariants
!
interface IsNaNmom
  module procedure IsNaNmom_mom
  module procedure IsNaNmom_part
end interface IsNaNmom
!
contains
!
subroutine CreateParts_mom_arr(p_arr,flv_arr,parts)
implicit none
!
  type(mom) , dimension(:) , intent(in) :: p_arr
  integer , dimension(:) , intent(in) :: flv_arr
  type(particle) , dimension(:) , intent(out) :: parts
!
  integer :: ipart
!
! We run through all the momenta:
  do ipart=1,size(p_arr)
    parts(ipart)%p = p_arr(ipart)
    parts(ipart)%flv = flv_arr(ipart)
  end do
!
end subroutine CreateParts_mom_arr
!
subroutine CreateParts_mom_sngle(p,flv,part)
implicit none
!
  type(mom) , intent(in) :: p
  integer , intent(in) :: flv
  type(particle) , intent(out) :: part
!
  part%p   = p
  part%flv = flv
!
end subroutine CreateParts_mom_sngle
!
subroutine CreateParts_sngle(pin,flv,part)
use momenta
implicit none
!
  real(kind(1d0)) , dimension(4) , intent(in) :: pin
  integer , intent(in) :: flv
  type(particle) , intent(out) :: part
!
  part%p   = pin
  part%flv = flv
!
end subroutine CreateParts_sngle
!
subroutine CreateParts_arr(pin,flv,part)
use momenta
implicit none
!
  real(kind(1d0)) , dimension(:,:) , intent(in) :: pin
  integer , dimension(:) , intent(in) :: flv
  type(particle) , dimension(:) , intent(out) :: part
!
  integer :: ipart
  type(mom) :: p
!
  do ipart=1,size(part)
    p = pin(1:4,ipart)
    part(ipart)%p   = p
    part(ipart)%flv = flv(ipart)
  end do
!
end subroutine CreateParts_arr
!
subroutine PrintParts_arr(p)
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
!
  integer :: ipart
!
  print 100,"The particle array holds",size(p),"particles"
  do ipart=1,size(p)
    call PrintParts_sngle(p(ipart))
  end do
!
  100 format(a,1x,I3,1x,a)
!
end subroutine PrintParts_arr
!
subroutine PrintParts_sngle(p)
use utils
implicit none
!
  type(particle) , intent(in) :: p
!
!
  write(*,fmt=100,advance='no') ConvertFromPDG(p%flv)
  call PrintMom(p%p)
  100 format ("PDG: ",a2,1x,',',1x)
!
end subroutine PrintParts_sngle
!
subroutine PrintSubproc_flv(flv)
use utils
implicit none
!
  integer , dimension(:) , intent(in) :: flv
!
  integer :: ipart
!
  do ipart=1,size(flv)
    write(*,'(A,1x)',advance='no') ConvertFromPDG(flv(ipart))
    if (ipart.eq.2) write(*,'(A)',advance='no') "-> " 
    if (ipart.eq.size(flv)) write(*,*)
  end do
!
end subroutine PrintSubproc_flv
!
subroutine PrintSubproc_file_flv(iun,flv)
use utils
implicit none
!
  integer , intent(in) :: iun
  integer , dimension(:) , intent(in) :: flv
!
  integer :: ipart
!
  do ipart=1,size(flv)
    write(iun,'(A,1x)',advance='no') ConvertFromPDG(flv(ipart))
    if (ipart.eq.2) write(iun,'(A)',advance='no') "-> "
    if (ipart.eq.size(flv)) write(iun,*)
  end do
!
end subroutine PrintSubproc_file_flv
!
subroutine PrintSubproc_parts(p)
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
!
  call PrintSubproc_flv(p(:)%flv)
!
end subroutine PrintSubproc_parts
!
function mult_part(x,p) result(q)
implicit none
!
  type(particle) :: q
!
  real(kind(1d0)) , intent(in) :: x
  type(particle) , intent(in) :: p
!
!
  q%flv = p%flv
  q%p   = x*p%p
  return
!
end function mult_part
!
function div_part(p,x) result(q)
implicit none
!
  type(particle) :: q
!
  type(particle) , intent(in) :: p
  real(kind(1d0)) , intent(in) :: x
!
!
  q%flv = p%flv
  q%p   = p%p/x
  return
!
end function div_part
!
! This is just a carbon-copy of gen_ranmom but it takes
! particles instead of momenta and at the very end fills
! up the flavor array with zeros:
subroutine gen_ranparts_mom(nleg,parr)
use process
use collider
use random
implicit none
!
  integer , intent(in) :: nleg
  type(particle) , dimension(:) , intent(out) :: parr
!
  integer :: ileg
  real(kind(1d0)) :: avgEn
  real(kind(1d0)) :: xi1,xi2,xi3
  real(kind(1d0)) :: px,py,pz,pabs
  real(kind(1d0)) :: beta
  real(kind(1d0)) , dimension(3) :: vec
  type(mom) :: p_tmp
!
! To generate a set of random momenta we need a scale,
! we use the hadronic CM energy
  avgEn = rstot / (nleg - 2)
!
  p_tmp = 0d0
  do ileg=3,nleg-1
    xi1 = gen_rand()
    xi2 = gen_rand()
    xi3 = gen_rand()
    px  = avgEn*(1d0 - 2d0*xi1)
    py  = avgEn*(1d0 - 2d0*xi2)
    pz  = avgEn*(1d0 - 2d0*xi3)
    pabs = sqrt(px**2 + py**2 + pz**2)
    parr(ileg)%p = createmom(sqrt(fs_masses(ileg)**2 + pabs**2),px,py,pz)
    p_tmp = p_tmp + parr(ileg)%p
  end do
! The transverse components of the last momentum are given 
! by conservation:
  parr(nleg)%p%px = -p_tmp%px
  parr(nleg)%p%py = -p_tmp%py
  parr(nleg)%p%pz = avgEn*(1d0 - 2d0*gen_rand())
  pabs = sqrt(parr(nleg)%p%px**2 + parr(nleg)%p%py**2 &
       +      parr(nleg)%p%pz**2)
  parr(nleg)%p%E = sqrt(fs_masses(nleg)**2 + pabs**2)
  p_tmp = 0d0
! In order to boost to the CM frame we calculate the total momentum:
  do ileg=3,nleg
    p_tmp = p_tmp + parr(ileg)%p
  end do
  beta = p_tmp%pz/p_tmp%E
  vec = 0d0
  vec(3) = -1d0
  do ileg=3,nleg
    call boost4vec(beta,vec,parr(ileg)%p,parr(ileg)%p)
  end do
  p_tmp = 0d0
  do ileg=3,nleg
    p_tmp = p_tmp + parr(ileg)%p
  end do
  parr(1)%p = createmom(p_tmp%E/2d0,0d0,0d0,p_tmp%E/2d0)
  parr(2)%p = createmom(p_tmp%E/2d0,0d0,0d0,-p_tmp%E/2d0)
!
! We simply fill up the flavor with zeros, that is with gluons:
  do ileg=1,nleg
    parr(ileg)%flv = 0
  end do
!
end subroutine gen_ranparts_mom
!
! This routine does the same as the previous one, but it
! takes as input a flavor list and assign it to p:
subroutine gen_ranparts_wflv(nleg,flv,p)
implicit none
!
  integer , intent(in) :: nleg
  integer , dimension(:) , intent(in) :: flv
  type(particle) , dimension(:) , intent(out) :: p
!
  integer :: ipart
!
!
  if (size(flv).ne.size(p)) then
    print *,"Size mismatch in gen_ranparts_wflv..."
    print *,"size(flv): ",size(flv)
    print *,"size(p):   ",size(p)
    stop
  end if
!
  call gen_ranparts_mom(nleg,p)
!
  do ipart=1,size(flv)
    p(ipart)%flv = flv(ipart)
  end do
!
end subroutine gen_ranparts_wflv
!
subroutine PrintInvariants_parts(p)
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
!
  integer :: ipart,jpart,npart
!
!
  npart = size(p)
  call PrintParts(p)
  print *,"The invariants: "
  do ipart=1,npart
    do jpart=ipart+1,npart
      write(*,'(a,2(I2,1x))') "sij, i,j: ",ipart,jpart
      write(*,*) (p(ipart)%p + p(jpart)%p)*(p(ipart)%p + p(jpart)%p)
    end do
  end do
!
end subroutine PrintInvariants_parts
!
function cross_pi_pj_mom_mom(pi,pj) result(pixpj)
implicit none
!
  real(kind(1d0)) , dimension(3) :: pixpj
!
  type(mom) , intent(in) :: pi,pj
!
!
  pixpj(1) = pi%py*pj%pz - pi%pz*pj%py
  pixpj(2) = pi%pz*pj%px - pi%px*pj%pz
  pixpj(3) = pi%px*pj%py - pi%py*pj%px
!
end function cross_pi_pj_mom_mom
!
function cross_pi_pj_part_part(pi,pj) result(pixpj)
implicit none
!
  real(kind(1d0)) , dimension(3) :: pixpj
!
  type(particle) , intent(in) :: pi,pj
!
!
  pixpj = cross_pi_pj_mom_mom(pi%p,pj%p)
!
end function cross_pi_pj_part_part
!
! This routine checks an array of momenta for NaNs, if inan is
! zero at least one NaN is found:
subroutine IsNaNmom_mom(p,inan)
use momenta
implicit none
!
  type(mom) , dimension(:) , intent(in) :: p
  integer , intent(out) :: inan
!
  integer ipart,npart
!
!
  inan = 0
!
  npart = size(p)
!
  do ipart=1,npart
    if ((p(ipart)%px.ne.p(ipart)%px).or. &
        (p(ipart)%py.ne.p(ipart)%py).or. &
        (p(ipart)%pz.ne.p(ipart)%pz).or. &
        (p(ipart)%E.ne.p(ipart)%E)) return
  end do
!
  inan = 1
!
end subroutine IsNaNmom_mom
!
subroutine IsNaNmom_part(p,inan)
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
  integer , intent(out) :: inan
!
!
!
!
  call IsNaNmom_mom(p(:)%p,inan)
!
end subroutine IsNaNmom_part
!
end module particles
!
module subprocesses
implicit none
!
  type eqv_cont
    integer :: eqv
    real(kind(1d0)) :: coeff
  end type eqv_cont
!
  integer :: num_flv_LO,num_flv_NLO_V,num_flv_NLO_R, &
             num_flv_NNLO_VV,num_flv_NNLO_RV,num_flv_NNLO_RR
!
  integer , allocatable , dimension(:,:) :: flv_LO
  integer , allocatable , dimension(:,:) :: flv_NLO_V
  integer , allocatable , dimension(:,:) :: flv_NLO_R
  integer , allocatable , dimension(:,:) :: flv_NNLO_VV
  integer , allocatable , dimension(:,:) :: flv_NNLO_RV
  integer , allocatable , dimension(:,:) :: flv_NNLO_RR
!
  integer :: num_flv_ch_Bkin,num_flv_ch_Rkin,num_flv_ch_RRkin
!
  integer , allocatable , dimension(:,:) :: flv_ch_Bkin
  integer , allocatable , dimension(:,:) :: flv_ch_Rkin
  integer , allocatable , dimension(:,:) :: flv_ch_RRkin
!
  integer :: num_flv_irr_LO,num_flv_irr_NLO_V,num_flv_irr_NLO_R, &
             num_flv_irr_NNLO_VV,num_flv_irr_NNLO_RV,num_flv_irr_NNLO_RR
!
  real(kind(1d0)) , allocatable , dimension(:) :: weights_Bkin
  real(kind(1d0)) , allocatable , dimension(:) :: weights_Rkin
  real(kind(1d0)) , allocatable , dimension(:) :: weights_RRkin
!
  type(eqv_cont) , allocatable , dimension(:) :: eqv_LO
  type(eqv_cont) , allocatable , dimension(:) :: eqv_NLO_V
  type(eqv_cont) , allocatable , dimension(:) :: eqv_NNLO_VV
  type(eqv_cont) , allocatable , dimension(:) :: eqv_NLO_R
  type(eqv_cont) , allocatable , dimension(:) :: eqv_NNLO_RR
  type(eqv_cont) , allocatable , dimension(:) :: eqv_NNLO_RV
!
contains
!
subroutine init_subprocesses()
use process
use input
implicit none
!
!
  integer :: istat
  integer :: iproc
!
!
  call init_processes('gen ')
!
  allocate(flv_LO(nleg_born,num_flv_LO),stat=istat)
  if (istat.ne.0) then
    print *,"Problem occured during allocating flv_LO..."
    stop
  end if
  allocate(flv_NLO_V(nleg_born,num_flv_NLO_V),stat=istat)
  if (istat.ne.0) then
    print *,"Problem occured during allocating flv_NLO_V..."
    stop
  end if
  allocate(flv_NLO_R(nleg_born + 1,num_flv_NLO_R),stat=istat)
  if (istat.ne.0) then
    print *,"Problem occured during allocating flv_NLO_R..."
    stop
  end if
  allocate(flv_NNLO_VV(nleg_born,num_flv_NNLO_VV),stat=istat)
  if (istat.ne.0) then
    print *,"Problem occured during allocating flv_NNLO_VV..."
    stop
  end if
  allocate(flv_NNLO_RV(nleg_born + 1,num_flv_NNLO_RV),stat=istat)
  if (istat.ne.0) then
    print *,"Problem occured during allocating flv_NNLO_RV..."
    stop
  end if
  allocate(flv_NNLO_RR(nleg_born + 2,num_flv_NNLO_RR),stat=istat)
  if (istat.ne.0) then
    print *,"Problem occured during allocating flv_NNLO_RR..."
    stop
  end if
!
  call init_processes('load')
!
  print *,"The following subprocesses were loaded: "
  write(*,'(A,I4)') "Number of subprocesses at the Born level: ", &
                    num_flv_LO
  do iproc=1,num_flv_LO
    write(*,'(20(I3,1x))') flv_LO(:,iproc)
  end do
!
  write(*,'(A,I4)') "Number of subprocesses at the NLO-V level: ", &
                    num_flv_NLO_V
  do iproc=1,num_flv_NLO_V
    write(*,'(20(I3,1x))') flv_NLO_V(:,iproc)
  end do
!
  write(*,'(A,I4)') "Number of subprocesses at the NLO-R level: ", &
                    num_flv_NLO_R
  do iproc=1,num_flv_NLO_R
    write(*,'(20(I3,1x))') flv_NLO_R(:,iproc)
  end do
!
  write(*,'(A,I4)') "Number of subprocesses at the NNLO-VV level: ", &
                    num_flv_NNLO_VV
  do iproc=1,num_flv_NNLO_VV
    write(*,'(20(I3,1x))') flv_NNLO_VV(:,iproc)
  end do
!
  write(*,'(A,I4)') "Number of subprocesses at the NNLO-RV level: ", &
                    num_flv_NNLO_RV
  do iproc=1,num_flv_NNLO_RV
    write(*,'(20(I3,1x))') flv_NNLO_RV(:,iproc)
  end do
!
  write(*,'(A,I4)') "Number of subprocesses at the NNLO-RR level: ", &
                    num_flv_NNLO_RR
  do iproc=1,num_flv_NNLO_RR
    write(*,'(20(I3,1x))') flv_NNLO_RR(:,iproc)
  end do
!
! Before we go to check the possible relations between the various
! subprocesses we offer the possibility to check the SMEs:
  if (nnloinput("#debugsme").eq.1) call DebugSME
!
! We try to obtain relations between the various contributions:
  call calc_relations()
!
! After setting up the possible relations between various SMEs
! we can check the insertion operators:
  if (nnloinput("#debugI1op").eq.1) call DebugI1op
!
end subroutine init_subprocesses
!
! This routine takes all the contributions and tries to deduce 
! relationships within:
subroutine calc_relations()
use flags
use process
use momenta
use particles
use scales
implicit none
!
!
  integer :: iproc
  integer :: imom,ileg
  integer :: istat
  integer , parameter :: nmom = 10
!
  real(kind(1d0)) :: SME
  real(kind(1d0)) , dimension(:,:) , allocatable :: numrels
  type(mom) , dimension(:) , allocatable :: parr_tmp
  type(particle) , dimension(:) , allocatable :: parts_tmp
!
  interface
    subroutine CalcB(parts,smeB)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: parts
      real(kind(1d0)) , intent(out) :: smeB
!
    end subroutine CalcB
!
    subroutine CalcV(parts,smeV,smeVLaurent)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: parts
      real(kind(1d0)) , intent(out) :: smeV
      real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: smeVLaurent
!
    end subroutine CalcV
!
    subroutine CalcVddim(parts,smeV,smeVLaurent)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: parts
      real(kind(1d0)) , intent(out) :: smeV
      real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: smeVLaurent
!
    end subroutine CalcVddim
!
    subroutine CalcVV(parts,smeVV,smeVVLaurent)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: parts
      real(kind(1d0)) , intent(out) :: smeVV
      real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: smeVVLaurent
!
    end subroutine CalcVV
!
    subroutine CalcR(parts,smeR)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: parts
      real(kind(1d0)) , intent(out) :: smeR
!
    end subroutine CalcR
!
    subroutine CalcRV(parts,smeRV,smeRVLaurent)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: parts
      real(kind(1d0)) , intent(out) :: smeRV
      real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: smeRVLaurent
!
    end subroutine CalcRV
!
    subroutine CalcRR(parts,smeRR)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: parts
      real(kind(1d0)) , intent(out) :: smeRR
!
    end subroutine CalcRR
  end interface
!
! *********
! **Born:**
! *********
  if (flg_LO) then
! We consider the Born subprocesses:
! We allocate an array to hold the momenta:
    allocate(parr_tmp(nleg_born),stat=istat)
    if (istat.ne.0) then; print *,"Problem with parr_tmp"; stop; end if
! We allocate an array for the particles:
    allocate(parts_tmp(nleg_born),stat=istat)
    if (istat.ne.0) then; print *,"Problem with parts_tmp"; stop; end if
! We allocate an array which holds the numerical values of the SMEs
! in the PS points:
    allocate(numrels(nmom,num_flv_LO))
    if (istat.ne.0) then; print *,"Problem with numrels"; stop; end if
! We also have to allocate an array which will hold the equivalence 
! relations:
    allocate(eqv_LO(num_flv_LO),stat=istat)
    if (istat.ne.0) then; print *,"Problem with eqv_LO"; stop; end if
! We consider nmom different PS points:
    do imom=1,nmom
      call gen_ranmom(nleg_born,parr_tmp)
!
!      print *,"The generated momenta: "
!      do ileg=1,nleg_born
!        call printmom(parr_tmp(ileg))
!      end do
!
! For each PS point we calculate all the contributions:
      do iproc=1,num_flv_LO
!        write(*,'(A,1x,I4,1x,20(I3,1x))') "iproc, flv: ",iproc, &
!                                          flv_LO(:,iproc)
        call CreateParts(parr_tmp,flv_LO(:,iproc),parts_tmp)
        call CalcB(parts_tmp,SME)
! And store them in an array accordingly:
        numrels(imom,iproc) = SME
      end do
    end do
! We try to deduce relations between various contributions:
    call AnalyzeConts(num_flv_LO,flv_LO,numrels,eqv_LO)
    deallocate(parr_tmp,parts_tmp,numrels)
    Print *,"We found the following relations at the Born level: "
    call PrintRelations(num_flv_LO,flv_LO,eqv_LO)
  end if
!
! ************
! **Virtual:**
! ************
  if (flg_NLO_V) then
! We consider the Virtual subprocesses:
! We allocate an array to hold the momenta:
    allocate(parr_tmp(nleg_born),stat=istat)
    if (istat.ne.0) then; print *,"Problem with parr_tmp"; stop; end if
! We allocate an array for the particles:
    allocate(parts_tmp(nleg_born),stat=istat)
    if (istat.ne.0) then; print *,"Problem with parts_tmp"; stop; end if
! We allocate an array which holds the numerical values of the SMEs
! in the PS points:
    allocate(numrels(nmom,num_flv_NLO_V))
    if (istat.ne.0) then; print *,"Problem with numrels"; stop; end if
! We also have to allocate an array which will hold the equivalence 
! relations:
    allocate(eqv_NLO_V(num_flv_NLO_V),stat=istat)
    if (istat.ne.0) then; print *,"Problem with eqv_NLO_V"; stop; end if
! We consider nmom different PS points:
    do imom=1,nmom
      call gen_ranmom(nleg_born,parr_tmp)
! When the first set of momenta is generated we calculate a fictitious scale,
! this is needed since the virtual part has a \mu_R dependence beside
! of the one coming from \alpha_S:
      if (imom.eq.1) mur = sqrt(2d0*parr_tmp(1)*parr_tmp(2))
!
!      print *,"The generated momenta: "
!      do ileg=1,nleg_born
!        call printmom(parr_tmp(ileg))
!      end do
!
! For each PS point we calculate all the contributions:
      do iproc=1,num_flv_NLO_V
!        write(*,'(A,1x,I4,1x,20(I3,1x))') "iproc, flv: ",iproc, &
!                                          flv_NLO_V(:,iproc)
        call CreateParts(parr_tmp,flv_NLO_V(:,iproc),parts_tmp)
        call CalcV(parts_tmp,SME)
!        call CalcVddim(parts_tmp,SME)
! And store them in an array accordingly:
        numrels(imom,iproc) = SME
      end do
    end do
! We try to deduce relations between various contributions:
    call AnalyzeConts(num_flv_NLO_V,flv_NLO_V,numrels,eqv_NLO_V)
    deallocate(parr_tmp,parts_tmp,numrels)
    Print *,"We found the following relations for the virtual: "
    call PrintRelations(num_flv_NLO_V,flv_NLO_V,eqv_NLO_V)
  end if
!
! *******
! **VV:**
! *******
  if (flg_NNLO_VV) then
! We consider the Virtual-Virtual subprocesses:
! We allocate an array to hold the momenta:
    allocate(parr_tmp(nleg_born),stat=istat)
    if (istat.ne.0) then; print *,"Problem with parr_tmp"; stop; end if
! We allocate an array for the particles:
    allocate(parts_tmp(nleg_born),stat=istat)
    if (istat.ne.0) then; print *,"Problem with parts_tmp"; stop; end if
! We allocate an array which holds the numerical values of the SMEs
! in the PS points:
    allocate(numrels(nmom,num_flv_NNLO_VV))
    if (istat.ne.0) then; print *,"Problem with numrels"; stop; end if
! We also have to allocate an array which will hold the equivalence 
! relations:
    allocate(eqv_NNLO_VV(num_flv_NNLO_VV),stat=istat)
    if (istat.ne.0) then; print *,"Problem with eqv_NNLO_VV"; stop; end if
! We consider nmom different PS points:
    do imom=1,nmom
      call gen_ranmom(nleg_born,parr_tmp)
! When the first set of momenta is generated we calculate a fictitious scale,
! this is needed since the virtual part has a \mu_R dependence beside
! of the one coming from \alpha_S:
      if (imom.eq.1) mur = sqrt(2d0*parr_tmp(1)*parr_tmp(2))
!
!      print *,"The generated momenta: "
!      do ileg=1,nleg_born
!        call printmom(parr_tmp(ileg))
!      end do
!
! For each PS point we calculate all the contributions:
      do iproc=1,num_flv_NNLO_VV
!        write(*,'(A,1x,I4,1x,20(I3,1x))') "iproc, flv: ",iproc, &
!                                          flv_NNLO_VV(:,iproc)
        call CreateParts(parr_tmp,flv_NNLO_VV(:,iproc),parts_tmp)
        call CalcVV(parts_tmp,SME)
! And store them in an array accordingly:
        numrels(imom,iproc) = SME
      end do
    end do
! We try to deduce relations between various contributions:
    call AnalyzeConts(num_flv_NNLO_VV,flv_NNLO_VV,numrels,eqv_NNLO_VV)
    deallocate(parr_tmp,parts_tmp,numrels)
    Print *,"We found the following relations for the virtual: "
    call PrintRelations(num_flv_NNLO_VV,flv_NNLO_VV,eqv_NNLO_VV)
  end if
!
! *********
! **Real:**
! *********
  if (flg_NLO_R.or.flg_NNLO_RRA1_A1) then
! We consider the Real subprocesses:
! We allocate an array to hold the momenta:
    allocate(parr_tmp(nleg_born+1),stat=istat)
    if (istat.ne.0) then; print *,"Problem with parr_tmp"; stop; end if
! We allocate an array for the particles:
    allocate(parts_tmp(nleg_born+1),stat=istat)
    if (istat.ne.0) then; print *,"Problem with parts_tmp"; stop; end if
! We allocate an array which holds the numerical values of the SMEs
! in the PS points:
    allocate(numrels(nmom,num_flv_NLO_R))
    if (istat.ne.0) then; print *,"Problem with numrels"; stop; end if
! We also have to allocate an array which will hold the equivalence 
! relations:
    allocate(eqv_NLO_R(num_flv_NLO_R),stat=istat)
    if (istat.ne.0) then; print *,"Problem with eqv_NLO_R"; stop; end if
! We consider nmom different PS points:
    do imom=1,nmom
      call gen_ranmom(nleg_born+1,parr_tmp)
!
!      print *,"The generated momenta: "
!      do ileg=1,nleg_born+1
!        call printmom(parr_tmp(ileg))
!      end do
!
! For each PS point we calculate all the contributions:
      do iproc=1,num_flv_NLO_R
!        write(*,'(A,1x,I4,1x,20(I3,1x))') "iproc, flv: ",iproc, &
!                                          flv_NLO_R(:,iproc)
        call CreateParts(parr_tmp,flv_NLO_R(:,iproc),parts_tmp)
        call CalcR(parts_tmp,SME)
! And store them in an array accordingly:
        numrels(imom,iproc) = SME
      end do
    end do
! We try to deduce relations between various contributions:
    call AnalyzeConts(num_flv_NLO_R,flv_NLO_R,numrels,eqv_NLO_R)
    deallocate(parr_tmp,parts_tmp,numrels)
    Print *,"We found the following relations for the real: "
    call PrintRelations(num_flv_NLO_R,flv_NLO_R,eqv_NLO_R)
  end if
!
! *****************
! **Real-Virtual:**
! *****************
  if (flg_NNLO_RV) then
! We consider the Real-Virtual subprocesses:
! We allocate an array to hold the momenta:
    allocate(parr_tmp(nleg_born+1),stat=istat)
    if (istat.ne.0) then; print *,"Problem with parr_tmp"; stop; end if
! We allocate an array for the particles:
    allocate(parts_tmp(nleg_born+1),stat=istat)
    if (istat.ne.0) then; print *,"Problem with parts_tmp"; stop; end if
! We allocate an array which holds the numerical values of the SMEs
! in the PS points:
    allocate(numrels(nmom,num_flv_NNLO_RV))
    if (istat.ne.0) then; print *,"Problem with numrels"; stop; end if
! We also have to allocate an array which will hold the equivalence 
! relations:
    allocate(eqv_NNLO_RV(num_flv_NNLO_RV),stat=istat)
    if (istat.ne.0) then; print *,"Problem with eqv_NNLO_RV"; stop; end if
! We consider nmom different PS points:
    do imom=1,nmom
      call gen_ranmom(nleg_born+1,parr_tmp)
! When the first set of momenta is generated we calculate a fictitious scale,
! this is needed since the virtual part has a \mu_R dependence beside
! of the one coming from \alpha_S:
      if (imom.eq.1) mur = sqrt(2d0*parr_tmp(1)*parr_tmp(2))
!
!      print *,"The generated momenta: "
!      do ileg=1,nleg_born+1
!        call printmom(parr_tmp(ileg))
!      end do
!
! For each PS point we calculate all the contributions:
      do iproc=1,num_flv_NNLO_RV
!        write(*,'(A,1x,I4,1x,20(I3,1x))') "iproc, flv: ",iproc, &
!                                          flv_NNLO_RV(:,iproc)
        call CreateParts(parr_tmp,flv_NNLO_RV(:,iproc),parts_tmp)
        call CalcRV(parts_tmp,SME)
! And store them in an array accordingly:
        numrels(imom,iproc) = SME
      end do
    end do
! We try to deduce relations between various contributions:
    call AnalyzeConts(num_flv_NNLO_RV,flv_NNLO_RV,numrels,eqv_NNLO_RV)
    deallocate(parr_tmp,parts_tmp,numrels)
    Print *,"We found the following relations for the real-virtual: "
    call PrintRelations(num_flv_NNLO_RV,flv_NNLO_RV,eqv_NNLO_RV)
  end if
!
! *********
! ** RR: **
! *********
  if (flg_NNLO_RR) then
! We consider the RR subprocesses:
! We allocate an array to hold the momenta:
    allocate(parr_tmp(nleg_born+2),stat=istat)
    if (istat.ne.0) then; print *,"Problem with parr_tmp"; stop; end if
! We allocate an array for the particles:
    allocate(parts_tmp(nleg_born+2),stat=istat)
    if (istat.ne.0) then; print *,"Problem with parts_tmp"; stop; end if
! We allocate an array which holds the numerical values of the SMEs
! in the PS points:
    allocate(numrels(nmom,num_flv_NNLO_RR))
    if (istat.ne.0) then; print *,"Problem with numrels"; stop; end if
! We also have to allocate an array which will hold the equivalence 
! relations:
    allocate(eqv_NNLO_RR(num_flv_NNLO_RR),stat=istat)
    if (istat.ne.0) then; print *,"Problem with eqv_NLO_RR"; stop; end if
! We consider nmom different PS points:
    do imom=1,nmom
      call gen_ranmom(nleg_born+2,parr_tmp)
!
!      print *,"The generated momenta: "
!      do ileg=1,nleg_born+2
!        call printmom(parr_tmp(ileg))
!      end do
!
! For each PS point we calculate all the contributions:
      do iproc=1,num_flv_NNLO_RR
!        write(*,'(A,1x,I4,1x,20(I3,1x))') "iproc, flv: ",iproc, &
!                                          flv_NNLO_RR(:,iproc)
        call CreateParts(parr_tmp,flv_NNLO_RR(:,iproc),parts_tmp)
        call CalcRR(parts_tmp,SME)
! And store them in an array accordingly:
        numrels(imom,iproc) = SME
      end do
    end do
! We try to deduce relations between various contributions:
    call AnalyzeConts(num_flv_NNLO_RR,flv_NNLO_RR,numrels,eqv_NNLO_RR)
    deallocate(parr_tmp,parts_tmp,numrels)
    Print *,"We found the following relations for the double real: "
    call PrintRelations(num_flv_NNLO_RR,flv_NNLO_RR,eqv_NNLO_RR)
  end if
!
! Finally we create the possible irreducible subprocesses for the PS
! generator:
  call SetupWeights
!
end subroutine calc_relations
!
subroutine SetupWeights
use process
use flags
use particles
implicit none
!
!
  integer :: iproc,jproc,ipart
  integer :: eqv_old
  integer :: numproc
  integer :: istat
  real(kind(1d0)) :: coeff_old
  real(kind(1d0)) :: weight
!
  type(mom) , dimension(:) , allocatable :: p_tmp
  type(particle) , dimension(:) , allocatable :: parts_tmp
  real(kind(1d0)) :: SME1,SME2
!
  real(kind(1d0)) , parameter :: eps = 1d-8
!
  interface
    subroutine CalcB(parts,smeB)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: parts
      real(kind(1d0)) , intent(out) :: smeB
!
    end subroutine CalcB
!
    subroutine CalcV(parts,smeV,smeVLaurent)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: parts
      real(kind(1d0)) , intent(out) :: smeV
      real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: smeVLaurent
!
    end subroutine CalcV
!
    subroutine CalcR(parts,smeR)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: parts
      real(kind(1d0)) , intent(out) :: smeR
!
    end subroutine CalcR
!
    subroutine CalcRV(parts,smeRV,smeRVLaurent)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: parts
      real(kind(1d0)) , intent(out) :: smeRV
      real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: smeRVLaurent
!
    end subroutine CalcRV
!
    subroutine CalcRR(parts,smeRR)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: parts
      real(kind(1d0)) , intent(out) :: smeRR
!
    end subroutine CalcRR
  end interface
!
  num_flv_ch_Bkin = 0
  num_flv_ch_Rkin = 0
  num_flv_ch_RRkin = 0
!
  num_flv_irr_LO = 0
  num_flv_irr_NLO_V = 0
  num_flv_irr_NLO_R = 0
  num_flv_irr_NNLO_VV = 0
  num_flv_irr_NNLO_RV = 0
  num_flv_irr_NNLO_RR = 0
!
! For each multiplicity we set up the arrays for the irreducible 
! subprocesses for the multichannel PS generator:
!=======================================================================
! Born-like final state (LO + NLO-V + NNLO-VV):
!=======================================================================
  if (flg_LO) then
    do iproc=1,num_flv_LO
! New configuration is found:
      if (eqv_LO(iproc)%eqv.eq.-1) then
! Counting all irreducibile contributions with born kinematics:
        num_flv_ch_Bkin = num_flv_ch_Bkin + 1
! counting only those coming from the Born:
        num_flv_irr_LO = num_flv_irr_LO + 1
      end if
    end do
  end if
  if (flg_NLO_V) then
! The number of irreducible Virtual subprocesses coincides with
! the number of irreducible Born ones, although we can have 
! additional channels opened up via loop corrections.
! Start with the number taken from the Born:
    num_flv_irr_NLO_V = num_flv_irr_LO
    do iproc=1,num_flv_NLO_V
! The ordering of subprocesses should be the same as for the Born
      if (iproc.le.num_flv_LO) then
        do ipart=1,size(flv_LO,1)
          if (flv_LO(ipart,iproc).ne.flv_NLO_V(ipart,iproc)) then
            print *,"Error..."
            print *,"The ordering for Born and Virtual subprocesses are different..."
            print *,"iproc: ",iproc
            print *,"Born: "
            call PrintSubproc(flv_LO(:,iproc))
            print *,"Virtual: "
            call PrintSubproc(flv_NLO_V(:,iproc))
            stop
          end if
        end do
! The ordering is the same there is no need to register the 
! irreducible subprocess for the phase space generator 
! if the Born is also calculated:
        if (flg_LO) cycle
      end if
! New configuration is found:
      if (eqv_NLO_V(iproc)%eqv.eq.-1) then
        num_flv_ch_Bkin = num_flv_ch_Bkin + 1
! counting only those coming from the Born:
        num_flv_irr_NLO_V = num_flv_irr_NLO_V + 1
      end if
    end do
  end if
  if (flg_NNLO_VV) then
! In principle the number of VV channels is the same as the 
! number of Virtual ones, since if a channel is opened it is
! already opened by the virtual.
    num_flv_irr_NNLO_VV = num_flv_irr_NLO_V
    do iproc=1,num_flv_NNLO_VV
! The ordering of subprocesses should be the same as for the Virtual
      if (iproc.le.num_flv_NLO_V) then
        do ipart=1,size(flv_NLO_V,1)
          if (flv_NLO_V(ipart,iproc).ne.flv_NNLO_VV(ipart,iproc)) then
            print *,"Error..."
            print *,"The ordering for Born and Double-Virtual ",&
                    "subprocesses are different..."
            print *,"iproc: ",iproc
            print *,"Virtual: "
            call PrintSubproc(flv_NLO_V(:,iproc))
            print *,"Double-Virtual: "
            call PrintSubproc(flv_NNLO_VV(:,iproc))
            stop
          end if
        end do
! The ordering is the same there is no need to register the 
! irreducible subprocess for the phase space generator 
! if the Virtual is also calculated:
        if (flg_NLO_V) cycle
      end if
! New configuration is found:
      if (eqv_NNLO_VV(iproc)%eqv.eq.-1) then
        num_flv_ch_Bkin = num_flv_ch_Bkin + 1
! counting only those coming from the Born:
        num_flv_irr_NNLO_VV = num_flv_irr_NNLO_VV + 1
      end if
    end do
  end if
!=======================================================================
! Real-like final state (NLO_R + NNLO_RV):
!=======================================================================
! Before counting down irreducible subprocesses for the Real 
! contribution a check is made whether a subprocess is present or not
! which is reducible for the Real but not for the Real-Virtual. This
! can easily happen when having loop contributions:
! Loop through all Real subprocesses if the Real-Virtual is also about
! to be calculated:
  if (flg_NNLO_RV.and.(flg_NLO_R.or.flg_NNLO_RRA1_A1)) then
    do iproc=1,num_flv_NLO_R
! The ordering of subprocesses should be the same as for the Real
      do ipart=1,size(flv_NLO_R,1)
        if (flv_NLO_R(ipart,iproc).ne.flv_NNLO_RV(ipart,iproc)) then
          print *,"Error..."
          print *,"The ordering for Real and Real-Virtual ", &
                  "subprocesses are different..."
          print *,"iproc: ",iproc
          print *,"Real: "
          call PrintSubproc(flv_NLO_R(:,iproc))
          print *,"Real-Virtual: "
          call PrintSubproc(flv_NNLO_RV(:,iproc))
          stop
        end if
      end do
! If a given subprocess is irreducible for the Real-Virtual but was
! reducible for the Real we make it irreducible even for the Real:
      if ((eqv_NLO_R(iproc)%eqv.ne.-1).and. &
          (eqv_NNLO_RV(iproc)%eqv.eq.-1)) then
! For the Real it was reduced to eqv_old with some coefficient:
        eqv_old = eqv_NLO_R(iproc)%eqv
        coeff_old = eqv_NLO_R(iproc)%coeff
! Made it irreducible:
        eqv_NLO_R(iproc)%eqv = eqv_NNLO_RV(iproc)%eqv
! Irreducible the coefficient is irrelevant:
        eqv_NLO_R(iproc)%coeff = 1d0
        print *,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!", &
                "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        print *,"A Real subprocess is made irreducible since it is ", &
                "irreducible for the Real-Virtual:"
        call PrintSubproc(flv_NLO_R(:,iproc))
        print *,"Originally it was reduced to: "
        call PrintSubproc(flv_NLO_R(:,eqv_old))
        print *,"With reduction coefficient ",coeff_old
! There can be other subprocesses which are reduced to the iproc-th
! one:
        do jproc=iproc+1,num_flv_NLO_R
! Found one among the Real-Virtual subprocesses:
          if (eqv_NNLO_RV(jproc)%eqv.eq.iproc) then
! This subprocess should be a reducible one even for the Real and 
! should have been expressed with the subprocess eqv_old
            if (eqv_NLO_R(jproc)%eqv.ne.eqv_old) then
              print *,"Error in ordering of subprocesses..."
              stop
            end if
! Change the reduction to the subprocess sitting at position iproc:
            eqv_NLO_R(jproc)%eqv = iproc
! The coefficient is changed to coincide with the coefficient taken
! from the Real-Virtual subprocess:
            eqv_NLO_R(jproc)%coeff = eqv_NNLO_RV(jproc)%coeff
            print *,"A reducible Real subprocess: "
            call PrintSubproc(flv_NLO_R(:,jproc))
            print *,"is reassigned to: "
            call PrintSubproc(flv_NLO_R(:,iproc))
            print *,"with coefficient of: ",eqv_NLO_R(jproc)%coeff
            print *,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!", &
                    "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
! This is risky as hell hence better check it:
! Allocate arrays to hold momenta and particles:
            allocate(p_tmp(nleg_born+1), &
                     parts_tmp(nleg_born+1), &
                     stat=istat)
            if (istat.ne.0) then
              print *,"Problem during allocating p_tmp..."
              stop
            end if
! Generate a random set of momenta:
            call gen_ranmom(nleg_born+1,p_tmp)
! Create particles for the irreducible subprocess:
            call CreateParts(p_tmp,flv_NLO_R(:,iproc),parts_tmp)
! Calculate the SME for it:
            call CalcR(parts_tmp,SME1)
! Create particles for the reducible one:
            call CreateParts(p_tmp,flv_NLO_R(:,jproc),parts_tmp)
! Calculate the SME for it:
            call CalcR(parts_tmp,SME2)
            if (abs(eqv_NNLO_RV(jproc)%coeff - SME2/SME1).gt.eps) then
              print *,"Coefficient mismatch in Real <-> Real-Virtual..."
              print *,"The original Real irreducible subprocess: "
              call PrintSubproc(flv_NLO_R(:,eqv_old))
              print *,"The irreducible Real-Virtual subprocess: "
              call PrintSubproc(flv_NNLO_RV(:,iproc))
              print *,"The reducible Real subprocess: "
              call PrintSubproc(flv_NLO_R(:,jproc))
              print *,"The reduction coefficient for Real-Virtual: ", &
                      eqv_NNLO_RV(jproc)%coeff
              print *,"The ratio obtained for Real subprocesses ",iproc,jproc
              print *,SME2/SME1
              stop
            end if
! Deallocate what is not needed:
            deallocate(p_tmp,parts_tmp)
          end if
        end do
      end if
    end do
  end if
  if (flg_NLO_R.or.flg_NNLO_RRA1_A1) then
    do iproc=1,num_flv_NLO_R
! New configuration is found:
      if (eqv_NLO_R(iproc)%eqv.eq.-1) then
! Counting all irreducibile contributions with born kinematics:
        num_flv_ch_Rkin = num_flv_ch_Rkin + 1
! counting only those coming from the Born:
        num_flv_irr_NLO_R = num_flv_irr_NLO_R + 1
      end if
    end do
  end if
  if (flg_NNLO_RV) then
! The number of irreducible Real-Virtual subprocesses coincides with
! the number of irreducible Real ones, although we can have 
! additional channels opened up via loop corrections.
! Start with the number taken from the Real:
    num_flv_irr_NNLO_RV = num_flv_irr_NLO_R
    do iproc=1,num_flv_NNLO_RV
      if (iproc.le.num_flv_NLO_R) then
! The ordering is the same there is no need to register the 
! irreducible subprocess for the phase space generator 
! if the Real is also calculated:
        if (flg_NLO_R.or.flg_NNLO_RRA1_A1) cycle
      end if
! New configuration is found:
      if (eqv_NNLO_RV(iproc)%eqv.eq.-1) then
        num_flv_ch_Rkin = num_flv_ch_Rkin + 1
! counting only those coming from the Real:
        num_flv_irr_NNLO_RV = num_flv_irr_NNLO_RV + 1
      end if
    end do
  end if
!=======================================================================
! Double-Real-like final state (NNLO_RR):
!=======================================================================
  if (flg_NNLO_RR) then
    do iproc=1,num_flv_NNLO_RR
! New configuration is found:
      if (eqv_NNLO_RR(iproc)%eqv.eq.-1) then
! Counting all irreducibile contributions with born kinematics:
        num_flv_ch_RRkin = num_flv_ch_RRkin + 1
! counting only those coming from the Born:
        num_flv_irr_NNLO_RR = num_flv_irr_NNLO_RR + 1
      end if
    end do
  end if
!
! We allocate an array which will hold the irreducible 
! subprocesses:
!=======================================================================
! Born-like kinematics:
!=======================================================================
  if (flg_LO.or.flg_NLO_V.or.flg_NNLO_VV) then
    allocate(flv_ch_Bkin(nleg_born,num_flv_ch_Bkin),stat=istat)
    if (istat.ne.0) then
      print *,"Problem with flv_ch_Bkin..."
      stop
    end if
  end if
!=======================================================================
! Real-like kinematics:
!=======================================================================
  if (flg_NLO_R.or.flg_NNLO_RV) then
    allocate(flv_ch_Rkin(nleg_born+1,num_flv_ch_Rkin),stat=istat)
    if (istat.ne.0) then
      print *,"Problem with flv_ch_Rkin..."
      stop
    end if
  end if
!=======================================================================
! Double-Real-like kinematics:
!=======================================================================
  if (flg_NNLO_RR) then
    allocate(flv_ch_RRkin(nleg_born+2,num_flv_ch_RRkin),stat=istat)
    if (istat.ne.0) then
      print *,"Problem with flv_ch_RRkin..."
      stop
    end if
  end if
! We allocate an array which will hold the weights:
!=======================================================================
! Born-like kinematics:
!=======================================================================
  if (flg_LO.or.flg_NLO_V.or.flg_NNLO_VV) then
    allocate(weights_Bkin(num_flv_ch_Bkin),stat=istat)
    if (istat.ne.0) then
      print *,"Problem with weights_Bkin..."
      stop
    end if
  end if
!=======================================================================
! Real-like kinematics:
!=======================================================================
  if (flg_NLO_R.or.flg_NNLO_RV) then
    allocate(weights_Rkin(num_flv_ch_Rkin),stat=istat)
    if (istat.ne.0) then
      print *,"Problem with weights_Rkin..."
      stop
    end if
  end if
!=======================================================================
! Double-Real-like kinematics:
!=======================================================================
  if (flg_NNLO_RR) then
    allocate(weights_RRkin(num_flv_ch_RRkin),stat=istat)
    if (istat.ne.0) then
      print *,"Problem with weights_RRkin..."
      stop
    end if
  end if
!
! We fill up the new arrays:
!=======================================================================
! Born:
!=======================================================================
  if (flg_LO) then
    numproc = 0
    do iproc=1,num_flv_LO
      if (eqv_LO(iproc)%eqv.eq.-1) then
        numproc = numproc + 1
        if (numproc.gt.size(flv_ch_Bkin,2)) then
          print *,"More irreducible processes in LO..."
          stop
        end if
        flv_ch_Bkin(:,numproc) = flv_LO(:,iproc)
! Obtaining the total weight associated to this subprocess:
        weight = 1d0
        do jproc=iproc+1,num_flv_LO
          if (eqv_LO(jproc)%eqv.eq.iproc) then
            weight = weight + eqv_LO(jproc)%coeff
          end if
        end do
        weights_Bkin(numproc) = weight
      end if
    end do
  end if
!
!=======================================================================
! Virtual:
!=======================================================================
!
  if (flg_NLO_V) then
    numproc = 0
    do iproc=1,num_flv_NLO_V
      if (eqv_NLO_V(iproc)%eqv.eq.-1) then
        numproc = numproc + 1
        if (numproc.gt.size(flv_ch_Bkin,2)) then
          print *,"More irreducible processes in the NLO V..."
          stop
        end if
! We have two possibilities: the subprocess is already added
! by the Born hence we only check for consistency:
        if (numproc.le.num_flv_irr_LO) then
          do ipart=1,size(flv_NLO_V,1)
            if ((flv_ch_Bkin(ipart,numproc).ne.flv_NLO_V(ipart,iproc)).and. &
                flg_LO) then
              print *,"Inconsistency detected betweeen Born and virtual"
              print *,"subprocesses in the equivalence array..."
              print *,"iproc: ",iproc
              print *,"Born: "
              call PrintSubproc(flv_ch_Bkin(:,numproc))
              print *,"Virtual: "
              call PrintSubproc(flv_NLO_V(:,iproc))
              stop
            end if
          end do
        else
! Otherwise we have to include the subprocess:
          flv_ch_Bkin(:,numproc) = flv_LO(:,iproc)
! Obtaining the total weight associated to this subprocess:
          weight = 1d0
          do jproc=iproc+1,num_flv_NLO_V
            if (eqv_NLO_V(jproc)%eqv.eq.iproc) then
              weight = weight + eqv_NLO_V(jproc)%coeff
            end if
          end do
          weights_Bkin(numproc) = weight
        end if
      end if
    end do
  end if
!
!=======================================================================
! Double-Virtual:
!=======================================================================
!
  if (flg_NNLO_VV) then
    numproc = 0
    do iproc=1,num_flv_NNLO_VV
      if (eqv_NNLO_VV(iproc)%eqv.eq.-1) then
        numproc = numproc + 1
        if (numproc.gt.size(flv_ch_Bkin,2)) then
          print *,"More irreducible processes in the NNLO VV..."
          stop
        end if
! We have two possibilities: the subprocess is already added
! by the Born or the virtual hence we only check for consistency:
        if (numproc.le.num_flv_irr_NLO_V) then
          do ipart=1,size(flv_NNLO_VV,1)
            if ((flv_ch_Bkin(ipart,numproc).ne.flv_NNLO_VV(ipart,iproc)).and. &
                (flg_LO.or.flg_NLO_V)) then
              print *,"Inconsistency detected betweeen Born and doublevirtual"
              print *,"subprocesses in the equivalence array..."
              print *,"iproc: ",iproc
              print *,"Virtual: "
              call PrintSubproc(flv_ch_Bkin(:,numproc))
              print *,"Double-Virtual: "
              call PrintSubproc(flv_NNLO_VV(:,iproc))
              stop
            end if
          end do
        else
! Otherwise we have to include the subprocess:
          if (.not.flg_NLO_V) then
            flv_ch_Bkin(:,numproc) = flv_LO(:,iproc)
          else
            flv_ch_Bkin(:,numproc) = flv_NLO_V(:,iproc)
          end if
! Obtaining the total weight associated to this subprocess:
          weight = 1d0
          do jproc=iproc+1,num_flv_NNLO_VV
            if (eqv_NNLO_VV(jproc)%eqv.eq.iproc) then
              weight = weight + eqv_NNLO_VV(jproc)%coeff
            end if
          end do
          weights_Bkin(numproc) = weight
        end if
      end if
    end do
  end if
!
!=======================================================================
! Real:
!=======================================================================
!
  if (flg_NLO_R.or.flg_NNLO_RRA1_A1) then
    numproc = 0
    do iproc=1,num_flv_NLO_R
      if (eqv_NLO_R(iproc)%eqv.eq.-1) then
        numproc = numproc + 1
        if (numproc.gt.size(flv_ch_Rkin,2)) then
          print *,"More irreducible processes in NLO R..."
          stop
        end if
        flv_ch_Rkin(:,numproc) = flv_NLO_R(:,iproc)
! Obtaining the total weight associated to this subprocess:
        weight = 1d0
        do jproc=iproc+1,num_flv_NLO_R
! It is possible that there is a simple relationship between
! certain real-emission subprocesses, but this does not hold 
! for the real-virtual contribution, if this happens the increment
! in the weight has to be avoided:
          if (flg_NNLO_RV) then
! The ordering of subprocesses is strictly the same for the real and the
! real-virtual, hence if the jproc-th subprocess for the real-virtual
! is an independent one we are not allowed to increment the weight:
            if (eqv_NNLO_RV(jproc)%eqv.eq.-1) cycle
          end if
          if (eqv_NLO_R(jproc)%eqv.eq.iproc) then
            weight = weight + eqv_NLO_R(jproc)%coeff
          end if
        end do
        weights_Rkin(numproc) = weight
      end if
    end do
  end if
!
!=======================================================================
! Real-Virtual:
!=======================================================================
!
  if (flg_NNLO_RV) then
    numproc = 0
    do iproc=1,num_flv_NNLO_RV
      if (eqv_NNLO_RV(iproc)%eqv.eq.-1) then
        numproc = numproc + 1
        if (numproc.gt.size(flv_ch_Rkin,2)) then
          print *,"More irreducible processes in the NNLO RV..."
          stop
        end if
! We have two possibilities: the subprocess is already added
! by the Real hence we only check for consistency:
        if (numproc.le.num_flv_irr_NLO_R) then
          do ipart=1,size(flv_NNLO_RV,1)
            if ((flv_ch_Rkin(ipart,numproc).ne.flv_NNLO_RV(ipart,iproc)).and. &
                (flg_NLO_R.or.flg_NNLO_RRA1_A1)) then
              print *,"Inconsistency detected betweeen Real and Real-Virtual"
              print *,"subprocesses in the equivalence array..."
              print *,"iproc: ",iproc
              print *,"Real: "
              call PrintSubproc(flv_ch_Rkin(:,numproc))
              print *,"Real-Virtual: "
              call PrintSubproc(flv_NNLO_RV(:,iproc))
              stop
            end if
          end do
        else
! Otherwise we have to include the subprocess:
          flv_ch_Rkin(:,numproc) = flv_NLO_R(:,iproc)
! Obtaining the total weight associated to this subprocess:
          weight = 1d0
          do jproc=iproc+1,num_flv_NNLO_RV
            if (eqv_NNLO_RV(jproc)%eqv.eq.iproc) then
              weight = weight + eqv_NNLO_RV(jproc)%coeff
            end if
          end do
          weights_Rkin(numproc) = weight
        end if
      end if
    end do
  end if
!
!=======================================================================
! Double-Real:
!=======================================================================
!
  if (flg_NNLO_RR) then
    numproc = 0
    do iproc=1,num_flv_NNLO_RR
      if (eqv_NNLO_RR(iproc)%eqv.eq.-1) then
        numproc = numproc + 1
        if (numproc.gt.size(flv_ch_RRkin,2)) then
          print *,"More irreducible processes in NNLO RR..."
          stop
        end if
        flv_ch_RRkin(:,numproc) = flv_NNLO_RR(:,iproc)
! Obtaining the total weight associated to this subprocess:
        weight = 1d0
        do jproc=iproc+1,num_flv_NNLO_RR
          if (eqv_NNLO_RR(jproc)%eqv.eq.iproc) then
            weight = weight + eqv_NNLO_RR(jproc)%coeff
          end if
        end do
        weights_RRkin(numproc) = weight
      end if
    end do
  end if
!
! Some printout:
!
  print *,"For the MC integration the following subprocesses are used:"
  if (flg_LO.or.flg_NLO_V.or.flg_NNLO_VV) then
    print *,"Born kinematics: "
    call PrintChannels(num_flv_ch_Bkin,flv_ch_Bkin,weights_Bkin)
  end if
  if (flg_NLO_R.or.flg_NNLO_RV) then
    print *,"Real kinematics: "
    call PrintChannels(num_flv_ch_Rkin,flv_ch_Rkin,weights_Rkin)
  end if
  if (flg_NNLO_RR) then
    print *,"Double Real kinematics: "
    call PrintChannels(num_flv_ch_RRkin,flv_ch_RRkin,weights_RRkin)
  end if
!
end subroutine SetupWeights
!
subroutine AnalyzeConts(num_flv,flv_arr,numrels,eqv_arr)
use particles
implicit none
!
  integer , intent(in) :: num_flv
  integer , dimension(:,:) , intent(in) :: flv_arr
  real(kind(1d0)) , dimension(:,:) , intent(in) :: numrels
  type(eqv_cont) , dimension(:) , intent(out) :: eqv_arr
!
  integer :: imom,iproc,jproc
  logical :: passed
  real(kind(1d0)) :: ratio,ratio_tmp
  real(kind(1d0)) , parameter :: eps = 1d-8
!
! We initialize the array holding the equivalence information to
! zero:
  eqv_arr(:)%eqv = 0
  eqv_arr(:)%coeff = 0d0
!
! If we have only one process there is no relation to be found:
  if (num_flv.eq.1) then
    eqv_arr(1)%eqv = -1
    eqv_arr(1)%coeff = 1d0
    return
  end if
!
! We loop over all the subprocesses but the last:
  do iproc=1,num_flv
! If it is already assigned a status we skip it:
    if (eqv_arr(iproc)%eqv.ne.0) cycle
! We consider all the other subprocesses which are higher in the
! list:
    do jproc=iproc+1,num_flv
! We consider all PS points:
      passed = .true.
      do imom=1,size(numrels,1)
        if (imom.eq.1) then
          ratio = numrels(imom,jproc) / numrels(imom,iproc)
        end if
        ratio_tmp = numrels(imom,jproc) / numrels(imom,iproc)
        if (abs(ratio - ratio_tmp).gt.eps) then
! Only considering the case when the ratio is one:
!        if (abs(ratio - 1d0).gt.eps) then
! If there is one for which the ratio is not correct we can
! leave the loop since there is no relation:
          passed = .false.
          exit
        end if
      end do
! If a relation is found we have to fill in the ratio: 
      if (passed) then
        eqv_arr(jproc)%eqv   = iproc
        eqv_arr(jproc)%coeff = ratio
      end if
    end do
! In any case the iproc subprocess is irreducible:
    eqv_arr(iproc)%eqv   = -1
    eqv_arr(iproc)%coeff = 1d0
  end do
!
!  print *,"The equivalence relations at the end are: "
!  do iproc=1,num_flv
!    call PrintSubproc(flv_arr(:,iproc))
!    print *,"eqv_arr: ",eqv_arr(iproc)%eqv,eqv_arr(iproc)%coeff
!  end do
!
end subroutine AnalyzeConts
!
subroutine PrintRelations(num_flv,flv_arr,eqv_arr)
use particles
use utils
implicit none
!
  integer , intent(in) :: num_flv
  integer , dimension(:,:) , intent(in) :: flv_arr
  type(eqv_cont) , dimension(:) , intent(in) :: eqv_arr
!
  integer :: iproc,ipart
!
  do iproc=1,num_flv
    do ipart=1,size(flv_arr,1)
      write(*,'(A,1x)',advance='no') &
        ConvertFromPDG(flv_arr(ipart,iproc))
      if (ipart.eq.2) write(*,'(A)',advance='no') "-> " 
    end do
    if (eqv_arr(iproc)%eqv.eq.-1) then
      write(*,'(A)') " : Irreducible"
      cycle
    end if
    write(*,'(A,F6.4,1x)',advance='no') " ~ ",eqv_arr(iproc)%coeff
    do ipart=1,size(flv_arr,1)
      write(*,'(A,1x)',advance='no') &
        ConvertFromPDG(flv_arr(ipart,eqv_arr(iproc)%eqv))
      if (ipart.eq.2) write(*,'(A)',advance='no') "-> " 
    end do
    write(*,*)
  end do
!
end subroutine PrintRelations
!
subroutine PrintChannels(num_ch,flv_ch,weights)
use utils
implicit none
!
  integer , intent(in) :: num_ch
  integer , dimension(:,:) , intent(in) :: flv_ch
  real(kind(1d0)) , dimension(:) , intent(in) :: weights
!
  integer :: iproc,ipart
!
  do iproc=1,num_ch
    do ipart=1,size(flv_ch,1)
      write(*,'(A,1x)',advance='no') &
        ConvertFromPDG(flv_ch(ipart,iproc))
      if (ipart.eq.2) write(*,'(A)',advance='no') "-> " 
    end do
    write(*,'(A,F8.4)') " * ",weights(iproc)
  end do
!
end subroutine PrintChannels
!
end module subprocesses
!
module coupling
implicit none
!
  integer :: orderas
  logical :: runningaEM
  real(kind(1d0)) :: alphaEM0,alphaS0
  real(kind(1d0)) :: alphaS
  real(kind(1d0)) :: alphaEM
!
! These should be set by the user, these are the
! exponents of alphaS and alphaEM at Born level:
  integer :: border_as
  integer :: border_aem
!
contains
!
subroutine init_couplings
use QCDparams
use input
implicit none
!
!
  real(kind(1d0)) :: muscale
!
  real(kind(1d0)) , external :: genericxlambdL, &
                                genericxlambdNL, &
                                genericxlambdNNL
!
  runningaEM = .false.
! If runaem is present in the input card \alpha_{EM}
! is running otherwise kept fixed
  if (nnloinput("#runaem").gt.0d0) runningaEM = .true.
!
! We set alphaEM:
  alphaEM0 = 1d0/137d0
! The inverse can be supplied through the input card:
  if (nnloinput("#alphaem").gt.0d0) then
    alphaEM0 = 1d0 / nnloinput("#alphaem")
  end if
!
! This is an NNLO code hence alphaS should be 3-loop running:
  orderas = 3
! Though, through input can be modified:
  if (nnloinput("#orderas").gt.0) then
    orderas = nnloinput("#orderas")
  end if
!
  print *,"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
  write(*,'(a,1x,I2,1x,a)') &
  " %%%%%%%%%%%%%%%%% alphaS loop order:",orderas,"%%%%%%%%%%%%%%%%%%%%"
  write(*,'(a,1x,F10.6,1x,a)') &
  " %%%%%%%%%%%%%%%%% 1/alphaEM:",1d0/alphaEM0,"%%%%%%%%%%%%%%%%%%%%"
  print *,"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
!
  print *,"alphaEM: ",alphaEM0
!
! If a line with alphas is in the input card the alphaS at muscale
! is read from there than \Lambda_{QCD} is calculated:
  if (nnloinput("#alphas").gt.0) then
    alphaS0 = nnloinput("alphas")
! We have to read in the central scale too, which cannot be dynamical:
    muscale = nnloinput("muscale")
    if (muscale.lt.0) then
      print *,"with fixed alphaS the muscale can only be fixed too!"
      stop
    end if
    write(*,'(a,F6.4,a,E13.6)') "alphaS is ",alphaS0," at ",muscale
    print *,muscale
! We calculate \Lambda_{QCD}:
    if (orderas.eq.1) then
      qcd_lambda = genericxlambdL(alphaS0,muscale,qcd_nf)
    elseif (orderas.eq.2) then
      qcd_lambda = genericxlambdNL(alphaS0,muscale,qcd_nf)
    elseif (orderas.eq.3) then
      qcd_lambda = genericxlambdNNL(alphaS0,muscale,qcd_nf)
    else
      print *,"Fixed alphaS is supplied, but the corresponding ", &
              "lambdaQCD routine is not present..."
      stop
    end if
  end if
!
end subroutine init_couplings
!
subroutine calc_couplings
use QCDparams
use scales
implicit none
!
!
  real(kind(1d0)) , external :: calc_alphas
!
! We calculate the couplings:
  alphas = calc_alphas(orderas,mur**2,qcd_lambda,qcd_nf)
  if (runningaem) then
    print *,"Running aEM is not implemented yet..."
    stop
  else
    alphaEM = alphaEM0
  end if
!
end subroutine calc_couplings
!
end module coupling
!
module misc
implicit none
!
!
interface CastHelicityToLorentz
  module procedure CastHelicityToLorentz_munu
  module procedure CastHelicityToLorentz_q1q2_munu
end interface CastHelicityToLorentz
!
interface CalcSpinorProd
  module procedure CalcSpinorProd_mom
  module procedure CalcSpinorProd_part
end interface CalcSpinorProd
!
interface CalcPolVector
  module procedure CalcPolVector_mom
  module procedure CalcPolVector_part
end interface CalcPolVector
!
interface Contract
  module procedure Contract_mom_mom
  module procedure Contract_part_part
  module procedure Contract_mom_part
  module procedure Contract_part_mom
  module procedure Contract_M_M
end interface Contract
!
interface LorentzLambda
  module procedure LorentzLambda_mom_mom
  module procedure LorentzLambda_part_part
end interface LorentzLambda
!
interface SwapColor
  module procedure SwapColorij_mom_mom
  module procedure SwapColorij_part_part
  module procedure SwapColorijkl_mom_mom
  module procedure SwapColorijkl_part_part
end interface SwapColor
!
interface SwapColorLaurent
  module procedure SwapColorijLaurent_mom_mom
  module procedure SwapColorijLaurent_part_part
  module procedure SwapColorijklLaurent_mom_mom
  module procedure SwapColorijklLaurent_part_part
end interface SwapColorLaurent
!
interface CalcFluxFact
  module procedure CalcFluxFact_mom
  module procedure CalcFluxFact_part
end interface CalcFluxFact
!
contains
!
subroutine CastHelicityToLorentz_munu(k,q,M2ij,M2munu)
use momenta
use particles
use utils
implicit none
!
  type(particle) , intent(in) :: k,q
  complex(kind(1d0)) , dimension(2,2) , intent(in) :: M2ij
  real(kind(1d0)) , dimension(0:3,0:3) , intent(out) :: M2munu
!
  integer :: mu,nu
  type(cmom) :: eplus,eminus
  complex(kind(1d0)) :: tmp
  complex(kind(1d0)) , dimension(0:3) :: eplus_arr,eminus_arr
!
!
! For debugging:
!  real(kind(1d0)) , dimension(0:3,0:3) :: amunu,bmunu
!  real(kind(1d0)) , dimension(0:3) :: k_arr,q_arr
!
! Sanity check on the particle: it should be a gluon:
  if (k%flv.ne.0) then
    print *,"Error in CastHelicityToLorentz_munu..."
    print *,"We expected a gluon but we got: ",ConvertFromPDG(k%flv)
  end if
!
! We have to construct the polarization vectors for p:
  call CalcPolVector(k,q,eplus,eminus)
!
  eplus_arr(0) = eplus%E  ; eplus_arr(1) = eplus%px
  eplus_arr(2) = eplus%py ; eplus_arr(3) = eplus%pz
  eminus_arr(0) = eminus%E  ; eminus_arr(1) = eminus%px
  eminus_arr(2) = eminus%py ; eminus_arr(3) = eminus%pz
!
! M2ij(1,1) = {|M^2|}^{++}
! M2ij(1,2) = {|M^2|}^{+-}
! M2ij(2,1) = {|M^2|}^{-+}
! M2ij(2,2) = {|M^2|}^{--}
!
  M2munu = 0d0
  do mu=0,3
    do nu=0,3
      tmp = eminus_arr(mu)*M2ij(1,1)*eplus_arr(nu)  &
          + eminus_arr(mu)*M2ij(1,2)*eminus_arr(nu) &
          + eplus_arr(mu)*M2ij(2,1)*eplus_arr(nu)   &
          + eplus_arr(mu)*M2ij(2,2)*eminus_arr(nu)
      M2munu(mu,nu) = M2munu(mu,nu) + tmp
    end do
  end do
!
!  print *,"Checking normalization: "
!  print *,eplus*eplus
!  print *,eminus*eminus
!  print *,eplus*eminus
!  print *,k%p*eplus
!  print *,q%p*eplus
!  print *,k%p*eminus
!  print *,q%p*eminus
!
!   amunu = 0d0
!   bmunu = 0d0
!   k_arr(0) = k%p%E
!   k_arr(1) = k%p%px
!   k_arr(2) = k%p%py
!   k_arr(3) = k%p%pz
!   q_arr(0) = q%p%E
!   q_arr(1) = q%p%px
!   q_arr(2) = q%p%py
!   q_arr(3) = q%p%pz
!!
!   do mu=0,3
!     do nu=0,3
!       if (mu.eq.nu) then
!         if (mu.eq.0) bmunu(mu,nu) = bmunu(mu,nu) - 1d0
!         if (mu.ne.0) bmunu(mu,nu) = bmunu(mu,nu) + 1d0
!       end if
!       amunu(mu,nu) = amunu(mu,nu) &
!             + eplus_arr(mu)*conjg(eplus_arr(nu))   &
!             + eminus_arr(mu)*conjg(eminus_arr(nu))
!       bmunu(mu,nu) = bmunu(mu,nu) &
!             + (k_arr(mu)*q_arr(nu) + k_arr(nu)*q_arr(mu))/(k%p*q%p)
!     end do
!   end do
!   print *,"amunu: ",amunu
!   print *,"bmunu: ",bmunu
!
end subroutine CastHelicityToLorentz_munu
!
! This routine is the same as the previous one but it takes
! two reference momenta instead of one:
subroutine CastHelicityToLorentz_q1q2_munu(k,q1,q2,M2ij,M2munu)
use momenta
use particles
use utils
implicit none
!
  type(particle) , intent(in) :: k
  type(mom) , intent(in) :: q1,q2
  complex(kind(1d0)) , dimension(2,2) , intent(in) :: M2ij
  real(kind(1d0)) , dimension(0:3,0:3) , intent(out) :: M2munu
!
  integer :: mu,nu
  type(cmom) :: eplus1,eminus1,eplus2,eminus2
  complex(kind(1d0)) :: tmp
  complex(kind(1d0)) , dimension(0:3) :: eplus1_arr,eminus1_arr
  complex(kind(1d0)) , dimension(0:3) :: eplus2_arr,eminus2_arr
!
!
! Sanity check on the particle: it should be a gluon:
  if (k%flv.ne.0) then
    print *,"Error in CastHelicityToLorentz_q1q2_munu..."
    print *,"We expected a gluon but we got: ",ConvertFromPDG(k%flv)
  end if
!
! We have to construct the polarization vectors for p:
  call CalcPolVector(k%p,q1,eplus1,eminus1)
  call CalcPolVector(k%p,q2,eplus2,eminus2)
!
  eplus1_arr(0) = eplus1%E  ; eplus1_arr(1) = eplus1%px
  eplus1_arr(2) = eplus1%py ; eplus1_arr(3) = eplus1%pz
  eminus1_arr(0) = eminus1%E  ; eminus1_arr(1) = eminus1%px
  eminus1_arr(2) = eminus1%py ; eminus1_arr(3) = eminus1%pz
!
  eplus2_arr(0) = eplus2%E  ; eplus2_arr(1) = eplus2%px
  eplus2_arr(2) = eplus2%py ; eplus2_arr(3) = eplus2%pz
  eminus2_arr(0) = eminus2%E  ; eminus2_arr(1) = eminus2%px
  eminus2_arr(2) = eminus2%py ; eminus2_arr(3) = eminus2%pz
!
! M2ij(1,1) = {|M^2|}^{++}
! M2ij(1,2) = {|M^2|}^{+-}
! M2ij(2,1) = {|M^2|}^{-+}
! M2ij(2,2) = {|M^2|}^{--}
!
  M2munu = 0d0
  do mu=0,3
    do nu=0,3
      tmp = eminus1_arr(mu)*M2ij(1,1)*eplus2_arr(nu)  &
          + eminus1_arr(mu)*M2ij(1,2)*eminus2_arr(nu) &
          + eplus1_arr(mu)*M2ij(2,1)*eplus2_arr(nu)   &
          + eplus1_arr(mu)*M2ij(2,2)*eminus2_arr(nu)
      M2munu(mu,nu) = M2munu(mu,nu) + tmp
    end do
  end do
!
end subroutine CastHelicityToLorentz_q1q2_munu
!
! This routine calculates the spinor products <12> and [12]
! for two real momenta the definition is in accord with
! arXiv:hep-ph/9601359
subroutine CalcSpinorProd_mom(p1,p2,a,b)
use momenta
implicit none
!
  type(mom) , intent(in) :: p1,p2
  complex(kind(1d0)) , intent(out) :: a,b
!
  real(kind(1d0)) :: p1t,p2t
  real(kind(1d0)) :: p1plus,p2plus
  real(kind(1d0)) :: p1minus,p2minus
  complex(kind(1d0)) :: expphase1p,expphase1m
  complex(kind(1d0)) :: expphase2p,expphase2m
!
!
! p1 has positive energy:
  if (p1%E.gt.0d0) then
    p1plus  = p1%E + p1%pz
    p1minus = p1%E - p1%pz
!
    p1t = sqrt(p1%px**2 + p1%py**2)
!
! the pt for p1 is non-zero:
    if (p1t.gt.0d0) then
      expphase1p = cmplx(p1%px, p1%py)/p1t
      expphase1m = cmplx(p1%px,-p1%py)/p1t
! the pt for p1 is zero:
    else
      expphase1p = cmplx(1d0,0d0)
      expphase1m = cmplx(1d0,0d0)
    end if
! p2 has positive energy:
    if (p2%E.gt.0d0) then
      p2plus  = p2%E + p2%pz
      p2minus = p2%E - p2%pz
!
      p2t = sqrt(p2%px**2 + p2%py**2)
! 
! the pt for p2 is non-zero:
      if (p2t.gt.0d0) then
        expphase2p = cmplx(p2%px, p2%py)/p2t
        expphase2m = cmplx(p2%px,-p2%py)/p2t
! the pt for p2 is zero:
      else
        expphase2p = cmplx(1d0,0d0)
        expphase2m = cmplx(1d0,0d0)
      end if
! the spinor products for the p_1^0 > 0 p_2^0 > 0 case:
      a =  sqrt(p1minus*p2plus)*expphase1p &
        -  sqrt(p1plus*p2minus)*expphase2p
      a = -a
      b = -sqrt(p1minus*p2plus)*expphase1m &
        +  sqrt(p1plus*p2minus)*expphase2m
      b = -b
! p2 has negative energy:
    else
      p2plus  = -(p2%E + p2%pz)
      p2minus = -(p2%E - p2%pz)
!
      p2t = sqrt(p2%px**2 + p2%py**2)
! 
! the pt for p2 is non-zero:
      if (p2t.gt.0d0) then
        expphase2p = -cmplx(p2%px, p2%py)/p2t
        expphase2m = -cmplx(p2%px,-p2%py)/p2t
! the pt for p2 is zero:
      else
        expphase2p = cmplx(1d0,0d0)
        expphase2m = cmplx(1d0,0d0)
      end if
! the spinor products for the p_1^0 > 0 p_2^0 < 0 case:
      a =  sqrt(p1minus*p2plus)*expphase1p &
        -  sqrt(p1plus*p2minus)*expphase2p
      a = cmplx(0d0,1d0) * a
      b = -sqrt(p1minus*p2plus)*expphase1m &
        +  sqrt(p1plus*p2minus)*expphase2m
      b = cmplx(0d0,1d0) * b
    end if
! p1 has negative energy:
  else
    p1plus  = -(p1%E + p1%pz)
    p1minus = -(p1%E - p1%pz)
!
    p1t = sqrt(p1%px**2 + p1%py**2)
!
! the pt for p1 is non-zero:
    if (p1t.gt.0d0) then
      expphase1p = -cmplx(p1%px, p1%py)/p1t
      expphase1m = -cmplx(p1%px,-p1%py)/p1t
! the pt for p1 is zero:
    else
      expphase1p = cmplx(1d0,0d0)
      expphase1m = cmplx(1d0,0d0)
    end if
! p2 has positive energy:
    if (p2%E.gt.0d0) then
      p2plus  = p2%E + p2%pz
      p2minus = p2%E - p2%pz
!
      p2t = sqrt(p2%px**2 + p2%py**2)
! 
! the pt for p2 is non-zero:
      if (p2t.gt.0d0) then
        expphase2p = cmplx(p2%px, p2%py)/p2t
        expphase2m = cmplx(p2%px,-p2%py)/p2t
! the pt for p2 is zero:
      else
        expphase2p = cmplx(1d0,0d0)
        expphase2m = cmplx(1d0,0d0)
      end if
! the spinor products for the p_1^0 < 0 p_2^0 > 0 case:
      a =  sqrt(p1minus*p2plus)*expphase1p &
        -  sqrt(p1plus*p2minus)*expphase2p
      a = cmplx(0d0,1d0) * a
      b = -sqrt(p1minus*p2plus)*expphase1m &
        +  sqrt(p1plus*p2minus)*expphase2m
      b = cmplx(0d0,1d0) * b
! p2 has negative energy:
    else
      p2plus  = -(p2%E + p2%pz)
      p2minus = -(p2%E - p2%pz)
!
      p2t = sqrt(p2%px**2 + p2%py**2)
! 
! the pt for p2 is non-zero:
      if (p2t.gt.0d0) then
        expphase2p = -cmplx(p2%px, p2%py)/p2t
        expphase2m = -cmplx(p2%px,-p2%py)/p2t
! the pt for p2 is zero:
      else
        expphase2p = cmplx(1d0,0d0)
        expphase2m = cmplx(1d0,0d0)
      end if
! the spinor products for the p_1^0 < 0 p_2^0 < 0 case:
      a =  sqrt(p1minus*p2plus)*expphase1p &
        -  sqrt(p1plus*p2minus)*expphase2p
      b = -sqrt(p1minus*p2plus)*expphase1m &
        +  sqrt(p1plus*p2minus)*expphase2m
    end if
  end if
!
end subroutine CalcSpinorProd_mom
!
subroutine CalcSpinorProd_part(p1,p2,a,b)
use particles
use momenta
implicit none
!
  type(particle) , intent(in) :: p1,p2
  complex(kind(1d0)) , intent(out) :: a,b
!
!
!
  call CalcSpinorProd_mom(p1%p,p2%p,a,b)
!
end subroutine CalcSpinorProd_part
!
subroutine CalcPolVector_mom(p,q,eplus,eminus)
use momenta
implicit none
!
  type(mom) , intent(in) :: p,q
  type(cmom) , intent(out) :: eplus,eminus
!
  real(kind(1d0)) :: pplus,pminus,qplus,qminus
  real(kind(1d0)) :: pt,qt
  complex(kind(1d0)) :: expphasepp,expphasepm
  complex(kind(1d0)) :: expphaseqp,expphaseqm
  complex(kind(1d0)) :: Apq,Bpq
!
!
! p^0 > 0:
  if (p%E.gt.0d0) then
    pplus  = p%E + p%pz
    pminus = p%E - p%pz
!
    pt = sqrt(p%px**2 + p%py**2)
!
    if (pt.gt.0d0) then
      expphasepp = cmplx(p%px, p%py)/pt
      expphasepm = cmplx(p%px,-p%py)/pt
    else
      expphasepp = cmplx(1d0,0d0)
      expphasepm = cmplx(1d0,0d0)
    end if
! q^0 > 0:
    if (q%E.gt.0d0) then
      qplus  = q%E + q%pz
      qminus = q%E - q%pz
!
      qt = sqrt(q%px**2 + q%py**2)
!
      if (qt.gt.0d0) then
        expphaseqp = cmplx(q%px, q%py)/qt
        expphaseqm = cmplx(q%px,-q%py)/qt
      else
        expphaseqp = cmplx(1d0,0d0)
        expphaseqm = cmplx(1d0,0d0)
      end if
! calculating Apq = <pq> and Bpq = [pq]
      Apq =  sqrt(abs(pminus*qplus))*expphasepp &
          -  sqrt(abs(pplus*qminus))*expphaseqp
      Bpq = -sqrt(abs(pminus*qplus))*expphasepm &
          +  sqrt(abs(pplus*qminus))*expphaseqm
!
      eplus%E  = sqrt(abs(pplus*qplus)) + sqrt(abs(pminus*qminus)) &
               * expphasepm * expphaseqp
      eplus%px = sqrt(abs(pplus*qminus))*expphaseqp &
               + sqrt(abs(pminus*qplus))*expphasepm
      eplus%py = cmplx(0d0,1d0)*sqrt(abs(pminus*qplus))*expphasepm &
               - cmplx(0d0,1d0)*sqrt(abs(pplus*qminus))*expphaseqp
      eplus%pz = sqrt(abs(pplus*qplus)) - sqrt(abs(pminus*qminus)) &
               * expphasepm * expphaseqp
      eplus = -eplus/sqrt(2d0)/Apq
!
      eminus%E  = conjg(eplus%E)
      eminus%px = conjg(eplus%px)
      eminus%py = conjg(eplus%py)
      eminus%pz = conjg(eplus%pz)
! q^0 < 0:
    else
      qplus  = -q%E - q%pz
      qminus = -q%E + q%pz
!
      qt = sqrt(q%px**2 + q%py**2)
!
      if (qt.gt.0d0) then
        expphaseqp = -cmplx(q%px, q%py)/qt
        expphaseqm = -cmplx(q%px,-q%py)/qt
      else
        expphaseqp = cmplx(1d0,0d0)
        expphaseqm = cmplx(1d0,0d0)
      end if
! calculating Apq = <pq> and Bpq = [pq]
      Apq =  sqrt(abs(pminus*qplus))*expphasepp &
          -  sqrt(abs(pplus*qminus))*expphaseqp
      Bpq = -sqrt(abs(pminus*qplus))*expphasepm &
          +  sqrt(abs(pplus*qminus))*expphaseqm
!
      eplus%E  = sqrt(abs(pplus*qplus)) + sqrt(abs(pminus*qminus)) &
               * expphasepm * expphaseqp
      eplus%px = sqrt(abs(pplus*qminus))*expphaseqp &
               + sqrt(abs(pminus*qplus))*expphasepm
      eplus%py = cmplx(0d0,1d0)*sqrt(abs(pminus*qplus))*expphasepm &
               - cmplx(0d0,1d0)*sqrt(abs(pplus*qminus))*expphaseqp
      eplus%pz = sqrt(abs(pplus*qplus)) - sqrt(abs(pminus*qminus)) &
               * expphasepm * expphaseqp
      eplus = -eplus/sqrt(2d0)/Apq
!
      eminus%E  = conjg(eplus%E)
      eminus%px = conjg(eplus%px)
      eminus%py = conjg(eplus%py)
      eminus%pz = conjg(eplus%pz)
    end if
! p^0 < 0:
  else
    print *,"CalcPolVector_mom..."
    print *,"p^0 < 0 case is not implemented yet..."
    stop
  end if
!
end subroutine CalcPolVector_mom
!
subroutine CalcPolVector_part(p,q,eplus,eminus)
use momenta
use particles
implicit none
!
  type(particle) , intent(in) :: p,q
  type(cmom) , intent(out) :: eplus,eminus
!
!
!
  call CalcPolVector_mom(p%p,q%p,eplus,eminus)
!
end subroutine CalcPolVector_part
!
! The input is given as: p^\mu, M^{\mu\nu} and q^\nu:
function Contract_mom_mom(p,M,q) result(ans)
use momenta
implicit none
!
  real(kind(1d0)) :: ans
!
  type(mom) , intent(in) :: p,q
  real(kind(1d0)) , dimension(0:3,0:3) , intent(in) :: M
!
  integer :: mu,nu
  real(kind(1d0)) , dimension(0:3) :: p_arr,q_arr
!
!
! To make the contraction it is better to convert the momenta into
! arrays and we also have to bring down the indices:
  p_arr(0) = p%E
  p_arr(1) = -p%px
  p_arr(2) = -p%py
  p_arr(3) = -p%pz
!
  q_arr(0) = q%E
  q_arr(1) = -q%px
  q_arr(2) = -q%py
  q_arr(3) = -q%pz
!
  ans = 0d0
!
  do mu=0,3
    do nu=0,3
      ans = ans + p_arr(mu)*M(mu,nu)*q_arr(nu)
    end do
  end do
!
end function Contract_mom_mom
!
function Contract_part_part(p,M,q) result(ans)
use momenta
use particles
implicit none
!
  real(kind(1d0)) :: ans
!
  type(particle) , intent(in) :: p,q
  real(kind(1d0)) , dimension(0:3,0:3) :: M
!
!
!
  ans = Contract_mom_mom(p%p,M,q%p)
!
end function Contract_part_part
!
function Contract_part_mom(p,M,q) result(ans)
use momenta
use particles
implicit none
!
  real(kind(1d0)) :: ans
!
  type(particle) , intent(in) :: p
  type(mom) , intent(in) :: q
  real(kind(1d0)) , dimension(0:3,0:3) :: M
!
!
!
  ans = Contract_mom_mom(p%p,M,q)
!
end function Contract_part_mom
!
function Contract_mom_part(p,M,q) result(ans)
use momenta
use particles
implicit none
!
  real(kind(1d0)) :: ans
!
  type(mom) , intent(in) :: p
  type(particle) , intent(in) :: q
  real(kind(1d0)) , dimension(0:3,0:3) :: M
!
!
!
  ans = Contract_mom_mom(p,M,q%p)
!
end function Contract_mom_part
!
function Contract_M_M(A,B) result(C)
implicit none
!
  real(kind(1d0)) , dimension(0:3,0:3) :: C
!
  real(kind(1d0)) , dimension(0:3,0:3) , intent(in) :: A,B
!
  integer :: i,j
!
  do i=0,3
    do j=0,3
      C(i,j) = A(0,i)*B(j,0) - A(1,i)*B(j,1) &
             - A(2,i)*B(j,2) - A(3,i)*B(j,3)
    end do
  end do
!
end function Contract_M_M
!
subroutine LorentzLambda_mom_mom(K,Khat,p,q)
use momenta
implicit none
!
  type(mom) , intent(in) :: K,Khat,p
  type(mom) , intent(out) :: q
!
!
!
  call Lambda4vec(K,Khat,p,q)
!
end subroutine LorentzLambda_mom_mom
!
subroutine LorentzLambda_part_part(K,Khat,p,q)
use particles
implicit none
!
  type(mom) , intent(in) :: K,Khat
  type(particle) , intent(in) :: p
  type(particle) , intent(out) :: q
!
!
!
  call Lambda4vec(K,Khat,p%p,q%p)
  q%flv = p%flv
!
end subroutine LorentzLambda_part_part
!
subroutine SwapColorij_mom_mom(nleg,pin,pout,Bij)
use momenta
implicit none
!
  integer , intent(in) :: nleg
  type(mom) , dimension(:) , intent(in) :: pin,pout
  real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
!
  integer :: ileg,jleg
  real(kind(1d0)) , dimension(nleg) :: tmp
!
!
  do ileg=1,nleg-1
    do jleg=ileg+1,nleg
      if ((pin(ileg)%E.eq.pout(jleg)%E).and. &
          (pin(ileg)%px.eq.pout(jleg)%px).and. &
          (pin(ileg)%py.eq.pout(jleg)%py).and. &
          (pin(ileg)%pz.eq.pout(jleg)%pz)) then
        tmp = Bij(ileg,:)
        Bij(ileg,:) = Bij(jleg,:)
        Bij(jleg,:) = tmp
        tmp = Bij(:,ileg)
        Bij(:,ileg) = Bij(:,jleg)
        Bij(:,jleg) = tmp
      end if
    end do
  end do
!
end subroutine SwapColorij_mom_mom
!
subroutine SwapColorij_part_part(nleg,pin,pout,Bij)
use particles
implicit none
!
  integer , intent(in) :: nleg
  type(particle) , dimension(:) , intent(in) :: pin,pout
  real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
!
  integer :: ileg
  type(mom) , dimension(nleg) :: pmomin,pmomout
!
!
  do ileg=1,nleg
    pmomin(ileg)  = pin(ileg)%p
    pmomout(ileg) = pout(ileg)%p
  end do
!
  call SwapColorij_mom_mom(nleg,pmomin,pmomout,Bij)
!
end subroutine SwapColorij_part_part
!
subroutine SwapColorijLaurent_mom_mom(nleg,pin,pout,BijLaurent)
use particles
implicit none
!
  integer , intent(in) :: nleg
  type(mom) , dimension(:) , intent(in) :: pin,pout
  real(kind(1d0)) , dimension(:,:,-4:) , intent(inout) :: BijLaurent
!
  integer :: i
!
!
  do i=-4,2
    call SwapColorij_mom_mom(nleg,pin,pout,BijLaurent(:,:,i))
  end do
!
end subroutine SwapColorijLaurent_mom_mom
!
subroutine SwapColorijLaurent_part_part(nleg,pin,pout,BijLaurent)
use particles
implicit none
!
  integer , intent(in) :: nleg
  type(particle) , dimension(:) , intent(in) :: pin,pout
  real(kind(1d0)) , dimension(:,:,-4:) , intent(inout) :: BijLaurent
!
!
!
  call SwapColorijLaurent_mom_mom(nleg,pin(:)%p,pout(:)%p,BijLaurent)
!
end subroutine SwapColorijLaurent_part_part
!
subroutine SwapColorijklLaurent_mom_mom(nleg,pin,pout,BijklLaurent)
use particles
implicit none
!
  integer , intent(in) :: nleg
  type(mom) , dimension(:) , intent(in) :: pin,pout
  real(kind(1d0)) , dimension(:,:,:,:,-4:) , intent(inout) :: BijklLaurent
!
  integer :: i
!
!
  do i=-4,2
    call SwapColorijkl_mom_mom(nleg,pin,pout,BijklLaurent(:,:,:,:,i))
  end do
!
end subroutine SwapColorijklLaurent_mom_mom
!
subroutine SwapColorijklLaurent_part_part(nleg,pin,pout,BijklLaurent)
use particles
implicit none
!
  integer , intent(in) :: nleg
  type(particle) , dimension(:) , intent(in) :: pin,pout
  real(kind(1d0)) , dimension(:,:,:,:,-4:) , intent(inout) :: BijklLaurent
!
!
!
  call SwapColorijklLaurent_mom_mom(nleg,pin(:)%p,pout(:)%p,BijklLaurent)
!
end subroutine SwapColorijklLaurent_part_part
!
! When the input array has rank 4 we can face a doubly
! color-correlated SME or a spin- and color-correlated SME, too:
subroutine SwapColorRank4_part_part(nleg,pin,pout,Bijkl)
use particles
implicit none
!
  integer , intent(in) :: nleg
  type(particle) , dimension(:) , intent(in) :: pin,pout
  real(kind(1d0)) , dimension(:,:,:,:) , intent(inout) :: Bijkl
!
  integer :: ileg
  type(mom) , dimension(nleg) :: pmomin,pmomout
!
!
  do ileg=1,nleg
    pmomin(ileg)  = pin(ileg)%p
    pmomout(ileg) = pout(ileg)%p
  end do
!
! If the lower bound in the first two dimension is zero
! we have a simultaneously spin- and color-correlated SME:
  if ((lbound(Bijkl,1).eq.0).and.(lbound(bijkl,2).eq.0)) then
    call SwapColormunuij_mom_mom(nleg,pmomin,pmomout,Bijkl)
  else
    call SwapColorijkl_mom_mom(nleg,pmomin,pmomout,Bijkl)
  end if
!
end subroutine SwapColorRank4_part_part
!
! When the input array has rank 4 we can face a doubly
! color-correlated SME or a spin- and color-correlated SME, too:
subroutine SwapColorRank4_mom_mom(nleg,pin,pout,Bijkl)
use particles
implicit none
!
  integer , intent(in) :: nleg
  type(mom) , dimension(:) , intent(in) :: pin,pout
  real(kind(1d0)) , dimension(:,:,:,:) , intent(inout) :: Bijkl
!
  integer :: ileg
!
!
! If the lower bound in the first two dimension is zero
! we have a simultaneously spin- and color-correlated SME:
  if ((lbound(Bijkl,1).eq.0).and.(lbound(bijkl,2).eq.0)) then
    call SwapColormunuij_mom_mom(nleg,pin,pout,Bijkl)
  else
    call SwapColorijkl_mom_mom(nleg,pin,pout,Bijkl)
  end if
!
end subroutine SwapColorRank4_mom_mom
!
subroutine SwapColorijkl_mom_mom(nleg,pin,pout,Bijkl)
use momenta
implicit none
!
  integer , intent(in) :: nleg
  type(mom) , dimension(:) , intent(in) :: pin,pout
  real(kind(1d0)) , dimension(:,:,:,:) , intent(inout) :: Bijkl
!
  integer :: ileg,jleg
  integer :: i,j,k
  real(kind(1d0)) :: tmp
!
!
  do ileg=1,nleg-1
    do jleg=ileg+1,nleg
      if ((pin(ileg)%E.eq.pout(jleg)%E).and. &
          (pin(ileg)%px.eq.pout(jleg)%px).and. &
          (pin(ileg)%py.eq.pout(jleg)%py).and. &
          (pin(ileg)%pz.eq.pout(jleg)%pz)) then
        do i=1,nleg
          do j=1,nleg
            do k=1,nleg
              tmp = Bijkl(ileg,i,j,k)
              Bijkl(ileg,i,j,k) = Bijkl(jleg,i,j,k)
              Bijkl(jleg,i,j,k) = tmp
            end do
          end do
        end do
        do i=1,nleg
          do j=1,nleg
            do k=1,nleg
              tmp = Bijkl(i,ileg,j,k)
              Bijkl(i,ileg,j,k) = Bijkl(i,jleg,j,k)
              Bijkl(i,jleg,j,k) = tmp
            end do
          end do
        end do
        do i=1,nleg
          do j=1,nleg
            do k=1,nleg
              tmp = Bijkl(i,j,ileg,k)
              Bijkl(i,j,ileg,k) = Bijkl(i,j,jleg,k)
              Bijkl(i,j,jleg,k) = tmp
            end do
          end do
        end do
        do i=1,nleg
          do j=1,nleg
            do k=1,nleg
              tmp = Bijkl(i,j,k,ileg)
              Bijkl(i,j,k,ileg) = Bijkl(i,j,k,jleg)
              Bijkl(i,j,k,jleg) = tmp
            end do
          end do
        end do
      end if
    end do
  end do
!
end subroutine SwapColorijkl_mom_mom
!
subroutine SwapColorijkl_part_part(nleg,pin,pout,Bijkl)
use particles
implicit none
!
  integer , intent(in) :: nleg
  type(particle) , dimension(:) , intent(in) :: pin,pout
  real(kind(1d0)) , dimension(:,:,:,:) , intent(inout) :: Bijkl
!
  integer :: ileg
  type(mom) , dimension(nleg) :: pmomin,pmomout
!
!
  do ileg=1,nleg
    pmomin(ileg)  = pin(ileg)%p
    pmomout(ileg) = pout(ileg)%p
  end do
!
  call SwapColorijkl_mom_mom(nleg,pmomin,pmomout,Bijkl)
!
end subroutine SwapColorijkl_part_part
!
subroutine SwapColormunuij_mom_mom(nleg,pin,pout,Bmunuij)
use momenta
implicit none
!
  integer , intent(in) :: nleg
  type(mom) , dimension(:) , intent(in) :: pin,pout
  real(kind(1d0)) , dimension(0:,0:,:,:) , intent(inout) :: Bmunuij
!
  integer :: ileg,jleg,mu,nu,i,j
  real(kind(1d0)) :: tmp
!
!
  do ileg=1,nleg-1
    do jleg=ileg+1,nleg
      if ((pin(ileg)%E.eq.pout(jleg)%E).and. &
          (pin(ileg)%px.eq.pout(jleg)%px).and. &
          (pin(ileg)%py.eq.pout(jleg)%py).and. &
          (pin(ileg)%pz.eq.pout(jleg)%pz)) then
        do mu=0,3
          do nu=0,3
            do j=1,nleg
              tmp = Bmunuij(mu,nu,ileg,j)
              Bmunuij(mu,nu,ileg,j) = Bmunuij(mu,nu,jleg,j)
              Bmunuij(mu,nu,jleg,j) = tmp
            end do
          end do
        end do
        do mu=0,3
          do nu=0,3
            do i=1,nleg
              tmp = Bmunuij(mu,nu,i,ileg)
              Bmunuij(mu,nu,i,ileg) = Bmunuij(mu,nu,i,jleg)
              Bmunuij(mu,nu,i,jleg) = tmp
            end do
          end do
        end do
      end if
    end do
  end do
!
end subroutine SwapColormunuij_mom_mom
!
subroutine SwapColormunuij_part_part(nleg,pin,pout,Bmunuij)
use particles
implicit none
!
  integer , intent(in) :: nleg
  type(particle) , dimension(:) , intent(in) :: pin,pout
  real(kind(1d0)) , dimension(0:,0:,:,:) , intent(inout) :: Bmunuij
!
  integer :: ileg
  type(mom) , dimension(nleg) :: pmomin,pmomout
!
!
  do ileg=1,nleg
    pmomin(ileg)  = pin(ileg)%p
    pmomout(ileg) = pout(ileg)%p
  end do
!
  call SwapColormunuij_mom_mom(nleg,pmomin,pmomout,Bmunuij)
!
end subroutine SwapColormunuij_part_part
!
function CalcFluxFact_mom(p) result(flux)
use momenta
implicit none
!
  real(kind(1d0)) :: flux
!
  type(mom) , dimension(:) , intent(in) :: p
!
!
!
  flux = 1d0/(4d0*p(1)*p(2))
!
end function CalcFluxFact_mom
!
function CalcFluxFact_part(p) result(flux)
use particles
implicit none
!
  real(kind(1d0)) :: flux
!
  type(particle) , dimension(:) , intent(in) :: p
!
!
!
  flux = 1d0/(4d0*p(1)%p*p(2)%p)
!
end function CalcFluxFact_part
!
! Trivial implementation of the NAND logic gate:
function nand(a,b) result(ans)
implicit none
!
  logical :: ans
!
  logical , intent(in) :: a,b
!
!
  ans = (.not.a).or.(.not.b)
!
end function nand
!
function SeriesProd(series1,series2) result(ans)
implicit none
!
  real(kind(1d0)) , dimension(-4:2) :: ans
!
  real(kind(1d0)) , dimension(-4:2) , intent(in) :: series1,series2
!
  integer :: iterm,jterm,kterm
!
  ans = 0d0
  do jterm=-4,2
    do kterm=2,-4,-1
      if (((kterm + jterm).ge.-4).and.((kterm + jterm).le.2)) then
        ans(kterm+jterm) = ans(kterm+jterm) + series1(jterm)*series2(kterm)
      end if
    end do
  end do
!
end function SeriesProd
!
subroutine CheckMom(p)
use momenta
implicit none
!
  type(mom) , dimension(:) , intent(in) :: p
!
  integer :: ipart
  real(kind(1d0)) , parameter :: eps = 1d-06
  type(mom) :: sum_p
!
!
  sum_p = 0d0
!
  do ipart=1,size(p)
    if (ipart.gt.2) sum_p = sum_p + p(ipart)
    if (ipart.lt.3) sum_p = sum_p - p(ipart)
  end do
!
  if (max(abs(sum_p%E), &
          abs(sum_p%px), &
          abs(sum_p%py), &
          abs(sum_p%pz)).gt.eps) then
    print *,"Error in momentum conservation..."
    call PrintMom(p)
    print *,"sum_p: "
    print *,"E,px,py,pz: ",sum_p%E,sum_p%px,sum_p%py,sum_p%pz
    stop "CheckMom"
  end if
!
end subroutine CheckMom
!
! This routine takes an array of particles and checks upon momentum
! and flavor conservation:
subroutine CheckParts(p)
use particles
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
!
  integer :: sum_flv
  integer :: ipart
!
!
  sum_flv = 0
  do ipart=1,size(p)
    if (ipart.gt.2) sum_flv = sum_flv + p(ipart)%flv
    if (ipart.lt.3) sum_flv = sum_flv - p(ipart)%flv
  end do
!
  if (sum_flv.ne.0) then
    print *,"Flavour conservation is violated..."
    call PrintParts(p)
    stop "CheckParts"
  end if
!
  call CheckMom(p(:)%p)
!
end subroutine CheckParts
!
! This routine can be used to obtain the SME from the color-correlated
! one:
subroutine CastSMEijToSME(p,SMEij,SME)
use particles
use QCDparams
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
  real(kind(1d0)) , dimension(:,:) , intent(in) :: SMEij
  real(kind(1d0)) , intent(out) :: SME
!
  integer :: ipart,npart
!
!
  npart = size(p)
! Looking for the first parton, if found divide by the quadratic 
! Casimir to obtain the SME:
  do ipart=1,npart
    if (abs(p(ipart)%flv).le.6) then
      if (p(ipart)%flv.ne.0) SME = SMEij(ipart,ipart)/qcd_cf
      if (p(ipart)%flv.eq.0) SME = SMEij(ipart,ipart)/qcd_ca
      return
    end if
  end do
!
  print *,"Error in CastSMEijToSME..."
  print *,"Unable to find a parton..."
  call PrintParts(p)
  stop "CastSMEijToSME"
!
end subroutine CastSMEijToSME
!
! This routine can be used to obtain the SME from the spin-correlated
! one:
subroutine CastSMEmunuToSME(p,SMEmunu,SME)
use particles
use QCDparams
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
  real(kind(1d0)) , dimension(0:3,0:3) , intent(in) :: SMEmunu
  real(kind(1d0)) , intent(out) :: SME
!
!
!
! The SME can be obtained by contracting the spin-correlated SME
! with the metric tensor times minus one:
  SME = SMEmunu(1,1) + SMEmunu(2,2) + SMEmunu(3,3) - SMEmunu(0,0)
!
end subroutine CastSMEmunuToSME
!
! This routine can be used to turn the simultaneously spin-
! and color-correlated SME into a spin-correlated one.
! If allpart is true then all partons are used to obtain the
! spin-correlated SME:
subroutine CastSMEmunuijToSMEmunu(allpart,p,SMEmunuij,SMEmunu)
use particles
use QCDparams
implicit none
!
  logical , intent(in) :: allpart
  type(particle) , dimension(:) , intent(in) :: p
  real(kind(1d0)) , dimension(0:,0:,:,:) , intent(in) :: SMEmunuij
  real(kind(1d0)) , dimension(0:3,0:3) , intent(out) :: SMEmunu
!
  integer :: ileg
  real(kind(1d0)) :: CiSum
!
!
  CiSum = 0d0
  SMEmunu = 0d0
  do ileg=1,size(p)
! Looking for partons:
    if (abs(p(ileg)%flv).le.6) then
! Gluon is found:
      if (abs(p(ileg)%flv).eq.0) CiSum = CiSum + qcd_ca
! Quark is found:
      if (abs(p(ileg)%flv).ne.0) CiSum = CiSum + qcd_cf
      SMEmunu(0:3,0:3) = SMEmunu(0:3,0:3) &
                       + SMEmunuij(0:3,0:3,ileg,ileg)
! If we only want to construct the spin-correlated SME in
! a fast way only one parton is considered:
      if (.not.allpart) exit
    end if
  end do
!
! The spin-correlated SME is scaled by a or the
! sum of quatdratic Casimirs:
  SMEmunu = SMEmunu / CiSum
!
end subroutine CastSMEmunuijToSMEmunu
!
! This routine can be used to turn the simultaneously spin-
! and color-correlated SME into a color-correlated one.
subroutine CastSMEmunuijToSMEij(p,SMEmunuij,SMEij)
use particles
use QCDparams
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
  real(kind(1d0)) , dimension(0:,0:,:,:) , intent(in) :: SMEmunuij
  real(kind(1d0)) , dimension(:,:) , intent(out) :: SMEij
!
!
!
  SMEij(:,:) = SMEmunuij(1,1,:,:) + SMEmunuij(2,2,:,:) &
             + SMEmunuij(3,3,:,:) - SMEmunuij(0,0,:,:)
!
end subroutine CastSMEmunuijToSMEij
!
! This routine can be used to turn the double color-correlated SME 
! into a color-correlated one.
subroutine CastSMEijklToSMEij(allpart,p,SMEijkl,SMEij)
use particles
use QCDparams
implicit none
!
  logical , intent(in) :: allpart
  type(particle) , dimension(:) , intent(in) :: p
  real(kind(1d0)) , dimension(:,:,:,:) , intent(in) :: SMEijkl
  real(kind(1d0)) , dimension(:,:) , intent(out) :: SMEij
!
  integer :: ileg,jleg,nleg
  real(kind(1d0)) :: CiSum
!
!
  SMEij = 0d0
  CiSum = 0d0
  nleg = size(p)
  do ileg=1,nleg
    if (abs(p(ileg)%flv).le.6) then
      if (p(ileg)%flv.eq.0) CiSum = CiSum + qcd_ca
      if (p(ileg)%flv.ne.0) CiSum = CiSum + qcd_cf
      SMEij(:,:) = SMEij(:,:) + SMEijkl(ileg,ileg,:,:)
      if (.not.allpart) exit
      do jleg=1,nleg
        if (ileg.eq.jleg) cycle
        SMEij(:,:) = SMEij(:,:) - SMEijkl(ileg,jleg,:,:)
      end do
    end if
  end do
  if (allpart) SMEij = SMEij / (4*CiSum)
  if (.not.allpart) SMEij = SMEij / (2*CiSum)
!
end subroutine CastSMEijklToSMEij
!
! This routine can be used to obtain the spin-correlated SME 
! from the double spin-correlated one:
subroutine CastSMEalbemunuToSMEmunu(p,SMEalbemunu,SMEmunu1,SMEmunu2)
use particles
use QCDparams
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
  real(kind(1d0)) , dimension(0:3,0:3,0:3,0:3) , intent(in) :: SMEalbemunu
  real(kind(1d0)) , dimension(0:3,0:3) , intent(out) :: SMEmunu1,SMEmunu2
!
!
!
! The spin correlated SME can be obtained by contracting
! certain indices:
  SMEmunu1 = SMEalbemunu(:,:,1,1) + SMEalbemunu(:,:,2,2) &
           + SMEalbemunu(:,:,3,3) - SMEalbemunu(:,:,0,0)
  SMEmunu2 = SMEalbemunu(1,1,:,:) + SMEalbemunu(2,2,:,:) &
           + SMEalbemunu(3,3,:,:) - SMEalbemunu(0,0,:,:)
!
end subroutine CastSMEalbemunuToSMEmunu
!
end module misc
!
module alphamax
implicit none
!
  real(kind(1d0)) :: alpha0,y0,alphamaxval
!
contains
!
subroutine init_alphamax
use input
implicit none
!
!
!
!
! Be default alpha0 and y0 coincide in value and made equal to
! alphamax:
  alphamaxval = 1d0
  if (nnloinput("#alphamax").gt.0d0) then
    alphamaxval = nnloinput("alphamax")
  end if
! We define alpha0 and y0 through alphamax, though we keep the
! possibility to give different values for them via the input card:
  alpha0 = alphamaxval
  y0     = alphamaxval
  if (nnloinput("#alpha0").gt.0d0) then
    alpha0 = nnloinput("alpha0")
  end if
  if (nnloinput("#y0").gt.0d0) then
    y0 = nnloinput("y0")
  end if
  write(*,*) "*********************************************************"
  write(*,*) "*                                                       *"
  write(*,"(1x,a,E11.4,a)") "*                 alpha0 is: ",alpha0, &
          "                 *"
  write(*,"(1x,a,E11.4,a)") "*                 y0 is:     ",y0, &
          "                 *"
  write(*,*) "*                                                       *"
  write(*,*) "*********************************************************"
!
end subroutine init_alphamax
!
end module alphamax
