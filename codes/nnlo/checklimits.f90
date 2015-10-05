module FKSutils
implicit none
!
  integer , parameter :: nexp = 5
  real(kind(1d0)) , parameter :: tiny_xi = 1d-6
  real(kind(1d0)) , parameter :: tiny_y  = 1d-6
!
contains
!
! Note that this condition is only valid if the emitted parton
! is massless:
function CalcximaxFSR(p,emit,rad)
use momenta
implicit none
!
  real(kind(1d0)) :: CalcximaxFSR
!
  type(mom) , dimension(:) , intent(in) :: p
  integer , intent(in) :: emit,rad
!
  integer :: ipart
  real(kind(1d0)) :: q2,Mrec2
  type(mom) :: q,krec
!
  q = 0d0
  krec = 0d0
  do ipart=3,size(p)
    q = q + p(ipart)
    if ((ipart.ne.emit).and.(emit.le.size(p))) krec = krec + p(ipart)
    if ((ipart.ne.rad).and.(emit.gt.size(p))) krec = krec + p(ipart)
  end do
!
  q2 = q*q
  Mrec2 = krec*krec
!  print *,"q2,Mrec2: ",q2,Mrec2
!
  CalcximaxFSR = 1d0 - Mrec2/q2
!
end function CalcximaxFSR
!
subroutine MapFSR(pbar,xitilde,y,azi,emit,rad,flv,p)
use math
use momenta
use particles
implicit none
!
  type(mom) , dimension(:) , intent(in) :: pbar
  real(kind(1d0)) , intent(in) :: xitilde,azi
  real(kind(1d0)) , intent(out) :: y
  integer , intent(in) :: emit,rad
  integer , dimension(:) , intent(in) :: flv
  type(particle) , dimension(:) , intent(out) :: p
!
  integer :: ipart,jpart
  integer :: n
  real(kind(1d0)) :: ximax,xi
  real(kind(1d0)) :: q0,q2,ki0,kr0,Mrec2,kabs,krecabs,beta
  real(kind(1d0)) :: cospsii,cospsir
  real(kind(1d0)) :: sinpsii,sinpsir
  real(kind(1d0)) , dimension(3) :: dir
  type(mom) :: k,ki,kr,krec,q
!
  type(mom) :: ptot
!
!
! xi is sampled between 0 and 1 which is non-physical,
! we have to map it into the physical region:
  ximax = CalcximaxFSR(pbar,emit,rad)
! The next line should be commented out if the MakeScatterPlots and
! the one after it should be uncommented:
  xi = ximax*(xitilde*(1d0 - 2d0*tiny_xi) + tiny_xi)
! The next line should be commented out if the MakeScatterPlots
! routine is used...
!  y  = 1.5d0*(y - y**3/3d0)*(1d0 - tiny_y)
  y  = y*(1d0 - tiny_y)
!
  n = size(pbar)
!  print *,"New point: "
!  print *,"emit, rad: ",emit,rad
!  print *,"n: ",n
!  print *,"pbar: "
!  call PrintMom(pbar)
  do ipart=1,n+1
    p(ipart)%p = 0d0
  end do
! We copy momenta from pbar to p omitting the affected ones: i and r:
  ipart = 1
  jpart = 1
  do while (.true.)
    if (ipart.gt.n) exit
! This is commented out because we decided to include the emitter:
! If we hit the position of the emitter we have to increase both
! counters because this momentum will be left out:
!    if (jpart.eq.emit) then
!      ipart = ipart + 1
!      jpart = jpart + 1
!      cycle
!    end if
! If we hit the position of the radiated one we increment only
! jpart since the momentum sitting at position ipart in pbar is not
! omittable 
    if (jpart.eq.rad) then
      jpart = jpart + 1
      cycle
    end if
    p(jpart)%p = pbar(ipart)
    ipart = ipart + 1
    jpart = jpart + 1
  end do
!  call PrintMom(p)
!  print *,"n: ",n
! We calculate q which is the sum of all final state momenta:
! and the recoiling momentum by leaving out the momentum sitting
! at the emit-th position:
  q = 0d0
  krec = 0d0
  do ipart=3,n+1
    q = q + p(ipart)%p
! This is needed to correctly sum up momenta to obtain the
! recoiling one if the emitter is the n+1st momentum:
    if (ipart.ne.emit) krec = krec + p(ipart)%p
  end do
  q0 = q%E
  q2 = q*q
! We determine the energies of the splitting pair:
! kr0 = k_r^0 ; ki0 = k_i^0:
  Mrec2 = krec*krec
  kr0 = xi*q0/2d0
  ki0 = (q2 - Mrec2 - 2d0*q0*kr0)/(2d0*(q0 - kr0*(1d0 - y)))
! We determine k which has the same direction as pbar(emit), but
! with a different normalization:
! Ek = k_i^0 + k_r^0
  kabs = sqrt(ki0**2 + kr0**2 + 2d0*ki0*kr0*y)
! k = k_i + k_r is defined to be parallel with \bar{k}_{ir}:
  k = p(emit)%p
  k%E = ki0 + kr0
  k%px = kabs*k%px/p(emit)%p%E
  k%py = kabs*k%py/p(emit)%p%E
  k%pz = kabs*k%pz/p(emit)%p%E
! We insert the splitting pair by simply inserting k, but
! rescaling accordingly, with this ki || k and kr || k,
! additional rotations will be needed:
  p(emit)%p = ki0*p(emit)%p/p(emit)%p%E
  p(rad)%p  = kr0*p(emit)%p/p(emit)%p%E
! We calculate the cosine of the angle enclosed by k and k_i:
  cospsii = (kabs**2 + ki0**2 - kr0**2)/(2d0*kabs*ki0)
  sinpsii = sqrt(abs((1d0 - cospsii)*(1d0 + cospsii)))
! We calculate the cosine of the angle enclosed by k and k_r:
  cospsir = (kabs**2 + kr0**2 - ki0**2)/(2d0*kabs*kr0)
  sinpsir = sqrt(abs((1d0 - cospsir)*(1d0 + cospsir)))
! we have to rotate both k_i and k_r:
  dir(3) = 0d0
  dir(1) = k%py
  dir(2) = -k%px
  dir = dir / sqrt(k%px**2 + k%py**2)
! We rotate k_i and k_r in opposite direction:
  call rotate4vec(dir,sinpsii,cospsii,p(emit)%p)
  call rotate4vec(dir,-sinpsir,cospsir,p(rad)%p)
! We hav to rotate the k_i - k_r system around k by the azimuth:
  dir(1) = k%px
  dir(2) = k%py
  dir(3) = k%pz
  dir = dir / sqrt(k%px**2 + k%py**2 + k%pz**2)
  call rotate4vec(dir,sin(azi),cos(azi),p(emit)%p)
  call rotate4vec(dir,sin(azi),cos(azi),p(rad)%p)
! The true k_{rec} is - k:
  krec%E =  q0 - k%E
  krec%px = -k%px
  krec%py = -k%py
  krec%pz = -k%pz
!  krecabs = sqrt(krec%px**2 + krec%py**2 + krec%pz**2)
! By construction:
  krecabs = kabs
! We have to boost the remainder momenta:
  beta = (q2 - (krec%E + krecabs)**2) &
       / (q2 + (krec%E + krecabs)**2)
! The boost should be performed on all final state
! momenta except for k_i and k_r residing at position emit and rad:
  do ipart=3,n+1
    if ((ipart.eq.rad).or.(ipart.eq.emit)) cycle
    call boost4vec(beta,dir,p(ipart)%p,p(ipart)%p)
  end do
!
! Finally, we include flavor information:
  do ipart=1,n+1
    p(ipart)%flv = flv(ipart)
  end do
!
end subroutine MapFSR
!
end module FKSutils
!
!
! This routine tests the various limits of the matrix elements:
subroutine TestLimits
use process
use flags
use subprocesses
use regions
implicit none
!
!
  integer , parameter :: iun = 99
  logical :: cuts_tmp,cutfunc_tmp,analysis_tmp
  character (len=128) :: fname
!
  interface 
    subroutine CheckNLOlimits(iun,nleg,p,ptilde,Bij, &
                              numflv,flv_arr,subtype, &
                              Cir,Sr,CSir, &
                              CalcSMEB,CalcSMEBmunu,CalcSMEBij,CalcSMER)
    use particles
    use observables
    use regions
    implicit none
!
      integer , intent(in) :: iun,nleg,numflv
      type(particle) , dimension(:) , intent(inout) :: p,ptilde
      real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
      integer , dimension(:,:) , allocatable , intent(in) :: flv_arr
      character (len=2) , intent(in) :: subtype
      type(subterms) , dimension(:) , allocatable , intent(in) :: Cir,Sr,CSir
!
      interface 
        subroutine CalcSMEB(p,smeB)
        use particles
        implicit none
!
          type(particle) , dimension(:) , intent(in) :: p
          real(kind(1d0)) , intent(out) :: smeB
!
        end subroutine CalcSMEB
!
        subroutine CalcSMEBmunu(ileg,p,Bmunu)
        use particles
        implicit none
!
          integer , intent(in) :: ileg
          type(particle) , dimension(:) , intent(in) :: p
          real(kind(1d0)) , dimension(0:3,0:3) , intent(out) :: Bmunu
!
        end subroutine CalcSMEBmunu
!
        subroutine CalcSMEBij(p,Bij)
        use particles
        implicit none
!
          type(particle) , dimension(:) , intent(in) :: p
          real(kind(1d0)) , dimension(:,:) , intent(out) :: Bij
!
        end subroutine CalcSMEBij
!
        subroutine CalcSMER(p,smeR)
        use particles
        implicit none
!
          type(particle) , dimension(:) , intent(in) :: p
          real(kind(1d0)) , intent(out) :: smeR
!
        end subroutine CalcSMER
      end interface
!
    end subroutine CheckNLOlimits
!
    subroutine CalcB(p,smeB)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: smeB
!
    end subroutine CalcB
!
    subroutine CalcR(p,smeR)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: smeR
!
    end subroutine CalcR
!
    subroutine CalcRR(p,smeRR)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: smeRR
!
    end subroutine CalcRR
!
    subroutine CalcBmunu(ileg,p,Bmunu)
    use particles
    implicit none
!
      integer , intent(in) :: ileg
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(0:3,0:3) , intent(out) :: Bmunu
!
    end subroutine CalcBmunu
!
    subroutine CalcRmunu(ileg,p,Rmunu)
    use particles
    implicit none
!
      integer , intent(in) :: ileg
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(0:3,0:3) , intent(out) :: Rmunu
!
    end subroutine CalcRmunu
!
    subroutine CalcBij(p,Bij)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Bij
!
    end subroutine CalcBij
!
    subroutine CalcRij(p,Rij)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Rij
!
    end subroutine CalcRij
!
    subroutine CheckNNLOlimits(iun,numflv,flv_arr,cont,subtype,terms)
    use particles
    use regions
    implicit none
!
      integer , intent(in) :: iun,numflv
      integer , dimension(:,:) , intent(in) :: flv_arr
      character (len=4) , intent(in) :: cont
      character (len=14) , intent(in) :: subtype
      type(subterms) , dimension(:) , intent(in) :: terms
!
    end subroutine CheckNNLOlimits
  end interface
!
!
! for the time duration of the limit check we turn off all
! possible cuts and histogramming:
  cuts_tmp     = flg_cuts
  cutfunc_tmp  = flg_cutfunc
  analysis_tmp = flg_analysis
  flg_cuts     = .false.
  flg_cutfunc  = .false.
  flg_analysis = .false.
!
  fname = "checklimits.txt"
  open(unit=iun,file=fname,status='unknown')
  write(iun,'(a)') &
  "********************************************************************"
  write(iun,'(a)') &
  "**************************** Limits ********************************"
  write(iun,'(a)') &
  "********************************************************************"
  write(iun,*) 
!
! Testing all the singular regions of the Real emission ME:
  if (flg_NLO_R.and.flg_NLO_R_A1) then
    write(iun,'(a)') &
  "####################################################################"
    write(iun,'(a)') &
  "############################# NLO R - A1 ###########################"
    write(iun,'(a)') &
  "####################################################################"
    write(iun,'(a)') &
  "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ R - Cir ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    call CheckNLOlimits(iun,nleg_born+1,psubR,psubB,subBij, &
                        num_flv_irr_NLO_R,flv_ch_Rkin,'C ', &
                        subterms_Cir_R,subterms_Sr_R,subterms_CSir_R, &
                        CalcB,CalcBmunu,CalcBij,CalcR)
    write(iun,'(a)') &
  "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ R - Sr ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    call CheckNLOlimits(iun,nleg_born+1,psubR,psubB,subBij, &
                        num_flv_irr_NLO_R,flv_ch_Rkin,'S ', &
                        subterms_Cir_R,subterms_Sr_R,subterms_CSir_R, &
                        CalcB,CalcBmunu,CalcBij,CalcR)
    write(iun,'(a)') &
  "~~~~~~~~~~~~~~~~~~~~~~~~~~~~ R - CirSr ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    call CheckNLOlimits(iun,nleg_born+1,psubR,psubB,subBij, &
                        num_flv_irr_NLO_R,flv_ch_Rkin,'CS', &
                        subterms_Cir_R,subterms_Sr_R,subterms_CSir_R, &
                        CalcB,CalcBmunu,CalcBij,CalcR)
  end if
!
! Testing all the singular regions of the Double-Real emission ME:
  if (flg_NNLO_RR.and.flg_NNLO_RR_A1) then
    write(iun,'(a)') &
  "####################################################################"
    write(iun,'(a)') &
  "############################ NNLO RR - A1 ##########################"
    write(iun,'(a)') &
  "####################################################################"
    write(iun,'(a)') &
  "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ RR - Cir ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    call CheckNLOlimits(iun,nleg_born+2,psubRR,psubR,subRij, &
                        num_flv_irr_NNLO_RR,flv_ch_RRkin,'C ', &
                        subterms_Cir_RR,subterms_Sr_RR,subterms_CSir_RR, &
                        CalcR,CalcRmunu,CalcRij,CalcRR)
    write(iun,'(a)') &
  "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ RR - Sr ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    call CheckNLOlimits(iun,nleg_born+2,psubRR,psubR,subRij, &
                        num_flv_irr_NNLO_RR,flv_ch_RRkin,'S ', &
                        subterms_Cir_RR,subterms_Sr_RR,subterms_CSir_RR, &
                        CalcR,CalcRmunu,CalcRij,CalcRR)
    write(iun,'(a)') &
  "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ RR - CirSr ~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    call CheckNLOlimits(iun,nleg_born+2,psubRR,psubR,subRij, &
                        num_flv_irr_NNLO_RR,flv_ch_RRkin,'CS', &
                        subterms_Cir_RR,subterms_Sr_RR,subterms_CSir_RR, &
                        CalcR,CalcRmunu,CalcRij,CalcRR)
  end if
!
! Testing all the singular regions of the Real-Virtual emission ME:
  if (flg_NNLO_RV.and.flg_NNLO_RV_A1) then
    write(iun,'(a)') &
  "####################################################################"
    write(iun,'(a)') &
  "############################ NNLO RV - A1 ##########################"
    write(iun,'(a)') &
  "####################################################################"
    write(iun,'(a)') &
  "~~~~~~~~~~~~~~~~~~~~~~~~~~~~ RV - Cir ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    call CheckNNLOlimits(iun,num_flv_irr_NNLO_RV,flv_ch_Rkin, &                         
                         'RV  ','C             ',subterms_Cir_R)
    write(iun,'(a)') &
  "~~~~~~~~~~~~~~~~~~~~~~~~~~~~ RV - Sr ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    call CheckNNLOlimits(iun,num_flv_irr_NNLO_RV,flv_ch_Rkin, &                         
                         'RV  ','S             ',subterms_Sr_R)
    write(iun,'(a)') &
  "~~~~~~~~~~~~~~~~~~~~~~~~~~~~ RV - CirSr ~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    call CheckNNLOlimits(iun,num_flv_irr_NNLO_RV,flv_ch_Rkin, &                         
                         'RV  ','CS            ',subterms_CSir_R)
  end if
!
! Testing all the singular regions of the A_1 subtractions for RR:
  if (flg_NNLO_R_I1.and.flg_NNLO_RRA1_A1) then
    write(iun,'(a)') &
  "####################################################################"
    write(iun,'(a)') &
  "########################## NNLO RRA1 - A1 ##########################"
    write(iun,'(a)') &
  "####################################################################"
    write(iun,'(a)') &
  "~~~~~~~~~~~~~~~~~~~~~~~~~~ RRA1 - Cir ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    call CheckNNLOlimits(iun,num_flv_irr_NLO_R,flv_ch_Rkin, &                         
                         'RRA1','C             ',subterms_Cir_R)
    write(iun,'(a)') &
  "~~~~~~~~~~~~~~~~~~~~~~~~~~~ RRA1 - Sr ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    call CheckNNLOlimits(iun,num_flv_irr_NLO_R,flv_ch_Rkin, &                         
                         'RRA1','S             ',subterms_Sr_R)
    write(iun,'(a)') &
  "~~~~~~~~~~~~~~~~~~~~~~~~~~ RRA1 - CirSr ~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    call CheckNNLOlimits(iun,num_flv_irr_NLO_R,flv_ch_Rkin, &                         
                         'RRA1','CS            ',subterms_CSir_R)
  end if
!
! Testing A2-type subtractions defined for RR:
  if (flg_NNLO_RR.and.flg_NNLO_RR_A2) then
    write(iun,'(a)') &
  "####################################################################"
    write(iun,'(a)') &
  "########################### NNLO RR - A2 ###########################"
    write(iun,'(a)') &
  "####################################################################"
    write(iun,'(a)') &
  "~~~~~~~~~~~~~~~~~~~~~~~~~~~~ RR - Cirs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    call CheckNNLOlimits(iun,num_flv_irr_NNLO_RR,flv_ch_RRkin, &                         
                         'RR  ','Cirs          ',subterms_Cirs_RR)
    write(iun,'(a)') &
  "~~~~~~~~~~~~~~~~~~~~~~~~~~~ RR - Cirjs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    call CheckNNLOlimits(iun,num_flv_irr_NNLO_RR,flv_ch_RRkin, &                         
                         'RR  ','Cirjs         ',subterms_Cirjs_RR)
    write(iun,'(a)') &
  "~~~~~~~~~~~~~~~~~~~~~~~~~~~ RR - CSirs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    call CheckNNLOlimits(iun,num_flv_irr_NNLO_RR,flv_ch_RRkin, &                         
                         'RR  ','CSirs         ',subterms_CSirs_RR)
    write(iun,'(a)') &
  "~~~~~~~~~~~~~~~~~~~~~~~~~ RR - CirsCSirs ~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    call CheckNNLOlimits(iun,num_flv_irr_NNLO_RR,flv_ch_RRkin, &                         
                         'RR  ','CirsCSirs     ',subterms_CSirs_RR)
    write(iun,'(a)') &
  "~~~~~~~~~~~~~~~~~~~~~~~~~ RR - CirjsCSirs ~~~~~~~~~~~~~~~~~~~~~~~~~~"
    call CheckNNLOlimits(iun,num_flv_irr_NNLO_RR,flv_ch_RRkin, &                         
                         'RR  ','CirjsCSirs    ',subterms_CSirs_RR)
    write(iun,'(a)') &
  "~~~~~~~~~~~~~~~~~~~~~~~~~~~~ RR - Srs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    call CheckNNLOlimits(iun,num_flv_irr_NNLO_RR,flv_ch_RRkin, &                         
                         'RR  ','Srs            ',subterms_Srs_RR)
    write(iun,'(a)') &
  "~~~~~~~~~~~~~~~~~~~~~~~~~~ RR - CirsSrs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    call CheckNNLOlimits(iun,num_flv_irr_NNLO_RR,flv_ch_RRkin, &                         
                         'RR  ','CirsSrs       ',subterms_Srs_RR)
    write(iun,'(a)') &
  "~~~~~~~~~~~~~~~~~~~~~~~~~~ RR - CSirsSrs ~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    call CheckNNLOlimits(iun,num_flv_irr_NNLO_RR,flv_ch_RRkin, &                         
                         'RR  ','CSirsSrs      ',subterms_Srs_RR)
    write(iun,'(a)') &
  "~~~~~~~~~~~~~~~~~~~~~~~~~~ RR - CirjsSrs ~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    call CheckNNLOlimits(iun,num_flv_irr_NNLO_RR,flv_ch_RRkin, &                         
                         'RR  ','CirjsSrs      ',subterms_Srs_RR)
    write(iun,'(a)') &
  "~~~~~~~~~~~~~~~~~~~~~~~~ RR - CirsCSirsSrs ~~~~~~~~~~~~~~~~~~~~~~~~~"
    call CheckNNLOlimits(iun,num_flv_irr_NNLO_RR,flv_ch_RRkin, &                         
                         'RR  ','CirsCSirsSrs  ',subterms_Srs_RR)
  end if
!
! Testing A12-type subtractions defined for RR:
  if (flg_NNLO_RR.and.flg_NNLO_RR_A12) then
    write(iun,'(a)') &
  "####################################################################"
    write(iun,'(a)') &
  "########################## NNLO RR - A12 ###########################"
    write(iun,'(a)') &
  "####################################################################"
    write(iun,'(a)') &
  "~~~~~~~~~~~~~~~~~~~~~~~~~~ RR - CktCktr ~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    call CheckNNLOlimits(iun,num_flv_irr_NNLO_RR,flv_ch_RRkin, &
                         'RR  ','CktCktr       ',subterms_Cirs_RR)
    write(iun,'(a)') &
  "~~~~~~~~~~~~~~~~~~~~~~~~~ RR - CktCirkt ~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    call CheckNNLOlimits(iun,num_flv_irr_NNLO_RR,flv_ch_RRkin, &
                         'RR  ','CktCirkt      ',subterms_Cirjs_RR)
    write(iun,'(a)') &
  "~~~~~~~~~~~~~~~~~~~~~~~~~ RR - CktCSktr ~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    call CheckNNLOlimits(iun,num_flv_irr_NNLO_RR,flv_ch_RRkin, &
                         'RR  ','CktCSktr      ',subterms_CSirs_RR)
    write(iun,'(a)') &
  "~~~~~~~~~~~~~~~~~~~~~~~ RR - CktCirktCSktr ~~~~~~~~~~~~~~~~~~~~~~~~~"
    call CheckNNLOlimits(iun,num_flv_irr_NNLO_RR,flv_ch_RRkin, &                         
                         'RR  ','CktCirktCSktr ',subterms_CSirs_RR)
    write(iun,'(a)') &
  "~~~~~~~~~~~~~~~~~~~~~~~ RR - CktCktrCSktr ~~~~~~~~~~~~~~~~~~~~~~~~~~"
    call CheckNNLOlimits(iun,num_flv_irr_NNLO_RR,flv_ch_RRkin, &
                         'RR  ','CktCktrCSktr  ',subterms_CSirs_RR)
    write(iun,'(a)') &
  "~~~~~~~~~~~~~~~~~~~~~~~~~~ RR - CktSkt ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    call CheckNNLOlimits(iun,num_flv_irr_NNLO_RR,flv_ch_RRkin, &                         
                         'RR  ','CktSkt        ',subterms_Srs_RR)
    write(iun,'(a)') &
  "~~~~~~~~~~~~~~~~~~~~~~~~ RR - CktCrktSkt ~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    call CheckNNLOlimits(iun,num_flv_irr_NNLO_RR,flv_ch_RRkin, &                         
                         'RR  ','CktCrktSkt    ',subterms_Srs_RR)
    write(iun,'(a)') &
  "~~~~~~~~~~~~~~~~~~~~~~~~ RR - StCirt     ~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    call CheckNNLOlimits(iun,num_flv_irr_NNLO_RR,flv_ch_RRkin, &                         
                         'RR  ','StCirt        ',subterms_Cirs_RR)
    write(iun,'(a)') &
  "~~~~~~~~~~~~~~~~~~~~~~~~ RR - StCSirt    ~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    call CheckNNLOlimits(iun,num_flv_irr_NNLO_RR,flv_ch_RRkin, &                         
                         'RR  ','StCSirt       ',subterms_CSirs_RR)
    write(iun,'(a)') &
  "~~~~~~~~~~~~~~~~~~~~~~~~ RR - StCirtCSirt ~~~~~~~~~~~~~~~~~~~~~~~~~~"
    call CheckNNLOlimits(iun,num_flv_irr_NNLO_RR,flv_ch_RRkin, &                         
                         'RR  ','StCirtCSirt   ',subterms_CSirs_RR)
    write(iun,'(a)') &
  "~~~~~~~~~~~~~~~~~~~~~~~~~ RR - StCirtSrt ~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    call CheckNNLOlimits(iun,num_flv_irr_NNLO_RR,flv_ch_RRkin, &                         
                         'RR  ','StCirtSrt     ',subterms_Srs_RR)
    write(iun,'(a)') &
  "~~~~~~~~~~~~~~~~~~~~~~~~ RR - StCSirtSrt ~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    call CheckNNLOlimits(iun,num_flv_irr_NNLO_RR,flv_ch_RRkin, &                         
                         'RR  ','StCSirtSrt    ',subterms_Srs_RR)
    write(iun,'(a)') &
  "~~~~~~~~~~~~~~~~~~~~~~ RR - StCirtCSirtSrt ~~~~~~~~~~~~~~~~~~~~~~~~~"
    call CheckNNLOlimits(iun,num_flv_irr_NNLO_RR,flv_ch_RRkin, &                         
                         'RR  ','StCirtCSirtSrt',subterms_Srs_RR)
    write(iun,'(a)') &
  "~~~~~~~~~~~~~~~~~~~~~~~~~~~~ RR - StSrt ~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    call CheckNNLOlimits(iun,num_flv_irr_NNLO_RR,flv_ch_RRkin, &                         
                         'RR  ','StSrt         ',subterms_Srs_RR)
    write(iun,'(a)') &
  "~~~~~~~~~~~~~~~~~~~~~~~~~~ RR - CitStCirt ~~~~~~~~~~~~~~~~~~~~~~~~~~"
    call CheckNNLOlimits(iun,num_flv_irr_NNLO_RR,flv_ch_RRkin, &                         
                         'RR  ','CitStCirt     ',subterms_Cirs_RR)
    write(iun,'(a)') &
  "~~~~~~~~~~~~~~~~~~~~~~~~~ RR - CktStCSirt ~~~~~~~~~~~~~~~~~~~~~~~~~~"
    call CheckNNLOlimits(iun,num_flv_irr_NNLO_RR,flv_ch_RRkin, &                         
                         'RR  ','CktStCSirt    ',subterms_CSirs_RR)
    write(iun,'(a)') &
  "~~~~~~~~~~~~~~~~~~~~~~~~ RR - CktStCSirtSrt ~~~~~~~~~~~~~~~~~~~~~~~~"
    call CheckNNLOlimits(iun,num_flv_irr_NNLO_RR,flv_ch_RRkin, &                         
                         'RR  ','CktStCSirtSrt ',subterms_Srs_RR)
    write(iun,'(a)') &
  "~~~~~~~~~~~~~~~~~~~~~~~~ RR - CktStCkrtSrt ~~~~~~~~~~~~~~~~~~~~~~~~"
    call CheckNNLOlimits(iun,num_flv_irr_NNLO_RR,flv_ch_RRkin, &                         
                         'RR  ','CktStCkrtSrt  ',subterms_Srs_RR)
    write(iun,'(a)') &
  "~~~~~~~~~~~~~~~~~~~~~~~~ RR - CrtStCkrtSrt ~~~~~~~~~~~~~~~~~~~~~~~~"
    call CheckNNLOlimits(iun,num_flv_irr_NNLO_RR,flv_ch_RRkin, &                         
                         'RR  ','CrtStCkrtSrt  ',subterms_Srs_RR)
    write(iun,'(a)') &
  "~~~~~~~~~~~~~~~~~~~~~~~~~~ RR - CktStSrt ~~~~~~~~~~~~~~~~~~~~~~~~~~"
    call CheckNNLOlimits(iun,num_flv_irr_NNLO_RR,flv_ch_RRkin, &                         
                         'RR  ','CktStSrt      ',subterms_Srs_RR)
    write(iun,'(a)') &
  "~~~~~~~~~~~~~~~~~~~~~~~~~~ RR - CrtStSrt ~~~~~~~~~~~~~~~~~~~~~~~~~~"
    call CheckNNLOlimits(iun,num_flv_irr_NNLO_RR,flv_ch_RRkin, &                         
                         'RR  ','CrtStSrt      ',subterms_Srs_RR)
  end if
!
!
  close(iun)
!
!  stop "TestLimits"
!
! We restore the original values:
  flg_cuts     = cuts_tmp
  flg_cutfunc  = cutfunc_tmp
  flg_analysis = analysis_tmp
!
end subroutine TestLimits
!
! This routines checks NLO-type subtraction terms to the type
! specified in subtype: C/S/CS
! level:
subroutine CheckNLOlimits(iun,nleg,p,ptilde,Bij, &
                          numflv,flv_arr,subtype, &
                          Cir,Sr,CSir, &
                          CalcSMEB,CalcSMEBmunu,CalcSMEBij,CalcSMER)
use particles
use observables
use regions
implicit none
!
  integer , intent(in) :: iun,nleg,numflv
  type(particle) , dimension(:) , intent(inout) :: p,ptilde
  real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
  integer , dimension(:,:) , allocatable , intent(in) :: flv_arr
  character (len=2) , intent(in) :: subtype
  type(subterms) , dimension(:) , allocatable , intent(in) :: Cir,Sr,CSir
!
  integer :: istat
  integer :: iproc,iterm,ipart,jpart
  type(mom) , dimension(:) , allocatable :: pborn
!
  real(kind(1d0)) , parameter :: PtMin = 10d0
  real(kind(1d0)) , parameter :: dRMin = 0.4d0
!
!
  interface 
    subroutine CirLimit(iun,pborn,p,ptilde,Bij,flv,emit,rad, &
                        Cir,Sr,CSir, &
                        CalcSMEB,CalcSMEBmunu,CalcSMEBij,CalcSMER)
    use particles
    use regions
    implicit none
!
      integer , intent(in) :: iun
      type(mom) , dimension(:) , intent(in) :: pborn
      type(particle) , dimension(:) , intent(inout) :: p,ptilde
      real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
      integer , dimension(:) , intent(in) :: flv
      integer , intent(in) :: emit,rad
      type(subterms) , intent(in) :: Cir,Sr,CSir
!
      interface 
        subroutine CalcSMEB(p,smeB)
        use particles
        implicit none
!
          type(particle) , dimension(:) , intent(in) :: p
          real(kind(1d0)) , intent(out) :: smeB
!
        end subroutine CalcSMEB
!
        subroutine CalcSMEBmunu(ileg,p,Bmunu)
        use particles
        implicit none
!
          integer , intent(in) :: ileg
          type(particle) , dimension(:) , intent(in) :: p
          real(kind(1d0)) , dimension(0:3,0:3) , intent(out) :: Bmunu
!
        end subroutine CalcSMEBmunu
!
        subroutine CalcSMEBij(p,Bij)
        use particles
        implicit none
!
          type(particle) , dimension(:) , intent(in) :: p
          real(kind(1d0)) , dimension(:,:) , intent(out) :: Bij
!
        end subroutine CalcSMEBij
!
        subroutine CalcSMER(p,smeR)
        use particles
        implicit none
!
          type(particle) , dimension(:) , intent(in) :: p
          real(kind(1d0)) , intent(out) :: smeR
!
        end subroutine CalcSMER
      end interface
!
    end subroutine CirLimit
!
    subroutine SrLimit(iun,pborn,p,ptilde,Bij,flv,rad,Cir,Sr,CSir, &
                       CalcSMEB,CalcSMEBmunu,CalcSMEBij,CalcSMER)
    use particles
    use regions
    implicit none
!
      integer , intent(in) :: iun
      type(mom) , dimension(:) , intent(in) :: pborn
      type(particle) , dimension(:) , intent(inout) :: p,ptilde
      real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
      integer , dimension(:) , intent(in) :: flv
      integer , intent(in) :: rad
      type(subterms) , intent(in) :: Cir,Sr,CSir
!
      interface 
        subroutine CalcSMEB(p,smeB)
        use particles
        implicit none
!
          type(particle) , dimension(:) , intent(in) :: p
          real(kind(1d0)) , intent(out) :: smeB
!
        end subroutine CalcSMEB
!
        subroutine CalcSMEBmunu(ileg,p,Bmunu)
        use particles
        implicit none
!
          integer , intent(in) :: ileg
          type(particle) , dimension(:) , intent(in) :: p
          real(kind(1d0)) , dimension(0:3,0:3) , intent(out) :: Bmunu
!
        end subroutine CalcSMEBmunu
!
        subroutine CalcSMEBij(p,Bij)
        use particles
        implicit none
!
          type(particle) , dimension(:) , intent(in) :: p
          real(kind(1d0)) , dimension(:,:) , intent(out) :: Bij
!
        end subroutine CalcSMEBij
!
        subroutine CalcSMER(p,smeR)
        use particles
        implicit none
!
          type(particle) , dimension(:) , intent(in) :: p
          real(kind(1d0)) , intent(out) :: smeR
!
        end subroutine CalcSMER
      end interface
!
    end subroutine SrLimit
!
    subroutine CirSrLimit(iun,pborn,p,ptilde,Bij,flv,emit,rad, &
                          Cir,Sr,CSir, &
                          CalcSMEB,CalcSMEBmunu,CalcSMEBij,CalcSMER)
    use particles
    use regions
    implicit none
!
      integer , intent(in) :: iun
      type(mom) , dimension(:) , intent(in) :: pborn
      type(particle) , dimension(:) , intent(inout) :: p,ptilde
      real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
      integer , dimension(:) , intent(in) :: flv
      integer , intent(in) :: emit,rad
      type(subterms) , intent(in) :: Cir,Sr,CSir
!
      interface 
        subroutine CalcSMEB(p,smeB)
        use particles
        implicit none
!
          type(particle) , dimension(:) , intent(in) :: p
          real(kind(1d0)) , intent(out) :: smeB
!
        end subroutine CalcSMEB
!
        subroutine CalcSMEBmunu(ileg,p,Bmunu)
        use particles
        implicit none
!
          integer , intent(in) :: ileg
          type(particle) , dimension(:) , intent(in) :: p
          real(kind(1d0)) , dimension(0:3,0:3) , intent(out) :: Bmunu
!
        end subroutine CalcSMEBmunu
!
        subroutine CalcSMEBij(p,Bij)
        use particles
        implicit none
!
          type(particle) , dimension(:) , intent(in) :: p
          real(kind(1d0)) , dimension(:,:) , intent(out) :: Bij
!
        end subroutine CalcSMEBij
!
        subroutine CalcSMER(p,smeR)
        use particles
        implicit none
!
          type(particle) , dimension(:) , intent(in) :: p
          real(kind(1d0)) , intent(out) :: smeR
!
        end subroutine CalcSMER
      end interface
!
    end subroutine CirSrLimit
!
    subroutine CalcSMEB(p,smeB)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: smeB
!
    end subroutine CalcSMEB
!
    subroutine CalcSMEBmunu(ileg,p,Bmunu)
    use particles
    implicit none
!
      integer , intent(in) :: ileg
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(0:3,0:3) , intent(out) :: Bmunu
!
    end subroutine CalcSMEBmunu
!
    subroutine CalcSMEBij(p,Bij)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Bij
!
    end subroutine CalcSMEBij
!
    subroutine CalcSMER(p,smeR)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: smeR
!
    end subroutine CalcSMER
  end interface
!
! Since we would like to use the same underlying Born phase space
! point for all flavor configs we create a separate array for it:
  allocate(pborn(nleg-1),stat=istat)
  if (istat.ne.0) then
    print *,"Pproblem with pborn allocation in CheckCir...."
    stop
  end if
!
! We generate a common Born-like kinematics:
! We strive for a phase space point where all momenta are 
! well separated:
  do while (.true.)
    call gen_ranmom(nleg-1,pborn)
    do ipart=3,nleg-1
! Checking pts:
      if (get_pt(pborn(ipart)).lt.PtMin) goto 100
      do jpart=ipart+1,nleg-1
! Checking separation:
        if (get_dR(pborn(ipart),pborn(jpart)).lt.dRMin) goto 100
      end do
    end do
    exit
100 continue
  end do
!  call PrintMom(pborn)
!
! Loop over subprocesses:
  do iproc=1,numflv
    call PrintSubProc(iun,flv_arr(:,iproc))
! C-type check
    if (subtype.eq.'C ') then
! For each subprocess several collinear regions can be present:
      do iterm=1,Cir(iproc)%numterm
        call PrintCir(iun,iterm,Cir(iproc)%term(iterm))
        call CirLimit(iun,pborn,p,ptilde,Bij,           &
                      flv_arr(:,iproc),                 &
                      Cir(iproc)%term(iterm)%rad(1)%i,  &
                      Cir(iproc)%term(iterm)%rad(2)%i,  &
                      Cir(iproc),Sr(iproc),CSir(iproc), &
                      CalcSMEB,CalcSMEBmunu,CalcSMEBij,CalcSMER)
      end do
    elseif (subtype.eq.'S ') then
      do iterm=1,Sr(iproc)%numterm
        call PrintSr(iun,iterm,Sr(iproc)%term(iterm))
        call SrLimit(iun,pborn,p,ptilde,Bij, &
                     flv_arr(:,iproc),   &
                     Sr(iproc)%term(iterm)%rad(1)%i,   &
                     Cir(iproc),Sr(iproc),CSir(iproc), &
                     CalcSMEB,CalcSMEBmunu,CalcSMEBij,CalcSMER)
      end do
    elseif (subtype.eq.'CS') then
      do iterm=1,CSir(iproc)%numterm
        call PrintCirSr(iun,iterm,CSir(iproc)%term(iterm))
        call CirSrLimit(iun,pborn,p,ptilde,Bij, &
                        flv_arr(:,iproc),   &
                        CSir(iproc)%term(iterm)%rad(1)%i,  &
                        CSir(iproc)%term(iterm)%rad(2)%i,  &
                        Cir(iproc),Sr(iproc),CSir(iproc), &
                        CalcSMEB,CalcSMEBmunu,CalcSMEBij,CalcSMER)
      end do
    end if
  end do
!
  deallocate(pborn)
!
end subroutine CheckNLOlimits
!
subroutine CheckNNLOlimits(iun,numflv,flv_arr,cont,subtype,terms)
use particles
use observables
use regions
use process
use QCDparams
use utils
implicit none
!
  integer , intent(in) :: iun,numflv
  integer , dimension(:,:) , intent(in) :: flv_arr
  character (len=4) , intent(in) :: cont
  character (len=14) , intent(in) :: subtype
  type(subterms) , dimension(:) , intent(in) :: terms
!
  integer :: i,iproc,iterm
  integer :: ipart,jpart
  integer :: istat
  integer :: nleg
  integer :: emitir,emitjs,emitirs,emitkt,emitirt,emitktr
  integer :: emit12,emit34
  integer :: radi,radj,radk,radr,rads,radt,rad,rad_pr
  integer :: rad1,rad2,rad3,rad4
  integer :: radiID,radrID,radsID,radkID,radtID
  integer :: rad1ID,rad2ID,rad3ID,rad4ID
  integer :: emitirID,emitktID,emitktrID,emitirsID
  integer :: emit12ID,emit34ID
  real(kind(1d0)) , dimension(:,:) , allocatable :: Bij
  real(kind(1d0)) , dimension(:,:) , allocatable :: Vij
  real(kind(1d0)) , dimension(:,:,:) , allocatable :: Bijk
  real(kind(1d0)) , dimension(:,:,:,:) , allocatable :: Bijkl
  real(kind(1d0)) , dimension(:,:,:,:) , allocatable :: Bmunuij
  real(kind(1d0)) , dimension(:,:) , allocatable :: Rij
  type(mom) , dimension(:) , allocatable :: pborn
  type(particle) , dimension(:) , allocatable :: p,ptilde,phat
!
  real(kind(1d0)) , parameter :: PtMin = 10d0
  real(kind(1d0)) , parameter :: dRMin = 0.4d0
!
!
  interface
    subroutine CirLimitRV(iun,pborn,p,ptilde,flv,emit,rad,Cirterm)
    use particles
    use regions
    implicit none
!
      integer , intent(in) :: iun
      type(mom) , dimension(:) , intent(in) :: pborn
      type(particle) , dimension(:) , intent(inout) :: p,ptilde
      integer , dimension(:) , intent(in) :: flv
      integer , intent(in) :: emit,rad
      type(subterms) , intent(in) :: Cirterm
!
    end subroutine CirLimitRV
!
    subroutine SrLimitRV(iun,pborn,p,ptilde,Bij,Vij,Bijk,flv,rad,Srterm)
    use particles
    use regions
    implicit none
!
      integer , intent(in) :: iun
      type(mom) , dimension(:) , intent(in) :: pborn
      type(particle) , dimension(:) , intent(inout) :: p,ptilde
      real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
      real(kind(1d0)) , dimension(:,:) , intent(inout) :: Vij
      real(kind(1d0)) , dimension(:,:,:) , intent(inout) :: Bijk
      integer , dimension(:) , intent(in) :: flv
      integer , intent(in) :: rad
      type(subterms) , intent(in) :: Srterm
!
    end subroutine SrLimitRV
!
    subroutine CirSrLimitRV(iun,pborn,p,ptilde,flv,emit,rad,CirSrterm)
    use particles
    use regions
    implicit none
!
      integer , intent(in) :: iun
      type(mom) , dimension(:) , intent(in) :: pborn
      type(particle) , dimension(:) , intent(inout) :: p,ptilde
      integer , dimension(:) , intent(in) :: flv
      integer , intent(in) :: emit,rad
      type(subterms) , intent(in) :: CirSrterm
!
    end subroutine CirSrLimitRV
!
    subroutine CirLimitRRA1(iun,pborn,p,ptilde,Bij,Bmunuij,Rij,flv,emit,rad,Cirterm)
    use particles
    use regions
    implicit none
!
      integer , intent(in) :: iun
      type(mom) , dimension(:) , intent(in) :: pborn
      type(particle) , dimension(:) , intent(inout) :: p,ptilde
      real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij,Rij
      real(kind(1d0)) , dimension(0:,0:,:,:) , intent(inout) :: Bmunuij
      integer , dimension(:) , intent(in) :: flv
      integer , intent(in) :: emit,rad
      type(subterms) , intent(in) :: Cirterm
!
    end subroutine CirLimitRRA1
!
    subroutine SrLimitRRA1(iun,pborn,p,ptilde,Bij,Bijkl,Rij,flv,rad,Srterm)
    use particles
    use regions
    implicit none
!
      integer , intent(in) :: iun
      type(mom) , dimension(:) , intent(in) :: pborn
      type(particle) , dimension(:) , intent(inout) :: p,ptilde
      real(kind(1d0)) , dimension(:,:) , allocatable , intent(inout) :: Bij,Rij
      real(kind(1d0)) , dimension(:,:,:,:) , allocatable , intent(inout) :: Bijkl
      integer , dimension(:) , intent(in) :: flv
      integer , intent(in) :: rad
      type(subterms) , intent(in) :: Srterm
!
    end subroutine SrLimitRRA1
!
    subroutine CirSrLimitRRA1(iun,pborn,p,ptilde,Bij,Rij,flv,emit,rad,CirSrterm)
    use particles
    use regions
    implicit none
!
      integer , intent(in) :: iun
      type(mom) , dimension(:) , intent(in) :: pborn
      type(particle) , dimension(:) , intent(inout) :: p,ptilde
      real(kind(1d0)) , dimension(:,:) , allocatable , intent(inout) :: Bij,Rij
      integer , dimension(:) , intent(in) :: flv
      integer , intent(in) :: emit,rad
      type(subterms) , intent(in) :: CirSrterm
!
    end subroutine CirSrLimitRRA1
!
    subroutine CirsLimit(iun,pborn,p,phat,ptilde,flv, &
                         emit,radi,radr,rads,Cirsterm)
    use particles
    use regions
    implicit none
!
      integer , intent(in) :: iun
      type(mom) , dimension(:) , intent(in) :: pborn
      type(particle) , dimension(:) , intent(inout) :: p,ptilde,phat
      integer , dimension(:) , intent(in) :: flv
      integer , intent(in) :: emit,radi,radr,rads
      type(subterms) , intent(in) :: Cirsterm
!
    end subroutine CirsLimit
!
    subroutine CirjsLimit(iun,pborn,p,phat,ptilde,flv, &
                          emiti,radi,radr,emitj,radj,rads, &
                          Cirjsterm)
    use particles
    use regions
    implicit none
!
      integer , intent(in) :: iun
      type(mom) , dimension(:) , intent(in) :: pborn
      type(particle) , dimension(:) , intent(inout) :: p,ptilde,phat
      integer , dimension(:) , intent(in) :: flv
      integer , intent(in) :: emiti,emitj,radi,radr,radj,rads
      type(subterms) , intent(in) :: Cirjsterm
!
    end subroutine CirjsLimit
!
    subroutine CSirsLimit(iun,pborn,p,phat,ptilde,Bij,Bmunuij, &
                          flv,emiti,radi,radr,rads,CSirsterm)
    use particles
    use regions
    implicit none
!
      integer , intent(in) :: iun
      type(mom) , dimension(:) , intent(in) :: pborn
      type(particle) , dimension(:) , intent(inout) :: p,ptilde,phat
      real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
      real(kind(1d0)) , dimension(0:,0:,:,:) , intent(inout) :: Bmunuij
      integer , dimension(:) , intent(in) :: flv
      integer , intent(in) :: emiti,radi,radr,rads
      type(subterms) , intent(in) :: CSirsterm
!
    end subroutine CSirsLimit
!
    subroutine CirsCSirsLimit(iun,pborn,p,phat,ptilde,  &
                              Bij,Bmunuij,              &
                              flv,emiti,radi,radr,rads, &
                              CSirsterm)
    use particles
    use regions
    implicit none
!
      integer , intent(in) :: iun
      type(mom) , dimension(:) , intent(in) :: pborn
      type(particle) , dimension(:) , intent(inout) :: p,ptilde,phat
      real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
      real(kind(1d0)) , dimension(0:,0:,:,:) , intent(inout) :: Bmunuij
      integer , dimension(:) , intent(in) :: flv
      integer , intent(in) :: emiti,radi,radr,rads
      type(subterms) , intent(in) :: CSirsterm
!
    end subroutine CirsCSirsLimit
!
    subroutine CirjsCSirsLimit(iun,pborn,p,phat,ptilde,  &
                               Bij,Bmunuij,              &
                               flv,emitir,radi,radr,     &
                               emitjs,radj,rads,         &
                               CSirsterm)
    use particles
    use regions
    implicit none
!
      integer , intent(in) :: iun
      type(mom) , dimension(:) , intent(in) :: pborn
      type(particle) , dimension(:) , intent(inout) :: p,ptilde,phat
      real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
      real(kind(1d0)) , dimension(0:,0:,:,:) , intent(inout) :: Bmunuij
      integer , dimension(:) , intent(in) :: flv
      integer , intent(in) :: emitir,emitjs,radi,radr,radj,rads
      type(subterms) , intent(in) :: CSirsterm
!
    end subroutine CirjsCSirsLimit
!
    subroutine SrsLimit(iun,pborn,p,phat,ptilde,Bij,Bijkl, &
                        flv,radr,rads,Srsterm)
    use particles
    use regions
    implicit none
!
      integer , intent(in) :: iun
      type(mom) , dimension(:) , intent(in) :: pborn
      type(particle) , dimension(:) , intent(inout) :: p,ptilde,phat
      real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
      real(kind(1d0)) , dimension(:,:,:,:) , intent(inout) :: Bijkl
      integer , dimension(:) , intent(in) :: flv
      integer , intent(in) :: radr,rads
      type(subterms) , intent(in) :: Srsterm
!
    end subroutine SrsLimit
!
    subroutine CirsSrsLimit(iun,pborn,p,phat,ptilde,Bij,Bijkl, &
                            flv,emitirs,radi,radr,rads,        &
                            Srsterm)
    use particles
    use regions
    implicit none
!
      integer , intent(in) :: iun
      type(mom) , dimension(:) , intent(in) :: pborn
      type(particle) , dimension(:) , intent(inout) :: p,ptilde,phat
      real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
      real(kind(1d0)) , dimension(:,:,:,:) , intent(inout) :: Bijkl
      integer , dimension(:) , intent(in) :: flv
      integer , intent(in) :: emitirs,radi,radr,rads
      type(subterms) , intent(in) :: Srsterm
!
    end subroutine CirsSrsLimit
!
    subroutine CSirsSrsLimit(iun,pborn,p,phat,ptilde,Bij,Bijkl, &
                             flv,emitir,radi,radr,rads,         &
                             Srsterm)
    use particles
    use regions
    implicit none
!
      integer , intent(in) :: iun
      type(mom) , dimension(:) , intent(in) :: pborn
      type(particle) , dimension(:) , intent(inout) :: p,ptilde,phat
      real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
      real(kind(1d0)) , dimension(:,:,:,:) , intent(inout) :: Bijkl
      integer , dimension(:) , intent(in) :: flv
      integer , intent(in) :: emitir,radi,radr,rads
      type(subterms) , intent(in) :: Srsterm
!
    end subroutine CSirsSrsLimit
!
    subroutine CirjsSrsLimit(iun,pborn,p,phat,ptilde,Bij,Bijkl, &
                             flv,                               &
                             emitir,radi,radr,                  &
                             emitjs,radj,rads,                  &
                             Srsterm)
    use particles
    use regions
    implicit none
!
      integer , intent(in) :: iun
      type(mom) , dimension(:) , intent(in) :: pborn
      type(particle) , dimension(:) , intent(inout) :: p,ptilde,phat
      real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
      real(kind(1d0)) , dimension(:,:,:,:) , intent(inout) :: Bijkl
      integer , dimension(:) , intent(in) :: flv
      integer , intent(in) :: emitir,emitjs,radi,radr,radj,rads
      type(subterms) , intent(in) :: Srsterm
!
    end subroutine CirjsSrsLimit
!
    subroutine CirsCSirsSrsLimit(iun,pborn,p,phat,ptilde,Bij, &
                                 flv,emitirs,radi,radr,rads,  &
                                 Srsterm)
    use particles
    use regions
    implicit none
!
      integer , intent(in) :: iun
      type(mom) , dimension(:) , intent(in) :: pborn
      type(particle) , dimension(:) , intent(inout) :: p,ptilde,phat
      real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
      integer , dimension(:) , intent(in) :: flv
      integer , intent(in) :: emitirs,radi,radr,rads
      type(subterms) , intent(in) :: Srsterm
!
    end subroutine CirsCSirsSrsLimit
!
! A12 related interface section starts here:
!
    subroutine CktCktrLimit(iun,pborn,p,phat,ptilde,    &
                            flv,emitktr,radk,radt,radr, &
                            Cirterm,Cirsterm)
    use particles
    use regions
    implicit none
!
      integer , intent(in) :: iun
      type(mom) , dimension(:) , intent(in) :: pborn
      type(particle) , dimension(:) , intent(inout) :: p,ptilde,phat
      integer , dimension(:) , intent(in) :: flv
      integer , intent(in) :: emitktr,radk,radt,radr
      type(subterms) , intent(in) :: Cirterm,Cirsterm
!
    end subroutine CktCktrLimit
!
    subroutine CktCirktLimit(iun,pborn,p,phat,ptilde, &
                             flv,                     &
                             emitir,radi,radr,        &
                             emitkt,radk,radt,        &
                             Cirterm,Cirktterm)
    use particles
    use regions
    implicit none
!
      integer , intent(in) :: iun
      type(mom) , dimension(:) , intent(in) :: pborn
      type(particle) , dimension(:) , intent(inout) :: p,ptilde,phat
      integer , dimension(:) , intent(in) :: flv
      integer , intent(in) :: emitir,radi,radr,emitkt,radk,radt
      type(subterms) , intent(in) :: Cirterm,Cirktterm
!
    end subroutine CktCirktLimit
!
    subroutine CktCSktrLimit(iun,pborn,p,phat,ptilde,Bij,Bmunuij,    &
                             flv,emitkt,radk,radt,radr, &
                             Cirterm,CSirsterm)
    use particles
    use regions
    implicit none
!
      integer , intent(in) :: iun
      type(mom) , dimension(:) , intent(in) :: pborn
      type(particle) , dimension(:) , intent(inout) :: p,ptilde,phat
      real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
      real(kind(1d0)) , dimension(0:,0:,:,:) , intent(inout) :: Bmunuij
      integer , dimension(:) , intent(in) :: flv
      integer , intent(in) :: emitkt,radk,radt,radr
      type(subterms) , intent(in) :: Cirterm,CSirsterm
!
    end subroutine CktCSktrLimit
!
    subroutine CktCirktCSktrLimit(iun,pborn,p,phat,ptilde, &
                                  flv,                     &
                                  emitir,radi,radr,        &
                                  emitkt,radk,radt,        &
                                  Cirterm,CSktrterm)
    use particles
    use regions
    implicit none
!
      integer , intent(in) :: iun
      type(mom) , dimension(:) , intent(in) :: pborn
      type(particle) , dimension(:) , intent(inout) :: p,ptilde,phat
      integer , dimension(:) , intent(in) :: flv
      integer , intent(in) :: emitir,radi,radr,emitkt,radk,radt
      type(subterms) , intent(in) :: Cirterm,CSktrterm
!
    end subroutine CktCirktCSktrLimit
!
    subroutine CktCktrCSktrLimit(iun,pborn,p,phat,ptilde, &
                                 flv,                     &
                                 emitktr,radk,radt,radr,  &
                                 Cirterm,CSktrterm)
    use particles
    use regions
    implicit none
!
      integer , intent(in) :: iun
      type(mom) , dimension(:) , intent(in) :: pborn
      type(particle) , dimension(:) , intent(inout) :: p,ptilde,phat
      integer , dimension(:) , intent(in) :: flv
      integer , intent(in) :: emitktr,radk,radt,radr
      type(subterms) , intent(in) :: Cirterm,CSktrterm
!
    end subroutine CktCktrCSktrLimit
!
    subroutine CktSktLimit(iun,pborn,p,phat,ptilde,Bij,Bijkl, &
                           flv,                               &
                           emitkt,radk,radt,                  &
                           Cirterm,Sktterm)
    use particles
    use regions
    implicit none
!
      integer , intent(in) :: iun
      type(mom) , dimension(:) , intent(in) :: pborn
      type(particle) , dimension(:) , intent(inout) :: p,ptilde,phat
      real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
      real(kind(1d0)) , dimension(:,:,:,:) , intent(inout) :: Bijkl
      integer , dimension(:) , intent(in) :: flv
      integer , intent(in) :: emitkt,radk,radt
      type(subterms) , intent(in) :: Cirterm,Sktterm
!
    end subroutine CktSktLimit
!
    subroutine CktCrktSktLimit(iun,pborn,p,phat,ptilde,   &
                           flv,                           &
                           emitktr,emitkt,radr,radk,radt, &
                           Cirterm,Sktterm)
    use particles
    use regions
    implicit none
!
      integer , intent(in) :: iun
      type(mom) , dimension(:) , intent(in) :: pborn
      type(particle) , dimension(:) , intent(inout) :: p,ptilde,phat
      integer , dimension(:) , intent(in) :: flv
      integer , intent(in) :: emitktr,emitkt,radr,radk,radt
      type(subterms) , intent(in) :: Cirterm,Sktterm
!
    end subroutine CktCrktSktLimit
!
    subroutine StCirtLimit(iun,pborn,p,phat,ptilde,Rij, &
                           flv,emitirt,radi,radr,radt,  &
                           Srterm,Cirsterm)
    use particles
    use regions
    implicit none
!
      integer , intent(in) :: iun
      type(mom) , dimension(:) , intent(in) :: pborn
      type(particle) , dimension(:) , intent(inout) :: p,ptilde,phat
      real(kind(1d0)) , dimension(:,:) , intent(inout) :: Rij
      integer , dimension(:) , intent(in) :: flv
      integer , intent(in) :: emitirt,radi,radr,radt
      type(subterms) , intent(in) :: Srterm,Cirsterm
!
    end subroutine StCirtLimit
!
    subroutine StCSirtLimit(iun,pborn,p,phat,ptilde,Bij,Rij,Bmunuij,    &
                            flv,emitir,radi,radr,radt, &
                            Stterm,CSirtterm)
    use particles
    use regions
    implicit none
!
      integer , intent(in) :: iun
      type(mom) , dimension(:) , intent(in) :: pborn
      type(particle) , dimension(:) , intent(inout) :: p,ptilde,phat
      real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij,Rij
      real(kind(1d0)) , dimension(0:,0:,:,:) , intent(inout) :: Bmunuij
      integer , dimension(:) , intent(in) :: flv
      integer , intent(in) :: emitir,radi,radr,radt
      type(subterms) , intent(in) :: Stterm,CSirtterm
!
    end subroutine StCSirtLimit
!
    subroutine StCirtCSirtLimit(iun,pborn,p,phat,ptilde, &
                                flv,                     &
                                emitirt,radi,radr,radt,  &
                                Srterm,CSirtterm)
    use particles
    use regions
    implicit none
!
      integer , intent(in) :: iun
      type(mom) , dimension(:) , intent(in) :: pborn
      type(particle) , dimension(:) , intent(inout) :: p,ptilde,phat
      integer , dimension(:) , intent(in) :: flv
      integer , intent(in) :: emitirt,radi,radr,radt
      type(subterms) , intent(in) :: Srterm,CSirtterm
!
    end subroutine StCirtCSirtLimit
!
    subroutine StCirtSrtLimit(iun,pborn,p,phat,ptilde, &
                              flv,                     &
                              emitirt,radi,radr,radt,  &
                              Srterm,Srtterm)
    use particles
    use regions
    implicit none
!
      integer , intent(in) :: iun
      type(mom) , dimension(:) , intent(in) :: pborn
      type(particle) , dimension(:) , intent(inout) :: p,ptilde,phat
      integer , dimension(:) , intent(in) :: flv
      integer , intent(in) :: emitirt,radi,radr,radt
      type(subterms) , intent(in) :: Srterm,Srtterm
!
    end subroutine StCirtSrtLimit
!
    subroutine StCSirtSrtLimit(iun,pborn,p,phat,ptilde,Bij, &
                               flv,                         &
                              emitir,radi,radr,radt,        &
                              Srterm,Srtterm)
    use particles
    use regions
    implicit none
!
      integer , intent(in) :: iun
      type(mom) , dimension(:) , intent(in) :: pborn
      type(particle) , dimension(:) , intent(inout) :: p,ptilde,phat
      real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
      integer , dimension(:) , intent(in) :: flv
      integer , intent(in) :: emitir,radi,radr,radt
      type(subterms) , intent(in) :: Srterm,Srtterm
!
    end subroutine StCSirtSrtLimit
!
    subroutine StCirtCSirtSrtLimit(iun,pborn,p,phat,ptilde, &
                                   flv,                     &
                                   emitirt,radi,radr,radt,  &
                                   Srterm,Srtterm)
    use particles
    use regions
    implicit none
!
      integer , intent(in) :: iun
      type(mom) , dimension(:) , intent(in) :: pborn
      type(particle) , dimension(:) , intent(inout) :: p,ptilde,phat
      integer , dimension(:) , intent(in) :: flv
      integer , intent(in) :: emitirt,radi,radr,radt
      type(subterms) , intent(in) :: Srterm,Srtterm
!
    end subroutine StCirtCSirtSrtLimit
!
    subroutine StSrtLimit(iun,pborn,p,phat,ptilde,Bij,Bijkl, &
                          flv,radr,radt,                     &
                          Srterm,Srtterm)
    use particles
    use regions
    implicit none
!
      integer , intent(in) :: iun
      type(mom) , dimension(:) , intent(in) :: pborn
      type(particle) , dimension(:) , intent(inout) :: p,ptilde,phat
      real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
      real(kind(1d0)) , dimension(:,:,:,:) , intent(inout) :: Bijkl
      integer , dimension(:) , intent(in) :: flv
      integer , intent(in) :: radr,radt
      type(subterms) , intent(in) :: Srterm,Srtterm
!
    end subroutine StSrtLimit
!
    subroutine CitStCirtLimit(iun,pborn,p,phat,ptilde, &
                              flv,radi,radr,radt,      &
                              CSirterm,Cirtterm)
    use particles
    use regions
    implicit none
!
      integer , intent(in) :: iun
      type(mom) , dimension(:) , intent(in) :: pborn
      type(particle) , dimension(:) , intent(inout) :: p,ptilde,phat
      integer , dimension(:) , intent(in) :: flv
      integer , intent(in) :: radi,radr,radt
      type(subterms) , intent(in) :: CSirterm,Cirtterm
!
    end subroutine CitStCirtLimit
!
    subroutine CktStCSirtLimit(iun,pborn,p,phat,ptilde,Bij,Bmunuij, &
                               flv,radi,radr,radk,radt,             &
                               CSirterm,CSirtterm)
    use particles
    use regions
    implicit none
!
      integer , intent(in) :: iun
      type(mom) , dimension(:) , intent(in) :: pborn
      type(particle) , dimension(:) , intent(inout) :: p,ptilde,phat
      real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
      real(kind(1d0)) , dimension(0:,0:,:,:) , intent(inout) :: Bmunuij
      integer , dimension(:) , intent(in) :: flv
      integer , intent(in) :: radi,radr,radk,radt
      type(subterms) , intent(in) :: CSirterm,CSirtterm
!
    end subroutine CktStCSirtLimit
!
    subroutine CktStCSirtSrtLimit(iun,pborn,p,phat,ptilde,Bij,    &
                                  flv,emitir,radi,radr,radk,radt, &
                                  CSirterm,Srtterm)
    use particles
    use regions
    implicit none
!
      integer , intent(in) :: iun
      type(mom) , dimension(:) , intent(in) :: pborn
      type(particle) , dimension(:) , intent(inout) :: p,ptilde,phat
      real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
      integer , dimension(:) , intent(in) :: flv
      integer , intent(in) :: emitir,radi,radr,radk,radt
      type(subterms) , intent(in) :: CSirterm,Srtterm
!
    end subroutine CktStCSirtSrtLimit
!
    subroutine CktStCkrtSrtLimit(iun,pborn,p,phat,ptilde,    &
                                 flv,emitkrt,radk,radr,radt, &
                                 CSirterm,Srtterm)
    use particles
    use regions
    implicit none
!
      integer , intent(in) :: iun
      type(mom) , dimension(:) , intent(in) :: pborn
      type(particle) , dimension(:) , intent(inout) :: p,ptilde,phat
      integer , dimension(:) , intent(in) :: flv
      integer , intent(in) :: emitkrt,radk,radr,radt
      type(subterms) , intent(in) :: CSirterm,Srtterm
!
    end subroutine CktStCkrtSrtLimit
!
    subroutine CrtStCkrtSrtLimit(iun,pborn,p,phat,ptilde,    &
                                 flv,emitkrt,radk,radr,radt, &
                                 CSirterm,Srtterm)
    use particles
    use regions
    implicit none
!
      integer , intent(in) :: iun
      type(mom) , dimension(:) , intent(in) :: pborn
      type(particle) , dimension(:) , intent(inout) :: p,ptilde,phat
      integer , dimension(:) , intent(in) :: flv
      integer , intent(in) :: emitkrt,radk,radr,radt
      type(subterms) , intent(in) :: CSirterm,Srtterm
!
    end subroutine CrtStCkrtSrtLimit
!
    subroutine CktStSrtLimit(iun,pborn,p,phat,ptilde,Bij,Bijkl, &
                                 flv,radk,radr,radt,            &
                                 CSirterm,Srtterm)
    use particles
    use regions
    implicit none
!
      integer , intent(in) :: iun
      type(mom) , dimension(:) , intent(in) :: pborn
      type(particle) , dimension(:) , intent(inout) :: p,ptilde,phat
      real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
      real(kind(1d0)) , dimension(:,:,:,:) , intent(inout) :: Bijkl
      integer , dimension(:) , intent(in) :: flv
      integer , intent(in) :: radk,radr,radt
      type(subterms) , intent(in) :: CSirterm,Srtterm
!
    end subroutine CktStSrtLimit
!
    subroutine CrtStSrtLimit(iun,pborn,p,phat,ptilde,Bij,Bijkl, &
                                 flv,radr,radt,                 &
                                 CSirterm,Srtterm)
    use particles
    use regions
    implicit none
!
      integer , intent(in) :: iun
      type(mom) , dimension(:) , intent(in) :: pborn
      type(particle) , dimension(:) , intent(inout) :: p,ptilde,phat
      real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
      real(kind(1d0)) , dimension(:,:,:,:) , intent(inout) :: Bijkl
      integer , dimension(:) , intent(in) :: flv
      integer , intent(in) :: radr,radt
      type(subterms) , intent(in) :: CSirterm,Srtterm
!
    end subroutine CrtStSrtLimit
  end interface
!
! To test the behavior of subtraction terms near in the
! singular regions an underlying Born phase space point has to be
! generated, an array is allocated to accomodate for the momenta:
  allocate(pborn(nleg_born),stat=istat)
  if (istat.ne.0) then
    print *,"Allocation went wrong for pborn in CheckNNLOlimits..."
    stop
  end if
!
  if ((cont.eq.'RV  ').or. &
      (cont.eq.'RRA1')) then
    nleg = nleg_born + 1
  elseif (cont.eq.'RR  ') then
    nleg = nleg_born + 2
  else
    print *,"Unknown contribution is given to CheckNNLOlimits..."
    print *,"cont: ",cont
    stop "CheckNNLOlimits"
  end if
!
! Allocation of further momentum arrays and an array for the
! color-correlated ME:
  allocate(p(nleg), &
           phat(nleg_born+1), &
           ptilde(nleg_born), &
           Bij(nleg_born,nleg_born), &
           Vij(nleg_born,nleg_born), &
           Bijk(nleg_born,nleg_born,nleg_born), &
           Bijkl(nleg_born,nleg_born,nleg_born,nleg_born), &
           Bmunuij(0:3,0:3,nleg_born,nleg_born), &
           Rij(nleg_born+1,nleg_born+1), &
           stat=istat)
  if (istat.ne.0) then
    print *,"Allocation went wrong in CheckNNLOlimits..."
    stop
  end if
!
! An underlying Born kinematics has to be built up this is
! done in such a fashion to have well-separated particles:
  do while (.true.)
    call gen_ranmom(nleg_born,pborn)
    do ipart=3,nleg_born
! Checking pts:
      if (get_pt(pborn(ipart)).lt.PtMin) goto 100
      do jpart=ipart+1,nleg_born
! Checking separation:
        if (get_dR(pborn(ipart),pborn(jpart)).lt.dRMin) goto 100
      end do
    end do
    exit
100 continue
  end do
!
! Loop over all subprocesses:
  do iproc=1,numflv
    call PrintSubProc(iun,flv_arr(:,iproc))
!    call PrintSubProc(flv_arr(:,iproc))
! For each subprocess several separate subtraction terms
! can be present:
!    print *,"subtype: ",subtype
!    print *,"number of subterms: ",terms(iproc)%numterm
    do iterm=1,terms(iproc)%numterm
! Having several subtraction terms, different in type:
      if (cont.eq.'RV  ') then
        if (subtype.eq.'C             ') then
          call PrintCir(iun,iterm,terms(iproc)%term(iterm))
!          call PrintCir(6,iterm,terms(iproc)%term(iterm))
          call CirLimitRV(iun,pborn,p,ptilde,                &
                          flv_arr(:,iproc),                  &
                          terms(iproc)%term(iterm)%rad(1)%i, &
                          terms(iproc)%term(iterm)%rad(2)%i, &
                          terms(iproc))
        elseif (subtype.eq.'S             ') then
          call PrintSr(iun,iterm,terms(iproc)%term(iterm))
!          call PrintSr(6,iterm,terms(iproc)%term(iterm))
          call SrLimitRV(iun,pborn,p,ptilde,Bij,Vij,Bijk,   &
                         flv_arr(:,iproc),                  &
                         terms(iproc)%term(iterm)%rad(1)%i, &
                         terms(iproc))
        elseif (subtype.eq.'CS            ') then
          call PrintCirSr(iun,iterm,terms(iproc)%term(iterm))
!          call PrintCirSr(6,iterm,terms(iproc)%term(iterm))
          call CirSrLimitRV(iun,pborn,p,ptilde,                &
                            flv_arr(:,iproc),                  &
                            terms(iproc)%term(iterm)%rad(1)%i, &
                            terms(iproc)%term(iterm)%rad(2)%i, &
                            terms(iproc))
        else
          print *,"Wrong subtraction term type is specified..."
          print *,"cont: ",cont
          print *,"subtype: ",subtype
          stop
        end if
      elseif (cont.eq.'RRA1') then
        if (subtype.eq.'C             ') then
          call PrintCir(iun,iterm,terms(iproc)%term(iterm))
!          call PrintCir(6,iterm,terms(iproc)%term(iterm))
          call CirLimitRRA1(iun,pborn,p,ptilde,Bij,Bmunuij,Rij, &
                            flv_arr(:,iproc),                   &
                            terms(iproc)%term(iterm)%rad(1)%i,  &
                            terms(iproc)%term(iterm)%rad(2)%i,  &
                            terms(iproc))
        elseif (subtype.eq.'S             ') then
          call PrintSr(iun,iterm,terms(iproc)%term(iterm))
!          call PrintSr(6,iterm,terms(iproc)%term(iterm))
          call SrLimitRRA1(iun,pborn,p,ptilde,Bij,Bijkl,Rij,  &
                           flv_arr(:,iproc),                  &
                           terms(iproc)%term(iterm)%rad(1)%i, &
                           terms(iproc))
        elseif (subtype.eq.'CS            ') then
          call PrintCirSr(iun,iterm,terms(iproc)%term(iterm))
!          call PrintCirSr(6,iterm,terms(iproc)%term(iterm))
          call CirSrLimitRRA1(iun,pborn,p,ptilde,Bij,Rij,        &
                              flv_arr(:,iproc),                  &
                              terms(iproc)%term(iterm)%rad(1)%i, &
                              terms(iproc)%term(iterm)%rad(2)%i, &
                              terms(iproc))
        else
          print *,"Wrong subtraction term type is specified..."
          print *,"cont: ",cont
          print *,"subtype: ",subtype
          stop "CheckNNLOlimits"
        end if
      elseif (cont.eq.'RR  ') then
        if (subtype.eq.'Cirs          ') then
          call PrintCirs(iun,iterm,terms(iproc)%term(iterm))
!          call PrintCirs(6,iterm,terms(iproc)%term(iterm))
          call CirsLimit(iun,pborn,p,phat,ptilde,            &
                         flv_arr(:,iproc),                   &
                         terms(iproc)%term(iterm)%emit(1)%i, &
                         terms(iproc)%term(iterm)%rad(1)%i,  &
                         terms(iproc)%term(iterm)%rad(2)%i,  &
                         terms(iproc)%term(iterm)%rad(3)%i,  &
                         terms(iproc))
        elseif (subtype.eq.'Cirjs         ') then
          call PrintCirjs(iun,iterm,terms(iproc)%term(iterm))
!          call PrintCirjs(6,iterm,terms(iproc)%term(iterm))
          call CirjsLimit(iun,pborn,p,phat,ptilde,            &
                          flv_arr(:,iproc),                   &
                          terms(iproc)%term(iterm)%emit(1)%i, &
                          terms(iproc)%term(iterm)%rad(1)%i,  &
                          terms(iproc)%term(iterm)%rad(2)%i,  &
                          terms(iproc)%term(iterm)%emit(2)%i, &
                          terms(iproc)%term(iterm)%rad(3)%i,  &
                          terms(iproc)%term(iterm)%rad(4)%i,  &
                          terms(iproc))
        elseif (subtype.eq.'CSirs         ') then
          call PrintCSirs(iun,iterm,terms(iproc)%term(iterm))
!          call PrintCSirs(6,iterm,terms(iproc)%term(iterm))
          call CSirsLimit(iun,pborn,p,phat,ptilde,Bij,Bmunuij, &
                          flv_arr(:,iproc),                    &
                          terms(iproc)%term(iterm)%emit(1)%i,  &
                          terms(iproc)%term(iterm)%rad(1)%i,   &
                          terms(iproc)%term(iterm)%rad(2)%i,   &
                          terms(iproc)%term(iterm)%rad(3)%i,   &
                          terms(iproc))
        elseif (subtype.eq.'CirsCSirs     ') then
          call PrintCirsCSirs(iun,iterm,terms(iproc)%term(iterm))
!          call PrintCirsCSirs(6,iterm,terms(iproc)%term(iterm))
          call CirsCSirsLimit(iun,pborn,p,phat,ptilde,Bij,Bmunuij, &
                              flv_arr(:,iproc),                    &
                              terms(iproc)%term(iterm)%emit(1)%i,  &
                              terms(iproc)%term(iterm)%rad(1)%i,   &
                              terms(iproc)%term(iterm)%rad(2)%i,   &
                              terms(iproc)%term(iterm)%rad(3)%i,   &
                              terms(iproc))
        elseif (subtype.eq.'CirjsCSirs    ') then
! The CirjsCSirs counterterm needs a bit different treatment provided
! by the extra j index which is not present in the CSirs term, hence
! a loop has to be defined all massless partons of the RR subprocess
! and pick all those indices for which a (js) -> j + s splitting is
! possible:
          emitir = terms(iproc)%term(iterm)%emit(1)%i
          radi   = terms(iproc)%term(iterm)%rad(1)%i
          radr   = terms(iproc)%term(iterm)%rad(2)%i
          rads   = terms(iproc)%term(iterm)%rad(3)%i
          do radj=1,nleg
! the additional radiated leg cannot coincide with other radiated 
! partons:
            if ((radj.eq.radi).or.(radj.eq.radr).or.(radj.eq.rads)) cycle
! Only massless partons count:
            if (abs(flv_arr(radj,iproc)).gt.qcd_nf) cycle
! Calculating the position of the emitter in the underlying Born
! configuration:
            emitjs = min(radj,rads)
! If the (ir) -> i + r splitting comes first the position has to be
! decreased by one unit:
            if (emitjs.gt.max(radi,radr)) emitjs = emitjs - 1
!
!            call PrintSubProc(flv_arr(:,iproc))
!
!            call PrintCirjsCSirs(6,iterm,emitjs,radj,flv_arr(radj,iproc), &
!                                 terms(iproc)%term(iterm))
            call PrintCirjsCSirs(iun,iterm,emitjs,radj,flv_arr(radj,iproc), &
                                 terms(iproc)%term(iterm))
            call CirjsCSirsLimit(iun,pborn,p,phat,ptilde,Bij,Bmunuij, &
                                 flv_arr(:,iproc),                    &
                                 terms(iproc)%term(iterm)%emit(1)%i,  &
                                 terms(iproc)%term(iterm)%rad(1)%i,   &
                                 terms(iproc)%term(iterm)%rad(2)%i,   &
                                 emitjs,radj,                         &
                                 terms(iproc)%term(iterm)%rad(3)%i,   &
                                 terms(iproc))
          end do
        elseif (subtype.eq.'Srs           ') then
          call PrintSrs(iun,iterm,terms(iproc)%term(iterm))
!          call PrintSrs(6,iterm,terms(iproc)%term(iterm))
          call SrsLimit(iun,pborn,p,phat,ptilde,Bij,Bijkl, &
                        flv_arr(:,iproc),                  &
                        terms(iproc)%term(iterm)%rad(1)%i, &
                        terms(iproc)%term(iterm)%rad(2)%i, &
                        terms(iproc))
        elseif (subtype.eq.'CirsSrs       ') then
          radr = terms(iproc)%term(iterm)%rad(1)%i
          rads = terms(iproc)%term(iterm)%rad(2)%i
! Determining the position for radi:
          do radi=1,nleg
! Cannot coincide with the soft partons:
            if ((radi.eq.radr).or.(radi.eq.rads)) cycle
! It has to be a massless parton:
            if (abs(flv_arr(radi,iproc)).gt.qcd_nf) cycle
! Determining the position for the emitter of the i-r-s triplet:
            do emitirs=3,nleg-2
              if (terms(iproc)%term(iterm)%UBorn(emitirs).eq. &
                  (flv_arr(radi,iproc) + &
                   flv_arr(radr,iproc) + &
                   flv_arr(rads,iproc))) exit
            end do
            call PrintCirsSrs(iun,iterm,emitirs,radi,flv_arr(radi,iproc), &
                              terms(iproc)%term(iterm))
!            call PrintCirsSrs(6,iterm,emitirs,radi,flv_arr(radi,iproc), &
!                              terms(iproc)%term(iterm))
            call CirsSrsLimit(iun,pborn,p,phat,ptilde,Bij,Bijkl, &
                              flv_arr(:,iproc),                  &
                              emitirs,radi,                      &
                              terms(iproc)%term(iterm)%rad(1)%i, &
                              terms(iproc)%term(iterm)%rad(2)%i, &
                              terms(iproc))
          end do
        elseif (subtype.eq.'CSirsSrs      ') then
          do rad=1,2
            radr = terms(iproc)%term(iterm)%rad(rad)%i
            rads = terms(iproc)%term(iterm)%rad(mod(rad,2)+1)%i
! This counter term is only possible if both r and s are gluons:
            if ((flv_arr(radr,iproc).ne.0).or.(flv_arr(rads,iproc).ne.0)) cycle
! Determining the position for radi:
            do radi=1,nleg
! Cannot coincide with the soft partons:
              if ((radi.eq.radr).or.(radi.eq.rads)) cycle
! It has to be a massless parton:
              if (abs(flv_arr(radi,iproc)).gt.qcd_nf) cycle
! Determining the position for the emitter of the i-r pair:
              do emitir=3,nleg-2
                if (terms(iproc)%term(iterm)%UBorn(emitir).eq. &
                    (flv_arr(radi,iproc) + &
                     flv_arr(radr,iproc))) exit
              end do
              call PrintCSirsSrs(iun,iterm,emitir,radi,flv_arr(radi,iproc), &
                                 radr,rads,terms(iproc)%term(iterm))
!              call PrintCSirsSrs(6,iterm,emitir,radi,flv_arr(radi,iproc), &
!                                 radr,rads,terms(iproc)%term(iterm))
              call CSirsSrsLimit(iun,pborn,p,phat,ptilde,Bij,Bijkl, &
                                 flv_arr(:,iproc),                  &
                                 emitir,radi,radr,rads,             &
                                 terms(iproc))
            end do
          end do
        elseif (subtype.eq.'CirjsSrs      ') then
          radr = terms(iproc)%term(iterm)%rad(1)%i
          rads = terms(iproc)%term(iterm)%rad(2)%i
! This counter term is only possible if both r and s are gluons:
          if ((flv_arr(radr,iproc).ne.0).or.(flv_arr(rads,iproc).ne.0)) cycle
! Determining the position for radi:
          do radi=1,nleg
! Cannot coincide with the soft partons:
            if ((radi.eq.radr).or.(radi.eq.rads)) cycle
! It has to be a massless parton:
            if (abs(flv_arr(radi,iproc)).gt.qcd_nf) cycle
! Determining the position for radj, note that radr and rads are
! ordered hence radi and radj can be anything:
            do radj=1,nleg
! radi and radj cannot coincide:
              if (radi.eq.radj) cycle
! radj cannot be the same as the soft legs:
              if ((radj.eq.radr).or.(radj.eq.rads)) cycle
! Must be a massless parton:
              if (abs(flv_arr(radj,iproc)).gt.qcd_nf) cycle
! Determining the position for the emitter of the i-r duplet:
              do emitir=3,nleg-2
                if (terms(iproc)%term(iterm)%UBorn(emitir).eq. &
                    (flv_arr(radi,iproc) + &
                     flv_arr(radr,iproc))) exit
              end do
              do emitjs=3,nleg-2
                if (emitjs.eq.emitir) cycle
                if (terms(iproc)%term(iterm)%UBorn(emitjs).eq. &
                    (flv_arr(radj,iproc) + &
                     flv_arr(rads,iproc))) exit
              end do
!              print *,"radi,radr,radj,rads: ",radi,radr,radj,rads
!              print *,"emitir,emitjs: ",emitir,emitjs
              call PrintCirjsSrs(iun,iterm, &
                                 emitir,radi,flv_arr(radi,iproc), &
                                 emitjs,radj,flv_arr(radj,iproc), &
                                 terms(iproc)%term(iterm))
!              call PrintCirjsSrs(6,iterm, &
!                                 emitir,radi,flv_arr(radi,iproc), &
!                                 emitjs,radj,flv_arr(radj,iproc), &
!                                 terms(iproc)%term(iterm))
              call CirjsSrsLimit(iun,pborn,p,phat,ptilde,Bij,Bijkl, &
                                 flv_arr(:,iproc),                  &
                                 emitir,radi,radr,emitjs,radj,rads,           &
                                 terms(iproc))
            end do
          end do
        elseif (subtype.eq.'CirsCSirsSrs  ') then
          do rad=1,2
            radr = terms(iproc)%term(iterm)%rad(rad)%i
            rads = terms(iproc)%term(iterm)%rad(mod(rad,2)+1)%i
! This counter term is only possible if both r and s are gluons:
            if ((flv_arr(radr,iproc).ne.0).or.(flv_arr(rads,iproc).ne.0)) cycle
! Determining the position for radi:
            do radi=1,nleg
! Cannot coincide with the soft partons:
              if ((radi.eq.radr).or.(radi.eq.rads)) cycle
! It has to be a massless parton:
              if (abs(flv_arr(radi,iproc)).gt.qcd_nf) cycle
! Determining the position for the emitter of the i-r-s triplet,
! note that s is also soft!
              do emitirs=3,nleg-2
                if (terms(iproc)%term(iterm)%UBorn(emitirs).eq. &
                    (flv_arr(radi,iproc) + &
                     flv_arr(radr,iproc) + &
                     flv_arr(rads,iproc))) exit
              end do
              call PrintCirsCSirsSrs(iun,iterm, &
                                     emitirs,radi,flv_arr(radi,iproc), &
                                     radr,rads,terms(iproc)%term(iterm))
!              call PrintCirsCSirsSrs(6,iterm, &
!                                     emitirs,radi,flv_arr(radi,iproc), &
!                                     radr,rads,terms(iproc)%term(iterm))
              call CirsCSirsSrsLimit(iun,pborn,p,phat,ptilde,Bij, &
                                     flv_arr(:,iproc),            &
                                     emitirs,radi,radr,rads,      &
                                     terms(iproc))
            end do
          end do
        elseif (subtype.eq.'CktCktr       ') then
          rad1      = terms(iproc)%term(iterm)%rad(1)%i
          rad2      = terms(iproc)%term(iterm)%rad(2)%i
          rad3      = terms(iproc)%term(iterm)%rad(3)%i
          emitktr   = terms(iproc)%term(iterm)%emit(1)%i
          rad1ID    = terms(iproc)%term(iterm)%rad(1)%ID
          rad2ID    = terms(iproc)%term(iterm)%rad(2)%ID
          rad3ID    = terms(iproc)%term(iterm)%rad(3)%ID
          emitktrID = terms(iproc)%term(iterm)%emit(1)%ID
! Cirs corresponds to one certain ordering among radiated partons.
! all possible combinations have to be generated to find
! all k-t combinations:
          if (((rad1ID+rad2ID).eq.emitktrID).or. &
              ((rad1ID+rad2ID).eq.0).or. &
              ((rad1ID+rad2ID).eq.-rad3ID)) then
            call OrderRadir(rad1,rad1ID,rad2,rad2ID, &
                            radk,radkID,radt,radtID)
! For A12 the mapping is iterated where the first one is the
! collinear pair involving the k-t pair, hence the position
! for the emitter just coincides with the position of the lowest
! lying radiated parton, since the gluon has the ID of zero
! the flavor of the emitter is just the sum of the flavors
! of the two radiated partons:
            emitkt   = min(radk,radt)
            emitktID = radkID + radtID
            radr     = rad3
            radrID   = rad3ID
            call PrintCktCktr(iun,iterm,terms(iproc)%term(iterm), &
                              radk,radkID,radt,radtID,radr,radrID)
!            call PrintCktCktr(6,iterm,terms(iproc)%term(iterm), &
!                              radk,radkID,radt,radtID,radr,radrID)
            call CktCktrLimit(iun,pborn,p,phat,ptilde, &
                              flv_arr(:,iproc),        &
                              emitktr,radk,radt,radr,  &
                              subterms_Cir_RR(iproc),  &
                              terms(iproc))
          end if
          if (((rad1ID+rad3ID).eq.emitktrID).or. &
              ((rad1ID+rad3ID).eq.0).or. &
              ((rad1ID+rad3ID).eq.-rad2ID)) then
            call OrderRadir(rad1,rad1ID,rad3,rad3ID, &
                            radk,radkID,radt,radtID)
            emitkt   = min(radk,radt)
            emitktID = radkID + radtID
            radr     = rad2
            radrID   = rad2ID
            call PrintCktCktr(iun,iterm,terms(iproc)%term(iterm), &
                              radk,radkID,radt,radtID,radr,radrID)
!            call PrintCktCktr(6,iterm,terms(iproc)%term(iterm), &
!                              radk,radkID,radt,radtID,radr,radrID)
            call CktCktrLimit(iun,pborn,p,phat,ptilde, &
                              flv_arr(:,iproc),        &
                              emitktr,radk,radt,radr,  &
                              subterms_Cir_RR(iproc),  &
                              terms(iproc))
          end if
          if (((rad2ID+rad3ID).eq.emitktrID).or. &
              ((rad2ID+rad3ID).eq.0).or. &
              ((rad2ID+rad3ID).eq.-rad1ID)) then
            call OrderRadir(rad2,rad2ID,rad3,rad3ID, &
                            radk,radkID,radt,radtID)
            emitkt   = min(radk,radt)
            emitktID = radkID + radtID
            radr     = rad1
            radrID   = rad1ID
            call PrintCktCktr(iun,iterm,terms(iproc)%term(iterm), &
                              radk,radkID,radt,radtID,radr,radrID)
!            call PrintCktCktr(6,iterm,terms(iproc)%term(iterm), &
!                              radk,radkID,radt,radtID,radr,radrID)
            call CktCktrLimit(iun,pborn,p,phat,ptilde, &
                              flv_arr(:,iproc),        &
                              emitktr,radk,radt,radr,  &
                              subterms_Cir_RR(iproc),  &
                              terms(iproc))
          end if
        elseif (subtype.eq.'CktCirkt      ') then
          rad1     = terms(iproc)%term(iterm)%rad(1)%i
          rad2     = terms(iproc)%term(iterm)%rad(2)%i
          rad3     = terms(iproc)%term(iterm)%rad(3)%i
          rad4     = terms(iproc)%term(iterm)%rad(4)%i
          emit12   = terms(iproc)%term(iterm)%emit(1)%i
          emit34   = terms(iproc)%term(iterm)%emit(2)%i
          rad1ID   = terms(iproc)%term(iterm)%rad(1)%ID
          rad2ID   = terms(iproc)%term(iterm)%rad(2)%ID
          rad3ID   = terms(iproc)%term(iterm)%rad(3)%ID
          rad4ID   = terms(iproc)%term(iterm)%rad(4)%ID
          emit12ID = terms(iproc)%term(iterm)%emit(1)%ID
          emit34ID = terms(iproc)%term(iterm)%emit(2)%ID
! In the case of the Cirjs counterterms the counterterms are defined
! such no symmetry factor is needed, hence each Cirjs term will generate
! two distinct CktCirkt term:
! Note that our counterterms are ordered hence no twiddling is needed
! with the order of the legs, in other words: we always have:
! q(~) -> q(~) g , g -> q q~ , g -> g(i) g(j) , such that i < j
          do i=1,2
            if (i.eq.1) then
! radi = rad1, radr = rad2, radk = rad3, radt = rad4:
              radi = rad1 ; radr = rad2 ; radk = rad3 ; radt = rad4
              radiID = rad1ID ; radrID = rad2ID ; radkID = rad3ID ; radtID = rad4ID
              emitir = emit12 ; emitkt = emit34
              emitirID = emit12ID ; emitktID = emit34ID
            else
! radi = rad3, radr = rad4, radk = rad1, radt = rad2:
              radi = rad3 ; radr = rad4 ; radk = rad1 ; radt = rad2
              radiID = rad3ID ; radrID = rad4ID ; radkID = rad1ID ; radtID = rad2ID
              emitir = emit34 ; emitkt = emit12
              emitirID = emit34ID ; emitktID = emit12ID
            end if
            call PrintCktCirkt(iun,iterm,terms(iproc)%term(iterm), &
                               emitir,emitirID,radi,radiID,radr,radrID, &
                               emitkt,emitktID,radk,radkID,radt,radtID)
!            call PrintCktCirkt(6,iterm,terms(iproc)%term(iterm), &
!                               emitir,emitirID,radi,radiID,radr,radrID, &
!                               emitkt,emitktID,radk,radkID,radt,radtID)
            call CktCirktLimit(iun,pborn,p,phat,ptilde, &
                               flv_arr(:,iproc),        &
                               emitir,radi,radr,        &
                               emitkt,radk,radt,        &
                               subterms_Cir_RR(iproc),  &
                               terms(iproc))
          end do
        elseif (subtype.eq.'CktCSktr      ') then
          radk   = terms(iproc)%term(iterm)%rad(1)%i
          radt   = terms(iproc)%term(iterm)%rad(2)%i
          radr   = terms(iproc)%term(iterm)%rad(3)%i
          emitkt = terms(iproc)%term(iterm)%emit(1)%i
          call PrintCktCSktr(iun,iterm,terms(iproc)%term(iterm))
!          call PrintCktCSktr(6,iterm,terms(iproc)%term(iterm))
          call CktCSktrLimit(iun,pborn,p,phat,ptilde,Bij,Bmunuij, &
                             flv_arr(:,iproc),                    &
                             emitkt,radk,radt,radr,               &
                             subterms_Cir_RR(iproc),              &
                             terms(iproc))
        elseif (subtype.eq.'CktCirktCSktr ') then
          radk   = terms(iproc)%term(iterm)%rad(1)%i
          radt   = terms(iproc)%term(iterm)%rad(2)%i
          radr   = terms(iproc)%term(iterm)%rad(3)%i
! For this counterterm one additional leg has to be found:
          do radi=1,nleg
            if (abs(flv_arr(radi,iproc)).gt.qcd_nf) cycle
            if ((radi.eq.radk).or.(radi.eq.radt).or.(radi.eq.radr)) cycle
! No additional constraint is present provided parton r can 
! only be a gluon.
! Obtaining the emitter for the i-r pair:
            emitir = min(radi,radr)
! The position for the emitter of the k-t pair can change
! hence the original position should be obtained for each
! counterterm:
            emitkt = terms(iproc)%term(iterm)%emit(1)%i
! The positions of the emitters can change depending upon the
! order of legs:
            if (emitir.gt.max(radk,radt)) emitir = emitir - 1
            if (emitkt.gt.max(radi,radr)) emitkt = emitkt - 1
            call PrintCktCirktCSktr(iun,iterm,emitir,radi, &
                                    flv_arr(radi,iproc), &
                                    terms(iproc)%term(iterm))
!            call PrintCktCirktCSktr(6,iterm,emitir,radi, &
!                                    flv_arr(radi,iproc), &
!                                    terms(iproc)%term(iterm))
            call CktCirktCSktrLimit(iun,pborn,p,phat,ptilde, &
                                    flv_arr(:,iproc),        &
                                    emitir,radi,radr,        &
                                    emitkt,radk,radt,        &
                                    subterms_Cir_RR(iproc),  &
                                    terms(iproc))
          end do
        elseif (subtype.eq.'CktCktrCSktr  ') then
          radk    = terms(iproc)%term(iterm)%rad(1)%i
          radt    = terms(iproc)%term(iterm)%rad(2)%i
          radr    = terms(iproc)%term(iterm)%rad(3)%i
          emitktr = terms(iproc)%term(iterm)%emit(1)%i
          call PrintCktCktrCSktr(iun,iterm,terms(iproc)%term(iterm))
!          call PrintCktCktrCSktr(6,iterm,terms(iproc)%term(iterm))
          call CktCktrCSktrLimit(iun,pborn,p,phat,ptilde, &
                                 flv_arr(:,iproc),        &
                                 emitktr,radk,radt,radr,  &
                                 subterms_Cir_RR(iproc),  &
                                 terms(iproc))
        elseif (subtype.eq.'CktSkt        ') then
          radk   = terms(iproc)%term(iterm)%rad(1)%i
          radt   = terms(iproc)%term(iterm)%rad(2)%i
          emitkt = min(radk,radt)
          call PrintCktSkt(iun,iterm,terms(iproc)%term(iterm))
!          call PrintCktSkt(6,iterm,terms(iproc)%term(iterm))
          call CktSktLimit(iun,pborn,p,phat,ptilde,Bij,Bijkl, &
                           flv_arr(:,iproc),                  &
                           emitkt,radk,radt,                  &
                           subterms_Cir_RR(iproc),            &
                           terms(iproc))
        elseif (subtype.eq.'CktCrktSkt    ') then
          radk   = terms(iproc)%term(iterm)%rad(1)%i
          radt   = terms(iproc)%term(iterm)%rad(2)%i
          emitkt = min(radk,radt)
          do radr=1,nleg
            if (abs(flv_arr(radr,iproc)).gt.qcd_nf) cycle
            if ((radr.eq.radk).or.(radr.eq.radt)) cycle
            emitktr = min(radr,radk,radt)
            radrID = flv_arr(radr,iproc)
            call PrintCktCrktSkt(iun,iterm,terms(iproc)%term(iterm), &
                                 emitktr,radr,radrID)
!            call PrintCktCrktSkt(6,iterm,terms(iproc)%term(iterm), &
!                                 emitktr,radr,radrID)
            call CktCrktSktLimit(iun,pborn,p,phat,ptilde,       &
                                 flv_arr(:,iproc),              &
                                 emitktr,emitkt,radr,radk,radt, &
                                 subterms_Cir_RR(iproc),        &
                                 terms(iproc))
          end do
        elseif (subtype.eq.'StCirt        ') then
          rad1      = terms(iproc)%term(iterm)%rad(1)%i
          rad2      = terms(iproc)%term(iterm)%rad(2)%i
          rad3      = terms(iproc)%term(iterm)%rad(3)%i
          rad1ID    = terms(iproc)%term(iterm)%rad(1)%ID
          rad2ID    = terms(iproc)%term(iterm)%rad(2)%ID
          rad3ID    = terms(iproc)%term(iterm)%rad(3)%ID
          emitktr   = terms(iproc)%term(iterm)%emit(1)%i
! Only gluon can become soft, hence leg t is selected according to this:
          if (rad1ID.eq.0) then
            radt = rad1
            radi = rad2
            radr = rad3
            radtID = rad1ID
            radiID = rad2ID
            radrID = rad3ID
            call PrintStCirt(iun,iterm,terms(iproc)%term(iterm), &
                             radi,radiID,radr,radrID,radt,radtID)
!            call PrintStCirt(6,iterm,terms(iproc)%term(iterm), &
!                             radi,radiID,radr,radrID,radt,radtID)
            call StCirtLimit(iun,pborn,p,phat,ptilde,Rij, &
                             flv_arr(:,iproc),            &
                             emitktr,radi,radr,radt,      &
                             subterms_Sr_RR(iproc),       &
                             terms(iproc))
          end if
          if (rad2ID.eq.0) then
            radt = rad2
            radi = rad1
            radr = rad3
            radtID = rad2ID
            radiID = rad1ID
            radrID = rad3ID
            call PrintStCirt(iun,iterm,terms(iproc)%term(iterm), &
                             radi,radiID,radr,radrID,radt,radtID)
!            call PrintStCirt(6,iterm,terms(iproc)%term(iterm), &
!                             radi,radiID,radr,radrID,radt,radtID)
            call StCirtLimit(iun,pborn,p,phat,ptilde,Rij, &
                             flv_arr(:,iproc),            &
                             emitktr,radi,radr,radt,      &
                             subterms_Sr_RR(iproc),       &
                             terms(iproc))
          end if
          if (rad3ID.eq.0) then
            radt = rad3
            radi = rad1
            radr = rad2
            radtID = rad3ID
            radiID = rad1ID
            radrID = rad2ID
            call PrintStCirt(iun,iterm,terms(iproc)%term(iterm), &
                             radi,radiID,radr,radrID,radt,radtID)
!            call PrintStCirt(6,iterm,terms(iproc)%term(iterm), &
!                             radi,radiID,radr,radrID,radt,radtID)
            call StCirtLimit(iun,pborn,p,phat,ptilde,Rij, &
                             flv_arr(:,iproc),            &
                             emitktr,radi,radr,radt,      &
                             subterms_Sr_RR(iproc),       &
                             terms(iproc))
          end if
        elseif (subtype.eq.'StCSirt       ') then
          radi   = terms(iproc)%term(iterm)%rad(1)%i
          radr   = terms(iproc)%term(iterm)%rad(2)%i
          radt   = terms(iproc)%term(iterm)%rad(3)%i
          emitir = terms(iproc)%term(iterm)%emit(1)%i
          call PrintStCSirt(iun,iterm,terms(iproc)%term(iterm))
!          call PrintStCSirt(6,iterm,terms(iproc)%term(iterm))
          call StCSirtLimit(iun,pborn,p,phat,ptilde,Bij,Rij,Bmunuij, &
                            flv_arr(:,iproc),                        &
                            emitir,radi,radr,radt,                   &
                            subterms_Sr_RR(iproc),                   &
                            terms(iproc))
        elseif (subtype.eq.'StCirtCSirt   ') then
          radi    = terms(iproc)%term(iterm)%rad(1)%i
          radr    = terms(iproc)%term(iterm)%rad(2)%i
          radt    = terms(iproc)%term(iterm)%rad(3)%i
          emitirt = terms(iproc)%term(iterm)%emit(1)%i
          call PrintStCirtCSirt(iun,iterm,terms(iproc)%term(iterm))
!          call PrintStCirtCSirt(6,iterm,terms(iproc)%term(iterm))
          call StCirtCSirtLimit(iun,pborn,p,phat,ptilde, &
                                flv_arr(:,iproc),        &
                                emitirt,radi,radr,radt,  &
                                subterms_Sr_RR(iproc),   &
                                terms(iproc))
        elseif (subtype.eq.'StCirtSrt     ') then
          do rad=1,2
            radr   = terms(iproc)%term(iterm)%rad(rad)%i
            radt   = terms(iproc)%term(iterm)%rad(mod(rad,2)+1)%i
            radrID = terms(iproc)%term(iterm)%rad(rad)%ID
            radtID = terms(iproc)%term(iterm)%rad(mod(rad,2)+1)%ID
! Contribution is only possible if parton t is a gluon,
! since a single soft emission is only possible then:
            if (flv_arr(radt,iproc).ne.0) cycle
            do radi=1,nleg
              if (abs(flv_arr(radi,iproc)).gt.qcd_nf) cycle
              if ((radi.eq.radr).or.(radi.eq.radt)) cycle
              emitirt = min(radi,radr,radt)
              radiID = flv_arr(radi,iproc)
              call PrintStCirtSrt(iun,iterm,terms(iproc)%term(iterm), &
                                  emitirt,radi,radiID, &
                                  radr,radrID,radt,radtID)
!              call PrintStCirtSrt(6,iterm,terms(iproc)%term(iterm), &
!                                  emitirt,radi,radiID, &
!                                  radr,radrID,radt,radtID)
              call StCirtSrtLimit(iun,pborn,p,phat,ptilde, &
                                  flv_arr(:,iproc),        &
                                  emitirt,radi,radr,radt,  &
                                  subterms_Sr_RR(iproc),   &
                                  terms(iproc))
            end do
          end do
        elseif (subtype.eq.'StCSirtSrt    ') then
          do rad=1,2
            radr   = terms(iproc)%term(iterm)%rad(rad)%i
            radt   = terms(iproc)%term(iterm)%rad(mod(rad,2)+1)%i
            radrID = terms(iproc)%term(iterm)%rad(rad)%ID
            radtID = terms(iproc)%term(iterm)%rad(mod(rad,2)+1)%ID
! Contribution is only possible if parton t is a gluon,
! since a single soft emission is only possible then:
            if (flv_arr(radt,iproc).ne.0) cycle
            do radi=1,nleg
              if (abs(flv_arr(radi,iproc)).gt.qcd_nf) cycle
              if ((radi.eq.radr).or.(radi.eq.radt)) cycle
              radiID = flv_arr(radi,iproc)
              emitir = min(radi,radr)
              if (emitir.gt.radt) emitir = emitir - 1
              call PrintStCSirtSrt(iun,iterm,terms(iproc)%term(iterm), &
                                   emitir,radi,radiID, &
                                   radr,radrID,radt,radtID)
!              call PrintStCSirtSrt(6,iterm,terms(iproc)%term(iterm), &
!                                   emitir,radi,radiID, &
!                                   radr,radrID,radt,radtID)
              call StCSirtSrtLimit(iun,pborn,p,phat,ptilde,Bij, &
                                   flv_arr(:,iproc),             &
                                   emitir,radi,radr,radt,        &
                                   subterms_Sr_RR(iproc),        &
                                   terms(iproc))
            end do
          end do
        elseif (subtype.eq.'StCirtCSirtSrt') then
          do rad=1,2
            radr = terms(iproc)%term(iterm)%rad(rad)%i
            radt = terms(iproc)%term(iterm)%rad(mod(rad,2)+1)%i
! This counter term is only possible if both r and t are gluons:
            if ((flv_arr(radr,iproc).ne.0).or.(flv_arr(radt,iproc).ne.0)) cycle
! Determining the position for radi:
            do radi=1,nleg
! Cannot coincide with the soft partons:
              if ((radi.eq.radr).or.(radi.eq.radt)) cycle
! It has to be a massless parton:
              if (abs(flv_arr(radi,iproc)).gt.qcd_nf) cycle
! Determining the position for the emitter of the i-r-s triplet,
! note that s is also soft!
              do emitirt=3,nleg-2
                if (terms(iproc)%term(iterm)%UBorn(emitirt).eq. &
                    (flv_arr(radi,iproc) + &
                     flv_arr(radr,iproc) + &
                     flv_arr(radt,iproc))) exit
              end do
              call PrintStCirtCSirtSrt(iun,iterm, &
                                       emitirt,radi,flv_arr(radi,iproc), &
                                       radr,radt,terms(iproc)%term(iterm))
!              call PrintStCirtCSirtSrt(6,iterm, &
!                                       emitirt,radi,flv_arr(radi,iproc), &
!                                       radr,radt,terms(iproc)%term(iterm))
              call StCirtCSirtSrtLimit(iun,pborn,p,phat,ptilde, &
                                       flv_arr(:,iproc),        &
                                       emitirt,radi,radr,radt,  &
                                       subterms_Sr_RR(iproc),   &
                                       terms(iproc))
            end do
          end do
        elseif (subtype.eq.'StSrt         ') then
          do rad=1,2
            radr = terms(iproc)%term(iterm)%rad(rad)%i
            radt = terms(iproc)%term(iterm)%rad(mod(rad,2)+1)%i
! This counter term is only possible if both r and t are gluons:
            if ((flv_arr(radr,iproc).ne.0).or.(flv_arr(radt,iproc).ne.0)) cycle
            call PrintStSrt(iun,iterm,radr,radt,terms(iproc)%term(iterm))
!            call PrintStSrt(6,iterm,radr,radt,terms(iproc)%term(iterm))
            call StSrtLimit(iun,pborn,p,phat,ptilde,Bij,Bijkl, &
                            flv_arr(:,iproc),                  &
                            radr,radt,                         &
                            subterms_Sr_RR(iproc),             &
                            terms(iproc))
          end do
        elseif (subtype.eq.'CitStCirt     ') then
! In the Cirs counterterm the legs are ordered: i < r < s, hence
! to have covered all CitStCirt terms the legs have to be permuted:
          do rad=1,3
            rad1 = terms(iproc)%term(iterm)%rad(rad)%i
            rad2 = terms(iproc)%term(iterm)%rad(mod(rad,3)+1)%i
            rad3 = terms(iproc)%term(iterm)%rad(mod(rad+1,3)+1)%i
! These are just half of the possible pairings, to have the other half
! we have to swap rad1 and rad3:
            do rad_pr=1,2
              if (rad_pr.eq.1) then
                radt = rad3
              else
                radt = rad1
                rad1 = rad3
              end if
! This counter term is only possible if t is gluon:
              if (flv_arr(radt,iproc).ne.0) cycle
              radi = rad1
              radr = rad2
              call PrintCitStCirt(iun,iterm,                &
                                  terms(iproc)%term(iterm), &
                                  radi,flv_arr(radi,iproc), &
                                  radr,flv_arr(radr,iproc), &
                                  radt,flv_arr(radt,iproc))
!               call PrintCitStCirt(6,iterm,                  &
!                                  terms(iproc)%term(iterm), &
!                                  radi,flv_arr(radi,iproc), &
!                                  radr,flv_arr(radr,iproc), &
!                                  radt,flv_arr(radt,iproc))
              call CitStCirtLimit(iun,pborn,p,phat,ptilde, &
                                  flv_arr(:,iproc),        &
                                  radi,radr,radt,          &
                                  subterms_CSir_RR(iproc), &
                                  terms(iproc))
            end do
          end do
        elseif (subtype.eq.'CktStCSirt    ') then
! The position for parton i, r and t is fixed due to the CS 
! counterterm:
          radi = terms(iproc)%term(iterm)%rad(1)%i
          radr = terms(iproc)%term(iterm)%rad(2)%i
          radt = terms(iproc)%term(iterm)%rad(3)%i
! Determining the position for radk:
          do radk=1,nleg
! Cannot coincide with the other rad legs:
            if ((radk.eq.radi).or. &
                (radk.eq.radr).or. &
                (radk.eq.radt)) cycle
! It has to be a massless parton:
            if (abs(flv_arr(radk,iproc)).gt.qcd_nf) cycle
            call PrintCktStCSirt(iun,iterm,                &
                                 terms(iproc)%term(iterm), &
                                 radk,flv_arr(radk,iproc), &
                                 radt,flv_arr(radt,iproc))
!            call PrintCktStCSirt(6,iterm,                  &
!                                 terms(iproc)%term(iterm), &
!                                 radk,flv_arr(radk,iproc), &
!                                 radt,flv_arr(radt,iproc))
            call CktStCSirtLimit(iun,pborn,p,phat,ptilde, &
                                 Bij,Bmunuij,             &
                                 flv_arr(:,iproc),        &
                                 radi,radr,radk,radt,     &
                                 subterms_CSir_RR(iproc), &
                                 terms(iproc))
          end do
        elseif (subtype.eq.'CktStCSirtSrt ') then
! The position for parton r and t is fixed due to the S 
! counterterm, but since the legs are ordered we have to permute them
! to get all the counterterms:
          do rad=1,2
            radr = terms(iproc)%term(iterm)%rad(rad)%i
            radt = terms(iproc)%term(iterm)%rad(mod(rad,2)+1)%i
! Only considering those Srt terms which are defined for
! gluons:
            if (flv_arr(radr,iproc).ne.0) cycle
! Determining the position for radi and radk:
            do radi=1,nleg
! Cannot coincide with the other rad legs:
              if ((radi.eq.radr).or. &
                  (radi.eq.radt)) cycle
! It has to be a massless parton:
              if (abs(flv_arr(radi,iproc)).gt.qcd_nf) cycle
              do radk=1,nleg
! Cannot coincide with the other rad legs:
                if ((radk.eq.radi).or. &
                    (radk.eq.radr).or. &
                    (radk.eq.radt)) cycle
! It has to be a massless parton:
                if (abs(flv_arr(radk,iproc)).gt.qcd_nf) cycle
! Determing the ir emitter on the underlying Born level:
                emitir = min(radi,radr)
                if (emitir.gt.radt) emitir = emitir - 1
                call PrintCktStCSirtSrt(iun,iterm,                &
                                        terms(iproc)%term(iterm), &
                                        emitir,                   &
                                        radi,flv_arr(radi,iproc), &
                                        radr,flv_arr(radr,iproc), &
                                        radk,flv_arr(radk,iproc), &
                                        radt,flv_arr(radt,iproc))
!                call PrintCktStCSirtSrt(6,iterm,                  &
!                                        terms(iproc)%term(iterm), &
!                                        emitir,                   &
!                                        radi,flv_arr(radi,iproc), &
!                                        radr,flv_arr(radr,iproc), &
!                                        radk,flv_arr(radk,iproc), &
!                                        radt,flv_arr(radt,iproc))
                call CktStCSirtSrtLimit(iun,pborn,p,phat,ptilde,Bij, &
                                        flv_arr(:,iproc),            &
                                        emitir,radi,radr,radk,radt,  &
                                        subterms_CSir_RR(iproc),     &
                                        terms(iproc))
              end do
            end do
          end do
        elseif (subtype.eq.'CktStCkrtSrt  ') then
! The position for parton r and t is fixed due to the S 
! counterterm, but a permutation is needed to cover all the 
! possible counterterms:
          do rad=1,2
            radr = terms(iproc)%term(iterm)%rad(rad)%i
            radt = terms(iproc)%term(iterm)%rad(mod(rad,2)+1)%i
! Only considering those Srt terms which are defined for
! gluons:
            if (flv_arr(radr,iproc).ne.0) cycle
! Determining the position for radk
            do radk=1,nleg
! Cannot coincide with the other rad legs:
              if ((radk.eq.radr).or. &
                  (radk.eq.radt)) cycle
! It has to be a massless parton:
              if (abs(flv_arr(radk,iproc)).gt.qcd_nf) cycle
! Determing the krt emitter on the underlying Born level:
! The triple collinear splitting remains the same no matter
! how we order the radiated partons, hence emitktr = emitkrt
              emitktr = min(radk,radr,radt)
              call PrintCktStCkrtSrt(iun,iterm,                &
                                     terms(iproc)%term(iterm), &
                                     emitktr,                  &
                                     radk,flv_arr(radk,iproc), &
                                     radr,flv_arr(radr,iproc), &
                                     radt,flv_arr(radt,iproc))
!              call PrintCktStCkrtSrt(6,iterm,                  &
!                                     terms(iproc)%term(iterm), &
!                                     emitktr,                  &
!                                     radk,flv_arr(radk,iproc), &
!                                     radr,flv_arr(radr,iproc), &
!                                     radt,flv_arr(radt,iproc))
              call CktStCkrtSrtLimit(iun,pborn,p,phat,ptilde, &
                                     flv_arr(:,iproc),        &
                                     emitktr,radk,radr,radt,  &
                                     subterms_CSir_RR(iproc), &
                                     terms(iproc))
            end do
          end do
        elseif (subtype.eq.'CrtStCkrtSrt  ') then
! The position for parton r and t is fixed due to the S 
! counterterm, a permutation is needed to get covered all the
! counterterms:
          do rad=1,2
            radr = terms(iproc)%term(iterm)%rad(rad)%i
            radt = terms(iproc)%term(iterm)%rad(mod(rad,2)+1)%i
! Only considering those Srt terms which are defined for
! gluons:
            if (flv_arr(radr,iproc).ne.0) cycle
! Determining the position for radk
            do radk=1,nleg
! Cannot coincide with the other rad legs:
              if ((radk.eq.radr).or. &
                  (radk.eq.radt)) cycle
! It has to be a massless parton:
              if (abs(flv_arr(radk,iproc)).gt.qcd_nf) cycle
! Determing the krt emitter on the underlying Born level:
! The triple collinear splitting remains the same no matter
! how we order the radiated partons, hence emitktr = emitkrt
              emitktr = min(radk,radr,radt)
              call PrintCrtStCkrtSrt(iun,iterm,                &
                                     terms(iproc)%term(iterm), &
                                     emitktr,                  &
                                     radk,flv_arr(radk,iproc), &
                                     radr,flv_arr(radr,iproc), &
                                     radt,flv_arr(radt,iproc))
!              call PrintCrtStCkrtSrt(6,iterm,                  &
!                                     terms(iproc)%term(iterm), &
!                                     emitktr,                  &
!                                     radk,flv_arr(radk,iproc), &
!                                     radr,flv_arr(radr,iproc), &
!                                     radt,flv_arr(radt,iproc))
              call CrtStCkrtSrtLimit(iun,pborn,p,phat,ptilde, &
                                     flv_arr(:,iproc),        &
                                     emitktr,radk,radr,radt,  &
                                     subterms_CSir_RR(iproc), &
                                     terms(iproc))
            end do
          end do
        elseif (subtype.eq.'CktStSrt      ') then
! The position for parton r and t is fixed due to the S 
! counterterm, but still have to permutate them to get all the
! possible counterterms:
          do rad=1,2
            radr = terms(iproc)%term(iterm)%rad(rad)%i
            radt = terms(iproc)%term(iterm)%rad(mod(rad,2)+1)%i
! Only considering those Srt terms which are defined for
! gluons:
            if (flv_arr(radr,iproc).ne.0) cycle
! Determining the position for radk
            do radk=1,nleg
! Cannot coincide with the other rad legs:
              if ((radk.eq.radr).or. &
                  (radk.eq.radt)) cycle
! It has to be a massless parton:
              if (abs(flv_arr(radk,iproc)).gt.qcd_nf) cycle
              call PrintCktStSrt(iun,iterm,                &
                                 terms(iproc)%term(iterm), &
                                 radk,flv_arr(radk,iproc), &
                                 radr,flv_arr(radr,iproc), &
                                 radt,flv_arr(radt,iproc))
!              call PrintCktStSrt(6,iterm,                  &
!                                 terms(iproc)%term(iterm), &
!                                 radk,flv_arr(radk,iproc), &
!                                 radr,flv_arr(radr,iproc), &
!                                 radt,flv_arr(radt,iproc))
              call CktStSrtLimit(iun,pborn,p,phat,ptilde,Bij,Bijkl, &
                                 flv_arr(:,iproc),                  &
                                 radk,radr,radt,                    &
                                 subterms_CSir_RR(iproc),           &
                                 terms(iproc))
            end do
          end do
        elseif (subtype.eq.'CrtStSrt      ') then
! The position for parton r and t is fixed due to the S 
! counterterm, but still need an additional permutation to
! get everything covered:
          do rad=1,2
            radr = terms(iproc)%term(iterm)%rad(rad)%i
            radt = terms(iproc)%term(iterm)%rad(mod(rad,2)+1)%i
! Only considering those Srt terms which are defined for
! gluons:
            if (flv_arr(radr,iproc).ne.0) cycle
            call PrintCrtStSrt(iun,iterm,                &
                               terms(iproc)%term(iterm), &
                               radr,flv_arr(radr,iproc), &
                               radt,flv_arr(radt,iproc))
!            call PrintCrtStSrt(6,iterm,                  &
!                               terms(iproc)%term(iterm), &
!                               radr,flv_arr(radr,iproc), &
!                               radt,flv_arr(radt,iproc))
            call CrtStSrtLimit(iun,pborn,p,phat,ptilde,Bij,Bijkl, &
                               flv_arr(:,iproc),                  &
                               radr,radt,                         &
                               subterms_CSir_RR(iproc),           &
                               terms(iproc))
          end do
        else
          print *,"Wrong subtraction term type is specified..."
          print *,"cont: ",cont
          print *,"subtype: ",subtype
          stop "CheckNNLOlimits"
        end if
      else
        print *,"Wrong contribution type is specified..."
        print *,"cont: ",cont
        stop "CheckNNLOlimits"
      end if
    end do
  end do
!
  deallocate(pborn,p,ptilde,phat,Bij,Bijk,Bijkl,Bmunuij,Rij)
!
!  stop "CheckNNLOlimits"
!
end subroutine CheckNNLOlimits
!
subroutine CirLimit(iun,pborn,p,ptilde,Bij,flv,emit,rad,Cir,Sr,CSir, &
                    CalcSMEB,CalcSMEBmunu,CalcSMEBij,CalcSMER)
use process
use particles
use math
use random
use FKSutils
use regions
use scales
use my_scale
use my_model
use coupling
implicit none
!
  integer , intent(in) :: iun
  type(mom) , dimension(:) , intent(in) :: pborn
  type(particle) , dimension(:) , intent(inout) :: p,ptilde
  real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
  integer , dimension(:) , intent(in) :: flv
  integer , intent(in) :: emit,rad
  type(subterms) , intent(in) :: Cir,Sr,CSir
!
  integer :: iexp,iscale
  integer :: istat
  integer :: nleg
  integer :: emit1,rad1
  real(kind(1d0)) :: xi,y,azi
  real(kind(1d0)) :: smeR,Cirterm
  real(kind(1d0)) , dimension(:) , allocatable :: A1term,smeR_arr
!
!
  interface 
!
    subroutine CalcSubInvariants(p)
    use regions
    use particles
!
      type(particle) , dimension(:) , intent(in) :: p
!
    end subroutine CalcSubInvariants

!
    subroutine CalcA1subtraction(p,ptilde,Bij,weight, &
                                 Cir,Sr,CSir, &
                                 CalcSMEB,CalcSMEBmunu,CalcSMEBij, &
                                 A1term)
    use regions
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      type(particle) , dimension(:) , intent(out) :: ptilde
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Bij
      real(kind(1d0)) , intent(in) :: weight
      type(subterms) , intent(in) :: Cir,Sr,CSir
      real(kind(1d0)) , dimension(:) , intent(out) :: A1term
!
      interface
        subroutine CalcSMEB(p,smeB)
        use particles
        implicit none
!
          type(particle) , dimension(:) , intent(in) :: p
          real(kind(1d0)) , intent(out) :: smeB
!
        end subroutine CalcSMEB
!
        subroutine CalcSMEBmunu(ileg,p,Bmunu)
        use particles
        implicit none
!
          integer , intent(in) :: ileg
          type(particle) , dimension(:) , intent(in) :: p
          real(kind(1d0)) , dimension(0:3,0:3) , intent(out) :: Bmunu
!
        end subroutine CalcSMEBmunu
!
        subroutine CalcSMEBij(p,Bij)
        use particles
        implicit none
!
          type(particle) , dimension(:) , intent(in) :: p
          real(kind(1d0)) , dimension(:,:) , intent(out) :: Bij
!
        end subroutine CalcSMEBij
      end interface
!
    end subroutine CalcA1subtraction
!
    subroutine PickCir(p,ptilde,emit,rad,Cir,CalcSMEB,CalcSMEBmunu, &
                       Cirterm)
    use regions
    use particles
    implicit none
!
      type(particle), dimension(:), intent(in) :: p
      type(particle), dimension(:), intent(out) :: ptilde
      integer , intent(in) :: emit,rad
      type(subterms) , intent(in) :: Cir
      real(kind(1d0)) , intent(out) :: Cirterm
!
      interface
        subroutine CalcSMEB(p,smeB)
        use particles
        implicit none
!
          type(particle) , dimension(:) , intent(in) :: p
          real(kind(1d0)) , intent(out) :: smeB
!
        end subroutine CalcSMEB
!
        subroutine CalcSMEBmunu(ileg,p,Bmunu)
        use particles
        implicit none
!
          integer , intent(in) :: ileg
          type(particle) , dimension(:) , intent(in) :: p
          real(kind(1d0)) , dimension(0:3,0:3) , intent(out) :: Bmunu
!
        end subroutine CalcSMEBmunu
      end interface
    end subroutine PickCir
!
    subroutine CalcSMEB(p,smeB)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: smeB
!
    end subroutine CalcSMEB
!
    subroutine CalcSMEBmunu(ileg,p,Bmunu)
    use particles
    implicit none
!
      integer , intent(in) :: ileg
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(0:3,0:3) , intent(out) :: Bmunu
!
    end subroutine CalcSMEBmunu
!
    subroutine CalcSMEBij(p,Bij)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Bij
!
    end subroutine CalcSMEBij
!
    subroutine CalcSMER(p,smeR)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: smeR
!
    end subroutine CalcSMER
  end interface
!
!  print *,"************************************************************"
!  print *,"We will approach the limit of: "
!  print *,"emit,rad: ",emit,rad
!  call PrintSubProc(flv)
!  print *,"************************************************************"
!  read(*,*)
!
  allocate(A1term(nscales), &
           smeR_arr(nscales), &
           stat=istat)
  if (istat.ne.0) then
    print *,"Problem with allocation of A1term in CirLimit..."
    stop
  end if
!
! The number of legs for the real process can be determined by
! examining the size of the momentum array holding the momenta:
  nleg = size(p)
! To go into the collinear limit we simply fix (randomly)
! xi and the azimuth:
  xi = gen_rand()
  azi = 2d0*pi*gen_rand()
!  print *,"New sequence is started: "
  emit1 = min(emit,rad)
  rad1  = max(emit,rad)
  do iexp=1,nexp
    y = 1d0 - 10d0**(-iexp)
! We construct the real momenta:
    call MapFSR(pborn,xi,y,azi,emit1,rad1,flv,p)
! We have to calculate the common invariants for the real
! kinematics:
    call CalcSubInvariants(p)
! We have to calculate the scales:
    call calcmyscales(p)
! we obtain the SME for the Real: 
    call CalcSMER(p,smeR)
! we calculate the real contribution with different scales:
    call StoreXis_scales
    do iscale=1,nscales
      call ChangeXis_scales(iscale)
      call calcmyscales(p)
      call calc_couplings
      smeR_arr(iscale) = smeR &
                       * alphas**((nleg - nleg_born)+border_as) &
                       * alphaEM**border_aEM
    end do
    call RestoreXis_scales
! We calculate the A1 type subtraction:
!    call CalcA1subtraction(p,ptilde,Bij,1d0, &
!                           Cir,Sr,CSir, &
!                           CalcSMEB,CalcSMEBmunu,CalcSMEBij, &
!                           A1term)
! We calculate only that particular Cir term which actually cancels the
! leading divergence:
    call PickCir(p,ptilde,emit,rad,Cir,CalcSMEB,CalcSMEBmunu, &
                 Cirterm)
!    print *,A1term / smeR
    write(iun,*) "iexp= ",iexp, &
                 ", Cir/M2= ",Cirterm/smeR
!
!    Print *,"Sequence starts here: "
!    print *,"orig. smeR: ",smeR
!    print *,"A1: ",A1term
!    print *,"smeR: ",smeR_arr
!    do iscale=1,nscales
!      print *,A1term(iscale)/smeR_arr(iscale)
!    end do
!
  end do
!
  deallocate(A1term,smeR_arr)
!
end subroutine CirLimit
!
subroutine SrLimit(iun,pborn,p,ptilde,Bij,flv,rad,Cir,Sr,CSir, &
                   CalcSMEB,CalcSMEBmunu,CalcSMEBij,CalcSMER)
use process
use particles
use math
use random
use FKSutils
use regions
use QCDparams
use scales
use my_scale
use coupling
use misc
implicit none
!
  integer , intent(in) :: iun
  type(mom) , dimension(:) , intent(in) :: pborn
  type(particle) , dimension(:) , intent(inout) :: p,ptilde
  real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
  integer , dimension(:) , intent(in) :: flv
  integer , intent(in) :: rad
  type(subterms) , intent(in) :: Cir,Sr,CSir
!
  integer :: iexp,iscale
  integer :: istat
  integer :: nleg
  integer :: emit
  real(kind(1d0)) :: xi,y,azi
  real(kind(1d0)) :: smeR,Srterm
  real(kind(1d0)) :: fluxfact
  real(kind(1d0)) , dimension(:) , allocatable :: A1term,smeR_arr
!
!
  interface 
!
    subroutine CalcSubInvariants(p)
    use regions
    use particles
!
      type(particle) , dimension(:) , intent(in) :: p
!
    end subroutine CalcSubInvariants

!
    subroutine CalcA1subtraction(p,ptilde,Bij,weight, &
                                 Cir,Sr,CSir, &
                                 CalcSMEB,CalcSMEBmunu,CalcSMEBij, &
                                 A1term)
    use regions
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      type(particle) , dimension(:) , intent(out) :: ptilde
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Bij
      real(kind(1d0)) , intent(in) :: weight
      type(subterms) , intent(in) :: Cir,Sr,CSir
      real(kind(1d0)) , dimension(:) , intent(out) :: A1term
!
      interface
        subroutine CalcSMEB(p,smeB)
        use particles
        implicit none
!
          type(particle) , dimension(:) , intent(in) :: p
          real(kind(1d0)) , intent(out) :: smeB
!
        end subroutine CalcSMEB
!
        subroutine CalcSMEBmunu(ileg,p,Bmunu)
        use particles
        implicit none
!
          integer , intent(in) :: ileg
          type(particle) , dimension(:) , intent(in) :: p
          real(kind(1d0)) , dimension(0:3,0:3) , intent(out) :: Bmunu
!
        end subroutine CalcSMEBmunu
!
        subroutine CalcSMEBij(p,Bij)
        use particles
        implicit none
!
          type(particle) , dimension(:) , intent(in) :: p
          real(kind(1d0)) , dimension(:,:) , intent(out) :: Bij
!
        end subroutine CalcSMEBij
      end interface
!
    end subroutine CalcA1subtraction
!
    subroutine PickSr(p,ptilde,Bij,radr,Sr,CalcSMEBij,Srterm)
    use regions
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      type(particle) , dimension(:) , intent(out) :: ptilde
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Bij
      integer , intent(in) :: radr
      type(subterms) , intent(in) :: Sr
      real(kind(1d0)) , intent(out) :: Srterm
!
      interface
        subroutine CalcSMEBij(p,Bij)
        use particles
        implicit none
!
          type(particle) , dimension(:) , intent(in) :: p
          real(kind(1d0)) , dimension(:,:) , intent(out) :: Bij
!
        end subroutine CalcSMEBij
      end interface
!
    end subroutine PickSr
!
    subroutine CalcSMEB(p,smeB)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: smeB
!
    end subroutine CalcSMEB
!
    subroutine CalcSMEBmunu(ileg,p,Bmunu)
    use particles
    implicit none
!
      integer , intent(in) :: ileg
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(0:3,0:3) , intent(out) :: Bmunu
!
    end subroutine CalcSMEBmunu
!
    subroutine CalcSMEBij(p,Bij)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Bij
!
    end subroutine CalcSMEBij
!
    subroutine CalcSMER(p,smeR)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: smeR
!
    end subroutine CalcSMER
  end interface
!
!  print *,"************************************************************"
!  print *,"We will approach the limit of: "
!  print *,"rad: ",rad
!  call PrintSubProc(flv)
!  print *,"************************************************************"
!  read(*,*)
!
  allocate(A1term(nscales), &
           smeR_arr(nscales), &
           stat=istat)
  if (istat.ne.0) then
    print *,"Problem with allocation of A1term and smeR_arr in SrLimit..."
    stop
  end if
!
! The number of legs for the real process can be determined by
! examining the size of the momentum array holding the momenta:
  nleg = size(p)
! To go into the soft limit we simply fix (randomly)
! y and the azimuth:
!  y = 1d0 - 2d0*gen_rand()
  y = 1d0
  azi = 2d0*pi*gen_rand()
!  print *,"New sequence is started: "
  do iexp=1,nexp
    xi = 10d0**(-iexp)
! In order to construct the real momenta we have to pick a final state
! emitter, to do so we select the first massless final state parton
! which does not coincide with the emitted one:
    do emit=3,size(pborn)
! massless parton:
      if (abs(flv(emit)).gt.qcd_nf) cycle
! does not coincide with the soft gluon:
      if (emit.eq.rad) cycle
      exit
    end do
! We construct the real momenta:
    call MapFSR(pborn,xi,y,azi,emit,rad,flv,p)
! We have to calculate the common invariants for the real
! kinematics:
    call CalcSubInvariants(p)
! We have to calculate the scales:
    call calcmyscales(p)
! we obtain the SME for the Real: 
    call CalcSMER(p,smeR)
! Calculating the flux factor:
    fluxfact = CalcFluxFact(p)
! we calculate the real contribution with different scales:
    call StoreXis_scales
    do iscale=1,nscales
      call ChangeXis_scales(iscale)
      call calcmyscales(p)
      call calc_couplings
      smeR_arr(iscale) = smeR &
                       * fluxfact &
                       * alphas**((nleg - nleg_born)*border_as) &
                       * alphaEM**border_aEM
    end do
    call RestoreXis_scales
! We calculate the A1 type subtraction:
!    call CalcA1subtraction(p,ptilde,Bij,1d0, &
!                           Cir,Sr,CSir, &
!                           CalcSMEB,CalcSMEBmunu,CalcSMEBij, &
!                           A1term)
! We calculate only that particular Sr term which actually cancels the
! leading divergence:
    call PickSr(p,ptilde,Bij,rad,Sr,CalcSMEBij,Srterm)
!    print *,A1term / smeR
    write(iun,*) "iexp= ",iexp, &
!                 ", A1/M2= ",A1term(centscale)/smeR_arr(centscale), &
                 ", Sr/M2= ",Srterm/smeR
!
!    Print *,"Sequence starts here: "
!    print *,"orig. smeR: ",smeR
!    print *,"A1: ",A1term
!    print *,"smeR: ",smeR_arr
!    do iscale=1,nscales
!      print *,A1term(iscale)/smeR_arr(iscale)
!    end do
  end do
!
  deallocate(A1term,smeR_arr)
!
end subroutine SrLimit
!
subroutine CirSrLimit(iun,pborn,p,ptilde,Bij,flv,emit,rad,Cir,Sr,CSir, &
                      CalcSMEB,CalcSMEBmunu,CalcSMEBij,CalcSMER)
use process
use particles
use math
use random
use FKSutils
use regions
use QCDparams
use scales
use my_scale
use coupling
implicit none
!
  integer , intent(in) :: iun
  type(mom) , dimension(:) , intent(in) :: pborn
  type(particle) , dimension(:) , intent(inout) :: p,ptilde
  real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
  integer , dimension(:) , intent(in) :: flv
  integer , intent(in) :: emit,rad
  type(subterms) , intent(in) :: Cir,Sr,CSir
!
  integer :: iexp,iscale
  integer :: istat
  integer :: nleg
  real(kind(1d0)) :: xi,y,azi
  real(kind(1d0)) :: smeR,CSirterm
  real(kind(1d0)) , dimension(:) , allocatable :: A1term,smeR_arr
!
!
  interface 
!
    subroutine CalcSubInvariants(p)
    use regions
    use particles
!
      type(particle) , dimension(:) , intent(in) :: p
!
    end subroutine CalcSubInvariants

!
    subroutine CalcA1subtraction(p,ptilde,Bij,weight, &
                                 Cir,Sr,CSir, &
                                 CalcSMEB,CalcSMEBmunu,CalcSMEBij, &
                                 A1term)
    use regions
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      type(particle) , dimension(:) , intent(out) :: ptilde
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Bij
      real(kind(1d0)) , intent(in) :: weight
      type(subterms) , intent(in) :: Cir,Sr,CSir
      real(kind(1d0)) , dimension(:) , intent(out) :: A1term
!
      interface
        subroutine CalcSMEB(p,smeB)
        use particles
        implicit none
!
          type(particle) , dimension(:) , intent(in) :: p
          real(kind(1d0)) , intent(out) :: smeB
!
        end subroutine CalcSMEB
!
        subroutine CalcSMEBmunu(ileg,p,Bmunu)
        use particles
        implicit none
!
          integer , intent(in) :: ileg
          type(particle) , dimension(:) , intent(in) :: p
          real(kind(1d0)) , dimension(0:3,0:3) , intent(out) :: Bmunu
!
        end subroutine CalcSMEBmunu
!
        subroutine CalcSMEBij(p,Bij)
        use particles
        implicit none
!
          type(particle) , dimension(:) , intent(in) :: p
          real(kind(1d0)) , dimension(:,:) , intent(out) :: Bij
!
        end subroutine CalcSMEBij
      end interface
!
    end subroutine CalcA1subtraction
!
    subroutine PickCirSr(p,ptilde,emit,rad,CSir,CalcSMEB,CSirterm)
    use regions
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      type(particle) , dimension(:) , intent(out) :: ptilde
      integer , intent(in) :: emit,rad
      type(subterms) , intent(in) :: CSir
      real(kind(1d0)) , intent(out) :: CSirterm
!
      interface
        subroutine CalcSMEB(p,smeB)
        use particles
        implicit none
!
          type(particle) , dimension(:) , intent(in) :: p
          real(kind(1d0)) , intent(out) :: smeB
!
        end subroutine CalcSMEB
      end interface
!
    end subroutine PickCirSr
!
    subroutine CalcSMEB(p,smeB)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: smeB
!
    end subroutine CalcSMEB
!
    subroutine CalcSMEBmunu(ileg,p,Bmunu)
    use particles
    implicit none
!
      integer , intent(in) :: ileg
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(0:3,0:3) , intent(out) :: Bmunu
!
    end subroutine CalcSMEBmunu
!
    subroutine CalcSMEBij(p,Bij)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Bij
!
    end subroutine CalcSMEBij
!
    subroutine CalcSMER(p,smeR)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: smeR
!
    end subroutine CalcSMER
  end interface
!
!  print *,"************************************************************"
!  print *,"We will approach the limit of: "
!  print *,"emit,rad: ",emit,rad
!  call PrintSubProc(flv)
!  print *,"************************************************************"
!  read(*,*)
!
  allocate(A1term(nscales), &
           smeR_arr(nscales), &
           stat=istat)
  if (istat.ne.0) then
    print *,"Problem with allocation of A1term in CirSrLimit..."
    stop
  end if
!
! The number of legs for the real process can be determined by
! examining the size of the momentum array holding the momenta:
  nleg = size(p)
! We only fix the azimuth randomly:
  azi = 2d0*pi*gen_rand()
!  print *,"New sequence is started: "
  do iexp=1,nexp
! This is the collinear-soft limit, hence we have to approach
! the boundary of both xi and y simultaneously:
    xi = 10d0**(-iexp)
    y  = 1d0 - 10d0**(-iexp)
! We construct the real momenta:
    call MapFSR(pborn,xi,y,azi,emit,rad,flv,p)
! We have to calculate the common invariants for the real
! kinematics:
    call CalcSubInvariants(p)
! We have to calculate the scales:
    call calcmyscales(p)
! we obtain the SME for the Real: 
    call CalcSMER(p,smeR)
! we calculate the real contribution with different scales:
    call StoreXis_scales
    do iscale=1,nscales
      call ChangeXis_scales(iscale)
      call calcmyscales(p)
      call calc_couplings
      smeR_arr(iscale) = smeR &
                       * alphas**((nleg - nleg_born)*border_as) &
                       * alphaEM**border_aEM
    end do
! We calculate the A1 type subtraction:
    call CalcA1subtraction(p,ptilde,Bij,1d0, &
                           Cir,Sr,CSir, &
                           CalcSMEB,CalcSMEBmunu,CalcSMEBij, &
                           A1term)
! We calculate only that particular Sr term which actually cancels the
! leading divergence:
    call PickCirSr(p,ptilde,emit,rad,CSir,CalcSMEB,CSirterm)
!    print *,A1term / smeR
    write(iun,*) "iexp= ",iexp, &
!                 ", A1/M2= ",A1term/smeR, &
                 ", CirSr/M2= ",CSirterm/smeR
!
!    Print *,"Sequence starts here: "
!    print *,"orig. smeR: ",smeR
!    print *,"A1: ",A1term
!    print *,"smeR: ",smeR_arr
!    do iscale=1,nscales
!      print *,A1term(iscale)/smeR_arr(iscale)
!    end do
  end do
!
  deallocate(A1term,smeR_arr)
!
end subroutine CirSrLimit
!
subroutine CirLimitRV(iun,pborn,p,ptilde,flv,emit,rad,Cir)
use particles
use random
use regions
use math
use FKSutils
use scales
use my_scale
use coupling
implicit none
!
  integer , intent(in) :: iun
  type(mom) , dimension(:) , intent(in) :: pborn
  type(particle) , dimension(:) , intent(inout) :: p,ptilde
  integer , dimension(:) , intent(in) :: flv
  integer , intent(in) :: emit,rad
  type(subterms) , intent(in) :: Cir
!
  integer :: iexp,iscale
  integer :: emit1,rad1
  real(kind(1d0)) :: xi,y,azi
  real(kind(1d0)) :: s_ir
  real(kind(1d0)) :: smeRV,Cirterm
!
!
  interface 
    subroutine CalcSubInvariants(p)
    use regions
    use particles
!
      type(particle) , dimension(:) , intent(in) :: p
!
    end subroutine CalcSubInvariants
!
    subroutine PickCirRV(p,ptilde,emit,rad,Cir,Cirterm)
    use regions
    use particles
    implicit none
!
      type(particle), dimension(:), intent(in) :: p
      type(particle), dimension(:), intent(out) :: ptilde
      integer , intent(in) :: emit,rad
      type(subterms) , intent(in) :: Cir
      real(kind(1d0)) , intent(out) :: Cirterm
!
    end subroutine PickCirRV
!
    subroutine CalcRV(parts,smeRV,smeRVLaurent)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: parts
      real(kind(1d0)) , intent(out) :: smeRV
      real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: &
        smeRVLaurent
!
    end subroutine CalcRV
  end interface
!
!  print *,"************************************************************"
!  print *,"We will approach the limit of: "
!  print *,"emit,rad: ",emit,rad
!  call PrintSubProc(flv)
!  print *,"************************************************************"
!  read(*,*)
!
! To go into the collinear limit we simply fix (randomly)
! xi and the azimuth:
  xi = gen_rand()
  azi = 2d0*pi*gen_rand()
!  print *,"New sequence is started: "
  emit1 = min(emit,rad)
  rad1  = max(emit,rad)
  do iexp=1,nexp
    y = 1d0 - 10d0**(-iexp)
! We construct the real momenta:
    call MapFSR(pborn,xi,y,azi,emit1,rad1,flv,p)
    s_ir = 2d0*p(emit)%p*p(rad)%p
! We have to calculate the common invariants for the real
! kinematics:
!    call PrintParts(p)
    call CalcSubInvariants(p)
! We have to calculate the scales:
    call calcmyscales(p)
! The next line is only a simple test implies \mu_R^2 = Q^2:
!    mur = sqrt(2d0*p(1)%p*p(2)%p)
! we obtain the SME for the real-virtual: 
    call CalcRV(p,smeRV)
! We calculate only that particular Cir term which actually cancels the
! leading divergence:
    call PickCirRV(p,ptilde,emit,rad,Cir,Cirterm)
!    write(iun,*) "iexp= ",int(log10(s_ir)), &
!                 ", Cir/M2= ",Cirterm/smeRV
    write(iun,fmt='(a,I3,1x)',advance='no') "iexp= ",int(log10(s_ir))
    write(iun,*) ", Cir/M2= ",Cirterm/smeRV
!
  end do
!
end subroutine CirLimitRV
!
subroutine SrLimitRV(iun,pborn,p,ptilde,Bij,Vij,Bijk,flv,rad,Sr)
use particles
use random
use regions
use math
use FKSutils
use scales
use my_scale
use coupling
use QCDparams
implicit none
!
  integer , intent(in) :: iun
  type(mom) , dimension(:) , intent(in) :: pborn
  type(particle) , dimension(:) , intent(inout) :: p,ptilde
  real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
  real(kind(1d0)) , dimension(:,:) , intent(inout) :: Vij
  real(kind(1d0)) , dimension(:,:,:) , intent(inout) :: Bijk
  integer , dimension(:) , intent(in) :: flv
  integer , intent(in) :: rad
  type(subterms) , intent(in) :: Sr
!
  integer :: iexp,iscale
  integer :: istat
  integer :: nleg
  integer :: emit
  real(kind(1d0)) :: xi,y,azi
  real(kind(1d0)) :: Eg
  real(kind(1d0)) :: smeRV,Srterm
  real(kind(1d0)) , dimension(-4:2) :: RVLaurent
!
!
  interface 
    subroutine CalcSubInvariants(p)
    use regions
    use particles
!
      type(particle) , dimension(:) , intent(in) :: p
!
    end subroutine CalcSubInvariants
!
    subroutine PickSrRV(p,ptilde,Bij,Vij,Bijk,rad,Sr,Srterm)
    use regions
    use particles
    implicit none
!
      type(particle), dimension(:), intent(in) :: p
      type(particle), dimension(:), intent(out) :: ptilde
      real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
      real(kind(1d0)) , dimension(:,:) , intent(inout) :: Vij
      real(kind(1d0)) , dimension(:,:,:) , intent(inout) :: Bijk
      integer , intent(in) :: rad
      type(subterms) , intent(in) :: Sr
      real(kind(1d0)) , intent(out) :: Srterm
!
    end subroutine PickSrRV
!
    subroutine CalcRV(parts,smeRV,smeRVLaurent)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: parts
      real(kind(1d0)) , intent(out) :: smeRV
      real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: &
        smeRVLaurent
!
    end subroutine CalcRV
  end interface
!
!  print *,"************************************************************"
!  print *,"We will approach the limit of: "
!  print *,"rad: ",rad
!  call PrintSubProc(flv)
!  print *,"************************************************************"
!  read(*,*)
!
! To go into the soft limit we simply fix 
! y and the azimuth, azimuth is fixed randomly:
  y = 0d0
  azi = 2d0*pi*gen_rand()
!  print *,"New sequence is started: "
  do iexp=1,nexp
    xi = 10d0**(-iexp)
! In order to construct the real momenta we have to pick a final state
! emitter, to do so we select the first massless final state parton
! which does not coincide with the emitted one:
    do emit=3,size(pborn)
! massless parton:
      if (abs(flv(emit)).gt.qcd_nf) cycle
! does not coincide with the soft gluon:
      if (emit.eq.rad) cycle
      exit
    end do
! We construct the real momenta:
    call MapFSR(pborn,xi,y,azi,emit,rad,flv,p)
    Eg = p(rad)%p%E
! We have to calculate the common invariants for the real
! kinematics:
!    call PrintParts(p)
    call CalcSubInvariants(p)
! We have to calculate the scales:
    call calcmyscales(p)
! The next line is only a simple test implies \mu_R^2 = Q^2:
!    mur = sqrt(2d0*p(1)%p*p(2)%p)
! we obtain the SME for the real-virtual: 
    call CalcRV(p,smeRV)
! We calculate only that particular Cir term which actually cancels the
! leading divergence:
    call PickSrRV(p,ptilde,Bij,Vij,Bijk,rad,Sr,Srterm)
!    write(iun,*) "iexp= ",int(log10(Eg)), &
!                 ", Sr/M2= ",Srterm/smeRV
    write(iun,fmt='(a,I3,1x)',advance='no') "iexp= ",int(log10(Eg))
    write(iun,*) ", Sr/M2= ",Srterm/smeRV
!
  end do
!
end subroutine SrLimitRV
!
subroutine CirSrLimitRV(iun,pborn,p,ptilde,flv,emit,rad,CirSr)
use particles
use math
use random
use FKSutils
use regions
use scales
use my_scale
implicit none
!
  integer , intent(in) :: iun
  type(mom) , dimension(:) , intent(in) :: pborn
  type(particle) , dimension(:) , intent(inout) :: p,ptilde
  integer , dimension(:) , intent(in) :: flv
  integer , intent(in) :: emit,rad
  type(subterms) , intent(in) :: CirSr
!
  integer :: iexp,iscale
  real(kind(1d0)) :: xi,y,azi
  real(kind(1d0)) :: s_ir,Eg
  real(kind(1d0)) :: smeRV,CirSrterm
!
!
  interface 
    subroutine CalcSubInvariants(p)
    use regions
    use particles
!
      type(particle) , dimension(:) , intent(in) :: p
!
    end subroutine CalcSubInvariants
!
    subroutine PickCirSrRV(p,ptilde,emit,rad,CirSr,CirSrterm)
    use regions
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      type(particle) , dimension(:) , intent(out) :: ptilde
      integer , intent(in) :: emit,rad
      type(subterms) , intent(in) :: CirSr
      real(kind(1d0)) , intent(out) :: CirSrterm
!
    end subroutine PickCirSrRV
!
    subroutine CalcRV(parts,smeRV,smeRVLaurent)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: parts
      real(kind(1d0)) , intent(out) :: smeRV
      real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: &
        smeRVLaurent
!
    end subroutine CalcRV
  end interface
!
!  print *,"************************************************************"
!  print *,"We will approach the limit of: "
!  print *,"emit,rad: ",emit,rad
!  call PrintSubProc(flv)
!  print *,"************************************************************"
!  read(*,*)
!
! We only fix the azimuth randomly:
  azi = 2d0*pi*gen_rand()
!  print *,"New sequence is started: "
  do iexp=1,nexp
! This is the collinear-soft limit, hence we have to approach
! the boundary of both xi and y simultaneously:
    xi = 10d0**(-iexp)
    y  = 1d0 - 10d0**(-iexp)
! We construct the real momenta:
    call MapFSR(pborn,xi,y,azi,emit,rad,flv,p)
    s_ir = 2d0*p(emit)%p*p(rad)%p
    Eg   = p(rad)%p%E
! We have to calculate the common invariants for the real
! kinematics:
    call CalcSubInvariants(p)
! We have to calculate the scales:
    call calcmyscales(p)
! The next line is only a simple test implies \mu_R^2 = Q^2:
!    mur = sqrt(2d0*p(1)%p*p(2)%p)
! we obtain the SME for the Real: 
    call CalcRV(p,smeRV)
! We calculate only that particular CirSr term which actually cancels the
! leading divergence:
    call PickCirSrRV(p,ptilde,emit,rad,CirSr,CirSrterm)
!    write(iun,*) "iexp= ",int(log10(s_ir)),int(log10(Eg)), &
!                 ", CirSr/M2= ",CirSrterm/smeRV
    write(iun,fmt='(a,I3,1x,I3,1x)',advance='no') "iexp= ",int(log10(s_ir)),int(log10(Eg))
    write(iun,*) ", CirSr/M2= ",CirSrterm/smeRV
!
  end do
!
end subroutine CirSrLimitRV
!
subroutine CirLimitRRA1(iun,pborn,p,ptilde,Bij,Bmunuij,Rij,flv,emit,rad,Cir)
use particles
use random
use regions
use math
use FKSutils
use scales
use my_scale
use coupling
implicit none
!
  integer , intent(in) :: iun
  type(mom) , dimension(:) , intent(in) :: pborn
  type(particle) , dimension(:) , intent(inout) :: p,ptilde
  real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij,Rij
  real(kind(1d0)) , dimension(0:,0:,:,:) , intent(inout) :: Bmunuij
  integer , dimension(:) , intent(in) :: flv
  integer , intent(in) :: emit,rad
  type(subterms) , intent(in) :: Cir
!
  integer :: iexp
  integer :: emit1,rad1
  real(kind(1d0)) :: xi,y,azi
  real(kind(1d0)) :: s_ir
  real(kind(1d0)) :: I1term,Cirterm
!
!
  interface 
    subroutine CalcSubInvariants(p)
    use regions
    use particles
!
      type(particle) , dimension(:) , intent(in) :: p
!
    end subroutine CalcSubInvariants
!
    subroutine PickCirRRA1(p,ptilde,Bij,Bmunuij,emit,rad,Cir,Cirterm)
    use regions
    use particles
    implicit none
!
      type(particle), dimension(:), intent(in) :: p
      type(particle), dimension(:), intent(out) :: ptilde
      real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
      real(kind(1d0)) , dimension(0:,0:,:,:) , intent(inout) :: Bmunuij
      integer , intent(in) :: emit,rad
      type(subterms) , intent(in) :: Cir
      real(kind(1d0)) , intent(out) :: Cirterm
!
    end subroutine PickCirRRA1
!
    subroutine CalcR(p,smeR)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: smeR
!
    end subroutine CalcR
!
    subroutine CalcRij(p,Rij)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Rij
!
    end subroutine CalcRij
!
    subroutine CalcI1(p,Bij,CalcSMEB,CalcSMEBij,I1term,I1Laurent)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
      real(kind(1d0)) , intent(out) :: I1term
      real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: I1Laurent
!
      interface
        subroutine CalcSMEB(p,smeB)
        use particles
        implicit none
!
          type(particle) , dimension(:) , intent(in) :: p
          real(kind(1d0)) , intent(out) :: smeB
!
        end subroutine CalcSMEB
!
        subroutine CalcSMEBij(p,Bij)
        use particles
        implicit none
!
          type(particle) , dimension(:) , intent(in) :: p
          real(kind(1d0)) , dimension(:,:) , intent(out) :: Bij
!
        end subroutine CalcSMEBij
      end interface
    end subroutine CalcI1
  end interface
!
!  print *,"************************************************************"
!  print *,"We will approach the limit of: "
!  print *,"emit,rad: ",emit,rad
!  call PrintSubProc(flv)
!  print *,"************************************************************"
!  read(*,*)
!
! To go into the collinear limit we simply fix (randomly)
! xi and the azimuth:
  xi = gen_rand()
  azi = 2d0*pi*gen_rand()
!  print *,"New sequence is started: "
  emit1 = min(emit,rad)
  rad1  = max(emit,rad)
  do iexp=1,nexp
    y = 1d0 - 10d0**(-iexp)
! We construct the real momenta:
    call MapFSR(pborn,xi,y,azi,emit1,rad1,flv,p)
    s_ir = 2d0*p(emit)%p*p(rad)%p
! We have to calculate the common invariants for the real
! kinematics:
!    call PrintParts(p)
    call CalcSubInvariants(p)
! We have to calculate the scales:
    call calcmyscales(p)
! The next line is only a simple test implies \mu_R^2 = Q^2:
!    mur = sqrt(2d0*p(1)%p*p(2)%p)
! We obtain the RR_A1 contribution for the subprocess in question:
    call CalcI1(p,Rij,CalcR,CalcRij,I1term)
! We calculate only that particular Cir term which actually cancels the
! leading divergence:
    call PickCirRRA1(p,ptilde,Bij,Bmunuij,emit,rad,Cir,Cirterm)
!    write(iun,*) "iexp= ",int(log10(s_ir)), &
!                 ", Cir/I1= ",Cirterm/I1term
    write(iun,fmt='(a,I3,1x)',advance='no') "iexp= ",int(log10(s_ir))
    write(iun,*) ", Cir/I1= ",Cirterm/I1term
!
  end do
!
end subroutine CirLimitRRA1
!
subroutine SrLimitRRA1(iun,pborn,p,ptilde,Bij,Bijkl,Rij,flv,rad,Sr)
use particles
use random
use regions
use math
use FKSutils
use scales
use my_scale
use QCDparams
implicit none
!
  integer , intent(in) :: iun
  type(mom) , dimension(:) , intent(in) :: pborn
  type(particle) , dimension(:) , intent(inout) :: p,ptilde
  real(kind(1d0)) , dimension(:,:) , allocatable , intent(inout) :: Bij,Rij
  real(kind(1d0)) , dimension(:,:,:,:) , allocatable , intent(inout) :: Bijkl
  integer , dimension(:) , intent(in) :: flv
  integer , intent(in) :: rad
  type(subterms) , intent(in) :: Sr
!
  integer :: iexp
  integer :: emit
  real(kind(1d0)) :: xi,y,azi
  real(kind(1d0)) :: Eg
  real(kind(1d0)) :: I1term,Srterm
!
!
  interface 
    subroutine CalcSubInvariants(p)
    use regions
    use particles
!
      type(particle) , dimension(:) , intent(in) :: p
!
    end subroutine CalcSubInvariants
!
    subroutine PickSrRRA1(p,ptilde,Bij,Bijkl,rad,Sr,Srterm)
    use regions
    use particles
    implicit none
!
      type(particle), dimension(:), intent(in) :: p
      type(particle), dimension(:), intent(out) :: ptilde
      real(kind(1d0)) , dimension(:,:) , allocatable , intent(inout) :: Bij
      real(kind(1d0)) , dimension(:,:,:,:) , allocatable , intent(inout) :: Bijkl
      integer , intent(in) :: rad
      type(subterms) , intent(in) :: Sr
      real(kind(1d0)) , intent(out) :: Srterm
!
    end subroutine PickSrRRA1
!
    subroutine CalcR(p,smeR)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: smeR
!
    end subroutine CalcR
!
    subroutine CalcRij(p,Rij)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Rij
!
    end subroutine CalcRij
!
    subroutine CalcI1(p,Bij,CalcSMEB,CalcSMEBij,I1term,I1Laurent)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
      real(kind(1d0)) , intent(out) :: I1term
      real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: I1Laurent
!
      interface
        subroutine CalcSMEB(p,smeB)
        use particles
        implicit none
!
          type(particle) , dimension(:) , intent(in) :: p
          real(kind(1d0)) , intent(out) :: smeB
!
        end subroutine CalcSMEB
!
        subroutine CalcSMEBij(p,Bij)
        use particles
        implicit none
!
          type(particle) , dimension(:) , intent(in) :: p
          real(kind(1d0)) , dimension(:,:) , intent(out) :: Bij
!
        end subroutine CalcSMEBij
      end interface
    end subroutine CalcI1
  end interface
!
!  print *,"************************************************************"
!  print *,"We will approach the limit of: "
!  print *,"rad: ",rad
!  call PrintSubProc(flv)
!  print *,"************************************************************"
!  read(*,*)
!
! To go into the soft limit we simply fix 
! y and the azimuth, the azimuth is fixed randomly:
  y = 1d0
  azi = 2d0*pi*gen_rand()
!  print *,"New sequence is started: "
  do iexp=1,nexp
    xi = 10d0**(-iexp)
! In order to construct the real momenta we have to pick a final state
! emitter, to do so we select the first massless final state parton
! which does not coincide with the emitted one:
    do emit=3,size(pborn)
! massless parton:
      if (abs(flv(emit)).gt.qcd_nf) cycle
! does not coincide with the soft gluon:
      if (emit.eq.rad) cycle
      exit
    end do
! We construct the real momenta:
    call MapFSR(pborn,xi,y,azi,emit,rad,flv,p)
    Eg = p(rad)%p%E
! We have to calculate the common invariants for the real
! kinematics:
!    call PrintParts(p)
    call CalcSubInvariants(p)
! We have to calculate the scales:
    call calcmyscales(p)
! The next line is only a simple test implies \mu_R^2 = Q^2:
!    mur = sqrt(2d0*p(1)%p*p(2)%p)
! We obtain the RR_A1 contribution for the subprocess in question:
    call CalcI1(p,Rij,CalcR,CalcRij,I1term)
! We calculate only that particular Sr term which actually cancels the
! leading divergence:
    call PickSrRRA1(p,ptilde,Bij,Bijkl,rad,Sr,Srterm)
!    write(iun,*) "iexp= ",int(log10(Eg)), &
!                 ", Sr/I1= ",Srterm/I1term
    write(iun,fmt='(a,I3,1x)',advance='no') "iexp= ",int(log10(Eg))
    write(iun,*) ", Sr/M2= ",Srterm/I1term
!
  end do
!
end subroutine SrLimitRRA1
!
subroutine CirSrLimitRRA1(iun,pborn,p,ptilde,Bij,Rij,flv,emit,rad,CirSr)
use particles
use random
use regions
use math
use FKSutils
use scales
use my_scale
use coupling
implicit none
!
  integer , intent(in) :: iun
  type(mom) , dimension(:) , intent(in) :: pborn
  type(particle) , dimension(:) , intent(inout) :: p,ptilde
  real(kind(1d0)) , dimension(:,:) , allocatable , intent(inout) :: Bij,Rij
  integer , dimension(:) , intent(in) :: flv
  integer , intent(in) :: emit,rad
  type(subterms) , intent(in) :: CirSr
!
  integer :: iexp
  real(kind(1d0)) :: xi,y,azi
  real(kind(1d0)) :: s_ir,Eg
  real(kind(1d0)) :: I1term,CirSrterm
!
!
  interface 
    subroutine CalcSubInvariants(p)
    use regions
    use particles
!
      type(particle) , dimension(:) , intent(in) :: p
!
    end subroutine CalcSubInvariants
!
    subroutine PickCirSrRRA1(p,ptilde,Bij,emit,rad,CirSr,CirSrterm)
    use regions
    use particles
    implicit none
!
      type(particle), dimension(:), intent(in) :: p
      type(particle), dimension(:), intent(out) :: ptilde
      real(kind(1d0)) , dimension(:,:) , allocatable , intent(inout) :: Bij
      integer , intent(in) :: emit,rad
      type(subterms) , intent(in) :: CirSr
      real(kind(1d0)) , intent(out) :: CirSrterm
!
    end subroutine PickCirSrRRA1
!
    subroutine CalcR(p,smeR)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: smeR
!
    end subroutine CalcR
!
    subroutine CalcRij(p,Rij)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Rij
!
    end subroutine CalcRij
!
    subroutine CalcI1(p,Bij,CalcSMEB,CalcSMEBij,I1term,I1Laurent)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
      real(kind(1d0)) , intent(out) :: I1term
      real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: I1Laurent
!
      interface
        subroutine CalcSMEB(p,smeB)
        use particles
        implicit none
!
          type(particle) , dimension(:) , intent(in) :: p
          real(kind(1d0)) , intent(out) :: smeB
!
        end subroutine CalcSMEB
!
        subroutine CalcSMEBij(p,Bij)
        use particles
        implicit none
!
          type(particle) , dimension(:) , intent(in) :: p
          real(kind(1d0)) , dimension(:,:) , intent(out) :: Bij
!
        end subroutine CalcSMEBij
      end interface
    end subroutine CalcI1
  end interface
!
!  print *,"************************************************************"
!  print *,"We will approach the limit of: "
!  print *,"emit,rad: ",emit,rad
!  call PrintSubProc(flv)
!  print *,"************************************************************"
!  read(*,*)
!
! Only fixing the azimuth randomly:
  azi = 2d0*pi*gen_rand()
!  print *,"New sequence is started: "
  do iexp=1,nexp
! This is the collinear-soft limit, hence we have to approach
! the boundary of both xi and y simultaneously:
    xi = 10d0**(-iexp)
    y = 1d0 - 10d0**(-iexp)
! We construct the real momenta:
    call MapFSR(pborn,xi,y,azi,emit,rad,flv,p)
    s_ir = 2d0*p(emit)%p*p(rad)%p
    Eg   = p(rad)%p%E
! We have to calculate the common invariants for the real
! kinematics:
!    call PrintParts(p)
    call CalcSubInvariants(p)
! We have to calculate the scales:
    call calcmyscales(p)
! The next line is only a simple test implies \mu_R^2 = Q^2:
!    mur = sqrt(2d0*p(1)%p*p(2)%p)
! We obtain the RR_A1 contribution for the subprocess in question:
    call CalcI1(p,Rij,CalcR,CalcRij,I1term)
! We calculate only that particular Cir term which actually cancels the
! leading divergence:
    call PickCirSrRRA1(p,ptilde,Bij,emit,rad,CirSr,CirSrterm)
!    write(iun,*) "iexp= ",int(log10(s_ir)),int(log10(Eg)), &
!                 ", CirSr/I1= ",CirSrterm/I1term
    write(iun,fmt='(a,I3,1x,I3,1x)',advance='no') "iexp= ",int(log10(s_ir)),int(log10(Eg))
    write(iun,*) ", CirSr/M2= ",CirSrterm/I1term
!
  end do
!
end subroutine CirSrLimitRRA1
!
subroutine CirsLimit(iun,pborn,p,phat,ptilde,flv,emit,radi,radr,rads,Cirs)
use particles
use random
use regions
use math
use FKSutils
use scales
use my_scale
use coupling
implicit none
!
  integer , intent(in) :: iun
  type(mom) , dimension(:) , intent(in) :: pborn
  type(particle) , dimension(:) , intent(inout) :: p,ptilde,phat
  integer , dimension(:) , intent(in) :: flv
  integer , intent(in) :: emit,radi,radr,rads
  type(subterms) , intent(in) :: Cirs
!
  integer :: iexp
  integer :: rad1,rad2
  real(kind(1d0)) :: xi1,xi2,azi1,azi2,y
  real(kind(1d0)) :: s_ir,s_is,s_irs
  real(kind(1d0)) :: smeRR,Cirsterm
!
!
  interface 
    subroutine CalcSubInvariants(p)
    use regions
    use particles
!
      type(particle) , dimension(:) , intent(in) :: p
!
    end subroutine CalcSubInvariants
!
    subroutine PickCirs(p,ptilde,emit,radi,radr,rads,Cirs,Cirsterm)
    use regions
    use particles
    implicit none
!
      type(particle), dimension(:), intent(in) :: p
      type(particle), dimension(:), intent(out) :: ptilde
      integer , intent(in) :: emit,radi,radr,rads
      type(subterms) , intent(in) :: Cirs
      real(kind(1d0)) , intent(out) :: Cirsterm
!
    end subroutine PickCirs
!
    subroutine CalcRR(p,smeRR)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: smeRR
!
    end subroutine CalcRR
  end interface
!
!  print *,"************************************************************"
!  print *,"We will approach the limit of: "
!  print *,"emit,radi,radr,rads: ",emit,radi,radr,rads
!  call PrintSubProc(flv)
!  print *,"************************************************************"
!  read(*,*)
!
! The triply collinear limit is approached by applying the FSR mapping
! twice in each step, hence two energy fractions and azimuths have to be
! specified (randomly):
  xi1  = gen_rand()
  xi2  = gen_rand()
  azi1 = 2d0*pi*gen_rand()
  azi2 = 2d0*pi*gen_rand()
! Having one emitter and three radiated partons it has to be specified
! which radiated one corresponds to the emitter:
  if (emit.eq.radi) then
    rad1 = min(radr,rads)
    rad2 = max(radr,rads)
  elseif (emit.eq.radr) then
    rad1 = min(radi,rads)
    rad2 = max(radi,rads)
  elseif (emit.eq.rads) then
    rad1 = min(radi,radr)
    rad2 = max(radi,radr)
  end if
!  print *,"New sequence is started: "
  do iexp=1,nexp
! The scaling is chosen to be the same:
    y = 1d0 - 10d0**(-iexp)
! Using the Born momenta to construct an intermediate PS (phat):
    call MapFSR(pborn,xi1,y,azi1,emit,rad1,flv,phat)
! Using the intermediate momenta to construct the final RR PS: 
    y = 1d0 - 10d0**(-iexp)
    call MapFSR(phat(:)%p,xi2,y,azi2,emit,rad2,flv,p)
    s_ir = 2d0*p(emit)%p*p(rad1)%p
    s_is = 2d0*p(emit)%p*p(rad2)%p
    s_irs = (p(radi)%p + p(radr)%p + p(rads)%p)**2
! We have to calculate the common invariants for the real
! kinematics:
!    call PrintParts(phat)
!    call PrintParts(p)
!    if ((radi.eq.5).and.(radr.eq.6).and.(rads.eq.7)) then
!      print *,"Cggg is reached..."
!      call DebugRRealPS(p(:)%p)
!      call PrintParts(p)
!    end if
    call CalcSubInvariants(p)
! We have to calculate the scales:
    call calcmyscales(p)
! The next line is only a simple test implies \mu_R^2 = Q^2:
!    mur = sqrt(2d0*p(1)%p*p(2)%p)
! We obtain the RR contribution for the subprocess in question:
    call CalcRR(p,smeRR)
!    print *,"smeRR: ",smeRR
! We calculate only that particular Cirs term which actually cancels the
! leading divergence:
    call PickCirs(p,ptilde,emit,radi,radr,rads,Cirs,Cirsterm)
    write(iun,fmt='(a,I3,1x)',advance='no') "iexp= ",int(log10(s_irs))
    write(iun,*) ", Cirs/RR= ",Cirsterm/smeRR
!
  end do
!
end subroutine CirsLimit
!
subroutine CirjsLimit(iun,pborn,p,phat,ptilde,flv, &
                      emiti,radi,radr,emitj,radj,rads, &
                      Cirjs)
use particles
use random
use regions
use math
use FKSutils
use scales
use my_scale
use coupling
implicit none
!
  integer , intent(in) :: iun
  type(mom) , dimension(:) , intent(in) :: pborn
  type(particle) , dimension(:) , intent(inout) :: p,ptilde,phat
  integer , dimension(:) , intent(in) :: flv
  integer , intent(in) :: emiti,emitj,radi,radr,radj,rads
  type(subterms) , intent(in) :: Cirjs
!
  integer :: iexp
  integer :: emit1,emit2,rad1,rad2
  real(kind(1d0)) :: xi1,xi2,azi1,azi2,y
  real(kind(1d0)) :: s_ir,s_js
  real(kind(1d0)) :: smeRR,Cirjsterm
!
!
  interface 
    subroutine CalcSubInvariants(p)
    use regions
    use particles
!
      type(particle) , dimension(:) , intent(in) :: p
!
    end subroutine CalcSubInvariants
!
    subroutine PickCirjs(p,ptilde,emiti,radi,radr,emitj,radj,rads, &
                         Cirjs,Cirjsterm)
    use regions
    use particles
    implicit none
!
      type(particle), dimension(:), intent(in) :: p
      type(particle), dimension(:), intent(out) :: ptilde
      integer , intent(in) :: emiti,radi,radr,emitj,radj,rads
      type(subterms) , intent(in) :: Cirjs
      real(kind(1d0)) , intent(out) :: Cirjsterm
!
    end subroutine PickCirjs
!
    subroutine CalcRR(p,smeRR)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: smeRR
!
    end subroutine CalcRR
  end interface
!
!  print *,"************************************************************"
!  print *,"We will approach the limit of: "
!  print *,"emiti,radi,radr: ",emiti,radi,radr
!  print *,"emitj,radj,rads: ",emitj,radj,rads
!  call PrintSubProc(flv)
!  print *,"************************************************************"
!  read(*,*)
!
! The double collinear limit is approached by applying the FSR mapping
! twice in each step, hence two energy fractions and azimuths have to be
! specified (randomly):
  xi1  = gen_rand()
  xi2  = gen_rand()
  azi1 = 2d0*pi*gen_rand()
  azi2 = 2d0*pi*gen_rand()
! In the case of the double collinear limit we have two
! emitters which are ordered, that is emiti < emitj, but
! not for the radiated ones. To build up phase space firstly 
! that branching has to be done where the radiated parton 
! comes first.
! Starting with the hypothesis that the first radiation is
! indeed coming from the first pair, the true radiated parton
! is that one for which the position is higher:
  rad1 = max(radi,radr)
  rad2 = max(radj,rads)
! The hypothesis was correct the first branching happens before
! the second one:
  if (rad1.lt.rad2) then
    emit1 = emiti
    emit2 = emitj
! The ordering of the branchings is just the other way around:
  else
    rad1 = max(radj,rads)
    rad2 = max(radi,radr)
    emit1 = emitj
    emit2 = emiti
  end if
! If the first radiated parton is before the second emitter the
! position of the second emitter has to be increased by one unit
! since the emitter numbering corresponds to the underlying Born:
  if (rad1.le.emit2) emit2 = emit2 + 1
!  print *,"New sequence is started: "
  do iexp=1,nexp
! The scaling is chosen to be the same:
    y = 1d0 - 10d0**(-iexp)
! Using the Born momenta to construct an intermediate PS (phat):
    call MapFSR(pborn,xi1,y,azi1,emit1,rad1,flv,phat)
! Using the intermediate momenta to construct the final RR PS: 
    y = 1d0 - 10d0**(-iexp)
    call MapFSR(phat(:)%p,xi2,y,azi2,emit2,rad2,flv,p)
    s_ir = 2*p(radi)%p*p(radr)%p
    s_js = 2*p(radj)%p*p(rads)%p
!    print *,"s_ir: ",s_ir
!    print *,"s_js: ",s_js
! We have to calculate the common invariants for the real
! kinematics:
!    call PrintParts(phat)
!    call PrintParts(p)
    call CalcSubInvariants(p)
! We have to calculate the scales:
    call calcmyscales(p)
! The next line is only a simple test implies \mu_R^2 = Q^2:
!    mur = sqrt(2d0*p(1)%p*p(2)%p)
! We obtain the RR contribution for the subprocess in question:
    call CalcRR(p,smeRR)
!    print *,"smeRR: ",smeRR
! We calculate only that particular Cirs term which actually cancels the
! leading divergence:
    call PickCirjs(p,ptilde,emiti,radi,radr,emitj,radj,rads, &
                   Cirjs,Cirjsterm)
    write(iun,fmt='(a,2(I3,1x))',advance='no') "iexp= ",int(log10(s_ir)),int(log10(s_js))
    write(iun,*) ", Cirjs/RR= ",Cirjsterm/smeRR
!
  end do
!
end subroutine CirjsLimit
!
subroutine CSirsLimit(iun,pborn,p,phat,ptilde,Bij,Bmunuij, &
                      flv,emit,radi,radr,rads,CSirs)
use particles
use random
use regions
use math
use FKSutils
use scales
use my_scale
use coupling
use QCDparams
use misc
implicit none
!
  integer , intent(in) :: iun
  type(mom) , dimension(:) , intent(in) :: pborn
  type(particle) , dimension(:) , intent(inout) :: p,ptilde,phat
  real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
  real(kind(1d0)) , dimension(0:,0:,:,:) , intent(inout) :: Bmunuij
  integer , dimension(:) , intent(in) :: flv
  integer , intent(in) :: emit,radi,radr,rads
  type(subterms) , intent(in) :: CSirs
!
  logical :: soft_first,coll_first
  integer :: iexp
  integer :: emit1,emit2,rad1,rad2
  real(kind(1d0)) :: xi1,xi2,azi1,azi2,y1,y2
  real(kind(1d0)) :: s_ir,Es
  real(kind(1d0)) :: smeRR,CSirsterm
!
!
  interface 
    subroutine CalcSubInvariants(p)
    use regions
    use particles
!
      type(particle) , dimension(:) , intent(in) :: p
!
    end subroutine CalcSubInvariants
!
    subroutine PickCSirs(p,phat,ptilde,Bij,Bmunuij,emit,radi,radr,rads, &
                         CSirs,CSirsterm)
    use regions
    use particles
    implicit none
!
      type(particle), dimension(:), intent(in) :: p
      type(particle), dimension(:), intent(out) :: phat,ptilde
      real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
      real(kind(1d0)) , dimension(0:,0:,:,:) , intent(inout) :: Bmunuij
      integer , intent(in) :: emit,radi,radr,rads
      type(subterms) , intent(in) :: CSirs
      real(kind(1d0)) , intent(out) :: CSirsterm
!
    end subroutine PickCSirs
!
    subroutine CalcRR(p,smeRR)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: smeRR
!
    end subroutine CalcRR
  end interface
!
!  print *,"************************************************************"
!  print *,"We will approach the limit of: "
!  print *,"emit,radi,radr,rads: ",emit,radi,radr,rads
!  call PrintSubProc(flv)
!  print *,"************************************************************"
!  read(*,*)
!
! In the case of the double soft-collinear limit azimuths are
! fixed randomly.
  azi1 = 2d0*pi*gen_rand()
  azi2 = 2d0*pi*gen_rand()
! The phase space construction depends upon the order of legs,
! if the soft leg comes first, the soft branching is considered first:
  soft_first = .false.
! Soft radiation comes first:
  if (max(radi,radr).gt.rads) then
    soft_first = .true.
! For the soft the enclosed angle is fixed such that it is perpendicular to
! the emitter:
    y1  = 1d0
! For the collinear pair the energy is randomly reinstated:
    xi2 = 0.5d0
    rad1 = rads
    rad2 = max(radi,radr)
    do emit1=3,size(pborn)
      if (emit1.eq.emit) cycle
      if (abs(flv(emit1)).gt.qcd_nf) cycle
      exit
    end do
! If the soft radiation happens first position of the emitter
! for the collinear pair can change, since if the soft parton
! is inserted before the emitter of the collinear pair the
! position of that is got shifted:
    if (rads.le.emit) then
      emit2 = emit + 1
    else 
      emit2 = emit
    end if
! Collinear radiation comes first:
  else
! For the soft the enclosed angle is fixed:
    y2  = 0d0
    xi1 = 0.5d0
    rad1 = max(radi,radr)
    rad2 = rads
    emit1 = emit
    do emit2=3,size(pborn)
      if (emit2.eq.radi) cycle
      if (emit2.eq.radr) cycle
      if (abs(flv(emit2)).gt.qcd_nf) cycle
      exit
    end do
  end if
!  print *,"New sequence is started: "
  do iexp=1,nexp
! Using the Born momenta to construct an intermediate PS (phat):
! If soft radiation comes first xi1 and y2 has to be scaled:
    if (soft_first) then
      xi1 = 10d0**(-iexp/2d0)
      xi1 = (xi1*(1 - 2*sqrt(tiny_y)) + sqrt(tiny_y) - tiny_xi) &
          / (1 - 2*tiny_xi)
      y2  = 1d0 - 10d0**(-iexp)
    else
      y1  = 1d0 - 10d0**(-iexp)
      xi2 = 10d0**(-iexp/2d0)
      xi2 = (xi2*(1 - 2*sqrt(tiny_y)) + sqrt(tiny_y) - tiny_xi) &
          / (1 - 2*tiny_xi)
    end if
    call MapFSR(pborn,xi1,y1,azi1,emit1,rad1,flv,phat)
    call MapFSR(phat(:)%p,xi2,y2,azi2,emit2,rad2,flv,p)
    s_ir = 2*p(radi)%p*p(radr)%p
    Es   = p(rads)%p%E
! We have to calculate the common invariants for the real
! kinematics:
    call CalcSubInvariants(p)
! We have to calculate the scales:
    call calcmyscales(p)
! We obtain the RR contribution for the subprocess in question:
    call CalcRR(p,smeRR)
! We calculate only that particular Cirs term which actually cancels the
! leading divergence:
    call PickCSirs(p,phat,ptilde,Bij,Bmunuij,emit,radi,radr,rads, &
                   CSirs,CSirsterm)
    write(iun,fmt='(a,2(I3,1x))',advance='no') "iexp= ", &
      int(log10(s_ir)),int(log10(Es))
    write(iun,*) ", CSirs/RR= ",CSirsterm/smeRR
!
  end do
!
end subroutine CSirsLimit
!
subroutine CirsCSirsLimit(iun,pborn,p,phat,ptilde,Bij,Bmunuij, &
                          flv,emit,radi,radr,rads,CSirs)
use particles
use random
use regions
use math
use FKSutils
use scales
use my_scale
use coupling
use QCDparams
use misc
use nnlo_mom_maps
implicit none
!
  integer , intent(in) :: iun
  type(mom) , dimension(:) , intent(in) :: pborn
  type(particle) , dimension(:) , intent(inout) :: p,ptilde,phat
  real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
  real(kind(1d0)) , dimension(0:,0:,:,:) , intent(inout) :: Bmunuij
  integer , dimension(:) , intent(in) :: flv
  integer , intent(in) :: emit,radi,radr,rads
  type(subterms) , intent(in) :: CSirs
!
  integer :: iexp
  integer :: emiti,rad1,rad2
  real(kind(1d0)) :: xi1,xi2,azi1,azi2,y
  real(kind(1d0)) :: s_ir,s_is,s_rs,s_irs
  real(kind(1d0)) :: CSirsTerm,CirsCSirsTerm
!
!
  interface 
    subroutine CalcSubInvariants(p)
    use regions
    use particles
!
      type(particle) , dimension(:) , intent(in) :: p
!
    end subroutine CalcSubInvariants
!
    subroutine PickCirsCSirs(p,phat,ptilde,emit,radi,radr,rads, &
                             CSirs,CirsCSirsterm)
    use regions
    use particles
    implicit none
!
      type(particle), dimension(:), intent(in) :: p
      type(particle), dimension(:), intent(out) :: phat,ptilde
      integer , intent(in) :: emit,radi,radr,rads
      type(subterms) , intent(in) :: CSirs
      real(kind(1d0)) , intent(out) :: CirsCSirsterm
!
    end subroutine PickCirsCSirs
!
    subroutine PickCSirs(p,phat,ptilde,Bij,Bmunuij,emit,radi,radr,rads, &
                         CSirs,CSirsterm)
    use regions
    use particles
    implicit none
!
      type(particle), dimension(:), intent(in) :: p
      type(particle), dimension(:), intent(out) :: phat,ptilde
      real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
      real(kind(1d0)) , dimension(0:,0:,:,:) , intent(inout) :: Bmunuij
      integer , intent(in) :: emit,radi,radr,rads
      type(subterms) , intent(in) :: CSirs
      real(kind(1d0)) , intent(out) :: CSirsterm
!
    end subroutine PickCSirs
  end interface
!
!  print *,"************************************************************"
!  print *,"We will approach the limit of: "
!  print *,"emit,radi,radr,rads: ",emit,radi,radr,rads
!  call PrintSubProc(flv)
!  print *,"************************************************************"
!  read(*,*)
! For the CirsCSirs counter term a triply collinear limit
! has to be approached, to do so the azimuths and energies
! are selected randomly:
  xi1  = gen_rand()
  xi2  = gen_rand()
  azi1 = 2d0*pi*gen_rand()
  azi2 = 2d0*pi*gen_rand()
! Since for CSirs the ordering among the positions for the
! radiated partons can be anything the emitter is defined
! by the position of the least high radiated parton:
  emiti = min(radi,radr,rads)
! In the triply collinear limit two partons are emitted from
! the same source hence to sequencially create the PS the
! position of the radiated partons has to be selected with
! care, that is the first radiation produces a parton at
! position which is the second smallest after the position
! of the emitter:
  if (emiti.eq.radi) then
    rad1 = min(radr,rads)
    rad2 = max(radr,rads)
  elseif (emiti.eq.radr) then
    rad1 = min(radi,rads)
    rad2 = max(radi,rads)
  elseif (emiti.eq.rads) then
    rad1 = min(radi,radr)
    rad2 = max(radi,radr)
  end if
!  print *,"New sequence is started: "
  do iexp=1,nexp
    y = 1d0 - 10d0**(-iexp)
    call MapFSR(pborn,xi1,y,azi1,emiti,rad1,flv,phat)
    y = 1d0 - 10d0**(-iexp)
    call MapFSR(phat(:)%p,xi2,y,azi2,emiti,rad2,flv,p)
    s_ir  = 2*p(radi)%p*p(radr)%p
    s_is  = 2*p(radi)%p*p(rads)%p
    s_rs  = 2*p(radr)%p*p(rads)%p
    s_irs = (p(radi)%p + p(radr)%p + p(rads)%p)**2
!    print *,"s_ir,s_is,s_rs,s_irs: ",s_ir,s_is,s_rs,s_irs
!    print *,"emiti: ",emiti
!    print *,"rad1,rad2: ",rad1,rad2
! We have to calculate the common invariants for the real
! kinematics:
!    call PrintParts(phat)
!    call PrintParts(p)
!    call CheckParts(phat)
!    call CheckParts(p)
    call CalcSubInvariants(p)
! We have to calculate the scales:
    call calcmyscales(p)
! The next line is only a simple test implies \mu_R^2 = Q^2:
!    mur = sqrt(2d0*p(1)%p*p(2)%p)
! Calculating the CSirs term:
    call PickCSirs(p,phat,ptilde,Bij,Bmunuij,emit,radi,radr,rads, &
                   CSirs,CSirsTerm)
! We calculate only that particular Cirs term which actually cancels the
! leading divergence:
    call PickCirsCSirs(p,phat,ptilde,emit,radi,radr,rads, &
                       CSirs,CirsCSirsTerm)
    write(iun,fmt='(a,I3,1x)',advance='no') "iexp= ", &
      int(log10(s_irs))
    write(iun,*) ", CirsCSirs/CSirs= ",CirsCSirsTerm/CSirsTerm
!
  end do
!
end subroutine CirsCSirsLimit
!
subroutine CirjsCSirsLimit(iun,pborn,p,phat,ptilde,Bij,Bmunuij, &
                           flv,emitir,radi,radr,emitjs,radj,rads,CSirs)
use particles
use random
use regions
use math
use FKSutils
use scales
use my_scale
use coupling
use QCDparams
use misc
use nnlo_mom_maps
implicit none
!
  integer , intent(in) :: iun
  type(mom) , dimension(:) , intent(in) :: pborn
  type(particle) , dimension(:) , intent(inout) :: p,ptilde,phat
  real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
  real(kind(1d0)) , dimension(0:,0:,:,:) , intent(inout) :: Bmunuij
  integer , dimension(:) , intent(in) :: flv
  integer , intent(in) :: emitir,emitjs,radi,radr,radj,rads
  type(subterms) , intent(in) :: CSirs
!
  integer :: iexp
  integer :: emit1,emit2,rad1,rad2
  real(kind(1d0)) :: xi1,xi2,azi1,azi2,y
  real(kind(1d0)) :: s_ir,s_js
  real(kind(1d0)) :: CSirsTerm,CirjsCSirsTerm
!
!
  interface 
    subroutine CalcSubInvariants(p)
    use regions
    use particles
!
      type(particle) , dimension(:) , intent(in) :: p
!
    end subroutine CalcSubInvariants
!
    subroutine PickCirjsCSirs(p,phat,ptilde, &
                              emitir,radi,radr, &
                              emitjs,radj,rads, &
                              CSirs,CirjsCSirsterm)
    use regions
    use particles
    implicit none
!
      type(particle), dimension(:), intent(in) :: p
      type(particle), dimension(:), intent(out) :: phat,ptilde
      integer , intent(in) :: emitir,emitjs,radi,radr,radj,rads
      type(subterms) , intent(in) :: CSirs
      real(kind(1d0)) , intent(out) :: CirjsCSirsterm
!
    end subroutine PickCirjsCSirs
!
    subroutine PickCSirs(p,phat,ptilde,Bij,Bmunuij,emit,radi,radr,rads, &
                         CSirs,CSirsterm)
    use regions
    use particles
    implicit none
!
      type(particle), dimension(:), intent(in) :: p
      type(particle), dimension(:), intent(out) :: phat,ptilde
      real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
      real(kind(1d0)) , dimension(0:,0:,:,:) , intent(inout) :: Bmunuij
      integer , intent(in) :: emit,radi,radr,rads
      type(subterms) , intent(in) :: CSirs
      real(kind(1d0)) , intent(out) :: CSirsterm
!
    end subroutine PickCSirs
  end interface
!
!  print *,"************************************************************"
!  print *,"We will approach the limit of: "
!  print *,"emitir,emitjs,radi,radr,radj,rads: ",emitir,emitjs,radi,radr,radj,rads
!  call PrintSubProc(flv)
!  print *,"************************************************************"
!  read(*,*)
! For the CirjsCSirs counter term a doubly collinear limit
! has to be approached, to do so the azimuths and energies
! are selected randomly:
  xi1  = gen_rand()
  xi2  = gen_rand()
  azi1 = 2d0*pi*gen_rand()
  azi2 = 2d0*pi*gen_rand()
! To construct a PS point in the doubly collinear limit the collinear
! emissions have to be ordered:
  if (max(radi,radr).lt.max(radj,rads)) then
    emit1 = emitir
    rad1  = max(radi,radr)
    emit2 = emitjs
    rad2  = max(radj,rads)
  else
    emit1 = emitjs
    rad1  = max(radj,rads)
    emit2 = emitir
    rad2  = max(radi,radr)
  end if
  if (rad1.le.emit2) emit2 = emit2 + 1
!  print *,"New sequence is started: "
  do iexp=1,nexp
    y = 1d0 - 10d0**(-iexp)
    call MapFSR(pborn,xi1,y,azi1,emit1,rad1,flv,phat)
    y = 1d0 - 10d0**(-iexp)
    call MapFSR(phat(:)%p,xi2,y,azi2,emit2,rad2,flv,p)
    s_ir  = 2*p(radi)%p*p(radr)%p
    s_js  = 2*p(radj)%p*p(rads)%p
!    print *,"s_ir,s_js: ",s_ir,s_js
!    print *,"emitir,emitjs: ",emitir,emitjs
!    print *,"rad1,rad2: ",rad1,rad2
! We have to calculate the common invariants for the real
! kinematics:
!    call PrintParts(phat)
!    call PrintParts(p)
!    call CheckParts(phat)
!    call CheckParts(p)
    call CalcSubInvariants(p)
! We have to calculate the scales:
    call calcmyscales(p)
! The next line is only a simple test implies \mu_R^2 = Q^2:
!    mur = sqrt(2d0*p(1)%p*p(2)%p)
! Calculating the CSirs term:
    call PickCSirs(p,phat,ptilde,Bij,Bmunuij,emitir,radi,radr,rads, &
                   CSirs,CSirsTerm)
! We calculate only that particular Cirs term which actually cancels the
! leading divergence:
    call PickCirjsCSirs(p,phat,ptilde, &
                        emitir,radi,radr,emitjs,radj,rads, &
                        CSirs,CirjsCSirsTerm)
    write(iun,fmt='(a,2(I3,1x))',advance='no') "iexp= ", &
      int(log10(s_ir)),int(log10(s_js))
    write(iun,*) ", CirjsCSirs/CSirs= ",CirjsCSirsTerm/CSirsTerm
!
  end do
!
end subroutine CirjsCSirsLimit
!
subroutine SrsLimit(iun,pborn,p,phat,ptilde,Bij,Bijkl, &
                    flv,radr,rads,Srs)
use particles
use random
use regions
use math
use FKSutils
use scales
use my_scale
use coupling
use QCDparams
implicit none
!
  integer , intent(in) :: iun
  type(mom) , dimension(:) , intent(in) :: pborn
  type(particle) , dimension(:) , intent(inout) :: p,ptilde,phat
  real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
  real(kind(1d0)) , dimension(:,:,:,:) , intent(inout) :: Bijkl
  integer , dimension(:) , intent(in) :: flv
  integer , intent(in) :: radr,rads
  type(subterms) , intent(in) :: Srs
!
  integer :: iexp
  integer :: emit1,emit2,rad1,rad2
  real(kind(1d0)) :: xi1,xi2,azi1,azi2,y1,y2
  real(kind(1d0)) :: Er,Es
  real(kind(1d0)) :: smeRR,Srsterm
!
!
  interface 
    subroutine CalcSubInvariants(p)
    use regions
    use particles
!
      type(particle) , dimension(:) , intent(in) :: p
!
    end subroutine CalcSubInvariants
!
    subroutine PickSrs(p,ptilde,Bij,Bijkl,radr,rads, &
                       Srs,Srsterm)
    use regions
    use particles
    implicit none
!
      type(particle), dimension(:), intent(in) :: p
      type(particle), dimension(:), intent(out) :: ptilde
      real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
      real(kind(1d0)) , dimension(:,:,:,:) , intent(inout) :: Bijkl
      integer , intent(in) :: radr,rads
      type(subterms) , intent(in) :: Srs
      real(kind(1d0)) , intent(out) :: Srsterm
!
    end subroutine PickSrs
!
    subroutine CalcRR(p,smeRR)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: smeRR
!
    end subroutine CalcRR
  end interface
!
!  print *,"************************************************************"
!  print *,"We will approach the limit of: "
!  print *,"radr,rads: ",radr,rads
!  call PrintSubProc(flv)
!  print *,"************************************************************"
!  read(*,*)
!
! The double soft limit is approached by applying the FSR mapping
! twice in each step, hence azimuths are ys have to be
! specified (randomly), since the mapping changes the y it is set
! to a fixed value in each iteration:
  azi1 = 2d0*pi*gen_rand()
  azi2 = 2d0*pi*gen_rand()
! The iterated mappings only work if the radiations are ordered
! in the sense that that radiation is produced first which position
! comes first:
  if (radr.lt.rads) then
    rad1 = radr
    rad2 = rads
  else
    rad1 = rads
    rad2 = radr
  end if
! For soft emissions no emitter is assigned, but to make the
! actual mappings emitters have to be defined, to do so two
! final state massless partons are selected such that they have to
! be different and cannot coincide with the radiated partons:
  do emit1=3,size(p) - 2
    if (abs(flv(emit1)).gt.qcd_nf) cycle
    if (emit1.eq.rad1) cycle
    exit
  end do
  do emit2=3,size(p) - 1
    if (emit1.eq.emit2) cycle
    if (abs(flv(emit2)).gt.qcd_nf) cycle
    if ((emit2.eq.rad1).or.(emit2.eq.rad2)) cycle
    exit
  end do
!  print *,"emit1,emit2: ",emit1,emit2
!  print *,"rad1,rad2: ",rad1,rad2
!  print *,"New sequence is started: "
  do iexp=1,nexp
! The scaling is chosen to be the same:
    y1 = 0
    y2 = 0
    xi1 = 10d0**(-iexp)
    xi2 = 10d0**(-iexp)
! Using the Born momenta to construct an intermediate PS (phat):
    call MapFSR(pborn,xi1,y1,azi1,emit1,rad1,flv,phat)
    call MapFSR(phat(:)%p,xi2,y2,azi2,emit2,rad2,flv,p)
    Er = p(radr)%p%E
    Es = p(rads)%p%E
!    print *,"Er: ",Er
!    print *,"Es: ",Es
! We have to calculate the common invariants for the real
! kinematics:
!    call PrintParts(phat)
!    call PrintParts(p)
    call CalcSubInvariants(p)
! We have to calculate the scales:
    call calcmyscales(p)
! The next line is only a simple test implies \mu_R^2 = Q^2:
!    mur = sqrt(2d0*p(1)%p*p(2)%p)
! We obtain the RR contribution for the subprocess in question:
    call CalcRR(p,smeRR)
!    print *,"smeRR: ",smeRR
! We calculate only that particular Cirs term which actually cancels the
! leading divergence:
    call PickSrs(p,ptilde,Bij,Bijkl,radr,rads,Srs,Srsterm)
    write(iun,fmt='(a,2(I3,1x))',advance='no') "iexp= ", &
      int(log10(Er)),int(log10(Es))
    write(iun,*) ", Srs/RR= ",Srsterm/smeRR
!
  end do
!
end subroutine SrsLimit
!
subroutine CirsSrsLimit(iun,pborn,p,phat,ptilde,Bij,Bijkl, &
                        flv,emitirs,radi,radr,rads,Srs)
use particles
use random
use regions
use math
use FKSutils
use scales
use my_scale
use coupling
use QCDparams
use misc
use nnlo_mom_maps
implicit none
!
  integer , intent(in) :: iun
  type(mom) , dimension(:) , intent(in) :: pborn
  type(particle) , dimension(:) , intent(inout) :: p,ptilde,phat
  real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
  real(kind(1d0)) , dimension(:,:,:,:) , intent(inout) :: Bijkl
  integer , dimension(:) , intent(in) :: flv
  integer , intent(in) :: emitirs,radi,radr,rads
  type(subterms) , intent(in) :: Srs
!
  integer :: iexp
  integer :: emit,rad1,rad2
  real(kind(1d0)) :: xi1,xi2,azi1,azi2,y
  real(kind(1d0)) :: s_ir,s_is,s_rs,s_irs
  real(kind(1d0)) :: SrsTerm,CirsSrsTerm
!
!
  interface 
    subroutine CalcSubInvariants(p)
    use regions
    use particles
!
      type(particle) , dimension(:) , intent(in) :: p
!
    end subroutine CalcSubInvariants
!
    subroutine PickSrs(p,ptilde,Bij,Bijkl,radr,rads, &
                       Srs,Srsterm)
    use regions
    use particles
    implicit none
!
      type(particle), dimension(:), intent(in) :: p
      type(particle), dimension(:), intent(out) :: ptilde
      real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
      real(kind(1d0)) , dimension(:,:,:,:) , intent(inout) :: Bijkl
      integer , intent(in) :: radr,rads
      type(subterms) , intent(in) :: Srs
      real(kind(1d0)) , intent(out) :: Srsterm
!
    end subroutine PickSrs
!
    subroutine PickCirsSrs(p,ptilde,emitirs,radi,radr,rads, &
                           Srs,CirsSrsterm)
    use regions
    use particles
    implicit none
!
      type(particle), dimension(:), intent(in) :: p
      type(particle), dimension(:), intent(out) :: ptilde
      integer , intent(in) :: emitirs,radi,radr,rads
      type(subterms) , intent(in) :: Srs
      real(kind(1d0)) , intent(out) :: CirsSrsterm
!
    end subroutine PickCirsSrs
  end interface
!
!  print *,"************************************************************"
!  print *,"We will approach the limit of: "
!  print *,"emitirs,radi,radr,rads: ",emitirs,radi,radr,rads
!  call PrintSubProc(flv)
!  print *,"************************************************************"
!  read(*,*)
! For the CirsSrs counter term a triply collinear limit
! has to be approached, to do so the azimuths and energies
! are selected randomly:
  xi1  = gen_rand()
  xi2  = gen_rand()
  azi1 = 2d0*pi*gen_rand()
  azi2 = 2d0*pi*gen_rand()
! Determine the position for the emitted partons:
  emit = min(radi,radr,rads)
  if (emit.eq.radi) then
    rad1 = min(radr,rads)
    rad2 = max(radr,rads)
  elseif (emit.eq.radr) then
    rad1 = min(radi,rads)
    rad2 = max(radi,rads)
  elseif (emit.eq.rads) then
    rad1 = min(radi,radr)
    rad2 = max(radi,radr)
  end if
!  print *,"New sequence is started: "
  do iexp=1,nexp
    y = 1d0 - 10d0**(-iexp)
    call MapFSR(pborn,xi1,y,azi1,emit,rad1,flv,phat)
    y = 1d0 - 10d0**(-iexp)
    call MapFSR(phat(:)%p,xi2,y,azi2,emit,rad2,flv,p)
    s_ir  = 2*p(radi)%p*p(radr)%p
    s_is  = 2*p(radi)%p*p(rads)%p
    s_rs  = 2*p(radr)%p*p(rads)%p
    s_irs = (p(radi)%p + p(radr)%p + p(rads)%p)**2
!    print *,"s_ir,s_is,s_rs,s_irs: ",s_ir,s_is,s_rs,s_irs
!    print *,"emit: ",emit
!    print *,"rad1,rad2: ",rad1,rad2
! We have to calculate the common invariants for the real
! kinematics:
!    call PrintParts(phat)
!    call PrintParts(p)
!    call CheckParts(phat)
!    call CheckParts(p)
    call CalcSubInvariants(p)
! We have to calculate the scales:
    call calcmyscales(p)
! The next line is only a simple test implies \mu_R^2 = Q^2:
!    mur = sqrt(2d0*p(1)%p*p(2)%p)
! Calculating the CSirs term:
    call PickSrs(p,ptilde,Bij,Bijkl,radr,rads, &
                 Srs,SrsTerm)
! We calculate only that particular Cirs term which actually cancels the
! leading divergence:
    call PickCirsSrs(p,ptilde,emitirs,radi,radr,rads, &
                     Srs,CirsSrsTerm)
    write(iun,fmt='(a,I3,1x)',advance='no') "iexp= ", &
      int(log10(s_irs))
    write(iun,*) ", CirsSrs/Srs= ",CirsSrsTerm/SrsTerm
!
  end do
!
end subroutine CirsSrsLimit
!
subroutine CSirsSrsLimit(iun,pborn,p,phat,ptilde,Bij,Bijkl, &
                         flv,emitir,radi,radr,rads,Srs)
use particles
use random
use regions
use math
use FKSutils
use scales
use my_scale
use coupling
use QCDparams
use misc
use nnlo_mom_maps
implicit none
!
  integer , intent(in) :: iun
  type(mom) , dimension(:) , intent(in) :: pborn
  type(particle) , dimension(:) , intent(inout) :: p,ptilde,phat
  real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
  real(kind(1d0)) , dimension(:,:,:,:) , intent(inout) :: Bijkl
  integer , dimension(:) , intent(in) :: flv
  integer , intent(in) :: emitir,radi,radr,rads
  type(subterms) , intent(in) :: Srs
!
  logical :: soft_first
  integer :: iexp
  integer :: emit1,emit2,rad1,rad2
  real(kind(1d0)) :: xi1,xi2,azi1,azi2,y1,y2
  real(kind(1d0)) :: s_ir,Es
  real(kind(1d0)) :: SrsTerm,CSirsSrsTerm
!
!
  interface 
    subroutine CalcSubInvariants(p)
    use regions
    use particles
!
      type(particle) , dimension(:) , intent(in) :: p
!
    end subroutine CalcSubInvariants
!
    subroutine PickSrs(p,ptilde,Bij,Bijkl,radr,rads, &
                       Srs,Srsterm)
    use regions
    use particles
    implicit none
!
      type(particle), dimension(:), intent(in) :: p
      type(particle), dimension(:), intent(out) :: ptilde
      real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
      real(kind(1d0)) , dimension(:,:,:,:) , intent(inout) :: Bijkl
      integer , intent(in) :: radr,rads
      type(subterms) , intent(in) :: Srs
      real(kind(1d0)) , intent(out) :: Srsterm
!
    end subroutine PickSrs
!
    subroutine PickCSirsSrs(p,ptilde,Bij,emitir,radi,radr,rads, &
                            Srs,CSirsSrsterm)
    use regions
    use particles
    implicit none
!
      type(particle), dimension(:), intent(in) :: p
      type(particle), dimension(:), intent(out) :: ptilde
      real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
      integer , intent(in) :: emitir,radi,radr,rads
      type(subterms) , intent(in) :: Srs
      real(kind(1d0)) , intent(out) :: CSirsSrsterm
!
    end subroutine PickCSirsSrs
  end interface
!
!  print *,"************************************************************"
!  print *,"We will approach the limit of: "
!  print *,"emitir,radi,radr,rads: ",emitir,radi,radr,rads
!  call PrintSubProc(flv)
!  print *,"************************************************************"
!  read(*,*)
! For the CSirsSrs counter term a collinear-soft limit
! has to be approached.
  azi1 = 2d0*pi*gen_rand()
  azi2 = 2d0*pi*gen_rand()
! Having a collinear and a soft emission we have to determine
! which comes first to correctly construct the phase space:
  soft_first = .false.
! Soft radiation comes first:
  if (max(radi,radr).gt.rads) then
    soft_first = .true.
! For the soft the enclosed angle is fixed such that it is perpendicular to
! the emitter:
    y1  = 0d0
! For the collinear pair the energy is randomly reinstated:
    xi2 = 0.5d0
    rad1 = rads
    rad2 = max(radi,radr)
    do emit1=3,size(pborn)
      if (emit1.eq.emitir) cycle
      if (abs(flv(emit1)).gt.qcd_nf) cycle
      exit
    end do
! If the soft radiation happens first position of the emitter
! for the collinear pair can change, since if the soft parton
! is inserted before the emitter of the collinear pair the
! position of that is got shifted:
    if (rads.le.emitir) then
      emit2 = emitir + 1
    else 
      emit2 = emitir
    end if
! Collinear radiation comes first:
  else
! For the soft the enclosed angle is fixed:
    y2  = 0d0
    xi1 = 0.5d0
    rad1 = max(radi,radr)
    rad2 = rads
    emit1 = emitir
    do emit2=3,size(pborn)
      if (emit2.eq.radi) cycle
      if (emit2.eq.radr) cycle
      if (abs(flv(emit2)).gt.qcd_nf) cycle
      exit
    end do
  end if
!  print *,"New sequence is started: "
  do iexp=1,nexp
! If soft radiation comes first xi1 and y2 has to be scaled:
    if (soft_first) then
      xi1 = 10d0**(-iexp/2d0)
      xi1 = (xi1*(1 - 2*sqrt(tiny_y)) + sqrt(tiny_y) - tiny_xi) &
          / (1 - 2*tiny_xi)
      y2  = 1d0 - 10d0**(-iexp)
    else
      y1  = 1d0 - 10d0**(-iexp)
      xi2 = 10d0**(-iexp/2d0)
      xi2 = (xi2*(1 - 2*sqrt(tiny_y)) + sqrt(tiny_y) - tiny_xi) &
          / (1 - 2*tiny_xi)
    end if
    call MapFSR(pborn,xi1,y1,azi1,emit1,rad1,flv,phat)
    call MapFSR(phat(:)%p,xi2,y2,azi2,emit2,rad2,flv,p)
    s_ir = 2*p(radi)%p*p(radr)%p
    Es   = p(rads)%p%E
!    print *,"s_ir,Es: ",s_ir,Es
!    print *,"emit1,emit2: ",emit1,emit2
!    print *,"rad1,rad2: ",rad1,rad2
! We have to calculate the common invariants for the real
! kinematics:
!    call PrintParts(phat)
!    call PrintParts(p)
!    call CheckParts(phat)
!    call CheckParts(p)
    call CalcSubInvariants(p)
! We have to calculate the scales:
    call calcmyscales(p)
! The next line is only a simple test implies \mu_R^2 = Q^2:
!    mur = sqrt(2d0*p(1)%p*p(2)%p)
! Calculating the CSirs term:
    call PickSrs(p,ptilde,Bij,Bijkl,radr,rads, &
                 Srs,SrsTerm)
! We calculate only that particular Cirs term which actually cancels the
! leading divergence:
    call PickCSirsSrs(p,ptilde,Bij,emitir,radi,radr,rads, &
                      Srs,CSirsSrsTerm)
    write(iun,fmt='(a,2(I3,1x))',advance='no') "iexp= ", &
      int(log10(s_ir)),int(log10(Es))
    write(iun,*) ", CSirsSrs/Srs= ",CSirsSrsTerm/SrsTerm
!
  end do
!
end subroutine CSirsSrsLimit
!
subroutine CirjsSrsLimit(iun,pborn,p,phat,ptilde,Bij,Bijkl, &
                         flv,emitir,radi,radr,emitjs,radj,rads,Srs)
use particles
use random
use regions
use math
use FKSutils
use scales
use my_scale
use coupling
use QCDparams
use misc
use nnlo_mom_maps
implicit none
!
  integer , intent(in) :: iun
  type(mom) , dimension(:) , intent(in) :: pborn
  type(particle) , dimension(:) , intent(inout) :: p,ptilde,phat
  real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
  real(kind(1d0)) , dimension(:,:,:,:) , intent(inout) :: Bijkl
  integer , dimension(:) , intent(in) :: flv
  integer , intent(in) :: emitir,emitjs,radi,radr,radj,rads
  type(subterms) , intent(in) :: Srs
!
  integer :: iexp
  integer :: emit1,emit2,rad1,rad2
  real(kind(1d0)) :: xi1,xi2,azi1,azi2,y
  real(kind(1d0)) :: s_ir,s_js
  real(kind(1d0)) :: SrsTerm,CirjsSrsTerm
!
!
  interface 
    subroutine CalcSubInvariants(p)
    use regions
    use particles
!
      type(particle) , dimension(:) , intent(in) :: p
!
    end subroutine CalcSubInvariants
!
    subroutine PickSrs(p,ptilde,Bij,Bijkl,radr,rads, &
                       Srs,Srsterm)
    use regions
    use particles
    implicit none
!
      type(particle), dimension(:), intent(in) :: p
      type(particle), dimension(:), intent(out) :: ptilde
      real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
      real(kind(1d0)) , dimension(:,:,:,:) , intent(inout) :: Bijkl
      integer , intent(in) :: radr,rads
      type(subterms) , intent(in) :: Srs
      real(kind(1d0)) , intent(out) :: Srsterm
!
    end subroutine PickSrs
!
    subroutine PickCirjsSrs(p,ptilde, &
                            emitir,radi,radr, &
                            emitjs,radj,rads, &
                            Srs,CirjsSrsterm)
    use regions
    use particles
    implicit none
!
      type(particle), dimension(:), intent(in) :: p
      type(particle), dimension(:), intent(out) :: ptilde
      integer , intent(in) :: emitir,emitjs,radi,radr,radj,rads
      type(subterms) , intent(in) :: Srs
      real(kind(1d0)) , intent(out) :: CirjsSrsterm
!
    end subroutine PickCirjsSrs
  end interface
!
!  print *,"************************************************************"
!  print *,"We will approach the limit of: "
!  print *,"emitir,emitjs,radi,radr,radj,rads: ",emitir,radi,radr,radj,rads
!  call PrintSubProc(flv)
!  print *,"************************************************************"
!  read(*,*)
! For the CirjsSrs counter term a double collinear limit
! has to be approached, to do so the azimuths and energies
! are selected randomly:
  xi1  = gen_rand()
  xi2  = gen_rand()
  azi1 = 2d0*pi*gen_rand()
  azi2 = 2d0*pi*gen_rand()
! In the case of the double collinear limit we have two
! emitters which are ordered, that is emiti < emitj, but
! not for the radiated ones. To build up phase space firstly 
! that branching has to be done where the radiated parton 
! comes first.
! Starting with the hypothesis that the first radiation is
! indeed coming from the first pair, the true radiated parton
! is that one for which the position is higher:
  rad1 = max(radi,radr)
  rad2 = max(radj,rads)
! The hypothesis was correct the first branching happens before
! the second one:
  if (rad1.lt.rad2) then
    emit1 = emitir
    emit2 = emitjs
! The ordering of the branchings is just the other way around:
  else
    rad1 = max(radj,rads)
    rad2 = max(radi,radr)
    emit1 = emitjs
    emit2 = emitir
  end if
! If the first radiated parton is before the second emitter the
! position of the second emitter has to be increased by one unit
! since the emitter numbering corresponds to the underlying Born:
  if (rad1.le.emit2) emit2 = emit2 + 1
!  print *,"New sequence is started: "
  do iexp=1,nexp
    y = 1d0 - 10d0**(-iexp)
    call MapFSR(pborn,xi1,y,azi1,emit1,rad1,flv,phat)
    y = 1d0 - 10d0**(-iexp)
    call MapFSR(phat(:)%p,xi2,y,azi2,emit2,rad2,flv,p)
    s_ir  = 2*p(radi)%p*p(radr)%p
    s_js  = 2*p(radj)%p*p(rads)%p
!    print *,"s_ir,s_js: ",s_ir,s_js
!    print *,"emit1,emit2: ",emit1,emit2
!    print *,"rad1,rad2: ",rad1,rad2
! We have to calculate the common invariants for the real
! kinematics:
!    call PrintParts(phat)
!    call PrintParts(p)
!    call CheckParts(phat)
!    call CheckParts(p)
    call CalcSubInvariants(p)
! We have to calculate the scales:
    call calcmyscales(p)
! The next line is only a simple test implies \mu_R^2 = Q^2:
!    mur = sqrt(2d0*p(1)%p*p(2)%p)
! Calculating the CSirs term:
    call PickSrs(p,ptilde,Bij,Bijkl,radr,rads, &
                 Srs,SrsTerm)
! We calculate only that particular Cirs term which actually cancels the
! leading divergence:
    call PickCirjsSrs(p,ptilde,emitir,radi,radr,emitjs,radj,rads, &
                      Srs,CirjsSrsTerm)
    write(iun,fmt='(a,2(I3,1x))',advance='no') "iexp= ", &
      int(log10(s_ir)),int(log10(s_js))
    write(iun,*) ", CirjsSrs/Srs= ",CirjsSrsTerm/SrsTerm
!
  end do
!
end subroutine CirjsSrsLimit
!
subroutine CirsCSirsSrsLimit(iun,pborn,p,phat,ptilde,Bij,    &
                             flv,emitirs,radi,radr,rads,Srs)
use particles
use random
use regions
use math
use FKSutils
use scales
use my_scale
use coupling
use QCDparams
use misc
use nnlo_mom_maps
implicit none
!
  integer , intent(in) :: iun
  type(mom) , dimension(:) , intent(in) :: pborn
  type(particle) , dimension(:) , intent(inout) :: p,ptilde,phat
  real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
  integer , dimension(:) , intent(in) :: flv
  integer , intent(in) :: emitirs,radi,radr,rads
  type(subterms) , intent(in) :: Srs
!
  integer :: iexp
  integer :: emit,rad1,rad2
  real(kind(1d0)) :: xi1,xi2,azi1,azi2,y
  real(kind(1d0)) :: s_ir,s_is,s_rs,s_irs
  real(kind(1d0)) :: CSirsSrsTerm,CirsCSirsSrsTerm
!
!
  interface 
    subroutine CalcSubInvariants(p)
    use regions
    use particles
!
      type(particle) , dimension(:) , intent(in) :: p
!
    end subroutine CalcSubInvariants
!
    subroutine PickSrs(p,ptilde,Bij,Bijkl,radr,rads, &
                       Srs,Srsterm)
    use regions
    use particles
    implicit none
!
      type(particle), dimension(:), intent(in) :: p
      type(particle), dimension(:), intent(out) :: ptilde
      real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
      real(kind(1d0)) , dimension(:,:,:,:) , intent(inout) :: Bijkl
      integer , intent(in) :: radr,rads
      type(subterms) , intent(in) :: Srs
      real(kind(1d0)) , intent(out) :: Srsterm
!
    end subroutine PickSrs
!
    subroutine PickCSirsSrs(p,ptilde,Bij,emitir,radi,radr,rads, &
                            Srs,CSirsSrsterm)
    use regions
    use particles
    implicit none
!
      type(particle), dimension(:), intent(in) :: p
      type(particle), dimension(:), intent(out) :: ptilde
      real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
      integer , intent(in) :: emitir,radi,radr,rads
      type(subterms) , intent(in) :: Srs
      real(kind(1d0)) , intent(out) :: CSirsSrsterm
!
    end subroutine PickCSirsSrs
!
    subroutine PickCirsCSirsSrs(p,ptilde,emitirs,radi,radr,rads, &
                                Srs,CirsCSirsSrsterm)
    use regions
    use particles
    implicit none
!
      type(particle), dimension(:), intent(in) :: p
      type(particle), dimension(:), intent(out) :: ptilde
      integer , intent(in) :: emitirs,radi,radr,rads
      type(subterms) , intent(in) :: Srs
      real(kind(1d0)) , intent(out) :: CirsCSirsSrsterm
!
    end subroutine PickCirsCSirsSrs
  end interface
!
!  print *,"************************************************************"
!  print *,"We will approach the limit of: "
!  print *,"emitirs,radi,radr,rads: ",emitirs,radi,radr,rads
!  call PrintSubProc(flv)
!  print *,"************************************************************"
!  read(*,*)
! For the CirsCSirsSrs counter term a triply collinear limit
! has to be approached
  xi1  = gen_rand()
  xi2  = gen_rand()
  azi1 = 2d0*pi*gen_rand()
  azi2 = 2d0*pi*gen_rand()
! Determine the position for the emitted partons:
  emit = min(radi,radr,rads)
  if (emit.eq.radi) then
    rad1 = min(radr,rads)
    rad2 = max(radr,rads)
  elseif (emit.eq.radr) then
    rad1 = min(radi,rads)
    rad2 = max(radi,rads)
  elseif (emit.eq.rads) then
    rad1 = min(radi,radr)
    rad2 = max(radi,radr)
  end if
!  print *,"New sequence is started: "
  do iexp=1,nexp
    y = 1d0 - 10d0**(-iexp)
    call MapFSR(pborn,xi1,y,azi1,emit,rad1,flv,phat)
    y = 1d0 - 10d0**(-iexp)
    call MapFSR(phat(:)%p,xi2,y,azi2,emit,rad2,flv,p)
    s_ir  = 2*p(radi)%p*p(radr)%p
    s_is  = 2*p(radi)%p*p(rads)%p
    s_rs  = 2*p(radr)%p*p(rads)%p
    s_irs = (p(radi)%p + p(radr)%p + p(rads)%p)**2
!    print *,"s_ir,s_is,s_rs,s_irs: ",s_ir,s_is,s_rs,s_irs
!    print *,"emit: ",emit
!    print *,"rad1,rad2: ",rad1,rad2
! We have to calculate the common invariants for the real
! kinematics:
!    call PrintParts(phat)
!    call PrintParts(p)
!    call CheckParts(phat)
!    call CheckParts(p)
    call CalcSubInvariants(p)
! We have to calculate the scales:
    call calcmyscales(p)
! The next line is only a simple test implies \mu_R^2 = Q^2:
!    mur = sqrt(2d0*p(1)%p*p(2)%p)
! Calculating the CSirs term:
    call PickCSirsSrs(p,ptilde,Bij,emitirs,radi,radr,rads, &
                      Srs,CSirsSrsTerm)
    call PickCirsCSirsSrs(p,ptilde,emitirs,radi,radr,rads, &
                          Srs,CirsCSirsSrsTerm)
    write(iun,fmt='(a,I3,1x)',advance='no') "iexp= ", &
      int(log10(s_irs))
    write(iun,*) ", CirsCSirsSrs/CSirsSrs= ",CirsCSirsSrsTerm/CSirsSrsTerm
!
  end do
!
end subroutine CirsCSirsSrsLimit
!
subroutine CktCktrLimit(iun,pborn,p,phat,ptilde,flv, &
                        emitktr,radk,radt,radr,Ckt,Cktr)
use particles
use random
use regions
use math
use FKSutils
use scales
use my_scale
use coupling
use misc
use QCDparams
implicit none
!
  integer , intent(in) :: iun
  type(mom) , dimension(:) , intent(in) :: pborn
  type(particle) , dimension(:) , intent(inout) :: p,ptilde,phat
  integer , dimension(:) , intent(in) :: flv
  integer , intent(in) :: emitktr,radk,radt,radr
  type(subterms) , intent(in) :: Ckt,Cktr
!
  integer :: iexp
  integer :: nleg
  integer :: emit1,emit2,rad1,rad2
  integer :: rad_i,rad_r,rad_s
  real(kind(1d0)) :: xi1,xi2,azi1,azi2,y,y1,y2
  real(kind(1d0)) :: y_kt,y_ktr
  real(kind(1d0)) :: Cirsterm,Cktterm,CktCktrterm
!
!
  interface 
    subroutine CalcSubInvariants(p)
    use regions
    use particles
!
      type(particle) , dimension(:) , intent(in) :: p
!
    end subroutine CalcSubInvariants
!
    subroutine CalcR(p,smeR)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: smeR
!
    end subroutine CalcR
!
    subroutine CalcRmunu(ileg,p,Rmunu)
    use particles
    implicit none
!
      integer , intent(in) :: ileg
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(0:3,0:3) , intent(out) :: Rmunu
!
    end subroutine CalcRmunu
!
    subroutine PickCir(p,ptilde,emit,rad,Cir,CalcSMEB,CalcSMEBmunu, &
                       Cirterm)
    use regions
    use particles
    implicit none
!
      type(particle), dimension(:), intent(in) :: p
      type(particle), dimension(:), intent(out) :: ptilde
      integer , intent(in) :: emit,rad
      type(subterms) , intent(in) :: Cir
      real(kind(1d0)) , intent(out) :: Cirterm
!
      interface
        subroutine CalcSMEB(p,smeB)
        use particles
        implicit none
!
          type(particle) , dimension(:) , intent(in) :: p
          real(kind(1d0)) , intent(out) :: smeB
!
        end subroutine CalcSMEB
!
        subroutine CalcSMEBmunu(ileg,p,Bmunu)
        use particles
        implicit none
!
          integer , intent(in) :: ileg
          type(particle) , dimension(:) , intent(in) :: p
          real(kind(1d0)) , dimension(0:3,0:3) , intent(out) :: Bmunu
!
        end subroutine CalcSMEBmunu
      end interface
    end subroutine PickCir
!
    subroutine PickCirs(p,ptilde,emit,radi,radr,rads,Cirs,Cirsterm)
    use regions
    use particles
    implicit none
!
      type(particle), dimension(:), intent(in) :: p
      type(particle), dimension(:), intent(out) :: ptilde
      integer , intent(in) :: emit,radi,radr,rads
      type(subterms) , intent(in) :: Cirs
      real(kind(1d0)) , intent(out) :: Cirsterm
!
    end subroutine PickCirs
!
    subroutine PickCktCktr(p,phat,ptilde, &
                           emitktr,radk,radt,radr,rad_i,rad_r,rad_s, &
                           Cktr,CktCktrterm)
    use regions
    use particles
    implicit none
!
      type(particle), dimension(:), intent(in) :: p
      type(particle), dimension(:), intent(out) :: phat,ptilde
      integer , intent(in) :: emitktr,radk,radt,radr,rad_i,rad_r,rad_s
      type(subterms) , intent(in) :: Cktr
      real(kind(1d0)) , intent(out) :: CktCktrterm
!
    end subroutine PickCktCktr
  end interface
!
!  print *,"************************************************************"
!  print *,"We will approach the limit of: "
!  print *,"emitktr,radk,radt,radr: ",emitktr,radk,radt,radr
!  call PrintSubProc(flv)
!  print *,"************************************************************"
!  read(*,*)
! These counter terms can be tested in two different limits:
! In the doubly unresolved region these counter terms regulate
! the A1-type ones while in the singly unresolved part of the
! phase space the spurious singularities of the A2-type ones:
! Starting with the singly unresolved one:
!
! First, building a PS point which has n-1 legs:
  xi1  = gen_rand()
  y1   = gen_rand()
  azi1 = gen_rand()
  nleg = size(p)
! To generate the intermediate PS point the emitter is chosen 
! to be the first massless parton and the position for the 
! radiated one is the last position:
  do emit1=3,nleg-2
    if (abs(flv(emit1)).gt.qcd_nf) cycle
    exit
  end do
  rad1 = nleg - 1
! As the emitter of the second branching that parton's
! position is used which comes first:
  emit2 = min(radk,radt)
! The position for the radiated one is just the other one:
  rad2  = max(radk,radt)
! The energy and azimuth are picked randomly:
  xi2   = gen_rand()
  azi2  = gen_rand()
! Obtaining the parton order for the Cirs counterterm: 
  call ReorderRadirs(radk,flv(radk), &
                     radt,flv(radt), &
                     radr,flv(radr), &
                     rad_i,rad_r,rad_s)
!  print *,"New sequence is started: "
  write(iun,"(/a/)") "Singly unresolved region: "
  do iexp=1,nexp
! The scaling is chosen to be the same:
    y2 = 1d0 - 10d0**(-iexp)
    call MapFSR(pborn,xi1,y1,azi1,emit1,rad1,flv,phat)
! Using the intermediate PS point to build the RR one:
    call MapFSR(phat(:)%p,xi2,y2,azi2,emit2,rad2,flv,p)
! We have to calculate the common invariants for the real
! kinematics:
!    call PrintParts(p)
    call CalcSubInvariants(p)
    y_kt = (p(radk)%p + p(radt)%p)**2 / Q2
!    print *,"y_kt: ",y_kt
! We have to calculate the scales:
    call calcmyscales(p)
    call PickCirs(p,ptilde,emitktr,rad_i,rad_r,rad_s,Cktr,Cirsterm)
    call PickCktCktr(p,phat,ptilde, &
                     emitktr,radk,radt,radr,rad_i,rad_r,rad_s, &
                     Cktr,CktCktrterm)
!    print *,Cirsterm,CktCktrterm
    write(iun,fmt='(a,I3,1x)',advance='no') "iexp= ",int(log10(y_kt))
    write(iun,*) ", CktCktr/Cirs= ",CktCktrterm/Cirsterm
!
  end do
!
! In the second case the counterterm is checked in a doubly 
! unresolved part of the phase space:
  xi1  = gen_rand()
  xi2  = gen_rand()
  azi1 = 2d0*pi*gen_rand()
  azi2 = 2d0*pi*gen_rand()
  if (emitktr.eq.radk) then
    rad1 = min(radt,radr)
    rad2 = max(radt,radr)
  elseif (emitktr.eq.radt) then
    rad1 = min(radk,radr)
    rad2 = max(radk,radr)
  elseif (emitktr.eq.radr) then
    rad1 = min(radk,radt)
    rad2 = max(radk,radt)
  end if
!  print *,"New sequence is started: "
  write(iun,"(/a/)") "Doubly unresolved region: "
  do iexp=1,nexp
! The scaling is chosen to be the same:
    y = 1d0 - 10d0**(-iexp)
! Using the Born momenta to construct an intermediate PS (phat):
    call MapFSR(pborn,xi1,y,azi1,emitktr,rad1,flv,phat)
! Using the intermediate momenta to construct the final RR PS: 
    y = 1d0 - 10d0**(-iexp)
    call MapFSR(phat(:)%p,xi2,y,azi2,emitktr,rad2,flv,p)
! We have to calculate the common invariants for the real
! kinematics:
!    call PrintParts(phat)
!    call PrintParts(p)
    call CalcSubInvariants(p)
    y_ktr = (p(radk)%p + p(radt)%p + p(radr)%p)**2/Q2
! We have to calculate the scales:
    call calcmyscales(p)
!    print *,"y_ktr: ",y_ktr
! Calculating the corresponding Ckt term:
    call PickCir(p,phat,radk,radt,Ckt,CalcR,CalcRmunu,Cktterm)
    call PickCktCktr(p,phat,ptilde, &
                     emitktr,radk,radt,radr,rad_i,rad_r,rad_s, &
                     Cktr,CktCktrterm)
    write(iun,fmt='(a,I3,1x)',advance='no') "iexp= ",int(log10(y_ktr))
    write(iun,*) ", CktCktr/Ckt= ",CktCktrterm/Cktterm
!
  end do
!
end subroutine CktCktrLimit
!
subroutine CktCirktLimit(iun,pborn,p,phat,ptilde,flv, &
                         emitir,radi,radr,emitkt,radk,radt, &
                         Ckt,Cirkt)
use particles
use random
use regions
use math
use FKSutils
use scales
use my_scale
use coupling
use misc
use QCDparams
implicit none
!
  integer , intent(in) :: iun
  type(mom) , dimension(:) , intent(in) :: pborn
  type(particle) , dimension(:) , intent(inout) :: p,ptilde,phat
  integer , dimension(:) , intent(in) :: flv
  integer , intent(in) :: emitir,emitkt,radi,radr,radk,radt
  type(subterms) , intent(in) :: Ckt,Cirkt
!
  integer :: iexp
  integer :: nleg
  integer :: emit1,emit2,rad1,rad2
  integer :: emit_ir,emit_js,rad_i,rad_r,rad_j,rad_s
  real(kind(1d0)) :: xi1,xi2,azi1,azi2,y,y1,y2
  real(kind(1d0)) :: y_kt,y_ir
  real(kind(1d0)) :: Cirjsterm,Cktterm,CktCirktterm
!
!
  interface 
    subroutine CalcSubInvariants(p)
    use regions
    use particles
!
      type(particle) , dimension(:) , intent(in) :: p
!
    end subroutine CalcSubInvariants
!
    subroutine CalcR(p,smeR)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: smeR
!
    end subroutine CalcR
!
    subroutine CalcRmunu(ileg,p,Rmunu)
    use particles
    implicit none
!
      integer , intent(in) :: ileg
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(0:3,0:3) , intent(out) :: Rmunu
!
    end subroutine CalcRmunu
!
    subroutine PickCir(p,ptilde,emit,rad,Cir,CalcSMEB,CalcSMEBmunu, &
                       Cirterm)
    use regions
    use particles
    implicit none
!
      type(particle), dimension(:), intent(in) :: p
      type(particle), dimension(:), intent(out) :: ptilde
      integer , intent(in) :: emit,rad
      type(subterms) , intent(in) :: Cir
      real(kind(1d0)) , intent(out) :: Cirterm
!
      interface
        subroutine CalcSMEB(p,smeB)
        use particles
        implicit none
!
          type(particle) , dimension(:) , intent(in) :: p
          real(kind(1d0)) , intent(out) :: smeB
!
        end subroutine CalcSMEB
!
        subroutine CalcSMEBmunu(ileg,p,Bmunu)
        use particles
        implicit none
!
          integer , intent(in) :: ileg
          type(particle) , dimension(:) , intent(in) :: p
          real(kind(1d0)) , dimension(0:3,0:3) , intent(out) :: Bmunu
!
        end subroutine CalcSMEBmunu
      end interface
    end subroutine PickCir
!
    subroutine PickCirjs(p,ptilde,emiti,radi,radr,emitj,radj,rads, &
                         Cirjs,Cirjsterm)
    use regions
    use particles
    implicit none
!
      type(particle), dimension(:), intent(in) :: p
      type(particle), dimension(:), intent(out) :: ptilde
      integer , intent(in) :: emiti,radi,radr,emitj,radj,rads
      type(subterms) , intent(in) :: Cirjs
      real(kind(1d0)) , intent(out) :: Cirjsterm
!
    end subroutine PickCirjs
!
    subroutine PickCktCirkt(p,phat,ptilde, &
                            emitir,radi,radr,emitkt,radk,radt, &
                            rad_i,rad_r,rad_j,rad_s, &
                            Cktr,CktCirktterm)
    use regions
    use particles
    implicit none
!
      type(particle), dimension(:), intent(in) :: p
      type(particle), dimension(:), intent(out) :: phat,ptilde
      integer , intent(in) :: emitir,emitkt,radi,radr,radk,radt, &
                              rad_i,rad_r,rad_j,rad_s
      type(subterms) , intent(in) :: Cktr
      real(kind(1d0)) , intent(out) :: CktCirktterm
!
    end subroutine PickCktCirkt
  end interface
!
!  print *,"************************************************************"
!  print *,"We will approach the limit of: "
!  print *,"emitir,radi,radr,emitkt,radk,radt: ",emitir,radi,radr,emitkt,radk,radt
!  call PrintSubProc(flv)
!  print *,"************************************************************"
!  read(*,*)
! These counter terms can be tested in two different limits:
! In the doubly unresolved region these counter terms regulate
! the A1-type ones while in the singly unresolved part of the
! phase space the spurious singularities of the A2-type ones:
! Starting with the singly unresolved one:
!
! First, building a PS point which has n-1 legs:
  xi1  = gen_rand()
  y1   = gen_rand()
  azi1 = gen_rand()
  nleg = size(p)
! To generate the intermediate PS point the emitter is chosen 
! to be the first massless parton and the position for the 
! radiated one is the last position:
  do emit1=3,nleg-2
    if (abs(flv(emit1)).gt.qcd_nf) cycle
    exit
  end do
  rad1 = nleg - 1
! As the emitter of the second branching that parton's
! position is used which comes first:
  emit2 = min(radk,radt)
! The position for the radiated one is just the other one:
  rad2  = max(radk,radt)
! The energy and azimuth are picked randomly:
  xi2   = gen_rand()
  azi2  = gen_rand()
! Obtaining the parton order for the Cirjs counterterm: 
  call ReorderRadirjs(emitir,radi,radr,emitkt,radk,radt, &
                      emit_ir,rad_i,rad_r,emit_js,rad_j,rad_s)
!  write(*,"(a,4(I0,1x))") "radi,radr,radk,radt: ",radi,radr,radk,radt
!  write(*,"(a,4(I0,1x))") "rad_i,rad_r,rad_j,rad_s: ",rad_i,rad_r,rad_j,rad_s
!  print *,"New sequence is started: "
  write(iun,"(/a/)") "Singly unresolved region: "
  do iexp=1,nexp
! The scaling is chosen to be the same:
    y2 = 1d0 - 10d0**(-iexp)
    call MapFSR(pborn,xi1,y1,azi1,emit1,rad1,flv,phat)
! Using the intermediate PS point to build the RR one:
    call MapFSR(phat(:)%p,xi2,y2,azi2,emit2,rad2,flv,p)
! We have to calculate the common invariants for the real
! kinematics:
    call CalcSubInvariants(p)
    y_kt = (p(radk)%p + p(radt)%p)**2 / Q2
!    print *,"y_kt: ",y_kt
! We have to calculate the scales:
    call calcmyscales(p)
    call PickCirjs(p,ptilde,emit_ir,rad_i,rad_r,emit_js,rad_j,rad_s,Cirkt,Cirjsterm)
    call PickCktCirkt(p,phat,ptilde, &
                      emitir,radi,radr,emitkt,radk,radt, &
                      rad_i,rad_r,rad_j,rad_s, &
                      Cirkt,CktCirktterm)
    write(iun,fmt='(a,I3,1x)',advance='no') "iexp= ",int(log10(y_kt))
    write(iun,*) ", CktCirkt/Cirjs= ",CktCirktterm/Cirjsterm
!
  end do
!
! In the second case the counterterm is checked in a doubly 
! unresolved part of the phase space:
  xi1  = gen_rand()
  xi2  = gen_rand()
  azi1 = 2d0*pi*gen_rand()
  azi2 = 2d0*pi*gen_rand()
  rad1 = max(radi,radr)
  rad2 = max(radk,radt)
! The hypothesis was correct the first branching happens before
! the second one:
  if (rad1.lt.rad2) then
    emit1 = min(radi,radr)
    emit2 = min(radk,radt)
! The ordering of the branchings is just the other way around:
  else
    rad1 = max(radk,radt)
    rad2 = max(radi,radr)
    emit1 = min(radk,radt)
    emit2 = min(radi,radr)
  end if
!  write(*,"(a,4(I0,1x))") "radi,radr,radk,radt: ",radi,radr,radk,radt
!  write(*,"(a,4(I0,1x))") "radi,radr,radk,radt: ",rad1,rad2,emit1,emit2
!  print *,"New sequence is started: "
  write(iun,"(/a/)") "Doubly unresolved region: "
  do iexp=1,nexp
! The scaling is chosen to be the same:
    y = 1d0 - 10d0**(-iexp)
! Using the Born momenta to construct an intermediate PS (phat):
    call MapFSR(pborn,xi1,y,azi1,emit1,rad1,flv,phat)
! Using the intermediate momenta to construct the final RR PS: 
    y = 1d0 - 10d0**(-iexp)
    call MapFSR(phat(:)%p,xi2,y,azi2,emit2,rad2,flv,p)
! We have to calculate the common invariants for the real
! kinematics:
    call CalcSubInvariants(p)
    y_ir = (p(radi)%p + p(radr)%p)**2/Q2
    y_kt = (p(radk)%p + p(radt)%p)**2/Q2
! We have to calculate the scales:
    call calcmyscales(p)
!    print *,"y_ir,y_kt: ",y_ir,y_kt
! Calculating the corresponding Ckt term:
    call PickCir(p,phat,radk,radt,Ckt,CalcR,CalcRmunu,Cktterm)
    call PickCktCirkt(p,phat,ptilde, &
                      emitir,radi,radr,emitkt,radk,radt, &
                      rad_i,rad_r,rad_j,rad_s, &
                      Cirkt,CktCirktterm)
    write(iun,fmt='(a,2(I3,1x))',advance='no') "iexp= ", &
      int(log10(y_ir)),int(log10(y_kt))
    write(iun,*) ", CktCirkt/Ckt= ",CktCirktterm/Cktterm
!
  end do
!
end subroutine CktCirktLimit
!
subroutine CktCSktrLimit(iun,pborn,p,phat,ptilde,Bij,Bmunuij, &
                         flv,emitkt,radk,radt,radr,Ckt,CSktr)
use particles
use random
use regions
use math
use FKSutils
use scales
use my_scale
use coupling
use misc
use QCDparams
implicit none
!
  integer , intent(in) :: iun
  type(mom) , dimension(:) , intent(in) :: pborn
  type(particle) , dimension(:) , intent(inout) :: p,ptilde,phat
  real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
  real(kind(1d0)) , dimension(0:,0:,:,:) , intent(inout) :: Bmunuij
  integer , dimension(:) , intent(in) :: flv
  integer , intent(in) :: emitkt,radk,radt,radr
  type(subterms) , intent(in) :: Ckt,CSktr
!
  logical :: soft_first
  integer :: iexp
  integer :: nleg
  integer :: emit1,emit2,rad1,rad2
  integer :: rad_i,rad_r,rad_s
  real(kind(1d0)) :: xi1,xi2,azi1,azi2,y,y1,y2
  real(kind(1d0)) :: y_kt,E_r
  real(kind(1d0)) :: CSirsterm,Cktterm,CktCSktrterm
!
!
  interface 
    subroutine CalcSubInvariants(p)
    use regions
    use particles
!
      type(particle) , dimension(:) , intent(in) :: p
!
    end subroutine CalcSubInvariants
!
    subroutine CalcR(p,smeR)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: smeR
!
    end subroutine CalcR
!
    subroutine CalcRmunu(ileg,p,Rmunu)
    use particles
    implicit none
!
      integer , intent(in) :: ileg
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(0:3,0:3) , intent(out) :: Rmunu
!
    end subroutine CalcRmunu
!
    subroutine PickCir(p,ptilde,emit,rad,Cir,CalcSMEB,CalcSMEBmunu, &
                       Cirterm)
    use regions
    use particles
    implicit none
!
      type(particle), dimension(:), intent(in) :: p
      type(particle), dimension(:), intent(out) :: ptilde
      integer , intent(in) :: emit,rad
      type(subterms) , intent(in) :: Cir
      real(kind(1d0)) , intent(out) :: Cirterm
!
      interface
        subroutine CalcSMEB(p,smeB)
        use particles
        implicit none
!
          type(particle) , dimension(:) , intent(in) :: p
          real(kind(1d0)) , intent(out) :: smeB
!
        end subroutine CalcSMEB
!
        subroutine CalcSMEBmunu(ileg,p,Bmunu)
        use particles
        implicit none
!
          integer , intent(in) :: ileg
          type(particle) , dimension(:) , intent(in) :: p
          real(kind(1d0)) , dimension(0:3,0:3) , intent(out) :: Bmunu
!
        end subroutine CalcSMEBmunu
      end interface
    end subroutine PickCir
!
    subroutine PickCSirs(p,phat,ptilde,Bij,Bmunuij,emit,radi,radr,rads, &
                         CSirs,CSirsterm)
    use regions
    use particles
    implicit none
!
      type(particle), dimension(:), intent(in) :: p
      type(particle), dimension(:), intent(out) :: phat,ptilde
      real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
      real(kind(1d0)) , dimension(0:,0:,:,:) , intent(inout) :: Bmunuij
      integer , intent(in) :: emit,radi,radr,rads
      type(subterms) , intent(in) :: CSirs
      real(kind(1d0)) , intent(out) :: CSirsterm
!
    end subroutine PickCSirs
!
    subroutine PickCktCSktr(p,phat,ptilde,Bij,Bmunuij, &
                            emitkt,radk,radt,radr, &
                            CSktr,CktCSktrterm)
    use regions
    use particles
    implicit none
!
      type(particle), dimension(:), intent(in) :: p
      type(particle), dimension(:), intent(out) :: phat,ptilde
      real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
      real(kind(1d0)) , dimension(0:,0:,:,:) , intent(inout) :: Bmunuij
      integer , intent(in) :: emitkt,radk,radt,radr
      type(subterms) , intent(in) :: CSktr
      real(kind(1d0)) , intent(out) :: CktCSktrterm
!
    end subroutine PickCktCSktr
  end interface
!
!  print *,"************************************************************"
!  print *,"We will approach the limit of: "
!  print *,"emitkt,radk,radt,radr: ",emitkt,radk,radt,radr
!  call PrintSubProc(flv)
!  print *,"************************************************************"
!  read(*,*)
! These counter terms can be tested in two different limits:
! In the doubly unresolved region these counter terms regulate
! the A1-type ones while in the singly unresolved part of the
! phase space the spurious singularities of the A2-type ones:
! Starting with the singly unresolved one:
!
! First, building a PS point which has n-1 legs:
  xi1  = gen_rand()
  y1   = gen_rand()
  azi1 = gen_rand()
  nleg = size(p)
! To generate the intermediate PS point the emitter is chosen 
! to be the first massless parton and the position for the 
! radiated one is the last position:
  do emit1=3,nleg-2
    if (abs(flv(emit1)).gt.qcd_nf) cycle
    exit
  end do
  rad1 = nleg - 1
! As the emitter of the second branching that parton's
! position is used which comes first:
  emit2 = min(radk,radt)
! The position for the radiated one is just the other one:
  rad2  = max(radk,radt)
! The energy and azimuth are picked randomly:
  xi2   = gen_rand()
  azi2  = gen_rand()
! Obtaining the parton order for the Cirs counterterm: 
  call ReorderRadirs(radk,flv(radk), &
                     radt,flv(radt), &
                     radr,flv(radr), &
                     rad_i,rad_r,rad_s)
!  print *,"New sequence is started: "
  write(iun,"(/a/)") "Singly unresolved region: "
  do iexp=1,nexp
! The scaling is chosen to be the same:
    y2 = 1d0 - 10d0**(-iexp)
    call MapFSR(pborn,xi1,y1,azi1,emit1,rad1,flv,phat)
! Using the intermediate PS point to build the RR one:
    call MapFSR(phat(:)%p,xi2,y2,azi2,emit2,rad2,flv,p)
! We have to calculate the common invariants for the real
! kinematics:
!    call PrintParts(p)
    call CalcSubInvariants(p)
    y_kt = (p(radk)%p + p(radt)%p)**2 / Q2
!    print *,"y_kt: ",y_kt
! We have to calculate the scales:
    call calcmyscales(p)
    call PickCSirs(p,phat,ptilde,Bij,Bmunuij,emitkt,radk,radt,radr,CSktr,CSirsterm)
    call PickCktCSktr(p,phat,ptilde,Bij,Bmunuij, &
                      emitkt,radk,radt,radr,     &
                      CSktr,CktCSktrterm)
    write(iun,fmt='(a,I3,1x)',advance='no') "iexp= ",int(log10(y_kt))
    write(iun,*) ", CktCSktr/CSirs= ",CktCSktrterm/CSirsterm
!
  end do
!
! In the second case the counterterm is checked in a doubly 
! unresolved part of the phase space:
!
! In the case of the double soft-collinear limit azimuths are
! fixed randomly.
  azi1 = 2d0*pi*gen_rand()
  azi2 = 2d0*pi*gen_rand()
! The phase space construction depends upon the order of legs,
! if the soft leg comes first, the soft branching is considered first:
  soft_first = .false.
! Soft radiation comes first:
  if (max(radk,radt).gt.radr) then
    soft_first = .true.
! For the soft the enclosed angle is fixed such that it is perpendicular to
! the emitter:
    y1  = 1d0
! For the collinear pair the energy is randomly reinstated:
    xi2 = 0.5d0
    rad1 = radr
    rad2 = max(radk,radt)
    do emit1=3,size(pborn)
      if (emit1.eq.emitkt) cycle
      if (abs(flv(emit1)).gt.qcd_nf) cycle
      exit
    end do
! If the soft radiation happens first position of the emitter
! for the collinear pair can change, since if the soft parton
! is inserted before the emitter of the collinear pair the
! position of that is got shifted:
    if (radr.le.emitkt) then
      emit2 = emitkt + 1
    else 
      emit2 = emitkt
    end if
! Collinear radiation comes first:
  else
! For the soft the enclosed angle is fixed:
    y2  = 0d0
    xi1 = 0.5d0
    rad1 = max(radk,radt)
    rad2 = radr
    emit1 = emitkt
    do emit2=3,size(pborn)
      if (emit2.eq.radk) cycle
      if (emit2.eq.radt) cycle
      if (abs(flv(emit2)).gt.qcd_nf) cycle
      exit
    end do
  end if
!  print *,"New sequence is started: "
  write(iun,"(/a/)") "Doubly unresolved region: "
  do iexp=1,nexp
! Using the Born momenta to construct an intermediate PS (phat):
! If soft radiation comes first xi1 and y2 has to be scaled:
    if (soft_first) then
      xi1 = 10d0**(-iexp/2d0)
      xi1 = (xi1*(1 - 2*sqrt(tiny_y)) + sqrt(tiny_y) - tiny_xi) &
          / (1 - 2*tiny_xi)
      y2  = 1d0 - 10d0**(-iexp)
    else
      y1  = 1d0 - 10d0**(-iexp)
      xi2 = 10d0**(-iexp/2d0)
      xi2 = (xi2*(1 - 2*sqrt(tiny_y)) + sqrt(tiny_y) - tiny_xi) &
          / (1 - 2*tiny_xi)
    end if
    call MapFSR(pborn,xi1,y1,azi1,emit1,rad1,flv,phat)
    call MapFSR(phat(:)%p,xi2,y2,azi2,emit2,rad2,flv,p)
! We have to calculate the common invariants for the real
! kinematics:
!    call PrintParts(phat)
!    call PrintParts(p)
!    call CheckParts(phat)
!    call CheckParts(p)
    call CalcSubInvariants(p)
    y_kt = 2*p(radk)%p*p(radt)%p / Q2
    E_r  = p(radr)%p%E / sqrt(Q2)
!    print *,"y_kt, E_r: ",y_kt,E_r
! We have to calculate the scales:
    call calcmyscales(p)
! Calculating the corresponding Ckt term:
    call PickCir(p,phat,radk,radt,Ckt,CalcR,CalcRmunu,Cktterm)
    call PickCktCSktr(p,phat,ptilde,Bij,Bmunuij, &
                      emitkt,radk,radt,radr,     &
                      CSktr,CktCSktrterm)
    write(iun,fmt='(a,2(I3,1x))',advance='no') "iexp= ",int(log10(y_kt)),int(log10(E_r))
    write(iun,*) ", CktCSktr/Ckt= ",CktCSktrterm/Cktterm
!
  end do
!
end subroutine CktCSktrLimit
!
subroutine CktCirktCSktrLimit(iun,pborn,p,phat,ptilde,flv, &
                              emitir,radi,radr,emitkt,radk,radt, &
                              Ckt,CSktr)
use particles
use random
use regions
use math
use FKSutils
use scales
use my_scale
use coupling
use misc
use QCDparams
implicit none
!
  integer , intent(in) :: iun
  type(mom) , dimension(:) , intent(in) :: pborn
  type(particle) , dimension(:) , intent(inout) :: p,ptilde,phat
  integer , dimension(:) , intent(in) :: flv
  integer , intent(in) :: emitir,emitkt,radi,radr,radk,radt
  type(subterms) , intent(in) :: Ckt,CSktr
!
  integer :: iexp
  integer :: nleg
  integer :: emit1,emit2,rad1,rad2
  real(kind(1d0)) :: xi1,xi2,azi1,azi2,y1,y2
  real(kind(1d0)) :: y_kt
  real(kind(1d0)) :: CirjsCSirsterm,Cktterm,CktCirktCSktrterm
!
!
  interface 
    subroutine CalcSubInvariants(p)
    use regions
    use particles
!
      type(particle) , dimension(:) , intent(in) :: p
!
    end subroutine CalcSubInvariants
!
    subroutine CalcR(p,smeR)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: smeR
!
    end subroutine CalcR
!
    subroutine CalcRmunu(ileg,p,Rmunu)
    use particles
    implicit none
!
      integer , intent(in) :: ileg
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(0:3,0:3) , intent(out) :: Rmunu
!
    end subroutine CalcRmunu
!
    subroutine PickCir(p,ptilde,emit,rad,Cir,CalcSMEB,CalcSMEBmunu, &
                       Cirterm)
    use regions
    use particles
    implicit none
!
      type(particle), dimension(:), intent(in) :: p
      type(particle), dimension(:), intent(out) :: ptilde
      integer , intent(in) :: emit,rad
      type(subterms) , intent(in) :: Cir
      real(kind(1d0)) , intent(out) :: Cirterm
!
      interface
        subroutine CalcSMEB(p,smeB)
        use particles
        implicit none
!
          type(particle) , dimension(:) , intent(in) :: p
          real(kind(1d0)) , intent(out) :: smeB
!
        end subroutine CalcSMEB
!
        subroutine CalcSMEBmunu(ileg,p,Bmunu)
        use particles
        implicit none
!
          integer , intent(in) :: ileg
          type(particle) , dimension(:) , intent(in) :: p
          real(kind(1d0)) , dimension(0:3,0:3) , intent(out) :: Bmunu
!
        end subroutine CalcSMEBmunu
      end interface
    end subroutine PickCir
!
    subroutine PickCirjsCSirs(p,phat,ptilde, &
                              emitir,radi,radr, &
                              emitjs,radj,rads, &
                              CSirs,CirjsCSirsterm)
    use regions
    use particles
    implicit none
!
      type(particle), dimension(:), intent(in) :: p
      type(particle), dimension(:), intent(out) :: phat,ptilde
      integer , intent(in) :: emitir,emitjs,radi,radr,radj,rads
      type(subterms) , intent(in) :: CSirs
      real(kind(1d0)) , intent(out) :: CirjsCSirsterm
!
    end subroutine PickCirjsCSirs
!
    subroutine PickCktCirktCSktr(p,phat,ptilde, &
                                 emitir,radi,radr,emitkt,radk,radt, &
                                 CSktr,CktCirktCSktrterm)
    use regions
    use particles
    implicit none
!
      type(particle), dimension(:), intent(in) :: p
      type(particle), dimension(:), intent(out) :: phat,ptilde
      integer , intent(in) :: emitir,emitkt,radi,radr,radk,radt
      type(subterms) , intent(in) :: CSktr
      real(kind(1d0)) , intent(out) :: CktCirktCSktrterm
!
    end subroutine PickCktCirktCSktr
  end interface
!
!  print *,"************************************************************"
!  print *,"We will approach the limit of: "
!  print *,"emitir,radi,radr,emitkt,radk,radt: ",emitir,radi,radr,emitkt,radk,radt
!  call PrintSubProc(flv)
!  print *,"************************************************************"
!  read(*,*)
!
! First, building a PS point which has n-1 legs:
  xi1  = gen_rand()
  y1   = gen_rand()
  azi1 = gen_rand()
  nleg = size(p)
! To generate the intermediate PS point the emitter is chosen 
! to be the first massless parton and the position for the 
! radiated one is the last position:
  do emit1=3,nleg-2
    if (abs(flv(emit1)).gt.qcd_nf) cycle
    exit
  end do
  rad1 = nleg - 1
! As the emitter of the second branching that parton's
! position is used which comes first:
  emit2 = min(radk,radt)
! The position for the radiated one is just the other one:
  rad2  = max(radk,radt)
! The energy and azimuth are picked randomly:
  xi2   = gen_rand()
  azi2  = gen_rand()
!  print *,"New sequence is started: "
  write(iun,"(/a/)") "Singly unresolved region: "
  do iexp=1,nexp
    y2 = 1d0 - 10d0**(-iexp)
    call MapFSR(pborn,xi1,y1,azi1,emit1,rad1,flv,phat)
! Using the intermediate PS point to build the RR one:
    call MapFSR(phat(:)%p,xi2,y2,azi2,emit2,rad2,flv,p)
! We have to calculate the common invariants for the real
! kinematics:
    call CalcSubInvariants(p)
    y_kt = (p(radk)%p + p(radt)%p)**2 / Q2
!    print *,"y_kt: ",y_kt
! We have to calculate the scales:
    call calcmyscales(p)
    call PickCirjsCSirs(p,phat,ptilde, &
                        emitkt,radk,radt,emitir,radi,radr, &
                        CSktr,CirjsCSirsterm)
    call PickCktCirktCSktr(p,phat,ptilde, &
                           emitir,radi,radr,emitkt,radk,radt, &
                           CSktr,CktCirktCSktrterm)
    write(iun,fmt='(a,I3,1x)',advance='no') "iexp= ",int(log10(y_kt))
    write(iun,fmt='(a,G0)') ", CktCirktCSktr/CirjsCSirs= ", &
                            CktCirktCSktrterm/CirjsCSirsterm
!
  end do
!
end subroutine CktCirktCSktrLimit
!
subroutine CktCktrCSktrLimit(iun,pborn,p,phat,ptilde,flv, &
                             emitktr,radk,radt,radr,Ckt,CSktr)
use particles
use random
use regions
use math
use FKSutils
use scales
use my_scale
use coupling
use misc
use QCDparams
implicit none
!
  integer , intent(in) :: iun
  type(mom) , dimension(:) , intent(in) :: pborn
  type(particle) , dimension(:) , intent(inout) :: p,ptilde,phat
  integer , dimension(:) , intent(in) :: flv
  integer , intent(in) :: emitktr,radk,radt,radr
  type(subterms) , intent(in) :: Ckt,CSktr
!
  integer :: iexp
  integer :: nleg
  integer :: emit1,emit2,rad1,rad2
  integer :: rad_i,rad_r,rad_s
  real(kind(1d0)) :: xi1,xi2,azi1,azi2,y1,y2
  real(kind(1d0)) :: y_kt
  real(kind(1d0)) :: CktrCSktrterm,CktCktrCSktrterm
!
!
  interface 
    subroutine CalcSubInvariants(p)
    use regions
    use particles
!
      type(particle) , dimension(:) , intent(in) :: p
!
    end subroutine CalcSubInvariants
!
    subroutine PickCirsCSirs(p,phat,ptilde,emit,radi,radr,rads, &
                             CSirs,CirsCSirsterm)
    use regions
    use particles
    implicit none
!
      type(particle), dimension(:), intent(in) :: p
      type(particle), dimension(:), intent(out) :: phat,ptilde
      integer , intent(in) :: emit,radi,radr,rads
      type(subterms) , intent(in) :: CSirs
      real(kind(1d0)) , intent(out) :: CirsCSirsterm
!
    end subroutine PickCirsCSirs
!
    subroutine PickCktCktrCSktr(p,phat,ptilde,          &
                                emitktr,radk,radt,radr, &
                                CSktr,CktCktrCSktrterm)
    use regions
    use particles
    implicit none
!
      type(particle), dimension(:), intent(in) :: p
      type(particle), dimension(:), intent(out) :: phat,ptilde
      integer , intent(in) :: emitktr,radk,radt,radr
      type(subterms) , intent(in) :: CSktr
      real(kind(1d0)) , intent(out) :: CktCktrCSktrterm
!
    end subroutine PickCktCktrCSktr
  end interface
!
!  print *,"************************************************************"
!  print *,"We will approach the limit of: "
!  print *,"emitktr,radk,radt,radr: ",emitktr,radk,radt,radr
!  call PrintSubProc(flv)
!  print *,"************************************************************"
!  read(*,*)
!
! First, building a PS point which has n-1 legs:
  xi1  = gen_rand()
  y1   = gen_rand()
  azi1 = gen_rand()
  nleg = size(p)
! To generate the intermediate PS point the emitter is chosen 
! to be the first massless parton and the position for the 
! radiated one is the last position:
  do emit1=3,nleg-2
    if (abs(flv(emit1)).gt.qcd_nf) cycle
    exit
  end do
  rad1 = nleg - 1
! As the emitter of the second branching that parton's
! position is used which comes first:
  emit2 = min(radk,radt)
! The position for the radiated one is just the other one:
  rad2  = max(radk,radt)
! The energy and azimuth are picked randomly:
  xi2   = gen_rand()
  azi2  = gen_rand()
!  print *,"New sequence is started: "
  write(iun,"(/a/)") "Singly unresolved region: "
  do iexp=1,nexp
    y2 = 1d0 - 10d0**(-iexp)
    call MapFSR(pborn,xi1,y1,azi1,emit1,rad1,flv,phat)
! Using the intermediate PS point to build the RR one:
    call MapFSR(phat(:)%p,xi2,y2,azi2,emit2,rad2,flv,p)
! We have to calculate the common invariants for the real
! kinematics:
!    call PrintParts(p)
    call CalcSubInvariants(p)
    y_kt = (p(radk)%p + p(radt)%p)**2 / Q2
!    print *,"y_kt: ",y_kt
! We have to calculate the scales:
    call calcmyscales(p)
    call PickCirsCSirs(p,phat,ptilde,emitktr,radk,radt,radr, &
                       CSktr,CktrCSktrterm)
    call PickCktCktrCSktr(p,phat,ptilde,          &
                          emitktr,radk,radt,radr, &
                          CSktr,CktCktrCSktrterm)
    write(iun,fmt='(a,I3,1x)',advance='no') "iexp= ",int(log10(y_kt))
    write(iun,*) ", CktCktrCSktr/CktrCSktr= ",CktCktrCSktrterm/CktrCSktrterm
!
  end do
!
end subroutine CktCktrCSktrLimit
!
subroutine CktSktLimit(iun,pborn,p,phat,ptilde,Bij,Bijkl,flv, &
                       emitkt,radk,radt,Ckt,Skt)
use particles
use random
use regions
use math
use FKSutils
use scales
use my_scale
use coupling
use misc
use QCDparams
implicit none
!
  integer , intent(in) :: iun
  type(mom) , dimension(:) , intent(in) :: pborn
  type(particle) , dimension(:) , intent(inout) :: p,ptilde,phat
  real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
  real(kind(1d0)) , dimension(:,:,:,:) , intent(inout) :: Bijkl
  integer , dimension(:) , intent(in) :: flv
  integer , intent(in) :: emitkt,radk,radt
  type(subterms) , intent(in) :: Ckt,Skt
!
  integer :: iexp
  integer :: nleg
  integer :: emit1,emit2,rad1,rad2
  integer :: rad_i,rad_r,rad_s
  real(kind(1d0)) :: xi1,xi2,azi1,azi2,y1,y2
  real(kind(1d0)) :: y_kt
  real(kind(1d0)) :: Sktterm,CktSktterm
!
!
  interface 
    subroutine CalcSubInvariants(p)
    use regions
    use particles
!
      type(particle) , dimension(:) , intent(in) :: p
!
    end subroutine CalcSubInvariants
!
    subroutine PickSrs(p,ptilde,Bij,Bijkl,radr,rads, &
                       Srs,Srsterm)
    use regions
    use particles
    implicit none
!
      type(particle), dimension(:), intent(in) :: p
      type(particle), dimension(:), intent(out) :: ptilde
      real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
      real(kind(1d0)) , dimension(:,:,:,:) , intent(inout) :: Bijkl
      integer , intent(in) :: radr,rads
      type(subterms) , intent(in) :: Srs
      real(kind(1d0)) , intent(out) :: Srsterm
!
    end subroutine PickSrs
!
    subroutine PickCktSkt(p,phat,ptilde,Bij, &
                          emitkt,radk,radt,  &
                          Skt,CktSktterm)
    use regions
    use particles
    implicit none
!
      type(particle), dimension(:), intent(in) :: p
      type(particle), dimension(:), intent(out) :: phat,ptilde
      real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
      integer , intent(in) :: emitkt,radk,radt
      type(subterms) , intent(in) :: Skt
      real(kind(1d0)) , intent(out) :: CktSktterm
!
    end subroutine PickCktSkt
  end interface
!
!  print *,"************************************************************"
!  print *,"We will approach the limit of: "
!  print *,"emitktr,radk,radt,radr: ",emitktr,radk,radt,radr
!  call PrintSubProc(flv)
!  print *,"************************************************************"
!  read(*,*)
!
! First, building a PS point which has n-1 legs:
  xi1  = gen_rand()
  y1   = gen_rand()
  azi1 = gen_rand()
  nleg = size(p)
! To generate the intermediate PS point the emitter is chosen 
! to be the first massless parton and the position for the 
! radiated one is the last position:
  do emit1=3,nleg-2
    if (abs(flv(emit1)).gt.qcd_nf) cycle
    exit
  end do
  rad1 = nleg - 1
! As the emitter of the second branching that parton's
! position is used which comes first:
  emit2 = min(radk,radt)
! The position for the radiated one is just the other one:
  rad2  = max(radk,radt)
! The energy and azimuth are picked randomly:
  xi2   = gen_rand()
  azi2  = gen_rand()
!  print *,"New sequence is started: "
  write(iun,"(/a/)") "Singly unresolved region: "
  do iexp=1,nexp
    y2 = 1d0 - 10d0**(-iexp)
    call MapFSR(pborn,xi1,y1,azi1,emit1,rad1,flv,phat)
! Using the intermediate PS point to build the RR one:
    call MapFSR(phat(:)%p,xi2,y2,azi2,emit2,rad2,flv,p)
! We have to calculate the common invariants for the real
! kinematics:
!    call PrintParts(p)
    call CalcSubInvariants(p)
    y_kt = (p(radk)%p + p(radt)%p)**2 / Q2
!    print *,"y_kt: ",y_kt
! We have to calculate the scales:
    call calcmyscales(p)
    call PickSrs(p,ptilde,Bij,Bijkl,radk,radt,Skt,Sktterm)
    call PickCktSkt(p,phat,ptilde,Bij, &
                    emitkt,radk,radt,  &
                    Skt,CktSktterm)
    write(iun,fmt='(a,I3,1x)',advance='no') "iexp= ",int(log10(y_kt))
    write(iun,*) ", CktSkt/Skt= ",CktSktterm/Sktterm
!
  end do
!
end subroutine CktSktLimit
!
subroutine CktCrktSktLimit(iun,pborn,p,phat,ptilde,flv, &
                           emitrkt,emitkt,radr,radk,radt,Ckt,Skt)
use particles
use random
use regions
use math
use FKSutils
use scales
use my_scale
use coupling
use misc
use QCDparams
implicit none
!
  integer , intent(in) :: iun
  type(mom) , dimension(:) , intent(in) :: pborn
  type(particle) , dimension(:) , intent(inout) :: p,ptilde,phat
  integer , dimension(:) , intent(in) :: flv
  integer , intent(in) :: emitrkt,emitkt,radr,radk,radt
  type(subterms) , intent(in) :: Ckt,Skt
!
  integer :: iexp
  integer :: nleg
  integer :: emit1,emit2,rad1,rad2
  integer :: rad_i,rad_r,rad_s
  real(kind(1d0)) :: xi1,xi2,azi1,azi2,y1,y2
  real(kind(1d0)) :: y_kt
  real(kind(1d0)) :: CrktSktterm,CktCrktSktterm
!
!
  interface 
    subroutine CalcSubInvariants(p)
    use regions
    use particles
!
      type(particle) , dimension(:) , intent(in) :: p
!
    end subroutine CalcSubInvariants
!
    subroutine PickCirsSrs(p,ptilde,emitirs,radi,radr,rads, &
                           Srs,CirsSrsterm)
    use regions
    use particles
    implicit none
!
      type(particle), dimension(:), intent(in) :: p
      type(particle), dimension(:), intent(out) :: ptilde
      integer , intent(in) :: emitirs,radi,radr,rads
      type(subterms) , intent(in) :: Srs
      real(kind(1d0)) , intent(out) :: CirsSrsterm
!
    end subroutine PickCirsSrs
!
    subroutine PickCktCrktSkt(p,phat,ptilde,      &
                              radr,radk,radt,     &
                              Skt,CktCrktSktterm)
    use regions
    use particles
    implicit none
!
      type(particle), dimension(:), intent(in) :: p
      type(particle), dimension(:), intent(out) :: phat,ptilde
      integer , intent(in) :: radr,radk,radt
      type(subterms) , intent(in) :: Skt
      real(kind(1d0)) , intent(out) :: CktCrktSktterm
!
    end subroutine PickCktCrktSkt
  end interface
!
!  print *,"************************************************************"
!  print *,"We will approach the limit of: "
!  print *,"emitktr,radk,radt,radr: ",emitktr,radk,radt,radr
!  call PrintSubProc(flv)
!  print *,"************************************************************"
!  read(*,*)
!
! First, building a PS point which has n-1 legs:
  xi1  = gen_rand()
  y1   = gen_rand()
  azi1 = gen_rand()
  nleg = size(p)
! To generate the intermediate PS point the emitter is chosen 
! to be the first massless parton and the position for the 
! radiated one is the last position:
  do emit1=3,nleg-2
    if (abs(flv(emit1)).gt.qcd_nf) cycle
    exit
  end do
  rad1 = nleg - 1
! As the emitter of the second branching that parton's
! position is used which comes first:
  emit2 = min(radk,radt)
! The position for the radiated one is just the other one:
  rad2  = max(radk,radt)
! The energy and azimuth are picked randomly:
  xi2   = gen_rand()
  azi2  = gen_rand()
!  print *,"New sequence is started: "
  write(iun,"(/a/)") "Singly unresolved region: "
  do iexp=1,nexp
    y2 = 1d0 - 10d0**(-iexp)
    call MapFSR(pborn,xi1,y1,azi1,emit1,rad1,flv,phat)
! Using the intermediate PS point to build the RR one:
    call MapFSR(phat(:)%p,xi2,y2,azi2,emit2,rad2,flv,p)
! We have to calculate the common invariants for the real
! kinematics:
!    call PrintParts(p)
    call CalcSubInvariants(p)
    y_kt = (p(radk)%p + p(radt)%p)**2 / Q2
!    print *,"y_kt: ",y_kt
! We have to calculate the scales:
    call calcmyscales(p)
    call PickCirsSrs(p,ptilde,emitrkt,radr,radk,radt, &
                     Skt,CrktSktterm)
    call PickCktCrktSkt(p,phat,ptilde,      &
                        radr,radk,radt,     &
                        Skt,CktCrktSktterm)
    write(iun,fmt='(a,I3,1x)',advance='no') "iexp= ",int(log10(y_kt))
    write(iun,*) ", CktCrktSkt/CrktSkt= ",CktCrktSktterm/CrktSktterm
!
  end do
!
end subroutine CktCrktSktLimit
!
subroutine StCirtLimit(iun,pborn,p,phat,ptilde,Rij,flv, &
                       emitirt,radi,radr,radt,St,Cirt)
use particles
use random
use regions
use math
use FKSutils
use scales
use my_scale
use coupling
use misc
use QCDparams
implicit none
!
  integer , intent(in) :: iun
  type(mom) , dimension(:) , intent(in) :: pborn
  type(particle) , dimension(:) , intent(inout) :: p,ptilde,phat
  real(kind(1d0)) , dimension(:,:) , intent(inout) :: Rij
  integer , dimension(:) , intent(in) :: flv
  integer , intent(in) :: emitirt,radi,radr,radt
  type(subterms) , intent(in) :: St,Cirt
!
  integer :: iexp
  integer :: nleg
  integer :: emit1,emit2,rad1,rad2
  integer :: rad_i,rad_r,rad_s
  real(kind(1d0)) :: xi1,xi2,azi1,azi2,y,y1,y2
  real(kind(1d0)) :: E_t,y_irt
  real(kind(1d0)) :: Cirsterm,Stterm,StCirtterm
!
!
  interface 
    subroutine CalcSubInvariants(p)
    use regions
    use particles
!
      type(particle) , dimension(:) , intent(in) :: p
!
    end subroutine CalcSubInvariants
!
    subroutine CalcRij(p,Rij)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Rij
!
    end subroutine CalcRij
!
    subroutine PickSr(p,ptilde,Bij,radr,Sr,CalcSMEBij,Srterm)
    use regions
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      type(particle) , dimension(:) , intent(out) :: ptilde
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Bij
      integer , intent(in) :: radr
      type(subterms) , intent(in) :: Sr
      real(kind(1d0)) , intent(out) :: Srterm
!
      interface
        subroutine CalcSMEBij(p,Bij)
        use particles
        implicit none
!
          type(particle) , dimension(:) , intent(in) :: p
          real(kind(1d0)) , dimension(:,:) , intent(out) :: Bij
!
        end subroutine CalcSMEBij
      end interface
!
    end subroutine PickSr
!
    subroutine PickCirs(p,ptilde,emit,radi,radr,rads,Cirs,Cirsterm)
    use regions
    use particles
    implicit none
!
      type(particle), dimension(:), intent(in) :: p
      type(particle), dimension(:), intent(out) :: ptilde
      integer , intent(in) :: emit,radi,radr,rads
      type(subterms) , intent(in) :: Cirs
      real(kind(1d0)) , intent(out) :: Cirsterm
!
    end subroutine PickCirs
!
    subroutine PickStCirt(p,phat,ptilde, &
                          emitirt,radi,radr,radt,rad_i,rad_r,rad_s, &
                          Cirt,StCirtterm)
    use regions
    use particles
    implicit none
!
      type(particle), dimension(:), intent(in) :: p
      type(particle), dimension(:), intent(out) :: phat,ptilde
      integer , intent(in) :: emitirt,radi,radr,radt,rad_i,rad_r,rad_s
      type(subterms) , intent(in) :: Cirt
      real(kind(1d0)) , intent(out) :: StCirtterm
!
    end subroutine PickStCirt
  end interface
!
!  print *,"************************************************************"
!  print *,"We will approach the limit of: "
!  print *,"emitirt,radi,radr,radt: ",emitirt,radi,radr,radt
!  call PrintSubProc(flv)
!  print *,"************************************************************"
!  read(*,*)
! These counter terms can be tested in two different limits:
! In the doubly unresolved region these counter terms regulate
! the A1-type ones while in the singly unresolved part of the
! phase space the spurious singularities of the A2-type ones:
! Starting with the singly unresolved one:
!
! First, building a PS point which has n-1 legs:
  xi1  = gen_rand()
  y1   = gen_rand()
  azi1 = gen_rand()
  nleg = size(p)
! To generate the intermediate PS point the emitter is chosen 
! to be the first massless parton and the position for the 
! radiated one is the last position:
  do emit1=3,nleg-2
    if (abs(flv(emit1)).gt.qcd_nf) cycle
    exit
  end do
  rad1 = nleg - 1
! Since the second emission only produces a soft parton its
! emitter can be freely chosen to be the same as for emitter 1:
  emit2 = emit1
  rad2  = radt
! The azimuth are picked randomly:
  azi2  = gen_rand()
! while the angle is chosen to be perpendicular to emitter:
  y2    = 0d0
! Obtaining the parton order for the Cirs counterterm: 
  call ReorderRadirs(radi,flv(radi), &
                     radr,flv(radr), &
                     radt,flv(radt), &
                     rad_i,rad_r,rad_s)
!  write(*,"(a,3(I0,1x))") "radi,radr,radt: ",radi,radr,radt
!  write(*,"(a,3(I0,1x))") "radiID,radrID,radtID: ",flv(radi),flv(radr),flv(radt)
!  write(*,"(a,3(I0,1x))") "radi,radr,rads: ",rad_i,rad_r,rad_s
!  write(*,"(a,3(I0,1x))") "radiID,radrID,radsID: ",flv(rad_i),flv(rad_r),flv(rad_s)
!  print *,"New sequence is started: "
  write(iun,"(/a/)") "Singly unresolved region: "
  do iexp=1,nexp
    xi2 = 10d0**(-iexp)
    call MapFSR(pborn,xi1,y1,azi1,emit1,rad1,flv,phat)
! Using the intermediate PS point to build the RR one:
    call MapFSR(phat(:)%p,xi2,y2,azi2,emit2,rad2,flv,p)
! We have to calculate the common invariants for the real
! kinematics:
!    call PrintParts(p)
    call CalcSubInvariants(p)
    E_t = p(radt)%p%E / sqrt(Q2)
!    print *,"E_t: ",E_t
! We have to calculate the scales:
    call calcmyscales(p)
    call PickCirs(p,ptilde,emitirt,rad_i,rad_r,rad_s,Cirt,Cirsterm)
    call PickStCirt(p,phat,ptilde, &
                    emitirt,radi,radr,radt,rad_i,rad_r,rad_s, &
                    Cirt,StCirtterm)
!    print *,Cirsterm,StCirtterm
    write(iun,fmt='(a,I3,1x)',advance='no') "iexp= ",int(log10(E_t))
    write(iun,*) ", StCirt/Cirs= ",StCirtterm/Cirsterm
!
  end do
!
! In the second case the counterterm is checked in a doubly 
! unresolved part of the phase space:
  xi1  = gen_rand()
  xi2  = gen_rand()
  azi1 = 2d0*pi*gen_rand()
  azi2 = 2d0*pi*gen_rand()
  if (emitirt.eq.radi) then
    rad1 = min(radr,radt)
    rad2 = max(radr,radt)
  elseif (emitirt.eq.radr) then
    rad1 = min(radi,radt)
    rad2 = max(radi,radt)
  elseif (emitirt.eq.radt) then
    rad1 = min(radi,radr)
    rad2 = max(radi,radr)
  end if
!  print *,"New sequence is started: "
  write(iun,"(/a/)") "Doubly unresolved region: "
  do iexp=1,nexp
! The scaling is chosen to be the same:
    y = 1d0 - 10d0**(-iexp)
! Using the Born momenta to construct an intermediate PS (phat):
    call MapFSR(pborn,xi1,y,azi1,emitirt,rad1,flv,phat)
! Using the intermediate momenta to construct the final RR PS: 
    y = 1d0 - 10d0**(-iexp)
    call MapFSR(phat(:)%p,xi2,y,azi2,emitirt,rad2,flv,p)
! We have to calculate the common invariants for the real
! kinematics:
!    call PrintParts(phat)
!    call PrintParts(p)
    call CalcSubInvariants(p)
    y_irt = (p(radi)%p + p(radr)%p + p(radt)%p)**2/Q2
! We have to calculate the scales:
    call calcmyscales(p)
!    print *,"y_irt: ",y_irt
! Calculating the corresponding Ckt term:
    call PickSr(p,phat,Rij,radt,St,CalcRij,Stterm)
    call PickStCirt(p,phat,ptilde, &
                    emitirt,radi,radr,radt,rad_i,rad_r,rad_s, &
                    Cirt,StCirtterm)
    write(iun,fmt='(a,I3,1x)',advance='no') "iexp= ",int(log10(y_irt))
    write(iun,*) ", StCirt/St= ",StCirtterm/Stterm
!
  end do
!
end subroutine StCirtLimit
!
subroutine StCSirtLimit(iun,pborn,p,phat,ptilde,Bij,Rij,Bmunuij, &
                        flv,emitir,radi,radr,radt,St,CSirt)
use particles
use random
use regions
use math
use FKSutils
use scales
use my_scale
use coupling
use misc
use QCDparams
implicit none
!
  integer , intent(in) :: iun
  type(mom) , dimension(:) , intent(in) :: pborn
  type(particle) , dimension(:) , intent(inout) :: p,ptilde,phat
  real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij,Rij
  real(kind(1d0)) , dimension(0:,0:,:,:) , intent(inout) :: Bmunuij
  integer , dimension(:) , intent(in) :: flv
  integer , intent(in) :: emitir,radi,radr,radt
  type(subterms) , intent(in) :: St,CSirt
!
  logical :: soft_first
  integer :: iexp
  integer :: nleg
  integer :: emit1,emit2,rad1,rad2
  real(kind(1d0)) :: xi1,xi2,azi1,azi2,y,y1,y2
  real(kind(1d0)) :: y_ir,E_t
  real(kind(1d0)) :: CSirtterm,Stterm,StCSirtterm
!
!
  interface 
    subroutine CalcSubInvariants(p)
    use regions
    use particles
!
      type(particle) , dimension(:) , intent(in) :: p
!
    end subroutine CalcSubInvariants
!
    subroutine CalcRij(p,Rij)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Rij
!
    end subroutine CalcRij
!
    subroutine PickSr(p,ptilde,Bij,radr,Sr,CalcSMEBij,Srterm)
    use regions
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      type(particle) , dimension(:) , intent(out) :: ptilde
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Bij
      integer , intent(in) :: radr
      type(subterms) , intent(in) :: Sr
      real(kind(1d0)) , intent(out) :: Srterm
!
      interface
        subroutine CalcSMEBij(p,Bij)
        use particles
        implicit none
!
          type(particle) , dimension(:) , intent(in) :: p
          real(kind(1d0)) , dimension(:,:) , intent(out) :: Bij
!
        end subroutine CalcSMEBij
      end interface
!
    end subroutine PickSr
!
    subroutine PickCSirs(p,phat,ptilde,Bij,Bmunuij,emit,radi,radr,rads, &
                         CSirs,CSirsterm)
    use regions
    use particles
    implicit none
!
      type(particle), dimension(:), intent(in) :: p
      type(particle), dimension(:), intent(out) :: phat,ptilde
      real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
      real(kind(1d0)) , dimension(0:,0:,:,:) , intent(inout) :: Bmunuij
      integer , intent(in) :: emit,radi,radr,rads
      type(subterms) , intent(in) :: CSirs
      real(kind(1d0)) , intent(out) :: CSirsterm
!
    end subroutine PickCSirs
!
    subroutine PickStCSirt(p,phat,ptilde,Bij,Bmunuij, &
                           emitir,radi,radr,radt, &
                           CSirt,StCSirtterm)
    use regions
    use particles
    implicit none
!
      type(particle), dimension(:), intent(in) :: p
      type(particle), dimension(:), intent(out) :: phat,ptilde
      real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
      real(kind(1d0)) , dimension(0:,0:,:,:) , intent(inout) :: Bmunuij
      integer , intent(in) :: emitir,radi,radr,radt
      type(subterms) , intent(in) :: CSirt
      real(kind(1d0)) , intent(out) :: StCSirtterm
!
    end subroutine PickStCSirt
  end interface
!
!  print *,"************************************************************"
!  print *,"We will approach the limit of: "
!  print *,"emitir,radi,radr,radt: ",emitir,radi,radr,radt
!  call PrintSubProc(flv)
!  print *,"************************************************************"
!  read(*,*)
! These counter terms can be tested in two different limits:
! In the doubly unresolved region these counter terms regulate
! the A1-type ones while in the singly unresolved part of the
! phase space the spurious singularities of the A2-type ones:
! Starting with the singly unresolved one:
!
! First, building a PS point which has n-1 legs:
  xi1  = gen_rand()
  y1   = gen_rand()
  azi1 = gen_rand()
  nleg = size(p)
! To generate the intermediate PS point the emitter is chosen 
! to be the first massless parton and the position for the 
! radiated one is the last position:
  do emit1=3,nleg-2
    if (abs(flv(emit1)).gt.qcd_nf) cycle
    exit
  end do
  rad1 = nleg - 1
! Since the second emission only produces a soft parton its
! emitter can be freely chosen to be the same as for emitter 1:
  emit2 = emit1
  rad2  = radt
! The azimuth are picked randomly:
  azi2  = gen_rand()
! while the angle is chosen to be perpendicular to emitter:
  y2    = 0d0
!  print *,"New sequence is started: "
  write(iun,"(/a/)") "Singly unresolved region: "
  do iexp=1,nexp
    xi2 = 10d0**(-iexp)
    call MapFSR(pborn,xi1,y1,azi1,emit1,rad1,flv,phat)
! Using the intermediate PS point to build the RR one:
    call MapFSR(phat(:)%p,xi2,y2,azi2,emit2,rad2,flv,p)
! We have to calculate the common invariants for the real
! kinematics:
!    call PrintParts(p)
    call CalcSubInvariants(p)
    E_t = p(radt)%p%E / sqrt(Q2)
!    print *,"E_t: ",E_t
! We have to calculate the scales:
    call calcmyscales(p)
    call PickCSirs(p,phat,ptilde,Bij,Bmunuij,emitir,radi,radr,radt,CSirt,CSirtterm)
!    print *,"CSirtterm: ",CSirtterm
    call PickStCSirt(p,phat,ptilde,Bij,Bmunuij, &
                     emitir,radi,radr,radt,     &
                     CSirt,StCSirtterm)
    write(iun,fmt='(a,I3,1x)',advance='no') "iexp= ",int(log10(E_t))
    write(iun,*) ", StCSirt/CSirt= ",StCSirtterm/CSirtterm
!
  end do
!
! In the second case the counterterm is checked in a doubly 
! unresolved part of the phase space:
!
! In the case of the double soft-collinear limit azimuths are
! fixed randomly.
  azi1 = 2d0*pi*gen_rand()
  azi2 = 2d0*pi*gen_rand()
! The phase space construction depends upon the order of legs,
! if the soft leg comes first, the soft branching is considered first:
  soft_first = .false.
! Soft radiation comes first:
  if (max(radi,radr).gt.radt) then
    soft_first = .true.
! For the soft the enclosed angle is fixed such that it is perpendicular to
! the emitter:
    y1  = 1d0
! For the collinear pair the energy is randomly reinstated:
    xi2 = 0.5d0
    rad1 = radt
    rad2 = max(radi,radr)
    do emit1=3,size(pborn)
      if (emit1.eq.emitir) cycle
      if (abs(flv(emit1)).gt.qcd_nf) cycle
      exit
    end do
! If the soft radiation happens first position of the emitter
! for the collinear pair can change, since if the soft parton
! is inserted before the emitter of the collinear pair the
! position of that is got shifted:
    if (radt.le.emitir) then
      emit2 = emitir + 1
    else 
      emit2 = emitir
    end if
! Collinear radiation comes first:
  else
! For the soft the enclosed angle is fixed:
    y2  = 0d0
    xi1 = 0.5d0
    rad1 = max(radi,radr)
    rad2 = radt
    emit1 = emitir
    do emit2=3,size(pborn)
      if (emit2.eq.radi) cycle
      if (emit2.eq.radr) cycle
      if (abs(flv(emit2)).gt.qcd_nf) cycle
      exit
    end do
  end if
!  print *,"New sequence is started: "
  write(iun,"(/a/)") "Doubly unresolved region: "
  do iexp=1,nexp
! Using the Born momenta to construct an intermediate PS (phat):
! If soft radiation comes first xi1 and y2 has to be scaled:
    if (soft_first) then
      xi1 = 10d0**(-iexp/2d0)
      xi1 = (xi1*(1 - 2*sqrt(tiny_y)) + sqrt(tiny_y) - tiny_xi) &
          / (1 - 2*tiny_xi)
      y2  = 1d0 - 10d0**(-iexp)
    else
      y1  = 1d0 - 10d0**(-iexp)
      xi2 = 10d0**(-iexp/2d0)
      xi2 = (xi2*(1 - 2*sqrt(tiny_y)) + sqrt(tiny_y) - tiny_xi) &
          / (1 - 2*tiny_xi)
    end if
    call MapFSR(pborn,xi1,y1,azi1,emit1,rad1,flv,phat)
    call MapFSR(phat(:)%p,xi2,y2,azi2,emit2,rad2,flv,p)
! We have to calculate the common invariants for the real
! kinematics:
!    call PrintParts(phat)
!    call PrintParts(p)
!    call CheckParts(phat)
!    call CheckParts(p)
    call CalcSubInvariants(p)
    y_ir = 2*p(radi)%p*p(radr)%p / Q2
    E_t  = p(radt)%p%E / sqrt(Q2)
!    print *,"y_ir, E_t: ",y_ir,E_t
! We have to calculate the scales:
    call calcmyscales(p)
! Calculating the corresponding Ckt term:
    call PickSr(p,phat,Rij,radt,St,CalcRij,Stterm)
    call PickStCSirt(p,phat,ptilde,Bij,Bmunuij, &
                     emitir,radi,radr,radt,     &
                     CSirt,StCSirtterm)
    write(iun,fmt='(a,2(I3,1x))',advance='no') "iexp= ",int(log10(y_ir)),int(log10(E_t))
    write(iun,*) ", StCSirt/St= ",StCSirtterm/Stterm
!
  end do
!
end subroutine StCSirtLimit
!
subroutine StCirtCSirtLimit(iun,pborn,p,phat,ptilde,flv, &
                            emitirt,radi,radr,radt,St,CSirt)
use particles
use random
use regions
use math
use FKSutils
use scales
use my_scale
use coupling
use misc
use QCDparams
implicit none
!
  integer , intent(in) :: iun
  type(mom) , dimension(:) , intent(in) :: pborn
  type(particle) , dimension(:) , intent(inout) :: p,ptilde,phat
  integer , dimension(:) , intent(in) :: flv
  integer , intent(in) :: emitirt,radi,radr,radt
  type(subterms) , intent(in) :: St,CSirt
!
  integer :: iexp
  integer :: nleg
  integer :: emit1,emit2,rad1,rad2
  integer :: rad_i,rad_r,rad_s
  real(kind(1d0)) :: xi1,xi2,azi1,azi2,y1,y2
  real(kind(1d0)) :: E_t
  real(kind(1d0)) :: CirtCSirtterm,StCirtCSirtterm
!
!
  interface 
    subroutine CalcSubInvariants(p)
    use regions
    use particles
!
      type(particle) , dimension(:) , intent(in) :: p
!
    end subroutine CalcSubInvariants
!
    subroutine PickCirsCSirs(p,phat,ptilde,emit,radi,radr,rads, &
                             CSirs,CirsCSirsterm)
    use regions
    use particles
    implicit none
!
      type(particle), dimension(:), intent(in) :: p
      type(particle), dimension(:), intent(out) :: phat,ptilde
      integer , intent(in) :: emit,radi,radr,rads
      type(subterms) , intent(in) :: CSirs
      real(kind(1d0)) , intent(out) :: CirsCSirsterm
!
    end subroutine PickCirsCSirs
!
    subroutine PickStCirtCSirt(p,phat,ptilde,          &
                               emitirt,radi,radr,radt, &
                               CSirt,StCirtCSirtterm)
    use regions
    use particles
    implicit none
!
      type(particle), dimension(:), intent(in) :: p
      type(particle), dimension(:), intent(out) :: phat,ptilde
      integer , intent(in) :: emitirt,radi,radr,radt
      type(subterms) , intent(in) :: CSirt
      real(kind(1d0)) , intent(out) :: StCirtCSirtterm
!
    end subroutine PickStCirtCSirt
  end interface
!
!  print *,"************************************************************"
!  print *,"We will approach the limit of: "
!  print *,"emitirt,radi,radr,radt: ",emitirt,radi,radr,radt
!  call PrintSubProc(flv)
!  print *,"************************************************************"
!  read(*,*)
!
! First, building a PS point which has n-1 legs:
  xi1  = gen_rand()
  y1   = gen_rand()
  azi1 = gen_rand()
  nleg = size(p)
! To generate the intermediate PS point the emitter is chosen 
! to be the first massless parton and the position for the 
! radiated one is the last position:
  do emit1=3,nleg-2
    if (abs(flv(emit1)).gt.qcd_nf) cycle
    exit
  end do
  rad1 = nleg - 1
! Since the second emission only produces a soft parton its
! emitter can be freely chosen to be the same as for emitter 1:
  emit2 = emit1
  rad2  = radt
! The azimuth are picked randomly:
  azi2  = gen_rand()
! while the angle is chosen to be perpendicular to emitter:
  y2    = 0d0
!  print *,"New sequence is started: "
  write(iun,"(/a/)") "Singly unresolved region: "
  do iexp=1,nexp
    xi2 = 10d0**(-iexp)
    call MapFSR(pborn,xi1,y1,azi1,emit1,rad1,flv,phat)
! Using the intermediate PS point to build the RR one:
    call MapFSR(phat(:)%p,xi2,y2,azi2,emit2,rad2,flv,p)
! We have to calculate the common invariants for the real
! kinematics:
!    call PrintParts(p)
    call CalcSubInvariants(p)
    E_t = p(radt)%p%E / sqrt(Q2)
!    print *,"E_t: ",E_t
! We have to calculate the scales:
    call calcmyscales(p)
    call PickCirsCSirs(p,phat,ptilde,emitirt,radi,radr,radt, &
                       CSirt,CirtCSirtterm)
    call PickStCirtCSirt(p,phat,ptilde,          &
                         emitirt,radi,radr,radt, &
                         CSirt,StCirtCSirtterm)
    write(iun,fmt='(a,I3,1x)',advance='no') "iexp= ",int(log10(E_t))
    write(iun,*) ", StCirtCSirt/CirtCSirt= ",StCirtCSirtterm/CirtCSirtterm
!
  end do
!
end subroutine StCirtCSirtLimit
!
subroutine StCirtSrtLimit(iun,pborn,p,phat,ptilde,flv, &
                          emitirt,radi,radr,radt,St,Srt)
use particles
use random
use regions
use math
use FKSutils
use scales
use my_scale
use coupling
use misc
use QCDparams
implicit none
!
  integer , intent(in) :: iun
  type(mom) , dimension(:) , intent(in) :: pborn
  type(particle) , dimension(:) , intent(inout) :: p,ptilde,phat
  integer , dimension(:) , intent(in) :: flv
  integer , intent(in) :: emitirt,radi,radr,radt
  type(subterms) , intent(in) :: St,Srt
!
  integer :: iexp
  integer :: nleg
  integer :: emit1,emit2,rad1,rad2
  integer :: rad_i,rad_r,rad_s
  real(kind(1d0)) :: xi1,xi2,azi1,azi2,y1,y2
  real(kind(1d0)) :: E_t
  real(kind(1d0)) :: CirtSrtterm,StCirtSrtterm
!
!
  interface 
    subroutine CalcSubInvariants(p)
    use regions
    use particles
!
      type(particle) , dimension(:) , intent(in) :: p
!
    end subroutine CalcSubInvariants
!
    subroutine PickCirsSrs(p,ptilde,emit,radi,radr,rads, &
                           Srs,CirsSrsterm)
    use regions
    use particles
    implicit none
!
      type(particle), dimension(:), intent(in) :: p
      type(particle), dimension(:), intent(out) :: ptilde
      integer , intent(in) :: emit,radi,radr,rads
      type(subterms) , intent(in) :: Srs
      real(kind(1d0)) , intent(out) :: CirsSrsterm
!
    end subroutine PickCirsSrs
!
    subroutine PickStCirtSrt(p,phat,ptilde,          &
                             emitirt,radi,radr,radt, &
                             Srt,StCirtSrtterm)
    use regions
    use particles
    implicit none
!
      type(particle), dimension(:), intent(in) :: p
      type(particle), dimension(:), intent(out) :: phat,ptilde
      integer , intent(in) :: emitirt,radi,radr,radt
      type(subterms) , intent(in) :: Srt
      real(kind(1d0)) , intent(out) :: StCirtSrtterm
!
    end subroutine PickStCirtSrt
  end interface
!
!  print *,"************************************************************"
!  print *,"We will approach the limit of: "
!  print *,"emitktr,radk,radt,radr: ",emitktr,radk,radt,radr
!  call PrintSubProc(flv)
!  print *,"************************************************************"
!  read(*,*)
!
! First, building a PS point which has n-1 legs:
  xi1  = gen_rand()
  y1   = gen_rand()
  azi1 = gen_rand()
  nleg = size(p)
  rad_i = radi
  rad_r = min(radr,radt)
  rad_s = max(radr,radt)
! To generate the intermediate PS point the emitter is chosen 
! to be the first massless parton and the position for the 
! radiated one is the last position:
  do emit1=3,nleg-2
    if (abs(flv(emit1)).gt.qcd_nf) cycle
    exit
  end do
  rad1 = nleg - 1
! Since the second emission only produces a soft parton its
! emitter can be freely chosen to be the same as for emitter 1:
  emit2 = emit1
  rad2  = radt
! The azimuth are picked randomly:
  azi2  = gen_rand()
! while the angle is chosen to be perpendicular to emitter:
  y2    = 0d0
!  print *,"New sequence is started: "
  write(iun,"(/a/)") "Singly unresolved region: "
  do iexp=1,nexp
    xi2 = 10d0**(-iexp)
    call MapFSR(pborn,xi1,y1,azi1,emit1,rad1,flv,phat)
! Using the intermediate PS point to build the RR one:
    call MapFSR(phat(:)%p,xi2,y2,azi2,emit2,rad2,flv,p)
! We have to calculate the common invariants for the real
! kinematics:
!    call PrintParts(p)
    call CalcSubInvariants(p)
    E_t = p(radt)%p%E / sqrt(Q2)
!    print *,"E_t: ",E_t
! We have to calculate the scales:
    call calcmyscales(p)
    call PickCirsSrs(p,ptilde,emitirt,rad_i,rad_r,rad_s, &
                     Srt,CirtSrtterm)
    call PickStCirtSrt(p,phat,ptilde,          &
                       emitirt,radi,radr,radt, &
                       Srt,StCirtSrtterm)
    write(iun,fmt='(a,I3,1x)',advance='no') "iexp= ",int(log10(E_t))
    write(iun,*) ", StCirtSrt/CirtSrt= ",StCirtSrtterm/CirtSrtterm
!
  end do
!
end subroutine StCirtSrtLimit
!
subroutine StCSirtSrtLimit(iun,pborn,p,phat,ptilde,Bij,flv, &
                           emitir,radi,radr,radt,St,Srt)
use particles
use random
use regions
use math
use FKSutils
use scales
use my_scale
use coupling
use misc
use QCDparams
implicit none
!
  integer , intent(in) :: iun
  type(mom) , dimension(:) , intent(in) :: pborn
  type(particle) , dimension(:) , intent(inout) :: p,ptilde,phat
  real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
  integer , dimension(:) , intent(in) :: flv
  integer , intent(in) :: emitir,radi,radr,radt
  type(subterms) , intent(in) :: St,Srt
!
  integer :: iexp
  integer :: nleg
  integer :: emit1,emit2,rad1,rad2
  real(kind(1d0)) :: xi1,xi2,azi1,azi2,y1,y2
  real(kind(1d0)) :: E_t
  real(kind(1d0)) :: CSirtSrtterm,StCSirtSrtterm
!
!
  interface 
    subroutine CalcSubInvariants(p)
    use regions
    use particles
!
      type(particle) , dimension(:) , intent(in) :: p
!
    end subroutine CalcSubInvariants
!
    subroutine PickCSirsSrs(p,ptilde,Bij,emitir,radi,radr,rads, &
                            Srs,CSirsSrsterm)
    use regions
    use particles
    implicit none
!
      type(particle), dimension(:), intent(in) :: p
      type(particle), dimension(:), intent(out) :: ptilde
      real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
      integer , intent(in) :: emitir,radi,radr,rads
      type(subterms) , intent(in) :: Srs
      real(kind(1d0)) , intent(out) :: CSirsSrsterm
!
    end subroutine PickCSirsSrs
!
    subroutine PickStCSirtSrt(p,phat,ptilde,Bij,     &
                              emitir,radi,radr,radt, &
                              Srt,StCSirtSrtterm)
    use regions
    use particles
    implicit none
!
      type(particle), dimension(:), intent(in) :: p
      type(particle), dimension(:), intent(out) :: phat,ptilde
      real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
      integer , intent(in) :: emitir,radi,radr,radt
      type(subterms) , intent(in) :: Srt
      real(kind(1d0)) , intent(out) :: StCSirtSrtterm
!
    end subroutine PickStCSirtSrt
  end interface
!
!  print *,"************************************************************"
!  print *,"We will approach the limit of: "
!  print *,"emitirt,radi,radr,radt: ",emitirt,radi,radr,radt
!  call PrintSubProc(flv)
!  print *,"************************************************************"
!  read(*,*)
!
! First, building a PS point which has n-1 legs:
  xi1  = gen_rand()
  y1   = gen_rand()
  azi1 = gen_rand()
  nleg = size(p)
! To generate the intermediate PS point the emitter is chosen 
! to be the first massless parton and the position for the 
! radiated one is the last position:
  do emit1=3,nleg-2
    if (abs(flv(emit1)).gt.qcd_nf) cycle
    exit
  end do
  rad1 = nleg - 1
! Since the second emission only produces a soft parton its
! emitter can be freely chosen to be the same as for emitter 1:
  emit2 = emit1
  rad2  = radt
! The azimuth are picked randomly:
  azi2  = gen_rand()
! while the angle is chosen to be perpendicular to emitter:
  y2    = 0d0
!  print *,"New sequence is started: "
  write(iun,"(/a/)") "Singly unresolved region: "
  do iexp=1,nexp
    xi2 = 10d0**(-iexp)
    call MapFSR(pborn,xi1,y1,azi1,emit1,rad1,flv,phat)
! Using the intermediate PS point to build the RR one:
    call MapFSR(phat(:)%p,xi2,y2,azi2,emit2,rad2,flv,p)
! We have to calculate the common invariants for the real
! kinematics:
!    call PrintParts(p)
    call CalcSubInvariants(p)
    E_t = p(radt)%p%E / sqrt(Q2)
!    print *,"E_t: ",E_t
! We have to calculate the scales:
    call calcmyscales(p)
    call PickCSirsSrs(p,ptilde,Bij,emitir,radi,radr,radt, &
                      Srt,CSirtSrtterm)
    call PickStCSirtSrt(p,phat,ptilde,Bij,     &
                        emitir,radi,radr,radt, &
                        Srt,StCSirtSrtterm)
    write(iun,fmt='(a,I3,1x)',advance='no') "iexp= ",int(log10(E_t))
    write(iun,*) ", StCSirtSrt/CSirtSrt= ",StCSirtSrtterm/CSirtSrtterm
!
  end do
!
end subroutine StCSirtSrtLimit
!
subroutine StCirtCSirtSrtLimit(iun,pborn,p,phat,ptilde,flv, &
                               emitirt,radi,radr,radt,St,Srt)
use particles
use random
use regions
use math
use FKSutils
use scales
use my_scale
use coupling
use misc
use QCDparams
implicit none
!
  integer , intent(in) :: iun
  type(mom) , dimension(:) , intent(in) :: pborn
  type(particle) , dimension(:) , intent(inout) :: p,ptilde,phat
  integer , dimension(:) , intent(in) :: flv
  integer , intent(in) :: emitirt,radi,radr,radt
  type(subterms) , intent(in) :: St,Srt
!
  integer :: iexp
  integer :: nleg
  integer :: emit1,emit2,rad1,rad2
  real(kind(1d0)) :: xi1,xi2,azi1,azi2,y1,y2
  real(kind(1d0)) :: E_t
  real(kind(1d0)) :: CirtCSirtSrtterm,StCirtCSirtSrtterm
!
!
  interface 
    subroutine CalcSubInvariants(p)
    use regions
    use particles
!
      type(particle) , dimension(:) , intent(in) :: p
!
    end subroutine CalcSubInvariants
!
    subroutine PickCirsCSirsSrs(p,ptilde,emitirs,radi,radr,rads, &
                                Srs,CirsCSirsSrsterm)
    use regions
    use particles
    implicit none
!
      type(particle), dimension(:), intent(in) :: p
      type(particle), dimension(:), intent(out) :: ptilde
      integer , intent(in) :: emitirs,radi,radr,rads
      type(subterms) , intent(in) :: Srs
      real(kind(1d0)) , intent(out) :: CirsCSirsSrsterm
!
    end subroutine PickCirsCSirsSrs
!
    subroutine PickStCirtCSirtSrt(p,phat,ptilde,          &
                                  emitirk,radi,radr,radt, &
                                  Srt,StCirtCSirtSrtterm)
    use regions
    use particles
    implicit none
!
      type(particle), dimension(:), intent(in) :: p
      type(particle), dimension(:), intent(out) :: phat,ptilde
      integer , intent(in) :: emitirk,radi,radr,radt
      type(subterms) , intent(in) :: Srt
      real(kind(1d0)) , intent(out) :: StCirtCSirtSrtterm
!
    end subroutine PickStCirtCSirtSrt
  end interface
!
!  print *,"************************************************************"
!  print *,"We will approach the limit of: "
!  print *,"emitirt,radi,radr,radt: ",emitirt,radi,radr,radt
!  call PrintSubProc(flv)
!  print *,"************************************************************"
!  read(*,*)
!
! First, building a PS point which has n-1 legs:
  xi1  = gen_rand()
  y1   = gen_rand()
  azi1 = gen_rand()
  nleg = size(p)
! To generate the intermediate PS point the emitter is chosen 
! to be the first massless parton and the position for the 
! radiated one is the last position:
  do emit1=3,nleg-2
    if (abs(flv(emit1)).gt.qcd_nf) cycle
    exit
  end do
  rad1 = nleg - 1
! Since the second emission only produces a soft parton its
! emitter can be freely chosen to be the same as for emitter 1:
  emit2 = emit1
  rad2  = radt
! The azimuth are picked randomly:
  azi2  = gen_rand()
! while the angle is chosen to be perpendicular to emitter:
  y2    = 0d0
!  print *,"New sequence is started: "
  write(iun,"(/a/)") "Singly unresolved region: "
  do iexp=1,nexp
    xi2 = 10d0**(-iexp)
    call MapFSR(pborn,xi1,y1,azi1,emit1,rad1,flv,phat)
! Using the intermediate PS point to build the RR one:
    call MapFSR(phat(:)%p,xi2,y2,azi2,emit2,rad2,flv,p)
! We have to calculate the common invariants for the real
! kinematics:
!    call PrintParts(p)
    call CalcSubInvariants(p)
    E_t = p(radt)%p%E / sqrt(Q2)
!    print *,"E_t: ",E_t
! We have to calculate the scales:
    call calcmyscales(p)
    call PickCirsCSirsSrs(p,ptilde,emitirt,radi,radr,radt, &
                          Srt,CirtCSirtSrtterm)
    call PickStCirtCSirtSrt(p,phat,ptilde,         &
                            emitirt,radi,radr,radt, &
                            Srt,StCirtCSirtSrtterm)
    write(iun,fmt='(a,I3,1x)',advance='no') "iexp= ",int(log10(E_t))
    write(iun,fmt='(a,G0)') ", StCirtCSirtSrt/CirtCSirtSrt= ", &
                            StCirtCSirtSrtterm/CirtCSirtSrtterm
!
  end do
!
end subroutine StCirtCSirtSrtLimit
!
subroutine StSrtLimit(iun,pborn,p,phat,ptilde,Bij,Bijkl,flv, &
                      radr,radt,St,Srt)
use particles
use random
use regions
use math
use FKSutils
use scales
use my_scale
use coupling
use misc
use QCDparams
implicit none
!
  integer , intent(in) :: iun
  type(mom) , dimension(:) , intent(in) :: pborn
  type(particle) , dimension(:) , intent(inout) :: p,ptilde,phat
  real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
  real(kind(1d0)) , dimension(:,:,:,:) , intent(inout) :: Bijkl
  integer , dimension(:) , intent(in) :: flv
  integer , intent(in) :: radr,radt
  type(subterms) , intent(in) :: St,Srt
!
  integer :: iexp
  integer :: nleg
  integer :: emit1,emit2,rad1,rad2
  integer :: rad_r,rad_s
  real(kind(1d0)) :: xi1,xi2,azi1,azi2,y1,y2
  real(kind(1d0)) :: E_t
  real(kind(1d0)) :: Srtterm,StSrtterm
!
!
  interface 
    subroutine CalcSubInvariants(p)
    use regions
    use particles
!
      type(particle) , dimension(:) , intent(in) :: p
!
    end subroutine CalcSubInvariants
!
    subroutine PickSrs(p,ptilde,Bij,Bijkl,radr,rads, &
                       Srs,Srsterm)
    use regions
    use particles
    implicit none
!
      type(particle), dimension(:), intent(in) :: p
      type(particle), dimension(:), intent(out) :: ptilde
      real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
      real(kind(1d0)) , dimension(:,:,:,:) , intent(inout) :: Bijkl
      integer , intent(in) :: radr,rads
      type(subterms) , intent(in) :: Srs
      real(kind(1d0)) , intent(out) :: Srsterm
!
    end subroutine PickSrs
!
    subroutine PickStSrt(p,phat,ptilde,Bij,Bijkl, &
                         radr,radt,Srt,StSrtterm)
    use regions
    use particles
    implicit none
!
      type(particle), dimension(:), intent(in) :: p
      type(particle), dimension(:), intent(out) :: phat,ptilde
      real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
      real(kind(1d0)) , dimension(:,:,:,:) , intent(inout) :: Bijkl
      integer , intent(in) :: radr,radt
      type(subterms) , intent(in) :: Srt
      real(kind(1d0)) , intent(out) :: StSrtterm
!
    end subroutine PickStSrt
  end interface
!
!  print *,"************************************************************"
!  print *,"We will approach the limit of: "
!  print *,"radr,radt: ",radr,radt
!  call PrintSubProc(flv)
!  print *,"************************************************************"
!  read(*,*)
!
! First, building a PS point which has n-1 legs:
  xi1  = gen_rand()
  y1   = gen_rand()
  azi1 = gen_rand()
  nleg = size(p)
  rad_r = min(radr,radt)
  rad_s = max(radr,radt)
! To generate the intermediate PS point the emitter is chosen 
! to be the first massless parton and the position for the 
! radiated one is the last position:
  do emit1=3,nleg-2
    if (abs(flv(emit1)).gt.qcd_nf) cycle
    exit
  end do
  rad1 = nleg - 1
! Since the second emission only produces a soft parton its
! emitter can be freely chosen to be the same as for emitter 1:
  emit2 = emit1
  rad2  = radt
! The azimuth are picked randomly:
  azi2  = gen_rand()
! while the angle is chosen to be perpendicular to emitter:
  y2    = 0d0
!  print *,"New sequence is started: "
  write(iun,"(/a/)") "Singly unresolved region: "
  do iexp=1,nexp
    xi2 = 10d0**(-iexp)
    call MapFSR(pborn,xi1,y1,azi1,emit1,rad1,flv,phat)
! Using the intermediate PS point to build the RR one:
    call MapFSR(phat(:)%p,xi2,y2,azi2,emit2,rad2,flv,p)
! We have to calculate the common invariants for the real
! kinematics:
!    call PrintParts(p)
    call CalcSubInvariants(p)
    E_t = p(radt)%p%E / sqrt(Q2)
!    print *,"E_t: ",E_t
! We have to calculate the scales:
    call calcmyscales(p)
    call PickSrs(p,ptilde,Bij,Bijkl,rad_r,rad_s, &
                 Srt,Srtterm)
    call PickStSrt(p,phat,ptilde,Bij,Bijkl, &
                   radr,radt,Srt,StSrtterm)
    write(iun,fmt='(a,I3,1x)',advance='no') "iexp= ",int(log10(E_t))
    write(iun,fmt='(a,G0)') ", StSrt/Srt= ", &
                            StSrtterm/Srtterm
!
  end do
!
end subroutine StSrtLimit
!
subroutine CitStCirtLimit(iun,pborn,p,phat,ptilde,flv, &
                          radi,radr,radt,CSit,Cirt)
use particles
use random
use regions
use math
use FKSutils
use scales
use my_scale
use coupling
use misc
use QCDparams
implicit none
!
  integer , intent(in) :: iun
  type(mom) , dimension(:) , intent(in) :: pborn
  type(particle) , dimension(:) , intent(inout) :: p,ptilde,phat
  integer , dimension(:) , intent(in) :: flv
  integer , intent(in) :: radi,radr,radt
  type(subterms) , intent(in) :: CSit,Cirt
!
  logical :: swap
  integer :: iexp
  integer :: nleg
  integer :: emit1,emit2,rad1,rad2
  integer :: emitirs,rad_i,rad_r,rad_s
  real(kind(1d0)) :: xi1,xi2,azi1,azi2,y1,y2
  real(kind(1d0)) :: y_it,E_t
  real(kind(1d0)) :: Cirtterm,CitStCirtterm
  type(mom) :: p_tmp
!
!
  interface 
    subroutine CalcSubInvariants(p)
    use regions
    use particles
!
      type(particle) , dimension(:) , intent(in) :: p
!
    end subroutine CalcSubInvariants
!
    subroutine PickCirs(p,ptilde,emit,radi,radr,rads,Cirs,Cirsterm)
    use regions
    use particles
    implicit none
!
      type(particle), dimension(:), intent(in) :: p
      type(particle), dimension(:), intent(out) :: ptilde
      integer , intent(in) :: emit,radi,radr,rads
      type(subterms) , intent(in) :: Cirs
      real(kind(1d0)) , intent(out) :: Cirsterm
!
    end subroutine PickCirs
!
    subroutine PickCitStCirt(p,phat,ptilde,      &
                             radi,radr,radt,     &
                             rad_i,rad_r,rad_s,  &
                             Cirt,CitStCirtterm)
    use regions
    use particles
    implicit none
!
      type(particle), dimension(:), intent(in) :: p
      type(particle), dimension(:), intent(out) :: phat,ptilde
      integer , intent(in) :: radi,radr,radt,rad_i,rad_r,rad_s
      type(subterms) , intent(in) :: Cirt
      real(kind(1d0)) , intent(out) :: CitStCirtterm
!
    end subroutine PickCitStCirt
  end interface
!
!  print *,"************************************************************"
!  print *,"We will approach the limit of: "
!  print *,"radi,radr,radt: ",radi,radr,radt
!  call PrintSubProc(flv)
!  print *,"************************************************************"
!  read(*,*)
!
! First, building a PS point which has n-1 legs:
  xi1  = gen_rand()
  y1   = gen_rand()
  azi1 = gen_rand()
  nleg = size(p)
! Obtaining the parton order for the Cirs counterterm: 
  call ReorderRadirs(radi,flv(radi), &
                     radr,flv(radr), &
                     radt,flv(radt), &
                     rad_i,rad_r,rad_s)
  emitirs = min(radi,radr,radt)
! To generate the intermediate PS point the emitter is chosen 
! to be the first massless parton and the position for the 
! radiated one is the last position:
  do emit1=3,nleg-2
    if (abs(flv(emit1)).gt.qcd_nf) cycle
    exit
  end do
  rad1 = nleg - 1
! legs i and t are not ordered hence it is possible that
! i > t which can produce a problem when phase space is
! constructed:
  swap = .false.
  if (radi.lt.radt) then
    emit2 = radi
    rad2  = radt
  else
    swap = .true.
    emit2 = radt
    rad2  = radi
  end if
! The azimuth are picked randomly:
  azi2  = gen_rand()
!  print *,"New sequence is started: "
  write(iun,"(/a/)") "Singly unresolved region: "
  do iexp=1,nexp
    xi2 = 10d0**(-iexp)
    xi2 = sqrt(xi2)
    xi2 = (xi2*(1 - 2*sqrt(tiny_y)) + sqrt(tiny_y) - tiny_xi) &
        / (1 - 2*tiny_xi)
    y2  = 1d0 - 10d0**(-iexp)
    call MapFSR(pborn,xi1,y1,azi1,emit1,rad1,flv,phat)
! Using the intermediate PS point to build the RR one:
    call MapFSR(phat(:)%p,xi2,y2,azi2,emit2,rad2,flv,p)
    if (swap) then
      p_tmp     = p(radi)%p
      p(radi)%p = p(radt)%p
      p(radt)%p = p_tmp
    end if
! We have to calculate the common invariants for the real
! kinematics:
!    call PrintParts(p)
    call CalcSubInvariants(p)
    y_it = (p(radi)%p + p(radt)%p)**2 / Q2
    E_t = p(radt)%p%E / sqrt(Q2)
!    print *,"E_t: ",E_t
! We have to calculate the scales:
    call calcmyscales(p)
    call PickCirs(p,ptilde,emitirs,rad_i,rad_r,rad_s, &
                  Cirt,Cirtterm)
    call PickCitStCirt(p,phat,ptilde,      &
                       radi,radr,radt,     &
                       rad_i,rad_r,rad_s,  &
                       Cirt,CitStCirtterm)
    write(iun,fmt='(a,2(I3,1x))',advance='no') "iexp= ",int(log10(y_it)),int(log10(E_t))
    write(iun,fmt='(a,G0)') ", CitStCirt/Cirt= ", &
                            CitStCirtterm/Cirtterm
!
  end do
!
end subroutine CitStCirtLimit
!
subroutine CktStCSirtLimit(iun,pborn,p,phat,ptilde,Bij,Bmunuij, &
                           flv,radi,radr,radk,radt,CSit,CSirt)
use particles
use random
use regions
use math
use FKSutils
use scales
use my_scale
use coupling
use misc
use QCDparams
implicit none
!
  integer , intent(in) :: iun
  type(mom) , dimension(:) , intent(in) :: pborn
  type(particle) , dimension(:) , intent(inout) :: p,ptilde,phat
  real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
  real(kind(1d0)) , dimension(0:,0:,:,:) , intent(inout) :: Bmunuij
  integer , dimension(:) , intent(in) :: flv
  integer , intent(in) :: radi,radr,radk,radt
  type(subterms) , intent(in) :: CSit,CSirt
!
  logical :: swap
  integer :: iexp
  integer :: nleg
  integer :: emit1,emit2,rad1,rad2
  integer :: emit_ir
  real(kind(1d0)) :: xi1,xi2,azi1,azi2,y1,y2
  real(kind(1d0)) :: y_kt,E_t
  real(kind(1d0)) :: CSirtterm,CktStCSirtterm
  type(mom) :: p_tmp
!
!
  interface 
    subroutine CalcSubInvariants(p)
    use regions
    use particles
!
      type(particle) , dimension(:) , intent(in) :: p
!
    end subroutine CalcSubInvariants
!
    subroutine PickCSirs(p,phat,ptilde,Bij,Bmunuij,emit,radi,radr,rads, &
                         CSirs,CSirsterm)
    use regions
    use particles
    implicit none
!
      type(particle), dimension(:), intent(in) :: p
      type(particle), dimension(:), intent(out) :: phat,ptilde
      real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
      real(kind(1d0)) , dimension(0:,0:,:,:) , intent(inout) :: Bmunuij
      integer , intent(in) :: emit,radi,radr,rads
      type(subterms) , intent(in) :: CSirs
      real(kind(1d0)) , intent(out) :: CSirsterm
!
    end subroutine PickCSirs
!
    subroutine PickCktStCSirt(p,phat,ptilde,       &
                              radi,radr,radk,radt, &
                              CSirt,CktStCSirtterm)
    use regions
    use particles
    implicit none
!
      type(particle), dimension(:), intent(in) :: p
      type(particle), dimension(:), intent(out) :: phat,ptilde
      integer , intent(in) :: radi,radr,radk,radt
      type(subterms) , intent(in) :: CSirt
      real(kind(1d0)) , intent(out) :: CktStCSirtterm
!
    end subroutine PickCktStCSirt
  end interface
!
!  print *,"************************************************************"
!  print *,"We will approach the limit of: "
!  print *,"radi,radr,radk,radt: ",radi,radr,radk,radt
!  call PrintSubProc(flv)
!  print *,"************************************************************"
!  read(*,*)
!
! First, building a PS point which has n-1 legs:
  xi1  = gen_rand()
  y1   = gen_rand()
  azi1 = gen_rand()
  nleg = size(p)
  emit_ir = min(radi,radr)
  if (emit_ir.gt.radt) emit_ir = emit_ir - 1
! To generate the intermediate PS point the emitter is chosen 
! to be the first massless parton and the position for the 
! radiated one is the last position:
  do emit1=3,nleg-2
    if (abs(flv(emit1)).gt.qcd_nf) cycle
    exit
  end do
  rad1 = nleg - 1
! legs k and t are not ordered hence it is possible that
! k > t which can produce a problem when phase space is
! constructed:
  swap = .false.
  if (radk.lt.radt) then
    emit2 = radk
    rad2  = radt
  else
    swap = .true.
    emit2 = radt
    rad2  = radk
  end if
! The azimuth are picked randomly:
  azi2  = gen_rand()
!  print *,"New sequence is started: "
  write(iun,"(/a/)") "Singly unresolved region: "
  do iexp=1,nexp
    xi2 = 10d0**(-iexp)
    xi2 = sqrt(xi2)
    xi2 = (xi2*(1 - 2*sqrt(tiny_y)) + sqrt(tiny_y) - tiny_xi) &
        / (1 - 2*tiny_xi)
    y2  = 1d0 - 10d0**(-iexp)
    call MapFSR(pborn,xi1,y1,azi1,emit1,rad1,flv,phat)
! Using the intermediate PS point to build the RR one:
    call MapFSR(phat(:)%p,xi2,y2,azi2,emit2,rad2,flv,p)
    if (swap) then
      p_tmp     = p(radk)%p
      p(radk)%p = p(radt)%p
      p(radt)%p = p_tmp
    end if
! We have to calculate the common invariants for the real
! kinematics:
!    call PrintParts(p)
    call CalcSubInvariants(p)
    y_kt = (p(radk)%p + p(radt)%p)**2 / Q2
    E_t = p(radt)%p%E / sqrt(Q2)
!    print *,"E_t: ",E_t
! We have to calculate the scales:
    call calcmyscales(p)
    call PickCSirs(p,phat,ptilde,Bij,Bmunuij, &
                   emit_ir,radi,radr,radt,    &
                   CSirt,CSirtterm)
    call PickCktStCSirt(p,phat,ptilde,       &
                        radi,radr,radk,radt, &
                        CSirt,CktStCSirtterm)
    write(iun,fmt='(a,2(I3,1x))',advance='no') "iexp= ",int(log10(y_kt)),int(log10(E_t))
    write(iun,fmt='(a,G0)') ", CktStCSirt/CSirt= ", &
                            CktStCSirtterm/CSirtterm
!
  end do
!
end subroutine CktStCSirtLimit
!
subroutine CktStCSirtSrtLimit(iun,pborn,p,phat,ptilde,Bij, &
                              flv,emitir,radi,radr,radk,radt,CSit,Srt)
use particles
use random
use regions
use math
use FKSutils
use scales
use my_scale
use coupling
use misc
use QCDparams
implicit none
!
  integer , intent(in) :: iun
  type(mom) , dimension(:) , intent(in) :: pborn
  type(particle) , dimension(:) , intent(inout) :: p,ptilde,phat
  real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
  integer , dimension(:) , intent(in) :: flv
  integer , intent(in) :: emitir,radi,radr,radk,radt
  type(subterms) , intent(in) :: CSit,Srt
!
  logical :: swap
  integer :: iexp
  integer :: nleg
  integer :: emit1,emit2,rad1,rad2
  real(kind(1d0)) :: xi1,xi2,azi1,azi2,y1,y2
  real(kind(1d0)) :: y_kt,E_t
  real(kind(1d0)) :: CSirtSrtterm,CktStCSirtSrtterm
  type(mom) :: p_tmp
!
!
  interface 
    subroutine CalcSubInvariants(p)
    use regions
    use particles
!
      type(particle) , dimension(:) , intent(in) :: p
!
    end subroutine CalcSubInvariants
!
    subroutine PickCSirsSrs(p,ptilde,Bij,emitir,radi,radr,rads, &
                            Srs,CSirsSrsterm)
    use regions
    use particles
    implicit none
!
      type(particle), dimension(:), intent(in) :: p
      type(particle), dimension(:), intent(out) :: ptilde
      real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
      integer , intent(in) :: emitir,radi,radr,rads
      type(subterms) , intent(in) :: Srs
      real(kind(1d0)) , intent(out) :: CSirsSrsterm
!
    end subroutine PickCSirsSrs
!
    subroutine PickCktStCSirtSrt(p,phat,ptilde,       &
                                 radi,radr,radk,radt, &
                                 Srt,CktStCSirtSrtterm)
    use regions
    use particles
    implicit none
!
      type(particle), dimension(:), intent(in) :: p
      type(particle), dimension(:), intent(out) :: phat,ptilde
      integer , intent(in) :: radi,radr,radk,radt
      type(subterms) , intent(in) :: Srt
      real(kind(1d0)) , intent(out) :: CktStCSirtSrtterm
!
    end subroutine PickCktStCSirtSrt
  end interface
!
!  print *,"************************************************************"
!  print *,"We will approach the limit of: "
!  print *,"emitir,radi,radr,radk,radt: ",emitir,radi,radr,radk,radt
!  call PrintSubProc(flv)
!  print *,"************************************************************"
!  read(*,*)
!
! First, building a PS point which has n-1 legs:
  xi1  = gen_rand()
  y1   = gen_rand()
  azi1 = gen_rand()
  nleg = size(p)
! To generate the intermediate PS point the emitter is chosen 
! to be the first massless parton and the position for the 
! radiated one is the last position:
  do emit1=3,nleg-2
    if (abs(flv(emit1)).gt.qcd_nf) cycle
    exit
  end do
  rad1 = nleg - 1
! legs k and t are not ordered hence it is possible that
! k > t which can produce a problem when phase space is
! constructed:
  swap = .false.
  if (radk.lt.radt) then
    emit2 = radk
    rad2  = radt
  else
    swap = .true.
    emit2 = radt
    rad2  = radk
  end if
! The azimuth are picked randomly:
  azi2  = gen_rand()
!  print *,"New sequence is started: "
  write(iun,"(/a/)") "Singly unresolved region: "
  do iexp=1,nexp
    xi2 = 10d0**(-iexp)
    xi2 = sqrt(xi2)
    xi2 = (xi2*(1 - 2*sqrt(tiny_y)) + sqrt(tiny_y) - tiny_xi) &
        / (1 - 2*tiny_xi)
    y2  = 1d0 - 10d0**(-iexp)
    call MapFSR(pborn,xi1,y1,azi1,emit1,rad1,flv,phat)
! Using the intermediate PS point to build the RR one:
    call MapFSR(phat(:)%p,xi2,y2,azi2,emit2,rad2,flv,p)
    if (swap) then
      p_tmp     = p(radk)%p
      p(radk)%p = p(radt)%p
      p(radt)%p = p_tmp
    end if
! We have to calculate the common invariants for the real
! kinematics:
!    call PrintParts(p)
    call CalcSubInvariants(p)
    y_kt = (p(radk)%p + p(radt)%p)**2 / Q2
    E_t = p(radt)%p%E / sqrt(Q2)
!    print *,"E_t: ",E_t
! We have to calculate the scales:
    call calcmyscales(p)
    call PickCSirsSrs(p,ptilde,Bij,          &
                      emitir,radi,radr,radt, &
                      Srt,CSirtSrtterm)
    call PickCktStCSirtSrt(p,phat,ptilde,       &
                           radi,radr,radk,radt, &
                           Srt,CktStCSirtSrtterm)
    write(iun,fmt='(a,2(I3,1x))',advance='no') "iexp= ",int(log10(y_kt)),int(log10(E_t))
    write(iun,fmt='(a,G0)') ", CktStCSirtSrt/CSirtSrt= ", &
                            CktStCSirtSrtterm/CSirtSrtterm
!
  end do
!
end subroutine CktStCSirtSrtLimit
!
subroutine CktStCkrtSrtLimit(iun,pborn,p,phat,ptilde, &
                             flv,emitkrt,radk,radr,radt,CSit,Srt)
use particles
use random
use regions
use math
use FKSutils
use scales
use my_scale
use coupling
use misc
use QCDparams
implicit none
!
  integer , intent(in) :: iun
  type(mom) , dimension(:) , intent(in) :: pborn
  type(particle) , dimension(:) , intent(inout) :: p,ptilde,phat
  integer , dimension(:) , intent(in) :: flv
  integer , intent(in) :: emitkrt,radk,radr,radt
  type(subterms) , intent(in) :: CSit,Srt
!
  logical :: swap
  integer :: iexp
  integer :: nleg
  integer :: emit1,emit2,rad1,rad2
  real(kind(1d0)) :: xi1,xi2,azi1,azi2,y1,y2
  real(kind(1d0)) :: y_kt,E_t
  real(kind(1d0)) :: CkrtSrtterm,CktStCkrtSrtterm
  type(mom) :: p_tmp
!
!
  interface 
    subroutine CalcSubInvariants(p)
    use regions
    use particles
!
      type(particle) , dimension(:) , intent(in) :: p
!
    end subroutine CalcSubInvariants
!
    subroutine PickCirsSrs(p,ptilde,emitirs,radi,radr,rads, &
                           Srs,CirsSrsterm)
    use regions
    use particles
    implicit none
!
      type(particle), dimension(:), intent(in) :: p
      type(particle), dimension(:), intent(out) :: ptilde
      integer , intent(in) :: emitirs,radi,radr,rads
      type(subterms) , intent(in) :: Srs
      real(kind(1d0)) , intent(out) :: CirsSrsterm
!
    end subroutine PickCirsSrs
!
    subroutine PickCktStCkrtSrt(p,phat,ptilde,  &
                                radk,radr,radt, &
                                Srt,CktStCkrtSrtterm)
    use regions
    use particles
    implicit none
!
      type(particle), dimension(:), intent(in) :: p
      type(particle), dimension(:), intent(out) :: phat,ptilde
      integer , intent(in) :: radk,radr,radt
      type(subterms) , intent(in) :: Srt
      real(kind(1d0)) , intent(out) :: CktStCkrtSrtterm
!
    end subroutine PickCktStCkrtSrt
  end interface
!
!  print *,"************************************************************"
!  print *,"We will approach the limit of: "
!  print *,"emitkrt,radk,radr,radt: ",emitkrt,radk,radr,radt
!  call PrintSubProc(flv)
!  print *,"************************************************************"
!  read(*,*)
!
! First, building a PS point which has n-1 legs:
  xi1  = gen_rand()
  y1   = gen_rand()
  azi1 = gen_rand()
  nleg = size(p)
! To generate the intermediate PS point the emitter is chosen 
! to be the first massless parton and the position for the 
! radiated one is the last position:
  do emit1=3,nleg-2
    if (abs(flv(emit1)).gt.qcd_nf) cycle
    exit
  end do
  rad1 = nleg - 1
! legs k and t are not ordered hence it is possible that
! k > t which can produce a problem when phase space is
! constructed:
  swap = .false.
  if (radk.lt.radt) then
    emit2 = radk
    rad2  = radt
  else
    swap = .true.
    emit2 = radt
    rad2  = radk
  end if
! The azimuth are picked randomly:
  azi2  = gen_rand()
!  print *,"New sequence is started: "
  write(iun,"(/a/)") "Singly unresolved region: "
  do iexp=1,nexp
    xi2 = 10d0**(-iexp)
    xi2 = sqrt(xi2)
    xi2 = (xi2*(1 - 2*sqrt(tiny_y)) + sqrt(tiny_y) - tiny_xi) &
        / (1 - 2*tiny_xi)
    y2  = 1d0 - 10d0**(-iexp)
    call MapFSR(pborn,xi1,y1,azi1,emit1,rad1,flv,phat)
! Using the intermediate PS point to build the RR one:
    call MapFSR(phat(:)%p,xi2,y2,azi2,emit2,rad2,flv,p)
    if (swap) then
      p_tmp     = p(radk)%p
      p(radk)%p = p(radt)%p
      p(radt)%p = p_tmp
    end if
! We have to calculate the common invariants for the real
! kinematics:
!    call PrintParts(p)
    call CalcSubInvariants(p)
    y_kt = (p(radk)%p + p(radt)%p)**2 / Q2
    E_t = p(radt)%p%E / sqrt(Q2)
!    print *,"E_t: ",E_t
! We have to calculate the scales:
    call calcmyscales(p)
    call PickCirsSrs(p,ptilde,               &
                     emitkrt,radk,radr,radt, &
                     Srt,CkrtSrtterm)
    call PickCktStCkrtSrt(p,phat,ptilde,  &
                          radk,radr,radt, &
                          Srt,CktStCkrtSrtterm)
    write(iun,fmt='(a,2(I3,1x))',advance='no') "iexp= ",int(log10(y_kt)),int(log10(E_t))
    write(iun,fmt='(a,G0)') ", CktStCkrtSrt/CkrtSrt= ", &
                            CktStCkrtSrtterm/CkrtSrtterm
!
  end do
!
end subroutine CktStCkrtSrtLimit
!
subroutine CrtStCkrtSrtLimit(iun,pborn,p,phat,ptilde, &
                             flv,emitkrt,radk,radr,radt,CSit,Srt)
use particles
use random
use regions
use math
use FKSutils
use scales
use my_scale
use coupling
use misc
use QCDparams
implicit none
!
  integer , intent(in) :: iun
  type(mom) , dimension(:) , intent(in) :: pborn
  type(particle) , dimension(:) , intent(inout) :: p,ptilde,phat
  integer , dimension(:) , intent(in) :: flv
  integer , intent(in) :: emitkrt,radk,radr,radt
  type(subterms) , intent(in) :: CSit,Srt
!
  logical :: swap
  integer :: iexp
  integer :: nleg
  integer :: emit1,emit2,rad1,rad2
  real(kind(1d0)) :: xi1,xi2,azi1,azi2,y1,y2
  real(kind(1d0)) :: y_rt,E_t
  real(kind(1d0)) :: CkrtSrtterm,CrtStCkrtSrtterm
  type(mom) :: p_tmp
!
!
  interface 
    subroutine CalcSubInvariants(p)
    use regions
    use particles
!
      type(particle) , dimension(:) , intent(in) :: p
!
    end subroutine CalcSubInvariants
!
    subroutine PickCirsSrs(p,ptilde,emitirs,radi,radr,rads, &
                           Srs,CirsSrsterm)
    use regions
    use particles
    implicit none
!
      type(particle), dimension(:), intent(in) :: p
      type(particle), dimension(:), intent(out) :: ptilde
      integer , intent(in) :: emitirs,radi,radr,rads
      type(subterms) , intent(in) :: Srs
      real(kind(1d0)) , intent(out) :: CirsSrsterm
!
    end subroutine PickCirsSrs
!
    subroutine PickCrtStCkrtSrt(p,phat,ptilde,  &
                                radk,radr,radt, &
                                Srt,CrtStCkrtSrtterm)
    use regions
    use particles
    implicit none
!
      type(particle), dimension(:), intent(in) :: p
      type(particle), dimension(:), intent(out) :: phat,ptilde
      integer , intent(in) :: radk,radr,radt
      type(subterms) , intent(in) :: Srt
      real(kind(1d0)) , intent(out) :: CrtStCkrtSrtterm
!
    end subroutine PickCrtStCkrtSrt
  end interface
!
!  print *,"************************************************************"
!  print *,"We will approach the limit of: "
!  print *,"emitkrt,radk,radr,radt: ",emitkrt,radk,radr,radt
!  call PrintSubProc(flv)
!  print *,"************************************************************"
!  read(*,*)
!
! First, building a PS point which has n-1 legs:
  xi1  = gen_rand()
  y1   = gen_rand()
  azi1 = gen_rand()
  nleg = size(p)
! To generate the intermediate PS point the emitter is chosen 
! to be the first massless parton and the position for the 
! radiated one is the last position:
  do emit1=3,nleg-2
    if (abs(flv(emit1)).gt.qcd_nf) cycle
    exit
  end do
  rad1 = nleg - 1
! legs k and t are not ordered hence it is possible that
! k > t which can produce a problem when phase space is
! constructed:
  swap = .false.
  if (radr.lt.radt) then
    emit2 = radr
    rad2  = radt
  else
    swap = .true.
    emit2 = radt
    rad2  = radr
  end if
! The azimuth are picked randomly:
  azi2  = gen_rand()
!  print *,"New sequence is started: "
  write(iun,"(/a/)") "Singly unresolved region: "
  do iexp=1,nexp
    xi2 = 10d0**(-iexp)
    xi2 = sqrt(xi2)
    xi2 = (xi2*(1 - 2*sqrt(tiny_y)) + sqrt(tiny_y) - tiny_xi) &
        / (1 - 2*tiny_xi)
    y2  = 1d0 - 10d0**(-iexp)
    call MapFSR(pborn,xi1,y1,azi1,emit1,rad1,flv,phat)
! Using the intermediate PS point to build the RR one:
    call MapFSR(phat(:)%p,xi2,y2,azi2,emit2,rad2,flv,p)
    if (swap) then
      p_tmp     = p(radr)%p
      p(radr)%p = p(radt)%p
      p(radt)%p = p_tmp
    end if
! We have to calculate the common invariants for the real
! kinematics:
!    call PrintParts(p)
    call CalcSubInvariants(p)
    y_rt = (p(radr)%p + p(radt)%p)**2 / Q2
    E_t = p(radt)%p%E / sqrt(Q2)
!    print *,"E_t: ",E_t
! We have to calculate the scales:
    call calcmyscales(p)
    call PickCirsSrs(p,ptilde,               &
                     emitkrt,radk,radr,radt, &
                     Srt,CkrtSrtterm)
    call PickCrtStCkrtSrt(p,phat,ptilde,  &
                          radk,radr,radt, &
                          Srt,CrtStCkrtSrtterm)
    write(iun,fmt='(a,2(I3,1x))',advance='no') "iexp= ",int(log10(y_rt)),int(log10(E_t))
    write(iun,fmt='(a,G0)') ", CrtStCkrtSrt/CkrtSrt= ", &
                            CrtStCkrtSrtterm/CkrtSrtterm
!
  end do
!
end subroutine CrtStCkrtSrtLimit
!
subroutine CktStSrtLimit(iun,pborn,p,phat,ptilde,Bij,Bijkl, &
                         flv,radk,radr,radt,CSit,Srt)
use particles
use random
use regions
use math
use FKSutils
use scales
use my_scale
use coupling
use misc
use QCDparams
implicit none
!
  integer , intent(in) :: iun
  type(mom) , dimension(:) , intent(in) :: pborn
  type(particle) , dimension(:) , intent(inout) :: p,ptilde,phat
  real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
  real(kind(1d0)) , dimension(:,:,:,:) , intent(inout) :: Bijkl
  integer , dimension(:) , intent(in) :: flv
  integer , intent(in) :: radk,radr,radt
  type(subterms) , intent(in) :: CSit,Srt
!
  logical :: swap
  integer :: iexp
  integer :: nleg
  integer :: emit1,emit2,rad1,rad2
  real(kind(1d0)) :: xi1,xi2,azi1,azi2,y1,y2
  real(kind(1d0)) :: y_kt,E_t
  real(kind(1d0)) :: Srtterm,CktStSrtterm
  type(mom) :: p_tmp
!
!
  interface 
    subroutine CalcSubInvariants(p)
    use regions
    use particles
!
      type(particle) , dimension(:) , intent(in) :: p
!
    end subroutine CalcSubInvariants
!
    subroutine PickSrs(p,ptilde,Bij,Bijkl,radr,rads, &
                       Srs,Srsterm)
    use regions
    use particles
    implicit none
!
      type(particle), dimension(:), intent(in) :: p
      type(particle), dimension(:), intent(out) :: ptilde
      real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
      real(kind(1d0)) , dimension(:,:,:,:) , intent(inout) :: Bijkl
      integer , intent(in) :: radr,rads
      type(subterms) , intent(in) :: Srs
      real(kind(1d0)) , intent(out) :: Srsterm
!
    end subroutine PickSrs
!
    subroutine PickCktStSrt(p,phat,ptilde,Bij, &
                            radk,radr,radt,    &
                            Srt,CktStSrtterm)
    use regions
    use particles
    implicit none
!
      type(particle), dimension(:), intent(in) :: p
      type(particle), dimension(:), intent(out) :: phat,ptilde
      real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
      integer , intent(in) :: radk,radr,radt
      type(subterms) , intent(in) :: Srt
      real(kind(1d0)) , intent(out) :: CktStSrtterm
!
    end subroutine PickCktStSrt
  end interface
!
!  print *,"************************************************************"
!  print *,"We will approach the limit of: "
!  print *,"radk,radr,radt: ",radk,radr,radt
!  call PrintSubProc(flv)
!  print *,"************************************************************"
!  read(*,*)
!
! First, building a PS point which has n-1 legs:
  xi1  = gen_rand()
  y1   = gen_rand()
  azi1 = gen_rand()
  nleg = size(p)
! To generate the intermediate PS point the emitter is chosen 
! to be the first massless parton and the position for the 
! radiated one is the last position:
  do emit1=3,nleg-2
    if (abs(flv(emit1)).gt.qcd_nf) cycle
    exit
  end do
  rad1 = nleg - 1
! legs k and t are not ordered hence it is possible that
! k > t which can produce a problem when phase space is
! constructed:
  swap = .false.
  if (radk.lt.radt) then
    emit2 = radk
    rad2  = radt
  else
    swap = .true.
    emit2 = radt
    rad2  = radk
  end if
! The azimuth are picked randomly:
  azi2  = gen_rand()
!  print *,"New sequence is started: "
  write(iun,"(/a/)") "Singly unresolved region: "
  do iexp=1,nexp
    xi2 = 10d0**(-iexp)
    xi2 = sqrt(xi2)
    xi2 = (xi2*(1 - 2*sqrt(tiny_y)) + sqrt(tiny_y) - tiny_xi) &
        / (1 - 2*tiny_xi)
    y2  = 1d0 - 10d0**(-iexp)
    call MapFSR(pborn,xi1,y1,azi1,emit1,rad1,flv,phat)
! Using the intermediate PS point to build the RR one:
    call MapFSR(phat(:)%p,xi2,y2,azi2,emit2,rad2,flv,p)
    if (swap) then
      p_tmp     = p(radk)%p
      p(radk)%p = p(radt)%p
      p(radt)%p = p_tmp
    end if
! We have to calculate the common invariants for the real
! kinematics:
!    call PrintParts(p)
    call CalcSubInvariants(p)
    y_kt = (p(radk)%p + p(radt)%p)**2 / Q2
    E_t = p(radt)%p%E / sqrt(Q2)
!    print *,"E_t: ",E_t
! We have to calculate the scales:
    call calcmyscales(p)
    call PickSrs(p,ptilde,Bij,Bijkl,radr,radt,Srt,Srtterm)
    call PickCktStSrt(p,phat,ptilde,Bij, &
                      radk,radr,radt,    &
                      Srt,CktStSrtterm)
    write(iun,fmt='(a,2(I3,1x))',advance='no') "iexp= ",int(log10(y_kt)),int(log10(E_t))
    write(iun,fmt='(a,G0)') ", CktStSrt/Srt= ", &
                            CktStSrtterm/Srtterm
!
  end do
!
end subroutine CktStSrtLimit
!
subroutine CrtStSrtLimit(iun,pborn,p,phat,ptilde,Bij,Bijkl, &
                         flv,radr,radt,CSit,Srt)
use particles
use random
use regions
use math
use FKSutils
use scales
use my_scale
use coupling
use misc
use QCDparams
implicit none
!
  integer , intent(in) :: iun
  type(mom) , dimension(:) , intent(in) :: pborn
  type(particle) , dimension(:) , intent(inout) :: p,ptilde,phat
  real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
  real(kind(1d0)) , dimension(:,:,:,:) , intent(inout) :: Bijkl
  integer , dimension(:) , intent(in) :: flv
  integer , intent(in) :: radr,radt
  type(subterms) , intent(in) :: CSit,Srt
!
  logical :: swap
  integer :: iexp
  integer :: nleg
  integer :: emit1,emit2,rad1,rad2
  real(kind(1d0)) :: xi1,xi2,azi1,azi2,y1,y2
  real(kind(1d0)) :: y_rt,E_t
  real(kind(1d0)) :: Srtterm,CrtStSrtterm
  type(mom) :: p_tmp
!
!
  interface 
    subroutine CalcSubInvariants(p)
    use regions
    use particles
!
      type(particle) , dimension(:) , intent(in) :: p
!
    end subroutine CalcSubInvariants
!
    subroutine PickSrs(p,ptilde,Bij,Bijkl,radr,rads, &
                       Srs,Srsterm)
    use regions
    use particles
    implicit none
!
      type(particle), dimension(:), intent(in) :: p
      type(particle), dimension(:), intent(out) :: ptilde
      real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
      real(kind(1d0)) , dimension(:,:,:,:) , intent(inout) :: Bijkl
      integer , intent(in) :: radr,rads
      type(subterms) , intent(in) :: Srs
      real(kind(1d0)) , intent(out) :: Srsterm
!
    end subroutine PickSrs
!
    subroutine PickCrtStSrt(p,phat,ptilde,Bij, &
                            radr,radt,         &
                            Srt,CrtStSrtterm)
    use regions
    use particles
    implicit none
!
      type(particle), dimension(:), intent(in) :: p
      type(particle), dimension(:), intent(out) :: phat,ptilde
      real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
      integer , intent(in) :: radr,radt
      type(subterms) , intent(in) :: Srt
      real(kind(1d0)) , intent(out) :: CrtStSrtterm
!
    end subroutine PickCrtStSrt
  end interface
!
!  print *,"************************************************************"
!  print *,"We will approach the limit of: "
!  print *,"radr,radt: ",radr,radt
!  call PrintSubProc(flv)
!  print *,"************************************************************"
!  read(*,*)
!
! First, building a PS point which has n-1 legs:
  xi1  = gen_rand()
  y1   = gen_rand()
  azi1 = gen_rand()
  nleg = size(p)
! To generate the intermediate PS point the emitter is chosen 
! to be the first massless parton and the position for the 
! radiated one is the last position:
  do emit1=3,nleg-2
    if (abs(flv(emit1)).gt.qcd_nf) cycle
    exit
  end do
  rad1 = nleg - 1
! legs k and t are not ordered hence it is possible that
! k > t which can produce a problem when phase space is
! constructed:
  swap = .false.
  if (radr.lt.radt) then
    emit2 = radr
    rad2  = radt
  else
    swap = .true.
    emit2 = radt
    rad2  = radr
  end if
! The azimuth are picked randomly:
  azi2  = gen_rand()
!  print *,"New sequence is started: "
  write(iun,"(/a/)") "Singly unresolved region: "
  do iexp=1,nexp
    xi2 = 10d0**(-iexp)
    xi2 = sqrt(xi2)
    xi2 = (xi2*(1 - 2*sqrt(tiny_y)) + sqrt(tiny_y) - tiny_xi) &
        / (1 - 2*tiny_xi)
    y2  = 1d0 - 10d0**(-iexp)
    call MapFSR(pborn,xi1,y1,azi1,emit1,rad1,flv,phat)
! Using the intermediate PS point to build the RR one:
    call MapFSR(phat(:)%p,xi2,y2,azi2,emit2,rad2,flv,p)
    if (swap) then
      p_tmp     = p(radr)%p
      p(radr)%p = p(radt)%p
      p(radt)%p = p_tmp
    end if
! We have to calculate the common invariants for the real
! kinematics:
!    call PrintParts(p)
    call CalcSubInvariants(p)
    y_rt = (p(radr)%p + p(radt)%p)**2 / Q2
    E_t = p(radt)%p%E / sqrt(Q2)
!    print *,"E_t: ",E_t
! We have to calculate the scales:
    call calcmyscales(p)
    call PickSrs(p,ptilde,Bij,Bijkl,radr,radt,Srt,Srtterm)
    call PickCrtStSrt(p,phat,ptilde,Bij, &
                      radr,radt,         &
                      Srt,CrtStSrtterm)
    write(iun,fmt='(a,2(I3,1x))',advance='no') "iexp= ",int(log10(y_rt)),int(log10(E_t))
    write(iun,fmt='(a,G0)') ", CrtStSrt/Srt= ", &
                            CrtStSrtterm/Srtterm
!
  end do
!
end subroutine CrtStSrtLimit
!
! To correctly use this routine a line has to be commented out 
! in MapFSR, above...
subroutine MakeScatterPlot
use process
use particles
use momenta
use collider
use subprocesses
use math
use random
use regions
use FKSutils
use scales
use my_scale
use observables
implicit none
!
!
  integer :: ipoint,istat,ibin,ipart,jpart
  integer :: nleg
  integer :: emit,rad
  real(kind(1d0)) :: rescale
  real(kind(1d0)) :: xi,y,azi
  real(kind(1d0)) :: s_ir,Eg
  real(kind(1d0)) :: sme,I1term,subterm,subterm1,subterm2,ratio
  real(kind(1d0)) :: dx,xpos
  integer , dimension(:) , allocatable :: histo
  real(kind(1d0)) , dimension(:) , allocatable :: histobins
  real(kind(1d0)) , dimension(:,:) , allocatable :: Bij
  real(kind(1d0)) , dimension(:,:) , allocatable :: Rij
  real(kind(1d0)) , dimension(:,:,:) , allocatable :: Bijk
  real(kind(1d0)) , dimension(:,:,:,:) , allocatable :: Bijkl
  real(kind(1d0)) , dimension(:,:,:,:) , allocatable :: Bmunuij
  type(mom) , dimension(:) , allocatable :: pborn
  type(particle) , dimension(:), allocatable :: p,ptilde
!
  character (len=4) , parameter :: cont = 'm+1 '
  character (len=2) , parameter :: subtype = 'C '
  character (len=128) , parameter :: fname = 'scatter.dat'
  integer , parameter :: iun = 99
  integer , parameter :: numpoint = 10000
! for ee -> q q~ g g:
  integer , parameter :: iproc = 4
! C_{gg}: 
  integer , parameter :: iterm = 5
! S_{g_2}:
!  integer , parameter :: iterm = 2
!
  real(kind(1d0)) , parameter :: yval = 1d0 - 1d-8
  real(kind(1d0)) , parameter :: xival = 1d-2
  integer , parameter :: nbin = 50
  real(kind(1d0)) , parameter :: xmin = 0.99d0
  real(kind(1d0)) , parameter :: xmax = 1.01d0
!
  logical , parameter :: check_on_spead = .false.
  real(kind(1d0)) , parameter :: PtMin = 10d0
  real(kind(1d0)) , parameter :: dRMin = 0.4d0
!
!
  interface
    subroutine CalcSubInvariants(p)
    use regions
    use particles
!
      type(particle) , dimension(:) , intent(in) :: p
!
    end subroutine CalcSubInvariants
!
    subroutine PickCirRV(p,ptilde,emit,rad,Cir,Cirterm)
    use regions
    use particles
    implicit none
!
      type(particle), dimension(:), intent(in) :: p
      type(particle), dimension(:), intent(out) :: ptilde
      integer , intent(in) :: emit,rad
      type(subterms) , intent(in) :: Cir
      real(kind(1d0)) , intent(out) :: Cirterm
!
    end subroutine PickCirRV
!
    subroutine PickCirRRA1(p,ptilde,Bij,Bmunuij,emit,rad,Cir,Cirterm)
    use regions
    use particles
    implicit none
!
      type(particle), dimension(:), intent(in) :: p
      type(particle), dimension(:), intent(out) :: ptilde
      real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
      real(kind(1d0)) , dimension(0:,0:,:,:) , intent(inout) :: Bmunuij
      integer , intent(in) :: emit,rad
      type(subterms) , intent(in) :: Cir
      real(kind(1d0)) , intent(out) :: Cirterm
!
    end subroutine PickCirRRA1
!
    subroutine PickSrRV(p,ptilde,Bij,Bijk,rad,Sr,Srterm)
    use regions
    use particles
    implicit none
!
      type(particle), dimension(:), intent(in) :: p
      type(particle), dimension(:), intent(out) :: ptilde
      real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
      real(kind(1d0)) , dimension(:,:,:) , intent(inout) :: Bijk
      integer , intent(in) :: rad
      type(subterms) , intent(in) :: Sr
      real(kind(1d0)) , intent(out) :: Srterm
!
    end subroutine PickSrRV
!
    subroutine PickSrRRA1(p,ptilde,Bij,Bijkl,rad,Sr,Srterm)
    use regions
    use particles
    implicit none
!
      type(particle), dimension(:), intent(in) :: p
      type(particle), dimension(:), intent(out) :: ptilde
      real(kind(1d0)) , dimension(:,:) , allocatable , intent(inout) :: Bij
      real(kind(1d0)) , dimension(:,:,:,:) , allocatable , intent(inout) :: Bijkl
      integer , intent(in) :: rad
      type(subterms) , intent(in) :: Sr
      real(kind(1d0)) , intent(out) :: Srterm
!
    end subroutine PickSrRRA1
!
    subroutine CalcRV(parts,smeRV,smeRVLaurent)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: parts
      real(kind(1d0)) , intent(out) :: smeRV
      real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: &
        smeRVLaurent
!
    end subroutine CalcRV
!
    subroutine CalcR(p,smeR)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: smeR
!
    end subroutine CalcR
!
    subroutine CalcRij(p,Rij)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Rij
!
    end subroutine CalcRij
!
    subroutine CalcI1(p,Bij,CalcSMEB,CalcSMEBij,I1term,I1Laurent)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
      real(kind(1d0)) , intent(out) :: I1term
      real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: I1Laurent
!
      interface
        subroutine CalcSMEB(p,smeB)
        use particles
        implicit none
!
          type(particle) , dimension(:) , intent(in) :: p
          real(kind(1d0)) , intent(out) :: smeB
!
        end subroutine CalcSMEB
!
        subroutine CalcSMEBij(p,Bij)
        use particles
        implicit none
!
          type(particle) , dimension(:) , intent(in) :: p
          real(kind(1d0)) , dimension(:,:) , intent(out) :: Bij
!
        end subroutine CalcSMEBij
      end interface
    end subroutine CalcI1
  end interface
!
  write(*,"(a)") "#####################################################"
  write(*,"(a)") "#####################################################"
  write(*,"(a)") "############## Generating scatter plot ##############"
  write(*,"(a)") "#####################################################"
  write(*,"(a)") "#####################################################"
!
  write(*,"(a,a)") "The scatter plot will be generated for ",cont
  write(*,"(a,a)") "for subtraction type ",subtype
  write(*,"(a)",advance='no') "and subprocess of "
!
  open(unit=iun,file=fname,status='unknown')
!
! Allocations for the histogram:
  allocate(histo(nbin), &
           histobins(0:nbin))
! Filling up x positions for the histogram:
  dx = (xmax - xmin)/nbin
  histo = 0
  do ibin=0,nbin
    histobins(ibin) = xmin + ibin*dx
  end do
!
  if ((cont.eq.'RV  ').or. &
      (cont.eq.'RRA1').or. &
      (cont.eq.'m+1 ')) then
    nleg = nleg_born+1
    call PrintSubProc(flv_ch_Rkin(:,iproc))
  else
    print *,"Wrong contribution type is given..."
    print *,"cont: ",cont
    stop "MakeScatterPlot"
  end if
!
! Allocation of momentum arrays:
  allocate(pborn(nleg-1), &
           p(nleg), &
           ptilde(nleg-1), &
           Bij(nleg-1,nleg-1), &
           Rij(nleg,nleg), &
           Bijk(nleg-1,nleg-1,nleg-1), &
           Bijkl(nleg-1,nleg-1,nleg-1,nleg-1), &
           Bmunuij(0:3,0:3,nleg-1,nleg-1), &
           stat=istat)
  if (istat.ne.0) then
    print *,"Problem with allocation..."
    stop "MakeScatterPlot"
  end if
!
  do ipoint=1,numpoint
!    print *,"ipoint: ",ipoint
    do while (.true.)
      call gen_ranmom(nleg-1,pborn)
      if (.not.check_on_spead) exit
      do ipart=3,nleg-1
! Checking pts:
        if (get_pt(pborn(ipart)).lt.PtMin) goto 100
        do jpart=ipart+1,nleg-1
! Checking separation:
          if (get_dR(pborn(ipart),pborn(jpart)).lt.dRMin) goto 100
        end do
      end do
      exit
100   continue
    end do
! The momenta are needed to be rescaled:
    rescale = ebeam1/pborn(1)%E
    pborn = rescale*pborn
!    call PrintMom(pborn)
    if (cont.eq.'RV  ') then
      if (subtype.eq.'C ') then
        if (ipoint.eq.1) call PrintCir(6,iterm,subterms_Cir_R(iproc)%term(iterm))
        emit = subterms_Cir_R(iproc)%term(iterm)%rad(1)%i
        rad  = subterms_Cir_R(iproc)%term(iterm)%rad(2)%i
        y   = yval
        xi  = gen_rand()
        azi = 2d0*pi*gen_rand()
        call MapFSR(pborn,xi,y,azi,emit,rad,flv_ch_Rkin(:,iproc),p)
        s_ir = 2d0*p(emit)%p*p(rad)%p
!        print *,"sir: ",s_ir
        call CalcSubInvariants(p)
        call calcmyscales(p)
        call CalcRV(p,sme)
        call PickCirRV(p,ptilde,emit,rad,subterms_Cir_R(iproc),subterm)
        ratio = subterm/sme
!        print *,"sme:     ",sme
!        print *,"subterm: ",subterm
!        print *,"ratio:   ",subterm/sme
        do ibin=1,nbin
          if ((ratio.gt.histobins(ibin-1)).and. &
              (ratio.lt.histobins(ibin))) then
            histo(ibin) = histo(ibin) + 1
          end if
        end do
      elseif (subtype.eq.'S ') then
        if (ipoint.eq.1) call PrintSr(6,iterm,subterms_Sr_R(iproc)%term(iterm))
        emit = 3
        rad  = subterms_Sr_R(iproc)%term(iterm)%rad(1)%i
        y   = 0d0
        xi  = xival
        azi = 2d0*pi*gen_rand()
        call MapFSR(pborn,xi,y,azi,emit,rad,flv_ch_Rkin(:,iproc),p)
!        call PrintParts(p)
        Eg = p(rad)%p%E
!        print *,"Eg: ",Eg
        call CalcSubInvariants(p)
        call calcmyscales(p)
        call CalcRV(p,sme)
        call PickSrRV(p,ptilde,Bij,Bijk,rad,subterms_Sr_R(iproc),subterm)
        ratio = subterm/sme
!        print *,"sme:     ",sme
!        print *,"subterm: ",subterm
!        print *,"ratio:   ",subterm/sme
        do ibin=1,nbin
          if ((ratio.gt.histobins(ibin-1)).and. &
              (ratio.lt.histobins(ibin))) then
            histo(ibin) = histo(ibin) + 1
          end if
        end do
      else
        print *,"Wrong subtype value is given..."
        stop "MakeScatterPlot"
      end if
    elseif (cont.eq.'RRA1') then
      if (subtype.eq.'C ') then
        if (ipoint.eq.1) call PrintCir(6,iterm,subterms_Cir_R(iproc)%term(iterm))
        emit = subterms_Cir_R(iproc)%term(iterm)%rad(1)%i
        rad  = subterms_Cir_R(iproc)%term(iterm)%rad(2)%i
        y   = yval
        xi  = gen_rand()
        azi = 2d0*pi*gen_rand()
        call MapFSR(pborn,xi,y,azi,emit,rad,flv_ch_Rkin(:,iproc),p)
        s_ir = 2d0*p(emit)%p*p(rad)%p
!        print *,"sir: ",s_ir
        call CalcSubInvariants(p)
        call calcmyscales(p)
        call CalcI1(p,Rij,CalcR,CalcRij,I1term)
        call PickCirRRA1(p,ptilde,Bij,Bmunuij,emit,rad,subterms_Cir_R(iproc),subterm)
        ratio = subterm/I1term
!        print *,"I1term:  ",I1term
!        print *,"subterm: ",subterm
!        print *,"ratio:   ",subterm/I1term
        do ibin=1,nbin
          if ((ratio.gt.histobins(ibin-1)).and. &
              (ratio.lt.histobins(ibin))) then
            histo(ibin) = histo(ibin) + 1
          end if
        end do
      elseif (subtype.eq.'S ') then
        if (ipoint.eq.1) call PrintSr(6,iterm,subterms_Sr_R(iproc)%term(iterm))
        emit = 3
        rad  = subterms_Sr_R(iproc)%term(iterm)%rad(1)%i
        y   = 0d0
        xi  = xival
        azi = 2d0*pi*gen_rand()
        call MapFSR(pborn,xi,y,azi,emit,rad,flv_ch_Rkin(:,iproc),p)
!        call PrintParts(p)
        Eg = p(rad)%p%E
!        print *,"Eg: ",Eg
        call CalcSubInvariants(p)
        call calcmyscales(p)
        call CalcI1(p,Rij,CalcR,CalcRij,I1term)
        call PickSrRRA1(p,ptilde,Bij,Bijkl,rad,subterms_Sr_R(iproc),subterm)
        ratio = subterm/I1term
!        print *,"I1term:  ",I1term
!        print *,"subterm: ",subterm
!        print *,"ratio:   ",subterm/I1term
        do ibin=1,nbin
          if ((ratio.gt.histobins(ibin-1)).and. &
              (ratio.lt.histobins(ibin))) then
            histo(ibin) = histo(ibin) + 1
          end if
        end do
      else
        print *,"Wrong subtype value is given..."
        stop "MakeScatterPlot"
      end if
    elseif (cont.eq.'m+1 ') then
      if (subtype.eq.'C ') then
        if (ipoint.eq.1) call PrintCir(6,iterm,subterms_Cir_R(iproc)%term(iterm))
        emit = subterms_Cir_R(iproc)%term(iterm)%rad(1)%i
        rad  = subterms_Cir_R(iproc)%term(iterm)%rad(2)%i
        y   = yval
        xi  = gen_rand()
        azi = 2d0*pi*gen_rand()
        call MapFSR(pborn,xi,y,azi,emit,rad,flv_ch_Rkin(:,iproc),p)
        s_ir = 2d0*p(emit)%p*p(rad)%p
!        print *,"sir: ",s_ir
        call CalcSubInvariants(p)
        call calcmyscales(p)
        call CalcRV(p,sme)
        call CalcI1(p,Rij,CalcR,CalcRij,I1term)
        call PickCirRV(p,ptilde,emit,rad,subterms_Cir_R(iproc),subterm1)
        call PickCirRRA1(p,ptilde,Bij,Bmunuij,emit,rad,subterms_Cir_R(iproc),subterm2)
        ratio = (subterm1+subterm2)/(I1term+sme)
!        print *,"I1term+sme:        ",sme+I1term
!        print *,"subterm1+subterm2: ",subterm1+subterm2
!        print *,"ratio:   ",ratio
        do ibin=1,nbin
          if ((ratio.gt.histobins(ibin-1)).and. &
              (ratio.lt.histobins(ibin))) then
            histo(ibin) = histo(ibin) + 1
          end if
        end do
      elseif (subtype.eq.'S ') then
        if (ipoint.eq.1) call PrintSr(6,iterm,subterms_Sr_R(iproc)%term(iterm))
        emit = 3
        rad  = subterms_Sr_R(iproc)%term(iterm)%rad(1)%i
        y   = 0d0
        xi  = xival
        azi = 2d0*pi*gen_rand()
        call MapFSR(pborn,xi,y,azi,emit,rad,flv_ch_Rkin(:,iproc),p)
!        call PrintParts(p)
        Eg = p(rad)%p%E
!        print *,"Eg: ",Eg
        call CalcSubInvariants(p)
        call calcmyscales(p)
        call CalcRV(p,sme)
        call CalcI1(p,Rij,CalcR,CalcRij,I1term)
        call PickSrRV(p,ptilde,Bij,Bijk,rad,subterms_Sr_R(iproc),subterm1)
        call PickSrRRA1(p,ptilde,Bij,Bijkl,rad,subterms_Sr_R(iproc),subterm2)
        ratio = (subterm1+subterm2)/(I1term+sme)
!        print *,"I1term+sme:        ",sme+I1term
!        print *,"subterm1+subterm2: ",subterm1+subterm2
!        print *,"ratio:   ",ratio
        do ibin=1,nbin
          if ((ratio.gt.histobins(ibin-1)).and. &
              (ratio.lt.histobins(ibin))) then
            histo(ibin) = histo(ibin) + 1
          end if
        end do
      else
        print *,"Wrong subtype value is given..."
        stop "MakeScatterPlot"
      end if
    end if
  end do
!
  do ibin=1,nbin
    xpos = (histobins(ibin-1) + histobins(ibin))/2d0
    write(iun,fmt="(E16.6,4x,I0)") xpos,histo(ibin)
  end do
!
  deallocate(pborn,p,ptilde,Bij,Bijk,Bijkl,Bmunuij,histo,histobins)
!
  close(iun)
!
  stop
!
end subroutine MakeScatterPlot
