! This source holds routines to check the pole cancellations
! among various subtractions terms defined for RV and int1RR,A1
subroutine CheckSubtractPoles
use regions
use subprocesses
implicit none
!
!
  integer , parameter :: iun = 99
  integer :: iproc
  character (len=128) :: fname
!
!
  interface
    subroutine CheckSubPoles(iun,num_flv,flv_arr,subtype,terms)
    use regions
    implicit none
!
      integer , intent(in) :: iun
      integer , intent(in) :: num_flv
      integer , dimension(:,:) , intent(in) :: flv_arr
      character (len=2) , intent(in) :: subtype
      type(subterms) , dimension(:) , intent(in) :: terms
!
    end subroutine CheckSubPoles
  end interface
!
  fname = "checkpoles.txt"
  open(unit=iun,file=fname,status='unknown')
  write(iun,'(a)') &
  "********************************************************************"
  write(iun,'(a)') &
  "***************************** Poles ********************************"
  write(iun,'(a)') &
  "********************************************************************"
  write(iun,*) 
  write(iun,'(a)') &
  "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Cir ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
  call CheckSubPoles(iun,num_flv_irr_NLO_R,flv_ch_Rkin, &
                     'C ',subterms_Cir_R)
  write(iun,*) 
  write(iun,'(a)') &
  "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Sr ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
  call CheckSubPoles(iun,num_flv_irr_NLO_R,flv_ch_Rkin, &
                     'S ',subterms_Sr_R)
  write(iun,*) 
  write(iun,'(a)') &
  "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CirSr ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
  call CheckSubPoles(iun,num_flv_irr_NLO_R,flv_ch_Rkin, &
                     'CS',subterms_CSir_R)
!
  close(iun)
!
end subroutine CheckSubtractPoles
!
subroutine CheckSubPoles(iun,num_flv,flv_arr,subtype,terms)
use process
use particles
use regions
use observables
use FKSutils
use math
use random
use my_scale
use scales
use QCDparams
use misc
use nlo_mom_maps
use nnlo_subtractions
implicit none
!
  integer , intent(in) :: iun
  integer , intent(in) :: num_flv
  integer , dimension(:,:) , intent(in) :: flv_arr
  character (len=2) , intent(in) :: subtype
  type(subterms) , dimension(:) , intent(in) :: terms
!
  integer :: iproc,ipoint,iterm,ipart,jpart
  integer :: istat
  integer :: emit,rad
  integer :: ipol
  real(kind(1d0)) :: xi,y,azi
  real(kind(1d0)) :: Cir01,Cir10,Cir00I,CirR00
  real(kind(1d0)) :: Sr01,Sr10,Sr00I,SrR00
  real(kind(1d0)) :: CirSr01,CirSr10,CirSr00I,CirSrR00
  real(kind(1d0)) , dimension(-4:2) :: Cir01Laurent, &
                                       Cir10Laurent, &
                                       Cir00ILaurent, &
                                       CirR00Laurent, &
                                       Sr01Laurent, &
                                       Sr10Laurent, &
                                       Sr00ILaurent, &
                                       SrR00Laurent, &
                                       CirSr01Laurent, &
                                       CirSr10Laurent, &
                                       CirSr00ILaurent, &
                                       CirSrR00Laurent
  real(kind(1d0)) :: smeB,smeV
  real(kind(1d0)) , dimension(0:3,0:3) :: Bmunu
  real(kind(1d0)) , dimension(:,:) , allocatable :: Bij
  real(kind(1d0)) , dimension(:,:) , allocatable :: Vij
  real(kind(1d0)) , dimension(:,:,:) , allocatable :: Bijk
  real(kind(1d0)) , dimension(:,:,:,:) , allocatable :: Bijkl
  real(kind(1d0)) , dimension(:,:,:,:) , allocatable :: Bmunuij
  type(mom) , dimension(:) , allocatable :: pborn
  type(particle) , dimension(:) , allocatable :: p,ptilde
!
  integer , parameter :: numpoint = 2
  real(kind(1d0)) , parameter :: PtMin = 10d0
  real(kind(1d0)) , parameter :: dRMin = 0.4d0
!
!
  interface
    subroutine CalcB(p,smeB)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: smeB
!
    end subroutine CalcB
!
    subroutine CalcV(p,smeV,smeVLaurent)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: smeV
      real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: smeVLaurent
!
    end subroutine CalcV
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
    subroutine CalcVij(p,Vij,VijLaurent)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Vij
      real(kind(1d0)) , optional , dimension(:,:,-4:) , intent(out) :: VijLaurent
!
    end subroutine CalcVij
!
    subroutine CalcBijkl(p,Bijkl)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:,:,:,:) , intent(out) :: Bijkl
!
    end subroutine CalcBijkl
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
    subroutine CalcBmunuij(ileg,p,Bmunuij)
    use particles
    implicit none
!
      integer , intent(in) :: ileg
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(0:,0:,:,:) , intent(out) :: Bmunuij
!
    end subroutine CalcBmunuij
!
    subroutine CalcSubInvariants(p)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
!
    end subroutine CalcSubInvariants
!
    subroutine CalcSubtildeInvariants(p)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
!
    end subroutine CalcSubtildeInvariants
!
    subroutine MapMomCir0FF(p,ptilde,emit,emitID,radi,radiID, &
                            radr,radrID)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      type(particle) , dimension(:) , intent(out) :: ptilde
      integer , intent(in) :: emit,radi,radr
      integer , intent(in) :: emitID,radiID,radrID
!
    end subroutine MapMomCir0FF
  end interface
!
! Allocations:
  allocate(pborn(nleg_born), &
           p(nleg_born+1), &
           ptilde(nleg_born), &
           Bij(nleg_born,nleg_born), &
           Vij(nleg_born,nleg_born), &
           Bijk(nleg_born,nleg_born,nleg_born), &
           Bijkl(nleg_born,nleg_born,nleg_born,nleg_born), &
           Bmunuij(0:3,0:3,nleg_born,nleg_born), &
           stat=istat)
  if (istat.ne.0) then
    print *,"Problem with allocations..."
    stop "CheckSubPoles"
  end if
!
! Generating a common Born PS point which is used to construct
! the real-emission-like one:
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
! Loop over all the subprocesses:
  do iproc=1,num_flv
    call PrintSubProc(iun,flv_arr(:,iproc))
! Loop over subtraction terms:
    do iterm=1,terms(iproc)%numterm
! There is the possibility of checking in multiple PS points:
      do ipoint=1,numpoint
        write(iun,fmt='(a,I0)') "Checking pole cancellation in point ",ipoint
!******************
!***** C-type *****
!******************
        if (subtype.eq.'C ') then
! Determining the emitter and radiated parton:
          emit = terms(iproc)%term(iterm)%rad(1)%i
          rad  = terms(iproc)%term(iterm)%rad(2)%i
! Construction of PS:
          xi  = gen_rand()
          y   = 1d0 - 2d0*gen_rand()
          azi = 2*pi*gen_rand()
          call MapFSR(pborn,xi,y,azi,emit,rad,flv_arr(:,iproc),p)
! Calculating invariants for the subtraction terms:
          call CalcSubInvariants(p)
! The Q2 has to be set up even in the I operator invariants:
          call SetupQ2(Q2)
! Setting up scales:
          call calcmyscales(p)
!          mur = sqrt(2d0*p(1)%p*p(2)%p)
          call PrintCir(iun,iterm,terms(iproc)%term(iterm))
          if (emit.gt.2) then
! Applying collinear momentum mapping:
            call MapMomCir0FF(p,ptilde, &
                              terms(iproc)%term(iterm)%emit(1)%i,  &
                              terms(iproc)%term(iterm)%emit(1)%ID, &
                              terms(iproc)%term(iterm)%rad(1)%i,   &
                              terms(iproc)%term(iterm)%rad(1)%ID,  &
                              terms(iproc)%term(iterm)%rad(2)%i,   &
                              terms(iproc)%term(iterm)%rad(2)%ID)
            call CalcSubtildeInvariants(ptilde)
!
            call CalcB(ptilde,smeB)
            call CalcBij(ptilde,Bij)
            if (terms(iproc)%term(iterm)%emit(1)%ID.eq.0) then
              call CalcBmunu(emit,ptilde,Bmunu)
              call CalcBmunuij(emit,ptilde,Bmunuij)
            end if
!
            call CalcCir01FF(p,ptilde, &
                             terms(iproc)%term(iterm)%emit(1)%i,  &
                             terms(iproc)%term(iterm)%rad(1)%i,   &
                             terms(iproc)%term(iterm)%rad(2)%i,   &
                             Cir01,Cir01Laurent)
            call CalcCir10FF(p,ptilde,smeB,Bmunu, &
                             terms(iproc)%term(iterm)%emit(1)%i,  &
                             terms(iproc)%term(iterm)%rad(1)%i,   &
                             terms(iproc)%term(iterm)%rad(2)%i,   &
                             Cir10,Cir10Laurent)
            call CalcCir00IFF(p,ptilde,smeB,Bij,Bmunu,Bmunuij,     &
                              terms(iproc)%term(iterm)%emit(1)%i,  &
                              terms(iproc)%term(iterm)%rad(1)%i,   &
                              terms(iproc)%term(iterm)%rad(2)%i,   &
                              Cir00I,Cir00ILaurent)
            call CalcCirR00FF(p,ptilde,smeB,Bmunu, &
                              terms(iproc)%term(iterm)%emit(1)%i,  &
                              terms(iproc)%term(iterm)%rad(1)%i,   &
                              terms(iproc)%term(iterm)%rad(2)%i,   &
                              CirR00,CirR00Laurent)
            write(iun,*)
            write(iun,fmt="(a)") "Cancellation for Cir00I + Cir01: "
            write(iun,fmt="(a,a,7x,a)") "             Cir00I    ", &
                                        "                Cir01 ", &
                                        "             norm. sum"
            do ipol=-2,-1
              write(iun,fmt="(a,I0,a)",advance="no") "O(e^",ipol,") : "
              write(iun,*) Cir00ILaurent(ipol),Cir01Laurent(ipol), &
                           (Cir00ILaurent(ipol) + Cir01Laurent(ipol)) &
                           / Cir00ILaurent(ipol)
            end do
            write(iun,*)
            write(iun,fmt="(a)") "Cancellation for CirR00 + Cir10: "
            write(iun,fmt="(a,a,7x,a)") "             CirR00    ", &
                                        "                Cir10 ", &
                                        "             norm. sum"
            do ipol=-2,-1
              write(iun,fmt="(a,I0,a)",advance="no") "O(e^",ipol,") : "
              write(iun,*) CirR00Laurent(ipol),Cir10Laurent(ipol), &
                           (CirR00Laurent(ipol) + Cir10Laurent(ipol)) &
                           / CirR00Laurent(ipol)
            end do
          else
            print *,"Checks for collinear subtractions for ISR are not ready..."
            stop "CheckSubPoles"
          end if
!******************
!***** S-type *****
!******************
        elseif (subtype.eq.'S ') then
! The phase space point is generated according to a FSR fomula,
! hence an emitter has to be selected:
          do emit=3,size(pborn)
! massless parton:
            if (abs(flv_arr(emit,iproc)).gt.qcd_nf) cycle
! does not coincide with the soft gluon:
            if (emit.eq.rad) cycle
            exit
          end do
          rad  = terms(iproc)%term(iterm)%rad(1)%i
! Construction of PS:
          xi  = gen_rand()
          y   = 1d0 - 2d0*gen_rand()
          azi = 2*pi*gen_rand()
          call MapFSR(pborn,xi,y,azi,emit,rad,flv_arr(:,iproc),p)
! Calculating invariants for the subtraction terms:
          call CalcSubInvariants(p)
! The Q2 has to be set up even in the I operator invariants:
          call SetupQ2(Q2)
! Setting up scales:
          call calcmyscales(p)
!          mur = sqrt(2d0*p(1)%p*p(2)%p)
          call PrintSr(iun,iterm,terms(iproc)%term(iterm))
          call MapMomSr(p,ptilde, &
                        terms(iproc)%term(iterm)%rad(1)%i, &
                        terms(iproc)%term(iterm)%rad(1)%ID)
          call CalcSubtildeInvariants(ptilde)
!
          call CalcBij(ptilde,Bij)
          call CalcVij(ptilde,Vij)
          call CalcBijkl(ptilde,Bijkl)
!
          call CalcSr01(p,ptilde,Vij, &
                        terms(iproc)%term(iterm)%rad(1)%i, &
                        Sr01,Sr01Laurent)
          call CalcSr10(p,ptilde,Bij,Bijk, &
                        terms(iproc)%term(iterm)%rad(1)%i, &
                        Sr10,Sr10Laurent)
          call CalcSr00I(p,ptilde,Bij,Bijkl, &
                         terms(iproc)%term(iterm)%rad(1)%i, &
                         Sr00I,Sr00ILaurent)
          call CalcSrR00(p,ptilde,Bij, &
                         terms(iproc)%term(iterm)%rad(1)%i, &
                         SrR00,SrR00Laurent)
          write(iun,*)
          write(iun,fmt="(a)") "Cancellation for Sr00I + Sr01: "
          write(iun,fmt="(a,a,8x,a)") "             Sr00I    ", &
                                      "                 Sr01 ", &
                                      "             norm. sum"
          do ipol=-2,-1
            write(iun,fmt="(a,I0,a)",advance="no") "O(e^",ipol,") : "
            write(iun,*) Sr00ILaurent(ipol),Sr01Laurent(ipol), &
                         (Sr00ILaurent(ipol) + Sr01Laurent(ipol)) &
                         / Sr00ILaurent(ipol)
          end do
          write(iun,*)
          write(iun,fmt="(a)") "Cancellation for SrR00 + Sr10: "
          write(iun,fmt="(a,a,8x,a)") "             SrR00    ", &
                                      "                 Sr10 ", &
                                      "             norm. sum"
          do ipol=-2,-1
            write(iun,fmt="(a,I0,a)",advance="no") "O(e^",ipol,") : "
            write(iun,*) SrR00Laurent(ipol),Sr10Laurent(ipol), &
                         (SrR00Laurent(ipol) + Sr10Laurent(ipol)) &
                         / SrR00Laurent(ipol)
          end do
!*******************
!***** CS-type *****
!*******************
        elseif (subtype.eq.'CS') then
! Determining the emitter and radiated parton:
          emit = terms(iproc)%term(iterm)%rad(1)%i
          rad  = terms(iproc)%term(iterm)%rad(2)%i
! Construction of PS:
          xi  = gen_rand()
          y   = 1d0 - 2d0*gen_rand()
          azi = 2*pi*gen_rand()
          call MapFSR(pborn,xi,y,azi,emit,rad,flv_arr(:,iproc),p)
! Calculating invariants for the subtraction terms:
          call CalcSubInvariants(p)
! The Q2 has to be set up even in the I operator invariants:
          call SetupQ2(Q2)
! Setting up scales:
          call calcmyscales(p)
!          mur = sqrt(2d0*p(1)%p*p(2)%p)
          call PrintCirSr(iun,iterm,terms(iproc)%term(iterm))
! Applying soft momentum mapping:
          call MapMomSr(p,ptilde, &
                        terms(iproc)%term(iterm)%rad(2)%i, &
                        terms(iproc)%term(iterm)%rad(2)%ID)
          call CalcSubtildeInvariants(ptilde)
! The Born SME is calculated for once and for all:
          call CalcBij(ptilde,Bij)
          call CalcV(ptilde,smeV)
          call CastSMEijToSME(ptilde,Bij,smeB)
          if (terms(iproc)%term(iterm)%emit(1)%i.gt.2) then
            call CalcCirFFSr01(p,ptilde,smeV, &
                               terms(iproc)%term(iterm)%emit(1)%i, &
                               terms(iproc)%term(iterm)%rad(1)%i,  &
                               terms(iproc)%term(iterm)%rad(2)%i,  &
                               CirSr01,CirSr01Laurent)
            call CalcCirFFSr10(p,ptilde,smeB, &
                               terms(iproc)%term(iterm)%emit(1)%i, &
                               terms(iproc)%term(iterm)%rad(1)%i,  &
                               terms(iproc)%term(iterm)%rad(2)%i,  &
                               CirSr10,CirSr10Laurent)
            call CalcCirFFSr00I(p,ptilde,smeB,Bij, &
                                terms(iproc)%term(iterm)%emit(1)%i,  &
                                terms(iproc)%term(iterm)%rad(1)%i,   &
                                terms(iproc)%term(iterm)%rad(2)%i,   &
                                CirSr00I,CirSr00ILaurent)
            call CalcCirFFSrR00(p,ptilde,smeB, &
                                terms(iproc)%term(iterm)%emit(1)%i,  &
                                terms(iproc)%term(iterm)%rad(1)%i,   &
                                terms(iproc)%term(iterm)%rad(2)%i,   &
                                CirSrR00,CirSrR00Laurent)
          write(iun,*)
          write(iun,fmt="(a)") "Cancellation for CirSr00I + CirSr01: "
          write(iun,fmt="(a,a,7x,a)") "             CirSr00I  ", &
                                      "                CirSr01", &
                                      "            norm. sum"
          do ipol=-2,-1
            write(iun,fmt="(a,I0,a)",advance="no") "O(e^",ipol,") : "
            write(iun,*) CirSr00ILaurent(ipol),CirSr01Laurent(ipol), &
                         (CirSr00ILaurent(ipol) + CirSr01Laurent(ipol)) &
                         / CirSr00ILaurent(ipol)
          end do
          write(iun,*)
          write(iun,fmt="(a)") "Cancellation for CirSrR00 + CirSr10: "
          write(iun,fmt="(a,a,7x,a)") "             CirSrR00  ", &
                                      "                CirSr10", &
                                      "            norm. sum"
          do ipol=-2,-1
            write(iun,fmt="(a,I0,a)",advance="no") "O(e^",ipol,") : "
            write(iun,*) CirSrR00Laurent(ipol),CirSr10Laurent(ipol), &
                         (CirSrR00Laurent(ipol) + CirSr10Laurent(ipol)) &
                         / CirSrR00Laurent(ipol)
          end do
          else
            print *,"ISR cannot be checked for CS..."
            stop "CheckSubPoles"
          end if
        else
          print *,"Wrong subtraction type is given..."
          print *,"subtype: ",subtype
          stop "CheckSubPoles"
        end if
      end do
    end do
  end do
!
  deallocate(p,pborn,ptilde,Bij,Bijk,Bijkl,Bmunuij)
!
end subroutine CheckSubPoles
