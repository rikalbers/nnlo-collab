! This source contains routines needed for code debugging:
subroutine DebugSME
use process
use momenta
use particles
use utils
use my_scale
use coupling
use math
implicit none
!
!
  integer :: ipart,istat,ipoint
  integer :: nleg
  integer , parameter :: numpoint = 10
  real(kind(1d0)) :: sme,smepr
  real(kind(1d0)) , dimension(-4:2) :: SMELaurent
  real(kind(1d0)) , dimension(5,20) :: p_helac
  integer , dimension(:) , allocatable :: flv
  type(particle) , dimension(:) , allocatable :: p
  character (2) :: cont = 'RV'
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
    subroutine CalcR(parts,smeR)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: parts
      real(kind(1d0)) , intent(out) :: smeR
!
    end subroutine CalcR
!
    subroutine CalcV(parts,smeV,smeVlaurent)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: parts
      real(kind(1d0)) , intent(out) :: smeV
      real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: smeVlaurent
!
    end subroutine CalcV
!
    subroutine CalcRV(parts,smeRV,smeRVlaurent)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: parts
      real(kind(1d0)) , intent(out) :: smeRV
      real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: smeRVlaurent
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
  print *,"############################################################"
  print *,"#                                                          #"
  print *,"#                  Debugging the SMEs...                   #"
  print *,"#                                                          #"
  print *,"############################################################"
!
  if (cont.eq.'B ') then
    print *,"Debugging the Born contributions..."
    nleg = nleg_born
    allocate(flv(nleg),p(nleg),stat=istat)
    if (istat.ne.0) then 
      print *,"Allocation error in DebugSME..." 
      stop
    end if
  end if
!
  if (cont.eq.'R ') then
    print *,"Debugging the Real contributions..."
    nleg = nleg_born+1
    allocate(flv(nleg),p(nleg),stat=istat)
    if (istat.ne.0) then 
      print *,"Allocation error in DebugSME..." 
      stop
    end if
  end if
!
  if (cont.eq.'V ') then
    print *,"Debugging the Virtual contributions..."
    nleg = nleg_born
    allocate(flv(nleg),p(nleg),stat=istat)
    if (istat.ne.0) then 
      print *,"Allocation error in DebugSME..." 
      stop
    end if
  end if
!
  if (cont.eq.'RV') then
    print *,"Debugging the Real-Virtual contributions..."
    nleg = nleg_born + 1
    allocate(flv(nleg),p(nleg),stat=istat)
    if (istat.ne.0) then 
      print *,"Allocation error in DebugSME..." 
      stop
    end if
  end if
!
  if (cont.eq.'RR') then
    print *,"Debugging the RR contributions..."
    nleg = nleg_born+2
    allocate(flv(nleg),p(nleg),stat=istat)
    if (istat.ne.0) then 
      print *,"Allocation error in DebugSME..." 
      stop
    end if
  end if
!
  flv(1) = -11 ; flv(2) =  11
  flv(3) =   2 ; flv(4) =  -2 ; flv(5) =  2
  flv(6) =  -2
!  flv(7) =   0
  print *,"The subprocess is: "
  call PrintSubproc(flv)
!
  do ipoint=1,numpoint
! We are always taking the momenta from fort.301
    read(301,*) p_helac(1:5,1:nleg)
    call CreateParts(p_helac(1:4,1:nleg),flv,p)
    print *,"The read in particles: "
    call PrintParts(p)
! We can have several different contributions:
    if (cont.eq.'B ') then
      call CalcB(p,sme)
      sme = sme / (4d0*pi)**(border_as+border_aEM)
    elseif (cont.eq.'R ') then
      call CalcR(p,sme)
      sme = sme / (4d0*pi)**(border_as+border_aEM+1)
    elseif (cont.eq.'V ') then
! For the virtual part we have to set up the scales too:
      call calcmyscales(p)
      call calc_couplings
! If the full laurent series is requested the full virtual
! is calculated to test whether the combination of the muR
! independent and depend parts returns the full virtual
! we calculate only the finite part, too. Since in that case
! the full virtual is constructed from the muR dependent and
! independent parts.
      call CalcV(p,sme,SMELaurent)
      call CalcV(p,smepr)
      if (abs(sme-smepr).gt.1d-12) then
        print *,"Inconsistency between sme and smepr: "
        print *,"sme:   ",sme
        print *,"smepr: ",smepr
      end if
      sme = sme / (4d0*pi)**(border_as+border_aEM+1)
      SMELaurent = SMELaurent / (4d0*pi)**(border_as+border_aEM+1)
      write(404,*) SMELaurent(-2)
      write(405,*) SMELaurent(-1)
    elseif (cont.eq.'RV') then
! For the real-virtual part we have to set up the scales too:
      call calcmyscales(p)
      call calc_couplings
! For more details why we are calling the real-virtual routine
! twice take a look at the part written for the virtual part...
      call CalcRV(p,sme,SMELaurent)
      call CalcRV(p,smepr)
      if (abs(sme-smepr).gt.1d-12) then
        print *,"Inconsistency between sme and smepr: "
        print *,"sme:   ",sme
        print *,"smepr: ",smepr
      end if
      sme = sme / (4d0*pi)**(border_as+border_aEM+2)
      SMELaurent = SMELaurent / (4d0*pi)**(border_as+border_aEM+2)
      write(404,*) SMELaurent(-2)
      write(405,*) SMELaurent(-1)
    elseif (cont.eq.'RR') then
      call CalcRR(p,sme)
      sme = sme / (4d0*pi)**(border_as+border_aEM+2)
    else
      print *,"Unknown type of contribution..."
      print *,cont
      stop
    end if
! We write out the result:
! The result is appearing in fort.401
    write(401,*) sme
  end do
!
  deallocate(flv,p)
  stop
!
end subroutine DebugSME
!
! This routine manages the debug of the spin-correlated 
! modulus squared amplitudes:
subroutine DebugSCSME
use process
use particles
use momenta
use subprocesses
use misc
implicit none
!
!
  integer :: ipart,istat,ipoint,iproc,ileg
  integer :: nleg
  integer , parameter :: numpoint = 2
  type(particle) , dimension(:) , allocatable :: p
  character (2) :: cont = 'RV'
  real(kind(1d0)) :: Born,Bornpr
  real(kind(1d0)) , dimension(0:3,0:3) :: Bmunu
  real(kind(1d0)) :: Virt,Virtpr
  real(kind(1d0)) , dimension(0:3,0:3) :: Vmunu
  real(kind(1d0)) :: Real,Realpr
  real(kind(1d0)) , dimension(0:3,0:3) :: Rmunu
!
!
  interface
    subroutine CalcB(p,Born)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: Born
!
    end subroutine calcB
!
    subroutine CalcBmunu(ileg,p,Bmunu)
    use particles
    implicit none
!
      integer , intent(in) :: ileg
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(0:3,0:3) , intent(out) :: Bmunu
!
    end subroutine calcBmunu
!
    subroutine CalcV(p,Virt,VirtLaurent)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: Virt
      real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: VirtLaurent
!
    end subroutine calcV
!
    subroutine CalcVmunu(ileg,p,Vmunu)
    use particles
    implicit none
!
      integer , intent(in) :: ileg
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(0:3,0:3) , intent(out) :: Vmunu
!
    end subroutine calcVmunu
!
    subroutine CalcR(p,Real)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: Real
!
    end subroutine calcR
!
    subroutine CalcRmunu(ileg,p,Rmunu)
    use particles
    implicit none
!
      integer , intent(in) :: ileg
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(0:3,0:3) , intent(out) :: Rmunu
!
    end subroutine calcRmunu
  end interface
!
  print *,"############################################################"
  print *,"#                                                          #"
  print *,"#                       Debugging the                      #"
  print *,"#                  Spin-correlated SMEs....                #"
  print *,"#                                                          #"
  print *,"############################################################"
!
  if (cont.eq.'R ') nleg = nleg_born
  if (cont.eq.'RV') nleg = nleg_born
  if (cont.eq.'RR') nleg = nleg_born+1
!
! We allocate array an array to hold the particles and:
  allocate(p(nleg),stat=istat)
  if (istat.ne.0) then
    print *,"Problem during allocation of p in DebugSCSME..."
    stop
  end if
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! NLO R - A1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  if (cont.eq.'R ') then
! We go through all the subprocesses:
! We only have to go through the independent ones:
    do iproc=1,num_flv_irr_LO
      call PrintSubProc(flv_ch_Bkin(:,iproc))
! We check the SCSME in several PS points:
      do ipoint=1,numpoint
        write(*,'(a,I0)') "Checking the SCSME in point ",ipoint
! We pick a randomly chosen PS point:
        call gen_ranparts(nleg,flv_ch_Bkin(:,iproc),p)
        call PrintParts(p)
! We calculate the spin-correlated SME for legs having a gluon:
        do ileg=1,nleg
          if (p(ileg)%flv.eq.0) then
            write(*,'(a,I0)') "Gluon found at leg ",ileg
! We calculate Bmunu:
            call CalcBmunu(ileg,p,Bmunu)
            call CalcB(p,Born)
! We calculate the Born from Bmunu, this is done by a trivial
! contraction with d_{\mu\nu}, since the underlying on-shell
! currents are conserved we only have to contract with -g_{\mu\nu}: 
            Bornpr = -Bmunu(0,0) + Bmunu(1,1) + Bmunu(2,2) + Bmunu(3,3)
            write(*,*) "The Born is:            ",Born
            write(*,*) "The Born from Bmunu is: ",Bornpr
            write(*,*) "The ratio is: ",Bornpr/Born
            write(*,*) "p.Bmunu.p: ",Contract(p(ileg),Bmunu,p(ileg))
            write(*,*) "p.Bmunu.p: ",p(ileg)%p*Bmunu*p(ileg)%p
          end if
        end do
      end do
    end do
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! NNLO RR - A1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  elseif (cont.eq.'RR') then
! We go through all the subprocesses:
! We only have to go through the independent ones:
    do iproc=1,num_flv_irr_NLO_R
      call PrintSubProc(flv_ch_Rkin(:,iproc))
! We check the SCSME in several PS points:
      do ipoint=1,numpoint
        write(*,'(a,I0)') "Checking the SCSME in point ",ipoint
! We pick a randomly chosen PS point:
        call gen_ranparts(nleg,flv_ch_Rkin(:,iproc),p)
        call PrintParts(p)
! We calculate the spin-correlated SME for legs having a gluon:
        do ileg=1,nleg
          if (p(ileg)%flv.eq.0) then
            write(*,'(a,I0)') "Gluon found at leg ",ileg
! We calculate Bmunu:
            call CalcRmunu(ileg,p,Rmunu)
            call CalcR(p,Real)
! We calculate the Born from Bmunu, this is done by a trivial
! contraction with d_{\mu\nu}, since the underlying on-shell
! currents are conserved we only have to contract with -g_{\mu\nu}: 
            Realpr = -Rmunu(0,0) + Rmunu(1,1) + Rmunu(2,2) + Rmunu(3,3)
            write(*,*) "The Real is:            ",Real
            write(*,*) "The Real from Rmunu is: ",Realpr
            write(*,*) "The ratio is: ",Realpr/Real
            write(*,*) "p.Rmunu.p: ",Contract(p(ileg),Rmunu,p(ileg))
            write(*,*) "p.Rmunu.p: ",p(ileg)%p*Rmunu*p(ileg)%p
          end if
        end do
      end do
    end do
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! NNLO RV - A1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  elseif (cont.eq.'RV') then
! We go through all the subprocesses:
! We only have to go through the independent ones:
    do iproc=1,num_flv_irr_NLO_V
      call PrintSubProc(flv_ch_Bkin(:,iproc))
! We check the SCSME in several PS points:
      do ipoint=1,numpoint
        write(*,'(a,I0)') "Checking the SCSME in point ",ipoint
! We pick a randomly chosen PS point:
        call gen_ranparts(nleg,flv_ch_Bkin(:,iproc),p)
        call PrintParts(p)
! We calculate the spin-correlated SME for legs having a gluon:
        do ileg=1,nleg
          if (p(ileg)%flv.eq.0) then
            write(*,'(a,I0)') "Gluon found at leg ",ileg
! We calculate Bmunu:
            call CalcVmunu(ileg,p,Vmunu)
            call CalcV(p,Virt)
! We calculate the Virt from Vmunu, this is done by a trivial
! contraction with d_{\mu\nu}, since the underlying on-shell
! currents are conserved we only have to contract with -g_{\mu\nu}: 
            Virtpr = -Vmunu(0,0) + Vmunu(1,1) + Vmunu(2,2) + Vmunu(3,3)
            write(*,*) "The Virt is:            ",Virt
            write(*,*) "The Virt from Vmunu is: ",Virtpr
            write(*,*) "The ratio is: ",Virtpr/Virt
            write(*,*) "p.Vmunu.p: ",Contract(p(ileg),Vmunu,p(ileg))
            write(*,*) "p.Vmunu.p: ",p(ileg)%p*Vmunu*p(ileg)%p
          end if
        end do
      end do
    end do
  end if
!
  stop
!
end subroutine DebugSCSME
!
! This routine manages the debug of the color-correlated 
! modulus squared amplitudes:
subroutine DebugCCSME
use process
use particles
use momenta
use subprocesses
use misc
use utils
use QCDparams
implicit none
!
!
  integer :: ipart,istat,ipoint,iproc,ileg,jleg
  integer :: nleg
  integer , parameter :: numpoint = 2
  type(particle) , dimension(:) , allocatable :: p
  character (2) :: cont = 'RR'
  real(kind(1d0)) :: Born,CBorn,CBornpr
  real(kind(1d0)) :: Real,CReal,CRealpr
  real(kind(1d0)) , dimension(:,:) , allocatable :: Bij
  real(kind(1d0)) , dimension(:,:) , allocatable :: Rij
!
!
  interface
    subroutine CalcB(p,Born)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: Born
!
    end subroutine calcB
!
    subroutine CalcBij(p,Bij)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Bij
!
    end subroutine calcBij
!
    subroutine CalcR(p,Real)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: Real
!
    end subroutine calcR
!
    subroutine CalcRij(p,Rij)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Rij
!
    end subroutine calcRij
  end interface
!
  print *,"############################################################"
  print *,"#                                                          #"
  print *,"#                       Debugging the                      #"
  print *,"#                   Color-correlated SMEs....              #"
  print *,"#                                                          #"
  print *,"############################################################"
!
  if (cont.eq.'R ') nleg = nleg_born
  if (cont.eq.'RV') nleg = nleg_born
  if (cont.eq.'RR') nleg = nleg_born+1
!
! We allocate array an array to hold the particles and color correlated
! SMEs:
  allocate(p(nleg), &
           Bij(nleg,nleg), &
           Rij(nleg,nleg), &
           stat=istat)
  if (istat.ne.0) then
    print *,"Problem during allocation of p and Bij,Rij in DebugCCSME..."
    stop
  end if
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! NLO R - A1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  if (cont.eq.'R ') then
! We go through all the subprocesses:
! We only have to go through the independent ones:
    do iproc=1,num_flv_irr_LO
      call PrintSubProc(flv_ch_Bkin(:,iproc))
! We check the CCSME in several PS points:
      do ipoint=1,numpoint
        write(*,'(a,I0)') "Checking the CCSME in point ",ipoint
! We pick a randomly chosen PS point:
        call gen_ranparts(nleg,flv_ch_Bkin(:,iproc),p)
        call PrintParts(p)
        call CalcB(p,Born)
        call CalcBij(p,Bij)
! Loop over all partons:
        do ileg=1,nleg
          if (abs(flv_ch_Bkin(ileg,iproc)).gt.6) cycle
          write(*,fmt='(a,I0,a,a)',advance='no') "ileg: ",ileg," , ", &
            ConvertFromPDG(flv_ch_Bkin(ileg,iproc))
          if (flv_ch_Bkin(ileg,iproc).ne.0) CBorn = qcd_cf*Born
          if (flv_ch_Bkin(ileg,iproc).eq.0) CBorn = qcd_ca*Born
          CBornpr = 0d0
          do jleg=1,nleg
            if (jleg.eq.ileg) cycle
            CBornpr = CBornpr - Bij(ileg,jleg) - Bij(jleg,ileg)
          end do
          CBornpr = CBornpr/2d0
          write(*,*) "C * Born: ",CBorn," C * Born': ",CBornpr
        end do
      end do
    end do
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! NNLO RR - A1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  elseif (cont.eq.'RR') then
! We go through all the subprocesses:
! We only have to go through the independent ones:
    do iproc=1,num_flv_irr_NLO_R
      call PrintSubProc(flv_ch_Rkin(:,iproc))
! We check the CCSME in several PS points:
      do ipoint=1,numpoint
        write(*,'(a,I0)') "Checking the CCSME in point ",ipoint
! We pick a randomly chosen PS point:
        call gen_ranparts(nleg,flv_ch_Rkin(:,iproc),p)
        call PrintParts(p)
        call CalcR(p,Real)
        call CalcRij(p,Rij)
! Loop over all partons:
        do ileg=1,nleg
          if (abs(flv_ch_Rkin(ileg,iproc)).gt.6) cycle
          write(*,fmt='(a,I0,a,a)',advance='no') "ileg: ",ileg," , ", &
            ConvertFromPDG(flv_ch_Rkin(ileg,iproc))
          if (flv_ch_Rkin(ileg,iproc).ne.0) CReal = qcd_cf*Real
          if (flv_ch_Rkin(ileg,iproc).eq.0) CReal = qcd_ca*Real
          CRealpr = 0d0
          do jleg=1,nleg
            if (jleg.eq.ileg) cycle
            CRealpr = CRealpr - Rij(ileg,jleg) - Rij(jleg,ileg)
          end do
          CRealpr = CRealpr/2d0
          write(*,*) "C * Real: ",CReal," C * Real': ",CRealpr
          print *,"ratio: ",CRealpr/CReal
        end do
      end do
    end do
  else
    print *,"Error occured in DebugCCSME..." 
    print *,"There is no contribution like: ",cont
    stop
  end if
!
  stop
!
end subroutine DebugCCSME
!
! This routine helps debugging the I1 type integrated subtraction
! terms by evaluating them up to O(eps^-2) and calculating the
! corresponding virtual part also up to O(eps^-2) in order to see
! the cancelation.
subroutine DebugI1op
use process
use flags
use subprocesses
use particles
use my_scale
use coupling
use scales
implicit none
!
!
  integer :: ipart,istat,ipoint,iproc,ipol
  integer :: nleg
  integer , parameter :: numpoint = 2
  type(particle) , dimension(:) , allocatable :: p
  character (2) :: cont = 'V '
  real(kind(1d0)) :: I1term,sme
  real(kind(1d0)) , dimension(-4:2) :: I1Laurent
  real(kind(1d0)) , dimension(-4:2) :: SMELaurent
  real(kind(1d0)) , dimension(:,:) , allocatable :: Bij
!
!
  interface
    subroutine CalcB(p,Born)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: Born
!
    end subroutine calcB
!
    subroutine CalcR(p,Real)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: Real
!
    end subroutine calcR
!
    subroutine CalcBij(p,Bij)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Bij
!
    end subroutine calcBij
!
    subroutine CalcRij(p,Rij)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Rij
!
    end subroutine calcRij
!
    subroutine CalcV(parts,smeV,smeVlaurent)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: parts
      real(kind(1d0)) , intent(out) :: smeV
      real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: smeVlaurent
!
    end subroutine CalcV
!
    subroutine CalcRV(parts,smeRV,smeRVlaurent)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: parts
      real(kind(1d0)) , intent(out) :: smeRV
      real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: smeRVlaurent
!
    end subroutine CalcRV
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
  print *,"############################################################"
  print *,"#                                                          #"
  print *,"#                       Debugging the                      #"
  print *,"#                I1-type insertion operators...            #"
  print *,"#                                                          #"
  print *,"############################################################"
!
! The I1 operator can turn up at different places:
  if ((cont.eq.'V ').and.flg_NLO_V) then 
! Allocating arrays holding the momenta and the color-correlated MEs:
    allocate(p(nleg_born), &
             Bij(nleg_born,nleg_born), &
             stat=istat)
    if (istat.ne.0) then
      print *,"Problem with allocation of p and Bij in DebugI1op..."
      stop
    end if
! We would like to see the cancellation of poles for each and
! every irreducible subprocess, loop over them:
    do iproc=1,num_flv_irr_LO
      do ipoint=1,numpoint
        write(*,'(a,I0)') "Checking the I1 operator in point ",ipoint
! We pick a randomly chosen PS point:
        call gen_ranparts(nleg_born,flv_ch_Bkin(:,iproc),p)
        call PrintParts(p)
! We have to set up the scales:
        call calcmyscales(p)
        call calc_couplings
        call CalcI1(p,Bij,CalcB,CalcBij,I1term,I1Laurent)
        call CalcV(p,sme,SMELaurent)
        write(*,*) "                      Virtual", &
                   "                      I1op", &
                   "                  norm. sum."
        do ipol=-4,-1
          if (SMELaurent(ipol).eq.0d0) cycle
          write(*,'(a,I2,a)',advance='no') "O(eps^(",ipol,"))"
          write(*,*) SMELaurent(ipol),I1Laurent(ipol), &
            (SMELaurent(ipol) + I1Laurent(ipol)) / &
            SMELaurent(ipol)
        end do
      end do
    end do
  elseif ((cont.eq.'RV').and.flg_NNLO_RV) then 
! Allocating arrays holding the momenta and the color-correlated MEs:
    allocate(p(nleg_born+1), &
             Bij(nleg_born+1,nleg_born+1), &
             stat=istat)
    if (istat.ne.0) then
      print *,"Problem with allocation of p and Bij in DebugI1op..."
      stop
    end if
! We would like to see the cancellation of poles for each and
! every irreducible subprocess, loop over them:
    do iproc=1,num_flv_irr_NLO_R
      do ipoint=1,numpoint
        write(*,'(a,I0)') "Checking the I1 operator in point ",ipoint
! We pick a randomly chosen PS point:
        call gen_ranparts(nleg_born+1,flv_ch_Rkin(:,iproc),p)
        call PrintParts(p)
! We have to set up the scales:
        call calcmyscales(p)
        call calc_couplings
        call CalcI1(p,Bij,CalcR,CalcRij,I1term,I1Laurent)
        call CalcRV(p,sme,SMELaurent)
        write(*,*) "                      Virtual", &
                   "                      I1op", &
                   "                  norm. sum."
        do ipol=-4,-1
          if (SMELaurent(ipol).eq.0d0) cycle
          write(*,'(a,I2,a)',advance='no') "O(eps^(",ipol,"))"
          write(*,*) SMELaurent(ipol),I1Laurent(ipol), &
            (SMELaurent(ipol) + I1Laurent(ipol)) / &
            SMELaurent(ipol)
        end do
      end do
    end do
  else
    print *,"I1 check is not implemented for ",cont
    stop
  end if
!
  stop
!
end subroutine DebugI1op
!
! This routine tests the simultaneously spin and color correlated SME:
! SCCSME = Spin- and Color-Correlated SME:
subroutine DebugSCCSME
use utils
use process
use particles
use momenta
use subprocesses
use misc
implicit none
!
!
  integer :: ipart,istat,ipoint,iproc,ileg
  integer :: i,j,mu,nu
  integer , parameter :: numpoint = 2
  type(particle) , dimension(:) , allocatable :: p
  real(kind(1d0)) :: Born,Bornpr
  real(kind(1d0)) :: diff
  real(kind(1d0)) , dimension(0:3,0:3) :: Bmunu,Bmunupr
  real(kind(1d0)) , dimension(:,:) , allocatable :: Bij,Bijpr
  real(kind(1d0)) , dimension(:,:,:,:) , allocatable :: Bmunuij
!
  real(kind(1d0)) , parameter :: eps = 1d-12
!
!
  interface
    subroutine CalcB(p,Born)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: Born
!
    end subroutine calcB
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
    subroutine CalcBij(p,Bij)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Bij
!
    end subroutine CalcBij
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
  end interface
!
  print *,"############################################################"
  print *,"#                                                          #"
  print *,"#                       Debugging the                      #"
  print *,"#            Spin- and Color correlated SMEs....           #"
  print *,"#                                                          #"
  print *,"############################################################"
!
! We allocate array an array to hold the particles and the various
! spin- and/or color-correlated SMEs:
  allocate(p(nleg_born), &
           Bij(nleg_born,nleg_born), &
           Bijpr(nleg_born,nleg_born), &
           Bmunuij(0:3,0:3,nleg_born,nleg_born), &
           stat=istat)
  if (istat.ne.0) then
    print *,"Problem with allocation of arrays..."
    stop "DebugSCCSME"
  end if
! We go through all the subprocesses:
! We only have to go through the independent ones:
  do iproc=1,num_flv_irr_LO
    call PrintSubProc(flv_ch_Bkin(:,iproc))
! We check the SCCSME in several PS points:
    do ipoint=1,numpoint
      write(*,'(a,I0)') "Checking the SCCSME in point ",ipoint
! We pick a randomly chosen PS point:
      call gen_ranparts(nleg_born,flv_ch_Bkin(:,iproc),p)
      call PrintParts(p)
! Loop over all particles selecting only gluons:
      do ileg=1,nleg_born
        if (abs(flv_ch_Bkin(ileg,iproc)).ne.0) cycle
        write(*,fmt='(a,I0,a,a)') "ileg: ",ileg," , ", &
          ConvertFromPDG(flv_ch_Bkin(ileg,iproc))
! Obtaining the simultaneously spin- and color-correlated SME:
        call CalcBmunuij(ileg,p,Bmunuij)
! Obtaining the spin-correlated SME:
        call CalcBmunu(ileg,p,Bmunu)
! Obtaining the color-correlated SME:
        call CalcBij(p,Bij)
! Obtaining the SME itself:
        call CalcB(p,Born)
! Constructing the spin-correlated SME from the SC-correlated SME:
        call CastSMEmunuijToSMEmunu(.true.,p,Bmunuij,Bmunupr)
! To compare we make the difference then add all the entries:
        diff = abs(sum(Bmunu - Bmunupr))
        if (diff.gt.eps) then
          print *,"Check failed with diff = ",diff
          print *,"Bmunu: "
          do mu=0,3
            print *,(Bmunu(mu,nu) , nu=0,3)
          end do
          print *,"Bmunupr: "
          do mu=0,3
            print *,(Bmunupr(mu,nu) , nu=0,3)
          end do
        else
          print *,"Bmunuij -> Bmunu check passed with diff = ",diff
        end if
! Constructing the color-correlated SME from the SC-correlated one:
        call CastSMEmunuijToSMEij(p,Bmunuij,Bijpr)
! To compare we make the difference then add all the entries:
        diff = abs(sum(Bij - Bijpr))
        if (diff.gt.eps) then
          print *,"Check failed with diff = ",diff
          print *,"Bij: "
          do i=1,nleg_born
            print *,(Bij(i,j) , j=1,nleg_born)
          end do
          print *,"Bijpr: "
          do i=1,nleg_born
            print *,(Bijpr(i,j) , j=1,nleg_born)
          end do
        else
          print *,"Bmunuij -> Bij check passed with diff = ",diff
        end if
! Constructing the SME from the SC-correlated one through the 
! spin-correlated intermediate result:
        call CastSMEmunuToSME(p,Bmunupr,Bornpr)
        diff = abs(Born - Bornpr)
        if (diff.gt.eps) then
          print *,"Check failed with diff = ",diff
          print *,"Born: ",Born
          print *,"Bijpr: ",Bornpr
        else
          print *,"Bmunuij -> Bmunu -> B check passed with diff = ",diff
        end if
! Constructing the SME from the SC-correlated one through the 
! color-correlated intermediate result:
        call CastSMEijToSME(p,Bijpr,Bornpr)
        diff = abs(Born - Bornpr)
        if (diff.gt.eps) then
          print *,"Check failed with diff = ",diff
          print *,"Born: ",Born
          print *,"Bijpr: ",Bornpr
        else
          print *,"Bmunuij -> Bij -> B check passed with diff = ",diff
        end if
      end do
    end do
  end do
!
  deallocate(p,Bij,Bijpr,Bmunuij)
!
  stop
!
end subroutine DebugSCCSME
!
! This routine tests the doubly color correlated SME:
! DCCSME = Doubly Color-Correlated SME:
subroutine DebugDCCSME
use utils
use process
use particles
use momenta
use subprocesses
use misc
implicit none
!
!
  integer :: ipart,istat,ipoint,iproc,ileg
  integer :: i,j,mu,nu
  integer , parameter :: numpoint = 2
  type(particle) , dimension(:) , allocatable :: p
  real(kind(1d0)) :: Born,Bornpr
  real(kind(1d0)) :: diff
  real(kind(1d0)) , dimension(:,:) , allocatable :: Bij,Bijpr
  real(kind(1d0)) , dimension(:,:,:,:) , allocatable :: Bijkl
!
  real(kind(1d0)) , parameter :: eps = 1d-12
!
!
  interface
    subroutine CalcB(p,Born)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: Born
!
    end subroutine calcB
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
    subroutine CalcBijkl(p,Bijkl)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:,:,:,:) , intent(out) :: Bijkl
!
    end subroutine CalcBijkl
  end interface
!
  print *,"############################################################"
  print *,"#                                                          #"
  print *,"#                       Debugging the                      #"
  print *,"#               Doubly Color correlated SMEs....           #"
  print *,"#                                                          #"
  print *,"############################################################"
!
! We allocate array an array to hold the particles and the various
! color-correlated SMEs:
  allocate(p(nleg_born), &
           Bij(nleg_born,nleg_born), &
           Bijpr(nleg_born,nleg_born), &
           Bijkl(nleg_born,nleg_born,nleg_born,nleg_born), &
           stat=istat)
  if (istat.ne.0) then
    print *,"Problem with allocation of arrays..."
    stop "DebugDCCSME"
  end if
! We go through all the subprocesses:
! We only have to go through the independent ones:
  do iproc=1,num_flv_irr_LO
    call PrintSubProc(flv_ch_Bkin(:,iproc))
! We check the DCCSME in several PS points:
    do ipoint=1,numpoint
      write(*,'(a,I0)') "Checking the DCCSME in point ",ipoint
! We pick a randomly chosen PS point:
      call gen_ranparts(nleg_born,flv_ch_Bkin(:,iproc),p)
      call PrintParts(p)
! Obtaining the doubly color-correlated SME:
      call CalcBijkl(p,Bijkl)
! Obtaining the color-correlated SME:
      call CalcBij(p,Bij)
! Obtaining the SME itself:
      call CalcB(p,Born)
! Constructing the color-correlated one:
      call CastSMEijklToSMEij(.true.,p,Bijkl,Bijpr)
! To compare we make the difference then add all the entries:
      diff = abs(sum(Bij - Bijpr))
      if (diff.gt.eps) then
        print *,"Check failed with diff = ",diff
        print *,"Bij: "
        do i=1,nleg_born
          print *,(Bij(i,j) , j=1,nleg_born)
        end do
        print *,"Bijpr: "
        do i=1,nleg_born
          print *,(Bijpr(i,j) , j=1,nleg_born)
        end do
      else
        print *,"Bijkl -> Bij check passed with diff = ",diff
      end if
! Constructing the SME out of the doubly color-correlated SME in an
! indirect way:
      call CastSMEijToSME(p,Bijpr,Bornpr)
      diff = abs(Born - Bornpr)
      if (diff.gt.eps) then
        print *,"Check failed with diff = ",diff
        print *,"Born: ",Born
        print *,"Bornpr: ",Bornpr
      else
        print *,"Bijkl -> Bij -> B check passed with diff = ",diff
      end if
    end do
  end do
!
  deallocate(p,Bij,Bijpr,Bijkl)
!
  stop
!
end subroutine DebugDCCSME
!
subroutine DebugRRealPS(p)
use momenta
implicit none
!
  type(mom) , dimension(:) , intent(out) :: p
!
  real(kind(1d0)) :: Q2
  type(mom) :: Q
!
!
  p(3) = (/22.83834791233909d0,-3.470281096743738d0,13.139110305162939d0,26.57572587657201d0/)
  p(4) = (/-1.4218257356086295d0,15.199285354437809d0,1.424200979175682d0,15.331934389931087d0/) 
  p(5) = (/-29.69802717349611d0,-7.550456122913248d0,-9.37501943394676d0,32.04486222545873d0/)
  p(6) = (/6.626183329263178d0,-3.3476045509362873d0,-4.138399523745368d0,8.499359526493194d0/)
  p(7) = (/1.6553216675024718d0,-0.8309435838445329d0,-1.049892326646499d0,2.129044612023611d0/)
!
  Q = p(3) + p(4) + p(5) + p(6) + p(7)
  Q2 = Q**2
!
  p(1) = 0d0
  p(2) = 0d0
  p(1)%pz = sqrt(Q2)/2d0
  p(2)%pz = -sqrt(Q2)/2d0
  p(1)%E = sqrt(Q2)/2d0
  p(2)%E = sqrt(Q2)/2d0
!
end subroutine DebugRRealPS
