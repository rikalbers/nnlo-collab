! This file holds the user analysis routine and related functions
! and subroutines:
module analysis_supp
implicit none
!
  integer , parameter :: nycut = 15
  real(kind(1d0)) , dimension(nycut) :: ycut = &
    (/0.3d0,0.15d0,0.1d0,0.06d0,0.03d0,0.015d0,0.01d0,0.006d0,0.003d0, &
      0.0015d0,0.001d0,0.0006d0,0.0003d0,0.00015d0,0.0001d0/)
  integer , parameter :: njetalgo = 5
  character (len=10), dimension(njetalgo) :: cjetalgo = &
    (/"Durham3jet","Geneva3jet","JadeE03jet","JadeE3jet ","Cmbrdg3jet"/)
  character (len=1) , dimension(9) :: cn = &
    (/'1','2','3','4','5','6','7','8','9'/)
!
end module analysis_supp
!
! This analysis corresponds to the one published in arXiv:0904.1077:
!
subroutine init_analysis
use histo
use analysis_supp
implicit none
!
!
  integer :: ialgo
!
!
! We have to initialize the histograms, we allow for 100 bins in each:
  call init_hist(100)
!
! 1 - Thrust:
  call bookup_hist("1-T",0.01d0,0d0,0.45d0)
! Heavy jet mass:
  call bookup_hist("rho",0.01d0,0d0,0.43d0)
! Wide jet broadenning:
  call bookup_hist("BW",0.01d0,0d0,0.34d0)
! Total jet broadenning:
  call bookup_hist("BT",0.01d0,0d0,0.42d0)
! C-parameter:
  call bookup_hist("Cpar",0.01d0,0d0,1d0)
! Three-to-two jet transition:
  call bookup_hist("y23",0.25d0,-10.0d0,-0.75d0)
! Jet-algos:
  do ialgo=1,njetalgo
    call bookup_hist(cjetalgo(ialgo),1d0,0.5d0,15.5d0)
  end do
!
end subroutine init_analysis
!
subroutine analysis(p,wgt)
use histo
use particles
use observables
use analysis_supp
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
  real(kind(1d0)) , dimension(:) , intent(in) :: wgt
!
  integer :: ipart,iycut,ialgo
  integer :: npart,ntrack,njet
  integer , parameter :: maxjet = 5
  real(kind(1d0)) :: Q2
  real(kind(1d0)) , dimension(maxjet) :: YnArr
  real(kind(1d0)) , dimension(4,maxjet) :: ptrack
  real(kind(1d0)) , dimension(20) :: EvntShapes
  real(kind(1d0)) , dimension(4,3) :: ThrustDirs
  integer :: err
!
  real(kind(1d0)) :: T,Tmaj,Tmin,Oblat,Cpar,Dpar,MH2,md2,BW,BT, &
                     y23,logy23,rho
!
!
  ntrack = 0
  njet   = 0
!
  Q2 = (p(1)%p + p(2)%p)**2
!
  npart = size(p)
! Constructing the tracks:
  do ipart=3,npart
    if (abs(p(ipart)%flv).gt.5) cycle
    ntrack = ntrack + 1
    ptrack(1,ntrack) = p(ipart)%p%px
    ptrack(2,ntrack) = p(ipart)%p%py
    ptrack(3,ntrack) = p(ipart)%p%pz
    ptrack(4,ntrack) = p(ipart)%p%E
  end do
!
  call EvShapes(ptrack,ntrack,.false.,EvntShapes,ThrustDirs,err)
!
  T     = EvntShapes(1)
  Tmaj  = EvntShapes(2)
  Tmin  = EvntShapes(3)
  Oblat = EvntShapes(4)
  Cpar  = EvntShapes(5)
  MH2   = EvntShapes(6)
  md2   = EvntShapes(7)
  BT    = EvntShapes(8)
  BW    = EvntShapes(9)
  y23   = EvntShapes(10)
  Dpar  = EvntShapes(11)
!
  rho   = MH2
!
  logy23 = log(y23)
!
! 1 - Thrust:
  call fill_hist("1-T",1d0 - T,(1-T)*wgt)
! Heavy jet mass:
  call fill_hist("rho",rho,rho*wgt)
! Wide jet broadening:
  call fill_hist("BW",BW,BW*wgt)
! Total jet broadening:
  call fill_hist("BT",BT,BT*wgt)
! C parameter:
  call fill_hist("Cpar",Cpar,Cpar*wgt)
! Three-to-two jet transition:
  call fill_hist("y23",logy23,wgt)
!
! Jet clustering:
!
! Old way, which is fine for NLO but for NNLO we can increase
! speed a bit:
! Loop over different jet algos:
!  do ialgo=1,njetalgo
! Loop over different ycut values:
!    do iycut=15,1,-1
!      call eeJetAlgo(ialgo,Q2,ntrack,ptrack,ycut(iycut),njet)
! If the event does not pass the jet cut all the others will not
! pass either:
!      if (njet.eq.3) then
!        call fill_hist(cjetalgo(ialgo),dble(iycut),wgt)
!      elseif (njet.lt.3) then
!        goto 100
!      end if
!    end do
!    100 continue
!  end do
!
! New way:
! For each jet algorithm we calculate the transition variables:
  do ialgo=1,njetalgo
    if (ialgo.ne.5) then
      call CalcYn(ialgo,Q2,ntrack,ptrack,YnArr)
    else
      call CalcCambridgeYn(ialgo,Q2,ntrack,ptrack, &
                           ycut(15),YnArr)
    end if
! Loop over all ycut values:
    do iycut=15,1,-1
! To have exactly 3 jets: Y34 < ycut < Y23:
! That is the resolution should be less than Y34 not to 
! resolve 4 jets, but greater than Y23 to resolve at least 3
! jets:
      if ((ycut(iycut).gt.YnArr(2)).and.(ycut(iycut).lt.YnArr(1))) then
        call fill_hist(cjetalgo(ialgo),dble(iycut),wgt)
      end if
    end do
  end do
!
end subroutine analysis
!
! This routine implements all the relevant jet algoithms for
! ee-type colliders:
! the operation can be controlled via algo:
! algo = 1 : Durham, E-scheme
! algo = 2 : Geneva, E-scheme
! algo = 3 : Jade, E-scheme
! algo = 4 : Jade, E0-scheme
! algo = 5 : Cambridge, E-scheme
! Q2 : CM energy squared
subroutine eeJetAlgo(algo,Q2,ntrack,ptrack,ycut,njet)
implicit none
!
  integer , intent(in) :: algo
  real(kind(1d0)) , intent(in) :: Q2
  integer , intent(in) :: ntrack
  real(kind(1d0)) , dimension(4,ntrack) , intent(in) :: ptrack
  real(kind(1d0)) , intent(in) :: ycut
  integer , intent(out) :: njet
!
  integer :: itrack,i,j,itrk,jtrk
  integer(kind(1d0)) , dimension(ntrack) :: jetvec
  real(kind(1d0)) :: vijmin,vij,yij,piabs,pjabs,pijabs,costhetaij,Ei,Ej
  real(kind(1d0)) , dimension(4) :: p_tmp
  real(kind(1d0)) , dimension(4,ntrack) :: pjet
!
!
  njet = 0
!
! Copy the tracks into the jets array:
  pjet = ptrack
! Defining a jetvec to keep track of which pseudojets are combined away
! and which are still present:
  jetvec = 1
! When jetvec equals -1 it means the particle is recombined with an
! other one, if it is 0 it is a resolved jet (used in the Cambridge algo.)
!
! Note that vijmin is the same as yijmin except for the Cambridge-algo:
  vijmin = 1d99
! Loop only terminates if yij < ycut for pair corresponding to vij = vijmin
  do while (.true.)
    do i=1,ntrack-1
! If the actual jetvec item is -1 the jet is already combined with
! another one, or if it is 0 it is a resolved jet in the Cambridge algo:
      if (jetvec(i).ne.1) cycle
      Ei = pjet(4,i)
      piabs = sqrt(sum(pjet(1:3,i)**2))
      do j=i+1,ntrack
        if (jetvec(j).ne.1) cycle
        Ej = pjet(4,j)
        pjabs = sqrt(sum(pjet(1:3,j)**2))
! The cosine of the enclosed angle is always needed, hence
! calculate it for once and for all:
        if ((piabs.eq.0).or.(pjabs.eq.0)) then
          costhetaij = 1
        else
          costhetaij = (pjet(1,i)*pjet(1,j) + pjet(2,i)*pjet(2,j) &
                     + pjet(3,i)*pjet(3,j))/(piabs*pjabs)
        end if
! Calculation of vij:
! Durham-algo, vij = yij:
        if (algo.eq.1) then
          vij = 2*min(Ei**2,Ej**2)*(1 - costhetaij)/Q2
! Geneva-algo, vij = yij:
        elseif (algo.eq.2) then
          vij = 8d0/9d0*Ei*Ej*(1 - costhetaij)/(Ei + Ej)**2
! Jade-algo with E0 or with E, vij = yij:
        elseif ((algo.eq.3).or.(algo.eq.4)) then
          vij = 2*Ei*Ej*(1 - costhetaij)/Q2
! Cambridge-algo, vij != yij:
        elseif (algo.eq.5) then
          vij = 2*(1 - costhetaij)
        else
          print *,"Error is eeJetAlgo, algorithm is not implemented yet"
          stop "eeJetAlgo..."
        end if
        if (vij.lt.vijmin) then
          itrk = i
          jtrk = j
          vijmin = vij
!          print *,"itrk,jtrk: ",itrk,jtrk
!          print *,"vijmin changed: ",vijmin
        end if
      end do
    end do
! If the algorithm was the Cambridge one vij is defined
! differently from yij hence yij has to be calculated:
    if (algo.eq.5) then
      Ei = pjet(4,itrk)
      Ej = pjet(4,jtrk)
      piabs = sqrt(sum(pjet(1:3,itrk)**2))
      pjabs = sqrt(sum(pjet(1:3,jtrk)**2))
      if ((piabs.eq.0).or.(pjabs.eq.0)) then
        costhetaij = 1
      else
        costhetaij = (pjet(1,itrk)*pjet(1,jtrk) + pjet(2,itrk)*pjet(2,jtrk) &
                   + pjet(3,itrk)*pjet(3,jtrk))/(piabs*pjabs)
      end if
      yij = 2*min(Ei**2,Ej**2)*(1 - costhetaij)/Q2
! Otherwise yij is just vijmin:
    else
      yij = vijmin
    end if
! If the smallest distance is larger than the ycut the loop is quit,
! unless having the Cambridge algorithm:
    if ((yij.gt.ycut).and.(algo.ne.5)) then
      exit
! In the case of the Cambridge algorithm that pseudojet which
! has the smaller energy is deleted from the table of pseudojets and
! advocated to a resolved jet:
    elseif ((yij.gt.ycut).and.(algo.eq.5)) then
      if (Ei.lt.Ej) then
        jetvec(itrk) = 0
      else
        jetvec(jtrk) = 0
      end if
! We have to count how many pseudojets are remaining:
      njet = 0
      do i=1,ntrack
        if (jetvec(i).eq.1) njet = njet + 1
      end do
! If more than one are present go back to the beginning:
      if (njet.gt.1) then
        vijmin = 1d99
        cycle
! Otherwise quit the loop:
      else
        exit
      end if
    end if
! The recombination scheme varies from algo to algo:
! For Durham, Geneva,Jade-E and Cambridge we use the E-scheme:
    if (algo.ne.3) then
      p_tmp = pjet(:,itrk) + pjet(:,jtrk)
! For Jade-E0 which is algo = 3 we use the E0-scheme:
    else
      p_tmp = pjet(:,itrk) + pjet(:,jtrk)
! The spatial part has to be rescaled:
      pijabs = sqrt(sum(p_tmp(1:3)**2))
      p_tmp(1:3) = p_tmp(4)/pijabs*p_tmp(1:3)
    end if
    pjet(:,itrk) = p_tmp
! Nullifying the jet momentum:
    pjet(:,jtrk) = 0
! And putting the actual jetvec item to -1:
    jetvec(jtrk) = -1
    vijmin = 1d99
  end do
!
  njet = 0
  do i=1,ntrack
    if (jetvec(i).ne.-1) njet = njet + 1
  end do
!
end subroutine eeJetAlgo
!
! This routine Calculates the Durham y_{n-1 n}:
subroutine CalcDurhamyn(ptrack,ntrack,Q2,njet,ymin)
implicit none
!
  integer , intent(in) :: ntrack
  real(kind(1d0)) , dimension(4,ntrack) , intent(in) :: ptrack
  real(kind(1d0)) , intent(in) :: Q2
  integer , intent(in) :: njet
  real(kind(1d0)) , intent(out) :: ymin
!
  integer :: i,j,itrk,jtrk,mjet
  integer , dimension(ntrack) :: jetvec
  real(kind(1d0)) :: Ei,Ej,piabs,pjabs,costhetaij
  real(kind(1d0)) :: yij
  real(kind(1d0)) , dimension(4,ntrack) :: pjet
!
!
  ymin = -1d99
!
! Copying the tracks into the jet array:
  pjet = ptrack
! Filling jetvec as well:
  jetvec = 1
!
  do while (.true.)
! counting jets:
    mjet = 0
    do i=1,ntrack
! Skipping already combined jets:
      if (jetvec(i).eq.-1) cycle
      mjet = mjet + 1
    end do
! When the number of jets is less than njet we return with 
! ymin
    if (mjet.lt.njet) return
    ymin = 1d99
    do i=1,ntrack-1
      if (jetvec(i).eq.-1) cycle
      Ei = pjet(4,i)
      piabs = sqrt(sum(pjet(1:3,i)**2))
      do j=i+1,ntrack
        if (jetvec(j).eq.-1) cycle
        Ej = pjet(4,j)
        pjabs = sqrt(sum(pjet(1:3,j)**2))
        if ((piabs.eq.0).or.(pjabs.eq.0)) then
          costhetaij = 1
        else
          costhetaij = (pjet(1,i)*pjet(1,j) + pjet(2,i)*pjet(2,j) &
                     +  pjet(3,i)*pjet(3,j))/(piabs*pjabs)
        end if
        yij = 2*min(Ei**2,Ej**2)*(1 - costhetaij)/Q2
        if (yij.lt.ymin) then
          itrk = i
          jtrk = j
          ymin = yij
        end if
      end do
    end do
! Recombine using the E-scheme:
    pjet(:,itrk) = pjet(:,itrk) + pjet(:,jtrk)
    pjet(:,jtrk) = 0
    jetvec(jtrk) = -1
  end do
!
end subroutine CalcDurhamyn
!
! This is a general routine calculating all the n-to-n-1
! jet transition variables:
subroutine CalcYn(algo,Q2,ntrack,ptrack,Yn)
implicit none
!
  integer , intent(in) :: algo
  real(kind(1d0)) , intent(in) :: Q2
  integer , intent(in) :: ntrack
  real(kind(1d0)) , dimension(4,ntrack) , intent(in) :: ptrack
  real(kind(1d0)) , dimension(ntrack) , intent(out) :: Yn
!
  integer :: itrack,i,j,itrk,jtrk,mjet
  integer(kind(1d0)) , dimension(ntrack) :: jetvec
  real(kind(1d0)) :: vijmin,vij,yij,piabs,pjabs,pijabs,costhetaij,Ei,Ej
  real(kind(1d0)) , dimension(4) :: p_tmp
  real(kind(1d0)) , dimension(4,ntrack) :: pjet
!
!
  Yn = 0
!
! Copy the tracks into the jets array:
  pjet = ptrack
! Defining a jetvec to keep track of which pseudojets are combined away
! and which are still present:
  jetvec = 1
! When jetvec equals -1 it means the particle is recombined with an
! other one, if it is 0 it is a resolved jet (used in the 
! Cambridge algo., though not used here...)
!
! In the zeroth iteration we have as many jets as tracks,
! hence if we have n tracks the distance to resolve n+1 jets
! is exactly zero:
  yij = 0
!
  do while (.true.)
! Counting the jets:
    mjet = 0
    do i=1,ntrack
! Skipping recombined pseudojets:
      if (jetvec(i).eq.-1) cycle
      mjet = mjet + 1
    end do
! If the number of remaining jets is less than 2 we quit:
    if (mjet.lt.2) return
! The y value where the transition happens from a mjet + 1
! jet configuration to a mjet configuration is stored in yij:
! Note, that the 3-to-2 jet transition is stored in the first
! position:
    Yn(mjet-1) = yij
! vijmin is the minimal distance found in the current 
! iteration which corresponds to the yij except for the 
! Cambridge algorithm:
    vijmin = 1d99
    do i=1,ntrack-1
      if (jetvec(i).eq.-1) cycle
      Ei = pjet(4,i)
      piabs = sqrt(sum(pjet(1:3,i)**2))
      do j=i+1,ntrack
        if (jetvec(j).eq.-1) cycle
        Ej = pjet(4,j)
        pjabs = sqrt(sum(pjet(1:3,j)**2))
        if ((piabs.eq.0).or.(pjabs.eq.0)) then
          costhetaij = 1
        else
          costhetaij = (pjet(1,i)*pjet(1,j) + pjet(2,i)*pjet(2,j) &
                     +  pjet(3,i)*pjet(3,j))/(piabs*pjabs)
        end if
! Calculation of vij:
! Durham-algo, vij = yij:
        if (algo.eq.1) then
          vij = 2*min(Ei**2,Ej**2)*(1 - costhetaij)/Q2
! Geneva-algo, vij = yij:
        elseif (algo.eq.2) then
          vij = 8d0/9d0*Ei*Ej*(1 - costhetaij)/(Ei + Ej)**2
! Jade-algo with E0 or with E, vij = yij:
        elseif ((algo.eq.3).or.(algo.eq.4)) then
          vij = 2*Ei*Ej*(1 - costhetaij)/Q2
! Cambridge-algo, vij != yij:
        elseif (algo.eq.5) then
          vij = 2*(1 - costhetaij)
          print *,"Error!!!! This routine cannot be used with this algo"
          stop "CalcYn"
        else
          print *,"Error is CalcYn, algorithm is not implemented yet"
          stop "CalcYn..."
        end if
        if (vij.lt.vijmin) then
          itrk = i
          jtrk = j
          vijmin = vij
!          print *,"itrk,jtrk: ",itrk,jtrk
!          print *,"vijmin changed: ",vijmin
        end if
      end do
    end do
! If the algorithm was the Cambridge one vij is defined
! differently from yij hence yij has to be calculated:
    if (algo.eq.5) then
      Ei = pjet(4,itrk)
      Ej = pjet(4,jtrk)
      piabs = sqrt(sum(pjet(1:3,itrk)**2))
      pjabs = sqrt(sum(pjet(1:3,jtrk)**2))
      if ((piabs.eq.0).or.(pjabs.eq.0)) then
        costhetaij = 1
      else
        costhetaij = (pjet(1,itrk)*pjet(1,jtrk) + pjet(2,itrk)*pjet(2,jtrk) &
                   + pjet(3,itrk)*pjet(3,jtrk))/(piabs*pjabs)
      end if
      yij = 2*min(Ei**2,Ej**2)*(1 - costhetaij)/Q2
! Otherwise yij is just vijmin:
    else
      yij = vijmin
    end if
! We calculate transition functions, hence the pair with the
! smallest separation is always recombined:
! The recombination scheme varies from algo to algo:
! For Durham, Geneva,Jade-E and Cambridge we use the E-scheme:
    if (algo.ne.3) then
      p_tmp = pjet(:,itrk) + pjet(:,jtrk)
! For Jade-E0 which is algo = 3 we use the E0-scheme:
    else
      p_tmp = pjet(:,itrk) + pjet(:,jtrk)
! The spatial part has to be rescaled:
      pijabs = sqrt(sum(p_tmp(1:3)**2))
      p_tmp(1:3) = p_tmp(4)/pijabs*p_tmp(1:3)
    end if
    pjet(:,itrk) = p_tmp
! Nullifying the jet momentum:
    pjet(:,jtrk) = 0
! And putting the actual jetvec item to -1:
    jetvec(jtrk) = -1
  end do
!
end subroutine CalcYn
!
subroutine CalcCambridgeYn(algo,Q2,ntrack,ptrack,ymin,Yn)
implicit none
!
  integer , intent(in) :: algo
  real(kind(1d0)) , intent(in) :: Q2
  integer , intent(in) :: ntrack
  real(kind(1d0)) , dimension(4,ntrack) , intent(in) :: ptrack
  real(kind(1d0)) , intent(in) :: ymin
  real(kind(1d0)) , dimension(ntrack) , intent(out) :: Yn
!
  integer :: i,j,itrk,jtrk,mjet
  integer , dimension(ntrack) :: jetvec
  real(kind(1d0)) :: yijmax,vijmin,ycut,yij,vij
  real(kind(1d0)) :: Ei,Ej,piabs,pjabs,costhetaij
  real(kind(1d0)) , dimension(4) :: p_tmp
  real(kind(1d0)) , dimension(4,ntrack) :: pjet
!
!
  Yn = 0
!
! Setting yijmax to the maximal value:
  yijmax = 1
!
  do while (.true.)
! Coping the content of ptrack into pjet:
    pjet = ptrack
! Defining a jetvec array, to hold info about recombined jets:
    jetvec = 1
! ycut is made equal to the yijmax found in the previous iteration:
    ycut = yijmax
    yijmax = 0
! If ycut reached the minimum we quit:
    if (ycut.lt.ymin) then
      exit
    end if
    do while (.true.)
      vijmin = 1d99
! Calculating jets:
      mjet = 0
      do i=1,ntrack
        if (jetvec(i).ne.-1) mjet = mjet + 1
      end do
! If only one jet is remaining we can quit the recombination loop:
      if (mjet.eq.1) exit
      do i=1,ntrack-1
! We consider the Cambridge algorithm, hence 0 is also allowed:
        if (jetvec(i).ne.1) cycle
        Ei = pjet(4,i)
        piabs = sqrt(sum(pjet(1:3,i)**2))
        do j=i+1,ntrack
          if (jetvec(j).ne.1) cycle
          Ej = pjet(4,j)
          pjabs = sqrt(sum(pjet(1:3,j)**2))
          if ((piabs.eq.0).or.(pjabs.eq.0)) then
            costhetaij = 1
          else
            costhetaij = (pjet(1,i)*pjet(1,j) + pjet(2,i)*pjet(2,j) &
                       +  pjet(3,i)*pjet(3,j))/(piabs*pjabs)
          end if
          vij = 2*(1 - costhetaij)
          yij = 2*min(Ei**2,Ej**2)*(1 - costhetaij)/Q2
! We have to keep the maximal value for yij which should be
! less than ycut since pairs with yij > ycut are already treated:
          if ((yij.gt.yijmax).and.(yij.lt.ycut)) yijmax = yij
          if (vij.lt.vijmin) then
            itrk = i
            jtrk = j
            vijmin = vij
!            print *,"itrk,jtrk: ",itrk,jtrk
!            print *,"vijmin changed: ",vijmin
          end if
        end do
      end do
! At this point the closest pair is identified, by 
! calculating yij and comparing it with ycut it can be
! determined wether they will be recombined or not:
! Calculating yij for the closest pair:
      Ei = pjet(4,itrk)
      Ej = pjet(4,jtrk)
      piabs = sqrt(sum(pjet(1:3,itrk)**2))
      pjabs = sqrt(sum(pjet(1:3,jtrk)**2))
      if ((piabs.eq.0).or.(pjabs.eq.0)) then
        costhetaij = 1
      else
        costhetaij = (pjet(1,itrk)*pjet(1,jtrk) + pjet(2,itrk)*pjet(2,jtrk) &
                   + pjet(3,itrk)*pjet(3,jtrk))/(piabs*pjabs)
      end if
      yij = 2*min(Ei**2,Ej**2)*(1 - costhetaij)/Q2
! Having a resolved jet:
      if (yij.ge.ycut) then
        if (Ei.lt.Ej) then
          jetvec(itrk) = 0
        else
          jetvec(jtrk) = 0
        end if
! Counting remaining pseudojets to decide what to do:
        mjet = 0
        do i=1,ntrack
          if (jetvec(i).eq.1) mjet = mjet + 1
        end do
! At least two pseudojets are present, continue...
        if (mjet.gt.1) then
          cycle
! Only one pseudojet remains quit the loop...
        else
          exit
        end if
      end if
      p_tmp = pjet(:,itrk) + pjet(:,jtrk)
      pjet(:,itrk) = p_tmp
! Nullifying the jet momentum:
      pjet(:,jtrk) = 0
! And putting the actual jetvec item to -1:
      jetvec(jtrk) = -1
    end do
! Counting resolved jets to determine which Yn is
! to be assigned the value of ycut:
    mjet = 0
    do i=1,ntrack
      if (jetvec(i).ne.-1) mjet = mjet + 1
    end do
    if (mjet.gt.1) Yn(mjet-1) = yijmax
  end do
!
end subroutine CalcCambridgeYn
