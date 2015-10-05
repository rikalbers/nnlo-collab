! This program is written to test the reshuffle algorithm:
module reshuffle_data
implicit none
!
  integer , parameter :: maxpart = 20
  integer , parameter :: maxpartons = 9
!
contains
!
  subroutine PrintPatterns(n,numptrn,ptrn)
  implicit none
!
    integer , intent(in) :: n,numptrn
    integer , intent(in) , dimension(:,:) :: ptrn
!
    integer :: iptrn,ipart
!
!
! We print out the patterns:
    print 100,"The number of patterns we have: ",numptrn
    100 format(A,I2)
    do iptrn=1,numptrn
      print 101,"The ",iptrn,"th pattern is: ",ptrn(:,iptrn)
      101 format(A,I2,A,20(I3,1X))
    end do
!
  end subroutine PrintPatterns
!
end module reshuffle_data
!
subroutine reshufflemom(n,p,pout,numptrn,ptrn,returnptrn,prefactor)
use momenta
use particles
use process
use QCDparams
use reshuffle_data
implicit none
!
  integer , intent(in) :: n,numptrn
  type(particle) , intent(in) , dimension(:) :: p
  type(particle) , intent(out) , dimension(:) :: pout
  integer , intent(in) , dimension(:,:) :: ptrn
  integer , intent(out) :: returnptrn
  real(kind(1d0)) , intent(out) :: prefactor
!
  integer , dimension(maxpartons) :: flavin
  integer , dimension(maxpartons) :: flav_ptrn
  integer , dimension(maxpartons) :: io,ordering
  integer :: ipart,iparton,jparton,iptrn
  integer :: nparton
  integer :: nextparton,partonid
  integer , dimension(maxpartons) :: partonlist
  logical :: found
  real(kind(1d0)) :: prefact
  type(mom) , dimension(maxpart) :: pmom_tmp
!
!
  nextparton = 1
  iparton = 0
  partonlist = 0
  returnptrn = -1000
  prefactor = 0d0
  flavin = 0
!  print *,p(:)%flv
! We encode the input flavour configuration with our scheme:
  do ipart=1,n
!    print *,"flav(ipart) is: ",p(ipart)%flv
! skip heavy and/or non-QCD stuff:
    if ((ipart.gt.2).and.(ipart.lt.proc_firstlight)) then
      flavin(ipart) = p(ipart)%flv
      cycle
    end if
    if (abs(p(ipart)%flv).gt.qcd_nf) then
      flavin(ipart) = p(ipart)%flv
      cycle
    end if
!    print *,"iparton is increased..."
    iparton = iparton + 1
! we select light, non-gluon partons:
    if ((p(ipart)%flv.ne.0).and.(abs(p(ipart)%flv).le.qcd_nf)) then
!      print *,"quark!"
! We presume that the partontype is not yet stored in partonlist:
      partonid = nextparton
! We search for it in the partonlist:
      do jparton=1,nextparton-1
! The presumption was wrong we found it, so use jparton for id:
        if (partonlist(jparton).eq.abs(p(ipart)%flv)) then
!          print *,"We found the quark in the list!"
          partonid = jparton
          exit
        end if
      end do
! If we did not find the parton in the list we have to insert it:
      if (partonid.eq.nextparton) then
        nextparton = nextparton + 1
        partonlist(partonid) = abs(p(ipart)%flv)
      end if
! In any case we insert it into the flavour list:
      flavin(ipart) = sign(partonid,p(ipart)%flv)
! If it is a gluon we simply include a 'g':
    elseif (p(ipart)%flv.eq.0) then
!      print *,"gluon!"
      flavin(ipart) = 0
    end if
  end do
  nparton = iparton
!  print 200,flavin(1:n)
!  print *,flavin(1:n)
!
! At this point we bear the flavor configuration in an internal 
! representation.
! We have to take all the available patterns and compare flavin with 
! them:
  do iptrn=1,numptrn
    flav_ptrn(1:n) = ptrn(1:n,iptrn)
    iparton = 0
    call CompFlavors(n,flav_ptrn,flavin,io,ordering,found)
!    print *,"iptrn, found: ",iptrn,found
    if (found) then
!      print 200,ordering(1:n)
!      print 200,io(1:n)
      exit
    end if
  end do
!
! For safety measures we check wether we indeed found the pattern:
  if (.not.found) then
    print *,"We failed to found the pattern for: ",flavin(1:n)
    stop
  end if
!
! Since we can have non-QCD and/or massive partons we have to
! modify the numbers in ordering if they are larger than 2:
  where (ordering > 2)
    ordering = ordering + proc_firstlight - 3
  end where
!  print *,"The ordering when nonQCD part is considered: "
!  print 200,ordering(1:n)
!
! We have to take the original pmom array and reshuffle the momenta
! using a previously specified convention:
  do ipart=1,n
    pmom_tmp(ipart) = p(ordering(ipart))%p
  end do
!
! We give back the momenta and the flavor in a different array:
  call CreateParts(pmom_tmp(1:n),flav_ptrn(1:n),pout(1:n))
!
! We calculate the needed new prefactor as well:
  call CalcPrefactor(n,flavin,prefact)
  prefactor = prefact
!
!
  returnptrn = iptrn
!
  200 format(20(I3,1x))
!
end subroutine reshufflemom
!
! This routine is the same as the previous one but it distinguishes
! between up and down quarks...
subroutine reshufflemomud(n,p,pout,numptrn,ptrn,returnptrn,prefactor)
use momenta
use particles
use process
use QCDparams
use reshuffle_data
implicit none
!
  integer , intent(in) :: n,numptrn
  type(particle) , intent(in) , dimension(:) :: p
  type(particle) , intent(out) , dimension(:) :: pout
  integer , intent(in) , dimension(:,:) :: ptrn
  integer , intent(out) :: returnptrn
  real(kind(1d0)) , intent(out) :: prefactor
!
  integer , dimension(maxpartons) :: flavin
  integer , dimension(maxpartons) :: flav_ptrn
  integer , dimension(maxpartons) :: io,ordering
  integer :: ipart,ipartonup,ipartondown,jparton,iptrn
  integer :: nparton
  integer :: nextparton,partonid
  integer :: nextup,nextdown
  integer , dimension(maxpartons) :: partonlist
  integer , dimension(maxpartons) :: upquarklist,downquarklist
  logical :: found
  real(kind(1d0)) :: prefact
  type(mom) , dimension(maxpart) :: pmom_tmp
!
!
  nextdown = 1
  nextup = 1
  partonlist = 0
  upquarklist = 0
  downquarklist = 0
  returnptrn = -1000
  prefactor = 0d0
  flavin = 0
!  print *,p(:)%flv
! We encode the input flavour configuration with our scheme:
  do ipart=1,n
!    print *,"flav(ipart) is: ",p(ipart)%flv
! skip heavy and/or non-QCD stuff:
    if ((ipart.gt.2).and.(ipart.lt.proc_firstlight)) then
      flavin(ipart) = p(ipart)%flv
      cycle
    end if
    if (abs(p(ipart)%flv).gt.qcd_nf) then
      flavin(ipart) = p(ipart)%flv
      cycle
    end if
! we select light, down quarks:
    if ((p(ipart)%flv.ne.0).and.(abs(p(ipart)%flv).le.qcd_nf).and. &
        (abs(mod(p(ipart)%flv,2)).eq.1)) then
!      print *,"down quark!"
! We presume that the partontype is not yet stored in downquarklist:
      partonid = nextdown
! We search for it in the partonlist:
      do jparton=1,nextdown-1
! The presumption was wrong we found it, so use jparton for id:
        if (downquarklist(jparton).eq.abs(p(ipart)%flv)) then
!          print *,"We found the down quark in the list!"
          partonid = jparton
          exit
        end if
      end do
! If we did not find the parton in the list we have to insert it:
      if (partonid.eq.nextdown) then
        nextdown = nextdown + 1
        downquarklist(partonid) = abs(p(ipart)%flv)
      end if
! In any case we insert it into the flavour list:
      flavin(ipart) = sign(2*partonid-1,p(ipart)%flv)
! If it is a gluon we simply include a 'g':
    elseif ((p(ipart)%flv.ne.0).and.(abs(p(ipart)%flv).le.qcd_nf).and. &
        (abs(mod(p(ipart)%flv,2)).eq.0)) then
!      print *,"up quark!"
! We presume that the partontype is not yet stored in upquarklist:
      partonid = nextup
! We search for it in the partonlist:
      do jparton=1,nextup-1
! The presumption was wrong we found it, so use jparton for id:
        if (upquarklist(jparton).eq.abs(p(ipart)%flv)) then
!          print *,"We found the up quark in the list!"
          partonid = jparton
          exit
        end if
      end do
! If we did not find the parton in the list we have to insert it:
      if (partonid.eq.nextup) then
        nextup = nextup + 1
        upquarklist(partonid) = abs(p(ipart)%flv)
      end if
! In any case we insert it into the flavour list:
      flavin(ipart) = sign(2*partonid,p(ipart)%flv)
! If it is a gluon we simply include a 'g':
    elseif (p(ipart)%flv.eq.0) then
!      print *,"gluon!"
      flavin(ipart) = 0
    else
!      print *,"Problem occured in reshufflemomud..."
!      print *,"Flavor: ",flavin(ipart)
      stop
    end if
  end do
!  print 200,flavin(1:n)
!  print *,flavin(1:n)
!
! At this point we bear the flavor configuration in an internal 
! representation.
! We have to take all the available patterns and compare flavin with 
! them:
  do iptrn=1,numptrn
    flav_ptrn(1:n) = ptrn(1:n,iptrn)
    call CompFlavors(n,flav_ptrn,flavin,io,ordering,found)
!    print *,"iptrn, found: ",iptrn,found
    if (found) then
!      print 200,ordering(1:n)
!      print 200,io(1:n)
      exit
    end if
  end do
!
! For safety measures we check wether we indeed found the pattern:
  if (.not.found) then
    print *,"We failed to find the pattern for: ",flavin(1:n)
    stop
  end if
!
! Since we can have non-QCD and/or massive partons we have to
! modify the numbers in ordering if they are larger than 2:
! FIXME at this point I am not sure whether the next 3 lines are
! needed or not...
!  where (ordering > 2)
!    ordering = ordering + proc_firstlight - 3
!  end where
!  print *,"The ordering when nonQCD part is considered: "
!  print 200,ordering(1:n)
!
! We have to take the original pmom array and reshuffle the momenta
! using a previously specified convention:
  do ipart=1,n
    pmom_tmp(ipart) = p(ordering(ipart))%p
  end do
!
! We give back the momenta and the flavor in a different array:
  call CreateParts(pmom_tmp(1:n),flav_ptrn(1:n),pout(1:n))
!
! We calculate the needed new prefactor as well:
  call CalcPrefactor(n,flavin,prefact)
  prefactor = prefact
!
!
  returnptrn = iptrn
!
  200 format(20(I3,1x))
!
end subroutine reshufflemomud
!
subroutine CompFlavors(npart,flavA,flavB,io,ordering,found)
use reshuffle_data
implicit none
!
  integer :: npart
  integer , dimension(maxpartons) , intent(in) :: flavA,flavB
  integer , dimension(maxpartons) , intent(out) :: io,ordering
  logical , intent(out) :: found
!
  integer :: ipart,jpart
  integer :: flv1,flv2
  integer , dimension(maxpartons) :: flav_tmp
!
!
  found = .false.
  io(1:2) = 1
  io(3:npart) = -1
  ordering = 0
!
!  print *,"The input flavor configs are: "
!  print 110,"flavA : ",flavA(1:npart)
!  print 110,"flavB : ",flavB(1:npart)
!
! We keep flavA and we try to shuffle the flavors in flavB to arrive at
! flavA:
  flav_tmp = flavB
!
! We go through all partons in flavA and try to find the same in flav_tmp:
  do ipart=1,npart
    do jpart=1,npart+1
      if (jpart.eq.npart+1) return
      if (flav_tmp(jpart).eq.99) cycle
      if (((ipart.lt.3).and.(jpart.gt.2)).or. &
          ((ipart.gt.2).and.(jpart.lt.3))) then
        flv1 = -flavA(ipart)
      else
        flv1 =  flavA(ipart)
      end if
      flv2 = flav_tmp(jpart)
      if (flv1.eq.flv2) then
        ordering(ipart) = jpart
        flav_tmp(jpart) = 99
        exit
      end if
    end do
  end do
!
! When we found a good pattern we have to save the io array as well,
! which is used in HELAC but can be useful even in NJETS too:
  do ipart=1,npart
    if (((ipart.lt.3).and.(ordering(ipart).ge.3)).or. &
        ((ipart.ge.3).and.(ordering(ipart).lt.3))) then
      io(ipart) = -io(ipart)
    end if
  end do
!
  found = .true.
!
  100 format(20(I3,1x))
  110 format(A,20(I3,1x))
!
end subroutine CompFlavors
!
! This routine presents a simple way how to recalculate the symmetry
! factors associated to the colour factor of initial states and the
! factor associated to identical partons in the final state. The 
! factor is calculated for each evaluation, so the efficiency can
! be at least argued, but this routine is intended to be used for
! complex processes, that is, the CPU time needed to perform these
! calculations is negligible to calculations related to the dynamics.
subroutine CalcPrefactor(npart,flav,prefact)
use QCDparams
use reshuffle_data
implicit none
!
  integer , intent(in) :: npart
  integer , intent(in) , dimension(maxpartons) :: flav
  real(kind(1d0)) , intent(out) :: prefact
!
  integer :: ipart
  integer , dimension(-maxpartons:maxpartons) :: numparton = 0
!
  integer :: ifactorial
  external :: ifactorial
!
  prefact = 1d0
  numparton = 0
! We loop through all final state particles and count the
! number of each flavour.
  do ipart=3,npart
    if (abs(flav(ipart)).gt.qcd_nf) cycle
    numparton(flav(ipart)) = numparton(flav(ipart)) + 1
  end do
!
!  print *,"Number of different partons: "
!  print *,numparton
!
  do ipart=-maxpartons,maxpartons
    if (numparton(ipart).ne.0) prefact = prefact*ifactorial(numparton(ipart))
  end do
!
! We also include the colour factor associated to the initial state:
  do ipart=1,2
    if (abs(flav(ipart)).gt.qcd_nf) cycle
    if (flav(ipart).ne.0) prefact = prefact*3d0
    if (flav(ipart).eq.0) prefact = prefact*8d0
  end do
!
!  print *,"The prefactor is: ",prefact
!
end subroutine CalcPrefactor
!
recursive integer function ifactorial(i) result(ans)
implicit none
!
  integer , intent(in) :: i
!
!
!
  if (i==0) then
    ans = 1
  else
    ans = i*ifactorial(i-1)
  end if
!
end function ifactorial
