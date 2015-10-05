! This source contains routines for checking the various ensembles
! of subtraction terms:
module RegCheck
implicit none
!
  logical , parameter :: checksingly = .true.
  logical , parameter :: checkdoubly = .true.
  integer , parameter :: niter = 15
  integer , parameter :: verbosity = 3
  real(kind(1d0)) :: xi1,xi2
  real(kind(1d0)) :: y1,y2
  real(kind(1d0)) :: azi1,azi2
!
  real(kind(1d0)) , dimension(:) , allocatable :: sme
  real(kind(1d0)) , dimension(:) , allocatable :: subs1
!
  interface FindDoublyUnresolved
    module procedure FindDoublyUnresolved_2parton
    module procedure FindDoublyUnresolved_3parton
    module procedure FindDoublyUnresolved_4parton
  end interface FindDoublyUnresolved
  interface GenerateDoublyUnresLimits
    module procedure GenerateDoublyUnresLimits_2parton
    module procedure GenerateDoublyUnresLimits_3parton
    module procedure GenerateDoublyUnresLimits_4parton
  end interface GenerateDoublyUnresLimits
!
contains
!
subroutine CheckValidity(iter,sme,subs,centscale,passed)
implicit none
!
  integer , intent(in) :: iter
  real(kind(1d0)) , dimension(:) , intent(in) :: sme
  real(kind(1d0)) , dimension(:) , intent(in) :: subs
  integer , intent(in) :: centscale
  logical , intent(out) :: passed
!
  character(len=15) :: errmsg
  integer :: iscale,nscale,scale1,scale2
  real(kind(1d0)) :: ratio,dist
!
!
  passed = .true.
!
! If verbosity is odd only checking the central scale:
  if (mod(verbosity,2).eq.1) then
    scale1 = centscale
    scale2 = scale1
! Otherwise we check for each scale choice:
  else
    scale1 = 1
    scale2 = size(sme)
  end if
!
  do iscale=scale1,scale2
    if ((sme(iscale).eq.0).or.(subs(iscale).eq.0)) then
      write(*,*) "Error in CheckValidity..."
      write(*,"(A,1x,I0)") "SME or Subtraction is zero for scale: ",iscale
      write(*,"(A,1x,I0)") "In iteration ",iter
      write(*,"(A,1x,G0)") "SME: ",sme(iscale)
      write(*,"(A,1x,G0)") "Subs: ",subs(iscale)
      write(*,*) "Better check the cuts..."
!      stop "CheckValidity"
!      read(*,*)
    end if
! To put out some error message we use the convention of 
! POWHEG-BOX:
    ratio = subs(iscale)/sme(iscale)
    dist = abs(ratio - 1)
    errmsg = " "
    if (dist.gt.0.01) then
      if (iter.eq.1) then
        if (dist.lt.0.1) then
          errmsg = "*-WARN-*"
        else
          errmsg = "*-WWARN-*"
        end if
      elseif (iter.eq.3) then
        if (dist.lt.0.1) then
          errmsg = "*-WWARN-*"
        else
          errmsg = "*-WWWARN-*"
        end if
      elseif (iter.ge.4) then
        if (dist.lt.0.1) then
          errmsg = "*-WWWARN-*"
        elseif (dist.lt.0.3) then
          errmsg = "*-WWWWARN-*"
        else
          errmsg = "*-WWWWWARN-*"
        end if
      end if
      passed = .false.
    end if
    if (((errmsg.ne." ").and.(verbosity.ne.0)).or. &
        (verbosity.gt.2)) then
      write(*,"(2(A,1X,I2),1X,G0,1X,A)") "iter no.",iter," scale no.",iscale,ratio,errmsg
    end if
  end do
!
end subroutine CheckValidity
!
subroutine FindSinglyUnresolved(flv_born_arr,flv_real,radi,radr,foundregion)
use regions
use utils
implicit none
!
  integer , dimension(:,:) , intent(in) :: flv_born_arr
  integer , dimension(:) , intent(in) :: flv_real
  integer , intent(inout) :: radi,radr
  logical , intent(out) :: foundregion
!
  integer :: ipart,rpart,npart
  integer , dimension(size(flv_real)-1) :: flv_tmp
!
  npart = size(flv_real)
!
! radi and radr are given by previous valid region
  do ipart=radi,npart-1
! if the ipart-th particle is not a parton simply skip:
    if (abs(flv_real(ipart)).gt.6) then
      radr = ipart + 1
      cycle
    end if
    do rpart=radr+1,npart
! Both particles cannot be in the initial state:
      if (max(ipart,rpart).lt.3) cycle
! if the rpart-th particle is not a parton simply skip:
      if (abs(flv_real(rpart)).gt.6) cycle
!      print 100,"ipart,rpart: ",ipart,rpart
!      100 format(A,1x,I0,1x,I0)
! To have a splitting it has to be valid QCD-wise and the
! corresponding underlying Born should exist.
! Distinguish between initial and final state splitting:
! ISR:
      if (ipart.lt.3) then
!        print *,"ISR: "
        stop "FindSinglyUnresolved not implemented yet..."
! FSR:
      else
!        print *,"FSR: "
! Is it valid QCD-wise?
        if (.not.((flv_real(ipart)+flv_real(rpart).eq.0).or. &
                  (flv_real(ipart).eq.0).or. &
                  (flv_real(rpart).eq.0))) cycle
! Is there a corresponding underlying Born?
        call UndoBranching(ipart,rpart,flv_real,flv_tmp)
        if (.not.CheckEquivProc(npart-1,size(flv_born_arr,2), &
                                flv_born_arr,flv_tmp)) cycle
!        write(*,"(A)",advance='no') "Underlying Born: "
!        call PrintSubproc(flv_tmp)
! Good branching is found return with legs:
        radi = ipart
        radr = rpart
        foundregion = .true.
        return
      end if
    end do
    radr = ipart + 1
  end do
!
  foundregion = .false.
!
end subroutine FindSinglyUnresolved
!
subroutine FindDoublyUnresolved_2parton(flv_born_arr,flv_rreal, &
                                        radr,rads,foundregion)
use regions
use utils
implicit none
!
  integer , dimension(:,:) , intent(in) :: flv_born_arr
  integer , dimension(:) , intent(in) :: flv_rreal
  integer , intent(inout) :: radr,rads
  logical , intent(out) :: foundregion
!
  integer :: ipart,ileg,rpart,spart,npart
  integer , dimension(size(flv_rreal)-2) :: flv_tmp
!
  npart = size(flv_rreal)
!
! radr and rads are given by previous valid region
  do rpart=radr,npart-1
! if the rpart-th particle is not a parton simply skip:
    if (abs(flv_rreal(rpart)).gt.6) then
      rads = rpart + 1
      cycle
    end if
    do spart=rads+1,npart
! Both particles cannot be in the initial state:
      if (max(rpart,spart).lt.3) cycle
! if the spart-th particle is not a parton simply skip:
      if (abs(flv_rreal(spart)).gt.6) cycle
!      print 100,"rpart,spart: ",rpart,spart
      100 format(A,2(1x,I0))
! This routine is only used to identify two soft legs, hence
! only checking for quark pairs and gluon pairs:
      if (flv_rreal(rpart)+flv_rreal(spart).ne.0) cycle
! The construction of the underlying Born is by simply dropping
! legs rpart and spart:
      ileg = 0
      do ipart=1,npart
        if ((ipart.eq.rpart).or.(ipart.eq.spart)) cycle
        ileg = ileg + 1
        flv_tmp(ileg) = flv_rreal(ipart)
      end do
      if (.not.CheckEquivProc(npart-2,size(flv_born_arr,2), &
          flv_born_arr,flv_tmp)) cycle
!      write(*,"(A)",advance='no') "Underlying Born: "
!      call PrintSubproc(flv_tmp)
! Good branching is found return with legs:
      radr = rpart
      rads = spart
      foundregion = .true.
      return
    end do
    rads = rpart + 1
  end do
!
  foundregion = .false.
!
end subroutine FindDoublyUnresolved_2parton
!
subroutine FindDoublyUnresolved_3parton(flv_born_arr,flv_rreal, &
                                        radi,radr,rads,foundregion)
use regions
use utils
implicit none
!
  integer , dimension(:,:) , intent(in) :: flv_born_arr
  integer , dimension(:) , intent(in) :: flv_rreal
  integer , intent(inout) :: radi,radr,rads
  logical , intent(out) :: foundregion
!
  integer :: ipart,rpart,spart,npart
  integer , dimension(size(flv_rreal)-2) :: flv_tmp
!
  npart = size(flv_rreal)
!
! radi, radr and rads are given by previous valid region
  do ipart=radi,npart-2
! if the ipart-th particle is not a parton simply skip:
    if (abs(flv_rreal(ipart)).gt.6) then
      radr = ipart + 2
      rads = ipart + 2
      cycle
    end if
    do rpart=radr,npart-1
! Both particles cannot be in the initial state:
      if (max(ipart,rpart).lt.3) then
        rads = rpart + 1
        cycle
      end if
! if the rpart-th particle is not a parton simply skip:
      if (abs(flv_rreal(rpart)).gt.6) then
        rads = rpart + 1
        cycle
      end if
      do spart=rads+1,npart
! This parton cannot be in the initial state.
! If not a parton simply cycle:
        if (abs(flv_rreal(spart)).gt.6) cycle
!        print 100,"ipart,rpart,spart: ",ipart,rpart,spart
!        100 format(A,3(1x,I0))
! To have a splitting it has to be valid QCD-wise and the
! corresponding underlying Born should exist.
! Distinguish between initial and final state splitting:
! ISR:
        if (ipart.lt.3) then
!          print *,"ISR: "
          stop "FindDoublyUnresolved not implemented yet..."
! FSR:
        else
!          print *,"FSR: "
! Is it valid QCD-wise?
          if (.not.((flv_rreal(ipart)+flv_rreal(rpart)+flv_rreal(spart).eq.0).or. &
                    (flv_rreal(ipart)+flv_rreal(rpart)+flv_rreal(spart).eq.flv_rreal(ipart)).or. &
                    (flv_rreal(ipart)+flv_rreal(rpart)+flv_rreal(spart).eq.flv_rreal(rpart)).or. &
                    (flv_rreal(ipart)+flv_rreal(rpart)+flv_rreal(spart).eq.flv_rreal(spart)))) cycle
! Is there a corresponding underlying Born?
          call UndoBranchingMod(flv_rreal,flv_tmp,ipart,rpart,spart)
          if (.not.CheckEquivProc(npart-2,size(flv_born_arr,2), &
                                  flv_born_arr,flv_tmp)) cycle
!          write(*,"(A)",advance='no') "Underlying Born: "
!          call PrintSubproc(flv_tmp)
! Good branching is found return with legs:
          radi = ipart
          radr = rpart
          rads = spart
          foundregion = .true.
          return
        end if
      end do
      rads = rpart + 1
    end do
    radr = ipart + 1
  end do
!
  foundregion = .false.
!
end subroutine FindDoublyUnresolved_3parton
!
subroutine FindDoublyUnresolved_4parton(flv_born_arr,flv_rreal, &
                                        radi,radr,radj,rads,foundregion)
use regions
use utils
implicit none
!
  integer , dimension(:,:) , intent(in) :: flv_born_arr
  integer , dimension(:) , intent(in) :: flv_rreal
  integer , intent(inout) :: radi,radr,radj,rads
  logical , intent(out) :: foundregion
!
  integer :: ipart,rpart,jpart,spart,npart
  integer , dimension(size(flv_rreal)-2) :: flv_tmp
!
  npart = size(flv_rreal)
!
! Two separate two-particle splittings have to be identified:
  do ipart=radi,npart-2
    if (abs(flv_rreal(ipart)).gt.6) then
      radr = ipart + 2
      cycle
    end if
    do rpart=radr,npart
      if (abs(flv_rreal(rpart)).gt.6) then
        cycle
      end if
      do jpart=radj,npart-1
        if ((abs(flv_rreal(jpart)).gt.6)) then
          rads = jpart + 1
          cycle
        end if
        do spart=rads+1,npart
          if ((abs(flv_rreal(spart)).gt.6)) then
            cycle
          end if
          if (ipart.eq.jpart) cycle
          if (rpart.eq.spart) cycle
          if (rpart.eq.jpart) cycle
!          write(*,"(A,4(I0,1x))") "i,r,j,s: ",ipart,rpart,jpart,spart
! ISR:
          if (ipart.lt.3) then
!            print *,"ISR: "
            stop "FindDoublyUnresolved not implemented yet..."
! FSR:
          else
!            print *,"FSR: "
! Is it valid QCD-wise?
            if (.not.((flv_rreal(ipart)+flv_rreal(rpart).eq.0).or. &
                      (flv_rreal(ipart).eq.0).or. &
                      (flv_rreal(rpart).eq.0))) cycle
          end if
! ISR:
          if (jpart.lt.3) then
!            print *,"ISR: "
            stop "FindDoublyUnresolved not implemented yet..."
! FSR:
          else
!            print *,"FSR: "
! Is it valid QCD-wise?
            if (.not.((flv_rreal(jpart)+flv_rreal(spart).eq.0).or. &
                      (flv_rreal(jpart).eq.0).or. &
                      (flv_rreal(spart).eq.0))) cycle
          end if
! Is there a corresponding underlying Born?
          call UndoBranchingMod(flv_rreal,flv_tmp,ipart,rpart,jpart,spart)
          if (.not.CheckEquivProc(npart-2,size(flv_born_arr,2), &
                                  flv_born_arr,flv_tmp)) cycle
          radi = ipart
          radr = rpart
          radj = jpart
          rads = spart
          foundregion = .true.
          return 
        end do
        rads = jpart + 1
      end do
      radj = ipart
    end do
    radr = ipart + 2
  end do
!
  foundregion = .false.
!
end subroutine FindDoublyUnresolved_4parton
!
subroutine SplittingCtype(lambda,p_ub,p,flv,radi,radr)
use FKSutils
use particles
implicit none
!
  real(kind(1d0)) , intent(in) :: lambda
  type(mom) , dimension(:) , intent(in) :: p_ub
  type(particle) , dimension(:) , intent(out) :: p
  integer , dimension(:) , intent(in) :: flv
  integer , intent(in) :: radi,radr
!
  integer :: rad1,rad2
  real(kind(1d0)) :: y,xi
  type(mom) :: p_tmp
!
!
! This is a collinear splitting hence y is lambda, the
! energy fraction is picked up randomly:
  y  = 1 - lambda
! Removing artificial cut-off:
  y = y/(1 - tiny_y)
  xi = 0.5d0!xi1
!
! radi and radr should be ordered but for safety measures
! reordering is done:
  rad1 = min(radi,radr)
  rad2 = max(radi,radr)
!
! ISR:
  if (rad1.lt.3) then
! FSR:
  else
    call MapFSR(p_ub,xi,y,azi1,rad1,rad2,flv,p)
  end if
!
! If the positions were swaped we have to swap them back:
  if (rad1.ne.radi) then
    p_tmp     = p(rad1)%p
    p(rad1)%p = p(rad2)%p
    p(rad2)%p = p_tmp
  end if
!
end subroutine SplittingCtype
!
subroutine SplittingStype(lambda,p_ub,p,flv,rads)
use FKSutils
use particles
implicit none
!
  real(kind(1d0)) , intent(in) :: lambda
  type(mom) , dimension(:) , intent(in) :: p_ub
  type(particle) , dimension(:) , intent(out) :: p
  integer , dimension(:) , intent(in) :: flv
  integer , intent(in) :: rads
!
  integer :: radr,rad1,rad2
  real(kind(1d0)) :: y,xi
  type(mom) :: p_tmp
!
!
! This is the soft emission of a gluon, hence lambda is the
! energy fraction while y is picked randomly:
  xi = (lambda - tiny_xi) / (1 - 2*tiny_xi)
  y = 0d0!y1
!
! When a soft radiation is requested only the emitted, soft
! parton is specified, an emitter has to be picked:
  do radr=1,size(flv)
    if (radr.eq.rads) cycle
    if (abs(flv(radr)).le.6) exit
  end do
!
! radr and rads should be ordered but for safety measures
! reordering is done:
  rad1 = min(radr,rads)
  rad2 = max(radr,rads)
!
! ISR:
  if (rad1.lt.3) then
! FSR:
  else
    call MapFSR(p_ub,xi,y,azi1,rad1,rad2,flv,p)
  end if
!
! If the positions were swaped we have to swap them back:
  if (rad1.ne.radr) then
    p_tmp     = p(rad1)%p
    p(rad1)%p = p(rad2)%p
    p(rad2)%p = p_tmp
  end if
!
end subroutine SplittingStype
!
subroutine SplittingCStype(lambda,p_ub,p,flv,radi,radr)
use FKSutils
use particles
implicit none
!
  real(kind(1d0)) , intent(in) :: lambda
  type(mom) , dimension(:) , intent(in) :: p_ub
  type(particle) , dimension(:) , intent(out) :: p
  integer , dimension(:) , intent(in) :: flv
  integer , intent(in) :: radi,radr
!
  integer :: rad1,rad2
  real(kind(1d0)) :: y,xi
  type(mom) :: p_tmp
!
!
! In this limit the collinear and soft region is approached
! in the same time, the soft scaling goes with the square root
! of the collinear scaling:
  y  = 1 - lambda
  xi = sqrt(lambda)
! Removing artificial cut-offs:
  y = y/(1 - tiny_y)
  xi = (xi - tiny_xi) / (1 - 2*tiny_xi)
!
! radi and radr should be ordered but for safety measures
! reordering is done:
  rad1 = min(radi,radr)
  rad2 = max(radi,radr)
!
! ISR:
  if (rad1.lt.3) then
! FSR:
  else
    call MapFSR(p_ub,xi,y,azi1,rad1,rad2,flv,p)
  end if
!
! If the positions were swaped we have to swap them back:
  if (rad1.ne.radi) then
    p_tmp     = p(rad1)%p
    p(rad1)%p = p(rad2)%p
    p(rad2)%p = p_tmp
  end if
!
end subroutine SplittingCStype
!
subroutine SplittingTripleCtype(lambda,p_ub,p,flv,emitirs,radi,radr,rads)
use particles
use random
use math
implicit none
!
  real(kind(1d0)) , intent(in) :: lambda
  type(mom) , dimension(:) , intent(in) :: p_ub
  type(particle) , dimension(:) , intent(out) :: p
  integer , dimension(:) , intent(in) :: flv
  integer , intent(in) :: emitirs,radi,radr,rads
!
  integer :: ipart,jpart
  real(kind(1d0)) :: s_irs
  real(kind(1d0)) :: E,cosgamma,theta,sinphi,cosphi,norm
  real(kind(1d0)) , dimension(3) :: dir
  type(mom) :: p_irs,p_i,p_r,p_s
  type(mom) :: p_tmp
!
!
!
!  write(*,"(A,4(I0,1x))") "emit,i,r,s: ",emitirs,radi,radr,rads
!  call PrintMom(p_ub)
!
! Calculating the overall kinematic invariant and obtaining
! the energy of the emitter:
  p_irs = p_ub(emitirs)
  s_irs = p_irs**2
  E = p_irs%E
! For the splitting we assume complete symmetry and three
! massless decay products that is s_ir = s_is = s_rs = s
! The enclosed angle between any final state particle is 
! \gamma, hence:
  cosgamma = 1 - 1.5*s_irs/E**2
! Due to symmetry the azimuth can be obtained:
  theta = 0.5*acos(1d0/3d0*(4*cosgamma - 1))
  p_i%E = E/3d0
  p_i%px = 0
  p_i%py = E/3d0*sin(theta)
  p_i%pz = E/3d0*cos(theta)
! All the other momenta are obtained from p_i by simply 
! rotating it along the z axis:
  sinphi = sqrt(3d0)/2d0
  cosphi = -0.5d0
  dir    = 0
  dir(3) = -1
  p_r = p_i
  call rotate4vec(dir,sinphi,cosphi,p_r)
  p_s = p_r
  call rotate4vec(dir,sinphi,cosphi,p_s)
! All three momenta are set up around the z axis, an overall
! azimuth is set randomly:
  sinphi = gen_rand()
  cosphi = sqrt(1 - sinphi**2)
  call rotate4vec(dir,sinphi,cosphi,p_i)
  call rotate4vec(dir,sinphi,cosphi,p_r)
  call rotate4vec(dir,sinphi,cosphi,p_s)
! The common axis is the z axis, we have to rotate all three
! vectors such that p_irs becomes the common axis:
  dir(1) = -p_irs%py
  dir(2) = p_irs%px
  dir(3) = 0
  norm = sqrt(dir(1)**2 + dir(2)**2 + dir(3)**2)
  dir = dir/norm
  cosphi = p_irs%pz/sqrt(p_irs%px**2 + p_irs%py**2 + p_irs%pz**2)
  sinphi = sqrt(1 - cosphi**2)
  call rotate4vec(dir,sinphi,cosphi,p_i)
  call rotate4vec(dir,sinphi,cosphi,p_r)
  call rotate4vec(dir,sinphi,cosphi,p_s)
!
  jpart = 0
  do ipart=1,size(p)
    p(ipart)%flv = flv(ipart)
    if ((ipart.eq.radi).or.(ipart.eq.radr).or.(ipart.eq.rads)) then
      if (ipart.eq.radi) p(ipart)%p = p_i
      if (ipart.eq.radr) p(ipart)%p = p_r
      if (ipart.eq.rads) p(ipart)%p = p_s
      if (ipart.eq.emitirs) jpart = jpart + 1
      cycle
    end if
    jpart = jpart + 1
    p(ipart)%p = p_ub(jpart)
  end do
!
end subroutine SplittingTripleCtype
!
subroutine GenerateSinglyUnresLimits(iproc,cont,flv,p_ub,p,phat,ptilde, &
                                     Bij,Bijkl,Bmunuij,Rij, &
                                     radi,radr)
use particles
use regions
use observables
use scales
use utils
implicit none
!
  integer , intent(in) :: iproc
  character(len=2) , intent(in) :: cont
  integer , dimension(:) , intent(in) :: flv
  type(mom) , dimension(:) , intent(inout) :: p_ub
  type(particle) , dimension(:) , intent(inout) :: p,phat,ptilde
  real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
  real(kind(1d0)) , dimension(:,:,:,:) , intent(inout) :: Bijkl
  real(kind(1d0)) , dimension(0:,0:,:,:) , intent(inout) :: Bmunuij
  real(kind(1d0)) , dimension(:,:) , intent(inout) :: Rij
  integer , intent(in) :: radi,radr
!
  logical :: passed
  integer :: ipart,jpart,npart
  integer :: iter,irad,rad1,rad2
  integer :: emit,emitID
  real(kind(1d0)) :: lambda
  real(kind(1d0)) :: s_ir,E_s
  real(kind(1d0)) :: cfunc
  type(particle) , dimension(size(p_ub)) :: parts_ub
!
  real(kind(1d0)) , parameter :: PtMin = 10d0
  real(kind(1d0)) , parameter :: dRMin = 0.4d0
!
  interface
    subroutine Cuts(p,cfunc)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: cfunc
!
    end subroutine Cuts
!
    subroutine sigmaRR(p,weightPS,dsigmaRRi)
    use particles
    implicit none
    !
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(in) :: weightPS
      real(kind(1d0)) , dimension(:) , intent(out) :: dsigmaRRi
    !
    end subroutine sigmaRR
!
    subroutine Calc_mp2_subs(p,phat,ptilde,         &
                             Bij,Bijkl,Bmunuij,Rij, &
                             weightPS,              &
                             Cir,Sr,CirSr,          &
                             Cirs,CSirs,Cirjs,Srs,  &
                             subsRR)
    use regions
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      type(particle) , dimension(:) , intent(out) :: phat,ptilde
      real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
      real(kind(1d0)) , dimension(:,:,:,:) , intent(inout) :: Bijkl
      real(kind(1d0)) , dimension(0:,0:,:,:) , intent(inout) :: Bmunuij
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Rij
      real(kind(1d0)) , intent(in) :: weightPS
      type(subterms) , intent(in) :: Cir,Sr,CirSr
      type(subterms) , intent(in) :: Cirs,CSirs,Cirjs,Srs
      real(kind(1d0)) , dimension(:) , intent(out) :: subsRR
!
    end subroutine Calc_mp2_subs
  end interface
!
  npart = size(p)
!
! Generating an underlying Born configuration to build upon it:
! The generation is such that we restrict to that phase space region
! where the final state particles are well-separated:
  do while (.true.)
    call gen_ranmom(npart-1,p_ub)
    do ipart=3,npart-1
! Checking pts:
      if (get_pt(p_ub(ipart)).lt.PtMin) goto 100
      do jpart=ipart+1,npart-1
! Checking separation:
        if (get_dR(p_ub(ipart),p_ub(jpart)).lt.dRMin) goto 100
      end do
    end do
! If we special cuts, or less stringent it is possible to left without
! any subtractions resulting in termination, to avoid it is desirable to
! check the underlying Born even with the cuts used in the real calculation,
! too:
    call CreateParts(p_ub,flv,parts_ub)
    call Cuts(parts_ub,cfunc)
    if (cfunc.eq.0) goto 100
    exit
100 continue
  end do
! Phase space point is ready:
!  call PrintMom(p_ub)
!
! Determining the position of the emitter and its kind:
  emit = min(radi,radr)
  if (emit.gt.2) then
    emitID = flv(radi) + flv(radr)
  else
  end if
!
! When dealing with singly unresolved regions three different
! cases can emerge: C,S,CS
! C can always, but S and CS only if there is a gluon:
!
!  goto 222
!
! C-case:
  do iter=1,niter
    lambda = 10d0**(-iter)
    call SplittingCtype(lambda,p_ub,p,flv,radi,radr)
    s_ir = (p(radi)%p + p(radr)%p)**2
! RR:
    if (cont.eq.'RR') then
      call sigmaRR(p,1d0,sme)
      call Calc_mp2_subs(p,phat,ptilde, &
                         Bij,Bijkl,Bmunuij,Rij, &
                         1d0, &
                         subterms_Cir_RR(iproc),   &
                         subterms_Sr_RR(iproc),    &
                         subterms_CSir_RR(iproc),  &
                         subterms_Cirs_RR(iproc),  &
                         subterms_CSirs_RR(iproc), &
                         subterms_Cirjs_RR(iproc), &
                         subterms_Srs_RR(iproc),   &
                         subs1)
!
      call CheckValidity(iter,sme,subs1,centscale,passed)
    end if
  end do
  write(*,"(3(3(a),I0),a)",advance='no') &
    "Cir: ",ConvertFromPDG(emitID),"(",emit,") -> ", &
    ConvertFromPDG(p(radi)%flv),"(",radi,") || ", &
    ConvertFromPDG(p(radr)%flv),"(",radr,") "
  if (passed) then
    write(*,*) "VALID"
  else
    write(*,*) "INVALID"
  end if
!
!  goto 111
!  222 continue 
!
! S-type:
  do irad=1,2
    if (irad.eq.1) then
      rad1 = radi
    else
      rad1 = radr
    end if
! Only consider the emitted parton if it is a gluon:
    if (flv(rad1).ne.0) cycle
!    print *,"Soft candidate is found: ",rad1
! Construct phase space:
    do iter=1,niter
      lambda = 10d0**(-iter)
      call SplittingStype(lambda,p_ub,p,flv,rad1)
      E_s = p(rad1)%p%E
!      print *,"E_s: ",E_s
! RR:
      if (cont.eq.'RR') then
        call sigmaRR(p,1d0,sme)
        call Calc_mp2_subs(p,phat,ptilde, &
                           Bij,Bijkl,Bmunuij,Rij, &
                           1d0, &
                           subterms_Cir_RR(iproc),   &
                           subterms_Sr_RR(iproc),    &
                           subterms_CSir_RR(iproc),  &
                           subterms_Cirs_RR(iproc),  &
                           subterms_CSirs_RR(iproc), &
                           subterms_Cirjs_RR(iproc), &
                           subterms_Srs_RR(iproc),   &
                           subs1)
!
        call CheckValidity(iter,sme,subs1,centscale,passed)
      end if
    end do
    write(*,"(3(a),I0,a)",advance='no') &
      "Sr: ",ConvertFromPDG(p(rad1)%flv),"(",rad1,") -> 0 "
    if (passed) then
      write(*,*) "VALID"
    else
      write(*,*) "INVALID"
    end if
  end do
!
!  goto 111
!
! CS-type:
  do irad=1,2
    if (irad.eq.1) then
      rad1 = radi
      rad2 = radr
    else
      rad1 = radr
      rad2 = radi
    end if
! Only consider the case if the second parton (rad2) is a gluon:
    if (flv(rad2).ne.0) cycle
!    print *,"Soft candidate is found: ",rad2
! Construct phase space:
    do iter=1,niter
      lambda = 10d0**(-iter)
      call SplittingCStype(lambda,p_ub,p,flv,rad1,rad2)
      s_ir = (p(radi)%p + p(radr)%p)**2
      E_s = p(rad2)%p%E
!      print *,"s_ir: ",s_ir
!      print *,"E_s: ",E_s
! RR:
      if (cont.eq.'RR') then
        call sigmaRR(p,1d0,sme)
        call Calc_mp2_subs(p,phat,ptilde, &
                           Bij,Bijkl,Bmunuij,Rij, &
                           1d0, &
                           subterms_Cir_RR(iproc),   &
                           subterms_Sr_RR(iproc),    &
                           subterms_CSir_RR(iproc),  &
                           subterms_Cirs_RR(iproc),  &
                           subterms_CSirs_RR(iproc), &
                           subterms_Cirjs_RR(iproc), &
                           subterms_Srs_RR(iproc),   &
                           subs1)
!
        call CheckValidity(iter,sme,subs1,centscale,passed)
      end if
    end do
    write(*,"(4(3(a),I0),a)",advance='no') &
      "CirSr: ",ConvertFromPDG(emitID),"(",emit,") -> ", &
                ConvertFromPDG(p(radi)%flv),"(",radi,") || ", &
                ConvertFromPDG(p(radr)%flv),"(",radr,") , ", &
                ConvertFromPDG(p(rad2)%flv),"(",rad2,") -> 0 "
    if (passed) then
      write(*,*) "VALID"
    else
      write(*,*) "INVALID"
    end if
  end do
!
  111 continue
!
end subroutine GenerateSinglyUnresLimits
!
subroutine GenerateDoublyUnresLimits_2parton(iproc,flv,p_ub,p,phat,ptilde, &
                                             Bij,Bijkl,Bmunuij,Rij, &
                                             radr,rads)
use particles
use regions
use observables
use scales
use utils
implicit none
!
  integer , intent(in) :: iproc
  integer , dimension(:) , intent(in) :: flv
  type(mom) , dimension(:) , intent(inout) :: p_ub
  type(particle) , dimension(:) , intent(inout) :: p,phat,ptilde
  real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
  real(kind(1d0)) , dimension(:,:,:,:) , intent(inout) :: Bijkl
  real(kind(1d0)) , dimension(0:,0:,:,:) , intent(inout) :: Bmunuij
  real(kind(1d0)) , dimension(:,:) , intent(inout) :: Rij
  integer , intent(in) :: radr,rads
!
  logical :: passed
  integer :: ipart,jpart,npart
  integer :: iter
  real(kind(1d0)) :: lambda
  real(kind(1d0)) :: E_r,E_s
  real(kind(1d0)) :: cfunc
  type(particle) , dimension(size(p_ub)) :: parts_ub
!
  real(kind(1d0)) , parameter :: PtMin = 10d0
  real(kind(1d0)) , parameter :: dRMin = 0.4d0
!
  interface
    subroutine Cuts(p,cfunc)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: cfunc
!
    end subroutine Cuts
!
    subroutine sigmaRR(p,weightPS,dsigmaRRi)
    use particles
    implicit none
    !
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(in) :: weightPS
      real(kind(1d0)) , dimension(:) , intent(out) :: dsigmaRRi
    !
    end subroutine sigmaRR
!
    subroutine Calc_mp2_subs(p,phat,ptilde,         &
                             Bij,Bijkl,Bmunuij,Rij, &
                             weightPS,              &
                             Cir,Sr,CirSr,          &
                             Cirs,CSirs,Cirjs,Srs,  &
                             subsRR)
    use regions
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      type(particle) , dimension(:) , intent(out) :: phat,ptilde
      real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
      real(kind(1d0)) , dimension(:,:,:,:) , intent(inout) :: Bijkl
      real(kind(1d0)) , dimension(0:,0:,:,:) , intent(inout) :: Bmunuij
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Rij
      real(kind(1d0)) , intent(in) :: weightPS
      type(subterms) , intent(in) :: Cir,Sr,CirSr
      type(subterms) , intent(in) :: Cirs,CSirs,Cirjs,Srs
      real(kind(1d0)) , dimension(:) , intent(out) :: subsRR
!
    end subroutine Calc_mp2_subs
  end interface
!
  npart = size(p)
!
! Generating an underlying Born configuration to build upon it:
! The generation is such that we restrict to that phase space region
! where the final state particles are well-separated:
  do while (.true.)
    call gen_ranmom(npart-2,p_ub)
    do ipart=3,npart-2
! Checking pts:
      if (get_pt(p_ub(ipart)).lt.PtMin) goto 100
      do jpart=ipart+1,npart-2
! Checking separation:
        if (get_dR(p_ub(ipart),p_ub(jpart)).lt.dRMin) goto 100
      end do
    end do
! If we have special cuts, or less stringent it is possible to left without
! any subtractions resulting in termination, to avoid it is desirable to
! check the underlying Born even with the cuts used in the real calculation,
! too:
    call CreateParts(p_ub,flv,parts_ub)
    call Cuts(parts_ub,cfunc)
    if (cfunc.eq.0) goto 100
    exit
100 continue
  end do
! Phase space point is ready:
!  call PrintMom(p_ub)
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@ Srs @@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  do iter=1,niter
    lambda = 10d0**(-iter)
    lambda = lambda
    call SplittingStype(lambda,p_ub,phat,flv,radr)
    call SplittingStype(lambda,phat(:)%p,p,flv,rads)
    E_r = p(radr)%p%E
    E_s = p(rads)%p%E
!    print *,"E_r: ",E_r
!    print *,"E_s: ",E_s
    call sigmaRR(p,1d0,sme)
    call Calc_mp2_subs(p,phat,ptilde, &
                       Bij,Bijkl,Bmunuij,Rij, &
                       1d0, &
                       subterms_Cir_RR(iproc),   &
                       subterms_Sr_RR(iproc),    &
                       subterms_CSir_RR(iproc),  &
                       subterms_Cirs_RR(iproc),  &
                       subterms_CSirs_RR(iproc), &
                       subterms_Cirjs_RR(iproc), &
                       subterms_Srs_RR(iproc),   &
                       subs1)
!
    call CheckValidity(iter,sme,subs1,centscale,passed)
  end do
  write(*,"(2(3(a),I0),a)",advance='no') &
    "Srs: ",ConvertFromPDG(p(radr)%flv),"(",radr,") -> 0, ", &
    ConvertFromPDG(p(rads)%flv),"(",rads,") -> 0 "
  if (passed) then
    write(*,*) "VALID"
  else
    write(*,*) "INVALID"
  end if
!
end subroutine GenerateDoublyUnresLimits_2parton
!
subroutine GenerateDoublyUnresLimits_3parton(iproc,flv,p_ub,p,phat,ptilde, &
                                             Bij,Bijkl,Bmunuij,Rij, &
                                             radi,radr,rads)
use particles
use regions
use observables
use scales
use utils
implicit none
!
  integer , intent(in) :: iproc
  integer , dimension(:) , intent(in) :: flv
  type(mom) , dimension(:) , intent(inout) :: p_ub
  type(particle) , dimension(:) , intent(inout) :: p,phat,ptilde
  real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
  real(kind(1d0)) , dimension(:,:,:,:) , intent(inout) :: Bijkl
  real(kind(1d0)) , dimension(0:,0:,:,:) , intent(inout) :: Bmunuij
  real(kind(1d0)) , dimension(:,:) , intent(inout) :: Rij
  integer , intent(in) :: radi,radr,rads
!
  logical :: passed
  integer :: ipart,jpart,npart
  integer :: swap
  integer :: iter,irad,rad1,rad2,rad3
  integer :: emit,emitID,emit1,emit2
  real(kind(1d0)) :: lambda
  real(kind(1d0)) :: cfunc
  real(kind(1d0)) :: s_ir,s_irs,E_s
  type(mom) :: p_tmp
  type(mom) , dimension(size(p_ub)) :: pborn
  type(particle) , dimension(size(p_ub)) :: parts_ub
!
  real(kind(1d0)) , parameter :: PtMin = 10d0
  real(kind(1d0)) , parameter :: dRMin = 0.4d0
!
!
  interface
    subroutine Cuts(p,cfunc)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: cfunc
!
    end subroutine Cuts
!
    subroutine sigmaRR(p,weightPS,dsigmaRRi)
    use particles
    implicit none
    !
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(in) :: weightPS
      real(kind(1d0)) , dimension(:) , intent(out) :: dsigmaRRi
    !
    end subroutine sigmaRR
!
    subroutine Calc_mp2_subs(p,phat,ptilde,         &
                             Bij,Bijkl,Bmunuij,Rij, &
                             weightPS,              & 
                             Cir,Sr,CirSr,          &
                             Cirs,CSirs,Cirjs,Srs,  &
                             subsRR)
    use regions
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      type(particle) , dimension(:) , intent(out) :: phat,ptilde
      real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
      real(kind(1d0)) , dimension(:,:,:,:) , intent(inout) :: Bijkl
      real(kind(1d0)) , dimension(0:,0:,:,:) , intent(inout) :: Bmunuij
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Rij
      real(kind(1d0)) , intent(in) :: weightPS
      type(subterms) , intent(in) :: Cir,Sr,CirSr
      type(subterms) , intent(in) :: Cirs,CSirs,Cirjs,Srs
      real(kind(1d0)) , dimension(:) , intent(out) :: subsRR
!
    end subroutine Calc_mp2_subs
  end interface
!
  npart = size(p)
!
! Generating an underlying Born configuration to build upon it:
! The generation is such that we restrict to that phase space region
! where the final state particles are well-separated:
  do while (.true.)
    call gen_ranmom(npart-2,p_ub)
    do ipart=3,npart-2
! Checking pts:
      if (get_pt(p_ub(ipart)).lt.PtMin) goto 100
      do jpart=ipart+1,npart-2
! Checking separation:
        if (get_dR(p_ub(ipart),p_ub(jpart)).lt.dRMin) goto 100
      end do
    end do
! If we have special cuts, or less stringent it is possible to left without
! any subtractions resulting in termination, to avoid it is desirable to
! check the underlying Born even with the cuts used in the real calculation,
! too:
    call CreateParts(p_ub,flv,parts_ub)
    call Cuts(parts_ub,cfunc)
    if (cfunc.eq.0) goto 100
    exit
100 continue
  end do
! Phase space point is ready:
!  call PrintMom(p_ub)
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@ Cirs @@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! Determining the position of the emitter and its kind:
  emit = min(radi,radr,rads)
  if (emit.gt.2) then
    emitID = flv(radi) + flv(radr) + flv(rads)
  else
  end if
  do iter=1,niter
    lambda = 10d0**(-iter)
! We create an auxiliary intermediate PS point:
    call SplittingCtype(lambda,p_ub,phat,flv,radi,radr)
! We replace the momentum at position emit with p(radi) + p(radr):
    pborn = p_ub
    pborn(emit) = phat(radi)%p + phat(radr)%p
    call SplittingTripleCtype(lambda,pborn,p,flv,emit,radi,radr,rads)
    s_irs = (p(radi)%p + p(radr)%p + p(rads)%p)**2
!    print *,"s_irs: ",s_irs
    call sigmaRR(p,1d0,sme)
    call Calc_mp2_subs(p,phat,ptilde, &
                       Bij,Bijkl,Bmunuij,Rij, &
                       1d0, &
                       subterms_Cir_RR(iproc),   &
                       subterms_Sr_RR(iproc),    &
                       subterms_CSir_RR(iproc),  &
                       subterms_Cirs_RR(iproc),  &
                       subterms_CSirs_RR(iproc), &
                       subterms_Cirjs_RR(iproc), &
                       subterms_Srs_RR(iproc),   &
                       subs1)
!
    call CheckValidity(iter,sme,subs1,centscale,passed)
  end do
  write(*,"(4(3(a),I0),a)",advance='no') &
    "Cirs: ",ConvertFromPDG(emitID),"(",emit,") -> ", &
    ConvertFromPDG(p(radi)%flv),"(",radi,") || ", &
    ConvertFromPDG(p(radr)%flv),"(",radr,") || ", &
    ConvertFromPDG(p(rads)%flv),"(",rads,") "
  if (passed) then
    write(*,*) "VALID"
  else
    write(*,*) "INVALID"
  end if
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@ CSirs @@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! For the CSirs limit a soft leg is distinguished, three 
! possible cases can emerge:
  do irad=1,3
! rad3 is the soft leg:
    if (irad.eq.1) then
      if (flv(rads).ne.0) cycle
      rad1 = radi
      rad2 = radr
      rad3 = rads
    elseif (irad.eq.2) then
      if (flv(radr).ne.0) cycle
      rad1 = radi
      rad2 = rads
      rad3 = radr
    else
      if (flv(radi).ne.0) cycle
      rad1 = radr
      rad2 = rads
      rad3 = radi
    end if
! Determining the position of the emitter and its kind:
    emit = min(rad1,rad2)
    if (emit.gt.2) then
      emitID = flv(rad1) + flv(rad2)
    else
    if (emit.gt.rad3) emit = emit - 1
    end if
    do iter=1,niter
      lambda = 10d0**(-iter)
! Unless the soft leg has the highest position we start
! with it:
      if (rad3.ne.rads) then
        call SplittingStype(sqrt(lambda),p_ub,phat,flv,rad3)
        call SplittingCtype(lambda,phat(:)%p,p,flv,rad1,rad2)
      else
        call SplittingCtype(lambda,p_ub,phat,flv,rad1,rad2)
        call SplittingStype(sqrt(lambda),phat(:)%p,p,flv,rad3)
      end if
      s_ir = (p(rad1)%p + p(rad2)%p)**2
      E_s = p(rad3)%p%E
!      print *,"s_ir: ",s_ir
!      print *,"E_s: ",E_s
      call sigmaRR(p,1d0,sme)
      call Calc_mp2_subs(p,phat,ptilde, &
                         Bij,Bijkl,Bmunuij,Rij, &
                         1d0, &
                         subterms_Cir_RR(iproc),   &
                         subterms_Sr_RR(iproc),    &
                         subterms_CSir_RR(iproc),  &
                         subterms_Cirs_RR(iproc),  &
                         subterms_CSirs_RR(iproc), &
                         subterms_Cirjs_RR(iproc), &
                         subterms_Srs_RR(iproc),   &
                         subs1)
!
      call CheckValidity(iter,sme,subs1,centscale,passed)
    end do
    write(*,"(4(3(a),I0),a)",advance='no') &
      "CSirs: ",ConvertFromPDG(emitID),"(",emit,") -> ", &
      ConvertFromPDG(p(rad1)%flv),"(",rad1,") || ", &
      ConvertFromPDG(p(rad2)%flv),"(",rad2,") , ", &
      ConvertFromPDG(p(rad3)%flv),"(",rad3,") -> 0 "
    if (passed) then
      write(*,*) "VALID"
    else
      write(*,*) "INVALID"
    end if
  end do
  goto 111
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@ CirsCSirs @@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! For the CirsCSirs limit a soft leg is distinguished, three 
! possible cases can emerge:
  do irad=1,3
! rad3 is the soft leg:
    if (irad.eq.1) then
      if (flv(rads).ne.0) cycle
      rad1 = radi
      rad2 = radr
      rad3 = rads
    elseif (irad.eq.2) then
      if (flv(radr).ne.0) cycle
      rad1 = radi
      rad2 = rads
      rad3 = radr
    else
      if (flv(radi).ne.0) cycle
      rad1 = radr
      rad2 = rads
      rad3 = radi
    end if
! Determining the position of the emitter and its kind:
    emit = min(rad1,rad2,rad3)
    if (emit.gt.2) then
      emitID = flv(rad1) + flv(rad2) + flv(rad3)
    else
    end if
    do iter=1,niter
      lambda = 10d0**(-iter)
! Unless the soft leg has the highest position we start
! with it:
      if (rad3.ne.rads) then
        call SplittingCStype(lambda,p_ub,phat,flv,rad1,rad3)
        call SplittingCtype(lambda,phat(:)%p,p,flv,rad1,rad2)
      else
        call SplittingCtype(lambda,p_ub,phat,flv,rad1,rad2)
        call SplittingCStype(lambda,phat(:)%p,p,flv,rad1,rad3)
      end if
      s_irs = (p(radi)%p + p(radr)%p + p(rads)%p)**2
      E_s = p(rad3)%p%E
!      print *,"s_irs: ",s_irs
!      print *,"E_s: ",E_s
      call sigmaRR(p,1d0,sme)
      call Calc_mp2_subs(p,phat,ptilde, &
                         Bij,Bijkl,Bmunuij,Rij, &
                         1d0, &
                         subterms_Cir_RR(iproc),   &
                         subterms_Sr_RR(iproc),    &
                         subterms_CSir_RR(iproc),  &
                         subterms_Cirs_RR(iproc),  &
                         subterms_CSirs_RR(iproc), &
                         subterms_Cirjs_RR(iproc), &
                         subterms_Srs_RR(iproc),   &
                         subs1)
!
      call CheckValidity(iter,sme,subs1,centscale,passed)
    end do
    write(*,"(5(3(a),I0),a)",advance='no') &
      "CirsCSirs: ",ConvertFromPDG(emitID),"(",emit,") -> ", &
      ConvertFromPDG(p(radi)%flv),"(",radi,") || ", &
      ConvertFromPDG(p(radr)%flv),"(",radr,") || ", &
      ConvertFromPDG(p(rads)%flv),"(",rads,") , ", &
      ConvertFromPDG(p(rad3)%flv),"(",rad3,") -> 0 "
    if (passed) then
      write(*,*) "VALID"
    else
      write(*,*) "INVALID"
    end if
  end do
  111 continue
!
end subroutine GenerateDoublyUnresLimits_3parton
!
subroutine GenerateDoublyUnresLimits_4parton(iproc,flv,p_ub,p,phat,ptilde, &
                                             Bij,Bijkl,Bmunuij,Rij, &
                                             radi,radr,radj,rads)
use particles
use regions
use observables
use scales
use utils
implicit none
!
  integer , intent(in) :: iproc
  integer , dimension(:) , intent(in) :: flv
  type(mom) , dimension(:) , intent(inout) :: p_ub
  type(particle) , dimension(:) , intent(inout) :: p,phat,ptilde
  real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
  real(kind(1d0)) , dimension(:,:,:,:) , intent(inout) :: Bijkl
  real(kind(1d0)) , dimension(0:,0:,:,:) , intent(inout) :: Bmunuij
  real(kind(1d0)) , dimension(:,:) , intent(inout) :: Rij
  integer , intent(in) :: radi,radr,radj,rads
!
  logical :: passed
  integer :: ipart,jpart,npart
  integer :: iter,irad,rad1,rad2,rad3,rad4
  integer :: emitir,emitirID
  integer :: emitjs,emitjsID
  real(kind(1d0)) :: lambda
  real(kind(1d0)) :: cfunc
  real(kind(1d0)) :: s_ir,s_js
  type(particle) , dimension(size(p_ub)) :: parts_ub
!
  real(kind(1d0)) , parameter :: PtMin = 10d0
  real(kind(1d0)) , parameter :: dRMin = 0.4d0
!
!
  interface
    subroutine Cuts(p,cfunc)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: cfunc
!
    end subroutine Cuts
!
    subroutine sigmaRR(p,weightPS,dsigmaRRi)
    use particles
    implicit none
    !
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(in) :: weightPS
      real(kind(1d0)) , dimension(:) , intent(out) :: dsigmaRRi
    !
    end subroutine sigmaRR
!
    subroutine Calc_mp2_subs(p,phat,ptilde,         &
                             Bij,Bijkl,Bmunuij,Rij, &
                             weightPS,              & 
                             Cir,Sr,CirSr,          &
                             Cirs,CSirs,Cirjs,Srs,  &
                             subsRR)
    use regions
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      type(particle) , dimension(:) , intent(out) :: phat,ptilde
      real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
      real(kind(1d0)) , dimension(:,:,:,:) , intent(inout) :: Bijkl
      real(kind(1d0)) , dimension(0:,0:,:,:) , intent(inout) :: Bmunuij
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Rij
      real(kind(1d0)) , intent(in) :: weightPS
      type(subterms) , intent(in) :: Cir,Sr,CirSr
      type(subterms) , intent(in) :: Cirs,CSirs,Cirjs,Srs
      real(kind(1d0)) , dimension(:) , intent(out) :: subsRR
!
    end subroutine Calc_mp2_subs
  end interface
!
  npart = size(p)
!
! Generating an underlying Born configuration to build upon it:
! The generation is such that we restrict to that phase space region
! where the final state particles are well-separated:
  do while (.true.)
    call gen_ranmom(npart-2,p_ub)
    do ipart=3,npart-2
! Checking pts:
      if (get_pt(p_ub(ipart)).lt.PtMin) goto 100
      do jpart=ipart+1,npart-2
! Checking separation:
        if (get_dR(p_ub(ipart),p_ub(jpart)).lt.dRMin) goto 100
      end do
    end do
! If we have special cuts, or less stringent it is possible to left without
! any subtractions resulting in termination, to avoid it is desirable to
! check the underlying Born even with the cuts used in the real calculation,
! too:
    call CreateParts(p_ub,flv,parts_ub)
    call Cuts(parts_ub,cfunc)
    if (cfunc.eq.0) goto 100
    exit
100 continue
  end do
! Phase space point is ready:
!  call PrintMom(p_ub)
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@ Cirjs @@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! Determining the position of the emitters:
  emitir = min(radi,radr)
  emitjs = min(radj,rads)
  if (max(radi,radr).lt.max(radj,rads)) then
    rad1 = radi
    rad2 = radr
    rad3 = radj
    rad4 = rads
  else
    rad1 = radj
    rad2 = rads
    rad3 = radi
    rad4 = radr
  end if
! The position can change:
  if (emitir.gt.max(radj,rads)) emitir = emitir - 1
  if (emitjs.gt.max(radi,radr)) emitjs = emitjs - 1
  if (emitir.gt.2) then
    emitirID = flv(radi) + flv(radr)
  else
  end if
  if (emitjs.gt.2) then
    emitjsID = flv(radj) + flv(rads)
  else
  end if
  do iter=1,niter
    lambda = 10d0**(-iter)
    call SplittingCtype(lambda,p_ub,phat,flv,rad1,rad2)
    call SplittingCtype(lambda,phat(:)%p,p,flv,rad3,rad4)
    s_ir = (p(radi)%p + p(radr)%p)**2
    s_js = (p(radj)%p + p(rads)%p)**2
!    print *,"s_ir, s_js: ",s_ir,s_js
    call sigmaRR(p,1d0,sme)
    call Calc_mp2_subs(p,phat,ptilde, &
                       Bij,Bijkl,Bmunuij,Rij, &
                       1d0, &
                       subterms_Cir_RR(iproc),   &
                       subterms_Sr_RR(iproc),    &
                       subterms_CSir_RR(iproc),  &
                       subterms_Cirs_RR(iproc),  &
                       subterms_CSirs_RR(iproc), &
                       subterms_Cirjs_RR(iproc), &
                       subterms_Srs_RR(iproc),   &
                       subs1)
!
    call CheckValidity(iter,sme,subs1,centscale,passed)
  end do
  write(*,"(6(3(a),I0),a)",advance='no') &
    "Cirjs: ",ConvertFromPDG(emitirID),"(",emitir,") -> ", &
    ConvertFromPDG(p(radi)%flv),"(",radi,") || ", &
    ConvertFromPDG(p(radr)%flv),"(",radr,") , ", &
    ConvertFromPDG(emitjsID),"(",emitjs,") -> ", &
    ConvertFromPDG(p(radj)%flv),"(",radj,") || ", &
    ConvertFromPDG(p(rads)%flv),"(",rads,") "
  if (passed) then
    write(*,*) "VALID"
  else
    write(*,*) "INVALID"
  end if
!
end subroutine GenerateDoublyUnresLimits_4parton
!
subroutine ReorderRads(flv,radi,radr,rads,emit1,rad1,emit2,rad2,swap)
implicit none
!
  integer , dimension(:) , intent(in) :: flv
  integer , intent(in) :: radi,radr,rads
  integer , intent(out) :: emit1,emit2,rad1,rad2
  integer , intent(out) :: swap
!
  integer :: irad
!
!
  print *,"i,r,s: ",radi,radr,rads
  print *,flv(radi),flv(radr),flv(rads)
!
  swap = 0
!
! Only reorder if quarks are present:
  if (min(flv(radi),flv(radr),flv(rads)).ne. &
      max(flv(radi),flv(radr),flv(rads))) then
    emit1 = min(radi,radr,rads)
! Finding the gluon or quark pair:
    if (flv(radr)+flv(rads).eq.0) then
      rad1  = max(radi,min(radr,rads))
      emit2 = min(radr,rads)
      rad2  = max(radr,rads)
    elseif (flv(radi)+flv(rads).eq.0) then
      rad1  = max(radr,min(radi,rads))
      emit2 = min(radi,rads)
      rad2  = max(radi,rads)
    elseif (flv(radi)+flv(radr).eq.0) then
      rad1  = max(rads,min(radi,radr))
      emit2 = min(radi,radr)
      rad2  = max(radi,radr)
    end if
! It is possible that the emitter of the first branching
! corresponds to the last parton in that case a swap is
! needed:
    if (rad1.eq.size(flv)) then
      do irad=3,size(flv)-1
        if ((irad.eq.radi).or.(irad.eq.radr).or.(irad.eq.rads)) cycle
        exit
      end do
      swap = rad1
      rad1 = irad
    end if
    print *,"emit1,rad1: ",emit1,rad1
    print *,"emit2,rad2: ",emit2,rad2
    print *,"swap : ",swap
  else
    emit1 = radi
    rad1  = radr
    emit2 = radr
    rad2  = rads
  end if
!
  read(*,*)
!
end subroutine ReorderRads
!
end module RegCheck
!
subroutine RegionCheck
use process
use flags
use subprocesses
use particles
use RegCheck
use random
use math
use scales
implicit none
!
!
  logical :: foundregion
  integer :: istat
  integer :: iproc
  integer :: radi,radr,radj,rads
  real(kind(1d0)) , dimension(:,:) , allocatable :: Bij,Rij
  real(kind(1d0)) , dimension(:,:,:,:) , allocatable :: Bijkl
  real(kind(1d0)) , dimension(:,:,:,:) , allocatable :: Bmunuij
  type(mom) , dimension(:) , allocatable :: pborn,preal
  type(particle) , dimension(:) , allocatable :: p_r,p_rr,ptilde,phat
!
!
!
  write(*,*)
  write(*,*) "*********************************************************"
  write(*,*) "********************** RegionCheck **********************"
  write(*,*) "***************** Checking regions for ******************"
  write(*,*) "************ proper subtraction term usage **************"
  write(*,*) "*********************************************************"
!
! If there is no NLO/NNLO R* piece activated there is nothing to do:
  if (.not.(flg_NLO_R.or.flg_NNLO_RV.or.flg_NNLO_RR)) return
!
! Allocating auxiliary arrays to hold particles:
  allocate(pborn(nleg_born),   &
           preal(nleg_born+1), &
           p_r(nleg_born+1),   &
           p_rr(nleg_born+2),  &
           ptilde(nleg_born),  &
           phat(nleg_born+1),  &
           stat=istat)
  if (istat.ne.0) then
    print *,"Allocation for pborn,preal,p_r,p_rr,ptilde and phat"
    print *,"went wrong in RegionCheck..."
    stop "RegionCheck"
  end if
!
! Allocating auxiliary arrays to hold SMEs and subtractions:
  allocate(sme(nscales),   &
           subs1(nscales), &
           stat=istat)
  if (istat.ne.0) then
    print *,"Allocation for subs* went wrong in RegionCheck..."
    stop "RegionCheck"
  end if
!
! Allocating arrays holding color-correlated SMEs:
  allocate(Bij(nleg_born,nleg_born),                       &
           Bijkl(nleg_born,nleg_born,nleg_born,nleg_born), &
           Bmunuij(0:3,0:3,nleg_born,nleg_born),           &
           Rij(nleg_born+1,nleg_born+1),                   &
           stat=istat)
  if (istat.ne.0) then
    print *,"Allocation for color-corr. SMEs went wrong in RegionCheck..."
    stop "RegionCheck"
  end if
!
! Giving random values for xi and y variables:
  xi1  = gen_rand()
  xi2  = gen_rand()
  y1   = gen_rand()
  y2   = gen_rand()
  azi1 = 2*pi*gen_rand()
  azi2 = 2*pi*gen_rand()
!
! Checking the entire RR line:
  if (flg_NNLO_RR) then
! Loop over all subprocesses:
    do iproc=1,num_flv_irr_NNLO_RR
      if (checksingly) then
        write(*,'(A)',advance='no') "Checking singly unresolved regions for "
        call PrintSubproc(flv_ch_RRkin(:,iproc))
! Giving initial values for radi radr:
        radi = 1
        radr = 1
        do while (.true.)
! Trying to find a region:
          call FindSinglyUnresolved(flv_NLO_R, &
                                    flv_ch_RRkin(:,iproc), &
                                    radi,radr,foundregion)
! If not found anything just jump out:
          if (.not.foundregion) exit
!          print *,"Region is found with : ",radi,radr
! Region found work out all possible limits:
! This is done for the RR contribution, the region is singly
! unresolved hence the underlying Born has Real-type kinematics:
          call GenerateSinglyUnresLimits(iproc,"RR",flv_ch_RRkin(:,iproc), &
                                         preal,p_rr,phat,ptilde,           &
                                         Bij,Bijkl,Bmunuij,Rij,            &
                                         radi,radr)
        end do
      end if
      if (checkdoubly) then
        write(*,'(A)',advance='no') "Checking doubly unresolved regions for "
        call PrintSubproc(flv_ch_RRkin(:,iproc))
! Giving initial values for radi radr and rads:
        radi = 1
        radr = 2
        rads = 2
        do while (.true.)
! Trying to find a region:
          call FindDoublyUnresolved(flv_LO, &
                                    flv_ch_RRkin(:,iproc), &
                                    radi,radr,rads,foundregion)
! If not found anything just jump out:
          if (.not.foundregion) exit
!          print *,"Region is found with : ",radi,radr,rads
! Region found work out all possible limits:
          call GenerateDoublyUnresLimits(iproc,flv_ch_RRkin(:,iproc), &
                                         pborn,p_rr,phat,ptilde,      &
                                         Bij,Bijkl,Bmunuij,Rij,       &
                                         radi,radr,rads)
        end do
! Giving initial values for radi, radr, radj and rads:
        radi = 1
        radr = 2
        radj = 2
        rads = 2
        do while (.true.)
! Trying to find a region:
          call FindDoublyUnresolved(flv_LO, &
                                    flv_ch_RRkin(:,iproc), &
                                    radi,radr,radj,rads,foundregion)
! If not found anything just jump out:
          if (.not.foundregion) exit
!          print *,"Region is found with : ",radi,radr,radj,rads
! Region found work out all possible limits:
          call GenerateDoublyUnresLimits(iproc,flv_ch_RRkin(:,iproc), &
                                         pborn,p_rr,phat,ptilde,      &
                                         Bij,Bijkl,Bmunuij,Rij,       &
                                         radi,radr,radj,rads)
        end do
! Giving initial values for radr and rads:
        radr = 1
        rads = 1
        do while (.true.)
! Trying to find a region:
          call FindDoublyUnresolved(flv_LO,                &
                                    flv_ch_RRkin(:,iproc), &
                                    radr,rads,foundregion)
! If not found anything just jump out:
          if (.not.foundregion) exit
!          print *,"Region is found with : ",radr,rads
! Region found work out all possible limits:
          call GenerateDoublyUnresLimits(iproc,flv_ch_RRkin(:,iproc), &
                                         pborn,p_rr,phat,ptilde,      &
                                         Bij,Bijkl,Bmunuij,Rij,       &
                                         radr,rads)
        end do
      end if
    end do
  end if
!
  deallocate(pborn,preal,p_r,p_rr,ptilde,phat)
  deallocate(sme,subs1)
  deallocate(Bij,Bijkl,Bmunuij,Rij)
!
  stop "RegionCheck"
!
end subroutine RegionCheck
!
! This routine can be used to produce ratio plots for the whole 
! partonic lines: lines containing m+1 and m+2 partons.
! For the m+2 line the possbile types of regions are:
! Cirs,CSirs,Cirjs,Srs;
! Cir,Sr,CirSr
subroutine MakeRatioPlot
use particles
use RegCheck
use process
use flags
use subprocesses
use scales
use regions
use observables
implicit none
!
!
  integer , parameter :: iproc = 4
  integer , parameter :: npoint = 10000
  character(len=3) , parameter :: row = 'm+2'
  character(len=5) , parameter :: limit = 'Sr   '
  real(kind(1d0)) , parameter :: lambda = 1d-5
  integer , parameter :: radi = 5
  integer , parameter :: radr = 6
  integer , parameter :: radj = 4
  integer , parameter :: rads = 7
!
  integer , parameter :: iun = 99
  integer , parameter :: nbin = 50
  real(kind(1d0)) , parameter :: xmin = 0.999d0
  real(kind(1d0)) , parameter :: xmax = 1.001d0
  character (len=128) :: fname = 'scatter'
!
  real(kind(1d0)) , parameter :: PtMin = 10d0
  real(kind(1d0)) , parameter :: dRMin = 0.4d0
!
  integer :: istat,ipoint,ibin,ipart,jpart
  integer :: emitirs
  integer , dimension(:) , allocatable :: flv
  real(kind(1d0)) :: cfunc
  real(kind(1d0)) :: s_irs,s_ir,s_js,E_r,E_s
  real(kind(1d0)) :: dx,xpos,ratio
  integer , dimension(:) , allocatable :: histo
  real(kind(1d0)) , dimension(:) , allocatable :: histobins
  type(mom) , dimension(:) , allocatable :: pborn,preal
  type(particle) , dimension(:) , allocatable :: parts_born,parts_real
  type(particle) , dimension(:) , allocatable :: p,phat,ptilde
  real(kind(1d0)) , dimension(:,:) , allocatable :: Bij
  real(kind(1d0)) , dimension(:,:,:,:) , allocatable :: Bijkl
  real(kind(1d0)) , dimension(:,:,:,:) , allocatable :: Bmunuij
  real(kind(1d0)) , dimension(:,:) , allocatable :: Rij
!
!
  interface
    subroutine Cuts(p,cfunc)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: cfunc
!
    end subroutine Cuts
!
    subroutine sigmaRR(p,weightPS,dsigmaRRi)
    use particles
    implicit none
    !
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(in) :: weightPS
      real(kind(1d0)) , dimension(:) , intent(out) :: dsigmaRRi
    !
    end subroutine sigmaRR
!
    subroutine Calc_mp2_subs(p,phat,ptilde,         &
                             Bij,Bijkl,Bmunuij,Rij, &
                             weightPS,              & 
                             Cir,Sr,CirSr,          &
                             Cirs,CSirs,Cirjs,Srs,  &
                             subsRR)
    use regions
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      type(particle) , dimension(:) , intent(out) :: phat,ptilde
      real(kind(1d0)) , dimension(:,:) , intent(inout) :: Bij
      real(kind(1d0)) , dimension(:,:,:,:) , intent(inout) :: Bijkl
      real(kind(1d0)) , dimension(0:,0:,:,:) , intent(inout) :: Bmunuij
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Rij
      real(kind(1d0)) , intent(in) :: weightPS
      type(subterms) , intent(in) :: Cir,Sr,CirSr
      type(subterms) , intent(in) :: Cirs,CSirs,Cirjs,Srs
      real(kind(1d0)) , dimension(:) , intent(out) :: subsRR
!
    end subroutine Calc_mp2_subs
  end interface
!
! Allocation of all momenta arrays needed by the Calc_mp2_subs 
! and the Calc_mp1_subs routines:
  allocate(pborn(nleg_born),        &
           preal(nleg_born+1),      &
           parts_born(nleg_born),   &
           parts_real(nleg_born+1), &
           ptilde(nleg_born),       &
           phat(nleg_born+1),       &
           p(nleg_born+2),          &
! allocating arrays to hold smes and subtractions:
           sme(nscales),   &
           subs1(nscales), &
! allocating arrays holding color-correlated amplitudes:
           Bij(nleg_born,nleg_born),                       &
           Bijkl(nleg_born,nleg_born,nleg_born,nleg_born), &
           Bmunuij(0:3,0:3,nleg_born,nleg_born),           &
           Rij(nleg_born+1,nleg_born+1),                   &
           stat=istat)
!
  if (istat.ne.0) then
    print *,"Problem occured with allocation in MakeRatioPlotLines"
    stop "MakeRatioPlotLines"
  end if
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
! Selection of the subprocess:
  if (row.eq.'m+2') then
    allocate(flv(nleg_born+2))
! Filling up the flv array:
    flv = flv_ch_RRkin(:,iproc)
    call PrintSubProc(flv)
  elseif (row.eq.'m+1') then
    allocate(flv(nleg_born+1))
    flv = flv_ch_Rkin(:,iproc)
    call PrintSubProc(flv)
  else
    print *,"Wrong row is selected..."
    print *,"row: ",row
    stop "MakeRatioPlotLines"
  end if
  do ipoint=1,npoint
    if (row.eq.'m+2') then
! Generating the underlying process, this can have Born or
! even real kinematics depending upon which limit is intended
! to be taken:
      if ((limit.eq.'Cirs ').or. &
          (limit.eq.'Cirjs').or. &
          (limit.eq.'CSirs').or. &
          (limit.eq.'Srs  ')) then
        do while (.true.)
          call gen_ranmom(nleg_born,pborn)
          call CreateParts(pborn,flv,parts_born)
          call Cuts(parts_born,cfunc)
          if (cfunc.ne.0) then
            exit
          end if
        end do
!        print *,"Generated underlying Born PS point: "
!        call PrintMom(pborn)
! Generation of the RR kinematics:
        if (limit.eq.'Cirs ') then
          emitirs = min(radi,radr)
! Generation of an auxiliary PS point:
          call SplittingCtype(lambda,pborn,phat,flv,radi,radr)
! Replacing massless momentum at the place of emitter to the
! recently created massive one (the sum of the two daughter 
! partons):
          pborn(emitirs) = phat(radi)%p + phat(radr)%p
          call SplittingTripleCtype(lambda,pborn,p,flv,emitirs,radi,radr,rads)
          s_irs = (p(radi)%p + p(radr)%p + p(rads)%p)**2
          Q2 = (p(1)%p + p(2)%p)**2
          if (ipoint.lt.10) then
            print *,"s_irs: ",s_irs
            print *,"y_irs: ",s_irs/Q2
          end if
        elseif (limit.eq.'CSirs') then
          call SplittingCtype(lambda,pborn,phat,flv,radi,radr)
          call SplittingStype(sqrt(lambda),phat(:)%p,p,flv,rads)
          Q2 = (p(1)%p + p(2)%p)**2
          s_ir = (p(radi)%p + p(radr)%p)**2
          E_s = p(rads)%p%E
          if (ipoint.lt.10) then
            print *,"s_ir: ",s_ir
            print *,"y_ir: ",s_ir/Q2
            print *,"E_s: ",E_s
            print *,"e_s: ",E_s/sqrt(Q2)
          end if
        elseif (limit.eq.'Cirjs') then
          call SplittingCtype(lambda,pborn,phat,flv,radi,radr)
          call SplittingCtype(lambda,phat(:)%p,p,flv,radj,rads)
          Q2 = (p(1)%p + p(2)%p)**2
          s_ir = (p(radi)%p + p(radr)%p)**2
          s_js = (p(radj)%p + p(rads)%p)**2
          if (ipoint.lt.10) then
            print *,"s_ir: ",s_ir
            print *,"s_js: ",s_js
            print *,"y_ir: ",s_ir/Q2
            print *,"y_js: ",s_js/Q2
          end if
        elseif (limit.eq.'Srs  ') then
          call SplittingStype(lambda,pborn,phat,flv,radr)
          call SplittingStype(lambda,phat(:)%p,p,flv,rads)
          Q2 = (p(1)%p + p(2)%p)**2
          E_r = p(radr)%p%E
          E_s = p(rads)%p%E
          if (ipoint.lt.10) then
            print *,"E_r: ",E_r
            print *,"E_s: ",E_s
            print *,"e_r: ",E_r/sqrt(Q2)
            print *,"e_s: ",E_s/sqrt(Q2)
          end if
        else
          print *,"Wrong limit is specified..."
          print *,"limit: ",limit
          stop "MakeRatioPlot"
        end if
        call sigmaRR(p,1d0,sme)
        call Calc_mp2_subs(p,phat,ptilde, &
                           Bij,Bijkl,Bmunuij,Rij, &
                           1d0, &
                           subterms_Cir_RR(iproc),   &
                           subterms_Sr_RR(iproc),    &
                           subterms_CSir_RR(iproc),  &
                           subterms_Cirs_RR(iproc),  &
                           subterms_CSirs_RR(iproc), &
                           subterms_Cirjs_RR(iproc), &
                           subterms_Srs_RR(iproc),   &
                           subs1)
        ratio = subs1(centscale)/sme(centscale)
      elseif ((limit.eq.'Cir  ').or. &
              (limit.eq.'Sr   ').or. &
              (limit.eq.'CirSr')) then
! When generating PS points in the singly unresolved region we make
! an attempt to really generate point there, to this end we cannot
! use the cut function since it has to be set such that the A_2 and
! A_{12} terms give contribution too:
        do while (.true.)
          call gen_ranmom(nleg_born+1,preal)
          do ipart=3,nleg_born+1
! Checking pts:
            if (get_pt(preal(ipart)).lt.PtMin) goto 100
            do jpart=ipart+1,nleg_born+1
! Checking separation:
              if (get_dR(preal(ipart),preal(jpart)).lt.dRMin) goto 100
            end do
          end do
          call CreateParts(preal,flv,parts_real)
          call Cuts(parts_real,cfunc)
          if (cfunc.ne.0) then
            exit
          end if
          100 continue
        end do
!        print *,"Generated underlying Real PS point: "
!        call PrintMom(preal)
! Generation of the RR kinematics:
        if (limit.eq.'Cir  ') then
          call SplittingCtype(lambda,preal,p,flv,radi,radr)
          Q2 = (p(1)%p + p(2)%p)**2
          s_ir = (p(radi)%p + p(radr)%p)**2
          if (ipoint.lt.10) then
            print *,"s_ir: ",s_ir
            print *,"y_ir: ",s_ir/Q2
          end if
        elseif (limit.eq.'Sr   ') then
          call SplittingStype(lambda,preal,p,flv,radr)
          Q2 = (p(1)%p + p(2)%p)**2
          E_r = p(radr)%p%E
          if (ipoint.lt.10) then
            print *,"E_r: ",E_r
            print *,"e_r: ",E_r/sqrt(Q2)
          end if
        elseif (limit.eq.'CirSr') then
          call SplittingCStype(lambda,preal,p,flv,radi,radr)
          Q2 = (p(1)%p + p(2)%p)**2
          s_ir = (p(radi)%p + p(radr)%p)**2
          E_r = p(radr)%p%E
          if (ipoint.lt.10) then
            print *,"s_ir: ",s_ir
            print *,"y_ir: ",s_ir/Q2
            print *,"E_r: ",E_r
            print *,"e_r: ",E_r/sqrt(Q2)
          end if
        else
          print *,"Wrong limit is specified..."
          print *,"limit: ",limit
          stop "MakeRatioPlot"
        end if
        call sigmaRR(p,1d0,sme)
        call Calc_mp2_subs(p,phat,ptilde, &
                           Bij,Bijkl,Bmunuij,Rij, &
                           1d0, &
                           subterms_Cir_RR(iproc),   &
                           subterms_Sr_RR(iproc),    &
                           subterms_CSir_RR(iproc),  &
                           subterms_Cirs_RR(iproc),  &
                           subterms_CSirs_RR(iproc), &
                           subterms_Cirjs_RR(iproc), &
                           subterms_Srs_RR(iproc),   &
                           subs1)
        ratio = subs1(centscale)/sme(centscale)
      else
        print *,"Wrong limit is specified..."
        print *,"row: ",row
        print *,"limit: ",limit
        stop "MakeRatioPlot"
      end if
    end if
    do ibin=1,nbin
      if ((ratio.gt.histobins(ibin-1)).and. &
          (ratio.lt.histobins(ibin))) then
        histo(ibin) = histo(ibin) + 1
      end if
    end do
  end do
!
  fname = trim(fname)//"-"//trim(limit)//".dat"
  open(unit=iun,file=fname,status='unknown')
!
  do ibin=1,nbin
    xpos = (histobins(ibin-1) + histobins(ibin))/2d0
    write(iun,fmt="(E16.6,4x,I0)") xpos,histo(ibin)
  end do
!
  deallocate(pborn,preal,ptilde,phat,p, &
             parts_born,parts_real,     &
             sme,subs1,                 &
             Bij,Bijkl,Bmunuij,Rij,     &
             flv,                       &
             histo,histobins            &
            )
!
  close(iun)
!
  stop "MakeRatioPlot"
!
end subroutine MakeRatioPlot
