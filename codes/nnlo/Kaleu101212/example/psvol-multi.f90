program montecarlo
  use particles
  use cuts_module
  use rambo_module
  implicit none
  type(rambo_type) :: rambo
  type(cuts_type) :: cuts
  integer :: iproc,nproc
  integer , allocatable :: process(:,:),nfinst(:)
  real(kind(1d0)) , allocatable , dimension(:) :: s0,s1,s2,se, &
                                                  ave,sig,rel, &
                                                  s0opt,seopt
!
! The shape (-2:17) is hard-wired in Kaleu. (-2,-1) are the initial states,
! (1,2,...) are the final states, the maximal number being 17. 0 is nothing.
! The same holds for the other arrays below with this shape.
!
  real(kind(1d0)) ,allocatable :: masses(:),widths(:)
  character(20)   ,allocatable :: labels(:)
  logical         :: onlyqcd,withqcd,withhiggs
  real(kind(1d0)) :: smin(-2:17,-2:17),ecm
  real(kind(1d0)) :: pkaleu(0:3,-2:17),x1,x2,weight,thrs
!
! As you can see above, momentum components are in the (0,1,2,3) convention.
!
  integer :: itmp(99),nev,nbatch,nstep,option,iev,noptim,igenerator
  logical :: discard,next,optim
  character(5) :: start
!
  real(kind(1d0)) , parameter :: pi = 3.14159265358979323846264338327d0
!
! Initialize particles
  call init_particles
!
! Skip header comments in input
  next = .true.
  do while (next)
    read(5,'(a5)') start
    next = start.ne.'start'
  enddo
!
! We have multiply processes:
  read(5,*) nproc
  allocate(process(-2:17,nproc))
  allocate(nfinst(nproc))
  allocate(s0(nproc),s1(nproc),s2(nproc),se(nproc))
  allocate(ave(nproc),sig(nproc),rel(nproc))
  allocate(s0opt(nproc),seopt(nproc))
!
! Read process. Numbering defined in particles.f90
  do iproc=1,nproc
    read(5,*) nfinst(iproc)
    read(5,*) itmp(1:nfinst(iproc)+2)
    process(-2:-1,iproc) = itmp(1:2)
    process( 1:nfinst(iproc),iproc) = itmp(3:nfinst(iproc)+2)
  end do
  print *,"The processes which are read in:"
  do iproc=1,nproc
    print *,"iproc: ",iproc
    print *,process(-2:nfinst(iproc),iproc)
  end do
!
! Read some options
  read(5,*) onlyqcd,withqcd,withhiggs
!
! Read center-of-mass energy
  read(5,*) ecm
!
! Read the number of events to be generated
  read(5,*) nev
!
! Read which phase space generator to use
  read(5,*) igenerator
!
! Read optimization parameters for Kaleu
  read(5,*) nbatch ! number of accepted ps-points before optimization step
  read(5,*) nstep  ! number of optimization steps
  read(5,*) thrs   ! threshold for channel-removal (after full optimization)
!
! Initialize particle labels, masses and widths as defined in "constants.h90"
  allocate(labels(nfirst:nlast)) ! "nfirst" and "nlast" are defined in
  allocate(masses(nfirst:nlast)) ! the module "particles"
  allocate(widths(nfirst:nlast)) !
  include 'constants.h90'
!
! Initialize cuts, and get smin to be delivered to Kaleu
! This should be modified...
!  call cuts_init(cuts,process(:,1),nfinst(1),masses,smin)
  smin = 1d0
!
! No pdfs in this example
  option = 0
!
! Initialize phase space generator
  do iproc=1,nproc
    call multi_kaleu_init(iproc,process(:,iproc),nfinst(iproc),ecm,option, &
                          masses,widths,labels, &
                          onlyqcd,withqcd,withhiggs, &
                          nbatch,nstep,thrs,smin)
  end do
!
! Number of events to be spent in optimization phase
  noptim = nbatch*nstep
  optim  = .true.
!
! Initialize Monte Carlo
  s0 = 0d0 ! number of generated events
  s1 = 0d0 ! sum of full weights
  s2 = 0d0 ! sum of squares of full weights
  se = 0d0 ! number of accepted events
!
  do iev=1,nev
! run through all the subprocesses:
    do iproc=1,nproc
      call multi_kaleu_gnrt(iproc,discard,x1,x2,pkaleu)
!
      if (discard) then
        weight = 0d0
        cycle
      end if
!  
      call multi_kaleu_wght(iproc,weight)
! We have to include a tower of 2*pis:
      weight = weight * (2d0*pi)**(4 - 3*nfinst(iproc))
      se(iproc) = se(iproc) + 1d0
!  
! The weight should be given back to optimize the PS generation:
      call multi_kaleu_collect(iproc,weight)
!
! Update statistics
      s0(iproc) = s0(iproc) + 1d0
      s1(iproc) = s1(iproc) + weight
      s2(iproc) = s2(iproc) + weight*weight
!  
! Print intermediate results
      if (mod(iev,100000).eq.1.and.iproc.eq.1) write(6,'(a8,a16,a12,a8,2a10,a9)') &
        'iproc','average','sigma','in %','accepted','generatd','iev/nev'
      if (mod(iev,10000).eq.0) then
        call statistics( ave(iproc),sig(iproc),rel(iproc) ,s0(iproc),s1(iproc),s2(iproc) )
        write(6,'(i3, e16.8, e12.4,    f8.4, 2f10.0,         f9.4 )') &
                  iproc, ave(iproc),   sig(iproc), 100*rel(iproc),  se(iproc),s0(iproc), dble(iev)/nev
        if (iproc.eq.nproc) write(6,*) "******************************"
      endif
!
! Optimization phase
      if (optim.and.int(se(iproc)).eq.noptim) then
        optim = .false.
        call statistics( ave(iproc),sig(iproc),rel(iproc) ,s0(iproc),s1(iproc),s2(iproc) )
        write(6,*)
        write(6,'(e16.8,e12.4,f8.4,2f10.0)') ave(iproc),sig(iproc),100*rel(iproc),se(iproc),s0(iproc)
        write(6,*)
        write(6,*) 'Throwing away estimates, ' &
                  ,'keeping optimized distributions, '
        write(6,*) 'and re-starting MC.'
        write(6,*)
        s0opt(iproc) = s0(iproc)
        seopt(iproc) = se(iproc)
        s0(iproc) = 0d0
        s1(iproc) = 0d0
        s2(iproc) = 0d0
        se(iproc) = 0d0
      endif
!  
    end do
  end do
!
contains
!
!
  subroutine statistics( ave,sig,rel ,s0,s1,s2 )
  implicit none
  real(kind(1d0)) ,intent(out) :: ave,sig,rel
  real(kind(1d0)) ,intent(in)  :: s0,s1,s2
  ave = 0d0
  sig = 0d0
  rel = 0d0
  if (s0.gt.0d0) then
    ave = s1/s0
    sig = abs(ave)
    if (s0.gt.1d0) then
      sig = ( s2/s0 - ave*ave )/( s0-1d0 )
      if (sig.gt.0d0) sig = sqrt(sig)
    endif
    if (ave.ne.0d0) rel = sig/abs(ave)
  endif
  end subroutine
!
!
end program
