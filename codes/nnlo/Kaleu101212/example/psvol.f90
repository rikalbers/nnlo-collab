program montecarlo
  use particles
  use cuts_module
  use rambo_module
  implicit none
  type(rambo_type) :: rambo
  type(cuts_type) :: cuts
  integer         :: process(-2:17),nfinst
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
  real(kind(1d0)) :: ave,sig,rel ,s0,s1,s2,se ,s0opt,seopt
  complex(kind(1d0)) :: amp
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
! Read process. Numbering defined in particles.f90
  read(5,*) nfinst
  read(5,*) itmp(1:nfinst+2)
  process(-2:-1    ) = itmp(1:2)
  process( 1:nfinst) = itmp(3:nfinst+2)
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
  allocate( labels(nfirst:nlast) ) ! "nfirst" and "nlast" are defined in
  allocate( masses(nfirst:nlast) ) ! the module "particles"
  allocate( widths(nfirst:nlast) ) !
  include 'constants.h90'
!
!! Initialize amplitude calculator
!  call user_toyamp_init( process ,nfinst           &
!                        ,masses,widths,labels      &
!                        ,onlyqcd,withqcd,withhiggs )
!
! Initialize cuts, and get smin to be delivered to Kaleu
  call cuts_init( cuts ,process,nfinst ,masses ,smin )
!
! No pdfs in this example
  option = 0
!
! Initialize phase space generator
  if (igenerator.eq.0) then
    call rambo_init( rambo ,process,nfinst ,ecm,option ,masses )
  else
    call user_kaleu_init( process ,nfinst           & ! Routines starting with
                         ,ecm ,option               & ! "user_kaleu" are defined
                         ,masses,widths,labels      & ! in the file
                         ,onlyqcd,withqcd,withhiggs & ! "kaleu_userfile.f90"
                         ,nbatch ,nstep ,thrs ,smin ) !
  endif
!
! Number of events to be spent in optimization phase
  if (igenerator.eq.0) then
    noptim = 0
    optim = .false.
  else
    noptim = nbatch*nstep
    optim  = .true.
  endif
!
! Initialize Monte Carlo
  s0 = 0d0 ! number of generated events
  s1 = 0d0 ! sum of full weights
  s2 = 0d0 ! sum of squares of full weights
  se = 0d0 ! number of accepted events
!
! Begin Monte Carlo loop
  do iev=1,nev
!
! Generate event
    if (igenerator.eq.0) then
      call rambo_gnrt( rambo ,pkaleu )
      discard = .false.
    else
      call user_kaleu_gnrt( discard ,x1,x2 ,pkaleu )
    endif
!
    if (discard) then
      weight = 0d0
!  
    elseif (cuts_pass(cuts,pkaleu)) then
! Calculate weight from generator
      if (igenerator.eq.0) then
        call rambo_wght( rambo ,weight )
      else
        call user_kaleu_wght( weight )
      endif
! Evaluate amplitude
!      call user_toyamp_calc( amp ,pkaleu )
! Full weight
!      weight = weight*dreal( amp*conjg(amp) )
! We have to include a tower of 2*pis:
      weight = weight * (2d0*pi)**(4 - 3*nfinst)
      se = se + 1d0
!  
    else
      weight = 0d0
!  
    endif
!
! Give full weight back to Kaleu
    if (igenerator.eq.0) then
    else
      call user_kaleu_collect( weight )
    endif
!
! Update statistics
    s0 = s0 + 1d0
    s1 = s1 + weight
    s2 = s2 + weight*weight
!  
! Print intermediate results
    if (mod(iev,100000).eq.1) write(6,'(a16,a12,a8,2a10,a9)') &
      'average','sigma','in %','accepted','generatd','iev/nev'
    if (mod(iev,10000).eq.0) then
      call statistics( ave,sig,rel ,s0,s1,s2 )
      write(6,'(e16.8, e12.4,    f8.4, 2f10.0,         f9.4 )') &
                  ave,   sig, 100*rel,  se,s0, dble(iev)/nev
    endif
!
! Optimization phase
    if (optim.and.int(se).eq.noptim) then
      optim = .false.
      call statistics( ave,sig,rel ,s0,s1,s2 )
      write(6,*)
      write(6,'(e16.8,e12.4,f8.4,2f10.0)') ave,sig,100*rel,se,s0
      write(6,*)
      write(6,*) 'Throwing away estimates, ' &
                ,'keeping optimized distributions, '
      write(6,*) 'and re-starting MC.'
      write(6,*)
      s0opt = s0
      seopt = se
      s0 = 0d0
      s1 = 0d0
      s2 = 0d0
      se = 0d0
    endif
!  
  enddo
! End Monte Carlo loop
!
! Print results
  call statistics( ave,sig,rel ,s0,s1,s2 )
  write(6,*) 'result:',ave,' +/-',sig
  write(6,*) 'relative error in %:',100*rel
  write(6,*) 'number of evaluations:',int(se),', for optimization:',int(seopt)
  write(6,*) 'number of generations:',int(s0),', for optimization:',int(s0opt)
!
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
