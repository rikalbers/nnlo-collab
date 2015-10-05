module user_kaleu
  use avh_kaleu
  use avh_kaleu_dipol
  use avh_kaleu_model
  use avh_kaleu_strf
  type(model_type) ,save :: mdl ! Masses and widths
  type(dipol_type) ,save :: obj ! Instance of the phase space generator
  type(strf_type)  ,save :: str ! Concerns the generation of x1,x2
  integer          ,save :: iinst ! How to deal with inst-momenta ?
end module
      

  subroutine user_kaleu_init( process ,nfinst            &
                             ,ecm ,option                &
                             ,masses,widths,labels       &
                             ,onlyqcd,withqcd,withhiggs  &
                             ,nbatch ,nstep ,thrs ,smin  &
                             ,dip,ndip                   &
                             ,aff,afi,aif,aii            )
!*********************************************************************
!* Everything is input.
!*  process: the particles involved in the process
!*   nfinst: the number of final-state particles
!*      ecm: the centre-of-mass energy
!*   option: how should Kaleu deal with the inst-momenta?
!*           option=0: inst-momenta are the same for each event
!*           option=1: Kaleu generates x1,x2 and constructs inst-momenta
!*           option=2: Kaleu takes inst-momenta as input
!*   masses: masses of particles in the model
!*   widths: widths of particles in the model
!*   labels: labels of particles in the model
!*   onlyqcd: are only QCD-particles involved?
!*   withqcd: should there be QCD-particles?
!* withhiggs: should there be a Higgs?
!*   nbatch: number of accepted ps-points before optimization step
!*    nstep: number of optimization steps
!*     thrs: threshold for channel-removal (after full optimization)
!*     smin: matrix of lower limits on positive 2-particle invariants,
!*           and upper limits on negative 2-particle invariants
!*      dip: list of dipoles
!*     ndip: the number of dipoles
!*     aff,afi,aif,aii: cut-offs for dipole phase space
!*********************************************************************
  use particles ! This is not part of Kaleu, defines "nfirst,nlast".
  use user_kaleu
  implicit none
  integer         ,intent(in) :: process(-2:17),nfinst
  real(kind(1d0)) ,intent(in) :: ecm
  integer         ,intent(in) :: option
  real(kind(1d0)) ,intent(in) :: masses(nfirst:nlast),widths(nfirst:nlast)
  character(20)   ,intent(in) :: labels(nfirst:nlast)
  logical         ,intent(in) :: onlyqcd,withqcd,withhiggs
  integer         ,intent(in) :: nbatch,nstep
  real(kind(1d0)) ,intent(in) :: thrs,smin(-2:17,-2:17)
  integer         ,intent(in) :: ndip
  integer         ,intent(in) :: dip(3,ndip)
  real(kind(1d0)) ,intent(in) :: aff,afi,aif,aii
  integer :: prcss(-2:17),noffset,nbatch_grid,nbatch_cut
  logical :: cancel
! The list vtx of possible interaction vertices is temporary. It is
! only used to create the tree stored inside obj.
  type(vertx_type) :: vtx
!
  call kaleu_hello
!
  call user_model( mdl,vtx ,masses,widths,labels       &
                           ,onlyqcd,withqcd,withhiggs  )
! Kaleu doesn't distinguish anti-particles from particles.
  prcss(-2:nfinst) = abs(process(-2:nfinst))
! Some more optimization parameters.
  noffset     = nbatch*nstep
  nbatch_grid = 1000
  nbatch_cut  = 1000
! Kaleu_dipol needs the integer value of the gluon label in the model,
! so it is input for the following routine.
  call dipol_init( mdl,vtx,obj ,prcss,nfinst ,ecm &
                  ,smin                           &
                  ,nbatch,nstep,thrs              &
                  ,noffset,nbatch_grid,nbatch_cut &
                  ,dip,ndip,gluon                 & ! gluon here
                  ,aff,afi,aif,aii                )
  iinst = option
  if (iinst.eq.1) call dipol_init_strf( str,obj )
  end


  subroutine user_kaleu_gnrt( discard ,x1kaleu,x2kaleu ,pkaleu )
!*********************************************************************
!* Everything is output, pkaleu is also input if you put option=2 at
!*   initialization. In the latter case, pkaleu(0:3,-2:-1) are taken
!*   as inst-momenta. It doesn't matter wether they are positive or
!*   negative energy on input, but realize that they are always
!*   negative-energy on output.
!* Event should get  weight=0  if discard=.true.
!* You will use x1kaleu x2kaleu only if you put option=1 at
!*   initialization, 
!*********************************************************************
  use user_kaleu
  implicit none
  logical          ,intent(out)   :: discard
  real(kind(1d0))  ,intent(out)   :: x1kaleu,x2kaleu
  real(kind(1d0))  ,intent(inout) :: pkaleu(0:3,-2:17)
  x1kaleu = 1d0
  x2kaleu = 1d0
  if     (iinst.eq.1) then
! Helac-convention: momenta are in cm-frame with shat=x1*x2*s
    call dipol_gnrt_strf( str,obj ,discard ,x1kaleu,x2kaleu )
! Alpgen-convention: momenta are in lab-frame,
!   with  p(0,-1)=-x1*sqrt(s)/2  and  p(0,-2)=-x2*sqrt(s)/2
!    NOT AVAILABLE
  elseif (iinst.eq.2) then
    call dipol_inst( obj ,pkaleu )
  endif
  call dipol_gnrt( mdl,obj ,discard ,pkaleu )
  end


  subroutine user_kaleu_wght( weight )
!*********************************************************************
!* Everything is output.
!* If you put option=2 at initialization, don't forget to put the
!* weight related to your generation of the inst-momenta yourself in
!* your Monte Carlo.
!*********************************************************************
  use user_kaleu
  implicit none
  real(kind(1d0)) ,intent(out) :: weight
  real(kind(1d0)) :: ww
  call dipol_wght( mdl,obj ,weight )
  if (iinst.eq.1) then
    call strf_wght( str ,ww )
    weight = weight*ww
  endif
  end


  subroutine user_kaleu_collect( weight )
!*********************************************************************
!* Everything is input.
!*********************************************************************
  use user_kaleu
  implicit none
  real(kind(1d0)) ,intent(in) :: weight
  if (iinst.eq.1.and.weight.gt.0d0) call strf_collect( str ,weight )
  call dipol_collect( obj ,weight )
  end


  subroutine user_kaleu_close
!*********************************************************************
!* If you do not need Kaleu anymore while your program continues,
!* you may call this routine to de-allocate anything related to Kaleu
!*********************************************************************
  use user_kaleu
  implicit none
  call dipol_close( obj )
  if (iinst.eq.1) call strf_close( str )
  end


  subroutine user_kaleu_plotgrids( iunit )
!*********************************************************************
!* Will pollute your directory with lots of files containing the
!* adapted distributions for the generation of the splitting
!* variables. These can be visualized with gnuplot.
!*********************************************************************
  use user_kaleu
  implicit none
  integer ,intent(in) :: iunit
  call dipol_plotgrids( obj ,iunit )
  if (iinst.eq.1) call strf_plot( str ,iunit )
  end 
