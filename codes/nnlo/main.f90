program nnlo
implicit none
!
!
!
  call NNLOinit()
!
  call NNLOcalc()
!
  stop
!
end program nnlo
!
subroutine NNLOinit()
use process
use collider
use scales
use flags
use QCDparams
use coupling
use my_model
use random
use alphamax
use subprocesses
use phasespace
implicit none
!
!
call init_process()
!
call init_collider()
!
call init_flags()
!
call init_seed()
!
call init_scales()
!
call init_QCDparams()
!
call init_couplings()
!
call init_mymodel()
!
call init_alphamax()
!
call init_contribs()
!
call init_subprocesses()
!
call init_PS()
!
call init_Regions()
!
end subroutine NNLOinit
