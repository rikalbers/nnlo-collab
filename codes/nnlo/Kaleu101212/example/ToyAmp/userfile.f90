module user_toyamp
  use avh_toyamp
  use avh_toyamp_model
  type(model_type)  ,save :: mdl ! Masses and widths
  type(toyamp_type) ,save :: obj ! Instance of the amplitude calculator
end module
      

  subroutine user_toyamp_init( process ,nfinst &
                              ,masses,widths,labels    &
                              ,onlyqcd,withqcd,withhiggs )
!*********************************************************************
!* Input:
!*   process: the particles involved in the process
!*    nfinst: the number of final-state particles
!*    masses: particle masses defined in model
!*    widths: particle widths defined in model
!*    labels: particle labels defined in model
!*   onlyqcd: is only QCD involved?
!*   withqcd: is      QCD involved?
!* withhiggs: is there a higgs?    
!*********************************************************************
  use particles ! This is not part of ToyAmp, should come from outside
  use user_toyamp
  implicit none
  integer          ,intent(in) :: process(-2:17),nfinst
!
  real(kind(1d0))  ,intent(in) :: masses(nfirst:nlast),widths(nfirst:nlast)
  character(20)    ,intent(in) :: labels(nfirst:nlast)
  logical          ,intent(in) :: onlyqcd,withqcd,withhiggs
  logical :: cancel
! The list vtx of possible interaction vertices is temporary. It is
! only used to create the tree stored inside obj
  type(vertx_type) :: vtx
!
  call simple_model( mdl,vtx ,masses,widths,labels &
                             ,onlyqcd,withqcd,withhiggs )
  call toyamp_put_process( mdl,vtx,obj ,process ,nfinst ,cancel=cancel )
  if (cancel) then
    write(*,*) 'ERROR in ToyAmp: impossible process'
    stop
  endif
!
  end

  subroutine user_toyamp_calc( amp ,pkaleu )
!*********************************************************************
!* Input:
!*   pkaleu: momenta
!* Output:
!*   amp: value of the amplitude
!*********************************************************************
  use user_toyamp
  implicit none
  complex(kind(1d0)) ,intent(out) :: amp 
  real(kind(1d0))    ,intent(in)  :: pkaleu(0:3,-2:17)
  call toyamp_put_mom( obj ,pkaleu )
  call toyamp_calc( mdl,obj ,amp )
  end

  subroutine user_toyamp_close
!*********************************************************************
!* If you do not need ToyAmp anymore while your program continues,
!* you may call this routine to de-allocate anything related to ToyAmp
!*********************************************************************
  use user_toyamp
  implicit none
  call toyamp_close( obj )
  end
