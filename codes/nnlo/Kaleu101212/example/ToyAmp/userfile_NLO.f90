module user_toyamp
  use avh_toyamp_model
  use avh_toyamp_dipol
  type(model_type) ,save :: mdl ! Masses and widths
  type(dipol_type) ,save :: obj ! Instance of the amplitude calculator
end module

  subroutine user_toyamp_init( process ,nfinst      &
                              ,masses,widths,labels &
                              ,onlyqcd,withhiggs    &
                              ,aff,afi,aif,aii      &
                              ,dip,ndip             )
!*********************************************************************
!* Input:
!*            process: the particles involved in the process
!*             nfinst: the number of final-state particles, ie "n+1"
!*             masses: particle masses defined in model
!*             widths: particle widths defined in model
!*             labels: particle labels defined in model
!*            onlyqcd: is only QCD involved?
!*          withhiggs: is there a higgs?
!*    aff,afi,aif,aii: cut-offs for the dipole phase space
!* Output:    
!*                dip: list of i,j,k of dipole terms
!*               ndip: number of dipole terms
!*********************************************************************
  use particles ! This is not part of ToyAmp, should come from outside
  use user_toyamp
  implicit none
  integer         ,intent(in) :: process(-2:17),nfinst
  real(kind(1d0)) ,intent(in) :: masses(nfirst:nlast),widths(nfirst:nlast)
  character(20)   ,intent(in) :: labels(nfirst:nlast)
  logical         ,intent(in) :: onlyqcd,withhiggs
  real(kind(1d0)) ,intent(in) :: aff,afi,aif,aii
  integer         ,intent(out) :: ndip
  integer         ,intent(out) :: dip(3,maxndip) ! maxndip defined in avh_toyamp_dipol
  logical :: withqcd,cancel
! The list vtx of possible interaction vertices is temporary. It is
! only used to create the tree stored inside obj
  type(vertx_type) :: vtx
!
  withqcd = .true.
  call simple_model( mdl,vtx ,masses,widths,labels &
                             ,onlyqcd,withqcd,withhiggs )
  call dipol_init( mdl,vtx,obj ,process,nfinst ,cancel &
                  ,aff,afi,aif,aii ,dip,ndip )
!
  end

  subroutine user_toyamp_mom( idip ,pkaleu ,pdip ,notzero )
!*********************************************************************
!* Input:
!*     idip: dipole term number. idip=0 refers to (n+1) real-rad.
!*   pkaleu: (n+1)-momentum configuration.
!* Output:
!*     pdip: mapped n-momentum configuration if dip>0,
!*           else it is just a copy of pkaleu.
!*  notzero: is false if the dipole term gives no contribution
!*           due to the cut-offs for the dipole phase space.
!*********************************************************************
  use user_toyamp
  implicit none
  real(kind(1d0))  ,intent(in)  :: pkaleu(0:3,-2:17)
  integer          ,intent(in)  :: idip
  real(kind(1d0))  ,intent(out) :: pdip(0:3,-2:17)
  logical          ,intent(out) :: notzero
  call dipol_mom( obj ,pkaleu ,idip ,pdip ,notzero )
  end

  subroutine user_toyamp_calc( idip ,pkaleu ,weight )
!*********************************************************************
!* Input:
!*     idip: dipole term number. idip=0 refers to (n+1) case.
!*   pkaleu: (n+1)-momentum configuration.
!* Output:
!*   weight: 
!*********************************************************************
  use user_toyamp
  implicit none
  integer          ,intent(in)  :: idip
  real(kind(1d0))  ,intent(in)  :: pkaleu(0:3,-2:17)
  real(kind(1d0))  ,intent(out) :: weight
  call dipol_calc( mdl,obj ,idip ,pkaleu ,weight )
  end

  subroutine user_toyamp_close
!*********************************************************************
!* If you do not need ToyAmp anymore while your program continues,
!* you may call this routine to de-allocate anything related to ToyAmp
!*********************************************************************
  use user_toyamp
  implicit none
  call dipol_close( obj )
  end
