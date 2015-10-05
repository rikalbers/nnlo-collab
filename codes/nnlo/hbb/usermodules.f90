module my_scale
implicit none
!
!
contains
!
subroutine calcmyscales(p)
use scales
use particles
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
!
  if (muscale.lt.0) then
    print *,"No dynamic scale is defined..."
    stop
  end if
!
  mur = xir*muscale
  muf = xif*muscale
!
end subroutine calcmyscales
!
end module my_scale
!
module my_model
implicit none
!
  real(kind(1d0)) :: phys_Zmass
  real(kind(1d0)) :: phys_Wmass
  real(kind(1d0)) :: phys_Hmass
  real(kind(1d0)) :: phys_Tmass
!
  real(kind(1d0)) :: phys_Zwidth
  real(kind(1d0)) :: phys_Wwidth
  real(kind(1d0)) :: phys_Hwidth
  real(kind(1d0)) :: phys_Twidth
!
  real(kind(1d0)) :: phys_yb
!
contains
!
subroutine init_mymodel()
use math
use collider
use input
use QCDparams
use phasespace
implicit none
!
!
!
!
!
  phys_Zmass = 91.1876d0
  phys_Wmass = 80.385d0
!
  phys_Zwidth = 2.4952d0
  phys_Wwidth = 2.085d0
!
  phys_Tmass = 172.5d0
  phys_Hmass = 120d0
!
  phys_Twidth = 0d0
  phys_Hwidth = 0d0
!
  if (nnloinput("#zmass").gt.0d0) then
    phys_Zmass = nnloinput("#zmass")
  end if
  if (nnloinput("#zwidth").gt.0d0) then
    phys_Zwidth = nnloinput("#zwidth")
  end if
  if (nnloinput("#wmass").gt.0d0) then
    phys_Wmass = nnloinput("#wmass")
  end if
  if (nnloinput("#wwidth").gt.0d0) then
    phys_Wwidth = nnloinput("#wwidth")
  end if
  if (nnloinput("#hmass").gt.0d0) then
    phys_Wmass = nnloinput("#hmass")
  end if
  if (nnloinput("#hwidth").gt.0d0) then
    phys_Wwidth = nnloinput("#hwidth")
  end if
!
  write(*,*) "********************************************************"
  write(*,*) "*                                                      *"
  write(*,*) "*           Physical parameteres used in               *"
  write(*,*) "*                this calculation:                     *"
  write(*,*) "*                mZ     = ",phys_Zmass
  write(*,*) "*                mW     = ",phys_Wmass
  write(*,*) "*                mH     = ",phys_Hmass
  write(*,*) "*                GammaZ = ",phys_Zwidth
  write(*,*) "*                GammaW = ",phys_Wwidth
  write(*,*) "*                GammaH = ",phys_Hwidth
  write(*,*) "*                                                      *"
  write(*,*) "********************************************************"
!
! We have to set up the couplings:
! We set up the Yukawa coupling such a way that the LO
! result is one:
  phys_yb = sqrt(4d0*pi/phys_Hmass/qcd_nc)/sqrt(0.5d0*gevm22pb/phys_Hmass)
!
end subroutine init_mymodel
!
end module my_model
