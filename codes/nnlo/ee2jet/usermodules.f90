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
  complex(kind(1d0)) , dimension(2,2,2) :: Couplings
  complex(kind(1d0)) , dimension(2) :: AxialCouplings
!
contains
!
subroutine init_mymodel()
use math
use collider
use input
implicit none
!
!
!
  real(kind(1d0)) :: qup,qdown
  real(kind(1d0)) :: sw,cw,scw
  real(kind(1d0)) :: veL,veR
  real(kind(1d0)) , dimension(2) :: vqL,vqR
  complex(kind(1d0)) :: Pz
!
!
  phys_Zmass = 91.1876d0
  phys_Wmass = 80.385d0
!
  phys_Zwidth = 2.4952d0
  phys_Wwidth = 2.085d0
!
  phys_Tmass = 172.5d0
  phys_Hmass = 125d0
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
!
  write(*,*) "********************************************************"
  write(*,*) "*                                                      *"
  write(*,*) "*           Physical parameters used in                *"
  write(*,*) "*                this calculation:                     *"
  write(*,*) "*                mZ     = ",phys_Zmass
  write(*,*) "*                mW     = ",phys_Wmass
  write(*,*) "*                GammaZ = ",phys_Zwidth
  write(*,*) "*                GammaW = ",phys_Wwidth
  write(*,*) "*                                                      *"
  write(*,*) "********************************************************"
!
! We have to set up the couplings:
  qup   = 2d0/3d0
  qdown = -1d0/3d0
!
  sw  = sqrt(1d0 - phys_Wmass**2/phys_Zmass**2)
  cw  = phys_Wmass/phys_Zmass
  scw = sw*cw
!
  veL = (sw**2 - 0.5d0)/scw
  veR = sw/cw
!
  vqL(1) = ( 0.5d0 - sw**2*qup)/scw
  vqL(2) = (-0.5d0 - sw**2*qdown)/scw
!
  vqR(1) = -sw/cw*qup
  vqR(2) = -sw/cw*qdown
!
  Pz = stot/(stot - phys_Zmass**2 &
     + cmplx(0d0,1d0)*phys_Zwidth*phys_Zmass)
!
  couplings(1,1,1) = - qup   + 0*veL*vqL(1)*Pz
  couplings(1,1,2) = - qdown + 0*veL*vqL(2)*Pz
  couplings(1,2,1) = - qup   + 0*veL*vqR(1)*Pz
  couplings(1,2,2) = - qdown + 0*veL*vqR(2)*Pz
  couplings(2,1,1) = - qup   + 0*veR*vqL(1)*Pz
  couplings(2,1,2) = - qdown + 0*veR*vqL(2)*Pz
  couplings(2,2,1) = - qup   + 0*veR*vqR(1)*Pz
  couplings(2,2,2) = - qdown + 0*veR*vqR(2)*Pz
!
  axialcouplings(1) = 0*veL*Pz/scw
  axialcouplings(2) = 0*veR*Pz/scw
!
  print *,"sw,cw: ",sw,cw
  print *,pi/nnloinput("alphaem")/sqrt(2d0)/phys_Wmass**2/sw**2
  print *,"stot: ",sqrt(stot)
!
end subroutine init_mymodel
!
end module my_model
