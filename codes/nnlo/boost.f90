! This routine performs a boost with beta in the direction of vec:
subroutine boost4vec(beta,vec,vin,vout)
use momenta
implicit none
!
  real(kind(1d0)) , intent(in) :: beta
  real(kind(1d0)) , dimension(3) , intent(in) :: vec
  type(mom) , intent(in) :: vin
  type(mom) , intent(out) :: vout
!
  real(kind(1d0)) :: gamma
  real(kind(1d0)) :: vb
!
!
  gamma   = 1d0/sqrt(1d0 - beta**2)
  vb      = vin%px*vec(1) + vin%py*vec(2) + vin%pz*vec(3)
  vout%px = vin%px + vec(1)*((gamma - 1d0)*vb + gamma*beta*vin%E)
  vout%py = vin%py + vec(2)*((gamma - 1d0)*vb + gamma*beta*vin%E)
  vout%pz = vin%pz + vec(3)*((gamma - 1d0)*vb + gamma*beta*vin%E)
  vout%E  = gamma*(vin%E + vb*beta)
!
end subroutine boost4vec
!
! This routine rotates a three-vector if the sine and cosine is given with a
! desired direction:
subroutine rotate3vec(dir,sinphi,cosphi,vec)
implicit none
!
  real(kind(1d0)) , dimension(3) , intent(in) :: dir
  real(kind(1d0)) , intent(in) :: sinphi,cosphi
  real(kind(1d0)) , dimension(3) , intent(inout) :: vec
!
  real(kind(1d0)) :: dirdotvec
  real(kind(1d0)) , dimension(3) :: dircrossvec
!
!
  dircrossvec(1)=dir(2)*vec(3)-dir(3)*vec(2)
  dircrossvec(2)=dir(3)*vec(1)-dir(1)*vec(3)
  dircrossvec(3)=dir(1)*vec(2)-dir(2)*vec(1)
!
  dirdotvec=dir(1)*vec(1)+dir(2)*vec(2)+dir(3)*vec(3)
!
  vec(:) = vec(:)+sinphi*dircrossvec(:) &
         - (1d0-cosphi)*(vec(:)-dir(:)*dirdotvec)
!
end subroutine rotate3vec
!
! This routine rotates the spatial part of a four-vector if the 
! sine and cosine is given with a desired direction:
subroutine rotate4vec(dir,sinphi,cosphi,p)
use momenta
implicit none
!
  real(kind(1d0)) , dimension(3) , intent(in) :: dir
  real(kind(1d0)) , intent(in) :: sinphi,cosphi
  type(mom), intent(inout) :: p
!
  real(kind(1d0)) :: dirdotvec
  real(kind(1d0)) , dimension(3) :: dircrossvec,vec
!
!
  vec(1) = p%px
  vec(2) = p%py
  vec(3) = p%pz
!
  dircrossvec(1)=dir(2)*vec(3)-dir(3)*vec(2)
  dircrossvec(2)=dir(3)*vec(1)-dir(1)*vec(3)
  dircrossvec(3)=dir(1)*vec(2)-dir(2)*vec(1)
!
  dirdotvec=dir(1)*vec(1)+dir(2)*vec(2)+dir(3)*vec(3)
!
  vec(:) = vec(:)+sinphi*dircrossvec(:) &
         - (1d0-cosphi)*(vec(:)-dir(:)*dirdotvec)
!
  p%px = vec(1)
  p%py = vec(2)
  p%pz = vec(3)
!
end subroutine rotate4vec
!
subroutine Lambda4vec(K,Khat,p,q)
use momenta
implicit none
!
  type(mom) , intent(in) :: K,Khat,p
  type(mom) , intent(out) :: q
!
  real(kind(1d0)) :: KKh2,K2,KKhp,Khp
  type(mom) :: KKh
!
!
  KKh  = K + Khat
  KKh2 = KKh*KKh
  K2   = K*K
  KKhp = KKh*p
  Khp  = Khat*p
!
  q = p - 2d0/KKh2*KKhp*KKh + 2d0/K2*Khp*K
!
end subroutine Lambda4vec
