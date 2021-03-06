c      program sqmatbeta
c      implicit none
c
c      integer i, ps, con
c      real*8 MMfull(-4:0), y13, y23, mu, Q, xi
c      real*8 MM0(0:2), MM1(-2:2), MM2(-4:0)
c
c      write(6,*) 'Input kinematic variable values (y13,y23)!'
c      read(5,*) y13, y23
c      write(6,*) 'Input renormalization scale and total energy! ',
c     $ '(not squared)'
c      read(5,*) mu, Q
c      xi=(mu/Q)**2
c      write(6,*) 'Choose convention!'
c      write(6,*) ' 0 - GR'
c      write(6,*) ' 1 - ST'
c      read(5,*) con
c      write(6,*) 'Do you want to calculate the poles?'
c      write(6,*) ' 0 - no'
c      write(6,*) ' 1 - yes'
c      read(5,*) ps
c
c      call amp2(MM0,MM1,MM2,y13,y23,xi,ps,con)
c
c      write(6,*) 'Born'
c      do i=0,2
c       write(6,*) i, MM0(i)
c      end do
c
c      write(6,*) 'One loop squared matrix element'
c      do i=-2,2
c       write(6,*) i, MM1(i)
c      end do
c
c      write(6,*) 'Two loop squared matrix element'
c      do i=-4,0
c       write(6,*) i, MM2(i)
c      end do
c
c      end


************************************************************************
****  doubly virtual squared matrix element for gamma -> q g qbar   ****
****            based on the work of T. Gehrmann et al.             ****
****              hep-ph/0112081v1, hep-ph/0710.0346v3               ****
! UV renormalized
!   yij - kinematic variable, 1:quark, 2:antiquark, 3:gluon
!   MM0(0:2) - Born
!   MM1(-2:2) - one loop squared matrix element
!   MM2full(-4:0) - two loop squared matrix element
!   xi - mu^2/Q^2,  mu: renormaization scale, Q: scattering energy
!   ps - pole switch, 0: no pole calculation, 1: include pole calculation
!   con - convention switch, 1: use ST convention, 0: use GR convention
************************************************************************
      subroutine amp2(MM0,MM1,MM2full,y13,y23,xi,ps,con)
      implicit none
      real*8 NC, nf, TR, nfg, pi, z3, eg
      common /piez/ pi, z3
      common /const/ NC, nf, TR, eg

      real*8 MM2full(-4:0), MM0(0:2), MM1(-2:2)
      real*8 M0M0(-4:2), M0M1(-4:2), MM2(-4:0)
      real*8 fin20, fin11, lof, pol2(-4:0), pol1(-4:0)
      real*8 I1qqg(-4:2), I1qqg2(-4:2), I1qqgsq(-4:2), H2qqgm1
      real*8 xi, y13, y23, y12
      real*8 beta0, beta1, K, CF, CA
      complex*16 f1yz(-4:2), f1zy(-4:2), f2yz(-4:2), f2zy(-4:2)
      integer i, j, ps, con

! checking kinematic variable range
      if(y13.le.0d0 .or. y13.ge.1d0) then
       write(6,*) 'Warning! y13 is outside of the ]0,1[ interval!'
      end if
      if(y23.le.0d0 .or. y23.ge.1d0) then
       write(6,*) 'Warning! y23 is outside of the ]0,1[ interval!'
      end if
! checking mu^2/Q^2
      if(xi.le.0d0) then
       write(6,*) 'Squared ratio of renormalization scale and ',
     $  'scattering energy must be nonzero positive!'
       stop
      end if
! checking ps
      if(ps.ne.1 .and. ps.ne.0) then
       write(6,*) 'ps must be 0 or 1!'
       stop
      end if
! checking con
      if(con.ne.1 .and. con.ne.0) then
       write(6,*) 'con must be 0 or 1!'
       stop
      end if


      NC=3d0
      nf=5d0
      TR=0.5d0
      if(nf.eq.6.0d0) then
       nfg=9/15d0
      else if(nf.eq.5.0d0) then
       nfg=1/11d0
      else if(nf.eq.4.0d0) then
       nfg=2/5d0
      else if(nf.eq.3.0d0) then
       nfg=0d0
      else if(nf.eq.2.0d0) then
       nfg=1/5d0
      else if(nf.eq.1.0d0) then
       nfg=1d0
      end if
      pi=3.14159265358979324d0
      z3=1.20205690315959429d0 !zeta(3)
      eg=0.57721566490153286d0 !eulergamma

      CF=TR*(NC**2-1)/NC
      CA=2*TR*NC
      beta0=(11*CA-4*TR*nf)/6
      beta1=(17*CA**2-10*CA*TR*nf-6*CF*TR*nf)/6
      K=(67.0d0/18-pi**2/6)*CA-10.0d0/9*TR*nf
      y12=1-y13-y23


! |M0|^2
      do i=-4,-1
       M0M0(i)=0.0d0
      end do
      M0M0(0)=4*(NC**2-1)*(y13/y23+y23/y13+2*y12/(y13*y23))
      M0M0(1)=-8*(NC**2-1)*(y13/y23+y23/y13+y12/(y13*y23)+1)
      M0M0(2)=4*(NC**2-1)*(y13/y23+y23/y13+2)


! 2Re<M0|M1>
      call f1(f1yz,y13,y23)
      call f1(f1zy,y23,y13)
      call f2(f2yz,y13,y23)
      call f2(f2zy,y23,y13)

      do i=-4,2
       M0M1(i)=2*(NC**2-1)*REAL(NC*(f1yz(i)+f1zy(i))
     $         +1/NC*(f2yz(i)+f2zy(i)))
       if(i.lt.2) M0M1(i)=M0M1(i)-beta0*M0M0(i+1)
      end do


      if(ps.eq.1) then

      I1qqg(-4)=0.0d0
      I1qqg(-3)=0.0d0
      I1qqg(-2)=-NC+1/(2*NC)

      I1qqg(-1)=NC/6d0*(-10+3*lof(y13,1)+3*lof(y23,1))-1/NC*(
     $ lof(y12,1)/2d0-3/4d0)+nf/6d0

      I1qqg(0)=NC/12d0*(pi**2+10*(lof(y13,1)+lof(y23,1))
     $ - 3*(lof(y13,2)+lof(y23,2)))-1/NC*(pi**2/24d0
     $ + 3*lof(y12,1)/4d0-lof(y12,2)/4d0)
     $ - nf/12d0*(lof(y13,1)+lof(y23,1))

      I1qqg(1)=NC/72d0*(10*pi**2+24*z3-3*pi**2*(lof(y13,1)+lof(y23,1))
     $ - 30*(lof(y13,2)+lof(y23,2))+6*(lof(y13,3)+lof(y23,3)) )
     $ - 1/(48*NC)*(3*pi**2+8*z3-2*pi**2*lof(y12,1)-18*lof(y12,2)
     $ + 4*lof(y12,3) )+nf/72d0*(3*(lof(y13,2)+lof(y23,2))-pi**2)

      I1qqg(2)=NC/1440d0*( -pi**4+800*z3+10*(
     $ - (10*pi**2+24*z3)*(lof(y13,1)+lof(y23,1))
     $ + 3*pi**2*(lof(y13,2)+lof(y23,2))+20*(lof(y13,3)+lof(y23,3))
     $ - 3*(lof(y13,4)+lof(y23,4))) )
     $ - 1/(2880*NC)*( -pi**4+720*z3-60*((3*pi**2+8*z3)*lof(y12,1)
     $ - pi**2*lof(y12,2)-6*lof(y12,3)+lof(y12,4)) )
     $ + nf/144d0*( -8*z3+pi**2*(lof(y13,1)+lof(y23,1))
     $ - 2*(lof(y13,3)+lof(y23,3)) )

      do i=-4,2
      I1qqg2(i)=I1qqg(i)*(2d0**i)
      end do

      I1qqgsq(-4)=I1qqg(-2)**2
      I1qqgsq(-3)=2*I1qqg(-2)*I1qqg(-1)
      I1qqgsq(-2)=2*I1qqg(-2)*I1qqg(0)+I1qqg(-1)**2
      I1qqgsq(-1)=2*(I1qqg(-2)*I1qqg(1)+I1qqg(-1)*I1qqg(0))
      I1qqgsq(0) =2*(I1qqg(-1)*I1qqg(1)+I1qqg(-2)*I1qqg(2))+I1qqg(0)**2


      H2qqgm1=( (4*z3+589/432d0-11*pi**2/72d0)*NC**2
     $ - z3/2d0-41/54d0-pi**2/48d0+(pi**2/4d0-3*z3-3/16d0)/NC**2
     $ + (pi**2/36d0-19/18d0)*NC*nf-(1/54d0+pi**2/24d0)*nf/NC
     $ + 5/27d0*nf**2 )/4d0


! poles for 2Re<M0|M2> + |M1|^2
! 2x0 part
      do i=-4,0
       pol2(i)=2*(-I1qqgsq(i)+beta0*(I1qqg2(i+1)-I1qqg(i+1))
     $         +K*I1qqg2(i))
       if(i.gt.-4) then
        pol2(i)=pol2(i)+2*(beta0*pi**2/4d0*I1qqg2(i-1))
       if(i.gt.-3) then
        pol2(i)=pol2(i)+2*(K*pi**2/4+7*beta0/3d0*z3)*I1qqg2(i-2)
       end if
       end if
      end do
      pol2(-1)=pol2(-1)+H2qqgm1*2

      do i=0,2
       pol2(-i)=pol2(-i)*M0M0(0)+pol2(-i-1)*M0M0(1)+pol2(-i-2)*M0M0(2)
      end do
      pol2(-3)=pol2(-3)*M0M0(0)+pol2(-4)*M0M0(1)
      pol2(-4)=pol2(-4)*M0M0(0)

! 1x1 part
      do i=-4,0
       pol1(i)=0.0d0
       do j=0,2
        pol1(i)=pol1(i)+2*I1qqg(-j)*M0M1(i+j)
       end do
       if(i.gt.-4) then
        pol1(i)=pol1(i)+2*I1qqg(1)*M0M1(i-1)
       if(i.gt.-3) then
        pol1(i)=pol1(i)+2*I1qqg(2)*M0M1(i-2)
       end if
       end if
      end do


      do i=-4,0
       MM2(i)=pol2(i)+pol1(i)
      end do

      else
      do i=-4,0
       MM2(i)=0d0
      end do
      end if


      MM2(0)=MM2(0) + (NC**2-1)*(
     $  NC**2*(fin20(1,y13,y23)+fin20(1,y23,y13))
     $ + fin20(2,y13,y23)+fin20(2,y23,y13)
     $ + (fin20(3,y13,y23)+fin20(3,y23,y13))/(NC**2)
     $ + NC*nf*(fin20(4,y13,y23)+fin20(4,y23,y13))
     $ + nf/NC*(fin20(5,y13,y23)+fin20(5,y23,y13))
     $ + nf**2*(fin20(6,y13,y23)+fin20(6,y23,y13))
     $ + nfg*(4/NC-NC)*(fin20(7,y13,y23)+fin20(7,y23,y13)) )
     $+(NC**2-1)*(
     $ NC**2*(fin11(1,y13,y23)+fin11(1,y23,y13))
     $ + fin11(2,y13,y23)+fin11(2,y23,y13)
     $ + (fin11(3,y13,y23)+fin11(3,y23,y13))/(NC**2)
     $ + NC*nf*(fin11(4,y13,y23)+fin11(4,y23,y13))
     $ + nf/NC*(fin11(5,y13,y23)+fin11(5,y23,y13))
     $ + nf**2*(fin11(6,y13,y23)+fin11(6,y23,y13)) )


! restoring renormalization scale dependency

      do i=-2,2
       MM1(i)=M0M1(i)+beta0*log(xi)*M0M0(i)
     $        -beta0/2d0*log(xi)**2*M0M0(i-1)
     $        +beta0/6d0*log(xi)**3*M0M0(i-2)
       if(i.ge.0) MM0(i)=M0M0(i)
      end do

      do i=-4,0
      MM2full(i)=MM2(i)
     $          +(2*beta0*M0M1(i) + beta1*M0M0(i))*log(xi)
     $          +M0M0(i)*(beta0*log(xi))**2
      if(i.gt.-4) then
       MM2full(i)=MM2full(i)
     $           -(beta0*M0M1(i-1)+beta1*M0M0(i-1))*log(xi)**2
     $           -beta0**2*M0M0(i-1)*log(xi)**3
      if(i.gt.-3) then
       MM2full(i)=MM2full(i)
     $           +((4*beta0*M0M1(i-2)+8*beta1*M0M0(i-2))*log(xi)**3 
     $            +7*beta0**2*M0M0(i-2)*log(xi)**4)/12d0
      end if
      end if
      end do


      if(con.eq.1) then
! multiplied by Sms/S
       MM1(2)=MM1(2)+MM1(0)*pi**2/12d0+MM1(-1)*z3/3d0
     $  +MM1(-2)*pi**4/160d0
       MM1(1)=MM1(1)+MM1(-1)*pi**2/12d0+MM1(-2)*z3/3d0
       MM1(0)=MM1(0)+MM1(-2)*pi**2/12d0
! multiplied by (Sms/S)^2
       MM2full(0)=MM2full(0)+MM2full(-2)*pi**2/6d0+MM2full(-3)*2*z3/3d0
     $  +MM2full(-4)*7*pi**4/360d0
       MM2full(-1)=MM2full(-1)+MM2full(-3)*pi**2/6d0
     $  +MM2full(-4)*2*z3/3d0
       MM2full(-2)=MM2full(-2)+MM2full(-4)*pi**2/6d0
      end if

      end subroutine


************************************************************************
****     functions f1 and f2 are necessary for the <M0|M1> term     ****
****                  coded from hep-ph/0112081v1                   ****
! y and z - kinematic variables y13 and y23 in arbitrary order
! i - epsilon exponent
************************************************************************
      subroutine f1(f,y,z)
      implicit none
      integer i
      real*8 y, z
      complex*16 f(-4:2)
      complex*16 Bub1(-8:3), Buby(-8:3), Bubz(-8:3), Box6(-8:2)

      call bub(Bub1,1d0)
      call bub(Buby,y)
      call bub(Bubz,z)
      call boxsix(Box6,y,z)

      do i=-4,2
      f(i)=(-3*Bub1(i)+Bub1(i-1)+2*Bub1(i-2)-4*Buby(i+1)
     $   +12*Buby(i)-8*Buby(i-1))/(y*z)
     $  + y/z*(-2*Bubz(i+1)+8*Bubz(i)-10*Bubz(i-1)+3*Bubz(i-2)
     $   +Bubz(i-3)-3*Bub1(i)+4*Bub1(i-1)+Bub1(i-2)
     $   -2*Bub1(i-3)-2*Buby(i+1)+8*Buby(i)-10*Buby(i-1)+4*Buby(i-2))
     $  + (4*Bubz(i+1)-12*Bubz(i)+9*Bubz(i-1)-Bubz(i-2)
     $   +6*Bub1(i)-2*Bub1(i-1)-4*Bub1(i-2)+4*Buby(i+1)
     $   -12*Buby(i)+8*Buby(i-1))/z
     $  + y/(1-z)**2*(Bubz(i)-Bubz(i-1)-Bub1(i)+Bub1(i-1))
     $  + y/(1-z)*(3*Bubz(i)-5*Bubz(i-1)+2*Bubz(i-3)-3*Bub1(i)
     $   +4*Bub1(i-1)+Bub1(i-2)-2*Bub1(i-3))
     $  + (4*Bub1(i)-3*Bub1(i-1)-3*Bub1(i-2)-2*Bub1(i-3)
     $   -4*Bubz(i)+3*Bubz(i-1)+3*Bubz(i-2)+2*Bubz(i-3))/(1-z)
     $  + 4*Bubz(i)-9*Bubz(i-1)+6*Bubz(i-2)-Bubz(i-3)
     $ + Box6(i)*(-8/z-2*y**2/z+6*y/z-2*z+4+2/(y*z))
     $ + Box6(i-1)*(24/z+8*y**2/z-20*y/z+6*z-11-6/(y*z))
     $ + Box6(i-2)*(-16/z-10*y**2/z+18*y/z-12*z+9+4/(y*z))
     $ + Box6(i-3)*(-6-4*y/z+4*y**2/z+16*z)
      end do

      end subroutine


      subroutine f2(f,y,z)
      implicit none
      integer i
      real*8 y, z
      complex*16 f(-4:2)
      complex*16 Bub1(-8:3), Bubx(-8:3), Bubz(-8:3), Box6(-8:2)

      call bub(Bub1,1d0)
      call bub(Bubx,1-y-z)
      call bub(Bubz,z)
      call boxsix(Box6,1-y-z,z)

      do i=-4,2
      f(i)=(3*Bub1(i)-Bub1(i-1)-2*Bub1(i-2)+2*Bubx(i+1)
     $  -6*Bubx(i)+4*Bubx(i-1))/(y*z) + y/z*(-Bubz(i-2)
     $  +Bubz(i-3)+3*Bub1(i)-4*Bub1(i-1)-Bub1(i-2)+2*Bub1(i-3)
     $  +2*Bubx(i+1)-8*Bubx(i)+10*Bubx(i-1)
     $  -4*Bubx(i-2)) + (Bubz(i-1)-Bubz(i-2)-6*Bub1(i)
     $  +2*Bub1(i-1)+4*Bub1(i-2)-4*Bubx(i+1)+12*Bubx(i)
     $  -8*Bubx(i-1))/z + 2/(y+z)**2*(Bubx(i)-Bub1(i))
     $  +2/(y+z)*(Bub1(i-1)-2*Bubx(i-1))+y/(1-z)**2*(Bub1(i)
     $   -Bub1(i-1)-Bubz(i)+Bubz(i-1)) + y/(1-z)*(3*Bub1(i)
     $  -4*Bub1(i-1)-Bub1(i-2)+2*Bub1(i-3)-3*Bubz(i)+5*Bubz(i-1)
     $  -2*Bubz(i-3))+(2*Bubz(i)+Bubz(i-1)-5*Bubz(i-2)-2*Bubz(i-3)
     $  -2*Bub1(i)-Bub1(i-1)+5*Bub1(i-2)+2*Bub1(i-3))/(1-z)
     $  +2*Bubz(i)-7*Bubz(i-1)+2*Bubz(i-2)+3*Bubz(i-3)
     $  -4*Bubx(i)+10*Bubx(i-1)-4*Bubx(i-2)
     $ +Box6(i)*(8-4*y-4/z+4*y/z-2*y**2/z-4*z)
     $ +Box6(i-1)*(12*y+12/z-12*y/z+8*y**2/z+12*z-20)
     $ +Box6(i-2)*(8-14*y-8/z+8*y/z-10*y**2/z-14*z)
     $ +Box6(i-3)*(10*y+4*y**2/z+10*z)
     $ +Box6(i-4)*4*(y+z)
      end do

      end subroutine


************************************************************************
****                    Bubble and Box6 itegrals                    ****
****                  coded from hep-ph/0112081v1                   ****
! Multiplication by S_e^-1 included. S_e=(4pi/Exp[EulerGamma])^e
! Corrections:
!  - a factor of i/(16Pi^2) had to be removed
!  - the coefficient index on l in Box6 had to be reversed
! Using T. Gehrmann and E. Remiddi's HPLs and 2dHPLs code0 from
! hep-ph/0107173 and hep-ph/0111255
************************************************************************
      subroutine Bub(buba,x)
      implicit none
      integer i
      real*8 x, pi, z3
      complex*16 buba(-8:3)
      common /piez/ pi, z3

      do i=-8,-2
       buba(i)=0.0d0
      end do

      buba(-1)=DCMPLX(1.0d0,0.0d0)

      buba(0)=DCMPLX((2-log(x)),-pi)

      buba(1)=DCMPLX(4-7*pi**2/12d0-2*log(x)+log(x)**2/2d0,
     $ -2*pi+pi*log(x))

      buba(2)=1/12d0*DCMPLX(96-14*pi**2-(48-7*pi**2)*log(x)+12*log(x)**2
     $ -2*log(x)**3-28*z3,-48*pi+3*pi**3+24*pi*log(x)-6*pi*log(x)**2)

      buba(3)=1/1440d0*DCMPLX(23040-3360*pi**2+73*pi**4-11520*log(x)
     $ +(1680*pi**2+3360*z3)*log(x)+(2880-420*pi**2)*log(x)**2
     $ -480*log(x)**3+60*log(x)**4-6720*z3,-11520*pi+720*pi**3
     $ +5760*pi*log(x)-360*pi**3*log(x)-1440*pi*log(x)**2
     $ +240*pi*log(x)**3+3360*pi*z3)

      end subroutine


      subroutine boxsix(Box6,y,z)
      implicit none
      integer i
      real*8 y, z, pi, z3, l(-2:2), GYZ1(0:3),GYZ2(0:3,0:3),
     $ GYZ3(0:3,0:3,0:3),GYZ4(0:3,0:3,0:3,0:3), Hr1(0:1),Hr2(0:1,0:1),
     $ Hr3(0:1,0:1,0:1),Hr4(0:1,0:1,0:1,0:1)
      complex*16 Bub1(-8:3), Buby(-8:3), Bubz(-8:3), c(0:4), Box6(-8:2)
      common /piez/ pi, z3

      call bub(Bub1,1d0)
      call bub(Buby,y)
      call bub(Bubz,z)

      call tdhpl(y,z,4,GYZ1,GYZ2,GYZ3,GYZ4,Hr1,Hr2,Hr3,Hr4)

      c(0)=-1/2d0*DCMPLX(1/(1-y-z),0.0d0)
      c(1)=-1/(2*(1-y-z))*DCMPLX(2.0d0,-pi)
      c(2)=-1/(24*(1-y-z))*DCMPLX(48-7*pi**2,-24*pi)
      c(3)=-1/(24*(1-y-z))*DCMPLX(96-14*pi**2-28*z3,3*pi**3-48*pi)
      c(4)=-1/(2880*(1-y-z))*DCMPLX(23040-3360*pi**2+73*pi**4
     $      -6720*z3,-11520*pi+720*pi**3+3360*pi*z3)


      l(2)=2*Hr1(0)*GYZ3(2,2,0)-2*Hr1(0)*GYZ3(2,0,0)
     $     -2*Hr1(0)*GYZ3(0,2,0)+2*Hr1(0)*GYZ3(0,0,0)
     $     -2*Hr2(0,0)*GYZ2(2,0)+2*Hr2(0,0)*GYZ2(0,0)
     $     +2*Hr3(0,0,0)*GYZ1(0)+2*Hr4(0,0,0,0)+2*Hr4(0,0,1,0)
     $     -2*Hr3(0,1,0)*GYZ1(2)+2*Hr3(0,1,0)*GYZ1(0)
     $     +2*Hr4(0,1,0,0)+2*Hr4(0,1,1,0)+2*Hr2(1,0)*GYZ2(2,2)
     $     -2*Hr2(1,0)*GYZ2(2,0)-2*Hr2(1,0)*GYZ2(0,2)
     $     +2*Hr2(1,0)*GYZ2(0,0)-2*Hr3(1,0,0)*GYZ1(2)
     $     +2*Hr3(1,0,0)*GYZ1(0)+2*Hr4(1,0,0,0)+2*Hr4(1,0,1,0)
     $     -2*Hr3(1,1,0)*GYZ1(2)+2*Hr3(1,1,0)*GYZ1(0)
     $     +2*Hr4(1,1,0,0)+2*Hr4(1,1,1,0)-2*GYZ4(2,2,1,0)
     $     +2*GYZ4(2,0,1,0)+2*GYZ4(2,1,0,0)+2*GYZ4(0,2,1,0)
     $     +2*GYZ4(0,0,0,0)-2*GYZ4(0,0,1,0)-2*GYZ4(0,1,0,0)
     $     -2*GYZ4(1,0,0,0)+7*pi**4/180 + pi**2/3*(-Hr1(0)*GYZ1(2)
     $     +Hr1(0)*GYZ1(0)+Hr2(0,0)+Hr2(0,1)-Hr1(1)*GYZ1(2)
     $     +Hr1(1)*GYZ1(0)+Hr2(1,0)+Hr2(1,1)+GYZ2(2,2)-GYZ2(2,0)
     $     -GYZ2(0,2)+GYZ2(0,0) )


      l(1)=2*Hr1(0)*GYZ2(2,0)-2*Hr1(0)*GYZ2(0,0)-2*Hr2(0,0)*GYZ1(0)
     $     -2*Hr3(0,0,0)-2*Hr3(0,1,0)+2*Hr2(1,0)*GYZ1(2)
     $     -2*Hr2(1,0)*GYZ1(0)-2*Hr3(1,0,0)-2*Hr3(1,1,0)-2*GYZ3(2,1,0)
     $     -2*GYZ3(0,0,0)+2*GYZ3(0,1,0)+2*GYZ3(1,0,0)
     $     +pi**2/3*(-Hr1(0)-Hr1(1)+GYZ1(2)-GYZ1(0) )


      l(0)=2*Hr1(0)*GYZ1(0)+2*Hr2(0,0)+2*Hr2(1,0)+2*GYZ2(0,0)
     $    -2*GYZ2(1,0)+pi**2/3

      l(-1)=-2*Hr1(0)-2*GYZ1(0)

      l(-2)=2.0d0


      do i=-8,-3
       Box6(i)=0.0d0
      end do

      Box6(-2)=c(0)*l(-2)

      Box6(-1)=c(0)*l(-1)+c(1)*l(-2)

      Box6(0)=c(0)*l(0)+c(1)*l(-1)+c(2)*l(-2)

      Box6(1)=c(0)*l(1)+c(1)*l(0)+c(2)*l(-1)+c(3)*l(-2)

      Box6(2)=c(0)*l(2)+c(1)*l(1)+c(2)*l(0)+c(3)*l(-1)+c(4)*l(-2)


      do i=-8,2
       Box6(i)=Box6(i)+(Buby(i+1)+Bubz(i+1)-Bub1(i+1))/(1-y-z)
      end do

      end subroutine


************************************************************************
****      2x0 two-loop finite part coded from hep-ph/0112081v1      ****
************************************************************************
      real*8 function fin20(i,y,z)
      implicit none

      integer i !A, B, C, D, E, F, G - (1-7)
      real*8 y, z, T
      real*8 GYZ1(0:3),GYZ2(0:3,0:3),GYZ3(0:3,0:3,0:3),
     $ GYZ4(0:3,0:3,0:3,0:3), HZr1(0:1),HZr2(0:1,0:1),HZr3(0:1,0:1,0:1),
     $ HZr4(0:1,0:1,0:1,0:1)
      real*8 NC, nf, TR, pi, eg, z3
      common /piez/ pi, z3
      common /const/ NC, nf, TR, eg

      T= y/z+z/y+2*(1/(y*z)-1/y-1/z)

      call tdhpl(y,z,4,GYZ1,GYZ2,GYZ3,GYZ4,HZr1,HZr2,HZr3,HZr4)

      if(i.eq.1) then
      fin20=
     $(z)/(12*y)*(
     $   2*pi**2 
     $ + 6*pi**2* HZr1(0)
     $ - 12*pi**2* GYZ1(1)
     $ - 72*z3
     $ + 8* HZr1(0)
     $ - 36* HZr1(0)* GYZ2(1,0)
     $ - 36* HZr3(0,1,0)
     $ + 39* HZr2(1,0)
     $ + 39* GYZ2(1,0)
     $ + 72* GYZ3(1,1,0)
     $ )
     $ + (1.0d0)/(2*y*(y+z))*(
     $   17* HZr2(1,0)
     $ + 17* GYZ2(1,0)
     $ )
     $ + (1.0d0)/(36*y)*(
     $ - 12*pi**2
     $ - 24*pi**2* HZr1(0)
     $ + 48*pi**2* GYZ1(1)
     $ + 288*z3
     $ + 457
     $ - 84* HZr1(0)
     $ - 36* HZr1(0)* GYZ1(0)
     $ + 144* HZr1(0)* GYZ2(1,0)
     $ + 144* HZr3(0,1,0)
     $ - 306* HZr2(1,0)
     $ - 192* GYZ1(0)
     $ - 234* GYZ2(1,0)
     $ - 288* GYZ3(1,1,0)
     $ )
     $ + (z)/(36*(1-y)**2)*(
     $ -pi**2
     $ + 6*pi**2* HZr1(0)
     $ + 6*pi**2* HZr1(1)
     $ - 6*pi**2* GYZ1(2)
     $ + 18*pi**2* GYZ1(0)
     $ - 12*pi**2* GYZ1(1)
     $ + 36*z3
     $ - 36* HZr1(0)* GYZ2(2,0)
     $ + 60* HZr1(0)* GYZ1(0)
     $ + 72* HZr1(0)* GYZ2(0,0)
     $ + 36* HZr3(0,1,0)
     $ - 36* HZr2(1,0)* GYZ1(2)
     $ + 36* HZr2(1,0)* GYZ1(0)
     $ + 36* HZr3(1,1,0)
     $ + 36* GYZ3(2,1,0)
     $ - 355* GYZ1(0)
     $ + 270* GYZ2(0,0)
     $ - 108* GYZ3(0,1,0)
     $ + 6* GYZ2(1,0)
     $ - 72* GYZ3(1,0,0)
     $ + 72* GYZ3(1,1,0)
     $ )
     $ + (z)/(36*(1-y))*(
     $ - 33*pi**2 
     $ + 18*pi**2* HZr1(0)
     $ + 18*pi**2* HZr1(1)
     $ - 18*pi**2* GYZ1(2)
     $ + 54*pi**2* GYZ1(0)
     $ - 36*pi**2* GYZ1(1)
     $ + 108*z3
     $ - 277
     $ + 60* HZr1(0)
     $ - 108* HZr1(0)* GYZ2(2,0)
     $ + 216* HZr1(0)* GYZ1(0)
     $ + 216* HZr1(0)* GYZ2(0,0)
     $ + 108* HZr3(0,1,0)
     $ + 36* HZr2(1,0)
     $ - 108* HZr2(1,0)* GYZ1(2)
     $ + 108* HZr2(1,0)* GYZ1(0)
     $ + 108* HZr3(1,1,0)
     $ + 108* GYZ3(2,1,0)
     $ - 615* GYZ1(0)
     $ + 594* GYZ2(0,0)
     $ - 324* GYZ3(0,1,0)
     $ + 198* GYZ2(1,0)
     $ - 216* GYZ3(1,0,0)
     $ + 216* GYZ3(1,1,0)
     $ )
     $ + (z)/((y+z)**3)*(
     $   (11*pi**2)/(2.0d0)* HZr1(1)
     $ - (11*pi**2)/(2.0d0)* GYZ1(2)
     $ - 33* HZr1(0)* GYZ2(2,0)
     $ - 33* HZr3(0,1,0)
     $ - 33* HZr2(1,0)
     $ - 33* HZr2(1,0)* GYZ1(2)
     $ + 33* HZr2(1,0)* GYZ1(0)
     $ + 33* HZr3(1,1,0)
     $ + 33* GYZ3(2,1,0)
     $ + 33* GYZ3(0,1,0)
     $ - 33* GYZ2(1,0)
     $ )
     $ + (z)/(2*(y+z)**2)*(
     $ - 11*pi**2
     $ - (22*pi**2)/(3.0d0)* HZr1(1)
     $ + (22*pi**2)/(3.0d0)* GYZ1(2)
     $ + 33* HZr1(0)
     $ + 44* HZr1(0)* GYZ2(2,0)
     $ - 66* HZr1(0)* GYZ1(0)
     $ + 44* HZr3(0,1,0)
     $ - 22* HZr2(1,0)
     $ + 44* HZr2(1,0)* GYZ1(2)
     $ - 44* HZr2(1,0)* GYZ1(0)
     $ - 44* HZr3(1,1,0)
     $ - 44* GYZ3(2,1,0)
     $ - 33* GYZ1(0)
     $ - 44* GYZ3(0,1,0)
     $ + 110* GYZ2(1,0)
     $ )
     $ + (z)/(2*(y+z))*(
     $   (11*pi**2)/(6.0d0)
     $ - 11* HZr1(0)
     $ + 11* HZr1(0)* GYZ1(0)
     $ + 11* HZr2(1,0)
     $ + 11* GYZ1(0)
     $ - 11* GYZ2(1,0)
     $ )
     $ + (z**2)/((y+z)**4)*(
     $ - (11*pi**2)/(2.0d0)* HZr1(1)
     $ + (11*pi**2)/(2.0d0)* GYZ1(2)
     $ + 33* HZr1(0)* GYZ2(2,0)
     $ + 33* HZr3(0,1,0)
     $ + 33* HZr2(1,0)* GYZ1(2)
     $ - 33* HZr2(1,0)* GYZ1(0)
     $ - 33* HZr3(1,1,0)
     $ - 33* GYZ3(2,1,0)
     $ - 33* GYZ3(0,1,0)
     $ )
     $ + (z**2)/((y+z)**3)*(
     $   (11*pi**2)/(2.0d0)
     $ + (11*pi**2)/(3.0d0)* HZr1(1)
     $ - (11*pi**2)/(3.0d0)* GYZ1(2)
     $ - 22* HZr1(0)* GYZ2(2,0)
     $ + 33* HZr1(0)* GYZ1(0)
     $ - 22* HZr3(0,1,0)
     $ + 33* HZr2(1,0)
     $ - 22* HZr2(1,0)* GYZ1(2)
     $ + 22* HZr2(1,0)* GYZ1(0)
     $ + 22* HZr3(1,1,0)
     $ + 22* GYZ3(2,1,0)
     $ + 22* GYZ3(0,1,0)
     $ - 33* GYZ2(1,0)
     $ )
     $ + (z**2)/(2*(y+z)**2)*(
     $ - (11*pi**2)/(6.0d0)
     $ - 11* HZr1(0)* GYZ1(0)
     $ - 11* HZr2(1,0)
     $ + 11* GYZ2(1,0)
     $ )
     $ + (1.0d0)/(18*(1-y))*(
     $ + 23*pi**2 
     $ - 12*pi**2* HZr1(0)
     $ - 12*pi**2* HZr1(1)
     $ + 12*pi**2* GYZ1(2)
     $ - 36*pi**2* GYZ1(0)
     $ + 24*pi**2* GYZ1(1)
     $ - 72*z3
     $ + 72* HZr1(0)* GYZ2(2,0)
     $ - 120* HZr1(0)* GYZ1(0)
     $ - 144* HZr1(0)* GYZ2(0,0)
     $ - 72* HZr3(0,1,0)
     $ - 18* HZr2(1,0)
     $ + 72* HZr2(1,0)* GYZ1(2)
     $ - 72* HZr2(1,0)* GYZ1(0)
     $ - 72* HZr3(1,1,0)
     $ - 72* GYZ3(2,1,0)
     $ + 515* GYZ1(0)
     $ - 432* GYZ2(0,0)
     $ + 216* GYZ3(0,1,0)
     $ - 138* GYZ2(1,0)
     $ + 144* GYZ3(1,0,0)
     $ - 144* GYZ3(1,1,0)
     $ )
     $ + (1.0d0)/((y+z)**2)*(
     $ - (7*pi**2)/(3.0d0)* HZr1(1)
     $ + (7*pi**2)/(3.0d0)* GYZ1(2)
     $ + 14* HZr1(0)* GYZ2(2,0)
     $ + 14* HZr3(0,1,0)
     $ + 14* HZr2(1,0)* GYZ1(2)
     $ - 14* HZr2(1,0)* GYZ1(0)
     $ - 14* HZr3(1,1,0)
     $ - 14* GYZ3(2,1,0)
     $ - 14* GYZ3(0,1,0)
     $ )
     $ + (1.0d0)/(4*(y+z))*(
     $   (28*pi**2)/(3.0d0)
     $ + (14*pi**2)/(3.0d0)* HZr1(1)
     $ - (14*pi**2)/(3.0d0)* GYZ1(2)
     $ - 22
     $ + 11* HZr1(0)
     $ - 28* HZr1(0)* GYZ2(2,0)
     $ + 56* HZr1(0)* GYZ1(0)
     $ - 28* HZr3(0,1,0)
     $ + 56* HZr2(1,0)
     $ - 28* HZr2(1,0)* GYZ1(2)
     $ + 28* HZr2(1,0)* GYZ1(0)
     $ + 28* HZr3(1,1,0)
     $ + 28* GYZ3(2,1,0)
     $ + 11* GYZ1(0)
     $ + 28* GYZ3(0,1,0)
     $ - 56* GYZ2(1,0)
     $ )
     $ + (T*pi**2)/(216.0d0)*(
     $ - 1045
     $ + 147* HZr1(0)
     $ + 36* HZr1(0)* GYZ1(2)
     $ + 108* HZr1(0)* GYZ1(0)
     $ - 36* HZr1(0)* GYZ1(1)
     $ + 72* HZr2(0,1)
     $ + 54* HZr1(1)
     $ + 72* HZr1(1)* GYZ1(2)
     $ - 72* HZr1(1)* GYZ1(3)
     $ - 36* HZr1(1)* GYZ1(1)
     $ + 72* HZr2(1,0)
     $ + 36* HZr2(1,1)
     $ - 186* GYZ1(2)
     $ + 36* GYZ2(2,0)
     $ - 72* GYZ2(2,1)
     $ + 72* GYZ2(3,2)
     $ - 72* GYZ2(0,2)
     $ + 147* GYZ1(0)
     $ - 72* GYZ2(0,1)
     $ + 36* GYZ2(1,2)
     $ + 132* GYZ1(1)
     $ - 108* GYZ2(1,0)
     $ + 72* GYZ2(1,1)
     $ )
     $ + (T)/(216.0d0)*(
     $ - (99931.0d0)/(12.0d0)
     $ + (132*pi**4)/(5.0d0)
     $ + 4776*z3
     $ - 216*z3* HZr1(0)
     $ + 1080*z3* HZr1(1)
     $ - 864*z3* GYZ1(2)
     $ - 216*z3* GYZ1(0)
     $ - 216*z3* GYZ1(1)
     $ + 304* HZr1(0)
     $ - 1116* HZr1(0)* GYZ2(2,0)
     $ - 216* HZr1(0)* GYZ3(2,1,0)
     $ + 432* HZr1(0)* GYZ3(3,2,0)
     $ - 432* HZr1(0)* GYZ3(0,2,0)
     $ - 144* HZr1(0)* GYZ1(0)
     $ + 1512* HZr1(0)* GYZ2(0,0)
     $ - 216* HZr1(0)* GYZ3(0,1,0)
     $ + 216* HZr1(0)* GYZ3(1,2,0)
     $ - 36* HZr1(0)* GYZ2(1,0)
     $ - 432* HZr1(0)* GYZ3(1,0,0)
     $ + 1920* HZr2(0,0)
     $ + 1512* HZr2(0,0)* GYZ1(0)
     $ + 432* HZr2(0,0)* GYZ2(0,0)
     $ + 864* HZr4(0,0,1,0)
     $ + 1008* HZr3(0,1,0)
     $ - 216* HZr3(0,1,0)* GYZ1(2)
     $ + 432* HZr3(0,1,0)* GYZ1(3)
     $ + 216* HZr3(0,1,0)* GYZ1(0)
     $ - 216* HZr3(0,1,0)* GYZ1(1)
     $ + 432* HZr4(0,1,1,0)
     $ - 1095* HZr2(1,0)
     $ - 1116* HZr2(1,0)* GYZ1(2)
     $ + 216* HZr2(1,0)* GYZ2(2,0)
     $ + 432* HZr2(1,0)* GYZ2(3,2)
     $ - 432* HZr2(1,0)* GYZ2(3,0)
     $ - 432* HZr2(1,0)* GYZ2(0,2)
     $ + 1152* HZr2(1,0)* GYZ1(0)
     $ + 216* HZr2(1,0)* GYZ2(1,2)
     $ - 216* HZr2(1,0)* GYZ2(1,0)
     $ + 1512* HZr3(1,0,0)
     $ + 432* HZr3(1,0,0)* GYZ1(0)
     $ + 864* HZr4(1,0,1,0)
     $ + 324* HZr3(1,1,0)
     $ + 432* HZr3(1,1,0)* GYZ1(2)
     $ - 432* HZr3(1,1,0)* GYZ1(3)

     $ - 216* HZr3(1,1,0)* GYZ1(1)
     $ + 432* HZr4(1,1,0,0)
     $ + 216* HZr4(1,1,1,0)
     $ + 216* GYZ4(2,0,1,0)

     $ + 1116* GYZ3(2,1,0)
     $ + 432* GYZ4(2,1,1,0)
     $ - 432* GYZ4(3,2,1,0)
     $ - 432* GYZ4(3,0,1,0)

     $ + 432* GYZ4(0,2,1,0)
     $ + 304* GYZ1(0)
     $ + 1920* GYZ2(0,0)
     $ - 432* GYZ4(0,0,1,0)
     $ - 1008* GYZ3(0,1,0)

     $ + 432* GYZ4(0,1,1,0)
     $ - 216* GYZ4(1,2,1,0)
     $ + 1095* GYZ2(1,0)
     $ - 1512* GYZ3(1,0,0)
     $ + 648* GYZ4(1,0,1,0)

     $ - 792* GYZ3(1,1,0)
     $ + 432* GYZ4(1,1,0,0)
     $ - 432* GYZ4(1,1,1,0)
     $ )
     $ + (1.0d0)/(2.0d0)*(
     $ - (7*pi**2)/(3.0d0)
     $ + 1
     $ - 14* HZr1(0)* GYZ1(0)

     $ - 14* HZr2(1,0)
     $ + 14* GYZ2(1,0)
     $ )


      else if(i.eq.2) then
      fin20=
     $(z)/(y**2)*(
     $ - 3* HZr1(0)* GYZ1(2)
     $ - 3* HZr1(1)* GYZ1(3)
     $ + 3* GYZ2(3,2)
     $)
     $ + (z**2)/(y**2)*(
     $   HZr1(0)* GYZ1(2)

     $ + HZr1(1)* GYZ1(3)
     $ - GYZ2(3,2)
     $)
     $ + (1.0d0)/(y**2)*(
     $   2* HZr1(0)* GYZ1(2)
     $ + 2* HZr1(1)* GYZ1(3)
     $ - 2* GYZ2(3,2)
     $)

     $ + (z*pi**2)/(18*y)*(
     $   3* HZr1(0)
     $ - 24* HZr1(0)* GYZ1(2)
     $ + 21* HZr1(1)
     $ - 12* HZr1(1)* GYZ1(2)
     $ + 12* GYZ2(2,2)

     $ + GYZ1(2)
     $ + 12* GYZ2(2,1)
     $ + 6* GYZ2(0,2)
     $ + 36* GYZ1(1)
     $)
     $ + (z)/(9*y)*(
     $ - 27*z3
     $ - 90*z3* GYZ1(2)

     $ - 36* HZr1(0)
     $ - 84* HZr1(0)* GYZ2(2,2)
     $ + 18* HZr1(0)* GYZ3(2,2,0)

     $ + 18* HZr1(0)* GYZ3(2,3,2)
     $ + 152* HZr1(0)* GYZ1(2)
     $ - 57* HZr1(0)* GYZ2(2,0)

     $ + 36* HZr1(0)* GYZ3(3,2,2)
     $ + 78* HZr1(0)* GYZ2(3,2)
     $ - 18* HZr1(0)* GYZ3(3,0,2)

     $ - 18* HZr1(0)* GYZ3(0,2,2)
     $ - 3* HZr1(0)* GYZ2(0,2)
     $ + 36* HZr1(0)* GYZ3(0,2,0)

     $ + 18* HZr1(0)* GYZ3(0,3,2)
     $ + 54* HZr1(0)* GYZ2(1,0)
     $ - 18* HZr2(0,0)* GYZ2(2,2)

     $ - 36* HZr2(0,0)* GYZ1(2)
     $ - 18* HZr2(0,0)* GYZ2(2,0)
     $ - 18* HZr2(0,0)* GYZ2(0,2)

     $ - 9* HZr3(0,0,1)
     $ + 54* HZr3(0,0,1)* GYZ1(2)
     $ - 72* HZr3(0,0,1)* GYZ1(3)
     $ - 9* HZr2(0,1)

     $ + 54* HZr2(0,1)* GYZ2(2,3)
     $ + 9* HZr2(0,1)* GYZ1(2)
     $ - 18* HZr2(0,1)* GYZ2(2,0)

     $ - 108* HZr2(0,1)* GYZ2(3,3)
     $ + 6* HZr2(0,1)* GYZ1(3)
     $ + 18* HZr2(0,1)* GYZ2(3,0)
     $ + 18* HZr2(0,1)* GYZ2(0,2)

     $ - 36* HZr2(0,1)* GYZ2(0,3)
     $ + 9* HZr3(0,1,0)
     $ - 36* HZr3(0,1,0)* GYZ1(2)
     $ + 18* HZr3(0,1,0)* GYZ1(3)

     $ + 72* HZr1(1)* GYZ3(2,3,3)
     $ - 75* HZr1(1)* GYZ2(2,3)
     $ - 18* HZr1(1)* GYZ3(2,3,0)

     $ - 18* HZr1(1)* GYZ3(2,0,3)
     $ + 36* HZr1(1)* GYZ3(3,2,3)
     $ - 84* HZr1(1)* GYZ2(3,2)

     $ + 36* HZr1(1)* GYZ3(3,3,2)
     $ - 108* HZr1(1)* GYZ3(3,3,3)
     $ + 84* HZr1(1)* GYZ2(3,3)

     $ + 36* HZr1(1)* GYZ3(3,3,0)
     $ + 143* HZr1(1)* GYZ1(3)
     $ + 24* HZr1(1)* GYZ2(3,0)
     $ - 18* HZr1(1)* GYZ3(0,3,2)

     $ - 18* HZr1(1)* GYZ3(0,3,3)
     $ - 3* HZr1(1)* GYZ2(0,3)
     $ + 36* HZr1(1)* GYZ3(0,3,0)
     $ + 9* HZr1(1)* GYZ1(0)

     $ + 27* HZr1(1)* GYZ2(1,0)
     $ + 36* HZr2(1,0)* GYZ2(2,2)
     $ - 18* HZr2(1,0)* GYZ2(2,3)

     $ + 57* HZr2(1,0)* GYZ1(2)
     $ - 18* HZr2(1,0)* GYZ2(2,0)
     $ - 18* HZr2(1,0)* GYZ2(3,2)

     $ - 78* HZr2(1,0)* GYZ1(3)
     $ + 36* HZr2(1,0)* GYZ2(0,2)
     $ - 18* HZr2(1,0)* GYZ2(0,3)
     $ - 18* HZr3(1,0,0)* GYZ1(2)

     $ - 36* HZr3(1,0,1)
     $ - 18* HZr3(1,0,1)* GYZ1(2)
     $ - 36* HZr2(1,1)* GYZ2(3,3)
     $ + 84* HZr2(1,1)* GYZ1(3)

     $ + 18* HZr2(1,1)* GYZ2(0,3)
     $ - 18* HZr3(1,1,0)
     $ - 54* HZr3(1,1,0)* GYZ1(2)
     $ + 18* HZr3(1,1,0)* GYZ1(3)

     $ - 18* GYZ4(2,2,1,0)
     $ + 75* GYZ3(2,3,2)
     $ + 18* GYZ4(2,3,2,0)

     $ - 72* GYZ4(2,3,3,2)
     $ + 18* GYZ4(2,3,0,2)
     $ + 18* GYZ4(2,0,3,2)
     $ - 9* GYZ2(2,0)

     $ + 18* GYZ4(2,0,1,0)
     $ + 54* GYZ3(2,1,0)
     $ - 36* GYZ4(2,1,1,0)
     $ + 84* GYZ3(3,2,2)

     $ - 36* GYZ4(3,2,3,2)
     $ - 143* GYZ2(3,2)
     $ - 24* GYZ3(3,2,0)
     $ + 18* GYZ4(3,2,1,0)

     $ - 36* GYZ4(3,3,2,2)
     $ - 84* GYZ3(3,3,2)
     $ - 36* GYZ4(3,3,2,0)

     $ + 108* GYZ4(3,3,3,2)
     $ - 36* GYZ4(3,3,0,2)
     $ - 24* GYZ3(3,0,2)
     $ + 18* GYZ4(3,0,1,0)

     $ - 9* GYZ2(0,2)
     $ + 18* GYZ4(0,3,2,2)
     $ + 3* GYZ3(0,3,2)
     $ - 36* GYZ4(0,3,2,0)

     $ + 18* GYZ4(0,3,3,2)
     $ - 36* GYZ4(0,3,0,2)
     $ + 27* GYZ3(0,1,0)
     $ - 27* GYZ3(1,2,0)

     $ - 27* GYZ3(1,0,2)
     $ + 9* GYZ2(1,0)
     $ - 108* GYZ3(1,1,0)
     $)
     $ + (z**2)/(y)*(
     $   2* HZr1(0)* GYZ1(2)
     $ + 2* HZr1(1)* GYZ1(3)

     $ - 2* GYZ2(3,2)
     $)
     $ + (1.0d0)/(y*(y+z))*(
     $   2* HZr1(0)* GYZ1(2)
     $ + 6* HZr1(0)* GYZ2(3,2)
     $ - 6* HZr1(0)* GYZ2(0,2)

     $ - 2* HZr2(0,1)
     $ - 6* HZr2(0,1)* GYZ1(3)
     $ + 6* HZr1(1)* GYZ2(3,0)
     $ - 6* HZr1(1)* GYZ2(0,3)
     $ + 2* HZr1(1)* GYZ1(0)

     $ + 6* HZr2(1,0)* GYZ1(2)
     $ - 6* HZr2(1,0)* GYZ1(3)
     $ - 6* HZr3(1,1,0)
     $ - 2* GYZ2(2,0)
     $ + 6* GYZ3(2,1,0)

     $ - 6* GYZ3(3,2,0)
     $ - 6* GYZ3(3,0,2)
     $ - 2* GYZ2(0,2)
     $ + 6* GYZ3(0,3,2)
     $ + 6* GYZ3(0,1,0)

     $ + 2* GYZ2(1,0)
     $)
     $ + (pi**2)/(9*y)*(
     $ - 3
     $ + 3* HZr1(0)
     $ + 24* HZr1(0)* GYZ1(2)
     $ - 15* HZr1(1)
     $ + 12* HZr1(1)* GYZ1(2)

     $ - 12* GYZ2(2,2)
     $ - 7* GYZ1(2)
     $ - 12* GYZ2(2,1)
     $ - 6* GYZ2(0,2)
     $ - 24* GYZ1(1)
     $)
     $ + (1.0d0)/(9*y)*(
     $   (45.0d0)/(2.0d0)* HZr1(1)

     $ - (45.0d0)/(2.0d0)* GYZ1(2)
     $ - 54*z3
     $ + 180*z3* GYZ1(2)
     $ + 139
     $ + 57* HZr1(0)
     $ + 132* HZr1(0)* GYZ2(2,2)

     $ - 36* HZr1(0)* GYZ3(2,2,0)
     $ - 36* HZr1(0)* GYZ3(2,3,2)
     $ - 250* HZr1(0)* GYZ1(2)

     $ + 96* HZr1(0)* GYZ2(2,0)
     $ - 72* HZr1(0)* GYZ3(3,2,2)
     $ - 210* HZr1(0)* GYZ2(3,2)

     $ + 36* HZr1(0)* GYZ3(3,0,2)
     $ + 36* HZr1(0)* GYZ3(0,2,2)
     $ + 96* HZr1(0)* GYZ2(0,2)

     $ - 72* HZr1(0)* GYZ3(0,2,0)
     $ - 36* HZr1(0)* GYZ3(0,3,2)
     $ - 9* HZr1(0)* GYZ1(0)
     $ - 72* HZr1(0)* GYZ2(1,0)

     $ + 36* HZr2(0,0)* GYZ2(2,2)
     $ + 144* HZr2(0,0)* GYZ1(2)
     $ + 36* HZr2(0,0)* GYZ2(2,0)

     $ + 36* HZr2(0,0)* GYZ2(0,2)
     $ - 18* HZr3(0,0,1)
     $ - 108* HZr3(0,0,1)* GYZ1(2)
     $ + 144* HZr3(0,0,1)* GYZ1(3)

     $ + 36* HZr2(0,1)
     $ - 108* HZr2(0,1)* GYZ2(2,3)
     $ + 36* HZr2(0,1)* GYZ2(2,0)
     $ + 216* HZr2(0,1)* GYZ2(3,3)

     $ + 42* HZr2(0,1)* GYZ1(3)
     $ - 36* HZr2(0,1)* GYZ2(3,0)
     $ - 36* HZr2(0,1)* GYZ2(0,2)
     $ + 72* HZr2(0,1)* GYZ2(0,3)

     $ + 18* HZr3(0,1,0)
     $ + 72* HZr3(0,1,0)* GYZ1(2)
     $ - 36* HZr3(0,1,0)* GYZ1(3)
     $ - 144* HZr1(1)* GYZ3(2,3,3)

     $ + 132* HZr1(1)* GYZ2(2,3)
     $ + 36* HZr1(1)* GYZ3(2,3,0)
     $ + 36* HZr1(1)* GYZ3(2,0,3)

     $ - 72* HZr1(1)* GYZ3(3,2,3)
     $ + 132* HZr1(1)* GYZ2(3,2)
     $ - 72* HZr1(1)* GYZ3(3,3,2)

     $ + 216* HZr1(1)* GYZ3(3,3,3)
     $ - 168* HZr1(1)* GYZ2(3,3)
     $ - 72* HZr1(1)* GYZ3(3,3,0)

     $ - 214* HZr1(1)* GYZ1(3)
     $ - 30* HZr1(1)* GYZ2(3,0)
     $ + 36* HZr1(1)* GYZ3(0,3,2)
     $ + 36* HZr1(1)* GYZ3(0,3,3)

     $ + 96* HZr1(1)* GYZ2(0,3)
     $ - 72* HZr1(1)* GYZ3(0,3,0)
     $ - 27* HZr1(1)* GYZ1(0)
     $ - 18* HZr1(1)* GYZ2(1,0)
     $ - 18* HZr2(1,0)

     $ - 72* HZr2(1,0)* GYZ2(2,2)
     $ + 36* HZr2(1,0)* GYZ2(2,3)
     $ - 168* HZr2(1,0)* GYZ1(2)

     $ + 36* HZr2(1,0)* GYZ2(2,0)
     $ + 36* HZr2(1,0)* GYZ2(3,2)
     $ + 210* HZr2(1,0)* GYZ1(3)

     $ - 72* HZr2(1,0)* GYZ2(0,2)
     $ + 36* HZr2(1,0)* GYZ2(0,3)
     $ + 36* HZr3(1,0,0)* GYZ1(2)
     $ + 36* HZr3(1,0,1)

     $ + 36* HZr3(1,0,1)* GYZ1(2)
     $ + 72* HZr2(1,1)* GYZ2(3,3)
     $ - 132* HZr2(1,1)* GYZ1(3)
     $ - 36* HZr2(1,1)* GYZ2(0,3)

     $ + 72* HZr3(1,1,0)
     $ + 108* HZr3(1,1,0)* GYZ1(2)
     $ - 36* HZr3(1,1,0)* GYZ1(3)
     $ + 36* GYZ4(2,2,1,0)

     $ - 132* GYZ3(2,3,2)
     $ - 36* GYZ4(2,3,2,0)
     $ + 144* GYZ4(2,3,3,2)

     $ - 36* GYZ4(2,3,0,2)
     $ - 36* GYZ4(2,0,3,2)
     $ + 27* GYZ2(2,0)
     $ - 36* GYZ4(2,0,1,0)

     $ - 108* GYZ3(2,1,0)
     $ + 72* GYZ4(2,1,1,0)
     $ - 132* GYZ3(3,2,2)
     $ + 72* GYZ4(3,2,3,2)

     $ + 214* GYZ2(3,2)
     $ + 30* GYZ3(3,2,0)
     $ - 36* GYZ4(3,2,1,0)
     $ + 72* GYZ4(3,3,2,2)

     $ + 168* GYZ3(3,3,2)
     $ + 72* GYZ4(3,3,2,0)
     $ - 216* GYZ4(3,3,3,2)

     $ + 72* GYZ4(3,3,0,2)
     $ + 30* GYZ3(3,0,2)
     $ - 36* GYZ4(3,0,1,0)
     $ + 27* GYZ2(0,2)

     $ - 36* GYZ4(0,3,2,2)
     $ - 96* GYZ3(0,3,2)
     $ + 72* GYZ4(0,3,2,0)
     $ - 36* GYZ4(0,3,3,2)

     $ + 72* GYZ4(0,3,0,2)
     $ - 48* GYZ1(0)
     $ - 54* GYZ3(0,1,0)
     $ + 18* GYZ3(1,2,0)
     $ + 18* GYZ3(1,0,2)

     $ - 18* GYZ2(1,0)
     $ + 144* GYZ3(1,1,0)
     $)
     $ + (z)/(9*(1-y)**2)*(
     $ - 11*pi**2
     $ - (27*pi**2)/(2.0d0)* HZr1(0)
     $ - (27*pi**2)/(2.0d0)* HZr1(1)
     $ + (27*pi**2)/(2.0d0)* GYZ1(2)

     $ - 6*pi**2* GYZ1(0)
     $ + 6*pi**2* GYZ1(1)
     $ + 225*z3
     $ + (81.0d0)/(2.0d0)* HZr1(1)* GYZ1(0)
     $ - (81.0d0)/(2.0d0)* GYZ2(2,0)
     $ - (81.0d0)/(2.0d0)* GYZ2(0,2)

     $ + 9* HZr1(0)* GYZ2(2,0)
     $ - 72* HZr1(0)* GYZ2(0,2)
     $ + 12* HZr1(0)* GYZ1(0)
     $ - 9* HZr1(0)* GYZ2(0,0)
     $ + 81* HZr3(0,0,1)

     $ + 36* HZr2(0,1)
     $ - 81* HZr2(0,1)* GYZ1(2)
     $ - 9* HZr2(0,1)* GYZ1(0)
     $ - 81* HZr3(0,1,0)
     $ - 81* HZr1(1)* GYZ2(2,3)

     $ + 36* HZr1(1)* GYZ1(3)
     $ - 81* HZr1(1)* GYZ2(0,3)
     $ + 27* HZr1(1)* GYZ2(0,0)
     $ - 36* HZr2(1,0)
     $ + 81* HZr2(1,0)* GYZ1(2)

     $ - 9* HZr2(1,0)* GYZ1(0)
     $ + 81* HZr3(1,0,1)
     $ - 81* HZr3(1,1,0)
     $ + 81* GYZ3(2,3,2)
     $ - 27* GYZ3(2,0,0)

     $ - 36* GYZ2(3,2)
     $ - 27* GYZ3(0,2,0)
     $ + 81* GYZ3(0,3,2)
     $ + 37* GYZ1(0)
     $ - 27* GYZ3(0,0,2)

     $ - 45* GYZ2(0,0)
     $ + 45* GYZ3(0,1,0)
     $ + 30* GYZ2(1,0)
     $ + 36* GYZ3(1,0,0)
     $ - 36* GYZ3(1,1,0)
     $)

     $ + (z)/(18*(1-y))*(
     $ - 9*pi**2* HZr1(0)
     $ - 9*pi**2* HZr1(1)
     $ + 9*pi**2* GYZ1(2)
     $ - 36*pi**2* GYZ1(0)
     $ + 36*pi**2* GYZ1(1)
     $ + 54*z3
     $ + 62

     $ + 60* HZr1(0)
     $ - 144* HZr1(0)* GYZ1(2)
     $ + 54* HZr1(0)* GYZ2(2,0)
     $ - 108* HZr1(0)* GYZ1(0)
     $ - 54* HZr1(0)* GYZ2(0,0)

     $ + 54* HZr3(0,0,1)
     $ - 18* HZr2(0,1)
     $ - 54* HZr2(0,1)* GYZ1(2)
     $ - 54* HZr2(0,1)* GYZ1(0)
     $ - 54* HZr3(0,1,0)

     $ + 117* HZr1(1)
     $ - 54* HZr1(1)* GYZ2(2,3)
     $ - 162* HZr1(1)* GYZ1(3)
     $ - 54* HZr1(1)* GYZ2(0,3)
     $ + 99* HZr1(1)* GYZ1(0)

     $ + 162* HZr1(1)* GYZ2(0,0)
     $ - 18* HZr2(1,0)
     $ + 54* HZr2(1,0)* GYZ1(2)
     $ - 54* HZr2(1,0)* GYZ1(0)
     $ + 54* HZr3(1,0,1)

     $ - 54* HZr3(1,1,0)
     $ + 54* GYZ3(2,3,2)
     $ - 117* GYZ1(2)
     $ - 99* GYZ2(2,0)
     $ - 162* GYZ3(2,0,0)

     $ + 162* GYZ2(3,2)
     $ - 99* GYZ2(0,2)
     $ - 162* GYZ3(0,2,0)
     $ + 54* GYZ3(0,3,2)
     $ + 132* GYZ1(0)

     $ - 162* GYZ3(0,0,2)
     $ - 198* GYZ2(0,0)
     $ + 270* GYZ3(0,1,0)
     $ + 18* GYZ2(1,0)
     $ + 216* GYZ3(1,0,0)
     $ - 216* GYZ3(1,1,0)
     $)

     $ + (z)/((y+z)**3)*(
     $ - 2*pi**2* HZr1(1)
     $ + 2*pi**2* GYZ1(2)
     $ + 12* HZr1(0)* GYZ2(2,0)
     $ + 12* HZr3(0,1,0)
     $ + 12* HZr2(1,0)

     $ + 12* HZr2(1,0)* GYZ1(2)
     $ - 12* HZr2(1,0)* GYZ1(0)
     $ - 12* HZr3(1,1,0)
     $ - 12* GYZ3(2,1,0)
     $ - 12* GYZ3(0,1,0)

     $ + 12* GYZ2(1,0)
     $)
     $ + (z)/((y+z)**2)*(
     $   2*pi**2
     $ + (4*pi**2)/(3.0d0)* HZr1(1)
     $ - (4*pi**2)/(3.0d0)* GYZ1(2)
     $ - 6* HZr1(0)
     $ - 8* HZr1(0)* GYZ2(2,0)

     $ + 12* HZr1(0)* GYZ1(0)
     $ - 8* HZr3(0,1,0)
     $ + 4* HZr2(1,0)
     $ - 8* HZr2(1,0)* GYZ1(2)
     $ + 8* HZr2(1,0)* GYZ1(0)

     $ + 8* HZr3(1,1,0)
     $ + 8* GYZ3(2,1,0)
     $ + 6* GYZ1(0)
     $ + 8* GYZ3(0,1,0)
     $ - 20* GYZ2(1,0)
     $)
     $ + (z)/(y+z)*(
     $ - (pi**2)/(3.0d0)

     $ + 2* HZr1(0)
     $ - 2* HZr1(0)* GYZ1(0)
     $ - 2* HZr2(1,0)
     $ - 2* GYZ1(0)
     $ + 2* GYZ2(1,0)
     $)
     $ + (z**2)/((1-y)**3)*(
     $ + (2*pi**2)/(3.0d0)* HZr1(0)

     $ + (2*pi**2)/(3.0d0)* HZr1(1)
     $ - (2*pi**2)/(3.0d0)* GYZ1(2)
     $ - 12*z3
     $ + 4* HZr1(0)* GYZ2(0,2)
     $ - 4* HZr3(0,0,1)
     $ + 4* HZr2(0,1)* GYZ1(2)

     $ + 4* HZr3(0,1,0)
     $ + 4* HZr1(1)* GYZ2(2,3)
     $ + 4* HZr1(1)* GYZ2(0,3)
     $ - 4* HZr2(1,0)* GYZ1(2)
     $ - 4* HZr3(1,0,1)

     $ + 4* HZr3(1,1,0)
     $ - 4* GYZ3(2,3,2)
     $ - 4* GYZ3(0,3,2)
     $)
     $ + (z**2)/((1-y)**2)*(
     $ + 4* HZr1(0)* GYZ1(2)

     $ + 4* HZr1(1)* GYZ1(3)
     $ - 4* GYZ2(3,2)
     $)
     $ + (z**2)/(1-y)*(
     $ + 2* HZr1(0)* GYZ1(2)
     $ + 2* HZr1(1)* GYZ1(3)

     $ - 2* GYZ2(3,2)
     $)
     $ + (z**2)/((y+z)**4)*(
     $   2*pi**2* HZr1(1)
     $ - 2*pi**2* GYZ1(2)
     $ - 12* HZr1(0)* GYZ2(2,0)
     $ - 12* HZr3(0,1,0)

     $ - 12* HZr2(1,0)* GYZ1(2)
     $ + 12* HZr2(1,0)* GYZ1(0)
     $ + 12* HZr3(1,1,0)
     $ + 12* GYZ3(2,1,0)
     $ + 12* GYZ3(0,1,0)
     $)

     $ + (z**2)/((y+z)**3)*(
     $ - 2*pi**2
     $ - (4*pi**2)/(3.0d0)* HZr1(1)
     $ + (4*pi**2)/(3.0d0)* GYZ1(2)
     $ + 8* HZr1(0)* GYZ2(2,0)
     $ - 12* HZr1(0)* GYZ1(0)

     $ + 8* HZr3(0,1,0)
     $ - 12* HZr2(1,0)
     $ + 8* HZr2(1,0)* GYZ1(2)
     $ - 8* HZr2(1,0)* GYZ1(0)
     $ - 8* HZr3(1,1,0)

     $ - 8* GYZ3(2,1,0)
     $ - 8* GYZ3(0,1,0)
     $ + 12* GYZ2(1,0)
     $)
     $ + (z**2)/((y+z)**2)*(
     $   (pi**2)/(3.0d0)
     $ + 2* HZr1(0)* GYZ1(0)
     $ + 2* HZr2(1,0)

     $ - 2* GYZ2(1,0)
     $)
     $ + (1.0d0)/(1-y-z)*(
     $   (pi**2)/(6.0d0)
     $ + HZr1(0)* GYZ1(0)
     $ + HZr2(1,0)
     $ - GYZ2(1,0)
     $)
     $ + (1.0d0)/(9*(1-y))*(
     $   10*pi**2
     $ + 12*pi**2* HZr1(0)

     $ + 9*pi**2* HZr1(1)
     $ - 9*pi**2* GYZ1(2)
     $ + 12*pi**2* GYZ1(0)
     $ - 21*pi**2* GYZ1(1)
     $ - 234*z3
     $ + 54* HZr1(0)* GYZ2(0,2)

     $ + 21* HZr1(0)* GYZ1(0)
     $ + 18* HZr1(0)* GYZ2(0,0)
     $ - 18* HZr1(0)* GYZ2(1,0)
     $ - 90* HZr3(0,0,1)
     $ - 63* HZr2(0,1)

     $ + 90* HZr2(0,1)* GYZ1(2)
     $ + 36* HZr2(0,1)* GYZ1(0)
     $ + 54* HZr3(0,1,0)
     $ + 90* HZr1(1)* GYZ2(2,3)

     $ - 63* HZr1(1)* GYZ1(3)
     $ + 90* HZr1(1)* GYZ2(0,3)
     $ - 18* HZr1(1)* GYZ1(0)
     $ - 72* HZr1(1)* GYZ2(0,0)
     $ + 63* HZr2(1,0)

     $ - 54* HZr2(1,0)* GYZ1(2)
     $ - 90* HZr3(1,0,1)
     $ + 54* HZr3(1,1,0)
     $ - 90* GYZ3(2,3,2)
     $ + 18* GYZ2(2,0)

     $ + 72* GYZ3(2,0,0)
     $ - 36* GYZ3(2,1,0)
     $ + 63* GYZ2(3,2)
     $ + 18* GYZ2(0,2)
     $ + 72* GYZ3(0,2,0)

     $ - 90* GYZ3(0,3,2)
     $ + 55* GYZ1(0)
     $ + 72* GYZ3(0,0,2)
     $ + 108* GYZ2(0,0)
     $ - 108* GYZ3(0,1,0)
     $ + 3* GYZ2(1,0)

     $ - 90* GYZ3(1,0,0)
     $ + 126* GYZ3(1,1,0)
     $)
      fin20 = fin20
     $ + (1.0d0)/(9*(y+z)**2)*(
     $ + 18* HZr1(0)* GYZ2(2,2)
     $ + 27* HZr1(0)* GYZ1(2)

     $ - 18* HZr1(0)* GYZ2(3,2)
     $ - 72* HZr3(0,0,1)
     $ + 39* HZr2(0,1)
     $ + 36* HZr2(0,1)* GYZ1(2)

     $ - 90* HZr2(0,1)* GYZ1(3)
     $ + 18* HZr2(0,1)* GYZ1(0)
     $ + 18* HZr3(0,1,0)
     $ - 18* HZr3(0,1,1)
     $ + 230* HZr1(1)

     $ + 54* HZr1(1)* GYZ2(2,3)
     $ - 102* HZr1(1)* GYZ1(2)
     $ - 18* HZr1(1)* GYZ2(2,0)
     $ + 36* HZr1(1)* GYZ2(3,2)

     $ - 108* HZr1(1)* GYZ2(3,3)
     $ + 66* HZr1(1)* GYZ1(3)
     $ + 18* HZr1(1)* GYZ2(3,0)
     $ - 18* HZr1(1)* GYZ2(0,2)

     $ + 18* HZr1(1)* GYZ2(0,3)
     $ - 27* HZr1(1)* GYZ1(0)
     $ - 27* HZr2(1,0)
     $ - 18* HZr2(1,0)* GYZ1(2)
     $ + 18* HZr2(1,0)* GYZ1(3)

     $ - 36* HZr3(1,0,1)
     $ + 102* HZr2(1,1)
     $ - 36* HZr2(1,1)* GYZ1(3)
     $ + 18* HZr2(1,1)* GYZ1(0)
     $ + 18* HZr3(1,1,0)

     $ + 102* GYZ2(2,2)
     $ + 18* GYZ3(2,2,0)
     $ - 54* GYZ3(2,3,2)
     $ - 230* GYZ1(2)

     $ + 18* GYZ3(2,0,2)
     $ + 27* GYZ2(2,0)
     $ - 36* GYZ3(3,2,2)
     $ - 66* GYZ2(3,2)

     $ - 18* GYZ3(3,2,0)
     $ + 108* GYZ3(3,3,2)
     $ - 18* GYZ3(3,0,2)
     $ + 18* GYZ3(0,2,2)

     $ + 27* GYZ2(0,2)
     $ - 18* GYZ3(0,3,2)
     $)
     $ + (1.0d0)/(9*(y+z))*(
     $ - (3*pi**2)/(2.0d0)* HZr1(1)
     $ + (3*pi**2)/(2.0d0)* GYZ1(2)
     $ - 170

     $ - 18* HZr1(0)* GYZ1(2)
     $ + 9* HZr1(0)* GYZ2(2,0)
     $ - 72* HZr2(0,1)
     $ + 9* HZr3(0,1,0)
     $ - 123* HZr1(1)

     $ - 90* HZr1(1)* GYZ1(3)
     $ + 18* HZr1(1)* GYZ1(0)
     $ + 18* HZr2(1,0)
     $ + 9* HZr2(1,0)* GYZ1(2)
     $ - 9* HZr2(1,0)* GYZ1(0)

     $ - 9* HZr3(1,1,0)
     $ + 123* GYZ1(2)
     $ - 18* GYZ2(2,0)
     $ - 9* GYZ3(2,1,0)
     $ + 90* GYZ2(3,2)

     $ - 18* GYZ2(0,2)
     $ - 9* GYZ3(0,1,0)
     $)
     $ + (T*pi**2)/(72.0d0)*(
     $ - 115
     $ - 24* HZr1(0)* GYZ1(2)
     $ - 12* HZr1(0)* GYZ1(1)
     $ - 12* HZr2(0,1)

     $ + 7* HZr1(1)
     $ - 48* HZr1(1)* GYZ1(2)
     $ + 36* HZr1(1)* GYZ1(0)
     $ - 12* HZr1(1)* GYZ1(1)
     $ - 12* HZr2(1,0)
     $ + 12* HZr2(1,1)

     $ + 48* GYZ2(2,2)
     $ + 19* GYZ1(2)
     $ - 24* GYZ2(2,0)
     $ - 36* GYZ2(0,2)
     $ + 48* GYZ2(0,1)

     $ + 12* GYZ2(1,2)
     $ - 26* GYZ1(1)
     $ + 36* GYZ2(1,0)
     $ - 48* GYZ2(1,1)
     $)
     $ + (T)/(54.0d0)*(
     $ -(3*pi**4)/(8.0d0)
     $ + (15251.0d0)/(12.0d0)
     $ - (357.0d0)/(2.0d0)*z3 
 
     $ + 108*z3* HZr1(1)
     $ - 270*z3* GYZ1(2)
     $ + 162*z3* GYZ1(1)
     $ - 360* HZr1(0)
     $ - 198* HZr1(0)* GYZ2(2,2)
 
     $ + 108* HZr1(0)* GYZ3(2,2,0)
     $ + 108* HZr1(0)* GYZ3(2,3,2)
     $ + 78* HZr1(0)* GYZ1(2)
 
     $ + 54* HZr1(0)* GYZ3(2,0,2)
     $ - 180* HZr1(0)* GYZ2(2,0)
     $ - 108* HZr1(0)* GYZ3(2,0,0)
 
     $ + 54* HZr1(0)* GYZ3(2,1,0)
     $ + 108* HZr1(0)* GYZ3(3,2,2)
     $ + 297* HZr1(0)* GYZ2(3,2)
 
     $ - 108* HZr1(0)* GYZ3(3,3,2)
     $ - 108* HZr1(0)* GYZ3(0,2,2)
     $ - 180* HZr1(0)* GYZ2(0,2)
 
     $ + 216* HZr1(0)* GYZ3(0,3,2)
     $ - 216* HZr1(0)* GYZ1(0)
     $ + 54* HZr1(0)* GYZ3(0,1,0)
     $ - 54* HZr1(0)* GYZ3(1,2,0)
 
     $ - 108* HZr1(0)* GYZ3(1,0,2)
     $ + 9* HZr1(0)* GYZ2(1,0)
     $ + 108* HZr1(0)* GYZ3(1,0,0)
     $ - 378* HZr2(0,0)* GYZ1(2)
 
     $ - 108* HZr2(0,0)* GYZ2(2,0)
     $ - 108* HZr2(0,0)* GYZ2(0,2)
     $ + 18* HZr3(0,0,1)
     $ + 216* HZr3(0,0,1)* GYZ1(2)
 
     $ - 432* HZr3(0,0,1)* GYZ1(3)
     $ - 108* HZr3(0,0,1)* GYZ1(0)
     $ + 54* HZr3(0,0,1)* GYZ1(1)
     $ + 216* HZr4(0,0,1,0)
 
     $ + 348* HZr2(0,1)
     $ + 324* HZr2(0,1)* GYZ2(2,3)
     $ - 279* HZr2(0,1)* GYZ1(2)
     $ - 162* HZr2(0,1)* GYZ2(2,0)
 
     $ + 108* HZr2(0,1)* GYZ2(3,2)
     $ - 540* HZr2(0,1)* GYZ2(3,3)
     $ - 63* HZr2(0,1)* GYZ1(3)
 
     $ + 108* HZr2(0,1)* GYZ2(3,0)
     $ + 54* HZr2(0,1)* GYZ2(0,2)
     $ - 108* HZr2(0,1)* GYZ2(0,3)
     $ - 117* HZr2(0,1)* GYZ1(0)
 
     $ + 108* HZr2(0,1)* GYZ2(0,0)
     $ - 54* HZr2(0,1)* GYZ2(1,2)
     $ + 54* HZr2(0,1)* GYZ2(1,0)
     $ - 252* HZr3(0,1,0)
 
     $ - 216* HZr3(0,1,0)* GYZ1(2)
     $ + 108* HZr3(0,1,0)* GYZ1(3)
     $ - 54* HZr3(0,1,0)* GYZ1(0)
     $ - 54* HZr3(0,1,0)* GYZ1(1)
 
     $ + 198* HZr3(0,1,1)
     $ - 108* HZr3(0,1,1)* GYZ1(3)
     $ + 108* HZr4(0,1,1,0)
     $ + 17* HZr1(1)
 
     $ + 432* HZr1(1)* GYZ3(2,3,3)
     $ - 477* HZr1(1)* GYZ2(2,3)
     $ - 108* HZr1(1)* GYZ3(2,3,0)
 
     $ + 297* HZr1(1)* GYZ1(2)
     $ - 108* HZr1(1)* GYZ3(2,0,3)
     $ + 198* HZr1(1)* GYZ2(2,0)
 
     $ + 162* HZr1(1)* GYZ3(2,1,0)
     $ + 216* HZr1(1)* GYZ3(3,2,3)
     $ - 396* HZr1(1)* GYZ2(3,2)
 
     $ - 108* HZr1(1)* GYZ3(3,2,0)
     $ + 216* HZr1(1)* GYZ3(3,3,2)
     $ - 648* HZr1(1)* GYZ3(3,3,3)
 
     $ + 234* HZr1(1)* GYZ2(3,3)
     $ + 108* HZr1(1)* GYZ3(3,3,0)
     $ + 426* HZr1(1)* GYZ1(3)
 
     $ - 108* HZr1(1)* GYZ3(3,0,2)
     $ + 108* HZr1(1)* GYZ3(3,0,3)
     $ - 297* HZr1(1)* GYZ2(3,0)
 
     $ - 54* HZr1(1)* GYZ3(0,2,3)
     $ + 198* HZr1(1)* GYZ2(0,2)
     $ - 108* HZr1(1)* GYZ3(0,3,2)
 
     $ + 108* HZr1(1)* GYZ3(0,3,3)
     $ - 297* HZr1(1)* GYZ2(0,3)
     $ + 216* HZr1(1)* GYZ3(0,3,0)
     $ - 78* HZr1(1)* GYZ1(0)
 
     $ + 108* HZr1(1)* GYZ3(0,0,3)
     $ + 378* HZr1(1)* GYZ2(0,0)
     $ - 54* HZr1(1)* GYZ3(0,1,0)
     $ - 54* HZr1(1)* GYZ3(1,2,3)
 
     $ - 54* HZr1(1)* GYZ3(1,0,3)
     $ - 81* HZr1(1)* GYZ2(1,0)
     $ - 108* HZr1(1)* GYZ3(1,0,0)
     $ - 81* HZr2(1,0)
 
     $ + 216* HZr2(1,0)* GYZ2(2,2)
     $ - 108* HZr2(1,0)* GYZ2(2,3)
     $ + 117* HZr2(1,0)* GYZ1(2)
 
     $ - 162* HZr2(1,0)* GYZ2(2,0)
     $ - 108* HZr2(1,0)* GYZ2(3,2)
     $ + 108* HZr2(1,0)* GYZ2(3,3)
 
     $ - 297* HZr2(1,0)* GYZ1(3)
     $ - 216* HZr2(1,0)* GYZ2(0,3)
     $ + 171* HZr2(1,0)* GYZ1(0)
 
     $ + 108* HZr2(1,0)* GYZ2(0,0)
     $ + 54* HZr2(1,0)* GYZ2(1,2)
     $ + 54* HZr2(1,0)* GYZ2(1,0)
     $ - 108* HZr3(1,0,0)* GYZ1(2)
 
     $ - 54* HZr4(1,0,0,1)
     $ + 360* HZr3(1,0,1)
     $ - 162* HZr3(1,0,1)* GYZ1(2)
     $ - 108* HZr3(1,0,1)* GYZ1(3)
 
     $ + 54* HZr3(1,0,1)* GYZ1(0)
     $ + 54* HZr3(1,0,1)* GYZ1(1)
     $ + 108* HZr4(1,0,1,0)
     $ - 297* HZr2(1,1)
 
     $ - 216* HZr2(1,1)* GYZ2(3,3)
     $ + 396* HZr2(1,1)* GYZ1(3)
     $ + 108* HZr2(1,1)* GYZ2(3,0)
     $ + 108* HZr2(1,1)* GYZ2(0,3)
 
     $ - 198* HZr2(1,1)* GYZ1(0)
     $ + 81* HZr3(1,1,0)
     $ - 378* HZr3(1,1,0)* GYZ1(2)
     $ + 108* HZr3(1,1,0)* GYZ1(3)
 
     $ + 108* HZr3(1,1,0)* GYZ1(0)
     $ - 54* HZr3(1,1,0)* GYZ1(1)
     $ + 108* HZr4(1,1,0,1)
     $ + 162* HZr4(1,1,1,0)
 
     $ - 297* GYZ2(2,2)
     $ - 198* GYZ3(2,2,0)
     $ - 216* GYZ4(2,2,1,0)
     $ + 477* GYZ3(2,3,2)
 
     $ + 108* GYZ4(2,3,2,0)
     $ - 432* GYZ4(2,3,3,2)
     $ + 108* GYZ4(2,3,0,2)
 
     $ - 17* GYZ1(2)
     $ - 198* GYZ3(2,0,2)
     $ + 108* GYZ4(2,0,3,2)
     $ + 78* GYZ2(2,0)
 
     $ - 378* GYZ3(2,0,0)
     $ + 108* GYZ4(2,0,1,0)
     $ - 162* GYZ4(2,1,2,0)
     $ - 162* GYZ4(2,1,0,2)
 
     $ + 81* GYZ3(2,1,0)
     $ + 108* GYZ4(2,1,0,0)
     $ + 396* GYZ3(3,2,2)
     $ + 108* GYZ4(3,2,2,0)
 
     $ - 216* GYZ4(3,2,3,2)
     $ - 426* GYZ2(3,2)
     $ + 108* GYZ4(3,2,0,2)
     $ + 297* GYZ3(3,2,0)
 
     $ - 216* GYZ4(3,3,2,2)
     $ - 234* GYZ3(3,3,2)
     $ - 108* GYZ4(3,3,2,0)
 
     $ + 648* GYZ4(3,3,3,2)
     $ - 108* GYZ4(3,3,0,2)
     $ + 108* GYZ4(3,0,2,2)
 
     $ + 297* GYZ3(3,0,2)
     $ - 108* GYZ4(3,0,3,2)
     $ - 198* GYZ3(0,2,2)
     $ + 54* GYZ4(0,2,3,2)
 
     $ + 78* GYZ2(0,2)
     $ - 378* GYZ3(0,2,0)
     $ + 162* GYZ4(0,2,1,0)
     $ + 108* GYZ4(0,3,2,2)
 
     $ + 297* GYZ3(0,3,2)
     $ - 216* GYZ4(0,3,2,0)
     $ - 108* GYZ4(0,3,3,2)
     $ - 216* GYZ4(0,3,0,2)
 
     $ - 360* GYZ1(0)
     $ - 378* GYZ3(0,0,2)
     $ - 108* GYZ4(0,0,3,2)
     $ + 108* GYZ4(0,0,1,0)
     $ + 54* GYZ4(0,1,2,0)
 
     $ + 54* GYZ4(0,1,0,2)
     $ + 333* GYZ3(0,1,0)
     $ - 216* GYZ4(0,1,1,0)
     $ + 54* GYZ4(1,2,3,2)
 
     $ + 81* GYZ3(1,2,0)
     $ + 108* GYZ4(1,2,0,0)
     $ + 81* GYZ3(1,0,2)
     $ + 108* GYZ4(1,0,2,0)
 
     $ + 54* GYZ4(1,0,3,2)
     $ + 3* GYZ2(1,0)
     $ + 108* GYZ4(1,0,0,2)
     $ + 378* GYZ3(1,0,0)
     $ - 216* GYZ4(1,0,1,0)
 
     $ + 117* GYZ3(1,1,0)
     $ - 216* GYZ4(1,1,0,0)
     $ + 216* GYZ4(1,1,1,0)
     $)
     $ + (pi**2)/(18.0d0)*(
     $   11
     $ + 9* HZr1(0)
     $ - 24* HZr1(0)* GYZ1(2)
 
     $ - 6* HZr2(0,1)
     $ + 8* HZr1(1)
     $ - 36* HZr1(1)* GYZ1(2)
     $ + 24* HZr1(1)* GYZ1(0)
     $ + 6* HZr1(1)* GYZ1(1)
     $ + 24* HZr2(1,1)
 
     $ + 24* GYZ2(2,2)
     $ + 2* GYZ1(2)
     $ - 24* GYZ2(2,0)
     $ + 12* GYZ2(2,1)
     $ - 6* GYZ2(0,2)
     $ + 9* GYZ1(0)
 
     $ + 12* GYZ2(0,1)
     $ - 6* GYZ2(1,2)
     $ - 10* GYZ1(1)
     $ + 24* GYZ2(1,0)
     $ - 6* GYZ2(1,1)
     $)
     $ + (1.0d0)/(18.0d0)*(
     $   288*z3
     $ + 180*z3* HZr1(1)
 
     $ - 360*z3* GYZ1(2)
     $ + 180*z3* GYZ1(1)
     $ - 188* HZr1(0)
     $ - 150* HZr1(0)* GYZ2(2,2)
 
     $ + 72* HZr1(0)* GYZ3(2,2,0)
     $ + 72* HZr1(0)* GYZ3(2,3,2)
     $ + 295* HZr1(0)* GYZ1(2)
 
     $ - 300* HZr1(0)* GYZ2(2,0)
     $ - 36* HZr1(0)* GYZ3(2,0,0)
     $ + 72* HZr1(0)* GYZ3(3,2,2)
 
     $ + 216* HZr1(0)* GYZ2(3,2)
     $ - 72* HZr1(0)* GYZ3(3,3,2)
     $ - 36* HZr1(0)* GYZ3(0,2,2)
 
     $ - 156* HZr1(0)* GYZ2(0,2)
     $ + 36* HZr1(0)* GYZ3(0,2,0)
     $ + 144* HZr1(0)* GYZ3(0,3,2)
     $ + 78* HZr1(0)* GYZ1(0)
 
     $ - 72* HZr1(0)* GYZ3(0,0,2)
     $ + 18* HZr1(0)* GYZ2(0,0)
     $ + 36* HZr1(0)* GYZ3(0,1,0)
     $ - 72* HZr1(0)* GYZ3(1,2,0)
 
     $ - 36* HZr1(0)* GYZ3(1,0,2)
     $ + 132* HZr1(0)* GYZ2(1,0)
     $ + 36* HZr1(0)* GYZ3(1,0,0)
     $ + 36* HZr1(0)* GYZ3(1,1,0)
 
     $ + 36* HZr2(0,0)
     $ - 36* HZr2(0,0)* GYZ2(2,2)
     $ - 108* HZr2(0,0)* GYZ1(2)
     $ - 36* HZr2(0,0)* GYZ2(2,0)
 
     $ - 36* HZr2(0,0)* GYZ2(0,2)
     $ + 18* HZr2(0,0)* GYZ1(0)
     $ - 108* HZr4(0,0,0,1)
     $ + 84* HZr3(0,0,1)
 
     $ + 180* HZr3(0,0,1)* GYZ1(2)
     $ - 288* HZr3(0,0,1)* GYZ1(3)
     $ + 36* HZr3(0,0,1)* GYZ1(1)
     $ + 108* HZr4(0,0,1,0)
 
     $ - 36* HZr4(0,0,1,1)
     $ + 289* HZr2(0,1)
     $ + 216* HZr2(0,1)* GYZ2(2,3)
     $ - 222* HZr2(0,1)* GYZ1(2)
 
     $ - 72* HZr2(0,1)* GYZ2(2,0)
     $ + 72* HZr2(0,1)* GYZ2(3,2)
     $ - 360* HZr2(0,1)* GYZ2(3,3)
 
     $ + 192* HZr2(0,1)* GYZ1(3)
     $ + 72* HZr2(0,1)* GYZ2(3,0)
     $ - 72* HZr2(0,1)* GYZ2(0,3)
     $ - 60* HZr2(0,1)* GYZ1(0)
 
     $ + 36* HZr2(0,1)* GYZ2(0,0)
     $ - 36* HZr2(0,1)* GYZ2(1,2)
     $ - 114* HZr3(0,1,0)
     $ - 144* HZr3(0,1,0)* GYZ1(2)
 
     $ + 72* HZr3(0,1,0)* GYZ1(3)
     $ - 72* HZr3(0,1,0)* GYZ1(0)
     $ + 36* HZr3(0,1,0)* GYZ1(1)
     $ - 72* HZr4(0,1,0,1)
 
     $ + 150* HZr3(0,1,1)
     $ - 72* HZr3(0,1,1)* GYZ1(3)
     $ + 36* HZr3(0,1,1)* GYZ1(0)
     $ + 36* HZr4(0,1,1,0)
     $ - 376* HZr1(1)
 
     $ + 288* HZr1(1)* GYZ3(2,3,3)
     $ - 372* HZr1(1)* GYZ2(2,3)
     $ - 72* HZr1(1)* GYZ3(2,3,0)
 
     $ + 204* HZr1(1)* GYZ1(2)
     $ - 72* HZr1(1)* GYZ3(2,0,3)
     $ + 150* HZr1(1)* GYZ2(2,0)
 
     $ + 36* HZr1(1)* GYZ3(2,0,0)
     $ + 36* HZr1(1)* GYZ3(2,1,0)
     $ + 144* HZr1(1)* GYZ3(3,2,3)
 
     $ - 300* HZr1(1)* GYZ2(3,2)
     $ - 72* HZr1(1)* GYZ3(3,2,0)
     $ + 144* HZr1(1)* GYZ3(3,3,2)
 
     $ - 432* HZr1(1)* GYZ3(3,3,3)
     $ + 408* HZr1(1)* GYZ2(3,3)
     $ + 72* HZr1(1)* GYZ3(3,3,0)
 
     $ + 584* HZr1(1)* GYZ1(3)
     $ - 72* HZr1(1)* GYZ3(3,0,2)
     $ + 72* HZr1(1)* GYZ3(3,0,3)
     $ - 216* HZr1(1)* GYZ2(3,0)
 
     $ - 36* HZr1(1)* GYZ3(0,2,3)
     $ + 150* HZr1(1)* GYZ2(0,2)
     $ + 36* HZr1(1)* GYZ3(0,2,0)
 
     $ - 72* HZr1(1)* GYZ3(0,3,2)
     $ + 72* HZr1(1)* GYZ3(0,3,3)
     $ - 216* HZr1(1)* GYZ2(0,3)
 
     $ + 144* HZr1(1)* GYZ3(0,3,0)
     $ - 295* HZr1(1)* GYZ1(0)
     $ + 36* HZr1(1)* GYZ3(0,0,2)
     $ - 36* HZr1(1)* GYZ3(0,0,3)
 
     $ + 108* HZr1(1)* GYZ2(0,0)
     $ - 36* HZr1(1)* GYZ3(1,2,3)
     $ - 36* HZr1(1)* GYZ3(1,0,3)
     $ - 72* HZr1(1)* GYZ3(1,0,0)
 
     $ + 15* HZr2(1,0)
     $ + 108* HZr2(1,0)* GYZ2(2,2)
     $ - 72* HZr2(1,0)* GYZ2(2,3)
     $ + 6* HZr2(1,0)* GYZ1(2)
 
     $ - 72* HZr2(1,0)* GYZ2(2,0)
     $ - 72* HZr2(1,0)* GYZ2(3,2)
     $ + 72* HZr2(1,0)* GYZ2(3,3)
 
     $ - 216* HZr2(1,0)* GYZ1(3)
     $ + 108* HZr2(1,0)* GYZ2(0,2)
     $ - 144* HZr2(1,0)* GYZ2(0,3)
     $ + 168* HZr2(1,0)* GYZ1(0)
 
     $ + 36* HZr2(1,0)* GYZ2(0,0)
     $ - 36* HZr2(1,0)* GYZ2(1,2)
     $ + 72* HZr2(1,0)* GYZ2(1,0)
     $ + 18* HZr3(1,0,0)
 
     $ - 36* HZr3(1,0,0)* GYZ1(2)
     $ - 72* HZr4(1,0,0,1)
     $ + 222* HZr3(1,0,1)
     $ - 36* HZr3(1,0,1)* GYZ1(2)
 
     $ - 72* HZr3(1,0,1)* GYZ1(3)
     $ + 36* HZr3(1,0,1)* GYZ1(0)
     $ + 36* HZr3(1,0,1)* GYZ1(1)
     $ + 72* HZr4(1,0,1,0)
 
     $ - 204* HZr2(1,1)
     $ - 144* HZr2(1,1)* GYZ2(3,3)
     $ + 300* HZr2(1,1)* GYZ1(3)
     $ + 72* HZr2(1,1)* GYZ2(3,0)
 
     $ + 72* HZr2(1,1)* GYZ2(0,3)
     $ - 150* HZr2(1,1)* GYZ1(0)
     $ - 36* HZr2(1,1)* GYZ2(0,0)
     $ + 54* HZr3(1,1,0)
 
     $ - 216* HZr3(1,1,0)* GYZ1(2)
     $ + 72* HZr3(1,1,0)* GYZ1(3)
     $ + 36* HZr3(1,1,0)* GYZ1(0)
     $ + 36* HZr3(1,1,0)* GYZ1(1)
 
     $ + 108* HZr4(1,1,1,0)
     $ - 204* GYZ2(2,2)
     $ - 150* GYZ3(2,2,0)
     $ - 36* GYZ4(2,2,0,0)
 
     $ - 108* GYZ4(2,2,1,0)
     $ + 372* GYZ3(2,3,2)
     $ + 72* GYZ4(2,3,2,0)
 
     $ - 288* GYZ4(2,3,3,2)
     $ + 72* GYZ4(2,3,0,2)
     $ + 376* GYZ1(2)
     $ - 150* GYZ3(2,0,2)
 
     $ - 36* GYZ4(2,0,2,0)
     $ + 72* GYZ4(2,0,3,2)
     $ + 295* GYZ2(2,0)
     $ - 36* GYZ4(2,0,0,2)
 
     $ - 108* GYZ3(2,0,0)
     $ + 108* GYZ4(2,0,1,0)
     $ - 36* GYZ4(2,1,2,0)
     $ - 36* GYZ4(2,1,0,2)
 
     $ + 144* GYZ3(2,1,0)
     $ + 72* GYZ4(2,1,0,0)
     $ - 72* GYZ4(2,1,1,0)
     $ + 300* GYZ3(3,2,2)
 
     $ + 72* GYZ4(3,2,2,0)
     $ - 144* GYZ4(3,2,3,2)
     $ - 584* GYZ2(3,2)
 
     $ + 72* GYZ4(3,2,0,2)
     $ + 216* GYZ3(3,2,0)
     $ - 144* GYZ4(3,3,2,2)
 
     $ - 408* GYZ3(3,3,2)
     $ - 72* GYZ4(3,3,2,0)
     $ + 432* GYZ4(3,3,3,2)
 
     $ - 72* GYZ4(3,3,0,2)
     $ + 72* GYZ4(3,0,2,2)
     $ + 216* GYZ3(3,0,2)
 
     $ - 72* GYZ4(3,0,3,2)
     $ - 150* GYZ3(0,2,2)
     $ - 36* GYZ4(0,2,2,0)
 
     $ + 36* GYZ4(0,2,3,2)
     $ + 295* GYZ2(0,2)
     $ - 36* GYZ4(0,2,0,2)
     $ - 108* GYZ3(0,2,0)
 
     $ + 72* GYZ4(0,2,1,0)
     $ + 72* GYZ4(0,3,2,2)
     $ + 216* GYZ3(0,3,2)
     $ - 144* GYZ4(0,3,2,0)
 
     $ - 72* GYZ4(0,3,3,2)
     $ - 144* GYZ4(0,3,0,2)
     $ - 188* GYZ1(0)
     $ - 36* GYZ4(0,0,2,2)
 
     $ - 108* GYZ3(0,0,2)
     $ + 36* GYZ4(0,0,3,2)
     $ + 36* GYZ2(0,0)
     $ + 108* GYZ4(0,0,1,0)
     $ + 6* GYZ3(0,1,0)
 
     $ - 72* GYZ4(0,1,1,0)
     $ + 36* GYZ4(1,2,3,2)
     $ + 72* GYZ4(1,2,0,0)
     $ + 72* GYZ4(1,2,1,0)
 
     $ + 72* GYZ4(1,0,2,0)
     $ + 36* GYZ4(1,0,3,2)
     $ - 310* GYZ2(1,0)
     $ + 72* GYZ4(1,0,0,2)
 
     $ + 90* GYZ3(1,0,0)
     $ - 144* GYZ4(1,0,1,0)
     $ + 60* GYZ3(1,1,0)
     $ - 108* GYZ4(1,1,0,0)
     $ + 36* GYZ4(1,1,1,0)
     $)


      else if(i.eq.3) then
      fin20=
     $  + (z)/(y**2)*(
     $   6* HZr1(0)* GYZ1(2)
     $ + 6* HZr1(1)* GYZ1(3)
     $ - 6* GYZ2(3,2)
     $ )
     $ + (z**2)/(y**2)*(
     $ - 2* HZr1(0)* GYZ1(2)
 
     $ - 2* HZr1(1)* GYZ1(3)
     $ + 2* GYZ2(3,2)
     $ )
     $ + (1.0d0)/(y**2)*(
     $ - 4* HZr1(0)* GYZ1(2)
     $ - 4* HZr1(1)* GYZ1(3)
 
     $ + 4* GYZ2(3,2)
     $ )
     $ + (z*pi**2)/(6*y)*(
     $ -1
     $ + 4* HZr1(0)
     $ + 10* HZr1(0)* GYZ1(2)
     $ + HZr1(1)
     $ + 6* HZr1(1)* GYZ1(2)
 
     $ + 10* HZr1(1)* GYZ1(3)
     $ - 16* GYZ2(2,2)
     $ + GYZ1(2)
     $ + 2* GYZ2(2,0)
     $ + 4* GYZ2(2,1)
 
     $ - 10* GYZ2(3,2)
     $ )
     $ + (z)/(4*y)*(
     $ - 36*z3
     $ - 16*z3* GYZ1(2)
     $ - 4* HZr1(0)
     $ + 6* HZr1(0)* GYZ2(2,2)
 
     $ - 40* HZr1(0)* GYZ3(2,2,0)
     $ + 9* HZr1(0)* GYZ1(2)
     $ - 8* HZr1(0)* GYZ3(2,0,2)
 
     $ - 24* HZr1(0)* GYZ2(2,0)
     $ - 16* HZr1(0)* GYZ3(3,2,2)
     $ + 46* HZr1(0)* GYZ2(3,2)
 
     $ - 40* HZr1(0)* GYZ3(3,2,0)
     $ - 10* HZr1(0)* GYZ2(0,2)
     $ + 12* HZr1(0)* GYZ2(1,0)
 
     $ + 16* HZr2(0,0)* GYZ2(2,2)
     $ + 40* HZr2(0,0)* GYZ2(2,0)
     $ - 4* HZr3(0,0,1)
     $ + 17* HZr2(0,1)
 
     $ + 24* HZr2(0,1)* GYZ2(2,2)
     $ - 14* HZr2(0,1)* GYZ1(2)
     $ - 8* HZr2(0,1)* GYZ2(2,0)
 
     $ - 8* HZr2(0,1)* GYZ2(3,2)
     $ - 24* HZr2(0,1)* GYZ2(3,3)
     $ + 14* HZr2(0,1)* GYZ1(3)
     $ - 8* HZr3(0,1,0)
 
     $ + 32* HZr3(0,1,0)* GYZ1(2)
     $ - 40* HZr3(0,1,0)* GYZ1(3)
     $ - 6* HZr3(0,1,1)
     $ + 8* HZr3(0,1,1)* GYZ1(3)
 
     $ + 24* HZr1(1)* GYZ3(2,2,3)
     $ - 8* HZr1(1)* GYZ2(2,3)
     $ - 16* HZr1(1)* GYZ3(2,0,3)
 
     $ - 6* HZr1(1)* GYZ2(2,0)
     $ - 24* HZr1(1)* GYZ3(3,2,3)
     $ + 12* HZr1(1)* GYZ2(3,2)
 
     $ + 8* HZr1(1)* GYZ3(3,2,0)
     $ - 24* HZr1(1)* GYZ3(3,3,2)
     $ - 24* HZr1(1)* GYZ3(3,3,3)
 
     $ + 60* HZr1(1)* GYZ2(3,3)
     $ + 24* HZr1(1)* GYZ3(3,3,0)
     $ + 26* HZr1(1)* GYZ1(3)
     $ + 8* HZr1(1)* GYZ3(3,0,2)
 
     $ + 18* HZr1(1)* GYZ2(3,0)
     $ - 6* HZr1(1)* GYZ2(0,2)
     $ - 10* HZr1(1)* GYZ2(0,3)
     $ - 17* HZr1(1)* GYZ1(0)
 
     $ + 24* HZr1(1)* GYZ2(1,0)
     $ - 13* HZr2(1,0)
     $ - 40* HZr2(1,0)* GYZ2(2,2)
     $ + 4* HZr2(1,0)* GYZ1(2)
 
     $ - 24* HZr2(1,0)* GYZ2(3,2)
     $ - 46* HZr2(1,0)* GYZ1(3)
     $ + 40* HZr2(1,0)* GYZ2(3,0)
     $ + 40* HZr3(1,0,0)* GYZ1(2)
 
     $ - 10* HZr3(1,0,1)
     $ - 16* HZr3(1,0,1)* GYZ1(2)
     $ + 8* HZr3(1,0,1)* GYZ1(3)
     $ + 24* HZr2(1,1)* GYZ2(3,3)
 
     $ - 12* HZr2(1,1)* GYZ1(3)
     $ - 8* HZr2(1,1)* GYZ2(3,0)
     $ + 6* HZr2(1,1)* GYZ1(0)
     $ - 20* HZr3(1,1,0)
 
     $ + 8* HZr3(1,1,0)* GYZ1(2)
     $ + 24* HZr3(1,1,0)* GYZ1(3)
     $ - 24* GYZ4(2,2,3,2)
 
     $ + 6* GYZ3(2,2,0)
     $ + 48* GYZ4(2,2,1,0)
     $ + 8* GYZ3(2,3,2)
     $ + 6* GYZ3(2,0,2)
 
     $ + 16* GYZ4(2,0,3,2)
     $ + 17* GYZ2(2,0)
     $ - 48* GYZ4(2,0,1,0)
     $ + 18* GYZ3(2,1,0)
 
     $ - 16* GYZ4(2,1,1,0)
     $ - 12* GYZ3(3,2,2)
     $ - 8* GYZ4(3,2,2,0)
 
     $ + 24* GYZ4(3,2,3,2)
     $ - 26* GYZ2(3,2)
     $ - 8* GYZ4(3,2,0,2)
     $ - 18* GYZ3(3,2,0)
 
     $ + 64* GYZ4(3,2,1,0)
     $ + 24* GYZ4(3,3,2,2)
     $ - 60* GYZ3(3,3,2)
 
     $ - 24* GYZ4(3,3,2,0)
     $ + 24* GYZ4(3,3,3,2)
     $ - 24* GYZ4(3,3,0,2)
 
     $ - 8* GYZ4(3,0,2,2)
     $ - 18* GYZ3(3,0,2)
     $ + 64* GYZ4(3,0,1,0)
     $ + 6* GYZ3(0,2,2)
 
     $ + 17* GYZ2(0,2)
     $ + 10* GYZ3(0,3,2)
     $ + 18* GYZ3(0,1,0)
     $ - 24* GYZ3(1,2,0)
     $ - 24* GYZ3(1,0,2)
 
     $ - 30* GYZ2(1,0)
     $ )
     $ + (z**2)/(y)*(
     $ - 2* HZr1(0)* GYZ1(2)
     $ - 2* HZr1(1)* GYZ1(3)
     $ + 2* GYZ2(3,2)
     $ )
 
     $ + (1.0d0)/(2*y*(y+z))*(
     $ - 6* HZr1(0)* GYZ2(2,2)
     $ - 17* HZr1(0)* GYZ1(2)
     $ + 18* HZr1(0)* GYZ2(3,2)
 
     $ - 18* HZr1(0)* GYZ2(0,2)
     $ + 17* HZr2(0,1)
     $ + 6* HZr2(0,1)* GYZ1(2)
     $ - 18* HZr2(0,1)* GYZ1(3)
     $ - 6* HZr3(0,1,1)
 
     $ - 6* HZr1(1)* GYZ2(2,0)
     $ + 18* HZr1(1)* GYZ2(3,0)
     $ - 6* HZr1(1)* GYZ2(0,2)
     $ - 18* HZr1(1)* GYZ2(0,3)
 
     $ - 17* HZr1(1)* GYZ1(0)
     $ - 17* HZr2(1,0)
     $ + 24* HZr2(1,0)* GYZ1(2)
     $ - 18* HZr2(1,0)* GYZ1(3)
     $ - 6* HZr3(1,0,1)
 
     $ + 6* HZr2(1,1)* GYZ1(0)
     $ - 24* HZr3(1,1,0)
     $ + 6* GYZ3(2,2,0)
     $ + 6* GYZ3(2,0,2)
     $ + 17* GYZ2(2,0)
 
     $ + 18* GYZ3(2,1,0)
     $ - 18* GYZ3(3,2,0)
     $ - 18* GYZ3(3,0,2)
     $ + 6* GYZ3(0,2,2)
 
     $ + 17* GYZ2(0,2)
     $ + 18* GYZ3(0,3,2)
     $ + 18* GYZ3(0,1,0)
     $ - 34* GYZ2(1,0)
     $ )
     $ + (pi**2)/(3*y)*(
     $ - 4* HZr1(0)* GYZ1(2)
 
     $ + HZr1(1)
     $ - 6* HZr1(1)* GYZ1(2)
     $ - 4* HZr1(1)* GYZ1(3)
     $ + 10* GYZ2(2,2)
     $ - 6* GYZ1(2)
 
     $ - 2* GYZ2(2,0)
     $ - 4* GYZ2(2,1)
     $ + 4* GYZ2(3,2)
     $ + GYZ1(1)
     $ )
     $ + (1.0d0)/(4*y)*(
     $ - 8*z3
     $ + 80*z3* GYZ1(2)
     $ + 19
 
     $ + 8* HZr1(0)
     $ - 12* HZr1(0)* GYZ2(2,2)
     $ + 32* HZr1(0)* GYZ3(2,2,0)
     $ + 26* HZr1(0)* GYZ1(2)
 
     $ + 16* HZr1(0)* GYZ3(2,0,2)
     $ + 40* HZr1(0)* GYZ2(2,0)
     $ + 32* HZr1(0)* GYZ3(3,2,2)
 
     $ - 108* HZr1(0)* GYZ2(3,2)
     $ + 32* HZr1(0)* GYZ3(3,2,0)
     $ + 68* HZr1(0)* GYZ2(0,2)
     $ - 16* HZr1(0)* GYZ2(1,0)
 
     $ - 32* HZr2(0,0)* GYZ2(2,2)
     $ - 32* HZr2(0,0)* GYZ2(2,0)
     $ - 8* HZr3(0,0,1)
     $ - 26* HZr2(0,1)
 
     $ - 48* HZr2(0,1)* GYZ2(2,2)
     $ + 36* HZr2(0,1)* GYZ1(2)
     $ + 16* HZr2(0,1)* GYZ2(2,0)
 
     $ + 16* HZr2(0,1)* GYZ2(3,2)
     $ + 48* HZr2(0,1)* GYZ2(3,3)
     $ - 28* HZr2(0,1)* GYZ1(3)
     $ + 16* HZr3(0,1,0)
 
     $ - 16* HZr3(0,1,0)* GYZ1(2)
     $ + 32* HZr3(0,1,0)* GYZ1(3)
     $ + 12* HZr3(0,1,1)
     $ - 16* HZr3(0,1,1)* GYZ1(3)
     $ - 2* HZr1(1)
 
     $ - 48* HZr1(1)* GYZ3(2,2,3)
     $ + 24* HZr1(1)* GYZ2(2,3)
     $ + 32* HZr1(1)* GYZ3(2,0,3)
 
     $ + 12* HZr1(1)* GYZ2(2,0)
     $ + 48* HZr1(1)* GYZ3(3,2,3)
     $ - 24* HZr1(1)* GYZ2(3,2)
 
     $ - 16* HZr1(1)* GYZ3(3,2,0)
     $ + 48* HZr1(1)* GYZ3(3,3,2)
     $ + 48* HZr1(1)* GYZ3(3,3,3)
 
     $ - 136* HZr1(1)* GYZ2(3,3)
     $ - 48* HZr1(1)* GYZ3(3,3,0)
     $ - 16* HZr1(1)* GYZ3(3,0,2)
     $ - 20* HZr1(1)* GYZ2(3,0)
 
     $ + 12* HZr1(1)* GYZ2(0,2)
     $ + 68* HZr1(1)* GYZ2(0,3)
     $ + 30* HZr1(1)* GYZ1(0)
     $ - 24* HZr1(1)* GYZ2(1,0)
     $ + 26* HZr2(1,0)
 
     $ + 32* HZr2(1,0)* GYZ2(2,2)
     $ - 56* HZr2(1,0)* GYZ1(2)
     $ + 108* HZr2(1,0)* GYZ1(3)
     $ - 32* HZr2(1,0)* GYZ2(3,0)
 
     $ - 32* HZr3(1,0,0)* GYZ1(2)
     $ + 4* HZr3(1,0,1)
     $ + 32* HZr3(1,0,1)* GYZ1(2)
     $ - 16* HZr3(1,0,1)* GYZ1(3)
 
     $ - 48* HZr2(1,1)* GYZ2(3,3)
     $ + 24* HZr2(1,1)* GYZ1(3)
     $ + 16* HZr2(1,1)* GYZ2(3,0)
     $ - 12* HZr2(1,1)* GYZ1(0)
 
     $ + 56* HZr3(1,1,0)
     $ - 16* HZr3(1,1,0)* GYZ1(2)
     $ + 48* GYZ4(2,2,3,2)
     $ - 12* GYZ3(2,2,0)
 
     $ - 48* GYZ4(2,2,1,0)
     $ - 24* GYZ3(2,3,2)
     $ + 2* GYZ1(2)
     $ - 12* GYZ3(2,0,2)
 
     $ - 32* GYZ4(2,0,3,2)
     $ - 30* GYZ2(2,0)
     $ + 48* GYZ4(2,0,1,0)
     $ - 36* GYZ3(2,1,0)
 
     $ + 32* GYZ4(2,1,1,0)
     $ + 24* GYZ3(3,2,2)
     $ + 16* GYZ4(3,2,2,0)
 
     $ - 48* GYZ4(3,2,3,2)
     $ + 16* GYZ4(3,2,0,2)
     $ + 20* GYZ3(3,2,0)
     $ - 80* GYZ4(3,2,1,0)
 
     $ - 48* GYZ4(3,3,2,2)
     $ + 136* GYZ3(3,3,2)
     $ + 48* GYZ4(3,3,2,0)
 
     $ - 48* GYZ4(3,3,3,2)
     $ + 48* GYZ4(3,3,0,2)
     $ + 16* GYZ4(3,0,2,2)
     $ + 20* GYZ3(3,0,2)
 
     $ - 80* GYZ4(3,0,1,0)
     $ - 12* GYZ3(0,2,2)
     $ - 30* GYZ2(0,2)
     $ - 68* GYZ3(0,3,2)
     $ - 36* GYZ3(0,1,0)
 
     $ + 24* GYZ3(1,2,0)
     $ + 24* GYZ3(1,0,2)
     $ + 60* GYZ2(1,0)
     $ - 8* GYZ3(1,1,0)
     $ )
     $ + (z)/(4*(1-y)**2)*(
     $  5*pi**2
     $ + (14*pi**2)/(3.0d0)* HZr1(0)
     $ + (14*pi**2)/(3.0d0)* HZr1(1)
     $ - (14*pi**2)/(3.0d0)* GYZ1(2)
     $ + (2*pi**2)/(3.0d0)* GYZ1(0)
     $ - (4*pi**2)/(3.0d0)* GYZ1(1)
     $ - 92*z3
 
     $ + 28* HZr1(0)* GYZ2(0,2)
     $ + 12* HZr1(0)* GYZ1(0)
     $ + 4* HZr1(0)* GYZ2(0,0)
     $ - 32* HZr3(0,0,1)
     $ - 16* HZr2(0,1)
 
     $ + 32* HZr2(0,1)* GYZ1(2)
     $ + 4* HZr2(0,1)* GYZ1(0)
     $ + 28* HZr3(0,1,0)
     $ + 32* HZr1(1)* GYZ2(2,3)
 
     $ - 16* HZr1(1)* GYZ1(3)
     $ + 32* HZr1(1)* GYZ2(0,3)
     $ + 6* HZr1(1)* GYZ1(0)
     $ - 4* HZr1(1)* GYZ2(0,0)
     $ + 16* HZr2(1,0)
 
     $ - 28* HZr2(1,0)* GYZ1(2)
     $ - 32* HZr3(1,0,1)
     $ + 28* HZr3(1,1,0)
     $ - 32* GYZ3(2,3,2)
     $ - 6* GYZ2(2,0)
 
     $ + 4* GYZ3(2,0,0)
     $ - 4* GYZ3(2,1,0)
     $ + 16* GYZ2(3,2)
     $ - 6* GYZ2(0,2)
     $ + 4* GYZ3(0,2,0)
 
     $ - 32* GYZ3(0,3,2)
     $ + 23* GYZ1(0)
     $ + 4* GYZ3(0,0,2)
     $ - 10* GYZ2(0,0)
     $ - 8* GYZ3(0,1,0)
     $ - 14* GYZ2(1,0)
 
     $ - 8* GYZ3(1,0,0)
     $ + 8* GYZ3(1,1,0)
     $ )
     $ + (z)/(4*(1-y))*(
     $   (11*pi**2)/(3.0d0) 
     $ + 6*pi**2* HZr1(0)
     $ + 6*pi**2* HZr1(1)
     $ - 6*pi**2* GYZ1(2)
 
     $ + 2*pi**2* GYZ1(0)
     $ - 4*pi**2* GYZ1(1)
     $ - 132*z3
     $ + 17
     $ + 4* HZr1(0)
     $ + 28* HZr1(0)* GYZ1(2)
     $ + 36* HZr1(0)* GYZ2(0,2)
 
     $ + 12* HZr1(0)* GYZ2(0,0)
     $ - 48* HZr3(0,0,1)
     $ + 4* HZr2(0,1)
     $ + 48* HZr2(0,1)* GYZ1(2)
     $ + 12* HZr2(0,1)* GYZ1(0)
 
     $ + 36* HZr3(0,1,0)
     $ - 2* HZr1(1)
     $ + 48* HZr1(1)* GYZ2(2,3)
     $ + 32* HZr1(1)* GYZ1(3)
     $ + 48* HZr1(1)* GYZ2(0,3)
 
     $ - 22* HZr1(1)* GYZ1(0)
     $ - 12* HZr1(1)* GYZ2(0,0)
     $ - 36* HZr2(1,0)* GYZ1(2)
     $ - 48* HZr3(1,0,1)
     $ + 36* HZr3(1,1,0)
 
     $ - 48* GYZ3(2,3,2)
     $ + 2* GYZ1(2)
     $ + 22* GYZ2(2,0)
     $ + 12* GYZ3(2,0,0)
     $ - 12* GYZ3(2,1,0)
 
     $ - 32* GYZ2(3,2)
     $ + 22* GYZ2(0,2)
     $ + 12* GYZ3(0,2,0)
     $ - 48* GYZ3(0,3,2)
     $ + 39* GYZ1(0)
 
     $ + 12* GYZ3(0,0,2)
     $ - 22* GYZ2(0,0)
     $ - 24* GYZ3(0,1,0)
     $ - 26* GYZ2(1,0)
     $ - 24* GYZ3(1,0,0)
     $ + 24* GYZ3(1,1,0)
     $ )
 
     $ + (z)/((y+z)**3)*(
     $ - 6*pi**2* HZr1(1)
     $ + 6*pi**2* GYZ1(2)
     $ + 36* HZr1(0)* GYZ2(2,0)
     $ + 36* HZr3(0,1,0)
     $ + 36* HZr2(1,0)
 
     $ + 36* HZr2(1,0)* GYZ1(2)
     $ - 36* HZr2(1,0)* GYZ1(0)
     $ - 36* HZr3(1,1,0)
     $ - 36* GYZ3(2,1,0)
     $ - 36* GYZ3(0,1,0)
 
     $ + 36* GYZ2(1,0)
     $ )
     $ + (z)/((y+z)**2)*(
     $  6*pi**2
     $ + 4*pi**2* HZr1(1)
     $ - 4*pi**2* GYZ1(2)
     $ - 18* HZr1(0)
     $ - 24* HZr1(0)* GYZ2(2,0)
 
     $ + 36* HZr1(0)* GYZ1(0)
     $ - 24* HZr3(0,1,0)
     $ + 12* HZr2(1,0)
     $ - 24* HZr2(1,0)* GYZ1(2)
     $ + 24* HZr2(1,0)* GYZ1(0)
 
     $ + 24* HZr3(1,1,0)
     $ + 24* GYZ3(2,1,0)
     $ + 18* GYZ1(0)
     $ + 24* GYZ3(0,1,0)
     $ - 60* GYZ2(1,0)
     $ )
     $ + (z)/(y+z)*(
     $  6* HZr1(0)
     $ - 6* HZr1(0)* GYZ1(0)
     $ - 6* HZr2(1,0)
     $ - 6* GYZ1(0)
     $ + 6* GYZ2(1,0)
     $ - pi**2)
     $ + (z**2)/((1-y)**3)*(
     $  (-2*pi**2)/(3.0d0)* HZr1(0)
 
     $ - (2*pi**2)/(3.0d0)* HZr1(1)
     $ + (2*pi**2)/(3.0d0)* GYZ1(2)
     $ + 12*z3
     $ - 4* HZr1(0)* GYZ2(0,2)
     $ + 4* HZr3(0,0,1)
     $ - 4* HZr2(0,1)* GYZ1(2)
 
     $ - 4* HZr3(0,1,0)
     $ - 4* HZr1(1)* GYZ2(2,3)
     $ - 4* HZr1(1)* GYZ2(0,3)
     $ + 4* HZr2(1,0)* GYZ1(2)
     $ + 4* HZr3(1,0,1)
 
     $ - 4* HZr3(1,1,0)
     $ + 4* GYZ3(2,3,2)
     $ + 4* GYZ3(0,3,2)
     $ )
     $ + (z**2)/((1-y)**2)*(
     $ - 4* HZr1(0)* GYZ1(2)
 
     $ - 4* HZr1(1)* GYZ1(3)
     $ + 4* GYZ2(3,2)
     $ )
     $ + (z**2)/(1-y)*(
     $ - 2* HZr1(0)* GYZ1(2)
     $ - 2* HZr1(1)* GYZ1(3)
 
     $ + 2* GYZ2(3,2)
     $ )
     $ + (z**2)/((y+z)**4)*(
     $  6*pi**2* HZr1(1)
     $ - 6*pi**2* GYZ1(2)
     $ - 36* HZr1(0)* GYZ2(2,0)
     $ - 36* HZr3(0,1,0)
 
     $ - 36* HZr2(1,0)* GYZ1(2)
     $ + 36* HZr2(1,0)* GYZ1(0)
     $ + 36* HZr3(1,1,0)
     $ + 36* GYZ3(2,1,0)
     $ + 36* GYZ3(0,1,0)
     $ )
 
     $ + (z**2)/((y+z)**3)*(
     $ - 6*pi**2
     $ - 4*pi**2* HZr1(1)
     $ + 4*pi**2* GYZ1(2)
     $ + 24* HZr1(0)* GYZ2(2,0)
     $ - 36* HZr1(0)* GYZ1(0)
 
     $ + 24* HZr3(0,1,0)
     $ - 36* HZr2(1,0)
     $ + 24* HZr2(1,0)* GYZ1(2)
     $ - 24* HZr2(1,0)* GYZ1(0)
     $ - 24* HZr3(1,1,0)
 
     $ - 24* GYZ3(2,1,0)
     $ - 24* GYZ3(0,1,0)
     $ + 36* GYZ2(1,0)
     $ )
     $ + (z**2)/((y+z)**2)*(
     $  pi**2
     $ + 6* HZr1(0)* GYZ1(0)
     $ + 6* HZr2(1,0)
 
     $ - 6* GYZ2(1,0)
     $ )
     $ + (1.0d0)/(1-y-z)*(
     $  (-pi**2)/(3.0d0)
     $ - 2* HZr1(0)* GYZ1(0)
     $ - 2* HZr2(1,0)
     $ + 2* GYZ2(1,0)
     $ )
     $ + (1.0d0)/(2*(1-y))*(
     $ - 2*pi**2
 
     $ - (8*pi**2)/(3.0d0)* HZr1(0)
     $ - (4*pi**2)/(3.0d0)* HZr1(1)
     $ + (4*pi**2)/(3.0d0)* GYZ1(2)
     $ + (2*pi**2)/(3.0d0)* GYZ1(0)
     $ + (8*pi**2)/(3.0d0)* GYZ1(1)
     $ + 64*z3
     $ - 8* HZr1(0)* GYZ2(2,0)
 
     $ - 16* HZr1(0)* GYZ2(0,2)
     $ + 6* HZr1(0)* GYZ1(0)
     $ + 4* HZr1(0)* GYZ2(0,0)
     $ + 8* HZr1(0)* GYZ2(1,0)
     $ + 20* HZr3(0,0,1)
 
     $ + 12* HZr2(0,1)
     $ - 20* HZr2(0,1)* GYZ1(2)
     $ - 4* HZr2(0,1)* GYZ1(0)
     $ - 8* HZr3(0,1,0)
     $ - 20* HZr1(1)* GYZ2(2,3)
 
     $ + 12* HZr1(1)* GYZ1(3)
     $ - 20* HZr1(1)* GYZ2(0,3)
     $ - 2* HZr1(1)* GYZ1(0)
     $ + 4* HZr1(1)* GYZ2(0,0)
     $ - 10* HZr2(1,0)
 
     $ + 8* HZr2(1,0)* GYZ1(2)
     $ + 8* HZr2(1,0)* GYZ1(0)
     $ + 20* HZr3(1,0,1)
     $ - 8* HZr3(1,1,0)
     $ + 20* GYZ3(2,3,2)
 
     $ + 2* GYZ2(2,0)
     $ - 4* GYZ3(2,0,0)
     $ + 12* GYZ3(2,1,0)
     $ - 12* GYZ2(3,2)
     $ + 2* GYZ2(0,2)
 
     $ - 4* GYZ3(0,2,0)
     $ + 20* GYZ3(0,3,2)
     $ - 5* GYZ1(0)
     $ - 4* GYZ3(0,0,2)
     $ - 4* GYZ2(0,0)
     $ - 16* GYZ3(1,1,0)
     $ )
 
      fin20 = fin20
     $ + (1.0d0)/(2*(y+z)**2)*(
     $   (22*pi**2)/(3.0d0)* HZr1(1)
     $ - (22*pi**2)/(3.0d0)* GYZ1(2)
     $ - 2* HZr1(0)* GYZ2(2,2)
     $ + 3* HZr1(0)* GYZ1(2)
 
     $ - 44* HZr1(0)* GYZ2(2,0)
     $ - 2* HZr1(0)* GYZ2(3,2)
     $ + 2* HZr1(0)* GYZ2(0,2)
     $ + 3* HZr2(0,1)
 
     $ - 2* HZr2(0,1)* GYZ1(2)
     $ - 2* HZr2(0,1)* GYZ1(3)
     $ - 44* HZr3(0,1,0)
     $ + 2* HZr3(0,1,1)
     $ + 27* HZr1(1)
 
     $ - 4* HZr1(1)* GYZ2(2,3)
     $ + 6* HZr1(1)* GYZ1(2)
     $ + 2* HZr1(1)* GYZ2(2,0)
     $ - 4* HZr1(1)* GYZ2(3,2)
 
     $ - 4* HZr1(1)* GYZ2(3,3)
     $ + 6* HZr1(1)* GYZ1(3)
     $ + 2* HZr1(1)* GYZ2(3,0)
     $ + 2* HZr1(1)* GYZ2(0,2)
 
     $ + 2* HZr1(1)* GYZ2(0,3)
     $ - 3* HZr1(1)* GYZ1(0)
     $ - 3* HZr2(1,0)
     $ - 44* HZr2(1,0)* GYZ1(2)
     $ + 2* HZr2(1,0)* GYZ1(3)
 
     $ + 44* HZr2(1,0)* GYZ1(0)
     $ + 2* HZr3(1,0,1)
     $ - 6* HZr2(1,1)
     $ + 4* HZr2(1,1)* GYZ1(3)
     $ - 2* HZr2(1,1)* GYZ1(0)
 
     $ + 44* HZr3(1,1,0)
     $ - 6* GYZ2(2,2)
     $ - 2* GYZ3(2,2,0)
     $ + 4* GYZ3(2,3,2)
     $ - 27* GYZ1(2)
 
     $ - 2* GYZ3(2,0,2)
     $ + 3* GYZ2(2,0)
     $ + 46* GYZ3(2,1,0)
     $ + 4* GYZ3(3,2,2)
     $ - 6* GYZ2(3,2)
 
     $ - 2* GYZ3(3,2,0)
     $ + 4* GYZ3(3,3,2)
     $ - 2* GYZ3(3,0,2)
     $ - 2* GYZ3(0,2,2)
 
     $ + 3* GYZ2(0,2)
     $ - 2* GYZ3(0,3,2)
     $ + 46* GYZ3(0,1,0)
     $ )
     $ + (1.0d0)/(2*(y+z))*(
     $ (-22*pi**2)/(3.0d0)
     $ - 2*pi**2* HZr1(1)
 
     $ + 2*pi**2* GYZ1(2)
     $ - 17
     $ - 5* HZr1(0)
     $ - 4* HZr1(0)* GYZ1(2)
     $ + 12* HZr1(0)* GYZ2(2,0)
     $ - 44* HZr1(0)* GYZ1(0)
 
     $ - 4* HZr2(0,1)
     $ + 12* HZr3(0,1,0)
     $ + 4* HZr1(1)
     $ - 4* HZr1(1)* GYZ1(2)
     $ - 8* HZr1(1)* GYZ1(3)
     $ + 4* HZr1(1)* GYZ1(0)
 
     $ - 42* HZr2(1,0)
     $ + 12* HZr2(1,0)* GYZ1(2)
     $ - 12* HZr2(1,0)* GYZ1(0)
     $ + 4* HZr2(1,1)
     $ - 12* HZr3(1,1,0)
 
     $ + 4* GYZ2(2,2)
     $ - 4* GYZ1(2)
     $ - 4* GYZ2(2,0)
     $ - 12* GYZ3(2,1,0)
     $ + 8* GYZ2(3,2)
 
     $ - 4* GYZ2(0,2)
     $ - 5* GYZ1(0)
     $ - 12* GYZ3(0,1,0)
     $ + 46* GYZ2(1,0)
     $ )
     $ + (T*pi**2)/(24.0d0)*(
     $ + 29
     $ + 6* HZr1(1)
 
     $ + 8* HZr1(1)* GYZ1(2)
     $ - 16* GYZ2(2,2)
     $ + 12* GYZ1(2)
     $ + 8* GYZ2(2,1)
     $ + 8* GYZ2(0,2)
 
     $ - 8* GYZ2(0,1)
     $ - 18* GYZ1(1)
     $ + 8* GYZ2(1,1)
     $ )
     $ + (T)/(8.0d0)*(
     $ - (22*pi**4)/(45.0d0)
     $ + (255.0d0)/(4.0d0)
     $ - 60*z3
     $ - 16*z3* HZr1(1)
     $ - 16*z3* GYZ1(2)
 
     $ + 32*z3* GYZ1(1)
     $ + 18* HZr1(0)* GYZ2(2,2)
     $ + 15* HZr1(0)* GYZ1(2)
     $ - 8* HZr1(0)* GYZ3(2,0,2)
 
     $ - 16* HZr1(0)* GYZ3(3,2,2)
     $ + 42* HZr1(0)* GYZ2(3,2)
     $ - 16* HZr1(0)* GYZ3(3,3,2)
 
     $ + 16* HZr1(0)* GYZ3(3,0,2)
     $ - 42* HZr1(0)* GYZ2(0,2)
     $ + 16* HZr1(0)* GYZ3(0,3,2)
 
     $ - 16* HZr1(0)* GYZ3(0,0,2)
     $ + 16* HZr2(0,0)* GYZ2(2,2)
     $ - 16* HZr3(0,0,1)* GYZ1(2)
 
     $ + 8* HZr3(0,0,1)* GYZ1(1)
     $ + 16* HZr4(0,0,1,1)
     $ + 15* HZr2(0,1)
     $ + 32* HZr2(0,1)* GYZ2(2,2)
 
     $ + 6* HZr2(0,1)* GYZ1(2)
     $ + 8* HZr2(0,1)* GYZ2(2,0)
     $ - 16* HZr2(0,1)* GYZ2(3,2)
 
     $ - 16* HZr2(0,1)* GYZ2(3,3)
     $ + 42* HZr2(0,1)* GYZ1(3)
     $ + 8* HZr2(0,1)* GYZ2(0,2)
     $ - 8* HZr2(0,1)* GYZ2(1,2)
 
     $ - 8* HZr2(0,1)* GYZ2(1,0)
     $ - 8* HZr3(0,1,0)* GYZ1(2)
     $ + 16* HZr4(0,1,0,1)
     $ - 18* HZr3(0,1,1)
 
     $ + 16* HZr3(0,1,1)* GYZ1(3)
     $ - 16* HZr3(0,1,1)* GYZ1(0)
     $ + 16* HZr4(0,1,1,0)
     $ - 52* HZr1(1)
 
     $ + 32* HZr1(1)* GYZ3(2,2,3)
     $ + 24* HZr1(1)* GYZ2(2,3)
     $ - 18* HZr1(1)* GYZ1(2)
 
     $ - 18* HZr1(1)* GYZ2(2,0)
     $ - 16* HZr1(1)* GYZ3(2,0,0)
     $ - 8* HZr1(1)* GYZ3(2,1,0)
 
     $ - 32* HZr1(1)* GYZ3(3,2,3)
     $ + 36* HZr1(1)* GYZ2(3,2)
     $ + 16* HZr1(1)* GYZ3(3,2,0)
 
     $ - 32* HZr1(1)* GYZ3(3,3,2)
     $ - 32* HZr1(1)* GYZ3(3,3,3)
     $ + 84* HZr1(1)* GYZ2(3,3)
 
     $ + 16* HZr1(1)* GYZ3(3,3,0)
     $ + 30* HZr1(1)* GYZ1(3)
     $ + 16* HZr1(1)* GYZ3(3,0,2)
 
     $ + 16* HZr1(1)* GYZ3(3,0,3)
     $ - 42* HZr1(1)* GYZ2(3,0)
     $ + 8* HZr1(1)* GYZ3(0,2,3)
 
     $ - 18* HZr1(1)* GYZ2(0,2)
     $ - 16* HZr1(1)* GYZ3(0,2,0)
     $ + 16* HZr1(1)* GYZ3(0,3,2)
 
     $ + 16* HZr1(1)* GYZ3(0,3,3)
     $ - 42* HZr1(1)* GYZ2(0,3)
     $ - 15* HZr1(1)* GYZ1(0)
     $ - 16* HZr1(1)* GYZ3(0,0,2)
 
     $ - 16* HZr1(1)* GYZ3(0,0,3)
     $ + 8* HZr1(1)* GYZ3(0,1,0)
     $ - 8* HZr1(1)* GYZ3(1,2,3)
     $ - 8* HZr1(1)* GYZ3(1,0,3)
 
     $ + 16* HZr1(1)* GYZ3(1,0,0)
     $ + 13* HZr2(1,0)
     $ + 24* HZr2(1,0)* GYZ1(2)
     $ + 16* HZr2(1,0)* GYZ2(3,3)
 
     $ - 42* HZr2(1,0)* GYZ1(3)
     $ + 16* HZr2(1,0)* GYZ2(0,2)
     $ - 16* HZr2(1,0)* GYZ2(0,3)
     $ + 8* HZr4(1,0,0,1)
 
     $ - 6* HZr3(1,0,1)
     $ - 24* HZr3(1,0,1)* GYZ1(2)
     $ + 16* HZr3(1,0,1)* GYZ1(3)
     $ - 8* HZr3(1,0,1)* GYZ1(0)
 
     $ + 8* HZr3(1,0,1)* GYZ1(1)
     $ + 18* HZr2(1,1)
     $ + 32* HZr2(1,1)* GYZ2(3,3)
     $ - 36* HZr2(1,1)* GYZ1(3)
 
     $ - 16* HZr2(1,1)* GYZ2(3,0)
     $ - 16* HZr2(1,1)* GYZ2(0,3)
     $ + 18* HZr2(1,1)* GYZ1(0)
     $ + 16* HZr2(1,1)* GYZ2(0,0)
 
     $ + 12* HZr3(1,1,0)
     $ - 8* HZr3(1,1,0)* GYZ1(2)
     $ + 16* HZr4(1,1,0,1)
     $ + 16* HZr4(1,1,1,0)
 
     $ - 32* GYZ4(2,2,3,2)
     $ + 18* GYZ2(2,2)
     $ + 18* GYZ3(2,2,0)
 
     $ + 16* GYZ4(2,2,0,0)
     $ - 24* GYZ3(2,3,2)
     $ + 52* GYZ1(2)
     $ + 18* GYZ3(2,0,2)
 
     $ + 16* GYZ4(2,0,2,0)
     $ + 15* GYZ2(2,0)
     $ + 16* GYZ4(2,0,0,2)
     $ - 8* GYZ4(2,0,1,0)
 
     $ + 8* GYZ4(2,1,2,0)
     $ + 8* GYZ4(2,1,0,2)
     $ - 42* GYZ3(2,1,0)
     $ - 16* GYZ4(2,1,0,0)
 
     $ - 16* GYZ4(2,1,1,0)
     $ - 36* GYZ3(3,2,2)
     $ - 16* GYZ4(3,2,2,0)
 
     $ + 32* GYZ4(3,2,3,2)
     $ - 30* GYZ2(3,2)
     $ - 16* GYZ4(3,2,0,2)
     $ + 42* GYZ3(3,2,0)
 
     $ + 16* GYZ4(3,2,1,0)
     $ + 32* GYZ4(3,3,2,2)
     $ - 84* GYZ3(3,3,2)
 
     $ - 16* GYZ4(3,3,2,0)
     $ + 32* GYZ4(3,3,3,2)
     $ - 16* GYZ4(3,3,0,2)
 
     $ - 16* GYZ4(3,0,2,2)
     $ + 42* GYZ3(3,0,2)
     $ - 16* GYZ4(3,0,3,2)
     $ + 16* GYZ4(3,0,1,0)
 
     $ + 18* GYZ3(0,2,2)
     $ + 16* GYZ4(0,2,2,0)
     $ - 8* GYZ4(0,2,3,2)
     $ + 15* GYZ2(0,2)
 
     $ + 16* GYZ4(0,2,0,2)
     $ - 8* GYZ4(0,2,1,0)
     $ - 16* GYZ4(0,3,2,2)
     $ + 42* GYZ3(0,3,2)
 
     $ - 16* GYZ4(0,3,3,2)
     $ + 16* GYZ4(0,0,2,2)
     $ + 16* GYZ4(0,0,3,2)
     $ - 8* GYZ4(0,1,2,0)
 
     $ - 8* GYZ4(0,1,0,2)
     $ - 42* GYZ3(0,1,0)
     $ + 16* GYZ4(0,1,1,0)
     $ + 8* GYZ4(1,2,3,2)
 
     $ - 16* GYZ4(1,2,0,0)
     $ + 8* GYZ4(1,2,1,0)
     $ - 16* GYZ4(1,0,2,0)
     $ + 8* GYZ4(1,0,3,2)
     $ - 28* GYZ2(1,0)
 
     $ - 16* GYZ4(1,0,0,2)
     $ + 8* GYZ4(1,0,1,0)
     $ + 36* GYZ3(1,1,0)
     $ + 16* GYZ4(1,1,0,0)
     $ - 16* GYZ4(1,1,1,0)
     $ )
 
     $ + (pi**2)/(6.0d0)*(
     $   9
     $ - 7* HZr1(0)
     $ + 6* HZr1(0)* GYZ1(2)
     $ - 2* HZr1(0)* GYZ1(1)
     $ + 4* HZr2(0,1)
     $ - 24* HZr1(1)
 
     $ + 12* HZr1(1)* GYZ1(2)
     $ + 8* HZr1(1)* GYZ1(3)
     $ - 4* HZr1(1)* GYZ1(0)
     $ - 6* HZr1(1)* GYZ1(1)
     $ - 2* HZr2(1,0)
     $ - 6* HZr2(1,1)
 
     $ - 20* GYZ2(2,2)
     $ + 34* GYZ1(2)
     $ + 6* GYZ2(2,0)
     $ + 8* GYZ2(2,1)
     $ - 8* GYZ2(3,2)
 
     $ + 8* GYZ2(0,2)
     $ - 7* GYZ1(0)
     $ - 4* GYZ2(0,1)
     $ + 6* GYZ2(1,2)
     $ - 10* GYZ1(1)
     $ - 4* GYZ2(1,0)
     $ )
 
     $ + (1.0d0)/(4.0d0)*( 
     $   40*z3
     $ + 40*z3* HZr1(1)
     $ - 80*z3* GYZ1(2)
     $ + 40*z3* GYZ1(1)
     $ - 2
     $ - 29* HZr1(0)
     $ + 20* HZr1(0)* GYZ2(2,2)
 
     $ - 32* HZr1(0)* GYZ3(2,2,0)
     $ - 4* HZr1(0)* GYZ1(2)
     $ + 32* HZr1(0)* GYZ2(2,0)
 
     $ + 16* HZr1(0)* GYZ3(2,0,0)
     $ + 16* HZr1(0)* GYZ3(2,1,0)
     $ - 24* HZr1(0)* GYZ3(3,2,2)
 
     $ + 52* HZr1(0)* GYZ2(3,2)
     $ - 32* HZr1(0)* GYZ3(3,2,0)
     $ - 24* HZr1(0)* GYZ3(3,3,2)
 
     $ + 24* HZr1(0)* GYZ3(3,0,2)
     $ + 8* HZr1(0)* GYZ3(0,2,2)
     $ - 44* HZr1(0)* GYZ2(0,2)
 
     $ + 16* HZr1(0)* GYZ3(0,2,0)
     $ + 24* HZr1(0)* GYZ3(0,3,2)
     $ + 40* HZr1(0)* GYZ1(0)
     $ - 24* HZr1(0)* GYZ3(0,0,2)
 
     $ - 20* HZr1(0)* GYZ2(0,0)
     $ + 16* HZr1(0)* GYZ3(1,2,0)
     $ - 8* HZr1(0)* GYZ3(1,0,2)
     $ + 4* HZr1(0)* GYZ2(1,0)
 
     $ - 16* HZr1(0)* GYZ3(1,0,0)
     $ - 16* HZr1(0)* GYZ3(1,1,0)
     $ + 20* HZr2(0,0)
     $ + 16* HZr2(0,0)* GYZ2(2,2)
 
     $ - 4* HZr2(0,0)* GYZ1(2)
     $ + 16* HZr2(0,0)* GYZ2(2,0)
     $ - 20* HZr2(0,0)* GYZ1(0)
     $ + 36* HZr3(0,0,1)
 
     $ - 16* HZr3(0,0,1)* GYZ1(2)
     $ + 16* HZr3(0,0,1)* GYZ1(1)
     $ - 16* HZr4(0,0,1,0)
     $ + 16* HZr4(0,0,1,1)
 
     $ + 48* HZr2(0,1)* GYZ2(2,2)
     $ - 52* HZr2(0,1)* GYZ1(2)
     $ - 24* HZr2(0,1)* GYZ2(3,2)
 
     $ - 24* HZr2(0,1)* GYZ2(3,3)
     $ + 84* HZr2(0,1)* GYZ1(3)
     $ + 8* HZr2(0,1)* GYZ2(0,2)
     $ - 8* HZr2(0,1)* GYZ1(0)
 
     $ - 16* HZr2(0,1)* GYZ2(1,2)
     $ - 8* HZr2(0,1)* GYZ2(1,0)
     $ - 28* HZr3(0,1,0)
     $ + 32* HZr3(0,1,0)* GYZ1(2)
 
     $ - 32* HZr3(0,1,0)* GYZ1(3)
     $ + 16* HZr3(0,1,0)* GYZ1(0)
     $ - 24* HZr3(0,1,0)* GYZ1(1)
     $ + 16* HZr4(0,1,0,1)
 
     $ - 20* HZr3(0,1,1)
     $ + 24* HZr3(0,1,1)* GYZ1(3)
     $ - 16* HZr3(0,1,1)* GYZ1(0)
     $ + 32* HZr4(0,1,1,0)
     $ - 58* HZr1(1)
 
     $ + 48* HZr1(1)* GYZ3(2,2,3)
     $ - 32* HZr1(1)* GYZ2(2,3)
     $ - 4* HZr1(1)* GYZ1(2)
 
     $ - 20* HZr1(1)* GYZ2(2,0)
     $ - 16* HZr1(1)* GYZ3(2,0,0)
     $ - 8* HZr1(1)* GYZ3(2,1,0)
 
     $ - 48* HZr1(1)* GYZ3(3,2,3)
     $ + 40* HZr1(1)* GYZ2(3,2)
     $ + 24* HZr1(1)* GYZ3(3,2,0)
 
     $ - 48* HZr1(1)* GYZ3(3,3,2)
     $ - 48* HZr1(1)* GYZ3(3,3,3)
     $ + 136* HZr1(1)* GYZ2(3,3)
 
     $ + 24* HZr1(1)* GYZ3(3,3,0)
     $ - 4* HZr1(1)* GYZ1(3)
     $ + 24* HZr1(1)* GYZ3(3,0,2)
     $ + 24* HZr1(1)* GYZ3(3,0,3)
 
     $ - 52* HZr1(1)* GYZ2(3,0)
     $ + 16* HZr1(1)* GYZ3(0,2,3)
     $ - 20* HZr1(1)* GYZ2(0,2)
 
     $ - 16* HZr1(1)* GYZ3(0,2,0)
     $ + 24* HZr1(1)* GYZ3(0,3,2)
     $ + 24* HZr1(1)* GYZ3(0,3,3)
 
     $ - 52* HZr1(1)* GYZ2(0,3)
     $ + 4* HZr1(1)* GYZ1(0)
     $ - 16* HZr1(1)* GYZ3(0,0,2)
     $ - 24* HZr1(1)* GYZ3(0,0,3)
 
     $ + 4* HZr1(1)* GYZ2(0,0)
     $ + 8* HZr1(1)* GYZ3(0,1,0)
     $ - 16* HZr1(1)* GYZ3(1,2,3)
     $ - 16* HZr1(1)* GYZ3(1,0,3)
 
     $ + 24* HZr1(1)* GYZ2(1,0)
     $ + 16* HZr1(1)* GYZ3(1,0,0)
     $ + 60* HZr2(1,0)
     $ - 40* HZr2(1,0)* GYZ2(2,2)
 
     $ + 80* HZr2(1,0)* GYZ1(2)
     $ + 16* HZr2(1,0)* GYZ2(2,0)
     $ - 32* HZr2(1,0)* GYZ2(3,2)
 
     $ + 24* HZr2(1,0)* GYZ2(3,3)
     $ - 52* HZr2(1,0)* GYZ1(3)
     $ + 32* HZr2(1,0)* GYZ2(3,0)
     $ + 32* HZr2(1,0)* GYZ2(0,2)
 
     $ - 24* HZr2(1,0)* GYZ2(0,3)
     $ - 36* HZr2(1,0)* GYZ1(0)
     $ - 16* HZr2(1,0)* GYZ2(0,0)
     $ + 24* HZr2(1,0)* GYZ2(1,2)
 
     $ - 16* HZr2(1,0)* GYZ2(1,0)
     $ - 20* HZr3(1,0,0)
     $ + 16* HZr3(1,0,0)* GYZ1(2)
     $ + 16* HZr4(1,0,0,1)
     $ + 28* HZr3(1,0,1)
 
     $ - 40* HZr3(1,0,1)* GYZ1(2)
     $ + 24* HZr3(1,0,1)* GYZ1(3)
     $ - 8* HZr3(1,0,1)* GYZ1(0)
     $ + 16* HZr3(1,0,1)* GYZ1(1)
 
     $ - 24* HZr4(1,0,1,0)
     $ + 4* HZr2(1,1)
     $ + 48* HZr2(1,1)* GYZ2(3,3)
     $ - 40* HZr2(1,1)* GYZ1(3)
 
     $ - 24* HZr2(1,1)* GYZ2(3,0)
     $ - 24* HZr2(1,1)* GYZ2(0,3)
     $ + 20* HZr2(1,1)* GYZ1(0)
     $ + 16* HZr2(1,1)* GYZ2(0,0)
 
     $ - 64* HZr3(1,1,0)
     $ + 16* HZr3(1,1,0)* GYZ1(2)
     $ + 32* HZr3(1,1,0)* GYZ1(3)
     $ - 16* HZr3(1,1,0)* GYZ1(0)
 
     $ - 24* HZr3(1,1,0)* GYZ1(1)
     $ + 24* HZr4(1,1,0,1)
     $ - 8* HZr4(1,1,1,0)
     $ - 48* GYZ4(2,2,3,2)
 
     $ + 4* GYZ2(2,2)
     $ + 20* GYZ3(2,2,0)
     $ + 16* GYZ4(2,2,0,0)
     $ + 40* GYZ4(2,2,1,0)
 
     $ + 32* GYZ3(2,3,2)
     $ + 58* GYZ1(2)
     $ + 20* GYZ3(2,0,2)
     $ + 16* GYZ4(2,0,2,0)
 
     $ - 4* GYZ2(2,0)
     $ + 16* GYZ4(2,0,0,2)
     $ - 4* GYZ3(2,0,0)
     $ - 48* GYZ4(2,0,1,0)
 
     $ + 8* GYZ4(2,1,2,0)
     $ + 8* GYZ4(2,1,0,2)
     $ - 100* GYZ3(2,1,0)
     $ - 32* GYZ4(2,1,0,0)
 
     $ - 32* GYZ4(2,1,1,0)
     $ - 40* GYZ3(3,2,2)
     $ - 24* GYZ4(3,2,2,0)
 
     $ + 48* GYZ4(3,2,3,2)
     $ + 4* GYZ2(3,2)
     $ - 24* GYZ4(3,2,0,2)
 
     $ + 52* GYZ3(3,2,0)
     $ + 56* GYZ4(3,2,1,0)
     $ + 48* GYZ4(3,3,2,2)
 
     $ - 136* GYZ3(3,3,2)
     $ - 24* GYZ4(3,3,2,0)
     $ + 48* GYZ4(3,3,3,2)
 
     $ - 24* GYZ4(3,3,0,2)
     $ - 24* GYZ4(3,0,2,2)
     $ + 52* GYZ3(3,0,2)
     $ - 24* GYZ4(3,0,3,2)
 
     $ + 56* GYZ4(3,0,1,0)
     $ + 20* GYZ3(0,2,2)
     $ + 16* GYZ4(0,2,2,0)
     $ - 16* GYZ4(0,2,3,2)
 
     $ - 4* GYZ2(0,2)
     $ + 16* GYZ4(0,2,0,2)
     $ - 4* GYZ3(0,2,0)
     $ - 24* GYZ4(0,2,1,0)
 
     $ - 24* GYZ4(0,3,2,2)
     $ + 52* GYZ3(0,3,2)
     $ - 24* GYZ4(0,3,3,2)
     $ - 29* GYZ1(0)
 
     $ + 16* GYZ4(0,0,2,2)
     $ - 4* GYZ3(0,0,2)
     $ + 24* GYZ4(0,0,3,2)
     $ + 20* GYZ2(0,0)
     $ - 16* GYZ4(0,0,1,0)
 
     $ - 8* GYZ4(0,1,2,0)
     $ - 8* GYZ4(0,1,0,2)
     $ - 20* GYZ3(0,1,0)
     $ + 16* GYZ4(0,1,1,0)
 
     $ + 16* GYZ4(1,2,3,2)
     $ - 24* GYZ3(1,2,0)
     $ - 16* GYZ4(1,2,0,0)
     $ - 8* GYZ4(1,2,1,0)
 
     $ - 24* GYZ3(1,0,2)
     $ - 16* GYZ4(1,0,2,0)
     $ + 16* GYZ4(1,0,3,2)
     $ - 56* GYZ2(1,0)
     $ - 16* GYZ4(1,0,0,2)
 
     $ + 24* GYZ3(1,0,0)
     $ + 24* GYZ4(1,0,1,0)
     $ + 40* GYZ3(1,1,0)
     $ + 32* GYZ4(1,1,0,0)
     $ ) 


      else if(i.eq.4) then
      fin20=
     $(z)/(3*y)*(
     $ - 2* HZr1(0)
     $ - 3* HZr2(1,0)
     $ - 3* GYZ2(1,0)
     $)
     $ + (1.0d0)/(y*(y+z))*(
     $ - 2* HZr2(1,0)
     $ - 2* GYZ2(1,0)
     $)
     $ + (1.0d0)/(18*y)*(
     $ - 74

     $ + 15* HZr1(0)
     $ + 36* HZr2(1,0)
     $ + 15* GYZ1(0)
     $ + 36* GYZ2(1,0)
     $)
     $ + (z)/(18*(1-y)**2)*(
     $   2*pi**2
     $ - 3* HZr1(0)* GYZ1(0)
     $ + 50* GYZ1(0)

     $ - 18* GYZ2(0,0)
     $ - 12* GYZ2(1,0)
     $)
     $ + (z)/(18*(1-y))*(
     $   6*pi**2
     $ + 38
     $ - 3* HZr1(0)
     $ - 9* HZr1(0)* GYZ1(0)
     $ + 87* GYZ1(0)

     $ - 54* GYZ2(0,0)
     $ - 36* GYZ2(1,0)
     $)
     $ + (z)/((y+z)**3)*(
     $ - 2*pi**2* HZr1(1)
     $ + 2*pi**2* GYZ1(2)
     $ + 12* HZr1(0)* GYZ2(2,0)

     $ + 12* HZr3(0,1,0)
     $ + 12* HZr2(1,0)
     $ + 12* HZr2(1,0)* GYZ1(2)
     $ - 12* HZr2(1,0)* GYZ1(0)
     $ - 12* HZr3(1,1,0)

     $ - 12* GYZ3(2,1,0)
     $ - 12* GYZ3(0,1,0)
     $ + 12* GYZ2(1,0)
     $)
     $ + (z)/((y+z)**2)*(
     $   2*pi**2 
     $ + (4*pi**2)/(3.0d0)* HZr1(1)
     $ - (4*pi**2)/(3.0d0)* GYZ1(2)

     $ - 6* HZr1(0)
     $ - 8* HZr1(0)* GYZ2(2,0)
     $ + 12* HZr1(0)* GYZ1(0)
     $ - 8* HZr3(0,1,0)
     $ + 4* HZr2(1,0)

     $ - 8* HZr2(1,0)* GYZ1(2)
     $ + 8* HZr2(1,0)* GYZ1(0)
     $ + 8* HZr3(1,1,0)
     $ + 8* GYZ3(2,1,0)
     $ + 6* GYZ1(0)

     $ + 8* GYZ3(0,1,0)
     $ - 20* GYZ2(1,0)
     $)
     $ + (z)/(y+z)*(
     $ - (pi**2)/(3.0d0)
     $ + 2* HZr1(0)
     $ - 2* HZr1(0)* GYZ1(0)
     $ - 2* HZr2(1,0)
     $ - 2* GYZ1(0)

     $ + 2* GYZ2(1,0)
     $)
     $ + (z**2)/((y+z)**4)*(
     $   2*pi**2* HZr1(1)
     $ - 2*pi**2* GYZ1(2)
     $ - 12* HZr1(0)* GYZ2(2,0)
     $ - 12* HZr3(0,1,0)

     $ - 12* HZr2(1,0)* GYZ1(2)
     $ + 12* HZr2(1,0)* GYZ1(0)
     $ + 12* HZr3(1,1,0)
     $ + 12* GYZ3(2,1,0)
     $ + 12* GYZ3(0,1,0)
     $)

     $ + (z**2)/((y+z)**3)*(
     $ - 2*pi**2 
     $ - (4*pi**2)/(3.0d0)* HZr1(1)
     $ + (4*pi**2)/(3.0d0)* GYZ1(2)
     $ + 8* HZr1(0)* GYZ2(2,0)
     $ - 12* HZr1(0)* GYZ1(0)

     $ + 8* HZr3(0,1,0)
     $ - 12* HZr2(1,0)
     $ + 8* HZr2(1,0)* GYZ1(2)
     $ - 8* HZr2(1,0)* GYZ1(0)
     $ - 8* HZr3(1,1,0)

     $ - 8* GYZ3(2,1,0)
     $ - 8* GYZ3(0,1,0)
     $ + 12* GYZ2(1,0)
     $)
     $ + (z**2)/(3*(y+z)**2)*(
     $   pi**2
     $ + 6* HZr1(0)* GYZ1(0)
     $ + 6* HZr2(1,0)

     $ - 6* GYZ2(1,0)
     $)
     $ + (1.0d0)/(9*(1-y))*(
     $ - 4*pi**2 
     $ + 6* HZr1(0)* GYZ1(0)
     $ - 70* GYZ1(0)
     $ + 36* GYZ2(0,0)
     $ + 24* GYZ2(1,0)
     $)

     $ + (1.0d0)/((y+z)**2)*(
     $   (2*pi**2)/(3.0d0)* HZr1(1)
     $ - (2*pi**2)/(3.0d0)* GYZ1(2)
     $ - 4* HZr1(0)* GYZ2(2,0)
     $ - 4* HZr3(0,1,0)
     $ - 4* HZr2(1,0)* GYZ1(2)

     $ + 4* HZr2(1,0)* GYZ1(0)
     $ + 4* HZr3(1,1,0)
     $ + 4* GYZ3(2,1,0)
     $ + 4* GYZ3(0,1,0)
     $)
     $ + (1.0d0)/(y+z)*(
     $ - (2*pi**2)/(3.0d0)
     $ - (pi**2)/(3.0d0)* HZr1(1)

     $ + (pi**2)/(3.0d0)* GYZ1(2)
     $ + 2
     $ - HZr1(0)
     $ + 2* HZr1(0)* GYZ2(2,0)
     $ - 4* HZr1(0)* GYZ1(0)
     $ + 2* HZr3(0,1,0)
     $ - 4* HZr2(1,0)

     $ + 2* HZr2(1,0)* GYZ1(2)
     $ - 2* HZr2(1,0)* GYZ1(0)
     $ - 2* HZr3(1,1,0)
     $ - 2* GYZ3(2,1,0)
     $ - GYZ1(0)
     $ - 2* GYZ3(0,1,0)

     $ + 4* GYZ2(1,0)
     $)
     $ + (T*pi**2)/(216.0d0)*(
     $ + 431
     $ - 12* HZr1(0)
     $ + 24* GYZ1(2)
     $ - 12* GYZ1(0)
     $ - 24* GYZ1(1)
     $)
     $ + (T)/(18.0d0)*(
     $   (4345.0d0)/(36.0d0)
     $ - 38*z3

     $ + 31* HZr1(0)
     $ + 12* HZr1(0)* GYZ2(2,0)
     $ + 10* HZr1(0)* GYZ1(0)
     $ - 18* HZr1(0)* GYZ2(0,0)
     $ + 3* HZr1(0)* GYZ2(1,0)

     $ - 41* HZr2(0,0)
     $ - 18* HZr2(0,0)* GYZ1(0)
     $ - 3* HZr3(0,1,0)
     $ + 29* HZr2(1,0)
     $ + 12* HZr2(1,0)* GYZ1(2)

     $ - 15* HZr2(1,0)* GYZ1(0)
     $ - 18* HZr3(1,0,0)
     $ - 12* GYZ3(2,1,0)
     $ + 31* GYZ1(0)
     $ - 41* GYZ2(0,0)
     $ + 3* GYZ3(0,1,0)

     $ - 29* GYZ2(1,0)
     $ + 18* GYZ3(1,0,0)
     $ + 12* GYZ3(1,1,0)
     $)
     $ + (pi**2)/(3.0d0)
     $ + 2* HZr1(0)* GYZ1(0)
     $ + 2* HZr2(1,0)
     $ - 2* GYZ2(1,0)


      else if(i.eq.5) then
      fin20=
     $(z)/(9*y)*(
     $ - 2*pi**2* GYZ1(2)
     $ + 12* HZr1(0)* GYZ2(2,2)
     $ - 47* HZr1(0)* GYZ1(2)
     $ + 3* HZr1(0)* GYZ2(2,0)

     $ - 15* HZr1(0)* GYZ2(3,2)
     $ + 3* HZr1(0)* GYZ2(0,2)
     $ + 18* HZr2(0,0)* GYZ1(2)
     $ + 9* HZr2(0,1)

     $ + 3* HZr2(0,1)* GYZ1(3)
     $ + 12* HZr1(1)* GYZ2(2,3)
     $ + 12* HZr1(1)* GYZ2(3,2)
     $ - 12* HZr1(1)* GYZ2(3,3)

     $ - 38* HZr1(1)* GYZ1(3)
     $ + 3* HZr1(1)* GYZ2(3,0)
     $ + 3* HZr1(1)* GYZ2(0,3)
     $ - 9* HZr1(1)* GYZ1(0)
     $ + 9* HZr2(1,0)

     $ - 12* HZr2(1,0)* GYZ1(2)
     $ + 15* HZr2(1,0)* GYZ1(3)
     $ - 12* HZr2(1,1)* GYZ1(3)
     $ - 12* GYZ3(2,3,2)

     $ + 9* GYZ2(2,0)
     $ - 12* GYZ3(3,2,2)
     $ + 38* GYZ2(3,2)
     $ - 3* GYZ3(3,2,0)

     $ + 12* GYZ3(3,3,2)
     $ - 3* GYZ3(3,0,2)
     $ + 9* GYZ2(0,2)
     $ - 3* GYZ3(0,3,2)
     $)

     $ + (1.0d0)/(y*(y+z))*(
     $ - 2* HZr1(0)* GYZ1(2)
     $ + 2* HZr2(0,1)
     $ - 2* HZr1(1)* GYZ1(0)
     $ + 2* HZr2(1,0)
     $ + 2* GYZ2(2,0)

     $ + 2* GYZ2(0,2)
     $)
     $ + (1.0d0)/(18*y)*(
     $   8*pi**2* GYZ1(2)
     $ - 38
     $ + 3* HZr1(0)
     $ - 48* HZr1(0)* GYZ2(2,2)

     $ + 188* HZr1(0)* GYZ1(2)
     $ - 12* HZr1(0)* GYZ2(2,0)
     $ + 60* HZr1(0)* GYZ2(3,2)
     $ - 12* HZr1(0)* GYZ2(0,2)

     $ - 72* HZr2(0,0)* GYZ1(2)
     $ - 36* HZr2(0,1)
     $ - 12* HZr2(0,1)* GYZ1(3)
     $ - 48* HZr1(1)* GYZ2(2,3)

     $ - 48* HZr1(1)* GYZ2(3,2)
     $ + 48* HZr1(1)* GYZ2(3,3)
     $ + 152* HZr1(1)* GYZ1(3)
     $ - 12* HZr1(1)* GYZ2(3,0)

     $ - 12* HZr1(1)* GYZ2(0,3)
     $ + 36* HZr1(1)* GYZ1(0)
     $ - 36* HZr2(1,0)
     $ + 48* HZr2(1,0)* GYZ1(2)

     $ - 60* HZr2(1,0)* GYZ1(3)
     $ + 48* HZr2(1,1)* GYZ1(3)
     $ + 48* GYZ3(2,3,2)
     $ - 36* GYZ2(2,0)

     $ + 48* GYZ3(3,2,2)
     $ - 152* GYZ2(3,2)
     $ + 12* GYZ3(3,2,0)
     $ - 48* GYZ3(3,3,2)

     $ + 12* GYZ3(3,0,2)
     $ - 36* GYZ2(0,2)
     $ + 12* GYZ3(0,3,2)
     $ + 15* GYZ1(0)
     $)
     $ + (z)/(18*(1-y)**2)*(
     $ - 2*pi**2

     $ + 3* HZr1(0)* GYZ1(0)
     $ - 50* GYZ1(0)
     $ + 18* GYZ2(0,0)
     $ + 12* GYZ2(1,0)
     $)
     $ + (z)/(18*(1-y))*(
     $ - 6*pi**2
     $ - 38
     $ + 3* HZr1(0)

     $ + 9* HZr1(0)* GYZ1(0)
     $ - 87* GYZ1(0)
     $ + 54* GYZ2(0,0)
     $ + 36* GYZ2(1,0)
     $)
     $ + (1.0d0)/(9*(1-y))*(
     $   2*pi**2
     $ - 3* HZr1(0)* GYZ1(0)

     $ + 38* GYZ1(0)
     $ - 18* GYZ2(0,0)
     $ - 12* GYZ2(1,0)
     $)
     $ + (1.0d0)/(9*(y+z)**2)*(
     $ - 9* HZr1(0)* GYZ1(2)
     $ - 3* HZr2(0,1)
     $ - 26* HZr1(1)

     $ + 12* HZr1(1)* GYZ1(2)
     $ - 12* HZr1(1)* GYZ1(3)
     $ + 9* HZr1(1)* GYZ1(0)
     $ + 9* HZr2(1,0)
     $ - 12* HZr2(1,1)

     $ - 12* GYZ2(2,2)
     $ + 26* GYZ1(2)
     $ - 9* GYZ2(2,0)
     $ + 12* GYZ2(3,2)
     $ - 9* GYZ2(0,2)
     $)

     $ + (1.0d0)/(9*(y+z))*(
     $ + 38
     $ - 9* HZr1(0)
     $ - 12* HZr1(1)
     $ + 12* GYZ1(2)
     $ - 9* GYZ1(0)
     $)
     $ + (T*pi**2)/(36.0d0)*(
     $ - 7
     $ +  HZr1(1)

     $ - 5* GYZ1(2)
     $ + 4* GYZ1(1)
     $)
     $ + (T)/(108.0d0)*(
     $ - (4085.0d0)/(6.0d0)
     $ + 6*z3
     $ + 72* HZr1(0)
     $ + 72* HZr1(0)* GYZ2(2,2)

     $ - 147* HZr1(0)* GYZ1(2)
     $ + 36* HZr1(0)* GYZ2(2,0)
     $ - 108* HZr1(0)* GYZ2(3,2)
     $ + 36* HZr1(0)* GYZ2(0,2)

     $ - 18* HZr1(0)* GYZ2(1,0)
     $ + 108* HZr2(0,0)* GYZ1(2)
     $ - 36* HZr3(0,0,1)
     $ - 201* HZr2(0,1)

     $ + 72* HZr2(0,1)* GYZ1(2)
     $ - 36* HZr2(0,1)* GYZ1(3)
     $ + 72* HZr2(0,1)* GYZ1(0)
     $ + 18* HZr3(0,1,0)
     $ - 72* HZr3(0,1,1)

     $ + 68* HZr1(1)
     $ + 144* HZr1(1)* GYZ2(2,3)
     $ - 108* HZr1(1)* GYZ1(2)
     $ - 72* HZr1(1)* GYZ2(2,0)

     $ + 144* HZr1(1)* GYZ2(3,2)
     $ - 144* HZr1(1)* GYZ2(3,3)
     $ - 348* HZr1(1)* GYZ1(3)
     $ + 108* HZr1(1)* GYZ2(3,0)

     $ - 72* HZr1(1)* GYZ2(0,2)
     $ + 108* HZr1(1)* GYZ2(0,3)
     $ + 147* HZr1(1)* GYZ1(0)
     $ - 108* HZr1(1)* GYZ2(0,0)

     $ - 81* HZr2(1,0)
     $ - 72* HZr2(1,0)* GYZ1(2)
     $ + 108* HZr2(1,0)* GYZ1(3)
     $ - 18* HZr2(1,0)* GYZ1(0)
     $ - 72* HZr3(1,0,1)

     $ + 108* HZr2(1,1)
     $ - 144* HZr2(1,1)* GYZ1(3)
     $ + 72* HZr2(1,1)* GYZ1(0)
     $ + 108* GYZ2(2,2)

     $ + 72* GYZ3(2,2,0)
     $ - 144* GYZ3(2,3,2)
     $ - 68* GYZ1(2)
     $ + 72* GYZ3(2,0,2)

     $ - 147* GYZ2(2,0)
     $ + 108* GYZ3(2,0,0)
     $ - 144* GYZ3(3,2,2)
     $ + 348* GYZ2(3,2)

     $ - 108* GYZ3(3,2,0)
     $ + 144* GYZ3(3,3,2)
     $ - 108* GYZ3(3,0,2)
     $ + 72* GYZ3(0,2,2)

     $ - 147* GYZ2(0,2)
     $ + 108* GYZ3(0,2,0)
     $ - 108* GYZ3(0,3,2)
     $ + 72* GYZ1(0)
     $ + 108* GYZ3(0,0,2)

     $ - 18* GYZ3(0,1,0)
     $ + 228* GYZ2(1,0)
     $ - 108* GYZ3(1,0,0)
     $ - 72* GYZ3(1,1,0)
     $)
     $ + (1.0d0)/(9.0d0)*(
     $   2*pi**2
     $ + 2*pi**2* HZr1(1)

     $ - 4*pi**2* GYZ1(2)
     $ + 2*pi**2* GYZ1(1)
     $ + 19* HZr1(0)
     $ + 12* HZr1(0)* GYZ2(2,2)
     $ - 29* HZr1(0)* GYZ1(2)

     $ + 6* HZr1(0)* GYZ2(2,0)
     $ - 18* HZr1(0)* GYZ2(3,2)
     $ + 6* HZr1(0)* GYZ2(0,2)
     $ - 3* HZr1(0)* GYZ1(0)

     $ - 3* HZr1(0)* GYZ2(1,0)
     $ - 9* HZr2(0,0)
     $ + 18* HZr2(0,0)* GYZ1(2)
     $ - 6* HZr3(0,0,1)
     $ - 35* HZr2(0,1)

     $ + 12* HZr2(0,1)* GYZ1(2)
     $ - 6* HZr2(0,1)* GYZ1(3)
     $ + 12* HZr2(0,1)* GYZ1(0)
     $ + 3* HZr3(0,1,0)
     $ - 12* HZr3(0,1,1)

     $ + 38* HZr1(1)
     $ + 24* HZr1(1)* GYZ2(2,3)
     $ - 12* HZr1(1)* GYZ1(2)
     $ - 12* HZr1(1)* GYZ2(2,0)

     $ + 24* HZr1(1)* GYZ2(3,2)
     $ - 24* HZr1(1)* GYZ2(3,3)
     $ - 64* HZr1(1)* GYZ1(3)
     $ + 18* HZr1(1)* GYZ2(3,0)

     $ - 12* HZr1(1)* GYZ2(0,2)
     $ + 18* HZr1(1)* GYZ2(0,3)
     $ + 29* HZr1(1)* GYZ1(0)
     $ - 18* HZr1(1)* GYZ2(0,0)
     $ - 3* HZr2(1,0)

     $ - 12* HZr2(1,0)* GYZ1(2)
     $ + 18* HZr2(1,0)* GYZ1(3)
     $ - 3* HZr2(1,0)* GYZ1(0)
     $ - 12* HZr3(1,0,1)
     $ + 12* HZr2(1,1)

     $ - 24* HZr2(1,1)* GYZ1(3)
     $ + 12* HZr2(1,1)* GYZ1(0)
     $ + 12* GYZ2(2,2)
     $ + 12* GYZ3(2,2,0)

     $ - 24* GYZ3(2,3,2)
     $ - 38* GYZ1(2)
     $ + 12* GYZ3(2,0,2)
     $ - 29* GYZ2(2,0)

     $ + 18* GYZ3(2,0,0)
     $ - 24* GYZ3(3,2,2)
     $ + 64* GYZ2(3,2)
     $ - 18* GYZ3(3,2,0)

     $ + 24* GYZ3(3,3,2)
     $ - 18* GYZ3(3,0,2)
     $ + 12* GYZ3(0,2,2)
     $ - 29* GYZ2(0,2)

     $ + 18* GYZ3(0,2,0)
     $ - 18* GYZ3(0,3,2)
     $ + 19* GYZ1(0)
     $ + 18* GYZ3(0,0,2)
     $ - 9* GYZ2(0,0)

     $ - 3* GYZ3(0,1,0)
     $ + 32* GYZ2(1,0)
     $ - 18* GYZ3(1,0,0)
     $ - 12* GYZ3(1,1,0)
     $)


      else if(i.eq.6) then
      fin20=T/108d0*(-17*pi**2-20*HZr1(0)+3*HZr1(0)*GYZ1(0)+15*HZr2(0,0)
     $ - 20*GYZ1(0)+15*GYZ2(0,0))

      else if(i.eq.7) then
      fin20=
     $(z)/(y**2)*(
     $ - 9* HZr1(0)* GYZ1(2)
     $ - 9* HZr1(1)* GYZ1(3)
     $ + 9* GYZ2(3,2)
     $)
     $ + (z**2)/(y**2)*(
     $   3* HZr1(0)* GYZ1(2)

     $ + 3* HZr1(1)* GYZ1(3)
     $ - 3* GYZ2(3,2)
     $)
     $ + (1.0d0)/(y**2)*(
     $   6* HZr1(0)* GYZ1(2)
     $ + 6* HZr1(1)* GYZ1(3)
     $ - 6* GYZ2(3,2)
     $)

     $ + (z*pi**2)/(9*y)*(
     $   12* HZr1(1)* GYZ1(2)
     $ - 12* GYZ2(2,2)
     $ - 2* GYZ1(2)
     $ + 12* GYZ2(2,0)
     $ - 2* GYZ1(1)
     $)

     $ + (z)/(3*y)*(
     $ - 72*z3* GYZ1(2)
     $ - 9* HZr1(0)
     $ + 12* HZr1(0)* GYZ2(2,2)
     $ - 18* HZr1(0)* GYZ1(2)

     $ - 4* HZr1(0)* GYZ2(2,0)
     $ - 8* HZr1(0)* GYZ2(3,2)
     $ - 4* HZr1(0)* GYZ2(1,0)
     $ + 24* HZr3(0,0,1)* GYZ1(2)

     $ + 24* HZr2(0,1)* GYZ2(2,2)
     $ - 8* HZr2(0,1)* GYZ1(2)
     $ - 24* HZr2(0,1)* GYZ2(2,0)

     $ - 24* HZr2(0,1)* GYZ2(3,2)
     $ - 24* HZr2(0,1)* GYZ2(3,3)
     $ + 4* HZr3(0,1,1)
     $ + 24* HZr3(0,1,1)* GYZ1(3)

     $ + 24* HZr1(1)* GYZ3(2,2,3)
     $ + 4* HZr1(1)* GYZ2(2,3)
     $ - 24* HZr1(1)* GYZ3(2,0,3)

     $ + 4* HZr1(1)* GYZ2(2,0)
     $ - 24* HZr1(1)* GYZ3(3,2,3)
     $ + 8* HZr1(1)* GYZ2(3,2)

     $ + 24* HZr1(1)* GYZ3(3,2,0)
     $ - 24* HZr1(1)* GYZ3(3,3,2)
     $ - 24* HZr1(1)* GYZ3(3,3,3)

     $ - 8* HZr1(1)* GYZ2(3,3)
     $ + 24* HZr1(1)* GYZ3(3,3,0)
     $ - 18* HZr1(1)* GYZ1(3)
     $ + 24* HZr1(1)* GYZ3(3,0,2)

     $ - 8* HZr1(1)* GYZ2(3,0)
     $ + 4* HZr1(1)* GYZ2(0,2)
     $ - 4* HZr1(1)* GYZ2(1,0)
     $ - 12* HZr2(1,0)* GYZ1(2)

     $ + 8* HZr2(1,0)* GYZ1(3)
     $ + 4* HZr3(1,0,1)
     $ - 24* HZr3(1,0,1)* GYZ1(2)
     $ + 24* HZr3(1,0,1)* GYZ1(3)

     $ + 24* HZr2(1,1)* GYZ2(3,3)
     $ - 8* HZr2(1,1)* GYZ1(3)
     $ - 24* HZr2(1,1)* GYZ2(3,0)
     $ - 4* HZr2(1,1)* GYZ1(0)

     $ + 4* HZr3(1,1,0)
     $ - 24* GYZ4(2,2,3,2)
     $ - 4* GYZ3(2,2,0)
     $ + 24* GYZ4(2,2,1,0)

     $ - 4* GYZ3(2,3,2)
     $ - 4* GYZ3(2,0,2)
     $ + 24* GYZ4(2,0,3,2)
     $ - 24* GYZ4(2,0,1,0)

     $ - 8* GYZ3(3,2,2)
     $ - 24* GYZ4(3,2,2,0)
     $ + 24* GYZ4(3,2,3,2)

     $ + 18* GYZ2(3,2)
     $ - 24* GYZ4(3,2,0,2)
     $ + 8* GYZ3(3,2,0)
     $ + 24* GYZ4(3,2,1,0)

     $ + 24* GYZ4(3,3,2,2)
     $ + 8* GYZ3(3,3,2)
     $ - 24* GYZ4(3,3,2,0)

     $ + 24* GYZ4(3,3,3,2)
     $ - 24* GYZ4(3,3,0,2)
     $ - 24* GYZ4(3,0,2,2)
     $ + 8* GYZ3(3,0,2)

     $ + 24* GYZ4(3,0,1,0)
     $ - 4* GYZ3(0,2,2)
     $ + 4* GYZ3(1,2,0)
     $ + 4* GYZ3(1,0,2)
     $ + 4* GYZ3(1,1,0)
     $)

     $ + (z**2)/(y)*(
     $   3* HZr1(0)* GYZ1(2)
     $ + 3* HZr1(1)* GYZ1(3)
     $ - 3* GYZ2(3,2)
     $)
     $ + (1.0d0)/(3*y*(1-y-z))*(
     $   (2*pi**2)/(3.0d0)* GYZ1(1)

     $ + 4* HZr1(0)* GYZ2(3,2)
     $ + 4* HZr1(0)* GYZ2(1,0)
     $ + 4* HZr1(1)* GYZ2(3,3)
     $ + 4* HZr1(1)* GYZ2(3,0)

     $ + 4* HZr1(1)* GYZ2(1,0)
     $ - 4* HZr2(1,0)* GYZ1(3)
     $ - 4* GYZ3(3,2,0)
     $ - 4* GYZ3(3,3,2)

     $ - 4* GYZ3(3,0,2)
     $ - 4* GYZ3(1,2,0)
     $ - 4* GYZ3(1,0,2)
     $ - 4* GYZ3(1,1,0)
     $)
     $ + (1.0d0)/(3*y*(1-z))*(
     $   (2*pi**2)/(3.0d0)* GYZ1(1)

     $ + 4* HZr1(0)* GYZ2(3,2)
     $ + 4* HZr1(0)* GYZ2(1,0)
     $ + 4* HZr1(1)* GYZ2(3,3)
     $ + 4* HZr1(1)* GYZ2(3,0)

     $ + 4* HZr1(1)* GYZ2(1,0)
     $ - 4* HZr2(1,0)* GYZ1(3)
     $ - 4* GYZ3(3,2,0)
     $ - 4* GYZ3(3,3,2)

     $ - 4* GYZ3(3,0,2)
     $ - 4* GYZ3(1,2,0)
     $ - 4* GYZ3(1,0,2)
     $ - 4* GYZ3(1,1,0)
     $)

     $ + (1.0d0)/(3*y*(y+z))*(
     $   8* HZr1(0)* GYZ2(2,2)
     $ - 8* HZr2(0,1)* GYZ1(2)
     $ + 8* HZr3(0,1,1)
     $ + 8* HZr1(1)* GYZ2(2,0)

     $ + 8* HZr1(1)* GYZ2(0,2)
     $ - 8* HZr2(1,0)* GYZ1(2)
     $ + 8* HZr3(1,0,1)
     $ - 8* HZr2(1,1)* GYZ1(0)
     $ + 8* HZr3(1,1,0)

     $ - 8* GYZ3(2,2,0)
     $ - 8* GYZ3(2,0,2)
     $ - 8* GYZ3(0,2,2)
     $)
     $ + (pi**2)/(9*y)*(
     $ - 6* HZr1(1)* GYZ1(2)

     $ + 6* GYZ2(2,2)
     $ + GYZ1(2)
     $ - 6* GYZ2(2,0)
     $ + GYZ1(1)
     $)
     $ + (1.0d0)/(3*y)*(
     $   36*z3* GYZ1(2)
     $ + 18* HZr1(0)

     $ - 12* HZr1(0)* GYZ2(2,2)
     $ + 11* HZr1(0)* GYZ1(2)
     $ - 4* HZr1(0)* GYZ2(2,0)
     $ - 2* HZr1(0)* GYZ2(3,2)

     $ + 6* HZr1(0)* GYZ2(0,2)
     $ + 2* HZr1(0)* GYZ2(1,0)
     $ - 12* HZr3(0,0,1)* GYZ1(2)
     $ - 12* HZr2(0,1)* GYZ2(2,2)

     $ + 4* HZr2(0,1)* GYZ1(2)
     $ + 12* HZr2(0,1)* GYZ2(2,0)
     $ + 12* HZr2(0,1)* GYZ2(3,2)

     $ + 12* HZr2(0,1)* GYZ2(3,3)
     $ + 6* HZr2(0,1)* GYZ1(3)
     $ - 8* HZr3(0,1,1)
     $ - 12* HZr3(0,1,1)* GYZ1(3)
     $ + 9* HZr1(1)

     $ - 12* HZr1(1)* GYZ3(2,2,3)
     $ - 8* HZr1(1)* GYZ2(2,3)
     $ + 12* HZr1(1)* GYZ3(2,0,3)

     $ - 8* HZr1(1)* GYZ2(2,0)
     $ + 12* HZr1(1)* GYZ3(3,2,3)
     $ - 4* HZr1(1)* GYZ2(3,2)

     $ - 12* HZr1(1)* GYZ3(3,2,0)
     $ + 12* HZr1(1)* GYZ3(3,3,2)
     $ + 12* HZr1(1)* GYZ3(3,3,3)

     $ + 4* HZr1(1)* GYZ2(3,3)
     $ - 12* HZr1(1)* GYZ3(3,3,0)
     $ + 11* HZr1(1)* GYZ1(3)
     $ - 12* HZr1(1)* GYZ3(3,0,2)

     $ - 2* HZr1(1)* GYZ2(3,0)
     $ - 8* HZr1(1)* GYZ2(0,2)
     $ + 6* HZr1(1)* GYZ2(0,3)
     $ + 2* HZr1(1)* GYZ2(1,0)

     $ + 6* HZr2(1,0)* GYZ1(2)
     $ + 2* HZr2(1,0)* GYZ1(3)
     $ - 8* HZr3(1,0,1)
     $ + 12* HZr3(1,0,1)* GYZ1(2)

     $ - 12* HZr3(1,0,1)* GYZ1(3)
     $ - 12* HZr2(1,1)* GYZ2(3,3)
     $ + 4* HZr2(1,1)* GYZ1(3)
     $ + 12* HZr2(1,1)* GYZ2(3,0)

     $ + 8* HZr2(1,1)* GYZ1(0)
     $ - 8* HZr3(1,1,0)
     $ + 12* GYZ4(2,2,3,2)
     $ + 8* GYZ3(2,2,0)

     $ - 12* GYZ4(2,2,1,0)
     $ + 8* GYZ3(2,3,2)
     $ - 9* GYZ1(2)
     $ + 8* GYZ3(2,0,2)

     $ - 12* GYZ4(2,0,3,2)
     $ + 12* GYZ4(2,0,1,0)
     $ + 4* GYZ3(3,2,2)

     $ + 12* GYZ4(3,2,2,0)
     $ - 12* GYZ4(3,2,3,2)
     $ - 11* GYZ2(3,2)

     $ + 12* GYZ4(3,2,0,2)
     $ + 2* GYZ3(3,2,0)
     $ - 12* GYZ4(3,2,1,0)

     $ - 12* GYZ4(3,3,2,2)
     $ - 4* GYZ3(3,3,2)
     $ + 12* GYZ4(3,3,2,0)

     $ - 12* GYZ4(3,3,3,2)
     $ + 12* GYZ4(3,3,0,2)
     $ + 12* GYZ4(3,0,2,2)

     $ + 2* GYZ3(3,0,2)
     $ - 12* GYZ4(3,0,1,0)
     $ + 8* GYZ3(0,2,2)
     $ - 6* GYZ3(0,3,2)

     $ - 2* GYZ3(1,2,0)
     $ - 2* GYZ3(1,0,2)
     $ - 2* GYZ3(1,1,0)
     $)
     $ + (z)/(3*(1-y-z))*(
     $ - (2*pi**2)/(3.0d0)* HZr1(1)
     $ - (2*pi**2)/(3.0d0)* GYZ1(1)

     $ - 8* HZr1(0)* GYZ2(3,2)
     $ + 4* HZr1(0)* GYZ2(0,2)
     $ - 4* HZr1(0)* GYZ2(1,0)
     $ + 8* HZr3(0,0,1)

     $ + 8* HZr2(0,1)* GYZ1(3)
     $ - 4* HZr2(0,1)* GYZ1(0)
     $ + 4* HZr3(0,1,0)
     $ - 8* HZr1(1)* GYZ2(3,0)
     $ - 4* HZr1(1)* GYZ2(1,0)

     $ + 8* HZr2(1,0)* GYZ1(3)
     $ - 4* HZr2(1,0)* GYZ1(0)
     $ - 4* HZr3(1,0,1)
     $ - 8* HZr3(1,1,0)
     $ + 8* GYZ3(3,2,0)

     $ + 8* GYZ3(3,0,2)
     $ - 4* GYZ3(0,1,0)
     $ + 4* GYZ3(1,2,0)
     $ + 4* GYZ3(1,0,2)
     $ + 4* GYZ3(1,1,0)
     $)

     $ + (z*pi**2)/(3*(1-y)**2)*(
     $ - 3
     $ - HZr1(0)
     $ - HZr1(1)
     $ + GYZ1(2)
     $)
     $ + (z)/((1-y)**2)*(
     $   6*z3
     $ - 2* HZr1(0)* GYZ2(0,2)

     $ + 2* HZr3(0,0,1)
     $ + 6* HZr2(0,1)
     $ - 2* HZr2(0,1)* GYZ1(2)
     $ - 2* HZr3(0,1,0)
     $ - 2* HZr1(1)* GYZ2(2,3)

     $ + 6* HZr1(1)* GYZ1(3)
     $ - 2* HZr1(1)* GYZ2(0,3)
     $ - 6* HZr2(1,0)
     $ + 2* HZr2(1,0)* GYZ1(2)
     $ + 2* HZr3(1,0,1)

     $ - 2* HZr3(1,1,0)
     $ + 2* GYZ3(2,3,2)
     $ - 6* GYZ2(3,2)
     $ + 2* GYZ3(0,3,2)
     $)

     $ + (z)/(1-y)*(
     $   3* HZr1(0)
     $ - 2* HZr1(0)* GYZ1(2)
     $ + 3* HZr1(1)
     $ - 2* HZr1(1)* GYZ1(3)
     $ - 3* GYZ1(2)
     $ + 2* GYZ2(3,2)
     $)

     $ + (z)/(3*(y+z)**3)*(
     $ - 8* HZr1(0)* GYZ2(2,2)
     $ - 24* HZr1(0)* GYZ1(2)
     $ + 24* HZr2(0,1)
     $ + 8* HZr2(0,1)* GYZ1(2)

     $ - 8* HZr3(0,1,1)
     $ - 8* HZr1(1)* GYZ2(2,0)
     $ - 8* HZr1(1)* GYZ2(0,2)
     $ - 24* HZr1(1)* GYZ1(0)
     $ + 24* HZr2(1,0)

     $ + 8* HZr2(1,0)* GYZ1(2)
     $ - 8* HZr3(1,0,1)
     $ + 8* HZr2(1,1)* GYZ1(0)
     $ - 8* HZr3(1,1,0)
     $ + 8* GYZ3(2,2,0)

     $ + 8* GYZ3(2,0,2)
     $ + 24* GYZ2(2,0)
     $ + 8* GYZ3(0,2,2)
     $ + 24* GYZ2(0,2)
     $)

     $ + (z)/(3*(y+z)**2)*(
     $ - 24* HZr1(0)
     $ + 4* HZr1(0)* GYZ2(2,2)
     $ + 16* HZr1(0)* GYZ1(2)
     $ - 16* HZr2(0,1)

     $ - 4* HZr2(0,1)* GYZ1(2)
     $ + 4* HZr3(0,1,1)
     $ + 4* HZr1(1)* GYZ2(2,0)
     $ + 4* HZr1(1)* GYZ2(0,2)

     $ + 16* HZr1(1)* GYZ1(0)
     $ - 16* HZr2(1,0)
     $ - 4* HZr2(1,0)* GYZ1(2)
     $ + 4* HZr3(1,0,1)
     $ - 4* HZr2(1,1)* GYZ1(0)

     $ + 4* HZr3(1,1,0)
     $ - 4* GYZ3(2,2,0)
     $ - 4* GYZ3(2,0,2)
     $ - 16* GYZ2(2,0)
     $ - 4* GYZ3(0,2,2)

     $ - 16* GYZ2(0,2)
     $ + 24* GYZ1(0)
     $)
     $ + (z)/(3*(y+z))*(
     $   8* HZr1(0)
     $ - 8* HZr1(0)* GYZ2(2,2)
     $ + 8* HZr2(0,1)* GYZ1(2)

     $ - 8* HZr3(0,1,1)
     $ - 8* HZr1(1)* GYZ2(2,0)
     $ - 8* HZr1(1)* GYZ2(0,2)
     $ + 8* HZr2(1,0)* GYZ1(2)
     $ - 8* HZr3(1,0,1)

     $ + 8* HZr2(1,1)* GYZ1(0)
     $ - 8* HZr3(1,1,0)
     $ + 8* GYZ3(2,2,0)
     $ + 8* GYZ3(2,0,2)
     $ + 8* GYZ3(0,2,2)

     $ - 8* GYZ1(0)
     $)
     $ + (z**2)/((1-y)**3)*(
     $  pi**2* HZr1(0)
     $ +pi**2* HZr1(1)
     $ -pi**2* GYZ1(2)
     $ - 18*z3
     $ + 6* HZr1(0)* GYZ2(0,2)

     $ - 6* HZr3(0,0,1)
     $ + 6* HZr2(0,1)* GYZ1(2)
     $ + 6* HZr3(0,1,0)
     $ + 6* HZr1(1)* GYZ2(2,3)
     $ + 6* HZr1(1)* GYZ2(0,3)

     $ - 6* HZr2(1,0)* GYZ1(2)
     $ - 6* HZr3(1,0,1)
     $ + 6* HZr3(1,1,0)
     $ - 6* GYZ3(2,3,2)
     $ - 6* GYZ3(0,3,2)
     $)

     $ + (z**2)/((1-y)**2)*(
     $   6* HZr1(0)* GYZ1(2)
     $ + 6* HZr1(1)* GYZ1(3)
     $ - 6* GYZ2(3,2)
     $)
     $ + (z**2)/(1-y)*(
     $   3* HZr1(0)* GYZ1(2)

     $ + 3* HZr1(1)* GYZ1(3)
     $ - 3* GYZ2(3,2)
     $)
     $ + (1.0d0)/(3*(1-y-z)*(1-y))*(
     $ - (2*pi**2)/(3.0d0)* HZr1(1)
     $ - 4* HZr1(0)* GYZ2(3,2)

     $ + 4* HZr1(0)* GYZ2(0,2)
     $ + 8* HZr3(0,0,1)
     $ + 8* HZr2(0,1)* GYZ1(3)
     $ - 4* HZr2(0,1)* GYZ1(0)
     $ + 4* HZr3(0,1,0)

     $ + 4* HZr1(1)* GYZ2(3,3)
     $ - 4* HZr1(1)* GYZ2(3,0)
     $ + 4* HZr2(1,0)* GYZ1(3)
     $ - 4* HZr2(1,0)* GYZ1(0)
     $ - 4* HZr3(1,0,1)

     $ - 8* HZr3(1,1,0)
     $ + 4* GYZ3(3,2,0)
     $ - 4* GYZ3(3,3,2)
     $ + 4* GYZ3(3,0,2)
     $ - 4* GYZ3(0,1,0)
     $)

     $ + (1.0d0)/(3*(1-y-z))*(
     $   (2*pi**2)/(3.0d0)
     $ + (pi**2)/(2.0d0)* HZr1(1)
     $ - (pi**2)/(2.0d0)* GYZ1(1)
     $ - 3* HZr1(0)* GYZ2(0,2)
     $ + 4* HZr1(0)* GYZ1(0)

     $ - 3* HZr1(0)* GYZ2(1,0)
     $ - 6* HZr3(0,0,1)
     $ - 6* HZr2(0,1)* GYZ1(3)
     $ + 3* HZr2(0,1)* GYZ1(0)
     $ - 3* HZr3(0,1,0)

     $ - 6* HZr1(1)* GYZ2(3,3)
     $ - 3* HZr1(1)* GYZ2(1,0)
     $ + 4* HZr2(1,0)
     $ + 3* HZr2(1,0)* GYZ1(0)
     $ + 3* HZr3(1,0,1)

     $ + 6* HZr3(1,1,0)
     $ + 6* GYZ3(3,3,2)
     $ + 3* GYZ3(0,1,0)
     $ + 3* GYZ3(1,2,0)
     $ + 3* GYZ3(1,0,2)

     $ - 4* GYZ2(1,0)
     $ + 3* GYZ3(1,1,0)
     $)
     $ + (1.0d0)/(1-y)*(
     $ - (pi**2)/(2.0d0)
     $ + (2*pi**2)/(3.0d0)* HZr1(0)
     $ + (2*pi**2)/(3.0d0)* HZr1(1)
     $ - (2*pi**2)/(3.0d0)* GYZ1(2)
     $ - 12*z3

     $ + 4* HZr1(0)* GYZ2(0,2)
     $ - HZr1(0)* GYZ1(0)
     $ - 4* HZr3(0,0,1)
     $ + HZr2(0,1)
     $ + 4* HZr2(0,1)* GYZ1(2)

     $ + 4* HZr3(0,1,0)
     $ + 4* HZr1(1)* GYZ2(2,3)
     $ + HZr1(1)* GYZ1(3)
     $ + 4* HZr1(1)* GYZ2(0,3)
     $ - HZr1(1)* GYZ1(0)

     $ - HZr2(1,0)
     $ - 4* HZr2(1,0)* GYZ1(2)
     $ - 4* HZr3(1,0,1)
     $ + 4* HZr3(1,1,0)
     $ - 4* GYZ3(2,3,2)

     $ + GYZ2(2,0)
     $ - GYZ2(3,2)
     $ + GYZ2(0,2)
     $ - 4* GYZ3(0,3,2)
     $ + 3* GYZ1(0)
     $ + 2* GYZ2(1,0)
     $)

     $ + (1.0d0)/(3*(y+z)**2)*(
     $ -pi**2* HZr1(1)
     $ +pi**2* GYZ1(2)
     $ - 6* HZr1(0)* GYZ2(2,2)
     $ - 6* HZr1(0)* GYZ1(2)

     $ + 6* HZr1(0)* GYZ2(2,0)
     $ - 6* HZr1(0)* GYZ2(3,2)
     $ + 6* HZr1(0)* GYZ2(0,2)
     $ - 6* HZr2(0,1)

     $ - 6* HZr2(0,1)* GYZ1(2)
     $ - 6* HZr2(0,1)* GYZ1(3)
     $ + 6* HZr3(0,1,0)
     $ + 6* HZr3(0,1,1)
     $ + 24* HZr1(1)

     $ - 12* HZr1(1)* GYZ2(2,3)
     $ - 8* HZr1(1)* GYZ1(2)
     $ + 6* HZr1(1)* GYZ2(2,0)
     $ - 12* HZr1(1)* GYZ2(3,2)

     $ - 12* HZr1(1)* GYZ2(3,3)
     $ - 12* HZr1(1)* GYZ1(3)
     $ + 6* HZr1(1)* GYZ2(3,0)
     $ + 6* HZr1(1)* GYZ2(0,2)

     $ + 6* HZr1(1)* GYZ2(0,3)
     $ + 6* HZr1(1)* GYZ1(0)
     $ + 6* HZr2(1,0)
     $ + 6* HZr2(1,0)* GYZ1(2)
     $ + 6* HZr2(1,0)* GYZ1(3)

     $ - 6* HZr2(1,0)* GYZ1(0)
     $ + 6* HZr3(1,0,1)
     $ + 8* HZr2(1,1)
     $ + 12* HZr2(1,1)* GYZ1(3)
     $ - 6* HZr2(1,1)* GYZ1(0)

     $ - 6* HZr3(1,1,0)
     $ + 8* GYZ2(2,2)
     $ - 6* GYZ3(2,2,0)
     $ + 12* GYZ3(2,3,2)
     $ - 24* GYZ1(2)

     $ - 6* GYZ3(2,0,2)
     $ - 6* GYZ2(2,0)
     $ + 12* GYZ3(3,2,2)
     $ + 12* GYZ2(3,2)

     $ - 6* GYZ3(3,2,0)
     $ + 12* GYZ3(3,3,2)
     $ - 6* GYZ3(3,0,2)
     $ - 6* GYZ3(0,2,2)

     $ - 6* GYZ2(0,2)
     $ - 6* GYZ3(0,3,2)
     $)
     $ + (1.0d0)/(3*(y+z))*(
     $  pi**2
     $ + (pi**2)/(2.0d0)* HZr1(1)
     $ - (pi**2)/(2.0d0)* GYZ1(2)
     $ - 12* HZr1(0)

     $ - 12* HZr1(0)* GYZ1(2)
     $ - 3* HZr1(0)* GYZ2(2,0)
     $ + 6* HZr1(0)* GYZ1(0)
     $ - 18* HZr2(0,1)
     $ - 3* HZr3(0,1,0)

     $ - 32* HZr1(1)
     $ - 8* HZr1(1)* GYZ1(2)
     $ - 30* HZr1(1)* GYZ1(3)
     $ + 12* HZr1(1)* GYZ1(0)
     $ + 12* HZr2(1,0)

     $ - 3* HZr2(1,0)* GYZ1(2)
     $ + 3* HZr2(1,0)* GYZ1(0)
     $ + 8* HZr2(1,1)
     $ + 3* HZr3(1,1,0)
     $ + 8* GYZ2(2,2)

     $ + 32* GYZ1(2)
     $ - 12* GYZ2(2,0)
     $ + 3* GYZ3(2,1,0)
     $ + 30* GYZ2(3,2)
     $ - 12* GYZ2(0,2)

     $ - 12* GYZ1(0)
     $ + 3* GYZ3(0,1,0)
     $)
      fin20 = fin20
     $ + (T)/(6.0d0)*(
     $ - (2*pi**2)/(3.0d0)* GYZ1(2)
     $ + (2*pi**2)/(3.0d0)* GYZ1(1)
     $ + 4* HZr1(0)* GYZ2(2,2)

     $ - 4* HZr1(0)* GYZ2(2,0)
     $ + 4* HZr1(0)* GYZ2(1,0)
     $ - 4* HZr3(0,1,1)
     $ - 9* HZr1(1)
     $ + 4* HZr1(1)* GYZ2(2,3)

     $ - 4* HZr1(1)* GYZ2(2,0)
     $ + 8* HZr1(1)* GYZ2(3,2)
     $ - 4* HZr1(1)* GYZ2(0,2)
     $ + 4* HZr1(1)* GYZ2(1,0)

     $ - 4* HZr2(1,0)* GYZ1(2)
     $ - 4* HZr3(1,0,1)
     $ - 8* HZr2(1,1)* GYZ1(3)
     $ + 4* HZr2(1,1)* GYZ1(0)
     $ - 4* HZr3(1,1,0)

     $ + 4* GYZ3(2,2,0)
     $ - 4* GYZ3(2,3,2)
     $ + 9* GYZ1(2)
     $ + 4* GYZ3(2,0,2)

     $ - 8* GYZ3(3,2,2)
     $ + 4* GYZ3(0,2,2)
     $ - 4* GYZ3(1,2,0)
     $ - 4* GYZ3(1,0,2)
     $ - 4* GYZ3(1,1,0)
     $)

     $ + (pi**2)/(18.0d0)*(
     $   25
     $ - 12* HZr1(0)
     $ + 6* HZr1(0)* GYZ1(2)
     $ - 6* HZr1(0)* GYZ1(1)
     $ - 23* HZr1(1)
     $ + 12* HZr1(1)* GYZ1(2)

     $ - 6* HZr1(1)* GYZ1(1)
     $ - 6* HZr2(1,0)
     $ - 6* HZr2(1,1)
     $ - 12* GYZ2(2,2)
     $ + 22* GYZ1(2)
     $ + 6* GYZ2(2,0)

     $ - 12* GYZ1(0)
     $ + 6* GYZ2(1,2)
     $ + GYZ1(1)
     $)
     $ + (1.0d0)/(3.0d0)*(
     $   72*z3 
     $ + 18*z3* HZr1(1)
     $ - 36*z3* GYZ1(2)
     $ + 18*z3* GYZ1(1)

     $ - 5* HZr1(0)
     $ + 11* HZr1(0)* GYZ2(2,2)
     $ + 5* HZr1(0)* GYZ1(2)
     $ + 6* HZr1(0)* GYZ3(2,0,2)

     $ + 4* HZr1(0)* GYZ2(2,0)
     $ - 6* HZr1(0)* GYZ3(3,2,2)
     $ + 9* HZr1(0)* GYZ2(3,2)

     $ - 6* HZr1(0)* GYZ3(3,3,2)
     $ + 6* HZr1(0)* GYZ3(3,0,2)
     $ + 6* HZr1(0)* GYZ3(0,2,2)

     $ - 22* HZr1(0)* GYZ2(0,2)
     $ + 6* HZr1(0)* GYZ3(0,3,2)
     $ + HZr1(0)* GYZ1(0)
     $ - 6* HZr1(0)* GYZ3(0,0,2)

     $ - 6* HZr1(0)* GYZ3(1,0,2)
     $ - 2* HZr1(0)* GYZ2(1,0)
     $ - 2* HZr3(0,0,1)
     $ + 6* HZr3(0,0,1)* GYZ1(1)
     $ - 11* HZr2(0,1)

     $ + 12* HZr2(0,1)* GYZ2(2,2)
     $ - 9* HZr2(0,1)* GYZ1(2)
     $ - 6* HZr2(0,1)* GYZ2(2,0)

     $ - 6* HZr2(0,1)* GYZ2(3,2)
     $ - 6* HZr2(0,1)* GYZ2(3,3)
     $ + 7* HZr2(0,1)* GYZ1(3)
     $ + 13* HZr2(0,1)* GYZ1(0)

     $ - 6* HZr2(0,1)* GYZ2(1,2)
     $ - 10* HZr3(0,1,0)
     $ + 6* HZr3(0,1,0)* GYZ1(2)
     $ - 6* HZr3(0,1,0)* GYZ1(1)

     $ - 11* HZr3(0,1,1)
     $ + 6* HZr3(0,1,1)* GYZ1(3)
     $ - 10* HZr1(1)
     $ + 12* HZr1(1)* GYZ3(2,2,3)

     $ + 2* HZr1(1)* GYZ2(2,3)
     $ + 16* HZr1(1)* GYZ1(2)
     $ - 11* HZr1(1)* GYZ2(2,0)

     $ - 12* HZr1(1)* GYZ3(3,2,3)
     $ + 22* HZr1(1)* GYZ2(3,2)
     $ + 6* HZr1(1)* GYZ3(3,2,0)

     $ - 12* HZr1(1)* GYZ3(3,3,2)
     $ - 12* HZr1(1)* GYZ3(3,3,3)
     $ + 16* HZr1(1)* GYZ2(3,3)

     $ + 6* HZr1(1)* GYZ3(3,3,0)
     $ - 6* HZr1(1)* GYZ1(3)
     $ + 6* HZr1(1)* GYZ3(3,0,2)
     $ + 6* HZr1(1)* GYZ3(3,0,3)

     $ - 9* HZr1(1)* GYZ2(3,0)
     $ + 6* HZr1(1)* GYZ3(0,2,3)
     $ - 11* HZr1(1)* GYZ2(0,2)
     $ + 6* HZr1(1)* GYZ3(0,3,2)

     $ + 6* HZr1(1)* GYZ3(0,3,3)
     $ - 9* HZr1(1)* GYZ2(0,3)
     $ - 5* HZr1(1)* GYZ1(0)
     $ - 6* HZr1(1)* GYZ3(0,0,3)

     $ - 6* HZr1(1)* GYZ3(1,2,3)
     $ - 6* HZr1(1)* GYZ3(1,0,3)
     $ - 2* HZr1(1)* GYZ2(1,0)
     $ + HZr2(1,0)

     $ - 6* HZr2(1,0)* GYZ2(2,2)
     $ + 13* HZr2(1,0)* GYZ1(2)
     $ + 6* HZr2(1,0)* GYZ2(3,3)
     $ - 9* HZr2(1,0)* GYZ1(3)

     $ - 6* HZr2(1,0)* GYZ2(0,3)
     $ - 2* HZr2(1,0)* GYZ1(0)
     $ + 6* HZr2(1,0)* GYZ2(1,2)
     $ + 6* HZr4(1,0,0,1)
     $ + 11* HZr3(1,0,1)

     $ - 12* HZr3(1,0,1)* GYZ1(2)
     $ + 6* HZr3(1,0,1)* GYZ1(3)
     $ + 6* HZr3(1,0,1)* GYZ1(1)
     $ - 6* HZr4(1,0,1,0)

     $ - 16* HZr2(1,1)
     $ + 12* HZr2(1,1)* GYZ2(3,3)
     $ - 22* HZr2(1,1)* GYZ1(3)
     $ - 6* HZr2(1,1)* GYZ2(3,0)

     $ - 6* HZr2(1,1)* GYZ2(0,3)
     $ + 11* HZr2(1,1)* GYZ1(0)
     $ - 12* HZr3(1,1,0)
     $ + 6* HZr3(1,1,0)* GYZ1(2)

     $ - 6* HZr3(1,1,0)* GYZ1(1)
     $ + 6* HZr4(1,1,0,1)
     $ - 6* HZr4(1,1,1,0)
     $ - 12* GYZ4(2,2,3,2)

     $ - 16* GYZ2(2,2)
     $ + 11* GYZ3(2,2,0)
     $ + 6* GYZ4(2,2,1,0)
     $ - 2* GYZ3(2,3,2)

     $ + 10* GYZ1(2)
     $ + 11* GYZ3(2,0,2)
     $ + 5* GYZ2(2,0)
     $ - 6* GYZ4(2,0,1,0)
     $ - 24* GYZ3(2,1,0)

     $ - 22* GYZ3(3,2,2)
     $ - 6* GYZ4(3,2,2,0)
     $ + 12* GYZ4(3,2,3,2)

     $ + 6* GYZ2(3,2)
     $ - 6* GYZ4(3,2,0,2)
     $ + 9* GYZ3(3,2,0)
     $ + 6* GYZ4(3,2,1,0)

     $ + 12* GYZ4(3,3,2,2)
     $ - 16* GYZ3(3,3,2)
     $ - 6* GYZ4(3,3,2,0)

     $ + 12* GYZ4(3,3,3,2)
     $ - 6* GYZ4(3,3,0,2)
     $ - 6* GYZ4(3,0,2,2)

     $ + 9* GYZ3(3,0,2)
     $ - 6* GYZ4(3,0,3,2)
     $ + 6* GYZ4(3,0,1,0)
     $ + 11* GYZ3(0,2,2)

     $ - 6* GYZ4(0,2,3,2)
     $ + 5* GYZ2(0,2)
     $ - 6* GYZ4(0,3,2,2)
     $ + 9* GYZ3(0,3,2)

     $ - 6* GYZ4(0,3,3,2)
     $ - 5* GYZ1(0)
     $ + 6* GYZ4(0,0,3,2)
     $ + GYZ3(0,1,0)
     $ + 6* GYZ4(1,2,3,2)

     $ + 2* GYZ3(1,2,0)
     $ + 2* GYZ3(1,0,2)
     $ + 6* GYZ4(1,0,3,2)
     $ - 6* GYZ2(1,0)
     $ - GYZ3(1,1,0)
     $)


!      else
      end if
      end function


************************************************************************
****      1x1 two-loop finite part coded from hep-ph/0112081v1      ****
************************************************************************
      real*8 function fin11(i,y,z)
      implicit none

      integer i ! A, B, C, D, E, F - (1-6)
      real*8 y, z, T
      real*8 GYZ1(0:3),GYZ2(0:3,0:3),GYZ3(0:3,0:3,0:3),
     $ GYZ4(0:3,0:3,0:3,0:3), HZr1(0:1),HZr2(0:1,0:1),HZr3(0:1,0:1,0:1),
     $ HZr4(0:1,0:1,0:1,0:1)
      real*8 NC, nf, TR, pi, z3, eg
      common /piez/ pi, z3
      common /const/ NC, nf, TR, eg

      T= y/z+z/y+2*(1/(y*z)-1/y-1/z)

      call tdhpl(y,z,4,GYZ1,GYZ2,GYZ3,GYZ4,HZr1,HZr2,HZr3,HZr4)


      if(i.eq.1) then
      fin11=
     $(1.0d0)/(2.0d0) +  (z)/(4*y)
     $+ (1.0d0)/(6*y)*(
     $ - pi**2
     $ - 24
     $ - 10* HZr1(0)
     $ - 6* HZr1(0)*GYZ1(0)
     $ - 6* HZr2(1,0)
     $ - 10* GYZ1(0)
     $ + 6* GYZ2(1,0)
     $)
     $+ (z)/(6*(1-y)**2)*(
     $  pi**2* GYZ1(0)
     $ + 10* HZr1(0)*GYZ1(0)
     $ + 12* HZr1(0)*GYZ2(0,0)
     $ + 6* HZr2(1,0)*GYZ1(0)
     $ + 21* GYZ1(0)
     $ + 23* GYZ2(0,0)
     $ - 6* GYZ3(0,1,0)
     $ - 12* GYZ3(1,0,0)
     $)
     $+ (z)/(6*(1-y))*(
     $  pi**2 
     $ + 3*pi**2* GYZ1(0)
     $ + 21
     $ + 10* HZr1(0)
     $ + 36*HZr1(0)*GYZ1(0)
     $ + 36*HZr1(0)*GYZ2(0,0)
     $ + 6* HZr2(1,0)
     $ + 18* HZr2(1,0)*GYZ1(0)
     $ + 73* GYZ1(0)
     $ + 33* GYZ2(0,0)
     $ - 18*GYZ3(0,1,0)
     $ - 6* GYZ2(1,0)
     $ - 36* GYZ3(1,0,0)
     $)
     $+ (1.0d0)/(3*(1-y))*(
     $ - 2*pi**2* GYZ1(0)
     $ - 20*HZr1(0)*GYZ1(0)
     $ - 24*HZr1(0)*GYZ2(0,0)
     $ - 12* HZr2(1,0)*GYZ1(0)
     $ - 42* GYZ1(0)
     $ - 28* GYZ2(0,0)
     $ + 12* GYZ3(0,1,0)
     $ + 24* GYZ3(1,0,0)
     $)
     $+ (T*pi**2 )/(72.0d0)*(
     $   pi**2
     $ + 169
     $ + 20* HZr1(0)
     $ + 12* HZr1(0)*GYZ1(0)
     $ + 12* HZr2(1,0)
     $ + 20* GYZ1(0)
     $ - 12* GYZ2(1,0)
     $)
     $+ (T)/(9.0d0)*(
     $ + 72
     $ + 60* HZr1(0)
     $ + 61* HZr1(0)*GYZ1(0)
     $ + 30*HZr1(0)*GYZ2(0,0)
     $ - 9* HZr1(0)*GYZ3(0,1,0)
     $ - 15* HZr1(0)*GYZ2(1,0)
     $ - 18* HZr1(0)*GYZ3(1,0,0)
     $ + 25*HZr2(0,0)
     $ + 30*HZr2(0,0)*GYZ1(0)
     $ + 18* HZr2(0,0)*GYZ2(0,0)
     $ + 15*HZr3(0,1,0)
     $ + 9* HZr3(0,1,0)*GYZ1(0)
     $ + 36* HZr2(1,0)
     $ + 15* HZr2(1,0)*GYZ1(0)
     $ - 9* HZr2(1,0)*GYZ2(1,0)
     $ + 30*HZr3(1,0,0)
     $ + 18*HZr3(1,0,0)*GYZ1(0)
     $ + 9* HZr4(1,0,1,0)
     $ + 18* HZr4(1,1,0,0)
     $ + 60*GYZ1(0)
     $ + 25*GYZ2(0,0)
     $ - 15*GYZ3(0,1,0)
     $ - 36* GYZ2(1,0)
     $ - 30*GYZ3(1,0,0)
     $ + 9* GYZ4(1,0,1,0)
     $ + 18*GYZ4(1,1,0,0)
     $)


      else if(i.eq.2) then
      fin11=
     $+ (z)/(6*y)*(
     $ - 2*pi**2* HZr1(0)*GYZ1(2)
     $ - 2*pi**2* HZr1(1)*GYZ1(3)
     $ + 2*pi**2* GYZ2(3,2)
     $ + 3
     $ - 42* HZr1(0)*GYZ1(2)

     $ - 20* HZr1(0)*GYZ2(2,0)
     $ + 12* HZr1(0)*GYZ3(2,1,0)
     $ + 2* HZr1(0)*GYZ2(3,2)

     $ + 12* HZr1(0)*GYZ3(3,2,0)
     $ + 12* HZr1(0)*GYZ3(3,0,2)
     $ - 20* HZr1(0)*GYZ2(0,2)

     $ + 12* HZr1(0)*GYZ3(0,3,2)
     $ + 12* HZr1(0)*GYZ3(1,2,0)
     $ + 12* HZr1(0)*GYZ3(1,0,2)

     $ - 4* HZr2(0,0)*GYZ1(2)
     $ - 24* HZr2(0,0)*GYZ2(2,0)
     $ - 24* HZr2(0,0)*GYZ2(0,2)
     $ - 2* HZr2(0,1)*GYZ1(3)

     $ - 12* HZr2(0,1)*GYZ2(3,0)
     $ - 12* HZr2(0,1)*GYZ2(0,3)
     $ - 12* HZr3(0,1,0)*GYZ1(2)
     $ - 42* HZr1(1)*GYZ1(3)

     $ - 20* HZr1(1)*GYZ2(3,0)
     $ + 12* HZr1(1)*GYZ3(3,1,0)
     $ - 20* HZr1(1)*GYZ2(0,3)
     $ + 12* HZr1(1)*GYZ3(1,3,0)

     $ + 12* HZr1(1)*GYZ3(1,0,3)
     $ + 12* HZr2(1,0)*GYZ2(3,2)
     $ - 2* HZr2(1,0)*GYZ1(3)
     $ - 12* HZr2(1,0)*GYZ2(3,0)

     $ - 12* HZr2(1,0)*GYZ2(0,3)
     $ - 24* HZr3(1,0,0)*GYZ1(2)
     $ - 12* HZr3(1,0,1)*GYZ1(3)
     $ - 24* HZr3(1,1,0)*GYZ1(3)

     $ + 42* GYZ2(3,2)
     $ + 20* GYZ3(3,2,0)
     $ - 12* GYZ4(3,2,1,0)
     $ + 20* GYZ3(3,0,2)

     $ - 12* GYZ4(3,1,2,0)
     $ - 12* GYZ4(3,1,0,2)
     $ + 20* GYZ3(0,3,2)
     $ - 12* GYZ4(1,3,2,0)

     $ - 12* GYZ4(1,3,0,2)
     $ - 12* GYZ4(1,0,3,2)
     $)
     $+ (1.0d0)/(6*y)*(
     $ - pi**2 
     $ + 4*pi**2* HZr1(0)*GYZ1(2)

     $ + 4*pi**2* HZr1(1)*GYZ1(3)
     $ - 4*pi**2* GYZ2(3,2)
     $ - 10* HZr1(0)
     $ + 90* HZr1(0)*GYZ1(2)
     $ + 40* HZr1(0)*GYZ2(2,0)

     $ - 24* HZr1(0)*GYZ3(2,1,0)
     $ - 16* HZr1(0)*GYZ2(3,2)
     $ - 24* HZr1(0)*GYZ3(3,2,0)

     $ - 24* HZr1(0)*GYZ3(3,0,2)
     $ + 40* HZr1(0)*GYZ2(0,2)
     $ - 24* HZr1(0)*GYZ3(0,3,2)
     $ - 6* HZr1(0)*GYZ1(0)

     $ - 24* HZr1(0)*GYZ3(1,2,0)
     $ - 24* HZr1(0)*GYZ3(1,0,2)
     $ + 32* HZr2(0,0)*GYZ1(2)

     $ + 48* HZr2(0,0)*GYZ2(2,0)
     $ + 48* HZr2(0,0)*GYZ2(0,2)
     $ + 6* HZr2(0,1)
     $ + 16* HZr2(0,1)*GYZ1(3)

     $ + 24* HZr2(0,1)*GYZ2(3,0)
     $ + 24* HZr2(0,1)*GYZ2(0,3)
     $ + 24* HZr3(0,1,0)*GYZ1(2)
     $ - 9* HZr1(1)

     $ + 96* HZr1(1)*GYZ1(3)
     $ + 40* HZr1(1)*GYZ2(3,0)
     $ - 24* HZr1(1)*GYZ3(3,1,0)
     $ + 40* HZr1(1)*GYZ2(0,3)

     $ - 6* HZr1(1)*GYZ1(0)
     $ - 24* HZr1(1)*GYZ3(1,3,0)
     $ - 24* HZr1(1)*GYZ3(1,0,3)
     $ - 6* HZr2(1,0)

     $ - 24* HZr2(1,0)*GYZ2(3,2)
     $ + 16* HZr2(1,0)*GYZ1(3)
     $ + 24* HZr2(1,0)*GYZ2(3,0)
     $ + 24* HZr2(1,0)*GYZ2(0,3)

     $ + 48* HZr3(1,0,0)*GYZ1(2)
     $ + 24* HZr3(1,0,1)*GYZ1(3)
     $ + 48* HZr3(1,1,0)*GYZ1(3)
     $ + 9* GYZ1(2)

     $ + 6* GYZ2(2,0)
     $ - 96* GYZ2(3,2)
     $ - 40* GYZ3(3,2,0)
     $ + 24* GYZ4(3,2,1,0)

     $ - 40* GYZ3(3,0,2)
     $ + 24* GYZ4(3,1,2,0)
     $ + 24* GYZ4(3,1,0,2)
     $ + 6* GYZ2(0,2)

     $ - 40* GYZ3(0,3,2)
     $ - 10* GYZ1(0)
     $ + 24* GYZ4(1,3,2,0)
     $ + 24* GYZ4(1,3,0,2)

     $ + 24* GYZ4(1,0,3,2)
     $)
     $+ (z)/(6*(1-y)**2)*(
     $ - pi**2* GYZ1(0)
     $ - 6* HZr1(0)*GYZ2(2,0)
     $ - 6* HZr1(0)*GYZ2(0,2)

     $ - 10* HZr1(0)*GYZ1(0)
     $ - 12* HZr1(0)*GYZ2(0,0)
     $ - 6* HZr2(0,1)*GYZ1(0)
     $ - 12* HZr1(1)*GYZ2(3,0)

     $ - 12* HZr1(1)*GYZ2(0,3)
     $ + 9* HZr1(1)*GYZ1(0)
     $ + 12* HZr1(1)*GYZ2(0,0)
     $ - 6* HZr2(1,0)*GYZ1(0)
     $ - 9* GYZ2(2,0)

     $ - 12* GYZ3(2,0,0)
     $ + 12* GYZ3(3,2,0)
     $ + 12* GYZ3(3,0,2)
     $ - 9* GYZ2(0,2)
     $ - 12* GYZ3(0,2,0)

     $ + 12* GYZ3(0,3,2)
     $ - 42* GYZ1(0)
     $ - 12* GYZ3(0,0,2)
     $ - 26* GYZ2(0,0)
     $ + 12* GYZ3(0,1,0)
     $ + 24* GYZ3(1,0,0)
     $)

     $+ (z)/(6*(1-y))*(
     $ - pi**2 
     $ - 3*pi**2* GYZ1(0)
     $ - 42
     $ - 10* HZr1(0)
     $ - 6* HZr1(0)*GYZ1(2)
     $ - 18* HZr1(0)*GYZ2(2,0)

     $ - 18* HZr1(0)*GYZ2(0,2)
     $ - 36* HZr1(0)*GYZ1(0)
     $ - 36* HZr1(0)*GYZ2(0,0)
     $ - 6* HZr2(0,1)
     $ - 18* HZr2(0,1)*GYZ1(0)

     $ + 9* HZr1(1)
     $ - 12* HZr1(1)*GYZ1(3)
     $ - 36* HZr1(1)*GYZ2(3,0)
     $ - 36* HZr1(1)*GYZ2(0,3)
     $ + 33* HZr1(1)*GYZ1(0)

     $ + 36* HZr1(1)*GYZ2(0,0)
     $ - 6* HZr2(1,0)
     $ - 18* HZr2(1,0)*GYZ1(0)
     $ - 9* GYZ1(2)
     $ - 33* GYZ2(2,0)

     $ - 36* GYZ3(2,0,0)
     $ + 12* GYZ2(3,2)
     $ + 36* GYZ3(3,2,0)
     $ + 36* GYZ3(3,0,2)
     $ - 33* GYZ2(0,2)

     $ - 36* GYZ3(0,2,0)
     $ + 36* GYZ3(0,3,2)
     $ - 136* GYZ1(0)
     $ - 36* GYZ3(0,0,2)
     $ - 6* GYZ2(0,0)

     $ + 36* GYZ3(0,1,0)
     $ + 12* GYZ2(1,0)
     $ + 72* GYZ3(1,0,0)
     $)
     $+ (1.0d0)/((1-y)**2*(y+z)**2)*(
     $ - HZr1(1)*GYZ1(0)
     $ + GYZ2(2,0)
 
     $+ GYZ2(0,2)
     $)
     $+ (1.0d0)/((1-y)**2*(y+z))*(
     $   2* HZr1(1)*GYZ1(0)
     $ - 2* GYZ2(2,0)
     $ - 2* GYZ2(0,2)
     $ + GYZ1(0)
     $)

     $+ (1.0d0)/((1-y)**2)*(
     $ - HZr1(1)*GYZ1(0)
     $ + GYZ2(2,0)
     $ + GYZ2(0,2)
     $ - GYZ1(0)
     $)
     $+ (1.0d0)/((1-y)*(y+z)**2)*(
     $ - HZr1(1)

     $ - 2* HZr1(1)*GYZ1(0)
     $ + GYZ1(2)
     $ + 2* GYZ2(2,0)
     $ + 2* GYZ2(0,2)
     $)
     $+ (1.0d0)/((1-y)*(y+z))*(
     $   1
     $ + 2* HZr1(1)

     $ + 2* HZr1(1)*GYZ1(0)
     $ - 2* GYZ1(2)
     $ - 2* GYZ2(2,0)
     $ - 2* GYZ2(0,2)
     $ + 2* GYZ1(0)
     $)
     $+ (1.0d0)/(3*(1-y))*(
     $   pi**2* GYZ1(0)

     $ - 3
     $ + 12* HZr1(0)*GYZ2(2,0)
     $ + 12* HZr1(0)*GYZ2(0,2)
     $ + 10* HZr1(0)*GYZ1(0)
     $ + 12* HZr1(0)*GYZ2(0,0)

     $ + 6* HZr2(0,1)*GYZ1(0)
     $ - 3* HZr1(1)
     $ + 18* HZr1(1)*GYZ2(3,0)
     $ + 18* HZr1(1)*GYZ2(0,3)
     $ - 18* HZr1(1)*GYZ1(0)

     $ - 12* HZr1(1)*GYZ2(0,0)
     $ + 6* HZr2(1,0)*GYZ1(0)
     $ + 3* GYZ1(2)
     $ + 18* GYZ2(2,0)
     $ + 12* GYZ3(2,0,0)

     $ - 18* GYZ3(3,2,0)
     $ - 18* GYZ3(3,0,2)
     $ + 18* GYZ2(0,2)
     $ + 12* GYZ3(0,2,0)

     $ - 18* GYZ3(0,3,2)
     $ + 66* GYZ1(0)
     $ + 12* GYZ3(0,0,2)
     $ + 14* GYZ2(0,0)
     $ - 12* GYZ3(0,1,0)
     $ - 24* GYZ3(1,0,0)
     $)

     $+ (1.0d0)/(6*(y+z)**2)*(
     $ - 2*pi**2* HZr1(1)
     $ + 2*pi**2* GYZ1(2)
     $ + 11* HZr1(0)*GYZ1(2)
     $ + 12* HZr1(0)*GYZ2(2,0)

     $ + 12* HZr1(0)*GYZ2(0,2)
     $ - 11* HZr2(0,1)
     $ - 12* HZr2(0,1)*GYZ1(0)
     $ - 42* HZr1(1)
     $ - 11* HZr1(1)*GYZ1(0)

     $ + 12* HZr1(1)*GYZ2(1,0)
     $ - 11* HZr2(1,0)
     $ + 12* HZr2(1,0)*GYZ1(2)
     $ - 12* HZr2(1,0)*GYZ1(0)
     $ - 12* HZr3(1,0,1)

     $ - 24* HZr3(1,1,0)
     $ + 42* GYZ1(2)
     $ + 11* GYZ2(2,0)
     $ - 12* GYZ3(2,1,0)
     $ + 11* GYZ2(0,2)

     $ - 12* GYZ3(1,2,0)
     $ - 12* GYZ3(1,0,2)
     $)
     $+ (1.0d0)/(6*(y+z))*(
     $   2*pi**2
     $ + 42
     $ + 11* HZr1(0)
     $ + 12* HZr1(0)*GYZ1(0)

     $ + 12* HZr2(1,0)
     $ + 11* GYZ1(0)
     $ - 12* GYZ2(1,0)
     $)
     $+ (T*pi**2)/(12.0d0)*(
     $ - 8
     $ - 2* HZr1(0)*GYZ1(2)
     $ - 2* HZr2(0,1)
     $ + 3* HZr1(1)

     $ - 4* HZr1(1)*GYZ1(3)
     $ + 2* HZr1(1)*GYZ1(0)
     $ - 3* GYZ1(2)
     $ - 2* GYZ2(2,0)
     $ + 4* GYZ2(3,2)

     $ - 2* GYZ2(0,2)
     $ + 2* GYZ2(1,0)
     $)
     $+ (T)/(6.0d0)*(
     $ - 96
     $ - 40* HZr1(0)
     $ - 39* HZr1(0)*GYZ1(2)
     $ - 29* HZr1(0)*GYZ2(2,0)

     $ - 12* HZr1(0)*GYZ3(2,0,0)
     $ + 6* HZr1(0)*GYZ3(2,1,0)
     $ + 20* HZr1(0)*GYZ2(3,2)

     $ + 12* HZr1(0)*GYZ3(3,2,0)
     $ + 12* HZr1(0)*GYZ3(3,0,2)
     $ - 29* HZr1(0)*GYZ2(0,2)

     $ - 12* HZr1(0)*GYZ3(0,2,0)
     $ + 12* HZr1(0)*GYZ3(0,3,2)
     $ - 24* HZr1(0)*GYZ1(0)
     $ - 12* HZr1(0)*GYZ3(0,0,2)

     $ + 6* HZr1(0)*GYZ3(0,1,0)
     $ + 6* HZr1(0)*GYZ3(1,2,0)
     $ + 6* HZr1(0)*GYZ3(1,0,2)
     $ + 10* HZr1(0)*GYZ2(1,0)

     $ + 12* HZr1(0)*GYZ3(1,0,0)
     $ - 20* HZr2(0,0)*GYZ1(2)
     $ - 12* HZr2(0,0)*GYZ2(2,0)
     $ - 12* HZr2(0,0)*GYZ2(0,2)

     $ - 20* HZr3(0,0,1)
     $ - 12* HZr3(0,0,1)*GYZ1(0)
     $ - 9* HZr2(0,1)
     $ - 20* HZr2(0,1)*GYZ1(3)
     $ - 12* HZr2(0,1)*GYZ2(3,0)

     $ - 12* HZr2(0,1)*GYZ2(0,3)
     $ + 9* HZr2(0,1)*GYZ1(0)
     $ + 12* HZr2(0,1)*GYZ2(0,0)
     $ + 6* HZr2(0,1)*GYZ2(1,0)

     $ - 10* HZr3(0,1,0)
     $ - 6* HZr3(0,1,0)*GYZ1(2)
     $ - 6* HZr3(0,1,0)*GYZ1(0)
     $ - 6* HZr4(0,1,0,1)
     $ - 12* HZr4(0,1,1,0)

     $ + 36* HZr1(1)
     $ - 48* HZr1(1)*GYZ1(3)
     $ - 20* HZr1(1)*GYZ2(3,0)
     $ + 12* HZr1(1)*GYZ3(3,1,0)
     $ - 20* HZr1(1)*GYZ2(0,3)

     $ + 39* HZr1(1)*GYZ1(0)
     $ + 20* HZr1(1)*GYZ2(0,0)
     $ - 6* HZr1(1)*GYZ3(0,1,0)
     $ + 12* HZr1(1)*GYZ3(1,3,0)

     $ + 12* HZr1(1)*GYZ3(1,0,3)
     $ - 9* HZr1(1)*GYZ2(1,0)
     $ - 12* HZr1(1)*GYZ3(1,0,0)
     $ - 9* HZr2(1,0)

     $ - 9* HZr2(1,0)*GYZ1(2)
     $ - 6* HZr2(1,0)*GYZ2(2,0)
     $ + 12* HZr2(1,0)*GYZ2(3,2)
     $ - 20* HZr2(1,0)*GYZ1(3)

     $ - 12* HZr2(1,0)*GYZ2(3,0)
     $ - 6* HZr2(1,0)*GYZ2(0,2)
     $ - 12* HZr2(1,0)*GYZ2(0,3)
     $ + 19* HZr2(1,0)*GYZ1(0)

     $ + 12* HZr2(1,0)*GYZ2(0,0)
     $ + 6* HZr2(1,0)*GYZ2(1,0)
     $ - 12* HZr3(1,0,0)*GYZ1(2)
     $ - 12* HZr4(1,0,0,1)

     $ + 9* HZr3(1,0,1)
     $ - 12* HZr3(1,0,1)*GYZ1(3)
     $ + 6* HZr3(1,0,1)*GYZ1(0)
     $ - 6* HZr4(1,0,1,0)
     $ + 18* HZr3(1,1,0)

     $ - 24* HZr3(1,1,0)*GYZ1(3)
     $ + 12* HZr3(1,1,0)*GYZ1(0)
     $ - 36* GYZ1(2)
     $ - 39* GYZ2(2,0)
     $ - 20* GYZ3(2,0,0)

     $ + 6* GYZ4(2,0,1,0)
     $ + 9* GYZ3(2,1,0)
     $ + 12* GYZ4(2,1,0,0)
     $ + 48* GYZ2(3,2)

     $ + 20* GYZ3(3,2,0)
     $ - 12* GYZ4(3,2,1,0)
     $ + 20* GYZ3(3,0,2)
     $ - 12* GYZ4(3,1,2,0)

     $ - 12* GYZ4(3,1,0,2)
     $ - 39* GYZ2(0,2)
     $ - 20* GYZ3(0,2,0)
     $ + 6* GYZ4(0,2,1,0)

     $ + 20* GYZ3(0,3,2)
     $ - 40* GYZ1(0)
     $ - 20* GYZ3(0,0,2)
     $ + 6* GYZ4(0,1,2,0)
     $ + 6* GYZ4(0,1,0,2)

     $ + 10* GYZ3(0,1,0)
     $ + 9* GYZ3(1,2,0)
     $ + 12* GYZ4(1,2,0,0)
     $ - 12* GYZ4(1,3,2,0)

     $ - 12* GYZ4(1,3,0,2)
     $ + 9* GYZ3(1,0,2)
     $ + 12* GYZ4(1,0,2,0)
     $ - 12* GYZ4(1,0,3,2)

     $ + 48* GYZ2(1,0)
     $ + 12* GYZ4(1,0,0,2)
     $ + 20* GYZ3(1,0,0)
     $ - 12* GYZ4(1,0,1,0)
     $ - 24* GYZ4(1,1,0,0)
     $)

     $+ (pi**2)/(6.0d0)*(
     $   HZr1(0)
     $ - 2* HZr1(0)*GYZ1(2)
     $ - 2* HZr2(0,1)
     $ + 2* HZr1(1)
     $ - 4* HZr1(1)*GYZ1(3)
     $ + 2* HZr1(1)*GYZ1(0)

     $ - 2* GYZ1(2)
     $ - 2* GYZ2(2,0)
     $ + 4* GYZ2(3,2)
     $ - 2* GYZ2(0,2)
     $ + GYZ1(0)
     $ + 2* GYZ2(1,0)
     $)

     $+ (1.0d0)/(6.0d0)*(
     $   21* HZr1(0)
     $ - 53* HZr1(0)*GYZ1(2)
     $ - 52* HZr1(0)*GYZ2(2,0)
     $ - 24* HZr1(0)*GYZ3(2,0,0)

     $ + 12* HZr1(0)*GYZ3(2,1,0)
     $ + 22* HZr1(0)*GYZ2(3,2)
     $ + 24* HZr1(0)*GYZ3(3,2,0)

     $ + 24* HZr1(0)*GYZ3(3,0,2)
     $ - 52* HZr1(0)*GYZ2(0,2)
     $ - 24* HZr1(0)*GYZ3(0,2,0)

     $ + 24* HZr1(0)*GYZ3(0,3,2)
     $ + 20* HZr1(0)*GYZ1(0)
     $ - 24* HZr1(0)*GYZ3(0,0,2)
     $ + 12* HZr1(0)*GYZ2(0,0)

     $ + 12* HZr1(0)*GYZ3(0,1,0)
     $ + 12* HZr1(0)*GYZ3(1,2,0)
     $ + 12* HZr1(0)*GYZ3(1,0,2)
     $ + 14* HZr1(0)*GYZ2(1,0)

     $ + 24* HZr1(0)*GYZ3(1,0,0)
     $ + 2* HZr2(0,0)
     $ - 4* HZr2(0,0)*GYZ1(2)
     $ - 24* HZr2(0,0)*GYZ2(2,0)

     $ - 24* HZr2(0,0)*GYZ2(0,2)
     $ + 12* HZr2(0,0)*GYZ1(0)
     $ - 40* HZr3(0,0,1)
     $ - 24* HZr3(0,0,1)*GYZ1(0)
     $ - 31* HZr2(0,1)

     $ - 22* HZr2(0,1)*GYZ1(3)
     $ - 24* HZr2(0,1)*GYZ2(3,0)
     $ - 24* HZr2(0,1)*GYZ2(0,3)
     $ + 30* HZr2(0,1)*GYZ1(0)

     $ + 24* HZr2(0,1)*GYZ2(0,0)
     $ + 12* HZr2(0,1)*GYZ2(1,0)
     $ - 14* HZr3(0,1,0)
     $ - 12* HZr3(0,1,0)*GYZ1(2)

     $ - 12* HZr3(0,1,0)*GYZ1(0)
     $ - 12* HZr4(0,1,0,1)
     $ - 24* HZr4(0,1,1,0)
     $ + 42* HZr1(1)
     $ - 84* HZr1(1)*GYZ1(3)

     $ - 22* HZr1(1)*GYZ2(3,0)
     $ + 24* HZr1(1)*GYZ3(3,1,0)
     $ - 22* HZr1(1)*GYZ2(0,3)
     $ + 53* HZr1(1)*GYZ1(0)

     $ + 4* HZr1(1)*GYZ2(0,0)
     $ - 12* HZr1(1)*GYZ3(0,1,0)
     $ + 24* HZr1(1)*GYZ3(1,3,0)
     $ + 24* HZr1(1)*GYZ3(1,0,3)

     $ - 12* HZr1(1)*GYZ2(1,0)
     $ - 24* HZr1(1)*GYZ3(1,0,0)
     $ + 11* HZr2(1,0)
     $ - 12* HZr2(1,0)*GYZ1(2)

     $ - 12* HZr2(1,0)*GYZ2(2,0)
     $ + 24* HZr2(1,0)*GYZ2(3,2)
     $ - 22* HZr2(1,0)*GYZ1(3)

     $ - 24* HZr2(1,0)*GYZ2(3,0)
     $ - 12* HZr2(1,0)*GYZ2(0,2)
     $ - 24* HZr2(1,0)*GYZ2(0,3)
     $ + 38* HZr2(1,0)*GYZ1(0)

     $ + 24* HZr2(1,0)*GYZ2(0,0)
     $ + 12* HZr2(1,0)*GYZ2(1,0)
     $ + 12* HZr3(1,0,0)
     $ - 24* HZr3(1,0,0)*GYZ1(2)

     $ - 24* HZr4(1,0,0,1)
     $ + 12* HZr3(1,0,1)
     $ - 24* HZr3(1,0,1)*GYZ1(3)
     $ + 12* HZr3(1,0,1)*GYZ1(0)
     $ - 12* HZr4(1,0,1,0)

     $ + 24* HZr3(1,1,0)
     $ - 48* HZr3(1,1,0)*GYZ1(3)
     $ + 24* HZr3(1,1,0)*GYZ1(0)
     $ - 42* GYZ1(2)
     $ - 53* GYZ2(2,0)

     $ - 4* GYZ3(2,0,0)
     $ + 12* GYZ4(2,0,1,0)
     $ + 12* GYZ3(2,1,0)
     $ + 24* GYZ4(2,1,0,0)

     $ + 84* GYZ2(3,2)
     $ + 22* GYZ3(3,2,0)
     $ - 24* GYZ4(3,2,1,0)
     $ + 22* GYZ3(3,0,2)

     $ - 24* GYZ4(3,1,2,0)
     $ - 24* GYZ4(3,1,0,2)
     $ - 53* GYZ2(0,2)
     $ - 4* GYZ3(0,2,0)

     $ + 12* GYZ4(0,2,1,0)
     $ + 22* GYZ3(0,3,2)
     $ + 21* GYZ1(0)
     $ - 4* GYZ3(0,0,2)
     $ + 2* GYZ2(0,0)

     $ + 12* GYZ4(0,1,2,0)
     $ + 12* GYZ4(0,1,0,2)
     $ - 4* GYZ3(0,1,0)
     $ + 12* GYZ3(1,2,0)
     $ + 24* GYZ4(1,2,0,0)

     $ - 24* GYZ4(1,3,2,0)
     $ - 24* GYZ4(1,3,0,2)
     $ + 12* GYZ3(1,0,2)
     $ + 24* GYZ4(1,0,2,0)

     $ - 24* GYZ4(1,0,3,2)
     $ + 42* GYZ2(1,0)
     $ + 24* GYZ4(1,0,0,2)
     $ - 8* GYZ3(1,0,0)
     $ - 24* GYZ4(1,0,1,0)

     $ - 48* GYZ4(1,1,0,0)
     $)


      else if(i.eq.3) then
      fin11=
     $(z*(1-z)**2)/(y**3)*(
     $ - 2* HZr1(0)*GYZ3(2,3,2)
     $ - 4* HZr1(0)*GYZ3(3,2,2)
     $ + 4* HZr2(0,0)*GYZ2(2,2)
     $ + 2* HZr2(0,1)*GYZ2(2,3)
     $ + 2* HZr2(0,1)*GYZ2(3,2)
     $ - 2* HZr1(1)*GYZ3(3,2,3)
     $ - 4* HZr1(1)*GYZ3(3,3,2)
     $ + 2* HZr2(1,0)*GYZ2(2,3)
     $ + 2* HZr2(1,0)*GYZ2(3,2)
     $ + 4* HZr2(1,1)*GYZ2(3,3)
     $ + 2* GYZ4(3,2,3,2)
     $ + 4* GYZ4(3,3,2,2)
     $)
     $+ (z)/(y**2)*(
     $ - 8* HZr1(0)*GYZ2(2,2)
     $ + 12* HZr1(0)*GYZ3(2,3,2)
     $ + 24* HZr1(0)*GYZ3(3,2,2)
     $ - 2* HZr1(0)*GYZ2(3,2)
     $ - 24* HZr2(0,0)*GYZ2(2,2)
     $ + 4* HZr2(0,0)*GYZ1(2)
     $ - 12* HZr2(0,1)*GYZ2(2,3)
     $ + 4* HZr2(0,1)*GYZ1(2)
     $ - 12* HZr2(0,1)*GYZ2(3,2)
     $ + 2* HZr2(0,1)*GYZ1(3)
     $ - 4* HZr1(1)*GYZ2(2,3)
     $ + 12* HZr1(1)*GYZ3(3,2,3)
     $ - 8* HZr1(1)*GYZ2(3,2)
     $ + 24* HZr1(1)*GYZ3(3,3,2)
     $ - 12* HZr2(1,0)*GYZ2(2,3)
     $ + 4* HZr2(1,0)*GYZ1(2)
     $ - 12* HZr2(1,0)*GYZ2(3,2)
     $ + 2* HZr2(1,0)*GYZ1(3)
     $ - 24* HZr2(1,1)*GYZ2(3,3)
     $ + 8* HZr2(1,1)*GYZ1(3)
     $ + 4* GYZ3(2,3,2)
     $ + 8* GYZ3(3,2,2)
     $ - 12* GYZ4(3,2,3,2)
     $ - 24* GYZ4(3,3,2,2)
     $)
     $+ (z**2)/(y**2)*(
     $   4* HZr1(0)*GYZ2(2,2)
     $ - 8* HZr1(0)*GYZ3(2,3,2)
     $ - 16* HZr1(0)*GYZ3(3,2,2)
     $ + 2* HZr1(0)*GYZ2(3,2)
     $ + 16* HZr2(0,0)*GYZ2(2,2)
     $ - 4* HZr2(0,0)*GYZ1(2)
     $ + 8* HZr2(0,1)*GYZ2(2,3)
     $ - 2* HZr2(0,1)*GYZ1(2)
     $ + 8* HZr2(0,1)*GYZ2(3,2)
     $ - 2* HZr2(0,1)*GYZ1(3)
     $ + 2* HZr1(1)*GYZ2(2,3)
     $ - 8* HZr1(1)*GYZ3(3,2,3)
     $ + 4* HZr1(1)*GYZ2(3,2)
     $ - 16* HZr1(1)*GYZ3(3,3,2)
     $ + 8* HZr2(1,0)*GYZ2(2,3)
     $ - 2* HZr2(1,0)*GYZ1(2)
     $ + 8* HZr2(1,0)*GYZ2(3,2)
     $ - 2* HZr2(1,0)*GYZ1(3)
     $ + 16* HZr2(1,1)*GYZ2(3,3)
     $ - 4* HZr2(1,1)*GYZ1(3)
     $ - 2* GYZ3(2,3,2)
     $ - 4* GYZ3(3,2,2)
     $ + 8* GYZ4(3,2,3,2)
     $ + 16* GYZ4(3,3,2,2)
     $)
     $+ (1.0d0)/(y**2)*(
     $ + 4* HZr1(0)*GYZ2(2,2)
     $ - 4* HZr1(0)*GYZ3(2,3,2)
     $ - 8* HZr1(0)*GYZ3(3,2,2)
     $ + 8* HZr2(0,0)*GYZ2(2,2)
     $ + 4* HZr2(0,1)*GYZ2(2,3)
     $ - 2* HZr2(0,1)*GYZ1(2)
     $ + 4* HZr2(0,1)*GYZ2(3,2)
     $ + 2* HZr1(1)*GYZ2(2,3)
     $ - 4* HZr1(1)*GYZ3(3,2,3)
     $ + 4* HZr1(1)*GYZ2(3,2)
     $ - 8* HZr1(1)*GYZ3(3,3,2)
     $ + 4* HZr2(1,0)*GYZ2(2,3)
     $ - 2* HZr2(1,0)*GYZ1(2)
     $ + 4* HZr2(1,0)*GYZ2(3,2)
     $ + 8* HZr2(1,1)*GYZ2(3,3)
     $ - 4* HZr2(1,1)*GYZ1(3)
     $ - 2* GYZ3(2,3,2)
     $ - 4* GYZ3(3,2,2)
     $ + 4* GYZ4(3,2,3,2)
     $ + 8* GYZ4(3,3,2,2)
     $)
     $+ (z)/(y)*(
     $   (1.0d0)/(4.0d0)
     $ + 14* HZr1(0)*GYZ2(2,2)
     $ + 4* HZr1(0)*GYZ3(2,2,0)
     $ - 14* HZr1(0)*GYZ3(2,3,2)
     $ + 5* HZr1(0)*GYZ1(2)
     $ + 4* HZr1(0)*GYZ3(2,0,2)
     $ - 2* HZr1(0)*GYZ3(2,1,0)
     $ - 28* HZr1(0)*GYZ3(3,2,2)
     $ + 7* HZr1(0)*GYZ2(3,2)
     $ + 4* HZr1(0)*GYZ3(0,2,2)
     $ - 2* HZr1(0)*GYZ3(1,2,0)
     $ - 2* HZr1(0)*GYZ3(1,0,2)
     $ + 2* HZr2(0,0)
     $ + 24* HZr2(0,0)*GYZ2(2,2)
     $ - 14* HZr2(0,0)*GYZ1(2)
     $ + 4* HZr3(0,0,1)*GYZ1(2)
     $ + 2* HZr2(0,1)
     $ + 14* HZr2(0,1)*GYZ2(2,3)
     $ - 7* HZr2(0,1)*GYZ1(2)
     $ - 2* HZr2(0,1)*GYZ2(2,0)
     $ + 12* HZr2(0,1)*GYZ2(3,2)
     $ - 7* HZr2(0,1)*GYZ1(3)
     $ - 2* HZr2(0,1)*GYZ2(0,2)
     $ + 2* HZr3(0,1,0)*GYZ1(2)
     $ + 4* HZr3(0,1,1)*GYZ1(3)
     $ + 7* HZr1(1)*GYZ2(2,3)
     $ + 2* HZr1(1)*GYZ3(2,3,0)
     $ - HZr1(1)*GYZ1(2)
     $ + 2* HZr1(1)*GYZ3(2,0,3)
     $ - 16* HZr1(1)*GYZ3(3,2,3)
     $ + 14* HZr1(1)*GYZ2(3,2)
     $ + 4* HZr1(1)*GYZ3(3,2,0)
     $ - 32* HZr1(1)*GYZ3(3,3,2)
     $ + 7* HZr1(1)*GYZ1(3)
     $ + 4* HZr1(1)*GYZ3(3,0,2)
     $ - 2* HZr1(1)*GYZ3(3,1,0)
     $ + 2* HZr1(1)*GYZ3(0,2,3)
     $ + 4* HZr1(1)*GYZ3(0,3,2)
     $ - 2* HZr1(1)*GYZ3(1,3,0)
     $ - 2* HZr1(1)*GYZ3(1,0,3)
     $ + 2* HZr2(1,0)
     $ + 14* HZr2(1,0)*GYZ2(2,3)
     $ - 7* HZr2(1,0)*GYZ1(2)
     $ - 2* HZr2(1,0)*GYZ2(2,0)
     $ + 14* HZr2(1,0)*GYZ2(3,2)
     $ - 7* HZr2(1,0)*GYZ1(3)
     $ - 2* HZr2(1,0)*GYZ2(0,2)
     $ + 2* HZr3(1,0,1)*GYZ1(3)
     $ + HZr2(1,1)
     $ + 32* HZr2(1,1)*GYZ2(3,3)
     $ - 14* HZr2(1,1)*GYZ1(3)
     $ - 4* HZr2(1,1)*GYZ2(3,0)
     $ - 4* HZr2(1,1)*GYZ2(0,3)
     $ + GYZ2(2,2)
     $ - 7* GYZ3(2,3,2)
     $ - 2* GYZ4(2,3,2,0)
     $ - 2* GYZ4(2,3,0,2)
     $ - 2* GYZ4(2,0,3,2)
     $ - 14* GYZ3(3,2,2)
     $ - 4* GYZ4(3,2,2,0)
     $ + 16* GYZ4(3,2,3,2)
     $ - 7* GYZ2(3,2)
     $ - 4* GYZ4(3,2,0,2)
     $ + 2* GYZ4(3,2,1,0)
     $ + 32* GYZ4(3,3,2,2)
     $ - 4* GYZ4(3,0,2,2)
     $ + 2* GYZ4(3,1,2,0)
     $ + 2* GYZ4(3,1,0,2)
     $ - 2* GYZ4(0,2,3,2)
     $ - 4* GYZ4(0,3,2,2)
     $ + 2* GYZ4(1,3,2,0)
     $ + 2* GYZ4(1,3,0,2)
     $ + 2* GYZ4(1,0,3,2)
     $)
     $+ (1.0d0)/(y)*(
     $ + 4
     $ - 20* HZr1(0)*GYZ2(2,2)
     $ - 8* HZr1(0)*GYZ3(2,2,0)
     $ + 16* HZr1(0)*GYZ3(2,3,2)
     $ - 13* HZr1(0)*GYZ1(2)
     $ - 8* HZr1(0)*GYZ3(2,0,2)
     $ + 4* HZr1(0)*GYZ3(2,1,0)
     $ + 32* HZr1(0)*GYZ3(3,2,2)
     $ - 4* HZr1(0)*GYZ2(3,2)
     $ - 8* HZr1(0)*GYZ3(0,2,2)
     $ + 4* HZr1(0)*GYZ3(1,2,0)
     $ + 4* HZr1(0)*GYZ3(1,0,2)
     $ - 24* HZr2(0,0)*GYZ2(2,2)
     $ + 8* HZr2(0,0)*GYZ1(2)
     $ - 8* HZr3(0,0,1)*GYZ1(2)
     $ - HZr2(0,1)
     $ - 16* HZr2(0,1)*GYZ2(2,3)
     $ + 10* HZr2(0,1)*GYZ1(2)
     $ + 4* HZr2(0,1)*GYZ2(2,0)
     $ - 12* HZr2(0,1)*GYZ2(3,2)
     $ + 4* HZr2(0,1)*GYZ1(3)
     $ + 4* HZr2(0,1)*GYZ2(0,2)
     $ - 4* HZr3(0,1,0)*GYZ1(2)
     $ - 8* HZr3(0,1,1)*GYZ1(3)
     $ - (3.0d0)/(2.0d0)*HZr1(1)
     $ - 10* HZr1(1)*GYZ2(2,3)
     $ - 4* HZr1(1)*GYZ3(2,3,0)
     $ + 2* HZr1(1)*GYZ1(2)
     $ - 4* HZr1(1)*GYZ3(2,0,3)
     $ + 20* HZr1(1)*GYZ3(3,2,3)
     $ - 20* HZr1(1)*GYZ2(3,2)
     $ - 8* HZr1(1)*GYZ3(3,2,0)
     $ + 40* HZr1(1)*GYZ3(3,3,2)
     $ - 14* HZr1(1)*GYZ1(3)
     $ - 8* HZr1(1)*GYZ3(3,0,2)
     $ + 4* HZr1(1)*GYZ3(3,1,0)
     $ - 4* HZr1(1)*GYZ3(0,2,3)
     $ - 8* HZr1(1)*GYZ3(0,3,2)
     $ - HZr1(1)*GYZ1(0)
     $ + 4* HZr1(1)*GYZ3(1,3,0)
     $ + 4* HZr1(1)*GYZ3(1,0,3)
     $ - 2* HZr2(1,0)
     $ - 16* HZr2(1,0)*GYZ2(2,3)
     $ + 10* HZr2(1,0)*GYZ1(2)
     $ + 4* HZr2(1,0)*GYZ2(2,0)
     $ - 16* HZr2(1,0)*GYZ2(3,2)
     $ + 4* HZr2(1,0)*GYZ1(3)
     $ + 4* HZr2(1,0)*GYZ2(0,2)
     $ - 4* HZr3(1,0,1)*GYZ1(3)
     $ - 2* HZr2(1,1)
     $ - 40* HZr2(1,1)*GYZ2(3,3)
     $ + 20* HZr2(1,1)*GYZ1(3)
     $ + 8* HZr2(1,1)*GYZ2(3,0)
     $ + 8* HZr2(1,1)*GYZ2(0,3)
     $ - 2* GYZ2(2,2)
     $ + 10* GYZ3(2,3,2)
     $ + 4* GYZ4(2,3,2,0)
     $ + 4* GYZ4(2,3,0,2)
     $ + (3.0d0)/(2.0d0)*GYZ1(2)
     $ + 4* GYZ4(2,0,3,2)
     $ + GYZ2(2,0)
     $ + 20* GYZ3(3,2,2)
     $ + 8* GYZ4(3,2,2,0)
     $ - 20* GYZ4(3,2,3,2)
     $ + 14* GYZ2(3,2)
     $ + 8* GYZ4(3,2,0,2)
     $ - 4* GYZ4(3,2,1,0)
     $ - 40* GYZ4(3,3,2,2)
     $ + 8* GYZ4(3,0,2,2)
     $ - 4* GYZ4(3,1,2,0)
     $ - 4* GYZ4(3,1,0,2)
     $ + 4* GYZ4(0,2,3,2)
     $ + GYZ2(0,2)
     $ + 8* GYZ4(0,3,2,2)
     $ - 4* GYZ4(1,3,2,0)
     $ - 4* GYZ4(1,3,0,2)
     $ - 4* GYZ4(1,0,3,2)
     $ - GYZ2(1,0)
     $)
     $+ (z)/(2*(1-y)**2)*(
     $   2*HZr1(0)*GYZ2(2,0)
     $ + 2*HZr1(0)*GYZ2(0,2)
     $ + 2*HZr2(0,1)*GYZ1(0)
     $ + 4* HZr1(1)*GYZ2(3,0)
     $ + 4* HZr1(1)*GYZ2(0,3)
     $ - 3* HZr1(1)*GYZ1(0)
     $ - 4* HZr1(1)*GYZ2(0,0)
     $ + 3* GYZ2(2,0)
     $ + 4* GYZ3(2,0,0)
     $ - 4* GYZ3(3,2,0)
     $ - 4* GYZ3(3,0,2)
     $ + 3* GYZ2(0,2)
     $ + 4* GYZ3(0,2,0)
     $ - 4* GYZ3(0,3,2)
     $ + 7* GYZ1(0)
     $ + 4* GYZ3(0,0,2)
     $ + 1* GYZ2(0,0)
     $ - 2* GYZ3(0,1,0)
     $ - 4* GYZ3(1,0,0)
     $)
     $+ (z)/(2*(1-y))*(
     $   7
     $ + 2* HZr1(0)*GYZ1(2)
     $ + 6* HZr1(0)*GYZ2(2,0)
     $ + 6* HZr1(0)*GYZ2(0,2)
     $ + 2* HZr2(0,1)
     $ + 6* HZr2(0,1)*GYZ1(0)
     $ - 3* HZr1(1)
     $ + 4* HZr1(1)*GYZ1(3)
     $ + 12* HZr1(1)*GYZ2(3,0)
     $ + 12* HZr1(1)*GYZ2(0,3)
     $ - 11* HZr1(1)*GYZ1(0)
     $ - 12* HZr1(1)*GYZ2(0,0)
     $ + 3* GYZ1(2)
     $ + 11* GYZ2(2,0)
     $ + 12* GYZ3(2,0,0)
     $ - 4* GYZ2(3,2)
     $ - 12* GYZ3(3,2,0)
     $ - 12* GYZ3(3,0,2)
     $ + 11* GYZ2(0,2)
     $ + 12* GYZ3(0,2,0)
     $ - 12* GYZ3(0,3,2)
     $ + 21* GYZ1(0)
     $ + 12* GYZ3(0,0,2)
     $ - 9* GYZ2(0,0)
     $ - 6* GYZ3(0,1,0)
     $ - 2* GYZ2(1,0)
     $ - 12* GYZ3(1,0,0)
     $)
     $+ (1.0d0)/((1-y)**2*(y+z)**2)*(
     $   HZr1(1)*GYZ1(0)
     $ - GYZ2(2,0)
     $ - GYZ2(0,2)
     $)
     $+ (1.0d0)/((1-y)**2*(y+z))*(
     $ - 2* HZr1(1)*GYZ1(0)
     $ + 2* GYZ2(2,0)
     $ + 2* GYZ2(0,2)
     $ - GYZ1(0)
     $)
     $+ (1.0d0)/((1-y)**2)*(
     $ + HZr1(1)*GYZ1(0)
     $ - GYZ2(2,0)
     $ - GYZ2(0,2)
     $ + GYZ1(0)
     $)
     $+ (1.0d0)/((1-y)*(y+z)**2)*(
     $ + HZr1(1)
     $ - GYZ1(2)
     $)
     $+ (1.0d0)/((1-y)*(y+z))*(
     $ - 1
     $ - 2* HZr1(1)
     $ + 2* GYZ1(2)
     $)
     $+ (1.0d0)/(1-y)*(
     $   1
     $ - 2* HZr1(0)*GYZ2(2,0)
     $ - 2* HZr1(0)*GYZ2(0,2)
     $ + HZr1(1)
     $ - 2* HZr1(1)*GYZ2(3,0)
     $ - 2* HZr1(1)*GYZ2(0,3)
     $ + 3* HZr1(1)*GYZ1(0)
     $ - GYZ1(2)
     $ - 3* GYZ2(2,0)
     $ + 2* GYZ3(3,2,0)
     $ + 2* GYZ3(3,0,2)
     $ - 3* GYZ2(0,2)
     $ + 2* GYZ3(0,3,2)
     $ - 8* GYZ1(0)
     $ - 2* GYZ2(0,0)
     $)
     $+ (1.0d0)/(2*(y+z)**2)*(
     $ - 4* HZr1(0)*GYZ2(2,2)
     $ +   HZr1(0)*GYZ1(2)
     $ -   HZr2(0,1)
     $ + 4* HZr3(0,1,1)
     $ + 14* HZr1(1)
     $ - 4* HZr1(1)*GYZ2(2,3)
     $ + 8* HZr1(1)*GYZ1(2)
     $ + 4* HZr1(1)*GYZ2(2,0)
     $ - 8* HZr1(1)*GYZ2(3,2)
     $ + 4* HZr1(1)*GYZ2(0,2)
     $ -  HZr1(1)*GYZ1(0)
     $ - 2* HZr1(1)*GYZ2(1,0)
     $ -  HZr2(1,0)
     $ + 2* HZr2(1,0)*GYZ1(2)
     $ + 2* HZr3(1,0,1)
     $ - 8* HZr2(1,1)
     $ + 8* HZr2(1,1)*GYZ1(3)
     $ - 4* HZr2(1,1)*GYZ1(0)
     $ - 8* GYZ2(2,2)
     $ - 4* GYZ3(2,2,0)
     $ + 4* GYZ3(2,3,2)
     $ - 14* GYZ1(2)
     $ - 4* GYZ3(2,0,2)
     $ +  GYZ2(2,0)
     $ + 2* GYZ3(2,1,0)
     $ + 8* GYZ3(3,2,2)
     $ - 4* GYZ3(0,2,2)
     $ +  GYZ2(0,2)
     $ + 2* GYZ3(1,2,0)
     $ + 2* GYZ3(1,0,2)
     $)
     $+ (1.0d0)/(2*(y+z))*(
     $ - 14
     $ +  HZr1(0)
     $ - 2* HZr1(0)*GYZ1(2)
     $ - 2* HZr2(0,1)
     $ + 6* HZr1(1)
     $ - 4* HZr1(1)*GYZ1(3)
     $ + 2* HZr1(1)*GYZ1(0)
     $ - 6* GYZ1(2)
     $ - 2* GYZ2(2,0)
     $ + 4* GYZ2(3,2)
     $ - 2* GYZ2(0,2)
     $ +  GYZ1(0)
     $ + 2* GYZ2(1,0)
     $)
     $+ T*(
     $   8
     $ + 3* HZr1(0)*GYZ2(2,2)
     $ + 2* HZr1(0)*GYZ3(2,2,0)
     $ - 2* HZr1(0)*GYZ3(2,3,2)
     $ + 4* HZr1(0)*GYZ1(2)
     $ + 2* HZr1(0)*GYZ3(2,0,2)
     $ - HZr1(0)*GYZ3(2,1,0)
     $ - 4* HZr1(0)*GYZ3(3,2,2)
     $ + 2* HZr1(0)*GYZ3(0,2,2)
     $ - HZr1(0)*GYZ3(1,2,0)
     $ - HZr1(0)*GYZ3(1,0,2)
     $ + 2* HZr2(0,0)*GYZ2(2,2)
     $ + 2* HZr3(0,0,1)*GYZ1(2)
     $ + 2* HZr4(0,0,1,1)
     $ + 4* HZr2(0,1)
     $ + 2* HZr2(0,1)*GYZ2(2,3)
     $ - HZr2(0,1)*GYZ2(1,0)
     $ + HZr3(0,1,0)*GYZ1(2)
     $ + HZr4(0,1,0,1)
     $ - 3* HZr3(0,1,1)
     $ + 4* HZr3(0,1,1)*GYZ1(3)
     $ - 2* HZr3(0,1,1)*GYZ1(0)
     $ - 6* HZr1(1)
     $ + 3* HZr1(1)*GYZ2(2,3)
     $ + 2* HZr1(1)*GYZ3(2,3,0)
     $ + 2* HZr1(1)*GYZ3(2,0,3)
     $ - 3* HZr1(1)*GYZ2(2,0)
     $ - 2* HZr1(1)*GYZ3(2,0,0)
     $ - 4* HZr1(1)*GYZ3(3,2,3)
     $ + 6* HZr1(1)*GYZ2(3,2)
     $ + 4* HZr1(1)*GYZ3(3,2,0)
     $ - 8* HZr1(1)*GYZ3(3,3,2)
     $ + 8* HZr1(1)*GYZ1(3)
     $ + 4* HZr1(1)*GYZ3(3,0,2)
     $ - 2* HZr1(1)*GYZ3(3,1,0)
     $ + 2* HZr1(1)*GYZ3(0,2,3)
     $ - 3* HZr1(1)*GYZ2(0,2)
     $ - 2* HZr1(1)*GYZ3(0,2,0)
     $ + 4* HZr1(1)*GYZ3(0,3,2)
     $ - 4* HZr1(1)*GYZ1(0)
     $ - 2* HZr1(1)*GYZ3(0,0,2)
     $ + HZr1(1)*GYZ3(0,1,0)
     $ - 2* HZr1(1)*GYZ3(1,3,0)
     $ - 2* HZr1(1)*GYZ3(1,0,3)
     $ + 2* HZr1(1)*GYZ3(1,0,0)
     $ + 2* HZr2(1,0)*GYZ2(2,3)
     $ - HZr2(1,0)*GYZ2(2,0)
     $ + 2* HZr2(1,0)*GYZ2(3,2)
     $ - HZr2(1,0)*GYZ2(0,2)
     $ + 2* HZr3(1,0,1)*GYZ1(3)
     $ - HZr3(1,0,1)*GYZ1(0)
     $ + 8* HZr2(1,1)*GYZ2(3,3)
     $ - 6* HZr2(1,1)*GYZ1(3)
     $ - 4* HZr2(1,1)*GYZ2(3,0)
     $ - 4* HZr2(1,1)*GYZ2(0,3)
     $ + 3* HZr2(1,1)*GYZ1(0)
     $ + 2* HZr2(1,1)*GYZ2(0,0)
     $ + 3* GYZ3(2,2,0)
     $ + 2* GYZ4(2,2,0,0)
     $ - 3* GYZ3(2,3,2)
     $ - 2* GYZ4(2,3,2,0)
     $ - 2* GYZ4(2,3,0,2)
     $ + 6* GYZ1(2)
     $ + 3* GYZ3(2,0,2)
     $ + 2* GYZ4(2,0,2,0)
     $ - 2* GYZ4(2,0,3,2)
     $ + 4* GYZ2(2,0)
     $ + 2* GYZ4(2,0,0,2)
     $ - GYZ4(2,0,1,0)
     $ - 2* GYZ4(2,1,0,0)
     $ - 6* GYZ3(3,2,2)
     $ - 4* GYZ4(3,2,2,0)
     $ + 4* GYZ4(3,2,3,2)
     $ - 8* GYZ2(3,2)
     $ - 4* GYZ4(3,2,0,2)
     $ + 2* GYZ4(3,2,1,0)
     $ + 8* GYZ4(3,3,2,2)
     $ - 4* GYZ4(3,0,2,2)
     $ + 2* GYZ4(3,1,2,0)
     $ + 2* GYZ4(3,1,0,2)
     $ + 3* GYZ3(0,2,2)
     $ + 2* GYZ4(0,2,2,0)
     $ - 2* GYZ4(0,2,3,2)
     $ + 4* GYZ2(0,2)
     $ + 2* GYZ4(0,2,0,2)
     $ - GYZ4(0,2,1,0)
     $ - 4* GYZ4(0,3,2,2)
     $ + 2* GYZ4(0,0,2,2)
     $ - GYZ4(0,1,2,0)
     $ - GYZ4(0,1,0,2)
     $ - 2* GYZ4(1,2,0,0)
     $ + 2* GYZ4(1,3,2,0)
     $ + 2* GYZ4(1,3,0,2)
     $ - 2* GYZ4(1,0,2,0)
     $ + 2* GYZ4(1,0,3,2)
     $ - 4* GYZ2(1,0)
     $ - 2* GYZ4(1,0,0,2)
     $ + GYZ4(1,0,1,0)
     $ + 2* GYZ4(1,1,0,0)
     $ + (11.0d0)/(4.0d0)*HZr2(1,1)
     $ - (11.0d0)/(4.0d0)*HZr1(1)*GYZ1(2)
     $ + (11.0d0)/(4.0d0)*GYZ2(2,2)
     $ + (3.0d0)/(2.0d0)*HZr1(1)*GYZ2(1,0)
     $ - (3.0d0)/(2.0d0)*HZr2(1,0)*GYZ1(2)
     $ - (3.0d0)/(2.0d0)*HZr3(1,0,1)
     $ - (3.0d0)/(2.0d0)*GYZ3(2,1,0)
     $ - (3.0d0)/(2.0d0)*GYZ3(1,2,0)
     $ - (3.0d0)/(2.0d0)*GYZ3(1,0,2)
     $)
      fin11 = fin11
     $ - (1.0d0)/(2.0d0)
     $ - (7.0d0)/(2.0d0)*HZr1(0)
     $ + 10* HZr1(0)*GYZ2(2,2)
     $ + 8* HZr1(0)*GYZ3(2,2,0)
     $ - 8* HZr1(0)*GYZ3(2,3,2)
     $ + 4* HZr1(0)*GYZ1(2)
     $ + 8* HZr1(0)*GYZ3(2,0,2)
     $ - 2* HZr1(0)*GYZ2(2,0)
     $ - 4* HZr1(0)*GYZ3(2,1,0)
     $ - 16* HZr1(0)*GYZ3(3,2,2)
     $ + 5* HZr1(0)*GYZ2(3,2)
     $ + 8* HZr1(0)*GYZ3(0,2,2)
     $ - 2* HZr1(0)*GYZ2(0,2)
     $ - 4* HZr1(0)*GYZ3(1,2,0)
     $ - 4* HZr1(0)*GYZ3(1,0,2)
     $ + HZr1(0)*GYZ2(1,0)
     $ + 3* HZr2(0,0)
     $ + 8* HZr2(0,0)*GYZ2(2,2)
     $ - 8* HZr2(0,0)*GYZ1(2)
     $ - 2* HZr3(0,0,1)
     $ + 8* HZr3(0,0,1)*GYZ1(2)
     $ + 8* HZr4(0,0,1,1)
     $ + 10* HZr2(0,1)
     $ + 8* HZr2(0,1)*GYZ2(2,3)
     $ - 5* HZr2(0,1)*GYZ1(3)
     $ - 3* HZr2(0,1)*GYZ1(0)
     $ - 4* HZr2(0,1)*GYZ2(1,0)
     $ - HZr3(0,1,0)
     $ + 4* HZr3(0,1,0)*GYZ1(2)
     $ + 4* HZr4(0,1,0,1)
     $ - 10* HZr3(0,1,1)
     $ + 16* HZr3(0,1,1)*GYZ1(3)
     $ - 8* HZr3(0,1,1)*GYZ1(0)
     $ - 7* HZr1(1)
     $ + 10* HZr1(1)*GYZ2(2,3)
     $ + 8* HZr1(1)*GYZ3(2,3,0)
     $ - 6* HZr1(1)*GYZ1(2)
     $ + 8* HZr1(1)*GYZ3(2,0,3)
     $ - 10* HZr1(1)*GYZ2(2,0)
     $ - 8* HZr1(1)*GYZ3(2,0,0)
     $ - 16* HZr1(1)*GYZ3(3,2,3)
     $ + 20* HZr1(1)*GYZ2(3,2)
     $ + 16* HZr1(1)*GYZ3(3,2,0)
     $ - 32* HZr1(1)*GYZ3(3,3,2)
     $ + 14* HZr1(1)*GYZ1(3)
     $ + 16* HZr1(1)*GYZ3(3,0,2)
     $ - 5* HZr1(1)*GYZ2(3,0)
     $ - 8* HZr1(1)*GYZ3(3,1,0)
     $ + 8* HZr1(1)*GYZ3(0,2,3)
     $ - 10* HZr1(1)*GYZ2(0,2)
     $ - 8* HZr1(1)*GYZ3(0,2,0)
     $ + 16* HZr1(1)*GYZ3(0,3,2)
     $ - 5* HZr1(1)*GYZ2(0,3)
     $ - 4* HZr1(1)*GYZ1(0)
     $ - 8* HZr1(1)*GYZ3(0,0,2)
     $ + 8* HZr1(1)*GYZ2(0,0)
     $ + 4* HZr1(1)*GYZ3(0,1,0)
     $ - 8* HZr1(1)*GYZ3(1,3,0)
     $ - 8* HZr1(1)*GYZ3(1,0,3)
     $ + 5* HZr1(1)*GYZ2(1,0)
     $ + 8* HZr1(1)*GYZ3(1,0,0)
     $ + 3* HZr2(1,0)
     $ + 8* HZr2(1,0)*GYZ2(2,3)
     $ - 5* HZr2(1,0)*GYZ1(2)
     $ - 4* HZr2(1,0)*GYZ2(2,0)
     $ + 8* HZr2(1,0)*GYZ2(3,2)
     $ - 5* HZr2(1,0)*GYZ1(3)
     $ - 4* HZr2(1,0)*GYZ2(0,2)
     $ + HZr2(1,0)*GYZ1(0)
     $ - 5* HZr3(1,0,1)
     $ + 8* HZr3(1,0,1)*GYZ1(3)
     $ - 4* HZr3(1,0,1)*GYZ1(0)
     $ + 6* HZr2(1,1)
     $ + 32* HZr2(1,1)*GYZ2(3,3)
     $ - 20* HZr2(1,1)*GYZ1(3)
     $ - 16* HZr2(1,1)*GYZ2(3,0)
     $ - 16* HZr2(1,1)*GYZ2(0,3)
     $ + 10* HZr2(1,1)*GYZ1(0)
     $ + 8* HZr2(1,1)*GYZ2(0,0)
     $ + 6* GYZ2(2,2)
     $ + 10* GYZ3(2,2,0)
     $ + 8* GYZ4(2,2,0,0)
     $ - 10* GYZ3(2,3,2)
     $ - 8* GYZ4(2,3,2,0)
     $ - 8* GYZ4(2,3,0,2)
     $ + 7* GYZ1(2)
     $ + 10* GYZ3(2,0,2)
     $ + 8* GYZ4(2,0,2,0)
     $ - 8* GYZ4(2,0,3,2)
     $ + 4* GYZ2(2,0)
     $ + 8* GYZ4(2,0,0,2)
     $ - 8* GYZ3(2,0,0)
     $ - 4* GYZ4(2,0,1,0)
     $ - 5* GYZ3(2,1,0)
     $ - 8* GYZ4(2,1,0,0)
     $ - 20* GYZ3(3,2,2)
     $ - 16* GYZ4(3,2,2,0)
     $ + 16* GYZ4(3,2,3,2)
     $ - 14* GYZ2(3,2)
     $ - 16* GYZ4(3,2,0,2)
     $ + 5* GYZ3(3,2,0)
     $ + 8* GYZ4(3,2,1,0)
     $ + 32* GYZ4(3,3,2,2)
     $ - 16* GYZ4(3,0,2,2)
     $ + 5* GYZ3(3,0,2)
     $ + 8* GYZ4(3,1,2,0)
     $ + 8* GYZ4(3,1,0,2)
     $ + 10* GYZ3(0,2,2)
     $ + 8* GYZ4(0,2,2,0)
     $ - 8* GYZ4(0,2,3,2)
     $ + 4* GYZ2(0,2)
     $ + 8* GYZ4(0,2,0,2)
     $ - 8* GYZ3(0,2,0)
     $ - 4* GYZ4(0,2,1,0)
     $ - 16* GYZ4(0,3,2,2)
     $ + 5* GYZ3(0,3,2)
     $ - (7.0d0)/(2.0d0)*GYZ1(0)
     $ + 8* GYZ4(0,0,2,2)
     $ - 8* GYZ3(0,0,2)
     $ + 3* GYZ2(0,0)
     $ - 4* GYZ4(0,1,2,0)
     $ - 4* GYZ4(0,1,0,2)
     $ + 4* GYZ3(0,1,0)
     $ - 5* GYZ3(1,2,0)
     $ - 8* GYZ4(1,2,0,0)
     $ + 8* GYZ4(1,3,2,0)
     $ + 8* GYZ4(1,3,0,2)
     $ - 5* GYZ3(1,0,2)
     $ - 8* GYZ4(1,0,2,0)
     $ + 8* GYZ4(1,0,3,2)
     $ - 7* GYZ2(1,0)
     $ - 8* GYZ4(1,0,0,2)
     $ + 8* GYZ3(1,0,0)
     $ + 4* GYZ4(1,0,1,0)
     $ + 8* GYZ4(1,1,0,0)


      else if(i.eq.4) then
      fin11=
     $(1.0d0)/(6*y)*(
     $      HZr1(0)
     $ +    GYZ1(0)
     $)
     $+ (z)/(6*(1-y)**2)*(
     $ -  HZr1(0)*GYZ1(0)
     $ - 2*GYZ2(0,0)
     $)
     $+ (z)/(6*(1-y))*(
     $ - HZr1(0)

     $ - 3*HZr1(0)*GYZ1(0)
     $ -  GYZ1(0)
     $ - 6*GYZ2(0,0)
     $)
     $+ (1.0d0)/(3*(1-y))*(
     $   2*HZr1(0)*GYZ1(0)
     $ + 4*GYZ2(0,0)
     $)

     $+ (T*pi**2)/(36.0d0)*(
     $ - 22
     $ - HZr1(0)
     $ - GYZ1(0)
     $)
     $+ (T)/(18.0d0)*(
     $ - 12*HZr1(0)
     $ - 10*HZr1(0)*GYZ1(0)
     $ - 6*HZr1(0)*GYZ2(0,0)

     $ + 3*HZr1(0)*GYZ2(1,0)
     $ - 10*HZr2(0,0)
     $ - 6*HZr2(0,0)*GYZ1(0)
     $ - 3*HZr3(0,1,0)
     $ - 3*HZr2(1,0)*GYZ1(0)

     $ - 6*HZr3(1,0,0)
     $ - 12*GYZ1(0)
     $ - 10*GYZ2(0,0)
     $ + 3*GYZ3(0,1,0)
     $ + 6*GYZ3(1,0,0)
     $)


      else if(i.eq.5) then
      fin11=
     $(z)/(3*y)*(
     $     HZr1(0)*GYZ2(2,0)
     $ -   HZr1(0)*GYZ2(3,2)
     $ +   HZr1(0)*GYZ2(0,2)
     $ + 2* HZr2(0,0)*GYZ1(2)

     $ +   HZr2(0,1)*GYZ1(3)
     $ +   HZr1(1)*GYZ2(3,0)
     $ +   HZr1(1)*GYZ2(0,3)
     $ +   HZr2(1,0)*GYZ1(3)
     $ -   GYZ3(3,2,0)

     $ -   GYZ3(3,0,2)
     $ -   GYZ3(0,3,2)
     $)
     $+ (1.0d0)/(6*y)*(
     $     HZr1(0)
     $ - 4*HZr1(0)*GYZ2(2,0)

     $ + 4*HZr1(0)*GYZ2(3,2)
     $ - 4*HZr1(0)*GYZ2(0,2)
     $ - 8*HZr2(0,0)*GYZ1(2)
     $ - 4*HZr2(0,1)*GYZ1(3)

     $ - 4*HZr1(1)*GYZ2(3,0)
     $ - 4*HZr1(1)*GYZ2(0,3)
     $ - 4*HZr2(1,0)*GYZ1(3)
     $ + 4*GYZ3(3,2,0)

     $ + 4*GYZ3(3,0,2)
     $ + 4*GYZ3(0,3,2)
     $ +     GYZ1(0)
     $ )
     $+ (z)/(6*(1-y)**2)*(
     $     HZr1(0)*GYZ1(0)
     $ + 2*GYZ2(0,0)
     $ )

     $+ (z)/(6*(1-y))*(
     $     HZr1(0)
     $ + 3*HZr1(0)*GYZ1(0)
     $ +     GYZ1(0)
     $ + 6*GYZ2(0,0)
     $ )
     $+ (1.0d0)/(3*(1-y))*(
     $ -   HZr1(0)*GYZ1(0)

     $ - 2* GYZ2(0,0)
     $ )
     $+ (1.0d0)/(3*(y+z)**2)*(
     $ - HZr1(0)*GYZ1(2)
     $ + HZr2(0,1)
     $ + HZr1(1)*GYZ1(0)
     $ + HZr2(1,0)

     $ - GYZ2(2,0)
     $ - GYZ2(0,2)
     $ )
     $+ (1.0d0)/(3*(y+z))*(
     $ - HZr1(0)
     $ - GYZ1(0)
     $ )
     $+ (T)/(12.0d0)*(
     $   8* HZr1(0)

     $ + 3* HZr1(0)*GYZ1(2)
     $ + 4* HZr1(0)*GYZ2(2,0)
     $ - 4* HZr1(0)*GYZ2(3,2)
     $ + 4* HZr1(0)*GYZ2(0,2)

     $ - 2* HZr1(0)*GYZ2(1,0)
     $ + 4* HZr2(0,0)*GYZ1(2)
     $ + 4* HZr3(0,0,1)
     $ - 3* HZr2(0,1)
     $ + 4* HZr2(0,1)*GYZ1(3)

     $ + 2* HZr3(0,1,0)
     $ + 4* HZr1(1)*GYZ2(3,0)
     $ + 4* HZr1(1)*GYZ2(0,3)
     $ - 3* HZr1(1)*GYZ1(0)
     $ - 4* HZr1(1)*GYZ2(0,0)

     $ - 3* HZr2(1,0)
     $ + 4* HZr2(1,0)*GYZ1(3)
     $ - 2* HZr2(1,0)*GYZ1(0)
     $ + 3* GYZ2(2,0)
     $ + 4* GYZ3(2,0,0)

     $ - 4* GYZ3(3,2,0)
     $ - 4* GYZ3(3,0,2)
     $ + 3* GYZ2(0,2)
     $ + 4* GYZ3(0,2,0)
     $ - 4* GYZ3(0,3,2)

     $ + 8* GYZ1(0)
     $ + 4* GYZ3(0,0,2)
     $ - 2* GYZ3(0,1,0)
     $ - 4* GYZ3(1,0,0)
     $ )
     $+ (1.0d0)/(3.0d0)*(
     $     HZr1(0)*GYZ1(2)

     $ + 2* HZr1(0)*GYZ2(2,0)
     $ - 2* HZr1(0)*GYZ2(3,2)
     $ + 2* HZr1(0)*GYZ2(0,2)
     $ -   HZr1(0)*GYZ1(0)

     $ -   HZr1(0)*GYZ2(1,0)
     $ -   HZr2(0,0)
     $ + 2* HZr2(0,0)*GYZ1(2)
     $ + 2* HZr3(0,0,1)
     $ -   HZr2(0,1)

     $ + 2* HZr2(0,1)*GYZ1(3)
     $ +   HZr3(0,1,0)
     $ + 2* HZr1(1)*GYZ2(3,0)
     $ + 2* HZr1(1)*GYZ2(0,3)
     $ -   HZr1(1)*GYZ1(0)

     $ - 2* HZr1(1)*GYZ2(0,0)
     $ -   HZr2(1,0)
     $ + 2* HZr2(1,0)*GYZ1(3)
     $ -   HZr2(1,0)*GYZ1(0)
     $ +   GYZ2(2,0)

     $ + 2* GYZ3(2,0,0)
     $ - 2* GYZ3(3,2,0)
     $ - 2* GYZ3(3,0,2)
     $ +   GYZ2(0,2)
     $ + 2* GYZ3(0,2,0)

     $ - 2* GYZ3(0,3,2)
     $ + 2* GYZ3(0,0,2)
     $ -   GYZ2(0,0)
     $ -   GYZ3(0,1,0)
     $ - 2* GYZ3(1,0,0)
     $ )


      else if(i.eq.6) then
      fin11=T/36d0*(2*pi**2+HZr1(0)*GYZ1(0)+HZr2(0,0)+GYZ2(0,0))

!      else
      end if

      end function


************************************************************************
****                       Li_2(x) = H(0,1;x)                       ****
****        using T. Gehrmann and E. Remiddi's hplog.f code         ****
****                         hep-ph/0107173                         ****
************************************************************************
      real*8 function li2(x)
      implicit none
      real*8 x, Hr1(0:1),Hr2(0:1,0:1),Hr3(0:1,0:1,0:1),
     $ Hr4(0:1,0:1,0:1,0:1), Hi1(0:1),Hi2(0:1,0:1),Hi3(0:1,0:1,0:1),
     $ Hi4(0:1,0:1,0:1,0:1)
      complex*16 Hc1(0:1),Hc2(0:1,0:1),Hc3(0:1,0:1,0:1),
     $ Hc4(0:1,0:1,0:1,0:1)

      call hplog(x,2,Hc1,Hc2,Hc3,Hc4,Hr1,Hr2,Hr3,Hr4,
     $                 Hi1,Hi2,Hi3,Hi4,0,1)

      li2=Hr2(0,1)

      end function


************************************************************************
****                      real part of log(-y)                      ****
************************************************************************
      real*8 function lof(x,i)
      implicit none
      integer i
      real*8 x, pi, z3
      common /piez/ pi, z3

      if(i.eq.1) then
        lof=log(ABS(x))
      else if(i.eq.2) then
        if(x.ge.0.0d0) then
         lof=log(x)**2-pi**2
        else
         lof=log(ABS(x))**2
        end if
      else if(i.eq.3) then
        if(x.ge.0.0d0) then
         lof=log(x)**3-3*log(x)*pi**2
        else
         lof=log(ABS(x))**3
        end if
      else if(i.eq.4) then
        if(x.ge.0.0d0) then
         lof=log(x)**4-6*(log(x)*pi)**2+pi**4
        else
         lof=log(ABS(x))**4
        end if
      end if

      end function
