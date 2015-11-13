c-----------------------------------------------------------------------

      real(kind(1d0)) function dotprod(p1,p2)
      implicit none
c --- RETURN TWICE THE DOT PRODUCT OF P1 AND P2 (S_{12})
      real(kind(1d0)) p1,p2
      dimension p1(4),p2(4)
      dotprod = p1(4)*p2(4) - p1(1)*p2(1) - p1(2)*p2(2) - p1(3)*p2(3)
      dotprod = 2d0*dotprod
      end

c-----------------------------------------------------------------------

      FUNCTION sdot(P,I,J)
c --- RETURN TWICE THE DOT PRODUCT OF P(*,I) AND P(*,J) (S_{IJ})
      IMPLICIT NONE
      INTEGER I,J
      REAL(KIND(1D0)) sdot,P(4,7)
      sdot=P(4,I)*P(4,J)-P(3,I)*P(3,J)-P(2,I)*P(2,J)-P(1,I)*P(1,J)
      sdot=2d0*sdot
      END

c-----------------------------------------------------------------------

      subroutine getdotproducts(plab,n,slin)
      implicit none
      integer i,j,n
      real(kind(1d0)) plab,slin,xp,pi,pj,dotprod
      dimension plab(4,11),xp(4,11),slin(11,11),pi(4),pj(4)

c if part is commented out by AK
c      if (x.lt.1) then
c       write(6,*)'x<1 not implemented'
c      else
       do i=1,4
        do j=1,n
         xp(i,j) = plab(i,j)
        enddo
       enddo
c      endif

      do i=1,n
       do j=i+1,n
        call selectmom(xp,pi,i)
        call selectmom(xp,pj,j)
        slin(i,j) = dotprod(pi,pj)
        slin(j,i) = slin(i,j)
       enddo
      enddo

      end

c-----------------------------------------------------------------------

      subroutine selectmom(plab,p,i)
      implicit none
      integer i,j
      real(kind(1d0)) plab,p
      dimension plab(4,11),p(4)

      do j=1,4
       p(j) = plab(j,i)
      enddo

      end

c-----------------------------------------------------------------------

      subroutine spinorprod(p1,p2,a,b)
      implicit none
      real(kind(1d0)) p1,p2
      real(kind(1d0)) p1t,p2t,p1p,p1m,p2p,p2m
      dimension p1(4),p2(4)
      complex(kind(1d0)) a,b,expphase1p,expphase1m,expphase2p,expphase2m

      p1t = sqrt(p1(1)**2+p1(2)**2)
      p2t = sqrt(p2(1)**2+p2(2)**2)

      if (p1(4).gt.0d0) then
       p1p = p1(4)+p1(3)
       p1m = p1(4)-p1(3)
       if (p1t.gt.0d0) then
        expphase1p = cmplx(p1(1), p1(2))/p1t
        expphase1m = cmplx(p1(1),-p1(2))/p1t
       else
        expphase1p = cmplx(1d0,0d0)
        expphase1m = cmplx(1d0,0d0)
       endif

       if (p2(4).gt.0d0) then
        p2p = p2(4)+p2(3)
        p2m = p2(4)-p2(3)
        if (p2t.gt.0d0) then
         expphase2p = cmplx(p2(1), p2(2))/p2t
         expphase2m = cmplx(p2(1),-p2(2))/p2t
        else
         expphase2p = cmplx(1d0,0d0)
         expphase2m = cmplx(1d0,0d0)
        endif

        if (p1p.lt.0) p1p = -p1p
        if (p1m.lt.0) p1m = -p1m
        if (p2p.lt.0) p2p = -p2p
        if (p2m.lt.0) p2m = -p2m

        a =-(sqrt(p1m*p2p)*expphase1p - sqrt(p1p*p2m)*expphase2p)
        b =  sqrt(p1m*p2p)*expphase1m - sqrt(p1p*p2m)*expphase2m
       else
        p2p =-(p2(4)+p2(3))
        p2m =-(p2(4)-p2(3))
        if (p2t.gt.0d0) then
         expphase2p =-cmplx(p2(1), p2(2))/p2t
         expphase2m =-cmplx(p2(1),-p2(2))/p2t
        else
         expphase2p = cmplx(1d0,0d0)
         expphase2m = cmplx(1d0,0d0)
        endif

        if (p1p.lt.0) p1p = -p1p
        if (p1m.lt.0) p1m = -p1m
        if (p2p.lt.0) p2p = -p2p
        if (p2m.lt.0) p2m = -p2m

        a = sqrt(p1m*p2p)*expphase1p - sqrt(p1p*p2m)*expphase2p
        a = cmplx(0d0,1d0)*a
        b = sqrt(p1m*p2p)*expphase1m - sqrt(p1p*p2m)*expphase2m
        b =-cmplx(0d0,1d0)*b
       endif
      else
       p1p =-(p1(4)+p1(3))
       p1m =-(p1(4)-p1(3))
       if (p1t.gt.0d0) then
        expphase1p =-cmplx(p1(1), p1(2))/p1t
        expphase1m =-cmplx(p1(1),-p1(2))/p1t
       else
        expphase1p = cmplx(1d0,0d0)
        expphase1m = cmplx(1d0,0d0)
       endif

       if (p2(4).gt.0d0) then
        p2p = p2(4)+p2(3)
        p2m = p2(4)-p2(3)
        if (p2t.gt.0d0) then
         expphase2p = cmplx(p2(1), p2(2))/p2t
         expphase2m = cmplx(p2(1),-p2(2))/p2t
        else
         expphase2p = cmplx(1d0,0d0)
         expphase2m = cmplx(1d0,0d0)
        endif

        if (p1p.lt.0) p1p = -p1p
        if (p1m.lt.0) p1m = -p1m
        if (p2p.lt.0) p2p = -p2p
        if (p2m.lt.0) p2m = -p2m

        a = sqrt(p1m*p2p)*expphase1p - sqrt(p1p*p2m)*expphase2p
        a = cmplx(0d0,1d0)*a
        b = sqrt(p1m*p2p)*expphase1m - sqrt(p1p*p2m)*expphase2m
        b =-cmplx(0d0,1d0)*b
       else
        p2p =-(p2(4)+p2(3))
        p2m =-(p2(4)-p2(3))
        if (p2t.gt.0d0) then
         expphase2p =-cmplx(p2(1), p2(2))/p2t
         expphase2m =-cmplx(p2(1),-p2(2))/p2t
        else
         expphase2p = cmplx(1d0,0d0)
         expphase2m = cmplx(1d0,0d0)
        endif

        if (p1p.lt.0) p1p = -p1p
        if (p1m.lt.0) p1m = -p1m
        if (p2p.lt.0) p2p = -p2p
        if (p2m.lt.0) p2m = -p2m

        a = sqrt(p1m*p2p)*expphase1p - sqrt(p1p*p2m)*expphase2p
        b = sqrt(p1m*p2p)*expphase1m - sqrt(p1p*p2m)*expphase2m
        b =-b
       endif
      endif

      end

c-----------------------------------------------------------------------

      subroutine getspinorproducts(plab,n,ain,bin)
      implicit none
      integer i,j,n
      real(kind(1d0)) plab,xp,pi,pj
      dimension plab(4,11),xp(4,11),pi(4),pj(4)
      complex(kind(1d0)) a,b,ain(11,11),bin(11,11)

c if part is commented out by AK
c      if (x.lt.1) then
c       write(6,*)'x<1 not implemented'
c      else
       do i=1,4
        do j=1,n
         xp(i,j) = plab(i,j)
        enddo
       enddo
c      endif

      do i=1,n
       do j=i+1,n
        call selectmom(xp,pi,i)
        call selectmom(xp,pj,j)
        call spinorprod(pi,pj,a,b)
        ain(i,j) = a
        bin(i,j) = b
        ain(j,i) = -a
        bin(j,i) = -b
       enddo
      enddo

      end

c-----------------------------------------------------------------------
      subroutine LtoH(K,P,KPP,KPM,KMP,KMM)

C Subroutine calculates
C
C     Eps(P)^r_mu (K^mu K^nu)/K^2 Eps(P)^s_nu
C
C for massive momentum K and P is massless

      implicit none

C Variables

      integer i

      double precision K(4),P(4),k1(4),k2(4),ptmp1(4),ptmp2(4),c1,c2
     >      ,cf,sf,ct,st
      double complex KPP,KPM,KMP,KMM
     >          ,Ak1P,Bk1P,Ak2P,Bk2P,Ak1k2,Bk1k2

      double precision dotprod
      external dotprod

C
C Begin code
C

C Set up two massless momenta, k1 and k2, such that K = +/-(k1 + k2)

      cf = K(2)/sqrt(K(1)**2 + K(2)**2)
      sf = K(1)/sqrt(K(1)**2 + K(2)**2)

      ct = K(3)/sqrt(K(1)**2 + K(2)**2 + K(3)**2)
      st = sqrt(K(1)**2 + K(2)**2)/sqrt(K(1)**2 + K(2)**2+K(3)**2)

      call rotatez(cf,-sf,K,ptmp1)
      call rotatex(ct,-st,ptmp1,ptmp2)

c      write(6,*)'K       = ',K,dotprod(K,K)
c      write(6,*)'ptmp1   = ',ptmp1,dotprod(ptmp1,ptmp1)
c      write(6,*)'ptmp2   = ',ptmp2,dotprod(ptmp2,ptmp2)
c      write(6,*)'************************************'

      if ( (ptmp2(4)+ptmp2(3).gt.0d0).and.(ptmp2(4)-ptmp2(3).gt.0d0) )
     >     then
         k1(4) = (ptmp2(4)+ptmp2(3))/2d0
         k1(3) = (ptmp2(4)+ptmp2(3))/2d0
         k1(1) = 0d0
         k1(2) = 0d0

         k2(4) = (ptmp2(4)-ptmp2(3))/2d0
         k2(3) = -(ptmp2(4)-ptmp2(3))/2d0
         k2(1) = 0d0
         k2(2) = 0d0

         c1 = 1d0
         c2 = 1d0
c         write(6,*)'k1 + k2 = ',k1(1)+k2(1),k1(2)+k2(2)
c     >        ,k1(3)+k2(3),k1(4)+k2(4)
c         write(6,*)'k1 = ',k1,dotprod(k1,k1)
c         write(6,*)'k2 = ',k2,dotprod(k2,k2)        
c         write(6,*)'c1 , c2 = ',c1,c2
c         write(6,*)'^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^'
      endif

      if ( (ptmp2(4)+ptmp2(3).gt.0d0).and.(ptmp2(4)-ptmp2(3).lt.0d0) )
     >     then
         k1(4) = (ptmp2(4)+ptmp2(3))/2d0
         k1(3) = (ptmp2(4)+ptmp2(3))/2d0
         k1(1) = 0d0
         k1(2) = 0d0

         k2(4) = -(ptmp2(4)-ptmp2(3))/2d0
         k2(3) = +(ptmp2(4)-ptmp2(3))/2d0
         k2(1) = 0d0
         k2(2) = 0d0

         c1 = 1d0
         c2 = -1d0
c         write(6,*)'k1 + k2 = ',k1(1)-k2(1),k1(2)-k2(2)
c     >        ,k1(3)-k2(3),k1(4)-k2(4)
c         write(6,*)'k1 = ',k1,dotprod(k1,k1)
c         write(6,*)'k2 = ',k2,dotprod(k2,k2)        
c         write(6,*)'c1 , c2 = ',c1,c2
c         write(6,*)'^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^'
      endif

      if ( (ptmp2(4)+ptmp2(3).lt.0d0).and.(ptmp2(4)-ptmp2(3).gt.0d0) )
     >     then
         k1(4) = -(ptmp2(4)+ptmp2(3))/2d0
         k1(3) = -(ptmp2(4)+ptmp2(3))/2d0
         k1(1) = 0d0
         k1(2) = 0d0

         k2(4) = (ptmp2(4)-ptmp2(3))/2d0
         k2(3) = -(ptmp2(4)-ptmp2(3))/2d0
         k2(1) = 0d0
         k2(2) = 0d0

         c1 = -1d0
         c2 = 1d0
c         write(6,*)'k1 + k2 = ',-k1(1)+k2(1),-k1(2)+k2(2)
c     >        ,-k1(3)+k2(3),-k1(4)+k2(4)
c         write(6,*)'k1 = ',k1,dotprod(k1,k1)
c         write(6,*)'k2 = ',k2,dotprod(k2,k2)        
c         write(6,*)'c1 , c2 = ',c1,c2
c         write(6,*)'^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^'
      endif

      if ( (ptmp2(4)+ptmp2(3).lt.0d0).and.(ptmp2(4)-ptmp2(3).lt.0d0) )
     >     then
         k1(4) = -(ptmp2(4)+ptmp2(3))/2d0
         k1(3) = -(ptmp2(4)+ptmp2(3))/2d0
         k1(1) = 0d0
         k1(2) = 0d0

         k2(4) = -(ptmp2(4)-ptmp2(3))/2d0
         k2(3) = +(ptmp2(4)-ptmp2(3))/2d0
         k2(1) = 0d0
         k2(2) = 0d0

         c1 = -1d0
         c2 = -1d0
c         write(6,*)'k1 + k2 = ',-k1(1)-k2(1),-k1(2)-k2(2)
c     >        ,-k1(3)-k2(3),-k1(4)-k2(4)
c         write(6,*)'k1 = ',k1,dotprod(k1,k1)
c         write(6,*)'k2 = ',k2,dotprod(k2,k2)        
c         write(6,*)'c1 , c2 = ',c1,c2
c         write(6,*)'^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^'
      endif

      call rotatex(ct,st,k1,ptmp1)
      call rotatez(cf,sf,ptmp1,k1)

      call rotatex(ct,st,k2,ptmp1)
      call rotatez(cf,sf,ptmp1,k2)

C Some tests

c      write(6,*)'P       = ',P
c      write(6,*)'K       = ',K,dotprod(K,K)/2d0
c      write(6,*)'k1 + k2 = ',c1*k1(1)+c2*k2(1),c1*k1(2)+c2*k2(2)
c     >                      ,c1*k1(3)+c2*k2(3),c1*k1(4)+c2*k2(4)
c      write(6,*)'k1 = ',k1,dotprod(k1,k1)
c      write(6,*)'k2 = ',k2,dotprod(k2,k2)
c      write(6,*)

C Calculate the spinorproducts

      call spinorprod(k1,P,Ak1P,Bk1P)
      call spinorprod(k2,P,Ak2P,Bk2P)
      call spinorprod(k1,k2,Ak1k2,Bk1k2)

      KPP = cmplx(-1d0/2d0,0)
      KMM = KPP
      KPM = 1d0/2d0 * (Bk1P * Ak1k2 * Bk2P)/(Ak1P * Bk1k2 * Ak2P)
      KMP = conjg(KPM)

C
C End code
C

      end
*************************************************************************
************************************************************************
C     For the rotation
************************************************************************
      subroutine rotatex(cosp,sinp,p,prot)
      double precision cosp,sinp,p(4),prot(4)

C Performs a rotation of four-vector p around x with angle phi,
C (cosp=cos(phi)) with convention p=(x,y,z,e)

      prot(4) = p(4)
      prot(1) = p(1)
      prot(2) = cosp*p(2)+sinp*p(3)
      prot(3) = cosp*p(3)-sinp*p(2)

      return
      entry rotatey(cosp,sinp,p,prot)

C Performs a rotation of four-vector p around y with angle phi,

      prot(4) = p(4)
      prot(2) = p(2)
      prot(3) = cosp*p(3)+sinp*p(1)
      prot(1) = cosp*p(1)-sinp*p(3)

      return
      entry rotatez(cosp,sinp,p,prot)

C Performs a rotation of four-vector p around z with angle phi,

      prot(4) = p(4)
      prot(3) = p(3)
      prot(1) = cosp*p(1)+sinp*p(2)
      prot(2) = cosp*p(2)-sinp*p(1)

      end

************************************************************************

c-----------------------------------------------------------------------
      subroutine MakeMasslessF77(K,K1,K2)
      implicit none

C Variables

      integer i

      double precision K(4),k1(4),k2(4),ptmp1(4),ptmp2(4),c1,c2
     >      ,cf,sf,ct,st

C
C Begin code
C

C Set up two massless momenta, k1 and k2, such that K = +/-(k1 + k2)

      cf = K(2)/sqrt(K(1)**2 + K(2)**2)
      sf = K(1)/sqrt(K(1)**2 + K(2)**2)

      ct = K(3)/sqrt(K(1)**2 + K(2)**2 + K(3)**2)
      st = sqrt(K(1)**2 + K(2)**2)/sqrt(K(1)**2 + K(2)**2+K(3)**2)

      call rotatez(cf,-sf,K,ptmp1)
      call rotatex(ct,-st,ptmp1,ptmp2)

c      write(6,*)'K       = ',K,dotprod(K,K)
c      write(6,*)'ptmp1   = ',ptmp1,dotprod(ptmp1,ptmp1)
c      write(6,*)'ptmp2   = ',ptmp2,dotprod(ptmp2,ptmp2)
c      write(6,*)'************************************'

      if ( (ptmp2(4)+ptmp2(3).gt.0d0).and.(ptmp2(4)-ptmp2(3).gt.0d0) )
     >     then
         k1(4) = (ptmp2(4)+ptmp2(3))/2d0
         k1(3) = (ptmp2(4)+ptmp2(3))/2d0
         k1(1) = 0d0
         k1(2) = 0d0

         k2(4) = (ptmp2(4)-ptmp2(3))/2d0
         k2(3) = -(ptmp2(4)-ptmp2(3))/2d0
         k2(1) = 0d0
         k2(2) = 0d0

         c1 = 1d0
         c2 = 1d0
c         write(6,*)'k1 + k2 = ',k1(1)+k2(1),k1(2)+k2(2)
c     >        ,k1(3)+k2(3),k1(4)+k2(4)
c         write(6,*)'k1 = ',k1,dotprod(k1,k1)
c         write(6,*)'k2 = ',k2,dotprod(k2,k2)        
c         write(6,*)'c1 , c2 = ',c1,c2
c         write(6,*)'^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^'
      endif

      if ( (ptmp2(4)+ptmp2(3).gt.0d0).and.(ptmp2(4)-ptmp2(3).lt.0d0) )
     >     then
         k1(4) = (ptmp2(4)+ptmp2(3))/2d0
         k1(3) = (ptmp2(4)+ptmp2(3))/2d0
         k1(1) = 0d0
         k1(2) = 0d0

         k2(4) = -(ptmp2(4)-ptmp2(3))/2d0
         k2(3) = +(ptmp2(4)-ptmp2(3))/2d0
         k2(1) = 0d0
         k2(2) = 0d0

         c1 = 1d0
         c2 = -1d0
c         write(6,*)'k1 + k2 = ',k1(1)-k2(1),k1(2)-k2(2)
c     >        ,k1(3)-k2(3),k1(4)-k2(4)
c         write(6,*)'k1 = ',k1,dotprod(k1,k1)
c         write(6,*)'k2 = ',k2,dotprod(k2,k2)        
c         write(6,*)'c1 , c2 = ',c1,c2
c         write(6,*)'^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^'
      endif

      if ( (ptmp2(4)+ptmp2(3).lt.0d0).and.(ptmp2(4)-ptmp2(3).gt.0d0) )
     >     then
         k1(4) = -(ptmp2(4)+ptmp2(3))/2d0
         k1(3) = -(ptmp2(4)+ptmp2(3))/2d0
         k1(1) = 0d0
         k1(2) = 0d0

         k2(4) = (ptmp2(4)-ptmp2(3))/2d0
         k2(3) = -(ptmp2(4)-ptmp2(3))/2d0
         k2(1) = 0d0
         k2(2) = 0d0

         c1 = -1d0
         c2 = 1d0
c         write(6,*)'k1 + k2 = ',-k1(1)+k2(1),-k1(2)+k2(2)
c     >        ,-k1(3)+k2(3),-k1(4)+k2(4)
c         write(6,*)'k1 = ',k1,dotprod(k1,k1)
c         write(6,*)'k2 = ',k2,dotprod(k2,k2)        
c         write(6,*)'c1 , c2 = ',c1,c2
c         write(6,*)'^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^'
      endif

      if ( (ptmp2(4)+ptmp2(3).lt.0d0).and.(ptmp2(4)-ptmp2(3).lt.0d0) )
     >     then
         k1(4) = -(ptmp2(4)+ptmp2(3))/2d0
         k1(3) = -(ptmp2(4)+ptmp2(3))/2d0
         k1(1) = 0d0
         k1(2) = 0d0

         k2(4) = -(ptmp2(4)-ptmp2(3))/2d0
         k2(3) = +(ptmp2(4)-ptmp2(3))/2d0
         k2(1) = 0d0
         k2(2) = 0d0

         c1 = -1d0
         c2 = -1d0
c         write(6,*)'k1 + k2 = ',-k1(1)-k2(1),-k1(2)-k2(2)
c     >        ,-k1(3)-k2(3),-k1(4)-k2(4)
c         write(6,*)'k1 = ',k1,dotprod(k1,k1)
c         write(6,*)'k2 = ',k2,dotprod(k2,k2)        
c         write(6,*)'c1 , c2 = ',c1,c2
c         write(6,*)'^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^'
      endif

      call rotatex(ct,st,k1,ptmp1)
      call rotatez(cf,sf,ptmp1,k1)

      call rotatex(ct,st,k2,ptmp1)
      call rotatez(cf,sf,ptmp1,k2)

C Some tests

c      write(6,*)'P       = ',P
c      write(6,*)'K       = ',K,dotprod(K,K)/2d0
c      write(6,*)'k1 + k2 = ',c1*k1(1)+c2*k2(1),c1*k1(2)+c2*k2(2)
c     >                      ,c1*k1(3)+c2*k2(3),c1*k1(4)+c2*k2(4)
c      write(6,*)'k1 = ',k1,dotprod(k1,k1)
c      write(6,*)'k2 = ',k2,dotprod(k2,k2)
c      write(6,*)
C
C End code
C

      end
