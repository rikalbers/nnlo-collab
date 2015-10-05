       subroutine physics
       include 'declare.h'

       include 'common_coup.h'
       include 'common_masses.h'
       include 'common_flags.h'
       include 'common_print.h'
       include 'common_strf.h'
       include 'common_lha.h'
!Kaleu       include 'common_qcdrun.h'
       include 'common_qcdrun_kaleu.h' !Kaleu

C  CHARGE and ISOSPIN
       dimension q(4),t(4)
C  CKM MATRIX
       dimension zckm(3,3)

!Kaleu       logical onlyqcd
!Kaleu       logical withqcd

c      save init,root2,pi,zi,zero
       save
       data init/0/
        
c  -----------------------------------------------------------------------
c  start initialize
c  -----------------------------------------------------------------------
       if(init.eq.0)then

! CKM initialize
c[
       zckm(1:3,1:3)=0
       do i=1,3
       zckm(i,i)=1
       enddo
c]

       read*,onlyqcd,withqcd,irun

       write(*,'(a30,2l3,i3)') '         onlyqcd,withqcd,irun:'
     &                     ,onlyqcd,withqcd,irun

       root2=dsqrt(dnou(2))
       call mypi(pi)
       zi=dcmplx(dnou(0),dnou(1))
       zero=dcmplx(dnou(0),dnou(0))
       
       if(iunitary.eq.0)iuni=1
       if(iunitary.eq.1)iuni=0
       if(ihiggs.eq.0)ihi=0
       if(ihiggs.eq.1)ihi=1
      
       zgv3(31:35,31:35,31:35)=zero
       zgv4(31:35,31:35,31:35,31:35)=zero
       zgvffl(31:35,-12:-1,1:12)=zero
       zgvffr(31:35,-12:-1,1:12)=zero
       zgvvs(31:35,31:35,41:44)=zero
       zgsvv(41:44,31:35,31:35)=zero
       zgvvss(31:35,31:35,41:44,41:44)=zero
       zgssvv(41:44,41:44,31:35,31:35)=zero
       zgvss(31:35,41:44,41:44)=zero
       zgssv(41:44,41:44,31:35)=zero
       zgsffl(41:44,-12:-1,1:12)=zero
       zgsffr(41:44,-12:-1,1:12)=zero
       zgs3(41:44,41:44,41:44)=zero
       zgs4(41:44,41:44,41:44,41:44)=zero
       parmas(-12:44)=0
       parwid(-12:44)=0
       
       include 'constants.h'
       aqedup=alpha
       aqcdup=gqcd**2/4/pi

       write(*,'(a16,e23.16)') '         gfermi:',gfermi
       write(*,'(a16,e23.16)') '           rsw2:'
     &                       ,dnou(1)-parmas(33)**2/parmas(32)**2
       write(*,'(a16,e23.16)') '          alpha:',alpha
       write(*,'(a16,e15.8,e23.16)') '   alphas, gqcd:'
     &                              ,gqcd**2/4d0/pi,gqcd
       write(*,'(a16,2e15.8)') ' mass,width   Z:',parmas(32),parwid(32)
       write(*,'(a16,2e15.8)') ' mass,width   W:',parmas(33),parwid(33)
       write(*,'(a16,2e15.8)') ' mass,width   H:',parmas(41),parwid(41)
       write(*,'(a16,2e15.8)') ' mass,width nel:',parmas( 1),parwid( 1)
       write(*,'(a16,2e15.8)') ' mass,width  el:',parmas( 2),parwid( 2)
       write(*,'(a16,2e15.8)') ' mass,width   u:',parmas( 3),parwid( 3)
       write(*,'(a16,2e15.8)') ' mass,width   d:',parmas( 4),parwid( 4)
       write(*,'(a16,2e15.8)') ' mass,width nmu:',parmas( 5),parwid( 5)
       write(*,'(a16,2e15.8)') ' mass,width  mu:',parmas( 6),parwid( 6)
       write(*,'(a16,2e15.8)') ' mass,width   c:',parmas( 7),parwid( 7)
       write(*,'(a16,2e15.8)') ' mass,width   s:',parmas( 8),parwid( 8)
       write(*,'(a16,2e15.8)') ' mass,width nta:',parmas( 9),parwid( 9)
       write(*,'(a16,2e15.8)') ' mass,width  ta:',parmas(10),parwid(10)
       write(*,'(a16,2e15.8)') ' mass,width   t:',parmas(11),parwid(11)
       write(*,'(a16,2e15.8)') ' mass,width   b:',parmas(12),parwid(12)

c HHH coupling
       chi=1
       
       mm=33
       include 'compl_mass.h'
       zwm=zmas
       mm=32
       include 'compl_mass.h'
       zzm=zmas
       parmas(42)=parmas(32)
       parwid(42)=parwid(32)
       parmas(43)=parmas(33)
       parwid(43)=parwid(33)
       parmas(44)=parmas(34)
       parwid(44)=parwid(34)
       
       do i=1,12
       write(nunit1,'(i8,5x,e10.4,5x,e10.4)')i,parmas(i),parwid(i)
       enddo
       do i=31,35
       write(nunit1,'(i8,5x,e10.4,5x,e10.4)')i,parmas(i),parwid(i)
       enddo
       do i=41,44
       write(nunit1,'(i8,5x,e10.4,5x,e10.4)')i,parmas(i),parwid(i)
       enddo

       hm=parmas(41)
        
        do l=1,12
        parmas(-l)=parmas(l)
        parwid(-l)=parwid(l)
        enddo
       
       init=1
       endif
c  -----------------------------------------------------------------------
c  end initialize
c  -----------------------------------------------------------------------
       
c  -----------------------------------------------------------------------
c  start QCD couplings
c  -----------------------------------------------------------------------

       if(init.eq.1)then
       if(gqcd.eq.0)print*,'QCD COUPLING SET TO ZERO'
       if(gqcd.gt.0)print*,'QCD INCLUDED, g=',gqcd,gqcd**2/4/pi
       init=2
       elseif(init.ge.2.and.irun.eq.1) then
           call qcdscale(scale)
           gqcd=dsqrt(alphas(scale)*4*pi)
       if(init.eq.2)print*,'QCD INCLUDED, g=',gqcd,gqcd**2/4/pi
          init=3
       endif
        
       zgvffl(35,-3,3)=gqcd/root2
       zgvffr(35,-3,3)=gqcd/root2
       zgvffl(35,-4,4)=gqcd/root2
       zgvffr(35,-4,4)=gqcd/root2
       zgvffl(35,-7,7)=gqcd/root2
       zgvffr(35,-7,7)=gqcd/root2
       zgvffl(35,-8,8)=gqcd/root2
       zgvffr(35,-8,8)=gqcd/root2
       zgvffl(35,-11,11)=gqcd/root2
       zgvffr(35,-11,11)=gqcd/root2
       zgvffl(35,-12,12)=gqcd/root2
       zgvffr(35,-12,12)=gqcd/root2
       zgv3(35,35,35)=-gqcd/4*dsqrt(dnou(2))**3
       zgv4(35,35,35,35)= gqcd**2/8*dsqrt(dnou(2))**4
       
! IF ONLY QCD THEN
       if(onlyqcd)return
       
c  -----------------------------------------------------------------------
c  end QCD couplings
c  -----------------------------------------------------------------------
c  no running of couplings: next line does only for qcd
c  -----------------------------------------------------------------------
       if(init.eq.3)return
c  -----------------------------------------------------------------------

       zsinw2=dnou(1)-zwm**2/zzm**2
       zsinw =cdsqrt(zsinw2)
       zcosw =cdsqrt(dcmplx(dnou(1),dnou(0))-zsinw2)
       
       e=dsqrt(dnou(4)*pi*alpha)
       zgw=e/zsinw
       q(1)= dnou(0)
       q(2)=-dnou(1)
       q(3)= dnou(2)/dnou(3)
       q(4)=-dnou(1)/dnou(3)
       
       t(1)=+dnou(1)/dnou(2)
       t(2)=-dnou(1)/dnou(2)
       t(3)=+dnou(1)/dnou(2)
       t(4)=-dnou(1)/dnou(2)
        
       zgv3(31,33,34)=-e 
       zgv3(33,31,33)=-e 
       zgv3(34,34,31)=-e 

       zgv3(32,33,34)= zgw*zcosw
       zgv3(33,32,33)= zgw*zcosw
       zgv3(34,34,32)= zgw*zcosw
    
       zgv4(31,31,33,34)=-e**2                   !AAW+W-
       zgv4(31,32,33,34)= e**2*zcosw/zsinw       !AZW+W-
              
       zgv4(32,32,33,34)=-e**2*zcosw**2/zsinw**2  !ZZW+W-
       zgv4(32,31,33,34)= e**2*zcosw/zsinw        !ZAW+W-
       
       zgv4(33,34,33,33)= zgw*zgw           !/dnou(2)!W+W-W+W+
       zgv4(33,33,32,32)=-zgw*zgw*zcosw**2   !/dnou(2)!W+W+ZZ
       zgv4(33,33,31,32)= e**2*zcosw/zsinw        !W+W+AZ
       zgv4(33,33,31,31)=-e**2            !/dnou(2)!W+W+AA
       
       zgv4(34,33,34,34)= zgw*zgw           !/dnou(2)!W-W+W-W-
       zgv4(34,34,32,32)=-zgw*zgw*zcosw**2   !/dnou(2)!W-W-ZZ
       zgv4(34,34,31,32)= e**2*zcosw/zsinw        !W-W-AZ
       zgv4(34,34,31,31)=-e**2            !/dnou(2)!W-W-AA
         
       zgvffl(31,-2,2)=-e*q(2)
       zgvffl(31,-3,3)=-e*q(3)
       zgvffl(31,-4,4)=-e*q(4)

       zgvffl(31,-6,6)=zgvffl(31,-2,2)
       zgvffl(31,-7,7)=zgvffl(31,-3,3)
       zgvffl(31,-8,8)=zgvffl(31,-4,4)

       zgvffl(31,-10,10)=zgvffl(31,-6,6)
       zgvffl(31,-11,11)=zgvffl(31,-7,7)
       zgvffl(31,-12,12)=zgvffl(31,-8,8)
      
       zgvffr(31,-2,2)=-e*q(2)
       zgvffr(31,-3,3)=-e*q(3)
       zgvffr(31,-4,4)=-e*q(4)

       zgvffr(31,-6,6)=zgvffr(31,-2,2)
       zgvffr(31,-7,7)=zgvffr(31,-3,3)
       zgvffr(31,-8,8)=zgvffr(31,-4,4)

       zgvffr(31,-10,10)=zgvffr(31,-6,6)
       zgvffr(31,-11,11)=zgvffr(31,-7,7)
       zgvffr(31,-12,12)=zgvffr(31,-8,8)
        
       zgvffl(32,-1,1)=(t(1)-q(1)*zsinw**2)/zcosw*zgw 
       zgvffl(32,-2,2)=(t(2)-q(2)*zsinw**2)/zcosw*zgw 
       zgvffl(32,-3,3)=(t(3)-q(3)*zsinw**2)/zcosw*zgw 
       zgvffl(32,-4,4)=(t(4)-q(4)*zsinw**2)/zcosw*zgw 

       zgvffl(32,-5,5)=zgvffl(32,-1,1)
       zgvffl(32,-6,6)=zgvffl(32,-2,2)
       zgvffl(32,-7,7)=zgvffl(32,-3,3)
       zgvffl(32,-8,8)=zgvffl(32,-4,4)
      
       zgvffl(32,-9,9)  =zgvffl(32,-1,1)
       zgvffl(32,-10,10)=zgvffl(32,-2,2)
       zgvffl(32,-11,11)=zgvffl(32,-3,3)
       zgvffl(32,-12,12)=zgvffl(32,-4,4)
      
       zgvffr(32,-1,1)=-q(1)*zsinw**2/zcosw*zgw
       zgvffr(32,-2,2)=-q(2)*zsinw**2/zcosw*zgw 
       zgvffr(32,-3,3)=-q(3)*zsinw**2/zcosw*zgw
       zgvffr(32,-4,4)=-q(4)*zsinw**2/zcosw*zgw

       zgvffr(32,-5,5)=zgvffr(32,-1,1)
       zgvffr(32,-6,6)=zgvffr(32,-2,2)
       zgvffr(32,-7,7)=zgvffr(32,-3,3)
       zgvffr(32,-8,8)=zgvffr(32,-4,4)
      
       zgvffr(32,-9,9)  =zgvffr(32,-1,1)
       zgvffr(32,-10,10)=zgvffr(32,-2,2)
       zgvffr(32,-11,11)=zgvffr(32,-3,3)
       zgvffr(32,-12,12)=zgvffr(32,-4,4)
         
       zgvffl(34,-1,2)=zgw/root2
       zgvffl(34,-3,4)=zgw/root2*dconjg(zckm(1,1))
       zgvffl(34,-3,8)=zgw/root2*dconjg(zckm(1,2))
       zgvffl(34,-3,12)=zgw/root2*dconjg(zckm(1,3))

       zgvffl(34,-5,6)=zgvffl(34,-1,2)
       zgvffl(34,-7,4)=zgw/root2*dconjg(zckm(2,1))
       zgvffl(34,-7,8)=zgw/root2*dconjg(zckm(2,2))
       zgvffl(34,-7,12)=zgw/root2*dconjg(zckm(2,3))
       
       zgvffl(34,-9,10)=zgvffl(34,-1,2) 
       zgvffl(34,-11,4)=zgw/root2*dconjg(zckm(3,1))
       zgvffl(34,-11,8)=zgw/root2*dconjg(zckm(3,2))
       zgvffl(34,-11,12)=zgw/root2*dconjg(zckm(3,3))
       
       zgvffl(33,-2,1)=zgw/root2
       zgvffl(33,-4,3)=zgw/root2*zckm(1,1) 
       zgvffl(33,-8,3)=zgw/root2*zckm(2,1) 
       zgvffl(33,-12,3)=zgw/root2*zckm(3,1) 
       
       zgvffl(33,-6,5)=zgvffl(33,-2,1)
       zgvffl(33,-4,7)=zgw/root2*zckm(1,2)
       zgvffl(33,-8,7)=zgw/root2*zckm(2,2)
       zgvffl(33,-12,7)=zgw/root2*zckm(3,2)
       
       zgvffl(33,-10,9)=zgvffl(33,-2,1)
       zgvffl(33,-4,11)=zgw/root2*zckm(1,3)
       zgvffl(33,-8,11)=zgw/root2*zckm(2,3)
       zgvffl(33,-12,11)=zgw/root2*zckm(3,3)
      
       zgsvv(41,33,34)= zgw*zwm*ihi                               !HW+W-
       zgsvv(41,32,32)= zgw*zwm/zcosw**2 *ihi                !/dnou(2)!HZZ
      
       zgsvv(43,33,32)=-e*zwm*zsinw/zcosw*iuni
       zgsvv(43,33,31)=-e*zwm*iuni
        
       zgsvv(44,34,32)=-e*zwm*zsinw/zcosw*iuni
       zgsvv(44,34,31)=-e*zwm*iuni
     
       zgvvs(31,33,44)=-e*zwm*iuni
       zgvvs(31,34,43)=-e*zwm*iuni
       
       zgvvs(32,32,41)= zgw*zwm/zcosw**2*ihi                        !ZZH
       zgvvs(32,33,44)=-e*zwm*zsinw/zcosw*iuni
       zgvvs(32,34,43)=-e*zwm*zsinw/zcosw*iuni
      
       zgvvs(33,33,41)= zgw*zwm*ihi                                !W+W+H
       zgvvs(33,31,43)=-e*zwm*iuni
       zgvvs(33,32,43)=-e*zwm*zsinw/zcosw*iuni
      
       zgvvs(34,34,41)= zgw*zwm*ihi                                !W-W-H
       zgvvs(34,31,44)=-e*zwm*iuni
       zgvvs(34,32,44)=-e*zwm*zsinw/zcosw*iuni
      
       zgvvss(34,34,41,41)=zgw**2/dnou(2)*ihi                    !/dnou(2)!W-W-HH
       zgvvss(34,34,42,42)=zgw**2/dnou(2)*iuni                   !/dnou(2)!W-W-XX
       zgvvss(34,34,43,44)=zgw**2/dnou(2)*iuni
       zgvvss(34,31,44,41)=-e**2/(dnou(2)*zsinw)*iuni
       zgvvss(34,31,44,42)=-zi*e**2/(dnou(2)*zsinw)*iuni
       zgvvss(34,32,44,41)=-e**2/(dnou(2)*zcosw)*iuni
       zgvvss(34,32,44,42)=-zi*e**2/(dnou(2)*zcosw)*iuni

       zgvvss(33,33,41,41)=zgw**2/dnou(2)*ihi                    !/dnou(2)!W+W+HH
       zgvvss(33,33,42,42)=zgw**2/dnou(2)*iuni                   !/dnou(2)!W+W+XX
       zgvvss(33,33,43,44)=zgw**2/dnou(2)*iuni
       zgvvss(33,31,43,41)=-e**2/(dnou(2)*zsinw)*iuni
       zgvvss(33,31,43,42)=zi*e**2/(dnou(2)*zsinw)*iuni
       zgvvss(33,32,43,41)=-e**2/(dnou(2)*zcosw)*iuni
       zgvvss(33,32,43,42)=zi*e**2/(dnou(2)*zcosw)*iuni

       zgvvss(32,32,41,41)=zgw**2/(dnou(2)*zcosw**2)*ihi          !/dnou(2)!ZZHH
       zgvvss(32,32,42,42)=zgw**2/(dnou(2)*zcosw**2)*iuni         !/dnou(2)!ZZXX
       zgvvss(32,32,43,44)=zgw**2*(zsinw**2-zcosw**2)**2
     . /(dnou(2)*zcosw**2)*iuni
       zgvvss(32,31,43,44)=e**2*(zsinw**2-zcosw**2)/(zsinw*zcosw)*iuni
       zgvvss(32,33,44,41)=-e**2/(dnou(2)*zcosw)*iuni
       zgvvss(32,34,43,41)=-e**2/(dnou(2)*zcosw)*iuni
       zgvvss(32,33,44,42)=-zi*e**2/(dnou(2)*zcosw)*iuni
       zgvvss(32,34,43,42)=zi*e**2/(dnou(2)*zcosw)*iuni

       zgvvss(31,31,43,44)=dnou(2)*e**2*iuni
       zgvvss(31,32,43,44)=e**2*(zsinw**2-zcosw**2)/(zsinw*zcosw)*iuni
       zgvvss(31,33,44,41)=-e**2/(dnou(2)*zsinw)*iuni
       zgvvss(31,34,43,41)=-e**2/(dnou(2)*zsinw)*iuni      
       zgvvss(31,33,44,42)=-zi*e**2/(dnou(2)*zsinw)*iuni    
       zgvvss(31,34,43,42)= zi*e**2/(dnou(2)*zsinw)*iuni          

       zgvss(31,43,44)=-e*iuni 

       zgvss(32,42,41)=-zi*zgw/(dnou(2)*zcosw)*iuni 
       zgvss(32,43,44)=-zgw*(zsinw**2-zcosw**2)/(dnou(2)*zcosw)*iuni 
      
       zgvss(33,43,41)= zgw/dnou(2)*iuni 
       zgvss(33,43,42)=-zi*zgw/dnou(2)*iuni 
     
       zgvss(34,44,41)=-zgw/dnou(2)*iuni 
       zgvss(34,44,42)=-zi*zgw/dnou(2)*iuni 
       
       zgssv(41,42,32)= zi*zgw/(dnou(2)*zcosw)*iuni  
       zgssv(41,43,34)=-zgw/dnou(2)*iuni 
       zgssv(41,44,33)= zgw/dnou(2)*iuni  
 
       zgssv(42,41,32)=-zi*zgw/(dnou(2)*zcosw)*iuni               !XHZ 
       zgssv(42,44,33)= zi*zgw/dnou(2)*iuni  
       zgssv(42,43,34)= zi*zgw/dnou(2)*iuni  
       
       zgssv(43,43,31)= e*iuni 
       zgssv(43,43,32)= zgw*(zsinw**2-zcosw**2)/(dnou(2)*zcosw)*iuni 
       zgssv(43,41,33)=-zgw/dnou(2)*iuni 
       zgssv(43,42,33)=-zi*zgw/dnou(2)*iuni 

       zgssv(44,44,31)=-e*iuni 
       zgssv(44,44,32)=-zgw*(zsinw**2-zcosw**2)/(dnou(2)*zcosw)*iuni 
       zgssv(44,41,34)= zgw/dnou(2)*iuni 
       zgssv(44,42,34)=-zi*zgw/dnou(2)*iuni 
       
       zgssvv(41,41,32,32)= zgw**2/(dnou(2)*zcosw**2)*ihi      !/dnou(2)!HHZZ
       zgssvv(41,41,33,34)= zgw**2/dnou(2)*ihi                      !HHW+W-
       zgssvv(41,43,34,31)=-e**2/(dnou(2)*zsinw)*iuni 
       zgssvv(41,44,33,31)=-e**2/(dnou(2)*zsinw)*iuni 
       zgssvv(41,43,34,32)=-e**2/(dnou(2)*zcosw)*iuni 
       zgssvv(41,44,33,32)=-e**2/(dnou(2)*zcosw)*iuni 
       
       zgssvv(42,42,32,32)= zgw**2/(dnou(2)*zcosw**2)*iuni 
       zgssvv(42,42,33,34)= zgw**2/dnou(2)*iuni 
       zgssvv(42,43,34,31)= zi*e**2/(dnou(2)*zsinw)*iuni 
       zgssvv(42,44,33,31)=-zi*e**2/(dnou(2)*zsinw)*iuni 
       zgssvv(42,43,34,32)= zi*e**2/(dnou(2)*zcosw)*iuni 
       zgssvv(42,44,33,32)=-zi*e**2/(dnou(2)*zcosw)*iuni 

       zgssvv(43,43,31,31)= dnou(2)*e**2*iuni 
       zgssvv(43,43,32,31)= e**2*(zsinw**2-zcosw**2)/(zsinw*zcosw)*iuni 
       zgssvv(43,43,32,32)= zgw**2*(zsinw**2-zcosw**2)**2
     . /(dnou(2)*zcosw**2)
     . *iuni 
       zgssvv(43,43,33,34)= zgw**2/dnou(2)*iuni 
       zgssvv(43,41,33,31)=-e**2/(dnou(2)*zsinw)*iuni 
       zgssvv(43,42,33,31)=-zi*e**2/(dnou(2)*zsinw)*iuni 
       zgssvv(43,41,33,32)=-e**2/(dnou(2)*zcosw)*iuni 
       zgssvv(43,42,33,32)=-zi*e**2/(dnou(2)*zcosw)*iuni 
       
       zgssvv(44,44,31,31)=dnou(2)*e**2*iuni 
       zgssvv(44,44,32,31)=e**2*(zsinw**2-zcosw**2)/(zsinw*zcosw)*iuni 
       zgssvv(44,44,32,32)=zgw**2*(zsinw**2-zcosw**2)**2
     . /(dnou(2)*zcosw**2)
     . *iuni 
       zgssvv(44,44,33,34)=zgw**2/dnou(2)*iuni 
       zgssvv(44,41,34,31)=-e**2/(dnou(2)*zsinw)*iuni 
       zgssvv(44,42,34,31)= zi*e**2/(dnou(2)*zsinw)*iuni 
       zgssvv(44,41,34,32)=-e**2/(dnou(2)*zcosw)*iuni 
       zgssvv(44,42,34,32)= zi*e**2/(dnou(2)*zcosw)*iuni 
      
       do i=1,12
       zgsffl(41,-i,i)=-e/(dnou(2)*zsinw)/zwm*parmas(i)*ihi
       zgsffr(41,-i,i)=-e/(dnou(2)*zsinw)/zwm*parmas(i)*ihi
       enddo
       
       do i=1,12
       j=mod(i-1,4)+1
       zgsffl(42,-i,i)=-zi*e/(dnou(2)*zsinw)
     . *(dnou(2)*t(j))/zwm*parmas(i)*iuni 
       zgsffr(42,-i,i)= zi*e/(dnou(2)*zsinw)
     . *(dnou(2)*t(j))/zwm*parmas(i)*iuni 
       enddo
       
       zgsffr(44,-1,2)= -zgw*parmas(2)/(root2*zwm)*iuni       
       zgsffr(44,-5,6)= -zgw*parmas(6)/(root2*zwm)*iuni 
       zgsffr(44,-9,10)=-zgw*parmas(10)/(root2*zwm)*iuni 

       zgsffl(44,-1,2)= zero
       zgsffl(44,-5,6)= zero
       zgsffl(44,-9,10)=zero
       
       zgsffr(44,-3,4)=  -zgw*parmas(4)/(root2*zwm)*iuni 
       zgsffr(44,-7,8)=  -zgw*parmas(8)/(root2*zwm)*iuni 
       zgsffr(44,-11,12)=-zgw*parmas(12)/(root2*zwm)*iuni 
       
       zgsffl(44,-3,4)=  zgw*parmas(3)/(root2*zwm)*iuni 
       zgsffl(44,-7,8)=  zgw*parmas(7)/(root2*zwm)*iuni 
       zgsffl(44,-11,12)=zgw*parmas(11)/(root2*zwm)*iuni 
       
       zgsffl(43,-2,1) =-zgw*parmas(2)/(root2*zwm)*iuni 
       zgsffl(43,-6,5) =-zgw*parmas(6)/(root2*zwm)*iuni 
       zgsffl(43,-10,9)=-zgw*parmas(10)/(root2*zwm)*iuni 
       
       zgsffr(43,-2,1)=zero
       zgsffr(43,-6,5)=zero
       zgsffr(43,-10,9)=zero
       
       zgsffl(43,-4,3)=  -zgw*parmas(4)/(root2*zwm)*iuni 
       zgsffl(43,-8,7)=  -zgw*parmas(8)/(root2*zwm)*iuni 
       zgsffl(43,-12,11)=-zgw*parmas(12)/(root2*zwm)*iuni 
       
       zgsffr(43,-4,3)=   zgw*parmas(3)/(root2*zwm)*iuni 
       zgsffr(43,-8,7)=   zgw*parmas(7)/(root2*zwm)*iuni 
       zgsffr(43,-12,11)= zgw*parmas(11)/(root2*zwm)*iuni 
       
       zgs3(41,41,41)=(-dnou(3)*zgw/dnou(2))*hm**2/zwm*ihi        !/dnou(2)!HHH
     .               *chi
       zgs3(41,42,42)=(-zgw/dnou(2))*hm**2/zwm*iuni            !/dnou(2)!HXX
       zgs3(41,43,44)=(-zgw/dnou(2))*hm**2/zwm*iuni 
       
       zgs3(42,42,41)=(-zgw/dnou(2))*hm**2/zwm*iuni 
       
       zgs3(43,43,41)=(-zgw/dnou(2))*hm**2/zwm*iuni 
       
       zgs3(44,44,41)=(-zgw/dnou(2))*hm**2/zwm*iuni 
       
       zgs4(41,41,41,41)=(-dnou(3)/dnou(4))*zgw**2*hm**2/zwm**2*ihi   !/6.d0!HHHH
       zgs4(41,41,42,42)=(-dnou(1)/dnou(4))*zgw**2*hm**2/zwm**2*iuni  !/dnou(2)!HHXX
       zgs4(41,41,43,44)=(-dnou(1)/dnou(4))*zgw**2*hm**2/zwm**2*iuni 

       zgs4(42,42,41,41)=(-dnou(1)/dnou(4))*zgw**2*hm**2/zwm**2*iuni   !/dnou(2)!XXHH
       zgs4(42,42,42,42)=(-dnou(3)/dnou(4))*zgw**2*hm**2/zwm**2*iuni   !/6.d0!XXXX 
       zgs4(42,42,43,44)=(-dnou(1)/dnou(4))*zgw**2*hm**2/zwm**2*iuni 

       zgs4(43,43,41,41)=(-dnou(1)/dnou(4))*zgw**2*hm**2/zwm**2*iuni   !/dnou(2)!F+F+HH
       zgs4(43,43,42,42)=(-dnou(1)/dnou(4))*zgw**2*hm**2/zwm**2*iuni   !/dnou(2)!F+F+XX 
       zgs4(43,43,44,44)=(-dnou(1)/dnou(2))*zgw**2*hm**2/zwm**2*iuni   !/dnou(2)!F+F+F-F-

       zgs4(44,44,41,41)=(-dnou(1)/dnou(4))*zgw**2*hm**2/zwm**2*iuni   !/dnou(2)!F-F-HH
       zgs4(44,44,42,42)=(-dnou(1)/dnou(4))*zgw**2*hm**2/zwm**2*iuni   !/dnou(2)!F-F-XX  
       zgs4(44,44,43,43)=(-dnou(1)/dnou(2))*zgw**2*hm**2/zwm**2*iuni   !/dnou(2)!F-F-F+F+
       
       iqed=0
       if(iqed.eq.1)then
       
       zgv3(32:34,32:34,32:34)=zero
       zgv4(32:34,32:34,32:34,32:34)=zero
       zgvffl(32:34,-12:-1,1:12)=zero
       zgvffr(32:34,-12:-1,1:12)=zero
       zgvvs(31:35,31:35,41:44)=zero
       zgsvv(41:44,31:35,31:35)=zero
       zgvvss(31:35,31:35,41:44,41:44)=zero
       zgssvv(41:44,41:44,31:35,31:35)=zero
       zgvss(31:35,41:44,41:44)=zero
       zgssv(41:44,41:44,31:35)=zero
       zgsffl(41:44,-12:-1,1:12)=zero
       zgsffr(41:44,-12:-1,1:12)=zero
       zgs3(41:44,41:44,41:44)=zero
       zgs4(41:44,41:44,41:44,41:44)=zero
       
       endif

       include 'special/physics.h'

       return
       end   
