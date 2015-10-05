       subroutine readcuts_kaleu
       include 'declare.h'
       include 'common_masses.h'
       include 'common_int.h'
       include 'common_cuts.h'
       character*24 file
     
       cutoff=1.0d-3
       print*,'WARNING CUTOFF SET:',cutoff
       ptc(3:20)=cutoff
       drc(3:20,3:20)=0
       etac(3:20)=20
       ec(3:20)=cutoff
       c1(3:20)=1
       c2(3:20)=1
       cc(3:20,3:20)=1
       gmas(3:20,3:20)=cutoff
       pi=dacos(-dnou(1))
       
c -- PP cuts
       include 'cuts.h'
c --

       do i=3,n
       if(ifl(i).lt.13)then
        i1=iabs(mod(ifl(i),4))
       else
        i1=ifl(i)
       endif
        if(i1.eq.35)i1=0
c- bottom quarks
        if(iabs(ifl(i)).eq.12)i1=12
c- top quarks
        if(iabs(ifl(i)).eq.11)i1=11
        goto 111
c -------
c special for mc4lhc 2003
c[
c       if(ifl(i).eq.11)i1=-1
c       if(ifl(i).eq.-11)i1=-1
c       if(ifl(i).eq.12)i1=-1
c       if(ifl(i).eq.-12)i1=-1
c]

 111    continue
c   photons
       if(i1.eq.31)then
        ptc(i)=ptg
        etac(i)=etag
       endif
c   neutrinos
       if(i1.eq.1)then
       endif
c   ch. leptons
       if(i1.eq.2)then
c      if(i1.le.2)then
        ptc(i)=ptl
        etac(i)=etal
       endif
c   quarks
       if(i1.eq.3.or.i1.eq.0)then
        ptc(i)=ptq
        etac(i)=etaq
       endif
c  bottom 
       if(i1.eq.12)then
        ptc(i)=ptbottom
        etac(i)=etabottom
       endif
c   top
       if(i1.eq.11)then
        ptc(i)=pttop
        etac(i)=etatop
       endif
       enddo
    
       do i=3,n-1
       if(ifl(i).lt.13)then
        i1=iabs(mod(ifl(i),4))
       else
        i1=ifl(i)
       endif
        if(i1.eq.35)i1=0
c- bottom quarks                                           
        if(iabs(ifl(i)).eq.12)i1=12
c- top quarks 
        if(iabs(ifl(i)).eq.11)i1=11
        goto 222
c -------
c special for mc4lhc 2003
c[
c       if(ifl(i).eq.11)i1=-1
c       if(ifl(i).eq.-11)i1=-1
c       if(ifl(i).eq.12)i1=-1
c       if(ifl(i).eq.-12)i1=-1
c]

 222    continue
        do j=i+1,n
         if(ifl(j).lt.13)then
          j1=iabs(mod(ifl(j),4))
         else
          j1=ifl(j)
         endif
         if(j1.eq.35)j1=0
c- bottom quarks 
         if(iabs(ifl(j)).eq.12)j1=12
c- top quarks
         if(iabs(ifl(j)).eq.11)j1=11
         goto 333

c -------
c special for mc4lhc 2003
c[
c       if(ifl(j).eq.11)j1=-1
c       if(ifl(j).eq.-11)j1=-1
c       if(ifl(j).eq.12)j1=-1
c       if(ifl(j).eq.-12)j1=-1
c]
 333    continue

c  n-n
         if(i1.eq.1.or.j1.eq.1)then
c  l-q
         elseif(i1.eq.2.and.(j1.eq.3.or.j1.eq.0.or.j1.eq.12))then
         drc(i,j)=drlq
c  l-l
         elseif(i1.eq.2.and.j1.eq.2)then
         drc(i,j)=drll
c  q-l
         elseif(j1.eq.2.and.(i1.eq.3.or.i1.eq.0.or.i1.eq.12))then
         drc(i,j)=drlq
c  g-x  photon-charged
         elseif(i1.eq.31.and.j1.ne.1)then
         drc(i,j)=drgx
c  x-g  charged-photon
         elseif(j1.eq.31.and.i1.ne.1)then
         drc(i,j)=drgx
c  q-q
         elseif((i1.eq.3.or.i1.eq.0).and.
     .          (j1.eq.3.or.j1.eq.0))then
         drc(i,j)=drqq
         gmas(i,j)=max(gqq,gmas(i,j))
c  q-b
         elseif((i1.eq.3.or.i1.eq.0).and.
     .          (j1.eq.12))then
         drc(i,j)=drqb
         gmas(i,j)=max(gqb,gmas(i,j))
c  b-q
         elseif((j1.eq.3.or.j1.eq.0).and.
     .          (i1.eq.12))then
         drc(i,j)=drqb
         gmas(i,j)=max(gqb,gmas(i,j))
c  b-b
         elseif((i1.eq.12).and.
     .          (j1.eq.12))then
         drc(i,j)=drbb
         gmas(i,j)=max(gbb,gmas(i,j))
c
         endif
        enddo
       enddo
   
       do i=3,n
c       ec(i)=max(ptc(i),parmas(ifl(i)))
        if(ptc(i).gt.0)then
         ec(i)=max(ec(i),ptc(i))
        else
         ec(i)=max(ec(i),parmas(ifl(i)))
        endif
        r=exp(2*etac(i))
        c1(i)=(r-1)/(r+1)
        c2(i)=c1(i)
       enddo 
   
       do i=3,n-1
        do j=i+1,n
c        ri=sqrt(1-c1(i)*c1(i))
c        rj=sqrt(1-c1(j)*c1(j))
c        cc(i,j)=c1(i)*c1(j)+ri*rj*cos(drc(i,j))
         cc(i,j)=cos(drc(i,j))
        enddo
       enddo

       do i=3,n-1
        do j=i+1,n
         gmas(i,j)=max(gmas(i,j),
     &             dsqrt( 2*ptc(i)*ptc(j)*(1-cos(drc(i,j))) ) )
        enddo
       enddo

       do l=4,20
        do k=3,l-1
          drc(l,k)= drc(k,l)
           cc(l,k)=  cc(k,l)
           gmas(l,k)=  gmas(k,l)
        enddo
       enddo
  
       print*,'---------------------------------------------------'
       print*,'        the cuts        '
       do i=3,n
       print*,'pt     of  ',i,'   particle   ',ptc(i)
       print*,'energy of  ',i,'   particle   ',ec(i)
       enddo
       do i=3,n
       print*,'rapidity of  ',i,'   particle   ',etac(i)
       enddo
       do i=3,n
       print*,'cos-beam1 of ',i,'   particle   ',c1(i)
       enddo
       do i=3,n
       print*,'cos-beam2 of ',i,'   particle   ',c2(i)
       enddo

       do i=3,n-1
        do j=i+1,n
         if(i.ne.j)then
         print*,'DR     ',i,'  with  ',j,drc(i,j)
         print*,'cos of ',i,'  with  ',j,cc(i,j)
         endif
        enddo
       enddo

       do i=3,n-1
        do j=i+1,n
         if(i.ne.j)
     .   print*,'mass of ',i,'  with  ',j,gmas(i,j)
        enddo
       enddo
       print*,'---------------------------------------------------'

       return
       end
