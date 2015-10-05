       subroutine integration_kaleu( w0 ,ikaleu,nbatch,nstep )
       include 'declare.h'

       include 'common_int.h'
!       include 'common_mom.h'
!       include 'common_masses.h'
       include 'common_feyn.h'
       include 'common_phegas_mom.h'
       include 'common_new.h'
       include 'common_print.h'
       include 'common_norm.h'
       include 'common_strf.h'
       include 'common_unweight.h'
!       include 'common_lha.h'
       include 'common_unw.h'
       include 'common_qcdrun_kaleu.h'
       include 'common_debug.h'
       include 'common_onep.h'

       integer ipdt(-100:100)

       dimension wmax(10)
 
       logical onep ,discard,initphase
       
c WRITE OUT
c START
       logical lunwei
c lunwei: for unweightinga   !1/3
       data lunwei/.false./  !2/3
c      data lunwei/.true./   !3/3
       data nunwei/1000/   !3/3
       logical lwmax,lwri
       data lwmax/.false./
       integer nwri,nwmax,nwarn
       data nwri/0/,nwmax/0/,nwarn/0/
       double precision umax,umax1
       data umax/0/
c END

       save init,initphase
       data init/0/,initphase/.true./


       print*,'NG=',ng
       if(ng.eq.0)then
          print*,'no FG contributions found: BYE'
          stop
       endif
c --------------------------------------------------------------
c         START OF INITIALIZATION
c---------------------------------------------------------------
       if(init.eq.0)then
          
        idebug=0
        print*,'FOR DEBUGGING',idebug
        onep=.false.
     
c unweighting
      read(*,*,end=1000)lunwei,nunwei,nwri_tot
c     read*,lunwei,nunwei,nwri_tot
      print*,'UNWEIGHTING IS', lunwei

          inonzero=0
          read(*,*,end=1000)ncha
          read*,nmc
          read*,idummy1,idummy2,xdummy3,idummy4,idummy5,idummy6
          print*,'nmc:',nmc
          print*,'nbatch,nstep:',nbatch,nstep
          wmax(1:10)=0
          wmemax=0
          wme=1
          imc=0
          init=1
       endif

       call kaleu2helac_init( ikaleu,nbatch,nstep ) ! This routine is defined below

c --------------------------------------------------------------
c         END  OF INITIALIZATION
c---------------------------------------------------------------
      call cpu_time( tstart)
      call vtime0(tcpu0,tclock0)
 4     continue
       ipass=0
       iev=0
       do while(iev.le.nmc)  ![iev
        iev=iev+1
          w=0
          wme=0

          call kaleu2helac_gnrt( discard ) ! This routine is defined below
          if (discard) goto 1

          call cuts(icut)
          if(icut.eq.0) goto 1

cdef include cuts from possible decays
c      call decays(0,icut,w)
c      if(icut.eq.0)goto1

          call helac_kaleu_wght( wpsp_t )
          if (wpsp_t.eq.0d0) goto 2

          inonzero=inonzero+1
          if(irun.eq.1)call physics

c         call cpu_time(tstart)
c         npoi=10000
c         do i=1,npoi
          call helac_master(wme)
!          wme = 1d0 !DEBUG
!          wme = wme**(0.5d0) !DEBUG
c         enddo
c         call cpu_time(tstop)
c     print*,'Elapsed time', (tstop-tstart)*1000/npoi
c         stop

cx     wrest=1
cx     call test_func(wme)
          w=wpsp_t*wme*wrest
          wsf=1
          if(istruc.eq.1)call strf(wsf)
          w=w*wsf
          ipass=ipass+1

          call helac_kaleu_collect( w )

       wr=1
       include 'adapt1.h'
       include 'ktreweight.h'

       w1=w

cdef include cuts just for histograms
c      call decays(1,icut,w)

          if(wme.gt.wmemax)then
             wmemax=wme
          endif
          if(w.gt.wmax(1))then
             do i=10,2,-1
                wmax(i)=wmax(i-1)
             enddo
             if(wmax(1).gt.0)then
                if(w/wmax(1).gt.10)then
                  if(idebug.eq.1)then
                   write(nunit2,*)'WARNING ABOUT WEIGHT'
                   write(nunit2,*)iev
                   write(nunit2,*)w,wme,wpsp_t
                   write(nunit2,*)w,wmax
                   write(nunit2,*)(i,pmom(i,1:4),i=1,n)
                   write(nunit2,*)'--------------------'
                  endif
                endif
             endif
             wmax(1)=w
          endif

C HISTOGRAMING

c          call tohisto(w)
      include 'toh.h'


c WRITE OUT EVENTS
c START
c --------------
          if(lunwei)then
c --------------
             if(.not.lwmax.and..not.initphase)then
                if(umax.lt.w1)then
                   umax=w1
c                  call put_unwei(umax)
                endif
                nwmax=nwmax+1
                umax1=umax
                call put_unwei(w1)
             endif
c            if(nwmax.eq.nopt_step) then
             if(nwmax.eq.nunwei) then
                nwmax=nwmax+1
                lwmax=.true.
                call get_unwei(umax)
                write(nunit1,*)'START UNWEIGHTING',iev,umax,nwmax,nunwei
             endif
             if(lwmax)then
                if(w1.gt.umax1)then
                   umax1=w1
                   write(nunit1,*)'iev,umax1,umax',iev,umax1,umax
                endif
                lwri=.false.
                if(umax*rnmy(0).lt.w1)lwri=.true.
                if(umax.lt.w1)nwarn=nwarn+1
               a=0.d0
               b=10*umax
               call histo1(199,100,a,b,w1,1.d0)

               if(lwri)then                   
                nwri=nwri+1
                idprup=81
                xwgtup=1 !w1*10**3 !1
                call qcdscale(scale1)
                scalup=scale1
                write(nunit3)N,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP

c               call cmstolab
                include 'pdt.h'
                do i=1,n
                   idup=ifl(i)
                   idup=ipdt(ifl(i))
                   istup=-io(i)

                   imothup1=0
                   imothup2=0
                   if(i.gt.2.and.i.le.n) then
                      imothup1=1
                      imothup2=2
                   endif

                   if(io(i).eq.1) then
                      icol1=icol_un(i,1)+100
                      if(icol1.eq.100)icol1=0
                      icol2=icol_un(i,2)+100
                      if(icol2.eq.100)icol2=0
                   else
                      icol1=icol_un(i,2)+100
                      if(icol1.eq.100)icol1=0
                      icol2=icol_un(i,1)+100
                      if(icol2.eq.100)icol2=0
                   endif

                   px=pmom(i,1)
                   py=pmom(i,2)
                   pz=pmom(i,3)
                   p0=pmom(i,4)
                   pmass=pmom(i,5)
                   vtime=0
                   vspin=9

                   write(nunit3)idup,istup,imothup1,imothup2,icol1,icol2
     &                  ,px,py,pz,p0
     &                  ,pmass,vtime,vspin

                enddo  
               endif
             endif
c     rewind(nunit3)

             if(nwri.eq.nwri_tot) go to 5

c --------------
          endif
c --------------
c END
 1        continue  

          call bookin(1,w)
          call bookin(2,wme)
c         if(w.gt.0)call bookin(3,dlog(w))
          if(lwri)call bookin(3,w)
 2        continue   
          call geti(1,0,w0)
          if(mod(iev,10000).eq.0)then
             call errest(1,x1,y1,0)
             print'(a6,2d15.6,3i10)','sigma=',
     .            x1,sqrt(y1)/x1,inonzero,int(w0),iev
             call errest(3,xx1,yy1,0)
             if(yy1.gt.0)print'(a7,2d15.6)','<w=/=0>',
     .            xx1,sqrt(yy1)/dabs(xx1)
             if(nwri.gt.0)print*,nwri
             write(17,*)iev,x1,sqrt(y1)
             print*,'----------------------------'
          endif

         if (initphase.and.ipass.eq.nbatch*nstep) then
           initphase = .false. ! necessary to avoid repetetion
           print*,''
           print'(a7,2d15.6,3i10)',' sigma=',
     .            x1,sqrt(y1)/x1,inonzero,int(w0),iev
           print*,'Throwing away estimates'
     &           ,', keeping optimized distributions'
     &           ,', and restarting MC'
           print*,''
           inonzero = 0
           call clear(1)
           call clear(2)
         endif

       enddo ![iev

      call vtime0(tcpu1,tclock1)

 5      continue
       call geti(1,0,w0)
       print*,'total number of acc. and gen. events:',ipass,iev
       print*,'number of events used for estimates :',inonzero,int(w0)
       print*,'                cpu-time, clock-time:'
     &        ,tcpu1-tcpu0,tclock1-tclock0
       call errest(1,x1,y1,1)

       call geti(3,0,ss0)
       print*,'lwri: points have used',ss0

        read*,alimit,nlimit

        read*,IDBMUP1,IDBMUP2,EBMUP1,EBMUP2,IDWTUP,NPRUP
        print*,IDBMUP1,IDBMUP2,EBMUP1,EBMUP2,IDWTUP,NPRUP
       WRITE(200,5100) IDBMUP1,IDBMUP2,EBMUP1,EBMUP2,
     &      iPDFGUP1,iPDFGUP2,iPDFSUP1,iPDFSUP2,IDWTUP,NPRUP
c    &             0,       0,   10042,   10042,IDWTUP,NPRUP
       
       WRITE(200,5200) x1*10d0**3,dsqrt(y1)*10d0**3,1.d0,81

       if(x1.gt.0)print*,'% error:',dsqrt(y1)/x1*100
       write(16,*)x1,sqrt(y1)
       if(wmax(1).gt.0)
     . print*,'<w>/w_max,w_max',x1/wmax(1),wmax(1)
       if(wmax(10).gt.0)
     . print*,'<w>/w_max,w_max',x1/wmax(10),wmax(10)
       print'(10(E14.6,2X))',wmax(1:10)
       call errest(2,x1,y1,0)
       if(wmemax.gt.0)
     . print*,'<me>/memax,memax',x1/wmemax,wmemax
     
       print*,'number of w=1 events',nwri
       print*,'number of w>1 events',nwarn
       print*,'maximum weight used for un:vs',umax,umax1

 5100 FORMAT(1P,2I8,2E14.6,6I6)
 5200 FORMAT(1P,3E14.6,I6)

      call histo2(199,199,1,1.)

1000   return
       end


      subroutine kaleu2helac_init( ikaleu,nbatch,nstep )
!* ********************************************************************
!* ********************************************************************
      include 'declare.h'
!       
      include 'common_flags.h'
      include 'common_int.h'
      include 'common_mom.h'
      include 'common_masses.h'
      include 'common_phegas_mom.h'
      include 'common_norm.h'
      include 'common_strf.h'
      include 'common_cuts.h'
      include 'common_qcdrun_kaleu.h'
!
      dimension p(2,4),smin(-2:17,-2:17)
!      
      parmas(0)=0
      parwid(0)=-1
!
      gevnb=dnou(38937966)/100
      print*,'gevnb:',gevnb
      read*,e
      print*,'energy:',e
      read*,istruc
      print*,'istruc:',istruc
!
      xp1=1
      xp2=1
      he=e
!
! Initialize initial-state momenta for Helac
      p(1,4)=(e**2+parmas(ifl(1))**2-parmas(ifl(2))**2)/e/2
      p(1,3)=p(1,4)*dsqrt( 1-(parmas(ifl(1))/p(1,4))**2 )
      p(1,2)=0
      p(1,1)=0
      p(2,4)=(e**2+parmas(ifl(2))**2-parmas(ifl(1))**2)/e/2
      p(2,3)=-p(2,4)*dsqrt( 1-(parmas(ifl(2))/p(2,4))**2 )
      p(2,2)=0
      p(2,1)=0
      do i=1,2
        pt=dsqrt(p(i,1)**2+p(i,2)**2)
        pq=dsqrt(pt**2+p(i,3)**2)
        zq(i,1)= dcmplx(p(i,4)+p(i,3),p(i,3))
        zq(i,2)= dcmplx(p(i,4)-p(i,3),pt)
        zq(i,3)= dcmplx(p(i,1), p(i,2))
        zq(i,4)= dcmplx(p(i,1),-p(i,2))
        zq(i,5)= dcmplx( parmas(ifl(i)) , pq )
        pmom(i,1:4)=p(i,1:4)
        pmom(i,5)=parmas(ifl(i))
      enddo  
!
! Constant factor in weight       
      call mypi(pifac)
      pifac = (2*pifac)**(4-3*(n-2))
      wrest=4*dsqrt( scalar_product(p(1,1:4),p(2,1:4))**2
     .               -parmas(ifl(1))**2*parmas(ifl(2))**2   )
      wrest=gevnb/wrest * pifac
!     
! Read cuts     
      call readcuts
!
! Derive kinematical limits on 2-particle invariants from cuts
      do j1=   3,n-1
      do j2=j1+1,n
        k1 = j1-2
        k2 = j2-2
        rm1 = parmas(ifl(j1))
        rm2 = parmas(ifl(j2))
        if(rm1.eq.0.and.rm2.eq.0)then
          e1 = ec(j1)
          e2 = ec(j2)
          p1 = 0
          p2 = 0
          if (e1.gt.0) p1 = e1*dsqrt( 1-(rm1/e1)**2 )
          if (e2.gt.0) p2 = e2*dsqrt( 1-(rm2/e2)**2 )
          xx = rm1**2 + rm2**2 + 2*e1*e2 - 2*p1*p2*cc(j1,j2)
          smin(k1,k2) = max( xx ,gmas(j1,j2)**2 ,(rm1+rm2)**2 )
          smin(k2,k1) = smin(k1,k2)
        else
          smin(k1,k2) = max( gmas(j1,j2)**2 ,(rm1+rm2)**2 )
          smin(k2,k1) = smin(k1,k2)
        endif
      enddo
      enddo
      do j1=3,n
        k1  = j1-2
        rm1 = parmas(ifl(j1))
        e1  = ec(j1)
        ptsq = ptc(j1)**2
        k2 = -1
          rm2 = parmas(ifl(1))
!          ptsq = ( 1d0-c1(j1)**2 )*( e1**2-rm1**2 )
          smin(k1,k2) = rm1**2 + rm2**2 - ptsq/2
          smin(k2,k1) = smin(k1,k2)
        k2 = -2
          rm2 = parmas(ifl(2))
!          ptsq = ( 1d0-c2(j1)**2 )*( e1**2-rm1**2 )
          smin(k1,k2) = rm1**2 + rm2**2 - ptsq/2
          smin(k2,k1) = smin(k1,k2)
      enddo
!
! Initialize Kaleu
      call helac_kaleu_init( parmas,parwid,ihiggs,onlyqcd,withqcd 
     &                      ,ifl,n
     &                      ,e,istruc
     &                      ,nbatch,nstep
     &                      ,smin )
      end
!
      subroutine kaleu2helac_gnrt( discard )
!* ********************************************************************
!* ********************************************************************
      include 'declare.h'
!       
      include 'common_int.h'
      include 'common_mom.h'
      include 'common_masses.h'
      include 'common_phegas_mom.h'
      include 'common_strf.h'
!
      double precision pkaleu(0:3,-2:17)
      logical discard
!       
! Generate momenta
      call helac_kaleu_gnrt( discard ,xp1,xp2,pkaleu )
      if (discard) return
!
! Put initial-state momenta to Helac
      if (istruc.eq.1) then
        ehat = -pkaleu(0,-1)-pkaleu(0,-2)
        wjac = 1 ! wjac is already in the weight of Kaleu
        do ii=1,2
          jj = -ii
          pt = pkaleu(1,jj)**2 + pkaleu(2,jj)**2
          pq = dsqrt( pt + pkaleu(3,jj)**2 )
          pt = dsqrt( pt )
          zq(ii,1)= dcmplx(-pkaleu(0,jj) -pkaleu(3,jj) ,-pkaleu(3,jj) )
          zq(ii,2)= dcmplx(-pkaleu(0,jj) +pkaleu(3,jj) ,pt ) 
          zq(ii,3)= dcmplx(-pkaleu(1,jj),-pkaleu(2,jj) )
          zq(ii,4)= dcmplx(-pkaleu(1,jj),+pkaleu(2,jj) )
          zq(ii,5)= dcmplx( parmas(ifl(ii)) , pq )
          pmom(ii,4) = -pkaleu(0,jj)
          pmom(ii,1) = -pkaleu(1,jj)
          pmom(ii,2) = -pkaleu(2,jj)
          pmom(ii,3) = -pkaleu(3,jj)
          pmom(ii,5) = parmas(ifl(ii))
        enddo
      endif !(istruc.eq.1)
!
! Put final-state momenta to Helac
      do ii=3,n
        jj = ii-2
        pmom(ii,4) = pkaleu(0,jj)
        pmom(ii,1) = pkaleu(1,jj)
        pmom(ii,2) = pkaleu(2,jj)
        pmom(ii,3) = pkaleu(3,jj)
        pmom(ii,5) = parmas(ifl(ii))
        pt = pkaleu(1,jj)**2 + pkaleu(2,jj)**2
        pq = dsqrt( pt + pkaleu(3,jj)**2 )
        pt = dsqrt( pt )
        zq(ii,1)= dcmplx( pkaleu(0,jj) +pkaleu(3,jj) ,pkaleu(3,jj) )
        zq(ii,2)= dcmplx( pkaleu(0,jj) -pkaleu(3,jj) ,pt ) 
        zq(ii,3)= dcmplx( pkaleu(1,jj), pkaleu(2,jj) )
        zq(ii,4)= dcmplx( pkaleu(1,jj),-pkaleu(2,jj) )
        zq(ii,5)= dcmplx( parmas(ifl(ii)) , pq )
      enddo
      end
