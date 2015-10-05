       subroutine integration_kaleu( w0 ,ikaleu,nbatch,nstep
     &                                  ,nbatch_grid,nbatch_cut )

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
       include 'common_warn.h'
       include 'common_unweight.h'
!       include 'common_lha.h'
       include 'common_unw.h'
       include 'common_psp.h'
       include 'common_qcdrun.h'
       include 'common_debug.h'
       include 'common_onep.h'

       integer ipdt(-100:100)

       dimension wmax(10)
 
       logical onep ,discard,initphase

       logical lunwei
       data lunwei/.false./  
       data nunwei/1000/  
       logical lwmax,lwri
       data lwmax/.false./
       integer nwri,nwmax,nwarn
       data nwri/0/,nwmax/0/,nwarn/0/
       double precision umax,umax1
       data umax/0/

       save init,iplot,initphase
       data init/0/,iplot/0/,initphase/.true./

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
      print*,'UNWEIGHTING IS', lunwei

          inonzero=0
          inonzeroplus=0
          inonzerominus=0
          read(*,*,end=1000)ncha
          print*,'Number of MC points'
          read*,nmc
          print*,' nmc  = ',nmc
          read*,idummy1,idummy2,xdummy3,idummy4,idummy5,idummy6
          print*,'nbatch,nstep,nbatch_grid,nbatch_cut:'
     &          ,nbatch,nstep,nbatch_grid,nbatch_cut
          nopt=100
          wmax(1:10)=0
          wmemax=0
          wme=1
          imc=0
          init=1
       endif

       call kaleu2helac_init( ikaleu,nbatch,nstep      ! This routine is
     &                       ,nbatch_grid,nbatch_cut ) ! defined below

       noptim = nbatch*nstep ! number of accepted events to be thrown away
       nplot  = 100000 ! number of accepted events before plotting grids

c --------------------------------------------------------------
c         END  OF INITIALIZATION
c---------------------------------------------------------------
 4     continue
       ipass=0
       iev=0
       do while(iev.le.nmc)  ![iev
        iev=iev+1
          w=0
          wme=0
          call kaleu2helac_gnrt( discard ) ! This routine is defined below
          if (discard) goto 1

          if (iev.gt.1) then
             call cuts(icut)
             if(icut.eq.0) goto 1
          endif

          call helac_kaleu_wght( wpsp_t )
          if (wpsp_t.eq.0d0) goto 1

c         call cpu_time(t2)
c         write(122,*)t2-t1
c         print*,'PS',t2-t1
          inonzero=inonzero+1
          if(irun.eq.1)call physics

          wsf=1
          if(istruc.eq.1)call strf(wsf)

c         call cpu_time(t1)
          call dipoles_I(wme,wrest*wsf*wpsp_t)
c         call cpu_time(t2)
c         write(123,*)t2-t1
c         print*,'ME',t2-t1
c        if(iev.gt.1000)stop
       if(iev.eq.1)goto 1

       w=wpsp_t*wme*wrest
       w=w*wsf
       ipass=ipass+1

       if (.not.w.le.0d0.and..not.w.ge.0d0) then     !DEBUG
         write(6,*) 'wpsp_t,wme,wsf:',wpsp_t,wme,wsf !DEBUG
         call printmomenta                           !DEBUG
         stop                                        !DEBUG
       endif                                         !DEBUG

       call helac_kaleu_collect( dabs(w) )

c      print*,'iev',iev,wme,wpsp_t

       wr=1
       include 'adapt1.h'
       include 'ktreweight.h'

       w1=w

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

c WRITE OUT EVENTS
c START
c --------------
          if(lunwei)then
c --------------
             if(.not.lwmax.and..not.initphase)then
                if(umax.lt.w1)then
                   umax=w1
                endif
                nwmax=nwmax+1
                umax1=umax
                call put_unwei(w1)
             endif
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
                call cmstolab
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

             if(nwri.eq.nwri_tot) go to 5

c --------------
          endif
c --------------
c END
 1        continue  

       include 'adapt2.h'

          if(w.gt.0)then
             inonzeroplus=inonzeroplus+1
           call bookin(1,w)
           call bookin(2,0d0)
          elseif(w.lt.0)then 
             inonzerominus=inonzerominus+1
           wmi=-w
           call bookin(2,wmi)
           call bookin(1,0d0)
          else
           call bookin(1,w)
           call bookin(2,w)
          endif
          call bookin(3,wme)
          call bookin(4,w)
          if(lwri)call bookin(4,w)
 2        continue   
          call geti(1,0,w0)
          call geti(2,0,w00)
          if(mod(iev,10000).eq.0)then
             call errest(1,x1,y1,0)
             if(x1.gt.0)
     .       print'(a6,2d15.6,3i10)','sigma=',
     .            x1,sqrt(y1)/x1,inonzeroplus,int(w0),iev
             call errest(2,x2,y2,0)
             if(x2.gt.0)
     .       print'(a6,2d15.6,3i10)','sigma=',
     .            x2,sqrt(y2)/x2,inonzerominus,int(w00),iev
             if(abs(x1-x2).gt.0)
     .       print'(a6,2d15.6,3i10)','sigma=',
     .            x1-x2,sqrt(y1+y2)/abs(x1-x2),
     .            inonzeroplus+inonzerominus,int(w00),iev
             call errest(4,xx1,yy1,0)
             if(abs(xx1).gt.0)print'(a7,2d15.6)','sigma=',
     .            xx1,sqrt(yy1)
c            if(yy1.gt.0)print'(a7,2d15.6)','<w=/=0>',
c    .            xx1,sqrt(yy1)/dabs(xx1)
             if(nwri.gt.0)print*,nwri
             write(17,*)iev,xx1,sqrt(yy1)
             print*,'----------------------------'
          endif

         if (ipass.ne.iplot.and.mod(ipass,nplot).eq.0) then
           iplot = ipass
           print*,''
           print*,'Plotting grids'
           print*,''
           call helac_kaleu_plotgrids( 21 ) ! unit
         endif

         if (initphase.and.ipass.eq.noptim) then
           initphase = .false.
           print*,''
           print*,'Throwing away estimates'
     &           ,', keeping optimized distributions'
     &           ,', and restarting MC'
           print*,''
           inonzero = 0
           inonzeroplus = 0
           inonzerominus = 0
           call clear(1)
           call clear(2)
           call clear(3)
           call clear(4)
         endif

       enddo ![iev

 5      continue
       call geti(1,0,w0)
       print*,'total number of acc. and gen. events:',ipass,iev
       print*,'number of events used for estimates :',inonzero,int(w0)
       call errest(1,x1,y1,1)
       call errest(2,x2,y2,1)
       print*,'total XS',x1-x2,sqrt(y1+y2)

       call geti(4,0,ss0)
       print*,'lwri: points have used',ss0

       read*,alimit,nlimit
 
        read*,IDBMUP1,IDBMUP2,EBMUP1,EBMUP2,IDWTUP,NPRUP
        print*,IDBMUP1,IDBMUP2,EBMUP1,EBMUP2,IDWTUP,NPRUP
       WRITE(200,5100) IDBMUP1,IDBMUP2,EBMUP1,EBMUP2,
     &      iPDFGUP1,iPDFGUP2,iPDFSUP1,iPDFSUP2,IDWTUP,NPRUP
       WRITE(200,5200) x1*10d0**3,dsqrt(y1)*10d0**3,1.d0,81

       if(x1.gt.0)print*,'% error:',dsqrt(y1)/x1*100
       write(16,*)x1,sqrt(y1)
       if(x2.gt.0)print*,'% error:',dsqrt(y2)/x2*100
       if(wmax(1).gt.0)
     . print*,'<w>/w_max,w_max',x1/wmax(1),wmax(1)
       if(wmax(10).gt.0)
     . print*,'<w>/w_max,w_max',x1/wmax(10),wmax(10)
       print'(10(E14.6,2X))',wmax(1:10)
       call errest(3,x1,y1,0)
       if(wmemax.gt.0)
     . print*,'<me>/memax,memax',x1/wmemax,wmemax
     
       do j=1,25
        if(iwarning(j).gt.0)print*,'iwarning(',j,') = ',iwarning(j)
        if(iwonders(j).gt.0)print*,'iwonders(',j,') = ',iwonders(j)
       enddo
       print*,'number of w=1 events',nwri
       print*,'number of w>1 events',nwarn
       print*,'maximum weight used for un:vs',umax,umax1

 5100 FORMAT(1P,2I8,2E14.6,6I6)
 5200 FORMAT(1P,3E14.6,I6)

      call histo2(199,199,1,1.)

1000   return
       end
