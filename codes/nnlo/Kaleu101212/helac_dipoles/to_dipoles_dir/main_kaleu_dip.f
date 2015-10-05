       program main

       include 'declare.h'
       
       include 'common_int.h'
       include 'common_masses.h'
       include 'common_flags.h'
       include 'common_feyn.h'
       include 'common_new.h'
       include 'common_print.h'

       character*24, histo,error
 
       character*2 char(20)
       character*100 process
      
       call zmset(34)

       nunit1=6
     
       print*,'histo file'
       read*,histo,error
       print*,'histo file',histo,error
       print*,'repeat,ranhel'
       read*,repeat,iranhel
       print*,'repeat',repeat,'  ranhel',iranhel
       print*,'Number of particles'
       read*,n
       print*,n
       print*,'Flavour of particles'
       read*,(ifl(i),i=1,n)
       print*,(ifl(i),i=1,n)
       print*,'iflag,iunitary,ihiggs,iwidth'
       read*,iflag,iunitary,ihiggs,iwidth
       print*,iflag,iunitary,ihiggs,iwidth
       
       call physics

       read*,ikaleu,nbatch,nstep,nbatch_grid,nbatch_cut !Kaleu

       do i=1,n
       if(i.le.2)then
        io(i)=1
       else
        io(i)=-1
       endif
       enddo

       do i=1,n
        if(ifl(i).eq. 2)char(i)='E-'
        if(ifl(i).eq.-2)char(i)='E+'
        if(ifl(i).eq. 1)char(i)='Ne'
        if(ifl(i).eq.-1)char(i)='N1'
        if(ifl(i).eq. 3)char(i)='Uq'
        if(ifl(i).eq.-3)char(i)='Ua'
        if(ifl(i).eq. 4)char(i)='Dq'
        if(ifl(i).eq.-4)char(i)='Da'
        if(ifl(i).eq. 6)char(i)='M-'
        if(ifl(i).eq.-6)char(i)='M+'
        if(ifl(i).eq. 5)char(i)='Nm'
        if(ifl(i).eq.-5)char(i)='N2'
        if(ifl(i).eq. 7)char(i)='Cq'
        if(ifl(i).eq.-7)char(i)='Ca'
        if(ifl(i).eq. 8)char(i)='Sq'
        if(ifl(i).eq.-8)char(i)='Sa'
        if(ifl(i).eq. 10)char(i)='T-'
        if(ifl(i).eq.-10)char(i)='T+'
        if(ifl(i).eq. 9)char(i)='Nt'
        if(ifl(i).eq.-9)char(i)='N3'
        if(ifl(i).eq. 11)char(i)='Tq'
        if(ifl(i).eq.-11)char(i)='Ta'
        if(ifl(i).eq. 12)char(i)='Bq'
        if(ifl(i).eq.-12)char(i)='Ba'
        if(ifl(i).eq.31)char(i)='A0'
        if(ifl(i).eq.32)char(i)='Z0'
        if(ifl(i).eq.33)char(i)='W+'
        if(ifl(i).eq.34)char(i)='W-'
        if(ifl(i).eq.35)char(i)='G0'
        if(ifl(i).eq.41)char(i)='H0'
       enddo
       write(process,*)(char(i),i=1,n)
       print*,process(2:2*n+1)

       nunit2=32
       open(nunit2,file='kine_'//process(2:2*n+1)//'.out')
       nunit3=30
       open(nunit3,file='even_'//process(2:2*n+1)//'.out'
     . ,form='unformatted')

c for LHA
       open(200,file='sample'//process(2:2*n+1)//'.init')

       if(repeat.eq.0) then
          call mtime
          call helac_init
       endif
       if(repeat.eq.1) then
        open(21,file='tree_'//process(2:2*n+1)//'.in')
          call mtime
        call helac_init
          call mtime
        close(0)
        stop
       endif
       if(repeat.eq.2)then
        open(21,file='tree_'//process(2:2*n+1)//'.in')
        call getlist
        read(21,*)ng
        print*,'ng',ng
        do i1=1,ng
        do j1=1,8
        read(21,*)(is(i1,j1,k1),k1=1,8)
        enddo
        enddo
       endif
       
c some initialization for histograms       
       include 'nh.h'
        open(16,file=histo)
        open(17,file=error)
       write(16,*)ifl(1:n)
       do i=1,nh
        call histo3(i)
       enddo
       call histo3(199)
c ---------

c ---------
c THE MAIN INTEGRATION ROUTINE: see below       
c ---------
       call mtime               ! timing
       if (ikaleu.eq.0) then                               !Kaleu
       call integration(w0)
       else                                                !Kaleu
       call integration_kaleu(w0 ,ikaleu,nbatch,nstep      !Kaleu
     &                           ,nbatch_grid,nbatch_cut ) !Kaleu
       endif                                               !Kaleu
       call mtime ! timing
       
c ---------
c output histograms
       do i=1,nh
        call histo2(i,16,0,w0)
       enddo
c ---------
       stop
       end


       subroutine integration(w0)
 
  
!ANDRE       use, intrinsic :: ieee_arithmetic, only :
!ANDRE     .      ieee_set_flag, ieee_get_flag, ieee_usual,
!ANDRE     .      ieee_is_finite


       include 'declare.h'
       include 'common_int.h'
       include 'common_mom.h'
       include 'common_masses.h'
       include 'common_feyn.h'
       include 'common_phegas_mom.h'
       include 'common_new.h'
       include 'common_print.h'
       include 'common_norm.h'
       include 'common_strf.h'
       include 'common_warn.h'
       include 'common_unweight.h'
       include 'common_lha.h'
       include 'common_unw.h'
       include 'common_psp.h'
       include 'common_qcdrun.h'
       include 'common_debug.h'
       include 'common_onep.h'

       logical flag_value(3)

       common/myparni/xx,r2

       integer ipdt(-100:100)

       dimension alpha(ng+1),beta(0:ng+1),wei(ng+1),wpsp(ng+1)
     .           ,p(20,4),igraph(ng+1)
     .           ,icount0(0:ng+1),icount1(0:ng+1)
     .           ,wmax(10),alpha_m(ng+1),x(ng+1),y(ng+1)
 
       dimension ig_m(ng+1)
       
       real*8 hlower,hupper,hvaria,hweigh

       character*2 file
       logical onep
       logical lunwei
       data lunwei/.false./  
       data nunwei/1000/  
       logical lwmax,lwri
       data lwmax/.false./
       integer nwri,nwmax,nwarn
       data nwri/0/,nwmax/0/,nwarn/0/
       double precision umax,umax1
       data umax/0/

       save init
       data init/0/

       integer ind(ng+1)
       
       icount0(0:ng+1)=0
       icount1(0:ng+1)=0
      
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

          beta(0)=0
          igraph(1:ng+1)=0
          inonzero=0
          inonzeroplus=0
          inonzerominus=0
          read(*,*,end=1000)ncha
          print*,'option for kinematical channels:',ncha
          if(ncha.eq.0)then
             nch=ng+1
             do i=1,nch
                igraph(i)=mod(i,ng+1)
             enddo
          elseif(ncha.ge.1)then
             nch=1
             igraph(1)=0
          elseif(ncha.eq.-1)then
             read*,nch
             do i=1,nch
                igraph(i)=i
             enddo
          endif
          if(ncha.le.1)i_psp=1
          if(ncha.eq.2)i_psp=2

          print*,'Number of channels',nch
          print*,'Number of MC points'
          read*,nmc
          print*,' nmc  = ',nmc
          alpha(1:nch)=dnou(1)/nch
          alpha_m(1:nch)=alpha(1:nch)
          do i=1,nch
             beta(i)=beta(i-1)+alpha(i)
          enddo
          read*,nopti,nopt_step,optf,maxopt,noptlim,iopt
          nopt=100
          number_opt=0
          wmax(1:10)=0
          wmemax=0
          wme=1
          if(nch.eq.1)then
             iopt=0
             maxopt=-1
          endif
          print*,'nopt,nopt_step,optf,maxopt,iopt'
          print*, nopt,nopt_step,optf,maxopt,iopt
          print*,'number of channels=',nch
          imc=0
          if(iopt.eq.0)imc=0
          init=1
       endif
c --------------------------------------------------------------
c         END  OF INITIALIZATION
c---------------------------------------------------------------
 4     continue
       print*,'NUMBER OF CHANNELS',nch
       ipass=0
       iev=0
       do while(iev.le.nmc)  ![iev
        iev=iev+1
          w=0
          wme=0
          wei(1:nch)=0
          call igen(nch,beta(1:nch),ig)
          icount0(igraph(ig))=icount0(igraph(ig))+1
          idir=0
          if(iev.eq.1)then 
           ig=nch
           print*,nch,ig,igraph(ig)
          endif
          call phegas(wpsp(ig),igraph(ig),idir,iweight)

          if(iweight.eq.0) goto 2
          if(iweight.eq.-1) goto 1
       
          if(wpsp(ig).le.0)then
             goto 1
          endif

          if (iev.gt.1) then
!ANDRE             call ieee_set_flag(ieee_usual,.false.)
             call checkpoint(icut)
!ANDRE             call ieee_get_flag(ieee_usual,flag_value)
!ANDRE             if(any(flag_value)) icut=0
             if(icut.eq.0) goto 1
          endif

          wpsp_t=0
          do i=1,nch
             if(i.ne.ig)then
                idir=1
                call phegas(wpsp(i),igraph(i),idir,iweight)
                if(iweight.eq.0)goto 2
                if(iweight.eq.-1)goto 1
                if(wpsp(i).gt.0)then
                   wpsp_t=wpsp_t+alpha(i)/wpsp(i)
                else 
                   print*,'wpsp2=',wpsp(i),i,wpsp(ig),ig,nch,igraph(ig)
                endif
             else
                wpsp_t=wpsp_t+alpha(i)/wpsp(i)
             endif
          enddo
          inonzero=inonzero+1
          icount1(igraph(ig))=icount1(igraph(ig))+1
          wpsp_t=wpsp_t**(-1)
          if(irun.eq.1)call physics

          wsf=1
          if(istruc.eq.1)call strf(wsf)

!ANDRE          if(.not.ieee_is_finite(wrest*wsf*wpsp_t)) goto 1

          if (iopt.eq.0.or.number_opt.gt.maxopt) then
             call dipoles(3,wme,wrest*wsf*wpsp_t)
          else
             call dipoles(2,wme,wrest*wsf*wpsp_t)
          endif

       if(iev.eq.1)goto 1

       w=wpsp_t*wme*wrest
       if (.not.w.le.0d0.and..not.w.ge.0d0) then !DEBUG
         write(6,*) 'wpsp_t,wme:',wpsp_t,wme !DEBUG
         call printmomenta
         stop !DEBUG
       endif !DEBUG
       w=w*wsf

       ipass=ipass+1

       wr=1
       include 'adapt1.h'
       include 'ktreweight.h'

       if(.not.lwmax)then
       include 'myparni.h'
       endif

       w1=w

       if(number_opt.le.maxopt)then
          wei(1:nch)=0
          do i=1,nch
             if(wpsp(i).gt.0) wei(i)=w*w*wpsp_t/wpsp(i)
          enddo
       endif

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
                   ii=igraph(ig)
                   write(nunit2,*)iev,igraph(ig),alpha(ig)
                   write(nunit2,*)w,wme,wpsp_t
                   if(ii.gt.0)then
                      k=1
                      do while(is(ii,k,1).ne.0)
                         write(nunit2,*)is(ii,k,1:8)
                         k=k+1
                      enddo
                   endif
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
             if(.not.lwmax.and.
     .            (number_opt.gt.maxopt.or.iopt.eq.0))then
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
          if(lwri)call bookin(4,w)
          if(number_opt.le.maxopt)then
             do i=1,nch
                call bookin(i+4,wei(i))
             enddo
          endif 
 2        continue   
          call geti(1,0,w0)
          call geti(2,0,w00)
          if(mod(iev,1000).eq.0)then
             call errest(1,x1,y1,0)
             if(x1.gt.0)
     .       print'(a6,2d15.6,3i10)','sigma=',
     .            x1,sqrt(y1)/x1,inonzeroplus,int(w0),iev
             call errest(2,x2,y2,0)
             if(x2.gt.0)
     .       print'(a6,2d15.6,3i10)','sigma=',
     .            x2,sqrt(y2)/x2,inonzerominus,int(w00),iev
             call errest(4,xx1,yy1,0)
             if(yy1.gt.0)print'(a7,2d15.6)','<w=/=0>',
     .            xx1,sqrt(yy1)/dabs(xx1)
             if(nwri.gt.0)print*,nwri
             write(17,*)iev,x1,sqrt(y1)
             print*,'----------------------------'
          endif
          if(iev.eq.nopt)then
             if(ipass.eq.0)nopt=nopt*2
          endif
c optimization   [
          if(iev.eq.nopt.and.iopt.eq.1
     .         .and.number_opt.le.maxopt.and.nch.gt.1)then
             sum=0
             do i=1,nch 
                call errest(i+4,x(i),y(i),0)
                sum=sum+alpha(i)*x(i)
             enddo
             nch0=nch
           call optimize(nch,igraph(1:nch),alpha,nm,ig_m,alpha_m,x,iend)
             if(nch.eq.1)iopt=0
             if(iend.eq.2)number_opt=maxopt
             if(number_opt.ge.1)then
c
       include 'adapt3.h'
c
                data iunit/80/
                write(nunit1,*)'UNIT=',iunit
                write(file,'(i2)') iunit
                open(iunit,file='data.'//file)
                write(iunit,'(i4)')nch
                do i=1,nch
                   write(iunit,'(3(e10.4,5x),i8)')
     .                  alpha(i),x(i),dsqrt(y(i)),igraph(i)
                enddo
                close(iunit)
                iunit=iunit+1
             endif
             if(number_opt.eq.0)then
                nopt=nopti
             else
                nopt=nopt+nopt_step
             endif
             write(nunit1,*)'----------------------------'
             write(nunit1,*)'iev,nopt,number_opt,maxopt,nch,nch0'
     .            , iev,nopt,number_opt,maxopt,nch,nch0
             if(nopt.gt.noptlim)number_opt=maxopt
             nopt_step=nopt_step*optf
             write(nunit1,*)'nopt_step',nopt_step
             number_opt=number_opt+1
             if(number_opt.gt.maxopt)then
                nch=nm
                igraph(1:nch)=ig_m(1:nch)
                print*,'# of graphs , # of channels',ng,nch
                do i=1,nch
                   alpha(i)=alpha_m(i)
                   if(10*alpha(i).gt.1)then
                      k=1
                      ii=igraph(i)
                      write(nunit1,*)'the graph',ii,alpha(i)
                      if(ii.eq.0)cycle
                      do while(is(ii,k,1).ne.0)
                         write(nunit1,*)is(ii,k,1:8)
                         k=k+1
                      enddo
                   endif
                enddo
                open(iunit,file='data.final')
                write(iunit,'(i4)')nch
                do i=1,nch
                   write(iunit,'(e10.4,i8)')alpha(i),igraph(i)
                enddo
                close(iunit)
                iunit=iunit+1
             endif
c optimization   ]
             do i=1,nch
                beta(i)=beta(i-1)+alpha(i)
             enddo
             if(number_opt.gt.maxopt.and.nch.gt.555)then
              call SOR(alpha,nch,ind,1)
              sum=0
              i=0
              do while(sum.lt.0.90d0)
               i=i+1
               sum=sum+alpha(i)
              enddo
              nch=i
              if(nch.gt.555) nch=555
               ig_m(1:nch)=igraph(1:nch)
               do i=1,nch
                igraph(i)=ig_m(ind(i))
               enddo
              sum=0
              do i=1,nch
               sum=sum+alpha(i)
              enddo

              print*,'ALPHA(I)',sum,nch
              print'(d12.4)',alpha(1:nch)
              print*,ind(1:nch)
              do i=1,nch
              print*,igraph(i)
              print'(8I6)',(is(igraph(i),j,1:8),j=1,n-2)
              enddo

              alpha(1:nch)=alpha(1:nch)/sum
              wmax(1:10)=0
              imc=0
              do i=1,nch
                beta(i)=beta(i-1)+alpha(i)
              enddo
             endif

             if(number_opt.eq.1.or.number_opt.gt.maxopt)then
                write(nunit1,*)'clear ',number_opt,' opt'
                do l=1,nch0
                   call clear(l+2)
                enddo
                inonzero=0
                inonzeroplus=0
                inonzerominus=0
                call clear(1)
                call clear(2)
                goto 4
             endif

 3           continue       
          endif
       enddo ![iev

 5      continue
       call geti(1,0,w0)
       print*,'out of ',nmc,'  ',int(w0),' points have been used'
       print*,'and  ',inonzero,' points resulted to =/= 0 weight'
       print*,'whereas  ',int(w0)-inonzero,' points to 0 weight'
       call errest(1,x1,y1,1)
       call errest(2,x2,y2,1)
       print*,'total XS',x1-x2,sqrt(y1+y2)

       call geti(4,0,ss0)
       print*,'lwri: points have used',ss0

       if(ncha.ne.0.or.number_opt.eq.0) read*,alimit,nlimit
 
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
      
     
       subroutine optimize(n,igraph,a,nm,ig_m,am,x,iend)
       include 'declare.h'

       include 'common_feyn.h'
       include 'common_debug.h'
       dimension a(n),am(n),x(n),igraph(n),ig_m(*)
       save
       data init/0/,ex/.25d0/,alimit/0/
       data icount/0/

       if(init.eq.0)then
                     if(n.le.100)nlim=10
        if(n.gt.100.and.n.le.200)nlim=20
        if(n.gt.200.and.n.le.500)nlim=40
                     if(n.gt.500)nlim=60
        d=0
        xmax=0
        iend=0
        read*,alimit,nlim
        print*,'LIMIT=',alimit,nlim
        print*,'DUMPFAC=',ex
        goto13
        return
       endif

       print*,'LIMIT=',alimit

       a1=0
       do i=1,n
	a1=a1+a(i)*x(i)**ex
       enddo

       xmax=0
       do i=1,n
        if(i.eq.1)then
         xmax=abs(a1-x(1)**ex)
        else
         xmax=max(xmax,abs(a1-x(i)**ex))
        endif
       enddo
        
       if(init.eq.1)then
        d=xmax
        init=2
       endif
       
       print*,'This opt. max=',xmax,n

       if(xmax.le.d)then
        ig_m(1:n)=igraph(1:n)
        nm=n
        am(1:n)=a(1:n)
        d=xmax
       else
       endif

       print*,'All opt.    d=',d
       icount=icount+1

 13    continue
       do i=1,n
        a(i)=a(i)*x(i)**ex
       enddo
       
       n1=n
 11    continue
       call normalize(a(1:n1),n1)

       if(alimit.eq.0.and.init.gt.0)return

       if(init.gt.0)then
        j=n1
        do k=1,n1
         if(a(k).lt.alimit)j=j-1
        enddo
        if(j.le.nlim)return
       endif
       k=0
       j=0
       istop=0
       do while(istop.eq.0)
        k=k+1
        j=j+1
        if(j.eq.n1)istop=1
	if(init.eq.0)then
        eps1=dnou(1)/dnou(10)**11
	do l=1,n1
         ratio=abs(a(k)-a(l))/(a(k)+a(l))
         if(k.ne.l.and.ratio.lt.eps1.and.a(k).gt.0)then
          if(idebug.eq.1)then
          data ix/0/
          ix=ix+1
          print*,igraph(k),a(k)
          print*,igraph(l),a(l)
          print*,is(igraph(k),1,1:8)
          print*,is(igraph(k),2,1:8)
          print*,is(igraph(l),1,1:8)
          print*,is(igraph(l),2,1:8)
          print*,ix,'----------------------------------'
          endif
          igraph(k:n1-1)=igraph(k+1:n1)
          igraph(n1)=0
          a(k:n1-1)=a(k+1:n1)
          a(n1)=0
          k=k-1
	  goto 12
         endif
	enddo
	endif
        if(a(k).lt.alimit.and.init.gt.0)then
         igraph(k:n1-1)=igraph(k+1:n1)
         igraph(n1)=0
         a(k:n1-1)=a(k+1:n1)
         a(n1)=0
         k=k-1
        endif

 12     continue
        enddo

       if(k.lt.n1)then
        n1=k
        call normalize(a(1:n1),n1)
       endif
       
       if(n1.le.n.and.init.eq.0)then
	d=0
	init=1
        do k=1,n1
         a(k)=dnou(1)/n1
        enddo
       endif
       n=n1
       
       return
       end

       subroutine igen(n,b,ig)
       include 'declare.h'

       dimension b(1:n)

       r=rnmy(0)
       do i=1,n
       ig=i
       if(r-b(i).lt.0)return
       enddo 

       end
      
       subroutine normalize(a,n)
       include 'declare.h'
       dimension a(n)
       
       atot=0
       do i=1,n
        atot=atot+a(i)
       enddo
       do i=1,n
        a(i)=a(i)/atot
       enddo
       
       return
       end

       subroutine put_unwei(x)
       include 'declare.h'
       parameter (n=500)
       dimension y(n)
       save y,init

c use weighted events
c      x=0
c      return

c use unweighted events
       if(init.eq.0)then
        y(1:n)=0
        init=1
       endif
 
       do i=1,n
        if(x.gt.y(i))then
         y(i+1:n)=y(i:n-1) 
         y(i)=x
         return
        endif
       enddo
       return
       
       entry get_unwei(t)
 
c use weighted events
c      t=0
c      return

c use unweighted events
c      sum=0
c      do i=1,n
c       sum=sum+y(i)
c      enddo
c      t=sum/dble(n)
c      t=y(5)
       
       call histo3(200)
       a=y(n)-1d-11
       b=y(1)+1d-11
       nb=100
       w=1.
       do i=1,n
        call histo1(200,nb,a,b,y(i),w)
       enddo
       open(198,file='tmp16')
       call histo2(200,198,0,w)
       rewind(198)
       sum=0
       do i=1,nb
        read(198,*)hx,hy,he,hn
        sum=sum+hn
        if(sum.gt.0.9*n)then
         t=hx
         close(198)
         return
        endif
       enddo

       return
       end


      SUBROUTINE SOR(A,N,K,IOPT)
C-----------------------------------------------------------------------
C     Sort A(N) into ascending order
C     IOPT = 1 : return sorted A and index array K
C     IOPT = 2 : return index array K only
C-----------------------------------------------------------------------
c     DOUBLE PRECISION A(N),B(5000)
c     INTEGER N,I,J,IOPT,K(N),K1(5000),IL(5000),IR(5000)
      DOUBLE PRECISION A(N),B(N)
      INTEGER N,I,J,IOPT,K(N),K1(N),IL(N),IR(N)
c     IF (N.GT.5000) then
c       write(*,*) 'Too many entries to sort in srt, stop'
c       stop
c     endif
      if(n.le.0) return
      IL(1)=0
      IR(1)=0
      DO 10 I=2,N
      IL(I)=0
      IR(I)=0
      J=1
   2  IF(A(I).GT.A(J)) GOTO 5
   3  IF(IL(J).EQ.0) GOTO 4
      J=IL(J)
      GOTO 2
   4  IR(I)=-J
      IL(J)=I
      GOTO 10
   5  IF(IR(J).LE.0) GOTO 6
      J=IR(J)
      GOTO 2
   6  IR(I)=IR(J)
      IR(J)=I
  10  CONTINUE
      I=1
      J=1
      GOTO 8
  20  J=IL(J)
   8  IF(IL(J).GT.0) GOTO 20
   9  K(I)=J
      B(I)=A(J)
      I=I+1
      IF(IR(J)) 12,30,13
  13  J=IR(J)
      GOTO 8
  12  J=-IR(J)
      GOTO 9
  30  CONTINUE
      DO I=1,N
      K1(I)=K(N+1-I)
      ENDDO
      DO I=1,N
      K(I)=K1(I)
      ENDDO
      IF(IOPT.EQ.2) RETURN
      DO 31 I=1,N
  31  A(I)=B(N+1-I)
 999  END
