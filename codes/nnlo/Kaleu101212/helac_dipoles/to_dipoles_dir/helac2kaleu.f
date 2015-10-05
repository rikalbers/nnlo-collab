      subroutine kaleu2helac_init( ikaleu,nbatch,nstep
     &                            ,nbatch_grid,nbatch_cut )
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
!      include 'common_cuts.h'
      include 'common_qcdrun.h'
      include 'common_dipoles.h'
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
! Read cuts. The following routine is defined below 
      call helac2kaleu_cuts( smin )
      small = 1d-6*e**2 ! minimum value for border below/above smin
      do i=-2,n-3
      if (i.eq.0) cycle
      do j=i+1,n-2
      if (j.le.0) cycle
        if (dabs(smin(i,j)).lt.small) then
          if (i.lt.0) smin(i,j) = -small
          if (i.gt.0) smin(i,j) =  small
          smin(j,i) = smin(i,j)
        endif
        write(*,'(a32,i3,a1,i3,a3,f14.4)')                 !DEBUG
     &    ' MESSAGE from helac2kaleu: smin(',i,',',j,') =' !DEBUG
     &                               ,smin(i,j)            !DEBUG
      enddo
      enddo
!
!          nbatch =             ! number of accepted ps-points before optimization step
!           nstep =             ! number of optimization steps
             thrs = 1d-3         ! threshold for channel removal after optimization
          noffset = nbatch*nstep ! number of accepted events before start grids
!     nbatch_grid =             ! for optimization grids
!     nbatch_cut  =             ! for optimization of choice below/above smin
!
! Initialize Kaleu
      call helac_kaleu_init( parmas,parwid,ihiggs,onlyqcd,withqcd 
     &                      ,ifl,n
     &                      ,e,istruc
     &                      ,nbatch,nstep,thrs
     &                      ,smin
     &                      ,noffset,nbatch_grid,nbatch_cut
     &                    ,ijkdip,ndip,mxdip
!     &                    ,alphaMaxFF,alphaMaxFI,alphaMaxIF,alphaMaxII )
     &                    ,1d0 ,1d0 ,1d0 ,1d0 )
      end
!
      subroutine helac2kaleu_cuts( smin )
!* ********************************************************************
!* * Derive kinematical limits on 2-particle invariants from cuts
!* ********************************************************************
      include 'declare.h'
!       
      include 'common_int.h'
      include 'common_masses.h'
      include 'common_cuts.h'
!
      dimension p(2,4),smin(-2:17,-2:17)
!      
! Read LO cuts. The following routine is defined in readcuts_kaleu_dip.f 
      call readcuts_kaleu
!
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
          smin(k1,k2) = rm1**2 + rm2**2 - ptsq/2
          smin(k2,k1) = smin(k1,k2)
        k2 = -2
          rm2 = parmas(ifl(2))
          smin(k1,k2) = rm1**2 + rm2**2 - ptsq/2
          smin(k2,k1) = smin(k1,k2)
      enddo
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
 
      subroutine printmomenta
!* ********************************************************************
!* ********************************************************************
      include 'declare.h'
      include 'common_int.h'
      include 'common_phegas_mom.h'
      do ii=1,n
        write(6,'(i2,4d24.16)')
     &  ii,pmom(ii,4),pmom(ii,1),pmom(ii,2),pmom(ii,3)
      enddo
      end
