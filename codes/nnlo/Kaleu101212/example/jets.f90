module jets_module
  use particles
  use cuts_module
  private
  public :: jets_type,jets_init,jets_pass
!
  real(kind(1d0)) ,parameter :: twopi=6.2831853071795864769252867665590d0
!
  type :: jets_type
    private
    integer :: ngluon=0,nparton=0,nfinst=0
    integer ,dimension(20) :: gluon,parton
    real(kind(1d0)) :: dR=0.4d0
    type(cuts_type) :: cuts
  end type
!
contains
!
!
  subroutine jets_init( jets ,cuts ,process,nfinst ,masses ,smin )
!***********************************************************************
!* Allows only gluons to become soft and/or collinear.
!* The jet algorithm finds the pair of partons with the smallest
!* value of the dmeasure in the kT-algorithm, and adds the momenta.
!* The thus obtained set of momenta (so with one less) is then subject
!* to LO cuts.
!***********************************************************************
  implicit none
  type(jets_type) ,intent(inout) :: jets
  type(cuts_type) ,intent(inout) :: cuts
  integer         ,intent(in)  :: process(-2:17),nfinst
  real(kind(1d0)) ,intent(in)  :: masses(nfirst:nlast)
  real(kind(1d0)) ,intent(out) :: smin(-2:17,-2:17)
  integer         :: ii,tmpproc(-2:17)
  type(cuts_type) :: tmpcuts
!
  jets%nfinst = nfinst
!
  jets%ngluon  = 0
  jets%nparton = 0
  do ii=-2,nfinst
    if (ii.eq.0) cycle
    if (is_qcd(abs(process(ii)))) then
      jets%nparton = jets%nparton+1
      jets%parton(jets%nparton) = ii
    endif
    if (ii.gt.0.and.abs(process(ii)).eq.gluon) then
      jets%ngluon = jets%ngluon+1
      jets%gluon(jets%ngluon) = ii
    endif
  enddo
!
  if (jets%ngluon.eq.0) then
    write(*,*) 'ERROR in jets_init: no final-state gluons present.'
    stop
  endif
!
! Cuts going to be applied after one recombination
  ii = jets%gluon(1)
  tmpproc(-2:ii-1    ) = process(  -2:ii-1)
  tmpproc(ii:nfinst-1) = process(ii+1:nfinst)
  call cuts_init( cuts ,tmpproc,nfinst-1 ,masses ,smin )
!
! Determine smin for LO cuts promoted to NLO phase space
  call cuts_init( jets%cuts ,process,nfinst ,masses ,smin )
!
  end subroutine
!
!
  function jets_pass( jets,cuts ,pkaleu ) result(value)
  implicit none
  type(jets_type) ,intent(in) :: jets
  type(cuts_type) ,intent(in) :: cuts
  real(kind(1d0)) ,intent(in) :: pkaleu(0:3,-2:17)
  logical :: value
  integer :: ii,jj,imin,jmin,l1,l2
  real(kind(1d0)) :: pt1,dmeasure,pt2,y1,y2,dphi,dmin,ptmp(0:3,-2:17),sgn
!
!  value = cuts_pass( jets%cuts ,pkaleu ); return !DEBUG
  imin = 0
  jmin = 0
  dmin = 9d99
  do l1=1,jets%ngluon
    jj = jets%gluon(l1)
    do l2=1,jets%nparton
      ii = jets%parton(l2)
      if (jj.eq.ii) cycle
!...  final-initial    
      if     (ii.lt.0) then
         pt1=pkaleu(1,jj)**2+pkaleu(2,jj)**2
         dmeasure=pt1
         sgn = sign(1d0,-pkaleu(3,ii)*pkaleu(3,jj))
!...  final-final   
      elseif (ii.gt.0) then
         pt1=pkaleu(1,jj)**2+pkaleu(2,jj)**2
         pt2=pkaleu(1,ii)**2+pkaleu(2,ii)**2
         y1=0.5d0*log((pkaleu(0,jj)+pkaleu(3,jj))/(pkaleu(0,jj)-pkaleu(3,jj)))
         y2=0.5d0*log((pkaleu(0,ii)+pkaleu(3,ii))/(pkaleu(0,ii)-pkaleu(3,ii)))
         dphi= ph4(pkaleu(1:3,jj)) - ph4(pkaleu(1:3,ii))
         dphi=min(dabs(dphi),twopi-dabs(dphi))
         dmeasure=min(pt1,pt2)*((y1-y2)**2+dphi**2)/jets%dR**2
         sgn = 1d0
      endif
      if (dmeasure.lt.dmin.and.sgn.eq.1d0) then
        dmin = dmeasure
        imin = ii
        jmin = jj
      endif
    enddo
  enddo
!
  ptmp(0:3,-2:jets%nfinst) = pkaleu(0:3,-2:jets%nfinst)
  ptmp(0:3,imin) = ptmp(0:3,imin)+ptmp(0:3,jmin)
  ptmp(0:3,jmin:jets%nfinst-1) = pkaleu(0:3,jmin+1:jets%nfinst)
!
  value = cuts_pass( cuts ,ptmp )  
!
  end function
!
!
  function ph4(p) result(value)
  implicit none
  real(kind(1d0)) ,intent(in) :: p(3)
  real(kind(1d0)) :: value
  real(kind(1d0)) :: ptsq,s2,s
  ptsq = p(1)**2+p(2)**2
  if (ptsq.le.0) then
    if (p(3).ge.0) value=0
    if (p(3).lt.0) value=twopi/2
    return
  endif
  s2=p(1)**2/ptsq
  s=sqrt(s2)
  if(s.gt.1)then
    write(*,*) 'WARNING from jets_ph4: s=',s
    value=0
    if (p(1).lt.0) value=twopi/2
    return
  endif
  if (p(1).lt.0) s=-sqrt(s2)
  value = acos(s)
  if (p(2).lt.0) value=twopi-value
  end function
!
!
end module
