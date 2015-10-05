module rambo_module
  use kaleu_particles
!
  private
  public :: rambo_type,rambo_init,rambo_gnrt,rambo_wght
!
  integer ,parameter :: nmax=17
!
  type :: rambo_type
    private
    real(kind(1d0)) :: w=0d0
    real(kind(1d0)) :: m(-2:nmax)=0d0
    real(kind(1d0)) :: e=0d0
    integer         :: n=0
    integer         :: iw(5)=0
  end type
    
!
contains
!
!
  subroutine rambo_init( obj ,process,nfinst ,ecm,option ,masses )
  implicit none
  type(rambo_type) ,intent(inout) :: obj
  integer         ,intent(in) :: process(-2:17),nfinst,option
  real(kind(1d0)) ,intent(in) :: ecm,masses(nfirst:nlast)
  real(kind(1d0)) :: pp(4,100),xm(100)
  integer :: ii
  logical :: firstcall=.true.
  if (firstcall) then
    firstcall = .false.
    write(*,'(a72)') '########################################################################'
    write(*,'(a72)') '#                                                                      #'
    write(*,'(a72)') '#                         You are using RAMBO                          #'
    write(*,'(a72)') '#                                                                      #'
    write(*,'(a72)') '#                     for phase space generation                       #'
    write(*,'(a72)') '#                                                                      #'
    write(*,'(a72)') '# AUTHORS:  S.D. ELLIS,  R. KLEISS,  W.J. STIRLING                     #'
    write(*,'(a72)') '# THIS IS VERSION 1.0 -  WRITTEN BY R. KLEISS                          #'
    write(*,'(a72)') '# Wrapped by A. van Hameren                                            #'
    write(*,'(a72)') '#                                                                      #'
    write(*,'(a72)') '########################################################################'
  endif
  if (option.ne.0) then
    write(*,*) 'ERROR in rambo: only option=0 possible for now.'
    stop
  endif
  obj%n = nfinst
  obj%e = ecm
  do ii=-2,nfinst
    if (ii.eq.0) cycle
    obj%m(ii) = masses(abs(process(ii)))
  enddo
  end subroutine
!  
!
  subroutine rambo_gnrt( obj, pkaleu )
  implicit none
  type(rambo_type) ,intent(inout) :: obj
  real(kind(1d0)) ,intent(out) :: pkaleu(0:3,-2:17)
  real(kind(1d0)) :: stot,s1,s2,xm(100),pp(4,100)
  stot = obj%e**2
  s1 = obj%m(-1)**2
  s2 = obj%m(-2)**2
  pkaleu(0,-1) = -( stot + s1 - s2 ) / (2*obj%e)
  pkaleu(3,-1) = -dsqrt( pkaleu(0,-1)**2 - s1 )
  pkaleu(2,-1) =  0d0
  pkaleu(1,-1) =  0d0
  pkaleu(0,-2) = -( stot + s2 - s1 ) / (2*obj%e)
  pkaleu(3,-2) = -pkaleu(3,-1)
  pkaleu(2,-2) = -pkaleu(2,-1)
  pkaleu(1,-2) = -pkaleu(1,-1)
!
  xm(1:obj%n) = obj%m(1:obj%n)
  call rambo( obj%n ,obj%e ,xm ,pp ,obj%w ,obj%iw )
!
  pkaleu(  0,1:obj%n) = pp(  4,1:obj%n)
  pkaleu(1:3,1:obj%n) = pp(1:3,1:obj%n)
  end subroutine  
!
!
  subroutine rambo_wght( obj ,weight )
  implicit none
  type(rambo_type) ,intent(in) :: obj
  real(kind(1d0)) ,intent(out) :: weight
  weight = obj%w
  end subroutine
!
!
      SUBROUTINE RAMBO(N,ET,XM,P,WT,IWARN)
!------------------------------------------------------
!
!                       RAMBO
!
!    RA(NDOM)  M(OMENTA)  B(EAUTIFULLY)  O(RGANIZED)
!
!    A DEMOCRATIC MULTI-PARTICLE PHASE SPACE GENERATOR
!    AUTHORS:  S.D. ELLIS,  R. KLEISS,  W.J. STIRLING
!    THIS IS VERSION 1.0 -  WRITTEN BY R. KLEISS
!
!    N  = NUMBER OF PARTICLES (>1, IN THIS VERSION <101)
!    ET = TOTAL CENTRE-OF-MASS ENERGY
!    XM = PARTICLE MASSES ( DIM=100 )
!    P  = PARTICLE MOMENTA ( DIM=(4,100) )
!    WT = WEIGHT OF THE EVENT
!
!------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION XM(100),P(4,100),Q(4,100),Z(100),R(4)  &
        ,B(3),P2(100),XM2(100),E(100),V(100),IWARN(5) &
        ,RN(4)
      real rho(4)
      DATA ACC/1.D-14/,ITMAX/6/,IBEGIN/0/!,IWARN/5*0/
      save ACC,ITMAX,IBEGIN,Z
      save TWOPI,PI2LOG
!------------------------------------------------------
!       VERSION DEFINITION
!------------------------------------------------------
!      DATA IVERSION/0/
!      IF(IVERSION.EQ.0)PRINT*,'VERSION OF 13/12/94'
!      IVERSION=1
!------------------------------------------------------
!         IBEGIN=0
! INITIALIZATION STEP: FACTORIALS FOR THE PHASE SPACE WEIGHT
      IF(IBEGIN.NE.0) GOTO 103
      IBEGIN=1
      TWOPI=8.D0*DATAN(1.D0)
      PI2LOG=DLOG(TWOPI/4.D0)
      Z(2)=PI2LOG
      DO 101 K=3,nmax!100
  101 Z(K)=Z(K-1)+PI2LOG-2.D0*DLOG(DBLE(K-2))
      DO 102 K=3,nmax!100
  102 Z(K)=(Z(K)-DLOG(DBLE(K-1)))
!
! CHECK ON THE NUMBER OF PARTICLES
  103 IF(N.GT.1.AND.N.LT.101) GOTO 104
      PRINT 1001,N
      STOP
!
! CHECK WHETHER TOTAL ENERGY IS SUFFICIENT; COUNT NONZERO MASSES
  104 XMT=0.D0
      NM=0
      DO 105 I=1,N
       IF((XM(I).GT.0.D0).OR.(XM(I).LT.0.D0)) NM=NM+1
  105 XMT=XMT+DABS(XM(I))
      IF(XMT.LE.ET) GOTO 201
      PRINT 1002,XMT,ET
      STOP
!
! THE PARAMETER VALUES ARE NOW ACCEPTED
!
! GENERATE N MASSLESS MOMENTA IN INFINITE PHASE SPACE
 201  CONTINUE
      DO 202 I=1,N
        call avh_randarr(rho,4)
        RN(1)=rho(1)
        RN(2)=rho(2)
        RN(3)=rho(3)
        RN(4)=rho(4)
      C=2.D0*RN(1)-1.D0
      S=DSQRT(1.D0-C*C)
      F=TWOPI*RN(2)
      Q(4,I)=-DLOG(RN(3)*RN(4))
      Q(3,I)=Q(4,I)*C
      Q(2,I)=Q(4,I)*S*DSIN(F)
  202 Q(1,I)=Q(4,I)*S*DCOS(F)
!
! CALCULATE THE PARAMETERS OF THE CONFORMAL TRANSFORMATION
      DO 203 I=1,4
  203 R(I)=0.D0
      DO 204 I=1,N
      DO 204 K=1,4
  204 R(K)=R(K)+Q(K,I)
      RMAS=DSQRT(R(4)**2-R(3)**2-R(2)**2-R(1)**2)
      DO 205 K=1,3
  205 B(K)=-R(K)/RMAS
      G=R(4)/RMAS
      A=1.D0/(1.D0+G)
      X=ET/RMAS
!
! TRANSFORM THE Q'S CONFORMALLY INTO THE P'S
      DO 207 I=1,N
      BQ=B(1)*Q(1,I)+B(2)*Q(2,I)+B(3)*Q(3,I)
      DO 206 K=1,3
  206 P(K,I)=X*(Q(K,I)+B(K)*(Q(4,I)+A*BQ))
  207 P(4,I)=X*(G*Q(4,I)+BQ)
!
! CALCULATE WEIGHT AND POSSIBLE WARNINGS
      WT=PI2LOG
      IF(N.NE.2) WT=(2.D0*N-4.D0)*DLOG(ET)+Z(N)
      IF(WT.GE.-180.D0) GOTO 208
      IF(IWARN(1).LE.5) PRINT 1004,WT
      IWARN(1)=IWARN(1)+1
  208 IF(WT.LE. 174.D0) GOTO 209
      IF(IWARN(2).LE.5) PRINT 1005,WT
      IWARN(2)=IWARN(2)+1
!
! RETURN FOR WEIGHTED MASSLESS MOMENTA
  209 IF(NM.NE.0) GOTO 210
      WT=DEXP(WT)
      RETURN
!
! MASSIVE PARTICLES: RESCALE THE MOMENTA BY A FACTOR X
  210 XMAX=DSQRT(1.D0-(XMT/ET)**2)
      DO 301 I=1,N
      XM2(I)=XM(I)**2
  301 P2(I)=P(4,I)**2
      ITER=0
      X=XMAX
      ACCU=ET*ACC
  302 F0=-ET
      G0=0.D0
      X2=X*X
      DO 303 I=1,N
      E(I)=DSQRT(XM2(I)+X2*P2(I))
      F0=F0+E(I)
  303 G0=G0+P2(I)/E(I)
      IF(DABS(F0).LE.ACCU) GOTO 305
      ITER=ITER+1
      IF(ITER.LE.ITMAX) GOTO 304
      PRINT 1006,ITMAX
      GOTO 305
  304 X=X-F0/(X*G0)
      GOTO 302
  305 DO 307 I=1,N
      V(I)=X*P(4,I)
      DO 306 K=1,3
  306 P(K,I)=X*P(K,I)
  307 P(4,I)=E(I)
!
! CALCULATE THE MASS-EFFECT WEIGHT FACTOR
      WT2=1.D0
      WT3=0.D0
      DO 308 I=1,N
      WT2=WT2*V(I)/E(I)
  308 WT3=WT3+V(I)**2/E(I)
      WTM=(2.D0*N-3.D0)*DLOG(X)+DLOG(WT2/WT3*ET)
!
! RETURN FOR  WEIGHTED MASSIVE MOMENTA
      WT=WT+WTM
      IF(WT.GE.-180.D0) GOTO 309
      IF(IWARN(3).LE.5) PRINT 1004,WT
      IWARN(3)=IWARN(3)+1
  309 IF(WT.LE. 174.D0) GOTO 310
      IF(IWARN(4).LE.5) PRINT 1005,WT
      IWARN(4)=IWARN(4)+1
  310 WT=DEXP(WT)
      RETURN
!
 1001 FORMAT(' RAMBO FAILS: # OF PARTICLES =',I5,' IS NOT ALLOWED')
 1002 FORMAT(' RAMBO FAILS: TOTAL MASS =',D15.6,' IS NOT', &
       ' SMALLER THAN TOTAL ENERGY =',D15.6)
 1004 FORMAT(' RAMBO WARNS: WEIGHT = EXP(',F20.9,') MAY UNDERFLOW')
 1005 FORMAT(' RAMBO WARNS: WEIGHT = EXP(',F20.9,') MAY  OVERFLOW')
 1006 FORMAT(' RAMBO WARNS:',I3,' ITERATIONS DID NOT GIVE THE', &
       ' DESIRED ACCURACY =',D15.6)
      END subroutine
!
end module
