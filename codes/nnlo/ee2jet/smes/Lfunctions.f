c     implicit real(kind(1d0)) (a-h,o-z)
c     real(kind(1d0)) I3m3,I3m3sym
c     complex(kind(1d0)) L0,L1,Lsm1,Lsm1tilde
c     common /qsquared/ q2
c     q2 = 1d0
c     goto 1
c1    read(5,*)x,y
c     write(6,*)L0(-x,-y)
c     write(6,*)L1(-x,-y)
c     goto 1
c2    read(5,*)x1,y1,x2,y2
c     write(6,*)Lsm1(-x1,-y1,-x2,-y2)
c     goto 2
c3    read(5,*)s,t,xm12,xm22
c     write(6,*)Lsm1tilde(s,t,xm12,xm22)
c     goto 3
c4    read(5,*)s12,s34,s56
c     write(6,*)I3m3(s12,s34,s56)
c     write(6,*)I3m3sym(s12,s34,s56)
c     goto 4
c     end

************************************************************************

      REAL(KIND(1D0)) FUNCTION THETA(BOOL)
      IMPLICIT NONE

      LOGICAL BOOL

      THETA = 0D0
      IF (BOOL) THETA = 1D0

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION Ln(X,Y)
      IMPLICIT NONE

      REAL(KIND(1D0)) X,Y
      COMPLEX(KIND(1D0)) LOGX,LOGY

      REAL(KIND(1D0)) Q2
      COMMON /QSQUARED/ Q2
      REAL(KIND(1D0)) PI,THETA
      PARAMETER ( PI  = 3.14159265358979 D0)

      LOGX = CMPLX(LOG(ABS(X)/Q2),-PI*THETA(-X.GT.0D0))
      LOGY = CMPLX(LOG(ABS(Y)/Q2),-PI*THETA(-Y.GT.0D0))
      Ln = (LOGX - LOGY)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION L0(X,Y)
      IMPLICIT NONE

      REAL(KIND(1D0)) X,Y
      COMPLEX(KIND(1D0)) LOGX,LOGY

      REAL(KIND(1D0)) Q2
      COMMON /QSQUARED/ Q2
      REAL(KIND(1D0)) PI,THETA
      PARAMETER ( PI  = 3.14159265358979 D0)

      LOGX = CMPLX(LOG(ABS(X)/Q2),-PI*THETA(-X.GT.0D0))
      LOGY = CMPLX(LOG(ABS(Y)/Q2),-PI*THETA(-Y.GT.0D0))
      L0 = (LOGX - LOGY)/(1D0-X/Y)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION L1(X,Y)
      IMPLICIT NONE

      REAL(KIND(1D0)) X,Y
      COMPLEX(KIND(1D0)) L0

      REAL(KIND(1D0)) Q2
      COMMON /QSQUARED/ Q2

      L1 = (L0(X,Y) + 1D0)/(1D0-X/Y)

      END


************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION LSM1(X1,Y1,X2,Y2)
      IMPLICIT NONE

      REAL(KIND(1D0)) X1,Y1,X2,Y2,R1,R2,DDILOG
      COMPLEX(KIND(1D0)) LOGX1,LOGY1,LOGX2,LOGY2,LI2R1,LI2R2

      REAL(KIND(1D0)) Q2
      COMMON /QSQUARED/ Q2
      REAL(KIND(1D0)) PI,PI2OVER6,THETA
      PARAMETER ( PI  = 3.14159265358979 D0)
      PARAMETER ( PI2OVER6 = 1.64493406684822 D0)
      EXTERNAL DDILOG

      LOGX1 = CMPLX(LOG(ABS(X1)/Q2),-PI*THETA(-X1.GT.0D0))
      LOGY1 = CMPLX(LOG(ABS(Y1)/Q2),-PI*THETA(-Y1.GT.0D0))
      LOGX2 = CMPLX(LOG(ABS(X2)/Q2),-PI*THETA(-X2.GT.0D0))
      LOGY2 = CMPLX(LOG(ABS(Y2)/Q2),-PI*THETA(-Y2.GT.0D0))
      R1 = X1/Y1
      IF (R1.GT.0.D0) THEN
         LI2R1 = CMPLX(DDILOG(1D0 - R1),0D0)
      ELSE
         LI2R1 = CMPLX(PI2OVER6 - DDILOG(R1),0D0)
     >         - CMPLX(LOG(1D0 - R1),0D0)*(LOGX1 - LOGY1)
      ENDIF
      R2 = X2/Y2
      IF (R2.GT.0.D0) THEN
         LI2R2 = CMPLX(DDILOG(1D0 - R2),0D0)
      ELSE
         LI2R2 = CMPLX(PI2OVER6 - DDILOG(R2),0D0)
     >         - CMPLX(LOG(1D0 - R2),0D0)*(LOGX2 - LOGY2)
      ENDIF

      LSM1 = LI2R1 + LI2R2 + (LOGX1-LOGY1)*(LOGX2-LOGY2) - PI2OVER6

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION LS1(X1,Y1,X2,Y2)
      IMPLICIT NONE

      REAL(KIND(1D0)) X1,Y1,X2,Y2,ONEMR
      COMPLEX(KIND(1D0)) LSM1,L0

      ONEMR = 1D0 - X1/Y1 - X2/Y2

      LS1 = (LSM1(X1,Y1,X2,Y2)/ONEMR + L0(X1,Y1) + L0(X2,Y2))/ONEMR

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION LSM1H(S,T,SM1,SM2)
      IMPLICIT NONE

      REAL(KIND(1D0)) S,T,SM1,SM2,I3M3
      COMPLEX(KIND(1D0)) LSM1TILDE

      LSM1H = LSM1TILDE(S,T,SM1,SM2)
     > + (0.5D0*(S - SM1 - SM2) + SM1*SM2/T)*I3M3(S,SM1,SM2)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION LSM1TILDE(S,T,SM1,SM2)
      IMPLICIT NONE

      REAL(KIND(1D0)) S,T,SM1,SM2,R1,R2,DDILOG
      COMPLEX(KIND(1D0)) LOGS,LOGT,LOGSM1,LOGSM2,LI2R1,LI2R2

      REAL(KIND(1D0)) Q2
      COMMON /QSQUARED/ Q2
      REAL(KIND(1D0)) PI,PI2OVER6,THETA
      PARAMETER ( PI  = 3.14159265358979 D0)
      PARAMETER ( PI2OVER6 = 1.64493406684822 D0)
      EXTERNAL DDILOG

      LOGS = CMPLX(LOG(ABS(S)/Q2),-PI*THETA(-S.GT.0D0))
      LOGT = CMPLX(LOG(ABS(T)/Q2),-PI*THETA(-T.GT.0D0))
      LOGSM1 = CMPLX(LOG(ABS(SM1)/Q2),-PI*THETA(-SM1.GT.0D0))
      LOGSM2 = CMPLX(LOG(ABS(SM2)/Q2),-PI*THETA(-SM2.GT.0D0))
      R1 = SM1/T
      IF (R1.GT.0.D0) THEN
         LI2R1 = CMPLX(DDILOG(1D0 - R1),0D0)
      ELSE
         LI2R1 = CMPLX(PI2OVER6 - DDILOG(R1),0D0)
     >         - CMPLX(LOG(1D0 - R1),0D0)*(LOGSM1 - LOGT)
      ENDIF
      R2 = SM2/T
      IF (R2.GT.0.D0) THEN
         LI2R2 = CMPLX(DDILOG(1D0 - R2),0D0)
      ELSE
         LI2R2 = CMPLX(PI2OVER6 - DDILOG(R2),0D0)
     >         - CMPLX(LOG(1D0 - R2),0D0)*(LOGSM2 - LOGT)
      ENDIF

      LSM1TILDE = - LI2R1 - LI2R2 - 0.5D0*(LOGS-LOGT)**2
     > + 0.5D0*(LOGS-LOGSM1)*(LOGS-LOGSM2)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION LSM1E(S,T,SM1,SM3)
      IMPLICIT NONE

      REAL(KIND(1D0)) S,T,SM1,SM3,R1S,R1T,R3S,R3T,DDILOG
      COMPLEX(KIND(1D0)) LOGS,LOGT,LOGSM1,LOGSM3
     >          ,LI2R1S,LI2R1T,LI2R3S,LI2R3T,LI2R1SR3T

      REAL(KIND(1D0)) Q2
      COMMON /QSQUARED/ Q2
      REAL(KIND(1D0)) PI,PI2OVER6,THETA
      PARAMETER ( PI  = 3.14159265358979 D0)
      PARAMETER ( PI2OVER6 = 1.64493406684822 D0)
      EXTERNAL DDILOG

      LOGS = CMPLX(LOG(ABS(S)/Q2),-PI*THETA(-S.GT.0D0))
      LOGT = CMPLX(LOG(ABS(T)/Q2),-PI*THETA(-T.GT.0D0))
      LOGSM1 = CMPLX(LOG(ABS(SM1)/Q2),-PI*THETA(-SM1.GT.0D0))
      LOGSM3 = CMPLX(LOG(ABS(SM3)/Q2),-PI*THETA(-SM3.GT.0D0))
      R1S = SM1/S
      IF (R1S.GT.0.D0) THEN
         LI2R1S = CMPLX(DDILOG(1D0 - R1S),0D0)
      ELSE
         LI2R1S = CMPLX(PI2OVER6 - DDILOG(R1S),0D0)
     >          - CMPLX(LOG(1D0 - R1S),0D0)*(LOGSM1 - LOGS)
      ENDIF
      R1T = SM1/T
      IF (R1T.GT.0.D0) THEN
         LI2R1T = CMPLX(DDILOG(1D0 - R1T),0D0)
      ELSE
         LI2R1T = CMPLX(PI2OVER6 - DDILOG(R1T),0D0)
     >          - CMPLX(LOG(1D0 - R1T),0D0)*(LOGSM1 - LOGT)
      ENDIF
      R3S = SM3/S
      IF (R3S.GT.0.D0) THEN
         LI2R3S = CMPLX(DDILOG(1D0 - R3S),0D0)
      ELSE
         LI2R3S = CMPLX(PI2OVER6 - DDILOG(R3S),0D0)
     >          - CMPLX(LOG(1D0 - R3S),0D0)*(LOGSM3 - LOGS)
      ENDIF
      R3T = SM3/T
      IF (R3T.GT.0.D0) THEN
         LI2R3T = CMPLX(DDILOG(1D0 - R3T),0D0)
      ELSE
         LI2R3T = CMPLX(PI2OVER6 - DDILOG(R3T),0D0)
     >          - CMPLX(LOG(1D0 - R3T),0D0)*(LOGSM3 - LOGT)
      ENDIF
      IF (R1S*R3T.GT.0.D0) THEN
         LI2R1SR3T = CMPLX(DDILOG(1D0 - R1S*R3T),0D0)
      ELSE
         LI2R1SR3T = CMPLX(PI2OVER6 - DDILOG(R1S*R3T),0D0)
     >             - CMPLX(LOG(1D0 - R1S*R3T),0D0)
     >              *(LOGSM1 + LOGSM3 - LOGS - LOGT)
      ENDIF

      LSM1E = - LI2R1S - LI2R1T - LI2R3S - LI2R3T + LI2R1SR3T
     > - 0.5D0*(LOGS-LOGT)**2

      END

************************************************************************

      REAL(KIND(1D0)) FUNCTION GRAM(S12,S34,S56)
      IMPLICIT NONE

      REAL(KIND(1D0)) S12,S34,S56

      GRAM = 2D0*(S12*S34+S34*S56+S56*S12) - (S12**2+S34**2+S56**2 )
      END

************************************************************************

      REAL(KIND(1D0)) FUNCTION RELI2(R)
      IMPLICIT NONE

      REAL(KIND(1D0)) R,DDILOG
      REAL(KIND(1D0)) PI2OVER6
      PARAMETER ( PI2OVER6 = 1.64493406684822 D0)
      EXTERNAL DDILOG

      IF (R.LT.1D0) THEN
       RELI2 = DDILOG(R)
      ELSE
       RELI2 =PI2OVER6 - DDILOG(1D0-R) - LOG(ABS(R))*LOG(ABS(1D0-R))
      ENDIF

      END

************************************************************************

      REAL(KIND(1D0)) FUNCTION I3M3(S12,S34,S56)
      IMPLICIT NONE

      REAL(KIND(1D0)) S12,S34,S56,D12,D34,D56,DELTA3,RD3
     > ,DCLAUS,RELI2,GRAM,X,Y,RHO

      REAL(KIND(1D0)) PI2,PI2OVER3,THETA
      PARAMETER ( PI2= 9.86960440108936 D0)
      PARAMETER ( PI2OVER3 = 3.28986813369645 D0)

      DELTA3 = -GRAM(S12,S34,S56)
      D12 = S12 - S34 - S56
      D56 = S56 - S12 - S34

      IF (DELTA3.GT.0) GOTO 2

      D34 = S34 - S56 - S12
      RD3 = SQRT(-DELTA3)
      I3M3 = 2D0/RD3
     > *( DCLAUS(2D0*ATAN(RD3/D12))
     >  + DCLAUS(2D0*ATAN(RD3/D34))
     >  + DCLAUS(2D0*ATAN(RD3/D56)) )

      RETURN

 2    RD3 = SQRT(DELTA3)
      X = S12/S56
      Y = S34/S56
      RHO = 2D0*S56/(D56+RD3)
      I3M3 = -(2D0*(RELI2(-RHO*X) + RELI2(-RHO*Y)) 
     > + LOG(ABS(RHO*X))*LOG(ABS(RHO*Y)) - PI2*THETA(RHO.LT.0D0)
     > + LOG(Y/X)*LOG(ABS((1D0+RHO*Y)/(1D0+RHO*X))) + PI2OVER3)/RD3

      END

************************************************************************

      REAL(KIND(1D0)) FUNCTION I3M3SYM(S12,S34,S56)
      IMPLICIT NONE

      REAL(KIND(1D0)) S12,S34,S56,D12,D34,D56,DELTA3,RD3,DCLAUS,GRAM,
     >                X,Y,RHO,I3M3

      DELTA3 = -GRAM(S12,S34,S56)
      D12 = S12 - S34 - S56
      D56 = S56 - S12 - S34

      IF (DELTA3.GT.0) GOTO 2

      D34 = S34 - S56 - S12
      RD3 = SQRT(-DELTA3)
      I3M3SYM = 2D0/RD3
     > *( DCLAUS(2D0*ATAN(RD3/D12))
     >  + DCLAUS(2D0*ATAN(RD3/D34))
     >  + DCLAUS(2D0*ATAN(RD3/D56)) )

      RETURN

 2    RD3 = SQRT(DELTA3)
      X = S12/S56
      Y = S34/S56
      RHO = 2D0*S56/(D56+RD3)
      IF (RHO.GT.0D0) THEN
       I3M3SYM = I3M3(S12,S34,S56)
      ELSE
       I3M3SYM = I3M3(S56,S34,S12)
      ENDIF

      END

