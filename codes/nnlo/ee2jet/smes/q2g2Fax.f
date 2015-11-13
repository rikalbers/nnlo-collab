************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbggFaxPMMP(ONE,TWO,THREE,FOUR,FIVE,
     >                                         SIX,A,B)
      IMPLICIT NONE 

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      REAL(KIND(1D0)) S(11,11),XMT,mtsq
      COMMON /DOTPRODUCTS/ S
     >       /TOPMASS/ XMT,mtsq
      REAL(KIND(1D0)) S12, S34, S56, T123, T124, T356, T456, I3m3
      COMPLEX(KIND(1D0)) A12,A13,A14,A15,A23,A24,A25,A34,A35,A36,A45,A46
     >          ,A56,B12,B13,B14,B16,B23,B24,B26,B34,B35,B36,B45,B46,B56
     >          ,L0, L1, Ln, Lsm1h

      A12 = A(ONE,TWO)
      A13 = A(ONE,THREE)
      A14 = A(ONE,FOUR)
      A15 = A(ONE,FIVE)
      A23 = A(TWO,THREE)
      A24 = A(TWO,FOUR)
      A25 = A(TWO,FIVE)
      A34 = A(THREE,FOUR)
      A35 = A(THREE,FIVE)
      A36 = A(THREE,SIX)
      A45 = A(FOUR,FIVE)
      A46 = A(FOUR,SIX)
      A56 = A(FIVE,SIX)
      B12 = B(ONE,TWO)
      B13 = B(ONE,THREE)
      B14 = B(ONE,FOUR)
      B16 = B(ONE,SIX)
      B23 = B(TWO,THREE)
      B24 = B(TWO,FOUR)
      B26 = B(TWO,SIX)
      B34 = B(THREE,FOUR)
      B35 = B(THREE,FIVE)
      B36 = B(THREE,SIX)
      B45 = B(FOUR,FIVE)
      B46 = B(FOUR,SIX)
      B56 = B(FIVE,SIX)
      S12 = S(ONE,TWO)
      S34 = S(THREE,FOUR)
      S56 = S(FIVE,SIX)
      T123 = S(ONE,TWO)+S(ONE,THREE)+S(TWO,THREE)
      T124 = S(ONE,TWO)+S(ONE,FOUR)+S(TWO,FOUR)
      T356 = S(THREE,FIVE)+S(THREE,SIX)+S(FIVE,SIX)
      T456 = S(FOUR,FIVE)+S(FOUR,SIX)+S(FIVE,SIX)

      qqbggFaxPMMP = 
     >  -((A23**2*A45*B16)/(A24*A34*(-(A14*B13)-A24*B23)*S56))+
     >  (A25*B14**2*B36)/(B13*(-(A14*B13)-A24*B23)*B34*S56)+
     >  (A23**2*A35*B16)/(12.*A24*A34*mtsq*S56)-
     >  (A25*B14**2*B46)/(12.*B13*B34*mtsq*S56)-
     >  (A23*B46*(-((A23*A45)/(A12*A34*S56))-
     >       (2*A45*B14)/
     >        (S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+S56**2)
     >       +(B14*B36*(S12-S34-S56))/
     >        (B34*B56*(S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+
     >            S56**2))+(A23*B36*(-S12+S34-S56))/
     >        (A12*B56*(S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+
     >            S56**2))+(A23*A45*(-S12-S34+S56))/
     >        (A12*A34*(S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+
     >            S56**2))))/(-(A14*B13)-A24*B23)-
     >  (A35*B14*(-((B14*B36)/(B12*B34*S56))-
     >       (2*A23*B36)/
     >        (S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+S56**2)
     >       +(A23*A45*(S12-S34-S56))/
     >        (A34*A56*(S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+
     >            S56**2))+(A45*B14*(-S12+S34-S56))/
     >        (A56*B12*(S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+
     >            S56**2))+(B14*B36*(-S12-S34+S56))/
     >        (B12*B34*(S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+
     >            S56**2))))/(-(A14*B13)-A24*B23)-
     >  (-((A23*A35*B14*B46)/
     >        (S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+S56**2))
     >     -(3*(-(A13*B14)-A23*B24)*
     >        (-(A12*A25*B12*B16)-A23*A45*B14*B36-A24*A35*B13*B46-
     >          A25*A56*B16*B56))/
     >      (2.*(-(A14*B13)-A24*B23)*
     >        (S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+S56**2))
     >     -(3*(A25*B12+A35*B13)*(-(A13*B14)-A23*B24)*
     >        (-S12+S34-S56)*
     >        (-(A12*B16*(S12-S34-S56))-A25*B56*(-S12-S34+S56)))/
     >      ((-(A14*B13)-A24*B23)*
     >        (S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+S56**2)**
     >         2)-(A23*(A25*B12+A35*B13)*B46)/
     >      (2.*(-(A14*B13)-A24*B23)*T123)+
     >     (A24*(A25*B12+A35*B13)*(-(A13*B14)-A23*B24)*B36*(T123-T124))/
     >      ((-(A14*B13)-A24*B23)**2*
     >        (S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+S56**2)))
     >    *I3m3(S12,S34,S56)-(-((A23*A35*B14*B46)/
     >        (S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+S56**2))
     >     -(3*(-(A13*B14)-A23*B24)*
     >        (-(A12*A25*B12*B16)-A23*A45*B14*B36-A24*A35*B13*B46-
     >          A25*A56*B16*B56))/
     >      (2.*(-(A14*B13)-A24*B23)*
     >        (S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+S56**2))
     >     -(3*(-(A13*B14)-A23*B24)*(-(A12*B16)+A24*B46)*
     >        (-S12+S34-S56)*
     >        (A25*B12*(S12-S34-S56)+A56*B16*(-S12-S34+S56)))/
     >      ((-(A14*B13)-A24*B23)*
     >        (S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+S56**2)**
     >         2)-(A35*B14*(-(A12*B16)+A24*B46))/
     >      (2.*(-(A14*B13)-A24*B23)*T124)+
     >     (A45*B13*(-(A13*B14)-A23*B24)*(-(A12*B16)+A24*B46)*
     >        (-T123+T124))/
     >      ((-(A14*B13)-A24*B23)**2*
     >        (S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+S56**2)))
     >    *I3m3(S12,S34,S56)-(A23*(A25*B12+A35*B13)*
     >     (-(A15*B13)-A25*B23)*L0(-T123,-S12))/
     >   (A56*(-(A14*B13)-A24*B23)**2*S12)+
     >  ((A24*B12+A34*B13)*(A25*B12+A35*B13)*B14*(-(A15*B13)-A25*B23)*
     >     L0(-T123,-S56))/(A56*B12*B13*(-(A14*B13)-A24*B23)**2*S56)+
     >  (B14*(-(A14*B16)-A24*B26)*(-(A12*B16)+A24*B46)*L0(-T124,-S12))/
     >   ((-(A14*B13)-A24*B23)**2*B56*S12)+
     >  (A23*(-(A14*B16)-A24*B26)*(-(A12*B13)-A24*B34)*
     >     (-(A12*B16)+A24*B46)*L0(-T124,-S56))/
     >   (A12*A24*(-(A14*B13)-A24*B23)**2*B56*S56)+
     >  (A45*(A25*B12+A35*B13)*B14**2*L1(-S56,-T123))/
     >   (A56*B12*B13*(-(A14*B13)-A24*B23)*T123)+
     >  (A23**2*B36*(-(A12*B16)+A24*B46)*L1(-S56,-T124))/
     >   (A12*A24*(-(A14*B13)-A24*B23)*B56*T124)-
     >  ((A23*A45*(A25*B12+A35*B13))/(A34*A56*(-(A14*B13)-A24*B23)**2)+
     >     (6*A12*(A25*B12+A35*B13)*(-(A13*B14)-A23*B24)*
     >        (-2*A25*B12*B56+B16*(-S12+S34-S56)))/
     >      ((-(A14*B13)-A24*B23)*
     >        (S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+S56**2)**
     >         2)+(A23*A35**2*B14)/
     >      (A34*A56*(S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+
     >          S56**2))-(A24*(A25*B12+A35*B13)*(-(A13*B14)-A23*B24)*
     >        (3*(-(A14*A35*B13)-A24*A35*B23)-A45*(T123-T124)))/
     >      (A34*A56*(-(A14*B13)-A24*B23)**2*
     >        (S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+S56**2)))
     >    *Ln(-S12,-S34)-((B14*B36*(-(A12*B16)+A24*B46))/
     >      ((-(A14*B13)-A24*B23)**2*B34*B56)-
     >     (6*B12*(-(A13*B14)-A23*B24)*(-(A12*B16)+A24*B46)*
     >        (-2*A12*A56*B16+A25*(-S12+S34-S56)))/
     >      ((-(A14*B13)-A24*B23)*
     >        (S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+S56**2)**
     >         2)+(A23*B14*B46**2)/
     >      (B34*B56*(S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+
     >          S56**2))-(B13*(-(A13*B14)-A23*B24)*(-(A12*B16)+A24*B46)*
     >        (3*(-(A14*B13*B46)-A24*B23*B46)-B36*(-T123+T124)))/
     >      ((-(A14*B13)-A24*B23)**2*B34*B56*
     >        (S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+S56**2)))
     >    *Ln(-S12,-S34)-((A25*B12+A35*B13)**2*Ln(-S56,-S34))/
     >   (A56*B12*(-(A14*B13)-A24*B23)**2)-
     >  ((-(A12*B16)+A24*B46)**2*Ln(-S56,-S34))/
     >   (A12*(-(A14*B13)-A24*B23)**2*B56)-
     >  (-((A24*A35*(A23*B36+A25*B56))/
     >        (A12*A34*(-(A45*B35)-A46*B36)**2))-
     >     (6*A56*(-(A35*B45)-A36*B46)*(A23*B36+A25*B56)*
     >        (2*A25*B12*B56-B16*(-S12+S34-S56)))/
     >      ((-(A45*B35)-A46*B36)*
     >        (S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+S56**2)**
     >         2)-(A23**2*A35*B46)/
     >      (A12*A34*(S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+
     >          S56**2))-(A45*(-(A35*B45)-A36*B46)*(A23*B36+A25*B56)*
     >        (3*(A23*A45*B35+A23*A46*B36)+A24*(T356-T456)))/
     >      (A12*A34*(-(A45*B35)-A46*B36)**2*
     >        (S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+S56**2)))
     >    *Ln(-S56,-S34)-(-((B13*(A45*B14-A56*B16)*B46)/
     >        (B12*B34*(-(A45*B35)-A46*B36)**2))+
     >     (6*(A45*B14-A56*B16)*(-(A35*B45)-A36*B46)*B56*
     >        (2*A12*A56*B16-A25*(-S12+S34-S56)))/
     >      ((-(A45*B35)-A46*B36)*
     >        (S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+S56**2)**
     >         2)-(A35*B14**2*B46)/
     >      (B12*B34*(S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+
     >          S56**2))-((A45*B14-A56*B16)*B36*(-(A35*B45)-A36*B46)*
     >        (3*(A45*B14*B35+A46*B14*B36)+B13*(-T356+T456)))/
     >      (B12*B34*(-(A45*B35)-A46*B36)**2*
     >        (S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+S56**2)))
     >    *Ln(-S56,-S34)+(((A24*B12+A34*B13)**2*
     >        (-(A15*B13)-A25*B23)**2-A45**2*B13**2*T123**2)*
     >    Lsm1h(S34,T123,S12,S56))/(2.*A56*B12*(-(A14*B13)-A24*B23)**4)+
     >  (((-(A14*B16)-A24*B26)**2*(-(A12*B13)-A24*B34)**2-
     >       A24**2*B36**2*T124**2)*Lsm1h(S34,T124,S12,S56))/
     >   (2.*A12*(-(A14*B13)-A24*B23)**4*B56)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbggFaxPMPM(ONE,TWO,THREE,FOUR,FIVE,
     >                                         SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      REAL(KIND(1D0)) S(11,11),XMT,mtsq
      COMMON /DOTPRODUCTS/ S
     >       /TOPMASS/ XMT,mtsq
      REAL(KIND(1D0)) S12, S34, S56, T123, T124, T356, T456, I3m3
      COMPLEX(KIND(1D0)) A12,A13,A14,A15,A23,A24,A25,A34,A35,A36,A45,A46
     >          ,A56,B12,B13,B14,B16,B23,B24,B26,B34,B35,B36,B45,B46,B56
     >          ,L0, L1, Ln, Lsm1h

      A12 = A(ONE,TWO)
      A13 = A(ONE,THREE)
      A14 = A(ONE,FOUR)
      A15 = A(ONE,FIVE)
      A23 = A(TWO,THREE)
      A24 = A(TWO,FOUR)
      A25 = A(TWO,FIVE)
      A34 = A(THREE,FOUR)
      A35 = A(THREE,FIVE)
      A36 = A(THREE,SIX)
      A45 = A(FOUR,FIVE)
      A46 = A(FOUR,SIX)
      A56 = A(FIVE,SIX)
      B12 = B(ONE,TWO)
      B13 = B(ONE,THREE)
      B14 = B(ONE,FOUR)
      B16 = B(ONE,SIX)
      B23 = B(TWO,THREE)
      B24 = B(TWO,FOUR)
      B26 = B(TWO,SIX)
      B34 = B(THREE,FOUR)
      B35 = B(THREE,FIVE)
      B36 = B(THREE,SIX)
      B45 = B(FOUR,FIVE)
      B46 = B(FOUR,SIX)
      B56 = B(FIVE,SIX)
      S12 = S(ONE,TWO)
      S34 = S(THREE,FOUR)
      S56 = S(FIVE,SIX)
      T123 = S(ONE,TWO)+S(ONE,THREE)+S(TWO,THREE)
      T124 = S(ONE,TWO)+S(ONE,FOUR)+S(TWO,FOUR)
      T356 = S(THREE,FIVE)+S(THREE,SIX)+S(FIVE,SIX)
      T456 = S(FOUR,FIVE)+S(FOUR,SIX)+S(FIVE,SIX)
      
c No top quark...AK
      qqbggFaxPMPM =
     >  (A24*A35*(-(A14*B16)-A34*B36))/
     >   (A13*A34*(-(A13*B14)-A23*B24)*S56)-
     >  (B13*(-(A25*B23)+A45*B34)*B46)/
     >   (B24*(-(A13*B14)-A23*B24)*B34*S56)+
     >  (B13*(-(A25*B23)+A45*B34)*B36)/
     >   (12.*B24*B34*mtsq*S56)-
     >  (A24*A45*(-(A14*B16)-A34*B36))/
     >   (12.*A13*A34*mtsq*S56)+
     >  (A24*B36*((A24*A35)/(A12*A34*S56)-
     >       (2*A35*B13)/
     >        (S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+S56**2)
     >       -(B13*B46*(S12-S34-S56))/
     >        (B34*B56*(S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+
     >            S56**2))+(A24*B46*(-S12+S34-S56))/
     >        (A12*B56*(S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+
     >            S56**2))-(A24*A35*(-S12-S34+S56))/
     >        (A12*A34*(S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+
     >            S56**2))))/(-(A13*B14)-A23*B24)+
     >  (A45*B13*((B13*B46)/(B12*B34*S56)-
     >       (2*A24*B46)/
     >        (S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+S56**2)
     >       -(A24*A35*(S12-S34-S56))/
     >        (A34*A56*(S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+
     >            S56**2))+(A35*B13*(-S12+S34-S56))/
     >        (A56*B12*(S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+
     >            S56**2))-(B13*B46*(-S12-S34+S56))/
     >        (B12*B34*(S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+
     >            S56**2))))/(-(A13*B14)-A23*B24)+
     >  (-((A24*A45*B13*B36)/
     >        (S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+S56**2))
     >     -(3*(-(A14*B13)-A24*B23)*
     >        (-(A12*A25*B12*B16)-A23*A45*B14*B36-A24*A35*B13*B46-
     >          A25*A56*B16*B56))/
     >      (2.*(-(A13*B14)-A23*B24)*
     >        (S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+S56**2))
     >     -(3*(-(A14*B13)-A24*B23)*(-(A12*B16)+A23*B36)*
     >        (-S12+S34-S56)*
     >        (A25*B12*(S12-S34-S56)+A56*B16*(-S12-S34+S56)))/
     >      ((-(A13*B14)-A23*B24)*
     >        (S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+S56**2)**
     >         2)-(A45*B13*(-(A12*B16)+A23*B36))/
     >      (2.*(-(A13*B14)-A23*B24)*T123)+
     >     (A35*B14*(-(A14*B13)-A24*B23)*(-(A12*B16)+A23*B36)*
     >        (T123-T124))/
     >      ((-(A13*B14)-A23*B24)**2*
     >        (S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+S56**2)))
     >    *I3m3(S12,S34,S56)+(-((A24*A45*B13*B36)/
     >        (S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+S56**2))
     >     -(3*(-(A14*B13)-A24*B23)*
     >        (-(A12*A25*B12*B16)-A23*A45*B14*B36-A24*A35*B13*B46-
     >          A25*A56*B16*B56))/
     >      (2.*(-(A13*B14)-A23*B24)*
     >        (S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+S56**2))
     >     -(3*(A25*B12+A45*B14)*(-(A14*B13)-A24*B23)*
     >        (-S12+S34-S56)*
     >        (-(A12*B16*(S12-S34-S56))-A25*B56*(-S12-S34+S56)))/
     >      ((-(A13*B14)-A23*B24)*
     >        (S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+S56**2)**
     >         2)-(A24*(A25*B12+A45*B14)*B36)/
     >      (2.*(-(A13*B14)-A23*B24)*T124)+
     >     (A23*(A25*B12+A45*B14)*(-(A14*B13)-A24*B23)*B46*
     >        (-T123+T124))/
     >      ((-(A13*B14)-A23*B24)**2*
     >        (S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+S56**2)))
     >    *I3m3(S12,S34,S56)-(B13*(-(A13*B16)-A23*B26)*
     >     (-(A12*B16)+A23*B36)*L0(-T123,-S12))/
     >   ((-(A13*B14)-A23*B24)**2*B56*S12)-
     >  (A24*(-(A13*B16)-A23*B26)*(A12*B24+A13*B34)*
     >     (-(A12*B16)+A23*B36)*L0(-T123,-S56))/
     >   (A12*A13*(-(A13*B14)-A23*B24)**2*B56*S56)+
     >  (A24*(A25*B12+A45*B14)*(-(A15*B14)-A25*B24)*L0(-T124,-S12))/
     >   (A56*(-(A13*B14)-A23*B24)**2*S12)-
     >  (B13*(A25*B12+A45*B14)*(-(A15*B14)-A25*B24)*
     >     (-(A13*B12)-A34*B24)*L0(-T124,-S56))/
     >   (A56*B12*B24*(-(A13*B14)-A23*B24)**2*S56)-
     >  (A14*A24*(-(A12*B16)+A23*B36)*B46*L1(-S56,-T123))/
     >   (A12*A13*(-(A13*B14)-A23*B24)*B56*T123)-
     >  (A35*B13*(A25*B12+A45*B14)*B23*L1(-S56,-T124))/
     >   (A56*B12*B24*(-(A13*B14)-A23*B24)*T124)+
     >  (-((B13*(-(A12*B16)+A23*B36)*B46)/
     >        ((-(A13*B14)-A23*B24)**2*B34*B56))-
     >     (6*B12*(-(A14*B13)-A24*B23)*(-(A12*B16)+A23*B36)*
     >        (-2*A12*A56*B16+A25*(-S12+S34-S56)))/
     >      ((-(A13*B14)-A23*B24)*
     >        (S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+S56**2)**
     >         2)-(A24*B13*B36**2)/
     >      (B34*B56*(S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+
     >          S56**2))+(B14*(-(A14*B13)-A24*B23)*(-(A12*B16)+A23*B36)*
     >        (3*(-(A13*B14*B36)-A23*B24*B36)-B46*(T123-T124)))/
     >      ((-(A13*B14)-A23*B24)**2*B34*B56*
     >        (S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+S56**2)))
     >    *Ln(-S12,-S34)+(-((A24*A35*(A25*B12+A45*B14))/
     >        (A34*A56*(-(A13*B14)-A23*B24)**2))+
     >     (6*A12*(A25*B12+A45*B14)*(-(A14*B13)-A24*B23)*
     >        (-2*A25*B12*B56+B16*(-S12+S34-S56)))/
     >      ((-(A13*B14)-A23*B24)*
     >        (S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+S56**2)**
     >         2)-(A24*A45**2*B13)/
     >      (A34*A56*(S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+
     >          S56**2))+(A23*(A25*B12+A45*B14)*(-(A14*B13)-A24*B23)*
     >        (3*(-(A13*A45*B14)-A23*A45*B24)-A35*(-T123+T124)))/
     >      (A34*A56*(-(A13*B14)-A23*B24)**2*
     >        (S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+S56**2)))
     >    *Ln(-S12,-S34)+((A25*B12+A45*B14)**2*Ln(-S56,-S34))/
     >   (A56*B12*(-(A13*B14)-A23*B24)**2)+
     >  ((-(A12*B16)+A23*B36)**2*Ln(-S56,-S34))/
     >   (A12*(-(A13*B14)-A23*B24)**2*B56)+
     >  ((B14*(A35*B13-A56*B16)*B36)/(B12*B34*(-(A35*B45)-A36*B46)**2)+
     >     (6*(A35*B13-A56*B16)*(-(A45*B35)-A46*B36)*B56*
     >        (2*A12*A56*B16-A25*(-S12+S34-S56)))/
     >      ((-(A35*B45)-A36*B46)*
     >        (S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+S56**2)**
     >         2)+(A45*B13**2*B36)/
     >      (B12*B34*(S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+
     >          S56**2))+((A35*B13-A56*B16)*(-(A45*B35)-A46*B36)*B46*
     >        (3*(A35*B13*B45+A36*B13*B46)+B14*(T356-T456)))/
     >      (B12*B34*(-(A35*B45)-A36*B46)**2*
     >        (S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+S56**2)))
     >    *Ln(-S56,-S34)+((A23*A45*(A24*B46+A25*B56))/
     >      (A12*A34*(-(A35*B45)-A36*B46)**2)-
     >     (6*A56*(-(A45*B35)-A46*B36)*(A24*B46+A25*B56)*
     >        (2*A25*B12*B56-B16*(-S12+S34-S56)))/
     >      ((-(A35*B45)-A36*B46)*
     >        (S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+S56**2)**
     >         2)+(A24**2*A45*B36)/
     >      (A12*A34*(S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+
     >          S56**2))+(A35*(-(A45*B35)-A46*B36)*(A24*B46+A25*B56)*
     >        (3*(A24*A35*B45+A24*A36*B46)+A23*(-T356+T456)))/
     >      (A12*A34*(-(A35*B45)-A36*B46)**2*
     >        (S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+S56**2)))
     >    *Ln(-S56,-S34)-(((-(A13*B16)-A23*B26)**2*
     >        (-(A12*B14)+A23*B34)**2-A23**2*B46**2*T123**2)*
     >    Lsm1h(S34,T123,S12,S56))/(2.*A12*(-(A13*B14)-A23*B24)**4*B56)-
     >  (((A23*B12-A34*B14)**2*(-(A15*B14)-A25*B24)**2-
     >       A35**2*B14**2*T124**2)*Lsm1h(S34,T124,S12,S56))/
     >   (2.*A56*B12*(-(A13*B14)-A23*B24)**4)

      END                                                   

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbggFaxPMPP(ONE,TWO,THREE,FOUR,FIVE,
     >                                         SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      REAL(KIND(1D0)) S(11,11),XMT,mtsq
      COMMON /DOTPRODUCTS/ S
     >       /TOPMASS/ XMT,mtsq
      REAL(KIND(1D0)) S12, S14, S34, S56, T123, T124, T134
      COMPLEX(KIND(1D0)) A12, A13, A23, A24, A25, A34, A35, A45, A56
     >          ,B13, B14, B36, B46
     >          ,L0, L1, Ln, Lsm1e

      A12 = A(ONE,TWO)
      A13 = A(ONE,THREE)
      A23 = A(TWO,THREE)
      A24 = A(TWO,FOUR)
      A25 = A(TWO,FIVE)
      A34 = A(THREE,FOUR)
      A35 = A(THREE,FIVE)
      A45 = A(FOUR,FIVE)
      A56 = A(FIVE,SIX)
      B13 = B(ONE,THREE)
      B14 = B(ONE,FOUR)
      B36 = B(THREE,SIX)
      B46 = B(FOUR,SIX)
      S12 = S(ONE,TWO)
      S14 = S(ONE,FOUR)
      S34 = S(THREE,FOUR)
      S56 = S(FIVE,SIX)
      T123 = S(ONE,TWO)+S(ONE,THREE)+S(TWO,THREE)
      T124 = S(ONE,TWO)+S(ONE,FOUR)+S(TWO,FOUR)
      T134 = S(ONE,THREE)+S(ONE,FOUR)+S(THREE,FOUR)

      qqbggFaxPMPP =
     -   -(A25*(((-(A23*B13) - A24*B14)*B36)/A24 + (B46*T134)/A13))/
     -    (12.*A34*mtsq*S56) - 
     -   (A25*A35*B13*L0(-T123,-S12))/(A34**2*A56*S12) + 
     -   (A25*A45*B14*L0(-T124,-S12))/(A34**2*A56*S12) + 
     -   (A25*B46*(S14 + S34)*L1(-T123,-S56))/(A13*A34*S56**2) - 
     -   (A24*A25*B46*(L0(-T123,-S56)/S56+(S34*L1(-T123,-S56))/S56**2))/
     -  (A12*A34**2)-(A23*A25*B13*B36*L1(-T124,-S56))/(A24*A34*S56**2) + 
     -   (A23*A25*B36*(L0(-T124,-S56)/S56+(S34*L1(-T124,-S56))/S56**2))/
     -    (A12*A34**2) + (A25**2*Ln(-T123,-S56))/(A12*A34**2*A56) - 
     -   (A25**2*Ln(-T124,-S56))/(A12*A34**2*A56) - 
     -   (A25*(A24*A35 + A23*A45)*Lsm1e(T123,T124,S12,S56))/
     -    (2.*A12*A34**3*A56)

      END                                                   

************************************************************************
      
      COMPLEX(KIND(1D0)) FUNCTION qqbggFaxPMMM(ONE,TWO,THREE,FOUR,FIVE,
     >                                         SIX,A,B) 
      IMPLICIT NONE               

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qqbggFaxPMPP
      
      qqbggFaxPMMM = -qqbggFaxPMPP(TWO,ONE,FOUR,THREE,SIX,FIVE,B,A)
      
      END              
      
************************************************************************

      COMPLEX(KIND(1D0)) FUNCTIONqqbggFaxslPMPM(ONE,TWO,THREE,FOUR,FIVE,
     >                                          SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      REAL(KIND(1D0)) S(11,11),XMT,mtsq
      COMMON /DOTPRODUCTS/ S
     >       /TOPMASS/ XMT,mtsq
      REAL(KIND(1D0)) S56, T123, T124
      COMPLEX(KIND(1D0)) A12, A13, A23, A24, A25, A45
     >          ,B12, B13, B14, B16, B24, B36
     >          ,L1

      A12 = A(ONE,TWO)
      A13 = A(ONE,THREE)
      A23 = A(TWO,THREE)
      A24 = A(TWO,FOUR)
      A25 = A(TWO,FIVE)
      A45 = A(FOUR,FIVE)
      B12 = B(ONE,TWO)
      B13 = B(ONE,THREE)
      B14 = B(ONE,FOUR)
      B16 = B(ONE,SIX)
      B24 = B(TWO,FOUR)
      B36 = B(THREE,SIX)
      S56 = S(FIVE,SIX)
      T123 = S(ONE,TWO)+S(ONE,THREE)+S(TWO,THREE)
      T124 = S(ONE,TWO)+S(ONE,FOUR)+S(TWO,FOUR)

      qqbggFaxslPMPM =
     -   (2*A24*A45*(-(A12*B16) + A23*B36)*
     -      (1/(24.*mtsq) - L1(-T123,-S56)/(2.*S56)))/
     -    (A13*A23*S56) + (2*B13*(A25*B12 + A45*B14)*B36*
     -      (1/(24.*mtsq) - L1(-T124,-S56)/(2.*S56)))/
     -    (B14*B24*S56)

      END                                                   

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTIONqqbggFaxslPMPP(ONE,TWO,THREE,FOUR,FIVE,
     >                                          SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      REAL(KIND(1D0)) S(11,11),XMT,mtsq
      COMMON /DOTPRODUCTS/ S
     >       /TOPMASS/ XMT,mtsq
      REAL(KIND(1D0)) S56, T123, T124
      COMPLEX(KIND(1D0)) A12, A13, A14, A23, A24, A25
     >          ,B13, B14, B34, B36, B46
     >          ,L1

      A12 = A(ONE,TWO)
      A13 = A(ONE,THREE)
      A14 = A(ONE,FOUR)
      A23 = A(TWO,THREE)
      A24 = A(TWO,FOUR)
      A25 = A(TWO,FIVE)
      B13 = B(ONE,THREE)
      B14 = B(ONE,FOUR)
      B34 = B(THREE,FOUR)
      B36 = B(THREE,SIX)
      B46 = B(FOUR,SIX)
      S56 = S(FIVE,SIX)
      T123 = S(ONE,TWO)+S(ONE,THREE)+S(TWO,THREE)
      T124 = S(ONE,TWO)+S(ONE,FOUR)+S(TWO,FOUR)

c No top quark...AK
      qqbggFaxslPMPP =
     -   (2*A25*(-(A12*B14) + A23*B34)*B46*
     -      (1/(24.*mtsq) - L1(-T123,-S56)/(2.*S56)))/
     -    (A13*A23*S56) + (2*A25*(-(A12*B13) - A24*B34)*B36*
     -      (1/(24.*mtsq) - L1(-T124,-S56)/(2.*S56)))/
     -    (A14*A24*S56)

      END                                                   

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTIONqqbggFaxslPMMP(ONE,TWO,THREE,FOUR,FIVE,
     >                                          SIX,A,B) 
      IMPLICIT NONE               

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qqbggFaxslPMPM
      
      qqbggFaxslPMMP = -qqbggFaxslPMPM(TWO,ONE,THREE,FOUR,SIX,FIVE,B,A)
      
      END              
      
************************************************************************

      COMPLEX(KIND(1D0)) FUNCTIONqqbggFaxslPMMM(ONE,TWO,THREE,FOUR,FIVE,
     >                                          SIX,A,B) 
      IMPLICIT NONE               

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qqbggFaxslPMPP
      
      qqbggFaxslPMMM = -qqbggFaxslPMPP(TWO,ONE,THREE,FOUR,SIX,FIVE,B,A)
      
      END              
      
************************************************************************
