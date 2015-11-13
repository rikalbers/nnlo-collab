************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbggFvsPMPM(ONE,TWO,THREE,FOUR,FIVE,
     >                                         SIX,A,B)
      IMPLICIT NONE 

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      REAL(KIND(1D0)) S(11,11)
      COMMON /DOTPRODUCTS/ S
      REAL(KIND(1D0)) S12, S34, S56, T123, T124, T356, T456, I3m3
      COMPLEX(KIND(1D0)) A12,A13,A14,A15,A23,A24,A25,A26,A34,A35,A36,A45
     >          ,B12,B13,B14,B15,B16,B23,B24,B26,B34,B35,B36,B45,B46,B56
     >          ,L0, L1, Ln, Lsm1h,A46,A56

      A12 = A(ONE,TWO)
      A13 = A(ONE,THREE)
      A14 = A(ONE,FOUR)
      A15 = A(ONE,FIVE)
      A23 = A(TWO,THREE)
      A24 = A(TWO,FOUR)
      A25 = A(TWO,FIVE)
      A26 = A(TWO,SIX)
      A34 = A(THREE,FOUR)
      A35 = A(THREE,FIVE)
      A36 = A(THREE,SIX)
      A45 = A(FOUR,FIVE)
      A46 = A(FOUR,SIX)
      A56 = A(FIVE,SIX)
      B12 = B(ONE,TWO)
      B13 = B(ONE,THREE)
      B14 = B(ONE,FOUR)
      B15 = B(ONE,FIVE)
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

      qqbggFvsPMPM = 
     >  ((A25*B12+A35*B13)**2-A25*A56*B12*B16)/
     >   (A56*B12*(-(A13*B14)-A23*B24)**2)+
     >  ((-(A12*B16)+A24*B46)**2-A12*A25*B16*B56)/
     >   (A12*(-(A13*B14)-A23*B24)**2*B56)+
     >  (A25*(-(A14*B13)-A24*B23)*
     >     (-2*A12*A56*B16+A25*(-S12+S34-S56)))/
     >   (A12*A56*(-(A13*B14)-A23*B24)*
     >     (S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+S56**2))+
     >  (B16*(-(A14*B13)-A24*B23)*
     >     (-2*A25*B12*B56+B16*(-S12+S34-S56)))/
     >   (B12*(-(A13*B14)-A23*B24)*B56*
     >     (S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+S56**2))+
     >  (-((A24*A45*B13*B36)/
     >        (S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+S56**2))
     >     -(A23*A35*B14*(-(A14*B13)-A24*B23)**2*B46)/
     >      ((-(A13*B14)-A23*B24)**2*
     >        (S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+S56**2))
     >     -(3*(-(A14*B13)-A24*B23)*
     >        (-(A24*A35*B14*B36)-A23*A45*B13*B46+
     >          ((-(A23*B13)-A24*B14)*(-(A35*B36)-A45*B46)*S34*
     >             (-S12+S34-S56))/
     >           (S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+
     >             S56**2)))/
     >      ((-(A13*B14)-A23*B24)*
     >        (S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+S56**2))
     >     +(2*A23*(-(A14*B13)-A24*B23)*B46*
     >        (-((A25*B12+A45*B14)*(-(A13*B14)-A23*B24))+
     >          A23*A56*B12*B46)*(T123-T124))/
     >      ((-(A13*B14)-A23*B24)**3*
     >        (S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+S56**2)))
     >    *I3m3(S12,S34,S56)+(-((A24*A45*B13*B36)/
     >        (S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+S56**2))
     >     -(A23*A35*B14*(-(A14*B13)-A24*B23)**2*B46)/
     >      ((-(A13*B14)-A23*B24)**2*
     >        (S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+S56**2))
     >     -(3*(-(A14*B13)-A24*B23)*
     >        (-(A24*A35*B14*B36)-A23*A45*B13*B46+
     >          ((-(A23*B13)-A24*B14)*(-(A35*B36)-A45*B46)*S34*
     >             (-S12+S34-S56))/
     >           (S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+
     >             S56**2)))/
     >      ((-(A13*B14)-A23*B24)*
     >        (S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+S56**2))
     >     +(2*A35*B14*(-(A14*B13)-A24*B23)*
     >        (-((-(A13*B14)-A23*B24)*(-(A12*B16)+A23*B36))+
     >          A12*A35*B14*B56)*(-T123+T124))/
     >      ((-(A13*B14)-A23*B24)**3*
     >        (S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+S56**2)))
     >    *I3m3(S12,S34,S56)+(2*A23*B12*B36*
     >     (((-(A13*B16)-A23*B26)*(-(A12*B14)+A23*B34)*L0(-T123,-S12))/
     >        ((-(A13*B14)-A23*B24)*S12)-
     >       A23*B36*(L0(-T123,-S12)/S12-L1(-T123,-S12)/(2.*S12))))/
     >   ((-(A13*B14)-A23*B24)**2*B56)+
     >  (2*A12*A45*B14*(((A23*B12-A34*B14)*(-(A15*B14)-A25*B24)*
     >          L0(-T124,-S12))/((-(A13*B14)-A23*B24)*S12)-
     >       A45*B14*(L0(-T124,-S12)/S12-L1(-T124,-S12)/(2.*S12))))/
     >   (A56*(-(A13*B14)-A23*B24)**2)+
     >  (2*A35*B13*B56*(((-(A35*B15)-A36*B16)*(-(A35*B34)-A56*B46)*
     >          L0(-T356,-S56))/((-(A35*B45)-A36*B46)*S56)-
     >       A35*B13*(L0(-T356,-S56)/S56-L1(-T356,-S56)/(2.*S56))))/
     >   (B12*(-(A35*B45)-A36*B46)**2)+
     >  (2*A24*A56*B46*(((-(A25*B45)-A26*B46)*(A34*B46+A35*B56)*
     >          L0(-T456,-S56))/((-(A35*B45)-A36*B46)*S56)-
     >       A24*B46*(L0(-T456,-S56)/S56-L1(-T456,-S56)/(2.*S56))))/
     >   (A12*(-(A35*B45)-A36*B46)**2)-
     >  ((3*(-(A23*B13)-A24*B14)*(-(A14*B13)-A24*B23)*
     >        (-(A35*B36)-A45*B46)*(-S12-S34+S56))/
     >      ((-(A13*B14)-A23*B24)*
     >        (S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+S56**2)**
     >         2)+((-(A14*B13)-A24*B23)*
     >        ((A25**2*B12)/A56-2*A25*B16+(A12*B16**2)/B56))/
     >      (2.*(-(A13*B14)-A23*B24)*
     >        (S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+S56**2))
     >     -(A25*B12*(((S12-S34-S56)*(A12*A56*B16+A25*T123))/A56+
     >          2*A23*B36*(T123-T124)))/
     >      ((-(A13*B14)-A23*B24)**2*
     >        (S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+S56**2))
     >     +(2*A23*B12*B46*(A12*A56*B16+A25*T123)*(T123-T124))/
     >      ((-(A13*B14)-A23*B24)**3*
     >        (S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+S56**2)))
     >    *Ln(-S12,-S34)-((3*(-(A23*B13)-A24*B14)*(-(A14*B13)-A24*B23)*
     >        (-(A35*B36)-A45*B46)*(-S12-S34+S56))/
     >      ((-(A13*B14)-A23*B24)*
     >        (S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+S56**2)**
     >         2)+((-(A14*B13)-A24*B23)*
     >        ((A25**2*B12)/A56-2*A25*B16+(A12*B16**2)/B56))/
     >      (2.*(-(A13*B14)-A23*B24)*
     >        (S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+S56**2))
     >     -(2*A12*A35*B14*(-T123+T124)*(A25*B12*B56+B16*T124))/
     >      ((-(A13*B14)-A23*B24)**3*
     >        (S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+S56**2))
     >     +(A12*B16*(2*A45*B14*(-T123+T124)-
     >          ((S12-S34-S56)*(A25*B12*B56+B16*T124))/B56))/
     >      ((-(A13*B14)-A23*B24)**2*
     >        (S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+S56**2)))
     >    *Ln(-S12,-S34)-((3*(-(A23*B13)-A24*B14)*(-(A45*B35)-A46*B36)*
     >        (-(A35*B36)-A45*B46)*(S12-S34-S56))/
     >      ((-(A35*B45)-A36*B46)*
     >        (S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+S56**2)**
     >         2)+((-(A45*B35)-A46*B36)*
     >        (-2*A25*B16+(A56*B16**2)/B12+(A25**2*B56)/A12))/
     >      (2.*(-(A35*B45)-A36*B46)*
     >        (S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+S56**2))
     >     -(A25*B56*(-(((-S12-S34+S56)*(-(A12*A56*B16)-A25*T356))/
     >             A12)+2*A35*B13*(T356-T456)))/
     >      ((-(A35*B45)-A36*B46)**2*
     >        (S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+S56**2))
     >     -(2*A35*B14*B56*(-(A12*A56*B16)-A25*T356)*(T356-T456))/
     >      ((-(A35*B45)-A36*B46)**3*
     >        (S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+S56**2)))
     >    *Ln(-S56,-S34)-((3*(-(A23*B13)-A24*B14)*(-(A45*B35)-A46*B36)*
     >        (-(A35*B36)-A45*B46)*(S12-S34-S56))/
     >      ((-(A35*B45)-A36*B46)*
     >        (S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+S56**2)**
     >         2)+((-(A45*B35)-A46*B36)*
     >        (-2*A25*B16+(A56*B16**2)/B12+(A25**2*B56)/A12))/
     >      (2.*(-(A35*B45)-A36*B46)*
     >        (S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+S56**2))
     >     +(2*A23*A56*B46*(-T356+T456)*(-(A25*B12*B56)-B16*T456))/
     >      ((-(A35*B45)-A36*B46)**3*
     >        (S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+S56**2))
     >     +(A56*B16*(2*A24*B46*(-T356+T456)+
     >          ((-S12-S34+S56)*(-(A25*B12*B56)-B16*T456))/B12))/
     >      ((-(A35*B45)-A36*B46)**2*
     >        (S12**2-2*S12*S34+S34**2-2*S12*S56-2*S34*S56+S56**2)))
     >    *Ln(-S56,-S34)+(-((2*A23*A24*B36*B46+(-(A12*B16)+A24*B46)**2)/
     >        (A12*(-(A13*B14)-A23*B24)**2*B56))+
     >     (2*A23*B46*(-(A12*B16)+A24*B46)*T123)/
     >      (A12*(-(A13*B14)-A23*B24)**3*B56))*Ln(-T123,-S34)+
     >  (-(((A25*B12+A35*B13)**2+2*A35*A45*B13*B14)/
     >        (A56*B12*(-(A13*B14)-A23*B24)**2))+
     >     (2*A35*(A25*B12+A35*B13)*B14*T124)/
     >      (A56*B12*(-(A13*B14)-A23*B24)**3))*Ln(-T124,-S34)-
     >  (2*A23*(-(A13*B16)-A23*B26)*B46*(-(A25*B45)-A26*B46)*T123*
     >     Lsm1h(S34,T123,S12,S56))/(A12*(-(A13*B14)-A23*B24)**4*B56)-
     >  (2*A35*B14*(-(A35*B15)-A36*B16)*(-(A15*B14)-A25*B24)*T124*
     >     Lsm1h(S34,T124,S12,S56))/(A56*B12*(-(A13*B14)-A23*B24)**4)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbggFvsPMPP(ONE,TWO,THREE,FOUR,FIVE,
     >                                         SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      REAL(KIND(1D0)) S(11,11)
      COMMON /DOTPRODUCTS/ S
      REAL(KIND(1D0)) S12, S56, T123, T124, T356, T456
      COMPLEX(KIND(1D0)) A12, A23, A24, A25, A34, A35, A45, A56
     >          ,B12, B13, B14, B16, B36, B46, B56
     >          ,L0, L1, Ln, Lsm1e
      
      A12 = A(ONE,TWO)
      A23 = A(TWO,THREE)
      A24 = A(TWO,FOUR)
      A25 = A(TWO,FIVE)
      A34 = A(THREE,FOUR)
      A35 = A(THREE,FIVE)
      A45 = A(FOUR,FIVE)
      A56 = A(FIVE,SIX)
      B12 = B(ONE,TWO)
      B13 = B(ONE,THREE)
      B14 = B(ONE,FOUR)
      B16 = B(ONE,SIX)
      B36 = B(THREE,SIX)
      B46 = B(FOUR,SIX)
      B56 = B(FIVE,SIX)
      S12 = S(ONE,TWO)
      S56 = S(FIVE,SIX)
      T123 = S(ONE,TWO)+S(ONE,THREE)+S(TWO,THREE)
      T124 = S(ONE,TWO)+S(ONE,FOUR)+S(TWO,FOUR)
      T356 = S(THREE,FIVE)+S(THREE,SIX)+S(FIVE,SIX)
      T456 = S(FOUR,FIVE)+S(FOUR,SIX)+S(FIVE,SIX)

      qqbggFvsPMPP =
     -   (A25**2/(A12*A56) - B16**2/(B12*B56))/A34**2 + 
     -   (A25*A35*B13*L0(-T123,-S12))/(A34**2*A56*S12) + 
     -   (A25*A45*B14*L0(-T124,-S12))/(A34**2*A56*S12) + 
     -   (A23*A25*B36*L0(-T356,-S56))/(A12*A34**2*S56) + 
     -   (A24*A25*B46*L0(-T456,-S56))/(A12*A34**2*S56) - 
     -   (A12*A35**2*B13**2*L1(-T123,-S12))/(A34**2*A56*S12**2) - 
     -   (A12*A45**2*B14**2*L1(-T124,-S12))/(A34**2*A56*S12**2) - 
     -   (A23**2*A56*B36**2*L1(-T356,-S56))/(A12*A34**2*S56**2) - 
     -   (A24**2*A56*B46**2*L1(-T456,-S56))/(A12*A34**2*S56**2) - 
     -   ((A24*A35 + A23*A45)*((A12*A35*B13*L0(-T123,-S12))/S12 + 
     -   (A24*A56*B46*L0(-T123,-S56))/S56 - (A25*Ln(-T123,-T124))/2.))/
     -    (A12*A34**3*A56) + ((A24*A35 + A23*A45)*
     -      ((A12*A45*B14*L0(-T124,-S12))/S12 + 
     -   (A23*A56*B36*L0(-T124,-S56))/S56 - (A25*Ln(-T124,-T123))/2.))/
     -    (A12*A34**3*A56) - (A23*A24*A35*A45*Lsm1e(T123,T124,S12,S56))/
     -    (A12*A34**4*A56) - (A23*A24*A35*A45*Lsm1e(T124,T123,S12,S56))/
     -    (A12*A34**4*A56)

      END                                                   

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbggFvsPMMP(ONE,TWO,THREE,FOUR,FIVE,
     >                                         SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qqbggFvsPMPM

      qqbggFvsPMMP = qqbggFvsPMPM(TWO,ONE,THREE,FOUR,SIX,FIVE,B,A)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbggFvsPMMM(ONE,TWO,THREE,FOUR,FIVE,
     >                                         SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qqbggFvsPMPP

      qqbggFvsPMMM = qqbggFvsPMPP(TWO,ONE,THREE,FOUR,SIX,FIVE,B,A)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbggFvfPMPM(ONE,TWO,THREE,FOUR,FIVE,
     >                                         SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      REAL(KIND(1D0)) S(11,11)
      COMMON /DOTPRODUCTS/ S
      REAL(KIND(1D0)) S12, S34, S56, T123, T124, T356, T456, I3m3
      COMPLEX(KIND(1D0)) A12, A13, A23, A24, A25, A35, A36, A45, A56
     >          ,B12, B13, B14, B16, B24, B36, B45, B46, B56
     >          ,Lsm1h

      A12 = A(ONE,TWO)
      A13 = A(ONE,THREE)
      A23 = A(TWO,THREE)
      A24 = A(TWO,FOUR)
      A25 = A(TWO,FIVE)
      A35 = A(THREE,FIVE)
      A36 = A(THREE,SIX)
      A45 = A(FOUR,FIVE)
      A56 = A(FIVE,SIX)
      B12 = B(ONE,TWO)
      B13 = B(ONE,THREE)
      B14 = B(ONE,FOUR)
      B16 = B(ONE,SIX)
      B24 = B(TWO,FOUR)
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

      qqbggFvfPMPM =
     -   -(A45*B13*(-(A12*B16) + A23*B36)*I3m3(S12,S34,S56))/
     -    (2.*(-(A13*B14) - A23*B24)*T123) - 
     -   (A24*(A25*B12 + A45*B14)*B36*I3m3(S12,S34,S56))/
     -    (2.*(-(A13*B14) - A23*B24)*T124) - 
     -   (A24*(A35*B13 - A56*B16)*B36*I3m3(S56,S34,S12))/
     -    (2.*(-(A35*B45) - A36*B46)*T356) - 
     -   (A45*B13*(A24*B46 + A25*B56)*I3m3(S56,S34,S12))/
     -    (2.*(-(A35*B45) - A36*B46)*T456) - 
     -   ((-(A12*B16) + A23*B36)**2*Lsm1h(S34,T123,S12,S56))/
     -    (A12*(-(A13*B14) - A23*B24)**2*B56) - 
     -   ((A25*B12 + A45*B14)**2*Lsm1h(S34,T124,S12,S56))/
     -    (A56*B12*(-(A13*B14) - A23*B24)**2)

      END                                                   

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbggFvfPMPP(ONE,TWO,THREE,FOUR,FIVE,
     >                                         SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      REAL(KIND(1D0)) S(11,11)
      COMMON /DOTPRODUCTS/ S
      REAL(KIND(1D0)) T123,T124,S12,S56
      COMPLEX(KIND(1D0)) A12, A25, A34, A56, Lsm1e

      A12 = A(ONE,TWO)
      A25 = A(TWO,FIVE)
      A34 = A(THREE,FOUR)
      A56 = A(FIVE,SIX)
      S12 = S(ONE,TWO)
      S56 = S(FIVE,SIX)
      T123 = S(ONE,TWO)+S(ONE,THREE)+S(TWO,THREE)
      T124 = S(ONE,TWO)+S(ONE,FOUR)+S(TWO,FOUR)

      qqbggFvfPMPP =
     >  -((A25**2*Lsm1e(T123,T124,S12,S56))/(A12*A34**2*A56))

      END                                                   

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbggFvfPMMP(ONE,TWO,THREE,FOUR,FIVE,
     >                                         SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qqbggFvfPMPM

      qqbggFvfPMMP = qqbggFvfPMPM(TWO,ONE,THREE,FOUR,SIX,FIVE,B,A)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbggFvfPMMM(ONE,TWO,THREE,FOUR,FIVE,
     >                                         SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qqbggFvfPMPP

      qqbggFvfPMMM = qqbggFvfPMPP(TWO,ONE,THREE,FOUR,SIX,FIVE,B,A)

      END

************************************************************************
