************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qggqbA60PMPM(ONE,TWO,THREE,FOUR,FIVE,
     >                                         SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      REAL(KIND(1D0)) S(11,11)
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMMON /DOTPRODUCTS/ S
      REAL(KIND(1D0)) S23,S56,T123,T234
      COMPLEX(KIND(1D0)) A12, A23, A24, A25, A34, A45
     >          ,B12, B13, B16, B23, B34, B36

      A12 = A(ONE,TWO)
      A23 = A(TWO,THREE)
      A24 = A(TWO,FOUR)
      A25 = A(TWO,FIVE)
      A34 = A(THREE,FOUR)
      A45 = A(FOUR,FIVE)
      B12 = B(ONE,TWO)
      B13 = B(ONE,THREE)
      B16 = B(ONE,SIX)
      B23 = B(TWO,THREE)
      B34 = B(THREE,FOUR)
      B36 = B(THREE,SIX)
      S23 = S(TWO,THREE)
      S56 = S(FIVE,SIX)
      T123 = S(ONE,TWO)+S(ONE,THREE)+S(TWO,THREE)
      T234 = S(TWO,THREE)+S(TWO,FOUR)+S(THREE,FOUR)

      qggqbA60PMPM = 
     -   (0,1)*((A24*A45*B13*B16)/(A34*B12*S23*S56) - 
     -     (A45*B13**2*(-(A12*B16) + A23*B36))/(B12*S23*S56*T123) + 
     -     (A24**2*B16*(-(A25*B23) + A45*B34))/(A34*S23*S56*T234))

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qggqbA60PPMM(ONE,TWO,THREE,FOUR,FIVE,
     >                                         SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      REAL(KIND(1D0)) S(11,11)
      COMMON /DOTPRODUCTS/ S
      REAL(KIND(1D0)) S23,S56,T123,T234
      COMPLEX(KIND(1D0)) A12, A13, A23, A34, A35, A45
     >          ,B12, B16, B23, B24, B26, B34

      A12 = A(ONE,TWO)
      A13 = A(ONE,THREE)
      A23 = A(TWO,THREE)
      A34 = A(THREE,FOUR)
      A35 = A(THREE,FIVE)
      A45 = A(FOUR,FIVE)
      B12 = B(ONE,TWO)
      B16 = B(ONE,SIX)
      B23 = B(TWO,THREE)
      B24 = B(TWO,FOUR)
      B26 = B(TWO,SIX)
      B34 = B(THREE,FOUR)
      S23 = S(TWO,THREE)
      S56 = S(FIVE,SIX)
      T123 = S(ONE,TWO)+S(ONE,THREE)+S(TWO,THREE)
      T234 = S(TWO,THREE)+S(TWO,FOUR)+S(THREE,FOUR)

      qggqbA60PPMM = 
     -   (0,1)*(-(((A35*B23 + A45*B24)*(-(A13*B16) - A23*B26))/
     -        (A12*B34*S23*S56)) - 
     -     (A13*A45*B12*(-(A13*B16) - A23*B26))/(A12*S23*S56*T123) + 
     -     (A34*B16*B24*(A35*B23 + A45*B24))/(B34*S23*S56*T234))

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qggqbA60PPPM(ONE,TWO,THREE,FOUR,FIVE,
     >                                         SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) A12, A23, A34, A45, A56

      A12 = A(ONE,TWO)
      A23 = A(TWO,THREE)
      A34 = A(THREE,FOUR)
      A45 = A(FOUR,FIVE)
      A56 = A(FIVE,SIX)

      qggqbA60PPPM = ((0,-1)*A45**2)/(A12*A23*A34*A56)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qgqbgA60PPMM(ONE,TWO,THREE,FOUR,FIVE,
     >                                         SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      REAL(KIND(1D0)) S(11,11)
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMMON /DOTPRODUCTS/ S
      REAL(KIND(1D0)) S56,T124,T234
      COMPLEX(KIND(1D0)) A12, A13, A14, A23, A24, A34, A35, A45
     >          ,B12, B13, B14, B16, B23, B24, B26, B34

      A12 = A(ONE,TWO)
      A13 = A(ONE,THREE)
      A14 = A(ONE,FOUR)
      A23 = A(TWO,THREE)
      A24 = A(TWO,FOUR)
      A34 = A(THREE,FOUR)
      A35 = A(THREE,FIVE)
      A45 = A(FOUR,FIVE)
      B12 = B(ONE,TWO)
      B13 = B(ONE,THREE)
      B14 = B(ONE,FOUR)
      B16 = B(ONE,SIX)
      B23 = B(TWO,THREE)
      B24 = B(TWO,FOUR)
      B26 = B(TWO,SIX)
      B34 = B(THREE,FOUR)
      S56 = S(FIVE,SIX)
      T124 = S(ONE,TWO)+S(ONE,FOUR)+S(TWO,FOUR)
      T234 = S(TWO,THREE)+S(TWO,FOUR)+S(THREE,FOUR)

      qgqbgA60PPMM = 
     -         (0,1)*(((A35*B13 + A45*B14)*(-(A13*B16) - A23*B26))/
     -      (A12*A23*B14*B34*S56) + 
     -     (A35*B12*(-(A14*B16) - A24*B26))/(A12*B14*S56*T124) - 
     -     (A34*B16*(A35*B23 + A45*B24))/(A23*B34*S56*T234))

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qgqbgA60PPMP(ONE,TWO,THREE,FOUR,FIVE,
     >                                         SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) A12, A13, A14, A23, A34, A35, A56

      A12 = A(ONE,TWO)
      A13 = A(ONE,THREE)
      A14 = A(ONE,FOUR)
      A23 = A(TWO,THREE)
      A34 = A(THREE,FOUR)
      A35 = A(THREE,FIVE)
      A56 = A(FIVE,SIX)

      qgqbgA60PPMP = ((0,-1)*A13*A35**2)/(A12*A14*A23*A34*A56)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbggA60PMMP(ONE,TWO,THREE,FOUR,FIVE,
     >                                         SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      REAL(KIND(1D0)) S(11,11)
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMMON /DOTPRODUCTS/ S
      REAL(KIND(1D0)) S34, S56, T134, T234
      COMPLEX(KIND(1D0)) A13, A14, A23, A25, A34, A35
     >          ,B14, B16, B23, B24, B34, B46

      A13 = A(ONE,THREE)
      A14 = A(ONE,FOUR)
      A23 = A(TWO,THREE)
      A25 = A(TWO,FIVE)
      A34 = A(THREE,FOUR)
      A35 = A(THREE,FIVE)
      B14 = B(ONE,FOUR)
      B16 = B(ONE,SIX)
      B23 = B(TWO,THREE)
      B24 = B(TWO,FOUR)
      B34 = B(THREE,FOUR)
      B46 = B(FOUR,SIX)
      S34 = S(THREE,FOUR)
      S56 = S(FIVE,SIX)
      T134 = S(ONE,THREE)+S(ONE,FOUR)+S(THREE,FOUR)
      T234 = S(TWO,THREE)+S(TWO,FOUR)+S(THREE,FOUR)

      qqbggA60PMMP =
     -   (0,1)*(((-(A25*B24) - A35*B34)*(-(A13*B16) + A34*B46))/
     -      (A14*B23*S34*S56) - 
     -     (A13*A25*B14*(-(A13*B16) + A34*B46))/(A14*S34*S56*T134) - 
     -     (A23*B16*B24*(-(A25*B24) - A35*B34))/(B23*S34*S56*T234))

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbggA60PMPM(ONE,TWO,THREE,FOUR,FIVE,
     >                                         SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      REAL(KIND(1D0)) S(11,11)
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMMON /DOTPRODUCTS/ S
      REAL(KIND(1D0)) S34, S56, T134, T234
      COMPLEX(KIND(1D0)) A14, A23, A24, A25, A34, A45
     >          ,B13, B14, B16, B23, B34, B36

      A14 = A(ONE,FOUR)
      A23 = A(TWO,THREE)
      A24 = A(TWO,FOUR)
      A25 = A(TWO,FIVE)
      A34 = A(THREE,FOUR)
      A45 = A(FOUR,FIVE)
      B13 = B(ONE,THREE)
      B14 = B(ONE,FOUR)
      B16 = B(ONE,SIX)
      B23 = B(TWO,THREE)
      B34 = B(THREE,FOUR)
      B36 = B(THREE,SIX)
      S34 = S(THREE,FOUR)
      S56 = S(FIVE,SIX)
      T134 = S(ONE,THREE)+S(ONE,FOUR)+S(THREE,FOUR)
      T234 = S(TWO,THREE)+S(TWO,FOUR)+S(THREE,FOUR)

      qqbggA60PMPM =
     -   (0,-1)*(-((A24*A25*B13*B16)/(A23*B14*S34*S56)) + 
     -     (A25*B13**2*(-(A14*B16) - A34*B36))/(B14*S34*S56*T134) + 
     -     (A24**2*B16*(-(A25*B23) + A45*B34))/(A23*S34*S56*T234))

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbggA60PMPP(ONE,TWO,THREE,FOUR,FIVE,
     >                                         SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      REAL(KIND(1D0)) S(11,11)
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMMON /DOTPRODUCTS/ S
      COMPLEX(KIND(1D0)) A14, A23, A25, A34, A56

      A14 = A(ONE,FOUR)
      A23 = A(TWO,THREE)
      A25 = A(TWO,FIVE)
      A34 = A(THREE,FOUR)
      A56 = A(FIVE,SIX)

      qqbggA60PMPP = ((0,-1)*A25**2)/(A14*A23*A34*A56)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qggqbA60PMMM(ONE,TWO,THREE,FOUR,FIVE,
     >                                         SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qggqbA60PPPM

      qggqbA60PMMM = qggqbA60PPPM(FOUR,THREE,TWO,ONE,SIX,FIVE,B,A)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qgqbgA60PMMP(ONE,TWO,THREE,FOUR,FIVE,
     >                                         SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qgqbgA60PPMM

      qgqbgA60PMMP = qgqbgA60PPMM(THREE,TWO,ONE,FOUR,SIX,FIVE,B,A)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qgqbgA60PMMM(ONE,TWO,THREE,FOUR,FIVE,
     >                                         SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qgqbgA60PPMP

      qgqbgA60PMMM = qgqbgA60PPMP(THREE,TWO,ONE,FOUR,SIX,FIVE,B,A)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbggA60PMMM(ONE,TWO,THREE,FOUR,FIVE,
     >                                         SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qqbggA60PMPP

      qqbggA60PMMM = qqbggA60PMPP(TWO,ONE,FOUR,THREE,SIX,FIVE,B,A)

      END

************************************************************************
