************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qggqbA6sPMPM(ONE,TWO,THREE,FOUR,FIVE,
     >                                         SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)

      qggqbA6sPMPM = (0,0)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qggqbA6sPPMM(ONE,TWO,THREE,FOUR,FIVE,
     >                                         SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)

      qggqbA6sPPMM = (0,0)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qggqbA6sPPPM(ONE,TWO,THREE,FOUR,FIVE,
     >                                         SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      REAL(KIND(1D0)) S(11,11)
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMMON /DOTPRODUCTS/ S
      REAL(KIND(1D0)) S56,T123,T234
      COMPLEX(KIND(1D0)) A13, A23, A25, A34, A45
     >          ,B13, B16, B23, B26, B34

      A13 = A(ONE,THREE)
      A23 = A(TWO,THREE)
      A25 = A(TWO,FIVE)
      A34 = A(THREE,FOUR)
      A45 = A(FOUR,FIVE)
      B13 = B(ONE,THREE)
      B16 = B(ONE,SIX)
      B23 = B(TWO,THREE)
      B26 = B(TWO,SIX)
      B34 = B(THREE,FOUR)
      S56 = S(FIVE,SIX)
      T123 = S(ONE,TWO)+S(ONE,THREE)+S(TWO,THREE)
      T234 = S(TWO,THREE)+S(TWO,FOUR)+S(THREE,FOUR)

      qggqbA6sPPPM = (0,1)*
     > (-A45*(B16*A13+B26*A23)*B13/T123 
     >  +B16*(A45*B34-A25*B23)*A34/T234)/(3.*A23**2*S56)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qggqbA6sPMMM(ONE,TWO,THREE,FOUR,FIVE,
     >                            SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qggqbA6sPPPM

      qggqbA6sPMMM = qggqbA6sPPPM(FOUR,THREE,TWO,ONE,SIX,FIVE,B,A)

      END

************************************************************************
************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qggqbA6tPMPM(ONE,TWO,THREE,FOUR,FIVE,
     >                                         SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)

      qggqbA6tPMPM = (0,0)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qggqbA6tPPMM(ONE,TWO,THREE,FOUR,FIVE,
     >                                         SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)

      qggqbA6tPPMM = (0,0)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qggqbA6tPPPM(ONE,TWO,THREE,FOUR,FIVE,
     >                                         SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      REAL(KIND(1D0)) S(11,11),XMT,mtsq
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMMON /DOTPRODUCTS/ S
     >       /TOPMASS/ XMT,mtsq
      REAL(KIND(1D0)) S23,S56,T123,T234
      COMPLEX(KIND(1D0)) A13, A23, A25, A34, A45
     >          ,B13, B16, B23, B26, B34

      A13 = A(ONE,THREE)
      A23 = A(TWO,THREE)
      A25 = A(TWO,FIVE)
      A34 = A(THREE,FOUR)
      A45 = A(FOUR,FIVE)
      B13 = B(ONE,THREE)
      B16 = B(ONE,SIX)
      B23 = B(TWO,THREE)
      B26 = B(TWO,SIX)
      B34 = B(THREE,FOUR)
      S23 = S(TWO,THREE)
      S56 = S(FIVE,SIX)
      T123 = S(ONE,TWO)+S(ONE,THREE)+S(TWO,THREE)
      T234 = S(TWO,THREE)+S(TWO,FOUR)+S(THREE,FOUR)

      qggqbA6tPPPM = 0d0

      qggqbA6tPPPM = (0,1)*
     > (-A45*(B16*A13+B26*A23)*B13/T123 
     >  +B16*(A45*B34-A25*B23)*A34/T234)/(3.*A23**2*S56)
      qggqbA6tPPPM = qggqbA6tPPPM*S23/(20.*mtsq)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qggqbA6tPMMM(ONE,TWO,THREE,FOUR,FIVE,
     >                                         SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qggqbA6tPPPM

      qggqbA6tPMMM = qggqbA6tPPPM(FOUR,THREE,TWO,ONE,SIX,FIVE,B,A)

      END

************************************************************************
