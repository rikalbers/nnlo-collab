************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qggqbA60PMPMNT
     >                                (ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
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

      qggqbA60PMPMNT = 
     >   (B13**2 * A45 * (-A12 * B16 + A23 * B36))/(B12 * S23 * T123)
     >  -(A24**2 * B16 * (-A25 * B23 + A45 * B34))/(A34 * S23 * T234)
     >  -(B13 * A24 * B16 * A45)/(B12 * A34 * S23)

      qggqbA60PMPMNT = (0,1)/S56 * qggqbA60PMPMNT

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qggqbA60PPMMNT
     >                                (ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
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

      qggqbA60PPMMNT = 
     >   (A13 * B12 * A45 * (-A13 * B16 - A23 * B26))/(A12 * S23 * T123)
     >  -(A34 * B24 * B16 * (A35 * B23 + A45 * B24))/(B34 * s23 * T234)
     >  +(A35 * B23 + A45 * B24) * (-A13 * B16 - A23 * B26)
     >  /(A12 * B34 * S23)

      qggqbA60PPMMNT = (0,1)/S56 * qggqbA60PPMMNT 

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qggqbA60PPPMNT
     >                                (ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      REAL(KIND(1D0)) S(11,11)
      COMMON /DOTPRODUCTS/ S
      REAL(KIND(1D0)) S56
      COMPLEX(KIND(1D0)) A12, A23, A34, A45, B56

      A12 = A(ONE,TWO)
      A23 = A(TWO,THREE)
      A34 = A(THREE,FOUR)
      A45 = A(FOUR,FIVE)
      B56 = B(FIVE,SIX)
      S56 = S(FIVE,SIX)

      qggqbA60PPPMNT = -(A45**2 * B56)/(A12 * A23 * A34)
      
      qggqbA60PPPMNT = (0,1)/S56 * qggqbA60PPPMNT

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qggqbA60PMMMNT
     >                                (ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      REAL(KIND(1D0)) S(11,11)
      COMMON /DOTPRODUCTS/ S
      REAL(KIND(1D0)) S56
      COMPLEX(KIND(1D0)) B12, B23, B34, B16, A56

      B12 = B(ONE,TWO)
      B23 = B(TWO,THREE)
      B34 = B(THREE,FOUR)
      B16 = B(ONE,SIX)
      A56 = A(FIVE,SIX)
      S56 = S(FIVE,SIX)

      qggqbA60PMMMNT = -(B16**2 * A56)/(B12 * B23 * B34)
      
      qggqbA60PMMMNT = (0,1)/S56 * qggqbA60PMMMNT

      END

************************************************************************
