************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qggqbA70PPPPMNT
     >                          (ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN
      REAL(KIND(1D0)) S(11,11)
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMMON /DOTPRODUCTS/ S
      REAL(KIND(1D0)) S67
      COMPLEX(KIND(1D0)) A12, A23, A34, A45, A56
     >          ,B67

      A12 = A(ONE,TWO)
      A23 = A(TWO,THREE)
      A34 = A(THREE,FOUR)
      A45 = A(FOUR,FIVE)
      A56 = A(FIVE,SIX)
      B67 = B(SIX,SEVEN)
      S67 = S(SIX,SEVEN)

      qggqbA70PPPPMNT = -(A56**2 * B67)/(A12 * A23 * A34 * A45) 

      qggqbA70PPPPMNT = (0,1)/S67 * qggqbA70PPPPMNT

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qggqbA70PPPMMNT
     >                          (ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN
      REAL(KIND(1D0)) S(11,11)
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMMON /DOTPRODUCTS/ S
      REAL(KIND(1D0)) S34,S67,T167,T234,T345,T567
      COMPLEX(KIND(1D0)) A12,A23,A24,A34,A45,A56
     >          ,A4123,A4127,A4231,A4567,A6172,A6173,A6543
     >          ,B17,B23,B34,B35,B45

      A12 = A(ONE,TWO)
      A23 = A(TWO,THREE)
      A24 = A(TWO,FOUR)
      A34 = A(THREE,FOUR)
      A45 = A(FOUR,FIVE)
      A56 = A(FIVE,SIX)

      A4123 = A(FOUR,ONE) * B(ONE,THREE) 
     >      + A(FOUR,TWO) * B(TWO,THREE)
      A4127 = A(FOUR,ONE) * B(ONE,SEVEN) 
     >      + A(FOUR,TWO) * B(TWO,SEVEN)
      A4231 = A(FOUR,TWO) * B(TWO,ONE) 
     >      + A(FOUR,THREE) * B(THREE,ONE)
      A4567 = A(FOUR,FIVE) * B(FIVE,SEVEN) 
     >      + A(FOUR,SIX) * B(SIX,SEVEN)
      A6172 = A(SIX,ONE) * B(ONE,TWO) 
     >      + A(SIX,SEVEN) * B(SEVEN,TWO)
      A6173 = A(SIX,ONE) * B(ONE,THREE) 
     >      + A(SIX,SEVEN) * B(SEVEN,THREE)
      A6543 = A(SIX,FIVE) * B(FIVE,THREE) 
     >      + A(SIX,FOUR) * B(FOUR,THREE)

      B17 = B(ONE,SEVEN)
      B23 = B(TWO,THREE)
      B34 = B(THREE,FOUR)
      B35 = B(THREE,FIVE)
      B45 = B(FOUR,FIVE)

      S34 = S(THREE,FOUR)
      S67 = S(SIX,SEVEN)
      
      T167 = S(ONE,SIX) + S(ONE,SEVEN) + S(SIX,SEVEN)
      T234 = S(TWO,THREE) + S(TWO,FOUR) + S(THREE,FOUR)
      T345 = S(THREE,FOUR) + S(THREE,FIVE) + S(FOUR,FIVE)
      T567 = S(FIVE,SIX) + S(FIVE,SEVEN) + S(SIX,SEVEN)

      qggqbA70PPPMMNT = 
     > -(A56 * A4567)/(A23 * A34 * B34 * T567)
     >  *((B23 * A4231)/T234 + A4123/A12)
     > +(A4567 * A6543)/(A12 * A23 * A34 * B34 * B45)
     > +(A4127 * A6543 * A45 * B35)/(A12 * A24 * B45 * S34 * T345)
     > +(B17 * A6172 * A45**2 * B35**2)/(B45 * A24 * S34 * T345 * T167)
     > -(B17 * A45 * B35)/(A23 * B34 * B45 * T167)
     >  *(A6172/A34 + A6173/A24)
     > -(B17 * A45 * B23)/(A23 * B34 * T234 * T167)
     >  *((A6172 * A24)/A34 + A6173)

      qggqbA70PPPMMNT = (0,1)/S67 * qggqbA70PPPMMNT

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qggqbA70PPMPMNT
     >                          (ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN
      REAL(KIND(1D0)) S(11,11)
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMMON /DOTPRODUCTS/ S
      REAL(KIND(1D0)) S23,S34,S67,T123,T167,T234,T345,T567
      COMPLEX(KIND(1D0)) A12,A13,A23,A34,A35,A45,A56
     >          ,A3124,A3127,A3241,A3542,A3567,A6172,A6534,A617243
     >          ,B12,B17,B23,B24,B45

      A12 = A(ONE,TWO)
      A13 = A(ONE,THREE)
      A23 = A(TWO,THREE)
      A34 = A(THREE,FOUR)
      A35 = A(THREE,FIVE)
      A45 = A(FOUR,FIVE)
      A56 = A(FIVE,SIX)

      A3124 = A(THREE,ONE) * B(ONE,FOUR)
     >      + A(THREE,TWO) * B(TWO,FOUR)
      A3127 = A(THREE,ONE) * B(ONE,SEVEN)
     >      + A(THREE,TWO) * B(TWO,SEVEN)
      A3241 = A(THREE,TWO) * B(TWO,ONE)
     >      + A(THREE,FOUR) * B(FOUR,ONE)
      A3542 = A(THREE,FIVE) * B(FIVE,TWO)
     >      + A(THREE,FOUR) * B(FOUR,TWO)
      A3567 = A(THREE,FIVE) * B(FIVE,SEVEN)
     >      + A(THREE,SIX) * B(SIX,SEVEN)      
      A6172 = A(SIX,ONE) * B(ONE,TWO)
     >      + A(SIX,SEVEN) * B(SEVEN,TWO)      
      A6534 = A(SIX,FIVE) * B(FIVE,FOUR)
     >      + A(SIX,THREE) * B(THREE,FOUR)
      A617243 = A(SIX,ONE) * B(ONE,TWO) * A(TWO,THREE)
     >        + A(SIX,ONE) * B(ONE,FOUR) * A(FOUR,THREE)
     >        + A(SIX,SEVEN) * B(SEVEN,TWO) * A(TWO,THREE)
     >        + A(SIX,SEVEN) * B(SEVEN,FOUR) * A(FOUR,THREE)

      B12 = B(ONE,TWO)
      B17 = B(ONE,SEVEN)
      B23 = B(TWO,THREE)
      B24 = B(TWO,FOUR)
      B45 = B(FOUR,FIVE)

      S23 = S(TWO,THREE)
      S34 = S(THREE,FOUR)
      S67 = S(SIX,SEVEN)

      T123 = S(ONE,TWO) + S(ONE,THREE) + S(TWO,THREE)
      T167 = S(ONE,SIX) + S(ONE,SEVEN) + S(SIX,SEVEN)
      T234 = S(TWO,THREE) + S(TWO,FOUR) + S(THREE,FOUR)
      T345 = S(THREE,FOUR) + S(THREE,FIVE) + S(FOUR,FIVE)
      T567 = S(FIVE,SIX) + S(FIVE,SEVEN) + S(SIX,SEVEN)

      qggqbA70PPMPMNT = 
     > (A13 * B12 * A3124 * A3567 * A56)/(A12 * A34 * S23 * T123 * T567)
     > +(A13 * B12 * A3127 * A56 * A35)/(A12 * A34 * A45 * S23 * T123)
     > -(A3127 * A6172 * A35**2)/(A12 * A23 * A34 * A45 * B23 * T345)
     > -(A3127 * A6534 * A35 * B24)/(B23 * A12 * A23 * S34 * T345)
     > -(B24**2 * B17 * A617243 * A35)/(S23 * S34 * T234 * T167)
     > +(B24 * A3567 * A56)/(S23 * S34 * T567)
     >  *(-(B24 * A3241)/T234 - A3124/A12)
     > +(B17 * A6172 * A35**2)/(S23 * T345 * T167)
     >  *(A3542/(A34 * A45) + (B24 * B45)/S34)

      qggqbA70PPMPMNT = (0,1)/S67 * qggqbA70PPMPMNT

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qggqbA70PMPPMNT
     >                          (ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN
      REAL(KIND(1D0)) S(11,11)
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMMON /DOTPRODUCTS/ S
      REAL(KIND(1D0)) S23,S34,S67,T123,T167,T234,T345,T567
      COMPLEX(KIND(1D0)) A23,A24,A25,A34,A45,A56
     >          ,A2134,A2137,A2341,A2543,A2534,A2567,A6173,A6534,A6543
     >          ,A617342
     >          ,B12,B13,B17,B34

      A23 = A(TWO,THREE)
      A24 = A(TWO,FOUR)
      A25 = A(TWO,FIVE)
      A34 = A(THREE,FOUR)
      A45 = A(FOUR,FIVE)
      A56 = A(FIVE,SIX)

      A2134 = A(TWO,ONE) * B(ONE,FOUR)
     >      + A(TWO,THREE) * B(THREE,FOUR)
      A2137 = A(TWO,ONE) * B(ONE,SEVEN)
     >      + A(TWO,THREE) * B(THREE,SEVEN)
      A2341 = A(TWO,THREE) * B(THREE,ONE)
     >      + A(TWO,FOUR) * B(FOUR,ONE)
      A2543 = A(TWO,FIVE) * B(FIVE,THREE)
     >      + A(TWO,FOUR) * B(FOUR,THREE)
      A2534 = A(TWO,FIVE) * B(FIVE,FOUR)
     >      + A(TWO,THREE) * B(THREE,FOUR)
      A2567 = A(TWO,FIVE) * B(FIVE,SEVEN)
     >      + A(TWO,SIX) * B(SIX,SEVEN)
      A6173 = A(SIX,ONE) * B(ONE,THREE)
     >      + A(SIX,SEVEN) * B(SEVEN,THREE)
      A6534 = A(SIX,FIVE) * B(FIVE,FOUR)
     >      + A(SIX,THREE) * B(THREE,FOUR)
      A6543 = A(SIX,FIVE) * B(FIVE,THREE)
     >      + A(SIX,FOUR) * B(FOUR,THREE)
      A617342 = A(SIX,ONE) * B(ONE,THREE) * A(THREE,TWO)
     >        + A(SIX,ONE) * B(ONE,FOUR) * A(FOUR,TWO)
     >        + A(SIX,SEVEN) * B(SEVEN,THREE) * A(THREE,TWO)
     >        + A(SIX,SEVEN) * B(SEVEN,FOUR) * A(FOUR,TWO)

      B12 = B(ONE,TWO)
      B13 = B(ONE,THREE)
      B17 = B(ONE,SEVEN)
      B34 = B(THREE,FOUR)

      S23 = S(TWO,THREE)
      S34 = S(THREE,FOUR)
      S67 = S(SIX,SEVEN)

      T123 = S(ONE,TWO) + S(ONE,THREE) + S(TWO,THREE)
      T167 = S(ONE,SIX) + S(ONE,SEVEN) + S(SIX,SEVEN)
      T234 = S(TWO,THREE) + S(TWO,FOUR) + S(THREE,FOUR)
      T345 = S(THREE,FOUR) + S(THREE,FIVE) + S(FOUR,FIVE)
      T567 = S(FIVE,SIX) + S(FIVE,SEVEN) + S(SIX,SEVEN)

      qggqbA70PMPPMNT =
     > -(B34**2 * A2341 * A2567 * A56)/(S23 * S34 * T234 * T567)
     > +(B13 * A2341 * A2567 * A56)/(B12 * A34 * A24 * S23 * T567)
     > +(B13**2 * A2134 * A2567 * A56)/(B12 * A24 * S23 * T123 * T567)
     > +(B13**2 * A2137 * A56 * A25)/(B12 * A24 * A45 * S23 * T123)
     > -(B17 * A25 * B34**2 * A617342)/(S23 * S34 * T234 * T167)
     > +(B13 * B17 * A25)/(B12 * A24 * S23 * T345)
     >  *(A6543 * (A25/A45 + A23/A34) + A6534 * A24/A34)
     > +(B17 * A6173 * A25)/(A24 * S23 * T345 * T167)
     >  *(A2543 * (A25/A45 + A23/A34) + A2534 * A24/A34)
      
      qggqbA70PMPPMNT = (0,1)/S67 * qggqbA70PMPPMNT

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qggqbA70PPMMMNT
     >                          (ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN
      REAL(KIND(1D0)) S(11,11)
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMMON /DOTPRODUCTS/ S
      REAL(KIND(1D0)) S23,S34,S67,T123,T167,T234,T567
      COMPLEX(KIND(1D0)) A12,A13,A34,A56
     >          ,A3127,A3542,A4567,A5342,A6172,A6542
     >          ,B12,B17,B24,B34,B45
     >          ,B234567

      A12 = A(ONE,TWO)
      A13 = A(ONE,THREE)
      A34 = A(THREE,FOUR)
      A56 = A(FIVE,SIX)

      A3127 = A(THREE,ONE) * B(ONE,SEVEN)
     >      + A(THREE,TWO) * B(TWO,SEVEN)
      A3542 = A(THREE,FIVE) * B(FIVE,TWO)
     >      + A(THREE,FOUR) * B(FOUR,TWO)
      A4567 = A(FOUR,FIVE) * B(FIVE,SEVEN)
     >      + A(FOUR,SIX) * B(SIX,SEVEN)
      A5342 = A(FIVE,THREE) * B(THREE,TWO)
     >      + A(FIVE,FOUR) * B(FOUR,TWO)
      A6172 = A(SIX,ONE) * B(ONE,TWO) 
     >      + A(SIX,SEVEN) * B(SEVEN,TWO)
      A6542 = A(SIX,FIVE) * B(FIVE,TWO)
     >      + A(SIX,FOUR) * B(FOUR,TWO)

      B12 = B(ONE,TWO)
      B17 = B(ONE,SEVEN)
      B24 = B(TWO,FOUR)
      B34 = B(THREE,FOUR)
      B45 = B(FOUR,FIVE)

      B234567 = B(TWO,THREE) * A(THREE,FIVE) * B(FIVE,SEVEN)
     >        + B(TWO,THREE) * A(THREE,SIX) * B(SIX,SEVEN)
     >        + B(TWO,FOUR) * A(FOUR,FIVE) * B(FIVE,SEVEN)
     >        + B(TWO,FOUR) * A(FOUR,SIX) * B(SIX,SEVEN)

      S23 = S(TWO,THREE)
      S34 = S(THREE,FOUR)
      S67 = S(SIX,SEVEN)

      T123 = S(ONE,TWO) + S(ONE,THREE) + S(TWO,THREE)
      T167 = S(ONE,SIX) + S(ONE,SEVEN) + S(SIX,SEVEN)
      T234 = S(TWO,THREE) + S(TWO,FOUR) + S(THREE,FOUR)
      T567 = S(FIVE,SIX) + S(FIVE,SEVEN) + S(SIX,SEVEN)

      qggqbA70PPMMMNT =
     > -(B12 * B234567 * A56)/(S23 * T567)
     >  *(A34**2/(S34 * T234) - A13/(B34 * B24 * A12))
     > +(A13**2 * B12**2 * A4567 * A56)/(A12 * B24 * S23 * T123 * T567)
     > -(A13 * B12 * A3127 * A6542)/(A12 * B24 * B45 * S23 * T123)
     > +(A3127 * A6172)/(A12 * B34 * B45 * S23)
     > +(B17 * A6172 * A3542)/(B34 * B45 * S23 * T167)
     > -(B17 * A6172 * A5342 * A34**2)/(S23 * S34 * T234 * T167)

      qggqbA70PPMMMNT = (0,1)/S67 * qggqbA70PPMMMNT

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qggqbA70PMMPMNT
     >                          (ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN
      REAL(KIND(1D0)) S(11,11)
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMMON /DOTPRODUCTS/ S
      REAL(KIND(1D0)) S23,S34,S67,T123,T167,T234,T345,T567
      COMPLEX(KIND(1D0)) A23,A35,A45,A56
     >          ,A2134,A2137,A2534,A3124,A3127,A3567,A5234,A6174,A6534
     >          ,B12,B14,B17,B23,B24,B34
     >          ,B423567

      A23 = A(TWO,THREE)
      A35 = A(THREE,FIVE)
      A45 = A(FOUR,FIVE)
      A56 = A(FIVE,SIX)

      A2134 = A(TWO,ONE) * B(ONE,FOUR)
     >      + A(TWO,THREE) * B(THREE,FOUR)
      A2137 = A(TWO,ONE) * B(ONE,SEVEN)
     >      + A(TWO,THREE) * B(THREE,SEVEN)
      A2534 = A(TWO,FIVE) * B(FIVE,FOUR)
     >      + A(TWO,THREE) * B(THREE,FOUR)
      A3124 = A(THREE,ONE) * B(ONE,FOUR)
     >      + A(THREE,TWO) * B(TWO,FOUR)
      A3127 = A(THREE,ONE) * B(ONE,SEVEN)
     >      + A(THREE,TWO) * B(TWO,SEVEN)
      A3567 = A(THREE,FIVE) * B(FIVE,SEVEN)
     >      + A(THREE,SIX) * B(SIX,SEVEN)      
      A5234 = A(FIVE,TWO) * B(TWO,FOUR)
     >      + A(FIVE,THREE) * B(THREE,FOUR)
      A6174 = A(SIX,ONE) * B(ONE,FOUR) 
     >      + A(SIX,SEVEN) * B(SEVEN,FOUR)
      A6534 = A(SIX,FIVE) * B(FIVE,FOUR)
     >      + A(SIX,THREE) * B(THREE,FOUR)

      B12 = B(ONE,TWO)
      B14 = B(ONE,FOUR)
      B17 = B(ONE,SEVEN)
      B23 = B(TWO,THREE)
      B24 = B(TWO,FOUR)
      B34 = B(THREE,FOUR)

      B423567 = B(FOUR,TWO) * A(TWO,FIVE) * B(FIVE,SEVEN)
     >        + B(FOUR,TWO) * A(TWO,SIX) * B(SIX,SEVEN)
     >        + B(FOUR,THREE) * A(THREE,FIVE) * B(FIVE,SEVEN)
     >        + B(FOUR,THREE) * A(THREE,SIX) * B(SIX,SEVEN)

      S23 = S(TWO,THREE)
      S34 = S(THREE,FOUR)
      S67 = S(SIX,SEVEN)

      T123 = S(ONE,TWO) + S(ONE,THREE) + S(TWO,THREE)
      T167 = S(ONE,SIX) + S(ONE,SEVEN) + S(SIX,SEVEN)
      T234 = S(TWO,THREE) + S(TWO,FOUR) + S(THREE,FOUR)
      T345 = S(THREE,FOUR) + S(THREE,FIVE) + S(FOUR,FIVE)
      T567 = S(FIVE,SIX) + S(FIVE,SEVEN) + S(SIX,SEVEN)

      qggqbA70PMMPMNT =
     > -(A23**2 * B14 * B423567 * A56)/(S23 * S34 * T234 * T567)
     > -(B14 * A3567 * A56)/(B24 * S34 * T123 * T567)
     >  *(-(B24 * A2134 + B34 * A3124)/B23 - (B14 * A3124)/B12)
     > -(B14 * A56 * A35)/(B24 * A45 * S34 * T123)
     >  *(-(B24 * A2137 + B34 * A3127)/B23 - (B14 * A3127)/B12)
     > +(B14 * B17 * A6534 * A35**2)/(B12 * B24 * A45 * S34 * T345)
     > +(B17 * A6174 * A2534 * A35**2)/(B24 * A45 * S34 * T345 * T167)
     > +(B17 * A6174 * A5234 * A35)/(A45 * B24 * B23 * S34 * T167)
     > -(B17 * A6174 * A5234 * A23**2)/(S23 * S34 * T234 * T167)

      qggqbA70PMMPMNT = (0,1)/S67 * qggqbA70PMMPMNT

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qggqbA70PMPMMNT
     >                          (ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qggqbA70PPMPMNT

      qggqbA70PMPMMNT = 
     > qggqbA70PPMPMNT(FIVE,FOUR,THREE,TWO,ONE,SEVEN,SIX,B,A)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qggqbA70PMMMMNT
     >                          (ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qggqbA70PPPPMNT

      qggqbA70PMMMMNT = 
     > qggqbA70PPPPMNT(FIVE,FOUR,THREE,TWO,ONE,SEVEN,SIX,B,A)

      END
