************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qgqbA50PPM(ONE,TWO,THREE,FOUR,FIVE,
     >                                       A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) A12, A23, A34, A45

      A12 = A(ONE,TWO)
      A23 = A(TWO,THREE)
      A34 = A(THREE,FOUR)
      A45 = A(FOUR,FIVE)

      qgqbA50PPM = -(0,1)*A34**2/(A12*A23*A45)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbgA50PMP(ONE,TWO,THREE,FOUR,FIVE,
     >                                       A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) A24, A23, A31, A45

      A23 = A(TWO,THREE)
      A24 = A(TWO,FOUR)
      A31 = A(THREE,ONE)
      A45 = A(FOUR,FIVE)

      qqbgA50PMP = (0,1)*A24**2/(A23*A31*A45)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qgqbA50PMM(ONE,TWO,THREE,FOUR,FIVE,
     >                                       A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qgqbA50PPM

! FIXME minus -> plus
      qgqbA50PMM = -qgqbA50PPM(THREE,TWO,ONE,FIVE,FOUR,B,A)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbgA50PMM(ONE,TWO,THREE,FOUR,FIVE,
     >                                       A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qqbgA50PMP

! FIXME minus -> plus
      qqbgA50PMM = -qqbgA50PMP(TWO,ONE,THREE,FIVE,FOUR,B,A)

      END

************************************************************************
