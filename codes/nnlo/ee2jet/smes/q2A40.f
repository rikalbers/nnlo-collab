************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbA40PM(ONE,TWO,FOUR,FIVE,
     >                                      A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,FOUR,FIVE
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) A12, A25, A45

      A12 = A(ONE,TWO)
      A25 = A(TWO,FIVE)
      A45 = A(FOUR,FIVE)
c Factor 2 difference!
      qqbA40PM = -(0,1)*A25**2/(A12*A45)

      END

************************************************************************
