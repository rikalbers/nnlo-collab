************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qgqbA51PPM(ONE,TWO,THREE,FOUR,FIVE,
     >                                       A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qgqbA5PPM,qqbgA5PMP
      REAL(KIND(1D0)) XNf,XNu,XNd,XNc,XNa                  
      COMMON /QCD/ XNf,XNu,XNd,XNc,XNa

      qgqbA51PPM = qgqbA5PPM(ONE,TWO,THREE,FOUR,FIVE,A,B)
     >+ 1D0/XNc**2*qqbgA5PMP(ONE,THREE,TWO,FOUR,FIVE,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qgqbA51PPMmudep(ONE,TWO,THREE,FOUR,
     >                                            FIVE,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qgqbA5PPMmudep,qqbgA5PMPmudep
      REAL(KIND(1D0)) XNf,XNu,XNd,XNc,XNa                  
      COMMON /QCD/ XNf,XNu,XNd,XNc,XNa

      qgqbA51PPMmudep = qgqbA5PPMmudep(ONE,TWO,THREE,FOUR,FIVE,A,B)
     >+ 1D0/XNc**2*qqbgA5PMPmudep(ONE,THREE,TWO,FOUR,FIVE,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qgqbA51PPMmuindep(ONE,TWO,THREE,FOUR,
     >                                              FIVE,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qgqbA5PPMmuindep,qqbgA5PMPmuindep
      REAL(KIND(1D0)) XNf,XNu,XNd,XNc,XNa                  
      COMMON /QCD/ XNf,XNu,XNd,XNc,XNa

      qgqbA51PPMmuindep = qgqbA5PPMmuindep(ONE,TWO,THREE,FOUR,FIVE,A,B)
     >+ 1D0/XNc**2*qqbgA5PMPmuindep(ONE,THREE,TWO,FOUR,FIVE,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qgqbA51PPMem2(ONE,TWO,THREE,FOUR,FIVE,
     >                                          A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qgqbA5PPMem2,qqbgA5PMPem2
      REAL(KIND(1D0)) XNf,XNu,XNd,XNc,XNa                  
      COMMON /QCD/ XNf,XNu,XNd,XNc,XNa

      qgqbA51PPMem2 = qgqbA5PPMem2(ONE,TWO,THREE,FOUR,FIVE,A,B)
     >+ 1D0/XNc**2*qqbgA5PMPem2(ONE,THREE,TWO,FOUR,FIVE,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qgqbA51PPMem1(ONE,TWO,THREE,FOUR,FIVE,
     >                                          A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qgqbA5PPMem1,qqbgA5PMPem1
      REAL(KIND(1D0)) XNf,XNu,XNd,XNc,XNa                  
      COMMON /QCD/ XNf,XNu,XNd,XNc,XNa

      qgqbA51PPMem1 = qgqbA5PPMem1(ONE,TWO,THREE,FOUR,FIVE,A,B)
     >+ 1D0/XNc**2*qqbgA5PMPem1(ONE,THREE,TWO,FOUR,FIVE,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qgqbA51PMM(ONE,TWO,THREE,FOUR,FIVE,
     >                                       A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qgqbA5PMM,qqbgA5PMM
      REAL(KIND(1D0)) XNf,XNu,XNd,XNc,XNa                  
      COMMON /QCD/ XNf,XNu,XNd,XNc,XNa

      qgqbA51PMM = qgqbA5PMM(ONE,TWO,THREE,FOUR,FIVE,A,B)
     >+ 1D0/XNc**2*qqbgA5PMM(ONE,THREE,TWO,FOUR,FIVE,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qgqbA51PMMmudep(ONE,TWO,THREE,FOUR,
     >                                            FIVE,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qgqbA5PMMmudep,qqbgA5PMMmudep
      REAL(KIND(1D0)) XNf,XNu,XNd,XNc,XNa                  
      COMMON /QCD/ XNf,XNu,XNd,XNc,XNa

      qgqbA51PMMmudep = qgqbA5PMMmudep(ONE,TWO,THREE,FOUR,FIVE,A,B)
     >+ 1D0/XNc**2*qqbgA5PMMmudep(ONE,THREE,TWO,FOUR,FIVE,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qgqbA51PMMmuindep(ONE,TWO,THREE,FOUR,
     >                                              FIVE,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qgqbA5PMMmuindep,qqbgA5PMMmuindep
      REAL(KIND(1D0)) XNf,XNu,XNd,XNc,XNa                  
      COMMON /QCD/ XNf,XNu,XNd,XNc,XNa

      qgqbA51PMMmuindep = qgqbA5PMMmuindep(ONE,TWO,THREE,FOUR,FIVE,A,B)
     >+ 1D0/XNc**2*qqbgA5PMMmuindep(ONE,THREE,TWO,FOUR,FIVE,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qgqbA51PMMem2(ONE,TWO,THREE,FOUR,FIVE,
     >                                          A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qgqbA5PMMem2,qqbgA5PMMem2
      REAL(KIND(1D0)) XNf,XNu,XNd,XNc,XNa                  
      COMMON /QCD/ XNf,XNu,XNd,XNc,XNa

      qgqbA51PMMem2 = qgqbA5PMMem2(ONE,TWO,THREE,FOUR,FIVE,A,B)
     >+ 1D0/XNc**2*qqbgA5PMMem2(ONE,THREE,TWO,FOUR,FIVE,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qgqbA51PMMem1(ONE,TWO,THREE,FOUR,FIVE,
     >                                          A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qgqbA5PMMem1,qqbgA5PMMem1
      REAL(KIND(1D0)) XNf,XNu,XNd,XNc,XNa                  
      COMMON /QCD/ XNf,XNu,XNd,XNc,XNa

      qgqbA51PMMem1 = qgqbA5PMMem1(ONE,TWO,THREE,FOUR,FIVE,A,B)
     >+ 1D0/XNc**2*qqbgA5PMMem1(ONE,THREE,TWO,FOUR,FIVE,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qgqbA5PPM(ONE,TWO,THREE,FOUR,FIVE,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qgqbA50PPM,qgqbV,qgqbFccPPM,qgqbFscPPM

      qgqbA5PPM = qgqbA50PPM(ONE,TWO,THREE,FOUR,FIVE,A,B)
     >            *qgqbV(ONE,TWO,THREE,FOUR,FIVE)
     > + CMPLX(0D0,1D0)*qgqbFccPPM(ONE,TWO,THREE,FOUR,FIVE,A,B)
     > + CMPLX(0D0,1D0)*qgqbFscPPM(ONE,TWO,THREE,FOUR,FIVE,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qgqbA5PPMmudep(ONE,TWO,THREE,FOUR,
     >                                           FIVE,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qgqbA50PPM,qgqbV

      qgqbA5PPMmudep = qgqbA50PPM(ONE,TWO,THREE,FOUR,FIVE,A,B)
     >               * qgqbV(ONE,TWO,THREE,FOUR,FIVE)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qgqbA5PPMmuindep(ONE,TWO,THREE,FOUR,
     >                                             FIVE,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qgqbFccPPM,qgqbFscPPM

      qgqbA5PPMmuindep = 
     > + CMPLX(0D0,1D0)*qgqbFccPPM(ONE,TWO,THREE,FOUR,FIVE,A,B)
     > + CMPLX(0D0,1D0)*qgqbFscPPM(ONE,TWO,THREE,FOUR,FIVE,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qgqbA5PPMem2(ONE,TWO,THREE,FOUR,FIVE,
     >                                         A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qgqbA50PPM,qgqbVem2

      qgqbA5PPMem2 = qgqbA50PPM(ONE,TWO,THREE,FOUR,FIVE,A,B)
     >             *qgqbVem2(ONE,TWO,THREE,FOUR,FIVE)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qgqbA5PPMem1(ONE,TWO,THREE,FOUR,FIVE,
     >                                         A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qgqbA50PPM,qgqbVem1

      qgqbA5PPMem1 = qgqbA50PPM(ONE,TWO,THREE,FOUR,FIVE,A,B)
     >             *qgqbVem1(ONE,TWO,THREE,FOUR,FIVE)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qgqbA5PMM(ONE,TWO,THREE,FOUR,FIVE,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qgqbA5PPM

      qgqbA5PMM = -qgqbA5PPM(THREE,TWO,ONE,FIVE,FOUR,B,A)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qgqbA5PMMmudep(ONE,TWO,THREE,FOUR,
     >                                           FIVE,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qgqbA5PPMmudep

      qgqbA5PMMmudep = -qgqbA5PPMmudep(THREE,TWO,ONE,FIVE,FOUR,B,A)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qgqbA5PMMmuindep(ONE,TWO,THREE,FOUR,
     >                                             FIVE,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qgqbA5PPMmuindep

      qgqbA5PMMmuindep = -qgqbA5PPMmuindep(THREE,TWO,ONE,FIVE,FOUR,B,A)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qgqbA5PMMem2(ONE,TWO,THREE,FOUR,FIVE,
     >                                         A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qgqbA5PPMem2

      qgqbA5PMMem2 = -qgqbA5PPMem2(THREE,TWO,ONE,FIVE,FOUR,B,A)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qgqbA5PMMem1(ONE,TWO,THREE,FOUR,FIVE,
     >                                         A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qgqbA5PPMem1

      qgqbA5PMMem1 = -qgqbA5PPMem1(THREE,TWO,ONE,FIVE,FOUR,B,A)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbgA5PMP(ONE,TWO,THREE,FOUR,FIVE,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qqbgA50PMP,qqbgV,qqbgFccPMP,qqbgFscPMP

      qqbgA5PMP = qqbgA50PMP(ONE,TWO,THREE,FOUR,FIVE,A,B)
     >            *qqbgV(ONE,TWO,THREE,FOUR,FIVE)
     > + CMPLX(0D0,1D0)*qqbgFccPMP(ONE,TWO,THREE,FOUR,FIVE,A,B)
     > + CMPLX(0D0,1D0)*qqbgFscPMP(ONE,TWO,THREE,FOUR,FIVE,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbgA5PMPmudep(ONE,TWO,THREE,FOUR,
     >                                           FIVE,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qqbgA50PMP,qqbgV

      qqbgA5PMPmudep = qqbgA50PMP(ONE,TWO,THREE,FOUR,FIVE,A,B)
     >               * qqbgV(ONE,TWO,THREE,FOUR,FIVE)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbgA5PMPmuindep(ONE,TWO,THREE,FOUR,
     >                                             FIVE,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qqbgFccPMP,qqbgFscPMP

      qqbgA5PMPmuindep = 
     > + CMPLX(0D0,1D0)*qqbgFccPMP(ONE,TWO,THREE,FOUR,FIVE,A,B)
     > + CMPLX(0D0,1D0)*qqbgFscPMP(ONE,TWO,THREE,FOUR,FIVE,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbgA5PMPem2(ONE,TWO,THREE,FOUR,FIVE,
     >                                         A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qqbgA50PMP,qqbgVem2

      qqbgA5PMPem2 = qqbgA50PMP(ONE,TWO,THREE,FOUR,FIVE,A,B)
     >             * qqbgVem2(ONE,TWO,THREE,FOUR,FIVE)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbgA5PMPem1(ONE,TWO,THREE,FOUR,FIVE,
     >                                         A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qqbgA50PMP,qqbgVem1

      qqbgA5PMPem1 = qqbgA50PMP(ONE,TWO,THREE,FOUR,FIVE,A,B)
     >             * qqbgVem1(ONE,TWO,THREE,FOUR,FIVE)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbgA5PMM(ONE,TWO,THREE,FOUR,FIVE,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qqbgA5PMP

      qqbgA5PMM = -qqbgA5PMP(TWO,ONE,THREE,FIVE,FOUR,B,A)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbgA5PMMmudep(ONE,TWO,THREE,FOUR,
     >                                           FIVE,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qqbgA5PMPmudep

      qqbgA5PMMmudep = -qqbgA5PMPmudep(TWO,ONE,THREE,FIVE,FOUR,B,A)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbgA5PMMmuindep(ONE,TWO,THREE,FOUR,
     >                                             FIVE,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qqbgA5PMPmuindep

      qqbgA5PMMmuindep = -qqbgA5PMPmuindep(TWO,ONE,THREE,FIVE,FOUR,B,A)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbgA5PMMem2(ONE,TWO,THREE,FOUR,FIVE,
     >                                         A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qqbgA5PMPem2

      qqbgA5PMMem2 = -qqbgA5PMPem2(TWO,ONE,THREE,FIVE,FOUR,B,A)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbgA5PMMem1(ONE,TWO,THREE,FOUR,FIVE,
     >                                         A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qqbgA5PMPem1

      qqbgA5PMMem1 = -qqbgA5PMPem1(TWO,ONE,THREE,FIVE,FOUR,B,A)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbgA5axPMP(ONE,TWO,THREE,FOUR,FIVE,
     >                                        A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qqbgFaxPMP

      qqbgA5axPMP =
     >   CMPLX(0D0,1D0)*qqbgFaxPMP(ONE,TWO,THREE,FOUR,FIVE,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbgA5axPMM(ONE,TWO,THREE,FOUR,FIVE,
     >                                        A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qqbgA5axPMP

      qqbgA5axPMM = qqbgA5axPMP(TWO,ONE,THREE,FIVE,FOUR,B,A)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qgqbV(ONE,TWO,THREE,FOUR,FIVE)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE
      REAL(KIND(1D0)) S(11,11)
      REAL(KIND(1D0)) MUSQ,COMU
      COMMON /DOTPRODUCTS/ S
     >       /scales/ COMU
      COMPLEX(KIND(1D0)) LOGS12,LOGS23
      REAL(KIND(1D0)) PI
      PARAMETER ( PI  = 3.14159265358979 D0)
      REAL(KIND(1D0)) THETA
      EXTERNAL THETA

      MUSQ = COMU**2
      LOGS12 = CMPLX(LOG(ABS(S(ONE,TWO))/MUSQ),
     >               -PI*THETA(S(ONE,TWO).GT.0D0))
      LOGS23 = CMPLX(LOG(ABS(S(TWO,THREE))/MUSQ),
     >               -PI*THETA(S(TWO,THREE).GT.0D0))

      qgqbV = -0.5D0*LOGS12**2 - 0.5D0*LOGS23**2
     >        +1.5D0*LOGS23 - 3D0

c eps^-2:
c      qgqbV = -2d0

c eps^-1:
c      qgqbV = 1d0*(LOGS12 + LOGS23) - 1d0*1.5d0

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qgqbVem2(ONE,TWO,THREE,FOUR,FIVE)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE

      qgqbVem2 = -2d0

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qgqbVem1(ONE,TWO,THREE,FOUR,FIVE)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE
      REAL(KIND(1D0)) S(11,11)
      REAL(KIND(1D0)) MUSQ,COMU
      COMMON /DOTPRODUCTS/ S
     >       /scales/ COMU
      COMPLEX(KIND(1D0)) LOGS12,LOGS23
      REAL(KIND(1D0)) PI
      PARAMETER ( PI  = 3.14159265358979 D0)
      REAL(KIND(1D0)) THETA
      EXTERNAL THETA

      MUSQ = COMU**2
      LOGS12 = CMPLX(LOG(ABS(S(ONE,TWO))/MUSQ),
     >               -PI*THETA(S(ONE,TWO).GT.0D0))
      LOGS23 = CMPLX(LOG(ABS(S(TWO,THREE))/MUSQ),
     >               -PI*THETA(S(TWO,THREE).GT.0D0))

c eps^-1:
      qgqbVem1 = 1d0*(LOGS12 + LOGS23) - 1d0*1.5d0

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbgV(ONE,TWO,THREE,FOUR,FIVE)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE
      REAL(KIND(1D0)) S(11,11)
      REAL(KIND(1D0)) MUSQ,COMU
      COMMON /DOTPRODUCTS/ S
     >       /scales/ COMU
      COMPLEX(KIND(1D0)) LOGS12,LOGS45
      REAL(KIND(1D0)) PI
      PARAMETER ( PI  = 3.14159265358979 D0)
      REAL(KIND(1D0)) THETA
      EXTERNAL THETA

      MUSQ = COMU**2

      LOGS12 = CMPLX(LOG(ABS(S(ONE,TWO))/MUSQ),
     >               -PI*THETA(S(ONE,TWO).GT.0D0))
      LOGS45 = CMPLX(LOG(ABS(S(FOUR,FIVE))/MUSQ),
     >               -PI*THETA(S(FOUR,FIVE).GT.0D0))

      qqbgV = -0.5D0*LOGS12**2 + 1.5D0*LOGS45 - 3.5D0

c eps^-2:
c      qqbgV = -1d0

C eps^-1:
c      qqbgV = 1d0*LOGS12 - 1d0*1.5d0

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbgVem2(ONE,TWO,THREE,FOUR,FIVE)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE

c eps^-2:
      qqbgVem2 = -1d0

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbgVem1(ONE,TWO,THREE,FOUR,FIVE)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE
      REAL(KIND(1D0)) S(11,11)
      REAL(KIND(1D0)) MUSQ,COMU
      COMMON /DOTPRODUCTS/ S
     >       /scales/ COMU
      COMPLEX(KIND(1D0)) LOGS12,LOGS45
      REAL(KIND(1D0)) PI
      PARAMETER ( PI  = 3.14159265358979 D0)
      REAL(KIND(1D0)) THETA
      EXTERNAL THETA

      MUSQ = COMU**2

      LOGS12 = CMPLX(LOG(ABS(S(ONE,TWO))/MUSQ),
     >               -PI*THETA(S(ONE,TWO).GT.0D0))
      LOGS45 = CMPLX(LOG(ABS(S(FOUR,FIVE))/MUSQ),
     >               -PI*THETA(S(FOUR,FIVE).GT.0D0))

C eps^-1:
      qqbgVem1 = 1d0*LOGS12 - 1d0*1.5d0

      END

************************************************************************
************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qgqbFccPPM(ONE,TWO,THREE,FOUR,FIVE,
     >                                       A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      REAL(KIND(1D0)) S(11,11)
      COMMON /DOTPRODUCTS/ S
      REAL(KIND(1D0)) S12, S23, S45
      COMPLEX(KIND(1D0)) A12, A13, A23, A34, A45 ,B15 ,L0, Lsm1
      
      A12 = A(ONE,TWO)
      A13 = A(ONE,THREE)
      A23 = A(TWO,THREE)
      A34 = A(THREE,FOUR)
      A45 = A(FOUR,FIVE)
      B15 = B(ONE,FIVE)
      S12 = S(ONE,TWO)
      S23 = S(TWO,THREE)
      S45 = S(FOUR,FIVE)

      qgqbFccPPM = A34**2/(A12*A23*A45)*
     - (Lsm1(-S12,-S45,-S23,-S45)-2D0*A13*B15*A45/A34*L0(-S23,-S45)/S45)

      END                                                   

************************************************************************
      
      COMPLEX(KIND(1D0)) FUNCTION qgqbFscPPM(ONE,TWO,THREE,FOUR,FIVE,
     >                                       A,B)
      IMPLICIT NONE 
      
      INTEGER ONE,TWO,THREE,FOUR,FIVE
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      REAL(KIND(1D0)) S(11,11)
      COMMON /DOTPRODUCTS/ S
      REAL(KIND(1D0)) S23, S45
      COMPLEX(KIND(1D0)) A12, A13, A23, A34, A45 ,B15 ,L0, L1
      
      A12 = A(ONE,TWO)
      A13 = A(ONE,THREE)
      A23 = A(TWO,THREE)
      A34 = A(THREE,FOUR)
      A45 = A(FOUR,FIVE)
      B15 = B(ONE,FIVE)
      S23 = S(TWO,THREE)
      S45 = S(FOUR,FIVE)

      qgqbFscPPM = A13*B15*A45/(A12*A23*A45)*
     - (A34*L0(-S23,-S45)/S45 + 0.5D0*A13*B15*A45*L1(-S23,-S45)/S45**2)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbgFccPMP(ONE,TWO,THREE,FOUR,FIVE,
     >                                       A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      REAL(KIND(1D0)) S(11,11)
      COMMON /DOTPRODUCTS/ S
      REAL(KIND(1D0)) S12, S13, S23, S45
      COMPLEX(KIND(1D0)) A12,A13,A14,A23,A24,A34,A45,B13,L0,Lsm1
      
      A12 = A(ONE,TWO)
      A13 = A(ONE,THREE)
      A14 = A(ONE,FOUR)
      A23 = A(TWO,THREE)
      A24 = A(TWO,FOUR)
      A34 = A(THREE,FOUR)
      A45 = A(FOUR,FIVE)
      B13 = B(ONE,THREE)
      S12 = S(ONE,TWO)
      S13 = S(ONE,THREE)
      S23 = S(TWO,THREE)
      S45 = S(FOUR,FIVE)

      qqbgFccPMP =
     -  A24**2/(A23*A13*A45)*Lsm1(-S12,-S45,-S13,-S45)
     - +A24*(A12*A34-A14*A23)/(A23*A13**2*A45)*Lsm1(-S12,-S45,-S23,-S45)
     - +2D0*B13*A14*A24/(A13*A45)*L0(-S23,-S45)/S45

      END                                                   

************************************************************************
      
      COMPLEX(KIND(1D0)) FUNCTION qqbgFscPMP(ONE,TWO,THREE,FOUR,FIVE, 
     >                                       A,B)
      IMPLICIT NONE 
      
      INTEGER ONE,TWO,THREE,FOUR,FIVE
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      REAL(KIND(1D0)) S(11,11)
      COMMON /DOTPRODUCTS/ S
      REAL(KIND(1D0)) S12, S23, S45
      COMPLEX(KIND(1D0)) A12, A13, A14, A23, A34, A45 
     >          ,B12, B13, B15, B23, B25, B35, B45
     >          ,L0, L1, Lsm1
      
      A12 = A(ONE,TWO)
      A13 = A(ONE,THREE)
      A14 = A(ONE,FOUR)
      A23 = A(TWO,THREE)
      A34 = A(THREE,FOUR)
      A45 = A(FOUR,FIVE)
      B12 = B(ONE,TWO)
      B13 = B(ONE,THREE)
      B15 = B(ONE,FIVE)
      B23 = B(TWO,THREE)
      B25 = B(TWO,FIVE)
      B35 = B(THREE,FIVE)
      B45 = B(FOUR,FIVE)
      S12 = S(ONE,TWO)
      S23 = S(TWO,THREE)
      S45 = S(FOUR,FIVE)

      qqbgFscPMP =
     -  A14**2*A23/(A13**3*A45)*Lsm1(-S12,-S45,-S23,-S45)
     - -0.5D0*(A14*B13)**2*A23/(A13*A45)*L1(-S45,-S23)/S23**2
     - -A14**2*A23*B13/(A13**2*A45)*L0(-S45,-S23)/S23
     - -A12*B13*A34*B35/A13*L1(-S45,-S12)/S12**2
     - +A12*B13*A34*A14/(A13**2*A45)*L0(-S45,-S12)/S12
     - -0.5D0*B35*(B13*B25+B23*B15)/(B12*B23*A13*B45)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbgFaxPMP(ONE,TWO,THREE,FOUR,FIVE,
     >                                       A,B)
      IMPLICIT NONE 

      INTEGER ONE,TWO,THREE,FOUR,FIVE
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      REAL(KIND(1D0)) S(11,11),XMT,mtsq
      COMMON /DOTPRODUCTS/ S
     >       /TOPMASS/ XMT,mtsq
      REAL(KIND(1D0)) S12, S45
      COMPLEX(KIND(1D0)) A24,B13,B35,L1

      A24 = A(TWO,FOUR)
      B13 = B(ONE,THREE)
      B35 = B(THREE,FIVE)
      S12 = S(ONE,TWO)
      S45 = S(FOUR,FIVE)

c No top contribution ... AK 
      qqbgFaxPMP = 
c     >  - A24*B13*B35*(L1(-S12,-S45)/S45**2-1D0/(12D0*S45*mtsq))
     >  - A24*B13*B35*(L1(-S12,-S45)/S45**2)

      END

************************************************************************
