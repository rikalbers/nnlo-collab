************************************************************************
C This file edited by SG 2005.03.09
C Introduced qggqbV, qgqbgV, qqbggV
************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qggqbA61PMPM(ONE,TWO,THREE,FOUR,FIVE,
     >                                         SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qggqbA6PMPM,qqbggA6PMPM,qggqbA6SPMPM,
     >                   qggqbA6TPMPM
      REAL(KIND(1D0)) XNf,XNu,XNd,XNc,XNa                      
      COMMON /QCD/ XNf,XNu,XNd,XNc,XNa

      qggqbA61PMPM = qggqbA6PMPM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >  - 1D0/XNc**2*qqbggA6PMPM(ONE,FOUR,THREE,TWO,FIVE,SIX,A,B)
     >     - XNf/XNc*qggqbA6SPMPM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >     + 1D0/XNc*qggqbA6TPMPM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qggqbA61PPMM(ONE,TWO,THREE,FOUR,FIVE,
     >                                         SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qggqbA6PPMM,qqbggA6PMMP,qggqbA6SPPMM,
     >                   qggqbA6TPPMM
      REAL(KIND(1D0)) XNf,XNu,XNd,XNc,XNa                      
      COMMON /QCD/ XNf,XNu,XNd,XNc,XNa

      qggqbA61PPMM = qggqbA6PPMM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >  - 1D0/XNc**2*qqbggA6PMMP(ONE,FOUR,THREE,TWO,FIVE,SIX,A,B)
     >     - XNf/XNc*qggqbA6SPPMM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >     + 1D0/XNc*qggqbA6TPPMM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qggqbA61PPPM(ONE,TWO,THREE,FOUR,FIVE,
     >                                         SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qggqbA6PPPM,qqbggA6PMPP,qggqbA6SPPPM,
     >                   qggqbA6TPPPM
      REAL(KIND(1D0)) XNf,XNu,XNd,XNc,XNa                      
      COMMON /QCD/ XNf,XNu,XNd,XNc,XNa

      qggqbA61PPPM = qggqbA6PPPM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >  - 1D0/XNc**2*qqbggA6PMPP(ONE,FOUR,THREE,TWO,FIVE,SIX,A,B)
     >     - XNf/XNc*qggqbA6SPPPM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >     + 1D0/XNc*qggqbA6TPPPM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qggqbA61PMMM(ONE,TWO,THREE,FOUR,FIVE,
     >                                         SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qggqbA6PMMM,qqbggA6PMMM,qggqbA6SPMMM,
     >                   qggqbA6TPMMM
      REAL(KIND(1D0)) XNf,XNu,XNd,XNc,XNa                      
      COMMON /QCD/ XNf,XNu,XNd,XNc,XNa

      qggqbA61PMMM = qggqbA6PMMM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >  - 1D0/XNc**2*qqbggA6PMMM(ONE,FOUR,THREE,TWO,FIVE,SIX,A,B)
     >     - XNf/XNc*qggqbA6SPMMM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >     + 1D0/XNc*qggqbA6TPMMM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)

      END

************************************************************************
************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbggA63PMPM(ONE,FOUR,TWO,THREE,FIVE,
     >                                         SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,FOUR,TWO,THREE,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qggqbA6PPMM,qggqbA6PMPM,qgqbgA6PPMM,qgqbgA6PMMP
     >          ,qqbggA6PMPM,qqbggA6PMMP

      qqbggA63PMPM = qggqbA6PPMM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >             + qggqbA6PMPM(ONE,THREE,TWO,FOUR,FIVE,SIX,A,B)
     >             + qgqbgA6PPMM(ONE,TWO,FOUR,THREE,FIVE,SIX,A,B)
     >             + qgqbgA6PMMP(ONE,THREE,FOUR,TWO,FIVE,SIX,A,B)
     >             + qqbggA6PMPM(ONE,FOUR,TWO,THREE,FIVE,SIX,A,B)
     >             + qqbggA6PMMP(ONE,FOUR,THREE,TWO,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbggA63PMMP(ONE,FOUR,TWO,THREE,FIVE,
     >                                         SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,FOUR,TWO,THREE,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qggqbA6PPMM,qggqbA6PMPM,qgqbgA6PPMM,qgqbgA6PMMP
     >          ,qqbggA6PMPM,qqbggA6PMMP

      qqbggA63PMMP = qggqbA6PMPM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >             + qggqbA6PPMM(ONE,THREE,TWO,FOUR,FIVE,SIX,A,B)
     >             + qgqbgA6PMMP(ONE,TWO,FOUR,THREE,FIVE,SIX,A,B)
     >             + qgqbgA6PPMM(ONE,THREE,FOUR,TWO,FIVE,SIX,A,B)
     >             + qqbggA6PMMP(ONE,FOUR,TWO,THREE,FIVE,SIX,A,B)
     >             + qqbggA6PMPM(ONE,FOUR,THREE,TWO,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbggA63PMPP(ONE,FOUR,TWO,THREE,FIVE,
     >                                         SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,FOUR,TWO,THREE,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qggqbA6PPPM,qgqbgA6PPMP,qqbggA6PMPP

      qqbggA63PMPP = qggqbA6PPPM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >             + qggqbA6PPPM(ONE,THREE,TWO,FOUR,FIVE,SIX,A,B)
     >             + qgqbgA6PPMP(ONE,TWO,FOUR,THREE,FIVE,SIX,A,B)
     >             + qgqbgA6PPMP(ONE,THREE,FOUR,TWO,FIVE,SIX,A,B)
     >             + qqbggA6PMPP(ONE,FOUR,TWO,THREE,FIVE,SIX,A,B)
     >             + qqbggA6PMPP(ONE,FOUR,THREE,TWO,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbggA63PMMM(ONE,FOUR,TWO,THREE,FIVE,
     >                                         SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,FOUR,TWO,THREE,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qggqbA6PMMM,qgqbgA6PMMM,qqbggA6PMMM

      qqbggA63PMMM = qggqbA6PMMM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >             + qggqbA6PMMM(ONE,THREE,TWO,FOUR,FIVE,SIX,A,B)
     >             + qgqbgA6PMMM(ONE,TWO,FOUR,THREE,FIVE,SIX,A,B)
     >             + qgqbgA6PMMM(ONE,THREE,FOUR,TWO,FIVE,SIX,A,B)
     >             + qqbggA6PMMM(ONE,FOUR,TWO,THREE,FIVE,SIX,A,B)
     >             + qqbggA6PMMM(ONE,FOUR,THREE,TWO,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbggA64vPMPM(ONE,FOUR,TWO,THREE,FIVE,
     >                                          SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,FOUR,TWO,THREE,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qqbggA6vsPMPM,qqbggA6vfPMPM

      qqbggA64vPMPM =- qqbggA6vsPMPM(ONE,FOUR,TWO,THREE,FIVE,SIX,A,B)
     >               - qqbggA6vfPMPM(ONE,FOUR,TWO,THREE,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbggA64vPMMP(ONE,FOUR,TWO,THREE,FIVE,
     >                                          SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,FOUR,TWO,THREE,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qqbggA6vsPMMP,qqbggA6vfPMMP

      qqbggA64vPMMP =- qqbggA6vsPMMP(ONE,FOUR,TWO,THREE,FIVE,SIX,A,B)
     >               - qqbggA6vfPMMP(ONE,FOUR,TWO,THREE,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbggA64vPMPP(ONE,FOUR,TWO,THREE,FIVE,
     >                                          SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,FOUR,TWO,THREE,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qqbggA6vsPMPP,qqbggA6vfPMPP

      qqbggA64vPMPP =- qqbggA6vsPMPP(ONE,FOUR,TWO,THREE,FIVE,SIX,A,B)
     >               - qqbggA6vfPMPP(ONE,FOUR,TWO,THREE,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbggA64vPMMM(ONE,FOUR,TWO,THREE,FIVE,
     >                                          SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,FOUR,TWO,THREE,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qqbggA6vsPMMM,qqbggA6vfPMMM

      qqbggA64vPMMM =- qqbggA6vsPMMM(ONE,FOUR,TWO,THREE,FIVE,SIX,A,B)
     >               - qqbggA6vfPMMM(ONE,FOUR,TWO,THREE,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION 
     >           qqbggA64axPMPM(ONE,FOUR,TWO,THREE,FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,FOUR,TWO,THREE,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qqbggA6axPMPM

      qqbggA64axPMPM = qqbggA6axPMPM(ONE,FOUR,TWO,THREE,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION 
     >           qqbggA64axPMMP(ONE,FOUR,TWO,THREE,FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,FOUR,TWO,THREE,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qqbggA6axPMMP

      qqbggA64axPMMP = qqbggA6axPMMP(ONE,FOUR,TWO,THREE,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION 
     >           qqbggA64axPMPP(ONE,FOUR,TWO,THREE,FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,FOUR,TWO,THREE,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qqbggA6axPMPP

      qqbggA64axPMPP = qqbggA6axPMPP(ONE,FOUR,TWO,THREE,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION 
     >           qqbggA64axPMMM(ONE,FOUR,TWO,THREE,FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,FOUR,TWO,THREE,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qqbggA6axPMMM

      qqbggA64axPMMM = qqbggA6axPMMM(ONE,FOUR,TWO,THREE,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION 
     >           qqbggA65axPMPM(ONE,FOUR,TWO,THREE,FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,FOUR,TWO,THREE,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qqbggA6axslPMPM

      qqbggA65axPMPM = qqbggA6axslPMPM(ONE,FOUR,TWO,THREE,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION 
     >           qqbggA65axPMMP(ONE,FOUR,TWO,THREE,FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,FOUR,TWO,THREE,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qqbggA6axslPMMP

      qqbggA65axPMMP = qqbggA6axslPMMP(ONE,FOUR,TWO,THREE,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION 
     >           qqbggA65axPMPP(ONE,FOUR,TWO,THREE,FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,FOUR,TWO,THREE,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qqbggA6axslPMPP

      qqbggA65axPMPP = qqbggA6axslPMPP(ONE,FOUR,TWO,THREE,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION 
     >           qqbggA65axPMMM(ONE,FOUR,TWO,THREE,FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,FOUR,TWO,THREE,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qqbggA6axslPMMM

      qqbggA65axPMMM = qqbggA6axslPMMM(ONE,FOUR,TWO,THREE,FIVE,SIX,A,B)

      END

************************************************************************
************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qggqbA6PMPM(ONE,TWO,THREE,FOUR,FIVE,
     >                                        SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qggqbA60PMPM,qggqbV,qggqbFccPMPM,qggqbFscPMPM

      qggqbA6PMPM = qggqbA60PMPM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >              *qggqbV(ONE,TWO,THREE,FOUR,FIVE,SIX)
     > + CMPLX(0D0,1D0)*qggqbFccPMPM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     > + CMPLX(0D0,1D0)*qggqbFscPMPM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qggqbA6PPMM(ONE,TWO,THREE,FOUR,FIVE,
     >                                        SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qggqbA60PPMM,qggqbV,qggqbFccPPMM,qggqbFscPPMM

      qggqbA6PPMM = qggqbA60PPMM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >              *qggqbV(ONE,TWO,THREE,FOUR,FIVE,SIX)
     > + CMPLX(0D0,1D0)*qggqbFccPPMM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     > + CMPLX(0D0,1D0)*qggqbFscPPMM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qggqbA6PPPM(ONE,TWO,THREE,FOUR,FIVE,
     >                                        SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qggqbA60PPPM,qggqbV,qggqbFccPPPM,qggqbFscPPPM

      qggqbA6PPPM = qggqbA60PPPM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >              *qggqbV(ONE,TWO,THREE,FOUR,FIVE,SIX)
     > + CMPLX(0D0,1D0)*qggqbFccPPPM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     > + CMPLX(0D0,1D0)*qggqbFscPPPM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qggqbA6PMMM(ONE,TWO,THREE,FOUR,FIVE,
     >                                        SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qggqbA60PMMM,qggqbV,qggqbFccPMMM,qggqbFscPMMM

      qggqbA6PMMM = qggqbA60PMMM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >              *qggqbV(ONE,TWO,THREE,FOUR,FIVE,SIX)
     > + CMPLX(0D0,1D0)*qggqbFccPMMM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     > + CMPLX(0D0,1D0)*qggqbFscPMMM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qgqbgA6PPMP(ONE,TWO,THREE,FOUR,FIVE,
     >                                        SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qgqbgA60PPMP,qgqbgV,qgqbgFccPPMP,qgqbgFscPPMP

      qgqbgA6PPMP = qgqbgA60PPMP(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >              *qgqbgV(ONE,TWO,THREE,FOUR,FIVE,SIX)
     > + CMPLX(0D0,1D0)*qgqbgFccPPMP(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     > + CMPLX(0D0,1D0)*qgqbgFscPPMP(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qgqbgA6PPMM(ONE,TWO,THREE,FOUR,FIVE,
     >                                        SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qgqbgA60PPMM,qgqbgV,qgqbgFccPPMM,qgqbgFscPPMM

      qgqbgA6PPMM = qgqbgA60PPMM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >              *qgqbgV(ONE,TWO,THREE,FOUR,FIVE,SIX)
     > + CMPLX(0D0,1D0)*qgqbgFccPPMM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     > + CMPLX(0D0,1D0)*qgqbgFscPPMM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qgqbgA6PMMP(ONE,TWO,THREE,FOUR,FIVE,
     >                                        SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qgqbgA60PMMP,qgqbgV,qgqbgFccPMMP,qgqbgFscPMMP

      qgqbgA6PMMP = qgqbgA60PMMP(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >              *qgqbgV(ONE,TWO,THREE,FOUR,FIVE,SIX)
     > + CMPLX(0D0,1D0)*qgqbgFccPMMP(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     > + CMPLX(0D0,1D0)*qgqbgFscPMMP(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qgqbgA6PMMM(ONE,TWO,THREE,FOUR,FIVE,
     >                                        SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qgqbgA60PMMM,qgqbgV,qgqbgFccPMMM,qgqbgFscPMMM

      qgqbgA6PMMM = qgqbgA60PMMM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >              *qgqbgV(ONE,TWO,THREE,FOUR,FIVE,SIX)
     > + CMPLX(0D0,1D0)*qgqbgFccPMMM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     > + CMPLX(0D0,1D0)*qgqbgFscPMMM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbggA6PMPP(ONE,TWO,THREE,FOUR,FIVE,
     >                                        SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qqbggA60PMPP,qqbggV,qqbggFccPMPP,qqbggFscPMPP

      qqbggA6PMPP = qqbggA60PMPP(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >              *qqbggV(ONE,TWO,THREE,FOUR,FIVE,SIX)
     > + CMPLX(0D0,1D0)*qqbggFccPMPP(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     > + CMPLX(0D0,1D0)*qqbggFscPMPP(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbggA6PMPM(ONE,TWO,THREE,FOUR,FIVE,
     >                                        SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qqbggA60PMPM,qqbggV,qqbggFccPMPM,qqbggFscPMPM

      qqbggA6PMPM = qqbggA60PMPM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >              *qqbggV(ONE,TWO,THREE,FOUR,FIVE,SIX)
     > + CMPLX(0D0,1D0)*qqbggFccPMPM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     > + CMPLX(0D0,1D0)*qqbggFscPMPM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbggA6PMMP(ONE,TWO,THREE,FOUR,FIVE,
     >                                        SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qqbggA60PMMP,qqbggV,qqbggFccPMMP,qqbggFscPMMP

      qqbggA6PMMP = qqbggA60PMMP(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >              *qqbggV(ONE,TWO,THREE,FOUR,FIVE,SIX)
     > + CMPLX(0D0,1D0)*qqbggFccPMMP(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     > + CMPLX(0D0,1D0)*qqbggFscPMMP(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbggA6PMMM(ONE,TWO,THREE,FOUR,FIVE,
     >                                        SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qqbggA60PMMM,qqbggV,qqbggFccPMMM,qqbggFscPMMM

      qqbggA6PMMM = qqbggA60PMMM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >              *qqbggV(ONE,TWO,THREE,FOUR,FIVE,SIX)
     > + CMPLX(0D0,1D0)*qqbggFccPMMM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     > + CMPLX(0D0,1D0)*qqbggFscPMMM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbggA6vsPMPP(ONE,TWO,THREE,FOUR,FIVE,
     >                                          SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qqbggFvsPMPP

      qqbggA6vsPMPP =
     >   CMPLX(0D0,1D0)*qqbggFvsPMPP(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbggA6vsPMPM(ONE,TWO,THREE,FOUR,FIVE,
     >                                          SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qqbggFvsPMPM

      qqbggA6vsPMPM =
     >   CMPLX(0D0,1D0)*qqbggFvsPMPM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbggA6vsPMMP(ONE,TWO,THREE,FOUR,FIVE,
     >                                          SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qqbggFvsPMMP

      qqbggA6vsPMMP =
     >   CMPLX(0D0,1D0)*qqbggFvsPMMP(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbggA6vsPMMM(ONE,TWO,THREE,FOUR,FIVE,
     >                                          SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qqbggFvsPMMM

      qqbggA6vsPMMM =
     >   CMPLX(0D0,1D0)*qqbggFvsPMMM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbggA6vfPMPP(ONE,TWO,THREE,FOUR,FIVE,
     >                                          SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qqbggFvfPMPP

      qqbggA6vfPMPP =
     >   CMPLX(0D0,1D0)*qqbggFvfPMPP(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbggA6vfPMPM(ONE,TWO,THREE,FOUR,FIVE,
     >                                          SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qqbggFvfPMPM

      qqbggA6vfPMPM =
     >   CMPLX(0D0,1D0)*qqbggFvfPMPM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbggA6vfPMMP(ONE,TWO,THREE,FOUR,FIVE,
     >                                          SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qqbggFvfPMMP

      qqbggA6vfPMMP =
     >   CMPLX(0D0,1D0)*qqbggFvfPMMP(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbggA6vfPMMM(ONE,TWO,THREE,FOUR,FIVE,
     >                                          SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qqbggFvfPMMM

      qqbggA6vfPMMM =
     >   CMPLX(0D0,1D0)*qqbggFvfPMMM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbggA6axPMPP(ONE,TWO,THREE,FOUR,FIVE,
     >                                          SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qqbggFaxPMPP

      qqbggA6axPMPP =
     >   CMPLX(0D0,1D0)*qqbggFaxPMPP(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbggA6axPMPM(ONE,TWO,THREE,FOUR,FIVE,
     >                                          SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qqbggFaxPMPM

      qqbggA6axPMPM =
     >   CMPLX(0D0,1D0)*qqbggFaxPMPM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbggA6axPMMP(ONE,TWO,THREE,FOUR,FIVE,
     >                                          SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qqbggFaxPMMP

      qqbggA6axPMMP =
     >   CMPLX(0D0,1D0)*qqbggFaxPMMP(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbggA6axPMMM(ONE,TWO,THREE,FOUR,FIVE,
     >                                          SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qqbggFaxPMMM

      qqbggA6axPMMM =
     >   CMPLX(0D0,1D0)*qqbggFaxPMMM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION
     > qqbggA6axslPMPP(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qqbggFaxslPMPP

      qqbggA6axslPMPP =
     >   CMPLX(0D0,1D0)*qqbggFaxslPMPP(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION
     > qqbggA6axslPMPM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qqbggFaxslPMPM

      qqbggA6axslPMPM =
     >   CMPLX(0D0,1D0)*qqbggFaxslPMPM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION
     > qqbggA6axslPMMP(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qqbggFaxslPMMP

      qqbggA6axslPMMP =
     >   CMPLX(0D0,1D0)*qqbggFaxslPMMP(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION
     > qqbggA6axslPMMM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qqbggFaxslPMMM

      qqbggA6axslPMMM =
     >   CMPLX(0D0,1D0)*qqbggFaxslPMMM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)

      END

************************************************************************
************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qggqbV(ONE,TWO,THREE,FOUR,FIVE,SIX)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      REAL(KIND(1D0)) S(11,11)
      REAL(KIND(1D0)) MUSQ,COMU
      COMMON /DOTPRODUCTS/ S
     >       /scales/ COMU
      COMPLEX(KIND(1D0)) LOGS12,LOGS23,LOGS34,LOGS56
      REAL(KIND(1D0)) THETA,PI
      PARAMETER ( PI  = 3.14159265358979 D0)

      MUSQ = COMU**2
      LOGS12 = CMPLX(LOG(ABS(S(ONE,TWO))/MUSQ),
     >               -PI*THETA(S(ONE,TWO).GT.0D0))

      LOGS23 = CMPLX(LOG(ABS(S(TWO,THREE))/MUSQ),
     >               -PI*THETA(S(TWO,THREE).GT.0D0))

      LOGS34 = CMPLX(LOG(ABS(S(THREE,FOUR))/MUSQ),
     >               -PI*THETA(S(THREE,FOUR).GT.0D0))

      LOGS56 = CMPLX(LOG(ABS(S(FIVE,SIX))/MUSQ),
     >               -PI*THETA(S(FIVE,SIX).GT.0D0))
      
C The e^-2 part

c      qggqbV = -3D0

C The e^-1 part

c      qggqbV = 1d0*(LOGS12 + LOGS23 + LOGS34) - 1d0*1.5d0

C The e^0 part

      qggqbV = -0.5D0*LOGS12**2 - 0.5*LOGS23**2 - 0.5D0*LOGS34**2
     >         +1.5D0*LOGS56 - 3.5D0
c      qggqbV  = -3.5D0

c      qgV = 1.5D0*LOGS56 - 3.5D0

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qgqbgV(ONE,TWO,THREE,FOUR,FIVE,SIX)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      REAL(KIND(1D0)) S(11,11)
      REAL(KIND(1D0)) MUSQ,COMU
      COMMON /DOTPRODUCTS/ S
     >       /scales/ COMU
      COMPLEX(KIND(1D0)) LOGS12,LOGS23,LOGS56
      REAL(KIND(1D0)) THETA,PI
      PARAMETER ( PI  = 3.14159265358979 D0)

      MUSQ = COMU**2
      LOGS12 = CMPLX(LOG(ABS(S(ONE,TWO))/MUSQ),
     >               -PI*THETA(S(ONE,TWO).GT.0D0))

      LOGS23 = CMPLX(LOG(ABS(S(TWO,THREE))/MUSQ),
     >               -PI*THETA(S(TWO,THREE).GT.0D0))

      LOGS56 = CMPLX(LOG(ABS(S(FIVE,SIX))/MUSQ),
     >               -PI*THETA(S(FIVE,SIX).GT.0D0))
      
C The e^-2 part

c      qgqbgV = -2D0

C The e^-1 part

c      qgqbgV = 1d0*(LOGS12 + LOGS23) - 1d0*1.5d0

C The e^0 part

      qgqbgV = -0.5D0*LOGS12**2 - 0.5*LOGS23**2
     >         +1.5D0*LOGS56 - 3.5D0
c      qgqbgV = -3.5D0

c      qgV = 1.5D0*LOGS56 - 3.5D0

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbggV(ONE,TWO,THREE,FOUR,FIVE,SIX)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      REAL(KIND(1D0)) S(11,11)
      REAL(KIND(1D0)) MUSQ,COMU
      COMMON /DOTPRODUCTS/ S
     >       /scales/ COMU
      COMPLEX(KIND(1D0)) LOGS12,LOGS56
      REAL(KIND(1D0)) THETA,PI
      PARAMETER ( PI  = 3.14159265358979 D0)

      MUSQ = COMU**2
      LOGS12 = CMPLX(LOG(ABS(S(ONE,TWO))/MUSQ),
     >               -PI*THETA(S(ONE,TWO).GT.0D0))

      LOGS56 = CMPLX(LOG(ABS(S(FIVE,SIX))/MUSQ),
     >               -PI*THETA(S(FIVE,SIX).GT.0D0))
      
C The e^-2 part

c      qqbggV = -1D0

C The e^-1 part

c      qqbggV = 1d0*LOGS12 - 1d0*1.5d0

C The e^0 part

      qqbggV = -0.5D0*LOGS12**2
     >         +1.5D0*LOGS56 - 3.5D0
c      qqbggV = -3.5D0

c      qgV = 1.5D0*LOGS56 - 3.5D0

      END

************************************************************************
c The next part contains routines to calculate the pole parts of the
c virtual contribution:
************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qggqbA61PMPMem2(ONE,TWO,THREE,FOUR,
     >                                            FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qggqbA6PMPMem2,qqbggA6PMPMem2
      REAL(KIND(1D0)) XNf,XNu,XNd,XNc,XNa                      
      COMMON /QCD/ XNf,XNu,XNd,XNc,XNa

      qggqbA61PMPMem2 = qggqbA6PMPMem2(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >     - 1D0/XNc**2*qqbggA6PMPMem2(ONE,FOUR,THREE,TWO,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qggqbA61PPMMem2(ONE,TWO,THREE,FOUR,
     >                                            FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qggqbA6PPMMem2,qqbggA6PMMPem2
      REAL(KIND(1D0)) XNf,XNu,XNd,XNc,XNa                      
      COMMON /QCD/ XNf,XNu,XNd,XNc,XNa

      qggqbA61PPMMem2 = qggqbA6PPMMem2(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >     - 1D0/XNc**2*qqbggA6PMMPem2(ONE,FOUR,THREE,TWO,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qggqbA61PPPMem2(ONE,TWO,THREE,FOUR,
     >                                            FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qggqbA6PPPMem2,qqbggA6PMPPem2
      REAL(KIND(1D0)) XNf,XNu,XNd,XNc,XNa                      
      COMMON /QCD/ XNf,XNu,XNd,XNc,XNa

      qggqbA61PPPMem2 = qggqbA6PPPMem2(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >     - 1D0/XNc**2*qqbggA6PMPPem2(ONE,FOUR,THREE,TWO,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qggqbA61PMMMem2(ONE,TWO,THREE,FOUR,
     >                                            FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qggqbA6PMMMem2,qqbggA6PMMMem2
      REAL(KIND(1D0)) XNf,XNu,XNd,XNc,XNa                      
      COMMON /QCD/ XNf,XNu,XNd,XNc,XNa

      qggqbA61PMMMem2 = qggqbA6PMMMem2(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >     - 1D0/XNc**2*qqbggA6PMMMem2(ONE,FOUR,THREE,TWO,FIVE,SIX,A,B)

      END

************************************************************************
************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbggA63PMPMem2(ONE,FOUR,TWO,THREE,
     >                                            FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,FOUR,TWO,THREE,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qggqbA6PPMMem2,qggqbA6PMPMem2,qgqbgA6PPMMem2,
     >                   qgqbgA6PMMPem2,qqbggA6PMPMem2,qqbggA6PMMPem2

      qqbggA63PMPMem2 = qggqbA6PPMMem2(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >                + qggqbA6PMPMem2(ONE,THREE,TWO,FOUR,FIVE,SIX,A,B)
     >                + qgqbgA6PPMMem2(ONE,TWO,FOUR,THREE,FIVE,SIX,A,B)
     >                + qgqbgA6PMMPem2(ONE,THREE,FOUR,TWO,FIVE,SIX,A,B)
     >                + qqbggA6PMPMem2(ONE,FOUR,TWO,THREE,FIVE,SIX,A,B)
     >                + qqbggA6PMMPem2(ONE,FOUR,THREE,TWO,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbggA63PMMPem2(ONE,FOUR,TWO,THREE,
     >                                            FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,FOUR,TWO,THREE,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qggqbA6PPMMem2,qggqbA6PMPMem2,qgqbgA6PPMMem2,
     >                   qgqbgA6PMMPem2,qqbggA6PMPMem2,qqbggA6PMMPem2

      qqbggA63PMMPem2 = qggqbA6PMPMem2(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >                + qggqbA6PPMMem2(ONE,THREE,TWO,FOUR,FIVE,SIX,A,B)
     >                + qgqbgA6PMMPem2(ONE,TWO,FOUR,THREE,FIVE,SIX,A,B)
     >                + qgqbgA6PPMMem2(ONE,THREE,FOUR,TWO,FIVE,SIX,A,B)
     >                + qqbggA6PMMPem2(ONE,FOUR,TWO,THREE,FIVE,SIX,A,B)
     >                + qqbggA6PMPMem2(ONE,FOUR,THREE,TWO,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbggA63PMPPem2(ONE,FOUR,TWO,THREE,
     >                                            FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,FOUR,TWO,THREE,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qggqbA6PPPMem2,qgqbgA6PPMPem2,qqbggA6PMPPem2

      qqbggA63PMPPem2 = qggqbA6PPPMem2(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >                + qggqbA6PPPMem2(ONE,THREE,TWO,FOUR,FIVE,SIX,A,B)
     >                + qgqbgA6PPMPem2(ONE,TWO,FOUR,THREE,FIVE,SIX,A,B)
     >                + qgqbgA6PPMPem2(ONE,THREE,FOUR,TWO,FIVE,SIX,A,B)
     >                + qqbggA6PMPPem2(ONE,FOUR,TWO,THREE,FIVE,SIX,A,B)
     >                + qqbggA6PMPPem2(ONE,FOUR,THREE,TWO,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbggA63PMMMem2(ONE,FOUR,TWO,THREE,
     >                                            FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,FOUR,TWO,THREE,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qggqbA6PMMMem2,qgqbgA6PMMMem2,qqbggA6PMMMem2

      qqbggA63PMMMem2 = qggqbA6PMMMem2(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >                + qggqbA6PMMMem2(ONE,THREE,TWO,FOUR,FIVE,SIX,A,B)
     >                + qgqbgA6PMMMem2(ONE,TWO,FOUR,THREE,FIVE,SIX,A,B)
     >                + qgqbgA6PMMMem2(ONE,THREE,FOUR,TWO,FIVE,SIX,A,B)
     >                + qqbggA6PMMMem2(ONE,FOUR,TWO,THREE,FIVE,SIX,A,B)
     >                + qqbggA6PMMMem2(ONE,FOUR,THREE,TWO,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qggqbA61PMPMem1(ONE,TWO,THREE,FOUR,
     >                                            FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qggqbA6PMPMem1,qqbggA6PMPMem1
      REAL(KIND(1D0)) XNf,XNu,XNd,XNc,XNa                      
      COMMON /QCD/ XNf,XNu,XNd,XNc,XNa

      qggqbA61PMPMem1 = qggqbA6PMPMem1(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >     - 1D0/XNc**2*qqbggA6PMPMem1(ONE,FOUR,THREE,TWO,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qggqbA61PPMMem1(ONE,TWO,THREE,FOUR,
     >                                            FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qggqbA6PPMMem1,qqbggA6PMMPem1
      REAL(KIND(1D0)) XNf,XNu,XNd,XNc,XNa                      
      COMMON /QCD/ XNf,XNu,XNd,XNc,XNa

      qggqbA61PPMMem1 = qggqbA6PPMMem1(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >     - 1D0/XNc**2*qqbggA6PMMPem1(ONE,FOUR,THREE,TWO,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qggqbA61PPPMem1(ONE,TWO,THREE,FOUR,
     >                                            FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qggqbA6PPPMem1,qqbggA6PMPPem1
      REAL(KIND(1D0)) XNf,XNu,XNd,XNc,XNa                      
      COMMON /QCD/ XNf,XNu,XNd,XNc,XNa

      qggqbA61PPPMem1 = qggqbA6PPPMem1(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >     - 1D0/XNc**2*qqbggA6PMPPem1(ONE,FOUR,THREE,TWO,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qggqbA61PMMMem1(ONE,TWO,THREE,FOUR,
     >                                            FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qggqbA6PMMMem1,qqbggA6PMMMem1
      REAL(KIND(1D0)) XNf,XNu,XNd,XNc,XNa                      
      COMMON /QCD/ XNf,XNu,XNd,XNc,XNa

      qggqbA61PMMMem1 = qggqbA6PMMMem1(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >     - 1D0/XNc**2*qqbggA6PMMMem1(ONE,FOUR,THREE,TWO,FIVE,SIX,A,B)

      END

************************************************************************
************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbggA63PMPMem1(ONE,FOUR,TWO,THREE,
     >                                            FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,FOUR,TWO,THREE,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qggqbA6PPMMem1,qggqbA6PMPMem1,qgqbgA6PPMMem1,
     >                   qgqbgA6PMMPem1,qqbggA6PMPMem1,qqbggA6PMMPem1

      qqbggA63PMPMem1 = qggqbA6PPMMem1(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >                + qggqbA6PMPMem1(ONE,THREE,TWO,FOUR,FIVE,SIX,A,B)
     >                + qgqbgA6PPMMem1(ONE,TWO,FOUR,THREE,FIVE,SIX,A,B)
     >                + qgqbgA6PMMPem1(ONE,THREE,FOUR,TWO,FIVE,SIX,A,B)
     >                + qqbggA6PMPMem1(ONE,FOUR,TWO,THREE,FIVE,SIX,A,B)
     >                + qqbggA6PMMPem1(ONE,FOUR,THREE,TWO,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbggA63PMMPem1(ONE,FOUR,TWO,THREE,
     >                                            FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,FOUR,TWO,THREE,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qggqbA6PPMMem1,qggqbA6PMPMem1,qgqbgA6PPMMem1,
     >                   qgqbgA6PMMPem1,qqbggA6PMPMem1,qqbggA6PMMPem1

      qqbggA63PMMPem1 = qggqbA6PMPMem1(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >                + qggqbA6PPMMem1(ONE,THREE,TWO,FOUR,FIVE,SIX,A,B)
     >                + qgqbgA6PMMPem1(ONE,TWO,FOUR,THREE,FIVE,SIX,A,B)
     >                + qgqbgA6PPMMem1(ONE,THREE,FOUR,TWO,FIVE,SIX,A,B)
     >                + qqbggA6PMMPem1(ONE,FOUR,TWO,THREE,FIVE,SIX,A,B)
     >                + qqbggA6PMPMem1(ONE,FOUR,THREE,TWO,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbggA63PMPPem1(ONE,FOUR,TWO,THREE,
     >                                            FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,FOUR,TWO,THREE,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qggqbA6PPPMem1,qgqbgA6PPMPem1,qqbggA6PMPPem1

      qqbggA63PMPPem1 = qggqbA6PPPMem1(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >                + qggqbA6PPPMem1(ONE,THREE,TWO,FOUR,FIVE,SIX,A,B)
     >                + qgqbgA6PPMPem1(ONE,TWO,FOUR,THREE,FIVE,SIX,A,B)
     >                + qgqbgA6PPMPem1(ONE,THREE,FOUR,TWO,FIVE,SIX,A,B)
     >                + qqbggA6PMPPem1(ONE,FOUR,TWO,THREE,FIVE,SIX,A,B)
     >                + qqbggA6PMPPem1(ONE,FOUR,THREE,TWO,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbggA63PMMMem1(ONE,FOUR,TWO,THREE,
     >                                            FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,FOUR,TWO,THREE,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qggqbA6PMMMem1,qgqbgA6PMMMem1,qqbggA6PMMMem1

      qqbggA63PMMMem1 = qggqbA6PMMMem1(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >                + qggqbA6PMMMem1(ONE,THREE,TWO,FOUR,FIVE,SIX,A,B)
     >                + qgqbgA6PMMMem1(ONE,TWO,FOUR,THREE,FIVE,SIX,A,B)
     >                + qgqbgA6PMMMem1(ONE,THREE,FOUR,TWO,FIVE,SIX,A,B)
     >                + qqbggA6PMMMem1(ONE,FOUR,TWO,THREE,FIVE,SIX,A,B)
     >                + qqbggA6PMMMem1(ONE,FOUR,THREE,TWO,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qggqbA6PMPMem2(ONE,TWO,THREE,FOUR,
     >                                           FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qggqbA60PMPM,qggqbVem2

      qggqbA6PMPMem2 = qggqbA60PMPM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >               * qggqbVem2(ONE,TWO,THREE,FOUR,FIVE,SIX)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qggqbA6PPMMem2(ONE,TWO,THREE,FOUR,
     >                                           FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qggqbA60PPMM,qggqbVem2

      qggqbA6PPMMem2 = qggqbA60PPMM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >               * qggqbVem2(ONE,TWO,THREE,FOUR,FIVE,SIX)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qggqbA6PPPMem2(ONE,TWO,THREE,FOUR,
     >                                           FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qggqbA60PPPM,qggqbVem2

      qggqbA6PPPMem2 = qggqbA60PPPM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >               * qggqbVem2(ONE,TWO,THREE,FOUR,FIVE,SIX)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qggqbA6PMMMem2(ONE,TWO,THREE,FOUR,
     >                                           FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qggqbA60PMMM,qggqbVem2

      qggqbA6PMMMem2 = qggqbA60PMMM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >               * qggqbVem2(ONE,TWO,THREE,FOUR,FIVE,SIX)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qgqbgA6PPMPem2(ONE,TWO,THREE,FOUR,
     >                                           FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qgqbgA60PPMP,qgqbgVem2

      qgqbgA6PPMPem2 = qgqbgA60PPMP(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >               * qgqbgVem2(ONE,TWO,THREE,FOUR,FIVE,SIX)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qgqbgA6PPMMem2(ONE,TWO,THREE,FOUR,
     >                                           FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qgqbgA60PPMM,qgqbgVem2

      qgqbgA6PPMMem2 = qgqbgA60PPMM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >               * qgqbgVem2(ONE,TWO,THREE,FOUR,FIVE,SIX)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qgqbgA6PMMPem2(ONE,TWO,THREE,FOUR,
     >                                           FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qgqbgA60PMMP,qgqbgVem2

      qgqbgA6PMMPem2 = qgqbgA60PMMP(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >               * qgqbgVem2(ONE,TWO,THREE,FOUR,FIVE,SIX)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qgqbgA6PMMMem2(ONE,TWO,THREE,FOUR,
     >                                           FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qgqbgA60PMMM,qgqbgVem2

      qgqbgA6PMMMem2 = qgqbgA60PMMM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >               * qgqbgVem2(ONE,TWO,THREE,FOUR,FIVE,SIX)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbggA6PMPPem2(ONE,TWO,THREE,FOUR,
     >                                           FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qqbggA60PMPP,qqbggVem2

      qqbggA6PMPPem2 = qqbggA60PMPP(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >               * qqbggVem2(ONE,TWO,THREE,FOUR,FIVE,SIX)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbggA6PMPMem2(ONE,TWO,THREE,FOUR,
     >                                           FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qqbggA60PMPM,qqbggVem2

      qqbggA6PMPMem2 = qqbggA60PMPM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >               * qqbggVem2(ONE,TWO,THREE,FOUR,FIVE,SIX)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbggA6PMMPem2(ONE,TWO,THREE,FOUR,
     >                                           FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qqbggA60PMMP,qqbggVem2

      qqbggA6PMMPem2 = qqbggA60PMMP(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >               * qqbggVem2(ONE,TWO,THREE,FOUR,FIVE,SIX)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbggA6PMMMem2(ONE,TWO,THREE,FOUR,
     >                                           FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qqbggA60PMMM,qqbggVem2

      qqbggA6PMMMem2 = qqbggA60PMMM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >               * qqbggVem2(ONE,TWO,THREE,FOUR,FIVE,SIX)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qggqbA6PMPMem1(ONE,TWO,THREE,FOUR,
     >                                           FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qggqbA60PMPM,qggqbVem1

      qggqbA6PMPMem1 = qggqbA60PMPM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >               * qggqbVem1(ONE,TWO,THREE,FOUR,FIVE,SIX)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qggqbA6PPMMem1(ONE,TWO,THREE,FOUR,
     >                                           FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qggqbA60PPMM,qggqbVem1

      qggqbA6PPMMem1 = qggqbA60PPMM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >               * qggqbVem1(ONE,TWO,THREE,FOUR,FIVE,SIX)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qggqbA6PPPMem1(ONE,TWO,THREE,FOUR,
     >                                           FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qggqbA60PPPM,qggqbVem1

      qggqbA6PPPMem1 = qggqbA60PPPM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >               * qggqbVem1(ONE,TWO,THREE,FOUR,FIVE,SIX)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qggqbA6PMMMem1(ONE,TWO,THREE,FOUR,
     >                                           FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qggqbA60PMMM,qggqbVem1

      qggqbA6PMMMem1 = qggqbA60PMMM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >               * qggqbVem1(ONE,TWO,THREE,FOUR,FIVE,SIX)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qgqbgA6PPMPem1(ONE,TWO,THREE,FOUR,
     >                                           FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qgqbgA60PPMP,qgqbgVem1

      qgqbgA6PPMPem1 = qgqbgA60PPMP(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >               * qgqbgVem1(ONE,TWO,THREE,FOUR,FIVE,SIX)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qgqbgA6PPMMem1(ONE,TWO,THREE,FOUR,
     >                                           FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qgqbgA60PPMM,qgqbgVem1

      qgqbgA6PPMMem1 = qgqbgA60PPMM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >               * qgqbgVem1(ONE,TWO,THREE,FOUR,FIVE,SIX)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qgqbgA6PMMPem1(ONE,TWO,THREE,FOUR,
     >                                           FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qgqbgA60PMMP,qgqbgVem1

      qgqbgA6PMMPem1 = qgqbgA60PMMP(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >               * qgqbgVem1(ONE,TWO,THREE,FOUR,FIVE,SIX)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qgqbgA6PMMMem1(ONE,TWO,THREE,FOUR,
     >                                           FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qgqbgA60PMMM,qgqbgVem1

      qgqbgA6PMMMem1 = qgqbgA60PMMM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >               * qgqbgVem1(ONE,TWO,THREE,FOUR,FIVE,SIX)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbggA6PMPPem1(ONE,TWO,THREE,FOUR,
     >                                           FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qqbggA60PMPP,qqbggVem1

      qqbggA6PMPPem1 = qqbggA60PMPP(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >               * qqbggVem1(ONE,TWO,THREE,FOUR,FIVE,SIX)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbggA6PMPMem1(ONE,TWO,THREE,FOUR,
     >                                           FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qqbggA60PMPM,qqbggVem1

      qqbggA6PMPMem1 = qqbggA60PMPM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >               * qqbggVem1(ONE,TWO,THREE,FOUR,FIVE,SIX)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbggA6PMMPem1(ONE,TWO,THREE,FOUR,
     >                                           FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qqbggA60PMMP,qqbggVem1

      qqbggA6PMMPem1 = qqbggA60PMMP(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >               * qqbggVem1(ONE,TWO,THREE,FOUR,FIVE,SIX)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbggA6PMMMem1(ONE,TWO,THREE,FOUR,
     >                                           FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qqbggA60PMMM,qqbggVem1

      qqbggA6PMMMem1 = qqbggA60PMMM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >               * qqbggVem1(ONE,TWO,THREE,FOUR,FIVE,SIX)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qggqbVem2(ONE,TWO,THREE,FOUR,FIVE,SIX)
      IMPLICIT NONE
      
      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      REAL(KIND(1D0)) S(11,11)
      REAL(KIND(1D0)) MUSQ,COMU
      COMMON /DOTPRODUCTS/ S
     >       /scales/ COMU
      COMPLEX(KIND(1D0)) LOGS12,LOGS23,LOGS34,LOGS56
      REAL(KIND(1D0)) THETA,PI
      PARAMETER ( PI  = 3.14159265358979 D0)

      MUSQ = COMU**2
      LOGS12 = CMPLX(LOG(ABS(S(ONE,TWO))/MUSQ),
     >               -PI*THETA(S(ONE,TWO).GT.0D0))

      LOGS23 = CMPLX(LOG(ABS(S(TWO,THREE))/MUSQ),
     >               -PI*THETA(S(TWO,THREE).GT.0D0))

      LOGS34 = CMPLX(LOG(ABS(S(THREE,FOUR))/MUSQ),
     >               -PI*THETA(S(THREE,FOUR).GT.0D0))

      LOGS56 = CMPLX(LOG(ABS(S(FIVE,SIX))/MUSQ),
     >               -PI*THETA(S(FIVE,SIX).GT.0D0))
      
C The e^-2 part

      qggqbVem2 = -3D0

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qgqbgVem2(ONE,TWO,THREE,FOUR,FIVE,SIX)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      REAL(KIND(1D0)) S(11,11)
      REAL(KIND(1D0)) MUSQ,COMU
      COMMON /DOTPRODUCTS/ S
     >       /scales/ COMU
      COMPLEX(KIND(1D0)) LOGS12,LOGS23,LOGS56
      REAL(KIND(1D0)) THETA,PI
      PARAMETER ( PI  = 3.14159265358979 D0)

      MUSQ = COMU**2
      LOGS12 = CMPLX(LOG(ABS(S(ONE,TWO))/MUSQ),
     >               -PI*THETA(S(ONE,TWO).GT.0D0))

      LOGS23 = CMPLX(LOG(ABS(S(TWO,THREE))/MUSQ),
     >               -PI*THETA(S(TWO,THREE).GT.0D0))

      LOGS56 = CMPLX(LOG(ABS(S(FIVE,SIX))/MUSQ),
     >               -PI*THETA(S(FIVE,SIX).GT.0D0))
      
C The e^-2 part

      qgqbgVem2 = -2D0

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbggVem2(ONE,TWO,THREE,FOUR,FIVE,SIX)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      REAL(KIND(1D0)) S(11,11)
      REAL(KIND(1D0)) MUSQ,COMU
      COMMON /DOTPRODUCTS/ S
     >       /scales/ COMU
      COMPLEX(KIND(1D0)) LOGS12,LOGS56
      REAL(KIND(1D0)) THETA,PI
      PARAMETER ( PI  = 3.14159265358979 D0)

      MUSQ = COMU**2
      LOGS12 = CMPLX(LOG(ABS(S(ONE,TWO))/MUSQ),
     >               -PI*THETA(S(ONE,TWO).GT.0D0))

      LOGS56 = CMPLX(LOG(ABS(S(FIVE,SIX))/MUSQ),
     >               -PI*THETA(S(FIVE,SIX).GT.0D0))
      
C The e^-2 part

      qqbggVem2 = -1D0

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qggqbVem1(ONE,TWO,THREE,FOUR,FIVE,SIX)
      IMPLICIT NONE
      
      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      REAL(KIND(1D0)) S(11,11)
      REAL(KIND(1D0)) MUSQ,COMU
      COMMON /DOTPRODUCTS/ S
     >       /scales/ COMU
      COMPLEX(KIND(1D0)) LOGS12,LOGS23,LOGS34,LOGS56
      REAL(KIND(1D0)) THETA,PI
      PARAMETER ( PI  = 3.14159265358979 D0)

      MUSQ = COMU**2
      LOGS12 = CMPLX(LOG(ABS(S(ONE,TWO))/MUSQ),
     >               -PI*THETA(S(ONE,TWO).GT.0D0))

      LOGS23 = CMPLX(LOG(ABS(S(TWO,THREE))/MUSQ),
     >               -PI*THETA(S(TWO,THREE).GT.0D0))

      LOGS34 = CMPLX(LOG(ABS(S(THREE,FOUR))/MUSQ),
     >               -PI*THETA(S(THREE,FOUR).GT.0D0))

      LOGS56 = CMPLX(LOG(ABS(S(FIVE,SIX))/MUSQ),
     >               -PI*THETA(S(FIVE,SIX).GT.0D0))
      
C The e^-1 part

      qggqbVem1 = 1d0*(LOGS12 + LOGS23 + LOGS34) - 1d0*1.5d0

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qgqbgVem1(ONE,TWO,THREE,FOUR,FIVE,SIX)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      REAL(KIND(1D0)) S(11,11)
      REAL(KIND(1D0)) MUSQ,COMU
      COMMON /DOTPRODUCTS/ S
     >       /scales/ COMU
      COMPLEX(KIND(1D0)) LOGS12,LOGS23,LOGS56
      REAL(KIND(1D0)) THETA,PI
      PARAMETER ( PI  = 3.14159265358979 D0)

      MUSQ = COMU**2
      LOGS12 = CMPLX(LOG(ABS(S(ONE,TWO))/MUSQ),
     >               -PI*THETA(S(ONE,TWO).GT.0D0))

      LOGS23 = CMPLX(LOG(ABS(S(TWO,THREE))/MUSQ),
     >               -PI*THETA(S(TWO,THREE).GT.0D0))

      LOGS56 = CMPLX(LOG(ABS(S(FIVE,SIX))/MUSQ),
     >               -PI*THETA(S(FIVE,SIX).GT.0D0))
      
C The e^-1 part

      qgqbgVem1 = 1d0*(LOGS12 + LOGS23) - 1d0*1.5d0

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbggVem1(ONE,TWO,THREE,FOUR,FIVE,SIX)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      REAL(KIND(1D0)) S(11,11)
      REAL(KIND(1D0)) MUSQ,COMU
      COMMON /DOTPRODUCTS/ S
     >       /scales/ COMU
      COMPLEX(KIND(1D0)) LOGS12,LOGS56
      REAL(KIND(1D0)) THETA,PI
      PARAMETER ( PI  = 3.14159265358979 D0)

      MUSQ = COMU**2
      LOGS12 = CMPLX(LOG(ABS(S(ONE,TWO))/MUSQ),
     >               -PI*THETA(S(ONE,TWO).GT.0D0))

      LOGS56 = CMPLX(LOG(ABS(S(FIVE,SIX))/MUSQ),
     >               -PI*THETA(S(FIVE,SIX).GT.0D0))
      
C The e^-1 part

      qqbggVem1 = 1d0*LOGS12 - 1d0*1.5d0

      END

************************************************************************
c The next part comprises routines decomposed into a muR dependent and
c an independent part:
************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qggqbA6PMPMmudep(ONE,TWO,THREE,
     >                                             FOUR,FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qggqbA60PMPM,qggqbV,qggqbFccPMPM,qggqbFscPMPM

      qggqbA6PMPMmudep = qggqbA60PMPM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >                 * qggqbV(ONE,TWO,THREE,FOUR,FIVE,SIX)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qggqbA6PMPMmuindep(ONE,TWO,THREE,
     >                                               FOUR,FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qggqbA60PMPM,qggqbV,qggqbFccPMPM,qggqbFscPMPM

      qggqbA6PMPMmuindep = 
     > + CMPLX(0D0,1D0)*qggqbFccPMPM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     > + CMPLX(0D0,1D0)*qggqbFscPMPM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qggqbA6PPMMmudep(ONE,TWO,THREE,
     >                                             FOUR,FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qggqbA60PPMM,qggqbV,qggqbFccPPMM,qggqbFscPPMM

      qggqbA6PPMMmudep = qggqbA60PPMM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >                 * qggqbV(ONE,TWO,THREE,FOUR,FIVE,SIX)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qggqbA6PPMMmuindep(ONE,TWO,THREE,
     >                                               FOUR,FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qggqbA60PPMM,qggqbV,qggqbFccPPMM,qggqbFscPPMM

      qggqbA6PPMMmuindep = 
     > + CMPLX(0D0,1D0)*qggqbFccPPMM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     > + CMPLX(0D0,1D0)*qggqbFscPPMM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qggqbA6PPPMmudep(ONE,TWO,THREE,
     >                                             FOUR,FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qggqbA60PPPM,qggqbV,qggqbFccPPPM,qggqbFscPPPM

      qggqbA6PPPMmudep = qggqbA60PPPM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >                 * qggqbV(ONE,TWO,THREE,FOUR,FIVE,SIX)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qggqbA6PPPMmuindep(ONE,TWO,THREE,
     >                                               FOUR,FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qggqbA60PPPM,qggqbV,qggqbFccPPPM,qggqbFscPPPM

      qggqbA6PPPMmuindep = 
     > + CMPLX(0D0,1D0)*qggqbFccPPPM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     > + CMPLX(0D0,1D0)*qggqbFscPPPM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qggqbA6PMMMmudep(ONE,TWO,THREE,
     >                                             FOUR,FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qggqbA60PMMM,qggqbV,qggqbFccPMMM,qggqbFscPMMM

      qggqbA6PMMMmudep = qggqbA60PMMM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >                 * qggqbV(ONE,TWO,THREE,FOUR,FIVE,SIX)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qggqbA6PMMMmuindep(ONE,TWO,THREE,
     >                                               FOUR,FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qggqbA60PMMM,qggqbV,qggqbFccPMMM,qggqbFscPMMM

      qggqbA6PMMMmuindep = 
     > + CMPLX(0D0,1D0)*qggqbFccPMMM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     > + CMPLX(0D0,1D0)*qggqbFscPMMM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qgqbgA6PPMPmudep(ONE,TWO,THREE,
     >                                             FOUR,FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qgqbgA60PPMP,qgqbgV,qgqbgFccPPMP,qgqbgFscPPMP

      qgqbgA6PPMPmudep = qgqbgA60PPMP(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >                 * qgqbgV(ONE,TWO,THREE,FOUR,FIVE,SIX)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qgqbgA6PPMPmuindep(ONE,TWO,THREE,
     >                                               FOUR,FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qgqbgA60PPMP,qgqbgV,qgqbgFccPPMP,qgqbgFscPPMP

      qgqbgA6PPMPmuindep = 
     > + CMPLX(0D0,1D0)*qgqbgFccPPMP(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     > + CMPLX(0D0,1D0)*qgqbgFscPPMP(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qgqbgA6PPMMmudep(ONE,TWO,THREE,
     >                                             FOUR,FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qgqbgA60PPMM,qgqbgV,qgqbgFccPPMM,qgqbgFscPPMM

      qgqbgA6PPMMmudep = qgqbgA60PPMM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >                 * qgqbgV(ONE,TWO,THREE,FOUR,FIVE,SIX)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qgqbgA6PPMMmuindep(ONE,TWO,THREE,
     >                                               FOUR,FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qgqbgA60PPMM,qgqbgV,qgqbgFccPPMM,qgqbgFscPPMM

      qgqbgA6PPMMmuindep = 
     > + CMPLX(0D0,1D0)*qgqbgFccPPMM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     > + CMPLX(0D0,1D0)*qgqbgFscPPMM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qgqbgA6PMMPmudep(ONE,TWO,THREE,
     >                                             FOUR,FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qgqbgA60PMMP,qgqbgV,qgqbgFccPMMP,qgqbgFscPMMP

      qgqbgA6PMMPmudep = qgqbgA60PMMP(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >                 * qgqbgV(ONE,TWO,THREE,FOUR,FIVE,SIX)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qgqbgA6PMMPmuindep(ONE,TWO,THREE,
     >                                               FOUR,FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qgqbgA60PMMP,qgqbgV,qgqbgFccPMMP,qgqbgFscPMMP

      qgqbgA6PMMPmuindep = 
     > + CMPLX(0D0,1D0)*qgqbgFccPMMP(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     > + CMPLX(0D0,1D0)*qgqbgFscPMMP(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qgqbgA6PMMMmudep(ONE,TWO,THREE,
     >                                             FOUR,FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qgqbgA60PMMM,qgqbgV,qgqbgFccPMMM,qgqbgFscPMMM

      qgqbgA6PMMMmudep = qgqbgA60PMMM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >                 * qgqbgV(ONE,TWO,THREE,FOUR,FIVE,SIX)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qgqbgA6PMMMmuindep(ONE,TWO,THREE,
     >                                               FOUR,FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qgqbgA60PMMM,qgqbgV,qgqbgFccPMMM,qgqbgFscPMMM

      qgqbgA6PMMMmuindep = 
     > + CMPLX(0D0,1D0)*qgqbgFccPMMM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     > + CMPLX(0D0,1D0)*qgqbgFscPMMM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbggA6PMPPmudep(ONE,TWO,THREE,
     >                                             FOUR,FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qqbggA60PMPP,qqbggV,qqbggFccPMPP,qqbggFscPMPP

      qqbggA6PMPPmudep = qqbggA60PMPP(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >                 * qqbggV(ONE,TWO,THREE,FOUR,FIVE,SIX)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbggA6PMPPmuindep(ONE,TWO,THREE,
     >                                               FOUR,FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qqbggA60PMPP,qqbggV,qqbggFccPMPP,qqbggFscPMPP

      qqbggA6PMPPmuindep = 
     > + CMPLX(0D0,1D0)*qqbggFccPMPP(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     > + CMPLX(0D0,1D0)*qqbggFscPMPP(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbggA6PMPMmudep(ONE,TWO,THREE,
     >                                             FOUR,FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qqbggA60PMPM,qqbggV,qqbggFccPMPM,qqbggFscPMPM

      qqbggA6PMPMmudep = qqbggA60PMPM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >                 * qqbggV(ONE,TWO,THREE,FOUR,FIVE,SIX)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbggA6PMPMmuindep(ONE,TWO,THREE,
     >                                               FOUR,FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qqbggA60PMPM,qqbggV,qqbggFccPMPM,qqbggFscPMPM

      qqbggA6PMPMmuindep = 
     > + CMPLX(0D0,1D0)*qqbggFccPMPM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     > + CMPLX(0D0,1D0)*qqbggFscPMPM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbggA6PMMPmudep(ONE,TWO,THREE,
     >                                             FOUR,FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qqbggA60PMMP,qqbggV,qqbggFccPMMP,qqbggFscPMMP

      qqbggA6PMMPmudep = qqbggA60PMMP(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >                 * qqbggV(ONE,TWO,THREE,FOUR,FIVE,SIX)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbggA6PMMPmuindep(ONE,TWO,THREE,
     >                                               FOUR,FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qqbggA60PMMP,qqbggV,qqbggFccPMMP,qqbggFscPMMP

      qqbggA6PMMPmuindep = 
     > + CMPLX(0D0,1D0)*qqbggFccPMMP(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     > + CMPLX(0D0,1D0)*qqbggFscPMMP(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbggA6PMMMmudep(ONE,TWO,THREE,
     >                                             FOUR,FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qqbggA60PMMM,qqbggV,qqbggFccPMMM,qqbggFscPMMM

      qqbggA6PMMMmudep = qqbggA60PMMM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >                 * qqbggV(ONE,TWO,THREE,FOUR,FIVE,SIX)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbggA6PMMMmuindep(ONE,TWO,THREE,
     >                                               FOUR,FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qqbggA60PMMM,qqbggV,qqbggFccPMMM,qqbggFscPMMM

      qqbggA6PMMMmuindep = 
     > + CMPLX(0D0,1D0)*qqbggFccPMMM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     > + CMPLX(0D0,1D0)*qqbggFscPMMM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)

      END

************************************************************************
************************************************************************
************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qggqbA61PMPMmudep(ONE,TWO,THREE,
     >                                              FOUR,FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qggqbA6PMPMmudep,qqbggA6PMPMmudep,
     >                   qggqbA6SPMPM,qggqbA6TPMPM
      REAL(KIND(1D0)) XNf,XNu,XNd,XNc,XNa                      
      COMMON /QCD/ XNf,XNu,XNd,XNc,XNa

      qggqbA61PMPMmudep = 
     >             + qggqbA6PMPMmudep(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >  - 1D0/XNc**2*qqbggA6PMPMmudep(ONE,FOUR,THREE,TWO,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qggqbA61PMPMmuindep(ONE,TWO,THREE,
     >                                                FOUR,FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qggqbA6PMPMmuindep,qqbggA6PMPMmuindep,
     >                   qggqbA6SPMPM,qggqbA6TPMPM
      REAL(KIND(1D0)) XNf,XNu,XNd,XNc,XNa                      
      COMMON /QCD/ XNf,XNu,XNd,XNc,XNa

      qggqbA61PMPMmuindep = 
     >             + qggqbA6PMPMmuindep(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >  - 1D0/XNc**2*qqbggA6PMPMmuindep(ONE,FOUR,THREE,TWO,FIVE,SIX,A,B)
     >     - XNf/XNc*qggqbA6SPMPM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >     + 1D0/XNc*qggqbA6TPMPM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qggqbA61PPMMmudep(ONE,TWO,THREE,
     >                                              FOUR,FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qggqbA6PPMMmudep,qqbggA6PMMPmudep,
     >                   qggqbA6SPPMM,qggqbA6TPPMM
      REAL(KIND(1D0)) XNf,XNu,XNd,XNc,XNa                      
      COMMON /QCD/ XNf,XNu,XNd,XNc,XNa

      qggqbA61PPMMmudep = 
     >             + qggqbA6PPMMmudep(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >  - 1D0/XNc**2*qqbggA6PMMPmudep(ONE,FOUR,THREE,TWO,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qggqbA61PPMMmuindep(ONE,TWO,THREE,
     >                                                FOUR,FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qggqbA6PPMMmuindep,qqbggA6PMMPmuindep,
     >                   qggqbA6SPPMM,qggqbA6TPPMM
      REAL(KIND(1D0)) XNf,XNu,XNd,XNc,XNa                      
      COMMON /QCD/ XNf,XNu,XNd,XNc,XNa

      qggqbA61PPMMmuindep = 
     >             + qggqbA6PPMMmuindep(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >  - 1D0/XNc**2*qqbggA6PMMPmuindep(ONE,FOUR,THREE,TWO,FIVE,SIX,A,B)
     >     - XNf/XNc*qggqbA6SPPMM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >     + 1D0/XNc*qggqbA6TPPMM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qggqbA61PPPMmudep(ONE,TWO,THREE,
     >                                              FOUR,FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qggqbA6PPPMmudep,qqbggA6PMPPmudep,
     >                   qggqbA6SPPPM,qggqbA6TPPPM
      REAL(KIND(1D0)) XNf,XNu,XNd,XNc,XNa                      
      COMMON /QCD/ XNf,XNu,XNd,XNc,XNa

      qggqbA61PPPMmudep = 
     >             + qggqbA6PPPMmudep(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >  - 1D0/XNc**2*qqbggA6PMPPmudep(ONE,FOUR,THREE,TWO,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qggqbA61PPPMmuindep(ONE,TWO,THREE,
     >                                                FOUR,FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qggqbA6PPPMmuindep,qqbggA6PMPPmuindep,
     >                   qggqbA6SPPPM,qggqbA6TPPPM
      REAL(KIND(1D0)) XNf,XNu,XNd,XNc,XNa                      
      COMMON /QCD/ XNf,XNu,XNd,XNc,XNa

      qggqbA61PPPMmuindep = 
     >             + qggqbA6PPPMmuindep(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >  - 1D0/XNc**2*qqbggA6PMPPmuindep(ONE,FOUR,THREE,TWO,FIVE,SIX,A,B)
     >     - XNf/XNc*qggqbA6SPPPM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >     + 1D0/XNc*qggqbA6TPPPM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qggqbA61PMMMmudep(ONE,TWO,THREE,
     >                                              FOUR,FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qggqbA6PMMMmudep,qqbggA6PMMMmudep,
     >                   qggqbA6SPMMM,qggqbA6TPMMM
      REAL(KIND(1D0)) XNf,XNu,XNd,XNc,XNa                      
      COMMON /QCD/ XNf,XNu,XNd,XNc,XNa

      qggqbA61PMMMmudep = 
     >             + qggqbA6PMMMmudep(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >  - 1D0/XNc**2*qqbggA6PMMMmudep(ONE,FOUR,THREE,TWO,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qggqbA61PMMMmuindep(ONE,TWO,THREE,
     >                                                FOUR,FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qggqbA6PMMMmuindep,qqbggA6PMMMmuindep,
     >                   qggqbA6SPMMM,qggqbA6TPMMM
      REAL(KIND(1D0)) XNf,XNu,XNd,XNc,XNa                      
      COMMON /QCD/ XNf,XNu,XNd,XNc,XNa

      qggqbA61PMMMmuindep = 
     >             + qggqbA6PMMMmuindep(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >  - 1D0/XNc**2*qqbggA6PMMMmuindep(ONE,FOUR,THREE,TWO,FIVE,SIX,A,B)
     >     - XNf/XNc*qggqbA6SPMMM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >     + 1D0/XNc*qggqbA6TPMMM(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)

      END

************************************************************************
************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbggA63PMPMmudep(ONE,FOUR,TWO,
     >                                              THREE,FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,FOUR,TWO,THREE,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qggqbA6PPMMmudep,qggqbA6PMPMmudep,
     >                   qgqbgA6PPMMmudep,qgqbgA6PMMPmudep,
     >                   qqbggA6PMPMmudep,qqbggA6PMMPmudep

      qqbggA63PMPMmudep = 
     >             + qggqbA6PPMMmudep(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >             + qggqbA6PMPMmudep(ONE,THREE,TWO,FOUR,FIVE,SIX,A,B)
     >             + qgqbgA6PPMMmudep(ONE,TWO,FOUR,THREE,FIVE,SIX,A,B)
     >             + qgqbgA6PMMPmudep(ONE,THREE,FOUR,TWO,FIVE,SIX,A,B)
     >             + qqbggA6PMPMmudep(ONE,FOUR,TWO,THREE,FIVE,SIX,A,B)
     >             + qqbggA6PMMPmudep(ONE,FOUR,THREE,TWO,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbggA63PMPMmuindep(ONE,FOUR,TWO,
     >                                               THREE,FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,FOUR,TWO,THREE,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qggqbA6PPMMmuindep,qggqbA6PMPMmuindep,
     >                   qgqbgA6PPMMmuindep,qgqbgA6PMMPmuindep,
     >                   qqbggA6PMPMmuindep,qqbggA6PMMPmuindep

      qqbggA63PMPMmuindep = 
     >             + qggqbA6PPMMmuindep(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >             + qggqbA6PMPMmuindep(ONE,THREE,TWO,FOUR,FIVE,SIX,A,B)
     >             + qgqbgA6PPMMmuindep(ONE,TWO,FOUR,THREE,FIVE,SIX,A,B)
     >             + qgqbgA6PMMPmuindep(ONE,THREE,FOUR,TWO,FIVE,SIX,A,B)
     >             + qqbggA6PMPMmuindep(ONE,FOUR,TWO,THREE,FIVE,SIX,A,B)
     >             + qqbggA6PMMPmuindep(ONE,FOUR,THREE,TWO,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbggA63PMMPmudep(ONE,FOUR,TWO,
     >                                              THREE,FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,FOUR,TWO,THREE,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qggqbA6PPMMmudep,qggqbA6PMPMmudep,
     >                   qgqbgA6PPMMmudep,qgqbgA6PMMPmudep,
     >                   qqbggA6PMPMmudep,qqbggA6PMMPmudep

      qqbggA63PMMPmudep = 
     >             + qggqbA6PMPMmudep(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >             + qggqbA6PPMMmudep(ONE,THREE,TWO,FOUR,FIVE,SIX,A,B)
     >             + qgqbgA6PMMPmudep(ONE,TWO,FOUR,THREE,FIVE,SIX,A,B)
     >             + qgqbgA6PPMMmudep(ONE,THREE,FOUR,TWO,FIVE,SIX,A,B)
     >             + qqbggA6PMMPmudep(ONE,FOUR,TWO,THREE,FIVE,SIX,A,B)
     >             + qqbggA6PMPMmudep(ONE,FOUR,THREE,TWO,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbggA63PMMPmuindep(ONE,FOUR,TWO,
     >                                               THREE,FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,FOUR,TWO,THREE,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qggqbA6PPMMmuindep,qggqbA6PMPMmuindep,
     >                   qgqbgA6PPMMmuindep,qgqbgA6PMMPmuindep,
     >                   qqbggA6PMPMmuindep,qqbggA6PMMPmuindep

      qqbggA63PMMPmuindep = 
     >             + qggqbA6PMPMmuindep(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >             + qggqbA6PPMMmuindep(ONE,THREE,TWO,FOUR,FIVE,SIX,A,B)
     >             + qgqbgA6PMMPmuindep(ONE,TWO,FOUR,THREE,FIVE,SIX,A,B)
     >             + qgqbgA6PPMMmuindep(ONE,THREE,FOUR,TWO,FIVE,SIX,A,B)
     >             + qqbggA6PMMPmuindep(ONE,FOUR,TWO,THREE,FIVE,SIX,A,B)
     >             + qqbggA6PMPMmuindep(ONE,FOUR,THREE,TWO,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbggA63PMPPmudep(ONE,FOUR,TWO,
     >                                              THREE,FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,FOUR,TWO,THREE,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qggqbA6PPPMmudep,qgqbgA6PPMPmudep,
     >                   qqbggA6PMPPmudep

      qqbggA63PMPPmudep = 
     >             + qggqbA6PPPMmudep(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >             + qggqbA6PPPMmudep(ONE,THREE,TWO,FOUR,FIVE,SIX,A,B)
     >             + qgqbgA6PPMPmudep(ONE,TWO,FOUR,THREE,FIVE,SIX,A,B)
     >             + qgqbgA6PPMPmudep(ONE,THREE,FOUR,TWO,FIVE,SIX,A,B)
     >             + qqbggA6PMPPmudep(ONE,FOUR,TWO,THREE,FIVE,SIX,A,B)
     >             + qqbggA6PMPPmudep(ONE,FOUR,THREE,TWO,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbggA63PMPPmuindep(ONE,FOUR,TWO,
     >                                               THREE,FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,FOUR,TWO,THREE,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qggqbA6PPPMmuindep,qgqbgA6PPMPmuindep,
     >                   qqbggA6PMPPmuindep

      qqbggA63PMPPmuindep = 
     >             + qggqbA6PPPMmuindep(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >             + qggqbA6PPPMmuindep(ONE,THREE,TWO,FOUR,FIVE,SIX,A,B)
     >             + qgqbgA6PPMPmuindep(ONE,TWO,FOUR,THREE,FIVE,SIX,A,B)
     >             + qgqbgA6PPMPmuindep(ONE,THREE,FOUR,TWO,FIVE,SIX,A,B)
     >             + qqbggA6PMPPmuindep(ONE,FOUR,TWO,THREE,FIVE,SIX,A,B)
     >             + qqbggA6PMPPmuindep(ONE,FOUR,THREE,TWO,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbggA63PMMMmudep(ONE,FOUR,TWO,
     >                                              THREE,FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,FOUR,TWO,THREE,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qggqbA6PMMMmudep,qgqbgA6PMMMmudep,
     >                   qqbggA6PMMMmudep

      qqbggA63PMMMmudep = 
     >             + qggqbA6PMMMmudep(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >             + qggqbA6PMMMmudep(ONE,THREE,TWO,FOUR,FIVE,SIX,A,B)
     >             + qgqbgA6PMMMmudep(ONE,TWO,FOUR,THREE,FIVE,SIX,A,B)
     >             + qgqbgA6PMMMmudep(ONE,THREE,FOUR,TWO,FIVE,SIX,A,B)
     >             + qqbggA6PMMMmudep(ONE,FOUR,TWO,THREE,FIVE,SIX,A,B)
     >             + qqbggA6PMMMmudep(ONE,FOUR,THREE,TWO,FIVE,SIX,A,B)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qqbggA63PMMMmuindep(ONE,FOUR,TWO,
     >                                               THREE,FIVE,SIX,A,B)
      IMPLICIT NONE

      INTEGER ONE,FOUR,TWO,THREE,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qggqbA6PMMMmuindep,qgqbgA6PMMMmuindep,
     >                   qqbggA6PMMMmuindep

      qqbggA63PMMMmuindep = 
     >             + qggqbA6PMMMmuindep(ONE,TWO,THREE,FOUR,FIVE,SIX,A,B)
     >             + qggqbA6PMMMmuindep(ONE,THREE,TWO,FOUR,FIVE,SIX,A,B)
     >             + qgqbgA6PMMMmuindep(ONE,TWO,FOUR,THREE,FIVE,SIX,A,B)
     >             + qgqbgA6PMMMmuindep(ONE,THREE,FOUR,TWO,FIVE,SIX,A,B)
     >             + qqbggA6PMMMmuindep(ONE,FOUR,TWO,THREE,FIVE,SIX,A,B)
     >             + qqbggA6PMMMmuindep(ONE,FOUR,THREE,TWO,FIVE,SIX,A,B)

      END
