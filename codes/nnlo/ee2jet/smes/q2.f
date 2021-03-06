************************************************************************

      REAL(KIND(1D0)) FUNCTION PSI2qBorn(ONE,TWO,FOUR,FIVE)
      IMPLICIT NONE

      INTEGER ONE,TWO,FOUR,FIVE,I,IHELE
      REAL(KIND(1D0)) PSI(2),XNf,XNu,XNd,XNc,XNa,XNc2,CC
      COMPLEX(KIND(1D0)) qgA40,C(2,2,2)
      COMMON /COUPLINGS/ C
     >       /QCD/XNf,XNu,XNd,XNc,XNa

C u type: I = 1
C d type: I = 2

      DO I = 1,2
      PSI(I) = 0D0
C IHELE = 1 -> veL
C IHELE = 2 -> veR
      DO IHELE = 1,2

C couplings(eL/eR,qL/qR,up/down)

      CC = C(IHELE,2,I)*CONJG(C(IHELE,2,I))
      CALL qgA40PM(ONE,TWO,FOUR,FIVE,IHELE,qgA40)
      PSI(I) = PSI(I) + CC*qgA40*CONJG(qgA40)

      CC = C(IHELE,1,I)*CONJG(C(IHELE,1,I))
      CALL qgA40MP(ONE,TWO,FOUR,FIVE,IHELE,qgA40)
      PSI(I) = PSI(I) + CC*qgA40*CONJG(qgA40)

      ENDDO

      ENDDO

      PSI2qBorn = XNu*PSI(1) + XNd*PSI(2)
c The averaging over initial state spin states are reinstated...AK
C      PSI2qBorn = (XNc2-1D0)*PSI2qBorn
      PSI2qBorn = XNc*PSI2qBorn

C---NORMALISATION OK, 4 PI ALPHA_EM = 4 PI ALPHA_S = 1

      END

************************************************************************

      REAL(KIND(1D0)) FUNCTION PSI2uBorn(ONE,TWO,FOUR,FIVE)
      IMPLICIT NONE

      INTEGER ONE,TWO,FOUR,FIVE,I,IHELE
      REAL(KIND(1D0)) PSI,XNf,XNu,XNd,XNc,XNa,XNc2,CC
      COMPLEX(KIND(1D0)) qgA40,C(2,2,2)
      COMMON /COUPLINGS/ C
     >       /QCD/XNf,XNu,XNd,XNc,XNa

      XNc2 = XNc**2


      PSI = 0D0
C IHELE = 1 -> veL
C IHELE = 2 -> veR
      DO IHELE = 1,2

      CC = C(IHELE,2,1)*CONJG(C(IHELE,2,1))
      CALL qgA40PM(ONE,TWO,FOUR,FIVE,IHELE,qgA40)
      PSI = PSI + CC*qgA40*CONJG(qgA40)

      CC = C(IHELE,1,1)*CONJG(C(IHELE,1,1))
      CALL qgA40MP(ONE,TWO,FOUR,FIVE,IHELE,qgA40)
      PSI = PSI + CC*qgA40*CONJG(qgA40)

      ENDDO
      PSI2uBorn = PSI
c The averaging over initial state spin states are reinstated...AK
c      PSI2uBorn = (XNc2-1D0)*PSI2uBorn
      PSI2uBorn = XNc*PSI2uBorn
C---NORMALISATION OK, 4 PI ALPHA_EM = 4 PI ALPHA_S = 1

      END

************************************************************************

      REAL(KIND(1D0)) FUNCTION PSI2dBorn(ONE,TWO,FOUR,FIVE)
      IMPLICIT NONE

      INTEGER ONE,TWO,FOUR,FIVE,I,IHELE
      REAL(KIND(1D0)) PSI,XNf,XNu,XNd,XNc,XNa,XNc2,CC
      COMPLEX(KIND(1D0)) qgA40,C(2,2,2)
      COMMON /COUPLINGS/ C
     >       /QCD/XNf,XNu,XNd,XNc,XNa

      XNc2 = XNc**2


      PSI = 0D0
C IHELE = 1 -> veL
C IHELE = 2 -> veR
      DO IHELE = 1,2

      CC = C(IHELE,2,2)*CONJG(C(IHELE,2,2))
      CALL qgA40PM(ONE,TWO,FOUR,FIVE,IHELE,qgA40)
      PSI = PSI + CC*qgA40*CONJG(qgA40)

      CC = C(IHELE,1,2)*CONJG(C(IHELE,1,2))
      CALL qgA40MP(ONE,TWO,FOUR,FIVE,IHELE,qgA40)
      PSI = PSI + CC*qgA40*CONJG(qgA40)

      ENDDO

      PSI2dBorn = PSI
c The averaging over initial state spin states are reinstated...AK
c      PSI2dBorn = (XNc2-1D0)*PSI2dBorn
      PSI2dBorn = XNc*PSI2dBorn
C---NORMALISATION OK, 4 PI ALPHA_EM = 4 PI ALPHA_S = 1

      END

************************************************************************

      SUBROUTINE qgA4tree(ONE,TWO,FOUR,FIVE,IHELE,OUT)
      IMPLICIT NONE

C MPHe -> 1-,2+,He
C PMHe -> 1+,2-,He
C IHELE = 1 for 4+
C IHELE = 2 for 4-

      INTEGER ONE,TWO,FOUR,FIVE,IHELE
      COMPLEX(KIND(1D0)) OUT,A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qqbA40PM
      COMMON /SPINORPRODUCTS/ A,B

      ENTRY qgA40PM(ONE,TWO,FOUR,FIVE,IHELE,OUT)
      IF (IHELE.EQ.2) THEN
       OUT = qqbA40PM(ONE,TWO,FIVE,FOUR,A,B)
      ELSE
       OUT = qqbA40PM(ONE,TWO,FOUR,FIVE,A,B)
      ENDIF
      RETURN

      ENTRY qgA40MP(ONE,TWO,FOUR,FIVE,IHELE,OUT)
      IF (IHELE.EQ.2) THEN
       OUT = qqbA40PM(ONE,TWO,FIVE,FOUR,B,A)
      ELSE
       OUT = qqbA40PM(ONE,TWO,FOUR,FIVE,B,A)
      ENDIF
      RETURN

      END

************************************************************************

      REAL(KIND(1D0)) FUNCTION PSI2qVirtNLO(ONE,TWO,FOUR,FIVE)
      IMPLICIT NONE

      INTEGER ONE,TWO,FOUR,FIVE
      REAL(KIND(1D0)) S(11,11)
      REAL(KIND(1D0)) PSI,XNf,XNu,XNd,XNc,XNa
      REAL(KIND(1D0)) MUSQ, COMU, EPSLOG
      COMMON /QCD/XNf,XNu,XNd,XNc,XNa
     >       /scales/ COMU
     >       /DOTPRODUCTS/ S
      REAL(KIND(1D0)) PI
      PARAMETER ( PI  = 3.14159265358979 D0)

      MUSQ = COMU**2
c LOG from epsilon factor
      EPSLOG = LOG(MUSQ/S(FOUR,FIVE))

c excluding r[e]^-1
      PSI = PI**2 -3D0*EPSLOG-EPSLOG**2

c including r[e]^-1
c      PSI = 7D0*PI**2 /6D0 -3D0*EPSLOG-EPSLOG**2

      PSI2qVirtNLO = PSI
      END
      
************************************************************************

      REAL(KIND(1D0)) FUNCTION PSI2qVirtNLOmuindep(ONE,TWO,FOUR,FIVE)
      IMPLICIT NONE

      INTEGER ONE,TWO,FOUR,FIVE
      REAL(KIND(1D0)) PSI
      REAL(KIND(1D0)) PI
      PARAMETER ( PI  = 3.14159265358979 D0)

c excluding r[e]^-1
      PSI = PI**2

c including r[e]^-1
c      PSI = 7D0/6D0 *PI**2

      PSI2qVirtNLOmuindep = PSI
      END
      
************************************************************************

      REAL(KIND(1D0)) FUNCTION PSI2qVirtNLOmudep(ONE,TWO,FOUR,FIVE)
      IMPLICIT NONE

      INTEGER ONE,TWO,FOUR,FIVE
      REAL(KIND(1D0)) S(11,11)
      REAL(KIND(1D0)) PSI,XNf,XNu,XNd,XNc,XNa
      REAL(KIND(1D0)) MUSQ, COMU, EPSLOG
      COMMON /QCD/XNf,XNu,XNd,XNc,XNa
     >       /scales/ COMU
     >       /DOTPRODUCTS/ S
      REAL(KIND(1D0)) PI
      PARAMETER ( PI  = 3.14159265358979 D0)

      MUSQ = COMU**2
c LOG from epsilon factor
      EPSLOG = LOG(MUSQ/S(FOUR,FIVE))

      PSI = -3D0*EPSLOG-EPSLOG**2

      PSI2qVirtNLOmudep = PSI
      END
      
************************************************************************

      REAL(KIND(1D0)) FUNCTION PSI2qVirtNLOem1(ONE,TWO,FOUR,FIVE)
      IMPLICIT NONE

      INTEGER ONE,TWO,FOUR,FIVE
      REAL(KIND(1D0)) S(11,11)
      REAL(KIND(1D0)) PSI,XNf,XNu,XNd,XNc,XNa
      REAL(KIND(1D0)) MUSQ, COMU, EPSLOG
      COMMON /QCD/XNf,XNu,XNd,XNc,XNa
     >       /scales/ COMU
     >       /DOTPRODUCTS/ S
      REAL(KIND(1D0)) PI
      PARAMETER ( PI  = 3.14159265358979 D0)

      MUSQ = COMU**2
c LOG from epsilon factor
      EPSLOG = LOG(MUSQ/S(FOUR,FIVE))

      PSI = -3D0 -2D0*EPSLOG

      PSI2qVirtNLOem1 = PSI
      END
      
************************************************************************

      REAL(KIND(1D0)) FUNCTION PSI2qVirtNLOem2(ONE,TWO,FOUR,FIVE)
      IMPLICIT NONE

      INTEGER ONE,TWO,FOUR,FIVE
      REAL(KIND(1D0)) PSI

      PSI = -2D0

      PSI2qVirtNLOem2 = PSI
      END
      
************************************************************************