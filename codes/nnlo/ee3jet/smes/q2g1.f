************************************************************************

      REAL(KIND(1D0)) FUNCTION PSI2q1gBorn(ONE,TWO,THREE,FOUR,FIVE)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,I,IHELE
      REAL(KIND(1D0)) PSI(2),XNf,XNu,XNd,XNc,XNa,XNc2,CC
      COMPLEX(KIND(1D0)) qgA50,C(2,2,2)
      COMMON /COUPLINGS/ C
     >       /QCD/XNf,XNu,XNd,XNc,XNa

      XNc2 = XNc**2

C u type: I = 1
C d type: I = 2

      DO I = 1,2
      PSI(I) = 0D0
C IHELE = 1 -> veL
C IHELE = 2 -> veR
      DO IHELE = 1,2

      CC = C(IHELE,2,I)*CONJG(C(IHELE,2,I))
      CALL qgA50PPHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50)
      PSI(I) = PSI(I) + CC*qgA50*CONJG(qgA50)
      CALL qgA50PMHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50)
      PSI(I) = PSI(I) + CC*qgA50*CONJG(qgA50)

      CC = C(IHELE,1,I)*CONJG(C(IHELE,1,I))
      CALL qgA50MPHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50)
      PSI(I) = PSI(I) + CC*qgA50*CONJG(qgA50)
      CALL qgA50MMHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50)
      PSI(I) = PSI(I) + CC*qgA50*CONJG(qgA50)

      ENDDO

      ENDDO

      PSI2q1gBorn = XNu*PSI(1) + XNd*PSI(2)
c The averaging over initial state spin states are reinstated...AK
      PSI2q1gBorn = (XNc2-1D0)*PSI2q1gBorn

C---NORMALISATION OK, 4 PI ALPHA_EM = 4 PI ALPHA_S = 1

      END

************************************************************************

      REAL(KIND(1D0)) FUNCTION PSI2u1gBorn(ONE,TWO,THREE,FOUR,FIVE)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,I,IHELE
      REAL(KIND(1D0)) PSI,XNf,XNu,XNd,XNc,XNa,XNc2,CC
      COMPLEX(KIND(1D0)) qgA50,C(2,2,2)
      COMMON /COUPLINGS/ C
     >       /QCD/XNf,XNu,XNd,XNc,XNa

      XNc2 = XNc**2


      PSI = 0D0
C IHELE = 1 -> veL
C IHELE = 2 -> veR
      DO IHELE = 1,2

      CC = C(IHELE,2,1)*CONJG(C(IHELE,2,1))
      CALL qgA50PPHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50)
      PSI = PSI + CC*qgA50*CONJG(qgA50)
      CALL qgA50PMHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50)
      PSI = PSI + CC*qgA50*CONJG(qgA50)

      CC = C(IHELE,1,1)*CONJG(C(IHELE,1,1))
      CALL qgA50MPHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50)
      PSI = PSI + CC*qgA50*CONJG(qgA50)
      CALL qgA50MMHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50)
      PSI = PSI + CC*qgA50*CONJG(qgA50)

      ENDDO

      PSI2u1gBorn = PSI
c The averaging over initial state spin states are reinstated...AK
      PSI2u1gBorn = (XNc2-1D0)*PSI2u1gBorn

C---NORMALISATION OK, 4 PI ALPHA_EM = 4 PI ALPHA_S = 1

      END

************************************************************************

      REAL(KIND(1D0)) FUNCTION PSI2d1gBorn(ONE,TWO,THREE,FOUR,FIVE)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,I,IHELE
      REAL(KIND(1D0)) PSI,XNf,XNu,XNd,XNc,XNa,XNc2,CC
      COMPLEX(KIND(1D0)) qgA50,C(2,2,2)
      COMMON /COUPLINGS/ C
     >       /QCD/XNf,XNu,XNd,XNc,XNa

      XNc2 = XNc**2


      PSI = 0D0
C IHELE = 1 -> veL
C IHELE = 2 -> veR
      DO IHELE = 1,2

      CC = C(IHELE,2,2)*CONJG(C(IHELE,2,2))
      CALL qgA50PPHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50)
      PSI = PSI + CC*qgA50*CONJG(qgA50)
      CALL qgA50PMHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50)
      PSI = PSI + CC*qgA50*CONJG(qgA50)

      CC = C(IHELE,1,2)*CONJG(C(IHELE,1,2))
      CALL qgA50MPHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50)
      PSI = PSI + CC*qgA50*CONJG(qgA50)
      CALL qgA50MMHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50)
      PSI = PSI + CC*qgA50*CONJG(qgA50)

      ENDDO

      PSI2d1gBorn = PSI
c The averaging over initial state spin states are reinstated...AK
      PSI2d1gBorn = (XNc2-1D0)*PSI2d1gBorn

C---NORMALISATION OK, 4 PI ALPHA_EM = 4 PI ALPHA_S = 1

      END

************************************************************************

************************************************************************
************************************************************************

      REAL(KIND(1D0)) FUNCTION PSI2q1gBornV(ONE,TWO,THREE,FOUR,FIVE
     >                            ,VPP,VPM,VMP,VMM)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,I,IHELE
      REAL(KIND(1D0)) PSI(2),XNf,XNu,XNd,XNc,XNa,XNc2
      COMPLEX(KIND(1D0)) qgA50,C(2,2,2),VPP,VPM,VMP,VMM
     >                         ,MPP,MPM,MMP,MMM
     >          ,PSII(2)
      COMMON /COUPLINGS/ C
     >       /QCD/XNf,XNu,XNd,XNc,XNa

      DO I=1,2
         PSII(I) = 0D0
      ENDDO
      XNc2 = XNc**2

C u type: I = 1
C d type: I = 2

      DO I = 1,2
      PSI(I) = 0D0
C IHELE = 1 -> veL
C IHELE = 2 -> veR
      DO IHELE = 1,2

c         write(6,*)'c1,c2 = ',C(IHELE,1,I),C(IHELE,2,I)
c         write(6,*)

         CALL qgA50PPHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,MPP)
         MPP = C(IHELE,2,I)*MPP
         
         CALL qgA50PMHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,MPM)
         MPM = C(IHELE,2,I)*MPM
         
         CALL qgA50MPHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,MMP)
         MMP = C(IHELE,1,I)*MMP
         
         CALL qgA50MMHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,MMM)
         MMM = C(IHELE,1,I)*MMM

c         write(6,*)'MPP , VPP = ',MPP,VPP
c         write(6,*)'MPM , VPM = ',MPM,VPM
c         write(6,*)'MMP , VMP = ',MMP,VMP
c         write(6,*)'MMM , VMM = ',MMM,VMM
         
         PSII(I) = PSII(I)+CONJG(MPP)*VPP*MPP
     >                  +CONJG(MPP)*VPM*MPM
     >                  +CONJG(MPM)*VMP*MPP
     >                  +CONJG(MPM)*VMM*MPM
c
     >                  +CONJG(MMP)*VPP*MMP
     >                  +CONJG(MMP)*VPM*MMM
     >                  +CONJG(MMM)*VMP*MMP
     >                  +CONJG(MMM)*VMM*MMM

c         write(6,*)'PSII = ',CONJG(MPP)*VPP*MPP
c     >                  +CONJG(MPP)*VPM*MPM
c     >                  +CONJG(MPM)*VMP*MPP
c     >                  +CONJG(MPM)*VMM*MPM
cc
c     >                  +CONJG(MMP)*VPP*MMP
c     >                  +CONJG(MMP)*VPM*MMM
c     >                  +CONJG(MMM)*VMP*MMP
c     >                  +CONJG(MMM)*VMM*MMM

c         write(6,*)PSII

c         write(6,*)CONJG(MPP)*VPP*MPP
c         write(6,*)CONJG(MPP)*VPM*MPM
c         write(6,*)CONJG(MPM)*VMP*MPP
c         write(6,*)CONJG(MPM)*VMM*MPM
c         write(6,*)CONJG(MMP)*VPP*MMP
c         write(6,*)CONJG(MMP)*VPM*MMM
c         write(6,*)CONJG(MMM)*VMP*MMP
c         write(6,*)CONJG(MMM)*VMM*MMM

         PSI(I) = PSII(I)
c         write(6,*)'PSI',PSI
    
C      CC = C(IHELE,2,I)*CONJG(C(IHELE,2,I))
C      CALL qgA50PPHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50)
C      PSI(I) = PSI(I) + CC*qgA50*CONJG(qgA50)
C      CALL qgA50PMHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50)
C      PSI(I) = PSI(I) + CC*qgA50*CONJG(qgA50)
C
C      CC = C(IHELE,1,I)*CONJG(C(IHELE,1,I))
C      CALL qgA50MPHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50)
C      PSI(I) = PSI(I) + CC*qgA50*CONJG(qgA50)
C      CALL qgA50MMHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50)
C      PSI(I) = PSI(I) + CC*qgA50*CONJG(qgA50)

      ENDDO

      ENDDO
      
c      write(6,*)'PSI = ',PSI

      PSI2q1gBornV = XNu*PSI(1) + XNd*PSI(2)
      PSI2q1gBornV = 4D0*(XNc2-1D0)*PSI2q1gBornV

c      write(6,*)'PSI2q1gBornV = ',PSI2q1gBornV

      END

************************************************************************
************************************************************************

************************************************************************
************************************************************************

      REAL(KIND(1D0)) FUNCTION PSI2d1gBornV(ONE,TWO,THREE,FOUR,FIVE,
     >                                      M2IJ)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,I,IHELE
      REAL(KIND(1D0)) XNf,XNu,XNd,XNc,XNa,XNc2
      COMPLEX(KIND(1D0)) M2IJ(2,2)
      COMPLEX(KIND(1D0)) qgA50,C(2,2,2),MPP,MPM,MMP,MMM
      COMMON /COUPLINGS/ C
     >       /QCD/XNf,XNu,XNd,XNc,XNa

      XNc2 = XNc**2

      M2IJ = 0d0
C IHELE = 1 -> veL
C IHELE = 2 -> veR
c
c The ordering is 1_{q},2_{g},3_{q~}
c PPHe : 1^+,2^+
c PMHe : 1^+,2^-
c MPHe : 1^-,2^+
c MMHe : 1^-,2^-
c
      DO IHELE = 1,2

         CALL qgA50PPHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,MPP)
         MPP = C(IHELE,2,2)*MPP
         
         CALL qgA50PMHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,MPM)
         MPM = C(IHELE,2,2)*MPM
         
         CALL qgA50MPHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,MMP)
         MMP = C(IHELE,1,2)*MMP
         
         CALL qgA50MMHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,MMM)
         MMM = C(IHELE,1,2)*MMM

c (1^{\lambda_i},2^+) * (1^{\lambda_i},2^+)^*
         M2IJ(1,1) = M2IJ(1,1) + MPP*CONJG(MPP)
         M2IJ(1,1) = M2IJ(1,1) + MMP*CONJG(MMP)
c (1^{\lambda_i},2^-) * (1^{\lambda_i},2^-)^*
         M2IJ(2,2) = M2IJ(2,2) + MPM*CONJG(MPM)
         M2IJ(2,2) = M2IJ(2,2) + MMM*CONJG(MMM)
c (1^{\lambda_i},2^+) * (1^{\lambda_i},2^-)^*
         M2IJ(1,2) = M2IJ(1,2) - MPP*CONJG(MPM)
         M2IJ(1,2) = M2IJ(1,2) - MMP*CONJG(MMM)
c (1^{\lambda_i},2^-) * (1^{\lambda_i},2^+)^*
         M2IJ(2,1) = M2IJ(2,1) - MMM*CONJG(MMP)
         M2IJ(2,1) = M2IJ(2,1) - MPM*CONJG(MPP)

      ENDDO

c In the normalization the averaging over initial state spin states are
c included
      M2IJ = (XNc2-1D0)*M2IJ

      PSI2d1gbornV = M2IJ(1,1) + M2IJ(2,2)

      END

************************************************************************
************************************************************************

************************************************************************
************************************************************************

      REAL(KIND(1D0)) FUNCTION PSI2u1gBornV(ONE,TWO,THREE,FOUR,FIVE,
     >                                      M2IJ)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,I,IHELE
      REAL(KIND(1D0)) XNf,XNu,XNd,XNc,XNa,XNc2
      COMPLEX(KIND(1D0)) M2IJ(2,2)
      COMPLEX(KIND(1D0)) qgA50,C(2,2,2),MPP,MPM,MMP,MMM
      COMMON /COUPLINGS/ C
     >       /QCD/XNf,XNu,XNd,XNc,XNa

      XNc2 = XNc**2

      M2IJ = 0d0
C IHELE = 1 -> veL
C IHELE = 2 -> veR
c
c The ordering is 1_{q},2_{g},3_{q~}
c PPHe : 1^+,2^+
c PMHe : 1^+,2^-
c MPHe : 1^-,2^+
c MMHe : 1^-,2^-
c
      DO IHELE = 1,2

         CALL qgA50PPHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,MPP)
         MPP = C(IHELE,2,1)*MPP
         
         CALL qgA50PMHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,MPM)
         MPM = C(IHELE,2,1)*MPM
         
         CALL qgA50MPHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,MMP)
         MMP = C(IHELE,1,1)*MMP
         
         CALL qgA50MMHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,MMM)
         MMM = C(IHELE,1,1)*MMM

c (1^{\lambda_i},2^+) * (1^{\lambda_i},2^+)^*
         M2IJ(1,1) = M2IJ(1,1) + MPP*CONJG(MPP)
         M2IJ(1,1) = M2IJ(1,1) + MMP*CONJG(MMP)
c (1^{\lambda_i},2^-) * (1^{\lambda_i},2^-)^*
         M2IJ(2,2) = M2IJ(2,2) + MPM*CONJG(MPM)
         M2IJ(2,2) = M2IJ(2,2) + MMM*CONJG(MMM)
c (1^{\lambda_i},2^+) * (1^{\lambda_i},2^-)^*
         M2IJ(1,2) = M2IJ(1,2) - MPP*CONJG(MPM)
         M2IJ(1,2) = M2IJ(1,2) - MMP*CONJG(MMM)
c (1^{\lambda_i},2^-) * (1^{\lambda_i},2^+)^*
         M2IJ(2,1) = M2IJ(2,1) - MMM*CONJG(MMP)
         M2IJ(2,1) = M2IJ(2,1) - MPM*CONJG(MPP)

      ENDDO

c In the normalization the averaging over initial state spin states are
c included
      M2IJ = (XNc2-1D0)*M2IJ

      PSI2u1gbornV = M2IJ(1,1) + M2IJ(2,2)

      END

************************************************************************
************************************************************************

************************************************************************

      REAL(KIND(1D0)) FUNCTION PSI2q1gvirtNLO(ONE,TWO,THREE,FOUR,FIVE)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE
      REAL(KIND(1D0)) PSI2q1gNLO,PSI2q1gaxNLO

      PSI2q1gvirtNLO = PSI2q1gNLO(ONE,TWO,THREE,FOUR,FIVE)
     >               + PSI2q1gaxNLO(ONE,TWO,THREE,FOUR,FIVE)

      END

************************************************************************

      REAL(KIND(1D0)) FUNCTION PSI2u1gvirtNLO(ONE,TWO,THREE,FOUR,FIVE)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE
      REAL(KIND(1D0)) PSI2u1gNLO,PSI2u1gaxNLO

      PSI2u1gvirtNLO = PSI2u1gNLO(ONE,TWO,THREE,FOUR,FIVE)
     >               + PSI2u1gaxNLO(ONE,TWO,THREE,FOUR,FIVE)

      END

************************************************************************

      REAL(KIND(1D0)) FUNCTION PSI2u1gvirtNLOmuindep(ONE,TWO,THREE,
     >                                               FOUR,FIVE)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE
      REAL(KIND(1D0)) PSI2u1gNLOmuindep,PSI2u1gaxNLO

      PSI2u1gvirtNLOmuindep = PSI2u1gNLOmuindep(ONE,TWO,THREE,FOUR,FIVE)
     >                      + PSI2u1gaxNLO(ONE,TWO,THREE,FOUR,FIVE)

      END

************************************************************************

      REAL(KIND(1D0)) FUNCTION PSI2u1gvirtNLOmudep(ONE,TWO,THREE,FOUR,
     >                                             FIVE)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE
      REAL(KIND(1D0)) PSI2u1gNLOmudep

      PSI2u1gvirtNLOmudep = PSI2u1gNLOmudep(ONE,TWO,THREE,FOUR,FIVE)

      END

************************************************************************

      REAL(KIND(1D0)) FUNCTION PSI2u1gvirtNLOem2(ONE,TWO,THREE,FOUR,
     >                                           FIVE)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE
      REAL(KIND(1D0)) PSI2u1gNLOem2

      PSI2u1gvirtNLOem2 = PSI2u1gNLOem2(ONE,TWO,THREE,FOUR,FIVE)

      END

************************************************************************

      REAL(KIND(1D0)) FUNCTION PSI2u1gvirtNLOem1(ONE,TWO,THREE,FOUR,
     >                                           FIVE)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE
      REAL(KIND(1D0)) PSI2u1gNLOem1

      PSI2u1gvirtNLOem1 = PSI2u1gNLOem1(ONE,TWO,THREE,FOUR,FIVE)

      END

************************************************************************

      REAL(KIND(1D0)) FUNCTION PSI2d1gvirtNLO(ONE,TWO,THREE,FOUR,FIVE)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE
      REAL(KIND(1D0)) PSI2d1gNLO,PSI2d1gaxNLO

      PSI2d1gvirtNLO = PSI2d1gNLO(ONE,TWO,THREE,FOUR,FIVE)
     >               + PSI2d1gaxNLO(ONE,TWO,THREE,FOUR,FIVE)

      END

************************************************************************

      REAL(KIND(1D0)) FUNCTION PSI2d1gvirtNLOmuindep(ONE,TWO,THREE,FOUR,
     >                                               FIVE)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE
      REAL(KIND(1D0)) PSI2d1gNLOmuindep,PSI2d1gaxNLO

      PSI2d1gvirtNLOmuindep = PSI2d1gNLOmuindep(ONE,TWO,THREE,FOUR,FIVE)
     >                      + PSI2d1gaxNLO(ONE,TWO,THREE,FOUR,FIVE)

      END

************************************************************************

      REAL(KIND(1D0)) FUNCTION PSI2d1gvirtNLOmudep(ONE,TWO,THREE,FOUR,
     >                                           FIVE)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE
      REAL(KIND(1D0)) PSI2d1gNLOmudep

      PSI2d1gvirtNLOmudep = PSI2d1gNLOmudep(ONE,TWO,THREE,FOUR,FIVE)

      END

************************************************************************

      REAL(KIND(1D0)) FUNCTION PSI2d1gvirtNLOem2(ONE,TWO,THREE,FOUR,
     >                                           FIVE)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE
      REAL(KIND(1D0)) PSI2d1gNLOem2

      PSI2d1gvirtNLOem2 = PSI2d1gNLOem2(ONE,TWO,THREE,FOUR,FIVE)

      END

************************************************************************

      REAL(KIND(1D0)) FUNCTION PSI2d1gvirtNLOem1(ONE,TWO,THREE,FOUR,
     >                                           FIVE)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE
      REAL(KIND(1D0)) PSI2d1gNLOem1

      PSI2d1gvirtNLOem1 = PSI2d1gNLOem1(ONE,TWO,THREE,FOUR,FIVE)

      END


************************************************************************

      REAL(KIND(1D0)) FUNCTION PSI2q1gNLO(ONE,TWO,THREE,FOUR,FIVE)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,I,IHELE
      REAL(KIND(1D0)) XNf,XNu,XNd,XNc,XNa,XNc2,CC
      COMPLEX(KIND(1D0)) PSI(2),qgA50,qgA51,C(2,2,2)
      COMMON /COUPLINGS/ C
     >       /QCD/XNf,XNu,XNd,XNc,XNa

      XNc2 = XNc**2

C u type: I = 1
C d type: I = 2

      DO I = 1,2
      PSI(I) = 0D0
C IHELE = 1 -> veL
C IHELE = 2 -> veR
      DO IHELE = 1,2

      CC = C(IHELE,2,I)*CONJG(C(IHELE,2,I))
      CALL qgA51PPHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA51)
      CALL  qgA50PPHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50)
      PSI(I) = PSI(I) + CC*qgA51*CONJG(qgA50)
      CALL qgA51PMHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA51)
      CALL qgA50PMHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50)
      PSI(I) = PSI(I) + CC*qgA51*CONJG(qgA50)

      CC = C(IHELE,1,I)*CONJG(C(IHELE,1,I))
      CALL qgA51MPHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA51)
      CALL qgA50MPHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50)
      PSI(I) = PSI(I) + CC*qgA51*CONJG(qgA50)
      CALL qgA51MMHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA51)
      CALL qgA50MMHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50)
      PSI(I) = PSI(I) + CC*qgA51*CONJG(qgA50)

      ENDDO

      ENDDO

      PSI2q1gNLO =
     >   XNu*(PSI(1)+CONJG(PSI(1))) + XNd*(PSI(2)+CONJG(PSI(2)))
      PSI2q1gNLO = 4D0*(XNc2-1D0)*XNc*PSI2q1gNLO

      END

************************************************************************

      REAL(KIND(1D0)) FUNCTION PSI2u1gNLO(ONE,TWO,THREE,FOUR,FIVE)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,I,IHELE
      REAL(KIND(1D0)) XNf,XNu,XNd,XNc,XNa,XNc2,CC
      COMPLEX(KIND(1D0)) PSI,qgA50,qgA51,C(2,2,2)
      COMMON /COUPLINGS/ C
     >       /QCD/XNf,XNu,XNd,XNc,XNa

      XNc2 = XNc**2

      PSI = 0D0
C IHELE = 1 -> veL
C IHELE = 2 -> veR
      DO IHELE = 1,2

      CC = C(IHELE,2,1)*CONJG(C(IHELE,2,1))
      CALL qgA51PPHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA51)
      CALL  qgA50PPHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50)
      PSI = PSI + CC*qgA51*CONJG(qgA50)
      CALL qgA51PMHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA51)
      CALL qgA50PMHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50)
      PSI = PSI + CC*qgA51*CONJG(qgA50)

      CC = C(IHELE,1,1)*CONJG(C(IHELE,1,1))
      CALL qgA51MPHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA51)
      CALL qgA50MPHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50)
      PSI = PSI + CC*qgA51*CONJG(qgA50)
      CALL qgA51MMHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA51)
      CALL qgA50MMHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50)
      PSI = PSI + CC*qgA51*CONJG(qgA50)

      ENDDO

      PSI2u1gNLO = PSI+CONJG(PSI)

c In the normalization the averaging over initial state spin states are
c included
      PSI2u1gNLO = (XNc2-1D0)*XNc*PSI2u1gNLO

      END

************************************************************************

      REAL(KIND(1D0)) FUNCTION PSI2d1gNLO(ONE,TWO,THREE,FOUR,FIVE)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,I,IHELE
      REAL(KIND(1D0)) XNf,XNu,XNd,XNc,XNa,XNc2,CC
      COMPLEX(KIND(1D0)) PSI,qgA50,qgA51,C(2,2,2)
      COMMON /COUPLINGS/ C
     >       /QCD/XNf,XNu,XNd,XNc,XNa

      XNc2 = XNc**2

      PSI = 0D0
C IHELE = 1 -> veL
C IHELE = 2 -> veR
      DO IHELE = 1,2

      CC = C(IHELE,2,2)*CONJG(C(IHELE,2,2))
      CALL qgA51PPHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA51)
      CALL  qgA50PPHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50)
      PSI = PSI + CC*qgA51*CONJG(qgA50)
      CALL qgA51PMHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA51)
      CALL qgA50PMHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50)
      PSI = PSI + CC*qgA51*CONJG(qgA50)

      CC = C(IHELE,1,2)*CONJG(C(IHELE,1,2))
      CALL qgA51MPHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA51)
      CALL qgA50MPHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50)
      PSI = PSI + CC*qgA51*CONJG(qgA50)
      CALL qgA51MMHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA51)
      CALL qgA50MMHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50)
      PSI = PSI + CC*qgA51*CONJG(qgA50)

      ENDDO

      PSI2d1gNLO = PSI+CONJG(PSI)

c In the normalization the averaging over initial state spin states are
c included
      PSI2d1gNLO = (XNc2-1D0)*XNc*PSI2d1gNLO

      END

************************************************************************

      REAL(KIND(1D0)) FUNCTION PSI2u1gNLOmudep(ONE,TWO,THREE,FOUR,FIVE)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,I,IHELE
      REAL(KIND(1D0)) XNf,XNu,XNd,XNc,XNa,XNc2,CC
      COMPLEX(KIND(1D0)) PSI,qgA50,qgA51,C(2,2,2)
      COMMON /COUPLINGS/ C
     >       /QCD/XNf,XNu,XNd,XNc,XNa

      XNc2 = XNc**2

      PSI = 0D0
C IHELE = 1 -> veL
C IHELE = 2 -> veR
      DO IHELE = 1,2

      CC = C(IHELE,2,1)*CONJG(C(IHELE,2,1))
      CALL qgA51PPHemudep(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA51)
      CALL  qgA50PPHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50)
      PSI = PSI + CC*qgA51*CONJG(qgA50)
      CALL qgA51PMHemudep(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA51)
      CALL qgA50PMHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50)
      PSI = PSI + CC*qgA51*CONJG(qgA50)

      CC = C(IHELE,1,1)*CONJG(C(IHELE,1,1))
      CALL qgA51MPHemudep(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA51)
      CALL qgA50MPHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50)
      PSI = PSI + CC*qgA51*CONJG(qgA50)
      CALL qgA51MMHemudep(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA51)
      CALL qgA50MMHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50)
      PSI = PSI + CC*qgA51*CONJG(qgA50)

      ENDDO

      PSI2u1gNLOmudep = PSI+CONJG(PSI)

c In the normalization the averaging over initial state spin states are
c included
      PSI2u1gNLOmudep = (XNc2-1D0)*XNc*PSI2u1gNLOmudep

      END

************************************************************************

      REAL(KIND(1D0)) FUNCTION PSI2u1gNLOmuindep(ONE,TWO,THREE,FOUR,
     >                                           FIVE)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,I,IHELE
      REAL(KIND(1D0)) XNf,XNu,XNd,XNc,XNa,XNc2,CC
      COMPLEX(KIND(1D0)) PSI,qgA50,qgA51,C(2,2,2)
      COMMON /COUPLINGS/ C
     >       /QCD/XNf,XNu,XNd,XNc,XNa

      XNc2 = XNc**2

      PSI = 0D0
C IHELE = 1 -> veL
C IHELE = 2 -> veR
      DO IHELE = 1,2

      CC = C(IHELE,2,1)*CONJG(C(IHELE,2,1))
      CALL qgA51PPHemuindep(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA51)
      CALL  qgA50PPHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50)
      PSI = PSI + CC*qgA51*CONJG(qgA50)
      CALL qgA51PMHemuindep(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA51)
      CALL qgA50PMHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50)
      PSI = PSI + CC*qgA51*CONJG(qgA50)

      CC = C(IHELE,1,1)*CONJG(C(IHELE,1,1))
      CALL qgA51MPHemuindep(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA51)
      CALL qgA50MPHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50)
      PSI = PSI + CC*qgA51*CONJG(qgA50)
      CALL qgA51MMHemuindep(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA51)
      CALL qgA50MMHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50)
      PSI = PSI + CC*qgA51*CONJG(qgA50)

      ENDDO

      PSI2u1gNLOmuindep = PSI+CONJG(PSI)

c In the normalization the averaging over initial state spin states are
c included
      PSI2u1gNLOmuindep = (XNc2-1D0)*XNc*PSI2u1gNLOmuindep

      END

************************************************************************

      REAL(KIND(1D0)) FUNCTION PSI2u1gNLOem2(ONE,TWO,THREE,FOUR,FIVE)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,I,IHELE
      REAL(KIND(1D0)) XNf,XNu,XNd,XNc,XNa,XNc2,CC
      COMPLEX(KIND(1D0)) PSI,qgA50,qgA51,C(2,2,2)
      COMMON /COUPLINGS/ C
     >       /QCD/XNf,XNu,XNd,XNc,XNa

      XNc2 = XNc**2

      PSI = 0D0
C IHELE = 1 -> veL
C IHELE = 2 -> veR
      DO IHELE = 1,2

      CC = C(IHELE,2,1)*CONJG(C(IHELE,2,1))
      CALL qgA51PPHeem2(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA51)
      CALL  qgA50PPHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50)
      PSI = PSI + CC*qgA51*CONJG(qgA50)
      CALL qgA51PMHeem2(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA51)
      CALL qgA50PMHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50)
      PSI = PSI + CC*qgA51*CONJG(qgA50)

      CC = C(IHELE,1,1)*CONJG(C(IHELE,1,1))
      CALL qgA51MPHeem2(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA51)
      CALL qgA50MPHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50)
      PSI = PSI + CC*qgA51*CONJG(qgA50)
      CALL qgA51MMHeem2(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA51)
      CALL qgA50MMHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50)
      PSI = PSI + CC*qgA51*CONJG(qgA50)

      ENDDO

      PSI2u1gNLOem2 = PSI+CONJG(PSI)

c In the normalization the averaging over initial state spin states are
c included
      PSI2u1gNLOem2 = (XNc2-1D0)*XNc*PSI2u1gNLOem2

      END

************************************************************************

      REAL(KIND(1D0)) FUNCTION PSI2u1gNLOem1(ONE,TWO,THREE,FOUR,FIVE)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,I,IHELE
      REAL(KIND(1D0)) XNf,XNu,XNd,XNc,XNa,XNc2,CC
      COMPLEX(KIND(1D0)) PSI,qgA50,qgA51,C(2,2,2)
      COMMON /COUPLINGS/ C
     >       /QCD/XNf,XNu,XNd,XNc,XNa

      XNc2 = XNc**2

      PSI = 0D0
C IHELE = 1 -> veL
C IHELE = 2 -> veR
      DO IHELE = 1,2

      CC = C(IHELE,2,1)*CONJG(C(IHELE,2,1))
      CALL qgA51PPHeem1(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA51)
      CALL  qgA50PPHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50)
      PSI = PSI + CC*qgA51*CONJG(qgA50)
      CALL qgA51PMHeem1(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA51)
      CALL qgA50PMHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50)
      PSI = PSI + CC*qgA51*CONJG(qgA50)

      CC = C(IHELE,1,1)*CONJG(C(IHELE,1,1))
      CALL qgA51MPHeem1(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA51)
      CALL qgA50MPHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50)
      PSI = PSI + CC*qgA51*CONJG(qgA50)
      CALL qgA51MMHeem1(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA51)
      CALL qgA50MMHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50)
      PSI = PSI + CC*qgA51*CONJG(qgA50)

      ENDDO

      PSI2u1gNLOem1 = PSI+CONJG(PSI)

c In the normalization the averaging over initial state spin states are
c included
      PSI2u1gNLOem1 = (XNc2-1D0)*XNc*PSI2u1gNLOem1

      END

************************************************************************

      REAL(KIND(1D0)) FUNCTION PSI2d1gNLOmudep(ONE,TWO,THREE,FOUR,FIVE)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,I,IHELE
      REAL(KIND(1D0)) XNf,XNu,XNd,XNc,XNa,XNc2,CC
      COMPLEX(KIND(1D0)) PSI,qgA50,qgA51,C(2,2,2)
      COMMON /COUPLINGS/ C
     >       /QCD/XNf,XNu,XNd,XNc,XNa

      XNc2 = XNc**2

      PSI = 0D0
C IHELE = 1 -> veL
C IHELE = 2 -> veR
      DO IHELE = 1,2

      CC = C(IHELE,2,2)*CONJG(C(IHELE,2,2))
      CALL qgA51PPHemudep(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA51)
      CALL  qgA50PPHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50)
      PSI = PSI + CC*qgA51*CONJG(qgA50)
      CALL qgA51PMHemudep(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA51)
      CALL qgA50PMHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50)
      PSI = PSI + CC*qgA51*CONJG(qgA50)

      CC = C(IHELE,1,2)*CONJG(C(IHELE,1,2))
      CALL qgA51MPHemudep(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA51)
      CALL qgA50MPHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50)
      PSI = PSI + CC*qgA51*CONJG(qgA50)
      CALL qgA51MMHemudep(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA51)
      CALL qgA50MMHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50)
      PSI = PSI + CC*qgA51*CONJG(qgA50)

      ENDDO

      PSI2d1gNLOmudep = PSI+CONJG(PSI)

c In the normalization the averaging over initial state spin states are
c included
      PSI2d1gNLOmudep = (XNc2-1D0)*XNc*PSI2d1gNLOmudep

      END

************************************************************************

      REAL(KIND(1D0)) FUNCTION PSI2d1gNLOmuindep(ONE,TWO,THREE,FOUR,
     >                                           FIVE)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,I,IHELE
      REAL(KIND(1D0)) XNf,XNu,XNd,XNc,XNa,XNc2,CC
      COMPLEX(KIND(1D0)) PSI,qgA50,qgA51,C(2,2,2)
      COMMON /COUPLINGS/ C
     >       /QCD/XNf,XNu,XNd,XNc,XNa

      XNc2 = XNc**2

      PSI = 0D0
C IHELE = 1 -> veL
C IHELE = 2 -> veR
      DO IHELE = 1,2

      CC = C(IHELE,2,2)*CONJG(C(IHELE,2,2))
      CALL qgA51PPHemuindep(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA51)
      CALL  qgA50PPHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50)
      PSI = PSI + CC*qgA51*CONJG(qgA50)
      CALL qgA51PMHemuindep(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA51)
      CALL qgA50PMHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50)
      PSI = PSI + CC*qgA51*CONJG(qgA50)

      CC = C(IHELE,1,2)*CONJG(C(IHELE,1,2))
      CALL qgA51MPHemuindep(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA51)
      CALL qgA50MPHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50)
      PSI = PSI + CC*qgA51*CONJG(qgA50)
      CALL qgA51MMHemuindep(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA51)
      CALL qgA50MMHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50)
      PSI = PSI + CC*qgA51*CONJG(qgA50)

      ENDDO

      PSI2d1gNLOmuindep = PSI+CONJG(PSI)

c In the normalization the averaging over initial state spin states are
c included
      PSI2d1gNLOmuindep = (XNc2-1D0)*XNc*PSI2d1gNLOmuindep

      END

************************************************************************

      REAL(KIND(1D0)) FUNCTION PSI2d1gNLOem2(ONE,TWO,THREE,FOUR,FIVE)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,I,IHELE
      REAL(KIND(1D0)) XNf,XNu,XNd,XNc,XNa,XNc2,CC
      COMPLEX(KIND(1D0)) PSI,qgA50,qgA51,C(2,2,2)
      COMMON /COUPLINGS/ C
     >       /QCD/XNf,XNu,XNd,XNc,XNa

      XNc2 = XNc**2

      PSI = 0D0
C IHELE = 1 -> veL
C IHELE = 2 -> veR
      DO IHELE = 1,2

      CC = C(IHELE,2,2)*CONJG(C(IHELE,2,2))
      CALL qgA51PPHeem2(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA51)
      CALL  qgA50PPHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50)
      PSI = PSI + CC*qgA51*CONJG(qgA50)
      CALL qgA51PMHeem2(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA51)
      CALL qgA50PMHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50)
      PSI = PSI + CC*qgA51*CONJG(qgA50)

      CC = C(IHELE,1,2)*CONJG(C(IHELE,1,2))
      CALL qgA51MPHeem2(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA51)
      CALL qgA50MPHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50)
      PSI = PSI + CC*qgA51*CONJG(qgA50)
      CALL qgA51MMHeem2(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA51)
      CALL qgA50MMHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50)
      PSI = PSI + CC*qgA51*CONJG(qgA50)

      ENDDO

      PSI2d1gNLOem2 = PSI+CONJG(PSI)

c In the normalization the averaging over initial state spin states are
c included
      PSI2d1gNLOem2 = (XNc2-1D0)*XNc*PSI2d1gNLOem2

      END

************************************************************************

      REAL(KIND(1D0)) FUNCTION PSI2d1gNLOem1(ONE,TWO,THREE,FOUR,FIVE)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,I,IHELE
      REAL(KIND(1D0)) XNf,XNu,XNd,XNc,XNa,XNc2,CC
      COMPLEX(KIND(1D0)) PSI,qgA50,qgA51,C(2,2,2)
      COMMON /COUPLINGS/ C
     >       /QCD/XNf,XNu,XNd,XNc,XNa

      XNc2 = XNc**2

      PSI = 0D0
C IHELE = 1 -> veL
C IHELE = 2 -> veR
      DO IHELE = 1,2

      CC = C(IHELE,2,2)*CONJG(C(IHELE,2,2))
      CALL qgA51PPHeem1(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA51)
      CALL  qgA50PPHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50)
      PSI = PSI + CC*qgA51*CONJG(qgA50)
      CALL qgA51PMHeem1(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA51)
      CALL qgA50PMHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50)
      PSI = PSI + CC*qgA51*CONJG(qgA50)

      CC = C(IHELE,1,2)*CONJG(C(IHELE,1,2))
      CALL qgA51MPHeem1(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA51)
      CALL qgA50MPHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50)
      PSI = PSI + CC*qgA51*CONJG(qgA50)
      CALL qgA51MMHeem1(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA51)
      CALL qgA50MMHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50)
      PSI = PSI + CC*qgA51*CONJG(qgA50)

      ENDDO

      PSI2d1gNLOem1 = PSI+CONJG(PSI)

c In the normalization the averaging over initial state spin states are
c included
      PSI2d1gNLOem1 = (XNc2-1D0)*XNc*PSI2d1gNLOem1

      END

************************************************************************

      REAL(KIND(1D0)) FUNCTION PSI2u1gaxNLO(ONE,TWO,THREE,FOUR,FIVE)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,I,IHELE
      REAL(KIND(1D0)) XNf,XNu,XNd,XNc,XNa,XNc2
      COMPLEX(KIND(1D0)) PSI,qgA50,qgA54,C(2,2,2),Cax(2),TEMP1,TEMP2
      COMMON /COUPLINGS/ C
     >       /AXIALCOUPLINGS/ Cax
     >       /QCD/XNf,XNu,XNd,XNc,XNa

      XNc2 = XNc**2

      PSI = 0D0
C IHELE = 1 -> veL
C IHELE = 2 -> veR
      DO IHELE = 1,2
      CALL qgA5axPMPHe(ONE,THREE,TWO,FOUR,FIVE,IHELE,qgA54)
      TEMP1 = Cax(IHELE)*qgA54
      CALL qgA50PPHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50)
      TEMP2 = C(IHELE,2,1)*qgA50
      PSI = PSI + TEMP1*CONJG(TEMP2)

      CALL qgA5axPMMHe(ONE,THREE,TWO,FOUR,FIVE,IHELE,qgA54)
      TEMP1 = Cax(IHELE)*qgA54
      CALL qgA50PMHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50)
      TEMP2 = C(IHELE,2,1)*qgA50
      PSI = PSI + TEMP1*CONJG(TEMP2)

      CALL qgA5axMPPHe(ONE,THREE,TWO,FOUR,FIVE,IHELE,qgA54)
      TEMP1 = Cax(IHELE)*qgA54
      CALL qgA50MPHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50)
      TEMP2 = C(IHELE,1,1)*qgA50
      PSI = PSI + TEMP1*CONJG(TEMP2)

      CALL qgA5axMPMHe(ONE,THREE,TWO,FOUR,FIVE,IHELE,qgA54)
      TEMP1 = Cax(IHELE)*qgA54
      CALL qgA50MMHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50)
      TEMP2 = C(IHELE,1,1)*qgA50
      PSI = PSI + TEMP1*CONJG(TEMP2)

      ENDDO

      PSI2u1gaxNLO = PSI + CONJG(PSI)

c In the normalization the averaging over initial state spin states are
c included
      PSI2u1gaxNLO = (XNc2-1D0)*PSI2u1gaxNLO

      END

************************************************************************

      REAL(KIND(1D0)) FUNCTION PSI2d1gaxNLO(ONE,TWO,THREE,FOUR,FIVE)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,I,IHELE
      REAL(KIND(1D0)) XNf,XNu,XNd,XNc,XNa,XNc2
      COMPLEX(KIND(1D0)) PSI,qgA50,qgA54,C(2,2,2),Cax(2),TEMP1,TEMP2
      COMMON /COUPLINGS/ C
     >       /AXIALCOUPLINGS/ Cax
     >       /QCD/XNf,XNu,XNd,XNc,XNa

      XNc2 = XNc**2

      PSI = 0D0
C IHELE = 1 -> veL
C IHELE = 2 -> veR
      DO IHELE = 1,2
      CALL qgA5axPMPHe(ONE,THREE,TWO,FOUR,FIVE,IHELE,qgA54)
      TEMP1 = Cax(IHELE)*qgA54
      CALL qgA50PPHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50)
      TEMP2 = C(IHELE,2,2)*qgA50
      PSI = PSI + TEMP1*CONJG(TEMP2)

      CALL qgA5axPMMHe(ONE,THREE,TWO,FOUR,FIVE,IHELE,qgA54)
      TEMP1 = Cax(IHELE)*qgA54
      CALL qgA50PMHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50)
      TEMP2 = C(IHELE,2,2)*qgA50
      PSI = PSI + TEMP1*CONJG(TEMP2)

      CALL qgA5axMPPHe(ONE,THREE,TWO,FOUR,FIVE,IHELE,qgA54)
      TEMP1 = Cax(IHELE)*qgA54
      CALL qgA50MPHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50)
      TEMP2 = C(IHELE,1,2)*qgA50
      PSI = PSI + TEMP1*CONJG(TEMP2)

      CALL qgA5axMPMHe(ONE,THREE,TWO,FOUR,FIVE,IHELE,qgA54)
      TEMP1 = Cax(IHELE)*qgA54
      CALL qgA50MMHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50)
      TEMP2 = C(IHELE,1,2)*qgA50
      PSI = PSI + TEMP1*CONJG(TEMP2)

      ENDDO

      PSI2d1gaxNLO = PSI + CONJG(PSI)

c In the normalization the averaging over initial state spin states are
c included
      PSI2d1gaxNLO = (XNc2-1D0)*PSI2d1gaxNLO

      END


************************************************************************

      REAL(KIND(1D0)) FUNCTION PSI2q1gaxNLO(ONE,TWO,THREE,FOUR,FIVE)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,I,IHELE
      REAL(KIND(1D0)) XNf,XNu,XNd,XNc,XNa,XNc2
      COMPLEX(KIND(1D0)) PSI(2),qgA50,qgA54,C(2,2,2),Cax(2),TEMP1,TEMP2
      COMMON /COUPLINGS/ C
     >       /AXIALCOUPLINGS/ Cax
     >       /QCD/XNf,XNu,XNd,XNc,XNa

      XNc2 = XNc**2

C u type: I = 1
C d type: I = 2

      DO I = 1,2
      PSI(I) = 0D0
C IHELE = 1 -> veL
C IHELE = 2 -> veR
      DO IHELE = 1,2
      CALL qgA5axPMPHe(ONE,THREE,TWO,FOUR,FIVE,IHELE,qgA54)
      TEMP1 = Cax(IHELE)*qgA54
      CALL qgA50PPHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50)
      TEMP2 = C(IHELE,2,I)*qgA50
      PSI(I) = PSI(I) + TEMP1*CONJG(TEMP2)

      CALL qgA5axPMMHe(ONE,THREE,TWO,FOUR,FIVE,IHELE,qgA54)
      TEMP1 = Cax(IHELE)*qgA54
      CALL qgA50PMHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50)
      TEMP2 = C(IHELE,2,I)*qgA50
      PSI(I) = PSI(I) + TEMP1*CONJG(TEMP2)

      CALL qgA5axMPPHe(ONE,THREE,TWO,FOUR,FIVE,IHELE,qgA54)
      TEMP1 = Cax(IHELE)*qgA54
      CALL qgA50MPHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50)
      TEMP2 = C(IHELE,1,I)*qgA50
      PSI(I) = PSI(I) + TEMP1*CONJG(TEMP2)

      CALL qgA5axMPMHe(ONE,THREE,TWO,FOUR,FIVE,IHELE,qgA54)
      TEMP1 = Cax(IHELE)*qgA54
      CALL qgA50MMHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50)
      TEMP2 = C(IHELE,1,I)*qgA50
      PSI(I) = PSI(I) + TEMP1*CONJG(TEMP2)

      ENDDO

      ENDDO

      PSI2q1gaxNLO =
     >   XNu*(PSI(1)+CONJG(PSI(1))) + XNd*(PSI(2)+CONJG(PSI(2)))
      PSI2q1gaxNLO = 4D0*(XNc2-1D0)*PSI2q1gaxNLO

      END

************************************************************************

      SUBROUTINE qgA5tree(ONE,TWO,THREE,FOUR,FIVE,IHELE,OUT)
      IMPLICIT NONE

C PPHe -> 1+,2+,He
C MPHe -> 1-,2+,He
C PMHe -> 1+,2-,He
C MMHe -> 1-,2-,He
C IHELE = 1 for 4+
C IHELE = 2 for 4-

      INTEGER ONE,TWO,THREE,FOUR,FIVE,IHELE
      COMPLEX(KIND(1D0)) OUT,A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qgqbA50PPM,qgqbA50PMM
      COMMON /SPINORPRODUCTS/ A,B

      ENTRY qgA50PPHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,OUT)
      IF (IHELE.EQ.2) THEN
       OUT = qgqbA50PPM(ONE,TWO,THREE,FOUR,FIVE,A,B)
      ELSE
       OUT = qgqbA50PPM(ONE,TWO,THREE,FIVE,FOUR,A,B)
      ENDIF
      RETURN

      ENTRY qgA50MMHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,OUT)
      IF (IHELE.EQ.2) THEN
       OUT = -qgqbA50PPM(ONE,TWO,THREE,FIVE,FOUR,B,A)
      ELSE
       OUT = -qgqbA50PPM(ONE,TWO,THREE,FOUR,FIVE,B,A)
      ENDIF
      RETURN

      ENTRY qgA50PMHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,OUT)
      IF (IHELE.EQ.2) THEN
       OUT = qgqbA50PMM(ONE,TWO,THREE,FOUR,FIVE,A,B)
      ELSE
       OUT = qgqbA50PMM(ONE,TWO,THREE,FIVE,FOUR,A,B)
      ENDIF
      RETURN

      ENTRY qgA50MPHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,OUT)
      IF (IHELE.EQ.2) THEN
       OUT = -qgqbA50PMM(ONE,TWO,THREE,FIVE,FOUR,B,A)
      ELSE
       OUT = -qgqbA50PMM(ONE,TWO,THREE,FOUR,FIVE,B,A)
      ENDIF
      RETURN

      END

************************************************************************

      SUBROUTINE qgA5loop(ONE,TWO,THREE,FOUR,FIVE,IHELE,OUT)
      IMPLICIT NONE

C PPHe -> 1+,2+,He
C MPHe -> 1-,2+,He
C PMHe -> 1+,2-,He
C MMHe -> 1-,2-,He
C IHELE = 1 for 4+
C IHELE = 2 for 4-

      INTEGER ONE,TWO,THREE,FOUR,FIVE,IHELE
      COMPLEX(KIND(1D0)) OUT,A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qgqbA51PPM,qgqbA51PMM,qqbgA5axPMP,qqbgA5axPMM
      COMMON /SPINORPRODUCTS/ A,B

      ENTRY qgA51PPHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,OUT)
      IF (IHELE.EQ.2) THEN
       OUT = qgqbA51PPM(ONE,TWO,THREE,FOUR,FIVE,A,B)
      ELSE
       OUT = qgqbA51PPM(ONE,TWO,THREE,FIVE,FOUR,A,B)
      ENDIF
      RETURN

      ENTRY qgA51MMHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,OUT)
      IF (IHELE.EQ.2) THEN
       OUT = -qgqbA51PPM(ONE,TWO,THREE,FIVE,FOUR,B,A)
      ELSE
       OUT = -qgqbA51PPM(ONE,TWO,THREE,FOUR,FIVE,B,A)
      ENDIF
      RETURN

      ENTRY qgA51PMHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,OUT)
      IF (IHELE.EQ.2) THEN
       OUT = qgqbA51PMM(ONE,TWO,THREE,FOUR,FIVE,A,B)
      ELSE
       OUT = qgqbA51PMM(ONE,TWO,THREE,FIVE,FOUR,A,B)
      ENDIF
      RETURN

      ENTRY qgA51MPHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,OUT)
      IF (IHELE.EQ.2) THEN
       OUT = -qgqbA51PMM(ONE,TWO,THREE,FIVE,FOUR,B,A)
      ELSE
       OUT = -qgqbA51PMM(ONE,TWO,THREE,FOUR,FIVE,B,A)
      ENDIF
      RETURN

      ENTRY qgA5axPMPHe(ONE,THREE,TWO,FOUR,FIVE,IHELE,OUT)
      IF (IHELE.EQ.2) THEN
       OUT = qqbgA5axPMP(ONE,THREE,TWO,FOUR,FIVE,A,B)
      ELSE
       OUT = qqbgA5axPMP(ONE,THREE,TWO,FIVE,FOUR,A,B)
      ENDIF
      RETURN

      ENTRY qgA5axMPMHe(ONE,THREE,TWO,FOUR,FIVE,IHELE,OUT)
      IF (IHELE.EQ.2) THEN
       OUT = -qqbgA5axPMP(ONE,THREE,TWO,FIVE,FOUR,B,A)
      ELSE
       OUT = -qqbgA5axPMP(ONE,THREE,TWO,FOUR,FIVE,B,A)
      ENDIF
      RETURN

      ENTRY qgA5axPMMHe(ONE,THREE,TWO,FOUR,FIVE,IHELE,OUT)
      IF (IHELE.EQ.2) THEN
       OUT = qqbgA5axPMM(ONE,THREE,TWO,FOUR,FIVE,A,B)
      ELSE
       OUT = qqbgA5axPMM(ONE,THREE,TWO,FIVE,FOUR,A,B)
      ENDIF
      RETURN

      ENTRY qgA5axMPPHe(ONE,THREE,TWO,FOUR,FIVE,IHELE,OUT)
      IF (IHELE.EQ.2) THEN
       OUT = -qqbgA5axPMM(ONE,THREE,TWO,FIVE,FOUR,B,A)
      ELSE
       OUT = -qqbgA5axPMM(ONE,THREE,TWO,FOUR,FIVE,B,A)
      ENDIF
      RETURN

      END

************************************************************************

      SUBROUTINE qgA5loopem2(ONE,TWO,THREE,FOUR,FIVE,IHELE,OUT)
      IMPLICIT NONE

C PPHe -> 1+,2+,He
C MPHe -> 1-,2+,He
C PMHe -> 1+,2-,He
C MMHe -> 1-,2-,He
C IHELE = 1 for 4+
C IHELE = 2 for 4-

      INTEGER ONE,TWO,THREE,FOUR,FIVE,IHELE
      COMPLEX(KIND(1D0)) OUT,A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qgqbA51PPMem2,qgqbA51PMMem2
      COMMON /SPINORPRODUCTS/ A,B

      ENTRY qgA51PPHeem2(ONE,TWO,THREE,FOUR,FIVE,IHELE,OUT)
      IF (IHELE.EQ.2) THEN
       OUT = qgqbA51PPMem2(ONE,TWO,THREE,FOUR,FIVE,A,B)
      ELSE
       OUT = qgqbA51PPMem2(ONE,TWO,THREE,FIVE,FOUR,A,B)
      ENDIF
      RETURN

      ENTRY qgA51MMHeem2(ONE,TWO,THREE,FOUR,FIVE,IHELE,OUT)
      IF (IHELE.EQ.2) THEN
       OUT = -qgqbA51PPMem2(ONE,TWO,THREE,FIVE,FOUR,B,A)
      ELSE
       OUT = -qgqbA51PPMem2(ONE,TWO,THREE,FOUR,FIVE,B,A)
      ENDIF
      RETURN

      ENTRY qgA51PMHeem2(ONE,TWO,THREE,FOUR,FIVE,IHELE,OUT)
      IF (IHELE.EQ.2) THEN
       OUT = qgqbA51PMMem2(ONE,TWO,THREE,FOUR,FIVE,A,B)
      ELSE
       OUT = qgqbA51PMMem2(ONE,TWO,THREE,FIVE,FOUR,A,B)
      ENDIF
      RETURN

      ENTRY qgA51MPHeem2(ONE,TWO,THREE,FOUR,FIVE,IHELE,OUT)
      IF (IHELE.EQ.2) THEN
       OUT = -qgqbA51PMMem2(ONE,TWO,THREE,FIVE,FOUR,B,A)
      ELSE
       OUT = -qgqbA51PMMem2(ONE,TWO,THREE,FOUR,FIVE,B,A)
      ENDIF
      RETURN

      END

************************************************************************

      SUBROUTINE qgA5loopem1(ONE,TWO,THREE,FOUR,FIVE,IHELE,OUT)
      IMPLICIT NONE

C PPHe -> 1+,2+,He
C MPHe -> 1-,2+,He
C PMHe -> 1+,2-,He
C MMHe -> 1-,2-,He
C IHELE = 1 for 4+
C IHELE = 2 for 4-

      INTEGER ONE,TWO,THREE,FOUR,FIVE,IHELE
      COMPLEX(KIND(1D0)) OUT,A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qgqbA51PPMem1,qgqbA51PMMem1
      COMMON /SPINORPRODUCTS/ A,B

      ENTRY qgA51PPHeem1(ONE,TWO,THREE,FOUR,FIVE,IHELE,OUT)
      IF (IHELE.EQ.2) THEN
       OUT = qgqbA51PPMem1(ONE,TWO,THREE,FOUR,FIVE,A,B)
      ELSE
       OUT = qgqbA51PPMem1(ONE,TWO,THREE,FIVE,FOUR,A,B)
      ENDIF
      RETURN

      ENTRY qgA51MMHeem1(ONE,TWO,THREE,FOUR,FIVE,IHELE,OUT)
      IF (IHELE.EQ.2) THEN
       OUT = -qgqbA51PPMem1(ONE,TWO,THREE,FIVE,FOUR,B,A)
      ELSE
       OUT = -qgqbA51PPMem1(ONE,TWO,THREE,FOUR,FIVE,B,A)
      ENDIF
      RETURN

      ENTRY qgA51PMHeem1(ONE,TWO,THREE,FOUR,FIVE,IHELE,OUT)
      IF (IHELE.EQ.2) THEN
       OUT = qgqbA51PMMem1(ONE,TWO,THREE,FOUR,FIVE,A,B)
      ELSE
       OUT = qgqbA51PMMem1(ONE,TWO,THREE,FIVE,FOUR,A,B)
      ENDIF
      RETURN

      ENTRY qgA51MPHeem1(ONE,TWO,THREE,FOUR,FIVE,IHELE,OUT)
      IF (IHELE.EQ.2) THEN
       OUT = -qgqbA51PMMem1(ONE,TWO,THREE,FIVE,FOUR,B,A)
      ELSE
       OUT = -qgqbA51PMMem1(ONE,TWO,THREE,FOUR,FIVE,B,A)
      ENDIF
      RETURN

      END

************************************************************************

      SUBROUTINE qgA5loopmudep(ONE,TWO,THREE,FOUR,FIVE,IHELE,OUT)
      IMPLICIT NONE

C PPHe -> 1+,2+,He
C MPHe -> 1-,2+,He
C PMHe -> 1+,2-,He
C MMHe -> 1-,2-,He
C IHELE = 1 for 4+
C IHELE = 2 for 4-

      INTEGER ONE,TWO,THREE,FOUR,FIVE,IHELE
      COMPLEX(KIND(1D0)) OUT,A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qgqbA51PPMmudep,qgqbA51PMMmudep
      COMMON /SPINORPRODUCTS/ A,B

      ENTRY qgA51PPHemudep(ONE,TWO,THREE,FOUR,FIVE,IHELE,OUT)
      IF (IHELE.EQ.2) THEN
       OUT = qgqbA51PPMmudep(ONE,TWO,THREE,FOUR,FIVE,A,B)
      ELSE
       OUT = qgqbA51PPMmudep(ONE,TWO,THREE,FIVE,FOUR,A,B)
      ENDIF
      RETURN

      ENTRY qgA51MMHemudep(ONE,TWO,THREE,FOUR,FIVE,IHELE,OUT)
      IF (IHELE.EQ.2) THEN
       OUT = -qgqbA51PPMmudep(ONE,TWO,THREE,FIVE,FOUR,B,A)
      ELSE
       OUT = -qgqbA51PPMmudep(ONE,TWO,THREE,FOUR,FIVE,B,A)
      ENDIF
      RETURN

      ENTRY qgA51PMHemudep(ONE,TWO,THREE,FOUR,FIVE,IHELE,OUT)
      IF (IHELE.EQ.2) THEN
       OUT = qgqbA51PMMmudep(ONE,TWO,THREE,FOUR,FIVE,A,B)
      ELSE
       OUT = qgqbA51PMMmudep(ONE,TWO,THREE,FIVE,FOUR,A,B)
      ENDIF
      RETURN

      ENTRY qgA51MPHemudep(ONE,TWO,THREE,FOUR,FIVE,IHELE,OUT)
      IF (IHELE.EQ.2) THEN
       OUT = -qgqbA51PMMmudep(ONE,TWO,THREE,FIVE,FOUR,B,A)
      ELSE
       OUT = -qgqbA51PMMmudep(ONE,TWO,THREE,FOUR,FIVE,B,A)
      ENDIF
      RETURN

      END

************************************************************************

      SUBROUTINE qgA5loopmuindep(ONE,TWO,THREE,FOUR,FIVE,IHELE,OUT)
      IMPLICIT NONE

C PPHe -> 1+,2+,He
C MPHe -> 1-,2+,He
C PMHe -> 1+,2-,He
C MMHe -> 1-,2-,He
C IHELE = 1 for 4+
C IHELE = 2 for 4-

      INTEGER ONE,TWO,THREE,FOUR,FIVE,IHELE
      COMPLEX(KIND(1D0)) OUT,A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qgqbA51PPMmuindep,qgqbA51PMMmuindep
      COMMON /SPINORPRODUCTS/ A,B

      ENTRY qgA51PPHemuindep(ONE,TWO,THREE,FOUR,FIVE,IHELE,OUT)
      IF (IHELE.EQ.2) THEN
       OUT = qgqbA51PPMmuindep(ONE,TWO,THREE,FOUR,FIVE,A,B)
      ELSE
       OUT = qgqbA51PPMmuindep(ONE,TWO,THREE,FIVE,FOUR,A,B)
      ENDIF
      RETURN

      ENTRY qgA51MMHemuindep(ONE,TWO,THREE,FOUR,FIVE,IHELE,OUT)
      IF (IHELE.EQ.2) THEN
       OUT = -qgqbA51PPMmuindep(ONE,TWO,THREE,FIVE,FOUR,B,A)
      ELSE
       OUT = -qgqbA51PPMmuindep(ONE,TWO,THREE,FOUR,FIVE,B,A)
      ENDIF
      RETURN

      ENTRY qgA51PMHemuindep(ONE,TWO,THREE,FOUR,FIVE,IHELE,OUT)
      IF (IHELE.EQ.2) THEN
       OUT = qgqbA51PMMmuindep(ONE,TWO,THREE,FOUR,FIVE,A,B)
      ELSE
       OUT = qgqbA51PMMmuindep(ONE,TWO,THREE,FIVE,FOUR,A,B)
      ENDIF
      RETURN

      ENTRY qgA51MPHemuindep(ONE,TWO,THREE,FOUR,FIVE,IHELE,OUT)
      IF (IHELE.EQ.2) THEN
       OUT = -qgqbA51PMMmuindep(ONE,TWO,THREE,FIVE,FOUR,B,A)
      ELSE
       OUT = -qgqbA51PMMmuindep(ONE,TWO,THREE,FOUR,FIVE,B,A)
      ENDIF
      RETURN

      END

************************************************************************
**** The next section contains routines needed for the spin correlated
**** virtual part:
************************************************************************

c An extra selector is included: I = 1 -> u-type, I = 2 -> d-type 
      REAL(KIND(1D0)) FUNCTION PSI2q1gVmunu(ONE,TWO,THREE,FOUR,
     >                                      FIVE,I,M2IJ)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,I
      REAL(KIND(1D0)) PSI2q1gNLOmunu,PSI2q1gaxNLOmunu
      COMPLEX(KIND(1D0)) M2IJ(2,2),M2IJ_tmp(2,2)

      M2IJ = 0d0

      PSI2q1gVmunu = 
     >    PSI2q1gNLOmunu(ONE,TWO,THREE,FOUR,FIVE,I,M2IJ_tmp)
      M2IJ = M2IJ_tmp
      PSI2q1gVmunu = PSI2q1gVmunu 
     >  + PSI2q1gaxNLOmunu(ONE,TWO,THREE,FOUR,FIVE,I,M2IJ_tmp)
      M2IJ = M2IJ + M2IJ_tmp

      END

c An extra selector is included: I = 1 -> u-type, I = 2 -> d-type 
c This routine calculates only the muR dependent part of the virtual,
c to obtain the full, correct answer the muR independent one has to
c be added too:
      REAL(KIND(1D0)) FUNCTION PSI2q1gVmunu_mudep(ONE,TWO,THREE,FOUR,
     >                                            FIVE,I,M2IJ)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,I
      REAL(KIND(1D0)) PSI2q1gNLOmunu_mudep
      COMPLEX(KIND(1D0)) M2IJ(2,2),M2IJ_tmp(2,2)

      M2IJ = 0d0

      PSI2q1gVmunu_mudep = 
     >    PSI2q1gNLOmunu_mudep(ONE,TWO,THREE,FOUR,FIVE,I,M2IJ_tmp)
      M2IJ = M2IJ_tmp

      END

c An extra selector is included: I = 1 -> u-type, I = 2 -> d-type 
c This routine gives back the muR independent part of the virtual, it has
c to be combined with the previous one to get a meaningful result.
      REAL(KIND(1D0)) FUNCTION PSI2q1gVmunu_muindep(ONE,TWO,THREE,FOUR,
     >                                              FIVE,I,M2IJ)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,I
      REAL(KIND(1D0)) PSI2q1gNLOmunu_muindep,PSI2q1gaxNLOmunu
      COMPLEX(KIND(1D0)) M2IJ(2,2),M2IJ_tmp(2,2)

      M2IJ = 0d0

      PSI2q1gVmunu_muindep = 
     >    PSI2q1gNLOmunu_muindep(ONE,TWO,THREE,FOUR,FIVE,I,M2IJ_tmp)
      M2IJ = M2IJ_tmp
c The axial contribution remains the same provided there is no
c muR dependence at all in it:
      PSI2q1gVmunu_muindep = PSI2q1gVmunu_muindep 
     >  + PSI2q1gaxNLOmunu(ONE,TWO,THREE,FOUR,FIVE,I,M2IJ_tmp)
      M2IJ = M2IJ + M2IJ_tmp

      END

************************************************************************

************************************************************************

      REAL(KIND(1D0)) FUNCTION PSI2q1gNLOmunu(ONE,TWO,THREE,FOUR,
     >                                        FIVE,I,M2IJ)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,I,IHELE
      REAL(KIND(1D0)) XNf,XNu,XNd,XNc,XNa,XNc2,CC
      COMPLEX(KIND(1D0)) C(2,2,2),M2IJ(2,2),
     >                   M2PP,M2MM,M2PM,M2MP
      COMPLEX(KIND(1D0)) qgA51PP,qgA51MM,qgA51PM,qgA51MP,
     >                   qgA50PP,qgA50MM,qgA50PM,qgA50MP 
      COMMON /COUPLINGS/ C
     >       /QCD/XNf,XNu,XNd,XNc,XNa

      XNc2 = XNc**2

      M2IJ = 0D0

c M2XY: X and Y are corresponding to the gluon polarization:
      M2PP = 0D0 ; M2MM = 0D0 ; M2PM = 0D0 ; M2MP = 0D0
c
c The ordering is 1_{q},2_{g},3_{q~}
c PPHe : 1^+,2^+
c PMHe : 1^+,2^-
c MPHe : 1^-,2^+
c MMHe : 1^-,2^-

c To ease up the computation we calculate everything than
c use whatever is needed:
      CALL qgA51PPHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA51PP)
      CALL qgA51MMHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA51MM)
      CALL qgA51PMHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA51PM)
      CALL qgA51MPHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA51MP)
c ^v
      CALL qgA50PPHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50PP)
      CALL qgA50MMHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50MM)
      CALL qgA50PMHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50PM)
      CALL qgA50MPHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50MP)

C IHELE = 1 -> veL
C IHELE = 2 -> veR

      DO IHELE = 1,2
c <M0+|M1+> + <M1+|M0+>:
      CC = C(IHELE,2,I)*CONJG(C(IHELE,2,I))
      M2PP = M2PP + CC*(qgA51PP*CONJG(qgA50PP) + qgA50PP*CONJG(qgA51PP))
      CC = C(IHELE,1,I)*CONJG(C(IHELE,1,I))
      M2PP = M2PP + CC*(qgA51MP*CONJG(qgA50MP) + qgA50MP*CONJG(qgA51MP))
c <M0-|M1-> + <M1-|M0->:
      CC = C(IHELE,2,I)*CONJG(C(IHELE,2,I))
      M2MM = M2MM + CC*(qgA51PM*CONJG(qgA50PM) + qgA50PM*CONJG(qgA51PM))
      CC = C(IHELE,1,I)*CONJG(C(IHELE,1,I))
      M2MM = M2MM + CC*(qgA51MM*CONJG(qgA50MM) + qgA50MM*CONJG(qgA51MM))
c <M0+|M1-> + <M1+|M0->:
      CC = C(IHELE,2,I)*CONJG(C(IHELE,2,I))
      M2PM = M2PM + CC*(qgA51PP*CONJG(qgA50PM) + qgA50PP*CONJG(qgA51PM))
      CC = C(IHELE,1,I)*CONJG(C(IHELE,1,I))
      M2PM = M2PM + CC*(qgA51MP*CONJG(qgA50MM) + qgA50MP*CONJG(qgA51MM))
c <M0-|M1+> + <M1-|M0+>:
      CC = C(IHELE,2,I)*CONJG(C(IHELE,2,I))
      M2MP = M2MP + CC*(qgA51PM*CONJG(qgA50PP) + qgA50PM*CONJG(qgA51PP))
      CC = C(IHELE,1,I)*CONJG(C(IHELE,1,I))
      M2MP = M2MP + CC*(qgA51MM*CONJG(qgA50MP) + qgA50MM*CONJG(qgA51MP))
      ENDDO

      M2IJ(1,1) = M2PP
      M2IJ(2,2) = M2MM
      M2IJ(1,2) = -M2PM
      M2IJ(2,1) = -M2MP

      M2IJ = (XNc2-1D0)*XNc*M2IJ

      PSI2q1gNLOmunu = M2IJ(1,1) + M2IJ(2,2)

      END

************************************************************************

      REAL(KIND(1D0)) FUNCTION PSI2q1gNLOmunu_mudep(ONE,TWO,THREE,FOUR,
     >                                              FIVE,I,M2IJ)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,I,IHELE
      REAL(KIND(1D0)) XNf,XNu,XNd,XNc,XNa,XNc2,CC
      COMPLEX(KIND(1D0)) C(2,2,2),M2IJ(2,2),
     >                   M2PP,M2MM,M2PM,M2MP
      COMPLEX(KIND(1D0)) qgA51PP,qgA51MM,qgA51PM,qgA51MP,
     >                   qgA50PP,qgA50MM,qgA50PM,qgA50MP 
      COMMON /COUPLINGS/ C
     >       /QCD/XNf,XNu,XNd,XNc,XNa

      XNc2 = XNc**2

      M2IJ = 0D0

c M2XY: X and Y are corresponding to the gluon polarization:
      M2PP = 0D0 ; M2MM = 0D0 ; M2PM = 0D0 ; M2MP = 0D0
c
c The ordering is 1_{q},2_{g},3_{q~}
c PPHe : 1^+,2^+
c PMHe : 1^+,2^-
c MPHe : 1^-,2^+
c MMHe : 1^-,2^-

c To ease up the computation we calculate everything than
c use whatever is needed:
      CALL qgA51PPHemudep(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA51PP)
      CALL qgA51MMHemudep(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA51MM)
      CALL qgA51PMHemudep(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA51PM)
      CALL qgA51MPHemudep(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA51MP)
c ^v
      CALL qgA50PPHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50PP)
      CALL qgA50MMHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50MM)
      CALL qgA50PMHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50PM)
      CALL qgA50MPHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50MP)

C IHELE = 1 -> veL
C IHELE = 2 -> veR

      DO IHELE = 1,2
c <M0+|M1+> + <M1+|M0+>:
      CC = C(IHELE,2,I)*CONJG(C(IHELE,2,I))
      M2PP = M2PP + CC*(qgA51PP*CONJG(qgA50PP) + qgA50PP*CONJG(qgA51PP))
      CC = C(IHELE,1,I)*CONJG(C(IHELE,1,I))
      M2PP = M2PP + CC*(qgA51MP*CONJG(qgA50MP) + qgA50MP*CONJG(qgA51MP))
c <M0-|M1-> + <M1-|M0->:
      CC = C(IHELE,2,I)*CONJG(C(IHELE,2,I))
      M2MM = M2MM + CC*(qgA51PM*CONJG(qgA50PM) + qgA50PM*CONJG(qgA51PM))
      CC = C(IHELE,1,I)*CONJG(C(IHELE,1,I))
      M2MM = M2MM + CC*(qgA51MM*CONJG(qgA50MM) + qgA50MM*CONJG(qgA51MM))
c <M0+|M1-> + <M1+|M0->:
      CC = C(IHELE,2,I)*CONJG(C(IHELE,2,I))
      M2PM = M2PM + CC*(qgA51PP*CONJG(qgA50PM) + qgA50PP*CONJG(qgA51PM))
      CC = C(IHELE,1,I)*CONJG(C(IHELE,1,I))
      M2PM = M2PM + CC*(qgA51MP*CONJG(qgA50MM) + qgA50MP*CONJG(qgA51MM))
c <M0-|M1+> + <M1-|M0+>:
      CC = C(IHELE,2,I)*CONJG(C(IHELE,2,I))
      M2MP = M2MP + CC*(qgA51PM*CONJG(qgA50PP) + qgA50PM*CONJG(qgA51PP))
      CC = C(IHELE,1,I)*CONJG(C(IHELE,1,I))
      M2MP = M2MP + CC*(qgA51MM*CONJG(qgA50MP) + qgA50MM*CONJG(qgA51MP))
      ENDDO

      M2IJ(1,1) = M2PP
      M2IJ(2,2) = M2MM
      M2IJ(1,2) = -M2PM
      M2IJ(2,1) = -M2MP

      M2IJ = (XNc2-1D0)*XNc*M2IJ

      PSI2q1gNLOmunu_mudep = M2IJ(1,1) + M2IJ(2,2)

      END

************************************************************************

      REAL(KIND(1D0)) FUNCTION PSI2q1gNLOmunu_muindep(ONE,TWO,THREE,
     >                                                FOUR,FIVE,I,M2IJ)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,I,IHELE
      REAL(KIND(1D0)) XNf,XNu,XNd,XNc,XNa,XNc2,CC
      COMPLEX(KIND(1D0)) C(2,2,2),M2IJ(2,2),
     >                   M2PP,M2MM,M2PM,M2MP
      COMPLEX(KIND(1D0)) qgA51PP,qgA51MM,qgA51PM,qgA51MP,
     >                   qgA50PP,qgA50MM,qgA50PM,qgA50MP 
      COMMON /COUPLINGS/ C
     >       /QCD/XNf,XNu,XNd,XNc,XNa

      XNc2 = XNc**2

      M2IJ = 0D0

c M2XY: X and Y are corresponding to the gluon polarization:
      M2PP = 0D0 ; M2MM = 0D0 ; M2PM = 0D0 ; M2MP = 0D0
c
c The ordering is 1_{q},2_{g},3_{q~}
c PPHe : 1^+,2^+
c PMHe : 1^+,2^-
c MPHe : 1^-,2^+
c MMHe : 1^-,2^-

c To ease up the computation we calculate everything than
c use whatever is needed:
      CALL qgA51PPHemuindep(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA51PP)
      CALL qgA51MMHemuindep(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA51MM)
      CALL qgA51PMHemuindep(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA51PM)
      CALL qgA51MPHemuindep(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA51MP)
c ^v
      CALL qgA50PPHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50PP)
      CALL qgA50MMHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50MM)
      CALL qgA50PMHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50PM)
      CALL qgA50MPHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50MP)

C IHELE = 1 -> veL
C IHELE = 2 -> veR

      DO IHELE = 1,2
c <M0+|M1+> + <M1+|M0+>:
      CC = C(IHELE,2,I)*CONJG(C(IHELE,2,I))
      M2PP = M2PP + CC*(qgA51PP*CONJG(qgA50PP) + qgA50PP*CONJG(qgA51PP))
      CC = C(IHELE,1,I)*CONJG(C(IHELE,1,I))
      M2PP = M2PP + CC*(qgA51MP*CONJG(qgA50MP) + qgA50MP*CONJG(qgA51MP))
c <M0-|M1-> + <M1-|M0->:
      CC = C(IHELE,2,I)*CONJG(C(IHELE,2,I))
      M2MM = M2MM + CC*(qgA51PM*CONJG(qgA50PM) + qgA50PM*CONJG(qgA51PM))
      CC = C(IHELE,1,I)*CONJG(C(IHELE,1,I))
      M2MM = M2MM + CC*(qgA51MM*CONJG(qgA50MM) + qgA50MM*CONJG(qgA51MM))
c <M0+|M1-> + <M1+|M0->:
      CC = C(IHELE,2,I)*CONJG(C(IHELE,2,I))
      M2PM = M2PM + CC*(qgA51PP*CONJG(qgA50PM) + qgA50PP*CONJG(qgA51PM))
      CC = C(IHELE,1,I)*CONJG(C(IHELE,1,I))
      M2PM = M2PM + CC*(qgA51MP*CONJG(qgA50MM) + qgA50MP*CONJG(qgA51MM))
c <M0-|M1+> + <M1-|M0+>:
      CC = C(IHELE,2,I)*CONJG(C(IHELE,2,I))
      M2MP = M2MP + CC*(qgA51PM*CONJG(qgA50PP) + qgA50PM*CONJG(qgA51PP))
      CC = C(IHELE,1,I)*CONJG(C(IHELE,1,I))
      M2MP = M2MP + CC*(qgA51MM*CONJG(qgA50MP) + qgA50MM*CONJG(qgA51MP))
      ENDDO

      M2IJ(1,1) = M2PP
      M2IJ(2,2) = M2MM
      M2IJ(1,2) = -M2PM
      M2IJ(2,1) = -M2MP

      M2IJ = (XNc2-1D0)*XNc*M2IJ

      PSI2q1gNLOmunu_muindep = M2IJ(1,1) + M2IJ(2,2)

      END

************************************************************************

      REAL(KIND(1D0)) FUNCTION PSI2q1gaxNLOmunu(ONE,TWO,THREE,FOUR,
     >                                          FIVE,I,M2IJ)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,I,IHELE
      REAL(KIND(1D0)) XNf,XNu,XNd,XNc,XNa,XNc2
      COMPLEX(KIND(1D0)) qgA50,qgA54,C(2,2,2),Cax(2),CC,M2IJ(2,2),
     >                   M2PP,M2MM,M2PM,M2MP
      COMPLEX(KIND(1D0)) qgA54PMP,qgA54PMM,qgA54MPP,qgA54MPM,
     >                   qgA50PP,qgA50MM,qgA50PM,qgA50MP 
      COMMON /COUPLINGS/ C
     >       /AXIALCOUPLINGS/ Cax
     >       /QCD/XNf,XNu,XNd,XNc,XNa

      XNc2 = XNc**2

      M2IJ = 0D0

c M2XY: X and Y are corresponding to the gluon polarization:
      M2PP = 0D0 ; M2MM = 0D0 ; M2PM = 0D0 ; M2MP = 0D0
c
c Tree-level:
c The ordering is 1_{q},2_{g},3_{q~}
c PPHe : 1^+,2^+
c PMHe : 1^+,2^-
c MPHe : 1^-,2^+
c MMHe : 1^-,2^-
      CALL qgA50PPHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50PP)
      CALL qgA50MMHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50MM)
      CALL qgA50PMHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50PM)
      CALL qgA50MPHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50MP)
c
c Axial:
c The ordering is 1_{q},3_{q~},2_{g}
c PMPHe : 1^+,3^-,2^+
c PMMHe : 1^+,3^-,2^-
c MPPHe : 1^-,3^+,2^+
c MMPHe : 1^-,3^+,2^-
      CALL qgA5axPMPHe(ONE,THREE,TWO,FOUR,FIVE,IHELE,qgA54PMP)
      CALL qgA5axPMMHe(ONE,THREE,TWO,FOUR,FIVE,IHELE,qgA54PMM)
      CALL qgA5axMPPHe(ONE,THREE,TWO,FOUR,FIVE,IHELE,qgA54MPP)
      CALL qgA5axMPMHe(ONE,THREE,TWO,FOUR,FIVE,IHELE,qgA54MPM)

C IHELE = 1 -> veL
C IHELE = 2 -> veR
      DO IHELE = 1,2
c <M0+|M1+> + <M1+|M0+>:
      CC = Cax(IHELE)*CONJG(C(IHELE,2,I))
      M2PP = M2PP + CC*qgA54PMP*CONJG(qgA50PP) 
     >     + CONJG(CC)*qgA50PP*CONJG(qgA54PMP)
      CC = Cax(IHELE)*CONJG(C(IHELE,1,I))
      M2PP = M2PP + CC*qgA54MPP*CONJG(qgA50MP) 
     >     + CONJG(CC)*qgA50MP*CONJG(qgA54MPP)
c <M0-|M1-> + <M1-|M0->:
      CC = Cax(IHELE)*CONJG(C(IHELE,2,I))
      M2MM = M2MM + CC*qgA54PMM*CONJG(qgA50PM) 
     >     + CONJG(CC)*qgA50PM*CONJG(qgA54PMM)
      CC = Cax(IHELE)*CONJG(C(IHELE,1,I))
      M2MM = M2MM + CC*qgA54MPM*CONJG(qgA50MM) 
     >     + CONJG(CC)*qgA50MM*CONJG(qgA54MPM)
c <M0+|M1-> + <M1+|M0->:
      CC = Cax(IHELE)*CONJG(C(IHELE,2,I))
      M2PM = M2PM + CC*qgA54PMP*CONJG(qgA50PM) 
     >     + CONJG(CC)*qgA50PP*CONJG(qgA54PMM)
      CC = Cax(IHELE)*CONJG(C(IHELE,1,I))
      M2PM = M2PM + CC*qgA54MPP*CONJG(qgA50MM) 
     >     + CONJG(CC)*qgA50MP*CONJG(qgA54MPM)
c <M0-|M1+> + <M1-|M0+>:
      CC = Cax(IHELE)*CONJG(C(IHELE,2,I))
      M2MP = M2MP + CC*qgA54PMM*CONJG(qgA50PP) 
     >     + CONJG(CC)*qgA50PM*CONJG(qgA54PMP)
      CC = Cax(IHELE)*CONJG(C(IHELE,1,I))
      M2MP = M2MP + CC*qgA54MPM*CONJG(qgA50MP) 
     >     + CONJG(CC)*qgA50MM*CONJG(qgA54MPP)

      ENDDO

      M2IJ(1,1) = M2PP
      M2IJ(2,2) = M2MM
      M2IJ(1,2) = -M2PM
      M2IJ(2,1) = -M2MP

      M2IJ = (XNc2-1D0)*XNc*M2IJ

      PSI2q1gaxNLOmunu = M2IJ(1,1) + M2IJ(2,2)

      END

************************************************************************
************************************************************************
* The next section deals with the poles of the spin-correlated Virtual:
************************************************************************
************************************************************************

c An extra selector is included: I = 1 -> u-type, I = 2 -> d-type 
      REAL(KIND(1D0)) FUNCTION PSI2q1gVmunuem1(ONE,TWO,THREE,FOUR,
     >                                         FIVE,I,M2IJ)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,I
      REAL(KIND(1D0)) PSI2q1gNLOmunuem1
      COMPLEX(KIND(1D0)) M2IJ(2,2)

      M2IJ = 0d0

      PSI2q1gVmunuem1 = 
     >    PSI2q1gNLOmunuem1(ONE,TWO,THREE,FOUR,FIVE,I,M2IJ)

      END

************************************************************************

c An extra selector is included: I = 1 -> u-type, I = 2 -> d-type 
      REAL(KIND(1D0)) FUNCTION PSI2q1gVmunuem2(ONE,TWO,THREE,FOUR,
     >                                         FIVE,I,M2IJ)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,I
      REAL(KIND(1D0)) PSI2q1gNLOmunuem2
      COMPLEX(KIND(1D0)) M2IJ(2,2)

      M2IJ = 0d0

      PSI2q1gVmunuem2 = 
     >    PSI2q1gNLOmunuem2(ONE,TWO,THREE,FOUR,FIVE,I,M2IJ)

      END

************************************************************************

      REAL(KIND(1D0)) FUNCTION PSI2q1gNLOmunuem1(ONE,TWO,THREE,FOUR,
     >                                           FIVE,I,M2IJ)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,I,IHELE
      REAL(KIND(1D0)) XNf,XNu,XNd,XNc,XNa,XNc2,CC
      COMPLEX(KIND(1D0)) C(2,2,2),M2IJ(2,2),
     >                   M2PP,M2MM,M2PM,M2MP
      COMPLEX(KIND(1D0)) qgA51PP,qgA51MM,qgA51PM,qgA51MP,
     >                   qgA50PP,qgA50MM,qgA50PM,qgA50MP 
      COMMON /COUPLINGS/ C
     >       /QCD/XNf,XNu,XNd,XNc,XNa

      XNc2 = XNc**2

      M2IJ = 0D0

c M2XY: X and Y are corresponding to the gluon polarization:
      M2PP = 0D0 ; M2MM = 0D0 ; M2PM = 0D0 ; M2MP = 0D0
c
c The ordering is 1_{q},2_{g},3_{q~}
c PPHe : 1^+,2^+
c PMHe : 1^+,2^-
c MPHe : 1^-,2^+
c MMHe : 1^-,2^-

c To ease up the computation we calculate everything than
c use whatever is needed:
      CALL qgA51PPHeem1(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA51PP)
      CALL qgA51MMHeem1(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA51MM)
      CALL qgA51PMHeem1(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA51PM)
      CALL qgA51MPHeem1(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA51MP)
c ^v
      CALL qgA50PPHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50PP)
      CALL qgA50MMHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50MM)
      CALL qgA50PMHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50PM)
      CALL qgA50MPHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50MP)

C IHELE = 1 -> veL
C IHELE = 2 -> veR

      DO IHELE = 1,2
c <M0+|M1+> + <M1+|M0+>:
      CC = C(IHELE,2,I)*CONJG(C(IHELE,2,I))
      M2PP = M2PP + CC*(qgA51PP*CONJG(qgA50PP) + qgA50PP*CONJG(qgA51PP))
      CC = C(IHELE,1,I)*CONJG(C(IHELE,1,I))
      M2PP = M2PP + CC*(qgA51MP*CONJG(qgA50MP) + qgA50MP*CONJG(qgA51MP))
c <M0-|M1-> + <M1-|M0->:
      CC = C(IHELE,2,I)*CONJG(C(IHELE,2,I))
      M2MM = M2MM + CC*(qgA51PM*CONJG(qgA50PM) + qgA50PM*CONJG(qgA51PM))
      CC = C(IHELE,1,I)*CONJG(C(IHELE,1,I))
      M2MM = M2MM + CC*(qgA51MM*CONJG(qgA50MM) + qgA50MM*CONJG(qgA51MM))
c <M0+|M1-> + <M1+|M0->:
      CC = C(IHELE,2,I)*CONJG(C(IHELE,2,I))
      M2PM = M2PM + CC*(qgA51PP*CONJG(qgA50PM) + qgA50PP*CONJG(qgA51PM))
      CC = C(IHELE,1,I)*CONJG(C(IHELE,1,I))
      M2PM = M2PM + CC*(qgA51MP*CONJG(qgA50MM) + qgA50MP*CONJG(qgA51MM))
c <M0-|M1+> + <M1-|M0+>:
      CC = C(IHELE,2,I)*CONJG(C(IHELE,2,I))
      M2MP = M2MP + CC*(qgA51PM*CONJG(qgA50PP) + qgA50PM*CONJG(qgA51PP))
      CC = C(IHELE,1,I)*CONJG(C(IHELE,1,I))
      M2MP = M2MP + CC*(qgA51MM*CONJG(qgA50MP) + qgA50MM*CONJG(qgA51MP))
      ENDDO

      M2IJ(1,1) = M2PP
      M2IJ(2,2) = M2MM
      M2IJ(1,2) = -M2PM
      M2IJ(2,1) = -M2MP

      M2IJ = (XNc2-1D0)*XNc*M2IJ

      PSI2q1gNLOmunuem1 = M2IJ(1,1) + M2IJ(2,2)

      END

************************************************************************

      REAL(KIND(1D0)) FUNCTION PSI2q1gNLOmunuem2(ONE,TWO,THREE,FOUR,
     >                                           FIVE,I,M2IJ)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,I,IHELE
      REAL(KIND(1D0)) XNf,XNu,XNd,XNc,XNa,XNc2,CC
      COMPLEX(KIND(1D0)) C(2,2,2),M2IJ(2,2),
     >                   M2PP,M2MM,M2PM,M2MP
      COMPLEX(KIND(1D0)) qgA51PP,qgA51MM,qgA51PM,qgA51MP,
     >                   qgA50PP,qgA50MM,qgA50PM,qgA50MP 
      COMMON /COUPLINGS/ C
     >       /QCD/XNf,XNu,XNd,XNc,XNa

      XNc2 = XNc**2

      M2IJ = 0D0

c M2XY: X and Y are corresponding to the gluon polarization:
      M2PP = 0D0 ; M2MM = 0D0 ; M2PM = 0D0 ; M2MP = 0D0
c
c The ordering is 1_{q},2_{g},3_{q~}
c PPHe : 1^+,2^+
c PMHe : 1^+,2^-
c MPHe : 1^-,2^+
c MMHe : 1^-,2^-

c To ease up the computation we calculate everything than
c use whatever is needed:
      CALL qgA51PPHeem2(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA51PP)
      CALL qgA51MMHeem2(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA51MM)
      CALL qgA51PMHeem2(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA51PM)
      CALL qgA51MPHeem2(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA51MP)
c ^v
      CALL qgA50PPHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50PP)
      CALL qgA50MMHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50MM)
      CALL qgA50PMHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50PM)
      CALL qgA50MPHe(ONE,TWO,THREE,FOUR,FIVE,IHELE,qgA50MP)

C IHELE = 1 -> veL
C IHELE = 2 -> veR

      DO IHELE = 1,2
c <M0+|M1+> + <M1+|M0+>:
      CC = C(IHELE,2,I)*CONJG(C(IHELE,2,I))
      M2PP = M2PP + CC*(qgA51PP*CONJG(qgA50PP) + qgA50PP*CONJG(qgA51PP))
      CC = C(IHELE,1,I)*CONJG(C(IHELE,1,I))
      M2PP = M2PP + CC*(qgA51MP*CONJG(qgA50MP) + qgA50MP*CONJG(qgA51MP))
c <M0-|M1-> + <M1-|M0->:
      CC = C(IHELE,2,I)*CONJG(C(IHELE,2,I))
      M2MM = M2MM + CC*(qgA51PM*CONJG(qgA50PM) + qgA50PM*CONJG(qgA51PM))
      CC = C(IHELE,1,I)*CONJG(C(IHELE,1,I))
      M2MM = M2MM + CC*(qgA51MM*CONJG(qgA50MM) + qgA50MM*CONJG(qgA51MM))
c <M0+|M1-> + <M1+|M0->:
      CC = C(IHELE,2,I)*CONJG(C(IHELE,2,I))
      M2PM = M2PM + CC*(qgA51PP*CONJG(qgA50PM) + qgA50PP*CONJG(qgA51PM))
      CC = C(IHELE,1,I)*CONJG(C(IHELE,1,I))
      M2PM = M2PM + CC*(qgA51MP*CONJG(qgA50MM) + qgA50MP*CONJG(qgA51MM))
c <M0-|M1+> + <M1-|M0+>:
      CC = C(IHELE,2,I)*CONJG(C(IHELE,2,I))
      M2MP = M2MP + CC*(qgA51PM*CONJG(qgA50PP) + qgA50PM*CONJG(qgA51PP))
      CC = C(IHELE,1,I)*CONJG(C(IHELE,1,I))
      M2MP = M2MP + CC*(qgA51MM*CONJG(qgA50MP) + qgA50MM*CONJG(qgA51MP))
      ENDDO

      M2IJ(1,1) = M2PP
      M2IJ(2,2) = M2MM
      M2IJ(1,2) = -M2PM
      M2IJ(2,1) = -M2MP

      M2IJ = (XNc2-1D0)*XNc*M2IJ

      PSI2q1gNLOmunuem2 = M2IJ(1,1) + M2IJ(2,2)

      END
