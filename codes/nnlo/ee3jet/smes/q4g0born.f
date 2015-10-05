************************************************************************

      REAL(KIND(1D0)) FUNCTION PSI4qBorn12(ONE,TWO,THREE,FOUR,FIVE,SIX,
     >                                     IFL)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX,IFL
      REAL(KIND(1D0)) PSI4qI,PSI4qslI
     >      ,XNf,XNu,XNd,XNc,XNa
      COMMON /QCD/XNf,XNu,XNd,XNc,XNa

      PSI4qBorn12 = - PSI4qI(ONE,TWO,THREE,FOUR,FIVE,SIX,IFL)
     >            + (XNa-1D0)*PSI4qI(ONE,FOUR,THREE,TWO,FIVE,SIX,IFL)
     >            - 1D0/XNc*PSI4qslI(ONE,TWO,THREE,FOUR,FIVE,SIX,IFL)
      PSI4qBorn12 = -1D0/XNc*PSI4qBorn12
c A factor of 1/2 is needed since we are having T_R = 1/2:
      PSI4qBorn12 = PSI4qBorn12 / 2d0

      END

************************************************************************

      REAL(KIND(1D0)) FUNCTION PSI4qBorn13(ONE,TWO,THREE,FOUR,FIVE,SIX,
     >                                     IFL)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX,IFL
      REAL(KIND(1D0)) PSI4qI,PSI4qslI
     >      ,XNf,XNu,XNd,XNc,XNa
      COMMON /QCD/XNf,XNu,XNd,XNc,XNa

      PSI4qBorn13 = 2D0*PSI4qI(ONE,TWO,THREE,FOUR,FIVE,SIX,IFL)
     >        + 2D0*PSI4qI(ONE,FOUR,THREE,TWO,FIVE,SIX,IFL)
     >        + (XNc+1D0/XNc)*PSI4qslI(ONE,TWO,THREE,FOUR,FIVE,SIX,IFL)
      PSI4qBorn13 = -1D0/XNc*PSI4qBorn13
c A factor of 1/2 is needed since we are having T_R = 1/2:
      PSI4qBorn13 = PSI4qBorn13 / 2d0

      END

************************************************************************

      REAL(KIND(1D0)) FUNCTION PSI4qBorn14(ONE,TWO,THREE,FOUR,FIVE,SIX,
     >                                     IFL)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX,IFL
      REAL(KIND(1D0)) PSI4qI,PSI4qslI
     >      ,XNf,XNu,XNd,XNc,XNa
      COMMON /QCD/XNf,XNu,XNd,XNc,XNa

      PSI4qBorn14 = (XNa-1D0)*PSI4qI(ONE,TWO,THREE,FOUR,FIVE,SIX,IFL)
     >            - PSI4qI(ONE,FOUR,THREE,TWO,FIVE,SIX,IFL)
     >            - 1D0/XNc*PSI4qslI(ONE,TWO,THREE,FOUR,FIVE,SIX,IFL)
      PSI4qBorn14 = -1D0/XNc*PSI4qBorn14
c A factor of 1/2 is needed since we are having T_R = 1/2:
      PSI4qBorn14 = PSI4qBorn14 / 2d0

      END

************************************************************************

      REAL(KIND(1D0)) FUNCTION PSI2u2dBorn12(ONE,TWO,THREE,FOUR,FIVE,
     >                                       SIX)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      REAL(KIND(1D0)) PSI2u2d
     >      ,XNf,XNu,XNd,XNc,XNa
      COMMON /QCD/XNf,XNu,XNd,XNc,XNa

      PSI2u2dBorn12 = - PSI2u2d(ONE,TWO,THREE,FOUR,FIVE,SIX)
c     >              + (XNa-1D0)*PSI2u2d(ONE,TWO,THREE,FOUR,FIVE,SIX)
      PSI2u2dBorn12 = -1D0/XNc*PSI2u2dBorn12
c A factor of 1/2 is needed since we are having T_R = 1/2:
      PSI2u2dBorn12 = PSI2u2dBorn12 / 2d0

      END

************************************************************************

      REAL(KIND(1D0)) FUNCTION PSI2u2dBorn13(ONE,TWO,THREE,FOUR,FIVE,
     >                                       SIX)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      REAL(KIND(1D0)) PSI2u2d
     >      ,XNf,XNu,XNd,XNc,XNa
      COMMON /QCD/XNf,XNu,XNd,XNc,XNa

      PSI2u2dBorn13 = 2D0*PSI2u2d(ONE,TWO,THREE,FOUR,FIVE,SIX)
c     >              + 2D0*PSI2u2d(ONE,TWO,THREE,FOUR,FIVE,SIX)
      PSI2u2dBorn13 = -1D0/XNc*PSI2u2dBorn13
c A factor of 1/2 is needed since we are having T_R = 1/2:
      PSI2u2dBorn13 = PSI2u2dBorn13 / 2d0

      END

************************************************************************

      REAL(KIND(1D0)) FUNCTION PSI2u2dBorn14(ONE,TWO,THREE,FOUR,FIVE,
     >                                       SIX)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX,IFL
      REAL(KIND(1D0)) X,PSI2u2d
     >      ,XNf,XNu,XNd,XNc,XNa
      COMMON /QCD/XNf,XNu,XNd,XNc,XNa

      PSI2u2dBorn14 = (XNa-1D0)*PSI2u2d(ONE,TWO,THREE,FOUR,FIVE,SIX)
c     >              - PSI2u2d(ONE,TWO,THREE,FOUR,FIVE,SIX)
      PSI2u2dBorn14 = -1D0/XNc*PSI2u2dBorn14
c A factor of 1/2 is needed since we are having T_R = 1/2:
      PSI2u2dBorn14 = PSI2u2dBorn14 / 2d0

      END

************************************************************************

      REAL(KIND(1D0)) FUNCTION PSI2u2cBorn12(ONE,TWO,THREE,FOUR,FIVE,
     >                                       SIX)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      REAL(KIND(1D0)) PSI2u2c
     >      ,XNf,XNu,XNd,XNc,XNa
      COMMON /QCD/XNf,XNu,XNd,XNc,XNa

      PSI2u2cBorn12 = - PSI2u2c(ONE,TWO,THREE,FOUR,FIVE,SIX)
      PSI2u2cBorn12 = -1D0/XNc*PSI2u2cBorn12
c A factor of 1/2 is needed since we are having T_R = 1/2:
      PSI2u2cBorn12 = PSI2u2cBorn12 / 2d0

      END

************************************************************************

      REAL(KIND(1D0)) FUNCTION PSI2u2cBorn13(ONE,TWO,THREE,FOUR,FIVE,
     >                                       SIX)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      REAL(KIND(1D0)) PSI2u2c
     >      ,XNf,XNu,XNd,XNc,XNa
      COMMON /QCD/XNf,XNu,XNd,XNc,XNa

      PSI2u2cBorn13 = 2D0*PSI2u2c(ONE,TWO,THREE,FOUR,FIVE,SIX)
      PSI2u2cBorn13 = -1D0/XNc*PSI2u2cBorn13
c A factor of 1/2 is needed since we are having T_R = 1/2:
      PSI2u2cBorn13 = PSI2u2cBorn13 / 2d0

      END

************************************************************************

      REAL(KIND(1D0)) FUNCTION PSI2u2cBorn14(ONE,TWO,THREE,FOUR,FIVE,
     >                                       SIX)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX,IFL
      REAL(KIND(1D0)) X,PSI2u2c
     >      ,XNf,XNu,XNd,XNc,XNa
      COMMON /QCD/XNf,XNu,XNd,XNc,XNa

      PSI2u2cBorn14 = (XNa-1D0)*PSI2u2c(ONE,TWO,THREE,FOUR,FIVE,SIX)
      PSI2u2cBorn14 = -1D0/XNc*PSI2u2cBorn14
c A factor of 1/2 is needed since we are having T_R = 1/2:
      PSI2u2cBorn14 = PSI2u2cBorn14 / 2d0

      END

************************************************************************

      REAL(KIND(1D0)) FUNCTION PSI2d2sBorn12(ONE,TWO,THREE,FOUR,FIVE,
     >                                       SIX)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      REAL(KIND(1D0)) PSI2d2s
     >      ,XNf,XNu,XNd,XNc,XNa
      COMMON /QCD/XNf,XNu,XNd,XNc,XNa

      PSI2d2sBorn12 = - PSI2d2s(ONE,TWO,THREE,FOUR,FIVE,SIX)
      PSI2d2sBorn12 = -1D0/XNc*PSI2d2sBorn12
c A factor of 1/2 is needed since we are having T_R = 1/2:
      PSI2d2sBorn12 = PSI2d2sBorn12 / 2d0

      END

************************************************************************

      REAL(KIND(1D0)) FUNCTION PSI2d2sBorn13(ONE,TWO,THREE,FOUR,FIVE,
     >                                       SIX)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      REAL(KIND(1D0)) PSI2d2s
     >      ,XNf,XNu,XNd,XNc,XNa
      COMMON /QCD/XNf,XNu,XNd,XNc,XNa

      PSI2d2sBorn13 = 2D0*PSI2d2s(ONE,TWO,THREE,FOUR,FIVE,SIX)
      PSI2d2sBorn13 = -1D0/XNc*PSI2d2sBorn13
c A factor of 1/2 is needed since we are having T_R = 1/2:
      PSI2d2sBorn13 = PSI2d2sBorn13 / 2d0

      END

************************************************************************

      REAL(KIND(1D0)) FUNCTION PSI2d2sBorn14(ONE,TWO,THREE,FOUR,FIVE,
     >                                       SIX)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX,IFL
      REAL(KIND(1D0)) X,PSI2d2s
     >      ,XNf,XNu,XNd,XNc,XNa
      COMMON /QCD/XNf,XNu,XNd,XNc,XNa

      PSI2d2sBorn14 = (XNa-1D0)*PSI2d2s(ONE,TWO,THREE,FOUR,FIVE,SIX)
      PSI2d2sBorn14 = -1D0/XNc*PSI2d2sBorn14
c A factor of 1/2 is needed since we are having T_R = 1/2:
      PSI2d2sBorn14 = PSI2d2sBorn14 / 2d0

      END
