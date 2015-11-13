************************************************************************
C Four-quark two-lepton squared matrix elements at tree level.
C Symmetry factor is not included!
C PSI4qBorn(ONE,TWO,THREE,FOUR,FIVE,SIX,IFL,X): squared matrix element
C PSI4qBornIK(I,K,ONE,TWO,THREE,FOUR,FIVE,SIX,IFL,X): 
C color correlated squared matrix element I and K are the connected legs
C ONE:   q
C TWO:   qb
C THREE: Q
C FOUR:  QB
C FIVE:  lb
C SIX:   l
C IFL = 0: flavor summed
C IFL = 1: QB is u flavor
C IFL = 2: QB is d flavor
C x is the scale of the cm energy
************************************************************************

      REAL(KIND(1D0)) FUNCTION PSI4qBorn(ONE,TWO,THREE,FOUR,FIVE,SIX,IFL,X)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX,IFL
      REAL(KIND(1D0)) X,PSI4q,PSI4qsl
     >      ,XNf,XNu,XNd,XNc,XNa
      COMMON /QCD/XNf,XNu,XNd,XNc,XNa

      PSI4qBorn = PSI4q(ONE,TWO,THREE,FOUR,FIVE,SIX,IFL,X)
     >          + PSI4q(ONE,FOUR,THREE,TWO,FIVE,SIX,IFL,X)
     >          + 1D0/XNc*PSI4qsl(ONE,TWO,THREE,FOUR,FIVE,SIX,IFL,X)

      END

************************************************************************

      REAL(KIND(1D0)) FUNCTION PSI4qBornIK(I,K,ONE,TWO,THREE,FOUR,FIVE,SIX,IFL,X)
      IMPLICIT NONE

      INTEGER I,K,ONE,TWO,THREE,FOUR,FIVE,SIX,IFL
      REAL(KIND(1D0)) X,PSI4qBorn12,PSI4qBorn13,PSI4qBorn14
     >      ,XNf,XNu,XNd,XNc,XNa
      COMMON /QCD/XNf,XNu,XNd,XNc,XNa

      IF (I.EQ.1) THEN
       IF (K.EQ.2) THEN
        PSI4qBornIK = PSI4qBorn12(ONE,TWO,THREE,FOUR,FIVE,SIX,IFL,X)
       ELSEIF (K.EQ.3) THEN
        PSI4qBornIK = PSI4qBorn13(ONE,TWO,THREE,FOUR,FIVE,SIX,IFL,X)
       ELSEIF (K.EQ.4) THEN
        PSI4qBornIK = PSI4qBorn14(ONE,TWO,THREE,FOUR,FIVE,SIX,IFL,X)
       ENDIF
      ELSEIF (I.EQ.2) THEN
       IF (K.EQ.1) THEN
        PSI4qBornIK = PSI4qBorn12(ONE,TWO,THREE,FOUR,FIVE,SIX,IFL,X)
       ELSEIF (K.EQ.3) THEN
        PSI4qBornIK = PSI4qBorn14(ONE,TWO,THREE,FOUR,FIVE,SIX,IFL,X)
       ELSEIF (K.EQ.4) THEN
        PSI4qBornIK = PSI4qBorn13(ONE,TWO,THREE,FOUR,FIVE,SIX,IFL,X)
       ENDIF
      ELSEIF (I.EQ.3) THEN
       IF (K.EQ.1) THEN
        PSI4qBornIK = PSI4qBorn13(ONE,TWO,THREE,FOUR,FIVE,SIX,IFL,X)
       ELSEIF (K.EQ.2) THEN
        PSI4qBornIK = PSI4qBorn14(ONE,TWO,THREE,FOUR,FIVE,SIX,IFL,X)
       ELSEIF (K.EQ.4) THEN
        PSI4qBornIK = PSI4qBorn12(ONE,TWO,THREE,FOUR,FIVE,SIX,IFL,X)
       ENDIF
      ELSEIF (I.EQ.4) THEN
       IF (K.EQ.1) THEN
        PSI4qBornIK = PSI4qBorn14(ONE,TWO,THREE,FOUR,FIVE,SIX,IFL,X)
       ELSEIF (K.EQ.2) THEN
        PSI4qBornIK = PSI4qBorn13(ONE,TWO,THREE,FOUR,FIVE,SIX,IFL,X)
       ELSEIF (K.EQ.3) THEN
        PSI4qBornIK = PSI4qBorn12(ONE,TWO,THREE,FOUR,FIVE,SIX,IFL,X)
       ENDIF
      ELSE
       WRITE(6,*)'Indeces I and K are equal in PSI4qBornIK'
       STOP
      ENDIF

      END

************************************************************************

      REAL(KIND(1D0)) FUNCTION PSI4qBorn12(ONE,TWO,THREE,FOUR,FIVE,SIX,IFL,X)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX,IFL
      REAL(KIND(1D0)) X,PSI4q,PSI4qsl
     >      ,XNf,XNu,XNd,XNc,XNa
      COMMON /QCD/XNf,XNu,XNd,XNc,XNa

      PSI4qBorn12 = - PSI4q(ONE,TWO,THREE,FOUR,FIVE,SIX,IFL,X)
     >            + (XNa-1D0)*PSI4q(ONE,FOUR,THREE,TWO,FIVE,SIX,IFL,X)
     >            - 1D0/XNc*PSI4qsl(ONE,TWO,THREE,FOUR,FIVE,SIX,IFL,X)
      PSI4qBorn12 = -1D0/XNc*PSI4qBorn12

      END

************************************************************************

      REAL(KIND(1D0)) FUNCTION PSI4qBorn13(ONE,TWO,THREE,FOUR,FIVE,SIX,IFL,X)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX,IFL
      REAL(KIND(1D0)) X,PSI4q,PSI4qsl
     >      ,XNf,XNu,XNd,XNc,XNa
      COMMON /QCD/XNf,XNu,XNd,XNc,XNa

      PSI4qBorn13 = 2D0*PSI4q(ONE,TWO,THREE,FOUR,FIVE,SIX,IFL,X)
     >        + 2D0*PSI4q(ONE,FOUR,THREE,TWO,FIVE,SIX,IFL,X)
     >        + (XNc+1D0/XNc)*PSI4qsl(ONE,TWO,THREE,FOUR,FIVE,SIX,IFL,X)
      PSI4qBorn13 = -1D0/XNc*PSI4qBorn13

      END

************************************************************************

      REAL(KIND(1D0)) FUNCTION PSI4qBorn14(ONE,TWO,THREE,FOUR,FIVE,SIX,IFL,X)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX,IFL
      REAL(KIND(1D0)) X,PSI4q,PSI4qsl
     >      ,XNf,XNu,XNd,XNc,XNa
      COMMON /QCD/XNf,XNu,XNd,XNc,XNa

      PSI4qBorn14 = (XNa-1D0)*PSI4q(ONE,TWO,THREE,FOUR,FIVE,SIX,IFL,X)
     >            - PSI4q(ONE,FOUR,THREE,TWO,FIVE,SIX,IFL,X)
     >            - 1D0/XNc*PSI4qsl(ONE,TWO,THREE,FOUR,FIVE,SIX,IFL,X)
      PSI4qBorn14 = -1D0/XNc*PSI4qBorn14

      END

************************************************************************

      REAL(KIND(1D0)) FUNCTION PSI4q(ONE,TWO,THREE,FOUR,FIVE,SIX,IFL,X)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX,IFL,IHEL5
      COMPLEX(KIND(1D0)) TEMP1,TEMP2
      REAL(KIND(1D0)) X,PSI(3),CG(2),XNf,XNu,XNd,XNc,XNa
      COMMON /QCD/XNf,XNu,XNd,XNc,XNa

      CG(1) = 0.6666666666667D0
      CG(2) =-0.3333333333333D0

      PSI(1) = 0D0
      PSI(2) = 0D0
      PSI(3) = 0D0
      DO IHEL5 = 1,2

      CALL qQA60PMPMH5(ONE,TWO,THREE,FOUR,FIVE,SIX,IHEL5,X,TEMP1)
      CALL qQA60PMPMH5(THREE,FOUR,ONE,TWO,FIVE,SIX,IHEL5,X,TEMP2)
      PSI(1) = PSI(1) + TEMP1*CONJG(TEMP1)
      PSI(2) = PSI(2) + TEMP2*CONJG(TEMP2)
      PSI(3) = PSI(3) + TEMP1*CONJG(TEMP2) + TEMP2*CONJG(TEMP1)

      CALL qQA60PMMPH5(ONE,TWO,THREE,FOUR,FIVE,SIX,IHEL5,X,TEMP1)
      CALL qQA60MPPMH5(THREE,FOUR,ONE,TWO,FIVE,SIX,IHEL5,X,TEMP2)
      PSI(1) = PSI(1) + TEMP1*CONJG(TEMP1) 
      PSI(2) = PSI(2) + TEMP2*CONJG(TEMP2)
      PSI(3) = PSI(3) + TEMP1*CONJG(TEMP2) + TEMP2*CONJG(TEMP1)

      CALL qQA60MPMPH5(ONE,TWO,THREE,FOUR,FIVE,SIX,IHEL5,X,TEMP1)
      CALL qQA60MPMPH5(THREE,FOUR,ONE,TWO,FIVE,SIX,IHEL5,X,TEMP2)
      PSI(1) = PSI(1) + TEMP1*CONJG(TEMP1) 
      PSI(2) = PSI(2) + TEMP2*CONJG(TEMP2)
      PSI(3) = PSI(3) + TEMP1*CONJG(TEMP2) + TEMP2*CONJG(TEMP1)

      CALL qQA60MPPMH5(ONE,TWO,THREE,FOUR,FIVE,SIX,IHEL5,X,TEMP1)
      CALL qQA60PMMPH5(THREE,FOUR,ONE,TWO,FIVE,SIX,IHEL5,X,TEMP2)
      PSI(1) = PSI(1) + TEMP1*CONJG(TEMP1) 
      PSI(2) = PSI(2) + TEMP2*CONJG(TEMP2)
      PSI(3) = PSI(3) + TEMP1*CONJG(TEMP2) + TEMP2*CONJG(TEMP1)

      ENDDO

      IF (IFL.EQ.0) THEN
       PSI4q = (CG(1)**2*XNu+CG(2)**2*XNd)*XNf*(PSI(1)+PSI(2))
     >       + (CG(1)*XNu+CG(2)*XNd)**2*PSI(3)
      ELSEIF(IFL.EQ.1) THEN
       PSI4q = (CG(1)**2*XNu+CG(2)**2*XNd)*PSI(1)
     >       + CG(1)**2*XNf*PSI(2)
     >       + (CG(1)*XNu+CG(2)*XNd)*CG(1)*PSI(3)
      ELSEIF(IFL.EQ.2) THEN
       PSI4q = (CG(1)**2*XNu+CG(2)**2*XNd)*PSI(1)
     >       + CG(2)**2*XNf*PSI(2)
     >       + (CG(1)*XNu+CG(2)*XNd)*CG(2)*PSI(3)
      ELSE
       WRITE(6,*)'Wrong flavour subscript in 4-q matrix element'
       STOP
      ENDIF
      PSI4q = 4D0*XNa*PSI4q

      END

************************************************************************

      REAL(KIND(1D0)) FUNCTION PSI4qsl(ONE,TWO,THREE,FOUR,FIVE,SIX,IFL,X)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX,IFL,IHEL5
      COMPLEX(KIND(1D0)) qQA60,TEMP1,TEMP2
      REAL(KIND(1D0)) X,PSI,CG(2),XNf,XNu,XNd,XNc,XNa
      COMMON /QCD/XNf,XNu,XNd,XNc,XNa

      CG(1) = 0.6666666666667D0
      CG(2) =-0.3333333333333D0

      PSI = 0D0
      DO IHEL5 = 1,2

      CALL qQA60PMPMH5(ONE,TWO,THREE,FOUR,FIVE,SIX,IHEL5,X,qQA60)
      TEMP1 = qQA60
      CALL qQA60PMPMH5(THREE,FOUR,ONE,TWO,FIVE,SIX,IHEL5,X,qQA60)
      TEMP1 = TEMP1 + qQA60
      CALL qQA60PMPMH5(ONE,FOUR,THREE,TWO,FIVE,SIX,IHEL5,X,qQA60)
      TEMP2 = qQA60
      CALL qQA60PMPMH5(THREE,TWO,ONE,FOUR,FIVE,SIX,IHEL5,X,qQA60)
      TEMP2 = TEMP2 + qQA60
      PSI = PSI + TEMP1*CONJG(TEMP2) + TEMP2*CONJG(TEMP1)

      CALL qQA60MPMPH5(ONE,TWO,THREE,FOUR,FIVE,SIX,IHEL5,X,qQA60)
      TEMP1 = qQA60
      CALL qQA60MPMPH5(THREE,FOUR,ONE,TWO,FIVE,SIX,IHEL5,X,qQA60)
      TEMP1 = TEMP1 + qQA60
      CALL qQA60MPMPH5(ONE,FOUR,THREE,TWO,FIVE,SIX,IHEL5,X,qQA60)
      TEMP2 = qQA60
      CALL qQA60MPMPH5(THREE,TWO,ONE,FOUR,FIVE,SIX,IHEL5,X,qQA60)
      TEMP2 = TEMP2 + qQA60
      PSI = PSI + TEMP1*CONJG(TEMP2) + TEMP2*CONJG(TEMP1)

      ENDDO

      IF (IFL.EQ.0) THEN
       PSI4qsl = (CG(1)**2*XNu + CG(2)**2*XNd)*PSI
      ELSEIF(IFL.EQ.1) THEN
       PSI4qsl = CG(1)**2*PSI
      ELSEIF(IFL.EQ.2) THEN
       PSI4qsl = CG(2)**2*PSI
      ELSE
       WRITE(6,*)'Wrong flavour subscript in 4-q matrix element'
       STOP
      ENDIF
      PSI4qsl = 4D0*XNa*PSI4qsl

      END

************************************************************************

      SUBROUTINE qQA6tree(ONE,TWO,THREE,FOUR,FIVE,SIX,IHEL5,X,OUT)
      IMPLICIT NONE

C PPH5 -> 1+,2+,H5
C MPH5 -> 1-,2+,H5
C PMH5 -> 1+,2-,H5
C MMH5 -> 1-,2-,H5
C IHEL5 = 1 for 5+
C IHEL5 = 2 for 5-

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX,IHEL5
      COMPLEX(KIND(1D0)) OUT,A(11,11),B(11,11),XA(11,11),XB(11,11)
     >          ,qQA60PMPM,qQA60PMMP
      REAL(KIND(1D0)) X,S(11,11),XS(11,11)
      COMMON /DOTPRODUCTS/ S,XS
     >       /SPINORPRODUCTS/ A,B,XA,XB

      ENTRY qQA60PMPMH5(ONE,TWO,THREE,FOUR,FIVE,SIX,IHEL5,X,OUT)
      IF (X.EQ.1) THEN
       IF (IHEL5.EQ.2) THEN
        OUT = qQA60PMPM(ONE,TWO,THREE,FOUR,FIVE,SIX,S,A,B)
       ELSE
        OUT = qQA60PMPM(ONE,TWO,THREE,FOUR,SIX,FIVE,S,A,B)
       ENDIF
      ELSE
       IF (IHEL5.EQ.2) THEN
        OUT = qQA60PMPM(ONE,TWO,THREE,FOUR,FIVE,SIX,XS,XA,XB)
       ELSE
        OUT = qQA60PMPM(ONE,TWO,THREE,FOUR,SIX,FIVE,XS,XA,XB)
       ENDIF
      ENDIF
      RETURN

      ENTRY qQA60MPMPH5(ONE,TWO,THREE,FOUR,FIVE,SIX,IHEL5,X,OUT)
      IF (X.EQ.1) THEN
       IF (IHEL5.EQ.2) THEN
        OUT = qQA60PMPM(ONE,TWO,THREE,FOUR,SIX,FIVE,S,B,A)
       ELSE
        OUT = qQA60PMPM(ONE,TWO,THREE,FOUR,FIVE,SIX,S,B,A)
       ENDIF
      ELSE
       IF (IHEL5.EQ.2) THEN
        OUT = qQA60PMPM(ONE,TWO,THREE,FOUR,SIX,FIVE,XS,XB,XA)
       ELSE
        OUT = qQA60PMPM(ONE,TWO,THREE,FOUR,FIVE,SIX,XS,XB,XA)
       ENDIF
      ENDIF
      RETURN

      ENTRY qQA60PMMPH5(ONE,TWO,THREE,FOUR,FIVE,SIX,IHEL5,X,OUT)
      IF (X.EQ.1) THEN
       IF (IHEL5.EQ.2) THEN
        OUT = qQA60PMMP(ONE,TWO,THREE,FOUR,FIVE,SIX,S,A,B)
       ELSE
        OUT = qQA60PMMP(ONE,TWO,THREE,FOUR,SIX,FIVE,S,A,B)
       ENDIF
      ELSE
       IF (IHEL5.EQ.2) THEN
        OUT = qQA60PMMP(ONE,TWO,THREE,FOUR,FIVE,SIX,XS,XA,XB)
       ELSE
        OUT = qQA60PMMP(ONE,TWO,THREE,FOUR,SIX,FIVE,XS,XA,XB)
       ENDIF
      ENDIF
      RETURN

      ENTRY qQA60MPPMH5(ONE,TWO,THREE,FOUR,FIVE,SIX,IHEL5,X,OUT)
      IF (X.EQ.1) THEN
       IF (IHEL5.EQ.2) THEN
        OUT = qQA60PMMP(ONE,TWO,THREE,FOUR,SIX,FIVE,S,B,A)
       ELSE
        OUT = qQA60PMMP(ONE,TWO,THREE,FOUR,FIVE,SIX,S,B,A)
       ENDIF
      ELSE
       IF (IHEL5.EQ.2) THEN
        OUT = qQA60PMMP(ONE,TWO,THREE,FOUR,SIX,FIVE,XS,XB,XA)
       ELSE
        OUT = qQA60PMMP(ONE,TWO,THREE,FOUR,FIVE,SIX,XS,XB,XA)
       ENDIF
      ENDIF
      RETURN

      END

************************************************************************
************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qQA60PMPM(ONE,TWO,THREE,FOUR,FIVE,SIX,S,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      REAL(KIND(1D0)) S(11,11)
      REAL(KIND(1D0)) S34,S56,T134,T234
      COMPLEX(KIND(1D0)) B13,A52,Sm41p36,A42,B61,Sm52p43,SmILpMJ
      INTEGER I,L,M,J
      SmILpMJ(I,L,M,J) = A(I,L)*B(L,J) + A(I,M)*B(M,J)

      S34 = S(THREE,FOUR)
      S56 = S(FIVE,SIX)
      T134 = S(ONE,THREE) + S(THREE,FOUR) + S(ONE,FOUR)
      T234 = S(TWO,THREE) + S(THREE,FOUR) + S(TWO,FOUR)
      B13 = B(ONE,THREE)
      A52 = A(FIVE,TWO)
      A42 = A(FOUR,TWO)
      B61 = B(SIX,ONE)
      Sm41p36 = SmILpMJ(FOUR,ONE,THREE,SIX)
      Sm52p43 = SmILpMJ(FIVE,TWO,FOUR,THREE)

      qQA60PMPM = CMPLX(0D0,1D0)
     > *(B13*A52*Sm41p36/T134+A42*B61*Sm52p43/T234)/(S34*S56)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qQA60PMMP(ONE,TWO,THREE,FOUR,FIVE,SIX,S,A,B)
      IMPLICIT NONE

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX
      REAL(KIND(1D0)) S(11,11)
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMPLEX(KIND(1D0)) qQA60PMPM

      qQA60PMMP = qQA60PMPM(ONE,TWO,FOUR,THREE,FIVE,SIX,S,A,B)

      END

************************************************************************
