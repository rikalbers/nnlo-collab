c --- Routines taken from NT, hep-ph/9806317

c --- Throughout this file the flavor ordering
c     
c     1 - q, 2 - qb, 3 - Q, 4 - Qb, 5 - g, 6 - e+, 7 - e-
c     
c     is always used! 


      SUBROUTINE qQgA7tree(ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,IHELE,OUT)
      IMPLICIT NONE

C PPPHe -> 1+,3+,5+,He
C PPMHe -> 1+,3+,5-,He
C PMPHe -> 1+,3-,5+,He
C PMMHe -> 1+,3-,5-,He
C MPPHe -> 1-,3+,5+,He
C MPMHe -> 1-,3+,5-,He
C MMPHe -> 1-,3-,5+,He
C MMMHe -> 1-,3-,5-,He
C IHELE = 1 for He+
C IHELE = 2 for He-

      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,IHELE
      COMPLEX(KIND(1D0)) OUT
      COMPLEX(KIND(1D0)) qQgA70PPP_A1,qQgA70PPM_A1,qQgA70PMP_A1,
     >                   qQgA70PMM_A1
      COMPLEX(KIND(1D0)) qQgA70PPP_A2,qQgA70PPM_A2,qQgA70PMP_A2,
     >                   qQgA70PMM_A2
      COMPLEX(KIND(1D0)) qQgA70PPP_A3,qQgA70PMM_A3
      COMPLEX(KIND(1D0)) qQgA70PPM_A4,qQgA70PMP_A4
      REAL(KIND(1D0)) S(11,11)
      COMPLEX(KIND(1D0)) A(11,11),B(11,11)
      COMMON /DOTPRODUCTS/S
      COMMON /SPINORPRODUCTS/A,B 

 
      ENTRY qQgA70PPPHe_A1(ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,IHELE,OUT)
      IF (IHELE.EQ.2) THEN
         OUT = qQgA70PPP_A1(ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,A,B)
      ELSE
         OUT = qQgA70PPP_A1(ONE,TWO,THREE,FOUR,FIVE,SEVEN,SIX,A,B)
      ENDIF
      RETURN
      ENTRY qQgA70MMMHe_A1(ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,IHELE,OUT)
      IF (IHELE.EQ.2) THEN
         OUT = -qQgA70PPP_A1(ONE,TWO,THREE,FOUR,FIVE,SEVEN,SIX,B,A)
      ELSE
         OUT = -qQgA70PPP_A1(ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,B,A)
      ENDIF
      RETURN
      ENTRY qQgA70PPMHe_A1(ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,IHELE,OUT)
      IF (IHELE.EQ.2) THEN
         OUT = qQgA70PPM_A1(ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,A,B)
      ELSE
         OUT = qQgA70PPM_A1(ONE,TWO,THREE,FOUR,FIVE,SEVEN,SIX,A,B)
      ENDIF
      RETURN
      ENTRY qQgA70MMPHe_A1(ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,IHELE,OUT)
      IF (IHELE.EQ.2) THEN
         OUT = -qQgA70PPM_A1(ONE,TWO,THREE,FOUR,FIVE,SEVEN,SIX,B,A)
      ELSE
         OUT = -qQgA70PPM_A1(ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,B,A)
      ENDIF
      RETURN
      ENTRY qQgA70PMPHe_A1(ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,IHELE,OUT)
      IF (IHELE.EQ.2) THEN
         OUT = qQgA70PMP_A1(ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,A,B)
      ELSE
         OUT = qQgA70PMP_A1(ONE,TWO,THREE,FOUR,FIVE,SEVEN,SIX,A,B)
      ENDIF
      RETURN
      ENTRY qQgA70MPMHe_A1(ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,IHELE,OUT)
      IF (IHELE.EQ.2) THEN
         OUT = -qQgA70PMP_A1(ONE,TWO,THREE,FOUR,FIVE,SEVEN,SIX,B,A)
      ELSE
         OUT = -qQgA70PMP_A1(ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,B,A)
      ENDIF
      RETURN
      ENTRY qQgA70PMMHe_A1(ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,IHELE,OUT)
      IF (IHELE.EQ.2) THEN
         OUT = qQgA70PMM_A1(ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,A,B)
      ELSE
         OUT = qQgA70PMM_A1(ONE,TWO,THREE,FOUR,FIVE,SEVEN,SIX,A,B)
      ENDIF
      RETURN
      ENTRY qQgA70MPPHe_A1(ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,IHELE,OUT)
      IF (IHELE.EQ.2) THEN
         OUT = -qQgA70PMM_A1(ONE,TWO,THREE,FOUR,FIVE,SEVEN,SIX,B,A)
      ELSE
         OUT = -qQgA70PMM_A1(ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,B,A)
      ENDIF
      RETURN
      
      ENTRY qQgA70PPPHe_A2(ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,IHELE,OUT)
      IF (IHELE.EQ.2) THEN
         OUT = qQgA70PPP_A2(ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,A,B)
      ELSE
         OUT = qQgA70PPP_A2(ONE,TWO,THREE,FOUR,FIVE,SEVEN,SIX,A,B)
      ENDIF
      RETURN
      ENTRY qQgA70MMMHe_A2(ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,IHELE,OUT)
      IF (IHELE.EQ.2) THEN
         OUT = -qQgA70PPP_A2(ONE,TWO,THREE,FOUR,FIVE,SEVEN,SIX,B,A)
      ELSE
         OUT = -qQgA70PPP_A2(ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,B,A)
      ENDIF
      RETURN
      ENTRY qQgA70PPMHe_A2(ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,IHELE,OUT)
      IF (IHELE.EQ.2) THEN
         OUT = qQgA70PPM_A2(ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,A,B)
      ELSE
         OUT = qQgA70PPM_A2(ONE,TWO,THREE,FOUR,FIVE,SEVEN,SIX,A,B)
      ENDIF
      RETURN
      ENTRY qQgA70MMPHe_A2(ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,IHELE,OUT)
      IF (IHELE.EQ.2) THEN
         OUT = -qQgA70PPM_A2(ONE,TWO,THREE,FOUR,FIVE,SEVEN,SIX,B,A)
      ELSE
         OUT = -qQgA70PPM_A2(ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,B,A)
      ENDIF
      RETURN
      ENTRY qQgA70PMPHe_A2(ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,IHELE,OUT)
      IF (IHELE.EQ.2) THEN
         OUT = qQgA70PMP_A2(ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,A,B)
      ELSE
         OUT = qQgA70PMP_A2(ONE,TWO,THREE,FOUR,FIVE,SEVEN,SIX,A,B)
      ENDIF
      RETURN
      ENTRY qQgA70MPMHe_A2(ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,IHELE,OUT)
      IF (IHELE.EQ.2) THEN
         OUT = -qQgA70PMP_A2(ONE,TWO,THREE,FOUR,FIVE,SEVEN,SIX,B,A)
      ELSE
         OUT = -qQgA70PMP_A2(ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,B,A)
      ENDIF
      RETURN
      ENTRY qQgA70PMMHe_A2(ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,IHELE,OUT)
      IF (IHELE.EQ.2) THEN
         OUT = qQgA70PMM_A2(ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,A,B)
      ELSE
         OUT = qQgA70PMM_A2(ONE,TWO,THREE,FOUR,FIVE,SEVEN,SIX,A,B)
      ENDIF
      RETURN
      ENTRY qQgA70MPPHe_A2(ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,IHELE,OUT)
      IF (IHELE.EQ.2) THEN
         OUT = -qQgA70PMM_A2(ONE,TWO,THREE,FOUR,FIVE,SEVEN,SIX,B,A)
      ELSE
         OUT = -qQgA70PMM_A2(ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,B,A)
      ENDIF
      RETURN


      ENTRY qQgA70PPPHe_A3(ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,IHELE,OUT)
      IF (IHELE.EQ.2) THEN
         OUT = qQgA70PPP_A3(ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,A,B)
      ELSE
         OUT = qQgA70PPP_A3(ONE,TWO,THREE,FOUR,FIVE,SEVEN,SIX,A,B)
      ENDIF
      RETURN
      ENTRY qQgA70MMMHe_A3(ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,IHELE,OUT)
      IF (IHELE.EQ.2) THEN
         OUT = -qQgA70PPP_A3(ONE,TWO,THREE,FOUR,FIVE,SEVEN,SIX,B,A)
      ELSE
         OUT = -qQgA70PPP_A3(ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,B,A)
      ENDIF
      RETURN
      ENTRY qQgA70PMMHe_A3(ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,IHELE,OUT)
      IF (IHELE.EQ.2) THEN
         OUT = qQgA70PMM_A3(ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,A,B)
      ELSE
         OUT = qQgA70PMM_A3(ONE,TWO,THREE,FOUR,FIVE,SEVEN,SIX,A,B)
      ENDIF
      RETURN
      ENTRY qQgA70MPPHe_A3(ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,IHELE,OUT)
      IF (IHELE.EQ.2) THEN
         OUT = -qQgA70PMM_A3(ONE,TWO,THREE,FOUR,FIVE,SEVEN,SIX,B,A)
      ELSE
         OUT = -qQgA70PMM_A3(ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,B,A)
      ENDIF
      RETURN
      ENTRY qQgA70PPMHe_A3(ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,IHELE,OUT)
      out = (0D0,0D0)
      return
      ENTRY qQgA70MMPHe_A3(ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,IHELE,OUT)
      out = (0D0,0D0)
      return
      ENTRY qQgA70PMPHe_A3(ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,IHELE,OUT)
      out = (0D0,0D0)
      return
      ENTRY qQgA70MPMHe_A3(ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,IHELE,OUT)
      out = (0D0,0D0)
      return

      ENTRY qQgA70PPMHe_A4(ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,IHELE,OUT)
      IF (IHELE.EQ.2) THEN
         OUT = qQgA70PPM_A4(ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,A,B)
      ELSE
         OUT = qQgA70PPM_A4(ONE,TWO,THREE,FOUR,FIVE,SEVEN,SIX,A,B)
      ENDIF
      RETURN
      ENTRY qQgA70MMPHe_A4(ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,IHELE,OUT)
      IF (IHELE.EQ.2) THEN
         OUT = -qQgA70PPM_A4(ONE,TWO,THREE,FOUR,FIVE,SEVEN,SIX,B,A)
      ELSE
         OUT = -qQgA70PPM_A4(ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,B,A)
      ENDIF
      RETURN
      ENTRY qQgA70PMPHe_A4(ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,IHELE,OUT)
      IF (IHELE.EQ.2) THEN
         OUT = qQgA70PMP_A4(ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,A,B)
      ELSE
         OUT = qQgA70PMP_A4(ONE,TWO,THREE,FOUR,FIVE,SEVEN,SIX,A,B)
      ENDIF
      RETURN
      ENTRY qQgA70MPMHe_A4(ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,IHELE,OUT)
      IF (IHELE.EQ.2) THEN
         OUT = -qQgA70PMP_A4(ONE,TWO,THREE,FOUR,FIVE,SEVEN,SIX,B,A)
      ELSE
         OUT = -qQgA70PMP_A4(ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,B,A)
      ENDIF
      RETURN
      ENTRY qQgA70PPPHe_A4(ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,IHELE,OUT)
      out = (0D0,0D0) 
      RETURN
      ENTRY qQgA70MMMHe_A4(ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,IHELE,OUT)
      out = (0D0,0D0) 
      RETURN
      ENTRY qQgA70PMMHe_A4(ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,IHELE,OUT)
      out = (0D0,0D0) 
      RETURN
      ENTRY qQgA70MPPHe_A4(ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,IHELE,OUT)
      out = (0D0,0D0) 
      RETURN
      END

************************************************************************

C Routines below translated form NLOJET++ 4.1.3

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qQgA70PPP_A1(ONE,TWO,THREE,FOUR,FIVE,
     $                                         SIX,SEVEN,Ain,Bin)
      IMPLICIT NONE
      COMPLEX(KIND(1D0)) ain(11,11),bin(11,11)
      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,I,J,K,L

      REAL(KIND(1D0)) S15,s34,S3,t267,t234,t167,t345
      COMPLEX(KIND(1D0)) b15,a62,a45,b17,a42,b23,a15,a54,b53,a34,
     >               b35
      COMPLEX(KIND(1D0)) C,c4153,c4267,c6175,c4157,c6243,c4351,c6173
      REAL(KIND(1D0)) S(11,11)
      COMMON /DOTPRODUCTS/S

      S3(I,J,K) = S(I,J) + S(I,K) + S(J,K)
      C(I,J,K,L) = Ain(I,J)*Bin(J,L) + Ain(I,K)*Bin(K,L)

      s15 = S(ONE,FIVE)
      s34 = S(THREE,FOUR)

      t267 = S3(TWO,SIX,SEVEN)
      t234 = S3(TWO,THREE,FOUR)
      t167 = S3(ONE,SIX,SEVEN)
      t345 = S3(THREE,FOUR,FIVE)
      
      b15 = Bin(ONE,FIVE)
      a62 = Ain(SIX,TWO)
      a45 = Ain(FOUR,FIVE)
      b17 = Bin(ONE,SEVEN)
      a42 = Ain(FOUR,TWO)
      b23 = Bin(TWO,THREE)
      a15 = Ain(ONE,FIVE)
      a54 = Ain(FIVE,FOUR)
      b53 = Bin(FIVE,THREE)
      a34 = Ain(THREE,FOUR)
      b35 = Bin(THREE,FIVE)

      c4153 = C(FOUR,ONE,FIVE,THREE)
      c4267 = C(FOUR,TWO,SIX,SEVEN)
      c6175 = C(SIX,ONE,SEVEN,FIVE) 
      c4157 = C(FOUR,ONE,FIVE,SEVEN) 
      c6243 = C(SIX,TWO,FOUR,THREE)
      c4351 = C(FOUR,THREE,FIVE,ONE) 
      c6173 = C(SIX,ONE,SEVEN,THREE)
    
      qQgA70PPP_A1 = - b15*c4153*c4267*a62/(a45*s15*s34*t267)
     > - b17*c6175*a42*a42*b23/(a45*s34*t234*t167)
     > - c4157*c6243*a42/(a15*a54*s34*t234)
     > + b53*c4351*c4267*a62/(a45*s34*t345*t267)
     > + b17*(c6173*a34 + c6175*a54)*b35*a42/(a45*s34*t345*t167)

      qQgA70PPP_A1 = (0,1)*qQgA70PPP_A1/S(SIX,SEVEN)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qQgA70PPM_A1(ONE,TWO,THREE,FOUR,FIVE,
     $                                         SIX,SEVEN,Ain,Bin)
      IMPLICIT NONE
      COMPLEX(KIND(1D0)) ain(11,11),bin(11,11)
      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,I,J,K,L
      REAL(KIND(1D0)) S3,s15,s34,t267,t234,t167,t345
      COMPLEX(KIND(1D0))  b13,a51,a62,b35,b17,a42,b15,b53,a54,b34
      COMPLEX(KIND(1D0)) C,c4267,c6243,c6173,c5243,c5267,c2453
      REAL(KIND(1D0)) S(11,11)
      COMMON /DOTPRODUCTS/S

      S3(I,J,K) = S(I,J) + S(I,K) + S(J,K)
      C(I,J,K,L) = Ain(I,J)*Bin(J,L) + Ain(I,K)*Bin(K,L)
 
      s15 = S(ONE,FIVE)
      s34 = S(THREE,FOUR)
      t267 = S3(TWO,SIX,SEVEN)
      t234 = S3(TWO,THREE,FOUR)
      t167 = S3(ONE,SIX,SEVEN)
      t345 = S3(THREE,FOUR,FIVE)
    
      b13 = Bin(ONE,THREE)
      a51 = Ain(FIVE,ONE)
      a62 = Ain(SIX,TWO)
      b35 = Bin(THREE,FIVE) 
      b17 = Bin(ONE,SEVEN)
      a42 = Ain(FOUR,TWO)
      b15 = Bin(ONE,FIVE)
      b53 = Bin(FIVE,THREE)
      a54 = Ain(FIVE,FOUR) 
      b34 = Bin(THREE,FOUR)
    
      c4267 = C(FOUR,TWO,SIX,SEVEN)
      c6243 = C(SIX,TWO,FOUR,THREE)
      c6173 = C(SIX,ONE,SEVEN,THREE) 
      c5243 = C(FIVE,TWO,FOUR,THREE)
      c5267 = C(FIVE,TWO,SIX,SEVEN)
      c2453 = C(TWO,FOUR,FIVE,THREE)  
    
      qQgA70PPM_A1 = b13*b13*a51*c4267*a62/(b35*s15*s34*t267)
     > + b17*c6173*c5243*a42/(b35*s34*t234*t167)
     > - b13*b17*c6243*a42/(b15*b53*s34*t234)
     > + b13*a54*(b34*c4267 + b35*c5267)*a62/(b35*s34*t345*t267)
     > - b17*c6173*a54*c2453/(b35*s34*t345*t167)
    
      qQgA70PPM_A1 = (0,1)*qQgA70PPM_A1/S(SIX,SEVEN)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qQgA70PMP_A1(ONE,TWO,THREE,FOUR,FIVE,
     $                                         SIX,SEVEN,Ain,Bin)
      IMPLICIT NONE
      COMPLEX(KIND(1D0)) ain(11,11),bin(11,11)
      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN
      COMPLEX(KIND(1D0)) qQgA70PPP_A1 
      REAL(KIND(1D0)) S(11,11)
      COMMON /DOTPRODUCTS/S
      

      qQgA70PMP_A1=qQgA70PPP_A1(ONE,TWO,FOUR,THREE,FIVE,SIX
     $     ,SEVEN,Ain,Bin)
      end

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qQgA70PMM_A1(ONE,TWO,THREE,FOUR,FIVE,
     $                                         SIX,SEVEN,Ain,Bin)
      IMPLICIT NONE
      COMPLEX(KIND(1D0)) ain(11,11),bin(11,11)
      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN
      COMPLEX(KIND(1D0)) qQgA70PPM_A1 
      REAL(KIND(1D0)) S(11,11)
      COMMON /DOTPRODUCTS/S
      

      qQgA70PMM_A1=qQgA70PPM_A1(ONE,TWO,FOUR,THREE,FIVE,SIX
     $     ,SEVEN,Ain,Bin)
      
      end

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qQgA70PPP_A2(ONE,TWO,THREE,FOUR,FIVE,
     $                                         SIX,SEVEN,Ain,Bin)
      IMPLICIT NONE
      COMPLEX(KIND(1D0)) ain(11,11),bin(11,11)
      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,I,J,K,L
      REAL(KIND(1D0)) S3,s34,s25,t345,t267,t167,t134
      COMPLEX(KIND(1D0)) b53,a62,a45,b17,a34,a54,b35,a42,b13,a52,b25
      COMPLEX(KIND(1D0)) C,c4351,c4267,c6173,c6175,c4135,c4137
      REAL(KIND(1D0)) S(11,11)
      COMMON /DOTPRODUCTS/S

      S3(I,J,K) = S(I,J) + S(I,K) + S(J,K)
      C(I,J,K,L) = Ain(I,J)*Bin(J,L) + Ain(I,K)*Bin(K,L)

      s34 = S(THREE,FOUR)
      s25 = S(TWO,FIVE)
      t345 = S3(THREE,FOUR,FIVE)
      t267 = S3(TWO,SIX,SEVEN)
      t167 = S3(ONE,SIX,SEVEN)
      t134 = S3(ONE,THREE,FOUR)
    
      b53 = Bin(FIVE,THREE)
      a62 = Ain(SIX,TWO)
      a45 = Ain(FOUR,FIVE)
      b17 = Bin(ONE,SEVEN)
      a34 = Ain(THREE,FOUR)
      a54 = Ain(FIVE,FOUR)
      b35 = Bin(THREE,FIVE)
      a42 = Ain(FOUR,TWO)
      b13 = Bin(ONE,THREE)
      a52 = Ain(FIVE,TWO)
      b25 = Bin(TWO,FIVE)
    
      c4351 = C(FOUR,THREE,FIVE,ONE)
      c4267 = C(FOUR,TWO,SIX,SEVEN)
      c6173 = C(SIX,ONE,SEVEN,THREE) 
      c6175 = C(SIX,ONE,SEVEN,FIVE)
      c4135 = C(FOUR,ONE,THREE,FIVE)
      c4137 = C(FOUR,ONE,THREE,SEVEN)
    
      qQgA70PPP_A2 =- b53*c4351*c4267*a62/(a45*s34*t345*t267)
     > - b17*(c6173*a34 + c6175*a54)*b35*a42/(a45*s34*t345*t167)
     > - b13*c4135*c4267*a62/(a45*s34*t134*t267)
     > - b13*c4137*a62*a42/(a45*a52*s34*t134)
     > - b17*c6173*b25*a42*a42/(a45*s25*s34*t167)

      qQgA70PPP_A2 = (0,1)*qQgA70PPP_A2/S(SIX,SEVEN)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qQgA70PPM_A2(ONE,TWO,THREE,FOUR,FIVE,
     $                                         SIX,SEVEN,Ain,Bin)
      IMPLICIT NONE
      COMPLEX(KIND(1D0)) ain(11,11),bin(11,11)
      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,I,J,K,L
      REAL(KIND(1D0)) S3,s34,s25,t345,t267,t167,t134
      COMPLEX(KIND(1D0)) b13,a54,b34,b35,a62,b17,a41,b52,a52 
      COMPLEX(KIND(1D0)) C,c4267,c5267,c6173,c2453,c4137,c6253,c4253 
      REAL(KIND(1D0)) S(11,11)
      COMMON /DOTPRODUCTS/S

      S3(I,J,K) = S(I,J) + S(I,K) + S(J,K)
      C(I,J,K,L) = Ain(I,J)*Bin(J,L) + Ain(I,K)*Bin(K,L)

      s34 = S(THREE,FOUR)
      s25 = S(TWO,FIVE)
      t345 = S3(THREE,FOUR,FIVE)
      t267 = S3(TWO,SIX,SEVEN)
      t167 = S3(ONE,SIX,SEVEN)
      t134 = S3(ONE,THREE,FOUR)
    
      b13 = Bin(ONE,THREE)
      a54 = Ain(FIVE,FOUR)
      b34 = Bin(THREE,FOUR)
      b35 = Bin(THREE,FIVE) 
      a62 = Ain(SIX,TWO)
      b17 = Bin(ONE,SEVEN)
      a41 = Ain(FOUR,ONE)
      b52 = Bin(FIVE,TWO)
      a52 = Ain(FIVE,TWO)
    
      c4267 = C(FOUR,TWO,SIX,SEVEN)
      c5267 = C(FIVE,TWO,SIX,SEVEN)
      c6173 = C(SIX,ONE,SEVEN,THREE)
      c2453 = C(TWO,FOUR,FIVE,THREE)
      c4137 = C(FOUR,ONE,THREE,SEVEN)
      c6253 = C(SIX,TWO,FIVE,THREE)
      c4253 = C(FOUR,TWO,FIVE,THREE)
    
      qQgA70PPM_A2 = 
     > - b13*a54*(b34*c4267 + b35*c5267)*a62/(b35*s34*t345*t267)
     > + b17*c6173*a54*c2453/(b35*s34*t345*t167)
     > + b13*b13*a41*c5267*a62/(b35*s34*t134*t267)
     > - b13*c4137*c6253/(b35*b52*s34*t134)
     > + b17*c6173*a52*c4253/(b35*s25*s34*t167)
    
      qQgA70PPM_A2 = (0,1)*qQgA70PPM_A2/S(SIX,SEVEN)

      END


************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qQgA70PMP_A2(ONE,TWO,THREE,FOUR,FIVE,
     $                                         SIX,SEVEN,Ain,Bin)
      IMPLICIT NONE
      COMPLEX(KIND(1D0)) ain(11,11),bin(11,11)
      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN
      COMPLEX(KIND(1D0)) qQgA70PPP_A2 
      REAL(KIND(1D0)) S(11,11)
      COMMON /DOTPRODUCTS/S
      

      qQgA70PMP_A2=qQgA70PPP_A2(ONE,TWO,FOUR,THREE,FIVE,SIX
     $     ,SEVEN,Ain,Bin)
      end

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qQgA70PMM_A2(ONE,TWO,THREE,FOUR,FIVE,
     $                                         SIX,SEVEN,Ain,Bin)
      IMPLICIT NONE
      COMPLEX(KIND(1D0)) ain(11,11),bin(11,11)
      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN
      COMPLEX(KIND(1D0)) qQgA70PPM_A2 
      REAL(KIND(1D0)) S(11,11)
      COMMON /DOTPRODUCTS/S
      

      qQgA70PMM_A2=qQgA70PPM_A2(ONE,TWO,FOUR,THREE,FIVE,SIX 
     $     ,SEVEN,Ain,Bin)
      
      end

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qQgA70PPP_A3(ONE,TWO,THREE,FOUR,FIVE,
     $                                         SIX,SEVEN,Ain,Bin)
      IMPLICIT NONE
      COMPLEX(KIND(1D0)) ain(11,11),bin(11,11)
      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,I,J,K,L
      REAL(KIND(1D0)) S3,s35,t345,t267,t167
      COMPLEX(KIND(1D0)) b53,a62,a45,b17,a34,a54,b35,a42
      COMPLEX(KIND(1D0)) C,c6175,c6173,c4351,c4267
      REAL(KIND(1D0)) S(11,11)
      COMMON /DOTPRODUCTS/S

      S3(I,J,K) = S(I,J) + S(I,K) + S(J,K)
      C(I,J,K,L) = Ain(I,J)*Bin(J,L) + Ain(I,K)*Bin(K,L)

      s35 = S(THREE,FIVE)
      t345 = S3(THREE,FOUR,FIVE)
      t267 = S3(TWO,SIX,SEVEN)
      t167 = S3(ONE,SIX,SEVEN)
    
      b53 = Bin(FIVE,THREE)
      a62 = Ain(SIX,TWO)
      a45 = Ain(FOUR,FIVE)
      b17 = Bin(ONE,SEVEN) 
      a34 = Ain(THREE,FOUR)
      a54 = Ain(FIVE,FOUR)
      b35 = Bin(THREE,FIVE)
      a42 = Ain(FOUR,TWO) 
    
      c6175 = C(SIX,ONE,SEVEN,FIVE)
      c6173 = C(SIX,ONE,SEVEN,THREE)
      c4351 = C(FOUR,THREE,FIVE,ONE)
      c4267 = C(FOUR,TWO,SIX,SEVEN)
    
      qQgA70PPP_A3 =- b53*c4351*c4267*a62/(a45*s35*t345*t267)
     > - b17*(c6173*a34 + c6175*a54)*b35*a42/(a45*s35*t345*t167)    
      
      qQgA70PPP_A3 = (0,1)*qQgA70PPP_A3/S(SIX,SEVEN)
      
      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qQgA70PMM_A3(ONE,TWO,THREE,FOUR,FIVE,
     $                                         SIX,SEVEN,Ain,Bin)
      IMPLICIT NONE
      COMPLEX(KIND(1D0)) ain(11,11),bin(11,11)
      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,I,J,K,L
      REAL(KIND(1D0)) S3,s35,t435,t267,t167
      COMPLEX(KIND(1D0)) b14,a53,b43,b45,a62,b17
      COMPLEX(KIND(1D0)) C,c3267,c5267,c6174,c2354
      REAL(KIND(1D0)) S(11,11)
      COMMON /DOTPRODUCTS/S

      S3(I,J,K) = S(I,J) + S(I,K) + S(J,K)
      C(I,J,K,L) = Ain(I,J)*Bin(J,L) + Ain(I,K)*Bin(K,L)

      s35 = S(THREE,FIVE)
      t435 = S3(FOUR,THREE,FIVE)
      t267 = S3(TWO,SIX,SEVEN)
      t167 = S3(ONE,SIX,SEVEN)
    
      b14 = Bin(ONE,FOUR)
      a53 = Ain(FIVE,THREE)
      b43 = Bin(FOUR,THREE)
      b45 = Bin(FOUR,FIVE)
      a62 = Ain(SIX,TWO)
      b17 = Bin(ONE,SEVEN)
    
      c3267 = C(THREE,TWO,SIX,SEVEN)
      c5267 = C(FIVE,TWO,SIX,SEVEN)
      c6174 = C(SIX,ONE,SEVEN,FOUR) 
      c2354 = C(TWO,THREE,FIVE,FOUR)
    
      qQgA70PMM_A3 = 
     > - b14*a53*(b43*c3267 + b45*c5267)*a62/(b45*s35*t435*t267)
     > + b17*c6174*a53*c2354/(b45*s35*t435*t167)
 
      qQgA70PMM_A3 = (0,1)*qQgA70PMM_A3/S(SIX,SEVEN)

      END


************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qQgA70PPM_A4(ONE,TWO,THREE,FOUR,FIVE,
     $                                         SIX,SEVEN,Ain,Bin)
      IMPLICIT NONE
      COMPLEX(KIND(1D0)) ain(11,11),bin(11,11)
      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,I,J,K,L
      REAL(KIND(1D0)) S3,s45,t345,t267,t167
      COMPLEX(KIND(1D0)) b13,a54,b34,b35,a62,b17
      COMPLEX(KIND(1D0)) C,c4267,c5267,c6173,c2453
      REAL(KIND(1D0)) S(11,11)
      COMMON /DOTPRODUCTS/S

      S3(I,J,K) = S(I,J) + S(I,K) + S(J,K)
      C(I,J,K,L) = Ain(I,J)*Bin(J,L) + Ain(I,K)*Bin(K,L)

      s45 = S(FOUR,FIVE)
      t345 = S3(THREE,FOUR,FIVE)
      t267 = S3(TWO,SIX,SEVEN)
      t167 = S3(ONE,SIX,SEVEN)
    
      b13 = Bin(ONE,THREE)
      a54 = Ain(FIVE,FOUR)
      b34 = Bin(THREE,FOUR)
      b35 = Bin(THREE,FIVE) 
      a62 = Ain(SIX,TWO)
      b17 = Bin(ONE,SEVEN) 
    
      c4267 = C(FOUR,TWO,SIX,SEVEN)
      c5267 = C(FIVE,TWO,SIX,SEVEN)
      c6173 = C(SIX,ONE,SEVEN,THREE) 
      c2453 = C(TWO,FOUR,FIVE,THREE)
    
      qQgA70PPM_A4 = 
     > b13*a54*(b34*c4267 + b35*c5267)*a62/(b35*s45*t345*t267)
     > - b17*c6173*a54*c2453/(b35*s45*t345*t167)
    
      qQgA70PPM_A4 = (0,1)*qQgA70PPM_A4/S(SIX,SEVEN)

      END

************************************************************************

      COMPLEX(KIND(1D0)) FUNCTION qQgA70PMP_A4(ONE,TWO,THREE,FOUR,FIVE,
     $                                         SIX,SEVEN,Ain,Bin)
      IMPLICIT NONE
      COMPLEX(KIND(1D0)) ain(11,11),bin(11,11)
      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,I,J,K,L
      REAL(KIND(1D0)) S3,s45,t345,t267,t167
      COMPLEX(KIND(1D0)) b54,a62,a35,b17,a43,a53,b45,a32
      COMPLEX(KIND(1D0)) C,c6175,c6174,c3451,c3267 
      REAL(KIND(1D0)) S(11,11)
      COMMON /DOTPRODUCTS/S

      S3(I,J,K) = S(I,J) + S(I,K) + S(J,K)
      C(I,J,K,L) = Ain(I,J)*Bin(J,L) + Ain(I,K)*Bin(K,L)

      s45 = S(FOUR,FIVE)
      t345 = S3(THREE,FOUR,FIVE)
      t267 = S3(TWO,SIX,SEVEN)
      t167 = S3(ONE,SIX,SEVEN)
    
      b54 = Bin(FIVE,FOUR)
      a62 = Ain(SIX,TWO)
      a35 = Ain(THREE,FIVE)
      b17 = Bin(ONE,SEVEN)
      a43 = Ain(FOUR,THREE)
      a53 = Ain(FIVE,THREE)
      b45 = Bin(FOUR,FIVE)
      a32 = Ain(THREE,TWO) 
    
      c6175 = C(SIX,ONE,SEVEN,FIVE)
      c6174 = C(SIX,ONE,SEVEN,FOUR)
      c3451 = C(THREE,FOUR,FIVE,ONE) 
      c3267 = C(THREE,TWO,SIX,SEVEN)
    
      qQgA70PMP_A4 = b54*c3451*c3267*a62/(a35*s45*t345*t267)
     > + b17*(c6174*a43 + c6175*a53)*b45*a32/(a35*s45*t345*t167)
 
      qQgA70PMP_A4 = (0,1)*qQgA70PMP_A4/S(SIX,SEVEN)

      END
