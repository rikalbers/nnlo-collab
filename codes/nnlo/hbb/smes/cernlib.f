      REAL FUNCTION DILOG(X)                                                    
      IMPLICIT NONE
      INTEGER I
      REAL X                                                                    
      REAL(KIND(1D0)) Y,T,S,A,PI3,PI6,ZERO,ONE,HALF,MALF,MONE,MTWO             
      REAL(KIND(1D0)) C(0:18),H,ALFA,B0,B1,B2                                  
                                                                                
      DATA ZERO /0.0D0/, ONE /1.0D0/                                            
      DATA HALF /0.5D0/, MALF /-0.5D0/, MONE /-1.0D0/, MTWO /-2.0D0/            
      DATA PI3 /3.28986 81336 96453D0/, PI6 /1.64493 40668 48226D0/             
                                                                                
      DATA C( 0) / 0.42996 69356 08137 0D0/                                     
      DATA C( 1) / 0.40975 98753 30771 1D0/                                     
      DATA C( 2) /-0.01858 84366 50146 0D0/                                     
      DATA C( 3) / 0.00145 75108 40622 7D0/                                     
      DATA C( 4) /-0.00014 30418 44423 4D0/                                     
      DATA C( 5) / 0.00001 58841 55418 8D0/                                     
      DATA C( 6) /-0.00000 19078 49593 9D0/                                     
      DATA C( 7) / 0.00000 02419 51808 5D0/                                     
      DATA C( 8) /-0.00000 00319 33412 7D0/                                     
      DATA C( 9) / 0.00000 00043 45450 6D0/                                     
      DATA C(10) /-0.00000 00006 05784 8D0/                                     
      DATA C(11) / 0.00000 00000 86121 0D0/                                     
      DATA C(12) /-0.00000 00000 12443 3D0/                                     
      DATA C(13) / 0.00000 00000 01822 6D0/                                     
      DATA C(14) /-0.00000 00000 00270 1D0/                                     
      DATA C(15) / 0.00000 00000 00040 4D0/                                     
      DATA C(16) /-0.00000 00000 00006 1D0/                                     
      DATA C(17) / 0.00000 00000 00000 9D0/                                     
      DATA C(18) /-0.00000 00000 00000 1D0/                                     
                                                                                
      IF(X .EQ. ONE) THEN                                                       
       DILOG=PI6                                                                
       RETURN                                                                   
      ELSE IF(X .EQ. MONE) THEN                                                 
       DILOG=MALF*PI6                                                           
       RETURN                                                                   
      END IF                                                                    
      T=-X                                                                      
      IF(T .LE. MTWO) THEN                                                      
       Y=MONE/(ONE+T)                                                           
       S=ONE                                                                    
       A=-PI3+HALF*(LOG(-T)**2-LOG(ONE+ONE/T)**2)                               
      ELSE IF(T .LT. MONE) THEN                                                 
       Y=MONE-T                                                                 
       S=MONE                                                                   
       A=LOG(-T)                                                                
       A=-PI6+A*(A+LOG(ONE+ONE/T))                                              
      ELSE IF(T .LE. MALF) THEN                                                 
       Y=(MONE-T)/T                                                             
       S=ONE                                                                    
       A=LOG(-T)                                                                
       A=-PI6+A*(MALF*A+LOG(ONE+T))                                             
      ELSE IF(T .LT. ZERO) THEN                                                 
       Y=-T/(ONE+T)                                                             
       S=MONE                                                                   
       A=HALF*LOG(ONE+T)**2                                                     
      ELSE IF(T .LE. ONE) THEN                                                  
       Y=T                                                                      
       S=ONE                                                                    
       A=ZERO                                                                   
      ELSE                                                                      
       Y=ONE/T                                                                  
       S=MONE                                                                   
       A=PI6+HALF*LOG(T)**2                                                     
      END IF                                                                    
                                                                                
      H=Y+Y-ONE                                                                 
      ALFA=H+H                                                                  
      B1=ZERO                                                                   
      B2=ZERO                                                                   
      DO 1 I = 18,0,-1                                                          
      B0=C(I)+ALFA*B1-B2                                                        
      B2=B1                                                                     
    1 B1=B0                                                                     
      DILOG=-(S*(B0-H*B2)+A)                                                    
      RETURN                                                                    
      END                                                                       
      REAL(KIND(1D0)) FUNCTION DDILOG(X)                                       
      IMPLICIT NONE
                                                                                
      INTEGER I
      REAL(KIND(1D0)) X,Y,T,S,A,PI3,PI6,ZERO,ONE,HALF,MALF,MONE,MTWO           
      REAL(KIND(1D0)) C(0:18),H,ALFA,B0,B1,B2                                  
                                                                                
      DATA ZERO /0.0D0/, ONE /1.0D0/                                            
      DATA HALF /0.5D0/, MALF /-0.5D0/, MONE /-1.0D0/, MTWO /-2.0D0/            
      DATA PI3 /3.28986 81336 96453D0/, PI6 /1.64493 40668 48226D0/             
                                                                                
      DATA C( 0) / 0.42996 69356 08137 0D0/                                     
      DATA C( 1) / 0.40975 98753 30771 1D0/                                     
      DATA C( 2) /-0.01858 84366 50146 0D0/                                     
      DATA C( 3) / 0.00145 75108 40622 7D0/                                     
      DATA C( 4) /-0.00014 30418 44423 4D0/                                     
      DATA C( 5) / 0.00001 58841 55418 8D0/                                     
      DATA C( 6) /-0.00000 19078 49593 9D0/                                     
      DATA C( 7) / 0.00000 02419 51808 5D0/                                     
      DATA C( 8) /-0.00000 00319 33412 7D0/                                     
      DATA C( 9) / 0.00000 00043 45450 6D0/                                     
      DATA C(10) /-0.00000 00006 05784 8D0/                                     
      DATA C(11) / 0.00000 00000 86121 0D0/                                     
      DATA C(12) /-0.00000 00000 12443 3D0/                                     
      DATA C(13) / 0.00000 00000 01822 6D0/                                     
      DATA C(14) /-0.00000 00000 00270 1D0/                                     
      DATA C(15) / 0.00000 00000 00040 4D0/                                     
      DATA C(16) /-0.00000 00000 00006 1D0/                                     
      DATA C(17) / 0.00000 00000 00000 9D0/                                     
      DATA C(18) /-0.00000 00000 00000 1D0/                                     
                                                                                
      IF(X .EQ. ONE) THEN                                                       
       DDILOG=PI6                                                               
       RETURN                                                                   
      ELSE IF(X .EQ. MONE) THEN                                                 
       DDILOG=MALF*PI6                                                          
       RETURN                                                                   
      END IF                                                                    
      T=-X                                                                      
      IF(T .LE. MTWO) THEN                                                      
       Y=MONE/(ONE+T)                                                           
       S=ONE                                                                    
       A=-PI3+HALF*(LOG(-T)**2-LOG(ONE+ONE/T)**2)                               
      ELSE IF(T .LT. MONE) THEN                                                 
       Y=MONE-T                                                                 
       S=MONE                                                                   
       A=LOG(-T)                                                                
       A=-PI6+A*(A+LOG(ONE+ONE/T))                                              
      ELSE IF(T .LE. MALF) THEN                                                 
       Y=(MONE-T)/T                                                             
       S=ONE                                                                    
       A=LOG(-T)                                                                
       A=-PI6+A*(MALF*A+LOG(ONE+T))                                             
      ELSE IF(T .LT. ZERO) THEN                                                 
       Y=-T/(ONE+T)                                                             
       S=MONE                                                                   
       A=HALF*LOG(ONE+T)**2                                                     
      ELSE IF(T .LE. ONE) THEN                                                  
       Y=T                                                                      
       S=ONE                                                                    
       A=ZERO                                                                   
      ELSE                                                                      
       Y=ONE/T                                                                  
       S=MONE                                                                   
       A=PI6+HALF*LOG(T)**2                                                     
      END IF                                                                    
                                                                                
      H=Y+Y-ONE                                                                 
      ALFA=H+H                                                                  
      B1=ZERO                                                                   
      B2=ZERO                                                                   
      DO 1 I = 18,0,-1                                                          
      B0=C(I)+ALFA*B1-B2                                                        
      B2=B1                                                                     
    1 B1=B0                                                                     
      DDILOG=-(S*(B0-H*B2)+A)                                                   
      RETURN                                                                    
      END                                                                       

      REAL(KIND(1D0)) FUNCTION DCLAUS(X)
      IMPLICIT NONE
 
      REAL(KIND(1D0)) A,B,X
      DIMENSION A(0:8),B(0:13)
      INTEGER I
 
      REAL(KIND(1D0)) B0,B1,B2,ALFA,H,U,S,V,R1,HF,PI,PI2,PIH,RPIH
      PARAMETER (R1 = 1, HF =R1/2)
      PARAMETER (PI = 3.14159 26535 89793 24D0)
      PARAMETER (PI2 = 2*PI, PIH = PI/2, RPIH = 2/PI)
 
      DATA A( 0) / 0.02795 28319 73575 6613D0/
      DATA A( 1) / 0.00017 63088 74389 8116D0/
      DATA A( 2) / 0.00000 12662 74146 1157D0/
      DATA A( 3) / 0.00000 00117 17181 8134D0/
      DATA A( 4) / 0.00000 00001 23006 4129D0/
      DATA A( 5) / 0.00000 00000 01395 2729D0/
      DATA A( 6) / 0.00000 00000 00016 6908D0/
      DATA A( 7) / 0.00000 00000 00000 2076D0/
      DATA A( 8) / 0.00000 00000 00000 0027D0/
 
      DATA B( 0) / 0.63909 70888 57265 341D0/
      DATA B( 1) /-0.05498 05693 01851 716D0/
      DATA B( 2) /-0.00096 12619 45950 606D0/
      DATA B( 3) /-0.00003 20546 86822 550D0/
      DATA B( 4) /-0.00000 13294 61695 426D0/
      DATA B( 5) /-0.00000 00620 93601 824D0/
      DATA B( 6) /-0.00000 00031 29600 656D0/
      DATA B( 7) /-0.00000 00001 66351 954D0/
      DATA B( 8) /-0.00000 00000 09196 527D0/
      DATA B( 9) /-0.00000 00000 00524 004D0/
      DATA B(10) /-0.00000 00000 00030 580D0/
      DATA B(11) /-0.00000 00000 00001 820D0/
      DATA B(12) /-0.00000 00000 00000 110D0/
      DATA B(13) /-0.00000 00000 00000 007D0/
 
      V=MOD(ABS(X),PI2)
      S=SIGN(R1,X)
      IF(V .GT. PI) THEN
       V=PI2-V
       S=-S
      ENDIF
      IF(V .EQ. 0 .OR. V .EQ. PI) THEN
       H=0
      ELSEIF(V .LT. PIH) THEN
       U=RPIH*V
       H=2*U**2-1
       ALFA=H+H
       B1=0
       B2=0
       DO 1 I = 8,0,-1
       B0=A(I)+ALFA*B1-B2
       B2=B1
    1  B1=B0
       H=V*(1-LOG(V)+HF*V**2*(B0-H*B2))
      ELSE
       U=RPIH*V-2
       H=2*U**2-1
       ALFA=H+H
       B1=0
       B2=0
       DO 2 I = 13,0,-1
       B0=B(I)+ALFA*B1-B2
       B2=B1
    2  B1=B0
       H=(PI-V)*(B0-H*B2)
      ENDIF
      DCLAUS=S*H
      RETURN
      END
