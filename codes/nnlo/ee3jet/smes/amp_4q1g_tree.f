c --- Routines taken from NT, hep-ph/9806317

c --- Throughout this file the flavor ordering
c     
c     1 - q, 2 - qb, 3 - q, 4 - qb, 5 - g, 6 - e+, 7 - e-
c     
c     is always used! 



c All amplitudes needed here are already in amp_2q2Q1g_tree.f
       
      SUBROUTINE qQgA7treeDUMMY(ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,IHELE
     $     ,OUT)
      IMPLICIT NONE
      INTEGER ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,TYPE,IHELE
      COMPLEX(KIND(1D0)) OUT

      stop "qQgA7treeDUMMY should never be called!"
      end

