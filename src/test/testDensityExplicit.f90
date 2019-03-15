PROGRAM testDensityExplicit

   USE constantProgram
   USE diffSolver

   INTEGER, PARAMETER          :: n = 7
   INTEGER                     :: i
   REAL(kind=xp), DIMENSION(n)          :: S,x,nu,sta

   CALL dtFix(0.1_xp)
   CALL dxFix(1._xp)

   DO i = 1,n
     CALL random_number(nu(i))
     CALL random_number(S(i))
     CALL random_number(sta(i))
     PRINT *, S(i)
   END DO

   x(1) = 1._xp
   DO i = 2, n
     x(i) = x(i-1) + 1._xp
   ENDDO

   CALL densityExplicit(S,nu)

   PRINT *, S

END PROGRAM testdensityExplicit
