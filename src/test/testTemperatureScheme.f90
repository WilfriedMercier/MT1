PROGRAM testTemperatureScheme

   USE constantProgram
   USE diffSolver

   INTEGER, PARAMETER          :: n = 7
   INTEGER                     :: i
   REAL(KIND=xp), DIMENSION(n) :: T1, T2, Qadv, Qp, Qm, Cv

   DO i = 1,n
      CALL RANDOM_NUMBER(Qadv)
      CALL RANDOM_NUMBER(Qp)
      CALL RANDOM_NUMBER(Qm)
      CALL RANDOM_NUMBER(Cv)
   END DO

   Qp(1)=0_xp
   Qm(1)=0_xp
   Qadv(1)=0_xp

   Qp(2)=1_xp
   Qm(2)=1_xp
   Qadv(2)=0_xp

   Qp(3)=1.1_xp
   Qm(3)=1_xp
   Qadv(3)=0_xp

   CALL dtFix(0.1_xp)

   CALL temperatureExplicitScheme(T1, Qadv, Qp, Qm, Cv)
   CALL temperatureImplicitScheme(T2, Qadv, Qp, Qm, Cv)

   PRINT *,T1
   PRINT *,T2
   PRINT *,(T1-T2)/T1

END PROGRAM
