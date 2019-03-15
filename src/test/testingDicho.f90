PROGRAM testingDicho
   !MERCIER Wilfried - 02/09/18 - 09/11/19

   USE constantProgram
   USE io
   USE simulationManager
   USE constantSimulation
   USE algebricEquations
   USE sCurve

   IMPLICIT NONE

   REAL(KIND = xp)      :: X0, Y0         !Initial position
   REAL(KIND = xp)      :: DeltaS, Delta  !Steps used by the dichotomy
   REAL(KIND = xp)      :: radius         !Radius at which we look for the sCurve
   INTEGER              :: flag           !Just an error flag
   CHARACTER(len = 300) :: filename

   !Temporary
   INTEGER                       :: i, j, unitout2, NMAX, iNStep, jNStep
   REAL(KIND = xp)               :: iStep, jStep, iMax, jMax,  iMin, jMin, accuracy
   REAL(KIND = xp), DIMENSION(1) :: Temp, Surf, DQ, RAD, TauEff, func, X0Mat, Y0Mat
   !------------------------------------------------------------------------

   ! ...Test1...
   IF (1==2) THEN
!      CALL Dichotomy(DeltaOfDims, -9.3_xp, -22.0_xp, X0, Y0, 0.01_xp, 1, 10.0_xp, flag)
      WRITE(*,*)flag
   ENDIF

   ! ...Test2...
   IF (1==2) THEN
      X0 = -10.0_xp
      Y0 = -250.0_xp
      DeltaS = 1.0_xp
      accuracy = 1.0_xp
      CALL sCurveFromDichotomy(CircleLand, X0, Y0, 10000.0_xp, -10000.0_xp, 10000.0_xp, -10000.0_xp, &
           accuracy, 0.3_xp, DeltaS, 0, 100, 1000000, 1.0_xp, "testi")
   ENDIF

   ! ...Test3...
   filename = "physicalsetting"
   CALL readInputFile(filename)
   CALL initOutput()
   CALL initParams()
   CALL initConstants()

   !Maximum number of steps before stopping the search of the first crossing
   NMAX = 1000000
   DeltaS   = 1.05_xp
   Delta    = 1.0_xp
   accuracy = 1.0E-13_xp
   radius   = 35
   RAD(1)   = radius

   Y0Mat(1) = 1.0E6_xp
   Y0Mat    = temperatureNormalize(Y0Mat)
   Y0       = Y0Mat(1)

   X0Mat(1) = 1.0E3_xp
   X0Mat    = SurfaceDensityNormalize(X0Mat, RAD)
   X0       = X0Mat(1)

!   !Print path name
   i = GETCWD(filename)
   WRITE(*,*)"Current path: ", filename

!   !Write phase space in file
   IF (1==2) THEN
      iMax     = 100000.0_xp
      iMin     = 3.0_xp
      iNStep   = 1000
      iStep    = (LOG10(iMax) - LOG10(iMin))/iNStep

      jMax     = 50.0_xp
      jMin     = 0.3_xp
      jNStep   = 1000
      jStep    = (LOG10(jMax) - LOG10(jMin))/jNStep
!
!      X0Mat(1) = 1.0E4_xp
!      X0Mat    = SurfaceDensityNormalize(X0Mat, RAD)
!      iMin     = X0Mat(1)
!      X0Mat(1) = 3.0E4_xp
!      X0Mat    = SurfaceDensityNormalize(X0Mat, RAD)
!      iMax     = X0Mat(1)
!      iNStep   = 2000
!      iStep    = (LOG10(iMax) - LOG10(iMin))/iNStep
!
!      Y0Mat(1) = 2.7E7_xp
!      Y0Mat    = temperatureNormalize(Y0Mat)
!      jMax     = Y0Mat(1)
!      Y0Mat(1) = 1.50E7_xp
!      Y0Mat    = temperatureNormalize(Y0Mat)
!      jMin     = Y0Mat(1)
!      jNStep   = 2000
!      jStep    = (LOG10(jMax) - LOG10(jMin))/jNStep
!
      OPEN(NEWUNIT=unitOut2, FILE="dichospace35", STATUS='REPLACE')
      WRITE(unitOut2, '(I15, I15)')iNStep, jNStep
      DO i=0, iNStep, 1
         DO j=0, jNStep, 1
            Temp(1)  = jMin*10**(j*jStep)
            Surf(1)  = iMin*10**(i*iStep)
            DQ(1) = ComputeDelQ(Surf(1), Temp(1), TauEff(1), RAD(1))

            !Convert in SI and write in file
            WRITE(unitOut2, '(E15.5, E15.5, E15.5, E15.5)')surfaceDensitySI(Surf, RAD), temperatureSI(Temp),&
                                                         massicEnergySI(DQ), TauEff
         ENDDO
      ENDDO
      CLOSE(unitOut2)
   ENDIF

   !Temporary
   IF (1==2) THEN
      X0Mat(1) = X0
      Y0Mat(1) = Y0
      CALL get_function_f(RAD, func)
      CALL get_init_sigma(RAD, func, Y0Mat)
      CALL get_init_temperature(RAD, func, X0Mat)

      WRITE(*,*)"f = ", func
      WRITE(*,*)"X = ", RAD
      WRITE(*,*)"T = ", Y0Mat
      WRITE(*,*)"S = ", X0Mat
   ENDIF

   filename = "dichodataGoodBranch"
!   CALL sCurveFromDichotomy(ComputeDelQ, X0, Y0, 1.0E10_xp, 1.0E10_xp, 0.0_xp, 0.0_xp,&
!                            accuracy, Delta, DeltaS, 1, 200*NMAX, NMAX, radius, filename)

!   CALL WriteMaxSlopePoint(accuracy, Delta, DeltaS, 1, 200*NMAX)

   CALL BuildAllSCurves(accuracy, Delta, DeltaS, 1, 200*NMAX, 200*NMAx)


   CONTAINS


!>Test functions for the dichotomy (can be anything)
!------------------------------------------------------------------------
   !>Straight line in X-Y plane
   REAL(KIND = xp) FUNCTION DeltaOfDims(X ,Y)

      REAL(KIND = xp), INTENT(IN) :: X, Y

      DeltaOfDims = X-Y
   END FUNCTION DeltaOfDims


!-------------------------------------------------------------------------
   !>Difference of squared coordinate components
   REAL(KIND = xp) FUNCTION DeltaOfSquares(X ,Y)

      REAL(KIND = xp), INTENT(IN) :: X, Y

      IF (X<0) THEN
         DeltaOfSquares = X**2.0_xp-Y**2.0_xp
      ELSE
         DeltaOfSquares = 10.0_xp
      ENDIF
   END FUNCTION DeltaOfSquares


!-------------------------------------------------------------------------
   !>Pretty nice circle
   REAL(KIND = xp) FUNCTION CircleLand(X ,Y, T, R)

      REAL(KIND = xp), INTENT(IN) :: X, Y, T, R

      CircleLand = X**2.0_xp+Y**2.0_xp - 30000.0_xp
   END FUNCTION CircleLand

END PROGRAM testingDicho
