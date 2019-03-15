PROGRAM testTransformation

   USE constantSimulation
   USE constantProgram
   USE io

   IMPLICIT NONE
   INTEGER, PARAMETER                               :: nTest=10000        !> number of step for surface density
   REAL(kind = xp), DIMENSION(nTest)                 :: randomVal, output, radius !> radius
   CHARACTER(len = 200)                             :: filename          !> Filename of simulation parameter
   INTEGER                                          :: i

   !Read filename from Bash input
   IF (COMMAND_ARGUMENT_COUNT() == 1) THEN
      CALL GET_COMMAND_ARGUMENT(1,filename)
   ELSE
      PRINT *,"No filename furnish for input file, use default file physicalsetting"
      filename = "physicalsetting"
   END IF

   !Read value frome file
   CALL readInputFile(filename)

   !initialize module
   CALL initParams()

   CALL random_number(randomVal)
   CALL random_number(radius)


   PRINT *,''
   PRINT *,'Test time'
   DO i = 1,nTest
      output(i) = timeSI(timeNormalize(randomVal(i)))
   END DO
   PRINT *,'Min',MINVAL(output/randomVal)
   PRINT *,'Max',MAXVAL(output/randomVal)
   PRINT *,'Moyenne',SUM(output/randomVal)/real(nTest,kind=xp)

   PRINT *,''
   PRINT *,'Test distance'
   output = distanceSI(distanceNormalize(randomVal))
   PRINT *,'Min',MINVAL(output/randomVal)
   PRINT *,'Max',MAXVAL(output/randomVal)
   PRINT *,'Moyenne',SUM(output/randomVal)/real(nTest,kind=xp)

   PRINT *,''
   PRINT *,'Test mass'
   DO i = 1,nTest
      output(i) = massSI(massNormalize(randomVal(i)))
   END DO
   PRINT *,'Min',MINVAL(output/randomVal)
   PRINT *,'Max',MAXVAL(output/randomVal)
   PRINT *,'Moyenne',SUM(output/randomVal)/real(nTest,kind=xp)

   PRINT *,''
   PRINT *,'Test temperature'
   output = temperatureSI(temperatureNormalize(randomVal))
   PRINT *,'Min',MINVAL(output/randomVal)
   PRINT *,'Max',MAXVAL(output/randomVal)
   PRINT *,'Moyenne',SUM(output/randomVal)/real(nTest,kind=xp)

   PRINT *,''
   PRINT *,'Test radius'
   output = radiusSI(radiusNormalize(randomVal))
   PRINT *,'Min',MINVAL(output/randomVal)
   PRINT *,'Max',MAXVAL(output/randomVal)
   PRINT *,'Moyenne',SUM(output/randomVal)/real(nTest,kind=xp)

   PRINT *,''
   PRINT *,'Test surface density'
   output = surfaceDensitySI(surfaceDensityNormalize(randomVal,radiusNormalize(radius)),radiusNormalize(radius))
   PRINT *,'Min',MINVAL(output/randomVal)
   PRINT *,'Max',MAXVAL(output/randomVal)
   PRINT *,'Moyenne',SUM(output/randomVal)/real(nTest,kind=xp)

   PRINT *,''
   PRINT *,'Test speed'
   output = speedSI(speedNormalize(randomVal))
   PRINT *,'Min',MINVAL(output/randomVal)
   PRINT *,'Max',MAXVAL(output/randomVal)
   PRINT *,'Moyenne',SUM(output/randomVal)/real(nTest,kind=xp)

   PRINT *,''
   PRINT *,'Test density'
   output = densitySI(densityNormalize(randomVal))
   PRINT *,'Min',MINVAL(output/randomVal)
   PRINT *,'Max',MAXVAL(output/randomVal)
   PRINT *,'Moyenne',SUM(output/randomVal)/real(nTest,kind=xp)

   PRINT *,''
   PRINT *,'Test pressure'
   output = pressureSI(pressureNormalize(randomVal))
   PRINT *,'Min',MINVAL(output/randomVal)
   PRINT *,'Max',MAXVAL(output/randomVal)
   PRINT *,'Moyenne',SUM(output/randomVal)/real(nTest,kind=xp)

   PRINT *,''
   PRINT *,'Test accretion rate'
   output = accretionRateSI(accretionRateNormalize(randomVal))
   PRINT *,'Min',MINVAL(output/randomVal)
   PRINT *,'Max',MAXVAL(output/randomVal)
   PRINT *,'Moyenne',SUM(output/randomVal)/real(nTest,kind=xp)

   PRINT *,''
   PRINT *,'Test massic energy'
   output = massicEnergySI(massicEnergyNormalize(randomVal))
   PRINT *,'Min',MINVAL(output/randomVal)
   PRINT *,'Max',MAXVAL(output/randomVal)
   PRINT *,'Moyenne',SUM(output/randomVal)/real(nTest,kind=xp)

   PRINT *,''
   PRINT *,'Test kinematic viscosity'
   output = kinematicViscositySI(kinematicViscosityNormalize(randomVal))
   PRINT *,'Min',MINVAL(output/randomVal)
   PRINT *,'Max',MAXVAL(output/randomVal)
   PRINT *,'Moyenne',SUM(output/randomVal)/real(nTest,kind=xp)

   PRINT *,''
   PRINT *,'Test heat capacity'
   output = heatCapacitySI(heatCapacityNormalize(randomVal))
   PRINT *,'Min',MINVAL(output/randomVal)
   PRINT *,'Max',MAXVAL(output/randomVal)
   PRINT *,'Moyenne',SUM(output/randomVal)/real(nTest,kind=xp)

   PRINT *,''
   PRINT *,'Test flux'
   output = fluxSI(fluxNormalize(randomVal))
   PRINT *,'Min',MINVAL(output/randomVal)
   PRINT *,'Max',MAXVAL(output/randomVal)
   PRINT *,'Moyenne',SUM(output/randomVal)/real(nTest,kind=xp)

END PROGRAM testTransformation
