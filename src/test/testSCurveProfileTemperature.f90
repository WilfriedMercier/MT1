PROGRAM testSCurveMesh

   USE constantSimulation
   USE constantProgram
   USE sCurve
   USE io

   IMPLICIT NONE
   REAL(kind = xp)                                  :: minTemperature, & !> minimum temperature of the grid
                                                       maxTemperature    !> maximum temperature of the grid
   INTEGER                                          :: nTemperature      !> number of step for the temperature
   REAL(kind = xp)                                  :: sigmaProfile      !>surface density
   REAL(kind = xp), DIMENSION(:, :), ALLOCATABLE :: temperature, &    !> temperature grid
                                                       sigma, &          !> surface density grid
                                                       Qp, &             !> heating grid
                                                       Qm, &             !> heat loss grid
                                                       Tau_eff,&         !> Opacity grid
                                                       xGrid
   CHARACTER(len = 200)                             :: filename, mainConfFilename          !> Filename of simulation parameter
   INTEGER                                          :: i, j, unitOut
   NAMELIST /physicalsettingSCurveProfileTemperature/ mainConfFilename,sigmaProfile,&
   minTemperature,maxTemperature,nTemperature

   !Read filename from Bash input
   IF (COMMAND_ARGUMENT_COUNT() == 1) THEN
      CALL GET_COMMAND_ARGUMENT(1,filename)
   ELSE
      PRINT *,"No filename furnish for input file, use default file physicalsettingSCurve"
      filename = "physicalsettingSCurveProfileTemperature"
   END IF

   !Load file setting
   OPEN (UNIT=10, File=trim(filename), STATUS='OLD')
   READ (10, NML=physicalsettingSCurveProfileTemperature)
   CLOSE (10)

   !Read value frome file
   CALL readInputFile(trim(mainConfFilename))

   !initialize module
   CALL initParams()
   CALL initOutput()
   CALL initConstants()

   !Allocate table
   ALLOCATE(temperature(nTemperature, nSpaceStep))
   ALLOCATE(sigma(nTemperature, nSpaceStep))
   ALLOCATE(Qp(nTemperature, nSpaceStep))
   ALLOCATE(Qm(nTemperature, nSpaceStep))
   ALLOCATE(Tau_eff(nTemperature, nSpaceStep))
   ALLOCATE(xGrid(nTemperature, nSpaceStep))

   !Initialization
   sigma = sigmaProfile
   DO j = 1, nTemperature
      temperature(j,:) = 10._xp**(log10(minTemperature)+&
      j*((log10(maxTemperature)-log10(minTemperature))/(nTemperature-1._xp)))
      xGrid(j,:) = xRadius
   END DO

   !Calculate data
   DO i = 1, nSpaceStep
      !Normalization
      temperature(:,i) = temperatureNormalize(temperature(:,i))
      sigma(:,i) = surfaceDensityNormalize(sigma(:,i), xGrid(:,i))

      !Calculation
      CALL heatingTermsFromMeshGrid(temperature(:,i), sigma(:,i),&
      Qp(:,i), Qm(:,i), xRadius(i), Tau_eff(:,i))

      !Back to SI
      temperature(:,i) = temperatureSI(temperature(:,i))
      sigma(:,i) = surfaceDensitySI(sigma(:,i), xGrid(:,i))
      Qp(:,i) = massicEnergySI(Qp(:,i))
      Qm(:,i) = massicEnergySI(Qm(:,i))
   END DO

   !Write data
   OPEN(NEWUNIT=unitOut, FILE = 'SMeshProfileTemperature.dat', STATUS='REPLACE')
   DO i = 1,nSpaceStep
      !Write mesh
      DO j = 1,nTemperature
            WRITE (unitOut, 101) i, j, rRadius(i), temperature(j,i), sigma(j,i), Qp(j,i), Qm(j,i), Tau_eff(j,i)
      END DO
   END DO
   CLOSE(unitOut)

   !Define the format of the output
   101 FORMAT (2(2x,I10),6(2x,1pg18.10e3))

   !Dealocate table
   DEALLOCATE(Tau_eff)
   DEALLOCATE(temperature)
   DEALLOCATE(sigma)
   DEALLOCATE(Qp)
   DEALLOCATE(Qm)
   DEALLOCATE(xGrid)

END PROGRAM testSCurveMesh
