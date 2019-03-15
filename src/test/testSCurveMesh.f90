PROGRAM testSCurveMesh

   USE constantSimulation
   USE constantProgram
   USE algebricEquations
   USE sCurve
   USE io

   IMPLICIT NONE
   REAL(kind = xp)                                  :: minTemperature, & !> minimum temperature of the grid
                                                       maxTemperature    !> maximum temperature of the grid
   INTEGER                                          :: nTemperature      !> number of step for the temperature
   REAL(kind = xp)                                  :: minSigma, &       !> minimum surface density
                                                       maxSigma          !> maximum surface density
   INTEGER                                          :: nSigma            !> number of step for surface density
   REAL(kind = xp), DIMENSION(:, :, :), ALLOCATABLE :: temperature, &    !> temperature grid
                                                       sigma, &          !> surface density grid
                                                       Qp, &             !> heating grid
                                                       Qm, &             !> heat loss grid
                                                       Tau_eff           !> Opacity grid
   REAL(kind = xp), DIMENSION(:)      , ALLOCATABLE :: highTemp, lowTemp, highSigma, lowSigma
   LOGICAL        , DIMENSION(:,:)    , ALLOCATABLE :: QTransition
   REAL(kind = xp), DIMENSION(:)      , ALLOCATABLE :: r                 !> radius
   CHARACTER(len = 200)                             :: filename, mainConfFilename          !> Filename of simulation parameter
   INTEGER                                          :: i, j, k, nS, unitOut, unitOutFrame
   NAMELIST /physicalsettingSCurve/ mainConfFilename,minSigma,maxSigma,nSigma,minTemperature,maxTemperature,nTemperature

   !Read filename from Bash input
   IF (COMMAND_ARGUMENT_COUNT() == 1) THEN
      CALL GET_COMMAND_ARGUMENT(1,filename)
   ELSE
      PRINT *,"No filename furnish for input file, use default file physicalsettingSCurve"
      filename = "physicalsettingSCurve"
   END IF

   !Load file setting
   OPEN (UNIT=10, File=trim(filename), STATUS='OLD')
   READ (10, NML=physicalsettingSCurve)
   CLOSE (10)

   !Read value frome file
   CALL readInputFile(trim(mainConfFilename))

   !initialize module
   CALL initParams()
   CALL initOutput()
   CALL initConstants()

   !Allocate table
   ALLOCATE(temperature(nSigma, nTemperature, nSpaceStep))
   ALLOCATE(sigma(nSigma, nTemperature, nSpaceStep))
   ALLOCATE(Qp(nSigma, nTemperature, nSpaceStep))
   ALLOCATE(Qm(nSigma, nTemperature, nSpaceStep))
   ALLOCATE(Tau_eff(nSigma, nTemperature, nSpaceStep))
   ALLOCATE(r(nSpaceStep))
   ALLOCATE(QTransition(nSigma-1, nTemperature-1))

   !Calculate data
   r = radiusSI(xRadius)
   !$OMP PARALLEL DO SCHEDULE(GUIDED,1)
   DO i = 1, nSpaceStep
      CALL meshGridScurve(minTemperature, maxTemperature, nTemperature, minSigma, maxSigma, nSigma, xradius(i), &
                          temperature(:,:,i), sigma(:,:,i), Qp(:,:,i), Qm(:,:,i), Tau_eff(:,:,i))
   END DO
   !$OMP END PARALLEL DO

   !Write data
   OPEN(NEWUNIT=unitOut, FILE = 'ScurveMesh.dat', STATUS='REPLACE')
   OPEN(NEWUNIT=unitOutFrame, FILE = 'ScurveFrame.dat', STATUS='REPLACE')
   DO i = 1,nSpaceStep
      !Write mesh
      DO k = 1,nTemperature
         DO j = 1,nSigma
            WRITE (unitOut, 101) i, j, k, r(i), temperature(j,k,i), sigma(j,k,i), Qp(j,k,i), Qm(j,k,i), Tau_eff(j,k,i)
         END DO
      END DO

      !Calculate and write frame
      CALL determineMeshTransition(sigma(:,:,i), temperature(:,:,i), nSigma, nTemperature, Qp(:,:,i), Qm(:,:,i),&
      QTransition,nS)
      ALLOCATE(highTemp(nS))
      ALLOCATE(lowTemp(nS))
      ALLOCATE(highSigma(nS))
      ALLOCATE(lowSigma(nS))
      CALL calculateSFrame(sigma(:,:,i), temperature(:,:,i), nSigma, nTemperature, Qp(:,:,i), Qm(:,:,i),&
      QTransition,lowSigma,highSigma,lowTemp,highTemp)
      DO j = 1,nS
         WRITE (unitOutFrame, 102) i, j, lowSigma(j), highSigma(j), lowTemp(j), highTemp(j)
      END DO
      DEALLOCATE(lowSigma)
      DEALLOCATE(highSigma)
      DEALLOCATE(lowTemp)
      DEALLOCATE(highTemp)
   END DO
   CLOSE(unitOut)
   CLOSE(unitOutFrame)

   !Define the format of the output
   101 FORMAT (3(2x,I10),6(2x,1pg18.10e3))
   102 FORMAT (2(2x,I10),4(2x,1pg18.10e3))

   !Dealocate table
   DEALLOCATE(QTransition)
   DEALLOCATE(Tau_eff)
   DEALLOCATE(temperature)
   DEALLOCATE(sigma)
   DEALLOCATE(Qp)
   DEALLOCATE(Qm)
   DEALLOCATE(r)

END PROGRAM testSCurveMesh
