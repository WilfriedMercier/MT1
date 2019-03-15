PROGRAM testInitialCondition

   USE constantSimulation
   USE algebricEquations
   USE sCurve
   USE io

   IMPLICIT NONE
   REAL(kind = xp), DIMENSION(:), ALLOCATABLE       :: temperature, &    !> temperature grid
                                                       sigma, &          !> surface density grid
                                                       Qp, &             !> heating grid
                                                       Qm, &             !> heat loss grid
                                                       S, &
                                                       f
   REAL(kind = xp), DIMENSION(:)      , ALLOCATABLE :: r                 !> radius
   CHARACTER(len = 200)                             :: filename          !> Filename of simulation parameter
   INTEGER                                          :: i, unitOut

   !Read filename from Bash input
   IF (COMMAND_ARGUMENT_COUNT() == 1) THEN
      CALL GET_COMMAND_ARGUMENT(1,filename)
   ELSE
      PRINT *,"No filename furnish for input file, use default file physicalsetting"
      filename = "physicalsetting"
   END IF

   !Read value frome file
   CALL readInputFile(trim(filename))

   !initialize module
   CALL initParams()
   CALL initOutput()
   CALL initConstants()

   !Allocate table
   ALLOCATE(temperature(nSpaceStep))
   ALLOCATE(sigma(nSpaceStep))
   ALLOCATE(Qp(nSpaceStep))
   ALLOCATE(Qm(nSpaceStep))
   ALLOCATE(r(nSpaceStep))
   ALLOCATE(S(nSpaceStep))

   !Calculate data
   CALL getInitDisc(xRadius,temperature,sigma,variationAccretionRateInit)
   f = getFunctionF(xRadius)

   !Write data
   OPEN(NEWUNIT=unitOut, FILE='discInitValue.dat', STATUS='REPLACE')
   DO i = 1,nSpaceStep
      WRITE (unitOut, 101) i, rRadius(i), temperature(i), sigma(i), f(i)
   END DO
   CLOSE(unitOut)

   !Define the format of the output
   101 FORMAT (1(2x,I10),4(2x,1pg18.10e3))

   !Dealocate table
   DEALLOCATE(temperature)
   DEALLOCATE(sigma)
   DEALLOCATE(Qp)
   DEALLOCATE(Qm)
   DEALLOCATE(r)
   DEALLOCATE(S)

END PROGRAM testInitialCondition
