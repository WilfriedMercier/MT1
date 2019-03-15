PROGRAM testSimulationExplicitScheme

   USE constantProgram
   USE constantSimulation
   USE simulationManager
   USE algebricEquations
   USE io

   IMPLICIT NONE
   CHARACTER(len = 200) :: filename

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
   CALL initOutput()

   CALL initParams()

   CALL initConstants()

   CALL initSimulation()

   !Run simulation
   CALL runExplicitSimulation()

   !Close simulation
   CALL CloseSimulation()
   CALL closeOutput()

END PROGRAM testSimulationExplicitScheme
