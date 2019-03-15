PROGRAM trouNoir

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
      PRINT *,"No filename given as input file. Using default file : physicalsetting"
      filename = "physicalsetting"
   END IF


   !Read value frome file (function in io)
   CALL readInputFile(filename)

   !Initializes output file and directory (function in io)
   CALL initOutput()

   !Initializes all parameters and constants in program units (function in constantSimulation)
   CALL initParams()

   !Initializes alpha (function in algebricEquations)
   CALL initConstants()

   !Initializes vectors values such as temperature (function in simulationManager)
   CALL initSimulation()

   !Run simulation
   CALL runSimulation()

   !Close simulation and ouput file
   CALL closeSimulation()
   CALL closeOutput()

END PROGRAM trouNoir
