MODULE simulationManager

   USE diffSolver
   USE constantSimulation
   USE constantProgram
   USE io
   !USE algebricEquations

   IMPLICIT NONE
   REAL(KIND = xp), PRIVATE                            :: time
   REAL(KIND = xp), PRIVATE, DIMENSION(:), ALLOCATABLE :: Temperature, &
                                                          Temperaturebis, &
                                                          pressureRad, &
                                                          pressureGas, &
                                                          pressureTotal, &
                                                          Sdensity, &
                                                          sigma, &
                                                          rho, &
                                                          Fz, &
                                                          accretionSpeed, &
                                                          accretionRate, &
                                                          opacityFree, &
                                                          ke, &
                                                          epsilonff, &
                                                          taueff, &
                                                          halflength, &
                                                          angularSpeed, &
                                                          soundSpeed, &
                                                          viscosity, &
                                                          Qp, &
                                                          Qm, &
                                                          Qadv, &
                                                          beta, &
                                                          heatCapacity, &
                                                          gamma3, &
                                                          stSdensity, &
                                                          stTemp
   INTEGER, PRIVATE                                    :: countOutput
   LOGICAL                                             :: superDt
   PRIVATE

   PUBLIC initSimulation, runSimulation, runExplicitSimulation, runImplicitSimulation, &
          closeSimulation, oneExplicitStep, oneImplicitStep

   CONTAINS

      !> initialises temperature and density first, and then algebraic values
      SUBROUTINE initSimulation()
        use algebricEquations
        IMPLICIT none

        PRINT *, 'Init Simulation'
        !Allocate vectors
        CALL allocation()

        !Temporary analytical solution
        CALL getInitDisc(xRadius , stTemp, sigma, variationAccretionRateInit)
        temperature = stTemp
        Sdensity = sigma * xRadius
        time = 0.0_xp
        ke = thomsonOpacityNorm
        superDt = .FALSE.
        CALL algebricStep(PressureRad, &
                                 angularSpeed, &
                                 halfLength, &
                                 soundSpeed, &
                                 viscosity, &
                                 Qp, &
                                 accretionSpeed, &
                                 accretionRate, &
                                 PressureGas, &
                                 opacityFree, &
                                 pressureTotal, &
                                 beta, &
                                 heatCapacity, &
                                 gamma3, &
                                 Qadv, &
                                 taueff, &
                                 epsilonff, &
                                 Fz, &
                                 Qm, &
                                 rho, &
                                 Sdensity, &
                                 Temperature)

        !Function defined in io (writes values in output file)
        CALL writeState(timeSI(time), &
                        radiusSI(xRadius), &
                        distanceSI(halflength), &
                        temperatureSI(Temperature), &
                        pressureSI(pressureTotal), &
                        pressureSI(pressureGas), &
                        pressureSI(pressureRad), &
                        beta, &
                        surfaceDensitySI(Sdensity, xRadius), &
                        speedSI(soundspeed), &
                        kinematicViscositySI(viscosity), &
                        speedSI(accretionSpeed), &
                        accretionRateSI(accretionRate), &
                        massicEnergySI(Qp), &
                        massicEnergySI(Qm), &
                        massicEnergySI(Qadv), &
                        heatCapacitySI(heatCapacity), &
                        fluxNormalize(Fz), &
                        opacityFree, &
                        ke, &
                        epsilonff, &
                        taueff, &
                        angularSpeed)

      END SUBROUTINE initSimulation


      !> run the simulation
      SUBROUTINE runSimulation()

         IMPLICIT NONE
         integer :: it

         PRINT *, 'Run main loop'
         countOutput = 1
         it=0

         DO WHILE (time < tmax .AND. it<itmax)  !>number of time steps during test phase
           CALL oneStep() !differential and algebraic equations

           time = time + dt
           it=it+1

           IF (countOutput == snapshotFrequency) THEN
              print*,'time = ', timeSI(time), timeSI(dt),&
              timeSI(viscT), timeSI(thermT), timeSI(dynT), timeSI(condStabilityAdv), timeSI(condStabilityDiff)
              CALL writeState(timeSI(time), &
                              radiusSI(xRadius), &
                              distanceSI(halflength), &
                              temperatureSI(Temperature), &
                              pressureSI(pressureTotal), &
                              pressureSI(pressureGas), &
                              pressureSI(pressureRad), &
                              beta, &
                              surfaceDensitySI(Sdensity, xRadius), &
                              speedSI(soundspeed), &
                              kinematicViscositySI(viscosity), &
                              speedSI(accretionSpeed), &
                              accretionRateSI(accretionRate), &
                              massicEnergySI(Qp), &
                              massicEnergySI(Qm), &
                              massicEnergySI(Qadv), &
                              heatCapacitySI(heatCapacity), &
                              fluxNormalize(Fz), &
                              opacityFree, &
                              ke, &
                              epsilonff, &
                              taueff, &
                              angularSpeed)
             countOutput = 1
          ELSE
             countOutput = countOutput + 1
          END IF
         END DO
         IF (countOutput /= 1) THEN
            print*,'time = ', timeSI(time), timeSI(dt),&
            timeSI(viscT), timeSI(thermT), timeSI(dynT), timeSI(condStabilityAdv), timeSI(condStabilityDiff)
            CALL writeState(timeSI(time), &
                            radiusSI(xRadius), &
                            distanceSI(halflength), &
                            temperatureSI(Temperature), &
                            pressureSI(pressureTotal), &
                            pressureSI(pressureGas), &
                            pressureSI(pressureRad), &
                            beta, &
                            surfaceDensitySI(Sdensity, xRadius), &
                            speedSI(soundspeed), &
                            kinematicViscositySI(viscosity), &
                            speedSI(accretionSpeed), &
                            accretionRateSI(accretionRate), &
                            massicEnergySI(Qp), &
                            massicEnergySI(Qm), &
                            massicEnergySI(Qadv), &
                            heatCapacitySI(heatCapacity), &
                            fluxNormalize(Fz), &
                            opacityFree, &
                            ke, &
                            epsilonff, &
                            taueff, &
                            angularSpeed)
         END IF
      END SUBROUTINE runSimulation


      !> run the simulation using only explicit scheme (for testing purpose only)
      SUBROUTINE runExplicitSimulation()
        IMPLICIT NONE
        integer :: it

        PRINT *, 'Run main loop'
        countOutput = 1
        it=0

        DO WHILE (time < tmax.and.it<itmax)  !>number of time steps during test phase
          CALL oneExplicitStep()
          time = time + dt
          it=it+1
          IF (countOutput == snapshotFrequency) THEN
             print*,'time = ', timeSI(time), timeSI(dt),&
             timeSI(viscT), timeSI(thermT), timeSI(dynT), timeSI(condStabilityAdv), timeSI(condStabilityDiff)
             CALL writeState(timeSI(time), &
                             radiusSI(xRadius), &
                             distanceSI(halflength), &
                             temperatureSI(Temperature), &
                             pressureSI(pressureTotal), &
                             pressureSI(pressureGas), &
                             pressureSI(pressureRad), &
                             beta, &
                             surfaceDensitySI(Sdensity, xRadius), &
                             speedSI(soundspeed), &
                             kinematicViscositySI(viscosity), &
                             speedSI(accretionSpeed), &
                             accretionRateSI(accretionRate), &
                             massicEnergySI(Qp), &
                             massicEnergySI(Qm), &
                             massicEnergySI(Qadv), &
                             heatCapacitySI(heatCapacity), &
                             fluxNormalize(Fz), &
                             opacityFree, &
                             ke, &
                             epsilonff, &
                             taueff, &
                             angularSpeed)
            countOutput = 1
         ELSE
            countOutput = countOutput + 1
         END IF
        END DO

        IF (countOutput /= 1) THEN
           print*,'time = ', timeSI(time), timeSI(dt),&
           timeSI(viscT), timeSI(thermT), timeSI(dynT), timeSI(condStabilityAdv), timeSI(condStabilityDiff)
           CALL writeState(timeSI(time), &
                           radiusSI(xRadius), &
                           distanceSI(halflength), &
                           temperatureSI(Temperature), &
                           pressureSI(pressureTotal), &
                           pressureSI(pressureGas), &
                           pressureSI(pressureRad), &
                           beta, &
                           surfaceDensitySI(Sdensity, xRadius), &
                           speedSI(soundspeed), &
                           kinematicViscositySI(viscosity), &
                           speedSI(accretionSpeed), &
                           accretionRateSI(accretionRate), &
                           massicEnergySI(Qp), &
                           massicEnergySI(Qm), &
                           massicEnergySI(Qadv), &
                           heatCapacitySI(heatCapacity), &
                           fluxNormalize(Fz), &
                           opacityFree, &
                           ke, &
                           epsilonff, &
                           taueff, &
                           angularSpeed)
        END IF
      END SUBROUTINE runExplicitSimulation


      SUBROUTINE runImplicitSimulation()
        IMPLICIT NONE
        integer :: it

        PRINT *, 'Run main loop'
        countOutput = 1
        it=0

        DO WHILE (time < tmax.and.it<itmax)  !>number of time steps during test phase
          CALL oneImplicitStep()
          time = time + dt
          it=it+1
          IF (countOutput == snapshotFrequency) THEN
             print*,'time = ', timeSI(time), timeSI(dt),&
             timeSI(viscT), timeSI(thermT), timeSI(dynT), timeSI(condStabilityAdv), timeSI(condStabilityDiff)
             CALL writeState(timeSI(time), &
                             radiusSI(xRadius), &
                             distanceSI(halflength), &
                             temperatureSI(Temperature), &
                             pressureSI(pressureTotal), &
                             pressureSI(pressureGas), &
                             pressureSI(pressureRad), &
                             beta, &
                             surfaceDensitySI(Sdensity, xRadius), &
                             speedSI(soundspeed), &
                             kinematicViscositySI(viscosity), &
                             speedSI(accretionSpeed), &
                             accretionRateSI(accretionRate), &
                             massicEnergySI(Qp), &
                             massicEnergySI(Qm), &
                             massicEnergySI(Qadv), &
                             heatCapacitySI(heatCapacity), &
                             fluxNormalize(Fz), &
                             opacityFree, &
                             ke, &
                             epsilonff, &
                             taueff, &
                             angularSpeed)
            countOutput = 1
         ELSE
            countOutput = countOutput + 1
         END IF
        END DO

        IF (countOutput /= 1) THEN
           print*,'time = ', timeSI(time), timeSI(dt),&
           timeSI(viscT), timeSI(thermT), timeSI(dynT), timeSI(condStabilityAdv), timeSI(condStabilityDiff)
           CALL writeState(timeSI(time), &
                           radiusSI(xRadius), &
                           distanceSI(halflength), &
                           temperatureSI(Temperature), &
                           pressureSI(pressureTotal), &
                           pressureSI(pressureGas), &
                           pressureSI(pressureRad), &
                           beta, &
                           surfaceDensitySI(Sdensity, xRadius), &
                           speedSI(soundspeed), &
                           kinematicViscositySI(viscosity), &
                           speedSI(accretionSpeed), &
                           accretionRateSI(accretionRate), &
                           massicEnergySI(Qp), &
                           massicEnergySI(Qm), &
                           massicEnergySI(Qadv), &
                           heatCapacitySI(heatCapacity), &
                           fluxNormalize(Fz), &
                           opacityFree, &
                           ke, &
                           epsilonff, &
                           taueff, &
                           angularSpeed)
        END IF
      END SUBROUTINE runImplicitSimulation

      SUBROUTINE oneImplicitStep()
        use algebricEquations
        IMPLICIT NONE
        !Valentin 2/12

        !> calculate the time step before calling the diffsolver code
        dt = timeStepCalculator(accretionSpeed, viscosity, xRadius, halflength, angularSpeed)

        !Verify stability condition and change dt if needed
        IF (dt > condStabilityAdv) THEN
           dt = condStabilityAdv
        END IF
        IF (dt > condStabilityDiff) THEN
           dt = condStabilityDiff
        END IF

        CALL temperatureExplicitScheme(Temperature, Qp, Qm, Qadv, heatCapacity)
        CALL densityImplicit(Sdensity,viscosity)

        !> all the algebraic equations are now computed inside 'algebricEquations.f90'
        CALL algebricStep(PressureRad, &
                                 angularSpeed, &
                                 halfLength, &
                                 soundSpeed, &
                                 viscosity, &
                                 Qp, &
                                 accretionSpeed, &
                                 accretionRate, &
                                 PressureGas, &
                                 opacityFree, &
                                 pressureTotal, &
                                 beta, &
                                 heatCapacity, &
                                 gamma3, &
                                 Qadv, &
                                 taueff, &
                                 epsilonff, &
                                 Fz, &
                                 Qm, &
                                 rho, &
                                 Sdensity, &
                                 Temperature)
      END SUBROUTINE oneImplicitStep



      !> SUBROUTINE wich execute one time step using only explicit scheme
      SUBROUTINE oneExplicitStep()
        use algebricEquations
        IMPLICIT NONE
        !Valentin 2/12

        !> calculate the time step before calling the diffsolver code
        dt = timeStepCalculator(accretionSpeed, viscosity, xRadius, halflength, angularSpeed)

        !Verify stability condition and change dt if needed
        IF (dt > condStabilityAdv) THEN
           dt = condStabilityAdv
        END IF
        IF (dt > condStabilityDiff) THEN
           dt = condStabilityDiff
        END IF

        CALL temperatureExplicitScheme(Temperature, Qp, Qm, Qadv, heatCapacity)
        CALL densityExplicit(Sdensity,viscosity)

        !> all the algebraic equations are now computed inside 'algebricEquations.f90'
        CALL algebricStep(PressureRad, &
                                 angularSpeed, &
                                 halfLength, &
                                 soundSpeed, &
                                 viscosity, &
                                 Qp, &
                                 accretionSpeed, &
                                 accretionRate, &
                                 PressureGas, &
                                 opacityFree, &
                                 pressureTotal, &
                                 beta, &
                                 heatCapacity, &
                                 gamma3, &
                                 Qadv, &
                                 taueff, &
                                 epsilonff, &
                                 Fz, &
                                 Qm, &
                                 rho, &
                                 Sdensity, &
                                 Temperature)
      END SUBROUTINE oneExplicitStep

      !> SUBROUTINE wich execute one time step
      SUBROUTINE oneStep()
        use algebricEquations
        IMPLICIT NONE
        !Valentin 16/11

        !> calculate the time step before calling the diffsolver code
        dt = timeStepCalculator(accretionSpeed, viscosity, xRadius, halflength, angularSpeed)

        !IF ((min(MinVAL(1._xp-ABS(Qp-Qm)/Qp),MinVAL(1._xp-ABS(Qp-Qm)/Qm)) > (1._xp-1e-6_xp)) .AND. (.NOT. superDt)) THEN
        IF ((max(MaxVAL(ABS(Qp-Qm)/Qp),MaxVAL(ABS(Qp-Qm)/Qm)) < 1e-4_xp) .AND. (.NOT. superDt)) THEN
           dt = viscT
           print*,'viscous time !', timeSI(time)
           CALL densitySolver(Sdensity,viscosity)
           superDt = .TRUE.
        ELSE
           !Verify stability condition and change dt if needed
           IF (dt > condStabilityAdv) THEN
              dt = condStabilityAdv
           END IF
           CALL temperatureSolver(Temperature, Qp, Qm, Qadv, heatCapacity)
           CALL densitySolver(Sdensity,viscosity)
           superDt = .FALSE.
        END IF


        !> all the algebraic equations are now computed inside 'algebricEquations.f90'
        CALL algebricStep(PressureRad, &
                                 angularSpeed, &
                                 halfLength, &
                                 soundSpeed, &
                                 viscosity, &
                                 Qp, &
                                 accretionSpeed, &
                                 accretionRate, &
                                 PressureGas, &
                                 opacityFree, &
                                 pressureTotal, &
                                 beta, &
                                 heatCapacity, &
                                 gamma3, &
                                 Qadv, &
                                 taueff, &
                                 epsilonff, &
                                 Fz, &
                                 Qm, &
                                 rho, &
                                 Sdensity, &
                                 Temperature)

      END SUBROUTINE oneStep

      !> Allocate every vector
      SUBROUTINE allocation()
        IMPLICIT NONE
        ALLOCATE(Temperature(nSpaceStep))
        ALLOCATE(Temperaturebis(nspacestep))
        ALLOCATE(pressureRad(nSpaceStep))
        ALLOCATE(pressureGas(nSpaceStep))
        ALLOCATE(pressureTotal(nSpaceStep))
        ALLOCATE(halflength(nSpaceStep))
        ALLOCATE(soundSpeed(nSpaceStep))
        ALLOCATE(viscosity(nSpaceStep))
        ALLOCATE(Qp(nSpaceStep))
        ALLOCATE(Qm(nSpaceStep))
        ALLOCATE(Qadv(nSpaceStep))
        ALLOCATE(beta(nSpaceStep))
        ALLOCATE(heatCapacity(nSpaceStep))
        ALLOCATE(Sdensity(nSpaceStep))
        ALLOCATE(Fz(nSpaceStep))
        ALLOCATE(accretionSpeed(nSpaceStep))
        ALLOCATE(accretionRate(nSpaceStep))
        ALLOCATE(opacityFree(nSpaceStep))
        ALLOCATE(ke(nSpaceStep))
        ALLOCATE(epsilonff(nSpaceStep))
        ALLOCATE(taueff(nSpaceStep))
        ALLOCATE(sigma(nspaceStep))
        ALLOCATE(stSdensity(nspaceStep))
        ALLOCATE(stTemp(nspaceStep))
        ALLOCATE(rho(nspaceStep))
        ALLOCATE(angularSpeed(nSpaceStep))
        ALLOCATE(gamma3(nSpaceStep))

      END SUBROUTINE allocation


      !> free memory used by the module
      SUBROUTINE closeSimulation()
         IMPLICIT NONE
         DEALLOCATE(Temperature)
         DEALLOCATE(Temperaturebis)
         DEALLOCATE(pressureRad)
         DEALLOCATE(pressureGas)
         DEALLOCATE(pressureTotal)
         DEALLOCATE(halflength)
         DEALLOCATE(soundSpeed)
         DEALLOCATE(viscosity)
         DEALLOCATE(Qp)
         DEALLOCATE(Qm)
         DEALLOCATE(Qadv)
         DEALLOCATE(beta)
         DEALLOCATE(heatCapacity)
         DEALLOCATE(Sdensity)
         DEALLOCATE(Fz)
         DEALLOCATE(accretionSpeed)
         DEALLOCATE(accretionRate)
         DEALLOCATE(opacityFree)
         DEALLOCATE(ke)
         DEALLOCATE(epsilonff)
         DEALLOCATE(taueff)
         DEALLOCATE(sigma)
         DEALLOCATE(stSdensity)
         DEALLOCATE(stTemp)
         DEALLOCATE(rho)
         DEALLOCATE(angularSpeed)
         DEALLOCATE(gamma3)
      END SUBROUTINE



END MODULE simulationManager
