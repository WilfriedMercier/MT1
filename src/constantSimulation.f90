MODULE constantSimulation
   USE constantPhysic
   USE constantProgram

   IMPLICIT NONE

   REAL(kind = xp), PUBLIC                            :: omegaMax, &          !> Maximum rotation speed in SI units
                                                         rs, &                !> Schwartzild radius in SI unit
                                                         T0, &                !> Typical temperature in SI unit
                                                         sigma0, &
                                                         mass0, &
                                                         rMin, &              !> Minimum radius in SI unit
                                                         accretionRate0,&     !> Initial accrestion rate in SI unit
                                                         blackHoleMassSI, &   !> Mass of the black hole in SI unit
                                                         blackHoleMass, &     !> Mass of the black hole in simulation unit
                                                         tmaxSI, &            !> Time at which stopping the simulation in SI unit
                                                         tmax, &              !> Time at which stopping the simulation in Simulation unit
                                                         itmax, &
                                                         rmax, &              !> maximal radius in SI unit
                                                         dx, &                !> Radius normalized step
                                                         xMin, &              !> Min radius normalized in simulation unit
                                                         xMax, &              !> Max radius normalized in simulation unit
                                                         GNormConst, &        !> Gravitational constant normalized
                                                         cNormConst, &        !> Speed of light normalized
                                                         kbNormConst, &       !> Stephan Boltzman constant normalized
                                                         hNormConst, &        !> Planck constant normalized
                                                         sigmaNormConst, &    !> Stephan constant normalized
                                                         mpNormConst, &       !> Mass proton normalized
                                                         thomsonOpacityNorm, &
                                                         aNormConst, &
                                                         RNormConst, &
                                                         opacityFreeNormConst, &
                                                         emissivityFreeNormConst, &
                                                         atomicMass           !> Mean atomic mass
   REAL(kind = xp), PUBLIC, DIMENSION(:), ALLOCATABLE :: xRadius, rRadius     !> List of the normalized radius position
   INTEGER        , PUBLIC                            :: nSpaceStep           !> Number of space step
   REAL(kind = xp), PRIVATE                           :: LTot                 !> Total Luminosity in SI unit

   REAL(kind = xp), PUBLIC                            :: chemX                !> Chemical composition
   REAL(kind = xp), PUBLIC                            :: chemY                !> Chemical composition
   REAL(kind = xp), PUBLIC                            :: chemZ                !> Chemical composition
   REAL(kind = xp), PUBLIC                            :: alpha                !> Alpha parameter

   INTEGER        , PUBLIC                            :: snapshotFrequency
   REAL(kind = xp), PUBLIC                            :: variationAccretionRateInit,&  !> parameter for init disc with a not full disc
                                                         minScaleTemperatureMesh, &    !>for mesh calculation
                                                         maxScaleTemperatureMesh, &    !>for mesh calculation
                                                         minScaleSigmaMesh, &          !>for mesh calculation
                                                         maxScaleSigmaMesh             !>for mesh calculation
   INTEGER        , PUBLIC                            :: nTemperatureMesh, nSigmaMesh  !>for mesh calculation
   LOGICAL        , PUBLIC                            :: outputMesh


   CONTAINS


!Routine for initializing value

      SUBROUTINE initParams()

        IMPLICIT NONE

        CALL calculateNormalizationFactors() !Program units computation
        CALL calculateParameter()            !Compute the atomic mass
        CALL normalizedConstants()           !Normalize all the constants
        CALL normalizedParameter()           !Normalize the black hole mass as well as the maximum time (found in ?)
        CALL calculateSimulationStep()       !Compute the spatial step as well as xRadius

      END SUBROUTINE initParams

      SUBROUTINE closeParams()

         IMPLICIT NONE

         DEALLOCATE(xRadius)
         DEALLOCATE(rRadius)
      END SUBROUTINE closeParams

      !> Subroutine which initializes main parameters with mass and initial accretion rate
      SUBROUTINE calculateNormalizationFactors()

         USE constantPhysic

         IMPLICIT NONE

         rs         = 2.0_xp*Gconst*blackHoleMassSI/(cConst*cConst)                       !Schwarzschild's radius
         rMin       = 3.0_xp*rs                                                           !Minimum radius of study
         rMax       = rMax*rs
         omegaMax   = sqrt(Gconst*blackHoleMassSI/(rMin**3._xp))                           !Maximum angular veclocity
         LTot       = accretionRate0*cConst*cConst/12.0_xp                                !Total disc's luminosity
         T0         = (LTot/(36.0_xp*piConst*sigmaConst*(rs**2._xp)))**(0.25_xp)           !Normalization temperature
         Mass0      = accretionRate0/omegaMax
         Sigma0     = Mass0/(rs*rs)
      END SUBROUTINE calculateNormalizationFactors


     !>Subroutine which takes physical constants and normalizes them
     !>denoted as 'constant_norm'. The '_norm' means 'normalized'.
     SUBROUTINE normalizedConstants()

        USE constantPhysic

        IMPLICIT NONE

        GNormConst     = Gconst * accretionRate0 / ((omegaMax*rs)**3._xp)
        cNormConst     = cConst / (rs*omegaMax)
        kbNormConst    = kbConst * T0 / (omegaMax * accretionRate0 * rs * rs)
        hNormConst     = hConst / (accretionRate0 * rs * rs)
        sigmaNormConst = sigmaConst * (T0**4._xp) / (omegaMax * omegaMax * accretionRate0)
        mpNormConst    = massNormalize(mpConst)
        thomsonOpacityNorm = thomsonOpacityConst * (1._xp+chemX) * accretionRate0 / (rs*rs*omegaMax)
        aNormConst = aConst * (T0**4._xp) * rs / (omegaMax * accretionRate0)
        RNormConst = kbNormConst / mpNormConst

        opacityFreeNormConst = opacityFreeConst * accretionRate0 * accretionRate0 /&
        (omegaMax * omegaMax * (rs**5._xp) * (T0**3.5_xp))

        emissivityFreeNormConst = emissivityFreeConst * accretionRate0 * (T0**0.5_xp) /&
        (rs**5._xp * omegaMax**4._xp)

     END SUBROUTINE normalizedConstants


     !>Subroutine which normalizes parameters
     SUBROUTINE normalizedParameter()

        IMPLICIT NONE

        blackHoleMass  = massNormalize(blackHoleMassSI)
        tmax = timeNormalize(tmaxSI)

     END SUBROUTINE normalizedParameter


     !> Subroutine which computes program parameters
     SUBROUTINE calculateParameter()

        IMPLICIT NONE

         IF (chemX*chemY*chemZ==0.0) THEN
            WRITE(*,*)"Error : one of the abudances is 0"
            RETURN
         ENDIF

        atomicMass = 1.0_xp/(2.0_xp*chemX+0.75_xp*chemY+0.5_xp*chemZ)

     END SUBROUTINE calculateParameter


     !>Subroutine which calculates spaceStep and builds the normalized distance vector
     SUBROUTINE calculateSimulationStep()

        IMPLICIT NONE
        REAL(kind = xp), DIMENSION(1) :: tempVar
        INTEGER                       :: i

        tempVar = radiusNormalize((/rMax/))
        xMax = tempVar(1)
        tempVar = radiusNormalize((/rMin/))
        xMin = tempVar(1)
        dx = (xMax-xMin)/REAL(nSpaceStep,xp)
        ALLOCATE(xRadius(nSpaceStep))
        ALLOCATE(rRadius(nSpaceStep))

        DO i=1,nSpaceStep
           xRadius(i) = xMin + dx*REAL(i,xp)
        END DO
        rRadius = radiusSI(xRadius)

     END SUBROUTINE calculateSimulationStep


!----------------------------------------------------------------------------------------
!Function to normalize value in simulation units

      !> Function which normalizes time in simulation unit
      FUNCTION timeNormalize(t)

         IMPLICIT NONE
         REAL(kind = xp), INTENT(in) :: t
         REAL(kind = xp)             :: timeNormalize

         timeNormalize = t*omegaMax

      END FUNCTION timeNormalize


      !> Function which normalizes distance
      FUNCTION distanceNormalize(dist)

         IMPLICIT NONE
         REAL(kind = xp), DIMENSION(:)         , INTENT(in) :: dist
         REAL(kind = xp), DIMENSION(size(dist))             :: distanceNormalize

         distanceNormalize = dist/rs

      END FUNCTION distanceNormalize


      !> Function which normalizes mass in simulation unit
      FUNCTION massNormalize(mass)

         IMPLICIT NONE
         REAL(kind = xp), INTENT(in) :: mass
         REAL(kind = xp)             :: massNormalize

         massNormalize = mass*omegaMax/accretionRate0

      END FUNCTION massNormalize


      !> Function which normalizes temperature in simulation unit
      FUNCTION temperatureNormalize(temperature)

         IMPLICIT NONE
         REAL(kind = xp), DIMENSION(:)                , INTENT(in) :: temperature
         REAL(kind = xp), DIMENSION(size(temperature))             :: temperatureNormalize

         temperatureNormalize = temperature/T0

      END FUNCTION temperatureNormalize


      !> Function which normalizes the radius
      FUNCTION radiusNormalize(r)

         IMPLICIT NONE
         REAL(kind = xp), DIMENSION(:)       , INTENT(in) :: r
         REAL(kind = xp), DIMENSION(size(r))              :: radiusNormalize

         radiusNormalize = sqrt(r/rs)

      END FUNCTION radiusNormalize


      !> Normalizes surface density
      FUNCTION surfaceDensityNormalize(sigma, x)

         IMPLICIT NONE
         REAL(kind = xp), DIMENSION(:)           , INTENT(in) :: sigma
         REAL(kind = xp), DIMENSION(:)           , INTENT(in) :: x
         REAL(kind = xp), DIMENSION(size(sigma))              :: surfaceDensityNormalize

         surfaceDensityNormalize = sigma * x / sigma0

      END FUNCTION surfaceDensityNormalize


      !> Function which normalizes speed
      FUNCTION speedNormalize(v)

         IMPLICIT NONE
         REAL(kind = xp), DIMENSION(:)       , INTENT(in) :: v
         REAL(kind = xp), DIMENSION(size(v))              :: speedNormalize

         speedNormalize = v/(rs*omegaMax)

      END FUNCTION speedNormalize


      !> Function which normalizes 3D density
      FUNCTION densityNormalize(rho)

         IMPLICIT NONE
         REAL(kind = xp), DIMENSION(:)         , INTENT(in) :: rho
         REAL(kind = xp), DIMENSION(size(rho))              :: densityNormalize

         densityNormalize = rho*(rs**3.0_xp)*omegaMax/accretionRate0

      END FUNCTION densityNormalize


      !> Function which normalizes pressure
      FUNCTION pressureNormalize(p)

         IMPLICIT NONE
         REAL(kind = xp), DIMENSION(:)       , INTENT(in) :: p
         REAL(kind = xp), DIMENSION(size(p))              :: pressureNormalize

         pressureNormalize = p*rs/(omegaMax*accretionRate0)

      END FUNCTION pressureNormalize


      !> Function which normalizes the accretion rate
      FUNCTION accretionRateNormalize(accretionRate)

         IMPLICIT NONE
         REAL(kind = xp), DIMENSION(:)                  , INTENT(in) :: accretionRate
         REAL(kind = xp), DIMENSION(size(accretionRate))             :: accretionRateNormalize

         accretionRateNormalize = accretionRate/accretionRate0

      END FUNCTION accretionRateNormalize


      !> Function which normalizes energy density
      FUNCTION massicEnergyNormalize(massicEnergy)

         IMPLICIT NONE
         REAL(kind = xp), DIMENSION(:)                 , INTENT(in) :: massicEnergy
         REAL(kind = xp), DIMENSION(size(massicEnergy))             :: massicEnergyNormalize

         massicEnergyNormalize = massicEnergy/((omegaMax**3._xp)*rs*rs)

      END FUNCTION massicEnergyNormalize


      !> Function which normalizes kinematic viscosity
      FUNCTION kinematicViscosityNormalize(kinematicViscosity)

         IMPLICIT NONE
         REAL(kind = xp), DIMENSION(:)                       , INTENT(in) :: kinematicViscosity
         REAL(kind = xp), DIMENSION(size(kinematicViscosity))             :: kinematicViscosityNormalize

         kinematicViscosityNormalize = 3._xp *kinematicViscosity / (4._xp * rs * rs * omegaMax)

      END FUNCTION kinematicViscosityNormalize


      !> function which normalizes heatCapacity
      FUNCTION heatCapacityNormalize(heatCapacity)

         IMPLICIT NONE
         REAL(kind = xp), DIMENSION(:)                 , INTENT(in) :: heatCapacity
         REAL(kind = xp), DIMENSION(size(heatCapacity))             :: heatCapacityNormalize

!         heatCapacityNormalize = heatCapacity * T0 / (omegaMax * rs * rs * accretionRate0)
         heatCapacityNormalize = heatCapacity * T0 / (omegaMax * rs * rs * omegaMax)

      END FUNCTION heatCapacityNormalize


      !> function which normalizes the flux
      FUNCTION fluxNormalize(flux)

         IMPLICIT NONE
         REAL(kind = xp), DIMENSION(:)         , INTENT(in) :: flux
         REAL(kind = xp), DIMENSION(size(flux))             :: fluxNormalize

         fluxNormalize = flux / (omegaMax * omegaMax * accretionRate0)
      END FUNCTION fluxNormalize


!-------------------------------------------------------------------------------------
!Functions used to go back to the International System


      !> Function which normalizes time in simulation unit
      FUNCTION timeSI(t)

         IMPLICIT NONE
         REAL(kind = xp), INTENT(in) :: t
         REAL(kind = xp)             :: timeSI

         timeSI = t/omegaMax

      END FUNCTION timeSI


      !> Function to normalize distance
      FUNCTION distanceSI(dist)

         IMPLICIT NONE
         REAL(kind = xp), DIMENSION(:)         , INTENT(in) :: dist
         REAL(kind = xp), DIMENSION(size(dist))             :: distanceSI

         distanceSI = dist*rs

      END FUNCTION distanceSI


      !> Function to normalize mass in simulation unit
      FUNCTION massSI(mass)

         IMPLICIT NONE
         REAL(kind = xp), INTENT(in) :: mass
         REAL(kind = xp)             :: massSI

         massSI = mass*accretionRate0/omegaMax

      END FUNCTION massSI


      !> Function for normalize temperature in simulation unit
      FUNCTION temperatureSI(temperature)

         IMPLICIT NONE
         REAL(kind = xp), DIMENSION(:)                , INTENT(in) :: temperature
         REAL(kind = xp), DIMENSION(size(temperature))             :: temperatureSI

         temperatureSI = temperature*T0

      END FUNCTION temperatureSI


      !> Function to normalize radius
      FUNCTION radiusSI(x)

         IMPLICIT NONE
         REAL(kind = xp), DIMENSION(:)       , INTENT(in) :: x
         REAL(kind = xp), DIMENSION(size(x))              :: radiusSI

         radiusSI = rs * x * x

      END FUNCTION radiusSI


      !> Takes S to return sigma value in SI
      FUNCTION surfaceDensitySI(S, x)

         IMPLICIT NONE
         REAL(kind = xp), DIMENSION(:)      , INTENT(in) :: S
         REAL(kind = xp), DIMENSION(:)      , INTENT(in) :: x
         REAL(kind = xp), DIMENSION(size(S))             :: surfaceDensitySI

         surfaceDensitySI = S * sigma0 / x

      END FUNCTION surfaceDensitySI


      !> Function to normalize speed
      FUNCTION speedSI(v)

         IMPLICIT NONE
         REAL(kind = xp), DIMENSION(:)       , INTENT(in) :: v
         REAL(kind = xp), DIMENSION(size(v))              :: speedSI

         speedSI = v * rs * omegaMax

      END FUNCTION speedSI


      !> Function to normalize density
      FUNCTION densitySI(rho)

         IMPLICIT NONE
         REAL(kind = xp), DIMENSION(:)         , INTENT(in) :: rho
         REAL(kind = xp), DIMENSION(size(rho))              :: densitySI

         densitySI = rho * accretionRate0  / ((rs**3.0_xp)*omegaMax)

      END FUNCTION densitySI


      !> Function to normalize pressure
      FUNCTION pressureSI(p)

         IMPLICIT NONE
         REAL(kind = xp), DIMENSION(:)       , INTENT(in) :: p
         REAL(kind = xp), DIMENSION(size(p))              :: pressureSI

         pressureSI = p * omegaMax * accretionRate0 / rs

      END FUNCTION pressureSI


      !> Function to normalize accretion rate
      FUNCTION accretionRateSI(accretionRate)

         IMPLICIT NONE
         REAL(kind = xp), DIMENSION(:)                  , INTENT(in) :: accretionRate
         REAL(kind = xp), DIMENSION(size(accretionRate))             :: accretionRateSI

         accretionRateSI = accretionRate*accretionRate0

      END FUNCTION accretionRateSI


      !> Function to normalize energy density
      FUNCTION massicEnergySI(massicEnergy)

         IMPLICIT NONE
         REAL(kind = xp), DIMENSION(:)                 , INTENT(in) :: massicEnergy
         REAL(kind = xp), DIMENSION(size(massicEnergy))             :: massicEnergySI

         massicEnergySI = massicEnergy * (omegaMax**3._xp) * rs * rs

      END FUNCTION massicEnergySI


      !> Function to normalize kinematic viscosity
      FUNCTION kinematicViscositySI(kinematicViscosity)

         IMPLICIT NONE
         REAL(kind = xp), DIMENSION(:)                       , INTENT(in) :: kinematicViscosity
         REAL(kind = xp), DIMENSION(size(kinematicViscosity))             :: kinematicViscositySI

         kinematicViscositySI = kinematicViscosity * rs * rs * omegaMax

      END FUNCTION kinematicViscositySI


      !> function to normalize heatCapacity
      FUNCTION heatCapacitySI(heatCapacity)

         IMPLICIT NONE
         REAL(kind = xp), DIMENSION(:)                 , INTENT(in) :: heatCapacity
         REAL(kind = xp), DIMENSION(size(heatCapacity))             :: heatCapacitySI

         heatCapacitySI = heatCapacity * omegaMax * rs * rs * omegaMax / T0

      END FUNCTION heatCapacitySI


      !> function to normalize the flux
      FUNCTION fluxSI(flux)

         IMPLICIT NONE
         REAL(kind = xp), DIMENSION(:)         , INTENT(in) :: flux
         REAL(kind = xp), DIMENSION(size(flux))             :: fluxSI

         fluxSI = flux * omegaMax * omegaMax * accretionRate0
      END FUNCTION fluxSI

END MODULE constantSimulation
