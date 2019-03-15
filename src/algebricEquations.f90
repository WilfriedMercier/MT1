MODULE algebricEquations

   USE constantPhysic
   USE constantProgram
   USE constantSimulation

! les noms des constantes dans ce fichier font references aux constantes normalisees


CONTAINS


!--------------------------------------------------------------------------------------------
   !>initialisation constante algebric equation
   SUBROUTINE initConstants()

      USE constantPhysic

      IMPLICIT NONE

      alpha = 1.0_xp

  END SUBROUTINE


  !> Equation  23, Page 15
  FUNCTION getFunctionF(x)
     REAL(kind = xp), INTENT(IN)   , DIMENSION(:)       :: x
     REAL(kind = xp)               , DIMENSION(size(x)) :: getFunctionF

     !----------------------------------

     getFunctionF=1.0_xp-sqrt(3.0_xp)/x
  END FUNCTION getFunctionF

  !> Equation 21, page 15
  FUNCTION getInitTemperature(x,f,accretionRateScale)
     REAL(kind = xp), INTENT(IN)   , DIMENSION(:)       :: x,f
     REAL(kind = xp), INTENT(IN)                        :: accretionRateScale
     REAL(kind = xp)               , DIMENSION(size(x)) :: getInitTemperature

     !----------------------------------

     getInitTemperature =  (1.4e4_xp/T0) * (alpha**(-0.2_xp)) * &
     ((accretionRateScale*accretionRate0/ 1.0e13_xp)**0.3_xp) * & !To change 1e-4
     ((blackHoleMassSI /  MSunSI)**(0.25_xp)) * &
     (((x**2.0_xp) * rs /  1.0e8_xp)**(-0.75_xp)) &
     * (f**0.3_xp)

  END FUNCTION getInitTemperature

  !> Equation 22, page 15
  FUNCTION getInitSigma(x,f,accretionRateScale)
     REAL(kind = xp), INTENT(IN)   , DIMENSION(:)       :: x,f
     REAL(kind = xp), INTENT(IN)                        :: accretionRateScale
     REAL(kind = xp)               , DIMENSION(size(x)) :: getInitSigma

     !----------------------------------

     getInitSigma =  (52._xp / sigma0) * (alpha**(-0.8_xp)) * &
     ((accretionRateScale*accretionrate0 / 1.0e13_xp)**0.7_xp) * & !To change 1e-4
     ((blackHoleMassSI /  MSunSI)**0.25_xp) * &
     (((x**2.0_xp) * rs /  1.0e8_xp)**(-0.75_xp)) &
     * (f**0.7_xp)
  END FUNCTION getInitSigma

  !> Initial disc state
  SUBROUTINE getInitDisc(x, T, sigma, accretionRateScale)
     REAL(kind = xp), INTENT(IN)    :: x(:), accretionRateScale
     REAL(kind = xp), INTENT(INOUT) :: T(:), sigma(:)
     ! .. local ..
     REAL(kind = xp)                :: f(size(T))
     !----------------------
     f = getFunctionF(x)
     T = getInitTemperature(x,f,accretionRateScale)
     sigma = getInitSigma(x,f,accretionRateScale)
  END SUBROUTINE getInitDisc

!--------------------------------------------------------------------------------------------
   !> calcule la vitesse angulaire Eq (2.3.1)
   FUNCTION getAngularSpeed(xradius)

      IMPLICIT NONE

      REAL(kind=xp), DIMENSION(:), INTENT(in) :: xradius
      REAL(kind=xp), DIMENSION(size(xradius))  :: getAngularSpeed

      getAngularSpeed = (sqrt(3.0_xp)/xradius)**3._xp

      IF (ANY(getAngularSpeed <= 0.0_xp)) THEN
         PRINT *, "negative angular speed"
      END IF

    END FUNCTION


!---------------------------------------------------------------------------------------------
   !> calcule la pression de radiation Eq(2.3.3)
   FUNCTION getPressureRad(Temp)

      IMPLICIT NONE

      REAL(kind=xp), DIMENSION(:)         , INTENT(in) :: Temp
      REAL(kind=xp), DIMENSION(size(Temp))             :: getPressureRad

      getPressureRad = (1.0_xp/3.0_xp) * aNormConst * Temp**4._xp

      IF (ANY(getPressureRad <= 0.0_xp)) THEN
         PRINT *, "negative radiation pressure"
      END IF

   END FUNCTION


!---------------------------------------------------------------------------------------------
   !> calcule la pression de radiation Eq(2.3.3)
   FUNCTION getPressureGas(rho, Temp)

      IMPLICIT NONE

      ! REAL(kind=xp), INTENT(in)                 :: atomicMass
      REAL(kind=xp), DIMENSION(:)         , INTENT(in)   :: Temp
      REAL(kind=xp), DIMENSION(:)         , INTENT(in)   :: rho
      REAL(kind=xp), DIMENSION(size(Temp))               :: getPressureGas


      getPressureGas =  rho * Temp *  kbNormConst / (atomicMass * MpNormConst)

      IF (ANY(getPressureGas <= 0.0_xp)) THEN
         PRINT *, "negative gas pressure"
      END IF

   END FUNCTION


!---------------------------------------------------------------------------------------------
   !> calcule la pression de radiation Eq(2.3.3)
   FUNCTION getPressureTotal(PressureGas, PressureRad)

      IMPLICIT NONE

      REAL(kind=xp), DIMENSION(:)                , INTENT(in) :: PressureGas
      REAL(kind=xp), DIMENSION(:)                , INTENT(in) :: PressureRad
      REAL(kind=xp), DIMENSION(size(pressureRad))             :: getPressureTotal

      getPressureTotal = PressureRad + PressureGas

   END FUNCTION


!---------------------------------------------------------------------------------------------
   !> calcule le parametre beta qui la pression du gaz sur la pression total Eq(2.3.4)
   FUNCTION getBeta(PressureTotal, PressureGas)

      IMPLICIT NONE

      REAL(kind=xp), DIMENSION(:)                  , INTENT(in) :: PressureTotal
      REAL(kind=xp), DIMENSION(:)                  , INTENT(in) :: PressureGas
      REAL(kind=xp), DIMENSION(size(PressureTotal))             :: getBeta


      getBeta = PressureGas / PressureTotal

  END FUNCTION


!---------------------------------------------------------------------------------------------
   !> calcule la vitesse son Eq(2.3.5)
   FUNCTION getSoundSpeed(HalfLenght,AngularSpeed)

      IMPLICIT NONE

      REAL(kind=xp), DIMENSION(:)           , INTENT(in) :: HalfLenght
      REAL(kind=xp), DIMENSION(:)           , INTENT(in) :: AngularSpeed
      REAL(kind=xp), DIMENSION(size(AngularSpeed))       :: getSoundSpeed


      getSoundSpeed = AngularSpeed * HalfLenght

      IF (ANY(getSoundSpeed <= 0.0_xp)) THEN
         PRINT *, "negative sound speed"
      END IF

   END FUNCTION


!---------------------------------------------------------------------------------------------
   !> calcule la demi hauteur du disque Eq (2.3.6)
   FUNCTION getHalfLength(radius, AngularSpeed, Temp, SurfaceDensity)

      USE constantPhysic
      IMPLICIT NONE

      REAL(kind=xp), DIMENSION(:)                 , INTENT(in) :: radius
      REAL(kind=xp), DIMENSION(:)                 , INTENT(in) :: AngularSpeed
      REAL(kind=xp), DIMENSION(:)                 , INTENT(in) :: Temp
      REAL(kind=xp), DIMENSION(:)                 , INTENT(in) :: SurfaceDensity
      REAL(kind=xp), DIMENSION(size(AngularSpeed))             :: getHalfLength
      REAL(kind=xp), DIMENSION(size(AngularSpeed))             :: A
      REAL(kind=xp), DIMENSION(size(AngularSpeed))             :: B
      REAL(kind=xp), DIMENSION(size(AngularSpeed))             :: C
      REAL(kind=xp), DIMENSION(size(AngularSpeed))             :: Delta

      INTEGER :: i

      IF ( ANY(radius==0.0_xp) .OR. AtomicMass==0.0_xp ) THEN
         WRITE(*,*)"ERROR : unexpected value for radius and/or Atomic mass"
         RETURN
      ENDIF

      A = 0.5_xp * (AngularSpeed**2._xp) * (surfaceDensity/radius)
      B = -(1.0_xp/3.0_xp) * aNormConst * (Temp**4._xp)
      C = -RNormConst * Temp * surfaceDensity / (2.0_xp * AtomicMass * radius)

      Delta = B**2._xp - 4._xp*A*C

      IF ( all(A>=0) .AND. all(B<=0) .AND. all(c<=0) )THEN
         IF (all(Delta >= 0)) Then
           getHalfLength = -0.5_xp * (B - sqrt(Delta)) / A
        ELSE
           print *, "ERROR negative DELTA"
        END IF
      ELSE
       DO i = 1,nspaceStep
         IF ( A(i)<=0 .OR. B(i)>=0 .OR. C(i)>=0 )THEN
           print *,'space step', i,"Error on A,B,C", A(i),B(i),C(i)
         END IF
       ENDDO
      END IF
   END FUNCTION


!---------------------------------------------------------------------------------------------
   FUNCTION getRho(radius,surfaceDensity, HalfLength)

      IMPLICIT NONE
      REAL(kind=xp), DIMENSION(:)              , INTENT(in)  :: radius
      REAL(kind=xp), DIMENSION(:)               , INTENT(in) :: HalfLength
      REAL(kind=xp), DIMENSION(:)               , INTENT(in) :: SurfaceDensity
      REAL(kind=xp), DIMENSION(size(HalfLength))             :: getRho

      getRho = surfaceDensity / (2.0_xp * radius * HalfLength)

      IF (ANY(getRho <= 0.0_xp)) THEN
         PRINT *, "negative rho"
      END IF

   END FUNCTION


!---------------------------------------------------------------------------------------------
   !> calcule de la viscosité cinematique Eq(2.3.8)
   FUNCTION getViscosity(SoundSpeed, HalfLength)

      IMPLICIT NONE
      !REAL(kind=xp), INTENT(in)                :: alpha
      REAL(kind=xp), DIMENSION(:)               , INTENT(in)     :: SoundSpeed
      REAL(kind=xp), DIMENSION(:)               , INTENT(in)     :: HalfLength
      REAL(kind=xp), DIMENSION(size(HalfLength))                 :: getViscosity

      getViscosity =  (2.0_xp/3.0_xp) * alpha * SoundSpeed * HalfLength

      IF (ANY(getViscosity <= 0.0_xp)) THEN
         PRINT *, "negative viscosity"
      END IF

   END FUNCTION


!---------------------------------------------------------------------------------------------
  !Solving accretion speed
  FUNCTION getAccretionSpeed(S,nu, radius)
    IMPLICIT NONE
    REAL(kind = xp), DIMENSION(:)         , INTENT(IN)    :: nu, &  !>viscosity list
                                                             S, &      !>surface density
                                                             radius
    REAL(kind = xp), DIMENSION(nspaceStep)                :: getAccretionSpeed
    INTEGER                                               :: i      !>some loop integers
    !------------------------------------------

    ! Boundary condition at rmin dirchlet with S = 0
    getAccretionSpeed(1) = -3._xp * (S(2) * nu(2) - S(1)*nu(1)) / (2._xp * radius(1) * dx * S(1))
    ! Boundary condition at rmax Newmann with Mdot = Mdot0
    getAccretionSpeed(nspaceStep) = - 1. / (2._xp * piConst * xRadius(nspaceStep) * S(nspaceStep))

    DO i = 2, (nSpaceStep-1)        !for all space step except first and last
      getAccretionSpeed(i) = -3._xp * (S(i+1)*nu(i+1)-S(i)*nu(i)) / (2._xp * radius(i) * dx * S(i))
    END DO

  END FUNCTION getAccretionSpeed


!---------------------------------------------------------------------------------------------
   !> calcule le taux accretion qui depends de r Eq(2.3.11)
   FUNCTION getAccretionRate(radius, SurfaceDensity, AccretionSpeed)

      IMPLICIT NONE

      REAL(kind=xp), DIMENSION(:)           , INTENT(in) :: radius
      REAL(kind=xp), DIMENSION(:)           , INTENT(in) :: SurfaceDensity
      REAL(kind=xp), DIMENSION(:)           , INTENT(in) :: AccretionSpeed
      REAL(kind=xp), DIMENSION(size(radius))             :: getAccretionRate


      getAccretionRate = -2._xp * piConst * radius * SurfaceDensity * AccretionSpeed

   END FUNCTION


!---------------------------------------------------------------------------------------------
   !> calcule opacité de opacityFree
   FUNCTION getOpacityFree(Temp, rho)

      IMPLICIT NONE

      REAL(kind=xp), DIMENSION(:)         , INTENT(in) :: Temp
      REAL(kind=xp), DIMENSION(:)         , INTENT(in) :: rho
      REAL(kind=xp), DIMENSION(size(Temp))             :: getOpacityFree

      getOpacityFree = OpacityFreeNormConst * rho / (Temp**3.5_xp)

      IF (ANY(getOpacityFree <= 0.0_xp)) THEN
         PRINT *, "negative opacity free"
      END IF

   END Function


!---------------------------------------------------------------------------------------------
   !> calcule emissivité free free
   FUNCTION getEmissivityFree(rho, temp)

      IMPLICIT NONE

      REAL(kind=xp), DIMENSION(:)         , INTENT(in) :: Temp
      REAL(kind=xp), DIMENSION(:)         , INTENT(in) :: rho
      REAL(kind=xp), DIMENSION(size(temp))             :: getEmissivityFree

      getEmissivityFree = emissivityFreeNormConst * (rho**2._xp) * sqrt(Temp)

      IF (ANY(getEmissivityFree <= 0.0_xp)) THEN
         PRINT *, "negative emissivity free"
      END IF

   END FUNCTION


!---------------------------------------------------------------------------------------------
!> calcule profondeur optique effective
   FUNCTION getTaueff( opacityFree, surfaceDensity, radius)
      IMPLICIT NONE

      REAL(kind=xp), DIMENSION(:)           , INTENT(in) :: SurfaceDensity
      REAL(kind=xp), DIMENSION(:)           , INTENT(in) :: radius
      REAL(kind=xp), DIMENSION(:)           , INTENT(in) :: opacityFree
      REAL(kind=xp), DIMENSION(size(radius))             :: getTaueff

      getTaueff = 0.5_xp * surfaceDensity * sqrt(thomsonOpacityNorm * opacityFree) / radius

      IF (ANY(getTaueff <= 0.0_xp)) THEN
         PRINT *, "negative taueff"
      END IF

   END FUNCTION


!---------------------------------------------------------------------------------------------
   !> equation du flux radiatif Eq(2.3.13)
   FUNCTION getRadiatifFlux(Temp, radius, opacityFree, sigma, HalfLenght, taueff, emissivityFree)

      IMPLICIT NONE
      REAL(kind=xp), DIMENSION(:)         , INTENT(in) :: Temp
      REAL(kind=xp), DIMENSION(:)         , INTENT(in) :: radius
      REAL(kind=xp), DIMENSION(:)         , INTENT(in) :: opacityFree
      REAL(kind=xp), DIMENSION(:)         , INTENT(in) :: sigma
      REAL(kind=xp), DIMENSION(:)         , INTENT(in) :: HalfLenght
      REAL(kind=xp), DIMENSION(:)         , INTENT(in) :: taueff
      REAL(kind=xp), DIMENSION(:)         , INTENT(in) :: emissivityFree
      REAL(kind=xp), DIMENSION(size(Temp))             :: getRadiatifFlux
      INTEGER                                          :: i

      DO i = 1, size(Temp)
         IF (taueff(i) >= 1._xp) THEN
            getRadiatifFlux(i) = 2.0_xp*aNormConst*cNormConst*radius(i)*(Temp(i)**4._xp)/&
            (sigma(i)*3._xp*(opacityFree(i)+thomsonOpacityNorm))
         ELSE
            getRadiatifFlux(i) = emissivityFree(i) *  HalfLenght(i)
         END IF
      END DO

      IF (ANY(getRadiatifFlux <= 0.0_xp)) THEN
         PRINT *, "negative radiative flow"
      END IF

    END FUNCTION


!---------------------------------------------------------------------------------------------
   !> calcule de Q plus Eq(2.3.13 14c)
   FUNCTION getQplus(viscosity, AngularSpeed)

      IMPLICIT NONE

      REAL(kind=xp), DIMENSION(:)                 , INTENT(in) :: viscosity
      REAL(kind=xp), DIMENSION(:)                 , INTENT(in) :: AngularSpeed
      REAL(kind=xp), DIMENSION(size(angularspeed))             :: getQplus

      getQplus = 2.25_xp * viscosity * (AngularSpeed**2._xp)

      IF (ANY(getQplus <= 0.0_xp)) THEN
         PRINT *, "negative Q+"
      END IF

   END FUNCTION


!---------------------------------------------------------------------------------------------
     !> calcule de Q moins Eq(2.3.13 14c)
    FUNCTION getQmoins(radius, RadiatifFlux, SurfaceDensity)

        IMPLICIT NONE

        REAL(kind=xp), DIMENSION(:)           , INTENT(in) :: radius
        REAL(kind=xp), DIMENSION(:)           , INTENT(in) :: RadiatifFlux
        REAL(kind=xp), DIMENSION(:)           , INTENT(in) :: SurfaceDensity
        REAL(kind=xp), DIMENSION(size(radius))             :: getQmoins

        getQmoins = 2._xp * radius * RadiatifFlux / SurfaceDensity

      IF (ANY(getQmoins <= 0.0_xp)) THEN
         PRINT *, "negative Q-"
      END IF

    END FUNCTION


!---------------------------------------------------------------------------------------------
   !>Compute the difference from Qplus, Qmoins and all their parameters
   FUNCTION getDeltaQ(viscosity, AngularSpeed, radius, RadiatifFlux, SurfaceDensity)
      !Wilfried - 09/11/18

      IMPLICIT NONE

      REAL(KIND = xp), DIMENSION(:), INTENT(IN) :: viscosity,  AngularSpeed, radius, RadiatifFlux, SurfaceDensity !Arguments of the functions
      REAL(KIND = xp), DIMENSION(SIZE(radius))  :: getDeltaQ

      getDeltaQ = getQplus(viscosity, AngularSpeed) - getQmoins(radius, RadiatifFlux, SurfaceDensity)

   END FUNCTION getDeltaQ


!---------------------------------------------------------------------------------------------
   FUNCTION getHeatCapacity(beta)

      IMPLICIT NONE

      REAL(kind=xp), DIMENSION(:)         , INTENT(in) :: Beta
      REAL(kind=xp), PARAMETER                         :: gammaG = 5./3.
      REAL(kind=xp), DIMENSION(size(beta))             :: getHeatCapacity

      getHeatCapacity = (RNormConst * (12.0_xp * (gammaG - 1._xp)&
      * ( 1._xp - beta ) + beta )) / ((gammaG - 1._xp) * beta * atomicMass)

      IF (ANY(getHeatCapacity <= 0.0_xp)) THEN
         PRINT *, "negative heat capacity"
      END IF

   END FUNCTION


!---------------------------------------------------------------------------------------------
   FUNCTION getGamma3(heatcapacity, beta)! la fonction calcule gamma3 -1

        IMPLICIT NONE

        REAL(kind=xp), DIMENSION(:)         , INTENT(in) :: Beta
        REAL(kind=xp), DIMENSION(:)         , INTENT(in) :: heatCapacity
        REAL(kind=xp), DIMENSION(size(beta))             :: getGamma3

        getGamma3 = (RNormConst*(4._xp - 3._xp*beta) / (beta*AtomicMass*heatCapacity))

    END FUNCTION


!---------------------------------------------------------------------------------------------
   !> Function solving Qadv
   FUNCTION getAdvectionHeatFlow(Temp, S, v, nu, Cv , gamma3, radius)
     !Valentin
     !14 novembre
     !basé sur le poly "equations.tex"

      IMPLICIT NONE

      REAL(kind = xp), DIMENSION(:)         , INTENT(IN) :: S, & !>surface density
                                                            v, &     !>speed list
                                                            nu, & !viscosity
                                                            Temp, &
                                                            Cv, &
                                                            gamma3, &
                                                            radius
      REAL(kind = xp)                                       dSoxdx, & !S/x space derivative
                                                            dTdx, & !Temperature space derivative
                                                            ddnuSdx2, & !second space derivative of nu*S
                                                            derivsum !just to make the line shorter
      REAL(kind = xp), DIMENSION(nspaceStep)             :: getAdvectionHeatFLow

      INTEGER                                            :: i        !> loop integer

      ! Boundary condition at rmax, we take value(imax+1) = 2value(imax) - value(imax-1)
      dTdx = (Temp(nSpaceStep)-Temp(nSpaceStep-1))/dx
      dSoxdx = (S(nSpaceStep)/radius(nSpaceStep)-S(nSpaceStep-1)/radius(nSpaceStep-1))/dx
      derivsum = gamma3(nSpaceStep)*Temp(nSpaceStep)/S(nSpaceStep)*v(nSpaceStep)*dSoxdx
      getAdvectionHeatFLow(nspaceStep) = Cv(nSpaceStep)/2._xp*(derivsum - v(nSpaceStep)/radius(nSpaceStep)*dTdx)

      ! Boundary condition at rmin, the second derivative is zero (see Valentin's Poly)
      dTdx = (Temp(2)-Temp(1))/dx
      dSoxdx = (S(2)/radius(2)-S(1)/radius(1))/dx
      derivsum = gamma3(1)*Temp(1)/S(1)*v(1)*dSoxdx
      getAdvectionHeatFLow(1) = Cv(1)/2._xp*(derivsum - v(1)/radius(1)*dTdx)

      DO i = 2, nSpaceStep-1        !for all space step except first and last
        dTdx = (Temp(i+1)-Temp(i))/dx
        dSoxdx = (S(i+1)/radius(i+1)-S(i)/radius(i))/dx
        ddnuSdx2 = (nu(i+1)*S(i+1) - 2._xp*nu(i)*S(i) + nu(i-1)*S(i-1))/dx**2_xp
        derivsum = gamma3(i)*Temp(i)/S(i)*(3._xp/(2._xp*radius(i)**2_xp)*ddnuSdx2 + v(i)*dSoxdx)
        getAdvectionHeatFLow(i) = Cv(i)/2._xp*(derivsum - v(i)/radius(i)*dTdx)
      END DO
      getAdvectionHeatFlow = 0._xp
   END FUNCTION getAdvectionHeatFlow


   SUBROUTINE algebricStep(PressureRad, &
                            angularSpeed, &
                            halfLength, &
                            soundSpeed, &
                            viscosity, &
                            Qplus, &
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
                            density, &
                            Temperature)
     IMPLICIT NONE

     REAL(KIND = xp), DIMENSION(nSpaceStep), INTENT(OUT) :: PressureRad, &
                                                            angularSpeed, &
                                                            halfLength, &
                                                            soundSpeed, &
                                                            viscosity, &
                                                            Qplus, &
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
                                                            rho
     REAL(KIND = xp), DIMENSION(nSpaceStep), INTENT(IN) :: Temperature, density


     !Calculate angular speed and halflength
     angularSpeed = getAngularspeed(xRadius)
     halfLength = getHalfLength(xRadius,angularSpeed,Temperature,density)

     !Calculate speed of sound
     soundSpeed = getSoundSpeed(halflength,angularSpeed)

     !Détermination de la masse volumique
     rho = getRho(xRadius,density,halfLength)

     !Calculate viscosity
     viscosity = getViscosity(soundSpeed, halflength)

     !Calculate Qplus
     Qplus = getQplus(viscosity, angularSpeed)

     !Calculate gas pressure
     pressureGas = getPressureGas(rho,Temperature)

     !calculate opacity free free
     opacityFree = getOpacityFree(Temperature, rho)

     !calculate radiation pressure
     PressureRad = getPressureRad(Temperature)

     !calculate total pressure
     pressureTotal = getPressureTotal(pressureGas,PressureRad)

     !Calculate beta the ration of pressure gas on total pressure
     beta = getBeta(pressureTotal,pressureGas)

     !Calculate heat capacity
     heatCapacity = getHeatCapacity(beta)

     !Calculate emissivity free free
     epsilonff = getEmissivityFree(density, Temperature)

     !Calculate optical depth
     taueff = getTaueff(opacityFree, density, xRadius)

     !Calculate radiation emission
     Fz = getRadiatifFlux(Temperature, xRadius, opacityFree, density, halflength, taueff, epsilonff)

     !Calculate Q moins
     Qm = getQmoins(xRadius, Fz, density)

     !calculate accretion speed
     accretionSpeed = getAccretionSpeed(density,viscosity,xRadius)

     !calculate accretion rate
     accretionRate = getAccretionRate(xRadius,density,accretionSpeed)

     !calculate gamma3-1
     gamma3 = getGamma3(heatCapacity, beta)

     !Calculate Qadv
     Qadv = getAdvectionHeatFlow(Temperature, density, accretionSpeed, viscosity, heatCapacity, gamma3, xRadius)

   END SUBROUTINE algebricStep

END MODULE
