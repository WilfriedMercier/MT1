MODULE diffSolver

   USE constantProgram
   USE constantSimulation

   IMPLICIT NONE

   REAL(kind=xp), PUBLIC :: dt, &!> TimeStep
                            condStabilityDiff, & !> Stability condition for diffusion explicit scheme
                            condStabilityAdv     !> Stability condition for advection explicit scheme
   REAL(kind=xp), PUBLIC :: viscT, &  !>Viscous time
                            thermT, & !> Thermic time
                            dynT      !> Dynamic time



   PRIVATE
   PUBLIC timeStepCalculator, &
          viscousTime, &
          thermicTime, &
          dynamicTime, &
          dtFix, &
          dxFix, &
          temperatureSolver, &
          temperatureExplicitScheme, &
          temperatureImplicitScheme, &
          densityExplicit, &
          densityImplicit, &
          densitySolver

   CONTAINS

      !> Determine time step (correspond to the lowest caracteristic time)
      FUNCTION timeStepCalculator(v, nu, x, h, omega)

         IMPLICIT NONE

         REAL(kind=xp), DIMENSION(:) , INTENT(IN)  :: v, &   !>Speed of accretion
                                                      nu, &  !>Viscosity
                                                      x, &   !>Normalised radius
                                                      h, &   !> Heigt of the disk
                                                      omega  !> Rotation Speed
         REAL(kind=xp)                             :: timeStepCalculator

         viscT = viscousTime(nu,x)
         thermT = thermicTime(x,h,viscT)
         dynT = dynamicTime(omega)
         condStabilityAdv = getStabilityAdvCond(v,dx)
         condStabilityDiff = getStabilityDiffCond(dx)

         timeStepCalculator = MIN(MIN(viscT,thermT),dynT)

      END FUNCTION


      !> Viscous time calculation (Page 17 3.3.1)
      FUNCTION viscousTime(nu, x)
         IMPLICIT NONE

         REAL(kind=xp), DIMENSION(:) , INTENT(IN) :: nu, & !>viscosity
                                                     x    !>Normalised radius
         REAL(KIND=xp)                            :: viscousTime
         INTEGER                                  :: n
         n = SIZE(nu)

         viscousTime = MINVAL(x**4_xp/nu)

      END FUNCTION viscousTime


      !> Thermic time calculation (Page 17 3.3.2)
      FUNCTION thermicTime(x, h, viscTime)

         IMPLICIT NONE

         REAL(kind=xp), DIMENSION(:) , INTENT(IN) :: x, &          !>Normalized radius
                                                     h             !> Heigt of the disk
         REAL(KIND=xp),                INTENT(IN) :: viscTime   !> Viscous Time
         REAL(KIND=xp)                            :: thermicTime


         thermicTime = MINVAL((h/x**2._xp)**2._xp*viscTime)

      END FUNCTION thermicTime


      !> Dynamic time calculation (Page 18 3.3.3)
      FUNCTION dynamicTime(omega)

         IMPLICIT NONE

         REAL(kind=xp), DIMENSION(:) , INTENT(IN) :: omega       !> Rotation Speed
         REAL(KIND=xp)                            :: dynamicTime !> Dynamic Time

         dynamicTime = MINVAL(1._xp/omega)

      END FUNCTION dynamicTime


      !>Calculate advection stability condition
      FUNCTION getStabilityAdvCond(v, dx)

         REAL(kind=xp), DIMENSION(:) , INTENT(IN) :: v !>Speed of accretion
         REAL(kind=xp),                INTENT(IN) :: dx
         REAL(kind=xp)                            :: getStabilityAdvCond

         getStabilityAdvCond = 2._xp*dx/(3._xp*maxval(abs(v)))

      END FUNCTION getStabilityAdvCond

      !>Calculate diffusion stability condition
      FUNCTION getStabilityDiffCond(dx)

         REAL(kind=xp),                INTENT(IN) :: dx
         REAL(kind=xp)                            :: getStabilityDiffCond

         getStabilityDiffCond = 2._xp*dx*dx

      END FUNCTION getStabilityDiffCond

      !>SUBROUTINE for testing purpose fix TimeStep
      SUBROUTINE dtFix(deltaT)

         IMPLICIT NONE

         REAL(kind=xp), INTENT(IN) :: deltaT !>Timestep

         dt = deltaT

      END SUBROUTINE dtFix


      SUBROUTINE dxFix(deltax)

         IMPLICIT NONE

         REAL(kind=xp), INTENT(IN) :: deltax !>Space step

         dx = deltax

      END SUBROUTINE dxFix

!---------------------------------------------------------------------------------
!Temperature solving

      !>subroutine which resolves the partial differential equation for the temperature
      SUBROUTINE temperatureSolver(T,Qp,Qm,Qadv,Cv)

         IMPLICIT NONE

         REAL(kind = xp), DIMENSION(:), INTENT(INOUT) :: T        !temperature list
         REAL(kind = xp), DIMENSION(:), INTENT(IN)    :: Qp, &   !incoming heat flow
                                                         Qm, &   !outcoming heat flow
                                                         Qadv, &    !advection heat flowdt, &    !time interval
                                                         Cv

         CALL temperatureExplicitScheme(T,Qp,Qm,Qadv,Cv)

      END SUBROUTINE temperatureSolver

      !> Subroutine for solving temperature with an explicit scheme
      SUBROUTINE temperatureExplicitScheme(T,Qp,Qm,Qadv,Cv)

         IMPLICIT NONE

         REAL(kind = xp), DIMENSION(:), INTENT(INOUT) :: T        !temperature list
         REAL(kind = xp), DIMENSION(:), INTENT(IN)    :: Qp, &   !incoming heat flow
                                                         Qm, &   !outcoming heat flow
                                                         Qadv, &    !advection heat flowdt, &    !time interval
                                                         Cv
         REAL(kind = xp), DIMENSION(SIZE(T))          :: oldT !temperature copy

         oldT = T
         T = oldT+dt*(Qp-Qm+Qadv)/Cv


      END SUBROUTINE temperatureExplicitScheme

      !> Subroutine for solving  temperature with an implicit scheme
      SUBROUTINE temperatureImplicitScheme(T,Qp,Qm,Qadv, Cv)

         IMPLICIT NONE

         REAL(kind = xp), DIMENSION(:), INTENT(INOUT) :: T        !temperature list
         REAL(kind = xp), DIMENSION(:), INTENT(IN)    :: Qp, &   !incoming heat flow
                                                         Qm, &   !outcoming heat flow
                                                         Qadv, &    !advection heat flowdt, &    !time interval
                                                         Cv

         CALL temperatureExplicitScheme(T,Qp,Qm,Qadv,Cv)

      END SUBROUTINE temperatureImplicitScheme



!---------------------------------------------------------------------------------
!Density solving


      SUBROUTINE densitySolver(S,nu)

         IMPLICIT NONE

         REAL(kind = xp), DIMENSION(:), INTENT(IN)    :: nu     !>viscosity list
         REAL(kind = xp), DIMENSION(:), INTENT(INOUT) :: S      !>density list to be changed

         !CALL densityCN(S,nu)
         IF (dt > condStabilityDiff) THEN
            CALL densityImplicit(S,nu)
         ELSE
            CALL densityExplicit(S,nu)
         END IF

      END SUBROUTINE densitySolver


      SUBROUTINE densityImplicit(S,nu)

        !> Valentin 3/12 for the latest modifs

         IMPLICIT NONE

         REAL(kind = xp), DIMENSION(:), INTENT(IN)    :: nu         !>viscosity list
         REAL(kind = xp), DIMENSION(:), INTENT(INOUT) :: S          !>density list to be changed
         REAL(kind = xp), DIMENSION(nspaceStep)       :: B, D       !>Element for lapack call
         REAL(kind = xp), DIMENSION(nspaceStep-1)     :: DL, DU     !>Element for lapack call
         REAL(kind = xp), DIMENSION(nspaceStep)       :: lmdExp     !>lambda factor (see method)
         INTEGER                                      :: N, NRHS, LDB, INFO !> lapack inout info

         !The general idea behind this function is to pose A*X = B with x the new value of S

         !> See the rapport (implicit Euler method) for the details of the B,D,DL,DU filling
         !Store the actual value into B
         B = S
          !Calculate lambda
         lmdExp = 3._xp*dt/(4._xp*dx**2._xp*xRadius**2._xp)

         !Fill D, Dl and Du
         D  = 1._xp+2._xp*lmdExp*nu
         DU = -lmdExp(1:(nSpaceStep-1))*nu(2:nSpaceStep)
         DL = -lmdExp(2:nSpaceStep)*nu(1:(nSpaceStep-1))

         !>Previous boundary conditions
         D(nSpaceStep)    = nu(nSpaceStep)/dx
         DL(nSpaceStep-1) = -nu(nSpaceStep-1)/dx
         B(nSpaceStep)    = 1._xp/(3._xp*piConst)

         !Init other value
         N = nSpaceStep
         NRHS = 1
         LDB = nSpaceStep
         INFO = 0
         !Solve the linear equation
         CALL DGTSV(N, NRHS, DL, D, DU, B, LDB, INFO)
         IF (INFO /= 0) THEN
            PRINT *, 'On density implicit, error ', INFO
         END IF
         S = B

      END SUBROUTINE densityImplicit


      SUBROUTINE densityCN(S,nu)


         IMPLICIT NONE

         REAL(kind = xp), DIMENSION(:), INTENT(IN)    :: nu     !>viscosity list
         REAL(kind = xp), DIMENSION(:), INTENT(INOUT) :: S      !>density list to be changed
         REAL(kind = xp), DIMENSION(nspaceStep)       :: B, D, L1, L2 !>Element for lapack call
         REAL(kind = xp), DIMENSION(nspaceStep-1)     :: DL, DU !>Element for lapack call
         REAL(kind = xp), DIMENSION(nspaceStep)       :: lmdExp !>lambda factor (see method)
         REAL(kind = xp)                              :: beta = 0.5_xp
         INTEGER                                      :: N, NRHS, LDB, INFO

         !The general idea behind this function is to pose A*X = B with x the new value of S

         !> See the rapport (implicit Euler method) for the details of the B,D,DL,DU filling
         !Store the actual value into B
         B = S
         !Calculate lambda
         lmdExp = 3._xp*dt/(4._xp*(dx**2._xp)*xRadius**2._xp)

         !Fill D, Dl and Du for explicit
         D  = 1._xp+2._xp*lmdExp*nu*(1._xp-beta)
         DU = -lmdExp(1:(nSpaceStep-1))*nu(2:nSpaceStep)*(1._xp-beta)
         DL = -lmdExp(2:nSpaceStep)*nu(1:(nSpaceStep-1))*(1._xp-beta)

         !Previous boundary conditions
         D(nSpaceStep)    = -lmdExp(nSpaceStep-1)*nu(nSpaceStep-1)*(1._xp-beta)/nu(nSpaceStep)
         DL(nSpaceStep-1) = (1._xp+2._xp*lmdExp(nSpaceStep-1)*(1._xp-beta))*nu(nSpaceStep-1)/nu(nSpaceStep)

         B = tridiagonalMVmul(DL,DU,D,B)
         B(nSpaceStep) = B(nSpaceStep) + (dx/(3._xp*piConst*nu(nSpaceStep)))*(1._xp-beta) &
         + beta*(1._xp/(3._xp*piConst)-S(nSpaceStep))

         !Fill D, Dl and Du
         D  = 1._xp+2._xp*lmdExp*nu*beta
         DU = -lmdExp(1:(nSpaceStep-1))*nu(2:nSpaceStep)*beta
         DL = -lmdExp(2:nSpaceStep)*nu(1:(nSpaceStep-1))*beta

         !Previous boundary conditions
         D(nSpaceStep)    = nu(nSpaceStep)/dx
         DL(nSpaceStep-1) = -nu(nSpaceStep-1)/dx
         !B(nSpaceStep)    = 1._xp/(3._xp*piConst)

         !Init other value
         N = nSpaceStep
         NRHS = 1
         LDB = nSpaceStep
         INFO = 0
         !Solve the linear equation
         CALL DGTSV(N, NRHS, DL, D, DU, B, LDB, INFO)
         IF (INFO /= 0) THEN
            PRINT *, 'On density implicit, error ', INFO
         END IF
         S = B

      END SUBROUTINE densityCN


      !>explicit method for resolving the partial differential equation on the density
      SUBROUTINE densityExplicit(S,nu)
         IMPLICIT NONE

         REAL(kind = xp), DIMENSION(:), INTENT(IN)    :: nu     !>viscosity list
         REAL(kind = xp), DIMENSION(:), INTENT(INOUT) :: S      !>density list to be changed
         REAL(kind = xp), DIMENSION(nspaceStep)       :: newS   !>S copy then new value
         REAL(kind = xp)                              :: lmdExp !>lambda factor (see method)
         INTEGER                                      :: i      !>some loop integers
         newS = S
         lmdExp = 3._xp/4._xp*dt/(dx**2._xp)
         ! Boundary condition at rmin
         newS(1) = S(1)+ (lmdExp/(xRadius(1)**2._xp))*(nu(2)*S(2)-2._xp*nu(1)*S(1))

         !for all space step except first and last
         DO i = 2, nSpaceStep-1
           newS(i) = S(i) + (lmdExp/(xRadius(i)**2._xp))*(nu(i-1)*S(i-1) - 2._xp*nu(i)*S(i) + nu(i+1)*S(i+1))
         ENDDO

         ! Boundary condition at rmax : dSnu is the condition so we differentiate dSnu once
         newS(nspaceStep) = (dx/(3._xp*piConst) + nu(nSpaceStep-1)*newS(nSpaceStep-1))/nu(nSpaceStep)


         S = newS ! New value affected to S, at the time t + dt
      END SUBROUTINE densityExplicit


      !> Tridiagonal multiplication
      FUNCTION tridiagonalMVmul(DL,DU,D,B)

         IMPLICIT NONE

         REAL(kind = xp), DIMENSION(:), INTENT(IN) :: DL,DU,D,B
         REAL(kind = xp), DIMENSION(size(B))       :: tridiagonalMVmul
         INTEGER                                   :: i

         tridiagonalMVmul(1) = D(1)*B(1) + DU(1)*B(2)
         DO i = 2,(size(B)-1)
            tridiagonalMVmul(i) = DL(i-1)*B(i-1) + D(i)*B(i) + DU(i)*B(i+1)
         END DO
         tridiagonalMVmul(size(B)) = DL(size(B)-1)*B(size(B)-1) + D(size(B))*B(size(B))

      END FUNCTION

END MODULE diffSolver
