MODULE sCurve

   USE constantProgram
   USE constantSimulation
   USE algebricEquations
   USE constantPhysic

   IMPLICIT NONE


   CONTAINS

   !--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
      !> Finds a zero by dichotomy within a certain range
   SUBROUTINE ComputeDicho(fonction, X0, Y0, X1, Y1, XOut, YOut, Eps, dimi, TauEff, radius, isOkFlag)

      REAL(KIND = xp), EXTERNAL        :: fonction                         !Function used
      REAL(KIND = xp), INTENT(IN)      :: Eps                              !Accuracy wanted for the dichotomy
      REAL(KIND = xp), INTENT(IN)      :: radius                           !radius
      INTEGER, INTENT(IN)              :: dimi                             !Dimension along which the dichotomy is performed
      REAL(KIND = xp), INTENT(INOUT)   :: X0, Y0, X1, Y1                   !Positions of the ends
      INTEGER, INTENT(INOUT)           :: isOkFlag                         !Flag to know if the dichotomy converged or not (0 if it did not, 1 if it did)
      REAL(KIND = xp), INTENT(OUT)     :: XOut, YOut                       !X, Y position of the zero
      REAL(KIND = xp), INTENT(OUT)     :: TauEff                           !Opacity

      ! ... local ...
      REAL(KIND = xp)                  :: XMean, YMean                     !Temporary positions
      REAL(KIND = xp)                  :: FMinus, FPlus                    !Values of the function at the two ends
      REAL(KIND = xp)                  :: FMean                            !Temporary variable for the dichotomy
      REAL(KIND = xp)                  :: accuracy_abs                     !Absolute accuracy of the dichotomy
      REAL(KIND = xp)                  :: accuracy_rel                     !Relative accuracy of the dichotomy
      INTEGER                          :: N                                !Number of iterations performed
      INTEGER, PARAMETER               :: NMax = 1000000                   !Maximum number of iterations
      !-------------------------------------------------------------------------------------------------------



      isOkFlag       = 1
      N              = 0
      accuracy_abs   = 1000.0_xp
      accuracy_rel   = accuracy_abs
      FMinus         = 0.0_xp
      FPlus          = 0.0_xp

      Loop_foundit : DO WHILE( (N<NMax) .AND. (accuracy_abs>Eps) .AND. (accuracy_rel>Eps) )  !As long as the accuracy is too high and we did not reach the maximum iterations, perform a new dichotomy within a more restricted range
         N = N+1
         IF (Eps>1.0E-2_xp) THEN
            WRITE(*,*)"Larger than 0.01 at ", N
         ENDIF

         !Compute two values of the function at (X,Y) and either (X,Y+DELTA) or (X+DELTA, Y) given the value of i
         FMinus   = fonction(X0, Y0, TauEff, radius)
         FPlus    = fonction(X1, Y1, TauEff, radius)

         !New accuracy
         accuracy_abs = ABS((Fplus-Fminus)/FPlus)
         accuracy_rel = ABS(Fplus-Fminus)

         !Compute the position and the value of the function at half distance
         XMean = (X0+X1)/2.0_xp
         YMean = (Y0+Y1)/2.0_xp
         FMean = fonction(XMean, YMean, TauEff, radius)

         !Find in which range lies the zero
         cond_limit_sol : IF (Fminus*FMean > 0._xp) THEN
            IF (dimi==0) THEN
               X0 = XMean
            ELSE
               Y0 = YMean
            END IF
         ELSE cond_limit_sol
            IF (dimi==0) THEN
               X1 = XMean
            ELSE
               Y1 = YMean
            END IF
         END IF cond_limit_sol

         XOut = X1
         YOut = Y1

!         WRITE(*,*)X0, Y0, X1, Y1, FPlus !, FMinus, accuracy

      END DO Loop_foundit

      IF ( (accuracy_rel>Eps .AND. accuracy_abs>Eps) .OR. (N>=NMax) ) THEN
         WRITE(*,*)"Not converging :", accuracy_abs, accuracy_rel
         isOkFlag = -2
      ENDIF

   END SUBROUTINE ComputeDicho


   !------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
      !>Look for the first crossing of the sCurve (hence find the range within which the zero lies) by moving two points (one upwards, one downward) and then performs a dichotomy
   SUBROUTINE Dichotomy(fonction, X, Y, TauEff, radius, XMin, YMin, XZero, YZero, Eps, dimi, Delta, NMax, isOkFlag)

      REAL(KIND = xp), EXTERNAL        :: fonction                         !Function used
      REAL(KIND = xp), INTENT(IN)      :: X, Y                             !Initial position
      REAL(KIND = xp), INTENT(IN)      :: Eps                              !Accuracy wanted for the dichotomy
      REAL(KIND = xp), INTENT(IN)      :: Delta                            !Step used in order to find the first crossing of the sCurve
      REAL(KIND = xp), INTENT(IN)      :: Xmin, Ymin                       !Minimum positions allowed
      REAL(KIND = xp), INTENT(IN)      :: radius                  !radius
      INTEGER, INTENT(IN)              :: dimi                             !Dimension along which the dichotomy is performed
      INTEGER, INTENT(IN)              :: NMax                             !Maximum number of iterations allowed
      REAL(KIND = xp), INTENT(INOUT)   :: XZero, YZero                     !Position of the zero found by dichotomy
      INTEGER, INTENT(OUT)             :: isOkFlag                         !Flag to know how well the dichotomy was performed (-2 = did not converge, -1 = wrong value for i/position, 0 = no first crossing of sCurve, 1 = ok)
      REAL(KIND = xp), INTENT(OUT)     :: TauEff                           !Opacity

      ! ...local...
      REAL(KIND = xp)                  :: XOldUp, YOldUp                   !Old positions for the point moving upward
      REAL(KIND = xp)                  :: XNewUp, YNewUp                   !New positions for the point moving upward

      REAL(KIND = xp)                  :: XOldDown, YOldDown               !Old positions for the point moving downward
      REAL(KIND = xp)                  :: XNewDown, YNewDown               !New positions for the point moving downward

      REAL(KIND = xp)                  :: XOld, YOld                       !Old positions for the lower bound within which the zero lies
      REAL(KIND = xp)                  :: XNew, YNew                       !New positions for the uper bound within which the zero lies

      REAL(KIND = xp)                  :: FMinus, FPlus                    !Values of the function at the two ends

      INTEGER                          :: N                                !Number of iterations performed
      INTEGER                          :: dichoFlag                        !Flag used to tell if the Dichotomy can be used
      ! ----------------------------------------------------------------------------

      N           = 0
      FMinus      = 0.0_xp
      FPlus       = 0.0_xp
      dichoFlag   = 0
      isOkFlag    = 1

      !Check we are in the correct ranges/have the correct values
      IF ( (dimi/=0) .AND. (dimi/=1) ) THEN !0->S 1->T
         WRITE(*,*)"Error : i should only be 0 (or 1) to perform dichotomy along X (or Y). Given value is ", dimi
         isOkFlag = -1
         RETURN
      END IF

      IF ( (X<XMin) .OR. (Y<YMin) ) THEN
         WRITE(*,*)"Error : one of the position coordinates was below the minimum allowed."
         isOkFlag = -1
         RETURN
      ENDIF

      !Initialize values
      XOldUp   = X
      YOldUp   = Y
      XOldDown = X
      YOldDown = Y
      IF (dimi==0) THEN
         XNewUp   = X + Delta
         YNewUp   = Y
         XNewDown = X - Delta
         YNewDown = Y
      ELSE IF (dimi==1) THEN
         YNewUp   = Y + Delta
         XNewUp   = X
         YNewDown = Y - Delta
         XNewDown = X
      ENDIF

      !Correct the values if we are out of bounds
      IF (XNewDown < XMin) THEN
         XNewDown = Xmin
      ENDIF
      IF (YNewDown < YMin) THEN
         YNewDown = YMin
      ENDIF

      loop_dichotomy : DO WHILE (N<NMax)
         N = N+1

         !Check for a first crossing between bounds and, if not found, move a step Delta upwards
         CALL checkFirstCrossing(fonction, XOldUp, YOldUp, XNewUp, YNewUp, XMin, YMin, Delta, dimi, &
                                 TauEff, radius, dichoFlag)

         If (dichoFlag==1) THEN
            XNew  = XNewUp
            YNew  = YNewUp
            XOld  = XOldUp
            YOld  = YOldUp
         ELSE
            !Check for a first crossing between bounds and, if not found, move a step Delta downwards (hence -Delta)
            CALL checkFirstCrossing(fonction, XOldDown, YOldDown, XNewDown, YNewDown, XMin, YMin, -Delta, dimi, &
                                    TauEff, radius, dichoFlag)
            If (dichoFlag==1) THEN
               XNew  = XNewDown
               YNew  = YNewDown
               XOld  = XOldDown
               YOld  = YOldDown
            ENDIF
         ENDIF

         !Once the first crossing is found (flag == 1) then call the dichotomy
         If (dichoFlag==1) THEN
            CALL ComputeDicho(fonction, XOld, YOld, XNew, YNew, XZero, YZero, Eps, dimi, TauEff, radius, isOkFlag) !Call the dichotomy within the range we found
            RETURN
         ENDIF

      END DO loop_dichotomy
      WRITE(*,*)"Error : maximum number (", NMax, &
                ") of iterations reached. No first crossing of the sCurve found from position (X, Y) = (", &
                X, Y, ") with total step ", NMax*Delta, " in both directions."
      isOkFlag = 0

   END SUBROUTINE Dichotomy

   SUBROUTINE DichotomyPositiveOnly(fonction, X, Y, TauEff, radius, XMin, YMin, XZero, YZero, Eps, dimi, Delta, NMax, isOkFlag)

      REAL(KIND = xp), EXTERNAL        :: fonction                         !Function used
      REAL(KIND = xp), INTENT(IN)      :: X, Y                             !Initial position
      REAL(KIND = xp), INTENT(IN)      :: Eps                              !Accuracy wanted for the dichotomy
      REAL(KIND = xp), INTENT(IN)      :: Delta                            !Step used in order to find the first crossing of the sCurve
      REAL(KIND = xp), INTENT(IN)      :: Xmin, Ymin                       !Minimum positions allowed
      REAL(KIND = xp), INTENT(IN)      :: radius                  !radius
      INTEGER, INTENT(IN)              :: dimi                             !Dimension along which the dichotomy is performed
      INTEGER, INTENT(IN)              :: NMax                             !Maximum number of iterations allowed
      REAL(KIND = xp), INTENT(INOUT)   :: XZero, YZero                     !Position of the zero found by dichotomy
      INTEGER, INTENT(OUT)             :: isOkFlag                         !Flag to know how well the dichotomy was performed (-2 = did not converge, -1 = wrong value for i/position, 0 = no first crossing of sCurve, 1 = ok)
      REAL(KIND = xp), INTENT(OUT)     :: TauEff                           !Opacity

      ! ...local...
      REAL(KIND = xp)                  :: XOldUp, YOldUp                   !Old positions for the point moving upward
      REAL(KIND = xp)                  :: XNewUp, YNewUp                   !New positions for the point moving upward

      REAL(KIND = xp)                  :: XOld, YOld                       !Old positions for the lower bound within which the zero lies
      REAL(KIND = xp)                  :: XNew, YNew                       !New positions for the uper bound within which the zero lies

      REAL(KIND = xp)                  :: FMinus, FPlus                    !Values of the function at the two ends

      INTEGER                          :: N                                !Number of iterations performed
      INTEGER                          :: dichoFlag                        !Flag used to tell if the Dichotomy can be used
      ! ----------------------------------------------------------------------------

      N           = 0
      FMinus      = 0.0_xp
      FPlus       = 0.0_xp
      dichoFlag   = 0
      isOkFlag    = 1

      !Check we are in the correct ranges/have the correct values
      IF ( (dimi/=0) .AND. (dimi/=1) ) THEN !0->S 1->T
         WRITE(*,*)"Error : i should only be 0 (or 1) to perform dichotomy along X (or Y). Given value is ", dimi
         isOkFlag = -1
         RETURN
      END IF

      IF ( (X<XMin) .OR. (Y<YMin) ) THEN
         WRITE(*,*)"Error : one of the position coordinates was below the minimum allowed."
         isOkFlag = -1
         RETURN
      ENDIF

      !Initialize values
      XOldUp   = X
      YOldUp   = Y
      IF (dimi==0) THEN
         XNewUp   = X + Delta
         YNewUp   = Y
      ELSE IF (dimi==1) THEN
         YNewUp   = Y + Delta
         XNewUp   = X
      ENDIF

      loop_dichotomy : DO WHILE (N<NMax)
         N = N+1

         !Check for a first crossing between bounds and, if not found, move a step Delta upwards
         CALL checkFirstCrossing(fonction, XOldUp, YOldUp, XNewUp, YNewUp, XMin, YMin, Delta, dimi, &
                                 TauEff, radius, dichoFlag)

         If (dichoFlag==1) THEN
            XNew  = XNewUp
            YNew  = YNewUp
            XOld  = XOldUp
            YOld  = YOldUp
         ENDIF

         !Once the first crossing is found (flag == 1) then call the dichotomy
         If (dichoFlag==1) THEN
!            WRITE(*,*)"First crossing found between (", XOld, ",", YOld, ") and (", XNew, ",", YNew, ")"
            CALL ComputeDicho(fonction, XOld, YOld, XNew, YNew, XZero, YZero, Eps, dimi, TauEff, radius, isOkFlag) !Call the dichotomy within the range we found
            !WRITE(*,*)"(XZero, YZero) = (", XZero, ",", YZero, ")"
            RETURN
         ENDIF

      END DO loop_dichotomy
      WRITE(*,*)"Error : maximum number (", NMax, &
                ") of iterations reached. No first crossing of the sCurve found from position (X, Y) = (", &
                X, Y, ") with total step ", NMax*Delta, " along axis ", dimi
      isOkFlag = 0

   END SUBROUTINE DichotomyPositiveOnly


   !------------------------------------------------------------------------------------------------------------------------------------------------------------------
      !>Check that the first crossing is in the range [XOld, XNew], [YOld, YNew]
   SUBROUTINE checkFirstCrossing(fonction, XOld, YOld, XNew, YNew, XMin, YMin, Delta, dimi, TauEff, radius, dichoFlag)

      REAL(KIND = xp), EXTERNAL        :: fonction                !Function used
      REAL(KIND = xp), INTENT(IN)      :: Delta                   !Step used to move upwards or downwards
      REAL(KIND = xp), INTENT(IN)      :: Xmin, Ymin              !Minimum positions allowed
      REAL(KIND = xp), INTENT(IN)      :: radius                  !radius
      INTEGER, INTENT(IN)              :: dimi                    !Dimension along which the dichotomy must be performed
      REAL(KIND = xp), INTENT(INOUT)   :: XOld, YOld, XNew, YNew  !Positions to test
      INTEGER, INTENT(INOUT)           :: dichoFlag               !Flag used to tell if the Dichotomy can be used
      REAL(KIND = xp), INTENT(OUT)     :: TauEff                  !Opacity

      ! ...local...
      REAL(KIND = xp)                  :: FMinus, FPlus           !Values of the function at the two ends
      ! ----------------------------------------------------------------------------

      !Compute two values of the function at (XOld,YOld) and (XNew,YNew)
      FMinus   = fonction(XOld, YOld, TauEff, radius)
      FPlus    = fonction(XNew, YNew, TauEff, radius)

      !Check whether we crossed the sCurve for the first time
      dichoFlag = 0
      IF (FMinus*FPlus > 0) THEN !If we did not cross the sCurve, we keep moving upwards/downwards
         IF (dimi==0) THEN
            XOld  = XNew
            XNew  = XNew+Delta
         ELSE
            YOld  = YNew
            YNew  = YNew+Delta
         ENDIF
      ELSE
         dichoFlag = 1
      END IF

      !Check we are not out of bounds
      IF (XNew < XMin) THEN
         XNew = XMin-Delta
         XOld = XNew-2*Delta
      ENDIF
      IF (YNew < YMin) THEN
         YNew = Ymin-Delta
         YOld = YNew-2*Delta
      ENDIF

   END SUBROUTINE checkFirstCrossing


   !----------------------------------------------------------------------------------------------------------------------------------------------------
      !> Find all the zeros on the sCurve from dichotomy up to a given temperature
   SUBROUTINE sCurveFromDichotomy(fonction, S, T, SMax, TMax, SMin, TMin, Eps, Delta, &
              DeltaS, LogStep, NPointMax, NMax, radius, outputfile)
   !MERCIER Wilfried - 02/09/18

      REAL(KIND = xp), EXTERNAL        :: fonction       !Function used
      REAL(KIND = xp), INTENT(IN)      :: TMax           !Maximum temperature at which we want to stop the search of the curve
      REAL(KIND = xp), INTENT(IN)      :: SMax           !Maximum surface density at which we want to stop the search of the curve
      REAL(KIND = xp), INTENT(IN)      :: TMin           !Minimum temperature at which we want to stop the search of the curve
      REAL(KIND = xp), INTENT(IN)      :: SMin           !Minimum surface density at which we want to stop the search of the curve
      REAL(KIND = xp), INTENT(IN)      :: Delta          !Step we use to find the first crossing of the curve
      REAL(KIND = xp), INTENT(IN)      :: radius         !Radius at which we are working
      INTEGER, INTENT(IN)              :: NMax           !Maximum number of iterations allowed for the dichotomy
      INTEGER, INTENT(IN)              :: NPointMax      !Maximum number of zeroes calculated if no boundary is reached
      INTEGER, INTENT(IN)              :: LogStep        !Boolean (1 = Logarithmic step, 0 = linear step)
      REAL(KIND = xp), INTENT(INOUT)   :: Eps            !Accuracy we want for the dichotomy
      REAL(KIND = xp), INTENT(INOUT)   :: DeltaS         !Step used to move (horizontically or vertically) to the next point on the sCurve
      REAL(KIND = xp), INTENT(INOUT)   :: S, T           !Initial position (bottom left corner) from which we start to look for the sCurve
      CHARACTER(LEN = *), intent(in)   :: outputfile     !Output file

      ! ... local...
      REAL(KIND = xp)                  :: SZero, TZero   !Position of the zeros
      REAL(KIND = xp)                  :: STemp, TTemp   !Temporary position
      REAL(KIND = xp)                  :: Slope          !Slope of the sCurve
      REAL(KIND = xp)                  :: TauEff         !Opacity at zeroes
      REAL(KIND = xp)                  :: TauEffTemp     !temporary opacity variable

      !Matrix representation (temporary)
      REAL(KIND = xp), DIMENSION(1)    :: SZeroMat, TZeroMat, DelQ, RAD

      INTEGER                          :: isOkFlag       !A flag to know if the dichotomy worked properly
      INTEGER                          :: i_refine       !Counter for the di4chotomy if it did not converge
      INTEGER                          :: N              !Number of iterations
      INTEGER                          :: dimi           !Indicates on which variable we should perform the dichotomy
      INTEGER                          :: unitOut        !Temporary unit used to write the zero positions in a file
      ! ----------------------------------------------------------------------------

      unitOut  = 100
      isOkFlag = 1
      N        = 0
      dimi     = 1
      Slope    = 1.0_xp

      STemp       = 0.0_xp
      TTemp       = 0.0_xp
      TauEff      = -1.0_xp
      TauEffTemp  = -1.0_xp
      RAD(1)      = radius

      WRITE(*,*)"LOG:", LogStep

      OPEN(NEWUNIT=unitOut, FILE=outputfile, STATUS='REPLACE')
      loop_buildsCurve : DO WHILE ( (T<TMax) .AND. (N<NPointMax) .AND. (isOkFlag==1) &
                                    .AND. (S<SMax) .AND. (T>TMin) .AND. (S>SMin) )

         !Perform the dichotomy
         Loop_refine_step : DO i_refine = 1,10

!            IF (ABS(TauEff-1)<0.1 ) THEN
!               WRITE(*,*)"I AM HERE"
!               CALL DichotomyPositiveOnly(fonction, S, T, TauEff, radius, Smin, TMin, SZero, &
!                    TZero, Eps, dimi, Delta, NMax, isOkFlag)
!            ELSE
               CALL Dichotomy(fonction, S, T, TauEff, radius, Smin, TMin, SZero, TZero, Eps, dimi, Delta, NMax, isOkFlag)
               cond_refine : IF ( (N>=NMax) .OR. (isOkflag==-2) ) THEN
                  WRITE(*,*)"REFINING for the ", i_refine, "th time."
                  cond_s : IF (dimi==1) THEN
                     S = STemp
                     CALL UpdateVar(S, DeltaS, REAL(2*i_refine,kind=xp), LogStep)
                  ELSE IF (dimi==0) THEN
                     T = TTemp
                     CALL UpdateVar(T, DeltaS, REAL(2*i_refine,kind=xp), LogStep)
                  END IF cond_s
               ELSE
                  EXIT
               END IF cond_refine
!            ENDIF
         END DO Loop_refine_step
         TauEffTemp = TauEff

         WRITE(*,*)"TAUEFF = ", TauEff

         !Check the dichotomy found a zero and that i as a good value
         IF ( isOkFlag==0 .OR. isOkFlag==-1 ) THEN
            WRITE(*,*)"Problem with dichotomy. Error flag", isOkFlag, " returned."
            RETURN
         ENDIF

         !If we reach the non-converging part we force the flag and the accuracy wanted
         IF ( ABS(TauEff-1)<0.1 ) THEN
!            WRITE(*,*)"proche de 1"
!            isOkFlag    = 1
!            Eps         = 1.0_xp
         ENDIF

         ! Store within a file the (SZero, TZero) coordinates
         SZeroMat(1) = SZero
         TZeroMat(1) = TZero
         IF (TZero-TTemp > 0) THEN
            N = N+1
            DelQ(1)     = fonction(SZero, TZero, TauEff, radius)
            WRITE(unitOut, '(E20.5, E20.5, E20.5, E20.5, E20.5, I2, I2)')surfaceDensitySI(SZeroMat, RAD), temperatureSI(TZeroMat), &
                            massicEnergySI(DelQ), delQ, TauEff, dimi, isOkFlag
            FLUSH(unitOut)

            WRITE(*,'(I5, E20.5, E20.5, I2, F20.5, F20.5, F20.5, E20.5, I2)')N, DelQ, massicEnergySI(DelQ), isOkFlag,&
            LOG10(surfaceDensitySI(SZeroMat, RAD)), LOG10(temperatureSI(TZeroMat)), Slope, DeltaS, dimi
         ENDIF

         !Compute the slope of the sCurve
         IF ( N>1) THEN
            IF (LogStep==0) THEN
               Slope = (TZero-TTemp)/(Szero-STemp)
            ELSE IF (LogStep==1) THEN
               Slope = (LOG10(TZero)-LOG10(TTemp))/(LOG10(SZero)-LOG10(STemp))
            ENDIF
         ENDIF

         !Find the dimension along which we must perform the next dichotomy given the slope
         IF ( (((ABS(Slope)>1) .OR. (TZero-TTemp<0)) .AND. (dimi/=0)) ) THEN
            dimi = 0
         ELSE IF ( (((ABS(Slope)<1) .OR. (TZero-TTemp<0)) .AND. (dimi/=1)) ) THEN
            dimi = 1
         ENDIF

         !If we reach the non-converging part we force the dichotomy axis
         IF ( ABS(TauEff-1)<0.1 ) THEN
!            WRITE(*,*)"CUCOUCOUUOUCOUCUOU"
!            dimi = 1
         ENDIF

         !Store previous position of the zero
         !If we went back in temperature (not possible for the SCurve) then go back to previous point
         !and change the dimension of the axis of the dichotomy
         IF (TZero-TTemp < 0) THEN
            WRITE(*,*)"WARNING : ", LOG10(surfaceDensitySI(SZeroMat, RAD)), LOG10(temperatureSI(TZeroMat)), dimi
            T     = TTemp
            S     = STemp
         ELSE
            T     = TZero
            S     = SZero
            TTemp = TZero
            STemp = SZero
         ENDIF

         !Move forward once a zero is found
         IF (dimi==1) THEN
            !Change the sign of the step if required
            IF ( (N>1) .AND. ((LogStep==0 .AND. Slope*DeltaS<0) .OR. &
               (LogStep==1 .AND. ((DeltaS>1 .AND. Slope<0) .OR. (DeltaS<1 .AND. Slope>0)))) ) THEN
               CALL InvertDelta(DeltaS, LogStep)
            ENDIF

            !If we reach the non-converging part we force it to go to the positive S
            IF ( ABS(TauEff-1)<0.1 ) THEN
!               DeltaS = ABS(DeltaS)
            ENDIF

            CALL UpdateVar(S, DeltaS, 1.0_xp, LogStep)
         ENDIF

         IF (dimi==0) THEN
            IF ( (N>1) .AND. ((LogStep==0 .AND. Slope*DeltaS<0) .OR. &
               (LogStep==1 .AND. (DeltaS<1 .AND. Slope>0))) ) THEN
               CALL InvertDelta(DeltaS, LogStep)
            ENDIF

            !If we reach the non-converging part we force it to go to the positive S
            IF ( ABS(TauEff-1)<0.1 ) THEN
!               DeltaS = ABS(DeltaS)
            ENDIF

            CALL UpdateVar(T, DeltaS, 1.0_xp, LogStep)
         ENDIF

      END DO loop_buildsCurve
      CLOSE(unitOut)

      IF (N>=NPointMax) THEN
         WRITE(*,*)"Maximum number of steps reached."
      ENDIF
      IF (T>=TMax) THEN
         WRITE(*,*)"The maximum temperature was reached."
      ENDIF
      IF (S>=SMax) THEN
         WRITE(*,*)"The maximum surface density was reached"
      ENDIF
      IF (T<=TMin) THEN
         WRITE(*,*)"The minimum temperature was reached."
      ENDIF
      IF (S<=SMin) THEN
         WRITE(*,*)"The minimum surface density was reached"
      ENDIF
      IF (isOkFlag .NE. 1) THEN
         WRITE(*,*)"Problem with dichotomy. Error flag", isOkFlag, " returned."
      ENDIF

   END SUBROUTINE sCurveFromDichotomy


!-------------------------------------------------------------------------------------------------------------------
   !> Inverts Delta given a linear or logarithmic step
   SUBROUTINE InvertDelta(Delta, LogStep)
      INTEGER, INTENT(IN)              :: LogStep  !1=log step, 0=linear step
      REAL(KIND = xp), INTENT(INOUT)   :: Delta    !Step to invert

      IF ( (LogStep==1) .AND. (Delta .NE. 0) ) THEN
         Delta = 1.0_xp/Delta
      ELSE IF (LogStep==0) THEN
         Delta = -Delta
      ELSE
         WRITE(*,*)"Could not invert the step as it is equal to 0."
      ENDIF
   END SUBROUTINE InvertDelta


!-------------------------------------------------------------------------------------------------------------------
   !> Update the value of a variable given a certain logarithmic/linear step
   SUBROUTINE UpdateVar(X, Delta, RefineCoeff, LogStep)
      INTEGER, INTENT(IN)              :: LogStep     !1=log step, 0=linear step
      REAL(KIND = xp), INTENT(IN)      :: Delta       !Step for the update
      REAL(KIND = xp), INTENT(IN)      :: RefineCoeff !Coefficient used to refine the step
      REAL(KIND = xp), INTENT(INOUT)   :: X           !Position coordinate to update

      IF (LogStep==1) THEN
         X = X*(Delta**(1.0_xp/RefineCoeff))
      ELSE
         X = X+(Delta/RefineCoeff)
      ENDIF
   END SUBROUTINE UpdateVar


!----------------------------------------------------------------------------------------------------------------------
   SUBROUTINE dichotomyFromGrid()

      USE io

      IMPLICIT NONE

      REAL(kind = xp), DIMENSION(:)      , ALLOCATABLE :: highTemp, lowTemp, highSigma, lowSigma
      LOGICAL        , DIMENSION(:,:)    , ALLOCATABLE :: QTransition
      REAL(kind = xp)                                  :: minTemperature,&
                                                          maxTemperature
      REAL(kind = xp)                                  :: minSigma, &
                                                          maxSigma
      REAL(kind = xp), DIMENSION(nSpaceStep)           :: initTemperature, initSigma, temp
      REAL(kind = xp), DIMENSION(:,:)    , ALLOCATABLE :: gridTemperature, &    !> temperature grid
                                                          gridSigma, &          !> surface density grid
                                                          gridQp, &             !> heating grid
                                                          gridQm, &             !> heat loss grid
                                                          gridTaueff           !> Opacity grid
      INTEGER                                          :: i, j, k, nS

      ! ...Dichotomy...
      REAL(KIND = xp)                                  :: Eps
      REAL(KIND = xp), DIMENSION(1)                    :: XZero, YZero, TauEff, DeltaQ, radius
      INTEGER                                          :: dim, isOkFlag, unitOut
      CHARACTER(LEN = 300)                             :: outputFile


      ALLOCATE(gridTemperature(nSigmaMesh,nTemperatureMesh))
      ALLOCATE(gridSigma(nSigmaMesh,nTemperatureMesh))
      ALLOCATE(gridQp(nSigmaMesh,nTemperatureMesh))
      ALLOCATE(gridQm(nSigmaMesh,nTemperatureMesh))
      ALLOCATE(gridTaueff(nSigmaMesh,nTemperatureMesh))

      !Init output
      IF (outputMesh) THEN
         CALL initOutputMeshFrame()
      END IF

      !Calculate init condition
      CALL getInitDisc(xRadius,initTemperature,initSigma,1.0_xp)

      !Initialize values for dichotomy
      Eps      = 1.0E-11_xp
      dim      = 0
      unitOut  = 11

      DO k = 1,nSpaceStep
         !Setup range in functions of initial condition
         minTemperature = initTemperature(k)*minScaleTemperatureMesh
         maxTemperature = initTemperature(k)*maxScaleTemperatureMesh
         minSigma = initSigma(k)*minScaleSigmaMesh
         maxSigma = initSigma(k)*maxScaleSigmaMesh

         DO i = 1,nSigmaMesh
            gridSigma(i,:) = 10._xp**(log10(minSigma) + REAL(i-1,xp)*(log10(maxSigma)-log10(minSigma))&
            /REAL(nSigmaMesh-1,xp))
         END DO
         DO j = 1,nTemperatureMesh
            gridTemperature(:,j) = 10._xp**(log10(minTemperature) + REAL(j-1,xp)*(log10(maxTemperature)-log10(minTemperature))&
            /REAL(nTemperatureMesh-1,xp))
         END DO

         !Calculate Grid
         DO j = 1,nTemperatureMesh
            CALL heatingTermsFromMeshGrid(gridTemperature(:,j), gridSigma(:,j), gridQp(:,j), gridQm(:,j),&
            xRadius(k), gridTaueff(:,j))
         END DO

         !Create Frame
         CALL determineMeshTransition(gridSigma, gridTemperature, nSigmaMesh, nTemperatureMesh, gridQp, gridQm,&
         QTransition,nS)
         ALLOCATE(highTemp(nS))
         ALLOCATE(lowTemp(nS))
         ALLOCATE(highSigma(nS))
         ALLOCATE(lowSigma(nS))
         CALL calculateSFrame(gridSigma, gridTemperature, nSigmaMesh, nTemperatureMesh, gridQp, gridQm,&
         QTransition,lowSigma,highSigma,lowTemp,highTemp)

         !>Dichotomy part
         !Creates new file for each radius
         WRITE(outputFile, '(a, I7, a)')"DichoFromMesh", k, ".dat"
         OPEN(NEWUNIT=unitOut, FILE=outputFile, STATUS='REPLACE')

         !For each range, perform a dichotomy to find the zero
         radius(1) = xradius(k)
         DO i = 1, nS
            CALL ComputeDicho(ComputeDelQ, lowSigma(i), highSigma(i), lowTemp(i), &
                 highTemp(i), XZero(1), YZero(1), Eps, dim, TauEff(1), xradius(k), isOkFlag)
            IF (isOkFlag==1) THEN
               DeltaQ(1) = ComputeDelQ(XZero(1), YZero(1), TauEff(1), xradius(k))
               WRITE(unitOut, '(E20.5, E20.5, E20.5)')surfaceDensitySI(XZero, radius), &
                                                      temperatureSI(YZero), massicEnergySI(DeltaQ)
            ENDIF
         END DO
         CLOSE(unitOut)

         !Output mesh
         IF (outputMesh) THEN
            CALL writeMesh(k, nTemperatureMesh, nSigmaMesh, gridTemperature, gridSigma, gridQp, gridQm, gridTaueff)
            CALL writeFrame(k,lowSigma,highSigma,lowTemp,highTemp)
         END IF

         !Deallocate memory use for frame
         DEALLOCATE(lowSigma)
         DEALLOCATE(highSigma)
         DEALLOCATE(lowTemp)
         DEALLOCATE(highTemp)
      END DO

      !Init output
      IF (outputMesh) THEN
         CALL closeOutputMeshFrame()
      END IF

      DEALLOCATE(gridTemperature)
      DEALLOCATE(gridSigma)
      DEALLOCATE(gridQp)
      DEALLOCATE(gridQm)
      DEALLOCATE(gridTaueff)

   END SUBROUTINE dichotomyFromGrid


   !------------------------------------------------------------------------------------------------------------------------------------
      !>Subroutine to create and calculate the meshgrid for s curve
   SUBROUTINE meshGridScurve(minTemperature, maxTemperature, nTemperature, minSigma, maxSigma, nSigma, x, &
                            temperature, sigma, Qp, Qm, Tau_eff)

      IMPLICIT NONE
      REAL(kind = xp),                                  INTENT(IN)  :: minTemperature, & !> minimum temperature of the grid
                                                                       maxTemperature    !> maximum temperature of the grid
      INTEGER        ,                                  INTENT(IN)  :: nTemperature      !> number of step for the temperature
      REAL(kind = xp),                                  INTENT(IN)  :: minSigma, &       !> minimum surface density
                                                                       maxSigma          !> maximum surface density
      INTEGER        ,                                  INTENT(IN)  :: nSigma            !> number of step for surface density
      REAL(kind = xp),                                  INTENT(IN)  :: x                 !> radius use for the calculation
      REAL(kind = xp), DIMENSION(nSigma, nTemperature), INTENT(OUT) :: temperature, &    !> temperature grid
                                                                       Tau_eff, &        !> Opacity grid
                                                                       sigma, &          !> surface density grid
                                                                       Qp, &             !> heating grid
                                                                       Qm                !> heat loss grid
      REAL(kind = xp), DIMENSION(nSigma, nTemperature)              :: xMatrix
      REAL(kind = xp), DIMENSION(nSigma*nTemperature)               :: temperatureVec, sigmaVec, QpVec, QmVec
      INTEGER                                                       :: i, j

      xMatrix = x
      !Create grid for surface density and temperature
      DO i = 1,nSigma
         sigma(i,:) = 10._xp**(log10(minSigma) + REAL(i-1,xp)*(log10(maxSigma)-log10(minSigma))&
         /REAL(nSigma-1,xp))
      END DO
      DO j = 1,nTemperature
         temperature(:,j) = 10._xp**(log10(minTemperature) + REAL(j-1,xp)*(log10(maxTemperature)-log10(minTemperature))&
         /REAL(nTemperature-1,xp))
      END DO

      !Normalize
      DO j = 1,nTemperature
         !xMatrix(:,j) = radiusNormalize(rMatrix(:,j))
         temperature(:,j) = temperatureNormalize(temperature(:,j))
         sigma(:,j) = surfaceDensityNormalize(sigma(:,j),xMatrix(:,1))
      END DO

      !Calculate grid
      DO j = 1,nTemperature
         CALL heatingTermsFromMeshGrid(temperature(:,j), sigma(:,j), Qp(:,j), Qm(:,j), xMatrix(1,1), Tau_eff(:,j))
      END DO

      !Denormalize
      DO j = 1,nTemperature
         temperature(:,j) = temperatureSI(temperature(:,j))
         sigma(:,j) = surfaceDensitySI(sigma(:,j), xMatrix(:,j))
         Qp(:,j) = massicEnergySI(Qp(:,j))
         Qm(:,j) = massicEnergySI(Qm(:,j))
      END DO

   END SUBROUTINE meshGridScurve


   !---------------------------------------------------------------------------------------------------------------------------
      !>Retourne le point de pente max pour chaque rayon
   SUBROUTINE WriteMaxSlopePoint(Eps, Delta, DeltaS, LogStep, NMax)

      REAL(KIND = xp), INTENT(IN)            :: Eps      !Accuracy for the dichotomy
      REAL(KIND = xp), INTENT(IN)            :: Delta    !Step for the dichotomy
      REAL(KIND = xp), INTENT(INOUT)         :: DeltaS   !Step used to find the first crossing
      INTEGER, INTENT(IN)                    :: LogStep  !0=Linear, 1=Logarithmic step
      INTEGER, INTENT(IN)                    :: NMax     !Max number of iterations for the dichotomy

      ! ...local variables...
      REAL(KIND = xp)                        :: TauEff
      REAL(KIND = xp), DIMENSION(nSpaceStep) :: initTemperature, initSigma, temp
      REAL(KIND = xp), DIMENSION(1)          :: S, T, SMin, TMin, radius
      INTEGER                                :: i, unitOut

      !Initialize x, T, S values
      CALL getInitDisc(xRadius,initTemperature,initSigma,1.0_xp)

      unitOut = 11
      OPEN(NEWUNIT=unitOut, FILE="MaxSlopePoint", STATUS='REPLACE')
      DO i = 1, nSpaceStep
         radius(1)   = xradius(i)
         SMin(1)     = initSigma(i)*minScaleSigmaMesh
         TMin(1)     = initTemperature(i)*minScaleTemperatureMesh
         CALL findInifniteSlopePoint(ComputeDelQ, SMin(1), TMin(1), &
              S(1), T(1), Eps, Delta, DeltaS, LogStep, NMax, TauEff, xRadius(i))

         WRITE(unitOut, '(I6, E20.5, E20.5, E20.5, E20.5, E20.5)')i, radiusSI(radius), S, &
                         surfaceDensitySI(S, radius), T, temperatureSI(T)
      ENDDO
      CLOSE(unitOut)

   END SUBROUTINE WriteMaxSlopePoint


   !---------------------------------------------------------------------------------------------------------------------------
      !>Build all sCurves
   SUBROUTINE BuildAllSCurves(Eps, Delta, DeltaS, LogStep, NPointMax, NMax)

      INTEGER, INTENT(IN)                    :: LogStep     !0=Linear, 1=Logarithmic step
      INTEGER, INTENT(IN)                    :: NMax        !Max number of iterations for the dichotomy
      INTEGER, INTENT(IN)                    :: NPointMax   !Max number of iterations for the dichotomy
      REAL(KIND = xp), INTENT(IN)            :: Delta       !Step for the dichotomy
      REAL(KIND = xp), INTENT(INOUT)         :: DeltaS      !Step used to find the first crossing
      REAL(KIND = xp), INTENT(INOUT)         :: Eps         !Accuracy for the dichotomy

      ! ...local variables...
      REAL(KIND = xp)                        :: TauEff
      REAL(KIND = xp), DIMENSION(nSpaceStep) :: initTemperature, initSigma
      REAL(KIND = xp), DIMENSION(1)          :: S, T, SMin, TMin, SMax, TMax
      REAL(KIND = xp), DIMENSION(1)          :: radius
      INTEGER                                :: i
      CHARACTER(LEN = 100)                     :: filename

      !Initialize x, T, S values
      CALL getInitDisc(xRadius,initTemperature,initSigma,1.0_xp)

      DO i = 35, 35 !nSpaceStep
         WRITE(filename, '(a, I6.6, a)')"SCurveFromFullDicho", i, ".dat"

         radius(1)   = xradius(i)
         SMin(1)     = initSigma(i)*minScaleSigmaMesh
         TMin(1)     = initTemperature(i)*minScaleTemperatureMesh

         SMin(1)     = 1.0E3_xp
         SMin        = SurfaceDensityNormalize(SMin, radius)
         TMin(1)     = 1.0E6_xp
         TMin        = temperatureNormalize(TMin)

         SMax(1)     = initSigma(i)*maxScaleSigmaMesh
         TMax(1)     = initTemperature(i)*maxScaleTemperatureMesh

         WRITE(*,*)surfaceDensitySI(SMin, radius), temperatureSI(TMin), radius 
         CALL sCurveFromDichotomy(ComputeDelQ, SMin(1), TMin(1), 1.0E10_xp, 1.0E10_xp, 0.0_xp, 0.0_xp, &
                                  Eps, Delta, DeltaS, LogStep, NPointMax, NMax, radius(1), filename)
      ENDDO

   END SUBROUTINE BuildAllSCurves


   !---------------------------------------------------------------------------------------------------------------------------
      !>Calcule les matrices Q+, Q- (et Qadv) à partir de T et S
   SUBROUTINE heatingTermsFromMeshGrid(Tgrid, Sgrid, Qp, Qm, x, Tau_eff)

      REAL(KIND=xp), INTENT(IN), DIMENSION(:)   :: Tgrid
      REAL(KIND=xp), INTENT(IN), DIMENSION(:)   :: Sgrid
      REAL(KIND=xp), INTENT(IN)                 :: x
      REAL(KIND=xp), INTENT(OUT), DIMENSION(:)  :: Qp
      REAL(KIND=xp), INTENT(OUT), DIMENSION(:)  :: Qm
      REAL(KIND=xp), INTENT(OUT), DIMENSION(:)  :: Tau_eff

      ! .. local ..
      INTEGER                                   :: n_size, i                        !Size of vectors and dummy variable
      REAL(KIND=xp)                             :: omega                            !Angular speed (scalar)
      REAL(KIND=xp), ALLOCATABLE, DIMENSION(:)  :: xx                               !Position vector (with same value everywhere)
      REAL(KIND=xp), ALLOCATABLE, DIMENSION(:)  :: xxx                              !Angular speed vector
      REAL(KIND=xp), ALLOCATABLE, DIMENSION(:)  :: H
      REAL(KIND=xp), ALLOCATABLE, DIMENSION(:)  :: cs
      REAL(KIND=xp), ALLOCATABLE, DIMENSION(:)  :: nu
      REAL(KIND=xp), ALLOCATABLE, DIMENSION(:)  :: masse_volumique
      REAL(KIND=xp), ALLOCATABLE, DIMENSION(:)  :: Pgaz
      REAL(KIND=xp), ALLOCATABLE, DIMENSION(:)  :: Prad
      REAL(KIND=xp), ALLOCATABLE, DIMENSION(:)  :: P
      REAL(KIND=xp), ALLOCATABLE, DIMENSION(:)  :: kappaff
      REAL(KIND=xp), ALLOCATABLE, DIMENSION(:)  :: bet
      REAL(KIND=xp), ALLOCATABLE, DIMENSION(:)  :: Cv
      REAL(KIND=xp), ALLOCATABLE, DIMENSION(:)  :: Fz
      REAL(KIND=xp), ALLOCATABLE, DIMENSION(:)  :: epsilonff
      !--------------------------------------------------------------

      !Allocating vectors
      n_size = size(Tgrid)
      ALLOCATE(H(n_size))
      ALLOCATE(cs(n_size))
      ALLOCATE(nu(n_size))
      ALLOCATE(masse_volumique(n_size))
      ALLOCATE(Pgaz(n_size))
      ALLOCATE(kappaff(n_size))
      ALLOCATE(Prad(n_size))
      ALLOCATE(P(n_size))
      ALLOCATE(bet(n_size))
      ALLOCATE(Cv(n_size))
      ALLOCATE(Fz(n_size))
      ALLOCATE(epsilonff(n_size))
      ALLOCATE(xx(n_size))
      ALLOCATE(xxx(n_size))

      DO i=1,n_size
         xx(i)    = x
         xxx(i)   = x
      END DO

      xxx = getAngularSpeed(xxx)
      omega = xxx(1)

      !Détermination de la hauteur du disque
      H = getHalfLength(xx,xxx,Tgrid,Sgrid)

      !print *, 'H =', H
      !Détermination de la vitesse du son
      cs = getSoundSpeed(H,xxx)

      !Détermination de la masse volumique
      masse_volumique = getRho(xx,Sgrid,H)

      !Détermination de la viscosité
      nu = getViscosity(cs,H)
      !Détermination du terme de chauffage
      Qp = getQplus(nu,xxx)

      !Détermination de la pression due au gaz
      Pgaz = getPressureGas(masse_volumique,Tgrid)

      !Détermination de l'opacité free-free
      kappaff = getOpacityFree(Tgrid,masse_volumique)

      !Détermination de la pression due à la radiation
      Prad = getPressureRad(Tgrid)

      !Détermination de la pression totale
      P = getPressureTotal(Pgaz,Prad)

      !Détermination du rapport de la pression du gaz à la pression totale
      bet = getBeta(P,Pgaz)

      !Détermination de la capacité calorifique à volume constant
      Cv = getHeatCapacity(bet)

      !Détermination de l'émissivité free-free
      epsilonff = getEmissivityFree(masse_volumique,Tgrid)

      !Détermination de la profondeur optique
      Tau_eff = getTaueff(kappaff,Sgrid,xx)

      !Détermination du flux radiatif
      Fz = getRadiatifFlux(Tgrid,xx,kappaff,Sgrid,H,Tau_eff,epsilonff)

      !Détermination du terme de refroidissement
      Qm = getQmoins(xx,Fz,Sgrid)

   END SUBROUTINE heatingTermsFromMeshGrid


   !>Calulate transitions in the mesh and their number
   SUBROUTINE determineMeshTransition(SGrid,TGrid,nSigma,nTemp,Qp,Qm,QTransition,nS)

      IMPLICIT NONE

      REAL(kind=xp), DIMENSION(:,:)              , INTENT(IN)  :: SGrid,TGrid,Qp,Qm
      INTEGER                                    , INTENT(IN)  :: nSigma,nTemp
      LOGICAL      , DIMENSION(nSigma-1,nTemp-1) , INTENT(OUT) :: QTransition
      INTEGER                                    , INTENT(OUT) :: nS
      REAL(kind=xp), DIMENSION(nSigma,nTemp)                   :: deltaQ
      INTEGER                                                  :: i,j
      nS = 0

      !Determine transition in Qp-Qm and the number of transition
      deltaQ = Qp-Qm
      DO j=1,(nTemp-1)
         DO i=1,(nSigma-1)
            QTransition(i,j) = (((deltaQ(i,j+1)*deltaQ(i,j))<0._xp) .OR. &
            ((deltaQ(i+1,j)*deltaQ(i,j))<0._xp))
            IF (QTransition(i,j)) THEN
               nS = nS+1
            END IF
         END DO
      END DO

   END SUBROUTINE determineMeshTransition


   !>Calulate value to frame the S curve
   SUBROUTINE calculateSFrame(SGrid,TGrid,nSigma,nTemp,Qp,Qm,QTransition,lowSigma,highSigma,lowTemp,highTemp)

      IMPLICIT NONE

      REAL(kind=xp), DIMENSION(:,:)           , INTENT(IN)     :: SGrid, TGrid,Qp,Qm
      INTEGER                                 , INTENT(IN)     :: nSigma,nTemp
      LOGICAL      , DIMENSION(nSigma-1,nTemp-1) , INTENT(IN)  :: QTransition
      REAL(kind=xp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)  :: highSigma, lowSigma, highTemp, lowTemp
      INTEGER                                                  :: i,j,k


      k=1
      DO j=1,(nTemp-1)
         DO i=1,(nSigma-1)
            IF (QTransition(i,j)) THEN
               lowTemp(k)   = TGrid(i,j)
               highTemp(k)  = TGrid(i+1,j+1)
               lowSigma(k)  = SGrid(i,j)
               highSigma(k) = SGrid(i+1,j+1)
               k = k+1
            END IF
         END DO
      END DO

   END SUBROUTINE calculateSFrame


   !----------------------------------------------------------------------------------------------------------------------------------------------------
      !> Look for the zeros up to the DeltaT=infinity point and return it
   SUBROUTINE findInifniteSlopePoint(fonction, S, T, SOut, TOut, Eps, Delta, DeltaS, LogStep, NMax, TauEff, radius)
   !MERCIER Wilfried - 02/09/18

      REAL(KIND = xp), EXTERNAL        :: fonction       !Function used
      REAL(KIND = xp), INTENT(IN)      :: Eps            !Accuracy we want for the dichotomy
      REAL(KIND = xp), INTENT(IN)      :: Delta          !Step we use to find the first crossing of the curve
      REAL(KIND = xp), INTENT(IN)      :: radius         !Radius at which we are working
      INTEGER, INTENT(IN)              :: NMax           !Maximum number of iterations allowed for the dichotomy
      INTEGER, INTENT(IN)              :: LogStep        !Boolean (1 = Logarithmic step, 0 = linear step)
      REAL(KIND = xp), INTENT(INOUT)   :: DeltaS         !Step used to move (horizontically or vertically) to the next point on the sCurve
      REAL(KIND = xp), INTENT(INOUT)   :: S, T           !Initial position (bottom left corner) from which we start looking for the sCurve

      ! ... local...
      REAL(KIND = xp)                  :: SZero, TZero   !Position of the zeros
      REAL(KIND = xp)                  :: STemp, TTemp   !Temporary position
      REAL(KIND = xp)                  :: Slope          !Slope of the sCurve
      REAL(KIND = xp)                  :: TauEff         !Opacity at zeroes
      REAL(KIND = xp)                  :: SMin, TMin     !Minimum values
      REAL(KIND = xp)                  :: SOut, TOut     !Output values

      INTEGER                          :: isOkFlag       !A flag to know if the dichotomy worked properly
      INTEGER                          :: i_refine       !Counter for the dichotomy if it did not converge
      INTEGER                          :: N              !Number of iterations
      INTEGER                          :: dimi           !Indicates on which variable we should perform the dichotomy
      ! ----------------------------------------------------------------------------

      isOkFlag    = 1
      N           = 0
      dimi        = 1
      Slope       = -11.0_xp
      SMin        = 0.0_xp
      TMin        = 0.0_xp

      STemp       = 0.0_xp
      TTemp       = 0.0_xp
      TauEff      = -1.0_xp

      loop_buildsCurve : DO WHILE ( (isOkFlag==1) .AND. (T>TMin) .AND. (S>SMin) .AND. (DeltaS>0) )

         !Perform the dichotomy
         Loop_refine_step : DO i_refine = 1,10
            CALL Dichotomy(fonction, S, T, TauEff, radius, Smin, TMin, SZero, TZero, Eps, dimi, Delta, NMax, isOkFlag)
            cond_refine : IF ( (N>=NMax) .OR. (isOkflag==-2) ) THEN
               cond_s : IF (dimi==1) THEN
                  S = STemp
                  CALL UpdateVar(S, DeltaS, REAL(2*i_refine,kind=xp), LogStep)
               ELSE IF (dimi==0) THEN
                  T = TTemp
                  CALL UpdateVar(T, DeltaS, REAL(2*i_refine,kind=xp), LogStep)
               END IF cond_s
            ELSE
               EXIT
            END IF cond_refine
         END DO Loop_refine_step

         !Check the dichotomy found a zero and that i as a good value
         IF ( isOkFlag==0 .OR. isOkFlag==-1 ) THEN
            WRITE(*,*)"Problem with dichotomy. Error flag", isOkFlag, " returned."
            RETURN
         ENDIF

         N = N+1

         !Compute the slope of the sCurve
         IF ( N>1) THEN
            IF (LogStep==0) THEN
               Slope = (TZero-TTemp)/(Szero-STemp)
            ELSE IF (LogStep==1) THEN
               Slope = (LOG10(TZero)-LOG10(TTemp))/(LOG10(SZero)-LOG10(STemp))
            ENDIF
         ENDIF

!         WRITE(*,*)N, Slope, TZero, SZero
         IF ( (Slope<0) .AND. (N>1) ) THEN
            SOut = STemp
            TOut = TTemp
            RETURN
         ENDIF

         !Find the dimension along which we must perform the next dichotomy given the slope
         IF (N>1) THEN
            IF ( (((ABS(Slope)>1) .OR. (TZero-TTemp<0)) .AND. (dimi/=0)) ) THEN
               dimi = 0
            ELSE IF ( (((ABS(Slope)<1) .OR. (TZero-TTemp<0)) .AND. (dimi/=1)) ) THEN
               dimi = 1
            ENDIF
         ENDIF

         !Store previous position of the zero
         T     = TZero
         S     = SZero
         TTemp = TZero
         STemp = SZero

         !Move forward once a zero is found
         IF (dimi==1) THEN
            !Change the sign of the step if required
            IF ( (N>1) .AND. ((LogStep==0 .AND. Slope*DeltaS<0) .OR. &
               (LogStep==1 .AND. ((DeltaS>1 .AND. Slope<0) .OR. (DeltaS<1 .AND. Slope>0)))) ) THEN
               CALL InvertDelta(DeltaS, LogStep)
            ENDIF
            CALL UpdateVar(S, DeltaS, 1.0_xp, LogStep)
         ENDIF

         IF (dimi==0) THEN
            IF ( (N>1) .AND. ((LogStep==0 .AND. Slope*DeltaS<0) .OR. &
               (LogStep==1 .AND. (DeltaS<1 .AND. Slope>0))) ) THEN
               CALL InvertDelta(DeltaS, LogStep)
            ENDIF
            CALL UpdateVar(T, DeltaS, 1.0_xp, LogStep)
         ENDIF

      END DO loop_buildsCurve

   END SUBROUTINE findInifniteSlopePoint


   !-------------------------------------------------------------------------------------------------------------------
      !>Qplus-Qminus at radius
   REAL(KIND = xp) FUNCTION ComputeDelQ(S, T, TauEff, rad)

      REAL(KIND = xp), INTENT(IN)   :: S, T           !Temperature and Surface density
      REAL(KIND = xp), INTENT(IN)   :: rad         !Radius
      REAL(KIND = xp), INTENT(OUT)  :: TauEff         !Opacity

      !...local...
      REAL(KIND = xp), DIMENSION(1) :: Qp, Qm, DelQ   !Heating terms
      REAL(KIND = xp), DIMENSION(1) :: SMat, TMat     !Matrix form of S and T
      REAL(KIND = xp), DIMENSION(1) :: TauEffMat
      !-----------------------------------------------------------------------------------------

      SMat(1)  = S
      TMat(1)  = T

      CALL heatingTermsFromMeshGrid(TMat, SMat, Qp, Qm, rad, TauEffMat)
      TauEff            = TauEffMat(1)
      DelQ              = Qp-Qm
      ComputeDelQ       = DelQ(1)

   END FUNCTION ComputeDelQ

END MODULE sCurve
