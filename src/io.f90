MODULE io

   USE constantProgram
   IMPLICIT NONE

   PRIVATE
   CHARACTER(len = 300) :: simulationName, path
   INTEGER              :: snapshotNum, unitOutMesh, unitOutFrame

   LOGICAL, PUBLIC      :: newSim

   PUBLIC readInputFile, initOutput, closeOutput, writeState, move2Dir,&
   initOutputMeshFrame, closeOutputMeshFrame, writeMesh, writeFrame

   CONTAINS

      !Subroutine for reading input file identified by filename
      SUBROUTINE readInputFile(filename)

         USE constantSimulation
         IMPLICIT NONE
         CHARACTER(len = *), INTENT(IN) :: filename
         NAMELIST /physicalsetting/ simulationName,path,blackHoleMass,accretionRate0,rmax,tmax,itmax, &
         nSpaceStep,alpha,chemX,chemY,chemZ,&
         variationAccretionRateInit,minScaleTemperatureMesh, &
         maxScaleTemperatureMesh,minScaleSigmaMesh,maxScaleSigmaMesh,&
         nTemperatureMesh,nSigmaMesh,outputMesh,snapshotFrequency
         !----------------------------------------


         OPEN (UNIT=10, File=filename , STATUS='OLD')
         READ (10, NML=physicalsetting)
         CLOSE (10)

         blackHoleMassSI = blackHoleMass
         tmaxSI = tmax

         newSim = .TRUE.

      END SUBROUTINE readInputFile

      !>Subroutine which opens a file and create a directory
      SUBROUTINE initOutput()

         IMPLICIT NONE

         CALL move2Dir()
         snapshotNum = 1

      END SUBROUTINE initOutput


      !> Used to move to the output directory, use it directly ONLY if you don't use initOutput
      SUBROUTINE move2Dir()

         IMPLICIT NONE

         !Create an output directory and move there
         CALL SYSTEM('mkdir -p ' // TRIM(path) // '/' // TRIM(simulationName))
         CALL CHDIR(TRIM(path) // '/' // TRIM(simulationName))

      END SUBROUTINE move2Dir


      !>Subroutine for opening file and create directory
      SUBROUTINE closeOutput()

         IMPLICIT NONE

      END SUBROUTINE closeOutput

      !>Subroutine for writing a statimee of code
      SUBROUTINE writeState(time,r,H,Temperature,P,Pgaz,Prad,beta,sigma,cs,nu,v,accretionRate, &
                            Qp,Qm,Qadv,Cv,Fz,kff,ke,epsilonff,tauff,omega)

         USE constantSimulation

         IMPLICIT NONE

         REAL (kind=xp),                INTENT(in)  :: time    !>time
         REAL (kind=xp), DIMENSION (:), INTENT(in)  :: r, &    !>Radius
                                                       H, &    !>heigt
                                                       Temperature, &    !>temperature
                                                       P, &    !>pressure
                                                       Pgaz, & !>pressure gaz
                                                       Prad, & !>pressure radiation
                                                       beta, & !>Pressure indicator
                                                       sigma, & !>matter surface density
                                                       cs, &
                                                       nu, &
                                                       v, &
                                                       accretionRate, &
                                                       Qp, &
                                                       Qm, &
                                                       Qadv, &
                                                       Cv, &
                                                       Fz, &
                                                       kff, &
                                                       ke, &
                                                       epsilonff, &
                                                       tauff, &
                                                       omega

         INTEGER                                    :: unitOut, i

         CHARACTER(len=40)                          :: filename !> time on string format

         !Determine filename
         WRITE (filename, "(a,i6.6,a)")  'snapshot_', snapshotNum, '.dat'

         !>Create output file and write variable
         OPEN(NEWUNIT=unitOut, FILE= filename, STATUS='REPLACE')
         WRITE (unitOut, *) '# time = ', time

         WRITE (unitOut, 103) 'r         ', &
                              'H         ', &
                              'T         ', &
                              'P         ', &
                              'Pgaz         ', &
                              'Prad         ', &
                              'beta         ', &
                              'sigma         ', &
                              'cs         ', &
                              'nu         ', &
                              'v         ', &
                              'accretionRate         ', &
                              'Qp         ', &
                              'Qm         ', &
                              'Qadv         ', &
                              'Cv         ', &
                              'Fz         ', &
                              'kff         ', &
                              'ke         ', &
                              'epsilonff         ', &
                              'tauff         ', &
                              'omega         '

         DO i = 1,nSpaceStep
            WRITE (unitOut, 103) r(i), &
                                 H(i), &
                                 Temperature(i), &
                                 P(i), &
                                 Pgaz(i), &
                                 Prad(i), &
                                 beta(i), &
                                 sigma(i), &
                                 cs(i), &
                                 nu(i), &
                                 v(i), &
                                 accretionRate(i), &
                                 Qp(i), &
                                 Qm(i), &
                                 Qadv(i), &
                                 Cv(i), &
                                 Fz(i), &
                                 kff(i), &
                                 ke(i), &
                                 epsilonff(i), &
                                 tauff(i), &
                                 omega(i)
         END DO
         CLOSE(unitOut)

         !Define the format of the output
         103 FORMAT (22(2x,1pg18.10e3))

         !Increase snapshot number
         snapshotNum = snapshotNum + 1

      END SUBROUTINE writeState

      SUBROUTINE initOutputMeshFrame()

         IMPLICIT NONE

         OPEN(NEWUNIT=unitOutMesh, FILE = 'ScurveMesh.dat', STATUS='REPLACE')
         OPEN(NEWUNIT=unitOutFrame, FILE = 'ScurveFrame.dat', STATUS='REPLACE')

      END SUBROUTINE initOutputMeshFrame

      SUBROUTINE closeOutputMeshFrame()

         IMPLICIT NONE

         CLOSE(unitOutMesh)
         CLOSE(unitOutFrame)

      END SUBROUTINE closeOutputMeshFrame

      !>Write mesh of the scurve
      SUBROUTINE writeMesh(iRadius, nTemperature, nSigma, temperature, sigma, Qp, Qm, Taueff)

         USE constantSimulation

         IMPLICIT NONE

         INTEGER, INTENT(IN)                         :: iRadius, nTemperature, nSigma
         REAL(kind = xp), DIMENSION(:,:), INTENT(IN) :: temperature, sigma, Qp, Qm, Taueff
         INTEGER                                     :: j, k

         DO k = 1,nTemperature
            DO j = 1,nSigma
               WRITE (unitOutMesh, 101) iRadius, j, k, rRadius(iRadius), temperature(j,k), sigma(j,k), Qp(j,k), Qm(j,k), TauEff(j,k)
            END DO
         END DO

         !Define the format of the output
         101 FORMAT (3(2x,I10),6(2x,1pg18.10e3))

      END SUBROUTINE writeMesh

      !>Write frame of the s curve
      SUBROUTINE writeFrame(iRadius,lowSigma,highSigma,lowTemp,highTemp)

         IMPLICIT NONE

         INTEGER, INTENT(IN)                       :: iRadius
         REAL(kind = xp), DIMENSION(:), INTENT(IN) :: lowSigma,highSigma,lowTemp,highTemp
         INTEGER                                   :: j

         DO j = 1,size(lowSigma)
            WRITE (unitOutFrame, 102) iRadius, j, lowSigma(j), highSigma(j), lowTemp(j), highTemp(j)
         END DO

         !Define the format of the output
         102 FORMAT (2(2x,I10),4(2x,1pg18.10e3))
      END SUBROUTINE writeFrame
END MODULE io
