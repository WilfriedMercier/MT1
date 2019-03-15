MODULE constantPhysic

   USE constantProgram

   REAL(KIND = xp), PARAMETER, PUBLIC :: GConst = 6.67408e-11_xp            !>Gravitational constant -> m^3/(kg)/(s^2)
   REAL(KIND = xp), PARAMETER, PUBLIC :: cConst = 299792458_xp              !>Speed of light -> m/s
   REAL(KIND = xp), PARAMETER, PUBLIC :: kbConst = 1.38064852e-23_xp        !>Boltzmann constant -> J/K
   REAL(KIND = xp), PARAMETER, PUBLIC :: hConst = 6.626070040e-34_xp        !>Planck constant -> J.s
   REAL(KIND = xp), PARAMETER, PUBLIC :: sigmaConst = 5.670367e-8_xp        !>Stefan -> W/(m^2)/(K^4)
   REAL(KIND = xp), PARAMETER, PUBLIC :: mpConst = 1.672621898e-27_xp       !>Proton mass -> kg
   REAL(KIND = xp), PARAMETER, PUBLIC :: muConst = 0.617283950617284_xp     !>Sans dimension

   REAL(KIND = xp), PARAMETER, PUBLIC :: aConst = 7.564e-16_xp              !>Radiation constant
   REAL(KIND = xp), PARAMETER, PUBLIC :: thomsonOpacityConst = 0.02_xp      !>Constant in thomson opacity
   REAL(KIND = xp), PARAMETER, PUBLIC :: opacityFreeConst = 6.13e18_xp      !>Constant in opacity freefree
   REAL(KIND = xp), PARAMETER, PUBLIC :: emissivityFreeConst = 6.22e13_xp   !>Constant in emossivity free free

   REAL(kind = xp), PARAMETER, PUBLIC :: MSunSI = 1.98847e30_xp             !> Sun mass -> kg

   REAL(KIND = xp), PARAMETER, PUBLIC :: piConst = 3.141592653589793238462643_xp !>Pi
   REAL(KIND = xp), PARAMETER, PUBLIC :: RConst = kbConst/mpConst !> m^2.s^(-1)  (m^2/(K.s^2) perfect gas constant

END MODULE constantPhysic
