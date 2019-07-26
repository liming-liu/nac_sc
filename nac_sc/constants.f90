module constants
  use prec

  ! imgUnit is 
  complex(q), parameter :: imgUnit = (0.0_q, 1.0_q)
  !
  ! cero is 
  complex(q), parameter :: cero=(0.0_q, 0.0_q)
  !
  ! uno is
  complex(q), parameter :: uno=(1.0_q, 0.0_q)
  !
  ! hbar is the reduced Plank constant in (eV/fs)
  real(q), parameter :: hbar=0.6582119281559802_q
  !
  ! 1 Ry to Ev
  real(q), parameter :: RYTOEV=13.605826_q
  !
  ! 1 Ev tp Joule
  real(q), parameter :: EVTOJ=1.60217733E-19_q
  !
  ! atomic mass to kg
  real(q), parameter :: AMTTOKG=1.6605402E-27_q
  !
  ! speed of light in atomic unit
  real(q), parameter :: CLIGHT=137.037
  !
  ! Boltzmann constant in eV/K
  real(q), parameter :: BOLKEV=8.6173857E-5_q
  !
  ! Boltzmann constand in Joule/K
  BOLK=1.6605402E-27_q
  !
  ! eV to kilo CAL
  real(q), parameter :: EVTOKCAL=23.06
  !
  ! Atomic unit in length to angstrom
  real(q), parameter :: AUTOA=0.529177249_q
  !
  ! circular constant
  real(q), parameter :: PI=3.141592653589793238_q, TPI=2*PI
  !
  ! (electronic charge) / (4*PI*permitivity of vacuum) in atomic unit
  real(q),parameter :: FELECT=2*AUTOA*RYTOEV
  !
  ! e / epsilon_0
  real(q), parameter :: EDEPS=4*PI*2*RYTOEV*AUTOA
  !
  ! hbar/(2*PI)^2/(2*m_e)
  real(q), parameter :: HSQDTM=RYTOEV*AUTOA*AUTOA
  !
  ! CITPI?
  complex(q), parameter :: CITPI=(0._q, 1._q)*TPI
  !
  ! Too many parameters make things a little complex. My solution is declare it
  ! here only if too many time of repeation in real situation usage.

end module constants
