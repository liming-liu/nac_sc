module datatype

    type psi
        integer :: iband, ikpts, ispin
    end type psi

    type waveinfo
        character(len=255) :: WAVECAR
        real*8             :: ENCUT
        integer            :: ISPIN, NKPTS, NBANDS, qw, IU, MAXPLWS
        integer, allocatable, dimension(:)     :: NPLWS
        real*8,  allocatable, dimension(:,:)   :: VKPTS
        real*8,  allocatable, dimension(:,:,:) :: BANDS
    end type waveinfo

    type overlap
        integer :: NBANDS
        integer :: TSTEPS
        real*8  :: dt
        real*8, allocatable, dimension(:,:,:) :: Dij
        real*8, allocatable, dimension(:,:)   :: Eig
    end type overlap
    
    real(kind=8), parameter :: c=299792458              ! light of speed.

end module datatype
