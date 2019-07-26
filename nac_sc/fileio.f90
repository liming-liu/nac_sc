module fileio
  use prec
  implicit none
  type namdInfo
      integer :: BMIN        ! Lower basis No.
      integer :: BMAX        ! Higher basis No.
      integer :: NBASIS      ! No. of adiabatic orbitals
      integer :: NBANDS      ! No. of total bands
      integer :: INIBAND     ! initial band for place e/h
      integer :: NSW         ! No. of MD steps
      integer :: NAMDTIME    ! No. of NAMD steps
      integer :: NAMDTINI    ! Initial time step of NAMD
      integer :: NTRAJ       ! No. of surface hopping tranectories
      integer :: NELM        ! No. of steps of electron wave propagation
      integer :: NSAMPLE     ! No. of samples for namd
      integer, allocatable, dimension(:) :: NAMDTINI_A  
      integer, allocatable, dimension(:) :: INIBAND_A
      real(kind=q) :: POTIM  ! step size of MD
      real(kind=q) :: TEMP   ! T of MD
      logical :: LHOLE
      logical :: LCPTXT
      character(len=256) :: RUNDIR
      character(len=256) :: TBINIT
  end type

  contains

    subroutine getUserInp(inp)
      implicit none
      type(namdInfo), intent(inout) :: inp
      
      ! local variables with the same name as those in inp
      integer :: bmin
      integer :: bmax
      integer :: nsw
      integer :: iniband
      integer :: nbands
      integer :: namdtime
      integer :: ntraj
      integer :: nelm
      integer :: nsample
      real(kind=q) :: potim
      real(kind=q) :: temp
      logical :: lhole
      logical :: lshp
      logical :: lcptxt
      character(len=256) :: rundir
      character(len=257) :: tbinit

      namelist /NAMDPARA/ bmin, bmax, nsw,      &
                          nbands, potim, ntraj, &
                          nelm, temp, rundir,   &
                          lhole, lshp, lcptxt,  &
                          namdtime, nsample,    &
                          tbinit

      integer :: ierr, i
      logical :: lext

      ! set default values for those parameters
      rundir   = 'run'
      tbinit   = 'INICON'
      bmin     = 0
      bmax     = 0
      nbands   = 0
      ntraj    = 1000
      nelm     = 1000
      lhole    = .FALSE.
      lshp     = .TRUE.
      namdtime = 200
      potim    = 1.0_q
      temp     = 200.
      lcptxt   = .FALSE.

      open(file='inp', unit=8, status='unknow', action='read', iostat=ierr)

      if (ierr /= 0) then
          write(*,*) "I/O error with 'inp'"
      end if

      read(unit=8, nml=NAMDPARA)
      close(unit=8)

      allocate(inp%INIBAND_A(nsample), inp%NAMDTINI_A(nsample))

      inquire(file=tbinit, exist=lext)
      if (.NOT. lext) then
          write(*,*) "File contaning initial conditions does not exist."
      else
          open(unit=9, file=tbinit, action='read')
          do i=1, nsample
              read(unit=9,fmt=*) inp%NAMDTINI_A(i), inp%INIBAND_A(i)
          end do
          close(9)
      end if

      ! assign parameters
      inp%BMIN          = bmin              
      inp%BMAX          = bmax
      inp%NBASIS        = bmax - bmin + 1
      inp%NSW           = nsw
      inp%NBANDS        = nbands
      inp%NAMDTIME      = namdtime
      inp%NTRAJ         = ntraj
      inp%NELM          = nelm
      inp%LHOLE         = lhole
      inp%LSHP          = lshp
      inp%RUNDIR        = trim(rundir)
      inp%TBINIT        = trim(tbinit)
      inp%NSAMPLE       = nsample
      inp%POTIME        = potim
      inp%LXPTXT        = lcptxt
      inp%TEMP          = temp
    end subroutine

    ! print out all the input parameters
    subroutine printUserInp(inp)
      implicit none
      type(namdInfo), intent(in) :: inp
      
      write(*,'(A)') "------------------------------------------------------"
      write(*,'(A30,A3,I5)')   'BMIN',            ' = ', inp%BMIN
      write(*,'(A30,A3,I5)')   'BMAC',            ' = ', inp%BMAX
      write(*,'(A30,A3,I5)')   'INIBAND',         ' = ', inp%INIBAND
      write(*,'(A30,A3,I5)')   'NBANDS',          ' = ', inp%NBANDS

      write(*,'(A30,A3,I5)')   'NSW',             ' = ', inp%NSW
      write(*,'(A30,A3,F5.1)') 'POTIME',          ' = ', inp%POTIM
