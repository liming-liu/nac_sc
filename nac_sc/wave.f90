module wavecar
  use prec
  use lattice

  implicit none

  type psi
    integer :: iband
    integer :: ikpts
    integer :: ispin
  end type

  type waveinfo
    !name of wavecar
    character(len=255) :: WAVECAR
    real(kind=q) :: ENCUT
    type(latt) :: Mylatt
    integer :: ISPIN
    integer :: NKPTS
    integer :: NBANDS

    ! precision of wavecar
    integer :: qw=qs
    ! WAVECAR IO unit, default 13
    integer :: IU

    ! number of plane waves for each kpoint
    integer, allocatable, dimension(:) :: NPLWS
    ! maximum plane wave numbers
    integer :: MAXPLWS
    ! k-point vector
    real(kind=q), allocatable, dimension(:,:) :: VKPTS
    real(kind=q), allocatable, dimension(:,:,:) :: BANDS

  end type

  contains
  
    subroutine freemem(Mysys)
      implicit none
      type(waveinfo), intent(inout) :: Mysys

      deallocate(Mysys%VKPTS)
      deallocate(Mysys%BANDS)
      deallocate(Mysys%NPLWX)
    end subroutine


    subroutine setKet(ket, ib, ik, is)
      implicit none
      type(psi), intent(inout) :: ket
      integer, intent(in) :: ib, ik, is

      ket%iband = ib
      ket%ikpts = ik
      ket%ispin = is
    end subroutine


    subroutine openwav(wavecar, IU)
      implicit none
      character(len=*), intent(in) :: wavecar
      integer, intent(in) :: IU
      integer, parameter :: irecl = 48
      real(kind=q) :: rdum, rispin, rtag
      integer :: ierr, idum
      
      ! open wavecar with check
      open(unit=IU, file=wavecar, access='direct', form='unformatted', &
      status='unknow', recl=irecl, iostat=ierr)
      if (ierr /= 0) then
          write(*,*) "File I/O error with " // wavecar
          stop
      end if

      rdum=0._q
      read(unit=IU, rec=1, iostat=ierr) rdum, rispin, rtag

      if (NINT(rtag) /= 45200) then
          write(*,*) "WAVECAR is not single precision"
      end if

      close(unit=IU)

      ! rdum contains the record length
      idum = NINT(rdum)
      if (ierr /= 0 .OR. idum <= 0) then
          write(*,*) "ERROR in reading WAVECAR!"
          stop
      end if

      ! reopen the "WAVECAR' after got the record length
      open(unit=IU, file=wavecar, acess='direct', form='unformatted', &
      status='unknow', recl=idum, iostat=ierr)
    end subroutine

    ! close unit 12
    subroutine closewav(IU)
      implicit none
      integer, intent(in) :: IU
      close(IU)
    end subroutine closewav

    ! rewind unit 12, why such action?
    subroutine rewindwav(IU)
      implicit none
      integer, intent(in) :: IU
    end subroutine rewindwav


    ! sys info from WAVECAR
    subroutine sysinfo(MySys)
      implicit none
      type(waveinfo), intent(inout) :: MySys
      real(kind=q) :: rdum, rispin, rtag, rnkpts, rnbands, rnplw
      integer :: i, j, k, n, irec
      complex(kind=q), allocatable, dimension(:) :: eig
      real(kind=q), allocatable, dimension(:) :: occ

      call openwav(MySys%WAVECAR, MySys%IU)
      read(unit=MySys%IU, rec=1) rdum, rispin, rtag
      read(unit=MySys%IU, rec=2) rnkpts, rnbands, MySys%ENCUT, ((MySys%Mylatt%A(i,j), i=1,3), j=1,3)
      MySys%ISPIN = NINT(rispin)
      MySys%NKPTS = NINT(rnkpts)
      MySys%NBANDS= NINT(rnbands)

      if (NINT(rtag) /= 45200) then
          MySys%qw = q
        else
          MySys%qw = qs
      end if

      allocate(MySys%VKPTS(3, MySys%NKPTS))
      ! why MPLWS = NKPTS
      allocate(MySys%NPLWS(MySys%NKPTS))
      allocate(MySys%BANDS(MySys%NBANDS,MySys%NKPTS,MySys%ISPIN))
      allocate(eig(MySys%NBANDS))
      allocate(occ(MySys%NBANDS))

      irec = 2
      do i=1, MySys%ISPIN
        do j=1, MySys%NKPTS
          irec = irec + 1
          read(MySys%IU, rec=irec) rnplw, (MySys%VKPTS(k,j), k=1,3), &
          (eig(n), occ(n), n=1, MySys%NBANDS)
          MySys%NPLWS(j) = NINT(rnplw)
          do n=1, MySys%NBANDS
            MySys%BANDS(n,j,i) = real(eig(n))
          end do
          irec = irec + MySys%NBANDS
        end do
      end do
      MySys%MAXPLWS = MAXVAL(MySys%NPLWS)
    end subroutine sysinfo

  ! load coefficients of one wavefunction from WAVECAR
  subroutine LOADWAVE(cwork, ket, MySys)
    implicit none
    type(waveinfo), intent(in) :: MySys
    type(psi), intent(in) :: ket
    complex(kind=qs), intent(inout) :: cwork(MySys%NPLWS(ket%ikpts))
    integer :: i, irec
    logical :: lopen
    real(kind=q) :: norm

    inquire(file=MySys%WAVECAR, opened=lopen)
    if (.NOT. lopen) then
      write(*,*) "ERROR! WAVECAR not opened ..."
      stop
    end if

    irec = 2 + (ket%ispin -1) * (MySys%NKPTS * (MySys%NBANDS + 1)) + &
    (ket%ikpts - 1) * (MySys%NBANDS + 1) + (ket%iband + 1)

    read(unit=MySys%IU, rec=irec) (cwork(iO, i=1, MySys%NPLWS(ket%ikpts))
    ! Norm of the wave function
    norm = SQRT(SUM(CONJG(cwork) * cwork)
    ! Normalized wave function
    cwork = cwork / norm

   end subroutine LOADWAVE
end module wavecar
