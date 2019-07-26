module wavecar
use data_type

  contains

  subroutine wavecoef(wave, coeff)
    implicit none
    type(waveinfo) :: wave
    real*8, dimension(3) :: a1, a2, a3, b1, b2, b3
    real*8, dimension(3) :: vtmp, sumkg, wk
    integer, allocatable :: igall(:,:) ! what is this?
    complex*8, allocatable :: coeff(:,:)
    complex*8, allocatable :: occ(:), cener(:)
    integer :: nrecl, irec, iost, nspin, nprec, nwk, nband, npmax, nplane
    integer :: ncnt, ig1, ig2, ig3, ig1p, ig2p, ig3p
    real*8  :: pi, ecut, Vcell, gtot, etot
    real*8  :: phi12, phi13, phi23, sinphi123
    real*8  :: nb1maxA, nb1maxB, nb1maxC, &
               nb2maxA, nb2maxB, nb2maxC, &
               nb3maxA, nb3maxB, nb3maxC, &
               nb1max, nb1max2max, nb3max, npmax

    ! constant c = 2m/hbar**2 in units of 1/eV Ang^2
    data c/0.262465831d0/

    pi=4.*atan(1.)
    nrecl=24

    ! open WAVECAR
    open(unit=10, file='WAVECAR', access='direct', recl=nrecl,&
         iostat=iost, status='old', form='unformatted')
    ! check iostate
    if (iost .ne. 0) write "Open WAVECAR error, iostat =", iost
    ! read in record length, spin index, precision info
    read(unit=10, rec=1) nrecl, nspin, nprec
!    read(unit=10, rec=1) xnrecl, xnspin, xnprec
!    nrecl=nint(xnrecl)
!    nspin=nint(xnspin)
!    nprec=nint(xnprec)
    ! check if WAVECAR is single precision (8 byte)
    if (nprec .ne. 45200) then
        write(*,*) 'ERROR: WAVECAR is not single precision.'
        stop
    endif
    close(unit=10)
    ! reopen WAVECAR with new record length
    open(unit=10, file='WAVECAR', access='direct', recl=nrecl,&
         iostat=iost, status='old', form='unformatted')
    ! read in No. of kpoints, No. of bands, planewave basis set &
    ! energy cut, lattice vector
    read(unit=10, rec=2) nwk, nband, ecut, (a1(j), j=1, 3), &
         (a2(j), j=1, 3), (a3(j), j=1,3)
!    read(unit=10, rec=2) xnwk, xnband, ecut, (a1(j), j=1, 3), &
!         (a2(j), j=1, 3), (a3(j), j=1,3)
!    nwk=nint(xnwk)
!    nband=nint(xnband)
    ! allocate occ and cener
    allocate(occ(nband))
    allocate(cener(nband))
    ! compute direct cell volume
    call vcross(vtmp, a2, a3)
    Vcell=sum(a1(j)*vtmp(j), j=1, 3)
    ! compute reciprocal lattice basis
    call vcross(b1, a2, a3)
    call vcross(b2, a3, a1)
    call vcross(b3, a1, a2)
    do j=1, 3
        b1(j)=2.*pi*b1(j)/Vcell
        b2(j)=2.*pi*b2(j)/Vcell
        b3(j)=2.*pi*b3(j)/Vcell
    enddo
    b1mag=dsqrt(b1(1)**2+b1(2)**2+b1(3)**2)
    b2mag=dsqrt(b2(1)**2+b2(2)**2+b2(3)**2)
    b3mag=dsqrt(b3(1)**2+b3(2)**2+b3(3)**2)

    ! compute planewave basis size from ecut
    phi12=acos((b1(1)*b2(1)+b1(2)*b2(2)+b1(3)*b2(3))/(b1mag*b2mag))
    call vcross(vtmp, b1, b2)
    vmag=dsqrt(vtmp(1)**2+vtmp(2)**2+vtmp(3)**2)
    sinphi123=(b3(1)*vtmp(1)+b3(2)*vtmp(2)+b3(3)*vtmp(3))/(b3mag*vmag)
    nb1maxA=(dsqrt(ecut*c)/(b1mag*abs(sin(phi12))))+1
    nb2maxA=(dsqrt(ecut*c)/(b2mag*abs(sin(phi12))))+1
    nb3maxA=(dsqrt(ecut*c)/(b3mag*abs(sinphi123)))+1
    npmaxA=nint(4.*pi*nb1maxA*nb2maxA*nb3maxA/3.)

    phi13=acos((b1(1)*b3(1)+b1(2)*b3(2)+b1(3)*b3(3))/(b1mag*b3mag))
    call vcross(vtmp, b1, b3)
    vmag=dsqrt(vtmp(1)**2+vtmp(2)**2+vtmp(3)**2)
    sinphi123=(b2(1)*vtmp(1)+b2(2)*vtmp(2)+b2(3)*vtmp(3))/(b2mag*vmag)
    nb1maxB=(dsqrt(ecut*c)/(b1mag*abs(sin(phi13))))+1
    nb2maxB=(dsqrt(ecut*c)/(b2mag*abs(sinphi123)))+1
    nb3maxB=(dsqrt(ecut*c)/(b3mag*abs(sin(phi13))))+1
    npmaxB=nint(4.*pi*nb1maxB*nb2maxB*nb3maxB/3.)
    
    phi23=acos((b3(1)*b2(1)+b3(2)*b2(2)+b3(3)*b2(3))/(b3mag*b2mag))
    call vcross(vtmp, b2, b3)
    vmag=dsqrt(vtmp(1)**2+vtmp(2)**2+vtmp(3)**2)
    sinphi123=(b1(1)*vtmp(1)+b1(2)*vtmp(2)+b1(3)*vtmp(3))/(b1mag*vmag)
    nb1maxC=(dsqrt(ecut*c)/(b1mag*abs(sinphi123)))+1
    nb2maxC=(dsqrt(ecut*c)/(b2mag*abs(sin(phi23))))+1
    nb3maxC=(dsqrt(ecut*c)/(b3mag*abs(sin(phi23))))+1
    npmaxC=nint(4.*pi*nb1maxC*nb2maxC*nb3maxC/3.)

    nb1max=max0(nb1maxA, nb1maxB, nb1maxC)
    nb2max=max0(nb2maxA, nb2maxB, nb2maxC)
    nb3max=max0(nb3maxA, nb3maxB, nb3maxC)
    npmax=min0(nb1max, nb2max, nb3max)

    allocate (igall(3, npmax))
    allocate (coeff(npmax, nband))

    ! initial record location
    irec=2
    ! loop over spin
    do isp=1, nspin
        ! loop over k-points
        do iwk=1, nwk
            irec=irec+1         ! skip first 2 records
            ! there is no recl. And read in energy and occ 
            read(unit=10, rec=irec) nplane, (wk(i), i=1, 3), &
                (cener(iband), occ(iband), iband=1, nband)
#            read(unit=10, rec=irec) xnplane, (wk(i), i=1, 3), &
#                (cener(iband), occ(iband), iband=1, nband)
#            nplane=nint(xnplane)
            ! plane waves coef
            ncnt=0
            ! loop over k_z
            do ig3=0, 2*nb3max
                ig3p=ig3
                if (ig3 .gt. nb3max) ig3p=ig3-2*nb3max-1
                do ig2=0, 2*nb2max
                    ig2p=ig2
                    if (ig2 .gt. nb2max) ig2p=ig2-2*nb2max-1
                    do ig1=0, 2*nb1max
                        ig1p=ig1
                        if (ig1 .gt. nb1max) ig1p=ig1-2*nb1max-1
                        do j=1, 3
                            ! calc new reciprocal vector
                            sumkg(j)=(wk(1)+ig1p)*b1(j)+&
                                     (wk(2)+ig2p)*b2(j)+&
                                     (wk(3)+ig3p)*b3(j)
                        enddo
                        ! calculate the mag of the reciprocal vector
                        gtot=sqrt(sumkg(1)**2+sumkg(2)**2+sumkg(3)**2)
                        ! calc the en(k)
                        etot=gtot**2/c
                        ! check if the k is in range of encut, if dose, save k
                        ! in igall
                        if (etot .lt. ecut) then
                            ncnt=ncnt+1
                            igall(1, ncnt)=ig1p
                            igall(2, ncnt)=ig2p
                            igall(3, ncnt)=ig3p
                        endif
                    enddo
                enddo
            enddo
            ! check basis size
            if (ncnt .ne. nplane) then
                write(6,*) 'ERROR, calculated No. of basis size != input one'
                stop
            endif
            ! read in plane wave coeffs for each band
            do iband=1, nband
                irec=irec+1
                read(unit=10, rec=irec) (coeff(iplane, iband), &
                                        iplane=1, nplane)
            enddo
        enddo
    enddo
  end subroutine wavecoef

  subroutine vcross(a, b, c)
    implicit none
    real*8, dimension(3) :: a, b, c

    a(1)=b(2)*c(3)-b(3)*c(2)
    a(2)=b(3)*c(1)-b(1)*c(3)
    a(3)=b(1)*c(2)-b(2)*c(1)
    
    return
  end subroutine vcross

! calculate nac matrix
! C_ij = <psi_i(t)|d/dt|psi_j(t+dt)> 
!      = (<psi_i(t)|psi_j(t+dt)> - <psi_i(t)|psi_j(t)>) / dt
!      = (<psi_i(t)|psi_j(t+dt)> - <psi_i(t)|psi_j(t)> + 
!         <psi_i(t)|psi_j(t)> - <psi_i(t+dt)|psi_j(t)> ) / (2*dt)
!      = (<psi_i(t)|psi_j(t+dt)> - <psi_j(t)|psi_i(t+dt)>) / (2*dt)
! 1. get the psi_i(t) i=1, nband as for psi_j(t+dt), j=1, nband
! 2. calculate nac matrix


    subroutine setKet()
        implicit none
        type(psi), intent(inout) :: ket
        integer, intent(in) :: ib, ik, is

        ket%iband=ib
        ket%ikpts=ik
        ket%ispin=is
        
    end subroutine setKet

    subroutine LOADWAVE()
        implicit none
        
        type(waveinfo), intent(in) :: MySys
        type(psi), intent(in) :: ket
        complex*8, intent(inout) :: cwork(MySys%NPLWS(ket%ikpts))

        integer :: i, irec
        logical :: lopen
        real*8  :: norm

        inquire(file=Mysys%WAVECAR, opened=lopen)
        if (.NOT. lopen) then
            write(*,*) "ERROR! WAVECAR not opened ..."
            stop
        end if

        irec=2+(ket%ispin-1)*(MySys%NKPTS*(MySys%NBNDS+1))+ &
               (ket%ikpts-1)*(MySys%NBANDS+1)+(ket%iband+1)

        ! write(*,*) rec, ket%ispin, ket%ikpts, ket%iband
        read(unit=MySys%IU, rec=irec) (cwork(i), i=1, MySys%NPLWS(ket%ikpts))
        ! Norm of the wave function
        norm=sqrt(sum(conjg(cwork)*cwork))
        ! Normalize wave function
        cwork=cwork/norm
        ! Unify the pase of the single-particle wave function
        ! cwork=cwork*conjg(cwork(1))/abs(cwork(1))
        
    end subroutine LOADWAVE

    subroutine calc_nac(waveA, waveB, Cij)
    implicit none
    type(waveinfo), intent(in) :: waveA, waveB

    integer :: i, j
    type(psi), allocatable, dimension(:) :: ket
    real*8, dimension(:,:), intent(inout) :: Cij
    complex*8, allocatable, dimension(:,:) :: crA, crB
    real*8 :: pij, pji, Oii
    integer, allocatable, dimesion(:) :: pha_correction

    allocate(crA(waveA%NPLWS(1), waveA%NBANDS))
    allocate(crB(waveB%NPLWS(1), waveB%NBANDS))
    
    allocate(ket(waveA%NBANDS))

    ! read in all the wavefunctions
    ! the gamma point wavecar has only one point

    ! $omp parallel do
    do i=1, waveA%NBANDS
        ! i-th band, first kpoint, first spin
        call setKet(ket(i), i, 1, 1)
        ! the coefficients are normalized in the LOADWAVE subroutine
        ! here, we don't have to worry about normalization problem
        call LOADWAVE(crA(:, i), ket(i), waveA)
        call LOADWAVE(crB(:, i), ket(i), waveB)
    end do
    ! $omp end parallel do

    ! Initialization
    Cij=0
    !  <psi_i(t)|d/dt|psi_j(t)>
    ! ~<psi_i(t)|psi_j(t+dt)>-<psi_i(t+dt)|psi_j(t)>/(2dt)
    write(*,*) "#", trim(waveA%WAVECAR)

    ! $omp parallel do
    do i=1, waveA%NBNDS
        ! write(*, *) "C_ZERO: ", REAL(crA(1, i), AIMAG(crA(1,i))
        do j=i+1, waveA%NBANDS
            ! <psi_i(t)|psi_j(t+dt)>
            pij=sum(conjg(crA(:, i)) * crB(:, j)) * pha_correction(j)
            pji=sum(conjg(crA(:, j)) * crB(:, i)) * pha_correction(i)

            ! not devided by 2*dt
            Cij(i, j)=pij-pji
            if (i/=j) then
                Cij(j, i)=-Cij(i, j)
            end if
        end do
    end do
    ! $omp end parallel do

    ! Cij=-0.658218*Cij/2.
    ! do i=1, waveA%NBNDS
    !     write(*, *) (Cij(i, j), j=1, waveA%NBANDS)
    ! end do

    deallocate(crA, crB, ket)
    deallocate(pha_correction)
    end subroutine calc_nac

end module wavecar 
