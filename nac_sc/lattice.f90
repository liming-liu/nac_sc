module lattice
  use prec
  type latt
    real(q) :: SCALE
    ! A: real space basis
    ! B: momentum space basis
    real(q) :: A(3,3), B(3,3)
    ! NORM of the basises
    real(q) :: ANORM(3), BNORM(3)
    ! volume of the cell
    real(q) :: OMEGA
  end type
  contains

    ! calculate cross production H = U1 X U2
    subroutine EXPRO(H, U1, U2)
      use prec
      implicit none
      real(kind=q), intent(in) :: U1, U2
      real(kind=q), intent(inout) :: H
      dimension H(3), U1(3), U2(3)

      H(1)=U1(2)*U2(3)-U1(3)*U2(2)
      H(2)=U1(3)*U2(1)-U1(1)*U2(3)
      H(3)=U1(1)*U2(2)-U1(2)*U2(1)
      
      return
    end subroutine

  ! calculating reciprocal lattice from the direct lattice
    subroutine lattic(Mylatt)
      use prec
      implicit none
      type(latt) :: Mylatt
      real(q) :: Omega
      integer :: i,j
      intrinsic SUM

      ! calculate reciprocal basis form real one
      CALL EXPRO(Mylatt%B(:1), Mylatt%A(:,2), Mylatt%A(:,3))
      CALL EXPRO(Mylatt%B(:2), Mylatt%A(:,3), Mylatt%A(:,1))
      CALL EXPRO(Mylatt%B(:3), Mylatt%A(:,1), Mylatt%A(:,2))

      ! calculate volume
      Omega=SUM(Mylatt%A(:,1) * Mylatt%B(:,1))

      ! normalize B
      Mylatt%B = Mylatt%B / Omega

      ! magnitude of basis vector
      do i=1,3
        Mylatt%ANORM(i)=SQRT(SUM(Mylatt%A(:,i)*Mylatt%A(:,i)))
        Mylatt%BNORM(i)=SQRT(SUM(Mylatt%B(:,i)*Mylatt%B(:,i)))
      end do

      Mylatt%OMEGA=Omega

      return
    end subroutine

end module lattice
