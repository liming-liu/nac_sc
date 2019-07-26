Program  main
  use prec
  use lattice
  use wavecar
  use nac

  implicit none

  type(namdInfo) :: inp
  type(TDKS) :: ks
  type(overlap) :: olap. olap_sec

  real(kind=q) :: start, fin
  integer :: ns

  call TDCoupIJ(trim(inp
