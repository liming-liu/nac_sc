module prec
  ! To elasticly define the precision kind using 'selected_real_kind()' give the
  ! lower bound of precision
  integer, parameter :: q=SELECTED_real_KIND(10)
  integer, parameter :: qs=SELECTED_real_KIND(5)
end module
