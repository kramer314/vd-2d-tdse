module potential
  use precision, only: ip, fp

  implicit none

  public :: potential_xyt

contains

  pure real(fp) function potential_xyt(i_x, i_y, i_t) result(val)
    integer(ip), intent(in) :: i_x, i_y, i_t

    val = 0.0_fp
  end function potential_xyt
end module potential
