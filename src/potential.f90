module potential

  use progvars

  use precision, only: ip, fp

  implicit none

  public :: potential_xyt

contains

  pure complex(fp) function potential_xyt(i_x, i_y, i_t) result(val)
    integer(ip), intent(in) :: i_x, i_y, i_t

    complex(fp) :: offset

    offset = 1.0_fp

    val = -1.0_fp / sqrt(x_range(i_x)**2 + y_range(i_x)**2 + offset**2)
    val = x_range(i_x)**2 * y_range(i_y)**2

  end function potential_xyt
end module potential
