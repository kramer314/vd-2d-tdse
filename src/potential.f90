module potential

  use progvars

  use precision, only: ip, fp

  implicit none

  public :: potential_xyt

contains

  pure complex(fp) function potential_xyt(i_x, i_y, i_t) result(val)
    integer(ip), intent(in) :: i_x, i_y, i_t

    real(fp) :: offset, depth

    depth = 10.0_fp
    offset = 1.0_fp

    val = -depth / sqrt(x_range(i_x)**2 + y_range(i_x)**2 + offset**2)

  end function potential_xyt
end module potential
