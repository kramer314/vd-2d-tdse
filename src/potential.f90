module potential

  use progvars

  use precision, only: ip, fp

  implicit none

  public :: potential_xyt

contains

  pure complex(fp) function potential_xyt(i_x, i_y, i_t) result(val)
    integer(ip), intent(in) :: i_x, i_y, i_t

    val = potential_sc_coulomb_2d(i_x, i_y)

    if (i_t > 0_ip) then
       val = val + potential_gobbler_quartic_2d(i_x, i_y)
    end if

  end function potential_xyt

  pure real(fp) function potential_sc_coulomb_2d(i_x, i_y) result(val)
    integer(ip), intent(in) :: i_x, i_y

    real(fp) :: offset, depth
    real(fp) :: x, y

    x = (x_range(i_x) - sc_coulomb_x0)
    y = (y_range(i_y) - sc_coulomb_y0)

    val = - sc_coulomb_depth / sqrt(x**2 + y**2 + sc_coulomb_offset**2)

  end function potential_sc_coulomb_2d

  ! Negative imaginary gobbler
  pure complex(fp) function potential_gobbler_quartic_2d(i_x, i_y) result(val)

    integer(ip), intent(in) :: i_x, i_y

    integer(fp) :: delta_x, delta_y

    val = 0.0_fp

    if (i_x > vd_xr_max) then
       delta_x = x_range(i_x) - x_range(vd_xr_max)
       val = val * gobbler_strength_x * (delta_x / gobbler_width_y)**4
    else if (i_x < vd_xl_min) then
       delta_x = x_range(vd_xl_min) - x_range(i_x)
       val = val * gobbler_strength_x * (delta_x / gobbler_width_y)**4
    end if

    if (i_y > vd_yr_max) then
       delta_y = y_range(i_y) - y_range(vd_yr_max)
       val = val * gobbler_strength_x * (delta_y / gobbler_width_y)**4
    else if (i_x < vd_yl_min) then
       delta_y = y_range(vd_yl_min) - y_range(i_y)
       val = val * gobbler_strength_x * (delta_y / gobbler_width_y)**4
    end if

    val = -j * val

  end function potential_gobbler_quartic_2d

end module potential
