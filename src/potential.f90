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
       val = val + potential_sin2_pulse_2d(i_x, i_y, i_t)
    end if

  end function potential_xyt

  pure real(fp) function potential_sin2_pulse_2d(i_x, i_y, i_t) result(val)
    integer(ip), intent(in) :: i_x, i_y, i_t

    real(fp) :: x, y, t

    real(fp) :: env, pol, int
    real(fp) :: t_ramp, t_plat, t_c
    real(fp) :: phi, omega

    real(fp) :: pulse

    x = x_range(i_x)
    y = y_range(i_y)
    t = t_range(i_t)

    int = pulse_intensity
    t_ramp = pulse_t_ramp
    t_plat = pulse_t_plateau

    phi = pulse_phi
    omega = pulse_omega

    ! Pulse envelope
    if (t .lt. t_ramp) then
       ! sin^2 envelope on up ramp
       env = sin(pi / 2.0_fp * t / t_ramp)**2
    else if (t .ge. t_ramp .and. t .le. t_ramp + t_plat) then
       ! constant envelope on plateau
       env = 1.0_fp
    else if (t .gt. t_ramp + t_plat .and. t .le. 2 * t_ramp + t_plat) then
       ! sin^2 envelope on down ramp
       env = sin(pi / 2.0_fp * (t - t_plat)/t_ramp)**2
    else
       env = 0.0_fp
    end if

    ! Right circular polarization
    pol = 1.0_fp

    t_c = t_ramp + t_plat / 2.0_fp

    pulse = env * sqrt(int)

    ! Dipole approximation: V = -q \vec{r} \dot \vec{E(t)}
    ! Here, q = -q_e
    val = pulse * ( x * cos(omega * (t - t_c) + pulse_phi) + &
         y * sin(omega * (t - t_c) + pulse_phi) )

  end function potential_sin2_pulse_2d

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
