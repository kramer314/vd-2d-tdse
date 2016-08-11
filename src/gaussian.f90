! Copyright (c) 2016 Alex Kramer <kramer.alex.kramer@gmail.com>
! See the LICENSE.txt file at the top-level of this distribution.

module gaussian
  ! 1D Gaussian module

  use progvars

  implicit none

  public :: gaussian_init
  public :: gaussian_cleanup
  public :: gaussian_pxpy
  public :: gaussian_xyt

  private

  ! Useful constants
  real(fp) :: p0_x_m, p0_x_2m
  real(fp) :: p0_y_m, p0_y_2m
  real(fp) :: sig_x2, sig_y2
  real(fp) :: sig_px, sig_px2
  real(fp) :: sig_py, sig_py2
  real(fp) :: norm_px, norm_py
  complex(fp) :: j_hb, jhb_m

contains
  subroutine gaussian_init()
    ! Precalculate useful constants
    p0_x_m = p0_x / m
    p0_x_2m = p0_x_m / 2.0_fp

    p0_y_m = p0_y / m
    p0_y_2m = p0_y_m / 2.0_fp

    sig_x2 = sig_x ** 2
    sig_y2 = sig_y ** 2

    j_hb = j / hbar
    jhb_m = j * hbar / m

    sig_px = hbar / sig_x
    sig_px2 = sig_px**2

    sig_py = hbar / sig_y
    sig_py2 = sig_py**2

    norm_px = sqrt(1.0_fp / (sig_px * sqrt_pi))
    norm_py = sqrt(1.0_fp / (sig_py * sqrt_pi))

  end subroutine gaussian_init

  subroutine gaussian_cleanup()
  end subroutine gaussian_cleanup

  real(fp) function gaussian_pxpy(px, py) result(val)
    ! Taken from Sakurai, p. 58
    real(fp), intent(in) :: px, py
    real(fp) :: exp_px, exp_py

    exp_px = - 0.5_fp * ( (px - p0_x) / sig_px )**2
    exp_py = - 0.5_fp * ( (py - p0_y) / sig_py )**2

    val = norm_px * norm_py * exp(exp_px + exp_py)
  end function gaussian_pxpy

  complex(fp) function gaussian_xyt(x, y, t) result(val)
    real(fp), intent(in) :: x, y
    real(fp), intent(in) :: t

    complex(fp) :: norm_x, norm_y
    complex(fp) :: exp_x1, exp_x2, exp_y1, exp_y2

    norm_x = 1.0_fp / sqrt( sqrt_pi * (sig_x + jhb_m / sig_x * t) )
    norm_y = 1.0_fp / sqrt( sqrt_pi * (sig_y + jhb_m / sig_y * t) )

    exp_x1 = - ( (x - x0) - p0_x_m * t )**2 / &
         ( 2.0_fp * sig_x2 * (1 + jhb_m / sig_x2 * t) )
    exp_y1 = - ( (y - y0) - p0_y_m * t )**2 / &
         ( 2.0_fp * sig_y2 * (1 + jhb_m / sig_y2 * t) )

    exp_x2 = j_hb * p0_x * ( (x - x0) - p0_x_2m * t )
    exp_y2 = j_hb * p0_y * ( (y - y0) - p0_y_2m * t )

    val = norm_x * norm_y * exp(exp_x1 + exp_x2 + exp_y1 + exp_y2)

  end function gaussian_xyt

end module gaussian
