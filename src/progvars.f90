! Copyright (c) 2016 Alex Kramer <kramer.alex.kramer@gmail.com>
! See the LICENSE.txt file at the top-level of this distribution.

module progvars
  ! Program variables module

  ! Imports -- library dependencies
  use globvars, only: pi_dp, e_dp, j_dp
  use config, only: config_get_param
  use numerics, only: numerics_linspace

  ! Imports -- program modules
  use precision
  use vd_2d_manager, only: vd_2d_manager_obj

  implicit none

  ! Numerical constants
  real(fp), parameter :: pi = real(pi_dp, kind=fp)
  real(fp), parameter :: sqrt_pi = sqrt(pi)
  real(fp), parameter :: e = real(e_dp, kind=fp)
  complex(fp), parameter :: j = cmplx(j_dp, kind=fp)

  ! Units
  real(fp) :: hbar

  ! Particle parameters
  real(fp) :: m

  ! Gaussian parameters
  real(fp) :: p0_x, sig_x, x0
  real(fp) :: p0_y, sig_y, y0

  ! Grid parameters
  real(fp) :: x_min, x_max, dx
  real(fp) :: y_min, y_max, dy
  integer(ip) :: nx, ny

  ! Virtual detector object along x spatial grid
  type(vd_2d_manager_obj) :: vd_xy

  ! Virtual detector parameters
  ! left/right number of grid points outside region of interested
  integer(ip) :: nxl_external, nxr_external
  integer(ip) :: nyl_external, nyr_external
  ! left/right number of virtual detector grid points in external grid
  integer(ip) :: nxl_vd, nxr_vd
  integer(ip) :: nyl_vd, nyr_vd

  ! disjoint VD switch
  logical :: vd_disjoint

  ! semi-classical VD switch
  logical :: vd_semi_classical

  ! Useful VD grid indices
  integer(ip) :: vd_xl_min, vd_xl_max, vd_xr_min, vd_xr_max
  integer(ip) :: vd_yl_min, vd_yl_max, vd_yr_min, vd_yr_max

  ! VD momentum bin parameters
  real(fp) :: vd_px_min, vd_px_max
  real(fp) :: vd_py_min, vd_py_max
  integer(ip) :: vd_npx, vd_npy

  ! Number of standard deviations to include in quantum VD binning
  integer(ip) :: vd_np_stdev

  ! Flux threshold for counting contributions
  real(fp) :: vd_j_eps

  ! Time grid parameters
  real(fp) :: t_min, t_max, dt
  integer(ip) :: nt

  ! Output parameters
  integer(ip) :: print_mod_x, print_mod_y, print_mod_t

  character(:), allocatable :: output_dir
  character(:), allocatable :: log_fname
  character(:), allocatable :: vd_p_fname
  character(:), allocatable :: vd_resid_fname
  character(:), allocatable :: vd_resid_cum_fname
  character(:), allocatable :: vd_resid_analysis_fname

  ! Residual analysis parameters
  real(fp) :: resid_p_eps

  ! Residual analysis variables
  real(fp) :: resid_mean, resid_var
  real(fp) :: resid_mean_sq_err, resid_mean_abs_err
  real(fp) :: resid_fivenum_arr(5)

  ! Arrays
  real(fp), allocatable :: x_range(:), y_range(:), t_range(:)
  complex(fp), allocatable :: psi_arr(:,:)

  real(fp), allocatable :: theor_np_arr(:,:), resid_np_arr(:,:)
  real(fp), allocatable :: resid_np_cum_arr(:,:)
  logical, allocatable :: resid_np_mask(:,:)

contains

  subroutine progvars_init()
    call progvars_read_params()
    call progvars_allocate_arrays()
    call progvars_set_arrays()

    call vd_xy%init(nx, ny, nxl_external, nxr_external, nxl_vd, nxr_vd, &
         nyl_external, nyr_external, nyl_vd, nyr_vd, dx, dy, vd_npx, vd_npy, &
         vd_px_min, vd_px_max, vd_py_min, vd_py_max, dt, vd_semi_classical, &
         vd_disjoint, hbar, m, vd_np_stdev=vd_np_stdev)

    vd_xl_min = vd_xy%xl_min
    vd_xl_max = vd_xy%xl_max
    vd_xr_min = vd_xy%xr_min
    vd_xr_max = vd_xy%xr_max

    vd_yl_min = vd_xy%yl_min
    vd_yl_max = vd_xy%yl_max
    vd_yr_min = vd_xy%yr_min
    vd_yr_max = vd_xy%yr_max


  end subroutine progvars_init

  subroutine progvars_cleanup()
    call vd_xy%cleanup()
    call progvars_deallocate_arrays()
    call progvars_deallocate_params()
  end subroutine progvars_cleanup

  subroutine progvars_allocate_arrays()
    allocate(psi_arr(nx, ny))
    allocate(t_range(nt))
    allocate(x_range(nx))
    allocate(y_range(ny))

    allocate(theor_np_arr(vd_npx, vd_npy))
    allocate(resid_np_arr(vd_npx, vd_npy))
    allocate(resid_np_cum_arr(vd_npx, vd_npy))
    allocate(resid_np_mask(vd_npx, vd_npy))
  end subroutine progvars_allocate_arrays

  subroutine progvars_deallocate_arrays()
    deallocate(x_range)
    deallocate(y_range)
    deallocate(t_range)
    deallocate(psi_arr)

    deallocate(theor_np_arr)
    deallocate(resid_np_arr)
    deallocate(resid_np_cum_arr)
    deallocate(resid_np_mask)
  end subroutine progvars_deallocate_arrays

  subroutine progvars_deallocate_params
    deallocate(output_dir)
    deallocate(log_fname)
    deallocate(vd_p_fname)
    deallocate(vd_resid_fname)
    deallocate(vd_resid_cum_fname)
    deallocate(vd_resid_analysis_fname)
  end subroutine progvars_deallocate_params

  subroutine progvars_set_arrays()
    ! Initialize numerical grids
    call numerics_linspace(t_min, t_max, t_range, dt)
    call numerics_linspace(x_min, x_max, x_range, dx)
    call numerics_linspace(y_min, y_max, y_range, dy)
  end subroutine progvars_set_arrays

  subroutine progvars_read_params()
    ! Read in parameters from input file
    logical :: success

    call config_get_param("hbar", hbar, success)

    call config_get_param("m", m, success)

    call config_get_param("p0_x", p0_x, success)
    call config_get_param("x0", x0, success)
    call config_get_param("sig_x", sig_x, success)

    call config_get_param("p0_y", p0_y, success)
    call config_get_param("y0", y0, success)
    call config_get_param("sig_y", sig_y, success)

    call config_get_param("x_min", x_min, success)
    call config_get_param("x_max", x_max, success)
    call config_get_param("nx", nx, success)

    call config_get_param("y_min", y_min, success)
    call config_get_param("y_max", y_max, success)
    call config_get_param("ny", ny, success)

    call config_get_param("t_min", t_min, success)
    call config_get_param("t_max", t_max, success)
    call config_get_param("nt", nt, success)

    call config_get_param("nxl_external", nxl_external, success)
    call config_get_param("nxr_external", nxr_external, success)
    call config_get_param("nxl_vd", nxl_vd, success)
    call config_get_param("nxr_vd", nxr_vd, success)

    call config_get_param("nyl_external", nyl_external, success)
    call config_get_param("nyr_external", nyr_external, success)
    call config_get_param("nyl_vd", nyl_vd, success)
    call config_get_param("nyr_vd", nyr_vd, success)

    call config_get_param("vd_px_min", vd_px_min, success)
    call config_get_param("vd_px_max", vd_px_max, success)
    call config_get_param("vd_npx", vd_npx, success)

    call config_get_param("vd_py_min", vd_py_min, success)
    call config_get_param("vd_py_max", vd_py_max, success)
    call config_get_param("vd_npy", vd_npy, success)

    call config_get_param("vd_j_eps", vd_j_eps, success)

    call config_get_param("vd_np_stdev", vd_np_stdev, success)

    call config_get_param("vd_semi_classical", vd_semi_classical, success)

    call config_get_param("vd_disjoint", vd_disjoint, success)

    call config_get_param("resid_p_eps", resid_p_eps, success)

    call config_get_param("output_dir", output_dir, success)
    call config_get_param("log_fname", log_fname, success)
    call config_get_param("vd_p_fname", vd_p_fname, success)
    call config_get_param("vd_resid_fname", vd_resid_fname, success)
    call config_get_param("vd_resid_cum_fname", vd_resid_cum_fname, success)
    call config_get_param("vd_resid_analysis_fname", vd_resid_analysis_fname, &
         success)

    call config_get_param("print_mod_x", print_mod_x, success)
    call config_get_param("print_mod_y", print_mod_y, success)
    call config_get_param("print_mod_t", print_mod_t, success)

  end subroutine progvars_read_params

end module progvars
