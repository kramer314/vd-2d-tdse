module propagate_cn2d
  use progvars

  use precision, only: fp, ip, fp_eps

  use potential, only: potential_xyt
  use tridiag, only: tridiag_sym_cnst, tridiag_sym_cnst_odiag

  implicit none

  private

  complex(dp), allocatable :: diag_arr_y(:)

  complex(dp), allocatable :: phi_arr_x(:), phi_arr_y(:)
  complex(dp), allocatable :: psi_arr_x(:), psi_arr_y(:)

  complex(dp), allocatable :: mat_coeff_arr_x(:), mat_coeff_arr_y(:)
  complex(dp), allocatable :: vec_coeff_arr_x(:), vec_coeff_arr_y(:)

  complex(fp) :: sym_cnst_dt, sym_cnst_dt_2
  complex(fp) :: sym_cnst_x_dt, sym_cnst_x_dt_2
  complex(fp) :: sym_cnst_y_dt, sym_cnst_y_dt_2
  complex(fp) :: diag_cnst_x_dt, diag_cnst_x_dt_2
  complex(fp) :: diag_cnst_y_dt, diag_cnst_y_dt_2
  complex(fp) :: diag_cnst_vy_dt

  public :: propagate_cn2d_init
  public :: propagate_cn2d_cleanup
  public :: propagate_cn2d_splitop

contains
  subroutine propagate_cn2d_init()
    integer(ip) :: i_x, i_y

    sym_cnst_dt = - (j * hbar**2 * dt) / (8.0_fp * m)
    sym_cnst_dt_2 = sym_cnst_dt / 2.0_fp

    sym_cnst_x_dt = sym_cnst_dt / dx**2
    sym_cnst_x_dt_2 = sym_cnst_dt_2 / dx **2

    sym_cnst_y_dt = sym_cnst_dt / dy**2
    sym_cnst_y_dt_2 = sym_cnst_dt_2 / dy **2

    diag_cnst_x_dt = (0.5_fp - 2.0_dp * sym_cnst_x_dt)
    diag_cnst_x_dt_2 = (0.5_fp - 2.0_dp * sym_cnst_x_dt_2)

    diag_cnst_y_dt = (0.5_fp - 2.0_dp * sym_cnst_y_dt)
    diag_cnst_y_dt_2 = (0.5_fp - 2.0_dp * sym_cnst_y_dt_2)

    diag_cnst_vy_dt = j * dt / 4.0_fp

    allocate(phi_arr_x(nx), phi_arr_y(ny))
    allocate(psi_arr_x(nx), psi_arr_y(ny))

    allocate(diag_arr_y(ny))

    allocate(mat_coeff_arr_x(nx - 1))
    allocate(mat_coeff_arr_y(ny - 1))
    allocate(vec_coeff_arr_x(nx))
    allocate(vec_coeff_arr_y(ny))

  end subroutine propagate_cn2d_init

  subroutine propagate_cn2d_cleanup()
    deallocate(phi_arr_x, phi_arr_y)
    deallocate(psi_arr_x, psi_arr_y)
    deallocate(diag_arr_y)
    deallocate(mat_coeff_arr_x, mat_coeff_arr_y)
    deallocate(vec_coeff_arr_x, vec_coeff_arr_y)
  end subroutine propagate_cn2d_cleanup

  subroutine propagate_cn2d_splitop(psi_arr, i_t)
    complex(dp), intent(inout) :: psi_arr(:,:)
    integer(ip), intent(in) :: i_t

    integer(ip) :: i_x, i_y

    ! Propagate dt/2 using T_x for each y
    !$omp parallel do &
    !$omp private(i_y, psi_arr_x, phi_arr_x, mat_coeff_arr_x, vec_coeff_arr_x)
    do i_y = 1, ny
       psi_arr_x(:) = psi_arr(:, i_y)
       call tridiag_sym_cnst(diag_cnst_x_dt_2, sym_cnst_x_dt_2, &
            psi_arr_x, mat_coeff_arr_x, vec_coeff_arr_x, &
            phi_arr_x)
       ! call tridiag_cnst(diag_cnst_x_dt_2, sym_cnst_x_dt_2, &
       !      sym_cnst_x_dt_2, psi_arr_x, mat_coeff_arr_x, vec_coeff_arr_x, &
       !      phi_arr_x)
       psi_arr(:, i_y) = phi_arr_x(:) - psi_arr(:, i_y)
    end do
    !$omp end parallel do

    ! Propagate dt using T_y + V for each x
    !$omp parallel do &
    !$omp private(i_x, psi_arr_y, phi_arr_y, mat_coeff_arr_y, vec_coeff_arr_y, diag_arr_y)
    do i_x = 1, nx
       psi_arr_y(:) = psi_arr(i_x, :)

       do i_y = 1, ny
          diag_arr_y(i_y) = diag_cnst_y_dt + diag_cnst_vy_dt * potential_xyt(i_x, i_y, i_t)
       end do

       call tridiag_sym_cnst_odiag(diag_arr_y, sym_cnst_y_dt, psi_arr_y, &
            mat_coeff_arr_y, vec_coeff_arr_y, phi_arr_y)
       psi_arr(i_x, :) = phi_arr_y(:) - psi_arr(i_x, :)
    end do
    !$omp end parallel do

    ! Propagate dt/2 using T_x for each yx
    !$omp parallel do &
    !$omp private(i_y, psi_arr_x, phi_arr_x, mat_coeff_arr_x, vec_coeff_arr_x)
    do i_y = 1, ny
       psi_arr_x(:) = psi_arr(:, i_y)
       call tridiag_sym_cnst(diag_cnst_x_dt_2, sym_cnst_x_dt_2, &
            psi_arr_x, mat_coeff_arr_x, vec_coeff_arr_x, &
            phi_arr_x)
       ! call tridiag_cnst(diag_cnst_x_dt_2, sym_cnst_x_dt_2, &
       !      sym_cnst_x_dt_2, psi_arr_x, mat_coeff_arr_x, vec_coeff_arr_x, &
       !      phi_arr_x)
       psi_arr(:, i_y) = phi_arr_x(:) - psi_arr(:, i_y)
    end do
    !$omp end parallel do

  end subroutine propagate_cn2d_splitop

end module propagate_cn2d
