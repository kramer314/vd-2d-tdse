module propagate_cn2d_itime
  use progvars

  use precision, only: fp, ip, fp_eps

  use potential, only: potential_xyt
  use tridiag, only: tridiag_sym_cnst, tridiag_sym_cnst_odiag

  use wfmath, only: wfmath_proj

  implicit none

  private

  complex(dp), allocatable :: diag_arr_y(:)

  complex(dp), allocatable :: phi_arr_x(:), phi_arr_y(:)
  complex(dp), allocatable :: psi_arr_x(:), psi_arr_y(:)

  complex(dp), allocatable :: mat_coeff_arr_x(:), mat_coeff_arr_y(:)
  complex(dp), allocatable :: vec_coeff_arr_x(:), vec_coeff_arr_y(:)

  complex(fp) :: sym_cnst_dtau, sym_cnst_dtau_2
  complex(fp) :: sym_cnst_x_dtau, sym_cnst_x_dtau_2
  complex(fp) :: sym_cnst_y_dtau, sym_cnst_y_dtau_2
  complex(fp) :: diag_cnst_x_dtau, diag_cnst_x_dtau_2
  complex(fp) :: diag_cnst_y_dtau, diag_cnst_y_dtau_2
  complex(fp) :: diag_cnst_vy_dtau

  complex(fp) :: dtau

  public :: propagate_cn2d_itime_init
  public :: propagate_cn2d_itime_cleanup
  public :: propagate_cn2d_itime_splitop
  public :: propagate_cn2d_itime_converged

contains
  subroutine propagate_cn2d_itime_init()
    integer(ip) :: i_x, i_y

    dtau = -j * dt

    sym_cnst_dtau = - (j * hbar**2 * dtau) / (8.0_fp * m)
    sym_cnst_dtau_2 = sym_cnst_dtau / 2.0_fp

    sym_cnst_x_dtau = sym_cnst_dtau / dx**2
    sym_cnst_x_dtau_2 = sym_cnst_dtau_2 / dx **2

    sym_cnst_y_dtau = sym_cnst_dtau / dy**2
    sym_cnst_y_dtau_2 = sym_cnst_dtau_2 / dy **2

    diag_cnst_x_dtau = (0.5_fp - 2.0_dp * sym_cnst_x_dtau)
    diag_cnst_x_dtau_2 = (0.5_fp - 2.0_dp * sym_cnst_x_dtau_2)

    diag_cnst_y_dtau = (0.5_fp - 2.0_dp * sym_cnst_y_dtau)
    diag_cnst_y_dtau_2 = (0.5_fp - 2.0_dp * sym_cnst_y_dtau_2)

    diag_cnst_vy_dtau = j * dtau / 4.0_fp

    allocate(phi_arr_x(nx), phi_arr_y(ny))
    allocate(psi_arr_x(nx), psi_arr_y(ny))

    allocate(diag_arr_y(ny))

    allocate(mat_coeff_arr_x(nx - 1))
    allocate(mat_coeff_arr_y(ny - 1))
    allocate(vec_coeff_arr_x(nx))
    allocate(vec_coeff_arr_y(ny))

  end subroutine propagate_cn2d_itime_init

  subroutine propagate_cn2d_itime_cleanup()
    deallocate(phi_arr_x, phi_arr_y)
    deallocate(psi_arr_x, psi_arr_y)
    deallocate(diag_arr_y)
    deallocate(mat_coeff_arr_x, mat_coeff_arr_y)
    deallocate(vec_coeff_arr_x, vec_coeff_arr_y)
  end subroutine propagate_cn2d_itime_cleanup

  subroutine propagate_cn2d_itime_splitop(psi_arr, i_t)
    complex(dp), intent(inout) :: psi_arr(:,:)
    integer(ip), intent(in) :: i_t

    integer(ip) :: i_x, i_y

    ! Propagate dtau/2 using T_x for each y
    !$omp parallel do &
    !$omp private(i_y, psi_arr_x, phi_arr_x, mat_coeff_arr_x, vec_coeff_arr_x)
    do i_y = 1, ny
       psi_arr_x(:) = psi_arr(:, i_y)
       call tridiag_sym_cnst(diag_cnst_x_dtau_2, sym_cnst_x_dtau_2, &
            psi_arr_x, mat_coeff_arr_x, vec_coeff_arr_x, &
            phi_arr_x)
       psi_arr(:, i_y) = phi_arr_x(:) - psi_arr(:, i_y)
    end do
    !$omp end parallel do
    
    ! Propagate dtau using T_y + V for each x
    !$omp parallel do &
    !$omp private(i_x, psi_arr_y, phi_arr_y, mat_coeff_arr_y, vec_coeff_arr_y, diag_arr_y)
    do i_x = 1, nx
       psi_arr_y(:) = psi_arr(i_x, :)

       do i_y = 1, ny
          diag_arr_y(i_y) = diag_cnst_y_dtau + diag_cnst_vy_dtau * potential_xyt(i_x, i_y, i_t)
       end do

       call tridiag_sym_cnst_odiag(diag_arr_y, sym_cnst_y_dtau, psi_arr_y, &
            mat_coeff_arr_y, vec_coeff_arr_y, phi_arr_y)
       psi_arr(i_x, :) = phi_arr_y(:) - psi_arr(i_x, :)
    end do
    !$omp end parallel do
    
    ! Propagate dtau/2 using T_x for each yx
    !$omp parallel do &
    !$omp private(i_y, psi_arr_x, phi_arr_x, mat_coeff_arr_x, vec_coeff_arr_x)
    do i_y = 1, ny
       psi_arr_x(:) = psi_arr(:, i_y)
       call tridiag_sym_cnst(diag_cnst_x_dtau_2, sym_cnst_x_dtau_2, &
            psi_arr_x, mat_coeff_arr_x, vec_coeff_arr_x, &
            phi_arr_x)
       psi_arr(:, i_y) = phi_arr_x(:) - psi_arr(:, i_y)
    end do
    !$omp end parallel do

  end subroutine propagate_cn2d_itime_splitop

  logical function propagate_cn2d_itime_converged(psi_arr, old_psi_arr, eps) result(val)

    complex(fp), intent(in) :: psi_arr(:,:)
    complex(fp), intent(in) :: old_psi_arr(:,:)
    real(fp), intent(in) :: eps

    real(fp) :: proj

    proj = wfmath_proj(psi_arr, old_psi_arr, dgrid)
    
    val = .false.
    if (abs(1.0_fp - proj) .lt. eps) then
       val = .true.
    end if

  end function propagate_cn2d_itime_converged

end module propagate_cn2d_itime
