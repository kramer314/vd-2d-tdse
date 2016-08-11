module propagate
  use progvars

  use gaussian, only: gaussian_init, gaussian_cleanup, gaussian_xyt
  
  use propagate_cn2d, only: propagate_cn2d_init, propagate_cn2d_cleanup, &
       propagate_cn2d_splitop

  use numerics, only: numerics_trapz
  
  use precision, only: fp, ip

  implicit none

  public :: propagate_init
  public :: propagate_cleanup
  public :: propagate_psi
  
contains
  
  subroutine propagate_init(psi_arr)
    ! Initialize wfunc
    complex(fp), intent(inout) :: psi_arr(:,:)

    integer(ip) :: i_x, i_y
    real(fp) :: x, y

    call gaussian_init()
    call propagate_cn2d_init()

    !$omp parallel do private(i_x, i_y, x, y) firstprivate(x_range, y_range)
    do i_y = 1, ny
       y = y_range(i_y)
       do i_x = 1, nx
          x = x_range(i_x)
          psi_arr(i_x, i_y) = gaussian_xyt(x, y, 0.0_fp)
       end do
    end do
    !$omp end parallel do
  end subroutine propagate_init
  
  subroutine propagate_cleanup()
    call gaussian_cleanup()
    call propagate_cn2d_cleanup()
  end subroutine propagate_cleanup

  subroutine propagate_psi(psi_arr, i_t)
    complex(fp), intent(inout) :: psi_arr(:,:)
    integer(ip), intent(in) :: i_t

    call propagate_cn2d_splitop(psi_arr, i_t)

  end subroutine propagate_psi

  real(dp) function propagate_norm(psi_arr) result(val)
    complex(fp), intent(inout) :: psi_arr(:,:)
    integer(ip) :: i_x
    real(fp) :: norm_integrand(nx)

    val = numerics_trapz(norm_integrand(:), dx)
    !$omp parallel do private(i_x)
    do i_x = 1, nx
       norm_integrand(i_x) = numerics_trapz(abs(psi_arr(i_x, :))**2, dy)
    end do
    !$omp end parallel do
    val = numerics_trapz(norm_integrand(:), dx)
    
  end function propagate_norm
end module propagate
