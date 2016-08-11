! Copyright (c) 2016 Alex Kramer <kramer.alex.kramer@gmail.com>
! See the LICENSE.txt file at the top-level of this distribution.

module propagate_gaussian
  use progvars
  use gaussian, only: gaussian_xyt, gaussian_init, gaussian_cleanup

  implicit none

  private

  public :: propagate_gaussian_psi
  public :: propagate_gaussian_init
  public :: propagate_gaussian_cleanup

contains

  subroutine propagate_gaussian_init()
    call gaussian_init()
  end subroutine propagate_gaussian_init

  subroutine propagate_gaussian_cleanup()
    call gaussian_cleanup()
  end subroutine propagate_gaussian_cleanup

  subroutine propagate_gaussian_psi(psi_arr, i_t)

    complex(fp), intent(inout) :: psi_arr(:,:)
    integer, intent(in) :: i_t

    real(fp) :: x, y, t
    integer(ip) :: i_x, i_y

    t = t_range(i_t)

    ! For testing efficiency, only fill VD region

    ! For parallel code, we assume that the number of grid points is much
    ! larger than the number of virtual detectors, and parallelize accordingly

    ! Efficiency improvements ->
    ! * Does improving stride access to psi_arr by splitting top / bottom
    !   loops beat overhead of putting more things on the stack, etc by
    !   having more omp calls?
    ! * Would assigning pointers to the correct points in psi_arr improve
    !   efficiency?

    ! Fill bottom array
    do i_y = vd_yl_min - 1, vd_yl_max + 1

       y = y_range(i_y)
       !$omp parallel do private(i_x, x) firstprivate(x_range)
       do i_x = vd_xl_min - 1, vd_xr_max + 1
          x = x_range(i_x)
          psi_arr(i_x, i_y) = gaussian_xyt(x, y, t)
       end do
       !$omp end parallel do
    end do

    ! Fill left/right arrays
    ! We loop in this way for optimal stride over psi_arr
    ! We only go from yl_max + 2 : yr_min - 2 so that we don't overlap

    !$omp parallel do private(i_x, i_y, x, y) firstprivate(x_range, y_range)
    do i_y = vd_yl_max + 2, vd_yr_min - 2
       y = y_range(i_y)
       ! Fill left array
       do i_x = vd_xl_min - 1, vd_xl_max + 1
          x = x_range(i_x)
          psi_arr(i_x, i_y) = gaussian_xyt(x, y, t)
       end do
    end do

    !$omp parallel do private(i_x, i_y, x, y) firstprivate(x_range, y_range)
    do i_y = vd_yl_max + 2, vd_yr_min - 2
       y = y_range(i_y)
       ! Fill right array
       do i_x = vd_xr_min - 1, vd_xr_max + 1
          x = x_range(i_x)
          psi_arr(i_x, i_y) = gaussian_xyt(x, y, t)
       end do
    end do
    !$omp end parallel do

    ! Fill top array
    do i_y = vd_yr_min - 1, vd_yr_max + 1

       y = y_range(i_y)
       !$omp parallel do private(i_x, x) firstprivate(x_range)
       do i_x = vd_xl_min - 1, vd_xr_max + 1
          x = x_range(i_x)
          psi_arr(i_x, i_y) = gaussian_xyt(x, y, t)
       end do
       !$omp end parallel do
    end do


  end subroutine propagate_gaussian_psi
end module propagate_gaussian
