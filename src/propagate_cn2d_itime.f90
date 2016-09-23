module propagate_cn2d_itime

  include "./propagate_cn2d_src/moduse.src"
  use wfmath, only: wfmath_proj

  implicit none

  private

  include "./propagate_cn2d_src/modvars.src"

  public :: propagate_cn2d_itime_init
  public :: propagate_cn2d_itime_cleanup
  public :: propagate_cn2d_itime_splitop
  public :: propagate_cn2d_itime_converged

contains
  subroutine propagate_cn2d_itime_init()
    integer(ip) :: i_x, i_y

    dtau = -j * dt

    include "./propagate_cn2d_src/init.src"

  end subroutine propagate_cn2d_itime_init

  subroutine propagate_cn2d_itime_cleanup()

    include "./propagate_cn2d_src/cleanup.src"
  end subroutine propagate_cn2d_itime_cleanup

  subroutine propagate_cn2d_itime_splitop(psi_arr)
    complex(dp), intent(inout) :: psi_arr(:,:)

    integer(ip), parameter :: i_t = 0_ip

    include "./propagate_cn2d_src/splitop.src"

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
