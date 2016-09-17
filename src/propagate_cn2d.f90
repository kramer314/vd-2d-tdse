module propagate_cn2d

  include "./propagate_cn2d_src/moduse.src"

  implicit none

  private

  include "./propagate_cn2d_src/modvars.src"

  public :: propagate_cn2d_init
  public :: propagate_cn2d_cleanup
  public :: propagate_cn2d_splitop

contains
  subroutine propagate_cn2d_init()
    integer(ip) :: i_x, i_y

    dtau = dt

    include "./propagate_cn2d_src/init.src"

  end subroutine propagate_cn2d_init

  subroutine propagate_cn2d_cleanup()

    include "./propagate_cn2d_src/cleanup.src"

  end subroutine propagate_cn2d_cleanup

  subroutine propagate_cn2d_splitop(psi_arr, i_t)
    complex(dp), intent(inout) :: psi_arr(:,:)
    integer(ip), intent(in) :: i_t

    include "./propagate_cn2d_src/splitop.src"

  end subroutine propagate_cn2d_splitop

end module propagate_cn2d
