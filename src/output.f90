! Copyright (c) 2016 Alex Kramer <kramer.alex.kramer@gmail.com>
! See the LICENSE.txt file at the top-level of this distribution.

module output
  ! Output module

  ! Imports -- library dependencies
  use log, only: log_log_info, log_log_critical, log_stderr
  use string, only: string_val

  ! Imports -- program modules
  use progvars

  implicit none

  private

  public :: output_init
  public :: output_cleanup

  public :: output_vd_counts
  public :: output_vd_residuals

  public :: output_logfile_unit

  ! Output file unit numbers
  integer(ip), parameter :: logfile_unit = 99
  integer(ip), parameter :: output_logfile_unit = logfile_unit

  integer(ip), parameter :: vd_p_unit = 98
  integer(ip), parameter :: vd_resid_unit = 97
  integer(ip), parameter :: vd_resid_cum_unit = 96
  integer(ip), parameter :: vd_resid_analysis_unit = 95


contains
  subroutine output_init()
    ! Initialize output module, including opening output files

    open(unit=logfile_unit, file=trim(output_dir)//trim(log_fname))
    open(unit=vd_p_unit, file=trim(output_dir)//trim(vd_p_fname))
    open(unit=vd_resid_unit, file=trim(output_dir)//trim(vd_resid_fname))
    open(unit=vd_resid_cum_unit, &
         file=trim(output_dir)//trim(vd_resid_cum_fname))
    open(unit=vd_resid_analysis_unit, file=trim(output_dir)//&
         trim(vd_resid_analysis_fname))

  end subroutine output_init

  subroutine output_cleanup()
    ! Cleanup output module, including closing output files

    close(unit=logfile_unit)
    close(unit=vd_p_unit)
    close(unit=vd_resid_unit)
    close(unit=vd_resid_cum_unit)
    close(unit=vd_resid_analysis_unit)

  end subroutine output_cleanup

  subroutine output_vd_counts()
    ! Output virtual detector momentum-binning counts
    integer(ip) :: i_py

    call log_log_info("Writing out VD counts", logfile_unit)

    do i_py = 1, vd_npy
       write(vd_p_unit, "("//string_val(vd_npx)//fp_format_raw//")") &
            vd_xy%vd_p_arr(:, i_py)
    end do

  end subroutine output_vd_counts

  subroutine output_vd_residuals()
    ! Output virtual detector residual comparisons
    integer(ip) :: i_py

    call log_log_info("Writing out VD residuals", logfile_unit)

    do i_py = 1, vd_npy
       write(vd_resid_unit, "("//string_val(vd_npx)//fp_format_raw//")") &
            resid_np_arr(:, i_py)
       write(vd_resid_cum_unit, "("//string_val(vd_npx)//fp_format_raw//")") &
            resid_np_cum_arr(:, i_py)
    end do

  end subroutine output_vd_residuals

end module output
