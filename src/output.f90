! Copyright (c) 2016 Alex Kramer <kramer.alex.kramer@gmail.com>
! See the LICENSE.txt file at the top-level of this distribution.

module output
  ! Output module

  ! Imports -- library dependencies
  use log, only: log_log_info, log_log_critical, log_stderr
  use string, only: string_val

  ! Imports -- program modules
  use progvars
  use potential, only: potential_xyt

  implicit none

  private

  public :: output_init
  public :: output_cleanup

  public :: output_vd_counts
  public :: output_psi
  public :: output_potential

  public :: output_logfile_unit

  ! Output file unit numbers
  integer(ip), parameter :: logfile_unit = 99
  integer(ip), parameter :: output_logfile_unit = logfile_unit

  integer(ip), parameter :: vd_p_unit = 98

contains
  subroutine output_init()
    ! Initialize output module, including opening output files

    open(unit=logfile_unit, file=trim(output_dir)//trim(log_fname))
    open(unit=vd_p_unit, file=trim(output_dir)//trim(vd_p_fname))

  end subroutine output_init

  subroutine output_cleanup()
    ! Cleanup output module, including closing output files

    close(unit=logfile_unit)
    close(unit=vd_p_unit)

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

  subroutine output_psi(fname, out_dir)
    ! Output abs(psi(x, y, t))**2 at the current time index

    character(:), allocatable, intent(in) :: fname
    character(:), allocatable, intent(in), optional :: out_dir

    character(:), allocatable :: path

    integer(ip) :: unit
    integer(ip) :: i_x, i_y

    if (present(out_dir)) then
       path = trim(out_dir)//trim(fname)
    else
       path = trim(output_dir)//trim(fname)
    end if

    open(newunit=unit, file=path)

    call log_log_info("Writing out abs(psi(x, y, t))^2 to "//path, &
         logfile_unit)

    do i_y = 1, ny / print_mod_y
       do i_x = 1, nx / print_mod_x
          write(unit, "("//string_val(nx)//fp_format_raw//")", advance="no") &
               abs(psi_arr(i_x * print_mod_x, i_y * print_mod_y))**2
       end do
       write(unit, *)
    end do

    close(unit)

  end subroutine output_psi

  subroutine output_potential(i_t, fname, out_dir)
    ! Output abs(psi(x, y, t))**2 at the current time index

    integer, intent(in) :: i_t
    character(:), allocatable, intent(in) :: fname
    character(:), allocatable, intent(in), optional :: out_dir

    character(:), allocatable :: path

    integer(ip) :: unit
    integer(ip) :: i_x, i_y

    if (present(out_dir)) then
       path = trim(out_dir)//trim(fname)
    else
       path = trim(output_dir)//trim(fname)
    end if

    open(newunit=unit, file=path)

    call log_log_info("Writing out V(x, y, t) to "//path, logfile_unit)

    do i_y = 1, ny / print_mod_y
       do i_x = 1, nx / print_mod_x
          write(unit, "("//string_val(nx)//fp_format_raw//")", advance="no") &
               real(potential_xyt(i_x * print_mod_x, i_y * print_mod_y, i_t))
       end do
       write(unit, *)
    end do

    close(unit)

  end subroutine output_potential

end module output
