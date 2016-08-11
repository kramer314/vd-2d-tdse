! Copyright (c) 2016 Alex Kramer <kramer.alex.kramer@gmail.com>
! See the LICENSE.txt file at the top-level of this distribution.

program main
  ! Main executable

  ! Imports -- library dependencies
  use log, only: log_log_info
  use string, only: string_val

  ! Imports -- program modules
  use progvars
  use setup, only: setup_init, setup_cleanup
  use output, only: output_vd_counts, output_vd_residuals, &
       logfile_unit=>output_logfile_unit
  use propagate, only: propagate_psi
  ! use resids, only: resids_calculate

  implicit none

  integer(ip) :: i_t

  call setup_init()

  call log_log_info("Beginning time propagation.", logfile_unit)
  do i_t = 1, nt
     
     call propagate_psi(psi_arr, i_t)
     call vd_xy%update(psi_arr)
     !call log_log_info("Done updating VDs", logfile_unit)

     if (mod(i_t, print_mod_t) .eq. 0) then
        call log_log_info("Timestep "//string_val(i_t)//" of "// &
             string_val(nt), logfile_unit)
     end if
     

  end do
  call log_log_info("Time propagation complete.", logfile_unit)

  call vd_xy%finalize()
  call output_vd_counts()

  ! call resids_calculate()
  ! call resids_stats()

  ! call output_vd_residuals()

  call setup_cleanup()
end program main