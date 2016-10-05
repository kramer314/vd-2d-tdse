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
  use output, only: output_vd_counts, output_psi, output_potential, &
       logfile_unit=>output_logfile_unit
  use propagate, only: propagate_psi, propagate_psi_itime
  use propagate_cn2d_itime, only: propagate_cn2d_itime_converged

  use wfmath, only: wfmath_norm, wfmath_normalize

  implicit none

  integer(ip) :: i_t
  logical :: converged

  call setup_init()

  call log_log_info("Beginning imaginary time propagation.", logfile_unit)

  converged = .false.
  i_t = 1_ip

  !$omp parallel workshare
  old_psi_arr(:,:) = psi_arr(:,:)
  !$omp end parallel workshare

  do while (.not. converged)

     call propagate_psi_itime(psi_arr)

     if (mod(i_t, print_mod_t) .eq. 0) then

        call log_log_info("Imaginary timestep "//string_val(i_t), logfile_unit)

        call wfmath_normalize(psi_arr, dgrid)

        if (i_t .gt. min_iter_itime) then

           converged = propagate_cn2d_itime_converged(psi_arr, old_psi_arr, epsilon(1.0_fp))

           !$omp parallel workshare
           old_psi_arr(:,:) = psi_arr(:,:)
           !$omp end parallel workshare
        end if

     end if

     i_t = i_t + 1

  end do

  call output_psi(psi_ground_fname)

  call log_log_info("Imaginary time propagation complete.", logfile_unit)

  call log_log_info("Beginning time propagation.", logfile_unit)
  do i_t = 1, nt

     call propagate_psi(psi_arr, i_t)
     call vd_xy%update(psi_arr)

     if (mod(i_t, print_mod_t) .eq. 0) then
        call log_log_info("Timestep "//string_val(i_t)//" of "// &
             string_val(nt), logfile_unit)
     end if

  end do
  call log_log_info("Time propagation complete.", logfile_unit)

  call vd_xy%finalize()
  call output_vd_counts()

  call setup_cleanup()
end program main
