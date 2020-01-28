!=======================================================================================================================
!
!  *******************
!  * pk_wrtie_output *
!  *******************
!
!  Purpose:
!  --------
!  Module to write the results to the output file.
!
!  Usage:
!  ------
!  CALL pk_write_init       to initiate the module and create the header in the output file.
!  CALL pk_write_to_file    to write the current time step to the output file.
!  CALL pk_write_time_steps to write the time step sizes to the output file
!
!=======================================================================================================================
MODULE pk_write_output
!=======================================================================================================================
! Record of revisions:
!       Date       Programmer               Description of changes
!       -----------------------------------------------------------
!       18/03/19   K. Luszczek              Original code. Aggregates functions previously in separate modules.
!
!=======================================================================================================================
!  Import subroutines
   USE pk_read_input , ONLY: newunit
   USE pk_data_types , ONLY: step_size_fetch
   USE ifport        , ONLY: makedirqq
!  Import variables
   USE pk_kinds      , ONLY: dp
   USE pk_data_types , ONLY: solution_at_step, fs_time_stamps, fs_n, fs_rho, fs_c, pc_group, exec_time
   USE pk_enumeration, ONLY: PLOT_RHO, PLOT_N, PLOT_C
!
   IMPLICIT NONE
!
   PRIVATE
!
!  Public subroutines
!  ------------------
   PUBLIC :: pk_write_init    ! Subroutine which process the input files.
   PUBLIC :: pk_write_to_file ! Function returns the available unit.
   PUBLIC :: pk_plot_select   ! Subroutine to plot selected, time-dependent final solution variables.
   PUBLIC :: pk_write_time_steps ! Subroutube which writes time step sizes to a file.
!
!=======================================================================================================================
CONTAINS
!=======================================================================================================================
SUBROUTINE pk_write_init(&
   current_step )          ! IN
!=======================================================================================================================
! Record of revisions:
!       Date       Programmer               Description of changes
!       -----------------------------------------------------------
!       11/03/18   K. Luszczek              Original code
!
!=======================================================================================================================
!
!  Arguments
!
   TYPE(solution_at_step), INTENT(IN) :: current_step
!
!  Locals
!
   INTEGER  :: unit_out
   INTEGER  :: ierror
   REAL(dp) :: current_time = 0.0_dp
!
!  Open an output file to print to
!  -----------------------
   unit_out = newunit()
!
   OPEN (UNIT = unit_out, FILE = 'pk_output.dat', STATUS = 'REPLACE', ACTION = 'WRITE', IOSTAT = ierror)
   IF (ierror /= 0) THEN
      WRITE(*,*) 'An error occurred when opening the output file. IOSTAT error: ', ierror
      WRITE(*,*) 'Execution stops.'
      STOP
   ENDIF
!  
   WRITE(unit_out,*, IOSTAT = ierror) ('      time(s)       rho($)   rel.n.dens.  rel.precursor con. ')
   IF (ierror /= 0) THEN
      WRITE(*,*) 'Error occurred when writing to output file. IOSTAT error: ', ierror
      WRITE(*,*) 'Execution stops.'
      STOP
   ENDIF
!
   WRITE(unit_out,110, IOSTAT = ierror) current_time, current_step%rho, current_step%n2p, current_step%c
   IF (ierror /= 0) THEN
      WRITE(*,*) 'An error occurred when writing to the output file. IOSTAT error: ', ierror
      WRITE(*,*) 'Execution stops.'
      STOP
   ENDIF
!
   110 FORMAT ( ES14.6, ES14.6, ES14.6, ES14.6, ES14.6,ES14.6, ES14.6, ES14.6, ES14.6 )
!
   CLOSE (UNIT = unit_out)
!
END SUBROUTINE pk_write_init
!=======================================================================================================================
!
!  ********************
!  * pk_write_to_file *
!  ********************
!
!  Purpose: Writes basic output to an output_data.dat file:
!           time stamps, rel. n. density, reactivity, delayed neutron precursor concnetration
!
!=======================================================================================================================
SUBROUTINE pk_write_to_file()
!=======================================================================================================================
! Record of revisions:
!       Date       Programmer               Description of changes
!       -----------------------------------------------------------
!       11/03/18   K. Luszczek              Original code
!
!=======================================================================================================================
!
!  Locals
!
   CHARACTER(120) :: output_file ! Output file name
   INTEGER :: unit_out
   INTEGER :: ierror
   INTEGER :: idx
   LOGICAL :: res ! Stores the result of makedirqq
!
   output_file = 'output/output_file.dat'
   res = makedirqq('output')
!
!  Open an output file to print to
!  -------------------------------
   unit_out = newunit()
!
   OPEN (UNIT=unit_out, FILE = output_file, STATUS = 'REPLACE', ACTION = 'WRITE', IOSTAT = ierror)
   IF (ierror /= 0) THEN
      WRITE(*,*) 'An error occurred when opening the output file. IOSTAT error: ', ierror
      WRITE(*,*) 'Execution stops.'
      STOP
   ENDIF
!  
   110 FORMAT ( ES14.6, ES14.6, ES14.6, ES14.6, ES14.6,ES14.6, ES14.6, ES14.6, ES14.6 )
!
!  Loop over the final solution arrays
!
   DO idx = 1, SIZE(fs_time_stamps)
      WRITE(unit_out,110, IOSTAT = ierror) fs_time_stamps(idx), fs_rho(idx), fs_n(idx), fs_c(idx,1:pc_group)
      IF (ierror /= 0) THEN
         WRITE(*,*) 'An error occurred when writing to the output file. IOSTAT error: ', ierror
         WRITE(*,*) 'Execution stops.'
         STOP
      ENDIF
   ENDDO
!
   CLOSE (UNIT = unit_out)
!
END SUBROUTINE pk_write_to_file
!=======================================================================================================================
!
!  ***********************
!  * pk_write_time_steps *
!  ***********************
!
!  Purpose: Writes time step sizes to step_sizes.dat file:
!
!=======================================================================================================================
SUBROUTINE pk_write_time_steps()
!=======================================================================================================================
! Record of revisions:
!       Date       Programmer               Description of changes
!       -----------------------------------------------------------
!       11/03/18   K. Luszczek              Original code
!
!=======================================================================================================================
!
!  Locals
!
   CHARACTER(120) :: output_file ! Output file name
   INTEGER  :: unit_out
   INTEGER  :: ierror
   INTEGER  :: idx
   LOGICAL  :: res ! Stores the result of makedirqq
   INTEGER  :: step_number
   REAL(dp) :: time_at_eos_old
   REAL(dp) :: time_at_eos_new
   REAL(dp) :: delta_t
   REAL(dp) :: omega
   REAL(dp) :: omega_rt
   REAL(dp) :: n_eos_rt
   LOGICAL  :: next_available      = .TRUE.
   INTEGER  :: implicit_iter_count = 0
   INTEGER  :: iteration_total     = 0
   INTEGER  :: time_reject_count   = 0
   INTEGER  :: rejected_total      = 0
   INTEGER  :: pos_in_file = 0 ! Stores position in the file
!
   output_file = 'output/step_sizes.dat'
   res = makedirqq('output')
!
!  Open an output file to print to
!  -------------------------------
   unit_out = newunit()
!
   OPEN (UNIT=unit_out, FILE = output_file, STATUS = 'REPLACE', ACTION = 'WRITE', IOSTAT = ierror)
   !ACTION = 'WRITE', IOSTAT = ierror)
   IF (ierror /= 0) THEN
      WRITE(*,*) 'An error occurred when opening the output file. IOSTAT error: ', ierror
      WRITE(*,*) 'Execution stops.'
      STOP
   ENDIF
!
!  Step number, Total time at end-of-step, Delta_t
   110 FORMAT ( I14, ES14.6, ES14.6, I14, I14, ES14.6, ES14.6, ES14.6)
!
!  Read the first entry from the step_size_store list
   CALL step_size_fetch(                         &
      step_number         = step_number        , & ! OUT
      time_at_eos         = time_at_eos_new    , & ! OUT
      time_reject_count   = time_reject_count  , & ! OUT
      implicit_iter_count = implicit_iter_count, & ! OUT
      next_available      = next_available     , & ! OUT
      omega               = omega              , & ! OUT
      omega_rt            = omega_rt           , & ! OUT
      n_eos_rt            = n_eos_rt )             ! OUT
!
!  Loop until there unread enties in the step_size_store list
!
   DO
      IF (.NOT. next_available) THEN
         EXIT
      ELSE
         pos_in_file =  pos_in_file + 1
      ENDIF
!
      time_at_eos_old = time_at_eos_new
!
      CALL step_size_fetch(                         &
         step_number         = step_number        , & ! OUT
         time_at_eos         = time_at_eos_new    , & ! OUT
         time_reject_count   = time_reject_count  , & ! OUT
         implicit_iter_count = implicit_iter_count, & ! OUT
         next_available      = next_available     , & ! OUT
         omega               = omega              , & ! OUT
         omega_rt            = omega_rt           , & ! OUT
         n_eos_rt            = n_eos_rt )             ! OUT
!
      delta_t = time_at_eos_new - time_at_eos_old
      rejected_total = rejected_total + time_reject_count
      iteration_total = iteration_total + implicit_iter_count
!
      WRITE(unit_out,110, IOSTAT = ierror) step_number      , time_at_eos_new, delta_t , implicit_iter_count, &
                                           time_reject_count, omega          , omega_rt, n_eos_rt
      IF (ierror /= 0) THEN
         WRITE(*,*) 'An error occurred when writing to the time step file. IOSTAT error: ', ierror
         WRITE(*,*) 'Execution stops.'
         STOP
      ENDIF
   END DO
!
!  Rewind to the first line and overwrite with some useful data
!
   WRITE(unit_out,*, IOSTAT = ierror) 'Average time step =', time_at_eos_new/step_number
   IF (ierror /= 0) THEN
      WRITE(*,*) 'An error occurred when writing to the output file. IOSTAT error: ', ierror
      WRITE(*,*) 'Execution stops.'
      STOP
   ENDIF
   WRITE(unit_out,*, IOSTAT = ierror) 'Main loop exec time =', exec_time
   IF (ierror /= 0) THEN
      WRITE(*,*) 'An error occurred when writing to the output file. IOSTAT error: ', ierror
      WRITE(*,*) 'Execution stops.'
      STOP
   ENDIF
   WRITE(unit_out,*, IOSTAT = ierror) 'Number of rejected steps =', rejected_total
   IF (ierror /= 0) THEN
      WRITE(*,*) 'An error occurred when writing to the output file. IOSTAT error: ', ierror
      WRITE(*,*) 'Execution stops.'
      STOP
   ENDIF
   WRITE(unit_out,*, IOSTAT = ierror) 'Number of iterations =', iteration_total
   IF (ierror /= 0) THEN
      WRITE(*,*) 'An error occurred when writing to the output file. IOSTAT error: ', ierror
      WRITE(*,*) 'Execution stops.'
      STOP
   ENDIF
   CLOSE (UNIT = unit_out)
!
END SUBROUTINE pk_write_time_steps
!=======================================================================================================================
!
!  *****************
!  * pk_plot_n_rho *
!  *****************
!
!  Purpose: Creates a temporary GNU plot file and calls GNUPLOT to create a plot showing relative neutron density
!           and reactivity versus time.
!
!=======================================================================================================================
SUBROUTINE pk_plot_n_rho()
!=======================================================================================================================
! Record of revisions:
!       Date       Programmer               Description of changes
!       -----------------------------------------------------------
!       31/05/19   K. Luszczek              Original code
!
!=======================================================================================================================
!
!  Locals
!
   CHARACTER(120) :: axis1_title
   CHARACTER(120) :: axis2_title
   CHARACTER(120) :: plot_name
   INTEGER :: unit_out ! stores a temporary gnuplot file unit number.
   INTEGER :: ierror
!
!  Set parameters
!
   plot_name = 'plot_N_RHO'
   axis1_title = 'relative neutron density (-)'
   axis2_title = 'reactivity ($)'
!
!  Write the temporary gnuplot file
!
   unit_out = newunit()
   OPEN (UNIT=unit_out, FILE = 'temp_plt', STATUS = 'REPLACE', ACTION = 'WRITE', IOSTAT = ierror)
   IF (ierror /= 0) THEN
      WRITE(*,*) 'An error occurred when opening the output file. IOSTAT error: ', ierror
      WRITE(*,*) 'Execution stops.'
      STOP
   ENDIF
!
!  Write to the dummy temrinal to get y1 and y2 axes bounds
!
   WRITE(unit_out,'(A)') 'set terminal unknown'
   WRITE(unit_out,'(A)') 'set xlabel "time (s)"'
   WRITE(unit_out,'(A)') 'set ylabel "'//trim(axis1_title)//'"'
   WRITE(unit_out,'(A)') 'set autoscale y'
   WRITE(unit_out,'(A)') 'set tics out'
   WRITE(unit_out,'(A)') 'set format x "%.2e"'
   WRITE(unit_out,'(A)') 'set format y "%.2e"'
   WRITE(unit_out,'(A)') 'set ytics nomirror'
   WRITE(unit_out,'(A)') 'set y2tics nomirror'
   WRITE(unit_out,'(A)') 'set y2label "'//trim(axis2_title)//'"'
   WRITE(unit_out,'(A)') 'set format y2 "%.2e"'
   WRITE(unit_out,'(A)') 'set autoscale y2'
   WRITE(unit_out,'(A)') 'plot \'
!
   WRITE(unit_out,'(A)') '"output/output_file.dat" using 1:3 axes x1y1'//&
                         ' with lines lt 1 lw 3'//' title "'//trim(axis1_title)//'"'//',\'
!
   WRITE(unit_out,'(A)') '"output/output_file.dat" using 1:2 axes x1y2'//&
                         ' with lines lt 2 lw 3'//' title "'//trim(axis2_title)//'"'
!
!  Get the y1 and y2 bounds for aligning y1 and y2 tics
!
   WRITE(unit_out,'(A)') 'min_y1 = GPVAL_Y_MIN'
   WRITE(unit_out,'(A)') 'max_y1 = GPVAL_Y_MAX'
!
   WRITE(unit_out,'(A)') 'min_y2 = GPVAL_Y2_MIN-0.1'
   WRITE(unit_out,'(A)') 'max_y2 = GPVAL_Y2_MAX+0.1'
   WRITE(unit_out,'(A)') 'if (abs(min_y2) > abs(max_y2)) max_y2 = abs(min_y2); else min_y2 = -abs(max_y2)'
!
   WRITE(unit_out,'(A)') 'N = 6'
   WRITE(unit_out,'(A)') 'dy1 = (max_y1 - min_y1) / N'
   WRITE(unit_out,'(A)') 'dy2 = (max_y2 - min_y2) / N'
   WRITE(unit_out,'(A)') 'set ytics min_y1, dy1'
   WRITE(unit_out,'(A)') 'set y2tics min_y2, dy2'
   WRITE(unit_out,'(A)') 'set yr [min_y1:max_y1]'
   WRITE(unit_out,'(A)') 'set y2r [min_y2:max_y2]'
   WRITE(unit_out,'(A)') 'set grid'
!
!  Replot
!
   WRITE(unit_out,'(A)') "set terminal png size 1224,768 font '/usr/share/fonts/truetype/FreeMono.ttf' 15"
   WRITE(unit_out,'(A)') "set tmargin 1"
   WRITE(unit_out,'(A)') 'set output "output/'//trim(plot_name)//'.png"'
   WRITE(unit_out,'(A)') 'replot'
!
!  Call gnuplot
!
   CALL system('gnuplot temp_plt')
! Erase the temp file
   IF (ierror == 0) CLOSE(UNIT=unit_out)
!
END SUBROUTINE pk_plot_n_rho
!=======================================================================================================================
!
!  *************
!  * pk_plot_c *
!  *************
!
!  Purpose: Creates a temporary GNU plot file and calls GNUPLOT to create a plot showing precursor densities over time.
!
!=======================================================================================================================
SUBROUTINE pk_plot_c()
!=======================================================================================================================
! Record of revisions:
!       Date       Programmer               Description of changes
!       -----------------------------------------------------------
!       27/05/19   K. Luszczek              Original code
!
!=======================================================================================================================
!
!  Locals
!
   CHARACTER(120) :: axis1_title
   CHARACTER(120) :: axis2_title
   CHARACTER(120) :: plot_name
   CHARACTER(120) :: line_title
   CHARACTER(120) :: string_to_write
   INTEGER :: unit_out ! stores a temporary gnuplot file unit number.
   INTEGER :: ierror
   INTEGER :: idy      ! Index variable.
   INTEGER :: x        ! Column in the file output_data.dat to write if different than y1 and y2.
!
!  Set parameters
!
   plot_name = 'plot_C'
   axis1_title = 'delayed neutron precusor concentrations'
!
!  Write the temporary gnuplot file
!
   unit_out = newunit()
   OPEN (UNIT=unit_out, FILE = 'temp_plt', STATUS = 'REPLACE', ACTION = 'WRITE', IOSTAT = ierror)
   IF (ierror /= 0) THEN
      WRITE(*,*) 'An error occurred when opening the output file. IOSTAT error: ', ierror
      WRITE(*,*) 'Execution stops.'
      STOP
   ENDIF
!
   WRITE(unit_out,'(A)') 'set terminal png size 1024,768'
   WRITE(unit_out,'(A)') 'set rmargin 10'
   WRITE(unit_out,'(A)') 'set xlabel "time (s)"'
   WRITE(unit_out,'(A)') 'set ylabel "'//trim(axis1_title)//'"'
   WRITE(unit_out,'(A)') 'set autoscale y'
   WRITE(unit_out,'(A)') 'set tics out'
   WRITE(unit_out,'(A)') 'set format x "%.2t*10^%+03T"'
   WRITE(unit_out,'(A)') 'set format y "%.2t*10^%+03T"'
   WRITE(unit_out,'(A)') 'set grid'
!
   WRITE(unit_out,'(A)') 'set output "output/'//trim(plot_name)//'.png"'
   WRITE(unit_out,'(A)') 'plot \'
   DO idy = 1, pc_group
      x = 3 + idy
      WRITE(string_to_write,'(I2)') x-3
      line_title = 'delayed neutron precusor concentration - family '//trim(adjustl(string_to_write))
      IF (idy == pc_group) THEN
         WRITE(unit_out,'(A)') '"output/output_file.dat" using 1:'//trim(adjustl(string_to_write))//' axes x1y1'//&
                               ' with lines'//' title "'//trim(line_title)//'"'
      ELSE
         WRITE(unit_out,'(A)') '"output/output_file.dat" using 1:'//trim(adjustl(string_to_write))//' axes x1y1'//&
                               ' with lines'//' title "'//trim(line_title)//'"'//',\'
      ENDIF
   END DO
!
!  Call gnuplot
!
   CALL system('gnuplot temp_plt')
! Erase the temp file
   IF (ierror == 0) CLOSE(UNIT=unit_out, status='delete')
END SUBROUTINE pk_plot_c
!=======================================================================================================================
!
!  ******************
!  * pk_plot_select *
!  ******************
!
!  Purpose: Matches the requested data with enumeration.
!
!=======================================================================================================================
SUBROUTINE pk_plot_select(&
   plot ) ! IN
!=======================================================================================================================
! Record of revisions:
!       Date       Programmer               Description of changes
!       -----------------------------------------------------------
!       27/05/19   K. Luszczek              Original code
!
!=======================================================================================================================
!
!  Arguments
!
   CHARACTER(*), INTENT(IN):: plot ! Type of plot to be created.
!
!
!
   SELECT CASE(plot)
   CASE('PLOT_N_RHO')
      CALL pk_plot_n_rho()
   CASE('PLOT_C')
      CALL pk_plot_c()
   CASE DEFAULT
      WRITE(*,*) 'No plot selected for plotting.'
   END SELECT
!
END SUBROUTINE pk_plot_select
!=======================================================================================================================
END MODULE pk_write_output
!=======================================================================================================================