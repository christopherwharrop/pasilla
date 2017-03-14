program QG

  ! Get unit numbers for stdin, stdout, stderr in a portable way
  use, intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
                                            stdout=>output_unit, &
                                            stderr=>error_unit

  use kind
  use QG_Config
  use QG_Model
  use QG_Writer
  use gptl

  implicit none

  integer :: start_step = 0
  integer :: spinup_steps = 720
  integer :: run_steps = 1440
  integer :: output_interval_steps = 3
  logical :: readstart = .true.

  namelist /runtime/ start_step, spinup_steps, run_steps, output_interval_steps, readstart

  type(qg_config_type) :: config
  type(qg_model_type)  :: model
  type(qg_writer_type) :: writer
  character(len=64)    :: filename  ! output/input filename

  integer step, ret

  ret = gptlsetutr (gptlnanotime)
  ret = gptlinitialize ()
  ret = gptlstart ('QG')

  ! Read namelist from stdin
  read(stdin,nml=runtime)
  rewind(stdin)

  config = qg_config_type(stdin)
  model = qg_model_type(config)
  writer = qg_writer_type('NETCDF')

  ! Spinup the model (if required)
  write(*,*) 'Integrating model spinup steps: ', spinup_steps
  call model%forward(spinup_steps)

  ! Output fields derived from initial state
  write(filename,'(A,I0.7)') 'qgout_', model%get_step()
  call writer%write(model, filename)

  ! Run the model
  write(*,*) 'Integrating model trajectory steps: ', run_steps
  do step = 1, run_steps, output_interval_steps

    ! Advance the model to next output interval
    call model%forward(min(output_interval_steps, run_steps))

    ! Output fields derived from current model state
    if (output_interval_steps <= run_steps) then
      write(filename,'(A,I0.7)') 'qgout_', model%get_step()
      call writer%write(model, filename)
    end if

  enddo

  ret = gptlstop ('QG') 
  ret = gptlpr (0) 
  ret = gptlfinalize ()

 
end program QG