program L96

  ! Get unit numbers for stdin, stdout, stderr in a portable way
  use, intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
                                            stdout=>output_unit, &
                                            stderr=>error_unit

  use kind,      only  : r8kind, i8kind
  use L96_Model, only  : l96_model_type
  use L96_Config, only : l96_config_type
  use L96_Reader, only : l96_reader_type
  use L96_Writer, only : l96_writer_type

  implicit none

  integer          :: start_step = 0
  integer          :: run_steps   = 1000
  integer          :: output_interval_steps = 100
  character(len=6) :: io_format = 'NETCDF'
  namelist /runtime/ start_step, run_steps, output_interval_steps, io_format

  type(L96_config_type) :: config ! Lorenz96 configuration
  type(L96_model_type)  :: model  ! Lorenz96 model
  type(L96_reader_type) :: reader ! Lorenz96 reader
  type(L96_writer_type) :: writer ! Lorenz96 writer

  integer :: t     ! Loop index
  integer :: ierr  ! Error code 
  character(len=64) :: filename  ! output/input filename
  character(len=16) :: format_extension

  ! Read namelist from stdin
  read(stdin,nml=runtime)
  rewind(stdin)

  ! If this is a restart, model from a restart file
  if (start_step /= 0) then

    ! Initialize a reader to read restart file
    reader = l96_reader_type(io_format)

    ! Construct name of restart file
    write(filename,'(A,I0.7)') 'lorenz96out_', start_step

    ! Load restart file into new model instance
    call reader%read(model, filename)

  ! Otherwise instantiate a new model
  else

    ! Create a Lorenz96 model configuration as specified in namelist
    config = l96_config_type(stdin)

    ! Create and configure a Lorenz96 model with the specified configuration
    model = l96_model_type(config)

  end if

!  call model%print()

  ! Construct name of output file
  write(filename,'(A,I0.7)') 'lorenz96out_', model%get_step()

  ! Create a writer to write model state
  writer = l96_writer_type(io_format)

  ! Write out model state if needed
  if (output_interval_steps <= run_steps) call writer%write(model, filename)

  ! Run the model
  do t = 1, run_steps, output_interval_steps

    ! Advance the model to next output interval
    call model%adv_nsteps(min(output_interval_steps, run_steps))

    ! Write out model state if needed
    if (output_interval_steps <= run_steps) then

      ! Construct name of output file
      write(filename,'(A,I0.7)') 'lorenz96out_', model%get_step()
      call writer%write(model, filename)

    end if

!    call model%print()

  end do

end program
