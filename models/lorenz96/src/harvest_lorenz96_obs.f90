program lorenz96Model

  ! Get unit numbers for stdin, stdout, stderr in a portable way
  use, intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
                                            stdout=>output_unit, &
                                            stderr=>error_unit

  use lorenz96, only : lorenz96_type
  use kind, only     : r8kind, i8kind

  implicit none

  integer :: step, i, mthd           ! Loop index
!  integer :: obs_idx(15)             ! Array of observation indices
!  integer :: obs_tim(15)             ! Array of observation times
  integer :: obs_idx(1000)             ! Array of observation indices
  integer :: obs_tim(1000)             ! Array of observation times
  type(lorenz96_type) :: obs_data(3) ! Array of the observation data
  integer :: ierr                    ! Error code 
  character(len=128) :: filename
  integer :: fileunit

  integer :: j, k

  ! Set the indices of the obs locations
!  obs_idx=(/4, 6, 8, 11, 12, 14, 18, 19, 25, 27, 29, 30, 33, 35, 40/)
!  obs_tim=(/1, 2, 3,  2,  1,  2,  3,  2,  1,  2,  3,  2,  1,  2,  3/)

  j=1
  do i = 1, 4000, 4
    obs_idx(j) = i
    k=modulo(j, 4)
    if (k ==0 ) k = 2
    obs_tim(j) = k
    j = j + 1
  end do

  ! Loop over experiment analysis times
  do step = 40060, 80000, 60

     ! Read in truth at analysis time and surrounding times
     obs_data(1) = lorenz96_type(step - 10, "NETCDF")
     obs_data(2) = lorenz96_type(step,"NETCDF")
     obs_data(3) = lorenz96_type(step + 10, "NETCDF")

    do mthd = 1, 4
    
      ! Construct name of obs input file
      write(filename, '(A,I0.7,A,I1,A)') 'lorenz96obs_', step, '_', mthd, '.txt'

      ! Open the obs file
      open(newunit=fileunit, file=trim(filename), form='formatted')

     ! Make observation file for method 1
     write(fileunit,'(I)') size(obs_idx)

     select case (mthd)
       case(1)
         do i=1, size(obs_idx)
           write(fileunit, '(2I,F7.3)') 1, obs_idx(i), obs_data(obs_tim(i))%state(obs_idx(i))
         end do
       case(2)
         do i=1, size(obs_idx)
           write(fileunit, '(2I,F7.3)') 1, obs_idx(i), obs_data(2)%state(obs_idx(i))
         end do
       case(3)
         do i=1, size(obs_idx)
           write(fileunit, '(2I,F7.3)') obs_tim(i), obs_idx(i), obs_data(obs_tim(i))%state(obs_idx(i))
         end do
       case(4)
         do i=1, size(obs_idx)
           write(fileunit, '(2I,F7.3)') obs_tim(i), obs_idx(i), obs_data(obs_tim(i))%state(obs_idx(i))
         end do
      end select

      close(fileunit)

    end do

  end do

end program
