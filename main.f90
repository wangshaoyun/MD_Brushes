program main
use global_variables
use initialize_update
use compute_acceleration
use input_output
implicit none

  !#################data#################!
  integer :: i
  !################begin#################!
  call cpu_time(started)
  call random_seed()
  !
  !input and initialize system, timing and histogram parameters.
  call initialize_parameters
  !
  !initialize force and parameters of potential
  call initialize_force_parameters
  !
  !
  if (restart_or_continue == 0) then
    i=1
    !
    !initialize position and velocity
    call initialize_position
    call initialize_velocity
    call write_pos
    call write_pos1
    call write_vel1(1)
    !
    !error analysis
    if ( qq /= 0 ) then
      call error_analysis
    end if
    !
    !compute force
    call compute_force
    call write_acc
    call write_hist
  else if (restart_or_continue /= 0) then
    !
    !read position and histogram data
    call continue_read_data(i)
    !
    !initialize lj_verlet_list
    call lj_verlet_list
    !
    !Construct the real verlet list and real_point vector
    if ( real_verlet == 1 ) then
      call real_verlet_list
    end if
    !
    !error analysis
    if (qq /= 0 ) then
      call error_analysis
    end if
    call rescale_velocity
    call compute_force
    !
    !write data
    call write_pos
    call write_vel
    call write_acc
  end if

  !##############preheating##############!
  if ( i <= StepNum0 ) then
    do step=i, StepNum0
      call new_position
      !
      !Rescale velocity to avoid the break of chemical bonds
      if ( mod( step, 20 ) == 0 .and. step < 300000 ) then
        call rescale_velocity
      end if
      call update_verlet_list
      if ( mod(step,DeltaStep1) == 0 ) then
        call height
        call write_height(step)
      end if
      if ( mod(step,DeltaStep3) == 0 ) then
        call error_analysis
        call write_pos1
        call write_vel1(step)
      end if
    end do
    i=step
  end if

  !################running###############!
  do step=i,StepNum0+StepNum            
    call new_position
    call update_verlet_list
    if ( mod(step,DeltaStep1) == 0 ) then
      call height
      call write_height(step)
    end if
    if ( mod(step,DeltaStep2) == 0 ) then     
      call histogram
    end if
    if ( mod(step,DeltaStep3) == 0 ) then
      call write_pos1                   
      call write_vel1(step)
      call write_hist
    end if
  end do

  !##################end#################!
  call cpu_time(finished)
  total_time=finished-started+total_time
  call write_time(total_time)
  write(*,*) 'Finished!'
  
end program






