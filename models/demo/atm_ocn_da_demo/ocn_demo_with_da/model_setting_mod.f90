module model_setting_mod

    !use mpi
    use parse_namelist_mod
    use CCPL_interface_mod
    use spmd_init_mod
    use coupling_atm_model_mod
    use grid_init_mod
    use decomp_init_mod
    use variable_mod

    implicit none

contains

  subroutine ocn_demo_init

      implicit none
      integer :: mpicom

      mpicom = CCPL_NULL_COMM

      call register_ocn_demo_component(mpicom)

      call parse_namelist
      call spmd_init(mpicom)
      call grid_init
      call decomp_init
      call variable_init

      call register_component_coupling_configuration

  end subroutine ocn_demo_init

  subroutine ocn_demo_step_on

      implicit none

      integer    :: i, j
      logical    :: interface_status
      
      do i=1,time_length/time_step
          
          sstm = sst_l
          sshm = ssh_l
          saltm = salt_l
          saltm_1 = salt_l_1

          do j=1, num_3d_vars
            vars_3d(:,:,j) = saltm_1
          end do
          
          call execute_ocn_coupling

      end do

  end subroutine ocn_demo_step_on

  subroutine finalize_ocn_demo

      implicit none
      integer :: ierr

      CALL MPI_BARRIER(ensemble_comm,ierr)

      call deallocate_fields
      if(allocated(lat)) deallocate(lat)
      if(allocated(lon)) deallocate(lon)
      if(allocated(local_grid_cell_index)) deallocate(local_grid_cell_index)

      call finalize_ocn_coupling
      call mpi_finalize(ier)

  end subroutine finalize_ocn_demo

end module model_setting_mod
