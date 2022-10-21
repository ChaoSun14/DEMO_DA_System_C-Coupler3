module model_setting_mod

    !use mpi
    use parse_namelist_mod
    use spmd_init_mod
    use coupling_ocn_model_mod
    use grid_init_mod
    use decomp_init_mod
    use variable_mod
    use ran_mod

    implicit none

    real(kind=4), public, allocatable :: ran_num(:)

contains

  subroutine atm_demo_init

      implicit none
      integer :: mpicom

      mpicom = CCPL_NULL_COMM

      call register_atm_demo_component(mpicom)

      call parse_namelist
      call spmd_init(mpicom)
      call grid_init
      call decomp_init
      call variable_init

      call register_component_coupling_configuration

  end subroutine atm_demo_init

  subroutine atm_demo_step_on

      implicit none

      integer    :: i

      allocate(ran_num(decomp_size(mytask_id+1)))

      do i=1, decomp_size(mytask_id+1)
        ran_num(i) = normal(0., 2.)
      end do

      do i=1,time_length/time_step
          
          TAUX_m  = TAUX_l
          TAUY_m  = TAUY_l
          PS_m = PS_l
          t_m   = T_l
         
      !    TAUX_m  = TAUX_l+ran_num
      !    TAUY_m  = TAUY_l+ran_num
      !    PS_m = PS_l+ran_num
      !   do j=1, levlen 
      !    T_m(:,j)  = T_l(:,j)+ran_num(:)
      !   end do
         
          call execute_atm_coupling

      end do

  end subroutine atm_demo_step_on

  subroutine finalize_atm_demo

      implicit none
      integer :: ierr

      CALL MPI_BARRIER(ensemble_comm,ierr)
      call deallocate_fields
      if(allocated(lat)) deallocate(lat)
      if(allocated(lon)) deallocate(lon)
      if(allocated(lev)) deallocate(lev)
      if(allocated(local_grid_cell_index)) deallocate(local_grid_cell_index)

      call finalize_atm_coupling
      call mpi_finalize(ier)

  end subroutine finalize_atm_demo

end module model_setting_mod
