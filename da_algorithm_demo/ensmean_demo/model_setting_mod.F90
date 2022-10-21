module model_setting_mod_1

    use parse_namelist_mod_1
    use spmd_init_mod_1
    use grid_init_mod_1
    use decomp_init_mod_1
    use variable_mod_1
    use ran_mod_1
#ifdef CCPL_DA
    use da_ensmean_demo_ccpl_coupling_mod
#endif
    implicit none

contains

  subroutine da_algorithm_demo_init

      use mpi
      implicit none
      integer :: mpicom
      integer                    :: start_time, end_time, ierr
      real                       :: total_time, max_total_time, max_total_time1, total_time2, max_total_time2, total_time3, total_time4, max_total_time3, max_total_time4
      
#ifndef CCPL_DA
      mpicom = MPI_COMM_WORLD
      call spmd_init_da(mpicom)
#endif

      call parse_namelist_da      
#ifdef CCPL_DA
      call system_clock(start_time)
      call get_grids_info
      nlat = ccpl_nlat
      nlon = ccpl_nlon
      nlev = ccpl_nlev
      call system_clock(end_time)
      total_time = real(end_time-start_time)/10000.
      call mpi_reduce(total_time, max_total_time, 1, MPI_REAL, MPI_MAX, 0, da_ensmean_demo_ccpl_comm, ierr)
#else
      call grid_init_all
#endif
      call decomp_init_da
#ifdef CCPL_DA
      call system_clock(start_time)
      call register_grids_decomps_fields
      call system_clock(end_time)
      total_time = real(end_time-start_time)/10000.
      call mpi_reduce(total_time, max_total_time1, 1, MPI_REAL, MPI_MAX, 0, da_ensmean_demo_ccpl_comm, ierr)
      if (masterpe) then
          write(*,*) "[CCPL <da_ensmean_demo>] ccpl time all: ", max_total_time+max_total_time1
          write(*,*) "[CCPL <da_ensmean_demo>] ccpl time grids: ", max_total_time
          write(*,*) "[CCPL <da_ensmean_demo>] ccpl time fields: ", max_total_time1
      end if
#endif
      call variable_init_da

      !call register_component_coupling_configuration

  end subroutine da_algorithm_demo_init

#ifdef ATM_DA
  subroutine da_algorithm_demo_step_on

      implicit none

      integer      :: i,j
      logical      :: interface_status
      real(kind=4) :: sum, ave, sum1, ave1
      real(kind=4), allocatable :: ran_num_T(:)
      real(kind=4), allocatable :: ran_num_Q(:)
      real(kind=4), allocatable :: ran_num_U(:)
      real(kind=4), allocatable :: ran_num_V(:)
  
      allocate(ran_num_T(NANALS))
      allocate(ran_num_Q(NANALS))
      allocate(ran_num_U(NANALS))
      allocate(ran_num_V(NANALS))
     

      
      do i=1, NANALS 
        ran_num_T(i) = normal(0., 2.)
      end do
      do i=1, NANALS
        ran_num_Q(i) = normal(0., 0.5)
      end do
      do i=1, NANALS
        ran_num_U(i) = normal(0., 2.)
      end do
      do i=1, NANALS
        ran_num_V(i) = normal(0., 2.)
      end do
      
      do i=1, NANALS
      do j =1, local_grid_size(mytask_id+1) 
        rTAUX_m(j,i)  = ran_num_U(i)+rTAUX_l(j)
        rTAUY_m(j,i)  = ran_num_U(i)+rTAUY_l(j)
        rPS_m(j,i) = ran_num_Q(i)+rPS_l(j)
        rT_m(j,:,i) = ran_num_T(i)+rT_l(j,:)
      end do
      end do
     
      sum = 0
      sum1 = 0
      do j=1, NANALS 
        sum =  sum + rTAUX_m(1,j)
        sum1 =  sum1 + ran_num_T(j)
      end do
      ave = sum/NANALS
      ave1 = sum1/NANALS
      IF (mype .EQ. 0) THEN
      write(*,*) "[CCPL <da_ensmean_demo>] date: ", ccpl_date
      write(*,*) "[CCPL <da_ensmean_demo>] DAFCC T: ", ccpl_t(1:10,1)
      write(*,*) "[CCPL <da_ensmean_demo>] READ  T: ", rt_l(1:10,1)
      END IF
      
      call check_diff_1D(1, "TAUX", rTAUX_l, ccpl_TAUX)
      call check_diff_1D(1, "TAUY", rTAUY_l, ccpl_TAUY)
      call check_diff_2D(1, "t", rt_l, ccpl_t)
      
      if(allocated(ran_num_T)) deallocate(ran_num_T)
      if(allocated(ran_num_Q)) deallocate(ran_num_Q)
      if(allocated(ran_num_U)) deallocate(ran_num_U)
      if(allocated(ran_num_V)) deallocate(ran_num_V)

  end subroutine da_algorithm_demo_step_on
#endif

#ifdef OCN_DA
  subroutine da_algorithm_demo_step_on

      implicit none

      integer      :: i,j
      logical      :: interface_status
      real(kind=4) :: sum, ave, sum1, ave1
      real(kind=4), allocatable :: ran_num_sst(:)
      real(kind=4), allocatable :: ran_num_ssh(:)
      real(kind=4), allocatable :: ran_num_salt(:)

      
      allocate(ran_num_SST(NANALS))
      allocate(ran_num_ssh(NANALS))
      allocate(ran_num_salt(NANALS))
     

      
      do i=1, NANALS 
        ran_num_sst(i) = normal(0., 2.)
      end do
      do i=1, NANALS
        ran_num_ssh(i) = normal(0., 0.5)
      end do
      do i=1, NANALS
        ran_num_salt(i) = normal(0., 2.)
      end do

      do i=1, NANALS
      do j =1, local_grid_size(mytask_id+1) 
        rsst_m(j,i)  = ran_num_sst(i)+rsst_l(j)
        rssh_m(j,i)  = ran_num_ssh(i)+rssh_l(j)
        rsalt_m(j,:,i) = ran_num_salt(i)+rsalt_l(j,:)
      end do
      end do
     
      sum = 0
      sum1 = 0
      do j=1, NANALS 
        sum =  sum + rsst_m(1,j)
        sum1 =  sum1 + ran_num_sst(j)
      end do
      ave = sum/NANALS
      ave1 = sum1/NANALS
      IF (mype .EQ. 0) THEN
      write(*,*) "[CCPL <da_ensmean_demo>] date: ", ccpl_date
      write(*,*) "[CCPL <da_ensmean_demo>] DAFCC SST: ", ccpl_sst(1:10)
      write(*,*) "[CCPL <da_ensmean_demo>] READ  SST: ", rsst_l(1:10)
      END IF
      
      call check_diff_1D(1, "SST", rsst_l, ccpl_sst)
      !call check_diff_1D(1, "LAT", ccpl_grid_lat, ccpl_lat)
      !call check_diff_1D(1, "SSH", rssh_l, ccpl_ssh)
      !call check_diff_1D(1, "SALT", rsalt_l, ccpl_salt)
      
      deallocate(ran_num_sst)
      deallocate(ran_num_ssh)
      deallocate(ran_num_salt)

  end subroutine da_algorithm_demo_step_on
#endif

  subroutine finalize_da_algorithm_demo

      implicit none
      integer  :: ierr
        
      CALL MPI_BARRIER(da_ensmean_demo_ccpl_comm,ierr)
      !call deallocate_fields_da
      if(allocated(xlat)) deallocate(xlat)
      if(allocated(xlon)) deallocate(xlon)
      if(allocated(xlev)) deallocate(xlev)
      if(allocated(local_grid_index)) deallocate(local_grid_index)
      
      !call mpi_finalize(ier)

  end subroutine finalize_da_algorithm_demo
end module model_setting_mod_1
