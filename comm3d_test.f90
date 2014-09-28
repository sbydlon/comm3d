program comm3d

  use mpi3dbasic, only : start_mpi,rank,size,is_master,finish_mpi,time_elapsed, &
       message,new_communicator,ps
  use mpi3dcomm, only : cartesian,decompose3d,allocate_array_body,exchange_all_neighbors, &
       test_exchange
  use mpi3dio, only : init_io,file_distributed,open_file_distributed,write_file_distributed,close_file_distributed

  implicit none

  character(256) :: str,method
  real :: walltime
  logical :: in_comm1,in_comm2,periodic_x,periodic_y,periodic_z
  integer :: comm1,comm2,size_x_in,size_y_in,size_z_in,i,j,k
  type(cartesian) :: C1,C2
  type(file_distributed) :: io1,io2
  integer,parameter :: nx1=21,ny1=11,nz1=21,nx2=21,ny2=11,nz2=21,nb=3,nF=6 ! could read these in from command line or input file

  real,dimension(:,:,:),allocatable :: F1,F2

  ! start MPI
  
  call start_mpi

  ! split MPI_COMM_WORLD into two communicators
  ! (for odd size, one process will be in both communicators)
  ! method uses integer division by truncation

  in_comm1 = (rank< (size+1)/2)
  in_comm2 = (rank>= size   /2)

  write(str,'(a,l1,a,l1)') &
       'in_comm1 = ',in_comm1,' in_comm2 = ',in_comm2
  call message(str)
  call new_communicator(in_comm1,comm1)
  call new_communicator(in_comm2,comm2)

  ! decompose each block, creating Cartesian communicators

  periodic_x = .false.; periodic_y = .false.; periodic_z = .false.
  method = '3D' ! 2D decomposition
  size_x_in = 0; size_y_in = 0; size_z_in = 0 ! number of processes in each direction (only used if method='manual')

  if (in_comm1) call decompose3d(C1,nx1,ny1,nz1,nb,nF,periodic_x,periodic_y,periodic_z,comm1,method,size_x_in,size_y_in,size_z_in)

  if (in_comm2) call decompose3d(C2,nx2,ny2,nz2,nb,nF,periodic_x,periodic_y,periodic_z,comm2,method,size_x_in,size_y_in,size_z_in)
    
  ! initialize parallel I/O

  if (in_comm1) call init_io(C1)
  if (in_comm2) call init_io(C2)

  ! initialize fields

  if (in_comm1) call allocate_array_body(F1,C1,ghost_nodes=.true.)
  if (in_comm2) call allocate_array_body(F2,C2,ghost_nodes=.true.)
              
  ! open files and set view for parallel I/O

  if (in_comm1) call open_file_distributed(io1,'F1.dat','write',C1%comm,C1%array_s,ps)
  if (in_comm2) call open_file_distributed(io2,'F2.dat','write',C2%comm,C2%array_s,ps)
  ! populate interior field values

  !if (in_comm1) ...
  !if (in_comm2) ...

  ! populate ghost nodes (MPI communication)

  !if (in_comm1) call exchange_all_neighbors(C1,F1)
  !if (in_comm2) call exchange_all_neighbors(C2,F2)

  ! or do all of above in testing subroutine


  if (in_comm1) call test_exchange(C1,F1)
  if (in_comm2) call test_exchange(C2,F2)

  ! output arrays (subroutine will convert F to single precision before writing)

  if (in_comm1) call write_file_distributed(io1,F1(C1%mx:C1%px,C1%my:C1%py,C1%mz:C1%pz))
  if (in_comm2) call write_file_distributed(io2,F2(C2%mx:C2%px,C2%my:C2%py,C2%mz:C2%pz))

  ! close parallel files

  if (in_comm1) call close_file_distributed(io1)
  if (in_comm2) call close_file_distributed(io2)

  ! to examine files in MATLAB:
  ! f = fopen('F1.dat','rb')
  ! F1 = fread(f,[nx ny],'real*4')
  ! fclose(f)

  ! timing information

  if (is_master()) then
     walltime = time_elapsed()
     write(str,'(a,i0,a,f0.8,a)') &
          'number of processes = ',size,'  wall clock time = ',walltime,' s'
     call message(str)
  end if


  ! finish MPI

  call finish_mpi

end program comm3d
