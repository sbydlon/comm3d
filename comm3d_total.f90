program comm3d_total

      use mpi3d_interface_test
      use mpi
      use mpi3d_interface
      use mpi3dcomm, only : cartesian,decompose3d,allocate_array_body,exchange_all_neighbors, &
       test_exchange  
      use mpi3dbasic, only : start_mpi,rank,size,is_master,finish_mpi,time_elapsed, &
       message,new_communicator,ps
      use mpi3dio, only : init_io,file_distributed,open_file_distributed,write_file_distributed,close_file_distributed
 

  implicit none
  
  character(256) :: str,method 
  real :: walltime 
  logical :: interface_comm_test_result,in_block1,in_block2,on_interface 
  integer :: n,cart_size(3),cart_rank,comm_cart,ierr,comm,comm1,comm2,commw
  integer :: size_x_in,size_y_in,size_z_in 
  integer :: coord(3),normal(3),comm_interface,i_rank
  logical,parameter :: reorder = .false., periodic(3) = &
  (/.false.,.false.,.false./)
  logical :: periodic_x,periodic_y,periodic_z,interface_exchange_neighbor_test_result
  type(interface3d) :: I
  integer,parameter :: nx1=61,ny1=61,nz1=61,nx2=61,ny2=61,nz2=61,nb=3,nF=6 ! could read these in from command line or input file
  type(cartesian) :: C1,C2
  type(file_distributed) :: io1,io2  
  real,dimension(:,:,:),allocatable :: F1,F2


  ! start MPI
  call start_mpi
  commw = MPI_COMM_WORLD
  i_rank = -1
  
  in_block1 = (rank< (size+1)/2)
  in_block2 = (rank>= size   /2)
  call new_communicator(in_block1,comm1)
  call new_communicator(in_block2,comm2)
  if(in_block1) comm = comm1
  if(in_block2) comm = comm2

  periodic_x = .false.; periodic_y = .false.; periodic_z = .false.
  method = '3D' ! 2D decomposition
  size_x_in = 0; size_y_in = 0; size_z_in = 0 ! number of processes in each direction (only used if method='manual')
 
  if (in_block1) call decompose3d(C1,nx1,ny1,nz1,nb,nF,periodic_x,periodic_y,periodic_z,comm1,method,size_x_in,size_y_in,size_z_in)
  if (in_block2) call decompose3d(C2,nx2,ny2,nz2,nb,nF,periodic_x,periodic_y,periodic_z,comm2,method,size_x_in,size_y_in,size_z_in)

  ! initialize parallel I/O

  if (in_block1) call init_io(C1)
  if (in_block2) call init_io(C2) 


  !! initialize fields
  if (in_block1) call &
    allocate_array_body(F1,C1,ghost_nodes=.true.,Fval=real(rank))
  if (in_block2) call &
    allocate_array_body(F2,C2,ghost_nodes=.true.,Fval=real(rank))

  ! open files and set view for parallel I/O

  if (in_block1) call open_file_distributed(io1,'F1.dat','write',C1%comm,C1%array_s,ps)
  if (in_block2) call open_file_distributed(io2,'F2.dat','write',C2%comm,C2%array_s,ps)


  ! populate interior field values

  !if (in_comm1) ...
  !if (in_comm2) ...

  ! populate ghost nodes (MPI communication)

  !if (in_comm1) call exchange_all_neighbors(C1,F1)
  !if (in_comm2) call exchange_all_neighbors(C2,F2)

  ! or do all of above in testing subroutine


  if (in_block1) call test_exchange(C1,F1)
  if (in_block2) call test_exchange(C2,F2)


  call MPI_Comm_size(comm,n,ierr)
  cart_size = (/0,0,0/)
  call MPI_Dims_create(n,3,cart_size,ierr)
  
  call MPI_Cart_create(comm,3,cart_size,periodic,reorder,comm_cart,ierr)
  call MPI_Comm_rank(comm_cart,cart_rank,ierr)
  call MPI_Cart_coords(comm_cart,cart_rank,3,coord,ierr)
  
  call MPI_Barrier(MPI_COMM_WORLD,ierr)

  normal = (/ 0,0,1 /)

  if(in_block2) normal = -normal

  if(in_block1) call new_interface(coord,cart_size,normal,commw,C1,I)
  if(in_block2) call new_interface(coord,cart_size,normal,commw,C2,I)

  
  ! Test parallel methods 
  call interface_communicator3d_test(coord,cart_size,normal,&
       commw,interface_comm_test_result,comm_interface,on_interface)

  if(in_block1) call &
   exchange_interface_neighbors3d_test(C1,I,interface_exchange_neighbor_test_result)
  if(in_block2) call &
   exchange_interface_neighbors3d_test(C2,I,interface_exchange_neighbor_test_result) 

  ! output arrays (subroutine will convert F to single precision before writing)

  if (in_block1) call write_file_distributed(io1,F1(C1%mx:C1%px,C1%my:C1%py,C1%mz:C1%pz))
  if (in_block2) call write_file_distributed(io2,F2(C2%mx:C2%px,C2%my:C2%py,C2%mz:C2%pz))

  ! close parallel files

  if (in_block1) call close_file_distributed(io1)
  if (in_block2) call close_file_distributed(io2)

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
   

  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  
  ! Test serial methods
  !if(on_interface)call MPI_Comm_rank(comm_interface,i_rank,ierr)
  !if(i_rank == 0) then
  !call init_fruit
  !call assert_true(interface_comm_test_result,'interface_communicator3d_test')
  !call assert_true(interface_exchange_neighbor_test_result, &
  !     'exchange_interface_neighbors3d_test')
  !call is_on_interface3d_test_all(cart_size)
  !call get_interface_coord3d_test_all(cart_size)
  !call get_dimensions3d_test_all(cart_size)
  !call fruit_summary
  !endif

  !call MPI_Barrier(MPI_COMM_WORLD,ierr)
  call finish_mpi


end program comm3d_total 
