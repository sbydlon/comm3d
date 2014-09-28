module mpi3d_interface

      implicit none

      public :: is_on_interface2d, get_interface_coord2d, &
          get_dimensions2d, interface_communicator2d, interface2d

      ! n             : number of processes on the interface
      ! coord         : global, cartesian coordinates of process
      ! icoord        : interface coordinate
      ! normal        : outward pointing unit normal, with respect to this
      ! side of the interface
      ! comm          : interface communicator (all processes on the interface)
      ! rank          : rank of process
      ! rank_neighbor : rank of process on opposite side of interface
      ! on_interface  : true, if this process is on the interface
      type interface2d
        integer :: n,coord(2),icoord,normal(2),comm,rank,rank_neighbor
        logical :: on_interface
      end type interface2d

      type interface3d
        integer :: n(2),coord(3),icoord(2),normal(3),comm,rank,rank_neighbor
        logical :: on_interface
        integer :: face_xy,face_xz,face_yz
      end type interface3d

      interface new_interface
          module procedure new_interface2d,new_interface3d
      end interface
     
      contains
      
      ! Construct a new interface between two blocks 
      ! processes on the interface are determined by the normal.
      ! The normal must be outward pointing with respect to a block side
      ! This means that n_1 = -n_2, where n_1 is the normal on one side
      ! of the interface, and n_2 is the normal on the opposite side.
      ! All data needed for future communication are contained in I
      subroutine new_interface2d(coord,cart_size,normal,comm,I)

          use mpi

          implicit none

          integer,intent(in) :: coord(2),cart_size(2),normal(2),comm
          type(interface2d),intent(out) :: I

          integer :: ierr

          call interface_communicator2d(coord,cart_size,normal,comm,&
                                        I%comm,I%on_interface)

          if(.not. I%on_interface) return

          call MPI_Comm_rank(I%comm,I%rank,ierr)

          call find_neighbor2d(coord,cart_size,normal,I%comm,&
                               I%rank_neighbor)

          
          I%n = get_dimensions2d(cart_size,normal)
          I%coord = coord
          I%icoord = get_interface_coord2d(coord,normal)
          I%normal = normal

      end subroutine new_interface2d

      subroutine new_interface3d(coord,cart_size,normal,comm,C,I)

          use mpi
          use mpi3dcomm, only : cartesian

          implicit none

          integer,intent(in) :: coord(3),cart_size(3),normal(3),comm
          type(cartesian),intent(in) :: C
          type(interface3d),intent(out) :: I

          integer :: ierr

          call interface_communicator3d(coord,cart_size,normal,comm,&
                                        I%comm,I%on_interface)

          if(.not. I%on_interface) return

          call MPI_Comm_rank(I%comm,I%rank,ierr)

          call find_neighbor3d(coord,cart_size,normal,I%comm,&
                               I%rank_neighbor)

          
          I%n = get_dimensions3d(cart_size,normal)
          I%coord = coord
          I%icoord = get_interface_coord3d(coord,normal)
          I%normal = normal
          
          call define_faces(C,I)

      end subroutine new_interface3d


      ! Check that a process is on the interface.
      function is_on_interface2d(coord,cart_size,normal) & 
          result(is_interface_node)
          
         use mpi3dcomm, only : cartesian
         use mpi, only : MPI_PROC_NULL

         implicit none

         integer :: coord(2),cart_size(2),normal(2)
         logical :: is_interface_node

         ! Check all four interface sides defined by the outward unit
         ! normal with respect to each side
         if (normal(1) == 1 .and. normal(2) == 0 ) then
             is_interface_node = (coord(1) == cart_size(1) - 1)
         elseif (normal(1) == 0 .and. normal(2) == 1 ) then
             is_interface_node = (coord(2) == cart_size(2) - 1)
         elseif (normal(1) == -1 .and. normal(2) == 0 ) then
             is_interface_node = (coord(1) == 0)
         elseif (normal(1) ==  0 .and. normal(2) == -1 ) then
             is_interface_node = (coord(2) == 0)
         endif

       end function is_on_interface2d
      
      ! Check that a process is on the interface.
      function is_on_interface3d(coord,cart_size,normal) & 
          result(is_interface_node)
          
         use mpi3dcomm, only : cartesian
         use mpi, only : MPI_PROC_NULL

         implicit none

         integer :: coord(3),cart_size(3),normal(3)
         logical :: is_interface_node

         ! Check all four interface sides defined by the outward unit
         ! normal with respect to each side
         if (normal(1) == 1 .and. normal(2) == 0 .and. normal(3) == 0) then
             is_interface_node = (coord(1) == cart_size(1) - 1)
         elseif (normal(1) == 0 .and. normal(2) == 1 .and. normal(3) == 0) then
             is_interface_node = (coord(2) == cart_size(2) - 1)
         elseif (normal(1) == -1 .and. normal(2) == 0 .and. normal(3) == 0 ) then
             is_interface_node = (coord(1) == 0)
         elseif (normal(1) ==  0 .and. normal(2) == -1 .and. normal(3) == 0 ) then
             is_interface_node = (coord(2) == 0)
         elseif (normal(1) ==  0 .and. normal(2) == 0 .and. normal(3) ==  1 ) then
             is_interface_node = (coord(3) == cart_size(3) - 1)
         elseif (normal(1) ==  0 .and. normal(2) == 0 .and. normal(3) == -1 ) then
             is_interface_node = (coord(3) == 0)
         endif

       end function is_on_interface3d

       ! Determine the number of processes on the interface
       function get_dimensions2d(cart_size,normal) &
           result(interface_size)

           implicit none 

           integer,intent(in) :: cart_size(2),normal(2)
           integer :: interface_size
           
           if (abs(normal(1)) == 1 .and. normal(2) == 0 ) then
               interface_size = cart_size(2)
           elseif (normal(1) == 0 .and. abs(normal(2)) == 1 ) then
               interface_size = cart_size(1)
           endif

       end function get_dimensions2d

       function get_dimensions3d(cart_size,normal) &
           result(interface_size)

           implicit none 

           integer,intent(in) :: cart_size(3),normal(3)
           integer :: interface_size(2)

           ! Plane convention
           ! x : (y,z)
           ! y : (x,z)
           ! z : (x,y)
           
           if (abs(normal(1)) == 1 .and. normal(2) == 0  .and. normal(3) == 0) then
               interface_size = (/ cart_size(2), cart_size(3) /)
           elseif (normal(1) == 0 .and. abs(normal(2)) == 1 .and. normal(3) == 0) then
               interface_size = (/ cart_size(1), cart_size(3) /) 
           elseif (normal(1) == 0 .and. normal(2) == 0 .and. abs(normal(3)) ==  1 ) then
               interface_size = (/ cart_size(1), cart_size(2) /)
           endif

       end function get_dimensions3d
       ! Convert between the global block coordinate system
       ! to a local, interface coordinate system
       function get_interface_coord2d(coord,normal) &
           result(interface_coord)

           implicit none

           integer,intent(in) :: coord(2), normal(2)
           integer :: interface_coord, tangent(2)

           ! Tangent vector with positive components
           tangent = (/ abs(normal(2)), abs(normal(1)) /) 

           interface_coord = tangent(1) * coord(1) + tangent(2) * coord(2)
       
       end function get_interface_coord2d
       

       function get_interface_coord3d(coord,normal) &
           result(interface_coord)

           implicit none

           integer,intent(in) :: coord(3), normal(3)
           integer :: interface_coord(2)

           if (abs(normal(1)) == 1 .and. normal(2) == 0  .and. normal(3) == 0) then
               interface_coord = (/ coord(2), coord(3) /)
           elseif (normal(1) == 0 .and. abs(normal(2)) == 1 .and. normal(3) == 0) then
               interface_coord = (/ coord(1), coord(3) /) 
           elseif (normal(1) == 0 .and. normal(2) == 0 .and. abs(normal(3)) ==  1 ) then
               interface_coord = (/ coord(1), coord(2) /)
           endif
       
       end function get_interface_coord3d

       ! Construct an interface communicator, 
       ! only processes on the interface with opposing normals
       ! and in 'comm' are included in 'comm_interface'
       !
       ! Example: Block1:  (0, 1) (1, 1)
       !                   (0, 0) (1, 0)
       !          Block2:  (0, 1) (1, 1)
       !                   (0, 0) (1, 0)
       ! The interface between blocks are defined with normal (0,-1), for Block1
       ! and normal (0,1), for Block2. 
       subroutine interface_communicator2d(coord,cart_size,normal,comm, &
                    comm_interface,on_interface)

           use mpi
           use mpi3dbasic, only : new_communicator

           implicit none

           integer,intent(in) :: coord(2),cart_size(2),normal(2),comm
           integer,intent(out) :: comm_interface
           logical,intent(out) :: on_interface

           on_interface = .false.

           ! Determine processes belonging to the interface
           if (is_on_interface2d(coord,cart_size,normal)) on_interface = .true.
           
           call new_communicator(on_interface,comm_interface,comm)

       end subroutine interface_communicator2d
       
       subroutine interface_communicator3d(coord,cart_size,normal,comm, &
                    comm_interface,on_interface)

           use mpi
           use mpi3dbasic, only : new_communicator

           implicit none

           integer,intent(in) :: coord(3),cart_size(3),normal(3),comm
           integer,intent(out) :: comm_interface
           logical,intent(out) :: on_interface

           on_interface = .false.

           ! Determine processes belonging to the interface
           if (is_on_interface3d(coord,cart_size,normal)) on_interface = .true.
           
           call new_communicator(on_interface,comm_interface,comm)

       end subroutine interface_communicator3d
       
       ! Find neighbors on opposite side of interface.
       ! This subroutine must only be called by processes in interface_comm
       subroutine find_neighbor2d(coord,cart_size,normal,comm,rank_neighbor)

           use mpi

           implicit none

           integer,intent(in) :: coord(2), cart_size(2), normal(2),comm 
           integer :: n,rank,interface_coord,ierr,i
           integer,dimension(:),allocatable :: ranks,x
           integer,intent(out) :: rank_neighbor 

           ! Ensure that only processes on the interface participate
           if (.not. is_on_interface2d(coord,cart_size,normal)) return
           
           n = get_dimensions2d(cart_size,normal)
           
           call MPI_Comm_rank(comm,rank,ierr)
           interface_coord = get_interface_coord2d(coord,normal)

           ! There are 2*n processes on the interface (n on each side)
           allocate(ranks(2*n),x(2*n))
           
           call MPI_Allgather(rank,1,MPI_INTEGER,ranks,1,MPI_INTEGER,comm,ierr)
           call MPI_Allgather(interface_coord,1,MPI_INTEGER,x,1,MPI_INTEGER,comm)

           ! Search for process with the same interface coordinates
           ! as self
           do i=1,2*n
              if (x(i) == interface_coord .and. ranks(i) /= rank) exit
           end do
           ! @todo should check that x(i) == interface_coord, if not true, 
           ! then the block sides are not adjacent 
           if (i > 2*n) i = 2*n

           rank_neighbor = ranks(i)
           
           deallocate(ranks,x)

       end subroutine find_neighbor2d

       subroutine find_neighbor3d(coord,cart_size,normal,comm,rank_neighbor)

           use mpi

           implicit none

           integer,intent(in) :: coord(3), cart_size(3), normal(3),comm 
           integer :: n(2),rank,interface_coord(2),ierr,i,nx,ny
           integer,dimension(:),allocatable :: ranks,x,y
           integer,intent(out) :: rank_neighbor 

           ! Ensure that only processes on the interface participate
           if (.not. is_on_interface3d(coord,cart_size,normal)) return
           
           n = get_dimensions3d(cart_size,normal)
           nx = n(1)
           ny = n(2)
           
           call MPI_Comm_rank(comm,rank,ierr)
           interface_coord = get_interface_coord3d(coord,normal)

           ! There are 2*n processes on the interface (n on each side)
           allocate(ranks(2*nx*ny),x(2*nx*ny),y(2*nx*ny))
           
           call MPI_Allgather(rank,1,MPI_INTEGER,ranks,1,MPI_INTEGER,comm,ierr)
           call MPI_Allgather(interface_coord(1),1,MPI_INTEGER,x,1,MPI_INTEGER,comm)
           call MPI_Allgather(interface_coord(2),1,MPI_INTEGER,y,1,MPI_INTEGER,comm)

           ! Search for process with the same interface coordinates
           ! as self
           do i=1,2*n(1)*n(2)
              if (x(i) == interface_coord(1) .and. &
                  y(i) == interface_coord(2) .and. ranks(i) /= rank) exit
           end do
           ! @todo should check that x(i) == interface_coord, if not true, 
           ! then the block sides are not adjacent 
           if (i > 2*n(1)*n(2)) i = 2*n(1)*n(2)

           rank_neighbor = ranks(i)
           
           deallocate(ranks,x,y)

       end subroutine find_neighbor3d

       subroutine exchange_interface_neighbors2d(C,F,I)

           use mpi
           use mpi3dcomm, only : cartesian

           implicit none

           type(cartesian),intent(in) :: C
           real,dimension(C%mbx:C%pbx,C%mby:C%pby),intent(inout) :: F
           type(interface2d),intent(in) :: I
           integer,parameter :: tag1 = 1,tag2 = 2,tag3 = 3,tag4 = 4
           integer :: ierr

           if(.not. I%on_interface) return
           
           if (I%normal(1) == 1 .and. I%normal(2) == 0) then
               call MPI_SendRecv( &
                    F(C%px,C%my),1,C%line_y,I%rank_neighbor,tag1, &
                    F(C%px+1,C%my),1,C%line_y,I%rank_neighbor,tag2, &
                    I%comm,MPI_STATUS_IGNORE,ierr)
           elseif (I%normal(1) == -1 .and. I%normal(2) == 0) then
               call MPI_SendRecv( &
                    F(C%mx,C%my),1,C%line_y,I%rank_neighbor,tag2, &
                    F(C%mx-1,C%my),1,C%line_y,I%rank_neighbor,tag1, &
                    I%comm,MPI_STATUS_IGNORE,ierr)          
           elseif (I%normal(1) == 0 .and. I%normal(2) == 1) then
               call MPI_SendRecv( &
                    F(C%mx,C%py),1,C%line_x,I%rank_neighbor,tag3, &
                    F(C%mx,C%py+1),1,C%line_x,I%rank_neighbor,tag4, &
                    I%comm,MPI_STATUS_IGNORE,ierr)
           elseif (I%normal(1) == 0 .and. I%normal(2) == -1) then
               call MPI_SendRecv( &
                    F(C%mx,C%my),1,C%line_x,I%rank_neighbor,tag4, &
                    F(C%mx,C%my-1),1,C%line_x,I%rank_neighbor,tag3, &
                    I%comm,MPI_STATUS_IGNORE,ierr)          
           endif

       end subroutine exchange_interface_neighbors2d
       
       subroutine exchange_interface_neighbors3d(C,F,I)

           use mpi
           use mpi3dcomm, only : cartesian

           implicit none

           type(cartesian),intent(in) :: C
           real,dimension(C%mbx:C%pbx,C%mby:C%pby,C%mbz:C%pbz),intent(inout) :: F
           type(interface3d),intent(in) :: I
           integer,parameter :: tag1 = 1,tag2 = 2,tag3 = 3,tag4 = 4, &
                                tag5 = 5,tag6 = 6
           integer :: ierr,nx,ny,nz

           if(.not. I%on_interface) return
           
           nx = I%normal(1)
           ny = I%normal(2)
           nz = I%normal(3)


           if (nx == 1 .and. ny ==  0 .and. nz == 0) then
               call MPI_SendRecv( &
                    F(C%px,C%my,C%mz),1,I%face_yz,I%rank_neighbor,tag1, &
                    F(C%px+1,C%my,C%mz),1,I%face_yz,I%rank_neighbor,tag2, &
                    I%comm,MPI_STATUS_IGNORE,ierr)          
           endif
           
           if (nx == -1 .and. ny ==  0 .and. nz == 0) then
               call MPI_SendRecv( &
                    F(C%mx,C%my,C%mz),1,I%face_yz,I%rank_neighbor,tag2, &
                    F(C%mx-1,C%my,C%mz),1,I%face_yz,I%rank_neighbor,tag1, &
                    I%comm,MPI_STATUS_IGNORE,ierr)          
           endif
           
           if (nx == 0 .and. ny ==  1 .and. nz == 0) then
               call MPI_SendRecv( &
                    F(C%mx,C%py,C%mz),1,I%face_xz,I%rank_neighbor,tag3, &
                    F(C%mx,C%py+1,C%mz),1,I%face_xz,I%rank_neighbor,tag4, &
                    I%comm,MPI_STATUS_IGNORE,ierr)          
           endif

           if (nx == 0 .and. ny == -1 .and. nz == 0) then
               call MPI_SendRecv( &
                    F(C%mx,C%my,C%mz),1,I%face_xz,I%rank_neighbor,tag4, &
                    F(C%mx,C%my-1,C%mz),1,I%face_xz,I%rank_neighbor,tag3, &
                    I%comm,MPI_STATUS_IGNORE,ierr)          
           endif

           if (nx == 0 .and. ny ==  0 .and. nz == 1) then
               call MPI_SendRecv( &
                    F(C%mx,C%my,C%pz),1,I%face_xy,I%rank_neighbor,tag5, &
                    F(C%mx,C%my,C%pz+1),1,I%face_xy,I%rank_neighbor,tag6, &
                    I%comm,MPI_STATUS_IGNORE,ierr)          
           endif
           
           if (nx == 0 .and. ny ==  0 .and. nz == -1) then
               call MPI_SendRecv( &
                    F(C%mx,C%my,C%mz),1,I%face_xy,I%rank_neighbor,tag6, &
                    F(C%mx,C%my,C%mz-1),1,I%face_xy,I%rank_neighbor,tag5, &
                    I%comm,MPI_STATUS_IGNORE,ierr)          
           endif
           
       end subroutine exchange_interface_neighbors3d

       ! Consider moving this to Cartesian
       subroutine define_faces(C,I)

           use mpi3dbasic, only : MPI_REAL_PW
           use mpi3dcomm, only : cartesian
           use mpi

           implicit none

           type(cartesian),intent(in) :: C
           type(interface3d),intent(inout) :: I

           integer :: count_xy,length_xy,stride_xy, &
                      count_yz,length_yz,stride_yz, &
                      count_xz,length_xz,stride_xz

           integer :: ierr,pw

           ! Define the face_xy, given by the plane with normal: (0 0 a)
           count_xy = C%lny
           stride_xy =  (C%lnx + 2*C%nb)
           length_xy = C%lnx
           call MPI_Type_vector(&
           count_xy,length_xy,stride_xy,MPI_REAL_PW,I%face_xy,ierr)
           call MPI_Type_commit(I%face_xy,ierr)

           ! Define the face_yz given by the plane with normal: (a 0 0)
           call MPI_Type_size(MPI_REAL_PW,pw,ierr)
           count_yz = C%lnz
           stride_yz = (C%lnx + 2*C%nb)*(C%lny + 2*C%nb)*pw
           length_yz = 1
           
           ! Note the use of hvector, which specifies the stride in bytes
           call MPI_Type_hvector(&
           count_yz,length_yz,stride_yz,C%line_y,I%face_yz,ierr)
           call MPI_Type_commit(I%face_yz,ierr)

           ! Define the face_xz given by the plane with normal: (0 a 0)
           count_xz = C%lnz
           stride_xz = (C%lnx + 2*C%nb)*(C%lny + 2*C%nb)*pw
           length_xz = 1
           
           call MPI_Type_hvector(&
           count_xz,length_xz,stride_xz,C%line_x,I%face_xz,ierr)
           call MPI_Type_commit(I%face_xz,ierr)
       end subroutine define_faces

end module mpi3d_interface
