module mpi3d_interface_test

  implicit none

  private

  public :: is_on_interface2d_test_all, &
            is_on_interface3d_test_all, &
            get_interface_coord2d_test_all, &
            get_interface_coord3d_test_all, &
            get_dimensions2d_test_all, &
            get_dimensions3d_test_all, &
            interface_communicator2d_test, &
            interface_communicator3d_test, &
            exchange_interface_neighbors2d_test, &
            exchange_interface_neighbors3d_test
  
contains
  subroutine is_on_interface2d_test_all(cart_size)

      implicit none

      integer,intent(in) :: cart_size(2)

      call is_on_interface2d_test((/  1,  0 /),cart_size) 
      call is_on_interface2d_test((/  0,  1 /),cart_size) 
      call is_on_interface2d_test((/ -1,  0 /),cart_size) 
      call is_on_interface2d_test((/  0, -1 /),cart_size) 
  end subroutine
  
  subroutine is_on_interface3d_test_all(cart_size)

      implicit none

      integer,intent(in) :: cart_size(3)

      call is_on_interface3d_test((/  1,  0, 0 /),cart_size) 
      call is_on_interface3d_test((/ -1,  0, 0 /),cart_size) 
      call is_on_interface3d_test((/  0,  1, 0 /),cart_size) 
      call is_on_interface3d_test((/  0, -1, 0 /),cart_size) 
      call is_on_interface3d_test((/  0,  0, 1 /),cart_size) 
      call is_on_interface3d_test((/  0,  0, -1 /),cart_size) 
  end subroutine


  ! Check that the interface nodes are found for a given interface
  ! defined by outward pointing unit normal 'normal'
  subroutine is_on_interface2d_test(normal,cart_size)

      use mpi3d_interface, only : is_on_interface2d
      use fruit

      implicit none

      integer,intent(in) :: normal(2),cart_size(2)
      integer,dimension(:),allocatable :: x,y
      integer :: coord(2),k
      logical :: test_passed
      test_passed = .true.

      call nodes_on_interface2d(normal,cart_size,x,y)

      do k=1,size(x)   
        coord = (/x(k), y(k)/)      
        if(.NOT. is_on_interface2d(coord,cart_size,normal)) & 
            test_passed = .false.
      end do

      deallocate(x,y)
      call assert_true(test_passed,'is_on_interface2d_test') 

  end subroutine is_on_interface2d_test
  
  subroutine is_on_interface3d_test(normal,cart_size)

      use mpi3d_interface, only : is_on_interface3d
      use fruit

      implicit none

      integer,intent(in) :: normal(3),cart_size(3)
      integer,dimension(:),allocatable :: x,y,z
      integer :: coord(3),k
      logical :: test_passed
      test_passed = .true.

      call nodes_on_interface3d(normal,cart_size,x,y,z)
      
      do k=1,size(x)   
        coord = (/ x(k), y(k), z(k) /)      
        if(.NOT. is_on_interface3d(coord,cart_size,normal)) & 
            test_passed = .false.
      end do

      deallocate(x,y,z)
      call assert_true(test_passed,'is_on_interface3d_test') 

  end subroutine is_on_interface3d_test

  subroutine nodes_on_interface2d(normal,cart_size,x,y)

      implicit none
      integer,intent(in) :: normal(2),cart_size(2)
      integer,intent(out),dimension(:),allocatable :: x,y

      if(normal(1) == 1 .and. normal(2) == 0) then 
         allocate(x(cart_size(2)),y(cart_size(2)))
         y = (/0,cart_size(2)-1/)
         x = 0*y + cart_size(1) - 1
      elseif(normal(1) == 0 .and. normal(2) == 1) then
         allocate(x(cart_size(1)),y(cart_size(1)))
         x = (/0,cart_size(1)-1/)
         y = 0*y + cart_size(2) - 1
      elseif(normal(1) == -1 .and. normal(2) == 0) then
         allocate(x(cart_size(2)),y(cart_size(2)))
         y = (/0,cart_size(2)-1/)
         x = 0*y
      elseif(normal(1) == 0 .and. normal(2) == -1) then
         allocate(x(cart_size(1)),y(cart_size(1)))
         x = (/0,cart_size(1)-1/)
         y = 0*y 
      endif
  end subroutine
  
  subroutine nodes_on_interface3d(normal,cart_size,x,y,z)

      implicit none
      integer,intent(in) :: normal(3),cart_size(3)
      integer,intent(out),dimension(:),allocatable :: x,y,z

      integer :: nx,ny,nz,N
         
      nx = cart_size(1)
      ny = cart_size(2)
      nz = cart_size(3)

      if(normal(1) == 1 .and. normal(2) == 0 .and. normal(3) == 0) then 
         N = nz*ny
         allocate(x(N),y(N),z(N))
         call grid(y,z,ny,nz)
         x = 0*y + nx - 1
      elseif(normal(1) == -1 .and. normal(2) == 0 .and. normal(3) == 0) then
         N = nz*ny
         allocate(x(N),y(N),z(N))
         call grid(y,z,ny,nz)
         x = 0*y
      elseif(normal(1) == 0 .and. normal(2) == 1 .and. normal(3) == 0) then
         N = nx*nz
         allocate(x(N),y(N),z(N))
         call grid(x,z,nx,nz)
         y = 0*x + ny - 1
      elseif(normal(1) == 0 .and. normal(2) == -1 .and. normal(3) == 0) then
         N = nx*nz
         allocate(x(N),y(N),z(N))
         call grid(x,z,nx,nz)
         y = 0*x
      elseif(normal(1) == 0 .and. normal(2) == 0 .and. normal(3) == 1) then
         N = nx*ny
         allocate(x(N),y(N),z(N))
         call grid(x,y,nx,ny)
         z = 0*x + nz - 1
      elseif(normal(1) == 0 .and. normal(2) == 0 .and. normal(3) == -1) then
         N = nx*ny
         allocate(x(N),y(N),z(N))
         call grid(x,y,nx,ny)
         z = 0*x
      endif
  end subroutine

  subroutine grid(x,y,nx,ny)

      implicit none

      integer,dimension(*),intent(inout) :: x,y
      integer,intent(in) :: nx,ny
      
      integer :: i,j,k

      k = 1
      do i=1,nx
         do j=1,ny
           x(k) = i + j*nx
           y(k) = j + i*ny
           k = k + 1
         end do
      end do

  end subroutine grid

  subroutine get_interface_coord2d_test_all(cart_size)

      implicit none

      integer,intent(in) :: cart_size(2)

      call get_interface_coord2d_test((/  1,  0 /),cart_size,2) 
      call get_interface_coord2d_test((/  0,  1 /),cart_size,1) 
      call get_interface_coord2d_test((/ -1,  0 /),cart_size,2) 
      call get_interface_coord2d_test((/  0, -1 /),cart_size,1) 

  end subroutine
  
  subroutine get_interface_coord3d_test_all(cart_size)

      implicit none

      integer,intent(in) :: cart_size(3)

      call get_interface_coord3d_test((/  1,  0,  0 /),cart_size,(/2,3/)) 
      call get_interface_coord3d_test((/ -1,  0,  0 /),cart_size,(/2,3/)) 
      call get_interface_coord3d_test((/  0,  1,  0 /),cart_size,(/1,3/)) 
      call get_interface_coord3d_test((/  0, -1,  0 /),cart_size,(/1,3/)) 
      call get_interface_coord3d_test((/  0,  0,  1 /),cart_size,(/1,2/)) 
      call get_interface_coord3d_test((/  0,  0, -1 /),cart_size,(/1,2/)) 

  end subroutine


  subroutine get_interface_coord2d_test(normal,cart_size,coord_comp) 
      
      use mpi3d_interface, only : is_on_interface2d, get_interface_coord2d
      use fruit

      implicit none

      integer,intent(in) :: normal(2),cart_size(2)
      integer,dimension(:),allocatable :: x,y
      integer :: coord(2),k,coord_comp,interface_coord
      logical :: test_passed
      
      test_passed = .true.
      
      call nodes_on_interface2d(normal,cart_size,x,y)
      
      do k=1,size(x)   
        coord = (/x(k), y(k)/)      
        if (is_on_interface2d(coord,cart_size,normal)) then 
            interface_coord = get_interface_coord2d(coord,normal)
            if (interface_coord /= coord(coord_comp)) test_passed = .false.
        end if
      end do

      deallocate(x,y)
      call assert_true(test_passed,'get_interface_coord2d_test') 

  end subroutine get_interface_coord2d_test
  
  subroutine get_interface_coord3d_test(normal,cart_size,coord_comp) 
      
      use mpi3d_interface, only : is_on_interface3d, get_interface_coord3d
      use fruit

      implicit none

      integer,intent(in) :: normal(3),cart_size(3)
      integer,dimension(:),allocatable :: x,y,z
      integer :: coord(3),k,coord_comp(2),interface_coord(2)
      logical :: test_passed
      
      test_passed = .true.
      
      call nodes_on_interface3d(normal,cart_size,x,y,z)
      
      do k=1,size(x)   
        coord = (/x(k), y(k), z(k) /)      
        if (is_on_interface3d(coord,cart_size,normal)) then 
            interface_coord = get_interface_coord3d(coord,normal)
            if (interface_coord(1) /= coord(coord_comp(1)) .or. &
                interface_coord(2) /= coord(coord_comp(2)) ) test_passed = .false.
        end if
      end do

      deallocate(x,y,z)
      call assert_true(test_passed,'get_interface_coord3d_test') 

  end subroutine get_interface_coord3d_test

  subroutine get_dimensions2d_test_all(cart_size)

      implicit none

      integer,intent(in) :: cart_size(2)
      
      call get_dimensions2d_test(cart_size,(/  1,  0 /),cart_size(2))
      call get_dimensions2d_test(cart_size,(/  0,  1 /),cart_size(1))
      call get_dimensions2d_test(cart_size,(/ -1,  0 /),cart_size(2))
      call get_dimensions2d_test(cart_size,(/  0, -1 /),cart_size(1))


  end subroutine get_dimensions2d_test_all

  subroutine get_dimensions3d_test_all(cart_size)

      implicit none

      integer,intent(in) :: cart_size(3)
      
      integer :: nx,ny,nz

      nx = cart_size(1)
      ny = cart_size(2)
      nz = cart_size(3)
      
      call get_dimensions3d_test(cart_size,(/  1,  0, 0 /),(/ny,nz/))
      call get_dimensions3d_test(cart_size,(/ -1,  0, 0 /),(/ny,nz/))
      call get_dimensions3d_test(cart_size,(/  0,  1, 0 /),(/nx,nz/))
      call get_dimensions3d_test(cart_size,(/  0, -1, 0 /),(/nx,nz/))
      call get_dimensions3d_test(cart_size,(/  0,  0, 1 /),(/nx,ny/))
      call get_dimensions3d_test(cart_size,(/  0,  0, -1 /),(/nx,ny/))


  end subroutine get_dimensions3d_test_all

  subroutine get_dimensions2d_test(cart_size,normal,interface_size)

      use mpi3d_interface, only : get_dimensions2d
      use fruit

      implicit none

      integer,intent(in) :: normal(2),cart_size(2),interface_size
      logical :: test_passed
      
      test_passed = (get_dimensions2d(cart_size,normal) == interface_size)
      
      call assert_true(test_passed,'get_interface_coord2d_test') 

  end subroutine get_dimensions2d_test
  
  subroutine get_dimensions3d_test(cart_size,normal,interface_size)

      use mpi3d_interface, only : get_dimensions3d
      use fruit

      implicit none

      integer,intent(in) :: normal(3),cart_size(3),interface_size(2)
      logical :: test_passed
      integer :: test_dim(2)
      
      test_dim = get_dimensions3d(cart_size,normal)

      test_passed = (test_dim(1) == interface_size(1) .and. &
                     test_dim(2) == interface_size(2) )
      
      call assert_true(test_passed,'get_interface_coord3d_test') 

  end subroutine get_dimensions3d_test

  subroutine interface_communicator2d_test(coord,cart_size,normal,&
             comm,test_result,comm_interface,on_interface)

      use mpi3d_interface, only : interface_communicator2d, &
                                  get_dimensions2d,find_neighbor2d, &
                                  get_interface_coord2d
      use mpi

      implicit none
      
      integer,intent(in) :: coord(2),normal(2),cart_size(2),comm
      logical,intent(out) :: test_result,on_interface
      integer,intent(out) :: comm_interface
      integer :: rank,rank_neighbor,ierr
      integer :: interface_coord,interface_neighbor_coord
      logical :: test_passed

      on_interface = .false.
      test_passed = .false.

      call interface_communicator2d(coord,cart_size,normal,comm,&
                                    comm_interface,on_interface)

      if(.not. on_interface) return

      call find_neighbor2d(coord,cart_size,normal,comm_interface,&
                           rank_neighbor)

      ! Send interface coordinates to neighbors
      interface_coord = get_interface_coord2d(coord,normal)

      call MPI_Comm_rank(comm,rank,ierr)
      call MPI_Comm_rank(comm_interface,rank_neighbor,ierr)

      call MPI_Sendrecv(interface_coord,1,MPI_INTEGER,rank_neighbor,&
                        1,interface_neighbor_coord,1,MPI_INTEGER,rank_neighbor,&
                        1,comm_interface,MPI_STATUS_IGNORE,ierr)
      test_passed = (interface_coord == interface_neighbor_coord)   

      call MPI_Reduce(test_passed,test_result,1,MPI_LOGICAL,MPI_LAND,0,comm_interface,ierr)


  end subroutine interface_communicator2d_test
  
  subroutine interface_communicator3d_test(coord,cart_size,normal,&
             comm,test_result,comm_interface,on_interface)

      use mpi3d_interface, only : interface_communicator3d, &
                                  get_dimensions3d,find_neighbor3d, &
                                  get_interface_coord3d
      use mpi

      implicit none
      
      integer,intent(in) :: coord(3),normal(3),cart_size(3),comm
      logical,intent(out) :: test_result,on_interface
      integer,intent(out) :: comm_interface
      integer :: rank,rank_neighbor,ierr,mpi_coord
      integer :: interface_coord(2),interface_neighbor_coord(2)
      logical :: test_passed

      on_interface = .false.
      test_passed = .false.

      call interface_communicator3d(coord,cart_size,normal,comm,&
                                    comm_interface,on_interface)

      if(.not. on_interface) return

      call find_neighbor3d(coord,cart_size,normal,comm_interface,&
                           rank_neighbor)

      ! Send interface coordinates to neighbors
      interface_coord = get_interface_coord3d(coord,normal)

      call MPI_Type_vector(2,1,1,MPI_INTEGER,mpi_coord,ierr)
      call MPI_Type_commit(mpi_coord,ierr)

      call MPI_Comm_rank(comm,rank,ierr)
      call MPI_Comm_rank(comm_interface,rank_neighbor,ierr)

      call MPI_Sendrecv(interface_coord,1,mpi_coord,rank_neighbor,&
                        1,interface_neighbor_coord,1,mpi_coord,rank_neighbor,&
                        1,comm_interface,MPI_STATUS_IGNORE,ierr)
      test_passed = (interface_coord(1) == interface_neighbor_coord(1) .and. &
                     interface_coord(2) == interface_neighbor_coord(2))   
      
      call MPI_Reduce(test_passed,test_result,1,MPI_LOGICAL,MPI_LAND,0,comm_interface,ierr)

  end subroutine interface_communicator3d_test

  subroutine exchange_interface_neighbors2d_test(C,I,test_result)

      use mpi
      use mpi3d_interface
      use mpi3dcomm, only : cartesian

      implicit none

      type(cartesian),intent(in) :: C
      real,dimension(C%mbx:C%pbx,C%mby:C%pby) :: F
      type(interface2d),intent(in) :: I
      logical,intent(out) :: test_result
      integer :: j,k,ierr
      logical :: test_passed


      test_passed = .true.

      if(.not. I%on_interface) return
      
      ! Fill all data cells excluding ghost cells
      do j=C%mx,C%px
        do k=C%my,C%py
            F(j,k) = real(I%rank)
        end do
      end do

      call exchange_interface_neighbors2d(C,F,I)

      ! Check that ghost cells match neighboring boundary

      if (I%normal(1) == 1 .and. I%normal(2) == 0) then
        do j=C%px+1,C%px+1
          do k=C%my,C%py
            if (abs(F(j,k) - I%rank_neighbor) > 1e-12  ) test_passed = .false.
          end do
        end do
      endif

      if (I%normal(1) == -1 .and. I%normal(2) == 0) then
        do j=C%mx-1,C%mx-1
          do k=C%my,C%py
            if (abs(F(j,k) - I%rank_neighbor) > 1e-12  ) test_passed = .false.
          end do
        end do
      endif
      
      if (I%normal(1) == 0 .and. I%normal(2) == 1) then
        do j=C%mx,C%px
          do k=C%py+1,C%py+1
            if (abs(F(j,k) - I%rank_neighbor) > 1e-12  ) test_passed = .false.
          end do
        end do
      endif
      
      if (I%normal(1) == 0 .and. I%normal(2) == -1) then
        do j=C%mx,C%px
          do k=C%my-1,C%my-1
            if (abs(F(j,k) - I%rank_neighbor) > 1e-12  ) test_passed = .false.
          end do
        end do
      endif

      call MPI_Reduce(test_passed,test_result,1,MPI_LOGICAL,MPI_LAND,0,I%comm,ierr)

  end subroutine exchange_interface_neighbors2d_test
  
  subroutine exchange_interface_neighbors3d_test(C,I,test_result)

      use mpi
      use mpi3d_interface
      use mpi3dcomm, only : cartesian

      implicit none

      type(cartesian),intent(in) :: C
      real,dimension(C%mbx:C%pbx,C%mby:C%pby,C%mbz:C%pbz) :: F
      type(interface3d),intent(in) :: I
      logical,intent(out) :: test_result
      integer :: j,k,l,nx,ny,nz,ierr
      logical :: test_passed


      test_passed = .true.


      F = real(0)



      if(.not. I%on_interface) return


      nx = I%normal(1)
      ny = I%normal(2)
      nz = I%normal(3)

      ! Fill all faces adjacent to ghost cells

      if (nx == 1 .and. ny == 0 .and. nz == 0) then
        do l=C%mz,C%pz
          do k=C%my,C%py
            do j=C%px,C%px
              F(j,k,l)  = 1.0 + real(k - C%my + (l - C%mz)*C%lny)
            end do
          end do
        end do
      endif  

      if (nx == -1 .and. ny == 0 .and. nz == 0) then
        do l=C%mz,C%pz
          do k=C%my,C%py
            do j=C%mx,C%mx
              F(j,k,l)  = 1.0 + real(k - C%my + (l - C%mz)*C%lny)
            end do
          end do
        end do
      endif  

      if (nx == 0 .and. ny == 1 .and. nz == 0) then
        do l=C%mz,C%pz
          do k=C%py,C%py
            do j=C%mx,C%px
              F(j,k,l)  = 1.0 + real(j - C%mx + (l - C%mz)*C%lnx)
            end do
          end do
        end do
      endif  

      if (nx == 0 .and. ny == -1 .and. nz == 0) then
        do l=C%mz,C%pz
          do k=C%my,C%my
            do j=C%mx,C%px
              F(j,k,l)  = 1.0 + real(j - C%mx + (l - C%mz)*C%lnx)
            end do
          end do
        end do
      endif  

      if (nx == 0 .and. ny == 0 .and. nz == 1) then
        do l=C%pz,C%pz
          do k=C%my,C%py
            do j=C%mx,C%px
              F(j,k,l)  = 1.0 + real(j - C%mx + (k - C%my)*C%lnx)
            end do
          end do
        end do
      endif                      
                    
      if (nx == 0 .and. ny == 0 .and. nz == -1) then
        do l=C%mz,C%mz
          do k=C%my,C%py
            do j=C%mx,C%px
              F(j,k,l)  = 1.0 + real(j - C%mx + (k - C%my)*C%lnx)
            end do
          end do
        end do
      endif                      

      call exchange_interface_neighbors3d(C,F,I)

      ! Check that ghost cells match neighboring boundary

      if (nx == 1 .and. ny == 0 .and. nz == 0) then
        do l=C%mz,C%pz
          do k=C%my,C%py
            do j=C%px,C%px
              if (F(j,k,l) /= F(j+1,k,l)) test_passed = .false. 
              print *,F(j,k,l),F(j+1,k,l)
            end do
          end do
        end do
      endif  
      
      if (nx == -1 .and. ny == 0 .and. nz == 0) then
        do l=C%mz,C%pz
          do k=C%my,C%py
            do j=C%mx,C%mx
              if (F(j,k,l) /= F(j-1,k,l)) test_passed = .false. 
              print *,F(j,k,l),F(j-1,k,l)
            end do
          end do
        end do
      endif  
      
      if (nx == 0 .and. ny == 1 .and. nz == 0) then
        do l=C%mz,C%pz
          do k=C%py,C%py
            do j=C%mx,C%px
              if (F(j,k,l) /= F(j,k+1,l)) test_passed = .false. 
            end do
          end do
        end do
      endif  

      if (nx == 0 .and. ny == -1 .and. nz == 0) then
        do l=C%mz,C%pz
          do k=C%my,C%my
            do j=C%mx,C%px
              if (F(j,k,l) /= F(j,k-1,l)) test_passed = .false. 
            end do
          end do
        end do
      endif  

      if (nx == 0 .and. ny == 0 .and. nz == 1) then
        do l=C%pz,C%pz
          do k=C%my,C%py
            do j=C%mx,C%px
              if (F(j,k,l) /= F(j,k,l+1)) test_passed = .false. 
            end do
          end do
        end do
      endif

      if (nx == 0 .and. ny == 0 .and. nz == -1) then
        do l=C%mz,C%mz
          do k=C%my,C%py
            do j=C%mx,C%px
              if (F(j,k,l) /= F(j,k,l-1))test_passed = .false. 
            end do
          end do
        end do
      endif

      call MPI_Reduce(test_passed,test_result,1,MPI_LOGICAL,MPI_LAND,0,I%comm,ierr)

  end subroutine exchange_interface_neighbors3d_test
end module mpi3d_interface_test
