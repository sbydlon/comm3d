module mpi3dcomm

  implicit none

  ! Cartesian communicator type, including
  ! COMM = communicator
  ! RANK = rank of process
  ! SIZE = size of communicator
  ! COORD = integer coordinates of process
  ! RANK_M = rank of neighbor in minus direction
  ! RANK_P = rank of neighbor in plus  direction

  type :: cartesian
     integer :: nb, &
          nx,mx,px,mbx,pbx,lnx, &
          ny,my,py,mby,pby,lny, &
          nz,mz,pz,mbz,pbz,lnz
     integer :: comm,rank,size_x,size_y,size_z,coord(3), & ! maybe turn size_x, etc. into vector: size(3)
          rank_mx,rank_px,rank_my,rank_py,rank_mz,rank_pz ! maybe turn rank_mx, etc. into vectors: rank_m(3), rank_p(3)
      integer :: line_x,line_y,line_z,block_x,block_y,block_z,block3d_xy,block3d_xz,block3d_yz ! MPI types for communication, again maybe organized in vectors?
     integer :: array_w,array_s ! MPI types for I/O (leave here, Eric will use in I/O module)
  end type cartesian

  ! naming convention on indices, example for x direction:
  ! nb = number of additional boundary points (ghost nodes)
  ! nx = total number of points (on all processes, excluding ghost nodes)
  ! mx = starting (minus) index on process
  ! px = ending   (plus ) index on process
  ! mbx = starting index on process, including ghost nodes
  ! pbx = ending   index on process, including ghost nodes
  ! lnx = total number of points (on process)

  ! module procedures for allocating arrays
  interface allocate_array_body
     module procedure allocate_array_body_2d,allocate_array_body_2d_nf,allocate_array_body_3d,allocate_array_body_3d_nF
  end interface


contains
  

! 3D decomposition subroutine

   subroutine decompose3d(C,nx,ny,nz,nb,nF,periodic_x,periodic_y,periodic_z,comm, &
       method,size_x_in,size_y_in,size_z_in)

    use mpi3dbasic, only : is_master,decompose1d,error,MPI_REAL_PW,pw
    use mpi

    implicit none

    type(cartesian),intent(out) :: C
    integer,intent(in) :: nx,ny,nz,nb,nF,comm,size_x_in,size_y_in,size_z_in
    logical,intent(in) :: periodic_x,periodic_y,periodic_z
    character(*),intent(inout) :: method

    integer,parameter :: dim=3
    integer :: ierr,size,index,shift,size_xyz(3),l,count_xy,stride_xy,length_xy
    logical :: periodic(3),reorder
    integer,dimension(:),allocatable :: blocklengths,types
    integer,dimension(:),allocatable :: displacements_mpi1 ! MPI-1
    !integer(MPI_ADDRESS_KIND),dimension(:),allocatable :: displacements ! MPI-2

    ! domain size

    C%nx = nx
    C%ny = ny
    C%nz = nz
    C%nb = nb

    ! properties of original communicator

    call MPI_Comm_size(comm,size,ierr)

    ! processor layout in Cartesian topology

    if (method=='manual') then

       ! manual decomposition, user specifies number of processes in x and y directions

       if (size/=size_x_in*size_y_in*size_z_in.and.is_master(comm)) &
            call error('Error: Incorrect number of processors for manual decomposition','decompose3d')
       size_xyz = (/ size_x_in,size_y_in,size_z_in /)

    else

       select case(method)
       case default ! 3D: MPI_Dims_create selects number of processes in all directions
          size_xyz = 0
       case('1Dx') ! 1D, x: MPI_Dims_create selects number of processes in x direction only
          size_xyz(1) = 0
          size_xyz(2) = 1
          size_xyz(3) = 1 
       case('1Dy') ! 1D, y: MPI_Dims_create selects number of processes in y direction only
          size_xyz(1) = 1
          size_xyz(2) = 0
          size_xyz(3) = 1 
       case('1Dz') ! 1D, z: MPI_Dims_create selects number of processes in z direction only
          size_xyz(1) = 1
          size_xyz(2) = 1 
          size_xyz(3) = 0 
       case('2Dxy') ! 2D, xy: MPI_Dims_create selects number of processes in x and y directions only
          size_xyz(1) = 0
          size_xyz(2) = 0 
          size_xyz(3) = 1    
       case('2Dxz') ! 2D, xz: MPI_Dims_create selects number of processes in x and z directions only
          size_xyz(1) = 0
          size_xyz(2) = 1 
          size_xyz(3) = 0 
       case('2Dyz') ! 2D, yz: MPI_Dims_create selects number of processes in y and z directions only
          size_xyz(1) = 1
          size_xyz(2) = 0 
          size_xyz(3) = 0 
       end select

       ! special cases
       if (C%nx==1) size_xyz(1) = 1
       if (C%ny==1) size_xyz(2) = 1
       if (C%nz==1) size_xyz(3) = 1
       
       call MPI_Dims_create(size,dim,size_xyz,ierr)

    end if

    C%size_x = size_xyz(1) ! number of processes in x direction
    C%size_y = size_xyz(2) ! number of processes in y direction
    C%size_z = size_xyz(3) ! number of processes in z direction

    ! 3D Cartesian communicator and coordinates

    periodic = (/ periodic_x,periodic_y,periodic_z /)
    reorder = .true.
    call MPI_Cart_create(comm,dim,size_xyz,periodic,reorder,C%comm,ierr)
    call MPI_Comm_rank(C%comm,C%rank,ierr)
    call MPI_Cart_coords(C%comm,C%rank,dim,C%coord,ierr)

    ! nearest neighbors in x-direction

    index = 0; shift = 1
    call MPI_Cart_shift(C%comm,index,shift,C%rank_mx,C%rank_px,ierr)

    ! nearest neighbors in y-direction

    index = 1; shift = 1
    call MPI_Cart_shift(C%comm,index,shift,C%rank_my,C%rank_py,ierr)

    ! nearest neighbors in z-direction

    index = 2; shift = 1
    call MPI_Cart_shift(C%comm,index,shift,C%rank_mz,C%rank_pz,ierr)


    ! initial data distribution on processors

    call decompose1d(C%nx,C%size_x,C%coord(1),C%mx,C%px,C%lnx)
    call decompose1d(C%ny,C%size_y,C%coord(2),C%my,C%py,C%lny)
    call decompose1d(C%nz,C%size_z,C%coord(3),C%mz,C%pz,C%lnz) 
    C%mbx = C%mx-C%nb
    C%pbx = C%px+C%nb
    C%mby = C%my-C%nb
    C%pby = C%py+C%nb
    C%mbz = C%mz-C%nb
    C%pbz = C%pz+C%nb


    ! MPI types containing lines of constant x, y, and z

    ! variable x, constant y, constant z
    call MPI_Type_vector(C%lnx,1,1,MPI_REAL_PW,C%line_x,ierr)
    call MPI_Type_commit(C%line_x,ierr)

    ! variable y, constant x, constant z
    call MPI_Type_vector(C%lny,1,C%lnx+2*C%nb,MPI_REAL_PW,C%line_y,ierr)
    call MPI_Type_commit(C%line_y,ierr)

    ! variable z, constant x, constant y
    call MPI_Type_vector(C%lnz,1,(C%lnx+2*C%nb)*(C%lny+2*C%nb),MPI_REAL_PW,C%line_z,ierr)
    call MPI_Type_commit(C%line_z,ierr) 


    ! MPI types containing 2D boundary blocks

    allocate(blocklengths(C%nb),types(C%nb))
    blocklengths = 1
    allocate(displacements_mpi1(C%nb)) ! MPI-1
    !allocate(displacements(C%nb)) ! MPI-2

    !types = C%line_y
    !! MPI-1
    !displacements_mpi1 = (/ (l,l=0,C%nb-1) /)*(C%lnx+2*C%nb)*(C%lny+2*C%nb)
    !call MPI_Type_struct(C%nb,blocklengths,displacements_mpi1*pw,types,C%block_z,ierr)
    !call MPI_Type_commit(C%block_z,ierr)
    !! MPI-2
    !!displacements = (/ (l,l=0,C%nb-1) /)*(C%lnx+2*C%nb)
    !!call MPI_Type_create_struct(C%nb,blocklengths,displacements*pw,types,C%block_y,ierr)
        
    count_xy = C%lny
    stride_xy =  (C%lnx + 2*C%nb)
    length_xy = C%lnx
    call MPI_Type_vector(&
        count_xy,length_xy,stride_xy,MPI_REAL_PW,C%block_z,ierr)
    call MPI_Type_commit(C%block_z,ierr)     

    types = C%line_y
    ! MPI-1
    displacements_mpi1 = (/ (l,l=0,C%nb-1) /)
    call MPI_Type_struct(C%nb,blocklengths,displacements_mpi1*pw,types,C%block_x,ierr)
    call MPI_Type_commit(C%block_x,ierr)
    ! MPI-2
    !displacements = (/ (l,l=0,C%nb-1) /)
    !call MPI_Type_create_struct(C%nb,blocklengths,displacements*pw,types,C%block_x,ierr)

    types = C%line_x
    ! MPI-1
    displacements_mpi1 = (/ (l,l=0,C%nb-1) /)*(C%lnx+2*C%nb)
    call MPI_Type_struct(C%nb,blocklengths,displacements_mpi1*pw,types,C%block_y,ierr)
    call MPI_Type_commit(C%block_y,ierr)
    ! MPI-2
    !displacements = (/ (l,l=0,C%nb-1) /)*(C%lnx+2*C%nb)
    !call MPI_Type_create_struct(C%nb,blocklengths,displacements*pw,types,C%block_y,ierr)

    deallocate(blocklengths,types)
    deallocate(displacements_mpi1) ! MPI-1
    !deallocate(displacements) ! MPI-2


    ! MPI Types containing 3D boundary blocks


    ! 3D MPI block covering faces in x-y plane

    allocate(blocklengths(C%lnx),types(C%lnx))
    blocklengths = 1
    allocate(displacements_mpi1(C%lnx)) ! MPI-1
    !allocate(displacements(C%nb)) ! MPI-2 

    types = C%block_z
    ! MPI-1
    displacements_mpi1 = (/ (l,l=0,C%lnx-1) /)*(C%lnx+2*C%nb)*(C%lny+2*C%nb)
    call MPI_Type_struct(C%nb,blocklengths,displacements_mpi1*pw,types,&
        C%block3d_xy,ierr)
    call MPI_Type_commit(C%block3d_xy,ierr)
    ! MPI-2
    !displacements = (/ (l,l=0,C%lnx-1) /)
    !call MPI_Type_create_struct(C%lnx,blocklengths,displacements_mpi1*pw,types,&
    !   C%block3d_xy,ierr) 

    !types = C%block_z
    !! MPI-1
    !displacements_mpi1 = (/ (l,l=0,C%lnx-1) /)
    !call MPI_Type_struct(C%lnx,blocklengths,displacements_mpi1*pw,types,&
    !    C%block3d_xy,ierr)
    !call MPI_Type_commit(C%block3d_xy,ierr)
    !! MPI-2
    !!displacements = (/ (l,l=0,C%lnx-1) /)
    !!call MPI_Type_create_struct(C%lnx,blocklengths,displacements_mpi1*pw,types,&
    !!   C%block3d_xy,ierr)

    deallocate(blocklengths,types)
    deallocate(displacements_mpi1) ! MPI-1
    !deallocate(displacements) ! MPI-2 


    ! 3D MPI block covering faces in x-z plane

    allocate(blocklengths(C%lnz),types(C%lnz))
    blocklengths = 1
    allocate(displacements_mpi1(C%lnz)) ! MPI-1
    !allocate(displacements(C%lnz)) ! MPI-2 

    types = C%block_y
    ! MPI-1
    displacements_mpi1 = (/ (l,l=0,C%lnz-1) /)*(C%lnx+2*C%nb)*(C%lny+2*C%nb) 
    call MPI_Type_struct(C%lnz,blocklengths,displacements_mpi1*pw,types,&
        C%block3d_xz,ierr)
    call MPI_Type_commit(C%block3d_xz,ierr)
    ! MPI-2
    !displacements = (/ (l,l=0,C%lnz-1) /)*(C%lnx+2*C%nb)*(C%lny+2*C%nb)
    !call MPI_Type_create_struct(C%lnz,blocklengths,displacements*pw,types,&
    !   C%block3d_xz,ierr)

    deallocate(blocklengths,types)
    deallocate(displacements_mpi1) ! MPI-1
    !deallocate(displacements) ! MPI-2   


    ! 3D MPI block covering faces in x-y plane 

    allocate(blocklengths(C%lnz),types(C%lnz))
    blocklengths = 1
    allocate(displacements_mpi1(C%lnz)) ! MPI-1
    !allocate(displacements(C%lnz)) ! MPI-2 

    types = C%block_x
    ! MPI-1
    displacements_mpi1 = (/ (l,l=0,C%lnz-1) /)*(C%lnx+2*C%nb)*(C%lny+2*C%nb)
    call MPI_Type_struct(C%lnz,blocklengths,displacements_mpi1*pw,types,&
        C%block3d_yz,ierr)
    call MPI_Type_commit(C%block3d_yz,ierr)
    ! MPI-2
    !displacements = (/ (l,l=0,C%lnz-1) /)*(C%lnx+2*C%nb)*(C%lny+2*C%nb)
    !call MPI_Type_create_struct(C%lnz,blocklengths,displacements*pw,types,&
    !   C%block3d_yz,ierr)

    deallocate(blocklengths,types)
    deallocate(displacements_mpi1) ! MPI-1
    !deallocate(displacements) ! MPI-2 


    ! TO DO: ADD NEW DATA TYPE THAT PACKS ALL NF FIELDS TOGETHER FOR MORE EFFICIENT COMMUNICATION

  end subroutine decompose3d 

!2D Decomposition subroutine

  subroutine decompose2d(C,nx,ny,nb,nF,periodic_x,periodic_y,comm, &
       method,size_x_in,size_y_in)

    use mpi3dbasic, only : is_master,decompose1d,error,MPI_REAL_PW,pw
    use mpi

    implicit none

    type(cartesian),intent(out) :: C
    integer,intent(in) :: nx,ny,nb,nF,comm,size_x_in,size_y_in
    logical,intent(in) :: periodic_x,periodic_y
    character(*),intent(inout) :: method

    integer,parameter :: dim=2
    integer :: ierr,size,index,shift,size_xy(2),l
    logical :: periodic(2),reorder
    integer,dimension(:),allocatable :: blocklengths,types
    integer,dimension(:),allocatable :: displacements_mpi1 ! MPI-1
    !integer(MPI_ADDRESS_KIND),dimension(:),allocatable :: displacements ! MPI-2

    ! domain size

    C%nx = nx
    C%ny = ny
    C%nb = nb

    ! properties of original communicator

    call MPI_Comm_size(comm,size,ierr)

    ! processor layout in Cartesian topology

    if (method=='manual') then

       ! manual decomposition, user specifies number of processes in x and y directions

       if (size/=size_x_in*size_y_in.and.is_master(comm)) &
            call error('Error: Incorrect number of processors for manual decomposition','decompose2d')
       size_xy = (/ size_x_in,size_y_in /)

    else

       select case(method)
       case default ! 2D: MPI_Dims_create selects number of processes in both directions
          size_xy = 0
       case('1Dx') ! 1D, x: MPI_Dims_create selects number of processes in x direction only
          size_xy(1) = 0
          size_xy(2) = 1
       case('1Dy') ! 1D, y: MPI_Dims_create selects number of processes in y direction only
          size_xy(1) = 1
          size_xy(2) = 0
       end select

       ! special cases
       if (C%nx==1) size_xy(1) = 1
       if (C%ny==1) size_xy(2) = 1
       
       call MPI_Dims_create(size,dim,size_xy,ierr)

    end if

    C%size_x = size_xy(1) ! number of processes in x direction
    C%size_y = size_xy(2) ! number of processes in y direction

    ! 2D Cartesian communicator and coordinates

    periodic = (/ periodic_x,periodic_y /)
    reorder = .true.
    call MPI_Cart_create(comm,dim,size_xy,periodic,reorder,C%comm,ierr)
    call MPI_Comm_rank(C%comm,C%rank,ierr)
    call MPI_Cart_coords(C%comm,C%rank,dim,C%coord,ierr)

    ! nearest neighbors in x-direction

    index = 0; shift = 1
    call MPI_Cart_shift(C%comm,index,shift,C%rank_mx,C%rank_px,ierr)

    ! nearest neighbors in y-direction

    index = 1; shift = 1
    call MPI_Cart_shift(C%comm,index,shift,C%rank_my,C%rank_py,ierr)

    ! initial data distribution on processors

    call decompose1d(C%nx,C%size_x,C%coord(1),C%mx,C%px,C%lnx)
    call decompose1d(C%ny,C%size_y,C%coord(2),C%my,C%py,C%lny)
    C%mbx = C%mx-C%nb
    C%pbx = C%px+C%nb
    C%mby = C%my-C%nb
    C%pby = C%py+C%nb

    ! MPI types containing lines of constant x and y

    ! variable x, constant y
    call MPI_Type_vector(C%lnx,1,1,MPI_REAL_PW,C%line_x,ierr)
    call MPI_Type_commit(C%line_x,ierr)

    ! variable y, constant x
    call MPI_Type_vector(C%lny,1,C%lnx+2*C%nb,MPI_REAL_PW,C%line_y,ierr)
    call MPI_Type_commit(C%line_y,ierr)

    ! MPI types containing boundary blocks

    allocate(blocklengths(C%nb),types(C%nb))
    blocklengths = 1
    allocate(displacements_mpi1(C%nb)) ! MPI-1
    !allocate(displacements(C%nb)) ! MPI-2


    types = C%line_y
    ! MPI-1
    displacements_mpi1 = (/ (l,l=0,C%nb-1) /)
    call MPI_Type_struct(C%nb,blocklengths,displacements_mpi1*pw,types,C%block_x,ierr)
    call MPI_Type_commit(C%block_x,ierr)
    ! MPI-2
    !displacements = (/ (l,l=0,C%nb-1) /)
    !call MPI_Type_create_struct(C%nb,blocklengths,displacements*pw,types,C%block_x,ierr)


    types = C%line_x
    ! MPI-1
    displacements_mpi1 = (/ (l,l=0,C%nb-1) /)*(C%lnx+2*C%nb)
    call MPI_Type_struct(C%nb,blocklengths,displacements_mpi1*pw,types,C%block_y,ierr)
    call MPI_Type_commit(C%block_y,ierr)
    ! MPI-2
    !displacements = (/ (l,l=0,C%nb-1) /)*(C%lnx+2*C%nb)
    !call MPI_Type_create_struct(C%nb,blocklengths,displacements*pw,types,C%block_y,ierr)

    deallocate(blocklengths,types)
    deallocate(displacements_mpi1) ! MPI-1
    !deallocate(displacements) ! MPI-2

    ! TO DO: ADD NEW DATA TYPE THAT PACKS ALL NF FIELDS TOGETHER FOR MORE EFFICIENT COMMUNICATION

  end subroutine decompose2d


  subroutine exchange_all_neighbors(C,F)

    implicit none

    type(cartesian),intent(in) :: C
    real,dimension(C%mbx:C%pbx,C%mby:C%pby,C%mbz:C%pbz),intent(inout) :: F

    call exchange_neighbors(C,F,'xm')
    call exchange_neighbors(C,F,'xp')
    call exchange_neighbors(C,F,'ym')
    call exchange_neighbors(C,F,'yp')
    call exchange_neighbors(C,F,'zm')
    call exchange_neighbors(C,F,'zp')

  end subroutine exchange_all_neighbors
  

  subroutine exchange_neighbors(C,F,side)

    use mpi

    implicit none

    type(cartesian),intent(in) :: C
    real,dimension(C%mbx:C%pbx,C%mby:C%pby,C%mbz:C%pbz),intent(inout) :: F
    character(*),intent(in) :: side

    integer :: ierr
    integer,parameter :: tagxm=1, tagxp=2, tagym=3, tagyp=4, tagzm=5, tagzp=6

    ! share data between processes using blocking sendrecv with 2D communicator
       
    select case(side)
    case('xm') ! px --> mbx
       call MPI_SendRecv( &
            F(C%px-C%nb+1,C%my,C%mz),1,C%block3d_yz,C%rank_px,tagxm, &
            F(C%mx-C%nb  ,C%my,C%mz),1,C%block3d_yz,C%rank_mx,tagxm, &
            C%comm,MPI_STATUS_IGNORE,ierr)
    case('xp') ! pbx <-- mx
       call MPI_SendRecv( &
            F(C%mx  ,C%my,C%mz),1,C%block3d_yz,C%rank_mx,tagxp, &
            F(C%px+1,C%my,C%mz),1,C%block3d_yz,C%rank_px,tagxp, &
            C%comm,MPI_STATUS_IGNORE,ierr)
    case('ym') ! py --> mby
       call MPI_SendRecv( &
            F(C%mx,C%py-C%nb+1,C%mz),1,C%block3d_xz,C%rank_py,tagym, &
            F(C%mx,C%my-C%nb  ,C%mz),1,C%block3d_xz,C%rank_my,tagym, &
            C%comm,MPI_STATUS_IGNORE,ierr)
    case('yp') ! pby <-- my
       call MPI_SendRecv( &
            F(C%mx,C%my  ,C%mz),1,C%block3d_xz,C%rank_my,tagyp, &
            F(C%mx,C%py+1,C%mz),1,C%block3d_xz,C%rank_py,tagyp, &
            C%comm,MPI_STATUS_IGNORE,ierr)
    case('zm') ! pz --> mbz
       call MPI_SendRecv( &
            F(C%mx,C%my,C%pz-C%nb+1),1,C%block3d_xy,C%rank_my,tagzp, &
            F(C%mx,C%py,C%pz-1     ),1,C%block3d_xy,C%rank_py,tagzp, &
            C%comm,MPI_STATUS_IGNORE,ierr)
    case('zp') ! pz --> mbz
       call MPI_SendRecv( &
            F(C%mx,C%my,C%pz  ),1,C%block3d_xy,C%rank_my,tagzp, &
            F(C%mx,C%py,C%pz+1),1,C%block3d_xy,C%rank_py,tagzp, &
            C%comm,MPI_STATUS_IGNORE,ierr) 
    end select

  end subroutine exchange_neighbors
    

  subroutine test_exchange(C,F)

    implicit none

    type(cartesian),intent(in) :: C
    real,dimension(C%mbx:C%pbx,C%mby:C%pby,C%mbz:C%pbz),intent(inout) :: F

    integer :: i,j,k
    real :: G

    ! initialize array


    do j = C%my,C%py
       do i = C%mx,C%px
          do k = C%mz,C%pz
          F(i,j,k) = real(i)+1000d0*real(j)+1000000d0*real(k)
          end do
       end do
    end do

    ! exchange

    call exchange_all_neighbors(C,F)

    ! check

    ! left
    if (C%mx/=1) then
       do j = C%my,C%py
          do i = C%mbx,C%mx-1 
            do k = C%mz,C%pz
             G = real(i)+1000d0*real(j)+1000000d0*real(K)
             if (F(i,j,k)/=G) print *, i,j,k,G,F(i,j,k)
            end do
          end do
       end do
    end if

    ! right
    if (C%px/=C%nx) then
       do j = C%my,C%py
          do i = C%px+1,C%pbx
            do k = C%mz,C%pz
             G = real(i)+1000d0*real(j)+1000000d0*real(k)
             if (F(i,j,k)/=G) print *, i,j,k,G,F(i,j,k)
             end do
          end do
       end do
    end if

    ! bottom
    if (C%my/=1) then
       do j = C%mby,C%my-1
          do i = C%mx,C%px
            do k = C%mz,C%pz
             G = real(i)+1000d0*real(j)+1000000d0*real(k)
             if (F(i,j,k)/=G) print *, i,j,k,G,F(i,j,k)
            end do
          end do
       end do
    end if

    ! top
    if (C%py/=C%ny) then
       do j = C%py+1,C%pby
          do i = C%mx,C%px
            do k = C%mz,C%pz
             G = real(i)+1000d0*real(j)+1000000d0*real(k)
             if (F(i,j,k)/=G) print *, i,j,k,G,F(i,j,k)
            end do
          end do
       end do
    end if

        ! up (z-direction)
    if (C%my/=1) then
       do j = C%my,C%py
          do i = C%mx,C%px
            do k = C%mbz,C%mz-1
             G = real(i)+1000d0*real(j)+1000000d0*real(k)
             if (F(i,j,k)/=G) print *, i,j,k,G,F(i,j,k)
            end do
          end do
       end do
    end if

    ! down (z-direction)
    if (C%py/=C%ny) then
       do j = C%my,C%py
          do i = C%mx,C%px
            do k = C%pz+1,C%pbz
             G = real(i)+1000d0*real(j)+1000000d0*real(k)
             if (F(i,j,k)/=G) print *, i,j,k,G,F(i,j,k)
            end do
          end do
       end do
    end if 
  end subroutine test_exchange

    subroutine allocate_array_body_2d(F,C,ghost_nodes,Fval)

    implicit none

    real,dimension(:,:),allocatable,intent(inout) :: F
    type(cartesian),intent(in) :: C
    logical,intent(in) :: ghost_nodes
    real,intent(in),optional :: Fval

    if (allocated(F)) return

    if (ghost_nodes) then
       allocate(F(C%mbx:C%pbx,C%mby:C%pby))
    else
       allocate(F(C%mx :C%px ,C%my :C%py))
    end if
    if (present(Fval)) then
       F = Fval
    else
       F = 1d20
    end if

  end subroutine allocate_array_body_2d


  subroutine allocate_array_body_2d_nF(F,C,nF,ghost_nodes,Fval)

    implicit none

    real,dimension(:,:,:),allocatable,intent(inout) :: F
    type(cartesian),intent(in) :: C
    integer,intent(in) :: nF
    logical,intent(in) :: ghost_nodes
    real,intent(in),optional :: Fval

    if (allocated(F)) return

    if (ghost_nodes) then
       allocate(F(C%mbx:C%pbx,C%mby:C%pby,nF))
    else
       allocate(F(C%mx :C%px ,C%my :C%py ,nF))
    end if
    if (present(Fval)) then
       F = Fval
    else
       F = 1d20
    end if
 
  end subroutine allocate_array_body_2d_nF 



  subroutine allocate_array_body_3d(F,C,ghost_nodes,Fval)

    implicit none

    real,dimension(:,:,:),allocatable,intent(inout) :: F
    type(cartesian),intent(in) :: C
    logical,intent(in) :: ghost_nodes
    real,intent(in),optional :: Fval

    if (allocated(F)) return

    if (ghost_nodes) then
       allocate(F(C%mbx:C%pbx,C%mby:C%pby,C%mbz:C%pbz))
    else
       allocate(F(C%mx :C%px ,C%my :C%py ,C%mz :C%pz ))
    end if
    if (present(Fval)) then
       F = Fval
    else
       F = 1d20
    end if

  end subroutine allocate_array_body_3d


  subroutine allocate_array_body_3d_nF(F,C,nF,ghost_nodes,Fval)

    implicit none

    real,dimension(:,:,:,:),allocatable,intent(inout) :: F
    type(cartesian),intent(in) :: C
    integer,intent(in) :: nF
    logical,intent(in) :: ghost_nodes
    real,intent(in),optional :: Fval

    if (allocated(F)) return

    if (ghost_nodes) then
       allocate(F(C%mbx:C%pbx,C%mby:C%pby,C%mbz:C%pbz,nF))
    else
       allocate(F(C%mx :C%px ,C%my :C%py ,C%mz :C%pz ,nF))
    end if
    if (present(Fval)) then
       F = Fval
    else
       F = 1d20
    end if
 
  end subroutine allocate_array_body_3d_nF


end module mpi3dcomm
