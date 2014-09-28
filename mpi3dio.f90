module mpi3dio

  implicit none

  ! derived type files distributed across multiple processes

  type :: file_distributed
     character(256) :: name
     integer :: fh,array,pr,MPI_REAL_PR
  end type file_distributed

  ! collective communication or not
  logical,parameter :: collective=.true.

  ! module procedures

  interface subarray
     module procedure subarray1d,subarray2d,subarray3d
  end interface

  interface write_file_distributed
     module procedure write_file_distributed_0d,write_file_distributed_1d, &
          write_file_distributed_2d,write_file_distributed_3d
  end interface

  interface read_file_distributed
     module procedure read_file_distributed_0d,read_file_distributed_1d, &
          read_file_distributed_2d,read_file_distributed_3d
  end interface


contains


  subroutine init_io(C)

    use mpi3dbasic, only : MPI_REAL_PW,MPI_REAL_PS,pw,ps
    use mpi3dcomm, only : cartesian

    implicit none

    type(cartesian),intent(inout) :: C

    ! MPI types used to set file views for data I/O;
    ! set processor-specific offset and extent

    ! 3D array section

    call subarray(C%nx,C%ny,C%nz,C%mx,C%px,C%my,C%py,C%mz,C%pz,MPI_REAL_PW,C%array_w)
    call subarray(C%nx,C%ny,C%nz,C%mx,C%px,C%my,C%py,C%mz,C%pz,MPI_REAL_PS,C%array_s)

  end subroutine init_io


  subroutine subarray1d(n,m,p,precision,array)

    use mpi

    implicit none

    integer,intent(in) :: n,m,p,precision
    integer,intent(out) :: array

    integer,parameter :: dim = 1
    integer,dimension(dim) :: sizes,subsizes,starts
    integer :: ierr

    sizes = n
    subsizes = p-m+1
    starts = m-1 ! subtract 1 since MPI starting index is 0

    call MPI_Type_create_subarray(dim,sizes,subsizes,starts, &
         MPI_ORDER_FORTRAN,precision,array,ierr)
    call MPI_Type_commit(array,ierr)

  end subroutine subarray1d


  subroutine subarray2d(nx,ny,mx,px,my,py,precision,array)

    use mpi

    implicit none

    integer,intent(in) :: nx,ny,mx,px,my,py,precision
    integer,intent(out) :: array

    integer,parameter :: dim = 2
    integer,dimension(dim) :: sizes,subsizes,starts
    integer :: ierr

    sizes = (/ nx,ny /)
    subsizes = (/ px-mx+1,py-my+1 /)
    starts = (/ mx,my /)-1 ! subtract 1 since MPI starting index is 0

    call MPI_Type_create_subarray(dim,sizes,subsizes,starts, &
         MPI_ORDER_FORTRAN,precision,array,ierr)
    call MPI_Type_commit(array,ierr)

  end subroutine subarray2d


  subroutine subarray3d(nx,ny,nz,mx,px,my,py,mz,pz,precision,array)

    use mpi

    implicit none

    integer,intent(in) :: nx,ny,nz,mx,px,my,py,mz,pz,precision
    integer,intent(out) :: array

    integer,parameter :: dim = 3
    integer,dimension(dim) :: sizes,subsizes,starts
    integer :: ierr

    sizes = (/ nx,ny,nz /)
    subsizes = (/ px-mx+1,py-my+1,pz-mz+1 /)
    starts = (/ mx,my,mz /)-1 ! subtract 1 since MPI starting index is 0

    
    call MPI_Type_create_subarray(dim,sizes,subsizes,starts, &
         MPI_ORDER_FORTRAN,precision,array,ierr)
    call MPI_Type_commit(array,ierr)

  end subroutine subarray3d


  subroutine open_file_distributed(fid,name,operation,comm,array,precision,offset_in)

    use mpi3dbasic, only : pw,ps,MPI_REAL_PW,MPI_REAL_PS,error
    use mpi

    implicit none

    type(file_distributed),intent(out) :: fid
    character(*),intent(in) :: name,operation
    integer,intent(in) :: comm,array,precision
    integer(MPI_OFFSET_KIND),intent(in),optional :: offset_in

    integer :: ierr
    integer(MPI_OFFSET_KIND) :: offset
    integer(MPI_OFFSET_KIND),parameter :: zero=0
    character(256) :: str

    ! offset (displacement) at beginning of file

    if (present(offset_in)) then
       offset = offset_in
    else
       offset = zero
    end if

    ! file name

    fid%name = name

    ! precision

    if (precision==pw) then
       fid%pr = pw
       fid%MPI_REAL_PR = MPI_REAL_PW
    else
       fid%pr = ps
       fid%MPI_REAL_PR = MPI_REAL_PS
    end if

    ! open file
    
    select case(operation)
    case('read')
       call MPI_File_open(comm,fid%name,MPI_MODE_RDONLY, &
            MPI_INFO_NULL,fid%fh,ierr)
    case('write','append')
       call MPI_File_open(comm,fid%name,MPI_MODE_CREATE+MPI_MODE_WRONLY, &
            MPI_INFO_NULL,fid%fh,ierr)
    case default
       call error('Invalid file operation','open_file_distributed') 
    end select

    if (ierr/=MPI_SUCCESS) then
        call io_error_message('MPI_File_open (on file ' // trim(name) // ')',ierr,str)
        call error(str,'open_file_distributed') 
    end if

    ! if writing, delete contents of file (by setting size to zero)

    if (operation=='write') then
       call MPI_File_set_size(fid%fh,zero,ierr)
       if (ierr/=MPI_SUCCESS) then
          call io_error_message('MPI_File_set_size',ierr,str)
          call error(str,'open_file_distributed') 
       end if
    end if

    ! set view

    call MPI_File_set_view(fid%fh,offset,fid%MPI_REAL_PR,array,'native',MPI_INFO_NULL,ierr)
    if (ierr/=MPI_SUCCESS) then
       call io_error_message('MPI_File_set_view',ierr,str)
       call error(str,'open_file_distributed') 
    end if

  end subroutine open_file_distributed

  
  subroutine write_file_distributed_0d(fid,data)

    use mpi3dbasic, only : pw,ps,error
    use mpi

    implicit none
    
    type(file_distributed),intent(in) :: fid
    real,intent(in) :: data

    integer :: ierr
    character(256) :: str

    if (fid%pr==pw) then
       if (collective) then
          call MPI_File_write_all(fid%fh,data,1, &
               fid%MPI_REAL_PR,MPI_STATUS_IGNORE,ierr)
       else
          call MPI_File_write    (fid%fh,data,1, &
               fid%MPI_REAL_PR,MPI_STATUS_IGNORE,ierr)
       end if
    else
       if (collective) then
          call MPI_File_write_all(fid%fh,real(data,ps),1, &
               fid%MPI_REAL_PR,MPI_STATUS_IGNORE,ierr)
       else
          call MPI_File_write    (fid%fh,real(data,ps),1, &
               fid%MPI_REAL_PR,MPI_STATUS_IGNORE,ierr)
       end if
    end if

    if (ierr/=MPI_SUCCESS) then
       call io_error_message('MPI_File_write_all',ierr,str)
       call error(str,'write_file_distributed_0d')
    end if

  end subroutine write_file_distributed_0d


  subroutine write_file_distributed_1d(fid,data)

    use mpi3dbasic, only : pw,ps,error
    use mpi

    implicit none
    
    type(file_distributed),intent(in) :: fid
    real,dimension(:),intent(in) :: data

    integer :: ierr
    character(256) :: str

    if (fid%pr==pw) then
       if (collective) then
          call MPI_File_write_all(fid%fh,data,size(data), &
               fid%MPI_REAL_PR,MPI_STATUS_IGNORE,ierr)
       else
          call MPI_File_write    (fid%fh,data,size(data), &
               fid%MPI_REAL_PR,MPI_STATUS_IGNORE,ierr)
       end if
    else
       if (collective) then
          call MPI_File_write_all(fid%fh,real(data,ps),size(data), &
               fid%MPI_REAL_PR,MPI_STATUS_IGNORE,ierr)
       else
          call MPI_File_write    (fid%fh,real(data,ps),size(data), &
               fid%MPI_REAL_PR,MPI_STATUS_IGNORE,ierr)
       end if
    end if

    if (ierr/=MPI_SUCCESS) then
       call io_error_message('MPI_File_write_all',ierr,str)
       call error(str,'write_file_distributed_1d')
    end if

  end subroutine write_file_distributed_1d


  subroutine write_file_distributed_2d(fid,data)

    use mpi3dbasic, only : pw,ps,error
    use mpi

    implicit none
    
    type(file_distributed),intent(in) :: fid
    real,dimension(:,:),intent(in) :: data

    integer :: ierr
    character(256) :: str

    if (fid%pr==pw) then
       if (collective) then
          call MPI_File_write_all(fid%fh,data,size(data), &
               fid%MPI_REAL_PR,MPI_STATUS_IGNORE,ierr)
       else
          call MPI_File_write    (fid%fh,data,size(data), &
               fid%MPI_REAL_PR,MPI_STATUS_IGNORE,ierr)
       end if
    else
       if (collective) then
          call MPI_File_write_all(fid%fh,real(data,ps),size(data), &
               fid%MPI_REAL_PR,MPI_STATUS_IGNORE,ierr)
       else
          call MPI_File_write    (fid%fh,real(data,ps),size(data), &
               fid%MPI_REAL_PR,MPI_STATUS_IGNORE,ierr)
       end if
    end if

    if (ierr/=MPI_SUCCESS) then
       call io_error_message('MPI_File_write_all',ierr,str)
       call error(str,'write_file_distributed_2d')
    end if

  end subroutine write_file_distributed_2d


  subroutine write_file_distributed_3d(fid,data)

    use mpi3dbasic, only : pw,ps,error
    use mpi

    implicit none
    
    type(file_distributed),intent(in) :: fid
    real,dimension(:,:,:),intent(in) :: data

    integer :: ierr
    character(256) :: str


    if (fid%pr==pw) then
       if (collective) then
          call MPI_File_write_all(fid%fh,data,size(data), &
               fid%MPI_REAL_PR,MPI_STATUS_IGNORE,ierr)
       else
          call MPI_File_write    (fid%fh,data,size(data), &
               fid%MPI_REAL_PR,MPI_STATUS_IGNORE,ierr)
       end if
    else
       if (collective) then
          call MPI_File_write_all(fid%fh,real(data,ps),size(data), &
               fid%MPI_REAL_PR,MPI_STATUS_IGNORE,ierr)
       else
          call MPI_File_write    (fid%fh,real(data,ps),size(data), &
               fid%MPI_REAL_PR,MPI_STATUS_IGNORE,ierr)
       end if
    end if

    if (ierr/=MPI_SUCCESS) then
       call io_error_message('MPI_File_write_all',ierr,str)
       call error(str,'write_file_distributed_3d')
    end if

  end subroutine write_file_distributed_3d


  subroutine read_file_distributed_0d(fid,data)

    use mpi3dbasic, only : MPI_REAL_PW,error
    use mpi

    implicit none
    
    type(file_distributed),intent(in) :: fid
    real,intent(out) :: data

    integer :: ierr
    character(256) :: str

    if (collective) then
       call MPI_File_read_all(fid%fh,data,1,MPI_REAL_PW,MPI_STATUS_IGNORE,ierr)
    else
       call MPI_File_read    (fid%fh,data,1,MPI_REAL_PW,MPI_STATUS_IGNORE,ierr)
    end if
    if (ierr/=MPI_SUCCESS) then
       call io_error_message('MPI_File_read_all',ierr,str)
       call error(str,'read_file_distributed_0d')
    end if

  end subroutine read_file_distributed_0d


  subroutine read_file_distributed_1d(fid,data)

    use mpi3dbasic, only : MPI_REAL_PW,error
    use mpi

    implicit none
    
    type(file_distributed),intent(in) :: fid
    real,dimension(:),intent(out) :: data

    integer :: ierr
    character(256) :: str

    if (collective) then
       call MPI_File_read_all(fid%fh,data,size(data),MPI_REAL_PW,MPI_STATUS_IGNORE,ierr)
    else
       call MPI_File_read    (fid%fh,data,size(data),MPI_REAL_PW,MPI_STATUS_IGNORE,ierr)
    end if
    if (ierr/=MPI_SUCCESS) then
       call io_error_message('MPI_File_read_all',ierr,str)
       call error(str,'read_file_distributed_1d')
    end if

  end subroutine read_file_distributed_1d


  subroutine read_file_distributed_2d(fid,data)

    use mpi3dbasic, only : MPI_REAL_PW,error
    use mpi

    implicit none
    
    type(file_distributed),intent(in) :: fid
    real,dimension(:,:),intent(out) :: data

    integer :: ierr
    character(256) :: str

    if (collective) then
       call MPI_File_read_all(fid%fh,data,size(data),MPI_REAL_PW,MPI_STATUS_IGNORE,ierr)
    else
       call MPI_File_read    (fid%fh,data,size(data),MPI_REAL_PW,MPI_STATUS_IGNORE,ierr)
    end if
    if (ierr/=MPI_SUCCESS) then
       call io_error_message('MPI_File_read_all',ierr,str)
       call error(str,'read_file_distributed_2d')
    end if

  end subroutine read_file_distributed_2d


  subroutine read_file_distributed_3d(fid,data)

    use mpi3dbasic, only : MPI_REAL_PW,error
    use mpi

    implicit none
    
    type(file_distributed),intent(in) :: fid
    real,dimension(:,:,:),intent(out) :: data

    integer :: ierr
    character(256) :: str

    if (collective) then
       call MPI_File_read_all(fid%fh,data,size(data),MPI_REAL_PW,MPI_STATUS_IGNORE,ierr)
    else
       call MPI_File_read    (fid%fh,data,size(data),MPI_REAL_PW,MPI_STATUS_IGNORE,ierr)
    end if
    if (ierr/=MPI_SUCCESS) then
       call io_error_message('MPI_File_read_all',ierr,str)
       call error(str,'read_file_distributed_3d')
    end if

  end subroutine read_file_distributed_3d


  subroutine close_file_distributed(fid)

    use mpi3dbasic, only : error
    use mpi

    implicit none

    type(file_distributed),intent(inout) :: fid

    integer :: ierr
    character(256) :: str

    call MPI_File_close(fid%fh,ierr)
    if (ierr/=MPI_SUCCESS) then
       call io_error_message('MPI_File_close',ierr,str)
       call error(str,'close_file_distributed')
    end if

  end subroutine close_file_distributed


  subroutine io_error_message(routine,ierr,str)

    !use mpi

    implicit none

    character(*),intent(in) :: routine
    integer,intent(in) :: ierr
    character(*),intent(out) :: str

    !select case(ierr)
    !case(MPI_ERR_IO)
    !   write(str,'(a)') 'Problem with ' // trim(routine) // ': MPI_ERR_IO'
    !case(MPI_ERR_NO_SPACE)
    !   write(str,'(a)') 'Problem with ' // trim(routine) // ': MPI_ERR_NO_SPACE'
    !case(MPI_ERR_NO_SUCH_FILE)
    !   write(str,'(a)') 'Problem with ' // trim(routine) // ': MPI_ERR_NO_SUCH_FILE'
    !case default
       write(str,'(a,i6)') 'Problem with ' // trim(routine) // ': ierr=',ierr
    !end select

  end subroutine io_error_message


end module mpi3dio
