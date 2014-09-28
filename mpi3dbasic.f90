module mpi3dbasic

  use mpi, only : MPI_REAL,MPI_DOUBLE_PRECISION

  implicit none

  integer,parameter ::  stdout=6, stderr=0 ! file units for stdout and stderr

  ! select precision for computation (usually double) and data output (usually single)
  !integer,parameter :: pw=4, MPI_REAL_PW=MPI_REAL             ! PW = working precision (single)
  integer,parameter :: pw=8, MPI_REAL_PW=MPI_DOUBLE_PRECISION ! PW = working precision (double)
  integer,parameter :: ps=4, MPI_REAL_PS=MPI_REAL             ! PS = saving  precision (single)
  !integer,parameter :: ps=8, MPI_REAL_PS=MPI_DOUBLE_PRECISION ! PS = saving  precision (double)

  ! properties of MPI_COMM_WORLD
  integer,save :: size,rank

  real,save :: time_start ! time at which MPI was initialized

contains

  
  subroutine start_mpi

    use mpi
    
    implicit none

    integer :: ierr
    
    ! initialize MPI, obtain size and rank in MPI_COMM_WORLD

    call MPI_Init(ierr)
    call MPI_Comm_size(MPI_COMM_WORLD,size,ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)

    time_start = MPI_WTime() ! start clock

  end subroutine start_mpi


  function time_elapsed() result(t)

    use mpi
    
    implicit none

    real :: t

    t = MPI_WTime()-time_start ! elapsed time since start_mpi() was called

  end function time_elapsed


  subroutine finish_mpi

    use mpi

    implicit none

    integer :: ierr

    if (is_master()) write(stdout,'(/,a,g20.10,a)') &
         'Total MPI time ',time_elapsed(),' s'
    
    call MPI_Finalize(ierr)

  end subroutine finish_mpi


  function is_master(comm)

    use mpi

    implicit none

    integer,intent(in),optional :: comm
    logical :: is_master

    integer :: rank,ierr

    if (present(comm)) then
       call MPI_Comm_rank(comm,rank,ierr)
    else
       call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
    end if

    is_master = (rank==0)

  end function is_master


  subroutine message(str,routine)
    ! MESSAGE writes informational message

    implicit none

    ! STR = string to be written to standard out
    ! ROUTINE =  subroutine name in which message originated

    character(*),intent(in) :: str
    character(*),intent(in),optional :: routine
   
    character(12) :: id
    character(256) :: str0,str1

    ! write message and subroutine name (if given)

    write(id,'(i10)') rank
    write(id,'(a)') '(' // trim(adjustl(id)) // ') '
    write(str1,'(a)') id // trim(adjustl(str))

    if (present(routine)) then
       write(str0,'(a)') id // 'Message from subroutine: ' &
            // trim(adjustl(routine))
       write(stdout,'(a,/,a)') trim(str0),trim(str1)
    else
       write(stdout,'(a)') trim(str1)
    end if

  end subroutine message


  subroutine warning(str,routine)
    ! WARNING writes error message, but does not terminate program

    implicit none

    ! STR = string to be written to standard error
    ! ROUTINE =  subroutine name in which error occurred

    character(*),intent(in) :: str
    character(*),intent(in),optional :: routine

    character(12) :: id
    character(256) :: str0,str1

    ! write error message and subroutine name (if given)
    
    write(id,'(i10)') rank
    write(id,'(a)') '(' // trim(adjustl(id)) // ') '
    write(str1,'(a)') id // trim(adjustl(str))

    if (present(routine)) then
       write(str0,'(a)') id // 'Warning in subroutine: ' &
            // trim(adjustl(routine))
       write(stderr,'(a,/,a)') trim(str0),trim(str1)
    else
       write(stderr,'(a)') trim(str1)
    end if

  end subroutine warning


  subroutine error(str,routine)
    ! ERROR writes error message and terminates program

    use mpi

    implicit none

    ! STR = string to be written to standard error
    ! ROUTINE =  subroutine name in which error occurred

    character(*),intent(in) :: str
    character(*),intent(in),optional :: routine
   
    ! IERR = MPI error flag

    integer :: ierr
    character(12) :: id
    character(256) :: str0,str1,str2

    ! write error message and subroutine name (if given)
    
    write(id,'(i10)') rank
    write(id,'(a)') '(' // trim(adjustl(id)) // ') '
    write(str1,'(a)') id // trim(adjustl(str))
    write(str2,'(a)') id // 'Terminating program'

    if (present(routine)) then
       write(str0,'(a)') id // 'Error in subroutine: ' &
            // trim(adjustl(routine))
       write(stderr,'(a,/,a,/,a)') trim(str0),trim(str1),trim(str2)
    else
       write(stderr,'(a,/,a)') trim(str1),trim(str2)
    end if

    ! terminate program

    call MPI_Abort(MPI_COMM_WORLD,0,ierr)

  end subroutine error


  subroutine new_communicator(in_new,comm_new,comm,new_rank)

    use mpi

    implicit none

    logical,intent(in) :: in_new
    integer,intent(out) :: comm_new
    integer,intent(in),optional :: comm, new_rank

    integer :: np,n,ierr,comm_old,group_old,group_new,size_new,rank
    integer,dimension(:),allocatable :: ranks,ranks_new

    ! creates new communicator (comm_new) from group of processes for which in_new=T,
    ! must be called by all processes in original communicator (comm)

    ! old communicator (defaults to MPI_COMM_WORLD)

    if (present(comm)) then
       comm_old = comm
    else
       comm_old = MPI_COMM_WORLD
    end if

    call MPI_Comm_size(comm_old,np,ierr)
    allocate(ranks(np),ranks_new(np))

    ! determine which processes will comprise new communicator

    if (in_new) then
        if (present(new_rank)) then
            rank = new_rank
        else
            call MPI_Comm_rank(comm_old,rank,ierr)
        end if
    else
       rank = -1
    end if

    call MPI_Allgather(rank,1,MPI_INTEGER,ranks,1,MPI_INTEGER,comm_old,ierr)

    size_new = 0
    do n = 1,np
       if (ranks(n)==-1) cycle
       size_new = size_new+1
       ranks_new(size_new) = ranks(n)
    end do

    ! create new communicator

    call MPI_Comm_group(comm_old,group_old,ierr)
    call MPI_Group_incl(group_old,size_new, &
         ranks_new(1:size_new),group_new,ierr)
    call MPI_Comm_create(comm_old,group_new,comm_new,ierr)
    if (in_new) call MPI_Group_free(group_new,ierr)
    call MPI_Group_free(group_old,ierr)

    deallocate(ranks,ranks_new)

  end subroutine new_communicator

  subroutine decompose1d(n,p,i,l,u,c)

    use mpi

    implicit none
    
    ! N = number of tasks
    ! P = number of processors
    ! I = rank of process (0<=i<=p-1)
    ! L = starting index of tasks (with tasks numbered 1:n not 0:n-1)
    ! U = ending index of tasks
    ! C = count, number of tasks assigned to process

    integer,intent(in) :: n,p,i
    integer,intent(out) :: l,u,c

    integer :: m,r

    if (n<p) call error('Error in decompose1d: number of tasks must not be less than number of processes','decompose1d')

    m = n/p ! integer division ignores any fractional part
    r = mod(n,p) ! remainder

    ! assign p-r processes m tasks each, r processes m+1 tasks each

    ! case 0, r=0

    if (r==0) then
       c = m
       l = 1+i*m
       u = l+c-1
       return
    end if

    ! case 1, p and r even OR p and r odd, symmetric task distribution
    
    if (mod(p,2)==mod(r,2)) then
       if (i<(p-r)/2) then ! low rank (0:(p-r)/2-1), m tasks
          c = m
          l = 1+i*m
       elseif (i<(p+r)/2) then ! intermediate rank ((p-r)/2:(p+r)/2-1), m+1 tasks
          c = m+1
          l = 1+(p-r)/2*m+(i-(p-r)/2)*(m+1)
       else ! high rank ((p+r)/2:p-1), m tasks
          c = m
          l = 1+(p-r)/2*m+r*(m+1)+(i-(p+r)/2)*m
       end if
    end if

    ! case 2, p odd and r even, symmetric task distribution
    
    if (mod(p,2)/=0.and.mod(r,2)==0) then
       if (i<r/2) then ! low rank (0:r/2-1), m+1 tasks
          c = m+1
          l = 1+i*(m+1)
       elseif (i<p-r/2) then ! intermediate rank (r/2:p-r/2-1), m tasks
          c = m
          l = 1+(r/2)*(m+1)+(i-r/2)*m
       else ! high rank (p-r/2:p-1), m+1 tasks
          c = m+1
          l = 1+(r/2)*(m+1)+(p-r)*m+(i-(p-r/2))*(m+1)
       end if
    end if

    ! case 3, p even and r odd, asymmetric task distribution

    if (mod(p,2)==0.and.mod(r,2)/=0) then
       if (i<p-r) then ! low rank (0:p-r-1), m tasks
          c = m
          l = 1+i*m
       else ! high rank (p-r:p-1), m+1 tasks
          c = m+1
          l = 1+(p-r)*m+(i-p+r)*(m+1)
       end if
    end if

    u = l+c-1

    if (c==0) call error('Error in decompose1d: at least one process has no data','decompose1d')

  end subroutine decompose1d


  subroutine test_decompose1d
    
    implicit none

    integer :: n,p,i,l,u,c

    do n = 14,14
       print *
       do p = 3,3
          print *
          do i = 0,p-1
             call decompose1d(n,p,i,l,u,c)
             write(6,'(6i6)') n,p,i,l,u,c
          end do
       end do
    end do

  end subroutine test_decompose1d


end module mpi3dbasic
