module MPI_parallel
    implicit none 

    integer :: MPI_rank, MPI_size
    integer :: ierr
    character(len=5)  :: chr_rank
    integer, parameter :: root_rank = 0

end module MPI_parallel