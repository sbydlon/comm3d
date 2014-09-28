module mpi3dcomm_test

      use mpi3dcomm, only : cartesian

      implicit none

contains 
    

      subroutine exchange_test

          implicit none

          real,dimension(:,:,:),allocatable :: F1,F2
          integer,parameter :: nx1=21,ny1=11,nz1=21,nx2=21,ny2=11,nz2=21,nb=3,nF=6


end module mpi3dcomm_test

