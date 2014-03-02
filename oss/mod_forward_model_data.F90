!
! Copyright (c) 2013 Paolo Antonelli, Tiziana Cherubini, Graziano Giuliani
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.
!
module mod_forward_model_data

  use mpi
  use iso_fortran_env

  implicit none

  private

  type atmospheric_data
    real :: skt
    real :: ps
    real , dimension(:) , allocatable :: t
    real , dimension(:) , allocatable :: q
    real , dimension(:) , allocatable :: o3
    real , dimension(:) , allocatable :: co2
  end type atmospheric_data

  type synthetic_measure
    real , dimension(:) , allocatable :: y
    real , dimension(:,:) , allocatable :: xkt
    real , dimension(:,:) , allocatable :: xkemrf
    real , dimension(:,:) , allocatable :: paxkemrf
  end type synthetic_measure

  public :: atmospheric_data
  public :: synthetic_measure

  contains

end module mod_forward_model_data

#ifdef TESTME

program test_forward_model_data
  use mod_forward_model
  use mpi
  implicit none

  integer :: icomm , myid , nproc

  call mpi_init(i)
  if ( i /= 0 ) then
    print *, 'Cannot test mpi part'
  end if

  icomm = MPI_COMM_WORLD

  call mpi_comm_rank(icomm,myid,i)
  call mpi_comm_size(icomm,nproc,i)

  ! Do something

  call mpi_finalize(i)

end program test_forward_model_data
#endif
