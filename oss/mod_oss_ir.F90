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
module mod_oss_ir

  use mpi
  use mod_solar
  use mod_hitran
  use iso_fortran_env

  implicit none

  private

  ! 2010 CODATA
  !real , parameter :: h = 6.62606957E-27
  !real , parameter :: c = 2.99792458E+10
  !real , parameter :: b = 1.380662E-16
  !real , parameter :: c1 = 2.0E0*h*c*c
  !real , parameter :: c2 = h*c/b
  !real , parameter :: pi = 3.1415926535897932384626433832795029E0
  !real , parameter :: deg2rad = pi/180.0E0
  !real , parameter :: grav = 9.80665E0
  !real , parameter :: rair = 10.0E0/grav
  ! Original Parameters
  real , parameter :: h = 6.626176E-27
  real , parameter :: c = 2.997925E+10
  real , parameter :: b = 1.380662E-16
  real , parameter :: c1 = 2.0E0*h*c*c
  real , parameter :: c2 = h*c/b
  real , parameter :: pi = 3.14159265
  real , parameter :: deg2rad = 0.017453293E0
  real , parameter :: rair = 1.020408163E0

  ! Weighting factor derivative wrt optical depth, for linear-in-tau,
  ! in low-OD limit
  real , parameter :: dvint = 10.0E0
  real , parameter :: epsiln = 1.0E-05

  type aoss_ir
    real , dimension(:) , allocatable :: pavlref
    real , dimension(:) , allocatable :: sunrad
    type(asolar) , pointer :: as
    type(ahitran) , pointer :: ah
    logical :: isinit = .false.
    contains
      procedure , pass :: init
      procedure , pass :: delete
  end type aoss_ir

  public :: aoss_ir

  contains

    integer function init(this,as,ah) result(iret)
      implicit none
      class(aoss_ir) :: this
      type(asolar) , target :: as
      type(ahitran) , target :: ah
      integer :: itmp

      if ( this%isinit ) then
        iret = this%delete()
      end if

      ! Check input

      if ( .not. as%isinit ) then
        write(error_unit,*) 'In file ',__FILE__,' at line ',__LINE__, &
          ' : Not initialized solar radiances given !'
        iret = -1
        return
      end if
      if ( .not. ah%isinit ) then
        write(error_unit,*) 'In file ',__FILE__,' at line ',__LINE__, &
          ' : Not initialized hitran coefficients given !'
        iret = -1
        return
      end if

      allocate(this%pavlref(ah%nlayod),this%sunrad(ah%nf),stat=iret)
      if ( iret /= 0 ) then
        write(error_unit,*) 'In file ',__FILE__,' at line ',__LINE__, &
          ' : Allocation error !'
        call zerodims(this)
        iret = -1
        return
      end if

      this%pavlref(1:ah%nlev-1) = 0.5E0 * &
                (ah%pref(1:ah%nlev-1)+ah%pref(2:ah%nlev))

      iret = as%interpolate(ah%vwvn)
      if ( iret /= 0 ) then
        write(error_unit,*) 'In file ',__FILE__,' at line ',__LINE__, &
          ' : interpolation error in asolar type !'
        itmp = this%delete( )
        return
      end if
      this%sunrad = as%ivals

      this%as => as
      this%ah => ah

      this%isinit = .true.
      iret = 0
    end function init

    integer function delete(this) result(iret)
      implicit none
      class(aoss_ir) :: this
      call zerodims(this)
      call cleanup(this)
      this%isinit = .false.
      iret = 0
    end function delete

    subroutine zerodims(ao)
      implicit none
      type(aoss_ir) :: ao
    end subroutine zerodims

    subroutine cleanup(ao)
      implicit none
      type(aoss_ir) :: ao
      if ( allocated(ao%pavlref) ) then
        deallocate(ao%pavlref)
        deallocate(ao%sunrad)
      end if
    end subroutine cleanup

end module mod_oss_ir

#ifdef TESTME

! mpif90 -O3 -g -mtune=native -c -I`nf-config --includedir` mod_solar.F90
! mpif90 -O3 -g -mtune=native -c -I`nf-config --includedir` mod_hitran.F90
! mpif90 -DTESTME -O3 -g -mtune=native -o test_oss_ir \
!        -I`nf-config --includedir` mod_oss_ir.F90 `nf-config --flibs` \
!        mod_solar.o mod_hitran.o

program test_oss_ir
  use mod_solar
  use mod_hitran
  use mod_oss_ir
  use mpi
  implicit none

  integer :: i , j , ni , icomm , myid , nproc
  real , dimension(1550) :: f
  type(asolar) , target :: as
  type(ahitran) , target :: ah
  type(aoss_ir) :: ao

  call mpi_init(i)
  if ( i /= 0 ) then
    print *, 'Cannot test mpi part'
  end if

  icomm = MPI_COMM_WORLD

  call mpi_comm_rank(icomm,myid,i)
  call mpi_comm_size(icomm,nproc,i)

  i = as%init_file_mpi('../data/solar_irradiances.nc',0,icomm)
  i = ah%init_file_mpi('../data/leo.iasi.0.05.nc',0,icomm)

  i = ao%init(as,ah)
  if ( i /= 0 ) then
    stop
  end if

  print *, 'MYID : ', myid , ', NF = ',ao%as%nf , 'IRR = ',ao%sunrad(ah%nf/2)

  i = as%delete()
  i = ah%delete()
  i = ao%delete()

  call mpi_finalize(i)

end program test_oss_ir
#endif
