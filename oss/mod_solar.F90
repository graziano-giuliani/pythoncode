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
module mod_solar

  use netcdf
  use mpi
  use iso_fortran_env

  implicit none

  private

  integer , parameter :: namelen = 16

  type asolar
    integer :: nf
    real , dimension(:) , allocatable :: frq
    real , dimension(:) , allocatable :: irr
    real , dimension(:) , allocatable :: ivals
    real :: f1
    real :: df
    logical :: isinit = .false.
    contains
      procedure , pass :: init_file
      procedure , pass :: init_file_mpi
      procedure , pass :: init_data
      procedure , pass :: interpolate
      procedure , pass :: delete
      ! final :: cleanup
  end type asolar

  public :: asolar

  character(len=namelen) , parameter :: fdim = 'nspectral'
  character(len=namelen) , parameter :: frq = 'FREQ'
  character(len=namelen) , parameter :: irr = 'IRRADIANCE'

  contains

    integer function init_file(this,datafile) result(iret)
      implicit none
      class(asolar) :: this
      character(len=*) , intent(in) :: datafile
      integer :: ncid , dimid , varid , itmp

      if ( this%isinit ) then
        iret = this%delete()
      end if

      iret = nf90_open(datafile,nf90_nowrite,ncid)
      if ( iret /= nf90_noerr ) then
        write(error_unit,*) 'In file ',__FILE__,' at line ',__LINE__, &
          ' : Cannot open file ',trim(datafile)
        write(error_unit,*) nf90_strerror(iret)
        return
      end if

      iret = nf90_inq_dimid(ncid,fdim,dimid)
      if ( iret /= nf90_noerr ) then
        write(error_unit,*) 'In file ',__FILE__,' at line ',__LINE__, &
          ' : No dimension ',trim(fdim),' in file ', trim(datafile)
        write(error_unit,*) nf90_strerror(iret)
        return
      end if

      iret = nf90_inquire_dimension(ncid,dimid,len=this%nf)
      if ( iret /= nf90_noerr ) then
        write(error_unit,*) 'In file ',__FILE__,' at line ',__LINE__, &
          ' : Cannot read dimension ',trim(fdim),' from file ', trim(datafile)
        write(error_unit,*) nf90_strerror(iret)
        return
      end if

      allocate(this%frq(this%nf), &
               this%irr(this%nf), stat=iret)
      if ( iret /= 0 ) then
        this%nf = 0
        write(error_unit,*) 'In file ',__FILE__,' at line ',__LINE__, &
          ' : Memory error allocating'
        return
      end if

      iret = nf90_inq_varid(ncid,frq,varid)
      if ( iret /= nf90_noerr ) then
        write(error_unit,*) 'In file ',__FILE__,' at line ',__LINE__, &
          ' : No variable ',trim(frq), ' in file ', trim(datafile)
        write(error_unit,*) nf90_strerror(iret)
        itmp = this%delete()
        return
      end if

      iret = nf90_get_var(ncid,varid,this%frq)
      if ( iret /= nf90_noerr ) then
        write(error_unit,*) 'In file ',__FILE__,' at line ',__LINE__, &
          ' : Cannot read variable ',trim(frq),' from file ', trim(datafile)
        write(error_unit,*) nf90_strerror(iret)
        itmp = this%delete()
        return
      end if

      iret = nf90_inq_varid(ncid,irr,varid)
      if ( iret /= nf90_noerr ) then
        write(error_unit,*) 'In file ',__FILE__,' at line ',__LINE__, &
          ' : No variable ',trim(irr), ' in file ', trim(datafile)
        write(error_unit,*) nf90_strerror(iret)
        itmp = this%delete()
        return
      end if

      iret = nf90_get_var(ncid,varid,this%irr)
      if ( iret /= nf90_noerr ) then
        write(error_unit,*) 'In file ',__FILE__,' at line ',__LINE__, &
          ' : Cannot read variable ',trim(irr),' from file ', trim(datafile)
        write(error_unit,*) nf90_strerror(iret)
        itmp = this%delete()
        return
      end if

      iret = nf90_close(ncid)
      if ( iret /= nf90_noerr ) then
        write(error_unit,*) 'In file ',__FILE__,' at line ',__LINE__, &
         ' : Cannot close file ', trim(datafile)
        write(error_unit,*) nf90_strerror(iret)
        itmp = this%delete()
        return
      end if

      this%f1 = this%frq(1)
      this%df = this%frq(2) - this%frq(1)
      this%isinit = .true.
      iret = 0
    end function init_file

    integer function init_data(this,frq,irr) result(iret)
      implicit none
      class(asolar) :: this
      real , dimension(:) , intent(in) :: frq
      real , dimension(:) , intent(in) :: irr

      if ( this%isinit ) then
        iret = this%delete()
      end if

      if ( size(frq) /= size(irr) ) then
        write(error_unit,*) 'In file ',__FILE__,' at line ',__LINE__, &
          ' : Dimension mismatch'
        iret = -1
        return
      end if

      this%nf = size(frq)
      allocate(this%frq(this%nf), this%irr(this%nf), stat=iret)
      if ( iret /= 0 ) then
        this%nf = 0
        write(error_unit,*) 'In file ',__FILE__,' at line ',__LINE__, &
          ' : Memory error allocating'
        return
      end if

      this%frq = frq
      this%irr = irr
      this%f1 = this%frq(1)
      this%df = this%frq(2) - this%frq(1)
      this%isinit = .true.
      iret = 0
    end function init_data

    integer function init_file_mpi(this,datafile,iocpu,icomm) result(iret)
      implicit none
      class(asolar) :: this
      character(len=*) , intent(in) :: datafile
      integer , intent(in) :: iocpu , icomm
      integer :: irank , mpierr

      if ( this%isinit ) then
        iret = this%delete()
      end if

      call mpi_comm_rank(icomm,irank,iret)
      if ( iret /= 0 ) then
        write(error_unit,*) 'In file ',__FILE__,' at line ',__LINE__, &
          ' : Error in MPI_COMM_RANK'
        return
      end if

      if ( irank == iocpu ) then
        iret = init_file(this,datafile)
        if ( iret /= 0 ) then
          write(error_unit,*) 'In file ',__FILE__,' at line ',__LINE__, &
            ' : IOCPU is not able to initialize'
          call mpi_abort(icomm,-1,iret)
          if ( iret /= 0 ) then
            write(error_unit,*) 'In file ',__FILE__,' at line ',__LINE__, &
              ' : THIS MUST NEVER HAPPEN !'
            call abort
          end if
        end if
      end if

      call mpi_bcast(this%nf,1,mpi_integer,iocpu,icomm,mpierr)
      if ( mpierr /= 0 ) then
        if ( irank == iocpu ) iret = this%delete( )
        iret = mpierr
        return
      end if

      if ( irank /= iocpu ) then
        allocate(this%frq(this%nf), this%irr(this%nf), stat=iret)
        if ( iret /= 0 ) then
          write(error_unit,*) 'In file ',__FILE__,' at line ',__LINE__, &
            ' : Allocation error on CPU ',irank
          call mpi_abort(icomm,-1,iret)
          if ( iret /= 0 ) then
            write(error_unit,*) 'In file ',__FILE__,' at line ',__LINE__, &
              ' : THIS MUST NEVER HAPPEN !'
            call abort
          end if
        end if
      end if

      call mpi_bcast(this%frq,this%nf,mpi_real,iocpu,icomm,mpierr)
      if ( mpierr /= 0 ) then
        iret = this%delete( )
        iret = mpierr
      end if

      call mpi_bcast(this%irr,this%nf,mpi_real,iocpu,icomm,mpierr)
      if ( mpierr /= 0 ) then
        iret = this%delete( )
        iret = mpierr
      end if

      if ( irank /= iocpu ) then
        this%f1 = this%frq(1)
        this%df = this%frq(2) - this%frq(1)
        this%isinit = .true.
      end if

      iret = 0
    end function init_file_mpi

    integer function interpolate(this,f) result(iret)
      implicit none
      class(asolar) :: this
      real , dimension(:) , intent(in) :: f
      real :: xf , w1 , w2
      integer :: i , i1 , i2 , nif

      if ( .not. this%isinit ) then
        write(error_unit,*) 'In file ',__FILE__,' at line ',__LINE__, &
          ' : Object not initialized'
        iret = -1
        return
      end if

      nif = size(f)
      if ( allocated(this%ivals) ) then
        if ( size(this%ivals) /= nif ) then
          deallocate(this%ivals)
          allocate(this%ivals(nif))
        end if
      else
        allocate(this%ivals(nif))
      end if
      do i = 1 , nif
        xf = (f(i)-this%f1)/this%df
        i1 = floor(xf)
        i2 = i1 + 1
        i1 = max(1,min(this%nf,i1))
        i2 = max(1,min(this%nf,i2))
        w1 = xf-real(i1)
        w2 = real(i2)-xf
        this%ivals(i) = this%irr(i1) * w2 + this%irr(i2) * w1
      end do
    end function interpolate

    integer function delete(this) result(iret)
      implicit none
      class(asolar) :: this
      if ( allocated(this%irr) ) then
        deallocate(this%frq)
        deallocate(this%irr)
        if ( allocated(this%ivals) ) then
          deallocate(this%ivals)
        end if
        this%nf = -1
        this%f1 = huge(this%f1)
        this%df = huge(this%df)
      end if
      this%isinit = .false.
      iret = 0
    end function delete

    subroutine cleanup(as)
      implicit none
      type(asolar) :: as
      if ( allocated(as%irr) ) then
        deallocate(as%frq)
        deallocate(as%irr)
        if ( allocated(as%ivals) ) then
          deallocate(as%ivals)
        end if
      end if
    end subroutine cleanup

end module mod_solar

#ifdef TESTME

! mpif90 -DTESTME -O3 -g -mtune=native -o test_solar \
!        -I`nf-config --includedir` mod_solar.F90 `nf-config --flibs`

program test_solar
  use mod_solar
  use mpi
  implicit none

  integer :: i , j , ni , icomm , myid , nproc
  real , dimension(1550) :: f
  type(asolar) :: as

  call mpi_init(i)
  if ( i /= 0 ) then
    print *, 'Cannot test mpi part'
  end if

  icomm = MPI_COMM_WORLD

  call mpi_comm_rank(icomm,myid,i)
  call mpi_comm_size(icomm,nproc,i)

  if ( myid == 0 ) then
    i = as%init_file('../data/solar_irradiances.nc')
    do ni = 0 , 999
      do j = 1 , 1550
        f(j) = 1550.0+ni + real(j-1)/(1550.0/1000.0)
      end do
      i = as%interpolate(f)
      print *, f(1) , f(1550)
      print *, as%ivals(1) , as%ivals(1550)
    end do
    i = as%delete()
  end if

  i = as%init_file_mpi('../data/solar_irradiances.nc',0,icomm)

  print *, 'MYID : ', myid , ', NF = ',as%nf , 'IRR = ',as%irr(as%nf/2)

  i = as%delete()

  call mpi_finalize(i)

end program test_solar
#endif
