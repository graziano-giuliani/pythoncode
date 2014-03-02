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
module mod_emissivity

  use netcdf
  use mpi
  use iso_fortran_env

  implicit none

  private

  integer , parameter :: namelen = 16

  type aemissivity
    integer :: nf
    real , dimension(:) , allocatable :: sfgrd
    real , dimension(:) , allocatable :: emrf
    logical :: isinit = .false.
    contains
      procedure , pass :: init_file
      procedure , pass :: init_file_mpi
      procedure , pass :: init_data
      procedure , pass :: delete
      ! final :: cleanup
  end type aemissivity

  public :: aemissivity

  character(len=namelen) , parameter :: fdim = 'nspectral'
  character(len=namelen) , parameter :: sfgrd = 'SfGrd'
  character(len=namelen) , parameter :: emrf = 'EmRf'

  contains

    integer function init_file(this,datafile) result(iret)
      implicit none
      class(aemissivity) :: this
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

      allocate(this%sfgrd(this%nf), &
               this%emrf(this%nf), stat=iret)
      if ( iret /= 0 ) then
        this%nf = 0
        write(error_unit,*) 'In file ',__FILE__,' at line ',__LINE__, &
          ' : Memory error allocating'
        return
      end if

      iret = nf90_inq_varid(ncid,sfgrd,varid)
      if ( iret /= nf90_noerr ) then
        write(error_unit,*) 'In file ',__FILE__,' at line ',__LINE__, &
          ' : No variable ',trim(sfgrd), ' in file ', trim(datafile)
        write(error_unit,*) nf90_strerror(iret)
        itmp = this%delete()
        return
      end if

      iret = nf90_get_var(ncid,varid,this%sfgrd)
      if ( iret /= nf90_noerr ) then
        write(error_unit,*) 'In file ',__FILE__,' at line ',__LINE__, &
          ' : Cannot read variable ',trim(sfgrd),' from file ', trim(datafile)
        write(error_unit,*) nf90_strerror(iret)
        itmp = this%delete()
        return
      end if

      iret = nf90_inq_varid(ncid,emrf,varid)
      if ( iret /= nf90_noerr ) then
        write(error_unit,*) 'In file ',__FILE__,' at line ',__LINE__, &
          ' : No variable ',trim(emrf), ' in file ', trim(datafile)
        write(error_unit,*) nf90_strerror(iret)
        itmp = this%delete()
        return
      end if

      iret = nf90_get_var(ncid,varid,this%emrf)
      if ( iret /= nf90_noerr ) then
        write(error_unit,*) 'In file ',__FILE__,' at line ',__LINE__, &
          ' : Cannot read variable ',trim(emrf),' from file ', trim(datafile)
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

      this%isinit = .true.
      iret = 0
    end function init_file

    integer function init_data(this,sfgrd,emrf) result(iret)
      implicit none
      class(aemissivity) :: this
      real , dimension(:) , intent(in) :: sfgrd
      real , dimension(:) , intent(in) :: emrf

      if ( this%isinit ) then
        iret = this%delete()
      end if

      if ( size(sfgrd) /= size(emrf) ) then
        write(error_unit,*) 'In file ',__FILE__,' at line ',__LINE__, &
          ' : Dimension mismatch'
        iret = -1
        return
      end if

      this%nf = size(sfgrd)
      allocate(this%sfgrd(this%nf), this%emrf(this%nf), stat=iret)
      if ( iret /= 0 ) then
        this%nf = 0
        write(error_unit,*) 'In file ',__FILE__,' at line ',__LINE__, &
          ' : Memory error allocating'
        return
      end if

      this%sfgrd = sfgrd
      this%emrf = emrf
      this%isinit = .true.
      iret = 0
    end function init_data

    integer function init_file_mpi(this,datafile,iocpu,icomm) result(iret)
      implicit none
      class(aemissivity) :: this
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
        allocate(this%sfgrd(this%nf), this%emrf(this%nf), stat=iret)
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

      call mpi_bcast(this%sfgrd,this%nf,mpi_real,iocpu,icomm,mpierr)
      if ( mpierr /= 0 ) then
        iret = this%delete( )
        iret = mpierr
      end if

      call mpi_bcast(this%emrf,this%nf,mpi_real,iocpu,icomm,mpierr)
      if ( mpierr /= 0 ) then
        iret = this%delete( )
        iret = mpierr
      end if

      if ( irank /= iocpu ) then
        this%isinit = .true.
      end if

      iret = 0
    end function init_file_mpi

    integer function delete(this) result(iret)
      implicit none
      class(aemissivity) :: this
      if ( allocated(this%emrf) ) then
        deallocate(this%sfgrd)
        deallocate(this%emrf)
        this%nf = -1
      end if
      this%isinit = .false.
      iret = 0
    end function delete

    subroutine cleanup(as)
      implicit none
      type(aemissivity) :: as
      if ( allocated(as%emrf) ) then
        deallocate(as%sfgrd)
        deallocate(as%emrf)
      end if
    end subroutine cleanup

end module mod_emissivity

#ifdef TESTME

! mpif90 -DTESTME -O3 -g -mtune=native -o test_emissivity \
!        -I`nf-config --includedir` mod_emissivity.F90 `nf-config --flibs`

program test_emissivity
  use mod_emissivity
  use mpi
  implicit none

  integer :: i , j , ni , icomm , myid , nproc
  real , dimension(1550) :: f
  type(aemissivity) :: as

  call mpi_init(i)
  if ( i /= 0 ) then
    print *, 'Cannot test mpi part'
  end if

  icomm = MPI_COMM_WORLD

  call mpi_comm_rank(icomm,myid,i)
  call mpi_comm_size(icomm,nproc,i)

  if ( myid == 0 ) then
    i = as%init_file('../data/emissivity.nc')
    print *, as%emrf(as%nf/2)
    i = as%delete()
  end if

  i = as%init_file_mpi('../data/emissivity.nc',0,icomm)

  print *, 'MYID : ', myid , ', NF = ',as%nf , 'EMSF = ',as%emrf(as%nf/2)

  i = as%delete()

  call mpi_finalize(i)

end program test_emissivity
#endif
