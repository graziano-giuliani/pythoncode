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
    contains
      procedure , pass :: init
      procedure , pass :: delete
      procedure , pass :: interpolate
  end type asolar

  public :: asolar

  character(len=namelen) , parameter :: fdim = 'nspectral'
  character(len=namelen) , parameter :: frq = 'FREQ'
  character(len=namelen) , parameter :: irr = 'IRRADIANCE'

  contains

    integer function init(this,datafile) result(iret)
      implicit none
      class(asolar) :: this
      character(len=*) , intent(in) :: datafile
      integer :: ncid , dimid , varid

      iret = nf90_open(datafile,nf90_nowrite,ncid)
      if ( iret /= nf90_noerr ) then
        write(error_unit,*) 'Cannot open file ',trim(datafile)
        write(error_unit,*) nf90_strerror(iret)
        return
      end if

      iret = nf90_inq_dimid(ncid,fdim,dimid)
      if ( iret /= nf90_noerr ) then
        write(error_unit,*) 'No dimension ',trim(fdim),' in file ', &
             trim(datafile)
        write(error_unit,*) nf90_strerror(iret)
        return
      end if

      iret = nf90_inquire_dimension(ncid,dimid,len=this%nf)
      if ( iret /= nf90_noerr ) then
        write(error_unit,*) 'Cannot read dimension ',trim(fdim), &
             ' from file ', trim(datafile)
        write(error_unit,*) nf90_strerror(iret)
        return
      end if

      allocate(this%frq(this%nf), &
               this%irr(this%nf), stat=iret)
      if ( iret /= 0 ) then
        write(error_unit,*) 'Memory error allocating'
        return
      end if

      iret = nf90_inq_varid(ncid,frq,varid)
      if ( iret /= nf90_noerr ) then
        write(error_unit,*) 'No variable ',trim(frq), &
             ' in file ', trim(datafile)
        write(error_unit,*) nf90_strerror(iret)
        return
      end if

      iret = nf90_get_var(ncid,varid,this%frq)
      if ( iret /= nf90_noerr ) then
        write(error_unit,*) 'Cannot read variable ',trim(frq), &
             ' from file ', trim(datafile)
        write(error_unit,*) nf90_strerror(iret)
        return
      end if

      iret = nf90_inq_varid(ncid,irr,varid)
      if ( iret /= nf90_noerr ) then
        write(error_unit,*) 'No variable ',trim(irr), &
             ' in file ', trim(datafile)
        write(error_unit,*) nf90_strerror(iret)
        return
      end if

      iret = nf90_get_var(ncid,varid,this%irr)
      if ( iret /= nf90_noerr ) then
        write(error_unit,*) 'Cannot read variable ',trim(irr), &
             ' from file ', trim(datafile)
        write(error_unit,*) nf90_strerror(iret)
        return
      end if

      iret = nf90_close(ncid)
      if ( iret /= nf90_noerr ) then
        write(error_unit,*) 'Cannot close file ', trim(datafile)
        write(error_unit,*) nf90_strerror(iret)
        return
      end if

      this%f1 = this%frq(1)
      this%df = this%frq(2) - this%frq(1)
      iret = 0
    end function init

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
      end if
      iret = 0
    end function delete

    integer function interpolate(this,f) result(iret)
      implicit none
      class(asolar) :: this
      real , dimension(:) , intent(in) :: f
      real :: xf , w1 , w2
      integer :: i , i1 , i2 , nif
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

end module mod_solar

#ifdef TESTME

! gfortran -DTESTME -O3 -g -mtune=native -o test_solar \
!           -I`nf-config --includedir` mod_solar.F90 `nf-config --flibs`

program test_solar
  use mod_solar
  implicit none

  integer :: i , j , ni
  real , dimension(1550) :: f
  type(asolar) :: as

  i = as%init('../data/solar_irradiances.nc')

  do ni = 0 , 999
    do j = 1 , 1550
      f(j) = 1550.0+ni + real(j-1)/(1550.0/1000.0)
    end do

    i = as%interpolate(f)

    print *, f(1) , f(1550)
    print *, as%ivals(1) , as%ivals(1550)
  end do

  i = as%delete()

end program test_solar
#endif
