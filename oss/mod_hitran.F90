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
module mod_hitran

  use netcdf
  use mpi
  use iso_fortran_env

  implicit none

  private

  integer , parameter :: namelen = 16

  real , parameter :: odfac = 0.0E0

  type ahitran
    integer :: nchan
    integer :: nf
    integer :: nchmax
    integer :: nmol
    integer :: nlev
    integer :: ntmpod
    integer :: nlayod
    integer :: nmoltab
    integer :: mxmols
    real , dimension(:) , allocatable :: cwvn
    real , dimension(:) , allocatable :: pref
    real , dimension(:,:) , allocatable :: tmptab
    real , dimension(:,:) , allocatable :: wvptab
    real , dimension(:,:) , allocatable :: coef
    real , dimension(:) , allocatable :: vwvn
    real , dimension(:,:,:) , allocatable :: kfix
    real , dimension(:,:,:) , allocatable :: dkh2o
    real , dimension(:,:,:) , allocatable :: kh2o
    real , dimension(:,:,:,:) , allocatable :: kvar
    integer , dimension(:,:) , allocatable :: ichmap
    integer , dimension(:) , allocatable :: nmols
    integer(2) , dimension(:) , allocatable :: molid
    integer(2) , dimension(:,:) , allocatable :: imols
    integer(2) , dimension(:) , allocatable :: nch
    logical :: isinit = .false.
    contains
      procedure , pass :: init_file
      procedure , pass :: init_file_mpi
      procedure , pass :: init_data
      procedure , pass :: delete
      ! final :: cleanup
  end type ahitran

  public :: ahitran

  character(len=namelen) , parameter :: chandim = 'nchan'
  character(len=namelen) , parameter :: fdim = 'nf'
  character(len=namelen) , parameter :: chmaxdim = 'nchmax'
  character(len=namelen) , parameter :: moldim = 'nmol'
  character(len=namelen) , parameter :: levdim = 'nlev'
  character(len=namelen) , parameter :: tmpoddim = 'ntmpod'
  character(len=namelen) , parameter :: layoddim = 'nlayod'
  character(len=namelen) , parameter :: moltabdim = 'nmoltab'
  character(len=namelen) , parameter :: mxmolsdim = 'mxmols'

  character(len=namelen) , parameter :: pref = 'pref'
  character(len=namelen) , parameter :: tmptab = 'tmptab'
  character(len=namelen) , parameter :: cwvn = 'cWvn'
  character(len=namelen) , parameter :: vwvn = 'vwvn'
  character(len=namelen) , parameter :: wvptab = 'wvptab'
  character(len=namelen) , parameter :: kfix = 'kfix'
  character(len=namelen) , parameter :: kh2o = 'kh2o'
  character(len=namelen) , parameter :: dkh2o = 'dkh2o'
  character(len=namelen) , parameter :: kvar = 'kvar'
  character(len=namelen) , parameter :: imols = 'imols'
  character(len=namelen) , parameter :: molid = 'molid'
  character(len=namelen) , parameter :: nch = 'nch'
  character(len=namelen) , parameter :: coef = 'coef'
  character(len=namelen) , parameter :: ichmap = 'ichmap'

  contains

    integer function init_file(this,datafile,imolind) result(iret)
      implicit none
      class(ahitran) :: this
      character(len=*) , intent(in) :: datafile
      integer , dimension(:) , intent(in) , optional :: imolind
      integer :: ncid , dimid , varid , itmp
      integer , dimension(3) :: default_imolind = (/1,2,3/)
      integer :: n , k

      if ( this%isinit ) then
        iret = this%delete()
      end if

      iret = nf90_open(datafile,nf90_nowrite,ncid)
      if ( iret /= nf90_noerr ) then
        write(error_unit,*) 'In file ',__FILE__,' at line ',__LINE__, &
          ': Cannot open file ',trim(datafile)
        write(error_unit,*) nf90_strerror(iret)
        return
      end if

      iret = get_dim_len(ncid,datafile,chandim,this%nchan)
      if ( iret /= 0 ) then
        call zerodims(this)
        return
      end if
      iret = get_dim_len(ncid,datafile,fdim,this%nf)
      if ( iret /= 0 ) then
        call zerodims(this)
        return
      end if
      iret = get_dim_len(ncid,datafile,chmaxdim,this%nchmax)
      if ( iret /= 0 ) then
        call zerodims(this)
        return
      end if
      iret = get_dim_len(ncid,datafile,moldim,this%nmol)
      if ( iret /= 0 ) then
        call zerodims(this)
        return
      end if
      iret = get_dim_len(ncid,datafile,levdim,this%nlev)
      if ( iret /= 0 ) then
        call zerodims(this)
        return
      end if
      iret = get_dim_len(ncid,datafile,tmpoddim,this%ntmpod)
      if ( iret /= 0 ) then
        call zerodims(this)
        return
      end if
      iret = get_dim_len(ncid,datafile,layoddim,this%nlayod)
      if ( iret /= 0 ) then
        call zerodims(this)
        return
      end if
      iret = get_dim_len(ncid,datafile,moltabdim,this%nmoltab)
      if ( iret /= 0 ) then
        call zerodims(this)
        return
      end if
      iret = get_dim_len(ncid,datafile,mxmolsdim,this%mxmols)
      if ( iret /= 0 ) then
        call zerodims(this)
        return
      end if

      allocate(this%cwvn(this%nchan), this%pref(this%nlev),  &
               this%tmptab(this%ntmpod,this%nlayod),  &
               this%wvptab(this%nmoltab,this%nlayod),  &
               this%coef(this%nchmax,this%nf), this%vwvn(this%nf),  &
               this%kfix(this%nlayod,this%ntmpod,this%nf),  &
               this%dkh2o(this%nlayod,this%ntmpod,this%nf),  &
               this%kh2o(this%nlayod,this%ntmpod,this%nf),  &
               this%kvar(this%mxmols,this%nlayod,this%ntmpod,this%nf),  &
               this%ichmap(this%nchmax,this%nf), &
               this%molid(this%nmol), &
               this%imols(this%mxmols,this%nf), &
               this%nch(this%nf), stat=iret)
      if ( iret /= 0 ) then
        call zerodims(this)
        write(error_unit,*) 'In file ',__FILE__,' at line ',__LINE__, &
          ' : Memory error allocating'
        return
      end if

      iret = nf90_inq_varid(ncid,cwvn,varid)
      if ( iret /= nf90_noerr ) then
        write(error_unit,*) 'In file ',__FILE__,' at line ',__LINE__, &
          ' : No variable ',trim(cwvn), ' in file ', trim(datafile)
        write(error_unit,*) nf90_strerror(iret)
        itmp = this%delete()
        return
      end if
      iret = nf90_get_var(ncid,varid,this%cwvn)
      if ( iret /= nf90_noerr ) then
        write(error_unit,*) 'In file ',__FILE__,' at line ',__LINE__, &
          ' : Cannot read variable ',trim(cwvn),' from file ', trim(datafile)
        write(error_unit,*) nf90_strerror(iret)
        itmp = this%delete()
        return
      end if
      iret = nf90_inq_varid(ncid,molid,varid)
      if ( iret /= nf90_noerr ) then
        write(error_unit,*) 'In file ',__FILE__,' at line ',__LINE__, &
          ' : No variable ',trim(molid),' in file ', trim(datafile)
        write(error_unit,*) nf90_strerror(iret)
        itmp = this%delete()
        return
      end if
      iret = nf90_get_var(ncid,varid,this%molid)
      if ( iret /= nf90_noerr ) then
        write(error_unit,*) 'In file ',__FILE__,' at line ',__LINE__, &
          ' : Cannot read variable ',trim(molid),' from file ', trim(datafile)
        write(error_unit,*) nf90_strerror(iret)
        itmp = this%delete()
        return
      end if
      iret = nf90_inq_varid(ncid,pref,varid)
      if ( iret /= nf90_noerr ) then
        write(error_unit,*) 'In file ',__FILE__,' at line ',__LINE__, &
          ' : No variable ',trim(pref), ' in file ', trim(datafile)
        write(error_unit,*) nf90_strerror(iret)
        itmp = this%delete()
        return
      end if
      iret = nf90_get_var(ncid,varid,this%pref)
      if ( iret /= nf90_noerr ) then
        write(error_unit,*) 'In file ',__FILE__,' at line ',__LINE__, &
          ' : Cannot read variable ',trim(pref), ' from file ', trim(datafile)
        write(error_unit,*) nf90_strerror(iret)
        itmp = this%delete()
        return
      end if
      iret = nf90_inq_varid(ncid,tmptab,varid)
      if ( iret /= nf90_noerr ) then
        write(error_unit,*) 'In file ',__FILE__,' at line ',__LINE__, &
          ' : No variable ',trim(tmptab), ' in file ', trim(datafile)
        write(error_unit,*) nf90_strerror(iret)
        itmp = this%delete()
        return
      end if
      iret = nf90_get_var(ncid,varid,this%tmptab)
      if ( iret /= nf90_noerr ) then
        write(error_unit,*) 'In file ',__FILE__,' at line ',__LINE__, &
          ' : Cannot read variable ',trim(tmptab),' from file ', trim(datafile)
        write(error_unit,*) nf90_strerror(iret)
        itmp = this%delete()
        return
      end if
      iret = nf90_inq_varid(ncid,wvptab,varid)
      if ( iret /= nf90_noerr ) then
        write(error_unit,*) 'In file ',__FILE__,' at line ',__LINE__, &
          ' : No variable ',trim(wvptab), ' in file ', trim(datafile)
        write(error_unit,*) nf90_strerror(iret)
        itmp = this%delete()
        return
      end if
      iret = nf90_get_var(ncid,varid,this%wvptab)
      if ( iret /= nf90_noerr ) then
        write(error_unit,*) 'In file ',__FILE__,' at line ',__LINE__, &
         ' : Cannot read variable ',trim(wvptab),' from file ', trim(datafile)
        write(error_unit,*) nf90_strerror(iret)
        itmp = this%delete()
        return
      end if
      iret = nf90_inq_varid(ncid,coef,varid)
      if ( iret /= nf90_noerr ) then
        write(error_unit,*) 'In file ',__FILE__,' at line ',__LINE__, &
         ' : No variable ',trim(coef),' in file ', trim(datafile)
        write(error_unit,*) nf90_strerror(iret)
        itmp = this%delete()
        return
      end if
      iret = nf90_get_var(ncid,varid,this%coef)
      if ( iret /= nf90_noerr ) then
        write(error_unit,*) 'In file ',__FILE__,' at line ',__LINE__, &
         ' : Cannot read variable ',trim(coef),' from file ', trim(datafile)
        write(error_unit,*) nf90_strerror(iret)
        itmp = this%delete()
        return
      end if
      iret = nf90_inq_varid(ncid,ichmap,varid)
      if ( iret /= nf90_noerr ) then
        write(error_unit,*) 'In file ',__FILE__,' at line ',__LINE__, &
          ' : No variable ',trim(ichmap), ' in file ', trim(datafile)
        write(error_unit,*) nf90_strerror(iret)
        itmp = this%delete()
        return
      end if
      iret = nf90_get_var(ncid,varid,this%ichmap)
      if ( iret /= nf90_noerr ) then
        write(error_unit,*) 'In file ',__FILE__,' at line ',__LINE__, &
          ' : Cannot read variable ',trim(ichmap),' from file ', trim(datafile)
        write(error_unit,*) nf90_strerror(iret)
        itmp = this%delete()
        return
      end if
      iret = nf90_inq_varid(ncid,nch,varid)
      if ( iret /= nf90_noerr ) then
        write(error_unit,*) 'In file ',__FILE__,' at line ',__LINE__, &
          ' : No variable ',trim(nch), ' in file ', trim(datafile)
        write(error_unit,*) nf90_strerror(iret)
        itmp = this%delete()
        return
      end if
      iret = nf90_get_var(ncid,varid,this%nch)
      if ( iret /= nf90_noerr ) then
        write(error_unit,*) 'In file ',__FILE__,' at line ',__LINE__, &
          ' : Cannot read variable ',trim(nch), ' from file ', trim(datafile)
        write(error_unit,*) nf90_strerror(iret)
        itmp = this%delete()
        return
      end if
      iret = nf90_inq_varid(ncid,vwvn,varid)
      if ( iret /= nf90_noerr ) then
        write(error_unit,*) 'In file ',__FILE__,' at line ',__LINE__, &
          ' : No variable ',trim(vwvn), ' in file ', trim(datafile)
        write(error_unit,*) nf90_strerror(iret)
        itmp = this%delete()
        return
      end if
      iret = nf90_get_var(ncid,varid,this%vwvn)
      if ( iret /= nf90_noerr ) then
        write(error_unit,*) 'In file ',__FILE__,' at line ',__LINE__, &
          ' : Cannot read variable ',trim(vwvn),' from file ', trim(datafile)
        write(error_unit,*) nf90_strerror(iret)
        itmp = this%delete()
        return
      end if
      iret = nf90_inq_varid(ncid,imols,varid)
      if ( iret /= nf90_noerr ) then
        write(error_unit,*) 'In file ',__FILE__,' at line ',__LINE__, &
          ' : No variable ',trim(imols), ' in file ', trim(datafile)
        write(error_unit,*) nf90_strerror(iret)
        itmp = this%delete()
        return
      end if
      iret = nf90_get_var(ncid,varid,this%imols)
      if ( iret /= nf90_noerr ) then
        write(error_unit,*) 'In file ',__FILE__,' at line ',__LINE__, &
          ' : Cannot read variable ',trim(imols),' from file ', trim(datafile)
        write(error_unit,*) nf90_strerror(iret)
        itmp = this%delete()
        return
      end if
      iret = nf90_inq_varid(ncid,kfix,varid)
      if ( iret /= nf90_noerr ) then
        write(error_unit,*) 'In file ',__FILE__,' at line ',__LINE__, &
          ' : No variable ',trim(kfix), ' in file ', trim(datafile)
        write(error_unit,*) nf90_strerror(iret)
        itmp = this%delete()
        return
      end if
      iret = nf90_get_var(ncid,varid,this%kfix)
      if ( iret /= nf90_noerr ) then
        write(error_unit,*) 'In file ',__FILE__,' at line ',__LINE__, &
          ' : Cannot read variable ',trim(kfix), ' from file ', trim(datafile)
        write(error_unit,*) nf90_strerror(iret)
        itmp = this%delete()
        return
      end if
      iret = nf90_inq_varid(ncid,dkh2o,varid)
      if ( iret /= nf90_noerr ) then
        write(error_unit,*) 'In file ',__FILE__,' at line ',__LINE__, &
          ' : No variable ',trim(dkh2o), ' in file ', trim(datafile)
        write(error_unit,*) nf90_strerror(iret)
        itmp = this%delete()
        return
      end if
      iret = nf90_get_var(ncid,varid,this%dkh2o)
      if ( iret /= nf90_noerr ) then
        write(error_unit,*) 'In file ',__FILE__,' at line ',__LINE__, &
          ' : Cannot read variable ',trim(dkh2o),' from file ', trim(datafile)
        write(error_unit,*) nf90_strerror(iret)
        itmp = this%delete()
        return
      end if
      iret = nf90_inq_varid(ncid,kh2o,varid)
      if ( iret /= nf90_noerr ) then
        write(error_unit,*) 'In file ',__FILE__,' at line ',__LINE__, &
          ' : No variable ',trim(kh2o), ' in file ', trim(datafile)
        write(error_unit,*) nf90_strerror(iret)
        itmp = this%delete()
        return
      end if
      iret = nf90_get_var(ncid,varid,this%kh2o)
      if ( iret /= nf90_noerr ) then
        write(error_unit,*) 'In file ',__FILE__,' at line ',__LINE__, &
          ' : Cannot read variable ',trim(kh2o), ' from file ', trim(datafile)
        write(error_unit,*) nf90_strerror(iret)
        itmp = this%delete()
        return
      end if
      iret = nf90_inq_varid(ncid,kvar,varid)
      if ( iret /= nf90_noerr ) then
        write(error_unit,*) 'In file ',__FILE__,' at line ',__LINE__, &
          ' : No variable ',trim(kvar), ' in file ', trim(datafile)
        write(error_unit,*) nf90_strerror(iret)
        itmp = this%delete()
        return
      end if
      iret = nf90_get_var(ncid,varid,this%kvar)
      if ( iret /= nf90_noerr ) then
        write(error_unit,*) 'In file ',__FILE__,' at line ',__LINE__, &
          ' : Cannot read variable ',trim(kvar), ' from file ', trim(datafile)
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

      do n = 1 , this%nf
        do k = 1 , this%ntmpod
          this%kfix(:,k,n) = (this%kfix(:,k,n) * this%wvptab(1,:))
        end do
      end do

      if ( present(imolind) ) then
        iret = select_species(this,imolind)
      else
        iret = select_species(this,default_imolind)
      end if
      if ( iret /= 0 ) then
        itmp = this%delete()
        return
      end if

      this%isinit = .true.
      iret = 0
    end function init_file

    integer function init_data(this,cwvn,pref,tmptab,wvptab,coef,vwvn, &
              kfix,dkh2o,kh2o,kvar,ichmap,molid,imols,nch) result(iret)
      implicit none
      class(ahitran) :: this
      real , dimension(:) :: cwvn
      real , dimension(:) :: pref
      real , dimension(:,:) :: tmptab
      real , dimension(:,:) :: wvptab
      real , dimension(:,:) :: coef
      real , dimension(:) :: vwvn
      real , dimension(:,:,:) :: kfix
      real , dimension(:,:,:) :: dkh2o
      real , dimension(:,:,:) :: kh2o
      real , dimension(:,:,:,:) :: kvar
      integer , dimension(:,:) :: ichmap
      integer(2) , dimension(:) :: molid
      integer(2) , dimension(:,:) :: imols
      integer(2) , dimension(:) :: nch

      if ( this%isinit ) then
        iret = this%delete()
      end if

      if ( size(wvptab,2) /= size(tmptab,2) .or. &
           size(vwvn) /= size(coef,2) .or.       &
           size(vwvn) /= size(kfix,3) .or.       &
           size(vwvn) /= size(dkh2o,3) .or.      &
           size(vwvn) /= size(kh2o,3) .or.       &
           size(vwvn) /= size(kvar,4) .or.       &
           size(vwvn) /= size(ichmap,2) .or.     &
           size(vwvn) /= size(imols,2) .or.      &
           size(vwvn) /= size(nch) .or.          &
           size(kfix,2) /= size(dkh2o,2) .or.    &
           size(kfix,2) /= size(kh2o,2) .or.     &
           size(kfix,2) /= size(kvar,2) .or.     &
           size(tmptab,2) /= size(kfix,1) .or.   &
           size(tmptab,2) /= size(dkh2o,1) .or.  &
           size(tmptab,2) /= size(kh2o,1) .or.   &
           size(tmptab,2) /= size(kvar,1) .or.   &
           size(imols,2) /= size(kvar,2) .or.    &
           size(ichmap,1) /= size(coef,1) ) then
        write(error_unit,*) 'In file ',__FILE__,' at line ',__LINE__, &
          ' : Dimension mismatch'
        iret = -1
        return
      end if

      this%nchan = size(cwvn)
      this%nf = size(vwvn)
      this%nchmax = size(coef,1)
      this%nmol = size(molid)
      this%nlev = size(pref)
      this%ntmpod = size(tmptab,1)
      this%nlayod = size(tmptab,2)
      this%nmoltab = size(wvptab,1)
      this%mxmols = size(imols,1)

      allocate(this%cwvn(this%nchan), this%pref(this%nlev),  &
               this%tmptab(this%ntmpod,this%nlayod),  &
               this%wvptab(this%nmoltab,this%nlayod),  &
               this%coef(this%nchmax,this%nf), this%vwvn(this%nf),  &
               this%kfix(this%nlayod,this%ntmpod,this%nf),  &
               this%dkh2o(this%nlayod,this%ntmpod,this%nf),  &
               this%kh2o(this%nlayod,this%ntmpod,this%nf),  &
               this%kvar(this%mxmols,this%nlayod,this%ntmpod,this%nf),  &
               this%ichmap(this%nchmax,this%nf), &
               this%molid(this%nmol), &
               this%imols(this%mxmols,this%nf), &
               this%nch(this%nf), stat=iret)
      if ( iret /= 0 ) then
        call zerodims(this)
        write(error_unit,*) 'In file ',__FILE__,' at line ',__LINE__, &
          ' : Memory error allocating'
        return
      end if

      this%cwvn = cwvn
      this%pref = pref
      this%tmptab = tmptab
      this%wvptab = wvptab
      this%coef = coef
      this%vwvn = vwvn
      this%kfix = kfix
      this%dkh2o = dkh2o
      this%kh2o = kh2o
      this%kvar = kvar
      this%ichmap = ichmap
      this%molid = molid
      this%imols = imols
      this%nch = nch
      this%isinit = .true.
      iret = 0
    end function init_data

    integer function init_file_mpi(this,datafile, &
                                   iocpu,icomm,imolind) result(iret)
      implicit none
      class(ahitran) :: this
      integer , dimension(:) , intent(in) , optional :: imolind
      character(len=*) , intent(in) :: datafile
      integer , intent(in) :: iocpu , icomm
      integer , dimension(9) :: dimarr
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
        if ( present(imolind) ) then
          iret = init_file(this,datafile,imolind)
        else
          iret = init_file(this,datafile)
        end if
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
        dimarr(1) = this%nchan
        dimarr(2) = this%nf
        dimarr(3) = this%nchmax
        dimarr(4) = this%nmol
        dimarr(5) = this%nlev
        dimarr(6) = this%ntmpod
        dimarr(7) = this%nlayod
        dimarr(8) = this%nmoltab
        dimarr(9) = this%mxmols
      end if

      call mpi_bcast(dimarr,9,mpi_integer,iocpu,icomm,mpierr)
      if ( mpierr /= 0 ) then
        if ( irank == iocpu ) iret = this%delete( )
        iret = mpierr
        return
      end if

      if ( irank /= iocpu ) then
        this%nchan = dimarr(1)
        this%nf = dimarr(2)
        this%nchmax = dimarr(3)
        this%nmol = dimarr(4)
        this%nlev = dimarr(5)
        this%ntmpod = dimarr(6)
        this%nlayod = dimarr(7)
        this%nmoltab = dimarr(8)
        this%mxmols = dimarr(9)
        allocate(this%cwvn(this%nchan), this%pref(this%nlev),  &
                 this%tmptab(this%ntmpod,this%nlayod),  &
                 this%wvptab(this%nmoltab,this%nlayod),  &
                 this%coef(this%nchmax,this%nf), this%vwvn(this%nf),  &
                 this%kfix(this%nlayod,this%ntmpod,this%nf),  &
                 this%dkh2o(this%nlayod,this%ntmpod,this%nf),  &
                 this%kh2o(this%nlayod,this%ntmpod,this%nf),  &
                 this%kvar(this%nmol,this%nlayod,this%ntmpod,this%nf),  &
                 this%ichmap(this%nchmax,this%nf), &
                 this%molid(this%nmol), &
                 this%imols(this%nmol,this%nf), &
                 this%nch(this%nf), stat=iret)
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

      call mpi_bcast(this%cwvn,this%nchan,mpi_real,iocpu,icomm,mpierr)
      if ( mpierr /= 0 ) then
        iret = this%delete( )
        iret = mpierr
      end if
      call mpi_bcast(this%pref,this%nlev,mpi_real,iocpu,icomm,mpierr)
      if ( mpierr /= 0 ) then
        iret = this%delete( )
        iret = mpierr
      end if
      call mpi_bcast(this%tmptab,product(shape(this%tmptab)), &
                     mpi_real,iocpu,icomm,mpierr)
      if ( mpierr /= 0 ) then
        iret = this%delete( )
        iret = mpierr
      end if
      call mpi_bcast(this%wvptab,product(shape(this%wvptab)), &
                     mpi_real,iocpu,icomm,mpierr)
      if ( mpierr /= 0 ) then
        iret = this%delete( )
        iret = mpierr
      end if
      call mpi_bcast(this%coef,product(shape(this%coef)), &
                     mpi_real,iocpu,icomm,mpierr)
      if ( mpierr /= 0 ) then
        iret = this%delete( )
        iret = mpierr
      end if
      call mpi_bcast(this%vwvn,this%nf,mpi_real,iocpu,icomm,mpierr)
      if ( mpierr /= 0 ) then
        iret = this%delete( )
        iret = mpierr
      end if
      call mpi_bcast(this%kfix,product(shape(this%kfix)), &
                     mpi_real,iocpu,icomm,mpierr)
      if ( mpierr /= 0 ) then
        iret = this%delete( )
        iret = mpierr
      end if
      call mpi_bcast(this%dkh2o,product(shape(this%dkh2o)), &
                     mpi_real,iocpu,icomm,mpierr)
      if ( mpierr /= 0 ) then
        iret = this%delete( )
        iret = mpierr
      end if
      call mpi_bcast(this%kh2o,product(shape(this%kh2o)), &
                     mpi_real,iocpu,icomm,mpierr)
      if ( mpierr /= 0 ) then
        iret = this%delete( )
        iret = mpierr
      end if
      call mpi_bcast(this%kvar,product(shape(this%kvar)), &
                     mpi_real,iocpu,icomm,mpierr)
      if ( mpierr /= 0 ) then
        iret = this%delete( )
        iret = mpierr
      end if
      call mpi_bcast(this%ichmap,product(shape(this%ichmap)), &
                     mpi_integer,iocpu,icomm,mpierr)
      if ( mpierr /= 0 ) then
        iret = this%delete( )
        iret = mpierr
      end if
      call mpi_bcast(this%molid,this%nmol,mpi_integer2,iocpu,icomm,mpierr)
      if ( mpierr /= 0 ) then
        iret = this%delete( )
        iret = mpierr
      end if
      call mpi_bcast(this%imols,product(shape(this%imols)), &
                     mpi_integer2,iocpu,icomm,mpierr)
      if ( mpierr /= 0 ) then
        iret = this%delete( )
        iret = mpierr
      end if
      call mpi_bcast(this%nch,this%nf,mpi_integer2,iocpu,icomm,mpierr)
      if ( mpierr /= 0 ) then
        iret = this%delete( )
        iret = mpierr
      end if

      if ( irank /= iocpu ) then
        this%isinit = .true.
      end if

      iret = 0
    end function init_file_mpi

    integer function select_species(this,molid) result(iret)
      implicit none
      type(ahitran) :: this
      integer , dimension(:) :: molid
      integer(2) , dimension(:) , allocatable :: new_molid
      integer(2) , dimension(:,:) , allocatable :: new_imols
      real , dimension(:,:,:,:) , allocatable :: new_kvar
      integer , dimension(this%mxmols) :: map
      integer , dimension(this%mxmols) :: maps
      integer , dimension(this%mxmols) :: iflag
      integer , dimension(this%mxmols) :: imols_indx
      real , dimension(this%mxmols) :: kbuf
      real :: ktot
      integer :: nmol , nm , iimol
      integer :: ismp , hinmol , nt , n , l , k , ks , kk , imol

      ! Check molecular selection requested by user.

      if ( count(molid > 0) > this%mxmols ) then
        write(error_unit,*) 'In file ',__FILE__,' at line ',__LINE__, &
          ' Database contains maximum ',this%mxmols,' species !'
        iret = -1
        return
      end if
      if ( molid(1) /= 1 ) then
        write(error_unit,*) 'In file ',__FILE__,' at line ',__LINE__, &
          ' : Water vapor must be specified in species selection.'
        iret = -1
        return
      end if
      nmol = count(molid > 0)
      nm = nmol
      do n = 2 , size(molid)
        if ( molid(n) == 0 ) then
          nm = n - 1
          exit
        end if
        if ( molid(n) <= molid(n-1) ) then
          write(error_unit,*) 'In file ',__FILE__,' at line ',__LINE__, &
            ' : Molecular selection of HITRAN ids must be in ascending order'
          iret = -1
          return
        end if
      end do
      if ( nm /= nmol ) then
        write(error_unit,*) 'In file ',__FILE__,' at line ',__LINE__, &
          ' : Spurious data in array molid : ', molid
        iret = -1
        return
      end if

      allocate(new_molid(nmol),new_imols(nmol,this%nf), &
               new_kvar(nmol,this%nlayod,this%ntmpod,this%nf),  &
               this%nmols(this%nf), stat=iret)
      if ( iret /= 0 ) then
        write(error_unit,*) 'In file ',__FILE__,' at line ',__LINE__, &
          ' : Allocation error !'
        iret = -1
        return
      end if

      new_molid = molid(1:nmol)
      ! Check and map the requested species
      map(:) = 0
      iflag(:) = 1
      kk = 0
      iimol = count(this%molid > 0)
      do k = 1 , iimol
        if ( any(new_molid(:) == this%molid(k)) ) then
          kk = kk + 1
          map(k) = kk
        end if
      end do
      if ( kk /= nmol ) then
        write(error_unit,*) 'Database species: ', this%molid(:)
        write(error_unit,*) 'Selected species: ', new_molid(:)
        write(error_unit,*) 'In file ',__FILE__,' at line ',__LINE__, &
          ' : Requested molecules not in the database file.'
        iret = -1
        return
      end if

      do ismp = 1 , this%nf
        hinmol = count(this%imols(:,ismp) > 0)
        iflag(2:iimol) = 0
        do nt = 1 , this%ntmpod
          do l = 1 , this%nlayod
            ktot = (this%kfix(l,nt,ismp))
            do ks = 2 , hinmol
              kbuf(ks-1) = this%kvar(ks-1,l,nt,ismp) * &
                           this%wvptab(this%imols(ks,ismp)+2,l)
              ktot = ktot + kbuf(ks-1)
            end do
            do ks = 2 , hinmol
              if ( kbuf(ks-1) > (odfac * ktot) ) then
                iflag(this%imols(ks,ismp)) = 1
              end if
            end do
          end do
        end do
        maps(1:iimol) = map(1:iimol) * iflag(1:iimol)
        kk = 0
        do ks = 1 , hinmol
          imol = this%imols(ks,ismp)
          if ( maps(imol) > 0 ) then
            kk = kk+1
            imols_indx(kk) = ks
          else
            do nt = 1 , this%ntmpod
              do l = 1 , this%nlayod
                this%kfix(l,nt,ismp) = this%kfix(l,nt,ismp) + &
                  (this%kvar(ks-1,l,nt,ismp) * this%wvptab(imol+2,l))
              end do
            end do
          end if
        end do
        this%nmols(ismp) = kk
        new_imols(1:kk,ismp) = maps(this%imols(imols_indx(1:kk),ismp))
        new_kvar(1:kk-1,:,:,ismp) = (this%kvar(imols_indx(2:kk)-1,:,:,ismp))
      end do

      this%nmol = nmol

      call move_alloc(new_kvar,this%kvar)
      call move_alloc(new_molid,this%molid)
      call move_alloc(new_imols,this%imols)

      iret = 0
    end function select_species

    integer function delete(this) result(iret)
      implicit none
      class(ahitran) :: this
      call cleanup(this)
      call zerodims(this)
      iret = 0
    end function delete

    subroutine cleanup(ah)
      implicit none
      type(ahitran) :: ah
      if ( allocated(ah%cwvn) ) then
        deallocate(ah%cwvn)
        deallocate(ah%pref)
        deallocate(ah%tmptab)
        deallocate(ah%wvptab)
        deallocate(ah%coef)
        deallocate(ah%vwvn)
        deallocate(ah%kfix)
        deallocate(ah%dkh2o)
        deallocate(ah%kh2o)
        deallocate(ah%kvar)
        deallocate(ah%ichmap)
        deallocate(ah%molid)
        deallocate(ah%imols)
        deallocate(ah%nch)
        if ( allocated(ah%nmols) ) deallocate(ah%nmols)
      end if
    end subroutine cleanup

    subroutine zerodims(ah)
      implicit none
      type(ahitran) :: ah
      ah%nchan = 0
      ah%nf = 0
      ah%nchmax = 0
      ah%nmol = 0
      ah%nlev = 0
      ah%ntmpod = 0
      ah%nlayod = 0
      ah%nmoltab = 0
      ah%mxmols = 0
    end subroutine zerodims

    integer function get_dim_len(ncid,datafile,dname,dlen) result(iret)
      implicit none
      integer , intent(in) :: ncid
      character(len=*) , intent(in) :: datafile
      character(len=*) , intent(in) :: dname
      integer , intent(out) :: dlen
      integer :: dimid
      iret = nf90_inq_dimid(ncid,dname,dimid)
      if ( iret /= nf90_noerr ) then
        write(error_unit,*) 'In file ',__FILE__,' at line ',__LINE__, &
          ' : No dimension ',trim(dname),' in file ', trim(datafile)
        write(error_unit,*) nf90_strerror(iret)
        return
      end if
      iret = nf90_inquire_dimension(ncid,dimid,len=dlen)
      if ( iret /= nf90_noerr ) then
        write(error_unit,*) 'In file ',__FILE__,' at line ',__LINE__, &
          ' : Cannot read dimension ',trim(dname),' from file ', trim(datafile)
        write(error_unit,*) nf90_strerror(iret)
        return
      end if
    end function get_dim_len

end module mod_hitran

#ifdef TESTME

! mpif90 -DTESTME -O3 -g -mtune=native -o test_hitran \
!        -I`nf-config --includedir` mod_hitran.F90 `nf-config --flibs`

program test_hitran
  use mod_hitran
  use mpi
  implicit none

  integer :: i , j , ni , icomm , myid , nproc
  real , dimension(1550) :: f
  type(ahitran) :: as

  call mpi_init(i)
  if ( i /= 0 ) then
    print *, 'Cannot test mpi part'
  end if

  icomm = MPI_COMM_WORLD

  call mpi_comm_rank(icomm,myid,i)
  call mpi_comm_size(icomm,nproc,i)

  if ( myid == 0 ) then
    i = as%init_file('../data/leo.iasi.0.05.nc')
    if ( i /= 0 ) then
      print *, 'Error init !'
      stop
    end if
    print *, as%pref(as%nlev/2)
    i = as%delete()
  end if

  i = as%init_file_mpi('../data/leo.iasi.0.05.nc',0,icomm,(/1,2,3,19/))

  print *, 'MYID : ', myid , ', NF = ',as%nf , 'PREF = ',as%pref(as%nlev/2)

  i = as%delete()

  call mpi_finalize(i)

end program test_hitran
#endif
