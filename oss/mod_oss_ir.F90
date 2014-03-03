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
  use mod_forward_model_data
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

  integer :: ipsfcg = -1
  integer :: itempg = -1
  integer :: itsking = -1

  logical , parameter :: lookup = .false.
  logical , parameter :: linintau = .false.

  ! 'variablegrid' specifies whether the input is on the absorption
  !  coefficient grid (i.e. normal operating mode) or on a user-specified grid.
  logical :: variablegrid = .false.
  logical :: is_planck_set = .false.

  type aoss_ir
    real , dimension(:) , allocatable :: pavlref
    real , dimension(:) , allocatable :: sunrad
    real , dimension(:) , allocatable :: ap1
    real , dimension(:) , allocatable :: ap2
    real , dimension(:) , allocatable :: at1_p1
    real , dimension(:) , allocatable :: at2_p1
    real , dimension(:) , allocatable :: at1_p2
    real , dimension(:) , allocatable :: at2_p2
    real , dimension(:) , allocatable :: adt1_p1
    real , dimension(:) , allocatable :: adt2_p1
    real , dimension(:) , allocatable :: adt1_p2
    real , dimension(:) , allocatable :: adt2_p2
    integer , dimension(:) , allocatable :: indxt_p1
    integer , dimension(:) , allocatable :: indxt_p2
    integer , dimension(:) , allocatable :: indxp
    real , dimension(:) , allocatable :: tavl
    real , dimension(:) , allocatable :: pavl
    real , dimension(:) , allocatable :: wfix
    real , dimension(:) , allocatable :: dtu
    real , dimension(:) , allocatable :: dtl
    real , dimension(:,:) , allocatable :: q
    real , dimension(:,:) , allocatable :: w
    real , dimension(:,:) , allocatable :: dwqu
    real , dimension(:,:) , allocatable :: dwql
    real , dimension(:) , allocatable :: qobs
    real , dimension(:) , allocatable :: wobs
    real , dimension(:) , allocatable :: dwquobs
    real , dimension(:) , allocatable :: dwqlobs
    real , dimension(:,:) , allocatable :: abso
    real , dimension(:) , allocatable :: tautot
    real , dimension(:) , allocatable :: dtaudtmp
    real , dimension(:) , allocatable :: absoobs
    real , dimension(:) , allocatable :: txdn
    real , dimension(:) , allocatable :: txup
    real , dimension(:) , allocatable :: bbar
    real , dimension(:) , allocatable :: dbbar
    real , dimension(:) , allocatable :: draddtmp
    real , dimension(:) , allocatable :: draddtau
    real , dimension(:) , allocatable :: drdw
    real , dimension(:) , allocatable :: draddtmpdw
    real , dimension(:) , allocatable :: draddtmpuw
    real , dimension(:) , allocatable :: dp
    real , dimension(:) , allocatable :: qcor
    real , dimension(:) , allocatable :: qr
    real , dimension(:) , allocatable :: ploc
    real , dimension(:) , allocatable :: b2
    real , dimension(:) , allocatable :: db2
    real , dimension(:) , allocatable :: ar
    real , dimension(:) , allocatable :: br
    real , dimension(:) , allocatable :: ad
    real , dimension(:) , allocatable :: bd
    real :: tautotobs
    real :: dtaudtmpobs
    real :: tavlobs
    real :: wfixobs
    real :: ap1obs
    real :: ap2obs
    real :: dtuobs
    real :: dtlobs
    real :: at1_p1obs
    real :: at2_p1obs
    real :: adt1_p1obs
    real :: adt2_p1obs
    real :: v2
    real :: bs2
    real :: dbs2
    real :: ars
    real :: brs
    real :: ads
    real :: bds
    integer :: indxt_p1obs = 0
    integer :: indxpobs
    type(asolar) , pointer :: as
    type(ahitran) , pointer :: ah
    integer :: mxlay
    integer :: mxhmol
    integer :: mxlev
    logical :: isinit = .false.
    contains
      procedure , pass :: init
      procedure , pass :: forward_model
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

      this%mxlev = ah%nlev
      this%mxlay = this%mxlev - 1
      this%mxhmol = ah%nmol

      allocate(this%ap1(this%mxlay), this%ap2(this%mxlay), &
               this%at1_p1(this%mxlay), this%at2_p1(this%mxlay), &
               this%at1_p2(this%mxlay), this%at2_p2(this%mxlay), &
               this%adt1_p1(this%mxlay), this%adt2_p1(this%mxlay), &
               this%adt1_p2(this%mxlay), this%adt2_p2(this%mxlay), &
               this%indxt_p1(this%mxlay), this%indxt_p2(this%mxlay), &
               this%indxp(this%mxlay), this%tavl(this%mxlay), &
               this%pavl(this%mxlay), this%wfix(this%mxlay), &
               this%dtu(this%mxlay), this%dtl(this%mxlay), &
               this%q(this%mxhmol,this%mxlay),  &
               this%w(this%mxhmol,this%mxlay), &
               this%dwqu(this%mxlay,this%mxhmol), &
               this%dwql(this%mxlay,this%mxhmol), &
               this%qobs(this%mxhmol), this%wobs(this%mxhmol), &
               this%dwquobs(this%mxhmol), this%dwqlobs(this%mxhmol), &
               this%abso(this%mxhmol,this%mxlay), &
               this%tautot(this%mxlay), this%dtaudtmp(this%mxlay), &
               this%absoobs(this%mxhmol), this%txdn(this%mxlev) , &
               this%txup(this%mxlev), this%bbar(this%mxlev), &
               this%dbbar(this%mxlev), this%draddtmp(this%mxlev), &
               this%draddtau(this%mxlev), this%drdw(this%mxlev), &
               this%draddtmpdw(this%mxlev), this%draddtmpuw(this%mxlev), &
               this%dp(this%mxlev), this%qcor(this%mxlev), &
               this%qr(this%mxlev), this%ploc(this%mxlay), &
               this%b2(this%mxlay), this%db2(this%mxlay), &
               this%ar(this%mxlay), this%br(this%mxlay), &
               this%ad(this%mxlay), this%bd(this%mxlay), &
               stat=iret)
      if ( iret /= 0 ) then
        write(error_unit,*) 'In file ',__FILE__,' at line ',__LINE__, &
          ' : Allocation error !'
        itmp = this%delete( )
        iret = -1
        return
      end if

      this%isinit = .true.
      iret = 0
    end function init

    subroutine forward_model(this,indata,outdata)
      implicit none
      class(aoss_ir) :: this
      type(atmospheric_data) , intent(in) :: indata
      type(synthetic_measure) , intent(out) :: outdata

    end subroutine forward_model

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
      ao%mxlay = 0
      ao%mxhmol = 0
      ao%mxlev = 0
    end subroutine zerodims

    subroutine cleanup(ao)
      implicit none
      type(aoss_ir) :: ao
      if ( allocated(ao%pavlref) ) then
        deallocate(ao%pavlref)
        deallocate(ao%sunrad)
        deallocate(ao%ap1)
        deallocate(ao%ap2)
        deallocate(ao%at1_p1)
        deallocate(ao%at2_p1)
        deallocate(ao%at1_p2)
        deallocate(ao%at2_p2)
        deallocate(ao%adt1_p1)
        deallocate(ao%adt2_p1)
        deallocate(ao%adt1_p2)
        deallocate(ao%adt2_p2)
        deallocate(ao%indxt_p1)
        deallocate(ao%indxt_p2)
        deallocate(ao%indxp)
        deallocate(ao%tavl)
        deallocate(ao%pavl)
        deallocate(ao%wfix)
        deallocate(ao%dtu)
        deallocate(ao%dtl)
        deallocate(ao%q)
        deallocate(ao%w)
        deallocate(ao%dwqu)
        deallocate(ao%dwql)
        deallocate(ao%qobs)
        deallocate(ao%wobs)
        deallocate(ao%dwquobs)
        deallocate(ao%dwqlobs)
        deallocate(ao%abso)
        deallocate(ao%tautot)
        deallocate(ao%dtaudtmp)
        deallocate(ao%absoobs)
        deallocate(ao%txdn)
        deallocate(ao%txup)
        deallocate(ao%bbar)
        deallocate(ao%dbbar)
        deallocate(ao%draddtmp)
        deallocate(ao%draddtau)
        deallocate(ao%drdw)
        deallocate(ao%draddtmpdw)
        deallocate(ao%draddtmpuw)
        deallocate(ao%dp)
        deallocate(ao%qcor)
        deallocate(ao%qr)
        deallocate(ao%ploc)
        deallocate(ao%b2)
        deallocate(ao%db2)
        deallocate(ao%ar)
        deallocate(ao%br)
        deallocate(ao%ad)
        deallocate(ao%bd)
      end if
    end subroutine cleanup

    subroutine planck_set(ao,vn,t,ts,nlay)
      type(aoss_ir) , intent(inout) :: ao
      integer , intent(in) :: nlay
      real , intent(in) :: vn , ts
      real , intent(in) , dimension(:) :: t
      integer :: l
      ao%v2 = vn
      do l = 1 , nlay
        call fplanck(vn,t(l),ao%b2(l),ao%db2(l))
      end do
      call fplanck(vn,ts,ao%bs2,ao%dbs2)
      is_planck_set = .true.
    end subroutine planck_set

    subroutine fplanck(vn,t,rad,draddt)
      !------------------------------------------------------------------------
      ! purpose: calculates planck function and its derivative  wrt temperature
      !------------------------------------------------------------------------
      implicit none
      real , intent(in) :: t , vn
      real , intent(out) :: rad , draddt
      real :: f1 , f2 , f3
      f1 = c1*vn**3
      f2 = c2*vn/t
      f3 = exp(f2)
      rad = f1/(f3-1.0E0)
      draddt =(rad*rad)*f2*f3/(f1*t)
    end subroutine fplanck

    subroutine planck_int(ao,vn,t,ts,nlay,radpl,db,bs,dbs)
      implicit none
      type(aoss_ir) , intent(inout) :: ao
      integer , intent(in) :: nlay
      real , intent(in) :: vn , ts
      real , intent(in) , dimension(:) :: t
      real , intent(out) :: bs , dbs
      real , intent(out) , dimension(:) :: radpl , db
      real :: bs1 , dbs1
      real :: v1 , v2n
      real , dimension(nlay) :: b1 , db1
      if ( .not. is_planck_set ) then
        call planck_set(ao,vn,t,ts,nlay)
      end if
      if ( vn >= ao%v2 ) then
        v2n = ao%v2 + dvint
        if ( vn >= v2n ) then
          call planck_set(ao,vn,t,ts,nlay)
          ao%v2 = vn
          v2n = ao%v2+dvint
        end if
        v1 = ao%v2
        b1(1:nlay)  = ao%b2(1:nlay)
        db1(1:nlay) = ao%db2(1:nlay)
        bs1 = ao%bs2
        dbs1 = ao%dbs2
        ao%v2 = v2n
        call planck_set(ao,ao%v2,t,ts,nlay)
        ao%ar(1:nlay) = (ao%b2(1:nlay)-b1(1:nlay))/dvint
        ao%br(1:nlay) = (b1(1:nlay)*ao%v2-ao%b2(1:nlay)*v1)/dvint
        ao%ad(1:nlay) = (ao%db2(1:nlay)-db1(1:nlay))/dvint
        ao%bd(1:nlay) = (db1(1:nlay)*ao%v2-ao%db2(1:nlay)*v1)/dvint
        ao%ars = (ao%bs2-bs1)/dvint
        ao%brs = (bs1*ao%v2-ao%bs2*v1)/dvint
        ao%ads = (ao%dbs2-dbs1)/dvint
        ao%bds = (dbs1*ao%v2-ao%dbs2*v1)/dvint
      end if
      radpl(1:nlay) = vn*ao%ar(1:nlay)+ao%br(1:nlay)
      db(1:nlay) = vn*ao%ad(1:nlay)+ao%bd(1:nlay)
      bs = vn*ao%ars+ao%brs
      dbs = vn*ao%ads+ao%bds
    end subroutine planck_int

    subroutine lint_log(xinp,pgrid,n0,p0,x0,x1)
      !---------------------------------------
      ! purpose: interpolation in log pressure
      !---------------------------------------
      integer , intent(in) :: n0
      real , intent(in) :: p0
      real , dimension(:) , intent(in) :: xinp
      real , dimension(:) , intent(in) :: pgrid
      real , intent(out) :: x0
      real , intent(out) , optional :: x1
      real :: xx
      xx = log(xinp(n0)/xinp(n0-1))/log(pgrid(n0)/pgrid(n0-1))
      x0 = xinp(n0-1)*(p0/pgrid(n0-1))**xx
      if ( present(x1) ) x1 = xinp(n0)*(p0/pgrid(n0))**xx
    end subroutine lint_log

    subroutine lpsum_log(pu,pl,xu,xl,scal,xint,dxu,dxl)
      !------------------------------------------------------------------
      ! purpose: computes average layer quantities (or integrated amount)
      !          using a log-x dependence on log-p.
      !------------------------------------------------------------------
      implicit none
      real , intent(in) :: pu , pl , scal
      real , intent(in) :: xu , xl
      real , intent(out) :: xint , dxu , dxl
      real :: hp , x0 , zeta , alza , alpha
      hp = log(pl/pu)
      x0 = pl*xl*hp
      zeta = pu*xu/(pl*xl)
      if ( abs(zeta-1.0E0) > epsiln ) then
        alza  = log(zeta)
        xint  = x0*(zeta-1.0E0)/alza
        alpha = zeta/(zeta-1.0E0)-1.0E0/alza
      else
        xint  = x0*2.0E0/(3.0E0-zeta)
        alpha = zeta/(3.0E0-zeta)
      end if
      xint = xint*scal
      dxu = xint*alpha/xu
      dxl = xint*(1.0E0-alpha)/xl
    end subroutine lpsum_log

end module mod_oss_ir

#ifdef TESTME

! mpif90 -O3 -g -mtune=native -c -I`nf-config --includedir` \
!         mod_forward_model_data.F90
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
