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

    !------------------------------------------------------------------------
    ! purpose: calculates planck function and its derivative wrt temperature
    !------------------------------------------------------------------------
    subroutine fplanck(vn,t,rad,draddt)
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

    !---------------------------------------
    ! purpose: interpolation in log pressure
    !---------------------------------------
    subroutine lint_log(xinp,pgrid,n0,p0,x0,x1)
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

    !------------------------------------------------------------------
    ! purpose: computes average layer quantities (or integrated amount)
    !          using a log-x dependence on log-p.
    !------------------------------------------------------------------
    subroutine lpsum_log(pu,pl,xu,xl,scal,xint,dxu,dxl)
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

    !---------------------------------------------------------------
    ! purpose:  given the luts, the function computes the layer
    !   optical depths. refer to ossrad_mw for the transmittance and
    !   brightness temperature calculation
    !---------------------------------------------------------------
    subroutine osstran(ao,n1,n2,nn)
      implicit none
      type(aoss_ir) , intent(inout) :: ao
      integer , intent(in) :: n1 , n2 , nn
      integer :: indxt1_p1 , indxt2_p1 , indxt3_p1
      integer :: l , ks , imol , indxp1
      real :: dft1_p1 , dft2_p1
      real :: abs0_p1 , abs0_p2 , dabs0_p1 , dabs0_p2
      real :: absh2ot1_p1 , absh2ot2_p1 , absh2ot3_p1
      real :: absh2o_p1 , absh2o_p2
      real :: dabsh2o_p1 , dabsh2o_p2
      real :: dabsh2odq_p1 , dabsh2odq_p2
      real :: at1_p1l , at2_p1l , adt1_p1l , adt2_p1l
      real :: adt1_p2l , adt2_p2l
      real :: ap1l , ap2l

      abs0_p2 = 0.0E0
      dabs0_p1 = 0.0E0
      dabs0_p2 = 0.0E0
      absh2o_p2 = 0.0E0
      dabsh2o_p2 = 0.0E0
      dabsh2odq_p2 = 0.0E0

      if ( n1 > 1 ) then
        indxt1_p1 = ao%indxt_p1obs-1
        indxt2_p1 = indxt1_p1+1
        indxt3_p1 = indxt1_p1+2
        indxp1 = ao%indxpobs
        at1_p1l = ao%at1_p1obs
        at2_p1l = ao%at2_p1obs
        adt1_p1l = ao%adt1_p1obs
        adt2_p1l = ao%adt2_p1obs
        ap1l = ao%ap1obs
        ap2l = ao%ap2obs
        !---fixed gases
        dft1_p1 = ao%ah%kfix(indxp1,indxt1_p1,nn) - &
                  ao%ah%kfix(indxp1,indxt3_p1,nn)
        dft2_p1 = ao%ah%kfix(indxp1,indxt2_p1,nn) - &
                  ao%ah%kfix(indxp1,indxt3_p1,nn)
        abs0_p1 = at1_p1l*dft1_p1 + at2_p1l*dft2_p1 + &
                  ao%ah%kfix(indxp1,indxt3_p1,nn)
        !---water vapor
        absh2ot1_p1 = ao%qobs(1) * ao%ah%dkh2o(indxp1,indxt1_p1,nn) + &
                      ao%ah%kh2o(indxp1,indxt1_p1,nn)
        absh2ot2_p1 = ao%qobs(1) * ao%ah%dkh2o(indxp1,indxt2_p1,nn) + &
                      ao%ah%kh2o(indxp1,indxt2_p1,nn)
        absh2ot3_p1 = ao%qobs(1) * ao%ah%dkh2o(indxp1,indxt3_p1,nn) + &
                      ao%ah%kh2o(indxp1,indxt3_p1,nn)
        dft1_p1 = absh2ot1_p1 - absh2ot3_p1
        dft2_p1 = absh2ot2_p1 - absh2ot3_p1
        absh2o_p1 = at1_p1l*dft1_p1 + at2_p1l*dft2_p1 + absh2ot3_p1
        dabsh2odq_p1 = at1_p1l*(ao%ah%dkh2o(indxp1,indxt1_p1,nn) -  &
                                ao%ah%dkh2o(indxp1,indxt3_p1,nn)) + &
                       at2_p1l*(ao%ah%dkh2o(indxp1,indxt2_p1,nn) -  &
                                ao%ah%dkh2o(indxp1,indxt3_p1,nn)) + &
                       ao%ah%dkh2o(indxp1,indxt3_p1,nn)
        dabsh2o_p1 = adt1_p1l*dft1_p1 + adt2_p1l*dft2_p1
        ao%absoobs(1) = -(abs0_p1*ap1l + abs0_p2*ap2l)
        ao%tautotobs = -ao%absoobs(1) * ao%wfixobs
        ao%dtaudtmpobs = (dabs0_p1*ap1l+dabs0_p2*ap2l) * ao%wfixobs
        ao%tautotobs  = ao%tautotobs + &
                        (absh2o_p1*ap1l+absh2o_p2*ap2l) * ao%wobs(1)
        ao%absoobs(1) = (dabsh2odq_p1*ap1l+dabsh2odq_p2*ap2l)*ao%qobs(1) + &
                     (absh2o_p1*ap1l+absh2o_p2*ap2l) + ao%absoobs(1)
        ao%dtaudtmpobs = ao%dtaudtmpobs + &
                      (dabsh2o_p1*ap1l+dabsh2o_p2*ap2l) * ao%wobs(1)
        !---variable gases
        do ks = 2 , ao%ah%nmols(nn)
          imol = ao%ah%imols(ks,nn)
          dft1_p1 = ao%ah%kvar(ks-1,indxp1,indxt1_p1,nn) - &
                    ao%ah%kvar(ks-1,indxp1,indxt3_p1,nn)
          dft2_p1 = ao%ah%kvar(ks-1,indxp1,indxt2_p1,nn) - &
                    ao%ah%kvar(ks-1,indxp1,indxt3_p1,nn)
          abs0_p1 = at1_p1l*dft1_p1+at2_p1l*dft2_p1 + &
                    ao%ah%kvar(ks-1,indxp1,indxt3_p1,nn)
          dabs0_p1 = adt1_p1l*dft1_p1 + adt2_p1l*dft2_p1
          ao%absoobs(ks) = (abs0_p1*ap1l + abs0_p2*ap2l)
          ao%tautotobs = ao%tautotobs + ao%absoobs(ks) * ao%wobs(imol)
          ao%dtaudtmpobs = ao%dtaudtmpobs + &
                    (dabs0_p1*ap1l+dabs0_p2*ap2l)*ao%wobs(imol)
          ao%absoobs(1) = ao%absoobs(1) - ao%absoobs(ks) * ao%qobs(imol)
        end do
      end if
      do l = 1 , n2
        indxt1_p1 = ao%indxt_p1(l)-1
        indxt2_p1 = indxt1_p1+1
        indxt3_p1 = indxt1_p1+2
        indxp1 = ao%indxp(l)
        at1_p1l = ao%at1_p1(l)
        at2_p1l = ao%at2_p1(l)
        adt1_p1l = ao%adt1_p1(l)
        adt2_p1l = ao%adt2_p1(l)
        adt1_p2l = ao%adt1_p2(l)
        adt2_p2l = ao%adt2_p2(l)
        ap1l = ao%ap1(l)
        ap2l = ao%ap2(l)
        !---fixed gases
        dft1_p1 = ao%ah%kfix(indxp1,indxt1_p1,nn) - &
                  ao%ah%kfix(indxp1,indxt3_p1,nn)
        dft2_p1 = ao%ah%kfix(indxp1,indxt2_p1,nn) - &
                  ao%ah%kfix(indxp1,indxt3_p1,nn)
        abs0_p1 = at1_p1l*dft1_p1 + at2_p1l*dft2_p1 + &
                  ao%ah%kfix(indxp1,indxt3_p1,nn)
        dabs0_p1 = adt1_p1l*dft1_p1 + adt2_p1l*dft2_p1
        !---water vapor
        absh2ot1_p1 = ao%q(1,l) * ao%ah%dkh2o(indxp1,indxt1_p1,nn) + &
                      ao%ah%kh2o(indxp1,indxt1_p1,nn)
        absh2ot2_p1 = ao%q(1,l) * ao%ah%dkh2o(indxp1,indxt2_p1,nn) + &
                      ao%ah%kh2o(indxp1,indxt2_p1,nn)
        absh2ot3_p1 = ao%q(1,l) * ao%ah%dkh2o(indxp1,indxt3_p1,nn) + &
                      ao%ah%kh2o(indxp1,indxt3_p1,nn)
        dft1_p1 = absh2ot1_p1 - absh2ot3_p1
        dft2_p1 = absh2ot2_p1 - absh2ot3_p1
        absh2o_p1 = at1_p1l*dft1_p1 + at2_p1l*dft2_p1 + absh2ot3_p1
        dabsh2odq_p1 = at1_p1l*(ao%ah%dkh2o(indxp1,indxt1_p1,nn) -  &
                                ao%ah%dkh2o(indxp1,indxt3_p1,nn)) + &
                       at2_p1l*(ao%ah%dkh2o(indxp1,indxt2_p1,nn) -  &
                                ao%ah%dkh2o(indxp1,indxt3_p1,nn)) + &
                       ao%ah%dkh2o(indxp1,indxt3_p1,nn)
        dabsh2o_p1 = adt1_p1l*dft1_p1 + adt2_p1l*dft2_p1
        ao%abso(1,l) = -(abs0_p1*ap1l + abs0_p2*ap2l)
        ao%tautot(l) = - ao%abso(1,l) * ao%wfix(l)
        ao%dtaudtmp(l) = (dabs0_p1*ap1l+dabs0_p2*ap2l) * ao%wfix(l)
        ao%tautot(l) = ao%tautot(l) + &
                    (absh2o_p1*ap1l+absh2o_p2*ap2l) * ao%w(1,l)
        ao%abso(1,l) = (dabsh2odq_p1*ap1l+dabsh2odq_p2*ap2l)*ao%q(1,l) + &
                    (absh2o_p1*ap1l+absh2o_p2*ap2l) + ao%abso(1,l)
        ao%dtaudtmp(l) = ao%dtaudtmp(l) + &
          (dabsh2o_p1*ap1l+dabsh2o_p2*ap2l)* ao%w(1,l)

        !---variable gases
        do ks = 2 , ao%ah%nmols(nn)
          imol = ao%ah%imols(ks,nn)
          dft1_p1 = ao%ah%kvar(ks-1,indxp1,indxt1_p1,nn) - &
                    ao%ah%kvar(ks-1,indxp1,indxt3_p1,nn)
          dft2_p1 = ao%ah%kvar(ks-1,indxp1,indxt2_p1,nn) - &
                    ao%ah%kvar(ks-1,indxp1,indxt3_p1,nn)
          abs0_p1 = at1_p1l*dft1_p1+at2_p1l*dft2_p1 + &
                    ao%ah%kvar(ks-1,indxp1,indxt3_p1,nn)
          dabs0_p1 = adt1_p1l*dft1_p1 + adt2_p1l*dft2_p1
          ao%abso(ks,l) = (abs0_p1*ap1l + abs0_p2*ap2l)
          ao%tautot(l) = ao%tautot(l) + ao%abso(ks,l) * ao%w(imol,l)
          ao%dtaudtmp(l) = ao%dtaudtmp(l) + &
                           (dabs0_p1*ap1l+dabs0_p2*ap2l)*ao%w(imol,l)
          ao%abso(1,l) = ao%abso(1,l) - ao%abso(ks,l) * ao%q(imol,l)
         end do
      end do
    end subroutine osstran

    !---------------------------------------------------------------
    ! purpose:  given the luts, the function computes the layer
    !   optical depths. refer to ossrad_mw for the transmittance and
    !   brightness temperature calculation
    !---------------------------------------------------------------
    subroutine osstran_vg(ao,n2,nn)
      implicit none
      type(aoss_ir) , intent(inout) :: ao
      integer , intent(in) :: nn , n2
      integer :: indxt1_p1 , indxt2_p1 , indxt3_p1
      integer :: indxt1_p2 , indxt2_p2 , indxt3_p2
      integer :: l , ks , imol , indxp1 , indxp2
      real :: dft1_p1 , dft2_p1 , dft1_p2 , dft2_p2
      real :: abs0_p1 , abs0_p2 , dabs0_p1 , dabs0_p2
      real :: absh2ot1_p1 , absh2ot2_p1 , absh2ot3_p1
      real :: absh2o_p1 , absh2o_p2
      real :: absh2ot1_p2 , absh2ot2_p2 , absh2ot3_p2
      real :: dabsh2o_p1 , dabsh2o_p2
      real :: dabsh2odq_p1 , dabsh2odq_p2
      real :: at1_p1l , at2_p1l , adt1_p1l , adt2_p1l
      real :: adt1_p2l , adt2_p2l
      real :: ap1l , ap2l
   
      do l = 1 , n2
        indxt1_p1 = ao%indxt_p1(l)-1
        indxt2_p1 = indxt1_p1+1
        indxt3_p1 = indxt1_p1+2
        indxp1 = ao%indxp(l)
        at1_p1l = ao%at1_p1(l)
        at2_p1l = ao%at2_p1(l)
        adt1_p1l = ao%adt1_p1(l)
        adt2_p1l = ao%adt2_p1(l)
        adt1_p2l = ao%adt1_p2(l)
        adt2_p2l = ao%adt2_p2(l)
        ap1l = ao%ap1(l)
        ap2l = ao%ap2(l)

        !---fixed gases
        dft1_p1 = ao%ah%kfix(indxp1,indxt1_p1,nn) - &
                  ao%ah%kfix(indxp1,indxt3_p1,nn)
        dft2_p1 = ao%ah%kfix(indxp1,indxt2_p1,nn) - &
                  ao%ah%kfix(indxp1,indxt3_p1,nn)
        abs0_p1 = at1_p1l*dft1_p1 + at2_p1l*dft2_p1 + &
                  ao%ah%kfix(indxp1,indxt3_p1,nn)
        dabs0_p1 = adt1_p1l*dft1_p1 + adt2_p1l*dft2_p1
        !---water vapor
        absh2ot1_p1 = ao%q(1,l) * ao%ah%dkh2o(indxp1,indxt1_p1,nn) + &
                      ao%ah%kh2o(indxp1,indxt1_p1,nn)
        absh2ot2_p1 = ao%q(1,l) * ao%ah%dkh2o(indxp1,indxt2_p1,nn) + &
                      ao%ah%kh2o(indxp1,indxt2_p1,nn)
        absh2ot3_p1 = ao%q(1,l) * ao%ah%dkh2o(indxp1,indxt3_p1,nn) + &
                      ao%ah%kh2o(indxp1,indxt3_p1,nn)
        dft1_p1 = absh2ot1_p1 - absh2ot3_p1
        dft2_p1 = absh2ot2_p1 - absh2ot3_p1
        absh2o_p1 = at1_p1l*dft1_p1 + at2_p1l*dft2_p1 + absh2ot3_p1
        dabsh2odq_p1 = at1_p1l*(ao%ah%dkh2o(indxp1,indxt1_p1,nn) -  &
                                ao%ah%dkh2o(indxp1,indxt3_p1,nn)) + &
                       at2_p1l*(ao%ah%dkh2o(indxp1,indxt2_p1,nn) -  &
                                ao%ah%dkh2o(indxp1,indxt3_p1,nn)) + &
                       ao%ah%dkh2o(indxp1,indxt3_p1,nn)
        dabsh2o_p1 = adt1_p1l*dft1_p1 + adt2_p1l*dft2_p1
        indxt1_p2 = ao%indxt_p2(l)-1
        indxt2_p2 = indxt1_p2+1
        indxt3_p2 = indxt1_p2+2
        indxp2 = indxp1+1
        !---fixed gases
        dft1_p2 = ao%ah%kfix(indxp2,indxt1_p2,nn) - &
                  ao%ah%kfix(indxp2,indxt3_p2,nn)
        dft2_p2 = ao%ah%kfix(indxp2,indxt2_p2,nn) - &
                  ao%ah%kfix(indxp2,indxt3_p2,nn)
        abs0_p2 = ao%at1_p2(l)*dft1_p2 + ao%at2_p2(l)*dft2_p2 + &
                  ao%ah%kfix(indxp2,indxt3_p2,nn)
        dabs0_p2 = adt1_p2l*dft1_p2 + adt2_p2l*dft2_p2
        !---water vapor
        absh2ot1_p2 = ao%q(1,l) * ao%ah%dkh2o(indxp2,indxt1_p2,nn) + &
                      ao%ah%kh2o(indxp2,indxt1_p2,nn)
        absh2ot2_p2 = ao%q(1,l) * ao%ah%dkh2o(indxp2,indxt2_p2,nn) + &
                      ao%ah%kh2o(indxp2,indxt2_p2,nn)
        absh2ot3_p2 = ao%q(1,l) * ao%ah%dkh2o(indxp2,indxt3_p2,nn) + &
                      ao%ah%kh2o(indxp2,indxt3_p2,nn)
        dft1_p2 = absh2ot1_p2 - absh2ot3_p2
        dft2_p2 = absh2ot2_p2 - absh2ot3_p2
        absh2o_p2 = ao%at1_p2(l)*dft1_p2 + ao%at2_p2(l)*dft2_p2 + absh2ot3_p2
        dabsh2odq_p2 = ao%at1_p2(l)*(ao%ah%dkh2o(indxp2,indxt1_p2,nn) -  &
                                  ao%ah%dkh2o(indxp2,indxt3_p2,nn)) + &
                       ao%at2_p2(l)*(ao%ah%dkh2o(indxp2,indxt2_p2,nn) -  &
                                  ao%ah%dkh2o(indxp2,indxt3_p2,nn)) + &
                       ao%ah%dkh2o(indxp2,indxt3_p2,nn)
        dabsh2o_p2 = adt1_p2l*dft1_p2 + adt2_p2l*dft2_p2
        ao%abso(1,l) = -(abs0_p1*ap1l + abs0_p2*ap2l)
        ao%tautot(l) = - ao%abso(1,l) * ao%wfix(l)
        ao%dtaudtmp(l) = (dabs0_p1*ap1l+dabs0_p2*ap2l) * ao%wfix(l)
        ao%tautot(l) = ao%tautot(l) + &
                       (absh2o_p1*ap1l+absh2o_p2*ap2l) * ao%w(1,l)
        ao%abso(1,l) = (dabsh2odq_p1*ap1l+dabsh2odq_p2*ap2l) * ao%q(1,l) + &
                    (absh2o_p1*ap1l+absh2o_p2*ap2l) + ao%abso(1,l)
        ao%dtaudtmp(l) = ao%dtaudtmp(l) + &
                        (dabsh2o_p1*ap1l+dabsh2o_p2*ap2l) * ao%w(1,l)

        !---variable gases
        do ks = 2 , ao%ah%nmols(nn)
          imol = ao%ah%imols(ks,nn)
          dft1_p1 = ao%ah%kvar(ks-1,indxp1,indxt1_p1,nn) - &
                    ao%ah%kvar(ks-1,indxp1,indxt3_p1,nn)
          dft2_p1 = ao%ah%kvar(ks-1,indxp1,indxt2_p1,nn) - &
                    ao%ah%kvar(ks-1,indxp1,indxt3_p1,nn)
          abs0_p1 = at1_p1l*dft1_p1+at2_p1l*dft2_p1 + &
                    ao%ah%kvar(ks-1,indxp1,indxt3_p1,nn)
          dabs0_p1 = adt1_p1l*dft1_p1 + adt2_p1l*dft2_p1
          dft1_p2 = ao%ah%kvar(ks-1,indxp2,indxt1_p2,nn) - &
                    ao%ah%kvar(ks-1,indxp2,indxt3_p2,nn)
          dft2_p2 = ao%ah%kvar(ks-1,indxp2,indxt2_p2,nn) - &
                    ao%ah%kvar(ks-1,indxp2,indxt3_p2,nn)
          abs0_p2 = ao%at1_p2(l)*dft1_p2+ao%at2_p2(l)*dft2_p2 + &
                    ao%ah%kvar(ks-1,indxp2,indxt3_p2,nn)
          dabs0_p2 = adt1_p2l*dft1_p2 + adt2_p2l*dft2_p2
          ao%abso(ks,l) = (abs0_p1*ap1l + abs0_p2*ap2l)
          ao%tautot(l) = ao%tautot(l) + ao%abso(ks,l) * ao%w(imol,l)
          ao%dtaudtmp(l) = ao%dtaudtmp(l) + &
                           (dabs0_p1*ap1l+dabs0_p2*ap2l) * ao%w(imol,l)
          ao%abso(1,l) = ao%abso(1,l) - ao%abso(ks,l) * ao%q(imol,l)
         end do
      end do
    end subroutine osstran_vg

    !---------------------------------------------------------------------
    ! purpose:  this subroutine calculates the average temperature and
    !   molecular amounts for all layers for given profiles of temperature
    !   and molecular concentrations. it also calculates the derivatives
    !   of tavl with respect to a change in the lower and upper boundary
    !   temperatures and the derivatives of the molecular amounts with
    !   respect to a change in the mixing ratios at the layer
    !   boundaries. molecular amounts are in molec./cm**2.
    !   integration assumes that t is linear in z (lnt linear in lnp)
    !   and lnq linear in lnp.
    !---------------------------------------------------------------------
    subroutine fpath_ir(ao,pref,xg,pobs,nobs,nlast,tsfc)
      implicit none
      type(aoss_ir) , intent(inout) :: ao
      integer , intent(in) :: nobs , nlast
      real , intent(in) :: pobs
      real , intent(in) , dimension(:) :: pref
      real , intent(in) , dimension(:) :: xg
      real , intent(out) :: tsfc
      integer :: nsurf , l , ih2o , k , ixoff , n , nend
      real :: psfc , tobs
      real :: scal , scal2 , qsfc , wtot , wdry
      real :: qtmp
      real :: dpobs = 0.0E0
      real :: scalobs2 = 0.0E0
      real :: scalobs
      real :: wdryobs , wtotobs
      real :: tobsdum , qtmpobs

      nsurf = nlast+1

      !------------------------------------------------------------------------
      !     compute average pressure for the layers
      !------------------------------------------------------------------------
      ao%ploc(1:nsurf) = pref(1:nsurf)
      nend = nlast-1
      psfc = xg(ipsfcg)
      !------------------------------------------------------------------------
      !     compute average temperature for the layers
      !------------------------------------------------------------------------
      do l = 1 , nend
        ao%dp(l) = (ao%ploc(l+1)-ao%ploc(l))
        scal = 1.0E0/ao%dp(l)
        call lpsum_log(ao%ploc(l),ao%ploc(l+1),xg(itempg+l-1),xg(itempg+l), &
                       scal,ao%tavl(l),ao%dtu(l),ao%dtl(l))
      end do
      tsfc = xg(itempg+nend)
      tobs = xg(itempg+nobs-1)
      !---surface layer
      ao%dp(nlast) = (psfc-ao%ploc(nsurf-1))
      scal = 1.0E0/ao%dp(nlast)
      scal2 = log(psfc/ao%ploc(nsurf))/log(ao%ploc(nsurf-1)/ao%ploc(nsurf))
      call lint_log(xg,ao%ploc,nsurf,psfc,tsfc)
      call lpsum_log(ao%ploc(nsurf-1),psfc,xg(nsurf-1),tsfc, &
                     scal,ao%tavl(nlast),ao%dtu(nlast),ao%dtl(nlast))
      ao%dtu(nlast) = ao%dtu(nlast)+ao%dtl(nlast)*scal2*tsfc/xg(nsurf-1)
      ao%dtl(nlast) = ao%dtl(nlast)*(1.0E0-scal2)*tsfc/xg(nsurf)
      if ( nobs > 1 ) then
        dpobs = ao%ploc(nobs)-pobs
        scalobs = 1.0E0/dpobs
        scalobs2 = log(pobs/ao%ploc(nobs-1))/log(ao%ploc(nobs)/ao%ploc(nobs-1))
        call lint_log(xg,ao%ploc,nobs,pobs,tobsdum,tobs)
        call lpsum_log(pobs,ao%ploc(nobs),tobs,xg(nobs),&
                       scalobs,ao%tavlobs,ao%dtuobs,ao%dtlobs)
        ao%dtuobs = ao%dtuobs*tobs/xg(nobs-1)*(1.0E0-scalobs2)
        ao%dtlobs = ao%dtlobs+ao%dtuobs*scalobs2*tobs/xg(nobs)
      end if
      !--------------------------------------------------------------------
      ! calculate amounts for individual species and derivatives wrt mixing
      ! ratios for retrieved constituents.
      !--------------------------------------------------------------------
      scal = rair
      ih2o = ao%ah%molid(1)-1
      ao%qcor(1:nsurf) = 1.0E0/(1.0E0+xg(ih2o+1:ih2o+nsurf))
      do k = 1 , ao%mxhmol
        ixoff = ao%ah%molid(k)-1
        !---transform mix. ratios into mass fractions
        !--- (relative to total air mass)
        ao%qr(1:nsurf) = xg(ixoff+1:ixoff+nsurf)*ao%qcor(1:nsurf)
        do l = 1 , nend
          call lpsum_log(ao%ploc(l),ao%ploc(l+1),ao%qr(l),ao%qr(l+1),&
                         scal,ao%w(k,l),ao%dwqu(l,k),ao%dwql(l,k))
        end do
        !---surface layer
        call lint_log(ao%qr,pref,nsurf,psfc,qsfc)
        call lpsum_log(pref(nsurf-1),psfc,ao%qr(nsurf-1),qsfc,&
                       scal,ao%w(k,nlast),ao%dwqu(nlast,k),ao%dwql(nlast,k))
        ao%dwqu(nlast,k) = ao%dwqu(nlast,k) + &
                           ao%dwql(nlast,k)*scal2*qsfc/ao%qr(nsurf-1)
        ao%dwql(nlast,k) = ao%dwql(nlast,k)*(1.0E0-scal2)*qsfc/ao%qr(nsurf)
        if ( nobs > 1 ) then
          call lint_log(ao%qr,ao%ploc,nobs,pobs,qtmpobs,qtmp)
          call lpsum_log(pobs,ao%ploc(nobs),qtmp,ao%qr(nobs),&
                         scal,ao%wobs(k),ao%dwquobs(k),ao%dwqlobs(k))
          ao%dwquobs(k) = ao%dwquobs(k)*qtmp/ao%qr(nobs-1)*(1.0E0-scalobs2)
          ao%dwqlobs(k) = ao%dwqlobs(k)+ao%dwquobs(k)*scalobs2*qtmp/ao%qr(nobs)
        end if
      end do
      !---derivatives of amount wrt dry mixing ratios
      do n = 1 , nlast
        ao%dwqu(n,1) = ao%dwqu(n,1)*ao%qcor(n)**2
        ao%dwql(n,1) = ao%dwql(n,1)*ao%qcor(n+1)**2
        ao%dwqu(n,2:ao%mxhmol) = ao%dwqu(n,2:ao%mxhmol)*ao%qcor(n)
        ao%dwql(n,2:ao%mxhmol) = ao%dwql(n,2:ao%mxhmol)*ao%qcor(n+1)
      end do
      ao%dwquobs(1) = ao%dwquobs(1)*ao%qcor(nobs-1)**2
      ao%dwqlobs(1) = ao%dwqlobs(1)*ao%qcor(nobs)**2
      ao%dwquobs(2:ao%mxhmol) = ao%dwquobs(2:ao%mxhmol)*ao%qcor(nobs-1)
      ao%dwqlobs(2:ao%mxhmol) = ao%dwqlobs(2:ao%mxhmol)*ao%qcor(nobs)
      !---compute fix gas amount and average layer mixing ratios
      do l = 1 , nlast
        wtot = ao%dp(l)*rair
        wdry = wtot-ao%w(1,l)
        ao%wfix(l) = wdry
        ao%q(1,l) = ao%w(1,l)/wtot
        ao%q(2:ao%mxhmol,l) = ao%w(2:ao%mxhmol,l)/wdry
      end do
      if ( nobs > 1 ) then
        wtotobs = dpobs*rair
        wdryobs = wtotobs-ao%wobs(1)
        ao%wfixobs = wdryobs
        ao%qobs(1) = ao%wobs(1)/wtotobs
        ao%qobs(2:ao%mxhmol) = ao%wobs(2:ao%mxhmol)/wdryobs
      end if
    end subroutine fpath_ir

    !---------------------------------------------------------------------
    ! purpose:  this subroutine calculates the average temperature and
    !   molecular amounts for all layers for given profiles of temperature
    !   and molecular concentrations. it also calculates the derivatives
    !   of tavl with respect to a change in the lower and upper boundary
    !   temperatures and the derivatives of the molecular amounts with
    !   respect to a change in the mixing ratios at the layer
    !   boundaries. molecular amounts are in molec./cm**2.
    !   integration assumes that t is linear in z (lnt linear in lnp)
    !   and lnq linear in lnp.
    !---------------------------------------------------------------------
    subroutine fpath_ir_vg(ao,xg,pobs,nobs,nlast,tsfc,puser)
      implicit none
      type(aoss_ir) , intent(inout) :: ao
      integer , intent(in) :: nobs , nlast
      real , intent(in) :: pobs
      real , intent(in) , dimension(:) :: xg
      real , intent(in) , dimension(:) :: puser
      real , intent(out) :: tsfc
      integer :: nsurf , l , ih2o , k , ixoff , n , nend
      real :: scal , wtot , wdry , tobs
      real :: wdryobs , wtotobs , dpobs

      ao%wfix = 0.0E0
      ao%q = 0.0E0
      ao%w = 0.0E0
      ao%dtu = 0.0E0
      ao%dtl = 0.0E0
      ao%dwqu = 0.0E0
      ao%dwql = 0.0E0

      nsurf = nlast+1
      !------------------------------------------------------------------------
      !     compute average pressure for the layers
      !------------------------------------------------------------------------
      ao%pavl(1:nlast) = 0.5E0*(puser(1:nlast)+puser(2:nsurf))
      ao%ploc(1:nsurf) = puser(1:nsurf)
      nend = nlast
      !------------------------------------------------------------------------
      !     compute average temperature for the layers
      !------------------------------------------------------------------------
      do l = 1 , nend
        ao%dp(l) = (ao%ploc(l+1)-ao%ploc(l))
        scal = 1.0E0/ao%dp(l)
        call lpsum_log(ao%ploc(l),ao%ploc(l+1),xg(itempg+l-1),xg(itempg+l), &
                       scal,ao%tavl(l),ao%dtu(l),ao%dtl(l))
      end do
      tsfc = xg(itempg+nend)
      tobs = xg(itempg+nobs-1)
      !--------------------------------------------------------------------
      ! calculate amounts for individual species and derivatives wrt mixing
      ! ratios for retrieved constituents.
      !--------------------------------------------------------------------
      scal = rair
      ih2o = ao%ah%molid(1)-1
      ao%qcor(1:nsurf) = 1.0E0/(1.0E0+xg(ih2o+1:ih2o+nsurf))
      do k = 1 , ao%mxhmol
        ixoff = ao%ah%molid(k)-1
        !---transform mix. ratios into mass fractions
        !--- (relative to total air mass)
        ao%qr(1:nsurf) = xg(ixoff+1:ixoff+nsurf)*ao%qcor(1:nsurf)
        do l = 1 , nend
          call lpsum_log(ao%ploc(l),ao%ploc(l+1),ao%qr(l),ao%qr(l+1),&
                         scal,ao%w(k,l),ao%dwqu(l,k),ao%dwql(l,k))
        end do
      end do
      !---derivatives of amount wrt dry mixing ratios
      do n = 1 , nlast
        ao%dwqu(n,1) = ao%dwqu(n,1)*ao%qcor(n)**2
        ao%dwql(n,1) = ao%dwql(n,1)*ao%qcor(n+1)**2
        ao%dwqu(n,2:ao%mxhmol) = ao%dwqu(n,2:ao%mxhmol)*ao%qcor(n)
        ao%dwql(n,2:ao%mxhmol) = ao%dwql(n,2:ao%mxhmol)*ao%qcor(n+1)
      end do
      if ( nobs > 1 ) then
        ao%dwquobs(1) = ao%dwquobs(1)*ao%qcor(nobs-1)**2
        ao%dwqlobs(1) = ao%dwqlobs(1)*ao%qcor(nobs)**2
        ao%dwquobs(2:ao%mxhmol) = ao%dwquobs(2:ao%mxhmol)*ao%qcor(nobs-1)
        ao%dwqlobs(2:ao%mxhmol) = ao%dwqlobs(2:ao%mxhmol)*ao%qcor(nobs)
      end if
      !---compute fix gas amount and average layer mixing ratios
      do l = 1 , nlast
        wtot = ao%dp(l)*rair
        wdry = wtot-ao%w(1,l)
        ao%wfix(l) = wdry
        ao%q(1,l) = ao%w(1,l)/wtot
        ao%q(2:ao%mxhmol,l) = ao%w(2:ao%mxhmol,l)/wdry
      end do
      if ( nobs > 1 ) then
        dpobs = ao%ploc(nobs)-pobs
        wtotobs = dpobs*rair
        wdryobs = wtotobs-ao%wobs(1)
        ao%wfixobs = wdryobs
        ao%qobs(1) = ao%wobs(1)/wtotobs
        ao%qobs(2:ao%mxhmol) = ao%wobs(2:ao%mxhmol)/wdryobs
      end if
    end subroutine fpath_ir_vg

    !---------------------------------------------------------
    ! purpose: computes items related to the viewing geometry.
    !---------------------------------------------------------
    subroutine setpath_ir(ao,pref,xg,pobs,obsang,sunang,nsurf,nobs, &
        n1,n2,sun,umu,umu0,delphi)
      implicit none
      type(aoss_ir) , intent(inout) :: ao
      real , intent(in) , dimension(:) :: pref , xg
      real , intent(in) :: obsang , sunang , pobs
      integer , intent(out) :: nobs , nsurf , n1 , n2
      logical , intent(out) :: sun
      real , intent(out) :: umu , umu0 , delphi
      real :: psurf , viewang
      integer :: i , icnt
      !---compute the surface level
      !nobs=1
      psurf = xg(ipsfcg)
      nlevloop: &
      do i = 2 , ao%mxlev-1
        if ( pref(i) >= psurf ) exit nlevloop
      end do nlevloop
      nsurf = i
      icnt = 1
      do while ( pobs > pref(icnt) )
        icnt = icnt+1
      end do
      nobs = icnt
      !---set viewing geometry parameters
      delphi = 0.0E0
      viewang = obsang
      if ( lookup ) then
        n1 = nsurf
        n2 = nobs-1
      else
        if ( viewang > 90.0E0 ) then
          call fatal(__FILE__,__LINE__,'eia must be less than 90')
        end if
        n1 = nobs
        n2 = nsurf-1
      end if
      umu = cos(viewang*deg2rad)
      if ( sunang < 80.0E0 ) then
        sun = .true.
        umu0 = cos(sunang*deg2rad)
      else
        sun = .false.
        umu0 = 1.0E0
      end if
    end subroutine setpath_ir

    !---------------------------------------------------------
    ! purpose: computes items related to the viewing geometry.
    !---------------------------------------------------------
    subroutine setpath_ir_vg(ao,pobs,obsang,sunang,nsurf,nobs, &
        n1,n2,sun,umu,umu0,delphi,puser)
      implicit none
      type(aoss_ir) , intent(inout) :: ao
      real , intent(in) :: obsang , sunang , pobs
      real , intent(in) , dimension(:) :: puser
      integer , intent(out) :: nobs , nsurf , n1 , n2
      logical , intent(out) :: sun
      real , intent(out) :: umu , umu0 , delphi
      real :: viewang
      integer :: icnt
      !---compute the surface level
      nsurf = size(puser)
      icnt = 1
      do while ( abs(pobs - puser(icnt)) > 1.0E-4 )
        icnt = icnt+1
      end do
      nobs = icnt
      !---set viewing geometry parameters
      delphi = 0.0E0
      viewang = obsang
      if ( lookup ) then
        n1 = nsurf
        n2 = nobs-1
      else
        if ( viewang > 90.0E0 ) then
          call fatal(__FILE__,__LINE__,'eia must be less than 90')
        end if
        n1 = nobs
        n2 = nsurf-1
      end if
      umu = cos(viewang*deg2rad)
      if ( sunang < 80.0E0 ) then
        sun = .true.
        umu0 = cos(sunang*deg2rad)
      else
        sun = .false.
        umu0 = 1.0E0
      end if
    end subroutine setpath_ir_vg

    subroutine settabindx_irobs(ao,n1)
      implicit none
      type(aoss_ir) , intent(inout) :: ao
      integer , intent(in) :: n1
      integer :: i , j1 , lp1
      real :: dent1_p1 , dent2_p1

      lp1 = n1-1
      ao%indxpobs = lp1
      ao%ap2obs = 0.0E0
      ao%ap1obs = 1.0E0
      temploop1: &
      do i = 3 , ao%ah%ntmpod-2
        if ( ao%ah%tmptab(i,lp1) >= ao%tavlobs ) exit temploop1
      end do temploop1
      if ( abs(ao%tavlobs-ao%ah%tmptab(i-1,lp1)) <= &
           abs(ao%tavlobs-ao%ah%tmptab(i,lp1)) ) then
        j1 = i-1
      else
        j1 = i
      end if
      ao%indxt_p1obs = j1
      dent1_p1 = (ao%ah%tmptab(j1-1,lp1)-ao%ah%tmptab(j1,lp1)) * &
                 (ao%ah%tmptab(j1-1,lp1)-ao%ah%tmptab(j1+1,lp1))
      dent2_p1 = (ao%ah%tmptab(j1,lp1)-ao%ah%tmptab(j1-1,lp1)) * &
                 (ao%ah%tmptab(j1,lp1)-ao%ah%tmptab(j1+1,lp1))
      ao%at1_p1obs = (ao%tavlobs-ao%ah%tmptab(j1,lp1)) * &
                  (ao%tavlobs-ao%ah%tmptab(j1+1,lp1))/dent1_p1
      ao%at2_p1obs = (ao%tavlobs-ao%ah%tmptab(j1-1,lp1)) * &
                  (ao%tavlobs-ao%ah%tmptab(j1+1,lp1))/dent2_p1
      ao%adt1_p1obs = (2.0E0*ao%tavlobs-ao%ah%tmptab(j1,lp1) - &
                       ao%ah%tmptab(j1+1,lp1))/dent1_p1
      ao%adt2_p1obs = (2.0E0*ao%tavlobs-ao%ah%tmptab(j1-1,lp1) - &
                       ao%ah%tmptab(j1+1,lp1))/dent2_p1
    end subroutine settabindx_irobs

    !--------------------------------------------------------------------
    ! purpose: computes temperature/water vapor indexes and interpolation
    !          coefficients.
    !--------------------------------------------------------------------
    subroutine settabindx_ir(ao,n2)
      implicit none
      type(aoss_ir) , intent(inout) :: ao
      integer , intent(in) :: n2
      integer :: l , i , j1 , lp1
      real :: dent1_p1 , dent2_p1

      !---compute coefficients for temperature interpolation of ods
      !!! indxp: is related to the index of od pressure layer located just above
      !!!  (by altitude) the current profile layer
      do l = 1 , n2
        lp1 = l
        ao%ap2(l) = 0.0E0
        ao%ap1(l) = 1.0E0
        ao%indxp(l) = lp1
        !!! indxt_p1: is related to the midst temperature index within the
        !!!           upper (by altitude) od pressure layer indxp
        !!! indxt_p2: is related to the midst temperature index within the
        !!!           lower (by altitude) od pressure layer indxp+1
        temploop1: &
        do i = 3 , ao%ah%ntmpod-2
          if ( ao%ah%tmptab(i,lp1) >= ao%tavl(l) ) exit temploop1
        end do temploop1
        if ( abs(ao%tavl(l)-ao%ah%tmptab(i-1,lp1)) <= &
             abs(ao%tavl(l)-ao%ah%tmptab(i,lp1)) ) then
          j1 = i-1
        else
          j1 = i
        end if
        ao%indxt_p1(l) = j1
        dent1_p1 = (ao%ah%tmptab(j1-1,lp1)-ao%ah%tmptab(j1,lp1)) * &
                   (ao%ah%tmptab(j1-1,lp1)-ao%ah%tmptab(j1+1,lp1))
        dent2_p1 = (ao%ah%tmptab(j1,lp1)-ao%ah%tmptab(j1-1,lp1)) * &
                   (ao%ah%tmptab(j1,lp1)-ao%ah%tmptab(j1+1,lp1))
        ao%at1_p1(l) = (ao%tavl(l)-ao%ah%tmptab(j1,lp1)) * &
                       (ao%tavl(l)-ao%ah%tmptab(j1+1,lp1))/dent1_p1
        ao%at2_p1(l) = (ao%tavl(l)-ao%ah%tmptab(j1-1,lp1)) * &
                       (ao%tavl(l)-ao%ah%tmptab(j1+1,lp1))/dent2_p1
        ao%adt1_p1(l) = (2.0E0*ao%tavl(l)-ao%ah%tmptab(j1,lp1) - &
                         ao%ah%tmptab(j1+1,lp1))/dent1_p1
        ao%adt2_p1(l) = (2.0E0*ao%tavl(l)-ao%ah%tmptab(j1-1,lp1) - &
                         ao%ah%tmptab(j1+1,lp1))/dent2_p1
      end do
    end subroutine settabindx_ir

    !--------------------------------------------------------------------
    ! purpose: computes temperature/water vapor indexes and interpolation
    !          coefficients.
    !--------------------------------------------------------------------
    subroutine settabindx_ir_vg(ao,n2)
      implicit none
      type(aoss_ir) , intent(inout) :: ao
      integer , intent(in) :: n2
      integer :: l , i , j1 , j2 , lp1 , lp2
      real :: dent1_p1 , dent2_p1 , dent1_p2 , dent2_p2

      !---compute coefficients for temperature interpolation of ods
      !!! indxp: is related to the index of od pressure layer located
      !!! just above (by altitude) the current profile layer
      do l = 1 , n2
        lp1 = 1
        do while ( ao%pavlref(lp1) < ao%pavl(l) )
          lp1 = lp1+1
        end do
        lp2 = lp1
        lp1 = lp1-1
        if ( lp1 == 0 ) then
          lp1 = 1
          lp2 = 2
        end if
        if ( lp2 >= ao%mxlev ) then
          lp1 = ao%mxlev-2
          lp2 = ao%mxlev-1
        end if
        !!! if ao%pavl(l) becomes located between lp1 and lp2:
        ao%ap2(l) = (ao%pavl(l)-ao%pavlref(lp1)) / &
                    (ao%pavlref(lp2)-ao%pavlref(lp1))
        ao%ap1(l) = (ao%pavlref(lp2)-ao%pavl(l)) / &
                    (ao%pavlref(lp2)-ao%pavlref(lp1))
        ao%indxp(l) = lp1
        !!! indxt_p1: is related to the midst temperature index within the
        !!!           upper (by altitude) od pressure layer indxp
        !!! indxt_p2: is related to the midst temperature index within the
        !!!           lower (by altitude) od pressure layer indxp+1
        temploop1: &
        do i = 3 , ao%ah%ntmpod-2
          if ( ao%ah%tmptab(i,lp1) >= ao%tavl(l) ) exit temploop1
        end do temploop1
        if ( abs(ao%tavl(l)-ao%ah%tmptab(i-1,lp1)) <= &
             abs(ao%tavl(l)-ao%ah%tmptab(i,lp1)) ) then
          j1 = i-1
        else
          j1 = i
        end if
        ao%indxt_p1(l) = j1
        dent1_p1 = (ao%ah%tmptab(j1-1,lp1)-ao%ah%tmptab(j1,lp1)) * &
                   (ao%ah%tmptab(j1-1,lp1)-ao%ah%tmptab(j1+1,lp1))
        dent2_p1 = (ao%ah%tmptab(j1,lp1)-ao%ah%tmptab(j1-1,lp1)) * &
                   (ao%ah%tmptab(j1,lp1)-ao%ah%tmptab(j1+1,lp1))
        ao%at1_p1(l) = (ao%tavl(l)-ao%ah%tmptab(j1,lp1)) * &
                       (ao%tavl(l)-ao%ah%tmptab(j1+1,lp1))/dent1_p1
        ao%at2_p1(l) = (ao%tavl(l)-ao%ah%tmptab(j1-1,lp1)) * &
                       (ao%tavl(l)-ao%ah%tmptab(j1+1,lp1))/dent2_p1
        ao%adt1_p1(l) = (2.0E0*ao%tavl(l)-ao%ah%tmptab(j1,lp1) - &
                        ao%ah%tmptab(j1+1,lp1))/dent1_p1
        ao%adt2_p1(l) = (2.0E0*ao%tavl(l)-ao%ah%tmptab(j1-1,lp1) - &
                        ao%ah%tmptab(j1+1,lp1))/dent2_p1

        temploop2: &
        do i = 3 , ao%ah%ntmpod-2
          if ( ao%ah%tmptab(i,lp2) >= ao%tavl(l) ) exit temploop2
        end do temploop2
        if ( abs(ao%tavl(l)-ao%ah%tmptab(i-1,lp2)) <= &
             abs(ao%tavl(l)-ao%ah%tmptab(i,lp2)) ) then
          j2 = i-1
        else
          j2 = i
        end if
        ao%indxt_p2(l) = j2
        !!! aty_px: interpolation coefficients, which are related to pressure
        !!!         layer (x=1,or,2, upper, or, lower) and temperature
        !!!         interpolation points (y=1,or,2, lower,or, midst)
        dent1_p2 = (ao%ah%tmptab(j2-1,lp2)-ao%ah%tmptab(j2,lp2)) * &
                   (ao%ah%tmptab(j2-1,lp2)-ao%ah%tmptab(j2+1,lp2))
        dent2_p2 = (ao%ah%tmptab(j2,lp2)-ao%ah%tmptab(j2-1,lp2)) * &
                   (ao%ah%tmptab(j2,lp2)-ao%ah%tmptab(j2+1,lp2))
        ao%at1_p2(l) = (ao%tavl(l)-ao%ah%tmptab(j2,lp2)) * &
                       (ao%tavl(l)-ao%ah%tmptab(j2+1,lp2))/dent1_p2
        ao%at2_p2(l) = (ao%tavl(l)-ao%ah%tmptab(j2-1,lp2)) * &
                       (ao%tavl(l)-ao%ah%tmptab(j2+1,lp2))/dent2_p2
        ao%adt1_p2(l) = (2.0E0*ao%tavl(l)-ao%ah%tmptab(j2,lp2) - &
                         ao%ah%tmptab(j2+1,lp2))/dent1_p2
        ao%adt2_p2(l) = (2.0E0*ao%tavl(l)-ao%ah%tmptab(j2-1,lp2) - &
                         ao%ah%tmptab(j2+1,lp2))/dent2_p2
      end do
    end subroutine settabindx_ir_vg

    !-------------------------------------
    ! purpose: interpolation in wavenumber
    !-------------------------------------
    subroutine vinterp(datin,vn,sfgrd,nsf,nxdim,ip0,coefint,datout)
      implicit none
      integer , intent(in) :: nxdim , nsf
      real , intent(in) , dimension(:,:) :: datin
      real , intent(in) :: vn
      real , intent(in) , dimension(:) :: sfgrd
      integer , intent(inout) :: ip0
      real , intent(out) :: coefint
      real , intent(out) , dimension(:) :: datout
      integer :: ip , i
      iploop: &
      do ip = ip0 , nsf-1
        if ( sfgrd(ip) > vn ) exit iploop
      end do iploop
      ip0 = ip
      coefint = (vn-sfgrd(ip0-1))/(sfgrd(ip0)-sfgrd(ip0-1))
      do i = 1 , nxdim
        datout(i) = coefint*(datin(i,ip0)-datin(i,ip0-1))+datin(i,ip0-1)
      end do
    end subroutine vinterp

    !------------------------------------------------------------------
    ! purpose: compute radiances (in mw/m2/str/cm-1) and derivatives of
    !   radiances with respect to atmospheric and surface parameters
    !------------------------------------------------------------------
    subroutine ossrad(ao,nn,xg,emrf,bs_sun,vn,n1,n2,sun,umu,umu0,   &
                      rad,xkt,xkemrf)
      implicit none
      type(aoss_ir) , intent(inout) :: ao
      logical , intent(in) :: sun
      integer , intent(in) :: n1 , n2 , nn
      real , intent(in) :: umu
      real , intent(in) :: umu0
      real , intent(in) :: bs_sun
      real , intent(in) :: vn
      real , dimension(:) , intent(in) :: xg
      real , dimension(:) , intent(in) :: emrf
      real , intent(out) :: rad
      real , intent(out) , dimension(:) :: xkt , xkemrf
      integer :: l , ks , k , ixoff
      real :: sec , secdif , sec0 , em , sun_rsfc , bs , dbs
      real :: tausfc , radsun , draddrsfc , draddtau_sun , draddemis
      real :: draddtskn , dtran , txsun , rsfc
      real :: bobs , dbobs
      real :: drdwobs2
      real :: draddtmpobs1 = 0.0E0
      real :: draddtmpobs2 = 0.0E0
      real :: draddtauobs1 = 0.0E0
      real :: draddtauobs2 = 0.0E0
      real :: sumtau0_up = 0.0E0
      real :: sumtau_dwn = 0.0E0
      sec = 1.0E0/umu
      sec0 = 1.0E0/umu0
      secdif = sec
      em = emrf(1)
      sun_rsfc = emrf(2)
      !-----------------------------------------------------------------------
      !     compute planck function and its derivative wrt temperature
      !-----------------------------------------------------------------------
      ! linintau is skipped in this version
      call planck_int(ao,vn,ao%tavl,xg(itsking),n2,ao%bbar,ao%dbbar,bs,dbs)
      !-----------------------------------------------------------------------
      !     compute transmittance profile along viewing path down to surface
      !-----------------------------------------------------------------------
      ao%txdn(1) = 1.0E0
      sumtau_dwn = 0.0E0
      if ( n1 > 1 ) then
        ao%txdn(1:n1-1) = 1.0E0
        sumtau_dwn = ao%tautotobs*sec
        ao%txdn(n1) = exp(-sumtau_dwn)
      end if
      do l = n1 , n2
        sumtau_dwn = sumtau_dwn+ao%tautot(l)*sec
        ao%txdn(l+1) = exp(-sumtau_dwn)
      end do
      tausfc = ao%txdn(n2+1)
      !---initialize radiance and derivative arrays
      rad = 0.0E0
      radsun = 0.0E0
      draddrsfc = 0.0E0
      draddtau_sun = 0.0E0
      ao%draddtmp(1:n2) = 0.0E0
      ao%draddtau(1:n2) = 0.0E0
      ao%draddtmpdw(1:n2+1)= 0.0E0
      ao%draddtmpuw(1:n2+1)= 0.0E0
      draddemis = 0.0E0
      draddtskn = 0.0E0
      !-----------------------------------------------------------------------
      !     1- downwelling thermal radiance calculation:
      !-----------------------------------------------------------------------
      if ( tausfc > 1.0E-06 ) then
        ao%txup(n2+1) = tausfc
        sumtau0_up = 0.0E0
        do l = n2 , 1 , -1
          sumtau0_up = sumtau0_up+ao%tautot(l)
          ao%txup(l) = exp(-(sumtau_dwn+sumtau0_up*secdif))
        end do
        do l = 1 , n2
          dtran = ao%txup(l+1)-ao%txup(l)
          !--- derivative of downwelling emission wrt to
          !--- temperature and constituents
          ao%draddtmp(l) = dtran*ao%dbbar(l)
          ao%draddtau(l) = (ao%txup(l)*ao%bbar(l)-rad)*secdif
          rad = rad+dtran*ao%bbar(l)
        end do
      end if
      !-----------------------------------------------------------------------
      !     2- add solar component
      !-----------------------------------------------------------------------
      if ( sun ) then
        txsun = exp(-(sumtau_dwn+sumtau0_up*sec0))
        !---derivatives wrt sfc solar reflectance
        draddrsfc = txsun*bs_sun*umu0/pi
        radsun = sun_rsfc*draddrsfc
        draddtau_sun = -radsun*(sec+sec0)
      end if
      !-----------------------------------------------------------------------
      !     3- add surface terms
      !-----------------------------------------------------------------------
      rsfc = (1.0E0-em)
      !---derivatives wrt emissivity and sfc skin temperature:
      draddemis = tausfc*bs-rad
      draddtskn = em*tausfc*dbs
      rad = rad*rsfc+em*tausfc*bs+radsun
      !-----------------------------------------------------------------------
      !     4- upwelling thermal radiance calculation
      !-----------------------------------------------------------------------
      do l = n2 , n1 , -1
        dtran = ao%txdn(l)-ao%txdn(l+1)
        ao%draddtau(l) = ao%draddtau(l)*rsfc + &
             (ao%txdn(l+1)*ao%bbar(l)-rad)*sec+draddtau_sun
        ao%draddtmp(l) = ao%draddtmp(l)*rsfc + &
              dtran*ao%dbbar(l)+ao%draddtau(l)*ao%dtaudtmp(l)
        rad = rad+dtran*ao%bbar(l)
      end do
      do l = 1 , n1-1
        ao%draddtau(l) = ao%draddtau(l)*rsfc
        ao%draddtmp(l) = ao%draddtmp(l)*rsfc
      end do
      if ( n1 > 1 ) then
        dtran = ao%txdn(n1-1)-ao%txdn(n1)
        call fplanck(vn,ao%tavlobs,bobs,dbobs)
        draddtauobs1 = ao%draddtau(l)*rsfc
        draddtauobs2 = (ao%txdn(n1)*bobs-rad)*sec+draddtau_sun
        draddtmpobs1 = ao%draddtmp(l)*rsfc+draddtauobs1*ao%dtaudtmp(n1-1)
        draddtmpobs2 = dtran*dbobs+draddtauobs2*ao%dtaudtmpobs
        rad = rad+dtran*bobs
      end if
      !-----------------------------------------------------------------------
      !     compute level derivatives and and map to array xkt:
      !-----------------------------------------------------------------------
      xkt = 0.0E0
      !  air temperature
      do l = n2 , 1 , -1
        xkt(itempg+l) = xkt(itempg+l)+ao%draddtmp(l)*ao%dtl(l)
        xkt(itempg+l-1) = xkt(itempg+l-1)+ao%draddtmp(l)*ao%dtu(l)
      end do
      if ( n1 > 1 ) then
        xkt(itempg+n1-1) = xkt(itempg+n1-1)+draddtmpobs2*ao%dtlobs
        xkt(itempg+n1-2) = xkt(itempg+n1-2)+draddtmpobs2*ao%dtuobs
      end if
      !---molecular concentrations
      do ks = 1 , ao%mxhmol
        k = ao%ah%imols(ks,nn)
        ixoff = ao%ah%molid(k)-1
        do l = n2 , 1 , -1
          ao%drdw(l) = ao%draddtau(l)*ao%abso(ks,l)
          xkt(ixoff+l+1) = xkt(ixoff+l+1)+ao%drdw(l)*ao%dwql(l,k)
          xkt(ixoff+l) = xkt(ixoff+l)+ao%drdw(l)*ao%dwqu(l,k)
        end do
        if ( n1 > 1 ) then
          ! drdwobs1 = draddtauobs1*ao%abso(ks,n1-1)
          drdwobs2 = draddtauobs2*ao%absoobs(ks)
          xkt(ixoff+n1) = xkt(ixoff+n1)+drdwobs2*ao%dwqlobs(ks)
          xkt(ixoff+n1-1) = xkt(ixoff+n1-1)+drdwobs2*ao%dwquobs(k)
        end if
      end do
      !---surface terms
      xkt(itsking) = draddtskn !--tskin
      xkemrf(1) = draddemis
      xkemrf(2) = draddrsfc
    end subroutine ossrad

    !------------------------------------------------------------------
    ! purpose: compute radiances (in mw/m2/str/cm-1) and derivatives of
    !   radiances with respect to atmospheric and surface parameters
    !------------------------------------------------------------------
    subroutine ossrad_vg(ao,nn,xg,emrf,bs_sun,vn,n1,n2, &
                         sun,umu,umu0,rad,xkt,xkemrf)
      implicit none
      type(aoss_ir) , intent(inout) :: ao
      logical , intent(in) :: sun
      integer , intent(in) :: n1 , n2 , nn
      real , intent(in) :: umu
      real , intent(in) :: umu0
      real , intent(in) :: bs_sun
      real , intent(in) :: vn
      real , dimension(:) , intent(in) :: xg
      real , dimension(:) , intent(in) :: emrf
      real , intent(out) :: rad
      real , intent(out) , dimension(:) :: xkt , xkemrf
      integer :: l , ks , k , ixoff
      real :: sec , secdif , sec0 , em , sun_rsfc , bs , dbs
      real :: tausfc , radsun , draddrsfc , draddtau_sun , draddemis
      real :: draddtskn , dtran , txsun , rsfc
      real :: sumtau0_up = 0.0E0
      real :: sumtau_dwn = 0.0E0
      sec = 1.0E0/umu
      sec0 = 1.0E0/umu0
      secdif = sec
      em = emrf(1)
      sun_rsfc = emrf(2)
      !-----------------------------------------------------------------------
      !     compute planck function and its derivative wrt temperature
      !-----------------------------------------------------------------------
      ! linintau is skipped in this version
      call planck_int(ao,vn,ao%tavl,xg(itsking),n2,ao%bbar,ao%dbbar,bs,dbs)
      !-----------------------------------------------------------------------
      !     compute transmittance profile along viewing path down to surface
      !-----------------------------------------------------------------------
      ao%txdn(1) = 1.0E0
      sumtau_dwn = 0.0E0
      if ( n1 > 1 ) then
        ao%txdn(1:n1-1) = 1.0E0
        sumtau_dwn = ao%tautotobs*sec
        ao%txdn(n1) = exp(-sumtau_dwn)
      end if
      do l = n1 , n2
        sumtau_dwn = sumtau_dwn+ao%tautot(l)*sec
        ao%txdn(l+1) = exp(-sumtau_dwn)
      end do
      tausfc = ao%txdn(n2+1)
      !---initialize radiance and derivative arrays
      rad = 0.0E0
      radsun = 0.0E0
      draddrsfc = 0.0E0
      draddtau_sun = 0.0E0
      ao%draddtmp(1:n2) = 0.0E0
      ao%draddtau(1:n2) = 0.0E0
      ao%draddtmpdw(1:n2+1)= 0.0E0
      ao%draddtmpuw(1:n2+1)= 0.0E0
      draddemis = 0.0E0
      draddtskn = 0.0E0
      !-----------------------------------------------------------------------
      !     1- downwelling thermal radiance calculation:
      !-----------------------------------------------------------------------
      if ( tausfc > 1.0E-06 ) then
        ao%txup(n2+1) = tausfc
        sumtau0_up = 0.0E0
        do l = n2 , 1 , -1
          sumtau0_up = sumtau0_up+ao%tautot(l)
          ao%txup(l) = exp(-(sumtau_dwn+sumtau0_up*secdif))
        end do
        do l = 1 , n2
          dtran = ao%txup(l+1)-ao%txup(l)
          !--- derivative of downwelling emission wrt to
          !--- temperature and constituents
          ao%draddtmp(l) = dtran*ao%dbbar(l)
          ao%draddtau(l) = (ao%txup(l)*ao%bbar(l)-rad)*secdif
          rad = rad+dtran*ao%bbar(l)
        end do
      end if
      !-----------------------------------------------------------------------
      !     2- add solar component
      !-----------------------------------------------------------------------
      if ( sun ) then
        txsun = exp(-(sumtau_dwn+sumtau0_up*sec0))
        !---derivatives wrt sfc solar reflectance
        draddrsfc = txsun*bs_sun*umu0/pi
        radsun = sun_rsfc*draddrsfc
        draddtau_sun = -radsun*(sec+sec0)
      end if
      !-----------------------------------------------------------------------
      !     3- add surface terms
      !-----------------------------------------------------------------------
      rsfc = (1.0E0-em)
      !---derivatives wrt emissivity and sfc skin temperature:
      draddemis = tausfc*bs-rad
      draddtskn = em*tausfc*dbs
      rad = rad*rsfc+em*tausfc*bs+radsun
      !-----------------------------------------------------------------------
      !     4- upwelling thermal radiance calculation
      !-----------------------------------------------------------------------
      do l = n2 , n1 , -1
        dtran = ao%txdn(l)-ao%txdn(l+1)
        ao%draddtau(l) = ao%draddtau(l)*rsfc + &
                 (ao%txdn(l+1)*ao%bbar(l)-rad)*sec+draddtau_sun
        ao%draddtmp(l) = ao%draddtmp(l)*rsfc + &
          dtran*ao%dbbar(l)+ao%draddtau(l)*ao%dtaudtmp(l)
        rad = rad+dtran*ao%bbar(l)
      end do
      do l = 1 , n1-1
        ao%draddtau(l) = ao%draddtau(l)*rsfc
        ao%draddtmp(l) = ao%draddtmp(l)*rsfc
      end do
      !-----------------------------------------------------------------------
      !     compute level derivatives and and map to array xkt:
      !-----------------------------------------------------------------------
      xkt = 0.0E0
      !  air temperature
      do l = n2 , 1 , -1
        xkt(itempg+l) = xkt(itempg+l)+ao%draddtmp(l)*ao%dtl(l)
        xkt(itempg+l-1) = xkt(itempg+l-1)+ao%draddtmp(l)*ao%dtu(l)
      end do
      !---molecular concentrations
      do ks = 1 , ao%ah%nmols(nn)
        k = ao%ah%imols(ks,nn)
        ixoff = ao%ah%molid(k)-1
        do l = n2 , 1 , -1
          ao%drdw(l) = ao%draddtau(l)*ao%abso(ks,l)
          xkt(ixoff+l+1) = xkt(ixoff+l+1)+ao%drdw(l)*ao%dwql(l,k)
          xkt(ixoff+l) = xkt(ixoff+l)+ao%drdw(l)*ao%dwqu(l,k)
        end do
      end do
      !---surface terms
      xkt(itsking) = draddtskn !--tskin
      xkemrf(1) = draddemis
      xkemrf(2) = draddrsfc
    end subroutine ossrad_vg

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
