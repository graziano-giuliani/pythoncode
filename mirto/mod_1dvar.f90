module mod_1DVar

  implicit none

  private

  public :: init_1DVar
  public :: compute_chi_square
  public :: update_solution_1DVar

  real , allocatable , dimension(:,:) :: fmKT
  real , allocatable , dimension(:,:) :: KtSeInv
  real , allocatable , dimension(:,:) :: KtSeInvK
  real , allocatable , dimension(:,:) :: A
  real , allocatable , dimension(:,:) :: L
  real , allocatable , dimension(:,:) :: U
  real , allocatable , dimension(:) :: dx
  real , allocatable , dimension(:) :: dot1
  real , allocatable , dimension(:) :: dot2
  real , allocatable , dimension(:) :: dot3
  real , allocatable , dimension(:) :: dot4
  real , allocatable , dimension(:) :: dot5
  real , allocatable , dimension(:) :: d
  real , allocatable , dimension(:) :: x
  real , allocatable , dimension(:) :: y
  real , allocatable , dimension(:) :: r
  real , allocatable , dimension(:) :: dz
  real , allocatable , dimension(:) :: ddx
  real , allocatable , dimension(:) :: totx
  logical , allocatable , dimension(:,:) :: utmask
  integer , allocatable , dimension(:) :: ipivot

  logical :: isinit = .false.

  real , parameter :: xgamma = 3.0

  interface init_1DVar
    module procedure init_1DVar_2arg
    module procedure init_1DVar_arr
  end interface init_1DVar

  interface dotprod
    module procedure dotprod_matrix_matrix
    module procedure dotprod_matrix_array
    module procedure dotprod_array_matrix
    module procedure dotprod_array_array
  end interface

  interface
    subroutine sgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
      character(len=1) , intent(in) :: transa , transb
      integer , intent(in) :: m , n , k , lda , ldb , ldc
      real , intent(in) :: alpha , beta
      real , intent(in) , dimension(m,k) :: a
      real , intent(in) , dimension(k,n) :: b
      real , intent(out) , dimension(m,n) :: c
    end subroutine sgemm
  end interface

  interface
    subroutine sgemv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
      implicit none
      character(len=1) , intent(in) :: trans
      integer , intent(in) :: m , n , lda , incx , incy
      real , intent(in) :: alpha , beta
      real , intent(in) , dimension(m,n) :: a
      real , intent(in) , dimension(n) :: x
      real , intent(out) , dimension(m) :: y
    end subroutine sgemv
  end interface

  interface
    subroutine sgetrf(m,n,a,lda,ipivot,info)
      implicit none
      integer , intent(in) :: m , n , lda
      integer , intent(inout) , dimension(min(m,n)) :: ipivot
      real , intent(inout) , dimension(lda,n) :: a
      integer , intent(out) :: info
    end subroutine sgetrf
  end interface

  interface
    subroutine sgesv(n,nrhs,a,lda,ipivot,b,ldb,info)
      implicit none
      integer , intent(in) :: lda , ldb , n , nrhs
      integer , intent(inout) , dimension(n) :: ipivot
      real , intent(in) , dimension(lda,n) :: a
      real , intent(inout) , dimension(ldb,nrhs) :: b
      integer , intent(out) :: info
    end subroutine sgesv
  end interface

  interface
    real function sdot(n,sx,incx,sy,incy)
      implicit none
      integer , intent(in) :: n , incx , incy
      real , intent(in) , dimension(n) :: sx
      real , intent(in) , dimension(n) :: sy
    end function sdot
  end interface

  contains

    subroutine dotprod_matrix_matrix(a,b,c)
      implicit none
      real , dimension(:,:) , intent(in) :: a
      real , dimension(:,:) , intent(in) :: b
      real , dimension(:,:) , intent(out) :: c
      integer :: dim1 , dim2 , dim3
      real , parameter :: alpha = 1.0
      real , parameter :: beta = 0.0
      dim1 = size(a,1)
      dim2 = size(a,2)
      dim3 = size(b,2)
      if ( dim1 /= size(c,1) ) then
        write(0,*) 'Dimension error in DOTPROD : A(dim1) /= C(dim1)'
        write(0,*) 'SHAPE A : ', shape(a)
        write(0,*) 'SHAPE B : ', shape(b)
        write(0,*) 'SHAPE C : ', shape(c)
        stop
      end if
      if ( dim2 /= size(b,1) ) then
        write(0,*) 'Dimension error in DOTPROD : A(dim2) /= B(dim1)'
        write(0,*) 'SHAPE A : ', shape(a)
        write(0,*) 'SHAPE B : ', shape(b)
        write(0,*) 'SHAPE C : ', shape(c)
        stop
      end if
      if ( dim3 /= size(c,2) ) then
        write(0,*) 'Dimension error in DOTPROD : B(dim2) /= C(dim(2)'
        write(0,*) 'SHAPE A : ', shape(a)
        write(0,*) 'SHAPE B : ', shape(b)
        write(0,*) 'SHAPE C : ', shape(c)
        stop
      end if
      call sgemm('N','N',dim1,dim3,dim2,alpha,a,dim1,b,dim2,beta,c,dim1)
    end subroutine dotprod_matrix_matrix

    subroutine dotprod_matrix_array(a,b,c)
      implicit none
      real , dimension(:,:) , intent(in) :: a
      real , dimension(:) , intent(in) :: b
      real , dimension(:) , intent(out) :: c
      integer :: dim1 , dim2
      real , parameter :: alpha = 1.0
      real , parameter :: beta = 0.0
      dim1 = size(a,1)
      dim2 = size(a,2)
      if ( dim1 /= size(c) ) then
        write(0,*) 'Dimension error in DOTPROD : A(dim1) /= size(C)'
        write(0,*) 'SHAPE A : ', shape(a)
        write(0,*) 'SIZE  B : ', size(b)
        write(0,*) 'SIZE  C : ', size(c)
        stop
      end if
      if ( dim2 /= size(b) ) then
        write(0,*) 'Dimension error in DOTPROD : A(dim2) /= B(dim1)'
        write(0,*) 'SHAPE A : ', shape(a)
        write(0,*) 'SIZE  B : ', size(b)
        write(0,*) 'SIZE  C : ', size(c)
        stop
      end if
      call sgemv('N',dim1,dim2,alpha,a,dim1,b,1,beta,c,1)
    end subroutine dotprod_matrix_array

    subroutine dotprod_array_matrix(a,b,c)
      implicit none
      real , dimension(:) , intent(in) :: a
      real , dimension(:,:) , intent(in) :: b
      real , dimension(:) , intent(out) :: c
      real , dimension(1,size(a)) :: tmp
      real , dimension(1,size(b,2)) :: tmp1
      integer :: dim1 , dim2 , dim3
      real , parameter :: alpha = 1.0
      real , parameter :: beta = 0.0
      tmp(1,:) = a
      dim1 = size(tmp,1)
      dim2 = size(tmp,2)
      dim3 = size(b,2)
      if ( dim1 /= size(tmp1,1) ) then
        write(0,*) 'Dimension error in DOTPROD : A(dim1) /= C(dim1)'
        write(0,*) 'SHAPE A : ', shape(tmp)
        write(0,*) 'SHAPE B : ', shape(b)
        write(0,*) 'SHAPE C : ', shape(tmp1)
        stop
      end if
      if ( dim2 /= size(b,1) ) then
        write(0,*) 'Dimension error in DOTPROD : A(dim2) /= B(dim1)'
        write(0,*) 'SHAPE A : ', shape(tmp)
        write(0,*) 'SHAPE B : ', shape(b)
        write(0,*) 'SHAPE C : ', shape(tmp1)
        stop
      end if
      if ( dim3 /= size(tmp1,2) ) then
        write(0,*) 'Dimension error in DOTPROD : B(dim2) /= C(dim(2)'
        write(0,*) 'SHAPE A : ', shape(tmp)
        write(0,*) 'SHAPE B : ', shape(b)
        write(0,*) 'SHAPE C : ', shape(tmp1)
        stop
      end if
      call sgemm('N','N',dim1,dim3,dim2,alpha,tmp,dim1,b,dim2,beta,tmp1,dim1)
      c = tmp1(:,1)
    end subroutine dotprod_array_matrix

    subroutine dotprod_array_array(a,b,d)
      implicit none
      real , intent(in) , dimension(:) :: a , b
      real , intent(out) :: d
      integer :: dim1
      dim1 = size(a)
      if ( dim1 /= size(b) ) then
        write(0,*) 'Dimension error in DOTPROD : A(dim1) /= B(dim1)'
        write(0,*) 'SIZE  A : ', size(a)
        write(0,*) 'SIZE  B : ', size(b)
        stop
      end if
      d = sdot(size(a),a,1,b,1)
    end subroutine dotprod_array_array

    subroutine ludecomp(a,l,u)
      implicit none
      real , intent(in) , dimension(:,:) :: a
      real , intent(out) , dimension(:,:) :: l
      real , intent(out) , dimension(:,:) :: u
      real , dimension(size(a,1),size(a,2)) :: lu
      integer :: i , dim1 , dim2 , info
      dim1 = size(a,1)
      dim2 = size(a,2)
      if ( dim1 /= size(l,1) .or. &
           dim2 /= size(l,2) .or. &
           dim1 /= size(u,1) .or. &
           dim2 /= size(u,2) ) then
        write(0,*) 'Dimension error in ludecomp'
        write(0,*) 'SHAPE A : ', shape(a)
        write(0,*) 'SHAPE L : ', shape(l)
        write(0,*) 'SHAPE U : ', shape(u)
        stop
      end if
      lu = a
      call sgetrf(dim1,dim2,lu,dim1,ipivot,info)
      if ( info /= 0 ) then
        if ( info < 0 ) then
          write(0,*) 'INVALID ARGUMENT NO ',-info
        else
          write(0,*) 'Singular matrix !!!'
          write(0,*) 'U(',info,',',info,') = ',lu(info,info)
        end if
        stop
      end if
      l = 0.0
      u = 0.0
      where (utmask)
        l = lu
      else where
        u = lu
      end where
      do i = 1 , dim2
        l(i,i) = 1.0
      end do
    end subroutine ludecomp

    subroutine solve(A,b,x)
      implicit none
      real , intent(in) , dimension(:,:) :: a
      real , intent(in) , dimension(:) :: b
      real , intent(out) , dimension(:) :: x
      real , dimension(size(a,1),size(a,2)) :: lu
      real , dimension(size(a,1),1) :: xb
      integer :: dim1 , dim2 , info
      dim1 = size(a,1)
      dim2 = size(a,2)
      lu = a
      xb(:,1) = b
      call sgesv(dim2,1,lu,dim1,ipivot,xb,dim1,info)
      x = xb(:,1)
    end subroutine solve

    subroutine init_1DVar_2arg(m,n)
      implicit none
      integer , intent(in) :: n , m
      integer :: i , j
      if ( isinit ) then
        if ( size(KtSeInv,1) /= n  .and. &
             size(KtSeInv,2) /= m ) then
          deallocate(KtSeInv)
          deallocate(KtSeInvK)
          deallocate(fmKT)
          deallocate(A)
          deallocate(utmask)
          deallocate(L)
          deallocate(U)
          deallocate(ipivot)
          deallocate(dx)
          deallocate(d)
          deallocate(x)
          deallocate(y)
          deallocate(r)
          deallocate(dz)
          deallocate(ddx)
          deallocate(totx)
          deallocate(dot1)
          deallocate(dot2)
          deallocate(dot3)
          deallocate(dot4)
          deallocate(dot5)
        else
          return
        end if
      end if
      allocate(KtSeInv(n,m))
      allocate(KtSeInvK(n,n))
      allocate(A(n,n))
      allocate(utmask(n,n))
      allocate(L(n,n))
      allocate(U(n,n))
      allocate(ipivot(n))
      allocate(dx(n))
      allocate(fmKT(n,m))
      allocate(dot1(n))
      allocate(dot2(n))
      allocate(dot3(n))
      allocate(d(n))
      allocate(x(n))
      allocate(y(n))
      allocate(r(n))
      allocate(dz(n))
      allocate(ddx(n))
      allocate(totx(n))
      allocate(dot4(m))
      allocate(dot5(m))
      do j = 1 , n
        do i = 1 , n
          utmask(i,j) = (j-i) > 0
        end do
      end do
      isinit = .true.
    end subroutine init_1DVar_2arg

    subroutine init_1DVar_arr(sh)
      implicit none
      integer , intent(in) , dimension(2) :: sh
      call init_1DVar_2arg(sh(1),sh(2))
    end subroutine init_1DVar_arr

    real function compute_chi_square(m,SEInv,yobs_minus_yhat) result(norm)
      implicit none
      integer , intent(in) :: m
      real , dimension(m,m) , intent(in) :: SEInv
      real , dimension(m) , intent(in) :: yobs_minus_yhat
      if ( .not. isinit ) return
      call dotprod(yobs_minus_yhat,SEInv,dot4)
      dot5 = yobs_minus_yhat/size(yobs_minus_yhat)
      call dotprod(dot4,dot5,norm)
    end function compute_chi_square

    subroutine update_solution_1DVar(m,n,fmK,SEInv,SaInv_ret, &
                                     yobs_minus_yhat,xhat,ap,xhat_new,d2)
      implicit none
      integer , intent(in) :: m , n
      real , dimension(m,n) , intent(in) :: fmK
      real , dimension(m,m) , intent(in) :: SEInv
      real , dimension(n,n) , intent(in) :: SaInv_ret
      real , dimension(m) , intent(in) :: yobs_minus_yhat
      real , dimension(n) , intent(in) :: xhat
      real , dimension(n) , intent(in) :: ap
      real , dimension(n) , intent(out) :: xhat_new
      real , intent(out) :: d2
      if ( .not. isinit ) return
      fmKT = transpose(fmK)
      call dotprod(fmKT,SEInv,KtSeInv)
      call dotprod(KtSeInv,fmK,KtSeInvK)
      A = KtSeInvK+(1.0+xgamma*SaInv_ret)
      dx = xhat - ap
      call dotprod(KtSeInv,yobs_minus_yhat,dot1)
      call dotprod(SaInv_ret,dx,dot2)
      d = dot1 - dot2
      call ludecomp(A,L,U)
      call solve(L,d,y)
      call solve(U,y,x)
      call dotprod(A,x,dot3)
      r = d - dot3
      call solve(L,r,dz)
      call solve(U,dz,ddx)
      totx = x+ddx
      xhat_new = xhat+totx
      call dotprod(totx,d,d2)
    end subroutine update_solution_1DVar

end module mod_1DVar
