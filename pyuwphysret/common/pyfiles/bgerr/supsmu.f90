module supsmu_mod

  private

  real(8) :: big , eps , sml
  real(8) , dimension(3) :: spans

  data spans , big , sml , eps /0.05 , 0.2 , 0.5 , 1.0E20 , 1.0E-7 , 1.0E-3/
  public :: supsmu , smooth

  contains

  subroutine supsmu(n,x,y,w,iper,span,alpha,smo)
    implicit none
    real(8) :: alpha , span
    integer :: iper , n
    real(8) , dimension(n) :: smo , w , x , y
    intent (in) alpha , iper
    real(8) , dimension(n) :: h
    real(8) , dimension(n,7) :: sc
    real(8) :: a , f , resmin , xscale , sw , sy , vsmlsq
    integer :: i , j , jper
!
!------------------------------------------------------------------
!
!     super-smoother.
!
!     Friedman J.H. (1984). A variable span smoother. Department of
!     Statistics, Stanford University, Technical Report LCS5.
!
!     version 10/10/84.
!
!     coded  and copyright (c) 1984 by:
!
!     Jerome H. Friedman
!     Department of Statistics
!     and
!     Stanford Linear Accelerator Center
!     Stanford University
!
!     all rights reserved.
!
!
!     input:
!     n : number of observations (x,y - pairs).
!     x(n) : ordered abscissa values.
!     y(n) : corresponding ordinate (response) values.
!     w(n) : weight for each (x,y) observation.
!     iper : periodic variable flag.
!     iper=1 => x is ordered interval variable.
!     iper=2 => x is a periodic variable with values
!     in the range (0.0,1.0) and peroid 1.0.
!     span : smoother span (fraction of observations in window).
!     span=0.0 => automatic (variable) span selection.
!     alpha : controles high frequency (small span) penality
!     used with automatic span selection (bass tone control).
!     (alpha.le.0.0 or alpha.gt.10.0 => no effect.)
!     output:
!     smo(n) : smoothed ordinate (response) values.
!     scratch:
!     sc(n,7) : internal working storage.
!
!     note:
!     for small samples (n < 40) or if there are substantial serial
!     correlations between obserations close in x - value, then
!     a prespecified fixed span smoother (span > 0) should be
!     used. reasonable span values are 0.2 to 0.4.
!
!------------------------------------------------------------------
!
    if ( x(n)>x(1) ) then
      i = n/4
      j = 3*i
      xscale = x(j) - x(i)
 50     continue
      if ( xscale>0.0D0 ) then
        vsmlsq = (eps*xscale)**2
        jper = iper
        if ( iper==2 .and. (x(1)<0.0D0 .or. x(n)>1.0D0) ) jper = 1
        if ( jper<1 .or. jper>2 ) jper = 1
        if ( span<=0.0D0 ) then
          do i = 1 , 3
            call smooth(n,x,y,w,spans(i),jper,vsmlsq,sc(1,2*i-1),sc(1,7))
            call smooth(n,x,sc(1,7),w,spans(2),-jper,vsmlsq,sc(1,2*i),h)
          end do
          do j = 1 , n
            resmin = big
            do i = 1 , 3
              if ( sc(j,2*i)<resmin ) then
                resmin = sc(j,2*i)
                sc(j,7) = spans(i)
              end if
            end do
            if ( alpha>0.0D0 .and. alpha<=10.0D0 .and. resmin<sc(j,6)     &
               & .and. resmin>0.0D0 ) sc(j,7) = sc(j,7)                 &
               & + (spans(3)-sc(j,7))*max(sml,resmin/sc(j,6))       &
               & **(10.0D0-alpha)
          end do
          call smooth(n,x,sc(1,7),w,spans(2),-jper,vsmlsq,sc(1,2),h)
          do j = 1 , n
            if ( sc(j,2)<=spans(1) ) sc(j,2) = spans(1)
            if ( sc(j,2)>=spans(3) ) sc(j,2) = spans(3)
            f = sc(j,2) - spans(2)
            if ( f>=0.0D0 ) then
              f = f/(spans(3)-spans(2))
              sc(j,4) = (1.0D0-f)*sc(j,3) + f*sc(j,5)
            else
              f = -f/(spans(2)-spans(1))
              sc(j,4) = (1.0D0-f)*sc(j,3) + f*sc(j,1)
            end if
          end do
          call smooth(n,x,sc(1,4),w,spans(1),-jper,vsmlsq,smo,h)
          return
        end if
      else
        if ( j<n ) j = j + 1
        if ( i>1 ) i = i - 1
        xscale = x(j) - x(i)
        go to 50
      end if
    else
      sy = 0.0D0
      sw = sy
      do j = 1 , n
        sy = sy + w(j)*y(j)
        sw = sw + w(j)
      end do
      a = 0.0D0
      if ( sw>0.0D0 ) a = sy/sw
      do j = 1 , n
        smo(j) = a
      end do
      return
    end if
    call smooth(n,x,y,w,span,jper,vsmlsq,smo,sc)
  end subroutine supsmu

  subroutine smooth(n,x,y,w,span,iper,vsmlsq,smo,acvr)
    implicit none
    integer :: iper , n
    real(8) :: span , vsmlsq
    real(8) , dimension(n) :: acvr , smo , w , x , y
    intent (in) iper , n , span , vsmlsq , w , x , y
    intent (inout) acvr , smo
    real(8) :: a , cvar , fbo , fbw , h , sy , tmp , var , wt ,  &
               xm , ym
    integer :: i , ibw , pin , it , j , j0 , jper , pout
    real(8) :: xti , xto
    xm = 0.0D0
    ym = xm
    var = ym
    cvar = var
    fbw = cvar
    jper = iabs(iper)
    ibw = int(0.5D0*span*n + 0.5D0)
    if ( ibw<2 ) ibw = 2
    it = 2*ibw + 1
    do i = 1 , it
      j = i
      if ( jper==2 ) j = i - ibw - 1
      xti = x(j)
      if ( j<1 ) then
        j = n + j
        xti = x(j) - 1.0D0
      end if
      wt = w(j)
      fbo = fbw
      fbw = fbw + wt
      if ( fbw>0.0D0 ) xm = (fbo*xm+wt*xti)/fbw
      if ( fbw>0.0D0 ) ym = (fbo*ym+wt*y(j))/fbw
      tmp = 0.0D0
      if ( fbo>0.0D0 ) tmp = fbw*wt*(xti-xm)/fbo
      var = var + tmp*(xti-xm)
      cvar = cvar + tmp*(y(j)-ym)
    end do
    do j = 1 , n
      pout = j - ibw - 1
      pin = j + ibw
      if ( .not.((jper/=2) .and. (pout<1 .or. pin>n)) ) then
        if ( pout<1 ) then
          pout = n + pout
          xto = x(pout) - 1.0D0
          xti = x(pin)
        else if ( pin<=n ) then
          xto = x(pout)
          xti = x(pin)
        else
          pin = pin - n
          xti = x(pin) + 1.0D0
          xto = x(pout)
        end if
        wt = w(pout)
        fbo = fbw
        fbw = fbw - wt
        tmp = 0.0D0
        if ( fbw>0.0D0 ) tmp = fbo*wt*(xto-xm)/fbw
        var = var - tmp*(xto-xm)
        cvar = cvar - tmp*(y(pout)-ym)
        if ( fbw>0.0D0 ) xm = (fbo*xm-wt*xto)/fbw
        if ( fbw>0.0D0 ) ym = (fbo*ym-wt*y(pout))/fbw
        wt = w(pin)
        fbo = fbw
        fbw = fbw + wt
        if ( fbw>0.0D0 ) xm = (fbo*xm+wt*xti)/fbw
        if ( fbw>0.0D0 ) ym = (fbo*ym+wt*y(pin))/fbw
        tmp = 0.0D0
        if ( fbo>0.0D0 ) tmp = fbw*wt*(xti-xm)/fbo
        var = var + tmp*(xti-xm)
        cvar = cvar + tmp*(y(pin)-ym)
      end if
      a = 0.0D0
      if ( var>vsmlsq ) a = cvar/var
      smo(j) = a*(x(j)-xm) + ym
      if ( iper>0 ) then
        h = 0.0D0
        if ( fbw>0.0D0 ) h = 1.0D0/fbw
        if ( var>vsmlsq ) h = h + (x(j)-xm)**2/var
        acvr(j) = 0.0D0
        a = 1.0D0 - w(j)*h
        if ( a>0.0D0 ) then
          acvr(j) = abs(y(j)-smo(j))/a
        else if ( j>1 ) then
          acvr(j) = acvr(j-1)
        end if
      end if
    end do
    j = 1
 100  continue
    j0 = j
    sy = smo(j)*w(j)
    fbw = w(j)
    if ( j<n ) then
 150  continue
      if ( x(j+1)<=x(j) ) then
        j = j + 1
        sy = sy + w(j)*smo(j)
        fbw = fbw + w(j)
        if ( j<n ) go to 150
      end if
    end if
    if ( j>j0 ) then
      a = 0.0D0
      if ( fbw>0.0D0 ) a = sy/fbw
      do i = j0 , j
        smo(i) = a
      end do
    end if
    j = j + 1
    if ( j<=n ) go to 100
  end subroutine smooth

end module supsmu_mod
