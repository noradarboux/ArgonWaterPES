!       program poth2oar
        module  Pot_h2oar
! Computes the 3D H2O-Ar potential in CM coordinates (R,theta,phi), with H2O rigid,
! as described in Section II of the paper.
! Distances should be given in bohrs, angles in degrees. The potential is returned in cm-1.
! The program requires the file v_lm.txt that contains the expansion coefficients
! v_lm(R) depending on the intermolecular distance R.
! 
! Test input: R=6.92 a0, theta=90 deg, phi=110 deg.
! Test output: V(R,theta,phi)=-139.4397 cm-1.

!     implicit none
!     double precision :: R,theta,phi,pot,potential,pot_h2oar
!       double precision :: pi
!       pi=dacos(-1.d0)
!       write(*,*) "enter R (a.u.), theta (deg), phi(deg)"
!       read(*,*) R, theta, phi
!       pot_h2oar=potential(R,theta,phi) !/219474.63137d0  ! potential in au
!       print*, R,theta*180d0/pi,phi*180d0/pi,pot_h2oar

!       end program

       contains
       function potHartree(sph)
       implicit none
       !input: spherical coordinate (theta, phi, r)
       double precision, dimension(3) :: sph
       !output: potential energy in hartree
       double precision               :: potHartree

       potHartree = potential(sph(3),sph(1),sph(2))/219474.63
       end function
       
       function potential(R,theta,phi)
       implicit none

       integer, parameter :: nr=22 ! number of internuclear distances R
       integer, parameter :: n1=11 ! number of values of l

       integer :: i,j,k,k2,l,m,il,im

       double precision :: potential
       double precision :: yp1,ypn
       double precision :: v_lm(n1,n1,nr)
       double precision :: pot(nr),pot2(nr)
       double precision :: rpot(nr)
       double precision :: x1,y,pi,energy
       double precision :: tesseral1,tes(n1*n1),theta,phi,R
       logical :: first
       data first/.true./
       save 
       pi=dacos(-1.d0)
       theta=theta*pi/180d0
       phi=phi*pi/180d0

      if (first) then
      open(unit=10,file='v_lm.txt',status='old')
      do l=0,n1-1
        do m=0,l
            do k=1,nr
              read(10,*) rpot(k),v_lm(l+1,m+1,k)
            enddo
!          pot(k,k2)=v_lm(l+1,m+1,k,k2)
        enddo
      enddo
      close(10)
      first=.false.
      endif

      call tesser(tes,n1-1,theta,phi)   ! tesseral harmonics with l=0...10
C       The tesseral t(l,m,theta,phi) is accessed via tes(l*(l+1)+m+1)

      x1=R
      energy=0d0
      do l=0,n1-1
        do m=0,l
          tesseral1=tes(l*(l+1)+m+1)
            pot(:)=v_lm(l+1,m+1,:)
        yp1=(pot(2)-pot(1))/(rpot(2)-rpot(1))
        ypn=(pot(nr)-pot(nr-1))/(rpot(nr)-rpot(nr-1))
        call spline(rpot,pot,nr,yp1,ypn,pot2)
        call splint(rpot,pot,pot2,nr,x1,y)
            energy=energy+y*tesseral1
        enddo
      enddo
      potential=energy

111   format(25(E14.6))
112   format(17(E16.8,1h,))

       end function
************************************************************************


C      ONE-DIMENSIONAL SPLINE
c
c======================================================================
c
      SUBROUTINE spline(x,y,n,yp1,ypn,y2)
      INTEGER n,NMAX
      DOUBLE PRECISION yp1,ypn,x(n),y(n),y2(n)
      PARAMETER (NMAX=500)
      INTEGER i,k
      DOUBLE PRECISION p,qn,sig,un,u(NMAX)
      if (yp1.gt..99e30) then
        y2(1)=0.
        u(1)=0.
      else
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+
     *1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*
     *u(i-1))/p
11    continue
      if (ypn.gt..99e30) then
        qn=0.
        un=0.
      else
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      return
      END


      SUBROUTINE splint(xa,ya,y2a,n,x,y)
      INTEGER n
      DOUBLE PRECISION x,y,xa(n),y2a(n),ya(n)
      INTEGER k,khi,klo
      DOUBLE PRECISION a,b,h
      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
!      if (h.eq.0.) pause 'bad xa input in splint'
      if (h.eq.0.) then
        print*, 'bad xa input in splint'
        stop
      endif
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**
     *2)/6.
      return
      END

c
c======================================================================
c


      SUBROUTINE TESSER(tes,n,theta,phi)
C
C     Subroutine computes tesseral harmonics as defined by J.L. Prather
C     [N.B.S. monograph 19 (1961)].
C
C     The tesserals are computed from l=m=0 upwards to l=m=n for the
C     polar angles theta and phi (radians).
C     The results are stored in the linear array tes in increasing order
C     of l and m.
C     Negative m-values refer to prather's s(l,m) (= p(l,m)*sin(m*phi))
C     positive m-values to prather's c(l,m) ( = p(l,m) * cos(m*phi) )
C     Note:
C         the tesseral t(l,m,theta,phi) is accessed via tes(l*(l+1)+m+1)
C
      implicit double precision (a-h,o-z)
      parameter ( nn = 50 )
      dimension tes( (n+1)*(n+1) )
      dimension p( (nn+1)*(nn+2)/2 ),sn(nn+1),cs(nn+1)
      data tpih/3.989422804014327d-1/
C     (   tpih = 1/dsqrt(2*pi)   )
      data twmh/7.071067811865475d-1/
C     (   twmh = 1./dsqrt(2)    )
C
      if (n .gt. nn) then
         write(*,'(''0n larger than:'',i3,'' increase nn in'',
     1              '' subroutine tesser '' )' ) nn
         stop 16
      endif
C
C     Compute associated Legendre functions
C
      cost = dcos(theta)
      Call assleg(p,n,cost)
C
      tes(1) = tpih * twmh * p(1)
      if ( n .eq. 0) return
C
C     Compute necessary sines and cosines
C
      cs(1) = 1.d0
      sn(1) = 0.d0
      cosp = dcos(phi)
      sinp = dsin(phi)
      do 10 m=1,n
      cs(m+1) = cs(m)*cosp - sn(m)*sinp
      sn(m+1) = sn(m)*cosp + cs(m)*sinp
   10 continue
C
      kp = 1
      ll = 1
      dl = 0.d0
C
      do 40 l=1,n
         ll = ll + l + l
C        ( ll = l(l+1) + 1   )
         kp = kp + 1
         dl = dl + 1.d0
         fac = dsqrt(dl+dl+1.d0) * tpih
         tes(ll) = fac*twmh * p(kp)
         kt1 = ll
         kt2 = ll
         dlm = dl
         dlm1 = dl + 1.d0
C
         do 30 m=1,l
            kt1 = kt1 + 1
            kt2 = kt2 - 1
C           ( kt1 = l(l+1) +1 + m  )
C           ( kt2 = l(l+1) +1 - m )
            kp = kp + 1
            dlm  = dlm  + 1.d0
            dlm1 = dlm1 - 1.d0
C           ( dlm = l+m      )
C           ( dlm1= l+1-m    )
            fac = fac / dsqrt(dlm*dlm1)
            tes(kt1) = fac * p(kp) * cs(m+1)
            tes(kt2) = fac * p(kp) * sn(m+1)
C           ( t(l,m) = fac * p(l,m) * cos(m*phi)  )
C           ( t(l,-m)= fac * p(l,m) * sin(m*phi)  )
C
   30    continue
   40 continue
C
      end

      SUBROUTINE ASSLEG(p,n,x)
C
C     The subroutine computes the associated Legendre polynomials as defined
C     by Edmonds (angular momentum in quantum mechanics).
C     x is the usual coordinate (-1 < x < +1 ), n is the maximum
C     quantum number. the associated legendre functions are computed
C     in the order (0,0),(1,0),(1,1),(2,0),(2,1),(2,2),.... ,(n,n)
C     and returned in the array p.
C     The associated Legendre function p(l,m,x) may be  accessed via
C     p( l*(l+1)/2 + m+1 )
C
      implicit double precision (a-h,o-z)
      dimension p((n+1)*(n+2)/2 )
C
      p(1) = 1.d0
      if (n .eq. 0) return
      sint = sqrt(1.d0 - x*x)
      p(2) = x
      p(3) = sint
      if (n .eq. 1) return
C
      lm1 = 1
      lm  = 3
      dl  = 1.d0
C
      do 20 l=2,n
      dl = dl + 1.d0
      lm1 = lm1 + 1
C     (   lm1 = l*(l-1)/2 + 1 )
      lm  = lm + 1
C     (   lm = l*(l+1)/2 + 1 )
C
      p(lm) = x*p(lm1) - sint*p(lm1+1)/dl
C     ( p(l,0) = x*p(l-1,0) - dsqrt(1-x*x)*p(l-1,1)/l   )
      mmax = l-1
      dlm = dl
      do 10 m=1,mmax
      lm1 = lm1 + 1
      lm  = lm + 1

      p(lm) = dlm*sint*p(lm1-1) + x*p(lm1)
C     (  p(l,m) = (l+m-1)*dsqrt(1-x*x)*p(l-1,m-1) + x*p(l-1,m)   )
      dlm = dlm + 1.d0
   10 continue
      lm = lm + 1
      p(lm) = (dl+dl-1.d0)*sint*p(lm1)
C     (  p(l,l) = (2*l-1)*dsqrt(1-x*x)*p(l-1,l-1)    )

   20 continue

      END

      end module
