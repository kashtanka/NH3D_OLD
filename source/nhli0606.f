      subroutine pos1(phib,n,dd,phib0,phib1,a1,a2,a3,c,d)
c
c 1d poisson solver with dirichlet b.c.
c
      implicit real*8(a-h,o-z)
      dimension phib(0:n+1),a1(n-1),a2(n-1),a3(n-1),c(n-1),d(n-1)
      save a21
c
      do 10 i=1,n-1
      a1(i)=1./dd
      a2(i)=-2./dd
      a3(i)=a1(i)
10    continue
c
      phib(2)=phib(2)-phib0/dd
      phib(n)=phib(n)-phib1/dd
c
      call tri1(phib,a1,a2,a3,n-1,c,d,a21)
c
      return
      end
c
c
c
      subroutine tri1(a4,a1,a2,a3,n,c,d,a21)
c-----------------------------------------------------------------------
c     to solve a system of linear tridigonal equations
c     using perchasing method.  the coefficients
c     of n equations are stored on arrays a1,a2,a3,a4
c     a1,a2,a3  -matrix of coefficients
c     a4        -r.h.s. of equation system
c     n         -number of eqations in the system
c     the solution vector is returned in a4.
c
c note the dimensions of a4!
c-----------------------------------------------------------------------
c
      implicit real*8(a-h,o-z)
      dimension c(n),d(n)
      dimension a1(n),a2(n),a3(n),a4(-1:n+1)
c
c     compute the recurrence coefficients
c
      c(1)=a3(1)/a2(1)
      a21=1./a2(1)
      do 10 j=2,n-1
      c(j)=a3(j)/(a2(j)-a1(j)*c(j-1))
10    continue
      c(n)=1./(a2(n)-a1(n)*c(n-1))
c
c coefficients denpendant on a a4 and solution:
c
      entry rtri1(a1,a3,a4,n,c,d,a21)
c
      d(1)=a4(1)*a21
      do 20 j=2,n-1
      d(j)=(a4(j)-a1(j)*d(j-1))*c(j)/a3(j)
20    continue
c
c     back substitution for  solution
c
      a4(n)=(a4(n)-a1(n)*d(n-1))*c(n)
      do 30 j=n-1,1,-1
      a4(j)=d(j)-c(j)*a4(j+1)
30    continue
      return
      end
c
c
c
      subroutine intri2(a1,a2,a3,nx,ny,c,d,a21)
c-----------------------------------------------------------------------
c     to solve a set of linear tridigonal equations
c     using perchasing method.  the coefficients
c     of n equations are stored on arrays a1,a2,a3,a4
c     a1,a2,a3  -matrix of coefficients
c     a4        -r.h.s. of equation system
c     ny   -number of eqations in each system
c     nx -number of equation sets to be solved in paralell
c     the solution vector is returned in a4.
c
c note the dimensions of a4!
c-----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      dimension c(nx,ny),d(nx,ny),a21(nx)
      dimension a1(ny),a2(nx,ny),a3(ny),a4(-1:nx+1,-1:ny+1)
c
c-----------------------------------------------------------------------
c     compute the recurrence coefficients
c-----------------------------------------------------------------------
      do 5 k=1,nx
      c(k,1)=a3(1)/a2(k,1)
      a21(k)=1./a2(k,1)
5     continue
      do 10 j=2,ny-1
      do 10 k=1,nx
      c(k,j)=a3(j)/(a2(k,j)-a1(j)*c(k,j-1))
10    continue
      do 11 k=1,nx
      c(k,ny)=1./(a2(k,ny)-a1(ny)*c(k,ny-1))
11    continue
      return
c
c coefficients denpendant on a a4 and solution:
c
      entry tri2(a1,a3,a4,nx,ny,c,d,a21)
c
      do 15 k=1,nx
      d(k,1)=a4(k,1)*a21(k)
15    continue
      do 20 j=2,ny-1
      do 20 k=1, nx
      d(k,j)=(a4(k,j)-a1(j)*d(k,j-1))*c(k,j)/a3(j)
20    continue
c-----------------------------------------------------------------------
c     back substitution for  solution
c-----------------------------------------------------------------------
      do 25 k=1,nx
      a4(k,ny)=(a4(k,ny)-a1(ny)*d(k,ny-1))*c(k,ny)
25    continue
      do 30 j=ny-1,1,-1
      do 30 k=1,nx
      a4(k,j)=d(k,j)-c(k,j)*a4(k,j+1)
30    continue
      return
      end
c
c
c
      subroutine iytri2(a1,a2,a3,nx,ny,c,d)
c-----------------------------------------------------------------------
c     to solve a set of linear tridigonal equations
c     using perchasing method.  the coefficients
c     of n equations are stored on arrays a1,a2,a3,a4
c     a1,a2,a3  -matrix of coefficients
c     a4        -r.h.s. of equation system
c     ny   -number of eqations in each system
c     nx -number of equation sets to be solved in paralell
c     the solution vector is returned in b.
c
c note the dimensions of b!
c-----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      dimension c(nx,ny),d(nx,ny)
      dimension a1(ny),a2(nx,ny),a3(ny)
      dimension b(-1:nx+1,-1:ny+1),x(-1:nx+1,-1:ny+1)
c
c compute the recurrence coefficients
c
      do 5 k=1,nx
      c(k,1)=a3(1)/a2(k,1)
5     continue
      do 10 j=2,ny-1
      do 10 k=1,nx
      c(k,j)=a3(j)/(a2(k,j)-a1(j)*c(k,j-1))
10    continue
      do 11 k=1,nx
      c(k,ny)=1./(a2(k,ny)-a1(ny)*c(k,ny-1))
11    continue
      return
c
c coefficients denpendant on a b and solution:
c
      entry ytri2(a1,a2,a3,b,x,nx,ny,c,d)
c
      do 15 k=1,nx
      d(k,1)=b(k,1)/a2(k,1)
15    continue
      do 20 j=2,ny-1
      do 20 k=1, nx
      d(k,j)=(b(k,j)-a1(j)*d(k,j-1))*c(k,j)/a3(j)
20    continue
c
c back substitution for  solution
c
      do 25 k=1,nx
      x(k,ny)=(b(k,ny)-a1(ny)*d(k,ny-1))*c(k,ny)
25    continue
      do 30 j=ny-1,1,-1
      do 30 k=1,nx
      x(k,j)=d(k,j)-c(k,j)*x(k,j+1)
30    continue
      return
      end
c
c ***
c
      subroutine intri3(a1,a2,a3,nx,ny,ns,c,a21)
c-----------------------------------------------------------------------
c     to solve nx*ny systems of ns linear tridiagonal equations
c     using perchasing method.  the coefficients
c     are stored on arrays a1,a2,a3,a4
c     a1,a2,a3  -matrix of coefficients
c     a4        -r.h.s. of equation system
c     ns    - number of equations in each system
c     nx*ny - number of equation sets to be solved in paralell
c     the solution vector is returned in a4.
c--------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      dimension c(nx,ny,ns),d(nx,ny,ns),a21(nx,ny)
      dimension a1(ns),a2(nx,ny,ns),a3(ns),a4(-1:nx+1,-1:ny+1,0:ns+2)
c
c     compute the recurrence coefficients
c
      do 5 iy=1,ny
      do 5 ix=1,nx
      c(ix,iy,1)=a3(1)/a2(ix,iy,1)
      a21(ix,iy)=1./a2(ix,iy,1)
5     continue
      do 10 is=2,ns-1
      do 10 iy=1,ny
      do 10 ix=1,nx
      c(ix,iy,is)=a3(is)/(a2(ix,iy,is)-a1(is)*c(ix,iy,is-1))
10    continue
      do 11 iy=1,ny
      do 11 ix=1,nx
      c(ix,iy,ns)=1./(a2(ix,iy,ns)-a1(ns)*c(ix,iy,ns-1))
11    continue
      return
c
c coefficients denpendant on a a4 and solution:
c
      entry tri3(a1,a3,a4,nx,ny,ns,c,d,a21)
c
      do 15 iy=1,ny
      do 15 ix=1,nx
      d(ix,iy,1)=a4(ix,iy,1)*a21(ix,iy)
15    continue


      do 20 is=2,ns-1
      do 20 iy=1,ny
      do 20 ix=1,nx
      d(ix,iy,is)=(a4(ix,iy,is)-a1(is)*d(ix,iy,is-1))*c(ix,iy,is)/a3(is)
20    continue
c
c back substitution for  solution
c
      do 25 iy=1,ny
      do 25 ix=1,nx
      a4(ix,iy,ns)=(a4(ix,iy,ns)-a1(ns)*d(ix,iy,ns-1))*c(ix,iy,ns)
25    continue
      do 30 is=ns-1,1,-1
      do 30 iy=1,ny
      do 30 ix=1,nx
      a4(ix,iy,is)=d(ix,iy,is)-c(ix,iy,is)*a4(ix,iy,is+1)
30    continue
      return
      end
c
c ***
c
      subroutine infft(ifax,ifay,trigsx,trigsy,sipcox,sipcoy,simcox
     :   ,simcoy,xlamb,ylamb,xylamb,nx,ny,dxx,dyy)
c----------------------------------------------------------------------
c preprocessing for fourier transforms and parameters for the poisson
c equation:
c----------------------------------------------------------------------
c
      implicit real*8(a-h,o-z)
      dimension ifax(*),ifay(*),sipcox(*),sipcoy(*),simcox(*),simcoy(*)
     :   ,xlamb(*),ylamb(*),xylamb(nx-1,ny-1)
      dimension trigsx(*),trigsy(*)
c
      double precision dpi
c      double precision cosinx(nx-1),cosiny(ny-1)
c
      dpi=4.d0*datan(1.d0)
c
c      call fax(ifax,nx-1,3)
c      call fax(ifay,ny-1,3)
c      call fftrig(trigsx,nx-1,3)
c      call fftrig(trigsy,ny-1,3)
c
      call set99(trigsx,ifax,nx-1)
      call set99(trigsy,ifay,ny-1)
c
      call insct(sipcox,simcox,nx-1)
      call insct(sipcoy,simcoy,ny-1)
c
      nxnx2=2*(nx-1)
      nyny2=2*(ny-1)
c
c      do 110 ix=1,nx-1
c      cosinx(ix)=dcos((ix-1)*dpi/dble(nx-1))
c110   continue
c      do 111 iy=1,ny-1
c      cosiny(iy)=dcos((iy-1)*dpi/dble(ny-1))
c111   continue
c
      do 120 ix=1,nx-1
      xlamb(ix)=4.d0*(dsin((ix-1)*dpi/dble(nxnx2)))**2/dxx
120   continue
      do 121 iy=1,ny-1
      ylamb(iy)=4.d0*(dsin((iy-1)*dpi/dble(nyny2)))**2/dyy
121   continue
c
c      write(6,*) 'xlamb'
c      write(6,1000) (xlamb(i),i=1,nx-1)
c      write(6,*) 'ylamb'
c      write(6,1000) (ylamb(i),i=1,ny-1)
c1000  format(3e22.15)
c
      do 130 iy=1,ny-1
      do 130 ix=1,nx-1
      xylamb(ix,iy)=xlamb(ix)+ylamb(iy)
130   continue
      return
      end
c
c
c
      subroutine insct(sipco,simco,n)
c
c preparation of coefficients for scfft
c
      implicit real*8(a-h,o-z)
      dimension sipco(n),simco(n)
      double precision dpi
c
      dpi=4.d0*datan(1.d0)
c
      do 20 k=1,n
      sipco(k)=dsin(dpi*dble(k-1)/(dble(2*n)))
     :        +dcos(dpi*dble(k-1)/(dble(2*n)))
      simco(k)=dsin(dpi*dble(k-1)/(dble(2*n)))
     :        -dcos(dpi*dble(k-1)/(dble(2*n)))
20    continue
c
c      write(6,*) 'sipco'
c      write(6,1000) (sipco(ix),ix=1,nx)
c      write(6,*) 'simco'
c      write(6,1000) (simco(ix),ix=1,nx)
c1000  format(1x,3e22.15)
c
      return
      end
c
c ***
c
      subroutine dscfft(f,wf,wx,wy,trigsx,ifax,sipcox,simcox,nx
     ;   ,trigsy,ifay,sipcoy,simcoy,ny,lot,isign)
c-----------------------------------------------------------------------
c double cosine transform on a staggered grid.
c
c for each is=1,ns performs the calculation of:
c
c    for i=1,nx;j=1,ny
c       sum(ix=1,nx;iy=1,ny) f(ix,iy,is)*cos((ix-0.5)*(i-1)*pi/nx)
c                            *cos((iy-0.5)*(j-1)*pi/ny)
c
c following wilhelmson and ericksen (j.c.p. 25,319-331) in the design
c of the pre and posprocessing (making a staggered cosine transform with
c a standard real transform).
c
c in the vectorization it has been assumed (in this routine, not in
c fft991) that nx>ny>ns.
c
c note the dimensions of f (they are appropriate for nh3d program)
c
c p. m. nov 1988
c-----------------------------------------------------------------------
c
      implicit real*8(a-h,o-z)
      dimension f(-1:nx+1,-1:ny+1,0:lot+2)
      dimension wf(-1:nx+1,-1:ny+1,0:lot+2)
c
      dimension wx(nx+2,ny,lot),wy(ny+2,nx,lot)
c
      dimension sipcox(nx),simcox(nx)
      dimension sipcoy(ny),simcoy(ny)
      dimension trigsx(*),trigsy(*)
      dimension ifax(10),ifay(10)
c      write(*,*) 'dscfft:',nx,ny,lot
c      write(*,'(''dscfft ifax:'',10i4)') (ifax(kkk),kkk=1,10)
c      write(*,'(''dscfft ifay:'',10i4)') (ifay(kkk),kkk=1,10)
c
      if(isign.eq.1) then
c
c x transform:
c
         do 10 iy=1,ny
         do 10 j=1,lot
         wx(1,iy,j)=f(1,iy,j)
         wx(2,iy,j)=0.
         wx(nx+1,iy,j)=f(nx,iy,j)
         wx(nx+2,iy,j)=0.
10       continue
         do 20 iy=1,ny
         do 20 j=1,lot
         do 20 ii=2,nx-2,2
         wx(ii+1,iy,j)=0.5*(f(ii,iy,j)+f(ii+1,iy,j))
         wx(ii+2,iy,j)=-0.5*(f(ii,iy,j)-f(ii+1,iy,j))
20       continue
c
         call fft991(wx,wf,trigsx,ifax,1,nx+2,nx,lot*ny,+1)
c
         do 30 iy=1,ny
         do 30 j=1,lot
         f(1,iy,j)=wx(1,iy,j)
30       continue
         do 40 iy=1,ny
         do 40 j=1,lot
         do 40 ix=2,nx
         f(ix,iy,j)=0.5*(sipcox(ix)*wx(ix,iy,j)-simcox(ix)
     :      *wx(nx-ix+2,iy,j))
40       continue
c
c y transform:
c
         do 50 j=1,lot
         do 50 ix=1,nx
         wy(1,ix,j)=f(ix,1,j)
         wy(2,ix,j)=0.
         wy(ny+1,ix,j)=f(ix,ny,j)
         wy(ny+2,ix,j)=0.
50       continue
         do 60 ii=2,ny-2,2
         do 60 j=1,lot
         do 60 ix=1,nx
         wy(ii+1,ix,j)=0.5*(f(ix,ii,j)+f(ix,ii+1,j))
         wy(ii+2,ix,j)=-0.5*(f(ix,ii,j)-f(ix,ii+1,j))
60       continue
c
         call fft991(wy,wf,trigsy,ifay,1,ny+2,ny,lot*nx,+1)
c
         do 70 j=1,lot
         do 70 ix=1,nx
         f(ix,1,j)=wy(1,ix,j)
70       continue
         do 80 iy=2,ny
         do 80 j=1,lot
         do 80 ix=1,nx
         f(ix,iy,j)=0.5*(sipcoy(iy)*wy(iy,ix,j)-simcoy(iy)
     :      *wy(ny-iy+2,ix,j))
80       continue
      else
c
c inverse x transform:
c
         do 110 iy=1,ny
         do 110 j=1,lot
         wx(1,iy,j)=f(1,iy,j)
110      continue
         do 120 j=1,lot
         do 120 iy=1,ny
         do 120 ix=2,nx
         wx(ix,iy,j)=sipcox(ix)*f(ix,iy,j)+simcox(ix)*f(nx-ix+2,iy,j)
120      continue
c
         call fft991(wx,wf,trigsx,ifax,1,nx+2,nx,lot*ny,-1)
c
         do 130 iy=1,ny
         do 130 j=1,lot
         f(1,iy,j)=wx(1,iy,j)
         f(nx,iy,j)=wx(nx+1,iy,j)
130      continue
         do 140 iy=1,ny
         do 140 j=1,lot
         do 140 ii=2,nx-2,2
         f(ii,iy,j)=wx(ii+1,iy,j)-wx(ii+2,iy,j)
         f(ii+1,iy,j)=wx(ii+1,iy,j)+wx(ii+2,iy,j)
140      continue
c
c inverse y transform:
c
         do 150 j=1,lot
         do 150 ix=1,nx
         wy(1,ix,j)=f(ix,1,j)
150      continue
         do 160 iy=2,ny
         do 160 j=1,lot
         do 160 ix=1,nx
         wy(iy,ix,j)=sipcoy(iy)*f(ix,iy,j)+simcoy(iy)*f(ix,ny-iy+2,j)
160      continue
c
         call fft991(wy,wf,trigsy,ifay,1,ny+2,ny,lot*nx,-1)
c
         do 170 j=1,lot
         do 170 ix=1,nx
         f(ix,1,j)=wy(1,ix,j)
         f(ix,ny,j)=wy(ny+1,ix,j)
170      continue
         do 180 ii=2,ny-2,2
         do 180 j=1,lot
         do 180 ix=1,nx
         f(ix,ii,j)=wy(ii+1,ix,j)-wy(ii+2,ix,j)
         f(ix,ii+1,j)=wy(ii+1,ix,j)+wy(ii+2,ix,j)
180      continue
      endif
      return
      end


c dscfft1 -- extended version of dscfft: computes level is = ns+1 also
c
      subroutine dscfft1(f,wf,wx,wy,trigsx,ifax,sipcox,simcox,nx
     ;   ,trigsy,ifay,sipcoy,simcoy,ny,lot,isign)
c-----------------------------------------------------------------------
c double cosine transform on a staggered grid.
c
c for each is=1,ns performs the calculation of:
c
c    for i=1,nx;j=1,ny
c       sum(ix=1,nx;iy=1,ny) f(ix,iy,is)*cos((ix-0.5)*(i-1)*pi/nx)
c                            *cos((iy-0.5)*(j-1)*pi/ny)
c
c following wilhelmson and ericksen (j.c.p. 25,319-331) in the design
c of the pre and posprocessing (making a staggered cosine transform with
c a standard real transform).
c
c in the vectorization it has been assumed (in this routine, not in
c fft991) that nx>ny>ns.
c
c note the dimensions of f (they are appropriate for nh3d program)
c
c p. m. nov 1988
c-----------------------------------------------------------------------
c
      implicit real*8(a-h,o-z)
C      dimension f(-1:nx+1,-1:ny+1,0:lot+2)
C      dimension wf(-1:nx+1,-1:ny+1,0:lot+2)
      dimension f(-1:nx+1,-1:ny+1,lot+1)
      dimension wf(-1:nx+1,-1:ny+1,lot+1)
c
      dimension wx(nx+2,ny,lot+1),wy(ny+2,nx,lot+1)
c
      dimension sipcox(nx),simcox(nx)
      dimension sipcoy(ny),simcoy(ny)
      dimension trigsx(*),trigsy(*)
      dimension ifax(10),ifay(10)
c
      lot1 = lot + 1
C
      if(isign.eq.1) then
c
c x transform:
c
         do 10 iy=1,ny
         do 10 j=1,lot1
         wx(1,iy,j)=f(1,iy,j)
         wx(2,iy,j)=0.
         wx(nx+1,iy,j)=f(nx,iy,j)
         wx(nx+2,iy,j)=0.
10       continue
         do 20 iy=1,ny
         do 20 j=1,lot1
         do 20 ii=2,nx-2,2
         wx(ii+1,iy,j)=0.5*(f(ii,iy,j)+f(ii+1,iy,j))
         wx(ii+2,iy,j)=-0.5*(f(ii,iy,j)-f(ii+1,iy,j))
20       continue
c
         call fft991(wx,wf,trigsx,ifax,1,nx+2,nx,lot1*ny,+1)
c
         do 30 iy=1,ny
         do 30 j=1,lot1
         f(1,iy,j)=wx(1,iy,j)
30       continue
         do 40 iy=1,ny
         do 40 j=1,lot1
         do 40 ix=2,nx
         f(ix,iy,j)=0.5*(sipcox(ix)*wx(ix,iy,j)-simcox(ix)
     :      *wx(nx-ix+2,iy,j))
40       continue
c
c y transform:
c
         do 50 j=1,lot1
         do 50 ix=1,nx
         wy(1,ix,j)=f(ix,1,j)
         wy(2,ix,j)=0.
         wy(ny+1,ix,j)=f(ix,ny,j)
         wy(ny+2,ix,j)=0.
50       continue
         do 60 ii=2,ny-2,2
         do 60 j=1,lot1
         do 60 ix=1,nx
         wy(ii+1,ix,j)=0.5*(f(ix,ii,j)+f(ix,ii+1,j))
         wy(ii+2,ix,j)=-0.5*(f(ix,ii,j)-f(ix,ii+1,j))
60       continue
c
         call fft991(wy,wf,trigsy,ifay,1,ny+2,ny,lot1*nx,+1)
c
         do 70 j=1,lot1
         do 70 ix=1,nx
         f(ix,1,j)=wy(1,ix,j)
70       continue
         do 80 iy=2,ny
         do 80 j=1,lot1
         do 80 ix=1,nx
         f(ix,iy,j)=0.5*(sipcoy(iy)*wy(iy,ix,j)-simcoy(iy)
     :      *wy(ny-iy+2,ix,j))
80       continue
      else
c
c inverse x transform:
c
         do 110 iy=1,ny
         do 110 j=1,lot1
         wx(1,iy,j)=f(1,iy,j)
110      continue
         do 120 j=1,lot1
         do 120 iy=1,ny
         do 120 ix=2,nx
         wx(ix,iy,j)=sipcox(ix)*f(ix,iy,j)+simcox(ix)*f(nx-ix+2,iy,j)
120      continue
c
         call fft991(wx,wf,trigsx,ifax,1,nx+2,nx,lot1*ny,-1)
c
         do 130 iy=1,ny
         do 130 j=1,lot1
         f(1,iy,j)=wx(1,iy,j)
         f(nx,iy,j)=wx(nx+1,iy,j)
130      continue
         do 140 iy=1,ny
         do 140 j=1,lot1
         do 140 ii=2,nx-2,2
         f(ii,iy,j)=wx(ii+1,iy,j)-wx(ii+2,iy,j)
         f(ii+1,iy,j)=wx(ii+1,iy,j)+wx(ii+2,iy,j)
140      continue
c
c inverse y transform:
c
         do 150 j=1,lot1
         do 150 ix=1,nx
         wy(1,ix,j)=f(ix,1,j)
150      continue
         do 160 iy=2,ny
         do 160 j=1,lot1
         do 160 ix=1,nx
         wy(iy,ix,j)=sipcoy(iy)*f(ix,iy,j)+simcoy(iy)*f(ix,ny-iy+2,j)
160      continue
c
         call fft991(wy,wf,trigsy,ifay,1,ny+2,ny,lot1*nx,-1)
c
         do 170 j=1,lot1
         do 170 ix=1,nx
         f(ix,1,j)=wy(1,ix,j)
         f(ix,ny,j)=wy(ny+1,ix,j)
170      continue
         do 180 ii=2,ny-2,2
         do 180 j=1,lot1
         do 180 ix=1,nx
         f(ix,ii,j)=wy(ii+1,ix,j)-wy(ii+2,ix,j)
         f(ix,ii+1,j)=wy(ii+1,ix,j)+wy(ii+2,ix,j)
180      continue
      endif
      return
      end
c




      subroutine raduv(nx,ny,ns
     :   ,u,um1,um2,ubx,ubx2,uby,uby2,ucc,ucc2,ds1y,yl,pp10p1
     :   ,v,vm1,vm2,vbx,vbx2,vby,vby2,vcc,vcc2,ds1x,xl,pp01p1
     :   ,ttmin,masout,uubout,prt,nchn1,ds0,ds1,s,grids0
     :   ,xlcor,xrcor,ylcor,yrcor,ifluxcor,dpmeddt,dx,dy)
c---------------------------------------------------------------------
c decide the modified lateral radiative b.c. for u
c
c in his version it is ensured that the total mass flux is zero.
c
c new 2006
c old scheme, put ifluxcor=1 Zero total mass flux  (OLD scheme)
c                         =0 No correction
c new scheme  put         =2 Adjust total mass flux to the given change in surface pressure
c
c In the OLD scheme the correction is done homogeneously at each boundary
c In the NEW scheme (2) the Inflow and Outflow are computed separately at each
c    boundary grid box (depending only on (xlcor, etc) and the corrections are
c    distributed by grid box in proportion of its relative contribution to the total flow.
c    One first computes the inbalance between the boundary mass flow and
c    the mean surface pressure tendency (residual) and then distribute it between
c    the boundary grid points.
c note: 1/g cancels
c
c---------------------------------------------------------------------
c
      implicit real*8(a-h,o-z)
      real(kind=8),allocatable::xflow(:,:,:),yflow(:,:,:)
      integer ifluxcor
      dimension u(0:nx+1,0:ny+1,0:ns+1),v(0:nx+1,0:ny+1,0:ns+1)
      dimension um1(0:nx+1,0:ny+1,0:ns+1),vm1(0:nx+1,0:ny+1,0:ns+1)
      dimension um2(0:nx+1,0:ny+1,0:ns+1),vm2(0:nx+1,0:ny+1,0:ns+1)
      dimension s(0:nx+1,0:ny+1,0:ns+1)
      dimension pp10p1(0:nx+1,0:ny+1),pp01p1(0:nx+1,0:ny+1)
      dimension ubx(2,0:ny+1,0:ns+1),ubx2(2,0:ny+1,0:ns+1)
      dimension vbx(2,0:ny+1,0:ns+1),vbx2(2,0:ny+1,0:ns+1)
      dimension uby(0:nx+1,2,0:ns+1),uby2(0:nx+1,2,0:ns+1)
      dimension vby(0:nx+1,2,0:ns+1),vby2(0:nx+1,2,0:ns+1)
      dimension ucc(2,2,0:ns+1),ucc2(2,2,0:ns+1)
      dimension vcc(2,2,0:ns+1),vcc2(2,2,0:ns+1)
      dimension ds1x(0:ns+1),ds1y(0:ns+1),ds0(0:ns+1),ds1(0:ns+1)
      logical masout,uubout,prt,grids0,fluxcor
c
c radiative boundary conditions for u and v:
c
c      write(*,*) 'raduv: call radbch u'
      call radbch(ucc,ubx,uby,u,um1,um2,nx+1,ny+1,ns+1
     :   ,1,ns-1,1,nx,1,ny+1)
c      write(*,*) 'raduv: call radbch v'
      call radbch(vcc,vbx,vby,v,vm1,vm2,nx+1,ny+1,ns+1
     :   ,1,ns-1,1,nx+1,1,ny)
c
c mass fluxes through boundaries:
c calculation of the integral of u dy dp/g = uu /g dy ds , etc.
c
      if(ifluxcor.eq.1) then
        xlmflx=0.0
        xlauu=0.0
        xrmflx=0.0
        xrauu=0.0
        do 10 iy=2,ny
        do 10 is=1,ns-1
        xlmflx=xlmflx+ubx(1,iy,is)*ds1y(is)*pp10p1(1,iy)
        xlauu=xlauu+abs(ubx(1,iy,is)*ds1y(is)*pp10p1(1,iy))
        xrmflx=xrmflx+ubx(2,iy,is)*ds1y(is)*pp10p1(nx,iy)
        xrauu=xrauu+abs(ubx(2,iy,is)*ds1y(is)*pp10p1(nx,iy))
10      continue
c
        ylmflx=0.0
        ylavv=0.0
        yrmflx=0.0
        yravv=0.0
        do 30 is=1,ns-1
        do 30 ix=2,nx
        ylmflx=ylmflx+vby(ix,1,is)*ds1x(is)*pp01p1(ix,1)
        ylavv=ylavv+abs(vby(ix,1,is)*ds1x(is)*pp01p1(ix,1))
        yrmflx=yrmflx+vby(ix,2,is)*ds1x(is)*pp01p1(ix,ny)
        yravv=yravv+abs(vby(ix,2,is)*ds1x(is)*pp01p1(ix,ny))
30      continue
c
        ttmflx=xlmflx-xrmflx+ylmflx-yrmflx
c
c ensuring zero mass flux:
c
c a correction is made to uu(ix=1,nx) and vv(iy=1,ny) to ensure zero
c mass flux. the correction is constant at each surface but varies
c between surfaces depending on the absolute value of the variables
c before correction.
c if xlcor eq 0 the correction is not made for xleft surface, and the
c same for other surfaces. this is to allow for the maintenance of
c upstream conditions.
c
        xlauu=xlcor*xlauu
        xrauu=xrcor*xrauu
        ylavv=ylcor*ylavv
        yravv=yrcor*yravv
c
        ttabs=xlauu+xrauu+ylavv+yravv
        if(ttabs.gt.ttmin) then
           duul=xlauu/ttabs*ttmflx/yl
           duur=xrauu/ttabs*ttmflx/yl
           dvvl=ylavv/ttabs*ttmflx/xl
           dvvr=yravv/ttabs*ttmflx/xl
c
           do 110 iy=2,ny
           do 110 is=1,ns-1
           ubx(1,iy,is)=ubx(1,iy,is)-duul/pp10p1(1,iy)
           ubx(2,iy,is)=ubx(2,iy,is)+duur/pp10p1(nx,iy)
110        continue
c
           do 130 is=1,ns-1
           do 130 ix=2,nx
           vby(ix,1,is)=vby(ix,1,is)-dvvl/pp01p1(ix,1)
           vby(ix,2,is)=vby(ix,2,is)+dvvr/pp01p1(ix,ny)
130        continue
        else
           duul=0.
           duur=0.
           dvvl=0.
           dvvr=0.
        endif
      elseif(ifluxcor.eq.2) then

! compute total flow accross boundaries and the inbalance
! (residual) between the integrated boundary fluxes and d pmed/dt
! Correct boundary velocity to impose mass balance, distribute corrections
! proportional to the initial grid point flow (if xlcor=1, etc).

        if(.not.allocated(xflow)) allocate(xflow(2,1:ny+1,1:ns-1))
        if(.not.allocated(yflow)) allocate(yflow(1:nx+1,2,1:ns-1))

        flowinx=0.
        flowoutx=0.
        flowincx=0.
        flowoutcx=0.

        flowiny=0.
        flowouty=0.
        flowincy=0.
        flowoutcy=0.

        do iy=2,ny
        do is=1,ns-1
          xflow(1,iy,is)=ubx(1,iy,is)*ds1y(is)*pp10p1(1,iy)
          xflow(2,iy,is)=ubx(2,iy,is)*ds1y(is)*pp10p1(nx,iy)
        enddo
        enddo

        do iy=2,ny
        do is=1,ns-1
          flowinx=flowinx+max(0.,xflow(1,iy,is))-min(0.,xflow(2,iy,is))
          flowoutx=flowoutx-min(0.,xflow(1,iy,is))
     :       +max(0.,xflow(2,iy,is))
        enddo
        enddo

        if(xlcor.eq.1.) then
          do iy=2,ny
          do is=1,ns-1
            flowincx=flowincx+max(0.,xflow(1,iy,is))
            flowoutcx=flowoutcx-min(0.,xflow(1,iy,is))
          enddo
          enddo
        endif

        if(xrcor.eq.1.) then
          do iy=2,ny
          do is=1,ns-1
            flowincx=flowincx-min(0.,xflow(2,iy,is))
            flowoutcx=flowoutcx+max(0.,xflow(2,iy,is))
          enddo
          enddo
        endif

        do is=1,ns-1
        do ix=2,nx
          yflow(ix,1,is)=vby(ix,1,is)*ds1x(is)*pp01p1(ix,1)
          yflow(ix,2,is)=vby(ix,2,is)*ds1x(is)*pp01p1(ix,ny)
        enddo
        enddo

        do is=1,ns-1
        do ix=2,nx
          flowiny=flowiny+max(0.,yflow(ix,1,is))
     :      -min(0.,yflow(ix,2,is))
          flowouty=flowouty-min(0.,yflow(ix,1,is))
     :      +max(0.,yflow(ix,2,is))
        enddo
        enddo

        if(ylcor.eq.1.) then
          do is=1,ns-1
          do ix=2,nx
            flowincy=flowincy+max(0.,yflow(ix,1,is))
            flowoutcy=flowoutcy-min(0.,yflow(ix,1,is))
          enddo
          enddo
        endif

        if(yrcor.eq.1.) then
          do is=1,ns-1
          do ix=2,nx
            flowincy=flowincy-min(0.,yflow(ix,2,is))
            flowoutcy=flowoutcy+max(0.,yflow(ix,2,is))
          enddo
          enddo
        endif

        flowin=flowinx+flowiny
        flowinc=flowincx+flowincy

        flowout=flowoutx+flowouty
        flowoutc=flowoutcx+flowoutcy

        area=(nx-1)*dx*(ny-1)*dy

c note: 1/g cancels

        residual=flowin-flowout-dpmeddt*area
        totalcf=flowinc+flowoutc
!       write(*,*) 'residual:',residual,flowin,flowout,dpmeddt*area
        if(residual.gt.1.e-6) then
          if(xlcor.eq.1.) then
            do iy=2,ny
            do is=1,ns-1
              ubx(1,iy,is)=ubx(1,iy,is)-residual*abs(xflow(1,iy,is))
     :          /(totalcf*pp10p1(1,iy)*ds1y(is))
            enddo
            enddo
          endif

          if(xrcor.eq.1.) then
            do iy=2,ny
            do is=1,ns-1
              ubx(2,iy,is)=ubx(2,iy,is)+residual*abs(xflow(2,iy,is))
     :          /(totalcf*pp10p1(nx,iy)*ds1y(is))
            enddo
            enddo
          endif
c
          if(ylcor.eq.1.) then
            do is=1,ns-1
            do ix=2,nx
              vby(ix,1,is)=vby(ix,1,is)-residual*abs(yflow(ix,1,is))
     :          /(totalcf*pp01p1(ix,1)*ds1x(is))
            enddo
            enddo
          endif

          if(yrcor.eq.1.) then
            do is=1,ns-1
            do ix=2,nx
              vby(ix,2,is)=vby(ix,2,is)+residual*abs(yflow(ix,2,is))
     :          /(totalcf*pp01p1(ix,ny)*ds1x(is))
            enddo
            enddo
          endif
        endif
c debug verify correction

!        flowinx=0.
!        flowoutx=0.
!
!        flowiny=0.
!        flowouty=0.
!
!        do iy=2,ny
!        do is=1,ns-1
!          xflow(1,iy,is)=ubx(1,iy,is)*ds1y(is)*pp10p1(1,iy)
!          xflow(2,iy,is)=ubx(2,iy,is)*ds1y(is)*pp10p1(nx,iy)
!        enddo
!        enddo
!
!        do iy=2,ny
!        do is=1,ns-1
!          flowinx=flowinx+max(0.,xflow(1,iy,is))-min(0.,xflow(2,iy,is))
!          flowoutx=flowoutx-min(0.,xflow(1,iy,is))
!     :       +max(0.,xflow(2,iy,is))
!        enddo
!        enddo
!
!        do is=1,ns-1
!        do ix=2,nx
!          yflow(ix,1,is)=vby(ix,1,is)*ds1x(is)*pp01p1(ix,1)
!          yflow(ix,2,is)=vby(ix,2,is)*ds1x(is)*pp01p1(ix,ny)
!        enddo
!        enddo
!
!        do is=1,ns-1
!        do ix=2,nx
!          flowiny=flowiny+max(0.,yflow(ix,1,is))
!     :      -min(0.,yflow(ix,2,is))
!          flowouty=flowouty-min(0.,yflow(ix,1,is))
!     :      +max(0.,yflow(ix,2,is))
!        enddo
!        enddo
!
!        flowin=flowinx+flowiny
!        flowout=flowoutx+flowouty
!
!        area=(nx-1)*dx*(ny-1)*dy
!
!c note: 1/g cancels
!
!        residual=flowin-flowout-dpmeddt*area
!        write(*,*) 'final residual:',residual,flowin,flowout
!     :    ,dpmeddt*area
      endif
c
c extend uu and vv to the outer ring:
c
c      write(*,*) 'raduv: call radbch u2'
      call radbch(ucc2,ubx2,uby2,u,um1,um2,nx+1,ny+1,ns+1
     :   ,1,ns-1,0,nx+1,0,ny+1)
c      write(*,*) 'raduv: call radbch v2'
      call radbch(vcc2,vbx2,vby2,v,vm1,vm2,nx+1,ny+1,ns+1
     :   ,1,ns-1,0,nx+1,0,ny+1)
c
      if(ifluxcor.eq.1) then
        xlmfl2=0.0
        xlauu=0.0
        xrmfl2=0.0
        xrauu=0.0
        do 210 iy=2,ny
        do 210 is=1,ns-1
        xlmfl2=xlmfl2+ubx2(1,iy,is)*ds1y(is)*pp10p1(0,iy)
        xlauu=xlauu+abs(ubx2(1,iy,is)*ds1y(is)*pp10p1(0,iy))
        xrmfl2=xrmfl2+ubx2(2,iy,is)*ds1y(is)*pp10p1(nx+1,iy)
        xrauu=xrauu+abs(ubx2(2,iy,is)*ds1y(is)*pp10p1(nx+1,iy))
210     continue
c
        ylmfl2=0.0
        ylavv=0.0
        yrmfl2=0.0
        yravv=0.0
        do 230 is=1,ns-1
        do 230 ix=2,nx
        ylmfl2=ylmfl2+vby2(ix,1,is)*ds1x(is)*pp01p1(ix,0)
        ylavv=ylavv+abs(vby2(ix,1,is)*ds1x(is)*pp01p1(ix,0))
        yrmfl2=yrmfl2+vby2(ix,2,is)*ds1x(is)*pp01p1(ix,ny+1)
        yravv=yravv+abs(vby2(ix,2,is)*ds1x(is)*pp01p1(ix,ny+1))
230     continue
c
        ttmfl2=xlmfl2-xrmfl2+ylmfl2-yrmfl2
c
c ensuring zero mass flux:
c
c a correction is made to uu(iy=1,ny) and vv(ix=1,nx) to ensure zero
c mass flux. the correction is constant at each surface but varies
c between surfaces depending on the absolute value of the variables
c before correction.
c
        xlauu=xlcor*xlauu
        xrauu=xrcor*xrauu
        ylavv=ylcor*ylavv
        yravv=yrcor*yravv
c
        ttabs=xlauu+xrauu+ylavv+yravv
        if(ttabs.gt.ttmin) then
           duul2=xlauu/ttabs*ttmfl2/yl
           duur2=xrauu/ttabs*ttmfl2/yl
           dvvl2=ylavv/ttabs*ttmfl2/xl
           dvvr2=yravv/ttabs*ttmfl2/xl
c
           do 310 iy=1,ny+1
           do 310 is=1,ns-1
           ubx2(1,iy,is)=ubx2(1,iy,is)-duul2/pp10p1(0,iy)
           ubx2(2,iy,is)=ubx2(2,iy,is)+duur2/pp10p1(nx+1,iy)
310        continue
c
           do 330 is=1,ns-1
           do 330 ix=1,nx+1
           vby2(ix,1,is)=vby2(ix,1,is)-dvvl2/pp01p1(ix,0)
           vby2(ix,2,is)=vby2(ix,2,is)+dvvr2/pp01p1(ix,ny+1)
330        continue
        else
           duul2=0.
           duur2=0.
           dvvl2=0.
           dvvr2=0.
        endif
c
        if(masout .and. prt) then
           write(nchn1,*)'xlcor=',xlcor,' xrcor=',xrcor
           write(nchn1,*)'ylcor=',ylcor,' yrcor=',yrcor
           write(nchn1,1) xlmflx,xrmflx,ylmflx,yrmflx,ttmflx
     :        ,duul,duur,dvvl,dvvr
           write(nchn1,2) xlmfl2,xrmfl2,ylmfl2,yrmfl2,ttmfl2
     :        ,duul2,duur2,dvvl2,dvvr2
        endif
      elseif(ifluxcor.eq.2) then
        flowinx=0.
        flowoutx=0.
        flowincx=0.
        flowoutcx=0.


        flowiny=0.
        flowouty=0.
        flowincy=0.
        flowoutcy=0.


        do iy=1,ny+1
        do is=1,ns-1
          xflow(1,iy,is)=ubx2(1,iy,is)*ds1y(is)*pp10p1(0,iy)
          xflow(2,iy,is)=ubx2(2,iy,is)*ds1y(is)*pp10p1(nx+1,iy)
        enddo
        enddo


        do iy=1,ny+1
        do is=1,ns-1
          flowinx=flowinx+max(0.,xflow(1,iy,is))-min(0.,xflow(2,iy,is))
          flowoutx=flowoutx-min(0.,xflow(1,iy,is))
     :       +max(0.,xflow(2,iy,is))
        enddo
        enddo

        if(xlcor.eq.1.) then
          do iy=1,ny+1
          do is=1,ns-1
            flowincx=flowincx+max(0.,xflow(1,iy,is))
            flowoutcx=flowoutcx-min(0.,xflow(1,iy,is))
          enddo
          enddo
        endif

        if(xrcor.eq.1.) then
          do iy=1,ny+1
          do is=1,ns-1
            flowincx=flowincx-min(0.,xflow(2,iy,is))
            flowoutcx=flowoutcx+max(0.,xflow(2,iy,is))
          enddo
          enddo
        endif


        do is=1,ns-1
        do ix=1,nx+1
          yflow(ix,1,is)=vby2(ix,1,is)*ds1x(is)*pp01p1(ix,0)
          yflow(ix,2,is)=vby2(ix,2,is)*ds1x(is)*pp01p1(ix,ny+1)
        enddo
        enddo


        do is=1,ns-1
        do ix=1,nx+1
          flowiny=flowiny+max(0.,yflow(ix,1,is))
     :      -min(0.,yflow(ix,2,is))
          flowouty=flowouty-min(0.,yflow(ix,1,is))
     :      +max(0.,yflow(ix,2,is))
        enddo
        enddo

        if(ylcor.eq.1.) then
          do is=1,ns-1
          do ix=1,nx+1
            flowincy=flowincy+max(0.,yflow(ix,1,is))
            flowoutcy=flowoutcy-min(0.,yflow(ix,1,is))
          enddo
          enddo
        endif

        if(yrcor.eq.1.) then
          do is=1,ns-1
          do ix=1,nx+1
            flowincy=flowincy-min(0.,yflow(ix,2,is))
            flowoutcy=flowoutcy+max(0.,yflow(ix,2,is))
          enddo
          enddo
        endif

        flowin=flowinx+flowiny
        flowinc=flowincx+flowincy

        flowout=flowoutx+flowouty
        flowoutc=flowoutcx+flowoutcy

        area=(nx+1)*dx*(ny+1)*dy

c note: 1/g cancels

        residual=flowin-flowout-dpmeddt*area
        totalcf=flowinc+flowoutc

        if(xlcor.eq.1.) then
          do iy=1,ny+1
          do is=1,ns-1
            ubx2(1,iy,is)=ubx2(1,iy,is)-residual*abs(xflow(1,iy,is))
     :        /(totalcf*pp10p1(0,iy)*ds1y(is))
          enddo
          enddo
        endif

        if(xrcor.eq.1.) then
          do iy=1,ny+1
          do is=1,ns-1
            ubx2(2,iy,is)=ubx2(2,iy,is)+residual*abs(xflow(2,iy,is))
     :        /(totalcf*pp10p1(nx+1,iy)*ds1y(is))
          enddo
          enddo
        endif
c
        if(ylcor.eq.1.) then
          do is=1,ns-1
          do ix=1,nx+1
            vby2(ix,1,is)=vby2(ix,1,is)-residual*abs(yflow(ix,1,is))
     :        /(totalcf*pp01p1(ix,0)*ds1x(is))
          enddo
          enddo
        endif

        if(yrcor.eq.1.) then
          do is=1,ns-1
          do ix=1,nx+1
            vby2(ix,2,is)=vby2(ix,2,is)+residual*abs(yflow(ix,2,is))
     :        /(totalcf*pp01p1(ix,ny+1)*ds1x(is))
          enddo
          enddo
        endif

      endif


c
      if(uubout .and. prt) then
        call wrigar(ubx,1,2,0,ny+1,0,ns+1,1,2,1,ny+1,1,ns-1
     :     ,'ubx     ',0.d0,3)
        call wrigar(ubx2,1,2,0,ny+1,0,ns+1,1,2,0,ny+1,1,ns-1
     :     ,'ubx2    ',0.d0,3)
        call wrigar(uby,0,nx+1,1,2,0,ns+1,1,nx+1,1,2,1,ns-1
     :     ,'uby     ',0.d0,2)
        call wrigar(uby2,0,nx+1,1,2,0,ns+1,0,nx+1,1,2,1,ns-1
     :     ,'uby2    ',0.d0,2)
        call wrigar(vbx,1,2,0,ny+1,0,ns+1,1,2,1,ny+1,1,ns-1
     :     ,'vbx     ',0.d0,3)
        call wrigar(vbx2,1,2,0,ny+1,0,ns+1,1,2,0,ny+1,1,ns-1
     :     ,'vbx2    ',0.d0,3)
        call wrigar(vby,0,nx+1,1,2,0,ns+1,1,nx+1,1,2,1,ns-1
     :     ,'vby     ',0.d0,2)
        call wrigar(vby2,0,nx+1,1,2,0,ns+1,0,nx+1,1,2,1,ns-1
     :     ,'vby2    ',0.d0,2)
c
        write(nchn1,*) 'ucc,ucc2'
        do 777 is=1,ns-1
        write(nchn1,*) is,ucc(1,1,is),ucc(2,1,is),ucc(2,2,is)
     :     ,ucc(1,2,is),ucc2(1,1,is),ucc2(2,1,is)
777     continue
        write(nchn1,*) 'vcc,vcc2'
        do 779 is=1,ns-1
        write(nchn1,*) is,vcc(1,1,is),vcc(2,1,is),vcc(2,2,is)
     :     ,vcc(1,2,is),vcc2(1,1,is),vcc2(1,2,is)
779     continue
      endif
c
      return
c
1     format(' xlmflx=',e15.8,' xrmflx=',e15.8,' ylmflx=',e15.8,
     :   ' yrmflx=',e15.8, ' ttmflx=',e15.8,/
     :   ,' duul=',e15.8,' duur=',e15.8
     :   /,' dvvl=',e15.8,' dvvr=',e15.8)
2     format(' xlmfl2=',e15.8,' xrmfl2=',e15.8,' ylmfl2=',e15.8,
     :   ' yrmfl2=',e15.8, ' ttmfl2=',e15.8,/
     :   ,' duul2=',e15.8,' duur2=',e15.8
     :   /,' dvvl2=',e15.8,' dvvr2=',e15.8)
      end
c
c
c
      subroutine radbco(abcc,abcx,abcy,a,am1,am2,nxm,nym,nsm
     :   ,k0,k1,ixl,ixr,iyl,iyr)
c
c-----------------------------------------------------------------------
c to compute radiative boundary conditions
c for horizontal boundaries of the field a.
c
c note that the values in the corner points depend on the order.
c
c this version is for pp only
c-----------------------------------------------------------------------
c
c
      implicit real*8(a-h,o-z)
      dimension a(0:nxm,0:nym,0:nsm),am1(0:nxm,0:nym,0:nsm)
     :   ,am2(0:nxm,0:nym,0:nsm)
      dimension abcx(2,0:nym,0:nsm),abcy(0:nxm,2,0:nsm),abcc(2,2,0:nsm)
      common/radpar/dx,dy,dt,dxdt,dydt,rxxxx,ryxxx
c
      ccmax=sqrt(dx*dx+dy*dy)/(2.*dt)
      rxmax=ccmax/dxdt
      rymax=ccmax/dydt
c
c x boundaries (ix=ixl,ixr)
c
      ib=ixl
      ib1=ib+1
      ib2=ib+2
      do 10 lr=1,2
      do 20 iy=iyl,iyr
      do 20 is=k0,k1
      rx=(am2(ib1,iy,is)-a(ib1,iy,is))
     :   /(a(ib1,iy,is)+am2(ib1,iy,is)
     :   -2.*am1(ib2,iy,is)+1.e-30)
      rx=dim(rx,0.d0)-dim(rx,rxmax)
      abcx(lr,iy,is)=((1.-rx)*am1(ib,iy,is)+2.*rx
     :   *a(ib1,iy,is))/(1.+rx)
20    continue
      ib=ixr
      ib1=ib-1
      ib2=ib-2
10    continue
c
c y boundaries (iy=iyl,iyr)
c
      ib=iyl
      ib1=ib+1
      ib2=ib+2
      do 110 lr=1,2
      do 120 is=k0,k1
      do 120 ix=ixl,ixr
      ry=(am2(ix,ib1,is)-a(ix,ib1,is))
     :   /(a(ix,ib1,is)+am2(ix,ib1,is)
     :   -2.*am1(ix,ib2,is)+1.e-30)
      ry=dim(ry,0.d0)-dim(ry,rymax)
      abcy(ix,lr,is)=((1.-ry)*am1(ix,ib,is)+2.*ry
     :   *a(ix,ib1,is))/(1.+ry)
120   continue
      ib=iyr
      ib1=ib-1
      ib2=ib-2
110   continue
c
c corner points:
c
      iby=iyl
      iby1=iby+1
      iby2=iby+2
      do 212 ly=1,2
      ibx=ixl
      ibx1=ibx+1
      ibx2=ibx+2
      do 217 lx=1,2
      do 222 is=k0,k1
      cx=(am2(ibx1,iby1,is)-a(ibx1,iby1,is))
     :   /(a(ibx1,iby1,is)+am2(ibx1,iby1,is)
     :   -2.*am1(ibx2,iby1,is)+1.e-30)*dxdt
      cx=dim(cx,0.d0)-dim(cx,ccmax)
      cy=(am2(ibx1,iby1,is)-a(ibx1,iby1,is))
     :   /(a(ibx1,iby1,is)+am2(ibx1,iby1,is)
     :   -2.*am1(ibx1,iby2,is)+1.e-30)*dydt
      cy=dim(cy,0.d0)-dim(cy,ccmax)
      cc=sqrt(cx*cx+cy*cy)
      rcc=cc/ccmax
      a13=0.25*(a(ibx1,iby,is)+a(ibx,iby1,is)+a(ibx2,iby1,is)
     :   +a(ibx1,iby2,is))
      abcc(lx,ly,is)=((1.-rcc)*am1(ibx,iby,is)+2.*rcc*a13)/(1.+rcc)
222   continue
      ibx=ixr
      ibx1=ibx-1
      ibx2=ibx-2
217   continue
      iby=iyr
      iby1=iby-1
      iby2=iby-2
212   continue
c
      return
      end
      subroutine ppbc(nx,ny,pp,ppbx,ppbx2,ppby,ppby2,ppcc,ppcc2
     :   ,iobppx,iobppy,ppbout,prt,nchn1)
c
c updates boundary variables from data obtained in the previous time
c step (radiative boundary conditions).
c
      implicit real*8(a-h,o-z)
      logical ppbout,prt
      dimension pp(0:nx+2,0:ny+2)
      dimension ppbx(2,0:ny+2),ppbx2(2,0:ny+2),ppby(0:nx+2,2)
     :   ,ppby2(0:nx+2,2)
      dimension ppcc(2,2),ppcc2(2,2)
c
      if(ppbout .and. prt) then
         write(nchn1,*) 'ppbx(1),ppbx(2),ppbx2(1),ppbx2(2)'
         do 700 iy=0,ny+2
         write(nchn1,*) iy,ppbx(1,iy),ppbx(2,iy),ppbx2(1,iy),ppbx2(2,iy)
700      continue
         write(nchn1,*) 'ppby(1),ppby(2),ppby2(1),ppby2(2)'
         do 701 ix=0,nx+2
         write(nchn1,*) ix,ppby(ix,1),ppby(ix,2),ppby2(ix,1),ppby2(ix,2)
701      continue
         write(nchn1,*) 'ppcc'
         write(nchn1,*) ppcc(1,1),ppcc(2,1),ppcc(2,2),ppcc(1,2)
         write(nchn1,*) ppcc2(1,1),ppcc2(2,1),ppcc2(2,2),ppcc2(1,2)
      endif
c
      if(iobppx.eq.1) then
         do 23 iy=0,ny+2
         pp(0,iy)=ppbx2(1,iy)
         pp(nx+2,iy)=ppbx2(2,iy)
23       continue
      else
         do 43 iy=0,ny+2
         pp(0,iy)=pp(1,iy)
         pp(nx+2,iy)=pp(nx+1,iy)
43       continue
      endif
c
      if(iobppy.eq.1) then
         do 22 ix=0,nx+2
         pp(ix,0)=ppby2(ix,1)
         pp(ix,ny+2)=ppby2(ix,2)
22       continue
      else
         do 42 ix=0,nx+2
         pp(ix,0)=pp(ix,1)
         pp(ix,ny+2)=pp(ix,ny+1)
42       continue
      endif
c
      if(iobppx.eq.1 .and. iobppy.eq.1) then
         pp(0,0)=ppcc2(1,1)
         pp(nx+2,0)=ppcc2(2,1)
         pp(nx+2,ny+2)=ppcc2(2,2)
         pp(0,ny+2)=ppcc2(1,2)
      endif
c
      return
      end




      subroutine hsmoot(f,work,nxm,nym,nsm,i0,i1,j0,j1,k0,k1,gama)
c
c 4th-order horizontal smoothing. 2nd order for grid-points adjacent
c to the boundary. no smoothing for boundary points.
c
c the second order coeficient is calculated from the fourth order one,
c imposing similar damping for 2 grid-length waves
c
c fourth order damping: response function
c r(lx,ly=)1-8*gama*(sin(pi*dx/lx)**4+sin(pi*dy/ly)**4)
c
      implicit real*8(a-h,o-z)
      dimension f(0:nxm,0:nym,0:nsm),work(0:nxm,0:nym)
      dimension gama(0:*)
      real(kind(0.d0)),allocatable,dimension(:)::
     :  gama05,um6g,gasco4,umgasc
c     dimension gama05(0:maxns),um6g(0:maxns),gasco4(0:maxns)
c    :   ,umgasc(0:maxns)
c
      if(.not.allocated(gama05)) then
        allocate(gama05(0:nsm))
        allocate(um6g(0:nsm))
        allocate(gasco4(0:nsm))
        allocate(umgasc(0:nsm))
      endif

      do is=k0,k1
        gama05(is)=0.5*gama(is)
      um6g(is)=1.-6.*gama(is)
      gasco4(is)=2.*gama(is)
      umgasc(is)=1.-8.*gama(is)
      enddo
c
!      do is=k0,k1
!        do iy=j0+2,j1-2
!        do ix=i0+2,i1-2
!          work(ix,iy)=um6g(is)*f(ix,iy,is)
!     :      -gama05(is)*(f(ix+2,iy,is)+f(ix-2,iy,is)
!     :      +f(ix,iy+2,is)+f(ix,iy-2,is)-4.*(f(ix+1,iy,is)+f(ix-1,iy,is)
!     :      +f(ix,iy+1,is)+f(ix,iy-1,is)))
!        enddo
!        enddo
c
!        do ix=i0+1,i1-1,i1-i0-2
!        do iy=j0+1,j1-1
!          work(ix,iy)=umgasc(is)*f(ix,iy,is)+gasco4(is)*(f(ix-1,iy,is)
!     :      +f(ix+1,iy,is)+f(ix,iy-1,is)+f(ix,iy+1,is))
!        enddo
!        enddo
!c
!        do iy=j0+1,j1-1,j1-j0-2
!        do ix=i0+1,i1-1
!          work(ix,iy)=umgasc(is)*f(ix,iy,is)+gasco4(is)*(f(ix-1,iy,is)
!     :      +f(ix+1,iy,is)+f(ix,iy-1,is)+f(ix,iy+1,is))
!        enddo
!        enddo
c
!        do iy=j0+1,j1-1
!        do ix=i0+1,i1-1
!          f(ix,iy,is)=work(ix,iy)
!        enddo
!        enddo
!      enddo
!---------------------version 2----------------------------------------!
      do is=k0,k1
      
      do iy=j0+1,j1-1
      do ix=i0+2,i1-2
      work(ix,iy)=um6g(is)*f(ix,iy,is)
     :      -gama05(is)*(f(ix+2,iy,is)+f(ix-2,iy,is)
     :      +f(ix,iy,is)+f(ix,iy,is)-4.*(f(ix+1,iy,is)+f(ix-1,iy,is)
     :      +f(ix,iy,is)+f(ix,iy,is)))
      enddo
      do ix=i0+1,i1-1,i1-i0-2
      work(ix,iy)=umgasc(is)*f(ix,iy,is)+gasco4(is)*(f(ix-1,iy,is)
     :      +f(ix+1,iy,is)+f(ix,iy,is)+f(ix,iy,is))
      enddo     
      enddo
      
      do iy=j0+1,j1-1
      do ix=i0+1,i1-1
      f(ix,iy,is)=work(ix,iy)
      enddo
      enddo
      
      do ix=i0+1,i1-1
      do iy=j0+2,j1-2
      work(ix,iy)=um6g(is)*f(ix,iy,is)
     :      -gama05(is)*(f(ix,iy,is)+f(ix,iy,is)
     :      +f(ix,iy+2,is)+f(ix,iy-2,is)-4.*(f(ix,iy,is)+f(ix,iy,is)
     :      +f(ix,iy+1,is)+f(ix,iy-1,is)))
      enddo
      do iy=j0+1,j1-1,j1-j0-2
       work(ix,iy)=umgasc(is)*f(ix,iy,is)+gasco4(is)*(f(ix,iy,is)
     :      +f(ix,iy,is)+f(ix,iy-1,is)+f(ix,iy+1,is))
      enddo
      enddo
      
      do iy=j0+1,j1-1
      do ix=i0+1,i1-1
      f(ix,iy,is)=work(ix,iy)
      enddo
      enddo
      enddo

      return
      end
c
c
c
      subroutine hsmop(f,fs,work,nxm,nym,nsm,i0,i1,j0,j1,k0,k1,gama
     :   ,perdam)
c
c only the perturbation is smoothed (if perdam is .true.)
c
c 4th-order horizontal smoothing. 2nd order for grid-points adjacent
c to the boundary. no smoothing for boundary points.
c
c the second order coeficient is calculated from the fourth order one,
c imposing similar damping for 2 grid-length waves
c
c fourth order damping: response function
c r(lx,ly=)1-8*gama*(sin(pi*dx/lx)**4+sin(pi*dy/ly)**4)
c
      implicit real*8(a-h,o-z)
      dimension f(0:nxm,0:nym,0:nsm),fs(0:nxm,0:nym,0:nsm)
     :   ,work(0:nxm,0:nym)
      logical perdam
      dimension gama(0:*)
c      dimension gama05(0:maxns),um6g(0:maxns),gasco4(0:maxns)
c     :   ,umgasc(0:maxns)
      real(kind(0.d0)),allocatable,dimension(:)::
     :  gama05,um6g,gasco4,umgasc
c
      if(.not.allocated(gama05)) then
        allocate(gama05(0:nsm))
        allocate(um6g(0:nsm))
        allocate(gasco4(0:nsm))
        allocate(umgasc(0:nsm))
      endif
c
      do 1 is=k0,k1
      gama05(is)=0.5*gama(is)
      um6g(is)=1.-6.*gama(is)
      gasco4(is)=2.*gama(is)
      umgasc(is)=1.-8.*gama(is)
1     continue
c
c
      if(perdam) then
         do 11 is=k0,k1
         do 11 iy=j0,j1
         do 11 ix=i0,i1
         f(ix,iy,is)=f(ix,iy,is)-fs(ix,iy,is)
11       continue
      endif
c
!      do 10 is=k0,k1
c
!      do 20 iy=j0+2,j1-2
!      do 20 ix=i0+2,i1-2
!      work(ix,iy)=um6g(is)*f(ix,iy,is)
!     :   -gama05(is)*(f(ix+2,iy,is)+f(ix-2,iy,is)
!     :   +f(ix,iy+2,is)+f(ix,iy-2,is)-4.*(f(ix+1,iy,is)+f(ix-1,iy,is)
!     :   +f(ix,iy+1,is)+f(ix,iy-1,is)))
!20    continue
c
!      do 30 ix=i0+1,i1-1,i1-i0-2
!      do 30 iy=j0+1,j1-1
!      work(ix,iy)=umgasc(is)*f(ix,iy,is)+gasco4(is)*(f(ix-1,iy,is)
!     :   +f(ix+1,iy,is)+f(ix,iy-1,is)+f(ix,iy+1,is))
!30    continue
c
!      do 40 iy=j0+1,j1-1,j1-j0-2
!      do 40 ix=i0+1,i1-1
!      work(ix,iy)=umgasc(is)*f(ix,iy,is)+gasco4(is)*(f(ix-1,iy,is)
!     :   +f(ix+1,iy,is)+f(ix,iy-1,is)+f(ix,iy+1,is))
!40    continue
c
!      do 50 iy=j0+1,j1-1
!      do 50 ix=i0+1,i1-1
!      f(ix,iy,is)=work(ix,iy)
!50    continue
c
!10    continue
c
!---------------------version 2----------------------------------------!
            do is=k0,k1
      
      do iy=j0+1,j1-1
      do ix=i0+2,i1-2
      work(ix,iy)=um6g(is)*f(ix,iy,is)
     :      -gama05(is)*(f(ix+2,iy,is)+f(ix-2,iy,is)
     :      +f(ix,iy,is)+f(ix,iy,is)-4.*(f(ix+1,iy,is)+f(ix-1,iy,is)
     :      +f(ix,iy,is)+f(ix,iy,is)))
      enddo
      do ix=i0+1,i1-1,i1-i0-2
      work(ix,iy)=umgasc(is)*f(ix,iy,is)+gasco4(is)*(f(ix-1,iy,is)
     :      +f(ix+1,iy,is)+f(ix,iy,is)+f(ix,iy,is))
      enddo     
      enddo
      
      do iy=j0+1,j1-1
      do ix=i0+1,i1-1
      f(ix,iy,is)=work(ix,iy)
      enddo
      enddo
      
      do ix=i0+1,i1-1
      do iy=j0+2,j1-2
      work(ix,iy)=um6g(is)*f(ix,iy,is)
     :      -gama05(is)*(f(ix,iy,is)+f(ix,iy,is)
     :      +f(ix,iy+2,is)+f(ix,iy-2,is)-4.*(f(ix,iy,is)+f(ix,iy,is)
     :      +f(ix,iy+1,is)+f(ix,iy-1,is)))
      enddo
      do iy=j0+1,j1-1,j1-j0-2
       work(ix,iy)=umgasc(is)*f(ix,iy,is)+gasco4(is)*(f(ix,iy,is)
     :      +f(ix,iy,is)+f(ix,iy-1,is)+f(ix,iy+1,is))
      enddo
      enddo
      
      do iy=j0+1,j1-1
      do ix=i0+1,i1-1
      f(ix,iy,is)=work(ix,iy)
      enddo
      enddo
      enddo
      
      if(perdam) then
         do 12 is=k0,k1
         do 12 iy=j0,j1
         do 12 ix=i0,i1
         f(ix,iy,is)=f(ix,iy,is)+fs(ix,iy,is)
12       continue
      endif
c
      return
      end
c
c
c
      subroutine vsmoot(f,work,nxm,nym,nsm,i0,i1,j0,j1,k0,k1,gama)
c
c 4th-order vertical smoothing. 2nd order for grid-points adjacent
c to the boundary. no smoothing for boundary points.
c
c the second order coeficient is calculated from the fourth order one,
c imposing similar damping for 2 grid-length waves
c
c fourth order damping: response function
c r(lx,ly=)1-16*gama*sin(pi*dx/lx)**4
c
      implicit real*8(a-h,o-z)
      dimension f(0:nxm,0:nym,0:nsm),work(0:nxm,0:nym,0:nsm)
c
      um6g=1.-6.*gama
      gasc=8.*gama
      gasc05=0.5*gasc
      umgasc=1.-gasc
c
      do 20 is=k0+2,k1-2
      do 20 iy=j0,j1
      do 20 ix=i0,i1
      work(ix,iy,is)=um6g*f(ix,iy,is)-gama*(f(ix,iy,is+2)+f(ix,iy,is-2)
     :   -4.*(f(ix,iy,is+1)+f(ix,iy,is-1)))
20    continue
c
      do 30 is=k0+1,k1-1,k1-k0-2
      do 30 iy=j0,j1
      do 30 ix=i0,i1
      work(ix,iy,is)=umgasc*f(ix,iy,is)+gasc05*(f(ix,iy,is-1)
     :   +f(ix,iy,is+1))
30    continue
c
c
      do 50 is=k0+1,k1-1
      do 50 iy=j0,j1
      do 50 ix=i0,i1
      f(ix,iy,is)=work(ix,iy,is)
50    continue
c
      return
      end
c
c
c
      subroutine vsmop(f,fs,work,nxm,nym,nsm,i0,i1,j0,j1,k0,k1,gama
     :   ,perdam)
c
c only the perturbation is smoothed (if perdam is .true.)
c
c 4th-order vertical smoothing. 2nd order for grid-points adjacent
c to the boundary. no smoothing for boundary points.
c
c the second order coeficient is calculated from the fourth order one,
c imposing similar damping for 2 grid-length waves
c
c fourth order damping: response function
c r(lx,ly=)1-16*gama*sin(pi*dx/lx)**4
c
      implicit real*8(a-h,o-z)
      dimension f(0:nxm,0:nym,0:nsm),fs(0:nxm,0:nym,0:nsm)
     :   ,work(0:nxm,0:nym,0:nsm)
      logical perdam
c
      um6g=1.-6.*gama
      gasc=8.*gama
      gasc05=0.5*gasc
      umgasc=1.-gasc
c
      if(perdam) then
         do 10 is=k0,k1
         do 10 iy=j0,j1
         do 10 ix=i0,i1
         f(ix,iy,is)=f(ix,iy,is)-fs(ix,iy,is)
10       continue
      endif
c
      do 20 is=k0+2,k1-2
      do 20 iy=j0,j1
      do 20 ix=i0,i1
      work(ix,iy,is)=um6g*f(ix,iy,is)-gama*(f(ix,iy,is+2)+f(ix,iy,is-2)
     :   -4.*(f(ix,iy,is+1)+f(ix,iy,is-1)))
20    continue
c
      do 30 is=k0+1,k1-1,k1-k0-2
      do 30 iy=j0,j1
      do 30 ix=i0,i1
      work(ix,iy,is)=umgasc*f(ix,iy,is)+gasc05*(f(ix,iy,is-1)
     :   +f(ix,iy,is+1))
30    continue
c
      do 50 is=k0+1,k1-1
      do 50 iy=j0,j1
      do 50 ix=i0,i1
      f(ix,iy,is)=work(ix,iy,is)
50    continue
c
      if(perdam) then
         do 60 is=k0,k1
         do 60 iy=j0,j1
         do 60 ix=i0,i1
         f(ix,iy,is)=f(ix,iy,is)+fs(ix,iy,is)
60       continue
      endif
c
      return
      end
c
c
c
      subroutine aselin(f,fm1,fm2,i0,i1,j0,j1,k0,k1,filt)
c-----------------------------------------------------------------------
c apply time filters on prognostic variables
c-----------------------------------------------------------------------
c
      implicit real*8(a-h,o-z)
      dimension f(i0:i1,j0:j1,k0:k1),fm1(i0:i1,j0:j1,k0:k1)
     :   ,fm2(i0:i1,j0:j1,k0:k1)
c
      do 20 k=k0,k1
      do 20 j=j0,j1
      do 20 i=i0,i1
      fm1(i,j,k)=fm1(i,j,k)+filt*(fm2(i,j,k)+f(i,j,k)-2.*fm1(i,j,k))
20    continue
      return
      end




      subroutine fct3d(f,fm1,fm2,u,um1,v,vm1,wsig,wsigm1,nx,ny,ns
     :   ,dtxs,dtys,dtxy,dxys,dxdt,dydt,dt,ds0
     :   ,flx,fly,fls,sx,sy,ss,wk7,wk8)
c-----------------------------------------------------------------------
c advection using 3-d flux-corrected transport method
c after zalesak (jcp, 1979).
c
c with leapfg=.false., the higher order scheme is the 2nd order
c leapfrog-trapezoidal scheme described in the apendice of the
c paper.
c
c uses work spaces 1-6, but only 1-3 need to be kept when call to fctlim
c
c note that it is assumed that variable f is on the 000 grid.
c-----------------------------------------------------------------------
c
      implicit real*8(a-h,o-z)
      dimension f(0:nx+1,0:ny+1,0:ns+1),fm1(0:nx+1,0:ny+1,0:ns+1)
     :   ,fm2(0:nx+1,0:ny+1,0:ns+1)
      dimension u(0:nx+1,0:ny+1,0:ns+1),um1(0:nx+1,0:ny+1,0:ns+1)
      dimension v(0:nx+1,0:ny+1,0:ns+1),vm1(0:nx+1,0:ny+1,0:ns+1)
      dimension wsig(0:nx+1,0:ny+1,0:ns+1),wsigm1(0:nx+1,0:ny+1,0:ns+1)
      dimension flx(0:nx+1,0:ny+1,0:ns+1),fly(0:nx+1,0:ny+1,0:ns+1)
     :   ,fls(0:nx+1,0:ny+1,0:ns+1)
      dimension sx(0:nx+1,0:ny+1,0:ns+1),sy(0:nx+1,0:ny+1,0:ns+1)
     :   ,ss(0:nx+1,0:ny+1,0:ns+1)
      dimension wk7(0:nx+1,0:ny+1,0:ns+1),wk8(0:nx+1,0:ny+1,0:ns+1)
      dimension dtxs(0:ns+1),dtys(0:ns+1),dxys(0:ns+1),ds0(0:ns+1)
c
c      equivalence(wk1,flx),(wk2,fly),(wk3,fls)
c      equivalence(wk4,sx),(wk5,sy),(wk6,ss)
c
      logical leapfg
      parameter(leapfg=.false.)
      parameter(dif=1.0,dif8=0.125*dif)
c
      do 10 is=1,ns-1
      do 10 iy=2,ny
      do 10 ix=1,nx
      flx(ix,iy,is)=dtys(is)*um1(ix,iy,is)*(fm1(ix+1,iy,is)
     :   +fm1(ix,iy,is))
10    continue
c
      do 20 is=1,ns-1
      do 20 iy=1,ny
      do 20 ix=2,nx
      fly(ix,iy,is)=dtxs(is)*vm1(ix,iy,is)*(fm1(ix,iy+1,is)
     :   +fm1(ix,iy,is))
20    continue
c
      do 30 is=0,ns-1
      do 30 iy=2,ny
      do 30 ix=2,nx
      fls(ix,iy,is)=dtxy*wsigm1(ix,iy,is+1)*(fm1(ix,iy,is+1)
     :   +fm1(ix,iy,is))
30    continue
c
      do 40 is=1,ns-1
      do 40 iy=2,ny
      do 40 ix=2,nx
      f(ix,iy,is)=fm2(ix,iy,is)-(flx(ix,iy,is)-flx(ix-1,iy,is)
     :   +fly(ix,iy,is)-fly(ix,iy-1,is)+fls(ix,iy,is)-fls(ix,iy,is-1))
     :   /dxys(is)
40    continue
c
      if( leapfg ) then
         call extra3(nx+1,ny+1,ns+1,f,1,nx+1,1,ny+1,0,ns)
      else
         do 50 is=1,ns-1
         do 50 iy=2,ny
         do 50 ix=2,nx
         f(ix,iy,is)=(f(ix,iy,is)+fm1(ix,iy,is))*0.5
50       continue
c
         call extra3(nx+1,ny+1,ns+1,f,1,nx+1,1,ny+1,0,ns)
c
         do 55 is=1,ns-1
         do 55 iy=2,ny
         do 55 ix=1,nx
         flx(ix,iy,is)=dtys(is)*(um1(ix,iy,is)+u(ix,iy,is))
     :      *(f(ix,iy,is)+f(ix+1,iy,is))*0.25
55       continue
c
         do 56 is=1,ns-1
         do 56 iy=1,ny
         do 56 ix=2,nx
         fly(ix,iy,is)=dtxs(is)*(vm1(ix,iy,is)+v(ix,iy,is))
     :      *(f(ix,iy,is)+f(ix,iy+1,is))*0.25
56       continue
c
         do 57 is=0,ns-1
         do 57 iy=2,ny
         do 57 ix=2,nx
         fls(ix,iy,is)=dtxy*(wsigm1(ix,iy,is+1)+wsig(ix,iy,is+1))
     :      *(f(ix,iy,is)+f(ix,iy,is+1))*0.25
57       continue
      endif
c
c low order scheme: donnor cell.
c
      do 60 is=1,ns-1
      do 60 iy=2,ny
      do 60 ix=1,nx
      sx(ix,iy,is)=dtys(is)*(
     :   ((um1(ix,iy,is)+abs(um1(ix,iy,is)))*fm1(ix,iy,is)
     :   +(um1(ix,iy,is)-abs(um1(ix,iy,is)))*fm1(ix+1,iy,is)))
  !   :   -dif8*dxdt*(fm1(ix+1,iy,is)-fm1(ix,iy,is))
      flx(ix,iy,is)=flx(ix,iy,is)-sx(ix,iy,is)
60    continue
c
      do 70 is=1,ns-1
      do 70 iy=1,ny
      do 70 ix=2,nx
      sy(ix,iy,is)=dtxs(is)*(
     :   ((vm1(ix,iy,is)+abs(vm1(ix,iy,is)))*fm1(ix,iy,is)
     :   +(vm1(ix,iy,is)-abs(vm1(ix,iy,is)))*fm1(ix,iy+1,is)))
 !    :   -dif8*dydt*(fm1(ix,iy+1,is)-fm1(ix,iy,is))
      fly(ix,iy,is)=fly(ix,iy,is)-sy(ix,iy,is)
70    continue
c
      do 80 is=0,ns-1
      do 80 iy=2,ny
      do 80 ix=2,nx
      ss(ix,iy,is)=dtxy*(
     :   ((wsigm1(ix,iy,is+1)+abs(wsigm1(ix,iy,is+1)))*fm1(ix,iy,is)
     :   +(wsigm1(ix,iy,is+1)-abs(wsigm1(ix,iy,is+1)))*fm1(ix,iy,is+1)))
 !    :   -dif8*ds0(is+1)/dt*(fm1(ix,iy,is+1)-fm1(ix,iy,is))
      fls(ix,iy,is)=fls(ix,iy,is)-ss(ix,iy,is)
80    continue
c
      do 90 is=1,ns-1
      do 90 iy=2,ny
      do 90 ix=2,nx
      f(ix,iy,is)=fm2(ix,iy,is)-(sx(ix,iy,is)-sx(ix-1,iy,is)
     :   +sy(ix,iy,is)-sy(ix,iy-1,is)+ss(ix,iy,is)-ss(ix,iy,is-1))
     :   /dxys(is)
90    continue
c
      call extra3(nx+1,ny+1,ns+1,f,1,nx+1,1,ny+1,0,ns)
      call extra3(nx+1,ny+1,ns+1,f,0,nx+1,0,ny+1,0,ns+1)
c
      do 800 iy=0,ny+1
      do 800 ix=0,nx+1
      fm1(ix,iy,ns+1)=fm1(ix,iy,ns)
800   continue
c
      call fctlim(f,fm2,flx,fly,fls
     :   ,sx,sy,ss,wk7,wk8,wk8,wk8,wk8
     :   ,nx,ny,ns,dxys)
c
      do 51 is=1,ns-1
      do 51 iy=2,ny
      do 51 ix=2,nx
      f(ix,iy,is)=f(ix,iy,is)-(flx(ix,iy,is)-flx(ix-1,iy,is)
     :   +fly(ix,iy,is)-fly(ix,iy-1,is)+fls(ix,iy,is)-fls(ix,iy,is-1))
     :   /dxys(is)
51    continue
c
      call extra3(nx+1,ny+1,ns+1,f,1,nx+1,1,ny+1,0,ns)
c
      return
      end




      subroutine drag0z(nchn1,ax,ay
     :   ,nx,ny,ns,psuf,hsuf,dx,dy,dt,nstep
     :   ,exptit,iomtyp,hmount,xmount,ymount,hxwid,hywid,hexp,pa,ptop
     :   ,ioreft,thdat,zthdat,usdat,vsdat,pts0,ivarhs
     :   ,ix00,ix01,ix10,ix11,iy00,iy01,iy10,iy11,verbose,dragx,dragy
     :   ,drag2d)
c
      implicit real*8(a-h,o-z)
      integer verbose
      parameter(r=287.d0,cp=1005.d0,akapa=r/cp,g=9.8066d0)
      dimension psuf(0:nx+1,0:ny+1),hsuf(0:nx+1,0:ny+1)
      dimension thdat(0:*),zthdat(0:*),usdat(0:*),vsdat(0:*)
      character*8 exptit
      logical gridx0,gridy0,firstcall
c
      save dxarea,dyarea,dx2,dy2,dx8,dy8,dx64,dy64,dx256,dy256
     :   ,gridx0,gridy0
      save firstcall
      data firstcall/.true./


c
      if(firstcall) then
        firstcall=.false.
        dx2=0.5*dx
        dx8=dx/8.
        dx64=dx/64.
        dx256=dx/256.
c
        dy2=0.5*dy
        dy8=dy/8.
        dy64=dy/64.
        dy256=dy/256.
c
        u0=usdat(0)
        v0=vsdat(0)
        th=thdat(0)
        if(ioreft.eq.3) then
          bv2=th
          t0=pts0*(pa/1.e5)**akapa
        else
          t0=th
          bv2=g/t0*(thdat(2)-thdat(1))/(zthdat(2)-zthdat(1))
          if(ioreft.eq.2) bv2=bv2+g*g/(cp*t0)
        endif
        bruntv=sqrt(max(0.d0,bv2))
c
        write(nchn1,*)'ax,ay,nx,ny,ns',ax,ay,nx,ny,ns
        dgarea=4.*ax*ay
        dxarea=dgarea
        dyarea=dgarea
c
        x0=xmount-ax
        x1=xmount+ax
        y0=ymount-ay
        y1=ymount+ay
c
        if(abs(x0-nint(x0/dx)*dx).lt..1) then
          gridx0=.false.
          ix01=int(x0/dx+1.1)
          ix11=int(x1/dx+1.1)
          ix00=ix01+1
          ix10=ix11
        else
          gridx0=.true.
          ix00=int(x0/dx+1.6)
          ix10=int(x1/dx+1.6)
          ix01=ix00
          ix11=ix10-1
        endif
c
        if(abs(y0-nint(y0/dy)*dy).lt..1) then
          gridy0=.false.
!indices for y grid limits
          iy01=int(y0/dy+1.1)
          iy11=int(y1/dy+1.1)
!indices for y grid limits
          iy00=iy01+1
          iy10=iy11
        else
          gridy0=.true.
          iy00=int(y0/dy+1.6)
          iy10=int(y1/dy+1.6)
          iy01=iy00
          iy11=iy10-1
        endif
c
        write(nchn1,1011) exptit,ivarhs
        write(nchn1,1010) x0/1000.,x1/1000.,y0/1000.,y1/1000.
        write(nchn1,*) dxarea,dyarea,ns,dx,dy
     :    ,ix00,ix01,ix10,ix11
     :    ,iy00,iy01,iy10,iy11
        write(nchn1,*)
     :    'dt,iomtyp,hmount,hxwid,hywid,u0,bruntv,pa,t0,hexp,v0'
        write(nchn1,1002) dt,iomtyp,hmount,hxwid,hywid,u0,bruntv,pa,t0
     :    ,hexp,v0
        write(nchn1,*) 'gridx0=',gridx0,ix00,ix01,ix10,ix11
        write(nchn1,*) 'gridy0=',gridy0,iy00,iy01,iy10,iy11
        write(50,*)
      endif


      dragx=0.


      do iy=iy01+1,iy11           !ok for .not.gridy0
      do ix=ix00,ix10-1
        dragx=dragx+(psuf(ix,iy)+psuf(ix+1,iy))
     :     *(hsuf(ix+1,iy)-hsuf(ix,iy))*dy2
      enddo
      enddo


      if(verbose.ge.3) write(*,*) 'Dragx(0):',dragx
      dragx0=dragx


      if(.not.gridx0) then
        do iy=iy01+1,iy11
          dragx=dragx+(psuf(ix00-1,iy)+3.*psuf(ix00,iy))
     :      *(hsuf(ix00,iy)-hsuf(ix00-1,iy))*dy8
          dragx=dragx+(3.*psuf(ix10,iy)+psuf(ix10+1,iy))
     :      *(hsuf(ix10+1,iy)-hsuf(ix10,iy))*dy8
        enddo
      endif
c
      if(verbose.ge.3) write(*,*) 'Dragx(1):',dragx
      dragx1=dragx


      if(gridy0) then
        do ix=ix00,ix10-1
          dragx=dragx+(3.*psuf(ix,iy00)+psuf(ix,iy00+1)
     :      +3.*psuf(ix+1,iy00)+psuf(ix+1,iy00+1))
     :      *(3.*hsuf(ix+1,iy00)+hsuf(ix+1,iy00+1)
     :      -3.*hsuf(ix,iy00)-hsuf(ix,iy00+1))*dy64
          dragx=dragx+(3.*psuf(ix,iy10)+psuf(ix,iy10-1)
     :      +3.*psuf(ix+1,iy10)+psuf(ix+1,iy10-1))
     :      *(3.*hsuf(ix+1,iy10)+hsuf(ix+1,iy10-1)
     :      -3.*hsuf(ix,iy10)-hsuf(ix,iy10-1))*dy64
        enddo
      endif
      if(verbose.ge.3) write(*,*) 'Dragx(2):',dragx
      dragx2=dragx
c
      if(.not.gridx0 .and. gridy0) then
         dragx=dragx+(9.*psuf(ix00,iy00)+3.*psuf(ix00,iy00+1)
     :      +3.*psuf(ix00-1,iy00)+psuf(ix00-1,iy00+1))
     :      *(3.*hsuf(ix00,iy00)+hsuf(ix00,iy00+1)-3.*hsuf(ix00-1,iy00)
     :      -hsuf(ix00-1,iy00+1))*dy256
         dragx=dragx+(9.*psuf(ix10,iy00)+3.*psuf(ix10,iy00+1)
     :      +3.*psuf(ix10+1,iy00)+psuf(ix10+1,iy00+1))
     :      *(3.*hsuf(ix10+1,iy00)+hsuf(ix10+1,iy00+1)
     :      -3.*hsuf(ix10,iy00)-hsuf(ix10,iy00+1))*dy256
         dragx=dragx+(9.*psuf(ix00,iy10)+3.*psuf(ix00,iy10-1)
     :      +3.*psuf(ix00-1,iy10)+psuf(ix00-1,iy10-1))
     :      *(3.*hsuf(ix00,iy10)+hsuf(ix00,iy10-1)-3.*hsuf(ix00-1,iy10)
     :      -hsuf(ix00-1,iy10-1))*dy256
         dragx=dragx+(9.*psuf(ix10,iy10)+3.*psuf(ix10,iy10-1)
     :      +3.*psuf(ix10+1,iy10)+psuf(ix10+1,iy10-1))
     :      *(3.*hsuf(ix10+1,iy10)+hsuf(ix10+1,iy10-1)
     :      -3.*hsuf(ix10,iy10)-hsuf(ix10,iy10-1))*dy256
c         arx=arx+dy
      endif
      if(verbose.ge.3) write(*,*) 'Dragx(3):',dragx
      dragx3=dragx
c
      dragy=0.
c      ary=0.
c
      do 20 iy=iy00,iy10-1
      do 20 ix=ix01+1,ix11
        dragy=dragy+(psuf(ix,iy)+psuf(ix,iy+1))
     :    *(hsuf(ix,iy+1)-hsuf(ix,iy))*dx2
c      ary=ary+dx
20    continue
c
c      write(50,*) 'ary=',ary
c
c      write(*,*) '20 cont'
c
      if(.not.gridy0) then
         do 21 ix=ix01+1,ix11
         dragy=dragy+(psuf(ix,iy00-1)+3.*psuf(ix,iy00))
     :      *(hsuf(ix,iy00)-hsuf(ix,iy00-1))*dx8
         dragy=dragy+(3.*psuf(ix,iy10)+psuf(ix,iy10+1))
     :      *(hsuf(ix,iy10+1)-hsuf(ix,iy10))*dx8
c         ary=ary+dx
21       continue
      endif
c
      if(gridx0) then
         do 22 iy=iy00,iy10-1
         dragy=dragy+(3.*psuf(ix00,iy)+psuf(ix00+1,iy)
     :      +3.*psuf(ix00,iy+1)+psuf(ix00+1,iy+1))
     :      *(3.*hsuf(ix00,iy+1)+hsuf(ix00+1,iy+1)
     :      -3.*hsuf(ix00,iy)-hsuf(ix00+1,iy))*dx64
         dragy=dragy+(3.*psuf(ix10,iy)+psuf(ix10-1,iy)
     :      +3.*psuf(ix10,iy+1)+psuf(ix10-1,iy+1))
     :      *(3.*hsuf(ix10,iy+1)+hsuf(ix10-1,iy+1)
     :      -3.*hsuf(ix10,iy)-hsuf(ix10-1,iy))*dx64
c         ary=ary+dx
22       continue
      endif
c
      if(.not.gridy0 .and. gridx0) then
         dragy=dragy+(9.*psuf(ix00,iy00)+3.*psuf(ix00+1,iy00)
     :      +3.*psuf(ix00,iy00-1)+psuf(ix00+1,iy00-1))
     :      *(3.*hsuf(ix00,iy00)+hsuf(ix00+1,iy00)-3.*hsuf(ix00,iy00-1)
     :      -hsuf(ix00+1,iy00-1))*dx256
         dragy=dragy+(9.*psuf(ix00,iy10)+3.*psuf(ix00+1,iy10)
     :      +3.*psuf(ix00,iy10+1)+psuf(ix00+1,iy10+1))
     :      *(3.*hsuf(ix00,iy10+1)+hsuf(ix00+1,iy10+1)
     :      -3.*hsuf(ix00,iy10)-hsuf(ix00+1,iy10))*dx256
         dragy=dragy+(9.*psuf(ix10,iy00)+3.*psuf(ix10-1,iy00)
     :      +3.*psuf(ix10,iy00-1)+psuf(ix10-1,iy00-1))
     :      *(3.*hsuf(ix10,iy00)+hsuf(ix10-1,iy00)-3.*hsuf(ix10,iy00-1)
     :      -hsuf(ix10-1,iy00-1))*dx256
         dragy=dragy+(9.*psuf(ix10,iy10)+3.*psuf(ix10-1,iy10)
     :      +3.*psuf(ix10,iy10+1)+psuf(ix10-1,iy10+1))
     :      *(3.*hsuf(ix10,iy10+1)+hsuf(ix10-1,iy10+1)
     :      -3.*hsuf(ix10,iy10)-hsuf(ix10-1,iy10))*dx256
c         ary=ary+dx
      endif


      drag2d=0.
      iy=ny/2+1
      do ix=ix00,ix10-1
        drag2d=drag2d+0.5*(psuf(ix,iy)+psuf(ix+1,iy))
     :     *(hsuf(ix+1,iy)-hsuf(ix,iy))
      enddo
c
c      write(*,*) 'dragx:',dragx/(arx*dx),arx*dx
c      write(*,*) 'dragy:',dragy/(ary*dy),ary*dy
      write(50,'(i5,3e15.7)') nstep,dragx,dragy,drag2d
c     write(50,'(i5,6e15.7)') nstep,dragx,dragy,dragx0,dragx1,dragx2
c    :  ,dragx3
      if(verbose.ge.1) then
        write(*,1001) nstep,dragx,dragy,drag2d
      endif
      return
1000  format(1x,'nstep',t15,'dragx',t30,'dragy')
c1001  format(1x,i5,t10,e14.7,t25,e14.7,f7.1)
1001  format(1x,i5,3e15.7)
1002  format(1x,f5.0,i2,f8.1,2f9.0,f7.2,f10.5,f9.1,f8.2,f5.1,f7.2)
1010  format(1x,'drag for x0=',f8.1,' x1=',f8.1,' y0=',f8.1,' y1=',f8.1
     :   ,' units: n/m**2')
1011  format('''',a8,'''',i2)
      end
c
c
c
      subroutine draguv(ix0,ix1,iy0,iy1,is0,is1
     :   ,nx,ny,ns,dx,dy,tems,psuf,u,v,w,nstep,us,vs,ptop,sigma0)
c
c calculates vertical momentum flux (no rotation).
c
c
      implicit real*8(a-h,o-z)
      parameter(r=287.d0,cp=1005.d0,akapa=r/cp)
      dimension tems(0:nx+1,0:ny+1,0:ns+1)
     :   ,psuf(0:nx+1,0:ny+1)
     :   ,u(0:nx+1,0:ny+1,0:ns+1),v(0:nx+1,0:ny+1,0:ns+1)
     :   ,us(0:nx+1,0:ny+1,0:ns+1),vs(0:nx+1,0:ny+1,0:ns+1)
     :   ,w(0:nx+1,0:ny+1,0:ns+1),sigma0(0:ns+1)
c      parameter(maxns=1000)
c      dimension duw(0:maxns),dvw(0:maxns)
      real(kind(0.d0)),allocatable,dimension(:):: duw,dvw

      if(.not.allocated(duw)) then
        allocate(duw(0:ns))
        allocate(dvw(0:ns))
      endif

c
      dxy025=0.25*dx*dy
c
      write(51,1020) nstep,ix0,ix1,iy0,iy1,is0,is1,dx,dy,ns
c
      do 40 is=is0,is1
c
      duw(is)=0.
      dvw(is)=0.
c
      do 30 iy=iy0+1,iy1
      do 30 ix=ix0+1,ix1
      ro=(sigma0(is)*psuf(ix,iy)+(1-sigma0(is))*ptop)
     :      /(r*tems(ix,iy,is))
      duw(is)=duw(is)-ro*w(ix,iy,is)
     :   *(u(ix,iy,is)+u(ix,iy,is-1)+u(ix-1,iy,is)
     :   +u(ix-1,iy,is-1)
     :   -us(ix,iy,is)-us(ix,iy,is-1)-us(ix-1,iy,is)
     :   -us(ix-1,iy,is-1))*dxy025
      dvw(is)=dvw(is)-ro*w(ix,iy,is)
     :   *(v(ix,iy,is)+v(ix,iy,is-1)+v(ix,iy-1,is)
     :   +v(ix,iy-1,is-1)
     :   -vs(ix,iy,is)-vs(ix,iy,is-1)-vs(ix,iy-1,is)
     :   -vs(ix,iy-1,is-1))*dxy025
30    continue
40    continue
c
      duwmmm=0.
      dvwmmm=0.
      do 90 is=is0,is1
      duwmmm=max(duwmmm,abs(duw(is)))
      dvwmmm=max(dvwmmm,abs(dvw(is)))
90    continue
c
      if(duwmmm.eq.0.) duwmmm=1.
      if(dvwmmm.eq.0.) dvwmmm=1.
c
      write(51,1001) duwmmm,dvwmmm
      do 100 is=is0,is1
      duw(is)=duw(is)/duwmmm
      dvw(is)=dvw(is)/dvwmmm
100   continue
      write(51,1030) (duw(is),is=is0,is1)
      write(51,1030) (dvw(is),is=is0,is1)
c
      return
1001  format(1x,2e15.7)
1020  format(i6,6i4,2f10.1,i4)
1030  format(1x,7f09.5)
      end

      subroutine drguv2(ix0,ix1,iy0,iy1,is0,is1
     :   ,nx,ny,ns,dx,dy,tems,psuf,u,v,w,nstep,ptop,sigma0)
c
c calculates vertical momentum flux (no rotation).
c
      implicit real*8(a-h,o-z)
      parameter(r=287.,cp=1005.,akapa=r/cp)
c
      dimension tems(0:nx+1,0:ny+1,0:ns+1)
     :   ,psuf(0:nx+1,0:ny+1)
     :   ,u(0:nx+1,0:ny+1,0:ns+1),v(0:nx+1,0:ny+1,0:ns+1)
     :   ,w(0:nx+1,0:ny+1,0:ns+1),sigma0(0:ns+1)
c      parameter(maxns=1000)
c      dimension duw(0:maxns),dvw(0:maxns)
      real(kind(0.d0)),allocatable,dimension(:):: duw,dvw

      if(.not.allocated(duw)) then
        allocate(duw(0:ns))
        allocate(dvw(0:ns))
      endif
c
      dxy=dx*dy
c
      write(51,1020) nstep,ix0,ix1,iy0,iy1,is0,is1,dx,dy,ns
c
      do 40 is=is0,is1
c
      duw(is)=0.
      dvw(is)=0.
c
      umed=0.
      vmed=0.
      nn=(nx-2)*(ny-2)
      do 25 iy=2,ny-1
      do 25 ix=2,nx-1
      umed=umed+u(ix,iy,is)
      vmed=vmed+v(ix,iy,is)
25    continue
      umed=umed/nn
      vmed=vmed/nn
c
      do 30 iy=iy0+1,iy1
      do 30 ix=ix0+1,ix1
      ro=(sigma0(is)*(psuf(ix,iy)-ptop)+ptop)/(r*tems(ix,iy,is))
      duw(is)=duw(is)-ro*w(ix,iy,is)
     :   *(0.25*(u(ix,iy,is)+u(ix,iy,is-1)+u(ix-1,iy,is)
     :   +u(ix-1,iy,is-1))-umed)*dxy
      dvw(is)=dvw(is)-ro*w(ix,iy,is)
     :   *(0.25*(v(ix,iy,is)+v(ix,iy,is-1)+v(ix,iy-1,is)
     :   +v(ix,iy-1,is-1))-vmed)*dxy
30    continue
40    continue
c
      duwmmm=0.
      dvwmmm=0.
      do 90 is=is0,is1
      duwmmm=max(duwmmm,abs(duw(is)))
      dvwmmm=max(dvwmmm,abs(dvw(is)))
90    continue
c
      if(duwmmm.eq.0.) duwmmm=1.
      if(dvwmmm.eq.0.) dvwmmm=1.
c
      write(51,1001) duwmmm,dvwmmm
      do 100 is=is0,is1
      duw(is)=duw(is)/duwmmm
      dvw(is)=dvw(is)/dvwmmm
100   continue
      write(51,1030) (duw(is),is=is0,is1)
      write(51,1030) (dvw(is),is=is0,is1)
c
      return
1001  format(1x,2e15.7)
1020  format(i6,6i4,2f10.1,i4)
1030  format(1x,7f09.5)
      end
c
c
c
      subroutine momflu(ix0,ix1,iy0,iy1,is0,is1,isover,z,zlev,nlev
     :   ,nx,ny,ns,dx,dy,tems,psuf,u,v,w
     :   ,nstep,us,vs,ptop,sigma0)
c
c calculates vertical momentum flux (no rotation).
c
      implicit real*8(a-h,o-z)
      parameter(r=287.,cp=1005.,akapa=r/cp,g=9.8066)
      dimension zlev(nlev)
c
      dimension isover(0:nx+1,0:ny+1,*)
      dimension tems(0:nx+1,0:ny+1,0:ns+1),psuf(0:nx+1,0:ny+1)
     :   ,u(0:nx+1,0:ny+1,0:ns+1),v(0:nx+1,0:ny+1,0:ns+1)
     :   ,us(0:nx+1,0:ny+1,0:ns+1),vs(0:nx+1,0:ny+1,0:ns+1)
     :   ,w(0:nx+1,0:ny+1,0:ns+1),z(0:nx+1,0:ny+1,0:ns+1)
     :   ,sigma0(0:ns+1)
c
      dxy=dx*dy
c
      write(53,1020) nstep,ix0,ix1,iy0,iy1,is0,is1,dx,dy,nlev
      write(53,1000)
c
      do 199 ilev=1,nlev
      do 199 iy=0,ny+1
      do 199 ix=0,nx+1
      isover(ix,iy,ilev)=0
199   continue
c
      do 200 ilev=1,nlev
      do 210 is=0,ns-1
      do 220 iy=iy0+1,iy1
      do 220 ix=ix0+1,ix1
      if(zlev(ilev).lt.z(ix,iy,is) .and. zlev(ilev).gt.z(ix,iy,is+1))
     :   isover(ix,iy,ilev)=is
220   continue
210   continue
200   continue
c
      do 40 ilev=1,nlev
c
      duw=0.
      dvw=0.
c
      do 30 iy=iy0+1,iy1
      do 30 ix=ix0+1,ix1
      isov=isover(ix,iy,ilev)
      rov=(sigma0(isov)*(psuf(ix,iy)-ptop)+ptop)
     :      /(r*tems(ix,iy,isov))
      rovp1=(sigma0(isov+1)*(psuf(ix,iy)-ptop)+ptop)
     :      /(r*tems(ix,iy,isov+1))
      delta=(zlev(ilev)-z(ix,iy,isov+1))/(z(ix,iy,isov)-z(ix,iy,isov+1))
      si=sign(1.d0,delta-0.5)
      isi=iabs(nint((si-1.)/2))
      wa=(delta-si*0.5)*w(ix,iy,isov+isi)
     :   +(1.-delta+si*0.5)*w(ix,iy,isov+1+isi)
      duw=duw-(delta*(rov
     :   *0.5*(u(ix,iy,isov)+u(ix-1,iy,isov)
     :   -us(ix,iy,isov)-us(ix-1,iy,isov)))
     :   +(1.-delta)*(rovp1
     :   *0.5*(u(ix,iy,isov+1)+u(ix-1,iy,isov+1)
     :   -us(ix,iy,isov+1)-us(ix-1,iy,isov+1))))*wa
      dvw=dvw-(delta*(rov
     :   *0.5*(v(ix,iy,isov)+v(ix,iy-1,isov)
     :   -vs(ix,iy,isov)-vs(ix,iy-1,isov)))
     :   +(1.-delta)*(rovp1
     :   *0.5*(v(ix,iy,isov+1)+v(ix,iy-1,isov+1)
     :   -vs(ix,iy,isov+1)-vs(ix,iy-1,isov+1))))*wa
30    continue
c
      write(53,1001) zlev(ilev),duw*dxy,dvw*dxy
40    continue
c
      return
1000  format(1x,'z',t15,'duw',t30,'dvw')
1001  format(1x,f9.1,2e15.7)
1020  format(i6,6i4,2f10.1,i4)
      end




      subroutine energy(ix0,ix1,iy0,iy1,is0,is1
     :   ,nx,ny,ns
     :   ,u,v,w,pt,pts,phi,phis,pp,pp10,pp01
     :   ,sigma0,sigma1,ds1x,ds1y,dxys,dxy,ds02a,ptop,nchn1,nchwr)
c
c calculates:
c
c kinetic and potential energy in a given box
c    ektot and eku2,ekv2,ekw2
c
c conversion from perturbation potential energy to kinetic energy:
c    cptk
c
c conversion from potential energy of the reference state to the potential
c energy of the perturbation:
c     cptspt
c
c boundary fluxes of kinetic energy and of pert.pot.energy:
c     flukx,fluky,fluphx,fluphy
c
      implicit real*8(a-h,o-z)
      dimension u(0:nx+1,0:ny+1,0:ns+1),w(0:nx+1,0:ny+1,0:ns+1)
     :   ,v(0:nx+1,0:ny+1,0:ns+1)
     :   ,pt(0:nx+1,0:ny+1,0:ns+1),pts(0:nx+1,0:ny+1,0:ns+1)
     :   ,phi(0:nx+1,0:ny+1,0:ns+1),phis(0:nx+1,0:ny+1,0:ns+1)
     :   ,pp(0:nx+2,0:ny+2),pp10(0:nx+1,0:ny+1)
     :   ,pp01(0:nx+1,0:ny+1)
      dimension sigma0(0:ns+1),sigma1(0:ns+1),ds1x(0:ns+1),ds1y(0:ns+1)
     :   ,dxys(0:ns+1),ds02a(0:ns+1)
      integer nchwr
c
      parameter(r=287.,cp=1005.,g=9.8066,akapa=r/cp)
      parameter(p00=1.e5)
c
c
      eku2=0.
      ekv2=0.
      ekw2=0.
c
      flu2xl=0.
      flu2yl=0.
      flv2xl=0.
      flv2yl=0.
      flw2xl=0.
      flw2yl=0.
      flu2xr=0.
      flu2yr=0.
      flv2xr=0.
      flv2yr=0.
      flw2xr=0.
      flw2yr=0.
c
      flphxl=0.
      flphyl=0.
      flphxr=0.
      flphyr=0.
c
      ept=0.
      cptk=0.
      cptspt=0.
c
      do 10 is=is0,is1-1
      do 10 iy=iy0+1,iy1
      do 10 ix=ix0+1,ix1
      eku2=eku2+u(ix,iy,is)*pp10(ix,iy)*u(ix,iy,is)*dxys(is)
      ekv2=ekv2+v(ix,iy,is)*pp01(ix,iy)*v(ix,iy,is)*dxys(is)
      ekw2=ekw2+pp(ix,iy)*w(ix,iy,is)*w(ix,iy,is)*dxys(is)
10    continue
c
      do 20 is=is0,is1-1
      do 20 iy=iy0+1,iy1
      eku2=eku2+(u(ix0,iy,is)*pp10(ix0,iy)*u(ix0,iy,is)
     :   +u(ix1,iy,is)*pp10(ix1,iy)*u(ix1,iy,is))*dxys(is)*0.5
20    continue
c
      do 30 is=is0,is1-1
      do 30 ix=ix0+1,ix1
      ekv2=ekv2+(v(ix,iy0,is)*pp01(ix,iy0)*v(ix,iy0,is)
     :   +v(ix,iy1,is)*pp01(ix,iy1)*v(ix,iy1,is))*dxys(is)*0.5
30    continue
c
      do 40 iy=iy0+1,iy1
      do 40 ix=ix0+1,ix1
      ekw2=ekw2+pp(ix,iy)*
     :   (w(ix,iy,is0)*w(ix,iy,is0)*(sigma0(is0)-sigma1(is0))
     :   +w(ix,iy,is1)*w(ix,iy,is1)*(sigma1(is1)-sigma0(is1-1)))
     :   *dxy
40    continue
c
      eku2=eku2/(2.*g)
      ekv2=ekv2/(2.*g)
      ekw2=ekw2/(2.*g)
c
      ektot=eku2+ekv2+ekw2
c
c total potential energy (perturbation only)
c and conversions:
c
      do 110 is=is0,is1-1
      do 110 iy=iy0+1,iy1
      do 110 ix=ix0+1,ix1
      p=sigma0(is)*pp(ix,iy)+ptop
      ept=ept+pt(ix,iy,is)*pp(ix,iy)*p**akapa*dxys(is)
      cptk=cptk+pt(ix,iy,is)*pp(ix,iy)/pts(ix,iy,is)
     :   *(w(ix,iy,is)+w(ix,iy,is+1))*dxys(is)
      cptspt=cptspt-p*(pts(ix,iy,is+1)-pts(ix,iy,is-1))/ds02a(is)
     :   /pts(ix,iy,is)*(w(ix,iy,is)+w(ix,iy,is+1))*dxys(is)
110   continue
c
      ept=ept*cp/(g*p00**akapa)
      cptk=cptk*0.5
      cptspt=cptspt*0.5*cp/r
c


      do 200 is=is0,is1-1
      do 200 iy=iy0+1,iy1
      flu2xl=flu2xl+u(ix0,iy,is)*pp10(ix0,iy)*u(ix0,iy,is)*u(ix0,iy,is)
     :   *ds1y(is)
      flu2xr=flu2xr+u(ix1,iy,is)*pp10(ix1,iy)*u(ix1,iy,is)*u(ix1,iy,is)
     :   *ds1y(is)
      flv2xl=flv2xl
     :   +(u(ix0,iy,is)*pp10(ix0,iy)*(v(ix0,iy,is)
     :   +v(ix0+1,iy,is)+v(ix0,iy-1,is)+v(ix0+1,iy-1,is))**2)
     :   *ds1y(is)
      flv2xr=flv2xr
     :   +(u(ix1,iy,is)*pp10(ix1,iy)*(v(ix1,iy,is)
     :   +v(ix1+1,iy,is)+v(ix1,iy-1,is)+v(ix1+1,iy-1,is))**2)
     :   *ds1y(is)
      flw2xl=flw2xl
     :   +(u(ix0,iy,is)*pp10(ix0,iy)*(w(ix0,iy,is)
     :   +w(ix0+1,iy,is)+w(ix0,iy,is+1)+w(ix0+1,iy,is+1))**2)
     :   *ds1y(is)
      flw2xr=flw2xr
     :   +(u(ix1,iy,is)*pp10(ix1,iy)*(w(ix1,iy,is)
     :   +w(ix1+1,iy,is)+w(ix1,iy,is+1)+w(ix1+1,iy,is+1))**2)
     :   *ds1y(is)
      flphxl=flphxl
     :   +u(ix0,iy,is)*pp10(ix0,iy)*(phi(ix0,iy,is)+phi(ix0+1,iy,is))
     :   *ds1y(is)
      flphxr=flphxr
     :   +u(ix1,iy,is)*pp10(ix1,iy)*(phi(ix1,iy,is)+phi(ix1+1,iy,is))
     :   *ds1y(is)
200   continue
      flu2xl=flu2xl/(2.*g)
      flv2xl=flv2xl/(32.*g)
      flw2xl=flw2xl/(32.*g)
      flphxl=flphxl/(2.*g)
      flu2xr=flu2xr/(2.*g)
      flv2xr=flv2xr/(32.*g)
      flw2xr=flw2xr/(32.*g)
      flphxr=flphxr/(2.*g)
c
      do 210 is=is0,is1-1
      do 210 ix=ix0+1,ix1
      flu2yl=flu2yl
     :   +(v(ix,iy0,is)*pp01(ix,iy0)*(u(ix,iy0,is)
     :   +u(ix,iy0+1,is)+u(ix-1,iy0,is)+u(ix-1,iy0+1,is))**2)
     :   *ds1x(is)
      flu2yr=flu2yr
     :   +(v(ix,iy1,is)*pp01(ix,iy1)*(u(ix,iy1,is)
     :   +u(ix,iy1+1,is)+u(ix-1,iy1,is)+u(ix-1,iy1-1,is))**2)
     :   *ds1x(is)
      flv2yl=flv2yl
     :   +v(ix,iy0,is)*pp01(ix,iy0)*v(ix,iy0,is)*v(ix,iy0,is)*ds1x(is)
      flv2yr=flv2yr
     :   +v(ix,iy1,is)*pp01(ix,iy1)*v(ix,iy1,is)*v(ix,iy,is)*ds1x(is)
      flw2yl=flw2yl
     :   +(v(ix,iy0,is)*pp01(ix,iy0)*(w(ix,iy0,is)
     :   +w(ix,iy0+1,is)+w(ix,iy0,is+1)+w(ix,iy0+1,is+1))**2)
     :   *ds1x(is)
      flw2yr=flw2yr
     :   +(v(ix,iy1,is)*pp01(ix,iy1)*(w(ix,iy1,is)
     :   +w(ix,iy1+1,is)+w(ix,iy1,is+1)+w(ix,iy1+1,is+1))**2)
     :   *ds1x(is)
      flphyl=flphyl
     :   +v(ix,iy0,is)*pp01(ix,iy0)*(phi(ix,iy0,is)+phi(ix,iy0+1,is))
     :   *ds1x(is)
      flphyr=flphyr
     :   +v(ix,iy1,is)*pp01(ix,iy1)*(phi(ix,iy1,is)+phi(ix,iy1+1,is))
     :   *ds1x(is)
210   continue
      flu2yl=flu2yl/(32.*g)
      flv2yl=flv2yl/(2.*g)
      flw2yl=flw2yl/(32.*g)
      flphyl=flphyl/(2.*g)
      flu2yr=flu2yr/(32.*g)
      flv2yr=flv2yr/(2.*g)
      flw2yr=flw2yr/(32.*g)
      flphyr=flphyr/(2.*g)
c
      fluk=flu2xl-flu2xr+flv2xl-flv2xr+flw2xl-flw2xr
     :   +flu2yl-flu2yr+flv2yl-flv2yr+flw2yl-flw2yr
c
      fluphi=flphxl-flphxr+flphyl-flphyr
c
      etende=fluk+fluphi+cptspt
c
      write(nchn1,1001)
      write(nchn1,1010) ektot,ept,ektot+ept
      write(nchn1,1002)
      write(nchn1,1010) eku2,ekv2,ekw2
      write(nchn1,1003)
      write(nchn1,1010) fluk,fluphi,etende
      write(nchn1,1004)
      write(nchn1,1010) cptk,cptspt
      write(nchn1,1005)
      write(nchn1,1010) flu2xl,flu2xr,flu2yl,flu2yr
      write(nchn1,1006)
      write(nchn1,1010) flv2xl,flv2xr,flv2yl,flv2yr
      write(nchn1,1007)
      write(nchn1,1010) flw2xl,flw2xr,flw2yl,flw2yr
      write(nchn1,1008)
      write(nchn1,1010) flphxl,flphxr,flphyl,flphyr
      write(nchwr,1067)ektot,ept,ektot+ept,cptspt
c
      return
1001  format(//1x,t13,'kinetic',t38,'perturb.pot.',t63,'total')
1002  format(1x,t13,'eku2',t38,'ekv2',t63,'ekw2')
1003  format(1x,t13,'fluk',t38,'fluphi',t63,'etende')
1004  format(1x,t13,'cptk',t38,'cptspt')
1005  format(1x,t13,'flu2xl',t38,'flu2xr',t63,'flu2yl',t88,'flu2yr')
1006  format(1x,t13,'flv2xl',t38,'flv2xr',t63,'flv2yl',t88,'flv2yr')
1007  format(1x,t13,'flw2xl',t38,'flw2xr',t63,'flw2yl',t88,'flw2yr')
1008  format(1x,t13,'flphxl',t38,'flphxr',t63,'flphyl',t88,'flphyr')
1010  format(1x,4e23.16)
1067  format(e23.16,3h   ,e23.16,3h   ,e23.16,3h   ,e23.16,3h   )
      end


      subroutine zero(x,lmn)
c
c puts a vector to zero
c
      implicit real*8(a-h,o-z)
      dimension x(lmn)
c
      do 10 ijk=1,lmn
      x(ijk)=0.
10    continue
      return
      end
c
c ***
c
      subroutine sintrx(f,w,sines,n,lot,isign)
c
c slow fourier transform
c
c   for j=1,lot
c       2/n*sum(i=1,n-1) f(i,j)*sin(i*pi/n)
c
c  in the inverse case (isign=-1) the factor 2/n is omitted.
c
      implicit real*8(a-h,o-z)
      dimension f(-1:n,-1:lot+1)
      dimension w(n,lot),sines(n-1,n-1)
c
      ton=2./float(n)
c
      do 10 j=1,lot
      do 10 i=1,n-1
      w(i,j)=0.
10    continue
c
      do 20 j=1,lot
      do 20 k=1,n-1
      do 20 i=1,n-1
      w(k,j)=w(k,j)+f(i,j)*sines(k,i)
20    continue
c
      if(isign.eq.1) then
         do 30 j=1,lot
         do 30 k=1,n-1
         f(k,j)=ton*w(k,j)
30       continue
      else
         do 40 j=1,lot
         do 40 k=1,n-1
         f(k,j)=w(k,j)
40       continue
      endif
c
      return
      end
c
c ***
c
      subroutine instra(sines,sine,n)
c
      implicit real*8(a-h,o-z)
      dimension sines(n-1,n-1)
      double precision sine(0:n-1),dpi
c
      dpi=4.d0*datan(1.d0)
c
      sine(0)=0.d0
      do 1 i=1,n-1
      sine(i)=dsin(dble(i)*dpi/dble(n))
1     continue
c
      nn=n+n
      nh=n/2+1
c
      do 10 k=1,n-1
      do 10 i=1,n-1
      ik=i*k
      nc=mod(ik,nn)
      ih=nc/n
      iih=mod(nc,n)
      iqh=iih/nh
      iq=iqh+1+2*ih
      iii=iqh*n+(-1)**iqh*iih
c      write(10,1000) i,k,ik,nc,ih,iih,iqh,iq,iii
c1000  format(10i4)
      sines(i,k)=(-1)**ih*sine(iii)
10    continue
c
c      call wrigar(sines,1,n-1,1,n-1,0,0,1,n-1,1,n-1,0,0,'sines   ',0.d0,1)
c
      return
      end
c
c ***
c
      subroutine scfft(f,w,trigs,ifax,sipco,simco,n,lot,isign)
c
c cosine transform on a staggered grid.
c following wilhelmson and ericksen (j.c.p. 25,319-331)
c
c performs:
c   for j=1,lot
c       sum(i=1,n) f(i,j)*cos((i-0.5)*(i-1)*pi/n)
c
c note the dimensions of f: they are appropriate for nh3d (phisuf solution)
c
c p.m. 1988/nov
c-----------------------------------------------------------------------
c
      implicit real*8(a-h,o-z)
      dimension f(-1:n+1,-1:lot+1)
      dimension w(n+2,lot)
      dimension sipco(n),simco(n)
      dimension trigs(*)
      dimension ifax(10)
c
      if(isign.eq.1) then
         do 10 j=1,lot
         w(1,j)=f(1,j)
         w(2,j)=0.
         w(n+1,j)=f(n,j)
         w(n+2,j)=0.
10       continue
         do 20 ii=2,n-2,2
         do 20 j=1,lot
         w(ii+1,j)=0.5*(f(ii,j)+f(ii+1,j))
         w(ii+2,j)=-0.5*(f(ii,j)-f(ii+1,j))
20       continue
c
         call fft991(w,f,trigs,ifax,1,n+2,n,lot,+1)
c
         do 30 j=1,lot
         f(1,j)=w(1,j)
30       continue
         do 40 i=2,n
         do 40 j=1,lot
         f(i,j)=0.5*(sipco(i)*w(i,j)-simco(i)*w(n-i+2,j))
40       continue
      else
         do 110 j=1,lot
         w(1,j)=f(1,j)
110      continue
         do 120 i=2,n
         do 120 j=1,lot
         w(i,j)=sipco(i)*f(i,j)+simco(i)*f(n-i+2,j)
120      continue
c
         call fft991(w,f,trigs,ifax,1,n+2,n,lot,-1)
c
         do 130 j=1,lot
         f(1,j)=w(1,j)
         f(n,j)=w(n+1,j)
130      continue
         do 140 ii=2,n-2,2
         do 140 j=1,lot
         f(ii,j)=w(ii+1,j)-w(ii+2,j)
         f(ii+1,j)=w(ii+1,j)+w(ii+2,j)
140      continue
      endif
      return
      end
c
c ***
c
c
c
c
      function abssuv(array,i0,i1,j0,j1,k0,k1)
c-----------------------------------------------------------------------
c abssum=sum of absolute value of array
c
c to speed up the process the 3-d array has been vectorized
c------------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      dimension array(*)
      lmn=(i1-i0+1)*(j1-j0+1)*(k1-k0+1)
      abssuv=0.0
      do 10 ijk=1,lmn
      abssuv=abssuv+abs(array(ijk))
10    continue
      return
      end
c
c
c
      function qsatf(pt,p)
c-----------------------------------------------------------------------
c compute saturation specific humidity using tetens' formula
c given potential temperature and pressure
c-----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      parameter(r=287.,cp=1005.,akapa=r/cp,p00=1.0e5)
      t= pt*(p/p00)**akapa
      esat=6.11*exp( 17.27*(t-273.)/(t-36.))
      qsatf=0.622*esat/(0.01*p- 0.378*esat )
      return
      end




      subroutine ptbc(nx,ny,ns,pt,ptbx,ptby,ptcc,iobptx,iobpty
     :   ,ptbout,prt,nchn1)
c
c updates boundary variables from data obtained in the previous time
c step (radiative boundary conditions).
c
      implicit real*8(a-h,o-z)
      logical ptbout,prt
      dimension pt(0:nx+1,0:ny+1,0:ns+1)
      dimension ptbx(2,0:ny+1,0:ns+1),ptby(0:nx+1,2,0:ns+1)
     :   ,ptcc(2,2,0:ns+1)
c
      if(ptbout .and. prt) then
         call wrigar(ptbx,1,2,0,ny+1,0,ns+1,1,2,1,ny+1,1,ns-1
     :      ,'ptbx    ',0.d0,3)
         call wrigar(ptby,0,nx+1,1,2,0,ns+1,1,nx+1,1,2,1,ns-1
     :      ,'ptby    ',0.d0,2)
         write(nchn1,*) 'ptcc'
         do 777 is=1,ns-1
         write(nchn1,*) is,ptcc(1,1,is),ptcc(2,1,is),ptcc(2,2,is)
     :      ,ptcc(1,2,is)
777      continue
      endif
c
      if(iobptx.eq.1) then
         do 21 iy=1,ny+1
         do 21 is=1,ns-1
         pt(1,iy,is)=ptbx(1,iy,is)
         pt(nx+1,iy,is)=ptbx(2,iy,is)
21       continue
      elseif (iobptx.eq.0) then
         do 22 iy=1,ny+1
         do 22 is=1,ns-1
         pt(1,iy,is)=pt(2,iy,is)
         pt(nx+1,iy,is)=pt(nx,iy,is)
22       continue
      elseif (iobptx.eq.2) then
         do 23 iy=1,ny+1
         do 23 is=1,ns-1
         pt(1,iy,is)=pt(nx,iy,is)
         pt(nx+1,iy,is)=pt(2,iy,is)
23       continue
      
      endif
c
      if(iobpty.eq.1) then
         do 11 is=1,ns-1
         do 11 ix=1,nx+1
         pt(ix,1,is)=ptby(ix,1,is)
         pt(ix,ny+1,is)=ptby(ix,2,is)
11       continue
      elseif(iobpty.eq.0) then
         do 12 is=1,ns-1
         do 12 ix=1,nx+1
         pt(ix,1,is)=pt(ix,2,is)
         pt(ix,ny+1,is)=pt(ix,ny,is)
12       continue
      elseif(iobpty.eq.2) then
         do 13 is=1,ns-1
         do 13 ix=1,nx+1
         pt(ix,1,is)=pt(ix,ny,is)
         pt(ix,ny+1,is)=pt(ix,2,is)
13       continue
      endif
c
      if(iobptx.eq.1 .and. iobpty.eq.1) then
         do 100 is=1,ns-1
         pt(1,1,is)=ptcc(1,1,is)
         pt(nx+1,1,is)=ptcc(2,1,is)
         pt(nx+1,ny+1,is)=ptcc(2,2,is)
         pt(1,ny+1,is)=ptcc(1,2,is)
100      continue
      endif
c
      return
      end
c
c
c
      subroutine uvbc(nx,ny,ns
     :   ,u,ubx,ubx2,uby,uby2,ucc,ucc2,iobuux,iobuuy
     :   ,v,vbx,vbx2,vby,vby2,vcc,vcc2,iobvvx,iobvvy)
c
c updates boundary variables from data obtained in the previous time
c step (radiative boundary conditions).
c
      implicit real*8(a-h,o-z)
      dimension u(0:nx+1,0:ny+1,0:ns+1),v(0:nx+1,0:ny+1,0:ns+1)
      dimension ubx(2,0:ny+1,0:ns+1),ubx2(2,0:ny+1,0:ns+1)
      dimension vbx(2,0:ny+1,0:ns+1),vbx2(2,0:ny+1,0:ns+1)
      dimension uby(0:nx+1,2,0:ns+1),uby2(0:nx+1,2,0:ns+1)
      dimension vby(0:nx+1,2,0:ns+1),vby2(0:nx+1,2,0:ns+1)
      dimension ucc(2,2,0:ns+1),ucc2(2,2,0:ns+1)
      dimension vcc(2,2,0:ns+1),vcc2(2,2,0:ns+1)
c
      if(iobuux.eq.1) then
         do 72 iy=2,ny
         do 72 is=1,ns-1
         u(1,iy,is)=ubx(1,iy,is)
         u(nx,iy,is)=ubx(2,iy,is)
72       continue
c
         do 75 iy=1,ny+1
         do 75 is=1,ns-1
         u(0,iy,is)=ubx2(1,iy,is)
         u(nx+1,iy,is)=ubx2(2,iy,is)
75       continue
      elseif (iobuux.eq.0) then
         do 90 iy=0,ny+1
         do 90 is=1,ns-1
         u(1,iy,is)=u(2,iy,is)
         u(nx,iy,is)=u(nx-1,iy,is)
         u(0,iy,is)=u(1,iy,is)
         u(nx+1,iy,is)=u(nx,iy,is)
90       continue
      elseif(iobuux.eq.2) then
         do 91 iy=0,ny+1
         do 91 is=1,ns-1
         u(1,iy,is)=u(nx-1,iy,is)
         u(nx,iy,is)=u(2,iy,is)
         u(0,iy,is)=u(nx-2,iy,is)
         u(nx+1,iy,is)=u(3,iy,is)
91       continue
      endif
c
      if(iobvvx.eq.1) then
         do 74 iy=2,ny-1
         do 74 is=1,ns-1
         v(1,iy,is)=vbx(1,iy,is)
         v(nx+1,iy,is)=vbx(2,iy,is)
74       continue
c
         do 76 iy=1,ny
         do 76 is=1,ns-1
         v(0,iy,is)=vbx2(1,iy,is)
76       continue
      elseif (iobvvx.eq.0) then
         do 94 iy=0,ny+1
         do 94 is=1,ns-1
         v(1,iy,is)=v(2,iy,is)
         v(nx+1,iy,is)=v(nx,iy,is)
         v(0,iy,is)=v(1,iy,is)
94       continue
      elseif (iobvvx.eq.2) then
         do 95 iy=0,ny+1
         do 95 is=1,ns-1
         v(1,iy,is)=v(nx-1,iy,is)
         v(nx+1,iy,is)=v(3,iy,is)
         v(0,iy,is)=v(nx-2,iy,is)
95       continue
      endif
c
      if(iobuuy.eq.1) then
         do 60 is=1,ns-1
         do 61 ix=2,nx-1
         u(ix,1,is)=uby(ix,1,is)
         u(ix,ny+1,is)=uby(ix,2,is)
61       continue
         do 65 ix=1,nx
         u(ix,0,is)=uby2(ix,1,is)
65       continue
60       continue
      else
         do 80 is=1,ns-1
         do 80 ix=0,nx+1
         u(ix,1,is)=u(ix,2,is)
         u(ix,0,is)=u(ix,1,is)
         u(ix,ny+1,is)=u(ix,ny,is)
80       continue
      endif
c
      if(iobvvy.eq.1) then
         do 160 is=1,ns-1
         do 63 ix=2,nx
         v(ix,1,is)=vby(ix,1,is)
         v(ix,ny,is)=vby(ix,2,is)
63       continue
         do 66 ix=1,nx+1
         v(ix,0,is)=vby2(ix,1,is)
         v(ix,ny+1,is)=vby2(ix,2,is)
66       continue
160      continue
      else
         do 180 is=1,ns-1
         do 180 ix=0,nx+1
         v(ix,1,is)=v(ix,2,is)
         v(ix,ny,is)=v(ix,ny-1,is)
         v(ix,0,is)=v(ix,1,is)
         v(ix,ny+1,is)=v(ix,ny,is)
180      continue
      endif
c
      if(iobuux.eq.1 .and. iobuuy.eq.1) then
         do 100 is=1,ns-1
         u(1,1,is)=ucc(1,1,is)
         u(nx,1,is)=ucc(2,1,is)
         u(nx,ny+1,is)=ucc(2,2,is)
         u(1,ny+1,is)=ucc(1,2,is)
         u(0,0,is)=ucc2(1,1,is)
         u(nx+1,0,is)=ucc2(2,1,is)
100      continue
      endif
c
      if(iobvvx.eq.1 .and. iobvvy.eq.1) then
         do 200 is=1,ns-1
         v(1,1,is)=vcc(1,1,is)
         v(nx+1,1,is)=vcc(2,1,is)
         v(nx+1,ny,is)=vcc(2,2,is)
         v(1,ny,is)=vcc(1,2,is)
         v(0,0,is)=vcc2(1,1,is)
         v(0,ny+1,is)=vcc2(1,2,is)
200      continue
      endif
c
      return
      end
c
c
c
      subroutine inipro(fdat,zdat,ndat,fpro,zpro,npro)
c
c interpolates from grid zdat to grid zpro (linear interpolation).
c
      implicit real*8(a-h,o-z)
      dimension fdat(0:ndat),zdat(0:ndat),zpro(0:npro),fpro(0:npro+1)
c
      j1=0
      do 10 i=1,ndat
      do 20 j=j1,npro

      if(zpro(j).gt.zdat(i)) then
         j1=j
         go to 10
      endif
      fpro(j)=fdat(i-1)+(fdat(i)-fdat(i-1))*(zpro(j)-zdat(i-1))
     :        /(zdat(i)-zdat(i-1))
!      write(*,*) 'inipro:',zpro(j),i,fdat(i),j1
 !     if(j.lt.100) then
  !    write(0,*) fpro(j),zpro(j),j
!	endif
20    continue
      if(j.eq.npro) exit
10    continue
      fpro(npro+1)=fpro(npro)
      return
      end
c
c
c
      subroutine extrah(nx1,ny1,var,i0,i1,j0,j1,k0,k1)
c
c extrapolates values of a variable to its outside horizontal border
c
      implicit real*8(a-h,o-z)
      dimension var(0:nx1,0:ny1,k0:k1)
c
      do 11 k=k0,k1
      do 11 ix=i0+1,i1-1
      var(ix,j0,k)=var(ix,j0+1,k)
      var(ix,j1,k)=var(ix,j1-1,k)
11    continue
c
      do 12 k=k0,k1
      do 12 iy=j0,j1
      var(i0,iy,k)=var(i0+1,iy,k)
      var(i1,iy,k)=var(i1-1,iy,k)
12    continue
c
      return
      end
c
c
c
      subroutine extrax(nx1,ny1,var,i0,i1,j0,j1,k0,k1)
c
c extrapolates values of a variable to its outside horizontal border
c
      implicit real*8(a-h,o-z)
      dimension var(0:nx1,0:ny1,k0:k1)
c
      do 12 k=k0,k1
      do 12 iy=j0,j1
      var(i0,iy,k)=var(i0+1,iy,k)
      var(i1,iy,k)=var(i1-1,iy,k)
12    continue
      return
      end
c
c
c
      subroutine extray(nx1,ny1,var,i0,i1,j0,j1,k0,k1)
c
c extrapolates values of a variable to its outside horizontal border
c
      implicit real*8(a-h,o-z)
      dimension var(0:nx1,0:ny1,k0:k1)
c
      do 11 k=k0,k1
      do 11 ix=i0,i1
      var(ix,j0,k)=var(ix,j0+1,k)
      var(ix,j1,k)=var(ix,j1-1,k)
11    continue
      return
      end
c
c
c
      subroutine extrpp(nx2,ny2,var,i0,i1,j0,j1)
c
c extrapolates values of a variable to its outside horizontal border
c
      implicit real*8(a-h,o-z)
      dimension var(0:nx2,0:ny2)
c
      do 11 ix=i0+1,i1-1
      var(ix,j0)=var(ix,j0+1)
      var(ix,j1)=var(ix,j1-1)
11    continue
c
      do 12 iy=j0,j1
      var(i0,iy)=var(i0+1,iy)
      var(i1,iy)=var(i1-1,iy)
12    continue
c
      return
      end
c
c
c
      subroutine extra3(nx1,ny1,ns1,var,i0,i1,j0,j1,k0,k1)
c
c extrapolates values of a variable to its outside border
c in all 3 directions.
c
      implicit real*8(a-h,o-z)
      dimension var(0:nx1,0:ny1,0:ns1)
c
      do 10 iy=j0+1,j1-1
      do 10 ix=i0+1,i1-1
      var(ix,iy,k0)=var(ix,iy,k0+1)
      var(ix,iy,k1)=var(ix,iy,k1-1)
10    continue
c
      do 20 is=k0,k1
      do 20 ix=i0+1,i1-1
      var(ix,j0,is)=var(ix,j0+1,is)
      var(ix,j1,is)=var(ix,j1-1,is)
20    continue
c
      do 30 is=k0,k1
      do 30 iy=j0,j1
      var(i0,iy,is)=var(i0+1,iy,is)
      var(i1,iy,is)=var(i1-1,iy,is)
30    continue
c
      return
      end
c
c
c
      subroutine radbch(abcc,abcx,abcy,a,am1,am2,nxm,nym,nsm
     :   ,k0,k1,ixl,ixr,iyl,iyr)
c
c-----------------------------------------------------------------------
c to compute radiative boundary conditions
c for horizontal boundaries of the field a.
c
c following raymond and kuo (qj 1984, 535-551)
c modified (p.miranda,1989)
c
c-----------------------------------------------------------------------
c
c
      implicit real*8(a-h,o-z)
      parameter (zero0=1.d-30)
      dimension a(0:nxm,0:nym,0:nsm),am1(0:nxm,0:nym,0:nsm)
     :   ,am2(0:nxm,0:nym,0:nsm)
      dimension abcx(2,0:nym,0:nsm),abcy(0:nxm,2,0:nsm),abcc(2,2,0:nsm)
      common/radpar/dx,dy,dt,dxdt,dydt,rxmax,rymax
c
c      call radrkh(abcc,abcx,abcy,a,am1,am2,nxm,nym,nsm                 #exp#
c     :   ,k0,k1,ixl,ixr,iyl,iyr)                                       #exp#
c      return                                                           #exp#
c      call radbco(abcc,abcx,abcy,a,am1,am2,nxm,nym,nsm                 #exp#
c     :   ,k0,k1,ixl,ixr,iyl,iyr)                                       #exp#
c      return                                                           #exp#
c
c      write(*,*) 'radpar1:',nxm,nym,nsm
c     :   ,k0,k1,ixl,ixr,iyl,iyr
c      write(*,*)'radpar:',dx,dy,dt,dxdt,dydt,rxmax,rymax
c      if(nsm.gt.100) then
c         write(6,*)  'fatal error in radbch'
c         stop
c      endif
c
c x boundaries (ix=ixl,ixr)
c
      ib=ixl
      ib1=ib+1
      ib2=ib+2
c     write(*,*) 'radbch:',dxdt
      do 10 lr=1,2
      do 20 iy=iyl+1,iyr-1
c
      do 20 is=k0,k1
c      dadx=(0.5*(a(ib1,iy,is)+am2(ib1,iy,is))-am1(ib2,iy,is))/dx
      dadx=(am1(ib,iy,is)-am1(ib2,iy,is))/(2.*dx)
      dady=(am1(ib1,iy+1,is)-am1(ib1,iy-1,is))/(2.*dy)
      gggg=dadx*dadx+dady*dady
      cc=(am2(ib1,iy,is)-a(ib1,iy,is))/(2.*dt*(gggg+zero0))
c     write(*,*) 'radbch:',cc,dadx,lr,iy,is


      rx=cc*dadx/dxdt
c     rx=dim(rx,0.)-dim(rx,rxmax)
      rx=max(0.d0,min(rx,rxmax))
      ry=cc*dady/dydt
c     ry=dim(ry,-rymax)-dim(ry,rymax)
      ry=max(-rymax,min(ry,rymax))
      abcx(lr,iy,is)=((1.-rx)*am1(ib,iy,is)+2.*rx*a(ib1,iy,is)
     :   -ry*(a(ib,iy+1,is)-a(ib,iy-1,is))
     :   )/(1.+rx)
c      rxx(lr,iy,is)=rx
c      ryx(lr,iy,is)=ry
20    continue
      ib=ixr
      ib1=ib-1
      ib2=ib-2
10    continue
c      write(6,*) 'r:'
c      write(6,1000)  (rxx(1,jj,15),jj=iyl+1,iyr-1)
c      write(6,1000)  (rxx(2,jj,15),jj=iyl+1,iyr-1)
c      write(6,1000)  (ryx(1,jj,15),jj=iyl+1,iyr-1)
c      write(6,1000)  (ryx(2,jj,15),jj=iyl+1,iyr-1)
1000  format(12f6.2)
c
c y boundaries (iy=iyl,iyr)
c
      ib=iyl
      ib1=ib+1
      ib2=ib+2
      do 110 lr=1,2
      do 120 ix=ixl+1,ixr-1
c
      do 120 is=k0,k1
      dadx=(am1(ix+1,ib1,is)-am1(ix-1,ib1,is))/(2.*dx)
      dady=(am1(ix,ib,is)-am1(ix,ib2,is))/(2.*dy)
c      dady=(0.5*(a(ix,ib1,is)+am2(ix,ib1,is))-am1(ix,ib2,is))/dy
      gggg=dadx*dadx+dady*dady
      cc=(am2(ix,ib1,is)-a(ix,ib1,is))/(2.*dt*(gggg+zero0))
      rx=cc*dadx/dxdt
c     rx=dim(rx,-rxmax)-dim(rx,rxmax)
      rx=max(-rxmax,min(rx,rxmax))
      ry=cc*dady/dydt
c     ry=dim(ry,0.)-dim(ry,rymax)
      ry=max(0.d0,min(ry,rymax))
      abcy(ix,lr,is)=((1.-ry)*am1(ix,ib,is)+2.*ry*a(ix,ib1,is)
     :   -rx*(a(ix+1,ib,is)-a(ix-1,ib,is))
     :   )/(1.+ry)
120   continue
      ib=iyr
      ib1=ib-1
      ib2=ib-2
110   continue
c
c corner points:
c
      iby=iyl
      iby1=iby+1
      iby2=iby+2
      do 310 ly=1,2
      ibx=ixl
      ibx1=ibx+1
      ibx2=ibx+2
      do 320 lx=1,2
c
      do 340 is=k0,k1
c      dadx=(0.5*(a(ibx1,iby1,is)+am2(ibx1,iby1,is))-am1(ibx2,iby1,is))
c     :   /dx
      dadx=(am1(ibx,iby1,is)-am1(ibx2,iby1,is))/(2.*dx)
c      dady=(0.5*(a(ibx1,iby1,is)+am2(ibx1,iby1,is))-am1(ibx1,iby2,is))
c     :   /dy
      dady=(am1(ibx1,iby,is)-am1(ibx1,iby2,is))/(2.*dy)
      gggg=dadx*dadx+dady*dady
      cc=(am2(ibx1,iby1,is)-a(ibx1,iby1,is))/(2.*dt*(gggg+zero0))
      rx=cc*dadx/dxdt
c     rx=dim(rx,0.)-dim(rx,rxmax)
      rx=max(0.d0,min(rx,rxmax))
      ry=cc*dady/dydt
c     ry=dim(ry,0.)-dim(ry,rymax)
      ry=max(0.d0,min(ry,rymax))
      abcc(lx,ly,is)=((1.-rx-ry)*am1(ibx,iby,is)
     :   +2.*rx*a(ibx1,iby,is)+2.*ry*a(ibx,iby1,is)
     :   )/(1.+rx+ry)
340   continue
      ibx=ixr
      ibx1=ibx-1
      ibx2=ibx-2
320   continue
      iby=iyr
      iby1=iby-1
      iby2=iby-2
310   continue
c
      do 400 is=k0,k1
      abcx(1,iyl,is)=abcc(1,1,is)
      abcx(1,iyr,is)=abcc(1,2,is)
      abcx(2,iyl,is)=abcc(2,1,is)
      abcx(2,iyr,is)=abcc(2,2,is)
      abcy(ixl,1,is)=abcx(1,iyl,is)
      abcy(ixr,1,is)=abcx(2,iyl,is)
      abcy(ixl,2,is)=abcx(1,iyr,is)
      abcy(ixr,2,is)=abcx(2,iyr,is)
400   continue
c
      return
      end




      subroutine fctlim(ftd,fa,flx,fly,fls
     :   ,sm,cx,cy,cs,rplus,rminus,pplus,pminus
     :   ,nx,ny,ns,dxys)
c
c calculates flux limit following zalesak, using strong limiter of
c boris and book.
c
c uses work spaces 4-9
c
c note the equivalences !
c
c
      implicit real*8(a-h,o-z)
      parameter (zero0=1.d-30)
      dimension ftd(0:nx+1,0:ny+1,0:ns+1),fa(0:nx+1,0:ny+1,0:ns+1)
      dimension flx(0:nx+1,0:ny+1,0:ns+1),fly(0:nx+1,0:ny+1,0:ns+1)
     :   ,fls(0:nx+1,0:ny+1,0:ns+1)
      dimension sm(0:nx+1,0:ny+1,-1:ns)
      dimension cx(0:nx+1,0:ny+1,0:ns+1),cy(0:nx+1,0:ny+1,0:ns+1)
     :   ,cs(0:nx+1,0:ny+1,0:ns+1)
      dimension pplus(0:nx+1,0:ny+1,0:ns+1),rplus(0:nx+1,0:ny+1,0:ns+1)
c     :   ,qplus(0:nx+1,0:ny+1,0:ns+1)
      dimension pminus(0:nx+1,0:ny+1,0:ns+1)
     :   ,rminus(0:nx+1,0:ny+1,0:ns+1)
c     :   ,qminus(0:nx+1,0:ny+1,0:ns+1)
      dimension dxys(0:ns+1)
c
c      equivalence (wk4,sm),(wk5,cx),(wk6,cy),(wk7,cs)
c      equivalence (wk8,rplus,rminus,pplus,pminus),(wk9,qplus,qminus)
c
      save prlim,iyprlim,fold
      logical prlim,iyprlim,fold
      data iyprlim,prlim,fold /.true.,.false., .true.  /
c
c in-line functions:
c
c     qplus(ix,iy,is)=(max(ftd(ix,iy,is),ftd(ix-1,iy,is),ftd(ix+1,iy,is)
c    :   ,ftd(ix,iy-1,is),ftd(ix,iy+1,is),ftd(ix,iy,is-1)
c    :   ,ftd(ix,iy,is+1)
c    :   ,fa(ix,iy,is),fa(ix-1,iy,is),fa(ix+1,iy,is)
c    :   ,fa(ix,iy-1,is),fa(ix,iy+1,is),fa(ix,iy,is-1),fa(ix,iy,is+1))
c    :   -ftd(ix,iy,is))*dxys(is)
c
c     qminus(ix,iy,is)=(ftd(ix,iy,is)-min(ftd(ix,iy,is)
c    :   ,ftd(ix-1,iy,is),ftd(ix+1,iy,is)
c    :   ,ftd(ix,iy-1,is),ftd(ix,iy+1,is)
c    :   ,ftd(ix,iy,is-1),ftd(ix,iy,is+1)
c    :   ,fa(ix,iy,is),fa(ix-1,iy,is),fa(ix+1,iy,is)
c    :   ,fa(ix,iy-1,is),fa(ix,iy+1,is)
c    :   ,fa(ix,iy,is-1),fa(ix,iy,is+1)))*dxys(is)
c
      if( iyprlim) then
c
c apply strong pre-limiter on flx
c
         do 10 is=1,ns-1
         do 10 iy=2,ny
         do 10 ix=1,nx
         sm(ix,iy,is)=(ftd(ix+1,iy,is)-ftd(ix,iy,is))
10       continue
c
c neuman b.c. in the horizontal
c
         do 11 ix=0,nx+1,nx+1
         do 11 is=1,ns-1
         do 11 iy=2,ny
         sm(ix,iy,is)=0.
11       continue
c
         do 20 is=1,ns-1
         do 20 iy=2,ny
         do 20 ix=1,nx
         a=flx(ix,iy,is)*sm(ix,iy,is)
         b=flx(ix,iy,is)*sm(ix+1,iy,is)
         c=flx(ix,iy,is)*sm(ix-1,iy,is)
c         if(a.lt.0. .and. (b.lt.0. .or. c.lt.0.)) flx(ix,iy,is)=0.
         flx(ix,iy,is)=flx(ix,iy,is)+flx(ix,iy,is)*(sign(0.5d0,a)-0.5)
     :      *min(1.d0,-sign(0.5d0,b)-sign(0.5d0,c)+1.d0)
20       continue
c
c use simplified flux pre-limiter for boundary fluxes if b.c. is not
c of the three basic types.
c
c         do 30 iy=1,ny+1
c         sm(1,iy,is)=flx(1,iy,is)*(ftd(2,iy,is)-ftd(1,iy,is))
c         flx(1,iy,is)=flx(1,iy,is)*max(0.,sign(1.,sm(1,iy,is)))
c30       continue
c         do 31 iy=1,ny+1
c         sm(ix,iy,is)=flx(ix,iy,is)*(ftd(ix+1,iy,is)-ftd(ix,iy,is))
c         flx(ix,iy,is)=flx(ix,iy,is)*max(0.,sign(1.,sm(ix,iy,is)))
c31       continue
c
c apply strong pre-limiter on fly
c
         do 110 is=1,ns-1
         do 110 iy=1,ny
         do 110 ix=2,nx
         sm(ix,iy,is)=(ftd(ix,iy+1,is)-ftd(ix,iy,is))
110      continue
c
c neuman b.c. in the horizontal
c
         do 111 iy=0,ny+1,ny+1
         do 111 is=1,ns-1
         do 111 ix=2,nx
         sm(ix,iy,is)=0.
111      continue
c
         do 120 is=1,ns-1
         do 120 iy=1,ny
         do 120 ix=2,nx
         a=fly(ix,iy,is)*sm(ix,iy,is)
         b=fly(ix,iy,is)*sm(ix,iy+1,is)
         c=fly(ix,iy,is)*sm(ix,iy-1,is)
c         if(a.lt.0. .and. (b.lt.0. .or. c.lt.0.)) fly(ix,iy,is)=0.
         fly(ix,iy,is)=fly(ix,iy,is)+fly(ix,iy,is)*(sign(0.5d0,a)-0.5)
     :      *min(1.d0,-sign(0.5d0,b)-sign(0.5d0,c)+1.d0)
120      continue
c
c apply strong pre-limiter on fls
c
         do 210 is=0,ns-1
         do 210 iy=2,ny
         do 210 ix=2,nx
         sm(ix,iy,is)=(ftd(ix,iy,is+1)-ftd(ix,iy,is))
210      continue
c
c vertical b.c.:
c
         do 211 is=-1,ns,ns+1
         do 211 iy=2,ny
         do 211 ix=2,nx
         sm(ix,iy,is)=0.
211      continue
c
         do 220 is=0,ns-1
         do 220 iy=2,ny
         do 220 ix=2,nx
         a=fls(ix,iy,is)*sm(ix,iy,is)
         b=fls(ix,iy,is)*sm(ix,iy,is+1)
         c=fls(ix,iy,is)*sm(ix,iy,is-1)
c         if(a.lt.0. .and. (b.lt.0. .or. c.lt.0.)) fls(ix,iy,is)=0.
         fls(ix,iy,is)=fls(ix,iy,is)+fls(ix,iy,is)*(sign(0.5d0,a)-0.5)
     :      *min(1.d0,-sign(0.5d0,b)-sign(0.5d0,c)+1.d0)
220      continue
c
      elseif(prlim) then
c
c use simplified pre-limiter on fluxes
c
         do 2 is=1,ns-1
         do 2 iy=2,ny
         do 2 ix=1,nx
         sm(ix,iy,is)=flx(ix,iy,is)*(ftd(ix+1,iy,is)-ftd(ix,iy,is))
         flx(ix,iy,is)=flx(ix,iy,is)*max(0.d0,sign(1.d0,sm(ix,iy,is)))
2        continue
c
         do 4 is=1,ns-1
         do 4 iy=1,ny
         do 4 ix=2,nx
         sm(ix,iy,is)=fly(ix,iy,is)*(ftd(ix,iy+1,is)-ftd(ix,iy,is))
         fly(ix,iy,is)=fly(ix,iy,is)*max(0.d0,sign(1.d0,sm(ix,iy,is)))
4        continue
c
         do 6 is=0,ns-1
         do 6 iy=2,ny
         do 6 ix=2,nx
         sm(ix,iy,is)=fls(ix,iy,is)*(ftd(ix,iy,is+1)-ftd(ix,iy,is))
         fls(ix,iy,is)=fls(ix,iy,is)*max(0.d0,sign(1.d0,sm(ix,iy,is)))
6        continue
      endif
c
c find upper bounds for f(ix,iy,is)
c
c      if(.not.fold) then
c         do 300 is=0,ns
c         do 300 iy=1,ny+1
c         do 300 ix=1,nx+1
c         sm(ix,iy,is)=ftd(ix,iy,is)
c300      continue
c      else
c         do 310 is=0,ns
c         do 310 iy=1,ny+1
c         do 310 ix=1,nx+1
c         sm(ix,iy,is)=max(ftd(ix,iy,is),fa(ix,iy,is))
c310      continue
c      endif
c
c calculate quantity q(ix,iy,is)+
c
c      do 400 is=1,ns-1
c      do 400 iy=2,ny
c      do 400 ix=2,nx
c      qplus(ix,iy,is)=(max(sm(ix,iy,is),sm(ix-1,iy,is),sm(ix+1,iy,is)
c     :   ,sm(ix,iy-1,is),sm(ix,iy+1,is),sm(ix,iy,is-1),sm(ix,iy,is+1))
c     :   -ftd(ix,iy,is))*dxys(is)
c400   continue
c
c calculate total flux directing towards point (ix,iy,is) p(ix,iy,is)+ .
c
      call extra3(nx+1,ny+1,ns+1,flx,0,nx+1,1,ny+1,0,ns)
      call extra3(nx+1,ny+1,ns+1,fly,1,nx+1,0,ny+1,0,ns)
      call extrah(nx+1,ny+1,fls,1,nx+1,1,ny+1,0,ns-1)
c
      do 800 iy=1,ny+1
      do 800 ix=1,nx+1
      fls(ix,iy,ns)=0.
800   continue
c
      do 410 is=1,ns
      do 410 iy=1,ny+1
      do 410 ix=1,nx+1
      pplus(ix,iy,is)=max(0.d0,flx(ix-1,iy,is))-min(0.d0,flx(ix,iy,is))
     :   +max(0.d0,fly(ix,iy-1,is))-min(0.d0,fly(ix,iy,is))
     :   +max(0.d0,fls(ix,iy,is-1))-min(0.d0,fls(ix,iy,is))
410   continue
c
c calculate inward flux amplification factor r(ix,iy,is)+
c r+=min(1.,q+/p+)
c
      do 420 is=1,ns
      do 420 iy=1,ny+1
      do 420 ix=1,nx+1
      qplus=(max(ftd(ix,iy,is),ftd(ix-1,iy,is),ftd(ix+1,iy,is)
     :   ,ftd(ix,iy-1,is),ftd(ix,iy+1,is),ftd(ix,iy,is-1)
     :   ,ftd(ix,iy,is+1)
     :   ,fa(ix,iy,is),fa(ix-1,iy,is),fa(ix+1,iy,is)
     :   ,fa(ix,iy-1,is),fa(ix,iy+1,is),fa(ix,iy,is-1),fa(ix,iy,is+1))
     :   -ftd(ix,iy,is))*dxys(is)
      rplus(ix,iy,is)=min(1.d0,qplus/(pplus(ix,iy,is)+zero0))
420   continue
c
c determine boundary values of r+ .
c
      call extra3(nx+1,ny+1,ns+1,rplus,1,nx+1,1,ny+1,0,ns)
c
c calculate correcting factor c(ix,iy,is) from inward flux amplification
c factor r+ .
c
      do 440 is=1,ns-1
      do 440 iy=2,ny
      do 440 ix=1,nx
      cx(ix,iy,is)=max(-sign(rplus(ix,iy,is),flx(ix,iy,is))
     :   ,sign(rplus(ix+1,iy,is),flx(ix,iy,is)))
440   continue
c
      do 450 is=1,ns-1
      do 450 iy=1,ny
      do 450 ix=2,nx
      cy(ix,iy,is)=max(-sign(rplus(ix,iy,is),fly(ix,iy,is))
     :   ,sign(rplus(ix,iy+1,is),fly(ix,iy,is)))
450   continue
c
      do 460 is=0,ns-1
      do 460 iy=2,ny
      do 460 ix=2,nx
      cs(ix,iy,is)=max(-sign(rplus(ix,iy,is),fls(ix,iy,is))
     :   ,sign(rplus(ix,iy,is+1),fls(ix,iy,is)))
460   continue
c
c find lower bounds for advected fields.
c
c      if(.not.fold) then
c         do 500 is=0,ns
c         do 500 iy=1,ny+1
c         do 500 ix=1,nx+1
c         sm(ix,iy,is)=ftd(ix,iy,is)
c500      continue
c      else
c         do 510 is=0,ns
c         do 510 iy=1,ny+1
c         do 510 ix=1,nx+1
c         sm(ix,iy,is)=min(ftd(ix,iy,is),fa(ix,iy,is))
c510      continue
c      endif
c
c calculate term q(ix,iy,is)-
c

c      do 520 is=1,ns-1
c      do 520 iy=2,ny
c      do 520 ix=2,nx
c      qminus(ix,iy,is)=(ftd(ix,iy,is)-min(sm(ix,iy,is)
c     :   ,sm(ix-1,iy,is),sm(ix+1,iy,is)
c     :   ,sm(ix,iy-1,is),sm(ix,iy+1,is)
c     :   ,sm(ix,iy,is-1),sm(ix,iy,is+1)))*dxys(is)
c520   continue
c
c apply flux correction on flx.
c
      do 530 is=1,ns-1
      do 530 iy=2,ny
      do 530 ix=1,nx
      sm(ix,iy,is)=cx(ix,iy,is)*flx(ix,iy,is)
530   continue
c
c calculate flux directing away  from grid point (ix,iy,is) p(ix,iy,is)- .
c
      do 540 is=1,ns-1
      do 540 iy=2,ny
      do 540 ix=2,nx
      pminus(ix,iy,is)=max(0.d0,sm(ix,iy,is))-min(0.d0,sm(ix-1,iy,is))
540   continue
c
c apply flux correction on fly.
c
      do 550 is=1,ns-1
      do 550 iy=1,ny
      do 550 ix=2,nx
      sm(ix,iy,is)=cy(ix,iy,is)*fly(ix,iy,is)
550   continue
c
c calculate total flux directing away  from grid point (ix,iy,is) p(ix,iy,is)- .
c
      do 560 is=1,ns-1
      do 560 iy=2,ny
      do 560 ix=2,nx
      pminus(ix,iy,is)=pminus(ix,iy,is)+max(0.d0,sm(ix,iy,is))
     :   -min(0.d0,sm(ix,iy-1,is))
560   continue
c
c apply flux correction on fls.
c
      do 570 is=0,ns-1
      do 570 iy=2,ny
      do 570 ix=2,nx
      sm(ix,iy,is)=cs(ix,iy,is)*fls(ix,iy,is)
570   continue
c
c calculate total flux directing away from grid point (ix,iy,is).
c
      do 580 is=1,ns-1
      do 580 iy=2,ny
      do 580 ix=2,nx
      pminus(ix,iy,is)=pminus(ix,iy,is)+max(0.d0,sm(ix,iy,is))
     :   -min(0.d0,sm(ix,iy,is-1))
580   continue
c
c calculate outgoing flux amplification factor r(ix,iy,is)- .
c
      do 600 is=1,ns-1
      do 600 iy=2,ny
      do 600 ix=2,nx
      qminus=(ftd(ix,iy,is)-min(ftd(ix,iy,is)
     :   ,ftd(ix-1,iy,is),ftd(ix+1,iy,is)
     :   ,ftd(ix,iy-1,is),ftd(ix,iy+1,is)
     :   ,ftd(ix,iy,is-1),ftd(ix,iy,is+1)
     :   ,fa(ix,iy,is),fa(ix-1,iy,is),fa(ix+1,iy,is)
     :   ,fa(ix,iy-1,is),fa(ix,iy+1,is)
     :   ,fa(ix,iy,is-1),fa(ix,iy,is+1)))*dxys(is)
      rminus(ix,iy,is)=min(1.d0,qminus
     :   /(pminus(ix,iy,is)+zero0))
600   continue
c
c determine boundary values of r- .
c
      call extra3(nx+1,ny+1,ns+1,rminus,1,nx+1,1,ny+1,0,ns)
c
c calculate the final correcting factor c(ix,iy,is).
c
      do 620 is=1,ns-1
      do 620 iy=2,ny
      do 620 ix=1,nx
      cx(ix,iy,is)=cx(ix,iy,is)*
     :   max(sign(rminus(ix,iy,is),flx(ix,iy,is)),
     :   -sign(rminus(ix+1,iy,is),flx(ix,iy,is)))
620   continue
c
      do 630 is=1,ns-1
      do 630 iy=1,ny
      do 630 ix=2,nx
      cy(ix,iy,is)=cy(ix,iy,is)*
     :   max(sign(rminus(ix,iy,is),fly(ix,iy,is)),
     :   -sign(rminus(ix,iy+1,is),fly(ix,iy,is)))
630   continue
c
      do 640 is=0,ns-1
      do 640 iy=2,ny
      do 640 ix=2,nx
      cs(ix,iy,is)=cs(ix,iy,is)*
     :   max(sign(rminus(ix,iy,is),fls(ix,iy,is)),
     :   -sign(rminus(ix,iy,is+1),fls(ix,iy,is)))
640   continue
c
c final flux corrections
c
      do 700 is=1,ns-1
      do 700 iy=2,ny
      do 700 ix=1,nx
      flx(ix,iy,is)=flx(ix,iy,is)*cx(ix,iy,is)
700   continue
c
      do 710 is=1,ns-1
      do 710 iy=1,ny
      do 710 ix=2,nx
      fly(ix,iy,is)=fly(ix,iy,is)*cy(ix,iy,is)
710   continue
c
      do 720 is=0,ns-1
      do 720 iy=2,ny
      do 720 ix=2,nx
      fls(ix,iy,is)=fls(ix,iy,is)*cs(ix,iy,is)
720   continue
c
      end




c     subroutine 'set99' - computes factors of n & trigonometric
c     functins required by fft99 & fft991
c
      subroutine set99(trigs,ifax,n)
c
      implicit real*8(a-h,o-z)
      common /qqqqfft/ ixxx
      dimension trigs(n),ifax(*),jfax(10),lfax(7)
      data lfax/6,8,5,4,3,2,1/
      ixxx=1
c
      del=4.0*asin(1.0)/float(n)
      nil=0
      nhl=(n/2)-1
      do 10 k=nil,nhl
      angle=float(k)*del
      trigs(2*k+1)=cos(angle)
      trigs(2*k+2)=sin(angle)
   10 continue
c
c     find factors of n (8,6,5,4,3,2; only one 8 allowed)
c     look for sixes first, store factors in descending order
      nu=n
      ifac=6
      k=0
      l=1
   20 continue
      if (mod(nu,ifac).ne.0) go to 30
      k=k+1
      jfax(k)=ifac
      if (ifac.ne.8) go to 25
      if (k.eq.1) go to 25
      jfax(1)=8
      jfax(k)=6
   25 continue
      nu=nu/ifac
      if (nu.eq.1) go to 50
      if (ifac.ne.8) go to 20
   30 continue
      l=l+1
      ifac=lfax(l)
      if (ifac.gt.1) go to 20
c
      write(6,'(''nh3dFatalError=Illegal number for fft'',i5)') n
      stop
      return
c
c     now reverse order of factors
   50 continue
      nfax=k
      ifax(1)=nfax
      do 60 i=1,nfax
      ifax(nfax+2-i)=jfax(i)
   60 continue
      return
      end
c
c
c
      subroutine fft991(a,work,trigs,ifax,inc,jump,n,lot,isign)
c
      implicit real*8(a-h,o-z)
      common /qqqqfft/ ixxx
      dimension a(*),work(*),trigs(*),ifax(*)
c
      if(ixxx.ne.1) call set99(trigs,ifax,n)
      nfax=ifax(1)
      nx=n+1
      if (mod(n,2).eq.1) nx=n
      nblox=1+(lot-1)/64
      nvex=lot-(nblox-1)*64
c      write(*,'(''(fft991:'',11i4)') nfax,(ifax(k19),k19=1,10)
      if (isign.eq.-1) go to 300
c
c     isign=+1, spectral to gridpoint transform
c     -----------------------------------------
  100 continue
      istart=1
      do 220 nb=1,nblox
      ia=istart
      i=istart
      do 110 j=1,nvex
      a(i+inc)=0.5*a(i)
      i=i+jump
  110 continue
      if (mod(n,2).eq.1) go to 130
      i=istart+n*inc
      do 120 j=1,nvex
      a(i)=0.5*a(i)
      i=i+jump
  120 continue
  130 continue
      ia=istart+inc
      la=1
      igo=+1
c
      do 160 k=1,nfax
      ifac=ifax(k+1)
      ierr=-1
      if (igo.eq.-1) go to 140
      call rpassm(a(ia),a(ia+la*inc),work(1),work(ifac*la+1),trigs,
     *    inc,1,jump,nx,nvex,n,ifac,la,ierr)
      go to 150
  140 continue
      call rpassm(work(1),work(la+1),a(ia),a(ia+ifac*la*inc),trigs,
     *    1,inc,nx,jump,nvex,n,ifac,la,ierr)
  150 continue
      if (ierr.ne.0) go to 500
      la=ifac*la
      igo=-igo
      ia=istart
  160 continue
c
c     if necessary, copy results back to a
c     ------------------------------------
      if (mod(nfax,2).eq.0) go to 190
      ibase=1
      jbase=ia
      do 180 jj=1,nvex
      i=ibase
      j=jbase
      do 170 ii=1,n
      a(j)=work(i)
      i=i+1
      j=j+inc
  170 continue
      ibase=ibase+nx
      jbase=jbase+jump
  180 continue
  190 continue
c
c     fill in zeros at end_
c     --------------------
      ix=istart+n*inc
      do 210 j=1,nvex
      a(ix)=0.0
      a(ix+inc)=0.0
      ix=ix+jump
  210 continue
c
      istart=istart+nvex*jump
      nvex=64
  220 continue
      return
c
c     isign=-1, gridpoint to spectral transform
c     -----------------------------------------
  300 continue
      istart=1
      do 410 nb=1,nblox
      ia=istart
      la=n
      igo=+1
c
      do 340 k=1,nfax
      ifac=ifax(nfax+2-k)
      la=la/ifac
      ierr=-1
      if (igo.eq.-1) go to 320
      call qpassm(a(ia),a(ia+ifac*la*inc),work(1),work(la+1),trigs,
     *    inc,1,jump,nx,nvex,n,ifac,la,ierr)
      go to 330
  320 continue
      call qpassm(work(1),work(ifac*la+1),a(ia),a(ia+la*inc),trigs,
     *    1,inc,nx,jump,nvex,n,ifac,la,ierr)
  330 continue
      if (ierr.ne.0) go to 500
      igo=-igo
      ia=istart+inc
  340 continue
c
c     if necessary, copy results back to a
c     ------------------------------------
      if (mod(nfax,2).eq.0) go to 370
      ibase=1
      jbase=ia
      do 360 jj=1,nvex
      i=ibase
      j=jbase
      do 350 ii=1,n
      a(j)=work(i)
      i=i+1
      j=j+inc
  350 continue
      ibase=ibase+nx
      jbase=jbase+jump
  360 continue
  370 continue
c
c     shift a(0) & fill in zero imag parts
c     ------------------------------------
      ix=istart
      do 380 j=1,nvex
      a(ix)=a(ix+inc)
      a(ix+inc)=0.0
      ix=ix+jump
  380 continue
      if (mod(n,2).eq.1) go to 400
      iz=istart+(n+1)*inc
      do 390 j=1,nvex
      a(iz)=0.0
      iz=iz+jump
  390 continue
  400 continue
c
      istart=istart+nvex*jump
      nvex=64
  410 continue
      return
c
c     error messages
c     --------------
  500 continue
      go to (510,530,550) ierr
  510 continue
      write(6,520) nvex
  520 format(16h1vector length =,i4,17h, greater than 64)
      go to 570
  530 continue
      write(6,540) ifac
  540 format( 9h1factor =,i3,17h, not catered for)
      go to 570
  550 continue
      write(6,560) ifac
  560 format(9h1factor =,i3,31h, only catered for if la*ifac=n)
  570 continue
      return
      end
c
c
c
      subroutine qpassm(a,b,c,d,trigs,inc1,inc2,inc3,inc4,lot,n,ifac,
     *    la,ierr)
c
      implicit real*8(a-h,o-z)
      dimension a(*),b(*),c(*),d(*),trigs(*)
c
      data sin36/0.587785252292473/,sin72/0.951056516295154/,
     *    qrt5/0.559016994374947/,sin60/0.866025403784437/
c
      m=n/ifac
      iink=la*inc1
      jink=la*inc2
      ijump=(ifac-1)*iink
      kstop=(n-ifac)/(2*ifac)
c
      ibad=1
      if (lot.gt.64) go to 910
      ibase=0
      jbase=0
      igo=ifac-1
      if (igo.eq.7) igo=6
      ibad=2
      if (igo.lt.1.or.igo.gt.6) go to 910
      go to (200,300,400,500,600,800),igo
c
c     coding for factor 2
c     -------------------
  200 continue
      ia=1
      ib=ia+iink
      ja=1
      jb=ja+(2*m-la)*inc2
c
      if (la.eq.m) go to 290
c
      do 220 l=1,la
      i=ibase
      j=jbase
cdir$ ivdep
      do 210 ijk=1,lot
      c(ja+j)=a(ia+i)+a(ib+i)
      c(jb+j)=a(ia+i)-a(ib+i)
      i=i+inc3
      j=j+inc4
  210 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  220 continue
      ja=ja+jink
      jink=2*jink
      jb=jb-jink
      ibase=ibase+ijump
      ijump=2*ijump+iink
      if (ja.eq.jb) go to 260
      do 250 k=la,kstop,la
      kb=k+k
      c1=trigs(kb+1)
      s1=trigs(kb+2)
      jbase=0
      do 240 l=1,la
      i=ibase
      j=jbase
cdir$ ivdep
      do 230 ijk=1,lot
      c(ja+j)=a(ia+i)+(c1*a(ib+i)+s1*b(ib+i))
      c(jb+j)=a(ia+i)-(c1*a(ib+i)+s1*b(ib+i))
      d(ja+j)=(c1*b(ib+i)-s1*a(ib+i))+b(ia+i)
      d(jb+j)=(c1*b(ib+i)-s1*a(ib+i))-b(ia+i)
      i=i+inc3
      j=j+inc4
  230 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  240 continue
      ibase=ibase+ijump
      ja=ja+jink
      jb=jb-jink
  250 continue
      if (ja.gt.jb) go to 900
  260 continue
      jbase=0
      do 280 l=1,la
      i=ibase
      j=jbase
cdir$ ivdep
      do 270 ijk=1,lot
      c(ja+j)=a(ia+i)
      d(ja+j)=-a(ib+i)
      i=i+inc3
      j=j+inc4
  270 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  280 continue
      go to 900
c
  290 continue
      z=1.0/float(n)
      do 294 l=1,la
      i=ibase
      j=jbase
cdir$ ivdep
      do 292 ijk=1,lot
      c(ja+j)=z*(a(ia+i)+a(ib+i))
      c(jb+j)=z*(a(ia+i)-a(ib+i))
      i=i+inc3
      j=j+inc4
  292 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  294 continue
      go to 900
c
c     coding for factor 3
c     -------------------
  300 continue
      ia=1
      ib=ia+iink
      ic=ib+iink
      ja=1
      jb=ja+(2*m-la)*inc2
      jc=jb
c
      if (la.eq.m) go to 390
c
      do 320 l=1,la
      i=ibase
      j=jbase
cdir$ ivdep
      do 310 ijk=1,lot
      c(ja+j)=a(ia+i)+(a(ib+i)+a(ic+i))
      c(jb+j)=a(ia+i)-0.5*(a(ib+i)+a(ic+i))
      d(jb+j)=sin60*(a(ic+i)-a(ib+i))
      i=i+inc3
      j=j+inc4
  310 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  320 continue
      ja=ja+jink
      jink=2*jink
      jb=jb+jink
      jc=jc-jink
      ibase=ibase+ijump
      ijump=2*ijump+iink
      if (ja.eq.jc) go to 360
      do 350 k=la,kstop,la
      kb=k+k
      kc=kb+kb
      c1=trigs(kb+1)
      s1=trigs(kb+2)
      c2=trigs(kc+1)
      s2=trigs(kc+2)
      jbase=0
      do 340 l=1,la
      i=ibase
      j=jbase
cdir$ ivdep
      do 330 ijk=1,lot
      a1=(c1*a(ib+i)+s1*b(ib+i))+(c2*a(ic+i)+s2*b(ic+i))
      b1=(c1*b(ib+i)-s1*a(ib+i))+(c2*b(ic+i)-s2*a(ic+i))
      a2=a(ia+i)-0.5*a1
      b2=b(ia+i)-0.5*b1
      a3=sin60*((c1*a(ib+i)+s1*b(ib+i))-(c2*a(ic+i)+s2*b(ic+i)))
      b3=sin60*((c1*b(ib+i)-s1*a(ib+i))-(c2*b(ic+i)-s2*a(ic+i)))
      c(ja+j)=a(ia+i)+a1
      d(ja+j)=b(ia+i)+b1
      c(jb+j)=a2+b3
      d(jb+j)=b2-a3
      c(jc+j)=a2-b3
      d(jc+j)=-(b2+a3)
      i=i+inc3
      j=j+inc4
  330 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  340 continue
      ibase=ibase+ijump
      ja=ja+jink
      jb=jb+jink
      jc=jc-jink
  350 continue
      if (ja.gt.jc) go to 900
  360 continue
      jbase=0
      do 380 l=1,la
      i=ibase
      j=jbase
cdir$ ivdep
      do 370 ijk=1,lot
      c(ja+j)=a(ia+i)+0.5*(a(ib+i)-a(ic+i))
      d(ja+j)=-sin60*(a(ib+i)+a(ic+i))
      c(jb+j)=a(ia+i)-(a(ib+i)-a(ic+i))
      i=i+inc3
      j=j+inc4
  370 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  380 continue
      go to 900
c
  390 continue
      z=1.0/float(n)
      zsin60=z*sin60
      do 394 l=1,la
      i=ibase
      j=jbase
cdir$ ivdep
      do 392 ijk=1,lot
      c(ja+j)=z*(a(ia+i)+(a(ib+i)+a(ic+i)))
      c(jb+j)=z*(a(ia+i)-0.5*(a(ib+i)+a(ic+i)))
      d(jb+j)=zsin60*(a(ic+i)-a(ib+i))
      i=i+inc3
      j=j+inc4
  392 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  394 continue
      go to 900
c
c     coding for factor 4
c     -------------------
  400 continue
      ia=1
      ib=ia+iink
      ic=ib+iink
      id=ic+iink
      ja=1
      jb=ja+(2*m-la)*inc2
      jc=jb+2*m*inc2
      jd=jb
c
      if (la.eq.m) go to 490
c
      do 420 l=1,la
      i=ibase
      j=jbase
cdir$ ivdep
      do 410 ijk=1,lot
      c(ja+j)=(a(ia+i)+a(ic+i))+(a(ib+i)+a(id+i))
      c(jc+j)=(a(ia+i)+a(ic+i))-(a(ib+i)+a(id+i))
      c(jb+j)=a(ia+i)-a(ic+i)
      d(jb+j)=a(id+i)-a(ib+i)
      i=i+inc3
      j=j+inc4
  410 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  420 continue
      ja=ja+jink
      jink=2*jink
      jb=jb+jink
      jc=jc-jink
      jd=jd-jink
      ibase=ibase+ijump
      ijump=2*ijump+iink
      if (jb.eq.jc) go to 460
      do 450 k=la,kstop,la
      kb=k+k
      kc=kb+kb
      kd=kc+kb
      c1=trigs(kb+1)
      s1=trigs(kb+2)
      c2=trigs(kc+1)
      s2=trigs(kc+2)
      c3=trigs(kd+1)
      s3=trigs(kd+2)
      jbase=0
      do 440 l=1,la
      i=ibase
      j=jbase
cdir$ ivdep
      do 430 ijk=1,lot
      a0=a(ia+i)+(c2*a(ic+i)+s2*b(ic+i))
      a2=a(ia+i)-(c2*a(ic+i)+s2*b(ic+i))
      a1=(c1*a(ib+i)+s1*b(ib+i))+(c3*a(id+i)+s3*b(id+i))
      a3=(c1*a(ib+i)+s1*b(ib+i))-(c3*a(id+i)+s3*b(id+i))
      b0=b(ia+i)+(c2*b(ic+i)-s2*a(ic+i))
      b2=b(ia+i)-(c2*b(ic+i)-s2*a(ic+i))
      b1=(c1*b(ib+i)-s1*a(ib+i))+(c3*b(id+i)-s3*a(id+i))
      b3=(c1*b(ib+i)-s1*a(ib+i))-(c3*b(id+i)-s3*a(id+i))
      c(ja+j)=a0+a1
      c(jc+j)=a0-a1
      d(ja+j)=b0+b1
      d(jc+j)=b1-b0
      c(jb+j)=a2+b3
      c(jd+j)=a2-b3
      d(jb+j)=b2-a3
      d(jd+j)=-(b2+a3)
      i=i+inc3
      j=j+inc4
  430 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  440 continue
      ibase=ibase+ijump
      ja=ja+jink
      jb=jb+jink
      jc=jc-jink
      jd=jd-jink
  450 continue
      if (jb.gt.jc) go to 900
  460 continue
      sin45=dsqrt(0.5d0)
      jbase=0
      do 480 l=1,la
      i=ibase
      j=jbase
cdir$ ivdep
      do 470 ijk=1,lot
      c(ja+j)=a(ia+i)+sin45*(a(ib+i)-a(id+i))
      c(jb+j)=a(ia+i)-sin45*(a(ib+i)-a(id+i))
      d(ja+j)=-a(ic+i)-sin45*(a(ib+i)+a(id+i))
      d(jb+j)=a(ic+i)-sin45*(a(ib+i)+a(id+i))
      i=i+inc3
      j=j+inc4
  470 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  480 continue
      go to 900
c
  490 continue
      z=1.0/float(n)
      do 494 l=1,la
      i=ibase
      j=jbase
cdir$ ivdep
      do 492 ijk=1,lot
      c(ja+j)=z*((a(ia+i)+a(ic+i))+(a(ib+i)+a(id+i)))
      c(jc+j)=z*((a(ia+i)+a(ic+i))-(a(ib+i)+a(id+i)))
      c(jb+j)=z*(a(ia+i)-a(ic+i))
      d(jb+j)=z*(a(id+i)-a(ib+i))
      i=i+inc3
      j=j+inc4
  492 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  494 continue
      go to 900
c
c     coding for factor 5
c     -------------------
  500 continue
      ia=1
      ib=ia+iink
      ic=ib+iink
      id=ic+iink
      ie=id+iink
      ja=1
      jb=ja+(2*m-la)*inc2
      jc=jb+2*m*inc2
      jd=jc
      je=jb
c
      if (la.eq.m) go to 590
c
      do 520 l=1,la
      i=ibase
      j=jbase
cdir$ ivdep
      do 510 ijk=1,lot
      a1=a(ib+i)+a(ie+i)
      a3=a(ib+i)-a(ie+i)
      a2=a(ic+i)+a(id+i)
      a4=a(ic+i)-a(id+i)
      a5=a(ia+i)-0.25*(a1+a2)
      a6=qrt5*(a1-a2)
      c(ja+j)=a(ia+i)+(a1+a2)
      c(jb+j)=a5+a6
      c(jc+j)=a5-a6
      d(jb+j)=-sin72*a3-sin36*a4
      d(jc+j)=-sin36*a3+sin72*a4
      i=i+inc3
      j=j+inc4
  510 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  520 continue
      ja=ja+jink
      jink=2*jink
      jb=jb+jink
      jc=jc+jink
      jd=jd-jink
      je=je-jink
      ibase=ibase+ijump
      ijump=2*ijump+iink
      if (jb.eq.jd) go to 560
      do 550 k=la,kstop,la
      kb=k+k
      kc=kb+kb
      kd=kc+kb
      ke=kd+kb
      c1=trigs(kb+1)
      s1=trigs(kb+2)
      c2=trigs(kc+1)
      s2=trigs(kc+2)
      c3=trigs(kd+1)
      s3=trigs(kd+2)
      c4=trigs(ke+1)
      s4=trigs(ke+2)
      jbase=0
      do 540 l=1,la
      i=ibase
      j=jbase
cdir$ ivdep
      do 530 ijk=1,lot
      a1=(c1*a(ib+i)+s1*b(ib+i))+(c4*a(ie+i)+s4*b(ie+i))
      a3=(c1*a(ib+i)+s1*b(ib+i))-(c4*a(ie+i)+s4*b(ie+i))
      a2=(c2*a(ic+i)+s2*b(ic+i))+(c3*a(id+i)+s3*b(id+i))
      a4=(c2*a(ic+i)+s2*b(ic+i))-(c3*a(id+i)+s3*b(id+i))
      b1=(c1*b(ib+i)-s1*a(ib+i))+(c4*b(ie+i)-s4*a(ie+i))
      b3=(c1*b(ib+i)-s1*a(ib+i))-(c4*b(ie+i)-s4*a(ie+i))
      b2=(c2*b(ic+i)-s2*a(ic+i))+(c3*b(id+i)-s3*a(id+i))
      b4=(c2*b(ic+i)-s2*a(ic+i))-(c3*b(id+i)-s3*a(id+i))
      a5=a(ia+i)-0.25*(a1+a2)
      a6=qrt5*(a1-a2)
      b5=b(ia+i)-0.25*(b1+b2)
      b6=qrt5*(b1-b2)
      a10=a5+a6
      a20=a5-a6
      b10=b5+b6
      b20=b5-b6
      a11=sin72*b3+sin36*b4
      a21=sin36*b3-sin72*b4
      b11=sin72*a3+sin36*a4
      b21=sin36*a3-sin72*a4
      c(ja+j)=a(ia+i)+(a1+a2)
      c(jb+j)=a10+a11
      c(je+j)=a10-a11
      c(jc+j)=a20+a21
      c(jd+j)=a20-a21
      d(ja+j)=b(ia+i)+(b1+b2)
      d(jb+j)=b10-b11
      d(je+j)=-(b10+b11)
      d(jc+j)=b20-b21
      d(jd+j)=-(b20+b21)
      i=i+inc3
      j=j+inc4
  530 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  540 continue
      ibase=ibase+ijump
      ja=ja+jink
      jb=jb+jink
      jc=jc+jink
      jd=jd-jink
      je=je-jink
  550 continue
      if (jb.gt.jd) go to 900
  560 continue
      jbase=0
      do 580 l=1,la
      i=ibase
      j=jbase
cdir$ ivdep
      do 570 ijk=1,lot
      a1=a(ib+i)+a(ie+i)
      a3=a(ib+i)-a(ie+i)
      a2=a(ic+i)+a(id+i)
      a4=a(ic+i)-a(id+i)
      a5=a(ia+i)+0.25*(a3-a4)
      a6=qrt5*(a3+a4)
      c(ja+j)=a5+a6
      c(jb+j)=a5-a6
      c(jc+j)=a(ia+i)-(a3-a4)
      d(ja+j)=-sin36*a1-sin72*a2
      d(jb+j)=-sin72*a1+sin36*a2
      i=i+inc3
      j=j+inc4
  570 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  580 continue
      go to 900
c
  590 continue
      z=1.0/float(n)
      zqrt5=z*qrt5
      zsin36=z*sin36
      zsin72=z*sin72
      do 594 l=1,la
      i=ibase
      j=jbase
cdir$ ivdep
      do 592 ijk=1,lot
      a1=a(ib+i)+a(ie+i)
      a3=a(ib+i)-a(ie+i)
      a2=a(ic+i)+a(id+i)
      a4=a(ic+i)-a(id+i)
      a5=z*(a(ia+i)-0.25*(a1+a2))
      a6=zqrt5*(a1-a2)
      c(ja+j)=z*(a(ia+i)+(a1+a2))
      c(jb+j)=a5+a6
      c(jc+j)=a5-a6
      d(jb+j)=-zsin72*a3-zsin36*a4
      d(jc+j)=-zsin36*a3+zsin72*a4
      i=i+inc3
      j=j+inc4
  592 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  594 continue
      go to 900
c
c     coding for factor 6
c     -------------------
  600 continue
      ia=1
      ib=ia+iink
      ic=ib+iink
      id=ic+iink
      ie=id+iink
      if=ie+iink
      ja=1
      jb=ja+(2*m-la)*inc2
      jc=jb+2*m*inc2
      jd=jc+2*m*inc2
      je=jc
      jf=jb
c
      if (la.eq.m) go to 690
c
      do 620 l=1,la
      i=ibase
      j=jbase
cdir$ ivdep
      do 610 ijk=1,lot
      a11=(a(ic+i)+a(if+i))+(a(ib+i)+a(ie+i))
      c(ja+j)=(a(ia+i)+a(id+i))+a11
      c(jc+j)=(a(ia+i)+a(id+i)-0.5*a11)
      d(jc+j)=sin60*((a(ic+i)+a(if+i))-(a(ib+i)+a(ie+i)))
      a11=(a(ic+i)-a(if+i))+(a(ie+i)-a(ib+i))
      c(jb+j)=(a(ia+i)-a(id+i))-0.5*a11
      d(jb+j)=sin60*((a(ie+i)-a(ib+i))-(a(ic+i)-a(if+i)))
      c(jd+j)=(a(ia+i)-a(id+i))+a11
      i=i+inc3
      j=j+inc4
  610 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  620 continue
      ja=ja+jink
      jink=2*jink
      jb=jb+jink
      jc=jc+jink
      jd=jd-jink
      je=je-jink
      jf=jf-jink
      ibase=ibase+ijump
      ijump=2*ijump+iink
      if (jc.eq.jd) go to 660
      do 650 k=la,kstop,la
      kb=k+k
      kc=kb+kb
      kd=kc+kb
      ke=kd+kb
      kf=ke+kb
      c1=trigs(kb+1)
      s1=trigs(kb+2)
      c2=trigs(kc+1)
      s2=trigs(kc+2)
      c3=trigs(kd+1)
      s3=trigs(kd+2)
      c4=trigs(ke+1)
      s4=trigs(ke+2)
      c5=trigs(kf+1)
      s5=trigs(kf+2)
      jbase=0
      do 640 l=1,la
      i=ibase
      j=jbase
cdir$ ivdep
      do 630 ijk=1,lot
      a1=c1*a(ib+i)+s1*b(ib+i)
      b1=c1*b(ib+i)-s1*a(ib+i)
      a2=c2*a(ic+i)+s2*b(ic+i)
      b2=c2*b(ic+i)-s2*a(ic+i)
      a3=c3*a(id+i)+s3*b(id+i)
      b3=c3*b(id+i)-s3*a(id+i)
      a4=c4*a(ie+i)+s4*b(ie+i)
      b4=c4*b(ie+i)-s4*a(ie+i)
      a5=c5*a(if+i)+s5*b(if+i)
      b5=c5*b(if+i)-s5*a(if+i)
      a11=(a2+a5)+(a1+a4)
      a20=(a(ia+i)+a3)-0.5*a11
      a21=sin60*((a2+a5)-(a1+a4))
      b11=(b2+b5)+(b1+b4)
      b20=(b(ia+i)+b3)-0.5*b11
      b21=sin60*((b2+b5)-(b1+b4))
      c(ja+j)=(a(ia+i)+a3)+a11
      d(ja+j)=(b(ia+i)+b3)+b11
      c(jc+j)=a20-b21
      d(jc+j)=a21+b20
      c(je+j)=a20+b21
      d(je+j)=a21-b20
      a11=(a2-a5)+(a4-a1)
      a20=(a(ia+i)-a3)-0.5*a11
      a21=sin60*((a4-a1)-(a2-a5))
      b11=(b5-b2)-(b4-b1)
      b20=(b3-b(ia+i))-0.5*b11
      b21=sin60*((b5-b2)+(b4-b1))
      c(jb+j)=a20-b21
      d(jb+j)=a21-b20
      c(jd+j)=a11+(a(ia+i)-a3)
      d(jd+j)=b11+(b3-b(ia+i))
      c(jf+j)=a20+b21
      d(jf+j)=a21+b20
      i=i+inc3
      j=j+inc4
  630 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  640 continue
      ibase=ibase+ijump
      ja=ja+jink
      jb=jb+jink
      jc=jc+jink
      jd=jd-jink
      je=je-jink
      jf=jf-jink
  650 continue
      if (jc.gt.jd) go to 900
  660 continue
      jbase=0
      do 680 l=1,la
      i=ibase
      j=jbase
cdir$ ivdep
      do 670 ijk=1,lot
      c(ja+j)=(a(ia+i)+0.5*(a(ic+i)-a(ie+i)))+ sin60*(a(ib+i)-a(if+i))
      d(ja+j)=-(a(id+i)+0.5*(a(ib+i)+a(if+i)))-sin60*(a(ic+i)+a(ie+i))
      c(jb+j)=a(ia+i)-(a(ic+i)-a(ie+i))
      d(jb+j)=a(id+i)-(a(ib+i)+a(if+i))
      c(jc+j)=(a(ia+i)+0.5*(a(ic+i)-a(ie+i)))-sin60*(a(ib+i)-a(if+i))
      d(jc+j)=-(a(id+i)+0.5*(a(ib+i)+a(if+i)))+sin60*(a(ic+i)+a(ie+i))
      i=i+inc3
      j=j+inc4
  670 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  680 continue
      go to 900
c
  690 continue
      z=1.0/float(n)
      zsin60=z*sin60
      do 694 l=1,la
      i=ibase
      j=jbase
cdir$ ivdep
      do 692 ijk=1,lot
      a11=(a(ic+i)+a(if+i))+(a(ib+i)+a(ie+i))
      c(ja+j)=z*((a(ia+i)+a(id+i))+a11)
      c(jc+j)=z*((a(ia+i)+a(id+i))-0.5*a11)
      d(jc+j)=zsin60*((a(ic+i)+a(if+i))-(a(ib+i)+a(ie+i)))
      a11=(a(ic+i)-a(if+i))+(a(ie+i)-a(ib+i))
      c(jb+j)=z*((a(ia+i)-a(id+i))-0.5*a11)
      d(jb+j)=zsin60*((a(ie+i)-a(ib+i))-(a(ic+i)-a(if+i)))
      c(jd+j)=z*((a(ia+i)-a(id+i))+a11)
      i=i+inc3
      j=j+inc4
  692 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  694 continue
      go to 900
c
c     coding for factor 8
c     -------------------
  800 continue
      ibad=3
      if (la.ne.m) go to 910
      ia=1
      ib=ia+iink
      ic=ib+iink
      id=ic+iink
      ie=id+iink
      if=ie+iink
      ig=if+iink
      ih=ig+iink
      ja=1
      jb=ja+la*inc2
      jc=jb+2*m*inc2
      jd=jc+2*m*inc2
      je=jd+2*m*inc2
      z=1.0/float(n)
      zsin45=z*dsqrt(0.5d0)
c
      do 820 l=1,la
      i=ibase
      j=jbase
cdir$ ivdep
      do 810 ijk=1,lot
      c(ja+j)=z*(((a(ia+i)+a(ie+i))+(a(ic+i)+a(ig+i)))+
     *    ((a(id+i)+a(ih+i))+(a(ib+i)+a(if+i))))
      c(je+j)=z*(((a(ia+i)+a(ie+i))+(a(ic+i)+a(ig+i)))-
     *    ((a(id+i)+a(ih+i))+(a(ib+i)+a(if+i))))
      c(jc+j)=z*((a(ia+i)+a(ie+i))-(a(ic+i)+a(ig+i)))
      d(jc+j)=z*((a(id+i)+a(ih+i))-(a(ib+i)+a(if+i)))
      c(jb+j)=z*(a(ia+i)-a(ie+i))
     *    +zsin45*((a(ih+i)-a(id+i))-(a(if+i)-a(ib+i)))
      c(jd+j)=z*(a(ia+i)-a(ie+i))
     *    -zsin45*((a(ih+i)-a(id+i))-(a(if+i)-a(ib+i)))
      d(jb+j)=zsin45*((a(ih+i)-a(id+i))+(a(if+i)-a(ib+i)))
     *    +z*(a(ig+i)-a(ic+i))
      d(jd+j)=zsin45*((a(ih+i)-a(id+i))+(a(if+i)-a(ib+i)))
     *    -z*(a(ig+i)-a(ic+i))
      i=i+inc3
      j=j+inc4
  810 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  820 continue
c
c     return
c     ------
  900 continue
      ibad=0
  910 continue
      ierr=ibad
      return
      end
*deck libfft02
c     subroutine 'rpassm' - performs one pass through data as part
c     of multiple real fft (fourier synthesis) routine
c
c     a is first real input vector
c         equivalence b(1) with a (la*inc1+1)
c     c is first real output vector
c         equivalence d(1) with c(ifac*la*inc2+1)
c     trigs is a precalculated list of sines & cosines
c     inc1 is the addressing increment for a
c     inc2 is the addressing increment for c
c     inc3 is the increment between input vectors a
c     inc4 is the increment between output vectors c
c     lot is the number of vectors
c     n is the length of the vectors
c     ifac is the current factor of n
c     la is the product of previous factors
c     ierr is an error indicator:
c              0 - pass completed without error
c              1 - lot greater than 64
c              2 - ifac not catered for
c              3 - ifac only catered for if la=n/ifac
c
c-----------------------------------------------------------------
c
      subroutine rpassm(a,b,c,d,trigs,inc1,inc2,inc3,inc4,lot,n,ifac
     *    ,la,ierr)
c
      implicit real*8(a-h,o-z)
      dimension a(*),b(*),c(*),d(*),trigs(*)
c
      dimension a10(64),a11(64),a20(64),a21(64),
     *    b10(64),b11(64),b20(64),b21(64)
      data sin36/0.587785252292473/,sin72/0.951056516295154/,
     *    qrt5/0.559016994374947/,sin60/0.866025403784437/
c
c      write(*,'(a8,3i8)')'rpassm:',lot,n,ifac
      m=n/ifac
      iink=la*inc1
      jink=la*inc2
      jump=(ifac-1)*jink
      kstop=(n-ifac)/(2*ifac)
c
      ibad=1
      if (lot.gt.64) go to 910
      ibase=0
      jbase=0
      igo=ifac-1
      if (igo.eq.7) igo=6
      ibad=2
      if (igo.lt.1.or.igo.gt.6) go to 910
      go to (200,300,400,500,600,800),igo
c
c     coding for factor 2
c     -------------------
  200 continue
      ia=1
      ib=ia+(2*m-la)*inc1
      ja=1
      jb=ja+jink
c
      if (la.eq.m) go to 290
c
      do 220 l=1,la
      i=ibase
      j=jbase
cdir$ ivdep
      do 210 ijk=1,lot
      c(ja+j)=a(ia+i)+a(ib+i)
      c(jb+j)=a(ia+i)-a(ib+i)
      i=i+inc3
      j=j+inc4
  210 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  220 continue
      ia=ia+iink
      iink=2*iink
      ib=ib-iink
      ibase=0
      jbase=jbase+jump
      jump=2*jump+jink
      if (ia.eq.ib) go to 260
      do 250 k=la,kstop,la
      kb=k+k
      c1=trigs(kb+1)
      s1=trigs(kb+2)
      ibase=0
      do 240 l=1,la
      i=ibase
      j=jbase
cdir$ ivdep
      do 230 ijk=1,lot
      c(ja+j)=a(ia+i)+a(ib+i)
      d(ja+j)=b(ia+i)-b(ib+i)
      c(jb+j)=c1*(a(ia+i)-a(ib+i))-s1*(b(ia+i)+b(ib+i))
      d(jb+j)=s1*(a(ia+i)-a(ib+i))+c1*(b(ia+i)+b(ib+i))
      i=i+inc3
      j=j+inc4
  230 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  240 continue
      ia=ia+iink
      ib=ib-iink
      jbase=jbase+jump
  250 continue
      if (ia.gt.ib) go to 900
  260 continue
      ibase=0
      do 280 l=1,la
      i=ibase
      j=jbase
cdir$ ivdep
      do 270 ijk=1,lot
      c(ja+j)=a(ia+i)
      c(jb+j)=-b(ia+i)
      i=i+inc3
      j=j+inc4
  270 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  280 continue
      go to 900
c
  290 continue
      do 294 l=1,la
      i=ibase
      j=jbase
cdir$ ivdep
      do 292 ijk=1,lot
      c(ja+j)=2.0*(a(ia+i)+a(ib+i))
      c(jb+j)=2.0*(a(ia+i)-a(ib+i))
      i=i+inc3
      j=j+inc4
  292 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  294 continue
      go to 900
c
c     coding for factor 3
c     -------------------
  300 continue
      ia=1
      ib=ia+(2*m-la)*inc1
      ic=ib
      ja=1
      jb=ja+jink
      jc=jb+jink
c
      if (la.eq.m) go to 390
c
      do 320 l=1,la
      i=ibase
      j=jbase
cdir$ ivdep
      do 310 ijk=1,lot
      c(ja+j)=a(ia+i)+a(ib+i)
      c(jb+j)=(a(ia+i)-0.5*a(ib+i))-(sin60*(b(ib+i)))
      c(jc+j)=(a(ia+i)-0.5*a(ib+i))+(sin60*(b(ib+i)))
      i=i+inc3
      j=j+inc4
  310 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  320 continue
      ia=ia+iink
      iink=2*iink
      ib=ib+iink
      ic=ic-iink
      jbase=jbase+jump
      jump=2*jump+jink
      if (ia.eq.ic) go to 360
      do 350 k=la,kstop,la
      kb=k+k
      kc=kb+kb
      c1=trigs(kb+1)
      s1=trigs(kb+2)
      c2=trigs(kc+1)
      s2=trigs(kc+2)
      ibase=0
      do 340 l=1,la
      i=ibase
      j=jbase
cdir$ ivdep
      do 330 ijk=1,lot
      c(ja+j)=a(ia+i)+(a(ib+i)+a(ic+i))
      d(ja+j)=b(ia+i)+(b(ib+i)-b(ic+i))
      c(jb+j)=
     * c1*((a(ia+i)-0.5*(a(ib+i)+a(ic+i)))-(sin60*(b(ib+i)+b(ic+i))))
     * -s1*((b(ia+i)-0.5*(b(ib+i)-b(ic+i)))+(sin60*(a(ib+i)-a(ic+i))))
      d(jb+j)=
     * s1*((a(ia+i)-0.5*(a(ib+i)+a(ic+i)))-(sin60*(b(ib+i)+b(ic+i))))
     * +c1*((b(ia+i)-0.5*(b(ib+i)-b(ic+i)))+(sin60*(a(ib+i)-a(ic+i))))
      c(jc+j)=
     * c2*((a(ia+i)-0.5*(a(ib+i)+a(ic+i)))+(sin60*(b(ib+i)+b(ic+i))))
     * -s2*((b(ia+i)-0.5*(b(ib+i)-b(ic+i)))-(sin60*(a(ib+i)-a(ic+i))))
      d(jc+j)=
     * s2*((a(ia+i)-0.5*(a(ib+i)+a(ic+i)))+(sin60*(b(ib+i)+b(ic+i))))
     * +c2*((b(ia+i)-0.5*(b(ib+i)-b(ic+i)))-(sin60*(a(ib+i)-a(ic+i))))
      i=i+inc3
      j=j+inc4
  330 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  340 continue
      ia=ia+iink
      ib=ib+iink
      ic=ic-iink
      jbase=jbase+jump
  350 continue
      if (ia.gt.ic) go to 900
  360 continue
      ibase=0
      do 380 l=1,la
      i=ibase
      j=jbase
cdir$ ivdep
      do 370 ijk=1,lot
      c(ja+j)=a(ia+i)+a(ib+i)
      c(jb+j)=(0.5*a(ia+i)-a(ib+i))-(sin60*b(ia+i))
      c(jc+j)=-(0.5*a(ia+i)-a(ib+i))-(sin60*b(ia+i))
      i=i+inc3
      j=j+inc4
  370 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  380 continue
      go to 900
c
  390 continue
      ssin60=2.0*sin60
      do 394 l=1,la
      i=ibase
      j=jbase
cdir$ ivdep
      do 392 ijk=1,lot
      c(ja+j)=2.0*(a(ia+i)+a(ib+i))
      c(jb+j)=(2.0*a(ia+i)-a(ib+i))-(ssin60*b(ib+i))
      c(jc+j)=(2.0*a(ia+i)-a(ib+i))+(ssin60*b(ib+i))
      i=i+inc3
      j=j+inc4
  392 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  394 continue
      go to 900
c
c     coding for factor 4
c     -------------------
  400 continue
      ia=1
      ib=ia+(2*m-la)*inc1
      ic=ib+2*m*inc1
      id=ib
      ja=1
      jb=ja+jink
      jc=jb+jink
      jd=jc+jink
c
      if (la.eq.m) go to 490
c
      do 420 l=1,la
      i=ibase
      j=jbase
cdir$ ivdep
      do 410 ijk=1,lot
      c(ja+j)=(a(ia+i)+a(ic+i))+a(ib+i)
      c(jb+j)=(a(ia+i)-a(ic+i))-b(ib+i)
      c(jc+j)=(a(ia+i)+a(ic+i))-a(ib+i)
      c(jd+j)=(a(ia+i)-a(ic+i))+b(ib+i)
      i=i+inc3
      j=j+inc4
  410 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  420 continue
      ia=ia+iink
      iink=2*iink
      ib=ib+iink
      ic=ic-iink
      id=id-iink
      jbase=jbase+jump
      jump=2*jump+jink
      if (ib.eq.ic) go to 460
      do 450 k=la,kstop,la
      kb=k+k
      kc=kb+kb
      kd=kc+kb
      c1=trigs(kb+1)
      s1=trigs(kb+2)
      c2=trigs(kc+1)
      s2=trigs(kc+2)
      c3=trigs(kd+1)
      s3=trigs(kd+2)
      ibase=0
      do 440 l=1,la
      i=ibase
      j=jbase
cdir$ ivdep
      do 430 ijk=1,lot
      c(ja+j)=(a(ia+i)+a(ic+i))+(a(ib+i)+a(id+i))
      d(ja+j)=(b(ia+i)-b(ic+i))+(b(ib+i)-b(id+i))
      c(jc+j)=
     *    c2*((a(ia+i)+a(ic+i))-(a(ib+i)+a(id+i)))
     *   -s2*((b(ia+i)-b(ic+i))-(b(ib+i)-b(id+i)))
      d(jc+j)=
     *    s2*((a(ia+i)+a(ic+i))-(a(ib+i)+a(id+i)))
     *   +c2*((b(ia+i)-b(ic+i))-(b(ib+i)-b(id+i)))
      c(jb+j)=
     *    c1*((a(ia+i)-a(ic+i))-(b(ib+i)+b(id+i)))
     *   -s1*((b(ia+i)+b(ic+i))+(a(ib+i)-a(id+i)))
      d(jb+j)=
     *    s1*((a(ia+i)-a(ic+i))-(b(ib+i)+b(id+i)))
     *   +c1*((b(ia+i)+b(ic+i))+(a(ib+i)-a(id+i)))
      c(jd+j)=
     *    c3*((a(ia+i)-a(ic+i))+(b(ib+i)+b(id+i)))
     *   -s3*((b(ia+i)+b(ic+i))-(a(ib+i)-a(id+i)))
      d(jd+j)=
     *    s3*((a(ia+i)-a(ic+i))+(b(ib+i)+b(id+i)))
     *   +c3*((b(ia+i)+b(ic+i))-(a(ib+i)-a(id+i)))
      i=i+inc3
      j=j+inc4
  430 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  440 continue
      ia=ia+iink
      ib=ib+iink
      ic=ic-iink
      id=id-iink
      jbase=jbase+jump
  450 continue
      if (ib.gt.ic) go to 900
  460 continue
      ibase=0
      sin45=dsqrt(0.5d0)
      do 480 l=1,la
      i=ibase
      j=jbase
cdir$ ivdep
      do 470 ijk=1,lot
      c(ja+j)=a(ia+i)+a(ib+i)
      c(jb+j)=sin45*((a(ia+i)-a(ib+i))-(b(ia+i)+b(ib+i)))
      c(jc+j)=b(ib+i)-b(ia+i)
      c(jd+j)=-sin45*((a(ia+i)-a(ib+i))+(b(ia+i)+b(ib+i)))
      i=i+inc3
      j=j+inc4
  470 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  480 continue
      go to 900
c
  490 continue
      do 494 l=1,la
      i=ibase
      j=jbase
cdir$ ivdep
      do 492 ijk=1,lot
      c(ja+j)=2.0*((a(ia+i)+a(ic+i))+a(ib+i))
      c(jb+j)=2.0*((a(ia+i)-a(ic+i))-b(ib+i))
      c(jc+j)=2.0*((a(ia+i)+a(ic+i))-a(ib+i))
      c(jd+j)=2.0*((a(ia+i)-a(ic+i))+b(ib+i))
      i=i+inc3
      j=j+inc4
  492 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  494 continue
      go to 900
c
c     coding for factor 5
c     -------------------
  500 continue
      ia=1
      ib=ia+(2*m-la)*inc1
      ic=ib+2*m*inc1
      id=ic
      ie=ib
      ja=1
      jb=ja+jink
      jc=jb+jink
      jd=jc+jink
      je=jd+jink
c
      if (la.eq.m) go to 590
c
      do 520 l=1,la
      i=ibase
      j=jbase
cdir$ ivdep
      do 510 ijk=1,lot
      c(ja+j)=a(ia+i)+(a(ib+i)+a(ic+i))
      c(jb+j)=((a(ia+i)-0.25*(a(ib+i)+a(ic+i)))+qrt5*(a(ib+i)-a(ic+i)))
     *    -(sin72*b(ib+i)+sin36*b(ic+i))
      c(jc+j)=((a(ia+i)-0.25*(a(ib+i)+a(ic+i)))-qrt5*(a(ib+i)-a(ic+i)))
     *    -(sin36*b(ib+i)-sin72*b(ic+i))
      c(jd+j)=((a(ia+i)-0.25*(a(ib+i)+a(ic+i)))-qrt5*(a(ib+i)-a(ic+i)))
     *    +(sin36*b(ib+i)-sin72*b(ic+i))
      c(je+j)=((a(ia+i)-0.25*(a(ib+i)+a(ic+i)))+qrt5*(a(ib+i)-a(ic+i)))
     *    +(sin72*b(ib+i)+sin36*b(ic+i))
      i=i+inc3
      j=j+inc4
  510 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  520 continue
      ia=ia+iink
      iink=2*iink
      ib=ib+iink
      ic=ic+iink
      id=id-iink
      ie=ie-iink
      jbase=jbase+jump
      jump=2*jump+jink
      if (ib.eq.id) go to 560
      do 550 k=la,kstop,la
      kb=k+k
      kc=kb+kb
      kd=kc+kb
      ke=kd+kb
      c1=trigs(kb+1)
      s1=trigs(kb+2)
      c2=trigs(kc+1)
      s2=trigs(kc+2)
      c3=trigs(kd+1)
      s3=trigs(kd+2)
      c4=trigs(ke+1)
      s4=trigs(ke+2)
      ibase=0
      do 540 l=1,la
      i=ibase
      j=jbase
cdir$ ivdep
      do 530 ijk=1,lot
c
      a10(ijk)=(a(ia+i)-0.25*((a(ib+i)+a(ie+i))+(a(ic+i)+a(id+i))))
     *    +qrt5*((a(ib+i)+a(ie+i))-(a(ic+i)+a(id+i)))
      a20(ijk)=(a(ia+i)-0.25*((a(ib+i)+a(ie+i))+(a(ic+i)+a(id+i))))
     *    -qrt5*((a(ib+i)+a(ie+i))-(a(ic+i)+a(id+i)))
      b10(ijk)=(b(ia+i)-0.25*((b(ib+i)-b(ie+i))+(b(ic+i)-b(id+i))))
     *    +qrt5*((b(ib+i)-b(ie+i))-(b(ic+i)-b(id+i)))
      b20(ijk)=(b(ia+i)-0.25*((b(ib+i)-b(ie+i))+(b(ic+i)-b(id+i))))
     *    -qrt5*((b(ib+i)-b(ie+i))-(b(ic+i)-b(id+i)))
      a11(ijk)=sin72*(b(ib+i)+b(ie+i))+sin36*(b(ic+i)+b(id+i))
      a21(ijk)=sin36*(b(ib+i)+b(ie+i))-sin72*(b(ic+i)+b(id+i))
      b11(ijk)=sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i))
      b21(ijk)=sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i))
c
      c(ja+j)=a(ia+i)+((a(ib+i)+a(ie+i))+(a(ic+i)+a(id+i)))
      d(ja+j)=b(ia+i)+((b(ib+i)-b(ie+i))+(b(ic+i)-b(id+i)))
      c(jb+j)=c1*(a10(ijk)-a11(ijk))-s1*(b10(ijk)+b11(ijk))
      d(jb+j)=s1*(a10(ijk)-a11(ijk))+c1*(b10(ijk)+b11(ijk))
      c(je+j)=c4*(a10(ijk)+a11(ijk))-s4*(b10(ijk)-b11(ijk))
      d(je+j)=s4*(a10(ijk)+a11(ijk))+c4*(b10(ijk)-b11(ijk))
      c(jc+j)=c2*(a20(ijk)-a21(ijk))-s2*(b20(ijk)+b21(ijk))
      d(jc+j)=s2*(a20(ijk)-a21(ijk))+c2*(b20(ijk)+b21(ijk))
      c(jd+j)=c3*(a20(ijk)+a21(ijk))-s3*(b20(ijk)-b21(ijk))
      d(jd+j)=s3*(a20(ijk)+a21(ijk))+c3*(b20(ijk)-b21(ijk))
c
      i=i+inc3
      j=j+inc4
  530 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  540 continue
      ia=ia+iink
      ib=ib+iink
      ic=ic+iink
      id=id-iink
      ie=ie-iink
      jbase=jbase+jump
  550 continue
      if (ib.gt.id) go to 900
  560 continue
      ibase=0
      do 580 l=1,la
      i=ibase
      j=jbase
cdir$ ivdep
      do 570 ijk=1,lot
      c(ja+j)=(a(ia+i)+a(ib+i))+a(ic+i)
      c(jb+j)=(qrt5*(a(ia+i)-a(ib+i))+(0.25*(a(ia+i)+a(ib+i))-a(ic+i)))
     *    -(sin36*b(ia+i)+sin72*b(ib+i))
      c(je+j)=-(qrt5*(a(ia+i)-a(ib+i))+(0.25*(a(ia+i)+a(ib+i))-a(ic+i)
     * ))
     *    -(sin36*b(ia+i)+sin72*b(ib+i))
      c(jc+j)=(qrt5*(a(ia+i)-a(ib+i))-(0.25*(a(ia+i)+a(ib+i))-a(ic+i)))
     *    -(sin72*b(ia+i)-sin36*b(ib+i))
      c(jd+j)=-(qrt5*(a(ia+i)-a(ib+i))-(0.25*(a(ia+i)+a(ib+i))-a(ic+i)
     * ))
     *    -(sin72*b(ia+i)-sin36*b(ib+i))
      i=i+inc3
      j=j+inc4
  570 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  580 continue
      go to 900
c
  590 continue
      qqrt5=2.0*qrt5
      ssin36=2.0*sin36
      ssin72=2.0*sin72
      do 594 l=1,la
      i=ibase
      j=jbase
cdir$ ivdep
      do 592 ijk=1,lot
      c(ja+j)=2.0*(a(ia+i)+(a(ib+i)+a(ic+i)))
      c(jb+j)=(2.0*(a(ia+i)-0.25*(a(ib+i)+a(ic+i)))
     *    +qqrt5*(a(ib+i)-a(ic+i)))-(ssin72*b(ib+i)+ssin36*b(ic+i))
      c(jc+j)=(2.0*(a(ia+i)-0.25*(a(ib+i)+a(ic+i)))
     *    -qqrt5*(a(ib+i)-a(ic+i)))-(ssin36*b(ib+i)-ssin72*b(ic+i))
      c(jd+j)=(2.0*(a(ia+i)-0.25*(a(ib+i)+a(ic+i)))
     *    -qqrt5*(a(ib+i)-a(ic+i)))+(ssin36*b(ib+i)-ssin72*b(ic+i))
      c(je+j)=(2.0*(a(ia+i)-0.25*(a(ib+i)+a(ic+i)))
     *    +qqrt5*(a(ib+i)-a(ic+i)))+(ssin72*b(ib+i)+ssin36*b(ic+i))
      i=i+inc3
      j=j+inc4
  592 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  594 continue
      go to 900
c
c     coding for factor 6
c     -------------------
  600 continue
      ia=1
      ib=ia+(2*m-la)*inc1
      ic=ib+2*m*inc1
      id=ic+2*m*inc1
      ie=ic
      if=ib
      ja=1
      jb=ja+jink
      jc=jb+jink
      jd=jc+jink
      je=jd+jink
      jf=je+jink
c
      if (la.eq.m) go to 690
c
      do 620 l=1,la
      i=ibase
      j=jbase
cdir$ ivdep
      do 610 ijk=1,lot
      c(ja+j)=(a(ia+i)+a(id+i))+(a(ib+i)+a(ic+i))
      c(jd+j)=(a(ia+i)-a(id+i))-(a(ib+i)-a(ic+i))
      c(jb+j)=((a(ia+i)-a(id+i))+0.5*(a(ib+i)-a(ic+i)))
     *    -(sin60*(b(ib+i)+b(ic+i)))
      c(jf+j)=((a(ia+i)-a(id+i))+0.5*(a(ib+i)-a(ic+i)))
     *    +(sin60*(b(ib+i)+b(ic+i)))
      c(jc+j)=((a(ia+i)+a(id+i))-0.5*(a(ib+i)+a(ic+i)))
     *    -(sin60*(b(ib+i)-b(ic+i)))
      c(je+j)=((a(ia+i)+a(id+i))-0.5*(a(ib+i)+a(ic+i)))
     *    +(sin60*(b(ib+i)-b(ic+i)))
      i=i+inc3
      j=j+inc4
  610 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  620 continue
      ia=ia+iink
      iink=2*iink
      ib=ib+iink
      ic=ic+iink
      id=id-iink
      ie=ie-iink
      if=if-iink
      jbase=jbase+jump
      jump=2*jump+jink
      if (ic.eq.id) go to 660
      do 650 k=la,kstop,la
      kb=k+k
      kc=kb+kb
      kd=kc+kb
      ke=kd+kb
      kf=ke+kb
      c1=trigs(kb+1)
      s1=trigs(kb+2)
      c2=trigs(kc+1)
      s2=trigs(kc+2)
      c3=trigs(kd+1)
      s3=trigs(kd+2)
      c4=trigs(ke+1)
      s4=trigs(ke+2)
      c5=trigs(kf+1)
      s5=trigs(kf+2)
      ibase=0
      do 640 l=1,la
      i=ibase
      j=jbase
cdir$ ivdep
      do 630 ijk=1,lot
c
      a11(ijk)= (a(ie+i)+a(ib+i))+(a(ic+i)+a(if+i))
      a20(ijk)=(a(ia+i)+a(id+i))-0.5*a11(ijk)
      a21(ijk)=sin60*((a(ie+i)+a(ib+i))-(a(ic+i)+a(if+i)))
      b11(ijk)= (b(ib+i)-b(ie+i))+(b(ic+i)-b(if+i))
      b20(ijk)=(b(ia+i)-b(id+i))-0.5*b11(ijk)
      b21(ijk)=sin60*((b(ib+i)-b(ie+i))-(b(ic+i)-b(if+i)))
c
      c(ja+j)=(a(ia+i)+a(id+i))+a11(ijk)
      d(ja+j)=(b(ia+i)-b(id+i))+b11(ijk)
      c(jc+j)=c2*(a20(ijk)-b21(ijk))-s2*(b20(ijk)+a21(ijk))
      d(jc+j)=s2*(a20(ijk)-b21(ijk))+c2*(b20(ijk)+a21(ijk))
      c(je+j)=c4*(a20(ijk)+b21(ijk))-s4*(b20(ijk)-a21(ijk))
      d(je+j)=s4*(a20(ijk)+b21(ijk))+c4*(b20(ijk)-a21(ijk))
c
      a11(ijk)=(a(ie+i)-a(ib+i))+(a(ic+i)-a(if+i))
      b11(ijk)=(b(ie+i)+b(ib+i))-(b(ic+i)+b(if+i))
      a20(ijk)=(a(ia+i)-a(id+i))-0.5*a11(ijk)
      a21(ijk)=sin60*((a(ie+i)-a(ib+i))-(a(ic+i)-a(if+i)))
      b20(ijk)=(b(ia+i)+b(id+i))+0.5*b11(ijk)
      b21(ijk)=sin60*((b(ie+i)+b(ib+i))+(b(ic+i)+b(if+i)))
c
      c(jd+j)=
     *  c3*((a(ia+i)-a(id+i))+a11(ijk))-s3*((b(ia+i)+b(id+i))-b11(ijk))
      d(jd+j)=
     *  s3*((a(ia+i)-a(id+i))+a11(ijk))+c3*((b(ia+i)+b(id+i))-b11(ijk))
      c(jb+j)=c1*(a20(ijk)-b21(ijk))-s1*(b20(ijk)-a21(ijk))
      d(jb+j)=s1*(a20(ijk)-b21(ijk))+c1*(b20(ijk)-a21(ijk))
      c(jf+j)=c5*(a20(ijk)+b21(ijk))-s5*(b20(ijk)+a21(ijk))
      d(jf+j)=s5*(a20(ijk)+b21(ijk))+c5*(b20(ijk)+a21(ijk))
c
      i=i+inc3
      j=j+inc4
  630 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  640 continue
      ia=ia+iink
      ib=ib+iink
      ic=ic+iink
      id=id-iink
      ie=ie-iink
      if=if-iink
      jbase=jbase+jump
  650 continue
      if (ic.gt.id) go to 900
  660 continue
      ibase=0
      do 680 l=1,la
      i=ibase
      j=jbase
cdir$ ivdep
      do 670 ijk=1,lot
      c(ja+j)=a(ib+i)+(a(ia+i)+a(ic+i))
      c(jd+j)=b(ib+i)-(b(ia+i)+b(ic+i))
      c(jb+j)=(sin60*(a(ia+i)-a(ic+i)))-(0.5*(b(ia+i)+b(ic+i))+b(ib+i))
      c(jf+j)=-(sin60*(a(ia+i)-a(ic+i)))-(0.5*(b(ia+i)+b(ic+i))+b(ib+i)
     *)
      c(jc+j)=sin60*(b(ic+i)-b(ia+i))+(0.5*(a(ia+i)+a(ic+i))-a(ib+i))
      c(je+j)=sin60*(b(ic+i)-b(ia+i))-(0.5*(a(ia+i)+a(ic+i))-a(ib+i))
      i=i+inc3
      j=j+inc4
  670 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  680 continue
      go to 900
c
  690 continue
      ssin60=2.0*sin60
      do 694 l=1,la
      i=ibase
      j=jbase
cdir$ ivdep
      do 692 ijk=1,lot
      c(ja+j)=(2.0*(a(ia+i)+a(id+i)))+(2.0*(a(ib+i)+a(ic+i)))
      c(jd+j)=(2.0*(a(ia+i)-a(id+i)))-(2.0*(a(ib+i)-a(ic+i)))
      c(jb+j)=(2.0*(a(ia+i)-a(id+i))+(a(ib+i)-a(ic+i)))
     *    -(ssin60*(b(ib+i)+b(ic+i)))
      c(jf+j)=(2.0*(a(ia+i)-a(id+i))+(a(ib+i)-a(ic+i)))
     *    +(ssin60*(b(ib+i)+b(ic+i)))
      c(jc+j)=(2.0*(a(ia+i)+a(id+i))-(a(ib+i)+a(ic+i)))
     *    -(ssin60*(b(ib+i)-b(ic+i)))
      c(je+j)=(2.0*(a(ia+i)+a(id+i))-(a(ib+i)+a(ic+i)))
     *    +(ssin60*(b(ib+i)-b(ic+i)))
      i=i+inc3
      j=j+inc4
  692 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  694 continue
      go to 900
c
c     coding for factor 8
c     -------------------
  800 continue
      ibad=3
      if (la.ne.m) go to 910
      ia=1
      ib=ia+la*inc1
      ic=ib+2*la*inc1
      id=ic+2*la*inc1
      ie=id+2*la*inc1
      ja=1
      jb=ja+jink
      jc=jb+jink
      jd=jc+jink
      je=jd+jink
      jf=je+jink
      jg=jf+jink
      jh=jg+jink
      ssin45=dsqrt(2.0d0)
c
      do 820 l=1,la
      i=ibase
      j=jbase
cdir$ ivdep
      do 810 ijk=1,lot
      c(ja+j)=2.0*(((a(ia+i)+a(ie+i))+a(ic+i))+(a(ib+i)+a(id+i)))
      c(je+j)=2.0*(((a(ia+i)+a(ie+i))+a(ic+i))-(a(ib+i)+a(id+i)))
      c(jc+j)=2.0*(((a(ia+i)+a(ie+i))-a(ic+i))-(b(ib+i)-b(id+i)))
      c(jg+j)=2.0*(((a(ia+i)+a(ie+i))-a(ic+i))+(b(ib+i)-b(id+i)))
      c(jb+j)=2.0*((a(ia+i)-a(ie+i))-b(ic+i))
     *    +ssin45*((a(ib+i)-a(id+i))-(b(ib+i)+b(id+i)))
      c(jf+j)=2.0*((a(ia+i)-a(ie+i))-b(ic+i))
     *    -ssin45*((a(ib+i)-a(id+i))-(b(ib+i)+b(id+i)))
      c(jd+j)=2.0*((a(ia+i)-a(ie+i))+b(ic+i))
     *    -ssin45*((a(ib+i)-a(id+i))+(b(ib+i)+b(id+i)))
      c(jh+j)=2.0*((a(ia+i)-a(ie+i))+b(ic+i))
     *    +ssin45*((a(ib+i)-a(id+i))+(b(ib+i)+b(id+i)))
      i=i+inc3
      j=j+inc4
  810 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  820 continue
c
c     return
c     ------
  900 continue
      ibad=0
  910 continue
      ierr=ibad
      return
      end


        subroutine calclmon(za,zo,zoh,ua,tvia,tvis,xlmon,ustar,tstar
     :    ,valnon)


c******************************************************************
c
c       Routine de calcul de la longueur de Monin-Obukhov, de ustar
c       tetastar par la methode itierative decrite dans Itier,1980
c
c       Arguments:
c               za      : niveau de mesure de reference
c               zo,zoh  : Longueurs de rugosite des moments et de
c                       la chaleur
c               ua      : Vent a za
c               tvia,tvis: Temperatures virtuelles a za et en surface
c               xlmon   : Longueur de Monin-Obukhov
c               ustar,tstar: Vitesse et temperature de frottement
c               valnon  : valeur manquante
c********************************************************************
        implicit real*8 (a-h,o-z)
        nitermax=30
        eps=1.e-5
c       zua=za+8.
        zua=za
        zm=zua/zo
        zh=za/zoh
        xlmon1=1e6
        xlmon=1e6
        xk=.4
        dift=tvia-tvis
        tmoy=(tvia+tvis)/2.
        gfxk=xk*9.80665


        do 100 n=1,nitermax
        xsi=zua/xlmon1
        ustar=(xk*ua)/(dlog(zm)-psim(xsi))
        xsi=za/xlmon1
        tstar=(xk*dift)/(dlog(zh)-psih(xsi))
        xlmon=(tmoy*ustar*ustar)/(gfxk*tstar)


c       Test de convergence


        if(abs(xlmon).lt.1.) then
                coef=abs(xlmon-xlmon1)
                if(coef.lt..05) then
                        if(abs(xlmon).gt.10000.) then
                                xlmon=dsign(1.d0,dift)*10000.
                        endif
                go to 110
                endif
                else
                coef=abs((xlmon-xlmon1)/xlmon)
                if(coef.lt..01) then
                        if(abs(xlmon).gt.10000.) then
                                xlmon=dsign(1.d0,dift)*10000.
                        endif
                go to 110
                endif
        endif
        xlmon1=xlmon
100     continue


c       Cas tres stable ou le processus ne converge pas


        xlmon=valnon
        ustar=valnon
        tstar=valnon
c       write(5,*) xlmon,ustar,tstar


110     continue
        if(abs(ustar-valnon).ge.eps) then
                ustar=max(ustar,.003d0)
c               if(abs(ustar-.003).lt.eps) write(5,*) 122
        endif
        return
        end






       subroutine rad(rg,rat,pl,xrg,xra,xpl
     :   ,nx,ny,ndados,nf,dtl,n,xugeos,xvgeos,ugeos,vgeos,ns,zonirr)
c
c updates radiation forcing at the surface (linear time interpolation)
c
       implicit double precision (a-h,o-z)
       dimension rg(0:nx+1,0:nx+1),rat(0:nx+1,0:ny+1),pl(0:nx+1,0:ny+1)
     :   ,xrg(0:ndados),xra(0:ndados),xpl(0:ndados)
     :   ,zonirr(0:nx+1,0:ny+1)
       dimension  xugeos(ns,5),xvgeos(ns,5),ugeos(ns),vgeos(ns)


        izx = int(n * dtl/ nf + 1.)
        izzx = int(n * dtl - nf * izx + nf)
c       write(*,*) 'izx=',izx,ndados
        if(izx.ge.ndados) then
          write(0,*) 'nh3dFatalError=Not enough radiation data !'
          stop
        endif
        do ix=1,nx
          do iy=1,ny
            rat(ix,iy)=xra(izx)+(xra(izx+1)-xra(izx))*izzx/nf
            pl(ix,iy)=(xpl(izx)+(xpl(izx+1)-xpl(izx))*izzx/nf)
     :         *zonirr(ix,iy)
            rg(ix,iy)=xrg(izx)+(xrg(izx+1)-xrg(izx))*izzx/nf
          enddo
        enddo


c
c  update of geostrophic wind (note only for fcori.ne.0)
c
c       izx = int(n * dtl/12./ nf + 1.)
c       izzx = int(n * dtl - nf*12 * izx + nf*12)
c       do 104 is=1,ns-1
c          ugeos(is)=xugeos(is,izx)+(xugeos(is,izx+1)-xugeos(is,izx))
c    &     *izzx/12./nf
c          vgeos(is)=xvgeos(is,izx)+(xvgeos(is,izx+1)-xvgeos(is,izx))
c    &     *izzx/12./nf
c       write(28,*) 'ug,uv',ugeos(is),vgeos(is)
c104    continue


c       call wrigar(rat,1,nx,1,ny,0,0,1,nx,1,ny,0,0,'rat     ',0.d0,1)
c       call wrigar(rg,1,nx,1,ny,0,0,1,nx,1,ny,0,0,'rg      ',0.d0,1)
c       call wrigar(pl,1,nx,1,ny,0,0,1,nx,1,ny,0,0,'pl      ',0.d0,1)
c
c      do igeo=1,ns-1
c          write(*,*) 'vgeo',igeo,ugeos(igeo),vgeos(igeo)
c       enddo


c
        return
        end

!-------------NEW subroutines (V. Stepanenko, May 2007): rad_par2, eps1, epsil, sin_sun

      SUBROUTINE rad_par2(taixiy,qa,pres,tclds,atm_rad,net_solar)
!-----Computes solar and atmospheric radiation according to simple
!-----parameterizations      
	use alloc
	implicit none
      real(8), external:: eps1,ndaym,epsil
	real(8) sigma,coef_l,tclds,atm_rad,net_solar,sol_const,
     & taixiy,qa,hume,pres,days,delta,tha,sinh0,phi1,h0,
     & coef_s,a_ang,d_ang,c_ang,a_brent,b_brent	 
      REAL, EXTERNAL:: sin_sun
	integer(4) i
	logical firstcall
	data firstcall /.true./
	SAVE
      	
	if (firstcall) then
	 sigma = 5.67d-8
	 coef_l = 0.22
	 sol_const = 1376.
	 coef_s = 1. - 0.4393
!	 tau_rg = 2.74 !3.5 !0.4 
	 
	 a_ang = 0.82
	 d_ang = 0.25
	 c_ang = 0.95/100.

       a_brent = 0.526
	 b_brent = 0.065/10.
      endif
		
	hume = qa*pres/0.622
	select case(rat_par)	
	case(1)
	 atm_rad = taixiy**4*sigma*eps1(taixiy,hume)
	case(2)
	 atm_rad=sigma*taixiy**4*(a_ang-d_ang*10.**(-c_ang*hume))
	case(3)
       atm_rad=sigma*taixiy**4*(a_brent+b_brent*dsqrt(hume))
	end select
	atm_rad = atm_rad*(1+coef_l*tclds*tclds)
	
	sinh0 = sin_sun(imonth,iday,ihour,iminu,iseco,xlatit)		
	h0 = dasin(sinh0)
	net_solar = sol_const*(1-coef_s*tclds)*sinh0/
     & (1+epsil(h0)*tau_rg/sinh0)
      if (sinh0 < 0) net_solar = 0
	
	if (firstcall) firstcall=.false.
	RETURN
	END SUBROUTINE rad_par2


	REAL(8) FUNCTION eps1(t,h)
	implicit none
	real(8) t,h,c
	if(t.ge.273.15) then
	 c = 0.14
	else
	 c = 0.17
	end if 
	if (h==0) h=0.1
!	eps1 = c*(h**(1./7.))*exp(350./t)
	eps1 = c*((h/100.)**(1./7.))*exp(350./t)
	RETURN
	END FUNCTION eps1


	REAL(8) FUNCTION epsil(h0)
	implicit none
	real(8) h0, pi
      pi = 3.1415926
	   
	if (h0 >= 15./180.*pi.and.h0 < 20./180.*pi) epsil = 0.24
	if (h0 >= 20./180.*pi.and.h0 < 25./180.*pi) epsil = 0.22
	if (h0 >= 25./180.*pi.and.h0 < 30./180.*pi) epsil = 0.20
	if (h0 >= 30./180.*pi.and.h0 < 40./180.*pi) epsil = 0.18
	if (h0 >= 40./180.*pi.and.h0 < 50./180.*pi) epsil = 0.16
	if (h0 >= 50./180.*pi.and.h0 <= 60./180.*pi) epsil = 0.14
	if (h0 < 15./180.*pi) epsil = 0.30
		   
	RETURN
	END FUNCTION epsil     


      REAL FUNCTION sin_sun(imonth,iday,ihour,iminu,iseco,xlatit)
!-----Computes sine of sun height angle (cosine of zenith angle)---------      
      implicit none
      real(8) month,day,hour,minute,second,ndaym(12),days,delta,pi,
     & tha,phi1,xlatit 
      integer i,imonth,iday,ihour,iminu,iseco
      data ndaym /31.,28.,31.,30.,31.,30.,31.,
     &            31.,30.,31.,30.,31./
      
      pi = 4.*atan(1.)
      phi1 = xlatit*pi/180.
      month = float(imonth)
      day   = float(iday)
      hour  = float(ihour)
      minute = float(iminu)
      second = float(iseco)
           
      if (month/=1.) then
	 days=0.
	 do i=1,int4(month)-1
        days = days + ndaym(i)
	 enddo
	 days = days + day
      else
       days = day
	endif
	days = days + int4(hour/24.)
	
	delta = 23.5*pi/180.*dCOS(2*pi*(days-173.)/365.)
	tha   = 2*pi*(hour+minute/60.+second/3600.-12.)/24.
      sin_sun=dSIN(phi1)*dSIN(delta)+dCOS(phi1)*dCOS(delta)*dCOS(tha)
      
!	continue
      END FUNCTION sin_sun



       subroutine ydrov (li,lf,dt,wg,w2,wr,pl,veg,wrmax,xl,tcdil
     :   ,lev,letr,leg,ruigc,ruitc,vsw,cw1,cw2,wgeq,xdd2)
c
      implicit double precision (a-h,o-z)
c
      real*8 leg,lev,letr
      dimension wg(2),w2(2),wr(2)
c
c       evolution of wr
c
!      write(*,*) 'ydrov>>',li,lf,dt,wg(li),w2(li),wr(li)
!    :   ,pl,veg,wrmax,xl,tcdil
!    :   ,lev,letr,leg,ruigc,ruitc,vsw,cw1,cw2,wgeq,xdd2

       zev = lev / xl
       zetr = letr / xl


       zeg = leg / xl
       wr(lf) = wr(li) - dt * (zev - zetr - veg * pl )
       wr(lf) =dmax1(0.d0, wr(lf))
       ztemp= (wr(lf) - wrmax) / dt
       ruir =dmax1(0.d0, ztemp)
       wr(lf) = min(wr(lf), wrmax)
       zpg = (1. - veg) * pl + ruir
c
c         evolution of wg
c       - - - - - - - - - - - - - -
c
c      write(*,*) 'ydrov:',wg(li),tcdil,cw2,wgeq,dt,cw1,zeg,zpg
       wg(lf) = (wg(li) -
     &    dt * (cw1 * (zeg - zpg) / 1000. - cw2 * wgeq / tcdil))
     &    /(1.+dt*cw2/tcdil)
c
       ztemp=wg(lf)-vsw
       ruig =dmax1(0.d0, ztemp)
       ruigc = ruigc + ruig * 10.

c        evolutlon of w2
c       _______________
c
      ruip = ruig * 10./(xdd2*1000.)
      w2(lf) = w2(li) - dt*(zeg + zetr - zpg + ruip) / (xdd2 * 1000.)

c     ---------------------------------------------------
c                Run-off calculations
c     ----------------------------------------------------

      ztemp = w2(lf) - vsw
      ruitc = ruitc +dmax1(0.d0,ztemp) * xdd2 * 1000. + ruig * 10.
      wg(lf) = min(wg(lf), vsw)
      w2(lf) = min(w2(lf), vsw)

      return
      end






      subroutine dscfft2(f,wf,wx,wy,trigsx,ifax,sipcox,simcox,nx
     ;   ,trigsy,ifay,sipcoy,simcoy,ny,lot,isign)
c-----------------------------------------------------------------------
c double cosine transform on a staggered grid.
c
c for each is=1,ns performs the calculation of:
c
c    for i=1,nx;j=1,ny
c       sum(ix=1,nx;iy=1,ny) f(ix,iy,is)*cos((ix-0.5)*(i-1)*pi/nx)
c                            *cos((iy-0.5)*(j-1)*pi/ny)
c
c following wilhelmson and ericksen (j.c.p. 25,319-331) in the design
c of the pre and posprocessing (making a staggered cosine transform with
c a standard real transform).
c
c in the vectorization it has been assumed (in this routine, not in
c fft991) that nx>ny>ns.
c
c note the dimensions of f (they are appropriate for nh3d program)
c
c p. m. nov 1994
c-----------------------------------------------------------------------
c
      implicit real*8(a-h,o-z)
      dimension f(-1:nx+1,-1:ny+1,lot)
      dimension wf(-1:nx+1,-1:ny+1,lot)
c
      dimension wx(nx+2,ny,lot),wy(ny+2,nx,lot)
c
      dimension sipcox(nx),simcox(nx)
      dimension sipcoy(ny),simcoy(ny)
      dimension trigsx(*),trigsy(*)
      dimension ifax(10),ifay(10)
c
      if(isign.eq.1) then
c
c x transform:
c
         do iy=1,ny
            do j=1,lot
               wx(1,iy,j)=f(1,iy,j)
               wx(2,iy,j)=0.
               wx(nx+1,iy,j)=f(nx,iy,j)
               wx(nx+2,iy,j)=0.
            enddo
         enddo


         do iy=1,ny
            do j=1,lot
               do ii=2,nx-2,2
                  wx(ii+1,iy,j)=0.5*(f(ii,iy,j)+f(ii+1,iy,j))
                  wx(ii+2,iy,j)=-0.5*(f(ii,iy,j)-f(ii+1,iy,j))
               enddo
            enddo
         enddo
c
         call fft991(wx,wf,trigsx,ifax,1,nx+2,nx,lot*ny,+1)
c
         do iy=1,ny
            do j=1,lot
               f(1,iy,j)=wx(1,iy,j)
            enddo
         enddo


         do iy=1,ny
            do j=1,lot
               do ix=2,nx
                  f(ix,iy,j)=0.5*(sipcox(ix)*wx(ix,iy,j)-simcox(ix)
     :               *wx(nx-ix+2,iy,j))
               enddo
            enddo
         enddo
c
c y transform:
c
         do j=1,lot
            do ix=1,nx
               wy(1,ix,j)=f(ix,1,j)
               wy(2,ix,j)=0.
               wy(ny+1,ix,j)=f(ix,ny,j)
               wy(ny+2,ix,j)=0.
            enddo
         enddo


         do ii=2,ny-2,2
            do j=1,lot
               do ix=1,nx
                  wy(ii+1,ix,j)=0.5*(f(ix,ii,j)+f(ix,ii+1,j))
                  wy(ii+2,ix,j)=-0.5*(f(ix,ii,j)-f(ix,ii+1,j))
               enddo
            enddo
         enddo
c
         call fft991(wy,wf,trigsy,ifay,1,ny+2,ny,lot*nx,+1)
c
         do j=1,lot
            do ix=1,nx
               f(ix,1,j)=wy(1,ix,j)
            enddo
         enddo


         do iy=2,ny
            do j=1,lot
               do ix=1,nx
                  f(ix,iy,j)=0.5*(sipcoy(iy)*wy(iy,ix,j)-simcoy(iy)
     :               *wy(ny-iy+2,ix,j))
               enddo
            enddo
         enddo
      else
c
c inverse x transform:
c
         do iy=1,ny
            do j=1,lot
               wx(1,iy,j)=f(1,iy,j)
            enddo
         enddo


         do j=1,lot
            do iy=1,ny
               do ix=2,nx
                  wx(ix,iy,j)=sipcox(ix)*f(ix,iy,j)
     :               +simcox(ix)*f(nx-ix+2,iy,j)
               enddo
            enddo
         enddo
c
         call fft991(wx,wf,trigsx,ifax,1,nx+2,nx,lot*ny,-1)
c
         do iy=1,ny
            do j=1,lot
               f(1,iy,j)=wx(1,iy,j)
               f(nx,iy,j)=wx(nx+1,iy,j)
            enddo
         enddo


         do iy=1,ny
            do j=1,lot
               do ii=2,nx-2,2
                  f(ii,iy,j)=wx(ii+1,iy,j)-wx(ii+2,iy,j)
                  f(ii+1,iy,j)=wx(ii+1,iy,j)+wx(ii+2,iy,j)
               enddo
            enddo
         enddo


c
c inverse y transform:
c
         do j=1,lot
            do ix=1,nx
               wy(1,ix,j)=f(ix,1,j)
            enddo
         enddo


         do iy=2,ny
            do j=1,lot
               do ix=1,nx
                  wy(iy,ix,j)=sipcoy(iy)*f(ix,iy,j)
     :               +simcoy(iy)*f(ix,ny-iy+2,j)
               enddo
            enddo
         enddo
c
         call fft991(wy,wf,trigsy,ifay,1,ny+2,ny,lot*nx,-1)
c
         do j=1,lot
            do ix=1,nx
               f(ix,1,j)=wy(1,ix,j)
               f(ix,ny,j)=wy(ny+1,ix,j)
            enddo
         enddo


         do ii=2,ny-2,2
            do j=1,lot
               do ix=1,nx
                  f(ix,ii,j)=wy(ii+1,ix,j)-wy(ii+2,ix,j)
                  f(ix,ii+1,j)=wy(ii+1,ix,j)+wy(ii+2,ix,j)
               enddo
            enddo
         enddo
      endif
      return
      end




      subroutine inisoil(tsmed,qvsoil,nx,ny,fnmap)
c
c initializes soil variables for a simple "soil"
c tsmed: parameter for the analytical definition of the soil
c        temperature (cf THERMO)
c qvsoil: constant near soil air specific humidity
c
      implicit real*8(a-h,o-z)
      dimension qvsoil(0:nx+1,0:ny+1)
      dimension tsmed(0:nx+1,0:ny+1)
      character*80 fnmap,fext


      open(1,file=fext(fnmap,'tsm'),status='old')
      call readgrd(1,tsmed,0,nx+1,0,ny+1)
      close(1)




      open(1,file=fext(fnmap,'qso'),status='old')
      call readgrd(1,qvsoil,0,nx+1,0,ny+1)
      close(1)


      return
      end




      subroutine ireadgrd(nunit,ivar,i0,i1,j0,j1)
c
c reads 2d array from a grd file
c
      implicit real*8(a-h,o-z)
      dimension ivar(i0:i1,j0:j1)
      character*4 dsaa


      read(nunit,'(a4)') dsaa
      read(nunit,*) nx,ny
      read(nunit,*) xmin,xmax
      read(nunit,*) ymin,ymax
      read(nunit,*) zmin,zmax


      do j=j0,j1
         read(nunit,*) (ivar(i,j),i=i0,i1)
      enddo


      return
      end




      subroutine outarr (a,l0,l1,i0,i1,is,ixp,m0,m1,j0,j1,js,iyp
     :   ,n0,n1,k0,k1,ks,isp,nplan,leng,title,nchan,scale,nexp
     :   ,factor,ndigit,tkoff,mode,nstep)
c
c-----------------------------------------------------------------------
c  write a section of array a(l,m,n) into channel nchan
c  data to be printed out :                                             s)
c  ((a(i,j,k),i=i0,i1,is),j=j0,j1,js),k=k0,k1,ks)
c
c  leng :   maximum line length required (number of characters)
c
c  this geometry is appropriate for sigma vertical coordinate. in spite
c  of the arrangement of the j (y) coordinate, in all cases i0,j0 and k0
c  are always displayed, but not necessarily i1,j1,k1, depending on
c  is,js,ks.
c
c  p.m. nov 1988
c-----------------------------------------------------------------------
c
      implicit real*8(a-h,o-z)
      integer lc(-200:200)
      dimension a(l0:l1,m0:m1,n0:*)
      character*8 title
      character*9 cmod(3)
      character*2 numero(39)
      character*20 forma1
      character*10 forma2
      character*13 forma3
      character*80 filgrd
      character*80 fnout
      logical dogrd
      common/cfnout/fnout
c
      dimension ixp(nplan),iyp(nplan),isp(nplan)
c
      save cmod,forma1,numero,forma2,forma3
c
      data cmod/'plane is=','plane iy=','plane ix='/
      data numero/'01','02','03','04','05','06','07','08','09'
     :   ,'10','11','12','13','14','15','16','17','18','19'
     :   ,'20','21','22','23','24','25','26','27','28','29'
     :   ,'30','31','32','33','34','35','36','37','38','39'/
c
      data forma1/'(1x,i2,''*'',t6,31i04)'/
      data forma2/'(t6,31i04)'/
      data forma3/'(i4,04e30.23)'/


      dogrd=.false.
      dx=4000.
      dy=4000.
c
c      if(ndigit.lt.1 .or. ndigit.gt.15) then
c         write(nchan,*) 'error in outarr: ndigit=',ndigit
c         return
c      endif
c
c      write(nchan,*) l0,l1,i0,i1,is,ixp
c      write(nchan,*) m0,m1,j0,j1,js,iyp
c      write(nchan,*) n0,n1,k0,k1,ks,isp
c      write(nchan,*) nplan,leng,title,nchan,scale,nexp
c      write(nchan,*) factor,ndigit,tkoff,mode
c
c
c maximum number of elements per line
c
         lengl=leng-6
         lfield=ndigit+2
         line=lengl/lfield
      if(.not.dogrd) then
         lengl=leng-6
         lfield=ndigit+2
         line=lengl/lfield
c
         forma1(18:19)=numero(ndigit+2)
         forma1(15:16)=numero(line)
         forma2(8:9)=numero(ndigit+2)
         forma2(5:6)=numero(line)
c
         ndec=max(1,ndigit+2-7)
         ndig=ndec+7
         forma3(8:9)=numero(ndig)
         forma3(11:12)=numero(ndec)
         forma3(5:6)=numero(line)
      endif
c
c      write(nchan,*) 'forma1:',forma1
c      write(nchan,*) 'forma2:',forma2
c      write(nchan,*) 'forma3:',forma3
c
c      write(nchan,2) nexp
c
c write horizontal plane (k=cst):
c
      if(mode.eq.0 .or. mode.eq.1) then
         i11=i0+min(i1-i0,(line-1)*is)
         imax=min(lengl,(lfield*((i1-i0)/is+1)))
         li=i0+(line-1)*is
         i12=min(i1,li)
         do 210 kp=1,nplan
         if(isp(kp).ge.n0) then
            k=isp(kp)
            if(dogrd) then
               write(filgrd,'(a40,''.'',i1,i2,i5,''.'',a2)') fnout
     :            ,1,k,nstep,title
               call packstr(filgrd,60)


               open(95,file=filgrd,status='unknown')
               write(95,'(a4)') 'DSAA'
               write(95,*) i1-i0+1,j1-j0+1
               xmin=0.
               xmax=(i1-i0)*dx/1000.
               ymin=0.
               ymax=(j1-j0)*dy/1000.
               write(95,*) xmin,xmax
               write(95,*) ymin,ymax
               do j=j0,j1
                  do i=i0,i1
                     amin=min(amin,a(i,j,k))
                     amax=max(amax,a(i,j,k))
                  enddo
               enddo
               write(95,*) amin,amax
               do 221 j=j0,j1
                  write(95,'(10e13.5)') (a(i,j,k),i=i0,i1)
221            continue
               close(95)
            else
               write(nchan,10) title,tkoff,cmod(1),k,nexp
               write(nchan,forma2) (i,i=i0,i11,is)
               write(nchan,4) ('*',i=2,imax)
               do 200 j=j1,j0,-js
               if(ndigit.lt.15) then
                  do 220 i=i0,i1,is
                  lc(i)=nint((a(i,j,k)-tkoff)*factor)
220               continue
                  write(nchan,forma1) j,(lc(i),i=i0,i12,is)
               else
                  write(nchan,forma3) j,(a(i,j,k),i=i0,i12,is)
               endif
200            continue
               write(nchan,4) ('*',i=2,imax)
            endif
         endif
210      continue
      endif
      if(mode.eq.0 .or. mode.eq.2) then
c
c write vertical plane (j=cst):
c
         i11=i0+min(i1-i0,(line-1)*is)
         imax=min(lengl,(lfield*((i1-i0)/is+1)))
         li=i0+(line-1)*is
         i12=min(i1,li)
         do 310 jp=1,nplan
         if(iyp(jp).ge.m0) then
            j=iyp(jp)
            if(dogrd) then
               write(filgrd,'(a40,''.'',i1,i2,i5,''.'',a2)') fnout
     :            ,2,j,nstep,title
               call packstr(filgrd,60)


               open(95,file=filgrd,status='unknown')
               write(95,'(a4)') 'DSAA'
               write(95,*) i1-i0+1,k1-k0+1
               xmin=0.
               xmax=(i1-i0)*dx/1000.
               write(95,*) xmin,xmax
               write(95,*) 0.,1.
               do j=j0,j1
                  do i=i0,i1
                     amin=min(amin,a(i,j,k))
                     amax=max(amax,a(i,j,k))
                  enddo
               enddo
               write(95,*) amin,amax
               do 321 k=k0,k1
                  write(95,'(10e13.5)') (a(i,j,k),i=i0,i1)
321            continue
               close(95)
            else
               write(nchan,10) title,tkoff,cmod(2),j,nexp
               write(nchan,forma2) (i,i=i0,i11,is)
               write(nchan,4) ('*',i=2,imax)
               do 300 k=k0,k1,ks
               if(ndigit.lt.15) then
                  do 320 i=i0,i1,is
                  lc(i)=nint((a(i,j,k)-tkoff)*factor)
320               continue
                  write(nchan,forma1) k,(lc(i),i=i0,i12,is)
               else
                  write(nchan,forma3) k,(a(i,j,k),i=i0,i12,is)
               endif
300            continue
               write(nchan,4) ('*',i=2,imax)
           endif
         endif
310      continue
      endif
      if(mode.eq.0 .or. mode.eq.3) then
c
c write vertical plane (i=cst):
c
         j11=j0+min(j1-j0,(line-1)*js)
         jmax=min(lengl,(lfield*((j1-j0)/js+1)))
         mi=j0+(line-1)*js
         j12=min(j1,mi)
         do 410 ip=1,nplan
         if(ixp(ip).ge.l0) then
            i=ixp(ip)
            if(dogrd) then
               write(filgrd,'(a40,''.'',i1,i2,i5,''.'',a2)') fnout
     :            ,3,i,nstep,title
               call packstr(filgrd,60)
               open(95,file=filgrd,status='unknown')
               write(95,'(a4)') 'DSAA'
               write(95,*) j1-j0+1,k1-k0+1
               ymin=0.
               ymax=(j1-j0)*dy/1000.
               write(95,*) ymin,ymax
               write(95,*) 0.,1.
               do j=j0,j1
                  do i=i0,i1
                     amin=min(amin,a(i,j,k))
                     amax=max(amax,a(i,j,k))
                  enddo
               enddo
               write(95,*) amin,amax
               do 421 k=k0,k1
                  write(95,'(10e13.5)') (a(i,j,k),j=j0,j1)
421            continue
               close(95)
            else
               write(nchan,10) title,tkoff,cmod(3),i,nexp
               write(nchan,forma2) (j,j=j11,j0,-js)
               write(nchan,4) ('*',j=2,jmax)
               do 400 k=k0,k1,ks
               if(ndigit.lt.15) then
                  do 420 j=j1,j0,-js
                  lc(j)=nint((a(i,j,k)-tkoff)*factor)
420               continue
                  write(nchan,forma1) k,(lc(j),j=j12,j0,-js)
               else
                  write(nchan,forma3) k,(a(i,j,k),j=j12,j0,-js)
               endif
400            continue
               write(nchan,4) ('*',j=2,jmax)
            endif
         endif
410      continue
      endif
c
      return
2     format(t30,'units of 10**',i3)
3     format(t6,31i4)
4     format(t7,125a1)
5     format(1x,i2,'*',t6,31i4)
10    format(//t7,a8,'(add:',e10.3,')',3x,a9,i4,5x,'units of 10**'
     :   ,i3,/)
      end




      subroutine packstr(string,n)
      character (len=*) string
      character*128 string0


      string0=string
      string=' '
      k=0
      do i=1,max(128,len(string))
        if(string0(i:i).ne.' ') then
           k=k+1
           string(k:k)=string0(i:i)
        endif
      enddo
      return
      end




      subroutine pois2n(phisurf,phsfco,i0,j0,wk1,wk2,dxx,dyy,eigxy
     :   ,trigsx,ifax,sipcox,simcox,nx
     :   ,trigsy,ifay,sipcoy,simcoy,ny)
c
c This routine solves a 2-d Poisson equation with
c homogeneous Neuman boundary conditions
c using a 2-d Staggered Cosine Fourier Transform
c The problem with the undetermined constant is solved
c by imposing the value of the first coefficient in the double
c Fourier series, which is related to the mean of the
c series.
c
c Note 2<=i0<=nx
c Note 2<=j0<=ny
c
c Note the dimensions of phisurf. They are appropriate for
c the nh3d program (for consistency with the 3d Poisson solver)
c
      implicit real*8(a-h,o-z)


      dimension phisurf(0:nx+1,0:ny+1)
     :   ,wk1(0:nx+1,0:ny+1),wk2(0:nx+1,0:ny+1)
     :   ,eigxy(2:nx,2:ny)
      dimension trigsx(*),trigsy(*)
      dimension sipcox(*),simcox(*),ifax(*)
      dimension sipcoy(*),simcoy(*),ifay(*)
      logical fstcall
      save fstcall


      data fstcall/.true./
c


      if(fstcall) then
         fstcall=.false.
         pi=4.*datan(1.d0)
         iy=2
         eigxy(2,2)=0.
         do ix=3,nx
            eigxy(ix,iy)=1.d0/((2./dxx*dcos((ix-2)*pi/(nx-1)))
     :         -(2./dxx+2./dyy)+(2./dyy*dcos((iy-2)*pi/(ny-1))))
         enddo
         do iy=3,ny
            do ix=2,nx
               eigxy(ix,iy)=1.d0/((2./dxx*dcos((ix-2)*pi/(nx-1)))
     :            -(2./dxx+2./dyy)+(2./dyy*dcos((iy-2)*pi/(ny-1))))
            enddo
         enddo
      endif


      call dscfft2(phisurf,wk1,wk2,wk2,trigsx,ifax,sipcox,simcox,nx-1
     :   ,trigsy,ifay,sipcoy,simcoy,ny-1,1,+1)


      do iy=2,ny
         do ix=2,nx
            phisurf(ix,iy)=phisurf(ix,iy)*eigxy(ix,iy)
         enddo
      enddo


      call dscfft2(phisurf,wk1,wk2,wk2,trigsx,ifax,sipcox,simcox,nx-1
     :   ,trigsy,ifay,sipcoy,simcoy,ny-1,1,-1)


c calculates additive constant:
c
      phicon=phisurf(i0,j0)-phsfco
c
      do iy=2,ny
         do ix=2,nx
            phisurf(ix,iy)=phisurf(ix,iy)-phicon
         enddo
      enddo
c
c fulfills boundary conditions:
c
      call extrah(nx+1,ny+1,phisurf,1,nx+1,1,ny+1,1,1)
c
      return
      end




      subroutine readgrd(nunit,var,i0,i1,j0,j1)
c
c reads 2d array from a grd file
c
      implicit real*8(a-h,o-z)
      dimension var(i0:i1,j0:j1)
      character*4 dsaa


      read(nunit,'(a4)') dsaa
      read(nunit,*) nx,ny
      read(nunit,*) xmin,xmax
      read(nunit,*) ymin,ymax
      read(nunit,*) zmin,zmax


      do j=j0,j1
         read(nunit,*) (var(i,j),i=i0,i1)
      enddo


      return
      end




      subroutine toupcase(string)
c
c converts string to uppercase
c
      character (len=*) string
      integer i,ismall,ibig
      ismall=ichar('a')
      ibig=ichar('A')
      do i=1,len(string)
        if(string(i:i).ge.'a' .and. string(i:i).le.'z') then
          string(i:i)=char(ichar(string(i:i))+ibig-ismall)
        endif
      enddo
      return
      end




      subroutine sigrid(sigma0,sigma1,ns)
c
c defines the vertical grid.
c
      implicit real*8(a-h,o-z)
      dimension sigma0(0:ns+1),sigma1(0:ns+1)
c
c equally spaced grid:
c
      dsp=1./(ns-1.)
      do is=0,ns+1
        sigma1(is)=(float(is)-1.0)*dsp
        sigma0(is)=(float(is)-0.5)*dsp
      enddo
      return
      end




      subroutine solonoi (iclay,isand,iveg,veg,xlai,rsmin,alb,emis,z0
     :   ,z0h,rra,ta,qa,ua,rg,rat,ps,pl,zref
     &   ,ts,t2,wg,w2,wr,rn,h,le,g,dt,li,lf,xdd2
     &   ,cdm,rhoa,tsfunc)
c
c---------------------------------------------------------------------
c
c  This routine contains the soil model of Noilhan e Planton (NP89)
c  it has been modified from code given by
c  Joel Noilhan (October 1994)
c
c iclay: percentage of clay in the soil
c isand: percentage of sand in the soil
c Last changed: December 1994
c
c
c
c---------------------------------------------------------------------
c
      implicit double precision (a-h,o-z)
      real*8 le,leg,lev,letr
      dimension ts(2),t2(2),wg(2),w2(2),wr(2)
      dimension csowi(-2:100),weai(-2:100),vswi(-2:100),bwi(-2:100)
     :   ,webi(-2:100),vfli(-2:100),vcci(-2:100),cw2ai(-2:100)
     :   ,cw1si(-2:100)
c    :   ,pswi(),coswi()
      dimension rsai(2),alfai(2),cvv(2)


c
c physical constants
c
      parameter(gg = 9.80665)
      parameter(ra = 287.05)
      parameter(rv = 461.51)
      parameter(cp = 1005.46)
      parameter(cpv = 1869.46)
      parameter(xl = 2500800.)
      parameter(stefan = 5.6697e-08)
      parameter(xpi = 3.141592654)
      parameter(too = 273.16)
      parameter(xpoo = 100000.)
      parameter(vkarmn = .4)
      parameter(xli = 2834500.)
      parameter(tcdil = 86400.)
      parameter(rascp = ra / cp)
      parameter(xlsg = xl / gg)
      parameter(xlscp = xl / cp)
      parameter(cpsl = cp / xl)
      parameter(cpsg = cp / gg)
      parameter(gscp = gg / cp)
      parameter(unsg = 1. / gg)
      parameter(gsra = gg / ra)
      parameter(rasl = ra / xl)
      parameter(rascp2 = ra / (2.*cp))
      parameter(rasrv = ra / rv)
      parameter(rvsra = rv / ra)
      parameter(etv = rvsra - 1.)
      parameter(ecph = cpv / cp - 1.)
      parameter(xlsrv = xl / rv)
      parameter(unscp = 1. / cp)
      parameter(cpvmcp = cpv - cp)
      parameter(etvq = 1. - rasrv)
      parameter(xlf = xli - xl)
      parameter(xliscp = xli / cp)
      parameter(xlisg = xli / gg)


c
c---------------------------------------------------------------------
c soil properties following Clapp e Hornberg, 1978, using parametrized
c formulas from Noilhan, as a function of the percentage of clay and
c sand
c
c note iclay.eq.-1 for imposed soil temperature
c note iclay.eq.-2 for water (lake)
c
c Note: initialization of the arrays with soils properties is done
c only in the first call to solonoi (when ifirst.eq.0)
c


      save ifirst
      save csowi,vswi,bwi,weai,webi,cw2ai,vcci,vfli,cw1si
      data ifirst/0/
c
cch   Clapp and Hornberg formulation
c
cch    data csowi/3.,3.222,3.057,3.56,4.418,4.111,3.67,3.593,3.995
cch   :   ,3.058,3.729,3.6/
cch    data pswi/.1,.121,.09,.218,.786,.478,.299,.356,.63,.153,.49,.405/
cch    data vswi/.4,.395,.41,.435,.485,.451,.42,.477,.476,.426,.492,.482/
cch    data bwi /4.,4.05,4.38,4.9,5.3,5.39,7.12,7.75,8.52,10.4,10.4,11.4/
cch    data weai /.4,.387,.404,.219,.105,.148,.135,.127,.084,.0139,.075,
cch   &          .083/
cch    data webi /4.,4.,4.,4.,6.,6.,6.,8.,10.,8.,10.,12./
cch    data cw2ai /4.,3.9,3.7,1.8,.8,.8,.8,.4,.6,.3,.3,.3/
cch    data vcci /.14,.135,.15,.195,.255,.24,.255,.322,.324,.31,.37,.367/
cch    data vfli/.1,.068,.075,.114,.179,.155,.175,.218,.25,.219,.283
cch   :   ,286./
cch    data coswi /176,176.,156.3,3.41,7.2,7.,6.3,1.7,2.5,2.2,1.,1.3/
c
c
c
c**********************************************************************
c properties of vegetation(Mattthews 1983)
c 2 classes
c
      data rsai,alfai,cvv /100.,30.,0.,.04,1.e-6,2.e-5/
c---------------------------------------------------------------------
c     for comparison with Giordani
c
c     data rsai,alfai,cvv /100.,30.,0.,.04,1.2e-5,2.e-5/
c---------------------------------------------------------------------
      if (ifirst.eq.0) then
         ifirst=1
         do icl=1,100
            bwi(icl)=0.137*icl+3.501
            weai(icl)=0.73242*icl**(-0.539)
            webi(icl)=0.134*icl+3.4
            cw1si(icl)=.558e-2*icl+.08488
            cw2ai(icl)=13.815*icl**(-0.954)
            csowi(icl)=(-0.008*icl+3.96)
            vswi(icl)=.494305-1.08e-3*icl
            vcci(icl)=.0890467*icl**.3496
            vfli(icl)=.0371342*(icl)**.5
         enddo
         bwi(0)=1.
         weai(0)=1.
         webi(0)=1.
         cw1si(0)=1.
         cw2ai(0)=1.
         csowi(0)=1.
         vswi(0)=1.
         vcci(0)=1.
         vfli(0)=1.
         bwi(-1)=1.
         weai(-1)=1.
         webi(-1)=1.
         cw1si(-1)=1.
         cw2ai(-1)=1.
         csowi(-1)=1.
         vswi(-1)=1.
         vcci(-1)=1.
         vfli(-1)=1.
         bwi(-2)=1.
         weai(-2)=1.
         webi(-2)=1.
         cw1si(-2)=1.
         cw2ai(-2)=1.
         csowi(-2)=1.
         vswi(-2)=1.
         vcci(-2)=1.
         vfli(-2)=1.
      endif




c
c soil properties
c
      csow=csowi(isand)
cch   psw= pswi(iclay)
      vsw= vswi(isand)
      bw = bwi(iclay)
      wea= weai(iclay)
      web= webi(iclay)
      cw2a=cw2ai(iclay)
      vcc= vcci(iclay)
      vfl= vfli(iclay)
cch   cosw=coswi(iclay)*.000001
      cw1s=cw1si(iclay)
c
c   wg inicialization (if wg(i,j)=0, in the begining)
c
c     if (wg(li).lt.1.e-5.and.iclay.ge.0) then
c        wg(li)=(dacos(1.-2.*qa/qsat(ts(li),ps)))*vcc/xpi
c     elseif(iclay.lt.0) then
c        wg(li)=1.
c        wg(lf)=1.
c     endif
      cv=cvv(iveg)
      wrmax = .2*veg*xlai
      rsa=rsai(iveg)
      alfa=alfai(iveg)
c
c**********************************************************************
c          c lculo das propriedades do solo (cso,cw1,cw2)
c               (corresponde … rotina ysol de Noilhan)
c
c     FOR EFEDA TEST (SAME CONDITIONS THAN GIORDANI ET AL.)
c     xlambda=.37
c     por=.477
c     wk=.05
c     ws=.59
      alfa2=.66
      gsrv=gg/rv
      valnon=-999.
c
c--------------------------------------------------------
c     for comparison with Giordani
c
c     cg=(1.4+wg(li)*4.18)*1.e6
c     cso=2.*sqrt(xpi/(xlambda*cg*tcdil))
c--------------------------------------------------------
c    if you dont know the thermal properties of the soil
c    then you have to use the following expressions
c-------------------------------------------------------
      cso = csow*((vsw/w2(li))**(bw/(2.*log(10.))))
     &      *1.e-6
c-------------------------------------------------------
      if (iclay.ge.0) then
cch      cw1s=2.*sqrt(xpi*vsw)/(bw*psw*cosw*tcdil)
         cw1=cw1s*((vsw/wg(li))**(.5*bw+1.))
         cw1=dmin1(45.d0,cw1)
         cw2 = cw2a * w2(li) / (vsw - w2(li) + .01)
         zx = w2(li) / vsw
         wgeq = zx - wea * zx**web* (1. - zx**(8. * web))
         wgeq = wgeq * vsw
      endif
c
c**********************************************************************
c
c        surface resistance calculation if veg>0
c
       if(veg.gt.0) then
         zf = .55 * rg / rsa * 2. / xlai
         zf1 = (1. + zf) / (zf + rsmin / 5000.)
         zf3 = 1. - alfa * (qsat(ta, ps) - qa) * 1000.
         zf4 = 1. - .0016 *((298. - ta)**2)
         zf4=max(1.d-03,zf4)
         if (w2(li).ge.vcc) then
             zf2=1.
         elseif (w2(li) .le. vfl) then
            zf2 =1.e-6
         else
            zf2 = (w2(li) - vfl) / (vcc - vfl)
         endif
         rs = rsmin / xlai * zf1 / zf2 / zf3 / zf4
         rs =dmin1(rs, 5000.d0)
c
         if (rs .lt. 0) then
            write(*,270)
            write(*,*) veg,iclay,tsfunc
            write (*,*) zf1,zf2,zf3,zf4,qa,ta,ps,qsat(ta,ps)
            write (*,*) 'rg,rsa,xlai,zf,rsmin',rg,rsa,xlai,zf,rsmin
270      format('negative surface resistance, progr stops')
            write(0,*) 'nh3dFatalError=ISBA negative surface resistance'
            stop
         endif
      else
         rs=0.
      endif
c
c**********************************************************************
c
      call ydifus (li,lf,dt,wg,wr,ts,ps,t2,qa,ta,rg,rat,ua
     &       ,veg,vcc,wrmax,cso,cv,alb,emis,rs,rra,z0,z0h
     &       ,le,leg,lev,letr,g,rn,h,hu
     &       ,ra,etv,stefan,xl,cp,tcdil,xpi,zref,ct,gscp,gg,vkarmn  !zref ->deltazt
     &       ,cdm,rhoa,iclay,tsfunc)


      if(iclay.ge.0) then
         call ydrov (li,lf,dt,wg,w2,wr,pl,veg,wrmax,xl,tcdil,
     &      lev,letr,leg,ruigc,ruitc,vsw,cw1,cw2,wgeq,xdd2)
      else
         wg(lf)=wg(li)
         w2(lf)=w2(li)
         wr(lf)=wr(li)
      endif
c
c
c next seems irrelevant (and would need initialization)
c
c     zplc=zplc+pl*dt
c     zle=zle+le*dt/xl
c     zlev=zlev+lev*dt/xl
c     zletr=zletr+letr*dt/xl
      return
      end




      subroutine uphsuf(hsuf,hsufmax,hnew,hmount,nx,ny,nstep,nstepgrow)
c
c updates hsuf : case of varying mountain height (no change in shape)
c
      implicit real*8(a-h,o-z)
c
      dimension hsuf(0:nx+1,0:ny+1)
      dimension hsufmax(0:nx+1,0:ny+1)
c
      hnew=hmount*min(1.d0,dble(nstep)/dble(nstepgrow))


      do iy=0,ny+1
         do ix=0,nx+1
            hsuf(ix,iy)=hsufmax(ix,iy)*hnew/hmount
         enddo
      enddo


      return
      end
      subroutine uphsuf2(hsuf,hsufmax,hnew,hmount,nx,ny,nstep,nstepgrow
     :  ,nstepflat,numsteps)
c
c updates hsuf : case of varying mountain height (no change in shape)
c
      implicit real*8(a-h,o-z)
c
      dimension hsuf(0:nx+1,0:ny+1)
      dimension hsufmax(0:nx+1,0:ny+1)
c
      istep=nstep/(nstepgrow+nstepflat)+1
      ifrac=mod(nstep,nstepgrow+nstepflat)
      hfrac=(float(istep-1)+max(1.d0,dble(ifrac)/nstepgrow))
     :  /float(numsteps)
      hnew=hmount*hfrac


      do iy=0,ny+1
         do ix=0,nx+1
            hsuf(ix,iy)=hsufmax(ix,iy)*hnew/hmount
         enddo
      enddo


      return
      end




      subroutine uprof(iorefu,usdat,zusdat,ndus,fpro,zpro,npro)
c
c creates reference state wind profile
c iorefu=1 then "interpolates data profile (usdat,vsdat)"
c else : analytical profile
c
      implicit real*8(a-h,o-z)
      dimension usdat(0:*),zusdat(0:*),fpro(0:*),zpro(0:*)
c
      if(iorefu.eq.0) then
         call inipro(usdat,zusdat,ndus,fpro(0),zpro,npro)
      elseif(iorefu.eq.1) then
         pi=4.*atan(1.)
         u0=10.
         bv=0.01
         hydlam=2.*pi*u0/bv
         un=u0/bv
         h0=3./4.*hydlam
         do 10 ipro=0,npro+1
         fpro(ipro)=-u0*tanh((zpro(ipro)-h0)/un)
10       continue
      elseif(iorefu.eq.2) then
         pi=4.*atan(1.)
         u0=10.
         bv=0.01
         hydlam=2.*pi*u0/bv
         un=u0/bv
         h0=1./4.*hydlam
         do 20 ipro=0,npro+1
         fpro(ipro)=-u0*tanh((zpro(ipro)-h0)/un)
20       continue
      elseif(iorefu.eq.3) then
         pi=4.*atan(1.)
         u0=10.
         bv=0.01
         hydlam=2.*pi*u0/bv
         un=u0/bv
         h0=2.*un
         do 30 ipro=0,npro+1
         fpro(ipro)=-u0*tanh((zpro(ipro)-h0)/un)
30       continue
      elseif(iorefu.eq.4) then
c
c for wind rotating with height (constant speed) choose iorefu=4, iorefv=5 (note uprof
c is called twice, first for u then for v)
c usdat(0) is the speed
c usdat(1) is the wind direction at z=0 (not the meteorological ddd)
c usdat(2) is the rate of change
c
         pi=4.*atan(1.)
         u0=usdat(0)
         alpha0=usdat(1)*pi/180.
         beta=usdat(2)
         do ipro=0,npro+1
           fpro(ipro)=u0*cos(alpha0+beta*zpro(ipro))
         enddo
      elseif(iorefu.eq.5) then
         pi=4.*atan(1.)
         u0=usdat(0)
         alpha0=usdat(1)*pi/180.
         beta=usdat(2)
         do ipro=0,npro+1
           fpro(ipro)=u0*sin(alpha0+beta*zpro(ipro))
         enddo
      elseif(iorefu.eq.6) then
         u0=usdat(0)
         alpha=usdat(1)
         hd=usdat(2)
         h=usdat(3)
         do ipro=0,npro+1
           z=zpro(ipro)
           if(zpro(ipro).le.hd) then
             fpro(ipro)=u0+alpha*z
           else
             fpro(ipro)=u0*(2*hd-h-z)*(z-h)/(hd-h)**2+alpha*(z-h)
     :         *(hd**2-h*z)/(h-hd)**2
           endif
         enddo
      elseif(iorefu.eq.7) then
         u0=usdat(0)
         alpha=usdat(1)
         hd=usdat(2)
         h=usdat(3)
         do ipro=0,npro+1
           z=zpro(ipro)
           if(zpro(ipro).le.hd) then
             fpro(ipro)=u0
           else
             fpro(ipro)=u0*(2*hd-h-z)*(z-h)/(hd-h)**2
           endif
         enddo
      elseif(iorefu.eq.8) then
         u0=usdat(0)
         hcrit=usdat(1)
         do ipro=0,npro+1
           z=zpro(ipro)
           fpro(ipro)=max(-30.d0,u0*(1.d0-(z/hcrit)**2))
         enddo
      elseif(iorefu.eq.9) then
         u0=usdat(0)
         dudz=usdat(1)
         umin=usdat(2)
         umax=usdat(3)
         do ipro=0,npro+1
           z=zpro(ipro)
           fpro(ipro)=max(umin,min(umax,u0+dudz*z))
         enddo
      elseif(iorefu.eq.10) then
c usdat(0)=u0
c usdat(1)=dudz0
c u=(a+bz)**(2/3)


         a=usdat(0)**(3./2.)
         b=3./2.*usdat(1)*usdat(0)**(0.5)
         umin=usdat(2)
         umax=usdat(3)
         do ipro=0,npro+1
           z=zpro(ipro)
           if(z.lt.-a/b) then
             fpro(ipro)=max(umin,min(umax,(a+b*z)**(2./3.)))
           else
             fpro(ipro)=0.
           endif
         enddo
      elseif(iorefu.eq.11) then
c piecewise linear profile u=u0,z<h1; u=u0-dudz*z, h1<z<h2; u=const, z>h2
c usdat(0)=u0
c usdat(1)=dudz
c usdat(2)=h1
c usdat(3)=h2
         u0=usdat(0)
         dudz=usdat(1)
         h1=usdat(2)
         h2=usdat(3)
         do ipro=0,npro+1
           z=zpro(ipro)
           if(z.le.h1) then
             fpro(ipro)=u0
           elseif(z.gt.h2) then
             fpro(ipro)=fpro(ipro-1)
           else
             fpro(ipro)=u0+dudz*(z-h1)
           endif
         enddo
      else
         write(0,*) 'nh3dFatalError=unknown u profile:',iorefu
         stop
      endif
      return
      end






      subroutine wri2ar(a,i0,i1,j0,j1,title,tkoff,nx,ny)
c-----------------------------------------------------------------------
c handler for wrigar (case of standard 2d arrays of nh3d)
c-----------------------------------------------------------------------
c
      implicit real*8(a-h,o-z)
c
      character*8 title
      dimension a(0:nx+1,0:ny+1)
c
      call wrigar(a,0,nx+1,0,ny+1,0,0,i0,i1,j0,j1,0,0,title,tkoff,1)
c
      return
      end




      subroutine wri3ar(a,i0,i1,j0,j1,k0,k1,title,tkoff,mode,nx,ny)
c-----------------------------------------------------------------------
c handler for wrigar (case of standard 3d arrays of nh3d)
c-----------------------------------------------------------------------
c
      implicit real*8(a-h,o-z)
c
      character*8 title
      dimension a(0:nx+1,0:ny+1,0:*)
c
      call wrigar(a,0,nx+1,0,ny+1,0,ns1,i0,i1,j0,j1,k0,k1
     :   ,title,tkoff,mode)
c
      return
      end




      subroutine ydifus (li,lf,dt,wg,wr,ts,ps,t2,qa,ta,rg,rat,ua
     &       ,veg,vcc,wrmax,cso,cv,alb,emis,rs,rra,z0,z0h
     &       ,le,leg,lev,letr,g,rn,h,hu
     &       ,ra,etv,stefan,xl,cp,tcdil,xpi,zref,ct,gscp,gg,vkarmn
     &       ,cdm,rhoa,iclay,tsfunc)
c
c Note tsfunc is the new value for ts, if iclay<0. It is imposed
c and should be defined outside, prior to this call. This allows for
c the imposition of a particular cycle for soil or lake temperature
c as a function of time.
c
      implicit double precision (a-h,o-z)
      real*8 le,leg,lev,letr
      parameter (ns=2)
      dimension ts(ns), t2(ns), wg(ns), wr(ns)
c
c              hu (surface humidity) calculation
c
      xhu=1.


      if (iclay.ge.-1.and.wg(li).le.vcc) then
         ztemp=xpi*wg(li)/vcc
         hu=.5*(1.-cos(ztemp))
c
c       specific treatment for dew (see mahfouf-noilhan,jam,91,tabl2)
c
        lmn91=1
	
        if (lmn91.eq.1) then
          if (hu*qsat(ts(li),ps).lt.qa) then
             if (qsat(ts(li),ps).gt.qa) then
               hu=qa/qsat(ts(li),ps)
               xhu=1.
             else
               hu=1.
             endif
           endif
         endif
      else
         hu=1.
      endif

c
c**********************************************************************
c           delta calculatlon
c
      delta=0.
      if (veg.gt.0.) delta=(wr(li)/wrmax)**(2./3.)
c
c**********************************************************************
c                         hv calculation
c                         ______________
c
       hv = 1. - dmax1(0.d0,dsign(1.d0,qsat(ts(li), ps) - qa)) * rs *
     &      (1. - delta) / (rra + rs )
		
c
c?    rra ‚ a resistencia aerodinamica = 1/ChVa  -volta
c        a ser calculada … frente (correctamente


c        Implicit resolution of the surface temperature equation
c        linearization of Ts

      
        ztvi = ta * (1. + etv * qa)
        rhoa = ps / (ra * ztvi)
        qs = (hu * (1. - veg) + veg * hv) * qsat (ts(li),ps) +
     &       (1. - hv) * veg * qa
	  
c
c**********************************************************************
c                CMWF, 1981/11/25-27, pp. 59-79)
c**********************************************************************
c
        zdepl = 0.
        !write(*,*) 'nhli9904zref',zref
        zua=zref + zdepl
c
c corrected in 2005 (Victor Stepanenko)
c
        ztvis = ts(li) *(1. + etv * qs)
c______________________________________________________________________
c       calculo do numero de richardson e da resistencia aerodinamica
c
      !write(*,*) 'zref',zref
      ua=dmax1(1.d0,ua)
      ri=2.*gg*zua*zua*(ztvi-ztvis+gscp*zref)/(ztvi+ztvis)/(ua*ua)/zref
c
c       calculo do comprimento de Monin -Obukhov,de ustar,tstar
c
        valnon=-999.
c        call calclmon(z,zO,zOh,ua,ztvi,ztvls,xlmon,ustar,tstar,valnon)
        zeps= 1.e-6
        !write(* ,*) 'zua,z0',zua,z0
        zch=(vkarmn/dlog(abs(zua/z0)))**2
        !write(*,*) 'zch',zch
        zsta=ri*ua*ua
c
c  iyra escolhe entre a nova e a velha versao da rotina yra
c         de NP89. Quando iyra (ficheiro de entrada(?!)) e maior
c         do que 0 utiliza-se a nova versao que destingue entre
c         os comprimentos de rugosidade
c
        iyra=1
		
        if(iyra.eq.0) then
c
          if(ri.lt.0.) then
           zdi=1./(ua+75.*zch*sqrt(-zref*zsta/z0h))
           rra=zch*(ua-15.*zsta*zdi)
           zdim=1./(1.+75.*zch*sqrt(-zref*ri/z0))
           cdm=zch*(1.-10.*ri*zdim)*ua
          else
           zds=sqrt(ua * ua + 5. * zsta + zeps)
           rra = zch*ua/(1.+15.*zsta*zds/ua/ua/ua)
           zdim=(1.+5.*ri)**2.
           cdm=zch*ua/zdim
          endif
c
        else
c
c nova rotina yramng
          ! write(*,*) 'estez0h',z0h
          !write(*,*) 'z,zref',z,zref
           z=zref
           xmu=dlog(z0/z0h)
           xfh=dlog(z/z0)/dlog(z/z0h)
	
          ! write(*,*) '--------------------'
          ! write(*,*) 'z0,z0h',z0,z0h
           if(ri.lt.0.) then
                zdi=1./(ua+chstar(xmu)*zch*15*(z/z0h)**ph(xmu)*xfh
     s              *sqrt(-zsta))
                rra=zch*(ua-15.*zsta*zdi)*xfh
c
c alteracoes para o calculo do coeficiente de dragg_______________
c        cdm=Cd*uvsuf
c
                cdh=rra
                zdim=1./(ua+cmstar(xmu)*zch*10*(z/z0)**pm(xmu)
     s              *sqrt(-zsta))
                cdm=zch*(ua-10.*zsta*zdim)
c_______________________________________________________________________
           else
                zds = sqrt(ua * ua + 5. * zsta + zeps)
                rra = zch*ua/(1.+15.*zsta*zds/ua/ua/ua)*xfh
c
                cdh=rra
                zdsm= sqrt(ua * ua + 5. * zsta + zeps)
                cdm = zch*ua/(1.+10.*zsta/zdsm/ua)
           endif
c
        endif

c
        rra=1./rra
c----------------------------------------------------------------------
c
      zrsra = rhoa / rra
      if (iclay.ge.0) then
         ct = 1. / ((1. - veg) / cso + veg / cv)
         hv=1.-dmax1(0.d0,dsign(1.d0,qsat(ts(li), ps) - qa)) * rs*
     &    (1. - delta) / (rra + rs)
c
         za=1. / dt + ct * (4. * emis * stefan * (ts(li)**3) +
     &      zrsra*xl*dqsat(ts(li),ps)*(veg*hv+(1.-veg)*hu*xhu)
     &      +zrsra*cp)+2.*xpi/tcdil
         zb=1./dt+ ct*(3. * emis * stefan * (ts(li)** 3) +
     &      zrsra*xl*dqsat(ts(li),ps)*(veg*hv+(1.-veg)*hu*xhu))
         zc=2.*xpi*t2(li)/tcdil+ct*(zrsra*cp*ta+rg*
     &      (1.-alb)+emis*rat-zrsra*xl*(veg*hv*(qsat(ts(li),
     &      ps)-qa)+(1.-veg)*xhu*(hu*qsat(ts(li),ps)-qa)))
c
         ts(lf) = (ts(li) * zb + zc) / za
c
c       resolution of t2 equation
c
         t2(lf) = (t2(li) + dt * ts(lf) / tcdil) / (1. + dt / tcdil)
      else
          ts(lf)=tsfunc
          t2(lf)=t2(li)
      endif
c
c**********************************************************************
c


      hqs2=hu*qsat(ts(lf),ps)
      rn = rg * (1. - alb) + emis * (rat - stefan * (ts(lf)**4))
      h = rhoa * cp * (ts(lf) - ta) / rra
      leg = rhoa * xl * (1.-veg)*xhu*(hu* qsat(ts(lf), ps)-qa)/rra
      lev = rhoa * xl * veg * hv * (qsat(ts(lf), ps) - qa) / rra
      zzhv = dmax1(0.d0,dsign(1.d0,qsat(ts(li), ps) - qa))
      letr = zzhv * (1. - delta) * rhoa * xl * veg * (qsat(ts(lf), ps)
     &       - qa) / (rra + rs)
      le = leg + lev
      g = rn - h - le
      return
      end

      subroutine wri3grd(var,varname,istep,i0,i1,j0,j1,k0,k1,mode
     :  ,ixgrid,iygrid,isgrid
     :  ,nx,ny,ns,nstep,fngrd,timestring
     :  ,xmount,ymount,dx,dy,phi,phis,ixp,iyp,isp,zplan,nplan
     :  ,sigma0,sigma1,zconst,zlev,nlev,hsuf
     :  ,z,z2,z3,varxy,varxs,varys,varxz,varyz,iwxy,iwxz,iwyz)


c-----------------------------------------------------------------------
c handler for wrigrd (case of standard 3d arrays of nh3d)
c-----------------------------------------------------------------------
c  write a section of array a(l,m,n) into channel nchan
c  data to be printed out :                                             s)
c  ((a(i,j,k),i=i0,i1,is),j=j0,j1,js),k=k0,k1,ks)
c
c  leng :   maximum line length required (number of characters)
c
c  this geometry is appropriate for sigma vertical coordinate. in spite
c  of the arrangement of the j (y) coordinate, in all cases i0,j0 and k0
c  are always displayed, but not necessarily i1,j1,k1, depending on
c  is,js,ks.
c
c  p.m. nov 1988
c-----------------------------------------------------------------------


      implicit real*8(a-h,o-z)
c
c maxgr: max of grid size on x,y direction
c
      parameter(maxgr=1000)
      parameter(g=9.8066)
      character*80 fngrd,fextgrdt
      character*2 varname
      character*14 timestring
      dimension var(0:nx+1,0:ny+1,0:ns+1)
      dimension z(0:nx+1,0:ny+1,0:ns+1)
      dimension ixp(nplan),iyp(nplan),isp(nplan)
      dimension zplan(nplan)
      dimension zlev(0:nlev)
      dimension zbot(0:maxgr)
      dimension sigma0(0:ns+1)
      dimension sigma1(0:ns+1)
      dimension phi(0:nx+1,0:ny+1,0:ns+1)
      dimension phis(0:nx+1,0:ny+1,0:ns+1)
      dimension hsuf(0:nx+1,0:ny+1)
      dimension varxy(0:nx+1,0:ny+1)
      dimension varxs(0:nx+1,0:ns+1)
      dimension z2(0:nx+1,0:ns+1)
      dimension varxz(0:nx+1,0:nlev)
      dimension varys(0:ny+1,0:ns+1)
      dimension z3(0:ny+1,0:ns+1)
      dimension varyz(0:ny+1,0:nlev)
      dimension iwxy(0:nx+1,0:ny+1)
      dimension iwxz(0:nx+1,0:nlev)
      dimension iwyz(0:ny+1,0:nlev)
      logical zconst


      if(k0.eq.k1) then
         modd=1
      elseif(j0.eq.j1) then
         modd=2
      elseif(i0.eq.i1) then
         modd=3
      else
         modd=mode
      endif


      if(ixgrid.eq.0) then
        xmin=((i0-1.5)*dx-xmount)
        xmax=((i1-1.5)*dx-xmount)
      else
        xmin=((i0-1.)*dx-xmount)
        xmax=((i1-1.)*dx-xmount)
      endif


      if(iygrid.eq.0) then
        ymin=((j0-1.5)*dy-ymount)
        ymax=((j1-1.5)*dy-ymount)
      else
        ymin=((j0-1.)*dy-ymount)
        ymax=((j1-1.)*dy-ymount)
      endif


      if(isgrid.eq.0) then
        smin=sigma0(k1)
        smax=sigma0(k0)
      else
        smin=sigma1(k1)
        smax=sigma1(k0)
      endif


      if(zplan(1).ge.0. .or. (zconst.and.mode.ne.1)) then
        z=(phi+phis)/g
        if(isgrid.ne.0) then
          do is=0,ns
            z(:,:,is)=0.5*(z(:,:,is)+z(:,:,is+1))
          enddo
        endif
      endif
c
c write horizontal plane (k=cst):
c






      if(modd.eq.0 .or. modd.eq.1) then
        if(k0.eq.k1) then
c
c case of horizontal arrays
c
          call wrigrd(var(0,0,k0),xmin,xmax,ymin,ymax,vmin,vmax
     :      ,0,nx+1,0,ny+1,i0,i1,j0,j1,fextgrdt(fngrd,varname,'s'
     :      ,k0,timestring))
          return
        endif
        do kp=1,nplan
          if(isp(kp).ge.0) then
             k=isp(kp)
             do iy=0,ny+1
               do ix=0,nx+1
                 varxy(ix,iy)=var(ix,iy,k)
               enddo
             enddo
             call wrigrd(varxy,xmin,xmax,ymin,ymax,vmin,vmax
     :         ,0,nx+1,0,ny+1,i0,i1,j0,j1,fextgrdt(fngrd,varname,'s'
     :         ,k,timestring))
          endif
        enddo
        if(k0.ne.k1) then
          do kp=1,nplan
            if(zplan(kp).ge.0.) then
               call inter3(i0,i1,j0,j1,iwxy,z,zplan(kp),nx+1,ny+1,ns
     :           ,var,varxy)
               iplan=nint(zplan(kp)/100.)
               call wrigrd(varxy,xmin,xmax,ymin,ymax,vmin,vmax
     :           ,0,nx+1,0,ny+1,i0,i1,j0,j1,fextgrdt(fngrd,varname,'z'
     :           ,iplan,timestring))
            endif
          enddo
        endif
      endif
c
c write vertical plane (j=cst):
c
      if(mode.eq.0 .or. mode.eq.2) then
        do jp=1,nplan
          if(iyp(jp).ge.0) then
            j=iyp(jp)
            do is=0,ns+1
              do ix=0,nx+1
                varxs(ix,is)=var(ix,j,is)
              enddo
            enddo
            do ix=0,nx+1
              zbot(ix)=hsuf(ix,j)
            enddo
            if(zconst) then
              do ix=0,nx+1
                do is=0,ns+1
                  z2(ix,is)=z(ix,j,is)
                enddo
              enddo
              call inter2i(nx+1,ns+1,i0,i1,k0,k1,iwxz,z2,zlev,nlev
     :          ,varxs,varxz,zbot)
              call wrigrd(varxz,xmin,xmax,zlev(0),zlev(nlev),vmin,vmax
     :          ,0,nx+1,0,nlev,i0,i1,0,nlev,fextgrdt(fngrd,varname,'y'
     :          ,j,timestring))
            else
              call wrigrd(varxs,xmin,xmax,smax,smin,vmin,vmax
     :          ,0,nx+1,0,ns+1,i0,i1,k0,k1,fextgrdt(fngrd,varname,'y'
     :          ,j,timestring))
            endif
          endif
        enddo
      endif
c
c write vertical plane (i=cst):
c
      if(mode.eq.0 .or. mode.eq.3) then
        do ip=1,nplan
          if(ixp(ip).ge.0) then
            i=ixp(ip)
            do is=0,ns+1
              do iy=0,ny+1
                varys(iy,is)=var(i,iy,is)
              enddo
            enddo
            do iy=0,ny+1
              zbot(iy)=hsuf(i,iy)
            enddo
            if(zconst) then
              do iy=0,ny+1
                do is=0,ns+1
                  z3(iy,is)=z(i,iy,is)
                enddo
              enddo
              call inter2i(ny+1,ns+1,j0,j1,k0,k1,iwyz,z3,zlev,nlev
     :          ,varys,varyz,zbot)
              call wrigrd(varyz,ymin,ymax,zlev(0),zlev(nlev),vmin,vmax
     :          ,0,ny+1,0,nlev,j0,j1,0,nlev,fextgrdt(fngrd,varname,'x'
     :          ,i,timestring))
            else
              call wrigrd(varys,xmin,xmax,smax,smin,vmin,vmax
     :          ,0,ny+1,0,ns+1,j0,j1,k0,k1,fextgrdt(fngrd,varname,'x'
     :          ,i,timestring))
            endif
          endif
        enddo
      endif
c
      return
      end




      subroutine writbln(fnbln,zbot,nx,dx,xmount)
      implicit real*8 (a-h,o-z)
      character*80 fnbln
      dimension zbot(0:nx+1)
      integer nx
      open(77,file=fnbln,status='unknown')
      write(77,'(2i6)') nx+3,2
      do ix=1,nx
        xis=((ix-1.5)*dx-xmount)
        yps=zbot(ix)
        write(77,'(2e15.7)') xis,yps
      enddo


      write(77,'(2e15.7)') ((nx-1.5)*dx-xmount),0.
     :  ,((1-1.5)*dx-xmount),0.
     :  ,((1-1.5)*dx-xmount),zbot(1)
      return
      end




      subroutine inter3(ix0,ix1,iy0,iy1,isover,z,zlev,nxm,nym,ns,f,fa)
c
c interpolates a variable for a z grid, as given by zlev.
c if a given point on the regular z grid is under the topography
c fa is put to zero.
c
      implicit real*8(a-h,o-z)
      dimension isover(0:nxm,0:nym)
      dimension z(0:nxm,0:nym,0:ns),f(0:nxm,0:nym,0:ns),fa(0:nxm,0:nym)
c
      do iy=iy0,iy1
      do ix=ix0,ix1
        isover(ix,iy)=-999
      enddo
      enddo
c
      do 210 is=0,ns-1
      do 220 iy=iy0,iy1
      do 220 ix=ix0,ix1
      if(zlev.lt.z(ix,iy,is).and.zlev.ge.z(ix,iy,is+1)) isover(ix,iy)=is
220   continue
210   continue
c
      do 30 iy=iy0,iy1
      do 30 ix=ix0,ix1
      isov=isover(ix,iy)
      if(isov.eq.-999) then
         fa(ix,iy)=0.
      else
         delta=max(0.d0,(zlev-z(ix,iy,isov+1))
     :      /(z(ix,iy,isov)-z(ix,iy,isov+1)+1.e-30))
         fa(ix,iy)=(delta*f(ix,iy,isov)+(1.-delta)*f(ix,iy,isov+1))
      endif
30    continue
c
      return
      end




      subroutine inter2i(nim,njm,i0,i1,j0,j1,jover,z,zlev,nlev
     :   ,f,fa,zmin)


c
c interpolates a variable for a z grid, as given by zlev.
c zmin(ix): surface height
c
c note: nlev.le.njm
c
      implicit real*8(a-h,o-z)
      dimension zlev(0:*),zmin(0:*)
      dimension jover(0:nim,0:nlev)
      dimension z(0:nim,0:njm),f(0:nim,0:njm),fa(0:nim,0:nlev)
c
c      if(nlev.gt.njm) then
c         write(*,*) 'error in int2i',nlev
c         stop
c      endif
c
      do ilev=0,nlev
      do i=i0,i1
        jover(i,ilev)=-999
      enddo
      enddo
c
      do ilev=0,nlev
      do j=j0,j1-1
      do i=i0,i1
        if(zlev(ilev).le.z(i,j) .and. zlev(ilev).gt.z(i,j+1))
     :    jover(i,ilev)=j
      enddo
      enddo
      enddo


      do ilev=0,nlev
      do i=i0,i1
        jov=jover(i,ilev)
        if(jov.eq.-999) then
          fa(i,ilev)=f(i,j1)
        else
          if(zlev(ilev).ge.zmin(i)) then
            delta=max(0.d0,(zlev(ilev)-z(i,jov+1)))/(z(i,jov)-z(i,jov+1)
     :        +1.e-30)
            fa(i,ilev)=(delta*f(i,jov)+(1.-delta)*f(i,jov+1))
          else
            fa(i,ilev)=f(i,j1)
          endif
        endif
      enddo
      enddo
c
c      do ilev=0,nlev
c      do i=i0,i1
c        f(i,ilev)=fa(i,ilev)
c      enddo
c      enddo
c
      return
      end




      subroutine wrigrd(z,xmin,xmax,ymin,ymax,zmin,zmax
     :   ,m0,m1,n0,n1,i0,i1,j0,j1,fname)


      implicit real*8(a-h,o-z)
      character*80 fname
      dimension z(m0:m1,n0:n1)


c     write(*,'(a40,8i4,4f10.0)') fname,i0,i1,j0,j1
c    :  ,m0,m1,n0,n1,xmin,xmax,ymin,ymax
      zmin=z(i0,j0)
      zmax=z(i0,j0)


      do j=j0,j1
         do i=i0,i1
           zmin=min(zmin,z(i,j))
           zmax=max(zmax,z(i,j))
         enddo
      enddo


      open(77,file=fname,status='unknown')
      write(77,'(a4)') 'DSAA'
      write(77,'(2i5)') (i1-i0+1),(j1-j0+1)
      write(77,'(2e15.7)') xmin,xmax
      write(77,'(2e15.7)') ymin,ymax
      write(77,'(2e15.7)') zmin,zmax


      do j=j0,j1
         write(77,'(5E15.7)') (z(i,j),i=i0,i1)
      enddo


      close(77)
      return
      end




      subroutine pstvec(x,y,u,v,nunit)


      implicit real*8(a-h,o-z)
      parameter(pi=3.14159265)
      parameter(zero=1.e-20)
c
c     write(*,*)x,y,u,v
      if(max(abs(u),abs(v)).le.zero) return
      angle=atan2(v,u)*180./pi
      size=sqrt(u*u+v*v)
      write(nunit,'(3f8.2,e15.7)') x/1000.,y/1000.,angle,size
      return
      end




      subroutine readgrdful(fname,var,i0,i1,j0,j1,xmin,xmax,ymin,ymax)
c
c reads 2d array from a grd file
c
      implicit real*8(a-h,o-z)
      dimension var(i0:i1,j0:j1)
      character*80 fname
      character*4 dsaa


      nunit=77
      read(nunit,'(a4)') dsaa
      read(nunit,*) nx,ny
      read(nunit,*) xmin,xmax
      read(nunit,*) ymin,ymax
      read(nunit,*) zmin,zmax


      do j=j0,j1
         read(nunit,*) (var(i,j),i=i0,i1)
      enddo


      close(nunit)


      return
      end
