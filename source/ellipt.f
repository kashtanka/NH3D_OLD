      subroutine ellipt3(iellip,
     :  f
     :  ,dphi,el07,wrk2
     :  ,u000,dphids,wrk3
     :  ,v000,sdphds,wrk4
     :  ,el04
     :  ,uflux
     :  ,vflux
     :  ,wflux)
c-----------------------------------------------------------------------
c solve elliptic equation for geopotential height perturbation
c uses wk4,5,6,7,8
c
c a version of p.m.s ellipt
c can be launced if adjust.=.true.  only
c
c -- uses pos3c and inver to solve the poisson equation
c -- uses fixed boundary pressure p_*
c -- employes detailed vertically integrated mass balance condition for phi
c    and eigenvectors in vertical
c
c-----------------------------------------------------------------------
c
      use alloc
      implicit real*8(a-h,o-z)
      character*80 crash,crash2
c
      parameter (itma=10)
      parameter (epsloc=1.e-4)
c
c work variables:
c
      real(kind(0.d0)),dimension(0:nx1,0:ny1,0:ns1)::uflux,vflux,wflux
      dimension f(0:nx1,0:ny1,0:ns1),dphi(0:nx1,0:ny1,0:ns1)
      dimension b(0:nx1,0:ny1),bi(0:nx1,0:ny1)
c      dimension pphimx(0:nx1,0:ny1),pphimy(0:nx1,0:ny1)
c      dimension tt1(0:nx1,0:ny1),tt2(0:nx1,0:ny1),tt3(0:nx1,0:ny1)
      real(kind(0.d0)),allocatable,dimension(:,:)::
     :  pphimx,pphimy,tt1,tt2,tt3

      dimension u000(0:nx1,0:ny1,0:ns1),v000(0:nx1,0:ny1,0:ns1)
      dimension el04(0:nx1,0:ny1,0:ns1),el07(0:nx1,0:ny1,0:ns1)
      dimension dphids(0:nx1,0:ny1,0:ns),sdphds(0:nx1,0:ny1,0:ns)
      dimension wrk2(0:nx1,0:ny1,0:ns1)
      dimension wrk3(0:nx1,0:ny1,0:ns1)
      dimension wrk4(0:nx1,0:ny1,0:ns1)
c      dimension el02(0:nx1,0:ny1,0:ns1)
c      dimension el04a(0:nx1,0:ny1,0:ns1)
c      dimension el07a(0:nx1,0:ny1,0:ns1)
c
c --------------------------------------------------------
c
c      equivalence (wk1,f),(wk2,el07,dphi),(wk3,u000,dphids)
c     :   ,(wk4,v000,sdphds),(wk5,el04)
c
c      equivalence(hk1,ppgra2),(hk2,pplap),(hk3,gpplap)
c
c in-line funtions:
c
cCR put to zero f and dphi: really not neccessary
c
CR      call zero(f,nall)
      allocate(pphimx(0:nx1,0:ny1))
      allocate(pphimy(0:nx1,0:ny1))
      allocate(tt1(0:nx1,0:ny1))
      allocate(tt2(0:nx1,0:ny1))
      allocate(tt3(0:nx1,0:ny1))

      dphi=0.       ! to get dphi = 0 at level is = 0 (optional)
c
c compute the r.h.s source term of the equation, f
c and steady r.h.s. part of integral condition, b = div int_0^1 (F_u,F_v) ds
c
      do is=0,ns
      do iy=2,ny
      do ix=2,nx
        u000(ix,iy,is)=0.5*(u(ix,iy,is,3)+u(ix-1,iy,is,3))
        v000(ix,iy,is)=0.5*(v(ix,iy,is,3)+v(ix,iy-1,is,3))
        el04(ix,iy,is)=(ppdx(ix,iy,3)*(uflux(ix,iy,is)
     :     +uflux(ix-1,iy,is))
     :     +ppdy(ix,iy,3)*(vflux(ix,iy,is)+vflux(ix,iy-1,is)))*0.5
        el07(ix,iy,is)=-(ppdx(ix,iy,3)*(difunu(ix,iy,is)
     :     +difunu(ix-1,iy,is))
     :     +ppdy(ix,iy,3)*(difunv(ix,iy,is)+difunv(ix,iy-1,is)))*0.5
      enddo
      enddo
      enddo

      do is=1,ns-1
      do iy=2,ny
      do ix=2,nx

        p=ptop+sigma0(is)*pp(ix,iy,3)

        dtsdsi=(tems(ix,iy,is+1)-tems(ix,iy,is-1))/ds02a(is)

        alfa=ptop/p+sigma0(is)*dtsdsi/tems(ix,iy,is)

        beta=-ptop*pp(ix,iy,3)/(p*p)+dtsdsi/tems(ix,iy,is)
     :     +sigma0(is)*((tems(ix,iy,is+1)-tems(ix,iy,is))/ds0(is+1)
     :     -(tems(ix,iy,is)-tems(ix,iy,is-1))/ds0(is))/(ds1(is)
     :     *tems(ix,iy,is))
        betas=beta*sigma0(is)

        upduds=u000(ix,iy,is)+s0ds2a(is)
     :     *(u000(ix,iy,is+1)-u000(ix,iy,is-1))
        vpdvds=v000(ix,iy,is)+s0ds2a(is)
     :     *(v000(ix,iy,is+1)-v000(ix,iy,is-1))
c
c f(ix  ,iy,is)=(el02+el04+el05+el06+el07)/pp(ix,iy,3)+el03
c
        f(ix,iy,is)=(
c        el02(ix,iy,is)=
     :     -(uflux(ix,iy,is)-uflux(ix-1,iy,is))/dx
     :     -(vflux(ix,iy,is)-vflux(ix,iy-1,is))/dy
     :     +(s(ix,iy,is+1)*wflux(ix,iy,is+1)-s(ix,iy,is)
     :     *wflux(ix,iy,is))/ds1(is)
c        el04a(ix,iy,is)=
     :     +el04(ix,iy,is)
     :     +s0ds2a(is)*(el04(ix,iy,is+1)-el04(ix,iy,is-1))
c        el07a(ix,iy,is)=
     :     +(difunu(ix,iy,is)-difunu(ix-1,iy,is))/dx
     :     +(difunv(ix,iy,is)-difunv(ix,iy-1,is))/dy
     :     -(s(ix,iy,is+1)*difunw(ix,iy,is+1)-s(ix,iy,is)
     :     *difunw(ix,iy,is))/ds1(is)
     :     +el07(ix,iy,is)
     :     +s0ds2a(is)*(el07(ix,iy,is+1)-el07(ix,iy,is-1))
c
c        f(ix,iy,is)=(el02(ix,iy,is)+el04a(ix,iy,is)
c     :   +el07a(ix,iy,is)
     :    )/pp(ix,iy,3)
c
        f(ix,iy,is)=f(ix,iy,is)
c el03  =
     :     -g*(s(ix,iy,is+1)*((pt(ix,iy,is,3)+pt(ix,iy,is+1,3))
     :     /(pts(ix,iy,is)+pts(ix,iy,is+1))
     :     )-s(ix,iy,is)*((pt(ix,iy,is,3)+pt(ix,iy,is-1,3))
     :     /(pts(ix,iy,is)+pts(ix,iy,is-1))
     :     ))/ds1(is)
      enddo
      enddo
      enddo

c      if(elfout.and.prt) then
c         tk=0.
c         call wri3ar(u,2,nx,2,ny,1,ns-1,'u       ',tk,iomod,nx,ny)
c         call wri3ar(v,2,nx,2,ny,1,ns-1,'v       ',tk,iomod,nx,ny)
c         call wri3ar(w,2,nx,2,ny,1,ns-1,'w       ',tk,iomod,nx,ny)
c         call wri3ar(wsig,2,nx,2,ny,1,ns-1,'wsig    ',tk,iomod,nx,ny)
c         call wri3ar(uflux,2,nx,2,ny,1,ns-1,'uflux   ',tk,iomod,nx,ny)
c         call wri3ar(vflux,2,nx,2,ny,1,ns-1,'vflux   ',tk,iomod,nx,ny)
c         call wri3ar(wflux,2,nx,2,ny,1,ns-1,'wflux   ',tk,iomod,nx,ny)
c         call wri3ar(s,2,nx,2,ny,1,ns-1,'s       ',tk,iomod,nx,ny)
c         call wri3ar(el02,2,nx,2,ny,1,ns-1,'el02    ',tk,iomod,nx,ny)
c         call wri3ar(el04,2,nx,2,ny,1,ns-1,'el04a   ',tk,iomod,nx,ny)
c         call wri3ar(el04a,2,nx,2,ny,1,ns-1,'el04    ',tk,iomod,nx,ny)
c         call wri3ar(el07,2,nx,2,ny,1,ns-1,'el07    ',tk,iomod,nx,ny)
c         call wri3ar(el07a,2,nx,2,ny,1,ns-1,'el07a   ',tk,iomod,nx,ny)
c         call wri3ar(f,2,nx,2,ny,1,ns-1,'f-1     ',tk,iomod,nx,ny)
c      endif
c
C b for detailed mass balance
c
      do ix=2,nx
      do iy=2,ny
        sumb = 0.
        do  is=1,ns-1
        sumb= sumb+
     :  (
     :    (uflux(ix,iy,is)-uflux(ix-1,iy,is))/dx
     :   +(vflux(ix,iy,is)-vflux(ix,iy-1,is))/dy
     :  )*ds1(is)
        enddo
       b(ix,iy) = sumb
      enddo
      enddo
c     write(*,*) 'detailed mass balance'
c
      if(elfout.and.prt) then
         tk=0.
         call wri3ar(f,2,nx,2,ny,1,ns-1,'f       ',tk,iomod,nx,ny)
      endif
c
      if(abssuv(f,0,nx1,0,ny1,0,ns-1) .eq. 0.) return
c
c calculate the phi gradients at boundaries
c
c x left and right boundaries:
c
      do iy=2,ny
      do is=1,ns-1
        sav=0.25*(s(1,iy,is)+s(2,iy,is)+s(1,iy,is+1)+s(2,iy,is+1))
        xgrads=(0.25*((wflux(1,iy,is)+wflux(1,iy,is+1))
     :    /pp(1,iy,3)
     :    +(wflux(2,iy,is)+wflux(2,iy,is+1))/pp(2,iy,3))
     :    -g*(pt(1,iy,is,3)+pt(2,iy,is,3))
     :    /(pts(1,iy,is)+pts(2,iy,is)))/sav
        pgradx(1,iy,is)=((-(pp10(1,iy,4)*ubx(1,iy,is)
     :    -u(1,iy,is,2)*pp10(1,iy,2))/dtl
     :    -uflux(1,iy,is)
     :    +sigma0(is)*dpdx10(1,iy,3)*xgrads)/pp10(1,iy,3))*dx
      enddo
      enddo
c     write(*,*) '20 cont'
c
c x right boundary
c
      do iy=2,ny
      do is=1,ns-1
        sav=0.25*(s(nx,iy,is)+s(nx1,iy,is)+s(nx,iy,is+1)+s(nx1,iy,is+1))
c       xgrads(2,iy,is)
        xgrads=(0.25*((wflux(nx,iy,is)+wflux(nx,iy,is+1))
     :    /pp(nx,iy,3)
     :    +(wflux(nx1,iy,is)+wflux(nx1,iy,is+1))/pp(nx1,iy,3))
     :    -g*(pt(nx,iy,is,3)+pt(nx1,iy,is,3))
     :    /(pts(nx,iy,is)+pts(nx1,iy,is)))/sav
        pgradx(2,iy,is)=((-(pp10(nx,iy,4)*ubx(2,iy,is)
     :    -u(nx,iy,is,2)*pp10(nx,iy,2))/dtl
     :    -uflux(nx,iy,is)
     :    +sigma0(is)*dpdx10(nx,iy,3)*xgrads)/pp10(nx,iy,3))*dx
      enddo
      enddo
c     write(*,*) '30 cont'
c
      do ixb=1,2
      do iy=2,ny
        pgradx(ixb,iy,0)=pgradx(ixb,iy,1)
        pgradx(ixb,iy,ns)=pgradx(ixb,iy,ns-1)
      enddo
      enddo
c
c y left boundary:
c
      do is=1,ns-1
      do ix=2,nx
        sav=0.25*(s(ix,1,is)+s(ix,2,is)+s(ix,1,is+1)+s(ix,2,is+1))
        grads=(0.25*((wflux(ix,1,is)+wflux(ix,1,is+1))/pp(ix,1,3)
     :    +(wflux(ix,2,is)+wflux(ix,2,is+1))/pp(ix,2,3))
     :    -g*(pt(ix,1,is,3)+pt(ix,2,is,3))
     :    /(pts(ix,1,is)+pts(ix,2,is)))/sav
        pgrady(ix,1,is)=((-(pp01(ix,1,4)*vby(ix,1,is)
     :    -v(ix,1,is,2)*pp01(ix,1,2))/dtl
     :    -vflux(ix,1,is)
     :    +sigma0(is)*dpdy01(ix,1,3)*grads)/pp01(ix,1,3))*dy
      enddo
      enddo
c
c y right boundary:
c
      do is=1,ns-1
      do ix=2,nx
        sav=0.25*(s(ix,ny,is)+s(ix,ny1,is)+s(ix,ny,is+1)
     :    +s(ix,ny1,is+1))
        grads=(0.25*((wflux(ix,ny,is)+wflux(ix,ny,is+1))/pp(ix,ny,3)
     :    +(wflux(ix,ny1,is)+wflux(ix,ny1,is+1))/pp(ix,ny1,3))
     :    -g*(pt(ix,ny,is,3)+pt(ix,ny1,is,3))
     :    /(pts(ix,ny,is)+pts(ix,ny1,is)))/sav
        pgrady(ix,2,is)=((-(pp01(ix,ny,4)*vby(ix,2,is)
     :    -v(ix,ny,is,2)*pp01(ix,ny,2))/dtl
     :    -vflux(ix,ny,is)
     :    +sigma0(is)*dpdy01(ix,ny,3)*grads)/pp01(ix,ny,3))*dy
      enddo
      enddo

      do iyb=1,2
      do ix=2,nx
        pgrady(ix,iyb,0)=pgrady(ix,iyb,1)
        pgrady(ix,iyb,ns)=pgrady(ix,iyb,ns-1)
      enddo
      enddo
c
c bottom boundary:
c
      do iy=2,ny
      do ix=2,nx
        ppgra2(ix,iy)=(
     :    ppdx10(ix,iy,2)*ppdx10(ix,iy,2)+
     :    ppdx10(ix-1,iy,2)*ppdx10(ix-1,iy,2)+
     :    ppdy01(ix,iy,2)*ppdy01(ix,iy,2)+
     :    ppdy01(ix,iy-1,2)*ppdy01(ix,iy-1,2)
     :    )/2.
        pplap(ix,iy)=((pp(ix+1,iy,3)+pp(ix-1,iy,3)-2.*pp(ix,iy,3))/dxx
     :    +(pp(ix,iy+1,3)+pp(ix,iy-1,3)-2.*pp(ix,iy,3))/dyy)/pp(ix,iy,3)
        gpplap(ix,iy)=(ppgra2(ix,iy)-pplap(ix,iy))
        alfa=s(ix,iy,ns)*s(ix,iy,ns)
        tt1(ix,iy)=1./(1.+ppgra2(ix,iy)/alfa)*ds0(ns)
        tt2(ix,iy)=tt1(ix,iy)/alfa
      enddo
      enddo

      do iy=2,ny
      do ix=2,nx
c      pgrads(ix,iy,1)=0.
c      pgrads(ix,iy,1)=((wflux(ix,iy,1)/pp(ix,iy,3)
c     :   -g*(pt(ix,iy,0,3)+pt(ix,iy,1,3))/(pts(ix,iy,0)
c     :   +pts(ix,iy,1)))/s(ix,iy,1))*ds0(1)
        pgrads(ix,iy,2)=((wflux(ix,iy,ns)/pp(ix,iy,3)
     :    -g*(pt(ix,iy,ns-1,3)
     :    +pt(ix,iy,ns,3))/(pts(ix,iy,ns-1)+pts(ix,iy,ns))
c ***************************************************************
c additional term, proportional to horizontal flux at the boundary
     :    +(
     :    ppdx10(ix,iy,2)*uflux(ix,iy,ns)/pp10(ix,iy,3)
     :    +ppdx10(ix-1,iy,2)*uflux(ix-1,iy,ns)/pp10(ix-1,iy,3)
     :    +ppdy01(ix,iy,2)*vflux(ix,iy,ns)/pp01(ix,iy,3)
     :    +ppdy01(ix,iy-1,2)*vflux(ix,iy-1,ns)/pp01(ix,iy-1,3)
     :    )/s(ix,iy,ns)/2.
c ---------------------------------------------------------------
     :    )/s(ix,iy,ns))*tt1(ix,iy)
      enddo
      enddo

      call extrah(nx1,ny1,pgrads,1,nx1,1,ny1,1,2)
c
      if(elgout.and.prt) then
         tk=0.
         call wrigar(pgradx,1,2,0,ny1,0,ns,1,2,2,ny,0,ns
     :      ,'pgradx  ',tk,3)
         call wrigar(pgrady,0,nx1,1,2,0,ns,2,nx,1,2,0,ns
     :      ,'pgrady  ',tk,2)
         call wrigar(pgrads,0,nx1,0,ny1,1,2,2,nx,2,ny,1,2
     :      ,'pgrads  ',tk,1)
      endif
c     write(*,*) '80 cont+2'
c
c initialize pos3c: calculation of constants which do
c not vary from step to step
c
      if(iellip.eq.0)then
         call inpos3c
      endif
c     write(*,*) '80 cont+3'
c
c beginning of iterations
c
c iterations:
c
c     write(*,*) 'beginning of iterations'
      do 1000 ite=1,itma
c      call wri3ar(phi,2,nx,2,ny,0,ns,'phi-1   ',tk,1,nx,ny)
c
c compute the value of phi on the boundaries
c
      do is=1,ns-1
        do iy=2,ny
          phi(1,iy,is)=phi(2,iy,is)-pgradx(1,iy,is)
          phi(nx1,iy,is)=phi(nx,iy,is)+pgradx(2,iy,is)
        enddo
        do ix=2,nx
          phi(ix,1,is)=phi(ix,2,is)-pgrady(ix,1,is)
          phi(ix,ny1,is)=phi(ix,ny,is)+pgrady(ix,2,is)
        enddo
        phi(1,1,is)=phi(1,2,is)-pgrady(2,1,is)
        phi(nx1,1,is)=phi(nx1,2,is)-pgrady(nx,1,is)
        phi(1,ny1,is)=phi(1,ny,is)+pgrady(2,2,is)
        phi(nx1,ny1,is)=phi(nx1,ny,is)+pgrady(nx,2,is)
      enddo
c     write(*,*) '121a'
c
c *******************************************************
c      call wri3ar(phi,2,nx,2,ny,0,ns,'phi-+1  ',tk,1,nx,ny)
c      call wri3ar(pgrads(0,0,2),2,nx,2,ny,0,0,'pgrads  ',tk,1,nx,ny)
c      call wri3ar(ppdx10(0,0,2),2,nx,2,ny,0,0,'ppdx10  ',tk,1,nx,ny)
c      call wri3ar(ppdy01(0,0,2),2,nx,1,ny,0,0,'ppdy01  ',tk,1,nx,ny)

C BUG CORRECTED PROBLEM: this WAS recursive and introduceD assymmetries!!!
      do iy=2,ny
      do ix=2,nx
        tt3(ix,iy)=
     :    +(ppdx10(ix,iy,2)*(phi(ix+1,iy,ns)-phi(ix,iy,ns))/dx
     :    +ppdx10(ix-1,iy,2)*(phi(ix,iy,ns)-phi(ix-1,iy,ns))/dx
     :    +ppdy01(ix,iy,2)*(phi(ix,iy+1,ns)-phi(ix,iy,ns))/dy
     :    +ppdy01(ix,iy-1,2)*(phi(ix,iy,ns)-phi(ix,iy-1,ns))/dy
     :    )*tt2(ix,iy)/2.
      enddo
      enddo

      do iy=2,ny
      do ix=2,nx
        phi(ix,iy,ns)=phi(ix,iy,ns-1)+pgrads(ix,iy,2)+tt3(ix,iy)
c + phi hor. gadiendist soltuv osa
c    :    +(ppdx10(ix,iy,2)*(phi(ix+1,iy,ns)-phi(ix,iy,ns))/dx
c    :    +ppdx10(ix-1,iy,2)*(phi(ix,iy,ns)-phi(ix-1,iy,ns))/dx
c    :    +ppdy01(ix,iy,2)*(phi(ix,iy+1,ns)-phi(ix,iy,ns))/dy
c    :    +ppdy01(ix,iy-1,2)*(phi(ix,iy,ns)-phi(ix,iy-1,ns))/dy
c    :    )*tt2(ix,iy)/2.
      enddo
      enddo

c      call wri3ar(phi,2,nx,2,ny,0,ns,'phi+1   ',tk,1,nx,ny)
c      call wri3ar(tt2,2,nx,2,ny,0,0,'tt2     ',tk,1,nx,ny)
c     write(*,*) '121'
      do iy = 2,ny
c         phi(1,iy,0)=phi(2,iy,0)
         phi(1,iy,ns)=phi(1,iy,ns-1)+pgrads(2,iy,2)+
     :   ( ppdx10(2,iy,2)*(phi(3,iy,ns)-phi(2,iy,ns))/dx
     :    +ppdx10(1,iy,2)*(phi(2,iy,ns)-phi(1,iy,ns))/dx
     :    +ppdy01(2,iy,2)*(phi(2,iy+1,ns)-phi(2,iy,ns))/dy
     :    +ppdy01(2,iy-1,2)*(phi(2,iy,ns)-phi(2,iy-1,ns))/dy
     :   )*tt2(2,iy)/2.
c
c         phi(nx1,iy,0)=phi(nx,iy,0)
         phi(nx1,iy,ns)=phi(nx1,iy,ns-1)+pgrads(nx,iy,2)+
     :   ( ppdx10(nx,iy,2)*(phi(nx1,iy,ns)-phi(nx,iy,ns))/dx
     :    +ppdx10(nx-1,iy,2)*(phi(nx,iy,ns)-phi(nx-1,iy,ns))/dx
     :    +ppdy01(nx,iy,2)*(phi(nx,iy+1,ns)-phi(nx,iy,ns))/dy
     :    +ppdy01(nx,iy-1,2)*(phi(nx,iy,ns)-phi(nx,iy-1,ns))/dy
     :   )*tt2(2,iy)/2.
      enddo
c     write(*,*) '122'
         do 124 ix=2,nx
c         phi(ix,1,0)=phi(ix,2,0)
         phi(ix,1,ns)=phi(ix,1,ns-1)+pgrads(ix,2,2)+
     :   ( ppdx10(ix,2,2)*(phi(ix+1,2,ns)-phi(ix,2,ns))/dx
     :    +ppdx10(ix-1,2,2)*(phi(ix,2,ns)-phi(ix-1,2,ns))/dx
     :    +ppdy01(ix,2,2)*(phi(ix,3,ns)-phi(ix,2,ns))/dy
     :    +ppdy01(ix,1,2)*(phi(ix,2,ns)-phi(ix,1,ns))/dy
     :   )*tt2(ix,2)/2.
c
c         phi(ix,ny1,0)=phi(ix,ny,0)
         phi(ix,ny1,ns)=phi(ix,ny1,ns-1)+pgrads(ix,ny,2)+
     :   ( ppdx10(ix,ny,2)*(phi(ix+1,ny,ns)-phi(ix,ny,ns))/dx
     :    +ppdx10(ix-1,ny,2)*(phi(ix,ny,ns)-phi(ix-1,ny,ns))/dx
     :    +ppdy01(ix,ny,2)*(phi(ix,ny1,ns)-phi(ix,ny,ns))/dy
     :    +ppdy01(ix,ny-1,2)*(phi(ix,ny,ns)-phi(ix,ny-1,ns))/dy
     :   )*tt2(ix,ny)/2.
124      continue
c     write(*,*) '124'
c      phi(1,1,0)= phi(2,1,0)
      phi(1,1,ns)= phi(2,1,ns)
c      phi(1,ny1,0)= phi(2,ny1,0)
      phi(1,ny1,ns)= phi(2,ny1,ns)
c
c      phi(nx1,1,0)= phi(nx,1,0)
      phi(nx1,1,ns)= phi(nx,1,ns)
c      phi(nx1,ny1,0)= phi(nx,ny1,0)
      phi(nx1,ny1,ns)= phi(nx,ny1,ns)
c
c evaluate the l.h.s. from phi value of previous iteration
c and compute the residual of equation
c
      do is=1,ns
      do iy=1,ny1
      do ix=1,nx1
        dphids(ix,iy,is)=(phi(ix,iy,is)-phi(ix,iy,is-1))/ds0(is)
      enddo
      enddo
      enddo
c     write(*,*) '130'
c
      do is=1,ns-1
      do iy=1,ny1
      do ix=1,nx1
        sdphds(ix,iy,is)=sigma0(is)*
     :    (
     :    (phi(ix,iy,is+1)-phi(ix,iy,is))/ds0(is+1)
     :    +(phi(ix,iy,is)-phi(ix,iy,is-1))/ds0(is)
     :    )/2.
      enddo
      enddo
      enddo
c     write(*,*) '140'
c      call wri3ar(phi,2,nx,2,ny,0,ns,'phi0    ',tk,1,nx,ny)
c      call wri3ar(dphids,2,nx,2,ny,1,ns-1,'dphids  ',tk,iomod,nx,ny)
c      call wri3ar(dphids,2,nx,2,ny,1,ns-1,'ppgra2  ',tk,iomod,nx,ny)
c      call wri3ar(sdphds,2,nx,2,ny,1,ns-1,'sdphds  ',tk,iomod,nx,ny)
c      call wri3ar(ppdx,2,nx,2,ny,0,0,'ppdx    ',tk,iomod,nx,ny)
c      call wri3ar(ppdy,2,nx,2,ny,0,0,'ppdy    ',tk,iomod,nx,ny)
c
      do is=1,ns-1
      do iy=2,ny
      do ix=2,nx
        dphi(ix,iy,is)=f(ix,iy,is)-(
     :    (phi(ix+1,iy,is)+phi(ix-1,iy,is)-2.*phi(ix,iy,is))/dxx
     :    +(phi(ix,iy+1,is)+phi(ix,iy-1,is)-2.*phi(ix,iy,is))/dyy
     :    +(s(ix,iy,is+1)*s(ix,iy,is+1)*dphids(ix,iy,is+1)
     :    -s(ix,iy,is)*s(ix,iy,is)*dphids(ix,iy,is))/ds1(is)
     :    +ppgra2(ix,iy)*s0ds1(is)
     :    *(sigma1(is+1)*dphids(ix,iy,is+1)-sigma1(is)*dphids(ix,iy,is))
     :    +gpplap(ix,iy)*sdphds(ix,iy,is)
     :    -2.*ppdx(ix,iy,3)*(sdphds(ix+1,iy,is)-sdphds(ix-1,iy,is))/dx2
     :    -2.*ppdy(ix,iy,3)*(sdphds(ix,iy+1,is)-sdphds(ix,iy-1,is))/dy2
     :    )
      enddo
      enddo
      enddo
c
C Right hand term bi for mass balance condition:
C ---------------------------------------------
c
c    1.  (pphimx,pphimy)=p_*[sum_1^(ns-1)nabla_p(phi)ds]
c
      do ix = 1,nx
      do iy = 1,ny
        sumx = 0.
        sumy = 0.
        do  is=1,ns-1
c
        sumx= sumx+ds1(is)*(
     +   (phi(ix+1,iy,is)-phi(ix,iy,is))/dx
     -   -s0ds4a(is)*ppdx10(ix,iy,2)*(phi(ix+1,iy,is+1)
     +   +phi(ix,iy,is+1)-phi(ix+1,iy,is-1)-phi(ix,iy,is-1))
     )                     )
c
        sumy= sumy+ds1(is)*(
     +    (phi(ix,iy+1,is)-phi(ix,iy,is))/dy
     -   -s0ds4a(is)*ppdy01(ix,iy,2)*(phi(ix,iy+1,is+1)
     +   +phi(ix,iy,is+1)-phi(ix,iy+1,is-1)-phi(ix,iy,is-1))
     )                     )
        enddo
       pphimx(ix,iy) = sumx*pp10(ix,iy,3)
       pphimy(ix,iy) = sumy*pp01(ix,iy,3)
      enddo
      enddo
c
c    2.   b + \nabla *pphim
c
      do ix = 2,nx
      do iy = 2,ny
       bi(ix,iy) = b(ix,iy) +
     + (pphimx(ix,iy)-pphimx(ix-1,iy))/dx
     + +
     + (pphimy(ix,iy)-pphimy(ix,iy-1))/dy
      enddo
      enddo
c
      do iy = 2,ny
       bi(1,iy) = bi(2,iy)
       bi(nx1,iy) = bi(nx,iy)
      enddo
      do ix = 1,nx1
       bi(ix,1) = bi(ix,2)
       bi(ix,ny1)=bi(ix,ny)
      enddo
c
c solve the standard poisson equation with homogeneous b.c.'s
c
      if(eldout.and.prt) then
         tk=0.
         call wri3ar(dphi,2,nx,2,ny,1,ns-1,'dphi0   ',tk,iomod,nx,ny)
      endif
c
c     write(*,*) 'call pos3c'
      call pos3c(bi,dphi,wrk3,wrk4,ite)
c
c update the phi value after one iteration and exam the convergence
c
      do is=0,ns1
      do iy=1,ny1
      do ix=1,nx1
        phi(ix,iy,is)=phi(ix,iy,is)+dphi(ix,iy,is)
      enddo
      enddo
      enddo
c
      phisum=abssuv(phi,0,nx1,0,ny1,0,ns)
      dpsum=abssuv(dphi,0,nx1,0,ny1,0,ns)
c
      if(eldout.and.prt .or. echeck) then
         write(nchn1,2002) nstep,ite,phisum,dpsum
         tk=0.
         call wri3ar(dphi,1,nx1,1,ny1,0,ns,'dphi    ',tk,iomod,nx,ny)
         call wri3ar(phi,1,nx1,1,ny1,0,ns,'phie    ',tk,iomod,nx,ny)
      endif
c
      if(phisum .eq. 0.) phisum=1.
CR      if(dpsum/phisum .le. eps) then
      if(dpsum/phisum .le. epsloc) then
CR         write(*,2000) nstep,ite
         write(*,2002) nstep,ite,phisum,dpsum
         return
      endif
c
1000  continue
c     write(*,*) '1000 cont'
c
      if(echeck) then
         write(nchn1,2001) nstep,itma
         tk=0.
         call wri3ar(f,2,nx,2,ny,1,ns-1,'f       ',tk,iomod,nx,ny)
         call wrigar(pgradx,1,2,0,ny1,0,ns,1,2,1,ny1,0,ns
     :      ,'pgradx  ',tk,3)
         call wrigar(pgrady,0,nx1,1,2,0,ns,1,nx1,1,2,0,ns
     :      ,'pgrady  ',tk,2)
         call wrigar(pgrads,0,nx1,0,ny1,1,2,2,nx,2,ny,1,2
     :      ,'pgrads  ',tk,1)
         crash='crash.out'
         crash2='crsh'
         call out6(0,3,crash,crash2)
         stop
      endif
c
CR      write(*,2000) nstep,ite
         write(*,2002) nstep,ite,phisum,dpsum


      deallocate(pphimx)
      deallocate(pphimy)
      deallocate(tt1)
      deallocate(tt2)
      deallocate(tt3)
      stop
      return
2000  format(' ellipt: nstep=   ',i4,'   ite=   ',i3)
2001  format(' ellipt2 did not converge at nstep=',i5,' after ',i3,
     :   ' iterations')
2002  format(' nh3dFatalError=ellipt1: nstep=',i4,' ite=',i3
     :   ,' phisum=',e14.7,' dpsum=',e14.7)
      end
      
      subroutine ellipt(
     :  f
     :  ,dphi,el07,wrk2
     :  ,u000,dphids,wrk3
     :  ,v000,sdphds,wrk4
     :  ,el04,wrk5
     :  ,uflux,vflux,wflux)

c      equivalence (wk1,f),(wk2,el07,dphi),(wk3,u000,dphids)
c     :   ,(wk4,v000,sdphds),(wk5,el04)
c-----------------------------------------------------------------------
c solve elliptic equation for geopotential height perturbation
c uses wk4,5,6,7,8
c-----------------------------------------------------------------------

      use alloc

      implicit real*8(a-h,o-z)
      real(kind(0.d0)),dimension(0:nx1,0:ny1,0:ns1)::uflux,vflux,wflux
      dimension f(0:nx1,0:ny1,0:ns1)
      dimension dphi(0:nx1,0:ny1,0:ns1)
      dimension el07(0:nx1,0:ny1,0:ns1)
      dimension u000(0:nx1,0:ny1,0:ns1)
      dimension dphids(0:nx1,0:ny1,0:ns1)
      dimension v000(0:nx1,0:ny1,0:ns1)
      dimension sdphds(0:nx1,0:ny1,0:ns1)
      dimension el04(0:nx1,0:ny1,0:ns1)
      dimension wrk2(0:nx1,0:ny1,0:ns1)
      dimension wrk3(0:nx1,0:ny1,0:ns1)
      dimension wrk4(0:nx1,0:ny1,0:ns1)
      dimension wrk5(0:nx1,0:ny1,0:ns1)
      character*80 crash,crash2
c      include 'nh3dpa10.inc'

c
      parameter (itma=20)
c
c work variables:
c
c      real lhs1(0:nx1,0:ny1,0:ns1),lhs2(0:nx1,0:ny1,0:ns1)
c      real lhs3(0:nx1,0:ny1,0:ns1),lhs4(0:nx1,0:ny1,0:ns1)
c      real lhs5(0:nx1,0:ny1,0:ns1),lhs6(0:nx1,0:ny1,0:ns1)
c      real el01(0:nx1,0:ny1,0:ns1),el02(0:nx1,0:ny1,0:ns1)
c      real el07a(0:nx1,0:ny1,0:ns1),el04a(0:nx1,0:ny1,0:ns1)
c      real el05(0:nx1,0:ny1,0:ns1),el06(0:nx1,0:ny1,0:ns1)
c
c      dimension f(0:nx1,0:ny1,0:ns1),dphi(0:nx1,0:ny1,0:ns1)
c      dimension xgrads(2,0:ny1,0:ns),ygrads(0:nx1,2,0:ns)
c      dimension u000(0:nx1,0:ny1,0:ns1),v000(0:nx1,0:ny1,0:ns1)
c      dimension el04(0:nx1,0:ny1,0:ns1),el07(0:nx1,0:ny1,0:ns1)
c      dimension dphids(0:nx1,0:ny1,0:ns),sdphds(0:nx1,0:ny1,0:ns)
c      dimension ppgra2(0:nx1,0:ny1),pplap(0:nx1,0:ny1)
c     :   ,gpplap(0:nx1,0:ny1)
c
c      equivalence (wk1,f),(wk2,el07,dphi),(wk3,u000,dphids)
c     :   ,(wk4,v000,sdphds),(wk5,el04)
c
c      equivalence(hk1,ppgra2),(hk2,pplap),(hk3,gpplap)
c
c put to zero f and dphi:
c
c      allocate(u000(0:nx1,0:ny1,0:ns1))
c      allocate(v000(0:nx1,0:ny1,0:ns1))
c      allocate(el04(0:nx1,0:ny1,0:ns1))
c      allocate(el07(0:nx1,0:ny1,0:ns1))
c
c compute the r.h.s source term of the equation
c
c       write(*,*) 'Enter ellipt',nstep
      do is=0,ns
      do iy=2,ny
      do ix=2,nx
        u000(ix,iy,is)=0.5*(u(ix,iy,is,3)+u(ix-1,iy,is,3))
        v000(ix,iy,is)=0.5*(v(ix,iy,is,3)+v(ix,iy-1,is,3))
        el04(ix,iy,is)=(ppdx(ix,iy,3)*(uflux(ix,iy,is)
     :    +uflux(ix-1,iy,is))
     :    +ppdy(ix,iy,3)*(vflux(ix,iy,is)+vflux(ix,iy-1,is)))*0.5
        el07(ix,iy,is)=-(ppdx(ix,iy,3)*(difunu(ix,iy,is)
     :    +difunu(ix-1,iy,is))
     :    +ppdy(ix,iy,3)*(difunv(ix,iy,is)+difunv(ix,iy-1,is)))*0.5
      enddo
      enddo
      enddo

c       write(*,*) 'Allocate f'
c      allocate(f(0:nx1,0:ny1,0:ns1))
      f=0.

      do is=1,ns-1
      do iy=2,ny
      do ix=2,nx
        p=ptop+sigma0(is)*pp(ix,iy,3)
        dtsdsi=(tems(ix,iy,is+1)-tems(ix,iy,is-1))/(ds02a(is)
     :    *tems(ix,iy,is))
        alfa=ptop/p+sigma0(is)*dtsdsi
c
c Corrected (version 9601)
c
        beta=-ptop*pp(ix,iy,3)/(p*p)+dtsdsi*(1.+sigma0(is)*dtsdsi)
     :    +sigma0(is)*((tems(ix,iy,is+1)-tems(ix,iy,is))/ds0(is+1)
     :    -(tems(ix,iy,is)-tems(ix,iy,is-1))/ds0(is))/(ds1(is)
     :    *tems(ix,iy,is))
        betas=beta*sigma0(is)
c
        upduds=u000(ix,iy,is)+s0ds2a(is)
     :    *(u000(ix,iy,is+1)-u000(ix,iy,is-1))
        vpdvds=v000(ix,iy,is)+s0ds2a(is)
     :    *(v000(ix,iy,is+1)-v000(ix,iy,is-1))
c
c f(ix,iy,is)=(el01+el02+el04+el05+el06+el07)/pp(ix,iy,3)+el03
c
        f(ix,iy,is)=
c el01(ix,iy,is)=
     :    ((pp(ix,iy,3)-pp(ix,iy,1))/dtl
     :    +(u(ix,iy,is,2)*pp10(ix,iy,2)-u(ix-1,iy,is,2)*pp10(ix-1,iy,2))
     :    /dx
     :    +(v(ix,iy,is,2)*pp01(ix,iy,2)-v(ix,iy-1,is,2)*pp01(ix,iy-1,2))
     :    /dy
     :    +pp(ix,iy,2)*(wsig(ix,iy,is+1,2)-wsig(ix,iy,is,2))/ds1(is))
     :    /dt2
c el02(ix,iy,is)=
     :    -(uflux(ix,iy,is)-uflux(ix-1,iy,is))/dx
     :    -(vflux(ix,iy,is)-vflux(ix,iy-1,is))/dy
     :    +(s(ix,iy,is+1)*wflux(ix,iy,is+1)-s(ix,iy,is)
     :    *wflux(ix,iy,is))/ds1(is)
c el04a(ix,iy,is)=
     :    +el04(ix,iy,is)
     :    +s0ds2a(is)*(el04(ix,iy,is+1)-el04(ix,iy,is-1))
c el05(ix,iy,is)=
     :    -upduds*(dpp(ix+1,iy)-dpp(ix-1,iy))/dx2
     :    -vpdvds*(dpp(ix,iy+1)-dpp(ix,iy-1))/dy2
c
        f(ix,iy,is)=(f(ix,iy,is)
c
c Corrected (version 9601)
c
c el06(ix,iy,is)=
     :    +dpp(ix,iy)
     :    *(-(alfa*(wsig(ix,iy,is+1,3)-wsig(ix,iy,is,3))/ds1(is)
     :    +beta*(wsig(ix,iy,is,3)+wsig(ix,iy,is+1,3))*0.5)
     :    -(alfa+betas)*dpp(ix,iy)/pp(ix,iy,3)
     :    +(upduds*(1.-alfa)-betas*u000(ix,iy,is))*ppdx(ix,iy,3)
     :    +(vpdvds*(1.-alfa)-betas*v000(ix,iy,is))*ppdy(ix,iy,3))
c el07a(ix,iy,is)=
     :    +(difunu(ix,iy,is)-difunu(ix-1,iy,is))/dx
     :    +(difunv(ix,iy,is)-difunv(ix,iy-1,is))/dy
     :    -(s(ix,iy,is+1)*difunw(ix,iy,is+1)-s(ix,iy,is)
     :    *difunw(ix,iy,is))/ds1(is)
     :    +el07(ix,iy,is)
     :    +s0ds2a(is)*(el07(ix,iy,is+1)-el07(ix,iy,is-1))
     :    )/pp(ix,iy,3)
c
c       if(qif.eq.0.) then
c         if(is.ge.34) then
c          write(*,'(3i4,6e10.3)') ix,iy,is,f(ix,iy,is)
c     :      ,s(ix,iy,is),s(ix,iy,is+1)
c         write(*,'(3i4,6e10.3)') ix,iy,is,f(ix,iy,is)
c         write(*,'(3i4,6e10.3)') ix,iy,is,s(ix,iy,is)
c         write(*,'(3i4,6e10.3)') ix,iy,is,s(ix,iy,is+1)
c         write(*,'(3i4,6e10.3)') ix,iy,is,thetav(ix,iy,is)
c         write(*,'(3i4,6e10.3)') ix,iy,is,thetav(ix,iy,is+1)
c         write(*,'(3i4,6e10.3)') ix,iy,is,ds1(is)
c         endif
          f(ix,iy,is)=f(ix,iy,is)
     :      -(s(ix,iy,is+1)*thetav(ix,iy,is+1)
     :      -s(ix,iy,is)*thetav(ix,iy,is))/ds1(is)

cc el03=
c     :      -g*(s(ix,iy,is+1)*((pt(ix,iy,is,3)+pt(ix,iy,is+1,3))
c     :      /(pts(ix,iy,is)+pts(ix,iy,is+1))
c     :      )-s(ix,iy,is)*((pt(ix,iy,is,3)+pt(ix,iy,is-1,3))
c     :      /(pts(ix,iy,is)+pts(ix,iy,is-1))
c     :      ))/ds1(is)
c        elseif(ifqc.eq.0) then
c          f(ix,iy,is)=f(ix,iy,is)
cc
c
cc el03=
c     :      -g*(s(ix,iy,is+1)*(
c     :      (pt(ix,iy,is,3)+pt(ix,iy,is+1,3))
c     :      /(pts(ix,iy,is)+pts(ix,iy,is+1))
c     :      +qdrag*0.5*(0.61*(qv(ix,iy,is,3)+qv(ix,iy,is+1,3)
c     :      -qvs(ix,iy,is)-qvs(ix,iy,is+1)))
c     :      )-s(ix,iy,is)*
c     :      ((pt(ix,iy,is,3)+pt(ix,iy,is-1,3))
c     :      /(pts(ix,iy,is)+pts(ix,iy,is-1))
c     :      +qdrag*0.5*(0.61*(qv(ix,iy,is,3)+qv(ix,iy,is-1,3)
c     :      -qvs(ix,iy,is)-qvs(ix,iy,is-1)))
c     :      ))/ds1(is)
c        elseif(ifqr.eq.0) then
c          f(ix,iy,is)=f(ix,iy,is)
cc el03=
c     :      -g*(s(ix,iy,is+1)*((pt(ix,iy,is,3)+pt(ix,iy,is+1,3))
c     :      /(pts(ix,iy,is)+pts(ix,iy,is+1))
c     :      +qdrag*0.5*(0.61*(qv(ix,iy,is,3)+qv(ix,iy,is+1,3)
c     :      -qvs(ix,iy,is)-qvs(ix,iy,is+1))
c     :      -(qc(ix,iy,is,3)+qc(ix,iy,is+1,3)))
c     :      )-s(ix,iy,is)*((pt(ix,iy,is,3)+pt(ix,iy,is-1,3))
c     :      /(pts(ix,iy,is)+pts(ix,iy,is-1))

c     :      +qdrag*0.5*(0.61*(qv(ix,iy,is,3)+qv(ix,iy,is-1,3)
c     :      -qvs(ix,iy,is)-qvs(ix,iy,is-1))
c     :      -(qc(ix,iy,is,3)+qc(ix,iy,is-1,3)))
c     :      ))/ds1(is)
c        else
c          f(ix,iy,is)=f(ix,iy,is)
cc el03=
c     :      -g*(s(ix,iy,is+1)*((pt(ix,iy,is,3)+pt(ix,iy,is+1,3))
c     :      /(pts(ix,iy,is)+pts(ix,iy,is+1))
c     :      +qdrag*0.5*(0.61*(qv(ix,iy,is,3)+qv(ix,iy,is+1,3)
c     :      -qvs(ix,iy,is)-qvs(ix,iy,is+1))
c     :      -(qc(ix,iy,is,3)+qc(ix,iy,is+1,3))
c     :      -(qr(ix,iy,is,3)+qr(ix,iy,is+1,3)))
c     :      )-s(ix,iy,is)*((pt(ix,iy,is,3)+pt(ix,iy,is-1,3))
c     :      /(pts(ix,iy,is)+pts(ix,iy,is-1))
c     :      +qdrag*0.5*(0.61*(qv(ix,iy,is,3)+qv(ix,iy,is-1,3)
c     :      -qvs(ix,iy,is)-qvs(ix,iy,is-1))
c     :      -(qc(ix,iy,is,3)+qc(ix,iy,is-1,3))
c     :      -(qr(ix,iy,is,3)+qr(ix,iy,is-1,3)))
c     :      ))/ds1(is)
c        endif
      enddo
      enddo
      enddo
c
c      deallocate(u000)
c      deallocate(v000)
c      deallocate(el04)
c      deallocate(el07)
c       write(*,*) '10 done'
c
      if(elfout.and.prt) then
         tk=0.
         call wri3ar(f,2,nx,2,ny,1,ns-1,'f       ',tk,iomod,nx,ny)
      endif
c
      if(abssuv(f,0,nx1,0,ny1,0,ns-1) .eq. 0.) then
c        deallocate(f)
        return
      endif
c
c calculate the phi gradients at boundaries
c
c
c x left and right boundaries:
c
      do iy=2,ny
      do is=1,ns-1
        sav=0.25*(s(1,iy,is)+s(2,iy,is)+s(1,iy,is+1)+s(2,iy,is+1))
c xgrads(1,iy,is)=
        xgrads=(0.25*((wflux(1,iy,is)+wflux(1,iy,is+1))
     :    /pp(1,iy,3)
     :    +(wflux(2,iy,is)+wflux(2,iy,is+1))/pp(2,iy,3))
     :    -g*(pt(1,iy,is,3)+pt(2,iy,is,3))
     :    /(pts(1,iy,is)+pts(2,iy,is)))/sav
        pgradx(1,iy,is)=((-(pp10(1,iy,4)*ubx(1,iy,is)
     :    -u(1,iy,is,2)*pp10(1,iy,2))/dtl-uflux(1,iy,is)
     :    +sigma0(is)*dpdx10(1,iy,3)*xgrads)/pp10(1,iy,3))*dx
      enddo
      enddo
c
c x right boundary
c
      do iy=2,ny
      do is=1,ns-1
        sav=0.25*(s(nx,iy,is)+s(nx1,iy,is)+s(nx,iy,is+1)+s(nx1,iy,is+1))
c xgrads(2,iy,is)
        xgrads=(0.25*((wflux(nx,iy,is)+wflux(nx,iy,is+1))
     :    /pp(nx,iy,3)
     :    +(wflux(nx1,iy,is)+wflux(nx1,iy,is+1))/pp(nx1,iy,3))
     :    -g*(pt(nx,iy,is,3)+pt(nx1,iy,is,3))
     :    /(pts(nx,iy,is)+pts(nx1,iy,is)))/sav
        pgradx(2,iy,is)=((-(pp10(nx,iy,4)*ubx(2,iy,is)
     :    -u(nx,iy,is,2)*pp10(nx,iy,2))/dtl-uflux(nx,iy,is)
     :    +sigma0(is)*dpdx10(nx,iy,3)*xgrads)/pp10(nx,iy,3))*dx
      enddo
      enddo
c
      do ixb=1,2
      do iy=2,ny
        pgradx(ixb,iy,0)=pgradx(ixb,iy,1)
        pgradx(ixb,iy,ns)=pgradx(ixb,iy,ns-1)
      enddo
      enddo
c
c y left boundary:
c
      do is=1,ns-1
      do ix=2,nx
        sav=0.25*(s(ix,1,is)+s(ix,2,is)+s(ix,1,is+1)+s(ix,2,is+1))
        grads=(0.25*((wflux(ix,1,is)+wflux(ix,1,is+1))/pp(ix,1,3)
     :    +(wflux(ix,2,is)+wflux(ix,2,is+1))/pp(ix,2,3))
     :    -g*(pt(ix,1,is,3)+pt(ix,2,is,3))
     :    /(pts(ix,1,is)+pts(ix,2,is)))/sav
        pgrady(ix,1,is)=((-(pp01(ix,1,4)*vby(ix,1,is)
     :    -v(ix,1,is,2)*pp01(ix,1,2))/dtl-vflux(ix,1,is)
     :    +sigma0(is)*dpdy01(ix,1,3)*grads)/pp01(ix,1,3))*dy
      enddo
      enddo
c
c y right boundary:
c
      do is=1,ns-1
      do ix=2,nx
        sav=0.25*(s(ix,ny,is)+s(ix,ny1,is)+s(ix,ny,is+1)+s(ix,ny1,is+1))
        grads=(0.25*((wflux(ix,ny,is)+wflux(ix,ny,is+1))/pp(ix,ny,3)
     :    +(wflux(ix,ny1,is)+wflux(ix,ny1,is+1))/pp(ix,ny1,3))
     :    -g*(pt(ix,ny,is,3)+pt(ix,ny1,is,3))
     :    /(pts(ix,ny,is)+pts(ix,ny1,is)))/sav
        pgrady(ix,2,is)=((-(pp01(ix,ny,4)*vby(ix,2,is)
     :    -v(ix,ny,is,2)*pp01(ix,ny,2))/dtl-vflux(ix,ny,is)
     :    +sigma0(is)*dpdy01(ix,ny,3)*grads)/pp01(ix,ny,3))*dy
      enddo
      enddo
      do iyb=1,2
      do ix=2,nx
        pgrady(ix,iyb,0)=pgrady(ix,iyb,1)
        pgrady(ix,iyb,ns)=pgrady(ix,iyb,ns-1)
      enddo
      enddo
c
c top and bottom boundaries:
c
      do iy=2,ny
      do ix=2,nx
        pgrads(ix,iy,1)=((wflux(ix,iy,1)/pp(ix,iy,3)
     :    -g*(pt(ix,iy,0,3)+pt(ix,iy,1,3))/(pts(ix,iy,0)
     :    +pts(ix,iy,1)))/s(ix,iy,1))*ds0(1)
        pgrads(ix,iy,2)=((wflux(ix,iy,ns)/pp(ix,iy,3)
     :    -g*(pt(ix,iy,ns-1,3)+pt(ix,iy,ns,3))
     :    /(pts(ix,iy,ns-1)+pts(ix,iy,ns)))/s(ix,iy,ns))*dsng
      enddo
      enddo
c
      call extrah(nx1,ny1,pgrads,1,nx1,1,ny1,1,2)
c
      if(elgout.and.prt) then
         tk=0.
         call wrigar(pgradx,1,2,0,ny1,0,ns,1,2,2,ny,0,ns
     :      ,'pgradx  ',tk,3)
         call wrigar(pgrady,0,nx1,1,2,0,ns,2,nx,1,2,0,ns
     :      ,'pgrady  ',tk,2)
         call wrigar(pgrads,0,nx1,0,ny1,1,2,2,nx,2,ny,1,2
     :      ,'pgrads  ',tk,1)
      endif
c
c begining of iterations
c
c initialize pos3 (calculation of temporary variables which do
c not vary between iterations and preparation of matrix solution)
c
c       write(*,*) 'Call inpos3'
c      allocate(work1(0:nx1,0:ny1,0:ns1))
c      allocate(work2(0:nx1,0:ny1,0:ns1))
      call inpos3(wrk5,wrk2)
c      deallocate(work2)
c      allocate(ppgra2(0:nx1,0:ny1))
c      allocate(pplap(0:nx1,0:ny1))
c      allocate(gpplap(0:nx1,0:ny1))
c       write(*,*) 'done inpos3'
c
      do iy=2,ny
      do ix=2,nx
        ppgra2(ix,iy)=ppdx(ix,iy,3)*ppdx(ix,iy,3)+ppdy(ix,iy,3)
     :    *ppdy(ix,iy,3)
        pplap(ix,iy)=((pp(ix+1,iy,3)+pp(ix-1,iy,3)-2.*pp(ix,iy,3))/dxx
     :   +(pp(ix,iy+1,3)+pp(ix,iy-1,3)-2.*pp(ix,iy,3))/dyy)/pp(ix,iy,3)
        gpplap(ix,iy)=(ppgra2(ix,iy)-pplap(ix,iy))
      enddo
      enddo
c
c      if(eldout.and.prt) then
c         tk=0.
c         call wri2ar(ppdx(0,0,3),1,nx1,1,ny1,'ppdx    ',tk,nx,ny)
c         call wri2ar(ppdy(0,0,3),1,nx1,1,ny1,'ppdy    ',tk,nx,ny)
c         call wri2ar(ppgra2,2,nx,2,ny,'ppgra2  ',tk,nx,ny)
c         call wri2ar(gpplap,2,nx,2,ny,'gpplap  ',tk,nx,ny)
c      endif
c
c iterations:
c
c      allocate(dphi(0:nx1,0:ny1,0:ns1))

c       write(*,*) 'Iterate'
      do 1000 ite=1,itma
c
c compute the value of phi on the boundaries
c
      do is=1,ns-1
        do iy=2,ny
          phi(1,iy,is)=phi(2,iy,is)-pgradx(1,iy,is)
          phi(nx1,iy,is)=phi(nx,iy,is)+pgradx(2,iy,is)
        enddo
        do ix=2,nx
          phi(ix,1,is)=phi(ix,2,is)-pgrady(ix,1,is)
          phi(ix,ny1,is)=phi(ix,ny,is)+pgrady(ix,2,is)
        enddo
        phi(1,1,is)=phi(1,2,is)-pgrady(2,1,is)
        phi(nx1,1,is)=phi(nx1,2,is)-pgrady(nx,1,is)
        phi(1,ny1,is)=phi(1,ny,is)+pgrady(2,2,is)
        phi(nx1,ny1,is)=phi(nx1,ny,is)+pgrady(nx,2,is)
      enddo
c
      do iy=1,ny1
      do ix=1,nx1
        phi(ix,iy,0)=phi(ix,iy,1)-pgrads(ix,iy,1)
        phi(ix,iy,ns)=phig(ix,iy)+pgrads(ix,iy,2)
      enddo
      enddo
c
c evaluate the l.h.s. from phi value of previous iteration
c and compute the residual of equation
c
c      allocate(dphids(0:nx1,0:ny1,0:ns1))
c      allocate(sdphds(0:nx1,0:ny1,0:ns1))
      do is=1,ns
      do iy=1,ny1
      do ix=1,nx1
        dphids(ix,iy,is)=(phi(ix,iy,is)-phi(ix,iy,is-1))/ds0(is)
      enddo
      enddo
      enddo
c
      do is=1,ns-1
      do iy=1,ny1
      do ix=1,nx1
        sdphds(ix,iy,is)=s0ds2a(is)*(phi(ix,iy,is+1)-phi(ix,iy,is-1))
      enddo
      enddo
      enddo
c
      do is=1,ns-1
      do iy=2,ny
      do ix=2,nx
      dphi(ix,iy,is)=f(ix,iy,is)-(
     :  (phi(ix+1,iy,is)+phi(ix-1,iy,is)-2.*phi(ix,iy,is))/dxx
     :  +(phi(ix,iy+1,is)+phi(ix,iy-1,is)-2.*phi(ix,iy,is))/dyy
     :  +(s(ix,iy,is+1)*s(ix,iy,is+1)*dphids(ix,iy,is+1)
     :  -s(ix,iy,is)*s(ix,iy,is)*dphids(ix,iy,is))/ds1(is)
c lhs4(ix,iy,is)=
     :  +ppgra2(ix,iy)*s0ds1(is)
     :  *(sigma1(is+1)*dphids(ix,iy,is+1)-sigma1(is)*dphids(ix,iy,is))
c lhs5(ix,iy,is)=
     :  +gpplap(ix,iy)*sdphds(ix,iy,is)
c lhs6(ix,iy,is)=
     :  -2.*ppdx(ix,iy,3)*(sdphds(ix+1,iy,is)-sdphds(ix-1,iy,is))/dx2
     :  -2.*ppdy(ix,iy,3)*(sdphds(ix,iy+1,is)-sdphds(ix,iy-1,is))/dy2
c dphi(ix,iy,is)=f(ix,iy,is)-(lhs1(ix,iy,is)+lhs2(ix,iy,is)
c +lhs3(ix,iy,is)+lhs4(ix,iy,is)+lhs5(ix,iy,is)+lhs6(ix,iy,is)
     :  )
      enddo
      enddo
      enddo
c
c solve the standard poisson equation with homogeneous b.c.'s
c
      if(eldout.and.prt) then
         tk=0.
         call wri3ar(dphi,2,nx,2,ny,1,ns-1,'dphi0   ',tk,iomod,nx,ny)
      endif
c      deallocate(dphids)
c      deallocate(sdphds)

c      allocate(work2(0:nx1,0:ny1,0:ns1))
c      allocate(work3(0:nx1,0:ny1,0:ns1))
c
c       write(*,*) 'call pos3'
      call pos3(dphi,wrk5,wrk3,wrk4)
c      deallocate(work2)
c      deallocate(work3)
c
c update the phi value after one iteration and exam the convergence
c
      do is=0,ns
      do iy=1,ny1
      do ix=1,nx1
        phi(ix,iy,is)=phi(ix,iy,is)+dphi(ix,iy,is)
      enddo
      enddo
      enddo
c
      phisum=abssuv(phi,0,nx1,0,ny1,0,ns)
      dpsum=abssuv(dphi,0,nx1,0,ny1,0,ns)
c
      if(eldout.and.prt .or. echeck) then
        write(nchn1,'(''Ellipt:'',2i5,2e15.7)') nstep,ite,phisum,dpsum
        tk=0.
        call wri3ar(dphi,1,nx1,1,ny1,0,ns,'dphi    ',tk,iomod,nx,ny)
      endif
c
      if(phisum .eq. 0.) phisum=1.
      if(dpsum/phisum .le. eps) then
c        deallocate(f)
c        deallocate(dphi)
c        deallocate(work1)
c        deallocate(ppgra2)
c        deallocate(pplap)
c        deallocate(gpplap)
         if(verbose.ge.3) write(*,*) 'ellipt iter=',ite
        return
      endif
1000  continue
c
      if(echeck) then
        write(nchn1,2001) nstep,itma
        tk=0.
        call wri3ar(f,2,nx,2,ny,1,ns-1,'f       ',tk,iomod,nx,ny)
        call wrigar(pgradx,1,2,0,ny1,0,ns,1,2,1,ny1,0,ns
     :    ,'pgradx  ',tk,3)
        call wrigar(pgrady,0,nx1,1,2,0,ns,1,nx1,1,2,0,ns
     :    ,'pgrady  ',tk,2)
        call wrigar(pgrads,0,nx1,0,ny1,1,2,2,nx,2,ny,1,2
     :    ,'pgrads  ',tk,1)
        crash='crash.out'
        crash2='crsh'
        call out6(iomod,2,crash,crash2)
      endif

      write(0,*) 'nh3dWarning=Did not converge in Ellipt after ',itma,
     :  'Rel Error=',dpsum/phisum

      return
2001  format(' ellipt did not converge at nstep=',i5,' after ',i3,
     :   ' iterations')
      end      