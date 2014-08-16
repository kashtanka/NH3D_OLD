      subroutine moment(um2,vm2,wm2,ptm2,wrk1
     :  ,uflux,vflux,wflux)
c-----------------------------------------------------------------------
c time integration of momentum equations
c-----------------------------------------------------------------------
c
      use alloc

      implicit real*8(a-h,o-z)

      dimension um2(0:nx1,0:ny1,0:ns1)
      dimension vm2(0:nx1,0:ny1,0:ns1)
      dimension wm2(0:nx1,0:ny1,0:ns1)
      dimension ptm2(0:nx1,0:ny1,0:ns1)
      dimension wrk1(0:nx1,0:ny1,0:ns1)
      dimension uflux(0:nx1,0:ny1,0:ns1)
      dimension vflux(0:nx1,0:ny1,0:ns1)
      dimension wflux(0:nx1,0:ny1,0:ns1)
      

      logical sim,nao
      parameter (sim=.true.,nao=.false.)
c
      do iy=0,ny1
      do ix=0,nx1
        ppdx10(ix,iy,2)=(pp(ix+1,iy,2)-pp(ix,iy,2))/(dx*pp10(ix,iy,2))
        ppdy01(ix,iy,2)=(pp(ix,iy+1,2)-pp(ix,iy,2))/(dy*pp01(ix,iy,2))
      enddo
      enddo
c
c x equation:
c
c      allocate(um2(0:nx1,0:ny1,0:ns1))
c      allocate(vm2(0:nx1,0:ny1,0:ns1))
c      allocate(wm2(0:nx1,0:ny1,0:ns1))
c      allocate(ptm2(0:nx1,0:ny1,0:ns1))
      if(verbose.ge.3) write(*,*) 'moment:0'
      do is=0,ns
      do iy=0,ny1
      do ix=0,nx1
        um2(ix,iy,is)=u(ix,iy,is,2)
        u(ix,iy,is,2)=u(ix,iy,is,3)
        vm2(ix,iy,is)=v(ix,iy,is,2)
        v(ix,iy,is,2)=v(ix,iy,is,3)
        wm2(ix,iy,is)=w(ix,iy,is,2)
        w(ix,iy,is,2)=w(ix,iy,is,3)
        ptm2(ix,iy,is)=pt(ix,iy,is,2)
        pt(ix,iy,is,2)=pt(ix,iy,is,3)
      enddo
      enddo
      enddo
     
c
      if(verbose.ge.3) write(*,*) 'moment:u'
      do is=1,ns-1
      do iy=2,ny
      do ix=2,nx-1
        fx=-(phi(ix+1,iy,is)-phi(ix,iy,is))/dx
     :    +s0ds4a(is)*ppdx10(ix,iy,2)*(phi(ix+1,iy,is+1)
     :    +phi(ix,iy,is+1)-phi(ix+1,iy,is-1)-phi(ix,iy,is-1))
        u(ix,iy,is,3)=(pp10(ix,iy,1)*um2(ix,iy,is)-dtl*(uflux(ix,iy,is)
     :    -difunu(ix,iy,is)-fx*pp10(ix,iy,2)))/pp10(ix,iy,3)
        if(fcori.ne.0.) then
         u(ix,iy,is,3)=u(ix,iy,is,3)+dtl*fcori*(v(ix,iy,is,3)-vgeos(is))
        endif
      enddo
      enddo
      enddo
      
c      deallocate(uflux)
c
c
c y equations:
c
      if(verbose.ge.3) write(*,*) 'moment:v'
      do is=1,ns-1
      do iy=2,ny-1
      do ix=2,nx
        fy=-(phi(ix,iy+1,is)-phi(ix,iy,is))/dy
     :    +s0ds4a(is)*ppdy01(ix,iy,2)*(phi(ix,iy+1,is+1)
     :    +phi(ix,iy,is+1)-phi(ix,iy+1,is-1)-phi(ix,iy,is-1))        
        v(ix,iy,is,3)=(pp01(ix,iy,1)*vm2(ix,iy,is)-dtl*(vflux(ix,iy,is)
     :    -difunv(ix,iy,is)-fy*pp01(ix,iy,2)))/pp01(ix,iy,3)
        if(fcori.ne.0.) then
         v(ix,iy,is,3)=v(ix,iy,is,3)-dtl*fcori*(u(ix,iy,is,3)-ugeos(is))
        endif
      enddo
      enddo
      enddo
      
c      deallocate(vflux)
c
c
c     call wri3ar(v(0,0,0,3),1,nx1,0,ny1,1,ns-1,'v-mom   ',tk,0,nx,ny)
c
c horizontal boundary conditions have been previously defined
c (in the previous time step)
c
      if(verbose.ge.3) write(*,*) 'moment:uvbc'
      call uvbc(nx,ny,ns
     :  ,u(0,0,0,3),ubx,ubx2,uby,uby2,ucc,ucc2,iobuux,iobuuy
     :  ,v(0,0,0,3),vbx,vbx2,vby,vby2,vcc,vcc2,iobvvx,iobvvy)
     
 !     do ix=0,nx+1
 !     do is=1,ns-1
 !     v(ix,ny+1,is,3)=vs(ix,ny+1,is)
 !     u(ix,ny+1,is,3)=us(ix,ny+1,is)      
 !     enddo
 !     enddo

c                                              _
c put in surface drag: surface stress=cd*ro*vs*vs
c
      if(ifsoil.ne.0.) then
      if(verbose.ge.3) write(*,*) 'moment:ifsoil.ne.0'
        do iy=2,ny
          do ix=2,nx-1
!            u(ix,iy,ns-1,3)=u(ix,iy,ns-1,3)-dtl*cdm(ix,iy)
!     :        *u(ix,iy,ns-1,2)/deltaz(ix,iy)
          enddo
        enddo
        do iy=2,ny-1
          do ix=2,nx
!            v(ix,iy,ns-1,3)=v(ix,iy,ns-1,3)-dtl*cdm(ix,iy)
!     :        *v(ix,iy,ns-1,2)/deltaz(ix,iy)
          enddo
        enddo
      elseif(cdcoef.gt.0.) then
        cddrag=cdcoef*dtl*0.5
        do iy=2,ny
        do ix=2,nx
          deltaz(ix,iy)=0.5*(phis(ix,iy,ns-1)-phis(ix,iy,ns)
     :      +phi(ix,iy,ns-1)-phi(ix,iy,ns))/g
c         deltaz(ix,iy)=(phis(ix,iy,ns-1)+phi(ix,iy,ns-1))/g-hsuf(ix,iy)
c          write(99,*) 'acol?eltaz(ix,iy)',deltaz(ix,iy), ix,iy
        enddo
        enddo
        do iy=2,ny
        do ix=2,nx
         uvsuf(ix,iy)=dsqrt((0.5*(u(ix,iy,ns-1,2)+u(ix-1,iy,ns-1,2)))**2
     :     +(0.5*(v(ix,iy,ns-1,2)+v(ix,iy-1,ns-1,2)))**2)
        enddo
        enddo
        do iy=2,ny
        do ix=2,nx-1
          u(ix,iy,ns-1,3)=u(ix,iy,ns-1,3)-cddrag
     :      *(uvsuf(ix,iy)+uvsuf(ix+1,iy))*u(ix,iy,ns-1,2)/deltaz(ix,iy)
        enddo
        enddo
        do iy=2,ny-1
        do ix=2,nx
          v(ix,iy,ns-1,3)=v(ix,iy,ns-1,3)-cddrag
     :      *(uvsuf(ix,iy)+uvsuf(ix,iy+1))*v(ix,iy,ns-1,2)/deltaz(ix,iy)
        enddo
        enddo
      endif
c
      if(raylei) then
      if(verbose.ge.3) write(*,*) 'moment:raylei'
         do is=1,idrmax
         do iy=0,ny1
         do ix=0,nx1
           u(ix,iy,is,3)=u(ix,iy,is,3)-(u(ix,iy,is,2)-us(ix,iy,is))
     :       *dtl/taudra(is)
           v(ix,iy,is,3)=v(ix,iy,is,3)-(v(ix,iy,is,2)-vs(ix,iy,is))
     :       *dtl/taudra(is)
         enddo
         enddo
         enddo
      endif

c
c Lateral absorption column
c
      if(nxsponge2.gt.0.or.nxsponge1.gt.0) then
      if(verbose.ge.3) write(*,*) 'moment:sponge'
        do is=0,ns
          do iy=0,ny1
            do ix=0,nxsponge2
              u(ix,iy,is,3)=u(ix,iy,is,3)-(u(ix,iy,is,2)-us(ix,iy,is))
     :          *dtl/tauspo(ix)
              v(ix,iy,is,3)=v(ix,iy,is,3)-(v(ix,iy,is,2)-vs(ix,iy,is))
     :          *dtl/tauspo(ix)
            enddo
            do ix=nx1-nxsponge1,nx1
              u(ix,iy,is,3)=u(ix,iy,is,3)-(u(ix,iy,is,2)-us(ix,iy,is))
     :          *dtl/tauspo(ix)
              v(ix,iy,is,3)=v(ix,iy,is,3)-(v(ix,iy,is,2)-vs(ix,iy,is))
     :          *dtl/tauspo(ix)
            enddo
          enddo
        enddo
      endif
c
c Wind - Lateral forcing column - ymin and y max (lateral sponge)
c
      if(nysponge2.gt.0.or.nysponge1.gt.0) then
        do is=0,ns
          do ix=0,nx1
            do iy=0,nysponge2
              u(ix,iy,is,3)=u(ix,iy,is,3)-(u(ix,iy,is,2)-us(ix,iy,is))
     :          *dtl/tauspoy(iy)
              v(ix,iy,is,3)=v(ix,iy,is,3)-(v(ix,iy,is,2)-vs(ix,iy,is))
     :          *dtl/tauspoy(iy)
            enddo
            do iy=ny1-nysponge1,ny1
              u(ix,iy,is,3)=u(ix,iy,is,3)-(u(ix,iy,is,2)-us(ix,iy,is))
     :          *dtl/tauspoy(iy)
              v(ix,iy,is,3)=v(ix,iy,is,3)-(v(ix,iy,is,2)-vs(ix,iy,is))
     :          *dtl/tauspoy(iy)
            enddo
          enddo
        enddo
      endif
c
      if(dohsmo .and. mod(nstep,numsmo).eq.0) then
      if(verbose.ge.3) write(*,*) 'moment:hsmop'
         call hsmop2(u(0,0,0,3),ufd,us,hk1,nx1,ny1,ns,1,nx,1,ny1,1,ns-1
     :      ,hdamp,perdam,dt)
         call hsmop2(v(0,0,0,3),vfd,vs,hk1,nx1,ny1,ns,1,nx1,1,ny,1,ns-1
     :      ,hdamp,perdam,dt)
      endif
      if(dovsmo .and. mod(nstep,numsmo).eq.0) then
c         allocate(work1(0:nx1,0:ny1,0:ns1))
         call vsmop(u(0,0,0,3),us,wrk1,nx1,ny1,ns,1,nx,1,ny1,0,ns
     :      ,vdamp,perdam)
         call vsmop(v(0,0,0,3),vs,wrk1,nx1,ny1,ns,1,nx1,1,ny,0,ns
     :      ,vdamp,perdam)
c         deallocate(work1)
      endif
c
c vertical boundary conditions
c (note that this cannot be done before uvbc)
c
      do iy=0,ny1
      do ix=0,nx1
        u(ix,iy,0,3)=u(ix,iy,1,3)
        u(ix,iy,ns,3)=u(ix,iy,ns-1,3)
        v(ix,iy,0,3)=v(ix,iy,1,3)
        v(ix,iy,ns,3)=v(ix,iy,ns-1,3)
      enddo
      enddo
c
      call aselin(u(0,0,0,3),u(0,0,0,2),um2
     :  ,0,nx1,0,ny1,0,ns,filt)
      call aselin(v(0,0,0,3),v(0,0,0,2),vm2
     :  ,0,nx1,0,ny1,0,ns,filt)
c
      if(verbose.ge.3) write(*,*) 'moment:raduv'
      call raduv(nx,ny,ns
     :  ,u(0,0,0,3),u(0,0,0,2),um2
     :  ,ubx,ubx2,uby,uby2,ucc,ucc2,ds1y,yl,pp10(0,0,4)
     :  ,v(0,0,0,3),v(0,0,0,2),vm2
     :  ,vbx,vbx2,vby,vby2,vcc,vcc2,ds1x,xl,pp01(0,0,4)
     :  ,ttmin,masout,uubout,prt,nchn1,ds0,ds1,s,nao
     :  ,xlcor,xrcor,ylcor,yrcor,ifluxcor,dpmeddt,dx,dy)
c      deallocate(um2)
c      deallocate(vm2)
c
c vertical equations (sigma):
c
      if(verbose.ge.3) write(*,*) 'moment:thetav'
      if(qif.eq.0.) then
        do is=0,ns
          do iy=1,ny1
            do ix=1,nx1
              thetav(ix,iy,is)=g*pt(ix,iy,is,2)/pts(ix,iy,is)
            enddo
          enddo
        enddo
      elseif(ifqc.eq.0) then
        do is=0,ns
          do iy=1,ny1
            do ix=1,nx1
              thetav(ix,iy,is)=g*(pt(ix,iy,is,2)/pts(ix,iy,is)
     :          +qdrag*0.61*(qv(ix,iy,is,2)-qvs(ix,iy,is)))
            enddo
          enddo
        enddo
      elseif(ifqr.eq.0) then
        do is=0,ns
          do iy=1,ny1
            do ix=1,nx1
              thetav(ix,iy,is)=g*(pt(ix,iy,is,2)/pts(ix,iy,is)
     :          +qdrag*0.61*(qv(ix,iy,is,2)-qvs(ix,iy,is)
     :          -qc(ix,iy,is,2)))
            enddo
          enddo
        enddo
      elseif(ifqi.eq.0) then                                                    DC,11.2009
        do is=0,ns
          do iy=1,ny1
            do ix=1,nx1
              thetav(ix,iy,is)=g*(pt(ix,iy,is,2)/pts(ix,iy,is)
     :          +qdrag*0.61*(qv(ix,iy,is,2)-qvs(ix,iy,is)
     :          -qc(ix,iy,is,2)-qr(ix,iy,is,2)))
            enddo
          enddo
        enddo
	else                                                                      DC,11.2009
	  do is=0,ns                                                              DC,11.2009
          do iy=1,ny1                                                           DC,11.2009
            do ix=1,nx1                                                         DC,11.2009
              thetav(ix,iy,is)=g*(pt(ix,iy,is,2)/pts(ix,iy,is)                  DC,11.2009
     :          +qdrag*0.61*(qv(ix,iy,is,2)-qvs(ix,iy,is)                       DC,11.2009
     :          -qc(ix,iy,is,2)-qr(ix,iy,is,2)-qci(ix,iy,is,2)                  DC,11.2009
     :          -qsn(ix,iy,is,2)))                                              DC,11.2009
            enddo                                                               DC,11.2009
          enddo                                                                 DC,11.2009
        enddo                                                                   DC,11.2009
      endif
      do is=ns,1,-1
        do iy=1,ny1
          do ix=1,nx1
            thetav(ix,iy,is)=0.5*(thetav(ix,iy,is)+thetav(ix,iy,is-1))
          enddo
        enddo
      enddo

      if(verbose.ge.3) write(*,*) 'moment:w'
      do is=1,ns-1
      do iy=1,ny1
      do ix=1,nx1
        fs=s(ix,iy,is)*(phi(ix,iy,is)-phi(ix,iy,is-1))/ds0(is)
     :    +thetav(ix,iy,is)

        w(ix,iy,is,3)=(wm2(ix,iy,is)*pp(ix,iy,1)-dtl*(wflux(ix,iy,is)
     :    -difunw(ix,iy,is)-pp(ix,iy,2)*fs))/pp(ix,iy,3)
      enddo
      enddo
      enddo
c      deallocate(wflux)
c
      if(raylei) then
         do is=1,idrmax
         do iy=1,ny1
         do ix=1,nx1
           w(ix,iy,is,3)=w(ix,iy,is,3)*(1.-dtl/(taudra(is)))
         enddo
         enddo
         enddo
      endif

      if(nxsponge2.gt.0.or.nxsponge1.gt.0) then
        do is=0,ns
          do iy=0,ny1
            do ix=0,nxsponge2
              w(ix,iy,is,3)=w(ix,iy,is,3)*(1.-dtl/tauspo(ix))
            enddo
            do ix=nx1-nxsponge1,nx1
              w(ix,iy,is,3)=w(ix,iy,is,3)*(1.-dtl/tauspo(ix))
            enddo
          enddo
        enddo
      endif
c
      if(nxsponge2.gt.0.or.nxsponge1.gt.0) then
        do is=0,ns
          do ix=0,nx1
            do iy=0,nysponge2
              w(ix,iy,is,3)=w(ix,iy,is,3)*(1.-dtl/tauspoy(iy))
            enddo
            do iy=ny1-nysponge1,ny1
              w(ix,iy,is,3)=w(ix,iy,is,3)*(1.-dtl/tauspoy(iy))
            enddo
          enddo
        enddo
      endif

!        do is=0,ns
!          do ix=0,nx1
!            do iy=0,10
!              w(ix,iy,is,3)=w(ix,iy,is,3)*(1.-dtl/1000.)
!            enddo
!          enddo
!         enddo

c
      if(dohsmo .and. mod(nstep,numsmo).eq.0) then
         call hsmoot2(w(0,0,0,3),wfd,hk1,nx1,
     :    ny1,ns,2,nx,2,ny,1,ns-1,hdamp,dt)
      endif
      if(dovsmo .and. mod(nstep,numsmo).eq.0) then
  !       call vsmoot(w(0,0,0,3),wrk1,nx1,ny1,ns,2,nx,2,ny,1,ns-1,vdamp)
      endif
c
c boundary conditions for w can only be calculated after both
c sufpp and refs, since they use dpp.
c
      return
      end
      
      subroutine moment1(ww1,ww2,ww3,ptm2,ww5
     :  ,uflux,vflux,wflux)
c-----------------------------------------------------------------------
c time integration of momentum equations in gsp--balanced model
c-----------------------------------------------------------------------
c
c this is a
c modification of moment, assumes that adjust  is true
c
      use alloc
      implicit real*8(a-h,o-z)
      dimension ptm2(0:nx1,0:ny1,0:ns1)
      dimension ww1(0:nx1,0:ny1,0:ns1)
      dimension ww2(0:nx1,0:ny1,0:ns1)
      dimension ww3(0:nx1,0:ny1,0:ns1)
      dimension ww5(0:nx1,0:ny1,0:ns1)
      dimension uflux(0:nx1,0:ny1,0:ns1)
      dimension vflux(0:nx1,0:ny1,0:ns1)
      dimension wflux(0:nx1,0:ny1,0:ns1)
c
      logical sim,nao
      parameter (sim=.true.,nao=.false.)
c
c
c x equation:
c
      do is=0,ns
      do iy=0,ny1
      do ix=0,nx1
        ww1(ix,iy,is)=u(ix,iy,is,2)
        u(ix,iy,is,2)=u(ix,iy,is,3)
        ww2(ix,iy,is)=v(ix,iy,is,2)
        v(ix,iy,is,2)=v(ix,iy,is,3)
        ptm2(ix,iy,is)=pt(ix,iy,is,2)
        pt(ix,iy,is,2)=pt(ix,iy,is,3)
      enddo
      enddo
      enddo
c
      do is=1,ns-1
      do iy=2,ny
      do ix=2,nx-1
        fx=-(phi(ix+1,iy,is)-phi(ix,iy,is))/dx
     :    +s0ds4a(is)*ppdx10(ix,iy,2)*(phi(ix+1,iy,is+1)
     :    +phi(ix,iy,is+1)-phi(ix+1,iy,is-1)-phi(ix,iy,is-1))
        u(ix,iy,is,3)=(pp10(ix,iy,1)*ww1(ix,iy,is)-dtl*(uflux(ix,iy,is)
     :    -difunu(ix,iy,is)-fx*pp10(ix,iy,2)))/pp10(ix,iy,3)
        if (fcori.ne.0.) then
         u(ix,iy,is,3)=u(ix,iy,is,3)+dtl*fcori*(v(ix,iy,is,3)-vgeos(is))
        endif
      enddo
      enddo
      enddo
c
c y equations:
c
      do is=1,ns-1
      do iy=2,ny-1
      do ix=2,nx
        fy=-(phi(ix,iy+1,is)-phi(ix,iy,is))/dy
     :    +s0ds4a(is)*ppdy01(ix,iy,2)*(phi(ix,iy+1,is+1)
     :    +phi(ix,iy,is+1)-phi(ix,iy+1,is-1)-phi(ix,iy,is-1))
        v(ix,iy,is,3)=(pp01(ix,iy,1)*ww2(ix,iy,is)-dtl*(vflux(ix,iy,is)
     :    -difunv(ix,iy,is)-fy*pp01(ix,iy,2)))/pp01(ix,iy,3)
        if (fcori.ne.0.) then
         v(ix,iy,is,3)=v(ix,iy,is,3)-dtl*fcori*(u(ix,iy,is,3)-ugeos(is))
        endif
      enddo
      enddo
      enddo
c
c horizontal boundary conditions have been previously defined
c (in the previous time step)
c
      call uvbc(nx,ny,ns
     :  ,u(0,0,0,3),ubx,ubx2,uby,uby2,ucc,ucc2,iobuux,iobuuy
     :  ,v(0,0,0,3),vbx,vbx2,vby,vby2,vcc,vcc2,iobvvx,iobvvy)
c                                              _
c put in surface drag: surface stress=cd*ro*vs*vs
c
      if(ifsoil.ne.0) then
        do iy=2,ny
          do ix=2,nx-1
!            u(ix,iy,ns-1,3)=u(ix,iy,ns-1,3)-dtl*cdm(ix,iy)
!     :        *u(ix,iy,ns-1,2)/deltaz(ix,iy)
          enddo
        enddo
        do iy=2,ny-1
          do ix=2,nx
 !           v(ix,iy,ns-1,3)=v(ix,iy,ns-1,3)-dtl*cdm(ix,iy)
 !    :        *v(ix,iy,ns-1,2)/deltaz(ix,iy)
          enddo
        enddo
      elseif(cdcoef.gt.0.) then
        cddrag=cdcoef*dtl*0.5
        do iy=2,ny
        do ix=2,nx
c        deltaz(ix,iy)=(phis(ix,iy,ns-1)+phi(ix,iy,ns-1))/g-hsuf(ix,iy)
          deltaz(ix,iy)=0.5*(phis(ix,iy,ns-1)-phis(ix,iy,ns)
     :      +phi(ix,iy,ns-1)-phi(ix,iy,ns))/g
c          write(99,*) 'maisuldeltaz',deltaz(ix,iy)
        enddo
        enddo
        do iy=2,ny
        do ix=2,nx
          uvsuf(ix,iy)=sqrt((0.5*(u(ix,iy,ns-1,2)+u(ix-1,iy,ns-1,2)))**2
     :      +(0.5*(v(ix,iy,ns-1,2)+v(ix,iy-1,ns-1,2)))**2)
        enddo
        enddo
        do iy=2,ny
        do ix=2,nx-1
          u(ix,iy,ns-1,3)=u(ix,iy,ns-1,3)-cddrag
     :      *(uvsuf(ix,iy)+uvsuf(ix+1,iy))*u(ix,iy,ns-1,2)/deltaz(ix,iy)
        enddo
        enddo
        do iy=2,ny-1
        do ix=2,nx
          v(ix,iy,ns-1,3)=v(ix,iy,ns-1,3)-cddrag
     :      *(uvsuf(ix,iy)+uvsuf(ix,iy+1))*v(ix,iy,ns-1,2)/deltaz(ix,iy)
        enddo
        enddo
      endif
c
      if(raylei) then
        do is=1,idrmax
        do iy=0,ny1
        do ix=0,nx1
          u(ix,iy,is,3)=u(ix,iy,is,3)-(u(ix,iy,is,2)-us(ix,iy,is))
     :      *dtl/taudra(is)
          v(ix,iy,is,3)=v(ix,iy,is,3)-(v(ix,iy,is,2)-vs(ix,iy,is))
     :      *dtl/taudra(is)
        enddo
        enddo
        enddo
      endif

c
c Lateral absorption column
c
      if(nxsponge2.gt.0.or.nxsponge1.gt.0) then
        do is=0,ns
          do iy=0,ny1
            do ix=0,nxsponge2
              u(ix,iy,is,3)=u(ix,iy,is,3)-(u(ix,iy,is,2)-us(ix,iy,is))
     :          *dtl/tauspo(ix)
              v(ix,iy,is,3)=v(ix,iy,is,3)-(v(ix,iy,is,2)-vs(ix,iy,is))
     :          *dtl/tauspo(ix)
            enddo
            do ix=nx1-nxsponge1,nx1
              u(ix,iy,is,3)=u(ix,iy,is,3)-(u(ix,iy,is,2)-us(ix,iy,is))
     :          *dtl/tauspo(ix)
              v(ix,iy,is,3)=v(ix,iy,is,3)-(v(ix,iy,is,2)-vs(ix,iy,is))
     :          *dtl/tauspo(ix)
            enddo
          enddo
        enddo
      endif

      if(nysponge2.gt.0.or.nysponge1.gt.0) then
        do is=0,ns
          do ix=0,nx1
            do iy=0,nysponge2
              u(ix,iy,is,3)=u(ix,iy,is,3)-(u(ix,iy,is,2)-us(ix,iy,is))
     :          *dtl/tauspoy(iy)
              v(ix,iy,is,3)=v(ix,iy,is,3)-(v(ix,iy,is,2)-vs(ix,iy,is))
     :          *dtl/tauspoy(iy)
            enddo
            do iy=ny1-nysponge1,ny1
              u(ix,iy,is,3)=u(ix,iy,is,3)-(u(ix,iy,is,2)-us(ix,iy,is))
     :          *dtl/tauspoy(iy)
              v(ix,iy,is,3)=v(ix,iy,is,3)-(v(ix,iy,is,2)-vs(ix,iy,is))
     :          *dtl/tauspoy(iy)
            enddo
          enddo
        enddo
      endif
c
      if(dohsmo .and. mod(nstep,numsmo).eq.0) then

        call hsmop(u(0,0,0,3),us,hk1,nx1,ny1,ns,1,nx,1,ny1,1,ns-1
     :    ,hdamp,perdam)
        call hsmop(v(0,0,0,3),vs,hk1,nx1,ny1,ns,1,nx1,1,ny,1,ns-1
     :    ,hdamp,perdam)
      endif
      if(dovsmo .and. mod(nstep,numsmo).eq.0) then
        call vsmop(u(0,0,0,3),us,ww5,nx1,ny1,ns,1,nx,1,ny1,0,ns
     :    ,vdamp,perdam)
        call vsmop(v(0,0,0,3),vs,ww5,nx1,ny1,ns,1,nx1,1,ny,0,ns
     :    ,vdamp,perdam)
      endif
c
      call aselin(u(0,0,0,3),u(0,0,0,2),ww1
     :  ,0,nx1,0,ny1,0,ns,filt)
      call aselin(v(0,0,0,3),v(0,0,0,2),ww2
     :  ,0,nx1,0,ny1,0,ns,filt)
c
c computation of vertical velocity and mass balance restoration
c      call nduvw(1,ww3,ww5)
      call nduvw(1)
c -------------------------------------
c
      call raduv(nx,ny,ns
     :  ,u(0,0,0,3),u(0,0,0,2),ww1
     :  ,ubx,ubx2,uby,uby2,ucc,ucc2,ds1y,yl,pp10(0,0,4)
     :  ,v(0,0,0,3),v(0,0,0,2),ww2
     :  ,vbx,vbx2,vby,vby2,vcc,vcc2,ds1x,xl,pp01(0,0,4)
     :  ,ttmin,masout,uubout,prt,nchn1,ds0,ds1,s,nao
     :  ,xlcor,xrcor,ylcor,yrcor,ifluxcor,dpmeddt,dx,dy)
c
      return
      end
