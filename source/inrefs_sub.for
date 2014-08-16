      subroutine inrefs

c-----------------------------------------------------------------------
c specification of reference state
c new version: if fcor.ne.0 the initial state is defined using
c geostrophic balance.
c-----------------------------------------------------------------------
c
      use alloc
      use refstate

      implicit real*8(a-h,o-z)

c
c resmax is the maximum relative error admited in all tems,pts and phis.
c
      parameter(nnomax=100,nitmax=1)
c
c adjust pp,tems with  phis in requirement of hydrostatic
c balance in the vertical.
c if fcor.ne.0 impose geostrophic balance for initial fields.
c
c note only profile 1 is used for the adjustment (nhad0520)
c
      allocate(errp(0:nx+1,0:ny+1))
      allocate(idath(0:nx+1,0:ny+1))
c
c initialize reference profiles:
c
      call refprofs
c     call refprofs(zpro,npro,dzpro,fpro,ndat,nprof,fnpro
c    :  ,iorefu,usdat,zusdat,uspro0,dusdz
c    :  ,iorefv,vsdat,zvsdat,vspro0,dvsdz
c    :  ,ioreft,thdat,zthdat,thpro0,dthdz,ptdat,pts0
c    :  ,iorefq,qvsdat,zqvsdat,qvspro0,dqvsdz
c    :  ,psdat,pa,verbose,nchn1)

      ip=1
c
      do id=ndat,1,-1
        do iy=0,ny+1
          do ix=0,nx+1
            if(hsuf(ix,iy).le.zthdat(id,ip)) then
              idath(ix,iy)=id
            endif
          enddo
        enddo
      enddo
c
      if(ioreft.eq.2) then
        do iy=0,ny1
        do ix=0,nx1
          id=idath(ix,iy)
c        tvgrad=(thdat(id,ip)-thdat(id-1,ip))/(zthdat(id,ip)-zthdat(id-1,ip))
c        psuf(ix,iy)=psdat(id-1,ip)*((1.-tvgrad*hsuf(ix,iy))
c           /(thdat(id-1,ip)+tvgrad*zthdat(id-1,ip)))**(g/(r*tvgrad))
          tsuf(ix,iy)=
     :      (thdat(id-1,ip)+(thdat(id,ip)-thdat(id-1,ip))*(hsuf(ix,iy)
     :      -zthdat(id-1,ip))/(zthdat(id,ip)-zthdat(id-1,ip)))
     :      *(1.+0.61*(qvsdat(id-1,ip)+(qvsdat(id,ip)-qvsdat(id-1,ip))
     :      *(hsuf(ix,iy)
     :      -zqvsdat(id-1,ip))/(zqvsdat(id,ip)-zqvsdat(id-1,ip))))
          tmed=0.5*(thdat(id-1,ip)*(1.+0.61*qvsdat(id-1,ip))
     :      +tsuf(ix,iy))
c         write(*,*) ix,iy,id,hsuf(ix,iy),tmed
c        ,psdat(id-1,ip),zthdat(id-1,ip)
          psuf(ix,iy)=psdat(id-1,ip)*exp(-g/(r*tmed)
     :      *(hsuf(ix,iy)-zthdat(id-1,ip)))
c         write(*,*) psuf(ix,iy)
        enddo
        enddo
      else
        do iy=0,ny1
        do ix=0,nx1
          id=idath(ix,iy)  !level of *dat arrays where the orography surface is, 1 in case of CAOs
          
          ptsuf=
     :      (ptdat(id-1,ip)+(ptdat(id,ip)-ptdat(id-1,ip))*(hsuf(ix,iy)      ! humidity is taken into account
     :      -zthdat(id-1,ip))/(zthdat(id,ip)-zthdat(id-1,ip)))              ! even if qif=0 !!! 
     :      *(1.+0.61*(qvsdat(id-1,ip)+(qvsdat(id,ip)-qvsdat(id-1,ip))  
     :      *(hsuf(ix,iy)
     :      -zqvsdat(id-1,ip))/(zqvsdat(id,ip)-zqvsdat(id-1,ip))))
           ptmed=0.5*(ptdat(id-1,ip)*(1.+0.61*qvsdat(id-1,ip))+ptsuf)
           psuf(ix,iy)=(psdat(id-1,ip)**akapa-g/cp*p00**akapa
     :     *(hsuf(ix,iy)-zthdat(id-1,ip))/ptmed)**(1./akapa)
           tsuf(ix,iy)=ptsuf*(psuf(ix,iy)/p00)**akapa
        enddo
        enddo
      endif

c        tk=1.e5
c        call wri2ar(psuf,0,nx1,0,ny1,'psuf-ini',tk,nx,ny)

      do is=0,ns+1
      do iy=0,ny1
      do ix=0,nx1
        tems(ix,iy,is)=tsuf(ix,iy)                  !ogo, ves massiv raven prisemnoy T
      enddo
      enddo
      enddo
c
      do iy=0,ny1
      do ix=0,nx1
        tsufi(ix,iy)=tsuf(ix,iy)
        psufi(ix,iy)=psuf(ix,iy)
      enddo
      enddo
c
c      delp=0.
      if(fcor.ne.0.) then
      !write(0,*)'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'
c         do 1777 iy=0,ny1
c            delp=p00*((1.-(fcor*ug*dy*(iy-ny1/2.))/
c     :         (cp*ptsuf))**(cp/r)-1.)
c            do 1776 ix=0,nx1
c               psuf(ix,iy)=psuf(ix,iy)+delp
c1776        continue
c1777     continue
c
        ixref=nx1/2 !3    
        iyref=ny1/2
        pref=psuf(ixref,iyref) 
        do iy=1,ny1
          tmed=(tsuf(0,iy-1)+tsuf(0,iy))/2.
          delpy=-psuf(0,iy-1)*fcor*ug*dy/(tmed*r+fcor*ug*dy/2.)
          psuf(0,iy)=psuf(0,iy-1)+delpy
        enddo
        do ix=1,nx1
          do iy=0,ny1
            tmed=(tsuf(ix-1,iy)+tsuf(ix,iy))/2.
            delpx=+psuf(ix-1,iy)*fcor*vg*dx/(tmed*r-fcor*vg*dx/2.)
            psuf(ix,iy)=psuf(ix-1,iy)+delpx
          enddo
        enddo
        pcor=pref-psuf(ixref,iyref)
        do ix=0,nx1
          do iy=0,ny1
            psuf(ix,iy)=psuf(ix,iy)+pcor
          enddo
        enddo
      endif
c
c identify layer in the reference profile to interpolate the
c reference pressure to the surface: (layer id-1/id)
c
c iteration for reference state variables
c
      if(verbose.ge.0) write(*,1000) nnomax,resmax
      if(verbose.ge.0) write(*,*)'nno,restem,respts,resphs,respsu,respt'
      if(iniout) write(nchn1,1000) nnomax,resmax
c
      allocate(work2(0:nx1,0:ny1,0:ns1))
      allocate(work3(0:nx1,0:ny1,0:ns1))
      allocate(work4(0:nx1,0:ny1,0:ns1))
      allocate(work5(0:nx1,0:ny1,0:ns1))

      do nno=1,nnomax
c
        if (nno.gt.1) then
          do iy=1,ny1
          do ix=1,nx1
            hk1(ix,iy)=psuf(ix,iy)
          enddo
          enddo
        endif
c
c if fcor.ne.0 interpolate surface pressure from reference profile
c assuming constant temperature gradient.
c
        if(fcor.ne.0..and.itepsuf) then
          if(ug.ne.0. .or. vg.ne.0.) then
            if(iophis.eq.3) then
              ixref=nx1/2 !3
              iyref=ny1/2
              pref=psuf(ixref,iyref)
              write(0,*) pref
              call surpsu(ixref,iyref,pref)
            else
              call inpos2(iobphy)
              call surpsd
            endif
          endif
        endif
c
        if(phsout) then
          tk=1.e5
          call wri2ar(psuf,0,nx1,0,ny1,'psuf-k1  ',tk,nx,ny)
        endif

        do iy=1,ny1
        do ix=1,nx1
          pp(ix,iy,3)=psuf(ix,iy)-ptop
        enddo
        enddo
c
        call extrpp(nx2,ny2,pp(0,0,3),0,nx2,0,ny2)
        do iy=0,ny1
        do ix=0,nx1
          dpdx10(ix,iy,3)=(pp(ix+1,iy,3)-pp(ix,iy,3))/dx
          dpdy01(ix,iy,3)=(pp(ix,iy+1,3)-pp(ix,iy,3))/dy
        enddo
        enddo
c
c internal iteration:
c
         do it=1,nitmax

          if(nno.gt.1) then
            work4=tems
            work2=pts
            work5=pt(:,:,:,3)
            work3=phis
c
c interpolate pts in the vertical, according to present values of phis:
c
            call upuvt
          endif
          
c
c integrate hydrostatic equation for reference state.
c
          if(iophis.eq.1) then
            call surhyd
          elseif(iophis.eq.2) then
            call inpos2(iobphy)
            call surphi
          else
            call surphn
            if(verbose.ge.3) write(*,*) 'surphn done'
          endif
c
c calculates phis integrating hydrostatic relation:
c
          do iy=0,ny1
          do ix=0,nx1
            phis(ix,iy,ns)=phisuf(ix,iy)-dsng*r*pp(ix,iy,3)
     :        *(.75*tems(ix,iy,ns)+.25*tems(ix,iy,ns-1))
     :        /(ptop+pp(ix,iy,3)*(1.+0.5*dsng))
          enddo
          enddo
c
          do is=ns-1,0,-1
          do iy=0,ny1
          do ix=0,nx1
            phis(ix,iy,is)=phis(ix,iy,is+1)+ds0(is+1)*r05*pp(ix,iy,3)
     :        *(tems(ix,iy,is+1)+tems(ix,iy,is))
     :        /(ptop+pp(ix,iy,3)*sigma1(is+1))
          enddo
          enddo
          enddo
c
c initialization of arrays for refs and preparation of the solution:
c
          if(inrout) then
            tk=1.e5
            call wri2ar(psuf,1,nx1,1,ny1,'psuf-k  ',tk,nx,ny)
            tk=273.
            call wri2ar(tsuf,1,nx1,1,ny1,'tsuf    ',tk,nx,ny)
           call wri3ar(tems,1,nx1,1,ny1,0,ns,'tems    ',tk,iomod,nx,ny)
            call wri3ar(pts,1,nx1,1,ny1,0,ns,'pts     ',tk,iomod,nx,ny)
c            call wri3ar(pt,1,nx1,1,ny1,0,ns,'pt      ',tk,iomod,nx,ny)
            tk=0.
            call wri3ar(phis,1,nx1,1,ny1,1,ns,'phis    ',tk,iomod,nx,ny)
          endif
c
          if(verbose.ge.3) write(*,*) 'inrefs 1'
          if(nno.gt.1) then
            restem=0.
            respts=0.
            respt=0.
            resphs=0.
            respsu=0.
            do is=0,ns
            do iy=1,ny1
            do ix=1,nx1
              restem=max(restem,abs(work4(ix,iy,is)-tems(ix,iy,is)))
              respts=max(respts,abs(work2(ix,iy,is)-pts(ix,iy,is)))
              respt=max(respt,abs(work5(ix,iy,is)-pt(ix,iy,is,3)))
              resphs=max(resphs,abs(work3(ix,iy,is)-phis(ix,iy,is)))
            enddo
            enddo
            enddo
            do iy=1,ny1
            do ix=1,nx1
              errp(ix,iy)=hk1(ix,iy)-psuf(ix,iy)
              respsu=max(respsu,abs(hk1(ix,iy)-psuf(ix,iy)))
            enddo
            enddo
            if(verbose.ge.0) then
              write(*,1) nno,restem,respts,resphs,respsu,respt
            endif
         if(iniout) write(nchn1,1) nno,restem,respts,resphs,respsu,respt
c            tk=0.
c            call wri2ar(errp,1,nx1,1,ny1,'errp    ',tk,nx,ny)
            if(max(restem/300.,respts/300.,resphs/1.e4,respsu/1.e5
     :        ,respt/300.).lt.resmax) go to 501
          endif
        enddo
      enddo
501   continue
      deallocate(work2,work3,work4,work5)

c      write(*,*) '501'
c
c      call extrpp(nx2,ny2,pp(0,0,3),1,nx1,1,ny1)
c      call extrpp(nx2,ny2,pp(0,0,3),0,nx2,0,ny2)
c
      do k=1,4
      do iy=0,ny2
      do ix=0,nx2
        pp(ix,iy,k)=pp(ix,iy,3)
      enddo
      enddo
      enddo
c
      do iy=0,ny1
      do ix=0,nx1
        psufi(ix,iy)=pp(ix,iy,3)+ptop
        tsufi(ix,iy)=tsuf(ix,iy)
        phisui(ix,iy)=phisuf(ix,iy)
      enddo
      enddo
         tk=1.e5
c         call wri2ar(psufi,1,nx1,1,ny1,'psufi-k ',tk,nx,ny)
c
c      do 220 iy=0,ny2
c      do 220 ix=0,nx2
c      pp(ix,iy,1)=pp(ix,iy,2)
c220   continue
c
      deallocate(errp)
      deallocate(idath)

      entry rinref
c
      do k=1,4
        do iy=0,ny1
        do ix=0,nx1
          pp10(ix,iy,k)=0.5*(pp(ix,iy,k)+pp(ix+1,iy,k))
          pp01(ix,iy,k)=0.5*(pp(ix,iy,k)+pp(ix,iy+1,k))
        enddo
        enddo
c
        do iy=1,ny1
        do ix=1,nx1
          ppdx(ix,iy,k)=(pp(ix+1,iy,k)-pp(ix-1,iy,k))/(dx2*pp(ix,iy,k))
          ppdy(ix,iy,k)=(pp(ix,iy+1,k)-pp(ix,iy-1,k))/(dy2*pp(ix,iy,k))
        enddo
        enddo
      enddo
c
      do iy=0,ny1
      do ix=0,nx1
        dpdx10(ix,iy,3)=(pp(ix+1,iy,3)-pp(ix,iy,3))/dx
        dpdy01(ix,iy,3)=(pp(ix,iy+1,3)-pp(ix,iy,3))/dy
        ppdx10(ix,iy,2)=(pp(ix+1,iy,2)-pp(ix,iy,2))/(dx*pp10(ix,iy,2))
        ppdy01(ix,iy,2)=(pp(ix,iy+1,2)-pp(ix,iy,2))/(dy*pp01(ix,iy,2))
      enddo
      enddo
c
      do iy=0,ny1
      do ix=0,nx1 
        psuf(ix,iy)=pp(ix,iy,3)+ptop
      enddo
      enddo
c
      do iy=1,ny1
      do ix=1,nx1
        tems(ix,iy,ns+1)=tems(ix,iy,ns)
      enddo
      enddo
c
      do is=1,ns
      do iy=1,ny1
      do ix=1,nx1
        s(ix,iy,is)=gr2*(ptop+sigma1(is)*pp(ix,iy,3))/
     :    ((tems(ix,iy,is)+tems(ix,iy,is-1))*pp(ix,iy,3))
      enddo
      enddo
      enddo
c
      do iy=1,ny1
      do ix=1,nx1
        s(ix,iy,0)=gr*(ptop+sigma1(0)*pp(ix,iy,3))/(tems(ix,iy,0)
     :    *pp(ix,iy,3))
        s(ix,iy,ns+1)=gr*(ptop+sigma1(ns1)*pp(ix,iy,3))/(tems(ix,iy,ns)
     :    *pp(ix,iy,3))
      enddo
      enddo
c
      do iy=1,ny1
      do ix=1,nx1
        phig(ix,iy)=g*hsuf(ix,iy)-phisuf(ix,iy)
      enddo
      enddo
c
c initialization of arrays for refs and preparation of the solution:
c
      if(iophis.eq.2) then
         call inpos2(iobphy)
      elseif(iophis.eq.3) then
      endif
c      write(*,*) 'done inrefs'
c
      return
1     format(' #',i5,5e12.3)
1000  format(//' iteration for reference state variables ',i3,e14.7)

      end