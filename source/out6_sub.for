      subroutine out6(imo,klev,fnout,fngrd)
c-----------------------------------------------------------------------
c write formatted data
c-----------------------------------------------------------------------
      use alloc
      use refstate
	

      implicit real*8(a-h,o-z)
c      include 'nh3dpa10.inc'

      character*80 fnout,fngrd, fext
      real*8 s_arr(ns),dsig(ns)
	real,allocatable:: z_s(:,:,:),v1(:,:,:),u1(:,:,:)
	real,allocatable:: prgrx(:,:,:),prgrxs(:,:,:)
     :                  ,prgry(:,:,:),prgrys(:,:,:)
	real koriolis
	integer xminen,xmaxen,yminen,ymaxen
	s_arr=sigma1
     	allocate (z_s(nx,ny,ns),v1(nx,ny,ns),u1(nx,ny,ns))
	do is=2,ns-1
      dsig(is)=s_arr(is+1)-s_arr(is-1)
	enddo

      if(proout) call wripro(fngrd)

c
      if(.not.grdout) then
        if(norout) then
           tk=0.
           call wri3ar(u(0,0,0,klev),0,nx1,1,ny1,1,ns-1,'u       ',tk
     :       ,imo,nx,ny)
           call wri3ar(v(0,0,0,klev),1,nx1,0,ny1,1,ns-1,'v       ',tk
     :       ,imo,nx,ny)
           call wri3ar(w(0,0,0,klev),1,nx1,1,ny1,1,ns,'w       ',tk
     :       ,imo,nx,ny)
           call wri3ar(pt(0,0,0,klev),1,nx1,1,ny1,1,ns-1,'pt       ',tk
     :        ,imo,nx,ny)
           if(qif.ne.0.) then
           call wri3ar(qv(0,0,0,klev),1,nx1,1,ny1,1,ns-1,'qv       ',tk
     :         ,imo,nx,ny)
             call wri3ar(cond,1,nx1,1,ny1,1,ns-1,'cond     ',tk
     :         ,imo,nx,ny)
             call wri3ar(evap,1,nx1,1,ny1,1,ns-1,'evap     ',tk
     :         ,imo,nx,ny)
           endif
           if(ifqc.ne.0) then
           call wri3ar(qc(0,0,0,klev),1,nx1,1,ny1,1,ns-1,'qc       ',tk
     :         ,imo,nx,ny)
           endif
           if(ifqr.ne.0) then
           call wri3ar(qr(0,0,0,klev),1,nx1,1,ny1,1,ns-1,'qr       ',tk
     :         ,imo,nx,ny)
           endif
        endif
        if(morout) then
          tk=0.
          call wri3ar(phi,1,nx1,1,ny1,0,ns,'phi     ',tk,imo,nx,ny)
          call wri3ar(wsig(0,0,0,klev),1,nx1,1,ny1,0,ns1
     :       ,'wsig    ',tk,imo,nx,ny)
          if(qif.ne.0.) then
            allocate(rh(0:nx1,0:ny1,0:ns1))
            do ix=1,nx1
              do iy=1,ny1
                do is=1,ns
                  p=sigma0(is)*pp(ix,iy,klev)+ptop
                  t=(pt(ix,iy,is,klev)+pts(ix,iy,is))*(p/p00)**akapa
                  rh(ix,iy,is)=qv(ix,iy,is,klev)/qsat(t,p)
                enddo
              enddo
            enddo
            call wri3ar(rh,1,nx1,1,ny1,1,ns-1,'RH      ',tk,imo,nx,ny)
            deallocate(rh)
          endif
        endif
        if(surout.and.ifsoil.ne.0) then
           tk=273.
           call wri2ar(tsuf,1,nx1,1,ny1,'tsuf    ',tk,nx,ny)
           call wri2ar(tsnoi,2,nx,2,ny,'tsnoi   ',tk,nx,ny)
           call wri2ar(tsurw,2,nx,2,ny,'tsurw   ',tk,nx,ny)
               tmsur=xlake*tsurw+(1.-xlake)*tsnoi
           call wri2ar(tmsur,2,nx,2,ny,'tmsur   ',tk,nx,ny)
           call wri2ar(t2noi,2,nx,2,ny,'t2noi   ',tk,nx,ny)
           tk=0.
           call wri2ar(wgnoi,2,nx,2,ny,'wgnoi   ',tk,nx,ny)
           call wri2ar(w2noi,2,nx,2,ny,'w2noi   ',tk,nx,ny)
           call wri2ar(wrnoi,2,nx,2,ny,'wrnoi   ',tk,nx,ny)
           call wri2ar(h,2,nx,2,ny,'entalp  ',tk,nx,ny)
           call wri2ar(le,2,nx,2,ny,'latent  ',tk,nx,ny)
           call wri2ar(gsolo,2,nx,2,ny,'gsol    ',tk,nx,ny)
           call wri2ar(rn,2,nx,2,ny,'netrad  ',tk,nx,ny)
           call wri2ar(cdm,2,nx,2,ny,'cdm     ',tk,nx,ny)
           tk=0.
c           call wri2ar(phig,1,nx1,1,ny1,'phig    ',tk,nx,ny)
           call wri2ar(psuf,1,nx1,1,ny1,'psuf    ',tk,nx,ny)
c           call wri2ar(dpdy01,1,nx1,1,ny1,'dpdy01 ',tk,nx,ny)
c           call wri2ar(dpdx10,1,nx1,1,ny1,'dpdx10 ',tk,nx,ny)
c
        endif

        call wri2ar(deltaz,1,nx1,1,ny1,'deltaz  ',tk,nx,ny)
        if(ifqr.ne.0) then
          call wri2ar(prec,1,nx1,1,ny1,'prec    ',tk,nx,ny)
        endif
c
        avpsuf=0.
        do iy=2,ny
        do ix=2,nx
          avpsuf=avpsuf+psuf(ix,iy)
        enddo
        enddo
        avpsuf=avpsuf/float((nx-1)*(ny-1))
        write(nchn1,1000) avpsuf/100.
c
        if(refout) then
          tk=0.
          call wri3ar(phis,1,nx1,1,ny1,0,ns,'phis    ',tk,imo,nx,ny)
          tk=273.
          call wri3ar(pts,1,nx1,1,ny1,0,ns,'pts     ',tk,imo,nx,ny)
          call wri3ar(tems,1,nx1,1,ny1,0,ns,'tems    ',tk,imo,nx,ny)
        endif
c
c       if(fluout) then
c         tk=0.
c         call wri3ar(wk6,1,nx,2,ny,0,ns,'uflux    ',tk,imo,nx,ny)
c         call wri3ar(wk7,2,nx,1,ny,0,ns,'vflux    ',tk,imo,nx,ny)
c         call wri3ar(wk8,1,nx1,1,ny1,1,ns,'wflux    ',tk,imo,nx,ny)
c       endif
c
        if(difout) then
          tk=0.
          call wri3ar(difunu,1,nx,2,ny,0,ns,'difunu  ',tk,imo,nx,ny)
          call wri3ar(difunv,2,nx,1,ny,0,ns,'difunv  ',tk,imo,nx,ny)
          call wri3ar(difunw,2,nx,2,ny,1,ns,'difunw  ',tk,imo,nx,ny)
c          call wri3ar(difunt,2,nx,2,ny,1,ns-1,'difunt  ',tk,imo,nx,ny)

          if(qif.ne.0.) then
           call wri3ar(difunq,2,nx,2,ny,1,ns-1,'difunq  ',tk,imo,nx,ny)
          endif
        endif
        if(desout) then
          tk=0.
          call wri3ar(s,1,nx1,1,ny1,1,ns,'s       ',tk,imo,nx,ny)
        endif
      else
c
c write surface line for vertical cross sections
c
        if(norout) then
          allocate(work1(0:nx1,0:ny1,0:ns1))
c
c output u,v and w at half grid points
c
          do ix=2,nx
          do iy=2,ny
          do is=1,ns-1
            work1(ix,iy,is)=0.5*(u(ix,iy,is,klev)+u(ix-1,iy,is,klev))
	      u1(ix,iy,is)=0.5*(u(ix,iy,is,klev)+u(ix-1,iy,is,klev))
          enddo
          enddo
          enddo
          call wgrids(work1,'um',2,nx,2,ny,1,ns-1,imo,0,0,0,fngrd)
          do ix=2,nx
          do iy=2,ny
          do is=1,ns-1
            work1(ix,iy,is)=0.5*(v(ix,iy,is,klev)+v(ix,iy-1,is,klev))
	      v1(ix,iy,is)=0.5*(v(ix,iy,is,klev)+v(ix,iy-1,is,klev))
          enddo
          enddo
          enddo
          call wgrids(work1,'vm',2,nx,2,ny,1,ns-1,imo,0,0,0,fngrd)
          do ix=2,nx
          do iy=2,ny
          do is=1,ns-1
            work1(ix,iy,is)=0.5*(w(ix,iy,is,klev)+w(ix,iy,is+1,klev))
	    enddo
          enddo
          enddo
          call wgrids(work1,'wm',2,nx,2,ny,1,ns-1,imo,0,0,0,fngrd)
!          deallocate(work1)
      !------------------output for divergence chlen---------------------------------------------------!
	     do ix=1,nx
          do iy=1,ny
          do is=1,ns
           z_s(ix,iy,is)=(phis(ix,iy,is)+phi(ix,iy,is))/g
	    enddo
	    enddo
	    enddo
!	koriolis=1.367/(10.**4) !koriolis parameter for lat=70
	    
!	     do ix=2,nx-1
!          do iy=2,ny-1
!          do is=2,ns-1
	    
!	    work1(ix,iy,is)=
!     &((v1(ix+1,iy,is)-v1(ix-1,iy,is))/(3*dx) -
!     &(u1(ix,iy+1,is)-u1(ix,iy-1,is))/(3*dy) +koriolis -
!     &(z_s(ix+1,iy,is)-z_s(ix-1,iy,is))/(2*dx)*
!     &dsig(is)/(z_s(ix,iy,is+1)-z_s(ix,iy,is-1))*
!     &(v1(ix,iy,is+1)-v1(ix,iy,is-1))/dsig(is)+
!     &(z_s(ix,iy+1,is)-z_s(ix,iy-1,is))/(2*dy)*
!     &dsig(is)/(z_s(ix,iy,is+1)-z_s(ix,iy,is-1))*
!     &(u1(ix,iy,is+1)-u1(ix,iy,is-1))/dsig(is))*	
!     &((u(ix,iy,is,klev)-u(ix-1,iy,is,klev))/dx+
!     &(v(ix,iy,is,klev)-v(ix,iy-1,is,klev))/dy-
!     &(z_s(ix+1,iy,is)-z_s(ix-1,iy,is))/(2*dx)*
!     &dsig(is)/(z_s(ix,iy,is+1)-z_s(ix,iy,is-1))*
!     &(u1(ix,iy,is+1)-u1(ix,iy,is-1))/dsig(is)-
!     &(z_s(ix,iy+1,is)-z_s(ix,iy-1,is))/(2*dy)*
!     &dsig(is)/(z_s(ix,iy,is+1)-z_s(ix,iy,is-1))*
!     &(v1(ix,iy,is+1)-v1(ix,iy,is-1))/dsig(is))	
            
!          enddo
!          enddo
!          enddo
!       call wgrids(work1,'di',3,nx-1,3,ny-1,2,ns-1,imo,0,0,0,fngrd)
!	    do ix=2,nx-1
!          do iy=2,ny-1
!          do is=2,ns-1
	    
!	    work1(ix,iy,is)=
!     &((v1(ix+1,iy,is)-v1(ix-1,iy,is))/(3*dx) -
!     &(u1(ix,iy+1,is)-u1(ix,iy-1,is))/(3*dy) -
!     &(z_s(ix+1,iy,is)-z_s(ix-1,iy,is))/(2*dx)*
!     &dsig(is)/(z_s(ix,iy,is+1)-z_s(ix,iy,is-1))*
!     &(v1(ix,iy,is+1)-v1(ix,iy,is-1))/dsig(is)+
!     &(z_s(ix,iy+1,is)-z_s(ix,iy-1,is))/(2*dy)*
!     &dsig(is)/(z_s(ix,iy,is+1)-z_s(ix,iy,is-1))*
!     &(u1(ix,iy,is+1)-u1(ix,iy,is-1))/dsig(is))	
     	
            
  !        enddo
  !        enddo
  !        enddo
!	call wgrids(work1,'vo',3,nx-1,3,ny-1,2,ns-1,imo,0,0,0,fngrd)
      !-----------------------------------------------------------------------------!

! following not needed any more (post files for surfer)

!         call wrivec(u(0,0,0,klev),v(0,0,0,klev),'uv',1,fngrd
!    :      ,2,nx,2,ny)
!         call wrivec(u(0,0,0,klev),w(0,0,0,klev),'uw',2,fngrd
!    :      ,2,nx,1,ns)
!         call wrivec(v(0,0,0,klev),w(0,0,0,klev),'vw',3,fngrd
!    :      ,2,ny,1,ns)
          
          call wgrids(pt(0,0,0,klev),'pt',2,nx,2,ny,1,ns-1,imo,0,0,0
     :      ,fngrd)
          call wgrids(pts(0,0,0),'rp',2,nx,2,ny,1,ns-1,imo,0,0,0
     :      ,fngrd)
          allocate(theta(0:nx1,0:ny1,0:ns1))
                    
          theta=pt(:,:,:,klev)+pts
          call wgrids(theta,'th',2,nx,2,ny,1,ns-1,imo,0,0,0
     :      ,fngrd)
     
          allocate(pressure(0:nx1,0:ny1,0:ns1))
          allocate(temperatur(0:nx1,0:ny1,0:ns1))
          do ix=1,nx1
            do iy=1,ny1
              do is=1,ns
                pressure(ix,iy,is)=sigma0(is)*pp(ix,iy,klev)+ptop
                temperatur(ix,iy,is)=theta(ix,iy,is)*
     :           (pressure(ix,iy,is)/p00)**akapa
              enddo
            enddo
          enddo
          call wgrids(pressure(0,0,0),'pres',2,nx,2,ny,1,ns-1,imo,0,0,0
     :      ,fngrd)
          call wgrids(temperatur(0,0,0),'tem',2,nx,2,ny,1,ns-1,imo,0,0,0
     :      ,fngrd)
          deallocate(pressure)
          deallocate(temperatur)

c         Dima@ 27.04.07
          if (impur.NE.0) then                                                  DM,05.2007
            call wgrids(qs_10(0,0,0,klev),'%qs',2,nx,2,ny,1,ns-1,imo,0,0        DM,05.2007
     :        ,0,fngrd)                                                         DM,05.2007
          endif                                                                 DM,05.2007
          if(ifqc.ne.0) then
            call wgrids(qc(0,0,0,klev),'qc',2,nx,2,ny,1,ns-1,imo,0,0,0
     :        ,fngrd)
          endif
          if(ifqr.ne.0) then
            call wgrids(qr(0,0,0,klev),'qr',2,nx,2,ny,1,ns-1,imo,0,0,0
     :        ,fngrd)
          endif
	    if(ifqi.ne.0) then                                                    DC,11.2009
            call wgrids(qci(0,0,0,klev),'q1',2,nx,2,ny,1,ns-1,imo,0,0,0         DC,11.2009
     :        ,fngrd)                                                           DC,11.2009
	      call wgrids(qsn(0,0,0,klev),'q2',2,nx,2,ny,1,ns-1,imo,0,0,0         DC,11.2009
     :        ,fngrd)                                                           DC,11.2009
          endif                                                                 DC,11.2009
        endif
        if(morout) then
          call wgrids(u(0,0,0,klev),'u_',1,nx,2,ny,1,ns-1,imo,1,0,0
     :      ,fngrd)
          call wgrids(v(0,0,0,klev),'v_',2,nx,1,ny,1,ns-1,imo,0,1,0
     :      ,fngrd)
          call wgrids(w(0,0,0,klev),'w_',2,nx,2,ny,1,ns,imo,0,0,1
     :      ,fngrd)
          tk=0.
          call wgrids(phi,'ph',2,nx,2,ny,1,ns-1,imo,0,0,0,fngrd)
          call wgrids(wsig,'ws',2,nx,2,ny,1,ns,imo,0,0,0,fngrd)
          if(qif.ne.0.) then
            call wgrids(qv(0,0,0,klev),'qv',2,nx,2,ny,1,ns-1,imo,0,0,0
     :        ,fngrd)
            call wgrids(cond,'co',2,nx,2,ny,1,ns-1,imo,0,0,0
     :        ,fngrd)
            call wgrids(evap,'ev',2,nx,2,ny,1,ns-1,imo,0,0,0
     :        ,fngrd)
	    if (ifqi.ne.0) then
	    call wgrids(vini,'f1',2,nx,2,ny,1,ns-1,imo,0,0,0
     :        ,fngrd)
	    call wgrids(vdeps,'f2',2,nx,2,ny,1,ns-1,imo,0,0,0
     :        ,fngrd)
	    call wgrids(vdepi,'f3',2,nx,2,ny,1,ns-1,imo,0,0,0
     :        ,fngrd)
	    call wgrids(hmfrz,'f4',2,nx,2,ny,1,ns-1,imo,0,0,0
     :        ,fngrd)
	    call wgrids(sacrw,'f5',2,nx,2,ny,1,ns-1,imo,0,0,0
     :        ,fngrd)
	    call wgrids(sacrwr,'f6',2,nx,2,ny,1,ns-1,imo,0,0,0
     :        ,fngrd)
	    call wgrids(saci,'f7',2,nx,2,ny,1,ns-1,imo,0,0,0
     :        ,fngrd)
	    call wgrids(sagg,'f8',2,nx,2,ny,1,ns-1,imo,0,0,0
     :        ,fngrd)
	    call wgrids(raci,'f9',2,nx,2,ny,1,ns-1,imo,0,0,0
     :        ,fngrd)
	    call wgrids(iacr,'g1',2,nx,2,ny,1,ns-1,imo,0,0,0
     :        ,fngrd)
	    call wgrids(sbercw,'g2',2,nx,2,ny,1,ns-1,imo,0,0,0
     :        ,fngrd)
	    call wgrids(sberci,'g3',2,nx,2,ny,1,ns-1,imo,0,0,0
     :        ,fngrd)
	    call wgrids(ibercw,'g4',2,nx,2,ny,1,ns-1,imo,0,0,0
     :        ,fngrd)
	    call wgrids(sacrr,'g5',2,nx,2,ny,1,ns-1,imo,0,0,0
     :        ,fngrd)
	    call wgrids(smlt,'g6',2,nx,2,ny,1,ns-1,imo,0,0,0
     :        ,fngrd)
	    call wgrids(col,'g7',2,nx,2,ny,1,ns-1,imo,0,0,0
     :        ,fngrd)
	    call wgrids(auto,'g8',2,nx,2,ny,1,ns-1,imo,0,0,0
     :        ,fngrd)
	    endif
            allocate(rh(0:nx1,0:ny1,0:ns1))
            do ix=1,nx1
              do iy=1,ny1
                do is=1,ns
                  p=sigma0(is)*pp(ix,iy,klev)+ptop
                  t=theta(ix,iy,is)*(p/p00)**akapa
                  rh(ix,iy,is)=qv(ix,iy,is,klev)/qsat(t,p)
                  
                enddo
              enddo
            enddo
            call wgrids(rh,'rh',2,nx,2,ny,1,ns-1,imo,0,0,0,fngrd)
            deallocate(rh)
          endif
        endif
        deallocate(theta)
     
!--------output of pressure gradient terms in momentum equations-------!        
        allocate(prgrx(0:nx1,0:ny1,0:ns1),
     :           prgry(0:nx1,0:ny1,0:ns1))
        prgrx=0.
        prgry=0.

     
        do is=1,ns-1
        do iy=2,ny
        do ix=2,nx-1
          prgrx(ix,iy,is)=(-(phi(ix+1,iy,is)-phi(ix,iy,is))/dx
     :    +s0ds4a(is)*ppdx10(ix,iy,2)*(phi(ix+1,iy,is+1)
     :    +phi(ix,iy,is+1)-phi(ix+1,iy,is-1)-phi(ix,iy,is-1)))
     :    *pp01(ix,iy,2)
       enddo
       enddo
       enddo
       
       do is=1,ns-1
       do iy=2,ny-1
       do ix=2,nx
          prgry(ix,iy,is)=(-(phi(ix,iy+1,is)-phi(ix,iy,is))/dy
     :    +s0ds4a(is)*ppdy01(ix,iy,2)*(phi(ix,iy+1,is+1)
     :    +phi(ix,iy,is+1)-phi(ix,iy+1,is-1)-phi(ix,iy,is-1)))
     :    *pp01(ix,iy,2)
       enddo
       enddo
       enddo
       
       call extrah(nx1,ny1,prgrx,1,nx,1,ny1,0,ns)
       call extrah(nx1,ny1,prgry,1,nx1,1,ny,0,ns)
       
       call wgrids(prgrx,'px',2,nx,2,ny,1,ns-1,imo,0,0,0,fngrd)
       call wgrids(prgry,'py',2,nx,2,ny,1,ns-1,imo,0,0,0,fngrd)
       
       deallocate(prgrx,prgry)
       call wgrids(ufd,'ud',2,nx,2,ny,1,ns-1,imo,0,0,0,fngrd)
       call wgrids(vfd,'vd',2,nx,2,ny,1,ns-1,imo,0,0,0,fngrd)
       call wgrids(wfd,'wd',2,nx,2,ny,1,ns-1,imo,0,0,0,fngrd)
!------------------------------------------------------------------------!        
!-------------------output for coriolis terms in momentum equations------!
      work1=0.
      if(fcor.ne.0.) then
        do is=0,ns
        do iy=2,ny
        do ix=1,nx
          work1(ix,iy,is)=-fcor*0.25*(v(ix,iy,is,3)
     :     *pp01(ix,iy,3)
     :     +v(ix+1,iy,is,3)*pp01(ix+1,iy,3)+v(ix,iy-1,is,3)
     :     *pp01(ix,iy-1,3)+v(ix+1,iy-1,is,3)*pp01(ix+1,iy-1,3))
        enddo
        enddo
        enddo
        call wgrids(work1,'cx',2,nx,2,ny,1,ns-1,imo,0,0,0,fngrd)
        do is=0,ns
        do iy=1,ny
        do ix=2,nx
            work1(ix,iy,is)=fcor*0.25*(u(ix,iy,is,3)
     :      *pp10(ix,iy,3)
     :      +u(ix-1,iy,is,3)*pp10(ix-1,iy,3)+u(ix,iy+1,is,3)
     :      *pp10(ix,iy+1,3)+u(ix-1,iy+1,is,3)*pp10(ix-1,iy+1,3))
        enddo
        enddo
        enddo
        call wgrids(work1,'cy',2,nx,2,ny,1,ns-1,imo,0,0,0,fngrd)
      endif


!------------------------------------------------------------------------!

!-------------output of heat budget terms--------------------------------!
         do is=1,ns-1
         do iy=2,ny
         do ix=2,nx
           work1(ix,iy,is)=((s(ix,iy,is)+s(ix,iy,is+1))*(pts(ix,iy,is+1)
     :     -pts(ix,iy,is-1))
     :     *(w(ix,iy,is+1,2)+w(ix,iy,is,2))/ds08a(is))*pp(ix,iy,2)
         enddo
         enddo
         enddo
         call wgrids(work1,'tq',2,nx,2,ny,1,ns-1,imo,0,0,0,fngrd)
         do is=1,ns-1
         do iy=2,ny
         do ix=2,nx
           work1(ix,iy,is)=(u(ix,iy,is,2)*pp10(ix,iy,2)
     :       *(pt(ix+1,iy,is,2)+pt(ix,iy,is,2))
     :       -u(ix-1,iy,is,2)*pp10(ix-1,iy,2)
     :       *(pt(ix,iy,is,2)+pt(ix-1,iy,is,2)))/dx2
         enddo
         enddo
         enddo
         call wgrids(work1,'tx',2,nx,2,ny,1,ns-1,imo,0,0,0,fngrd)
         do is=1,ns-1
         do iy=2,ny
         do ix=2,nx
           work1(ix,iy,is)=(v(ix,iy,is,2)*pp01(ix,iy,2)
     :       *(pt(ix,iy+1,is,2)+pt(ix,iy,is,2))
     :       -v(ix,iy-1,is,2)*pp01(ix,iy-1,2)
     :       *(pt(ix,iy,is,2)+pt(ix,iy-1,is,2)))/dy2
         enddo
         enddo
         enddo
         call wgrids(work1,'ty',2,nx,2,ny,1,ns-1,imo,0,0,0,fngrd)
         do is=1,ns-1
         do iy=2,ny
         do ix=2,nx
           work1(ix,iy,is)=pp(ix,iy,2)*
     :      ((wsig(ix,iy,is+1,2)*(pt(ix,iy,is+1,2)+pt(ix,iy,is,2))
     :       -wsig(ix,iy,is,2)*(pt(ix,iy,is,2)+pt(ix,iy,is-1,2)))
     :       /ds12(is))
          enddo
          enddo
          enddo
          call wgrids(work1,'tc',2,nx,2,ny,1,ns-1,imo,0,0,0,fngrd)
          call wgrids(ptfd,'td',2,nx,2,ny,1,ns-1,imo,0,0,0,fngrd)

!--------------output of surface stuff in case of fixed surf temp--------!        
	  if (tfix.eq.1) then                                                          DC,07.2010
	     call wgrids(h,'hh',2,nx,2,ny,0,0,1,0,0,0,fngrd)
	     call wgrids(tsurf,'ts',2,nx,2,ny,0,0,1,0,0,0,fngrd)                       DC,07.2010
	     call wgrids(imp,'im',2,nx,2,ny,0,0,1,0,0,0,fngrd)
	     call wgrids(deltaz,'dz',2,nx,2,ny,0,0,1,0,0,0,fngrd)                      DC,07.2010
           call wgrids(le,'lh',2,nx,2,ny,0,0,1,0,0,0,fngrd)                          DC,07.2010
 !          call wgrids(cdm,'cd',2,nx,2,ny,0,0,1,0,0,0,fngrd)                        DC,07.2010
           call wgrids(ust_s,'uz',2,nx,2,ny,0,0,1,0,0,0,fngrd)                       DC,07.2010
	  endif                                                                


        if(surout) then
          call wgrids(tsuf,'t0',2,nx,2,ny,0,0,1,0,0,0,fngrd)
          if(ifsoil.ne.0) then
            call wgrids(tsnoi,'ts',2,nx,2,ny,0,0,1,0,0,0,fngrd)
            tmsur(:,:)=xlake(:,:)*tsurw+(1.-xlake)*tsnoi
            call wgrids(tsurw,'tw',2,nx,2,ny,0,0,1,0,0,0,fngrd)
            call wgrids(tmsur,'tm',2,nx,2,ny,0,0,1,0,0,0,fngrd)
            call wgrids(t2noi,'t2',2,nx,2,ny,0,0,1,0,0,0,fngrd)
            call wgrids(wgnoi,'wg',2,nx,2,ny,0,0,1,0,0,0,fngrd)
            call wgrids(w2noi,'w2',2,nx,2,ny,0,0,1,0,0,0,fngrd)
            call wgrids(wrnoi,'wr',2,nx,2,ny,0,0,1,0,0,0,fngrd)
            call wgrids(h,'hh',2,nx,2,ny,0,0,1,0,0,0,fngrd)
            call wgrids(le,'lh',2,nx,2,ny,0,0,1,0,0,0,fngrd)
             call wgrids(watice,'wi',2,nx,2,ny,0,0,1,0,0,0,fngrd)               DC,04.2009
            call wgrids(gsolo,'gs',2,nx,2,ny,0,0,1,0,0,0,fngrd)
            call wgrids(rn,'rn',2,nx,2,ny,0,0,1,0,0,0,fngrd)
            call wgrids(cdm,'cd',2,nx,2,ny,0,0,1,0,0,0,fngrd)
            call wgrids(deltaz,'dz',2,nx,2,ny,0,0,1,0,0,0,fngrd)
            call wgrids(rat,   'ra',2,nx,2,ny,0,0,1,0,0,0,fngrd)                VS,06.2007
	      call wgrids(rg,    'rg',2,nx,2,ny,0,0,1,0,0,0,fngrd)                VS,06.2007
          endif
          call wgrids(psuf,'ps',2,nx,2,ny,0,0,1,0,0,0,fngrd)
          panom=psuf-pini
          call wgrids(panom,'px',2,nx,2,ny,0,0,1,0,0,0,fngrd)
          if(ifqr.ne.0) then
            call wgrids(prec,'pr',2,nx,2,ny,0,0,1,0,0,0,fngrd)
          endif
        endif

        if(norout.or.surout) then
          call wgrids(hsuf,'hs',1,nx+1,1,ny+1,0,0,1,0,0,0,fngrd)
        endif

        if(refout) then
          call wgrids(phis,'fs',2,nx,2,ny,1,ns-1,imo,0,0,0,fngrd)
          call wgrids(pts,'os',2,nx,2,ny,1,ns-1,imo,0,0,0,fngrd)
          call wgrids(tems,'tt',2,nx,2,ny,1,ns-1,imo,0,0,0,fngrd)
          call wgrids(psufref,'pp',2,nx,2,ny,0,0,1,0,0,0,fngrd)
c
c output us,vs  at half grid points
          if(.not.allocated(work1)) allocate(work1(0:nx1,0:ny1,0:ns1))
c
          do ix=2,nx
          do iy=2,ny
          do is=1,ns-1
            work1(ix,iy,is)=0.5*(us(ix,iy,is)+us(ix-1,iy,is))
          enddo
          enddo
          enddo
          call wgrids(work1,'us',2,nx,2,ny,1,ns-1,imo,0,0,0,fngrd)
c
c output us,vs  at half grid points
c
          do ix=2,nx
          do iy=2,ny
          do is=1,ns-1
            work1(ix,iy,is)=0.5*(vs(ix,iy,is)+vs(ix,iy-1,is))
          enddo
          enddo
          enddo
          call wgrids(work1,'vs',2,nx,2,ny,1,ns-1,imo,0,0,0,fngrd)
        endif
c
       if(fluout) then
         call wgrids(wk6,'uf',1,nx,2,ny,1,ns-1,imo,1,0,0,fngrd)
         call wgrids(wk7,'vf',1,nx,2,ny,1,ns-1,imo,0,1,0,fngrd)
         call wgrids(wk8,'wf',1,nx,2,ny,1,ns-1,imo,0,0,1,fngrd)
       endif
c
        if(difout) then
          call wgrids(difunu,'du',1,nx,2,ny,1,ns-1,imo,1,0,0,fngrd)
          call wgrids(difunv,'dv',2,nx,1,ny,1,ns-1,imo,0,1,0,fngrd)
          call wgrids(difunw,'dw',2,nx,2,ny,1,ns,imo,0,0,1,fngrd)
          call wgrids(difunt,'dt',2,nx,2,ny,1,ns-1,imo,0,0,0,fngrd)
          if(qif.ne.0.) then
            call wgrids(difunq,'dq',2,nx,2,ny,1,ns-1,imo,0,0,0,fngrd)
            if(ifqc.ne.0) then
             call wgrids(difunqc,'dc',2,nx,2,ny,1,ns-1,imo,0,0,0,fngrd)
              if(ifqr.ne.0) then
                call wgrids(difunqr,'dr',2,nx,2,ny,1,ns-1,imo,0,0,0
     :            ,fngrd)

              endif
            endif
          endif
        endif
        if(desout) then
          call wgrids(s,'s_',1,nx,2,ny,1,ns-1,imo,0,0,0,fngrd)
        endif
      endif

      if(allocated(work1)) deallocate(work1)
      if(allocated(u1)) deallocate(u1)
      if(allocated(v1)) deallocate(v1)
      if(allocated(z_s)) deallocate(z_s)
      

      if(eneout) call energy(1,nx,1,ny,1,ns                       
     :  ,nx,ny,ns
     :  ,u(0,0,0,3),v(0,0,0,3),w(0,0,0,3),pt(0,0,0,3)
     :  ,pts,phi,phis,pp(0,0,3),pp10(0,0,3),pp01(0,0,3)
     :  ,sigma0,sigma1,ds1x,ds1y,dxys,dxy,ds02a,ptop,nchn1,67)    !write kinetic energy, integr. over the whole domain into .en1
      if(eneout) call energy(1,nx,ny/2,ny,1,ns
     :  ,nx,ny,ns
     :  ,u(0,0,0,3),v(0,0,0,3),w(0,0,0,3),pt(0,0,0,3)
     :  ,pts,phi,phis,pp(0,0,3),pp10(0,0,3),pp01(0,0,3)
     :  ,sigma0,sigma1,ds1x,ds1y,dxys,dxy,ds02a,ptop,nchn1,68)
      if(eneout) call energy(1,nx,1,ny/2,1,ns
     :  ,nx,ny,ns
     :  ,u(0,0,0,3),v(0,0,0,3),w(0,0,0,3),pt(0,0,0,3)
     :  ,pts,phi,phis,pp(0,0,3),pp10(0,0,3),pp01(0,0,3)
     :  ,sigma0,sigma1,ds1x,ds1y,dxys,dxy,ds02a,ptop,nchn1,69)
      if(eneout) call energy(1,nx,1,ny,ns/2,ns
     :  ,nx,ny,ns
     :  ,u(0,0,0,3),v(0,0,0,3),w(0,0,0,3),pt(0,0,0,3)
     :  ,pts,phi,phis,pp(0,0,3),pp10(0,0,3),pp01(0,0,3)
     :  ,sigma0,sigma1,ds1x,ds1y,dxys,dxy,ds02a,ptop,nchn1,70)
       xminen=(nx+1)/2-1
       xmaxen=(nx+1)/2
       yminen=xice-5
       ymaxen=xice+1  
       write(0,*) 'boundaries for energy',xminen,xmaxen,yminen,ymaxen 
      if(eneout) call energy(xminen,xmaxen,yminen,ymaxen,1,ns
     :  ,nx,ny,ns
     :  ,u(0,0,0,3),v(0,0,0,3),w(0,0,0,3),pt(0,0,0,3)
     :  ,pts,phi,phis,pp(0,0,3),pp10(0,0,3),pp01(0,0,3)
     :  ,sigma0,sigma1,ds1x,ds1y,dxys,dxy,ds02a,ptop,nchn1,71)
c
      return
c
1000  format(//1x,'average psuf=',f15.8,' hpa')
      end
