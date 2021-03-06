       subroutine thermo(timeho,ptm2,wrk,fngrd)
c-----------------------------------------------------------------------
c     time integration of potential temperature equation
c-----------------------------------------------------------------------
c
      use alloc
     

      implicit real*8(a-h,o-z)
      dimension ptm2(0:nx1,0:ny1,0:ns1)
      dimension wrk(0:nx1,0:ny1,0:ns1)
      character*80 fngrd


      if(tfct) then
         do k=2,3
         do is=0,ns
         do iy=1,ny1
         do ix=1,nx1
           pt(ix,iy,is,k)=pt(ix,iy,is,k)*pp(ix,iy,k)
         enddo
         enddo
         enddo
         enddo
         do is=0,ns
         do iy=1,ny1
         do ix=1,nx1
           ptm2(ix,iy,is)=ptm2(ix,iy,is)*pp(ix,iy,1)
         enddo
         enddo
         enddo
         allocate(work1(0:nx1,0:ny1,0:ns1))
         allocate(work2(0:nx1,0:ny1,0:ns1))
         allocate(work3(0:nx1,0:ny1,0:ns1))
         allocate(work4(0:nx1,0:ny1,0:ns1))
         allocate(work5(0:nx1,0:ny1,0:ns1))
         allocate(work6(0:nx1,0:ny1,0:ns1))
         allocate(work7(0:nx1,0:ny1,0:ns1))
         allocate(work8(0:nx1,0:ny1,0:ns1))
         call fct2d(pt(0,0,0,3),pt(0,0,0,2),ptm2
     :     ,u(0,0,0,3),u(0,0,0,2)
     :     ,v(0,0,0,3),v(0,0,0,2)
     :     ,wsig(0,0,0,3),wsig(0,0,0,2)
     :     ,nx,ny,ns
     :     ,dtxs,dtys,dtxy,dxys,dxdt,dydt,dt,ds0
     :     ,work1,work2,work3,work4,work5,work6,work7,work8,prt,grdout,
     :      iomod,fngrd)
         deallocate(work1)
         deallocate(work2)
         deallocate(work3)
         deallocate(work4)
         deallocate(work5)
         deallocate(work6)
         deallocate(work7)
         deallocate(work8)

         do k=2,3
         do is=0,ns
         do iy=1,ny1
         do ix=1,nx1
           pt(ix,iy,is,k)=pt(ix,iy,is,k)/pp(ix,iy,k)
         enddo
         enddo
         enddo
         enddo
         do is=0,ns
         do iy=1,ny1
         do ix=1,nx1
           ptm2(ix,iy,is)=ptm2(ix,iy,is)/pp(ix,iy,1)
         enddo
         enddo
         enddo
	   do is=1,ns-1
         do iy=2,ny
         do ix=2,nx
           q=(s(ix,iy,is)+s(ix,iy,is+1))*(pts(ix,iy,is+1)
     :       -pts(ix,iy,is-1))*(w(ix,iy,is+1,2)+w(ix,iy,is,2))/ds08a(is)
           p=sigma0(is)*pp(ix,iy,2)+ptop
           constt=hlatcp*(p00/p)**akapa
!	           if(ifqc.ne.0.and.ifqi.eq.0) then                                                   DC,11.2009
!             pt(ix,iy,is,3)=pt(ix,iy,is,3)+dtl*(                                DC,11.2009
!     :       pp(ix,iy,2)*(q+difunt(ix,iy,is)                                    DC,11.2009
!     :         +constt*(cond(ix,iy,is)-evap(ix,iy,is))))/pp(ix,iy,3) 
!	     elseif(ifqi.ne.0) then           
!		 pt(ix,iy,is,3)=pt(ix,iy,is,3)+dtl*(                                  DC,11.2009
!     :       pp(ix,iy,2)*(q+difunt(ix,iy,is)                                    DC,11.2009
!     :         +constt*(cond(ix,iy,is)-evap(ix,iy,is))                                       DC,11.2009
!     :         +hsub/cp*(p00/p)**akapa*(vini(ix,iy,is)+vdeps(ix,iy,is)+         DC,11.2009
!     :         vdepi(ix,iy,is))+hfus/cp*(p00/p)**akapa*(                        DC,11.2009
!     :         hmfrz(ix,iy,is)+sacrw(ix,iy,is)+iacr(ix,iy,is)+                  DC,11.2009
!     :         sbercw(ix,iy,is)+ibercw(ix,iy,is)+sacrr(ix,iy,is)-               DC,11.2009
!     :         smlt(ix,iy,is)-imlt(ix,iy,is))                                   DC,11.2009
!     :         ))/pp(ix,iy,3)                                                   DC,11.2009
!           else
!  this tavs is only used here as a test for when FCT scheme used only in horizontal
            tavs=(wsig(ix,iy,is+1,2)*(pt(ix,iy,is+1,2)+pt(ix,iy,is,2))
     :       -wsig(ix,iy,is,2)*(pt(ix,iy,is,2)+pt(ix,iy,is-1,2)))
     :       /ds12(is)
             if(qif.eq.0) then
               pt(ix,iy,is,3)=pt(ix,iy,is,3)+dtl*(
     :         pp(ix,iy,2)*(-tavs+q+difunt(ix,iy,is)))/pp(ix,iy,3)
             elseif(qif.ne.0.and.ifqc.eq.0) then
               pt(ix,iy,is,3)=pt(ix,iy,is,3)+dtl*(
     :         pp(ix,iy,2)*(-tavs+q+difunt(ix,iy,is)
     :         -hlat/cp*difunq(ix,iy,is)))/pp(ix,iy,3)         !-hlat/cp*difunq(ix,iy,is)
             elseif (ifqc.ne.0) then 
               pt(ix,iy,is,3)=pt(ix,iy,is,3)+dtl*(
     :         pp(ix,iy,2)*(-tavs+q+difunt(ix,iy,is)
     :         -hlat/cp*difunq(ix,iy,is)
     :         +constt*(cond(ix,iy,is)-evap(ix,iy,is))))/pp(ix,iy,3)
             endif
!           endif
         enddo
         enddo
         enddo
        
c
      elseif(centered) then
         
        
      !simple centered difference scheme for advection, upstream in horizontal as an option
         do is=1,ns-1
         do iy=2,ny
         do ix=2,nx
           q=(s(ix,iy,is)+s(ix,iy,is+1))*(pts(ix,iy,is+1)
     :       -pts(ix,iy,is-1))
     :       *(w(ix,iy,is+1,2)+w(ix,iy,is,2))/ds08a(is)

           tavx=(u(ix,iy,is,2)*pp10(ix,iy,2)
     :       *(pt(ix+1,iy,is,2)+pt(ix,iy,is,2))
     :       -u(ix-1,iy,is,2)*pp10(ix-1,iy,2)
     :       *(pt(ix,iy,is,2)+pt(ix-1,iy,is,2)))/dx2
           advfx(ix,iy,is)=u(ix,iy,is,2)*pp10(ix,iy,2)*0.5
     :       *(pt(ix+1,iy,is,2)+pt(ix,iy,is,2))
! ------simple upstream scheme for x direction---------------------!   
 !           if(u(ix,iy,is,2).ge.0) then
 !              tavx=
 !              (u(ix-1,iy,is,2)*pt(ix-1,iy,is,2)*pp10(ix-1,iy,2)
 !    :          -u(ix,iy,is,2)*pt(ix,iy,is,2)*pp10(ix,iy,2))/dx
 !           else
 !              tavx=
 !             (u(ix,iy,is,2)*pp10(ix,iy,2)*pt(ix+1,iy,is,2)
 !    :         -u(ix-1,iy,is,2)*pp10(ix-1,iy,2)*pt(ix,iy,is,2))/dx       
 !           endif

!----------------normal centered difference scheme-----------------!
           tavy=(v(ix,iy,is,2)*pp01(ix,iy,2)
     :       *(pt(ix,iy+1,is,2)+pt(ix,iy,is,2))
     :       -v(ix,iy-1,is,2)*pp01(ix,iy-1,2)
     :       *(pt(ix,iy,is,2)+pt(ix,iy-1,is,2)))/dy2
     
            advfy(ix,iy,is)=v(ix,iy,is,2)*pp01(ix,iy,2)*0.5
     :       *(pt(ix,iy+1,is,2)+pt(ix,iy,is,2))
 !--------------simple upstream scheme for y direction-------------!    
 !           if(v(ix,iy,is,2).ge.0) then
 !              tavy=
 !    :          (pp01(ix,iy-1,2)*v(ix,iy-1,is,2)*pt(ix,iy-1,is,2)-
 !    :          pp01(ix,iy,2)*v(ix,iy,is,2)*pt(ix,iy,is,2))/dy     
 !           else
 !              tavy=
 !    :          (pp01(ix,iy,2)*v(ix,iy,is,2)*pt(ix,iy+1,is,2)-
 !    :         pp01(ix,iy-1,2)*v(ix,iy-1,is,2)*pt(ix,iy,is,2))/dy      
 !           endif
 !------------------------------------------------------------------!
           tavs=(wsig(ix,iy,is+1,2)*(pt(ix,iy,is+1,2)+pt(ix,iy,is,2))
     :       -wsig(ix,iy,is,2)*(pt(ix,iy,is,2)+pt(ix,iy,is-1,2)))
     :       /ds12(is)
           p=sigma0(is)*pp(ix,iy,2)+ptop
 !          if(ifqc.ne.0.) then
 !            condens=hlatcp*(p00/p)**akapa*(cond(ix,iy,is)
 !    :         -evap(ix,iy,is))
 !          else
 !            condens=0.
 !          endif
!	     if(ifqi.ne.0.) then
!             transform=
!     :         hsub/cp*(p00/p)**akapa*(vini(ix,iy,is)+vdeps(ix,iy,is)+          DC,11.2009
 !    :         vdepi(ix,iy,is))+hfus/cp*(p00/p)**akapa*(                        DC,11.2009
 !    :         hmfrz(ix,iy,is)+sacrw(ix,iy,is)+iacr(ix,iy,is)+                  DC,11.2009
 !    :         sbercw(ix,iy,is)+ibercw(ix,iy,is)+sacrr(ix,iy,is)-               DC,11.2009
 !   :         smlt(ix,iy,is)-imlt(ix,iy,is))                                   DC,11.2009
  !         else
  !           transform=0.
  !         endif
  !        if(qif.ne.0) then
!	     pt(ix,iy,is,3)=(ptm2(ix,iy,is)*pp(ix,iy,1)                           DC,11,2009
!     :       -dtl*(tavx+tavy+pp(ix,iy,2)                                        DC,11,2009
!     :       *(tavs-q-(difunt(ix,iy,is)                                          DC,11,2009
!     :       -hlat/cp*difunq(ix,iy,is))
!     :       -condens -transform                                                 DC,11,2009
!     :        )))/pp(ix,iy,3)                                                   DC,11,2009
!           endif
 !          if (qif.eq.0) then
           pt(ix,iy,is,3)=(ptm2(ix,iy,is)*pp(ix,iy,1)                           DC,11,2009
     :       -dtl*(tavx+tavy+pp(ix,iy,2)                                        DC,11,2009
     :       *(tavs-q-difunt(ix,iy,is)
     :        )))/pp(ix,iy,3) 
  !         endif
         enddo
         enddo
         enddo
         
         
         if(prt) then
          if(grdout) then
            call wgrids(advfx,'a1',2,nx,2,ny,2,ns-1,iomod
     :        ,0,0,0,fngrd)
	      call wgrids(advfy,'a2',2,nx,2,ny,2,ns-1,iomod
     :        ,0,0,0,fngrd)
	    endif
        endif
         elseif(dst3) then
         
         call dst3adv(pp10(0,0,2),pp01(0,0,2)
     :    ,pp(0,0,1),pp(0,0,3),u(0,0,0,2),v(0,0,0,2)
     :    ,ptm2(0,0,0),ptm2(0,0,0)
     :    ,adv_flx,adv_fly
     :    ,ptnew(0,0,0))
     
         do is=1,ns-1
         do iy=2,ny
         do ix=2,nx
           q=(s(ix,iy,is)+s(ix,iy,is+1))*(pts(ix,iy,is+1)
     :       -pts(ix,iy,is-1))
     :       *(w(ix,iy,is+1,2)+w(ix,iy,is,2))/ds08a(is)
           tavs=(wsig(ix,iy,is+1,2)*(pt(ix,iy,is+1,2)+pt(ix,iy,is,2))
     :       -wsig(ix,iy,is,2)*(pt(ix,iy,is,2)+pt(ix,iy,is-1,2)))
     :       /ds12(is)
           tavx=(u(ix,iy,is,2)*pp10(ix,iy,2)
     :       *(pt(ix+1,iy,is,2)+pt(ix,iy,is,2))
     :       -u(ix-1,iy,is,2)*pp10(ix-1,iy,2)
     :       *(pt(ix,iy,is,2)+pt(ix-1,iy,is,2)))/dx2
           pt(ix,iy,is,3)=ptnew(ix,iy,is)
     :       -dtl*(tavx+pp(ix,iy,2)   
     :       *(tavs-q-difunt(ix,iy,is)
     :        ))/pp(ix,iy,3) 
         enddo
         enddo
         enddo 
         
         endif

 !                 do ix = 1,nx1
 !           write(0,*) ix,w(ix,ny,ns-1,2),u(ix,ny,ns-1,2)
 !           enddo
         
	end subroutine thermo
	
	

	subroutine thermo_part2(timeho,ptm2,wrk)
	
	use alloc
      use optic_parameters                                                      VS,05.2007

      implicit real*8(a-h,o-z)
      dimension ptm2(0:nx1,0:ny1,0:ns1)
      dimension wrk(0:nx1,0:ny1,0:ns1)

!-----Input arrays into sorad---------------------
      real pl_rad(ny,ns),ta_rad(ny,ns-1),                                       VS,05.2007
     & wa_rad(ny,ns-1),zlayer(ny,ns),fcld(ny,ns-1),                             VS,05.2007
     & zdel(ny,ns-1)                                                            VS,05.2007
      real cwc(ny,ns-1,3)
!-----Output arrays into sorad---------------------           
      real flx  (ny,ns), flc  (ny,ns)                                           VS,05.2007
      REAL flx_d(ny,ns), flx_u(ny,ns)                                           VS,05.2007
      REAL flc_d(ny,ns), flc_u(ny,ns)                                           VS,05.2007
      real fdiruv  (ny), fdifuv  (ny)                                           VS,05.2007
      real fdirpar (ny), fdifpar (ny)                                           VS,05.2007
      real fdirir  (ny), fdifir  (ny)                                           VS,05.2007
!-----Input arrays to irrad------------------------                             VS,06.2007
      real fs(ny,2), tg(ny,2), tv(ny,2)                                         VS,06.2007
      real tsurfs (ny,2)                                                        VS,06.2007
!-----Output arrays from irrad---------------------                 
      real flx_lw(ny,ns), flc_lw(ny,ns)                                         VS,06.2007
      real dfdts(ny,ns), sfcem(ny)                                              VS,06.2007
      
      real, external:: sin_sun                                                  VS,05.2007
      real cosz(ny)                                                             VS,05.2007
      
      integer ict, icb                                                          VS,05.2007
       
      logical firstcall, overcast, cldwater,                                    VS,05.2007
     & overcast_lw, cldwater_lw, high, trace, vege, aerosol                     VS,06.2007
      data firstcall /.true./                                                   VS,05.2007
c
      logical sim,nao
      parameter (sim=.true.,nao=.false.)
      
!------Heating due to shortwave and longwave radiation absorption-------------         
      if (radpar==3) then                                                       VS,05.2007
       if (firstcall) then                                                      VS,05.2007
        call alloc_optics(ny,ns-1)                                              VS,05.2007
        cldwater = .true.                                                       VS,05.2007
        overcast = .true.                                                       VS,05.2007
	  cldwater_lw = .true.                                                    VS,06.2007
        overcast_lw = .true.                                                    VS,06.2007
	  high = .false.                                                          VS,06.2007
        trace = .true.                                                          VS,06.2007
        vege = .true.                                                           VS,06.2007
        aerosol = .true.                                                        VS,06.2007 
        firstcall = .false.                                                     VS,05.2007
       endif                                                                    VS,05.2007
!------Clirad subroutines SORAD and IRRAD are called once   ------
!------per nradcall time steps,                             ------
!------in between Radheat is not updated -------------------------   
       if (mod(nstep,nradcall)==0.or.nstep==1) then                             VS,05.2007
        do ix = 1, nx                                                           VS,05.2007
         cosz = sin_sun(imonth,iday,ihour,iminu,iseco,xlatit)                   VS,05.2007
         do iy = 1, ny                                                          VS,05.2007
          do is = 1, ns                                                         VS,05.2007
           pl_rad(iy,is) = 1.d-2* (sigma1(is)*pp(ix,iy,2)+ptop)                 VS,05.2007
           if (is>=2.and.iy==int(ny/2)) then                                    VS,05.2007
            if (pl_rad(iy,is)>400.and.pl_rad(iy,is-1)<400.) ict = is            VS,05.2007
            if (pl_rad(iy,is)>700.and.pl_rad(iy,is-1)<700.) icb = is            VS,05.2007
           endif                                                                VS,05.2007
           zlayer(iy,is) = (phis(ix,iy,is)+phi(ix,iy,is))/g                     VS,05.2007
           if (is==ns) zlayer(iy,is) = hsuf(ix,iy)                              VS,05.2007
!          Converting zlayer to km, as needed to input in sorad
           zlayer(iy,is) = zlayer(iy,is)/1000.                                  VS,05.2007
          enddo                                                                 VS,05.2007
         enddo                                                                  VS,05.2007
         do iy = 1, ny                                                          VS,05.2007
          do is = 1, ns-1                                                       VS,05.2007
           t1_rad = (pts(ix,iy,is)+pt(ix,iy,is,2))/                             VS,05.2007
     &      (p00/(sigma1(is)*pp(ix,iy,2)+ptop))**akapa                          VS,05.2007
           t2_rad = (pts(ix,iy,is+1)+pt(ix,iy,is+1,2))/                         VS,05.2007
     &      (p00/(sigma1(is+1)*pp(ix,iy,2)+ptop))**akapa                        VS,05.2007
           ta_rad(iy,is) = (t1_rad + t2_rad)/2.                                 VS,05.2007
           wa_rad(iy,is) = (qv(ix,iy,is,2)+qv(ix,iy,is+1,2))/2.                 VS,05.2007
           zdel(iy,is)   = zlayer(iy,is)-zlayer(iy,is+1)                        VS,05.2007
!--------Notes:
!--------Ice particles are not computed by atmospheric model----------         
!--------fcld should be calculated un atmospheric model---------------
           fcld(iy,is)   = 0.                                                   VS,05.2007
           cwc(iy,is,1)  = 0.                                                   VS,05.2007
           cwc(iy,is,2)  = 0.                                                   VS,05.2007
           cwc(iy,is,3)  = 0.                                                   VS,05.2007
           if (ifqc==1)                                                         VS,05.2007
     &      cwc(iy,is,2) = (qc(ix,iy,is,2)+qc(ix,iy,is+1,2))/2.                 VS,05.2007
           if (ifqr==1)                                                         VS,05.2007
     &      cwc(iy,is,3) = (qr(ix,iy,is,2)+qr(ix,iy,is+1,2))/2.                 VS,05.2007
          enddo                                                                 VS,05.2007
         enddo                                                                  VS,05.2007
         do iy = 1, ny                                                          VS,06.2007
          tsurfs(iy,:) = sngl(tsurf(ix,iy))                                     VS,06.2007 
         enddo                                                                  VS,06.2007
         call opticparset(ny,ns-1,alb(ix,1:ny),emis(ix,1:ny),zlayer)            VS,05.2007
        if (cosz(1)>0) then                                                     VS,11.07.2007
         call sorad     (ny,ns-1,crel,pl_rad,ta_rad,wa_rad,oa,co2,zdel,         VS,05.2007
     &                  overcast,cldwater,cwc,taucld,reff,fcld,ict,icb,         VS,05.2007
     &                  taual,ssaal,asyal,cosz,                                 VS,05.2007
     &                  rsuvbm,rsuvdf,rsirbm,rsirdf,                            VS,05.2007
     &                  flx,flc,fdiruv,fdifuv,fdirpar,fdifpar,                  VS,05.2007
     &                  fdirir,fdifir,                                          VS,05.2007
     &                  flx_d,flx_u,flc_d,flc_u)                                VS,05.2007
         if (shortwave==1) then                                                 VS,11.07.2007
          do is = 1, ns                                                         VS,05.2007
           Srad(ix,1:ny,is) = s0*flx(1:ny,is)*amax1(cosz(1:ny),0.)              VS,05.2007
          enddo                                                                 VS,05.2007
         else                                                                   VS,11.07.2007
          Srad(ix,1:ny,1:ns) = 0.                                               VS,11.07.2007
         endif                                                                  VS,11.07.2007
         Srad_surf(ix,1:ny) = flx(1:ny,ns)*                                     VS,06.2007 
!	   (fdiruv(1:ny)+fdifuv(1:ny)+fdirpar(1:ny)                               VS,05.2007
!     &   +fdifpar(1:ny)+fdirir(1:ny)+fdifir(1:ny))*                            VS,05.2007
     &   s0*amax1(cosz(1:ny),0.)                                                VS,05.2007
        else                                                                    VS,05.2007
         Srad(ix,1:ny,1:ns) = 0.                                                VS,05.2007
         Srad_surf(ix,1:ny) = 0.                                                VS,05.2007
        endif                                                                   VS,05.2007
        fs(:,1) = 0.5                                                           VS,06.2007 
        fs(:,2) = 0.5                                                           VS,06.2007 
!         print*, 'nsur, na', nsur, na
        call irrad(ny,ns-1,pl_rad,ta_rad,wa_rad,oa,                             VS,06.2007 
     &             0.5*(ta_rad(1:ny,ns-1)+tsurfs(1:ny,1)),                      VS,06.2007 
     &             co2,high,trace,n2o,ch4,cfc11,cfc12,cfc22,                    VS,06.2007 
     &             vege,nsur,fs,tsurfs(1:ny,1:2),                               VS,06.2007 
     &             eg,tsurfs(1:ny,1:2),ev,rvir,                                 VS,06.2007  
     &             overcast_lw,cldwater_lw,cwc,taucl_lw,fcld,ict,icb,           VS,06.2007 
     &             aerosol,na,taual_lw,ssaal_lw,asyal_lw,                       VS,06.2007 
     &             flx_lw,flc_lw,dfdts,sfcem)                                   VS,06.2007  
        Lrad_surf(ix,1:ny) = flx_lw(1:ny,ns)                                    VS,11.07.2007 
        if (longwave==1) then                                                   VS,11.07.2007
         Lrad(ix,1:ny,1:ns) = flx_lw(1:ny,1:ns)                                 VS,11.07.2007  
        else                                                                    VS,11.07.2007
         Lrad(ix,1:ny,1:ns) = 0.                                                VS,11.07.2007
        endif                                                                   VS,11.07.2007
        enddo                                                                   VS,05.2007
       endif                                                                    VS,05.2007
       
!       print*, 'Srad',Srad(15,15,1:ns)                                          VS,05.2007
!       print*, 'date',imonth,iday,ihour
!       read*
!       print*, 'Srad_sum', Srad_surf(15,15)
       continue
       
!      Radiational fluxes are assumed to be downward

!       Lrad = 0.                                                                VS,05.2007
       do is = 2, ns-1                                                          VS,05.2007
        do ix = 2, nx                                                           VS,05.2007
         do iy = 2, ny                                                          VS,05.2007
          p=sigma0(is)*pp(ix,iy,2)+ptop                                         VS,05.2007
          temperature=(pts(ix,iy,is)+pt(ix,iy,is,2))/(p00/p)**akapa             VS,05.2007
          Radheat(ix,iy,is) = - 0.5*(s(ix,iy,is-1)+s(ix,iy,is))*                VS,05.2007
     &     (Srad(ix,iy,is+1)-Srad(ix,iy,is-1)+                                  VS,05.2007
     &     Lrad(ix,iy,is+1)-Lrad(ix,iy,is-1))/                                  VS,05.2007
     &     ds02(is)/cp/p*temperature*r                                          VS,05.2007
           pt(ix,iy,is,3) = pt(ix,iy,is,3)+dtl*                                 VS,05.2007
     &     pp(ix,iy,2)*Radheat(ix,iy,is)/pp(ix,iy,3)                            VS,05.2007
         enddo                                                                  VS,05.2007
        enddo                                                                   VS,05.2007
       enddo                                                                    VS,05.2007
      endif                                                                     VS,05.2007

!	print*, Srad_surf
!	read*

!      print*, Lrad_surf
!	read*
      
c
c surface heat flux (bulk parametrization):
c
      if(ifsoil.ne.0.or.ifhle.eq.1.or.tfix.ne.0) then
         do iy=2,ny
            do ix=2,nx
	         p=sigma0(ns-1)*pp(ix,iy,2)+ptop
!               pt(ix,iy,ns-1,3)=pt(ix,iy,ns-1,3)
!     :            +dtl*h(ix,iy)*(1.e5/p)**akapa
!     :            /(deltaz(ix,iy)*cp*rhoa)
            enddo
         enddo
	
      elseif(cdcoef.gt.0.) then
        ctbulk=cdcoef*dtl
        do iy=2,ny
        do ix=2,nx
          p=sigma0(ns-1)*pp(ix,iy,2)+ptop
          if(xlake(ix,iy).eq.0.) then
            tsfunc=tsoil(ix,iy,timeho)
          else
            tslake(ix,iy)=tlake(ix,iy,timeho)
            tsfunc=tslake(ix,iy)
          endif
          pt(ix,iy,ns-1,3)=pt(ix,iy,ns-1,3)+ctbulk*(1.e5/p)**akapa
     :      *uvsuf(ix,iy)*(tsfunc-pts(ix,iy,ns-1)
     :      -pt(ix,iy,ns-1,2))/deltaz(ix,iy)
        enddo
        enddo
      endif
c
c horizontal boundary conditions:
c
      call ptbc(nx,ny,ns,pt(0,0,0,3),ptbx,ptby,ptcc,iobptx,iobpty
     :   ,ptbout,prt,nchn1)
!         do is=1,ns-1
 !        do ix=1,nx+1
  !       pts(ix,ny+1,is)=ini_pts(is)
  !       enddo
  !       enddo
     
c
      if(raylei) then
        do is=1,idrmax
        do iy=1,ny1
        do ix=1,nx1
          pt(ix,iy,is,3)=pt(ix,iy,is,3)-dtl/taudra(is)*pt(ix,iy,is,2)
        enddo
        enddo
        enddo
      endif
c
      do iy=1,ny1
      do ix=1,nx1
!       p=sigma0(ns-1)*pp(ix,iy,2)+ptop
!       thsurf=tsurf(ix,iy)*(p00/p)**akapa
!       pt(ix,iy,ns,3)=2.*thsurf-pt(ix,iy,ns-1,3)
!     :  -2.*pts(ix,iy,ns-1)
!      write(0,*)ix,iy,pt(ix,iy,ns,3)
        pt(ix,iy,ns,3)=pt(ix,iy,ns-1,3)
        pt(ix,iy,0,3)=pt(ix,iy,1,3)
      enddo
      enddo
c
c smoothing:
c
      if(dohsmo .and. mod(nstep,numsmo).eq.0) then
        call hsmoot2(pt(0,0,0,3),ptfd,hk1,nx1,ny1,ns,1,nx1,1,ny1,1,ns-1
     :    ,hdamp,dt)
      endif
 !     if(dovsmt .and. mod(nstep,numsmo).eq.0) then
 !       call vsmoot(pt(0,0,0,3),wrk,nx1,ny1,ns,1,nx1,1,ny1,1,ns-1
 !    :    ,vdampt)
 !     endif
c
      !if(.not.tfct) then
        call aselin(pt(0,0,0,3),pt(0,0,0,2),ptm2
     :    ,0,nx1,0,ny1,0,ns,filt)
      !endif
c
      call radbch(ptcc,ptbx,ptby,pt(0,0,0,3),pt(0,0,0,2),ptm2
     :  ,nx1,ny1,ns1,1,ns-1,1,nx1,1,ny1)
c      deallocate(ptm2)
c
c vertical boundary conditions:
c
      do iy=1,ny1
      do ix=1,nx1
 !      p=sigma0(ns-1)*pp(ix,iy,2)+ptop
 !      thsurf=tsurf(ix,iy)*(p00/p)**akapa
 !      pt(ix,iy,ns,3)=2.*thsurf-pt(ix,iy,ns-1,3)
 !    :  -2.*pts(ix,iy,ns-1)
       pt(ix,iy,ns,3)=pt(ix,iy,ns-1,3)
        pt(ix,iy,0,3)=pt(ix,iy,1,3)
      enddo
      enddo
c
      return
      end
      
      
      subroutine humid(wrk)
c-----------------------------------------------------------------------
c     time integration of specific humidity equation
c-----------------------------------------------------------------------
c
      use alloc

      implicit real*8(a-h,o-z)
      dimension wrk(0:nx1,0:ny1,0:ns1)

      if(tfct) then

         do k=1,3
         do is=0,ns
         do iy=1,ny1
         do ix=1,nx1
            qv(ix,iy,is,k)=qv(ix,iy,is,k)*pp(ix,iy,k)
         enddo
         enddo
         enddo
         enddo

         allocate(work1(0:nx1,0:ny1,0:ns1))
         allocate(work2(0:nx1,0:ny1,0:ns1))
         allocate(work3(0:nx1,0:ny1,0:ns1))
         allocate(work4(0:nx1,0:ny1,0:ns1))
         allocate(work5(0:nx1,0:ny1,0:ns1))
         allocate(work6(0:nx1,0:ny1,0:ns1))
         allocate(work7(0:nx1,0:ny1,0:ns1))
         allocate(work8(0:nx1,0:ny1,0:ns1))
         call fct2d(qv(0,0,0,3),qv(0,0,0,2),qv(0,0,0,1)
     :     ,u(0,0,0,3),u(0,0,0,2)
     :     ,v(0,0,0,3),v(0,0,0,2)
     :     ,wsig(0,0,0,3),wsig(0,0,0,2)
     :     ,nx,ny,ns
     :     ,dtxs,dtys,dtxy,dxys,dxdt,dydt,dt,ds0
     :     ,work1,work2,work3,work4,work5,work6,work7,work8,prt,grdout,
     :      iomod,fngrd)
 !        call fct3d(qv(0,0,0,3),qv(0,0,0,2),qv(0,0,0,1)
 !    :      ,u(0,0,0,3),u(0,0,0,2)
 !    :      ,v(0,0,0,3),v(0,0,0,2)
 !    :      ,wsig(0,0,0,3),wsig(0,0,0,2)
 !    :      ,nx,ny,ns
 !    :      ,dtxs,dtys,dtxy,dxys,dxdt,dydt,dt,ds0
 !    :      ,work1,work2,work3,work4,work5,work6,work7,work8)
         deallocate(work1)
         deallocate(work2)
         deallocate(work3)
         deallocate(work4)
         deallocate(work5)
         deallocate(work6)
         deallocate(work7)
         deallocate(work8)
         do k=1,3
         do is=0,ns
         do iy=1,ny1
         do ix=1,nx1
           qv(ix,iy,is,k)=qv(ix,iy,is,k)/pp(ix,iy,k)
         enddo
         enddo
         enddo
         enddo
         do is=1,ns-1
         do iy=2,ny
         do ix=2,nx
!	   if(ifqi.ne.0.) then
!           qv(ix,iy,is,3)=qv(ix,iy,is,3)+dtl*(difunq(ix,iy,is)
!     :        -cond(ix,iy,is)+evap(ix,iy,is)-vdepi(ix,iy,is)-                      DC.10,2009
!     :        vdeps(ix,iy,is)-vini(ix,iy,is))/pp(ix,iy,3)
!	   else
!	        qv(ix,iy,is,3)=qv(ix,iy,is,3)+dtl*(difunq(ix,iy,is)
!     :        -cond(ix,iy,is)+evap(ix,iy,is))/pp(ix,iy,3)
!	  endif
         tavs=(wsig(ix,iy,is+1,2)*(qv(ix,iy,is+1,2)+qv(ix,iy,is,2))
     :       -wsig(ix,iy,is,2)*(qv(ix,iy,is,2)+qv(ix,iy,is-1,2)))
     :       /ds12(is)
!         if(ifqc.eq.0) then
 !        qv(ix,iy,is,3)=qv(ix,iy,is,3)+dtl*pp(ix,iy,2)*                            ! DC.07.2012 added pp(ix,iy,2)
 !    :      (-tavs+difunq(ix,iy,is))/pp(ix,iy,3)
 !        else
         qv(ix,iy,is,3)=qv(ix,iy,is,3)+dtl*pp(ix,iy,2)*                 
     :      (-tavs+difunq(ix,iy,is)
     :        -cond(ix,iy,is)+evap(ix,iy,is))/pp(ix,iy,3)
 !        endif
         
         enddo
         enddo
         enddo
      else
        do is=1,ns-1
        do iy=2,ny
        do ix=2,nx
          tavx=(u(ix,iy,is,2)*pp10(ix,iy,2)*(qv(ix+1,iy,is,2)
     :      +qv(ix,iy,is,2))
     :      -u(ix-1,iy,is,2)*pp10(ix-1,iy,2)*(qv(ix,iy,is,2)
     :      +qv(ix-1,iy,is,2)))/dx2
          tavy=(v(ix,iy,is,2)*pp01(ix,iy,2)*(qv(ix,iy+1,is,2)
     :      +qv(ix,iy,is,2))
     :      -v(ix,iy-1,is,2)*pp01(ix,iy-1,2)*(qv(ix,iy,is,2)
     :      +qv(ix,iy-1,is,2)))/dy2
          tavs=(wsig(ix,iy,is+1,2)*(qv(ix,iy,is+1,2)+qv(ix,iy,is,2))
     :      -wsig(ix,iy,is,2)*(qv(ix,iy,is,2)+qv(ix,iy,is-1,2)))
     :      /ds12(is)
!	    if(ifqi.ne.0.) then 
!          qv(ix,iy,is,3)=(qv(ix,iy,is,1)*pp(ix,iy,1)
!     :      -dtl*(tavx+tavy+pp(ix,iy,2)
!     :      *(tavs-difunq(ix,iy,is)+cond(ix,iy,is)-evap(ix,iy,is)                    DC,10.2009
!     :      +vdepi(ix,iy,is)+vdeps(ix,iy,is)+vini(ix,iy,is))))
!     :      /pp(ix,iy,3)
!	    else
!	      qv(ix,iy,is,3)=(qv(ix,iy,is,1)*pp(ix,iy,1)
!     :      -dtl*(tavx+tavy+pp(ix,iy,2)
!     :      *(tavs-difunq(ix,iy,is)+cond(ix,iy,is)-evap(ix,iy,is) !+tavs                    DC,10.2009
!     :      )))
!     :      /pp(ix,iy,3)
           ! if(ix.eq.3.and.is.eq.35) write(0,*) qv(ix,iy,is,3)
!	    endif

            qv(ix,iy,is,3)=(qv(ix,iy,is,1)*pp(ix,iy,1)                                DC,7.2010
     :      -dtl*(tavx+tavy+pp(ix,iy,2)
     :      *(tavs-difunq(ix,iy,is)                 
     :      )))
     :      /pp(ix,iy,3)
        enddo
        enddo
        enddo
      endif

	return
	end

	subroutine humid_part2(wrk)
	use alloc

      implicit real*8(a-h,o-z)
	dimension wrk(0:nx1,0:ny1,0:ns1)
c
c
c surface water vapor flux (bulk parametrization):
c
c
c&solo
c
      if(ifsoil.ne.0.or.tfix.ne.0) then
!        do iy=2,ny
!        do ix=2,nx
!          qv(ix,iy,ns-1,3)=qv(ix,iy,ns-1,3)+dtl*le(ix,iy)
!     :      /deltaz(ix,iy)/hlat/rhoa
!        enddo
!        enddo
      elseif(cdcoef.gt.0.) then
         ctbulk=cdcoef*dtl
         do iy=2,ny
           do ix=2,nx
             p=sigma0(ns-1)*pp(ix,iy,2)+ptop
             qv(ix,iy,ns-1,3)=qv(ix,iy,ns-1,3)+ctbulk
     :         *uvsuf(ix,iy)*(qvsoil(ix,iy)-qvs(ix,iy,ns-1)
     :         -qv(ix,iy,ns-1,2))/deltaz(ix,iy)
           enddo
         enddo
      endif
c
c horizontal boundary conditions:
c
      call ptbc(nx,ny,ns,qv(0,0,0,3),qvbx,qvby,qvcc,iobptx,iobpty
     :   ,ptbout,prt,nchn1)
c
      if(raylei) then
         do is=1,idrmax
         do iy=1,ny1
         do ix=1,nx1
           qv(ix,iy,is,3)=qv(ix,iy,is,3)-dtl/taudra(is)
     :        *(qv(ix,iy,is,2)-qvs(ix,iy,is))
         enddo
         enddo
         enddo
      endif
c
      do iy=1,ny1
      do ix=1,nx1
        qv(ix,iy,ns,3)=qv(ix,iy,ns-1,3)
        qv(ix,iy,0,3)=qv(ix,iy,1,3)
      enddo
      enddo
c
      if(dohsmo .and. mod(nstep,numsmo).eq.0) then
        call hsmoot2(qv(0,0,0,3),ptfd,hk1,nx1qv,ny1qv,ns,1,nx1qv
     :   ,1,ny1qv,1,ns-1
     :    ,2*hdamp,dt)
      endif
!      if(dovsmt .and. mod(nstep,numsmo).eq.0) then
!         call vsmoot(qv(0,0,0,3),wrk,nx1,ny1,ns,1,nx1,1,ny1,1,ns-1
!     :      ,vdampt)
!      endif

!      if(.not.tfct) then
         call aselin(qv(0,0,0,3),qv(0,0,0,2),qv(0,0,0,1)
     :   ,0,nx1,0,ny1,0,ns,filt)
 !     endif

      call radbch(qvcc,qvbx,qvby,qv(0,0,0,3),qv(0,0,0,2),qv(0,0,0,1)
     :   ,nx1,ny1,ns1,1,ns-1,1,nx1,1,ny1)

      return
      end
      
      
      subroutine cloudw(wrk)
c-----------------------------------------------------------------------
c
c     time integration of the cloud water equation
c
c uses 1 wk
c-----------------------------------------------------------------------
      use alloc

      implicit real*8 (a-h,o-z)
      dimension wrk(0:nx1,0:ny1,0:ns1)

      if(tfct) then
         do k=1,3
         do is=0,ns
         do iy=1,ny1
         do ix=1,nx1
            qc(ix,iy,is,k)=qc(ix,iy,is,k)*pp(ix,iy,k)
         enddo
         enddo
         enddo
         enddo
         allocate(work1(0:nx1,0:ny1,0:ns1))
         allocate(work2(0:nx1,0:ny1,0:ns1))
         allocate(work3(0:nx1,0:ny1,0:ns1))
         allocate(work4(0:nx1,0:ny1,0:ns1))
         allocate(work5(0:nx1,0:ny1,0:ns1))
         allocate(work6(0:nx1,0:ny1,0:ns1))
         allocate(work7(0:nx1,0:ny1,0:ns1))
         allocate(work8(0:nx1,0:ny1,0:ns1))
         call fct2d(qc(0,0,0,3),qc(0,0,0,2),qc(0,0,0,1)
     :     ,u(0,0,0,3),u(0,0,0,2)
     :     ,v(0,0,0,3),v(0,0,0,2)
     :     ,wsig(0,0,0,3),wsig(0,0,0,2)
     :     ,nx,ny,ns
     :     ,dtxs,dtys,dtxy,dxys,dxdt,dydt,dt,ds0
     :     ,work1,work2,work3,work4,work5,work6,work7,work8,prt,grdout,
     :      iomod,fngrd)
!         call fct3d(qc(0,0,0,3),qc(0,0,0,2),qc(0,0,0,1)
!     :      ,u(0,0,0,3),u(0,0,0,2)
!     :      ,v(0,0,0,3),v(0,0,0,2)
!     :      ,wsig(0,0,0,3),wsig(0,0,0,2)
!     :      ,nx,ny,ns
!     :      ,dtxs,dtys,dtxy,dxys,dxdt,dydt,dt,ds0
!     :      ,work1,work2,work3,work4,work5,work6,work7,work8)
         deallocate(work1)
         deallocate(work2)
         deallocate(work3)
         deallocate(work4)
         deallocate(work5)
         deallocate(work6)
         deallocate(work7)
         deallocate(work8)

         do k=1,3
         do is=0,ns
         do iy=1,ny1
         do ix=1,nx1
           qc(ix,iy,is,k)=qc(ix,iy,is,k)/pp(ix,iy,k)
         enddo
         enddo
         enddo
         enddo
         do is=1,ns-1
         do iy=2,ny
         do ix=2,nx
         tavs=(wsig(ix,iy,is+1,2)*(qc(ix,iy,is+1,2)+qc(ix,iy,is,2))
     :       -wsig(ix,iy,is,2)*(qc(ix,iy,is,2)+qc(ix,iy,is-1,2)))
     :       /ds12(is)
!	   if (ifqi.ne.0.) then
!           qc(ix,iy,is,3)=qc(ix,iy,is,3)+
!     :        dtl*pp(ix,iy,2)*(-tavs+difunqc(ix,iy,is)                                                 DC,07.2012 added pp(ix,iy,2)
!     :       +cond(ix,iy,is)-auto(ix,iy,is)-col(ix,iy,is)+imlt(ix,iy,is)                         DC,10.2009
!     :       -hmfrz(ix,iy,is)-sacrw(ix,iy,is)-sacrwr(ix,iy,is)                                   DC,10.2009
!     :       -sbercw(ix,iy,is)-ibercw(ix,iy,is))/pp(ix,iy,3)                                     DC,10.2009
!	   else
	     qc(ix,iy,is,3)=qc(ix,iy,is,3)+
     :	     dtl*pp(ix,iy,2)*(+difunqc(ix,iy,is)                      !-tavs                            DC,07.2012 added pp(ix,iy,2)
     :       +cond(ix,iy,is)-auto(ix,iy,is)-col(ix,iy,is))/pp(ix,iy,3)                           DC,10.2009
!	   endif
!         qc(ix,iy,is,3)=qc(ix,iy,is,3)+dtl*     
!     :    pp(ix,iy,2)*difunqc(ix,iy,is)/pp(ix,iy,3)                               !-tavs
         enddo
         enddo
         enddo
      else
c      if(verbose.ge.) write(*,*) 'do 300'
        do is=1,ns-1
        do iy=2,ny
        do ix=2,nx
        tavx=(u(ix,iy,is,2)*pp10(ix,iy,2)*(qc(ix+1,iy,is,2)
     :    +qc(ix,iy,is,2))
     :    -u(ix-1,iy,is,2)*pp10(ix-1,iy,2)*(qc(ix,iy,is,2)
     :    +qc(ix-1,iy,is,2)))/dx2
        tavy=(v(ix,iy,is,2)*pp01(ix,iy,2)*(qc(ix,iy+1,is,2)
     :    +qc(ix,iy,is,2))
     :    -v(ix,iy-1,is,2)*pp01(ix,iy-1,2)*(qc(ix,iy,is,2)
     :    +qc(ix,iy-1,is,2)))/dy2
        tavs=(wsig(ix,iy,is+1,2)*(qc(ix,iy,is+1,2)+qc(ix,iy,is,2))
     :    -wsig(ix,iy,is,2)*(qc(ix,iy,is,2)+qc(ix,iy,is-1,2)))
     :    /ds12(is)
!	  if(ifqi.ne.0.) then
!        qc(ix,iy,is,3)=(qc(ix,iy,is,1)*pp(ix,iy,1)                                             DC,10.2009
!     :    -dtl*(tavx+tavy+pp(ix,iy,2)                                                          DC,10.2009
!     :    *(tavs-difunqc(ix,iy,is)-cond(ix,iy,is)+auto(ix,iy,is)                               DC,10.2009
!     :    +col(ix,iy,is)-imlt(ix,iy,is)+hmfrz(ix,iy,is)                                        DC,10.2009
!     :    +sacrw(ix,iy,is)+sacrwr(ix,iy,is)+sbercw(ix,iy,is)                                   DC,10,2009
!     :    +ibercw(ix,iy,is))))/pp(ix,iy,3)                                                     DC,10.2009
!	  else
!        qc(ix,iy,is,3)=(qc(ix,iy,is,1)*pp(ix,iy,1)                                             DC,10.2009
!     :    -dtl*(tavx+tavy+pp(ix,iy,2)                                                          DC,10.2009
!     :    *(-difunqc(ix,iy,is)-cond(ix,iy,is)+auto(ix,iy,is)    !tavs                            DC,10.2009
!     :    +col(ix,iy,is))))/pp(ix,iy,3)                                                     DC,10.2009
!	  endif
        qc(ix,iy,is,3)=(qc(ix,iy,is,1)*pp(ix,iy,1)                              DC,7.2010                                           DC,10.2009
     :    -dtl*(tavx+tavy+pp(ix,iy,2)                                                          DC,10.2009
     :    *(tavs-difunqc(ix,iy,is))))/pp(ix,iy,3)
        enddo
        enddo
        enddo
c        if(verbose.ge.) write(*,*) '300'
      endif

      return
	end

	subroutine cloudw_part2(wrk)
	use alloc

      implicit real*8 (a-h,o-z)
	dimension wrk(0:nx1,0:ny1,0:ns1)
c
c
c surface cloud water flux (bulk parametrization):
c ?????
c
c horizontal boundary conditions:
c
c      if(verbose.ge.) write(*,*) 'call ptbc'
      call ptbc(nx,ny,ns,qc(0,0,0,3),qcbx,qcby,qccc,iobptx,iobpty
     :  ,ptbout,prt,nchn1)
c
      if(raylei) then
c      if(verbose.ge.1) write(*,*) 'raylei'
        do is=1,idrmax
        do iy=1,ny1
        do ix=1,nx1
          qc(ix,iy,is,3)=qc(ix,iy,is,3)-dtl/taudra(is)*qc(ix,iy,is,2)
        enddo
        enddo
        enddo
      endif
c
      do iy=1,ny1
      do ix=1,nx1
        qc(ix,iy,ns,3)=qc(ix,iy,ns-1,3)
        qc(ix,iy,0,3)=qc(ix,iy,1,3)
      enddo
      enddo
c
 !     if(dohsmo .and. mod(nstep,numsmo).eq.0) then
c      if(verbose.ge.1) write(*,*) 'hsmoot'
 !      call hsmoot2(qc(0,0,0,3),ptfd,hk1,nx1qc,ny1qc,ns,1,nx1qc
 !    :   ,1,ny1qc,1,ns-1
 !    :    ,hdamp,dt)
 !     endif
!      if(dovsmt .and. mod(nstep,numsmo).eq.0) then
c      if(verbose.ge.) write(*,*) 'vsmoot'
!         call vsmoot(qc(0,0,0,3),wrk,nx1,ny1,ns,1,nx1,1,ny1,1,ns-1
!     :      ,vdampt)
!      endif

!      if(.not.tfct) then
c      if(verbose.ge.) write(*,*) 'aselin'
         call aselin(qc(0,0,0,3),qc(0,0,0,2),qc(0,0,0,1)
     :   ,0,nx1qc,0,ny1qc,0,ns,filt)
!      endif

c      if(verbose.ge.) write(*,*) 'radbch'
      call radbch(qccc,qcbx,qcby,qc(0,0,0,3),qc(0,0,0,2),qc(0,0,0,1)
     :   ,nx1qc,ny1qc,ns1,1,ns-1,1,nx1qc,1,ny1qc)

      return
      end
      
	subroutine cloudice(wrk)                                                  DC,10.2009
c-----------------------------------------------------------------------
c     time integration of the cloud ice equation
c-----------------------------------------------------------------------
      use alloc                                                                 DC,10.2009
	implicit none                                                             DC,10.2009
	real(8) wrk(0:nx1,0:ny1,0:ns1)                                            DC,10.2009
	integer k,ix,iy,is                                                        DC,10.2009
	real(8) tavs,tavy,tavx
	

	if(tfct) then
         do k=1,3
         do is=0,ns
         do iy=1,ny1
         do ix=1,nx1
            qci(ix,iy,is,k)=qci(ix,iy,is,k)*pp(ix,iy,k)
         enddo
         enddo
         enddo
         enddo
         allocate(work1(0:nx1,0:ny1,0:ns1))
         allocate(work2(0:nx1,0:ny1,0:ns1))
         allocate(work3(0:nx1,0:ny1,0:ns1))
         allocate(work4(0:nx1,0:ny1,0:ns1))
         allocate(work5(0:nx1,0:ny1,0:ns1))
         allocate(work6(0:nx1,0:ny1,0:ns1))
         allocate(work7(0:nx1,0:ny1,0:ns1))
         allocate(work8(0:nx1,0:ny1,0:ns1))
         call fct3d(qci(0,0,0,3),qci(0,0,0,2),qci(0,0,0,1)
     :      ,u(0,0,0,3),u(0,0,0,2)
     :      ,v(0,0,0,3),v(0,0,0,2)
     :      ,wsig(0,0,0,3),wsig(0,0,0,2)
     :      ,nx,ny,ns
     :      ,dtxs,dtys,dtxy,dxys,dxdt,dydt,dt,ds0
     :      ,work1,work2,work3,work4,work5,work6,work7,work8)
         deallocate(work1)
         deallocate(work2)
         deallocate(work3)
         deallocate(work4)
         deallocate(work5)
         deallocate(work6)
         deallocate(work7)
         deallocate(work8)

         do k=1,3
         do is=0,ns
         do iy=1,ny1
         do ix=1,nx1
           qci(ix,iy,is,k)=qci(ix,iy,is,k)/pp(ix,iy,k)
         enddo
         enddo
         enddo
         enddo
         do is=1,ns-1
         do iy=2,ny
         do ix=2,nx
           qci(ix,iy,is,3)=qci(ix,iy,is,3)+
     :      dtl*pp(ix,iy,2)*(difunqci(ix,iy,is)                                 DC,07.2012 added pp(ix,iy,2)
     :		+vini(ix,iy,is)+vdepi(ix,iy,is)+hmfrz(ix,iy,is)                         DC,11.2009
     :      +ibercw(ix,iy,is)-imlt(ix,iy,is)-sberci(ix,iy,is)                   DC,11.2009
     :      -sagg(ix,iy,is)-saci(ix,iy,is)-raci(ix,iy,is))/pp(ix,iy,3)          DC,11.2009    
!             qci(ix,iy,is,3)=qci(ix,iy,is,3)+
!     :       dtl*pp(ix,iy,2)*difunqci(ix,iy,is)/pp(ix,iy,3)
         enddo 
         enddo
         enddo  
      else                                                                      DC,11.2009
        do is=1,ns-1                                                            DC,11.2009
        do iy=2,ny                                                              DC,11.2009
        do ix=2,nx                                                              DC,11.2009
        tavx=(u(ix,iy,is,2)*pp10(ix,iy,2)*(qci(ix+1,iy,is,2)                    DC,11.2009
     :    +qci(ix,iy,is,2))                                                     DC,11.2009
     :    -u(ix-1,iy,is,2)*pp10(ix-1,iy,2)*(qci(ix,iy,is,2)                     DC,11.2009
     :    +qci(ix-1,iy,is,2)))/dx2                                              DC,11.2009
        tavy=(v(ix,iy,is,2)*pp01(ix,iy,2)*(qci(ix,iy+1,is,2)                    DC,11.2009
     :    +qci(ix,iy,is,2))                                                     DC,11.2009
     :    -v(ix,iy-1,is,2)*pp01(ix,iy-1,2)*(qci(ix,iy,is,2)                     DC,11.2009
     :    +qci(ix,iy-1,is,2)))/dy2                                              DC,11.2009
        tavs=(wsig(ix,iy,is+1,2)*(qci(ix,iy,is+1,2)+qci(ix,iy,is,2))            DC,11.2009
     :    -wsig(ix,iy,is,2)*(qci(ix,iy,is,2)+qci(ix,iy,is-1,2)))                DC,11.2009
     :    /ds12(is)                                                             DC,11.2009
        qci(ix,iy,is,3)=(qci(ix,iy,is,1)*pp(ix,iy,1)                            DC,11.2009                
     :    -dtl*(tavx+tavy+pp(ix,iy,2)                                           DC,11.2009               
     :    *(tavs-difunqci(ix,iy,is)-vini(ix,iy,is)                              DC,11.2009
     :      -vdepi(ix,iy,is)-hmfrz(ix,iy,is)-ibercw(ix,iy,is)                   DC,11.2009
     :      +imlt(ix,iy,is)+sberci(ix,iy,is)+sagg(ix,iy,is)                     DC,11.2009
     :	  +saci(ix,iy,is)+raci(ix,iy,is))))/pp(ix,iy,3)                       DC,11.2009 
!        qci(ix,iy,is,3)=(qci(ix,iy,is,1)*pp(ix,iy,1)                                          
!     :    -dtl*(tavx+tavy+pp(ix,iy,2)                                                          
!     :    *(tavs-difunqci(ix,iy,is))))/pp(ix,iy,3)                             
        enddo
        enddo
        enddo	

      endif

	return
	end

	subroutine cloudice_part2(wrk)
	use alloc
	implicit none                                                             DC,10.2009
	integer k,ix,iy,is 
	real(8):: wrk(0:nx1,0:ny1,0:ns1)


c horizontal boundary conditions:
c
c      if(verbose.ge.) write(*,*) 'call ptbc'
      call ptbc(nx,ny,ns,qci(0,0,0,3),qcibx,qciby,qcicc,iobptx,iobpty           DC,11.2009
     :  ,ptbout,prt,nchn1)                                                      DC,11.2009
c                     
      if(raylei) then                                                           DC,11.2009
c      if(verbose.ge.1) write(*,*) 'raylei'      
        do is=1,idrmax                                                          DC,11.2009
        do iy=1,ny1                                                             DC,11.2009
        do ix=1,nx1
          qci(ix,iy,is,3)=qci(ix,iy,is,3)-dtl/taudra(is)*qci(ix,iy,is,2)        DC,11.2009
        enddo                                                                   DC,11.2009
        enddo                                                                   DC,11.2009
        enddo                                                                   DC,11.2009
      endif                                                                     DC,11.2009
c 
      do iy=1,ny1                                                               DC,11.2009
      do ix=1,nx1                                                               DC,11.2009
        qci(ix,iy,ns,3)=qci(ix,iy,ns-1,3)                                       DC,11.2009
        qci(ix,iy,0,3)=qci(ix,iy,1,3)                                           DC,11.2009
      enddo                                                                     DC,11.2009
      enddo                                                                     DC,11.2009
c 
      if(dohsmo .and. mod(nstep,numsmo).eq.0) then                              DC,11.2009
c      if(verbose.ge.1) write(*,*) 'hsmoot'
         call hsmoot(qci(0,0,0,3),hk1,nx1,ny1,ns,1,nx1,1,ny1,1,ns-1             DC,11.2009
     :      ,hdamp)                                                             DC,11.2009
      endif                                                                     DC,11.2009
      if(dovsmt .and. mod(nstep,numsmo).eq.0) then                              DC,11.2009
c      if(verbose.ge.) write(*,*) 'vsmoot'        
         call vsmoot(qci(0,0,0,3),wrk,nx1,ny1,ns,1,nx1,1,ny1,1,ns-1             DC,11.2009
     :      ,vdampt)                                                            DC,11.2009
      endif                                                                     DC,11.2009
 
 !     if(.not.tfct) then                                                        DC,11.2009
c      if(verbose.ge.) write(*,*) 'aselin' 
         call aselin(qci(0,0,0,3),qci(0,0,0,2),qci(0,0,0,1)                     DC,11.2009
     :   ,0,nx1,0,ny1,0,ns,filt)                                                DC,11.2009
  !    endif                                                                     DC,11.2009

c      if(verbose.ge.) write(*,*) 'radbch'
      call radbch(qcicc,qcibx,qciby,qci(0,0,0,3),qci(0,0,0,2)                   DC,11.2009
     :	,qci(0,0,0,1),nx1,ny1,ns1,1,ns-1,1,nx1,1,ny1)                         DC,11.2009

	return                                                                    DC,10.2009
	end                                                                       DC,10.2009
      
 
                                                                                DC,10.2009
      subroutine rainw(wrk) 
c-----------------------------------------------------------------------
c     time integration of the rain water equation
c-----------------------------------------------------------------------
      use alloc

      implicit real*8(a-h,o-z)
      dimension wrk(0:nx1,0:ny1,0:ns1)

c
      if(tfct) then
        do k=1,3
        do is=0,ns
        do iy=1,ny1
        do ix=1,nx1
           qr(ix,iy,is,k)=qr(ix,iy,is,k)*pp(ix,iy,k)
        enddo
        enddo
        enddo
        enddo
        allocate(work1(0:nx1,0:ny1,0:ns1))
        allocate(work2(0:nx1,0:ny1,0:ns1))
        allocate(work3(0:nx1,0:ny1,0:ns1))
        allocate(work4(0:nx1,0:ny1,0:ns1))
        allocate(work5(0:nx1,0:ny1,0:ns1))
        allocate(work6(0:nx1,0:ny1,0:ns1))
        allocate(work7(0:nx1,0:ny1,0:ns1))
        allocate(work8(0:nx1,0:ny1,0:ns1))
        call fct3d(qr(0,0,0,3),qr(0,0,0,2),qr(0,0,0,1)
     :     ,u(0,0,0,3),u(0,0,0,2)
     :     ,v(0,0,0,3),v(0,0,0,2)
     :     ,wsig(0,0,0,3),wsig(0,0,0,2)
     :     ,nx,ny,ns
     :     ,dtxs,dtys,dtxy,dxys,dxdt,dydt,dt,ds0
     :     ,work1,work2,work3,work4,work5,work6,work7,work8)
        deallocate(work1)
        deallocate(work2)
        deallocate(work3)
        deallocate(work4)
        deallocate(work5)
        deallocate(work6)
        deallocate(work7)
        deallocate(work8)

        do k=1,3
        do is=0,ns
        do iy=1,ny1
        do ix=1,nx1
          qr(ix,iy,is,k)=qr(ix,iy,is,k)/pp(ix,iy,k)
        enddo
        enddo
        enddo
        enddo
        do is=1,ns-1
        do iy=2,ny
        do ix=2,nx
          p=sigma0(is)*pp(ix,iy,2)+ptop
          rho1=(sigma0(is-1)*pp(ix,iy,2)+ptop)/(r*(pts(ix,iy,is-1)
     :      +pt(ix,iy,is-1,2))*(p/p00)**akapa*(1+0.61*qv(ix,iy,is-1,2)
     :      -qc(ix,iy,is-1,2)))
          rho2=(sigma0(is+1)*pp(ix,iy,2)+ptop)/(r*(pts(ix,iy,is+1)
     :      +pt(ix,iy,is+1,2))*(p/p00)**akapa*(1+0.61*qv(ix,iy,is+1,2)
     :      -qc(ix,iy,is+1,2)))
	if (ifqi.ne.0.) then
          qr(ix,iy,is,3)=qr(ix,iy,is,3)+
     :      dtl*pp(ix,iy,2)*(difunqr(ix,iy,is)                                  DC,07.2012 added pp(ix,iy,2)
     :      +auto(ix,iy,is)+col(ix,iy,is)-evap(ix,iy,is)+                       DC.11,2009
     :      smlt(ix,iy,is)+sacrwr(ix,iy,is)-sacrr(ix,iy,is)                     DC,11.2009
     :      -iacr(ix,iy,is)                                                     DC,11.2009
     :      -2./(rho1+rho2)*(rho2*vrain(ix,iy,is+1)*qr(ix,iy,is+1,2)             DC,07.2012 replaced g/pp with 1/rho
     :      -rho1*vrain(ix,iy,is-1)*qr(ix,iy,is-1,2))/ds02a(is))
     :      /pp(ix,iy,3)
	 else
	    qr(ix,iy,is,3)=qr(ix,iy,is,3)+dtl*(difunqr(ix,iy,is)
     :      +auto(ix,iy,is)+col(ix,iy,is)-evap(ix,iy,is)+                       DC.11,2009
     :      -2./(rho1+rho2)*(rho2*vrain(ix,iy,is+1)*qr(ix,iy,is+1,2)            DC,07.2012 replaced g/pp with 1/rho
     :      -rho1*vrain(ix,iy,is-1)*qr(ix,iy,is-1,2))/ds02a(is))
     :      /pp(ix,iy,3)
	endif
!         qr(ix,iy,is,3)=qr(ix,iy,is,3)+
!     :      dtl*pp(ix,iy,2)*(difunqr(ix,iy,is)                      
!     :      -2./(rho1+rho2)*(rho2*vrain(ix,iy,is+1)*qr(ix,iy,is+1,2)
!     :      -rho1*vrain(ix,iy,is-1)*qr(ix,iy,is-1,2))/ds02a(is))
!     :      /pp(ix,iy,3)
        enddo
        enddo
        enddo
      else
        do is=1,ns-1
        do iy=2,ny
        do ix=2,nx
          tavx=(u(ix,iy,is,2)*pp10(ix,iy,2)*(qr(ix+1,iy,is,2)
     :      +qr(ix,iy,is,2))
     :      -u(ix-1,iy,is,2)*pp10(ix-1,iy,2)*(qr(ix,iy,is,2)
     :      +qr(ix-1,iy,is,2)))/dx2
          tavy=(v(ix,iy,is,2)*pp01(ix,iy,2)*(qr(ix,iy+1,is,2)
     :      +qr(ix,iy,is,2))
     :      -v(ix,iy-1,is,2)*pp01(ix,iy-1,2)*(qr(ix,iy,is,2)
     :      +qr(ix,iy-1,is,2)))/dy2
          tavs=(wsig(ix,iy,is+1,2)*(qr(ix,iy,is+1,2)+qr(ix,iy,is,2))
     :      -wsig(ix,iy,is,2)*(qr(ix,iy,is,2)+qr(ix,iy,is-1,2)))
     :      /ds12(is)
          p=sigma0(is)*pp(ix,iy,2)+ptop
          rho1=(sigma0(is-1)*pp(ix,iy,2)+ptop)/(r*(pts(ix,iy,is-1)
     :      +pt(ix,iy,is-1,2))*(p/p00)**akapa*(1+0.61*qv(ix,iy,is-1,2)
     :      -qc(ix,iy,is-1,2)))
          rho2=(sigma0(is+1)*pp(ix,iy,2)+ptop)/(r*(pts(ix,iy,is+1)
     :      +pt(ix,iy,is+1,2))*(p/p00)**akapa*(1+0.61*qv(ix,iy,is+1,2)
     :      -qc(ix,iy,is+1,2)))
	 if(ifqi.ne.0.) then
          qr(ix,iy,is,3)=(qr(ix,iy,is,1)*pp(ix,iy,1)
     :      -dtl*(tavx+tavy+pp(ix,iy,2)
     :      *(tavs-difunqr(ix,iy,is)-auto(ix,iy,is)-col(ix,iy,is)
     :      -smlt(ix,iy,is)-sacrwr(ix,iy,is)+sacrr(ix,iy,is)+                   DC,11.2009
     :      iacr(ix,iy,is)                                                      DC,11.2009
     :      +evap(ix,iy,is)+2./(rho1+rho2)*(rho2*vrain(ix,iy,is+1)              DC,07.2012 replaced g/pp with 1/rho
     :      *qr(ix,iy,is+1,2)-rho1*vrain(ix,iy,is-1)*qr(ix,iy,is-1,2))
     :      /ds02a(is))))/pp(ix,iy,3)
	 else
	qr(ix,iy,is,3)=(qr(ix,iy,is,1)*pp(ix,iy,1)
     :      -dtl*(tavx+tavy+pp(ix,iy,2)
     :      *(tavs-difunqr(ix,iy,is)-auto(ix,iy,is)-col(ix,iy,is)
     :      +evap(ix,iy,is)+2./(rho1+rho2)*(rho2*vrain(ix,iy,is+1)              DC,07.2012 replaced g/pp with 1/rho
     :      *qr(ix,iy,is+1,2)-rho1*vrain(ix,iy,is-1)*qr(ix,iy,is-1,2))
     :      /ds02a(is))))/pp(ix,iy,3)
	endif
!        qr(ix,iy,is,3)=(qr(ix,iy,is,1)*pp(ix,iy,1)
!     :      -dtl*(tavx+tavy+pp(ix,iy,2)
!     :      *(tavs-difunqr(ix,iy,is)
!     :      +g/pp(ix,iy,2)*(rho2*vrain(ix,iy,is+1)
!     :      *qr(ix,iy,is+1,2)-rho1*vrain(ix,iy,is-1)*qr(ix,iy,is-1,2))
!     :      /ds02a(is))))/pp(ix,iy,3)
        enddo
        enddo
        enddo
      endif

	return
	end

	subroutine rainw_part2(wrk)
	use alloc
	implicit none
	integer k,ix,iy,is 
	real(8):: wrk(0:nx1,0:ny1,0:ns1) 


c
c horizontal boundary conditions:
c
      call ptbc(nx,ny,ns,qr(0,0,0,3),qrbx,qrby,qrcc,iobptx,iobpty
     :   ,ptbout,prt,nchn1)
c
      if(raylei) then
        do is=1,idrmax
        do iy=1,ny1
        do ix=1,nx1
          qr(ix,iy,is,3)=qr(ix,iy,is,3)-dtl/taudra(is)*qr(ix,iy,is,2)
        enddo
        enddo
        enddo
      endif
c
      do iy=1,ny1
      do ix=1,nx1
        qr(ix,iy,ns,3)=qr(ix,iy,ns-1,3)
        qr(ix,iy,0,3)=qr(ix,iy,1,3)
      enddo
      enddo
c
      if(dohsmo .and. mod(nstep,numsmo).eq.0) then
        call hsmoot(qr(0,0,0,3),hk1,nx1,ny1,ns,1,nx1,1,ny1,1,ns-1
     :    ,hdamp)
      endif
      if(dovsmt .and. mod(nstep,numsmo).eq.0) then
        call vsmoot(qr(0,0,0,3),wrk,nx1,ny1,ns,1,nx1,1,ny1,1,ns-1
     :    ,vdampt)
      endif

 !     if(.not.tfct) then
        call aselin(qr(0,0,0,3),qr(0,0,0,2),qr(0,0,0,1)
     :    ,0,nx1,0,ny1,0,ns,filt)
 !     endif

      call radbch(qrcc,qrbx,qrby,qr(0,0,0,3),qr(0,0,0,2),qr(0,0,0,1)
     :  ,nx1,ny1,ns1,1,ns-1,1,nx1,1,ny1)

      return
      end

	subroutine snow(wrk) 
c-----------------------------------------------------------------------
c     time integration of the snow equation
c-----------------------------------------------------------------------
      use alloc

      implicit none                                                             DC,10.2009
	real(8) wrk(0:nx1,0:ny1,0:ns1)                                            DC,10.2009
	integer k,ix,iy,is                                                        DC,10.2009
	real(8) tavs,tavy,tavx, rho1,rho2,p                                       DC,11.2009

c
      if(tfct) then                                                             DC,11.2009
        do k=1,3
        do is=0,ns                                                              DC,11.2009
        do iy=1,ny1                                                             DC,11.2009
        do ix=1,nx1                                                             DC,11.2009
           qsn(ix,iy,is,k)=qsn(ix,iy,is,k)*pp(ix,iy,k)                          DC,11.2009
        enddo                                                                   DC,11.2009
        enddo                                                                   DC,11.2009
        enddo                                                                   DC,11.2009
        enddo                                                                   DC,11.2009
        allocate(work1(0:nx1,0:ny1,0:ns1))                                      DC,11.2009
        allocate(work2(0:nx1,0:ny1,0:ns1))                                      DC,11.2009
        allocate(work3(0:nx1,0:ny1,0:ns1))                                      DC,11.2009
        allocate(work4(0:nx1,0:ny1,0:ns1))                                      DC,11.2009
        allocate(work5(0:nx1,0:ny1,0:ns1))                                      DC,11.2009
        allocate(work6(0:nx1,0:ny1,0:ns1))                                      DC,11.2009
        allocate(work7(0:nx1,0:ny1,0:ns1))                                      DC,11.2009
        allocate(work8(0:nx1,0:ny1,0:ns1))                                      DC,11.2009
        call fct3d(qsn(0,0,0,3),qsn(0,0,0,2),qsn(0,0,0,1)                       DC,11.2009
     :     ,u(0,0,0,3),u(0,0,0,2)                                               DC,11.2009
     :     ,v(0,0,0,3),v(0,0,0,2)                                               DC,11.2009
     :     ,wsig(0,0,0,3),wsig(0,0,0,2)                                         DC,11.2009
     :     ,nx,ny,ns                                                            DC,11.2009
     :     ,dtxs,dtys,dtxy,dxys,dxdt,dydt,dt,ds0                                DC,11.2009
     :     ,work1,work2,work3,work4,work5,work6,work7,work8)                    DC,11.2009
        deallocate(work1)                                                       DC,11.2009
        deallocate(work2)                                                       DC,11.2009
        deallocate(work3)                                                       DC,11.2009
        deallocate(work4)                                                       DC,11.2009
        deallocate(work5)                                                       DC,11.2009
        deallocate(work6)                                                       DC,11.2009
        deallocate(work7)                                                       DC,11.2009
        deallocate(work8)                                                       DC,11.2009
                             
        do k=1,3                                                                DC,11.2009
        do is=0,ns                                                              DC,11.2009
        do iy=1,ny1                                                             DC,11.2009
        do ix=1,nx1                                                             DC,11.2009
          qsn(ix,iy,is,k)=qsn(ix,iy,is,k)/pp(ix,iy,k)                           DC,11.2009
        enddo                                                                   DC,11.2009
        enddo                                                                   DC,11.2009
        enddo                                                                   DC,11.2009
        enddo                                                                   DC,11.2009
        do is=1,ns-1                                                            DC,11.2009
        do iy=2,ny                                                              DC,11.2009
        do ix=2,nx                                                              DC,11.2009
          p=sigma0(is)*pp(ix,iy,2)+ptop                                         DC,11.2009
          rho1=(sigma0(is-1)*pp(ix,iy,2)+ptop)/(r*(pts(ix,iy,is-1)              DC,11.2009
     :      +pt(ix,iy,is-1,2))*(p/p00)**akapa*(1+0.61*qv(ix,iy,is-1,2)          DC,11.2009
     :      -qc(ix,iy,is-1,2)))                                                 DC,11.2009
          rho2=(sigma0(is+1)*pp(ix,iy,2)+ptop)/(r*(pts(ix,iy,is+1)              DC,11.2009
     :      +pt(ix,iy,is+1,2))*(p/p00)**akapa*(1+0.61*qv(ix,iy,is+1,2)          DC,11.2009
     :      -qc(ix,iy,is+1,2)))                                                 DC,11.2009
          qsn(ix,iy,is,3)=qsn(ix,iy,is,3)+
     :        dtl*pp(ix,iy,2)*(difunqsn(ix,iy,is)+                              DC,07.2012 added pp
     :      vdeps(ix,iy,is)+sbercw(ix,iy,is)+sberci(ix,iy,is)+                  DC,11.2009
     :      sacrw(ix,iy,is)+sagg(ix,iy,is)+saci(ix,iy,is)+                      DC,11.2009
     :      raci(ix,iy,is)+iacr(ix,iy,is)+sacrr(ix,iy,is)-                      DC,11.2009
     :      smlt(ix,iy,is)                                                      DC,11.2009
     :      -2./(rho1+rho2)*(rho2*vsnow(ix,iy,is+1)*qsn(ix,iy,is+1,2)            DC,07.2012 replaced g/pp with 1/rho
     :      -rho1*vsnow(ix,iy,is-1)*qsn(ix,iy,is-1,2))/ds02a(is))               DC,11.2009
     :      /pp(ix,iy,3)                                                        DC,11.2009
!        qsn(ix,iy,is,3)=qsn(ix,iy,is,3)+dtl*(difunqsn(ix,iy,is)                              DC,11.2009
!     :      -2./(rho1+rho2)*(rho2*vsnow(ix,iy,is+1)*qsn(ix,iy,is+1,2)            DC,11.2009
!     :      -rho1*vsnow(ix,iy,is-1)*qsn(ix,iy,is-1,2))/ds02a(is))               DC,11.2009
!     :      /pp(ix,iy,3)
	  enddo                                                                   DC,11.2009
        enddo                                                                   DC,11.2009
        enddo                                                                   DC,11.2009
      else                                                                      DC,11.2009
        do is=1,ns-1                                                            DC,11.2009
        do iy=2,ny                                                              DC,11.2009
        do ix=2,nx                                                              DC,11.2009
          tavx=(u(ix,iy,is,2)*pp10(ix,iy,2)*(qsn(ix+1,iy,is,2)                  DC,11.2009
     :      +qsn(ix,iy,is,2))                                                   DC,11.2009
     :      -u(ix-1,iy,is,2)*pp10(ix-1,iy,2)*(qsn(ix,iy,is,2)                   DC,11.2009
     :      +qsn(ix-1,iy,is,2)))/dx2                                            DC,11.2009
          tavy=(v(ix,iy,is,2)*pp01(ix,iy,2)*(qsn(ix,iy+1,is,2)                  DC,11.2009
     :      +qsn(ix,iy,is,2))                                                   DC,11.2009
     :      -v(ix,iy-1,is,2)*pp01(ix,iy-1,2)*(qsn(ix,iy,is,2)                   DC,11.2009
     :      +qsn(ix,iy-1,is,2)))/dy2                                            DC,11.2009
          tavs=(wsig(ix,iy,is+1,2)*(qsn(ix,iy,is+1,2)+qsn(ix,iy,is,2))          DC,11.2009
     :      -wsig(ix,iy,is,2)*(qsn(ix,iy,is,2)+qsn(ix,iy,is-1,2)))              DC,11.2009
     :      /ds12(is)                                                           DC,11.2009
          p=sigma0(is)*pp(ix,iy,2)+ptop                                         DC,11.2009
          rho1=(sigma0(is-1)*pp(ix,iy,2)+ptop)/(r*(pts(ix,iy,is-1)              DC,11.2009
     :      +pt(ix,iy,is-1,2))*(p/p00)**akapa*(1+0.61*qv(ix,iy,is-1,2)          DC,11.2009 
     :      -qc(ix,iy,is-1,2)))                                                 DC,11.2009
          rho2=(sigma0(is+1)*pp(ix,iy,2)+ptop)/(r*(pts(ix,iy,is+1)              DC,11.2009
     :      +pt(ix,iy,is+1,2))*(p/p00)**akapa*(1+0.61*qv(ix,iy,is+1,2)          DC,11.2009
     :      -qc(ix,iy,is+1,2)))                                                 DC,11.2009
          qsn(ix,iy,is,3)=(qsn(ix,iy,is,1)*pp(ix,iy,1)                          DC,11.2009
     :      -dtl*(tavx+tavy+pp(ix,iy,2)                                         DC,11.2009
     :      *(tavs-difunqsn(ix,iy,is)-                                          DC,11.2009
     :      vdeps(ix,iy,is)-sbercw(ix,iy,is)-sberci(ix,iy,is)-                  DC,11.2009
     :      sacrw(ix,iy,is)-sagg(ix,iy,is)-saci(ix,iy,is)                       DC,11.2009
     :      -raci(ix,iy,is)-iacr(ix,iy,is)-sacrr(ix,iy,is)+                     DC,11.2009
     :      smlt(ix,iy,is)                                                      DC,11.2009
     :      +2./(rho1+rho2)*(rho2*vsnow(ix,iy,is+1)                              DC,07.2012 replaced g/pp with 1/rho
     :      *qsn(ix,iy,is+1,2)-rho1*vsnow(ix,iy,is-1)*qsn(ix,iy,is-1,2))        DC,11.2009
     :      /ds02a(is))))/pp(ix,iy,3)                                           DC,11.2009
!        qsn(ix,iy,is,3)=(qsn(ix,iy,is,1)*pp(ix,iy,1)                             DC,11.2009
!     :      -dtl*(tavx+tavy+pp(ix,iy,2)                                         DC,11.2009
!     :      *(tavs-difunqsn(ix,iy,is)                                           DC,11.2009
!     :      +2./(rho1+rho2)*(rho2*vsnow(ix,iy,is+1)                              DC,11.2009
!     :      *qsn(ix,iy,is+1,2)-rho1*vsnow(ix,iy,is-1)*qsn(ix,iy,is-1,2))        DC,11.2009
!     :      /ds02a(is))))/pp(ix,iy,3)
	  enddo                                                                   DC,11.2009
        enddo                                                                   DC,11.2009
        enddo                                                                   DC,11.2009
      endif                                                                     DC,11.2009

	return
	end

	subroutine snow_part2(wrk)
	use alloc
	implicit none
	integer k,ix,iy,is 
	real(8):: wrk(0:nx1,0:ny1,0:ns1) 
                                   
 
c
c horizontal boundary conditions:
c
      call ptbc(nx,ny,ns,qsn(0,0,0,3),qsnbx,qsnby,qsncc,iobptx,iobpty           DC,11.2009
     :   ,ptbout,prt,nchn1)                                                     DC,11.2009
c
      if(raylei) then                                                           DC,11.2009
        do is=1,idrmax                                                          DC,11.2009
        do iy=1,ny1                                                             DC,11.2009
        do ix=1,nx1                                                             DC,11.2009
          qsn(ix,iy,is,3)=qsn(ix,iy,is,3)-dtl/taudra(is)*qsn(ix,iy,is,2)        DC,11.2009
        enddo                                                                   DC,11.2009
        enddo                                                                   DC,11.2009
        enddo                                                                   DC,11.2009
      endif                                                                     DC,11.2009
c 
      do iy=1,ny1                                                               DC,11.2009
      do ix=1,nx1                                                               DC,11.2009
        qsn(ix,iy,ns,3)=qsn(ix,iy,ns-1,3)                                       DC,11.2009
        qsn(ix,iy,0,3)=qsn(ix,iy,1,3)                                           DC,11.2009
      enddo                                                                     DC,11.2009
      enddo                                                                     DC,11.2009
c 
      if(dohsmo .and. mod(nstep,numsmo).eq.0) then                              DC,11.2009
        call hsmoot(qsn(0,0,0,3),hk1,nx1,ny1,ns,1,nx1,1,ny1,1,ns-1              DC,11.2009
     :    ,hdamp)                                                               DC,11.2009
      endif                                                                     DC,11.2009
      if(dovsmt .and. mod(nstep,numsmo).eq.0) then                              DC,11.2009
        call vsmoot(qsn(0,0,0,3),wrk,nx1,ny1,ns,1,nx1,1,ny1,1,ns-1              DC,11.2009
     :    ,vdampt)                                                              DC,11.2009
      endif                                                                     DC,11.2009

!      if(.not.tfct) then                                                        DC,11.2009
        call aselin(qsn(0,0,0,3),qsn(0,0,0,2),qsn(0,0,0,1)                      DC,11.2009
     :    ,0,nx1,0,ny1,0,ns,filt)                                               DC,11.2009
 !     endif                                                                     DC,11.2009

      call radbch(qsncc,qsnbx,qsnby,qsn(0,0,0,3),qsn(0,0,0,2),                  DC,11.2009
     :   qsn(0,0,0,1),nx1,ny1,ns1,1,ns-1,1,nx1,1,ny1)                           DC,11.2009

      return                                                                    DC,11.2009
      end                                                                       DC,11.2009
