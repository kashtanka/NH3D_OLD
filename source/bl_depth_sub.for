	subroutine bl_depth(fngrd,nonloc_tm86,difloc)
	use alloc
	implicit none
	integer ix,iy,is
	real,allocatable:: z_s(:,:,:)  ! ,thetas(:,:)
	real(8) Ri_cr,t,p,b,c1,thv(nx1,ny1),hbl2,hflmin(nx1,ny1)
     :,heatflux,th_h(nx1,ny1)
	integer nn,ms(nx1,ny1)
	logical nonloc_tm86,difloc
	character*80 fngrd
	real ws_m,th_m,wth_h,w_m3
	real,parameter:: karm=0.4,b_t=46.,A=4.5,B2=5.
	allocate (z_s(nx1,ny1,ns))
!	allocate (thetas(nx1,ny1))
        b=8.5 
	Ri_cr=0.25
        c1=0.28
	nn=1
	ms=0
	
	do is=1,ns
	do iy=2,ny
	do ix=2,nx

!	z_s(ix,iy,is)=0.5*(phis(ix,iy,is)+phi(ix,iy,is)
!     : + phis(ix,iy,is-1)+phi(ix,iy,is-1))/g

	   z_s(ix,iy,is)=(phi(ix,iy,is)+phis(ix,iy,is))/g
	enddo
	enddo
	enddo


	do iy=2,ny
	do ix=2,nx
	p=sigma0(ns)*pp(ix,iy,2)+ptop
	t=(pts(ix,iy,ns-1)+pt(ix,iy,ns-1,2))*(p/p00)**akapa
	
	if(qif.gt.0) then
	  thv(ix,iy)=(pts(ix,iy,ns-1)+pt(ix,iy,ns-1,2))
     :	  *(1+0.61*qv(ix,iy,ns-1,2))
	if(0.5*(t+tsfix(ix,iy)).gt.273.15) then
          Fv(ix,iy)=-ust_s(ix,iy)*tst_s(ix,iy)
     :    +0.61*t*le(ix,iy)/hlat
	else
          Fv(ix,iy)=-ust_s(ix,iy)*tst_s(ix,iy)
     :       +0.61*t*le(ix,iy)/hsub
	endif
	else
          thv(ix,iy)=(pts(ix,iy,ns-1)+pt(ix,iy,ns-1,2))
	  Fv(ix,iy)=-ust_s(ix,iy)*tst_s(ix,iy)
	endif

	
	
	 if(Fv(ix,iy).gt.0) then
	wstar(ix,iy)=(g/thv(ix,iy)*Fv(ix,iy)*hbl(ix,iy))**(1./3.)
	wm(ix,iy)=(ust_s(ix,iy)**3.+7.*deltaz(ix,iy)/hbl(ix,iy)*0.4
     :	*wstar(ix,iy)**3.)**(1./3.)
	!thetas(ix,iy)=thv(ix,iy)+b*Fv(ix,iy)/wm(ix,iy)
	 else
	!thetas(ix,iy)=thv(ix,iy)
	    wstar(ix,iy)=0.
	endif
	enddo
	enddo

	do is=2,ns-1
	do iy=2,ny
	do ix=2,nx
	thgrad(ix,iy,is)=-(pt(ix,iy,is+1,2)
     :    +pts(ix,iy,is+1)-pt(ix,iy,is-1,2)-pts(ix,iy,is-1))/
     :   (z_s(ix,iy,is-1)-z_s(ix,iy,is+1))
	enddo
	enddo
!	write(0,*) is, pts(3,30,is),thgrad(3,30,is)
	enddo

	do is=ns-2,2,-1
	do iy=2,ny
	do ix=2,nx
      if (thgrad(ix,iy,is).gt.0.003.and.
     :   thgrad(ix,iy,is+1).le.0.003.and.nn.eq.1) then
      	hbl(ix,iy) = max(50.,z_s(ix,iy,is)-(thgrad(ix,iy,is)-0.003)/
     :	(thgrad(ix,iy,is)-thgrad(ix,iy,is+1))
     :    *(z_s(ix,iy,is)-z_s(ix,iy,is+1)))
	ms(ix,iy)=min(ns-2,int(is+0.5*(ns-is)))
!	write(0,*)ms(ix,iy),pt(ix,iy,ms(ix,iy),2)+pts(ix,iy,ms(ix,iy))
!       hbl(ix,iy)=max(50.,z_s(ix,iy,is))

	 nn=0
	endif
	enddo
	nn=1
        enddo
	enddo

!	write(0,*) thgrad
!	write(0,*) hbl
!	stop

	do iy=2,ny
	do ix=2,nx
	   if(Fv(ix,iy).gt.0) then
	  w_m3=wstar(ix,iy)**3.+B2*ust_s(ix,iy)**3.
	  wth_h=-A*w_m3/hbl(ix,iy)
	  ws_m=(ust_s(ix,iy)**3.
     :    +7.*karm*wstar(ix,iy)**3.*0.5)**(1./3.)
	  th_m=min(2.,-b_t/ws_m*wth_h)
	  th_h(ix,iy)=pt(ix,iy,ms(ix,iy),2)+pts(ix,iy,ms(ix,iy))+th_m
!	  write(0,*) th_h(ix,iy),th_m
	  endif
	enddo
	enddo

        do is=ns,2,-1
	do iy=2,ny
	do ix=2,nx
	    if(h(ix,iy).gt.20) then
	if(pt(ix,iy,is+1,2)+pts(ix,iy,is+1).lt.th_h(ix,iy)
     :  .and.pt(ix,iy,is,2)+pts(ix,iy,is).ge.th_h(ix,iy)) then
          hbl(ix,iy)=z_s(ix,iy,is+1)+(z_s(ix,iy,is)-z_s(ix,iy,is+1))
     :   *(th_h(ix,iy)-pt(ix,iy,is+1,2)-pts(ix,iy,is+1))
     :  /(pt(ix,iy,is,2)+pts(ix,iy,is)-pt(ix,iy,is+1,2)-pts(ix,iy,is+1))	
        endif
	endif
	enddo
!	 write(0,*)iy, hbl(3,iy)
	enddo
	enddo

 !     write(0,*) Fv(18,2),Fv(18,3),Fv(18,4)

	call extrpp(nx1,ny1,hbl,1,nx1,1,ny1)
      
	
        if (nonloc_tm86) then
        do iy=1,ny
	do ix=1,nx
	if(Fv(ix,iy).gt.0) then
        gammah(ix,iy)=6.5*h(ix,iy)/cp/rhoa/wm(ix,iy)/hbl(ix,iy)
	endif
	enddo
	enddo
	endif

	if(prt.and.driout) then
          if(grdout) then
            call wgrids(hbl,'bl',2,nx,2,ny,0,0,1,0,0,0,fngrd)
	  endif
	endif
        deallocate (z_s)
!	deallocate (thetas)
	return
	end
