	
	module surfprop
!SURFACE CHARACTERISITICS FOR DIFFERENT LAND USE CATHEGORIES
!0 - 
!1 - water
!2 -
!3 -
!4 -
!5 - 
!6 -
      real(4), dimension (0:6)::veg_,wr_,iclay_,w2_,wg_,irr_,xlai_,top_
     &   ,alb_,isa_,z0_,ive_,xdd_,ts_,t2_,rsm_,tsw_,z0w_,tsm_

	contains
	SUBROUTINE surfprop_init

!0: STANDARD SURFACE

      iclay_(0)=30.
      veg_(0)=0.5
      alb_(0)=0.2
      z0_(0)=0.2
      z0w_(0)=0.01
      ts_(0)=291.
      tsm_(0)=291.
      tsw_(0)=291.
      t2_(0)=291.
      wr_(0)=0.
      wg_(0)=.30
      w2_(0)=.30
      rsm_(0)=30
      irr_(0)=1.
      ive_(0)=1
      xdd_(0)=.8

!1: BARE SOIL

	iclay_(1)=20.
      veg_(1)=0.
      alb_(1)=0.4
      z0_(1)=0.01
      z0w_(1)=0.01
      ts_(1)=300.
      tsm_(1)=300.
      tsw_(1)=300.
      t2_(1)=300.
      wr_(1)=0.
      wg_(1)=.12
      w2_(1)=.12
      rsm_(1)=30.
      irr_(1)=1.
      ive_(1)=1
      xdd_(1)=.8

!2: FOREST

	iclay_(2)=30.
      veg_(2)=1.
      alb_(2)=0.2
      z0_(2)=0.2
      z0w_(2)=0.01
      ts_(2)=291.
      tsm_(2)=291.
      tsw_(2)=291.
      t2_(2)=291.
      wr_(2)=0.
      wg_(2)=.30
      w2_(2)=.30
      rsm_(2)=30
      irr_(2)=1.
      ive_(2)=1
      xdd_(2)=.8

	END SUBROUTINE surfprop_init
	end module surfprop

      program masksoil
 	use surfprop
c
c prepares surface grds (including top)
c
      implicit double precision(a-h,o-z)
      parameter(nx=37,ny=37)
      dimension xlake(0:nx+1,0:ny+1),veg(0:nx+1,0:ny+1)
     :   ,wr(0:nx+1,0:ny+1),iclay(0:nx+1,0:ny+1)
     :   ,w2(0:nx+1,0:ny+1),wg(0:nx+1,0:ny+1),irr(0:nx+1,0:ny+1)
     :   ,xlai(0:nx+1,0:ny+1),top(0:nx+1,0:ny+1)
     :   ,alb(0:nx+1,0:ny+1),isa(0:nx+1,0:ny+1)
     :   ,z0(0:nx+1,0:ny+1),ive(0:nx+1,0:ny+1),xdd(0:nx+1,0:ny+1)
     :   ,ts(0:nx+1,0:ny+1),t2(0:nx+1,0:ny+1),rsm(0:nx+1,0:ny+1)
     :   ,tsw(0:nx+1,0:ny+1),z0w(0:nx+1,0:ny+1)
     :   ,tsm(0:nx+1,0:ny+1),dep(0:nx+1,0:ny+1)
      character*30 fnmap,fext,fnunm,fname
      parameter(wgsat=0.33,w2sat=0.35)
	real(4) dx,dy
	real(4) xx
      real(4), dimension(0:6,0:nx+1,0:ny+1):: landuse
	real(4)  landuse_w(0:6,0:2*nx,0:2*ny), top_w(0:2*nx,0:2*ny)
	integer ii,jj,ilake,itop,i0,j0,ellx,elly,surftyp

      fnmap='map/lake_map'

      
	dx=4400 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	dy=6000 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!dx=1500 !Это для сетки 0,5 на 0,5 градусов!
	!dy=1500


	ilake=1
c ilake=1 earth surface with lakes
c ilake=0 no lakes
	itop=0
c itop=1 relief is on
c itop=0 no relief (surface height = 0)

!_________________________________
!  surftyp   | surface landscape  |
!____________|____________________|
!     0      |     "standard"     |
!     1      |      bare soil     |
!     2      |       forest       |
!____________|____________________|

      surftyp = 1
      xlake=0.
	dep=-1.
	top=0.

      call surfprop_init
      landuse(:,:,:)=0.;landuse(surftyp,:,:)=1.
	


!READING FILE WITH TOPOGRAPHY
      open(UNIT=111, FILE='E:\Practice2008\Modell_gelend_august_2008\DEM
     &\Gelend_prorej.oro', status='old')
      if (itop==1) then
       
      do i=1,39  
	 do j=1,39  
	  read(111,*) ii, jj, top(i-1,j-1)
	  
	  if (top(i-1,j-1)==-200) then !!! Здесь было == 0
         xlake(i-1,j-1)=1
	   dep(i-1,j-1)= 200
	  else
	   xlake(i-1,j-1)=0
	   dep(i-1,j-1)= -1

	  endif
       enddo
      enddo
      
      else

      do i=1,39  
	 do j=1,39  
	  read(111,*) ii, jj, top(i-1,j-1)
	  if (top(i-1,j-1)==-200) then !!! Здесь было == 0
         xlake(i-1,j-1)=1
	   dep(i-1,j-1)= 200
         top(i-1,j-1) = 0
	  else
	   xlake(i-1,j-1)=0
	   dep(i-1,j-1)= -1
         top(i-1,j-1) = 0

	  endif
       enddo
      enddo


	!goto 1
	endif
 
1      call wrigrd(xlake,xmi,xma,ymi,yma,vmi,vma,dx,dy
     :   ,0,nx+1,0,ny+1,0,nx+1,0,ny+1,fext(fnmap,'lak'))
	call wrigrd(dep,xmi,xma,ymi,yma,vmi,vma,dx,dy
     :   ,0,nx+1,0,ny+1,0,nx+1,0,ny+1,fext(fnmap,'dep'))
c
      albs=0.2
      albw=0.05
      do j=0,ny+1
         do i=0,nx+1
            xlai(i,j)=2.
            isa(i,j)=35
c           top=0.
            if(xlake(i,j).eq.1) then
               iclay(i,j)=-2
               veg(i,j)=0.
               alb(i,j)=0.05
               z0w(i,j)=0.01
               z0(i,j)=0.01
               ts(i,j)=291.
               tsm(i,j)=291.
               tsw(i,j)=291.
               t2(i,j)=291.
               wr(i,j)=0.
               wg(i,j)=.15
               w2(i,j)=.15
               rsm(i,j)=30
               irr(i,j)=0
               ive(i,j)=1
               xdd(i,j)=.8
             else
               iclay(i,j)=sum(landuse(:,i,j)*iclay_(:))
               veg(i,j)=  sum(landuse(:,i,j)*veg_(:))
               alb(i,j)=  sum(landuse(:,i,j)*alb_(:))
               z0(i,j)=   sum(landuse(:,i,j)*z0_(:))
               z0w(i,j)=  sum(landuse(:,i,j)*z0w_(:))
               ts(i,j)=   sum(landuse(:,i,j)*ts_(:))
               tsm(i,j)=  sum(landuse(:,i,j)*tsm_(:))
               tsw(i,j)=  sum(landuse(:,i,j)*tsw_(:))
               t2(i,j)=   sum(landuse(:,i,j)*t2_(:))
               wr(i,j)=   sum(landuse(:,i,j)*wr_(:))
               wg(i,j)=   sum(landuse(:,i,j)*wg_(:))
               w2(i,j)=   sum(landuse(:,i,j)*w2_(:))
               rsm(i,j)=  sum(landuse(:,i,j)*rsm_(:))
               irr(i,j)=  sum(landuse(:,i,j)*irr_(:))
               ive(i,j)=  sum(landuse(:,i,j)*ive_(:))
               xdd(i,j)=  sum(landuse(:,i,j)*xdd_(:))
            endif
         enddo
      enddo

!	16 WET SOIL PATCHES
!	do ii=1,7,2
!	 do i=ii*4,(ii+1)*4
!	  do jj = 1,7,2
!	   do j=jj*4,(jj+1)*4
!	    veg(i,j)=0.
!	    w2(i,j)=.30
!          wg(i,j)=.30
!         enddo
!	  enddo
!       enddo
!	enddo   

      call iwrigrd(iclay,xmi,xma,ymi,yma,icmi,icma,dx,dy
     :   ,0,nx+1,0,ny+1,0,nx+1,0,ny+1,fext(fnmap,'icl'))
      call wrigrd(veg,xmi,xma,ymi,yma,vmi,vma,dx,dy
     :   ,0,nx+1,0,ny+1,0,nx+1,0,ny+1,fext(fnmap,'veg'))
      call wrigrd(w2,xmi,xma,ymi,yma,w2mi,w2ma,dx,dy
     :   ,0,nx+1,0,ny+1,0,nx+1,0,ny+1,fext(fnmap,'w2'))
      call wrigrd(wg,xmi,xma,ymi,yma,wgmi,wgma,dx,dy
     :   ,0,nx+1,0,ny+1,0,nx+1,0,ny+1,fext(fnmap,'wg'))
      call wrigrd(xlai,xmi,xma,ymi,yma,wgmi,wgma,dx,dy
     :   ,0,nx+1,0,ny+1,0,nx+1,0,ny+1,fext(fnmap,'lai'))
      call wrigrd(top,xmi,xma,ymi,yma,wgmi,wgma,dx,dy
     :   ,0,nx+1,0,ny+1,0,nx+1,0,ny+1,fext(fnmap,'top'))
      call wrigrd(alb,xmi,xma,ymi,yma,albmi,albma,dx,dy
     :   ,0,nx+1,0,ny+1,0,nx+1,0,ny+1,fext(fnmap,'alb'))
      call wrigrd(z0,xmi,xma,ymi,yma,z0mi,z0ma,dx,dy
     :   ,0,nx+1,0,ny+1,0,nx+1,0,ny+1,fext(fnmap,'z0'))
      call wrigrd(z0w,xmi,xma,ymi,yma,z0wmi,z0wma,dx,dy
     :   ,0,nx+1,0,ny+1,0,nx+1,0,ny+1,fext(fnmap,'z0w'))
      call wrigrd(ts,xmi,xma,ymi,yma,tsmi,tsma,dx,dy
     :   ,0,nx+1,0,ny+1,0,nx+1,0,ny+1,fext(fnmap,'ts'))
      call wrigrd(tsm,xmi,xma,ymi,yma,tsmi,tsma,dx,dy
     :   ,0,nx+1,0,ny+1,0,nx+1,0,ny+1,fext(fnmap,'tsm'))
      call wrigrd(tsw,xmi,xma,ymi,yma,tswmi,tswma,dx,dy
     :   ,0,nx+1,0,ny+1,0,nx+1,0,ny+1,fext(fnmap,'tsw'))
      call wrigrd(t2,xmi,xma,ymi,yma,t2mi,t2ma,dx,dy
     :   ,0,nx+1,0,ny+1,0,nx+1,0,ny+1,fext(fnmap,'t2'))
      call wrigrd(wr,xmi,xma,ymi,yma,wrmi,wrma,dx,dy
     :   ,0,nx+1,0,ny+1,0,nx+1,0,ny+1,fext(fnmap,'wr'))
      call iwrigrd(irr,xmi,xma,ymi,yma,irrmi,irrma,dx,dy
     :   ,0,nx+1,0,ny+1,0,nx+1,0,ny+1,fext(fnmap,'irr'))
      call iwrigrd(isa,xmi,xma,ymi,yma,isami,isama,dx,dy
     :   ,0,nx+1,0,ny+1,0,nx+1,0,ny+1,fext(fnmap,'isa'))
      call iwrigrd(ive,xmi,xma,ymi,yma,ivemi,ivema,dx,dy
     :   ,0,nx+1,0,ny+1,0,nx+1,0,ny+1,fext(fnmap,'ive'))
      call wrigrd(rsm,xmi,xma,ymi,yma,rsmmi,rsmma,dx,dy
     :   ,0,nx+1,0,ny+1,0,nx+1,0,ny+1,fext(fnmap,'rsm'))
      call wrigrd(xdd,xmi,xma,ymi,yma,xddmi,xddma,dx,dy
     :   ,0,nx+1,0,ny+1,0,nx+1,0,ny+1,fext(fnmap,'xdd'))
      end

      subroutine wrigrd(z,xmin,xmax,ymin,ymax,zmin,zmax,dx,dy
     :   ,m0,m1,n0,n1,i0,i1,j0,j1,fname)

      implicit double precision (a-h,o-z)
      character*30 fname
      dimension z(m0:m1,n0:n1)
	real(4) dx,dy

      zmin=z(i0,j0)
       zmax=z(i0,j0)

      do j=j0,j1
         do i=i0,i1
           zmin=min(zmin,z(i,j))
           zmax=max(zmax,z(i,j))
         enddo
      enddo

	xmin=-float(i1)/2.*dx/1000.
	xmax=float(i1)/2.*dx/1000.
	ymin=-float(j1)/2.*dy/1000.
	ymax=float(j1)/2.*dy/1000.

      open(77,file=fname)
      write(77,'(a4)') 'DSAA'
     

      write(77,*) (i1-i0+1),(j1-j0+1)
      write(77,*) xmin,xmax
      write(77,*) ymin,ymax
      write(77,*) zmin,zmax

      do j=j0,j1
         write(77,'(5E15.7)') (z(i,j),i=i0,i1)
      enddo

      close(77)
      return
      end

      subroutine iwrigrd(iz,xmin,xmax,ymin,ymax,izmin,izmax,dx,dy
     :   ,m0,m1,n0,n1,i0,i1,j0,j1,fname)

      implicit double precision (a-h,o-z)
      character*30 fname
      dimension iz(m0:m1,n0:n1)
	real(4) dx,dy

      izmin=iz(i0,j0)
      izmax=iz(i0,j0)

      do j=j0,j1
         do i=i0,i1
           izmin=min(izmin,iz(i,j))
           izmax=max(izmax,iz(i,j))
         enddo
      enddo

	xmin=-float(i1)/2.*dx/1000.
	xmax=float(i1)/2.*dx/1000.
	ymin=-float(j1)/2.*dy/1000.
	ymax=float(j1)/2.*dy/1000.

      open(77,file=fname)
      write(77,'(a4)') 'DSAA'
      write(77,*) (i1-i0+1),(j1-j0+1)
      write(77,*) xmin,xmax
      write(77,*) ymin,ymax
      write(77,*) izmin,izmax

      do j=j0,j1
         write(77,'(8i8)') (iz(i,j),i=i0,i1)
      enddo

      close(77)
      return
      end

      subroutine readgrd(nunit,var,i0,i1,j0,j1,x0,x1,y0,y1,z0,z1)
c
c reads 2d array from a grd file
c
      implicit double precision(a-h,o-z)
      dimension var(i0:i1,j0:j1)
      character*4 dsaa

      read(nunit,'(a4)') dsaa
      read(nunit,*) nx,ny
      read(nunit,*) x0,x1
      read(nunit,*) y0,y1
      read(nunit,*) z0,z1

      do j=j0,j1
         read(nunit,*) (var(i,j),i=i0,i1)
      enddo

      return
      end

      subroutine ireadgrd(nunit,ivar,i0,i1,j0,j1,x0,x1,y0,y1,iz0,iz1)
c
c reads 2d array from a grd file
c
      implicit double precision(a-h,o-z)
      dimension ivar(i0:i1,j0:j1)
      character*4 dsaa

      read(nunit,'(a4)') dsaa
      read(nunit,*) nx,ny
      read(nunit,*) x0,x1
      read(nunit,*) y0,y1
      read(nunit,*) iz0,iz1

      do j=j0,j1
         read(nunit,*) (ivar(i,j),i=i0,i1)
      enddo

      return
      end

      function fext(fname,ext)
      character*40 fname,fext
      character*3 ext
      parameter(n=40)
c
c elimina extensao se existir
c
      fext=fname

      do 10 i=30,1,-1
         if(fext(i:i).eq.'.') then
            fext(i:n)=char(0)
            go to 11
         endif
10    continue
11    continue

      do 20 i=30,1,-1
         if(fext(i:i).ne.' ' .and. fext(i:i).ne.char(0)) then
            fext(i+1:i+4)='.'//ext
            go to 21
         endif
20    continue
      write(*,*) 'erro em fext'
      pause
21    continue
      return
      end
