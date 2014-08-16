!     Dima@ 13.04.07 
!     Summary@ Ground impurity initialization routine
      subroutine impursetup(nx, ny, ns, ix_impur, iy_impur, impur_value,
     :  immediate_impur, constant_impur,
     :  le_impur, qs_10, qs_10s)
      implicit none
      
!     Parameters@
!     (nx,ny,ns) - grid dimensions
!     (ix_impur,iy_imput) - (x,y) coordinate of the impurity signal.
!     immediate_impur, constant_impur - flags indicating the type of impurity to setup.
!     qs_10, qs_10s - variables storing volume impurity distribution.
!       qs_10 - current level
!       qs_10s - leapfrog three levels pattern
!     le_impur - impurity flow
!     radius - the range of the impurity point to disribute exponentially
      
      integer(4):: ix, iy, nx, ny, ns, ix_impur,iy_impur, radius
      real(8)::impur_value
      logical::immediate_impur, constant_impur
      real(8)::le_impur(0:nx,0:ny),
     :  qs_10(0:nx,0:ny,0:ns,1:3), qs_10s(0:nx,0:ny,0:ns)
c
      radius = 1
c
      ix_impur = 8
      iy_impur = 30
      impur_value = 1000000
      immediate_impur = .TRUE.
      constant_impur = .FALSE.
c
      qs_10 = 1
      qs_10s = 1
      le_impur = 0
      if (immediate_impur) then
c      
c       Setup impurity values
c
        do ix=ix_impur-radius,ix_impur+radius
          do iy=iy_impur-radius, iy_impur+radius
            qs_10(ix,iy,20,1) = impur_value
            qs_10s(ix,iy,20) = impur_value
          enddo
        enddo
c
      elseif (constant_impur) then
c      
c       Setup impurity values
c
        do ix=ix_impur-radius,ix_impur+radius
          do iy=iy_impur-radius, iy_impur+radius      
            le_impur(ix,iy) = impur_value
          enddo
        enddo
c                   
      endif
      end
      
      subroutine humid2(wrk)                                                    DM,05.2007
c-----------------------------------------------------------------------
c     time integration of s using the 3d leap-frog scheme
c     Dima@ 04.04.07
c-----------------------------------------------------------------------
c 
      use alloc                                                                 DM,05.2007

      implicit real*8(a-h,o-z)                                                  DM,05.2007
      dimension wrk(0:nx1,0:ny1,0:ns1)                                          DM,05.2007
      
      print *, "humid2"                                                         DM,05.2007
      
      if(tfct) then                                                             DM,05.2007

         do k=1,3                                                               DM,05.2007
         do is=0,ns                                                             DM,05.2007
         do iy=1,ny1                                                            DM,05.2007
         do ix=1,nx1                                                            DM,05.2007
            qs_10(ix,iy,is,k)=qs_10(ix,iy,is,k)*pp(ix,iy,k)                     DM,05.2007
         enddo                                                                  DM,05.2007
         enddo                                                                  DM,05.2007
         enddo                                                                  DM,05.2007
         enddo                                                                  DM,05.2007

         allocate(work1(0:nx1,0:ny1,0:ns1))                                     DM,05.2007
         allocate(work2(0:nx1,0:ny1,0:ns1))                                     DM,05.2007
         allocate(work3(0:nx1,0:ny1,0:ns1))                                     DM,05.2007
         allocate(work4(0:nx1,0:ny1,0:ns1))                                     DM,05.2007
         allocate(work5(0:nx1,0:ny1,0:ns1))                                     DM,05.2007
         allocate(work6(0:nx1,0:ny1,0:ns1))                                     DM,05.2007
         allocate(work7(0:nx1,0:ny1,0:ns1))                                     DM,05.2007
         allocate(work8(0:nx1,0:ny1,0:ns1))                                     DM,05.2007
         call fct3d(qs_10(0,0,0,3),qs_10(0,0,0,2),qs_10(0,0,0,1)                DM,05.2007
     :      ,u(0,0,0,3),u(0,0,0,2)                                              DM,05.2007
     :      ,v(0,0,0,3),v(0,0,0,2)                                              DM,05.2007
     :      ,wsig(0,0,0,3),wsig(0,0,0,2)                                        DM,05.2007
     :      ,nx,ny,ns                                                           DM,05.2007
     :      ,dtxs,dtys,dtxy,dxys,dxdt,dydt,dt,ds0                               DM,05.2007 
     :      ,work1,work2,work3,work4,work5,work6,work7,work8)                   DM,05.2007
         deallocate(work1)                                                      DM,05.2007 
         deallocate(work2)                                                      DM,05.2007
         deallocate(work3)                                                      DM,05.2007
         deallocate(work4)                                                      DM,05.2007
         deallocate(work5)                                                      DM,05.2007
         deallocate(work6)                                                      DM,05.2007
         deallocate(work7)                                                      DM,05.2007
         deallocate(work8)                                                      DM,05.2007
         do k=1,3                                                               DM,05.2007
         do is=0,ns                                                             DM,05.2007
         do iy=1,ny1                                                            DM,05.2007 
         do ix=1,nx1                                                            DM,05.2007
           qs_10(ix,iy,is,k)=qs_10(ix,iy,is,k)/pp(ix,iy,k)                      DM,05.2007 
         enddo                                                                  DM,05.2007
         enddo                                                                  DM,05.2007
         enddo                                                                  DM,05.2007 
         enddo                                                                  DM,05.2007
         do is=1,ns-1                                                           DM,05.2007
         do iy=2,ny                                                             DM,05.2007  
         do ix=2,nx                                                             DM,05.2007
           qs_10(ix,iy,is,3)=qs_10(ix,iy,is,3)+dtl*difunq_10(ix,iy,is)          DM,05.2007  
     :        /pp(ix,iy,3)                                                      DM,05.2007
         enddo                                                                  DM,05.2007
         enddo                                                                  DM,05.2007
         enddo                                                                  DM,05.2007
      else                                                                      DM,05.2007 
      
        !print *,"Start check 3"

        !do i = 0, nx1
        !    do j = 0, ny1
        !    do k = 0, ns+1
        !        
        !        do lev = 1,3
        !        if (qs_10(i,j,k,lev).NE.1d-5) then
        !            print *,qs_10(i,j,k,lev)
        !        endif
        !        enddo
        !        
        !    enddo
        !    enddo
        !    read *
        !enddo       
      
        do is=1,ns-1                                                            DM,05.2007
        do iy=2,ny                                                              DM,05.2007
        do ix=2,nx                                                              DM,05.2007
          tavx=(u(ix,iy,is,2)*pp10(ix,iy,2)*(qs_10(ix+1,iy,is,2)                DM,05.2007
     :      +qs_10(ix,iy,is,2))                                                 DM,05.2007
     :      -u(ix-1,iy,is,2)*pp10(ix-1,iy,2)*(qs_10(ix,iy,is,2)                 DM,05.2007
     :      +qs_10(ix-1,iy,is,2)))/dx2                                          DM,05.2007
          tavy=(v(ix,iy,is,2)*pp01(ix,iy,2)*(qs_10(ix,iy+1,is,2)                DM,05.2007
     :      +qs_10(ix,iy,is,2))                                                 DM,05.2007
     :      -v(ix,iy-1,is,2)*pp01(ix,iy-1,2)*(qs_10(ix,iy,is,2)                 DM,05.2007
     :      +qs_10(ix,iy-1,is,2)))/dy2                                          DM,05.2007
          tavs=(wsig(ix,iy,is+1,2)*                                             DM,05.2007
     :      (qs_10(ix,iy,is+1,2)+qs_10(ix,iy,is,2))                             DM,05.2007  
     :      -wsig(ix,iy,is,2)*(qs_10(ix,iy,is,2)+qs_10(ix,iy,is-1,2)))          DM,05.2007
     :      /ds12(is)                                                           DM,05.2007
          qs_10(ix,iy,is,3)=(qs_10(ix,iy,is,1)*pp(ix,iy,1)                      DM,05.2007
     :      -dtl*(tavx+tavy+pp(ix,iy,2)                                         DM,05.2007
     :      *(tavs-difunq_10(ix,iy,is)))) ! ))) is replacing +cond(ix,iy,is)-evap(ix,iy,is))))
     :      /pp(ix,iy,3)                                                        DM,05.2007  
        enddo                                                                   DM,05.2007
        enddo                                                                   DM,05.2007
        enddo                                                                   DM,05.2007
        
        !print *,"Check diffunq"
        !read *
        !
        !do i = 0, nx1
        !    do j = 0, ny1
        !    do k = 0, ns+1
        !        print *,difunq_10(i,j,k)
        !    enddo
        !    enddo
        !    read *
        !enddo        
        
        !print *,"Start check 4"

        !do i = 0, nx1
        !    do j = 0, ny1
        !    do k = 0, ns+1
        !        
        !        do lev = 1,3
        !        if (qs_10(i,j,k,lev).NE.1d-5) then
        !            print *,qs_10(i,j,k,lev)
        !        endif
        !        enddo
        !        
        !    enddo
        !    enddo
        !    read *,skip
        !    if (skip.EQ.0) goto 999
        !enddo         
        
        endif                                                                   DM,05.2007
c
c
c surface water vapor flux (bulk parametrization):
c
c
c&solo
c
      if(ifsoil.ne.0) then                                                      DM,05.2007
        do iy=2,ny                                                              DM,05.2007
        do ix=2,nx                                                              DM,05.2007
          qs_10(ix,iy,ns-1,3)=qs_10(ix,iy,ns-1,3)+dtl*le_impur(ix,iy)           DM,05.2007
     :      /deltaz(ix,iy)/hlat/rhoa                                            DM,05.2007
        enddo                                                                   DM,05.2007

        enddo                                                                   DM,05.2007
      elseif(cdcoef.gt.0.) then                                                 DM,05.2007 
         ctbulk=cdcoef*dtl                                                      DM,05.2007
         do iy=2,ny                                                             DM,05.2007
           do ix=2,nx                                                           DM,05.2007
             p=sigma0(ns-1)*pp(ix,iy,2)+ptop                                    DM,05.2007
             qs_10(ix,iy,ns-1,3)=qs_10(ix,iy,ns-1,3)+ctbulk                     DM,05.2007
     :         *uvsuf(ix,iy)*(qvsoil(ix,iy)-qvs(ix,iy,ns-1)                     DM,05.2007
     :         -qs_10(ix,iy,ns-1,2))/deltaz(ix,iy)                              DM,05.2007 
           enddo                                                                DM,05.2007
         enddo                                                                  DM,05.2007
      endif                                                                     DM,05.2007
         
c
c horizontal boundary conditions:
c
      call ptbc(nx,ny,ns,qs_10(0,0,0,3),qsbx,qsby,qscc,iobptx,iobpty            DM,05.2007
     :   ,ptbout,prt,nchn1)                                                     DM,05.2007
c
      if(raylei) then                                                           DM,05.2007
         do is=1,idrmax                                                         DM,05.2007 
         do iy=1,ny1                                                            DM,05.2007
         do ix=1,nx1                                                            DM,05.2007
           qs_10(ix,iy,is,3)=qs_10(ix,iy,is,3)-dtl/taudra(is)                   DM,05.2007
     :        *(qs_10(ix,iy,is,2)-qs_10s(ix,iy,is))                             DM,05.2007
         enddo                                                                  DM,05.2007
         enddo                                                                  DM,05.2007
         enddo                                                                  DM,05.2007
      endif                                                                     DM,05.2007
c
      do iy=1,ny1                                                               DM,05.2007
      do ix=1,nx1                                                               DM,05.2007
        qs_10(ix,iy,ns,3)=qs_10(ix,iy,ns-1,3)                                   DM,05.2007
        qs_10(ix,iy,0,3)=qs_10(ix,iy,1,3)                                       DM,05.2007
      enddo                                                                     DM,05.2007
      enddo                                                                     DM,05.2007
c
      if(dohsmo .and. mod(nstep,numsmo).eq.0) then                              DM,05.2007
         call hsmoot(qs_10(0,0,0,3),hk1,nx1,ny1,ns,1,nx1,1,ny1,1,ns-1           DM,05.2007
     :      ,hdamp)                                                             DM,05.2007 
      endif                                                                     DM,05.2007
      if(dovsmt .and. mod(nstep,numsmo).eq.0) then                              DM,05.2007
         !print *,"humid2"                                                      DM,05.2007
         call vsmoot(qs_10(0,0,0,3),wrk,nx1,ny1,ns,1,nx1,1,ny1,1,ns-1           DM,05.2007
     :      ,vdampt)                                                            DM,05.2007 
      endif                                                                     DM,05.2007
c
      if(.not.tfct) then                                                        DM,05.2007
         call aselin(qs_10(0,0,0,3),qs_10(0,0,0,2),qs_10(0,0,0,1)               DM,05.2007
     :   ,0,nx1,0,ny1,0,ns,filt)                                                DM,05.2007 
      endif                                                                     DM,05.2007
     
      call radbch(qscc,qsbx,qsby,qs_10(0,0,0,3),qs_10(0,0,0,2),                 DM,05.2007
     : qs_10(0,0,0,1),nx1,ny1,ns1,1,ns-1,1,nx1,1,ny1)                           DM,05.2007

      return                                                                    DM,05.2007
      end                                                                       DM,05.2007
