      Subroutine diffu(fngrd
     :  ,def11,h1,q1
     :  ,def12,h2,q2
     :  ,def13,h3,q3
     :  ,def22
     :  ,def23
     :  ,def33
     :  ,dif000
     :  ,dif010
     :  ,dif001,dif100)

      use alloc

      implicit real*8 (a-h,o-z)
      real(kind(0.d0)),dimension(0:nx1,0:ny1,0:ns1)::
     :  def11,def12,def13,def22,def23,def33
     :  ,h1,h2,h3,q1,q2,q3,dif000,dif010,dif001,dif100
c
c 9 work spaces
c each line one wk
c-----------------------------------------------------------------------
c  compute diffusion terms which are deformation and richardson number
c  dependent after lilly. if rikey=0.0,the later dependency is dropped.
c  enhanced upper-level dampping is determined by dragt,dragm,etc.
c
c  difunu,difunv and difunw are multiplied for pp, to save in ellipt.
c
c rayleigh damping passed to prognostic equations.
c
c-----------------------------------------------------------------------
c

      character*80 fngrd

      parameter (third4=4./3.,third2=2./3.)

c
c calculation of the deformation tensor:
c

c      allocate(def11(0:nx1,0:ny1,0:ns1))
c      allocate(def12(0:nx1,0:ny1,0:ns1))
c      allocate(def13(0:nx1,0:ny1,0:ns1))
c      allocate(def22(0:nx1,0:ny1,0:ns1))
c      allocate(def23(0:nx1,0:ny1,0:ns1))
c      allocate(def33(0:nx1,0:ny1,0:ns1))
c      allocate(dif000(0:nx1,0:ny1,0:ns1))
        if(verbose.ge.2) write(*,*) 'enter diffu'
        
        do ix=1,nx1
	w(ix,ny1,ns+1,2)=w(ix,ny1-1,ns+1,2)
	enddo

	do iy=1,ny1
	w(nx1,iy,ns+1,2)=w(nx1-1,iy,ns+1,2)
	enddo
      
	
      do is=1,ns-1
      do iy=1,ny1
      do ix=1,nx1
        def33(ix,iy,is)=-(s(ix,iy,is+1)+s(ix,iy,is))
     :    *(w(ix,iy,is+1,2)-w(ix,iy,is,2))/ds12(is)
	enddo
      enddo
      enddo
      

      do is=1,ns-1
      do iy=1,ny1
      do ix=1,nx1
        zdudx=(u(ix,iy,is,2)-u(ix-1,iy,is,2))/dx
     :    -s0ds4a(is)*ppdx(ix,iy,2)*(u(ix-1,iy,is+1,2)
     :    +u(ix,iy,is+1,2)-u(ix-1,iy,is-1,2)-u(ix,iy,is-1,2))
        zdvdy=(v(ix,iy,is,2)-v(ix,iy-1,is,2))/dy
     :    -s0ds4a(is)*ppdy(ix,iy,2)*(v(ix,iy-1,is+1,2)
     :    +v(ix,iy,is+1,2)-v(ix,iy-1,is-1,2)-v(ix,iy,is-1,2))
        zdwdz=-(s(ix,iy,is+1)+s(ix,iy,is))
     :    *(w(ix,iy,is+1,2)-w(ix,iy,is,2))/ds12(is)
        def11(ix,iy,is)=third4*zdudx-third2*(zdvdy+zdwdz)
        def22(ix,iy,is)=third4*zdvdy-third2*(zdudx+zdwdz)
        def33(ix,iy,is)=third4*zdwdz-third2*(zdudx+zdvdy)
      enddo
      enddo
      enddo
c
c boundary conditions d/ds=0
c
      do is=0,ns,ns
      do iy=1,ny1
      do ix=1,nx1
        zdudx=(u(ix,iy,is,2)-u(ix-1,iy,is,2))/dx
        zdvdy=(v(ix,iy,is,2)-v(ix,iy-1,is,2))/dy
        zdwdz=-(s(ix,iy,is+1)+s(ix,iy,is))
     :    *(w(ix,iy,is+1,2)-w(ix,iy,is,2))/ds12(is)
        def11(ix,iy,is)=third4*zdudx-third2*(zdvdy+zdwdz)
        def22(ix,iy,is)=third4*zdvdy-third2*(zdudx+zdwdz)
        def33(ix,iy,is)=third4*zdwdz-third2*(zdudx+zdvdy)
      enddo
      enddo
      enddo
c
      do is=1,ns-1
      do iy=1,ny
      do ix=1,nx
        def12(ix,iy,is)=
     :    (u(ix,iy+1,is,2)-u(ix,iy,is,2))/dy
     :    -s0ds4a(is)*0.5*(ppdy01(ix,iy,2)+ppdy01(ix+1,iy,2))
     :    *(u(ix,iy,is+1,2)
     :    +u(ix,iy+1,is+1,2)-u(ix,iy,is-1,2)-u(ix,iy+1,is-1,2))
     :    +(v(ix+1,iy,is,2)-v(ix,iy,is,2))/dx
     :    -s0ds4a(is)*0.5*(ppdx10(ix,iy,2)+ppdx10(ix,iy+1,2))
     :    *(v(ix+1,iy,is+1,2)
     :    +v(ix,iy,is+1,2)-v(ix+1,iy,is-1,2)-v(ix,iy,is-1,2))
      enddo
      enddo
      enddo
c
      do is=0,ns,ns
      do iy=1,ny
      do ix=1,nx
        def12(ix,iy,is)=
     :    (u(ix,iy+1,is,2)-u(ix,iy,is,2))/dy
     :    +(v(ix+1,iy,is,2)-v(ix,iy,is,2))/dx
      enddo
      enddo
      enddo
c
      do is=1,ns
      do iy=1,ny1
      do ix=1,nx
        def13(ix,iy,is)=-(s(ix+1,iy,is)+s(ix,iy,is))
     :    *(u(ix,iy,is,2)-u(ix,iy,is-1,2))/ds02(is)
     :    +(w(ix+1,iy,is,2)-w(ix,iy,is,2))/dx-s1ds4a(is)*ppdx10(ix,iy,2)
     :    *(w(ix+1,iy,is+1,2)+w(ix,iy,is+1,2)-w(ix+1,iy,is-1,2)
     :    -w(ix,iy,is-1,2))
      enddo
      enddo
      enddo
c
      do is=1,ns
      do iy=1,ny
      do ix=1,nx1
        def23(ix,iy,is)=-(s(ix,iy,is)+s(ix,iy+1,is))
     :    *(v(ix,iy,is,2)-v(ix,iy,is-1,2))/ds02(is)
     :    +(w(ix,iy+1,is,2)-w(ix,iy,is,2))/dy-s1ds4a(is)*ppdy01(ix,iy,2)
     :    *(w(ix,iy+1,is+1,2)+w(ix,iy,is+1,2)-w(ix,iy+1,is-1,2)
     :    -w(ix,iy,is-1,2))
      enddo
      enddo
      enddo
c
c debug
c
c      call extrah(nx1,ny1,def11,1,nx,1,ny,0,ns)
c      call extrah(nx1,ny1,def22,1,nx,1,ny,0,ns)
c      call extrah(nx1,ny1,def33,1,nx,1,ny,0,ns)
c      call extrah(nx1,ny1,def12,1,nx,1,ny,0,ns)
c      call extrah(nx1,ny1,def13,1,nx,1,ny,0,ns)
c      call extrah(nx1,ny1,def23,1,nx,1,ny,0,ns)
c
      if(defout.and.prt) then
        tk=0.
        call wri3ar(def11,1,nx1,1,ny1,1,ns-1,'def11   ',tk,iomod,nx,ny)
        call wri3ar(def12,1,nx,1,ny,1,ns-1,'def12   ',tk,iomod,nx,ny)
        call wri3ar(def13,1,nx,1,ny1,1,ns-1,'def13   ',tk,iomod,nx,ny)
        call wri3ar(def22,1,nx1,1,ny1,1,ns-1,'def22   ',tk,iomod,nx,ny)
        call wri3ar(def23,1,nx1,1,ny,1,ns-1,'def23   ',tk,iomod,nx,ny)
        call wri3ar(def33,1,nx1,1,ny1,1,ns-1,'def33   ',tk,iomod,nx,ny)
      endif
c
c calculation of diffusion coefficients:
c
      if(iodif.eq.1) then
c-
c calculate linear diffusion coefficients.
c-
        do is=0,ns
        do iy=1,ny1
        do ix=1,nx1
          dif000(ix,iy,is)=difcof
        enddo
        enddo
        enddo
      else
c-
c calculate deformation/richardson number dependent diffusion coefficient
c-
        if(prt.and.driout) allocate(ri(0:nx1,0:ny1,0:ns1))

        do is=1,ns-1
        do iy=2,ny
        do ix=2,nx
c
          def=dsqrt(0.5*(def11(ix,iy,is)*def11(ix,iy,is)+def22(ix,iy,is)
     :      *def22(ix,iy,is)+def33(ix,iy,is)*def33(ix,iy,is))
     :      +0.25*(def12(ix,iy,is)*def12(ix,iy,is)
     :      +def12(ix-1,iy,is)*def12(ix-1,iy,is)
     :      +def12(ix,iy-1,is)*def12(ix,iy-1,is)
     :      +def12(ix-1,iy-1,is)*def12(ix-1,iy-1,is)
     :      +def13(ix,iy,is)*def13(ix,iy,is)
     :      +def13(ix-1,iy,is)*def13(ix-1,iy,is)
     :      +def13(ix,iy,is+1)*def13(ix,iy,is+1)
     :      +def13(ix-1,iy,is+1)*def13(ix-1,iy,is+1)
     :      +def23(ix,iy,is)*def23(ix,iy,is)
     :      +def23(ix,iy-1,is)*def23(ix,iy-1,is)
     :      +def23(ix,iy,is+1)*def23(ix,iy,is+1)
     :      +def23(ix,iy-1,is+1)*def23(ix,iy-1,is+1)))
c

          if(qif.ne.0. .and. ifqc.ne.0) then
            p=sigma0(is)*pp(ix,iy,2)+ptop
            p2=sigma0(is+1)*pp(ix,iy,2)+ptop
            p1=sigma0(is-1)*pp(ix,iy,2)+ptop
            t=(pts(ix,iy,is)+pt(ix,iy,is,2))*(p/p00)**akapa
            tp2=(pts(ix,iy,is+1)+pt(ix,iy,is+1,2))*(p2/p00)**akapa
            tp1=(pts(ix,iy,is-1)+pt(ix,iy,is-1,2))*(p1/p00)**akapa
            qsatur=qsat(t,p)
            qsat2=qsat(tp2,p2)
            qsat1=qsat(tp1,p1)
            xnsq=-g*(s(ix,iy,is+1)+s(ix,iy,is))*
     :      ((1+hlat*qsatur/(rv*t))/(1+0.622*hlat*hlatcp*qsatur/
     :      (rv*t*t))
     :      *((pt(ix,iy,is+1,2)
     :      +pts(ix,iy,is+1)-pt(ix,iy,is-1,2)-pts(ix,iy,is-1))
     :      /(pt(ix,iy,is,2)+pts(ix,iy,is))+hlatcp/t*(qsat2-qsat1))-
     :      (qv(ix,iy,is+1,2)+qc(ix,iy,is+1,2)-qv(ix,iy,is-1,2)-
     :      qc(ix,iy,is-1,2)))
     :      /ds04a(is)
          else
            xnsq=-g*(s(ix,iy,is+1)+s(ix,iy,is))*(pt(ix,iy,is+1,2)
     :        +pts(ix,iy,is+1)-pt(ix,iy,is-1,2)-pts(ix,iy,is-1))
     :        /(pt(ix,iy,is,2)+pts(ix,iy,is))
     :        /ds04a(is)
          endif
          rich=xnsq/(def*def+zero0)
          if(prt.and.driout) ri(ix,iy,is)=rich
c
          z_dif=(phis(ix,iy,ns-1)+phi(ix,iy,ns-1))/g-hsuf(ix,iy)
          delka3=(1./delka2(is)+1/((0.4*z_dif)**2))**(-1)
          difk=delka2(is)*def*dsqrt(dim(1.0d0,rich*rikey))
!          difk=delka3*def*dsqrt(max(0.d0,1.0d0-rich*rikey))
! for horizontal eddy-diffusivity horizontal grid spacing is used as the length scale:
          difk_h=(dx*dkapa)**2*def*dsqrt(max(0.d0,1.0d0-rich*rikey))
c
          dif000(ix,iy,is)=difk+difl
          dif000_h(ix,iy,is)=difk_h+difl
        enddo
        enddo
        enddo
c
        call extra3(nx1,ny1,ns1,dif000,1,nx1,1,ny1,0,ns)
        call extra3(nx1,ny1,ns1,dif000_h,1,nx1,1,ny1,0,ns)
c
        if(prt.and.driout) then
          if(grdout) then
            call wgrids(ri,'ri',2,nx,2,ny,2,ns-1,iomod
     :        ,0,0,0,fngrd)
            call wgrids(dif000,'kk',2,nx,2,ny,2,ns-1,iomod
     :        ,0,0,0,fngrd)
!            call wgrids(dif000_h,'kh',2,nx,2,ny,2,ns-1,iomod
!     :        ,0,0,0,fngrd)
          else
            tk=0.
            call wri3ar(ri,2,nx,2,ny,2,ns-1,'ri      ',tk,iomod,nx,ny)
            call wri3ar(dif000,2,nx,2,ny,2,ns-1,'dif000  ',tk,iomod
     :        ,nx,ny)
          endif
          deallocate(ri)
        endif
c

      endif
c
c top boundary enhanced diffusion:
c
      if(olddra) then
        ixm=nx1/2
        iym=ny1/2
        ztop=phis(ixm,iym,1)/g
        do is=0,ns
        do iy=1,ny1
        do ix=1,nx1
          z=(phi(ix,iy,is)+phis(ix,iy,is))/g
          drag=(dragm+(dragt-dragm)*cos(pi05*abs((z-ztop)/(zm-ztop)))
     :      **2)*dim(z,zm)/(z-zm+zero0)
          dif000(ix,iy,is)=dif000(ix,iy,is)+drag
          dif000_h(ix,iy,is)=dif000_h(ix,iy,is)+drag    ! horizontal diffusivity added
        enddo
        enddo
        enddo
      elseif(.not.raylei) then
        do is=0,idrmax
        do iy=1,ny1
        do ix=1,nx1
          dif000(ix,iy,is)=dif000(ix,iy,is)+endrag(is)
          dif000_h(ix,iy,is)=dif000_h(ix,iy,is)+endrag(is)   ! horizontal diffusivity added
        enddo
        enddo
        enddo
      endif
c
c
c impose linear stability for diffusion coefficient
c
      do is=0,ns
      do iy=1,ny1
      do ix=1,nx1
        dif000(ix,iy,is)=min(deltaz(ix,iy)**2/(8.*dt),dif000(ix,iy,is))
      enddo
      enddo
      enddo

! impose linear stability for horizontal diffusion coefficient
!
      do is=0,ns
      do iy=1,ny1
      do ix=1,nx1
        dif000_h(ix,iy,is)=min(dx**2/(8.*dt),dif000_h(ix,iy,is))
      enddo
      enddo
      enddo

      if(prt.and.ddiout) then
         tk=0.
         call wri3ar(dif000,1,nx1,1,ny1,0,ns,'difcof  ',tk,iomod,nx,ny)
      endif
c
c the upper layer drag is non-linear. so maybe it should be calculated
c separatelly for the other grid.
c
c diffusion coefficients on other grid locations:
c

c      allocate(dif010(0:nx1,0:ny1,0:ns1))
c      allocate(dif100(0:nx1,0:ny1,0:ns1))

      do is=0,ns
      do iy=1,ny
      do ix=1,nx
        dif010(ix,iy,is)=0.5*(dif000(ix,iy,is)+dif000(ix,iy+1,is))
        dif100(ix,iy,is)=0.5*(dif000(ix,iy,is)+dif000(ix+1,iy,is))
        dif010_h(ix,iy,is)=0.5*(dif000_h(ix,iy,is)+dif000_h(ix,iy+1,is))  ! horizontal diffusivity
        dif100_h(ix,iy,is)=0.5*(dif000_h(ix,iy,is)+dif000_h(ix+1,iy,is))  ! horizontal diffusivity
      enddo
      enddo
      enddo
c
      call extrah(nx1,ny1,dif010,0,nx1,0,ny1,0,ns)
      call extrah(nx1,ny1,dif100,0,nx1,0,ny1,0,ns)
      call extrah(nx1,ny1,dif010_h,0,nx1,0,ny1,0,ns)    ! horizontal diffusivity
      call extrah(nx1,ny1,dif100_h,0,nx1,0,ny1,0,ns)    ! horizontal diffusivity
c
c diffuse only the perturbated field:
c
      if(uvdif.ne.1.0) then    ! different horizontal diffusivity is not included for this option
        do is=1,ns-1
        do iy=1,ny1
        do ix=1,nx1
          zdusdx=(us(ix,iy,is)-us(ix-1,iy,is))/dx
     :      -s0ds4a(is)*ppdx(ix,iy,2)*(us(ix-1,iy,is+1)
     :      +us(ix,iy,is+1)-us(ix-1,iy,is-1)-us(ix,iy,is-1))
          zdvsdy=(vs(ix,iy,is)-vs(ix,iy-1,is))/dy
     :      -s0ds4a(is)*ppdy(ix,iy,2)*(vs(ix,iy-1,is+1)
     :      +vs(ix,iy,is+1)-vs(ix,iy-1,is-1)-vs(ix,iy,is-1))
c
          def11s=third4*zdusdx-third2*zdvsdy
          def11(ix,iy,is)=(def11(ix,iy,is)-(1.0-uvdif)*def11s)
     :      *dif000(ix,iy,is)
          def22s=third4*zdvsdy-third2*zdusdx
          def22(ix,iy,is)=(def22(ix,iy,is)-(1.0-uvdif)*def22s)
     :      *dif000(ix,iy,is)
c
          def33s=-third2*(zdusdx+zdvsdy)
          def33(ix,iy,is)=(def33(ix,iy,is)-(1.0-uvdif)*def33s)
     :      *dif000(ix,iy,is)
        enddo
        enddo
        enddo
c
        do is=0,ns,ns
        do iy=1,ny1
        do ix=1,nx1
          zdusdx=(us(ix,iy,is)-us(ix-1,iy,is))/dx
          zdvsdy=(vs(ix,iy,is)-vs(ix,iy-1,is))/dy
c
          def11s=third4*zdusdx-third2*zdvsdy
          def11(ix,iy,is)=(def11(ix,iy,is)-(1.0-uvdif)*def11s)
     :      *dif000(ix,iy,is)
c
          def22s=third4*zdvsdy-third2*zdusdx
          def22(ix,iy,is)=(def22(ix,iy,is)-(1.0-uvdif)*def22s)
     :      *dif000(ix,iy,is)
c
          def33s=-third2*(zdusdx+zdvsdy)
          def33(ix,iy,is)=(def33(ix,iy,is)-(1.0-uvdif)*def33s)
     :      *dif000(ix,iy,is)
c
        enddo
        enddo
        enddo
c
        do is=1,ns-1
        do iy=1,ny
        do ix=1,nx
          def12s=(us(ix,iy+1,is)-us(ix,iy,is))/dy
     :      -s0ds4a(is)*0.5*(ppdy01(ix,iy,2)+ppdy01(ix+1,iy,2))
     :      *(us(ix,iy+1,is+1)
     :      +us(ix,iy,is+1)-us(ix,iy+1,is-1)-us(ix,iy,is-1))
     :      +(vs(ix+1,iy,is)-vs(ix,iy,is))/dx
     :      -s0ds4a(is)*0.5*(ppdx10(ix,iy,2)+ppdx10(ix,iy+1,2))
     :      *(vs(ix+1,iy,is+1)
     :      +vs(ix,iy,is+1)-vs(ix+1,iy,is-1)-vs(ix,iy,is-1))
          def12(ix,iy,is)=(def12(ix,iy,is)-(1.0-uvdif)*def12s)
     :      *0.5*(dif100(ix,iy,is)+dif100(ix,iy+1,is))
        enddo
        enddo
        enddo
c
        do is=0,ns,ns
        do iy=1,ny
        do ix=1,nx
          def12s=(us(ix,iy+1,is)-us(ix,iy,is))/dy
     :      +(vs(ix+1,iy,is)-vs(ix,iy,is))/dx
          def12(ix,iy,is)=(def12(ix,iy,is)-(1.0-uvdif)*def12s)
     :      *0.5*(dif100(ix,iy,is)+dif100(ix,iy+1,is))
        enddo
        enddo
        enddo
c
        do is=1,ns
        do iy=1,ny1
        do ix=1,nx
          def13s=-(s(ix+1,iy,is)+s(ix,iy,is))
     :      *(us(ix,iy,is)-us(ix,iy,is-1))/ds02(is)
          def13(ix,iy,is)=(def13(ix,iy,is)-(1.0-uvdif)*def13s)
     :      *0.5*(dif100(ix,iy,is)+dif100(ix,iy,is-1))
        enddo
        enddo
        enddo
c
        do is=1,ns
        do iy=1,ny
        do ix=1,nx1
          def23s=-(s(ix,iy+1,is)+s(ix,iy,is))
     :      *(vs(ix,iy,is)-vs(ix,iy,is-1))/ds02(is)
          def23(ix,iy,is)=(def23(ix,iy,is)-(1.0-uvdif)*def23s)
     :      *0.5*(dif010(ix,iy,is)+dif010(ix,iy,is-1))
        enddo
        enddo
        enddo
c
      else                   !here different horizontal diffusivity is already included
         do is=0,ns
         do iy=1,ny1
         do ix=1,nx1
          ! def11(ix,iy,is)=def11(ix,iy,is)*dif000(ix,iy,is)
            def11(ix,iy,is)=def11(ix,iy,is)*dif000_h(ix,iy,is)
          ! def22(ix,iy,is)=def22(ix,iy,is)*dif000(ix,iy,is)
            def22(ix,iy,is)=def22(ix,iy,is)*dif000_h(ix,iy,is)
           def33(ix,iy,is)=def33(ix,iy,is)*dif000(ix,iy,is)
         enddo
         enddo
         enddo
c
         do is=0,ns
         do iy=1,ny
         do ix=1,nx
!           def12(ix,iy,is)=def12(ix,iy,is)*0.5*(dif100(ix,iy,is)
!     :       +dif100(ix,iy+1,is))
            def12(ix,iy,is)=def12(ix,iy,is)*0.5*(dif100_h(ix,iy,is)
     :       +dif100_h(ix,iy+1,is))
         enddo
         enddo
         enddo
c
         do is=1,ns
         do iy=1,ny1
         do ix=1,nx
           def13(ix,iy,is)=def13(ix,iy,is)*0.5*(dif100(ix,iy,is)
     :       +dif100(ix,iy,is-1))
         enddo
         enddo
         enddo
c
         do is=1,ns
         do iy=1,ny
         do ix=1,nx1
           def23(ix,iy,is)=def23(ix,iy,is)*0.5*(dif010(ix,iy,is)
     :       +dif010(ix,iy,is-1))
         enddo
         enddo
         enddo
c
      endif
      
      do iy=1,ny1
	do ix=1,nx
	def13(ix,iy,ns)=0.5*(cdm(ix,iy)+cdm(ix+1,iy))*u(ix,iy,ns-1,2)
	enddo
	enddo
	do iy=1,ny
         do ix=1,nx1
      def23(ix,iy,ns)=0.5*(cdm(ix,iy)+cdm(ix,iy+1))*v(ix,iy,ns-1,2)
	enddo
	enddo
c--
c difunu
c--
      do is=1,ns-1
      do iy=2,ny
      do ix=1,nx
        difunu(ix,iy,is)=((def11(ix+1,iy,is)-def11(ix,iy,is))/dx
     :    -s0ds4a(is)*ppdx10(ix,iy,2)*(def11(ix+1,iy,is+1)
     :    +def11(ix,iy,is+1)-def11(ix+1,iy,is-1)-def11(ix,iy,is-1))
     :    +(def12(ix,iy,is)-def12(ix,iy-1,is))/dy
     :    -s0ds4a(is)*0.5*(ppdy(ix,iy,2)+ppdy(ix+1,iy,2))
     :    *(def12(ix,iy-1,is+1)
     :    +def12(ix,iy,is+1)-def12(ix,iy-1,is-1)-def12(ix,iy,is-1))
     :    -(s(ix+1,iy,is+1)+s(ix+1,iy,is)+s(ix,iy,is+1)+s(ix,iy,is))
     :    *(def13(ix,iy,is+1)-def13(ix,iy,is))/ds14(is))*pp10(ix,iy,2)
      enddo
      enddo
      enddo
c
      do iy=2,ny
      do ix=1,nx
        difunu(ix,iy,0)=difunu(ix,iy,1)
        difunu(ix,iy,ns)=difunu(ix,iy,ns-1)
      enddo
      enddo
c
c--
c difunv
c--
      do is=1,ns-1
      do iy=1,ny
      do ix=2,nx
        difunv(ix,iy,is)=((def22(ix,iy+1,is)-def22(ix,iy,is))/dy
     :    -s0ds4a(is)*ppdy01(ix,iy,2)*(def22(ix,iy+1,is+1)
     :    +def22(ix,iy,is+1)-def22(ix,iy+1,is-1)-def22(ix,iy,is-1))
     :    +(def12(ix,iy,is)-def12(ix-1,iy,is))/dx
     :    -s0ds4a(is)*0.5*(ppdx(ix,iy,2)+ppdx(ix,iy+1,2))
     :    *(def12(ix-1,iy,is+1)
     :    +def12(ix,iy,is+1)-def12(ix-1,iy,is-1)-def12(ix,iy,is-1))
     :    -(s(ix,iy+1,is+1)+s(ix,iy+1,is)+s(ix,iy,is+1)+s(ix,iy,is))
     :    *(def23(ix,iy,is+1)-def23(ix,iy,is))/ds14(is))*pp01(ix,iy,2)
      enddo
      enddo
      enddo
c
      do iy=1,ny
      do ix=2,nx
        difunv(ix,iy,0)=difunv(ix,iy,1)
        difunv(ix,iy,ns)=difunv(ix,iy,ns-1)
      enddo
      enddo
c--
c  difunw
c--
      do is=2,ns-1
      do iy=2,ny
      do ix=2,nx
        difunw(ix,iy,is)=((def13(ix,iy,is)-def13(ix-1,iy,is))/dx
     :    -s1ds4a(is)*ppdx(ix,iy,2)*(def13(ix,iy,is+1)
     :    +def13(ix-1,iy,is+1)-def13(ix,iy,is-1)-def13(ix-1,iy,is-1))
     :    +(def23(ix,iy,is)-def23(ix,iy-1,is))/dy
     :    -s1ds4a(is)*ppdy(ix,iy,2)*(def23(ix,iy,is+1)
     :    +def23(ix,iy-1,is+1)-def23(ix,iy,is-1)-def23(ix,iy-1,is-1))
     :    -s(ix,iy,is)*(def33(ix,iy,is)-def33(ix,iy,is-1))/ds0(is))
     :    *pp(ix,iy,2)
      enddo
      enddo
      enddo
c
      do iy=2,ny
      do ix=2,nx
        difunw(ix,iy,ns)=difunw(ix,iy,ns-1)
        difunw(ix,iy,1)=difunw(ix,iy,2)
      enddo
      enddo
c-
c difunt
c-
c      deallocate(def11)
c      deallocate(def12)
c      deallocate(def13)
c      deallocate(def22)
c      deallocate(def23)
c      deallocate(def33)
c      allocate(h1(0:nx1,0:ny1,0:ns1))
c      allocate(h2(0:nx1,0:ny1,0:ns1))
c      allocate(h3(0:nx1,0:ny1,0:ns1))

      do is=1,ns-1
      do iy=1,ny
      do ix=1,nx
        h1(ix,iy,is)=(pt(ix+1,iy,is,2)-pt(ix,iy,is,2)
     :    +tsdif*(pts(ix+1,iy,is)-pts(ix,iy,is)))/dx
     :    -s0ds4a(is)*ppdx10(ix,iy,2)*(pt(ix+1,iy,is+1,2)
     :    +pt(ix,iy,is+1,2)-pt(ix+1,iy,is-1,2)-pt(ix,iy,is-1,2)
     :    +tsdif*(pts(ix+1,iy,is+1)+pts(ix,iy,is+1)
     :    -pts(ix+1,iy,is-1)-pts(ix,iy,is-1)))
        !h1(ix,iy,is)=h1(ix,iy,is)*dif100(ix,iy,is)*tdif
         h1(ix,iy,is)=h1(ix,iy,is)*dif100_h(ix,iy,is)*tdif     !horizontal diffusivity is used
      enddo
      enddo
      enddo
c
      do iy=1,ny
      do ix=1,nx
c       h1(ix,iy,0)=-h1(ix,iy,1)
c       h1(ix,iy,ns)=-h1(ix,iy,ns-1)
        h1(ix,iy,ns)=(pt(ix+1,iy,ns,2)-pt(ix,iy,ns,2)
     :    +tsdif*(pts(ix+1,iy,ns)-pts(ix,iy,ns)))/dx
        h1(ix,iy,ns)=h1(ix,iy,ns)*dif100(ix,iy,ns)*tdif
        h1(ix,iy,0)=(pt(ix+1,iy,0,2)-pt(ix,iy,0,2)
     :    +tsdif*(pts(ix+1,iy,0)-pts(ix,iy,0)))/dx
        !h1(ix,iy,0)=h1(ix,iy,0)*dif100(ix,iy,0)*tdif
         h1(ix,iy,0)=h1(ix,iy,0)*dif100_h(ix,iy,0)*tdif       !horizontal diffusivity is used
      enddo
      enddo
c
      do is=1,ns-1
      do iy=1,ny
      do ix=1,nx
        h2(ix,iy,is)=(pt(ix,iy+1,is,2)-pt(ix,iy,is,2)
     :    +tsdif*(pts(ix,iy+1,is)-pts(ix,iy,is)))/dy
     :    -s0ds4a(is)*ppdy01(ix,iy,2)*(pt(ix,iy+1,is+1,2)
     :    +pt(ix,iy,is+1,2)-pt(ix,iy+1,is-1,2)-pt(ix,iy,is-1,2)
     :    +tsdif*(pts(ix,iy+1,is+1)+pts(ix,iy,is+1)
     :    -pts(ix,iy+1,is-1)-pts(ix,iy,is-1)))
        !h2(ix,iy,is)=h2(ix,iy,is)*dif010(ix,iy,is)*tdif
         h2(ix,iy,is)=h2(ix,iy,is)*dif010_h(ix,iy,is)*tdif       !horizontal diffusivity is used
      enddo
      enddo
      enddo
c
      do iy=1,ny
      do ix=1,nx
c      h2(ix,iy,0)=-h2(ix,iy,1)
c      h2(ix,iy,ns)=-h2(ix,iy,ns-1)
        h2(ix,iy,ns)=(pt(ix,iy+1,ns,2)-pt(ix,iy,ns,2)
     :    +tsdif*(pts(ix,iy+1,ns)-pts(ix,iy,ns)))/dy
        h2(ix,iy,ns)=h2(ix,iy,ns)*dif010(ix,iy,ns)*tdif
        h2(ix,iy,0)=(pt(ix,iy+1,0,2)-pt(ix,iy,0,2)
     :    +tsdif*(pts(ix,iy+1,0)-pts(ix,iy,0)))/dy
        !h2(ix,iy,0)=h2(ix,iy,0)*dif010(ix,iy,0)*tdif
        h2(ix,iy,0)=h2(ix,iy,0)*dif010_h(ix,iy,0)*tdif             !horizontal diffusivity is used
      enddo
      enddo
c
      do is=2,ns-1
      do iy=1,ny
      do ix=1,nx
        h3(ix,iy,is)=-s(ix,iy,is)*(pt(ix,iy,is,2)-pt(ix,iy,is-1,2)
     :    +tsdif*(pts(ix,iy,is)-pts(ix,iy,is-1)))/ds0(is)
        h3(ix,iy,is)=h3(ix,iy,is)*0.5*(dif000(ix,iy,is)
     :    +dif000(ix,iy,is-1))*tdif
      enddo
      enddo
      enddo
c
      do iy=1,ny1
      do ix=1,nx1
        h3(ix,iy,1)=0.0
        !h3(ix,iy,ns)=0.0
        h3(ix,iy,ns)=ust_s(ix,iy)*tst_s(ix,iy)
        !write(0,*) ix,iy,h3(ix,iy,ns)
      enddo
      enddo
c
      do is=1,ns-1
      do iy=2,ny
      do ix=2,nx
        difunt(ix,iy,is)=(h1(ix,iy,is)-h1(ix-1,iy,is))/dx-s0ds4a(is)
     :    *ppdx(ix,iy,2)*(h1(ix,iy,is+1)+h1(ix-1,iy,is+1)-h1(ix,iy,is-1)
     :    -h1(ix-1,iy,is-1))
     :    +(h2(ix,iy,is)-h2(ix,iy-1,is))/dy-s0ds4a(is)
     :    *ppdy(ix,iy,2)*(h2(ix,iy,is+1)+h2(ix,iy-1,is+1)-h2(ix,iy,is-1)
     :    -h2(ix,iy-1,is-1))
     :    -(s(ix,iy,is+1)+s(ix,iy,is))*(h3(ix,iy,is+1)-h3(ix,iy,is))
     :    /ds12(is)
      enddo
      enddo
      enddo
c
c      do 3051 iy=2,ny
c      do 3051 ix=2,nx
c      difunt(ix,iy,ns-1)=difunt(ix,iy,ns-2)
c      difunt(ix,iy,1)=difunt(ix,iy,2)
c3051  continue
c
c
c      deallocate(h1)
c      deallocate(h2)
c      deallocate(h3)

      if(qif.gt.0.) then
c-
c difunq
c-
c      allocate(q1(0:nx1,0:ny1,0:ns1))
c      allocate(q2(0:nx1,0:ny1,0:ns1))
c      allocate(q3(0:nx1,0:ny1,0:ns1))

      do is=1,ns-1
      do iy=1,ny

      do ix=1,nx
        q1(ix,iy,is)=(qv(ix+1,iy,is,2)-qv(ix,iy,is,2)
     :    +qsdif*(qvs(ix+1,iy,is)-qvs(ix,iy,is)))/dx
     :    -s0ds4a(is)*ppdx10(ix,iy,2)*(qv(ix+1,iy,is+1,2)
     :    +qv(ix,iy,is+1,2)-qv(ix+1,iy,is-1,2)-qv(ix,iy,is-1,2)
     :    +qsdif*(qvs(ix+1,iy,is+1)+qvs(ix,iy,is+1)
     :    -qvs(ix+1,iy,is-1)-qvs(ix,iy,is-1)))
        q1(ix,iy,is)=q1(ix,iy,is)*dif100(ix,iy,is)*qdif
      enddo
      enddo
      enddo
c
      do iy=1,ny
      do ix=1,nx
c      q1(ix,iy,0)=-q1(ix,iy,1)
c      q1(ix,iy,ns)=-q1(ix,iy,ns-1)
        q1(ix,iy,ns)=(qv(ix+1,iy,ns,2)-qv(ix,iy,ns,2)
     :    +qsdif*(qvs(ix+1,iy,ns)-qvs(ix,iy,ns)))/dx
        q1(ix,iy,ns)=q1(ix,iy,ns)*dif100(ix,iy,ns)*qdif
        q1(ix,iy,0)=(qv(ix+1,iy,0,2)-qv(ix,iy,0,2)
     :    +qsdif*(qvs(ix+1,iy,0)-qvs(ix,iy,0)))/dx
        q1(ix,iy,0)=q1(ix,iy,0)*dif100(ix,iy,0)*qdif
      enddo
      enddo
c
      do is=1,ns-1
      do iy=1,ny
      do ix=1,nx
        q2(ix,iy,is)=(qv(ix,iy+1,is,2)-qv(ix,iy,is,2)
     :    +qsdif*(qvs(ix,iy+1,is)-qvs(ix,iy,is)))/dy
     :    -s0ds4a(is)*ppdy01(ix,iy,2)*(qv(ix,iy+1,is+1,2)
     :    +qv(ix,iy,is+1,2)-qv(ix,iy+1,is-1,2)-qv(ix,iy,is-1,2)
     :    +qsdif*(qvs(ix,iy+1,is+1)+qvs(ix,iy,is+1)
     :    -qvs(ix,iy+1,is-1)-qvs(ix,iy,is-1)))
        q2(ix,iy,is)=q2(ix,iy,is)*dif010(ix,iy,is)*qdif
      enddo
      enddo
      enddo
c
      do iy=1,ny
      do ix=1,nx
c     q2(ix,iy,0)=-q2(ix,iy,1)
c     q2(ix,iy,ns)=-q2(ix,iy,ns-1)
        q2(ix,iy,ns)=(qv(ix,iy+1,ns,2)-qv(ix,iy,ns,2)
     :    +qsdif*(qvs(ix,iy+1,ns)-qvs(ix,iy,ns)))/dy
        q2(ix,iy,ns)=q2(ix,iy,ns)*dif010(ix,iy,ns)*qdif
        q2(ix,iy,0)=(qv(ix,iy+1,0,2)-qv(ix,iy,0,2)
     :    +qsdif*(qvs(ix,iy+1,0)-qvs(ix,iy,0)))/dy
        q2(ix,iy,0)=q2(ix,iy,0)*dif010(ix,iy,0)*qdif
      enddo
      enddo
c
      do is=2,ns-1
      do iy=1,ny
      do ix=1,nx
        q3(ix,iy,is)=-s(ix,iy,is)*(qv(ix,iy,is,2)-qv(ix,iy,is-1,2)
     :    +qsdif*(qvs(ix,iy,is)-qvs(ix,iy,is-1)))/ds0(is)
        q3(ix,iy,is)=q3(ix,iy,is)*0.5*(dif000(ix,iy,is)
     :    +dif000(ix,iy,is-1))*qdif
      enddo
      enddo
      enddo
c
      do iy=1,ny1
      do ix=1,nx1
        q3(ix,iy,1)=0.0
        q3(ix,iy,ns)=0.0
      enddo
      enddo
c
      do is=1,ns-1
      do iy=2,ny
      do ix=2,nx
        difunq(ix,iy,is)=(q1(ix,iy,is)-q1(ix-1,iy,is))/dx-s0ds4a(is)
     :    *ppdx(ix,iy,2)*(q1(ix,iy,is+1)+q1(ix-1,iy,is+1)-q1(ix,iy,is-1)
     :    -q1(ix-1,iy,is-1))
     :    +(q2(ix,iy,is)-q2(ix,iy-1,is))/dy-s0ds4a(is)
     :    *ppdy(ix,iy,2)*(q2(ix,iy,is+1)+q2(ix,iy-1,is+1)-q2(ix,iy,is-1)
     :    -q2(ix,iy-1,is-1))
     :    -(s(ix,iy,is+1)+s(ix,iy,is))*(q3(ix,iy,is+1)-q3(ix,iy,is))
     :    /ds12(is)
      enddo
      enddo
      enddo
c
c      do 4051 iy=2,ny
c      do 4051 ix=2,nx
c      difunq(ix,iy,ns-1)=difunq(ix,iy,ns-2)
c      difunq(ix,iy,1)=difunq(ix,iy,2)
c4051  continue
c
      endif
c-
c difunqc
c-
      if(ifqc.gt.0.) then

      do is=1,ns-1
      do iy=1,ny
      do ix=1,nx
        q1(ix,iy,is)=(qc(ix+1,iy,is,2)-qc(ix,iy,is,2))/dx
     :     -s0ds4a(is)*ppdx10(ix,iy,2)*(qc(ix+1,iy,is+1,2)
     :     +qc(ix,iy,is+1,2)-qc(ix+1,iy,is-1,2)-qc(ix,iy,is-1,2))
        q1(ix,iy,is)=q1(ix,iy,is)*dif100(ix,iy,is)*qdif
      enddo
      enddo
      enddo
c
      do iy=1,ny
      do ix=1,nx
        q1(ix,iy,ns)=(qc(ix+1,iy,ns,2)-qc(ix,iy,ns,2))/dx
        q1(ix,iy,ns)=q1(ix,iy,ns)*dif100(ix,iy,ns)*qdif
        q1(ix,iy,0)=(qc(ix+1,iy,0,2)-qc(ix,iy,0,2))/dx
        q1(ix,iy,0)=q1(ix,iy,0)*dif100(ix,iy,0)*qdif
      enddo
      enddo
c
      do is=1,ns-1
      do iy=1,ny
      do ix=1,nx
        q2(ix,iy,is)=(qc(ix,iy+1,is,2)-qc(ix,iy,is,2))/dy
     :     -s0ds4a(is)*ppdy01(ix,iy,2)*(qc(ix,iy+1,is+1,2)
     :     +qc(ix,iy,is+1,2)-qc(ix,iy+1,is-1,2)-qc(ix,iy,is-1,2))
        q2(ix,iy,is)=q2(ix,iy,is)*dif010(ix,iy,is)*qdif
      enddo
      enddo
      enddo
c
      do iy=1,ny
      do ix=1,nx
        q2(ix,iy,ns)=(qc(ix,iy+1,ns,2)-qc(ix,iy,ns,2))/dy
        q2(ix,iy,ns)=q2(ix,iy,ns)*dif010(ix,iy,ns)*qdif
        q2(ix,iy,0)=(qc(ix,iy+1,0,2)-qc(ix,iy,0,2))/dy
        q2(ix,iy,0)=q2(ix,iy,0)*dif010(ix,iy,0)*qdif
      enddo
      enddo
c
      do is=2,ns-1
      do iy=1,ny
      do ix=1,nx
        q3(ix,iy,is)=-s(ix,iy,is)*(qc(ix,iy,is,2)-qc(ix,iy,is-1,2))
     :     /ds0(is)
        q3(ix,iy,is)=q3(ix,iy,is)*0.5*(dif000(ix,iy,is)
     :     +dif000(ix,iy,is-1))*qdif
      enddo
      enddo
      enddo
c
      do iy=1,ny1
      do ix=1,nx1
        q3(ix,iy,1)=0.0
        q3(ix,iy,ns)=0.0
      enddo
      enddo
c
      do is=1,ns-1
      do iy=2,ny
      do ix=2,nx
        difunqc(ix,iy,is)=(q1(ix,iy,is)-q1(ix-1,iy,is))/dx-s0ds4a(is)
     :    *ppdx(ix,iy,2)*(q1(ix,iy,is+1)+q1(ix-1,iy,is+1)-q1(ix,iy,is-1)
     :    -q1(ix-1,iy,is-1))
     :    +(q2(ix,iy,is)-q2(ix,iy-1,is))/dy-s0ds4a(is)
     :    *ppdy(ix,iy,2)*(q2(ix,iy,is+1)+q2(ix,iy-1,is+1)-q2(ix,iy,is-1)
     :    -q2(ix,iy-1,is-1))
     :    -(s(ix,iy,is+1)+s(ix,iy,is))*(q3(ix,iy,is+1)-q3(ix,iy,is))
     :    /ds12(is)
      enddo
      enddo
      enddo

      endif

c--------------------------------------------------------------------------------------
c  difunqci
c--------------------------------------------------------------------------------------
      if(ifqi.gt.0.) then                                                       DC,11.2009

      do is=1,ns-1                                                              DC,11.2009
      do iy=1,ny                                                                DC,11.2009
      do ix=1,nx                                                                DC,11.2009
        q1(ix,iy,is)=(qci(ix+1,iy,is,2)-qci(ix,iy,is,2))/dx                     DC,11.2009
     :     -s0ds4a(is)*ppdx10(ix,iy,2)*(qci(ix+1,iy,is+1,2)                     DC,11.2009
     :     +qci(ix,iy,is+1,2)-qci(ix+1,iy,is-1,2)-qci(ix,iy,is-1,2))            DC,11.2009
        q1(ix,iy,is)=q1(ix,iy,is)*dif100(ix,iy,is)*qdif                         DC,11.2009
      enddo                                                                     DC,11.2009
      enddo                                                                     DC,11.2009
      enddo                                                                     DC,11.2009
c                             
      do iy=1,ny                                                                DC,11.2009
      do ix=1,nx                                                                DC,11.2009
        q1(ix,iy,ns)=(qci(ix+1,iy,ns,2)-qci(ix,iy,ns,2))/dx                     DC,11.2009
        q1(ix,iy,ns)=q1(ix,iy,ns)*dif100(ix,iy,ns)*qdif                         DC,11.2009
        q1(ix,iy,0)=(qci(ix+1,iy,0,2)-qci(ix,iy,0,2))/dx                        DC,11.2009
        q1(ix,iy,0)=q1(ix,iy,0)*dif100(ix,iy,0)*qdif                            DC,11.2009
      enddo                                                                     DC,11.2009
      enddo                                                                     DC,11.2009
c
      do is=1,ns-1                                                              DC,11.2009
      do iy=1,ny                                                                DC,11.2009
      do ix=1,nx                                                                DC,11.2009
        q2(ix,iy,is)=(qci(ix,iy+1,is,2)-qci(ix,iy,is,2))/dy                     DC,11.2009
     :     -s0ds4a(is)*ppdy01(ix,iy,2)*(qci(ix,iy+1,is+1,2)                     DC,11.2009
     :     +qci(ix,iy,is+1,2)-qci(ix,iy+1,is-1,2)-qci(ix,iy,is-1,2))            DC,11.2009
        q2(ix,iy,is)=q2(ix,iy,is)*dif010(ix,iy,is)*qdif                         DC,11.2009
      enddo                                                                     DC,11.2009
      enddo                                                                     DC,11.2009
      enddo                                                                     DC,11.2009
c  
      do iy=1,ny                                                                DC,11.2009
      do ix=1,nx                                                                DC,11.2009
        q2(ix,iy,ns)=(qci(ix,iy+1,ns,2)-qci(ix,iy,ns,2))/dy                     DC,11.2009
        q2(ix,iy,ns)=q2(ix,iy,ns)*dif010(ix,iy,ns)*qdif                         DC,11.2009
        q2(ix,iy,0)=(qci(ix,iy+1,0,2)-qci(ix,iy,0,2))/dy                        DC,11.2009
        q2(ix,iy,0)=q2(ix,iy,0)*dif010(ix,iy,0)*qdif                            DC,11.2009
      enddo                                                                     DC,11.2009
      enddo                                                                     DC,11.2009
c 
      do is=2,ns-1                                                              DC,11.2009
      do iy=1,ny                                                                DC,11.2009
      do ix=1,nx                                                                DC,11.2009
        q3(ix,iy,is)=-s(ix,iy,is)*(qci(ix,iy,is,2)-qci(ix,iy,is-1,2))           DC,11.2009
     :     /ds0(is)                                                             DC,11.2009
        q3(ix,iy,is)=q3(ix,iy,is)*0.5*(dif000(ix,iy,is)                         DC,11.2009
     :     +dif000(ix,iy,is-1))*qdif                                            DC,11.2009
      enddo                                                                     DC,11.2009
      enddo                                                                     DC,11.2009
      enddo                                                                     DC,11.2009
c                                          
      do iy=1,ny1                                                               DC,11.2009
      do ix=1,nx1                                                               DC,11.2009
        q3(ix,iy,1)=0.0                                                         DC,11.2009
        q3(ix,iy,ns)=0.0                                                        DC,11.2009
      enddo                                                                     DC,11.2009
      enddo                                                                     DC,11.2009
c
      do is=1,ns-1                                                              DC,11.2009
      do iy=2,ny                                                                DC,11.2009
      do ix=2,nx                                                                DC,11.2009
        difunqci(ix,iy,is)=(q1(ix,iy,is)-q1(ix-1,iy,is))/dx-s0ds4a(is)          DC,11.2009
     :    *ppdx(ix,iy,2)*(q1(ix,iy,is+1)+q1(ix-1,iy,is+1)-q1(ix,iy,is-1)        DC,11.2009
     :    -q1(ix-1,iy,is-1))                                                    DC,11.2009
     :    +(q2(ix,iy,is)-q2(ix,iy-1,is))/dy-s0ds4a(is)                          DC,11.2009
     :    *ppdy(ix,iy,2)*(q2(ix,iy,is+1)+q2(ix,iy-1,is+1)-q2(ix,iy,is-1)        DC,11.2009
     :    -q2(ix,iy-1,is-1))                                                    DC,11.2009
     :    -(s(ix,iy,is+1)+s(ix,iy,is))*(q3(ix,iy,is+1)-q3(ix,iy,is))            DC,11.2009
     :    /ds12(is)                                                             DC,11.2009
      enddo                                                                     DC,11.2009
      enddo                                                                     DC,11.2009
      enddo                                                                     DC,11.2009

      endif                                                                     DC,11.2009

c--------------------------------------------------------------------------------------
c  difunqsn
c--------------------------------------------------------------------------------------
      if(ifqi.gt.0.) then                                                       DC,11.2009

      do is=1,ns-1                                                              DC,11.2009
      do iy=1,ny                                                                DC,11.2009
      do ix=1,nx                                                                DC,11.2009
        q1(ix,iy,is)=(qsn(ix+1,iy,is,2)-qsn(ix,iy,is,2))/dx                     DC,11.2009
     :     -s0ds4a(is)*ppdx10(ix,iy,2)*(qsn(ix+1,iy,is+1,2)                     DC,11.2009
     :     +qsn(ix,iy,is+1,2)-qsn(ix+1,iy,is-1,2)-qsn(ix,iy,is-1,2))            DC,11.2009
        q1(ix,iy,is)=q1(ix,iy,is)*dif100(ix,iy,is)*qdif                         DC,11.2009
      enddo                                                                     DC,11.2009
      enddo                                                                     DC,11.2009
      enddo                                                                     DC,11.2009
c                             
      do iy=1,ny                                                                DC,11.2009
      do ix=1,nx                                                                DC,11.2009
        q1(ix,iy,ns)=(qsn(ix+1,iy,ns,2)-qsn(ix,iy,ns,2))/dx                     DC,11.2009
        q1(ix,iy,ns)=q1(ix,iy,ns)*dif100(ix,iy,ns)*qdif                         DC,11.2009
        q1(ix,iy,0)=(qsn(ix+1,iy,0,2)-qsn(ix,iy,0,2))/dx                        DC,11.2009
        q1(ix,iy,0)=q1(ix,iy,0)*dif100(ix,iy,0)*qdif                            DC,11.2009
      enddo                                                                     DC,11.2009
      enddo                                                                     DC,11.2009
c
      do is=1,ns-1                                                              DC,11.2009
      do iy=1,ny                                                                DC,11.2009
      do ix=1,nx                                                                DC,11.2009
        q2(ix,iy,is)=(qsn(ix,iy+1,is,2)-qsn(ix,iy,is,2))/dy                     DC,11.2009
     :     -s0ds4a(is)*ppdy01(ix,iy,2)*(qsn(ix,iy+1,is+1,2)                     DC,11.2009
     :     +qsn(ix,iy,is+1,2)-qsn(ix,iy+1,is-1,2)-qsn(ix,iy,is-1,2))            DC,11.2009
        q2(ix,iy,is)=q2(ix,iy,is)*dif010(ix,iy,is)*qdif                         DC,11.2009
      enddo                                                                     DC,11.2009
      enddo                                                                     DC,11.2009
      enddo                                                                     DC,11.2009
c  
      do iy=1,ny                                                                DC,11.2009
      do ix=1,nx                                                                DC,11.2009
        q2(ix,iy,ns)=(qsn(ix,iy+1,ns,2)-qsn(ix,iy,ns,2))/dy                     DC,11.2009
        q2(ix,iy,ns)=q2(ix,iy,ns)*dif010(ix,iy,ns)*qdif                         DC,11.2009
        q2(ix,iy,0)=(qsn(ix,iy+1,0,2)-qsn(ix,iy,0,2))/dy                        DC,11.2009
        q2(ix,iy,0)=q2(ix,iy,0)*dif010(ix,iy,0)*qdif                            DC,11.2009
      enddo                                                                     DC,11.2009
      enddo                                                                     DC,11.2009
c 
      do is=2,ns-1                                                              DC,11.2009
      do iy=1,ny                                                                DC,11.2009
      do ix=1,nx                                                                DC,11.2009
        q3(ix,iy,is)=-s(ix,iy,is)*(qsn(ix,iy,is,2)-qsn(ix,iy,is-1,2))           DC,11.2009
     :     /ds0(is)                                                             DC,11.2009
        q3(ix,iy,is)=q3(ix,iy,is)*0.5*(dif000(ix,iy,is)                         DC,11.2009
     :     +dif000(ix,iy,is-1))*qdif                                            DC,11.2009
      enddo                                                                     DC,11.2009
      enddo                                                                     DC,11.2009
      enddo                                                                     DC,11.2009
c                                          
      do iy=1,ny1                                                               DC,11.2009
      do ix=1,nx1                                                               DC,11.2009
        q3(ix,iy,1)=0.0                                                         DC,11.2009
        q3(ix,iy,ns)=0.0                                                        DC,11.2009
      enddo                                                                     DC,11.2009
      enddo                                                                     DC,11.2009
c
      do is=1,ns-1                                                              DC,11.2009
      do iy=2,ny                                                                DC,11.2009
      do ix=2,nx                                                                DC,11.2009
        difunqsn(ix,iy,is)=(q1(ix,iy,is)-q1(ix-1,iy,is))/dx-s0ds4a(is)          DC,11.2009
     :    *ppdx(ix,iy,2)*(q1(ix,iy,is+1)+q1(ix-1,iy,is+1)-q1(ix,iy,is-1)        DC,11.2009
     :    -q1(ix-1,iy,is-1))                                                    DC,11.2009
     :    +(q2(ix,iy,is)-q2(ix,iy-1,is))/dy-s0ds4a(is)                          DC,11.2009
     :    *ppdy(ix,iy,2)*(q2(ix,iy,is+1)+q2(ix,iy-1,is+1)-q2(ix,iy,is-1)        DC,11.2009
     :    -q2(ix,iy-1,is-1))                                                    DC,11.2009
     :    -(s(ix,iy,is+1)+s(ix,iy,is))*(q3(ix,iy,is+1)-q3(ix,iy,is))            DC,11.2009
     :    /ds12(is)                                                             DC,11.2009
      enddo                                                                     DC,11.2009
      enddo                                                                     DC,11.2009
      enddo                                                                     DC,11.2009

      endif                                                                     DC,11.2009

c-
c difunqr
c-
      if(ifqr.gt.0.) then

      do is=1,ns-1
      do iy=1,ny
      do ix=1,nx
        q1(ix,iy,is)=(qr(ix+1,iy,is,2)-qr(ix,iy,is,2))/dx
     :     -s0ds4a(is)*ppdx10(ix,iy,2)*(qr(ix+1,iy,is+1,2)
     :     +qr(ix,iy,is+1,2)-qr(ix+1,iy,is-1,2)-qr(ix,iy,is-1,2))
        q1(ix,iy,is)=q1(ix,iy,is)*dif100(ix,iy,is)*qdif
      enddo
      enddo
      enddo
c
      do iy=1,ny
      do ix=1,nx
        q1(ix,iy,ns)=(qr(ix+1,iy,ns,2)-qr(ix,iy,ns,2))/dx
        q1(ix,iy,ns)=q1(ix,iy,ns)*dif100(ix,iy,ns)*qdif
        q1(ix,iy,0)=(qr(ix+1,iy,0,2)-qr(ix,iy,0,2))/dx
        q1(ix,iy,0)=q1(ix,iy,0)*dif100(ix,iy,0)*qdif
      enddo
      enddo
c
      do is=1,ns-1
      do iy=1,ny
      do ix=1,nx
        q2(ix,iy,is)=(qr(ix,iy+1,is,2)-qr(ix,iy,is,2))/dy
     :     -s0ds4a(is)*ppdy01(ix,iy,2)*(qr(ix,iy+1,is+1,2)
     :     +qr(ix,iy,is+1,2)-qr(ix,iy+1,is-1,2)-qr(ix,iy,is-1,2))
        q2(ix,iy,is)=q2(ix,iy,is)*dif010(ix,iy,is)*qdif
      enddo
      enddo
      enddo
c
      do iy=1,ny
      do ix=1,nx
        q2(ix,iy,ns)=(qr(ix,iy+1,ns,2)-qr(ix,iy,ns,2))/dy
        q2(ix,iy,ns)=q2(ix,iy,ns)*dif010(ix,iy,ns)*qdif
        q2(ix,iy,0)=(qr(ix,iy+1,0,2)-qr(ix,iy,0,2))/dy
        q2(ix,iy,0)=q2(ix,iy,0)*dif010(ix,iy,0)*qdif
      enddo
      enddo
c
      do is=2,ns-1
      do iy=1,ny
      do ix=1,nx
        q3(ix,iy,is)=-s(ix,iy,is)*(qr(ix,iy,is,2)-qr(ix,iy,is-1,2))
     :     /ds0(is)
        q3(ix,iy,is)=q3(ix,iy,is)*0.5*(dif000(ix,iy,is)
     :     +dif000(ix,iy,is-1))*qdif
      enddo
      enddo
      enddo
c
      do iy=1,ny1
      do ix=1,nx1
        q3(ix,iy,1)=0.0
        q3(ix,iy,ns)=0.0
      enddo
      enddo
c
      do is=1,ns-1
      do iy=2,ny
      do ix=2,nx
        difunqr(ix,iy,is)=(q1(ix,iy,is)-q1(ix-1,iy,is))/dx-s0ds4a(is)
     :    *ppdx(ix,iy,2)*(q1(ix,iy,is+1)+q1(ix-1,iy,is+1)-q1(ix,iy,is-1)
     :    -q1(ix-1,iy,is-1))
     :    +(q2(ix,iy,is)-q2(ix,iy-1,is))/dy-s0ds4a(is)
     :    *ppdy(ix,iy,2)*(q2(ix,iy,is+1)+q2(ix,iy-1,is+1)-q2(ix,iy,is-1)
     :    -q2(ix,iy-1,is-1))
     :    -(s(ix,iy,is+1)+s(ix,iy,is))*(q3(ix,iy,is+1)-q3(ix,iy,is))
     :    /ds12(is)
      enddo
      enddo
      enddo

      endif

c      if(qif.ne.0.) then
c        deallocate(q1)
c        deallocate(q2)
c        deallocate(q3)
c      endif

c
c debug:
c
      if(iodify.eq.0) then
        call extray(nx1,ny1,difunu,1,nx,2,ny,0,ns)
        call extray(nx1,ny1,difunv,2,nx,2,ny-1,0,ns)
        call extray(nx1,ny1,difunw,1,nx1,2,ny,0,ns)
      endif
c
c these extrapolations get rid of unused boundary values for diffusion:
c
      call extrah(nx1,ny1,difunu,1,nx,1,ny1,0,ns)
      call extrah(nx1,ny1,difunv,1,nx1,1,ny,0,ns)
      call extrah(nx1,ny1,difunw,1,nx1,1,ny1,1,ns)
c
c      deallocate(dif000)
c      deallocate(dif010)
c      deallocate(dif100)

      return
      end

      Subroutine diffuv(fngrd
     :  ,def11,h1,q1
     :  ,def12,h2,q2
     :  ,def13,h3,q3
     :  ,def22
     :  ,def23
     :  ,def33
     :  ,dif000
     :  ,dif010
     :  ,dif001,dif100)

      use alloc

      implicit real*8 (a-h,o-z)
      real(kind(0.d0)),dimension(0:nx1,0:ny1,0:ns1)::
     :  def11,def12,def13,def22,def23,def33
     :  ,h1,h2,h3,q1,q2,q3,dif000,dif010,dif001,dif100
      
c-----------------------------------------------------------------------
c  compute vertical diffusion terms which are deformation and richardson number
c  dependent after lilly. if rikey=0.0,the later dependency is dropped.
c  enhanced upper-level dampping is determined by dragt,dragm,etc.
c
c  difunu,difunv and difunw are multiplied for pp, to save in ellipt.
c
c rayleigh damping passed to prognostic equations.
c
c-----------------------------------------------------------------------
c
      character*80 fngrd

      parameter (third4=4./3.,third2=2./3.)

c
c calculation of the deformation tensor:
c

      if(verbose.ge.2) write(*,*) 'enter diffu-vertical'

      do is=1,ns-1
      do iy=1,ny1
      do ix=1,nx1
        zdwdz=-(s(ix,iy,is+1)+s(ix,iy,is))
     :    *(w(ix,iy,is+1,2)-w(ix,iy,is,2))/ds12(is)
        def11(ix,iy,is)=-third2*zdwdz
        def22(ix,iy,is)=-third2*zdwdz
        def33(ix,iy,is)=third4*zdwdz
      enddo
      enddo
      enddo
c
c boundary conditions d/ds=0
c
      do is=0,ns,ns
      do iy=1,ny1
      do ix=1,nx1
        zdwdz=-(s(ix,iy,is+1)+s(ix,iy,is))
     :    *(w(ix,iy,is+1,2)-w(ix,iy,is,2))/ds12(is)
        def11(ix,iy,is)=-third2*zdwdz
        def22(ix,iy,is)=-third2*zdwdz
        def33(ix,iy,is)=third4*zdwdz
      enddo
      enddo
      enddo
c
      do is=1,ns-1
      do iy=1,ny
      do ix=1,nx
        def12(ix,iy,is)=0.
      enddo
      enddo
      enddo
c
      do is=0,ns,ns
      do iy=1,ny
      do ix=1,nx
        def12(ix,iy,is)=0.
c     :    (u(ix,iy+1,is,2)-u(ix,iy,is,2))/dy
c     :    +(v(ix+1,iy,is,2)-v(ix,iy,is,2))/dx
      enddo
      enddo
      enddo
c
      do is=1,ns
      do iy=1,ny1
      do ix=1,nx
        def13(ix,iy,is)=-(s(ix+1,iy,is)+s(ix,iy,is))
     :    *(u(ix,iy,is,2)-u(ix,iy,is-1,2))/ds02(is)
      enddo
      enddo
      enddo
c
      do is=1,ns
      do iy=1,ny
      do ix=1,nx1
        def23(ix,iy,is)=-(s(ix,iy,is)+s(ix,iy+1,is))
     :    *(v(ix,iy,is,2)-v(ix,iy,is-1,2))/ds02(is)
      enddo
      enddo
      enddo
c
      if(defout.and.prt) then
        tk=0.
        call wri3ar(def11,1,nx1,1,ny1,1,ns-1,'def11   ',tk,iomod,nx,ny)
        call wri3ar(def12,1,nx,1,ny,1,ns-1,'def12   ',tk,iomod,nx,ny)
        call wri3ar(def13,1,nx,1,ny1,1,ns-1,'def13   ',tk,iomod,nx,ny)
        call wri3ar(def22,1,nx1,1,ny1,1,ns-1,'def22   ',tk,iomod,nx,ny)
        call wri3ar(def23,1,nx1,1,ny,1,ns-1,'def23   ',tk,iomod,nx,ny)
        call wri3ar(def33,1,nx1,1,ny1,1,ns-1,'def33   ',tk,iomod,nx,ny)
      endif
c
c calculation of diffusion coefficients:
c
      if(iodif.eq.1) then
c-
c calculate linear diffusion coefficients.
c-
        do is=0,ns
        do iy=1,ny1
        do ix=1,nx1
          dif000(ix,iy,is)=difcof
        enddo
        enddo
        enddo
      else
c-
c calculate deformation/richardson number dependent diffusion coefficient
c-

        if(prt.and.driout) allocate(ri(0:nx1,0:ny1,0:ns1))
	if(prt.and.driout) allocate(bosa(0:nx1,0:ny1,0:ns1))
	  do is=1,ns-1
        do iy=2,ny
        do ix=2,nx
	
c
          def=dsqrt(0.5*(def11(ix,iy,is)*def11(ix,iy,is)+def22(ix,iy,is)
     :      *def22(ix,iy,is)+def33(ix,iy,is)*def33(ix,iy,is))
     :      +0.25*(def12(ix,iy,is)*def12(ix,iy,is)
     :      +def12(ix-1,iy,is)*def12(ix-1,iy,is)
     :      +def12(ix,iy-1,is)*def12(ix,iy-1,is)
     :      +def12(ix-1,iy-1,is)*def12(ix-1,iy-1,is)
     :      +def13(ix,iy,is)*def13(ix,iy,is)
     :      +def13(ix-1,iy,is)*def13(ix-1,iy,is)
     :      +def13(ix,iy,is+1)*def13(ix,iy,is+1)
     :      +def13(ix-1,iy,is+1)*def13(ix-1,iy,is+1)
     :      +def23(ix,iy,is)*def23(ix,iy,is)
     :      +def23(ix,iy-1,is)*def23(ix,iy-1,is)
     :      +def23(ix,iy,is+1)*def23(ix,iy,is+1)
     :      +def23(ix,iy-1,is+1)*def23(ix,iy-1,is+1)))
	    if(prt.and.driout) bosa(ix,iy,is)=def
c
               
          if(qif.ne.0. .and. ifqc.ne.0) then
            p=sigma0(is)*pp(ix,iy,2)+ptop
            p2=sigma0(is+1)*pp(ix,iy,2)+ptop
            p1=sigma0(is-1)*pp(ix,iy,2)+ptop
            t=(pts(ix,iy,is)+pt(ix,iy,is,2))*(p/p00)**akapa
       
           tp2=(pts(ix,iy,is+1)+pt(ix,iy,is+1,2))*(p2/p00)**akapa
            tp1=(pts(ix,iy,is-1)+pt(ix,iy,is-1,2))*(p1/p00)**akapa
	      qsatur=qsat(t,p)
            qsat2=qsat(tp2,p2)
            qsat1=qsat(tp1,p1)
            xnsq=-g*(s(ix,iy,is+1)+s(ix,iy,is))*
     :      ((1+hlat*qsatur/(rv*t))/(1+0.622*hlat*hlatcp*qsatur/
     :      (rv*t*t))
     :      *((pt(ix,iy,is+1,2)
     :      +pts(ix,iy,is+1)-pt(ix,iy,is-1,2)-pts(ix,iy,is-1))
     :      /(pt(ix,iy,is,2)+pts(ix,iy,is))+hlatcp/t*(qsat2-qsat1))-
     :      (qv(ix,iy,is+1,2)+qc(ix,iy,is+1,2)-qv(ix,iy,is-1,2)-
     :      qc(ix,iy,is-1,2)))
     :      /ds04a(is)
          else
            xnsq=-g*(s(ix,iy,is+1)+s(ix,iy,is))*(pt(ix,iy,is+1,2)
     :        +pts(ix,iy,is+1)-pt(ix,iy,is-1,2)-pts(ix,iy,is-1))
     :        /(pt(ix,iy,is,2)+pts(ix,iy,is))
     :        /ds04a(is)
          endif
          rich=xnsq/(def*def+zero0)
          if(prt.and.driout) ri(ix,iy,is)=rich
	    !write(0,*) ix,iy,is, deltaz(ix,iy)
c
          difk=delka2(is)*def*dsqrt(dim(1.0d0,rich*rikey))

	    
c
          dif000(ix,iy,is)=difk+difl
        enddo
        enddo
        enddo
c
        call extra3(nx1,ny1,ns1,dif000,1,nx1,1,ny1,0,ns)
c
        if(prt.and.driout) then
          if(grdout) then
            call wgrids(ri,'ri',2,nx,2,ny,2,ns-1,iomod
     :        ,0,0,0,fngrd)
            call wgrids(dif000,'kk',2,nx,2,ny,2,ns-1,iomod
     :        ,0,0,0,fngrd)
	      call wgrids(bosa,'de',2,nx,2,ny,2,ns-1,iomod
     :        ,0,0,0,fngrd)
          else
            tk=0.
            call wri3ar(ri,2,nx,2,ny,2,ns-1,'ri      ',tk,iomod,nx,ny)
            call wri3ar(dif000,2,nx,2,ny,2,ns-1,'dif000  ',tk,iomod
     :        ,nx,ny)
          endif
          deallocate(ri)
	    deallocate(bosa)
        endif
      endif
c
c top boundary enhanced diffusion:
c
      if(olddra) then
        ixm=nx1/2
        iym=ny1/2
        ztop=phis(ixm,iym,1)/g
        do is=0,ns
        do iy=1,ny1
        do ix=1,nx1
          z=(phi(ix,iy,is)+phis(ix,iy,is))/g
          drag=(dragm+(dragt-dragm)*cos(pi05*abs((z-ztop)/(zm-ztop)))
     :      **2)*dim(z,zm)/(z-zm+zero0)
          dif000(ix,iy,is)=dif000(ix,iy,is)+drag
	  enddo
        enddo
        enddo
      elseif(.not.raylei) then
         do is=0,idrmax
         do iy=1,ny1
         do ix=1,nx1
           dif000(ix,iy,is)=dif000(ix,iy,is)+endrag(is)
	    enddo
         enddo
         enddo
      endif
c
c
c impose linear stability for diffusion coefficient
c
c      do is=0,ns
c      do iy=1,ny1
c      do ix=1,nx1
c        dif000(ix,iy,is)=min(deltaz(ix,iy)**2/(8.*dt),dif000(ix,iy,is))
c      enddo
c      enddo
c      enddo

      if(prt.and.ddiout) then
         tk=0.
         call wri3ar(dif000,1,nx1,1,ny1,0,ns,'difcof  ',tk,iomod,nx,ny)
      endif
c
c the upper layer drag is non-linear. so maybe it should be calculated
c separatelly for the other grid.
c
c diffusion coefficients on other grid locations:
c

c      allocate(dif010(0:nx1,0:ny1,0:ns1))
c      allocate(dif100(0:nx1,0:ny1,0:ns1))

      do is=0,ns
      do iy=1,ny
      do ix=1,nx
        dif010(ix,iy,is)=0.5*(dif000(ix,iy,is)+dif000(ix,iy+1,is))
        dif100(ix,iy,is)=0.5*(dif000(ix,iy,is)+dif000(ix+1,iy,is))
      enddo
      enddo
      enddo
c
      call extrah(nx1,ny1,dif010,0,nx1,0,ny1,0,ns)
      call extrah(nx1,ny1,dif100,0,nx1,0,ny1,0,ns)
c
c diffuse only the perturbated field:
c
      if(uvdif.ne.1.0) then
c
         do is=0,ns
         do iy=1,ny1
         do ix=1,nx1
           def33(ix,iy,is)=def33(ix,iy,is)*dif000(ix,iy,is)
         enddo
         enddo
         enddo
c
         do is=1,ns
         do iy=1,ny1
         do ix=1,nx
           def13s=-(s(ix+1,iy,is)+s(ix,iy,is))
     :       *(us(ix,iy,is)-us(ix,iy,is-1))/ds02(is)
           def13(ix,iy,is)=(def13(ix,iy,is)-(1.0-uvdif)*def13s)
     :       *0.5*(dif100(ix,iy,is)+dif100(ix,iy,is-1))
         enddo
         enddo
         enddo
c
         do is=1,ns
         do iy=1,ny
         do ix=1,nx1
           def23s=-(s(ix,iy+1,is)+s(ix,iy,is))
     :       *(vs(ix,iy,is)-vs(ix,iy,is-1))/ds02(is)
           def23(ix,iy,is)=(def23(ix,iy,is)-(1.0-uvdif)*def23s)
     :       *0.5*(dif010(ix,iy,is)+dif010(ix,iy,is-1))
         enddo
         enddo
         enddo
c
      else
         do is=0,ns
         do iy=1,ny1
         do ix=1,nx1
           def33(ix,iy,is)=def33(ix,iy,is)*dif000(ix,iy,is)
         enddo
         enddo
         enddo
c
         do is=1,ns
         do iy=1,ny1
         do ix=1,nx
           def13(ix,iy,is)=def13(ix,iy,is)*0.5*(dif100(ix,iy,is)
     :       +dif100(ix,iy,is-1))
         enddo
         enddo
         enddo
c
         do is=1,ns
         do iy=1,ny
         do ix=1,nx1
           def23(ix,iy,is)=def23(ix,iy,is)*0.5*(dif010(ix,iy,is)
     :       +dif010(ix,iy,is-1))
         enddo
         enddo
         enddo

	do iy=1,ny1
	do ix=1,nx
	def13(ix,iy,ns)=cdm(ix,iy)*u(ix,iy,ns-1,2)
	enddo
	enddo
	do iy=1,ny
         do ix=1,nx1
      def23(ix,iy,ns)=cdm(ix,iy)*v(ix,iy,ns-1,2)
	enddo
	enddo
c
      endif
c--
c difunu
c--
      do is=1,ns-1
      do iy=2,ny
      do ix=1,nx
        difunu(ix,iy,is)=(
     :    -(s(ix+1,iy,is+1)+s(ix+1,iy,is)+s(ix,iy,is+1)+s(ix,iy,is))
     :    *(def13(ix,iy,is+1)-def13(ix,iy,is))/ds14(is))*pp10(ix,iy,2)
      enddo
      enddo
      enddo
c
      do iy=2,ny
      do ix=1,nx
        difunu(ix,iy,0)=difunu(ix,iy,1)
        difunu(ix,iy,ns)=difunu(ix,iy,ns-1)
      enddo
      enddo
c
c--
c difunv
c--
      do is=1,ns-1
      do iy=1,ny
      do ix=2,nx
        difunv(ix,iy,is)=(
     :    -(s(ix,iy+1,is+1)+s(ix,iy+1,is)+s(ix,iy,is+1)+s(ix,iy,is))
     :    *(def23(ix,iy,is+1)-def23(ix,iy,is))/ds14(is))*pp01(ix,iy,2)
      enddo
      enddo
      enddo
c
      do iy=1,ny
      do ix=2,nx
        difunv(ix,iy,0)=difunv(ix,iy,1)
        difunv(ix,iy,ns)=difunv(ix,iy,ns-1)
      enddo
      enddo

c--
c  difunw
c--
      do is=2,ns-1
      do iy=2,ny
      do ix=2,nx
        difunw(ix,iy,is)=
     :   ( -s(ix,iy,is)*(def33(ix,iy,is)-def33(ix,iy,is-1))/ds0(is))
     :    *pp(ix,iy,2)
      enddo
      enddo
      enddo
c
      do iy=2,ny
      do ix=2,nx
        difunw(ix,iy,ns)=difunw(ix,iy,ns-1)
        difunw(ix,iy,1)=difunw(ix,iy,2)
      enddo
      enddo
c-
c difunt
c-

 !     do is=2,ns-1
 !     do iy=1,ny
 !     do ix=1,nx
 !       h3(ix,iy,is)=-s(ix,iy,is)*(pt(ix,iy,is,2)-pt(ix,iy,is-1,2)
 !    :    +tsdif*(pts(ix,iy,is)-pts(ix,iy,is-1)))/ds0(is)
 !       h3(ix,iy,is)=h3(ix,iy,is)*0.5*(dif000(ix,iy,is)
 !    :    +dif000(ix,iy,is-1))*tdif
 !     enddo
 !     enddo
 !     enddo

      do is=2,ns-1
      do iy=1,ny
      do ix=1,nx
        h3(ix,iy,is)=-s(ix,iy,is)*(pt(ix,iy,is,2)-pt(ix,iy,is-1,2)
     :    +tsdif*(pts(ix,iy,is)-pts(ix,iy,is-1)))/ds0(is)
        h3(ix,iy,is)=h3(ix,iy,is)*0.5*(dif000(ix,iy,is)
     :    +dif000(ix,iy,is-1))*tdif
      enddo
      enddo
      enddo
c
      do iy=1,ny1
      do ix=1,nx1
        h3(ix,iy,1)=0.0
        !h3(ix,iy,ns)=0.0
        h3(ix,iy,ns)=-ust_s(ix,iy)*tst_s(ix,iy)
      enddo
      enddo
c
	if(prt.and.driout) then
          if(grdout) then
            call wgrids(h3,'h3',2,nx,2,ny,2,ns-1,iomod
     :        ,0,0,0,fngrd)
	    endif
	endif
c
c
      do is=1,ns-1
      do iy=2,ny
      do ix=2,nx
        difunt(ix,iy,is)=
     :    -(s(ix,iy,is+1)+s(ix,iy,is))*(h3(ix,iy,is+1)-h3(ix,iy,is))
     :    /ds12(is)
      enddo
      enddo
      enddo

c
c     Dima@ 13.04.07
c-
c difunq_10
c-     
      if (impur.gt.0) then
      do is=2,ns-1                                                              DM,05.2007
      do iy=1,ny                                                                DM,05.2007
      do ix=1,nx                                                                DM,05.2007
        q3(ix,iy,is)=-s(ix,iy,is)*(qs_10(ix,iy,is,2)-qs_10(ix,iy,is-1,2)        DM,05.2007
     :    )/ds0(is) ! +qsdif*(qvs(ix,iy,is)-qvs(ix,iy,is-1)))                   DM,05.2007  
        q3(ix,iy,is)=q3(ix,iy,is)*0.5*(dif000(ix,iy,is)                         DM,05.2007
     :    +dif000(ix,iy,is-1))*qdif                                             DM,05.2007 
      enddo                                                                     DM,05.2007
      enddo                                                                     DM,05.2007
      enddo                                                                     DM,05.2007 
c
      do iy=1,ny1                                                               DM,05.2007
      do ix=1,nx1                                                               DM,05.2007
        q3(ix,iy,1)=0.0                                                         DM,05.2007
        q3(ix,iy,ns)=0.0                                                        DM,05.2007
      enddo                                                                     DM,05.2007
      enddo                                                                     DM,05.2007
      
      do is=1,ns-1                                                              DM,05.2007
      do iy=2,ny                                                                DM,05.2007
      do ix=2,nx                                                                DM,05.2007
        difunq_10(ix,iy,is)=                                                    DM,05.2007
     :    -(s(ix,iy,is+1)+s(ix,iy,is))*(q3(ix,iy,is+1)-q3(ix,iy,is))            DM,05.2007
     :    /ds12(is)                                                             DM,05.2007
      enddo                                                                     DM,05.2007
      enddo                                                                     DM,05.2007
      enddo                                                                     DM,05.2007
	endif

      if(qif.gt.0.) then
c-
c difunq
c-
      do is=2,ns-1
      do iy=1,ny
      do ix=1,nx
        q3(ix,iy,is)=-s(ix,iy,is)*(qv(ix,iy,is,2)-qv(ix,iy,is-1,2)
     :    +qsdif*(qvs(ix,iy,is)-qvs(ix,iy,is-1)))/ds0(is)
        q3(ix,iy,is)=q3(ix,iy,is)*0.5*(dif000(ix,iy,is)
     :    +dif000(ix,iy,is-1))*qdif
      enddo
      enddo
      enddo
c
      do iy=1,ny1
      do ix=1,nx1
        q3(ix,iy,1)=0.0
        !q3(ix,iy,ns)=0.0
        q3(ix,iy,ns)=-ust_s(ix,iy)*qst_s(ix,iy)
      enddo
      enddo
c
      do is=1,ns-1
      do iy=2,ny
      do ix=2,nx
        difunq(ix,iy,is)=
     :    -(s(ix,iy,is+1)+s(ix,iy,is))*(q3(ix,iy,is+1)-q3(ix,iy,is))
     :    /ds12(is)
      enddo
      enddo
      enddo
c
      endif
c-
c difunqc
c-
      if(ifqc.gt.0.) then

c
      do is=2,ns-1
      do iy=1,ny
      do ix=1,nx
        q3(ix,iy,is)=-s(ix,iy,is)*(qc(ix,iy,is,2)-qc(ix,iy,is-1,2))
     :     /ds0(is)
        q3(ix,iy,is)=q3(ix,iy,is)*0.5*(dif000(ix,iy,is)
     :     +dif000(ix,iy,is-1))*qdif
      enddo
      enddo
      enddo
c
      do iy=1,ny1
      do ix=1,nx1
        q3(ix,iy,1)=0.0
        q3(ix,iy,ns)=0.0
      enddo
      enddo
c
      do is=1,ns-1
      do iy=2,ny
      do ix=2,nx
        difunqc(ix,iy,is)=
     :    -(s(ix,iy,is+1)+s(ix,iy,is))*(q3(ix,iy,is+1)-q3(ix,iy,is))
     :    /ds12(is)
      enddo
      enddo
      enddo

      endif

c-----------------------------------------------------------------------------
c difunqci
c-----------------------------------------------------------------------------
      if(ifqi.gt.0.) then                                                       DC,11.2009

c
      do is=2,ns-1                                                              DC,11.2009
      do iy=1,ny                                                                DC,11.2009
      do ix=1,nx                                                                DC,11.2009
        q3(ix,iy,is)=-s(ix,iy,is)*(qci(ix,iy,is,2)-qci(ix,iy,is-1,2))           DC,11.2009
     :     /ds0(is)                                                             DC,11.2009
        q3(ix,iy,is)=q3(ix,iy,is)*0.5*(dif000(ix,iy,is)                         DC,11.2009
     :     +dif000(ix,iy,is-1))*qdif                                            DC,11.2009
      enddo                                                                     DC,11.2009
      enddo                                                                     DC,11.2009
      enddo                                                                     DC,11.2009
c                  
      do iy=1,ny1                                                               DC,11.2009
      do ix=1,nx1                                                               DC,11.2009
        q3(ix,iy,1)=0.0                                                         DC,11.2009
        q3(ix,iy,ns)=0.0                                                        DC,11.2009
      enddo                                                                     DC,11.2009
      enddo                                                                     DC,11.2009
c
      do is=1,ns-1                                                              DC,11.2009
      do iy=2,ny                                                                DC,11.2009
      do ix=2,nx                                                                DC,11.2009
        difunqci(ix,iy,is)=                                                     DC,11.2009
     :    -(s(ix,iy,is+1)+s(ix,iy,is))*(q3(ix,iy,is+1)-q3(ix,iy,is))            DC,11.2009
     :    /ds12(is)                                                             DC,11.2009
      enddo                                                                     DC,11.2009
      enddo                                                                     DC,11.2009
      enddo                                                                     DC,11.2009
                             
      endif                                                                     DC,11.2009

c-----------------------------------------------------------------------------
c difunqsn
c-----------------------------------------------------------------------------
      if(ifqi.gt.0.) then                                                       DC,11.2009

c
      do is=2,ns-1                                                              DC,11.2009
      do iy=1,ny                                                                DC,11.2009
      do ix=1,nx                                                                DC,11.2009
        q3(ix,iy,is)=-s(ix,iy,is)*(qsn(ix,iy,is,2)-qsn(ix,iy,is-1,2))           DC,11.2009
     :     /ds0(is)                                                             DC,11.2009
        q3(ix,iy,is)=q3(ix,iy,is)*0.5*(dif000(ix,iy,is)                         DC,11.2009
     :     +dif000(ix,iy,is-1))*qdif                                            DC,11.2009
      enddo                                                                     DC,11.2009
      enddo                                                                     DC,11.2009
      enddo                                                                     DC,11.2009
c                  
      do iy=1,ny1                                                               DC,11.2009
      do ix=1,nx1                                                               DC,11.2009
        q3(ix,iy,1)=0.0                                                         DC,11.2009
        q3(ix,iy,ns)=0.0                                                        DC,11.2009
      enddo                                                                     DC,11.2009
      enddo                                                                     DC,11.2009
c
      do is=1,ns-1                                                              DC,11.2009
      do iy=2,ny                                                                DC,11.2009
      do ix=2,nx                                                                DC,11.2009
        difunqsn(ix,iy,is)=                                                     DC,11.2009
     :    -(s(ix,iy,is+1)+s(ix,iy,is))*(q3(ix,iy,is+1)-q3(ix,iy,is))            DC,11.2009
     :    /ds12(is)                                                             DC,11.2009
      enddo                                                                     DC,11.2009
      enddo                                                                     DC,11.2009
      enddo                                                                     DC,11.2009
                             
      endif                                                                     DC,11.2009

      
c-
c difunqr
c-
      if(ifqr.gt.0.) then

      do is=2,ns-1
      do iy=1,ny
      do ix=1,nx
        q3(ix,iy,is)=-s(ix,iy,is)*(qr(ix,iy,is,2)-qr(ix,iy,is-1,2))
     :     /ds0(is)
        q3(ix,iy,is)=q3(ix,iy,is)*0.5*(dif000(ix,iy,is)
     :     +dif000(ix,iy,is-1))*qdif
      enddo
      enddo
      enddo
c
      do iy=1,ny1
      do ix=1,nx1
        q3(ix,iy,1)=0.0
        q3(ix,iy,ns)=0.0
      enddo
      enddo
c
      do is=1,ns-1
      do iy=2,ny
      do ix=2,nx
        difunqr(ix,iy,is)=
     :    -(s(ix,iy,is+1)+s(ix,iy,is))*(q3(ix,iy,is+1)-q3(ix,iy,is))
     :    /ds12(is)
      enddo
      enddo
      enddo

      endif

c
c debug:
c
      if(iodify.eq.0) then
        call extray(nx1,ny1,difunu,1,nx,2,ny,0,ns)
        call extray(nx1,ny1,difunv,2,nx,2,ny-1,0,ns)
        call extray(nx1,ny1,difunw,1,nx1,2,ny,0,ns)
      endif
c
c these extrapolations get rid of unused boundary values for diffusion:
c
      call extrah(nx1,ny1,difunu,1,nx,1,ny1,0,ns)
      call extrah(nx1,ny1,difunv,1,nx1,1,ny,0,ns)
      call extrah(nx1,ny1,difunw,1,nx1,1,ny1,1,ns)
c
      return
      end

      Subroutine diffuvt(fngrd
     :  ,def11,h1,q1
     :  ,def12,h2,q2
     :  ,def13,h3,q3
     :  ,def22
     :  ,def23
     :  ,def33
     :  ,dif000
     :  ,dif010
     :  ,dif001,dif100)

      use alloc

      implicit real*8 (a-h,o-z)
      real(kind(0.d0)),dimension(0:nx1,0:ny1,0:ns1)::
     :  def11,def12,def13,def22,def23,def33
     :  ,h1,h2,h3,q1,q2,q3,dif000,dif010,dif001,dif100
	
c-----------------------------------------------------------------------
c  compute vertical diffusion terms which are deformation and richardson number
c  dependent after lilly. if rikey=0.0,the later dependency is dropped.
c  enhanced upper-level dampping is determined by dragt,dragm,etc.
c
c  difunu,difunv and difunw are multiplied for pp, to save in ellipt.
c
c rayleigh damping passed to prognostic equations.
c
c teixeira's method (jcp, 1999)
c
c-----------------------------------------------------------------------
c

      character*80 fngrd

      parameter (third4=4./3.,third2=2./3.)

c
c calculation of the deformation tensor:
c

      if(verbose.ge.1) write(*,*) 'enter diffu-vertical teixeira'

      do is=1,ns-1
      do iy=1,ny1
      do ix=1,nx1
        zdwdz=-(s(ix,iy,is+1)+s(ix,iy,is))
     :    *(w(ix,iy,is+1,2)-w(ix,iy,is,2))/ds12(is)
        def11(ix,iy,is)=-third2*zdwdz
        def22(ix,iy,is)=-third2*zdwdz
        def33(ix,iy,is)=third4*zdwdz
      enddo
      enddo
      enddo
c
c boundary conditions d/ds=0
c
      do is=0,ns,ns
      do iy=1,ny1
      do ix=1,nx1
        zdwdz=-(s(ix,iy,is+1)+s(ix,iy,is))
     :    *(w(ix,iy,is+1,2)-w(ix,iy,is,2))/ds12(is)
        def11(ix,iy,is)=-third2*zdwdz
        def22(ix,iy,is)=-third2*zdwdz
        def33(ix,iy,is)=third4*zdwdz
      enddo
      enddo
      enddo
c
      do is=1,ns-1
      do iy=1,ny
      do ix=1,nx
        def12(ix,iy,is)=0.
      enddo
      enddo
      enddo
c
      do is=0,ns,ns
      do iy=1,ny
      do ix=1,nx
        def12(ix,iy,is)=0.
c     :    (u(ix,iy+1,is,2)-u(ix,iy,is,2))/dy
c     :    +(v(ix+1,iy,is,2)-v(ix,iy,is,2))/dx
      enddo
      enddo
      enddo
c
      do is=1,ns
      do iy=1,ny1
      do ix=1,nx
        def13(ix,iy,is)=-(s(ix+1,iy,is)+s(ix,iy,is))
     :    *(u(ix,iy,is,2)-u(ix,iy,is-1,2))/ds02(is)
      enddo
      enddo
      enddo
c
      do is=1,ns
      do iy=1,ny
      do ix=1,nx1
        def23(ix,iy,is)=-(s(ix,iy,is)+s(ix,iy+1,is))
     :    *(v(ix,iy,is,2)-v(ix,iy,is-1,2))/ds02(is)
      enddo
      enddo
      enddo
c
      if(defout.and.prt) then
        tk=0.
        call wri3ar(def11,1,nx1,1,ny1,1,ns-1,'def11   ',tk,iomod,nx,ny)
        call wri3ar(def12,1,nx,1,ny,1,ns-1,'def12   ',tk,iomod,nx,ny)
        call wri3ar(def13,1,nx,1,ny1,1,ns-1,'def13   ',tk,iomod,nx,ny)
        call wri3ar(def22,1,nx1,1,ny1,1,ns-1,'def22   ',tk,iomod,nx,ny)
        call wri3ar(def23,1,nx1,1,ny,1,ns-1,'def23   ',tk,iomod,nx,ny)
        call wri3ar(def33,1,nx1,1,ny1,1,ns-1,'def33   ',tk,iomod,nx,ny)
      endif
c
c calculation of diffusion coefficients:
c
      if(iodif.eq.1) then
c-
c calculate linear diffusion coefficients.
c-
        do is=0,ns
        do iy=1,ny1
        do ix=1,nx1
          dif000(ix,iy,is)=difcof
        enddo
        enddo
        enddo
      elseif(iodif.eq.2) then
c-
c calculate deformation/richardson number dependent diffusion coefficient
c-
        if(prt.and.driout) allocate(ri(0:nx1,0:ny1,0:ns1))
	  do is=1,ns-1
        do iy=2,ny
        do ix=2,nx
c
          def=dsqrt(0.5*(def11(ix,iy,is)*def11(ix,iy,is)+def22(ix,iy,is)
     :      *def22(ix,iy,is)+def33(ix,iy,is)*def33(ix,iy,is))
     :      +0.25*(def12(ix,iy,is)*def12(ix,iy,is)
     :      +def12(ix-1,iy,is)*def12(ix-1,iy,is)
     :      +def12(ix,iy-1,is)*def12(ix,iy-1,is)
     :      +def12(ix-1,iy-1,is)*def12(ix-1,iy-1,is)
     :      +def13(ix,iy,is)*def13(ix,iy,is)
     :      +def13(ix-1,iy,is)*def13(ix-1,iy,is)
     :      +def13(ix,iy,is+1)*def13(ix,iy,is+1)
     :      +def13(ix-1,iy,is+1)*def13(ix-1,iy,is+1)
     :      +def23(ix,iy,is)*def23(ix,iy,is)
     :      +def23(ix,iy-1,is)*def23(ix,iy-1,is)
     :      +def23(ix,iy,is+1)*def23(ix,iy,is+1)
     :      +def23(ix,iy-1,is+1)*def23(ix,iy-1,is+1)))
	    
c
          if(qif.ne.0. .and. ifqc.ne.0) then
            p=sigma0(is)*pp(ix,iy,2)+ptop
            p2=sigma0(is+1)*pp(ix,iy,2)+ptop
            p1=sigma0(is-1)*pp(ix,iy,2)+ptop
            t=(pts(ix,iy,is)+pt(ix,iy,is,2))*(p/p00)**akapa

           tp2=(pts(ix,iy,is+1)+pt(ix,iy,is+1,2))*(p2/p00)**akapa
            tp1=(pts(ix,iy,is-1)+pt(ix,iy,is-1,2))*(p1/p00)**akapa
            qsatur=qsat(t,p)
            qsat2=qsat(tp2,p2)
            qsat1=qsat(tp1,p1)
            xnsq=-g*(s(ix,iy,is+1)+s(ix,iy,is))*
     :      ((1+hlat*qsatur/(rv*t))/(1+0.622*hlat*hlatcp*qsatur/
     :      (rv*t*t))
     :      *((pt(ix,iy,is+1,2)
     :      +pts(ix,iy,is+1)-pt(ix,iy,is-1,2)-pts(ix,iy,is-1))
     :      /(pt(ix,iy,is,2)+pts(ix,iy,is))+hlatcp/t*(qsat2-qsat1))-
     :      (qv(ix,iy,is+1,2)+qc(ix,iy,is+1,2)-qv(ix,iy,is-1,2)-
     :      qc(ix,iy,is-1,2)))
     :      /ds04a(is)
	    else
            xnsq=-g*(s(ix,iy,is+1)+s(ix,iy,is))*(pt(ix,iy,is+1,2)
     :        +pts(ix,iy,is+1)-pt(ix,iy,is-1,2)-pts(ix,iy,is-1))
     :        /(pt(ix,iy,is,2)+pts(ix,iy,is))
     :        /ds04a(is)
          endif
          rich=xnsq/(def*def+zero0)
          if(prt.and.driout) ri(ix,iy,is)=rich
c
          difk=delka2(is)*def*dsqrt(dim(1.0d0,rich*rikey))
          dif000(ix,iy,is)=difk+difl
        enddo
        enddo
        enddo
c
        call extra3(nx1,ny1,ns1,dif000,1,nx1,1,ny1,0,ns)
        do iy=2,ny
        do ix=2,nx
          dif000(ix,iy,ns)=dif000(ix,iy,ns-1)
        enddo
        enddo
c
        if(prt.and.driout) then
          if(grdout) then
            call wgrids(ri,'ri',2,nx,2,ny,2,ns-1,iomod
     :        ,0,0,0,fngrd)
            call wgrids(dif000,'kk',2,nx,2,ny,2,ns-1,iomod
     :        ,0,0,0,fngrd)
	    else
            tk=0.
            call wri3ar(ri,2,nx,2,ny,2,ns-1,'ri      ',tk,iomod,nx,ny)
            call wri3ar(dif000,2,nx,2,ny,2,ns-1,'dif000  ',tk,iomod
     :        ,nx,ny)
          endif
          deallocate(ri)
        endif

      endif
c
c top boundary enhanced diffusion:
c
      if(olddra) then
        ixm=nx1/2
        iym=ny1/2
        ztop=phis(ixm,iym,1)/g
        do is=0,ns
        do iy=1,ny1
        do ix=1,nx1
          z=(phi(ix,iy,is)+phis(ix,iy,is))/g
          drag=(dragm+(dragt-dragm)*cos(pi05*abs((z-ztop)/(zm-ztop)))
     :      **2)*dim(z,zm)/(z-zm+zero0)
          dif000(ix,iy,is)=dif000(ix,iy,is)+drag
        enddo
        enddo
        enddo
      elseif(.not.raylei) then
         do is=0,idrmax
         do iy=1,ny1
         do ix=1,nx1
           dif000(ix,iy,is)=dif000(ix,iy,is)+endrag(is)
         enddo
         enddo
         enddo
      endif
c
c
c impose linear stability for diffusion coefficient
c
c      do is=0,ns
c      do iy=1,ny1
c      do ix=1,nx1
c        dif000(ix,iy,is)=min(deltaz(ix,iy)**2/(8.*dt),dif000(ix,iy,is))
c      enddo
c      enddo
c      enddo

      if(prt.and.ddiout) then
         tk=0.
         call wri3ar(dif000,1,nx1,1,ny1,0,ns,'difcof  ',tk,iomod,nx,ny)
      endif
c
c the upper layer drag is non-linear. so maybe it should be calculated
c separatelly for the other grid.
c
c diffusion coefficients on other grid locations:
c

c      allocate(dif010(0:nx1,0:ny1,0:ns1))
c      allocate(dif100(0:nx1,0:ny1,0:ns1))

      do is=0,ns
      do iy=1,ny
      do ix=1,nx
        dif010(ix,iy,is)=0.5*(dif000(ix,iy,is)+dif000(ix,iy+1,is))
        dif100(ix,iy,is)=0.5*(dif000(ix,iy,is)+dif000(ix+1,iy,is))
      enddo
      enddo
      enddo
c
      call extrah(nx1,ny1,dif010,0,nx1,0,ny1,0,ns)
      call extrah(nx1,ny1,dif100,0,nx1,0,ny1,0,ns)
c
c diffuse only the perturbated field:
c
      if(uvdif.ne.1.0) then
c
         do is=0,ns
         do iy=1,ny1
         do ix=1,nx1
           def33(ix,iy,is)=def33(ix,iy,is)*dif000(ix,iy,is)
         enddo
         enddo
         enddo
c
         do is=1,ns
         do iy=1,ny1
         do ix=1,nx
           def13s=-(s(ix+1,iy,is)+s(ix,iy,is))
     :       *(us(ix,iy,is)-us(ix,iy,is-1))/ds02(is)
           def13(ix,iy,is)=(def13(ix,iy,is)-(1.0-uvdif)*def13s)
     :       *0.5*(dif100(ix,iy,is)+dif100(ix,iy,is-1))
         enddo
         enddo
         enddo
c
         do is=1,ns
         do iy=1,ny
         do ix=1,nx1
           def23s=-(s(ix,iy+1,is)+s(ix,iy,is))
     :       *(vs(ix,iy,is)-vs(ix,iy,is-1))/ds02(is)
           def23(ix,iy,is)=(def23(ix,iy,is)-(1.0-uvdif)*def23s)
     :       *0.5*(dif010(ix,iy,is)+dif010(ix,iy,is-1))
         enddo
         enddo
         enddo
c
      else
         do is=0,ns
         do iy=1,ny1
         do ix=1,nx1
           def33(ix,iy,is)=def33(ix,iy,is)*dif000(ix,iy,is)
         enddo
         enddo
         enddo
c
         do is=1,ns
         do iy=1,ny1
         do ix=1,nx
           def13(ix,iy,is)=def13(ix,iy,is)*0.5*(dif100(ix,iy,is)
     :       +dif100(ix,iy,is-1))
         enddo
         enddo
         enddo
c
         do is=1,ns
         do iy=1,ny
         do ix=1,nx1
           def23(ix,iy,is)=def23(ix,iy,is)*0.5*(dif010(ix,iy,is)
     :       +dif010(ix,iy,is-1))
         enddo
         enddo
         enddo
c
      endif
c--
c difunu - teixeira
c--
      do iy=2,ny
      do ix=1,nx
        do is=1,ns
          deltazt(is)=-ds02(is+1)/(s(ix,iy,is+1)+s(ix+1,iy,is+1)) !(estava -)
c          write(*,*) 'deltazt - difunu', deltazt(is)
        enddo

        if(verbose.ge.2) write(*,*) 'Call teixeira/difunu'
        call teixeira(0,nx+1,0,ny+1,0,ns+1,ix,iy,2,ns-1
     :    ,u(0,0,0,2),difunu,dif100,1,w1d1,w1d2,deltazt,dtl) !2->3

        do is=1,ns
          difunu(ix,iy,is)=difunu(ix,iy,is)*pp10(ix,iy,2)
        enddo
      enddo
      enddo
c
      do iy=2,ny
      do ix=1,nx
        difunu(ix,iy,0)=difunu(ix,iy,1)
        difunu(ix,iy,ns)=difunu(ix,iy,ns-1)
      enddo
      enddo

      if(prt.and.ddiout) then
         tk=0.
         call wri3ar(difunu,1,nx1,1,ny1,0,ns,'difunu',tk,iomod,nx,ny) !estava difunu  '
      endif
c
c--
c difunv - teixeira
c--
      do iy=1,ny
      do ix=2,nx
        do is=1,ns
          deltazt(is)=-ds02(is+1)/(s(ix,iy,is+1)+s(ix,iy+1,is+1))!estava -
c          write(*,*) 'deltazt -difunv',deltazt(is)
        enddo

        if(verbose.ge.2) write(*,*) 'Call teixeira/difunv'
        call teixeira(0,nx+1,0,ny+1,0,ns+1,ix,iy,2,ns-1
     :    ,v(0,0,0,2),difunv,dif010,1,w1d1,w1d2,deltazt,dtl) !2->3

        do is=1,ns
          difunv(ix,iy,is)=difunv(ix,iy,is)*pp01(ix,iy,2)
        enddo
      enddo
      enddo
c
      do iy=1,ny
      do ix=2,nx
        difunv(ix,iy,0)=difunv(ix,iy,1)
        difunv(ix,iy,ns)=difunv(ix,iy,ns-1)
      enddo
      enddo
c--
c  difunw - teixeira
c--
      do is=1,ns
      do iy=1,ny
      do ix=1,nx
        dif001(ix,iy,is)=0.5d0*(dif000(ix,iy,is)+dif000(ix,iy,is-1))
      enddo
      enddo
      enddo

      do iy=2,ny
      do ix=2,nx
        do is=2,ns-1
          deltazt(is)=-ds12(is)/(s(ix,iy,is)+s(ix,iy,is+1)) !estava -
c          write(*,*) 'deltazt -difunw',deltazt(is)
        enddo
        if(verbose.ge.2) write(*,*) 'Call teixeira/difunw'
        call teixeira(0,nx+1,0,ny+1,0,ns+1,ix,iy,2,ns-1
     :    ,w(0,0,0,2),difunw,dif001,1,w1d1,w1d2,deltazt,dtl) !2->3

        do is=2,ns-1
          difunw(ix,iy,is)=difunw(ix,iy,is)*pp(ix,iy,2)
c          write(*,*) 'difunw(x,y,is),pp',difunw(ix,iy,is),pp(ix,iy,2),is
        enddo
      enddo
      enddo
c
      do iy=2,ny
      do ix=2,nx
        difunw(ix,iy,ns)=difunw(ix,iy,ns-1)
        difunw(ix,iy,1)=difunw(ix,iy,2)
      enddo
      enddo
c-
c difunt - teixeira         _>Acho que o problema com a estabilidade est?aqui
c-
      do is=1,ns
      do iy=1,ny
      do ix=1,nx
        h3(ix,iy,is)=pt(ix,iy,is,2)+tsdif*pts(ix,iy,is)
      enddo
      enddo
      enddo

      do iy=2,ny
      do ix=2,nx
        do is=2,ns-1
          deltazt(is)=-ds0(is+1)/s(ix,iy,is+1) !estava -
c          write(*,*) 'deltazt -difunt',deltazt(is)
        enddo
        if(verbose.ge.2) write(*,*) 'Call teixeira/difunt'
        call teixeira(0,nx+1,0,ny+1,0,ns+1,ix,iy,2,ns-1
     :    ,h3,difunt,dif000,1,w1d1,w1d2,deltazt,dtl)

c        do is=2,ns-1                                    !Esta express? n? deve estar correcta,
c          difunt(ix,iy,is)=difunt(ix,iy,is)/20000.d0    !foi acrescentada por mim para evitar a instabilidade.
c        enddo                                           !A quest? ? qual a real defini?o de difunt? Cf. com difunt
                                                         !do m?odo cl?sico, na qual entra h3(ix,iy,is)

      enddo
      enddo
c
      if(qif.gt.0.) then
c-
c difunq - teixeira
c-
        do is=1,ns
        do iy=1,ny
        do ix=1,nx
          h3(ix,iy,is)=qv(ix,iy,is,2)+qsdif*qvs(ix,iy,is)
        enddo
        enddo
        enddo

        do iy=1,ny
        do ix=1,nx
          do is=1,ns
            deltazt(is)=-ds0(is+1)/s(ix,iy,is+1)
          enddo
          if(verbose.ge.2) write(*,*) 'Call teixeira/difunq'
          call teixeira(0,nx+1,0,ny+1,0,ns+1,ix,iy,2,ns-1
     :      ,h3,difunq,dif000,1,w1d1,w1d2,deltazt,dtl)
        enddo
        enddo
c
      endif
c-
c difunqc - teixeira
c-
      if(ifqc.gt.0.) then

c
        do iy=1,ny
        do ix=1,nx
          do is=1,ns
            deltazt(is)=-ds0(is+1)/s(ix,iy,is+1)
          enddo
          if(verbose.ge.0) write(*,*) 'Call teixeira/difunqc'
          call teixeira(0,nx+1,0,ny+1,0,ns+1,ix,iy,2,ns-1
     :      ,qc,difunqc,dif000,1,w1d1,w1d2,deltazt,dtl)
        enddo
        enddo
c
      endif

c-
c difunqci - teixeira
c-
      if(ifqi.gt.0.) then                                                       DC,11.2009

c
        do iy=1,ny                                                              DC,11.2009
        do ix=1,nx                                                              DC,11.2009
          do is=1,ns                                                            DC,11.2009
            deltazt(is)=-ds0(is+1)/s(ix,iy,is+1)                                DC,11.2009
          enddo                                                                 DC,11.2009
          if(verbose.ge.0) write(*,*) 'Call teixeira/difunqci'                  DC,11.2009
          call teixeira(0,nx+1,0,ny+1,0,ns+1,ix,iy,2,ns-1                       DC,11.2009
     :      ,qci,difunqci,dif000,1,w1d1,w1d2,deltazt,dtl)                       DC,11.2009
        enddo                                                                   DC,11.2009
        enddo                                                                   DC,11.2009
c
      endif                                                                     DC,11.2009

c-
c difunqsn - teixeira
c-
      if(ifqi.gt.0.) then                                                       DC,11.2009

c
        do iy=1,ny                                                              DC,11.2009
        do ix=1,nx                                                              DC,11.2009
          do is=1,ns                                                            DC,11.2009
            deltazt(is)=-ds0(is+1)/s(ix,iy,is+1)                                DC,11.2009
          enddo                                                                 DC,11.2009
          if(verbose.ge.0) write(*,*) 'Call teixeira/difunqci'                  DC,11.2009
          call teixeira(0,nx+1,0,ny+1,0,ns+1,ix,iy,2,ns-1                       DC,11.2009
     :      ,qsn,difunqsn,dif000,1,w1d1,w1d2,deltazt,dtl)                       DC,11.2009
        enddo                                                                   DC,11.2009
        enddo                                                                   DC,11.2009
c
      endif                                                                     DC,11.2009
c-
c difunqr - teixeira
c-
      if(ifqr.gt.0.) then
        do iy=1,ny
        do ix=1,nx
          do is=1,ns
            deltazt(is)=-ds0(is+1)/s(ix,iy,is+1)
          enddo
          if(verbose.ge.0) write(*,*) 'Call teixeira/difunqr'
          call teixeira(0,nx+1,0,ny+1,0,ns+1,ix,iy,2,ns-1
     :      ,qr,difunqr,dif000,1,w1d1,w1d2,deltazt,dtl)
        enddo
        enddo
      endif

c
c debug:
c
      if(iodify.eq.0) then
        call extray(nx1,ny1,difunu,1,nx,2,ny,0,ns)
        call extray(nx1,ny1,difunv,2,nx,2,ny-1,0,ns)
        call extray(nx1,ny1,difunw,1,nx1,2,ny,0,ns)
      endif
c
c these extrapolations get rid of unused boundary values for diffusion:
c
      call extrah(nx1,ny1,difunu,1,nx,1,ny1,0,ns)
      call extrah(nx1,ny1,difunv,1,nx1,1,ny,0,ns)
      call extrah(nx1,ny1,difunw,1,nx1,1,ny1,1,ns)
c
      return
      end

	subroutine diffu_local(fngrd
     :  ,def11,h1,q1
     :  ,def12,h2,q2
     :  ,def13,h3,q3
     :  ,def22
     :  ,def23
     :  ,def33
     :  ,dif000
     :  ,dif010
     :  ,dif001,dif100)

      use alloc
	implicit none
	character*80 fngrd
      real(kind(0.d0)),dimension(0:nx1,0:ny1,0:ns1)::
     :  def11,def12,def13,def22,def23,def33
     :  ,h1,h2,h3,q1,q2,q3,dif000,dif010,dif001,dif100,dift000
	real(kind(0.d0)),dimension(0:nx1,0:ny1):: lmax
	real*8 z,mixl,def,rich,difk,xnsq,qsat1,qsat2,
     :  qsatur,tp1,tp2,t,p1,p2,p,drag,ztop,tk,dift,
     :  difk_h
	integer is,iy,ix,ixm,iym
	real,parameter:: karm=0.4
	real,external:: qsat
!-----------------------------------------------------------------------------------------------------!
!      write(0,*) 'BROCKEN SOCIAL SCENE 2010'
      do ix=1,nx1
	do iy=1,ny1
      if (ust_s(ix,iy)*tst_s(ix,iy).gt.0) then
!	  lmax(ix,iy)=max(50.,0.007*ust_s(ix,iy)/fcor)
!        lmax(ix,iy)=max(50.,0.15*hbl(ix,iy))
        lmax(ix,iy)=40.
	else
	lmax(ix,iy)= 200. !max(50.,0.15*hbl(ix,iy))
	endif
!      lmax(ix,iy)=150.
	enddo
	enddo
	

	do ix=1,nx1
	w(ix,ny1,ns+1,2)=w(ix,ny1-1,ns+1,2)
	enddo

	do iy=1,ny1
	w(nx1,iy,ns+1,2)=w(nx1-1,iy,ns+1,2)
	enddo
      
	
      do is=1,ns-1
      do iy=1,ny1
      do ix=1,nx1
        def33(ix,iy,is)=-(s(ix,iy,is+1)+s(ix,iy,is))
     :    *(w(ix,iy,is+1,2)-w(ix,iy,is,2))/ds12(is)
	enddo
      enddo
      enddo
c
c boundary conditions d/ds=0
c
      do is=0,ns,ns
      do iy=1,ny1
      do ix=1,nx1
        def33(ix,iy,is)=-(s(ix,iy,is+1)+s(ix,iy,is))
     :    *(w(ix,iy,is+1,2)-w(ix,iy,is,2))/ds12(is)
      enddo
      enddo
      enddo

!	do is=1,ns
!     do iy=1,ny1
!      do ix=1,nx
!        def13(ix,iy,is)=-(s(ix+1,iy,is)+s(ix,iy,is))
!     :    *(u(ix,iy,is,2)-u(ix,iy,is-1,2))/ds02(is)
!      enddo
!      enddo
!      enddo
c
!      do is=1,ns
!      do iy=1,ny
!      do ix=1,nx1
!        def23(ix,iy,is)=-(s(ix,iy,is)+s(ix,iy+1,is))
!     :    *(v(ix,iy,is,2)-v(ix,iy,is-1,2))/ds02(is)
!      enddo
!      enddo
!      enddo

	do is=1,ns
      do iy=1,ny1
      do ix=1,nx
        def13(ix,iy,is)=-(s(ix+1,iy,is)+s(ix,iy,is))
     :    *(u(ix,iy,is,2)-u(ix,iy,is-1,2))/ds02(is)
     :    +(w(ix+1,iy,is,2)-w(ix,iy,is,2))/dx-s1ds4a(is)*ppdx10(ix,iy,2)
     :    *(w(ix+1,iy,is+1,2)+w(ix,iy,is+1,2)-w(ix+1,iy,is-1,2)
     :    -w(ix,iy,is-1,2))
	 enddo
      enddo
      enddo
      do is=1,ns
      do iy=1,ny
      do ix=1,nx1
        def23(ix,iy,is)=-(s(ix,iy,is)+s(ix,iy+1,is))
     :    *(v(ix,iy,is,2)-v(ix,iy,is-1,2))/ds02(is)
     :    +(w(ix,iy+1,is,2)-w(ix,iy,is,2))/dy-s1ds4a(is)*ppdy01(ix,iy,2)
     :    *(w(ix,iy+1,is+1,2)+w(ix,iy,is+1,2)-w(ix,iy+1,is-1,2)
     :    -w(ix,iy,is-1,2))
      enddo
      enddo
      enddo

	if(iodif.eq.1) then
c-
c calculate linear diffusion coefficients.
c-
        do is=0,ns
        do iy=1,ny1
        do ix=1,nx1
          dif000(ix,iy,is)=difcof
        enddo
        enddo
        enddo
      else
c-
c calculate deformation/richardson number dependent diffusion coefficient
c-
        if(prt.and.driout) allocate(ri(0:nx1,0:ny1,0:ns1))
	if(prt.and.driout) allocate(bosa(0:nx1,0:ny1,0:ns1))
	
	do is=1,ns-1
        do iy=2,ny
        do ix=2,nx
c
          def=dsqrt((0.25*(s(ix,iy,is+1)+s(ix,iy,is))*((u(ix,iy,is-1,2)
     :	-u(ix,iy,is+1,2))/ds02(is)+(u(ix-1,iy,is-1,2)
     :    -u(ix-1,iy,is+1,2))
     :    /ds02(is)))**2+(0.25*(s(ix,iy,is+1)+s(ix,iy,is))*
     :    ((v(ix,iy,is-1,2)-v(ix,iy,is+1,2))/ds02(is)+
     :    (v(ix,iy-1,is-1,2)-v(ix,iy-1,is+1,2))/ds02(is)))**2)

	 !
	 ! if(is.eq.35) write(0,*) is,def
	    if(prt.and.driout) bosa(ix,iy,is)=def
c
    !       if(qif.ne.0. .and. ifqc.ne.0) then
              
           if(qif.ne.0.and.ifqc.ge.0) then
             p=sigma0(is)*pp(ix,iy,2)+ptop
              t=(pts(ix,iy,is)+pt(ix,iy,is,2))*(p/p00)**akapa
             if(qv(ix,iy,is,2).ge.0.9*qsat(t,p)) then
            
 !           p2=sigma0(is+1)*pp(ix,iy,2)+ptop
 !           p1=sigma0(is-1)*pp(ix,iy,2)+ptop
            
 !           tp2=(pts(ix,iy,is+1)+pt(ix,iy,is+1,2))*(p2/p00)**akapa
 !           tp1=(pts(ix,iy,is-1)+pt(ix,iy,is-1,2))*(p1/p00)**akapa
 !           qsatur=qsat(t,p)
 !           qsat2=qsat(tp2,p2)
 !           qsat1=qsat(tp1,p1)
 !           xnsq=-g*(s(ix,iy,is+1)+s(ix,iy,is))*
 !    :      ((1+hlat*qsatur/(rv*t))/(1+0.622*hlat*hlatcp*qsatur/
 !    :      (rv*t*t))
 !    :      *((pt(ix,iy,is+1,2)
 !    :      +pts(ix,iy,is+1)-pt(ix,iy,is-1,2)-pts(ix,iy,is-1))
 !    :      /(pt(ix,iy,is,2)+pts(ix,iy,is))+hlatcp/t*(qsat2-qsat1))-
 !    :      (qv(ix,iy,is+1,2)+qc(ix,iy,is+1,2)-qv(ix,iy,is-1,2)-
 !    :      qc(ix,iy,is-1,2)))
 !    :      /ds04a(is)
               
            xnsq=-g*(s(ix,iy,is+1)+s(ix,iy,is))*(pt(ix,iy,is+1,2)
     :        +pts(ix,iy,is+1)+hlat/cp*qv(ix,iy,is+1,2)
     :        -pt(ix,iy,is-1,2)-pts(ix,iy,is-1)
     :        -hlat/cp*qv(ix,iy,is-1,2))
     :        /(pt(ix,iy,is,2)+pts(ix,iy,is)+
     :          +hlat/cp*qv(ix,iy,is,2))
     :        /ds04a(is)
               
          else
            xnsq=-g*(s(ix,iy,is+1)+s(ix,iy,is))*(pt(ix,iy,is+1,2)
     :        +pts(ix,iy,is+1)-pt(ix,iy,is-1,2)-pts(ix,iy,is-1))
     :        /(pt(ix,iy,is,2)+pts(ix,iy,is))
     :        /ds04a(is)
          endif
           else
             xnsq=-g*(s(ix,iy,is+1)+s(ix,iy,is))*(pt(ix,iy,is+1,2)
     :        +pts(ix,iy,is+1)-pt(ix,iy,is-1,2)-pts(ix,iy,is-1))
     :        /(pt(ix,iy,is,2)+pts(ix,iy,is))
     :        /ds04a(is)
          endif
          
          rich=xnsq/(def*def+zero0)
          if(prt.and.driout) ri(ix,iy,is)=rich
c

          
		  z=(phis(ix,iy,is)+phi(ix,iy,is))/g
          mixl=karm*z/(1.+karm*z/lmax(ix,iy))
	    if(rich.ge.0) then
          difk=mixl**2.*def*(max(0.,(1.-5.*rich))**2.)
	    dift=difk
	    endif
	    if(rich.lt.0) then
	    difk=mixl**2.*def*(min(9.,sqrt(1.-16.*rich)))
	    dift=difk*(min(3.,(1.-16.*rich)**0.25))
	 
	     !difk=delka2(is)*def*dsqrt(dim(1.0d0,rich*rikey))
	     !dift=difk
	
          endif
  !        z=(phi(ix,iy,is)+phis(ix,iy,is))/g
  !         if(z.gt.hbl(ix,iy)) then
  !         difk=0
  !         dift=0
  !         endif
 	   dif000(ix,iy,is)=difk+difl
	    dift000(ix,iy,is)=dift
        enddo
        enddo
        enddo
	call extra3(nx1,ny1,ns1,dif000,1,nx1,1,ny1,0,ns)
	call extra3(nx1,ny1,ns1,dift000,1,nx1,1,ny1,0,ns)
c
        if(prt.and.driout) then
          if(grdout) then
            call wgrids(ri,'ri',2,nx,2,ny,2,ns-1,iomod
     :        ,0,0,0,fngrd)
            call wgrids(dif000,'kk',2,nx,2,ny,2,ns-1,iomod
     :        ,0,0,0,fngrd)
	      call wgrids(bosa,'de',2,nx,2,ny,2,ns-1,iomod
     :        ,0,0,0,fngrd)
	  else
            tk=0.
            call wri3ar(ri,2,nx,2,ny,2,ns-1,'ri      ',tk,iomod,nx,ny)
            call wri3ar(dif000,2,nx,2,ny,2,ns-1,'dif000  ',tk,iomod
     :        ,nx,ny)
          endif
          deallocate(ri)
	    deallocate(bosa)
        endif
        endif
c
c top boundary enhanced diffusion:
c
      if(olddra) then
        ixm=nx1/2
        iym=ny1/2
        ztop=phis(ixm,iym,1)/g
        do is=0,ns
        do iy=1,ny1
        do ix=1,nx1
          z=(phi(ix,iy,is)+phis(ix,iy,is))/g
          drag=(dragm+(dragt-dragm)*cos(pi05*abs((z-ztop)/(zm-ztop)))
     :      **2)*dim(z,zm)/(z-zm+zero0)
          dif000(ix,iy,is)=dif000(ix,iy,is)+drag
	    dift000(ix,iy,is)=dift000(ix,iy,is)+drag
	    
        enddo
        enddo
        enddo
      elseif(.not.raylei) then
        do is=0,idrmax
        do iy=1,ny1
        do ix=1,nx1
          dif000(ix,iy,is)=dif000(ix,iy,is)+endrag(is)
	    dift000(ix,iy,is)=dift000(ix,iy,is)+endrag(is)
        enddo
        enddo
        enddo
      endif
       
	do is=0,ns
      do iy=1,ny
      do ix=1,nx
        dif010(ix,iy,is)=0.5*(dif000(ix,iy,is)+dif000(ix,iy+1,is))
        dif100(ix,iy,is)=0.5*(dif000(ix,iy,is)+dif000(ix+1,iy,is))
      enddo
      enddo
      enddo
c
      call extrah(nx1,ny1,dif010,0,nx1,0,ny1,0,ns)
      call extrah(nx1,ny1,dif100,0,nx1,0,ny1,0,ns)
c
      do is=0,ns
         do iy=1,ny1
         do ix=1,nx1
           def33(ix,iy,is)=def33(ix,iy,is)*dif000(ix,iy,is)
         enddo
         enddo
         enddo
c
c
         do is=1,ns
         do iy=1,ny1
         do ix=1,nx
           def13(ix,iy,is)=def13(ix,iy,is)*0.5*(dif100(ix,iy,is)
     :       +dif100(ix,iy,is-1))
         enddo
         enddo
         enddo
c
         do is=1,ns
         do iy=1,ny
         do ix=1,nx1
           def23(ix,iy,is)=def23(ix,iy,is)*0.5*(dif010(ix,iy,is)
     :       +dif010(ix,iy,is-1))
         enddo
         enddo
         enddo

	do iy=1,ny1
	do ix=1,nx
	def13(ix,iy,ns)=0.5*(cdm(ix,iy)+cdm(ix+1,iy))*u(ix,iy,ns-1,2)
	enddo
	enddo
	do iy=1,ny
         do ix=1,nx1
      def23(ix,iy,ns)=0.5*(cdm(ix,iy)+cdm(ix,iy+1))*v(ix,iy,ns-1,2)
	enddo
	enddo
c--
c difunu
c--
      do is=1,ns-1
      do iy=2,ny
      do ix=1,nx
        difunu(ix,iy,is)=(
     :    -(s(ix+1,iy,is+1)+s(ix+1,iy,is)+s(ix,iy,is+1)+s(ix,iy,is))
     :    *(def13(ix,iy,is+1)-def13(ix,iy,is))/ds14(is))*pp10(ix,iy,2)
      enddo
      enddo
      enddo
c
      do iy=2,ny
      do ix=1,nx
        difunu(ix,iy,0)=difunu(ix,iy,1)
        difunu(ix,iy,ns)=difunu(ix,iy,ns-1)
      enddo
      enddo

c--
c difunv
c--
      do is=1,ns-1
      do iy=1,ny
      do ix=2,nx
        difunv(ix,iy,is)=(
     :    -(s(ix,iy+1,is+1)+s(ix,iy+1,is)+s(ix,iy,is+1)+s(ix,iy,is))
     :    *(def23(ix,iy,is+1)-def23(ix,iy,is))/ds14(is))*pp01(ix,iy,2)
      enddo
      enddo
      enddo
c
      do iy=1,ny
      do ix=2,nx
        difunv(ix,iy,0)=difunv(ix,iy,1)
        difunv(ix,iy,ns)=difunv(ix,iy,ns-1)
      enddo
      enddo


c--
c  difunw
c--
      do is=2,ns-1
      do iy=2,ny
      do ix=2,nx
        difunw(ix,iy,is)=(
     :    -s(ix,iy,is)*(def33(ix,iy,is)-def33(ix,iy,is-1))/ds0(is))
     :    *pp(ix,iy,2)
      enddo
      enddo
      enddo
c
      do iy=2,ny
      do ix=2,nx
        difunw(ix,iy,ns)=difunw(ix,iy,ns-1)
        difunw(ix,iy,1)=difunw(ix,iy,2)
      enddo
      enddo
      
 !     do is=1,ns-1
 !     do iy=1,ny1
 !     do ix=1,nx-1
 !       umfl(ix,iy,is)=0.25*(def13(ix,iy,is)+def13(ix+1,iy,is)+
 !    :  def13(ix,iy,is+1)+def13(ix+1,iy,is+1))
 !     enddo
 !     enddo
 !     enddo
c
 !     do is=1,ns-1
 !     do iy=1,ny-1
 !     do ix=1,nx1
 !       vmfl(ix,iy,is)=0.25*(def23(ix,iy,is)+def23(ix,iy+1,is)+
 !    :  def23(ix,iy,is+1)+def23(ix,iy+1,is+1))
 !     enddo
 !     enddo
 !     enddo
!--------------------------------------------------------------------------------!
!                            DIFUNT                                              !
!--------------------------------------------------------------------------------!
	if (qif.ne.0.and.ifqc.ge.0) then
	  do is=2,ns-1
        do iy=1,ny
        do ix=1,nx
        h3(ix,iy,is)=-s(ix,iy,is)*(pt(ix,iy,is,2)-pt(ix,iy,is-1,2)
     :    +hlat/cp*qv(ix,iy,is,2)
     :    -hlat/cp*qv(ix,iy,is-1,2)   
     :    +tsdif*(pts(ix,iy,is)-pts(ix,iy,is-1)))/ds0(is)
        h3(ix,iy,is)=h3(ix,iy,is)*0.5*(dift000(ix,iy,is)
     :    +dift000(ix,iy,is-1))*tdif	
        enddo
        enddo
        enddo
	else
	  do is=2,ns-1
        do iy=1,ny
        do ix=1,nx
        h3(ix,iy,is)=-s(ix,iy,is)*(pt(ix,iy,is,2)-pt(ix,iy,is-1,2)
     :    +tsdif*(pts(ix,iy,is)-pts(ix,iy,is-1)))/ds0(is)
        h3(ix,iy,is)=h3(ix,iy,is)*0.5*(dift000(ix,iy,is)
     :    +dift000(ix,iy,is-1))*tdif	   
        enddo
        enddo
        enddo
      endif
	
	!write(0,*)pts(18,18,47)+pt(18,18,47,2),u(18,18,47,2),v(18,18,47,2)
c
      do iy=1,ny1
      do ix=1,nx1
        h3(ix,iy,1)=0.0
        !h3(ix,iy,ns)=0.0
	  if(qif.ge.0.and.ifqc.ge.0) then
	    h3(ix,iy,ns)=(ust_s(ix,iy)*tst_s(ix,iy)+
     :    hlat/cp*qst_s(ix,iy)*ust_s(ix,iy))
!	    hfl(ix,iy,ns)=h(ix,iy)+le(ix,iy)
	  else
	    h3(ix,iy,ns)=ust_s(ix,iy)*tst_s(ix,iy)
!	    hfl(ix,iy,ns)=h(ix,iy)
	  endif
	  
      enddo
      enddo

	do is=1,ns-1
      do iy=2,ny
      do ix=2,nx
        difunt(ix,iy,is)=
     :    -(s(ix,iy,is+1)+s(ix,iy,is))*(h3(ix,iy,is+1)-h3(ix,iy,is))
     :    /ds12(is)
      enddo
      enddo
      enddo
      
      if(prt.and.driout) then
          if(grdout) then
            call wgrids(h3,'he',2,nx,2,ny,2,ns-1,iomod
     :        ,0,0,0,fngrd)
	    endif
	endif
	
      
!--------------------------------------------------------------------------------!
!                            DIFUNQ                                              !
!--------------------------------------------------------------------------------!
	if(qif.gt.0.) then

	do is=2,ns-1
      do iy=1,ny
      do ix=1,nx
        q3(ix,iy,is)=-s(ix,iy,is)*(qv(ix,iy,is,2)-qv(ix,iy,is-1,2)
     :    )/ds0(is)   !
        q3(ix,iy,is)=q3(ix,iy,is)*0.5*(dift000(ix,iy,is)
     :    +dift000(ix,iy,is-1))
      enddo
      enddo
      enddo
c
      do iy=1,ny1
      do ix=1,nx1
        q3(ix,iy,1)=0.0
        !q3(ix,iy,ns)=0.0
        q3(ix,iy,ns)=qst_s(ix,iy)*ust_s(ix,iy)
      enddo
      enddo
c
      do is=1,ns-1
      do iy=2,ny
      do ix=2,nx
        difunq(ix,iy,is)=
     :    -(s(ix,iy,is+1)+s(ix,iy,is))*(q3(ix,iy,is+1)-q3(ix,iy,is))
     :    /ds12(is)
      enddo
      enddo
      enddo
      endif
      
      if(prt.and.driout) then
          if(grdout) then
            call wgrids(q3,'qf',2,nx,2,ny,2,ns-1,iomod
     :        ,0,0,0,fngrd)
	    endif
	endif
!--------------------------------------------------------------------------------!
!                            DIFUNQC                                              !
!--------------------------------------------------------------------------------!
      if(ifqc.gt.0.) then
      do is=2,ns-1
      do iy=1,ny
      do ix=1,nx
        q3(ix,iy,is)=-s(ix,iy,is)*(qc(ix,iy,is,2)-qc(ix,iy,is-1,2))
     :     /ds0(is)
        q3(ix,iy,is)=q3(ix,iy,is)*0.5*(dift000(ix,iy,is)
     :     +dift000(ix,iy,is-1))
      enddo
      enddo
      enddo
c
      do iy=1,ny1
      do ix=1,nx1
        q3(ix,iy,1)=0.0
        q3(ix,iy,ns)=0.0
      enddo
      enddo
c
      do is=1,ns-1
      do iy=2,ny
      do ix=2,nx
        difunqc(ix,iy,is)=
     :    -(s(ix,iy,is+1)+s(ix,iy,is))*(q3(ix,iy,is+1)-q3(ix,iy,is))
     :    /ds12(is)
      enddo
      enddo
      enddo

      endif
!--------------------------------------------------------------------------------!
!                            DIFUNQCI                                              !
!--------------------------------------------------------------------------------!		
      if(ifqi.gt.0.) then
	do is=2,ns-1                                                              
      do iy=1,ny                                                                
      do ix=1,nx                                                                
        q3(ix,iy,is)=-s(ix,iy,is)*(qci(ix,iy,is,2)-qci(ix,iy,is-1,2))           
     :     /ds0(is)                                                            
        q3(ix,iy,is)=q3(ix,iy,is)*0.5*(dift000(ix,iy,is)                         
     :     +dift000(ix,iy,is-1))*qdif                                          
      enddo                                                                    
      enddo                                                                  
      enddo                                                                
c                                          
      do iy=1,ny1                                                        
      do ix=1,nx1                                                           
        q3(ix,iy,1)=0.0                                                     
        q3(ix,iy,ns)=0.0                                                   
      enddo                                                                  
      enddo                                                                   
c
      do is=1,ns-1                                                         
      do iy=2,ny                                                               
      do ix=2,nx                                                              
        difunqci(ix,iy,is)=                                                 
     :    -(s(ix,iy,is+1)+s(ix,iy,is))*(q3(ix,iy,is+1)-q3(ix,iy,is))           
     :    /ds12(is)                                                            
      enddo                                                                   
      enddo                                                                   
      enddo                                                                   

      endif         
!--------------------------------------------------------------------------------!
!                            DIFUNQSN                                              !
!--------------------------------------------------------------------------------!

      if(ifqi.gt.0.) then
	do is=2,ns-1                                                           
      do iy=1,ny                                                           
      do ix=1,nx                                                              
        q3(ix,iy,is)=-s(ix,iy,is)*(qsn(ix,iy,is,2)-qsn(ix,iy,is-1,2))         
     :     /ds0(is)                                                            
        q3(ix,iy,is)=q3(ix,iy,is)*0.5*(dift000(ix,iy,is)                      
     :     +dift000(ix,iy,is-1))*qdif                                           
      enddo                                                                  
      enddo                                                                 
      enddo                                                                   
c                                          
      do iy=1,ny1                                                             
      do ix=1,nx1                                                             
        q3(ix,iy,1)=0.0                                                       
        q3(ix,iy,ns)=0.0                                                    
      enddo                                                                   
      enddo                                                                    
c
      do is=1,ns-1                                                              
      do iy=2,ny                                                               
      do ix=2,nx                                                              
        difunqsn(ix,iy,is)=                                                    
     :    -(s(ix,iy,is+1)+s(ix,iy,is))*(q3(ix,iy,is+1)-q3(ix,iy,is))            
     :    /ds12(is)                                                             
      enddo                                                                     
      enddo                                                                  
      enddo                                                                    

      endif 
!--------------------------------------------------------------------------------!
!                            DIFUNQR                                              !
!--------------------------------------------------------------------------------!
	
	if(ifqr.gt.0.) then
	do is=2,ns-1
      do iy=1,ny
      do ix=1,nx
        q3(ix,iy,is)=-s(ix,iy,is)*(qr(ix,iy,is,2)-qr(ix,iy,is-1,2))
     :     /ds0(is)
        q3(ix,iy,is)=q3(ix,iy,is)*0.5*(dift000(ix,iy,is)
     :     +dift000(ix,iy,is-1))*qdif
      enddo
      enddo
      enddo
c
      do iy=1,ny1
      do ix=1,nx1
        q3(ix,iy,1)=0.0
        q3(ix,iy,ns)=0.0
      enddo
      enddo
c
      do is=1,ns-1
      do iy=2,ny
      do ix=2,nx
        difunqr(ix,iy,is)=
     :    -(s(ix,iy,is+1)+s(ix,iy,is))*(q3(ix,iy,is+1)-q3(ix,iy,is))
     :    /ds12(is)
      enddo
      enddo
      enddo	
      endif
	
	if(iodify.eq.0) then
        call extray(nx1,ny1,difunu,1,nx,2,ny,0,ns)
        call extray(nx1,ny1,difunv,2,nx,2,ny-1,0,ns)
        call extray(nx1,ny1,difunw,1,nx1,2,ny,0,ns)
      endif
c
c these extrapolations get rid of unused boundary values for diffusion:
c
      call extrah(nx1,ny1,difunu,1,nx,1,ny1,0,ns)
      call extrah(nx1,ny1,difunv,1,nx1,1,ny,0,ns)
      call extrah(nx1,ny1,difunw,1,nx1,1,ny1,1,ns)                                                                  

	return
	end
