	
	module alloc
c
c this module includes all public varianles
c
c grid size:
c
      integer::
     :  nx,ny,ns
     :  ,nx1,ny1,ns1,nx2,ny2
     :  ,nxm1,nym1,nall,nxy,nxs,nys,jy,js,nall3,nallm,nallmm
      integer::
     :  isqv,isqc,isqr,nx1qv,ny1qv,ns1qv
     :  ,nx1qc,ny1qc,ns1qc
     :  ,nx1qr,ny1qr,ns1qr
     :  ,nx1qi,ny1qi,ns1qi

c     Dima@ 02.05.07                                                           
c     variable for impurity switcher
      integer::                                                                 DM,05.2007
     :  impur                                                                   DM,05.2007
c     flags indicating the immediate and constant impurity
      logical::                                                                 DM,05.2007
     :   immediate_impur, constant_impur                                        DM,05.2007
c     coordinates of the impurity signal point
      integer::                                                                 DM,05.2007
     :   ix_impur, iy_impur                                                     DM,05.2007
c     impurity le
      real(kind(0.d0)),allocatable,dimension(:,:)::le_impur                     DM,05.2007
c     impurity value
      real(kind(0.0D0)) :: impur_value                                          DM,05.2007
c
c those lines can be changed
c
c      parameter(nx=49,ny=49,ns=40)
c      parameter(isqv=0,isqc=0,isqr=0)
c
c      parameter(nx1=nx+1,ny1=ny+1,ns1=ns+1)
c      parameter(nx2=nx+2,ny2=ny+2)
c      parameter(nxm1=nx-1,nym1=ny-1)
c      parameter(nall=(nx1+1)*(ny1+1)*(ns1+1),nxy=(nx1+1)*(ny1+1))
c      parameter(nxs=(nx1+1)*(ns1+1),nys=(ny1+1)*(ns1+1))
c      parameter(jy=nx1+1,js=nxy,nall3=3*nall)
c      parameter(nallm=nall-nxy,nallmm=nallm-nxy)

c      parameter(nx1qv=nx1*isqv,ny1qv=ny1*isqv,ns1qv=ns1)
c      parameter(nx1qc=nx1*isqc,ny1qc=ny1*isqc,ns1qc=ns1)
c      parameter(nx1qr=nx1*isqr,ny1qr=ny1*isqr,ns1qr=ns1)

c
c u,v,w : wind components
c wsig:   d/dt (sigma) component of velocity in sigma coordinates
c us,vs: reference state wind
c dppref: reference state dp/dt
c psufref: reference state psuf

c     real(kind(0.d0)),dimension(0:nx1,0:ny1,0:ns1,2:3)::
c    :  u,v,w,wsig,pt
      real(kind(0.d0)),allocatable,dimension(:,:,:,:)::
     :  u,v,w,wsig,pt
      logical nonloc_Noh

c     real(kind(0.0D0)),dimension(0:nx1qv,0:ny1qv,0:ns1qv,1:3)::qv
c     real(kind(0.0D0)),dimension(0:nx1qv,0:ny1qv,0:ns1qv)::
c    :  difunq,qvs,cond,evap

c      real(kind(0.0D0)),dimension(0:nx1qc,0:ny1qc,0:ns1qc,1:3)::qc
c      real(kind(0.0D0)),dimension(0:nx1qc,0:ny1qc,0:ns1qc)::
c     :  difunqc,auto,col

c      real(kind(0.0D0)),dimension(0:nx1qr,0:ny1qr,0:ns1qr,1:3)::qr
c      real(kind(0.0D0)),dimension(0:nx1qr,0:ny1qr,0:ns1qr)::
c     :  vrain,difunqr
c      real(kind(0.0D0)),dimension(0:nx1,0:ny1)::prec

      real(kind(0.0D0)),allocatable,dimension(:,:,:,:)::qv,qc,qr,qci,
     :  qsn                                                                     DC,10.2009
      real(kind(0.0D0)),allocatable,dimension(:,:,:)::
     :  difunq,qvs,cond,evap,difunqc,auto,col,vrain,difunqr,vdepi,              DC,10.2009
     :  vini,vdeps,imlt,hmfrz,sacrw,sacrwr,sbercw,sberci,ibercw,                DC,10.2009
     :  sagg,saci,raci,iacr,difunqci,smlt,sacrr,difunqsn,vsnow,qs,
     :  thgrad,adv_flx,adv_fly,tendx,tendy,ptnew,advfx,advfy                    DC,02,2012
      real(kind(0.0D0)),allocatable,dimension(:,:)::prec,precsnow               DC,11.2009

c     Dima@ 04.04.07
      real(kind(0.0D0)),allocatable,dimension(:,:,:,:)::qs_10                   DM,05.2007
      real(kind(0.0D0)),allocatable,dimension(:,:,:)::qs_10s                    DM,05.2007
      
c     Dima@ 13.04.07
      real(kind(0.0D0)),allocatable,dimension(:,:,:)::difunq_10                 DM,05.2007

c     real(kind(0.d0)),dimension(0:nx1,0:ny1,0:ns1)::
c    :  us,vs,pts,tems,difunu,difunv,difunw,difunt,phi,phis,s
c    :  ,wk1,wk2,wk3,wk4,wk5,wk6,wk7,wk8,wk9,thetav,difcoft

      real(kind(0.d0)),allocatable,dimension(:,:,:)::
     :  wk1,wk2,wk3,wk4,wk5,wk6,wk7,wk8,wk9

!     Dima@ 13.05.07     
      real(kind(0.d0)),allocatable,dimension(:,:,:)::                           DM,05.2007
     :  wk_impur                                                                DM,05.2007
      logical wk_impur_allocated                                                DM,05.2007 

      real(kind(0.d0)),allocatable,dimension(:,:,:)::
     :  us,vs,pts,tems,difunu,difunv,difunw,difunt,phi,phis,s
     :  ,thetav,difcoft,hfl,Ribulk,hflc,hfls,h3c,thpert,
     :  qflc,qfl,qfls,q3c,umfl,vmfl,mfl,ptfd,wfd,ufd,vfd,
     :  upert,vpert,rho,h3e
c     :  ,wk1,wk2,wk3,wk4,wk5,wk6,wk7,wk8,wk9
      real(kind(0.d0)),allocatable,dimension(:,:)::dppref,psufref
c      real(kind(0.d0)),dimension(0:nx1,0:ny1)::
c     :  tsuf,psufi,tsufi,phisui,deltaz,ptsoil,uvsuf,qvsoil,phisuf
c     :  ,phig,hk1,hk2,hk3,ppgra2,pplap,gpplap
      real(kind(0.d0)),allocatable,dimension(:,:)::
     :  tsuf,psufi,tsufi,phisui,deltaz,ptsoil,uvsuf,qvsoil,phisuf
     :  ,phig,hk1,hk2,hk3,ppgra2,pplap,gpplap,hbl
c      real(kind(0.d0))::
c     :  psuf(0:nx1,0:ny1)
c     :  ,pp(0:nx2,0:ny2,4),pp10(0:nx1,0:ny1,4),pp01(0:nx1,0:ny1,4)
c     :  ,ppdx10(0:nx1,0:ny1,2:2),ppdx(0:nx1,0:ny1,4)
c     :  ,dpdx10(0:nx1,0:ny1,3:3)
c     :  ,ppdy01(0:nx1,0:ny1,2:2),ppdy(0:nx1,0:ny1,4)
c     :  ,dpdy01(0:nx1,0:ny1,3:3)
c     :  ,hsuf(0:nx1,0:ny1),dpp(0:nx1,0:ny1)
c     :  ,hsufmax(0:nx1,0:ny1)
      real(kind(0.d0)),allocatable::
     :  psuf(:,:),pp(:,:,:),pp10(:,:,:),pp01(:,:,:),ppdx10(:,:,:)
     :  ,ppdx(:,:,:)
     :  ,dpdx10(:,:,:),ppdy01(:,:,:),ppdy(:,:,:),dpdy01(:,:,:)
     :  ,hsuf(:,:),dpp(:,:),hsufmax(:,:)
c
c boundary values of u,v,pp,pt - obtained by extrapolation (radiative
c boundary conditions)
c
      real(kind(0.d0)),allocatable::
     :   ubx(:,:,:),vbx(:,:,:)
     :  ,uby(:,:,:),vby(:,:,:)
     :  ,ucc(:,:,:),vcc(:,:,:)
     :  ,ubx2(:,:,:),vbx2(:,:,:)
     :  ,uby2(:,:,:),vby2(:,:,:)
     :  ,ucc2(:,:,:),vcc2(:,:,:)
     :  ,ppbx(:,:),ppby(:,:),ppcc(:,:)
     :  ,ppbx2(:,:),ppby2(:,:),ppcc2(:,:)
     :  ,ptbx(:,:,:),ptby(:,:,:),ptcc(:,:,:)
     :  ,qvbx(:,:,:),qvby(:,:,:),qvcc(:,:,:)
     :  ,qcbx(:,:,:),qcby(:,:,:),qccc(:,:,:)
     :  ,qrbx(:,:,:),qrby(:,:,:),qrcc(:,:,:)
     :  ,pgradx(:,:,:),pgrads(:,:,:)
     :  ,pgrady(:,:,:)
    
      real(kind(0.d0)),allocatable::                                            DM,05.2007
     & qscc(:,:,:),qsbx(:,:,:),qsby(:,:,:)                                      DM,05.2007 
	
	real(kind(0.d0)),allocatable::                                            DC,11.2009
     & qcicc(:,:,:),qcibx(:,:,:),qciby(:,:,:),                                  DC,11.2009
     & qsncc(:,:,:),qsnbx(:,:,:),qsnby(:,:,:)                                   DC,11.2009  
c
c pt: perturbation potential temperature
c pts: reference state potential temperature
c tems: reference state temperature
c tsuf: temperature at the surface
c psuf: surface pressure
c ptop: pressure at the top of the model
c pp=psuf-ptop
c hsuf: height of the lower surface
c
c psufi,tsufi: psuf,tsuf in the initial state
c
c data for calculation of bulk transfer in the lowest hakf layer
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c diffusion terms
c
c parameters for rayleigh relaxation in the absorption layer
c
c      real(kind(0.d0)),dimension(0:ns)::
c     :  endrag,taudra,deltazt
      real(kind(0.d0)),dimension(:),allocatable::
     :  endrag,taudra,deltazt
      integer:: idrmax
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c phi: perturbation geopotential
c phis: reference state geopotential
c phisuf: phis at the surface
c phig:
c
c pgradx,pgrady,pgrads: boundary values of pp gradient (auxiliary var)
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c grid variables
c
c      real(kind(0.d0)),dimension(0:ns1)::
c     :  sigma0,sigma1
c     :  ,ds0,ds02,ds04
c     :  ,ds02a,ds04a,ds08a
c     :  ,ds1,ds12,ds14
c     :  ,s1ds0,s0ds1,s0ds2a,s0ds4a
c     :  ,s1ds4a
c     :  ,ds1x,ds1dx,ds1y,ds1dy
c     :  ,dtxs,dtys,dxys
      real(kind(0.d0)),dimension(:),allocatable::
     :  sigma0,sigma1
     :  ,ds0,ds02,ds04
     :  ,ds02a,ds04a,ds08a
     :  ,ds1,ds12,ds14
     :  ,s1ds0,s0ds1,s0ds2a,s0ds4a
     :  ,s1ds4a
     :  ,ds1x,ds1dx,ds1y,ds1dy
     :  ,dtxs,dtys,dxys

      real(kind(0.d0))::
     :  dsng,dsr05
     :  ,dx2,dx4,dx8,dx05,dx025,dxx,dxdt2,dxt2
     :  ,dy2,dy4,dy8,dy05,dy025,dyy,dydt2,dyt2
     :  ,dxy,dxy025,dtxy
     :  ,dtl
     
c     parameters regarding wind speed perturbation

      logical uvpert    !   switcher
      real f_chord, xm_pert  ! focal chord length (km) and perturbation scale (ms-1)
      real*8 tau_nudg   ! characteristic time scale of nudging
      integer:: nperturb,pnudg    ! start and length of nudging
c

      integer:: nstep
c      real(kind(0.d0)):: delka2(ns)
      real(kind(0.d0)),allocatable:: delka2(:)
c
c fcor: coriolis parameter
c
      real(kind(0.d0)):: fcor

c
c utm coordinates
c
      real(kind(0.d0)) xminutm,xmaxutm,yminutm,ymaxutm
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c parameters for fft
c
      integer::ifax(10),ifay(10)
c      real(kind(0.d0))::trigsx(nx-1),trigsy(ny-1)
c      real(kind(0.d0))::
c     :  sipcox(nx-1),simcox(nx-1),sipcoy(ny-1),simcoy(ny-1)
c      real(kind(0.d0))::
c     :  xlamb(nx-1),ylamb(ny-1),xylamb(nx-1,ny-1)
c     :  ,eigxy(2:nx,2:ny)
      real(kind(0.d0)),allocatable::trigsx(:),trigsy(:)
     :  ,sipcox(:),simcox(:),sipcoy(:),simcoy(:)
     :  ,xlamb(:),ylamb(:),xylamb(:,:)
     :  ,eigxy(:,:)
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c work space
c
c note: wm2 needed from moment to wsigbc
c note: uflux,vflux,wflux needed from flux to moment next time step
c
      real(kind(0.d0)),allocatable::
     :  theta(:,:,:)
     :  ,hkk1(:,:),hkk2(:,:)
     :  ,hkk3(:,:),hkk4(:,:),hkk5(:,:),
     :  pressure(:,:,:),temperatur(:,:,:)
      integer,allocatable::ihkk1(:,:),ihkk2(:,:)

      real(kind(0.0d0)),allocatable::
     :  rh(:,:,:),ri(:,:,:),bosa(:,:,:)
     :  ,work1(:,:,:),work2(:,:,:),work3(:,:,:),work4(:,:,:)
     :  ,work5(:,:,:),work6(:,:,:),work7(:,:,:),work8(:,:,:)

c      real(kind(0.0d0)),dimension(0:nx1,0:ny1,0:ns1)::
c     :  uflux,vflux,wflux
c     :  ,def11,def12,def13,def22,def23,def33
c     :  ,h1,h2,h3,q1,q2,q3,dif000,dif010,dif100,ri,zed,dif001
c      real(kind(0.d0)),dimension(0:ns1):: w1d1,w1d2
      real(kind(0.0d0)),allocatable::
     :  dif000_h(:,:,:),dif010_h(:,:,:)
     :  ,dif100_h(:,:,:),dif001_h(:,:,:)
      real(kind(0.d0)),dimension(:),allocatable:: w1d1,w1d2

c     integer,dimension(0:nx1,0:ny1,0:ns1):: iwk1

c      equivalence(wk1,vfy,wfx,def11,h1,q1,zed)
c     equivalence(wk1,ufx,vfx,wfx,def11,h1,q1,zed)
c     equivalence(wk2,ufy,vfy,wfy,def12,h2,q2)
c     equivalence(wk3,ufs,vfs,wfs,def13,h3,q3)
c     equivalence(wk4,wsigp,def22)
c      equivalence(wk6,uflux)
c      equivalence(wk7,vflux)
c      equivalence(wk8,wflux)
c     equivalence(wk3,wm2),(wk9,ptm2)

      integer,allocatable::iworkk(:,:,:)
c     dimension iwk1(:,:,:)
c     equivalence(wk3,wm2),(wk9,ptm2),(wk1,iwk1)
c     equivalence(wk6,uflux),(wk7,vflux),(wk8,wflux)
c
c qv - water vapour
c qc - cloud water (liquid)
c
       real(kind(0.d0))::qif,qdrag
       integer::ifqc,ifqr,ifblob,ifqi,ifperturb,                         DC,10.2009
     : ifdirichlet
c
c soil
c fields for solonoi
c also radiation and precipitation and geostrophic wind
c rainacc: accumulated precipitation during the run
c
c      real(kind(0.d0)),dimension(0:nx1,0:ny1)::
c     :  veg
c     :  ,xlai,rsmin,alb
c     :   ,emis,z0,z0h
c     :   ,rra,z0hw,z0w
c     :   ,rg,rat,xlake
c     :   ,pl,zonirr
c     :   ,tsmed,tslake
c     :   ,tsnoi,t2noi,wgnoi
c     :   ,w2noi,wrnoi
c     :   ,tsurw,tmsur
c     :   ,rn,h,le
c     :   ,gsolo,xdd2
c     :   ,cdm,rainacc
c
c      integer,dimension(0:nx1,0:ny1)::
c     :  iclay,isand,iveg

      real(kind(0.d0)),allocatable,dimension(:,:)::
     :  veg,xlai,rsmin,alb,emis,z0,z0h,rra,z0hw,z0w,rg,rat,xlake
     :   ,pl,zonirr,tsmed,tslake,tsnoi,t2noi,wgnoi,w2noi,wrnoi
     :   ,tsurw,tmsur,rn,h,le,gsolo,xdd2,cdm,rainacc,watice,                    DC,03.2009
     :    snowacc,xsea,tsfix,imp                                                    DC,11.2009
      integer xice

      integer,dimension(:,:),allocatable::
     :  iclay,isand,iveg

      real(kind(0.d0)),allocatable::
     :  xpl(:),xra(:),xrg(:)

      integer::
     :   li,lf,nf,ndados,ifsoil,ifhle
	integer::                                                                 DC,06.2010
     :   tfix                                                                   DC,06.2010

      real(kind(0.d0))::
     :   ts(2),t2(2),wg(2),w2(2),wr(2),rhoa,z0hz0,timeini

      real(kind(0.d0)),allocatable,dimension(:,:)::
     :  temmax,temmin,evapor,ptforc
     :  ,phi00s

      real(kind(0.d0)),allocatable::
     :  ugeos(:),vgeos(:),xugeos(:,:),xvgeos(:,:)

      real(kind(0.d0)):: fcori

c      real(kind(0.d0)),dimension(0:nx1)::
c     :  tauspo
c      real(kind(0.d0)),dimension(0:ny1)::
c     :  tauspoy

      real(kind(0.d0)),allocatable,dimension(:)::
     :  tauspo,tauspoy

      real(kind(0.d0)), allocatable::
     :  a1phi(:),a2phi(:,:),a3phi(:),ctphi(:,:),dtphi(:,:)
     :  ,sine(:),sines(:,:)

c      real(kind(0.d0)),dimension(ns-1)::
c     :  a1pos3,a3pos3

c      real(kind(0.d0)),dimension(nx-1,ny-1)::
c     :  a21pos3
      real(kind(0.d0)),allocatable,dimension(:)::
     :  a1pos3,a3pos3
      real(kind(0.d0)),allocatable,dimension(:,:)::
     :  a21pos3

c      real(kind(0.d0)),dimension(0:ns)::
c     :  aves
      real(kind(0.d0)),allocatable,dimension(:)::
     :  aves

      real(kind(0.d0)), allocatable::
     :  errp(:,:)
      integer, allocatable::
     :  idath(:,:)

c      real(kind(0.d0))::
c     :  duul(ny1),duur(ny1),dvvl(nx1),dvvr(nx1)
c
c      real(kind(0.d0))::
c     :  phibx1(0:ny1),phiby1(0:nx1),phibx2(0:ny1),phiby2(0:nx1)
      real(kind(0.d0)),allocatable,dimension(:)::
     :  duul,duur,dvvl,dvvr
     :  ,phibx1,phiby1,phibx2,phiby2

      real(kind(0.d0)),allocatable::
     :  a1sphi(:),a2sphi(:),a3sphi(:),csphi(:),dsphi(:)
     :  ,hk11(:,:),hk12(:,:)

c      real(kind(0.d0))::
c     :  dphsdx(2,0:ny1),dphsdy(0:nx1,2),hdamp(0:ns1)
      real(kind(0.d0)),allocatable::
     :  dphsdx(:,:),dphsdy(:,:),hdamp(:)

      real(kind(0.d0)), allocatable::
     :  ff(:,:),tbot(:,:),dlnpdx(:,:),dlnpdy(:,:)
     :  ,dlntdx(:,:),dlntdy(:,:)
     :  ,dpsudx(:,:),dpsudy(:,:),hk21(:,:,:),hk22(:,:,:)
c
      real(kind(0.d0))::
     :  dragm,dragt
      real(kind(0.d0)), allocatable::
     :  zlevm(:),zlev(:)
     
      integer(4) iflake,radpar,rat_par,nradcall,nlakecall                       VS,05.2007
     & ,shortwave,longwave                                                      VS,26.06.2007
      real(8) month,day,tau_rg                                                  VS,05.2007
      real(kind(0.d0)), dimension(:,:,:), allocatable::                         VS,05.2007
     & Radheat,Srad,Lrad                                                        VS,05.2007
      real(kind(0.d0)), dimension(:,:), allocatable::                           VS,05.2007
     & Srad_surf,Lrad_surf                                                      VS,05.2007
      real(kind(0.d0)), dimension(:,:), allocatable::                           VS,06.2007
     & tsurf                                                                    VS,06.2007
      real(kind(0.d0)), dimension(:,:), allocatable::                           DC,06.2010
     & ust_s,dzits,Fv,wm,wstar,gammah,tst_s,qst_s,
     & gammaq                                                                   DC,06.2010
      integer imonth,iday,ihour,iminu,iseco                                     VS,05.2007

      real(kind(0.d0))::                                                        VS,12.08.2007
     & t1_cpu,t2_cpu,t_exec(1:100)                                              VS,12.08.2007


c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c options:
c
c hydro : if not.true. phi is obtained from the elliptic equation.
c
c iophis = 1(hydrostatic) 2(poisson eq dirichlet) 3(poisson eq neuman)
c ioppuv = 1(u v define boundary dpp/dt) 2(rad.b.c. for pp define u v)
c iobxxx =
c iowlbc = 0(lateral b.c. for w from wsig, and wsig b.c. from contin.)
c          1(lateral b.c. for wsig from w, and w from u and v b.c.)
c
      logical hydro,tfct,inipt,raylei,olddra,mxspan,adjust,centered,dst3

      real(kind(0.d0))::
     :   dx,dy,dt,xl,yl,dt2
     :  ,xlcor,xrcor,ylcor,yrcor
      logical fluxxlcor,fluxylcor,fluxxrcor,fluxyrcor
      integer::
     :  ntime

      integer::
     :  nupsur,nuppts,nupdif
     :  ,iophis,ioppuv,iowlbc
     :  ,iodify,iobphy,iobuux,iobuuy,iobvvx,iobvvy,iobptx,iobpty
     :  ,iobppx,iobppy

c
c mountain: iomtyp:1- 2d mountain; 2- 3d mountain
c ivarhs=0 mountain height is established at nstep=0 and is fixed
c ivarhs=1 mountain height grows linearly from nstep=0 to nstepgrow and stays fixed afterwards
c ivarhs=2 mountain height grows in steps
c
      integer::
     :  iomtyp,ivarhs,nstepgrow,nstepflat,numsteps
      real(kind(0.d0))::
     :  hmount,xmount,ymount,xspan,yspan,hxwid,hywid,hexp
     :  ,xmountc,ymountc
c
c aselin time filter
c
      real(kind(0.d0))::
     :  filt,ppfilt
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      real(kind(0.d0)):: r,cp,g,hlat,p00,rv,hlatcp,xk4,qco,xnor,rhol
     :  ,xk1,dv,xkt,xniu,rho0,xa,xb,sch,gm48,gm38,gm29,akapa,gr,g2,r05
     :  ,gr2,g4,g05,omega,eps,ttmin,nso,hsub,dsnow,csnow,nro,arain
     :  ,hfus                                                                   DC,10.2009
      parameter(r=287.,cp=1005.,g=9.8066,hlat=2.501e6,p00=1.e5,
     :hsub=2.837e6,nso=3.e6,dsnow=0.25,csnow=11.72,hfus=3.336e5
     :,brain=0.8)                                                              DC,10.2009
      parameter(rv=461.51,hlatcp=hlat/cp,xk4=1.e-3,qco=1.e-3,xnor=0.8e7)
      parameter(rhol=1.e3,xk1=1.506,dv=2.26e-5,xkt=2.43e-2
     : ,xniu=1.51e-5,rho0=1.23,xa=842.,xb=0.8,sch=xniu/dv
     : ,gm48=17.837870,gm38=4.694174,gm29=1.827355)
c
c note: ha is the height corresponding to the surface pressure pa.
c
      parameter(akapa=r/cp,gr=g/r,g2=2.*g,r05=0.5*r,gr2=gr*2.,g4=4.*g)
      parameter(g05=0.5*g)
      parameter(omega=7.292e-5)
      parameter(eps=1.e-6,ttmin=1.e-5)

      real(kind(0.d0)):: xlatit,ptop,pa,ha
      common/par301/xlatit
      common/par302/ptop,pa,ha
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c iodif = 1(constant dif.coef.) 2(deformation dependent)
c
      logical:: dohsmo,dovsmo,perdam,dovsmt
      real(kind(0.d0))::
     :  hdampm,hdampt,vdamp,vdampt
      integer::numsmo

      real(kind(0.d0))::
     : rikey,difcof,uvdif,tsdif,difl
      integer::
     :  iodif
     :  ,nxsponge1,nxsponge2,nysponge1,nysponge2
     :  ,ntaut,ntauspon
      real(kind(0.d0))::
     :  cdcoef
     :  ,taum,taut,zm,zt,tauspon
      parameter(tdif=1.d0,qsdif=0.d0,qdif=1.d0,dkapa=0.8d0)

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     aux
      real*8 hh, declive, planalts, planalti


c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      integer npro
      real(kind(0.d0))
     :   t0ini,ug,vg,dugdz,dvgdz,psref,dzpro

      real(kind(0.d0)),allocatable, dimension(:)::
     :   zpro,fpro,fnpro,press0,press1,dpressdt,ini_pts
      real(kind(0.d0)),allocatable, dimension(:,:)::
     :   thpro0,dthdz,dthdt,dthdzdt
     :   ,uspro0,dusdz,dusdt,dusdzdt
     :   ,vspro0,dvsdz,dvsdt,dvsdzdt
     :   ,qvspro0,dqvsdz,dqvsdt,dqvsdzdt
     :   ,thpro1,dthdz1,uspro1,dusdz1
     :   ,vspro1,dvsdz1,qvspro1,dqvsdz1
     :   ,tepro,presspro,rhopro,rhopro1,
     :   perpro0,dperdz
c    :   ,thpro2,dthdz2,uspro2,dusdz2
c    :   ,vspro2,dvsdz2,qvspro2,dqvsdz2
c
c ndat is the maximum number of points in the profile given
c ndat=max(ndth,ndus,ndvs,ndqvs)
c
      integer ndat

      real(kind(0.d0)),allocatable,dimension(:,:)::
     :   thdat,usdat,vsdat,qvsdat
     :   ,zthdat,zusdat,zvsdat,zqvsdat,psdat,ptdat,pressdat
     :   ,tedat,per_dat

      character*50:: exptit
      character*50:: comlin(3)
      logical:: itepsuf
c
      integer::
     :  ioreft,iorefq,iorefu,iorefv,ndth,ndus,ndvs,ndqvs
      real(kind(0.d0)):: pts0

c      common/dat02/pts0,ioreft,iorefq,iorefu,iorefv,ndth,ndus,ndvs,ndqvs
      real(kind(0.d0))::
     :  ug0,dugdz0,vg0,dvgdz0
      common/dat04/ug0,dugdz0,vg0,dvgdz0,itepsuf
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c input/output control:
c iotfor: period of formatted output (unit nchn1);
c iotunf: period of unformatted output (unit 31, with initial formatted
c data in unit 30);
c iotrst: period of output of data for restart (units 40+).
c
c restar: if .true. the model starts from data of a previous run
c which must be provided on unit 21 (transfer to file ft21):
c activate corresponding card !!!                                       #exp#
c
      logical restar
c
c if rstout is true then:
c restart data output every iotrst time steps and at ntime.
c
      logical rstout
c
c forout: if true there is formatted output each iotfor time steps
c and at nstep=ntime.
c
      logical forout
c
c grdout if true the formatted output will include a grid file
c
      logical grdout
c
c unfout: if true there is unformatted output each iotunf time steps
c and at nstep=ntime. this output, on unit 30,31 must be transfered
c for a permament disk file, after being packed.
c activate corresponding cards !!!
c
      logical unfout
c
c output of wave drag:
c
      logical drgprt,domom
c
c choice of formatted output (see out6):
c
      logical norout,morout,refout,fluout,difout,surout
     :   ,iniout,masout,driout,ddiout,defout,eldout,elfout,elgout
     :   ,phsout,desout,echeck,eneout,ppbout,uubout,ptbout
     :   ,duvout,inrout,iniunf,proout
c
c unit for formatted output:
c
      parameter(nchn1=28)
c
c verbose: -1(no screen output), 0(a few), 1(much), 2(all
c
      integer verbose
c
      logical prt,unfprt
      real(kind(0.d0))::
     :  dzlev
      integer::
     :  iotrst,iotfor,iotunf,iotdrg,iotduv,nlev,iomod
     :  ,iostart,ioend,iodelta,iotscrout

c      parameter(noutum=100)
c
c local output
c
c      integer::
c     :  noutu,ixoutu(noutum),iyoutu(noutum),isoutu(noutum)
      integer noutu
      integer,allocatable,dimension(:)::ixoutu,iyoutu,isoutu
c
c profile output
c
c      parameter(noutpm=100)
c      integer:: noutp,ixoutp(noutpm),iyoutp(noutpm)
      integer noutp
      integer,allocatable,dimension(:):: ixoutp,iyoutp
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c parameters for radiative boundary conditions:
c
      real(kind(0.d0)):: rdx,rdy,rdt,dxdt,dydt,rxmax,rymax
      common/radpar/rdx,rdy,rdt,dxdt,dydt,rxmax,rymax
      real(kind(0.d0)):: pi,pi05
      common /consta/ pi,pi05

c
c old nh3dprec.inc
c

      real(kind(0.d0)):: zero0,resmax
      parameter (zero0=1.d-30,resmax=1.d-12)

      integer::istep,kstep
c      real(kind(0.d0)),dimension(0:nx1,0:ny1)::
c     :  pequiv
      real(kind(0.d0)),dimension(:,:),allocatable::
     :  pequiv
c
c work variables for inpos3c/pos3c
c
      real(kind(0.d0)),allocatable::emo(:,:),emi(:,:),vlam(:),e(:)

c
c variables for out6
c
      real(kind(0.d0)),allocatable,dimension(:,:)::panom,pini

      integer ifluxcor
      real(kind(0.d0)) dpmeddt

      end module alloc


      module nh3dpa2
c
c data for wrigar and wgrids
c
      parameter(nplan=9)
      logical fixexp
c
      integer ixp(nplan),iyp(nplan),isp(nplan)
      real(kind(0.d0))::zplan(nplan)
      real(kind(0.d0))::scalew
      parameter(maxxy=1026)
      real(kind(0.d0))::zbot(0:maxxy)

      logical zconst
c
      parameter(ndigit=07,nunit=28,leng=132)
      common/wri01/ixp,iyp,isp,fixexp,nexfix,ixss,iyss,isss

      end module nh3dpa2

      module refstate

c
c icorprof selects coordinates of the given profiles
c =0 model coordinates
c =1 utm coordinates
c =2 lon,lat - to convert to utm

      integer nprof,intrefst,mrefst,ifuprefs,icoorprof
      real(kind(0.d0)),dimension(:), allocatable::
     :  xprof,yprof,distinv
      real(kind(0.d0)) tprof1,tprof2,starttime

      real(kind(0.d0)),dimension(:,:,:), allocatable:: pweight
     :  ,dweightx,dweighty
      logical realtime,forcedpp,forcepsuf,forcemflx,divrefstate

      end module refstate
      
      MODULE optic_parameters
!-----Stores optic parameters, required for shortwave and longwave 
!-----radiation parameterization;
!-----They are updated each timestep in opticparset, 
!-----and used as input in SORAD and IRRAD subroutines    
      real, dimension(:,:,:), allocatable:: taual,ssaal,asyal,taucld,
     & reff
      real, dimension(:,:,:),allocatable::taucl_lw, eg,ev, rvir
      real,dimension(:,:,:,:),allocatable::taual_lw,ssaal_lw,asyal_lw
      real, dimension(:,:), allocatable:: oa
      real, dimension(:), allocatable:: rsuvbm,rsuvdf,rsirbm,rsirdf
      real s0,co2,n2o,ch4,cfc11,cfc12,cfc22
      integer, parameter:: na = 3, nsur = 2
      character(len=1) crel
!      logical high
      
      contains
       subroutine alloc_optics(m,np)
       integer(4) m, np
       
       allocate (oa(m,np))
       allocate (taual(m,np,8),ssaal(m,np,8),asyal(m,np,8),
     &  taucld(m,np,3), reff(m,np,3))
       allocate (rsuvbm(m),rsuvdf(m),rsirbm(m),rsirdf(m))
       allocate (taucl_lw(m,np,3), taual_lw(m,np,10,na),
     & ssaal_lw(m,np,10,na),asyal_lw(m,np,10,na))
       allocate (eg(m,nsur,10),ev(m,nsur,10),rvir(m,nsur,10))
       
!       s0 = 1367.
!       s0 = 1372.45      
       
       end subroutine alloc_optics
      END MODULE optic_parameters      
