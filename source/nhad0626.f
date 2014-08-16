
C THIS IS THE MAIN FILE FOR THE NH3D MODEL
C VERSION NHaD06 BUILD 12
c
c To use the model:
c . Compile (f90)
c . Edit file nh3d.dat
c . run
c
c To change the output for particular ploting packages:
c    . post process unformatted files (UNFOUT=.T.)
c or . post process grd and vector files (GRDOUT=.T.)
c or . change subrout. WRIGRD and WRIVEC to suit your needs
c
c GOOD LUCK !
c
c
c Pedro Miranda 1988, May 1995,August 1995, December 1995, June 1997
c Last changed 2005/10
c Rein Room 1997, 2000
c Developed at: University of Reading, University of Lisbon,
c Tartu Observatory
c
c History:
c sep 1988 onwards:
c
c    .3d code
c    .variable ds grid
c    .option of dirichlet or neuman boundary conditions for phisuf
c    .double staggered cosine transform for solution of poisson eq.
c    .restart facility
c
c new 1994/95:
c    .(optional) water vapour
c    .(optional) surface scheme of Noilhan et Planton (1989)
c    .(optional) radiation forcing at the surface (imposed)
c
c    ..several important input/output changes
c    ..all capital letters have been eliminated
c    ..initialization with f works alright (any wind direction)
c    ..(optional) lateral sponge
c
c new 1996
c    .condensation (water vapour and liquid water)
c     (follows Klemp and Whilhelmson (1978), Yau and Austin (1979))
c
c new 1997
c    . dynamical allocation for some variables
c    . better memory management (water varibles can be switched off easily
c                                without memory allocation)
c    . other fortran90 features
c    . user friendly input file (nh3d.dat)
c    . direct grid output from the model
c
c new 1997 - nhad model
c
c   Revised by R. Room (15.03.97):
c   (1) lower boundary p_* fixed;
c   (2) diagnostical computation of vertical velocities W and WSIG
c   (3) nondivergent initial velocity
c
c   Revised by R. Room (19.04.2000):
c    (1) inversion of ellipt using eigenvector techique in vertical:
c      ellipt3, inpos3c, inver, tqli
c    (2) no boundary conditions for NH geopotential at bottom and top,
c        rather a detailed balance condition for vertically integrated mass
c        employed
c
c If ADJUST=.T., modified version by R. R., else original version by P. M.
c
c 1998
c   Included an option for vertical diffusion only: in that case
c   IODIF=11 (vertical diffusion with constant dif. coef)
c   IODIF=12 (vertical diffusion, deformation dependent coef.)
c   IODIF=10 is equiv. to IODIF=0 (No diffusion)
c   This option was included because in the case of high mountains
c      there was an instability in the diffusion computation
c      due to the metric terms in the horizontal derivatives
c
c 2000 - November-December
c    Included an output for a numerical lagrangean tracer
c    New mountain type, sinusoidal sequence, iomtyp=7
c
c 2001 - January
c    ifsoil=2 new option included: turns off heat and latent heat fluxes, while keeping
c             the momentum flux
c    New mountain types, upwind ramp (iomtyp=8), downwind ramp, (iomtyp=9)
c    symetric ramp (iomtyp=10)
c
c 2001 - February
c    Included an option for vertical diffusion using Teixeira's explicit stable method
c    IODIF=20 is equivalent to IODIF=0 (No diffusion)
c    IODIF=21 (Teixeira's vertical diffusion with constant dif coef)
c    IODIF=22 (Teixeira's vertical diffusion with deformation dependent dif coef)
c
c
c Jun 2005
c    isqv=qif, isqr=ifqr, isqc=ifqc imposed on readpa (only one is needed)
c
c
c 2005
c nhad0517
c    full dynamic allocation: no need for recompilation between runs,
c      everything may be controlled by nh3d.dat
c    various minor bug corrections
c nhad0523
c    3D time varying forcing (u,v,theta): "xprofile"
c    simplified time varying profile (xprofile in nh3d.dat, the profile is given in *.ref)
c    the 'ref' file also includes the nprof profiles in different locations at different times,
c    the model interpolates those profiles linearly in time and by inverse distance weighting
c    in the horizontal
c    to activate the time varying profile, it is necessary to turn on the lateral sponge
c       put nxsponge1,nxsponge2,nysponge1,nysponge2 >0 (e.g. 5-10)
c       and define ifuprefs=1
c
c
c FILES:
c
c Data files:
c
c nh3d.dat (Input file)
c     This file includes all parameters which may be changed for an
c     experiment.
c     See the sample file for more information
c     The file includes:
c         - NH3D (first line header)
c         - comments>> lines starting by #
c         - single parameter setting>> keyword (variable name)
c             followed by value in '*' format (comments may follow)
c         - group of data>> keyword, then several lines of data
c             - no comment lines are allowed inside groups
c         - END (marks last line)
c      The order of parameters is irrelevant, except within a group
c      Most parameters have default values and may be omitted
c      The model checks if all needed parameters were given
c      See sample file (samp.dat) for an explanation of all parameters
c
c initial profiles and boundary forcing (nhad0520 onwards)
c   the name of this file is specified in nh3d.dat (keyword fnprofile, see below)
c
c other files include surface information and depend on switchs.
c In the simplest case (analyical
c orography, no surface scheme) no other files are needed.
c
c All input/output file names are organized in 5 groups with
c automatically made names:
c
c Group 1: surface files: Name format [fnmap].ext
c                   e.g. "samp_map.???"
c                   samp_map.top (topography, if iomod=4)
c                   samp_map.xdd (soil depth, if ifsoil=1)
c                   samp_map.ts  (soil initial temperature, if ifsoil=1)
c                   samp_map.tsm (soil temperature for function tsoil: "analytical" soil)
c                   samp_map.t2  (deep soil initial temperature, if ifsoil=1)
c                   samp_map.wg  (soil initial liquid water, if ifsoil=1)
c                   samp_map.w2  (deep soil initial liquid water, if ifsoil=1)
c                   samp_map.wr  (initial liquid water in the leaves, if ifsoil=1)
c                   samp_map.z0  (surface rugosity length, if ifsoil=1)
c                   samp_map.veg (fraction of vegetation, if ifsoil=1)
c                   samp_map.ive (vegetation type, if ifsoil=1)
c                   samp_map.lai (leaf area index, if ifsoil=1)
c                   samp_map.rsm (stomatal resistance, if ifsoil=1)
c                   samp_map.alb (surface albedo, if ifsoil=1)
c                   samp_map.icl (% of clay in the soil , if ifsoil=1)
c                   samp_map.isa (% of sand in the soil , if ifsoil=1)
c                   samp_map.irr (precipitation map to be multiplied by a constant at each time step , if ifsoil=1)
c
c Group 2: files with external forcing: Name format "samp_for.???"
c          currently this includes only reference state information
c                   samp_for.ref  (this is the fnprofile file)
c
c Group 3: files with output: Name format [fnout].ext
c                   e.g. "samp_001.???"
c                   samp_001.u10 (unformatted, hsuf,psuf)
c                   samp_001.u11 (unformatted, u)
c                   samp_001.u12 (unformatted, v)
c                   samp_001.u13 (unformatted, w)
c                   samp_001.u14 (unformatted, pt)
c                   samp_001.u15 (unformatted, phi)
c                   samp_001.u16 (unformatted, phis)
c                   samp_001.u17 (unformatted, wsig)
c                   samp_001.u18 (unformatted, difunu)
c                   samp_001.u19 (unformatted, difunv)
c                   samp_001.u20 (unformatted, qv) (if qif.ne.0.)
c                   samp_001.usf (unformatted, surface fields) (if.ifsoil.eq.1)
c                   samp_001.fmd (formatted, parameters)
c                   samp_001.cdm (unformatted, drag coefficient) (if.ifsoil.eq.1)
c                   samp_001.dxy (formatted, drag)
c                   samp_001.duv (formatted, momentum flux, constant sigma)
c                   samp_001.mom (formatted, momentum flux, constant z)
c                   samp_001.rst (unformatted, output file for restart, if rstout=.true.)
c                   samp_001.out (formatted, general output)
c                   samp_001.loc (formatted, pressure in selected points (for "microbaropraphy")"
c                   samp_001.sur (formatted, all surface values in selected points)
c                   samp_001.pro (formatted, all data in selected profiles)
c                   samp_001.eva (formatted, total evaporation during run)
c                   samp_001.tmi (formatted, min temperature during run)
c                   samp_001.tma (formatted, max temperature during run)
c                   samp_001.rai (formatted, accumulated surface precipitation during run)
c
c Group 4: restart file: Name format [fnrst].rst
c                   e.g. "samp_rst.???"
c                   samp_rst.rst (restart from, if restart=.true)
c   Note: the restart system needes some rewriting after the major change to dynamical allocation
c
c
c Group 5: grid files: Name format [fngrd]vvmpp.ttt in GRD format(see below)
c                   where: vv - variable
c                          m  - mode (s-const sigma,y,x,z -interpolated to a z surface)
c                          pp - plane (ix,iy,is or z in 100m)
c                   variable list: vv=
c                                     u_ (u)
c                                     v_ (v)
c                                     w_ (w)
c                                     pt (pt)
c                                     th (pt+pts=theta)
c                                     qv (qv)
c                                     qc (qc)
c                                     qr (qr)
c                                     co (cond)
c                                     ev (evap)
c                                     ph (phi)
c                                     ws (wsig)
c                                     rh (relative  humidity)
c                                     hh (sensible heat flux)
c                                     lh (latent heat flux)
c                                     gs (soil heat flux)
c                                     rn (net radiaton flux)
c                                     cd (transfer coefficient)
c                                     ps (surface pressure)
c                                     pr (precipitation)
c                                     fs (phis)
c                                     os (pts)
c                                     tt (tems)
c                                     uf (uflux)
c                                     vf (vflux)
c                                     wf (wflux)
c                                     s_ (s)
c                                     t0 (tsuf)
c                                     ts (tsnoi)
c                                     t2 (t2noi)
c                                     wg (wgnoi)
c                                     w2 (w2noi)
c                                     wr (wrnoi)
c                                     du (difunu)
c                                     dv (difunv)
c                                     dw (difunw)
c                                     dt (difunt)
c                                     dq (difunq)
c                                     dc (difunqc)
c                                     dr (difunqr)
c                                     ri (richardson number)
c
c The suggestion for the choice of the names (the extensions are
c automatic) is the following:
c       "samp" is to be substituted by a general name of the experiment
c       "001" identifies a particular run
c
c Note that all input surface files consist of xy grids which are
c read by subroutine readgrd. They must be written in the "grd"
c format (see one from the sample) which is:
c
c DSAA (first line)
c nx,ny (number of columns, number of lines)
c xmin,xmax
c ymin,ymax
c zmin,zmax
c then data values by line (constant y) from lower left corner.
c
c All grid formatted output is in the same format (see subrout. wrigrd)
c-----------------------------------------------------------------------
c
      
C
      program nhad
c
c staggered grid:
c
c physical domain: x1(1 to nx)        x0(1.5 to nx+0.5)
c                  y1(1 to ny)        y0(1.5 to ny+0.5)
c                  sigma1(1 to ns)    sigma0(0.5 to ns-0.5)
c
c ********* plane y1(iy) *********
c
c sigma0(is-1) -
c
c sigma1(is) -        011
c                     vfs,wfy
c                     def23,dif011
c
c
c
c sigma0(is) :        010                                     110
c                     v,difunv,vflux,vs                       ufy,vfx,
c                     dif010,h2                               def12
c                                                             dif110
c
c sigma1(is+1)
c
c                     |                                       |
c                     x0(i)                                   x1(i)
c                     pp01,ppdy01
c
c ********* plane y0(iy) ********* (iy=1,ny1)
c
c sigma0(is-1) -
c
c sigma1(is) -        001                                 101
c                     w,s,wsig,difunw,                        ufs,wfx,
c                     wflux,h3,                               def13
c                     dif001                                  dif101
c
c
c sigma0(is) :        000                                     100
c                     phi,phis,pt,pts,tems,qv,qvs,qc,qr,      u,difunu,
c                     difunt,wfs,def11,def22,def33,def,ri     uflux,us,
c                     ufx,vfy,wfs,dif000                      h1,flx
c                                                             dif100
c sigma1(is+1)        fly
c
c                     |                                       |
c                     x0(i)                                  x1(i)
c                     pp,dpp,psuf,hsuf,ppdx,ppdy             pp10,ppdx10
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      use alloc
      use refstate

      implicit real*8 (a-h,o-z)
c      dimension zup(0:nx+1,0:ny+1)
      real(kind(0.d0)),dimension(:,:),allocatable:: zup

      logical scrout,dif3d,diftx
	logical difloc,nonloc_tm86,nonloc_ls96                                    DC 06.2010
      character*80 fnmap,fnfor,fext,fnout,fngrd
      character*14 timestring
      common/cfnout/fnout
      real time_begin, time_end

      pi=4.d0*datan(1.d0)
      call cpu_time(time_begin)

      verbose=0
      istep=0
      kstep=-1
      write(*,*) 'Begin NHAD'
      write(*,*) 'Vers. nhad06, Build 21, P.Miranda, R. Room 2000ff'
      write(*,*)
     :  'Developed at the Univ.Reading, Lisbon and Tartu'

      if(verbose.ge.1) write(*,*) 'call readpa'
      call readpa(fnmap,fnfor,fnout,itheta
     :   ,iwind,scrout,fngrd,dif3d,diftx,difloc,nonloc_tm86,nonloc_ls96)
      if(forout) prt=.true.

      if(verbose.ge.0) write(*,*) 'ADJUST=',adjust
      write(nchn1,*) ' output from nhad'
      write(nchn1,*) 'Vers. nhad06, Build 21, P.Miranda, R. Room 2000ff'
      write(nchn1,*) 'experiment: ',exptit
      write(nchn1,*) comlin(1)
      write(nchn1,*) comlin(2)
      write(nchn1,*) comlin(3)
      call allocvar(nonloc_tm86,nonloc_ls96)
      allocate(zup(0:nx1,0:ny1))

!     Dima@ 13.05.07
      if (impur.gt.0) then                                                      DM,05.2007  
        call impursetup(nx1qv, ny1qv, ns1qv,                                    DM,05.2007
     :    ix_impur, iy_impur, impur_value,                                      DM,05.2007
     :    immediate_impur, constant_impur,                                      DM,05.2007
     :    le_impur, qs_10, qs_10s)                                              DM,05.2007
      endif                                                                     DM,05.2007

      if(verbose.ge.1) write(*,*) ifuprefs,itheta,iwind

      zlevm(1)=hmount
      do ilev=2,nlev
        zlevm(ilev)=zlevm(ilev-1)+dzlev
      enddo

      zlev(0)=0.
      do ilev=1,nlev
        zlev(ilev)=zlev(ilev-1)+dzlev
      enddo

      ug=ug0
      vg=vg0
      dugdz=dugdz0
      dvgdz=dvgdz0
      dragx=0.
      dragy=0.

      if(verbose.ge.1) write(*,*) 'init0'
      call init0
c
c prepare constants to be used
c
      if(verbose.ge.1) write(*,*) 'const'
      call const
      if(verbose.ge.1) write(*,*) 'infft'
      call infft(ifax,ifay,trigsx,trigsy,sipcox,sipcoy,simcox
     :  ,simcoy,xlamb,ylamb,xylamb,nx,ny,dxx,dyy)
c
      if(restar) then                                                  !*
c                                                                      !*
c note iomtyp and hexp are not being verified                          !*
c                                                                      !*
        call inprst                                                    !*
c                                                                      !*
        write(nchn1,*) '... restarting ...'                            !*
        write(nchn1,*) 'restarting from a previous run at nstep='      !*
     :    ,nstep                                                       !*
        write(nchn1,*) '...'                                           !*
                                                                       !*
        if(verbose.ge.0) then                                          !*
          write(*,*) '... restarting ...'                              !*
          write(*,*) 'restarting from a previous run at nstep='        !8
     :      ,nstep                                                     !*
          write(*,*) '...'                                             !*
        endif                                                          !*
c                                                                      !*
        call refprofs                                                  !*
c       call refprofs(zpro,npro,dzpro,fpro,ndat,nprof,fnpro            !*
c    :    ,iorefu,usdat,zusdat,uspro0,dusdz                            !*
c    :    ,iorefv,vsdat,zvsdat,vspro0,dvsdz                            !*
c    :    ,ioreft,thdat,zthdat,thpro0,dthdz,ptdat,pts0                 !*
c    :    ,iorefq,qvsdat,zqvsdat,qvspro0,dqvsdz                        !8
c    :    ,psdat,pa,verbose,nchn1)                                     !8
        call upuvt                                                     !*
        call rinref                                                    !*
        call indrag                                                    !*
      else
        nstep=0
c
c surface height  hsuf
c
        if(verbose.ge.1) write(*,*) 'sufhit'
        call sufhit(fnmap)
c
c specification of reference state
c
c        call wri3ar(w,2,nx,2,ny,1,ns-1,'w-inref ',tk,iomod,nx,ny)
       
        if(verbose.ge.1) write(*,*) 'inrefs'
        call inrefs
        
       
c
c top enhanced drag or rayleigh friction
c
        if(verbose.ge.1) write(*,*) 'indrag'
        call indrag

c
c specification of initial fields
c
        if(verbose.ge.1) write(*,*) 'initia'
c        call wri3ar(w,2,nx,2,ny,1,ns-1,'w-init0 ',tk,iomod,nx,ny)
c       call initia(wk1,wk2)
        call initia
      endif
c        call wri3ar(w,2,nx,2,ny,1,ns-1,'w-init  ',tk,iomod,nx,ny)
      if(verbose.ge.1) write(*,*) 'done initia'

      tauspo=1.e30
      tauspoy=1.e30

      dragspt=(2.*dx)**2/(4.*pi*pi*tauspon)
      dragspm=(2.*dx)**2/(4.*pi*pi*1.e8)
      do ix=nx1-nxsponge1,nx1
        endrg=dragspm+(dragspt-dragspm)
     :    *cos(0.5*pi*abs(float(nx1-ix)/float(nxsponge1)))**2
        tauspo(ix)=max(2.*dt,(2.*dx)**2/(4.*pi*pi*endrg))
        if(verbose.ge.2) write(*,*)
     :    'endrg1,sponge1:',endrg,tauspo(ix)/dt
        write(nchn1,*) 'Tauspox/dt',ix,tauspo(ix)/(2*dt)
      enddo

      do ix=0,nxsponge2
        endrg=dragspm+(dragspt-dragspm)
     :    *cos(0.5*pi*abs(float(ix)/float(nxsponge2)))**2
        tauspo(ix)=max(dt,(2.*dx)**2/(4.*pi*pi*endrg))
        if(verbose.ge.2) write(*,*)
     :    'endrg2,sponge2:',endrg,tauspo(ix)/dt
        write(nchn1,*) 'Tauspox/dt',ix,tauspo(ix)/(2*dt)
      enddo

      do iy=ny1-nysponge1,ny1
        endrg=dragspm+(dragspt-dragspm)
     :    *cos(0.5*pi*abs(float(ny1-iy)/float(nysponge1)))**2
        tauspoy(iy)=max(dt,(2.*dy)**2/(4.*pi*pi*endrg))

        if(verbose.ge.2) write(*,*)
     :    'endrg1,ysponge1:',endrg,tauspoy(iy)/dt
        write(nchn1,*) 'Tauspoy/dt',iy,tauspoy(iy)/(2*dt)
      enddo

      do iy=0,nysponge2
        endrg=dragspm+(dragspt-dragspm)
     :    *cos(0.5*pi*abs(float(iy)/float(nysponge2)))**2
        tauspoy(iy)=max(dt,(2.*dy)**2/(4.*pi*pi*endrg))

        if(verbose.ge.2) write(*,*)
     :    'endrg2,ysponge2:',endrg,tauspoy(iy)/dt
        write(nchn1,*) 'Tauspoy/dt',iy,tauspoy(iy)/(2*dt)
      enddo
c
c not updated for soil and humidity !!!
c
      if(unfout) then
        call outdat(fnout)
      endif
c
      hnew=hmount
      if(drgprt) then
        if(mxspan) then
          axspan=min(xmount-2.*dx,(xl-xmount-2.*dx))
          ayspan=min(ymount-2.*dy,(yl-ymount-2.*dy))
        else
          axspan=xspan
          ayspan=yspan
        endif
        if(verbose.ge.2) write(*,*) 'call drag0z'
        call drag0z(nchn1,axspan,ayspan
     :    ,nx,ny,ns,psuf,hsuf,dx,dy,dt,nstep
     :    ,exptit,iomtyp,hnew,xmount,ymount,hxwid,hywid,hexp,pa,ptop
     :    ,ioreft,thdat,zthdat,usdat,vsdat,pts0,ivarhs
     :    ,ix00,ix01,ix10,ix11,iy00,iy01,iy10,iy11,verbose,dragx,dragy
     :    ,drag2d)
        if(verbose.ge.2) write(*,*) 'call draguv 0'
        call draguv(ix01,ix11,iy01,iy11,1,ns
     :    ,nx,ny,ns,dx,dy,tems,psuf,u(0,0,0,3),v(0,0,0,3)
     :    ,w(0,0,0,3),nstep,us,vs,ptop,sigma0)
        if(domom) then
          if(verbose.ge.2) write(*,*) 'call momflu 0'
         allocate(iworkk(0:nx1,0:ny1,0:nlev))
c          allocate(zed(0:nx1,0:ny1,0:ns1)) >>wk1
          do is=0,ns
          do iy=0,ny+1
          do ix=0,nx+1
c            zed(ix,iy,is)=(phi(ix,iy,is)+phis(ix,iy,is))/g
            wk1(ix,iy,is)=(phi(ix,iy,is)+phis(ix,iy,is))/g
          enddo
          enddo
          enddo
          call momflu(ix01,ix11,iy01,iy11,1,ns,iworkk,wk1,zlevm,nlev
     :      ,nx,ny,ns,dx,dy,tems,psuf,u,v,w
     :      ,nstep,us,vs,ptop,sigma0)
          deallocate(iworkk)
c          deallocate(zed)
        endif
      endif
c
      write(nchn1,'(//a/)') ' parameters:'
      write(nchn1,*) 'ADJUST=',adjust
      write(nchn1,'(''ntime='',i7)') ntime
      write(nchn1,'(''dt='',f10.1,''s'')') dt
      write(nchn1,*) 'nx=',nx,' ny=',ny,' ns=',ns
      write(nchn1,'(''nx='',i3,'' ny='',i3,'' ns='',i3)') nx,ny,ns
      write(nchn1,1) xl/1000.,yl/1000.,ptop/100.
      write(nchn1,*) 'dx=',dx,' m;  dy=',dy,' m'
      ix=int(xmount/xl*nx)+1
      iy=int(ymount/yl*ny)+1
c      if(verbose.ge.2) write(*,*) ix,iy
      write(nchn1,7) (phis(ix,iy,ns-1)-phis(ix,iy,ns))/g
     :   ,(phis(ix,iy,0)-phis(ix,iy,1))/g
      write(nchn1,6) ds1(ns),ds1(1)
      write(nchn1,*) 'Mountain type:',iomtyp
      write(nchn1,9) hmount/1000.,hxwid/1000.,hywid/1000.,hexp
     :   ,xmount/1000.,ymount/1000.
      write(nchn1,'(''Diffusion: iodif='',i2)') iodif
      write(nchn1,*) ' difcof=',difcof,' difl=',difl,' rikey=',rikey
     :   ,' uvdif=',uvdif,' tsdif=',tsdif
      write(nchn1,*) 'absortion layer:'
     :   ,' raylei=',raylei,' olddra=',olddra
      write(nchn1,*) ' zm=',zm,' taum=',taum,' dragm=',dragm
     :   ,' taut=',taut,' dragt=',dragt
      write(nchn1,*) 'Latitude=',xlatit,' fcor=',fcor
      write(nchn1,*) 'ug=',ug,' dugdz=',dugdz,' vg=',vg,' dvgdz=',dvgdz
      write(nchn1,*) 'water vapor:',' qif=',qif,' qdrag=',qdrag
      write(nchn1,*) 'water:',' ifqc=',ifqc,' ifqr=',ifqr,'ifqi=',ifqi                        DC,10.2009
      write(nchn1,*) 'other options:'
      write(nchn1,*) 'tfct=',tfct,' hydro=',hydro,' ivarhs=',ivarhs
      write(nchn1,*) 'dohsmo=',dohsmo,' hdampm=',hdampm
     :   ,' hdampt=',hdampt,' perdam=',perdam,' numsmo=',numsmo
      write(nchn1,*) 'dovsmo=',dovsmo,' vdamp=',vdamp
      write(nchn1,*) 'dovsmt=',dovsmt,' vdampt=',vdampt
      write(nchn1,*) 'iobppx=',iobppx,' iobppy=',iobppy
     :              ,' iobptx=',iobptx,' iobpty=',iobpty
      write(nchn1,*) 'iobuux=',iobuux,' iobuuy=',iobuuy
     :              ,' iobvvx=',iobvvx,' iobvvy=',iobvvy
      write(nchn1,*) 'ioppuv=',ioppuv,' iophis=',iophis
     :   ,' iowlbc=',iowlbc
     :   ,' iobphy=',iobphy,' iodify=',iodify,' inipt=',inipt
      write(nchn1,*) 'time filter=',filt,ppfilt
      write(nchn1,13) eps
      write(nchn1,*) '>'
      if(grdout) then

        if(iomtyp.eq.4) then
          xmountc=-xminutm-1.5*dx
          ymountc=-yminutm-1.5*dx
        else
          xmountc=xmount
          ymountc=ymount
        endif

        call wribln(hsuf,nx,ny,dx,dy,xmountc,ymountc,fngrd)

      endif
c
      if(refout) then
        write( nchn1,*) 'reference state'
        tkoff=0.
        call wri2ar(hsuf,1,nx1,1,ny1,'hsuf    ',tkoff,nx,ny)
        tkoff=1.e5
        call wri2ar(psufi,1,nx1,1,ny1,'psufi   ',tkoff,nx,ny)
        tkoff=273.
        call wri2ar(tsufi,1,nx1,1,ny1,'tsufi   ',tkoff,nx,ny)
c        call wri3ar(tems,1,nx1,1,ny1,0,ns,'tems    ',tkoff,iomod,nx,ny)
c        call wri3ar(pts,1,nx1,1,ny1,0,ns,'pts     ',tkoff,iomod,nx,ny)
c        tkoff=0.
c        call wri3ar(phis,1,nx1,1,ny1,1,ns,'phis    ',tkoff,iomod,nx,ny)
      endif

      if(iniout) then
         write(nchn1,'(a/)') ' The initial variable fields'
      endif

      if(restar) then
        nstep1=nstep+1
        l=2
        dtl=l*dt
      else
        nstep=0
        nstep1=1
        l=1
        dtl=l*dt
       
        if(verbose.ge.1) write(*,*) 'refs'
        call refs
      endif
c
      if(.not.adjust) then
        if(verbose.ge.1) write(*,*) 'wsigbc'
        call wsigbc
      if(verbose.ge.3) write(*,'(''if:'',10i4)') (ifax(k19),k19=1,10)
      endif
    
      do iy=1,ny1
        do ix=1,nx1
         deltaz(ix,iy)=(phis(ix,iy,ns-1)+phi(ix,iy,ns-1))/g-hsuf(ix,iy)
c          deltaz(ix,iy)=0.5*(phis(ix,iy,ns-1)-phis(ix,iy,ns)
c     :     +phi(ix,iy,ns-1)-phi(ix,iy,ns))/g
c          write(99,*) 'aquideltaz',deltaz(ix,iy)
c          write(99,*) 'phi(ix,iy,ns-1)',phi(ix,iy,ns-1)
        enddo
      enddo
c
      if(iodif.ne.0) then
        if(dif3d) then
          if(verbose.ge.1) write(*,*) 'diffu'
!          call diffu(fngrd
!     :      ,wk1,wk1,wk1
!     :      ,wk2,wk2,wk2
!     :      ,wk3,wk3,wk3
!     :      ,wk4
!     :      ,wk5
!     :      ,wk6
!     :      ,wk7
!     :      ,wk8
!     :      ,wk9,wk9)

        elseif(diftx) then
          if(verbose.ge.1) write(*,*) 'diffuvt'
!          call diffuvt(fngrd
!     :      ,wk1,wk1,wk1
!     :      ,wk2,wk2,wk2
!     :      ,wk3,wk3,wk3
!     :      ,wk4
!     :      ,wk5
!     :      ,wk6
!     :      ,wk7
!     :      ,wk8
!     :      ,wk9,wk9)

	  elseif (difloc) then 
	    call diffu_local(fngrd
     :      ,wk1,wk1,wk1
     :      ,wk2,wk2,wk2
     :      ,wk3,wk3,wk3
     :      ,wk4
     :      ,wk5
     :      ,wk6
     :      ,wk7
     :      ,wk8
     :      ,wk9,wk9)
	  elseif (nonloc_tm86) then
	   call nonlocal_diffu(fngrd
     :      ,wk1,wk1,wk1
     :      ,wk2,wk2,wk2
     :      ,wk3,wk3,wk3
     :      ,wk4
     :      ,wk5
     :      ,wk6
     :      ,wk7
     :      ,wk8
     :      ,wk9,wk9)
	    elseif (nonloc_Noh) then
	   call dif_Noh(fngrd
     :      ,wk1,wk1,wk1
     :      ,wk2,wk2,wk2
     :      ,wk3,wk3,wk3
     :      ,wk4
     :      ,wk5
     :      ,wk6
     :      ,wk7
     :      ,wk8
     :      ,wk9,wk9)
	elseif (nonloc_ls96) then
!	   call diffu_LS(fngrd
!     :      ,wk1,wk1,wk1
!     :      ,wk2,wk2,wk2
!     :      ,wk3,wk3,wk3
!     :      ,wk4
!     :      ,wk5
!     :      ,wk6
!     :      ,wk7
!     :      ,wk8
!     :      ,wk9,wk9)
        else

	
          if(verbose.ge.1) write(*,*) 'diffuv'
!          call diffuv(fngrd
!     :      ,wk1,wk1,wk1
!     :      ,wk2,wk2,wk2
!     :      ,wk3,wk3,wk3
!     :      ,wk4
!     :      ,wk5
!     :      ,wk6
!     :      ,wk7
!     :      ,wk8
!     :      ,wk9,wk9)
        endif
      endif
  
c        call wri3ar(w,2,nx,2,ny,1,ns-1,'w-flux1 ',tk,iomod,nx,ny)
      if(verbose.ge.1) write(*,*) 'flux'
      if(verbose.ge.3) write(*,'(''if:'',10i4)') (ifax(k19),k19=1,10)
      call flux(wk1,wk1,wk1,wk2,wk2,wk2,wk3,wk3,wk3,wk4
     :  ,wk6,wk7,wk8)
      if(verbose.ge.3) write(*,'(''if:'',10i4)') (ifax(k19),k19=1,10)
c
      if(adjust.and.hydro) then
        write(0,*)'nh3dFatalError=Incompatible options:adjust.and.hydro'
        stop
      endif
      
      if(adjust) then
          if(verbose.ge.1) write(*,*) 'ellipt3'
          call ellipt3(0,
     :      wk1
     :      ,wk2,wk2,wk2
     :      ,wk3,wk3,wk3
     :      ,wk4,wk4,wk4
     :      ,wk5
     :      ,wk6
     :      ,wk7
     :      ,wk8)
      elseif(.not.hydro) then
        if(verbose.ge.1) write(*,*) 'ellipt'
        write(0,*)'767 ELLIPT'
        call ellipt(
     :      wk1
     :      ,wk2,wk2,wk2
     :      ,wk3,wk3,wk3
     :      ,wk4,wk4,wk4
     :      ,wk5,wk5
     :      ,wk6
     :      ,wk7
     :      ,wk8)
        if(verbose.ge.2) write(*,*) 'ellipt done'
      if(verbose.ge.3) write(*,'(''if:'',10i4)') (ifax(k19),k19=1,10)
      else
        call hydros
      endif
      
c
c soil parameters initialization
c
      if(ifsoil.ne.0) then
        call solonoii(fnmap)
      elseif(cdcoef.gt.0.) then
        call inisoil(tsmed,qvsoil,nx,ny,fnmap)
	elseif (tfix.ne.0) then
	  call inisurf(fnmap)                                                     DC,06.2010
      endif

      if(verbose.ge.3) write(*,'(''if:'',10i4)') (ifax(k19),k19=1,10)
      if(iniout.and.forout) call out6(iomod,2,fnout,fngrd)
      if(iniunf.and.unfout) call out2(2,fnout)
c
      if(restar) then
        call update
      endif
c
c time integration
c
      if(verbose.ge.0) write(*,*) 'Begin integration',starttime
c
c keep initial surface pressure for the computation of psuf-anomaly (px)
c
      pini=psuf
c
      call cpu_time(t2_cpu)                                                     VS,12.08.2007
	t_exec(1) = t2_cpu - time_begin                                           VS,12.08.2007

      do 9000 nstep=nstep1,ntime
	call cpu_time(t1_cpu)                                                     VS,12.08.20

      call timestamp(nstep,dt,starttime,timestring)
      read(timestring,'(4x,5i2)') imonth,iday,ihour,iminu,iseco                 VS,05.2007
      timeho=ihour+float(iminu)/60.+float(iseco)/3600.
      if(verbose.ge.2) write(*,*) 'nstep=',nstep,timestring,timeho


      instep=nstep
      l=1
      if( nstep.ne.1) l=2
      dtl=l*dt

      if (uvpert) then
        if (nstep.eq.nperturb) call iniperturb
        if (nstep.gt.nperturb.and.nstep.lt.nperturb+pnudg) then
          call pert_nudging
        endif
      endif
c
c verify the output status:

c
      if(forout) then
        if(mod(nstep,iotfor).eq.0 .or. nstep.eq.ntime
     :    .or. (nstep.ge.iostart .and. nstep.le. ioend
     :    .and. mod(nstep-iostart,iodelta).eq.0) ) then
          do i=1,14
             if(timestring(i:i).eq." ")timestring(i:i)="0"
          enddo
          prt=.true.
          write(nchn1,'('' Time:'',i6,2(i3,'':''),i3,1x,a14)')
     :      nstep,ihour,iminu,iseco,timestring
          do i=1,14
             if(timestring(i:i).eq." ")timestring(i:i)="0"
          enddo
          write(0,"(a9,i10,1x,a14)") "nh3dtime=",nstep,timestring
        else
          prt=.false.
          if(mod(nstep,10).eq.0) then
           do i=1,14
             if(timestring(i:i).eq." ")timestring(i:i)="0"
           enddo
           write(0,"(a9,i10,1x,a14)") "nh3dtime=",nstep,timestring
          endif
        endif
      else
        if(mod(nstep,10).eq.0 .or.nstep.eq.ntime) then
          do i=1,14
             if(timestring(i:i).eq." ")timestring(i:i)="0"
          enddo
          write(0,"(a9,i10,1x,a14)") "nh3dtime=",nstep,timestring
        endif
        prt=.false.
      endif

      if(unfout) then
        if(mod(nstep,iotunf).eq.0 .or. nstep.eq.ntime
     :    .or. (nstep.ge.iostart .and. nstep.le. ioend
     :    .and. mod(nstep-iostart,iodelta).eq.0) ) then
          unfprt=.true.
        else
          unfprt=.false.
        endif
      else
        unfprt=.false.
      endif

	call cpu_time(t2_cpu)                                                     VS,12.08.2007
	t_exec(2) = t_exec(2) + t2_cpu - t1_cpu                                   VS,12.08.2007
c
c update variable forcing (boundary conditions)

      if(ifuprefs.ne.0) then
        istepp=nstep-nstep1
        if(verbose.ge.2) write(*,*) 'call uprefst'
        call uprefst(istepp)
        call upuvt
      endif

	call cpu_time(t1_cpu)                                                     VS,12.08.2007
	t_exec(3) = t_exec(3) + t1_cpu - t2_cpu                                   VS,12.08.2007
c
c soil balance equations
c
      if(ifsoil.ne.0) then
        if(verbose.ge.1) write(*,*) 'call soil'
        if(verbose.ge.3) write(*,'(''if:'',10i4)') (ifax(k19),k19=1,10)
        call soil(fnout,timeho)
        if(verbose.ge.2) write(*,*) 'soil done'
        if(verbose.ge.3) write(*,'(''if:'',10i4)') (ifax(k19),k19=1,10)
      elseif(ifhle.eq.1) then
        if(verbose.ge.3) write(*,*) 'ifhle activated'
        le=0
        h=50
        do ix=nx/2,nx+1
          h(ix,:)=-50
        enddo
        if(verbose.ge.3) write(*,*) 'ifhle done'
	elseif(tfix.eq.1) then                                                    DC,06.2010
	 call surf_layer_t                                                       DC,06.2010
      endif
      call cpu_time(t2_cpu)                                                     VS,12.08.2007
	t_exec(4) = t_exec(4) + t2_cpu - t1_cpu                                   VS,12.08.2007
c
c Time integration of momentum equations, including boundary cond.
c
c   moment1 includes vertical winds w and wsig
c
c 
      if(adjust) then
        if(verbose.ge.1) write(*,*) 'moment1'
        call moment1(wk1,wk2,wk3,wk4,wk5,wk6,wk7,wk8)
c        call wri3ar(w,2,nx,2,ny,1,ns-1,'w-mom1  ',tk,iomod,nx,ny)
      else
        if(verbose.ge.1) write(*,*) 'moment'
        call moment(wk1,wk2,wk3,wk4,wk5,wk6,wk7,wk8)
      endif
      if(verbose.ge.3) write(*,'(''if:'',10i4)') (ifax(k19),k19=1,10)

	call cpu_time(t1_cpu)                                                     VS,12.08.2007
	t_exec(5) = t_exec(5) + t1_cpu - t2_cpu                                   VS,12.08.2007
c
c Note wk3<=>wm2 needed until aselin(w)
c Note wk4<=>ptm2 needed until thermo
       
      if(.not.adjust) then
c
c compute surface pressure and its tendency
c
        if(verbose.ge.1) write(*,*) 'sufpp'
      if(verbose.ge.3) write(*,'(''if:'',10i4)') (ifax(k19),k19=1,10)
        call sufpp
c
c re-evaluate reference state variables and compute phi at surface
c from new surface pressure
c
        if(verbose.ge.1) write(*,*) 'refs'
        call refs
c
c calculates wsig and w on boundaries:
c
        if(verbose.ge.1) write(*,*) 'wsigbc'
        call wsigbc
c
        if(verbose.ge.2) write(*,*) 'asselin w'
        call aselin(w(0,0,0,3),w(0,0,0,2),wk3
     :    ,0,nx1,0,ny1,0,ns1,filt)
c
c deallocate(wm2)--wk3 is available
c
      endif

	call cpu_time(t2_cpu)                                                     VS,12.08.2007
	t_exec(6) = t_exec(6) + t2_cpu - t1_cpu                                   VS,12.08.2007
c
c integration of potential temperature equation
c followed by imposed forcing.
c     	
      if(verbose.ge.1) write(*,*) 'thermo'
      call thermo(timeho,wk4,wk1,fngrd)
      call thermo_part2(timeho,wk4,wk1)
	call cpu_time(t1_cpu)                                                     VS,12.08.2007
	t_exec(7) = t_exec(7) + t1_cpu - t2_cpu                                   VS,12.08.2007
c
c wk4 is available
c
c Lateral absorption column
c
      if(nxsponge1.gt.0 .or. nxsponge2.gt.0) then
        do is=0,ns
          do iy=0,ny1
            do ix=0,nxsponge2
              pt(ix,iy,is,3)=pt(ix,iy,is,3)-pt(ix,iy,is,2)
     :          *dtl/tauspo(ix)
            enddo
            do ix=nx1-nxsponge1,nx1
              pt(ix,iy,is,3)=pt(ix,iy,is,3)-pt(ix,iy,is,2)
     :          *dtl/tauspo(ix)
            enddo
          enddo
        enddo
      endif

      if(nysponge1.gt.0 .or. nysponge2.gt.0) then
        do is=0,ns
          do ix=0,nx1
            do iy=0,nysponge2
              pt(ix,iy,is,3)=pt(ix,iy,is,3)-(pt(ix,iy,is,2))    !-ini_pts(is))
     :          *dtl/tauspoy(iy)
            enddo
            do iy=ny1-nysponge1,ny1
              pt(ix,iy,is,3)=pt(ix,iy,is,3)-(pt(ix,iy,is,2))      !-ini_pts(is))
     :          *dtl/tauspoy(iy)
            enddo
          enddo
        enddo
      endif

	call cpu_time(t2_cpu)                                                     VS,12.08.2007
	t_exec(8) = t_exec(8) + t2_cpu - t1_cpu                                   VS,12.08.2007
c
c thermal forcing
c
      call tforci

	call cpu_time(t1_cpu)                                                     VS,12.08.2007
	t_exec(9) = t_exec(9) + t1_cpu - t2_cpu                                   VS,12.08.2007
c
c integrate water vapor equation:
c
      if(qif.ne.0.) then
        if(verbose.ge.1) write(*,*) 'humid'
        call humid(wk1)
        call humid_part2(wk1)
        
        
	  if (impur.NE.0) then                                                    DM,05.2007
          if (.NOT.wk_impur_allocated) then                                     DM,05.2007 
            allocate(wk_impur(0:nx1,0:ny1,0:ns1))                               DM,05.2007
            wk_impur_allocated = .TRUE.                                         DM,05.2007
          endif                                                                 DM,05.2007
          call humid2(wk_impur)                                                 DM,05.2007
        endif                                                                   DM,05.2007
c
c clip
c
        do k=1,3
          do is=0,ns
            do iy=0,ny1
              do ix=0,nx1
               qv(ix,iy,is,k)=max(0.d0,qv(ix,iy,is,k))
	         qs_10(ix,iy,is,k)=max(0.d0,qs_10(ix,iy,is,k))                    DM,05.2007
              enddo
            enddo
          enddo
        enddo
      endif
      
       if(qif.ne.0) then
      if(nysponge2.gt.0.or.nysponge1.gt.0) then
        do is=1,ns-1
          do ix=1,nx1-1
            do iy=1,nysponge2
            qv(ix,iy,is,3)=qv(ix,iy,is,3)-(qv(ix,iy,is,2)-qvs(ix,iy,is))
     :          *dtl/tauspoy(iy)
            enddo
            do iy=ny1-nysponge1,ny1
            qv(ix,iy,is,3)=qv(ix,iy,is,3)-(qv(ix,iy,is,3)-qvs(ix,iy,is))
     :          *dtl/tauspoy(iy)
             !write(0,*) qvs(ix,iy,is), qv(ix,iy,is,3), iy
            enddo
          enddo
        enddo
      endif
       endif
	call cpu_time(t2_cpu)                                                     VS,12.08.2007
	t_exec(10) = t_exec(10) + t2_cpu - t1_cpu                                 VS,12.08.2007
c
c integrate cloud water equation:
c
      if(ifqc.ne.0) then
        if(verbose.ge.1) write(*,*) 'cloudw'
        call cloudw(wk1)
        call cloudw_part2(wk1)
c
c clip
c
        do k=1,3
          do is=0,ns
            do iy=0,ny1
              do ix=0,nx1
                qc(ix,iy,is,k)=max(0.d0,qc(ix,iy,is,k))
              enddo
            enddo
          enddo
        enddo
      endif

	call cpu_time(t1_cpu)                                                     VS,12.08.2007
	t_exec(11) = t_exec(11) + t1_cpu - t2_cpu                                 VS,12.08.2007

	if(ifqi.ne.0) then                                                        DC,11.2009
	 if(verbose.ge.1) write(*,*) 'cloudice'                                   DC,11.2009
        call cloudice(wk1)                                                      DC,11.2009
	  call snow(wk1)                                                          DC,11.2009
	do k=1,3                                                                  DC,11.2009
          do is=0,ns                                                            DC,11.2009
            do iy=0,ny1                                                         DC,11.2009
              do ix=0,nx1                                                       DC,11.2009
                qci(ix,iy,is,k)=max(0.d0,qci(ix,iy,is,k))                       DC,11.2009
	          qsn(ix,iy,is,k)=max(0.d0,qsn(ix,iy,is,k))                       DC,11.2009
              enddo                                                             DC,11.2009
            enddo                                                               DC,11.2209
          enddo                                                                 DC,11.2009
        enddo                                                                   DC,11.2209
      endif                                                                     DC,11.2009

c
c integrate rain water equation:
c
      if(ifqr.ne.0) then
        if(verbose.ge.1) write(*,*) 'rainw'
        call rainw(wk1)
c
c clip
c
        do k=1,3
          do is=0,ns
            do iy=0,ny1
              do ix=0,nx1
                qr(ix,iy,is,k)=max(0.d0,qr(ix,iy,is,k))
              enddo
            enddo
          enddo
        enddo
        pl=prec
      endif

	call cpu_time(t2_cpu)                                                     VS,12.08.2007
	t_exec(12) = t_exec(12) + t2_cpu - t1_cpu                                 VS,12.08.2007
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     equation of state
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        call eq_state
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     microphysical calls
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (ifqc.ne.0.and.ifqi.eq.0.and.ifqr.eq.0) then
        if(verbose.ge.1) write(*,*) 'condens'
        call condens
      !  call microphys_update
	endif     
c
	if (ifqr.ne.0.and.ifqi.eq.0) then
        if(verbose.ge.1) write(*,*) 'convert-rain'
        call microphys_update
        call condens
        call convert
        call collect
        call evaporate
        call rain
      endif

	if (ifqi.ne.0) then
        call saturation 
        call convert
       call collect
       call evaporate
       call rain                                                                DC,11.2009
	 call iceini                                                              DC,11.2009
	 call homofreeze                                                          DC,11.2009
	 call icemelt                                                             DC,11.2009
	 call accretcwsn                                                          DC,11.2009
	 call accretcisn                                                          DC,11.2009
	 call saggregation                                                        DC,11.2009
	 call accret3comp                                                         DC,11.2009
	 call bergeron                                                            DC,11.2009
	 call accretsr                                                            DC,11.2009
	 call snowmelt                                                            DC,11.2009
	 call snowprec
	 call microphys_update
        call saturation 
	 call adjustment
	endif

!	call thermo_part2(timeho,wk4,wk1)
	if (qif.ne.0)	then 
!	call humid_part2(wk1)
      do k=1,3
          do is=0,ns
            do iy=0,ny1
              do ix=0,nx1
               qv(ix,iy,is,k)=max(0.d0,qv(ix,iy,is,k))
	         enddo
            enddo
          enddo
        enddo
      endif
	if (ifqc.ne.0)  then
!	call cloudw_part2(wk1)
	do k=1,3
          do is=0,ns
            do iy=0,ny1
              do ix=0,nx1
                qc(ix,iy,is,k)=max(0.d0,qc(ix,iy,is,k))
              enddo
            enddo
          enddo
        enddo
      endif
	if (ifqi.ne.0) then
	  call cloudice_part2(wk1)
        call snow_part2(wk1)
	do k=1,3                                                                  DC,11.2009
          do is=0,ns                                                            DC,11.2009
            do iy=0,ny1                                                         DC,11.2009
              do ix=0,nx1                                                       DC,11.2009
                qci(ix,iy,is,k)=max(0.d0,qci(ix,iy,is,k))                       DC,11.2009
	          qsn(ix,iy,is,k)=max(0.d0,qsn(ix,iy,is,k))                       DC,11.2009
              enddo                                                             DC,11.2009
            enddo                                                               DC,11.2209
          enddo                                                                 DC,11.2009
        enddo     
	  endif        
	if (ifqr.ne.0) then
	  call rainw_part2(wk1)
	do k=1,3
          do is=0,ns
            do iy=0,ny1
              do ix=0,nx1
                qr(ix,iy,is,k)=max(0.d0,qr(ix,iy,is,k))
              enddo
            enddo
          enddo
        enddo
        pl=prec
	endif 
c------------------------------------------------c
c   BOUNDARY LAYER DEPTH diagnostic
c------------------------------------------------c
	call bl_depth(fngrd,nonloc_tm86,difloc)
c
c update diffusion terms
c
      if(.not.raylei .or. iodif.ne.0) then
        if(mod(nstep,nupdif).eq.0) then
          if(dif3d) then
            if(verbose.ge.1) write(*,*) 'diffu'
!            call diffu(fngrd
!     :        ,wk1,wk1,wk1
!     :        ,wk2,wk2,wk2
!     :        ,wk3,wk3,wk3
!     :        ,wk4
!     :        ,wk5
!     :        ,wk6
!     :        ,wk7
!     :        ,wk8
!     :        ,wk9,wk9)
          elseif(diftx) then
            if(verbose.ge.1) write(*,*) 'diffuvt'
!            call diffuvt(fngrd
!     :        ,wk1,wk1,wk1
!     :        ,wk2,wk2,wk2
!     :        ,wk3,wk3,wk3
!     :        ,wk4
!     :        ,wk5
!     :        ,wk6
!     :        ,wk7
!     :        ,wk8
!     :        ,wk9,wk9)
	    elseif (difloc) then 
	    call diffu_local(fngrd
     :      ,wk1,wk1,wk1
     :      ,wk2,wk2,wk2
     :      ,wk3,wk3,wk3
     :      ,wk4
     :      ,wk5
     :      ,wk6
     :      ,wk7
     :      ,wk8
     :      ,wk9,wk9)
	    elseif (nonloc_tm86) then
	   call nonlocal_diffu(fngrd
     :      ,wk1,wk1,wk1
     :      ,wk2,wk2,wk2
     :      ,wk3,wk3,wk3
     :      ,wk4
     :      ,wk5
     :      ,wk6
     :      ,wk7
     :      ,wk8
     :      ,wk9,wk9)
	    elseif (nonloc_Noh) then
	   call dif_Noh(fngrd
     :      ,wk1,wk1,wk1
     :      ,wk2,wk2,wk2
     :      ,wk3,wk3,wk3
     :      ,wk4
     :      ,wk5
     :      ,wk6
     :      ,wk7
     :      ,wk8
     :      ,wk9,wk9)
	elseif (nonloc_ls96) then
	   call diffu_LS(fngrd
     :      ,wk1,wk1,wk1
     :      ,wk2,wk2,wk2
     :      ,wk3,wk3,wk3
     :      ,wk4
     :      ,wk5
     :      ,wk6
     :      ,wk7
     :      ,wk8
     :      ,wk9,wk9)
          else
            if(verbose.ge.1) write(*,*) 'diffuv'
!            call diffuv(fngrd
!     :        ,wk1,wk1,wk1
!     :        ,wk2,wk2,wk2
!     :        ,wk3,wk3,wk3
!     :        ,wk4
!     :        ,wk5
!     :        ,wk6
!     :        ,wk7
!     :        ,wk8
!     :        ,wk9,wk9)
          endif
        endif
      endif

	call cpu_time(t1_cpu)                                                     VS,12.08.2007
	t_exec(13) = t_exec(13) + t1_cpu - t2_cpu                                 VS,12.08.2007
c
c update fluxes
c
      if(verbose.ge.1) write(*,*) 'flux'
c        call wri3ar(w,2,nx,2,ny,1,ns-1,'w-flux  ',tk,iomod,nx,ny)
      call flux(wk1,wk1,wk1,wk2,wk2,wk2,wk3,wk3,wk3,wk4
     :  ,wk6,wk7,wk8)

	call cpu_time(t2_cpu)                                                     VS,12.08.2007
	t_exec(14) = t_exec(14) + t2_cpu - t1_cpu                                 VS,12.08.2007
c
c
      
	call cpu_time(t1_cpu)                                                     VS,12.08.2007
	t_exec(15) = t_exec(15) + t1_cpu - t2_cpu                                 VS,12.08.2007
c
c
c output data for restart
c
      if(rstout .and. (mod(nstep,iotrst).eq.0 .or. nstep.eq.ntime)) then
         if(verbose.ge.1) write(*,*) 'outrst'
         call outrst(fnout)
      endif

	call cpu_time(t2_cpu)                                                     VS,12.08.2007
	t_exec(16) = t_exec(16) + t2_cpu - t1_cpu                                 VS,12.08.2007
c
c solve for hydrostatic or non-hydrostatic phi field
c
      if(adjust) then
          if(verbose.ge.1) write(*,*) 'ellipt3'
          call ellipt3(1,
     :      wk1
     :      ,wk2,wk2,wk2
     :      ,wk3,wk3,wk3
     :      ,wk4,wk4,wk4
     :      ,wk5
     :      ,wk6
     :      ,wk7
     :      ,wk8)
      elseif(.not.hydro ) then
        if(verbose.ge.1) write(*,*) 'ellipt'
        call ellipt(
     :      wk1
     :      ,wk2,wk2,wk2
     :      ,wk3,wk3,wk3
     :      ,wk4,wk4,wk4
     :      ,wk5,wk5
     :      ,wk6
     :      ,wk7
     :      ,wk8)
         if(verbose.ge.2) write(*,*) 'ellipt done'
      else
         call hydros
      endif

	call cpu_time(t1_cpu)                                                          VS,12.08.2007
      t_exec(17) = t_exec(17) + t1_cpu - t2_cpu                                 VS,12.08.2007
!      print*, 't_exec = ', t_exec(17)
!      stop
c
c calculates and prints surface drag
c
      if(drgprt .and. hmount.gt.0.) then
        if(mod(nstep,iotdrg).eq.0 .or. nstep.eq.ntime)then
          if(verbose.ge.2) write(*,*)
     :      'call dragz',ix00,ix01,ix10,ix11,iy00,iy01,iy10,iy11
          if(adjust) then
            do iy=0,ny1
            do ix=0,nx1
              pequiv(ix,iy)=psuf(ix,iy)*(1.+0.5*(phi(ix,iy,ns)
     :          +phi(ix,iy,ns-1))
     :          /(r*tsuf(ix,iy)))
            enddo
            enddo
            call drag0z(nchn1,axspan,ayspan
     :      ,nx,ny,ns,pequiv,hsuf,dx,dy,dt,nstep
     :      ,exptit,iomtyp,hnew,xmount,ymount,hxwid,hywid,hexp,pa,ptop
     :      ,ioreft,thdat,zthdat,usdat,vsdat,pts0,ivarhs
     :      ,ix00,ix01,ix10,ix11,iy00,iy01,iy10,iy11,verbose
     :      ,dragx,dragy,drag2d)
          else

            call drag0z(nchn1,axspan,ayspan
     :      ,nx,ny,ns,psuf,hsuf,dx,dy,dt,nstep
     :      ,exptit,iomtyp,hnew,xmount,ymount,hxwid,hywid,hexp,pa,ptop
     :      ,ioreft,thdat,zthdat,usdat,vsdat,pts0,ivarhs
     :      ,ix00,ix01,ix10,ix11,iy00,iy01,iy10,iy11,verbose
     :      ,dragx,dragy,drag2d)
          endif
        endif
      endif
      if(drgprt) then
        if(mod(nstep,iotduv).eq.0 .or. nstep.eq.ntime)then
          if(verbose.ge.2) write(*,*)
     :      'call draguv',ix01,ix11,iy01,iy11
          call draguv(ix01,ix11,iy01,iy11,1,ns
     :      ,nx,ny,ns,dx,dy,tems,psuf,u(0,0,0,3),v(0,0,0,3)
     :      ,w(0,0,0,3),nstep,us,vs,ptop,sigma0)
          if(domom) then
c            allocate(zed(0:nx1,0:ny1,0:ns1)) >>wk1
            allocate(iworkk(0:nx1,0:ny1,0:nlev))
            do is=0,ns
            do iy=0,ny+1
            do ix=0,nx+1
c              zed(ix,iy,is)=(phi(ix,iy,is)+phis(ix,iy,is))/g
              wk1(ix,iy,is)=(phi(ix,iy,is)+phis(ix,iy,is))/g
            enddo
            enddo
            enddo
            call momflu(ix01,ix11,iy01,iy11,1,ns,iworkk,wk1,zlevm,nlev
     :        ,nx,ny,ns,dx,dy,tems,psuf,u,v,w
     :        ,nstep,us,vs,ptop,sigma0)
            deallocate(iworkk)
c            deallocate(zed)
          endif
        endif
      endif

	call cpu_time(t2_cpu)                                                     VS,12.08.2007
      t_exec(18) = t_exec(18) + t2_cpu - t1_cpu                                 VS,12.08.2007
c
c write formated data into channel 6
c
      if(prt) then
        call out6(iomod,2,fnout,fngrd)
      endif
c
c write unformated data into channel 2

      if(unfprt) then
        call out2(2,fnout)
      endif

	call cpu_time(t1_cpu)                                                     VS,12.08.2007
      t_exec(19) = t_exec(19) + t1_cpu - t2_cpu                                 VS,12.08.2007
c
c check the stability of integration
c
      if(u(3,3,ns-1,3).gt.100.) then
        write(nchn1,*) 'integration is unstable -aborted. nstep='
     :    ,nstep
        if(.not.prt) call out6(2,2,fnout,fngrd)
        if(unfout .and. .not.unfprt) call out2(2,fnout)
        go to 9001
      endif
c
c output of data at several points
c
             if(noutu.gt.0) then
        write(52,'(i8,20(3i5,f10.2,3f7.2,f8.2))') nstep
     :  ,(ixoutu(ku),iyoutu(ku),isoutu(ku)
     :  ,psuf(ixoutu(ku),iyoutu(ku)),u(ixoutu(ku),iyoutu(ku)
     :  ,isoutu(ku),3)
     :  ,v(ixoutu(ku),iyoutu(ku),isoutu(ku),3)
     :  ,w(ixoutu(ku),iyoutu(ku),isoutu(ku),3)
     :  ,pt(ixoutu(ku),iyoutu(ku),isoutu(ku),3)
     :  +pts(ixoutu(ku),iyoutu(ku),isoutu(ku))
     :  ,ku=1,noutu)
      endif

	call cpu_time(t2_cpu)                                                     VS,12.08.2007
      t_exec(20) = t_exec(20) + t2_cpu - t1_cpu                                 VS,12.08.2007
c
c Screen output of variables
c
      if(mod(nstep,iotscrout).eq.0 .and.scrout.and. noutu.ge.1) then
        ix=ixoutu(1)
        iy=iyoutu(1)
        is=isoutu(1)

        if(drgprt) then
          write(*,'(a8,2x,a14,2x,11(a8,2X))')
     :      'nstep','Time','t_air','u','v','w','z','dragx','dragy'
     :      ,'drag2d','psuf'
          write(*,'(i8,2x,a14,2x,4(f8.3,2x),f9.1,3e10.3,f10.1)')
     :      nstep,timestring
     :      ,pt(ix,iy,is,3)+pts(ix,iy,is)
     :      ,u(ix,iy,is,3),v(ix,iy,is,3),w(ix,iy,is,3)
     :      ,(phis(ix,iy,is)+phi(ix,iy,is))/g,dragx,dragy
     :      ,drag2d,psuf(ix,iy)
        else
          write(*,'(8(a9,2X))')
     :      'nstep','Time','t_air','u','v','w','z','ps'
          write(*,'(i8,2x,a14,2x,4(f8.3,2x),f9.1,f10.1)')
     :      nstep,timestring
     :      ,pt(ix,iy,is,3)+pts(ix,iy,is)
     :      ,u(ix,iy,is,3),v(ix,iy,is,3),w(ix,iy,is,3)
     :      ,(phis(ix,iy,is)+phi(ix,iy,is))/g,psuf(ix,iy)
        endif
        if(ifsoil.ne.0) then
          write(*,'(10(a8,2X))')
     :      ' tsurf','t2','wg','w2','wr','rn','h','le','g','q'
          write(*,'(1x,4(f8.2,2X),f8.5,2X,4(f7.1,2X),f9.4)')
     :      (1.-xlake(ix,iy))*tsnoi(ix,iy)+xlake(ix,iy)*tslake(ix,iy)
     :      ,t2noi(ix,iy),wgnoi(ix,iy),w2noi(ix,iy)
     :      ,wrnoi(ix,iy),rn(ix,iy),h(ix,iy),le(ix,iy),gsolo(ix,iy)
     :      ,qv(ix,iy,is,3)
        endif
        if(ifqr.eq.1) then
          p=sigma0(is)*pp(ix,iy,3)+ptop
          t=(pt(ix,iy,is,3)+pts(ix,iy,is))*(p/p00)**akapa
          relh=qv(ix,iy,is,3)/qsat(t,p)
          write(*,'(1x,2(A8,2X),4(A12,2X))')
     :      'qv','rh','qc','qr','cond','evap'
          write(*,'(1x,f8.5,2x,f8.5,2x,4(f12.9,2x))')
     :      qv(ix,iy,is,3)*1000,relh,qc(ix,iy,is,3)*1000
     :      ,qr(ix,iy,is,3)*1000,cond(ix,iy,is)*1000
     :      ,evap(ix,iy,is)*1000
        endif
      endif

	call cpu_time(t1_cpu)                                                     VS,12.08.2007
      t_exec(21) = t_exec(21) + t1_cpu - t2_cpu                                 VS,12.08.2007
c
c update variable fields
c
      if(.not.adjust) then
        call update
c
        if(ivarhs.eq.1 .and. hnew.lt.hmount) then
          zup(:,:)=hsuf(:,:)
          call uphsuf(hsuf,hsufmax,hnew,hmount,nx,ny,nstep,nstepgrow)
          zup(:,:)=hsuf(:,:)-zup(:,:)
          if(verbose.ge.1) write(*,*) 'hnew:',hnew, zup(nx/2,ny/2)
          call newrefs(zup)
        elseif(ivarhs.eq.2) then
          call uphsuf2(hsuf,hsufmax,hnew,hmount,nx,ny,nstep,nstepgrow
     :      ,nstepflat,numsteps)
          if(verbose.ge.1) write(*,*) 'hnew:',hnew
c         call newrefs(zup)
        endif
      endif

	call cpu_time(t2_cpu)                                                     VS,12.08.2007
      t_exec(22) = t_exec(22) + t2_cpu - t1_cpu                                 VS,12.08.2007
c
9000  continue
c
9001  continue
      if(unfout) then
        write(30,'(i4)') -999
        close(30)
        close(10)
        close(11)
        close(12)
        close(13)
        close(14)
        close(15)
        close(16)
        close(17)
        close(18)
        close(19)
        close(20)
        close(21)
      endif
c
      if(ifsoil.ne.0) then
        xmin=(2-1.5)*dx-xmount
        ymin=(2-1.5)*dy-ymount
        xmax=(nx-1.5)*dx-xmount
        ymax=(ny-1.5)*dy-ymount
        call wrigrd(temmin,xmin,xmax,ymin,ymax,ttt1,ttt2
     :     ,0,nx1,0,ny1,2,nx,2,ny,fext(fnout,'tmi'))
        call wrigrd(temmax,xmin,xmax,ymin,ymax,ttt1,ttt2
     :     ,0,nx1,0,ny1,2,nx,2,ny,fext(fnout,'tma'))
        call wrigrd(evapor,xmin,xmax,ymin,ymax,ttt1,ttt2
     :     ,0,nx1,0,ny1,2,nx,2,ny,fext(fnout,'eva'))
      endif
      if(ifqr.eq.1) then
        call wrigrd(rainacc,xmin,xmax,ymin,ymax,ttt1,ttt2
     :     ,0,nx1,0,ny1,2,nx,2,ny,fext(fnout,'rai'))
      endif

      call cpu_time(time_end)
      print*, "Time of execution of code segments:"                             VS,12.08.2007
      do i=1,22                                                                 VS,12.08.2007
       print*, i, t_exec(i), '  sec'                                            VS,12.08.2007
      enddo                                                                     VS,12.08.2007

      write(*,*) 'Total cpu time:',time_end-time_begin,' sec'
      stop
1     format(' xl=',f9.2,' km',' yl=',f9.2,' km',' ptop=', f9.2,' mb ')
6     format( ' ds range from  ',f9.5,' to ',f9.5)
7     format( ' dz range from ',f9.1,'  to ',f9.1,' m')
9     format(' hmount=',f9.4,' km ','hxwid=',f6.1,' km '
     :   ,'hywid=',f6.1,' km ',' exp=',f4.1,
     :   /,' xmount=',f6.1,' km  ','ymount=',f6.1,' km ')
13    format(' the relative error in phi during iteration is required
     :to be <=' ,e7.1///)
      end

      subroutine allocvar(nonloc_tm86,nonloc_ls96)

c
c allocates arrays
c

      use alloc
      implicit real*8 (a-h,o-z)
      logical nonloc_ls96,nonloc_tm86

c      allocate(ptforc(0:nx1,0:ns))
c
      allocate(wk1(0:nx1,0:ny1,0:ns1))
      allocate(wk2(0:nx1,0:ny1,0:ns1))
      allocate(wk3(0:nx1,0:ny1,0:ns1))
      allocate(wk4(0:nx1,0:ny1,0:ns1))
      allocate(wk5(0:nx1,0:ny1,0:ns1))
      allocate(wk6(0:nx1,0:ny1,0:ns1))
      allocate(wk7(0:nx1,0:ny1,0:ns1))
      allocate(wk8(0:nx1,0:ny1,0:ns1))
      allocate(wk9(0:nx1,0:ny1,0:ns1))

      allocate(u(0:nx1,0:ny1,0:ns1,2:3))
      allocate(v(0:nx1,0:ny1,0:ns1,2:3))
      allocate(w(0:nx1,0:ny1,0:ns1,2:3))
	allocate(thgrad(0:nx1,0:ny1,0:ns1))
	thgrad=0.
      allocate(wsig(0:nx1,0:ny1,0:ns1,2:3))
      allocate(pt(0:nx1,0:ny1,0:ns1,2:3))

      if (radpar == 3) then                                                     VS,05.2007
       allocate(Radheat(0:nx1,0:ny1,0:ns1))                                     VS,05.2007
       allocate(Srad(0:nx1,0:ny1,0:ns1))                                        VS,05.2007
       allocate(Lrad(0:nx1,0:ny1,0:ns1))                                        VS,05.2007 
       allocate(Srad_surf(0:nx1,0:ny1))                                         VS,05.2007
       allocate(Lrad_surf(0:nx1,0:ny1))                                         VS,05.2007
      endif                                                                     VS,05.2007

	allocate(tsurf(0:nx1,0:ny1))                                              VS,06.2007
      allocate(ust_s(0:nx1,0:ny1))                                              DC,06.2010
	allocate(tst_s(0:nx1,0:ny1))                                              DC,06.2010
	allocate(qst_s(0:nx1,0:ny1))                                              DC,06.2010
      allocate(imp(0:nx1,0:ny1))
      ust_s=0.
      tst_s=0.
      allocate(us(0:nx1,0:ny1,0:ns1))
      allocate(vs(0:nx1,0:ny1,0:ns1))
      allocate(upert(0:nx1,0:ny1,0:ns1))
      allocate(vpert(0:nx1,0:ny1,0:ns1))
      allocate(pts(0:nx1,0:ny1,0:ns1))
      allocate(ini_pts(0:ns1))
      allocate(thpert(0:nx1,0:ny1,0:ns1))
      allocate(tems(0:nx1,0:ny1,0:ns1))
      allocate(difunu(0:nx1,0:ny1,0:ns1))
      allocate(difunv(0:nx1,0:ny1,0:ns1))
      allocate(difunw(0:nx1,0:ny1,0:ns1))
      allocate(difunt(0:nx1,0:ny1,0:ns1))
      allocate(rho(0:nx1,0:ny1,0:ns1))
      rho=0.
      allocate(dif000_h(0:nx1,0:ny1,0:ns1))
      dif000_h=0.
      allocate(dif100_h(0:nx1,0:ny1,0:ns1))
      dif100_h=0.
      allocate(dif010_h(0:nx1,0:ny1,0:ns1))
      dif010_h=0.
      allocate(dif001_h(0:nx1,0:ny1,0:ns1))
      dif001_h=0.
      allocate(ptfd(0:nx1,0:ny1,0:ns1))
      allocate(wfd(0:nx1,0:ny1,0:ns1))
      allocate(ufd(0:nx1,0:ny1,0:ns1))
      allocate(vfd(0:nx1,0:ny1,0:ns1))
      if(nonloc_Noh.or.nonloc_ls96) then
        allocate(h3e(0:nx1,0:ny1,0:ns1))
        h3e=0.
      endif
      if (nonloc_tm86.or.nonloc_ls96.or.nonloc_Noh) then
	allocate(q3c(0:nx1,0:ny1,0:ns1))
        allocate(h3c(0:nx1,0:ny1,0:ns1))
        q3c=0.
        h3c=0.
      endif
      allocate(adv_flx(0:nx1,0:ny1,0:ns1))
      allocate(adv_fly(0:nx1,0:ny1,0:ns1))
      allocate(ptnew(0:nx1,0:ny1,0:ns1))
      adv_flx=0.
      adv_fly=0.
      ptnew=0.
      allocate(phi(0:nx1,0:ny1,0:ns1))
      allocate(phis(0:nx1,0:ny1,0:ns1))
      allocate(s(0:nx1,0:ny1,0:ns1))
      allocate(thetav(0:nx1,0:ny1,0:ns1))
      allocate(difcoft(0:nx1,0:ny1,0:ns1))

      allocate(tsuf(0:nx1,0:ny1))
      allocate(psufi(0:nx1,0:ny1))
      allocate(tsufi(0:nx1,0:ny1))
      allocate(phisui(0:nx1,0:ny1))
      allocate(deltaz(0:nx1,0:ny1))
      deltaz=0.
      allocate(ptsoil(0:nx1,0:ny1))
      allocate(uvsuf(0:nx1,0:ny1))
      allocate(qvsoil(0:nx1,0:ny1))
      allocate(phisuf(0:nx1,0:ny1))
      allocate(phig(0:nx1,0:ny1))
      allocate(hk1(0:nx1,0:ny1))
      allocate(hk2(0:nx1,0:ny1))
      allocate(hk3(0:nx1,0:ny1))
      allocate(ppgra2(0:nx1,0:ny1))
      allocate(pplap(0:nx1,0:ny1))
      allocate(gpplap(0:nx1,0:ny1))

      allocate(iclay(0:nx1,0:ny1))
      allocate(isand(0:nx1,0:ny1))
      allocate(iveg(0:nx1,0:ny1))

      allocate(psuf(0:nx1,0:ny1))
      allocate(pp(0:nx2,0:ny2,4))
      allocate(pp10(0:nx1,0:ny1,4))
      allocate(pp01(0:nx1,0:ny1,4))
      allocate(ppdx10(0:nx1,0:ny1,2:2))
      allocate(ppdx(0:nx1,0:ny1,4))
      allocate(dpdx10(0:nx1,0:ny1,3:3))
      allocate(ppdy01(0:nx1,0:ny1,2:2))
      allocate(ppdy(0:nx1,0:ny1,4))
      allocate(dpdy01(0:nx1,0:ny1,3:3))
      allocate(hsuf(0:nx1,0:ny1))
      allocate(dpp(0:nx1,0:ny1))
      allocate(hsufmax(0:nx1,0:ny1))


      allocate(dppref(0:nx1,0:ny1))
      allocate(psufref(0:nx1,0:ny1))

      allocate(ubx(2,0:ny1,0:ns1))
      allocate(vbx(2,0:ny1,0:ns1))
      allocate(uby(0:nx1,2,0:ns1))
      allocate(vby(0:nx1,2,0:ns1))
      allocate(ucc(2,2,0:ns1))
      allocate(vcc(2,2,0:ns1))
      allocate(ubx2(2,0:ny1,0:ns1))
      allocate(vbx2(2,0:ny1,0:ns1))
      allocate(uby2(0:nx1,2,0:ns1))
      allocate(vby2(0:nx1,2,0:ns1))
      allocate(ucc2(2,2,0:ns1))
      allocate(vcc2(2,2,0:ns1))
      allocate(ppbx(2,0:ny2))
      allocate(ppby(0:nx2,2))
      allocate(ppcc(2,2))
      allocate(ppbx2(2,0:ny2))
      allocate(ppby2(0:nx2,2))
      allocate(ppcc2(2,2))
      allocate(ptbx(2,0:ny1,0:ns1))
      allocate(ptby(0:nx1,2,0:ns1))
      allocate(ptcc(2,2,0:ns1))
      allocate(qvbx(2,0:ny1,0:ns1))
      allocate(qvby(0:nx1,2,0:ns1))
      allocate(qvcc(2,2,0:ns1))
      allocate(qcbx(2,0:ny1,0:ns1))
      allocate(qcby(0:nx1,2,0:ns1))
      allocate(qccc(2,2,0:ns1))
      allocate(qrbx(2,0:ny1,0:ns1))
      allocate(qrby(0:nx1,2,0:ns1))
      allocate(qrcc(2,2,0:ns1))
      allocate(pgradx(2,0:ny1,0:ns))
      allocate(pgrads(0:nx1,0:ny1,2))
      allocate(pgrady(0:nx1,2,0:ns))
    
      allocate(qsbx(2,0:ny1,0:ns1))                                             DM,05.2007
      allocate(qsby(0:nx1,2,0:ns1))                                             DM,05.2007 
      allocate(qscc(2,2,0:ns1))                                                 DM,05.2007

	allocate(qcibx(2,0:ny1,0:ns1))                                            DC,11.2009
      allocate(qciby(0:nx1,2,0:ns1))                                            DC,11.2009 
      allocate(qcicc(2,2,0:ns1))                                                DC,11.2009
	allocate(qsnbx(2,0:ny1,0:ns1))                                            DC,11.2009
      allocate(qsnby(0:nx1,2,0:ns1))                                            DC,11.2009 
      allocate(qsncc(2,2,0:ns1))                                                DC,11.2009

      allocate(endrag(0:ns))
      allocate(taudra(0:ns))
      allocate(deltazt(0:ns))

      allocate(ds0(0:ns1))
      allocate(ds02(0:ns1))
      allocate(ds04(0:ns1))
      allocate(ds02a(0:ns1))
      allocate(ds04a(0:ns1))
      allocate(ds08a(0:ns1))
      allocate(ds1(0:ns1))
      allocate(ds12(0:ns1))
      allocate(ds14(0:ns1))
      allocate(s1ds0(0:ns1))
      allocate(s0ds1(0:ns1))
      allocate(s0ds2a(0:ns1))
      allocate(s0ds4a(0:ns1))
      allocate(s1ds4a(0:ns1))
      allocate(ds1x(0:ns1))
      allocate(ds1dx(0:ns1))
      allocate(ds1y(0:ns1))
      allocate(ds1dy(0:ns1))
      allocate(dtxs(0:ns1))
      allocate(dtys(0:ns1))
      allocate(dxys(0:ns1))

      allocate(delka2(ns))

      allocate(trigsx(nx-1))
      allocate(trigsy(ny-1))
      allocate(sipcox(nx-1))
      allocate(simcox(nx-1))
      allocate(sipcoy(ny-1))
      allocate(simcoy(ny-1))
      allocate(xlamb(nx-1))
      allocate(ylamb(ny-1))
      allocate(xylamb(nx-1,ny-1))
      allocate(eigxy(2:nx,2:ny))

c      if(ifsoil.ne.0) then
        allocate(veg(0:nx1,0:ny1))
        allocate(xlai(0:nx1,0:ny1))
        allocate(rsmin(0:nx1,0:ny1))
        allocate(alb(0:nx1,0:ny1))
        allocate(emis(0:nx1,0:ny1))
        allocate(z0(0:nx1,0:ny1))
        allocate(z0h(0:nx1,0:ny1))
        allocate(rra(0:nx1,0:ny1))
        allocate(z0hw(0:nx1,0:ny1))
        allocate(z0w(0:nx1,0:ny1))
        allocate(rg(0:nx1,0:ny1))
        allocate(rat(0:nx1,0:ny1))
        allocate(xlake(0:nx1,0:ny1))
	  allocate(watice(0:nx1,0:ny1))                                           DC,07.2010
	  allocate(xsea(0:nx1,0:ny1))                                             DC,07.2010
	  allocate(tsfix(0:nx1,0:ny1))                                            DC,07.2010
	  allocate(hbl(0:nx1,0:ny1))                                              DC,07.2010
	  hbl=1.                                                                  DC,07.2010
	  allocate(dzits(0:nx1,0:ny1))                                            DC,07.2010
	  dzits=0.                                                                DC,07.2010
	  allocate(Fv(0:nx1,0:ny1))                                               DC,07.2010
	  Fv=0.                                                                   DC,07.2010
	  allocate(wm(0:nx1,0:ny1))                                               DC,07.2010
	  wm=0.                                                                   DC,07.2010
	  allocate(gammah(0:nx1,0:ny1))                                           DC,07.2010
	  gammah=0.                                                               DC,07.2010
	  allocate(gammaq(0:nx1,0:ny1))                                           DC,03.2011
	  gammaq=0.                                                               DC,03.2011
	  allocate(wstar(0:nx1,0:ny1))                                            DC,07.2010
	  wstar=0.                                                                DC,07.2010
	  allocate(Ribulk(0:nx1,0:ny1,0:ns1))                                     DC,07.2010
	  Ribulk=0.                                                                    DC,11.2009
	
        allocate(pl(0:nx1,0:ny1))
        allocate(zonirr(0:nx1,0:ny1))
        allocate(tsmed(0:nx1,0:ny1))
        allocate(tslake(0:nx1,0:ny1))
        tslake=0.
        allocate(tsnoi(0:nx1,0:ny1))
        allocate(t2noi(0:nx1,0:ny1))
        allocate(wgnoi(0:nx1,0:ny1))
        allocate(w2noi(0:nx1,0:ny1))
        allocate(wrnoi(0:nx1,0:ny1))
        allocate(tsurw(0:nx1,0:ny1))
        allocate(tmsur(0:nx1,0:ny1))
        allocate(rn(0:nx1,0:ny1))
        rn=0.
        allocate(h(0:nx1,0:ny1))
        h=0.
        allocate(le(0:nx1,0:ny1))
        le=0.
        allocate(gsolo(0:nx1,0:ny1))
        gsolo=0.
        allocate(xdd2(0:nx1,0:ny1))
        xdd2=0.
        allocate(cdm(0:nx1,0:ny1))
        cdm=0.
        allocate(rainacc(0:nx1,0:ny1))
        rainacc=0.
	  allocate(snowacc(0:nx1,0:ny1))                                          DC,11.2009
	  

!     Dima@ 110507 - Allocate buffer for impur
        allocate(le_impur(0:nx1,0:ny1))                                         DM,05.2007      
c      endif

      allocate(temmax(0:nx1,0:ny1))
      allocate(temmin(0:nx1,0:ny1))
      allocate(evapor(0:nx1,0:ny1))
      allocate(ptforc(0:nx1,0:ny1))
      allocate(phi00s(0:nx1,0:ny1))

      allocate(tauspo(0:nx1))
      allocate(tauspoy(0:ny1))
      allocate(a1pos3(ns-1))
      allocate(a3pos3(ns-1))
      allocate(a21pos3(nx-1,ny-1))
      allocate(aves(0:ns))


      allocate(duul(ny1))
      allocate(duur(ny1))
      allocate(dvvl(nx1))
      allocate(dvvr(nx1))
      allocate(phibx1(0:ny1))
      allocate(phiby1(0:nx1))
      allocate(phibx2(0:ny1))
      allocate(phiby2(0:nx1))

      allocate(dphsdx(2,0:ny1))
      allocate(dphsdy(0:nx1,2))
      allocate(hdamp(0:ns1))
      
      allocate(advfx(0:nx1,0:ny1,0:ns1))
	allocate(advfy(0:nx1,0:ny1,0:ns1))

      if(qif.ne.0.) then

c       Dima@ 13.04.07
        allocate(qs_10(0:nx1qv,0:ny1qv,0:ns1qv,1:3))                            DM,05.2007
        allocate(qs_10s(0:nx1qv,0:ny1qv,0:ns1qv))                               DM,05.2007 
        allocate(difunq_10(0:nx1qv,0:ny1qv,0:ns1qv))                            DM,05.2007

        allocate(qv(0:nx1qv,0:ny1qv,0:ns1qv,1:3))
        allocate(difunq(0:nx1qv,0:ny1qv,0:ns1qv))
        allocate(qvs(0:nx1qv,0:ny1qv,0:ns1qv))
        allocate(cond(0:nx1qv,0:ny1qv,0:ns1qv))
        allocate(evap(0:nx1qv,0:ny1qv,0:ns1qv))
        qv=0.
        qvs=0.
        cond=0.
        evap=0.
        difunq=0.
        if(ifqc.ne.0) then
          allocate(qc(0:nx1qc,0:ny1qc,0:ns1qc,1:3))
          allocate(difunqc(0:nx1qc,0:ny1qc,0:ns1qc))
          allocate(auto(0:nx1qc,0:ny1qc,0:ns1qc))
          allocate(col(0:nx1qc,0:ny1qc,0:ns1qc))

          qc=0.
          difunqc=0.
          auto=0.
          col=0.
          if(ifqr.ne.0) then
            allocate(qr(0:nx1qr,0:ny1qr,0:ns1qr,1:3))
            allocate(vrain(0:nx1qr,0:ny1qr,0:ns1qr))
            allocate(difunqr(0:nx1qr,0:ny1qr,0:ns1qr))
            allocate(prec(0:nx1,0:ny1))
            qr=0.
            difunqr=0.
            vrain=0.
            prec=0.
            rainacc=0.
          endif
	    if(ifqi.ne.0) then                                                       DC,10.2009
	allocate(qci(0:nx1qi,0:ny1qi,0:ns1qi,1:3))                                   DC,10.2009
	allocate(qsn(0:nx1qi,0:ny1qi,0:ns1qi,1:3))                                   DC,10.2009
	allocate(vdepi(0:nx1qi,0:ny1qi,0:ns1qi))                                     DC,10.2009
	allocate(vdeps(0:nx1qi,0:ny1qi,0:ns1qi))                                     DC,10.2009
	allocate(vini(0:nx1qi,0:ny1qi,0:ns1qi))                                      DC,10.2009
	allocate(imlt(0:nx1qi,0:ny1qi,0:ns1qi))                                      DC,10.2009
	allocate(hmfrz(0:nx1qi,0:ny1qi,0:ns1qi))                                     DC,10.2009
	allocate(sacrw(0:nx1qi,0:ny1qi,0:ns1qi))                                     DC,10.2009
	allocate(sacrwr(0:nx1qi,0:ny1qi,0:ns1qi))                                    DC,10.2009
	allocate(sbercw(0:nx1qi,0:ny1qi,0:ns1qi))                                    DC,10.2009
	allocate(sberci(0:nx1qi,0:ny1qi,0:ns1qi))                                    DC,10.2009
	allocate(ibercw(0:nx1qi,0:ny1qi,0:ns1qi))                                    DC,10.2009
	allocate(sagg(0:nx1qi,0:ny1qi,0:ns1qi))                                      DC,10.2009
	allocate(saci(0:nx1qi,0:ny1qi,0:ns1qi))                                      DC,11.2009
	allocate(raci(0:nx1qi,0:ny1qi,0:ns1qi))                                      DC,11.2009
	allocate(iacr(0:nx1qi,0:ny1qi,0:ns1qi))                                      DC,11.2009
	allocate(smlt(0:nx1qi,0:ny1qi,0:ns1qi))                                      DC,11.2009
      allocate(sacrr(0:nx1qi,0:ny1qi,0:ns1qi))                                     DC,11.2009
	allocate(difunqci(0:nx1qi,0:ny1qi,0:ns1qi))                                  DC,11.2009
	allocate(difunqsn(0:nx1qi,0:ny1qi,0:ns1qi))                                  DC,11.2009
	allocate(vsnow(0:nx1qi,0:ny1qi,0:ns1qi)) 
	allocate(qs(0:nx1qi,0:ny1qi,0:ns1qi))
	allocate(tendx(0:nx1,0:ny1,0:ns1))
	allocate(tendy(0:nx1,0:ny1,0:ns1))
	allocate(precsnow(0:nx1,0:ny1))                                              DC,11.2009
	
	
	qci=0.
	qsn=0.
	vdepi=0.
	vdeps=0.
	vini=0.
	imlt=0.
	hmfrz=0.
	sacrw=0.
	sacrwr=0.
	sbercw=0.
	sberci=0.
	ibercw=0.
	sagg=0.
	saci=0.
	raci=0.
	iacr=0.
	smlt=0.
	sacrr=0.
	difunqci=0.
	difunqsn=0.
	vsnow=0.
	precsnow=0.
	qs=0.
	    endif
        endif
      endif

      maxnxny=max(nx1,ny1)
      allocate(a1sphi(0:maxnxny))
      allocate(a2sphi(0:maxnxny))
      allocate(a3sphi(0:maxnxny))
      allocate(csphi(0:maxnxny))
      allocate(dsphi(0:maxnxny))

      allocate(zlevm(nlev),zlev(0:nlev))

      allocate(w1d1(0:ns1))
      allocate(w1d2(0:ns1))
      allocate(pequiv(0:nx1,0:ny1))

      allocate(panom(0:nx1,0:ny1))
      allocate(pini(0:nx1,0:ny1))

      allocate(emo(ns,ns))
      allocate(emi(ns,ns))
      allocate(vlam(ns))
      allocate(e(ns+1))
      return
      end



c**********************************************************************







      subroutine blob

      use alloc

      implicit real*8 (a-h,o-z)
c
c initial temperature perturbation (depends on experiment!)             #exp#
c
c delpt: max temperature inversion (at x=0.)
c depinv: depth of inversion
c
      amppt=3.
      ampqv=2.e-3
      centrx=0.
      centrz=1500.
      widtxr=10800
      widtxl=10800.
      widtz=2000.

      do iy=2,ny
        do ix=2,nx
          do is=0,ns-1
            z=phis(ix,iy,is)/g
            x=(ix-1.5)*dx-xmount
            if(x.ge.0.) then
              beta=(widtxr-x)/widtxr
              beta=max(0.d0,beta)
              beta=beta*(widtz-z)/widtz
            else
              beta=(x+widtxl)/widtxl
              beta=max(0.d0,beta)
              beta=beta*(widtz-z)/widtz
            endif
            beta=max(0.d0,beta)
            pt(ix,iy,is,3)=pt(ix,iy,is,3)+amppt*beta
            qv(ix,iy,is,3)=qv(ix,iy,is,3)+ampqv*beta
          enddo
        enddo
      enddo
      qv(:,:,:,1)=qv(:,:,:,3)
      qv(:,:,:,2)=qv(:,:,:,3)
      return
      end

      subroutine collect
c
c Evaluation of the collection of raindrops by cloud (cloud water->rain water)
c
c col=const*rho*(3/8)*qc*qr**(7/8)
c
c after Emanuel (1994)
c
      use alloc

      implicit real*8 (a-h,o-z)
c      include 'nh3dpa10.inc'


      do is=0,ns
        do iy=0,ny1
          do ix=0,nx1
!            p=sigma0(is)*pp(ix,iy,2)+ptop
!            rho=p/(r*(pts(ix,iy,is)
!     :     +pt(ix,iy,is,2))*(p/p00)**akapa*(1+0.61*qv(ix,iy,is,2)
!     :     -qc(ix,iy,is,2)))
            xlr=(rho(ix,iy,is)*qr(ix,iy,is,2)/(pi*rhol*xnor))**0.25
	      col(ix,iy,is)=pi*xnor*xa*qc(ix,iy,is,2)*gm38*xlr**(xb+3)
     :        *0.25*sqrt(rho0/rho(ix,iy,is))
          enddo
        enddo
      enddo
      return
      end


      subroutine condens
c
c Evaluation of the condensation/evaporation (vapour<->cloud water)
c
c cond=d qv/dt=(qsat-qv)/(dt*(1+lv**2*qsat/cp*r*t))
c
c after Yau and Austin (1979), Wilhelmson and Klemp (1978)
c
      use alloc

      implicit real*8 (a-h,o-z)

      do is=1,ns-1
        do iy=1,ny1-1
          do ix=1,nx1-1
            p=sigma0(is)*pp(ix,iy,2)+ptop
            t=(pts(ix,iy,is)+pt(ix,iy,is,2))*(p/p00)**akapa
            rasrv=287.05/461.51
            qsatur=qsat(t,p)
!            if(ix.eq.3.and.is.eq.35) 
!     :      write(0,*) is,qsatur,qv(ix,iy,is,2)
!            if(ix.eq.3.and.is.eq.40) write(0,*)
!     :       (1.+4097.93*hlat*qsatur/(cp*(t-35.86)**2)/(1.-
!     :          (1.-rasrv)*esat(t)/p))
            if(qv(ix,iy,is,2).gt.qsatur.or.qc(ix,iy,is,2).gt.0.) then
              cond(ix,iy,is)=-min(qc(ix,iy,is,2)
     :          ,(qsatur-qv(ix,iy,is,2))/
     :          (1.+4097.93*hlat*qsatur/(cp*(t-35.86)**2)/(1.-
     :          (1.-rasrv)*esat(t)/p)))/dtl
            else
              cond(ix,iy,is)=0.
            endif
          enddo
        enddo
      enddo
      return
      end

      subroutine const
c----------------------------------------------------------------------
c prepare the often-required constants
c-----------------------------------------------------------------------
c
      use alloc
      implicit real*8 (a-h,o-z)
c      include 'nh3dpa10.inc'
c
      parameter(third2=2./3.)
c
      dxx=dx*dx
      dx05=0.5*dx
      dx025=0.25*dx
      dx2=2.*dx
      dx4=4.*dx
      dx8=8.*dx
      dxdt2=dx/dt2
      dxdt=dx/dt
      dxt2=2.*dx*dt
      rdx=dx
      rdy=dy
      rdt=dt
c
      dyy=dy*dy
      dy05=0.5*dy
      dy025=0.25*dy
      dy2=2.*dy
      dy4=4.*dy
      dy8=8.*dy
      dydt2=dy/dt2
      dydt=dy/dt
      dyt2=2.*dy*dt
c
      dxy=dx*dy
      dxy025=0.25*dxy
c
      dtxy=dt*dx*dy
c
      ds1(0)=sigma1(1)-sigma1(0)
      ds0(0)=sigma0(1)-sigma0(0)
      do is=1,ns
        ds1(is)=sigma1(is+1)-sigma1(is)
        ds0(is)=sigma0(is)-sigma0(is-1)
      enddo

      ds1(ns1)=ds1(ns)
      ds0(ns1)=sigma0(ns1)-sigma0(ns)
c
      pi=4.d0*datan(1.0d0)
      pi05=0.5*pi
c
c dsng is the distance (in sigma coord.) between ground level (sigma(ns)) and
c the lower level in ssigma coordinates
c
      dsng=sigma0(ns)-sigma1(ns)
      dsr05=dsng*r05
c
      do is=0,ns1
        ds12(is)=2.*ds1(is)
        ds14(is)=4.*ds1(is)
        ds02(is)=2.*ds0(is)
        ds04(is)=4.*ds0(is)
        ds1x(is)=dx*ds1(is)
        ds1y(is)=dy*ds1(is)
        dtxs(is)=dt*ds1x(is)
        dtys(is)=dt*ds1y(is)
        dxys(is)=dx*ds1y(is)
        ds1dx(is)=ds1(is)/dx
        ds1dy(is)=ds1(is)/dy
        s0ds1(is)=sigma0(is)/ds1(is)
        s1ds0(is)=sigma1(is)/ds0(is)
      enddo
      do is=0,ns
        ds02a(is)=ds0(is+1)+ds0(is)
        s0ds2a(is)=sigma0(is)/ds02a(is)
        ds04a(is)=2.*ds02a(is)
        ds08a(is)=4.*ds02a(is)
      enddo

      do is=1,ns1
        s1ds4a(is)=sigma1(is)/(2.*(ds1(is-1)+ds1(is)))
      enddo

      do is=0,ns
        s0ds4a(is)=0.5*s0ds2a(is)
      enddo
c
c
c coriolis parameter:
c moved to readpa
c
c     fcor=2.*omega*sin(xlatit*pi/180.)
c


      return
      end


      subroutine convert
c
c Evaluation of the autoconversion term (cloud water->rain water)
c
c auto=k1(qc-qc0)
c
c after Kessler(1969)
c
      use alloc

      implicit real*8 (a-h,o-z)
c      include 'nh3dpa10.inc'


      do is=0,ns
        do iy=0,ny1
          do ix=0,nx1
            if (qc(ix,iy,is,3).gt.qco) then
              auto(ix,iy,is)=xk4*(qc(ix,iy,is,3)-qco)
            else
              auto(ix,iy,is)=0.
            endif
          enddo
        enddo
      enddo
      return
      end

      subroutine teixeira(nx0,nx1,ny0,ny1,ns0,ns1,ix,iy,is0,is1
     : ,var,dvar,difcoft,interpol,difcofm,difcofp,deltaz,dt)

      !diffuses column ix,iy of variable var with teixeira's method

      !nx0:nx1,ny0:ny1,ns0:ns1 domain size
      !ix0:ix1,iy0:iy1,is0:is1 interior points
      !call in nhad:
      !  call teixeira(0,nx+1,0,ny+1,0,ns+1,ix,iy,2,ns-1
      !    ,u(0,0,0,2),diffunu,dif100,1,w1d1,w1d2,deltaz,2.*dt)
      !var: variable to diffuse
      !deltaz: grid increment
      !dvar: change in variable due to diffusion (dt*k*d**2 var/dz**2)
      !difcoft: diffusion coefficient
      !interpol: switch to choose interpolation method
      !difcofm, difcofp: work arrays

      implicit real*8(a-h,o-z)
      real*8 var(nx0:nx1,ny0:ny1,ns0:ns1+1)
     :  ,difcoft(nx0:nx1,ny0:ny1,ns0:ns1)
     :  ,dvar(nx0:nx1,ny0:ny1,ns0:ns1)

      integer jiter, niter,interpol,km,kp, km1,kp1
      real*8 dsminus,dsplus,betap,betam,deltazm,deltazp
     : ,difcofmn, difcofmnm,difcofpnp,dsmed
      real*8 difcofm(ns0:ns1),difcofp(ns0:ns1),deltaz(ns0:ns1)
     :      ,deltaz1(ns0:ns1),difcofp1(ns0:ns1),difcofm1(ns0:ns1)

      niter=2
      gama=1.d0/6.d0!1.d0/6.d0

      do is=is0-1,is1+1
         difcofp(is)=0.5d0*(difcoft(ix,iy,is)+difcoft(ix,iy,is+1))
         difcofm(is)=0.5d0*(difcoft(ix,iy,is-1)+difcoft(ix,is,is))
c          difcofm(is)=difcoft(ix,iy,is)
c          difcofp(is)=difcoft(ix,iy,is+1)
         deltaz1(is)=-deltaz(is)
      enddo

      ismax=is1+1
      ismin=is0-1

      do jiter=1,niter
        do is=is0,is1

         dsminus=sqrt(2.d0*dt/((difcofp(is)+difcofm(is))*gama))
     :           *difcofm(is)
         dsplus=sqrt(2.d0*dt/((difcofp(is)+difcofm(is))*gama))
     :           *difcofp(is)

         km=is-1
         deltazm=deltaz1(km)
         do while(dsminus.gt.deltazm .and. km.ge.is0)
           km=km-1
           deltazm=deltaz1(km)+deltazm
         enddo
         betam=1.d0-(deltazm-dsminus)/deltaz1(km)

         kp=is
         deltazp=deltaz1(kp)
         do while(dsplus.gt.deltazp .and. kp.le.is1)
           kp=kp+1
           deltazp=deltaz1(kp)+deltazp
         enddo
         betap=1.d0-(deltazp-dsplus)/deltaz1(kp)

         difcofmnm=difcoft(ix,iy,min(ismax,km))
         difcofmn=difcoft(ix,iy,min(ismax,km+1))

         difcofpn=difcoft(ix,iy,max(ismin,kp))
         difcofpnp=difcoft(ix,iy,max(ismin,kp+1))

         difcofm(is)=(1.d0-betam)*difcofmn+betam*difcofmnm
         difcofp(is)=(1.d0-betap)*difcofpn+betap*difcofpnp

       enddo
      enddo

      do is=is0,is1
c        write(*,*) 'is 2:',is
c         dsminus=sqrt((difcofm(is)+difcofp(is))*difcofm(is)*
c     :     dt/(2.d0*difcofp(is)*gama))
c
c         dsplus=sqrt((difcofm(is)+difcofp(is))*difcofp(is)*
c     :       dt/(2.d0*difcofm(is)*gama))

         dsminus=sqrt(2.d0*dt/((difcofp(is)+difcofm(is))*gama))
     :           *difcofm(is)
         dsplus=sqrt(2.d0*dt/((difcofp(is)+difcofm(is))*gama))
     :           *difcofp(is)

         km=is-1
         deltazm=deltaz1(km)
         do while(dsminus.gt.deltazm .and. km.ge.is0)
           km=km-1
           deltazm=deltaz1(km)+deltazm
         enddo
         betam=1.d0-(deltazm-dsminus)/deltaz1(km)

         kp=is
         deltazp=deltaz1(kp)
         do while(dsplus.gt.deltazp .and. kp.le.is1)
           kp=kp+1
           deltazp=deltaz1(kp)+deltazp
         enddo
         betap=1.d0-(deltazp-dsplus)/deltaz1(kp)

         if(interpol.eq.1) then

c linear interpolation
           varminus2=var(ix,iy,min(ismax,km))!km
           varminus1=var(ix,iy,min(ismax,km+1))!is

           varplus1=var(ix,iy,max(ismin,kp))   !is
           varplus2=var(ix,iy,max(ismin,kp+1)) !kp

           varminus=(1.d0-betam)*varminus1+betam*varminus2
           varplus=(1.d0-betap)*varplus1+betap*varplus2

         elseif(interpol.eq.2) then

           gamam=1.d0-betam
           gamap=1.d0-betap

           varmin3=var(ix,iy,max(ismin,km+1))
           varmin2=var(ix,iy,max(ismin,is))
           varmin1=var(ix,iy,max(ismin,is-1))

           varplus1=var(ix,iy,min(ismax,is+1))
           varplus2=var(ix,iy,min(ismax,is))
           varplus3=var(ix,iy,min(ismax,kp-1))

           varminus=varmin1*(betam)*(-gamam)/((-2.d0*(-1.d0)))
     :  +varmin2*(1+betam)*(-gamam)/(1.d0*(-1.d0))
     :  +varmin3*(1+betam)*(betam)/(2.d0*1.d0)

           varplus=varplus1*(betap)*(-gamap)/((-2.d0*(-1.d0)))
     :   +varplus2*(1+betap)*(-gamap)/(1.d0*(-1.d0))
     :   +varplus3*(1+betap)*(betap)/(2.d0*1.d0)

c          pause 'Estou aqui'

           else
             write(0,*) 'nh3dFatalError= in interpol teix'
               stop 'error interpol'
           endif

c rate of change of dar due to diffusion (dvar=d(var)/dt)

c           dvar(ix,iy,is)=0.5d0*difcofp(is)*(1.d0/(dsplus*dsplus)
c     :  +1.d0/(dsminus*dsplus))*(varplus-var(ix,iy,is))
c     :  -0.5d0*difcofm(is)*(1.d0/(dsminus*dsminus)
c     :  +1.d0/(dsplus*dsminus))*(var(ix,iy,is)-varminus)
           dsmed=(dsminus+dsplus)/2.d0
           dvar(ix,iy,is)=((difcofp(is)*(varplus-var(ix,iy,is))/dsplus)
     :            -difcofm(is)*(var(ix,iy,is)-varminus)/dsminus)/dsmed

      enddo
      return
      end


      subroutine evaporate
c
c Evaluation of the evaporation of rainwater (rain water<->vapour)
c
c evap=a rather complex expression
c
c after Emanuel (1994) and Ogura and Takahashi (1973)
c
      use alloc

      implicit real*8(a-h,o-z)


      do is=0,ns
        do iy=0,ny1
          do ix=0,nx1
            p=sigma0(is)*pp(ix,iy,2)+ptop
            t=(pts(ix,iy,is)+pt(ix,iy,is,2))*(p/p00)**akapa
!	      rho=p/(r*t*(1+0.61*qv(ix,iy,is,2)-qc(ix,iy,is,2)))
            xlr=(rho(ix,iy,is)*qr(ix,iy,is,2)/(pi*rhol*xnor))**0.25
	      qsatur=qsat(t,p)
	      if (ifqi.eq.0) then 
            evap(ix,iy,is)=2*pi*dv*xnor*(qsatur-qv(ix,iy,is,2))/(1+
     :        dv*hlat*hlat*qsatur*rho(ix,iy,is)/(xkt*rv*t*t))*
     :        (0.78*xlr*xlr+0.31*sch**(1./3.)*gm29*sqrt(xa/xniu)*
     :        (rho0/rho(ix,iy,is))**0.25*xlr**((xb+5.)/2.))
	      else
	      evap(ix,iy,is)=2*pi*dv*xnor*(qs(ix,iy,is)-qv(ix,iy,is,3))
     :	  /(1+dv*hlat*hlat*qs(ix,iy,is)*rho(ix,iy,is)/(xkt*rv*t*t))*
     :        (0.78*xlr*xlr+0.31*sch**(1./3.)*gm29*sqrt(xa/xniu)*
     :        (rho0/rho(ix,iy,is))**0.25*xlr**((xb+5.)/2.))
	      endif
            
          enddo
        enddo
      enddo
      return
      end



      subroutine flux
     :  (ufx,vfx,wfx,ufy,vfy,wfy,ufs,vfs,wfs,wsigp
     :  ,uflux,vflux,wflux)

c
c 4 worksp
c
c wk1(ufx,vfx,wfx), wk2(ufy,vfy,wfy), wk3(ufs,vfs,wfs), wk4(wsigp)
c
c-----------------------------------------------------------------------
c calculate momentum fluxes
c-----------------------------------------------------------------------
c
      use alloc
      implicit real*8(a-h,o-z)
c
      dimension ufx(0:nx1,0:ny1,0:ns1),ufs(0:nx1,0:ny1,0:ns1)
      dimension ufy(0:nx1,0:ny1,0:ns1)
      dimension vfx(0:nx1,0:ny1,0:ns1),vfs(0:nx1,0:ny1,0:ns1)
      dimension vfy(0:nx1,0:ny1,0:ns1)
      dimension wfx(0:nx1,0:ny1,0:ns1),wfs(0:nx1,0:ny1,0:ns1)
      dimension wfy(0:nx1,0:ny1,0:ns1)
      dimension wsigp(0:nx1,0:ny1,0:ns1)
      dimension uflux(0:nx1,0:ny1,0:ns1)
      dimension vflux(0:nx1,0:ny1,0:ns1)
      dimension wflux(0:nx1,0:ny1,0:ns1)
cc
c      equivalence (wk1,ufx,vfy,wfx)
c      equivalence (wk2,ufy,vfx,wfy)
c      equivalence (wk3,ufs,vfs,wfs)
c      equivalence (wk4,wsigp)
c
c ufx - u momentum flux in x-direction
c ufy - u momentum flux in y-direction
c ufs - u momentum flux in sigma-direction
c vfx - v momentum flux in x-direction
c vfy - v momentum flux in y-direction
c vfs - v momentum flux in sigma-direction
c wfx - w momentum flux in x-direction
c wfy - w momentum flux in y-direction
c wfs - w momentum flux in sigma-direction
c
c      allocate(uflux(0:nx1,0:ny1,0:ns1))
c      allocate(vflux(0:nx1,0:ny1,0:ns1))
c      allocate(wflux(0:nx1,0:ny1,0:ns1))
c      allocate(wsigp(0:nx1,0:ny1,0:ns1))
      if(verbose.ge.3) write(*,'(''fl1:'',10i4)') (ifax(k19),k19=1,10)
      do is=0,ns1
      do iy=1,ny1
      do ix=1,nx1
        wsigp(ix,iy,is)=pp(ix,iy,3)*wsig(ix,iy,is,3)
      enddo
      enddo
      enddo

c
c update momentum fluxes
c
c      allocate(ufx(0:nx1,0:ny1,0:ns1))
c      allocate(ufy(0:nx1,0:ny1,0:ns1))
c      allocate(ufs(0:nx1,0:ny1,0:ns1))

      if(verbose.ge.3) write(*,'(''ufx:'',10i4)') (ifax(k19),k19=1,10)
      do is=0,ns
        do iy=2,ny
        do ix=1,nx1
          ufx(ix,iy,is)=(u(ix-1,iy,is,3)+u(ix,iy,is,3))
     :      *(u(ix-1,iy,is,3)*pp10(ix-1,iy,3)+u(ix,iy,is,3)
     :      *pp10(ix,iy,3))/dx4
        enddo
        enddo
      if(verbose.ge.3) write(*,'(''uf:'',11i4)')is,(ifax(k19),k19=1,10)
        do iy=1,ny
        do ix=1,nx
          ufy(ix,iy,is)=(v(ix,iy,is,3)+v(ix+1,iy,is,3))
     :      *(u(ix,iy,is,3)*pp10(ix,iy,3)+u(ix,iy+1,is,3)
     :      *pp10(ix,iy+1,3))/dy4
        enddo
        enddo
      if(verbose.ge.3) write(*,'(''uf:'',11i4)')is,(ifax(k19),k19=1,10)
      enddo
c
      if(verbose.ge.3) write(*,'(''ufs:'',10i4)') (ifax(k19),k19=1,10)
      do is=1,ns
      do iy=2,ny
      do ix=1,nx
        ufs(ix,iy,is)=(u(ix,iy,is-1,3)+u(ix,iy,is,3))
     :    *(wsigp(ix,iy,is)+wsigp(ix+1,iy,is))/ds14(is)
      enddo
      enddo
      enddo
c
      do iy=2,ny
      do ix=1,nx
        ufs(ix,iy,0)=u(ix,iy,0,3)
     :    *(wsigp(ix,iy,0)+wsigp(ix+1,iy,0))/ds12(0)
        ufs(ix,iy,ns+1)=u(ix,iy,ns,3)
     :    *(wsigp(ix,iy,ns+1)+wsigp(ix+1,iy,ns+1))/ds12(ns+1)
      enddo
      enddo
c
      if(verbose.ge.3) write(*,'(''uflux:'',10i4)') (ifax(k19),k19=1,10)
      do is=0,ns
      do iy=2,ny
      do ix=1,nx
        uflux(ix,iy,is)=
     :    +(ufs(ix,iy,is+1)-ufs(ix,iy,is))
     :    +(ufx(ix+1,iy,is)-ufx(ix,iy,is))
     :    +(ufy(ix,iy,is)-ufy(ix,iy-1,is))
      enddo
      enddo
      enddo
c
c      deallocate(ufx)
c      deallocate(ufy)
c      deallocate(ufs)

      if(fcor.ne.0.) then
        do is=0,ns
        do iy=2,ny
        do ix=1,nx
          uflux(ix,iy,is)=uflux(ix,iy,is)-fcor*0.25*(v(ix,iy,is,3)
     :     *pp01(ix,iy,3)
     :     +v(ix+1,iy,is,3)*pp01(ix+1,iy,3)+v(ix,iy-1,is,3)
     :     *pp01(ix,iy-1,3)+v(ix+1,iy-1,is,3)*pp01(ix+1,iy-1,3))
        enddo
        enddo
        enddo
      endif
c
c vflux:
c
c      allocate(vfx(0:nx1,0:ny1,0:ns1))
c      allocate(vfy(0:nx1,0:ny1,0:ns1))
c      allocate(vfs(0:nx1,0:ny1,0:ns1))
      if(verbose.ge.3) write(*,'(''vfx:'',10i4)') (ifax(k19),k19=1,10)
      do is=0,ns
        do iy=1,ny
        do ix=1,nx
          vfx(ix,iy,is)=(u(ix,iy,is,3)+u(ix,iy+1,is,3))
     :       *(v(ix,iy,is,3)*pp01(ix,iy,3)+v(ix+1,iy,is,3)
     :       *pp01(ix+1,iy,3))/dx4
        enddo
        enddo
c
        do iy=1,ny1
        do ix=2,nx
          vfy(ix,iy,is)=(v(ix,iy-1,is,3)+v(ix,iy,is,3))
     :       *(v(ix,iy-1,is,3)*pp01(ix,iy-1,3)+v(ix,iy,is,3)
     :       *pp01(ix,iy,3))/dy4
        enddo
        enddo
      enddo
c
      do is=1,ns
      do iy=1,ny
      do ix=2,nx
        vfs(ix,iy,is)=(v(ix,iy,is-1,3)+v(ix,iy,is,3))
     :     *(wsigp(ix,iy,is)+wsigp(ix,iy+1,is))/ds14(is)
      enddo
      enddo
      enddo
c
      do iy=1,ny
      do ix=2,nx
        vfs(ix,iy,0)=v(ix,iy,0,3)
     :     *(wsigp(ix,iy,0)+wsigp(ix,iy+1,0))/ds12(0)
        vfs(ix,iy,ns+1)=v(ix,iy,ns,3)
     :     *(wsigp(ix,iy,ns+1)+wsigp(ix,iy+1,ns+1))/ds12(ns+1)
      enddo
      enddo
c
      if(verbose.ge.3) write(*,'(''vflux:'',10i4)') (ifax(k19),k19=1,10)
      do is=0,ns
      do iy=1,ny
      do ix=2,nx
        vflux(ix,iy,is)=(vfx(ix,iy,is)-vfx(ix-1,iy,is))
     :     +(vfs(ix,iy,is+1)-vfs(ix,iy,is))
     :     +(vfy(ix,iy+1,is)-vfy(ix,iy,is))
      enddo
      enddo
      enddo
c      deallocate(vfx)
c      deallocate(vfy)
c      deallocate(vfs)
c
      if(fcor.ne.0.) then
        do is=0,ns
        do iy=1,ny
        do ix=2,nx
          vflux(ix,iy,is)=vflux(ix,iy,is)+fcor*0.25*(u(ix,iy,is,3)
     :      *pp10(ix,iy,3)
     :      +u(ix-1,iy,is,3)*pp10(ix-1,iy,3)+u(ix,iy+1,is,3)
     :      *pp10(ix,iy+1,3)+u(ix-1,iy+1,is,3)*pp10(ix-1,iy+1,3))
        enddo
        enddo
        enddo
      endif
c
c      allocate(wfx(0:nx1,0:ny1,0:ns1))
c      allocate(wfy(0:nx1,0:ny1,0:ns1))
c      allocate(wfs(0:nx1,0:ny1,0:ns1))
      if(verbose.ge.3) write(*,'(''wfx:'',10i4)') (ifax(k19),k19=1,10)
      do is=1,ns
        do iy=1,ny1
        do ix=1,nx
          wfx(ix,iy,is)=(w(ix+1,iy,is,3)+w(ix,iy,is,3))
     :       *(u(ix,iy,is,3)*pp10(ix,iy,3)+u(ix,iy,is-1,3)
     :       *pp10(ix,iy,3))/dx4
        enddo
        enddo
        do iy=1,ny
        do ix=1,nx1
        wfy(ix,iy,is)=(w(ix,iy+1,is,3)+w(ix,iy,is,3))
     :     *(v(ix,iy,is,3)*pp01(ix,iy,3)+v(ix,iy,is-1,3)
     :     *pp01(ix,iy,3))/dy4
        enddo
        enddo
      enddo
c
      call extrax(nx1,ny1,wfx,0,nx1,1,ny1,0,ns)
      call extray(nx1,ny1,wfy,1,nx1,0,ny1,0,ns)
c
c or: extrapolate wfx(i=b+1) from wfx(i=b+0.5) and wfx(i=b)
c     equivalent to one-sided difference for w-flux convergence term
c     do 206 is=2,n-1
c     wfx05= w(1,is,3)*(uu(0,is,3)+uu(0,is-1,3)+uu(1,is  ,3)+
c    :       uu(1,is-1,3))*0.25
c     wfx(0,is)= 2.*wfx(1,is)-wfx05
c     wfx05=w(m1,is,3)*(uu(m1,is,3)+uu(m1,is-1,3)+uu(m,is,3)
c    :      +uu(m,is-1,3))*0.25
c     wfx(m1,is)=2.*wfx(m,is)-wfx05
c206  continue
c
      do is=0,ns
      do iy=1,ny1
      do ix=1,nx1
        wfs(ix,iy,is)=(w(ix,iy,is+1,3)+w(ix,iy,is,3))
     :    *(wsigp(ix,iy,is+1)+wsigp(ix,iy,is))/ds04(is)
      enddo
      enddo
      enddo
c
      if(verbose.ge.3) write(*,'(''wflux:'',10i4)') (ifax(k19),k19=1,10)
      do is=1,ns
      do iy=1,ny1
      do ix=1,nx1
        wflux(ix,iy,is)=(wfx(ix,iy,is)-wfx(ix-1,iy,is))
     :    +(wfy(ix,iy,is)-wfy(ix,iy-1,is))
     :    +(wfs(ix,iy,is)-wfs(ix,iy,is-1))
      enddo
      enddo
      enddo

c      deallocate(wfx)
c      deallocate(wfy)
c      deallocate(wfs)
c      deallocate(wsigp)
c
c this extrapolation gets rid of secondary boundary values:
c unused as such in the momentum equation.
c
      call extrah(nx1,ny1,uflux,1,nx,1,ny1,0,ns)
      call extrah(nx1,ny1,vflux,1,nx1,1,ny,0,ns)
c      call extrah(nx1,ny1,wflux,1,nx1,1,ny1,0,ns)
c
      if(verbose.ge.3) write(*,'(''flux::'',10i4)') (ifax(k19),k19=1,10)
      return
      end


      subroutine  hydros
c-----------------------------------------------------------------------
c integrate hydrostatic relation for phi field
c-----------------------------------------------------------------------
c
      use alloc
      implicit real*8(a-h,o-z)
c      include 'nh3dpa10.inc'
c
      do iy=1,ny1
      do ix=1,nx1
        phi(ix,iy,ns)=phig(ix,iy)-dsng*g*(pt(ix,iy,ns,3)
     :    +pt(ix,iy,ns-1,3))/((pts(ix,iy,ns)+pts(ix,iy,ns-1))
     :    *s(ix,iy,ns))
      enddo
      enddo
      do is=ns-1,0,-1
      do iy=1,ny1
      do ix=1,nx1
        phi(ix,iy,is)=phi(ix,iy,is+1)+ds0(is+1)*g
     :    *(pt(ix,iy,is,3)+pt(ix,iy,is+1,3))
     :    /((pts(ix,iy,is)+pts(ix,iy,is+1))*s(ix,iy,is+1))
c       write(99,*) 'meuphi',phi(ix,iy,is)
      enddo
      enddo
      enddo
      return
      end



      subroutine indrag
c
      use alloc
      implicit real*8(a-h,o-z)
c      include 'nh3dpa10.inc'
c
c calculate dragm and dragt (for use in diffu if(.not.raylei)) using
c taum and taut as the linear dissipation time scales for 2*dx waves.
c
      do is=0,ns+1
        hdamp(is)=hdampm
      enddo
c
      if(taum.ge.1.e8) then
         dragm=0.
      else
         dragm=(2.*dx)**2/(4.*pi*pi*taum)
      endif
c
      if(taut.ge.1.e8) then
         dragt=0.
      else
         dragt=(2.*dx)**2/(4.*pi*pi*taut)
      endif
c
      ixm=nx1/2
      iym=ny1/2
c
      ztop=(phis(ixm,iym,1)+phi(ixm,iym,1))/g
c
      do is=1,ns
        z=(phis(ixm,iym,is)+phi(ixm,iym,is))/g
        if(z.lt.zm) then
          idrmax=is-1
          go to 11
        endif
      enddo
11    continue
c
c
      do is=1,idrmax
        z=(phis(ixm,iym,is)+phi(ixm,iym,is))/g
        endrag(is)=(dragm+(dragt-dragm)
     :     *cos(pi05*abs((z-ztop)/(zm-ztop)))**2)
        if(endrag(is).eq.0.)then
          taudra(is) = 1.d10
        else
          taudra(is)=max(dt,(2.*dx)**2/(4.*pi*pi*endrag(is)))
        endif
        hdamp(is)=(hdampm+(hdampt-hdampm)
     :     *cos(pi05*abs((z-ztop)/(zm-ztop)))**2)
      enddo
c
!      write(nchn1,1000)
!      do is=1,idrmax
!        write(nchn1,1001) is,phis(ixm,iym,is)/g,endrag(is),taudra(is)
!     :    ,hdamp(is)
 !     enddo
!      do is=1,idrmax
!        write(0,*) is,phis(ixm,iym,is)/g,endrag(is),taudra(is)
!     :    ,hdamp(is)
!      enddo
c
c constant for diffu (deformation dependent diffusion)
c
      dkapa2=dkapa*dkapa
      do is=1,ns-1
        dd=min(dx,dy,0.5*(phis(ixm,iym,is-1)-phis(ixm,iym,is+1)
     :   +phi(ixm,iym,is-1)-phi(ixm,iym,is+1))/g)
        delka2(is)=dd*dd*dkapa2
      enddo
      delka2(ns)=delka2(ns-1)
c
c constants for radbch
c
c      rxmax=dsqrt(2.)/2.
c      rymax=dsqrt(2.)/2.
      rxmax=1.
      rymax=1.
c
1000  format(//,' is',t10,'z',t30,'enhanced drag',t50,'rayleigh times'
     :   ,t70,'hdamp')
1001  format(1x,i3,t10,f8.0,t30,f8.1,t50,e14.7,t66,f8.4)
      return
      end






      subroutine init0
c-----------------------------------------------------------------------
c puts all arrays to zero
c-----------------------------------------------------------------------
c
      use alloc
      implicit real*8(a-h,o-z)
c      include 'nh3dpa10.inc'
c
      u=0.
      v=0.
      w=0.
      wsig=0.
      pt=0.
      thetav=0.                                                                 VS,05.2007
      phi=0.
      difunu=0.
      difunv=0.
      difunw=0.
      difunt=0.
      if(qif.ne.0.) then
        difunq=0.
        qv=0.
        qvs=0.
        cond=0.
        evap=0.
        if(ifqc.ne.0) then
          qc=0.
          if(ifqr.ne.0) then
            qr=0.
            auto=0.
            col=0.
            vrain=0.
          endif
        endif
      endif
      s=0.
      us=0.
      vs=0.
      pts=0.
      phis=0.
      tems=0.
      tsuf=0.
      dpp=0.
      phig=0.
c
      dpdx10=0.
      dpdy01=0.
      ppdx10=0.
      ppdy01=0.
      psuf=0.
      ppdx=0.
      ppdy=0.
c
      pgradx=0.
      pgrady=0.
      pgrads=0.
c
      ubx=0.
      ubx=0.
      vbx=0.
      vbx=0.
      uby=0.
      uby=0.
      vby=0.
      vby=0.
      ptbx=0.
      ptby=0.
      ppbx=0.
      ppbx=0.
      ppby=0.
      ppby=0.
      
	if (allocated(Srad_surf)) then
       Srad_surf = 0.                                                           VS,06,2007
       Lrad_surf = 0.                                                           VS,06,2007
	endif
	tsurf     = 273.                                                          VS,06,2007

      qsbx = 0.                                                                 DM,05.2007
      qsby = 0.                                                                 DM,05.2007
      qscc = 0.                                                                 DM,05.2007
c
      return
      end



      subroutine initia
c-----------------------------------------------------------------------
c specification of initial fields
c-----------------------------------------------------------------------
c
      use alloc
      use refstate

      implicit real*8(a-h,o-z)
	integer ix,iy,is
	real*8, external:: qsati
c      include 'nh3dpa10.inc'

c
      logical sim,nao
      parameter (sim=.true.,nao=.false.)
c
c linear interpolation of initial wind profile (given) in z:
c
      call upuvt
c
      u(:,:,:,2)=us
      u(:,:,:,3)=us
      v(:,:,:,2)=vs
      v(:,:,:,3)=vs
      pt=0.
      if(qif.ne.0.) then
        qv(:,:,:,1)=qvs
        qv(:,:,:,2)=qvs
        qv(:,:,:,3)=qvs
        cond=0.
        evap=0.
	!do is=0,ns
        !do iy=0,ny1
        !  do ix=0,nx1
	  !    p=sigma0(is)*pp(ix,iy,2)+ptop
       !     t=(pts(ix,iy,is)+pt(ix,iy,is,2))*(p/p00)**akapa
	 !if(t.lt.273.15) then
	! if(qv(ix,iy,is,1).gt.qsati(t,p)) qv(ix,iy,is,1)=0.40*qsati(t,p)
	! if(qv(ix,iy,is,2).gt.qsati(t,p)) qv(ix,iy,is,2)=0.40*qsati(t,p)
	 !if(qv(ix,iy,is,3).gt.qsati(t,p)) qv(ix,iy,is,3)=0.40*qsati(t,p)
	!endif
	  
	!enddo
	!enddo
	!enddo
	
        if(ifqc.ne.0) then
          qc=0.
          if(ifqr.ne.0) then
            qr=0.
            auto=0.
            col=0.
            vrain=0.
            prec=0.
          endif
	    if (ifqi.ne.0) then
	qci=0.
	qsn=0.
	vdepi=0.
	vdeps=0.
	vini=0.
	imlt=0.
	hmfrz=0.
	sacrw=0.
	sacrwr=0.
	sbercw=0.
	sberci=0.
	ibercw=0.
	sagg=0.
	saci=0.
	raci=0.
	iacr=0.
	smlt=0.
	sacrr=0.
	difunqci=0.
	difunqsn=0.
	vsnow=0.
	precsnow=0.

	    endif
        endif
      endif

c
c initial mass fluxes through boundaries:
c calculation of the integral of u dy dp/g = uu /g dy ds
c
      xlmflx=0.0
      do is=1,ns-1
      do iy=2,ny
        xlmflx=xlmflx+u(1,iy,is,2)*pp10(1,iy,2)*ds1y(is)
      enddo
      enddo
      xlmflx=xlmflx/g
c
      xrmflx=0.0
      do is=1,ns-1
      do iy=2,ny
        xrmflx=xrmflx+u(nx,iy,is,2)*pp10(nx,iy,2)*ds1y(is)
      enddo
      enddo
      xrmflx=xrmflx/g
c
      ylmflx=0.0
      do is=1,ns-1
      do ix=2,nx
        ylmflx=ylmflx+v(ix,1,is,2)*pp01(ix,1,2)*ds1x(is)
      enddo
      enddo
      ylmflx=ylmflx/g
c
      yrmflx=0.0
      do is=1,ns-1
      do ix=2,nx
        yrmflx=yrmflx+v(ix,ny,is,2)*pp01(ix,ny,2)*ds1x(is)
      enddo
      enddo
      yrmflx=yrmflx/g
c
      ttmflx=xlmflx-xrmflx+ylmflx-yrmflx
c
      if(iniout) then
         write(nchn1,4) xlmflx,xrmflx,ylmflx,yrmflx,ttmflx
      endif
c
c transformation of the initial field p*(u,v) to a new one with
c the nondivergent vertically averaged component
c
c initialization of w and wsig.
c
      if(adjust) then
      call disc(us,vs,disb,ii,jj)
c      write(nchn1,*)'initial unbalance (=-dp/dt):',disb,' ix=',ii,
c     :' iy=',jj
c
c      call wri3ar(w(0,0,0,3),2,nx,2,ny,1,ns-1,'w-nduvw-',tk,iomod,nx,ny)
c      call nduvw(0,ww1,ww2)
      call nduvw(0)
c      call wri3ar(w(0,0,0,3),2,nx,2,ny,1,ns-1,'w-nduvw+',tk,iomod,nx,ny)
c
      do is=0,ns1
      do iy=0,ny1
      do ix=0,nx1
        u(ix,iy,is,2)=u(ix,iy,is,3)
        v(ix,iy,is,2)=v(ix,iy,is,3)
        us(ix,iy,is)=u(ix,iy,is,3)
        vs(ix,iy,is)=v(ix,iy,is,3)
        wsig(ix,iy,is,2)=wsig(ix,iy,is,3)
        w(ix,iy,is,2)=w(ix,iy,is,3)
      enddo
      enddo
      enddo
c
c      call wri3ar(w(0,0,0,3),2,nx,2,ny,1,ns-1,'w-disc- ',tk,iomod,nx,ny)
      call disc(us,vs,disb,ii,jj)
c      call wri3ar(w(0,0,0,3),2,nx,2,ny,1,ns-1,'w-disc+ ',tk,iomod,nx,ny)
      write(nchn1,*)'unbalance after nduvw=', disb,' ix=',ii,' iy=',jj
cc*
c      write(98,*)'u algvali'
c      call test2_0(us,nx1,ny1,ns1)
c      write(98,*)'v algvali'
c      call test2_0(vs,nx1,ny1,ns1)
c      stop
cc*
      endif
c
c this is for a particular experiment:                                  #exp#
c and ensures constant mass flux at each sigma layer (x direction).     #exp#
c
c      do 60 iy=0,ny+1
c      do 60 is=1,ns-1
c      do 60 ix=2,nx+1
c      u(ix,iy,is,3)=uu(1,iy,is,3)/pp10(ix,iy,3)
c60    continue
c      do 61 iy=0,ny1
c      do 61 ix=0,nx1
c      u(ix,iy,ns,3)=u(ix,iy,ns-1,3)
c      u(ix,iy,0,3)=u(ix,iy,1,3)
c61    continue
c      call extrah(nx1,ny1,u(0,0,0,3),1,nx1,1,ny1,0,ns-1)
c      do 62 is=0,ns
c      do 62 iy=0,ny1
c      do 62 ix=0,nx1
c      u(ix,iy,is,2)=u(ix,iy,is,3)
c62    continue
cc
c      do 160 ix=0,nx+1
c      do 160 is=1,ns-1
c      do 160 iy=2,ny+1
c      v(ix,iy,is,3)=vv(ix,1,is,3)/pp01(ix,iy,3)
c160   continue
c      do 161 iy=0,ny1
c      do 161 ix=0,nx1
c      v(ix,iy,ns,3)=v(ix,iy,ns-1,3)
c      v(ix,iy,0,3)=v(ix,iy,1,3)
c161   continue
c      call extrah(nx1,ny1,v(0,0,0,3),1,nx1,1,ny1,0,ns-1)
c      do 162 is=0,ns
c      do 162 iy=0,ny1
c      do 162 ix=0,nx1
c      v(ix,iy,is,2)=v(ix,iy,is,3)
c162   continue
c
c initial temperature perturbation (depends on experiment!)             #exp#
c
c delpt: max temperature inversion (at x=0.)
c depinv: depth of inversion
c
c      delpt=-37.
c      depinv=200.
cc
c      iy=ny1/2
c      do 70 ix=2,nx
c      z0=0.5*(phis(ix,iy,ns-1)+phis(ix,iy,ns))/g
c      do 69 is=0,ns-1
c      z=phis(ix,iy,is)/g
c      x=(ix-1.5)*dx
c      ptforc(ix,is)=delpt*(xl-x)/xl*exp(-(z-z0)/depinv)
c69    continue
c      ptforc(ix,ns)=ptforc(ix,ns-1)
c      ptforc(ix,0)=ptforc(ix,1)
c70    continue
c      do 68 is=0,ns
c      ptforc(1,is)=ptforc(2,is)
c      ptforc(nx1,is)=ptforc(nx,is)
c68    continue
c
c
c      do 71 iy=2,ny
c      do 71 ix=2,nx
c      do 71 is=0,ns-1
cc      z=phis(ix,iy,is)/g
cc      x=(ix-1.5)*dx
cc      y=(iy-1.5)*dy
cc      rad=dsqrt(((x-centrx)/widtx)**2
cc     :   +((y-centry)/widty)**2
cc     :   +((z-centrz)/widtz)**2)
cc      pt(ix,iy,is,3)=pt(ix,iy,is,3)+amppt*cos(pi05*rad)
cc     :   *dim(1.,rad)/(1.-rad+zero0)
c      pt(ix,iy,is,3)=ptforc(ix,is)
c71    continue
c
c Introduce initial perturbation on pt
c
      if(ifblob.eq.1) call blob
      
!      call perturb_profile
c
      call extrah(nx1,ny1,pt(0,0,0,3),1,nx1,1,ny1,0,ns-1)
      do iy=0,ny1
      do ix=0,nx1
        pt(ix,iy,ns,3)=pt(ix,iy,ns-1,3)
        pt(ix,iy,0,3)=pt(ix,iy,1,3)
      enddo
      enddo
c
      do is=0,ns1
      do iy=0,ny1
      do ix=0,nx1
        pt(ix,iy,is,2)=pt(ix,iy,is,3)
      enddo
      enddo
      enddo
c
c integrate continuity equation
c
      dpp=0.
      do is=1,ns-1
      do iy=1,ny1
      do ix=1,nx1
        dpp(ix,iy)=dpp(ix,iy)-(u(ix,iy,is,2)*pp10(ix,iy,2)
     :    -u(ix-1,iy,is,2)*pp10(ix-1,iy,2))*ds1dx(is)
     :    -(v(ix,iy,is,2)*pp01(ix,iy,2)
     :    -v(ix,iy-1,is,2)*pp01(ix,iy-1,2))*ds1dy(is)
      enddo
      enddo
      enddo
c
c       write(*,*) 'initia: call raduv'
      call raduv(nx,ny,ns
     :   ,u(0,0,0,3),u(0,0,0,2),u(0,0,0,2)
     :   ,ubx,ubx2,uby,uby2,ucc,ucc2,ds1y,yl,pp10(0,0,4)
     :   ,v(0,0,0,3),v(0,0,0,2),v(0,0,0,2)
     :   ,vbx,vbx2,vby,vby2,vcc,vcc2,ds1x,xl,pp01(0,0,4)
     :   ,ttmin,masout,uubout,prt,nchn1,ds0,ds1,s,nao
     :   ,xlcor,xrcor,ylcor,yrcor,ifluxcor,dpmeddt,dx,dy)
c       write(*,*) 'initia: call radbch pt'
      call radbch(ptcc,ptbx,ptby,pt(0,0,0,3),pt(0,0,0,2),pt(0,0,0,2)
     :   ,nx1,ny1,ns1,1,ns-1,1,nx1,1,ny1)
      if(qif.ne.0.) then
        call radbch(qvcc,qvbx,qvby,qv(0,0,0,3),qv(0,0,0,2),qv(0,0,0,2)
     :    ,nx1,ny1,ns1,1,ns-1,1,nx1,1,ny1)
        if(ifqc.ne.0) then
          call radbch(qccc,qcbx,qcby,qc(0,0,0,3),qc(0,0,0,2),qc(0,0,0,2)
     :      ,nx1,ny1,ns1,1,ns-1,1,nx1,1,ny1)
          if(ifqr.ne.0) then
            call radbch(qrcc,qrbx,qrby,qr(0,0,0,3),qr(0,0,0,2)
     :        ,qr(0,0,0,2),nx1,ny1,ns1,1,ns-1,1,nx1,1,ny1)
	      if(ifqi.ne.0) then                                                  DC,11.2009
              call radbch(qcicc,qcibx,qciby,qci(0,0,0,3),qci(0,0,0,2)           DC,11.2009
     :          ,qci(0,0,0,2),nx1,ny1,ns1,1,ns-1,1,nx1,1,ny1)                   DC,11.2009
	        call radbch(qsncc,qsnbx,qsnby,qsn(0,0,0,3),qsn(0,0,0,2)           DC,11.2009
     :          ,qsn(0,0,0,2),nx1,ny1,ns1,1,ns-1,1,nx1,1,ny1)                   DC,11.2009
	      endif                                                               DC,11.2009
          endif
        endif
      endif

c       write(*,*) 'initia: call radbch pp'
      call radbch(ppcc,ppbx,ppby,pp(0,0,4),pp(0,0,3),pp(0,0,2)
     :   ,nx2,ny2,0,0,0,1,nx1,1,ny1)
c       write(*,*) 'initia: call radbch pp2'
      call radbch(ppcc2,ppbx2,ppby2,pp(0,0,4),pp(0,0,3),pp(0,0,2)
     :   ,nx2,ny2,0,0,0,0,nx2,0,ny2)
c
c print out initial sounding profile
c
c      if(iniout) then
c        write(nchn1,'(//a/)') ' initial sounding profile'
c        call wripro
c        tk=0.
c        if(phsout) call wri2ar(dpp,1,nx1,1,ny1,'dpp      ',tk,nx,ny)
c      endif
c
c initialize work variables:
c
      return
4     format(//' xlmflx=',e15.8,' xrmflx=',e15.8,' ylmflx=',e15.8,
     :   ' yrmflx=',e15.8, ' ttmflx=',e15.8)
      end



      subroutine inpos2(jobphy)
c
c-----------------------------------------------------------------------
c solves poisson equation:
c  2      2    2      2
c d phi/dx  + d phi/dy  = f(x,y)
c
c with dirichlet boundary conditions (phi=0), using fourier method.
c
c note the dimensions of phi00 (only 2:nx and 2:ny are
c inside the physical domain).
c-----------------------------------------------------------------------
c
c
      use alloc
      implicit real*8(a-h,o-z)
c      include 'nh3dpa10.inc'
c
c      common/triphi/a1phi(ny-1),a2phi(nx-1,ny-1),a3phi(ny-1)
c     :   ,ctphi(nx-1,ny-1),dtphi(nx-1,ny-1)
c
      dimension phi00(0:nx1,0:ny1)
c
c      dimension sines(nx-1,nx-1)
c      double precision dpi,sine(0:nx-1)
c      save sines
c
      if(allocated(a1phi)) then
        deallocate(a1phi,a2phi,a3phi,ctphi,dtphi,sine,sines)
      endif
      allocate(a1phi(ny-1),a2phi(nx-1,ny-1),a3phi(ny-1)
     :  ,ctphi(nx-1,ny-1),dtphi(nx-1,ny-1),sine(0:nx-1)
     :  ,sines(nx-1,ny-1))
      dpi=4.d0*datan(1.d0)
c
c initialization of arrays and preparation of the solution:
c
      do iy=1,ny-1
        a1phi(iy)=1./dyy
        a3phi(iy)=a1phi(iy)
      enddo
c
      do ix=1,nx-1
        sine(ix)=-2./dyy-4./dxx*dsin(dble(ix)*dpi/dble(2*nx))**2
      enddo
c
      do iy=1,ny-1
      do ix=1,nx-1
        a2phi(ix,iy)=sine(ix)
      enddo
      enddo
c
c case of homogenous neuman b.c. (non-homogeneity to be treated outside)
c
      if(jobphy.eq.0) then
         do ix=1,nx-1
           a2phi(ix,1)=a2phi(ix,1)+a1phi(1)
           a2phi(ix,ny-1)=a2phi(ix,ny-1)+a3phi(ny-1)
         enddo
      endif
c
      call iytri2(a1phi,a2phi,a3phi,nx-1,ny-1,ctphi,dtphi)
c
c initialization of parameters for sintra:
c
      call instra(sines,sine,nx)
c
      return
c
c
c
      entry pos2(phi00,jobphy)
c
      allocate(work1(0:nx1,0:ny1,0:ns1))

      call sintrx(phi00,work1,sines,nx,ny-1,+1)

c      tk=0.
c      call wri2ar(phi00 ,1,nx1,1,ny1,'phi00a  ',tk,nx,ny)
      call ytri2(a1phi,a2phi,a3phi,phi00,phi00,nx-1,ny-1,ctphi,dtphi)
c      call wri2ar(phi00 ,1,nx1,1,ny1,'phi00b  ',tk,nx,ny)
      call sintrx(phi00,work1,sines,nx,ny-1,-1)
c      call wri2ar(phi00 ,1,nx1,1,ny1,'phi00c  ',tk,nx,ny)
c
c impose b.c.
c
      if(jobphy.eq.0) then
        do ix=1,nx1
          phi00(ix,1)=phi00(ix,2)
          phi00(ix,ny1)=phi00(ix,ny)
        enddo
      else
        do ix=1,nx1
          phi00(ix,1)=0.
          phi00(ix,ny1)=0.
        enddo
      endif
c
      do iy=1,ny1
        phi00(1,iy)=0.
        phi00(nx1,iy)=0.
      enddo

      deallocate(work1)
c
      return
      end



      subroutine inpos3(c,a2)
c-----------------------------------------------------------------------
c
c note that the dimensions of wdphi do not coincide with those in the
c subrout. dscfft and tri3
c as it is assumed that the physical domain starts ,there, at
c wdphi(0.5x,0.5y) whereas it is wdphi(1.5x,1.5y)
c
c note that c must be reserved between sucessive calls to
c pos3, after beeing defined in inpos3.
c-----------------------------------------------------------------------
c
      use alloc
      implicit real*8(a-h,o-z)
c      include 'nh3dpa10.inc'
c
c      dimension a1(ns-1),a2(nx-1,ny-1,ns-1),a3(ns-1),a21(nx-1,ny-1)
c     :   ,aves(ns)

      dimension a2(nx-1,ny-1,ns-1)
      dimension c(nx-1,ny-1,ns-1)
      dimension wdphi(0:nx1,0:ny1,0:ns1)
      dimension ww1(0:nx1,0:ny1,0:ns1),ww2(0:nx1,0:ny1,0:ns1)
c
c      save a1,a3,a21
c
c initialization of temporary variables:
c
      do is=1,ns
        aves(is)=0.0
        do iy=1,ny1
        do ix=1,nx1
          aves(is)=aves(is)+s(ix,iy,is)*s(ix,iy,is)
        enddo
        enddo
        aves(is)=aves(is)/(nx1*ny1*ds0(is))
      enddo
c
c calculate l.h.s. coefficients of the transformed equation
c as the elements of a tridiagonal array.
c
      do is=1,ns-1
        a1pos3(is)=aves(is)/ds1(is)
        a3pos3(is)=aves(is+1)/ds1(is)
      enddo
      do is=1,ns-1
      do iy=1,ny-1
      do ix=1,nx-1
        a2(ix,iy,is)=-xylamb(ix,iy)-(aves(is+1)+aves(is))/ds1(is)
      enddo
      enddo
      enddo
c
c boundary conditions (neuman in the top, dirichlet at the bottom)
c
      do iy=1,ny-1
      do ix=1,nx-1
        a2(ix,iy,1)=a2(ix,iy,1)+a1pos3(1)
      enddo
      enddo
c
c preparation of the solution of the tridiagonal systems of equations:
c
      call intri3(a1pos3,a2,a3pos3,nx-1,ny-1,ns-1,c,a21pos3)
c
      return
c
c solve poisson equation:
c
      entry pos3(wdphi,c,ww1,ww2)
c
c perform double cosine transform of the r.h.s (forcing) of
c the poisson eq.
c
      call dscfft(wdphi,ww1,ww2,ww2,trigsx,ifax,sipcox,simcox,nx-1
     :  ,trigsy,ifay,sipcoy,simcoy,ny-1,ns-1,+1)
c
c on exit wdphi contains the transformed r.h.s. (in spectrum space).
c
c solve the linear tridiagonal system of equations in parallel for
c (nx-1)*(ny-1) wave components .
c
      call tri3(a1pos3,a3pos3,wdphi,nx-1,ny-1,ns-1,c,ww1,a21pos3)
c
c transform back to physical space for solution of poisson equation.
c
      call dscfft(wdphi,ww1,ww2,ww2,trigsx,ifax,sipcox,simcox,nx-1
     :  ,trigsy,ifay,sipcoy,simcoy,ny-1,ns-1,-1)
c
c fulfil the boundary condition requirement.
c
      do iy=2,ny
      do ix=2,nx
        wdphi(ix,iy,0)=wdphi(ix,iy,1)
        wdphi(ix,iy,ns)=0.0
      enddo
      enddo
c
      call extrah(nx1,ny1,wdphi,1,nx1,1,ny1,0,ns)
c
      do is=0,ns
        do iy=0,ny1
          wdphi(0,iy,is)=0.
        enddo
        do ix=0,nx1
          wdphi(ix,0,is)=0.
        enddo
      enddo
c
      return
      end




      subroutine tri3b(a1,a2,a3,a4,a4w,nx,ny,ns,nc)
c
c   to solve nx*ny systems of ns linear tridiagonal equations
c   using direct reccurences
c   nc = 1 -- homogeneous conditions for both initial values
c             at the bottom (hom. neumann and hom. dirichlet)
c   nc = 2 -- h.neumann at the bottom, h.dirichlet in the top
c   nc = 3 -- h.neumann at both sides, exept ix=iy=1, when both
c             hn and hd are at the bottom
c
c   a1, a2,a3 -- coefficients of the reccurence (a1=a1_miranda,
c                a2=(a2/a1)_miranda, a3=(a3/a1)_miranda)
c   a4 --  r.h.s. (input) and solution (output)
c
      implicit real*8(a-h,o-z)
      dimension a1(ns),a2(nx,ny,ns),a3(ns),a4(-1:nx+1,-1:ny+1,0:ns+2)
     :         ,a4w(0:ns+2)
c
      nss=ns+1
      if(nc.ne.1) go to 1000
      do iy=1,ny
      do ix=1,nx
        gisp1=a4(ix,iy,nss-1)
        a4(ix,iy,nss)=0.0
        a4(ix,iy,nss-1)=0.0
        do is=nss-2,0,-1
          gisp=a4(ix,iy,is)
          a4(ix,iy,is)=gisp1/a1(is+1)-a2(ix,iy,is+1)*a4(ix,iy,is+1)-
     :      a3(is+1)*a4(ix,iy,is+2)
          gisp1=gisp
        enddo
      enddo
      enddo
c
      return
c
1000  continue
c
c case nc=2: homog. neumann at the bottom, hom. dirichlet in the top
c
      if(nc.ne.2) go to 2000
c
      do iy=1,ny
      do ix=1,nx
        gisp1=a4(ix,iy,nss-1)
        a4(ix,iy,nss)=0.0
        a4(ix,iy,nss-1)=0.0
        a4w(nss)=1.0
        a4w(nss-1)=1.0
        do is=nss-2,0,-1
          gisp=a4(ix,iy,is)
          a4w(is)=-a2(ix,iy,is+1)*a4w(is+1)-a3(is+1)*a4w(is+2)
          a4(ix,iy,is)=gisp1/a1(is+1)-a2(ix,iy,is+1)*a4(ix,iy,is+1)-
     :      a3(is+1)*a4(ix,iy,is+2)
          gisp1=gisp
        enddo
        qq=-a4(ix,iy,0)/a4w(0)
        do is=0,nss
          a4(ix,iy,is)=a4(ix,iy,is)+a4w(is)*qq
        enddo
      enddo
      enddo
c
      return
c
2000  continue
c
c case nc=3: homog. neumann at the bottom, hom. neumann in the top
c exept ix=iy=1, when both hom neumann and hom dirichlet are employed
c at the bottom (because hom neumann at both ends has no solution in
c this case
c
      if(nc.ne.3) go to 3000
c
      do iy=1,ny
      do ix=1,nx
        gisp1=a4(ix,iy,nss-1)
        a4(ix,iy,nss)=0.0
        a4(ix,iy,nss-1)=0.0
        a4w(nss)=1.0
        a4w(nss-1)=1.0
        do is=nss-2,0,-1
          gisp=a4(ix,iy,is)
          a4w(is)=-a2(ix,iy,is+1)*a4w(is+1)-a3(is+1)*a4w(is+2)
          a4(ix,iy,is)=gisp1/a1(is+1)-a2(ix,iy,is+1)*a4(ix,iy,is+1)-
     :      a3(is+1)*a4(ix,iy,is+2)
          gisp1=gisp
        enddo
        if(ix*iy.ne.1) then
          qq=a4w(1)-a4w(0)
          if(abs(qq).le.1.e-25) then
            write(0,*)'nh3dFatalError=tri3b:singular matrix at ix='
     :        ,ix,'  iy=',iy
            stop
          endif
          qq=-(a4(ix,iy,1)-a4(ix,iy,0))/qq
          do is=0,nss
            a4(ix,iy,is)=a4(ix,iy,is)+a4w(is)*qq
          enddo
        endif
      enddo
      enddo
c
      return
c
3000  continue
      if(nc.ne.4) go to 4000
c
c case nc=4: like nc=3, exept ix=iy=1, where hom neumann is employed
c at the bottom and hom dirichlet is employed in the top (like in nc=2)
c
      do iy=1,ny
      do ix=1,nx
        gisp1=a4(ix,iy,nss-1)
        a4(ix,iy,nss)=0.0
        a4(ix,iy,nss-1)=0.0
        a4w(nss)=1.0
        a4w(nss-1)=1.0
        do is=nss-2,0,-1
          gisp=a4(ix,iy,is)
          a4w(is)=-a2(ix,iy,is+1)*a4w(is+1)-a3(is+1)*a4w(is+2)
          a4(ix,iy,is)=gisp1/a1(is+1)-a2(ix,iy,is+1)*a4(ix,iy,is+1)-
     :      a3(is+1)*a4(ix,iy,is+2)
          gisp1=gisp
        enddo
        if(ix*iy.ne.1) then
          qq=a4w(1)-a4w(0)
          if(abs(qq).le.1.e-25) then
            write(0,*)'nh3dFatalError=tri3b:singular matrix at ix='
     :        ,ix,'  iy=',iy
            stop
          endif
          qq=-(a4(ix,iy,1)-a4(ix,iy,0))/qq
        else
          qq=-a4(ix,iy,0)/a4w(0)
        endif
        do is=0,nss
          a4(ix,iy,is)=a4(ix,iy,is)+a4w(is)*qq
        enddo
      enddo
      enddo
      return
c
4000  continue
      write(0,*)'nh3dFatalError= tri3b: nc has wrong value'
      stop
      end


      subroutine nduvw(i)
c
c turns velocity u(ix,iy,is,3),v(ix,iy,is,3) into a field wih a
c nondivergent vertically integrated
c component extracting from initial vector (u,v) its irrotational
c vertically integrated component.
c computes vertical winds w and wsig from diagnostic relations
c
c i=0 - initial initialization in the extended area 0:nx1,0:ny1
c i>0 - alters u,v in the internal points u: 2,nx-1:2:ny
c       v: 2,nx, 2,ny-1 and preservers in boundary points where
c       u,v are updated by the uvbc
c employes the routine dscfft to solve the 2d poisson equation
c
      use alloc
      implicit real*8(a-h,o-z)

c      dimension um(0:nx1,0:ny1), vm(0:nx1,0:ny1), f(0:nx1,0:ny1,0:3),
c     :          ww1(0:nx1,0:ny1,0:3),ww2(0:nx1,0:ny1,0:3)
      real(kind(0.d0)),allocatable:: um(:,:),vm(:,:),f(:,:,:)
      real(kind(0.d0)),allocatable,dimension(:,:,:):: ww1,ww2

      allocate(ww1(0:nx1,0:ny1,0:3))
      allocate(ww2(0:nx1,0:ny1,0:3))

      allocate(um(0:nx1,0:ny1))
      allocate(vm(0:nx1,0:ny1))
      allocate(f(0:nx1,0:ny1,0:3))
c
c third index for fields f ww1 and ww2 is introduced to satisfy
c requirements set by the routine dscfft
c
c updating vertical velocity field
c
      do ix=0,nx1
      do iy=0,ny1
      do is=0,ns1
        w(ix,iy,is,2)=w(ix,iy,is,3)
        wsig(ix,iy,is,2)=wsig(ix,iy,is,3)
      enddo
      enddo
      enddo
c
c new wsig and r.h.s. for 2d poisson eq.
c
      do ix=1,nx1
      do iy=1,ny1
        f(ix,iy,1)=0.
        wsig(ix,iy,1,3)=0.
        do is=2,ns
          im=is-1
          f(ix,iy,1)=f(ix,iy,1)+
     :      ((u(ix,iy,im,3)*pp10(ix,iy,3)-u(ix-1,iy,im,3)
     :      *pp10(ix-1,iy,3))/dx
     :      +(v(ix,iy,im,3)*pp01(ix,iy,3)-v(ix,iy-1,im,3)
     :      *pp01(ix,iy-1,3))/dy)*ds1(is)
          wsig(ix,iy,is,3)=-f(ix,iy,1)/pp(ix,iy,3)
        enddo
        do is=2,ns
          wsig(ix,iy,is,3)=wsig(ix,iy,is,3)-sigma1(is)*wsig(ix,iy,ns,3)
        enddo
        wsig(ix,iy,0,3)=-wsig(ix,iy,2,3)
        wsig(ix,iy,ns1,3)=-wsig(ix,iy,ns-1,3)
      enddo
      enddo
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc
cc   internal points of the domain are (2:nx , 2:ny)
cc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c perform double cosine transform of the r.h.s (forcing) of
c the poisson eq.
c
      call dscfft(f,ww1,ww2,ww2,trigsx,ifax,sipcox,simcox,nx-1
     :   ,trigsy,ifay,sipcoy,simcoy,ny-1,1,+1)
c
c f contains on exit the transformed r.h.s.
c
c solution of the poisson equation in spectral space for
c (nx-1)*(ny-1) wave components .
c
c value of f(2,2,1) does not affect wind field
c
      f(2,2,1) = 0.0
c      f(2,2,1) = 1.0
      do ix = 3,nx
        f(ix,2,1)=-f(ix,2,1)/xylamb(ix-1,1)
      enddo
      do iy = 3,ny
        f(2,iy,1)=-f(2,iy,1)/xylamb(1,iy-1)
      enddo
      do ix=3,nx
      do iy=3,ny
        f(ix,iy,1)=-f(ix,iy,1)/xylamb(ix-1,iy-1)
      enddo
      enddo
c
c
c transformation of the solution of poisson equation
c back to the physical space.
c
      call dscfft(f,ww1,ww2,ww2,trigsx,ifax,sipcox,simcox,nx-1
     :   ,trigsy,ifay,sipcoy,simcoy,ny-1,1,-1)
c
c continuation of the solution to the extended domain (1:nx1,1:ny1)
c using hom. neumann bc for the potential f of correction field:
c correction field grad f must not affect the calculated from previous
c time-step boundary values of (u,v) (uvbc). that means, grad f
c must be zero on boundaries.
c
      do iy=2,ny
        f(nx1,iy,1)=f(nx,iy,1)
        f(1,iy,1)=f(2,iy,1)
      enddo
      do ix=1,nx1
        f(ix,ny1,1)=f(ix,ny,1)
        f(ix,1,1)=f(ix,2,1)
      enddo
c
c irrotational component in the extended area (0:nx1,0:ny1)
c
      do ix=1,nx
      do iy=1,ny
        um(ix,iy)=(f(ix+1,iy,1)-f(ix,iy,1))/pp10(ix,iy,3)/dx
        vm(ix,iy)=(f(ix,iy+1,1)-f(ix,iy,1))/pp01(ix,iy,3)/dy
      enddo
      enddo
      do iy =1,ny
        um(0,iy)=um(1,iy)
        um(nx1,iy)=um(nx,iy)
        vm(0,iy)=vm(1,iy)
        vm(nx1,iy)=vm(nx,iy)
      enddo
      do ix=0,nx1
         um(ix,0)=um(ix,1)
         um(ix,ny1)=um(ix,ny)
         vm(ix,0)=vm(ix,1)
         vm(ix,ny1)=vm(ix,ny)
      enddo
c
c extraction of 2d irrotational component from the initial wind field
c (this does not change horizonalt wind, if it is initially 2d-nondivergent)
c
c
      if(i.eq.0)then
c
        do ix=1,nx1
        do iy=1,ny1
        do is=0,ns-1
          u(ix,iy,is,3)=u(ix,iy,is,3)-um(ix,iy)
          v(ix,iy,is,3)=v(ix,iy,is,3)-vm(ix,iy)
        enddo
        enddo
        enddo
      else
        do is=0,ns-1
          do ix=2,nx-1
          do iy=2,ny-1
            u(ix,iy,is,3)=u(ix,iy,is,3)-um(ix,iy)
            v(ix,iy,is,3)=v(ix,iy,is,3)-vm(ix,iy)
          enddo
          enddo
          do ix=2,nx-1
            u(ix,ny,is,3)=u(ix,ny,is,3)-um(ix,ny)
          enddo
          do iy=2,ny-1
            v(nx,iy,is,3)=v(nx,iy,is,3)-vm(nx,iy)
          enddo
        enddo
      endif
c
c vertical wind computation
c
      do ix=1,nx1
      do iy=1,ny1
c impossible to put ix,iy=0, because s is not
c defined at these surfaces
c bug CORRECTED 2000/07/20
        do is=2,ns
          w(ix,iy,is,3)=-(
     :      wsig(ix,iy,is,3)+sigma1(is)*
     :      (
     :      (u(ix,iy,is,3)+u(ix,iy,is-1,3))*ppdx10(ix,iy,2)+
     :      (v(ix,iy,is,3)+v(ix,iy,is-1,3))*ppdy01(ix,iy,2)+
     :      (u(ix-1,iy,is,3)+u(ix-1,iy,is-1,3))*ppdx10(ix-1,iy,2)+
     :      (v(ix,iy-1,is,3)+v(ix,iy-1,is-1,3))*ppdy01(ix,iy-1,2)
     :      )/4.
     :      )/s(ix,iy,is)
        enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c original bcs, which yield wrong results, if dsigma*p_* >> p_top
c
        w(ix,iy,1,3)=0.
        w(ix,iy,0,3)=-w(ix,iy,2,3)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      w(ix,iy,0,3)=w(ix,iy,2,3)
c      w(ix,iy,1,3)=w(ix,iy,2,3)
      enddo
      enddo
c      call wri3ar(ppdx10(0,0,2),2,nx,2,ny,0,0,'ppdx10  ',tk,iomod,nx,ny)
c      call wri3ar(ppdy01(0,0,2),2,nx,2,ny,0,0,'ppdy01  ',tk,iomod,nx,ny)
c      call wri3ar(s,2,nx,2,ny,1,ns-1,'spp     ',tk,iomod,nx,ny)
c      call wri3ar(u(0,0,0,3),2,nx,2,ny,1,ns-1,'upp     ',tk,iomod,nx,ny)
c      call wri3ar(v(0,0,0,3),2,nx,2,ny,1,ns-1,'vpp     ',tk,iomod,nx,ny)
c      call wri3ar(wsig(0,0,0,3),2,nx,2,ny,1,ns-1,'wsigpp  ',tk,iomod,nx
c     : ,ny)
c      call wri3ar(w(0,0,0,3),2,nx,2,ny,1,ns-1,'wpp     ',tk,iomod,nx,ny)
c
c continuation of fields w, wsig to left
c
c
      do is=0,ns1
        do ix=1,nx1
          wsig(ix,0,is,3)=wsig(ix,1,is,3)
          w(ix,0,is,3)=w(ix,1,is,3)
        enddo
        do iy=0,ny1
          wsig(0,iy,is,3)=wsig(1,iy,is,3)
          w(0,iy,is,3)=w(1,iy,is,3)
        enddo
      enddo
c

c
c ***************************************************************
c      write(nchn1,*)'  um'
c      do 200 ix=0,nx1
c      write(nchn1,250)(um(ix,iy),iy=0,ny1)
c200   continue
c      write(nchn1,*)'  vm'
c      do 210 ix=0,nx1
c      write(nchn1,250)(vm(ix,iy),iy=0,ny1)
c210   continue
c      write(nchn1,*)'  u(ix,iy,2)'
c      do 211 ix=0,nx1
c      write(nchn1,250)(u(ix,iy,2,3),iy=0,ny1)
c211   continue
c      write(nchn1,*)'  v(ix,iy,2)'
c      do 212 ix=0,nx1
c      write(nchn1,250)(v(ix,iy,2,3),iy=0,ny1)
c212   continue
c      write(nchn1,*)'  pp10(ix,iy)'
c      do 213 ix=0,nx1
c      write(nchn1,250)(pp10(ix,iy,2),iy=0,ny1)
c213   continue
c      write(nchn1,*)'  pp01(ix,iy)'
c      do 214 ix=0,nx1
c      write(nchn1,250)(pp01(ix,iy,2),iy=0,ny1)
c214   continue
c      write(nchn1,*)'  s(ix,iy,2)'
c      do 215 ix=0,nx1
c      write(nchn1,250)(s(ix,iy,2),iy=0,ny1)
c215   continue
c250   format(1x,7e11.3)
c
c
      deallocate(um)
      deallocate(vm)
      deallocate(f)
      return
      end

      subroutine disc(uu,vv,disb,ii,jj)
c
c calculates maximum discontinuity
c
c                       max(abs(div(p_s*vm)))
c
c where p_s is gspressure and vector vm is vertically averaged wind.
c
      use alloc
      implicit real*8(a-h,o-z)

c      dimension uu(0:nx1,0:ny1,0:ns1), vv(0:nx1,0:ny1,0:ns1),
c     :          um(0:nx1,0:ny1), vm(0:nx1,0:ny1)
      dimension uu(0:nx1,0:ny1,0:ns1), vv(0:nx1,0:ny1,0:ns1)
      real(kind(0.d0)),dimension(:,:),allocatable:: um,vm

      allocate(um(0:nx1,0:ny1))
      allocate(vm(0:nx1,0:ny1))
c
c vertically integrated wind * pp
c
      do ix=0,nx1
      do iy=0,ny1
        sumu = 0.
        sumv = 0.
        do is=1,ns-1
          sumu=sumu+uu(ix,iy,is)*ds1(is+1)
          sumv=sumv+vv(ix,iy,is)*ds1(is+1)
        enddo
        um(ix,iy)=sumu*pp10(ix,iy,3)
        vm(ix,iy)=sumv*pp01(ix,iy,3)
      enddo
      enddo
c
c maximum discontinuity (continuity disbalance)
c
      disb=0.
      do ix=2,nx-1
      do iy=2,ny-1
        sumu=(um(ix,iy)-um(ix-1,iy))/dx+(vm(ix,iy)-vm(ix,iy-1))/dy
        sumu=abs(sumu)
        if(sumu.gt.disb) then
          disb=sumu
          ii=ix
          jj=iy
        endif
      enddo
      enddo
c
      deallocate(um)
      deallocate(vm)

      return
      end



      subroutine inprst
c
      use alloc

      implicit real*8(a-h,o-z)
c
c inputs data for restart.
c
      read(21) nxw,nyw,nsw,dxw,dyw,dsw1,hmw,xmw,ymw,hxwidw,hywidw,xlatw
     :   ,qifw,nstep,ifsoil
c
      if((nxw.ne.nx).or.(nyw.ne.ny).or.(nsw.ne.ns).or.(dxw.ne.dx)
     :   .or.(dyw.ne.dy).or.(dsw1.ne.ds0(1)).or.(hmw.ne.hmount)
     :   .or.(xmw.ne.xmount).or.(ymw.ne.ymount).or.(hxwidw.ne.hxwid)
     :   .or.(hywidw.ne.hywid).or.(xlatw.ne.xlatit).or.(qifw.ne.qif))
     :   then
         write(6,*) 'nx,ny,ns,dx,dy,ds,hm,xm,ym,hxwid,hywid,xlat,qif'
         write(6,*) nx,ny,ns,dx,dy,ds0(1),hmount,xmount,ymount
     :      ,hxwid,hywid,xlatit,qif
         write(6,*) 'here:'
         write(6,*) nxw,nyw,nsw,dxw,dyw,dsw1,hmw,xmw,ymw
     :      ,hxwidw,hywidw,xlatw,qifw
         write(6,*) 'cf restart parameters'
      endif
c
      read(21) ((hsuf(ix,iy),ix=0,nx+1),iy=0,ny+1)
      read(21) ((psufi(ix,iy),ix=0,nx+1),iy=0,ny+1)
      read(21) ((phisui(ix,iy),ix=0,nx+1),iy=0,ny+1)
      read(21) ((tsufi(ix,iy),ix=0,nx+1),iy=0,ny+1)
c
      read(21) ((phisuf(ix,iy),ix=0,nx+1),iy=0,ny+1)
c
      read(21) ((pp(ix,iy,4),ix=0,nx+2),iy=0,ny+2)
      read(21) ((pp(ix,iy,3),ix=0,nx+2),iy=0,ny+2)
      read(21) ((pp(ix,iy,2),ix=0,nx+2),iy=0,ny+2)
      read(21) ((pp(ix,iy,1),ix=0,nx+2),iy=0,ny+2)
      read(21) ((dpp(ix,iy),ix=0,nx+1),iy=0,ny+1)
c
      do is=0,ns+1
        read(21) ((u(ix,iy,is,3),ix=0,nx+1),iy=0,ny+1)
        read(21) ((u(ix,iy,is,2),ix=0,nx+1),iy=0,ny+1)
        read(21) ((v(ix,iy,is,3),ix=0,nx+1),iy=0,ny+1)
        read(21) ((v(ix,iy,is,2),ix=0,nx+1),iy=0,ny+1)
        read(21) ((w(ix,iy,is,3),ix=0,nx+1),iy=0,ny+1)
        read(21) ((w(ix,iy,is,2),ix=0,nx+1),iy=0,ny+1)
        read(21) ((pt(ix,iy,is,3),ix=0,nx+1),iy=0,ny+1)
        read(21) ((pt(ix,iy,is,2),ix=0,nx+1),iy=0,ny+1)
        read(21) ((wsig(ix,iy,is,2),ix=0,nx+1),iy=0,ny+1)
        read(21) ((phis(ix,iy,is),ix=0,nx+1),iy=0,ny+1)
        read(21) ((phi(ix,iy,is),ix=0,nx+1),iy=0,ny+1)
        if(qif.ne.0.) read(21) ((qv(ix,iy,is,3),ix=0,nx+1),iy=0,ny+1)
        if(qif.ne.0.) read(21) ((qv(ix,iy,is,2),ix=0,nx+1),iy=0,ny+1)
        if(qif.ne.0.) read(21) ((qv(ix,iy,is,1),ix=0,nx+1),iy=0,ny+1)
        if(ifqc.ne.0) read(21) ((qc(ix,iy,is,3),ix=0,nx+1),iy=0,ny+1)
        if(ifqc.ne.0) read(21) ((qc(ix,iy,is,2),ix=0,nx+1),iy=0,ny+1)
        if(ifqc.ne.0) read(21) ((qc(ix,iy,is,1),ix=0,nx+1),iy=0,ny+1)
        if(ifqr.ne.0) read(21) ((qr(ix,iy,is,3),ix=0,nx+1),iy=0,ny+1)
        if(ifqr.ne.0) read(21) ((qr(ix,iy,is,2),ix=0,nx+1),iy=0,ny+1)
        if(ifqr.ne.0) read(21) ((qr(ix,iy,is,1),ix=0,nx+1),iy=0,ny+1)
	  if(ifqi.ne.0) read(21) ((qci(ix,iy,is,3),ix=0,nx+1),iy=0,ny+1)          DC,11.2009
        if(ifqi.ne.0) read(21) ((qci(ix,iy,is,2),ix=0,nx+1),iy=0,ny+1)          DC,11.2009
        if(ifqi.ne.0) read(21) ((qci(ix,iy,is,1),ix=0,nx+1),iy=0,ny+1)          DC,11.2009
	  if(ifqi.ne.0) read(21) ((qsn(ix,iy,is,3),ix=0,nx+1),iy=0,ny+1)          DC,11.2009
        if(ifqi.ne.0) read(21) ((qsn(ix,iy,is,2),ix=0,nx+1),iy=0,ny+1)          DC,11.2009
        if(ifqi.ne.0) read(21) ((qsn(ix,iy,is,1),ix=0,nx+1),iy=0,ny+1)          DC,11.2009
      enddo
c
      do ib=1,2
        read(21) (ppbx(ib,iy),iy=0,ny+2)
        read(21) (ppbx2(ib,iy),iy=0,ny+2)
        read(21) (ppby(ix,ib),ix=0,nx+2)
        read(21) (ppby2(ix,ib),ix=0,nx+2)
      enddo
c
      do ib=1,2
        read(21) ((ubx(ib,iy,is),iy=0,ny+1),is=0,ns+1)
        read(21) ((ubx2(ib,iy,is),iy=0,ny+1),is=0,ns+1)
        read(21) ((uby(ix,ib,is),ix=0,nx+1),is=0,ns+1)
        read(21) ((uby2(ix,ib,is),ix=0,nx+1),is=0,ns+1)

        read(21) ((vbx(ib,iy,is),iy=0,ny+1),is=0,ns+1)
        read(21) ((vbx2(ib,iy,is),iy=0,ny+1),is=0,ns+1)
        read(21) ((vby(ix,ib,is),ix=0,nx+1),is=0,ns+1)
        read(21) ((vby2(ix,ib,is),ix=0,nx+1),is=0,ns+1)

        read(21) ((ptbx(ib,iy,is),iy=0,ny+1),is=0,ns+1)
        read(21) ((ptby(ix,ib,is),ix=0,nx+1),is=0,ns+1)
        if(qif.ne.0.) read(21) ((qvbx(ib,iy,is),iy=0,ny+1),is=0,ns+1)
        if(qif.ne.0.) read(21) ((qvby(ix,ib,is),ix=0,nx+1),is=0,ns+1)
        if(ifqc.ne.0) read(21) ((qcbx(ib,iy,is),iy=0,ny+1),is=0,ns+1)
        if(ifqc.ne.0) read(21) ((qcby(ix,ib,is),ix=0,nx+1),is=0,ns+1)
        if(ifqr.ne.0) read(21) ((qrbx(ib,iy,is),iy=0,ny+1),is=0,ns+1)
        if(ifqr.ne.0) read(21) ((qrby(ix,ib,is),ix=0,nx+1),is=0,ns+1)
	  if(ifqi.ne.0) read(21) ((qcibx(ib,iy,is),iy=0,ny+1),is=0,ns+1)          DC,11.2009
        if(ifqi.ne.0) read(21) ((qciby(ix,ib,is),ix=0,nx+1),is=0,ns+1)          DC,11.2009
	  if(ifqi.ne.0) read(21) ((qsnbx(ib,iy,is),iy=0,ny+1),is=0,ns+1)          DC,11.2009
        if(ifqi.ne.0) read(21) ((qsnby(ix,ib,is),ix=0,nx+1),is=0,ns+1)          DC,11.2009

      enddo
c
      read(21) ((ppcc(i,j),i=1,2),j=1,2)
      read(21) ((ppcc2(i,j),i=1,2),j=1,2)
      read(21) (((ucc(i,j,is),is=0,ns),j=1,2),i=1,2)
      read(21) (((ucc2(i,j,is),is=0,ns),j=1,2),i=1,2)
      read(21) (((vcc(i,j,is),is=0,ns),j=1,2),i=1,2)
      read(21) (((vcc2(i,j,is),is=0,ns),j=1,2),i=1,2)
      read(21) (((ptcc(i,j,is),is=0,ns),j=1,2),i=1,2)
      if(qif.ne.0.) read(21) (((qvcc(i,j,is),is=0,ns),j=1,2),i=1,2)
      if(ifqc.ne.0) read(21) (((qccc(i,j,is),is=0,ns),j=1,2),i=1,2)
      if(ifqr.ne.0) read(21) (((qrcc(i,j,is),is=0,ns),j=1,2),i=1,2)
	if(ifqi.ne.0) read(21) (((qcicc(i,j,is),is=0,ns),j=1,2),i=1,2)            DC,11.2009
	if(ifqi.ne.0) read(21) (((qsncc(i,j,is),is=0,ns),j=1,2),i=1,2)            DC,11.2009

      if(ifsoil.ne.0) then
         read(21) ((tsnoi(ix,iy),ix=0,nx+1),iy=0,ny+1)
         read(21) ((t2noi(ix,iy),ix=0,nx+1),iy=0,ny+1)
         read(21) ((wgnoi(ix,iy),ix=0,nx+1),iy=0,ny+1)
         read(21) ((w2noi(ix,iy),ix=0,nx+1),iy=0,ny+1)
         read(21) ((wrnoi(ix,iy),ix=0,nx+1),iy=0,ny+1)
         read(21) ((h(ix,iy),ix=0,nx+1),iy=0,ny+1)
         read(21) ((le(ix,iy),ix=0,nx+1),iy=0,ny+1)
         read(21) ((gsolo(ix,iy),ix=0,nx+1),iy=0,ny+1)
         read(21) ((rn(ix,iy),ix=0,nx+1),iy=0,ny+1)
      endif
c
      close(21)
c
      do iy=0,ny+1
      do ix=0,nx+1
        psuf(ix,iy)=pp(ix,iy,3)+ptop
      enddo
      enddo
c
      return
      end

      subroutine newrefs(zup)

c-----------------------------------------------------------------------
c specification of reference state
c new version: if fcor.ne.0 the initial state is defined using
c geostrophic balance.
c-----------------------------------------------------------------------
c
      use alloc

      implicit real*8(a-h,o-z)
      dimension zup(0:nx+1,0:ny+1)

c
c if fcor.ne.0 interpolate surface pressure from reference profile
c assuming constant temperature gradient.
c
c     do iy=0,ny+1
c     do ix=0,nx+1
c       psuf(ix,iy)=psuf(ix,iy)*dexp(-g*zup(ix,iy)/(r*tsuf(ix,iy)))
c     enddo
c     enddo

      if(fcor.ne.0..and.itepsuf) then
        if(ug.ne.0. .or. vg.ne.0.) then
          if(iophis.eq.3) then
            ixref=nx1/2 !3
            iyref=ny1/2
            pref=psuf(ixref,iyref)
            call surpsu(ixref,iyref,pref)
          else
            call inpos2(iobphy)
            call surpsd
          endif
        endif
      endif

      do iy=1,ny1
      do ix=1,nx1
        pp(ix,iy,3)=psuf(ix,iy)-ptop
      enddo
      enddo

      call extrpp(nx2,ny2,pp(0,0,3),0,nx2,0,ny2)
      do iy=0,ny1
      do ix=0,nx1
        dpdx10(ix,iy,3)=(pp(ix+1,iy,3)-pp(ix,iy,3))/dx
        dpdy01(ix,iy,3)=(pp(ix,iy+1,3)-pp(ix,iy,3))/dy
      enddo
      enddo
c
c interpolate pts in the vertical, according to present values of phis:
c
      call upuvt
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
      endif
c
c calculates phis integrating hydrostatic relation:
c
      do iy=0,ny1
      do ix=0,nx1
        phis(ix,iy,ns)=phisuf(ix,iy)-dsng*r*pp(ix,iy,3)
     :    *(.75*tems(ix,iy,ns)+.25*tems(ix,iy,ns-1))
     :    /(ptop+pp(ix,iy,3)*(1.+0.5*dsng))
      enddo
      enddo
c
      do is=ns-1,0,-1
      do iy=0,ny1
      do ix=0,nx1
        phis(ix,iy,is)=phis(ix,iy,is+1)+ds0(is+1)*r05*pp(ix,iy,3)
     :    *(tems(ix,iy,is+1)+tems(ix,iy,is))
     :    /(ptop+pp(ix,iy,3)*sigma1(is+1))
      enddo
      enddo
      enddo

c     do k=1,4
c     do iy=0,ny2
c     do ix=0,nx2
c       pp(ix,iy,k)=pp(ix,iy,3)
c     enddo
c     enddo
c     enddo
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
      do k=3,3
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

      return
      end

      subroutine out2(klev,fnout)
c
c unformatted output
c
      use alloc
      use refstate

      implicit real*8(a-h,o-z)
      character*80 fnout,fext
      integer jstep
      save jstep
c      include 'nh3dpa10.inc'


      save nrec,nrecp,nrecs
      data jstep/0/
      data nrec/0/,nrecp/1/,nrecs/0/,nrecprec/0/
c

c
      write(*,'(1x,a8,i3,a40)') 'out2:',jstep,fnout
      if(jstep.eq.0) then
        open(78,file=fext(fnout,'lgr'))
        write(78,'(1x,a40)') fnout
        write(78,'(1x,4i7,4f17.2)') ntime,nx,ny,ns,dx,dy,timeini,dt
        jstep=jstep+1
        write(78,'(i10)') nstep
      else
        jstep=jstep+1
        write(78,'(i10)') nstep
      endif





      write(30,1000) nstep
c      write(*,*) nstep,thdat(0),thdat(29)
      write(36,rec=nrecp) ndat,ndth,ndus,ndvs,ndqvs,nprof
     :  ,((zthdat(i,ip),thdat(i,ip),zusdat(i,ip),usdat(i,ip)
     :   ,zvsdat(i,ip),vsdat(i,ip),zqvsdat(i,ip),qvsdat(i,ip),i=1,ndat)
     :  ,ip=1,nprof)
c
      do is=0,ns
        nrec=nrec+1
        write(11,rec=nrec) ((u(ix,iy,is,klev),ix=0,nx1),iy=0,ny1)
        write(12,rec=nrec) ((v(ix,iy,is,klev),ix=0,nx1),iy=0,ny1)
        write(13,rec=nrec) ((w(ix,iy,is,klev),ix=0,nx1),iy=0,ny1)
        write(14,rec=nrec) ((pt(ix,iy,is,klev),ix=0,nx1),iy=0,ny1)
        if(qif.ne.0.) then
          write(20,rec=nrec) ((qv(ix,iy,is,klev),ix=0,nx1),iy=0,ny1)
        endif
        if(ifqc.ne.0) then
          write(37,rec=nrec) ((qc(ix,iy,is,klev),ix=0,nx1),iy=0,ny1)
          write(39,rec=nrec) ((cond(ix,iy,is),ix=0,nx1),iy=0,ny1)
          write(56,rec=nrec) ((evap(ix,iy,is),ix=0,nx1),iy=0,ny1)
        endif
        if(ifqr.ne.0) then
          write(38,rec=nrec) ((qr(ix,iy,is,klev),ix=0,nx1),iy=0,ny1)
          write(54,rec=nrec) ((auto(ix,iy,is),ix=0,nx1),iy=0,ny1)
          write(55,rec=nrec) ((col(ix,iy,is),ix=0,nx1),iy=0,ny1)
          write(57,rec=nrec) ((vrain(ix,iy,is),ix=0,nx1),iy=0,ny1)
        endif
        write(15,rec=nrec) ((phi(ix,iy,is),ix=0,nx1),iy=0,ny1)
        write(16,rec=nrec) ((phis(ix,iy,is),ix=0,nx1),iy=0,ny1)
c        write(17,rec=nrec) ((wsig(ix,iy,is,klev),ix=0,nx1),iy=0,ny1)
        if(iodif.ne.0) then
          write(18,rec=nrec) ((difunu(ix,iy,is)/pp10(ix,iy,klev)
     :      ,ix=0,nx1),iy=0,ny1)
          write(19,rec=nrec) ((difunv(ix,iy,is)/pp01(ix,iy,klev)
     :      ,ix=0,nx1),iy=0,ny1)
          write(47,rec=nrec) ((difunt(ix,iy,is),ix=0,nx1),iy=0,ny1)
        endif
      enddo
      nrecp=nrecp+1
      nrecprec=nrecprec+1
      if(ifqr.ne.0) then
        write(58,rec=nrecprec) ((prec(ix,iy),ix=0,nx1),iy=0,ny1)
      endif
      write(10,rec=nrecp) ((psuf(ix,iy),ix=0,nx1),iy=0,ny1)
      if(ifsoil.ne.0) then
         nrecs=nrecs+1
         write(33,rec=nrecs) (((1.-xlake(ix,iy))*tsnoi(ix,iy)+
     :     xlake(ix,iy)*tslake(ix,iy)
     :     ,ix=0,nx1),iy=0,ny1)
         nrecs=nrecs+1
         write(33,rec=nrecs) ((t2noi(ix,iy),ix=0,nx1),iy=0,ny1)
         nrecs=nrecs+1
         write(33,rec=nrecs) ((wgnoi(ix,iy),ix=0,nx1),iy=0,ny1)
         nrecs=nrecs+1
         write(33,rec=nrecs) ((w2noi(ix,iy),ix=0,nx1),iy=0,ny1)
         nrecs=nrecs+1
         write(33,rec=nrecs) ((wrnoi(ix,iy),ix=0,nx1),iy=0,ny1)
         nrecs=nrecs+1
         write(33,rec=nrecs) ((h(ix,iy),ix=0,nx1),iy=0,ny1)
         nrecs=nrecs+1
         write(33,rec=nrecs) ((le(ix,iy),ix=0,nx1),iy=0,ny1)
         nrecs=nrecs+1
         write(33,rec=nrecs) ((gsolo(ix,iy),ix=0,nx1),iy=0,ny1)
         nrecs=nrecs+1
         write(33,rec=nrecs) ((rn(ix,iy),ix=0,nx1),iy=0,ny1)
         write(34,rec=nrecp-1) ((cdm(ix,iy),ix=0,nx1),iy=0,ny1)
      endif
c
      return
1000  format(1x,i5)
      end



      
      subroutine outdat(fnout)
c
c outputs model parameters:
c
      use alloc
      use refstate
      implicit real*8(a-h,o-z)
c      include 'nh3dpa10.inc'
      character*80 fnout
c
      write(30,*) 'experiment: ',exptit
c
      write(30,*) 'ntime'
      write(30,*) ntime
c
      write(30,*) 'nx,ny,ns,dx,dy,dt,pa,ptop'
      write(30,1000) nx,ny,ns,dx,dy,dt,pa,ptop
c
      write(30,*) 'sigma0,sigma1'
      write(30,1001) (sigma0(is),sigma1(is),is=0,ns1)
c
      write(30,*) 'hmount,xmount,ymount,hxwid,hywid,iomtyp'
      write(30,1002) hmount,xmount,ymount,hxwid,hywid,iomtyp
c
      write(10,rec=1) ((hsuf(ix,iy),ix=0,nx1),iy=0,ny1),0.
c
      write(30,*) ' nuppts,ioreft,pts0,iorefq,nupsur,iophis
     :   ,iorefu,iorefv'
      write(30,1003) nuppts,ioreft,pts0,iorefq,nupsur,iophis
     :   ,iorefu,iorefv
c
      write(30,*) 'ndat,ndth,ndus,ndvs,ndqvs'
      write(30,1004) ndat,ndth,ndus,ndvs,ndqvs,nprof
c
      do ip=1,nprof
      write(30,*)'zthdat,thdat,zusdat,usdat,zvsdat,vsdat,zqvsdat,qvsdat'
      do i=0,ndat
        write(30,1005) zthdat(i,ip),thdat(i,ip),zusdat(i,ip),usdat(i,ip)
     :    ,zvsdat(i,ip),vsdat(i,ip),zqvsdat(i,ip),qvsdat(i,ip)
      enddo
      enddo
c
      write(30,*) 'xlatit,fcor'
      write(30,1006) xlatit,fcor
c
      write(30,*) 'qif,qdrag'
      write(30,1007) qif,qdrag
c
      write(30,*) 'iodif,difcof,difl,uvdif,tsdif,qsdif,nupdif'
      write(30,1008) iodif,difcof,difl,uvdif,tsdif,qsdif,nupdif
c
      write(30,*) 'taut,taum,dkapa,rikey,zm,zt'
      write(30,1009) taut,taum,dkapa,rikey,zm,zt
c
      write(30,*) 'hydro,tfct,raylei,olddra'
      write(30,*) hydro,tfct,raylei,olddra
c
      write(30,*) 'bc.:uux,uuy,vvx,vvy,ptx,pty,ppx,ppy,phy,dfy,puv'
      write(30,1010) iobuux,iobuuy,iobvvx,iobvvy,iobptx,iobpty
     :   ,iobppx,iobppy,iobphy,iodify,ioppuv
c
      return
1000  format(1x,3i5,2f10.1,f6.0,2f8.0)
1001  format(1x,2f12.7)
1002  format(1x,5f10.1,i3)
1003  format(1x,2i7,f8.2,5i5)
1004  format(1x,5i5)
1005  format(1x,f8.0,e12.5,3(f8.0,f10.3))
1006  format(1x,f6.1,e14.7)
1007  format(1x,2f4.1)
1008  format(1x,i2,5f6.1,i3)
1009  format(1x,6e14.7)
1010  format(3x,11i4)
      end



      subroutine outrst(fnout)
c
      use alloc

      implicit real*8(a-h,o-z)
      character*80 fnout,fext
c      include 'nh3dpa10.inc'

c
c outputs data for restart.
c
c     save norst
c     data norst/0/
c
c     norst=norst+1
c
c     n40=40+norst
      n40=41

      open(n40,file=fext(fnout,'rst'),form='unformatted'
     :  ,status='unknown')
c
c
      write(n40) nx,ny,ns,dx,dy,ds0(1),hmount,xmount,ymount,hxwid,hywid
     :   ,xlatit,qif,nstep,ifsoil
c
      write(n40) ((hsuf(ix,iy),ix=0,nx+1),iy=0,ny+1)
      write(n40) ((psufi(ix,iy),ix=0,nx+1),iy=0,ny+1)
      write(n40) ((phisui(ix,iy),ix=0,nx+1),iy=0,ny+1)
      write(n40) ((tsufi(ix,iy),ix=0,nx+1),iy=0,ny+1)
c
      write(n40) ((phisuf(ix,iy),ix=0,nx+1),iy=0,ny+1)
c
      write(n40) ((pp(ix,iy,4),ix=0,nx+2),iy=0,ny+2)
      write(n40) ((pp(ix,iy,3),ix=0,nx+2),iy=0,ny+2)
      write(n40) ((pp(ix,iy,2),ix=0,nx+2),iy=0,ny+2)
      write(n40) ((pp(ix,iy,1),ix=0,nx+2),iy=0,ny+2)
      write(n40) ((dpp(ix,iy),ix=0,nx+1),iy=0,ny+1)
c
      do is=0,ns+1
        write(n40) ((u(ix,iy,is,3),ix=0,nx+1),iy=0,ny+1)
        write(n40) ((u(ix,iy,is,2),ix=0,nx+1),iy=0,ny+1)
        write(n40) ((v(ix,iy,is,3),ix=0,nx+1),iy=0,ny+1)
        write(n40) ((v(ix,iy,is,2),ix=0,nx+1),iy=0,ny+1)
        write(n40) ((w(ix,iy,is,3),ix=0,nx+1),iy=0,ny+1)
        write(n40) ((w(ix,iy,is,2),ix=0,nx+1),iy=0,ny+1)
        write(n40) ((pt(ix,iy,is,3),ix=0,nx+1),iy=0,ny+1)
        write(n40) ((pt(ix,iy,is,2),ix=0,nx+1),iy=0,ny+1)
        write(n40) ((wsig(ix,iy,is,2),ix=0,nx+1),iy=0,ny+1)
        write(n40) ((phis(ix,iy,is),ix=0,nx+1),iy=0,ny+1)
        write(n40) ((phi(ix,iy,is),ix=0,nx+1),iy=0,ny+1)
        if(qif.ne.0.) write(n40) ((qv(ix,iy,is,3),ix=0,nx+1),iy=0,ny+1)
        if(qif.ne.0.) write(n40) ((qv(ix,iy,is,2),ix=0,nx+1),iy=0,ny+1)
        if(qif.ne.0.) write(n40) ((qv(ix,iy,is,1),ix=0,nx+1),iy=0,ny+1)
        if(ifqc.ne.0) write(n40) ((qc(ix,iy,is,3),ix=0,nx+1),iy=0,ny+1)
        if(ifqc.ne.0) write(n40) ((qc(ix,iy,is,2),ix=0,nx+1),iy=0,ny+1)
        if(ifqc.ne.0) write(n40) ((qc(ix,iy,is,1),ix=0,nx+1),iy=0,ny+1)
        if(ifqr.ne.0) write(n40) ((qr(ix,iy,is,3),ix=0,nx+1),iy=0,ny+1)
        if(ifqr.ne.0) write(n40) ((qr(ix,iy,is,2),ix=0,nx+1),iy=0,ny+1)
        if(ifqr.ne.0) write(n40) ((qr(ix,iy,is,1),ix=0,nx+1),iy=0,ny+1)
	  if(ifqi.ne.0) write(n40) ((qci(ix,iy,is,3),ix=0,nx+1),iy=0,ny+1)        DC,11.2009
        if(ifqi.ne.0) write(n40) ((qci(ix,iy,is,2),ix=0,nx+1),iy=0,ny+1)        DC,11.2009
        if(ifqi.ne.0) write(n40) ((qci(ix,iy,is,1),ix=0,nx+1),iy=0,ny+1)        DC,11.2009
	  if(ifqi.ne.0) write(n40) ((qsn(ix,iy,is,3),ix=0,nx+1),iy=0,ny+1)        DC,11.2009
        if(ifqi.ne.0) write(n40) ((qsn(ix,iy,is,2),ix=0,nx+1),iy=0,ny+1)        DC,11.2009
        if(ifqi.ne.0) write(n40) ((qsn(ix,iy,is,1),ix=0,nx+1),iy=0,ny+1)        DC,11.2009
      enddo
c
      do ib=1,2
        write(n40) (ppbx(ib,iy),iy=0,ny+2)
        write(n40) (ppbx2(ib,iy),iy=0,ny+2)
        write(n40) (ppby(ix,ib),ix=0,nx+2)
        write(n40) (ppby2(ix,ib),ix=0,nx+2)
      enddo
c
      do ib=1,2
        write(n40) ((ubx(ib,iy,is),iy=0,ny+1),is=0,ns+1)
        write(n40) ((ubx2(ib,iy,is),iy=0,ny+1),is=0,ns+1)
        write(n40) ((uby(ix,ib,is),ix=0,nx+1),is=0,ns+1)
        write(n40) ((uby2(ix,ib,is),ix=0,nx+1),is=0,ns+1)
c
        write(n40) ((vbx(ib,iy,is),iy=0,ny+1),is=0,ns+1)
        write(n40) ((vbx2(ib,iy,is),iy=0,ny+1),is=0,ns+1)
        write(n40) ((vby(ix,ib,is),ix=0,nx+1),is=0,ns+1)
        write(n40) ((vby2(ix,ib,is),ix=0,nx+1),is=0,ns+1)
c
        write(n40) ((ptbx(ib,iy,is),iy=0,ny+1),is=0,ns+1)
        write(n40) ((ptby(ix,ib,is),ix=0,nx+1),is=0,ns+1)
        if(qif.ne.0.) write(n40) ((qvbx(ib,iy,is),iy=0,ny+1),is=0,ns+1)
        if(qif.ne.0.) write(n40) ((qvby(ix,ib,is),ix=0,nx+1),is=0,ns+1)
        if(ifqc.ne.0) write(n40) ((qcbx(ib,iy,is),iy=0,ny+1),is=0,ns+1)
        if(ifqc.ne.0) write(n40) ((qcby(ix,ib,is),ix=0,nx+1),is=0,ns+1)
        if(ifqr.ne.0) write(n40) ((qrbx(ib,iy,is),iy=0,ny+1),is=0,ns+1)
        if(ifqr.ne.0) write(n40) ((qrby(ix,ib,is),ix=0,nx+1),is=0,ns+1)
	  if(ifqi.ne.0) write(n40) ((qcibx(ib,iy,is),iy=0,ny+1),is=0,ns+1)        DC,11.2009
        if(ifqi.ne.0) write(n40) ((qciby(ix,ib,is),ix=0,nx+1),is=0,ns+1)        DC,11.2009
	  if(ifqi.ne.0) write(n40) ((qsnbx(ib,iy,is),iy=0,ny+1),is=0,ns+1)        DC,11.2009
        if(ifqi.ne.0) write(n40) ((qsnby(ix,ib,is),ix=0,nx+1),is=0,ns+1)        DC,11.2009
      enddo
c
      write(n40) ((ppcc(i,j),i=1,2),j=1,2)
      write(n40) ((ppcc2(i,j),i=1,2),j=1,2)
      write(n40) (((ucc(i,j,is),is=0,ns),j=1,2),i=1,2)
      write(n40) (((ucc2(i,j,is),is=0,ns),j=1,2),i=1,2)
      write(n40) (((vcc(i,j,is),is=0,ns),j=1,2),i=1,2)
      write(n40) (((vcc2(i,j,is),is=0,ns),j=1,2),i=1,2)
      write(n40) (((ptcc(i,j,is),is=0,ns),j=1,2),i=1,2)
      if(qif.ne.0.) write(n40) (((qvcc(i,j,is),is=0,ns),j=1,2),i=1,2)
      if(ifqc.ne.0) write(n40) (((qccc(i,j,is),is=0,ns),j=1,2),i=1,2)
      if(ifqr.ne.0) write(n40) (((qrcc(i,j,is),is=0,ns),j=1,2),i=1,2)
	if(ifqi.ne.0) write(n40) (((qcicc(i,j,is),is=0,ns),j=1,2),i=1,2)          DC,11.2009
	if(ifqi.ne.0) write(n40) (((qsncc(i,j,is),is=0,ns),j=1,2),i=1,2)          DC,11.2009
      if(ifsoil.ne.0) then
        write(n40) ((tsnoi(ix,iy),ix=0,nx+1),iy=0,ny+1)
        write(n40) ((t2noi(ix,iy),ix=0,nx+1),iy=0,ny+1)
        write(n40) ((wgnoi(ix,iy),ix=0,nx+1),iy=0,ny+1)
        write(n40) ((w2noi(ix,iy),ix=0,nx+1),iy=0,ny+1)
        write(n40) ((wrnoi(ix,iy),ix=0,nx+1),iy=0,ny+1)
        write(n40) ((h(ix,iy),ix=0,nx+1),iy=0,ny+1)
        write(n40) ((le(ix,iy),ix=0,nx+1),iy=0,ny+1)
        write(n40) ((gsolo(ix,iy),ix=0,nx+1),iy=0,ny+1)
        write(n40) ((rn(ix,iy),ix=0,nx+1),iy=0,ny+1)
      endif
c
      close(n40)
c
      return
      end


      subroutine rain
c
c Evaluation of the liquid precipitation velocity (rain water->soil)
c
c vrain=const*rho**(-3/8)*qr**(1/8)
c
c after Emanuel(1994)
c
      use alloc

      implicit real*8(a-h,o-z)
c      include 'nh3dpa10.inc'


      do is=0,ns
        do iy=0,ny1
          do ix=0,nx1
!            p=sigma0(is)*pp(ix,iy,2)+ptop
!            rho=p/(r*(pts(ix,iy,is)+pt(ix,iy,is,2))*(p/p00)**akapa
!     :       *(1+0.61*qv(ix,iy,is,2)-qc(ix,iy,is,2)))
            xlr=(rho(ix,iy,is)*qr(ix,iy,is,2)/(pi*rhol*xnor))**0.25
            vrain(ix,iy,is)=xa*gm48/6.*xlr**xb*sqrt(rho0/rho(ix,iy,is))
          enddo
        enddo
      enddo

      is=ns
      do iy=0,ny1
        do ix=0,nx1
!          p=sigma0(is)*pp(ix,iy,2)+ptop
!          rho=p/(r*(pts(ix,iy,is)
!     :    +pt(ix,iy,is,2))*(p/p00)**akapa*(1+0.61*qv(ix,iy,is,2)
!     :    -qc(ix,iy,is,2)))
          prec(ix,iy)=(0.5*(vrain(ix,iy,is)
     :      +vrain(ix,iy,is-1))-w(ix,iy,is,2))
     :      *rho(ix,iy,is)*qr(ix,iy,is,2)
        enddo
      enddo
      rainacc=rainacc+prec
      return
      end



      subroutine readpa(fnmap,fnfor,fnout
     :   ,itheta,iwind,scrout,fngrd,dif3d,diftx,difloc,nonloc_tm86,
     : nonloc_ls96)
c
c Input parameters, data and option switchs
c
      use alloc
      use nh3dpa2
      use refstate
      implicit real*8(a-h,o-z)

      character*80 fnrst,fnmap,fnfor,fnout,fext,fngrd,fnrad
      character*80 line,form1
      character*80 fileext,lstname
      character*14 timestring
      character*20 keyword,uppercase
      logical scrout,autogrid,anamount,fluxsel,new_prof,new_point
      external fext,uppercase
      logical ok_nx,ok_ny,ok_ns,ok_dx,ok_dy,ok_dt,ok_ntime,ok_profile
     :  ,ok_fnmap,ok_fnout,ok_fngrd,ok_fnrst,ok_fnfor,ok_plan
     :  ,ok_rad,dif3d, diftx
      logical difloc,nonloc_tm86,nonloc_ls96                                            DC 06.2010
c
c      include 'nh3dpa10.inc'
c
      ifhle=0
      taum=0.
      tauspon=0.
      taut=0.
      autogrid=.true.
      anamount=.false.
      fluxsel=.false.
      new_prof=.false.
      new_point=.false.
      npro=400
      dzpro=50.
      fcor=-999.
      open(29,file='nh3d.dat',status='old')
c     open(99,file='deltaz.dat',status='unknown')
      read(29,'(a)') line
c
c check if new input file format (nh3d9715 onwards)
c
      if(uppercase(line(1:4)).ne.'NH3D') then
        write(0,*) 'nh3dFatalError=Wrong Input file nh3d.dat !'
        stop
      endif
c
c define default values for all adjustable parameters
c
      call defaults(fnmap,fnfor,fnout
     :   ,itheta,iwind,scrout,fngrd,nbyte_acc)

c
      open(27,status='scratch')

      timeini0=timeini
      ok_nx=.false.
      ok_ny=.false.
      ok_ns=.false.
      ok_dx=.false.
      ok_dy=.false.
      ok_dt=.false.
      ok_ntime=.false.
      ok_profile=.false.
      ok_fnmap=.false.
      ok_fnout=.false.
      ok_fngrd=.false.
      ok_fnrst=.false.
      ok_fnfor=.false.
      ok_rad=.false.
      icha=ichar('0')
      ichz=ichar('z')

      itext=0
      inperr=0
      do while(inperr.eq.0)
        read(29,'(a)') line
c        write(*,*) line
c
c isolate keyword
c
        ichlef=1
        do i=1,80
          if(line(i:i).ne.' ') exit
          ichlef=i
        enddo
        if(ichlef.eq.80) cycle
        ichrig=ichlef
        do i=ichlef+1,80
          if(ichar(line(i:i)).lt.icha.or.ichar(line(i:i)).gt.ichz) exit
          ichrig=i
        enddo
        keyword=line(ichlef:ichrig)
        nchark=ichrig-ichlef+1
        call toupcase(keyword)
        if(verbose.ge.2) write(*,*) keyword,nchark,ichlef,ichrig
        write(27,'(a)') line(ichrig+1:80)
c       write(*,'(a)') line(ichrig+1:80)
        backspace(27)
        if(keyword(1:3).eq.uppercase('end')) then
          exit
        elseif(keyword(1:1).eq.'#') then
          cycle
        elseif(keyword(1:nchark).eq.uppercase('verbose')) then
          read(27,*,iostat=inperr) verbose
          if(verbose.ge.1) write(*,*) 'verbose=',verbose
        elseif(keyword(1:nchark).eq.uppercase('nx')) then
          read(27,*,iostat=inperr) nx
          if(verbose.ge.0) write(*,*) 'nx=',nx
          ok_nx=.true.
        elseif(keyword(1:nchark).eq.uppercase('ny')) then
          read(27,*,iostat=inperr) ny
          if(verbose.ge.0) write(*,*) 'ny=',ny
          ok_ny=.true.
        elseif(keyword(1:nchark).eq.uppercase('ns')) then
          read(27,*,iostat=inperr) ns
          if(verbose.ge.0) write(*,*) 'ns=',ns
          ok_ns=.true.
        elseif(keyword(1:nchark).eq.uppercase('isqv')) then
          read(27,*,iostat=inperr) isqv
          if(verbose.ge.0) write(*,*) 'isqv=',isqv
          qif=isqv
          write(*,*) 'qif=>',qif
        elseif(keyword(1:nchark).eq.uppercase('isqc')) then
          read(27,*,iostat=inperr) isqc
          if(verbose.ge.0) write(*,*) 'isqc=',isqc
          ifqc=isqc
          write(*,*) 'ifqc=>',ifqc
        elseif(keyword(1:nchark).eq.uppercase('isqr')) then
          read(27,*,iostat=inperr) isqr
          if(verbose.ge.0) write(*,*) 'isqr=',isqr
          ifqr=isqr
          write(*,*) 'ifqr=>',ifqr
	  elseif(keyword(1:nchark).eq.uppercase('isqi')) then                     DC,10.2009
          read(27,*,iostat=inperr) isqi                                         DC,10.2009
          if(verbose.ge.0) write(*,*) 'isqi=',isqi                              DC,10.2009
          ifqi=isqi                                                             DC,10.2009
          write(*,*) 'ifqi=>',ifqi                                              DC,10.2009

c       Dima@ 020507          
        elseif(keyword(1:nchark).eq.uppercase('impur')) then                    DM,05.2007  
          read(27,*,iostat=inperr) impur                                        DM,05.2007
          if (impur.NE.0) write(*,*) 'impur =',impur,' (Impurity is ON)'        DM,05.2007 
        
        elseif(keyword(1:nchark).eq.uppercase('dx')) then
          read(27,*,iostat=inperr) dx
          if(verbose.ge.0) write(*,*) 'dx=',dx
          ok_dx=.true.
        elseif(keyword(1:nchark).eq.uppercase('dy')) then
          read(27,*,iostat=inperr) dy
          if(verbose.ge.0) write(*,*) 'dy=',dy
          ok_dy=.true.
        elseif(keyword(1:nchark).eq.uppercase('scrout')) then
          read(27,*,iostat=inperr) scrout
          if(verbose.ge.1) write(*,*) 'scrout=',scrout
        elseif(keyword(1:nchark).eq.uppercase('iotscrout')) then
          read(27,*,iostat=inperr) iotscrout
          if(verbose.ge.1) write(*,*) 'iotscrout=',iotscrout
        elseif(keyword(1:nchark).eq.uppercase('dt')) then
          read(27,*,iostat=inperr) dt
          if(verbose.ge.0) write(*,*) 'dt=',dt
          ok_dt=.true.
        elseif(keyword(1:nchark).eq.uppercase('ntime')) then
          read(27,*,iostat=inperr) ntime
          if(verbose.ge.0) write(*,*) 'ntime=',ntime
          ok_ntime=.true.
        elseif(keyword(1:nchark).eq.uppercase('timeini')) then
          read(27,*,iostat=inperr) timeini
          if(verbose.ge.1) write(*,*) 'timeini=',timeini,inperr
        elseif(keyword(1:nchark).eq.uppercase('MOUNTAIN')) then
          anamount=.true.
          read(29,*,iostat=inperr) iomtyp,hmount,xmount,ymount
     :      ,hxwid,hywid,hexp
	    write(*,'(i2,6f10.2)') iomtyp,hmount,xmount,ymount
     :       ,hxwid,hywid,hexp
          if(verbose.ge.1) then
            write(*,'(i2,6f10.2)') iomtyp,hmount,xmount,ymount
     :       ,hxwid,hywid,hexp
          endif
        elseif(keyword(1:nchark).eq.uppercase('ivarhs')) then
          read(27,*,iostat=inperr) ivarhs
          if(verbose.ge.1) write(*,*) 'ivarhs=',ivarhs
        elseif(keyword(1:nchark).eq.uppercase('nstepgrow')) then
          read(27,*,iostat=inperr) nstepgrow
          if(verbose.ge.1) write(*,*) 'nstepgrow=',nstepgrow
        elseif(keyword(1:nchark).eq.uppercase('nstepflat')) then
          read(27,*,iostat=inperr) nstepflat
          if(verbose.ge.1) write(*,*) 'nstepflat=',nstepflat
        elseif(keyword(1:nchark).eq.uppercase('numsteps')) then
          read(27,*,iostat=inperr) numsteps
          if(verbose.ge.1) write(*,*) 'numsteps=',numsteps
        elseif(keyword(1:nchark).eq.uppercase('xlatit')) then
          read(27,*,iostat=inperr) xlatit
          if(verbose.ge.1) write(*,*) 'xlatit=',xlatit
        elseif(keyword(1:nchark).eq.uppercase('fcor')) then
          read(27,*,iostat=inperr) fcor
          if(verbose.ge.1) write(*,*) 'fcor=',fcor
        elseif(keyword(1:nchark).eq.uppercase('fcori')) then
          read(27,*,iostat=inperr) fcori
          if(verbose.ge.1) write(*,*) 'fcori=',fcori
        elseif(keyword(1:nchark).eq.uppercase('ptop')) then
          read(27,*,iostat=inperr) ptop
          if(verbose.ge.1) write(*,*) 'ptop=',ptop
        elseif(keyword(1:nchark).eq.uppercase('pa')) then
          read(27,*,iostat=inperr) pa
          if(verbose.ge.1) write(*,*) 'pa=',pa
	         pressdat(0,1)=pa
        elseif(keyword(1:nchark).eq.uppercase('radpar')) then
          read(27,*,iostat=inperr) radpar
          if(verbose.ge.1) write(*,*) 'radpar=',radpar
        elseif(keyword(1:nchark).eq.uppercase('shortwave')) then                VS,26.06.2007
          read(27,*,iostat=inperr) shortwave                                    VS,26.06.2007 
          if(verbose.ge.1) write(*,*) 'shortwave=',shortwave                    VS,26.06.2007
        elseif(keyword(1:nchark).eq.uppercase('longwave')) then                 VS,26.06.2007
          read(27,*,iostat=inperr) longwave                                     VS,26.06.2007 
          if(verbose.ge.1) write(*,*) 'longwave=',longwave                      VS,26.06.2007
        elseif(keyword(1:nchark).eq.uppercase('rat_par')) then
          read(27,*,iostat=inperr) rat_par
          if(verbose.ge.1) write(*,*) 'rat_par=',rat_par  
        elseif(keyword(1:nchark).eq.uppercase('nradcall')) then
          read(27,*,iostat=inperr) nradcall
          if(verbose.ge.1) write(*,*) 'nradcall=', nradcall
        elseif(keyword(1:nchark).eq.uppercase('tau_rg')) then
          read(27,*,iostat=inperr) tau_rg
          if(verbose.ge.1) write(*,*) 'tau_rg=',tau_rg        
        elseif(keyword(1:nchark).eq.uppercase('month')) then
          read(27,*,iostat=inperr) month
          if(verbose.ge.1) write(*,*) 'month=',month
        elseif(keyword(1:nchark).eq.uppercase('day')) then
          read(27,*,iostat=inperr) day
          if(verbose.ge.1) write(*,*) 'day=', day                  
        elseif(keyword(1:nchark).eq.uppercase('iflake')) then
          read(27,*,iostat=inperr) iflake
          if(verbose.ge.1) write(*,*) 'iflake=', iflake
        elseif(keyword(1:nchark).eq.uppercase('nlakecall')) then
          read(27,*,iostat=inperr) nlakecall
          if(verbose.ge.1) write(*,*) 'nlakecall=', nlakecall  
        elseif(keyword(1:nchark).eq.uppercase('nupsur')) then
          read(27,*,iostat=inperr) nupsur
          if(verbose.ge.1) write(*,*) 'nupsur=',nupsur
        elseif(keyword(1:nchark).eq.uppercase('nuppts')) then
          read(27,*,iostat=inperr) nuppts
          if(verbose.ge.1) write(*,*) 'nuppts=',nuppts
        elseif(keyword(1:nchark).eq.uppercase('nupdif')) then
          read(27,*,iostat=inperr) nupdif
          if(verbose.ge.1) write(*,*) 'nupdif=',nupdif
        elseif(keyword(1:nchark).eq.uppercase('hydro')) then
          read(27,*,iostat=inperr) hydro
          if(verbose.ge.1) write(*,*) 'hydro=',hydro
        elseif(keyword(1:nchark).eq.uppercase('adjust')) then
          read(27,*,iostat=inperr) adjust
          if(verbose.ge.1) write(*,*) 'adjust=',adjust
        elseif(keyword(1:nchark).eq.uppercase('tfct')) then
          read(27,*,iostat=inperr) tfct
          if(verbose.ge.1) write(*,*) 'tfct=',tfct
        elseif(keyword(1:nchark).eq.uppercase('centered')) then
          read(27,*,iostat=inperr) centered
          if(verbose.ge.1) write(*,*) 'centered=',centered
        elseif(keyword(1:nchark).eq.uppercase('dst3')) then
          read(27,*,iostat=inperr) dst3
          if(verbose.ge.1) write(*,*) 'dst3=',dst3
        elseif(keyword(1:nchark).eq.uppercase('raylei')) then
          read(27,*,iostat=inperr) raylei
          if(verbose.ge.1) write(*,*) 'raylei=',raylei
        elseif(keyword(1:nchark).eq.uppercase('olddra')) then
          read(27,*,iostat=inperr) olddra
          if(verbose.ge.1) write(*,*) 'olddra=',olddra
        elseif(keyword(1:nchark).eq.uppercase('fluxxlcor')) then
          read(27,*,iostat=inperr) fluxxlcor
          if(verbose.ge.1) write(*,*) 'fluxxlcor=',fluxxlcor
          if(fluxxlcor) then
            xlcor=1.
          else
            xlcor=0.
          endif
        elseif(keyword(1:nchark).eq.uppercase('fluxylcor')) then
          read(27,*,iostat=inperr) fluxylcor
          if(verbose.ge.1) write(*,*) 'fluxylcor=',fluxylcor
          if(fluxylcor) then
            ylcor=1.
          else
            ylcor=0.
          endif
        elseif(keyword(1:nchark).eq.uppercase('fluxxrcor')) then
          read(27,*,iostat=inperr) fluxxrcor
          if(verbose.ge.1) write(*,*) 'fluxxrcor=',fluxxrcor
          if(fluxxrcor) then
            xrcor=1.
          else
            xrcor=0.
          endif
        elseif(keyword(1:nchark).eq.uppercase('fluxyrcor')) then
          read(27,*,iostat=inperr) fluxyrcor
          if(verbose.ge.1) write(*,*) 'fluxyrcor=',fluxyrcor
          if(fluxyrcor) then
            yrcor=1.
          else
            yrcor=0.
          endif
        elseif(keyword(1:nchark).eq.uppercase('iobuux')) then
          read(27,*,iostat=inperr) iobuux
          if(verbose.ge.1) write(*,*) 'iobuux=',iobuux
        elseif(keyword(1:nchark).eq.uppercase('iobuuy')) then
          read(27,*,iostat=inperr) iobuuy
          if(verbose.ge.1) write(*,*) 'iobuuy=',iobuuy
        elseif(keyword(1:nchark).eq.uppercase('iobvvx')) then
          read(27,*,iostat=inperr) iobvvx
          if(verbose.ge.1) write(*,*) 'iobvvx=',iobvvx
        elseif(keyword(1:nchark).eq.uppercase('iobvvy')) then
          read(27,*,iostat=inperr) iobvvy
          if(verbose.ge.1) write(*,*) 'iobvvy=',iobvvy
        elseif(keyword(1:nchark).eq.uppercase('iobptx')) then
          read(27,*,iostat=inperr) iobptx
          if(verbose.ge.1) write(*,*) 'iobptx=',iobptx
        elseif(keyword(1:nchark).eq.uppercase('iobpty')) then
          read(27,*,iostat=inperr) iobpty
          if(verbose.ge.1) write(*,*) 'iobpty=',iobpty
        elseif(keyword(1:nchark).eq.uppercase('iobppx')) then
          read(27,*,iostat=inperr) iobppx
          if(verbose.ge.1) write(*,*) 'iobppx=',iobppx
        elseif(keyword(1:nchark).eq.uppercase('iobppy')) then
          read(27,*,iostat=inperr) iobppy
          if(verbose.ge.1) write(*,*) 'iobppy=',iobppy
        elseif(keyword(1:nchark).eq.uppercase('perdam')) then
          read(27,*,iostat=inperr) perdam
          if(verbose.ge.1) write(*,*) 'perdam=',perdam
        elseif(keyword(1:nchark).eq.uppercase('dohsmo')) then
          read(27,*,iostat=inperr) dohsmo
          if(verbose.ge.1) write(*,*) 'dohsmo=',dohsmo
        elseif(keyword(1:nchark).eq.uppercase('dovsmo')) then
          read(27,*,iostat=inperr) dovsmo
          if(verbose.ge.1) write(*,*) 'dovsmo=',dovsmo
        elseif(keyword(1:nchark).eq.uppercase('dovsmt')) then
          read(27,*,iostat=inperr) dovsmt
          if(verbose.ge.1) write(*,*) 'dovsmt=',dovsmt
        elseif(keyword(1:nchark).eq.uppercase('numsmo')) then
          read(27,*,iostat=inperr) numsmo
          if(verbose.ge.1) write(*,*) 'numsmo=',numsmo
        elseif(keyword(1:nchark).eq.uppercase('hdampm')) then
          read(27,*,iostat=inperr) hdampm
          if(verbose.ge.1) write(*,*) 'hdampm=',hdampm
        elseif(keyword(1:nchark).eq.uppercase('hdampt')) then
          read(27,*,iostat=inperr) hdampt
          if(verbose.ge.1) write(*,*) 'hdampt=',hdampt
        elseif(keyword(1:nchark).eq.uppercase('vdamp')) then
          read(27,*,iostat=inperr) vdamp
          if(verbose.ge.1) write(*,*) 'vdamp=',vdamp
        elseif(keyword(1:nchark).eq.uppercase('vdampt')) then
          read(27,*,iostat=inperr) vdampt
          if(verbose.ge.1) write(*,*) 'vdampt=',vdampt
        elseif(keyword(1:nchark).eq.uppercase('iodif')) then
          read(27,*,iostat=inperr) iodif
          if(verbose.ge.0) write(*,*) 'iodif=',iodif
c     reading parameters for wind perturbation          
        elseif(keyword(1:nchark).eq.uppercase('uvpert')) then
          read(27,*,iostat=inperr) uvpert
          if(verbose.ge.0) write(*,*) 'uvpert=',uvpert
        elseif(keyword(1:nchark).eq.uppercase('f_chord')) then
          read(27,*,iostat=inperr) f_chord
          if(verbose.ge.0) write(*,*) 'f_chord=',f_chord
        elseif(keyword(1:nchark).eq.uppercase('tau_nudg')) then
          read(27,*,iostat=inperr) tau_nudg
          if(verbose.ge.0) write(*,*) 'tau_nudg=',tau_nudg
        elseif(keyword(1:nchark).eq.uppercase('nperturb')) then
          read(27,*,iostat=inperr) nperturb
          if(verbose.ge.0) write(*,*) 'nperturb=',nperturb
        elseif(keyword(1:nchark).eq.uppercase('pnudg')) then
          read(27,*,iostat=inperr) pnudg
          if(verbose.ge.0) write(*,*) 'pnudg=',pnudg
          tau_nudg=1./pnudg
        elseif(keyword(1:nchark).eq.uppercase('xm_pert')) then
          read(27,*,iostat=inperr) xm_pert
          if(verbose.ge.0) write(*,*) 'xm_pert=',xm_pert
        elseif(keyword(1:nchark).eq.uppercase('xice')) then
          read(27,*,iostat=inperr) xice
          if(verbose.ge.0) write(*,*) 'xice=',xice
        elseif(keyword(1:nchark).eq.uppercase('rikey')) then
          read(27,*,iostat=inperr) rikey
          if(verbose.ge.1) write(*,*) 'rikey=',rikey
        elseif(keyword(1:nchark).eq.uppercase('difcof')) then
          read(27,*,iostat=inperr) difcof
          if(verbose.ge.1) write(*,*) 'difcof=',difcof
        elseif(keyword(1:nchark).eq.uppercase('uvdif')) then
          read(27,*,iostat=inperr) uvdif
          if(verbose.ge.1) write(*,*) 'uvdif=',uvdif
        elseif(keyword(1:nchark).eq.uppercase('tsdif')) then
          read(27,*,iostat=inperr) tsdif
          if(verbose.ge.1) write(*,*) 'tsdif=',tsdif
        elseif(keyword(1:nchark).eq.uppercase('cdcoef')) then
          read(27,*,iostat=inperr) cdcoef
          if(verbose.ge.1) write(*,*) 'cdcoef=',cdcoef
        elseif(keyword(1:nchark).eq.uppercase('taum')) then
          read(27,*,iostat=inperr) taum
          if(verbose.ge.1) write(*,*) 'taum=',taum
        elseif(keyword(1:nchark).eq.uppercase('taut')) then
          read(27,*,iostat=inperr) taut
          if(verbose.ge.1) write(*,*) 'taut=',taut
        elseif(keyword(1:nchark).eq.uppercase('ntaut')) then
          read(27,*,iostat=inperr) ntaut
          if(verbose.ge.1) write(*,*) 'ntaut=',ntaut
          taut=ntaut*dt
        elseif(keyword(1:nchark).eq.uppercase('ntauspon')) then
          read(27,*,iostat=inperr) ntauspon
          if(verbose.ge.1) write(*,*) 'ntauspon=',ntauspon
          tauspon=ntauspon*dt
        elseif(keyword(1:nchark).eq.uppercase('zm')) then
          read(27,*,iostat=inperr) zm
          if(verbose.ge.1) write(*,*) 'zm=',zm
        elseif(keyword(1:nchark).eq.uppercase('zt')) then
          read(27,*,iostat=inperr) zt
          if(verbose.ge.1) write(*,*) 'zt=',zt
        elseif(keyword(1:nchark).eq.uppercase('nxsponge1')) then
          read(27,*,iostat=inperr) nxsponge1
          if(verbose.ge.1) write(*,*) 'nxsponge1=',nxsponge1
        elseif(keyword(1:nchark).eq.uppercase('nxsponge2')) then
          read(27,*,iostat=inperr) nxsponge2
          if(verbose.ge.1) write(*,*) 'nxsponge2=',nxsponge2
        elseif(keyword(1:nchark).eq.uppercase('nysponge1')) then
          read(27,*,iostat=inperr) nysponge1
          if(verbose.ge.1) write(*,*) 'nysponge1=',nysponge1
        elseif(keyword(1:nchark).eq.uppercase('nysponge2')) then
          read(27,*,iostat=inperr) nysponge2
          if(verbose.ge.1) write(*,*) 'nysponge2=',nysponge2
        elseif(keyword(1:nchark).eq.uppercase('text')) then
          itext=itext+1
          if(itext.le.3) then
            read(27,*,iostat=inperr) comlin(itext)
            if(verbose.ge.1) write(*,*) 'comlin=',comlin(itext)
          endif
        elseif(keyword(1:nchark).eq.uppercase('ifluxcor')) then
          read(27,*,iostat=inperr) ifluxcor
          if(verbose.ge.1) write(*,*) 'ifluxcor=',ifluxcor
          fluxsel=.true.
        elseif(keyword(1:nchark).eq.uppercase('restar')) then
          read(27,*,iostat=inperr) restar
          if(verbose.ge.1) write(*,*) 'restar=',restar
        elseif(keyword(1:nchark).eq.uppercase('rstout')) then
          read(27,*,iostat=inperr) rstout
          if(verbose.ge.1) write(*,*) 'rstout=',rstout
        elseif(keyword(1:nchark).eq.uppercase('iotrst')) then
          read(27,*,iostat=inperr) iotrst
          if(verbose.ge.1) write(*,*) 'iotrst=',iotrst
        elseif(keyword(1:nchark).eq.uppercase('forout')) then
          read(27,*,iostat=inperr) forout
          if(verbose.ge.1) write(*,*) 'forout=',forout
        elseif(keyword(1:nchark).eq.uppercase('iotfor')) then
          read(27,*,iostat=inperr) iotfor
          if(verbose.ge.1) write(*,*) 'iotfor=',iotfor
        elseif(keyword(1:nchark).eq.uppercase('unfout')) then
          read(27,*,iostat=inperr) unfout
          if(verbose.ge.1) write(*,*) 'unfout=',unfout
        elseif(keyword(1:nchark).eq.uppercase('iotunf')) then
          read(27,*,iostat=inperr) iotunf
          if(verbose.ge.1) write(*,*) 'iotunf=',iotunf
        elseif(keyword(1:nchark).eq.uppercase('iostart')) then
          read(27,*,iostat=inperr) iostart
          if(verbose.ge.1) write(*,*) 'iotstart=',iostart
        elseif(keyword(1:nchark).eq.uppercase('ioend')) then
          read(27,*,iostat=inperr) ioend
          if(verbose.ge.1) write(*,*) 'ioend=',ioend
        elseif(keyword(1:nchark).eq.uppercase('iodelta')) then
          read(27,*,iostat=inperr) iodelta
          if(verbose.ge.1) write(*,*) 'iodelta=',iodelta
        elseif(keyword(1:nchark).eq.uppercase('drgprt')) then
          read(27,*,iostat=inperr) drgprt
          if(verbose.ge.1) write(*,*) 'drgprt=',drgprt
        elseif(keyword(1:nchark).eq.uppercase('domom')) then
          read(27,*,iostat=inperr) domom
          if(verbose.ge.1) write(*,*) 'domom=',domom
        elseif(keyword(1:nchark).eq.uppercase('iotdrg')) then
          read(27,*,iostat=inperr) iotdrg
          if(verbose.ge.1) write(*,*) 'iotdrg=',iotdrg
        elseif(keyword(1:nchark).eq.uppercase('iotduv')) then
          read(27,*,iostat=inperr) iotduv
          if(verbose.ge.1) write(*,*) 'iotduv=',iotduv
        elseif(keyword(1:nchark).eq.uppercase('nlev')) then
          read(27,*,iostat=inperr) nlev
          if(verbose.ge.0) write(*,*) 'nlev=',nlev
        elseif(keyword(1:nchark).eq.uppercase('dzlev')) then
          read(27,*,iostat=inperr) dzlev
          if(verbose.ge.0) write(*,*) 'dzlev=',dzlev
        elseif(keyword(1:nchark).eq.uppercase('mxspan')) then
          read(27,*,iostat=inperr) mxspan
          if(verbose.ge.1) write(*,*) 'mxspan=',mxspan
        elseif(keyword(1:nchark).eq.uppercase('iomod')) then
          read(27,*,iostat=inperr) iomod
          if(verbose.ge.1) write(*,*) 'iomod=',iomod
        elseif(keyword(1:nchark).eq.uppercase('proout')) then
          read(27,*,iostat=inperr) proout
          if(verbose.ge.1) write(*,*) 'proout=',proout
        elseif(keyword(1:nchark).eq.uppercase('norout')) then
          read(27,*,iostat=inperr) norout
          if(verbose.ge.1) write(*,*) 'norout=',norout
        elseif(keyword(1:nchark).eq.uppercase('surout')) then
          read(27,*,iostat=inperr) surout
          if(verbose.ge.1) write(*,*) 'surout=',surout
        elseif(keyword(1:nchark).eq.uppercase('masout')) then
          read(27,*,iostat=inperr) masout
          if(verbose.ge.1) write(*,*) 'masout=',masout
        elseif(keyword(1:nchark).eq.uppercase('iniout')) then
          read(27,*,iostat=inperr) iniout
          if(verbose.ge.1) write(*,*) 'iniout=',iniout
        elseif(keyword(1:nchark).eq.uppercase('iniunf')) then
          read(27,*,iostat=inperr) iniunf
          if(verbose.ge.1) write(*,*) 'iniunf=',iniunf
        elseif(keyword(1:nchark).eq.uppercase('echeck')) then
          read(27,*,iostat=inperr) echeck
          if(verbose.ge.1) write(*,*) 'echeck=',echeck
        elseif(keyword(1:nchark).eq.uppercase('eneout')) then
          read(27,*,iostat=inperr) eneout
          if(verbose.ge.1) write(*,*) 'eneout=',eneout
        elseif(keyword(1:nchark).eq.uppercase('morout')) then
          read(27,*,iostat=inperr) morout
          if(verbose.ge.1) write(*,*) 'morout=',morout
        elseif(keyword(1:nchark).eq.uppercase('refout')) then
          read(27,*,iostat=inperr) refout
          if(verbose.ge.1) write(*,*) 'refout=',refout
        elseif(keyword(1:nchark).eq.uppercase('inrout')) then
          read(27,*,iostat=inperr) inrout
          if(verbose.ge.1) write(*,*) 'inrout=',inrout
        elseif(keyword(1:nchark).eq.uppercase('fluout')) then
          read(27,*,iostat=inperr) fluout
          if(verbose.ge.1) write(*,*) 'fluout=',fluout
        elseif(keyword(1:nchark).eq.uppercase('ppbout')) then
          read(27,*,iostat=inperr) ppbout
          if(verbose.ge.1) write(*,*) 'ppbout=',ppbout
        elseif(keyword(1:nchark).eq.uppercase('duvout')) then
          read(27,*,iostat=inperr) duvout
          if(verbose.ge.1) write(*,*) 'duvout=',duvout
        elseif(keyword(1:nchark).eq.uppercase('uubout')) then
          read(27,*,iostat=inperr) uubout
          if(verbose.ge.1) write(*,*) 'uubout=',uubout
        elseif(keyword(1:nchark).eq.uppercase('ptbout')) then
          read(27,*,iostat=inperr) ptbout
          if(verbose.ge.1) write(*,*) 'ptbout=',ptbout
        elseif(keyword(1:nchark).eq.uppercase('difout')) then
          read(27,*,iostat=inperr) difout
          if(verbose.ge.1) write(*,*) 'difout=',difout
        elseif(keyword(1:nchark).eq.uppercase('driout')) then
          read(27,*,iostat=inperr) driout
          if(verbose.ge.1) write(*,*) 'driout=',driout
        elseif(keyword(1:nchark).eq.uppercase('ddiout')) then
          read(27,*,iostat=inperr) ddiout
          if(verbose.ge.1) write(*,*) 'ddiout=',ddiout
        elseif(keyword(1:nchark).eq.uppercase('defout')) then
          read(27,*,iostat=inperr) defout
          if(verbose.ge.1) write(*,*) 'defout=',defout
        elseif(keyword(1:nchark).eq.uppercase('eldout')) then
          read(27,*,iostat=inperr) eldout
          if(verbose.ge.1) write(*,*) 'eldout=',eldout
        elseif(keyword(1:nchark).eq.uppercase('elfout')) then
          read(27,*,iostat=inperr) elfout
          if(verbose.ge.1) write(*,*) 'elfout=',elfout
        elseif(keyword(1:nchark).eq.uppercase('elgout')) then
          read(27,*,iostat=inperr) elgout
          if(verbose.ge.1) write(*,*) 'elgout=',elgout
        elseif(keyword(1:nchark).eq.uppercase('phsout')) then
          read(27,*,iostat=inperr) phsout
          if(verbose.ge.1) write(*,*) 'phsout=',phsout
        elseif(keyword(1:nchark).eq.uppercase('desout')) then
          read(27,*,iostat=inperr) desout
          if(verbose.ge.1) write(*,*) 'desout=',desout
        elseif(keyword(1:nchark).eq.uppercase('ixss')) then
          read(27,*,iostat=inperr) ixss
          if(verbose.ge.1) write(*,*) 'ixss=',ixss
        elseif(keyword(1:nchark).eq.uppercase('iyss')) then
          read(27,*,iostat=inperr) iyss
          if(verbose.ge.1) write(*,*) 'iyss=',iyss
        elseif(keyword(1:nchark).eq.uppercase('isss')) then
          read(27,*,iostat=inperr) isss
          if(verbose.ge.1) write(*,*) 'isss=',isss
        elseif(keyword(1:nchark).eq.uppercase('ixp')) then
          read(27,*,iostat=ixperr) (ixp(ii),ii=1,nplan)
          if(verbose.ge.1) write(*,'(a5,10i3)') 'ixp=',ixp
        elseif(keyword(1:nchark).eq.uppercase('iyp')) then
          read(27,*,iostat=ixperr) (iyp(ii),ii=1,nplan)
          if(verbose.ge.1) write(*,'(a5,10i3)') 'iyp=',iyp
        elseif(keyword(1:nchark).eq.uppercase('isp')) then
          read(27,*,iostat=ixperr) (isp(ii),ii=1,nplan)
          if(verbose.ge.1) write(*,'(a5,10i3)') 'isp=',isp
        elseif(keyword(1:nchark).eq.uppercase('zplan')) then
          read(27,*,iostat=ixperr) (zplan(ii),ii=1,nplan)
          if(verbose.ge.1) write(*,'(a5,10f7.0)') 'zplan=',zplan
        elseif(keyword(1:nchark).eq.uppercase('scalew')) then
          read(27,*,iostat=inperr) scalew
          if(verbose.ge.1) write(*,'(a5,10f7.0)') 'scalew=',scalew
        elseif(keyword(1:nchark).eq.uppercase('zconst')) then
          read(27,*,iostat=inperr) zconst
          if(verbose.ge.1) write(*,*) 'zconst=',zconst
        elseif(keyword(1:nchark).eq.uppercase('ifhle')) then
          read(27,*,iostat=inperr) ifhle
          if(verbose.ge.1) write(*,*) 'ifhle=',ifhle
	  elseif(keyword(1:nchark).eq.uppercase('tfix')) then                     DC 06.2010
          read(27,*,iostat=inperr) tfix                                         DC 06.2010
          if(verbose.ge.1) write(*,*) 'tfix=',tfix                              DC 06.2010
        elseif(keyword(1:nchark).eq.uppercase('profile')) then
          read(29,*,iostat=inperr) pts0,ioreft,iorefq,iorefu,iorefv
          read(29,*,iostat=inperr) ndth,ndus,ndvs,ndqvs
          ndat=ndth
          nprof=1
          call allocref
          do ii=0,ndat
            read (29,*,iostat=inperr) zthdat(ii,1),thdat(ii,1)
     :        ,zusdat(ii,1),usdat(ii,1)
     :        ,zvsdat(ii,1),vsdat(ii,1),zqvsdat(ii,1),qvsdat(ii,1)
!     :         ,zthdat(ii,1),per_dat(ii,1)
          enddo
          pressdat(0,1)=pa
          dzpro=zthdat(ndat,1)/npro
          if(verbose.ge.1) write(*,*) 'Profile given'
          ok_profile=.true.
          starttime = timeini                                                   VS,05.2007

c
c new option nhad0518 'fnprofile' the profile is given in the 'ref' file
c
        elseif(keyword(1:nchark).eq.uppercase('fnprofile')) then
          read(27,*,iostat=inperr) fnfor
          open(35,file=fnfor,status='old')
          read(35,'(a8)') line
          read(35,*) ndatp1,nprof,realtime
     :      ,ioreft,iorefq,iorefu,iorefv,pts0,pa
          ok_fnfor=.true.
          ndat=ndatp1
          ndth=ndat
          ndus=ndat
          ndvs=ndat
          ndqvs=ndat
          call allocref
          read(35,'(a8)') line
          read(35,*) tprof1
          starttime=tprof1
          do ip=1,nprof
            if(verbose.ge.2) write(*,*) 'Profile:',ip,tprof1
            read(35,'(a8)') line
            read(35,*) xprof(ip),yprof(ip)
            read(35,'(a8)') line
            do ii=0,ndat-1
              read (35,*) zthdat(ii,ip),thdat(ii,ip)
     :          ,usdat(ii,ip),vsdat(ii,ip),qvsdat(ii,ip),pressdat(ii,ip)
              if(verbose.ge.3)
     :         write (*,'(6e15.7)') zthdat(ii,ip),thdat(ii,ip)
     :          ,usdat(ii,ip),vsdat(ii,ip),qvsdat(ii,ip),pressdat(ii,ip)
            enddo
c
c extrapolate profile to 20km
c
            zthdat(ndat,ip)=max(20000.,zthdat(ndat-1,ip))
            usdat(ndat,ip)=usdat(ndat-1,ip)
            vsdat(ndat,ip)=vsdat(ndat-1,ip)
            qvsdat(ndat,ip)=qvsdat(ndat-1,ip)
            thdat(ndat,ip)=thdat(ndat-1,ip)
            pressdat(ndat,ip)=100.  ! not relevant?
            pa=pressdat(0,ip)
          enddo
c
c nhad0523: if one needs different z's for the different variables, they must be
c   read independently
c
          zusdat=zthdat
          zvsdat=zthdat
          zqvsdat=zthdat

c         rewind(35)
          if(verbose.ge.1) write(*,*) 'Profiles given',nprof
          ok_profile=.true.
        elseif(keyword(1:nchark).eq.uppercase('forcedpp')) then
          read(27,*,iostat=inperr) forcedpp
          if(verbose.ge.1) write(*,*) 'forcedpp=',forcedpp
        elseif(keyword(1:nchark).eq.uppercase('forcepsuf')) then
          read(27,*,iostat=inperr) forcepsuf
          if(verbose.ge.1) write(*,*) 'forcepsuf=',forcepsuf
        elseif(keyword(1:nchark).eq.uppercase('npro')) then
          read(27,*,iostat=inperr) npro
          if(verbose.ge.1) write(*,*) 'npro=',npro
        elseif(keyword(1:nchark).eq.uppercase('dzpro')) then
          read(27,*,iostat=inperr) dzpro
          if(verbose.ge.1) write(*,*) 'dzpro=',dzpro
        elseif(keyword(1:nchark).eq.uppercase('fnrad')) then
          read(27,*,iostat=inperr) fnrad
          if(verbose.ge.1) write(*,*) 'fnrad=',fnrad
          open(40,file=fnrad,status='old')
c
c xrg - global solar radiation (W m**-2)
c xra - atmospheric radiation (downward)
c xpl - precipitation rate (kg/m**2/s)
c nf - time in seconds between radiation data
c
          read(40,'(a)') line
          read(40,*,iostat=inperr) ndados,nf
          allocate(xrg(0:ndados))
          allocate(xra(0:ndados))
          allocate(xpl(0:ndados))
          read(40,'(a)') line
          do ii=0,ndados
            read (40,*,iostat=inperr) xra(ii),xrg(ii),xpl(ii)                   VS,06.2007
          enddo
          close(40)
          if(verbose.ge.1) write(*,*) 'Radiation given'
          ok_rad=.true.
        elseif(keyword(1:nchark).eq.uppercase('ug0')) then
          read(27,*,iostat=inperr) ug0
          if(verbose.ge.1) write(*,*) 'ug0=',ug0
        elseif(keyword(1:nchark).eq.uppercase('vg0')) then
          read(27,*,iostat=inperr) vg0
          if(verbose.ge.1) write(*,*) 'vg0=',vg0
        elseif(keyword(1:nchark).eq.uppercase('itepsuf')) then
          read(27,*,iostat=inperr) itepsuf
          if(verbose.ge.1) write(*,*) 'itepsuf=',itepsuf
        elseif(keyword(1:nchark).eq.uppercase('sigmagrid')) then
          if(.not.ok_ns) then
            write(0,*) 'nh3dFatalError='
     :      ,'NS must be given before sigmagrid'
            stop
          endif
          autogrid=.false.
          if(verbose.ge.1) write(*,*) 'User defined sigma grid'
          allocate(sigma0(0:ns+1))
          allocate(sigma1(0:ns+1))
          read (29,*) (sigma1(is),is=1,ns)
          do is=1,ns
          write(0,*) is,sigma1(is)
          enddo
          if(sigma1(1).ne.0. .or. sigma1(ns).ne.1.) then
            write(0,*) 'nh3dFatalError=Inconsistent sigma grid'
            write(*,*) sigma1(1),sigma1(ns)
            stop
          endif
          sigma1(0)=-sigma1(2)
          sigma1(ns+1)=1.+(sigma1(ns)-sigma1(ns-1))
          do is=1,ns+1
            sigma0(is-1)=(sigma1(is)+sigma1(is-1))/2.
            if(verbose.ge.1) write(*,*) sigma1(is-1)
          enddo
          sigma0(ns+1)=2.*sigma0(ns)-sigma0(ns-1)
        elseif(keyword(1:nchark).eq.uppercase('fnmap')) then
          read(27,*,iostat=inperr) fnmap
          if(verbose.ge.0) write(*,*) 'fnmap=',fnmap
          ok_fnmap=.true.
        elseif(keyword(1:nchark).eq.uppercase('fnout')) then
          read(27,*,iostat=inperr) fnout
          if(verbose.ge.0) write(*,*) 'fnout=',fnout
          ok_fnout=.true.
        elseif(keyword(1:nchark).eq.uppercase('fngrd')) then
          read(27,*,iostat=inperr) fngrd
          if(verbose.ge.0) write(*,*) 'fngrd=',fngrd
          ok_fngrd=.true.
        elseif(keyword(1:nchark).eq.uppercase('fnrst')) then
          read(27,*,iostat=inperr) fnrst
          if(verbose.ge.0) write(*,*) 'fnrst=',fnrst
          ok_fnrst=.true.
        elseif(keyword(1:nchark).eq.uppercase('presout')) then
          read(29,*,iostat=inperr) noutu
          if(allocated(ixoutu)) deallocate(ixoutu)
          allocate(ixoutu(noutu))
          if(allocated(iyoutu)) deallocate(iyoutu)
          allocate(iyoutu(noutu))
          if(allocated(isoutu)) deallocate(isoutu)
          allocate(isoutu(noutu))
          read(29,*) (ixoutu(k),iyoutu(k),isoutu(k),k=1,noutu)
          if(verbose.ge.1) then
            write(*,*) 'Local Pressure Output ON'
            write(*,'(3i3)') (ixoutu(k),iyoutu(k),isoutu(k),k=1,noutu)
          endif
        elseif(keyword(1:nchark).eq.uppercase('profout')) then
          read(29,*,iostat=inperr) noutp
          if(allocated(ixoutp)) deallocate(ixoutp)
          allocate(ixoutp(noutp))
          if(allocated(iyoutp)) deallocate(iyoutp)
          allocate(iyoutp(noutp))
          do k=1,noutp
            read(29,*) ixoutp(k),iyoutp(k)
          enddo
          if(verbose.ge.1) then
            write(*,*) 'Profile Output ON'
            write(*,'(2i3)') (ixoutp(k),iyoutp(k),k=1,noutp)
          endif
        elseif(keyword(1:nchark).eq.uppercase('ifsoil')) then
          read(27,*,iostat=inperr) ifsoil
          if(verbose.ge.1) write(*,*) 'ifsoil=',ifsoil
          if(verbose.ge.0 .and. ifsoil.ne.0) write(*,*) 'Soil ON'
        elseif(keyword(1:nchark).eq.uppercase('z0hz0')) then
          read(27,*,iostat=inperr) zohz0
          if(verbose.ge.1) write(*,*) 'zohz0=',zohz0
        elseif(keyword(1:nchark).eq.uppercase('qif')) then
          read(27,*,iostat=inperr) qif
          if(verbose.ge.1) write(*,*) 'qif=',qif
          isqv=nint(qif)
          write(*,*) 'isqv=>',isqv
        elseif(keyword(1:nchark).eq.uppercase('qdrag')) then
          read(27,*,iostat=inperr) qdrag
          if(verbose.ge.1) write(*,*) 'qdrag=',qdrag
        elseif(keyword(1:nchark).eq.uppercase('ifqc')) then
          read(27,*,iostat=inperr) ifqc
          if(verbose.ge.1) write(*,*) 'ifqc=',ifqc
          isqc=ifqc
          write(*,*) 'isqc=>',isqc
        elseif(keyword(1:nchark).eq.uppercase('ifqr')) then
          read(27,*,iostat=inperr) ifqr
          if(verbose.ge.1) write(*,*) 'ifqr=',ifqr
          isqr=ifqr
          write(*,*) 'isqr=>',isqr
	 elseif(keyword(1:nchark).eq.uppercase('ifqi')) then                        DC,10.2009
          read(27,*,iostat=inperr) ifqi                                           DC.10,2009
	    if(verbose.ge.1) write(*,*) 'ifqi=',ifqi                                DC,10.2009
          isqi=ifqi                                                               DC,10.2009
          write(*,*) 'isqi=>',isqi                                                DC,10.2009
        elseif(keyword(1:nchark).eq.uppercase('ifblob')) then
          read(27,*,iostat=inperr) ifblob
          if(verbose.ge.1) write(*,*) 'ifblob=',ifblob
        elseif(keyword(1:nchark).eq.uppercase('ifperturb')) then
          read(27,*,iostat=inperr) ifperturb
          if(verbose.ge.1) write(*,*) 'ifperturb=',ifperturb
        elseif(keyword(1:nchark).eq.uppercase('ifdirichlet')) then
          read(27,*,iostat=inperr) ifdirichlet
          if(verbose.ge.1) write(*,*) 'ifdirichlet=',ifdirichlet
        elseif(keyword(1:nchark).eq.uppercase('iophis')) then
          read(27,*,iostat=inperr) iophis
          if(verbose.ge.1) write(*,*) 'iophis=',iophis
        elseif(keyword(1:nchark).eq.uppercase('ifuprefs')) then
          read(27,*,iostat=inperr) ifuprefs
          if(verbose.ge.1) write(*,*) 'ifuprefs=',ifuprefs
        elseif(keyword(1:nchark).eq.uppercase('divrefstate')) then
          read(27,*,iostat=inperr) divrefstate
          if(verbose.ge.1) write(*,*) 'divrefstate=',divrefstate
        elseif(keyword(1:nchark).eq.uppercase('iwind')) then
          read(27,*,iostat=inperr) iwind
          if(verbose.ge.1) write(*,*) 'iwind=',iwind
        elseif(keyword(1:nchark).eq.uppercase('itheta')) then
          read(27,*,iostat=inperr) itheta
          if(verbose.ge.1) write(*,*) 'itheta=',itheta
        elseif(keyword(1:nchark).eq.uppercase('filt')) then
          read(27,*,iostat=inperr) filt
          if(verbose.ge.1) write(*,*) 'filt=',filt
        elseif(keyword(1:nchark).eq.uppercase('ppfilt')) then
          read(27,*,iostat=inperr) ppfilt
          if(verbose.ge.1) write(*,*) 'ppfilt=',ppfilt
        elseif(keyword(1:nchark).eq.uppercase('nbyte_acc')) then
          read(27,*,iostat=inperr) nbyte_acc
          if(verbose.ge.1) write(*,*) 'nbyte_acc=',nbyte_acc
        elseif(keyword(1:nchark).eq.uppercase('grdout')) then
          read(27,*,iostat=inperr) grdout
          if(verbose.ge.1) write(*,*) 'grdout=',grdout

! new way of getting 1 line input for selection of profiles and points
! also new way for topography (mountain) data
! also new simple profile

        elseif(keyword(1:nchark).eq.uppercase('outprof')) then
          new_prof=.true.
          noutp=noutp+1
        elseif(keyword(1:nchark).eq.uppercase('outpoint')) then
          new_point=.true.
          noutu=noutu+1
        elseif(keyword(1:nchark).eq.uppercase('topography')) then
          read(27,*,iostat=inperr) iomtyp
          if(iomtyp.ne.4) then
            anamount=.true.
            backspace(27)
            read(27,*,iostat=inperr) iomtyp,hmount,xmount,ymount
     :       ,hxwid,hywid,hexp
          endif
          if(verbose.ge.1) then
            write(*,'(i2,6f10.2)') iomtyp,hmount,xmount,ymount
     :       ,hxwid,hywid,hexp
          endif
        elseif(keyword(1:nchark).eq.uppercase('simpleatmos')) then
          if(ok_profile) write(*,*) 'Warning: Profile redefined!'
          read(27,*) pts0,bv0,u0,v0,rh0,pa
          if(rh0.lt.0. .or. rh0.gt.1.) then
            write(0,*) 'nh3dFatalError='
     :        ,'Invalid relative humidity in profile'
            stop
          endif
          ifuprefs=0
          ndatp1=31
          nprof=1
          realtime=.false.
          ioreft=3
          iorefq=1
          iorefu=0
          iorefv=0
          ok_fnfor=.false.
          ndat=ndatp1-1
          ndth=ndat
          ndus=ndat
          ndvs=ndat
          ndqvs=ndat
          call allocref
          starttime=timeini
          ip=1
          do ii=0,ndat
            zthdat(ii,ip)=ii*1000.
            thdat(ii,ip)=bv0
            pressdat(ii,ip)=pa
            usdat(ii,ip)=u0
            vsdat(ii,ip)=v0
	      qvsdat(ii,ip)=rh0                                                        DC,11.2009
            !if(ii.gt.2.and.ii.lt.5) qvsdat(ii,ip)=0.98                              DC,11.2009
	      !if(ii.gt.4.and.ii.lt.7) qvsdat(ii,ip)=0.98                              DC,11.2009
	      !if(ii.gt.6.and.ii.lt.9) qvsdat(ii,ip)=0.9                              DC,11.2009
	      !if(ii.gt.8) qvsdat(ii,ip)=0.2                                        DC,11.2009
          enddo
          pa=pressdat(0,ip)
          zusdat=zthdat
          zvsdat=zthdat
          zqvsdat=zthdat
          if(verbose.ge.1) write(*,*) 'Simple Profile given'
          ok_profile=.true.

        else
          write(0,*) 'nh3dFatalError=Unrecognized keyword !!:',keyword
          stop
        endif
        if(inperr.gt.0) then
          write(0,*) 'nh3dFatalError=Input error from nh3d.dat',keyword
          stop
        endif
      enddo

      if(new_prof .or. new_point) then
        if(new_prof) then
          if(allocated(ixoutp)) deallocate(ixoutp)
          allocate(ixoutp(noutp))
          if(allocated(iyoutp)) deallocate(iyoutp)
          allocate(iyoutp(noutp))
          ioutp=0
        endif
        if(new_point) then
          if(allocated(ixoutu)) deallocate(ixoutu)
          allocate(ixoutu(noutu))
          if(allocated(iyoutu)) deallocate(iyoutu)
          allocate(iyoutu(noutu))
          if(allocated(isoutu)) deallocate(isoutu)
          allocate(isoutu(noutu))
          ioutu=0
          if(noutu.gt.20) then
            write(0,*) 'nh3dFatalError=','More than 20 outpoints'
            write(0,*) 'If you need more outpoints the code must be'
            write(0,*) 'changed look for "write(52," '
            stop
          endif
        endif
        rewind(29)
        inperr=0
        do while(inperr.eq.0)
          read(29,'(a)') line
          ichlef=1
          do i=1,80
            if(line(i:i).ne.' ') exit
            ichlef=i
          enddo
          if(ichlef.eq.80) cycle
          ichrig=ichlef
          do i=ichlef+1,80
            if(ichar(line(i:i)).lt.icha.or.ichar(line(i:i)).gt.ichz)exit
            ichrig=i
          enddo
          keyword=line(ichlef:ichrig)
          nchark=ichrig-ichlef+1
          call toupcase(keyword)
          if(verbose.ge.2) write(*,*) keyword,nchark,ichlef,ichrig
          write(27,'(a)') line(ichrig+1:80)
c         write(*,'(a)') line(ichrig+1:80)
          backspace(27)
          if(keyword(1:3).eq.uppercase('end')) then
            exit
          elseif(keyword(1:1).eq.'#') then
            cycle
          elseif(keyword(1:nchark).eq.uppercase('outprof')) then
            ioutp=ioutp+1
            read(27,*) ixoutp(ioutp),iyoutp(ioutp)
            if(ixoutp(ioutp).lt.1 .or. ixoutp(ioutp).gt.nx .or.
     :        iyoutp(ioutp).lt.1. .or. iyoutp(ioutp).gt.ny) then
              write(0,*) 'nh3dFatalError=Profile out of the grid:'
     :          ,ixoutp(ioutp),iyoutp(ioutp)
              stop
            endif
          elseif(keyword(1:nchark).eq.uppercase('outpoint')) then
            ioutu=ioutu+1
            read(27,*) ixoutu(ioutu),iyoutu(ioutu),isoutu(ioutu)
            if(ixoutu(ioutu).lt.1 .or. ixoutu(ioutu).gt.nx .or.
     :        iyoutu(ioutu).lt.1. .or. iyoutu(ioutu).gt.ny .or.
     :        isoutu(ioutu).lt.1. .or. isoutu(ioutu).gt.ns) then
              write(0,*) 'nh3dFatalError=Point out of the grid:'
     :          ,ixoutu(ioutu),iyoutu(ioutu),isoutu(ioutu)
              stop
            endif
          endif
          if(inperr.gt.0) then
            write(0,*) 'nh3dFatalError=Input error from nh3d.dat'
     :        ,keyword
            stop
          endif
        enddo
        if(new_point .and. verbose.ge.1) then
          write(*,*) 'Local Output ON'
          write(*,'(3i3)') (ixoutu(k),iyoutu(k),isoutu(k),k=1,noutu)
        endif
        if(new_prof .and. verbose.ge.1) then
          write(*,*) 'Profile Output ON'
          write(*,'(2i3)') (ixoutp(k),iyoutp(k),k=1,noutp)
        endif
      endif

      if(verbose.ge.0.and.qif.ne.0.) write(*,*) 'Water vapour ON'

      if(verbose.ge.0.and.ifqc.ne.0) write(*,*) 'Cloud Water ON'

      if(verbose.ge.0.and.ifqr.ne.0) write(*,*) 'Rain Water ON'

	if(verbose.ge.0.and.ifqi.ne.0) write(*,*) 'Cloud Ice and Snow ON'                   DC,10.2009

c
c check data
c
      if(.not.ok_nx) then
        write(0,*) 'nh3dFatalError=Aborted due to missing data: NX'
        stop
      endif
      if(.not.ok_ny) then
        write(0,*) 'nh3dFatalError=Aborted due to missing data: NY'
        stop
      endif
      if(.not.ok_ns) then
        write(0,*) 'nh3dFatalError=Aborted due to missing data: NS'
        stop
      endif
      if(.not.ok_dx) then
        write(0,*) 'nh3dFatalError=Aborted due to missing data: DX'
        stop
      endif
      if(.not.ok_dy) then
        write(0,*) 'nh3dFatalError=Aborted due to missing data: DY'
        stop
      endif
      if(.not.ok_ntime) then
        write(0,*) 'nh3dFatalError=Aborted due to missing data: NTIME'
        stop
      endif
      if(.not.ok_profile) then
        write(0,*) 'nh3dFatalError=Aborted due to missing data: PROFILE'
        stop
      endif
      if(.not.ok_fnmap .and. ((ifsoil.ne.0).or.(iomtyp.eq.4)))then
        write(0,*) 'nh3dFatalError=Aborted due to missing data: FNMAP'
        stop
      endif
      if(.not.ok_fnout) then
        if(ok_fngrd) then
          fnout=fngrd
        else
          write(0,*) 'nh3dFatalError=missing data: FNOUT or fngrd'
          stop
        endif
      endif
      if(.not.ok_fngrd) then
        if(ok_fnout) then
          fngrd=fnout
        else
          write(0,*) 'nh3dFatalError=missing data: FNGRD or fnout'
          stop
        endif
      endif
c     if(.not.ok_fnrst) then
c        write(*,*) 'Aborted due to missing data: FNRST'
c        stop
c      endif

      exptit=fnout(1:8)
      if(.not.ok_rad  .and. ifsoil.eq.1) then
        write(0,*) 'nh3dFatalError=missing radiation data: RADIATION'
        write(0,*) 'And Soil Model ON'
        stop
      endif
	if (radpar==1) then                                                       VS,29.02.2008
        if(ifsoil.eq.1 .and. ndados*nf.lt.ntime*dt) then
          write(0,*) 'nh3dFatalError=not enough radiation data'
          write(0,*) 'And Soil Model ON'
          stop
        endif
	endif                                                                     VS,29.02.2008

c not needed anymore

      if(qif.ne.0. .and. isqv.eq.0) then
        write(0,*) 'nh3dFatalError=Aborted due to isqv.ne.qif'
        write(0,*) 'No allocation for water vapor (qv)'
        stop
      endif

c not needed anymore

      if(isqc.eq.0 .and. ifqc.ne.0) then
        write(0,*) 'nh3dFatalError=Aborted due to isqc.ne.ifqc'
        write(0,*) 'No allocation for cloud water (qc)'
        stop
      endif

c not needed anymore

      if(isqr.eq.0 .and. ifqr.ne.0) then
        write(0,*) 'nh3dFatalError=Aborted due to isqr.ne.ifqr'
        write(0,*) 'No allocation for rain water (qc)'
        stop
      endif

c not needed anymore

      if(isqc.ne.ifqc .or. isqr.ne.ifqr .or. nint(qif).ne.isqv) then
        write(0,*) 'nh3dWarning= You may save memory requirements if'
     :    ,'you turn off some water variables and recompile'
      endif
c
c coriolis parameter:
c
      if(fcor.ne.-999. .and. xlatit.ne.0.) then
        write(0,*) 'nh3dFatalError= conflicting options for fcor'
        stop
      endif
      if(fcor.eq.-999.) then
        fcor=2.*omega*sin(xlatit*pi/180.)
        if(verbose.ge.2) write(*,*) 'fcor=',fcor,omega,xlatit,pi
      endif

      dif3d=.false.
      diftx=.false.
	difloc=.false.
	nonloc_tm86=.false.
	nonloc_ls96=.false.
        nonloc_Noh=.false.
      if(iodif.ge.20) then
        diftx=.true.
        iodif=iodif-20
      elseif(iodif.ge.10) then
        iodif=iodif-10
	   if(iodif.eq.3) then
	   difloc=.true.
	   endif
	   if(iodif.eq.4) then
	   nonloc_tm86=.true.
	   endif
	   if(iodif.eq.5) then
	   nonloc_ls96=.true.
	   endif
           if(iodif.eq.6) then
	   nonloc_Noh=.true.
	   endif
      else
        dif3d=.true.
      endif
      if(iodif.lt.0 .or. iodif.gt.6) then
        write(0,*) 'nh3dFatalError=Aborted due wrong iodif',iodif
        stop
      endif

      ok_plan=.true.
      do ii=1,nplan
        if(ixp(ii).gt.nx) then
          ok_plan=.false.
          exit
        endif
        if(iyp(ii).gt.ny) then
          ok_plan=.false.
          exit
        endif
        if(isp(ii).gt.ns) then
          ok_plan=.false.
          exit
        endif
      enddo

      if(.not.ok_plan) then
        write(0,*) 'nh3dFatalError= Data out of the grid: IXP,IYP.ISP'
        stop
      endif

      if(nlev.le.0) nlev=ns

      if(tauspon.eq.0.) then
        tauspon=ntauspon*dt
      endif
      if(taut.eq.0.) then
        taut=ntaut*dt
      endif
      if(taum.eq.0.) taum=1.e8

      xl=(nx-1)*dx
      yl=(ny-1)*dy

      nx1=nx+1
      ny1=ny+1
      ns1=ns+1
      nx2=nx+2
      ny2=ny+2
      nxm1=nx-1
      nym1=ny-1
      nall=(nx1+1)*(ny1+1)*(ns1+1)
      nxy=(nx1+1)*(ny1+1)
      nxs=(nx1+1)*(ns1+1)
      nys=(ny1+1)*(ns1+1)
      jy=nx1+1
      js=nxy
      nall3=3*nall
      nallm=nall-nxy
      nallmm=nallm-nxy

      nx1qv=nx1*isqv
      ny1qv=ny1*isqv
      ns1qv=ns1
      nx1qc=nx1*isqc
      ny1qc=ny1*isqc
      ns1qc=ns1
      nx1qr=nx1*isqr
      ny1qr=ny1*isqr
      ns1qr=ns1
	nx1qi=nx1*isqi
      ny1qi=ny1*isqi
      ns1qi=ns1

      if(.not.anamount) then
        xmount=xl/2.
        ymount=yl/2.
      endif

      dt2=2.*dt
      difl=difcof
      xspan=xmount-2.*dx
      yspan=ymount-2.*dy
      ug=ug0
      vg=vg0

      allocate(ugeos(ns))
      allocate(xugeos(ns,5))
      allocate(vgeos(ns))
      allocate(xvgeos(ns,5))

      do is=1,ns-1
         ugeos(is)=ug
         vgeos(is)=vg
      enddo

      if(autogrid) then
        allocate(sigma0(0:ns+1))
        allocate(sigma1(0:ns+1))
        call sigrid(sigma0,sigma1,ns)
      endif

c     dzpro=zthdat(ndat)/npro

      if(restar) then
         open(21,file=fext(fnrst,'rst'),status='old'
     :      ,form='unformatted')
      endif

      if(numsmo.eq.0) then
        dohsmo=.false.
        dovsmo=.false.
        dovsmt=.false.
        numsmo=999999
      endif
c
c output files:
c
c (3) 10 - hsuf,psuf
c (4) 11 - u
c (5) 12 - v
c (6) 13 - w
c (7) 14 - pt
c (8) 15 - phi
c (9) 16 - phis
c (10)17 - wsig
c (11)18 - difunu
c (12)19 - difunv
c (13)30 - fmdata
c     36 - variable upstream profile (pts,etc)
c (14)50 - dragxy
c (15)51 - draguv
c (16)53 - momflu
c (17)41 - rsdata
c (18)98 - output formatted
c (19)52 - local pressure output
c (20)20 - q (specific humidity)
c (21)27 - vertical profiles
c (22)23 - soil parameters (input matrices) solnoii
c (23)31 - surface data (output) vs time
c     32 - profile data (output)
c     33 - tsnoi
c     33 - t2noi
c     33 - wgnoi
c     33 - w2noi
c     33 - wrnoi
c     33 - h (entalpy flux)
c     33 - le (latent heat flux)
c     33 - g  (heat flux in the soil)
c     34 - cdm (momentum transfer coeficient)
c     37 - qc (cloud water)
c     38 - qr (rain water)
c INPUT
c     35 - reference state
c
      if(unfout) then
         lengrc=((nx1+1)*(ny1+1)*8+16)/nbyte_acc
         open(10,file=fext(fnout,'u10'),form='unformatted'
     :      ,access='direct',recl=lengrc,status='unknown')
         open(11,file=fext(fnout,'u11'),form='unformatted'
     :      ,access='direct',recl=lengrc,status='unknown')
         open(12,file=fext(fnout,'u12'),form='unformatted'
     :      ,access='direct',recl=lengrc,status='unknown')
         open(13,file=fext(fnout,'u13'),form='unformatted'
     :      ,access='direct',recl=lengrc,status='unknown')
         open(14,file=fext(fnout,'u14'),form='unformatted'
     :      ,access='direct',recl=lengrc,status='unknown')
         open(15,file=fext(fnout,'u15'),form='unformatted'
     :      ,access='direct',recl=lengrc,status='unknown')
         open(16,file=fext(fnout,'u16'),form='unformatted'
     :      ,access='direct',recl=lengrc,status='unknown')
         open(17,file=fext(fnout,'u17'),form='unformatted'
     :      ,access='direct',recl=lengrc,status='unknown')
         open(18,file=fext(fnout,'u18'),form='unformatted'
     :      ,access='direct',recl=lengrc,status='unknown')
         open(19,file=fext(fnout,'u19'),form='unformatted'
     :      ,access='direct',recl=lengrc,status='unknown')
         open(30,file=fext(fnout,'fmd'),status='unknown')
         open(36,file=fext(fnout,'pts'),status='unknown'
     :     ,access='direct',form='unformatted',recl=2000)
         open(20,file=fext(fnout,'u20'),form='unformatted'
     :      ,access='direct',recl=lengrc,status='unknown')
         open(33,file=fext(fnout,'usf'),form='unformatted'
     :      ,access='direct',recl=lengrc,status='unknown')
         open(34,file=fext(fnout,'cdm'),form='unformatted'
     :      ,access='direct',recl=lengrc,status='unknown')
         open(37,file=fext(fnout,'qc '),form='unformatted'
     :      ,access='direct',recl=lengrc,status='unknown')
         open(38,file=fext(fnout,'qr '),form='unformatted'
     :      ,access='direct',recl=lengrc,status='unknown')
         open(39,file=fext(fnout,'con'),form='unformatted'
     :      ,access='direct',recl=lengrc,status='unknown')
         open(54,file=fext(fnout,'aut'),form='unformatted'
     :      ,access='direct',recl=lengrc,status='unknown')
         open(55,file=fext(fnout,'col'),form='unformatted'
     :      ,access='direct',recl=lengrc,status='unknown')
         open(56,file=fext(fnout,'evp'),form='unformatted'
     :      ,access='direct',recl=lengrc,status='unknown')
         open(57,file=fext(fnout,'vrn'),form='unformatted'
     :      ,access='direct',recl=lengrc,status='unknown')
         open(58,file=fext(fnout,'pre'),form='unformatted'
     :      ,access='direct',recl=lengrc,status='unknown')
         open(47,file=fext(fnout,'dfc'),form='unformatted'
     :      ,access='direct',recl=lengrc,status='unknown')
      endif
c
      if(drgprt) then
         open(50,file=fext(fnout,'dxy'),status='unknown')
         open(51,file=fext(fnout,'duv'),status='unknown')
         open(53,file=fext(fnout,'mom'),status='unknown')
      endif
c
      if(rstout) then
         open(41,file=fext(fnout,'rst'),form='unformatted'
     :      ,status='unknown')
      endif
c
      open(28,file=fext(fnout,'out'),status='unknown')
      open(67,file=fext(fnout,'en1'),status='unknown')
      open(68,file=fext(fnout,'en2'),status='unknown')
      open(69,file=fext(fnout,'en3'),status='unknown')
      open(70,file=fext(fnout,'en4'),status='unknown')
      open(71,file=fext(fnout,'en5'),status='unknown')
      open(52,file=fext(fnout,'loc'),status='unknown')
      write(52,'(a5,20(3a5,a10,3a7,a8))') 'step'
     :  ,('ix','iy','is'
     :  ,'ps','u','v','w','theta',ioutu=1,noutu)
      open(31,file=fext(fnout,'sur'),status='unknown')
!     open(32,file=fext(fnout,'pro'),status='unknown')
c
c Soil parameters:
c ifsoil=0 no soil
c ifsoil=2 h=le=0 momentum flux is computed
c ifsoil=3 h and le are given, momentum flux is computed
c note that the soil model needs qif=1. (water vapour in the air)
c
c
c Water switch (qif=0 no water; qdrag=0. no effect on density)
c else         (qif=1    water; qdrag=1.    effect on density)
c Cloud Water switch (ifqc=0 no cloud water)
c Rain Water switch (ifqr=0 no rain water)
c Initial temperature perturbation switch
c         (ifblob=1 put in initial perturbation see sub. blob)
c
      if(ifsoil.eq.1 .and. qif.eq.0.) then
        write(0,*) 'nh3dFatalError=Aborted: ifsoil=1 needs qif=1.'
        write(nchn1,*) 'Aborted: ifsoil=1 needs qif=1.'
        stop
      endif
      if(ifqr.eq.1 .and. ifqc.ne.1) then
        ifqc=1
        isqc=1
        write(0,*) 'nh3dWarning=Changed: ifqc=1 ! (Because ifqr=1)'
        write(nchn1,*) 'Changed: ifqc=1 ! (Because ifqr=1)'
      endif
      if(ifqr.eq.1 .and. qif.ne.1) then
        qif=1.
        isqv=1
        write(0,*) 'nh3dWarning=Changed: qif=1. ! (Because ifqr=1)'
        write(nchn1,*) 'Changed: qif=1. ! (Because ifqr=1)'
      endif
      if(ifqc.eq.1 .and. qif.ne.1) then
        qif=1.
        isqv=1
        write(0,*) 'nh3dWarning=Changed: qif=1. ! (Because ifqc=1)'
        write(nchn1,*) 'Changed: qif=1. ! (Because ifqc=1)'
      endif
	if(ifqi.eq.1.and.qif.ne.1) then                                           DC,11.2009
	  qif=1.                                                                  DC,11.2009
        isqv=1                                                                  DC,11.2009
        write(0,*) 'nh3dWarning=Changed: qif=1. ! (Because ifqi=1)'             DC,11.2009
        write(nchn1,*) 'Changed: qif=1. ! (Because ifqi=1)'                     DC,11.2009
      endif                                                                     DC,11.2009
	if(ifqi.eq.1.and.ifqc.ne.1) then                                          DC,11.2009
	  ifqc=1.                                                                 DC,11.2009
        isqc=1                                                                  DC,11.2009
        write(0,*) 'nh3dWarning=Changed: ifqc=1. ! (Because ifqi=1)'            DC,11.2009
        write(nchn1,*) 'Changed: ifqc=1. ! (Because ifqi=1)'                    DC,11.2009
      endif                                                                     DC,11.2009
	if(ifqi.eq.1.and.ifqr.ne.1) then                                          DC,11.2009
	  ifqr=1.                                                                 DC,11.2009
        isqr=1                                                                  DC,11.2009
        write(0,*) 'nh3dWarning=Changed: ifqr=1. ! (Because ifqi=1)'            DC,11.2009
        write(nchn1,*) 'Changed: ifqr=1. ! (Because ifqi=1)'                    DC,11.2009
      endif                                                                     DC,11.2009
c
c Switch for update of reference state
c    (ifuprefs=1 update; ifuprefs=0 no update; ifuprefs=2 update 3D)
c
c    (iwind=1 update uv; itheta=1 update temp)
c
c data for analytic temperature cycles at the surface
c

      if(ifuprefs.ne.0) then
        if(.not.ok_fnfor) then
          write(0,*) 'nh3dFatalError=Aborted due to missing data: FNFOR'
          stop
        endif
        if(max(nxsponge1,nxsponge2,nysponge1,nysponge2).lt.0) then
          write(0,*) 'nh3dWarning: the lateral sponge must be turned on'
     :      ,' to force a time-varying profile (nxsponge1....)'
          write(nchn1,*) 'Warning: the lateral sponge must be turned on'
          write(nchn1,*)
     :      'to force a time-varying profile (nxsponge1....)'
        endif

c initial value for interval between forcing profiles: the actual value is
c   given in the profile data

        intrefst=0
c        open(35,file=fext(fnfor,'ref'),status='old')
c        read(35,'(a8)') line
c        read(35,*) ndat,nprof,intrefst
      endif

      if(timeini.ne.starttime .and. ok_fnfor .and.
     :  timeini.ne.timeini0) then
        do while(starttime.lt.timeini)
          write(*,*) 'Looking for initial profiles'
          write(*,*) timeini,starttime
          line='NO'
          do while(line(1:4).ne.'yyyy')
            read(35,'(a)',err=9001,end=9001) line
          enddo
          read(35,*) starttime
          write(*,*) starttime
        enddo
        write(*,*) 'Starting time:',starttime
        tprof1=starttime
        do ip=1,nprof
          if(verbose.ge.2) write(*,*) '1st Profile:',ip,tprof1
          read(35,'(a8)') line
          read(35,*) xprof(ip),yprof(ip)
          read(35,'(a8)') line
          do ii=0,ndat-1
            read (35,*) zthdat(ii,ip),thdat(ii,ip)
     :        ,usdat(ii,ip),vsdat(ii,ip),qvsdat(ii,ip),pressdat(ii,ip)
            if(verbose.ge.3)
     :       write (*,'(6e15.7)') zthdat(ii,ip),thdat(ii,ip)
     :        ,usdat(ii,ip),vsdat(ii,ip),qvsdat(ii,ip),pressdat(ii,ip)
          enddo
          zthdat(ndat,ip)=max(20000.,zthdat(ndat-1,ip))
          usdat(ndat,ip)=usdat(ndat-1,ip)
          vsdat(ndat,ip)=vsdat(ndat-1,ip)
          qvsdat(ndat,ip)=qvsdat(ndat-1,ip)
          thdat(ndat,ip)=thdat(ndat-1,ip)
          pressdat(ndat,ip)=100.  ! not relevant?

        enddo
      endif



c
c nbyte_acc is the number of bytes per unit of direct access length
c it is equal to 1 in the IBM fortran and it equals 4 in the
c DEC Fortran
c This is only relevant for unformatted 3D output!
c
      if(unfout) then
        WRITE(*,*) 'Important !'
        write(*,*) 'nbyte_acc=',nbyte_acc
        write(*,*) 'VERIFY THAT THIS IS OK:'
        write(*,*) 'DIGITAL(4),IBM(1),MICROSOFT(1)'
        WRITE(nchn1,*) 'Important !'
        write(nchn1,*) 'nbyte_acc=',nbyte_acc
        write(nchn1,*) 'VERIFY THAT THIS IS OK:'
        write(nchn1,*) 'DIGITAL(4),IBM(1),MICROSOFT(1)'
      endif
      if(numsmo.eq.0) then
        numsmo=9999999
      endif




c verify need for mass-flux correction at the lateral boundaries (raduv)
!
!     if(max(nxsponge1,nxsponge2,nysponge1,nysponge2).gt.0
!    :  .and. .not.fluxsel) then
!       ifluxcor=0
!       if(verbose.ge.1) then
!         write(*,*) 'No flux correction at boundaries (sponge active)'
!       endif
!     endif

      dpmeddt=0.

      if(ifluxcor.ne.0 .and.max(xlcor,xrcor,ylcor,yrcor).lt.1.) then
        write(0,*) 'nh3dFatalError=Fatal error ifluxcor'
        stop
      endif

! prepare list of available grd files, considering the output options selected

      if(grdout) then
        lstname='var_out'
        open(77,file=fileext(fngrd,lstname))
        if(norout) then
          form1='(a15,1x,a2,1x,a30,1x,a10,a5)'
          write(77,form1) 'theta','th','Potential temperature','K','*'
          write(77,form1) 'theta"','pt','Potential temp perturb','K','*'
          write(77,form1) 'u','um','x-comp of wind','m/s','*'
          write(77,form1) 'v','vm','y-comp of wind','m/s','*'
          write(77,form1) 'w','wm','z-comp of wind','m/s','*'
          if(isqv.eq.1) then
            write(77,form1) 'qv','qv','specific humidity','kg/kg','*'
            write(77,form1) 'co','co','condensation','kg/kg','*'
            write(77,form1) 'ev','ev','evaporation','kg/kg','*'
          endif
          if(isqc.eq.1) then
            write(77,form1) 'qc','qc','cloud water','kg/kg','*'
          endif
          if(isqr.eq.1) then
            write(77,form1) 'qr','qr','rain water','kg/kg','*'
          endif
        endif
        if(morout) then

          write(77,form1) 'u','u_','x-comp of wind (stag)','m/s','*'
          write(77,form1) 'v','v_','y-comp of wind (stag)','m/s','*'
          write(77,form1) 'v','v_','y-comp of wind (stag)','m/s','*'
          write(77,form1) 'phi','ph','geopotential','J/kg','*'
          write(77,form1) 'dsigma/dt','ws','sigma velocity','1/s','*'
          if(isqv.eq.1) then
            write(77,form1) 'rh','rh','relative humidity','NA','*'
          endif
        endif
        if(surout) then
          write(77,form1) 'T_surf','t0','surface air temp','K','s000'
          write(77,form1) 'p_surf','ps','surface pressure','Pa','s000'
          write(77,form1) 'p_anom','px','surface pressure anomaly','Pa'
     :      ,'s000'
          write(77,form1) 'h_surf','hs','surface height','m','s000'
          if(isqr.eq.1) then
            write(77,form1) 'precipitation','pr','precipitation','mm'
     :        ,'s000'
          endif
          if(ifsoil.ne.0) then
            write(77,form1) 'tm_sur','tm','top s+w temperature','K'
     :        ,'s000'
            write(77,form1) 't_soil','ts','top soil temperature','K'
     :        ,'s000'
            write(77,form1) 't_2','t2','deeper soil temperature','K'
     :        ,'s000'
            write(77,form1) 'w_ground','wq','top soil moisture','m/m'
     :        ,'s000'
            write(77,form1) 'w_2','w2','deeper soil moisture','m/m'
     :        ,'s000'
            write(77,form1) 'SensHeat','hh','sensible heat flux','W/m^2'
     :        ,'s000'
            write(77,form1) 'LatHeat','lh','latent heat flux','W/m^2'
     :        ,'s000'
            write(77,form1) 'GrdHeat','gs','ground heat flux','W/m^2'
     :        ,'s000'
            write(77,form1) 'NetRad','rn','net radiation flux','W/m^2'
     :        ,'s000'
            write(77,form1) 'DragCoef','cd','drag coef computed','SI'
     :        ,'s000'
          endif
        endif
        if(refout) then
          write(77,form1) 'phi_s','fs','ref state geopotential','J/kg'
     :        ,'*'
          write(77,form1) 'u_s','us','ref state u-wind','m/s','*'
          write(77,form1) 'v_s','vs','ref state v-wind','m/s','*'
          write(77,form1) 'theta_s','os','ref state pot temp','K','*'
          write(77,form1) 't_ref','tt','ref state temp','K','*'
          write(77,form1) 'p_s_ref','pp','ref state pressure','K','*'
        endif
        if(difout) then
          write(77,form1) 'diff_u','du','diff term in u','m/s^2','*'
          write(77,form1) 'diff_v','dv','diff term in v','m/s^2','*'
          write(77,form1) 'diff_w','dw','diff term in w','m/s^2','*'
          write(77,form1) 'diff_theta','dt','diff term in theta','K/s'
     :        ,'*'
          if(isqv.eq.1) then
            write(77,form1) 'diff_qv','dq','diff term in qv','1/s','*'
          endif
          if(isqc.eq.1) then
            write(77,form1) 'diff_qc','dc','diff term in qc','1/s','*'
          endif
          if(isqr.eq.1) then
            write(77,form1) 'diff_qr','dr','diff term in qr','1/s','*'
          endif
          if(driout) then
            write(77,form1) 'Richardson','ri','Rich number','NA','*'
          endif
          write(77,form1) 'Diff_coef','kk','Diffusion coeficient','NA'
     :        ,'*'
        endif
        close(77)
!        write(0,*) 'var_out.lst ready'


! prepare list of available cross-sections
        lstname='sec_out'
        open(77,file=fileext(fngrd,lstname))
        do iplan=1,nplan
          if(ixp(iplan).gt.0) then
            write(77,'(a1,3i1)') 'x',ixp(iplan)/100
     :        ,mod(ixp(iplan),100)/10,mod(ixp(iplan),10)
          endif
        enddo
        do iplan=1,nplan
          if(iyp(iplan).gt.0) then
            write(77,'(a1,3i1)') 'y',iyp(iplan)/100
     :        ,mod(iyp(iplan),100)/10,mod(iyp(iplan),10)
          endif
        enddo
        do iplan=1,nplan
          if(isp(iplan).gt.0) then
            write(77,'(a1,3i1)') 's',isp(iplan)/100
     :        ,mod(isp(iplan),100)/10,mod(isp(iplan),10)
          endif
        enddo
        do iplan=1,nplan
          if(zplan(iplan).gt.0.) then
            izp=nint(zplan(iplan)/100)
            write(77,'(a1,3i1)') 'z',izp/100
     :        ,mod(izp,100)/10,mod(izp,10)
          endif
        enddo
        close(77)
!        write(0,*) 'sec_out.lst ready'

        lstname='prf_out'
        open(77,file=fileext(fngrd,lstname))
        do iii=1,noutp
          write(77,'(2(a1,3i1))')
     :      'x',ixoutp(iii)/100
     :      ,mod(ixoutp(iii),100)/10,mod(ixoutp(iii),10)
     :      ,'y',iyoutp(iii)/100
     :      ,mod(iyoutp(iii),100)/10,mod(iyoutp(iii),10)
        enddo
        close(77)
!        write(0,*) 'prf_out.lst ready'

        lstname='loc_out'
        open(77,file=fileext(fngrd,lstname))
        do iii=1,noutu
          write(77,'(3(a1,3i1))')
     :      'x',ixoutu(iii)/100
     :      ,mod(ixoutu(iii),100)/10,mod(ixoutu(iii),10)
     :      ,'y',iyoutu(iii)/100
     :      ,mod(iyoutu(iii),100)/10,mod(iyoutu(iii),10)
     :      ,'s',isoutu(iii)/100
     :      ,mod(isoutu(iii),100)/10,mod(isoutu(iii),10)
        enddo
        close(77)
!        write(0,*) 'loc_out.lst ready'

! prepare list of available times at output

        lstname='tim_out'
        open(77,file=fileext(fngrd,lstname))
        write(*,*) 'starttime:',starttime
        do iss=0,ntime,iotfor
          call timestamp(iss,dt,starttime,timestring)
          do iii=1,14
            if(timestring(iii:iii).eq.' ') timestring(iii:iii)='0'
          enddo
          write(77,'(i10,1x,a14)') iss,timestring
        enddo

        if(mod(ntime,iotfor).ne.0) then
          iss=ntime
          call timestamp(iss,dt,starttime,timestring)
          do iii=1,14
            if(timestring(iii:iii).eq.' ') timestring(iii:iii)='0'
          enddo
          write(77,'(i10,1x,a14)') iss,timestring
        endif

        close(77)
        write(0,*) 'tim_grd.lst ready'

      endif

      if(iomtyp.eq.1) then
        yymount=(ny-1)/2.*dy
        if(ymount.ne.yymount) then
          write(*,*) 'warning: ymount adjusted'
          ymount=yymount
        endif
      endif

      dpressdt=0.

      pminimal=sigma0(0)*(pa-ptop)+ptop
      if(pminimal.le.5.) then
        write(0,*) 'nh3dFatalError=ptop is too small',ptop,pminimal
        stop
      endif

      return

9001  write(0,*) 'nh3dFatalError=Profile not found'

      stop
      end

      subroutine weiprof(fngrd)

c compute weighting factors for the different forcing profiles
c case nprof>1

c neglecting staggering of the grid

c using carteian coordinates given in the topography file if available

c if nprof=1 no weighting
c if divrefstate=.false. no horizontal structure, the profiles are interpolated to the central grid point

      use alloc
      use refstate
      character*80 fngrd

      if(allocated(pweight)) deallocate(pweight)
      allocate(pweight(0:nx+1,0:ny+1,nprof))

      if(verbose.ge.2) write(*,*) 'weiprof',nprof,divrefstate

      if(iomtyp.eq.4) then
        xmountc=-xminutm
        ymountc=-yminutm
      else
        xmountc=xmount+1.5*dx
        ymountc=ymount+1.5*dy
      endif

      if(nprof.gt.1 .and. divrefstate) then
        do ix=0,nx+1
        do iy=0,ny+1
          sumdisti=0.
          do ip=1,nprof
            xis=ix*dx-xmountc
            yps=iy*dy-ymountc
            distinv(ip)=1./sqrt((xis-xprof(ip))**2+(yps-yprof(ip))**2)
            sumdisti=sumdisti+distinv(ip)
          enddo
          do ip=1,nprof
            pweight(ix,iy,ip)=distinv(ip)/sumdisti
          enddo
        enddo
        enddo
      elseif(nprof.eq.1) then
        pweight=1.
      else
        xis=(nx/2)*dx-xmountc
        yps=(ny/2)*dy-ymountc
        sumdisti=0.
        do ip=1,nprof
c         write(*,*) ip,xis,xprof(ip),yps,yprof(ip)
          distinv(ip)=1./sqrt((xis-xprof(ip))**2+(yps-yprof(ip))**2)
          sumdisti=sumdisti+distinv(ip)
        enddo
        do ip=1,nprof
          pweight(:,:,ip)=distinv(ip)/sumdisti
        enddo
      endif

      if(refout) then
        do ip=1,nprof
          call wgrids(pweight,'pw',2,nx,2,ny,ip-1,ip-1,1,0,0,0,fngrd)
        enddo
      endif
      return
      end

      subroutine allocref
      use alloc
      use refstate

          if(allocated(zthdat)) deallocate(zthdat)
          allocate(zthdat(0:ndat,nprof))
          if(allocated(thdat)) deallocate(thdat)
          allocate(thdat(0:ndat,nprof))
!          if(allocated(per_dat)) deallocate(per_dat) 
!          allocate(per_dat(0:ndat,nprof))
          if(allocated(zusdat)) deallocate(zusdat)
          allocate(zusdat(0:ndat,nprof))
          if(allocated(usdat)) deallocate(usdat)
          allocate(usdat(0:ndat,nprof))
          if(allocated(zvsdat)) deallocate(zvsdat)
          allocate(zvsdat(0:ndat,nprof))
          if(allocated(vsdat)) deallocate(vsdat)
          allocate(vsdat(0:ndat,nprof))
          if(allocated(psdat)) deallocate(psdat)
          allocate(psdat(0:ndat,nprof))
          if(allocated(ptdat)) deallocate(ptdat)
          allocate(ptdat(0:ndat,nprof))
          if(allocated(zqvsdat)) deallocate(zqvsdat)
          allocate(zqvsdat(0:ndat,nprof))
          if(allocated(qvsdat)) deallocate(qvsdat)
          allocate(qvsdat(0:ndat,nprof))
          if(allocated(pressdat)) deallocate(pressdat)
          allocate(pressdat(0:ndat,nprof))
          
          if(allocated(xprof)) deallocate(xprof)
          allocate(xprof(nprof))

          if(allocated(yprof)) deallocate(yprof)
          allocate(yprof(nprof))

          if(allocated(distinv)) deallocate(distinv)
          allocate(distinv(nprof))

          if(allocated(fpro)) deallocate(fpro)
          allocate(fpro(0:npro+1))

          if(allocated(zpro)) deallocate(zpro)
          allocate(zpro(0:npro))

          if(allocated(fnpro)) deallocate(fnpro)
          allocate(fnpro(0:npro+1))

          if(allocated(press0)) deallocate(press0)
          allocate(press0(nprof))

          if(allocated(press1)) deallocate(press1)
          allocate(press1(nprof))

          if(allocated(dpressdt)) deallocate(dpressdt)
          allocate(dpressdt(nprof))

          if(allocated(thpro0)) deallocate(thpro0)
          allocate(thpro0(0:npro,nprof))
          
          if(allocated(perpro0)) deallocate(perpro0)
          allocate(perpro0(0:npro,nprof))

          if(allocated(thpro1)) deallocate(thpro1)
          allocate(thpro1(0:npro,nprof))

c         if(allocated(thpro2)) deallocate(thpro2)
c         allocate(thpro2(0:npro,nprof))

          if(allocated(dthdz)) deallocate(dthdz)
          allocate(dthdz(0:npro,nprof))
          
          if(allocated(dperdz)) deallocate(dperdz)
          allocate(dperdz(0:npro,nprof))

          if(allocated(dthdz1)) deallocate(dthdz1)
          allocate(dthdz1(0:npro,nprof))

c         if(allocated(dthdz2)) deallocate(dthdz2)
c         allocate(dthdz2(0:npro,nprof))

          if(allocated(dthdt)) deallocate(dthdt)
          allocate(dthdt(0:npro,nprof))

          if(allocated(uspro0)) deallocate(uspro0)
          allocate(uspro0(0:npro,nprof))

          if(allocated(uspro1)) deallocate(uspro1)
          allocate(uspro1(0:npro,nprof))

c         if(allocated(uspro2)) deallocate(uspro2)
c         allocate(uspro2(0:npro,nprof))

          if(allocated(dusdz)) deallocate(dusdz)
          allocate(dusdz(0:npro,nprof))

          if(allocated(dusdz1)) deallocate(dusdz1)
          allocate(dusdz1(0:npro,nprof))

c         if(allocated(dusdz2)) deallocate(dusdz2)
c         allocate(dusdz2(0:npro,nprof))

          if(allocated(dusdt)) deallocate(dusdt)
          allocate(dusdt(0:npro,nprof))

          if(allocated(vspro0)) deallocate(vspro0)
          allocate(vspro0(0:npro,nprof))

          if(allocated(vspro1)) deallocate(vspro1)
          allocate(vspro1(0:npro,nprof))

c         if(allocated(vspro2)) deallocate(vspro2)
c         allocate(vspro2(0:npro,nprof))

          if(allocated(dvsdz)) deallocate(dvsdz)
          allocate(dvsdz(0:npro,nprof))

          if(allocated(dvsdz1)) deallocate(dvsdz1)
          allocate(dvsdz1(0:npro,nprof))

c         if(allocated(dvsdz2)) deallocate(dvsdz2)
c         allocate(dvsdz2(0:npro,nprof))

          if(allocated(dvsdt)) deallocate(dvsdt)
          allocate(dvsdt(0:npro,nprof))

          if(allocated(qvspro0)) deallocate(qvspro0)
          allocate(qvspro0(0:npro,nprof))

          if(allocated(qvspro1)) deallocate(qvspro1)
          allocate(qvspro1(0:npro,nprof))

c         if(allocated(qvspro2)) deallocate(qvspro2)
c         allocate(qvspro2(0:npro,nprof))

          if(allocated(dqvsdz)) deallocate(dqvsdz)
          allocate(dqvsdz(0:npro,nprof))

          if(allocated(dqvsdz1)) deallocate(dqvsdz1)
          allocate(dqvsdz1(0:npro,nprof))

c         if(allocated(dqvsdz2)) deallocate(dqvsdz2)
c         allocate(dqvsdz2(0:npro,nprof))

          if(allocated(dqvsdt)) deallocate(dqvsdt)
          allocate(dqvsdt(0:npro,nprof))


          if(allocated(dthdzdt)) deallocate(dthdzdt)
          allocate(dthdzdt(0:npro,nprof))

          if(allocated(dusdzdt)) deallocate(dusdzdt)
          allocate(dusdzdt(0:npro,nprof))

          if(allocated(dvsdzdt)) deallocate(dvsdzdt)
          allocate(dvsdzdt(0:npro,nprof))

          if(allocated(dqvsdzdt)) deallocate(dqvsdzdt)
          allocate(dqvsdzdt(0:npro,nprof))

          if(allocated(presspro)) deallocate(presspro)
          allocate(presspro(0:npro,nprof))

          if(allocated(tepro)) deallocate(tepro)
          allocate(tepro(0:npro,nprof))

          if(allocated(rhopro)) deallocate(rhopro)
          allocate(rhopro(0:npro,nprof))
          if(allocated(rhopro1)) deallocate(rhopro1)
          allocate(rhopro1(0:npro,nprof))

          if(allocated(tedat)) deallocate(tedat)
          allocate(tedat(0:ndat,nprof))

      return
      end

      subroutine defaults(fnmap,fnfor,fnout
     :   ,itheta,iwind,scrout,fngrd,nbyte_acc)
c
c default values for parameters
c
      use alloc
      use nh3dpa2
      use refstate
      implicit real*8(a-h,o-z)

      character*80 fnrst,fnmap,fnfor,fnout,fext,fngrd
      logical scrout
      external fext
c
c      include 'nh3dpa10.inc'
c
      isqv=0
      isqc=0
      isqr=0
      isqi=0

      timeini=00000101000000
      starttime=timeini
      realtime=.false.
      intrefst=1
      scrout=.true.
      iotscrout=1
      ifluxcor=1

      iomtyp=4
      hmount=0.
      ivarhs=0
      hexp=1.
      nstepgrow=0
      xlatit=0.
      fcori=0.
      ptop=200.
      pa=100000. !1.e5        ! dc 02.10
      ha=0.

      nupsur=1
      nuppts=1
      nupdif=1

      adjust=.false.
      hydro=.false.
      tfct=.false.
      inipt=.false.
      raylei=.true.
      olddra=.false.
      xlcor=1.
      xrcor=1.
      ylcor=1.
      yrcor=1.

      iophis=3
      ioppuv=2
      iowlbc=0
      iodify=1
      iobphy=1
      iobuux=1
      iobuuy=1
      iobvvx=1
      iobvvy=1
      iobptx=1
      iobpty=1
      iobppx=1
      iobppy=1

      perdam=.true.
      dohsmo=.true.
      dovsmo=.true.
      dovsmt=.true.

      numsmo=5
      hdampm=0.015
      hdampt=0.0625
      vdamp=0.001
      vdampt=0.001

      iodif=0
      rikey=1.
      difcof=0.
      uvdif=1.
      tsdif=1.

      cdcoef=0.

      taum=1.e8
      ntaut=10
      zm=8000.
      zt=18000.
      nxsponge1=-1
      nxsponge2=-1
      nysponge1=-1
      nysponge2=-1
      ntauspon=ntaut

      exptit='test_001'
      comlin(1)='NH3D TEST'
      comlin(2)=' '
      comlin(3)=' '

      restar=.false.

      rstout=.false.
      iotrst=999999

      forout=.true.
      iotfor=999999

      unfout=.false.
      iotunf=999999
      iotsart=0
      ioend=-1
      iodelta=150

      drgprt=.false.
      domom=.false.
      iotdrg=1
      iotduv=20
      nlev=-999
      dzlev=500.

      mxspan=.false.

      iomod=0
      proout=.true.
      norout=.true.
      surout=.true.
      masout=.false.
      iniout=.true.
      iniunf=.false.

      echeck=.false.
      eneout=.false.

      morout=.false.
      refout=.false.
      inrout=.false.
      fluout=.false.
      ppbout=.false.
      duvout=.false.
      uubout=.false.
      ptbout=.false.

      difout=.false.
      driout=.false.
      ddiout=.false.
      defout=.false.

      eldout=.false.
      elfout=.false.
      elgout=.false.
      phsout=.false.
      desout=.false.

      ixss=1
      iyss=1
      isss=1
      fixexp=.false.
      nexfix=0

      ixp=-1
      iyp=-1
      isp=-1
      zplan=-1.

      ug0=0.
      dugdz0=0.
      vg0=0.
      dvgdz0=0.
      itepsuf=.true.

      fnmap='test_map'
      fnfor='test_for'
      fnout='test_001'
      fngrd='test'
      fnrst='test_rst'

      ifsoil=0
      z0hz0=1.

      qif=0.
      qdrag=1.
      ifqc=0
      ifqr=0
	ifqi=0                                                                     DC,10.2009
      ifblob=0

      ifuprefs=0
      divrefstate=.false.
      iwind=0
      itheta=0

      filt=0.1
      pfilt=0.1

      nbyte_acc=1

      grdout=.true.
      scalew=1.
      zconst=.true.
      verbose=0
      dzpro=100.

      forcedpp=.false.
      forcepsuf=.false.
      return
      end


      function uppercase(string)
c
c converts string to uppercase
c
      character (len=*) string
      character*20 uppercase
      integer i,ismall,ibig

      uppercase=string
      ismall=ichar('a')
      ibig=ichar('A')
      do i=1,len(string)
        if(string(i:i).ge.'a' .and. string(i:i).le.'z') then
          uppercase(i:i)=char(ichar(string(i:i))+ibig-ismall)
        endif
      enddo
      return
      end




      subroutine refprofs
c-----------------------------------------------------------------------
c  specify vertical profile of reference variables
c  this subrout interpolates a profile. To speed
c  up the process the data is interpolated to a regularly spaced
c  grid. to avoid errors, all the points in the given profile should
c  belong to the regular grid.
c
c ioreft: 1 - interpolate thdat as pts;
c         2 - interpolate thdat as tems;
c         3 - interpolate thdat as N**2 and calculates pts profile.
c
c iorefq  1 - interpolate qvdat as relative humidity (0..1);
c         2 - interpolate qvdat as specific humidity (kg/kg);
c
c ptdat is a vertical profile of potential temperature calculated
c in the case ioreft.ne.2 to compute the pressure profile.
c
c note: the interpolation is done using z as the vertical coordinate.
c in fact, reference profiles should be a function of p.
c-----------------------------------------------------------------------
c p.m. 2005/10/11
c-----------------------------------------------------------------------

      use alloc
      use refstate
      implicit none
      integer i,ip,id,iz
      double precision bvf2m,ptmed,qsat,tid,tmed
      character*80 fndebug
      integer istepp
      save istepp
      data istepp/0/

c initialization of the profile: zdat is (in principle) irregularly
c spaced; zpro is regularly spaced. this should be done only once.

      if(verbose.ge.2) write(*,*) 'Refprofs ',npro,dzpro,nprof,ndat
     :  ,pts0
      if(istepp.eq.0) then
        do i=0,npro
          zpro(i)=i*dzpro
        enddo
c compute weighting factors for profiles

        fndebug='weiprof.grd'
        call weiprof(fndebug)
        istepp=1
      endif

      do ip=1,nprof
c NEW: the surface pressure is read from the ref file

      psref=pressdat(0,ip)
      press0(ip)=pressdat(0,ip)
      if(ioreft.eq.1) then
        call inipro(thdat(0,ip),zthdat(0,ip),ndat,fpro(0),zpro,npro)
        t0ini=thdat(0,ip)*(psref/p00)**akapa
        if(verbose.ge.2) write(*,*) 't0ini=',t0ini
        do id=0,ndat
          ptdat(id,ip)=thdat(id,ip)
        enddo
      elseif(ioreft.eq.2) then
        call inipro(thdat(0,ip),zthdat(0,ip),ndat,fpro(0),zpro,npro)
	
        t0ini=thdat(0,ip)
        if(verbose.ge.3) then
          do i=0,npro
            write(*,*) 'th:',i,fpro(i)
          enddo
        endif
        do id=0,ndat
          ptdat(id,ip)=thdat(id,ip)*(pressdat(id,ip)/1.e5)**(-akapa)
		enddo
	
      else
        call inipro(thdat(0,ip),zthdat(0,ip),ndat,fnpro(0),zpro,npro)
        t0ini=pts0*(psref/p00)**akapa
        if(verbose.ge.2) write(*,*) 't0ini=',t0ini,fnpro(1)
c
        fpro(0)=pts0
        do iz=1,npro
          bvf2m=0.5*(fnpro(iz-1)+fnpro(iz))
          fpro(iz)=fpro(iz-1)*dexp(bvf2m*dzpro/g)
c         if(verbose.ge.2) write(*,*) 'th:',iz,bvf2m,fpro(iz)
        enddo
        fpro(npro+1)=fpro(npro)
c
        ptdat(0,ip)=t0ini
        do id=1,ndat
          bvf2m=0.5*(thdat(id-1,ip)+thdat(id,ip))
          ptdat(id,ip)=ptdat(id-1,ip)
     :      *dexp(bvf2m*(zthdat(id,ip)-zthdat(id-1,ip))/g)
        enddo
      endif

      do i=0,npro
        dthdz(i,ip)=(fpro(i+1)-fpro(i))/dzpro
      enddo
c
      
	do i=0,npro
        thpro0(i,ip)=fpro(i)-dthdz(i,ip)*zpro(i)
      enddo

        if(verbose.ge.3) then
          do i=0,npro
            write(*,*) 'th:',i,thpro0(i,ip),dthdz(i,ip),zpro(i)
          enddo
        endif
c
      call uprof(iorefu,usdat(0,ip),zusdat(0,ip),ndat,fpro(0),zpro,npro)
c
      do i=0,npro
        dusdz(i,ip)=(fpro(i+1)-fpro(i))/dzpro
      enddo
c
      do i=0,npro
        uspro0(i,ip)=fpro(i)-dusdz(i,ip)*zpro(i)
      enddo
c
      call uprof(iorefv,vsdat(0,ip),zvsdat(0,ip),ndat,fpro(0),zpro,npro)
c
      do i=0,npro
        dvsdz(i,ip)=(fpro(i+1)-fpro(i))/dzpro
      enddo
c
      do i=0,npro
        vspro0(i,ip)=fpro(i)-dvsdz(i,ip)*zpro(i)
      enddo
      
c compute the pressure profile (integrate hydrostatic relation):

      psdat(0,ip)=psref
      if(verbose.ge.1) write(*,*) 'psref=',psref
	
      if(ioreft.eq.2) then

        if(iorefq.eq.1) then
          qvsdat(0,ip)=qvsdat(0,ip)*qsat(thdat(0,ip),psdat(0,ip))
        endif
        do id=1,ndat
          tmed=0.5*(thdat(id,ip)+thdat(id-1,ip))
          psdat(id,ip)=psdat(id-1,ip)*exp(-g/(r*tmed)
     :      *(zthdat(id,ip)-zthdat(id-1,ip)))

          if(iorefq.eq.1) then
            qvsdat(id,ip)=qvsdat(id,ip)
     :        *qsat(thdat(id,ip),psdat(id,ip))
          endif
          if(verbose.ge.2) write(*,*)
     :      id,' psdat=',psdat(id,ip),tmed,qvsdat(id,ip)
        enddo

      else
        tid=ptdat(0,ip)*(psdat(0,ip)/1.e5)**akapa
        if(iorefq.eq.1) then
          qvsdat(0,ip)=qvsdat(0,ip)*qsat(tid,psdat(0,ip))
        endif
        do id=1,ndat
          ptmed=0.5*(ptdat(id,ip)+ptdat(id-1,ip))
          psdat(id,ip)=(psdat(id-1,ip)**akapa-g/cp*p00**akapa
     :       *(zthdat(id,ip)-zthdat(id-1,ip))/ptmed)**(1./akapa)
          tid=ptdat(id,ip)*(psdat(id,ip)/1.e5)**akapa
          if(iorefq.eq.1) then
            qvsdat(id,ip)=qvsdat(id,ip)*qsat(tid,psdat(id,ip))
          endif
          if(verbose.ge.2) then
            write(*,*) id,' psdat=',psdat(id,ip),ptmed,qvsdat(id,ip)
     :        ,tid
          endif
        enddo
      endif

! compute temperature profile (needed for bc flux correction nhad0613ff)

      do id=0,ndat
        tedat(id,ip)=ptdat(id,ip)*(psdat(id,ip)/1.e5)**(akapa)
      enddo

      call inipro(tedat(0,ip),zthdat(0,ip),ndat,tepro(0,ip),zpro,npro)

      presspro(0,ip)=psref
	
      do i=1,npro
        presspro(i,ip)=presspro(i-1,ip)*exp(-g*(zpro(i)-zpro(i-1))
     :    /(r*0.5*(tepro(i,ip)+tepro(i-1,ip))))
	
      enddo
      do i=0,npro
        rhopro(i,ip)=presspro(i,ip)/(r*tepro(i,ip))
      enddo

      	
      call inipro(qvsdat(0,ip),zqvsdat(0,ip),ndat,fpro(0),zpro,npro)
c
      do i=0,npro
        dqvsdz(i,ip)=(fpro(i+1)-fpro(i))/dzpro
      enddo
c
      do i=0,npro
        qvspro0(i,ip)=fpro(i)-dqvsdz(i,ip)*zpro(i)
      enddo
c
      write(nchn1,1000) ioreft,pts0,iorefu,iorefv,psref,ip
      do i=0,ndat
        write(nchn1,1001) zthdat(i,ip),thdat(i,ip)
     :    ,zusdat(i,ip),usdat(i,ip)
     :    ,zvsdat(i,ip),vsdat(i,ip),ptdat(i,ip)
     :    ,zqvsdat(i,ip),qvsdat(i,ip)
      enddo

      if(verbose.ge.2) then
        write(*,*) 'Interpolated profile ',ip
        write(*,'(9a15)') 'zpro','th','us','vs','qvs','dthdz'
     :    ,'temp','press','rho'
        do i=0,npro
          write(*,'(i5,9e15.7)')
     :      i,zpro(i),thpro0(i,ip),uspro0(i,ip),vspro0(i,ip)
     :        ,qvspro0(i,ip),dthdz(i,ip)
     :        ,tepro(i,ip),presspro(i,ip),rhopro(i,ip)
        enddo
      endif

      enddo
c
      return
1000  format(1x,'reference profiles: ioreft=',i2,' pts0=',f8.2
     :   ,' iorefu=',i2,' iorefv=',i2,' psref=',f7.0,' prof=',i3//
     :   ,t2,'zthdat',t10,'thdat',t25,'zusdat',t35,'usdat'
     :   ,t50,'zvsdat',t60,'vsdat',t75,'ptdat'//)
1001  format(1x,t2,f8.0,t10,e12.5,t25,f8.0,t35,f10.3,t50,f8.0,t60,f10.3
     :   ,t75,f10.3,t90,f10.3,t112,f10.5)
      end



      subroutine refs
c-----------------------------------------------------------------------
c re-evaluate reference state variables and compute phi at surface
c with pts prespecified
c
c dirichlet lateral boundary condition for phisuf.
c-----------------------------------------------------------------------
c
      use alloc
      use refstate

      implicit real*8(a-h,o-z)
c      include 'nh3dpa10.inc'

c
c reevaluate variables of reference state from new pp(i,3)
c
      if (mod(nstep,nupsur).eq.0 .or. ifuprefs.ne.0) then
c
         if(verbose.ge.2) write(*,*) 'refs:sur'
         if(verbose.ge.3) write(*,'(10i4)') (ifax(k19),k19=1,10)
         if(iophis.eq.1) then
            call surhyd
         elseif(iophis.eq.2) then
            call surphi
         else
            call surphn
         endif
c
         do iy=1,ny1
         do ix=1,nx1
           phig(ix,iy)=g*hsuf(ix,iy)-phisuf(ix,iy)
         enddo
         enddo
c
c calculates phis integrating hydrostatic relation:
c
         do iy=0,ny1
         do ix=0,nx1
             phis(ix,iy,ns)=phisuf(ix,iy)-dsng*r*pp(ix,iy,3)
     :       *(.75*tems(ix,iy,ns)+.25*tems(ix,iy,ns-1))
     :       /(ptop+pp(ix,iy,3)*(1.+0.5*dsng))
         enddo
         enddo
         do is=ns-1,0,-1
         do iy=0,ny1
         do ix=0,nx1
           phis(ix,iy,is)=phis(ix,iy,is+1)+ds0(is+1)*r05*pp(ix,iy,3)
     :       *(tems(ix,iy,is+1)+tems(ix,iy,is))
     :       /(ptop+pp(ix,iy,3)*sigma1(is+1))
         enddo
         enddo
         enddo
      endif
c
c updates reference temperature (interpolating for actual grid):
c
      if (mod(nstep,nuppts).eq.0 .and. ifuprefs.eq.0) then
         if(verbose.ge.2) write(*,*) 'refs:upuvt'
         call upuvt
      endif
c-----------------------------------------------------------------------
c     prepare parameter field 's'
c-----------------------------------------------------------------------
      do is=1,ns
      do iy=1,ny1
      do ix=1,nx1
        s(ix,iy,is)=gr2*(ptop+sigma1(is)*pp(ix,iy,3))
     :    /((tems(ix,iy,is)+tems(ix,iy,is-1))*pp(ix,iy,3))
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
      return
      end






      subroutine soil(fnout,timeho)
      
      use alloc
	implicit real*8(a-h,o-z)
	real(8):: h10_1,l10_1,ls10_1,hs10_1,Ts0_1
     &          Tb0_1, 
     &          Tbb0_1,h_ML0_1,extwat_1,extice_1,kor_1,trib_inflow_1,
     &         Sals0_1,Salb0_1,fetch_1,phi_1,lam_1,us0_1,vs0_1,Tm_1,
     &         alphax_1,alphay_1,a_veg_1,c_veg_1,h_veg_1,area_lake_1 
	 integer(8):: init_T_1
	
      
c      include 'nh3dpa10.inc'

      real(8), dimension(nx1,ny1):: usuf,vsuf                                   VS,05.2007
      character*80 fnout,fext,fnmap
      dimension wgw(2),w2w(2),wrw(2),tsw(2),t2w(2)
      real(8), dimension(:,:),allocatable::                                     VS,05.2007
     & hw_2d,xlew_2d,cdmw_2d,tsw_2d                                             VS,05.2007
	real(8), parameter:: stefan = 5.6697e-08                                  VS,06.2007
	
      character(len=60), save :: outpath	

      save ifirst
      save hw_2d,xlew_2d,cdmw_2d,tsw_2d                                         VS,05.2007

      data ifirst/0/
      data vegw/0./,isandw/1/,ivegw/1/,xlaiw/1./,rsminw/30./
     :   ,emisw/1./,rraw/1./,xdd2w/1./,w2w/1.,1./,wgw/1.,1./

      if(ifirst.eq.0) then
        pi180=4.*atan(1.)/180.
        allocate (hw_2d(nx1,ny1),xlew_2d(nx1,ny1),                              VS,05.2007
     &   cdmw_2d(nx1,ny1),tsw_2d(nx1,ny1))                                      VS,05.2007
!        call INIT_LAKE(nx+2,ny+2,fnmap,nlakecall*dt,
!     &                 h10_1,l10_1,ls10_1,hs10_1,Ts0_1,
!     &          Tb0_1, 
!     &          Tbb0_1,h_ML0_1,extwat_1,extice_1,kor_1,trib_inflow_1,
!     &         Sals0_1,Salb0_1,fetch_1,phi_1,lam_1,us0_1,vs0_1,Tm_1,
!     &         alphax_1,alphay_1,a_veg_1,c_veg_1,h_veg_1,area_lake_1,
!     &         init_T_1	  )
	  
!        outpath = 'results/lake/'
        print*, 'The soil is initialized'
      endif
	  
c
c soil balance
c
      li=3-li
      lf=3-lf

c     do iy=2,ny
c     do ix=2,nx
c        write(99,111) 'fi,fis',phi(ix,iy,ns-1),phis(ix,iy,ns-1),ix,iy
c 111    format(a6,3x,f15.8,3x,f15.8,4x,2i3)
c     enddo
c     enddo

      if(verbose.ge.2) write(*,*) 'Call rad'
	if (radpar==1) then
        call rad(rg,rat,pl,xrg,xra,xpl,
     :    nx,ny,ndados,nf,dt,nstep,
     :    xugeos,xvgeos,ugeos,vgeos,ns,zonirr)
      else
	  pl = 0.d0
	endif  

      do iy=2,ny
      do ix=2,nx
c       deltaz(ix,iy)=(phis(ix,iy,ns-1)+phi(ix,iy,ns-1))/g-hsuf(ix,iy)     estava +phis
         deltaz(ix,iy)=0.5*(phis(ix,iy,ns-1)-phis(ix,iy,ns)
     :    +phi(ix,iy,ns-1)-phi(ix,iy,ns))/g
c        write(99,*) 'deltazmain=',deltaz(ix,iy)
c        write(99,*)'fi(ix,iy,ns-1)fis',phi(ix,iy,ns-1),phis(ix,iy,ns-1)
         if(deltaz(ix,iy).lt.0.1) then
          write(*,*) 'DELTAZ:',ix,iy,deltaz(ix,iy)
     :    ,phis(ix,iy,ns-1),phi(ix,iy,ns-1),hsuf(ix,iy),phisuf(ix,iy)
     :    ,phis(ix,iy,ns),phi(ix,iy,ns)
     :    ,phig(ix,iy)
          write(0,*) 'nh3dFatalError=Deltaz<0.1 Unstable!'
          stop
        endif
c
c   uvsuf: wind speed at the first model level
c
        uvsuf(ix,iy)=dsqrt((0.5*(u(ix,iy,ns-1,3)+u(ix-1,iy,ns-1,3)))**2
     :     +(0.5*(v(ix,iy,ns-1,3)+v(ix,iy-1,ns-1,3)))**2)
     
        usuf(ix,iy) = 0.5*(u(ix,iy,ns-1,3)+u(ix-1,iy,ns-1,3))                   VS,05.2007
        vsuf(ix,iy) = 0.5*(v(ix,iy,ns-1,3)+v(ix,iy-1,ns-1,3))                   VS,05.2007
c
c taixiy: air temperature on the first model level
c
        p=sigma0(ns-1)*pp(ix,iy,2)+ptop
        taixiy=(p/1.e5)**akapa*(pts(ix,iy,ns-1)+pt(ix,iy,ns-1,3))
c       write(*,'(2i5,3e14.7)') ix,iy,taixiy,p,pts(ix,iy,ns-1)
        if(ifsoil.eq.1) then
          qa=qv(ix,iy,ns-1,2)
        else
          qa=1.e-3
        endif
        
c radiation, geostrophic wind and precipitation update
c
      if(ifsoil.eq.1) then
        if (radpar==2) then                                                     VS,05.2007
         tclds = 0.d0                                                           VS,05.2007
         call rad_par2(taixiy,qa,psuf(ix,iy),tclds,                             VS,05.2007
     &    rat(ix,iy),rg(ix,iy))                                                 VS,05.2007
!        tclds is the cloud fraction in the grid cell: temporarily 0.d0        
        elseif (radpar==3) then                                                 VS,05.2007
         rg(ix,iy)  = Srad_surf(ix,iy)                                          VS,26.06.2007
         rat(ix,iy) = Lrad_surf(ix,iy)                                          VS,26.06.2007
     &    + emis(ix,iy)*stefan*tsurf(ix,iy)**4                                  VS,06.2007
        endif                                                                   VS,05.2007
        if(verbose.ge.2) write(*,*) 'Rad done'
      endif
c
        ts(li)=tsnoi(ix,iy)
        tsw(li)=tsurw(ix,iy)
        t2(li)=t2noi(ix,iy)
        wg(li)=wgnoi(ix,iy)
        w2(li)=w2noi(ix,iy)
        wr(li)=wrnoi(ix,iy)
c
c
        if (ifirst.eq.0) then
          ts(lf)=tsnoi(ix,iy)
          tsw(lf)=tsurw(ix,iy)
          t2(lf)=t2noi(ix,iy)
          wg(lf)=wgnoi(ix,iy)
          w2(lf)=w2noi(ix,iy)
          wr(lf)=wrnoi(ix,iy)
        endif

        if (xlake(ix,iy).gt.0.) then
         if (iflake == 0) then                                                  VS,05.2007
          tslake(ix,iy)=tlake(ix,iy,timeho)
c         write(*,*) 'Call solonoi_0',ix,iy
          call solonoi (-2,isandw,ivegw,vegw
     &      ,xlaiw,rsminw,alb(ix,iy)
     &      ,emisw,z0w(ix,iy),z0hw(ix,iy),rraw
     &      ,taixiy,qa,uvsuf(ix,iy),rgw,ratw,psuf(ix,iy)
     &      ,pl(ix,iy),deltaz(ix,iy)
     &      ,tsw,t2w,wgw,w2w,wrw
     &      ,rnw,hw,xlew,gsolow,dtl,li,lf
     &      ,xdd2w,cdmw,rhoa,tslake(ix,iy))
     
         elseif (iflake == 1) then                                              VS,05.2007
          if (mod(nstep,nlakecall)==0.or.nstep==1) then                         VS,05.2007
           trib_inflow = 0.                                                     VS,05.2007
!           call Lake (taixiy,qa,psuf(ix,iy),usuf(ix,iy),                        VS,05.2007
!     &      vsuf(ix,iy),rat(ix,iy),rg(ix,iy),pl(ix,iy),deltaz(ix,iy),           VS,05.2007
!     &      tsw_2d(ix,iy),hw_2d(ix,iy),xlew_2d(ix,iy),cdmw_2d(ix,iy),           VS,05.2007
!     &      trib_inflow,nlakecall*dt,ix,iy,nx+2,ny+2)                           VS,05.2007
            fictfalse = .FALSE.
            stop
!	      call LAKE (taixiy, qa, psuf(ix,iy), usuf(ix,iy),
!     &      vsuf(ix,iy), rat(ix,iy), rg(ix,iy), pl(ix,iy), 1.d-7,
!     :      deltaz(ix,iy),nlakecall*dt,
!     :      h10_1, watice(ix,iy), ls10_1, hs10_1, Ts0_1, Tb0_1, Tbb0_1,
!     :      h_ML0_1, extwat_1, extice_1,kor_1,trib_inflow_1,
!     :      Sals0_1, Salb0_1, fetch_1, phi_1, lam_1, us0_1, vs0_1,
!     :      Tm_1, alphax_1, alphay_1, a_veg_1, c_veg_1, h_veg_1,
!     :      area_lake_1,	
!     :      tsw_2d(ix,iy), hw_2d(ix,iy), xlew_2d(ix,iy),
!     :      cdmw_2d(ix,iy),  ix, iy, nx + 2, ny + 2, 2008,
!     :      imonth, iday, real(ihour), init_T_1,fictfalse, fictfalse,
!     :      outpath)

           endif                                                                 VS,05.2007
          hw =      hw_2d(ix,iy)                                                VS,05.2007
          xlew =    xlew_2d(ix,iy)                                              VS,05.2007
          cdmw =    cdmw_2d(ix,iy)                                              VS,05.2007
          tsw(lf) = tsw_2d(ix,iy)                                               VS,05.2007
         else                                                                   VS,05.2007
          print*, 'Invalid lake model identifier!'                              VS,05.2007
          stop                                                                  VS,05.2007
         endif                                                                  VS,05.2007
	  else
         hw=0.
         xlew=0.
         cdmw=0.
        endif
	 
c
c       write(*,*) ix,iy
c     & ,tsw,t2w,wgw
c     &    ,rnw,hw,xlew
c     &    ,cdmw,tslake(ix,iy)
c
        
        if (xlake(ix,iy).lt.1.) then
          if(iclay(ix,iy).eq.-1) then
            tsfunc=tsoil(ix,iy,timeho)
          endif
c         write(*,*) 'Call solonoi_1',ix,iy,iveg(ix,iy)
          call solonoi(iclay(ix,iy),isand(ix,iy),iveg(ix,iy),veg(ix,iy)
     &      ,xlai(ix,iy),rsmin(ix,iy),alb(ix,iy)
     &      ,emis(ix,iy),z0(ix,iy),z0h(ix,iy),rra(ix,iy)
     &      ,taixiy,qa,uvsuf(ix,iy),rg(ix,iy),rat(ix,iy),psuf(ix,iy)
     &      ,pl(ix,iy),deltaz(ix,iy)
     &      ,ts,t2,wg,w2,wr
     &      ,rn(ix,iy),h(ix,iy),le(ix,iy),gsolo(ix,iy),dtl,li,lf
     &      ,xdd2(ix,iy),cdm(ix,iy),rhoa,tsfunc)
        else
          h(ix,iy)=0.
          le(ix,iy)=0.
          cdm(ix,iy)=0.
        endif
        h(ix,iy)=h(ix,iy)*(1.-xlake(ix,iy))+hw*xlake(ix,iy)
        le(ix,iy)=le(ix,iy)*(1.-xlake(ix,iy))+xlew*xlake(ix,iy)
        cdm(ix,iy)=cdm(ix,iy)*(1.-xlake(ix,iy))+cdmw*xlake(ix,iy)
	  ust_s(ix,iy)=sqrt(cdm(ix,iy))*sqrt(uvsuf(ix,iy))
	 
c
c ifsoil=2 means no heat or moisture flux, but there is momentum flux.
c

        if(ifsoil.eq.2) then
          h(ix,iy)=0.
          le(ix,iy)=0.
        elseif(ifsoil.eq.3) then
          if(ix.gt.nx/2) then
            h(ix,iy)=100.
          else
            h(ix,iy)=-50
          endif
          le(ix,iy)=0.
        else

c  cdm=Cd*u
          tsnoi(ix,iy)=ts(lf)
          tsurw(ix,iy)=tsw(lf)
          tsurf(ix,iy)=(1.-xlake(ix,iy))*ts(lf)+                                VS,06.2007
     &                    xlake(ix,iy)*tsw(lf)                                  VS,06.2007
          t2noi(ix,iy)=t2(lf)
          wgnoi(ix,iy)=wg(lf)
          w2noi(ix,iy)=w2(lf)
          wrnoi(ix,iy)=wr(lf)
          temmin(ix,iy)=min(temmin(ix,iy),tsnoi(ix,iy))
          temmax(ix,iy)=max(temmax(ix,iy),tsnoi(ix,iy))
          evapor(ix,iy)=evapor(ix,iy)+le(ix,iy)/2.5e6*dt
        endif
      enddo
      enddo

!	print*, 'rg',rg
!	print*, 'ra',rat
!	read*
!      print*, 'The first call of soil finished: STOP'                           VS,15.08.2007
!      STOP                                                                      VS,15.08.2007
      continue

      if(ifirst.eq.0) then
        xmin=(2-1.5)*dx-xmount
        ymin=(2-1.5)*dy-ymount
        xmax=(nx-1.5)*dx-xmount
        ymax=(ny-1.5)*dy-ymount

        call wrigrd(wgnoi,xmin,xmax,ymin,ymax,ttt1,ttt2
     :     ,0,nx1,0,ny1,2,nx,2,ny,fext(fnout,'wga'))
        call wrigrd(h,xmin,xmax,ymin,ymax,ttt1,ttt2
     :     ,0,nx1,0,ny1,2,nx,2,ny,fext(fnout,'h_a'))
        call wrigrd(le,xmin,xmax,ymin,ymax,ttt1,ttt2
     :     ,0,nx1,0,ny1,2,nx,2,ny,fext(fnout,'lea'))
        ifirst=1
      endif

c      write(*,'(3(f10.3,1x))') tsnoi(5,5),wgnoi(5,5),rg(5,5)

      return
      end






      subroutine solonoii(fnmap)
      use alloc
      implicit real*8(a-h,o-z)
c      include 'nh3dpa10.inc'
      character*80 fnmap,fext
c
c   esta rotina tem que ser reescrita para dar conta das heterogeneidades
c    da superficie!!!!!
c
c     open (1,file=fext(fnmap,'sol'),status='old')
c     read (1,*) xdd2i
c     read (1,*) ivegi
c     read (1,*) isoli
c     read (1,*) ts0
c     read (1,*) t20
c     read (1,*) wg0
c     read (1,*) w20
c     read (1,*) vegi
c     read (1,*) z0i
c     read (1,*) z0wi
c     read (1,*) z0hi
c     read (1,*) xlaii
c     read (1,*) rsmini
c     read (1,*) albi
c     read (1,*) emisi
c     close (1)
c----------------------------------------------------------------------
c
c
      open (1,file=fext(fnmap,'xdd'),status='old')
      call readgrd(1,xdd2,0,nx+1,0,ny+1)
      close(1)
      open (1,file=fext(fnmap,'ts '),status='old')
      call readgrd(1,tsnoi,0,nx+1,0,ny+1)
      close(1)
      open (1,file=fext(fnmap,'tsw'),status='old')
      call readgrd(1,tsurw,0,nx+1,0,ny+1)
      close(1)
      open (1,file=fext(fnmap,'t2 '),status='old')
      call readgrd(1,t2noi,0,nx+1,0,ny+1)
      close(1)
      open (1,file=fext(fnmap,'wg '),status='old')
      call readgrd(1,wgnoi,0,nx+1,0,ny+1)
      close(1)
      open (1,file=fext(fnmap,'w2 '),status='old')
      call readgrd(1,w2noi,0,nx+1,0,ny+1)
      close(1)
      open (1,file=fext(fnmap,'wr '),status='old')
      call readgrd(1,wrnoi,0,nx+1,0,ny+1)
      close(1)
      open (1,file=fext(fnmap,'veg'),status='old')
      call readgrd(1,veg,0,nx+1,0,ny+1)
      close(1)
      open (1,file=fext(fnmap,'z0 '),status='old')
      call readgrd(1,z0,0,nx+1,0,ny+1)
      close(1)
      open (1,file=fext(fnmap,'z0w'),status='old')
      call readgrd(1,z0w,0,nx+1,0,ny+1)
      close(1)
      open (1,file=fext(fnmap,'lai'),status='old')
      call readgrd(1,xlai,0,nx+1,0,ny+1)
      close(1)
      open (1,file=fext(fnmap,'rsm'),status='old')
      call readgrd(1,rsmin,0,nx+1,0,ny+1)
      close(1)
      open (1,file=fext(fnmap,'alb'),status='old')
      call readgrd(1,alb,0,nx+1,0,ny+1)
      close(1)
      open (1,file=fext(fnmap,'icl'),status='old')
      call ireadgrd(1,iclay,0,nx+1,0,ny+1)
      close(1)
      open (1,file=fext(fnmap,'isa'),status='old')
      call ireadgrd(1,isand,0,nx+1,0,ny+1)
      close(1)
      open (1,file=fext(fnmap,'ive'),status='old')
      call ireadgrd(1,iveg,0,nx+1,0,ny+1)
      close(1)
      open (1,file=fext(fnmap,'irr'),status='old')
      call readgrd(1,zonirr,0,nx+1,0,ny+1)
      close(1)
      open (1,file=fext(fnmap,'lak'),status='old')
      call readgrd(1,xlake,0,nx+1,0,ny+1)
      close(1)
      open (1,file=fext(fnmap,'tsm'),status='old')
      call readgrd(1,tsmed,0,nx+1,0,ny+1)
      close(1)
	open (1,file=fext(fnmap,'icw'),status='old')
      call readgrd(1,watice,0,nx+1,0,ny+1)
      close(1)
	

      do ix=0,nx1
      do iy=0,ny1
         rra(ix,iy)=1.
         z0h(ix,iy)=z0(ix,iy)*z0hz0
         z0hw(ix,iy)=z0w(ix,iy)*z0hz0
         if(iclay(ix,iy).lt.0) then
            emis(ix,iy)=1.
         else
            emis(ix,iy)=0.98
         endif
         temmin(ix,iy)=tsnoi(ix,iy)
         temmax(ix,iy)=tsnoi(ix,iy)
         evapor(ix,iy)=0.
      enddo
      enddo


c
      zplc=0.
      zle=0.
      zletr=0.
      ruigc=0.
      ruitc=0.
      li=2
      lf=1
      return
      end



      subroutine sufhit(fnmap)
c----------------------------------------------------------------------
c surface height  hsuf
c-----------------------------------------------------------------------
c
      use alloc

      implicit real*8(a-h,o-z)
c      include 'nh3dpa10.inc'
      double precision x,y,dpi
      character*4 dsaa
      character*80 fnmap,fext
      external fext
c
      dpi=4.d0*datan(1.d0)
c
      do iy=0,ny1
      do ix=0,nx1
        hsuf(ix,iy)=0.
      enddo
      enddo
c
      if(iomtyp.eq.1) then
c
c bell-shaped 2d mountain:
c
         do iy=0,ny1
         do ix=0,nx1
           x=(ix-1.5)*dx
           y=(iy-1.5)*dy
           hsuf(ix,iy)=hmount/(1.d0+((x-xmount)/dble(hxwid))**2)
         enddo
         enddo
      elseif(iomtyp.eq.2) then
c
c bell-shaped 3d mountain:
c
         do iy=0,ny1
         do ix=0,nx1
           x=(ix-1.5)*dx
           y=(iy-1.5)*dy
           hsuf(ix,iy)=hmount/((1.d0+((x-xmount)/dble(hxwid))**2
     :       +((y-ymount)/dble(hywid))**2)**hexp)
         enddo
         enddo
      elseif(iomtyp.eq.3) then
c
c 3d mountain (thorsteinsson)
c
c         xywid=2.*hxwid
c         do 30 iy=0,ny1
c         do 30 ix=0,nx1
c         x=(ix-1.5)*dx-xmount
c         y=(iy-1.5)*dy-ymount
c         xy=dsqrt(x*x+y*y)
c         if(xy.lt.xywid) then
c            hsuf(ix,iy)=hmount*0.5d0*(1.d0+dcos(xy*dpi/xywid))
c         else
c            hsuf(ix,iy)=0.
c         endif
c30       continue

c thorpian sausage alps

       htemp=hmount+hexp
       do iy=0,ny1
       do ix=0,nx1
         x=(ix-1.5)*dx-xmount
         y=(iy-1.5)*dy-ymount
         hsuf(ix,iy)=dmax1((htemp/(1.e0+((dmax1((abs(x)
     :     -hywid),0.d0))**2
     :     +y*y)/(hxwid*hxwid))**2)-hexp,0.d0)
       enddo
       enddo

      elseif(iomtyp.eq.4) then
c
c reads orography from file ft22
c
         open(22,file=fext(fnmap,'top'),status='old')
         read(22,'(a4)') dsaa
         read(22,*) nnx,nny
         if(nnx.ne.nx+2 .or. nny.ne.ny+2) then
            write(0,*) 'nnx',nnx,nny
            write(0,*) 'nh3dFatalError=error in grd height file'
            stop
         endif
         read(22,*) xminutm,xmaxutm
         read(22,*) yminutm,ymaxutm
         read(22,*) hhhmin,hhhmax
c
c         do 40 ix=nx1,0,-1
c         read(22,*) (hsuf(ix,iy),iy=ny1,0,-1)
c40       continue
         do iy=0,ny+1
           read(22,*) (hsuf(ix,iy),ix=0,nx+1)
         enddo
         close(22)

         hmount=-1.e30
         do iy=0,ny+1
         do ix=0,nx+1
           hmount=max(hmount,hsuf(ix,iy))
         enddo
         enddo
      elseif(iomtyp.eq.5) then
c
c bell-shaped 2d valley:
c
         do iy=0,ny1
         do ix=0,nx1
           x=(ix-1.5)*dx
           y=(iy-1.5)*dy
           hsuf(ix,iy)=hmount+(-hmount/(1.d0+((x-xmount)
     :       /dble(hxwid))**2))
         enddo
         enddo
      elseif(iomtyp.eq.6) then
c
c bell-shaped 3d valley:
c
         do iy=0,ny1
         do ix=0,nx1
           x=(ix-1.5)*dx
           y=(iy-1.5)*dy
           hsuf(ix,iy)=hmount+(-hmount/((1.d0+((x-xmount)/dble(hxwid))
     :       **2+((y-ymount)/dble(hywid))**2)**hexp))
         enddo
         enddo

      elseif(iomtyp.eq.7) then
c
c sinusoidal sequence of mountains with buffer region
c amplitude: hmount // "wave lenght: hxwid
c buffer region  hxwid < x < nx1-hxwid
c
       domain=(nx1-1)*dx
       !buffer=domain*dx-hxwid
       space=(1+int(hexp))*(hxwid)

       !region=2*hxwid
       if (domain .le. space) then
          write(0,*) 'nh3dFatalError=error in sufhit'
          stop "Wrong values for hxwid or hexp"
       endif

       do iy=0,ny1
        do ix=0,nx1
          x=(ix-1.5)*dx
          y=(iy-1.5)*dy
          if (x .le. (hxwid-1000.) .or. x .ge. (space-1000.)) then
            hsuf(ix,iy)=0.
          else
c           hsuf(ix,iy)=(hmount+hmount*sin(2.*dpi*((2.*x/hxwid)-23./20.)
c    &                                                            ))/2.
            hsuf(ix,iy)=(hmount+hmount*sin(2.*dpi*((x/hxwid)-0.4)))/2.
          endif
        enddo
       enddo

      elseif(iomtyp.eq.8) then

c                                                           _______
c 1-D ramp of "upwind" inclination beta=hmount/hxwid  _____/
c

        domain=(nx-1)*dx
c       write(*,*) domain
        hh=hexp*xmount
        buffer=domain-hh
        declive=hmount/(hxwid)
        planalti=hh+hxwid  !inicio do planalto
        planalts=(buffer-hxwid)
c       write(*,*) 'planalti, planalts, hh'
c       write(*,*) planalti, planalts, hh
c       write(*,*) planalts-planalti-hh
        if ((domain-hh-hxwid).le.0.) then
          write(0,*) 'nh3dFatalError=error in sufhit 2'
           stop 'Erro na escolha de hexp, ou hxwid, ou hmount'
        endif

        do iy=0, ny1
         do ix=0, nx1
           x=(ix-1.5)*dx
           y=(iy-1.5)*dy
           if(x.lt.hh ) then
             hsuf(ix,iy)=0.
           else
             if(x.ge.(hh).and.x.lt.planalti) then
               hsuf(ix,iy)=declive*(x-hh)
             else
               if(x.ge.planalti.and.x.lt.planalts) then
                 hsuf(ix,iy)=hmount
               else
                 if(x.ge.planalts.and.x.lt.(planalts+hxwid)) then
                   hsuf(ix,iy)=hmount  !-declive*(x-(planalts+hxwid))
                 else
                   hsuf(ix,iy)=hmount  !0.
                 endif
               endif
             endif
           endif
         enddo
        enddo


      elseif(iomtyp.eq.9) then

c                                                      ___
c 1-D ramp of "downwind" inclination beta=hmount/hxwid    \_____
c

        domain=(nx-1)*dx
        hh=hexp*xmount
        buffer=domain-hh
        declive=hmount/(hxwid)
        planalti=hh+hxwid  !inicio do planalto
        planalts=(buffer-hxwid)
        if ((planalts-planalti-hh).le.0.) then
          write(0,*) 'nh3dFatalError=error in sufhit 3'
           stop 'Erro na escolha de hexp, ou hxwid, ou hmount'
        endif

        do iy=0, ny1
         do ix=0, nx1
           x=(ix-1.5)*dx
           y=(iy-1.5)*dy
           if(x.lt.hh ) then
             hsuf(ix,iy)=0.
           else
             if(x.ge.(hh).and.x.lt.planalti) then
               hsuf(ix,iy)=hmount
             else
               if(x.ge.planalti.and.x.lt.planalts) then
                 hsuf(ix,iy)=hmount
               else
                 if(x.ge.planalts.and.x.lt.(planalts+hxwid)) then
                   hsuf(ix,iy)=-declive*(x-(planalts+hxwid))
                 else
                   hsuf(ix,iy)=0.
                 endif
               endif
             endif
           endif
         enddo
        enddo


      elseif (iomtyp.eq.10) then

c                                                          ____
c symetrical 1-D ramp of inclination beta=hmount/hxwid  __/    \__
c

        domain=(nx-1)*dx
        hh=hexp*xmount
        buffer=domain-hh
        declive=hmount/(hxwid)
        planalti=hh+hxwid  !inicio do planalto
        planalts=(buffer-hxwid)
        if ((planalts-planalti-hh).le.0.) then
          write(0,*) 'nh3dFatalError=error in sufhit 4'
           stop 'Erro na escolha de hexp, ou hxwid, ou hmount'
        endif

        do iy=0, ny1
         do ix=0, nx1
           x=(ix-1.5)*dx
           y=(iy-1.5)*dy
           if(x.lt.hh ) then
             hsuf(ix,iy)=0.
           else
             if(x.ge.(hh).and.x.lt.planalti) then
               hsuf(ix,iy)=declive*(x-hh)
             else
               if(x.ge.planalti.and.x.lt.planalts) then
                 hsuf(ix,iy)=hmount
               else
                 if(x.ge.planalts.and.x.lt.(planalts+hxwid)) then
                   hsuf(ix,iy)=-declive*(x-(planalts+hxwid))
                 else
                   hsuf(ix,iy)=0.
                 endif
               endif
             endif
           endif
         enddo
        enddo

      elseif (iomtyp.eq.11) then
c three witches of agnesi same height
        do iy=0,ny1
         do ix=0,nx1
           x=(ix-1.5)*dx
           y=(iy-1.5)*dy
           hsuf(ix,iy)=hmount/((1.d0+((x-(xmount/2.d0))/dble(hxwid))**2
     :       +((y-ymount)/dble(hywid))**2)**hexp)         !
     :       +hmount/((1.d0+((x-xmount)/dble(hxwid))**2   !2nd mountain
     :       +((y-ymount)/dble(hywid))**2)**hexp)
     :       +hmount/((1.d0+((x-1.5d0*xmount)/dble(hxwid))**2   !3rd mountain
     :       +((y-ymount)/dble(hywid))**2)**hexp)
         enddo
        enddo

      elseif (iomtyp.eq.12) then
c three witches of agnesi central hill higher then the others
        do iy=0,ny1
         do ix=0,nx1
           x=(ix-1.5)*dx
           y=(iy-1.5)*dy
           hsuf(ix,iy)=(hmount/((1.d0+((x-(xmount/2.d0))/dble(hxwid))**2
     :       +((y-ymount)/dble(hywid))**2)**hexp))*0.5d0         !
     :       +hmount/((1.d0+((x-xmount)/dble(hxwid))**2   !2nd mountain
     :       +((y-ymount)/dble(hywid))**2)**hexp)
     :       +(hmount/((1.d0+((x-1.5d0*xmount)/dble(hxwid))**2   !3rd mountain
     :       +((y-ymount)/dble(hywid))**2)**hexp))*0.5d0
         enddo
        enddo

      endif


c
c Case with growing mountain
c
      if(ivarhs.ne.0) then
         do ix=0,nx1
            do iy=0,ny1
               hsufmax(ix,iy)=hsuf(ix,iy)
               hsuf(ix,iy)=0.
            enddo
         enddo
      endif

      do iy=0,ny1
      do ix=0,nx1
        phisuf(ix,iy)=g*hsuf(ix,iy)
      enddo
      enddo

c
      return
      end



      subroutine sufpp
c-----------------------------------------------------------------------
c     integration of continuity equation for surface pressure
c-----------------------------------------------------------------------
c
      use alloc
      use refstate
      implicit real*8(a-h,o-z)
c      include 'nh3dpa10.inc'
c
      logical sim,nao
      parameter (sim=.true.,nao=.false.)
c
c      dimension duul(ny1),duur(ny1),dvvl(nx1),dvvr(nx1)
c
      if( nstep.le.1) then
         do iy=1,ny1
         do ix=1,nx1
           pp(ix,iy,3)=pp(ix,iy,1)+dpp(ix,iy)*dtl
         enddo
         enddo
c
         call extrpp(nx2,ny2,pp(0,0,3),0,nx2,0,ny2)
      endif
c-----------------------------------------------------------------------
c     compute surface pressure tendency  (dpp)
c-----------------------------------------------------------------------
!     call zero(dpp,nxy)

      dpp=0.
      do is=1,ns-1
      do iy=1,ny1
      do ix=1,nx1
        dpp(ix,iy)=dpp(ix,iy)-(u(ix,iy,is,3)*pp10(ix,iy,3)
     :    -u(ix-1,iy,is,3)*pp10(ix-1,iy,3))*ds1dx(is)
     :    -(v(ix,iy,is,3)*pp01(ix,iy,3)
     :    -v(ix,iy-1,is,3)*pp01(ix,iy-1,3))*ds1dy(is)
      enddo
      enddo
      enddo
c
c         tk=0.
c         call wri2ar(dpp,1,nx1,1,ny1,'dpp1     ',tk,nx,ny)
c
      do iy=1,ny1
      do ix=1,nx1
        pp(ix,iy,4)=pp(ix,iy,2)+dpp(ix,iy)*dt2
      enddo
      enddo
c
      if(ioppuv.ne.1) then
         if(iobppx.eq.0) then
            do iy=2,ny
              ppbx(1,iy)=pp(2,iy,4)
              ppbx(2,iy)=pp(nx,iy,4)
            enddo
            ppcc(1,1)=ppby(2,1)
            ppcc(2,1)=ppby(nx,1)
            ppcc(2,2)=ppby(nx,2)
            ppcc(1,2)=ppby(2,2)
         endif
c
         if(iobppy.eq.0) then
            do ix=2,nx
              ppby(ix,1)=pp(ix,2,4)
              ppby(ix,2)=pp(ix,ny,4)
            enddo
            ppcc(1,1)=ppbx(1,2)
            ppcc(2,1)=ppbx(2,2)
            ppcc(2,2)=ppbx(2,ny)
            ppcc(1,2)=ppbx(1,ny)
         endif
c
c calculate correction to secondary boundary points for uu and vv
c
         do ix=2,nx
           dvvl(ix)=(ppby(ix,1)-pp(ix,1,4))*dydt2
           dvvr(ix)=-(ppby(ix,2)-pp(ix,ny1,4))*dydt2
         enddo
c
         do iy=2,ny
           duul(iy)=(ppbx(1,iy)-pp(1,iy,4))*dxdt2
           duur(iy)=-(ppbx(2,iy)-pp(nx1,iy,4))*dxdt2
         enddo
c
         sduul=zero0
         sduur=zero0
         sdvvl=zero0
         sdvvr=zero0
         do ix=2,nx
           sdvvl=sdvvl+abs(dvvl(ix))
           sdvvr=sdvvr+abs(dvvr(ix))
         enddo
         do iy=2,ny
           sduul=sduul+abs(duul(iy))
           sduur=sduur+abs(duur(iy))
         enddo
c
         aduu11=sduul/(sduul+sdvvl)
         advv11=1.-aduu11
         aduu12=sduul/(sduul+sdvvr)
         advv12=1.-aduu12
         aduu21=sduur/(sduur+sdvvl)
         advv21=1.-aduu21
         aduu22=sduur/(sduur+sdvvr)
         advv22=1.-aduu22
c
         dpp11=(ppcc(1,1)-pp(1,1,4))/dt2
         duul(1)=dpp11*dx*aduu11
         dvvl(1)=dpp11*dy*advv11
         dpp21=(ppcc(2,1)-pp(nx1,1,4))/dt2
         duur(1)=-dpp21*dx*aduu21
         dvvl(nx1)=dpp21*dy*advv21
         dpp22=(ppcc(2,2)-pp(nx1,ny1,4))/dt2
         duur(ny1)=-dpp22*dx*aduu22
         dvvr(nx1)=-dpp22*dy*advv22
         dpp12=(ppcc(1,2)-pp(1,ny1,4))/dt2
         duul(ny1)=dpp12*dx*aduu12
         dvvr(1)=-dpp12*dy*advv12
c
         if(duvout.and.prt) then
           write(nchn1,*) ' aduu11=',aduu11,' aduu21=',aduu21
     :       ,' aduu22=',aduu22,' aduu12=',aduu12
     :       ,' advv11=',advv11,' advv21=',advv21,' advv22=',advv22
     :       ,' advv12=',advv12
           write(nchn1,*) 'duul'
           write(nchn1,*) (duul(iy),iy=1,ny1)
           write(nchn1,*) 'duur'
           write(nchn1,*) (duur(iy),iy=1,ny1)
           write(nchn1,*) 'dvvl'
           write(nchn1,*) (dvvl(ix),ix=1,nx1)
           write(nchn1,*) 'dvvr'
           write(nchn1,*) (dvvr(ix),ix=1,nx1)
         endif
c
         if(iobvvy.eq.1) then
            do ix=1,nx1
            do is=0,ns
              v(ix,0,is,3)=v(ix,0,is,3)+dvvl(ix)/pp01(ix,0,3)
              v(ix,ny1,is,3)=v(ix,ny1,is,3)+dvvr(ix)/pp01(ix,ny1,3)
            enddo
            enddo
         endif
c
         if(iobuux.eq.1) then
            do iy=1,ny1
            do is=0,ns
              u(0,iy,is,3)=u(0,iy,is,3)+duul(iy)/pp10(0,iy,3)
              u(nx1,iy,is,3)=u(nx1,iy,is,3)+duur(iy)/pp10(nx1,iy,3)
            enddo
            enddo
         endif
c
         do iy=2,ny
           pp(1,iy,4)=ppbx(1,iy)
           pp(nx1,iy,4)=ppbx(2,iy)
         enddo
         do ix=2,nx
           pp(ix,1,4)=ppby(ix,1)
           pp(ix,ny1,4)=ppby(ix,2)
         enddo
         pp(1,1,4)=ppcc(1,1)
         pp(nx1,1,4)=ppcc(2,1)
         pp(nx1,ny1,4)=ppcc(2,2)
         pp(1,ny1,4)=ppcc(1,2)
c
         do ix=1,nx1
           dpp(ix,1)=(pp(ix,1,4)-pp(ix,1,2))/dt2
           dpp(ix,ny1)=(pp(ix,ny1,4)-pp(ix,ny1,2))/dt2
         enddo
c
         do iy=1,ny1
           dpp(1,iy)=(pp(1,iy,4)-pp(1,iy,2))/dt2
           dpp(nx1,iy)=(pp(nx1,iy,4)-pp(nx1,iy,2))/dt2
         enddo
c
      endif

c lateral sponges for dpp Use only if reference state is being forced
c varying in time (as in that case there is no lateral flux correction)

      if(ifuprefs.ne.0 .and. forcedpp) then
        if(nxsponge2.gt.0 .or. nxsponge1.gt.0) then
        if(verbose.ge.3) write(*,*) 'sufpp dpp:sponge'
          do iy=0,ny1
            do ix=0,nxsponge2
              dpp(ix,iy)=dpp(ix,iy)-(dpp(ix,iy)-dppref(ix,iy))
     :          *dtl/tauspo(ix)
            enddo
            do ix=nx1-nxsponge1,nx1
              dpp(ix,iy)=dpp(ix,iy)-(dpp(ix,iy)-dppref(ix,iy))
     :          *dtl/tauspo(ix)
            enddo
          enddo
        endif
c
c Wind   - Lateral forcing column - ymin and y max (lateral sponge)
c
        if(nysponge2.gt.0 .or. nysponge1.gt.0) then
          do ix=0,nx1
            do iy=0,nysponge2
              dpp(ix,iy)=dpp(ix,iy)-(dpp(ix,iy)-dppref(ix,iy))
     :          *dtl/tauspo(ix)
            enddo
            do iy=ny1-nysponge1,ny1
              dpp(ix,iy)=dpp(ix,iy)-(dpp(ix,iy)-dppref(ix,iy))
     :          *dtl/tauspo(ix)
            enddo
          enddo
        endif
      endif
c
      if(phsout.and.prt) then
         tk=0.
         call wri2ar(dpp,1,nx1,1,ny1,'dpp      ',tk,nx,ny)
      endif
c
      call ppbc(nx,ny,pp(0,0,4),ppbx,ppbx2,ppby,ppby2,ppcc,ppcc2
     :   ,iobppx,iobppy,ppbout,prt,nchn1)
c
      call aselin(pp(0,0,4),pp(0,0,3),pp(0,0,2),0,nx2,0,ny2,0,0,filt)
c

      if(ifuprefs.ne.0 .and. forcepsuf) then
        if(nxsponge2.gt.0 .or. nxsponge1.gt.0) then
        if(verbose.ge.3) write(*,*) 'sufpp dpp:sponge'
          do iy=0,ny1
            do ix=0,nxsponge2
              psuf(ix,iy)=psuf(ix,iy)-(psuf(ix,iy)-psufref(ix,iy))
     :          *dtl/tauspo(ix)
            enddo
            do ix=nx1-nxsponge1,nx1
              psuf(ix,iy)=psuf(ix,iy)-(psuf(ix,iy)-psufref(ix,iy))
     :          *dtl/tauspo(ix)
            enddo
          enddo
        endif
c
c Wind   - Lateral forcing column - ymin and y max (lateral sponge)
c
        if(nysponge2.gt.0 .or. nysponge1.gt.0) then
          do ix=0,nx1
            do iy=0,nysponge2
              psuf(ix,iy)=psuf(ix,iy)-(psuf(ix,iy)-psufref(ix,iy))
     :          *dtl/tauspo(ix)
            enddo
            do iy=ny1-nysponge1,ny1
              psuf(ix,iy)=psuf(ix,iy)-(psuf(ix,iy)-psufref(ix,iy))
     :          *dtl/tauspo(ix)
            enddo
          enddo
        endif
        do iy=0,ny1
        do ix=0,nx1
          pp(ix,iy,3)=psuf(ix,iy)-ptop
        enddo
        enddo
      endif


      do iy=0,ny1
      do ix=0,nx1
        psuf(ix,iy)=pp(ix,iy,3)+ptop
        pp10(ix,iy,3)=0.5*(pp(ix,iy,3)+pp(ix+1,iy,3))
        pp10(ix,iy,4)=0.5*(pp(ix,iy,4)+pp(ix+1,iy,4))
        dpdx10(ix,iy,3)=(pp(ix+1,iy,3)-pp(ix,iy,3))/dx
        pp01(ix,iy,3)=0.5*(pp(ix,iy,3)+pp(ix,iy+1,3))
        pp01(ix,iy,4)=0.5*(pp(ix,iy,4)+pp(ix,iy+1,4))
        dpdy01(ix,iy,3)=(pp(ix,iy+1,3)-pp(ix,iy,3))/dy
      enddo
      enddo
c
      do iy=1,ny1
      do ix=1,nx1
        ppdx(ix,iy,3)=(pp(ix+1,iy,3)-pp(ix-1,iy,3))/(dx2*pp(ix,iy,3))
        ppdx(ix,iy,4)=(pp(ix+1,iy,4)-pp(ix-1,iy,4))/(dx2*pp(ix,iy,4))
        ppdy(ix,iy,3)=(pp(ix,iy+1,3)-pp(ix,iy-1,3))/(dy2*pp(ix,iy,3))
        ppdy(ix,iy,4)=(pp(ix,iy+1,4)-pp(ix,iy-1,4))/(dy2*pp(ix,iy,4))
      enddo
      enddo
c
      call radbch(ppcc,ppbx,ppby,pp(0,0,4),pp(0,0,3),pp(0,0,2)
     :   ,nx2,ny2,0,0,0,1,nx1,1,ny1)
      call radbch(ppcc2,ppbx2,ppby2,pp(0,0,4),pp(0,0,3),pp(0,0,2)
     :   ,nx2,ny2,0,0,0,0,nx2,0,ny2)
c
      return
      end



      subroutine surhyd
c
      use alloc
      implicit real*8(a-h,o-z)
c      include 'nh3dpa10.inc'
c
      do iy=1,ny1
      do ix=1,nx1
        phisuf(ix,iy)=g*hsuf(ix,iy)-r05*(tsuf(ix,iy)+tsufi(ix,iy))
     :    *dlog(psuf(ix,iy)/psufi(ix,iy))
      enddo
      enddo
c
      do iy=1,ny1
        phisuf(0,iy)=phisuf(1,iy)
      enddo
c
      do ix=0,nx1
        phisuf(ix,0)=phisuf(ix,1)
      enddo
c
      if(phsout.and.prt) then
         tk=0.
         call wri2ar(phisuf,0,nx1,0,ny1,'phisuf  ',tk,nx,ny)
      endif
      return
      end




      subroutine surphi
c
c the following code performs the solution of the elliptic problem
c for the surface boundary condition for phis, using dirichlet non-
c -homogeneous boundary conditions, and integrates the hydrostatic
c equation to get phis everywhere.
c
c lateral boundary conditions for phis at the surface are obtained
c from the solution of 1 d similar equations. the value at the four
c corners is imposed (using hydrostatic relation)
c
c-----------------------------------------------------------------------
c
c dirichlet lateral boundary condition for phisuf.
c-----------------------------------------------------------------------
c
      use alloc
      implicit real*8(a-h,o-z)
c      include 'nh3dpa10.inc'
c
c      dimension phi00(0:nx1,0:ny1)
c      dimension phibx1(0:ny1),phiby1(0:nx1)
c      dimension phibx2(0:ny1),phiby2(0:nx1)
c      dimension a1(nx-1),a2(nx-1),a3(nx-1),c(nx-1),d(nx-1)
c
c reevaluate variables of reference state from new pp(i,3)
c
c
c lateral b.c. at the surface
c
      do ix=2,nx
        phiby1(ix)=-r05*(
     :    (dpdx10(ix,1,3)*(tsuf(ix+1,1)/psuf(ix+1,1)
     :    +tsuf(ix,1)/psuf(ix,1))
     :    -dpdx10(ix-1,1,3)*(tsuf(ix,1)/psuf(ix,1)
     :    +tsuf(ix-1,1)/psuf(ix-1,1)))/dx)
        phiby2(ix)=-r05*(
     :    (dpdx10(ix,ny1,3)*(tsuf(ix+1,ny1)/psuf(ix+1,ny1)
     :    +tsuf(ix,ny1)/psuf(ix,ny1))
     :    -dpdx10(ix-1,ny1,3)*(tsuf(ix,ny1)/psuf(ix,ny1)
     :    +tsuf(ix-1,ny1)/psuf(ix-1,ny1)))/dx)
      enddo
c
      do iy=2,ny
        phibx1(iy)=-r05*(
     :    (dpdy01(1,iy,3)*(tsuf(1,iy+1)/psuf(1,iy+1)
     :    +tsuf(1,iy)/psuf(1,iy))
     :    -dpdy01(1,iy-1,3)*(tsuf(1,iy)/psuf(1,iy)
     :    +tsuf(1,iy-1)/psuf(1,iy-1)))/dy)
        phibx2(iy)=-r05*(
     :    (dpdy01(nx1,iy,3)*(tsuf(nx1,iy+1)/psuf(nx1,iy+1)
     :    +tsuf(nx1,iy)/psuf(nx1,iy))
     :    -dpdy01(nx1,iy-1,3)*(tsuf(nx1,iy)/psuf(nx1,iy)
     :    +tsuf(nx1,iy-1)/psuf(nx1,iy-1)))/dy)
      enddo
c
c uses hydrostatic relation to calculate phisuf at the four corner points:
c
      phisuf(1,1)=phisui(1,1)-r05*(tsuf(1,1)+tsufi(1,1))
     :   *dlog(psuf(1,1)/psufi(1,1))
      phisuf(1,ny1)=phisui(1,ny1)-r05*(tsuf(1,ny1)+tsufi(1,ny1))
     :   *dlog(psuf(1,ny1)/psufi(1,ny1))
      phisuf(nx1,1)=phisui(nx1,1)-r05*(tsuf(nx1,1)+tsufi(nx1,1))
     :   *dlog(psuf(nx1,1)/psufi(nx1,1))
      phisuf(nx1,ny1)=phisui(nx1,ny1)-r05*(tsuf(nx1,ny1)+tsufi(nx1,ny1))
     :   *dlog(psuf(nx1,ny1)/psufi(nx1,ny1))
c
c      phisuf(1,1)=phisuf(1,1)-r*tsuf(1,1)
c     :   *dlog(psuf(1,1)/psuf(1,1,2))
c      phisuf(1,ny1)=phisuf(1,ny1)-r*tsuf(1,ny1)
c     :   *dlog(psuf(1,ny1)/psuf(1,ny1,2))
c      phisuf(nx1,1)=phisuf(nx1,1)-r*tsuf(nx1,1)
c     :   *dlog(psuf(nx1,1)/psuf(nx1,1,2))
c      phisuf(nx1,ny1)=phisuf(nx1,ny1)-r*tsuf(nx1,ny1)
c     :   *dlog(psuf(nx1,ny1)/psuf(nx1,ny1,2))
c
      if(phsout .and. prt) then
        write(nchn1,3000) phisuf(1,1),phisuf(nx1,1),phisuf(1,ny1)
     :     ,phisuf(nx1,ny1)
        write(nchn1,3001) (phibx1(i),i=2,ny)
        write(nchn1,3001) (phibx2(i),i=2,ny)
        write(nchn1,3001) (phiby1(i),i=2,nx)
        write(nchn1,3001) (phiby2(i),i=2,nx)
3000    format(' phisuf(1,1)=',e22.15,' phisuf(nx1,1)=',e22.15
     :     /,' phisuf(1,ny1)=',e22.15,' phisuf(nx1,ny1)=',e22.15)
3001    format(1x,3e22.15)
      endif
c
c solve 1d poisson equation
c
      call pos1(phibx1,ny,dyy,phisuf(1,1),phisuf(1,ny1)
     :  ,a1sphi,a2sphi,a3sphi,csphi,dsphi)
      call pos1(phibx2,ny,dyy,phisuf(nx1,1),phisuf(nx1,ny1)
     :  ,a1sphi,a2sphi,a3sphi,csphi,dsphi)
      call pos1(phiby1,nx,dxx,phisuf(1,1),phisuf(nx1,1)
     :  ,a1sphi,a2sphi,a3sphi,csphi,dsphi)
      call pos1(phiby2,nx,dxx,phisuf(1,ny1),phisuf(nx1,ny1)
     :  ,a1sphi,a2sphi,a3sphi,csphi,dsphi)
c
c take solution:
c
      do iy=2,ny
        phisuf(1,iy)=phibx1(iy)
        phisuf(nx1,iy)=phibx2(iy)
      enddo
      do ix=2,nx
        phisuf(ix,1)=phiby1(ix)
        phisuf(ix,ny1)=phiby2(ix)
      enddo
c
c calculate bottom b.c. for phis from new surface pressure.
c solves poisson equation for surface phis, with dirichlet lateral
c b.c. , i.e., fixed phisuf at ix=1,nx1;iy=1,ny1 (otherwise it can
c be imposed). the method consists in transforming the problem in one
c with homogeneous b.c. in a way similar to williams (jfm...).
c
c due to the fact that the poisson solution is made by use of a fourier
c transform in the x direction followed by explicit solution in the y
c direction , the routine takes care of the staggering of the grid
c only the x direction (a sttagered transform is performed). the boundary
c condition in x is then applied at x=1.5dx;nx+.5dx, where phi00 is 0.
c
      do iy=2,ny
      do ix=2,nx
        phi00s(ix,iy)=-r05*(
     :    (dpdx10(ix,iy,3)*(tsuf(ix+1,iy)/psuf(ix+1,iy)+tsuf(ix,iy)
     :    /psuf(ix,iy))
     :    -dpdx10(ix-1,iy,3)*(tsuf(ix,iy)/psuf(ix,iy)+tsuf(ix-1,iy)
     :    /psuf(ix-1,iy)))/dx
     :    +(dpdy01(ix,iy,3)*(tsuf(ix,iy+1)/psuf(ix,iy+1)+tsuf(ix,iy)
     :    /psuf(ix,iy))
     :    -dpdy01(ix,iy-1,3)*(tsuf(ix,iy)/psuf(ix,iy)+tsuf(ix,iy-1)
     :    /psuf(ix,iy-1)))/dy)
      enddo
      enddo
c
c modify forcing to impose homogeneous dirichlet b.c. (only in x):
c
      do iy=2,ny
        phi00s(2,iy)=phi00s(2,iy)-phisuf(1,iy)/dxx
        phi00s(nx,iy)=phi00s(nx,iy)-phisuf(nx1,iy)/dxx
      enddo
c
      if(iobphy.eq.1) then
         do ix=2,nx
           phi00s(ix,2)=phi00s(ix,2)-phisuf(ix,1)/dyy
           phi00s(ix,ny)=phi00s(ix,ny)-phisuf(ix,ny1)/dyy
         enddo
      endif
c
      if(phsout.and.prt) then
         tk=0.
         call wri2ar(phi00s,2,nx,2,ny,'phi00   ',tk,nx,ny)
      endif
c
c solve homogenous poisson equation
c
      call pos2(phi00s,iobphy)
c
c take solution (b.c. are implicit):
c
      do iy=2,ny
      do ix=2,nx
        phisuf(ix,iy)=phi00s(ix,iy)
      enddo
      enddo
c
      do iy=1,ny1
        phisuf(0,iy)=phisuf(1,iy)
      enddo
c
      do ix=0,nx1
        phisuf(ix,0)=phisuf(ix,1)
      enddo
c
      if(phsout.and.prt) then
         tk=0.
         call wri2ar(phisuf,0,nx1,0,ny1,'phisuf  ',tk,nx,ny)
      endif
c
      return
      end



      subroutine surphn
c
c ***
c
c
c the following code performs the solution of the elliptic problem
c for the surface boundary condition for phis, using neuman non-
c -homogeneous boundary conditions, and integrates the hydrostatic
c equation to get phis everywhere.
c
c
c ***
c
c-----------------------------------------------------------------------
c re-evaluate reference state variables and compute phi at surface
c with pts prespecified
c-----------------------------------------------------------------------
c
      use alloc
      implicit real*8(a-h,o-z)
c      include 'nh3dpa10.inc'
c      common/bphi/ dphsdx(2,0:ny1),dphsdy(0:nx1,2)
c      dimension hk11(0:nx1,0:ny1),hk12(0:nx1,0:ny1)
c
c-----------------------------------------------------------------------
c     reevaluate variables of reference state from new pp(i,3)
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c calculate bottom b.c. for phis from new surface pressure.
c solves poisson equation for surface phis, with neuman b.c.
c following williams (jfm,1969). both the diagnostic equation for
c phisuf and its lateral boundary conditions are imposed by pp.
c-----------------------------------------------------------------------
      do iy=2,ny
      do ix=2,nx
        phisuf(ix,iy)=-r05*(
     :    (dpdx10(ix,iy,3)*(tsuf(ix+1,iy)/psuf(ix+1,iy)+tsuf(ix,iy)
     :    /psuf(ix,iy))
     :    -dpdx10(ix-1,iy,3)*(tsuf(ix,iy)/psuf(ix,iy)+tsuf(ix-1,iy)
     :    /psuf(ix-1,iy)))/dx
     :    +(dpdy01(ix,iy,3)*(tsuf(ix,iy+1)/psuf(ix,iy+1)+tsuf(ix,iy)
     :    /psuf(ix,iy))
     :    -dpdy01(ix,iy-1,3)*(tsuf(ix,iy)/psuf(ix,iy)+tsuf(ix,iy-1)
     :    /psuf(ix,iy-1)))/dy)
      enddo
      enddo
c
c lateral boundary conditions:
c
      do iy=2,ny
        dphsdx(1,iy)=-r05*dpdx10(1,iy,3)*(tsuf(1,iy)
     :    /psuf(1,iy)+tsuf(2,iy)/psuf(2,iy))
        dphsdx(2,iy)=-r05*dpdx10(nx,iy,3)*(tsuf(nx,iy)
     :    /psuf(nx,iy)+tsuf(nx1,iy)/psuf(nx1,iy))
      enddo
      dphsdx(1,1)=dphsdx(1,2)
      dphsdx(1,ny1)=dphsdx(1,ny)
      dphsdx(2,1)=dphsdx(2,2)
      dphsdx(2,ny1)=dphsdx(2,ny)
      dphsdx(1,0)=dphsdx(1,1)
c
      do ix=2,nx
        dphsdy(ix,1)=-r05*dpdy01(ix,1,3)*(tsuf(ix,1)
     :    /psuf(ix,1)+tsuf(ix,2)/psuf(ix,2))
        dphsdy(ix,2)=-r05*dpdy01(ix,ny,3)*(tsuf(ix,ny)
     :    /psuf(ix,ny)+tsuf(ix,ny1)/psuf(ix,ny1))
      enddo
      dphsdy(1,1)=dphsdy(2,1)
      dphsdy(nx1,1)=dphsdy(nx,1)
      dphsdy(1,2)=dphsdy(2,2)
      dphsdy(nx1,2)=dphsdy(nx,2)
      dphsdy(0,1)=dphsdy(1,1)
c
      do iy=2,ny
        phisuf(2,iy)=phisuf(2,iy)+dphsdx(1,iy)/dx
        phisuf(nx,iy)=phisuf(nx,iy)-dphsdx(2,iy)/dx
      enddo
c
      do ix=2,nx      
        phisuf(ix,2)=phisuf(ix,2)+dphsdy(ix,1)/dy
        phisuf(ix,ny)=phisuf(ix,ny)-dphsdy(ix,2)/dy
      enddo
c
      i0=nx-1
      j0=2
c      write(6,*) 'surphn:',i0,j0,psuf(i0,j0),psufi(i0,j0)
c
      phsfco=g*hsuf(i0,j0)-r05*(tsuf(i0,j0)+tsufi(i0,j0))
     :   *dlog(psuf(i0,j0)/psufi(i0,j0))
c
      if(phsout.and.prt) then
        tk=0.
        call wrigar(dphsdx,1,2,1,ny1,0,0,1,2,1,ny1,0,0,'dphsdx   ',tk,1)
        call wrigar(dphsdy,1,nx1,1,2,0,0,1,nx1,1,2,0,0,'dphsdy   ',tk,1)
        write(nchn1,*) 'phsfco=',phsfco
        call wri2ar(phisuf,2,nx,2,ny,'phisuf00',tk,nx,ny)
      endif
c
c solve homogenous poisson equation
c
      allocate(hk11(0:nx1,0:ny1),hk12(0:nx1,0:ny1))

      if(verbose.ge.2) write(*,*) 'Surphn: Call pois2n'
      if(verbose.ge.3) write(*,'(10i4)') (ifax(k19),k19=1,10)
      call pois2n(phisuf,phsfco,i0,j0,hk11,hk12,dxx,dyy,eigxy
     :  ,trigsx,ifax,sipcox,simcox,nx
     :  ,trigsy,ifay,sipcoy,simcoy,ny)
      if(verbose.ge.3) write(*,*) 'pois2 done'

      deallocate(hk11,hk12)

c
c impose b.c.:
c
      do ix=2,nx
        phisuf(ix,1)=phisuf(ix,2)-dphsdy(ix,1)*dy
        phisuf(ix,ny1)=phisuf(ix,ny)+dphsdy(ix,2)*dy
      enddo
c
      do iy=1,ny1
        phisuf(1,iy)=phisuf(2,iy)-dphsdx(1,iy)*dx
        phisuf(nx1,iy)=phisuf(nx,iy)+dphsdx(2,iy)*dx
      enddo
c
      do iy=1,ny1
        phisuf(0,iy)=phisuf(1,iy)
      enddo
c
      do ix=0,nx1
        phisuf(ix,0)=phisuf(ix,1)
      enddo
c
      if(phsout.and.prt) then
        tk=0.
        call wri2ar(phisuf,0,nx1,0,ny1,'phisuf  ',tk,nx,ny)
      endif
c
      return
      end



      subroutine surpsd
c-----------------------------------------------------------------------
c initial psuf (case fcor.ne.0)
c boundary conditions dirichlet (psuf on four corners imposed)
c-----------------------------------------------------------------------
c
      use alloc
      implicit real*8(a-h,o-z)
c      include 'nh3dpa10.inc'
c
      parameter(nttt=5)

c      dimension ff(0:nx1,0:ny1),tbot(0:nx1,0:ny1)
c      dimension dlnpdx(0:nx1,0:ny1),dlnpdy(0:nx1,0:ny1)
c      dimension dlntdx(0:nx1,0:ny1),dlntdy(0:nx1,0:ny1)
c
      allocate(ff(0:nx1,0:ny1),tbot(0:nx1,0:ny1))
      allocate(dlnpdx(0:nx1,0:ny1),dlnpdy(0:nx1,0:ny1))
      allocate(dlntdx(0:nx1,0:ny1),dlntdy(0:nx1,0:ny1))

      ixref=nx1/2 !2
      iyref=ny1/2
      pref=psuf(ixref,iyref)
c
      do 9000 ittt=1,nttt

      do iy=1,ny1
      do ix=1,nx1
        tbot(ix,iy)=tsuf(ix,iy)
c     :   +(pt(ix,iy,ns,3)+pt(ix,iy,ns-1,3))
c     :   *0.5*(psuf(ix,iy)/p00)**akapa
      enddo
      enddo
c
c integration (euler)
c
      iy=1
      do ix=2,nx1
        psuf(ix,iy)=psuf(ix-1,iy)+(psuf(ix,iy)+psuf(ix-1,iy))
     :    /(r*(tbot(ix,iy)+tbot(ix-1,iy)))
     :    *(fcor*vg-g*(hsuf(ix,iy)-hsuf(ix-1,iy))/dx)*dx
      enddo
c
      do iy=2,ny1
        ix=1
        psuf(ix,iy)=psuf(ix,iy-1)+(psuf(ix,iy)+psuf(ix,iy-1))
     :    /(r*(tbot(ix,iy)+tbot(ix,iy-1)))
     :    *(-fcor*ug-g*(hsuf(ix,iy)-hsuf(ix,iy-1))/dy)*dy
        ix=nx1
        psuf(ix,iy)=psuf(ix,iy-1)+(psuf(ix,iy)+psuf(ix,iy-1))
     :    /(r*(tbot(ix,iy)+tbot(ix,iy-1)))
     :    *(-fcor*ug-g*(hsuf(ix,iy)-hsuf(ix,iy-1))/dy)*dy
      enddo
c
      iy=ny1
      do ix=2,nx1
        psuf(ix,iy)=psuf(ix-1,iy)+(psuf(ix,iy)+psuf(ix-1,iy))
     :    /(r*(tbot(ix,iy)+tbot(ix-1,iy)))
     :    *(fcor*vg-g*(hsuf(ix,iy)-hsuf(ix-1,iy))/dx)*dx
      enddo
9000  continue
c
      delprf=pref-psuf(ixref,iyref)
      do iy=1,ny1
      do ix=1,nx1
        psuf(ix,iy)=psuf(ix,iy)+delprf
      enddo
      enddo
c
      if(phsout .and. prt) then
         write(nchn1,3000) psuf(1,1),psuf(nx1,1),psuf(1,ny1)
     :      ,psuf(nx1,ny1)
3000     format(' psuf(1,1)=',e14.7,' psuf(nx1,1)=',e14.7
     :      ,' psuf(1,ny1)=',e14.7,' psuf(nx1,ny1)=',e14.7)
      endif
c
      do iy=2,ny
      do ix=2,nx
        dlnpdx(ix,iy)=(fcor*vg-g*(hsuf(ix+1,iy)-hsuf(ix-1,iy))/dx2)
     :    /(r*tbot(ix,iy))
        dlnpdy(ix,iy)=(-fcor*ug-g*(hsuf(ix,iy+1)-hsuf(ix,iy-1))/dy2)
     :    /(r*tbot(ix,iy))
        dlntdx(ix,iy)=(tbot(ix+1,iy)-tbot(ix-1,iy))/(dx2*tbot(ix,iy))
        dlntdy(ix,iy)=(tbot(ix,iy+1)-tbot(ix,iy-1))/(dy2*tbot(ix,iy))
      enddo
      enddo
c
      do iy=2,ny
      do ix=2,nx
        ff(ix,iy)=
     :    psuf(ix,iy)*(dlnpdx(ix,iy)*(dlnpdx(ix,iy)-dlntdx(ix,iy))
     :    +dlnpdy(ix,iy)*(dlnpdy(ix,iy)-dlntdy(ix,iy))
     :    -g/(r*tbot(ix,iy))
     :    *((hsuf(ix+1,iy)+hsuf(ix-1,iy)-2.*hsuf(ix,iy))/dxx
     :    +(hsuf(ix,iy+1)+hsuf(ix,iy-1)-2.*hsuf(ix,iy))/dyy))
      enddo
      enddo
c
c lateral boundary conditions:
c
      do iy=2,ny
        ff(2,iy)=ff(2,iy)-psuf(1,iy)/dxx
        ff(nx,iy)=ff(nx,iy)-psuf(nx1,iy)/dxx
      enddo
c
      if(iobphy.eq.1) then
        do ix=2,nx
          ff(ix,2)=ff(ix,2)-psuf(ix,1)/dyy
          ff(ix,ny)=ff(ix,ny)-psuf(ix,ny1)/dyy
        enddo
      endif
c
c solve homogenous poisson equation
c
      call pos2(ff,iobphy)
c
c take solution (b.c. are implicit):
c
      do iy=2,ny
      do ix=2,nx
        psuf(ix,iy)=ff(ix,iy)
      enddo
      enddo
c
      do iy=1,ny1
        psuf(0,iy)=psuf(1,iy)
      enddo
c
      do ix=0,nx1
        psuf(ix,0)=psuf(ix,1)
      enddo
c
      do iy=1,ny
      do ix=1,nx1
        ff(ix,iy)=-1./fcor*((r*(tbot(ix,iy)+tbot(ix,iy+1))
     :    /(psuf(ix,iy)+psuf(ix,iy+1)))
     :    *(psuf(ix,iy+1)-psuf(ix,iy))/dy
     :    +g*(hsuf(ix,iy+1)-hsuf(ix,iy))/dy)
      enddo
      enddo

      if(phsout.and.prt) then
         tk=1.e5
         call wri2ar(psuf,0,nx1,0,ny1,'psuf-p5  ',tk,nx,ny)
         tk=0.
         call wri2ar(ff,1,nx1,1,ny,'ug      ',tk,nx,ny)
      endif

      deallocate(ff,tbot,dlnpdx,dlnpdy,dlntdx,dlntdy)
c
      return
      end



      subroutine surpsu(ixref,iyref,pref)
c-----------------------------------------------------------------------
c initial psuf (case fcor.ne.0)
c-----------------------------------------------------------------------
c
      use alloc
      implicit real*8(a-h,o-z)
c
      allocate(dpsudx(2,ny1),dpsudy(nx1,2))
      allocate(ff(0:nx1,0:ny1),tbot(0:nx1,0:ny1))
      allocate(dlnpdx(0:nx1,0:ny1),dlnpdy(0:nx1,0:ny1))
      allocate(dlntdx(0:nx1,0:ny1),dlntdy(0:nx1,0:ny1))
      allocate(hk21(0:nx1,0:ny1,2),hk22(0:nx1,0:ny1,2))
c
      do iy=1,ny1
      do ix=1,nx1
        tbot(ix,iy)=tsuf(ix,iy)+(pt(ix,iy,ns,3)+pt(ix,iy,ns-1,3))
     :    *0.5*(psuf(ix,iy)/p00)**akapa
      enddo
      enddo
      do iy=2,ny
      do ix=2,nx
        dlnpdx(ix,iy)=(fcor*vg-g*(hsuf(ix+1,iy)-hsuf(ix-1,iy))/dx2)
     :    /(r*tbot(ix,iy))
        dlnpdy(ix,iy)=(-fcor*ug-g*(hsuf(ix,iy+1)-hsuf(ix,iy-1))/dy2)
     :    /(r*tbot(ix,iy))
        dlntdx(ix,iy)=(tbot(ix+1,iy)-tbot(ix-1,iy))/(dx2*tbot(ix,iy))
        dlntdy(ix,iy)=(tbot(ix,iy+1)-tbot(ix,iy-1))/(dy2*tbot(ix,iy))
      enddo
      enddo
c
      do iy=2,ny
      do ix=2,nx
        ff(ix,iy)=
     :    psuf(ix,iy)*(dlnpdx(ix,iy)*(dlnpdx(ix,iy)-dlntdx(ix,iy))
     :    +dlnpdy(ix,iy)*(dlnpdy(ix,iy)-dlntdy(ix,iy))-g/(r*tbot(ix,iy))
     :    *((hsuf(ix+1,iy)+hsuf(ix-1,iy)-2.*hsuf(ix,iy))/dxx
     :    +(hsuf(ix,iy+1)+hsuf(ix,iy-1)-2.*hsuf(ix,iy))/dyy))
      enddo
      enddo
c
c lateral boundary conditions:
c
      do iy=2,ny
        dpsudx(1,iy)=0.5*(psuf(1,iy)/tbot(1,iy)+psuf(2,iy)/tbot(2,iy))/r
     :    *(fcor*vg-g*(hsuf(2,iy)-hsuf(1,iy))/dx)
        dpsudx(2,iy)=0.5*(psuf(nx,iy)/tbot(nx,iy)+psuf(nx1,iy)
     :    /tbot(nx1,iy))/r
     :    *(fcor*vg-g*(hsuf(nx1,iy)-hsuf(nx,iy))/dx)
      enddo
      dpsudx(1,1)=dpsudx(1,2)
      dpsudx(1,ny1)=dpsudx(1,ny)
      dpsudx(2,1)=dpsudx(2,2)
      dpsudx(2,ny1)=dpsudx(2,ny)
c
      do ix=2,nx
        dpsudy(ix,1)=0.5*(psuf(ix,1)/tbot(ix,1)+psuf(ix,2)/tbot(ix,2))/r
     :    *(-fcor*ug-g*(hsuf(ix,2)-hsuf(ix,1))/dy)
        dpsudy(ix,2)=0.5*(psuf(ix,ny)/tbot(ix,ny)+psuf(ix,ny1)
     :    /tbot(ix,ny1))/r
     :    *(-fcor*ug-g*(hsuf(ix,ny1)-hsuf(ix,ny))/dy)
      enddo
      dpsudy(1,1)=dpsudy(2,1)
      dpsudy(nx1,1)=dpsudy(nx,1)
      dpsudy(1,2)=dpsudy(2,2)
      dpsudy(nx1,2)=dpsudy(nx,2)
c
      do iy=2,ny
        ff(2,iy)=ff(2,iy)+dpsudx(1,iy)/dx
        ff(nx,iy)=ff(nx,iy)-dpsudx(2,iy)/dx
      enddo
c
      do ix=2,nx
        ff(ix,2)=ff(ix,2)+dpsudy(ix,1)/dy
        ff(ix,ny)=ff(ix,ny)-dpsudy(ix,2)/dy
      enddo
c
c      ixref=2
c      iyref=ny1/2
c
c solve homogenous poisson equation
c
      if(phsout.and.prt) then
         tk=0.
         call wri2ar(ff,1,nx1,1,ny1,'ff      ',tk,nx,ny)
      endif
c

      call pois2n(ff,pref,ixref,iyref,hk21,hk22,dxx,dyy,eigxy
     :  ,trigsx,ifax,sipcox,simcox,nx
     :  ,trigsy,ifay,sipcoy,simcoy,ny)


c
c impose b.c.:
c
      do ix=2,nx
        ff(ix,1)=ff(ix,2)-dpsudy(ix,1)*dy
        ff(ix,ny1)=ff(ix,ny)+dpsudy(ix,2)*dy
      enddo
c
      do iy=1,ny1
        ff(1,iy)=ff(2,iy)-dpsudx(1,iy)*dx
        ff(nx1,iy)=ff(nx,iy)+dpsudx(2,iy)*dx
      enddo
c
      if(phsout.and.prt) then
         tk=0.
         call wri2ar(ff,1,nx1,1,ny1,'f       ',tk,nx,ny)
      endif
c
      do iy=1,ny1
      do ix=1,nx1
        psuf(ix,iy)=ff(ix,iy)
      enddo
      enddo

      deallocate(dpsudx,dpsudy)
      deallocate(ff,tbot)
      deallocate(dlnpdx,dlnpdy)
      deallocate(dlntdx,dlntdy)
      deallocate(hk21,hk22)
c
      return
      end



      subroutine tforci
c
c thermal forcing.
c
      use alloc
      implicit real*8(a-h,o-z)
c      include 'nh3dpa10.inc'
c
      parameter(relax1=0.05,relax2=0.5,ixle=8)
c
      return
c
      if(nstep.lt.100) then
         do iy=1,ny1
         do is=0,ns
         do ix=1,nx1
           pt(ix,iy,is,3)=pt(ix,iy,is,3)
     :       +relax1*(ptforc(ix,is)-pt(ix,iy,is,3))
         enddo
         enddo
         enddo
      elseif(nstep.lt.500) then
         do iy=1,ny1
         do ix=1,ixle
         do is=0,ns
           pt(ix,iy,is,3)=pt(ix,iy,is,3)
     :       +relax2*(ptforc(ix,is)-pt(ix,iy,is,3))
         enddo
         enddo
         enddo
      endif
c
      return
      end


      subroutine update
c-----------------------------------------------------------------------
c update variable fields
c-----------------------------------------------------------------------
      use alloc

      implicit real*8(a-h,o-z)
c      include 'nh3dpa10.inc'

c
c pp and related fields:
c
      do iy=0,ny2
      do ix=0,nx2
        pp(ix,iy,1)=pp(ix,iy,2)
        pp(ix,iy,2)=pp(ix,iy,3)
        pp(ix,iy,3)=pp(ix,iy,4)
      enddo
      enddo
c
      do iy=0,ny1
      do ix=0,nx1
        psuf(ix,iy)=pp(ix,iy,3)+ptop
        pp10(ix,iy,1)=pp10(ix,iy,2)
        pp10(ix,iy,2)=pp10(ix,iy,3)
        pp10(ix,iy,3)=pp10(ix,iy,4)
        pp01(ix,iy,1)=pp01(ix,iy,2)
        pp01(ix,iy,2)=pp01(ix,iy,3)
        pp01(ix,iy,3)=pp01(ix,iy,4)
        ppdx(ix,iy,1)=ppdx(ix,iy,2)
        ppdx(ix,iy,2)=ppdx(ix,iy,3)
        ppdx(ix,iy,3)=ppdx(ix,iy,4)
        ppdy(ix,iy,1)=ppdy(ix,iy,2)
        ppdy(ix,iy,2)=ppdy(ix,iy,3)
        ppdy(ix,iy,3)=ppdy(ix,iy,4)
      enddo
      enddo
c
c     Dima@ 27.04.07 
      if (impur.gt.0) then    
      do is=0,ns                                                                DM,05.2007
      do iy=0,ny1                                                               DM,05.2007
      do ix=0,nx1                                                               DM,05.2007
        qs_10(ix,iy,is,1)=qs_10(ix,iy,is,2)                                     DM,05.2007
        qs_10(ix,iy,is,2)=qs_10(ix,iy,is,3)                                     DM,05.2007
      enddo                                                                     DM,05.2007
      enddo                                                                     DM,05.2007
      enddo  
	endif                                                                      DM,05.2007

	

      if(qif.ne.0.) then
        do is=0,ns
        do iy=0,ny1
        do ix=0,nx1
          qv(ix,iy,is,1)=qv(ix,iy,is,2)
          qv(ix,iy,is,2)=qv(ix,iy,is,3)
          if(ifqc.ne.0) then
            qc(ix,iy,is,1)=qc(ix,iy,is,2)
            qc(ix,iy,is,2)=qc(ix,iy,is,3)
            if(ifqr.ne.0) then
              qr(ix,iy,is,1)=qr(ix,iy,is,2)
              qr(ix,iy,is,2)=qr(ix,iy,is,3)
	        if(ifqi.ne.0) then                                                DC,11.2009
	          qci(ix,iy,is,1)=qci(ix,iy,is,2)                                 DC,11.2009
                qci(ix,iy,is,2)=qci(ix,iy,is,3)                                 DC,11.2009
	          qsn(ix,iy,is,1)=qsn(ix,iy,is,2)                                 DC,11.2009
                qsn(ix,iy,is,2)=qsn(ix,iy,is,3)                                 DC,11.2009
	        endif
            endif
          endif
        enddo
        enddo
        enddo
      endif
c
c w and wsig:
c
      do is=0,ns1
      do iy=0,ny1
      do ix=0,nx1
        wsig(ix,iy,is,2)=wsig(ix,iy,is,3)
      enddo
      enddo
      enddo
c
      return
      end

      subroutine update1
c-----------------------------------------------------------------------
c update variable fields
c-----------------------------------------------------------------------
      use alloc

      implicit real*8(a-h,o-z)
c      include 'nh3dpa10.inc'

c
c pp and related fields:
c
      do iy=0,ny2
      do ix=0,nx2
        pp(ix,iy,1)=pp(ix,iy,2)
        pp(ix,iy,2)=pp(ix,iy,3)
        pp(ix,iy,3)=pp(ix,iy,4)
      enddo
      enddo
c
      do iy=0,ny1
      do ix=0,nx1
        psuf(ix,iy)=pp(ix,iy,3)+ptop
        pp10(ix,iy,1)=pp10(ix,iy,2)
        pp10(ix,iy,2)=pp10(ix,iy,3)
        pp10(ix,iy,3)=pp10(ix,iy,4)
        pp01(ix,iy,1)=pp01(ix,iy,2)
        pp01(ix,iy,2)=pp01(ix,iy,3)
        pp01(ix,iy,3)=pp01(ix,iy,4)
        ppdx(ix,iy,1)=ppdx(ix,iy,2)
        ppdx(ix,iy,2)=ppdx(ix,iy,3)
        ppdx(ix,iy,3)=ppdx(ix,iy,4)
        ppdy(ix,iy,1)=ppdy(ix,iy,2)
        ppdy(ix,iy,2)=ppdy(ix,iy,3)
        ppdy(ix,iy,3)=ppdy(ix,iy,4)
      enddo
      enddo
c
      if(qif.ne.0.) then
        do is=0,ns
        do iy=0,ny1
        do ix=0,nx1
          qv(ix,iy,is,1)=qv(ix,iy,is,2)
          qv(ix,iy,is,2)=qv(ix,iy,is,3)
          if(ifqc.ne.0) then
            qc(ix,iy,is,1)=qc(ix,iy,is,2)
            qc(ix,iy,is,2)=qc(ix,iy,is,3)
            if(ifqr.ne.0) then
              qr(ix,iy,is,1)=qr(ix,iy,is,2)
              qr(ix,iy,is,2)=qr(ix,iy,is,3)
	        if(ifqi.ne.0) then                                                DC,11.2009
	          qci(ix,iy,is,1)=qci(ix,iy,is,2)                                 DC,11.2009
                qci(ix,iy,is,2)=qci(ix,iy,is,3)                                 DC,11.2009
	          qsn(ix,iy,is,1)=qsn(ix,iy,is,2)                                 DC,11.2009
                qsn(ix,iy,is,2)=qsn(ix,iy,is,3)                                 DC,11.2009
	        endif
            endif
          endif
        enddo
        enddo
        enddo
      endif
c
      return
      end









      subroutine uprefst(istepp)

      use alloc
      use refstate

      implicit real*8(a-h,o-z)
      character*80 line,fndebug
      integer year1,month1,day1,hour1,min1,sec1
      integer year2,month2,day2,hour2,min2,sec2
      real(kind(0.d0)),allocatable,dimension(:,:,:)::
     :  xflop,yflop,ubxpro1,vbypro1,rhox,rhoy
     :  ,ubxpro2,vbypro2,dubxprodt,dvbyprodt
     :  ,dubxprodz1,dubxprodz2,dvbyprodz1,dvbyprodz2
     :  ,dubxprodzdt,dvbyprodzdt

      if(verbose.ge.2)
     :   write(*,*) 'uprefst:',istepp,intrefst

      if(istepp.eq.intrefst) then
        if(verbose.ge.1) write(*,*) 'uprefst: cycle',istepp

        read(35,'(a)',end=8888) line
        read(35,*) tprof2
        do ip=1,nprof
          read(35,'(a)',end=8888) line
          if(verbose.ge.2) write(*,*) 'Read prof:',ip,line
          read(35,*) xprof(ip),yprof(ip)
          read(35,'(a)',end=8888) line
          do ii=0,ndat-1
            read(35,*) zthdat(ii,ip),thdat(ii,ip)
     :        ,usdat(ii,ip),vsdat(ii,ip),qvsdat(ii,ip),pressdat(ii,ip)
          enddo
          zthdat(ndat,ip)=max(20000.,zthdat(ndat-1,ip))
          usdat(ndat,ip)=usdat(ndat-1,ip)
          vsdat(ndat,ip)=vsdat(ndat-1,ip)
          qvsdat(ndat,ip)=qvsdat(ndat-1,ip)
          thdat(ndat,ip)=thdat(ndat-1,ip)
          pressdat(ndat,ip)=100.  ! not relevant?
        enddo

        zusdat=zthdat
        zvsdat=zthdat
        zqvsdat=zthdat

        fndebug='NewWeipro.grd'
        call weiprof(fndebug)
c
c mrefst=timesteps until next update
c intrefst=nstep of update
c
        if(realtime) then
c         write(*,*) 'uprefst:',tprof1,tprof2
          write(line,'(f15.0)') tprof1
          write(*,*) line
          read(line,'(i4,5i2)') year1,month1,day1,hour1,min1,sec1
          write(line,'(f15.0)') tprof2
          read(line,'(i4,5i2)') year2,month2,day2,hour2,min2,sec2
          h1=julian(year1,month1,day1)
          h2=julian(year2,month2,day2)
          stepseconds=((h2*24+hour2-h1*24-hour1)*3600.
     :      +(min2*60+sec2-min1*60-sec1))
          mrefst=nint(stepseconds/dt)
        else
          mrefst=nint((tprof2-tprof1)/dt)
        endif

        intrefst=istepp+mrefst
        if(verbose.ge.1) write(*,*) 'uprefst: intrefst',intrefst
        write(nchn1,'(a20,2f15.0,2i8)')
     :    'Next bc profile:',tprof1,tprof2,istepp,intrefst

c save current ref profiles

        thpro1=thpro0
        uspro1=uspro0
        vspro1=vspro0
        qvspro1=qvspro0
        press1=press0
        rhopro1=rhopro

        dthdz1=dthdz
        dusdz1=dusdz
        dvsdz1=dvsdz
        dqvsdz1=dqvsdz

        call refprofs

        do ip=1,nprof
          dpressdt(ip)=(press0(ip)-press1(ip))/stepseconds
          do iz=0,npro
            dthdt(iz,ip)=(thpro0(iz,ip)-thpro1(iz,ip))/mrefst
            dthdzdt(iz,ip)=(dthdz(iz,ip)-dthdz1(iz,ip))/mrefst
            dusdt(iz,ip)=(uspro0(iz,ip)-uspro1(iz,ip))/mrefst
            dusdzdt(iz,ip)=(dusdz(iz,ip)-dusdz1(iz,ip))/mrefst
            dvsdt(iz,ip)=(vspro0(iz,ip)-vspro1(iz,ip))/mrefst
            dvsdzdt(iz,ip)=(dvsdz(iz,ip)-dvsdz1(iz,ip))/mrefst
            dqvsdt(iz,ip)=(qvspro0(iz,ip)-qvspro1(iz,ip))/mrefst
            dqvsdzdt(iz,ip)=(dqvsdz(iz,ip)-dqvsdz1(iz,ip))/mrefst
          enddo
        enddo

! verify mass balance in the reference state profiles
! computations at time zero (of the update cycle)

        if(divrefstate.and. nprof.gt.1) then

          dppref=0.
          dpmeddt=0.
          do ip=1,nprof
            do ix=0,nx1
            do iy=0,ny1
              dppref(ix,iy)=dppref(ix,iy)+pweight(ix,iy,ip)*dpressdt(ip)
            enddo
            enddo
            do ix=2,nx
            do iy=2,ny
              dpmeddt=dpmeddt+dppref(ix,iy)
            enddo
            enddo
          enddo

          if(.not.allocated(ubxpro1)) allocate(ubxpro1(2,2:ny,0:npro))
          if(.not.allocated(ubxpro2)) allocate(ubxpro2(2,2:ny,0:npro))
          if(.not.allocated(dubxprodt))
     :      allocate(dubxprodt(2,2:ny,0:npro))
          if(.not.allocated(dubxprodz1))
     :      allocate(dubxprodz1(2,2:ny,0:npro))
          if(.not.allocated(dubxprodz2))
     :      allocate(dubxprodz2(2,2:ny,0:npro))
          if(.not.allocated(dubxprodzdt))
     :      allocate(dubxprodzdt(2,2:ny,0:npro))
          if(.not.allocated(rhox)) allocate(rhox(2,2:ny,0:npro))
          if(.not.allocated(xflop)) allocate(xflop(2,2:ny,0:npro))

          if(.not.allocated(vbypro1)) allocate(vbypro1(2:nx,2,0:npro))
          if(.not.allocated(vbypro2)) allocate(vbypro2(2:nx,2,0:npro))
          if(.not.allocated(dvbyprodt))
     :      allocate(dvbyprodt(2:nx,2,0:npro))
          if(.not.allocated(dvbyprodz1))
     :      allocate(dvbyprodz1(2:nx,2,0:npro))
          if(.not.allocated(dvbyprodz2))
     :      allocate(dvbyprodz2(2:nx,2,0:npro))
          if(.not.allocated(dvbyprodzdt))
     :      allocate(dvbyprodzdt(2:nx,2,0:npro))
          if(.not.allocated(rhoy)) allocate(rhoy(2:nx,2,0:npro))
          if(.not.allocated(yflop)) allocate(yflop(2:nx,2,0:npro))


          ubxpro1=0.
          rhox=0.
          vbypro1=0.
          rhoy=0.

          do ip=1,nprof
            do iz=0,npro
            do iy=2,ny
              ubxpro1(1,iy,iz)=ubxpro1(1,iy,iz)
     :          +pweight(1,iy,ip)*(uspro1(iz,ip)+dusdz1(iz,ip)*zpro(iz))
              ubxpro1(2,iy,iz)=ubxpro1(2,iy,iz)
     :          +pweight(nx,iy,ip)*(uspro1(iz,ip)
     :          +dusdz1(iz,ip)*zpro(iz))

              rhox(1,iy,iz)=rhox(1,iy,iz)
     :          +pweight(1,iy,ip)*rhopro1(iz,ip)
              rhox(2,iy,iz)=rhox(2,iy,iz)
     :          +pweight(nx,iy,ip)*rhopro1(iz,ip)
            enddo
            enddo

            do iz=0,npro
            do ix=2,nx
              vbypro1(ix,1,iz)=vbypro1(ix,1,iz)
     :          +pweight(ix,1,ip)*(vspro1(iz,ip)+dvsdz1(iz,ip)*zpro(iz))
              vbypro1(ix,2,iz)=vbypro1(ix,2,iz)
     :          +pweight(ix,ny,ip)*(vspro1(iz,ip)
     :          +dvsdz1(iz,ip)*zpro(iz))

              rhoy(ix,1,iz)=rhoy(ix,1,iz)
     :          +pweight(ix,1,ip)*rhopro1(iz,ip)
              rhoy(ix,2,iz)=rhoy(ix,2,iz)
     :          +pweight(ix,ny,ip)*rhopro1(iz,ip)
            enddo
            enddo
          enddo

          flowinx=0.
          flowoutx=0.
          flowiny=0.
          flowouty=0.

          do iy=2,ny
          do iz=0,npro
            xflop(1,iy,iz)=ubxpro1(1,iy,iz)*rhox(1,iy,iz)
            xflop(2,iy,iz)=ubxpro1(2,iy,iz)*rhox(2,iy,iz)
          enddo
          enddo

          do iy=2,ny
          do iz=0,npro
            flowinx=flowinx+max(0.d0,xflop(1,iy,iz))
     :         -min(0.d0,xflop(2,iy,iz))
            flowoutx=flowoutx-min(0.d0,xflop(1,iy,iz))
     :         +max(0.d0,xflop(2,iy,iz))
          enddo
          enddo

          do iz=0,npro
          do ix=2,nx
            yflop(ix,1,iz)=vbypro1(ix,1,iz)*rhoy(ix,1,iz)
            yflop(ix,2,iz)=vbypro1(ix,2,iz)*rhoy(ix,2,iz)
          enddo
          enddo

          do iz=0,npro
          do ix=2,nx
            flowiny=flowiny+max(0.d0,yflop(ix,1,iz))
     :        -min(0.d0,yflop(ix,2,iz))
            flowouty=flowouty-min(0.d0,yflop(ix,1,iz))
     :        +max(0.d0,yflop(ix,2,iz))
          enddo
          enddo

          flowin=flowinx+flowiny
          flowout=flowoutx+flowouty
          area=(nx-1)*dx*(ny-1)*dy

          residual=flowin-flowout-dpmeddt*area
          totalcf=flowin+flowout

          if(residual.gt.1.e-6) then
            do iy=2,ny
            do iz=0,npro
              ubxpro1(1,iy,iz)=ubxpro1(1,iy,iz)
     :          -residual*abs(xflop(1,iy,iz))/(totalcf*rhox(1,iy,iz))
              ubxpro1(2,iy,iz)=ubxpro1(2,iy,iz)
     :          -residual*abs(xflop(2,iy,iz))/(totalcf*rhox(2,iy,iz))
            enddo
            enddo

            do iz=0,npro
            do ix=2,nx
              vbypro1(ix,1,iz)=vbypro1(ix,1,iz)
     :          -residual*abs(yflop(ix,1,iz))/(totalcf*rhoy(ix,1,iz))
              vbypro1(ix,2,iz)=vbypro1(ix,2,iz)
     :          -residual*abs(yflop(ix,2,iz))/(totalcf*rhoy(ix,2,iz))
            enddo
            enddo
          endif


! veri  fy mass balance in the reference state profiles
! comp  utations at time max (of the update cycle)

          ubxpro2=0.
          vbypro2=0.
          rhox=0.
          rhoy=0.

          do ip=1,nprof

            do iz=0,npro
            do iy=2,ny
              ubxpro2(1,iy,iz)=ubxpro2(1,iy,iz)
     :          +pweight(1,iy,ip)*(uspro0(iz,ip)+dusdz(iz,ip)*zpro(iz))
              ubxpro2(2,iy,iz)=ubxpro2(2,iy,iz)
     :          +pweight(nx,iy,ip)*(uspro0(iz,ip)+dusdz(iz,ip)*zpro(iz))

              rhox(1,iy,iz)=rhox(1,iy,iz)
     :          +pweight(1,iy,ip)*rhopro(iz,ip)
              rhox(2,iy,iz)=rhox(2,iy,iz)
     :          +pweight(nx,iy,ip)*rhopro(iz,ip)
            enddo
            enddo

            do iz=0,npro
            do ix=2,nx
              vbypro2(ix,1,iz)=vbypro2(ix,1,iz)
     :          +pweight(ix,1,ip)*(vspro0(iz,ip)+dvsdz(iz,ip)*zpro(iz))
              vbypro2(ix,2,iz)=vbypro2(ix,2,iz)
     :          +pweight(ix,ny,ip)*(vspro0(iz,ip)+dvsdz(iz,ip)*zpro(iz))

              rhoy(ix,1,iz)=rhoy(ix,1,iz)
     :          +pweight(ix,1,ip)*rhopro(iz,ip)
              rhoy(ix,2,iz)=rhoy(ix,2,iz)
     :          +pweight(ix,ny,ip)*rhopro(iz,ip)
            enddo
            enddo
          enddo

          flowinx=0.
          flowoutx=0.
          flowiny=0.
          flowouty=0.

          do iy=2,ny
          do iz=0,npro
            xflop(1,iy,iz)=ubxpro2(1,iy,iz)*rhox(1,iy,iz)
            xflop(2,iy,iz)=ubxpro2(2,iy,iz)*rhox(2,iy,iz)
          enddo
          enddo

          do iy=2,ny
          do iz=0,npro
            flowinx=flowinx+max(0.,xflop(1,iy,iz))
     :         -min(0.,xflop(2,iy,iz))
            flowoutx=flowoutx-min(0.,xflop(1,iy,iz))
     :         +max(0.,xflop(2,iy,iz))
          enddo
          enddo

          do iz=0,npro
          do ix=2,nx
            yflop(ix,1,iz)=vbypro2(ix,1,iz)*rhoy(ix,1,iz)
            yflop(ix,2,iz)=vbypro2(ix,2,iz)*rhoy(ix,2,iz)
          enddo
          enddo

          do iz=0,npro
          do ix=2,nx
            flowiny=flowiny+max(0.,yflop(ix,1,iz))
     :        -min(0.,yflop(ix,2,iz))
            flowouty=flowouty-min(0.,yflop(ix,1,iz))
     :        +max(0.,yflop(ix,2,iz))
          enddo
          enddo

          flowin=flowinx+flowiny
          flowout=flowoutx+flowouty
          area=(nx-1)*dx*(ny-1)*dy

          dppref=0.
          psufref=0.
          dpmeddt=0.

c note  : 1/g cancels

          residual=flowin-flowout-dpmeddt*area
          totalcf=flowin+flowout

c         write(*,*) 'residual:',residual,flowin,flowout,dpmeddt*area

          if(residual.gt.1.e-6) then
            do iy=2,ny
            do iz=0,npro
              ubxpro2(1,iy,iz)=ubxpro2(1,iy,iz)
     :          -residual*abs(xflop(1,iy,iz))/(totalcf*rhox(1,iy,iz))
              ubxpro2(2,iy,iz)=ubxpro2(2,iy,iz)
     :          -residual*abs(xflop(2,iy,iz))/(totalcf*rhox(2,iy,iz))
            enddo
            enddo

            do ix=2,nx
            do iz=0,npro
              vbypro2(ix,1,iz)=vbypro2(ix,1,iz)
     :          -residual*abs(yflop(ix,1,iz))/(totalcf*rhoy(ix,1,iz))
              vbypro2(ix,2,iz)=vbypro2(ix,2,iz)
     :          -residual*abs(yflop(ix,2,iz))/(totalcf*rhoy(ix,2,iz))
            enddo
            enddo
          endif

          do ii=1,2
          do iy=2,ny
          do iz=0,npro
            dubxprodt(1,iy,iz)=(ubxpro2(1,iy,iz)-ubxpro1(1,iy,iz))
     :        /mrefst
          enddo
          enddo
          enddo

          do ii=1,2
          do iy=2,ny
          do iz=0,npro-1
            dubxprodz1(ii,iy,iz)=(ubxpro1(ii,iy,iz+1)-ubxpro1(ii,iy,iz))
     :        /dzpro
            dubxprodz2(ii,iy,iz)=(ubxpro2(ii,iy,iz+1)-ubxpro2(ii,iy,iz))
     :        /dzpro
          enddo
          enddo
          enddo

          do ii=1,2
          do iy=2,ny
          do iz=0,npro-1
            dubxprodzdt(ii,iy,iz)=(dubxprodz2(ii,iy,iz)
     :        -dubxprodz1(ii,iy,iz))/mrefst
          enddo
          enddo
          enddo

          do ii=1,2
          do ix=2,nx
          do iz=0,npro
            dvbyprodt(ix,ii,iz)=(vbypro2(ix,ii,iz)-vbypro1(ix,ii,iz))
     :        /mrefst
          enddo
          enddo
          enddo

          do ii=1,2
          do ix=2,nx
          do iz=0,npro-1
            dvbyprodz1(ix,ii,iz)=(vbypro1(ix,ii,iz+1)
     :        -vbypro1(ix,ii,iz))/mrefst
            dvbyprodz2(ix,ii,iz)=(vbypro2(ix,ii,iz+1)
     :        -vbypro2(ix,ii,iz))/mrefst
          enddo
          enddo
          enddo

          do ii=1,2
          do ix=2,nx
          do iz=0,npro-1
            dvbyprodzdt(ix,ii,iz)=(dvbyprodz2(ix,ii,iz)
     :        -dvbyprodz1(ix,ii,iz))/dzpro
          enddo
          enddo
          enddo

          if(.not.allocated(dweightx)) then
            allocate(dweightx(0:nx+1,0:ny+1,2))
            allocate(dweighty(0:nx+1,0:ny+1,2))
            do ix=0,nx+1
              distx2=(ix-nx)*dx
              distx1=(ix-2)*dx
              do iy=0,ny+1
                disty2=(iy-ny)*dy
                disty1=(iy-2)*dy
                distinv11=1./sqrt(distx1**2+disty1**2)
                distinv12=1./sqrt(distx1**2+disty2**2)
                distinv21=1./sqrt(distx2**2+disty1**2)
                distinv22=1./sqrt(distx2**2+disty2**2)
                sumdisti=distinv11+distinv12+distinv21+distinv22
                dweightx(ix,iy,1)=distinv11/sumdisti
                dweightx(ix,iy,2)=distinv21/sumdisti
                dweighty(ix,iy,1)=distinv12/sumdisti
                dweighty(ix,iy,2)=distinv22/sumdisti
              enddo
            enddo
          endif
        endif

! end_divrefstate

        if(verbose.ge.2) then
          do ip=1,nprof
            write(*,*) 'New profile',ip,dpressdt(ip)
            do iz=1,npro
            write(*,'(i5,8e15.7)')
     :        iz,thpro0(iz,ip)
     :        ,uspro0(iz,ip),vspro0(iz,ip)
     :        ,qvspro0(iz,ip)
            enddo
          enddo
        endif

c compute tendencies

        do ip=1,nprof
          do iz=0,npro
            dusdt(iz,ip)=(uspro0(iz,ip)-uspro1(iz,ip))/mrefst
            dusdzdt(iz,ip)=(dusdz(iz,ip)-dusdz1(iz,ip))/mrefst
            dvsdt(iz,ip)=(vspro0(iz,ip)-vspro1(iz,ip))/mrefst
            dvsdzdt(iz,ip)=(dvsdz(iz,ip)-dvsdz1(iz,ip))/mrefst
          enddo
        enddo

c recover saved profiles

        thpro0=thpro1
        uspro0=uspro1
        vspro0=vspro1
        qvspro0=qvspro1
        press0=press1

        dthdz=dthdz1
        dusdz=dusdz1
        dvsdz=dvsdz1
        dqvsdz=dqvsdz1

        if(verbose.ge.2) then
          do ip=1,nprof
            write(*,'(a30,i5,2e15.7)') 'Uprefst Next current profile'
     :        ,ip,press0(ip),dpressdt(ip)
            do iz=1,npro
            write(*,'(i5,8e15.7)')
     :        iz,thpro0(iz,ip),dthdt(iz,ip)
     :        ,uspro0(iz,ip),dusdt(iz,ip),vspro0(iz,ip),dvsdt(iz,ip)
     :        ,qvspro0(iz,ip),dqvsdt(iz,ip)
            enddo
          enddo
        endif

      endif


      do ip=1,nprof
        press0(ip)=press0(ip)+dpressdt(ip)*dt
        do iz=0,npro
          thpro0(iz,ip)=thpro0(iz,ip)+dthdt(iz,ip)
          dthdz(iz,ip)=dthdz(iz,ip)+dthdzdt(iz,ip)
          uspro0(iz,ip)=uspro0(iz,ip)+dusdt(iz,ip)
          dusdz(iz,ip)=dusdz(iz,ip)+dusdzdt(iz,ip)
          vspro0(iz,ip)=vspro0(iz,ip)+dvsdt(iz,ip)
          dvsdz(iz,ip)=dvsdz(iz,ip)+dvsdzdt(iz,ip)
          qvspro0(iz,ip)=qvspro0(iz,ip)+dqvsdt(iz,ip)
          dqvsdz(iz,ip)=dqvsdz(iz,ip)+dqvsdzdt(iz,ip)
        enddo
      enddo

      if(verbose.ge.3) then
        do ip=1,nprof
          write(*,'(a20,i6,e15.7)') 'New profile',ip,press0(ip)
          do iz=1,npro
          write(*,'(i5,4f10.3)')
     :      iz,thpro0(iz,ip),uspro0(iz,ip),vspro0(iz,ip),qvspro0(iz,ip)
          enddo
        enddo
      endif

      return
8888  continue
      write(0,*) 'nh3dFatalError=ERROR in uprefst (end of file)'
      stop
      end


      function julian(year,month,day)
      integer julian,year,month,day
      integer nm(12)
      data nm/31,28,31,30,31,30,31,31,30,31,30,31/
      if(mod(year,4).eq.0) then
        nm(2)=29
      else
        nm(2)=28
      endif

      julian=day
      do im=1,month-1
        julian=julian+nm(im)
      enddo

      return
      end

      subroutine timeofyear(styear,stmonth,stday,sthour,stmin,stsec
     :  ,elapsedsec,cuyear,cumonth,cuday,cuhour,cumin,cusec)

! compute clock time from number of seconds in the run (for output only)

      implicit none
      integer styear,stmonth,stday,sthour,stmin,stsec
      integer cuyear,cumonth,cuday,cuhour,cumin,cusec
      integer elapseddays,elapsedsec,cujulian,julianstart,newyear
     :  ,newmonth,newday,julian,elapsedhours

      elapseddays=elapsedsec/86400
      elapsedhours=(elapsedsec-elapseddays*86400)/3600
!      elapsedmin=(elapsedsec-elapseddays*86400-elapsedhours*3600)/60
      if(sthour+elapsedhours.ge.24) elapseddays=elapseddays+1

      julianstart=julian(styear,stmonth,stday)
      cujulian=elapseddays+julianstart

c verify cross of year boundary

      if(elapseddays.le.julian(styear,12,31)-julian(styear
     :  ,stmonth,stday)) then
!       write(*,'(a10,3i8)') 'Timeof:',julianstart,elapseddays,cujulian
        call julian2ymd(cujulian,styear,cumonth,cuday)
        cusec=elapsedsec-elapseddays*86400
        cuhour=cusec/3600+sthour
        cumin=mod(elapsedsec,3600)/60
        cusec=mod(elapsedsec,60)
        cuyear=styear
      else
        newyear=styear+1
        newmonth=stmonth
        newday=stday
        elapseddays=elapseddays-julian(styear,12,31)+julian(styear
     :    ,stmonth,stday)
        do while(elapseddays.gt.julian(newyear,12,31))
          elapseddays=elapseddays-julian(newyear,12,31)
          elapsedsec=elapsedsec-julian(newyear,12,31)*86400
          newyear=newyear+1
        enddo
        call julian2ymd(elapseddays,newyear,cumonth,cuday)
        cusec=elapsedsec-elapseddays*86400
        cuhour=cusec/3600+sthour
        cumin=mod(elapsedsec,3600)/60
        cusec=mod(elapsedsec,60)
        cuyear=newyear
      endif
!     write(*,'(''time:'',12i4,i10)')
!    :  styear,stmonth,stday,sthour,stmin,stsec
!    :  ,cuyear,cumonth,cuday,cuhour,cumin,cusec
!    :  ,elapsedsec
      return
      end

      subroutine julian2ymd(julian,year,month,day)
      integer days,year,month,day
      integer nm(12)
      data nm/31,28,31,30,31,30,31,31,30,31,30,31/
      if(mod(year,4).eq.0) then
        nm(2)=29
      else
        nm(2)=28
      endif

      days=julian

      month=1
      do while(days.gt.nm(month))
        days=days-nm(month)
        month=month+1
      enddo

      day=days

!     write(*,*) 'julian2:',julian,year,month,day
      return
      end

      subroutine upuvt
c
c update all reference fields :
c

      use alloc
      use refstate

      implicit real*8(a-h,o-z)
c
      us=0.
      vs=0.
      if(qif.ne.0.) qvs=0.
      pts=0.
      tems=0.
      dppref=0.
      psufref=0.
      dpmeddt=0.

      if(verbose.ge.2) write(*,*) 'Upuvt:',nprof
      do ip=1,nprof

        do ix=0,nx1
        do iy=0,ny1
          dppref(ix,iy)=dppref(ix,iy)+pweight(ix,iy,ip)*dpressdt(ip)
          psufref(ix,iy)=psufref(ix,iy)+pweight(ix,iy,ip)*press0(ip)
        enddo
        enddo

        do ix=2,nx
        do iy=2,ny
          dpmeddt=dpmeddt+dppref(ix,iy)
        enddo
        enddo


        do is=1,ns-1
        do iy=0,ny1
        do ix=1,nx
          zsia=0.5*(phis(ix,iy,is)+phis(ix+1,iy,is))/g-
     :      0.25*(phis(ix,iy,ns)+phis(ix,iy,ns-1)+phis(ix+1,iy,ns)
     :      +phis(ix+1,iy,ns-1))/g
          izia=zsia/dzpro
          izia=max(0,min(izia,npro))
          us(ix,iy,is)=us(ix,iy,is)
     :      +pweight(ix,iy,ip)*(uspro0(izia,ip)+dusdz(izia,ip)*zsia)
        enddo
        enddo
        enddo
c
        if(ip.eq.nprof) then
          do is=1,ns-1
          do iy=0,ny1
            us(0,iy,is)=us(1,iy,is)
            us(nx1,iy,is)=us(nx,iy,is)
          enddo
          enddo
        endif
c
        do is=1,ns-1
        do iy=1,ny
        do ix=0,nx1
          zsia=0.5*(phis(ix,iy,is)+phis(ix,iy+1,is))/g-
     :      0.5*(phis(ix,iy,ns)+phis(ix,iy,ns-1)+phis(ix,iy+1,ns)
     :      +phis(ix,iy+1,ns-1))/g
          izia=zsia/dzpro
          izia=max(0,min(izia,npro))
          vs(ix,iy,is)=vs(ix,iy,is)
     :      +pweight(ix,iy,ip)*(vspro0(izia,ip)+dvsdz(izia,ip)*zsia)
        enddo
        enddo
        enddo
c
        if(ip.eq.nprof) then
          do is=1,ns-1
          do ix=0,nx1
            vs(ix,0,is)=vs(ix,1,is)
            vs(ix,ny1,is)=vs(ix,ny,is)
          enddo
          enddo
c
          do iy=0,ny1
          do ix=0,nx1
            us(ix,iy,0)=us(ix,iy,1)
            vs(ix,iy,0)=vs(ix,iy,1)
            us(ix,iy,ns)=us(ix,iy,ns-1)
            vs(ix,iy,ns)=vs(ix,iy,ns-1)
          enddo
          enddo
        endif

        if(qif.ne.0.) then
          do is=0,ns
          do iy=0,ny1
          do ix=0,nx1
    !        zsia=phis(ix,iy,is)/g
            zsia=phis(ix,iy,is)/g-
     :      0.5*(phis(ix,iy,ns)+phis(ix,iy,ns-1))/g
            izia=zsia/dzpro
            izia=max(0,min(izia,npro))
            qvs(ix,iy,is)=qvs(ix,iy,is)
     :        +pweight(ix,iy,ip)*(qvspro0(izia,ip)+dqvsdz(izia,ip)*zsia)
          enddo
          enddo
          enddo
        endif
c
c pts   only:
c unne  cessary (disabled)
c        write(*,*) 'entry'
c        entry upt
c
        if(ioreft.ne.2) then
          do is=0,ns
          do iy=0,ny1
          do ix=0,nx1
 !           zsia=phis(ix,iy,is)/g
            
  !          zsia=phis(ix,iy,is)/g-
  !   :      0.5*(phis(ix,iy,ns)+phis(ix,iy,ns-1))/g
           
            zsia=(phis(ix,iy,is)+phi(ix,iy,is))/g
            
            izia0=zsia/dzpro
            izia=max(0,min(izia0,npro))
            pts(ix,iy,is)=pts(ix,iy,is)
     :        +pweight(ix,iy,ip)*(thpro0(izia,ip)+dthdz(izia,ip)*zsia)
          !if(ix.eq.nx/2.and.is.eq.ns-10)write(0,*)zsia,pts(ix,iy,is)
            if(verbose.ge.4) write(*,'(a10,4i4,f8.0,i4,3e15.7)')
     :        'inter:',
     :        ip,ix,iy,is,zsia,izia,thpro0(izia,ip),dthdz(izia,ip)
     :        ,pts(ix,iy,is)
          enddo
          enddo
          enddo
c          do 120 iy=0,ny1
c          do 120 ix=0,nx1
c          pts(ix,iy,ns)=pts(ix,iy,ns-1)
c          pts(ix,iy,0)=pts(ix,iy,1)
c120       continue
          if(ip.eq.nprof) then
            do is=0,ns
            do iy=0,ny1
            do ix=0,nx1
              tems(ix,iy,is)=pts(ix,iy,is)*((ptop+pp(ix,iy,3)
     :           *sigma0(is))/p00)**akapa
              if(qif.ne.0.) then
                tems(ix,iy,is)=tems(ix,iy,is)*(1.+0.61*qvs(ix,iy,is))
              endif
            enddo
            enddo
            enddo
            do iy=0,ny1
            do ix=0,nx1
c            tsuf(ix,iy)=(pts(ix,iy,ns)+pts(ix,iy,ns-1))*0.5
c         :     *(psuf(ix,iy)/p00)**akapa
              tsuf(ix,iy)=(tems(ix,iy,ns)+tems(ix,iy,ns-1))*0.5              
            enddo
            enddo
          endif
c
        else
          do is=0,ns
          do iy=0,ny1
          do ix=0,nx1
            zsia=phis(ix,iy,is)/g-
     :      0.5*(phis(ix,iy,ns)+phis(ix,iy,ns-1))/g
            izia=zsia/dzpro
            izia=max(0,min(izia,npro))
            tems(ix,iy,is)=tems(ix,iy,is)
     :        +pweight(ix,iy,ip)*(thpro0(izia,ip)+dthdz(izia,ip)*zsia)
          enddo
          enddo
          enddo

          if(ip.eq.nprof) then
            do is=0,ns
            do iy=0,ny1
            do ix=0,nx1
              p=(ptop+pp(ix,iy,3)*sigma0(is))
              if(p.le.0) write(*,*) ix,iy,is,pp(ix,iy,3),p,sigma0(is)
              pts(ix,iy,is)=tems(ix,iy,is)*(p00/p)**akapa
              if(qif.ne.0.) then
                tems(ix,iy,is)=tems(ix,iy,is)*(1.+0.61*qvs(ix,iy,is))
              endif
            enddo
            enddo
            enddo
c            do 121 iy=0,ny1
c            do 121 ix=0,nx1
c            pts(ix,iy,ns)=pts(ix,iy,ns-1)
c            pts(ix,iy,0)=pts(ix,iy,1)
c121         continue
c            do 122 iy=0,ny1
c            do 122 ix=0,nx1
c            tems(ix,iy,ns)=pts(ix,iy,ns)*((ptop+pp(ix,iy,3)*sigma0(ns))
c         :     /p00)**akapa
c            tems(ix,iy,0)=pts(ix,iy,0)*((ptop+pp(ix,iy,3)*sigma0(0))
c         :     /p00)**akapa
c122         continue
            do iy=0,ny1
            do ix=0,nx1
              tsuf(ix,iy)=0.5*(tems(ix,iy,ns-1)+tems(ix,iy,ns))
            enddo
            enddo
          endif
        endif
      enddo

      dpmeddt=dpmeddt/((nx-1)*(ny-1))

      if(verbose.ge.2) write(*,*) 'dpmetdt=',dpmeddt
c
c
c     tkoff=0.
c     call wri3ar(pts,1,nx1,1,ny1,0,ns,'upt-pts ',tkoff,0,nx,ny)
c     call wri3ar(tems,1,nx1,1,ny1,0,ns,'upt-tems',tkoff,0,nx,ny)
c
      return
c
      end









      subroutine wrigar(a,mx0,mx1,my0,my1,ms0,ms1,i0,i1,j0,j1,k0,k1
     :   ,title,tkoff,mode)
c
c-----------------------------------------------------------------------
c handler for outarr
c also grid output
c-----------------------------------------------------------------------
c
      use alloc
      use nh3dpa2
      implicit real*8(a-h,o-z)
c      include 'nh3dpa10.inc'
      character*8 title
      dimension a(mx0:mx1,my0:my1,ms0:*)
      parameter(maxnxny=1026)
      integer iap(0:maxnxny)
      logical firstcall
      data firstcall/.true./

      if(firstcall) then
        firstcall=.false.
        do i=0,maxnxny
          iap(i)=i
        enddo
      endif

c
c calculates max and min for scaling:
c
      armin=a(i0,j0,k0)
      armax=armin
      do k=k0,k1
      do j=j0,j1
      do i=i0,i1
        armin=min(a(i,j,k),armin)
        armax=max(a(i,j,k),armax)
      enddo
      enddo
      enddo
c
c calculates the scaling:
c
      write(nunit,1001) title,armin,armax
      if(armax.eq.armin .or. mode.eq.4) return
c
      armax=armax-tkoff
      armin=armin-tkoff
c
      scale=zero0+max(abs(armax),abs(armin))
      nexp=int(dlog10(scale))-ndigit
      if(nint(scale*10.**(-nexp)) .eq. 10.**(-nexp+1)) nexp=nexp+1
      if(dlog10(scale) .le. 0.)  nexp=nexp-1
      if(fixexp) nexp=nexfix
      factor=10.**(-nexp)
c
      if(k0.eq.k1) then
         modd=1
      elseif(j0.eq.j1) then
         modd=2
      elseif(i0.eq.i1) then
         modd=3
      else
         modd=mode
      endif
c
      if(modd.eq.0) then
         call outarr(a,mx0,mx1,i0,i1,ixss,ixp,my0,my1,j0,j1,iyss,iyp
     :      ,ms0,ms1,k0,k1,isss,isp,nplan,leng,title,nunit,scale,nexp
     :      ,factor,ndigit,tkoff,modd,nstep)
      else
         if(modd.eq.1) then
            nap=k1-k0+1
         elseif(modd.eq.2) then
            nap=j1-j0+1
         else
            nap=i1-i0+1
         endif
         call outarr(a,mx0,mx1,i0,i1,ixss,iap(i0)
     :      ,my0,my1,j0,j1,iyss,iap(j0)
     :      ,ms0,ms1,k0,k1,isss,iap(k0),nap,leng,title,nunit,scale,nexp
     :      ,factor,ndigit,tkoff,modd,nstep)
      endif
c
      return
1001  format(/t2,a8,t30,'min=',e14.7,' max=',e14.7)
      end



      subroutine wripro(fngrd)
c-----------------------------------------------------------------------
c write profiles
c-----------------------------------------------------------------------
c
      use alloc
      use refstate

      implicit real*8(a-h,o-z)
      character*14 timestring
      character*80 fextprf,fngrd
      integer ix,iy
c      include 'nh3dpa10.inc'


      save ifirst
      data ifirst/0/
c
      call timestamp(nstep,dt,starttime,timestring)
      if(verbose.ge.2) write(*,*) 'Profiles at:',timestring
      do ip=1,noutp
         ix=ixoutp(ip)
         iy=iyoutp(ip)
         open(32,file=fextprf(fngrd,'pf',ix,iy,timestring))
         write(32,2)
         do is=1,ns
            p=ptop+pp(ix,iy,2)*sigma0(is)
            ptem=pts(ix,iy,is)+pt(ix,iy,is,2)
            t=ptem*(p/p00)**akapa
            if(qif.eq.1.) then
             if(ifqc.eq.1) then
              if(ifqr.eq.1) then
               write(32,3) is,sigma0(is),(phis(ix,iy,is)
     :          +phi(ix,iy,is))/g
     :          ,p/100.,t,ptem,u(ix,iy,is,2),v(ix,iy,is,2)
     :          ,qv(ix,iy,is,2)*1000.,qv(ix,iy,is,2)/qsat(t,p)*100.
     :          ,w(ix,iy,is,2),qc(ix,iy,is,2)*1000.
     :          ,qr(ix,iy,is,2)*1000.,hfl(ix,iy,is),hflc(ix,iy,is)
     :          ,hfls(ix,iy,is)
              else
               write(32,3) is,sigma0(is),(phis(ix,iy,is)
     :          +phi(ix,iy,is))/g
     :          ,p/100.,t,ptem,u(ix,iy,is,2),v(ix,iy,is,2)
     :          ,qv(ix,iy,is,2)*1000.,qv(ix,iy,is,2)/qsat(t,p)*100.
     :          ,w(ix,iy,is,2),qc(ix,iy,is,2)*1000.,0.,
     :           hfl(ix,iy,is),
     :           hflc(ix,iy,is)
     :          ,hfls(ix,iy,is)
              endif
             else
              write(32,3) is,sigma0(is),(phis(ix,iy,is)
     :          +phi(ix,iy,is))/g
     :          ,p/100.,t,ptem,u(ix,iy,is,2),v(ix,iy,is,2)
     :          ,qv(ix,iy,is,2)*1000.,qv(ix,iy,is,2)/qsat(t,p)*100.
     :          ,w(ix,iy,is,2),0.,0.,hfl(ix,iy,is),hflc(ix,iy,is)
     :          ,hfls(ix,iy,is)
             endif
            else
              write(32,3) is,sigma0(is),(phis(ix,iy,is)
     :          +phi(ix,iy,is))/g
     :          ,p/100.,t,ptem,u(ix,iy,is,2),v(ix,iy,is,2),0.,0.
     :          ,w(ix,iy,is,2),0.,0.,hfl(ix,iy,is),hflc(ix,iy,is)
     :          ,hfls(ix,iy,is)
            endif
         enddo
         if(ifirst.eq.0) then
           write(31,'(a8,a14,22a8)')
     :       'nstep,','Time,','ix,','iy,','ts,','tlak','t2,'
     :       ,'wg,','w2,','wr,','rn,'
     :       ,'h,','le,','gs,','cd,','u,','v,','w,','th,','qv,','psuf,'
     :       ,'dz,'
           ifirst=1
         endif
         if(qif.eq.1..and.ifsoil.ne.0) then
           write(31,'(i6,'','',(a14,'',''),2(i8,'',''),3(f8.3,'','')
     :      ,3(f8.4,'',''),4(f8.1,'','')
     :      ,(f8.4,'',''),4(f8.2,'',''),(f8.4,'',''),(f8.0,'','')
     :      ,(f8.1,'',''))')
     :      nstep,timestring,ix,iy
     :      ,tsnoi(ix,iy),tslake(ix,iy)
     :      ,t2noi(ix,iy),wgnoi(ix,iy),w2noi(ix,iy)
     :      ,wrnoi(ix,iy),rn(ix,iy),h(ix,iy),le(ix,iy),gsolo(ix,iy)
     :      ,cdm(ix,iy),u(ix,iy,ns-1,2),v(ix,iy,ns-1,2),w(ix,iy,ns-1,2)
     :      ,pt(ix,iy,ns-1,2)+pts(ix,iy,ns-1),qv(ix,iy,ns-1,2)
     :      ,psuf(ix,iy)
     :      ,(phi(ix,iy,ns-1)+phis(ix,iy,ns-1))/g-hsuf(ix,iy)
         elseif(ifsoil.ne.0) then
           write(31,'(i6,'','',(a14,'',''),2(i8,'',''),3(f8.3,'','')
     :      ,3(f8.4,'',''),4(f8.1,'','')
     :      ,(f8.4,'',''),4(f8.2,'',''),(f8.4,'',''),(f8.0,'','')
     :      ,(f8.1,'',''))')
     :      nstep,timestring,ix,iy
     :      ,tsnoi(ix,iy),tslake(ix,iy)
     :      ,t2noi(ix,iy),wgnoi(ix,iy),w2noi(ix,iy)
     :      ,wrnoi(ix,iy),rn(ix,iy),h(ix,iy),le(ix,iy),gsolo(ix,iy)
     :      ,cdm(ix,iy),u(ix,iy,ns-1,2),v(ix,iy,ns-1,2),w(ix,iy,ns-1,2)
     :      ,pt(ix,iy,ns-1,2)+pts(ix,iy,ns-1),0.
     :      ,psuf(ix,iy)
     :      ,(phi(ix,iy,ns-1)+phis(ix,iy,ns-1))/g-hsuf(ix,iy)
         endif
      enddo
      close(32)
      return
2     format(t2,'is,',t7,'sigma,',t15,'z(m),',t25,'p(mb),',t36,'t(k),'
     :   ,t47,'theta(k),',t58,'u(m/s),',t69,'v(m/s),',t79,'q(g/kg),'
     :   ,t87,'rh,',t95,'w',t105,'qc',t115,'qr',t125,'Hg(w/m2),',t135,
     :   'Hc(w/m2),',t145,'Ht(w/m2),')
3     format(t2,i3,',',t7,f6.3,',',t15,f7.0,',',t25,f8.2,','
     :   ,t36,f8.2,',',t47,f8.2
     :   ,',',t58,f8.2,',',t69,f8.2,',',t79,f7.3,',',t87,f6.1
     :   ,t95,f9.5,t106,f9.4,t116,f9.6,',',t126,f9.2,',',t136,f9.2,','
     :   ,t146,f9.2)
      end


      subroutine wsigbc
c-----------------------------------------------------------------------
c boundary conditions for w and compute wsig.
c-----------------------------------------------------------------------
c
      use alloc
      implicit real*8(a-h,o-z)
c      include 'nh3dpa10.inc'
c
c solve diagnostic eq. for wsig from w-wsig relation
c
      do is=2,ns-1
      do iy=1,ny1
      do ix=1,nx1
        wsig(ix,iy,is,3)=
     :    -(w(ix,iy,is,3)*s(ix,iy,is)
     :    +sigma1(is)*(dpp(ix,iy)/pp(ix,iy,3)
     :    +0.25*(u(ix,iy,is,3)+u(ix,iy,is-1,3)+u(ix-1,iy,is,3)
     :    +u(ix-1,iy,is-1,3))*ppdx(ix,iy,3)
     :    +0.25*(v(ix,iy,is,3)+v(ix,iy,is-1,3)+v(ix,iy-1,is,3)
     :    +v(ix,iy-1,is-1,3))*ppdy(ix,iy,3)))
      enddo
      enddo
      enddo
      do iy=1,ny1
      do ix=1,nx1
        wsig(ix,iy,1,3)=0.0
        wsig(ix,iy,ns,3)=0.0
      enddo
      enddo
c
c compute wsig on lateral boundaries from continuity
c

      if(iowlbc.eq.0) then
        do ix=1,nx1,nx
c wsig(ix,iy,ns,3)=0.
        do is=ns-1,1,-1
        do iy=1,ny1
          wsig(ix,iy,is,3)=wsig(ix,iy,is+1,3)-(dpp(ix,iy)*ds1(is)
     :      +(u(ix,iy,is,3)*pp10(ix,iy,3)
     :      -u(ix-1,iy,is,3)*pp10(ix-1,iy,3))*ds1dx(is)
     :      +(v(ix,iy,is,3)*pp01(ix,iy,3)
     :      -v(ix,iy-1,is,3)*pp01(ix,iy-1,3))*ds1dy(is))/pp(ix,iy,3)
        enddo
        enddo
        enddo

        do is=ns-1,1,-1
        do iy=1,ny1
           wsig(1,iy,is,3)=wsig(2,iy,is,3)
           wsig(nx1,iy,is,3)=wsig(nx,iy,is,3)
        enddo
        enddo

c
c ensuring wsig(ix,iy,1,3)=0
c
        do ix=1,nx1,nx
        do iy=1,ny1
          if(wsig(ix,iy,1,3) .ne. 0.) then
            do is=2,ns-1
              wsig(ix,iy,is,3)=wsig(ix,iy,is,3)-(ns-is)/(ns-1.0)
     :          *wsig(ix,iy,1,3)
            enddo
            wsig(ix,iy,1,3)=0.0
          endif
        enddo
        enddo
c
        do iy=1,ny1,ny
c wsig(ix,iy,ns,3)=0.
        do is=ns-1,1,-1
        do ix=1,nx1
          wsig(ix,iy,is,3)=wsig(ix,iy,is+1,3)-(dpp(ix,iy)*ds1(is)
     :      +(u(ix,iy,is,3)*pp10(ix,iy,3)
     :      -u(ix-1,iy,is,3)*pp10(ix-1,iy,3))*ds1dx(is)
     :      +(v(ix,iy,is,3)*pp01(ix,iy,3)
     :      -v(ix,iy-1,is,3)*pp01(ix,iy-1,3))*ds1dy(is))/pp(ix,iy,3)
        enddo
        enddo
        enddo

        do is=ns-1,1,-1
        do ix=1,nx1
           wsig(ix,1,is,3)=wsig(ix,2,is,3)
           wsig(ix,ny1,is,3)=wsig(ix,ny,is,3)
        enddo
        enddo

        
c
c ensuring wsig(ix,iy,1,3)=0
c
        do iy=1,ny1,ny
        do ix=1,nx1
          if(wsig(ix,iy,1,3) .ne. 0.) then
            do is=2,ns-1
              wsig(ix,iy,is,3)=wsig(ix,iy,is,3)-(ns-is)/(ns-1.0)
     :          *wsig(ix,iy,1,3)
            enddo
            wsig(ix,iy,1,3)=0.0
          endif
        enddo
        enddo
      endif
c
c      do 703 iy=1,ny1
c      do 703 ix=1,nx1
c      wsig(ix,iy,ns1,3)=-ds1(ns+1)*(dpp(ix,iy)
c     :   +(uu(ix,iy,ns,3)-uu(ix-1,iy,ns,3))/dx
c     :   +(vv(ix,iy,ns,3)-vv(ix,iy-1,ns,3))/dy)/pp(ix,iy,3)
c      wsig(ix,iy,0,3)=+ds1(0)*(dpp(ix,iy)
c     :   +(uu(ix,iy,0,3)-uu(ix-1,iy,0,3))/dx
c     :   +(vv(ix,iy,0,3)-vv(ix,iy-1,0,3))/dy)/pp(ix,iy,3)
c703   continue
c
      do iy=1,ny1
      do ix=1,nx1
        wsig(ix,iy,0,3)=-wsig(ix,iy,2,3)
        wsig(ix,iy,ns1,3)=-wsig(ix,iy,ns-1,3)
c      wsig(ix,iy,0,3)=0.
c      wsig(ix,iy,ns1,3)=0.
      enddo
      enddo
c
c compute w on boundaries from wsig
c
      if(iowlbc.eq.0) then
        do is=2,ns-1
          do ix=1,nx1,nx
          do iy=1,ny1
            w(ix,iy,is,3)=-(wsig(ix,iy,is,3)
     :        +sigma1(is)*(dpp(ix,iy)/pp(ix,iy,3)
     :        +0.25*(u(ix,iy,is,3)+u(ix,iy,is-1,3)+u(ix-1,iy,is,3)
     :        +u(ix-1,iy,is-1,3))*ppdx(ix,iy,3)
     :        +0.25*(v(ix,iy,is,3)+v(ix,iy,is-1,3)+v(ix,iy-1,is,3)
     :        +v(ix,iy-1,is-1,3))*ppdy(ix,iy,3)))/s(ix,iy,is)
          enddo
          enddo
          do iy=1,ny1,ny
          do ix=1,nx1
            w(ix,iy,is,3)=-(wsig(ix,iy,is,3)
     :         +sigma1(is)*(dpp(ix,iy)/pp(ix,iy,3)
     :         +0.25*(u(ix,iy,is,3)+u(ix,iy,is-1,3)+u(ix-1,iy,is,3)
     :         +u(ix-1,iy,is-1,3))*ppdx(ix,iy,3)
     :         +0.25*(v(ix,iy,is,3)+v(ix,iy,is-1,3)+v(ix,iy-1,is,3)
     :         +v(ix,iy-1,is-1,3))*ppdy(ix,iy,3)))/s(ix,iy,is)
          enddo
          enddo
        enddo
      endif

      do is=ns-1,1,-1
      do iy=1,ny1
           w(1,iy,is,3)=w(2,iy,is,3)
           w(nx1,iy,is,3)=w(nx,iy,is,3)
      enddo
      enddo
        do is=ns-1,1,-1
        do ix=1,nx1
           w(ix,1,is,3)=w(ix,2,is,3)
           w(ix,ny1,is,3)=w(ix,ny,is,3)
        enddo
        enddo

c
c   compute w on and outside the top and bottom boudaries
c   to implement  mirror bundary conditions
c
      do iy=1,ny1
      do ix=1,nx1
c
c note wsig=0 and sigma=0 then
c
        w(ix,iy,1,3)=0.0
c
c note: i am making tems(ix,iy,-1)=tems(ix,iy,0)
c
c      stop=gr*(ptop+sigma1(0)*pp(ix,iy,3))/(tems(ix,iy,0)*pp(ix,iy,3))
        w(ix,iy,0,3)=(-wsig(ix,iy,0,3)-sigma1(0)*(dpp(ix,iy)
     :    /pp(ix,iy,3)+(u(ix-1,iy,1,3)+u(ix,iy,1,3))*0.5*ppdx(ix,iy,3)
     :    +(v(ix,iy-1,1,3)+v(ix,iy,1,3))*0.5*ppdy(ix,iy,3)))/s(ix,iy,0)
c
c      w (ix,iy,0,3)=-w(ix,iy,2,3)
      enddo
      enddo
c
      do iy=1,ny1
      do ix=1,nx1
        w(ix,iy,ns,3)=-(dpp(ix,iy)/pp(ix,iy,3)
     :    +(u(ix,iy,ns,3)+u(ix,iy,ns-1,3)+u(ix-1,iy,ns-1,3)
     :    +u(ix-1,iy,ns,3))*0.25*ppdx(ix,iy,3)
     :    +(v(ix,iy,ns,3)+v(ix,iy,ns-1,3)+v(ix,iy-1,ns-1,3)
     :    +v(ix,iy-1,ns,3))*0.25*ppdy(ix,iy,3))/s(ix,iy,ns)
      enddo
      enddo
c
      do iy=1,ny1
      do ix=1,nx1
c      sbot=gr*(ptop+sigma1(ns1)*pp(ix,iy,3))/(tems(ix,iy,ns)
c     :   *pp(ix,iy,3))
        w(ix,iy,ns1,3)=(-wsig(ix,iy,ns1,3)
     :    -sigma1(ns1)*(dpp(ix,iy)/pp(ix,iy,3)
     :    +(u(ix,iy,ns,3)+u(ix-1,iy,ns,3))*0.5*ppdx(ix,iy,3)
     :    +(v(ix,iy,ns,3)+v(ix,iy-1,ns,3))*0.5*ppdy(ix,iy,3)))
     :    /s(ix,iy,ns+1)
      enddo
      enddo
c
      if(raylei .and. iowlbc.eq.0) then
         do is=1,idrmax
         do iy=1,ny1,ny
         do ix=1,nx1
           w(ix,iy,is,3)=w(ix,iy,is,3)*(1.-dtl/taudra(is))
         enddo
         enddo
         enddo
         do is=1,idrmax
         do ix=1,nx1,nx
         do iy=1,ny1
           w(ix,iy,is,3)=w(ix,iy,is,3)*(1.-dtl/taudra(is))
         enddo
         enddo
         enddo
      endif
c
      return
      end






      function chstar(xmu)

      implicit real*8(a-h,o-z)
      chstar=3.2165+4.3431*xmu+.536*xmu*xmu-.0781*xmu*xmu*xmu
      return
      end



      function cmstar(xmu)

      implicit real*8(a-h,o-z)
      cmstar=6.8741+2.6933*xmu+.3601*xmu*xmu-.0154*xmu*xmu*xmu
      return
      end



      function dqsat (t,p)

      implicit double precision (a-h,o-z)
c     write(*,*) 'Qsat:',t,p
      rasrv = 287.05/461.51
      etvq = 1. - rasrv
      esat=618.78*exp(17.269*(t-273.16)/(t-35.86))
      desat=4097.93*esat/((t-35.86)*(t-35.86))
      dqsat=rasrv*desat*p/((p-etvq*esat)*(p-etvq*esat))
      end



      function esat (t)
c        tensao de saturacao
      implicit double precision (a-h,o-z)
      esat=618.78*exp(17.269*(t-273.16)/(t-35.86))
	end

      function esati (t)
c        tensao de saturacao with respect to ice     DC, 10.2009
      implicit double precision (a-h,o-z)
      esati=611.2*10.**(0.43*22.46*(t-273.15)/t)
      end


      function fext(fname,ext)
      character*80 fname,fext
      character*3 ext
      parameter(n=80)
c
c elimina extensao se existir
c
      fext=fname




      do 10 i=n,1,-1
         if(fext(i:i).eq.'.') then
            fext(i:n)=char(0)
            go to 11
         endif
10    continue
11    continue

      do 20 i=n,1,-1
         if(fext(i:i).ne.' ' .and. fext(i:i).ne.char(0)) then
            fext(i+1:i+4)='.'//ext
            go to 21
         endif
20    continue
      write(*,*) 'error in fext'

21    continue
      return
      end

      function fileext(fname,ext)
      character*80 fname,fileext
      character*80 ext
!      parameter(n=80,next=80)                      
      parameter(n=80)                                                           VS,05.2007
      integer next                                                              VS,05.2007

      fileext=fname
      next = len_trim(ext)                                                      VS,05.2007

c appends ext to fname

      do i=n,1,-1
        if(fileext(i:i).ne.' ' .and. fileext(i:i).ne.char(0)) then
          iright=max(n,i+1+next)
          fileext(i+1:iright)='.'//ext
          exit
        endif
      enddo
      return
      end






c******************************************************************
       function ph(xmu)

       implicit real*8(a-h,o-z)
       ph    =0.5802-0.1571*xmu+.0327*xmu*xmu-.0026*xmu*xmu*xmu
       return
       end
c
c******************************************************************
        function pm(xmu)
c****************************************************************
c
       implicit real*8(a-h,o-z)
       pm    =0.5233-0.0815*xmu+.0135*xmu*xmu-.001*xmu*xmu*xmu
       return
       end



       function psih(xsi)
        implicit real*8(a-h,o-z)

        if(xsi) 10,10,20
c       cas instable
10      x=(1.-16.*xsi)**.25
        a=(1.+x*x)/2.
        psih=2.*dlog(a)
        return
c       cas stable
20      continue
        if (xsi.gt.1.) then
                psih=-5.*(1.+dlog(xsi))
                else
                psih=-5.*xsi
        endif
        return
        end


c*********************************************************
c                                                        *
c       Fonctions de stabilite psim et psih d'apres      *
c       les formules de Paulson                          *
c*********************************************************

        function psim(xsi)
        implicit real*8(a-h,o-z)

        if(xsi) 10,10,20
c       cas instable
10      x=(1.-16.*xsi)**.25
        a=(1.+x)/2.
        b=(1.+x*x)/2.
        psim=2.*dlog(a)+dlog(b)-2.*atan(x)+asin(1.)
        return
c       cas stable
20      continue
        if (xsi.gt.1.) then
                psim=-5.*(1.+dlog(xsi))
                else
                psim=-5.*xsi
        endif
        return
        end
c
c**********************************************************************
      function qsat (t,p)
c          humidade especifica de saturacao
c
      implicit double precision (a-h,o-z)
      rasrv = 287.05/461.51
      etvq = 1. - rasrv
      espf=esat(t)/p
      qsat=rasrv*espf/(1.-etvq*espf)
	!qsat=380./p*exp(17.27*(t-273.15)/(t-35.86))
      end
  
	function qsati (t,p)                                                       DC,10.2009
c          humidade especifica de saturacao with respect to ice
c
      implicit double precision (a-h,o-z)
      rasrv = 287.05/461.51
      etvq = 1. - rasrv
      espf=esati(t)/p
      qsati=rasrv*espf/(1.-etvq*espf)
	!qsati=380./p*exp(21.88*(t-273.15)/(t-7.66))
      end

      function lvap(t)                                                            DC,10.2009
	implicit double precision (a-h,o-z)
	alv=2.4*10**3
	if (t.lt.253.15) blv=2631.5-10.*t
	if (t.ge.253.15) blv=0.
	lvap=2.549*10**6+(alv+blv)*(253.15-t)
	end

	
      function tlake(ix,iy,timeho)
      use alloc
      implicit real*8(a-h,o-z)
c      include 'nh3dpa10.inc'

      pi180=4.*datan(1.d0)/180.

      tlake=tsmed(ix,iy)

      return
      end



      function tsoil(ix,iy,timeho)
      use alloc
      implicit real*8(a-h,o-z)
c      include 'nh3dpa10.inc'

      pi180=4.*datan(1.d0)/180.

c     tsoil=tsmed(ix,iy)

        tsoil=462.86962-62.3327*timeho+8.04565*(timeho**2)-
     :      0.4015*(timeho**3)+0.00678*(timeho**4)

	pi=4.*atan(1.)                                                            VS,06.2007
	A=10.                                                                     VS,06.2007
	tsoil=tsmed(ix,iy)+A*sin(2*pi*timeho/24.-pi/2.)                           VS,06.2007

      return
      end


      subroutine wgrids(var,varname,i0,i1,j0,j1,k0,k1,mode
     :  ,ixgrid,iygrid,isgrid,fngrd)

      use alloc
      use refstate
      use nh3dpa2
      implicit real*8(a-h,o-z)
      dimension var(0:nx1,0:ny1,0:*)
      character*2 varname
      character*80 fngrd
      character*20 line
      character*14 timestring
      integer i0,i1,j0,j1,k0,k1,mode,ixgrid,iygrid,isgrid

      if(kstep.ne.nstep) then
        kstep=nstep
        istep=istep+1
      endif

      if(verbose.ge.2) then
        write(*,'(a8,1x,a2,12i4)')
     :   ' wgrids:',varname,i0,i1,j0,j1,k0,k1,mode,ixgrid,iygrid,isgrid
      endif

c      allocate(zed(0:nx1,0:ny1,0:ns1))
      maxh=max(nx+1,ny+1)
      maxv=max(ny+1,ns+1,nlev)
      allocate(hkk1(0:maxh,0:maxv))
      allocate(hkk2(0:maxh,0:maxv))
      allocate(hkk3(0:maxh,0:maxv))
      allocate(ihkk1(0:nx1,0:ny1))
      allocate(ihkk2(0:maxh,0:maxv))

      if(iomtyp.eq.4) then
        xmountc=-xminutm-1.5*dx
        ymountc=-yminutm-1.5*dx
      else
        xmountc=xmount
        ymountc=ymount
      endif

      call timestamp(nstep,dt,starttime,timestring)

      call wri3grd(var,varname,nstep,i0,i1,j0,j1,k0,k1,mode
     :  ,ixgrid,iygrid,isgrid
     :  ,nx,ny,ns,nstep,fngrd,timestring
     :  ,xmountc,ymountc,dx,dy,phi,phis,ixp,iyp,isp,zplan,nplan
     :  ,sigma0,sigma1,zconst,zlev,nlev,hsuf
     :  ,wk1,hkk3,hkk3,hkk1,hkk1,hkk1,hkk2,hkk2,ihkk1,ihkk2,ihkk2)
c     :  ,zed,hkk3,hkk3,hkk1,hkk1,hkk1,hkk2,hkk2,ihkk1,ihkk2,ihkk2)

      deallocate(hkk1,hkk2,hkk3,ihkk1,ihkk2)

      return
      end

      subroutine timestamp(nstep,dt,starttime,timestring)

      implicit none
      character*20 line
      character*14 timestring
      integer styear,stmonth,stday,sthour,stmin,stsec
      integer cuyear,cumonth,cuday,cuhour,cumin,cusec
      integer elapseddays,elapsedsec,verbose,nstep
      real*8 dt,starttime

      write(line,'(f15.0)') starttime
      read(line,'(i4,5i2)') styear,stmonth,stday,sthour,stmin,stsec
      elapsedsec=nstep*dt

!     write(*,'(''timestamp:'',f15.0,6i5)')
!    :  starttime,styear,stmonth,stday,sthour,stmin,stsec

      call timeofyear(styear,stmonth,stday,sthour,stmin,stsec
     :  ,elapsedsec,cuyear,cumonth,cuday,cuhour,cumin,cusec)

      write(timestring,'(i4,5i2)')
     :  cuyear,cumonth,cuday,cuhour,cumin,cusec

      return
      end





      function fextgrd(fngrd,varname,cmode,iplan,istep)
c
c Automatic creation of file names for grd files
c Rule: fngrd+varname+mode+iplan.istep
c       (Ncar)+(2car)+(1car)+(2int)+(3int)
c Limitations: 99 planes per mode,999 different time steps
c size of fngrd depends on file system
c
      character*80 fngrd,fextgrd
      integer iplan,istep
      character*2 varname
      character*1 cmode
      character*13 suff
      parameter(n=80)
c
c clears file extension:
c
      fextgrd=fngrd

      do i=n,1,-1
         if(fextgrd(i:i).eq.'.') then
            fextgrd(i:n)=char(0)
            exit
         endif
      enddo

      write(suff,'(a2,a1,i3,''.'',i6)') varname,cmode,iplan,istep
      do i=1,13
         if(suff(i:i).eq.' '.or.suff(i:i).eq.char(0)) suff(i:i)='0'
      enddo

      do i=n,1,-1
         if(fextgrd(i:i).ne.' ' .and. fextgrd(i:i).ne.char(0)) then
            fextgrd(i+1:i+13)=suff
            exit
         endif
      enddo

      return
      end
      function fextgrdt(fngrd,varname,cmode,iplan,timestring)
c
c Automatic creation of file names for grd files
c Rule: fngrd+varname+mode+iplan.timestring
c       (Ncar)+(2car)+(1car)+(2int)+(3int)
c Limitations: 99 planes per mode,999 different time steps
c size of fngrd depends on file system
c
      character*80 fngrd,fextgrdt
      character*14 timestring
      integer iplan,istep
      character*2 varname
      character*1 cmode
      character*21 suff
      parameter(n=80)
c
c clears file extension:
c
      fextgrdt=fngrd

      do i=n,1,-1
         if(fextgrdt(i:i).eq.'.') then
            fextgrdt(i:n)=char(0)
            exit
         endif
      enddo

      write(suff,'(a2,a1,i3,''.'',a14)') varname,cmode,iplan,timestring
      do i=1,21
         if(suff(i:i).eq.' '.or.suff(i:i).eq.char(0)) suff(i:i)='0'
      enddo

      do i=n,1,-1
         if(fextgrdt(i:i).ne.' ' .and. fextgrdt(i:i).ne.char(0)) then
            fextgrdt(i+1:i+21)=suff
            exit
         endif
      enddo

      return
      end

      function fextprf(fngrd,proname,ix,iy,timestring)
c
c Automatic creation of file names for profiles
c Rule: fngrd+proname+ix+iy.timestring
c       (Ncar)+(2car)+(3int)+(3int)+(14a)
c size of fngrd depends on file system
c
      character*80 fngrd,fextprf
      character*14 timestring
      integer ix,iy
      character*2 proname
      character*25 suff
      parameter(n=80)
c
c clears file extension:
c
      fextprf=fngrd

      do i=n,1,-1
         if(fextprf(i:i).eq.'.') then
            fextprf(i:n)=char(0)
            exit
         endif
      enddo

      write(suff,'(a2,a1,i3,a1,i3,''.'',a14)')
     :  proname,'x',ix,'y',iy,timestring
      do i=1,25
         if(suff(i:i).eq.' '.or.suff(i:i).eq.char(0)) suff(i:i)='0'
      enddo

      do i=n,1,-1
         if(fextprf(i:i).ne.' ' .and. fextprf(i:i).ne.char(0)) then
            fextprf(i+1:i+25)=suff
            exit
         endif
      enddo

      return
      end

      function fextbln(fngrd,varname,cmode,iplan)
c
c Automatic creation of file names for bln files
c Rule: fngrd+varname+mode+iplan.bln
c       (Ncar)+(2car)+(1car)+(2int)+(3int)
c Limitations: 99 planes per mode,999 different time steps
c size of fngrd depends on file system
c
      character*80 fngrd,fextbln
      integer iplan,istep
      character*2 varname
      character*1 cmode
      character*10 suff
      parameter(n=80)

      nsuff=10
c
c clears file extension:
c
      fextbln=fngrd

      do i=n,1,-1
         if(fextbln(i:i).eq.'.') then
            fextbln(i:n)=char(0)
            exit
         endif
      enddo

      write(suff,'(a2,a1,i3,''.bln'')') varname,cmode,iplan
      do i=1,nsuff
         if(suff(i:i).eq.' '.or.suff(i:i).eq.char(0)) suff(i:i)='0'
      enddo

      do i=n,1,-1
         if(fextbln(i:i).ne.' ' .and. fextbln(i:i).ne.char(0)) then
            fextbln(i+1:i+nsuff)=suff
            exit
         endif
      enddo

      return
      end










      function fextvi(fname,var,iext)
      character*80 fname,fextvi
      integer iext
      character*3 ext,var
      parameter(n=80)
c
c elimina extensao se existir
c
      fextvi=fname

      do i=n,1,-1
         if(fextvi(i:i).ne.' '.and.fextvi(i:i).ne.char(0)) then
            in=i+1
            exit
         endif
      enddo

      fextvi(in:in)='_'
      fextvi(in+1:in+3)=var

      write(ext,'(i3)') iext
      do ii=1,3
         if(ext(ii:ii).eq.' '.or. ext(ii:ii).eq.char(0)) ext(ii:ii)='0'
      enddo

      do 20 i=n,1,-1
         if(fextvi(i:i).ne.' ' .and. fextvi(i:i).ne.char(0)) then
            fextvi(i+1:i+4)='.'//ext
            go to 21
         endif
20    continue
      write(*,*) 'error in fextvi'
21    continue
      return
      end


      subroutine wrivec(vx,vy,vecname,mode,fngrd,i0,i1,j0,j1)

      use alloc
      use nh3dpa2
      implicit real*8(a-h,o-z)
      dimension vx(0:nx+1,0:ny+1,0:ns+1),vy(0:nx+1,0:ny+1,0:ns+1)
      character*80 fngrd,fextgrd
      character*2 vecname
      integer mode,i0,i1,j0,j1

      if(kstep.ne.nstep) then
        kstep=nstep
        istep=istep+1
      endif

      if(verbose.ge.2) then
        write(*,'(a7,a3,4i4,2i6)') ' wrivec',vecname,i0,i1,j0,j1,nstep
     :    ,istep
      endif
c
c write horizontal plane (k=cst):
c
      if(mode.eq.1) then
        do kp=1,nplan
          if(isp(kp).ge.0) then
             k=isp(kp)
             open(77,file=fextgrd(fngrd,vecname,'s',k,nstep)
     :         ,status='unknown')
             do iy=j0,j1
               do ix=i0,i1
                 xvec=(ix-1.5)*dx-xmount
                 yvec=(iy-1.5)*dy-ymount
                 uvec=0.5*(vx(ix,iy,k)+vx(ix-1,iy,k))
                 vvec=0.5*(vy(ix,iy,k)+vy(ix,iy-1,k))
                 call pstvec(xvec,yvec,uvec,vvec,77)
               enddo
             enddo
             close(77)
          endif
        enddo
        if(zplan(1).ge.0.) then
c          allocate(zed(0:nx1,0:ny1,0:ns1))
          allocate(hkk1(0:nx1,0:ny1))
          allocate(hkk2(0:nx1,0:ny1))
          allocate(ihkk1(0:nx1,0:ny1))
c          zed=(phi+phis)/g
          wk1=(phi+phis)/g
          do kp=1,nplan
            if(zplan(kp).ge.0.) then
              call inter3(i0-1,i1,j0,j1,ihkk1,wk1,zplan(kp),nx+1,ny+1,ns
     :          ,vx,hkk1)
              call inter3(i0,i1,j0-1,j1,ihkk1,wk1,zplan(kp),nx+1,ny+1,ns
     :          ,vy,hkk2)
              iplan=nint(zplan(kp)/100.)
              open(77,file=fextgrd(fngrd,vecname,'z',iplan,nstep)
     :          ,status='unknown')
              do iy=j0,j1
                do ix=i0,i1
                  xvec=(ix-1.5)*dx-xmount
                  yvec=(iy-1.5)*dy-ymount
                  uvec=0.5*(hkk1(ix,iy)+hkk1(ix-1,iy))
                  vvec=0.5*(hkk2(ix,iy)+hkk2(ix,iy-1))
                  call pstvec(xvec,yvec,uvec,vvec,77)
                enddo
              enddo
              close(77)
            endif
          enddo
          deallocate(hkk1,hkk2,ihkk1)
        endif
      endif
c
c write vertical plane (j=cst):
c
      if(mode.eq.2.and.iyp(1).ge.0) then
        allocate(hkk1(0:nx1,0:ns1))
        allocate(hkk2(0:nx1,0:ns1))
        allocate(hkk3(0:nx1,0:nlev))
        allocate(hkk4(0:nx1,0:nlev))
        allocate(hkk5(0:nx1,0:ns1))
        allocate(ihkk1(0:nx1,0:nlev))
        do jp=1,nplan
          if(iyp(jp).ge.0) then
            j=iyp(jp)
            do is=0,ns+1
              do ix=0,nx+1
                hkk1(ix,is)=vx(ix,j,is)
              enddo
            enddo
            do is=0,ns+1
              do ix=0,nx+1
                hkk2(ix,is)=vy(ix,j,is)
              enddo
            enddo
            do ix=0,nx+1
              zbot(ix)=hsuf(ix,j)
            enddo
            if(zconst) then
              do is=0,ns+1
              do ix=0,nx+1
                hkk5(ix,is)=(phi(ix,j,is)+phis(ix,j,is))/g
              enddo
              enddo
              call inter2i(nx+1,ns+1,i0-1,i1,j0,j1,ihkk1,hkk5,zlev,nlev
     :          ,hkk1,hkk3,zbot)
              do is=1,ns-1
                hkk5(:,is)=0.5*(hkk5(:,is)+hkk5(:,is-1))
              enddo
              call inter2i(nx+1,ns+1,i0,i1,j0,j1,ihkk1,hkk5,zlev,nlev
     :          ,hkk2,hkk4,zbot)
              open(77,file=fextgrd(fngrd,vecname,'y',j,nstep)
     :          ,status='unknown')
              do iz=0,nlev
                do ix=i0,i1
                  xvec=(ix-1.5)*dx-xmount
                  yvec=iz*dzlev
                  uvec=0.5*(hkk3(ix,iz)+hkk3(ix-1,iz))
                  vvec=hkk4(ix,iz)*scalew
                  call pstvec(xvec,yvec,uvec,vvec,77)
                enddo
              enddo
              close(77)
            endif
          endif
        enddo
        deallocate(hkk1,hkk2,hkk3,hkk4,hkk5,ihkk1)
      endif
c
c write vertical plane (i=cst):
c
      if(mode.eq.3.and.ixp(1).ge.0) then
        allocate(hkk1(0:ny1,0:ns1))
        allocate(hkk2(0:ny1,0:ns1))
        allocate(hkk3(0:ny1,0:nlev))
        allocate(hkk4(0:ny1,0:nlev))
        allocate(hkk5(0:ny1,0:ns1))
        allocate(ihkk1(0:ny1,0:nlev))
c        write(*,*) 'mode3'
        do ip=1,nplan
          if(ixp(ip).ge.0) then
            i=ixp(ip)
c        write(*,*) 'i=',i
            do is=0,ns+1
              do iy=0,ny+1
                hkk1(iy,is)=vx(i,iy,is)
              enddo
            enddo
c        write(*,*) 'hkk1'
            do is=0,ns+1
              do iy=0,ny+1
                hkk2(iy,is)=vy(i,iy,is)
              enddo
            enddo
c        write(*,*) 'hkk2'
            do iy=0,ny+1
              zbot(iy)=hsuf(i,iy)
            enddo
c        write(*,*) 'zbot'
            if(zconst) then
c        write(*,*) 'zconst'
              do is=0,ns+1
              do iy=0,ny+1
                hkk5(iy,is)=(phi(i,iy,is)+phis(i,iy,is))/g
              enddo
              enddo
c        write(*,*) 'hkk5'
              call inter2i(ny+1,ns+1,i0-1,i1,j0,j1,ihkk1,hkk5,zlev,nlev
     :          ,hkk1,hkk3,zbot)
c        write(*,*) 'inter2i'
              do is=1,ns-1
                do iy=0,ny+1
                  hkk5(iy,is)=0.5*(hkk5(iy,is)+hkk5(iy,is-1))
                enddo
              enddo
c        write(*,*) 'hkk5'
              call inter2i(ny+1,ns+1,i0,i1,j0,j1,ihkk1,hkk5,zlev,nlev
     :          ,hkk2,hkk4,zbot)
c        write(*,*) 'inter2i'
              open(77,file=fextgrd(fngrd,vecname,'x',j,nstep)
     :          ,status='unknown')
              do iz=0,nlev
                do iy=i0,i1
c                  write(*,*) iz,iy
                  xvec=(iy-1.5)*dy-ymount
                  yvec=iz*dzlev
                  uvec=0.5*(hkk3(iy,iz)+hkk3(iy-1,iz))
                  vvec=hkk4(iy,iz)*scalew
                  call pstvec(xvec,yvec,uvec,vvec,77)
                enddo
              enddo
              close(77)
            endif
          endif
        enddo
        deallocate(hkk1,hkk2,hkk3,hkk4,hkk5,ihkk1)
      endif
c
      return
      end


      subroutine wribln(hsuf,nx,ny,dx,dy,xmount,ymount,fngrd)
      use nh3dpa2
      implicit real*8 (a-h,o-z)
      real*8 hsuf(0:nx+1,0:ny+1)
      character*80 fextbln,fngrd
c
c Output mountain line for cross sections
c
      do ip=1,nplan
        if(ixp(ip).ge.0) then
          i=ixp(ip)
          do iy=0,ny+1
            zbot(iy)=hsuf(i,iy)
          enddo
          call writbln(fextbln(fngrd,'h_','x',i),zbot,ny
     :      ,dy,ymount)
        endif
        if(iyp(ip).ge.0) then
          j=iyp(ip)
          do ix=0,nx+1
            zbot(ix)=hsuf(ix,j)
          enddo
          call writbln(fextbln(fngrd,'h_','y',j),zbot,nx
     :      ,dx,xmount)
        endif
      enddo
      return
      end

C
      subroutine inpos3c
c
c revised version of ipos3, pos3:
c employs inver (eigenvectors in vertical)
c-----------------------------------------------------------------------
c
c note that the dimensions of dphi do not coincide with those in the
c subroutines dscfft and tri3
c as it is assumed that the physical domain starts ,there, at
c dphi(0.5x,0.5y) whereas it is dphi(1.5x,1.5y)
c
c-----------------------------------------------------------------------
c
c in this version eigensystem is defined in the domain is = (0,ns-1),
c i.e., it  includes level is=0;
c due to use of tqli, emo, emi, and vlam are saved in shifted
c domain (1,..,ns)
c
c ----------------------------------------------------------------------
      use alloc
      implicit real*8(a-h,o-z)
c
      dimension bi(0:nx1,0:ny1)
      dimension dphi(0:nx1,0:ny1,0:ns1)
      dimension ww1(0:nx1,0:ny1,0:ns1),ww2(0:nx1,0:ny1,0:ns1)
c moved to alloc
c      dimension emo(ns,ns),emi(ns,ns),vlam(ns),e(ns+1)
c
c      save emo, emi, vlam
      save ppmean
c
c initialization of steady constants of elliptic operator
c ======================================================
c
c     write(*,*) 'inpos3c'
      do is=0,ns
        aves(is)=0.0
        do iy=1,ny1
        do ix=1,nx1
          aves(is)=aves(is)+s(ix,iy,is)*s(ix,iy,is)
        enddo
        enddo
        aves(is)=aves(is)/(nx1*ny1*ds0(is))
      enddo
c     write(*,*) 'inpos3c cont 10'
c
C*  Preparation of matrix elemnts
C
      do i=1,ns
      do j=1,ns
        emo(i,j)=0.0
      enddo
      enddo
c     write(*,*) 'inpos3c+1'

      do i=0,ns-1
        emo(i+1,i+1)=1.
        vlam(i+1)=-(aves(i)+aves(i+1))/ds1(i)
        e(i+2)=aves(i+1)/sqrt(ds1(i+1)*ds1(i))
      enddo
c     write(*,*) 'inpos3c+2'
C       e(1)=0.
C
C*  Solving eigenvalue problem
C
      call tqli(vlam,e,ns,ns,emo)
c     write(*,*) 'inpos3c+3'
C
C*  Inverse matrix emi = (emo)^transposed
C   And renormation of emo and emi
C
      do i=1,ns
      do j=1,ns
       emi(i,j)=emo(j,i)*sqrt(ds1(j-1))
       emo(j,i)=emo(j,i)/sqrt(ds1(j-1))
      enddo
      enddo
c     write(*,*) 'inpos3c+4'
C
C*    PRINT EIGEN-VALUES AND EIGEN-VECTORS
C     ------------------------------------
C
C      WRITE(*,'(/,1X,''EIGVAL='',10(6E10.3,/,8X))')(vlam(J),J=1,ns)
C
C      DO  K=1,klev
C      WRITE(6,'(/,1X,''EIGVEC='',10(6E10.3,/,8X))')(EMO(J,K),J=1,ns)
C      ENDDO
C
C      call plotgrid("EIGV",EMO,ns,ns)
C      call plotgrid("VLAM",vlam,ns,1)
C      stop
C
c
c Mean effective surface pressure ppmean
c
      ppmean=0.
      do ix=1,nx1
      do iy=1,ny1
       ppmean=ppmean+pp(ix,iy,3)
      enddo
      enddo
      ppmean=ppmean/nx1/ny1
c     write(*,*) 'inpos3c+5'
c
      return
c
c solve poisson equation:
c
      entry pos3c(bi,dphi,ww1,ww2,ite)
c
c perform double cosine transform of of bi and r.h.s of
c the poisson eq.
c
      call dscfft1(bi,ww1,ww2,ww2,trigsx,ifax,sipcox,simcox,nx-1
     :   ,trigsy,ifay,sipcoy,simcoy,ny-1,0,+1)
c
      call dscfft1(dphi,ww1,ww2,ww2,trigsx,ifax,sipcox,simcox,nx-1
     :   ,trigsy,ifay,sipcoy,simcoy,ny-1,ns-1,+1)
c
c on exit dphi contains the transformed r.h.s. (in spectral space).
c
c solve Poisson equation in vertical with the help of eigenvectors
c
      call inver(nx-1,ny-1,ns-1
     R  , dphi,xylamb,ds1
     R  , bi,ppmean,emo,emi,vlam)
c
c transform solution back to physical space
c
      call dscfft1(dphi,ww1,ww2,ww2,trigsx,ifax,sipcox,simcox,nx-1
     :   ,trigsy,ifay,sipcoy,simcoy,ny-1,ns-1,-1)
c
c Meet boundary condition requirement.
c
      call extrah(nx1,ny1,dphi,1,nx1,1,ny1,0,ns)
c
      do is=0,ns
        do iy=0,ny1
          dphi(0,iy,is)=0.
        enddo
        do ix=0,nx1
          dphi(ix,0,is)=0.
        enddo
      enddo
c
      return
      end


c *** inver ***
C inversion in F. space, making use of eigenvectors,
C with mass balance in the role of BC
C
C dimension of eigensystem is enhanced to nsxns to include phi at level is=0
C
      SUBROUTINE inver(klon, klat, klev
     R , skf,xylamb,ds1
     R , b,pmean,EMO,EMI,vlam )
C
C     Input: skf: source function, in F space
C            b: right hand side for mass balanse condition, in F space
c     output: skf: solution phi, in F space
c
C*    Declaration of global parameters
C     --------------------------------
C
      IMPLICIT real*8(a-h,o-z)
C
C
      dimension  skf(-1:klon+1,-1:klat+1,0:klev+2)
     R , ds1(0:klev+2)
     R , xylamb(klon,klat)
     R , b(-1:klon+1,-1:klat+1)
     R , EMO(0:klev,0:klev),EMI(0:klev,0:klev),vlam(0:klev)
C
C*    Declaration of local workspace
C     ------------------------------
C
      dimension  wr1(0:klev), whom(0:klev), wnon(0:klev)
C
C    =============================================
C
      do j=1,klat
      do i=1,klon
       if(i*j.ne.1)
     &     b(i,j)= b(i,j)/xylamb(i,j)/pmean
      enddo
      enddo
C      b(1,1) = 0.0 ! no effect (short check with 5 time-steps)
C
C*    Elliptic solver in Fourier space
C     --------------------------------
C
      do j=1,klat
      do i=1,klon
C
C*    Special solutions of nonhomogeneous equation *wnon*
C      and homogenous equation  *whom*
C
       do k=0,klev
        wr1(k)=0.
        do l=1,klev
         wr1(k) = wr1(k)+EMI(k,l)*skf(i,j,l)
        enddo
       enddo
       do k=0,klev
        wnons=0.
        whoms=0.
        do l=0,klev
         wnons = wnons+EMO(k,l)*wr1(l)/(vlam(l) - xylamb(i,j))
         whoms = whoms+EMO(k,l)*EMI(l,klev)/(vlam(l) - xylamb(i,j))
        enddo
        wnon(k)=wnons
        whom(k)=whoms
       enddo
C
C*  Weight whoms of homogeneous solution
C
       wnons=0.0
       whoms=0.0
       do k=1,klev
          wnons = wnons+wnon(k)*ds1(k)
          whoms = whoms+whom(k)*ds1(k)
       enddo
       whoms=(-wnons+b(i,j))/whoms
C
C*  Final mass-balanced solution
C
       do k=0,klev
          skf(i,j,k)=wnon(k)+whoms*whom(k)
       enddo
c
      enddo
      enddo
c
      do k=0,klev
       skf(1,1,k) = 0.0
      enddo
C
      RETURN
      END


C *tqli* from: NUM RECIPES in Fortran + pythac
C
      SUBROUTINE tqli(d,e,n,np,z)
c
      implicit real*8(a-h,o-z)
C      INTEGER n,np
      dimension d(np),e(np),z(np,np)
CU    USES pythag
C      INTEGER i,iter,k,l,m
C      REAL b,c,dd,f,g,p,r,s,pythag
      do 11 i=2,n
        e(i-1)=e(i)
11    continue
      e(n)=0.
      do 15 l=1,n
        iter=0
1       do 12 m=l,n-1
          dd=abs(d(m))+abs(d(m+1))
          if (abs(e(m))+dd.eq.dd) goto 2
12      continue
        m=n
2       if(m.ne.l)then
          if(iter.eq.30) write(*,*)'too many iterations in tqli'
          iter=iter+1
          g=(d(l+1)-d(l))/(2.*e(l))
          r=pythag(g,1.d0)
          g=d(m)-d(l)+e(l)/(g+sign(r,g))
          s=1.
          c=1.
          p=0.
          do 14 i=m-1,l,-1
            f=s*e(i)
            b=c*e(i)
            r=pythag(f,g)
            e(i+1)=r
            if(r.eq.0.)then
              d(i+1)=d(i+1)-p
              e(m)=0.
              goto 1
            endif
            s=f/r
            c=g/r
            g=d(i+1)-p
            r=(d(i)-g)*s+2.*c*b
            p=s*r
            d(i+1)=g+p
            g=c*r-b
C     Omit lines from here ...
            do 13 k=1,n
              f=z(k,i+1)
              z(k,i+1)=s*z(k,i)+c*f
              z(k,i)=c*z(k,i)-s*f
13          continue
C     ... to here when finding only eigenvalues.
14        continue
          d(l)=d(l)-p
          e(l)=g
          e(m)=0.
          goto 1
        endif
15    continue
      return
      END
C
C  function pythag
C
      FUNCTION pythag(a,b)
      implicit real*8(a-h,o-z)
C      REAL a,b,pythag
C      REAL absa,absb
      absa=abs(a)
      absb=abs(b)
      if(absa.gt.absb)then
        pythag=absa*sqrt(1.+(absb/absa)**2)
      else
        if(absb.eq.0.)then
          pythag=0.
        else
          pythag=absb*sqrt(1.+(absa/absb)**2)
        endif
      endif
      return
      END

      
