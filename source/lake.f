!-----ONE-DIMENSIONAL LAKE MODEL Lake1.0beta
!-----Moscow State University, Scientific Research Computing Center
!-----Moscow, 119992, GSP-2, Vorobievy Gory, SRCC MSU	
!-----Victor M. Stepanenko, stepanen@srcc.msu.ru
!-----Release: May 2007

      SUBROUTINE LAKE(tempair1,humair1,pressure1,uwind1,vwind1,
     &longwave1,sw,prec,zref1,tsw,hw1,xlew1,cdmw1,trib_inflow1,
     &dtl,ix,iy,nx,ny)
     
!     ENTRY SUBROUTINE OF THE MODEL 

!     In atmospheric model must be called once each timestep,
!     or once each N timesteps, where N = dt_lake/dt_atmos 
!     (dt_lake - timestep in the lake model, dt_atmos - timestep the atmospheric model) 

!     INPUT VARIABLES:
!     tempair1      --- air temperature, K;
!     humair1       --- specific humidity of air, kg/kg;
!     pressure1     --- air pressure, Pa;
!     uwind1        --- x-component of wind, m/s;
!     vwind1        --- y-component of wind, m/s;
!     longwave1     --- longwave downward radiation, W/m**2;
!     sw            --- net solar radiation, W/m**2;
!     prec          --- precipitation intensity, m/s;
!     xref1         --- the height (m) of level in atmosphere, where temperature, 
!                       humidity and wind are measured (calculated by atmospheric model);
!     trib_inflow1  --- tributary inflow rate, m/s;
!     dtl           --- timestep (for stability must be not larger than 10-15 min), s;
!     ix            --- current number of grid point in X-direction;
!     iy            --- current number of grid point in Y-direction;
!     nx            --- total number of grid points in X-direction of atmospheric model;
!     ny            --- total number of grid points in Y-direction of atmospheric model;

!     OUTPUT VARIABLES:   
!     tsw           --- surface temperature of a lake, K;
!     hw1           --- sensible heat flux from lake, W/m**2; 
!     xlew1         --- latent heat flux from lake, W/m**2; 
!     cdmw          --- exchange coefficient, m/s ( == wind_speed*exchange_nondimensional_coeficient)

      USE phys_parameters
	USE phys_constants 
	USE numeric_params
	USE phys_constants2
	use driving_params
	use atmos
	use arrays
	implicit none
	REAL*8 tempair1,humair1,pressure1,uwind1,vwind1,longwave1,
     & sw,prec,zref1,tsw,hw1,xlew1,cdmw1,dtl,trib_inflow1
	integer(4) nx,ny,ix,iy
	logical firstcall
	data firstcall /.true./

	SAVE

	tempair = tempair1-273.2
	humair = humair1
	pressure = pressure1
	trib_inflow = trib_inflow1
	uwind = uwind1
	if (uwind==0) uwind=0.1
	vwind = vwind1
      if (vwind==0) vwind=0.1
	wind = dsqrt(uwind**2+vwind**2)
	longwave = longwave1
	zref = zref1

      if (firstcall) then
	 print*, 'Lake: first call'
	 allocate (init(nx,ny))
	 init=0
	 call DEFINE_phys_parameters
	 call DEFINE_phys_constants
	 call DEFINE_phys_constants2
	 call DEFINE_driving_params
	 call allocarrays
	endif
			      
	CALL MAIN(sw,prec,tsw,dtl,ix,iy,nx,ny,firstcall)
      hw1=hw;xlew1=xlew;cdmw1=cdmw

	END SUBROUTINE LAKE


      SUBROUTINE MAIN(sw,prec,tsw,dtl,ix,iy,nx,ny,firstcall)
      
!     MAIN SUBROUTINE OF THE MODEL 
!     Calculates all the parameters of state of a lake (temperature, currents,
!     eddy diffusivity coefficients, snow, ice, soil characteristics, etc.) 
!     at the next timestep (i.e. implements transition of variables 
!     from i-th time step to (i+1)-th )      
      
      use numeric_params    
	use driving_params
	use atmos
	use phys_constants2
	use arrays
	implicit none


!------------------------ MAIN VARIABLES ------------------------
!  arrays:    Tw1 and Tw2 - temperature of water, C
!             Ti1 and Ti2 - temperature of ice, C
!             T - temperature of snow, C
!  functions of time:h1 and h2 - thickness of water, m
!             l1 and l2 - thickness of ice, m
!             hs1 - thickness of snow, m                    
!             flag_ice_surf - shows if ice is melting on the top surface, n/d        
!  constants: cw - specific heat of water, J/(kg*K) 
!             ci - specific heat of ice, J/(kg*K)
!             row - density of water, kg/m**3
!             roi - density of ice, kg/m**3
!             lamw - thermal conductivity of water, J*sec/(m**3*K)
!             lami - thermal conductivity of ice, J*sec/(m**3*K)
!             L - specific heat of melting, J/kg
!             Lwv - specific heat of evaporation, J/kg
!             Liv - specific heat of sublimation, J/kg    
!             ddz - size of grid element (space), m
!             dt - size of grid element (time), sec
!             kstratwater - coefficient for linear initial conditions in water
!             kstratice - coefficient for linear initial conditions in ice     	     
!  boundary conditions:
!             eFlux - total heat flux on the upper surface,
!      		  shortwave(1-A)-T**4+longwave-H-LE, J/(m**2*sec)
!             Elatent - latent heat flux, J/(m**2*sec) 
!             Erad - shortwave radiation, penetrated below a surface, J/(m**2*sec)
!             tempair - air temperature (2 m), C
!             precip - precipitation, m/sec
!            (M) number of layers in water and snow 	
	            
	real(8), allocatable::u_2d(:,:,:),v_2d(:,:,:),Tsoil1_2d(:,:,:),
     & wi1_2d(:,:,:),wl1_2d(:,:,:),Tw1_2d(:,:,:),Ti1_2d(:,:,:),
     & Tis1_2d(:,:,:),dz_2d(:,:,:),T_2d(:,:,:),wl_2d(:,:,:),
     & dens_2d(:,:,:),E_2d(:,:,:),eps_2d(:,:,:),dhwfsoil_2d(:,:),
     & Sal1_2d(:,:,:)
	real(8), allocatable::l1_2d(:,:),h1_2d(:,:),ls1_2d(:,:),
     & hs1_2d(:,:),time_2d(:,:) 
      real(kind(0.d0)), allocatable:: dep_2d(:,:),cdm2(:,:)
	integer(4), allocatable::fl_sn_2d(:,:),fl_sn_init_2d(:,:),
     & itop_2d(:,:)
      real(8) l2,h2,h1,dhwhigh,dhwlow,dhwp,dhwe,
     & dhihigh,dhilow,dhip,
     & totalevap,dhis,totalprecip,totalmelt,
     & dhwfsoil,totalpen,totalwat,xx,time_fut,dhwf,dhif,
     & ls1, ls2, ec, addt, snowmass, snowmass_init,
     & month,totalsurfrad,totallongwave,dhi0
     & totalerad,totalhflux,totalpen1,Tsnow,
     & Tf_old1,Tf_old2,Tf_old12,Tf_old22,
     & ma,mi,pi,totevap_real,toteeff_real,tothflux_real,
     & toterad_real,t1,cv,ts,cloud,wind10,day,
     & hflux_obs, Elatent_obs, transf,
     & wl_max,unfrwat,wi_max,dhwls,
     & totalerad,dep_av,totalsw,nyear,nmonth,nday,nhour	 
	real(8), dimension(1:350)::a,b,c,d 
	real(8), dimension(1:50):: z_sol,T_sol,Tsum_sol 

	real(8) cs(ms),lams(ms) 
	real(8) dz
      real(8) ST,PGR,TGROLD,QGROLD,RADIAT,PRECIP,WSOLD,SNOLD,
     & ZS,thsoil,whsoil,bOLD,RF1,RF2,SNMELT,HSNOW,
     & HS,ES,TGRNEW,QGRNEW,WSNEW,SNNEW,RUNOFF,
     & ElatOld,HFold,PRSold,shortwave,extinct
      real(8) AL,DLT,DVT,ALLL,DL,
     & ALV,DV,Z,T,WL,WV,WI,dens
      real(8) dt,l1,hs1,ddz
      real(8) q(110),alp(10),psi(10) 
	real(8) roughness,emissivity,albedo,aM,bM,relhums
      real(8), dimension(1:350):: Temp
     	real(8) sw,prec,tsw,dtl
	real(8) x1,x2,x3,x4,x5
	
	real(8), external:: meltpnt

      integer(4) i,j,p,r,y,ki,i_ML,
     & year,flag_snow,itop,mon, ac,bc,cc,dc,
     & spin_up,mix,reflect,nmonth_old, nhour_old,iter,mntsw,
     & mntsi,mnts,xday,xhour,i1,mnKT,ix,iy,nx,ny,
     & n_1cm,n_5cm,layer_case
      
	character infile,outfile,dir_in*4,fext2*40
	character*20 timedata
	character*40 fnmap
	     
      logical flux_calc,indstab,ind_bound,ts_iter_more,
     & soil_ts_ext_iter,firstcall,uniform_depth

      common /watericesnowarr/ lams,q
	common /SOILDAT/ dz(ms),itop
      common /bL/ ST,PGR,TGROLD,QGROLD,RADIAT,WSOLD,SNOLD,
     & ZS,thsoil,whsoil,bOLD,RF1,RF2,SNMELT,HSNOW,
     & HS,ES,TGRNEW,QGRNEW,WSNEW,SNNEW,RUNOFF,
     & ElatOld,HFold,PRSold,extinct(ms)
      common /SOILSOL/ AL(ML),DLT(ML),DVT(ML),ALLL(ML),DL(ML),
     & ALV(ML),DV(ML),Z(ML),T(ML),WL(ML),WV(ML),WI(ML),
     & dens(ms)
	common /ymdh/ nyear,nmonth,nday,nhour
	common /ts_ext_iter/ Tf_old1,soil_ts_ext_iter 
	common /layers/ h1,l1,hs1,ls1
	common /atmos2/ shortwave, precip
	common /surface/ roughness,emissivity,albedo,aM,bM,relhums
	
	data uniform_depth /.true./
	 
      SAVE
  
!BODY OF PROGRAM

	rosoil = 1200.	
	extwat = 1. !0.7 !0.9 !0.524  
	extice = 2.
	pi = 4.*atan(1.)

	year  = 60*60*24*365
	month = 30*24*60*60
	day   = 24*60*60
	         
	fnmap = 'sam_map'		      
      ddz = 1./float(M)
	
	shortwave=sw
c converting of precipitation units: from kg/m**s/sec to m/sec
	precip=prec !*row

	Erad = shortwave*(1-albedo)
	dt=dtl
	
	if (firstcall) then	
	 call pbldat
	 call comsoilforlake      
	 allocate (l1_2d(nx,ny))
	 allocate (h1_2d(nx,ny))
	 allocate (ls1_2d(nx,ny))
       allocate (dep_2d(0:nx-1,0:ny-1))
       allocate (cdm2(0:nx,0:ny))
	 allocate (u_2d(M+1,nx,ny))
	 allocate (v_2d(M+1,nx,ny))
	 allocate (E_2d(M+1,nx,ny))
	 allocate (eps_2d(M+1,nx,ny))
	 allocate (Tsoil1_2d(ns+1,nx,ny))
	 allocate (wi1_2d(ns+1,nx,ny))
	 allocate (wl1_2d(ns+1,nx,ny))
	 allocate (Tw1_2d(M+1,nx,ny))
	 allocate (Sal1_2d(M+1,nx,ny)) 
	 allocate (Ti1_2d(M+1,nx,ny))
       allocate (Tis1_2d(M+1,nx,ny))
	 allocate (fl_sn_2d(nx,ny))
	 allocate (fl_sn_init_2d(nx,ny))
	 allocate (itop_2d(nx,ny))
	 allocate (hs1_2d(nx,ny))
	 allocate (dz_2d(ms,nx,ny))
	 allocate (T_2d(ml,nx,ny))
	 allocate (wl_2d(ml,nx,ny))
	 allocate (dens_2d(ms,nx,ny))
	 allocate (time_2d(nx,ny))
	 allocate (dhwfsoil_2d(nx,ny))
	 totalsw = 0.
	 if (runmode == 2.and.uniform_depth == .false.) then
	  open (999, file=path(1:len_trim(path))//fext2(fnmap,'dep'),
     &     	     status='old')
	  call readgrd2 (999,dep_2d,0,nx-1,0,ny-1)
	  close (999)
	  dep_av=0.;x1=0.
	  do i=1,nx-1
	   do j=1,ny-1
          if (dep_2d(i,j)>0.) then
		   dep_av=dep_av+dep_2d(i,j)
	     x1=x1+1.
          endif
	   enddo
        enddo
	  dep_av=dep_av/x1
	 endif
       firstcall=.false.
      endif
      
	if (init(ix,iy)==0) then
!     Specification of initial profiles	
       h1 = h10
	 l1 = l10
	 hs1 = hs10
	 ls1 = ls10
       do i=1,M+1
        u1(i)=0.1*(1-float(i)/float(M+1))
        v1(i)=0.1*(1-float(i)/float(M+1))
       enddo 
       do i=1,M 
	  E1(i) = 1.d-5 + float(i)/float(M)*(1.d-8  -  1.d-5)
	  eps1(i)= 1.d-9 + float(i)/float(M)*(1.d-12  -  1.d-9)
	 enddo
	 do i=1,ns
	  if (i<0.7*ns) then 
         Tsoil1(i) = Tb0 + float(i-1)/(0.7*float(ns))*(Tbb0-Tb0)
        else
	   Tsoil1(i) = Tbb0
	  endif 
        if (Tsoil1(i)>0) then
	   wi1(i) = 0.
	   wl1(i) = wl_max(por(i),rosoil(i),0.d0)-0.01
	  else
	   wl1(i) = unfrwat(Tsoil1(i),i)
         wi1(i) = wi_max(por(i),rosoil(i))-wl1(i)-0.01
	  endif
       enddo
		 
	 i_ML = int4(M*h_ML0/h1)
	 do i=1, max(i_ML,1)
	  Tw1(i) = Ts0
	  Sal1(i) = Sals0
       enddo
	 do i=max(i_ML+1,2),M+1
	  Tw1(i) = Ts0 + (Tb0-Ts0)*float(i-i_ML)/float(M+1-i_ML)
	  Sal1(i) = Sals0 + (Salb0-Sals0)*float(i-i_ML)/float(M+1-i_ML)
	 enddo 
	 Ti1 = 0.
	 Tis1 = 0 

	 if (l1/=0) then
	  do i=1,M+1
         Ti1(i)=-5.+5.*float(i-1)/float(M)
	  enddo
	 endif
	  		
	 if (hs1==0.) then	 
	  flag_snow = 0
	  flag_snow_init = 1
       else
	  flag_snow = 1
	  flag_snow_init = 0
        hs1 = dmax1(2.d-2, hs1)
	  hs1 = int4(hs1/0.01)*0.01
	  n_5cm = int4(hs1/0.05)
	  n_1cm = int4((hs1-n_5cm*0.05)/0.01)
	  if (n_1cm==0) then
	   n_1cm=1; hs1=hs1+0.01
        endif
	  itop = ms - (n_1cm+n_5cm)
 	  do i=itop,itop+n_1cm-1
         dz(i) = 0.01
         T(i) = -5.
	   wl(i) = 0
	   dens(i) = 150.
        enddo
	  do i=itop+n_1cm,ms-1
	   dz(i) = 0.05
         T(i) = -5.
	   wl(i) = 0
	   dens(i) = 150.
        enddo
       endif
  
        Tsum_sol = 0
	  mnts = 0 
        cdmw2 = 1.d-3 

	  totalevap = 0 
	  totalmelt = 0 
	  totalprecip = 0. 
	  totalwat = 0
	  totalpen = 0
	  totalpen1 = 0
	  time = 0
	  dhwfsoil = 0.

	  par = 1
        mon = 1  
	  iter = 0 
	  nmonth_old = 1
	  nhour_old = 21
	  xday=1
        xhour=0
       if (runmode==2.and.uniform_depth==.false.) then
	  if (ix<=nx-1.and.iy<=ny-1) then
	   h1 = dep_2d(ix,iy)
        else
	   h1 = dep_av
	  endif  
	  if (h1<=0) then
	   print*, 'Negative or zero initial depth at ix=',ix,
     &   'iy=',iy,'h1=',h1,':terminating program'
         STOP 
        endif
	 endif
	else
!ALIGNING VALUES FROM PREVIOUS TIME STEP (time) IN CURRENT POINT (ix,iy)
       l1 = l1_2d(ix,iy) 
       h1 = h1_2d(ix,iy)
       ls1 = ls1_2d(ix,iy) 

	 do i=1,M+1
	  u1(i)=u_2d(i,ix,iy)
	  v1(i)=v_2d(i,ix,iy)
	 enddo

	 do i=1,M+1
        E1(i)=E_2d(i,ix,iy)
	  eps1(i)=eps_2d(i,ix,iy)
	 enddo
	  	            
       do i=1,ns 
        Tsoil1(i)=Tsoil1_2d(i,ix,iy)
        wi1(i)=wi1_2d(i,ix,iy)
	  wl1(i)=wl1_2d(i,ix,iy)
       enddo
		 	  
	 do i=1,M+1 
	  Tw1(i) = Tw1_2d(i,ix,iy)
	  Sal1(i) = Sal1_2d(i,ix,iy)
        Ti1(i) = Ti1_2d(i,ix,iy)
	  Tis1(i) = Tis1_2d(i,ix,iy)
       enddo 
      			 
	 flag_snow = fl_sn_2d(ix,iy)
	 flag_snow_init = fl_sn_init_2d(ix,iy)

       itop = itop_2d(ix,iy)
	 hs1 = hs1_2d(ix,iy)

       do i=max(1,itop),ms-1
        dz(i) = dz_2d(i,ix,iy)
        T(i) = T_2d(i,ix,iy)
	  wl(i) = wl_2d(i,ix,iy)
	  dens(i) = dens_2d(i,ix,iy)
       enddo
  
        Tsum_sol = 0
	  mnts = 0 
        cdmw2 = cdm2(ix,iy)

	  totalevap = 0 
	  totalmelt = 0 
	  totalprecip = 0. 
	  totalwat = 0
	  totalpen = 0
	  totalpen1 = 0
	  
        time = time_2d(ix,iy)
	  dhwfsoil = dhwfsoil_2d(ix,iy)
 	endif
     			                
!THE MAIN CYCLE		
         
      time = time + dt
      nstep = nstep + 1
      
      layer_case = 1
      if (h1 > 0  .and. l1 > 0) layer_case = 2 
      if (h1 ==0  .and. l1 > 0) layer_case = 3
      if (h1 ==0  .and. l1 ==0) layer_case = 4
      
      if (layer_case == 1.and.
     & Tw1(1)<Meltpnt(Sal1(1))-0.1.and.h1-0.01*roi/row>0.01) 
     & layer_case = 2
     
      if (layer_case == 3.and.
     & Ti1(M+1)>Meltpnt(0.d0)+0.1.and.l1-0.01*row/roi>0.01) 
     & layer_case = 2
      
	if (layer_case == 3.and.
     & Ti1(1)>Meltpnt(0.d0)+0.1.and.l1-0.01*row/roi>0.01.
     & and.flag_snow==0) then
	 h1 = 0.01
	 ls1 = l1 - 0.01*row/roi
	 Tis1 = Ti1
	 l2 = 0
	 Ti2 = 0
	 if (ls1<0.009999) print*, 'Thin layer! - ls'
	 layer_case = 1 
      endif
  
      if1: IF (layer_case == 1) THEN 
      
!---------------------------- CASE 1: LIQUID WATER LAYER AT THE TOP ("Summer")-----------------------

      if (Tw1(M+1)<Meltpnt(Sal1(M+1))-0.1.and.ls1 == 0.and.
     & h1-0.01*roi/row>0.01) then
	 ls1 = 0.01
	 h1 = h1 - 0.01*roi/row
	endif	
      
	do i=2,M
	 if (Tw1(i)<Meltpnt(Sal1(i))) Tw1(i)=Meltpnt(Sal1(i))
	enddo
	
      call KT_eq(dt)
	do i=1,M+1
	 lamw(i)=lamw0+KT(i) !+KC(i) !+KLengT(i)
	enddo
	lamsal = lamw/(cw*row) * alsal
      
	do i=1,M
       SR(i)=Erad*(1-sabs)*dexp(-extwat*(i-0.5)*ddz*h1)
	enddo
	do i=2,M
	 RS(i)=-(SR(i)-SR(i-1))/(ddz*h1)
      enddo
	if (ls1==0) then
	 ice=0;snow=0;water=1;deepice=0
       dhwp = precip*dt
	 dhwe = - Elatent/(Lwv*row)*dt
	 dhw = dhwp + dhwe + dhwfsoil
	 dhw0 = dhwp + dhwe
	else
	 ice=0;snow=0;water=1;deepice=1
	 dhwls = -lamw(M)*dt/(ddz*h1*row*Lwi)*(Tw1(M+1) - Tw1(M))+
     & lami*dt/(ddz*ls1*row*Lwi)*(Tis1(2) - Tis1(1))+
     & (1-albedoofice)*Erad*(1-sabs)*dexp(-extwat*h1)/(row*Lwi)*dt 
       dls = -dhwls*row/roi
	 dls0 = dls
	 dhwp = precip*dt
	 dhwe = - Elatent/(Lwv*row)*dt
	 dhw = dhwp + dhwe + dhwfsoil + dhwls
	 dhw0 = dhwp + dhwe
	endif
	
	call S_diff(dt)
	 
	call soilforlake(dt, h1, l1, ls1, hs1,dhwfsoil,
     & a,b,c,d,Temp,M+1)
						
	h2 = h1 + dhw
	ls2 = ls1 + dls
	
      ENDIF if1               

      if2: IF (layer_case==2) THEN
      
!---------------------------CASE 2: WATER AND ICE ("Winter")------------------------------

       if (l1 == 0.) then
       l1 = 0.01 
	 h1 = h1-0.01*roi/row
	 if (tempair<0.) then
	  do i=1, M+1
	   Ti1(i) = tempair - float(i-1)/float(M)*tempair
	  enddo 
       else 
	  Ti1=0.
       endif
	end if 
	
	do i=2,M
	 if (Tw1(i)<Meltpnt(Sal1(i))) Tw1(i)=Meltpnt(Sal1(i))+1.E-5
	enddo
		  
	if (h1 == 0.) then
	 h1 = 0.01
	 l1 = l1-0.01*row/roi
	end if

	!if (h1<0.009999.or.l1<0.009999) print*, 'Thin layer!-l or h'
		  
	if (Tw1(M+1)<Meltpnt(Sal1(M+1))-0.1.and.ls1 == 0.and.
     &	h1-0.01*roi/row>0.01) then
	 ls1 = 0.01
	 h1 = h1 - 0.01*roi/row
	endif

      if (Turbpar/=1) then
    	call KT_eq(dt)
    	 do i=1,M+1
	  lamw(i)=lamw0+KT(i) !+KC(i)
	 enddo	
      else
	 lamw=lamw0*3.
      endif
      lamsal = lamw/(cw*row) * alsal
   
!CASE 2.1: WATER, ICE AND SNOW 

	if (flag_snow == 1) then
	 SR=0;RS=0
	 dhwlow = lamw(1)*dt/(ddz*h1*row*Lwi)*(Tw1(2) - Tw1(1))-lami*dt/
     & (ddz*l1*row*Lwi)*(Ti1(M+1) - Ti1(M))
       dhilow = -dhwlow*row/roi
	 if (Ti1(1)>Meltpnt(0.d0)+0.1) then
	  dhwhigh = (Ti1(1)-Meltpnt(0.d0))*ci*roi*l1*ddz/2/(Lwi*row)
	  dhihigh = -dhwhigh*row/roi
       else
	  dhwhigh = 0
	  dhihigh = 0
       endif
	 if (ls1/=0) then
	  ice=1;snow=1;water=1;deepice=1
	  dhwls=-lamw(M)*dt/(ddz*h1*row*Lwi)*(Tw1(M+1) - Tw1(M))+lami*
     &  dt/(ddz*ls1*row*Lwi)*(Tis1(2) - Tis1(1))
        dls = -dhwls*row/roi
	  dhw = dhwhigh + dhwlow + snmelt*dt + dhwfsoil + dhwls
	  dhw0 = dhwlow+dhwhigh+snmelt*dt
	  dls0 = dls
	 else
	  ice=1;snow=1;water=1;deepice=0
        dls = 0
	  dhwls = 0
	  dhw = dhwhigh + dhwlow + snmelt*dt + dhwfsoil
	  dhw0 = dhwlow+dhwhigh+snmelt*dt
       endif
	 dhi = dhilow + dhihigh 
	 dhi0 = dhihigh
	 
	 call S_diff(dt)
	 
	 call snowtemp(snowmass,snowmass_init,a,b,c,d,Temp,dt)

       h2 = h1 + dhw
	 l2 = l1 + dhi 
	 ls2 = ls1 + dls
	 
      else
	     
!CASE 2.2: WATER AND ICE WITHOUT SNOW

       do i=1,M
        SR(i)=Erad*dexp(-extice*l1)*dexp(-extwat*(i-0.5)*ddz*h1)
	 enddo
	 do i=2,M
	  RS(i)=-(SR(i)-SR(i-1))/(ddz*h1)
       enddo
	 if (Ti1(1)>Meltpnt(0.d0)+1.E-5) then
        dhwhigh = (Ti1(1)-Meltpnt(0.d0))*ci*roi*ddz*l1/2/(Lwi*row)
!	  Ti1(1) = -1.E-5
!	  do i=2,M
!	   if (Ti1(i)>1.E-5) then
!	    dhwhigh = dhwhigh+Ti1(i)*ci*roi*l1*ddz/(Lwi*row)
!          Ti1(i)=-1.E-5
!         endif
!	  enddo
	  dhis = 0 
	  dhihigh = -dhwhigh*row/roi
       endif	
       if (precip>0.and. tempair<0.and.flag_snow /= 1.and.l1>0.05) then
	  flag_snow = 1
	  hs1 = 0.02
	  itop = ms-2
	  dhip = 0
       else
	  dhip = precip*row/roi*dt
	 end if      
	 dhwlow = lamw(1)*dt/(ddz*h1*row*Lwi)*(Tw1(2) - Tw1(1))-lami*dt/
     & (ddz*l1*row*Lwi)*(Ti1(M+1) - Ti1(M))
       dhilow = -dhwlow*row/roi
       if (Ti1(1)<=Meltpnt(0.d0)+0.1) then
	  dhis = - Elatent/(roi*Liv)*dt
	  dhwhigh = 0
	  dhihigh = 0
       endif
	 dhi = dhihigh + dhilow + dhis + dhip
	 dhi0 = dhis + dhihigh + dhip
	 if (ls1/=0) then
	  ice=1;snow=0;water=1;deepice=1 
	  dhwls =-lamw(M)*dt/(ddz*h1*row*Lwi)*(Tw1(M+1) - Tw1(M))+lami*
     &  dt/(ddz*ls1*row*Lwi)*(Tis1(2) - Tis1(1))
        dls = -dhwls*row/roi
	  dls0 = dls
	  dhw = dhwhigh + dhwfsoil + dhwls + dhwlow
	  dhw0 = dhwhigh+dhwlow
	 else
	  ice=1;snow=0;water=1;deepice=0
        dhw = dhwhigh + dhwfsoil + dhwlow
	  dhw0 = dhwhigh+dhwlow
	  dls = 0.
	 endif
	
	 call S_diff(dt)
	
	 call soilforlake(dt, h1, l1, ls1, hs1,
     &  dhwfsoil,a,b,c,d,Temp,M+1)
    	
       h2 = h1 + dhw
	 l2 = l1 + dhi
	 ls2 = ls1 + dls
	 
       endif
       
       ENDIF if2

!CASE 3: ICE WITHOUT WATER

      if3: IF (layer_case==3) THEN
          
!CASE 3.1: ICE WITH SNOW
      
	if (flag_snow == 1) then
	 if (Ti1(1)>Meltpnt(0.d0)+0.1) then
	  dhwhigh = (Ti1(1)-Meltpnt(0.d0))*ci*roi*l1*ddz/2/(Lwi*row)
	  dhif = -dhwhigh*row/roi
       else
	  dhwhigh = 0
	  dhif = 0
       endif
	 dhw = snmelt*dt + dhwhigh
	 dhi = dhif
	 dhi0 = dhi
	 ice=1;snow=1;water=0;deepice=0
	 call snowtemp(snowmass,snowmass_init,a,b,c,d,Temp,dt)
       l2 = l1 + dhi 
	 h2 = h1 + dhw
      
      else
!CASE 3.2: ICE WITHOUT SNOW !
       if (precip>0.and. tempair<0) then
	  flag_snow = 1
	  hs1 = 0.02
        dhip = 0
       else
	  dhip = precip*row/roi*dt
	 end if
	 if (Ti1(1)>Meltpnt(0.d0)+0.1) then
	  dhw = (Ti1(1)-Meltpnt(0.d0))*ci*roi*l1*ddz/2/(Lwi*row)
	  dhihigh = - dhw*row/roi
	  dhis = 0
       else
        dhis = - Elatent/(roi*Liv)*dt
	  dhw = 0
	  dhihigh = 0
       endif
	 dhi = dhis + dhihigh + dhip
	 dhi0 = dhi
	 ice=1;snow=0;water=0;deepice=0 
	 call soilforlake(dt, h1, l1, ls1, hs1,
     & dhwfsoil,a,b,c,d,Temp,M+1)
	 			        	 		   	 				
	 l2 = l1 + dhi
	 h2 = h1 + dhw
	 
      endif
      
      ENDIF if3 

! CASE 4:SNOW ANS SOIL WITHOUT ICE AND WATER
			  
 5    if4: IF (layer_case==4) THEN
       if (flag_snow==1) then
    	  ice=0;snow=1;water=0;deepice=0
	  call snowtemp(snowmass,snowmass_init,a,b,c,d,Temp,dt)
       else
	  ice=0;snow=0;water=0;deepice=0
	  call soilforlake(dt, h1, l1, ls1, hs1,
     &   dhwfsoil,a,b,c,d,Temp,M+1)
	 endif 
	ENDIF if4

!     TRIBUTARY INFLOW to a lake
      h2 = h2 + trib_inflow/year*dt 
	      
	if (h2 < 0.005.and.ls2/=0.and.l2==0) then
       ls2 = ls2 + h2*row/roi
	 if (ls2 < 0.005) then 
	  ls2 = 0; Tis2 = 0
	 endif
	 h2 = 0.
	 Tw2 = 0.
	 Sal2 = 0.
      endif

	if ((h2<0.005.and.h2>0).and.l2 /= 0) then
	 l2 = l2 + h2*row/roi
	 if (l2 < 0.005) then 
	  l2 = 0; Ti2 = 0
	 endif
	 h2 = 0.
	 Tw2 = 0.
	 Sal2 = 0. 
	end if
			  
	if (h2 < 0.005.and.l2/=0.and.ls2/=0) then
	 l2 = l2 + h2*row/roi + ls2 
	 addt=0
	 do i=1,M+1
	  if ((i==1).or.(i==M+1)) then
         addt = addt + Tis2(i)*ls2/l2*ddz/2+Ti2(i)*(l1/l2-1)*ddz/2
	  else
	   addt = addt + Tis2(i)*ls2/l2*ddz+Ti2(i)*(l1/l2-1)*ddz
	  endif  
       enddo
	 Ti2 = Ti2 + addt 
	 ls2 = 0
	 Tis2 = 0.
	 h2 = 0
	 Tw2 = 0.
	 Sal2 = 0.
      endif
			 
      if (l2 < 0.005.and.h2 /= 0) then
	 h2 = h2 + l2*roi/row
	 if (h2 < 0.005) then 
	  h2 = 0; Tw2 = 0; Sal2 = 0.
	 endif
	 l2 = 0.
	 Ti2 = 0.
	end if

	if (h2 < 0.005.and. l2 == 0) then
	 h2 = 0.
	 Tw2 = 0.
	 Sal2 = 0. 
	end if
	
	if (l2 < 0.005.and. h2 == 0) then
	 l2 = 0.
	 Ti2 = 0.
	end if
			 
	if ((hs1 < 0.005.and.hs1>0).or.(hs1>0.and.l2==0.and.h2/=0)) then
	 h2 = h2 + (totalprecips-totalevaps-totalmelts)
	 if (h2 < 0.005) then
	  l2 = l2 + row*h2/roi
	  h2 = 0; Tw2 = 0; Sal2 = 0.
	  if (l2 < 0.005) then 
	   l2 = 0; Ti2 = 0
	  endif
	 endif
	 hs1 = 0.
	 flag_snow = 0. 
	 flag_snow_init = 1
      end if
		 			         
	if (ls2 < 0.005.and.h2/=0) then
	 h2 = h2 + ls2*roi/row
	 if (h2 < 0.005) then 
	  h2 = 0; Tw2 = 0; Sal2 = 0.
	 endif
	 ls2 = 0
	 Tis2 = 0.
      endif 

	if (ls2>0.and.h2 == 0) then
	 l2 = ls2
	 Ti2 = Tis2
	 if (l2 < 0.005) then 
	  l2 = 0; Ti2 = 0
	 endif
	 ls2 = 0
	 Tis2 = 0.
      endif
		      
	h1 = h2
	l1 = l2
	ls1 = ls2
	Tw1 = Tw2
	Ti1 = Ti2
	Tis1 = Tis2
	Sal1 = Sal2
   	
	! CALCULATION OF SUMMARY FLUXES !			
      totalpen = totalpen+dhwfsoil 
	totalpen1 = totalpen1+dhwfsoil
	totalevap = totalevap+Elatent/(row*Lwv)*dt
	totalprecip = totalprecip+precip*dt
!	totalsurfrad = totalsurfrad + surfrad*dt
!	totallongwave = totallongwave + longwave*dt
	totalhflux = totalhflux + hflux*dt
	totalerad = totalerad + erad*dt
	
      
	ki=22
	do i=2,22
	 if (z_sol(i)>l1+h1) then
	  ki=i-1
        exit
       endif
      enddo
      
	do i=2,ki
       if (l1>z_sol(i)) then
	  j=1
1122	  j=j+1
	  if ((j-1)*ddz*l1>z_sol(i)) then
	   T_sol(i)=Ti1(j-1)+(Ti1(j)-Ti1(j-1))/(ddz*l1)*
     &   (z_sol(i)-(j-2)*ddz*l1)
	  else
	   goto 1122
	  endif
       else
	  j=1
1123	  j=j+1 
	  if ((j-1)*ddz*h1+l1>z_sol(i)) then
	   T_sol(i)=Tw1(j-1)+(Tw1(j)-Tw1(j-1))/(ddz*h1)*
     &   (z_sol(i)-(j-2)*ddz*h1-l1)
	  else
	   goto 1123
	  endif
	 endif    
	enddo
	T_sol(1) = Tw1(1)
	if (l1/=0) T_sol(1) = Ti1(1)
	T_sol(ki+1) = Tw1(M+1)
	Tsum_sol = Tsum_sol + T_sol
	mnts = mnts + 1

      mnKT=mnKT+1
	
	if (l1==0) tsw=Tw1(1)+273.2
	if (l1/=0.and.flag_snow==0) tsw=Ti1(1)+273.2
	if (flag_snow==1) tsw=T(itop)+273.2
	
	if (init(ix,iy)==0) then
	 init(ix,iy)=1	
      endif
	 
	
! VALUES AT NEXT TIME STEP (time + dt) IN CURRENT POINT (ix,iy)
      
      l1_2d(ix,iy)  = l1 
	h1_2d(ix,iy)  = h1
      ls1_2d(ix,iy)  = ls1

	do i=1,M+1
	 u_2d(i,ix,iy)=u1(i)
	 v_2d(i,ix,iy)=v1(i)
	enddo

	do i=1,M+1
	 E_2d(i,ix,iy)=E1(i)
	 eps_2d(i,ix,iy)=eps1(i)
	enddo
	  	            
      do i=1,ns 
       Tsoil1_2d(i,ix,iy)=Tsoil1(i)
       wi1_2d(i,ix,iy)=wi1(i)
	 wl1_2d(i,ix,iy)=wl1(i)
      enddo
		 	  
	do i=1,M+1 
	 Tw1_2d(i,ix,iy) = Tw1(i)
	 Sal1_2d(i,ix,iy) = Sal1(i) 
       Ti1_2d(i,ix,iy) = Ti1(i)
	 Tis1_2d(i,ix,iy) = Tis1(i)
      enddo 
      			 
	fl_sn_2d(ix,iy) = flag_snow
	fl_sn_init_2d(ix,iy) = flag_snow_init

      itop_2d(ix,iy) = itop
	hs1_2d(ix,iy) = hs1

      do i=max(1,itop),ms-1
       dz_2d(i,ix,iy) = dz(i)
       T_2d(i,ix,iy) = T(i)
	 wl_2d(i,ix,iy) = wl(i)
	 dens_2d(i,ix,iy) = dens(i)
      enddo
      cdm2(ix,iy) = cdmw
      time_2d(ix,iy) = time
      dhwfsoil_2d(ix,iy) = dhwfsoil

!	if (ix==10.and.iy==10) print*, 'Lake',
!     & 'T1=',Tw1(1),'T2=',Tw1(2),'h1=',h1,'S=',shortwave,
!     & 'hh=',hflux,'LE=',Elatent,
!     & 'eflux=',eflux, 'Longwave=', longwave,
!     & 'Bal=',(shortwave*(1-albedoofwater)*(sabs)+longwave-
!     &  surfrad-hflux-Elatent-eflux)/(cw*row*ddz*h1/2)*dt
     
      if (runmode==1) then
       if (monthly_out ==1) call mon_out(dt) 
	 if (daily_out   ==1) call day_out(dt) 
	 if (hourly_out  ==1) call hour_out(dt)
	 if (everystep   ==1) call everystep_out
	 if (time_series ==1) call series_out(tsw)
	endif
      	      
	!FORMATS!
  7   format (f7.3, 4i5,37f7.3) 
  60  format (f5.2, 2e16.5, 2f8.3, e15.5, f11.5, f7.1)
  !61  format (i5, f6.2, 2e16.5, 2f8.3, e15.5, f11.5, f7.1)
  61  format (f6.2, f5.2, 2e16.5, 2f8.3, e15.5, f11.5, f7.1)
  62  format (f6.2, f5.2, 2f12.3, 2f8.3, f11.5, f7.1,3f12.4,2f7.2,3f7.1)
  80  format (f7.3, 10f9.2)
  90  format (f7.3, 12f9.2)   
  100 format (a5, 2a16, 2a8, a15, a11, a7) 
      !FINISHING PROGRAM!       
  
      END SUBROUTINE MAIN        

!-----------------------------------------------------------------------------------
!                  TEMPERATURE EQUATION SOLVER
!-----------------------------------------------------------------------------------

      SUBROUTINE T_SOLVER(dt)
      
!     T_SOLVER implements iterations to find water surface temperature 
!     at the next time step    

	use driving_params
	use arrays
      implicit none

	real(8) dt_scan,Tsurf1,Tsurf2,Tsurf3,Bal1,Bal2,Bal3,dt
	integer(4) iter,maxiter

	SAVE

      open (666,file=path(1:len_trim(path))//
     & 'results/debug/iter_T.dat')

      maxiter=10
      
!     SCANNING INTERVAL OF SURFACE TEMPERATURE [-90 C, ...]

      dt_scan = 5.
      Tsurf2=-90.
	Bal1=1.;Bal2=1.
	do while (Bal1*Bal2>=0)
	 Tsurf1=Tsurf2
	 call Bal_calc(Bal1,Tsurf1,dt)  
	 Tsurf2=Tsurf1+dt_scan
	 if (Tsurf2>50.) then
	  print*, 'Temperature limit in scan process is exceeded: STOP'
	  STOP
       endif
       call Bal_calc(Bal2,Tsurf2,dt)
	 continue
      enddo
      
	iter=0
	Bal3=1.
	
!     CHORDE METHOD TO FIND SURFACE TEMPERATURE

      c1:do while (dabs(Bal3)>0.1)
	 iter=iter+1
	 call Bal_calc(Bal1,Tsurf1,dt)
	 call Bal_calc(Bal2,Tsurf2,dt)    
	 Tsurf3=(Tsurf1*Bal2-Tsurf2*Bal1)/(Bal2-Bal1)
	 call Bal_calc(Bal3,Tsurf3,dt)
	 if (iter>maxiter) then
	  write (666,*) Bal3
	  exit c1
	 endif
	 if (Bal1*Bal3<0.) then
	  Tsurf2=Tsurf3
       elseif (Bal2*Bal3<0.) then
	  Tsurf1=Tsurf3
       endif
	enddo c1     

	call T_diff(0,Tsurf1,dt)
	
      END SUBROUTINE T_SOLVER
      

	SUBROUTINE BAL_CALC(Bal,Tsurf,dt)
	
!     BAL_CALC calculates the components of energy balance at the water surface	

	use atmos
	use driving_params
	use numeric_params
	use phys_constants2
	use arrays
	implicit none

	real(8) Bal,Tsurf,dt

      real(8), dimension(1:ms):: Tsn,cs,lams,q
	real(8) shortwave,precip 
	real(8) TET2,TET1,ro,esatsurf,humsurf,bx(7),bix(11),
     &c_u,c_t,ksw,hwave,SHF,LHF,dHdt,ddz
      real(8) h1,l1,hs1,ls1
	real(8) wind10,urel,vrel,u,v,wr
	real(8) Hm, LEm
	real(8) roughness,emissivity,albedo,aM,bM,relhums
	real(8) AL,DLT,DVT,ALLL,DL,
     & ALV,DV,Z,T,WL,WV,WI,dens,dz
	integer(4) itdrag,surftyp,itop

	common /layers/ h1,l1,hs1,ls1
	common /wind/ wind10,urel,vrel,u,v
	common /measflux/ Hm, LEm
	common /surface/ roughness,emissivity,albedo,aM,bM,relhums
	common /snow_char/ Tsn,cs
	common /atmos2/ shortwave, precip
	common /SOILSOL/ AL(ML),DLT(ML),DVT(ML),ALLL(ML),DL(ML),
     & ALV(ML),DV(ML),Z(ML),T(ML),WL(ML),WV(ML),WI(ML),dens(ms)
	common /SOILDAT/ dz(ms),itop
	common /watericesnowarr/ lams,q
	
	SAVE

      
      ddz=1./float(M)
      if (snow==1) then
	 surftyp=3
      elseif (ice==1) then
	 surftyp=2
      elseif (water==1) then
	 surftyp=4
	else
	 surftyp=1
      endif
	call SURF_CHAR(surftyp)
	call T_diff(1,Tsurf,dt)
	SELECT CASE (surftyp)   
	 case (1)
	  eflux = -(lamsoil(1)+lamsoil(2))/2.*(Tsoil2(2)-Tsoil2(1))/dzs(1)
	 case (2)
	  eflux = -lami*(Ti2(2)-Ti2(1))/(ddz*l1) 
	  !if (l1<0.1) eflux = -lami*(Ti2(2)-Ti2(1))/(ddz*0.1) 
	 case (3)
	  eflux = -(lams(itop)+lams(itop+1))/2.*(Tsn(itop+1)-Tsn(itop))
     &  /dz(itop)
       case (4)  
	  eflux = -lamw(1)*(Tw2(2)-Tw2(1))/(ddz*h1) 
	  !if (h1<0.1) eflux = -lamw(1)*(Tw2(2)-Tw2(1))/(ddz*0.1) 
      END SELECT
	
!    (urel,vrel) is wind, relative to lake currents

      if (relwind==2) then
	 urel=uwind-u
	 vrel=vwind-v
	elseif (relwind==1) then
       urel=uwind
	 vrel=vwind
      endif
	wr=dmax1(dsqrt(urel**2+vrel**2),0.1d0)
      wind10 = wr*log(10./0.01)/log(2./0.01) !0.01 m is water surface roughness
      TET2 = (tempair+273.15)*(100000./pressure)**0.286
	ro = pressure/(Rd*(tempair+273.15))  
      TET1 = (Tsurf+273.15)*(100000./pressure)**0.286  
	esatsurf = 610.7*10.**(aM*Tsurf/(bM+Tsurf))
	humsurf = 0.622/pressure*esatsurf*relhums
	if (PBLpar==1.or.PBLpar==11) then
	
c     Businger-Dayer, Beljaars parameterization for exchange coefficients

       itdrag=10
       bx(1) = wr !wind
	 bx(2) = TET2
	 bx(3) = TET1
	 bx(4) = humair
	 bx(5) = humsurf
	 bx(6) = zref
	 bx(7) = roughness 
       call dragvl (bx, bix, itdrag)
	 c_u = bix(8)
	 c_t = bix(9)
	 hflux = -cp*ro*c_u*c_t*wr*(TET2-TET1)
	 Elatent =  -Lwv*ro*c_u*c_t*wr*(humair-humsurf) 
	 tau = ro*c_u**2*wr**2
	endif

c     Parameterization of exchange coefficients acording to (Louis, 1979)

	if (PBLpar==2.or.PBLpar==21) 
     &call RichPBL(hflux,Elatent,tau,Tsurf,humsurf,roughness)
	
c     SHALLOW WATER EFFECT ON HEAT, MOISTURE AND MOMENTUM FLUXES (PANIN ET. AL., 2006)

      hwave = 0.07*wind10**2*(g*h1/dMAX1(wind10,1.d0)**2)**(3./5.)/g
	if (l1==0.and.h1/=0.and.(PBLpar==11.or.PBLpar==21)) then
	 ksw=2. 
	 hflux=hflux+hflux*ksw*hwave/h1
	 Elatent=Elatent+Elatent*ksw*hwave/h1
	 tau=tau+tau*ksw*hwave/h1
      endif
	cdmw=tau/(ro*wr);xlew=Elatent;hw=hflux
		
	velfrict = dsqrt(tau/ro)	
      if (l1 /= 0) Elatent = -Liv*ro*c_u*c_t*wr*(humair-humsurf) 

c CALCULATED sensible and latent heat fluxes are used in heat balance

	SHF=hflux
	LHF=Elatent
	
c MEASURED sensible and latent heat fluxes are used in heat balance
	!SHF = Hm
	!LHF = LEm
	
	if (snow==1) then
       Radbal = longwave-emissivity*sigma*(Tsurf+273.15)**4
      elseif (ice==0.and.water==1) then
	 Radbal = shortwave*(1-albedo)*(1-(1-sabs)*
     & exp(-extwat*0.5*ddz*h1))+longwave-emissivity*sigma*
     &(Tsurf+273.15)**4
	elseif (ice==1.and.snow==0) then
	 Radbal = shortwave*(1-albedo)*(1-exp(-extice*l1))+longwave-
     & emissivity*sigma*(Tsurf+273.15)**4
	else
       Radbal = shortwave*(1-albedo)+longwave-emissivity*sigma*
     &(Tsurf+273.15)**4  
	endif 
	
      SELECT CASE (surftyp)   
	 case (1)
	  dHdt=csoil(1)*rosoil(1)*dzs(1)/2.*(Tsoil2(1)-Tsoil1(1))/dt
	 case (2)
	  dHdt=ci*roi*ddz*l1/2.*(Ti2(1)-Ti1(1))/dt
	 case (3)
	  dHdt=cs(1)*dens(1)*dz(1)/2.*(Tsn(1)-T(1))/dt
       case (4)  
        dHdt=cw*row*ddz*h1/2.*(Tw2(1)-Tw1(1))/dt
	END SELECT

	surfrad = emissivity*sigma*(Tsurf+273.15)**4  
	botflux = -lamw(M)*(Tw2(M+1)-Tw2(M))/(ddz*h1)

      Bal=dHdt-(Radbal-SHF-LHF)+eflux
	
	END SUBROUTINE BAL_CALC


      SUBROUTINE T_diff(surf,Tsurf,dt)
      
!     T_diff calculates temperature profile in water, soil, ice and snow,
!     with known temparature at the surface 

	use numeric_params
	use phys_constants2
	use driving_params
	use arrays
      implicit none

	real(8) dt,Tsurf,lam2,lam1,dzmean,Hflow,ddz
	real(8), dimension(1:ms):: Tsn,lams,q,cs
	real(8), dimension(1:350):: a,b,c,d,Temp
	real(8) h1,l1,hs1,ls1
	real(8) AL,DLT,DVT,ALLL,DL,ALV,DV,Z,T,WL,WV,WI,dens
	real(8) dz
	real(8), external:: meltpnt
	integer(4) num,i,j,itop,surf
   

      common /snow_char/ Tsn,cs
	common /layers/ h1,l1,hs1,ls1
	common /SOILSOL/ AL(ML),DLT(ML),DVT(ML),ALLL(ML),DL(ML),
     & ALV(ML),DV(ML),Z(ML),T(ML),WL(ML),WV(ML),WI(ML),dens(ms)
	common /SOILDAT/ dz(ms),itop
	common /watericesnowarr/ lams,q
	data num /1/
	
	SAVE

	Hflow=0.
!     Meltpnt - Melting point temperature, C degrees 
!	Meltpnt=0.
	ddz=1./float(M)

	if (surf==1) then
	 if (snow==1) then
	  c(itop)=1.
	  b(itop)=0.
	  d(itop)=Tsurf
	  if (itop<=ms-2) then
	   call diff_coef(a,b,c,d,itop+1,ms-1,itop+1,ms-1,3,dt) 
	  endif
	  if (ice==1) then
!----------------ICE-SNOW INTERFACE
	   a(ms)=-(lams(ms-1)+lams(ms))/(2.*dz(ms-1))
	   b(ms)=-lami/(ddz*l1)
	   c(ms)=a(ms)+b(ms)-(cs(ms)*dens(ms)*dz(ms)/(2.*dt)+
     &   ci*roi*ddz*l1/(2.*dt))
	   d(ms)=-Ti1(1)*(cs(ms)*dens(ms)*dz(ms)/(2.*dt)+
     &   ci*roi*ddz*l1/(2.*dt))
!-----------------------------------
	   call diff_coef(a,b,c,d,2,M,ms+1,ms+M-1,2,dt) 
	   if (water==1) then
	    a(ms+M)=0.
	    c(ms+M)=1.
	    d(ms+M)=Meltpnt(Sal2(1))
	    call progonka2 (a,b,c,d,Temp,itop,ms+M)
          do i=itop, ms
	     Tsn(i)=Temp(i)
          enddo
	    do i=ms, ms+M
	     Ti2(i-ms+1)=Temp(i)
	    enddo 
         else
!----------------ICE-SOIL INTERFACE-------------------------
	    a(ms+M)=-lami/(ddz*l1)
	    b(ms+M)=-(lamsoil(1)+lamsoil(2))/(2*dzs(1))
	    c(ms+M)=a(ms+M)+b(ms+M)-(csoil(1)*rosoil(1)*dzs(1)/(2.*dt)+
     &    ci*roi*ddz*l1/(2.*dt))
	    d(ms+M)=-Ti1(M+1)*(csoil(1)*rosoil(1)*dzs(1)/(2.*dt)+
     &    ci*roi*ddz*l1/(2.*dt))
!-----------------------------------------------------------       
          call diff_coef(a,b,c,d,2,ns-1,ms+M+1,ms+M+ns-2,1,dt) 
	    c(ms+M+ns-1)=1.
	    a(ms+M+ns-1)=1.
	    d(ms+M+ns-1)=Hflow*dzs(ns-1)/(lamsoil(ns-1)+lamsoil(ns))
	    call progonka2 (a,b,c,d,Temp,itop,ms+M+ns-1)
	    do i=itop, ms
	     Tsn(i)=Temp(i)
          enddo
	    do i=ms, ms+M
	     Ti2(i-ms+1)=Temp(i)
	    enddo 
	    do i=ms+M, ms+M+ns-1
	     Tsoil2(i-ms-M+1)=Temp(i)
          enddo
         endif
        else
!-------------------SNOW-SOIL INTERFACE-------------------------
         a(ms)=-(lams(ms-1)+lams(ms))/(2.*dz(ms-1))
	   b(ms)=-(lamsoil(1)+lamsoil(2))/(2*dzs(1))
	   c(ms)=a(ms)+b(ms)-(cs(ms)*dens(ms)*dz(ms)/(2.*dt)+
     &   csoil(1)*rosoil(1)*dzs(1)/(2.*dt))
	   d(ms)=-Tsoil1(1)*(cs(ms)*dens(ms)*dz(ms)/(2.*dt)+
     &   csoil(1)*rosoil(1)*dzs(1)/(2.*dt))
!---------------------------------------------------------------
         call diff_coef(a,b,c,d,2,ns-1,ms+1,ms+ns-2,1,dt) 
	   c(ms+ns-1)=1.
	   a(ms+ns-1)=1.
	   d(ms+ns-1)=Hflow*dzs(ns-1)/(lamsoil(ns-1)+lamsoil(ns))
	   call progonka2 (a,b,c,d,Temp,itop,ms+ns-1)
         do i=itop, ms
	    Tsn(i)=Temp(i)
         enddo
	   do i=ms, ms+ns-1
	    Tsoil2(i-ms+1)=Temp(i)
         enddo
        endif
       elseif (ice==1) then
	  b(1)=0.
	  c(1)=1.
	  d(1)=Tsurf
        call diff_coef(a,b,c,d,2,M,2,M,2,dt)  
	  if (water==1) then
	   a(M+1)=0.
	   c(M+1)=1.
	   d(M+1)=Meltpnt(Sal2(1))
	   call progonka2 (a,b,c,d,Temp,1,M+1)
	   do i=1,M+1
	    Ti2(i)=Temp(i)
         enddo
        else
!----------------ICE-SOIL INTERFACE-------------------------
	   a(M+1)=-lami/(ddz*l1)
	   b(M+1)=-(lamsoil(1)+lamsoil(2))/(2*dzs(1))
	   c(M+1)=a(M+1)+b(M+1)-(csoil(1)*rosoil(1)*dzs(1)/(2.*dt)+
     &   ci*roi*ddz*l1/(2.*dt))
	   d(M+1)=-Ti1(M+1)*(csoil(1)*rosoil(1)*dzs(1)/(2.*dt)+
     &   ci*roi*ddz*l1/(2.*dt))
!-----------------------------------------------------------
         call diff_coef(a,b,c,d,2,ns-1,M+2,M+ns-1,1,dt) 
	   c(M+ns)=1.
	   a(M+ns)=1.
	   d(M+ns)=Hflow*dzs(ns-1)/(lamsoil(ns-1)+lamsoil(ns))
         call progonka2 (a,b,c,d,Temp,1,M+ns)
	   do i=1,M+1
	    Ti2(i)=Temp(i)
         enddo
	   do i=M+1,M+ns
	    Tsoil2(i-M)=Temp(i)
         enddo
	  endif
	 elseif (water==1) then  
	  b(1)=0.
	  c(1)=1.
	  d(1)=Tsurf
	  call diff_coef(a,b,c,d,2,M,2,M,4,dt) 
	  if (deepice==1) then
	   a(M+1)=0.
	   c(M+1)=1.
	   d(M+1)=Meltpnt(Sal2(M+1))
	   call progonka2(a,b,c,d,Temp,1,M+1) 
	   do i=1,M+1
	    Tw2(i)=Temp(i)
         enddo
        else
!---------------------WATER-SOIL INTERFACE------------------
	   a(M+1)=-lamw(M)/(ddz*h1)
	   b(M+1)=-(lamsoil(1)+lamsoil(2))/(2*dzs(1))
	   c(M+1)=a(M+1)+b(M+1)-(csoil(1)*rosoil(1)*dzs(1)/(2.*dt)+
     &   cw*row*ddz*h1/(2.*dt))
	   d(M+1)=-Tw1(M+1)*(csoil(1)*rosoil(1)*dzs(1)/(2.*dt)+
     &   cw*row*ddz*h1/(2.*dt))-SR(M)
!-----------------------------------------------------------
	   call diff_coef(a,b,c,d,2,ns-1,M+2,M+ns-1,1,dt)
	   c(M+ns)=1.
	   a(M+ns)=1.
	   d(M+ns)=Hflow*dzs(ns-1)/(lamsoil(ns-1)+lamsoil(ns))
	   call progonka2 (a,b,c,d,Temp,1,M+ns)
	   do i=1,M+1
	    Tw2(i)=Temp(i)
         enddo
	   do i=M+1,M+ns
	    Tsoil2(i-M)=Temp(i)
         enddo
        endif
	 else
	  b(1)=0.
	  c(1)=1.
	  d(1)=Tsurf  
        call diff_coef(a,b,c,d,2,ns-1,2,ns-1,1,dt) 
	  c(ns)=1.
	  a(ns)=1.
	  d(ns)=Hflow*dzs(ns-1)/(lamsoil(ns-1)+lamsoil(ns))
	  call progonka2 (a,b,c,d,Temp,1,ns)
	  do i=1,ns
	   Tsoil2(i)=Temp(i)
        enddo
       endif
	else   
       if (water==0) then
	  RETURN
       else
	  if (ice==1) then
	   b(1)=0.
	   c(1)=1.
	   d(1)=Meltpnt(Sal2(1))
         call diff_coef(a,b,c,d,2,M,2,M,4,dt)  
	   if (deepice==1) then
	    a(M+1)=0.
	    c(M+1)=1.
	    d(M+1)=Meltpnt(Sal2(M+1))
          call progonka2 (a,b,c,d,Temp,1,M+1)
          do i=1,M+1
	     Tw2(i)=Temp(i)
          enddo
          b(1)=0.
	    c(1)=1.
	    d(1)=Meltpnt(Sal2(M+1))
	    call diff_coef(a,b,c,d,2,M,2,M,5,dt)
!----------------DEEPICE-SOIL INTERFACE-------------------------
	    a(M+1)=-lami/(ddz*ls1)
	    b(M+1)=-(lamsoil(1)+lamsoil(2))/(2*dzs(1))
	    c(M+1)=a(M+1)+b(M+1)-(csoil(1)*rosoil(1)*dzs(1)/(2.*dt)+
     &    ci*roi*ddz*ls1/(2.*dt))
	    d(M+1)=-Tis1(M+1)*(csoil(1)*rosoil(1)*dzs(1)/(2.*dt)+
     &    ci*roi*ddz*ls1/(2.*dt))
!-----------------------------------------------------------	    
          call diff_coef(a,b,c,d,2,ns-1,M+2,M+ns-1,1,dt)  
	    c(M+ns)=1.
	    a(M+ns)=1.
	    d(M+ns)=Hflow*dzs(ns-1)/(lamsoil(ns-1)+lamsoil(ns))
	    call progonka2 (a,b,c,d,Temp,1,M+ns)
	    do i=1,M+1
	     Tis2(i)=Temp(i)
	    enddo
	    do i=M+1,M+ns
	     Tsoil2(i-M)=Temp(i)
	    enddo
         else
!---------------------WATER-SOIL INTERFACE------------------
         a(M+1)=-lamw(M)/(ddz*h1)
	   b(M+1)=-(lamsoil(1)+lamsoil(2))/(2*dzs(1))
	   c(M+1)=a(M+1)+b(M+1)-(csoil(1)*rosoil(1)*dzs(1)/(2.*dt)+
     &   cw*row*ddz*h1/(2.*dt))
	   d(M+1)=-Tw1(M+1)*(csoil(1)*rosoil(1)*dzs(1)/(2.*dt)+
     &   cw*row*ddz*h1/(2.*dt))
!-----------------------------------------------------------
         call diff_coef(a,b,c,d,2,ns-1,M+2,M+ns-1,1,dt)
	   c(M+ns)=1.
	   a(M+ns)=1.
	   d(M+ns)=Hflow*dzs(ns-1)/(lamsoil(ns-1)+lamsoil(ns))
	   call progonka2 (a,b,c,d,Temp,1,M+ns)
	   do i=1,M+1
	    Tw2(i)=Temp(i)
         enddo
	   do i=M+1,M+ns
	    Tsoil2(i-M)=Temp(i)
         enddo
        endif
	 elseif (deepice==1) then   	    
        b(1)=0.
	  c(1)=1.
	  d(1)=Meltpnt(Sal2(M+1))
	  call diff_coef(a,b,c,d,2,M,2,M,5,dt)
!----------------DEEPICE-SOIL INTERFACE-------------------------
	  a(M+1)=-lami/(ddz*ls1)
	  b(M+1)=-(lamsoil(1)+lamsoil(2))/(2*dzs(1))
	  c(M+1)=a(M+1)+b(M+1)-(csoil(1)*rosoil(1)*dzs(1)/(2.*dt)+
     &  ci*roi*ddz*ls1/(2.*dt))
	  d(M+1)=-Tis1(M+1)*(csoil(1)*rosoil(1)*dzs(1)/(2.*dt)+
     &  ci*roi*ddz*ls1/(2.*dt))
!-----------------------------------------------------------	    
        call diff_coef(a,b,c,d,2,ns-1,M+2,M+ns-1,1,dt)  
	  c(M+ns)=1.
	  a(M+ns)=1.
	  d(M+ns)=Hflow*dzs(ns-1)/(lamsoil(ns-1)+lamsoil(ns))
	  call progonka2 (a,b,c,d,Temp,1,M+ns)
	  do i=1,M+1
	   Tis2(i)=Temp(i)
	  enddo
	  do i=M+1,M+ns
	   Tsoil2(i-M)=Temp(i)
	  enddo
       endif
      endif
	endif

	END SUBROUTINE T_diff

	
	SUBROUTINE S_diff(dt)
	use arrays
	use driving_params
	implicit none
	
!     S_diff solves salinity (mineralization) diffusion equation	

      real(8) dt,Sflux1,ddz
      real(8), dimension(1:350):: a,b,c,d,Sal
      real(8) h1,l1,hs1,ls1
      
!     Sflux1 --- salinity flux at the bottom boundary       
!     Sflux0 --- salinity flux at the top boundary       

!     Sflux0 = -dh0/dt*Sal2(1)

      common /layers/ h1,l1,hs1,ls1
      
      ddz = 1./float(M)
      Sflux1 = 0.
      
      call diff_coef(a,b,c,d,2,M,2,M,6,dt)
      c(1)   = 1. + dhw0*ddz*h1/(dt*lamsal(1))
	b(1)   = 1.
	d(1)   = 0. ! Sflux0*ddz*h1/lamsal(1)
	c(M+1) = 1.
	a(M+1) = 1.
	d(M+1) = Sflux1*ddz*h1/lamsal(M)
	call progonka2 (a,b,c,d,Sal,1,M+1)
	Sal2(1:M+1)=Sal(1:M+1)
	
!	print*, 'Tw1', Tw1	
!	print*, 'Sal2', Sal2
!     print*, sum(Sal2(2:M)-Sal1(2:M)) 
!	read*
		
	END SUBROUTINE S_diff
	
	
!---------------------------------------------------------------------------------
      

	SUBROUTINE diff_coef(a,b,c,d,n0,n1,m0,m1,substr,dt)
	
!-------------------DEFINES COEFFICIENTS FOR SOLVING A THERMAL DIFFUSIVITY-----------------
!-------------------EQUATION BY FACTORIXATION METHOD----------------------------------------

      use numeric_params
	use phys_constants2	 
	use driving_params
	use arrays
	implicit none

	real(8), dimension(1:ms):: Tsn,lams,cs,q
	real(8), dimension(1:350):: a,b,c,d,Temp
	real(8) ddz,lam1,lam2,dzmean,dt
	real(8) h1,l1,hs1,ls1
	real(8) AL,DLT,DVT,ALLL,DL,ALV,DV,Z,T,WL,WV,WI,dens
	real(8) dz
	integer(4) i,j,n0,n1,m0,m1,substr,itop

      common /layers/ h1,l1,hs1,ls1
	common /snow_char/ Tsn,cs
	common /SOILSOL/ AL(ML),DLT(ML),DVT(ML),ALLL(ML),DL(ML),
     & ALV(ML),DV(ML),Z(ML),T(ML),WL(ML),WV(ML),WI(ML),dens(ms)
	common /SOILDAT/ dz(ms),itop
	common /watericesnowarr/ lams,q

	SAVE

      ddz=1./float(M)
      if (m1-m0/=n1-n0) then
	 print*, 'Error: diff_coef'
	 STOP
      endif

!     substr is a type of media : 1 - soil
!                                 2 - ice
!                                 3 - snow
!                                 4 - water
!                                 5 - deepice

!                                 6 - salinity in water

	SELECT CASE (substr)
!     CASE 1: SOIL
	 case (1)
       do i=m0,m1
	  j=i-m0+n0
	  lam2 = (lamsoil(j)+lamsoil(j+1))/2
	  lam1 = (lamsoil(j-1)+lamsoil(j))/2
	  dzmean = (dzs(j-1)+dzs(j))/2
	  a(i)=-lam1/dzs(j-1)
	  b(i)=-lam2/dzs(j)
	  c(i)=-(lam1/dzs(j-1)+lam2/dzs(j))-
     &   1./dt*(csoil(j)*rosoil(j)*dzmean)
	  d(i)=-Tsoil1(j)/dt*(csoil(j)*rosoil(j)*dzmean)
     	 enddo 
!      CASE 2: ICE
	 case (2)
       do i=m0,m1
	  j=i-m0+n0
	  c(i)=-(lami/ddz+lami/ddz)/l1-1./dt*ddz*l1*ci*roi
	  a(i)=-lami/(ddz*l1)+ci*roi*dhi*(j-1)*ddz/(2*dt)-
     &  ci*roi*dhi0/(2*dt)  
        b(i)=-lami/(ddz*l1)-ci*roi*dhi*(j-1)*ddz/(2*dt)+
     &  ci*roi*dhi0/(2*dt)  
	  d(i)=-Ti1(j)/dt*ddz*l1*ci*roi !dt*Erad*(1-sabs)*extwat/(cw*row)*
     	 enddo
!      CASE	3: SNOW  
	 case (3)
       do i=m0,m1
	  j=i-m0+n0
	  lam2 = (lams(j)+lams(j+1))/2
	  lam1 = (lams(j-1)+lams(j))/2
	  dzmean = (dz(j-1)+dz(j))/2
	  a(i)=-lam1/dz(j-1)
	  b(i)=-lam2/dz(j)
	  c(i)=-(lam1/dz(j-1)+lam2/dz(j))-
     &   1./dt*(cs(j)*dens(j)*dzmean)
	  d(i)=-T(j)/dt*(cs(j)*dens(j)*dzmean)
     	 enddo 
!      CASE 4: WATER
	 case (4)
	 do i=m0,m1
	  j=i-m0+n0
	  c(i)=-(lamw(j)/ddz+lamw(j-1)/ddz)/h1-1./dt*ddz*h1*cw*row
	  a(i)=-lamw(j-1)/(ddz*h1)+cw*row*dhw*(j-1)*ddz/(2*dt)-
     &  cw*row*dhw0/(2*dt)  
        b(i)=-lamw(j)/(ddz*h1)-cw*row*dhw*(j-1)*ddz/(2*dt)+
     &  cw*row*dhw0/(2*dt)  
	  d(i)=-Tw1(j)/dt*ddz*h1*cw*row-RS(j)*ddz*h1 !dt*Erad*(1-sabs)*extwat/(cw*row)*
     	 enddo 
!      CASE 5: DEEP ICE
       case (5)
	 do i=m0,m1
	  j=i-m0+n0
	  c(i)=-(lami/ddz+lami/ddz)/ls1-1./dt*ddz*ls1*ci*roi
	  a(i)=-lami/(ddz*ls1)+ci*roi*dls*(j-1)*ddz/(2*dt)-
     &  ci*roi*dls0/(2*dt)  
        b(i)=-lami/(ddz*ls1)-ci*roi*dls*(j-1)*ddz/(2*dt)+
     &  ci*roi*dls0/(2*dt)  
	  d(i)=-Tis1(j)/dt*ddz*ls1*ci*roi !dt*Erad*(1-sabs)*extwat/(cw*row)*
     	 enddo
!      CASE 6: SALINITY IN WATER
	 case (6)
	 do i=m0,m1
	  j=i-m0+n0
	  c(i)=-(lamsal(j)/ddz+lamsal(j-1)/ddz)/h1-1./dt*ddz*h1 
	  a(i)=-lamsal(j-1)/(ddz*h1)+dhw*(j-1)*ddz/(2*dt)-dhw0/(2*dt)  
        b(i)=-lamsal(j)/(ddz*h1)-dhw*(j-1)*ddz/(2*dt)+dhw0/(2*dt)  
	  d(i)=-Sal1(j)/dt*ddz*h1 !-RS(j)*ddz*h1 !dt*Erad*(1-sabs)*extwat/(cw*row)*
     	 enddo 
      END SELECT
  
	END SUBROUTINE diff_coef


      SUBROUTINE SURF_CHAR(surftyp)
	use driving_params
	use phys_constants2
	
!----------------DEFINES SURFACE CHARACTERISTICS:-------------------
!----------------roughness,emissivity,albedo,relative humidity,-----
!----------------coefficients in Magnus formula---------------------

	implicit none
      real(8) roughness,emissivity,albedo,aM,bM,velfrict,niu,relhums,
     & dirdif,sinh0
      integer(4) surftyp
      common /surface/ roughness,emissivity,albedo,aM,bM,relhums
       
	SAVE 
	 	
!     niu is air molecular viscosity
      niu = 1.5d-5 !1.007E-6
      select case (surftyp)
	 case (1)
        albedo = albedoofsoil
	  roughness = 0.05
	  aM = aMagw 
	  bM = bMagw
	  emissivity = emissivityofsoil
	  relhums = 0.5
       case (2) 
        albedo = albedoofice
	  emissivity = emissivityofice
	  aM = aMagi
	  bM = bMagi 
	  roughness = 0.00001
	  relhums = 0.7
	 case (3)
	  albedo = albedoofsnow
	  emissivity = emissivityofsnow
        aM = aMagi
	  bM = bMagi 
	  roughness = 0.001
	  relhums = 0.7
       case (4)
	  velfrict = dmax1(velfrict, 1.d-2)            
	  roughness = dmax1(0.111*niu/velfrict+
     &  0.0144*velfrict**2/g, 1.d-5) !for water 
	  aM = aMagw
	  bM = bMagw
	  if (varalb==0) then
	   albedo = albedoofwater
	  elseif (varalb==1) then
	   albedo = 0.05/(sinh0()+0.15)   
	   !albedo = (dirdif()*0.05/(sinh0()+0.15)+0.05)/
     &   !(1+dirdif())
	  endif
	  emissivity = emissivityofwater 
	  relhums = 1.
	end select	

	END SUBROUTINE SURF_CHAR

      REAL(8) FUNCTION dirdif()
      implicit none
	real(8) cloud
	real(8) dirdif0,b,S0,shortwave,precip,sinh0
	common /cloud/ cloud
	common /atmos2/ shortwave,precip
!	data cloud /0./

	SAVE

      b=1./3.
	S0=1367.
	cloud = 0.
    
      dirdif0 = (shortwave-b*S0*sinh0())/(b*(S0*sinh0()-shortwave))
	dirdif = dirdif0*(1.-sngl(cloud))
	dirdif = dmax1(dirdif,0.d0)
      
	END FUNCTION dirdif

	REAL(8) FUNCTION sinh0()
	use driving_params
!     sinh0 is sine of solar height (or cosine of zenith angle)	
	
      implicit none
	
	real(8) delta,theta,pi
	real(8) hour,phi_rad,lam_rad
	real(8) year,month,day,nday
	
      common /ymdh/ year,month,day,hour

      pi=dble(4.*atan(1.))		

	nday = day+(month-1)*30
	delta = 23.5*pi/180.*dCOS(2*pi*(nday-173)/365.)
	theta = pi*(hour-12.)/12.
	phi_rad = phi_g*pi/180.
	!lam_rad = lam_g*pi/180.
      sinh0 = dSIN(phi_rad)*dSIN(delta)+
     & dCOS(phi_rad)*dCOS(delta)*dCOS(theta)
	sinh0=dmax1(sinh0,dble(0.d0)) 

	END FUNCTION sinh0

c_____________________________________________________________________________

	SUBROUTINE RichPBL(h,le,tau,Tsurf1,humsurf,z0)
	
!     RichPBL calculates heat fluxes according to Louis parameterization (Louis, 1079)	
	
	use atmos
	implicit none
      
	real(8),parameter:: gg = 9.80665
      real(8),parameter:: ra = 287.05
      real(8),parameter:: rv = 461.51
      real(8),parameter:: cp = 1005.46
      real(8),parameter:: cpv = 1869.46
      real(8),parameter:: xl = 2500800.
      real(8),parameter:: stefan = 5.6697d-08
      real(8),parameter:: xpi = 3.141592654
	real(8),parameter:: too = 273.16
      real(8),parameter:: xpoo = 100000.
      real(8),parameter:: vkarmn = .4
      real(8),parameter:: xli = 2834500.
      real(8),parameter:: tcdil = 86400.
      real(8),parameter:: rascp = ra / cp
      real(8),parameter:: xlsg = xl / gg
      real(8),parameter:: xlscp = xl / cp
      real(8),parameter:: cpsl = cp / xl
      real(8),parameter:: cpsg = cp / gg
      real(8),parameter:: gscp = gg / cp
      real(8),parameter:: unsg = 1. / gg
      real(8),parameter:: gsra = gg / ra
      real(8),parameter:: rasl = ra / xl
      real(8),parameter:: rascp2 = ra / (2.*cp)
      real(8),parameter:: rasrv = ra / rv
      real(8),parameter:: rvsra = rv / ra
      real(8),parameter:: etv = rvsra - 1.
      real(8),parameter:: ecph = cpv / cp - 1.
      real(8),parameter:: xlsrv = xl / rv
      real(8),parameter:: unscp = 1. / cp
      real(8),parameter:: cpvmcp = cpv - cp
      real(8),parameter:: etvq = 1. - rasrv
      real(8),parameter:: xlf = xli - xl
      real(8),parameter:: xliscp = xli / cp
      real(8),parameter:: xlisg = xli / gg

	real(8) z0hz0,a,b,esatsurf,Tsurf1,humsurf,z,xmu,xfh,ts,h,le,leg,
     & urel,vrel,ta,qa,ps,ua,z0h,z0,hu,veg,ztvi,rhoa,zdepl,zua,ztvis,Ri,
     & valnon,zeps,zch,zsta,iyra,zdi,rra,zdim,cdm,zds,chstar2,qs,ph2,
     & cdh,cmstar2,pm2,zdsm,xhu,tau,wind10,u,v
	common /wind/ wind10,urel,vrel,u,v

	SAVE
      
      ta=tempair+273.15
	qa=humair
	ps=pressure
	ua=dsqrt(urel**2+vrel**2)
	ts=Tsurf1+273.15
	qs=humsurf

      z0hz0=0.1
	z0h=z0*z0hz0

c
c**********************************************************************
c           delta calculatlon
c
!      delta=0.
!      if (veg.gt.0.) delta=(wr(li)/wrmax)**(2./3.)
c
c**********************************************************************
c                         hv calculation
c                         ______________
c
      ! hv = 1. - dmax1(0.d0,dsign(1.d0,qsat(ts(li), ps) - qa)) * rs *
      !&      (1. - delta) / (rra + rs )

	  hu=1. !air above water is saturated
	  veg=0. !no vegetation on the lake
c
c?    rra ‚ a resistencia aerodinamica = 1/ChVa  -volta
c        a ser calculada … frente (correctamente

c        Implicit resolution of the surface temperature equation
c        linearization of Ts

        ztvi = ta * (1. + etv * qa)
        rhoa = ps / (ra * ztvi)
        !qs = (hu * (1. - veg) + veg * hv) * qsat (ts(li),ps) +
     &  !     (1. - hv) * veg * qa
c
c**********************************************************************
c                CMWF, 1981/11/25-27, pp. 59-79)
c**********************************************************************
c
        zdepl = 0.
        !write(*,*) 'nhli9904zref',zref
        zua=zref + zdepl
        ztvis = ts * (1. + etv * qs)
c______________________________________________________________________
c       calculo do numero de richardson e da resistencia aerodinamica
c
      !write(*,*) 'zref',zref
       ua=dmax1(1.d0,ua)
       ri=2.*gg*zua*zua*(ztvi-ztvis+gscp*zref)/(ztvi+ztvis)/(ua*ua)/zref
c
c       calculo do comprimento de Monin -Obukhov,de ustar,tstar
c
        valnon=-999.
c        call calclmon(z,zO,zOh,ua,ztvi,ztvls,xlmon,ustar,tstar,valnon)
        zeps= 1.e-6
        !write(* ,*) 'zua,z0',zua,z0
        zch=(vkarmn/dlog(dabs(zua/z0)))**2
        !write(*,*) 'zch',zch
        zsta=ri*ua*ua
c
c  iyra escolhe entre a nova e a velha versao da rotina yra
c         de NP89. Quando iyra (ficheiro de entrada(?!)) e maior
c         do que 0 utiliza-se a nova versao que destingue entre
c         os comprimentos de rugosidade
c
        iyra=1
        if(iyra.eq.0) then
c
          if(ri.lt.0.) then
           zdi=1./(ua+75.*zch*dsqrt(-zref*zsta/z0h))
           rra=zch*(ua-15.*zsta*zdi)
           zdim=1./(1.+75.*zch*dsqrt(-zref*ri/z0))
           cdm=zch*(1.-10.*ri*zdim)*ua
          else
           zds=dsqrt(ua * ua + 5. * zsta + zeps)
           rra = zch*ua/(1.+15.*zsta*zds/ua/ua/ua)
           zdim=(1.+5.*ri)**2.
           cdm=zch*ua/zdim
          endif
c
        else
c
c nova rotina yramng
          ! write(*,*) 'estez0h',z0h
          !write(*,*) 'z,zref',z,zref
           z=zref
           xmu=dlog(z0/z0h)
           xfh=dlog(z/z0)/dlog(z/z0h)
          ! write(*,*) '--------------------'
          ! write(*,*) 'z0,z0h',z0,z0h
           if(ri.lt.0.) then
                zdi=1./(ua+chstar2(xmu)*zch*15.*(z/z0h)**
     &		  ph2(xmu)*xfh
     s              *dsqrt(-zsta))
                rra=zch*(ua-15.*zsta*zdi)*xfh
c
c alteracoes para o calculo do coeficiente de dragg_______________
c        cdm=Cd*uvsuf
c
                cdh=rra
                zdim=1./(ua+cmstar2(xmu)*zch*10.*(z/z0)**
     &		  pm2(xmu)
     s              *dsqrt(-zsta))
                cdm=zch*(ua-10.*zsta*zdim)
c_______________________________________________________________________
           else
                zds = dsqrt(ua * ua + 5. * zsta + zeps)
                rra = zch*ua/(1.+15.*zsta*zds/ua/ua/ua)*xfh
c
                cdh=rra
                zdsm= dsqrt(ua * ua + 5. * zsta + zeps)
                cdm = zch*ua/(1.+10.*zsta/zdsm/ua)
	
           endif
c
        endif
c
        rra=1./rra
c----------------------------------------------------------------------
c
!      zrsra = rhoa / rra
!      if (iclay.ge.0) then
!         ct = 1. / ((1. - veg) / cso + veg / cv)
!         hv=1.-dmax1(0.d0,dsign(1.d0,qsat(ts(li), ps) - qa)) * rs*
!     &    (1. - delta) / (rra + rs)
c
!         za=1. / dt + ct * (4. * emis * stefan * (ts(li)**3) +
!     &      zrsra*xl*dqsat(ts(li),ps)*(veg*hv+(1.-veg)*hu*xhu)
!     &      +zrsra*cp)+2.*xpi/tcdil
!         zb=1./dt+ ct*(3. * emis * stefan * (ts(li)** 3) +
!     &      zrsra*xl*dqsat(ts(li),ps)*(veg*hv+(1.-veg)*hu*xhu))
!         zc=2.*xpi*t2(li)/tcdil+ct*(zrsra*cp*ta+rg*
!     &      (1.-alb)+emis*rat-zrsra*xl*(veg*hv*(qsat(ts(li),
!     &      ps)-qa)+(1.-veg)*xhu*(hu*qsat(ts(li),ps)-qa)))
c
!         ts(lf) = (ts(li) * zb + zc) / za
c
c       resolution of t2 equation
c
!         t2(lf) = (t2(li) + dt * ts(lf) / tcdil) / (1. + dt / tcdil)
!      else
!          ts(lf)=tsfunc
!          t2(lf)=t2(li)
!      endif
c
c**********************************************************************
c
      xhu=1.
      !hqs2=hu*qsat(ts(lf),ps)
      !rn = rg * (1. - alb) + emis * (rat - stefan * (ts(lf)**4))
      h = rhoa * cp * (ts - ta) / rra
      leg = rhoa * xl*(1.-veg) *xhu*(qs - qa)/ rra 
      !lev = rhoa * xl * veg * hv * (qsat(ts(lf), ps) - qa) / rra
      !zzhv = dmax1(0.d0,dsign(1.d0,qsat(ts(li), ps) - qa))
      !letr = zzhv * (1. - delta) * rhoa * xl * veg * (qsat(ts(lf), ps)
     &!       - qa) / (rra + rs)
	le = leg !+ lev
	tau = rhoa*cdm*ua
      !hw=h;xlew=le;cdmw=cdm

      return
	END SUBROUTINE RichPBL
      	
      function chstar2(xmu)

      implicit real*8(a-h,o-z)
      chstar2=3.2165+4.3431*xmu+.536*xmu*xmu-.0781*xmu*xmu*xmu
      return
      end

      function cmstar2(xmu)

      implicit real*8(a-h,o-z)
      cmstar2=6.8741+2.6933*xmu+.3601*xmu*xmu-.0154*xmu*xmu*xmu
      return
      end

      function ph2(xmu)

      implicit real*8(a-h,o-z)
      ph2    =0.5802-0.1571*xmu+.0327*xmu*xmu-.0026*xmu*xmu*xmu
      return
      end
c
c******************************************************************
      function pm2(xmu)
c****************************************************************
c
      implicit real*8(a-h,o-z)
      pm2    =0.5233-0.0815*xmu+.0135*xmu*xmu-.001*xmu*xmu*xmu
      return
      end
	
      SUBROUTINE dragvl(bx2,bix2,itdrag)
      
!     dragvl calculates exchange coefficients in aerodynamic formulas for surface
!     sensible heat, latent heat and momentum fluxes at the surface
!     following Businger-Dayer interpolation formulas, Beljaars parameterization and others      
      
      real bx(7),bix(11)
	real(8) bix2(11),bx2(7)
	INTEGER KL,ML,MS,nt,nspinmax,Num_Soil,Num_Veget

      common /p_drag/ vkc,anu,z0min,chc,p_m,p_t,q_m,q_t,
     *                alpha_m,alpha_t,beta_m,beta_t,gamma_m,gamma_t,
     *                zeta_st,a_m,a_t,b_m,b_t
	real vkc,anu,z0min,chc,p_m,p_t,q_m,q_t,
     *                alpha_m,alpha_t,beta_m,beta_t,gamma_m,gamma_t,
     *                zeta_st,a_m,a_t,b_m,b_t,grav
	real u_a,tpot_a,tpot_s,q_a,q_s,z,z0m,z0t,tvir_a,
     &tvir_s,dtvir,alam,c_u,c_t,ust,hfl,re,xx,zeta,zeta0m,
     &zeta0t,smo,psi_m,psi_t,phi_m,phi_t
      INTEGER iter,itop,itdrag
     
	   
      grav = 9.81
	
	bx=sngl(bx2)  

      u_a = amax1(bx(1),0.1) ! speed 
      tpot_a = bx(2) ! potential temperature - air
      tpot_s = bx(3) ! potential temperature - surface
      q_a = bx(4)
      q_s = bx(5)
      z = bx(6)
      if(bx(7) .gt. 0.) then
		z0m = bx(7)
      else
		z0m = z0min
      end if
      tvir_a = tpot_a * (1. + 0.61 * q_a) ! virtual temp - air
      tvir_s = tpot_s * (1. + 0.61 * q_s) ! virtual temp - surf
      dtvir = tvir_a - tvir_s  ! for smo
      alam = grav/(0.5*(tpot_a + tpot_s))  ! for smo
      c_u = vkc/log(z/z0m) ! for smo
      c_t = c_u
      do iter = 1,itdrag
		ust = c_u * u_a  ! for smo
		hfl =  - c_u * c_t * u_a * dtvir ! for smo Hflux
		if(bx(7) .lt. 0.) z0m = amax1(chc*ust**2/grav, z0min)
		re = ust*z0m/anu
		if(re .le. 0.111) then
			xx = - 2.43
		else
			if(0.111 .lt. re .and. re .le. 16.3) then
				xx = 0.83*log(re) - 0.6
			else
				xx = 0.49 * re**0.45
			end if
		end if
		z0t = amax1(z0m*exp(-xx), z0min)
c         z0t = z0m
		if(abs(dtvir) .lt. 1.e-3) then
			zeta = 1.
			zeta0m = 1.
			zeta0t = 1.
		else
			smo = - ust**3/(alam*vkc*hfl) ! vkc = 0.4
			smo = sign(1.,smo)*amax1(abs(smo),1.e-3) !
			zeta = z/smo
			zeta0m = z0m/smo
			zeta0t = z0t/smo
		end if
		c_u = vkc/amax1( (log(z/z0m) - psi_m(zeta,zeta0m)), 0.1)
		c_t = vkc/amax1( (log(z/z0t) - psi_t(zeta,zeta0t)), 0.1)
      end do
CC         c_t = c_u
      bix(1) = zeta
      bix(2) = alam*z*dtvir/u_a**2
      bix(3) = re
      bix(4) = log(z0m/z0t)
      bix(5) = z0m
      bix(6) = z0t
      bix(7) = 0.
      bix(8) = c_u
      bix(9) = c_t
      bix(10) = vkc*ust*z/phi_m(zeta)
      bix(11) = phi_m(zeta) / phi_t(zeta)

	bix2=dble(bix)

      return
      END                       

      real function phi_m(zeta)
      common /p_drag/ vkc,anu,z0min,chc,p_m,p_t,q_m,q_t,
     *                alpha_m,alpha_t,beta_m,beta_t,gamma_m,gamma_t,
     *                zeta_st,a_m,a_t,b_m,b_t
	real vkc,anu,z0min,chc,p_m,p_t,q_m,q_t,
     *                alpha_m,alpha_t,beta_m,beta_t,gamma_m,gamma_t,
     *                zeta_st,a_m,a_t,b_m,b_t
      real zeta,y,a
	a=1.0

c*   universal function for momentum: businger et al. (1971) formulation
c*   for zeta .ge. zeta_st and large et al. (1994) formulation otherwise
c*   for zeta .ge. 0. : beljaars and holtslag (1991)
      if(zeta .ge. 0.) then
c         y = alpha_m + beta_m*zeta     !  businger et al.
		y = 1. + zeta*(a+0.667*(6.-0.35*zeta)*exp(-0.35*zeta)) ! b&h
      else
		if(zeta. ge. zeta_st) then
			y = (1. - gamma_m*zeta)**(p_m)
		else
			y = (a_m - b_m*zeta)**(q_m)
		end if
      end if
      phi_m = y

      return
      end

      real function phi_t(zeta)
      common /p_drag/ vkc,anu,z0min,chc,p_m,p_t,q_m,q_t,
     *                alpha_m,alpha_t,beta_m,beta_t,gamma_m,gamma_t,
     *                zeta_st,a_m,a_t,b_m,b_t
	real vkc,anu,z0min,chc,p_m,p_t,q_m,q_t,
     *                alpha_m,alpha_t,beta_m,beta_t,gamma_m,gamma_t,
     *                zeta_st,a_m,a_t,b_m,b_t
      real zeta,y,a
	a=1.0

c*   universal function for heat: businger et al. (1971) formulation
c*   for zeta .ge. zeta_st and large et al. (1994) formulation otherwise
c*   for zeta .ge. 0. : beljaars and holtslag (1991)
      if(zeta .ge. 0.) then
c         y = alpha_t + beta_t*zeta
		y = 1. + zeta*(
     &    a*sqrt(1.+2.*zeta*a/3.)+0.667*(6.-0.35*zeta)*exp(-0.35*zeta)) ! b&h
      else
		if(zeta .ge. zeta_st) then
			y = (1. - gamma_t*zeta)**(p_t)
		else
			y = (a_t - b_t*zeta)**(q_t)
		end if
      end if
      phi_t = y

      return
      end

      real function psi_m1(x)

	real x,pi
      pi = 3.1415926
c*   integrated universal function for momentum based on businger et al.
c*   (1971) formulation for zeta_st .le. zeta .le. 0.
c*    p_m = - 1./4. (paulson, 1970)
      psi_m1 = log(0.125*(1. + x**2)*(1. + x)**2) -2.*atan(x) + pi/2.
      return
      end

      real function psi_m2(a,x)
	real x,a
c*   integrated universal function for momentum based on large et al.
c*   (1994) formulation for zeta .le. zeta_st
c*    only for q_m = -1./3.
      psi_m2 = (0.5*log(x**2 + x + 1.)
     &         + sqrt(3.)*atan((2.*x + 1.)/sqrt(3.)))
     &         / sign(1.,a)*abs(a)**(1./3.)
      return
      end

      real function psi_m(zeta,zeta0)
      common /p_drag/ vkc,anu,z0min,chc,p_m,p_t,q_m,q_t,
     *                alpha_m,alpha_t,beta_m,beta_t,gamma_m,gamma_t,
     *                zeta_st,a_m,a_t,b_m,b_t
	real vkc,anu,z0min,chc,p_m,p_t,q_m,q_t,
     *                alpha_m,alpha_t,beta_m,beta_t,gamma_m,gamma_t,
     *                zeta_st,a_m,a_t,b_m,b_t
	real zeta,y,zeta0,x,x_st,psi_m2,psi_m1,x0,x0_st,
     &psi_t1,a
	a=1.0
c*   for zeta .ge. 0. : beljaars and holtslag (1991)
      if(zeta .ge. 0.) then
c         y = (1.-alpha_m)*log(zeta/zeta0) - beta_m*(zeta-zeta0)
		y = - (a*zeta+
     &      0.667*(zeta-(5./0.35))*exp(-0.35*zeta)+(0.667*5./0.35))  ! b&h
      else
		if(zeta .ge. zeta_st) then
			x = (1. - gamma_m*zeta)**(- p_m)
			y = psi_m1(x)
		else
			x0 = 1. - (b_m/a_m)*zeta
			x = sign(1.,x0)*abs(x0)**(- q_m)
			x0_st = 1. - (b_m/a_m)*zeta_st
			x_st = sign(1.,x0_st)*abs(x0_st)**(- q_m)
			y = psi_m1(zeta_st) + psi_m2(a_m,x) - psi_m2(a_m,x_st)
		end if
      end if
      psi_m = y
      return
      end

      real function psi_t1(x)
	real x
c*   integrated universal function for heat based on businger et al.
c*   (1971) formulation for zeta_st .le. zeta .le. 0.
c*    p_t = - 1./2. (paulson, 1970)
      psi_t1 = log(0.25*(1. + x)**2) !in (paulsen, 1970) log(0.25*(1. + x**2)**2)!
      return
      end

      real function psi_t(zeta,zeta0)
      common /p_drag/ vkc,anu,z0min,chc,p_m,p_t,q_m,q_t,
     *                alpha_m,alpha_t,beta_m,beta_t,gamma_m,gamma_t,
     *                zeta_st,a_m,a_t,b_m,b_t
	real vkc,anu,z0min,chc,p_m,p_t,q_m,q_t,
     *                alpha_m,alpha_t,beta_m,beta_t,gamma_m,gamma_t,
     *                zeta_st,a_m,a_t,b_m,b_t
	real zeta,zeta0,y,x,x0,psi_t1,x0_st,x_st,psi_m2,a
	a=1.0
c*   for zeta .ge. 0. : beljaars and holtslag (1991)
      if(zeta .ge. 0.) then
c         y = (1.-alpha_t)*log(zeta/zeta0) - beta_t*(zeta-zeta0)
		y = - (a*(1.+(2./3.)*zeta*a)**(1.5)
     &   + 0.667*(zeta-(5./0.35))*exp(-0.35*zeta)+((0.667*5./0.35)-1.))
      else
		if(zeta .ge. zeta_st) then
			x = (1. - gamma_t*zeta)**(- p_t)
			y = psi_t1(x)
		else
			x0 = 1. - (b_t/a_t)*zeta
			x = sign(1.,x0)*abs(x0)**(- q_t)
			x0_st = 1. - (b_t/a_t)*zeta_st
			x_st = sign(1.,x0_st)*abs(x0_st)**(- q_t)
			y = psi_t1(x_st) + psi_m2(a_t,x) - psi_m2(a_t,x_st)
		end if
      end if
      psi_t = y
      return
      end
      
      SUBROUTINE KT_eq(dt)
      
!     KT_eq calculates eddy diffusivity (ED) in water coloumn following parameterizations:

!     1) empirical profile of ED;
!     2) Semi-empirical formulations for ED;
!     2) k-epsilon parameterization for ED, including kolmogorov relation.
      
	use atmos
    	use driving_params
	use arrays
      implicit none

	real(8), allocatable::E2(:),eps2(:),eps3(:),
     & u2(:),v2(:),k1(:),k3(:),k4(:),k5(:),u3(:),
     & v3(:),C1(:),Re(:),row(:),Ri(:),uv(:),E_it1(:),
     & E_it2(:),C_num(:),E_it3(:),Feps(:),E_it21(:),
     & eps_it21(:),eps_it1(:),eps_it2(:),
     & l(:),k2_mid(:),k3_mid(:),k4_mid(:),k5_mid(:),
     & Um(:),Vm(:),E12(:),eps12(:),knum(:),ceps3(:),k2t(:),
     & dzeta1(:),Eeps(:),Ecent(:),Epscent(:),res_E(:),
     & res_eps(:),dresd(:,:,:,:),dres2dE(:),dres2deps(:) 
	real(8), allocatable:: WU(:),WV(:),GAMT(:),GAMU(:),GAMV(:),
     & TF(:),KCu(:),KLengu(:)
	integer(4), allocatable::num(:)
	real(8), dimension(350):: a,b,c,d,y
	real(8), dimension(350,2,2)::am,bm,cm
	real(8), dimension(350,2)::ym,dm
	real(8) alp(10)
!	real(8) sbix(11)
	real(8) nyear,nmonth,nday,nhour
	real(8) swdir,stau,dtdz0,KC0,wind10,kw,C10,co1,
     &kws,agw,urel,vrel,u,v,CL1,CL2,roughness,k21s,l0,b0,wr
	real(8) h1,l1,hs1,ls1
	real(8) pi,CE,ddz,row0,cw,g,lam_eps,Ck,ufr,minK,l_mid,
     & niu,lo,ACC,month,wdir,dt,kar,Cs,Cs1,dudz2dvdz2,
     & lam_E,Co,taux,tauy,lam_gen,dE_it,a0,a1,a2,a3,lam_T,drodt,
     & deps_it,tau_gr,tau_i,mi,ma,E_nom,drodt1,drodt2,dT4,E_mid,eps_mid,
     & dT_sb,kor,kor2,Cz_gr,coef1,coef2,coef3,rp,Cz_i,Ce_a,Ce_i,Ce_gr,
     & C2_a,C2_i,C2_gr,tau_sbl,tau_wav,tau_ix,tau_iy,tau_grx,tau_gry,
     & xx,yy,coef,kwe,k21,tau2,lm,ext_lamw,lamw0,day,hour,ceps1,ceps2,
     & E_nom_s,E_nom_gen,Conuv,al_it,kron,FS,FB,al,alft,GAMUN,zup,lamTM,
     & CEt,sigmae,sigmaeps,urels,vrels,minwind,cepsh,ceh,delta,delta0,
     & resid, resid0
        
	real(8) AG(5),GRADU,GRADV,GRADT,GAMVN,CON0,TFR,COGRTN,COGRUN,
     &COGRVN,DTDZH,CON1,CON2
	integer(4) ind_par,flag_num,iter,iter_t,
     & nl,i,j,k,M1,numEeps
     	integer(4) psi(10)
	logical indstab,ind_bound,iterat,firstcall
	character outpath*60,tp_name*10,numt*1
	
  	common /ymdh/ nyear,nmonth,nday,nhour
	common /wind/  wind10,urel,vrel,u,v
	common /out/   outpath
	common /layers/ h1,l1,hs1,ls1

	data firstcall /.true./
	data psi/6,3,7,2,5,4,8,1,0,0/ !/3,2,4,1,0,0,0,0,0,0/  

	SAVE

      if (firstcall) then
	if (turb_out==1) then
	 write (numt,'(i1)') Turbpar
	 open(2114,file=path(1:len_trim(path))//'results/'//
     &                         'err_progon.dat',  status = 'unknown')
	 open(2115,file=path(1:len_trim(path))//'results/'//'E-eps.dat', 
     &										    status = 'unknown')
	 open(2116,file=path(1:len_trim(path))//'results/'//'E-eps1.dat',
     &										    status = 'unknown')
	 open (2117,file=path(1:len_trim(path))//'results/'
     &  //'turb'//numt//'.dat',	    status = 'unknown')
  	 open (2118,file=path(1:len_trim(path))//'results/'
     &  //'temp'//numt//'.dat',	    status = 'unknown')
       open (2119,file=path(1:len_trim(path))//'results/'
     &  //'uv'//numt//'.dat',	    status = 'unknown')
       open (2120,file=path(1:len_trim(path))//'results/'
     &  //'KT'//numt//'.dat',	    status = 'unknown') 
	 open (2121,file=path(1:len_trim(path))//'results/'
     &  //'dE1'//numt//'.dat',	    status = 'unknown') 
	 open (2122,file=path(1:len_trim(path))//'results/'
     &  //'dE2'//numt//'.dat',	    status = 'unknown') 
	 open (2123,file=path(1:len_trim(path))//'results/
     &  err_Eeps_iter.dat',	        status = 'unknown')  
	endif
	allocate (WU(1:M+1),WV(1:M+1),GAMT(1:M+1),GAMU(1:M+1),
     & GAMV(1:M+1),TF(1:M+1),KCu(1:M+1),KLengu(1:M+1)) 
	allocate (num(1:M+1)) 
	allocate ( E2(1:M+1),eps2(1:M+1),eps3(1:M+1),
     & u2(1:M+1),v2(1:M+1),k1(1:M+1),k3(1:M+1),k4(1:M+1),k5(1:M+1),
     & u3(1:M+1),v3(1:M+1),C1(1:M+1),Re(1:M+1),row(1:M+1),Ri(1:M+1),
     & uv(1:M+1),E_it1(1:M+1),dres2dE(1:M+1),dres2deps(1:M+1),  
     & E_it2(1:M+1),C_num(1:M+1),E_it3(1:M+1),Feps(1:M+1),
     & E_it21(1:M+1),eps_it21(1:M+1),eps_it1(1:M+1),
     & eps_it2(1:M+1),l(1:M+1),k2_mid(1:M+1),k3_mid(1:M+1),
     & k4_mid(1:M+1),k5_mid(1:M+1),Um(1:M+1),Vm(1:M+1),E12(1:M+1),
     & eps12(1:M+1),knum(1:M+1),ceps3(1:M+1),k2t(1:M+1) )
	allocate ( dzeta1(1:M+1),Eeps(1:M+1),Ecent(1:M+1),
     & Epscent(1:M+1),res_E(1:M+1),res_eps(1:M+1),
     & dresd(2,2,1:M+1,1:M+1) ) 

	month = 30.*24.*60.*60.
	day = 24*60.*60.
	hour = 60.*60.
	ddz = 1./float(M)
	pi = 4.*dtan(1.d0)
!PHYSICAL CONSTANTS
      row0 = 1000.
	lamw0 = 0.561
	cw = 3990.
	g = 9.814 
	kar = 0.4
	kor = 1.26*10.**(-4) !(-4) !Koriolis parameter for lat=60 deg 
	niu = 1.007d-6
	roughness = 0.01 !water roughness	       
!DIMENSIONLESS CONSTANTS IN E-EPS PARAMETERIZATION, according to Goudsmit et al., 2002
      CE = 0.09 !0.09
	CEt = 0.072
	ceps1 = 1.44
	ceps2 = 1.92
	sigmae = 1.
	sigmaeps = 1.3
!	Ce_a=1. !4. !4.
!	Ce_i=1. !4. !1.
!	Ce_gr=1. !4. !3.
!	C2_a=Ce_a**(3./2.)  !*CE
!	C2_i=Ce_i**(3./2.)  !*CE
!	C2_gr=Ce_gr**(3./2.) !*CE 
	lam_E=CE/sigmaE ! 0.73
	lam_T=1. !0.05
	lam_eps=CE/sigmaeps ! 0.73
	lam_gen= 1. 
	kwe = 10. !100.
	lamTM = 0.14
!Co is according to Satyanarayana et al., 1999
	Co=1.9 !? 1.9 
!CONSTANTS IN CONVECTION PARAMETERIZATION
	KC0=0. !1.d+5 !500.
	dtdz0=2.
!CONSTANTS IN LENGMUIR PARAMETERIZATION
      CL1=60. !120./5. !120.
	CL2=2./5.
!CONSTANTS OF SIMOES'S PARAMETERIZATION
	Cs = 1.5  !0.15
	Cs1 = 100.	
	
	lo=0.4*ddz*h1*3.
	L0=0.1
		
!AUXILIARY PHYSICAL CONSTANTS	
	drodt1=-0.112
	drodt2=0.027 !-2.*drodt1 !0.02
	CONUV = 1.
!OTHER CONSTANTS
!	kws=0.1
		
!PARAMETERS OF NUMERICAL SCHEME
	AG(1)=0. !0.
	AG(2)=0. !0.
	AG(3)=0. !0.
	!AG(4)=1. !0.
	!AG(5)=1. !0.
	CON0 = 0.05
	CON1 = 0.11
      CON2 = 1.56
	al_it=0.2
	ma=20.5
	mi=1.
	ACC=1.d-10
	knum = 0.
	minwind = 1.

	numEeps = 1
	delta0 = 2.d-2
	
	do i=1,8 
	 alp(i)=2*(mi+ma-(ma-mi)*dcos((2*psi(i)-1)/16*pi))**(-1)
      enddo

	a0=800.969d-7
	a1=588.194d-7
	a2=811.465d-8
	a3=476.600d-10
      endif

	iterat=.true.
	kor2=kor*dt/2

	u=u1(1)
	v=v1(1)
	
	if (relwind==2) then
	 urel=uwind-u
	 vrel=vwind-v
	elseif (relwind==1) then
       urel=uwind
	 vrel=vwind
	else
	 print*, 'relwind=',relwind,' - incorrect meaning: STOP' 
	 stop
      endif
	
	if (dabs(urel)<minwind) then
	 urels=urel/dabs(urel)*minwind
!      elseif (dabs(urel)>1.) then
!       urels=urel/dabs(urel)*1.
	else
	 urels=urel
      endif

	if (dabs(vrel)<minwind) then
	 vrels=vrel/dabs(vrel)*minwind
!      elseif (dabs(vrel)>1.) then
!       vrels=vrel/dabs(vrel)*1.
	else   
       vrels=vrel
	endif
	
	wr=dmax1(dsqrt(urels**2+vrels**2),minwind) !,1.d-1)
	tau2=1.273*cdmw2/dmax1(dsqrt(urel**2+vrel**2),0.1d0)
     & *wr**2 !wind      
      ufr=dsqrt(tau2/row0)	
	do i = 1, M+1
!	 dT4=Tw1(i)-4.
!	 if (dT4/=0.) then 
!	  drodt=(0.5*drodt1*(dabs(dT4)+dT4)+0.5*drodt2*(dT4-dabs(dT4)))
!     &  /dT4
!	 row(i)=1000.+drodt*dT4 
!	 else
!	  row(i)=1000. 
!       endif
!	 T(i)=20.
       row(i) = row0*(1+a0+a1*Tw1(i)-a2*(Tw1(i)**2)+a3*(Tw1(i)**3))
	enddo  

      i=1
	KC=0
	KCu=0
! PARAMETERIZATION OF CONVECTION
	do while (i<M+1)
	 j=i
       c2:do while (j<M+1)
	    if (row(j)>row(j+1)) then
	     j=j+1
	    else
	     exit c2
          endif
         enddo c2
       if (i/=j) then
	  do k=i,j-1
	   KC(k)=dabs(Tw1(j)-Tw1(i))/((j-i)*ddz*h1)/dtdz0*KC0
        enddo
	  do k=i,j
	   KCu(k)=dabs(Tw1(j)-Tw1(i))/((j-i)*ddz*h1)/dtdz0*KC0/(cw*row0)
        enddo
       endif
	 i=j+1
      enddo     

	!PARAMETERIZATION OF WIND STRESS SPLIT-UP INTO TWO PARTS:
! tau=tau_sbl+tau_wave
      ! if (wind10>5.) then
      !  C10=1.26*10**(-3.)
      ! else
	!  C10=0.0044*wind10**(-1.15)
      ! endif
!a.PARAMETERIZATION USNG SIGNIFICANT WAVE HEIGHT (NB1, pp. 8)
	!co1=4.*kws/sqrt(C10*1000.*g) !1000 m is "average fetch" of lake Syrdakh
	!kw=min(max(co1*wind10,kws),0.75)
!b.PARAMETERIZATION USING WAVE AGE (agw) (NB1, pp. 9)
	!agw=0.14*(g*1000.)**(1./3.)*C10**(1./6.)*wind10**(-2./3.) 
	!kw=min(max(kws*1.15/agw,kws),0.75)
	!KW=0.75
	tau_wav=tau2 !*kw
	tau_sbl=tau2 !*(1.-kw)
!CALCULATION SHEZY COEFFICIENT
	rp=0.1 !height of roughness elements
     	if (l1/=0) then
	 Cz_gr=7.7*dsqrt(g)*(h1/2/rp)**(1./6.) !50. !ground Shezy coefficient
	else
	 Cz_gr=7.7*dsqrt(g)*(h1/rp)**(1./6.)
	endif
	rp=0.01 
      Cz_i=7.7*dsqrt(g)*(h1/2/rp)**(1./6.)

!FRICTION:GROUND
	tau_gr=row0*g*Cz_gr**(-2)*(u1(M)**2+v1(M)**2)
	tau_grx=row0*g*Cz_gr**(-2)*u1(M)*dsqrt((u1(M)**2+v1(M)**2)) !*0.005
	tau_gry=row0*g*Cz_gr**(-2)*v1(M)*dsqrt((u1(M)**2+v1(M)**2)) !*0.005
!FRICTION:ICE
	tau_i=row0*g*Cz_i**(-2)*(u1(1)**2+v1(1)**2)
	tau_ix=-row0*g*Cz_i**(-2)*u1(1)*dsqrt(u1(1)**2+v1(1)**2)
	tau_iy=row0*g*Cz_i**(-2)*v1(1)*dsqrt(u1(1)**2+v1(1)**2)
	coef1=-g*Cz_i**(-2)*dsqrt(u1(1)**2+v1(1)**2)
	coef2=g*Cz_gr**(-2)*dsqrt(u1(M)**2+v1(M)**2) !*0.005 !*row0
      
!________EDDY VISCOSITY PARAMETERIZATION__________________________________
      do i=1,M+1
	 uv(i)=dsqrt(u1(i)**2+v1(i)**2)
	enddo
	do i=1,M
!	 Ri(i)=dmin1(dmax1(g/row0/(((uv(i+1)-uv(i))/
!     & (ddz*h1))**2)*(row(i+1)-row(i))/(ddz*h1),-10.d0),10.d0)
      dudz2dvdz2 = dmax1( ( (u1(i+1)-u1(i))/(ddz*h1) )**2 +
     & ( (v1(i+1)-v1(i))/(ddz*h1) )**2, 1.d-1)
	 Ri(i) = g/row0*(row(i+1)-row(i))/(ddz*h1)/dudz2dvdz2
      enddo

      select case (Turbpar)
	 
!     1. "Empirical" parametrization: Stepanenko, Lykossov (2005)
	 case (1)
	 tp_name='SL_Empir' 
       k2(1) =(lamw0*10.+(wind*2./zref/20.)*(lamw0*1000.-lamw0*10.))/
     & (cw*row0)
	 ext_lamw = dlog(k2(1)*cw*row0/(lamw0*10.))/h1
	 do i=2,M+1
	  k2(i) = k2(1)*dexp(-ext_lamw*(i-1)*ddz*h1)+niu
	 enddo
	 k2(1)=k2(1)+niu

!     2. "E-epsilon" parameterization: k=E**2/eps with 
!         prognostic equations for E and eps
	 case (2)
	 do i=1,M+1
	  if (E1(i)<=0) E1(i)=10.**(-16)
	  if (eps1(i)<=0) eps1(i)=10.**(-18)
	 enddo
       do i = 1, M
	  k2(i)=dmax1(CE*E1(i)**2/(eps1(i)+ACC),1.d-6)+niu+knum(i)
	 end do

!     3. Nickuradze (NICK) formulation: Rodi (1993)
	 case (3)
	 tp_name='NICK'
	 do i=1,M
	  zup=h1-(i-0.5)*ddz*h1
	  lm=h1*(0.14-0.08*(1-zup/h1)**2-0.06*(1-zup/h1)**4)
	  k2(i)=lm**2*dabs((uv(i+1)-uv(i))/(ddz*h1))*dexp(-Cs*Ri(i))+niu
	 enddo
	
!     4. Parabolic (PARAB) formulation: Engelund (1976)
	 case (4)
	 tp_name='PARAB'
       do i=1,M
	  zup=h1-(i-0.5)*ddz*h1
	  k2(i)=kar*ufr*zup*(1-zup/h1)*dexp(-Cs*Ri(i))
     &  +niu
	 enddo 

!     5. W2 (used in Version 2 of CE-QUAL-W2 model): Cole and Buchak (1995)	
	 case (5)

!     6. W2N (W2 with mixing length of Nickuradze): Cole and Buchak (1995) and Rodi (1993)
	 case (6)
    
!     7. RNG (re-normalization group) formulation: Simoes (1998)
       case (7)
	 tp_name='RNG'
	 do i=1,M
	  zup=h1-(i-0.5)*ddz*h1
	  k2(i)=niu*(1+dmax1(3*kar*(zup/niu*ufr)**3*
     &  (1-zup/h1)**3-Cs1,0.d0))**(1./3.)*dexp(-Cs*Ri(i))+niu  
       enddo
	end select
	  
	 k2=k2+KCu !+KLengu
	
!(U,V) - EQUATION
       !xx=1+ddz**2*h1*dhw/(k2(1)*dt)
	 xx=0.5
	 yy=(ddz*h1)**2/(k2(1)*dt)
!UPPER BOUNDARY CONDITION FOR (U,V)-EQUATION.
!1 CASE:WATER SURFACE FREE FROM ICE (ATMOSPHERIC MOMENTUM FLUX ONLY)
	 if (l1==0) then
	  taux=urels/wr*tau_sbl/row0 !wdir+pi
	  tauy=vrels/wr*tau_sbl/row0 !wdir+pi
!	  if(dabs(xx+yy)<dabs(xx).and.debug_out==1)
!     &   write(2114,*) '1',nyear,nmonth,nday,time/month
	  do i=1,2
	   do j=1,2
	    cm(1,i,j)=kron(i,j)*(xx+yy) !xx
	    bm(1,i,j)=kron(i,j)*xx
         enddo
        enddo
	  cm(1,1,2)=-kor*(h1*ddz)**2/(2*k2(1))
        cm(1,2,1)=kor*(h1*ddz)**2/(2*k2(1))
	  dm(1,1)=(taux*ddz*h1/k2(1)+yy*u1(1))
	  dm(1,1)=dm(1,1)+(u1(2)-u1(1))*xx
     &  +kor*(h1*ddz)**2*v1(1)/(2*k2(1))
   	  dm(1,2)=(tauy*ddz*h1/k2(1)+yy*v1(1))
	  dm(1,2)=dm(1,2)+(v1(2)-v1(1))*xx
     &  -kor*(h1*ddz)**2*u1(1)/(2*k2(1))
       endif
!2 CASE: THIN ICE (COMBINATION OF ATMOSPHERIC MOMENTUM FLUX AND WATER-ICE FRICTION)
	 if (l1>0.and.l1<=L0) then 
	  b0=l1/L0
	  taux=urels/wr*tau_sbl/row0 !wdir+pi
	  tauy=vrels/wr*tau_sbl/row0 !wdir+pi
!	  if(dabs(xx+yy-2*b0*coef1*ddz*h1/k2(2))<dabs(xx)
!     &  .and.debug_out==1)
!     &   write(2114,*) '2',nyear,nmonth,nday,time/month 
	  do i=1,2
	   do j=1,2
	    cm(1,i,j)=kron(i,j)*(xx+yy-b0*coef1*ddz*h1/k2(1))
	    bm(1,i,j)=kron(i,j)*xx
	   enddo
        enddo
	  cm(1,1,2)=-kor*(h1*ddz)**2/(2*k2(1))
        cm(1,2,1)=kor*(h1*ddz)**2/(2*k2(1))
	  dm(1,1)=yy*u1(1) !-tau_ix*ddz*h1/k2(1)
        dm(1,1)=dm(1,1)+(u1(2)-u1(1))*xx+2*taux*ddz*h1/k2(1)*(1-b0)
     &  +kor*(h1*ddz)**2*v1(1)/(2*k2(1))
	  dm(1,2)=yy*v1(1) !-tau_iy*ddz*h1/k2(1)
        dm(1,2)=dm(1,2)+(v1(2)-v1(1))*xx+2*tauy*ddz*h1/k2(1)*(1-b0)
     &  -kor*(h1*ddz)**2*u1(1)/(2*k2(1))
	! do i=1,2
	!  do j=1,2
	!   cm(1,i,j)=dble(kron(i,j))
	!   bm(1,i,j)=0
	!  enddo
      ! enddo
	! dm(1,1)=0
	! dm(1,2)=0
       endif
!3 CASE: THICK ICE (WATER-ICE FRICTION ONLY)
	 if (l1>L0) then
!	  if(dabs(xx+yy-2*coef1*ddz*h1/k2(2))<dabs(xx).and.debug_out==1)
!     &   write(2114,*) '3',nyear,nmonth,nday,time/month 
	  do i=1,2
	   do j=1,2
	    cm(1,i,j)=kron(i,j)*(xx+yy-coef1*ddz*h1/k2(1))
	    bm(1,i,j)=kron(i,j)*xx
	   enddo
        enddo
	  cm(1,1,2)=-kor*(h1*ddz)**2/(2*k2(1))
        cm(1,2,1)=kor*(h1*ddz)**2/(2*k2(1))
	  dm(1,1)=yy*u1(1) !-tau_ix*ddz*h1/k2(1)
        dm(1,1)=dm(1,1)+(u1(2)-u1(1))*xx
     &  +kor*(h1*ddz)**2*v1(1)/(2*k2(1))
	  dm(1,2)=yy*v1(1) !-tau_iy*ddz*h1/k2(1)
        dm(1,2)=dm(1,2)+(v1(2)-v1(1))*xx
     &  -kor*(h1*ddz)**2*u1(1)/(2*k2(1))
	 ! do i=1,2
	 !  do j=1,2
	 !   cm(1,i,j)=dble(kron(i,j))
	 !   bm(1,i,j)=0
	 !  enddo
       ! enddo
	 ! dm(1,1)=0
	 ! dm(1,2)=0
	 endif
	 xx=0.5*(1-ddz*h1*dhw/(k2(M)*dt))
	 yy=(ddz*h1)**2/(k2(M)*dt)
!	 if(dabs(xx+yy+2*coef2*ddz*h1/k2(M))<dabs(xx).and.debug_out==1)
!     & write(2114,*) '3',nyear,nmonth,nday,time/month
	  do i=1,2
	   do j=1,2
	    !if (i==1) then
	    ! coef=coef2*(u1(M))/dabs(u1(M))
          !else
	    ! coef=coef2*(v1(M))/dabs(v1(M))
          !endif
	    cm(M+1,i,j)=kron(i,j)*(xx+yy+coef2*ddz*h1/k2(M))
	    am(M+1,i,j)=kron(i,j)*xx
	   enddo
        enddo
	  cm(M+1,1,2)=-kor*(h1*ddz)**2/(2*k2(M))
        cm(M+1,2,1)=kor*(h1*ddz)**2/(2*k2(M))
	  dm(M+1,1)=yy*u1(M+1) !-tau_grx*ddz*h1/k2(M)
	  dm(M+1,1)=dm(M+1,1)-(u1(M+1)-u1(M))*xx
     &  +kor*(h1*ddz)**2*v1(M+1)/(2*k2(M))
	  dm(M+1,2)=yy*v1(M+1) !-tau_gry*ddz*h1/k2(M)
	  dm(M+1,2)=dm(M+1,2)-(v1(M+1)-v1(M))*xx
     &  -kor*(h1*ddz)**2*u1(M+1)/(2*k2(M))
	 ! do i=1,2
	 !  do j=1,2
	 !   cm(M,i,j)=dble(kron(i,j))
	 !   am(M,i,j)=0
	 !  enddo
       ! enddo
	 ! dm(M,1)=0
	 ! dm(M,2)=0
	 do k=2,M
	  do i=1,2
	   do j=1,2
	    !a(i)=(i-1)*dhw/(2*h1)-k2_mid(i-1)*dt/(h1**2*ddz**2)
	    !b(i)=-(i-1)*dhw/(2*h1)-k2_mid(i)*dt/(h1**2*ddz**2)
	    !c(i)=a(i)+b(i)-1
	    !d(i)=-u1(i)
	    am(k,i,j)=kron(i,j)*((k-1)*dhw/(4*h1)-k2(k-1)*dt/
     &    (2*h1**2*ddz**2))
	    bm(k,i,j)=-kron(i,j)*((k-1)*dhw/(4*h1)+k2(k)*dt/
     &    (2*h1**2*ddz**2))
	   enddo
	  enddo
	  cm(k,1,1)=-(k2(k)+k2(k-1))*dt/(2*h1**2*ddz**2)-1
	  cm(k,1,2)=kor2 !Krank-Nickolson
!       cm(k,1,2)=kor2*2. !implicit scheme
!	  cm(k,1,2)=0. !explicit scheme
	  cm(k,2,1)=-cm(k,1,2)
	  cm(k,2,2)=cm(k,1,1)
	  dm(k,1)=-u1(k)-kor2*v1(k) !Krank-Nickolson
	  dm(k,1)=dm(k,1)-dt*(2*h1**2*ddz**2)**(-1)*(k2(k)*(u1(k+1)-
     &  u1(k)))
	  dm(k,1)=dm(k,1)+dt*(2*h1**2*ddz**2)**(-1)*(k2(k-1)*(u1(k)-
     &  u1(k-1)))
	  dm(k,1)=dm(k,1)-(k-1)*dhw*(4*h1)**(-1)*(u1(k+1)-u1(k-1))
	  dm(k,2)=-v1(k)+kor2*u1(k) !Krank-Nickolson
	  dm(k,2)=dm(k,2)-dt*(2*h1**2*ddz**2)**(-1)*(k2(k)*(v1(k+1)-
     &  v1(k)))
	  dm(k,2)=dm(k,2)+dt*(2*h1**2*ddz**2)**(-1)*(k2(k-1)*(v1(k)-
     &  v1(k-1)))
	  dm(k,2)=dm(k,2)-(k-1)*dhw*(4*h1)**(-1)*(v1(k+1)-v1(k-1))
!       dm(k,1)=-u1(k)  !implicit scheme
!	  dm(k,2)=-v1(k)  !implicit scheme
!	  dm(k,1)=-u1(k)-kor2*2*v1(k)  !explicit scheme
!	  dm(k,2)=-v1(k)+kor2*2*u1(k)  !explicit scheme
	 enddo
!      ind_bound=.false.
!	 call ind_stab_fact_db (a,b,c,1,M+1,indstab,ind_bound)
!	 if (indstab==.false.) then
!	  do i=2,M
!	   a(i)=-k2(i)*dt/(h1**2*ddz**2)
!	   b(i)=-k2(i+1)*dt/(h1**2*ddz**2)
!	  enddo
!	 endif  
       call MATRIXPROGONKA(am,bm,cm,dm,ym,M+1)
	 do k=1,M+1
	  u2(k)=ym(k,1)
	  v2(k)=ym(k,2)
	 enddo
	 
!	 print*, (u2(1)-u1(1))/dt*(h1*ddz)**2/k2(1)
!     & -coef1*u2(1)*h1*ddz/k2(1)
!     & - 0.5*(u2(2)-u2(1))-0.5*(u1(2)-u1(1)) -     
!     & 0.5*dhw*h1*ddz/(k2(1)*dt)*(u2(2)-u2(1)) -
!     & 0.5*dhw*h1*ddz/(k2(1)*dt)*(u1(2)-u1(1)) -
!     & kor*(ddz*h1)**2/k2(1)*(v2(1)+v1(1))/2.	 
      
!      do i=2,M
!      print*, (v2(i)-v1(i))/dt
!     &   -0.5/h1/ddz*(k2(i)*(v2(i+1)-v2(i))/h1/ddz+
!     &   k2(i)*(v1(i+1)-v1(i))/h1/ddz -
!     &   k2(i-1)*(v2(i)-v2(i-1))/h1/ddz -
!     &   k2(i-1)*(v1(i)-v1(i-1))/h1/ddz)-
!     &   0.25*(i-1)/h1*dhw/dt*(v2(i+1)-v2(i-1))-
!     &   0.25*(i-1)/h1*dhw/dt*(v1(i+1)-v1(i-1))+
!     &   kor*(u1(i)+u2(i))/2.
!      enddo
	 
!	 print*, (u2(M+1)-u1(M+1))/dt*(h1*ddz)**2/k2(M)
!     & +coef2*u2(M+1)*h1*ddz/k2(M)
!     & + 0.5*(u2(M+1)-u2(M))+ 0.5*(u1(M+1)-u1(M)) -
!     & 0.5*dhw*h1*ddz/(k2(M)*dt)*(u2(M+1)-u2(M)) -
!     & 0.5*dhw*h1*ddz/(k2(M)*dt)*(u1(M+1)-u1(M)) -
!     & kor*(ddz*h1)**2/k2(M)*(v2(M+1)+v1(M+1))/2.	 
	 

      if (nmonth==7.and.nyear==977) then
	 do i=1,M
        write (2118,*) time, (i-1)*ddz*h1, Tw1(i)
	  write (2119,*) time, (i-1)*ddz*h1, dsqrt(u2(i)**2+v2(i)**2)
        write (2120,'(f12.1,f6.2,3f16.5)') 
     &   time, (i-0.5)*ddz*h1, KT(i), Ri(i), 
     &   row(i+1)-row(i)
	 enddo
	endif

	select case (Turbpar)
	 case (2)
	 case default
	 goto 1
      end select

	if (numEeps==1) then

!Krank-Nikolson numerical scheme for k-eps parameterization
!using simple iterations

	Um=(u2+u1)/2
      Vm=(v2+v1)/2
      
!SOLUTION OF 'E-EPS' PARAMETERIZATION EQUATIONS	
       E_it1=E1
	 eps_it1=eps1
	 iter_t=0
	 k2_mid=0.

	 AL=g/row0 !*drodt
	 DTDZH=(Tw1(2)-Tw1(1))/(h1*ddz)
	 !dT_sb=T(M+1)-T(1)	 
       !ALFT=(0.5*30.*(dabs(dT_sb)+dT_sb)+0.5*30.*(dT_sb-dabs(dT_sb)))
     & !/dT_sb
	 !if (l1/=0) ALFT=1.
	 ALFT=1. !1.
 		
10     E12=(E1+E_it1)/2
	 eps12=(eps1+eps_it1)/2
      		
	 do i=1,M
        !E_mid=(E_it1(i)+E_it1(i+1))/2.
	  E_mid=(E12(i)+E12(i+1))/2.
	  !eps_mid=(eps_it1(i)+eps_it1(i+1))/2. 
	  eps_mid=(eps12(i)+eps12(i+1))/2. 
        !k2_mid(i)=dmax1(CE*E_mid*E_mid/(eps_mid+ACC),1.E-6)
	  k2_mid(i)=k2_mid(i)*al_it+
     &  dmax1(E_mid*E_mid/(eps_mid+ACC),1.d-6)*(1-al_it)
	  k5_mid(i)=k2_mid(i)*lam_eps
	  k3_mid(i)=k2_mid(i)*lam_E
!	  k2(i)=dmax1(CE*E12(i)*E12(i)/(eps12(i)+ACC),1.d-6)
	  k2(i)=k2(i)*al_it+
     &  dmax1(E12(i)*E12(i)/(eps12(i)+ACC),1.d-6)*(1-al_it)
!	  k2t(i)=k2(i)/CE*CEt
	  k5(i)=k2(i)*lam_eps
	  k3(i)=k2(i)*lam_E
	  !l(i)=CE*E_it1(i)**(3./2.)/(eps_it1(i)+ACC)
	  l(i)=E12(i)**(3./2.)/(eps12(i)+ACC)
	  !TF(i)=dSQRT(E_it1(i))/(l(i)+ACC)
	  TF(i)=dSQRT(E12(i))/(l(i)+ACC)
	 enddo
	  !k2(M+1)=dmax1(CE*E_it1(M+1)*E_it1(M+1)/(eps_it1(M+1)+ACC),1.E-8)
	  !k2(M)=dmax1(CE*E12(M)*E12(M)/(eps12(M)+ACC),1.E-8)
	  !k5(M)=k2(M)
	  !k3(M)=k2(M)

	 do i=2,M
        WU(i)=(E_it1(i+1)-E_it1(i-1))*(U2(i)-U2(i-1))/(2.*h1**2*ddz**2)
        WV(i)=(E_it1(i+1)-E_it1(i-1))*(V2(i)-V2(i-1))/(2.*h1**2*ddz**2)
	 enddo
	 
       do i=2,M
	  GRADT = (Tw1(i+1)-Tw1(i-1))/(2.*h1*ddz)
        GRADU = (U2(i)-U2(i-1))/(h1*ddz)
        GRADV = (V2(i)-V2(i-1))/(h1*ddz)
        GAMUN = - CONUV*(WU(i+1)-WU(i-1))/(2.*h1*ddz)
        GAMVN = - CONUV*(WV(i+1)-WV(i-1))/(2.*h1*ddz)
	 COGRTN = DTDZH * 0.1
	 TFR=E_it1(i)/(l(i)*l(i)+ACC) + acc
	 COGRUN = GAMUN/TFR
       COGRVN = GAMVN/TFR 
	 if(gradT < 0.) then
        gamt(i) = 0.
        gamu(i) = 0.
        gamv(i) = 0.
	  gamt(i) = - AG(3)*cogrtn
        gamu(i) = - AG(1)*cogrun
        gamv(i) = - AG(2)*cogrvn
       else
	    GAMT(i) = - AG(3)*(AL*GRADT**2+CON0*TFR*COGRTN)
     &              /(AL*GRADT+CON0*TFR)
          GAMU(i) = - AG(1)*(AL*(CON2*GRADT+(CON2-1.)*GAMT(i))*GRADU +
     &                CON1*TFR*COGRUN)/(AL*GRADT+CON1*TFR)
          GAMV(i) = - AG(2)*(AL*(CON2*GRADT+(CON2-1.)*GAMT(i))*GRADV +
     &                CON1*TFR*COGRVN)/(AL*GRADT+CON1*TFR)
       end if
	 enddo
	 GAMT(1)=0.;GAMU(1)=0.;GAMV(1)=0.

	 do i=1,M
	  Gen(i)=((Um(i+1)-Um(i))*(Um(i+1)-Um(i)+2.*h1*ddz*GAMU(i))+
     &  (Vm(i+1)-Vm(i))*(Vm(i+1)-Vm(i)+2.*h1*ddz*GAMV(i)))/
     &  (h1**2*ddz**2)*CE
	  S(i)=-((row(i+1)-row(i))/(h1*ddz)+GAMT(i))*ALFT*AL*CEt
	  F(i)=Gen(i)+S(i)
	 enddo
      
!E-EQUATION SOLUTION      
	 if (l1==0) then 
	  FS=kwe*(tau2/row0)**(3./2.)
	 elseif (l1>0.and.l1<=L0) then
	  FS=((b0*tau_i+(1-b0)*tau2)/row0)**(3./2.)
	 elseif (l1>L0) then
	  FS=(tau_i/row0)**(3./2.)
       endif
	 xx=0.5*(k3_mid(1)*dt*(h1*ddz)**(-2)+0.5*dhw*h1**(-1))
	 yy=dt/(h1*ddz)
	 b(1)=-xx
	 c(1)=-(xx+1)-(dABS(S(1))-S(1))/(TF(1)+ACC)*dt/2-
     & TF(1)*dt	  
	 d(1)=-E1(1)-xx*(E1(2)-E1(1))-yy*FS-
     & k2(1)*(S(1)+dabs(S(1)))*dt/2-k2(1)*Gen(1)*dt !+eps_it1(i)*dt
	 
	 FB=-(tau_gr/row0)**(3./2.)
	 xx=0.5*(k3_mid(M-1)*dt*(h1*ddz)**(-2)-
     & dhw*(h1*ddz)**(-1)*(1-ddz/2.))
	 yy=dt/(h1*ddz)
	 a(M)=-xx
	 c(M)=-(xx+1)-(dABS(S(M))-S(M))/(TF(M)+ACC)*dt/2-
     & TF(M)*dt
	 d(M)=-E1(M)-xx*(E1(M-1)-E1(M))+yy*FB
     & -k2(M)*(S(M)+dabs(S(M)))*dt/2-k2(M)*Gen(M)*dt !+eps_it1(i)*dt
	 do i=2,M-1
	  a(i)=(i-0.5)*dhw/(4*h1)-k3_mid(i-1)*dt/(2*h1**2*ddz**2)
	  b(i)=-(i-0.5)*dhw/(4*h1)-k3_mid(i)*dt/(2*h1**2*ddz**2)
	  c(i)=a(i)+b(i)-1-(dABS(S(i))-S(i))/(TF(i)+ACC)*dt/2-
     &  TF(i)*dt	  
     	  d(i)=-E1(i)-k2(i)*(S(i)+dabs(S(i)))*dt/2-k2(i)*Gen(i)*dt !+eps_it1(i)*dt
	  d(i)=d(i)-dt*(2*h1**2*ddz**2)**(-1)*(k3_mid(i)*(E1(i+1)-
     &  E1(i)))
	  d(i)=d(i)+dt*(2*h1**2*ddz**2)**(-1)*(k3_mid(i-1)*(E1(i)-
     &  E1(i-1)))
	  d(i)=d(i)-(i-0.5)*dhw*(4*h1)**(-1)*(E1(i+1)-E1(i-1))
	 enddo
	 ind_bound=.false.
	 call ind_stab_fact_db (a,b,c,1,M,indstab,ind_bound)
	 if (indstab==.false.) then
	  do i=2,M-1
	   a(i)=-k3_mid(i-1)*dt/(2*h1**2*ddz**2)
	   b(i)=-k3_mid(i)*dt/(2*h1**2*ddz**2)
	  enddo
    	 endif 
       call progonka2_db (a,b,c,d,E_it2,1,M)
	 !E_it2(1)=E_it2(2)
	 E_it2(M+1)=0.

		      
	 call CUT_E_DB(E_it2,M,0.d0)

	 dE_it=maxval(dabs(E_it1-E_it2))

      !if (iterat) then
	! if (dE_it>1.E-10.and.iter_t<35) then !.or.deps_it>10.**(-13)
	  !E_it1=E_it1*al_it+E_it2*(1-al_it)
	!  E_it1=E_it1+alp(iter+1)*(E_it2-E_it1)
	! endif
      !endif

!C1 is given following Satyanarayana et al., 1999; Aupoix et al., 1989
       do i=1,M !+1
        Re(i)=ACC+(2*E_it1(i)/3.)**2/(eps_it1(i)*niu+ACC)
	  C1(i)=Co/(1+0.69*(2-Co)/dsqrt(Re(i)))
	  if (S(i)<0.) then
	   ceps3(i)=-0.4
        else
         ceps3(i)= 1.
	  endif
	 enddo

!EPS - EQUATION
20	 if (l1==0) then
	  FS=k5(1)*E_it1(1)**(3./2.)/kar*(roughness+0.2)**(-2) !roughness+ddz/2*h1
	 else
	  FS=k5(1)*E_it1(1)**(3./2.)/kar*(0.01+ddz/2*h1)**(-2) !0.01 - roughness of ice !0.01+ddz/2*h1
       endif 
	 xx=0.5*(k5_mid(1)*dt*(h1*ddz)**(-2)+0.5*dhw*h1**(-1))
	 yy=dt/(h1*ddz)
	 b(1)=-xx
	 c(1)=-(xx+1)-C1(1)*TF(1)*dt
     & +C1(1)*0.5*(-dABS(S(1))+S(1))/(TF(1)+ACC)*dt
	 d(1)=-eps1(1)-xx*(eps1(2)-eps1(1))-yy*FS-0.5*C1(1)*
     & (dABS(S(1))+S(1))*TF(1)*k2(1)*dt-
     & C1(1)*TF(1)*k2(1)*Gen(1)*dt !+eps_it1(i)*dt
       
	 FB=-k5(M)*E_it1(M)**(3./2.)/kar*(0.01+ddz/2*h1)**(-2) !0.01+ddz/2*h1
	 xx=0.5*(k3_mid(M-1)*dt*(h1*ddz)**(-2)-
     & dhw*(h1*ddz)**(-1)*(1-ddz/2.))
	 yy=dt/(h1*ddz)
	 a(M)=-xx
	 c(M)=-(xx+1)-C1(M)*TF(M)*dt
     & +C1(M)*0.5*(-dABS(S(M))+S(M))/(TF(M)+ACC)*dt
	 d(M)=-eps1(M)-xx*(eps1(M-1)-eps1(M))+yy*FB-0.5*C1(M)*
     & (dABS(S(M))+S(M))*TF(M)*k2(M)*dt-
     & C1(M)*TF(M)*k2(M)*Gen(M)*dt
	 do i=2,M-1
	  a(i)=-k5_mid(i-1)*dt/(2*h1**2*ddz**2)+(i-0.5)*dhw/(4*h1)
	  b(i)=-k5_mid(i)*dt/(2*h1**2*ddz**2)-(i-0.5)*dhw/(4*h1)
        c(i)=a(i)+b(i)-1-C1(i)*TF(i)*dt
     &  +C1(i)*0.5*(-dABS(S(i))+S(i))/(TF(i)+ACC)*dt
	  d(i)=-eps1(i)-0.5*C1(i)*(dABS(S(i))+S(i))*TF(i)*k2(i)*dt-
     &   C1(i)*TF(i)*k2(i)*Gen(i)*dt
	  d(i)=d(i)-dt*(2*h1**2*ddz**2)**(-1)*(k5_mid(i)*(eps1(i+1)-
     &  eps1(i)))
	  d(i)=d(i)+dt*(2*h1**2*ddz**2)**(-1)*(k5_mid(i-1)*(eps1(i)-
     &  eps1(i-1)))
	  d(i)=d(i)-(i-0.5)*dhw*(4*h1)**(-1)*(eps1(i+1)-eps1(i-1))
       enddo
	 ind_bound=.false.
	 call ind_stab_fact_db (a,b,c,1,M,indstab,ind_bound)
	 if (indstab==.false.) then
	  do i=2,M-1
	   a(i)=-k5_mid(i-1)*dt/(2*h1**2*ddz**2)
	   b(i)=-k5_mid(i)*dt/(2*h1**2*ddz**2)
	  enddo
    	 endif  
	
	 call progonka2_db (a,b,c,d,eps_it2,1,M)
	 eps_it2(M+1)=0.
	 
       call CUT_EPS_DB(eps_it2,M,0.d0)
	
	 deps_it=maxval(dabs(eps_it1-eps_it2))
	
	 if (iterat) then
        if ((dE_it>1.E-10.or.deps_it>1.E-13).and.iter_t<35) then !.or.deps_it>10.**(-13)
	   E_it1=E_it1*al_it+E_it2*(1-al_it)
	   eps_it1=eps_it1*al_it+eps_it2*(1-al_it)
         !eps_it1=eps_it1+alp(iter+1)*(eps_it2-eps_it1)
	   !E_it1=E_it1+alp(iter+1)*(E_it2-E_it1)
	   !E_it1=E_it1+alp(iter+1)*(E_it2-E_it1)
	   iter=iter+1
	   iter_t=iter_t + 1
	   if (iter==8) iter=0
	   goto 10
	  else
	   if (dE_it<1.E-6.and.deps_it<1.E-7) then
	    E2=E_it2
	    eps2=eps_it2
	    iter=0
	    iter_t=0
	    nl=2
	   else
	    eps_it1=eps1
	    E_it1=E1
	    iterat=.false.
	    nl=1
	    goto 10
         endif
        endif
       else
        E2=E_it2
	  eps2=eps_it2
	  iter=0
	  iter_t=0
	 endif
	 if (dE_it<1.E-10.or.deps_it<1.E-13) nl=3
	 
!END OF ITERATIONAL PROCESS

      elseif (numEeps==2) then
	
	print*, 'This method is not operational at the moment'
	stop

      AL = g/row0
      do i=1,M
	 Gen(i)=( (U2(i+1)-U2(i))*(U2(i+1)-U2(i))+
     & (V2(i+1)-V2(i))*(V2(i+1)-V2(i)) )/
     & (h1**2*ddz**2)*CE
	 S(i)=-( (row(i+1)-row(i))/(h1*ddz) )*AL*CEt
	 Re(i)=ACC+(2*E1(i)/3.)**2/(eps1(i)*niu+ACC)
	 C1(i)=Co/(1+0.69*(2-Co)/dsqrt(Re(i)))
	 if (i==1.or.i==M) then
	  dzeta1(i) = (i-0.5)*dhw/(h1*dt)
       else
	  dzeta1(i) = (i-0.5)*dhw/(2*h1*dt)
	 endif
!	  if (S(i)<0.) then
!	   ceps3(i)=-0.4
!        else
!         ceps3(i)= 1.
!	  endif
      enddo

	cepsh = lam_eps/(ddz**2*h1**2)
	ceh = lam_E/(ddz**2*h1**2)
	if (l1==0) then 
	 FS=kwe*(tau2/row0)**(3./2.)
	 xx=ddz*h1*kar*(roughness+0.5*ddz*h1)**2
	elseif (l1>0.and.l1<=L0) then
	 FS=((b0*tau_i+(1-b0)*tau2)/row0)**(3./2.)
	 xx=ddz*h1*kar*(0.01+0.5*ddz*h1)**2
	elseif (l1>L0) then
	 FS=(tau_i/row0)**(3./2.)
	 xx=ddz*h1*kar*(0.01+0.5*ddz*h1)**2
      endif
	FB=-(tau_gr/row0)**(3./2.)
	yy=ddz*h1*kar*(0.01+0.5*ddz*h1)**2
      
	eps1(M+1)=0.
	E1(M+1)=0.
	E_it2(M+1)=0.
	eps_it2(M+1)=0.
	res_E(M+1)=0.
	res_eps(M+1)=0.
	dres2dE(M+1)=0.
	dres2deps(M+1)=0.

      eps_it1 = eps1
	E_it1 = E1
      res_E(1:M) = 1.
	res_eps(1:M) = 1.
!	res_Ep = res_E
!	res_epsp = res_eps

	iter_t = 0
	
! Counter-gradient iterations to find the minimum of quadratic residence
! of E-eps equations
      do while (maxval(dabs(res_E))  >1.d-7 .or.
     &	      maxval(dabs(res_eps))>1.d-7 )

      iter_t = iter_t + 1
!	resid0 = resid

      do i=2,M 
       Eeps(i) = 0.5*(E_it1(i-1)+E_it1(i))**2/   
     &  (eps_it1(i-1)+eps_it1(i)+ACC)
	 Ecent(i) = E_it1(i-1)+E_it1(i)
	 Epscent(i) = eps_it1(i-1)+eps_it1(i) + ACC
	enddo
	 
	do i = 2, M-1
	 res_E(i)=(E_it1(i)-E1(i))/dt-ceh*
     & (E_it1(i+1)*Eeps(i+1)-E_it1(i)*(Eeps(i+1)+Eeps(i))
     & +E_it1(i-1)*Eeps(i))-dzeta1(i)*(E_it1(i+1)-E_it1(i-1))
     & -E_it1(i)**2*(Gen(i)+S(i))/(eps_it1(i)+ACC)+eps_it1(i)
	 res_eps(i)=(eps_it1(i)-eps1(i))/dt-cepsh*
     & (eps_it1(i+1)*Eeps(i+1)-eps_it1(i)*(Eeps(i+1)+Eeps(i))
     & +eps_it1(i-1)*Eeps(i))-dzeta1(i)*(eps_it1(i+1)-eps_it1(i-1))
     & -C1(i)*E_it1(i)*(Gen(i)+S(i))+C1(i)*eps_it1(i)**2/(E_it1(i)+ACC)
      enddo

	res_E(1)=(E_it1(1)-E1(1))/dt-ceh*
     & (Eeps(2)*(E_it1(2)-E_it1(1))+FS*ddz*h1/lam_E)
     & -dzeta1(1)*(E_it1(2)-E_it1(1))
     & -E_it1(1)**2*(Gen(1)+S(1))/(eps_it1(1)+ACC)+eps_it1(1)  
      res_E(M)=(E_it1(M)-E1(M))/dt-ceh*
     & (-Eeps(M)*(E_it1(M)-E_it1(M-1))-FB*ddz*h1/lam_E)
     & -dzeta1(M)*(E_it1(M)-E_it1(M-1))
     & -E_it1(M)**2*(Gen(M)+S(M))/(eps_it1(M)+ACC)+eps_it1(M)  

	res_eps(1)=(eps_it1(1)-eps1(1))/dt-cepsh*
     & Eeps(2)*(eps_it1(2)-eps_it1(1))
     & -E_it1(1)**(7./2.)/(eps_it1(1)+ACC)/xx
     & -dzeta1(1)*(eps_it1(2)-eps_it1(1))
     & -C1(1)*E_it1(1)*(Gen(1)+S(1))+C1(1)*eps_it1(1)**2/(E_it1(1)+ACC)
      res_eps(M)=(eps_it1(M)-eps1(M))/dt+cepsh*
     & Eeps(M)*(eps_it1(M)-eps_it1(M-1))
     & -E_it1(M)**(7./2.)/(eps_it1(M)+ACC)/yy
     & -dzeta1(M)*(eps_it1(M)-eps_it1(M-1))
     & -C1(M)*E_it1(M)*(Gen(M)+S(M))+C1(M)*eps_it1(M)**2/(E_it1(M)+ACC)
      
      resid0 = dmax1(maxval(dabs(res_E)), maxval(dabs(res_eps)))

	do i=2,M-1
	 dresd(1,1,i,i+1) =
     & -ceh*( Eeps(i+1)+2*Eeps(i+1)/Ecent(i+1)*
     & (E_it1(i+1)-E_it1(i)) ) - dzeta1(i)

	 dresd(1,1,i,i) = 
     & 1./dt-ceh*(E_it1(i+1)*2*Eeps(i+1)/Ecent(i+1)-
     & ( Eeps(i+1)+Eeps(i))-2*E_it1(i)*
     & (Eeps(i+1)/Ecent(i+1)+Eeps(i)/Ecent(i))+
     & E_it1(i-1)*2*Eeps(i)/Ecent(i) )-
     & 2*E_it1(i)/(eps_it1(i)+ACC)*(Gen(i)+S(i))	 

	 dresd(1,1,i,i-1) = 
     & -ceh*( Eeps(i) + 2*Eeps(i)/Ecent(i)*
     & (E_it1(i-1)-E_it1(i)) ) + dzeta1(i)

	 dresd(1,2,i,i+1) = 
     & +ceh*( Eeps(i+1)/epscent(i+1)*(E_it1(i+1)-E_it1(i)) )

	 dresd(1,2,i,i) =
     & +ceh*( E_it1(i+1)*Eeps(i+1)/epscent(i+1)-E_it1(i)*
     & (Eeps(i+1)/epscent(i+1)+Eeps(i)/epscent(i))+
     & E_it1(i-1)*Eeps(i)/epscent(i) ) +
     & (E_it1(i)/(eps_it1(i)+ACC))**2*(Gen(i)+S(i)) + 1

	 dresd(1,2,i,i-1) =
     & +ceh*( Eeps(i)/epscent(i)*(E_it1(i-1)-E_it1(i)) )

	 dresd(2,2,i,i+1)=
     & -cepsh*( Eeps(i+1)+Eeps(i+1)/epscent(i+1)*
     & (eps_it1(i)-eps_it1(i+1)) ) - dzeta1(i)

	 dresd(2,2,i,i)=
     & 1./dt+cepsh*( eps_it1(i+1)*Eeps(i+1)/epscent(i+1)
     & +Eeps(i+1)+Eeps(i)-eps_it1(i)*(Eeps(i+1)/epscent(i+1)
     & +Eeps(i)/epscent(i))+Eeps(i)/epscent(i)*eps_it1(i-1) )
     & +2*C1(i)*eps_it1(i)/(E_it1(i)+ACC)

	 dresd(2,2,i,i-1)=
     & -cepsh*( Eeps(i)+Eeps(i)/epscent(i)*
     & (eps_it1(i)-eps_it1(i-1)) ) + dzeta1(i)

	 dresd(2,1,i,i+1)=
     & -cepsh*( 2*Eeps(i+1)/Ecent(i+1)*( eps_it1(i+1)-eps_it1(i) ) )

	 dresd(2,1,i,i)=
     & -cepsh*( 2*Eeps(i+1)/Ecent(i+1)*eps_it1(i+1)-
     & eps_it1(i)*(2*Eeps(i+1)/Ecent(i+1)+2*Eeps(i)/Ecent(i))
     & +eps_it1(i-1)*2*Eeps(i)/Ecent(i) ) -
     & C1(i)*( Gen(i)+S(i)+(eps_it1(i)/(E_it1(i)+ACC))**2 )

	 dresd(2,1,i,i-1)=
     & -cepsh*( 2*Eeps(i)/Ecent(i)*( eps_it1(i-1)-eps_it1(i) ) )
	enddo 

      dresd(1,1,1,1)=
     & 1/dt-ceh*( 2*Eeps(2)/Ecent(2)*(E_it1(2)-E_it1(1)) - Eeps(2) )
     & +dzeta1(1)-2*E_it1(1)/(eps_it1(1)+ACC)*(Gen(1)+S(1)) 
      dresd(1,1,1,2)=
     & -ceh*( 2*Eeps(2)/Ecent(2)*(E_it1(2)-E_it1(1)) + Eeps(2) ) 
     & -dzeta1(1)

	dresd(1,2,1,1)=
     & -ceh*( -Eeps(2)/epscent(2)*(E_it1(2)-E_it1(1)) )
     & +(E_it1(1)/(eps_it1(1)+ACC))**2*(Gen(1)+S(1)) + 1
	dresd(1,2,1,2)=
     & -ceh*( -Eeps(2)/epscent(2)*(E_it1(2)-E_it1(1)) )	

      dresd(1,1,M,M)=
     & 1./dt-ceh*( -2*Eeps(M)/Ecent(M)*(E_it1(M)-E_it1(M-1))-Eeps(M) )
     & - dzeta1(M)-2*E_it1(M)/(eps_it1(M)+ACC)*(Gen(M)+S(M))
      dresd(1,1,M,M-1)=
     & -ceh*( -2*Eeps(M)/Ecent(M)*(E_it1(M)-E_it1(M-1))+Eeps(M) ) 
     & + dzeta1(M)

	dresd(1,2,M,M)=
     & -ceh*( +Eeps(M)/epscent(M)*(E_it1(M)-E_it1(M-1)) )
     & +(E_it1(M)/(eps_it1(M)+ACC))**2*(Gen(M)+S(M)) + 1	
	dresd(1,2,M,M-1)=
     & -ceh*( +Eeps(M)/epscent(M)*(E_it1(M)-E_it1(M-1)) )	
        
      dresd(2,1,1,1)=
     & -cepsh*(eps_it1(2)-eps_it1(1))*2*Eeps(2)/Ecent(2)-
     & 7./2.*E_it1(1)**2.5/xx/(eps_it1(1)+ACC)-
     & C1(1)*(Gen(1)+S(1)+(eps_it1(1)/(E_it1(1)+ACC))**2)   
!      print*,'dresd', !-cepsh*(eps_it1(2)-eps_it1(1))*2*Eeps(2)/Ecent(2),
!     & 7./2.*E_it1(1)**2.5/xx/eps_it1(1),
!     & C1(1)*(Gen(1)+S(1)),
!     & eps_it1(1),E_it1(1)   
      dresd(2,1,1,2)=
     & -cepsh*2*Eeps(2)/Ecent(2)*(eps_it1(2)-eps_it1(1))  

	dresd(2,2,1,1)=
     & 1/dt+cepsh*Eeps(2)/epscent(2)*(eps_it1(2)-eps_it1(1))+
     & cepsh*Eeps(2)+E_it1(1)**(7./2.)/(eps_it1(1)+ACC)**2/xx+
     & dzeta1(1)+2*C1(1)*eps_it1(1)/(E_it1(1)+ACC)  	
	dresd(2,2,1,2)=
     & cepsh*Eeps(2)/epscent(2)*(eps_it1(2)-eps_it1(1))-
     & cepsh*Eeps(2)-dzeta1(1)

      dresd(2,1,M,M)=
     & cepsh*2*Eeps(M)/Ecent(M)*(eps_it1(M)-eps_it1(M-1))-
     & 7./2.*E_it1(M)**2.5/yy/(eps_it1(M)+ACC)-
     & C1(M)*(Gen(M)+S(M)+(eps_it1(M)/(E_it1(M)+ACC))**2)
      dresd(2,1,M,M-1)=
     & cepsh*2*Eeps(M)/Ecent(M)*(eps_it1(M)-eps_it1(M-1)) 
     
	dresd(2,2,M,M)=
     & 1/dt-cepsh*Eeps(M)/epscent(M)*( eps_it1(M)-eps_it1(M-1) )
     & +cepsh*Eeps(M)+E_it1(M)**(7./2.)/(eps_it1(M)+ACC)**2/yy-dzeta1(M)
     & +2*C1(M)*eps_it1(M)/(E_it1(M)+ACC)
	dresd(2,2,M,M-1)=
     & -cepsh*Eeps(M)/epscent(M)*( eps_it1(M)-eps_it1(M-1) )
     & -cepsh*Eeps(M)+dzeta1(M)	 

	
	dres2dE(1) = res_E  (1) * dresd(1,1,1,1) +
     &             res_E  (2) * dresd(1,1,2,1) +
     &             res_eps(1) * dresd(2,1,1,1) +
     &             res_eps(2) * dresd(2,1,2,1) 

	dres2dE(M) = res_E  (M)   * dresd(1,1,M,  M) +
     &             res_E  (M-1) * dresd(1,1,M-1,M) +
     &             res_eps(M)   * dresd(2,1,M,  M) +
     &             res_eps(M-1) * dresd(2,1,M-1,M)

      dres2deps(1) = res_E  (1) * dresd(1,2,1,1) +
     &               res_E  (2) * dresd(1,2,2,1) +
     &               res_eps(1) * dresd(2,2,1,1) +
     &               res_eps(2) * dresd(2,2,2,1)

	dres2deps(M) = res_E(M)    * dresd(1,2,M,  M) +
     &               res_E(M-1)  * dresd(1,2,M-1,M) +
     &               res_eps(M)  * dresd(2,2,M,  M) +
     &               res_eps(M-1)* dresd(2,2,M-1,M)

	do i=2, M-1
	dres2dE(i) =   res_eps(i-1) * dresd(2,1,i-1,i) +
     &               res_eps(i)   * dresd(2,1,i,  i) +
     &               res_eps(i+1) * dresd(2,1,i+1,i) +  
     &               res_E  (i-1) * dresd(1,1,i-1,i) +
     &               res_E  (i)   * dresd(1,1,i,  i) +
     &               res_E  (i+1) * dresd(1,1,i+1,i) 

	dres2deps(i) = res_eps(i-1) * dresd(2,2,i-1,i) +
     &               res_eps(i)   * dresd(2,2,i,  i) +
     &               res_eps(i+1) * dresd(2,2,i+1,i) +  
     &               res_E  (i-1) * dresd(1,2,i-1,i) +
     &               res_E  (i)   * dresd(1,2,i,  i) +
     &               res_E  (i+1) * dresd(1,2,i+1,i) 
	enddo

	dres2dE = dres2dE*2
	dres2deps = dres2deps*2

	if (iter_t==1) then
       delta = delta0  
      else
!      delta = delta0*(1.+3.*dmin1(resid/resid0,1.d0))
       delta = delta0*(1+float(iter_t)/1.)
!	 if (resid0<1.2d-6) delta=delta0 
!      delta = delta0
	endif

	iterat = .true.

	do while (iterat)
	
	E_it2 = E_it1 - 1.d-7*delta*dres2dE/maxval(dabs(dres2dE))
	eps_it2 = eps_it1 - 1.d-7*delta*dres2deps/maxval(dabs(dres2deps))
      
      do i=2,M 
       Eeps(i) = 0.5*(E_it2(i-1)+E_it2(i))**2/   
     &  (eps_it2(i-1)+eps_it2(i)+ACC)
	 Ecent(i) = E_it2(i-1)+E_it2(i)
	 Epscent(i) = eps_it2(i-1)+eps_it2(i)+ACC
	enddo

	do i = 2, M-1
	 res_E(i)=(E_it2(i)-E1(i))/dt-ceh*
     & (E_it2(i+1)*Eeps(i+1)-E_it2(i)*(Eeps(i+1)+Eeps(i))
     & +E_it2(i-1)*Eeps(i))-dzeta1(i)*(E_it2(i+1)-E_it2(i-1))
     & -E_it2(i)**2*(Gen(i)+S(i))/(eps_it2(i)+ACC)+eps_it2(i)
	 res_eps(i)=(eps_it2(i)-eps1(i))/dt-cepsh*
     & (eps_it2(i+1)*Eeps(i+1)-eps_it2(i)*(Eeps(i+1)+Eeps(i))
     & +eps_it2(i-1)*Eeps(i))-dzeta1(i)*(eps_it2(i+1)-eps_it2(i-1))
     & -C1(i)*E_it2(i)*(Gen(i)+S(i))+C1(i)*eps_it2(i)**2/(E_it2(i)+ACC)
      enddo

	res_E(1)=(E_it2(1)-E1(1))/dt-ceh*
     & (Eeps(2)*(E_it2(2)-E_it2(1))+FS*ddz*h1/lam_E)
     & -dzeta1(1)*(E_it2(2)-E_it2(1))
     & -E_it2(1)**2*(Gen(1)+S(1))/(eps_it2(1)+ACC)+eps_it2(1)  
      res_E(M)=(E_it2(M)-E1(M))/dt-ceh*
     & (-Eeps(M)*(E_it2(M)-E_it2(M-1))-FB*ddz*h1/lam_E)
     & -dzeta1(M)*(E_it2(M)-E_it2(M-1))
     & -E_it2(M)**2*(Gen(M)+S(M))/(eps_it2(M)+ACC)+eps_it2(M)  

	res_eps(1)=(eps_it2(1)-eps1(1))/dt-cepsh*
     & Eeps(2)*(eps_it2(2)-eps_it2(1))
     & -E_it2(1)**(7./2.)/(eps_it2(1)+ACC)/xx
     & -dzeta1(1)*(eps_it2(2)-eps_it2(1))
     & -C1(1)*E_it2(1)*(Gen(1)+S(1))+C1(1)*eps_it2(1)**2/(E_it2(1)+ACC)
      res_eps(M)=(eps_it2(M)-eps1(M))/dt+cepsh*
     & Eeps(M)*(eps_it2(M)-eps_it2(M-1))
     & -E_it2(M)**(7./2.)/(eps_it2(M)+ACC)/yy
     & -dzeta1(M)*(eps_it2(M)-eps_it2(M-1))
     & -C1(M)*E_it2(M)*(Gen(M)+S(M))+C1(M)*eps_it2(M)**2/(E_it2(M)+ACC)

      resid = dmax1(maxval(dabs(res_E)), maxval(dabs(res_eps)))	

      if (resid>resid0) then
	 delta=delta/2.
      else
	 iterat = .false.
	endif

      if (iter_t>10000) then
	 print*, 'resid=', resid, resid0
	 read*
	endif
	 	
	enddo

!	E_it2 = (E_it2+E_it1)/2.
!	Eps_it2 = (Eps_it2+Eps_it1)/2.

	    	
	call CUT_EPS_DB(eps_it2, M, 1.d-15) !1.d-11
	call CUT_E_DB(E_it2, M, 1.d-11) !1.d-9

!     print*, "The residual is"
!	print*, maxval(dabs(res_E)), maxval(dabs(res_eps)),iter_t !maxval(dabs(E_it2-E_it1))
!	if (iter_t>25000) then
!	 print*, (dsign(1.d0,dres2dE(i)), i=1,M )
!      endif
!	print*, maxval(dabs(dres2dE)), maxval(dabs(dres2deps))
!	print*, maxval(dabs(res_eps)), maxval(dabs(eps_it2-eps_it1))
!	print*, !dres2dE(1:M), dres2deps(1:M)
!     &	      res_E  (1) * dresd(1,1,1,1),
!     &        res_E  (2) * dresd(1,1,2,1),
!     &        res_eps(1) * dresd(2,1,1,1),
!     &        res_eps(2) * dresd(2,1,2,1)
!      read*
!      if (maxval(dabs(res_E))>10.) then
!       print*, 'The iterations are diverging', time/dt
!       print*, maxval(dabs(res_E)), resid, resid0
!       stop      
!	 read*
!	endif

      if (resid0==resid) then
	 if (maxval(dabs(res_E))>10.) then
        print*, 'The iterations are diverging', time/dt
        print*, maxval(dabs(res_E)), resid, resid0
        stop      
!	  read*
	 endif
	 if (resid>1.d-6) write(2123,*) resid	 
	 exit
      endif

      if (iter_t>1000) print*, iter_t, resid

	E_it1   = E_it2
	eps_it1 = eps_it2
!	res_Ep = res_E
!	res_epsp = res_eps

      enddo

	E2 = E_it2
	eps2 = eps_it2

      else

	 print*, 'No numerical method associated with numEeps=', numEeps
	 stop

      endif
	

      E12=(E1+E2)/2
	eps12=(eps1+eps2)/2
	do i=1,M
	 !k2(i)=dmin1(Ck*E2(i)**2/(eps3(i)+ACC),0.01)
	 !E_mid=(E12(i)+E12(i+1))/2.
	 !E_mid=(E2(i)+E2(i+1))/2.
	 !Eps_mid=(Eps12(i)+Eps12(i+1))/2.
       !Eps_mid=(Eps2(i)+Eps2(i+1))/2.
	 !KT(i)=lam_T*CE*E2(i)**2/(eps2(i)+ACC)*row0*cw
	  KT(i)=lam_T*dmax1(CEt*E2(i)**2/(eps2(i)+ACC),1.d-6)*row0*cw
	 !KT(i)=lam_T*dmax1(CE*E_mid**2/(eps_mid+ACC),1.E-6)*row0*cw
	enddo
	
!OUTPUT
      if (turb_out==1) then
	 !if (time/hour+12.-3.*24.>0.and.time/hour+12.-3.*24.<12.) then
	 do i=1,M+1
	  write (2115,7) time/hour+12.-3*24, (i-0.5)*ddz*h1, E2(i)*1.E9,
     &  eps2(i)*1.E9, KT(i),tau2*1.E3,Tw1(i)
	 enddo
	 !endif
       write (2116,8) time/hour+12.-3*24, E2(1)*1.E9,
     &  eps2(1)*1.E9, KT(1)
	endif
7 	format (f10.4,f8.2,2f12.1,3f8.2)      
8     format (f10.4,2f12.1,f8.2)      
      			
	 
      E1=E2
	eps1=eps2
1     u1=u2
	v1=v2
	
	select case (Turbpar)
	 case (2)
	 case (1)
	 do i=1,M+1
	  KT(i)=(k2(i)-niu)*cw*row0 !*0.01
       enddo
	 case default
	  KT=lamTM*k2*cw*row0 !*0.01
        if (turb_out==1) then
	   if (firstcall) write (2117,'(a)') tp_name
	   do i=1,M+1
          write (2117,'(F8.2, F5.2, E10.3, F8.1, F5.2, F6.2)')
     &    time/hour+12.-3*24, (i-0.5)*ddz*h1, Ri(i),
     &    KT(i),(i-1)*ddz*h1,uv(i)
	   enddo
	  endif
      endselect 
     
      if (firstcall) firstcall=.false.

	END SUBROUTINE KT_eq


      SUBROUTINE soilforlake(dt,h1,l1,ls1,hs1,dhwfsoil,a,b,c,d,Temp,
     & num)
     
!     SOILFORLAKE calculates profiles of temperature, liquid and solid water 
!     content in the soil under a lake     
     
	use phys_constants2 
	use driving_params
	use arrays
      implicit none

	      
	real(8), allocatable:: Tsoil4(:), Tsoil3(:), 
     & wl2(:), wl3(:), wl4(:), wl5(:),  
     & wi2(:), wi3(:), lammoist(:), rosdry(:), filtr(:)
	real(8), dimension(1:350):: a,b,c,d,Temp  
	real(8) Hfhigh,Hflow,dt,lam1,lam2,dzmean,ma,mi,pi,
     & wlmax_i,bb1,arg,psi,pf,cr,wflow,wfhigh,lammoist1,lammoist2,h1,
     & l1,dhwfsoil,surplus,lack,Potphenergy,watavailtofreeze,unfrwat,
     & penet,wlmax_i1,wlmax_i2,Tsoil_1,Tsoil_2,filtr_low,filtr_high,
     & time_fut, ls1,evol,evol1,hs1,month,aa,bb,cc,cc1,Tf_old,
     & dhwfsoil1,dhwfsoil2,dhwfsoil3,dhwfsoil4,wl_max,max1,cv
	real(8) alp(10),psiit(10)
	integer(4) i, j, k, num,iter2,iter,iter3
	logical flux_calc,indstab,soil_ts_ext_iter,ts_iter_more,
     & firstcall

	common/ts_ext_iter/ Tf_old,soil_ts_ext_iter
	data psiit/6,3,7,2,5,4,8,1,0,0/ !/3,2,4,1,0,0,0,0,0,0/  
	data firstcall /.true./
      SAVE

!	open (133, file=dir_out//'\dhwfsoil.dat', status = 'unknown')	
!	 num - number of point on the top boundary of soil
     
	if (firstcall) then
	allocate (Tsoil4(1:ns+2), Tsoil3(1:ns+2), 
     & wl2(1:ns+2), wl3(1:ns+2), wl4(1:ns+2), wl5(1:ns+2),  
     & wi2(1:ns+2), wi3(1:ns+2), lammoist(1:ns+2), rosdry(1:ns+2),
     & filtr(1:ns+2))
	endif

!	lam = 0.00041*418
	pi = 3.141592654
	cr =  0.2*4180
	rosdry = 1200.
	month = 24*30*60*60
!PARAMETERS OF ITERATIONAL PROCESS
	ma=10.5 !5.5
	mi=4.5
	do i=1,8 
	 alp(i)=2*(mi+ma-(ma-mi)*dcos((2*psiit(i)-1)/16*pi))**(-1)
      enddo
			       
     	do i=1, ns
       wlmax_i = POR(i)*ROW/(rosdry(i)*(1-por(i))) 
	 BB1 = BH(i) + 2.
	 csoil(i) = cr+WL1(i)*CW+WI1(i)*CI
	 rosoil(i) = rosdry(i)*(1-por(i))*(1+wi1(i)+wl1(i)) 
	 ARG = (WL1(i)+WI1(i))/wlmax_i
	 ARG = dmax1(ARG,1.d-2)
	 PSI = PSIMAX(i)*(ARG)**(-BH(i))
	 PF = dlog(-PSI*100.)/dLOG(10.d0)
	 IF(PF>5.1) THEN
	  lamsoil(i) = 4.1E-4*418.
	 ELSE
	  lamsoil(i) = dexp(-PF-2.7)*418.
	 END IF
	 wlmax_i = POR(i)*ROW/(rosdry(i)*(1-por(i))) - wi1(i)*row/roi
	 ARG = (WL1(i))/wlmax_i
	 ARG = dmax1(ARG,1.d-2)
	 lammoist(i) = dlmax(i)*ARG**BB1
	end do
      
      Hflow = 0
	
!SPLITTING-UP METHOD FOR TEMPERATURE EQUATION:
!STEP 1: HEAT DIFFUSION

      CALL T_SOLVER (dt)
      
!SPLITTING-UP METHOD FOR MOISTURE EQUATION:
!STEP 1: MOISTURE DIFFUSION
	
	Wflow = 0.
	dhwfsoil=0.
	dhwfsoil1=0.
	dhwfsoil2=0.
	dhwfsoil3=0.
	dhwfsoil4=0.
	
      if ((l1 /= 0.and. h1 == 0.).or.ls1/=0.or.(hs1/=0.and.l1==0.and.
     &h1==0)) then
	 Wfhigh = 0
	 c(1)=1
	 b(1)=1
	 d(1)=0	
	endif
	   
	if (h1 /= 0.and.ls1 == 0) then
	 c(1)=1
	 b(1)=0
	 d(1)=WL_MAX(por(1),rosdry(1),wi1(1))
	 dhwfsoil1 = dt*(wl1(2)-wl1(1))*rosdry(1)*(1-por(1))/row
     & /dzs(1)*(lammoist(1)+lammoist(2))/2-
     &  (WL_MAX(por(1),rosdry(1),wi1(1))-wl1(1))
     & *rosdry(1)*(1-por(1))*dzss(1)/row
      endif	 
	  
	lammoist1 = (lammoist(ns-1)+lammoist(ns))/2
       c(ns)=-lammoist1/dzs(ns-1)
	 a(ns)=-lammoist1/dzs(ns-1)
	 d(ns)=Wflow
           
      do i=2,ns-1
       lammoist2 = (lammoist(i)+lammoist(i+1))/2
       lammoist1 = (lammoist(i-1)+lammoist(i))/2
       dzmean = (dzs(i-1)+dzs(i))/2
       a(i)=-dt/dzmean*lammoist1/dzs(i-1)
	 b(i)=-dt/dzmean*lammoist2/dzs(i)
	 c(i)=-dt/dzmean*(lammoist1/dzs(i-1)+
     & lammoist2/dzs(i))-1
	 d(i)=-wl1(i)
       !wl2(i) = wl1(i)+dt*(lammoist2/dzs(i)*(wl1(i+1)-wl1(i))-
      !& lammoist1/dzs(i-1)*(wl1(i)-wl1(i-1)))/dzmean
	enddo   
		
	call progonka (a,b,c,d,wl2,ns) 	
		
      
!STEP 2 OF SPLITTING-UP METHOD FOR MOISTURE EQUATION: 
!EVOLUTION OF MOISTURE DUE TO GRAVITATIONAL INFILTRATION 	     	    

      do i=1, ns-1
	 wlmax_i=WL_MAX(por(i+1),rosdry(i+1),wi1(i+1))
	 if (wl2(i+1)/wlmax_i>0.98.or.wlmax_i<0.01) then
	  filtr(i) = 0
       else
	  wlmax_i = (POR(i)+POR(i+1))*ROW/((rosdry(i)+rosdry(i+1))*
     &  (1-(por(i)+por(i+1))/2)) - (wi1(i)+wi1(i+1))/2*row/roi      
	  filtr(i) = (flwmax(i+1)+flwmax(i))/2*((wl2(i+1)+
     &  wl2(i))/2/wlmax_i)**(bh(i+1)+bh(i)+3)
	 endif
      enddo
	 
      do i=2,ns-1
       wl3(i) = wl2(i) + dt*(filtr(i-1)-filtr(i))/dzss(i)
      enddo
     	 
      if ((h1==0.and.l1/=0).or.ls1 /= 0) then
       wl3(1) = wl2(1) + dt*(-filtr(1)+0)/dzss(1)
      endif

	if (h1/=0.and.ls1 == 0) then
	 wl3(1) = WL_MAX(por(1),rosdry(1),wi1(1))
	 dhwfsoil2 = - filtr(1)*rosdry(1)*(1-por(1))/row*dt 
      endif
	
	filtr_low = 0
	wl3(ns) = wl2(ns) + dt*(-filtr_low+filtr(ns-1))/dzss(ns)
	
     	do i=1,ns
	 if (wl3(i)>WL_MAX(por(i),rosdry(i),wi1(i))) then
	  surplus=wl3(i)-WL_MAX(por(i),rosdry(i),wi1(i))
	c1:do j=1,ns
	    if (WL_MAX(por(j),rosdry(j),wi1(j))-wl3(j)>0) then
	     if (WL_MAX(por(j),rosdry(j),wi1(j))
     &	  -wl3(j)>surplus*rosdry(i)*(1-por(i))*dzss(i)/
     &      (rosdry(j)*(1-por(j))*dzss(j))) then
       	  wl3(j)=wl3(j)+surplus*rosdry(i)*(1-por(i))*dzss(i)/
     &	  (rosdry(j)*(1-por(j))*dzss(j))
		  surplus=0
            exit c1
	     else
	      surplus=surplus-(WL_MAX(por(j),rosdry(j),wi1(j))-wl3(j))*
     &	  (rosdry(j)*(1-por(j))*dzss(j))/
     &      (rosdry(i)*(1-por(i))*dzss(i))
	   	  wl3(j)=WL_MAX(por(j),rosdry(j),wi1(j)) 
           endif
          endif
	   enddo c1
	  dhwfsoil3 = dhwfsoil3 + surplus*rosdry(i)*(1-por(i))*dzss(i)/row
	  wl3(i) = wl_max(por(i),rosdry(i),wi1(i))
	 endif
	 if(wl3(i)<0) then
	  if (i/=1) then 
	   do j=i-1,1,-1
	    if (wl3(j)>-wl3(i)*rosdry(i)*(1-por(i))*dzss(i)/
     &  	(rosdry(j)*(1-por(j))*dzss(j))) then
	     wl3(j)=wl3(j)+wl3(i)*rosdry(i)*(1-por(i))*dzss(i)/
     &     (rosdry(j)*(1-por(j))*dzss(j))
	    goto 1
	    endif
         enddo
        endif
	  if (i/=ns) then 
	   do j=i+1,ns
	    if (wl3(j)>-wl3(i)*rosdry(i)*(1-por(i))*dzss(i)/
     &    (rosdry(j)*(1-por(j))*dzss(j))) then
	     wl3(j)=wl3(j)+wl3(i)*rosdry(i)*(1-por(i))*dzss(i)/
     &     (rosdry(j)*(1-por(j))*dzss(j))
	    goto 1
	    endif
         enddo
        endif
       !dhwfsoil = dhwfsoil - abs(wl4(i))*rosdry(i)*(1-por(i))*dzss(i)
     & ! /row
1	  wl3(i)=0.
	 endif
      enddo
      
! STEP 3 OF SPLITTING-UP METHOD : PHASE PROCESSES

20    c2:do i=1, ns
       if (wi1(i)>0.and.Tsoil2(i)>0.1) then
        Potphenergy = (Tsoil2(i)-0.1)*csoil(i)*rosoil(i)*dzss(i)
        if (Potphenergy>=wi1(i)*rosdry(i)*dzss(i)*(1-por(i))*Lwi) then
         wl4(i) = wl3(i) + wi1(i)
	   wi2(i) = 0
	   Tsoil3(i) = 0.1+(Potphenergy-(wi1(i)*rosdry(i)*dzss(i)*(
     &   1-por(i))*Lwi))/(csoil(i)*rosoil(i)*dzss(i))
        else
	   wl4(i) = wl3(i) + Potphenergy/(rosdry(i)*dzss(i)*(1-por(i))*Lwi)
	   wi2(i) = wi1(i) - Potphenergy/(rosdry(i)*dzss(i)*(1-por(i))*Lwi)
	   Tsoil3(i) = 0.1
	  endif
       else
        if (wl3(i)>=0.and.Tsoil2(i)<-0.1) then
         Potphenergy = -(Tsoil2(i)+0.1)*csoil(i)*rosoil(i)*dzss(i)
         watavailtofreeze = (wl3(i)-unfrwat(Tsoil2(i),i))*rosdry(i)*
     &   dzss(i)*(1-por(i))*Lwi
	   if (watavailtofreeze<0) then 
	    !cc = csoil(i)*rosoil(i)*dzss(i)
	    !cc1 = rosdry(i)*dzss(i)*(1-por(i))*Lwi
	    !bb = (-Potphenergy-0.1*cc-wl3(i)*cc1)/cc
	    !aa = cc1/cc
	    !call phase_iter(aa,bb,i,Tsoil3(i))
	    Tsoil3(i) = Tsoil2(i) +  
     &   (wl3(i)-unfrwat(Tsoil2(i),i))*rosdry(i)*dzss(i)*(1-por(i))*Lwi/
     &    (csoil(i)*rosoil(i)*dzss(i))
	    wi2(i) = wi1(i) + (wl3(i)-unfrwat(Tsoil2(i),i))
          wl4(i) = unfrwat(Tsoil2(i),i)
	    cycle c2  
         endif
	   if (Potphenergy>=watavailtofreeze) then
	    Tsoil3(i) = -0.1 - (Potphenergy-watavailtofreeze)/
     &    (csoil(i)*rosoil(i)*dzss(i))
	    wi2(i) = wi1(i) + (wl3(i)-unfrwat(Tsoil2(i),i))
	    wl4(i) = unfrwat(Tsoil2(i),i)
	   else
	    wl4(i) = wl3(i)-Potphenergy/(rosdry(i)*dzss(i)*(1-por(i))*Lwi)
	    wi2(i) = wi1(i)+Potphenergy/(rosdry(i)*dzss(i)*(1-por(i))*Lwi)
	    Tsoil3(i) = - 0.1
	   endif
	  else
	   wl4(i) = wl3(i)
	   wi2(i) = wi1(i)
	   Tsoil3(i) = Tsoil2(i)
	  endif
       endif
      enddo c2
      
      max1=0.
      do i=1,ns
       if (Tsoil3(i)<-0.1.and.dabs(wl4(i)-unfrwat(Tsoil3(i),i))>max1)
     & then
	   max1=dabs(wl4(i)-unfrwat(Tsoil3(i),i))
	 endif
	enddo    
	  
      if (max1>0.01.and.iter3<20) then
	 !wl3=(wl4+wl3)/2.
	 !wi1=(wi2+wi1)/2.
	 !Tsoil2=(Tsoil3+Tsoil2)/2.
	 wl3=wl3+alp(iter+1)*(wl4-wl3)
	 wi1=wi1+alp(iter+1)*(wi2-wi1)
	 Tsoil2=Tsoil2+alp(iter+1)*(Tsoil3-Tsoil2)
	 iter = iter + 1
	 iter3 = iter3 + 1
	 if(iter==8) iter = 0
	 goto 20
      else
	 iter = 0
	 iter3 = 0 
	endif   
	
	!evol=0

      do i=1,ns
       if (wi2(i)<0) then
	  wl4(i)=wl4(i)+wi2(i)
	  wi2(i)=0
       endif
	enddo

	do i=1,ns
	 if (wl4(i)>POR(i)*ROW/(rosdry(i)*(1-por(i)))-wi2(i)*row/roi) then
	  surplus = wl4(i)-(POR(i)*ROW/(rosdry(i)*(1-por(i))) - wi2(i)*
     &  row/roi)
	c3:do j=1,ns
	    if (POR(j)*ROW/(rosdry(j)*(1-por(j))) - wi2(j)*row/roi
     &	-wl4(j)>0) then
	     if (POR(j)*ROW/(rosdry(j)*(1-por(j))) - wi2(j)*row/roi
     &	  -wl4(j)>surplus*rosdry(i)*(1-por(i))*dzss(i)/
     &      (rosdry(j)*(1-por(j))*dzss(j))) then
       	  wl4(j)=wl4(j)+surplus*rosdry(i)*(1-por(i))*dzss(i)/
     &	  (rosdry(j)*(1-por(j))*dzss(j))
		  surplus=0
            exit c3
	     else
	      surplus=surplus-(POR(j)*ROW/(rosdry(j)*(1-por(j))) -
     & 	  wi2(j)*row/roi-wl4(j))*(rosdry(j)*(1-por(j))*dzss(j))/
     &      (rosdry(i)*(1-por(i))*dzss(i))
		  wl4(j)=POR(j)*ROW/(rosdry(j)*(1-por(j))) - wi2(j)*row/roi 
           endif
          endif
	   enddo c3
	  dhwfsoil4 = dhwfsoil4 + surplus*rosdry(i)*(1-por(i))*dzss(i)/row
	  wl4(i) = POR(i)*ROW/(rosdry(i)*(1-por(i))) - wi2(i)*row/roi 
	 endif
	 if(wl4(i)<0) then
	  if (i/=1) then 
	   do j=i-1,1,-1
	    if (wl4(j)>-wl4(i)*rosdry(i)*(1-por(i))*dzss(i)/
     &	(rosdry(j)*(1-por(j))*dzss(j))) then
	     wl4(j)=wl4(j)+wl4(i)*rosdry(i)*(1-por(i))*dzss(i)/
     &     (rosdry(j)*(1-por(j))*dzss(j))
	    goto 2
	    endif
         enddo
        endif
	  if (i/=ns) then 
	   do j=i+1,ns
	    if (wl4(j)>-wl4(i)*rosdry(i)*(1-por(i))*dzss(i)/
     &	(rosdry(j)*(1-por(j))*dzss(j))) then
	     wl4(j)=wl4(j)+wl4(i)*rosdry(i)*(1-por(i))*dzss(i)/
     &     (rosdry(j)*(1-por(j))*dzss(j))
	    goto 2
	    endif
         enddo
        endif
       !dhwfsoil = dhwfsoil - abs(wl4(i))*rosdry(i)*(1-por(i))*dzss(i)
     & ! /row
2	  wl4(i)=0
	 endif
      enddo
            
  !!!!! II. FREEZING AND MELTING IN CASE OF TSOIL<MELTING POINT (O C)!!!!!

		!do i=1, ns
		! if (Tsoil3(i)<-1.) then
		!   wi3(i) = wi2(i) + (wl4(i)-unfrwat(Tsoil3(i),i))
        !   wl5(i) = unfrwat(Tsoil3(i),i)
		!   Tsoil4(i) = Tsoil3(i) + & 
		!   (wl4(i)-unfrwat(Tsoil3(i),i))*rosdry(i)*dzss(i)*(1-por(i))*Lwi/&
		!   (csoil(i)*rosoil(i)*dzss(i))
		! else
		!   wi3(i) = wi2(i)
		!   wl5(i) = wl4(i)
		!   Tsoil4(i) = Tsoil3(i)
        ! end if
		! if (wl5(i)>POR(i)*ROW/(rosdry(i)*(1-por(i))) - wi3(i)*row/roi) then
	    !  surplus = wl5(i)-(POR(i)*ROW/(rosdry(i)*(1-por(i))) - wi3(i)*row/roi)
	    !  dhwfsoil = dhwfsoil + surplus*rosdry(i)*(1-por(i))*dzss(i)/row
		!  wl5(i) = POR(i)*ROW/(rosdry(i)*(1-por(i))) - wi3(i)*row/roi 
        ! endif
	    ! if(wl5(i)<0) then 
		!  dhwfsoil = dhwfsoil - abs(wl5(i))*rosdry(i)*(1-por(i))*dzss(i)/row
		!  wl5(i)=0
	    ! endif
        !enddo
   
	   ! evol=0
	!	evol1=0

	!	do i=1,ns
	!	 evol=evol+(Tsoil3(i)-Tsoil1(i))*rosoil(i)*csoil(i)*dzss(i)
	!	 evol1=evol1+(wi2(i)-wi1(i))*rosdry(i)*(1-por(i))*dzss(i)*Lwi
      !   enddo

      dhwfsoil=dhwfsoil1+dhwfsoil2+dhwfsoil3+dhwfsoil4

	!if (time/month>133.133.and.time/month<140.) then
	! write (133,'(f8.3,2e12.3)') time/month,dhwfsoil1,dhwfsoil2
	! do i=1,ns
	!  write  (133,'(5f7.2)') wl1(i),wl2(i),wl3(i),wl4(i),
     &!  WL_MAX(por(i),rosdry(i),wi1(i))
      ! enddo
      !endif

     	  
  	Tsoil1 = Tsoil3
	wl1 = wl4
	wi1 = wi2

      if (firstcall) firstcall=.false.
      RETURN
	END SUBROUTINE soilforlake 


      SUBROUTINE snowtemp(snowmass,snowmass_init,a,b,c,d,Temp,dt)
      
!     SNOWTEMP calculates temperature profile in the snow cover      
      
      use numeric_params    
	use atmos
	use phys_constants2
	use driving_params
	use arrays
	implicit none


!-----------------------------MAIN VARIABLES-------------------------------
!         arrays: Tsn - snow temperature, C
!                 lams - thermal conductivity of snow
!                
!            rofresh - density of fresh snow, kg/m**3           
!            snowmass - snow mass, kg/m**2

!            Boundary conditions:
!            t2 - air temperature (2 m), C        
       
      real(8)   xx, yy, dzsn(ms),alp(10),psi(10)
	real(8)   rofresh, t2, snowmass, densav, pheat, pmelt, snowmass_init,
     & csurf,Ts_old,Tf_old1, Tf_old2, Tf_old12,mi,ma,month,pi,
     & dhwfsoil,CCT(ML)
      real(8)  cv
	real(8) q(110)
	real(8), dimension(1:ms)::lams,cs,Tsn
      real(8) dt, h1, l1, hs1, ddz, ls1
      real(8) dz
      real(8) ST,PGR,TGROLD,QGROLD,RADIAT,PRECIP,WSOLD,SNOLD,
     & ZS,thsoil,whsoil,BOLD,RF1,RF2,SNMELT,HSNOW,
     & HS,ES,TGRNEW,QGRNEW,WSNEW,SNNEW,RUNOFF,
     & ElatOld,HFold,PRSold,shortwave,extinct
      real(8) AL,DLT,DVT,ALLL,DL,ALV,DV,Z,T,WL,WV,WI,dens,Tf(2)
	real(8), dimension(1:350)::a,b,c,d,Temp,a1,b1,c1,d1,Temp1

	integer(4) itop, itop_old, num, itop1,itop2,iter,iter2,i,
     & iyear,imonth,iday  

	logical flux_calc,ts_iter_more,soil_ts_ext_iter

      common /watericesnowarr/ lams, q
	common /layers/ h1,l1,hs1,ls1
	common /BL/ ST,PGR,TGROLD,QGROLD,RADIAT,WSOLD,SNOLD,
     & ZS,thsoil,whsoil,BOLD,RF1,RF2,SNMELT,HSNOW,
     & HS,ES,TGRNEW,QGRNEW,WSNEW,SNNEW,RUNOFF,
     & ElatOld,HFold,PRSold,extinct(ms)
	common /atmos2/ shortwave, precip
	common /SOILSOL/ AL(ML),DLT(ML),DVT(ML),ALLL(ML),DL(ML),
     & ALV(ML),DV(ML),Z(ML),T(ML),WL(ML),WV(ML),WI(ML),dens(ms)
	common /SOILDAT/ dz(ms),itop
	common/ts_ext_iter/ Tf_old1,soil_ts_ext_iter 
	common /snow_char/ Tsn,cs
	
	SAVE
	 
!     BODY	

      pi=3.141592654
	month=30*24*60*60
	ddz = 1./float(M)
    
	if (flag_snow_init == 1) then
	 itop = ms-2 
	 hs1=0.02
	 klmn=0	
       rofresh = 67.9 + 51.3*dexp(tempair/2.6)
	! snowmass = rofresh*0.012
	 do i = itop, ms-1
	  WL(i) = 0.0
	  dens(i) = rofresh 
	  dz(i) = 0.01
	 end do
	 wl(ms)=0.
	 dens(ms)=rofresh
	 T(ms) = Ti1(1)
	 do i = itop, ms-1
	  T(i) = tempair
       end do
	 do i=1,itop-1
	  T(i)=0
	  wl(i)=0
	  dens(i)=0
	  dz(i)=0
       enddo
       !T(ms)=0
	 !wl(ms)=0
	 !dens(ms)=0
	 !dz(ms)=0
	 totalmelts=0.
	 totalprecips=0.
	 totalevaps=0.
      end if

	
	!ma=8.5
	!mi=1.5
	!data psi/6,3,7,2,5,4,8,1,0,0/   !/3,2,4,1,0,0,0,0,0,0/
	!do i=1,8 
	! alp(i)=2*(mi+ma-(ma-mi)*cos((2*psi(i)-1)/8*pi))**(-1)
      !enddo
	
       do i = itop, ms-2
        dzsn(i) = (dz(i)+dz(i+1))/2
       enddo
       do i = itop, ms-1
        !densav = (dens(i) + dens(i+1))/2
        !lams(i) =  0.419*(6.*(densav/row)**4+1.9*(densav/row)+0.05)
	  lams(i) =  0.419*(6.*(dens(i)/row)**4+1.9*(dens(i)/row)+0.05)
        xx = wl(i)*row/(dens(i)*dz(i)) 
	  cs(i) = ci*(1-xx)+cw*xx
	 end do

      call soilforlake(dt, h1, l1, ls1, hs1, 
     & dhwfsoil,a,b,c,d,Temp,num)
	do i=itop,ms
	 T(i)=Tsn(i)	
      enddo
	t2=tempair
	Erad=shortwave*(1-albedoofsnow)
	Elatent=xlew
	call addPrecipToSnow(t2, T(itop), snowmass, iyear, imonth, iday, 
     & CCT, pmelt, pheat,dt)
	call snow_calc(t2,Erad,T(itop),snowmass,iyear, 
     &imonth, iday, CCt, pmelt, pheat,dt)
	
	if (flag_snow_init == 1) then
	! snowmass_init = snowmass  - (precip*dt*row-snmelt*dt*row-
      !& Elatent/Liv*dt)
       flag_snow_init = 0
	endif 
      
	totalevaps = totalevaps+Elatent/(row*Lwv)*dt
	totalprecips = totalprecips+precip*dt
	totalmelts = totalmelts + snmelt*dt
	      
      hs1=hsnow
	if (hs1==0) hs1=0.00001
            
      END SUBROUTINE snowtemp

	SUBROUTINE snow_calc(t2,Erad_sc,ts,snow_sc,iyear,imonth,
     & iday,CCT,pmelt,pheat,dt)
      use phys_constants
      use phys_parameters
	use numeric_params
	use atmos
	use driving_params
	use arrays
	implicit none

CCCCCCCCCCCCCCC input parameters:
CC      e, cal/(cm2*s), --- surface heat balance: e = Radiation-SIGMA*(TS**4)-Latent-Sensible-InSoil
CC      Elatent, cal/(cm2*s), --- latent heat flux.
CC      ts, deg C, --- surface temperature.
CC      precip, cm --- precipitation during dt.
CC
CCCCCCCCC  some other parameters: CCCCCCCCCCC
CC      densny, g/cm3       --- density of fresh snow, depends on temperature at 2m
CC      dznorm, cm          --- initial thickness of snow layers
CC      dzmin, cm           --- minimal thickness of snow layers
CC      dt, sec             --- timestep
CC      ms                  --- max. number of layers in snowpack
CC      cwat, cal/(gr K)    --- the water specific heat content
CC      csnow, cal/(gr K)   --- the snow specific heat content
CC      rhow, gr/cm3        --- the liquid water density
CC      PLI,  cal/g         --- latent heat for freezing/melting
CC      PL,  cal/g          --- latent heat for evapor./condens.
CC      whc                 --- hydraulic constant
CC      hcond, cm/sec       --- hydraulic conductivity
CC
CCCCCCCCCCC  values to be calculate after every timestep:
CC      wl(i), cm       --- moisture content in snow layer No.i
CC      t(i), deg.C     --- mean temperature of snow layer No.i
CC      dz(i), cm       --- hickness of snow layer No.i
CC      dens(i), g/cm3  --- mean density of snow layer No.i
CC      hsnow, cm       --- snow depth
CC      swe, cm         --- water equivalent snow depth
CC      snmelt, cm/sec  --- mean snow melting speed during dt
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	
	INTEGER(4) j,i,k,itop,iyear,imonth,iday,count_perc,count_melt
     &,nspin
	 
	integer(4) trtr, nmonth,nyear,nday,nhour
	
	common /time/ nyear,nmonth,nday,nhour
	common /layers/ h1,l1,hs1,ls1
      COMMON /SOILDAT/ dz(ms),itop
	COMMON /BL/       ST,PGR,TGROLD,QGROLD,RADIAT,WSOLD,SNOLD,
     &                  ZS,thsoil,whsoil,BOLD,RF1,RF2,SNMELT,HSNOW,
     &                  HS,ES,TGRNEW,QGRNEW,WSNEW,SNNEW,RUNOFF,
     &                  ElatOld,HFold,PRSold,extinct(ms)
	common /atmos2/ shortwave, precip
	COMMON /SOILSOL/ AL(ML),DLT(ML),DVT(ML),ALLL(ML),DL(ML),
     &                  ALV(ML),DV(ML),Z(ML),T(ML),WL(ML),WV(ML),WI(ML),
     &                  dens(ms)
	COMMON /spin/ Cond,Phase

      common /watericesnowarr/ lams, q
      real(8) dt,L,dhs,ddz,xxxx
	real(8) dz,q(110)
	real(8) cs(ms),lams(ms)
	real(8) ST,PGR,TGROLD,QGROLD,RADIAT,PRECIP,WSOLD,SNOLD, 
     &ZS,thsoil,whsoil,BOLD,RF1,RF2,SNMELT,HSNOW,HS,ES,TGRNEW,QGRNEW,
     &WSNEW,SNNEW,RUNOFF,ElatOld,HFold,PRSold
	real(8) AL,DLT,DVT,ALLL,DL,ALV,DV,Z,T,WL,WV,WI,dens
	real(8) extinct,e,erad_sc,ts,snow_sc,dz0,z0z,
     &evap,dzitop,erain,countrad,qbase,dmass,densev,counter,qin,rf,
     &smelti,ein,dzi,ziz,porosityOfSnow,w0,fukt,dztop,smelt,eout,p1,q0,
     &dzmasi,swe,densny,shortwave,erad1,eflux1,Phase,Ei,sigma,h1,l1,hs1,
     &ls1,Cond
	real(8) EInRad,altop,ColdContent,ColdContent1,Eold,T2,
     &rho_dry,precip_sum,bmelt_sum,CCT(ML),delta,a2,eta,P,
     &dens_old,ti(ml)
      real(8) totalevap, pheat,pmelt

	SAVE
       
	precip = precip*dt
	!hsnow = 0.
	!do i = itop, ms-1
      ! hsnow = hsnow + dz(i)
	!end do
	!eFlux = eFlux - Erad
	eflux1 = eflux !/2 
	erad1 =  erad_sc !/2     
	!eflux1 = eflux1*(1.-exp(-0.02*snow))  ! not -0.02*snow*10.                  
	do i=itop,ms-1
	!	wl(i) = wl(i)*dz(i)*dens(i)/rhow    
		ti(i) = t(i)
	end do
	
	snmelt = 0.
	swe = 0.0
	Phase = 0.
c	if(t2 .le. -10) then
c		densny = 0.05
c	elseif(t2 .le. -2) then
c		densny = 5./800.*(t2+18.)
c		densny = 0.1
c	else
c		densny = (t2+4.)/20.
c	end if
	densny = (67.9 + 51.3 * dexp(t2/2.6))
	if(ts .gt. 0.00001) then
		Elatent = 0
	elseif(snow_sc .lt. 0.000001 .and. Elatent .gt. 0.0) then    
		Elatent = 0
	end if

	dz(ms-1) = dz(ms-1) + dz(ms)
	evap=Elatent/(PLI+PL) 
	!totalevap=totalevap+evap/rhow*dt                              
	hsnow=0

	if(itop.eq.ms) then

	else
CC   		--- first top snow layer ---
		i=itop
		dztop=dz(i)
		Erain = pheat                           
		if(extinct(itop).eq.0.0) then
			Eflux1 = Eflux1 + Erad1 + Erain                            
		else	
			if(itop.eq.ms-1) then
				Eflux1 = Eflux1 + Erain + Erad1            
			else
			   Eflux1 = Eflux1 + Erain + Erad1*(1 - extinct(i)**dz(i))      
			end if
			countrad = Erad1*(1 - extinct(i)**dz(i))
		end if               
CC rain 	heat from rain is added as pheat
CC rain 	the rain water itself is assumed to freeze giving heat to top layer
CC rain 	then this heat will melt snow
		if(Eflux1.gt.0.0)then
			if(itop.eq.ms-1) then
        		    if(Ts.ge.-0.0001)then
				  if(Eflux1/PLI*dt .le. DENS(i)*DZ(i)-WL(i)*rhow) then
						smelt = Eflux1/PLI/rhow                          
						Eout=0.
					else
						smelt = (DZ(i)*DENS(i)/rhow - WL(i))/dt
						Eout = Eflux1 - smelt*rhow*PLI    
					end if
				else
					smelt = 0.0
					Eout = 0.
				end if
			else
				smelt = 0.0
				Eout = Eflux1
			end if
		else
			smelt = 0.
			Eout = 0.
		end if

	
Cc
Cc	PERCOLATION
Cc
		smelt = smelt + wl(i)/dt
		wl(i) = 0.
		qbase = smelt

CC
CC	TOP  LAYER THICKNESS AND DENSITY
CC
		dmass=dens(i)*dztop-(evap+smelt*rhow)*dt     
		if(evap.gt.0.0)then
			densev=dens(i)
		else
			densev=rhoi
		endif
		dztop=dztop-(smelt*rhow/dens(i)+evap/densev)*dt   
	  if(dztop.eq.0.) dztop = 0.000001	
		dens(i)= dmax1(dmass/dztop,0.d0)
		if(dens(i).gt.rhoi) then
	      	dztop = dztop*dens(i)/rhoi
	         	dens(i)=rhoi
			dmass=dens(i)*dztop
		endif
		dz(i)=dztop
		if(itop .eq.ms-1) snmelt=qbase
	end if 


	counter = dz(itop)
	if(itop.lt.ms-1) then
CC
CC INTERNAL LAYERS
CC
		do i=itop+1,ms-1
			qin = qbase                                     
			rf = 0.0                                        
     			smelti = 0.0                                    
			Ein=Eout                                        
			T(i) = T(i)*dz(i)*dens(i)/(qbase*rhow*dt+dz(i)*dens(i))  
			dzi=dz(i)                                       
			wl(i) = wl(i) + qbase*dt
			if(wl(i).lt.0.) wl(i) = 0.
			w0=wl(i)
			if(extinct(i).eq.0.0) then
				Ei = Ein                
			else
				if(i.eq.ms-1) then
					Ei=Ein + (Erad1 - countrad)
				else
					Ei=Ein + Erad1*(extinct(i)**(cou
     &					nter)-extinct(i)**(counter+dz(i)) )
					EInRad = Erad1*(extinct(i)**(counter) -
     &					extinct(i)**(counter+dz(i)) )
				end if                            
				countrad = countrad + Erad1*(extinct(i)**(cou
     &              nter)-extinct(i)**(counter+dz(i)) )
				counter = counter+dz(i)
			end if

			dz(i)=dz(i)+qbase*dt  
			ziz=dzi/dz(i)
			dens(i)=ziz*dens(i)+(1.0-ziz)*rhow

            
			if(Ei.gt.0.0) then
				if(T(i).ge.0.0001) then  ! 0.0001
					Ei = Ei + csnow*dens(i)*T(i)*dz(i)/dt
					T(i) = 0.0
					if(Ei/PLI*dt .le. (DZ(i)*DENS(i)-WL(i)*rhow)) then
						smelti = Ei/PLI/rhow 
						Eout = 0.
					else
						smelti = (DZ(i)*DENS(i)/rhow - WL(i))/dt
						Eout = Ei - smelti*rhow*PLI 
						dzi = 0.
					end if
					wl(i) = wl(i) + smelti*dt

				else
Cc		T<0
					T(i) = T(i) + Ei*dt/(csnow*dens(i)*dz(i))
					if(T(i).lt.-0.0001) then   !0.0001
					if(wl(i) .gt. -T(i)*csnow*dz(i)*dens(i)/PLI/rhow) 
     &							then
							rf = -T(i)*csnow*dz(i)*dens(i)/PLI/dt/rhow 
							T(i) = 0.0
							wl(i) = wl(i) - rf*dt
						else
							rf = wl(i)/dt
			       T(i) = T(i) + rf*rhow*dt*PLI/(csnow*dz(i)*dens(i))    
							wl(i) = 0.
						end if
						Eout = 0.
					end if
					if(T(i).ge.0.0001)then  !0.0001
						if(csnow*dens(i)*T(i)*dz(i)/PLI .le.
     &						dz(i)*dens(i)) then
						smelti = csnow*dens(i)*T(i)*dz(i)/dt/PLI/rhow 
							Eout = 0.
						else
							smelti = dz(i)*dens(i)/dt/rhow  
							Eout = csnow*dens(i)*T(i)*dz(i)/dt
     &                              - smelti*PLI*rhow 
							dzi = 0.
						end if
						T(i) = 0.0
						wl(i) = wl(i) + smelti*dt
					endif
				endif
              else

Cc	E = 0 (E<0 is impossible case)
				if(T(i).ge.0.0001) then     !0.0001
					if(csnow*dens(i)*T(i)*dz(i)/PLI .le.
     &                        dz(i)*dens(i)) then
						smelti = csnow*dens(i)*T(i)*dz(i)/dt/PLI/rhow  
						T(i) = 0.0
						Eout = 0.
						wl(i) = wl(i) + smelti*dt
					else
						smelti = dz(i)*dens(i)/dt/rhow  
						Eout = csnow*dens(i)*T(i)*dz(i)/dt
     &                            - smelti*PLI*rhow  
						T(i) = 0.0
						dzi = 0.
						wl(i) = wl(i) + smelti*dt
					end if
				else
Cc	T<0
				if(wl(i) .gt. -csnow*dens(i)*T(i)*dz(i)/PLI/rhow) then  
						rf = -csnow*dens(i)*T(i)*dz(i)/dt/PLI/rhow   
						T(i) = 0.0
						wl(i) = wl(i) - rf*dt
					else
						rf = wl(i)/dt
						wl(i) = 0.
					T(i) = T(i) + rf*rhow*dt*PLI/(csnow*dz(i)*dens(i))
					end if
					Eout = 0.
					if(T(i).ge.0.0001) then   !0.0001
						if(csnow*dens(i)*T(i)*dz(i)/PLI .le.
     &                         dz(i)*dens(i)) then
						smelti = csnow*dens(i)*T(i)*dz(i)/dt/PLI/rhow 
							Eout = 0.0
						else
							smelti = dz(i)*dens(i)/dt/rhow
							Eout = csnow*dens(i)*T(i)*dz(i)/dt
     &                            - smelti*PLI*rhow 
							dzi = 0.
						end if
						T(i) = 0.0
						wl(i) = wl(i) + smelti*dt
					endif
				endif
			endif
	
Cc
Cc
Cc	internal q
Cc
			rho_dry = (dens(i)-w0/(dz(i))*RHOW)/
     *			(1-w0/(dz(i)))
			if(rho_dry.lt.0.) WRITE(*,*) 'rho_dry = ', rho_dry

!     PorosityOfSnow = 1 - rho_dry/rhoi - wl(i)/dz(i)*(1-rhow/rhoi) wl(i)/rhow/dz(i)
	       PorosityOfSnow = 1 - rho_dry/rhoi - wl(i)/dz(i)
			PorosityOfSnow = dmax1(PorosityOfSnow,whc + 0.1)
			p1 = PorosityOfSnow - whc   ! whc = water holding capacity
			fukt = wl(i)/dz(i)
			if(fukt.gt.whc)then
				fukt = (fukt-whc)/p1
				q0 = hcond*fukt**3.0
	
				qbase = dmin1(q0,wl(i)/dt)
				wl(i) = wl(i) - qbase*dt
			else
				qbase = 0.0
				q0 = 0.0
			end if
			if(dzi .le. 0.000001) then   
				qbase = qbase + wl(i)/dt
				dz(i)=dzi
			else
				dzmasi=dz(i)*dens(i)-qbase*dt*rhow   
	            dzi=dz(i)+rhow*(-qbase/RHOW-rf/RHOW+rf/rhoi-  
     &                    smelti/rho_dry+smelti/RHOW)*dt     
				dens(i)=dmax1(dzmasi/dzi,0.d0)
				dz(i)=dzi
				if(dens(i).gt.rhoi)then
					dz(i)=dz(i)*dens(i)/rhoi
					dens(i)=rhoi
				end if
			end if
			Phase = Phase                                 !diagnostic variable
     &    			+ (t(i)-ti(i))*dens(i)*dz(i)/snow_sc/dt
	
			if(i.eq.ms-1) snmelt=qbase
		end do
      end if
	if (itop==ms-1.and.dz(itop)<0.000001) then
	 dz(itop)=0
	 dens(itop)=0
	 goto 123
	endif
      if(itop .lt. ms .and. dz(itop) .le. 0.000001) then   
		dz(itop)=0.
		dens(itop)=0.
		itop=itop+1
      end if
		
222   do i=ms-1,itop,-1
		if(dz(i).le.0.000001) then   
			do k=i-1,itop,-1
				if(dz(k).gt.0.000001) goto 111   
			end do
	
111			do j=i,itop,-1
				dz(j)=dz(k+j-i)
				dens(j)=dens(k+j-i)
				t(j)=t(k+j-i)
				wl(j)=wl(k+j-i)
			end do
			itop=itop-k+i
			goto 222
		end if
      end do
123   continue
	 
333   do i=itop,ms-1,1
		if(dz(i).lt.dzmin.and.itop.lt.ms-1) then
		if(i.lt.ms-1) then
			dens(i+1)=dens(i+1)*dz(i+1)/(dz(i)+dz(i+1))+
     &				dens(i)*dz(i)/(dz(i)+dz(i+1))
			t(i+1)=t(i+1)*dz(i+1)/(dz(i)+dz(i+1))+
     &				t(i)*dz(i)/(dz(i)+dz(i+1))
			wl(i+1)=wl(i+1)+wl(i)
			dz(i+1)=dz(i+1)+dz(i)
			do j=i,itop+1,-1
				dens(j)=dens(j-1)
				t(j)=t(j-1)
				wl(j)=wl(j-1)
				dz(j)=dz(j-1)
			end do
			itop=itop + 1
			goto 333
		end if
		end if
      end do
	
      swe = 0.
      hsnow = 0.
      do j=ms-1,itop,-1
          swe = swe + dz(j)*dens(j)
          hsnow = hsnow + dz(j)
      end do
       
cccccccccccccc snow densification due to gravity and metamorphism cccccccccccccc
	a2 = 0.00000066
	sigma = 75.                             ! metamorphic processes, Pa
	do j = itop+1,ms-1
		if(dens(j).lt.3*densny) then
			dens_old = dens(j)
			eta =                           ! compactive viscosity of snow 
     &			a2*dexp(19.3*dens(j)/rhoi)*
     &			dexp(67300/8.314/(t(j)+273.15))
			P = 0.0                         ! gravity, Pa
			do i = j,itop,-1
				!P = P + dens(i)/1000.*grav*(dz(i))*10000.
	            P = P + dens(i)*grav*dz(i)
			end do
			dens(j) = dens(j) + 0.001*(dt*dens(j)*1000.*(sigma+P)/eta)
                          dens(j) = dmin1(dens(j),rhoi)
			dz(j) = dz(j) * dens_old/dens(j)
		end if
	end do
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      swe = 0.
      hsnow = 0.
      do j=ms-1,itop,-1
          swe = swe + dz(j)*dens(j)
          hsnow = hsnow + dz(j)
      end do
      !do m = itop,ms-1
	!	por(m) = 1 - dens(m)/rhoi
      !end do
	
      snow_sc = swe

        !! dz(ms) !!
	if(itop.ne.ms-1) then
		dz(ms) = dz(ms-1)/3.
		dz(ms-1) = dz(ms-1) - dz(ms)
	else
		dz(ms) = dz(ms-1)/2.
		dz(ms-1) = dz(ms-1) - dz(ms)
	end if
	if(itop.eq.ms-1) then
		z(ms-1) = z(ms) - dz(ms-1) - dz(ms)
	else
		z(ms-1) = z(ms) - dz(ms-1)/2. - dz(ms)
	end if
	DO j = MS-2,itop+1,-1
		Z(j) = z(j+1)-dz(j)/2.-dz(j+1)/2.
      END DO
	if(itop.lt.ms-1) z(itop) = z(itop+1) - dz(itop+1)/2. - dz(itop)

	!do i=itop,ms-1
!		wl(i) = wl(i)/dz(i)/dens(i)*rhow
!	end do
	!extinct(itop)=0
	do i = itop,ms-1
		!extinct(i) = (0.9 - 0.4/0.45*(dens(i) - 0.05))*100.
	 extinct(i) = dexp(-(0.13*dens(i)+3.4))
	 !extinct(i) = 0.
	end do
	
	do i = 1, itop-1
		dz(i) = 0.
      end do

      
      precip = precip/dt
	!eflux = eflux+Erad
	
      return
	END

      SUBROUTINE addPrecipToSnow(T2,ts,snow,iyear,imonth,
     &	iday,CCT,pmelt,pheat,dt)
      use phys_constants    
	use phys_parameters
	use numeric_params
	use driving_params
	use atmos
	
	implicit none
CCCCCCCCCCCCCCC input parameters:
CC      e, cal/(cm2*s), --- surface heat balance: e = Radiation-SIGMA*(TS**4)-Latent-Sensible-InSoil
CC      Elatent, cal/(cm2*s), --- latent heat flux.
CC      ts, deg C, --- surface temperature.
CC      precip, cm --- precipitation during dt.
CC
CCCCCCCCC  some other parameters: CCCCCCCCCCC
CC      densny, g/cm3       --- density of fresh snow, depends on surface temperature
CC      dznorm, cm          --- initial thickness of snow layers
CC      dzmin, cm           --- minimal thickness of snow layers
CC      dt, sec             --- timestep
CC      ms                  --- max. number of layers in snowpack
CC      cwat, cal/(gr K)    --- the water specific heat content
CC      csnow, cal/(gr K)   --- the snow specific heat content
CC      rhow, gr/cm3        --- the liquid water density
CC      PLI,  cal/g          --- latent heat for freezing/melting
CC      PL,  cal/g          --- latent heat for evapor./condens.
CC      whc                 --- hydraulic constant
CC      hcond, cm/sec       --- hydraulic conductivity
CC
CCCCCCCCCCC  values to be calculate after every timestep:
CC      wl(i), cm       --- moisture content in snow layer No.i
CC      t(i), deg.C     --- mean temperature of snow layer No.i
CC      dz(i), cm       --- thickness of snow layer No.i
CC      dens(i), g/cm3  --- mean density of snow layer No.i
CC      hsnow, cm       --- snow depth
CC      swe, cm         --- water equivalent snow depth
CC      snmelt, cm/sec  --- mean snow melting speed during dt
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      COMMON /SOILDAT/ dz(ms),itop
	COMMON /BL/       ST,PGR,TGROLD,QGROLD,RADIAT,WSOLD,SNOLD,
     &                  ZS,thsoil,whsoil,BOLD,RF1,RF2,SNMELT,HSNOW,
     &                  HS,ES,TGRNEW,QGRNEW,WSNEW,SNNEW,RUNOFF,
     &                  ElatOld,HFold,PRSold,extinct(ms)
	common /atmos2/ shortwave, precip
	COMMON /SOILSOL/ AL(ML),DLT(ML),DVT(ML),ALLL(ML),DL(ML),
     &                  ALV(ML),DV(ML),Z(ML),T(ML),WL(ML),WV(ML),WI(ML),
     &                  dens(ms)
	COMMON /spin/ Cond,Phase
	common /layers/ h1,l1,hs1,ls1
     	real(8) dz
	real(8) ST,PGR,TGROLD,QGROLD,RADIAT,PRECIP,WSOLD,SNOLD,
     &ZS,thsoil,whsoil,BOLD,RF1,RF2,SNMELT,HSNOW,HS,ES,TGRNEW,QGRNEW,
     &WSNEW,SNNEW,RUNOFF,ElatOld,HFold,PRSold
	real(8) AL,DLT,DVT,ALLL,DL,ALV,DV,Z,T,WL,WV,WI,dens
	real(8) extinct,e,erad,ts,snow,dz0,pmelt,pheat,
     &z0z,evap,dzitop,erain,countrad,qbase,dmass,densev,counter,qin,rf,
     &smelti,ein,dzi,ziz,porosityOfSnow,w0,fukt,dztop,smelt,eout,p1,q0,
     &dzmasi,swe,densny,shortwave,Cond,Phase,h1,l1,hs1,ls1
	INTEGER(4) j,i,k,mn,itop,iyear,imonth,iday,count_perc,count_melt
     &,nspin
	real(8) EInRad,altop,T2,rho_dry,
     &precip_sum,bmelt_sum,CCT(ML),delta
      real(8) dt 

	SAVE
      
      precip = precip*dt
		
            
	      !if (snow > 0.) then
			! write (*, *) snow
			! STOP
			!end if
	!T2 = T2 - 273.15
c	if(T2 .le. -10) then
c		densny = 0.05
c	elseif(ts .le. -2) then
c		densny = 5./800.*(t2+18.)
cc		densny = 0.1
c	else
c		densny = (t2+4.)/20.
c	end if
	densny = (67.9 + 51.3 * dexp(t2/2.6))
	if(itop.lt.ms) 	dz(ms-1) = dz(ms-1) + dz(ms) 
	
	snow = 0.
	do mn = itop,ms-1
		snow = snow + dz(mn)*dens(mn)
	end do
	
	
	j=itop
555	continue
	if(dz(itop).lt.dznorm+dzmin) goto 444
	dz(j-1) = dz(j) - dznorm
	dz(j) = dznorm
	z(j-1)= z(itop) 
	z(j)= z(j-1) + dz(j-1)*2.
	wl(j-1)=0.0
	dens(j-1)=dens(j)
	t(j-1)=t(j)
	j=j-1
	itop=j
	goto 555
444	continue
      
	if(T2.gt.0.) then
		pmelt = precip/dt                       
		!pmelt = 0.                       
		pheat = T2*precip*cwat*rhow/dt	    
	else	   
	    pmelt=0.0
		pheat=0.0
		if(itop.eq.ms .and. precip.gt.0.0) then
			itop=itop-1
			z(itop) = z(itop+1)
			t(itop)=t(itop+1)
			dens(itop)=densny
			wl(itop)=0
		end if
		if(itop.ne.ms) then
			dz0 = dz(itop)
			z(itop) = z(itop) - precip*rhow/densny   
			dz(itop) = dz(itop) + precip*rhow/densny   
			if(dz(itop) .lt. dznorm+dzmin) then
				z0z=dz0/dz(itop)
				dens(itop)=z0z*dens(itop)+(1.0-z0z)*densny
			end if
	        j=itop
55			continue
			if(dz(itop) .lt. dznorm+dzmin) goto 54
			dz(j-1) = dz(j) - dznorm
			dz(j) = dznorm
			z(j-1)= z(itop) 
			z(j)= z(itop) - dz(itop)*0.5
			z0z=dz0/(dz(itop)+dz(itop-1))
			dens(itop)=z0z*dens(itop)+(1.0-z0z)*densny
			wl(j-1)=0.0
			dz0 = 0.0
			dens(j-1)=densny
			t(j-1)=t(j)
			j=j-1
			itop=j
			goto 55
54			continue
      	end if
	end if
		
	if(pheat.gt.0.0) then
		if(itop.eq.ms-1) then
			dmass = dens(itop)*dz(itop) + precip*rhow    
			dz(itop) = dz(itop) + precip*rhow/rhoi   
			if(dz(itop).eq.0.) dz(itop) = 0.000001   
			dens(itop)= dmax1(dmass/dz(itop),0.d0)
			if(dens(itop).gt.rhoi) then
				dz(itop) = dz(itop)*dens(itop)/rhoi
				dens(itop)=rhoi
			end if
		else
			if(itop.lt.ms-1) then
				mn = itop+1
				dmass = dens(mn)*dz(mn) + precip*rhow  
				dz(mn) = dz(mn) + precip*rhow/rhoi 
				if(dz(mn).eq.0.) dz(mn) = 0.000001    
				dens(mn)= dmax1(dmass/dz(mn),0.d0)
				if(dens(mn).gt.rhoi) then
	                       		dz(mn) = dz(mn)*dens(mn)/rhoi
					dens(mn)=rhoi
				end if
				wl(mn) = wl(mn) + precip
		T(mn) = T(mn)*dz(mn)*dens(mn)/(precip*rhow+dz(mn)*dens(mn)) 
			end if
		end if
	end if

        !!  dz(ms) !!!
	if(itop.ne.ms-1) then
		dz(ms) = dz(ms-1)/3.
		dz(ms-1) = dz(ms-1) - dz(ms)
	else
		dz(ms) = dz(ms-1)/2.
		dz(ms-1) = dz(ms-1) - dz(ms)
	end if
	if(itop.eq.ms-1) then
		z(ms-1) = z(ms) - dz(ms-1) - dz(ms)
	else 
		z(ms-1) = z(ms) - dz(ms-1)/2. - dz(ms)
	end if
	DO Mn = MS-2,itop+1,-1
		Z(Mn) = z(mn+1)-dz(mn)/2.-dz(mn+1)/2.
      END DO
	if(itop.lt.ms-1) z(itop) = z(itop+1) - dz(itop+1)/2. - dz(itop)

	!if (dz(itop).lt.dznorm+dzmin) goto 777 
	
	swe = 0.0
	do mn=ms-1,itop,-1
		swe=swe+dz(mn)*dens(mn)
	end do
	swe = swe + dz(ms)*dens(ms-1)
	snow = swe
	
	hsnow = 0.
	do i = itop, ms
       hsnow = hsnow + dz(i)
	end do
		   
	precip = precip/dt
	
	return
	END


      SUBROUTINE comsoilforlake
      
!     COMSOILFORLAKE specifies parameters of soil according to soil type      
      
      USE numeric_params
	use driving_params
	use arrays
	implicit none


      real(8),dimension(1:Num_Soil):: b,psi_max, 
     & Porosity,gamma_max,lambda_max,W_0,W_m
      real(8), allocatable:: AL(:),DLT(:),DVT(:),DL(:),ALV(:),
     & DV(:),Z(:),T(:),WL(:),WV(:),WI(:),dens(:)
	real(8) pow, BH_soil,PSIMAX_soil,POR_soil,FLWMAX_soil,
     & DLMAX_soil, WLM0_soil, WLM7_soil
	  
      integer(4) i,j,itop,SoilType
	logical firstcall

	data firstcall /.true./

      SAVE
      if (firstcall) then
       allocate (AL(1:ns+2),DLT(1:ns+2),DVT(1:ns+2),DL(1:ns+2),
     & ALV(1:ns+2),DV(1:ns+2),Z(1:ns+2),T(1:ns+2),WL(1:ns+2),
     & WV(1:ns+2),WI(1:ns+2),dens(1:ns+2))
	endif
	
!     I=1       ! SAND
!     I=2       ! LOAMY SAND
!     I=3       ! SANDY LOAM
!     I=4       ! LOAM
!     I=5       ! SILT LOAM
!     I=6       ! SANDY CLAY LOAM
!     I=7       ! CLAY LOAM
!     I=8       ! SILTY CLAY LOAM
!     I=9       ! SANDY CLAY
!     I=10      ! SILTY CLAY
!     I=11      ! CLAY
!----------------------------------------------------------------
!  NN    b     psi_max  Por    gamma_max lambda_max W_0     W_m
!----------------------------------------------------------------
!   1   4.05    3.50    0.395   0.01760   0.29800   0.02    0.01
!   2   4.38    1.78    0.410   0.01560   0.13900   0.05    0.02
!   3   4.90    7.18    0.435   0.00340   0.13700   0.08    0.03
!   4   5.39    14.6    0.451   0.00069   0.06190   0.20    0.08
!   5   5.30    56.6    0.485   0.00072   0.20400   0.18    0.07
!   6   7.12    8.63    0.420   0.00063   0.03070   0.13    0.06
!   7   8.52    36.1    0.476   0.00024   0.03720   0.27    0.12
!   8   7.75    14.6    0.477   0.00017   0.01240   0.24    0.11
!   9   10.4    6.16    0.426   0.00021   0.00641   0.23    0.10
!  10   10.4    17.4    0.492   0.00010   0.00755   0.30    0.15
!  11   11.4    18.6    0.482   0.00013   0.00926   0.40    0.20
!-----------------------------------------------------------------

      z(1) = 0.

	if(UpperLayer > 0.) then
	 pow = dlog(UpperLayer/depth)/dlog(1.d0/float(ns-1))
       do i=2,ns
	  z(i) = (1.*i/(ns-1))**pow*depth
	 end do
	else
	 do i=2,ns
	  z(i) = (1.*i/(ns-1))*depth
	 end do
	end if
	
	do i=1, ns-1
	 dzs(i)=z(i+1)-z(i)
      end do

	do i=1, ns
	 if (i==1) then
        dzss(i) = dzs(i)/2
       elseif (i==ns) then
	  dzss(i) = dzs(i-1)/2
       else
	  dzss(i) = (dzs(i-1)+dzs(i))/2
       endif 
      end do 
	
      SoilType = 7

      b(1) = 4.05
      b(2) = 4.38
      b(3) = 4.90
      b(4) = 5.39
      b(5) = 5.30
      b(6) = 7.12
      b(7) = 8.52
      b(8) = 7.75
      b(9) = 10.4
      b(10)= 10.4
      b(11)= 11.4

      psi_max( 1) = 3.50/100. 
      psi_max( 2) = 1.78/100. 
      psi_max( 3) = 7.18/100. 
      psi_max( 4) = 14.6/100. 
      psi_max( 5) = 56.6/100. 
      psi_max( 6) = 8.63/100. 
      psi_max( 7) = 36.1/100. 
      psi_max( 8) = 14.6/100. 
      psi_max( 9) = 6.16/100. 
      psi_max(10) = 17.4/100. 
      psi_max(11) = 18.6/100. 

      Porosity( 1) = 0.395
      Porosity( 2) = 0.410
      Porosity( 3) = 0.435
      Porosity( 4) = 0.451
      Porosity( 5) = 0.485
      Porosity( 6) = 0.420
      Porosity( 7) = 0.476
      Porosity( 8) = 0.477
      Porosity( 9) = 0.426
      Porosity(10) = 0.492
      Porosity(11) = 0.482

      gamma_max( 1) = 0.01760/100. 
      gamma_max( 2) = 0.01560/100. 
      gamma_max( 3) = 0.00340/100. 
      gamma_max( 4) = 0.00069/100. 
      gamma_max( 5) = 0.00072/100. 
      gamma_max( 6) = 0.00063/100. 
      gamma_max( 7) = 0.00024/100. 
      gamma_max( 8) = 0.00017/100. 
      gamma_max( 9) = 0.00021/100. 
      gamma_max(10) = 0.00010/100. 
      gamma_max(11) = 0.00013/100. 

      lambda_max( 1) = 0.29800/10000. 
      lambda_max( 2) = 0.13900/10000.
      lambda_max( 3) = 0.13700/10000.
      lambda_max( 4) = 0.06190/10000.
      lambda_max( 5) = 0.20400/10000.
      lambda_max( 6) = 0.03070/10000.
      lambda_max( 7) = 0.03720/10000.
      lambda_max( 8) = 0.01240/10000.
      lambda_max( 9) = 0.00641/10000.
      lambda_max(10) = 0.00755/10000.
      lambda_max(11) = 0.00926/10000.

      W_0( 1) = 0.02
      W_0( 2) = 0.05
      W_0( 3) = 0.08
      W_0( 4) = 0.20
      W_0( 5) = 0.18
      W_0( 6) = 0.13
      W_0( 7) = 0.27
      W_0( 8) = 0.24
      W_0( 9) = 0.23
      W_0(10) = 0.30
      W_0(11) = 0.40

      W_m( 1) = 0.01
      W_m( 2) = 0.02
      W_m( 3) = 0.03
      W_m( 4) = 0.08
      W_m( 5) = 0.07
      W_m( 6) = 0.06
      W_m( 7) = 0.12
      W_m( 8) = 0.11
      W_m( 9) = 0.10
      W_m(10) = 0.15
      W_m(11) = 0.20

	if(SoilType>0) then
	 do j = 1, ns
	  BH(j)     = b(SoilType)               ! PARAMETER B, DIMENSOINLESS
	  PSIMAX(j) = - psi_max(SoilType)       ! SAT. WATER POTENTIAL, M.          
    	  POR(j)    = Porosity(SoilType)        ! POROSITY, DIMENSIONLESS
	  FLWMAX(j) = gamma_max(SoilType)       ! SAT. HYDR. CONDUCTIVITY, M/S
	  DLMAX(j)  = lambda_max(SoilType)      ! SAT. WATER DIFFUSIVITY, ?KG/(M*SEC)
	  WLM0(j)   = W_0(SoilType)             ! MAXIMAL UNFREEZING WATER AT 0C
	  WLM7(j)   = W_m(SoilType)             ! MAXIMAL UNFREEZING WATER AT T<<0C
	 end do
	else
!	 do j = 1, ns
!	  BH(j)     = BH_soil           ! PARAMETER B, DIMENSOINLESS
!	  PSIMAX(j) = PSIMAX_soil       ! SAT. WATER POTENTIAL, M.          
!    	  POR(j)    = POR_soil          ! POROSITY, DIMENSIONLESS
!	  FLWMAX(j) = FLWMAX_soil       ! SAT. HYDR. CONDUCTIVITY, M/S
!	  DLMAX(j)  = DLMAX_soil        ! SAT. WATER DIFFUSIVITY, ?KG/(M*SEC)
!	  WLM0(j)   = WLM0_soil         ! MAXIMAL UNFREEZING WATER AT 0C
!	  WLM7(j)   = WLM7_soil         ! MAXIMAL UNFREEZING WATER AT T<<0C
!	 end do
	end if
      
	if (firstcall) firstcall=.false.   
      RETURN
      END

      SUBROUTINE PBLDAT

!     PBLDAT assigns values to parameters used in surface layer parameterization
!     (see subr. dragvl and subr., called by dragvl)

      IMPLICIT real (A-H,O-Z)
C
C*=====================================================================
C*       INITIALIZATION OF PBL COMMON BLOCKS (FBL,DRAG)               =
C*    FBL(1)     :  FOR LINEAR EXTRAPOLATION OF WIND ONTO SURFACE     =
C*    FBL(2)     :  ---------------------------------------------     =
C*    FBL(3)     :  PRESCRIBED MINIMAL VALUE OF SURFACE WIND MODULE   =
C*    FBL(4)     :  HEIGHT OF CONSTANT FLUX LAYER ( M )               =
C*    FBL(5)     :  CRITICAL VALUE OF RELATIVE HUMIDITY IN            =
C*               :  CALCULATION OF EQUIVALENT POTENTIAL TEMPERATURE   =
C*    FBL(6)     :  SEA WATER TEMPERATURE UNDER ICE ( DEG. K )        =
C*    FBL(7)     :  PRECRIBED MINIMAL VALUE OF SEA ICE SURFACE        =
C*                  TEMPERATURE ( DEG. K )                            =
C*    FBL(8)     :  RESERVED                                          =
C*    FBL(9)     :  PBL WIND TURNING ANGLE OVER OCEAN ( DEG )         =
C*    FBL(10)    :  --------------------------- ICE ----------------- =
C*    FBL(11)    :  --------------------------- SNOW COVERED SOIL --- =
C*    FBL(12)    :  --------------------------- BARED SOIL ---------- =
C*    FBL(13)    :  --------------------------- IN TROPICS ---------- =
C*    FBL(14)    :  BOUNDARIES OF TROPICS ( RADIANS )                 =
C*    FBL(15)    :  RESERVED
C*    FBL(16)    :  (SNOW COVERED SOIL HEAT CONDUCTIVITY)/DEPTH ----- =
C*    FBL(17)    :  (SEA ICE HEAT CONDUCTIVITY)/(SEA ICE DEPTH) ----- =
C*    FBL(18)    :  (SOIL HEAT CONDUCTIVITY)/DEPTH, CAL/(DEG.*M**2*S) =
C*    FBL(19)    :  BULK SOIL HEAT CAPACITY, CAL/(DEG.*M**2)          =
C*    FBL(20)    :  (SOIL MOISTURE CONDUCTIVITY)/DEPTH,SEC**-1        =
C*=====================================================================
      COMMON /FBL/ FBL(20)
      COMMON /BL1/ AZ0P,VEG,WL,DZZ
      COMMON /HYDR/ SNCR,WLMMX,WSMAX,WSL,WSG,CEFF,
     &              CA,DZL,DZG,BMIN,WSDL,WSDG,DMIN,DMAX,D,HR,WIINF,
     &              ZRM,ZRMM,TRM,FLXMIN,TOMIN,HSNold
      COMMON /VEGSW/ CCQ,CK,SWW,TL
      common /p_drag/ vkc,anu,z0min,chc,p_m,p_t,q_m,q_t,
     &                alpha_m,alpha_t,beta_m,beta_t,gamma_m,gamma_t,
     &                zeta_st,a_m,a_t,b_m,b_t
	real FBL,AZ0P,VEG,WL,DZZ,SNCR,WLMMX,WSMAX,WSL,WSG,CEFF,
     &CA,DZL,DZG,BMIN,WSDL,WSDG,DMIN,DMAX,D,HR,WIINF,ZRM,ZRMM,TRM,
     &FLXMIN,TOMIN,CCQ,CK,SWW,TL,
     &vkc,anu,z0min,chc,p_m,p_t,q_m,q_t,alpha_m,alpha_t,beta_m,
     &beta_t,gamma_m,gamma_t,zeta_st,a_m,a_t,b_m,b_t

      FBL(1) = 1.
      FBL(2) = 1.
      FBL(3) = 1.
      FBL(4) = 70.
      FBL(5) = 1.
      FBL(6) = 271.5
      FBL(7) = 243.2
      FBL(8) = 0.
      FBL(9) = 20.
      FBL(10) = 10.
      FBL(11) = 30.
      FBL(12) = 30.
      FBL(13) = 0.
      FBL(14) = 20.*3.1415926536/180.
C     FBL(8)...FBL(14) WILL BE NOT USED IN THE MULTILEVEL PBL
      FBL(15) = 0.
      FBL(16) = 0.04/2.
      FBL(17) = 0.5/3.
      FBL(18) = 0.6/2.
      FBL(19) = 4.0E+04
      FBL(20) = 1.0E-6/5.
c------------------------------------------------------------------
c*    subroutine dragvl constants
c------------------------------------------------------------------
      vkc = 0.4
      anu = 0.000015
      z0min = 1.5e-5
      chc = 0.0132
      zeta_st = -2.
      p_m = - 0.25
      p_t = - 0.5
      q_m = - 1./3.
      q_t = - 1./3.
      alpha_m = 1.
      alpha_t = 1.
      beta_m = 4.7
      beta_t = 4.7
      gamma_m = 16.
      gamma_t = 16.
      a_m = (1. + gamma_m*(p_m/q_m - 1.)*zeta_st)
     &      *(1. - gamma_m*zeta_st)**(p_m/q_m - 1.)
      b_m = a_m*gamma_m*(p_m/q_m)/(1. + gamma_m*(p_m/q_m - 1.)*zeta_st)
      a_t = (1. + gamma_t*(p_t/q_t - 1.)*zeta_st)
     &      *(1. - gamma_t*zeta_st)**(p_t/q_t - 1.)
      b_t = a_t*gamma_t*(p_t/q_t)/(1. + gamma_t*(p_t/q_t - 1.)*zeta_st)

      RETURN
      END


      REAL(8) FUNCTION unfrwat(T,i)
	use driving_params
	use arrays
	
!    CALCULATION OF LIQUID WATER CONTENT IN FREEZING SOIL
!    T - DEG. C; WLM0,WLM7,WLMAX - KG/KG

      implicit none
	real(8) T
	integer(4) i

      unfrwat = (WLM0(i)-WLM7(i))*dexp(T/3.) + WLM7(i)
      !unfrwat = 0 
      
      RETURN
      END


	REAL(8) FUNCTION QS(phase,t,p)
	
!     QS - specific humidity, kg/kg, for saturated water vapour

	implicit none
	real(8) t,p,aw,bw,ai,bi,a,b,es
	integer(4) phase
!phase = 1 - liquid, 0 - ice
      aw = 7.6326    
      bw = 241.9    
      ai = 9.5  
      bi = 265.5 
	a  = phase*aw+(1-phase)*ai
	b  = phase*bw+(1-phase)*bi
      es=610.7*10.**(a*t/(b+t))
      QS=0.622*es/p
	END FUNCTION


	REAL(8) FUNCTION ES(phase,t,p)
	
!     ES is pressure of saturated water vapour, Pa	

	implicit none
	real(8) t,p,aw,bw,ai,bi,a,b
	integer(4) phase
!phase = 1 - liquid, 0 - ice
      aw = 7.6326    
      bw = 241.9    
      ai = 9.5  
      bi = 265.5 
	a  = phase*aw+(1-phase)*ai
	b  = phase*bw+(1-phase)*bi
      ES = 610.7*10.**(a*t/(b+t))
      END FUNCTION
      
      REAL(8) FUNCTION MELTPNT(C)
      implicit none
      real(8) C, dtdc
      
      SAVE
      
      dtdc = 66.7
      
      MELTPNT = 0. - C*dtdc
      END FUNCTION MELTPNT

      REAL(8) FUNCTION WL_MAX(por,rosdry,wi)
	implicit none
	real(8) por,rosdry,wi,row,roi
	row=1000.
	roi=917.
	wl_max=por*row/(rosdry*(1-por)) - wi*row/roi
	END FUNCTION WL_MAX

      REAL(8) FUNCTION WI_MAX(por,rosdry)
	implicit none
	real(8) por,rosdry,roi
	roi=917.
	wi_max=por*roi/(rosdry*(1-por))
	END FUNCTION WI_MAX
	
	SUBROUTINE progonka (a, b, c, f, y, N)
	implicit none
	
!      FACTORIZATION METHOD FOR THE FOLLOWING SYSTEM OF LINEAR EQUATIONS:
!     -a(i)*y(i-1)+c(i)*y(i)-b(i)*y(i+1)=f(i) i=2,N-1
!      c(1)*y(1)-b(1)*y(2)=f(1)
!     -a(N)*y(N-1)+c(N)*y(N)=f(N)
!
      integer(4), parameter:: M=350
      integer(4) N, i
      real(8) a(M), b(M), c(M), f(M), y(M)
	real(8) alpha(M+2), beta(M+2) 

	SAVE
	 	
	!y=0
      alpha(2) = b(1)/c(1)
	beta(2) = f(1)/c(1)
	do i = 3, N+1 
       alpha(i) = b(i-1)/(c(i-1)-a(i-1)*alpha(i-1))
	 beta(i) = (f(i-1)+a(i-1)*beta(i-1))/ 
     & (c(i-1)-a(i-1)*alpha(i-1))
      end do
	y(N) = beta(N+1)
	do i = N-1, 1, -1
	 y(i) = alpha(i+1)*y(i+1)+beta(i+1)
      end do
	!write (*, *) (c(i)*y(i)-b(i)*y(i+1)-a(i)*y(i-1)-f(i),i=2,N-1) !/f(10)
 
      END 
      
	SUBROUTINE progonka2_db (a, b, c, f, y, K, N)
	implicit none
!FACTORIZATION METHOD FOR THE FOLLOWING SYSTEM OF LINEAR EQUATIONS:
!     -a(i)*y(i-1)+c(i)*y(i)-b(i)*y(i+1)=f(i) i=K+1,N-1
!      c(K)*y(K)-b(K)*y(K+1)=f(K)
!     -a(N)*y(N-1)+c(N)*y(N)=f(N)
!
      integer(4), parameter:: M=350
      integer(4) K, N, i
      real(8) a(M), b(M), c(M), f(M), y(M)
	real(8) alpha(M+2), beta(M+2) 
	SAVE
				
      alpha(K+1) = b(K)/c(K)
	beta(K+1) = f(K)/c(K)
	do i = K+2, N+1 
       alpha(i) = b(i-1)/(c(i-1)-a(i-1)*alpha(i-1))
	 beta(i) = (f(i-1)+a(i-1)*beta(i-1))/ 
     & (c(i-1)-a(i-1)*alpha(i-1))
      end do
	y(N) = beta(N+1)
	do i = N-1, K, -1
	 y(i) = alpha(i+1)*y(i+1)+beta(i+1)
      end do
	 
      END 

	SUBROUTINE progonka2(a, b, c, f, y, K, N)
	implicit none
!FACTORIZATION METHOD FOR THE FOLLOWING SYSTEM OF LINEAR EQUATIONS:
!     -a(i)*y(i-1)+c(i)*y(i)-b(i)*y(i+1)=f(i) i=K+1,N-1
!      c(K)*y(K)-b(K)*y(K+1)=f(K)
!     -a(N)*y(N-1)+c(N)*y(N)=f(N)
!
      integer(4), parameter:: M=350
      integer(4) K, N, i
      real(8) a(M), b(M), c(M), f(M), y(M)
	real(8) alpha(M+2), beta(M+2) 
	SAVE
				
      alpha(K+1) = b(K)/c(K)
	beta(K+1) = f(K)/c(K)
	do i = K+2, N+1 
       alpha(i) = b(i-1)/(c(i-1)-a(i-1)*alpha(i-1))
	 beta(i) = (f(i-1)+a(i-1)*beta(i-1))/ 
     & (c(i-1)-a(i-1)*alpha(i-1))
      end do
	y(N) = beta(N+1)
	do i = N-1, K, -1
	 y(i) = alpha(i+1)*y(i+1)+beta(i+1)
      end do
	 
      END 


      SUBROUTINE ind_stab_fact (a,b,c,N,M,ind_stab,ind_bound)
      
!     IND_STAB_FACT checks, if the stability condition for factorization method
!     is satisfied      
	
	implicit none
	real(8), dimension(1:350):: a,b,c
      integer(4) M,i,N
      logical ind_stab, ind_bound 

	SAVE

	ind_stab=.true.
 	if (ind_bound==.true.) then 
	 if (dabs(b(N))>=dabs(c(N)).or.dabs(a(M))>=dabs(c(M))) then
	  ind_stab=.false.
	  RETURN
       endif
      endif
	do i=N+1,M-1
	 if (dabs(a(i))+dabs(b(i))>=dabs(c(i))) then
	  ind_stab=.false.
	  RETURN
	 endif
	enddo
	
	END

	
	SUBROUTINE ind_stab_fact_db (a,b,c,N,M,ind_stab,ind_bound)

	implicit none
	real(8), dimension(1:350):: a,b,c
      integer(4) M,i,N
      logical ind_stab, ind_bound 

	SAVE

	ind_stab=.true.
	if (ind_bound==.true.) then 
	 if (dabs(b(N))>=dabs(c(N)).or.dabs(a(M))>=dabs(c(M))) then
	  ind_stab=.false.
	  RETURN
       endif
      endif
	do i=N+1,M-1
	 if (dabs(a(i))+dabs(b(i))>=dabs(c(i))) then
	  ind_stab=.false.
	  RETURN
	 endif
	enddo
	
	END

	SUBROUTINE CUT_E_DB(E,N,porog)
	implicit none
      real(8) E(350),porog
	integer(4) N,i
!	porog = 10.**(-9) 
	do i=1,N 
	 if (E(i)<porog) E(i)=porog
  	enddo
      END SUBROUTINE CUT_E_DB

      SUBROUTINE CUT_EPS_DB(eps,N,porog)
	implicit none
      real(8) eps(350),porog
	integer(4) N,i
!	porog = 10.**(-11) 
	do i=1,N
	 if (eps(i)<porog) eps(i)=porog
  	enddo
      END SUBROUTINE CUT_EPS_DB
 
      
      SUBROUTINE MATRIXSUM(a,b,c,k)
	implicit none
	
!     MATRIXES: C=A+B
	real(8), dimension (350,2,2)::a,b,c
	integer(4) k,j,i

      do i=1,2
	 do j=1,2
	  c(k,i,j)=a(k,i,j)+b(k,i,j)
       enddo
      enddo
	
	END SUBROUTINE MATRIXSUM

	SUBROUTINE MATRIXMULT(a,b,c,k)
	implicit none
	
!     MATRIXES: C=A*B      

	real(8), dimension (350,2,2)::a,b,c
	integer(4) k

	c(k,1,1)=a(k,1,1)*b(k,1,1)+a(k,1,2)*b(k,2,1)
	c(k,1,2)=a(k,1,1)*b(k,1,2)+a(k,1,2)*b(k,2,2)
	c(k,2,1)=a(k,2,1)*b(k,1,1)+a(k,2,2)*b(k,2,1)
	c(k,2,2)=a(k,2,1)*b(k,1,2)+a(k,2,2)*b(k,2,2)

	END SUBROUTINE MATRIXMULT


	SUBROUTINE MATRIXMULTVECTOR(a,g,f,k)
	implicit none
	
!     MATRIX A, VECTORS g, f: Ag=f

	real(8) a(350,2,2),f(350,2),g(350,2)
	integer(4) k

      f(k,1)=a(k,1,1)*g(k,1)+a(k,1,2)*g(k,2)
	f(k,2)=a(k,2,1)*g(k,1)+a(k,2,2)*g(k,2)

      END SUBROUTINE MATRIXMULTVECTOR

	SUBROUTINE VECTORSUM(a,b,c,k)
	implicit none
	
!     VECTORS: C=A+B

      real(8), dimension(350,2)::a,b,c
	integer(4) k
	c(k,1)=a(k,1)+b(k,1)
      c(k,2)=a(k,2)+b(k,2)
	END SUBROUTINE VECTORSUM


	SUBROUTINE INVERSMATRIX(a,a1,k)
	implicit none 
	
!     MATRIXES: A1=A*(-1)

      real(8), dimension(350,2,2)::a,a1
	integer(4) k
	
	a1(k,1,1)=a(k,2,2)/(a(k,1,1)*a(k,2,2)-a(k,1,2)*a(k,2,1)) 
	a1(k,1,2)=-a(k,1,2)/(a(k,1,1)*a(k,2,2)-a(k,1,2)*a(k,2,1))
	a1(k,2,1)=-a(k,2,1)/(a(k,1,1)*a(k,2,2)-a(k,1,2)*a(k,2,1))
	a1(k,2,2)=a(k,1,1)/(a(k,1,1)*a(k,2,2)-a(k,1,2)*a(k,2,1))  

	END SUBROUTINE INVERSMATRIX

	
	SUBROUTINE MATRIXPROGONKA(a,b,c,d,y,N)
	
!     MATRIXPROGONKA solves the set of MATRIX three-point diference equations	
	
	implicit none
	real(8), dimension(350,2,2):: a,b,c,x3,x32,x31,alpha,x33
	real(8), dimension(350,2):: y,d,x2,beta,x21
	real(8), dimension(2,2):: a2,b2,c2
	integer(4) N,i,j,k

      call INVERSMATRIX(c,x32,1)
	call MATRIXMULT(x32,b,x3,1)
	do i=1,2
	 do j=1,2
	  alpha(2,i,j)=x3(1,i,j)
	 enddo
      enddo
      call MATRIXMULTVECTOR(x32,d,x2,1)
	do i=1,2
	 beta(2,i)=x2(1,i)
      enddo

	do k=3,N
	 CALL MATRIXMULT(A,ALPHA,X3,k-1)
	 CALL MATRIXSUM(C,-X3,X31,k-1)
	 CALL INVERSMATRIX(X31,X32,k-1)
	 CALL MATRIXMULT(X32,B,X3,k-1)
	 do i=1,2
	  do j=1,2
	   alpha(k,i,j)=X3(k-1,i,j)
	  enddo
       enddo
	 !call matrixmult(x3,x31,x33,k-1)
	 !call matrixsum(-b,x33,x3,k-1)
	 CALL MATRIXMULTVECTOR(A,BETA,X2,K-1)
	 CALL VECTORSUM(D,X2,X21,K-1)
	 CALL MATRIXMULTVECTOR(X32,X21,X2,K-1)
	 do i=1,2
	  beta(k,i)=X2(k-1,i)
       enddo
      enddo

      CALL MATRIXMULT(A,ALPHA,X3,N)
	CALL MATRIXSUM(C,-X3,X31,N)
	CALL INVERSMATRIX(X31,X32,N)
	CALL MATRIXMULTVECTOR(A,BETA,X2,N)
	CALL VECTORSUM(D,X2,X21,N)
	CALL MATRIXMULTVECTOR(X32,X21,X2,N)
	do i=1,2
	 Y(N,i)=X2(N,i)
      enddo

	do k=N-1,1,-1
	 CALL MATRIXMULTVECTOR(ALPHA,Y,X2,K+1)
	 CALL VECTORSUM(X2,BETA,X21,K+1)
	 Y(K,1)=X21(K+1,1)
	 Y(K,2)=X21(K+1,2)
      enddo

	END SUBROUTINE MATRIXPROGONKA


      REAL(8) FUNCTION KRON(i,j)
      implicit none
	integer(4) i,j
	kron=0.
	if(i==j) kron=1.
	END FUNCTION 
	
      
!	SUBROUTINE WaterMix(flag, T1, T2, ro1)
!	use driving_params
!	use arrays
!	implicit none
!                 
!      real(8), allocatable:: T1(:), T2(:), ro1(:), 
!     & ro2(:)
!      integer(4) n, i, k, p, flag
!	logical firstcall

!	data firstcall /.true./

!	SAVE

!	if (firstcall) then
!       allocate (T1(1:M+1), T2(1:M+1), ro1(1:M+1), 
!     & ro2(1:M+1))
!	endif
                            
!      do i = 1, M+1
!       ro1(i) = 800.969d-7+588.194d-7*T1(i)- 
!     & 811.465d-8*(T1(i)**2)+476.600d-10*(T1(i)**3)
!       ro2(i) = ro1(i)
!      end do    
                       
!      do i = 1, M+1
!       T2(i) = T1(i)
!      end do        
                   
!	!if (flag==0) then			      
!      ! do n = 1, M-1 
!      !  k = 1
!      !  do i = 1, M+1-n
!      !   if (ro1(i) > ro1(k)) then
!      !    k = i 
!      !   end if
!      !  end do
!      !  T2(M+1-n) = T1(k)
!      !  T2(k) = T1(M+1-n)
!      !  ro2(M+1-n) = ro1(k)
!      !  ro2(k) = ro1(M+1-n) 
!      !  do p = 1, M+1-n
!      !   ro1(p) = ro2(p)
!      !   T1(p) = T2(p)
!      !  end do
!      ! end do
!	!else       
!       do n = 1, M-2 
!        k = 2
!        do i = 2, M+1-n
!         if (ro1(i) > ro1(k)) then
!          k = i 
!         end if
!        end do
!        T2(M+1-n) = T1(k)
!        T2(k) = T1(M+1-n)
!        ro2(M+1-n) = ro1(k)
!        ro2(k) = ro1(M+1-n) 
!        do p = 1, M+1-n
!         ro1(p) = ro2(p)
!         T1(p) = T2(p)
!        end do
!       end do      
!      !endif
	
!      if (flag==0) then
!       ro1(1) = 800.969d-7+588.194d-7*(T1(1)+T1(2))/2.- 
!     & 811.465d-8*((T1(1)+T1(2))**2)/4.+476.600d-10*
!     & ((T1(1)+T1(2))**3)/8.
!       ro2(1) = ro1(1)
!      c1:do i=M,2,-1
!	  if (ro2(1)>ro2(i)) then
!	   T2(i)=(T1(1)+T1(2))/2
!	   ro2(i)=ro1(1)
!	   do k=2,i
!	    T2(k-1)=T1(k)
!	    ro2(k-1)=ro1(k)
!         enddo
!	   T1=T2
!	   ro1=ro2
!	   exit c1
!	  endif
!	 enddo c1    
!      endif

!	if (firstcall) firstcall=.false.                                       
!      END SUBROUTINE WaterMix 

	function fext2(fname,ext)
      character*40 fname,fext2
      character*3 ext
      parameter(n=40)

      fext2=fname

      do 10 i=30,1,-1
         if(fext2(i:i).eq.'.') then
            fext2(i:n)=char(0)
            go to 11
         endif
10    continue
11    continue

      do 20 i=30,1,-1
         if(fext2(i:i).ne.' ' .and. fext2(i:i).ne.char(0)) then
            fext2(i+1:i+4)='.'//ext
            go to 21
         endif
20    continue
      write(*,*) 'erro em fext'

21    continue
      return
      end
                 
      subroutine readgrd2(nunit,var,i0,i1,j0,j1)
c
c reads 2d array from a grd file
c
      implicit real*8(a-h,o-z)
      dimension var(i0:i1,j0:j1)
      character*4 dsaa

      read(nunit,'(a4)') dsaa
      read(nunit,*) nx,ny
      read(nunit,*) xmin,xmax
      read(nunit,*) ymin,ymax
      read(nunit,*) zmin,zmax

      
      do j=j0,j1
         read(nunit,*) (var(i,j),i=i0,i1)
	enddo
      

      return
      end


	SUBROUTINE MON_OUT(dt)
      use arrays
      implicit none
      real(8) year,month,day,hour,dt
      real(8) h1,l1,hs1,ls1
      real(8), allocatable:: Tw_m(:),Tw_sum(:),Sal_m(:),Sal_sum(:),
     & TKE_m(:),TKE_sum(:),eps_m(:),eps_sum(:)
      integer(4) ndm1,ndm2,i,ndaym(12),month_old
      character month1*2,year1*4,day1*2,hour1*2,outpath*60
      character*6 timestring
      logical firstcall
      common /ymdh/ year,month,day,hour
      common /layers/ h1,l1,hs1,ls1
      common /out/ outpath
      data firstcall /.true./
      data ndaym /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/

      SAVE

      if (firstcall) then
       month_old=int4(month)
       ndm1 = int4(day)
       allocate (Tw_m(1:M+1),Tw_sum(1:M+1))
       allocate (Sal_m(1:M+1),Sal_sum(1:M+1))
	 allocate (TKE_m(1:M+1),TKE_sum(1:M+1))
	 allocate (eps_m(1:M+1),eps_sum(1:M+1))
       Tw_sum = 0.
       Sal_sum = 0.
	 TKE_sum = 0.
	 eps_sum = 0.
      endif

      Tw_sum=Tw_sum+Tw1
      Sal_sum=Sal_sum+Sal1
	TKE_sum=TKE_sum+E1
	eps_sum=eps_sum+eps1
      if (month/=float(month_old)) then
       Tw_m=Tw_sum/(((ndaym(month_old)-ndm1+1)*24.*60.*60.)/dt)
       Sal_m=Sal_sum/(((ndaym(month_old)-ndm1+1)*24.*60.*60.)/dt)
	 TKE_m=TKE_sum/(((ndaym(month_old)-ndm1+1)*24.*60.*60.)/dt)
	 eps_m=eps_sum/(((ndaym(month_old)-ndm1+1)*24.*60.*60.)/dt)
	 call dataminus (1,year1,month1,day1,hour1)
	 call timestr(6,year1,month1,day1,hour1,timestring)
       open(600,file=outpath(1:len_trim(outpath))//'monthly/'//
     &  'Profiles'//timestring//'.dat',status='unknown')
	 write (600,*) '1 - depth, m' 
	 write (600,*) '2 - temperature, C'
       write (600,*) '3 - salinity, kg/kg' 
       write (600,*) '4 - turbulent kinetic energy, m**2/s**2'
	 write (600,*) '5 - disspation rate, m**2/s**3'
       do i=1,M+1 
        write (600,'(f6.2,f9.3,3e12.4)') 
     &   (i-1)/float(M)*h1,Tw_m(i),Sal_m(i),TKE_m(i),eps_m(i)
       enddo
       close(600)
       month_old=int4(month)
       ndm1=int4(day) 
       Tw_sum = 0.
       Sal_sum = 0.
	 TKE_sum = 0.
	 eps_sum = 0.
      endif

      if (firstcall) firstcall=.false.
      END SUBROUTINE MON_OUT


      SUBROUTINE DAY_OUT(dt)
      use arrays
      implicit none
      real(8) year,month,day,hour,dt,day_dur
      real(8) h1,l1,hs1,ls1
      real(8), allocatable:: Tw_m(:),Tw_sum(:),Sal_m(:),Sal_sum(:),
     & TKE_m(:),TKE_sum(:),eps_m(:),eps_sum(:)
      integer(4) i,nday_old
      character month1*2,year1*4,day1*2,hour1*2,outpath*60
      character*8 timestring
      logical firstcall, firstwrite
      common /ymdh/ year,month,day,hour
      common /layers/ h1,l1,hs1,ls1
      common /out/ outpath
      data firstcall, firstwrite, day_dur /.true., .true., 86400./
      
      SAVE

      if (firstcall) then
       nday_old = int4(day)
       allocate (Tw_m(1:M+1),Tw_sum(1:M+1))
       allocate (Sal_m(1:M+1),Sal_sum(1:M+1))
	 allocate (TKE_m(1:M+1),TKE_sum(1:M+1))
	 allocate (eps_m(1:M+1),eps_sum(1:M+1))
       Tw_sum = 0.
       Sal_sum = 0.
	 TKE_sum = 0.
	 eps_sum = 0.
      endif

      Tw_sum=Tw_sum+Tw1
      Sal_sum=Sal_sum+Sal1
	TKE_sum=TKE_sum+E1
	eps_sum=eps_sum+eps1
      if (int4(day)/=nday_old) then
       Tw_m =Tw_sum /(day_dur/dt)
       Sal_m=Sal_sum/(day_dur/dt)
	 TKE_m=TKE_sum/(day_dur/dt)
	 eps_m=eps_sum/(day_dur/dt)
       if (.not. firstwrite) then
        call dataminus (2,year1,month1,day1,hour1)
        call timestr(8,year1,month1,day1,hour1,timestring)
        open(600,file=outpath(1:len_trim(outpath))//'daily/'// 
     &   'Profiles'//timestring//'.dat',status='unknown')
        write (600,*) '1 - depth, m' 
	  write (600,*) '2 - temperature, C'
        write (600,*) '3 - salinity, kg/kg' 
        write (600,*) '4 - turbulent kinetic energy, m**2/s**2'
	  write (600,*) '5 - disspation rate, m**2/s**3'
        do i=1,M+1 
         write (600,'(f6.2,f9.3,3e12.4)') 
     &    (i-1)/float(M)*h1,Tw_m(i),Sal_m(i),TKE_m(i),eps_m(i)
        enddo
        close(600)
       endif 
       nday_old=int4(day) 
       Tw_sum = 0.
       Sal_sum = 0.
	 TKE_sum = 0.
	 eps_sum = 0.
       if (firstwrite) firstwrite=.false.
      endif

      if (firstcall) firstcall=.false.
      END SUBROUTINE DAY_OUT


	SUBROUTINE HOUR_OUT(dt)
      use arrays
      use driving_params
      implicit none
      real(8) year,month,day,hour,dt,hour_dur
      real(8) h1,l1,hs1,ls1
      real(8), allocatable:: Tw_m(:),Tw_sum(:),Sal_m(:),Sal_sum(:),
     & TKE_m(:),TKE_sum(:),eps_m(:),eps_sum(:)
      integer(4) i,nhour_old
      character month1*2,year1*4,day1*2,hour1*2,outpath*60
      character*10 timestring
      logical firstcall, firstwrite
      common /ymdh/ year,month,day,hour
      common /layers/ h1,l1,hs1,ls1
      common /out/ outpath
      data firstcall, firstwrite, hour_dur /.true., .true., 3600./
      
      SAVE

      if (firstcall) then
       nhour_old = int4(hour)
       allocate (Tw_m(1:M+1),Tw_sum(1:M+1))
       allocate (Sal_m(1:M+1),Sal_sum(1:M+1))
	 allocate (TKE_m(1:M+1),TKE_sum(1:M+1))
	 allocate (eps_m(1:M+1),eps_sum(1:M+1))
       Tw_sum = 0.
       Sal_sum = 0.
	 TKE_sum = 0.
	 eps_sum = 0.
      endif

      Tw_sum=Tw_sum+Tw1
      Sal_sum=Sal_sum+Sal1
	TKE_sum=TKE_sum+E1
	eps_sum=eps_sum+eps1
      if (int4(hour)/=nhour_old) then
       Tw_m=Tw_sum/(interval*hour_dur/dt)
       Sal_m=Sal_sum/(interval*hour_dur/dt)
	 TKE_m=TKE_sum/(interval*hour_dur/dt)
	 eps_m=eps_sum/(interval*hour_dur/dt)
       if (.not. firstwrite) then
        call dataminus (3,year1,month1,day1,hour1)
        call timestr(10,year1,month1,day1,hour1,timestring)
        open(600, file=outpath(1:len_trim(outpath))//'hourly/'// 
     &   'Profiles'//timestring//'.dat',
     &    status='unknown')
	  write (600,*) '1 - depth, m' 
	  write (600,*) '2 - temperature, C'
        write (600,*) '3 - salinity, kg/kg' 
        write (600,*) '4 - turbulent kinetic energy, m**2/s**2'
	  write (600,*) '5 - disspation rate, m**2/s**3'
        do i=1,M+1 
         write (600,'(f6.2,f9.3,3e12.4)') 
     &    (i-1)/float(M)*h1,Tw_m(i),Sal_m(i),TKE_m(i),eps_m(i)
        enddo
        close(600)
       endif 
       nhour_old=int4(hour) 
       Tw_sum = 0.
       Sal_sum = 0.
	 TKE_sum = 0.
	 eps_sum = 0.
       if (firstwrite) firstwrite=.false.
      endif

      if (firstcall) firstcall=.false.
      END SUBROUTINE HOUR_OUT          
      
      SUBROUTINE EVERYSTEP_OUT
      use arrays
      implicit none
      real(8) h1,l1,hs1,ls1
      integer(4) i
      character outpath*60
      logical firstcall
      data firstcall /.true./
      common /out/ outpath
      common /layers/ h1,l1,hs1,ls1
      
      SAVE
      
      if (firstcall) then
       open (605,file=outpath(1:len_trim(outpath))//'everystep/'// 
     &  'Profiles.dat',  status='unknown')
       write (605,*) '1 - depth, m' 
	 write (605,*) '2 - temperature, C'
       write (605,*) '3 - salinity, kg/kg' 
       write (605,*) '4 - turbulent kinetic energy, m**2/s**2'
	 write (605,*) '5 - disspation rate, m**2/s**3'
      endif
      
      write (605,*) 'nstep = ', nstep
      do i=1,M+1 
       write (605,'(f6.2,f9.3,3e12.4)') 
     &  (i-1)/float(M)*h1,Tw1(i),Sal1(i),E1(i),eps1(i)
      enddo
      
      if (firstcall) firstcall=.false.
      END SUBROUTINE EVERYSTEP_OUT
      
      
      SUBROUTINE series_out(tsw)
      use driving_params
      use atmos
      implicit none
      real(8) tsw,year,month,day,hour,h1,l1,hs1,ls1
      integer(4) n_out,n_count,nhour_old
      character*60 outpath
      logical firstcall
      data firstcall /.true./
      common /out/ outpath
      common /ymdh/ year,month,day,hour
      common /layers/ h1,l1,hs1,ls1
      
      SAVE
            
      if (firstcall) then
!       if (dt_out/interval/=int4(dt_out/interval)) then
!        print*, 'Illegal dt_out: STOP'
!        stop
!       endif
!       n_out = int4(dt_out/interval) 
       open (603,file=outpath(1:len_trim(outpath))//'time_series/'// 
     &  'layers.dat',  status='unknown')
	 write (603,*) '1 - year'
       write (603,*) '2 - month'
	 write (603,*) '3 - day'
	 write (603,*) '4 - hour'
	 write (603,*) '5 - water layer thickness, m'
	 write (603,*) '6 - ice layer thickness,   m'
	 write (603,*) '7 - snow layer thickness,  m'
	 write (603,*) '8 - bottom ice thickness,  m'
       open (604,file=outpath(1:len_trim(outpath))//'time_series/'// 
     &  'T_fluxes.dat',status='unknown')
	 write (604,*) '1 - year'
       write (604,*) '2 - month'
	 write (604,*) '3 - day'
	 write (604,*) '4 - hour'
	 write (604,*) '5 - surface temperature,  C'
	 write (604,*) '6 - sensible heat flux,   W/m**2'
	 write (604,*) '7 - latent heat flux,     W/m**2'
	 nhour_old = int4(hour)
      endif
      
      if (int4(hour)/=nhour_old) then
       write (603,'(4f5.0,4f9.4)') year,month,day,hour,h1,l1,hs1,ls1
       write (604,'(4f5.0,3f9.4)') year,month,day,hour,tsw,hw,xlew
       nhour_old = int4(hour)
      endif 
      
      if (firstcall) firstcall=.false.
      END SUBROUTINE series_out           
      
      
      SUBROUTINE dataminus(n,year1,month1,day1,hour1)
      use driving_params
      implicit none
      real(8) year,month,day,hour
      integer(4) n,ndaym(12)
      character year1*4,month1*2,day1*2,hour1*2
      common /ymdh/ year,month,day,hour
      data ndaym /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
      
      SAVE
      
      if (n==1) then
!     Minus 1 month      
       if (int4(month)==1) then
        write (month1,'(i2)') 12
        write (year1, '(i4)') int4(year)-1
       else
        write (month1,'(i2)') int4(month)-1
        write (year1, '(i4)') int4(year)
       endif 
       elseif (n==2) then
!     Minus 1 day      
       if (int4(day)==1) then
        if (int4(month)==1) then
         write (month1,'(i2)') 12
         write (year1, '(i4)') int4(year)-1
         write (day1, '(i2)') 31
        else
         write (month1,'(i2)') int4(month)-1
         write (year1, '(i4)') int4(year)
         write (day1, '(i2)') ndaym(int4(month)-1) 
        endif
       else
        write (month1,'(i2)') int4(month)
        write (year1, '(i4)') int4(year)
        write (day1, '(i2)')  int4(day)-1 
       endif
      elseif (n==3) then
!     Minus 1 hour      
       if (int4(hour)==0) then
        if (int4(day)==1) then
         if (int4(month)==1) then
          write (month1,'(i2)') 12
          write (year1, '(i4)') int4(year)-1
          write (day1, '(i2)')  31
          write (hour1, '(i2)') 24-int4(interval)
         else
          write (month1,'(i2)') int4(month)-1
          write (year1, '(i4)') int4(year)
          write (day1, '(i2)')  ndaym(int4(month)-1) 
          write (hour1, '(i2)') 24-int4(interval)
         endif
        else
         write (month1,'(i2)') int4(month)
         write (year1, '(i4)') int4(year)
         write (day1,  '(i2)') int4(day)-1
         write (hour1, '(i2)') 24-int4(interval)
        endif
       else
        write (month1,'(i2)') int4(month)
        write (year1, '(i4)') int4(year)
        write (day1,  '(i2)') int4(day)
        write (hour1, '(i2)') int4(hour)-int4(interval)
       endif
      else
       print*, 'Invalid meaning of n in dataminus: STOP'
       stop
      endif
      
      END SUBROUTINE dataminus
      
      SUBROUTINE timestr(n,y,m,d,h,timestring)
      implicit none
      integer n,i
      character y*4,m*2,d*2,h*2
      character(len=n) timestring 
      
      if (n==10) then
       write (timestring, '(a4,3a2)') y,m,d,h 
      elseif (n==8) then
       write (timestring, '(a4,2a2)') y,m,d 
      elseif (n==6) then
       write (timestring, '(a4,a2)') y,m 
      endif 
      
      do i = 1, n
       if (timestring(i:i)==' '.or.timestring(i:i)==char(0))
     &  timestring(i:i)='0'
      enddo
      
      END SUBROUTINE timestr