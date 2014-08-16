 SUBROUTINE ALLOCARRAYS(nx,ny)
 use DRIVING_PARAMS
 use ARRAYS
 use NUMERIC_PARAMS
 use ATMOS, only: &
 & velfrict_prev
 implicit none

 integer(4) :: nx
 integer(4) :: ny
 integer(4) :: i !Loop index
 SAVE

 allocate (ddz(1:M), ddzi(1:Mice), dz_full(1:M), &
 & z_full(1:M+1), z_half (1:M))
 allocate (Tw1(1:M+1),Tw2(1:M+1),Ti1(1:Mice+1),Ti2(1:Mice+1), &
 & Tis1(1:Mice+1),Tis2(1:Mice+1),RS(1:M+1),SR(1:M+1),lamw(1:M), &
 & Sal1(1:M+1),Sal2(1:M+1),lamsal(1:M), SRi(1:Mice), &
 & SRdi(1:Mice))
 allocate (Tsoil1(1:ns+2), rosoil(1:ns+2), csoil(1:ns+2), &
 & lamsoil(1:ns+2), dzs(1:ns+2), dzss(1:ns+2), wi1(1:ns+2), &
 & wl1(1:ns+2), Tsoil2(1:ns+2), zsoil(1:ns+2)) 
 allocate (S(1:M),Gen(1:M),F(1:M),TKE_turb_trans(1:M))
 allocate (KT(1:M+1),k2(1:M+1),u1(1:M+1),v1(1:M+1), w1(1:M+1), &
 & E1(1:M+1), eps1(1:M+1))
 allocate (WLM0(1:ns+2),WLM7(1:ns+2),bH(1:ns+2),PSIMAX(1:ns+2), &
 & POR(1:ns+2),FLWMAX(1:ns+2),DLMAX(1:ns+2))
 allocate (Tsoil4(1:ns+2), Tsoil3(1:ns+2), &
 & wl2(1:ns+2), wl3(1:ns+2), wl4(1:ns+2), wl5(1:ns+2),  &
 & wi2(1:ns+2), wi3(1:ns+2), lammoist(1:ns+2), rosdry(1:ns+2), &
 & filtr(1:ns+2), wsoil(1:ns), Sals1(1:ns), Sals2(1:ns) )
 allocate (WU_(1:M+1),WV_(1:M+1),GAMT(1:M+1),GAMU(1:M+1), &
 &  GAMV(1:M+1),TF(1:M+1),KLengu(1:M+1)) 
 allocate (num(1:M+1)) 
 allocate ( E2(1:M+1),eps2(1:M+1),eps3(1:M+1), &
 & u2(1:M+1),v2(1:M+1),w2(1:M+1), &
 & k1(1:M+1),k3(1:M+1),k4(1:M+1),k5(1:M+1), &
 & u3(1:M+1),v3(1:M+1),C1aup(1:M+1),Re(1:M+1),row(1:M+1),Ri(1:M+1), &
 & uv(1:M+1),E_it1(1:M+1),dres2dE(1:M+1),dres2deps(1:M+1), & 
 & E_it2(1:M+1),C_num(1:M+1),E_it3(1:M+1),Feps(1:M+1), &
 & E_it21(1:M+1),eps_it21(1:M+1),eps_it1(1:M+1), &
 & eps_it2(1:M+1),l(1:M+1),k2_mid(1:M+1),k3_mid(1:M+1), &
 & k4_mid(1:M+1),k5_mid(1:M+1),Um(1:M+1),Vm(1:M+1),E12(1:M+1), &
 & eps12(1:M+1),knum(1:M+1),k2t(1:M+1),veg_sw(1:M+1) )
 allocate ( Eeps(1:M+1),Ecent(1:M+1), &
 & Epscent(1:M+1),res_E(1:M+1),res_eps(1:M+1), &
 & dresd(2,2,1:M+1,1:M+1) ) 
 allocate (KC(1:M+1),KLengT(1:M+1),RSR(1:M+1), &
 & l1_2d(nx,ny), h1_2d(nx,ny), ls1_2d(nx,ny), &
 & dep_2d(0:nx-1,0:ny-1), cdm2(0:nx,0:ny))
 allocate (u_2d(M+1,nx,ny),v_2d(M+1,nx,ny),E_2d(M+1,nx,ny), &
 & eps_2d(M+1,nx,ny),Tsoil1_2d(ns+1,nx,ny),wi1_2d(ns+1,nx,ny), &
 & wl1_2d(ns+1,nx,ny),Tw1_2d(M+1,nx,ny),Sal1_2d(M+1,nx,ny), &
 & Ti1_2d(Mice+1,nx,ny), Tis1_2d(Mice+1,nx,ny),Sals1_2d(ns,nx,ny), &
 & hs1_2d(nx,ny),dz_2d(ms,nx,ny),T_2d(ml,nx,ny),wl_2d(ml,nx,ny), &
 & dens_2d(ms,nx,ny),time_2d(nx,ny),dhwfsoil_2d(nx,ny), &
 & Elatent_2d(nx,ny), dhw_2d(nx,ny), dhw0_2d(nx,ny), &
 & dhi_2d(nx,ny), dhi0_2d(nx,ny), dls0_2d(nx,ny), &
 & velfrict_2d(nx,ny),eflux0_kinem_2d(nx,ny), &
 & lamw_2d(nx,ny,M),Radbal_2d(nx,ny),hflux_2d(nx,ny))
 allocate ( snmelt_2d(nx,ny) )
 allocate (fl_sn_2d(nx,ny),fl_sn_init_2d(nx,ny), init(nx,ny), itop_2d(nx,ny))
 allocate (nstep_2d(nx,ny))
 allocate ( PEMF    (1:M)   , PDENS_DOWN (1:M+1) , &
 & PT_DOWN (1:M+1) , PSAL_DOWN  (1:M+1) , &
 & pt_down_f (1:M)  )
 allocate ( k_turb_T_flux(1:M), T_massflux(1:M) )
 allocate ( Tskin(1:2), Tskin_2d(nx,ny) )

 Tw1= 0. ; Tw2    = 0.
 Ti1= 0. ; Ti2    = 0.
 Tis1    = 0. ; Tis2   = 0.
 Sal1    = 0. ; Sal2   = 0.
 RS = 0. ; SR= 0.
 SRi= 0. ; SRdi   = 0.
 SR_botsnow = 0.d0
 lamw    = 0. ; lamsal = 0. ; lamsoil = 0. 
 Tsoil1  = 0. ; Tsoil2 = 0.
 rosoil  = 0. ; csoil  = 0. 
 dzs= 0. ; dzss   = 0. 
 wi1= 0. ; 
 wl1= 0. ; 
 S  = 0. ; Gen    = 0  ; F  = 0. ; TKE_turb_trans  = 0.
 KT = 0. ; k2= 0.
 u1 = 0. ; v1= 0. ; w1 = 0.
 u2 = 0. ; v2= 0. ; w2 = 0.
 E1 = 0. ; eps1   = 0.
 WLM0    = 0. ; WLM7   = 0. 
 bH = 0. ; POR    = 0.
 PSIMAX  = 0. ; FLWMAX = 0. ; DLMAX   = 0.
 KC = 0. ; KLengT = 0.
 dhw= 0. ; dhw0   = 0. 
 dhi= 0. ; dhi0   = 0.
 dls= 0. ; dls0   = 0.
 
 PEMF    = 0.d0 ; PDENS_DOWN = 0.d0
 PT_DOWN = 0.d0 ; PSAL_DOWN  = 0.d0
 pt_down_f = 0.d0
 zinv    = 0.d0
 
 k_turb_T_flux = 0.d0
 T_massflux    = 0.d0
  
 init    = 0
  
 !ddz = 1.d0/float(M)
 
! print*, ddz(1:M), sum(ddz(1:M))

! Specification of the dzeta-grid (liquid water layer)
 call GRID_CREATE

! Specification of the dzetai-grid (deep ice and upper ice layers)
 ddzi(:) = 1.d0/float(Mice)
 
 write(*,*) 'The dzeta-grid in water coloumn'
 do i = 1, M
   write(*,*) i, ddz(i)
 enddo
 !write(*,*) 'The sum is ', sum(ddz(1:M))
! STOP

 END SUBROUTINE ALLOCARRAYS


 SUBROUTINE GRID_CREATE
 use ARRAYS
 implicit none

 integer(4) :: i
 real(8) :: dzeta1,dzeta2

!The grid according to Burchard, 2002, pp. 97-98
 do i = 1, M
   dzeta1 = &
   & (dtanh(d_surf)+dtanh((d_surf+d_bot)*real(i-1)/real(M)-d_surf))/ &
   & (dtanh(d_surf)+dtanh(d_bot))
   dzeta2 = &
   & (dtanh(d_surf)+dtanh((d_surf+d_bot)*real(i)/real(M)-d_surf))/ &
   & (dtanh(d_surf)+dtanh(d_bot))
   ddz(i) = dzeta2 - dzeta1
 enddo

END SUBROUTINE GRID_CREATE


SUBROUTINE LININTERPOL (z1,f1,N1,z2,f2,N2,flag)

! Input variables
! Integers
integer(4), intent(in) :: N1
integer(4), intent(in) :: N2

! Reals
real(8), intent(in) :: z1(1:N1)
real(8), intent(in) :: z2(1:N2)

real(8), intent(in) :: f1(1:N1)

! Output variables
! Reals
real(8), intent(out) :: f2(1:N2)

logical, intent(out) :: flag

! Local variables
integer(4) :: i ! Loop index
integer(4) :: j ! Loop index

flag = .true.

do1 :do i = 1, N2
  do2: do j = 1, N1
    if (j == N1) exit do2
    if (z2(i)>=z1(j).and.z2(i)<z1(j+1)) exit do2
  enddo do2
  if (j == N1) then
!    flag = .false.
    f2(i) = f1(N1)
  else
    f2(i) = f1(j) + ( f1(j+1) - f1(j) )*( z2(i) - z1(j) )/( z1(j+1) - z1(j) )
  endif
enddo do1

END SUBROUTINE LININTERPOL
