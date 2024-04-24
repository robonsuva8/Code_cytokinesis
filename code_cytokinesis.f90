!**********************************************************************
! Program associated with to the manuscript:
! Cell intrinsic mechanical regulation of plasma membrane accumulation at the cytokinetic furrow, PNAS (2024)
! Roberto Alonso-Matilla, Alice Lam, and Teemu P. Miettinen
! Code developed by Roberto Alonso-Matilla Dec, 2023



!************************************************************************
!   MODULE CONSTANTS
!************************************************************************
module constants

double precision, parameter  :: pi = 3.1415926535897d0 ! Number pi

double precision, parameter  ::	R0      = 7.d0                      !Initial sphere radius
double precision, parameter  ::	L0      = pi*R0                     !Initial arc-length domain
double precision, parameter  :: Vcell0  = (4.d0/3.d0)*pi*R0**3.d0   !Cell volume
double precision, parameter  :: Acell0  = 4.d0*pi*R0**2.d0          !Cell surface area

!CORTEX
double precision, parameter  :: eta_c       = 50.d0  ! 2-D cortex viscosity
double precision, parameter  :: drag_mc     = 30.d0  ! Drag coefficient
double precision, parameter  :: bend_mod    = 0.d0   ! Bending modulus
double precision, parameter  :: sigma_a_max = 800.d0 ! Maximum cortical tension


!MYOSIN 
double precision, parameter  :: Nmyo        = 100.d4             !Total number of myosin molecules
double precision, parameter  :: konmyo      = 3636.d-2           !Myosin association rate
double precision, parameter  :: koffmyo     = 7.d-2              !Myosin dissociation rate
double precision, parameter  :: rho_myo_00  = Nmyo/(Acell0+6022.d-1*koffmyo*Vcell0/konmyo)
double precision, parameter  :: rho_sat_myo = 25.d-1*rho_myo_00  !Myosin saturation density 
double precision, parameter  :: diffmyo     = 1.d0               !Myosin diffusion coefficient


!ACTIN 
double precision, parameter  :: omega       = 2.d0               !Thickness anisotropic parameter
double precision, parameter  :: chi_ani     = 0.d0               !Strength anisotropic parameter
double precision, parameter  :: Nact        = 11447.d3           !Total number of actin units
double precision, parameter  :: konact      = 83864.d-2          !Actin association rate
double precision, parameter  :: koffact     = 3.d-1              !Actin dissociation rate
double precision, parameter  :: rho_act_00  = Nact/(Acell0+6022.d-1*koffact*Vcell0/konact)
double precision, parameter  :: rho_sat_act = 25.d-1*rho_act_00  !Actin saturation density 
double precision, parameter  :: diffact     = 1.d0               !Actin diffusion coefficient


!MEMBRANE
double precision, parameter  :: eta_m       = 3.d-3  !2-D membrane viscosity
double precision, parameter  :: E_m         = 40     !Effective area expansion modulus of the membrane
double precision, parameter  :: drag_m      = 100    !Drag force between membrane lipids and proteins 
double precision, parameter  :: rho_mem0    = 5.d6   !Membrane lipid density 

!MEMBRANE-CORTEX LINKERS
double precision, parameter  :: Nlink       = 0.d0   !Total number of linkers
double precision, parameter  :: konlink     = 0.d0   !Linkers association rate
double precision, parameter  :: kofflink    = 0.d0   !Linkers dissociation rate
double precision, parameter  :: difflink    = 0.d0   !Linkers diffusion coefficient


integer, parameter           :: Nu = 140         ! Number of equally spaced collocation points
double precision, parameter  :: du = dble(L0/Nu) ! Distance between neighboring collocation points

end module 


!************************************************************************
!   MAIN PROGRAM
!************************************************************************

program tubular
 use constants
 implicit none

 double precision                            :: dt,rho_myo0,rho_act0,rho_free_myo,rho_free_act,dd,pressure
 double precision                            :: N_myo_cortex,N_act_cortex,ddm,s_eq
 integer                                     :: cycles,i,j
 double precision,dimension(:),allocatable   :: psi,coordinatex,rho_myo,rho_act,r,alpha,h,css,ctt,gamma,f
 double precision,dimension(:),allocatable   :: du_f,du_h,du_css,duu_css,a1,a2,a3,a4,a5,b1,b2,b3,RHS,vs,vn
 double precision,dimension(:),allocatable   :: ss,z,du_rhomyo,duu_rhomyo,rhomyo_p,du_rhoact,duu_rhoact,rhoact_p
 double precision,dimension(:),allocatable   :: r_p,h_p,css_p,ctt_p,gamma_p,du_vs,du_vn,duu_vn,integrand
 double precision,dimension(:),allocatable   :: rho_mem,a1m,a2m,a3m,a4m,RHSm,vmem,rhomem_p,du_vmem,du_rhomem
 double precision,dimension(:),allocatable   :: surf_a,surf_a_p,surf_a_css,surf_a_p_css,surf_a_ctt,surf_a_p_ctt
 double precision,dimension(:),allocatable   :: chi_vector,f_act_ss,f_act_tt,du_f_act_ss,du_f_act_tt
 double precision,dimension(:,:),allocatable :: matrix,matrixm
 integer,dimension(:),allocatable            :: indice,indicem
 integer*8                                   :: itime,icycles,ifile
 character*32 filename1

 allocate(psi(Nu),coordinatex(Nu),rho_myo(Nu),rho_act(Nu),r(Nu),alpha(Nu),h(Nu),css(Nu),ctt(Nu),gamma(Nu),f(Nu))
 allocate(du_f(Nu),du_h(Nu),du_css(Nu),duu_css(Nu),a1(Nu),a2(Nu),a3(Nu),a4(Nu),a5(Nu),b1(Nu),b2(Nu),b3(Nu))
 allocate(RHS(2*Nu+1),matrix(2*Nu+1,2*Nu+1),indice(2*Nu+1),vs(Nu),vn(Nu),ss(Nu),z(Nu))
 allocate(du_rhomyo(Nu),duu_rhomyo(Nu),rhomyo_p(Nu),du_rhoact(Nu),duu_rhoact(Nu),rhoact_p(Nu))
 allocate(r_p(Nu),h_p(Nu),css_p(Nu),ctt_p(Nu),gamma_p(Nu),du_vs(Nu),du_vn(Nu),duu_vn(Nu),integrand(Nu))
 allocate(rho_mem(Nu),vmem(Nu),rhomem_p(Nu),du_vmem(Nu),du_rhomem(Nu))
 allocate(a1m(Nu),a2m(Nu),a3m(Nu),a4m(Nu),matrixm(Nu,Nu),indicem(Nu),RHSm(Nu),surf_a(Nu),surf_a_p(Nu))
 allocate(surf_a_css(Nu),surf_a_p_css(Nu),surf_a_ctt(Nu),surf_a_p_ctt(Nu))
 allocate(chi_vector(Nu),f_act_ss(Nu),f_act_tt(Nu),du_f_act_ss(Nu),du_f_act_tt(Nu))

! Time step
dt     = 1.d-6
cycles = 2

! Coordinates of collocation points (polar angle in spherical coordinates)
do i=1,Nu
 psi(i) = (i-5.d-1)*(pi/Nu)-(pi/2.d0)
 coordinatex(i) = R0*sin(psi(i))
enddo

! Initial values: Myosin density field rho_myo(u,t=0), actin density field rho_act(u,t=0), 
! Membrane lipid density field rho_mem(u,t=0), radial coordinate r(u,t=0)
! tangent angle alpha(u,t=0), parameter transformation h(u,t=0), curvature css(u,t=0), curvature ctt(u,t=0)
! christoffel symbol gamma(u,t=0)
rho_myo0     = Nmyo/(Acell0+6022.d-1*koffmyo*Vcell0/konmyo)
rho_free_myo = koffmyo*rho_myo0/konmyo
rho_act0     = Nact/(Acell0+6022.d-1*koffact*Vcell0/konact)
rho_free_act = koffact*rho_act0/konact

do i=1,Nu
 rho_myo(i) = rho_myo0+1.d-2*exp(-(coordinatex(i)/3.d-1)**2.d0)
 rho_act(i) = rho_act0+1.d-2*exp(-(coordinatex(i)/3.d-1)**2.d0)
 rho_mem(i) = rho_mem0
 r(i)       = R0*cos(psi(i))
 alpha(i)   = acos(-sin(psi(i)))
 h(i)       = 1.d0
 css(i)     = 1.d0/R0
 ctt(i)     = 1.d0/R0
 gamma(i)   = -tan(psi(i))/R0
enddo

! Obtain the arc-length 
ss(1) = h(1)*du/2.d0
do i=2,Nu
  ss(i) = ss(i-1)+(h(i-1)+h(i))*du/2.d0
enddo
s_eq = (ss(70)+ss(71))/2.d0

! Saturation function for the active stress
do i=1,Nu
  chi_vector(i) = chi_ani*exp(-5.d-1*((ss(i)-s_eq)/omega)**2.d0)
  f(i) = sigma_a_max*erf(rho_myo(i)/rho_sat_myo)*erf(rho_act(i)/rho_sat_act)
  f_act_ss(i) = f(i)
  f_act_tt(i) = f(i)*(1+chi_vector(i))
enddo

! Initialize membrane lipid velocity
vmem = 0.d0

! Derivatives of f(u,t),h(u,t),css(u,t)
call first_deriv_type1(f_act_ss,du_f_act_ss)
call first_deriv_type1(f_act_tt,du_f_act_tt)
call first_deriv_type1(h,du_h)
call first_deriv_type1(css,du_css)
call second_deriv_type1(css,duu_css)

! Parameters used to construct the hydrodynamic matrix
do i=1,Nu
 a1(i) = 2.d0*eta_c/(h(i)**2.d0)
 a2(i) = (2.d0*eta_c/h(i))*(gamma(i)-(du_h(i)/(h(i)**2.d0)))
 a3(i) = -(2.d0*eta_c*gamma(i)**2.d0)-drag_mc
 a4(i) = 2.d0*eta_c*(du_css(i)/h(i)+gamma(i)*css(i)-ctt(i)*gamma(i))
 a5(i) = 2.d0*eta_c*css(i)/h(i)
 b1(i) = 2.d0*eta_c*css(i)/h(i)
 b2(i) = 2.d0*eta_c*ctt(i)*gamma(i)
 b3(i) = 2.d0*eta_c*(css(i)**2.d0+ctt(i)**2.d0) 
enddo

RHS    = 0.d0
matrix = 0.d0
! Coefficients of the hydrodynamic matrix. Linear problem solving cortex tangential and normal force balance equations
do i=1,Nu
  if(i.eq.1) then 
   matrix(i,i)=-(3.d0*a1(i)/du**2.d0)+(a2(i)/(2.d0*du))+a3(i)
   matrix(i,i+1)=(a1(i)/du**2.d0)+(a2(i)/(2.d0*du))
   matrix(i,Nu+i)=-(a5(i)/(2.d0*du))+a4(i)
   matrix(i,Nu+i+1)=(a5(i)/(2.d0*du))
   RHS(i)=-(1.d0/h(i))*du_f_act_ss(i)-gamma(i)*(f_act_ss(i)-f_act_tt(i))-drag_mc*vmem(i) 
   matrix(Nu+i,i)=(b1(i)/(2.d0*du))+b2(i)   
   matrix(Nu+i,i+1)=b1(i)/(2.d0*du)   
   matrix(Nu+i,Nu+i)=b3(i)   
   RHS(Nu+i)=-f_act_ss(i)*css(i)-f_act_tt(i)*ctt(i)+bend_mod*(ctt(i)+css(i))*((ctt(i)-css(i))**2.d0)+&
   2.d0*bend_mod*((duu_css(i)/(h(i)**2.d0))-&
   (du_css(i)*du_h(i)/(h(i)**3.d0))+(2.d0*gamma(i)*du_css(i)/h(i))-(ctt(i)*css(i)+gamma(i)**2.d0)*(css(i)-ctt(i)))
  elseif(i.eq.Nu) then
   matrix(i,i-1)=(a1(i)/du**2.d0)-(a2(i)/(2.d0*du))     
   matrix(i,i)=-(3.d0*a1(i)/du**2.d0)-(a2(i)/(2.d0*du))+a3(i)
   matrix(i,Nu+i-1)=-(a5(i)/(2.d0*du))
   matrix(i,Nu+i)=(a5(i)/(2.d0*du))+a4(i)
   RHS(i)=-(1.d0/h(i))*du_f_act_ss(i)-gamma(i)*(f_act_ss(i)-f_act_tt(i))-drag_mc*vmem(i) 
   matrix(Nu+i,i-1)=-b1(i)/(2.d0*du)   
   matrix(Nu+i,i)=-(b1(i)/(2.d0*du))+b2(i)  
   matrix(Nu+i,Nu+i)=b3(i)
   RHS(Nu+i)=-f_act_ss(i)*css(i)-f_act_tt(i)*ctt(i)+bend_mod*(ctt(i)+css(i))*((ctt(i)-css(i))**2.d0)+&
   2.d0*bend_mod*((duu_css(i)/(h(i)**2.d0))-&
   (du_css(i)*du_h(i)/(h(i)**3.d0))+(2.d0*gamma(i)*du_css(i)/h(i))-(ctt(i)*css(i)+gamma(i)**2.d0)*(css(i)-ctt(i)))
  else     
   matrix(i,i-1)=(a1(i)/du**2.d0)-(a2(i)/(2.d0*du))      
   matrix(i,i)=-(2.d0*a1(i)/du**2.d0)+a3(i)
   matrix(i,i+1)=(a1(i)/du**2.d0)+(a2(i)/(2.d0*du))
   matrix(i,Nu+i-1)=-(a5(i)/(2.d0*du))
   matrix(i,Nu+i)=a4(i)
   matrix(i,Nu+i+1)=(a5(i)/(2.d0*du))
   RHS(i)=-(1.d0/h(i))*du_f_act_ss(i)-gamma(i)*(f_act_ss(i)-f_act_tt(i))-drag_mc*vmem(i) 
   matrix(Nu+i,i-1)=-b1(i)/(2.d0*du)   
   matrix(Nu+i,i)=b2(i)  
   matrix(Nu+i,i+1)=b1(i)/(2.d0*du)   
   matrix(Nu+i,Nu+i)=b3(i) 
   RHS(Nu+i)=-f_act_ss(i)*css(i)-f_act_tt(i)*ctt(i)+bend_mod*(ctt(i)+css(i))*((ctt(i)-css(i))**2.d0)+&
   2.d0*bend_mod*((duu_css(i)/(h(i)**2.d0))-&
   (du_css(i)*du_h(i)/(h(i)**3.d0))+(2.d0*gamma(i)*du_css(i)/h(i))-(ctt(i)*css(i)+gamma(i)**2.d0)*(css(i)-ctt(i)))
  endif
enddo

RHS(2*Nu+1)=0.d0
do i=1,Nu
 matrix(Nu+i,2*Nu+1)=-1.d0  
 matrix(2*Nu+1,Nu+i)=r(i)*h(i)*du   
enddo

! Solution of the system of equations. Matrix is destroyed. Solution is stored in RHS vector
 call ludcmp(matrix,2*Nu+1,2*Nu+1,indice,dd)
 call lubksb(matrix,2*Nu+1,2*Nu+1,indice,RHS)

! We store the values of vs, vn and p.
 do i=1,Nu
   vs(i) = RHS(i)
   vn(i) = RHS(Nu+i)
 enddo
 pressure = RHS(2*Nu+1)

! Obtain the arc-length 
ss(1) = h(1)*du/2.d0
do i=2,Nu
  ss(i) = ss(i-1)+(h(i-1)+h(i))*du/2.d0
enddo
s_eq = (ss(70)+ss(71))/2.d0

! Obtain z coordinate
z(1) = r(1)*ctt(1)*h(1)*du/2.d0
do i=2,Nu
  z(i) = z(i-1)+(r(i-1)*ctt(i-1)*h(i-1)+r(i)*ctt(i)*h(i))*du/2.d0
enddo  

! Initialize surface area factor
surf_a     = 1.d0
surf_a_css = 1.d0
surf_a_ctt = 1.d0

! Time marching
do icycles = 1,cycles
 do itime = 1,200000000

  ! Some derivatives
  call first_deriv_type2(vs,du_vs)
  call first_deriv_type1(vn,du_vn)
  call second_deriv_type1(vn,duu_vn)
  call first_deriv_type1(h,du_h)
  call first_deriv_type1(rho_myo,du_rhomyo)
  call first_deriv_type1(rho_act,du_rhoact)
  call first_deriv_type1(css,du_css)
  call second_deriv_type1(rho_myo,duu_rhomyo)
  call second_deriv_type1(rho_act,duu_rhoact)

  call first_deriv_type1(rho_mem,du_rhomem)

  ! Solve tangential force balance on membrane 
  do i=1,Nu
    a1m(i) = 2.d0*eta_m/(h(i)**2.d0)
    a2m(i) = (2.d0*eta_m/h(i))*(gamma(i)-(du_h(i)/(h(i)**2.d0)))
    a3m(i) = -(2.d0*eta_m*gamma(i)**2.d0)-drag_mc-drag_m
    a4m(i) = -(2.d0*eta_m/h(i))*(du_css(i)*vn(i)+css(i)*du_vn(i))-2.d0*eta_m*gamma(i)*vn(i)*(css(i)-ctt(i))+&
    (E_m/(rho_mem0*h(i)))*du_rhomem(i)-drag_mc*vs(i)
  enddo
  matrixm = 0.d0
  RHSm    = 0.d0
  do i = 1,Nu
    if(i.eq.1) then 
      matrixm(i,i)   = -(3.d0*a1m(i)/du**2.d0)+(a2m(i)/(2.d0*du))+a3m(i)
      matrixm(i,i+1) = (a1m(i)/du**2.d0)+(a2m(i)/(2.d0*du))
      RHSm(i)        = a4m(i)
    elseif(i.eq.Nu) then
      matrixm(i,i-1) = (a1m(i)/du**2.d0)-(a2m(i)/(2.d0*du))    
      matrixm(i,i)   = -(3.d0*a1m(i)/du**2.d0)-(a2m(i)/(2.d0*du))+a3m(i)
      RHSm(i)        = a4m(i)
    else     
      matrixm(i,i-1) = (a1m(i)/du**2.d0)-(a2m(i)/(2.d0*du))    
      matrixm(i,i)   = -(2.d0*a1m(i)/du**2.d0)+a3m(i)
      matrixm(i,i+1) = (a1m(i)/du**2.d0)+(a2m(i)/(2.d0*du))
      RHSm(i)        = a4m(i)
    endif
  enddo
  call ludcmp(matrixm,Nu,Nu,indicem,ddm)
  call lubksb(matrixm,Nu,Nu,indicem,RHSm)

  do i=1,Nu
    vmem(i) = RHSm(i)
  enddo
  call first_deriv_type2(vmem,du_vmem)

  ! Obtain values for the next time step
  do i=1,Nu
    r_p(i)      = r(i)+dt*vn(i)*r(i)*ctt(i)
    h_p(i)      = h(i)+dt*h(i)*css(i)*vn(i)
    css_p(i)    = css(i)+dt*(-(css(i)**2.d0)*vn(i)-(duu_vn(i)/(h(i)**2.d0))+(du_h(i)*du_vn(i)/h(i)**3.d0))
    ctt_p(i)    = ctt(i)+dt*(-(ctt(i)**2.d0)*vn(i)-gamma(i)*du_vn(i)/h(i))
    gamma_p(i)  = gamma(i)+dt*(ctt(i)*du_vn(i)/h(i)-ctt(i)*gamma(i)*vn(i))
 
    rhomyo_p(i) = rho_myo(i)+dt*(-(rho_myo(i)*du_vs(i)/h(i))-(vs(i)*du_rhomyo(i)/h(i))-gamma(i)*rho_myo(i)*vs(i)-&
    rho_myo(i)*vn(i)*(css(i)+ctt(i))+diffmyo*((duu_rhomyo(i)/h(i)**2.d0)-&
    du_h(i)*du_rhomyo(i)/h(i)**3.d0+gamma(i)*du_rhomyo(i)/h(i))+konmyo*rho_free_myo-koffmyo*rho_myo(i))

    rhoact_p(i) = rho_act(i)+dt*(-(rho_act(i)*du_vs(i)/h(i))-(vs(i)*du_rhoact(i)/h(i))-gamma(i)*rho_act(i)*vs(i)-&
    rho_act(i)*vn(i)*(css(i)+ctt(i))+diffact*((duu_rhoact(i)/h(i)**2.d0)-&
    du_h(i)*du_rhoact(i)/h(i)**3.d0+gamma(i)*du_rhoact(i)/h(i))+konact*rho_free_act-koffact*rho_act(i))

    rhomem_p(i) = rho_mem(i)+dt*(-(rho_mem(i)*du_vmem(i)/h(i))-(vmem(i)*du_rhomem(i)/h(i))-&
    gamma(i)*rho_mem(i)*vmem(i)-rho_mem(i)*vn(i)*(css(i)+ctt(i)))

    surf_a_p(i)     = surf_a(i)*(1.d0+dt*vn(i)*(css(i)+ctt(i)))
    surf_a_p_css(i) = surf_a_css(i)*(1.d0+dt*vn(i)*(css(i)))
    surf_a_p_ctt(i) = surf_a_ctt(i)*(1.d0+dt*vn(i)*(ctt(i)))

  enddo

  ! Update values
  r          = r_p
  h          = h_p
  css        = css_p
  ctt        = ctt_p
  gamma      = gamma_p 
  rho_myo    = rhomyo_p
  rho_act    = rhoact_p
  rho_mem    = rhomem_p
  surf_a     = surf_a_p
  surf_a_css = surf_a_p_css
  surf_a_ctt = surf_a_p_ctt

  do i=1,Nu
    if (rho_myo(i) .lt. 0.d0) then
      rho_myo(i) = 0.d0
    endif
    if (rho_act(i) .lt. 0.d0) then
      rho_act(i) = 0.d0
    endif
    f(i) = sigma_a_max*erf(rho_myo(i)/rho_sat_myo)*erf(rho_act(i)/rho_sat_act)
    chi_vector(i) = chi_ani*exp(-5.d-1*((ss(i)-s_eq)/omega)**2.d0)
    f_act_ss(i) = f(i)
    f_act_tt(i) = f(i)*(1+chi_vector(i))
  enddo

  !We compute cortex amounts 
  do i=1,Nu
    integrand(i) = rho_myo(i)*r(i)*h(i)
  enddo
  N_myo_cortex = 2*pi*du*(5.d-1*(integrand(1)+integrand(Nu))+sum(integrand(2:Nu-1)))
 
  do i=1,Nu
    integrand(i) = rho_act(i)*r(i)*h(i)
  enddo
  N_act_cortex = 2*pi*du*(5.d-1*(integrand(1)+integrand(Nu))+sum(integrand(2:Nu-1)))

  rho_free_myo = (Nmyo-N_myo_cortex)/(Vcell0*602.2)
  rho_free_act = (Nact-N_act_cortex)/(Vcell0*602.2)

  ! Derivatives of f(u,t),h(u,t),css(u,t)
  call first_deriv_type1(f_act_ss,du_f_act_ss)
  call first_deriv_type1(f_act_tt,du_f_act_tt)
  call first_deriv_type1(h,du_h)
  call first_deriv_type1(css,du_css)
  call second_deriv_type1(css,duu_css)

  ! Parameters used to construct the hydrodynamic matrix
  do i=1,Nu
   a1(i) = 2.d0*eta_c/(h(i)**2.d0)
   a2(i) = (2.d0*eta_c/h(i))*(gamma(i)-(du_h(i)/(h(i)**2.d0)))
   a3(i) = -(2.d0*eta_c*gamma(i)**2.d0)-drag_mc
   a4(i) = 2.d0*eta_c*(du_css(i)/h(i)+gamma(i)*css(i)-ctt(i)*gamma(i))
   a5(i) = 2.d0*eta_c*css(i)/h(i)
   b1(i) = 2.d0*eta_c*css(i)/h(i)
   b2(i) = 2.d0*eta_c*ctt(i)*gamma(i)
   b3(i) = 2.d0*eta_c*(css(i)**2.d0+ctt(i)**2.d0) 
  enddo

  RHS    = 0.d0
  matrix = 0.d0

  ! Coefficients of the hydrodynamic matrix. Linear problem solving tangential and normal force balance equations
  do i=1,Nu
    if(i.eq.1) then 
     matrix(i,i)=-(3.d0*a1(i)/du**2.d0)+(a2(i)/(2.d0*du))+a3(i)
     matrix(i,i+1)=(a1(i)/du**2.d0)+(a2(i)/(2.d0*du))
     matrix(i,Nu+i)=-(a5(i)/(2.d0*du))+a4(i)
     matrix(i,Nu+i+1)=(a5(i)/(2.d0*du))
     RHS(i)=-(1.d0/h(i))*du_f_act_ss(i)-gamma(i)*(f_act_ss(i)-f_act_tt(i))-drag_mc*vmem(i) 
     matrix(Nu+i,i)=(b1(i)/(2.d0*du))+b2(i)   
     matrix(Nu+i,i+1)=b1(i)/(2.d0*du)   
     matrix(Nu+i,Nu+i)=b3(i)   
     RHS(Nu+i)=-f_act_ss(i)*css(i)-f_act_tt(i)*ctt(i)+bend_mod*(ctt(i)+css(i))*((ctt(i)-css(i))**2.d0)+&
     2.d0*bend_mod*((duu_css(i)/(h(i)**2.d0))-&
     (du_css(i)*du_h(i)/(h(i)**3.d0))+(2.d0*gamma(i)*du_css(i)/h(i))-(ctt(i)*css(i)+gamma(i)**2.d0)*(css(i)-ctt(i)))
    elseif(i.eq.Nu) then
     matrix(i,i-1)=(a1(i)/du**2.d0)-(a2(i)/(2.d0*du))     
     matrix(i,i)=-(3.d0*a1(i)/du**2.d0)-(a2(i)/(2.d0*du))+a3(i)
     matrix(i,Nu+i-1)=-(a5(i)/(2.d0*du))
     matrix(i,Nu+i)=(a5(i)/(2.d0*du))+a4(i)
     RHS(i)=-(1.d0/h(i))*du_f_act_ss(i)-gamma(i)*(f_act_ss(i)-f_act_tt(i))-drag_mc*vmem(i) 
     matrix(Nu+i,i-1)=-b1(i)/(2.d0*du)   
     matrix(Nu+i,i)=-(b1(i)/(2.d0*du))+b2(i)  
     matrix(Nu+i,Nu+i)=b3(i)
     RHS(Nu+i)=-f_act_ss(i)*css(i)-f_act_tt(i)*ctt(i)+bend_mod*(ctt(i)+css(i))*((ctt(i)-css(i))**2.d0)+&
     2.d0*bend_mod*((duu_css(i)/(h(i)**2.d0))-&
     (du_css(i)*du_h(i)/(h(i)**3.d0))+(2.d0*gamma(i)*du_css(i)/h(i))-(ctt(i)*css(i)+gamma(i)**2.d0)*(css(i)-ctt(i)))
    else     
     matrix(i,i-1)=(a1(i)/du**2.d0)-(a2(i)/(2.d0*du))      
     matrix(i,i)=-(2.d0*a1(i)/du**2.d0)+a3(i)
     matrix(i,i+1)=(a1(i)/du**2.d0)+(a2(i)/(2.d0*du))
     matrix(i,Nu+i-1)=-(a5(i)/(2.d0*du))
     matrix(i,Nu+i)=a4(i)
     matrix(i,Nu+i+1)=(a5(i)/(2.d0*du))
     RHS(i)=-(1.d0/h(i))*du_f_act_ss(i)-gamma(i)*(f_act_ss(i)-f_act_tt(i))-drag_mc*vmem(i) 
     matrix(Nu+i,i-1)=-b1(i)/(2.d0*du)   
     matrix(Nu+i,i)=b2(i)  
     matrix(Nu+i,i+1)=b1(i)/(2.d0*du)   
     matrix(Nu+i,Nu+i)=b3(i) 
     RHS(Nu+i)=-f_act_ss(i)*css(i)-f_act_tt(i)*ctt(i)+bend_mod*(ctt(i)+css(i))*((ctt(i)-css(i))**2.d0)+&
     2.d0*bend_mod*((duu_css(i)/(h(i)**2.d0))-&
     (du_css(i)*du_h(i)/(h(i)**3.d0))+(2.d0*gamma(i)*du_css(i)/h(i))-(ctt(i)*css(i)+gamma(i)**2.d0)*(css(i)-ctt(i)))
    endif
  enddo

  RHS(2*Nu+1) = 0.d0
  do i=1,Nu
    matrix(Nu+i,2*Nu+1) = -1.d0  
    matrix(2*Nu+1,Nu+i) = r(i)*h(i)*du   
  enddo

  ! Solution of the system of equations. Matrix is destroyed. Solution is stored in RHS vector
  call ludcmp(matrix,2*Nu+1,2*Nu+1,indice,dd)
  call lubksb(matrix,2*Nu+1,2*Nu+1,indice,RHS)

  ! We store the values of vs, vn and p.
  do i=1,Nu
    vs(i) = RHS(i)
    vn(i) = RHS(Nu+i)
  enddo
  pressure = RHS(2*Nu+1)

  ! Obtain the arc-length 
  ss(1) = h(1)*du/2.d0
  do i=2,Nu
    ss(i) = ss(i-1)+(h(i-1)+h(i))*du/2.d0
  enddo
  s_eq = (ss(70)+ss(71))/2.d0  


  ! Obtain z coordinate
  z(1) = r(1)*ctt(1)*h(1)*du/2.d0
  do i=2,Nu
    z(i) = z(i-1)+(r(i-1)*ctt(i-1)*h(i-1)+r(i)*ctt(i)*h(i))*du/2.d0
  enddo   


  ! Save variables from time to time
  if (2000.0*(itime/2000).eq.1.0*itime) then
    ifile = itime/2000
    write(filename1,'(A4,I6,A3,I2,A4)') 'data',ifile,'cyc',icycles,'.txt'
    open(unit=5,file=filename1,status='unknown')
    do i=1,Nu
      write(5,*) vs(i),vn(i),r(i),ss(i),rho_myo(i),rho_act(i),z(i),rho_mem(i),vmem(i),surf_a(i),surf_a_css(i),surf_a_ctt(i)
    enddo
    close(5)
  endif
    
 enddo
enddo

end

!-------------------------------------------------------------------------------------------------------------------------------------------------
! Subroutines that solve the hydrodynamic linear system of equations. 
! This routine is used in combination with lubksb to solve linear equations or invert a matrix. We solve the linear set of equations A · x = b
! call ludcmp(a,n,np,indx,d)
! call lubksb(a,n,np,indx,b)
! The answer x will be returned in b. Your original matrix A will be destroyed.
!-------------------------------------------------------------------------------------------------------------------------------------------------
subroutine ludcmp(a,n,np,indx,d)

 integer n,np,indx(n),NMAX
 double precision d,a(np,np),TINY
 PARAMETER (NMAX=40000,TINY=1.0e-20) ! Largest expected n, and a small number.
 integer i,imax,j,k
 double precision aamax,dum,sum,vv(NMAX)

 d=1.d0

 do i=1,n
  aamax=0.d0
  do j=1,n
   if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
  enddo
  if (aamax.eq.0.)  then
!   write (*,*) 333	
  endif		
  vv(i)=1.d0/aamax
 enddo
 do j=1,n
  do i=1,j-1
   sum=a(i,j)
   do k=1,i-1
    sum=sum-a(i,k)*a(k,j)
   enddo
   a(i,j)=sum
  enddo
  aamax=0.d0
  do i=j,n
   sum=a(i,j)
   do k=1,j-1
    sum=sum-a(i,k)*a(k,j)
   enddo
   a(i,j)=sum
   dum=vv(i)*abs(sum)
   if (dum.ge.aamax) then
    imax=i
    aamax=dum
   endif
  enddo
  if (j.ne.imax)then
   do k=1,n
    dum=a(imax,k)
    a(imax,k)=a(j,k)
    a(j,k)=dum
   enddo
   d=-d
   vv(imax)=vv(j)
  endif
  indx(j)=imax
  if(a(j,j).eq.0.) a(j,j)=TINY
  if(j.ne.n)then
   dum=1./a(j,j)
   do i=j+1,n
    a(i,j)=a(i,j)*dum
   enddo
  endif
 enddo

 return

end subroutine ludcmp

subroutine lubksb(a,n,np,indx,b)

 integer n,np,indx(n)
 double precision a(np,np),b(n)
 integer i,ii,j,ll
 double precision sum

 ii=0

 do i=1,n
  ll=indx(i)
  sum=b(ll)
  b(ll)=b(i)
  if (ii.ne.0)then
   do j=ii,i-1
    sum=sum-a(i,j)*b(j)
   enddo 
  else if (sum.ne.0.) then
   ii=i
  endif
  b(i)=sum
 enddo
 do i=n,1,-1
  sum=b(i)
  do j=i+1,n
   sum=sum-a(i,j)*b(j)
  enddo
  b(i)=sum/a(i,i)
 enddo

 return

end subroutine lubksb

!-------------------------------------------------------------------------------------------------------------------------------------------------
! It calculates the first derivative (type 1) of vect and store the results in dx_vect. 2nd order accurate
!-------------------------------------------------------------------------------------------------------------------------------------------------
subroutine first_deriv_type1(vect,dx_vect)
 use constants
 implicit none

 integer i
 double precision vect(Nu), dx_vect(Nu) 

 do i=1,Nu
  if(i.eq.1) then
   dx_vect(i)=(vect(i+1)-vect(i))/(2.d0*du)
  elseif(i.eq.Nu) then
   dx_vect(i)=(vect(i)-vect(i-1))/(2.d0*du)
  else
   dx_vect(i)=(vect(i+1)-vect(i-1))/(2.d0*du)
  endif
 enddo

end subroutine first_deriv_type1

!-------------------------------------------------------------------------------------------------------------------------------------------------
! It calculates the first derivative (type 2) of vect and store the results in dx_vect. 2nd order accurate
!-------------------------------------------------------------------------------------------------------------------------------------------------
subroutine first_deriv_type2(vect,dx_vect)
 use constants
 implicit none

 integer i
 double precision vect(Nu), dx_vect(Nu) 

 do i=1,Nu
  if(i.eq.1) then
   dx_vect(i)=(vect(i+1)+vect(i))/(2*du)
  elseif(i.eq.Nu) then
   dx_vect(i)=(-vect(i)-vect(i-1))/(2*du)
  else
   dx_vect(i)=(vect(i+1)-vect(i-1))/(2*du)
  endif
 enddo

end subroutine first_deriv_type2

!-------------------------------------------------------------------------------------------------------------------------------------------------
! It calculates the second derivative of vect and store the results in dxx_vect. 2nd order accurate
!-------------------------------------------------------------------------------------------------------------------------------------------------
subroutine second_deriv_type1(vect,dxx_vect)
 use constants
 implicit none

 integer i
 double precision vect(Nu), dxx_vect(Nu) 

 do i=1,Nu
  if(i.eq.1) then
   dxx_vect(i)=(vect(i+1)-vect(i))/(du**2.d0)
  elseif(i.eq.Nu) then
   dxx_vect(i)=(-vect(i)+vect(i-1))/(du**2.d0)
  else
   dxx_vect(i)=(vect(i+1)-2*vect(i)+vect(i-1))/(du**2.d0)
  endif
 enddo

end subroutine second_deriv_type1







