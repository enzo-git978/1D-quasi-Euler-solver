module Module_solve

 ! modules utilisés
  use Module_sortie
  use Module_preprocess
  use Module_riemann

contains

! ----------------------------------------------------------------------
! METHODES
! ----------------------------------------------------------------------

subroutine Solve_Eqn_Euler_MUSCL(dx,rho,u,P,E)
! ======================================================================================= 
! Subroutine pour resoudre le systeme d'Euler fluide parfait 1D compressible
! avec MUSCL
! ======================================================================================= 
! Parametres globaux 
use Module_parametres, only : N, maxiter, tol_rho, tol_u, tol_e

! Déclarations
  implicit none
  double precision, intent(in) 	 			   	:: dx								    	! Double from the main
  double precision, dimension(:,:),allocatable 	:: w,w_s,w_L,w_R,St							! Vecteur w au temps t
  double precision, dimension(:,:),allocatable  :: w_mid									! Vecteur w au temps t{n+1/2}
  double precision, dimension(:,:),allocatable  :: w_new									! Vecteur w au temps t{n+1}
  double precision								:: rho_s, u_s, E_s,dt,t,res_rho,res_u,res_e	! Double 	
  integer 										:: iter										! Integer for loop
  double precision, dimension(:),allocatable	:: x_i,An,Anf	    				    	! Maillage et Geometrie de la tuyere
  double precision, dimension(:,:),allocatable 	:: F_L,F_R									! Vecteur Flux
  double precision, dimension(:,:),allocatable 	:: F_num									! Vecteur Flux numériques 
  double precision, dimension(:),allocatable 	:: Ma,s,Temp								! Mach number, entropy and Temperature
  double precision, dimension(:),intent(inout) 	:: rho, u, P, E		    					! Caractéristiques physique du fluide
  
  ! Fin des declarations
 !=======================================================================================
	! Allocations dynamiques des variables locales de la subroutine
	allocate(w(3,-1:N+2),w_new(3,-1:N+2))    ! intern cells : 1, N ; ghost cells : -1, 0 and N+1, N+2
	allocate(w_s(3,-1:N+2))					 ! Vecteur d'etat initial (sert pour subroutine BC)
	allocate(w_L(3,N+1),w_R(3,N+1))          ! Interpolation des etats gauches et droites du vecteur conservatif
	allocate(F_L(3,N+1), F_R(3,N+1))	     ! Flux 
	allocate(F_num(3,N+1))				     ! Flux numeriques calculés avec Roes scheme 
	allocate(St(3,N))				         ! Source term
	allocate(x_i(-1:N+2),An(-1:N+2),Anf(0:N))! Maillage et geometrie de tuyere	
	allocate(Ma(N), s(N), Temp(N))		     ! Variables thermophysiques Mach number, entropy, temperature
	
	t = 0.0D0   ! temps de la simulation 
	res_rho = 1.0d9; res_u = 1.0d9; res_e = 1.0d9 ! Residus
	
	! Mesh creation
	call Mesh(dx,x_i)
	
	! Nozzle geometry creation
	call Nozzle_geom(dx,x_i,An,Anf)
	
	! Calcul des Conditions initiales 
	call Condition_initiale(x_i,An,w)
	
	! Creer le vecteur init w_s (sert pour bc)
	w_s = w
	
	! Boucle temporelle
	do iter=1,maxiter
		call Compute_time_step(dx,rho,u,E,P,An,w,dt) ! Calcul du pas de temps avec CFL
		t = t + dt
		!print*,'t = ', t, ' 	| dt =', dt
		
		! Interpolation (ordre 1 ou reconstruction MUSCL ordre 2 ou plus)
		call Interp(rho,u,P,E,An,w,w_L,w_R)
		
		! Calcul des flux a gauche pour chaque cellule
		call Calcul_Flux(w_L,F_L) 
		
		! Calcul des flux a droite pour chaque cellule
		call Calcul_Flux(w_R,F_R) 
		
		! Calcul des Flux numérique aux interfaces avec un schéma numerique 
		call Num_flux(w_L,w_R,F_L,F_R,F_num)
		
		! Calcul du source term
		call Source_term(dx,w,An,Anf,St)
		
		! Calcul des variables au temps t{n+1}
		call Update_variables(dx,dt,w,F_num,St,rho_s,u_s,E_s,rho,u,P,E,An,Anf,w_new)
		
		! Appliquer les conditions aux bords
		call Boundary_conditions(w_s,An,w_new)
		
		! Calcul des residus
		call Compute_residus(An,w,w_new,res_rho,res_u,res_e)
		
		! Test de Convergence
		if (res_rho <= tol_rho .and. res_u <= tol_u .and. res_e <= tol_e) then
			write(*, '(A,I0,A)') "Convergence atteinte en ", iter, " iterations !"
			write(*, '(A,1PE10.4,A,1PE10.4,A,1PE10.4)') &
				"Residus: rho = ", res_rho, " | u = ", res_u, " | E = ", res_e
			exit
		end if
		
		! Affichage toutes les 100 itérations pour monitorer
		if (mod(iter, 100) == 0) then
			write(*, '(A,I6,A,1PE10.4,A,1PE10.4,A,1PE10.4,A,1PE10.4)') &
				"| It:", iter, &
				" | Res_rho=", res_rho, &
				" | Res_u=",   res_u, &
				" | Res_E=",   res_e, &
				" | dt=",      dt
		end if
		
		! Update du vecteur conservatif
		w = w_new
		
	end do
	
	if (iter >= maxiter) then
		print*, "ATTENTION : Pas de convergence apres ", maxiter, " iterations."
		write(*, '(A,1PE10.4,A,1PE10.4,A,1PE10.4)')"Residus: rho = ", res_rho, " | u = ", res_u, " | E = ", res_e
	end if
	
	! Sauvegarde
	call Convert_conservative_to_primal(An,w_new,rho,u,P,E) 		! Convert to primary variables
	call Obtain_data(rho, P, u, Ma, s, Temp)              			! Compute Mach number and entropy and Temperature
	call save_results("data_final",dx,t,An,rho,u,P,E, Ma,s,Temp)    ! Save results
	call create_gnuplot_script(t,"data_final") 						! Creates GNUPLOT script to plot
	
	! Liberation de mémoire
	deallocate(w,w_L,w_R)
	deallocate(w_new)
	deallocate(F_L,F_R)
	deallocate(F_num)
	deallocate(x_i,An,Anf,St)
	deallocate(Ma,s,Temp)
	
    
end subroutine Solve_Eqn_Euler_MUSCL

! -----------------------------
 
 subroutine Compute_time_step(dx,rho,u,E,P,An,w,dt)
! ==================================================
! Calcul le pas de temps pour assurer cdt stabilite
! ==================================================
! Parametres globaux 
use Module_parametres, only : N, gamma, CFL

! Déclarations
	implicit none
    double precision,dimension(:,-1:),intent(in)  :: w						! Vecteur etat
	double precision, intent(in) 				  :: dx						! paramètre flottant d'entrée
    double precision,dimension(-1:),intent(in)    :: An		 			    ! Geometrie de la tuyere					
    double precision,dimension(-1:),intent(inout) :: rho, u, E, P			! Variables primaires
	double precision,dimension(N)			      :: c						! Vitesse du son locale
    integer 									  :: i						! Entier
    double precision 							  :: max_speed				! Calcul de la vit max
	double precision, intent(inout) 			  :: dt						! Pas de tps
! =========================================================================================================================
! Fin des declarations 

   ! Conversion des variables conservées en primitives
    call Convert_conservative_to_primal(An,w,rho,u,P,E)
	
	! Calcul de la vitesse du son
	do i = 1, N
        c(i) = dsqrt(gamma * P(i) / max(rho(i), 1.0d-12)) ! Célérité des ondes acoustiques [m/s]
    end do			
	
	!dt=CFL*dx/max(abs(u)+a);
	max_speed = maxval(dabs(u(1:N)) + c) ! on utilise que les mailles internes
	! Calcul du pas de temps
    dt = CFL * dx / max_speed
	
 end subroutine Compute_time_step
 
  ! -----------------------------
 
 subroutine Calcul_Flux(w_r,F)
! ================================================
! Calcul du vecteur flux F :
! F = {rho*u,rho*u² + P, (rho*E+P)*u}
! ================================================
! Parametres globaux 
use Module_parametres, only : N, gamma

! Déclarations
  implicit none
  double precision, dimension(:,:),intent(in)    :: w_r 		 ! Composantes du vecteur w au temps t{n}
  double precision 					 			 :: P			 ! Pression
  integer 			 							 :: i 			 ! entier pour calcul de boucles
  double precision, dimension(:,:),intent(out)   :: F   		 ! Composantes du vecteur Flux
! =======================================================================================================
! Fin des declarations 
	
	! Calcul des flux
	do i = 1, N+1 
		F(1,i) = w_r(2,i)
		F(2,i) = 0.5d0*(3.0d0-gamma)*(w_r(2,i)**2 / w_r(1,i)) + (gamma-1.0d0)*(w_r(3,i))
		P = (w_r(3,i)-0.5d0*(w_r(2,i)**2 / w_r(1,i)))*(gamma-1.0d0)
		F(3,i) = (w_r(2,i)*w_r(3,i)) / w_r(1,i) + P*w_r(2,i)/w_r(1,i)
	end do
	
 end subroutine Calcul_Flux
 
 ! -----------------------------
 
 subroutine Update_variables(dx,dt,w,F_num,St,rho_s,u_s,E_s,rho,u,P,E,An,Anf,w_new)
! ================================================
! Calcul des variables au prochain temps 
! discrétisation avec Euler explicite ordre 1
! Possibilite d'implicitation : i_impl = 1
! ================================================
! Parametres globaux 
use Module_parametres, only : N, i_impl

! Déclarations
  implicit none
  double precision,intent(in) 					   :: dx,dt,rho_s, u_s, E_s	! paramètres d'entrées de la subroutine
  double precision, dimension(:,-1:),intent(in)    :: w 					! Composantes du vecteur w au temps t{n}
  double precision, dimension(-1:),  intent(in)    :: An   		 			! geom de la tuyere au cells
  double precision,dimension(0:),    intent(in)	   :: Anf	 				! Geometrie de la tuyere au interface
  double precision, dimension(:,:),  intent(in)    :: St   	 				! Source term
  double precision, dimension(:,:),  intent(in)    :: F_num					! Flux numériques
  double precision, dimension(:),    intent(inout) :: rho, u, P, E		    ! Caractéristiques physique du fluide
  integer 			 							   :: i 					! entier pour calcul de boucles
  double precision, dimension(:,-1:),intent(inout) :: w_new 				! Composantes du vecteur w au temps t{n+1}
! =========================================================================================================================
! Fin des declarations

! Validate input
  if (i_impl < 0 .or. i_impl > 1) then
     print *, "Error: Invalid i_impl value. Must be between 0 and 1."
     stop
  end if

! Nozzle geometry creation
select case(i_impl)
	case(0)	
		! EXPLICITE
		! Calcul des états au temps {t+1} (que les mailles internes) avec Euler ordre 1
		do i = 1, N 					
			w_new(:,i) = w(:,i) - (dt / dx) * (Anf(i)*F_num(:,i+1) - Anf(i-1)*F_num(:,i))+ dt*St(:,i)
		end do
	case(1)
		! IMPLICITE : matrix free method (voir papier de Goncalves) 
		call Implicite(dt,dx,rho_s,u_s,E_s,rho,u,P,E,F_num,An,Anf,St,w,w_new)
		
end select
	
 end subroutine Update_variables

! -----------------------------

subroutine Num_flux(w_L,w_R,F_L,F_R,F_num)
! ======================================================================================= 
! Subroutine pour calculer les flux numériques 
! 1) Schéma de Roe
! 2) Schema HLL Davis
! 3) Schema HLL Roe
! 4) Schema HLLE
! 5) Schema HLLC-ANRS (Adaptive Non-iterative Riemann Solver (ANRS))
! ===================================================================================================================================================== 
! Parametres globaux 
use Module_parametres, only : i_sc

  ! Déclaration des variables
	implicit none
	double precision,dimension(:,:),intent(in) 	:: w_R,w_L		! États conservés (rho, rho*u, rho*E) a doite et a gauche de l'interface
	double precision,dimension(:,:),intent(in)	:: F_L,F_R		! Flux
	double precision,dimension(:,:),intent(out) :: F_num 		! Flux num à retourner
! =======================================================================================================================================================
! Fin des declarations 

    ! Initialisation du flux numerique
    F_num = 0.0d0
	
	select case(i_sc)
		case(1)
			! Roe SCHEME
			call ROE_flux(w_L,w_R,F_L,F_R,F_num)
			
		case(2)
			! HLL-Davis SCHEME
			call HLL_flux(w_L,w_R,F_L,F_R,F_num)
			
		case(3)
			! HLLR SCHEME (Hll-Roe)
			call HLLR_flux(w_L,w_R,F_L,F_R,F_num)
			
		case(4)
			! HLLE SCHEME (Hll-Einfield)
			call HLLE_flux(w_L,w_R,F_L,F_R,F_num)
			
		case(5)
			! HLLC-ANRS SCHEME
			call HLLC_ANRS_flux(w_L,w_R,F_L,F_R,F_num)
			
		case(6)
			! HLLC SCHEME
			call HLLC_flux(w_L,w_R,F_L,F_R,F_num)
	end select
	
end subroutine Num_flux
 
! ---------------------------------------------------------------------------------------------------------------------------------------------  
 
subroutine Boundary_conditions(w_s,An,w_new)
! ================================================
! Calcul des BCs
! ================================================
! Parametres globaux 
use Module_parametres, only : N, gamma, i_BC, o_BC, Pext
! Declarations
  implicit none
  double precision,dimension(-1:),intent(in)        :: An		 			  ! Geometrie de la tuyere
  double precision,dimension(:,-1:),intent(in) 		:: w_s					  ! Valeurs initiales pour entree supersonique
  double precision 									:: u_N,rho_N,E_N,c_N,P_N  ! Valeurs primaires a la maille N
  double precision 									:: u_b,rho_b,E_b    	  ! Valeurs primaires au bord (ghost cell)
  double precision 									:: u_s,rho_s,E_s    	  ! Valeurs primaires au bord (ghost cell)
  integer 											:: i 					  ! Integer for loop
  double precision, dimension(:,-1:),intent(inout)  :: w_new 				  ! Variables des lois de conservations au temps t{n+1}
! =========================================================================================================================
! Fin des declarations 

	! Validate input
	if (i_BC < 1 .or. i_BC > 2) then
		print *, "Error: Invalid i_BC value. Must be between 1 and 2."
		stop
	end if 
    ! Validate input
	if (o_BC < 1 .or. o_BC > 2) then
		print *, "Error: Invalid o_BC value. Must be between 1 and 2."
		stop
	end if 

	! Inlet BC
	select case(i_BC)
		case(1)
			! Entree subsonique (cdt compabilite + 2 var imposée)
			call Newton(An,w_new,u_b,rho_b,E_b) ! Newton pour determiner u_b
			w_new(1,0) = rho_b*An(0)
			w_new(2,0) = rho_b*u_b*An(0)
			w_new(3,0) = rho_b*E_b*An(0)
			! Interpolation ordre 0 pour le ghost cell -1
			w_new(1,-1) = rho_b*An(-1)
			w_new(2,-1) = rho_b*u_b*An(-1)
			w_new(3,-1) = rho_b*E_b*An(-1)
		case(2)
			! Entree supersonique (Tout est imposée)
			do i = -1,0
				rho_s = w_s(1,i) / An(i)
				u_s   = w_s(2,i) / (An(i) * rho_s)
				E_s   = w_s(3,i) / (An(i) * rho_s)
				w_new(1,i) = rho_s*An(i)
				w_new(2,i) = rho_s*u_s*An(i)
				w_new(3,i) = rho_s*E_s*An(i)
			enddo
	end select
	
	! Outlet BC
	select case(o_BC)
		case(1)
			! Sortie subsonique (cdt compabilite + Pext imposée)
			! --- Extraction des variables à la première maille réelle (i=N) ---
			rho_N = w_new(1, N)/An(N)
			u_N   = w_new(2, N) / max(An(N)*rho_N, 1.0d-12) 		        ! Sécurité contre la division par zéro ou densité négative
			E_N   = w_new(3, N) / max(An(N)*rho_N, 1.0d-12) 		  		! E total par unité de masse
			P_N   = (gamma - 1.0d0) * rho_N * (E_N - 0.5d0 *u_N**2)  		! Calcul de la pression (Loi des gaz parfaits)
			c_N   = dsqrt(gamma * P_N / max(rho_N, 1.0d-12)) 		  	    ! Célérité des ondes acoustiques
			
			! Calcul des variables au bord N+1
			u_b   = u_N - (Pext-P_N)/(rho_N*c_N) ! Cdt de compatibilite
			rho_b = rho_N + (Pext-P_N)/(c_N**2)  ! Compute Energy
			E_b = Pext / (rho_b * (gamma - 1.0d0)) + 0.5d0 * u_b**2
			
			! Update the w_new vector at the ghost cells
			w_new(1, N+1) = rho_b*An(N+1)
			w_new(2, N+1) = rho_b*u_b*An(N+1)
			w_new(3, N+1) = rho_b*E_b*An(N+1)
			! Interpolation ordre 0 pour le ghost cell N+2
			w_new(1,N+2) = rho_b*An(N+2)
			w_new(2,N+2) = rho_b*u_b*An(N+2)
			w_new(3,N+2) = rho_b*E_b*An(N+2)
		
		case(2)
			! Sortie supersonique (Rien n'est imposée : Neumann)
			rho_N = w_new(1, N)/An(N)
			u_N   = w_new(2, N) / max(An(N)*rho_N, 1.0d-12) 		  ! Sécurité contre la division par zéro ou densité négative
			E_N   = w_new(3, N) / max(An(N)*rho_N, 1.0d-12) 		  ! E total par unité de masse
			
			! Update the w_new vector at the ghost cells
			w_new(1, N+1) = rho_N*An(N+1)
			w_new(2, N+1) = rho_N*u_N*An(N+1)
			w_new(3, N+1) = rho_N*E_N*An(N+1) 
			! Interpolation ordre 0 pour le ghost cell N+2
			w_new(1,N+2) = rho_N*An(N+2)
			w_new(2,N+2) = rho_N*u_N*An(N+2)
			w_new(3,N+2) = rho_N*E_N*An(N+2)
			
	end select

 end subroutine Boundary_conditions
 
 ! ---------------------------------------------------------------------------------------------------------------------------------------------  
 
subroutine Convert_conservative_to_primal(An,w_new,rho,u,P,E)
! =========================================================
! Convertit les variables d'état en variables primaires
! =========================================================
! Parametres globaux 
use Module_parametres, only : N, gamma

 ! Déclaration des variables
  implicit none
  double precision, dimension(:,-1:),intent(in) :: w_new 			! Variables des lois de conservations au temps t{n+1}
  double precision,dimension(-1:),intent(in)	:: An	 			! Geometrie de la tuyere au cell
  integer										:: j				! loop variables
  double precision, dimension(-1:),intent(out)  :: rho,u,P,E		! Variables primitives
! =========================================================================================================================
! Fin des declarations 

  ! Conversion des variables conservées en primitives
  do j = -1, N+2
     rho(j) = w_new(1, j)/An(j)
     u(j)   = w_new(2, j) / max(An(j)*rho(j), 1.0d-12) ! Sécurité contre la division par zéro ou densité négative
     ! E total par unité de masse
     E(j)   = w_new(3, j) / max(An(j)*rho(j), 1.0d-12) 
     ! Calcul de la pression (Loi des gaz parfaits)
     P(j)   = (gamma - 1.0d0) * rho(j) * (E(j) - 0.5d0 *u(j)**2)
     
     ! Alerte si pression négative 
     if (P(j) < 0.0d0) then
		print*, "Erreur physique majeure : Pression negative j =", j, " P =", P(j)
		stop "Calcul interrompu"
	 endif
  end do

 end subroutine Convert_conservative_to_primal
 
 ! -----------------------------
 
 
 subroutine Interp(rho,u,P,E,An,w,w_L,w_R)
! ======================================================================================
! Calcul les valeurs des variables aux interfaces
! avec reconstruction MUSCL : 
! ghost cell | Cell i-1 | Cell i | Cell i+1 | Cell i+2 | Cell i+3 | ghost cell
!     ...        ...       ...        ...        ...
!    w(i-2)     w(i-1)     w(i)     w(i+1)     w(i+2)
!
! Interfaces (flux calculés ici) : 
!       i-1/2     i+1/2      i+3/2      i+5/2
! Pour chaque interface on construit w_L et w_R qui sont ensuite donné au schéma de Roe
! ======================================================================================
! Parametres globaux 
use Module_parametres, only : i_ord, N, gamma, b, phi
! Declarations
  implicit none
  double precision, dimension(:,-1:),intent(in)   :: w 				     ! Composantes du vecteur w au temps t{n}
  double precision,dimension(-1:),intent(in)	  :: An	 				 ! Geometrie de la tuyere au cell
  double precision							      :: eps,c1,c2,psi		 ! Tolerance to avoid division by 0, coefficients
  integer										  :: j					 ! loop variables 
  double precision								  :: rho_L,u_L,P_L,E_L 	 ! Left reconstruction of primary variables
  double precision								  :: rho_r,u_R,P_R,E_R 	 ! Right reconstruction of primary variables
  double precision								  :: r1,r2,r3,r4 	 	 ! Local slope of a cell
  double precision								  :: l1,l2,l3,l4	     ! final slope after limiter : X*PSI(Y/X) with Y/X : slope    
  double precision, dimension(0:N+2)			  :: dif_r,dif_u,dif_P	 ! Difference of neighbours cells
  double precision, dimension(-1:),intent(inout)  :: rho, u, P, E		 ! Caractéristiques physique du fluide
  double precision, dimension(:,:),intent(inout)  :: w_L,w_R 		     ! Composantes du vecteur w au temps t{n+1}
! =========================================================================================================================
! Fin des declarations 

   ! Validate input
  if (i_ord < 1 .or. i_ord > 2) then
    print *, "Error: Invalid i_ord value. Must be between 1 and 2."
    stop
  end if
   
    select case(i_ord)
	case(1)
		! ordre 1 
		call Convert_conservative_to_primal(An,w,rho,u,P,E) ! conversion 
		! Interpolation aux interfaces
		do j=1,N+1
			! Right interface
			w_R(1,j) = rho(j)
			w_R(2,j) = rho(j)*u(j)
			w_R(3,j) = rho(j)*E(j)
			! Left interface
			w_L(1,j) = rho(j-1)
			w_L(2,j) = rho(j-1)*u(j-1)
			w_L(3,j) = rho(j-1)*E(j-1)
		enddo
		
	case(2)
		! Reconstruction MUSCL (Ref : chap2 MNA ensma)
		
		! Initialisation
		r1 = 0.0d0; r2 = 0.0d0; r3 = 0.0d0; r4 = 0.0d0
		l1 = 0.0d0; l2 = 0.0d0; l3 = 0.0d0; l4 = 0.0d0
		w_L = 0.0d0;w_R = 0.0d0; eps = 1.0d-8; psi = 0.0d0
		
		! Conversion en variables primaires
		call Convert_conservative_to_primal(An,w,rho,u,P,E) 
		
		! Compute coefficients
		c1 = (1.0d0-phi)/4.0d0
		c2 = (1.0d0+phi)/4.0d0
		
		! Compute difference of each cells
		do j = 0,N+2
			dif_r(j) = rho(j)-rho(j-1)
			dif_u(j) = u(j)-u(j-1)
			dif_P(j) = P(j)-P(j-1)
		enddo
		
		! Compute Left and right state for each interface and each variable
		do j = 1,N+1
			! Compute slopes for rho
			r1=  dif_r(j-1)/(b*dif_r(j)+eps) ! eps : avoid 0 at denom
			r2 = dif_r(j)/(b*dif_r(j-1)+eps)
			r3 = dif_r(j)/(b*dif_r(j+1)+eps)
			r4 = dif_r(j+1)/(b*dif_r(j)+eps)
			! Compute limiter
			call limiter_function(r1,psi)
			l1 = (b*dif_r(j))*psi
			call limiter_function(r2,psi)
			l2 = (b*dif_r(j-1))*psi
			call limiter_function(r3,psi)
			l3 = (b*dif_r(j+1))*psi
			call limiter_function(r4,psi)
			l4 = (b*dif_r(j))*psi
			! Compute Left and right states of the interface (rho)
			rho_L = rho(j-1) + c1*l1 + c2*l2
			rho_R = rho(j) - c2*l3 - c1*l4
			! ---------------------------------------------
			! Compute slopes for u
			r1=  dif_u(j-1)/(b*dif_u(j)+eps) ! eps : avoid 0 at denom
			r2 = dif_u(j)/(b*dif_u(j-1)+eps)
			r3 = dif_u(j)/(b*dif_u(j+1)+eps)
			r4 = dif_u(j+1)/(b*dif_u(j)+eps)
			! Compute limiter
			call limiter_function(r1,psi)
			l1 = (b*dif_u(j))*psi
			call limiter_function(r2,psi)
			l2 = (b*dif_u(j-1))*psi
			call limiter_function(r3,psi)
			l3 = (b*dif_u(j+1))*psi
			call limiter_function(r4,psi)
			l4 = (b*dif_u(j))*psi
			! Compute Left and right states of the interface (u)
			u_L = u(j-1) + c1*l1 + c2*l2
			u_R = u(j) - c2*l3 - c1*l4
			! ---------------------------------------------
			! Compute slopes for P
			r1=  dif_P(j-1)/(b*dif_P(j)+eps) ! eps : avoid 0 at denom
			r2 = dif_P(j)/(b*dif_P(j-1)+eps)
			r3 = dif_P(j)/(b*dif_P(j+1)+eps)
			r4 = dif_P(j+1)/(b*dif_P(j)+eps)
			! Compute limiter
			call limiter_function(r1,psi)
			l1 = (b*dif_P(j))*psi
			call limiter_function(r2,psi)
			l2 = (b*dif_P(j-1))*psi
			call limiter_function(r3,psi)
			l3 = (b*dif_P(j+1))*psi
			call limiter_function(r4,psi)
			l4 = (b*dif_P(j))*psi
			! Compute Left and right states of the interface (P)
			P_L = P(j-1) + c1*l1 + c2*l2
			P_R = P(j) - c2*l3 - c1*l4
			! ---------------------------------------------
			! Safety (No negative value)
			rho_L = max(rho_L, 1.d-12)
			rho_R = max(rho_R, 1.d-12)
			P_L   = max(P_L, 1.d-12)
			P_R   = max(P_R, 1.d-12)
			! ---------------------------------------------
			! Compute Left and right states of the interface (P)
			E_L = (P_L/(rho_L*(gamma-1.0d0)))+0.5d0*u_L**2.0d0
			E_R = (P_R/(rho_R*(gamma-1.0d0)))+0.5d0*u_R**2.0d0
			! ---------------------------------------------
			! Compute Left interface w_L
			w_L(1,j) = rho_L
			w_L(2,j) = rho_L*u_L
			w_L(3,j) = rho_L*E_L
			! Compute Right interface w_L
			w_R(1,j) = rho_R
			w_R(2,j) = rho_R*u_R
			w_R(3,j) = rho_R*E_R
			
		enddo
		
    end select
  
 end subroutine Interp


! -----------------------------

subroutine limiter_function(r,psi)
! ================================================
! Limiteurs 
! ================================================
! Parametres globaux 
use Module_parametres, only : i_fct, beta

  implicit none
  double precision,INTENT(IN)	   :: r    		 	 ! Vector of the slope
  double precision 			   	   :: ii,jj,kk,eps	 ! Local values
  double precision,intent(inout)   :: psi   		 ! Return value
! =========================================================================================================================
! Fin des declarations 

! Validate input
  if (i_fct < 1 .or. i_fct > 5) then
     print *, "Error: Invalid i_fct value. Must be between 1 and 5."
     stop
  end if

! Select the limiter function, 1: Minmod, 2 : VanLeer, 3 : VanAlbada, 4 : Superbee, 5 : Chakravarthy
select case(i_fct)
	case(1)
		! Minmod (most robust but most dissipative)
		ii = min(1.0d0,r)
		psi = max(0.0d0,ii)
	case(2)
		! VanLeer
		eps = 1.0d-8
		psi = (r+dabs(r))/(1+r+eps) ! avoid 0 division
	case(3)
		! VanAlbada
		ii = (r+r**2)/(1+r**2)
		psi = max(0.0d0,ii)
	case(4)
		! Superbee (captures stiff front shock but less robust)
		ii = min(2.0d0,r)
		jj = min(1.0d0,2.0d0*r)
		kk = max(ii,jj)
		psi = max(0.0d0,kk)
	case(5)
		! Chakravarthy
		ii = min(beta,r)
		psi = max(ii,0.0d0)
end select

END subroutine limiter_function

! -----------------------------

subroutine Obtain_data(rho, P, u, Ma, s, Temp)
! ================================================
! Calcul du Mach, Entropie et Température (SI)
! ================================================
use Module_parametres, only : N, gamma, r, Cp
! Déclaration 
  implicit none
  double precision, dimension(-1:),intent(in)  :: rho,u,P
  double precision                             :: c_sound
  integer                                      :: i
  double precision, dimension(:),intent(out)   :: Ma, s, Temp
! ================================================
! Fin des declarations 

do i = 1, N
    ! Température réelle en Kelvins : T = P / (rho * R)
    Temp(i) = P(i) / (max(rho(i), 1.0d-12) * r)
    
    ! Vitesse du son locale
    c_sound = dsqrt(gamma * r * Temp(i))
    
    ! 3. Nombre de Mach
    Ma(i) = abs(u(i)) / max(c_sound, 1.0d-12)
    
    ! "Entropie" (Indicateur isentropique : P/rho^gamma)
    s(i) = P(i) / (max(rho(i), 1.0d-12)**gamma)
end do

END subroutine Obtain_data

! -----------------------------

subroutine Source_term(dx,w,An,Anf,St)
! ================================================
! Compute source term
! ================================================
! Parametres globaux 
use Module_parametres, only : N, gamma, Cf, h, Tp, r
! Declarations
  implicit none
  double precision,intent(in)                  :: dx 	    				! Pas de discretisation
  double precision,dimension(-1:),intent(in)   :: An		 				! Geometrie de la tuyere au centre
  double precision,dimension(0:),intent(in)	   :: Anf		 	    		! Geometrie de la tuyere au interface
  double precision, dimension(:,-1:),intent(in):: w		 					! Variables des lois de conservations au temps t{n}
  double precision							   :: dlnA_dx,pi,d_sig,T		! Variables double locale
  double precision							   :: u2_u1,J_aire,J_frot		! Variables double locale
  double precision 							   :: P_local 					! Pression locale
  integer 			   		                   :: i			        		! Indice boucle
  double precision, dimension(:,:),intent(out) :: St   			    		! Source term
! =========================================================================================================================
! Fin des declarations 

! --- Définition de constantes  ---
pi = 4.0d0*datan(1.0d0)

! --- St(1) : Continuité ---
St(1,:) = 0.0d0

do i=1,N

	! --- Pré-calculs  ---
	dlnA_dx = (dlog(Anf(i)) - dlog(Anf(i-1)))/dx
	d_sig = 2.0d0*dsqrt(pi*An(i))
	u2_u1  = (w(2,i)**2) / w(1,i)  ! Terme rho*u^2*A
	! P = (gamma-1) * [ w3/A - 0.5 * (w2^2 / (w1*A)) ]
	P_local = (gamma - 1.0d0) * (w(3,i)/An(i) - 0.5d0 * (u2_u1/ An(i)))
	
	! --- St(2) : Quantité de mouvement ---
	! Partie Aire (Forme log-linéaire)
	J_aire = P_local * An(i) * dlnA_dx
	! Partie Frottement
	J_frot = -0.5d0 * (u2_u1 / An(i)) * Cf * d_sig
	St(2,i) = J_aire + J_frot

	! --- St(3) : Énergie ---
	T = P_local / ( (w(1,i)/An(i)) * r ) ! T = P / (rho * R)  
	St(3,i) = -h * (T - Tp) * d_sig
	
enddo
END subroutine Source_term

! -----------------------------


subroutine Compute_residus(An,w,w_new,res_rho,res_u,res_e)
! ================================================
! Compute residus
! ================================================
! Parametres globaux 
use Module_parametres, only : N
! Declarations
 implicit none
 double precision,dimension(-1:),intent(in)	     :: An		 	        	! Geometrie de la tuyere
 double precision, dimension(:,-1:),intent(in)   :: w 				    	! Composantes du vecteur w au temps t{n}
 double precision, dimension(:,-1:),intent(in)   :: w_new 					! Composantes du vecteur w au temps t{n}
 double precision                                :: norm_rho,norm_u,norm_e  ! Normes
 integer										 :: i						! Calcul de boucle
 double precision				                 :: rho_np1,rho_n  			! Variables locales
 double precision				                 :: u_np1,u_n  				! Variables locales
 double precision				                 :: e_np1,e_n  				! Variables locales
 double precision,intent(out)                    :: res_rho,res_u,res_e 	! Residus
! =========================================================================================================================
! Fin des declarations 

! Init
res_rho = 0.0d0 
res_u   = 0.0d0 
res_e   = 0.0d0 
norm_rho = 0.0d0
norm_u   = 0.0d0
norm_e   = 0.0d0

do i = 1, N
    ! --- Temps n ---
    rho_n = w(1,i) / An(i)
    u_n   = w(2,i) / max(An(i) * rho_n, 1.0d-12)
    e_n   = w(3,i) / max(An(i) * rho_n, 1.0d-12)

    ! --- Temps n+1 ---
    rho_np1 = w_new(1,i) / An(i)
    u_np1   = w_new(2,i) / max(An(i) * rho_np1, 1.0d-12)
    e_np1   = w_new(3,i) / max(An(i) * rho_np1, 1.0d-12)

    ! --- Somme des carrés des écarts ---
    res_rho = res_rho + (rho_np1 - rho_n)**2
    res_u   = res_u   + (u_np1   - u_n)**2
    res_e   = res_e   + (e_np1   - e_n)**2

    ! --- Somme des carrés pour normalisation ---
    norm_rho = norm_rho + rho_np1**2
    norm_u   = norm_u   + u_np1**2
    norm_e   = norm_e   + e_np1**2
end do

! --- Norme L2 Relative (Sans dimension) ---
res_rho = dsqrt(res_rho / max(norm_rho, 1.0d-12))
res_u   = dsqrt(res_u   / max(norm_u,   1.0d-12))
res_e   = dsqrt(res_e   / max(norm_e,   1.0d-12))

END subroutine Compute_residus

! -----------------------------


subroutine Newton(An,w_new,u_b,rho_b,E_b)
! ================================================
! Algo de Newton pour resoudre f(ub) = 0 
! non lineaire
! ================================================
! Parametres globaux 
use Module_parametres, only : N, gamma, Cp, r, Pt0, Tt0 
! Declarations
 implicit none
 double precision, dimension(:,-1:),intent(in) :: w_new 				! Variables des lois de conservations au temps t{n+1}
 double precision,dimension(-1:),intent(in)	   :: An		 	        ! Geometrie de la tuyere
 double precision		      				   :: rho1,c1,u1,P1,E1		! Variables primaires a la maille i=1
 double precision 							   :: u_next,error,tol,g_b  ! Doubles locaux
 double precision 							   :: f,f_prime,P_b			! Doubles locaux
 integer 							  		   :: max_newton,iter		! Entier locaux
 double precision,INTENT(out)         		   :: u_b,rho_b,E_b		    ! Variables primaires a la maille i=0
! =========================================================================================================================
! Fin des declarations 

	! --- Extraction des variables à la première maille réelle (i=1) ---
	rho1 = w_new(1, 1)/An(1)
	u1   = w_new(2, 1) / max(An(1)*rho1, 1.0d-12) 			   ! Sécurité contre la division par zéro ou densité négative
	E1   = w_new(3, 1) / max(An(1)*rho1, 1.0d-12) 			   ! E total par unité de masse
	P1   = (gamma - 1.0d0) * rho1 * (E1 - 0.5d0*u1**2)         ! Calcul de la pression (Loi des gaz parfaits)
	c1 = dsqrt(gamma * P1 / max(rho1, 1.0d-12)) 			   ! Célérité des ondes acoustiques
   
	! Initialisation du Newton
	u_b = u1 
	iter = 0
	tol = 1.0d-8
	max_newton = 10
	error = 1.0d0

	! --- 2. Boucle de Newton-Raphson ---
	do while (error > tol .and. iter < max_newton)
		iter = iter + 1
		
		! Calcul du ratio de température (g_b est sans unité, entre 0 et 1)
		g_b = 1.0d0 - (u_b**2 / (2.0d0 * Cp * Tt0))
		
		! Securite
		if (g_b <= 1.0d-6) then
			write(*,*) 'ERREUR FATALE : g_b est negatif ou nul (vitesse ub excessive)'
			write(*,*) 'Divergence probable. Verifiez votre CFL ou votre initialisation.'
			stop
		endif
		
		! Pression et densité au bord respectivement en Pa et (en kg/m3)
		P_b = Pt0 * (g_b**(gamma/(gamma-1.0d0)))
		rho_b = (Pt0 / (r * Tt0)) * (g_b**(1.0d0/(gamma-1.0d0)))
		
		! Fonction f(ub) et sa dérivée f'(ub)  = -rho_b * u_b - (rho1 * c1)
		f = P_b - P1 - (rho1 * c1) * (u_b - u1)
		f_prime = -rho_b * u_b - (rho1 * c1)
		
		! Mise à jour
		u_next = u_b - f / f_prime
		error = dabs(f) / (rho1 * c1 + 1.0d-12) ! normalise
		u_b = u_next
		
	end do
	! Si non convergence
	if(iter==max_newton)then
		write(*,*)'ATTENTION NEWTON N A PAS CONVERGE'
	endif
	
	! Après la boucle Newton, recalculer les valeurs finales
    g_b = 1.0d0 - (u_b**2 / (2.0d0 * Cp * Tt0))
    P_b = Pt0 * (g_b**(gamma/(gamma-1.0d0)))
    rho_b = (Pt0 / (r * Tt0)) * (g_b**(1.0d0/(gamma-1.0d0)))
	E_b = (P_b/(rho_b*(gamma-1.0d0)))+0.5d0*u_b**2.0d0

END subroutine Newton

! -----------------------------

subroutine Implicite(dt,dx,rho_s,u_s,E_s,rho,u,P,E,F_num,An,Anf,St,w,w_new)
! ================================================
! Phase d'implicitation (matrix free)
! voir EulerQuasi1D_Goncalves.pdf
! ================================================
! Parametres globaux 
use Module_parametres, only : N, gamma 
! Declarations
 implicit none
 ! Variables d'entree
 double precision,intent(in) 					:: dx,dt,rho_s,u_s,E_s						! paramètres d'entrées de la subroutine
 double precision, dimension(:,-1:),intent(in)  :: w 								 		! Composantes du vecteur w au temps t{n}
 double precision,dimension(-1:),intent(in)	    :: An		 	    						! Geometrie de la tuyere aux centre des faces
 double precision, dimension(0:), intent(in)    :: Anf  		 							! geom de la tuyere au Interfaces
 double precision, dimension(:,:),intent(in)    :: St   			    					! Source term
 double precision, dimension(:,:),intent(in)    :: F_num									! Flux numériques
 double precision, dimension(-1:),intent(inout) :: rho, u, P, E								! Caractéristiques physique du fluide
 ! Variables locales 																		! ---------------
 double precision, dimension(3,N+1)		        :: dF,F_imp,w_L_imp,w_R_imp,F_R_imp,F_L_imp	! Matrices locales
 double precision, dimension(3,-1:N+2)		    :: w_imp									! Vecteur d'etat implicite 
 double precision, dimension(-1:N+2)		    :: c										! Vitesse du son
 double precision, dimension(N)		   	        :: C_minus,C_plus,C_0						! Vecteurs locaux
 double precision, dimension(3)		   	        :: diffusion,transport					   	! Vecteurs locaux
 double precision, dimension(3,N)		   	    :: dw_exp,dw_impl,dw_old					! Matrices locales 
 double precision 							    :: u_i_d,c_i_d,u_i_g,c_i_g,sigma,lambda,u_inter,omega	    	! Doubles locaux
 integer 							  		    :: i, l, n_imp, j							! Entier locaux
 ! Variables de sortie																		! ---------------
 double precision, dimension(:,-1:),intent(out) :: w_new 									! Variables des lois de conservations au temps t{n+1}
! =========================================================================================================================
! Fin des declarations 
	
	! Conversion variables conservative en primaires
	call Convert_conservative_to_primal(An,w,rho,u,P,E) 
	! Calcul des célérités des ondes acoustiques
	c = dsqrt(gamma * P / max(rho, 1.0d-12)) 			   
	
	! -- Calcul des coefficients C -- 
	sigma = dt/dx
	do i = 1, N ! boucle sur mailles
		u_i_g = 0.5d0*(u(i-1)+u(i)) 	  ! vitesse materielle à l'interface gauche
		c_i_g = 0.5d0*(c(i-1)+c(i)) 	  ! vitesse acoustique à l'interface gauche
		u_i_d = 0.5d0*(u(i+1)+u(i)) 	  ! vitesse materielle à l'interface droite
		c_i_d = 0.5d0*(c(i+1)+c(i)) 	  ! vitesse acoustique à l'interface droite
		C_minus(i) = -0.5d0*sigma*(dabs(u_i_g)+c_i_g)
		C_plus(i) = -0.5d0*sigma*(dabs(u_i_d)+c_i_d)
		C_0(i) = 1.0d0 - C_minus(i) - C_plus(i)
	enddo
	
	! -- Calcul increment explicite -- 
	do i = 1, N 					
		dw_exp(:,i) =  - sigma * (Anf(i)*F_num(:,i+1) - Anf(i-1)*F_num(:,i))+ dt*St(:,i)
	enddo
	
	dw_impl = dw_exp 
	n_imp = 20 ! nb de loop
	! -- Boucle implicitation -- 
	do l= 1, n_imp
		dw_old = dw_impl
		!w_imp(:,1:N) = w(:,1:N) + dw_impl(:,1:N) 			 ! Construction du w implicite
		!call Boundary_conditions(rho_s,u_s,E_s,An,w_imp) 	 ! BC
		! Calcul approche jacobienne du flux
		!~do j = 1,N+1									     ! Boucle sur interfaces
		!~	lambda = dmax1(dabs(u(j-1)) + c(j-1), dabs(u(j)) + c(j))
		!~	u_inter = 0.5d0 * ( (dabs(u(j-1)) + c(j-1)) + (dabs(u(j)) + c(j)) )
		!~	if (j == 1) then 
		!~		dF(:,j) = u_inter * 0.5d0 * (dw_impl(:,1) + dw_impl(:,j)) 
		!~	else if (j == N+1) then
		!~		dF(:,j) = u_inter * 0.5d0 * (dw_impl(:,j-1) + dw_impl(:,N)) 
		!~	else
		!~		dF(:,j) = u_inter * 0.5d0 * (dw_impl(:,j-1) + dw_impl(:,j))
		!~	endif
		!~	!~write(*,*)'dF 1 : ',dF(1,j)
		!~	!~write(*,*)'dF 2 : ',dF(2,j)
		!~	!~write(*,*)'dF 3 : ',dF(3,j)
		!~enddo
		
		! --- ZONE DE TEST POUR ISOLER LE BUG ---
        ! On commente temporairement les gros appels de flux et on calcule un dF doux
        do j = 2, N
            ! Vitesse moyenne à l'interface j (entre maille j-1 et j)
            u_inter = 0.5d0 * (u(j-1) + u(j))
            ! Approximon dF par un transport linéaire simple
            dF(:,j) = u_inter * 0.5d0 * (dw_impl(:,j-1) + dw_impl(:,j))
        end do
        ! Conditions aux limites basiques pour dF aux frontières j=1 et j=N+1
        dF(:,1)   = dF(:,2)
        dF(:,N+1) = dF(:,N)
        ! -----------------------------------------
		
		do i = 1, N
			! Terme de transport (flux centrés aux interfaces de la maille i)
			!transport(:) = 0.5d0 * sigma * ( dF(:, i+1) - dF(:, i) )
			transport(:) = 0.0d0
			!~write(*,*)'Transport 1 : ',transport(1)
			!~write(*,*)'Transport 2 : ',transport(2)
			!~write(*,*)'Transport 3 : ',transport(3)
			! Terme de diffusion 
			diffusion(:) = 0.0d0
			if (i > 1) then
				diffusion(:) = diffusion(:) + C_minus(i) * dw_old(:, i-1)
			end if
			if (i < N) then
				diffusion(:) = diffusion(:) + C_plus(i) * dw_old(:, i+1)
			end if
			!~write(*,*)'diffusion 1 : ',diffusion(1)
			!~write(*,*)'diffusion 2 : ',diffusion(2)
			!~write(*,*)'diffusion 3 : ',diffusion(3)
			! Formule finale de mise à jour de l'incrément pour le tour l+1
			dw_impl(:, i) = ( dw_exp(:, i) - transport(:) - diffusion(:) ) / C_0(i)
			!~write(*,*)'dw_impl 1,i : ',dw_impl(1,i)
			!~write(*,*)'dw_impl 2,i : ',dw_impl(2,i)
			!~write(*,*)'dw_impl 3,i : ',dw_impl(3,i)
			
			! test
			!write(*,*) l,maxval(abs(dw_impl-dw_old))
			
		enddo
    enddo
	!~write(*,*) 'C0 min/max',minval(C_0),maxval(C_0)
	!~write(*,*) 'max rhoE=',maxval(abs(w(3,1:N)))
	!~write(*,*) 'max dwE=',maxval(abs(dw_impl(3,1:N)))
	!~write(*,*) 'max dw rho = ', maxval(abs(dw_impl(1,:)))
	!~write(*,*) 'max dw mom = ', maxval(abs(dw_impl(2,:)))
	!~write(*,*) 'max dw E   = ', maxval(abs(dw_impl(3,:)))
	omega = 0.3d0 ! facteur de relaxation
	w_new(:,1:N) = w(:,1:N) + omega*dw_impl(:,1:N) ! Construction du vecteur final

END subroutine Implicite

! -----------------------------


end module Module_solve
