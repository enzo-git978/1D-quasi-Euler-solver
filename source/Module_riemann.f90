module Module_riemann

contains

! ----------------------------------------------------------------------
! METHODES
! ----------------------------------------------------------------------


subroutine Convert_conservative_to_primal_v2(w_r,rho,u,P,E)
! =========================================================
! Convertit les variables d'état en variables primaires
! =========================================================
! Parametres globaux 
use Module_parametres, only : N, gamma

! Déclarations
  implicit none
  double precision, dimension(:,:),intent(in)   :: w_r 				! Variables des lois de conservations au temps t{n+1}
  integer										:: j				! loop variables
  double precision, dimension(:),intent(out)    :: rho,u,P,E		! Variables primitives
! =========================================================================================================================
! Fin des declarations 
 

! Conversion des variables conservées en primitives
  do j = 1, N+1
     rho(j) = w_r(1, j)
     u(j)   = w_r(2, j) / max(rho(j), 1.0d-12) ! Sécurité contre la division par zéro ou densité négative
     ! E total par unité de masse
     E(j)   = w_r(3, j) / max(rho(j), 1.0d-12) 
     ! Calcul de la pression (Loi des gaz parfaits)
     P(j)   = (gamma - 1.0d0) * rho(j) * (E(j) - 0.5d0 * u(j)**2)
  end do
  
 end subroutine Convert_conservative_to_primal_v2
 
 ! -----------------------------

subroutine ROE_flux(w_L,w_R,F_L,F_R,F_num)
!=====================================================================
!  Compute Numerical fluxes with Roe SCHEME
!=====================================================================
! Parametres globaux 
use Module_parametres, only : N, gamma, delta_star

! Déclarations
  implicit none
  double precision,dimension(:,:),intent(in) 	  :: w_R,w_L		    	! États conservés (rho, rho*u, rho*E) a doite et a gauche de l'interface
  double precision,dimension(:,:),intent(in)	  :: F_L,F_R	        	! Flux
  integer 										  :: j				    	! Entier pour boucle
  double precision,dimension(N+1)	              :: rho_L,u_L,P_L,E_L		! Variables primitives a gauche
  double precision,dimension(N+1)	              :: rho_R,u_R,P_R,E_R		! Variables primitives a droite
  double precision,dimension(N+1)	              :: H_R,H_L			    ! Enthalpie totale a gauche et droite
  double precision 								  :: u_moy, H_moy, c_moy,a	! Valeurs moyennes de Roe
  double precision,dimension(3)					  :: wdif			 		! Vecteurs de 3 composants
  double precision,dimension(3,3)				  :: A_tilde				! Matrice A_tilde
  double precision,dimension(:,:),intent(inout)   :: F_num 				    ! Flux num à retourner
 ! =========================================================================================================================
! Fin des declarations 

	! Calcul du left state
	call Convert_conservative_to_primal_v2(w_L,rho_L,u_L,P_L,E_L)
	
	! Calcul du right state
	call Convert_conservative_to_primal_v2(w_R,rho_R,u_R,P_R,E_R)
	
	! Calcul de l'enthalpie totale a gauche
	H_L = (gamma/(gamma-1))*(P_L/rho_L)+0.5d0*u_L**2
	
	! Calcul de l'enthalpie totale a droite
	H_R = (gamma/(gamma-1))*(P_R/rho_R)+0.5d0*u_R**2
	
	! Boucle sur les interfaces
	do j = 1, N + 1
		! Calcul des moyennes de Roe
		a = dsqrt(rho_R(j)/ rho_L(j))									! rapport des masses volumiques
		u_moy = (a*u_R(j)+u_L(j))/(1.0d0+a)								! Moyenne de Roe de la vitesse matérielle
		H_moy = (a*H_R(j)+H_L(j))/(1.0d0+a)								! Moyenne de Roe de l'enthalpie	
		c_moy = dsqrt((gamma - 1.0d0) * (H_moy - 0.5d0 * u_moy**2))		! Moyenne de Roe de la vitesse de l'onde acoustique
		
		! Construction des matrices (voir excellent cours : P.S volpiani Roe's scheme)
		call Construct_matrix(u_moy, H_moy, c_moy, A_tilde)
		
		! Différence des états 
		wdif(:) = w_R(:,j) - w_L(:,j) 
	 
		! Roe : 0.5*[F(w_L)+F(w_R) - |A_tilde| * (W_R - W_L)]
		F_num(:,j) = 0.5d0*(F_L(:,j) + F_R(:,j)) - 0.5d0*matmul(A_tilde, wdif(:))

	end do


 end subroutine ROE_flux
 
 ! -----------------------------
 
 subroutine Construct_matrix(u_moy, H_moy, c_moy, A_tilde)
! =====================================================================================================
! Subroutine pour calculer les composantes des matrices A_tilde = P_tilde * Delta_tilde * P_tilde ^-1
! ======================================================================================================
! Parametres globaux 
use Module_parametres, only : i_cor, delta_star, gamma

! Déclarations
	implicit none
    double precision,intent(in)	 				   :: u_moy, H_moy, c_moy		        ! Valeurs moyennes de Roe
	double precision							   :: alph1, alph2, delta				! Variables locales
	double precision							   :: lambda_1,lambda_2,lambda_3		! Variables locales
	double precision							   :: psi_1,psi_2,psi_3					! Variables locales
    double precision,dimension(3,3) 			   :: Delta_tilde, Pinv, P_tilde		! Matrices pour Diagonalisation
	double precision,dimension(:,:),intent(inout)  :: A_tilde							! Matrice A_tilde a retourner
! =========================================================================================================================
! Fin des declarations   

    ! Variables auxiliaires pour P^{-1}
    alph1 = (gamma - 1.0d0) * u_moy**2 / (2.0d0 * c_moy**2)
    alph2 = (gamma - 1.0d0) / (c_moy**2)
	
    ! Matrice P^{-1}
    Pinv(1, 1) = 0.5d0 * (alph1 + u_moy / c_moy)
    Pinv(1, 2) = -0.5d0 * (alph2 * u_moy + 1.0d0 / c_moy)
    Pinv(1, 3) = alph2 / 2.0d0
    Pinv(2, 1) = 1.0d0 - alph1
    Pinv(2, 2) = alph2 * u_moy
    Pinv(2, 3) = -alph2
    Pinv(3, 1) = 0.5d0 * (alph1 - u_moy / c_moy)
    Pinv(3, 2) = -0.5d0 * (alph2 * u_moy - 1.0d0 / c_moy)
    Pinv(3, 3) = alph2 / 2.0d0

    ! Matrice P
    P_tilde(1, 1) = 1.0d0
	P_tilde(1, 2) = 1.0d0
    P_tilde(1, 3) = 1.0d0
    P_tilde(2, 1) = u_moy - c_moy
    P_tilde(2, 2) = u_moy
    P_tilde(2, 3) = u_moy + c_moy
    P_tilde(3, 1) = H_moy - c_moy * u_moy
    P_tilde(3, 2) = 0.5d0 * u_moy**2
    P_tilde(3, 3) = H_moy + c_moy * u_moy
	
	! Initialisation de la matrice Delta
	Delta_tilde = 0.0d0

	! Remplissage de la matrice Delta
	select case(i_cor)
		case(0)
			! Sans correction
			! Matrice Delta (valeurs propres absolues)
			Delta_tilde(1, 1) = dabs(u_moy - c_moy)
			Delta_tilde(2, 2) = dabs(u_moy)
			Delta_tilde(3, 3) = dabs(u_moy + c_moy)
		case(1)
			! Entropie seule
			delta = delta_star*(dabs(u_moy)+c_moy)
			lambda_1 = dabs(u_moy - c_moy) ! vitesse acoustique 1 
			lambda_3 = dabs(u_moy + c_moy) ! vitesse acoustique 2 
			if(lambda_1 <= delta) then
				psi_1 = 0.5d0*(lambda_1**2+delta**2)/delta
			else 
				psi_1 = lambda_1
			endif
			if(lambda_3 <= delta) then
				psi_3 = 0.5d0*(lambda_3**2+delta**2)/delta
			else 
				psi_3 = lambda_3
			endif
			! Matrice Delta 
			Delta_tilde(1, 1) = psi_1
			Delta_tilde(2, 2) = dabs(u_moy)
			Delta_tilde(3, 3) = psi_3
		case(2)
			! Complete
			delta = delta_star*(dabs(u_moy)+c_moy)
			lambda_1 = dabs(u_moy - c_moy) ! vitesse acoustique 1 
			lambda_2  = dabs(u_moy)		   ! vitesse materielle 
			lambda_3 = dabs(u_moy + c_moy) ! vitesse acoustique 2
			if(lambda_1 <= delta) then
				psi_1 = 0.5d0*(lambda_1**2+delta**2)/delta
			else 
				psi_1 = lambda_1
			endif
			if(lambda_2 <= delta) then
				psi_2 = 0.5d0*(lambda_2**2+delta**2)/delta
			else 
				psi_2 = lambda_2
			endif
			if(lambda_3 <= delta) then
				psi_3 = 0.5d0*(lambda_3**2+delta**2)/delta
			else 
				psi_3 = lambda_3
			endif
			! Matrice Delta 
			Delta_tilde(1, 1) = psi_1
			Delta_tilde(2, 2) = psi_2
			Delta_tilde(3, 3) = psi_3
	end select

    ! Calcul de la matrice |A_tilde|
	A_tilde = matmul(P_tilde, Delta_tilde)
    A_tilde = matmul(A_tilde, Pinv)
	
end subroutine Construct_matrix
 
 ! -----------------------------
 
 subroutine HLL_flux(w_L,w_R,F_L,F_R,F_num)
!=====================================================================
!  Compute Numerical fluxes with HLL SCHEME
!=====================================================================
! Parametres globaux 
use Module_parametres, only : N, gamma

! Déclarations
	implicit none
	double precision,dimension(:,:),intent(in) 		:: w_R,w_L					! États conservés (rho, rho*u, rho*E) a doite et a gauche de l'interface
	double precision,dimension(:,:),intent(in)		:: F_L,F_R					! Flux
	double precision, dimension(N+1)				:: rho_L,u_L,P_L,E_L		! Variables primitives a gauche
	double precision, dimension(N+1)				:: rho_R,u_R,P_R,E_R		! Variables primitives a droite
	double precision, dimension(N+1)				:: H_R,H_L					! Enthalpie totale a gauche et droite
	integer 										:: j						! Entier pour boucle
	double precision 								:: S_R,S_L,a_R,a_L,eps		! Double pour calcul de HLL
	double precision,dimension(3,N+1)			    :: F_hll					! Local vector for computing HLL flux 
    double precision,dimension(3)					:: wdif,num					! Vecteurs de 3 composants
	double precision,dimension(:,:), intent(inout) 	:: F_num 					! Flux num à retourner
! =======================================================================================================================================================
! Fin des declarations 	
	
	eps = 1.0d-8           ! tol for denom
	F_hll = 0.0d0		   ! initialisation

	! Calcul du left state
	call Convert_conservative_to_primal_v2(w_L,rho_L,u_L,P_L,E_L)
	
	! Calcul du right state
	call Convert_conservative_to_primal_v2(w_R,rho_R,u_R,P_R,E_R)
	
	! Calcul de l'enthalpie totale a gauche
	H_L = (gamma/(gamma-1))*(P_L/rho_L)+0.5d0*u_L**2
	
	! Calcul de l'enthalpie totale a droite
	H_R = (gamma/(gamma-1))*(P_R/rho_R)+0.5d0*u_R**2
	
	! Compute S_L and S_R
	do j = 1, N + 1 			            ! Boucle sur les interfaces
		a_L = dsqrt(gamma*P_L(j)/ rho_L(j)) ! Vitesse du son 
		a_R = dsqrt(gamma*P_R(j)/ rho_R(j)) ! Vitesse du son 
		S_L = min(u_L(j) - a_L, u_R(j) - a_R)
		S_R = max(u_L(j) + a_L, u_R(j) + a_R)
	 
		! State difference
		wdif(:) = w_R(:,j) - w_L(:,j)
		
		! Compute the numerical flux
		num(:) = S_R*F_L(:,j)-S_L*F_R(:,j)+S_R*S_L*wdif(:)
		F_hll(:,j) = num(:) / max(dabs(S_R - S_L), eps) 
		if(0.0d0 <= S_L) then
			F_num(:,j) = F_L(:,j)
		elseif(S_L < 0.0d0 .and. 0.0d0 < S_R) then
			F_num(:,j) = F_hll(:,j)
		elseif(0.0d0 >= S_R) then
			F_num(:,j) = F_R(:,j)
		endif
	enddo

 end subroutine HLL_flux
 
 
 ! -----------------------------
 
  subroutine HLLR_flux(w_L,w_R,F_L,F_R,F_num)
!=====================================================================
!  Compute Numerical fluxes with HLLR SCHEME
!=====================================================================
! Parametres globaux 
use Module_parametres, only : N, gamma

  ! Déclaration des variables
	implicit none
	double precision,dimension(:,:),intent(in) 		:: w_R,w_L					! États conservés (rho, rho*u, rho*E) a doite et a gauche de l'interface
	double precision,dimension(:,:),intent(in)		:: F_L,F_R					! Flux
	double precision, dimension(N+1)				:: rho_L,u_L,P_L,E_L		! Variables primitives a gauche
	double precision, dimension(N+1)				:: rho_R,u_R,P_R,E_R		! Variables primitives a droite
	double precision, dimension(N+1)				:: H_R,H_L					! Enthalpie totale a gauche et droite
	double precision 								:: u_moy, H_moy, c_moy,a	! Valeurs moyennes de Roe
	integer 										:: j						! Entier pour boucle
	double precision 								:: S_R,S_L,a_R,a_L,eps		! Double pour calcul de HLL
	double precision,dimension(3,N+1)			    :: F_hll					! Local vector for computing HLL flux 
    double precision,dimension(3)					:: wdif,num					! Vecteurs de 3 composants
	double precision,dimension(:,:), intent(inout) 	:: F_num 					! Flux num à retourner
! =======================================================================================================================================================
! Fin des declarations 	
	
	eps = 1.0d-8           ! tol for denom
	F_hll = 0.0d0		   ! initialisation

	! Calcul du left state
	call Convert_conservative_to_primal_v2(w_L,rho_L,u_L,P_L,E_L)
	
	! Calcul du right state
	call Convert_conservative_to_primal_v2(w_R,rho_R,u_R,P_R,E_R)
	
	! Calcul de l'enthalpie totale a gauche
	H_L = (gamma/(gamma-1))*(P_L/rho_L)+0.5d0*u_L**2
	
	! Calcul de l'enthalpie totale a droite
	H_R = (gamma/(gamma-1))*(P_R/rho_R)+0.5d0*u_R**2
	
	! Compute S_L and S_R
	do j = 1, N + 1 													! Boucle sur les interfaces		            							
		a = dsqrt(rho_R(j)/ rho_L(j))									! rapport des masses volumiques
		u_moy = (a*u_R(j)+u_L(j))/(1.0d0+a)								! Moyenne de Roe de la vitesse matérielle
		H_moy = (a*H_R(j)+H_L(j))/(1.0d0+a)								! Moyenne de Roe de l'enthalpie	
		c_moy = dsqrt((gamma - 1.0d0) * (H_moy - 0.5d0 * u_moy**2))		! Moyenne de Roe de la vitesse de l'onde acoustique
		S_L = u_moy - c_moy
		S_R = u_moy + c_moy
	 
		! State difference
		wdif(:) = w_R(:,j) - w_L(:,j)
	 
		! Compute the numerical flux
		num(:) = S_R*F_L(:,j)-S_L*F_R(:,j)+S_R*S_L*wdif(:)
		F_hll(:,j) = num(:) / max(dabs(S_R - S_L), eps) 
		if(0.0d0 <= S_L) then
			F_num(:,j) = F_L(:,j)
		elseif(S_L < 0.0d0 .and. 0.0d0 < S_R) then
			F_num(:,j) = F_hll(:,j)
		elseif(0.0d0 >= S_R) then
			F_num(:,j) = F_R(:,j)
		endif
	enddo
	

 end subroutine HLLR_flux
 
 
 ! -----------------------------
 
 
  subroutine HLLE_flux(w_L,w_R,F_L,F_R,F_num)
!=====================================================================
!  Compute Numerical fluxes with HLLE SCHEME
!=====================================================================
! Parametres globaux 
use Module_parametres, only : N, gamma

  ! Déclaration des variables
	implicit none
	double precision,dimension(:,:),intent(in) 		:: w_R,w_L					! États conservés (rho, rho*u, rho*E) a doite et a gauche de l'interface
	double precision,dimension(:,:),intent(in)		:: F_L,F_R					! Flux
	double precision, dimension(N+1)				:: rho_L,u_L,P_L,E_L		! Variables primitives a gauche
	double precision, dimension(N+1)				:: rho_R,u_R,P_R,E_R		! Variables primitives a droite
	double precision, dimension(N+1)				:: H_R,H_L					! Enthalpie totale a gauche et droite
	double precision 								:: u_moy, H_moy, c_moy,a	! Valeurs moyennes de Roe
	double precision 								:: eps,eta_1,eta_2,d,d_b    ! Double pour calcul 
	integer 										:: j						! Entier pour boucle
	double precision 								:: S_R,S_L,a_R,a_L			! Double pour calcul de HLL
	double precision,dimension(3,N+1)			    :: F_hll					! Local vector for computing HLL flux 
    double precision,dimension(3)					:: wdif,num					! Vecteurs de 3 composants
	double precision,dimension(:,:), intent(inout) 	:: F_num 					! Flux num à retourner
! =======================================================================================================================================================
! Fin des declarations 	
	
	eps = 1.0d-8           ! tol for denom
	F_hll = 0.0d0		   ! initialisation

	! Calcul du left state
	call Convert_conservative_to_primal_v2(w_L,rho_L,u_L,P_L,E_L)
	
	! Calcul du right state
	call Convert_conservative_to_primal_v2(w_R,rho_R,u_R,P_R,E_R)
	
	! Calcul de l'enthalpie totale a gauche
	H_L = (gamma/(gamma-1))*(P_L/rho_L)+0.5d0*u_L**2
	
	! Calcul de l'enthalpie totale a droite
	H_R = (gamma/(gamma-1))*(P_R/rho_R)+0.5d0*u_R**2
	
	do j = 1, N + 1 			            ! Boucle sur les interfaces
		! Compute velocities 
		a = dsqrt(rho_R(j)/ rho_L(j))									
		u_moy = (a*u_R(j)+u_L(j))/(1.0d0+a)	! Moyenne de Roe
		a_L = dsqrt(gamma*P_L(j)/ rho_L(j)) ! Vitesse du son 
		a_R = dsqrt(gamma*P_R(j)/ rho_R(j)) ! Vitesse du son
		
		! Compute d_barre : d_b
		eta_1 = dsqrt(rho_L(j))+dsqrt(rho_R(j))
		eta_2 = 0.5d0*((dsqrt(rho_L(j))*dsqrt(rho_R(j)))/eta_1**2)
		d = dsqrt(rho_L(j))*a_L**2 + dsqrt(rho_R(j))*a_R**2
		d_b = d/eta_1  + eta_2*(u_R(j) - u_L(j))**2
		
		! Compute S_L and S_R
		S_L = u_moy - dsqrt(d_b)
		S_R = u_moy + dsqrt(d_b)
		
		! State difference
		wdif(:) = w_R(:,j) - w_L(:,j) 
		
		! Compute the numerical flux
	    num(:) = S_R*F_L(:,j)-S_L*F_R(:,j)+S_R*S_L*wdif(:)
		F_hll(:,j) = num(:) / max(dabs(S_R - S_L), eps) 
		if(0.0d0 <= S_L) then
			F_num(:,j) = F_L(:,j)
		elseif(S_L < 0.0d0 .and. 0.0d0 < S_R) then
			F_num(:,j) = F_hll(:,j)
		elseif(0.0d0 >= S_R) then
			F_num(:,j) = F_R(:,j)
		endif
	enddo
	
 end subroutine HLLE_flux
 
 
 ! -----------------------------

subroutine HLLC_ANRS_flux(w_L,w_R,F_L,F_R,F_num)
! ===========================================================================================
! Compute the HLLC-ANRS flux
! ===========================================================================================
! Parametres globaux 
use Module_parametres, only : N, gamma

! Déclarations
 implicit none
	double precision,dimension(:,:),intent(in) 		:: w_R,w_L								 				! États conservés (rho, rho*u, rho*E) a doite et a gauche de l'interface
	double precision,dimension(:,:),intent(in)		:: F_L,F_R								 				! Flux
	double precision, dimension(N+1)				:: rho_L,u_L,P_L,E_L			 		 				! Variables primitives a gauche
	double precision, dimension(N+1)				:: rho_R,u_R,P_R,E_R			 		 				! Variables primitives a droite
	double precision, dimension(N+1)				:: H_R,H_L								 				! Enthalpie totale a gauche et droite
	integer 										:: j									 				! Entier pour boucle
	double precision 								:: S_R,S_L,a_R,a_L						 				! Double pour calcul de HLL
	double precision 								:: rho_c,a_c,p_star,u_star,q_L,q_R,rho_starR			! Double pour calcul de HLLC
	double precision 								:: var,var1,var2,var3,var4,var5,var6,var7,var8			! Double pour calcul de Q_star
	double precision 								:: p_max,p_min,Q,z,PLR,p_0,AA_L,AA_R,B_L,B_R,g_L,g_R	! Double pour calcul de p_star
	double precision 								:: den1,s1,s2,S_star,E_star,rho_starL				    ! Double pour denominateur
    double precision,dimension(3)					:: wdif,num,Q_starL,Q_starR				 				! Vecteurs de 3 composants
	double precision,dimension(:,:), intent(inout) 	:: F_num 								 				! Flux num à retourner
! ===========================================================================================
! Fin des declarations 	
	
	! Calcul du left state
	call Convert_conservative_to_primal_v2(w_L,rho_L,u_L,P_L,E_L)
	
	! Calcul du right state
	call Convert_conservative_to_primal_v2(w_R,rho_R,u_R,P_R,E_R)
	
	! Calcul de l'enthalpie totale a gauche
	H_L = (gamma/(gamma-1))*(P_L/rho_L)+0.5d0*u_L**2
	
	! Calcul de l'enthalpie totale a droite
	H_R = (gamma/(gamma-1))*(P_R/rho_R)+0.5d0*u_R**2
	
	! Boucle sur les interfaces
	do j = 1, N + 1
      
		! vitesses du son
		a_L = dsqrt(gamma*P_L(j)/ rho_L(j)) 
		a_R = dsqrt(gamma*P_R(j)/ rho_R(j)) 
		rho_c = 0.5d0*(rho_L(j)+rho_R(j))   ! rho chapeau
		a_c = 0.5d0*(a_L + a_R)				! a chapeau
		
		! Compute p_star (PVRS)
		p_star = 0.5d0*(P_L(j)+P_R(j))+0.5d0*(u_L(j)-u_R(j))*rho_c*a_c
		p_star = max(0.0d0,p_star)
		
		! Choix du solver PVRS ou TRRS ou TSRS pour compute p_star and u_star
		p_max = max(P_L(j),P_R(j))
		p_min = min(P_L(j),P_R(j))
		Q = p_max / p_min
		if (Q <= 2.0d0 .and. p_min <= p_star .and. p_star <= p_max) then
			! PVRS
			!write(*,*)'PVRS'
		else
			if(p_star <= p_min) then
				! TRRS : Two–Rarefaction Riemann solver (version optimisee)
				!write(*,*)'TRRS'
				z = (gamma-1)/(2.0d0*gamma)
				PLR = (P_L(j)/P_R(j))**z
				var5 = PLR*(u_L(j)/a_L)+u_R(j)/a_R + 2.0d0*(PLR-1.0d0)/(gamma-1.0d0)
				var6 = PLR/a_L + 1.0d0/a_R
				u_star = var5 / var6
				var7 = 1.0d0 + ((gamma-1.0d0)/(2.0d0*a_L))*(u_L(j)-u_star)
				var8 = 1.0d0 + ((gamma-1.0d0)/(2.0d0*a_R))*(u_star-u_R(j))
				p_star = 0.5d0*(P_L(j)*var7**(1.0d0/z) + P_R(j)*var8**(1.0d0/z))
			else 
				! TSRS : Two–Shock Riemann Solver
				!write(*,*)'TSRS'
				p_0 = max(0.0d0,p_star)
				AA_L = 2.0d0/((gamma+1.0d0)*rho_L(j))
				AA_R = 2.0d0/((gamma+1.0d0)*rho_R(j))
				B_L = ((gamma-1.0d0)/(gamma+1.0d0))*P_L(j)
				B_R = ((gamma-1.0d0)/(gamma+1.0d0))*P_R(j)
				g_L = dsqrt(AA_L/(p_0+B_L))
				g_R = dsqrt(AA_R/(p_0+B_R))
				var5 = g_L*P_L(j)+g_R*P_R(j)-(u_R(j)-u_L(j))
				var6 = g_L + g_R
				p_star = var5 / var6
			endif
		endif
		
		! Compute q_L
		if(p_star <= P_L(j)) then
			q_L = 1.0d0
		else
			var = 1+((gamma+1)/(2.0d0*gamma))*(p_star/P_L(j) - 1.0d0)
			q_L = dsqrt(var)
 		endif
		! Compute q_R
		if(p_star <= P_R(j)) then
			q_R = 1.0d0
		else
			var = 1+((gamma+1)/(2.0d0*gamma))*(p_star/P_R(j) - 1.0d0)
			q_R = dsqrt(var)
 		endif
		
		! Compute S_L, S_star and S_R and S_star 
		S_L = u_L(j)-a_L*q_L
		S_R = u_R(j)+a_R*q_R
	 
		! 4) calcul S_star 
		s1 = S_L - u_L(j)
		s2 = S_R - u_R(j)
		den1 = rho_L(j)*s1 - rho_R(j)*s2
	 
		if (dabs(den1) < 1.d-12) then
			S_star = 0.d0
		else
			S_star = ( P_R(j) - P_L(j) + rho_L(j)*s1*u_L(j) - rho_R(j)*s2*u_R(j) ) / den1
		endif
	 
		! Etats étoile conservatifs
		! gauche
		if (dabs(S_L - S_star) < 1.d-12) then
			Q_starL(:) = w_L(:,j)
		else
			rho_starL = rho_L(j)*(S_L - u_L(j)) / (S_L - S_star)
			u_star = S_star
			E_star = (S_L - u_L(j))*rho_L(j)*(w_L(3,j)/rho_L(j)+(S_star+P_L(j)/(rho_L(j)*(S_L-u_L(j))))*(S_star-u_L(j)))/(S_L-S_star)
			Q_starL(1) = rho_starL
			Q_starL(2) = rho_starL*u_star
			Q_starL(3) = E_star
		endif
		
		! droite
		if (dabs(S_R - S_star) < 1.d-12) then
			Q_starR(:) = w_R(:,j)
		else
			rho_starR = rho_R(j)*(S_R - u_R(j)) / (S_R - S_star)
			u_star = S_star
			E_star = (S_R-u_R(j))*rho_R(j)*(w_R(3,j)/rho_R(j)+(S_star+P_R(j)/(rho_R(j)*(S_R-u_R(j))))*(S_star-u_R(j)))/(S_R-S_star)
			Q_starR(1) = rho_starR
			Q_starR(2) = rho_starR*u_star
			Q_starR(3) = E_star
		endif
		
		! Flux HLLC final
		if (0.0d0 <= S_L) then
			F_num(:,j) = F_L(:,j)
		elseif (S_L <= 0.0d0 .and. 0.0d0 <= S_star) then
			F_num(:,j) = F_L(:,j) + S_L*(Q_starL(:) - w_L(:,j))
		elseif (S_star <= 0.0d0 .and. 0.0d0 <= S_R) then
			F_num(:,j) = F_R(:,j) + S_R*(Q_starR(:) - w_R(:,j))
		else
			F_num(:,j) = F_R(:,j)
		endif
		
	enddo

 end subroutine HLLC_ANRS_flux
 
 ! -----------------------------
 
 
  subroutine HLLC_flux(w_L,w_R,F_L,F_R,F_num)
!=====================================================================
!  Compute Numerical fluxes with HLLC SCHEME
!=====================================================================
! Parametres globaux 
use Module_parametres, only : N, gamma

! Déclarations
 implicit none
	double precision,dimension(:,:),intent(in) 		:: w_R,w_L								 				! États conservés (rho, rho*u, rho*E) a doite et a gauche de l'interface
	double precision,dimension(:,:),intent(in)		:: F_L,F_R								 				! Flux
	double precision, dimension(N+1)				:: rho_L,u_L,P_L,E_L			 		 				! Variables primitives a gauche
	double precision, dimension(N+1)				:: rho_R,u_R,P_R,E_R			 		 				! Variables primitives a droite
	double precision, dimension(N+1)				:: H_R,H_L								 				! Enthalpie totale a gauche et droite
	 double precision 								:: u_moy, H_moy, c_moy,a							    ! Valeurs moyennes de Roe
	integer 										:: j									 				! Entier pour boucle
	double precision 								:: S_R,S_L,a_R,a_L						 				! Double pour calcul de HLL
	double precision 								:: u_star,rho_starR										! Double pour calcul de HLLC
	double precision 								:: den1,s1,s2,S_star,E_star,rho_starL				    ! Double pour denominateur
    double precision,dimension(3)					:: wdif,num,Q_starL,Q_starR				 				! Vecteurs de 3 composants
	double precision,dimension(:,:), intent(inout) 	:: F_num 								 				! Flux num à retourner
! =======================================================================================================================================================
! Fin des declarations 	

	! Calcul du left state
	call Convert_conservative_to_primal_v2(w_L,rho_L,u_L,P_L,E_L)
	
	! Calcul du right state
	call Convert_conservative_to_primal_v2(w_R,rho_R,u_R,P_R,E_R)
	
	! Calcul de l'enthalpie totale a gauche
	H_L = (gamma/(gamma-1))*(P_L/rho_L)+0.5d0*u_L**2
	
	! Calcul de l'enthalpie totale a droite
	H_R = (gamma/(gamma-1))*(P_R/rho_R)+0.5d0*u_R**2

	! Boucle sur les interfaces
	do j = 1, N+1
     
		! Vitesses du son
		a_L = dsqrt(gamma*P_L(j)/rho_L(j))
		a_R = dsqrt(gamma*P_R(j)/rho_R(j))
	 
		! Vitesses du son de Roe
		a = dsqrt(rho_R(j)/ rho_L(j))									! rapport des masses volumiques
		u_moy = (a*u_R(j)+u_L(j))/(1.0d0+a)								! Moyenne de Roe de la vitesse matérielle
		H_moy = (a*H_R(j)+H_L(j))/(1.0d0+a)								! Moyenne de Roe de l'enthalpie	
		c_moy = dsqrt((gamma - 1.0d0) * (H_moy - 0.5d0 * u_moy**2))		! Moyenne de Roe de la vitesse de l'onde acoustique
		
		! Vitesses d’onde HLLC 
		S_L = min(u_L(j) - a_L, u_moy - c_moy)
		S_R = max(u_R(j) + a_R, u_moy + c_moy)
     
		! Calcul S_star 
		s1 = S_L - u_L(j)
		s2 = S_R - u_R(j)
		den1 = rho_L(j)*s1 - rho_R(j)*s2
     
		if (abs(den1) < 1.d-12) then
			S_star = 0.d0
		else
			S_star = ( P_R(j) - P_L(j) + rho_L(j)*s1*u_L(j) - rho_R(j)*s2*u_R(j) ) / den1
		endif
    
		! Etats étoile conservatifs
		! Gauche
		if (abs(S_L - S_star) < 1.d-12) then
			Q_starL(:) = w_L(:,j)
		else
			rho_starL = rho_L(j)*(S_L - u_L(j)) / (S_L - S_star)
			u_star = S_star
			E_star = (S_L - u_L(j))*rho_L(j)*(w_L(3,j)/rho_L(j)+(S_star+P_L(j)/(rho_L(j)*(S_L-u_L(j))))*(S_star-u_L(j)))/(S_L-S_star)
			Q_starL(1) = rho_starL
			Q_starL(2) = rho_starL*u_star
			Q_starL(3) = E_star
		endif
     
		! Droite
		if (abs(S_R - S_star) < 1.d-12) then
			Q_starR(:) = w_R(:,j)
		else
			rho_starR = rho_R(j)*(S_R - u_R(j)) / (S_R - S_star)
			u_star = S_star
			E_star = (S_R-u_R(j))*rho_R(j)*(w_R(3,j)/rho_R(j)+(S_star+P_R(j)/(rho_R(j)*(S_R-u_R(j))))*(S_star-u_R(j)))/(S_R-S_star)
			Q_starR(1) = rho_starR
			Q_starR(2) = rho_starR*u_star
			Q_starR(3) = E_star
		endif
		
		! Flux HLLC final
		if (0.0d0 <= S_L) then
			F_num(:,j) = F_L(:,j)
		elseif (S_L <= 0.0d0 .and. 0.0d0 <= S_star) then
			F_num(:,j) = F_L(:,j) + S_L*(Q_starL(:) - w_L(:,j))
		elseif (S_star <= 0.0d0 .and. 0.0d0 <= S_R) then
			F_num(:,j) = F_R(:,j) + S_R*(Q_starR(:) - w_R(:,j))
		else
			F_num(:,j) = F_R(:,j)
		endif
		
	enddo

 end subroutine HLLC_flux
 
 ! -----------------------------

end module Module_riemann