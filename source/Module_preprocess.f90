module Module_preprocess

contains

! ----------------------------------------------------------------------
! METHODES
! ----------------------------------------------------------------------

subroutine Condition_initiale(x_i,An,w)
! ===========================================================================================
! Calcul des conditions initiales avec Pi et Ti 
! ===========================================================================================
! Parametres globaux 
use Module_parametres, only : N, gamma, Pt0, Tt0, r, L, i_BC, o_BC
! Declarations	  
  implicit none
  double precision,dimension(-1:),intent(in)    :: x_i  								   ! Variables d'entree
  double precision,dimension(-1:),intent(in)    :: An 									   ! Variables d'entree
  double precision                              :: rho_loc, T_loc, P_loc, u_loc, E_tot_loc ! Variables locales 
  double precision                              :: rapport, A_star,  M_loc, x_c       	   ! Variables locales 
  integer                                       :: i 									   ! Variables locales 
  double precision,dimension(:,-1:),intent(out) :: w 									   ! Variable de sortie
! ===========================================================================================
! Fin des declarations 
	
A_star = minval(An(0:N))  ! section au col

do i = -1, N+2
	! Calcul du rapport A/A*
    rapport = An(i) / A_star
    
    select case(i_BC)
        case(1)  ! subso/superso classique
            if (x_i(i) <= x_c) then
                call Newton_aire_mach(rapport, 0.5d0, M_loc)
            else
                call Newton_aire_mach(rapport, 2.0d0, M_loc)
            endif
        case(2)  ! superso pur
            call Newton_aire_mach(rapport, 2.0d0, M_loc)
    end select
    
    ! --- Relations Isentropiques ---
    T_loc     = Tt0 / (1.0d0 + 0.5d0*(gamma-1.0d0)*M_loc**2)
    P_loc     = Pt0 * (T_loc/Tt0)**(gamma/(gamma-1.0d0))
    rho_loc   = P_loc / (r*T_loc)
    u_loc     = M_loc * dsqrt(gamma*r*T_loc)
    E_tot_loc = P_loc/(rho_loc*(gamma-1.0d0)) + 0.5d0*u_loc**2
    
	! --- Remplissage des variables conservatives ---
    w(1,i) = rho_loc * An(i)
    w(2,i) = rho_loc * u_loc * An(i)
    w(3,i) = rho_loc * E_tot_loc * An(i)
enddo


 end subroutine Condition_initiale
 
  ! -----------------------------
 
 subroutine Mesh(dx,x_i)
! ===========================================================================================
! Creation du maillage VF cell-centered
! ===========================================================================================
! Parametres globaux 
use Module_parametres, only : N
  implicit none
  double precision,intent(in) 					  :: dx			 ! pas de discretisation
  integer 			 							  :: i 			 ! entier pour calcul de boucles
  double precision,dimension(-1:),intent(out)	  :: x_i		 ! maillage
! ===========================================================================================
! Fin des declarations  
  
  ! Mesh creation
	do i = -1, N+2
		x_i(i) = (i - 0.5D0) * dx		! position de la cell
	enddo

 end subroutine Mesh
 
 ! -----------------------------
 
subroutine Nozzle_geom(dx,x_i,An,Anf)
! ===========================================================================================
! Discretisation de la geometrie de la tuyere
! ===========================================================================================
! Parametres globaux 
use Module_parametres, only : N, i_geo, L
  implicit none
  double precision,intent(in) 					  :: dx			 ! pas de discretisation
  double precision,dimension(-1:),intent(in)	  :: x_i		 ! maillage
  double precision						     	  :: x_f	 	 ! maillage interface
  integer 			 							  :: i 			 ! entier pour calcul de boucles
  double precision, dimension(-1:), intent(out)   :: An   		 ! geom de la tuyere au cells
  double precision, dimension(0:), intent(out)    :: Anf  		 ! geom de la tuyere au Interfaces 
! ===========================================================================================
! Fin des declarations 

! Validate input
  if (i_geo < 1 .or. i_geo > 5) then
     print *, "Error: Invalid i_geo value. Must be between 1 and 7."
     stop
  end if

! Nozzle geometry creation
select case(i_geo)
	case(1)	
		! divergent 
		! --- Aires aux cellules (An) ---
		do i = -1, N+2
			An(i) = 1.398d0 + 0.347d0*dtanh((8.0d0/L)*(x_i(i)-0.5d0*L))
		enddo
		! --- Aires aux INTERFACES (Anf) ---
		do i = 0, N
			x_f = dble(i) * dx  ! Position de l'interface
			Anf(i) = 1.398d0 + 0.347d0*dtanh((8.0d0/L)*(x_f-0.5d0*L))
		enddo
	case(2)
		! Tuyère de Laval 
		! --- Aires aux cellules (An) ---
		do i = -1, N+2
			if (x_i(i) <= 0.5d0*L) then
				An(i) = 1.0d0 + 6.0d0*(x_i(i) - 0.5d0*L)**2
			else
				An(i) = 1.0d0 + 2.0d0*(x_i(i) - 0.5d0*L)**2
			endif
		enddo
		! --- Aires aux Interfaces (Anf) ---
		do i = 0, N
			x_f = dble(i) * dx
			if (x_f <= 0.5d0*L) then
				Anf(i) = 1.0d0 + 6.0d0*(x_f - 0.5d0*L)**2
			else
				Anf(i) = 1.0d0 + 2.0d0*(x_f - 0.5d0*L)**2
			endif
		enddo
		
	case(3)
		! Tuyère CV/DV
		! --- Aires aux cellules (An) ---
		do i = -1, N+2
			if (x_i(i) <= 0.25d0*L) then
				An(i) = 1.0d0 + 0.1d0*(x_i(i) - 0.25d0*L)**2
			else
				An(i) = 1.0d0 + 2.0d-2*(x_i(i) - 0.25d0*L)**2
			endif
		enddo
		! --- Aires aux Interfaces (Anf) ---
		do i = 0, N
			x_f = dble(i) * dx
			if (x_f <= 0.5d0*L) then
				Anf(i) = 1.0d0 + 0.1d0*(x_f - 0.5d0*L)**2
			else
				Anf(i) = 1.0d0 + 2.0d-2*(x_f - 0.5d0*L)**2
			endif
		enddo
	
	case(4)
		! Soufflerie a deux col 2eme ferme
		! --- Aires aux cellules (An) ---
		do i = -1, N+2
			if (x_i(i)<= 1.0d0) then
				An(i) = 0.093d0 + 3.628d0*(x_i(i)-0.5d0)**2
			else if (x_i(i)<= 2.0d0) then
				An(i) = 1.0d0 
			else
				An(i) = 1.0d0 + 4.0d0*(x_i(i) -2.0d0)**2
			endif
		enddo
		! --- Aires aux Interfaces (Anf) ---
		do i = 0, N
			x_f = dble(i) * dx
			if (x_f<= 1.0d0) then
				Anf(i) = 0.093d0 + 3.628d0*(x_f-0.5d0)**2
			else if (x_f<= 2.0d0) then
				Anf(i) = 1.0d0 
			else
				Anf(i) = 1.0d0 + 4.0d0*(x_f-2.0d0)**2
			endif
		enddo
		
	case(6)
		! Parabolic nozzle
		! --- Aires aux cellules (An) ---
		do i = -1, N+2
			An(i) = 1.0d0 + 2.2d0*(x_i(i) - 1.5d0)**2
		enddo
		!Acol = minval(An(1:N))  ! Référence unique pour tout le domaine
		! --- Aires aux Interfaces (Anf) ---
		do i = 0, N
			x_f = dble(i) * dx
			Anf(i) = 1.0d0 + 2.2d0*(x_f - 1.5d0)**2
		enddo
	
	
end select
  
 end subroutine Nozzle_geom
 
! -----------------------------
  
  
 subroutine Newton_aire_mach(rapport, M_init, M_sol)
! ===========================================================================================
! Résous la relation aire-Mach par Newton pour chaque maille :
! f(Mach) = A/A*
! ===========================================================================================
! Parametres globaux 
use Module_parametres, only : gamma
    implicit none
    double precision, intent(in)  :: rapport, M_init  ! A(x)/A*
    double precision 			  :: M, f, f_prime, expo
    double precision 			  :: tol, error
    integer          			  :: iter, max_iter
	double precision, intent(out) :: M_sol
! ===========================================================================================
! Fin des declarations    
    tol      = 1.0d-8
    max_iter = 50
    M        = M_init
    expo     = (gamma+1.0d0)/(2.0d0*(gamma-1.0d0))
    
	! Boucle du Newton
    do iter = 1, max_iter
		! Calcul de f
        f = rapport - (1.0d0/M) * &
            ((2.0d0/(gamma+1.0d0)) * &
            (1.0d0 + 0.5d0*(gamma-1.0d0)*M**2))**expo
		! Calcul de f prime
        f_prime = (1.0d0/M**2) * &
            ((2.0d0/(gamma+1.0d0)) * &
            (1.0d0 + 0.5d0*(gamma-1.0d0)*M**2))**expo &
            - (1.0d0/M) * expo * &
            ((2.0d0/(gamma+1.0d0)) * &
            (1.0d0 + 0.5d0*(gamma-1.0d0)*M**2))**(expo-1.0d0) &
            * (2.0d0/(gamma+1.0d0)) * (gamma-1.0d0)*M
		! Calcul de l'erreur	
        error = dabs(f)
        if (error < tol) exit
		! Calcul du mach
        M = M - f/f_prime
        M = max(M, 1.0d-3)  ! sécurité M > 0
    enddo
    
    M_sol = M
end subroutine Newton_aire_mach
 

end module Module_preprocess