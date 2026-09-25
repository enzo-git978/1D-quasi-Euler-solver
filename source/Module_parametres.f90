module Module_parametres

implicit none

    ! --- Choix des schémas ---
	integer	:: i_init ! choix de la cdt initiale voulue 
    integer :: i_sc   ! Schéma spatial (1: Roe, 2: HLL, 3: HLLC)
    integer :: i_ord  ! Ordre spatial (1: Ordre 1, 2: MUSCL)
    integer :: i_fct  ! Limiteur de pente (1: Minmod, 2: Superbee, etc.)
    integer :: i_BC   ! Conditions aux limites d'entree (1: entree subso, 2: Sortie superso)
	integer :: o_BC   ! Conditions aux limites de sortie (1: sortie subso, 2: Sortie superso)
    integer :: i_cor  ! Correction d'entropie (0: Non, 1: Oui)
	integer :: i_geo  ! Choix de la geometrie de la tuyere
	integer :: i_gaz  ! Choix du gaz
	integer :: i_impl ! Choix de l'implicitation
	
	! --- Parametres pour schema ---
	double precision :: delta_star	! Coefficient pour correction entropique de Harten
	double precision :: beta     	! Coefficient pour limiter Chakravarthy
    double precision :: b,phi		! Parametre de compression, Parametre ordre de reconstruction
	
	! --- Parametres du gaz considere ---
	double precision :: gamma  		  ! Constante Adiabatique = CP/Cv
	double precision :: Cp			  ! Capacite calorifique a pression constante
	double precision :: r			  ! Cp-Cv J/(kg.K)
	
	! --- Parametres physiques ---
	double precision :: Pt0			  ! Pression generatrice
	double precision :: Tt0			  ! Température generatrice
	double precision :: Pext		  ! Pression de sortie
    double precision :: Cf			  ! Coefficient de frottement
    double precision :: h 			  ! Coefficient de convection
    double precision :: Tp 			  ! Température paroi
	
	! --- Parametres numeriques ---
	integer	         :: N       ! Nombres de mailles
	integer          :: maxiter ! Nb max d iterations
	double precision :: tol_rho ! Precision des residus de rho
	double precision :: tol_u 	! Precision des residus de u
	double precision :: tol_e   ! Precision des residus de e
	double precision :: L   	! Longueur du domaine
    double precision :: CFL 	! Condition CFL
	
	
end module Module_parametres
