program main
! ========================================================================================
! Programme de resolution du systeme d'Euler constitu� des lois de conservation : 
! conservation de la masse, de la quantit� de mvt et de l'energie
! en 1D instationnaire : d_t w + d_x F(w) = 0
! ou w = {rho,rho u, rho E} et F = {rho u, rho u�+P, (rho E+P) u}
! Condition initiale : w(x,0) = {w_R si x< x_d w_L sinon}
! Discr�tisation spatiale choisit : n pts; x_i = (i-0.5)*h; i=1,n (cell-centered mesh)
!
! Domain cells :
!
!              |           |   u(i)    |           |
!              |  u(i-1)   |___________|           |
!              |___________|           |   u(i+1)  |
!              |           |           |___________|
!           ...|-----x-----|-----x-----|-----x-----|...
!              |    i-1    |     i     |    i+1    |
!              |-         +|-         +|-         +|
!            i-3/2       i-1/2       i+1/2       i+3/2

! Avec Reconstruction MUSCL : les variables d �tats ne sont plus cstes mais 
! �voluent de fa�on lin�aires dans chaque mailles
! ========================================================================================


!!!!!!!!!!!!!!!!!!!!!!
! to create executable : make 
! to execute prog 	   : ./solve
! to clean 			   : make clean
!!!!!!!!!!!!!!!!!!!!!!
	
! Modules utilises
use Module_solve
use Module_lecture
use Module_parametres										

! Declarations
  implicit none
  double precision 			  	 		     :: dx,t1,t2		! reels double precision
  double precision, dimension(:),allocatable :: rho, u, P, E	! Variables des lois de conservations du fluide 
! Fin des declarations 
! =========================================================================================================================
! Main
! =========================================================================================================================

! Lecture of menu file
call read_parameters(10)   ! reads variables 

! Computation of space step dx 
dx = L /(N-1)		    ! dx : Taille de la cell ; (N-1) : nombre de cells
	
! Allocation dynamique de memoire
allocate(rho(-1:N+2))   ! Masse volumique           [kg/m3]
allocate(P(-1:N+2))     ! Pression 	                [Pa]
allocate(u(-1:N+2))	    ! Vitesse 		            [m/s]
allocate(E(-1:N+2))	    ! Energie totale sp�cifique [J/kg]

! Resolution des equations d'Euler avec MUSCL pour un fluide parfait unidimensionnel compressible
call CPU_time(t1)
call Solve_Eqn_Euler_MUSCL(dx,rho,u,P,E)
call CPU_time(t2)
write(*,*)'Temps de resolution :',t2-t1 ! tps de calcul de la subroutine
	
! Trace avec GNUPLOT les graphes au temps final 
call system('gnuplot data_final/plotdata_nozzle.plt') 
	
! Ecris les parametres calcul a l'ecran
call display_parameters()  ! display on screen the variables selected 
	
! Desallocation dynamique de memoire
deallocate(rho, u, P, E)


!--------------fin
write(*,*)
write(*,*)'Execution avec succes !'
!===============
end program main
!===============