module Module_lecture

contains

! ----------------------------------------------------------------------
! METHODES
! ----------------------------------------------------------------------

subroutine read_parameters(nfich)
!=====================================================================
!  Fonction : Lecture du fichier menu pour les parametres de calcul
!=====================================================================
! Parametres globaux
use Module_parametres

! Déclarations
  implicit none
  integer,          intent(in)   :: nfich		    ! Numero du fichier pour ouverture
! =========================================================================================================================
! Fin des declarations
  
! Main
WRITE(*,*) 'Debut lecture du fichier menu :'
open(nfich,file='menu', form='formatted')
! Parametre physiques
read(nfich,'(17/)') ! skip header
read(nfich,*) i_gaz
select case(i_gaz)
	case(1)
		! Gaz frais
		gamma = 1.4d0
		Cp = 1004.64d0
	case(2)
		! Gaz moyennement chaud
		gamma = 1.35d0
		Cp = 1150d0
	case(3)
		! Gaz tres chaud
		gamma = 1.32d0
		Cp = 1150d0
end select
! Calcul de r pour le gaz choisi
r = Cp*(gamma-1.0d0)/gamma
read(nfich,*) ! Saute la ligne vide
read(nfich,*) ! Saute Description
read(nfich,*) Pt0
read(nfich,*) ! Saute la ligne vide
read(nfich,*) ! Saute Description
read(nfich,*) Tt0
read(nfich,*) ! Saute la ligne vide
read(nfich,*) ! Saute Description
read(nfich,*) Pext
read(nfich,*) ! Saute la ligne vide
read(nfich,*) ! Saute Description
read(nfich,*) Cf
read(nfich,*) ! Saute la ligne vide
read(nfich,*) ! Saute Description
read(nfich,*) h
read(nfich,*) ! Saute la ligne vide
read(nfich,*) ! Saute Description
read(nfich,*) Tp
read(nfich,*) ! Saute la ligne vide
read(nfich,*) ! Saute Description
read(nfich,*) ! Saute Description
read(nfich,*) N
read(nfich,*) ! Saute la ligne vide
read(nfich,*) ! Saute Description
read(nfich,*) L
read(nfich,*) ! Saute la ligne vide
read(nfich,*) ! Saute Description
read(nfich,*) CFL
read(nfich,*) ! Saute la ligne vide
read(nfich,*) ! Saute Description
read(nfich,*) maxiter
read(nfich,*) ! Saute la ligne vide
read(nfich,*) ! Saute Description
read(nfich,*) tol
read(nfich,*) ! Saute la ligne vide
read(nfich,*) ! Saute Description
read(nfich,*) i_geo
read(nfich,*) ! Saute la ligne vide
read(nfich,*) ! Saute Description
read(nfich,*) i_impl
read(nfich,*) ! Saute la ligne vide
read(nfich,*) ! Saute Description
read(nfich,*) ! Saute Description
! BC
read(nfich,*) i_BC
read(nfich,*) ! Saute la ligne vide
read(nfich,*) ! Saute Description
read(nfich,*) o_BC
read(nfich,*) ! Saute la ligne vide
read(nfich,*) ! Saute Description
read(nfich,*) ! Saute Description
! Simulation
read(nfich,*) i_sc
read(nfich,*) ! Saute la ligne vide
read(nfich,*) ! Saute Description
read(nfich,*) i_ord
read(nfich,*) ! Saute la ligne vide
read(nfich,*) ! Saute Description
read(nfich,*) i_cor
read(nfich,*) ! Saute la ligne vide
read(nfich,*) ! Saute Description
read(nfich,*) delta_star
read(nfich,*) ! Saute la ligne vide
read(nfich,*) ! Saute Description
read(nfich,*) phi
if(phi==3) then
	phi = 1.0d0/3.0d0
endif
read(nfich,*) ! Saute la ligne vide
read(nfich,*) ! Saute Description
read(nfich,*) b
read(nfich,*) ! Saute la ligne vide
read(nfich,*) ! Saute Description
read(nfich,*) i_fct
read(nfich,*) ! Saute la ligne vide
read(nfich,*) ! Saute Description
read(nfich,*) beta

close(nfich)

end subroutine read_parameters


! -----------------------------


subroutine display_parameters()

!=====================================================================
!  Fonction : Affichage des parametres de calcul selectionnes
!=====================================================================
! Parametres globaux
use Module_parametres
 ! =========================================================================================================================
! Fin des declarations 
  
! Main
write(*,*) '### PARAMETRES DU CALCUL :  ###'
write(*,*)'-- Physiques -- : '
write(*, '(A25, F8.3)')  'gamma = ', gamma
write(*, '(A25, F8.3)')  'Cp = ', Cp
write(*, '(A25, F8.1, A8)') 'Stagnation pressure = ', Pt0/1.0d5, ' bar'
write(*, '(A25, F8.1, A8)') 'Stagnation temperature = ', Tt0, ' Kelvin'
write(*, '(A25, F8.1, A8)') 'Exit pressure = ', Pext/1.0d5, ' bar'
if (Cf /= 0.0d0) then 
	write(*, '(A25)') 'Friction model activated ! '
	write(*, '(A25, F8.2, A8)') 'Cf = ', Cf, ' '
endif 
if (h /= 0.0d0) then 
	write(*, '(A25)') 'Heat model activated ! '
	write(*, '(A25, F8.2, A8)') 'h = ', h, ' W/(K.m2)'
	write(*, '(A25, F8.2, A8)') 'Tp = ', Tp, ' Kelvin'
endif 
write(*,*)'*****************************************'
write(*,*)'-- Simulation -- : '
write(*, '(A25, I10)')  'Nombre de points N = ', N
write(*, '(A25, F8.3)')  'Longueur de la tuyere = ', L
write(*, '(A25, F8.3)')  'cfl = ', cfl
write(*, '(A25, I10)')  'maxiter = ', maxiter
write(*, *)  'tol = ', tol


write(*,*)'*****************************************'
write(*,*)'-- Geometrie de la Tuyere -- : '
select case(i_geo)
	case(1)
	    write(*,*) 'Divergent'
	case(2)
		write(*,*) ' Tuyere de Laval'
	case(3)
		write(*,*) 'Tuyere CV/DV'
	case(4)
		write(*,*) 'Soufflerie double col (2eme col ferme)'
	case(5)
		write(*,*) 'Soufflerie double col (2eme col ouvert)'
end select

write(*,*)'*****************************************'
write(*,*)'-- BC -- : '
select case(i_BC)
	case(1)
	    write(*,*) 'Subsonic inlet'
	case(2)
		write(*,*) 'Supersonic inlet'
end select

select case(o_BC)
	case(1)
	    write(*,*) 'Subsonic outlet'
	case(2)
		write(*,*) 'Supersonic outlet'
end select

write(*,*)'*****************************************'
write(*,*)'Choice of the flux scheme :'
select case(i_sc)
	case(1)
		write(*,*) 'Roe'
	case(2)
	    write(*,*) 'HLL (Davis estimation)'
	case(3)
		write(*,*) 'HLL (Mean estimation)'
	case(4)
		write(*,*) 'HLLE (Einfeldt estimation)'
	case(5)
		write(*,*) 'HLLC-ANRS'
	case(6)
		write(*,*) 'HLLC robuste'
end select


select case(i_ord)
	case(1)
	    write(*,*) 'ordre 1 en espace'
	case(2)
		if(dabs(phi) - 1.0d0/3.0d0 <= 1.0d-6) then
			write(*,*) 'ordre 3 en espace'
		else
			write(*,*) 'ordre 2 en espace'
		endif
		write(*, '(A, F8.3)')  'PHI   = ', phi
		write(*, '(A, F8.3)')  'parametre de compression b = ', b
		select case(i_fct)
			case(1)
				write(*,*) 'Choix du limiter : Minmod'
			case(2)
				write(*,*) 'Choix du limiter : VanLeer'
			case(3)
				write(*,*) 'Choix du limiter : VanAlbada'
			case(4)
				write(*,*) 'Choix du limiter : Superbee'
			case(5)
				write(*,*) 'Choix du limiter : Chakravarthy'
				write(*, '(A, F8.3)')  'beta     = ', beta
		end select
end select

select case(i_cor)
	case(0)
	    write(*,*) 'Sans correction entropique'
		write(*, *)
	case(1)
		write(*,*) 'Correction entropique'
		write(*,'(A, F8.3)') 'Coefficient delta_star=',delta_star
		write(*, *)
	case(2)
		write(*,*) 'Correction entropique complete'
		write(*, '(A, F8.3)') 'Coefficient delta_star=',delta_star
		write(*, *)
end select


end subroutine display_parameters

end module Module_lecture