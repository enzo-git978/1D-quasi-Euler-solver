module Module_sortie

contains

! ----------------------------------------------------------------------
! ECRITURE VECTEURS ET MATRICES
! ----------------------------------------------------------------------

subroutine ecrit_vec(x,n) ! Ecrit le vecteur coeff par coeff
    implicit none
    real(kind=kind(0.d0)),dimension(:)			            :: x    
    integer                                    				:: i,n
    do i=1,n
        write(*,fmt = '(I4,1X,F7.3)')i,x(i)
    enddo
	write(*,*)
end subroutine ecrit_vec

! ----------------------------------------------------------------------

subroutine ecrit_mat(A,n) ! Ecrit la matrice à l'ecran sous un format 'reconnaissable'
    implicit none
    real(kind=kind(0.d0)),dimension(:,:)			     :: A    
    integer                                  			 :: i,j,n
    do i=1,3
        do j=1,n
            write(*,'(F8.4, "  ")', advance = 'no')A(i,j)
        enddo
      	write(*,*)
    enddo
	write(*,*)
end subroutine ecrit_mat

! ----------------------------------------------------------------------


subroutine save_results(Save_choice,dx,Tf,An,rho,u,P,E,Ma,s,Temp)
! ================================================
! write results in a file in the data folder
! ================================================
! Parametres globaux 
use Module_parametres, only : N, L

 ! Déclaration des variables
  implicit none
  double precision, dimension(-1:),intent(in)   :: rho, u, P, E						! Caractéristiques physique du fluide
  double precision, dimension(:),intent(in)     :: Ma, s, Temp					   	! Caractéristiques physique du fluide
  double precision,intent(in) 					:: dx,Tf							! parametres d'entrée de la subroutine
  character(len=*),intent(in) 				    :: Save_choice						! Pick either data_final or film
  double precision,dimension(-1:),intent(in)    :: An		 						! Geometrie de la tuyere au centre
  character(len=10) 							:: time_str							! char for time
  character(len=50) 							:: filename							! name of the file
  integer 			 							:: i,ierr 			 				! entier pour calcul de boucles
! =========================================================================================================================
! Fin des declarations 	  
	
	! Convert time into character
	write(time_str, '(F5.3)') Tf

	! Construct filename without spaces
	filename = trim(Save_choice) // '/results_at_t=' // trim(adjustl(time_str)) // 's.dat'
	
	! Verification 
	print *, "Save file : ", trim(filename)
    
    ! Open the file with the generated filename
    open(unit=10, file=trim(filename), status="unknown", form="formatted", iostat=ierr)
    if (ierr /= 0) then
        print*, 'Error opening file:', trim(filename)
        stop
    end if
	
	! writing
	do i = 1, N-1
		write(10, *) (i-0.5d0)*dx, An(i), rho(i), Ma(i), P(i), Temp(i), rho(i)*u(i)*An(i), s(i)
	end do
	write(10, *) L, An(N), rho(N), Ma(N), P(N), Temp(N), rho(N)*u(N)*An(N), s(N)
	close(10)
  
  end subroutine save_results
  
  ! ----------------------------------------------------------------------
  
subroutine create_gnuplot_script(Tf,Save_choice)
! =============================================================
! Creates a gnuplot script in the data folder
! This script plots the results for t=Tf in a multiplot layout
! =============================================================
  implicit none
  character(len=100)          :: gnuplot_file        ! Name of the gnuplot script
  integer                     :: iunit               ! For handling errors
  double precision,intent(in) :: Tf                  ! Simulation time
  character(len=*),intent(in) :: Save_choice		 ! Pick either data_final or film
  character(len=100) 		  :: output_png		     ! Name of output
	
	! name of the script
  gnuplot_file = trim(Save_choice) // '/plotdata_nozzle.plt'

  ! open file
  open(newunit=iunit, file=trim(gnuplot_file), status="replace")
  
  ! Temps final et fichier de donnees
  write(iunit, '(A,F6.3)') "t_final = ", Tf
  write(iunit, '(A,A,A)') "filename_final = sprintf('", trim(Save_choice), "/results_at_t=%.3fs.dat', t_final)"
  write(iunit, *) ""

  ! =============================================================
  ! PLOT 1 : Geometrie seule
  ! =============================================================
  write(iunit, '(A)') "set terminal pngcairo enhanced font 'arial,10' size 1200, 400"
  write(iunit, '(A,A,A)') "set output '", trim(Save_choice), "/Geometrie_tuyere.png'"
  write(iunit, *) ""
  write(iunit, '(A)') "set xlabel 'x (m)'"
  write(iunit, '(A)') "set ylabel 'A(x)'"
  write(iunit, '(A)') "set title 'Geometrie de la tuyere'"
  write(iunit, '(A)') "set grid"
  write(iunit, '(A,A,A)') "plot filename_final using 1:2 with lines lw 2 lc rgb 'black' title 'A(x)'"
  write(iunit, *) ""

  ! =============================================================
  ! PLOT 2 : 6 grandeurs physiques en multiplot 3x2
  ! =============================================================
  write(iunit, '(A)') "set terminal pngcairo enhanced font 'arial,10' size 1400, 900"
  write(iunit, '(A,A,A)') "set output '", trim(Save_choice), "/Solution_Numerique.png'"
  write(iunit, *) ""
  
  write(iunit, '(A)') "set multiplot layout 3,2 title 'Solution Numerique Tuyere'"
  write(iunit, *) ""

  ! --- 1. Densite ---
  write(iunit, '(A)') "set xlabel 'x (m)'"
  write(iunit, '(A)') "set ylabel 'rho'"
  write(iunit, '(A)') "set title 'Masse volumique [kg/m3]'"
  write(iunit, '(A)') "set grid"
  write(iunit, '(A)') "plot filename_final using 1:3 with lines lw 2 lc rgb 'blue' title sprintf('t=%.3fs', t_final)"
  write(iunit, *) ""

  ! --- 2. Mach ---
  write(iunit, '(A)') "set xlabel 'x (m)'"
  write(iunit, '(A)') "set ylabel 'Mach '"
  write(iunit, '(A)') "set title 'Nombre de Mach'"
  write(iunit, '(A)') "set grid"
  write(iunit, '(A)') "plot filename_final using 1:4 with lines lw 2 lc rgb 'red' title sprintf('t=%.3fs', t_final)"
  write(iunit, *) ""

  ! --- 3. Pression ---
  write(iunit, '(A)') "set xlabel 'x (m)'"
  write(iunit, '(A)') "set ylabel 'P'"
  write(iunit, '(A)') "set title 'Pression [Pa]'"
  write(iunit, '(A)') "set grid"
  write(iunit, '(A)') "plot filename_final using 1:5 with lines lw 2 lc rgb 'green' title sprintf('t=%.3fs', t_final)"
  write(iunit, *) ""

  ! --- 4. Temperature ---
  write(iunit, '(A)') "set xlabel 'x (m)'"
  write(iunit, '(A)') "set ylabel 'T'"
  write(iunit, '(A)') "set title 'Temperature [K]'"
  write(iunit, '(A)') "set grid"
  write(iunit, '(A)') "plot filename_final using 1:6 with lines lw 2 lc rgb 'orange' title sprintf('t=%.3fs', t_final)"
  write(iunit, *) ""

  ! --- 5. Debit ---
  write(iunit, '(A)') "set xlabel 'x (m)'"
  write(iunit, '(A)') "set ylabel 'rho_u_A'"
  write(iunit, '(A)') "set title 'Debit [kg/s]'"
  write(iunit, '(A)') "set grid"
  write(iunit, '(A)') "plot filename_final using 1:7 with lines lw 2 lc rgb 'purple' title sprintf('t=%.3fs', t_final)"
  write(iunit, *) ""

  ! --- 6. Entropie ---
  write(iunit, '(A)') "set xlabel 'x (m)'"
  write(iunit, '(A)') "set ylabel 's'"
  write(iunit, '(A)') "set title 'Entropie'"
  write(iunit, '(A)') "set grid"
  write(iunit, '(A)') "plot filename_final using 1:8 with lines lw 2 lc rgb 'brown' title sprintf('t=%.3fs', t_final)"
  write(iunit, *) ""

  write(iunit, '(A)') "unset multiplot"
  write(iunit, *) ""
  write(iunit, '(A)') "print 'Plots generes : Geometrie_tuyere.png et Solution_Numerique.png'"

  close(iunit)

end subroutine create_gnuplot_script



  ! ----------------------------------------------------------------------
end module Module_sortie


