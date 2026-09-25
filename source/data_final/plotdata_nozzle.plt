t_final =  0.007
filename_final = sprintf('data_final/results_at_t=%.3fs.dat', t_final)
 
set terminal pngcairo enhanced font 'arial,10' size 1200, 400
set output 'data_final/Geometrie_tuyere.png'
 
set xlabel 'x (m)'
set ylabel 'A(x)'
set title 'Geometrie de la tuyere'
set grid
plot filename_final using 1:2 with lines lw 2 lc rgb 'black' title 'A(x)'
 
set terminal pngcairo enhanced font 'arial,10' size 1400, 900
set output 'data_final/Solution_Numerique.png'
 
set multiplot layout 3,2 title 'Solution Numerique Tuyere'
 
set xlabel 'x (m)'
set ylabel 'rho'
set title 'Masse volumique [kg/m3]'
set grid
plot filename_final using 1:3 with lines lw 2 lc rgb 'blue' title sprintf('t=%.3fs', t_final)
 
set xlabel 'x (m)'
set ylabel 'Mach '
set title 'Nombre de Mach'
set grid
plot filename_final using 1:4 with lines lw 2 lc rgb 'red' title sprintf('t=%.3fs', t_final)
 
set xlabel 'x (m)'
set ylabel 'P'
set title 'Pression [Pa]'
set grid
plot filename_final using 1:5 with lines lw 2 lc rgb 'green' title sprintf('t=%.3fs', t_final)
 
set xlabel 'x (m)'
set ylabel 'T'
set title 'Temperature [K]'
set grid
plot filename_final using 1:6 with lines lw 2 lc rgb 'orange' title sprintf('t=%.3fs', t_final)
 
set xlabel 'x (m)'
set ylabel 'rho_u_A'
set title 'Debit [kg/s]'
set grid
plot filename_final using 1:7 with lines lw 2 lc rgb 'purple' title sprintf('t=%.3fs', t_final)
 
set xlabel 'x (m)'
set ylabel 's'
set title 'Entropie'
set grid
plot filename_final using 1:8 with lines lw 2 lc rgb 'brown' title sprintf('t=%.3fs', t_final)
 
unset multiplot
 
print 'Plots generes : Geometrie_tuyere.png et Solution_Numerique.png'
