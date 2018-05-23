gnuplot -e "
set term dumb enhanced; 
set key bottom;
set notics;
set xlabel 'Homogenized Strain XX';
set ylabel 'Homogenized Stress XX';
plot 'homogenized_stress.dat' u 0:1 title 'case' w l;
"

gnuplot -e "
set term pngcairo;
set output 'displacement-reaction_xx.png';
set key bottom;
set xlabel 'Homogenized Strain XX';
set ylabel 'Homogenized Stress XX';
plot 'homogenized_stress.dat' u 0:1 title 'case' w lp,
     './MatlabResults/Micro2D-Plasticity.cur001' u 0:2 title 'MatlabReference' w lp;
"

gnuplot -e "
set term pngcairo;
set output 'displacement-reaction_yy.png';
set key bottom;
set xlabel 'Homogenized Strain YY';
set ylabel 'Homogenized Stress YY';
plot 'homogenized_stress.dat' u 0:2 title 'case' w lp,
     './MatlabResults/Micro2D-Plasticity.cur002' u 0:2 title 'MatlabReference' w lp;
"

gnuplot -e "
set term pngcairo;
set output 'displacement-reaction_xy.png';
set key bottom;
set xlabel 'Homogenized Strain XY';
set ylabel 'Homogenized Stress XY';
plot 'homogenized_stress.dat' u 0:4 title 'case' w lp,
     './MatlabResults/Micro2D-Plasticity.cur003' u 0:2 title 'MatlabReference' w lp;
"
