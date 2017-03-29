gnuplot -e "
set term dumb enhanced; 
set key bottom;
set notics;
set xlabel 'displacement X';
set ylabel 'reaction X';
plot 'displacement-reaction.dat' u 1:2 title 'kratos test' w l;
"

gnuplot -e "
set term pngcairo;
set output 'displacement-reaction.png';
set key bottom;
set xlabel 'displacement X';
set ylabel 'reaction X';
plot 'displacement-reaction.dat' u 1:2 title 'kratos test' w lp,
     './MatlabResults/Plasticity2D.dat' u 1:2 title 'Matlab Reference' w lp;
"
