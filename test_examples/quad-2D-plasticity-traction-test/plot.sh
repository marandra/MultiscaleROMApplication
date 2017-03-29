gnuplot -e "
set term dumb enhanced; 
set key bottom;
set notics;
set xlabel 'displacement X';
set ylabel 'reaction X';
plot 'displacement-reaction.dat' u 1:2 title 'case' w l;
"

gnuplot -e "
set term pngcairo;
set output 'displacement-reaction.png';
set key bottom;
set xlabel 'displacement X';
set ylabel 'reaction X';
plot 'displacement-reaction.dat' u 1:2 title 'case' w lp,
     './MatlabResults/quad01-traction_bbar_Matlab.dat' u 1:2 title 'MatlabReference' w lp;
"
