gnuplot -e "
set term dumb enhanced; 
set key bottom;
set notics;
set xlabel 'displacement X';
set ylabel 'reaction X';
plot 'displacement-reaction.dat' u 3:4 title 'case' w l;
"

gnuplot -e "
set term pngcairo;
set output 'displacement-reaction.png';
set key bottom;
set xlabel 'displacement Y';
set ylabel 'reaction Y';
plot 'displacement-reaction.dat' u 3:4 title 'case' w lp,
     'displacement-reaction.dat' u 3:4 title 'goldstandard' w lp;
"
