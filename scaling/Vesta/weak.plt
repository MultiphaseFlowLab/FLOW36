set terminal epslatex standalone
set output 'ws_64.eps'
set output 'ws_64.tex'

set xlabel "$N_{task}$"
set ylabel "dimensionless time per time step"

set yrange [0.5:15]
f(x)=1
set logscale x 2
set logscale y

plot "./gnuplot_data/ws_512.dat" using ($1*$2):($4) with linespoints ls 4 lt rgb "red" dt 1 lw 2 title "512x8x8",\
   "./gnuplot_data/ws_1024.dat" using ($1*$2):($4) with linespoints ls 6 lt rgb "blue" dt 1 lw 2 title "1024x16x16", \
   "./gnuplot_data/ws_64.dat" using ($1*$2):($4) with linespoints ls 8 lt rgb "green" dt 1 lw 2 title "$64^3$", \
   f(x) with lines ls -1 dt 2 lw 2 title "ideal"


set output
system('latex ws_64.tex && dvips ws_64.dvi && ps2pdf ws_64.ps')

system('rm ws_64.aux')
system('rm ws_64.dvi')
system('rm ws_64.eps')
system('rm ws_64.log')
system('rm ws_64.ps')
system('rm ws_64-inc.eps')
system('rm ws_64.tex')

