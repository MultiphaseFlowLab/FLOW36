set terminal epslatex standalone
set output 'ws_64.eps'
set output 'ws_64.tex'

set xlabel "$N_{task}$"
set ylabel "dimensionless time per time step"

#set yrange [0.5:5.5]
f(x)=1
set logscale x 2
set logscale y
set key top left

plot "./data/ws_512.dat" using ($1*$2):($4) with linespoints ls 4 lt rgb "red" dt 1 lw 2 title "512x8x8",\
   "./data/ws_1024.dat" using ($1*$2):($4) with linespoints ls 6 lt rgb "blue" dt 1 lw 2 title "1024x16x16", \
   "./data/ws_64.dat" using ($1*$2):($4) with linespoints ls 8 lt rgb "green" dt 1 lw 2 title "$64^3$", \
   f(x) with lines ls -1 dt 2 lw 2 title "ideal"


set output
system('latex ws_64.tex && dvips ws_64.dvi && ps2pdf ws_64.ps')

system('rm ws_64.aux ws_64.dvi ws_64.eps ws_64.log ws_64.ps ws_64-inc.eps ws_64.tex')
