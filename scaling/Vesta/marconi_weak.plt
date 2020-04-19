set terminal epslatex standalone
set output 'marconi_ws_64.eps'
set output 'marconi_ws_64.tex'

set xlabel "$N_{task}$"
set ylabel "dimensionless time per time step"

set yrange [0.5:60]
f(x)=1
set logscale x 2
set logscale y

set key top left

plot "./gnuplot_data/ws_marconi512.dat" using ($1*$2):($4) with linespoints ls 4 lt rgb "red" dt 1 lw 2 title "512x8x8",\
   "./gnuplot_data/ws_marconi1024.dat" using ($1*$2):($4) with linespoints ls 6 lt rgb "blue" dt 1 lw 2 title "1024x16x16", \
   "./gnuplot_data/ws_marconi64.dat" using ($1*$2):($4) with linespoints ls 8 lt rgb "green" dt 1 lw 2 title "$64^3$", \
   f(x) with lines ls -1 dt 2 lw 2 title "ideal"


set output
system('latex marconi_ws_64.tex && dvips marconi_ws_64.dvi && ps2pdf marconi_ws_64.ps')

system('rm marconi_ws_64.aux')
system('rm marconi_ws_64.dvi')
system('rm marconi_ws_64.eps')
system('rm marconi_ws_64.log')
system('rm marconi_ws_64.ps')
system('rm marconi_ws_64-inc.eps')
system('rm marconi_ws_64.tex')

