set terminal epslatex standalone
set output 'ssv_512.eps'
set output 'ssv_512.tex'

set xlabel "$N_{task}$"
set ylabel "speed-up"

#set xrange [6:14]
f(x)=2**(log(x)/log(2))/256
set logscale x 2
set logscale y

set key bottom right
set title "$512\\times512\\times513$ Vesta"

plot "./gnuplot_strong/512_16.dat" using ($1*$2):($4) with linespoints ls 4 lt rgb "red" dt 1 lw 2 title "16 tasks/node",\
   "./gnuplot_strong/512_32.dat" using ($1*$2):($4) with linespoints ls 4 lt rgb "blue" dt 1 lw 2 title "32 tasks/node",\
   "./gnuplot_strong/512_64.dat" using ($1*$2):($4) with linespoints ls 4 lt rgb "black" dt 1 lw 2 title "64 tasks/node",\
   f(x) with lines ls -1 dt 2 lw 2 title "ideal"

#plot "./gnuplot_strong/512_16.dat" using ($1*$2):($3) with linespoints ls 4 lt rgb "red" dt 1 lw 2 title "16 task/node",\
#   "./gnuplot_strong/512_32.dat" using ($1*$2):($3) with linespoints ls 4 lt rgb "blue" dt 1 lw 2 title "32 task/node",\
#   "./gnuplot_strong/512_64.dat" using ($1*$2):($3) with linespoints ls 4 lt rgb "black" dt 1 lw 2 title "64 task/node"


set output
system('latex ssv_512.tex && dvips ssv_512.dvi && ps2pdf ssv_512.ps')

system('rm ssv_512.aux')
system('rm ssv_512.dvi')
system('rm ssv_512.eps')
system('rm ssv_512.log')
system('rm ssv_512.ps')
system('rm ssv_512-inc.eps')
system('rm ssv_512.tex')