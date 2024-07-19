set terminal epslatex standalone
set output 'tv_1024.eps'
set output 'tv_1024.tex'

set xlabel "$N_{task}$"
set ylabel "Time/time step [s]"

set logscale x 2
set logscale y

set key top right
set title "$1024\\times1024\\times1025$ Vesta"

plot "./gnuplot_strong/1024_16.dat" using ($1*$2):($3) with linespoints ls 4 lt rgb "red" dt 1 lw 2 title "16 task/node",\
   "./gnuplot_strong/1024_32.dat" using ($1*$2):($3) with linespoints ls 4 lt rgb "blue" dt 1 lw 2 title "32 task/node",\
   "./gnuplot_strong/1024_64.dat" using ($1*$2):($3) with linespoints ls 4 lt rgb "black" dt 1 lw 2 title "64 task/node"

    
set output
system('latex tv_1024.tex && dvips tv_1024.dvi && ps2pdf tv_1024.ps')

system('rm tv_1024.aux')
system('rm tv_1024.dvi')
system('rm tv_1024.eps')
system('rm tv_1024.log')
system('rm tv_1024.ps')
system('rm tv_1024-inc.eps')
system('rm tv_1024.tex')


