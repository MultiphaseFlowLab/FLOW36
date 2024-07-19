set terminal epslatex standalone
set output 'tv_512.eps'
set output 'tv_512.tex'

set xlabel "$N_{task}$"
set ylabel "Time/time step [s]"

set logscale x 2
set logscale y

set key top right
set title "$512\\times512\\times513$ Vesta"

plot "./gnuplot_strong/512_16.dat" using ($1*$2):($3) with linespoints ls 4 lt rgb "red" dt 1 lw 2 title "16 task/node",\
   "./gnuplot_strong/512_32.dat" using ($1*$2):($3) with linespoints ls 4 lt rgb "blue" dt 1 lw 2 title "32 task/node",\
   "./gnuplot_strong/512_64.dat" using ($1*$2):($3) with linespoints ls 4 lt rgb "black" dt 1 lw 2 title "64 task/node"

    
set output
system('latex tv_512.tex && dvips tv_512.dvi && ps2pdf tv_512.ps')

system('rm tv_512.aux')
system('rm tv_512.dvi')
system('rm tv_512.eps')
system('rm tv_512.log')
system('rm tv_512.ps')
system('rm tv_512-inc.eps')
system('rm tv_512.tex')


