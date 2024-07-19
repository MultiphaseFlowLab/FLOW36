set terminal epslatex standalone
set output 'ssv_64.eps'
set output 'ssv_64.tex'

set xlabel "$N_{task}$"
set ylabel "speed-up"

#set xrange [6:14]
f(x)=2**(log(x)/log(2))/64
set logscale x 2
set logscale y

set key bottom right

plot "./gnuplot_data/ss_1024.dat" using ($1*$2):($4) with linespoints ls 4 lt rgb "red" dt 1 lw 2 title "1024x1024x1025 Vesta",\
   f(x) with lines ls -1 dt 2 lw 2 title "ideal"

    
set output
system('latex ssv_64.tex && dvips ssv_64.dvi && ps2pdf ssv_64.ps')

system('rm ssv_64.aux')
system('rm ssv_64.dvi')
system('rm ssv_64.eps')
system('rm ssv_64.log')
system('rm ssv_64.ps')
system('rm ssv_64-inc.eps')
system('rm ssv_64.tex')


