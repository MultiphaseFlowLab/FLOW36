set terminal epslatex standalone
set output 'ss_64.eps'
set output 'ss_64.tex'

set xlabel "$N_{task}$"
set ylabel "speed-up"

#set xrange [6:14]
f(x)=2**(log(x)/log(2))/64
set logscale x 2
set logscale y

set key bottom right

plot "./gnuplot_data/ss_1024.dat" using ($1*$2):($4) with linespoints ls 4 lt rgb "red" dt 1 lw 2 title "1024x1024x1025 Vesta",\
   "./gnuplot_data/ss_marconi.dat" using ($1*$2):($4) with linespoints ls 8 lt rgb "blue" dt 1 lw 2 title "1024x1024x1025 Marconi", \
   f(x) with lines ls -1 dt 2 lw 2 title "ideal"

    
set output
system('latex ss_64.tex && dvips ss_64.dvi && ps2pdf ss_64.ps')

system('rm ss_64.aux')
system('rm ss_64.dvi')
system('rm ss_64.eps')
system('rm ss_64.log')
system('rm ss_64.ps')
system('rm ss_64-inc.eps')
system('rm ss_64.tex')


