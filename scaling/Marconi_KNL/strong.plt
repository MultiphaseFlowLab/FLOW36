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

plot "./data/ss_512x512x513.dat" using ($1*$2):($4) with linespoints ls 4 lt rgb "red" dt 1 lw 2 title "512x512x513 Marconi KNL",\
   "./data/ss_1024x1024x1025.dat" using ($1*$2):($4) with linespoints ls 8 lt rgb "blue" dt 1 lw 2 title "1024x1024x1025 Marconi KNL", \
   f(x) with lines ls -1 dt 2 lw 2 title "ideal"


set output
system('latex ss_64.tex && dvips ss_64.dvi && ps2pdf ss_64.ps')

system('rm ss_64.aux ss_64.dvi ss_64.eps ss_64.log ss_64.ps ss_64-inc.eps ss_64.tex')
