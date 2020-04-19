set terminal epslatex standalone
set output 'ss_marconi_64.eps'
set output 'ss_marconi_64.tex'

set xlabel "$N_{task}$"
set ylabel "speed-up"

#set xrange [6:14]
f(x)=2**(log(x)/log(2))/64
set logscale x 2
set logscale y

set key bottom right

plot "./gnuplot_data/marconi/ss_512x256x257s.dat" using ($1*$2):($4) with linespoints ls 6 lt rgb "red" dt 1 lw 2 title "512x256x257 slab",\
   "./gnuplot_data/marconi/ss_512x256x257p.dat" using ($1*$2):($4) with linespoints ls 8 lt rgb "red" dt 1 lw 2 title "512x256x257 pencil",\
   "./gnuplot_data/marconi/ss_512x512x513s.dat" using ($1*$2):($4) with linespoints ls 6 lt rgb "blue" dt 1 lw 2 title "512x512x513 slab",\
   "./gnuplot_data/marconi/ss_512x512x513p.dat" using ($1*$2):($4) with linespoints ls 8 lt rgb "blue" dt 1 lw 2 title "512x512x513 pencil",\
   "./gnuplot_data/marconi/ss_1024x1024x1025s.dat" using ($1*$2):($4) with linespoints ls 6 lt rgb "black" dt 1 lw 2 title "1024x1024x1025 slab",\
   "./gnuplot_data/marconi/ss_1024x1024x1025p.dat" using ($1*$2):($4) with linespoints ls 8 lt rgb "black" dt 1 lw 2 title "1024x1024x1025 pencil",\
   f(x) with lines ls -1 dt 2 lw 2 title "ideal"

    
set output
system('latex ss_marconi_64.tex && dvips ss_marconi_64.dvi && ps2pdf ss_marconi_64.ps')

system('rm ss_marconi_64.aux')
system('rm ss_marconi_64.dvi')
system('rm ss_marconi_64.eps')
system('rm ss_marconi_64.log')
system('rm ss_marconi_64.ps')
system('rm ss_marconi_64-inc.eps')
system('rm ss_marconi_64.tex')


