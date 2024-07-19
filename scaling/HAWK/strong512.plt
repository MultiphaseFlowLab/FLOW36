reset
set terminal epslatex standalone linewidth 2
set output 'p3.eps'
set output 'p3.tex'
set size ratio 0.9
set lmargin at screen 0.2
set rmargin at screen 0.8
set xrange [60:4100]
set x2range [60:4100]
set yrange [0.5:50]
set xtics ('64' 64,'128' 128,'256' 256,'512' 512,'1024' 1024,'2048' 2048,'4096' 4096,'8192' 8192,'16384' 16384)
set x2tics ('1' 64,'1' 128,'2' 256,'4' 512,'8' 1024,'16' 2048,'32' 4096,'64' 8192,'128' 16384)
#set ytics 0,0.2,1
set xlabel '$N_{tasks}$' offset 0,0.
set x2label '$N_{nodes}$' offset 0,0.
set ylabel 'speed-up' offset 5,0 rotate by 0
set key font ",20"
set key top left
set logscale xyx2


plot 'strongscaling.dat' u 3:4 title '$512^3$ HT 1' with lines dt 1 lc rgb '#ff0000', \
     'strongscaling.dat' u 3:5 title '$512^3$ HT 2' with lines dt 2 lc rgb '#ff0000', \
     'strongscaling.dat' u 3:($3/128) title 'Ideal' with lines dt 2 lc rgb '#000000'

# notitle with lines lt 6 lw 3 dt 1 lc rgb '#3081B7'

set output
system('latex p3.tex && dvips p3.dvi && ps2pdf p3.ps')
system('rm p3.tex p3.log p3.aux p3.dvi p3.eps p3-inc.eps p3.ps')
system('mv p3.pdf strong512.pdf')
system('open strong512.pdf')
