reset
set terminal epslatex standalone linewidth 2
set output 'p3.eps'
set output 'p3.tex'
set size ratio 0.9
set lmargin at screen 0.2
set rmargin at screen 0.8
set xrange [60:66000]
set yrange [0.5:180]
set xtics ('64' 64,'256' 256,'1024' 1024,'4096' 4096,'16384' 16384,'65536' 65536)
#set ytics 0,0.2,1
set xlabel '$N_{tasks}$' offset 0,0.
set ylabel 'speed-up' offset 5,0 rotate by 0
set key font ",20"
set key top left
set logscale xy

plot 'strongscalingflow.dat' u 3:4 title '$512^3$ HT 1' with lines dt 1 lc rgb '#ff0000', \
     'strongscalingflow.dat' u 3:5 title '$512^3$ HT 2' with lines dt 2 lc rgb '#ff0000', \
     'strongscalingflow.dat' u 3:6 title '$1024^3$ HT 1' with lines dt 1 lc rgb '#0008ff', \
     'strongscalingflow.dat' u 3:($3/128) title 'Ideal' with lines dt 2 lc rgb '#000000'

# notitle with lines lt 6 lw 3 dt 1 lc rgb '#3081B7'

set output
system('latex p3.tex && dvips p3.dvi && ps2pdf p3.ps')
system('rm p3.tex p3.log p3.aux p3.dvi p3.eps p3-inc.eps p3.ps')
system('mv p3.pdf flowstrong.pdf')
system('open flowstrong.pdf')
