make clean
#rm -r *.mod
#make &> /dev/null


rm -r output/*

make

NTASK="1"

mpirun -n $NTASK ./curvature


# gnuplot --persist gnuplot.plt
