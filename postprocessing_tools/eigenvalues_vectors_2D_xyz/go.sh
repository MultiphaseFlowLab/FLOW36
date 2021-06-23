make clean
#rm -r *.mod
#make &> /dev/null


rm -r output/*

make

NTASK="1"

# for large grids
ulimit -s unlimited

mpirun -n $NTASK ./mass_center

# gnuplot --persist gnuplot.plt
