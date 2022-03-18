make clean
rm -r *.mod
make &> /dev/null


rm -r output/*

make

NTASK="1"

mpirun -n $NTASK ./change_grid

