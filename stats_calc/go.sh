make clean
rm -r *.mod
make &> /dev/null


rm -r output
mkdir output

make

# rank 0 calculates u stats, rank 1 v stats, rank 2 w stats
NTASK="3"

mpirun -n $NTASK ./stats_calc
