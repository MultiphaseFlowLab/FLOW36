make clean
# rm -r *.mod
# make &> /dev/null

rm -r output
mkdir output

make
#sleep 3600
NTASK="1"

mpirun -n $NTASK ./invariants
