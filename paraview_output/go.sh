make clean
rm -r *.mod
make &> /dev/null


rm -r output
mkdir output

make

NTASK="5"

mpirun -n $NTASK ./read_paraview
