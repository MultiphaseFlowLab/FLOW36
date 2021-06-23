make clean
rm -r *.mod
make &> /dev/null

rm -r output
mkdir output

make

NTASK="1"

mpirun -n $NTASK ./marangoni
