rm -r output
mkdir ./output

make clean
make &> /dev/null
make

#for greater grids
#unlimits the stack for the calculation time
#otherwise : segmentation fault
# ligne à décommenter quand on passera sur le serveur
ulimit -s unlimited

mpirun -n 1 ./coalescence_breakup

