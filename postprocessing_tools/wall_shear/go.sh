rm -r output
mkdir ./output

#make clean
make &> /dev/null
make

mpirun -n 1 ./wall_shear
