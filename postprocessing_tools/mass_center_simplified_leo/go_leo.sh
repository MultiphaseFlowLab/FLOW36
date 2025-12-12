module purge
module load gcc/11.3.0
module load openmpi/4.1.4--gcc--11.3.0-cuda-11.8
module load fftw/3.3.10--openmpi--4.1.4--gcc--11.3.0

make clean
rm -r *.mod
make &> /dev/null


rm -r output
mkdir output

make

./mass_center
