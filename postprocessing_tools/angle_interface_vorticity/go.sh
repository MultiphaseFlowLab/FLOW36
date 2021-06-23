make clean

make


rm -r ./output/
mkdir output

NTASK="1"

mpirun -n "$NTASK" ./get_tke


gfortran create_pdf.f90 -o pdf_generation

./pdf_generation

