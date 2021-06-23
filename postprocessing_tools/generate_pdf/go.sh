rm get_pdf
rm -rf output

mkdir output

gfortran -Wall main.f90 -o get_pdf

./get_pdf
