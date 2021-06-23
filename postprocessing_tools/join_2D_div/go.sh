rm join_data
rm -rf output

mkdir output

gfortran -Wall main.f90 -o join_data

./join_data
