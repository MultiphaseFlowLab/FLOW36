

####Ulisses setup, only one simulation at time for the moment
##### SET UP of MULTIBLOCK RUN
#Ulisses

TAG="DROPSCO"

#### WORKING DIRECTORIES & FILENAMES
      NAME="xforplot"
      wdr=$(pwd) ######<-----Name of folder for post-process
      cd $wdr
      chmod 777 compile.sh
      ./compile.sh






echo "Creating the serial-launcher "
echo $'\360\237\232\200 \360\237\232\200 \360\237\232\200'
##### EDIT THE LAUNCHER FILE GO_ULISSES
touch launcher   #create the file go_carole
cat > launcher << "EOF"
#!/bin/bash
EOF
echo "#PBS -l nodes=1:ppn=1" >>launcher
echo "#PBS -N"$TAG >> launcher  #option
cat >> launcher << "EOF"
#PBS -l walltime=11:59:00
#PBS -q regular
#PBS -V
#PBS -A cmarchio
EOF
echo "cd "$wdr >> launcher
echo "module load fftw/3.3.4/gnu/4.9.2 " >> launcher
echo "./"$NAME >>launcher
cat >> launcher << "EOF"
EOF

qsub launcher



