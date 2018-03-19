# if backup in spectral space
spectral="1"

# 0 to recover from iteration, 1 from iteration_old
old="0"


if [ "$old" == "0" ]; then
  iteration="`cat iteration.dat`"
  cp ./time_check.dat ../time_check.dat
  cp ./stats.dat ../stats.dat
  cp ./budget.dat ../budget.dat
  cp ./power_xspectra.dat ../power_xspectra.dat
  cp ./power_yspectra.dat ../power_yspectra.dat
elif [ "$old" == "1" ]; then
  iteration="`cat iteration_old.dat`"
  cp ./time_check_old.dat ../time_check.dat
  cp ./stats_old.dat ../stats.dat
  cp ./budget_old.dat ../budget.dat
  cp ./power_xspectra_old.dat ../power_xspectra.dat
  cp ./power_yspectra_old.dat ../power_yspectra.dat
fi

echo " "
# read
read it <<<${iteration//[^0-9]/ }
printf -v it "%08d" $it




if [ "$spectral" == "1" ]; then
 echo "Restarting file from iteration $it complex space"
 printf -v it "c_$it" $it
 sp="c"
elif [ "$spectral" == "0" ]; then
 echo "Restarting file from iteration $it real space"
 printf -v it "_$it" $it
sp=""
fi

echo " "

cp "./u"$sp".dat" "../u"$it".dat"
cp "./v"$sp".dat" "../v"$it".dat"
cp "./w"$sp".dat" "../w"$it".dat"
cp "./phi"$sp".dat" "../phi"$it".dat"
cp "./psi"$sp".dat" "../psi"$it".dat"
cp "./T"$sp".dat" "../T"$it".dat"
