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
 # printf -v it "$it" $it
 sp="c"
elif [ "$spectral" == "0" ]; then
 echo "Restarting file from iteration $it real space"
 # printf -v it "_$it" $it
sp=""
fi

echo " "

if [ "$old" == "0" ]; then
  cp "./u"$sp".dat" "../u"$sp"_"$it".dat"
  cp "./v"$sp".dat" "../v"$sp"_"$it".dat"
  cp "./w"$sp".dat" "../w"$sp"_"$it".dat"
  cp "./phi"$sp".dat" "../phi"$sp"_"$it".dat"
  cp "./psi"$sp"_fg.dat" "../psi"$sp"_fg_"$it".dat"
  cp "./T"$sp".dat" "../T"$sp"_"$it".dat"
  cp "./pos"$sp".dat" "../pos"$sp"_"$it".dat"
  cp "./vel"$sp".dat" "../vel"$sp"_"$it".dat"
elif [ "$old" == "1" ]; then
  cp "./u"$sp"_old.dat" "../u"$sp"_"$it".dat"
  cp "./v"$sp"_old.dat" "../v"$sp"_"$it".dat"
  cp "./w"$sp"_old.dat" "../w"$sp"_"$it".dat"
  cp "./phi"$sp"_old.dat" "../phi"$sp"_"$it".dat"
  cp "./psi"$sp"_fg_old.dat" "../psi"$sp"_fg_"$it".dat"
  cp "./T"$sp"_old.dat" "../T"$sp"_"$it".dat"
  cp "./pos"$sp"_old.dat" "../pos"$sp"_"$it".dat"
  cp "./vel"$sp"_old.dat" "../vel"$sp"_"$it".dat"
fi
