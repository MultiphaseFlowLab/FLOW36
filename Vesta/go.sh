qsub -t 2:00 -n 64 --proccount 1024 --mode c16 -o output.out ./sc_compiled/flow36

# run script to submit

# -t : time (minutes)
# -n : nodes
# --proccount : number of ranks (total), equal to mode*nodes
# --mode : number of rank per node/core
# -A : account
# -o : file where to write output
