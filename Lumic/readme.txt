EUROHPC Project CARAVEL
1 node has 2 x AMD with 64 cores (physical 128 cores per node)
1 node has 256 GB of memory 


Load manager is SLURM.

Libraries linking.
If using CRAY, not necessary to specify the path (using -crayFFTW)

Documentation on Compiler, Libraries
https://docs.lumi-supercomputer.eu/firststeps/getstarted/#access-to-lumi

Dahsboard of the project:
https://my.lumi-supercomputer.eu/

Login using RSA pub key.
Upload key at:
https://mms.myaccessid.org/fed-apps/profile/#personal


So far, no problems with Cray (compiling, linking and running is fine).

To run:
-Specify the memory (if not enough)
-Use srun (no mpirun or mpiexec), no need to specify the number of MPI task (automatically recognised by the SLURM)