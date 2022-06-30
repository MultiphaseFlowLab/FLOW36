EUROHPC Project BLUMEN
1 node has 2 x AMD with 64 cores (physical 128 cores per node)
1 node has 256 GB of memory 


Load manager is SLURM.

Important:
When using more than 1 node, include:
export UCX_NET_DEVICES=mlx5_0:1

to your Slurm job file when using more than 1 node for MPI run. 
Otherwise the node-2-node communication uses TCP, not InfiniBand.



To check the comp. Time consumed:

Run the tool get_breakdonw.py in the following manner (replace the dates if you like to get accounting breakdown for a different period, different than in the example bellow):
/opt/software/tools/bin/get_breakdown.py aroccon 2022-06-01 2022-07-01
The output is CVS formatted. The positions in the output line are as follows:
-username
-threadseconds
-threadhours
-corehours
-ave_num_threads
-ave_num_cores
-ave_num_nodes
-corehours_per_node
-coredays_per_node


Connection to the machine is via VPN + openSSH key.
Documentation:
https://docs.discoverer.bg/resource_overview.html
