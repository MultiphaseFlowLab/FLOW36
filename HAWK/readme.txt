Prace Project RUBIN
Storage is organized in workspace
Take care of the following things: IT EXPIRES, stanadrad time is 10 days. 
After this period of time, the ws is deleted.
To access the workspace, type ws_list and then cd to the workspace path.

Minimun size of the job is 8192 MPI threads (e.g ny x nz = 64 * 128), which corresponds to 64 nodes.
System Architecture:
AMD Rome Epyc: 2X AMD Epyc Rome 7742, 64 cores/socket, 128 cores/node, support hyperthreading 2x
Check job details with qstat -a or batchstat

Details for job submission at:
https://kb.hlrs.de/platforms/index.php/Batch_System_PBSPro_(Hawk)



