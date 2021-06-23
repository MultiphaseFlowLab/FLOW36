0               ! 1: reads files in spectral space, 0: reads files in physical space

0               ! x_start: first index to be saved in x direction, leave zero to read from 1st index
0               ! x_end: last index to be saved in x direction, leave zero to read up to last index
1               ! Delta_x: save a point in x direction every Delta_x

0               ! y_start: first index to be saved in y direction, leave zero to read from 1st index
0               ! y_end: last index to be saved in y direction, leave zero to read up to last index
1               ! Delta_y: save a point in y direction every Delta_y

0               ! z_start: first index to be saved in z direction, leave zero to read from 1st index
0               ! z_end: last index to be saved in z direction, leave zero to read up to last index
1               ! Delta_z: save a point in z direction every Delta_z

-1              ! start n step, leave -1 to use the first time step
-1              ! end n step, leave -1 to use the last time step

1               ! 0 to skip phase variable, 1 to include it
0               ! 0 to skip surfactant, 1 to include it
0               ! 0 to skip temperature, 1 to include it
1               ! 0 to skip velocity fluctuations, 1 to include them
1               ! 0 to skip vorticity, 1 to include it
1               ! 0 to skip strain rate, 1 to include it
