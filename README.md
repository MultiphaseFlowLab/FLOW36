# FLOW36
Github Repository of FLOW36. \
FLOW36 is a pseudo-spectral code for direct numerical simulation of multiphase turbulence based on a phase-field model approach.


If you use this code, please cite the following work: 

```bibtex
@article{roccon2025,
title = {FLOW36: A spectral solver for phase-field based multiphase turbulence simulations on heterogeneous computing architectures},
journal = {Computer Physics Communications},
volume = {313},
pages = {109640},
year = {2025},
issn = {0010-4655},
author = {Roccon, A. and Soligo, G. and Soldati, A.},
}
```

Webminar on FLOW36 available on [Cassyni](https://cassyni.com/events/QYj3h1ohsebWbK8AFdiPJd)


## Main Developers:
G. Soligo (https://github.com/giovannisoligo) \
A. Roccon (https://github.com/aroccon)


![](flow36_handbook/cop.jpeg)



## Modules available in FLOW36:
 - Single-phase flow (NS): Turbulent channel flow, close and open configurations 
 - Phase-field method (NS + CH): Clean intergaces (surfactant-free): Standard CH + Profile and Flux-corrected 
 - Phase-field method (NS + CH1+ CH2): Surfactant-laden interfaces (possible use of the dual grid) 
 - Phase-field method + scalar (NS + CH + EE): Heat/mass transfer in drop- and bubble-laden flows 
 - Phase-field method + particles (NS + CH + LPT): Interface-particle interactions  
 - Single-phase + passive scalar (NS + EE): Heat transfer in single-phase turbulence
 - Single-phase + particles (NS + LPT): Particle-laden turbulent flows 
 - Single-phase + temperature + particles (NS + EE + LPT): Particle-laden turbulent flows with temperature 


## Published works:
Click [here](http://calliope.dem.uniud.it) for a list of the published works


## Validation data and initial fields:
Click [here](https://doi.org/10.6084/m9.figshare.26232683) for the validation data shown in the manuscript

## How to run the code:
- **Clone the repository.**
- **Edit the `compile.sh` file:**  
  Select the configurations you want to run and the modules to load. Many EUROHPC supercomputer configurations are included, along with options for different compilers such as GNU, NVIDIA, AMD, Intel, IBM, etc. Simulation parameters can be defined after the machine configuration section
- **Run the `compile.sh` script:**  
  This will create the `set_run` folder. Inside, you will find the `results` folder and the `sc_compiled` folder, which contains the compiled source code and the main executable (`flow36`)
- **Run the simulation:**  
  Use `mpirun` or edit/create a SLURM file for your target machine
- **Access the results:**  
  The output fields will be located in `set_run/results` and can be visualized using ParaView. ParaView files can be generated with the `paraview_output_fg` post-processing code

## Packages/libraries required:
 - For CPU runs: 
   - Fortran compiler (tested with gfortran, ifort, nvfortran, ftn and xlf)
   - MPI Library (tested with MPICH, Spectrum, IntelMPI, OpenMPI)
   - FFTW library or MKL library (if using intel)
- For GPU runs:
   - Nvidia HPC-SDK (> 22.X), this by default contains the compiler, the cuFFT and the openMPI libraries.

## Contributing
We welcome all contributions that can enhance FLOW36, including bug fixes, performance improvements, and new features. 
If you would like to contribute, please contact aroccon or open an Issue in the repository.
