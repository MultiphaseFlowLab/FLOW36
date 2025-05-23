# FLOW36
Github Repository of FLOW36. \
FLOW36 is a pseudo-spectral code for direct numerical simualtion of multiphase turbulence based on a phase-field model approach.


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


## How to:
See the handbook in the flow36_handbook folder


## Packages/libraries required:
 - Fortran compiler (tested wtth gfortran, ifort, nvfortran, ftn and xlf)
 - MPI Library (tested with MPICH, Spectrum, IntelMPI, OpenMPI)
 - FFTW library (for CPU runs)
 - cuFFT library (for GPU runs, please install the Nvidia HPC-SDK)


