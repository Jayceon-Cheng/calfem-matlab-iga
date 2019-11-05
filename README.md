# MultiGrid-Isogeometric Reanalysis (MG-IGR)

MATLAB code for paper "An efficient auxiliary projection based multigrid-isogeometric reanalysis method and its application in optimization framework".

## Installation instructions

1. Make sure you have CALFEM for MATLAB installed, see https://github.com/CALFEM/calfem-matlab.

2. To install the IGA toolbox, click "Download as zip" to download the package, and then unpack it. 

3. Add the directories to the MATLAB path by clicking "Set path" in MATLAB, then "Add Folder...", select "calfem-matlab-iga/IGA" and then "Save".

4. Using the same method as in 3. add the following folders too: "calfem-matlab-iga/IGAplot", "calfem-matlab-iga/IGAutil", "calfem-matlab-iga/NURBS".

## Files description

### `Tubular_initial.m`

Case I: Tubular structural Initialization code

### `GA_nsga3.m`

Code of Genetic Algorithm optimization

### `pso.m`

Code of Particle Swarm Optimization

### `IGA_static_optimization.m`

Code for IGA static optimization

## Note

This code sets different parameters according to different computer needs, such as the path of the file addpath (genpath('..\calfem-matlab-iga-master-screwdriver'));
