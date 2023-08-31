Code is tested on MATLAB R2022a, the Tensor Toolbox used is version 3.2.1.

Run setPath.m to set up the paths involving the RedbKit package (add dependent paths), the Tensor Toolbox, necessary to run the codes inside the folder Thermal_nonlinear, can be downloaded from: https://www.tensortoolbox.org/. The code is tested with the version of the Tensor Toolbox: v3.2.1.

To reproduce the results go to Thermal_nonlinear folder and use the following scripts:

Script assemble_matrices.m  assemble necessary FEM matrices and put data into mesh_data  

Script main_passive.m run the optimization algorithm.

For the code of the linear thermal cloak paper take a look at: https://github.com/CaliShulz/PDE_ROM_thermal_cloak_PRSA
