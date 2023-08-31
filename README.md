Code is tested on MATLAB R2022a, the Tensor Toolbox used is version 3.2.1.

Run setPath.m to set up the paths involving the RedbKit package (https://redbkit.github.io/redbKIT/) (add dependent paths), the Tensor Toolbox, necessary to run the codes inside the folder Thermal_nonlinear, can be downloaded from: https://www.tensortoolbox.org/. The code is tested with the version of the Tensor Toolbox: v3.2.1 (https://gitlab.com/tensors/tensor_toolbox/-/releases/v3.2.1).

To reproduce the results go to Thermal_nonlinear folder and use the following scripts:

Script assemble_matrices.m  assemble necessary FEM matrices and put data into mesh_data  

Script main_passive_stationary.m run the optimization algorithm in stationary conditions.

The folder ufv_time_dependent.m contains analogous functions to run the optimization in time-dependent conditions.

That two (2) .mat files: output_ufv_with_Hessian.mat and output_ufv_time_dep_hessian_2.mat have already been created to store the results for the paper.

For the code of the linear thermal cloak paper take a look at: https://github.com/CaliShulz/PDE_ROM_thermal_cloak_PRSA
