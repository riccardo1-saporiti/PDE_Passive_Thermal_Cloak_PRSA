clearvars
clc
close all
%Solve the stationary optimization problem and plot Figure 2, Figure 3,
%Figure 4, Figure 5, Figure 6 according to the mesh selected by the user at
%line 32-33 

if ispc
    sslash = '\';
elseif isunix
    sslash = '/';
end

addpath("useful_functions" , "results_folder");


% Setup fonts for plots
font_label = 18;
font_title = 19;
font_legend = 10;

fonts_data.font_title = font_title;
fonts_data.font_label = font_label;

% Define input data
mu_0 = 1; % background diffusibity
param.mu_0 = mu_0;
T_dir = 0;
I_s   = 100; %source intensity


data_name = "data_setup_coarse_grid";               % Change here for different mesh 
% data_name = "data_setup_coarse_boar_cloak_shape";               % Change here for different mesh 

data_name = strcat("mesh_data",sslash,data_name);    
load(data_name);
B_dd_u = FOM.B_dd_u;
B_dd_f = FOM.B_dd_f;
B_dl = FOM.B_dl;

A = mu_0*FOM.A_d_ocp+FOM.A_d_robin_ocp;  % + FOM.M_ocp;
F = FOM.F_ocp;
E = FOM.E;
z  = (mu_0*FOM.A_d+FOM.A_d_robin) \ (I_s*FOM.F);                           % Solve for reference dynamics (no obstacle)
%solve for z at beginning using the matrix not projected
y  = (mu_0*FOM.A_d_ocp+FOM.A_d_robin_ocp) \ (I_s*FOM.F_ocp);               % Solve for uncontrolled dynamics (diffusivity inside the cloak is mu_0)
Nz = size(z,1);

fig = gobjects(0);
set(0,'DefaultFigureVisible','on');

beta        = 1e-9;     
param.beta  = beta;
beta_g = 7e-6;
xi = 1e-9;
xi_g = 7e-6;
gamma = 1e-9;
gamma_g = 5e-5;
param.beta_g = beta_g;
param.xi = xi;
param.xi_g = xi_g;
param.gamma = gamma;
param.gamma_g = gamma_g;
param.A     = A;
param.B_dd_u  = B_dd_u;
param.B_dd_f = B_dd_f;
param.B_dl = B_dl;
param.M_o   = FOM.M_obs;         
M_o         = param.M_o;
param.M_u   = FOM.M_u;
M_u         = param.M_u;
param.A_u   = FOM.A_u;
A_u         = param.A_u;
param.E     = E;
param.z     = z;
param.T_dir = T_dir;
param.F     = F;
param.I_s   = I_s;

Nu    = size(B_dd_u,3);
param.Nu = Nu;
decrease    = 0.1;
tol_mu = 0.1 * mu_0 ;
mu_min  = 0 + tol_mu; 
param.mu_min = mu_min;
param.tol_mu = tol_mu;               
B_dd_u_dir = FOM.B_dd_u_dir;
param.B_dd_u_dir = B_dd_u_dir;
B_dd_f_dir = FOM.B_dd_f_dir;
param.B_dd_f_dir = B_dd_f_dir;
B_dl_dir = FOM.B_dl_dir;
param.B_dl_dir = B_dl_dir;
A_d_dir = FOM.A_d_dir;
param.A_d_dir = A_d_dir;boundary_dof_ocp = FOM.boundary_dof_ocp;
A_d_robin_dir = FOM.A_d_robin_dir;
param.A_d_robin_dir = A_d_robin_dir;
T_dir_vect = T_dir * sparse(ones(length(boundary_dof_ocp),1)); 
param.T_dir_vect = T_dir_vect;
param.T_dir = T_dir;
A_in_dir = mu_0 * A_d_dir + A_d_robin_dir;
param.A_in_d = A_in_dir;

%% Optimization 
z_0 = zeros( 3 * Nu,1 );
lb = [];
tol_constraint = 1e-3;
param.tol_constraint = tol_constraint;


Acon = [];
b = [];
Aeq = [];
beq = []; 
ub= [];
param.constraint_type = 'u_f_v';
param.control_type = 'u_f_v';
param.b_c = 'dirichlet';

tic
[xsol,fval,history,grad_norm,lambda] = runfmincon_ufv_Hessian(@( x ) cost_with_grad_optimization( x , param) , @( x ) constraint( x , param ) , param , z_0 , lb );
toc
fval_history = [ history.fval ];

u_opt = xsol( 1 : Nu ,1 );
f_opt = xsol( Nu + 1 : 2 * Nu , 1 );
v_opt = xsol( 2 * Nu + 1 : end , 1 );
B_u = spmatrix(ttv(B_dd_u,u_opt,3));                      % tensor toobox functions to multiply tensor B_dd and a vector u_opt, returns a sparse matrix
B_f = spmatrix( ttv( B_dd_f , f_opt , 3 ) );
B_v = spmatrix( ttv( B_dl , v_opt , 3 ) );
if param.T_dir~= 0
        A_in_dir = param.A_in_d;
        B_dd_u_dir = param.B_dd_u_dir;
        B_dd_dir_u = sparse( double( ttv( B_dd_u_dir , u , 3 ) ) );
        B_dd_f_dir = param.B_dd_f_dir;
        B_dd_dir_f = sparse( double( ttv( B_dd_f_dir , f , 3 ) ) );
        B_dl_dir = param.B_dl_dir;
        B_dl_v = spmatrix( ttv(B_dl_dir,v,3));
        T_dir_vect = param.T_dir_vect;
        F = F - ( A_in_dir + B_dd_dir_u + B_dd_dir_f+B_dl_v  ) * T_dir_vect / I_s;
end
y_opt      = (A+B_u+B_v+B_f)  \ (I_s*F);                          % optimal state 
p_opt      = (A+B_u+B_v+B_f)' \ (M_o*(y_opt-E*z));                % optimal adjoint
loglog( 1 : length( fval_history ) , fval_history , 'Linewidth' , 1 ) 
title( "fmincon, \beta_{g} = " + num2str( beta_g ) + ", \gamma_{g} = " + num2str( gamma_g ) + ", \xi_{g} = " + num2str( xi_g ) )
ylabel( '$J_{history}$' , 'Interpreter' , 'Latex' , 'FontSize',20 )
xlabel( '$Iterations$' , 'Interpreter' , 'Latex' , 'FontSize',20 )
figure()
norm_grad = grad_norm;
isnotzerograd = find( norm_grad < 1.5 );
loglog( 1 :  length( isnotzerograd ), norm_grad( isnotzerograd ), 'Linewidth' , 1)
ylabel( '$\nabla{J}_{history}$' , 'Interpreter' , 'Latex' , 'FontSize',20 )
xlabel( '$Iterations$' , 'Interpreter' , 'Latex' , 'FontSize',20 )



  
num_mte_star = sum( ( FOM.E_obs * y_opt - FOM.E_obs * param.E * param.z ) .^ 2 );
num_mte = sum( ( FOM.E_obs * y - FOM.E_obs * param.E * param.z  ) .^ 2 );
den_mte = size( FOM.E_obs , 1 );
MTE_star = sqrt( num_mte_star / den_mte );
MTE = sqrt( num_mte / den_mte );
eta_val =abs( MTE - MTE_star ) / ( MTE );



%% Computation of the eigenvalues of the diffusivity matrix
delta = ( u_opt - f_opt ) .^ 2 + 4 * v_opt .^ 2;
lambda_1 = ( 2 * mu_0 + u_opt + f_opt + sqrt( delta ) ) * 0.5;
lambda_2 = ( 2 * mu_0 + u_opt + f_opt - sqrt( delta ) ) * 0.5;
if lambda_1 > 0 
    disp( 'The first Eigenvalue is positive' )
end
if lambda_2 > 0 
    disp( 'The second Eigenvalue is positive' )
end


%% Plot of the eigenvalues 
ctrl_data.name = "passive_control_u_fom";
ctrl_data.y    =lambda_2;
ctrl_data.mesh = FOM.MESH;

ctrl_plot_data.title = "$Eigenvalue\:\lambda_{2}$";
ctrl_plot_data.limits = [min(ctrl_data.y) max(ctrl_data.y)];

% Control lives on a reduced mesh, get reduced mesh
[reduced_control_elements,~] = get_reduced_mesh(FOM.MESH,FOM.control_basis_index);
ctrl_data.reduced.vertices     = FOM.MESH.vertices(:,FOM.control_basis_index);
ctrl_data.reduced.elements     = reduced_control_elements; 
ctrl_data.reduced.indexes      = 1:length(FOM.control_basis_index);

fig(length(fig)+1)  = figure;
[fig] = plot_field(fig,ctrl_data,ctrl_plot_data,fonts_data);

%% Plot of the eigenvectors with the diffusivity matrix built by the 3 controls
alpha_e_1 = - v_opt ./ ( mu_0 + u_opt - lambda_1 );
beta_e_1 = ones( length( alpha_e_1 ) , 1 );
norm_1 = sqrt( beta_e_1 .^ 2 + alpha_e_1 .^ 2 );
alpha_e_1_norm = alpha_e_1 ./ norm_1;
beta_e_1_norm = beta_e_1 ./ norm_1;
v_1 = [ alpha_e_1 ./ norm_1 ; beta_e_1 ./ norm_1 ] ;
ang_1 = rad2deg( atan( beta_e_1_norm ./ alpha_e_1_norm ) );
vert_cont = FOM.vertices_control;
ctrl_data.name = "passive_control_u_fom";
ctrl_data.y    =ang_1;
ctrl_data.mesh = FOM.MESH;

ctrl_plot_data.title = "$\angle(\textbf{w}_{1},x_{1}) [^{\circ}]$";
ctrl_plot_data.limits = [min(ctrl_data.y) max(ctrl_data.y)];

% Control lives on a reduced mesh, get reduced mesh
[reduced_control_elements,~] = get_reduced_mesh(FOM.MESH,FOM.control_basis_index);
ctrl_data.reduced.vertices     = FOM.MESH.vertices(:,FOM.control_basis_index);
ctrl_data.reduced.elements     = reduced_control_elements; 
ctrl_data.reduced.indexes      = 1:length(FOM.control_basis_index);

fig(length(fig)+1)  = figure;
[fig] = plot_field(fig,ctrl_data,ctrl_plot_data,fonts_data);


alpha_e_2 = - v_opt ./ ( mu_0 + u_opt - lambda_2 );
beta_e_2 = ones( length( alpha_e_1 ) , 1 );
norm_2 = sqrt( beta_e_2 .^ 2 + alpha_e_2 .^ 2 );
alpha_e_2_norm = alpha_e_2 ./ norm_2;
beta_e_2_norm = beta_e_2 ./ norm_2;
v_2 = [ alpha_e_2 ./ norm_2 ; beta_e_2 ./ norm_2 ]; 
ang_2 = rad2deg( atan( beta_e_2_norm ./ alpha_e_2_norm ) );
fprintf( 'scalar product: %f' , v_1' * v_2 )
x_vert_cont = vert_cont( 1 , 1 : 1 : end )';
y_vert_cont = vert_cont( 2 , 1 : 1 : end )';
ctrl_data.name = "passive_control_u_fom";
ctrl_data.y    =ang_2;
ctrl_data.mesh = FOM.MESH;

ctrl_plot_data.title = "$\angle(\textbf{w}_{2},x_{1}) [^{\circ}]$";
ctrl_plot_data.limits = [min(ctrl_data.y) max(ctrl_data.y)];

% Control lives on a reduced mesh, get reduced mesh
[reduced_control_elements,~] = get_reduced_mesh(FOM.MESH,FOM.control_basis_index);
ctrl_data.reduced.vertices     = FOM.MESH.vertices(:,FOM.control_basis_index);
ctrl_data.reduced.elements     = reduced_control_elements; 
ctrl_data.reduced.indexes      = 1:length(FOM.control_basis_index);

fig(length(fig)+1)  = figure;
plot_field(fig,ctrl_data,ctrl_plot_data,fonts_data);


%%
% Setup figures
fig = gobjects(0);
set(0,'DefaultFigureVisible','on');


% Control field


ctrl_data.name = "passive_control_f_fom";
ctrl_data.y    = f_opt;
ctrl_data.mesh = FOM.MESH;

ctrl_plot_data.title = "Control f";
ctrl_plot_data.limits = [min(ctrl_data.y) max(ctrl_data.y)];

% Control lives on a reduced mesh, get reduced mesh
[reduced_control_elements,~] = get_reduced_mesh(FOM.MESH,FOM.control_basis_index);
ctrl_data.reduced.vertices     = FOM.MESH.vertices(:,FOM.control_basis_index);
ctrl_data.reduced.elements     = reduced_control_elements; 
ctrl_data.reduced.indexes      = 1:length(FOM.control_basis_index);

fig(length(fig)+1)  = figure;
[fig] = plot_field(fig,ctrl_data,ctrl_plot_data,fonts_data);





ctrl_data.name = "passive_control_u_fom";
ctrl_data.y    = u_opt;
ctrl_data.mesh = FOM.MESH;

ctrl_plot_data.title = "Control u ";
ctrl_plot_data.limits = [min(ctrl_data.y) max(ctrl_data.y)];

% Control lives on a reduced mesh, get reduced mesh
[reduced_control_elements,~] = get_reduced_mesh(FOM.MESH,FOM.control_basis_index);
ctrl_data.reduced.vertices     = FOM.MESH.vertices(:,FOM.control_basis_index);
ctrl_data.reduced.elements     = reduced_control_elements; 
ctrl_data.reduced.indexes      = 1:length(FOM.control_basis_index);

fig(length(fig)+1)  = figure;
[fig] = plot_field(fig,ctrl_data,ctrl_plot_data,fonts_data);


ctrl_data.name = "passive_control_v_fom";
ctrl_data.y    = v_opt;
ctrl_data.mesh = FOM.MESH;

ctrl_plot_data.title = "Control v ";
ctrl_plot_data.limits = [min(ctrl_data.y) max(ctrl_data.y)];

% Control lives on a reduced mesh, get reduced mesh
[reduced_control_elements,~] = get_reduced_mesh(FOM.MESH,FOM.control_basis_index);
ctrl_data.reduced.vertices     = FOM.MESH.vertices(:,FOM.control_basis_index);
ctrl_data.reduced.elements     = reduced_control_elements; 
ctrl_data.reduced.indexes      = 1:length(FOM.control_basis_index);

fig(length(fig)+1)  = figure;
plot_field(fig,ctrl_data,ctrl_plot_data,fonts_data);

%%
% Setup figures
fig = gobjects(0);
set(0,'DefaultFigureVisible','on');



% Scattered field
ref_data.name = "reference";
ref_data.y    = full(z);
ref_data.mesh = FOM.MESH;

ref_plot_data.limits = [min(ref_data.y) max(ref_data.y)];
ref_plot_data.title  = "Reference";
fig(length(fig)+1)  = figure;
[fig] = plot_field(fig,ref_data,ref_plot_data,fonts_data);



sc_data.name = "uncloaked";
yf = zeros(Nz,1);
yf(FOM.nodes_ocp_in,1) = y;
%yf( FOM.observation_basis_index,1 ) = FOM.E_obs * y;
yf( boundary_dof_ocp ) = T_dir;
sc_data.y    = full(yf);
sc_data.mesh = FOM.MESH;

[state_elements,state_boundaries] = get_reduced_mesh(FOM.MESH,FOM.nodes_ocp);
sc_data.reduced.vertices     = FOM.MESH.vertices(:,FOM.nodes_ocp);
sc_data.reduced.elements     = state_elements; 
sc_data.reduced.indexes      = FOM.nodes_ocp;

sc_plot_data.limits = [min(z) max(z)];
sc_plot_data.title  = "Uncloaked";
fig(length(fig)+1)  = figure;
[fig] = plot_field(fig,sc_data,sc_plot_data,fonts_data);


cl_data      = sc_data;
cl_data.name = "cloaked";
yf = zeros(Nz,1);
yf(FOM.nodes_ocp_in,1) = y_opt;
%yf( FOM.observation_basis_index,1 ) = FOM.E_obs * y_opt;
cl_data.y    = full(yf);

cl_plot_data = sc_plot_data;
cl_plot_data.title  = "Cloaked";
fig(length(fig)+1)  = figure;
[fig] = plot_field(fig,cl_data,cl_plot_data,fonts_data);



%% Tracking error field 

tracking_difference = FOM.E_obs * y_opt - FOM.E_obs * param.E * param.z;
 
yf = zeros(Nz,1);
yf( FOM.observation_basis_index,1 ) = abs(tracking_difference);
sc_data.y    = full(yf);
sc_data.mesh = FOM.MESH;

[state_elements,~] = get_reduced_mesh(FOM.MESH,FOM.nodes_ocp);
sc_data.reduced.vertices     = FOM.MESH.vertices(:,FOM.nodes_ocp);
sc_data.reduced.elements     = state_elements; 
sc_data.reduced.indexes      = FOM.nodes_ocp;

min_plot = min( min( tracking_difference ) , min( tracking_difference )  );
max_plot = max( max( tracking_difference ) , max( tracking_difference ) );
max_max = max( abs( min_plot ) , max_plot );
sc_plot_data.limits = [-max_max , max_max];
sc_plot_data.title  = "tracking error";
 fig(length(fig)+1)  = figure;
 [fig] = plot_field(fig,sc_data,sc_plot_data,fonts_data);
  hold on 

switch data_name
    case{'mesh_data\data_setup_coarse_1_refinement_2_ufv'}
        
        r = 0.2;
        th = 0:pi/50:2*pi;
        x = 0;
        y=0;
        x_circle = r * cos(th) + x;
        y_circle = r * sin(th) + y;
        circles = plot(x_circle, y_circle);
        fill(x_circle, y_circle, 'white')

        r=0.2;
        R=0.4;
        xf = 0;
        yf=0;
        Xf=0;
        Yf=0;
        t = linspace(0,2*pi,200);
        x = xf + r*cos(t);
        y = yf + r*sin(t);
        X = Xf + R*cos(t);
        Y = Yf + R*sin(t);
        fill(X,Y, 'black');
        fill(x,y, 'white');
        
        % L(1) = line(x,y,'color','w');
        % L(2) = line(X,Y,'color','w');

        
    case{'mesh_data\data_setup_coarse_boar_cloak_shape'}
        
        load('mesh_boar_cloak')
        get_grey_obs_scrofa( to_save.shape.outer_vert , to_save.shape.obs_vert )  
        hold off
        
end



axis equal





