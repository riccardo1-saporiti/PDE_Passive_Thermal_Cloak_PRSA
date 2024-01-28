clc
clearvars
close all

%Plot Figure 15

if ispc
    sslash = '\';
elseif isunix
    sslash = '/';
end

addpath("useful_functions" , "results_time_dependent");


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


%% Assemble matrices
data_name = "data_setup_coarse_grid";             
data_name = strcat("mesh_data",sslash,data_name);    
load(data_name);

B_dd_u = FOM.B_dd_u;
B_dd_f = FOM.B_dd_f;
B_dl = FOM.B_dl;
param.B_dd_u = B_dd_u;
param.B_dd_f = B_dd_f;
param.B_dl = B_dl;
A = mu_0*FOM.A_d_ocp+FOM.A_d_robin_ocp;  
A_z = FOM.A_d + FOM.A_d_robin;
F_z = FOM.F;
F = FOM.F_ocp;
E = FOM.E;
M = FOM.M;
M_ocp = FOM.M_ocp;
param.M_ocp = M_ocp;
boundary_dof_ocp = FOM.boundary_dof_ocp;
beta        = 1e-9;     
param.beta  = beta;
beta_g = 1e-5;   
xi = 1e-9;
xi_g = 1e-5;
gamma = 1e-9;
gamma_g = 1e-5;
param.beta_g = beta_g;
param.xi = xi;
param.xi_g = xi_g;
param.gamma = gamma;
param.gamma_g = gamma_g;
param.A     = A;
param.M_o   = FOM.M_obs;         
M_o         = param.M_o;
param.M_u   = FOM.M_u;
M_u         = param.M_u;
param.A_u   = FOM.A_u;
A_u         = param.A_u;
param.E = E;
param.F     = F;
param.I_s   = I_s;
param.T_dir = T_dir;
mu_min  = 0 + eps; 
param.mu_0 = mu_0;
param.mu_min = mu_min;
tol_constraint = 1e-7;
param.tol_constraint = tol_constraint;

t0 = 0;
tf = 2;
Nt = 14;
param.Nt = Nt;
dt = ( tf - t0 ) / Nt;
param.t0 = t0;
param.tf = tf;
param.dt = dt;
mu_min = 1e-1;
param.mu_min = mu_min;
theta = 1;
param.theta = theta;
M_LHS_z = ( M / dt ) + theta * A_z;
M_RHS_z = ( M / dt ) + ( theta - 1 ) * A_z;
z_ss = A_z \ ( I_s * F_z ); 
Nz = length( z_ss );
param.Nz = Nz;
z_times = zeros( Nz , 1 );
z_vect = zeros( Nz * ( Nt + 1 ) , 1 );
z_i = z_times( : , 1 );
z_vect( 1 : Nz , 1 ) = z_i; 
Ny = size( A , 1 );
z_vect_cost = zeros( Ny * ( Nt + 1 ) , 1 );
z_vect_cost( 1 : Ny , 1 ) = E * z_i;
Nu = size( B_dl , 3 );
param.Nu = Nu;
param.Ny = Ny;
M_LHS_y = ( M_ocp / dt ) + A * theta;
M_RHS_y = ( M_ocp / dt ) + A * ( theta - 1 );
y_unc = ( A ) \ ( F * I_s ) ;   
y_unc_times = zeros( Ny , 1 );
y_vect = zeros( Ny * ( Nt + 1 ) , 1 );
y_i = y_unc_times( : , 1 );
y_vect( 1 : Ny , 1 ) = y_i;
tt = 1;


for t = ( t0 + dt ) : dt : tf 

    RHS_z = M_RHS_z * z_i + I_s * F_z;
    RHS_y = M_RHS_y * y_i + I_s * F;
    z_new = M_LHS_z \ RHS_z;
    z_i = z_new;
    y_i = M_LHS_y \ RHS_y;
    z_times( : , tt + 1 ) = z_new;
    z_vect( Nz * tt + 1 : Nz * ( tt + 1 ) ) = z_new;
    z_vect_cost( Ny * tt + 1 : Ny * ( tt + 1 ) ) = E * z_new;
    y_unc_times( : , tt + 1 ) = y_i;
    y_vect( Ny * tt + 1 : Ny * ( tt + 1 ) ) = y_i;
    tt = tt + 1;

end

param.z_vect = z_vect;
param.z_vect_cost = z_vect_cost;
param.z_times = z_times;


y_opt = zeros( ( Nt + 1 ) * Ny , 1 );
p_opt = zeros( ( Nt + 1 ) * Ny , 1 );
y_old = y_opt( 1 : Ny );          %initial condition on the state

load( 'output_ufv_stationary.mat' , 'u_opt','f_opt','v_opt' )
u_opt = [ zeros( Nu , 1 ) ; repmat( u_opt , ( Nt ) , 1 ) ];
f_opt = [ zeros( Nu , 1 ) ; repmat( f_opt , ( Nt ) , 1 ) ];
v_opt = [ zeros( Nu , 1 ) ; repmat( v_opt , ( Nt ) , 1 ) ];


lext =length( u_opt( Nu + 1 : end ) );


ii = 1; 
B_u_old = sparse( double( ttv( B_dd_u , u_opt( 1 : Nu ) , 3 ) ) );
B_f_old = sparse( double( ttv( B_dd_f , f_opt( 1 : Nu ) , 3 ) ) );
B_v_old = sparse( double( ttv( B_dl , v_opt( 1 : Nu ) , 3 ) ) );

for t = ( t0 + dt ) : dt : tf

        B_u_new = sparse( double( ttv( B_dd_u , u_opt( Nu * ii + 1 :  Nu * ( ii + 1 ) ) , 3 ) ) );
        B_f_new = sparse( double( ttv( B_dd_f , f_opt( Nu * ii + 1 : Nu * ( ii + 1 ) ) , 3 ) ) );
        B_v_new = sparse( double( ttv( B_dl , v_opt( Nu * ii + 1 :  Nu * ( ii + 1 )  ) , 3 ) ) );
        M_state_LHS = ( M_ocp / dt + theta * ( A + B_u_new + B_f_new + B_v_new ) );
        M_state_RHS = ( M_ocp / dt + ( theta - 1 ) * ( A + B_u_old + B_f_old + B_v_old ) );
        RHS_y = M_state_RHS * y_old + I_s * F;
        y_opt( Ny * ii + 1 :  Ny * ( ii + 1 )  )     = M_state_LHS \ RHS_y;
        y_old = y_opt(  Ny * ii + 1 :  Ny * ( ii + 1 )  );
        B_u_old = B_u_new;
        B_f_old = B_f_new;
        B_v_old = B_v_new;
        ii = ii + 1;
                                
end

       


u_times = zeros( Nu , 1 );
f_times = zeros( Nu , 1 );
v_times = zeros( Nu , 1 );
y_times = zeros( Ny , 1 );
lambda_1_times = zeros( Nu , 1 );
lambda_2_times = zeros( Nu , 1 );

for ii = 1 : Nt + 1

        u_times( : , ii ) = u_opt( Nu * ( ii - 1 ) + 1 : Nu * ii );
        f_times( : , ii ) = f_opt( Nu * ( ii - 1 ) + 1 : Nu * ii );
        v_times(:,ii)=v_opt(Nu*(ii-1)+1:Nu*ii );
        y_times( : , ii ) = y_opt( Ny * ( ii - 1 ) + 1 : Ny * ii );
        delta = ( u_times( : , ii ) - f_times( : , ii ) ) .^ 2 + 4 * v_times( : , ii ) .^ 2;
        lambda_1_times( : , ii ) =  ( 2 * mu_0 + u_times( : , ii ) + f_times( : , ii ) + sqrt( delta ) ) * 0.5;
        lambda_2_times( : , ii ) = ( 2 * mu_0 + u_times( : , ii ) + f_times( : , ii ) - sqrt( delta ) ) * 0.5;
end


%%
if min(min(lambda_1_times)) > 0 
    disp( 'The first Eigenvalue is positive' )
end
if min(min(lambda_2_times)) > 0 
    disp( 'The second Eigenvalue is positive' )
end


%%
%Select with i_times the time instant at which display the plot, select
%from i_times = 2 to 15, since the controls are initialized with a value
%equal to 0 and the plot can not be displayed in such a case

fig = gobjects(0);
set(0,'DefaultFigureVisible','on');
i_times = 15;
ctrl_data.name = "passive_control_u_fom";
ctrl_data.y    = full( u_times );
ctrl_data.mesh = FOM.MESH;

ctrl_plot_data.title = "Control u ";
%ctrl_plot_data.limits = [min( u_times( : , i_times ) )  max( u_times( : , i_times ) ) ];
ctrl_plot_data.limits = [min( min( u_times ) )  max( max( u_times ) ) ];
ctrl_plot_data.tt = i_times;
% Control lives on a reduced mesh, get reduced mesh
ctrl_plot_data.dt = dt;
[reduced_control_elements,~] = get_reduced_mesh(FOM.MESH,FOM.control_basis_index);
ctrl_data.reduced.vertices     = FOM.MESH.vertices(:,FOM.control_basis_index);
ctrl_data.reduced.elements     = reduced_control_elements; 
ctrl_data.reduced.indexes      = 1:length(FOM.control_basis_index);
fig(length(fig)+1)  = figure;
plot_field_td(fig,ctrl_data,ctrl_plot_data,fonts_data);

%%

fig = gobjects(0);
set(0,'DefaultFigureVisible','on');
i_times = 15;
ctrl_data.name = "passive_control_u_fom";
ctrl_data.y    = full( f_times );
ctrl_data.mesh = FOM.MESH;

ctrl_plot_data.title = "Control f ";
%ctrl_plot_data.limits = [min( f_times( : , i_times ) )  max( f_times( : , i_times ) ) ];
ctrl_plot_data.limits = [min( min( f_times ) )  max( max( f_times ) ) ];
ctrl_plot_data.tt = i_times;
% Control lives on a reduced mesh, get reduced mesh
ctrl_plot_data.dt = dt;
[reduced_control_elements,~] = get_reduced_mesh(FOM.MESH,FOM.control_basis_index);
ctrl_data.reduced.vertices     = FOM.MESH.vertices(:,FOM.control_basis_index);
ctrl_data.reduced.elements     = reduced_control_elements; 
ctrl_data.reduced.indexes      = 1:length(FOM.control_basis_index);
fig(length(fig)+1)  = figure;
plot_field_td(fig,ctrl_data,ctrl_plot_data,fonts_data);

%%

fig = gobjects(0);
set(0,'DefaultFigureVisible','on');
i_times = 15;
ctrl_data.name = "passive_control_u_fom";
ctrl_data.y    = full( v_times );
ctrl_data.mesh = FOM.MESH;

ctrl_plot_data.title = "Control v ";
%ctrl_plot_data.limits = [min( v_times( : , i_times ) )  max( v_times( : , i_times ) ) ];
ctrl_plot_data.limits = [min( min( v_times ) )  max( max( v_times ) ) ];
ctrl_plot_data.tt = i_times;
% Control lives on a reduced mesh, get reduced mesh
ctrl_plot_data.dt = dt;
[reduced_control_elements,~] = get_reduced_mesh(FOM.MESH,FOM.control_basis_index);
ctrl_data.reduced.vertices     = FOM.MESH.vertices(:,FOM.control_basis_index);
ctrl_data.reduced.elements     = reduced_control_elements; 
ctrl_data.reduced.indexes      = 1:length(FOM.control_basis_index);
fig(length(fig)+1)  = figure;
plot_field_td(fig,ctrl_data,ctrl_plot_data,fonts_data);

% % 

%% Plot of the eigenvalues 


fig = gobjects(0);
set(0,'DefaultFigureVisible','on');
i_times = 15;
ctrl_data.name = "passive_control_u_fom";
ctrl_data.y    = full( lambda_1_times );
ctrl_data.mesh = FOM.MESH;

ctrl_plot_data.title = "$\lambda_{1}$";
ctrl_plot_data.limits = [min( min( lambda_1_times ) )  max( max( lambda_1_times  ) ) ];
ctrl_plot_data.tt = i_times;
% Control lives on a reduced mesh, get reduced mesh
ctrl_plot_data.dt = dt;
[reduced_control_elements,~] = get_reduced_mesh(FOM.MESH,FOM.control_basis_index);
ctrl_data.reduced.vertices     = FOM.MESH.vertices(:,FOM.control_basis_index);
ctrl_data.reduced.elements     = reduced_control_elements; 
ctrl_data.reduced.indexes      = 1:length(FOM.control_basis_index);
fig(length(fig)+1)  = figure;
plot_field_td(fig,ctrl_data,ctrl_plot_data,fonts_data);


%%

Ny = size( y_times , 1 );
i_times = 15;
fig = gobjects(0);
set(0,'DefaultFigureVisible','on');

sc_data.name = "cloaked";
yf = zeros( size( y_times ) );
yf( boundary_dof_ocp , : ) = T_dir;
yf( FOM.nodes_ocp_in , i_times ) = y_opt( Ny * ( i_times - 1 ) + 1 : Ny * i_times );
%yf( FOM.observation_basis_index, i_times ) = FOM.E_obs * y_times( : , i_times );
sc_data.y    = yf;
sc_data.mesh = FOM.MESH;

[state_elements,~] = get_reduced_mesh(FOM.MESH,FOM.nodes_ocp);
sc_data.reduced.vertices     = FOM.MESH.vertices(:,FOM.nodes_ocp);
sc_data.reduced.elements     = state_elements; 
sc_data.reduced.indexes      = FOM.nodes_ocp;

sc_plot_data.limits = [min(z_ss) max(z_ss)];
sc_plot_data.limits = [min(z_ss) , max(z_ss)];
sc_plot_data.title  = "Cloaked";
sc_plot_data.tt = i_times;
% Control lives on a reduced mesh, get reduced mesh
sc_plot_data.dt = dt;
fig(length(fig)+1)  = figure;
plot_field_td(fig,sc_data,sc_plot_data,fonts_data);





fig = gobjects(0);
set(0,'DefaultFigureVisible','on');



ref_data.name = "reference";
ref_data.y    = full(z_times);
ref_data.mesh = FOM.MESH;

ref_plot_data.limits = [min(min(z_ss)) max(max(z_ss))];
ref_plot_data.title  = "Reference";
ref_plot_data.tt = i_times;
ref_plot_data.dt = dt;
fig(length(fig)+1)  = figure;
 plot_field_td(fig,ref_data,ref_plot_data,fonts_data);

%%

i_times = 10;
fig = gobjects(0);
set(0,'DefaultFigureVisible','on');


sc_data.name = "uncloaked";
yf = zeros( size( y_unc_times ) );
yf( boundary_dof_ocp , : ) = T_dir;
yf( FOM.nodes_ocp_in , i_times ) = y_unc_times(:,i_times);
%yf( FOM.observation_basis_index, i_times ) = FOM.E_obs * y_times( : , i_times );
sc_data.y    = yf;
sc_data.mesh = FOM.MESH;

[state_elements,state_boundaries] = get_reduced_mesh(FOM.MESH,FOM.nodes_ocp);
sc_data.reduced.vertices     = FOM.MESH.vertices(:,FOM.nodes_ocp);
sc_data.reduced.elements     = state_elements; 
sc_data.reduced.indexes      = FOM.nodes_ocp;

sc_plot_data.limits = [min(z_ss) max(z_ss)];
sc_plot_data.limits = [min(z_ss) , max(z_ss)];
sc_plot_data.title  = "Uncloaked";
sc_plot_data.tt = i_times;
% Control lives on a reduced mesh, get reduced mesh
sc_plot_data.dt = dt;
fig(length(fig)+1)  = figure;
[fig] = plot_field_td(fig,sc_data,sc_plot_data,fonts_data);

%%
%Plot of the first eigenvectors
i_times = 5;
u_opt = u_times( : , i_times );
f_opt = f_times( : , i_times );
v_opt = v_times( : , i_times );
lambda_1 = lambda_1_times( : , i_times );
lambda_2 = lambda_2_times( : , i_times );
alpha_e_1 = - v_opt ./ ( mu_0 + u_opt - lambda_1 );
beta_e_1 = ones( length( alpha_e_1 ) , 1 );
norm_1 = sqrt( beta_e_1 .^ 2 + alpha_e_1 .^ 2 );
alpha_e_1_norm = alpha_e_1 ./ norm_1;
beta_e_1_norm = beta_e_1 ./ norm_1;
v_1 = [ alpha_e_1 ./ norm_1 ; beta_e_1 ./ norm_1 ] ;
ang_1 = rad2deg( atan( beta_e_1_norm ./ alpha_e_1_norm ) );
vert_cont = FOM.vertices_control;
x_vert_cont = vert_cont( 1 , 1 : 1 : end )';
y_vert_cont = vert_cont( 2 , 1 : 1 : end )';
quiver( x_vert_cont , y_vert_cont , alpha_e_1_norm( 1 : 1 : end , 1 ) , beta_e_1_norm( 1 : 1 : end , 1 ) )
title( '$Eigenvector\:w_{1}$' , 'Interpreter' , 'Latex' )




%%
%tracking error
y_track = y_times( : , end );
z_track = z_times( : , end );
num_mte_star = sum( ( FOM.E_obs * y_track - FOM.E_obs * param.E * z_track ) .^ 2 );
num_mte = sum( ( FOM.E_obs * y_unc_times( : , end ) - FOM.E_obs * param.E * z_track  ) .^ 2 );
den_mte = size( FOM.E_obs , 1 );
MTE_star = sqrt( num_mte_star / den_mte );
MTE = sqrt( num_mte / den_mte );
eta = abs( MTE - MTE_star ) / ( MTE );
tracking_difference = FOM.E_obs * y_track - FOM.E_obs * param.E * z_track;
sc_data.name = "tracking error Reference vs Cloaked";
yf = zeros(Nz,1);
yf( FOM.observation_basis_index,1 ) = tracking_difference;
sc_data.y    = full(yf);
sc_data.mesh = FOM.MESH;

[state_elements,~] = get_reduced_mesh(FOM.MESH,FOM.nodes_ocp);
sc_data.reduced.vertices     = FOM.MESH.vertices(:,FOM.nodes_ocp);
sc_data.reduced.elements     = state_elements; 
sc_data.reduced.indexes      = FOM.nodes_ocp;

sc_plot_data.limits = [min(tracking_difference) max(tracking_difference)];
sc_plot_data.title  = "tracking error";
sc_plot_data.tt = 1;
 fig(length(fig)+1)  = figure;
 [fig] = plot_field_td(fig,sc_data,sc_plot_data,fonts_data);

%% Tracking error field 
i_times = 15;
y_times_st = y_times; 

load('output_ufv_time_depependent.mat' , 'y_times' )
y_times_td = y_times ;

tracking_difference = FOM.E_obs * y_times_td( : , i_times ) - FOM.E_obs * param.E * z_times( : , i_times );
tracking_difference_st = FOM.E_obs * y_times_st( : , i_times ) - FOM.E_obs * param.E  * z_times( : , i_times );

yf = zeros( size( y_times_td ) );
yf( FOM.observation_basis_index, i_times ) = abs(tracking_difference_st);
sc_data.y    = full(yf);
sc_data.mesh = FOM.MESH;

[state_elements,~] = get_reduced_mesh(FOM.MESH,FOM.nodes_ocp);
sc_data.reduced.vertices     = FOM.MESH.vertices(:,FOM.nodes_ocp);
sc_data.reduced.elements     = state_elements; 
sc_data.reduced.indexes      = FOM.nodes_ocp;

min_plot = min( min( tracking_difference ) , min( tracking_difference_st )  );
max_plot = max( max( tracking_difference ) , max( tracking_difference_st ) );
max_max = max( abs( min_plot ) , max_plot );
sc_plot_data.limits = [-max_max , max_max];
sc_plot_data.title  = "tracking error";
% Control lives on a reduced mesh, get reduced mesh
sc_plot_data.tt = i_times;
sc_plot_data.dt = dt;
fig(length(fig)+1)  = figure;
plot_field_td(fig,sc_data,sc_plot_data,fonts_data);

  hold on 
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
axis equal





 