%% Assemble finite element matrices 

% Reference   domain without obstacle
% Control     domain with obstacle
% Observation domain outside of the cloak  


% Necessary imports
clear all
clc
close all

if ispc
    sslash = '\';
elseif isunix
    sslash = '/';
end

addpath("useful_functions");

% Load mesh as data file with vertices, boundaries alements and labels for
% the various domains

mesh_name = "coarse_grid";            %the first number is referred to the number of refinements, the second to the location of the control

set(0,'DefaultFigureVisible','on');
load(strcat("mesh_data",sslash,"mesh_obstacle_coarse_grid")); 

vertices         = mesh_pet.p;    %2xNumNodes matrix with nodal coordinates
boundaries       = mesh_pet.e;     
elements         = mesh_pet.t;

control_dom_labels = mesh_pet.control_dom_labels;
obs_dom_labels     = mesh_pet.obs_dom_labels;      %1
ocp_dom_labels     = mesh_pet.ocp_dom_labels;      %1,2,3
boundary_ocp_labels= mesh_pet.boundary_ocp_labels;       %5,6,7,8: are the ones around obstacle
outer_boundary_labels = mesh_pet.outer_boundary_labels;  %1,2,3,4: are the most external ones

pdeplot(vertices,boundaries,elements);                       % Check that everything makes sense

% edges i,  index_vertices  e(1,i) --- e(2,i) , e(3,i) edge measure,
% e(4,i)?   e(5,i) edge label , e(6,i) edge external, internal , e(7,i) all
% ones


% Assemble matrices full domain

dim_problem  = 2;
fem = 'P1';
quad_order = 4;

DATA       = read_DataFile('datafile_empty', dim_problem);
DATA.param = [];

DATA.flag_robin   = outer_boundary_labels; %1,2,3,4
DATA.bcRob_alpha    = @(x,y,t,param)(1+ 0.*x.*y);
 
mu_diff        = 1;
DATA.diffusion = @(x,y,t,param)(mu_diff + 0.*x.*y);

%% Create and fill the MESH data structure
[ MESH ] = buildMESH( dim_problem, elements, vertices, boundaries, fem, quad_order, DATA );

% OCP mesh contains subdomain 1 2 3
elements_ocp  = [];
for dom_label = ocp_dom_labels
   
    elements_ocp = unique([elements_ocp find(elements(4,:)==dom_label)]);
    
end

nodes_ocp     = unique( [ elements(1,elements_ocp) elements(2,elements_ocp) elements(3,elements_ocp) ] ) ;

boundary_dof_ocp = 0;
for bound_label = unique(boundary_ocp_labels)
    boundary_dof_ocp = boundary_dof_ocp + (boundaries(5,:)==bound_label);
end
boundary_dof_ocp = find(boundary_dof_ocp==1);
boundary_dof_ocp = unique([boundaries(1,boundary_dof_ocp) boundaries(2,boundary_dof_ocp)]);
 

% Remove internal boundary nodes from nodes to keep
nodes_ocp_in  = setdiff(nodes_ocp,boundary_dof_ocp);


vertices_ocp = vertices(:,nodes_ocp);
vertices_ocp_in = vertices( : , nodes_ocp_in );
[elements_ocp,boundaries_ocp] = get_reduced_mesh(MESH,nodes_ocp);

MESH_OCP.vertices   = vertices_ocp;
MESH_OCP.elements   = elements_ocp;
MESH_OCP.boundaries = boundaries_ocp;

pdeplot(vertices_ocp,boundaries_ocp,elements_ocp)

hold off

%%
scatter(vertices_ocp(1,:),vertices_ocp(2,:))

scatter(vertices(1,boundary_dof_ocp),vertices(2,boundary_dof_ocp))

%scatter(vertices(1,nodes_to_keep_in),vertices(2,nodes_to_keep_in),'r','MarkerSize',2)

%nodes_to_keep_internal    = unique( [nodes_to_keep internal_boundary_dof_ocp] ) ;

% Extraction matrix, maps full nodes to ocp nodes
E = nodes_ocp_in' ./ [1:length(vertices)];
E(E~=1) = 0;
 
E       = sparse(E);

% Extraction matrix, maps full nodes to dirichlet nodes
E_dir           = boundary_dof_ocp' ./ [1:length(vertices)];
E_dir(E_dir~=1) = 0;
E_dir           = sparse(E_dir);


[ FE_SPACE ]   = buildFESpace( MESH, fem, 1, quad_order );
param.FE_SPACE = FE_SPACE;

DATA.force = @(x,y,t,param)( 1 + 0.*x.*y );
[~, F_dom, ~] = ADR_Assembler(MESH, DATA, FE_SPACE, [], [], [], [], 0);

%%
% Define source paramters

if strcmp( mesh_name , 'coarse_bottom_source' )

    x0 = 0 ; 
    y0 = -0.6 ; 
    
else

    x0 =  0.6;
    y0 =  0;
end
r  = 0.025;
source_strength = 1;

id_source = @(x,y) ( ((x-x0).^2 + (y-y0).^2 - r) <=0 ) * source_strength;
 
DATA.force = @(x,y,t,param)( id_source(x,y) );


%% Create and fill the FE_SPACE data structure

[A_d, F, M] = ADR_Assembler(MESH, DATA, FE_SPACE,'all');

% % % % [T_x] = ADR_Assembler(MESH, DATA, FE_SPACE, 'transport', [], [1], []);
% % % % [T_y] = ADR_Assembler(MESH, DATA, FE_SPACE, 'transport', [], [2], []);

%[M,~,A_d,~] = compute_state_input_matrices_simple(MESH, DATA, FE_SPACE); % assemble matrices
%[~, F, ~] = ADR_Assembler(MESH, DATA, FE_SPACE, [], [], [], [], 0);(can be
%used even to evaluate matrices at time t)


[A_d_robin, F, ~] =  ADR_ApplyBC(0*A_d, F, FE_SPACE, MESH, DATA, 0);  % Should save A_d (affine in mu) + A_robin (not mu dependent)

A_d_dir       = E * A_d       * transpose(E_dir);          % Use projection matrices E and E_dir to split matrix for Dirichlet boundary conditions
A_d_robin_dir = E * A_d_robin * transpose(E_dir);
A_d_OCP       = E * A_d       * transpose(E);
A_d_robin_OCP = E * A_d_robin * transpose(E);
F_OCP         = E * F;
M_OCP         = E * M         * transpose(E);


% Set T boundary dirichlet
T_dir      = 1;
T_dir_vect = T_dir * sparse(ones(length(boundary_dof_ocp),1)); 

% Assemble B control matrix

% r_control_in  = 0.25;
% r_control_ext = 0.35;
% r_vert                  = sqrt(  vertices(1,:).^2 + vertices(2,:).^2 ); 
% control_basis_index     = find(r_vert > r_control_in & r_vert < r_control_ext);
% control basis nodes comes from mesh subdomain

elements_control  = [];
for dom_label = control_dom_labels
   
    elements_control = unique([elements_control find(elements(4,:)==dom_label)]);
    
end

control_basis_index = unique( [ elements(1,elements_control) elements(2,elements_control) elements(3,elements_control) ] ) ;

 
elements_obs  = [];
for dom_label = obs_dom_labels
   
    elements_obs = unique([elements_obs find(elements(4,:)==dom_label)]);
    
end

observation_basis_index = unique( [ elements(1,elements_obs) elements(2,elements_obs) elements(3,elements_obs) ] ) ;


% change with r>0 for full observation

% Observation matrix, nodes with r > r_control_ext

% Visualize control nodes

%TODO mesh should be refined so that the control region is nice and smooth
plot(vertices(1,control_basis_index),vertices(2,control_basis_index),'ok','MarkerFaceColor','c') 
% 
if strcmp( mesh_name , 'fine_grid' ) 
    load('control_basis_index_coarse_1_refinement.mat')
    Interp_mat = zeros( length( control_basis_index ) , length( control_basis_index_coarse ) );
    find_diag_el = find( control_basis_index_coarse == control_basis_index_coarse );
    Interp_mat( find_diag_el ,find_diag_el ) = eye( length( find_diag_el ));
    extra_nodes_interp = setdiff( control_basis_index , control_basis_index_coarse );
    vertices_control = reshape( vertices( : , extra_nodes_interp ).' , length( extra_nodes_interp ) , [] , 2 );
    vertices_control_intepr = reshape( vertices( : , control_basis_index_coarse ).' , [] , length( control_basis_index_coarse ) , 2 );
    l2_dist_vertices_control = vecnorm( vertices_control - vertices_control_intepr , 2 , 3  );
    
    for iii = 1 : length( extra_nodes_interp )
        num_avg = 2;
        [ ~ , idx_small ] = mink( l2_dist_vertices_control( iii , : ) , num_avg ) ; 
        Interp_mat( iii+length( control_basis_index_coarse) , idx_small ) = 1 / num_avg;
    end
    
    FOM.Interp_mat = Interp_mat; %matrix that extends a coarse control to a fine control.
end


B_dd_u             = Bilinear_Assembler(MESH,FE_SPACE,'diffusion_i_c_1');    %tensor: Nq x Nq x Nu, with 
%Nq = numer of rows of  E, Nu by control_basis_index.
B_dd_u             = B_dd_u(nodes_ocp_in,nodes_ocp_in,control_basis_index);      % Extract necessary index from assembly
B_dd_u_dir = B_dd_u( nodes_ocp_in , boundary_dof_ocp , control_basis_index );


B_dd_f             = Bilinear_Assembler(MESH,FE_SPACE,'diffusion_i_c_2');    %tensor: Nq x Nq x Nu, with 
%Nq = numer of rows of  E, Nu by control_basis_index.
B_dd_f             = B_dd_f(nodes_ocp_in,nodes_ocp_in,control_basis_index);      % Extract necessary index from assembly
B_dd_f_dir = B_dd_f( nodes_ocp_in , boundary_dof_ocp , control_basis_index );

% Diffusion control tensor: dd isotropic; diagonal terms B_dl the opposite
B_dd             = Bilinear_Assembler(MESH,FE_SPACE,'diffusion');    %tensor: Nq x Nq x Nu, with 
%Nq = numer of rows of  E, Nu by control_basis_index.

B_dd             = B_dd(nodes_ocp_in,nodes_ocp_in,control_basis_index);      % Extract necessary index from assembly
B_dd_dir = B_dd( nodes_ocp_in , boundary_dof_ocp , control_basis_index );

B_dl             = Bilinear_Assembler(MESH,FE_SPACE,'diffusion_dl');
B_dl             = B_dl(nodes_ocp_in,nodes_ocp_in,control_basis_index);
B_dl_dir = B_dl( nodes_ocp_in , boundary_dof_ocp , control_basis_index );


M_u              = M(control_basis_index,control_basis_index);
A_u              = A_d(control_basis_index,control_basis_index);


%Control matrix 

E_con = control_basis_index' ./ [1:length(vertices)];
E_con(E_con ~= 1 )=0;
E_con = sparse( E_con );
B_con = ( E_con * M * E' )';
% Observation matrix construction

%M_obs = M(nodes_ocp_in,observation_basis_index);  % TODO cancellare, solo da un lato..

% Extraction matrix, maps ocp_nodes to observation nodes
E_obs           = observation_basis_index' ./ nodes_ocp_in;
E_obs(E_obs~=1) = 0;
E_obs           = sparse(E_obs);

vertices_control = vertices( : , control_basis_index );
%%
%vertices_ocp_in = vertices(:,nodes_to_keep_in);
plot(vertices(1,observation_basis_index),vertices(2,observation_basis_index),'ok','MarkerFaceColor','g') 

M_obs_test = sparse(size(M,1),size(M,2));
M_obs_test(observation_basis_index,observation_basis_index) = M(observation_basis_index,observation_basis_index);
M_obs_test = M_obs_test(nodes_ocp_in,nodes_ocp_in);
M_obs      = M_obs_test;


% Save necessary matrices (FOM = full order model)
FOM.DATA            = DATA;
FOM.FE_SPACE        = FE_SPACE;
FOM.MESH            = MESH;
FOM.A_d_robin       = A_d_robin;
FOM.A_d_robin_dir   = A_d_robin_dir;
FOM.A_d             = A_d;
FOM.A_d_dir         = A_d_dir;

FOM.A_d_ocp         = A_d_OCP;
FOM.A_d_robin_ocp   = A_d_robin_OCP;

FOM.M = M;
% % % % FOM.T_x = T_x;
% % % % FOM.T_y = T_y;
FOM.E               = E;
FOM.E_obs           = E_obs;
FOM.E_dir           = E_dir;
FOM.F               = F;
FOM.F_dom           = F_dom;
FOM.F_ocp           = F_OCP;
FOM.M_ref           = M_obs*E;
FOM.B_dd            = B_dd;
FOM.B_dl            = B_dl;
FOM.B_dd_dir = B_dd_dir;
FOM.B_dl_dir = B_dl_dir;
FOM.M_u             = M_u;
FOM.A_u             = A_u;
FOM.B_con       = B_con;

FOM.M_ocp           = M_OCP;
FOM.M_obs           = M_obs;
FOM.B_dd_u            = B_dd_u;
FOM.B_dd_f            = B_dd_f;
FOM.B_dl          = B_dl;
FOM.B_dd_u_dir = B_dd_u_dir;
FOM.B_dd_f_dir = B_dd_f_dir;
FOM.boundary_dof_ocp = boundary_dof_ocp;
FOM.control_basis_index     = control_basis_index; 
FOM.observation_basis_index = observation_basis_index;
FOM.nodes_ocp               = nodes_ocp;
FOM.nodes_ocp_in            = nodes_ocp_in;
FOM.vertices_ocp_in = vertices_ocp_in;
FOM.vertices_ocp = vertices_ocp;
FOM.vertices_control = vertices_control;
FOM.MESH_OCP = MESH_OCP;

% Save FOM data to be used in "main_passive.m"
data_set_name = strcat('mesh_data',sslash,'data_setup_',mesh_name);
save(data_set_name,'FOM','-v7.3');


 
% spy( M_OCP )
% spy( M_obs )