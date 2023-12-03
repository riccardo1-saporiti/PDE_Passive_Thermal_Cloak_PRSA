function [ J , grad_x , y , p] = cost_with_grad_u_f_v_implicit_euler( x , param )

beta = param.beta;
beta_g = param.beta_g;
gamma = param.gamma;
gamma_g = param.gamma_g;
xi = param.xi;
xi_g = param.xi_g;
M_ocp = param.M_ocp;
B_dd_u = param.B_dd_u;
B_dd_f = param.B_dd_f;
B_dl = param.B_dl;
F = param.F;
I_s = param.I_s;
dt = param.dt;
theta = param.theta;
A = param.A;
M_o = param.M_o;
E = param.E;
A_G_cont_elem = param.A_G_cont_elem;
M_G_cont_elem = param.M_G_cont_elem;
M_G_elem = param.M_G_elem;

Nu = param.Nu;
Ny = param.Ny;
Nt = param.Nt;
t0 = param.t0;
tf = param.tf;
z_times = param.z_times;
z_vect_cost = param.z_vect_cost;
u_opt = x( 1 : Nu * ( Nt ) , 1 );
u_opt = [ zeros( Nu , 1 );u_opt];
f_opt = x(Nu * ( Nt  )+1:2*Nu * ( Nt ),1);
f_opt=[zeros(Nu,1);f_opt];
v_opt= x(2 * Nu * ( Nt ) + 1 : end , 1 );
v_opt=[zeros(Nu,1);v_opt];
y_old = zeros( Ny , 1 );
y_opt = zeros( ( Nt + 1 ) * Ny , 1 );
p_opt = zeros( ( Nt + 1 ) * Ny , 1 );

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

ii = 1;
B_u_n = sparse( double( ttv( B_dd_u , u_opt( end - Nu + 1 : end ) , 3 ) ) );
B_f_n = sparse( double( ttv( B_dd_f , f_opt( end - Nu + 1 : end ) , 3 ) ) );
B_v_n = sparse( double( ttv( B_dl , v_opt( end - Nu + 1 : end ) , 3 ) ) );
A_plus = ( M_ocp / dt ) + theta * ( A + B_u_n + B_f_n + B_v_n );
p_old = A_plus \ ( theta * M_o * ( y_opt( end - Ny + 1 : end ) - E * z_times( : , end ) ) );
p_opt( end - Ny + 1 : end , 1 ) = p_old;

for t = t0 : dt : ( tf - dt )
    y_new = y_opt( end - Ny * ( ii + 1 ) + 1 : end - Ny * ii , 1 );
    B_u_i = sparse( double( ttv( B_dd_u , u_opt( end - Nu * ( ii + 1 ) + 1 : end - Nu * ii , 1 ) , 3 ) ) );
    B_f_i = sparse( double( ttv( B_dd_f , f_opt( end - Nu * ( ii + 1 ) + 1 : end - Nu * ii , 1 ) , 3 ) ) );
    B_v_i = sparse( double( ttv( B_dl , v_opt( end - Nu * ( ii + 1 ) + 1 : end - Nu * ii , 1 ) , 3 ) ) );
    A_plus = ( M_ocp / dt ) + theta * ( A + B_u_i + B_f_i + B_v_i );
    A_minus = ( M_ocp / dt ) + ( theta - 1 ) * ( A + B_u_i + B_f_i + B_v_i );
    RHS_adj = A_minus * p_old + M_o * ( y_new - E * z_times( : , end - ii ) );
    p_opt( end - ( ii + 1 ) * Ny + 1 : end - ii * Ny ) = A_plus \ RHS_adj;
    p_old = p_opt( end - ( ii + 1 ) * Ny + 1 : end - ii * Ny );
    ii = ii + 1;
end


                   
B_vector_u = zeros( length( u_opt ) - Nu , 1 );
B_vector_f = zeros( length( u_opt ) - Nu , 1 );                   
B_vector_v = zeros( length( u_opt ) - Nu , 1 );

grad_u_elem = zeros( Nt * Nu ,1 );
grad_v_elem = grad_u_elem;
grad_f_elem = grad_v_elem;
J = 0 ;

for j = 2 : Nt + 1 

    if j < Nt + 1 
                 
        B_add_1_u = theta * sparse( double( ttv( ttv( B_dd_u , y_opt( Ny * ( j - 1 ) + 1 : Ny * j  , 1 ) , 1 ) , p_opt( Ny * ( j - 1 ) + 1 : Ny * j , 1 ) , 1 ) ) );
        B_add_2_u = - ( theta - 1 ) * sparse( double( ttv( ttv( B_dd_u , y_opt( Ny * ( j - 1 ) + 1 : Ny * j  , 1 ) , 1 ) , p_opt( Ny * j + 1 : Ny * ( j + 1 ) ) , 1 ) ) );
        B_vector_u( ( j - 2 ) * Nu + 1 : ( j - 1 ) * Nu ) = B_add_1_u + B_add_2_u;

        B_add_1_f = theta * sparse( double( ttv( ttv( B_dd_f , y_opt( Ny * ( j - 1 ) + 1 : Ny * j  , 1 ) , 1 ) , p_opt( Ny * ( j - 1 ) + 1 : Ny * j , 1 ) , 1 ) ) );
        B_add_2_f = - ( theta - 1 ) * sparse( double( ttv( ttv( B_dd_f , y_opt( Ny * ( j - 1 ) + 1 : Ny * j  , 1 ) , 1 ) , p_opt( Ny * j + 1 : Ny * ( j + 1 ) ) , 1 ) ) );
        B_vector_f( ( j - 2 ) * Nu + 1 : ( j - 1 ) * Nu ) = B_add_1_f + B_add_2_f;

        B_add_1_v = theta * sparse( double( ttv( ttv( B_dl , y_opt( Ny * ( j - 1 ) + 1 : Ny * j  , 1 ) , 1 ) , p_opt( Ny * ( j - 1 ) + 1 : Ny * j , 1 ) , 1 ) ) );
        B_add_2_v = - ( theta - 1 ) * sparse( double( ttv( ttv( B_dl , y_opt( Ny * ( j - 1 ) + 1 : Ny * j  , 1 ) , 1 ) , p_opt( Ny * j + 1 : Ny * ( j + 1 ) ) , 1 ) ) );
        B_vector_v( ( j - 2 ) * Nu + 1 : ( j - 1 ) * Nu ) = B_add_1_v + B_add_2_v;


    else

        y_ttv = y_opt( Ny * ( j - 1 ) + 1 : Ny * j , 1 );
        p_ttv = p_opt( Ny * ( j - 1 ) + 1 : Ny * j , 1 );
        B_vector_u( ( j - 2 ) * Nu + 1 : ( j - 1 ) * Nu ) = theta * sparse( double( ttv( ttv( B_dd_u , y_ttv , 1 ) , p_ttv , 1 ) ) );
        B_vector_f( ( j - 2 ) * Nu + 1 : ( j - 1 ) * Nu ) = theta * sparse( double( ttv( ttv( B_dd_f , y_ttv , 1 ) , p_ttv , 1 ) ) );
        B_vector_v( ( j - 2 ) * Nu + 1 : ( j - 1 ) * Nu ) = theta * sparse( double( ttv( ttv( B_dl , y_ttv , 1 ) , p_ttv , 1 ) ) );
        
    end


    grad_u_elem( ( j- 2 ) * Nu + 1 : ( j - 1 ) * Nu ) = ...
        beta * dt * M_G_cont_elem  * u_opt( ( j - 1 ) * Nu + 1  :  ( j  ) * Nu , 1 ) + beta_g * dt * A_G_cont_elem * u_opt(  ( j - 1) * Nu + 1  :  ( j ) * Nu , 1  ) - dt * B_vector_u(  ( j - 2 ) * Nu + 1  :  ( j  -1 ) * Nu , 1  );
    grad_f_elem( ( j- 2 ) * Nu + 1 : ( j - 1 ) * Nu ) = ...
            xi * dt * M_G_cont_elem  * f_opt( ( j - 1 ) * Nu + 1  :  ( j  ) * Nu , 1 ) + xi_g * dt * A_G_cont_elem * f_opt(  ( j - 1) * Nu + 1  :  ( j ) * Nu , 1  ) - dt * B_vector_f(  ( j - 2 ) * Nu + 1  :  ( j  -1 ) * Nu , 1  );
    grad_v_elem( ( j- 2 ) * Nu + 1 : ( j - 1 ) * Nu ) = ...
            gamma * dt * M_G_cont_elem  * v_opt( ( j - 1 ) * Nu + 1  :  ( j  ) * Nu , 1 ) + gamma_g * dt * A_G_cont_elem * v_opt(  ( j - 1) * Nu + 1  :  ( j ) * Nu , 1  ) - dt * B_vector_v(  ( j - 2 ) * Nu + 1  :  ( j  -1 ) * Nu , 1  );
    
    y_opt_iter = y_opt( ( j - 1 ) * Ny + 1 : j * Ny );
    z_vect_iter = z_vect_cost( ( j - 1 ) * Ny + 1 : j * Ny );
    u_opt_iter = u_opt(( j - 1) * Nu + 1  :  ( j ) * Nu , 1  );
    f_opt_iter = f_opt( ( j - 1) * Nu + 1  :  ( j ) * Nu , 1  );
    v_opt_iter = v_opt( ( j - 1) * Nu + 1  :  ( j ) * Nu , 1  );

    J = J + 0.5 * dt * transpose( y_opt_iter - z_vect_iter ) * M_G_elem * ( y_opt_iter - z_vect_iter ) + 0.5 * dt * beta * transpose( u_opt_iter ) * M_G_cont_elem * u_opt_iter + ... 
    0.5 * dt * beta_g * transpose( u_opt_iter ) * A_G_cont_elem * u_opt_iter + 0.5 * dt * xi * transpose( f_opt_iter ) * M_G_cont_elem * f_opt_iter + ... 
    0.5 * dt * xi_g * transpose( f_opt_iter ) * A_G_cont_elem * f_opt_iter + 0.5 * dt * gamma * transpose( v_opt_iter ) * M_G_cont_elem * v_opt_iter + ... 
    0.5 * dt * gamma_g * transpose( v_opt_iter ) * A_G_cont_elem * v_opt_iter;
end


grad_x = [ grad_u_elem ; grad_f_elem ; grad_v_elem ];



y = y_opt;
p=p_opt;


return 