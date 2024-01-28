function [ y , p] = get_adjoint_and_state( x , param )
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
Nu = param.Nu;
Ny = param.Ny;
Nt = param.Nt;
t0 = param.t0;
tf = param.tf;
z_times = param.z_times;
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

y = y_opt;
p=p_opt;


return 