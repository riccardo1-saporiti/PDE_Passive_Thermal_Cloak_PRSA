function [ J , grad_x ] = cost_with_grad_optimization( x , param )

beta   = param.beta;
beta_g = param.beta_g;
A      = param.A;
B_dd_u   = param.B_dd_u;
M_o    = param.M_o;
M_u    = param.M_u;
A_u    = param.A_u;
F      = param.F;
I_s    = param.I_s;
z      = param.z;
Nu   = size(B_dd_u,3);
E_z_cost = param.E;


B_dd_f = param.B_dd_f;
xi = param.xi;
xi_g = param.xi_g;
B_dl = param.B_dl;
gamma = param.gamma;
gamma_g = param.gamma_g;

u    = x(1:Nu,1);
f     = x( Nu + 1 : 2 * Nu , 1 );
v    = x( 2 * Nu + 1 : end , 1 );

B_u = sparse( double( ttv( B_dd_u , u , 3 ) ) );
B_f = sparse( double( ttv( B_dd_f , f , 3 ) ) );
B_v = sparse( double( ttv( B_dl , v , 3 ) ) );

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

y = ( A + B_u + B_f + B_v ) \ ( F * I_s );

J      = 0.5*((y-E_z_cost*z)'*M_o*(y-E_z_cost*z) + beta * u' * ( M_u ) * u + beta_g * u' * ( A_u ) * u + xi * f' * ( M_u ) * f + ...
                xi_g * f' * ( A_u ) * f + gamma * v' * ( M_u ) * v + gamma_g * v' * ( M_u ) * v );

p = ( A + B_u + B_f + B_v )' \ ( M_o * ( y - E_z_cost * z ) );
grad_u_old = ( beta * M_u * u + beta_g * A_u * u - sparse( double( ttv( ttv( B_dd_u , y , 1 ) , p , 1 ) ) ) );
grad_f_old = ( xi * M_u * f + xi_g * A_u * f - sparse( double( ttv( ttv( B_dd_f , y , 1 ) , p ,1 ) ) ) );
grad_v_old = ( gamma * M_u * v + gamma_g * A_u * v - sparse( double( ttv( ttv( B_dl , y , 1 ) , p , 1 ) ) ) );

grad_x = [ grad_u_old ; grad_f_old  ; grad_v_old ];

                


return 