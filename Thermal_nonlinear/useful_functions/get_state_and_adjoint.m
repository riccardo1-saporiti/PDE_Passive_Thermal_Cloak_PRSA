function [y,p] = get_state_and_adjoint(x,param)

B_dd_u = param.B_dd_u;
B_dd_f = param.B_dd_f;
A    = param.A;
M_o  = param.M_o;
A_u  = param.A_u;
M_u  = param.M_u;
F    = param.F;
I_s  = param.I_s;
z    = param.z;
beta = param.beta;
beta_g   = param.beta_g;
xi = param.xi;
xi_g = param.xi_g;


Nu   = size(B_dd_u,3);

E_z_cost = param.E;



gamma = param.gamma;
gamma_g = param.gamma_g;

B_dl = param.B_dl;
u    = x(1:Nu,1);
f    = x(Nu+1:2*Nu,1);
v   = x( 2 * Nu + 1 : end , 1 );

B_u = sparse(double(ttv(B_dd_u,u,3)));
B_f = sparse(double(ttv(B_dd_f,f,3)));
B_v = sparse(double(ttv(B_dl,v,3)));

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
 
y      = (A+B_u+B_v+B_f)  \ (I_s*F);
p      = transpose(A+B_u+B_v+B_f)\ (M_o*(y-E_z_cost*z));


end