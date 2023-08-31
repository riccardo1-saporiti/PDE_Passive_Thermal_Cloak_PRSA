function [ c ] = constraint_ipopt( x , param )

mu_0 = param.mu_0;
Nu = size( param.B_dl , 3  );
Nt = param.Nt;
t0 = param.t0;
dt = param.dt;
tf=param.tf;
c_v = zeros( Nu * Nt , 1 );
c_linear = zeros( Nu*Nt,1);
%grad =[];
%g = zeros( 3 * Nu *Nt , Nu*Nt);
jj = 1;
for t = t0 : dt : ( tf - dt )
    u = x(Nu*(jj-1)+1:Nu*jj);
    f =x( Nu*Nt+Nu*(jj-1)+1 : Nu*Nt+Nu*jj);
    v= x( 2*Nu*Nt+Nu*(jj-1)+1 : 2*Nu*Nt+Nu*jj);
    c_v( Nu * ( jj - 1 ) + 1 : Nu * jj ) = v.^2-(mu_0+u).*(mu_0+f);
    c_linear(Nu*(jj-1)+1:Nu*jj) = -(2*mu_0+u+f);
%     der_con_u=-(mu_0+f);
%     g(Nu*(jj-1)+1:Nu*jj, Nu*(jj-1)+1:Nu*jj) = diag( der_con_u);
%     der_con_f=-(mu_0 + u );
%     g(Nt*Nu+ Nu * ( jj - 1 ) + 1 :  Nu*Nt+Nu*jj , Nu*(jj-1)+1:Nu*jj ) = diag( der_con_f);
%     der_con_v= 2 * v;
%     g(2*Nt*Nu+Nu*(jj-1)+1 : 2*Nu*Nt+Nu*jj , Nu * (jj-1)+1:Nu*jj) = diag( der_con_v);
    jj=jj+1;
end
c= [ c_v;c_linear];

return 
