function [ c , ceq , g , grad ] = constraint( x , param )

mu_0 = param.mu_0;
tol_constraint = param.tol_constraint;
Nu = size( param.B_dl , 3  );
ceq = [];
Nt = param.Nt;
t0 = param.t0;
dt = param.dt;
tf=param.tf;
c = zeros( Nu * Nt , 1 );
grad =[];
g = zeros( 3 * Nu *Nt , Nu*Nt);
jj = 1;
for t = t0 : dt : ( tf - dt )
    u = x(Nu*(jj-1)+1:Nu*jj);
    f =x( Nu*Nt+Nu*(jj-1)+1 : Nu*Nt+Nu*jj);
    v= x( 2*Nu*Nt+Nu*(jj-1)+1 : 2*Nu*Nt+Nu*jj);
    c( Nu * ( jj - 1 ) + 1 : Nu * jj ) = v.^2-(mu_0+u).*(mu_0+f)+tol_constraint;
    der_con_u=-(mu_0+f);
    g(Nu*(jj-1)+1:Nu*jj, Nu*(jj-1)+1:Nu*jj) = diag( der_con_u);
    der_con_f=-(mu_0 + u );
    g(Nt*Nu+ Nu * ( jj - 1 ) + 1 :  Nu*Nt+Nu*jj , Nu*(jj-1)+1:Nu*jj ) = diag( der_con_f);
    der_con_v= 2*v;
    g(2*Nt*Nu+Nu*(jj-1)+1 : 2*Nu*Nt+Nu*jj , Nu * (jj-1)+1:Nu*jj) = diag( der_con_v);
    jj=jj+1;
end
return 
