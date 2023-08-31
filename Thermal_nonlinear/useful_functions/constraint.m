function [ c , ceq , g , grad ] = constraint( x , param )

mu_0 = param.mu_0;
tol_constraint = param.tol_constraint;


Nu = size( param.B_dl , 3  );
ceq = [];
c = zeros( Nu , 1 );
c( 1 : Nu ) = x( 2 * Nu + 1 : end , 1 ) .^ 2 - ( mu_0 + x( 1 : Nu , 1 ) ) .* ( mu_0 + x( Nu + 1 : 2 * Nu , 1 ) ) + tol_constraint;
grad = [];
g = zeros( 3 * Nu , Nu );
der_con_u = - ( mu_0 + x( Nu + 1 : 2 * Nu , 1 ) );
g( 1 : Nu , 1 : Nu ) = diag( der_con_u );
der_con_f = - ( mu_0 + x( 1 : Nu , 1 ) );
g( Nu + 1 : 2 * Nu , 1 : Nu ) = diag( der_con_f );
der_con_v = + 2 * x( 2 * Nu + 1 : end , 1 );
g( 2 * Nu + 1 : end , 1 : Nu ) = diag( der_con_v );




return 
