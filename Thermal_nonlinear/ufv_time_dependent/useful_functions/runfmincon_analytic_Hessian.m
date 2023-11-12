function [xsol,fval,history,grad_norm,lambda] = runfmincon_analytic_Hessian( cost_with_grad_param , constraint ,  param , z_0 , lb )

% % history.x = [];
history.fval = [];
% % searchdir = [];
grad_norm = [];

Nu = param.Nu;
% % % % Acon = [ eye( Nu , Nu ) , eye( Nu,Nu) , zeros(Nu,Nu)];
% % % % b = ( 2 *  mu_0 - tol_constraint)*ones(Nu,1);
Acon = [];
b=[];
Aeq = [];
beq = []; 
ub= [];
NL_con = constraint; 
options.HessianFcn = [];
objfun = cost_with_grad_param;
       options = optimoptions( 'fmincon','SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,'CheckGradients',false,'Display','iter', ...
           'OutputFcn',@outfun,'OptimalityTolerance', 0, ...
    'StepTolerance', 0, ...
    'MaxFunctionEvaluations', inf,...
    'MaxIterations', 1000 , 'FiniteDifferenceType' , 'central' , 'SubproblemAlgorithm' , 'cg' , 'HessianMultiplyFcn', @HessMultFcn );

[xsol,fval,~,~,lambda] = fmincon(objfun,z_0,Acon,b,Aeq,beq,lb,ub,NL_con,options);

    function stop = outfun(~,optimValues,state)
     stop = false;
 
     switch state
         case 'init'
             hold on
         case 'iter'
         % Concatenate current point and objective function
         % value with history. x must be a row vector.
           history.fval = [history.fval; optimValues.fval];
            grad_norm = [grad_norm;... 
                        optimValues.firstorderopt];
%            history.x = [history.x; x]; no need to save even values of the
%            control
         % Concatenate current search direction with 
         % searchdir.
% %            searchdir = [searchdir;... 
% %                         optimValues.searchdirection'];
% %            plot(x(1),x(2),'o');
         % Label points with iteration number and add title.
         % Add .15 to x(1) to separate label from plotted 'o'.
%            text(x(1)+.15,x(2),... 
%                 num2str(optimValues.iteration));
%            title('Sequence of Points Computed by fmincon');
         case 'done'
             hold off
         otherwise
     end
    end

    function d = HessMultFcn( x , lamdba , v )


        
      
        [y_opt ,p_opt ] = get_adjoint_and_state( x ,param );

        B_dd_u = param.B_dd_u;
        B_dd_f = param.B_dd_f;
        B_dl=param.B_dl;
        Nu = param.Nu;
        Ny = param.Ny;
        Nt = param.Nt;
        t0 =  param.t0;
        dt = param.dt;
        tf = param.tf;
%         u_opt = x( 1 : Nu , 1 );
%         f_opt = x( Nu + 1 : 2 * Nu , 1 );
%         v_opt = x( 2 * Nu + 1 : end , 1 );
        u_opt_total = x( 1 : Nu * ( Nt ) , 1 );
        u_opt_total = [ zeros( Nu , 1 );u_opt_total];
        v_u = [ zeros( Nu , 1 ) ; v( 1 : Nu * ( Nt ) , 1 ) ];
        f_opt_total = x(Nu * ( Nt  )+1:2*Nu * ( Nt ),1);
        f_opt_total=[zeros(Nu,1);f_opt_total];
        v_f = [zeros(Nu,1);v( Nu * ( Nt  )+1:2*Nu * ( Nt ),1 )];
        v_opt_total= x(2 * Nu * ( Nt ) + 1 : end , 1 );
        v_opt_total=[zeros(Nu,1);v_opt_total];
        v_v = [zeros(Nu,1);v( 2 * Nu * ( Nt ) + 1 : end , 1 )];
        act_on_v_of_u = zeros( length( u_opt_total ) , 1 );
        act_on_v_of_f = zeros( length( f_opt_total ) , 1 );
        act_on_v_of_v = zeros( length( v_opt_total ) , 1 );
        ii = 1;

for t = ( t0 + dt ) : dt : tf  
        
        u_opt = u_opt_total( Nu * ii + 1 :  Nu * ( ii + 1 ) );
        f_opt = f_opt_total( Nu * ii + 1 :  Nu * ( ii + 1 ) );
        v_opt = v_opt_total( Nu * ii + 1 :  Nu * ( ii + 1 ) );
        y = y_opt( Ny * ( ii - 1 ) + 1 : Ny * ii , 1 );
        p = p_opt( Ny * ( ii - 1 ) + 1 : Ny * ii , 1 );
        B_u = spmatrix( ttv( B_dd_u , u_opt , 3 ) );
        B_f = spmatrix( ttv( B_dd_f , f_opt , 3 ) );
        B_v = spmatrix( ttv( B_dl , v_opt , 3 ) );
        B_u_q = spmatrix( ttv( B_dd_u ,y , 1 ) );
        B_f_q = spmatrix( ttv( B_dd_f , y , 1 ) );
        B_v_q = spmatrix( ttv( B_dl , y , 1 ) );
        B_u_l = spmatrix( ttv( B_dd_u , p , 1 ) );
        B_f_l = spmatrix( ttv( B_dd_f , p , 1 ) );
        B_v_l = spmatrix( ttv( B_dl , p , 1 ) );
        
        A = param.A;
        A_z = A + B_u + B_f + B_v;
        M_o = param.M_o;
        M_u = param.M_u;
        A_u = param.A_u;
        beta = param.beta;
        beta_g   = param.beta_g;
        xi = param.xi;
        xi_g = param.xi_g;
        gamma = param.gamma;
        gamma_g = param.gamma_g;
        % Big matrix to get hessian direction
        
        
%         K = [  A_z , zeros( Ny , Ny ) , B_u_q , B_f_q , B_v_q ;
%                M_o , A_z' , -B_u_l , - B_f_l , - B_v_l ;
%                -B_u_l' , B_u_q' , beta * M_u + beta_g * A_u , zeros( Nu , Nu ) , zeros( Nu , Nu ) ; 
%                -B_f_l' , B_f_q' , zeros( Nu , Nu ) , xi * M_u + xi_g * A_u , zeros( Nu , Nu ) ;
%                -B_v_l' , B_v_q' , zeros( Nu , Nu ) , zeros( Nu , Nu ) , gamma * M_u + gamma_g * A_u ];
%         
%         F = [ zeros(Ny,1) ; zeros(Ny,1) ; -grad_x];
%         
        delta_v_y = A_z \ ( -B_u_q * v_u( Nu * ii + 1 :  Nu * ( ii + 1 ) ) - B_f_q * v_f( Nu * ii + 1 :  Nu * ( ii + 1 ) ) - B_v_q * v_v( Nu * ii + 1 :  Nu * ( ii + 1 ) ) );
        w_1 = M_o * delta_v_y - B_u_l * v_u( Nu * ii + 1 :  Nu * ( ii + 1 ) ) - B_f_l * v_f( Nu * ii + 1 :  Nu * ( ii + 1 ) ) - B_v_l * v_v( Nu * ii + 1 :  Nu * ( ii + 1 ) );
        p_w_1 = A_z' \ ( - w_1 );
        act_on_v_of_u( Nu * ii + 1 :  Nu * ( ii + 1 ) , 1 ) = -B_u_l' * delta_v_y + B_u_q' * p_w_1 + ( beta * M_u + beta_g * A_u ) * v_u( Nu * ii + 1 :  Nu * ( ii + 1 ) );
        act_on_v_of_f( Nu * ii + 1 :  Nu * ( ii + 1 ) , 1 ) = -B_f_l' * delta_v_y + B_f_q' * p_w_1 + ( xi * M_u + xi_g * A_u ) * v_f( Nu * ii + 1 :  Nu * ( ii + 1 ) );
        act_on_v_of_v( Nu * ii + 1 :  Nu * ( ii + 1 ) , 1 ) = -B_v_l' * delta_v_y + B_v_q' * p_w_1 + ( gamma * M_u + gamma_g * A_u ) * v_v( Nu * ii + 1 :  Nu * ( ii + 1 ) );      
        ii = ii + 1;
end
        act_on_v_of_u = act_on_v_of_u( Nu + 1 : end );
        act_on_v_of_f = act_on_v_of_f( Nu + 1 : end );
        act_on_v_of_v = act_on_v_of_v( Nu + 1 : end );
        d = [ act_on_v_of_u ; act_on_v_of_f ; act_on_v_of_v ];

            
    

    end

end