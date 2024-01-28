function [xsol,fval,history,grad_norm,lambda] = runfmincon_ufv_Hessian( cost_with_grad_param , constraint ,  param , z_0 , lb )

% % history.x = [];
history.fval = [];
% % searchdir = [];
grad_norm = [];

Nu = param.Nu;
mu_0=param.mu_0;
tol_constraint=param.tol_constraint;
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
    'MaxIterations', 2000 , 'FiniteDifferenceType' , 'central' , 'SubproblemAlgorithm' , 'cg' , 'HessianMultiplyFcn', @HessMultFcn );
%options=optimset('Algorithm','sqp','ScaleProblem','obj-and-constr','DerivativeCheck','off','GradObj','on','TolX',1e-16,'TolFun',1e-16,'MaxFunEvals',inf,'MaxIter',6000,'Display','iter');

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


        
        switch param.constraint_type
            case{'u_f_v'}
        [y,p] = get_state_and_adjoint(x,param);
        B_dd_u = param.B_dd_u;
        B_dd_f = param.B_dd_f;
        B_dl=param.B_dl;
        Nu = param.Nu;
        
        u_opt = x( 1 : Nu , 1 );
        f_opt = x( Nu + 1 : 2 * Nu , 1 );
        v_opt = x( 2 * Nu + 1 : end , 1 );
        
        B_u = spmatrix( ttv( B_dd_u , u_opt , 3 ) );
        B_f = spmatrix( ttv( B_dd_f , f_opt , 3 ) );
        B_v = spmatrix( ttv( B_dl , v_opt , 3 ) );
        B_u_q = spmatrix( ttv( B_dd_u , y , 1 ) );
        B_f_q = spmatrix( ttv( B_dd_f , y , 1 ) );
        B_v_q = spmatrix( ttv( B_dl , y , 1 ) );
        B_u_l = spmatrix( ttv( B_dd_u , p , 1 ) );
        B_f_l = spmatrix( ttv( B_dd_f , p , 1 ) );
        B_v_l = spmatrix( ttv( B_dl , p , 1 ) );
        
        A = param.A;
        A_z = A + B_u + B_f + B_v;
        Ny = length( y );
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
        delta_v_y = A_z \ ( -B_u_q * v( 1 : Nu ) - B_f_q * v( Nu + 1 : 2 * Nu ) - B_v_q * v( 2 * Nu + 1 : end ) );
        w_1 = M_o * delta_v_y - B_u_l * v( 1 : Nu ) - B_f_l * v( Nu + 1 : 2 * Nu ) - B_v_l * v( 2 * Nu + 1 : end );
        p_w_1 = A_z' \ ( - w_1 );
        act_on_v_of_u = -B_u_l' * delta_v_y + B_u_q' * p_w_1 + ( beta * M_u + beta_g * A_u ) * v( 1 : Nu , 1 );
        act_on_v_of_f = -B_f_l' * delta_v_y + B_f_q' * p_w_1 + ( xi * M_u + xi_g * A_u ) * v( Nu + 1 : 2 * Nu , 1 );
        act_on_v_of_v = -B_v_l' * delta_v_y + B_v_q' * p_w_1 + ( gamma * M_u + gamma_g * A_u ) * v( 2 * Nu + 1 : end , 1 );      
        d = [ act_on_v_of_u ; act_on_v_of_f ; act_on_v_of_v ];

            case{'uv'}
        [~,y,p] = compute_gradient( x ,param );
        B_dd = param.B_dd;
        B_dl=param.B_dl;
        Nu = param.Nu;
        
        u_opt = x( 1 : Nu , 1 );
        v_opt = x( Nu + 1 : 2 * Nu , 1 );
       
        
        B_u = spmatrix( ttv( B_dd , u_opt , 3 ) );
        B_v = spmatrix( ttv( B_dl , v_opt , 3 ) );
        B_u_q = spmatrix( ttv( B_dd , y , 1 ) );
        B_v_q = spmatrix( ttv( B_dl , y , 1 ) );
        B_u_l = spmatrix( ttv( B_dd, p , 1 ) );
        B_v_l = spmatrix( ttv( B_dl , p , 1 ) );
        
        A = param.A;
        A_z = A + B_u + B_v;
        Ny = length( y );
        M_o = param.M_o;
        M_u = param.M_u;
        A_u = param.A_u;
        beta = param.beta;
        beta_g   = param.beta_g;
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
        delta_v_y = A_z \ ( -B_u_q * v( 1 : Nu ) - B_v_q * v( Nu + 1 : 2 * Nu ) );
        w_1 = M_o * delta_v_y - B_u_l * v( 1 : Nu ) - B_v_l * v( Nu + 1 : 2 * Nu );
        p_w_1 = A_z' \ ( - w_1 );
        act_on_v_of_u = -B_u_l' * delta_v_y + B_u_q' * p_w_1 + ( beta * M_u + beta_g * A_u ) * v( 1 : Nu , 1 );
       
        act_on_v_of_v = -B_v_l' * delta_v_y + B_v_q' * p_w_1 + ( gamma * M_u + gamma_g * A_u ) * v( Nu + 1 : 2 * Nu );      
        d = [ act_on_v_of_u ; act_on_v_of_v ];


        end
    end



end