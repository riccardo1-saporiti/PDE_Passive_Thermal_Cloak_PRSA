function [xsol,fval,history,grad_norm,lambda] = runfmincon_optimization( cost_with_grad_param , constraint , z_0 , lb )

% % history.x = [];
history.fval = [];
% % searchdir = [];
grad_norm = [];

Acon = [];
b=[];
Aeq = [];
beq = []; 
ub= [];
NL_con = constraint; 
objfun = cost_with_grad_param;
       options = optimoptions( 'fmincon','SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,'CheckGradients',false,'Display','iter', ...
           'OutputFcn',@outfun,'OptimalityTolerance', 0, ...
    'StepTolerance', 0, ...
    'MaxFunctionEvaluations', inf,...
    'MaxIterations', 2000 , 'FiniteDifferenceType' , 'central' );


[xsol,fval,~,~,lambda] = fmincon(objfun,z_0,Acon,b,Aeq,beq,lb,ub,NL_con,options);

    function stop = outfun(~,optimValues,state)
     stop = false;
 
     switch state
         case 'init'
             hold on
         case 'iter'
           history.fval = [history.fval; optimValues.fval];
            grad_norm = [grad_norm;... 
                        optimValues.firstorderopt];
         case 'done'
             hold off
         otherwise
     end
end
end