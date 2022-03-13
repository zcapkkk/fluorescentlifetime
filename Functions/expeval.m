function [monodata, bidata, tridata, err1, err2, err3, ft12, ft23] = expeval(amplitudes, lifetimes, acquisitiontime)
    [data0] = groupproject_simulate(amplitudes, lifetimes, acquisitiontime);
    meas = data0(:,2);

    %drop zero term
    data0=data0(data0(:,2)~=0,:);
    
    % define lambda function as evaluation/optimisation function 
    fun1 = @(x)expdec1eval(x, data0);
    fun2 = @(x)expdec2eval(x, data0);
    fun3 = @(x)expdec3eval(x, data0);
    
    % predict some preliminary guess values
    x01 = [meas(1), 1+rand(1,1)*3, meas(end) ];
    x02 = [0.5*meas(1)*ones(1,2), 1+rand(1,2)*3, meas(end)];
    x03 = [0.33*meas(1)*ones(1,3), 1+rand(1,3)*3, meas(end)];
    
    % perform optimisation
    %options = optimset('MaxFunEval', 1e10, 'MaxIter', 1e10);
    A=[]; b=[];
    [bestx1,fval1,exitflag1,output,lambda,grad,hessian] = fmincon(fun1, x01, A, b);
    besty1 = expdec(bestx1, data0(:,1));
    err1 = sqrt(diag(inv(hessian)));
    [bestx2,fval2,exitflag2,output,lambda,grad,hessian] = fmincon(fun2, x02, A, b);
    besty2 = expdec(bestx2, data0(:,1));
    err2 = sqrt(diag(inv(hessian)));
    [bestx3,fval3,exitflag3,output,lambda,grad,hessian] = fmincon(fun3, x03, A, b);
    besty3 = expdec(bestx3, data0(:,1));
    err3 = sqrt(diag(inv(hessian)));
    c = length(data0);
    ft12 = fcdf(fval1/fval2,c-3,c-5);
    ft23 = fcdf(fval2/fval3,c-5,c-7);

    monodata = [bestx1,fval1, exitflag1];
    bidata = [bestx2, fval2, exitflag2];
    tridata = [bestx3, fval3, exitflag3];
 
end
