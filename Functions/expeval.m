function [monodata, bidata, tridata] = expeval(amplitudes, lifetimes, acquisitiontime)
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
    options = optimset('MaxFunEval', 1e10, 'MaxIter', 1e10);
    [bestx1,fval1,exitflag1] = fminsearch(fun1, x01, options);
    besty1 = expdec(bestx1, data0(:,1));
    [bestx2,fval2,exitflag2] = fminsearch(fun2, x02, options);
    besty2 = expdec(bestx2, data0(:,1));
    [bestx3,fval3,exitflag3] = fminsearch(fun3, x03, options);
    besty3 = expdec(bestx3, data0(:,1));
      
    monodata = [bestx1, fval1, exitflag1];
    bidata = [bestx2, fval2, exitflag2];
    tridata = [bestx3, fval3, exitflag3];
 
end
