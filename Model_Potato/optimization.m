% Circadian flowering model
% Alberto Gonzalez Delgado
%Centro de Biotecnologia y Genomica de Plantas (UPM/CSIC-INIA)
%04/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p,cost]= optimization(exp,CO,GI,TOC1,LHY,Rep,x0,night,PRR5,CDF3)
      coder.extrinsic('simulannealbnd');
      coder.extrinsic('optimoptions');
        % Define parammeters --------------------------------------------------
        lb = repmat(0.1, 1, 6); %min limit
        ub = repmat(10, 1, 6); %max limit
        
        %1: aCO
        %2: kaCO
        %3: aGI
        %4: kaGI
        %5: ka TOC1
        %6: rTOC1
        %7: kaLHY
        %8: rLHY
        %9: kaPRR5
        %10: rPRR5
        %11: kaCDF3
        %12: rCDF3
    
        % Cost function
         fun = @(x) cost_function(exp, x(1), x(2), CO, x(3), x(4), GI, x(5), x(6), TOC1, x(7), x(8), LHY, 0,0.05 ,Rep, ...
             night,x(9),x(10),PRR5,3,CDF3,x(11),x(12));
        
        % Optimization options
        opts = optimoptions(@simulannealbnd, 'MaxIterations', 1e24,'MaxFunctionEvaluations', 1e24,'TolFun',1e-24);
    
        % Run optimization ----------------------------------------------------
        [p,cost] = simulannealbnd(fun,x0, lb, ub, opts);
    end