% Circadian flowering model
% Alberto Gonzalez Delgado
%Centro de Biotecnologia y Genomica de Plantas (UPM/CSIC-INIA)
%04/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Rep_sim= model(aCO,kaCO,CO,aGI,kaGI,GI,kaTOC1,rTOC1,TOC1,kaLHY,rLHY,LHY,b,d,Rep,night,kaPRR5,rPRR5,PRR5,n,CDF3,kaCDF3,rCDF3)
      coder.extrinsic('ode23s');
% Solve ODE ---------------------------------------------------------------

% Initial conditions
Rep0 = Rep(1);

% Time span
tspan = 1:1:133;
% Solve ODE using ode15s
[~, Rep_sim] = ode23s(@(t, Rep) equation(t, Rep, aCO,kaCO,CO,aGI,kaGI,GI,kaTOC1,rTOC1,TOC1,kaLHY,rLHY,LHY,b,d,night,kaPRR5,rPRR5,PRR5,n,CDF3,kaCDF3,rCDF3), tspan, Rep0);

end
