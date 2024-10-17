% Circadian flowering model
% Alberto Gonzalez Delgado
%Centro de Biotecnologia y Genomica de Plantas (UPM/CSIC-INIA)
%04/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dRepdt= equation(t,Rep,aCO,kaCO,CO,aGI,kaGI,GI,kaTOC1,rTOC1,TOC1,kaLHY,rLHY,LHY,b,d,night,kaPRR5,rPRR5,PRR5,n,CDF3,kaCDF3,rCDF3)
    % Previous equations ----------------------------------------------
t=round(t);
n=3;
b=0;
d=0.05;

%Activation term CO
ACO=aCO*((kaCO*CO(t)).^n);

%Activation term GI
AGI=aGI*((kaGI*GI(t)).^n);

%Repression term TOC1
RTOC1=rTOC1*((kaTOC1*TOC1(t)).^n);

%Repression term LHY
RLHY=rLHY*((kaLHY*LHY(t)).^n);
%Repression term PRR5
RPRR5=rPRR5*((kaPRR5*PRR5(t)).^n);
%Repression term CDF3
RCDF3=rCDF3*((kaCDF3*CDF3(t)).^n);

% Model ------------------------------------------------------------------

dRepdt=(b+ ACO + AGI + ACO*AGI)./(0.1 + (ACO/aCO)+(AGI/aGI)+ RLHY + RTOC1 + RPRR5 + RTOC1*(ACO)+ RPRR5*(ACO)+RCDF3+(RCDF3*(AGI)))-d.*Rep;

end




