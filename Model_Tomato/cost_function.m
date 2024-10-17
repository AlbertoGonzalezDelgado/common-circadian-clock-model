% Circadian flowering model
% Alberto Gonzalez Delgado
%Centro de Biotecnologia y Genomica de Plantas (UPM/CSIC-INIA)
%04/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Cost= cost_function(exp,aCO,kaCO,CO,aGI,kaGI,GI,kaTOC1,rTOC1,TOC1,kaLHY,rLHY,LHY,b,d,Rep,night,kaPRR5,rPRR5,PRR5,n,CDF3,kaCDF3,rCDF3)
coder.extrinsic('str2double')
coder.extrinsic('trapz')
coder.extrinsic('sum')
coder.extrinsic('max')
coder.extrinsic('min')

%Obtain simulation --------------------------------------------------------
Rep_sim=model(aCO,kaCO,CO,aGI,kaGI,GI,kaTOC1,rTOC1,TOC1,kaLHY,rLHY,LHY,b,d,Rep,night,kaPRR5,rPRR5,PRR5,n,CDF3,kaCDF3,rCDF3);
%Calculate cost -----------------------------------------------------------
%Observed data
Rep=Rep(:); 
Rep_sim=Rep_sim(:);
[Amp_obs,amp_index]=max(Rep);

expStruct = table2struct(exp(:,2:134));
fieldNames = fieldnames(expStruct);
ZT= fieldNames{amp_index};ZT = strrep(ZT, '_', '.');
Phase_obs = str2double(ZT(find(ZT=='.', 1, 'first')+1:end));

AUC_obs = trapz(Rep);
DEI_obs=sum(Rep);
FC_obs = max(Rep) - min(Rep);
%Predicted data
[Amp_pred,amp_index]=max(Rep_sim);
ZT= fieldNames{amp_index};ZT = strrep(ZT, '_', '.');
Phase_pred = str2double(ZT(find(ZT=='.', 1, 'first')+1:end));
AUC_pred = trapz(Rep_sim);
DEI_pred=sum(Rep_sim);
FC_pred = max(Rep_sim) - min(Rep_sim);
Esum=0;
if length(Rep) ~= length(Rep_sim)
    Esum=10000000;
else
    for i = 1:length(Rep)
        Ediff=(Rep(i)-Rep_sim(i))^2;
        Esum=Esum+Ediff;
    end
end

%Cost
Cost=150*Esum+(Amp_obs/Amp_pred).^2+0.05*(AUC_obs-AUC_pred).^2+(DEI_obs-DEI_pred).^2+150*(FC_obs-FC_pred).^2+0.05*(Phase_obs-Phase_pred).^2;

end
