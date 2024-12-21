% Circadian flowering model
% Alberto Gonzalez Delgado
%Centro de Biotecnologia y Genomica de Plantas (UPM/CSIC-INIA)
%04/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function potato_model(filename,clock_file)
rng(123);
% import data -------------------------------------
coder.extrinsic('detectImportOptions');
coder.extrinsic('readtable');
opts = detectImportOptions(filename, 'Delimiter', '\t', 'FileType', 'text');
data = readtable(filename, opts);
opts = detectImportOptions(clock_file, 'Delimiter', ',', 'FileType', 'text');
cic_data = readtable(clock_file, opts);
%merge data
exp = innerjoin(data, cic_data(:, {'ID', 'Abbreviation'}), 'Keys', 'ID');

variables = {'CO', 'GI', 'TOC1', 'LHY', 'SP5G','PRR5','CDF'};
names = {'CO', 'GI', 'TOC1', 'LHY', 'Rep','PRR5','CDF'};
data = struct();

for i = 1:length(variables)
    temp = table2array(exp(contains(exp.Abbreviation, variables{i}), 2:134));
    temp_normalized = (temp - min(temp)) / (max(temp) - min(temp));

    data.(names{i}) = temp_normalized;
end
disp(data)
% Initial values
x0=[6,0.2,... %CO
    1.6,0.185,... %GI
   0.1,1.5,... %TOC1
    0.5,3,... % LHY
    0.2,1.5,... %PRR5
    0.5,0.5]
[p,cost] = optimization(exp,data.CO,data.GI,data.TOC1,data.LHY,data.Rep,x0,72,data.PRR5,data.CDF);
 %Save results
%Plot model --------------------------------------------------------------
Rep_sim= model(p(1),p(2),data.CO,p(3),p(4),data.GI,p(5),p(6),data.TOC1,p(7),p(8),data.LHY,0,0.05,data.Rep,72,p(9),p(10),data.PRR5,3,data.CDF,p(11),p(12));

str = sprintf('%.4f, ', p); 
str = str(1:end-2); 
disp(['p=', str])
subplot(3,1,1)
plot(data.Rep, 'r')
hold on
plot(Rep_sim, 'r--') 
title("LD")

variables = {'CO', 'GI', 'TOC1', 'LHY', 'SP5G','PRR5','CDF'};
names = {'CO', 'GI', 'TOC1', 'LHY', 'Rep','PRR5','CDF'};
data = struct();

for i = 1:length(variables)
    temp = table2array(exp(contains(exp.Abbreviation, variables{i}), 135:267));
    temp_normalized = (temp - min(temp)) / (max(temp) - min(temp));

    data.(names{i}) = temp_normalized;
end
plot(data.Rep, 'k') 
hold on
plot(Rep_sim, 'k--') 
title("ND")



end
