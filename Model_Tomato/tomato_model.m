% Circadian flowering model
% Alberto Gonzalez Delgado
%Centro de Biotecnologia y Genomica de Plantas (UPM/CSIC-INIA)
%04/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tomato_model(filename,clock_file)
rng(123);
% import data -------------------------------------
coder.extrinsic('detectImportOptions');
coder.extrinsic('readtable');
opts = detectImportOptions(filename, 'Delimiter', '\t', 'FileType', 'text');
data = readtable(filename, opts);
opts = detectImportOptions(clock_file, 'Delimiter', ';', 'FileType', 'text');
cic_data = readtable(clock_file, opts);
%merge data
exp = innerjoin(data, cic_data(:, {'ID', 'Abbreviation'}), 'Keys', 'ID');

%LD
variables = {'CO', 'GI', 'TOC1', 'LHY', 'SP5G','PRR5','CDF3'};
names = {'CO', 'GI', 'TOC1', 'LHY', 'Rep','PRR5','CDF3'};
data = struct();

for i = 1:length(variables)
    temp = table2array(exp(contains(exp.Abbreviation, variables{i}), 2:134));
    temp_normalized = (temp - min(temp)) / (max(temp) - min(temp));

    data.(names{i}) = temp_normalized;
end

x_values = 0.1:0.05:3;
cost_values=1:20:100;

%Preallocating
P = zeros(num_steps, 13); 
Costs = zeros(num_steps, 1);

x0=[8,0.1,5,0.1,0.1,1,0.1,1,0.1,1,0.5,0.5];
[p,cost] = optimization(exp,data.CO,data.GI,data.TOC1,data.LHY,data.Rep,x0,89,data.PRR5,data.CDF3);
disp(cost)


% Find miniumum cost

Rep_sim= model(p(1),p(2),data.CO,p(3),p(4),data.GI,p(5),p(6),data.TOC1,p(7),p(8),data.LHY,0,0.05,data.Rep,89,p(9),p(10),data.PRR5,3,data.CDF3,p(11),p(12));

%Plot model --------------------------------------------------------------
%LD
plot(data.Rep, 'r')
hold on
plot(Rep_sim, 'r--') 
title("LD")

end
