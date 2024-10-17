% Circadian flowering model
% Alberto Gonzalez Delgado
%Centro de Biotecnologia y Genomica de Plantas (UPM/CSIC-INIA)
%04/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function print_result(filename,clock_file)
% import data -------------------------------------
coder.extrinsic('detectImportOptions');
coder.extrinsic('readtable');
opts = detectImportOptions(filename, 'Delimiter', '\t', 'FileType', 'text');
data = readtable(filename, opts);
opts = detectImportOptions(clock_file, 'Delimiter', ',', 'FileType', 'text');
cic_data = readtable(clock_file, opts);
%merge data
exp = innerjoin(data, cic_data(:, {'ID', 'Abbreviation'}), 'Keys', 'ID');

variables = {'CO', 'GI', 'TOC1', 'LHY', 'FT','PRR5','CDF'};
names = {'CO', 'GI', 'TOC1', 'LHY', 'FT','PRR5','CDF'};
data = struct();
for i = 1:length(variables)
    temp = table2array(exp(contains(exp.Abbreviation, variables{i}), 2:134));
    temp_normalized = (temp - min(temp)) / (max(temp) - min(temp));
    data.(names{i}) = temp_normalized;
end

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
    
p=[ 2.8736  ,  0.1021   , 1.7656    ,0.3792  , 1.5689  ,  1.5590  ,  1.0859,    3.2218  ,  0.4593  ,  1.8292    ,0.0146  ,  0.7674];


%Plot model --------------------------------------------------------------
%LD

disp("p")
disp(p)
FT_sim= model(p(1),p(2),data.CO,p(3),p(4),data.GI,p(5),p(6),data.TOC1,p(7),p(8),data.LHY,0,0.05,data.FT,93,p(9),p(10),data.PRR5,3,data.CDF,p(11),p(12));
writematrix(data.FT, 'FT_LD.tsv', 'Delimiter', '\t', 'FileType', 'text');
writematrix(FT_sim, 'FT_simulated_LD.tsv', 'Delimiter', '\t', 'FileType', 'text');
subplot(2,1,1)
plot(data.FT, 'r')  
hold on
plot(FT_sim, 'r--') 
title("LD")
x_ticks = linspace(0, 133, 25);  
x_labels = 0:24; 

set(gca, 'XTick', x_ticks);  
set(gca, 'XTickLabel', x_labels); 

%ND
disp(exp)
data = struct();

for i = 1:length(variables)
    temp = table2array(exp(contains(exp.Abbreviation, variables{i}), 135:267));
    temp_normalized = (temp - min(temp)) / (max(temp) - min(temp));
    data.(names{i}) = temp_normalized;
end
p=[ 2.8736  ,  0.1021   , 1.7656    ,0.3792  , 1.5689  ,  1.5590  ,  1.0859,    3.2218  ,  0.4593  ,  1.8292    ,4.0146  ,  4.7674];

FT_sim= model(p(1),p(2),data.CO,p(3),p(4),data.GI,p(5),p(6),data.TOC1,p(7),p(8),data.LHY,0,0.05,data.FT,93,p(9),p(10),data.PRR5,3,data.CDF,p(11),p(12));
subplot(2,1,2)
writematrix(data.FT, 'FT_ND.tsv', 'Delimiter', '\t', 'FileType', 'text');
writematrix(FT_sim, 'FT_simulated_ND.tsv', 'Delimiter', '\t', 'FileType', 'text');
plot(data.FT, 'b') 
hold on
plot(FT_sim, 'b--') 
title("ND")
x_ticks = linspace(0, 133, 25);  
x_labels = 0:24; 

set(gca, 'XTick', x_ticks);  
set(gca, 'XTickLabel', x_labels); 
saveas(gcf, 'Simulated_FT.pdf')
hold on
end
