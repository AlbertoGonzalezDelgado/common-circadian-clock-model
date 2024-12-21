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

variables = {'CO', 'GI', 'TOC1', 'LHY', 'SP5G','PRR5','CDF'};
names = {'CO', 'GI', 'TOC1', 'LHY', 'Rep','PRR5','CDF'};
data = struct();

for i = 1:length(variables)
    temp = table2array(exp(contains(exp.Abbreviation, variables{i}), 2:145));
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
    
 p=[7.7342, 0.1162, 4.8091, 0.1030, 0.2723, 1.7014, 0.0036, 0.9614, 0.1543, 1.9161, 0.4581, 0.5816 ];
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

plot(data.Rep, 'k') 
hold on
plot(Rep_sim, 'k--') 
title("ND")
x_ticks = linspace(0, 144, 25);  
x_labels = 0:24; 

set(gca, 'XTick', x_ticks);  
set(gca, 'XTickLabel', x_labels); 

saveas(gcf, 'Simulated_SP5G.pdf')
hold on
end
