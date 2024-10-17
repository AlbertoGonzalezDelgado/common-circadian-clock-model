% Circadian flowering model
% Alberto Gonzalez Delgado
%Centro de Biotecnologia y Genomica de Plantas (UPM/CSIC-INIA)
%04/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Poplar_model(filename,clock_file)
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
    

x0=[3.0693, 0.1000, 1.8252, 0.2979, 1.0092, 1.7348, 0.5209, 3.0386, 0.5485, 0.7727, 0.0753, 0.0793];
[p,cost] = optimization(exp,data.CO,data.GI,data.TOC1,data.LHY,data.FT,x0,16,data.PRR5,data.CDF);
disp(p)
Rep_sim= model(p(1),p(2),data.CO,p(3),p(4),data.GI,p(5),p(6),data.TOC1,p(7),p(8),data.LHY,0,0.05,data.FT,93,p(9),p(10),data.PRR5,3,data.CDF,p(11),p(12));
plot(data.FT, 'r')  
hold on
plot(Rep_sim, 'r--') 
title("LD")
x_ticks = linspace(0, 133, 22);  
x_labels = 0:22; 

set(gca, 'XTick', x_ticks);  
set(gca, 'XTickLabel', x_labels); 


end
