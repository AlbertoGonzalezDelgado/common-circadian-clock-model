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

variables = {'CO', 'GI', 'TOC1', 'LHY', 'FT','PRR5','CDF1'};
names = {'CO', 'GI', 'TOC1', 'LHY', 'FT','PRR5','CDF'};
data = struct();
for i = 1:length(variables)
    temp = table2array(exp(contains(exp.Abbreviation, variables{i}), 111:219));
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
    
p=[3.9693, 0.1000, 1.3252, 0.1979, 0.5092, 0.1348, 0.5209, 3.0386, 0.5485, 0.4727, 0.0753, 0.0793];

%Plot model --------------------------------------------------------------
Rep_sim= model(p(1),p(2),data.CO,p(3),p(4),data.GI,p(5),p(6),data.TOC1,p(7),p(8),data.LHY,0,0.05,data.FT,72,p(9),p(10),data.PRR5,3,data.CDF,p(11),p(12));

str = sprintf('%.4f, ', p); 
str = str(1:end-2); 
disp(['p=', str])
subplot(3,1,1)
plot(data.FT, 'r')
hold on
plot(Rep_sim, 'r--') 
title("LD")
x_ticks = linspace(0, 111, 21);  
x_labels = 2:22; 

set(gca, 'XTick', x_ticks);  
set(gca, 'XTickLabel', x_labels); 

writematrix(data.FT, 'FT_LD.tsv', 'Delimiter', '\t', 'FileType', 'text');
writematrix(Rep_sim, 'FT_simulated_LD.tsv', 'Delimiter', '\t', 'FileType', 'text');

%-------------------------------------------
%SD
p=[3.9693, 0.1000, 1.3252, 0.1979, 0.5092, 0.1348, 0.5209, 3.0386, 0.5485, 0.4727, 4.0753, 4.0793];

variables = {'CO', 'GI', 'TOC1', 'LHY', 'FT','PRR5','CDF1'};
names = {'CO', 'GI', 'TOC1', 'LHY', 'FT','PRR5','CDF'};
data = struct();
for i = 1:length(variables)
    temp = table2array(exp(contains(exp.Abbreviation, variables{i}), 220:328));
    temp_normalized = (temp - min(temp)) / (max(temp) - min(temp));

    data.(names{i}) = temp_normalized;
end

%Plot model --------------------------------------------------------------
Rep_sim= model(p(1),p(2),data.CO,p(3),p(4),data.GI,p(5),p(6),data.TOC1,p(7),p(8),data.LHY,0,0.05,data.FT,36,p(9),p(10),data.PRR5,3,data.CDF,p(11),p(12));

str = sprintf('%.4f, ', p); 
str = str(1:end-2); 
disp(['p=', str])
subplot(3,1,2)
plot(data.FT, 'r')
hold on
plot(Rep_sim, 'r--') 
title("SD")
writematrix(data.FT, 'FT_SD.tsv', 'Delimiter', '\t', 'FileType', 'text');
writematrix(Rep_sim, 'FT_simulated_SD.tsv', 'Delimiter', '\t', 'FileType', 'text');
x_ticks = linspace(0, 111, 21);  
x_labels = 2:22; 

set(gca, 'XTick', x_ticks);  
set(gca, 'XTickLabel', x_labels); 


saveas(gcf, 'Simulated_FT.pdf')
hold on
end
