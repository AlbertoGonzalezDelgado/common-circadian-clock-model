% Circadian flowering model
% Alberto Gonzalez Delgado
%Centro de Biotecnologia y Genomica de Plantas (UPM/CSIC-INIA)
%04/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function print_result(filename,clock_file)
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

% Obtain gene exp ----------------------------------------------------
%LD
variables = {'CO', 'GI', 'TOC1', 'LHY', 'SP5G','PRR5','CDF3'};
names = {'CO', 'GI', 'TOC1', 'LHY', 'Rep','PRR5','CDF3'};
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
    

%Plot model --------------------------------------------------------------
%LD
variables = {'CO', 'GI', 'TOC1', 'LHY', 'SP5G','PRR5','CDF3'};
names = {'CO', 'GI', 'TOC1', 'LHY', 'Rep','PRR5','CDF3'};
data = struct();

for i = 1:length(variables)
    temp = table2array(exp(contains(exp.Abbreviation, variables{i}), 2:134));
    temp_normalized = (temp - min(temp)) / (max(temp) - min(temp));

    data.(names{i}) = temp_normalized;
end
disp("p")
p=[7.9682, 0.1088, 4.9738, 0.1115, 0.1229, 0.9929, 0.1533, 0.9470, 0.1289, 1.0599, 0.4584, 0.5154];
disp(p)
Rep_sim= model(p(1),p(2),data.CO,p(3),p(4),data.GI,p(5),p(6),data.TOC1,p(7),p(8),data.LHY,0,0.05,data.Rep,93,p(9),p(10),data.PRR5,3,data.CDF3,p(11),p(12));
subplot(2,1,1)
plot(data.Rep, 'r')  
hold on
plot(Rep_sim, 'r--') 
title("LD")
x_ticks = linspace(0, 133, 25);  
x_labels = 0:24; 
writematrix(data.Rep, 'SP5G_LD.tsv', 'Delimiter', '\t', 'FileType', 'text');
writematrix(Rep_sim, 'SP5G_simulated_LD.tsv', 'Delimiter', '\t', 'FileType', 'text');
set(gca, 'XTick', x_ticks);  
set(gca, 'XTickLabel', x_labels); 

%SD
variables = {'CO', 'GI', 'TOC1', 'LHY', 'SP5G','PRR5','CDF3'};
data = struct();
p=[7.9682, 0.1088, 4.9738, 0.1115, 0.1229, 0.9929, 0.1533, 0.9470, 0.1289, 1.0599, 4.4584, 4.5154];

for i = 1:length(variables)
    temp = table2array(exp(contains(exp.Abbreviation, variables{i}), 268:400));
    temp_normalized = (temp - min(temp)) / (max(temp) - min(temp));
    data.(names{i}) = temp_normalized;
end

Rep_sim= model(p(1),p(2),data.CO,p(3),p(4),data.GI,p(5),p(6),data.TOC1,p(7),p(8),data.LHY,0,0.05,data.Rep,46,p(9),p(10),data.PRR5,3,data.CDF3,p(11),p(12));
subplot(2,1,2)
plot(data.Rep, 'b') 
hold on
plot(Rep_sim, 'b--') 
title("SD")
x_ticks = linspace(0, 133, 25);  
x_labels = 0:24; 

set(gca, 'XTick', x_ticks);  
set(gca, 'XTickLabel', x_labels); 
saveas(gcf, 'Simulated_SP5G.pdf')
writematrix(data.Rep, 'SP5G_SD.tsv', 'Delimiter', '\t', 'FileType', 'text');
writematrix(Rep_sim, 'SP5G_simulated_SD.tsv', 'Delimiter', '\t', 'FileType', 'text');
hold on
end
