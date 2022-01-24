% Este programa compara los LFP de cada canal de un tetrodo con el LFP
% promedio del tetrodo (TODA LA SEÑAL)

close all
clear all

% Defino directorio
directorio = input('Directorio: ','s');
directorio = horzcat(directorio , '/');

% Defino tetrodo (cuatro canales) a levantar con nombre custom name
peine = input('\nDefino canal a levantar (custom name)\n \nPeine (X): ');
tetrodo = input('\nTetrodo (X): ');
puerto_canal_custom = horzcat('P',num2str(peine),'-','T', ... 
    num2str(tetrodo));

clear peine tetrodo

% Guardo figuras?
guardar = input('\n¿Guardo? (1 = SI / 0 = NO) : '); 

% Leer info INTAN
read_Intan_RHD2000_file(horzcat(directorio, 'info.rhd'));
clear notes spike_triggers supply_voltage_channels aux_input_channels 
    
% Calcula LFP de cada canal del tetrodo especificado y el promedio de los 4
[LFP_tetrodo, LFP_canales, spikes_canales] = LFP_1tetrode(directorio, ...
    amplifier_channels, frequency_parameters, puerto_canal_custom);
 
% Calculo correlacion del LFP de cada canal con el LFP promedio

corr = zeros(4,1); % Matriz fila donde guardo correlaciones

% Para cada canal calculo la correlación con el LFP promedio
for i = (1:1:4)
    corr_aux = corrcoef(LFP_tetrodo(:,1), LFP_canales(:,i));
    corr(i) = corr_aux(2);
end


% PLOTEA LFP promedio y de cada canal individual

figure('units','normalized','outerposition',[0 0 1 1])

h(1) = subplot(5,1,1);
plot(LFP_tetrodo)
title('promedio tetrodo','FontSize',22)

h(2) = subplot(5,1,2);
plot(LFP_canales(:,1))
legend( strcat('corr: ',num2str(corr(1))) ,'FontSize',22 )
title('canal 1','FontSize',22)

h(3) = subplot(5,1,3);
plot(LFP_canales(:,2))
legend( strcat('corr: ',num2str(corr(2))) ,'FontSize',22)
title('canal 2','FontSize',22)

h(4) = subplot(5,1,4);
plot(LFP_canales(:,3))
legend( strcat('corr: ',num2str(corr(3))) ,'FontSize',22)
title('canal 3','FontSize',22)

h(5) = subplot(5,1,5);
plot(LFP_canales(:,4))
legend( strcat('corr: ',num2str(corr(4))) ,'FontSize',22)
title('canal 4','FontSize',22)

sgtitle({strcat('LFP_', datestr(now, 'yyyy-mm-dd')); ...
string(directorio) ; ...
strcat(string(puerto_canal_custom)) }, ...
'Interpreter','None')

ylim(h, [-2000 2000]);

linkaxes(h, 'x');


if guardar == 1
    
     print_png(1, directorio, strcat('_', puerto_canal_custom, ...
         '_comparacion_LFP'))
    
end

