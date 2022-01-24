% Calcula PSTH y LFP de todos los canales de un tetrodo especificado y la
% señal promedio (solo funciona con NNx)
% Hace varios graficos

close all
clear all


% Cargo y defino parametros %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Defino directorio
directorio = input('Directorio: ','s');
directorio = horzcat(directorio , '/');

% Leer info INTAN
read_Intan_RHD2000_file(horzcat(directorio, 'info.rhd'));
clear notes spike_triggers supply_voltage_channels aux_input_channels 

% Carga vector con parametros del protocolo experimental
params = readtable(horzcat(directorio,'/parametros_protocolo.txt'), ...
    'Delimiter','\t');

% Carga vector con parametros para el analisis de datos
params_analisis = readtable(...
    horzcat(directorio,'/parametros_analisis.txt'), 'Delimiter','\t');

% Defino tetrodo a levantar con nombre custom name
peine = input('\nDefino tetrodo a levantar (custom name) \n\nPeine (X): ');
tetrodo = input('\nTetrodo (X): ');
puerto_canal_custom = horzcat('P',num2str(peine),'-','T',num2str(tetrodo));

% Pregunto si ploteo toda la grilla o solo algunos estimulos
plot_grilla = input('\nPloteo toda la grilla? (1 = SI / 0 = NO): ');
if plot_grilla == 0
    grilla_psth = input('\nMatris lineal con numero ID estimulos : ');
end

% Pregunto si el umbral se determina manualmente o automaticamente
thr_automatico = input('\n¿Busqueda thr automatica? (1 = SI / 0 = NO) : ');

% Definimos manualmente un umbral para deteccion de spikes (en uV)
if thr_automatico == 0 
    thr = input('\nThreshold para el threshold cutting (en uV):  ');
end

% Guardo figuras?
guardar = input('\n¿Guardo? (1 = SI / 0 = NO) : ');

% Cargamos cantidad de trials y tiempo que dura cada uno
ntrials = params.Ntrials
tiempo_file = params.tiempo_entre_estimulos

% Especifico numero de id del BOS y REV
id_BOS = params_analisis.id_bos(1)

% Cargo orden de la grilla
if plot_grilla == 1
    grilla_psth = str2num(string(params_analisis.grilla_psth(1)))
end

% Genero songs.mat a partir de las canciones
estimulos = carga_songs(directorio);    

% Cargo id_estimulos 
for i = (1:1:length(estimulos))
    estimulos(i).id = params_analisis.orden(i);
    estimulos(i).frec_corte = params_analisis.freq_corte(i);
    estimulos(i).tipo = categorical(params_analisis.tipo_estimulo(i));
end
clear i 


% Levanto señal neuronal y analizo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Levanta señal neuronal y la filtra para obtener: LFP de cada canal del
% tetrodo , LFP promediando todos los canales y SPIKES de cada canal
[LFP_tetrodo, LFP_canales, spikes_canales]= LFP_1tetrode(directorio,...
    amplifier_channels, frequency_parameters, puerto_canal_custom);

% Genero struct con nombre de los estimulos y el momento de presentacion
estimulos = find_t0s(estimulos, ntrials, tiempo_file, ...
    board_adc_channels, frequency_parameters, directorio, false);

% Genero struct donde guardo datos de todos los canales
estimulos_tetrodos = struct();

% Para cada canal del tetrodo
for c = ( 1:1: size(spikes_canales,2) )
    
    raw_filtered = spikes_canales(:,c);
    LFP = LFP_canales(:,c);

    % Calculamos el umbral de deteccion de spikes
    if thr_automatico == 1
        thr = find_thr(raw_filtered, estimulos, tiempo_file, ...
        frequency_parameters);
    end

    % Buscamos spike times (en samples, NO unidades tiempo) por umbral
    spike_times =find_spike_times(raw_filtered, thr, frequency_parameters);

    % Genero objeto con RASTERS de todos los estimulos
    estimulos = generate_raster(spike_times, estimulos , tiempo_file, ... 
        ntrials, frequency_parameters);

    % Calculo LFP promediado por estimulo todos los trials
    estimulos = trialAverage_LFP(LFP, estimulos, tiempo_file, ntrials, ...
        frequency_parameters);

    % Calculo scores
    estimulos = score_calculator(id_BOS, estimulos, ...
        frequency_parameters, spike_times, ntrials);

    % Guardo resultados de este canal en una struct con todos los datos
    estimulos_tetrodos(c).canal = estimulos;
end


% PLOTEO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Ploteo Grilla PSTH
plot_some_raster_LFP_1tetrode(grilla_psth, id_BOS, estimulos_tetrodos, ...
    frequency_parameters, tiempo_file, ntrials,thr,directorio,spike_times);
sgtitle({datestr(now, 'yyyy-mm-dd'); ...
string(directorio) ; ...
strcat('tetrodo = ',string(puerto_canal_custom),"  |  ", ...
string(thr), "uV", "  |  ", "ntrials:", string(ntrials), "  |  ", ...
"t_inter_estimulo:", string(tiempo_file)) }, 'Interpreter','None',...
'FontSize',20)

%%%%%%%%%%%%%%%%% desarrolle hasta aca.
% Falta ver si puedo calcular corr de PSTH_sw y LFP de estimulos vs BOS

% % Selecciono datos de ese protocolo 
% estimulos_table = struct2table(estimulos);
% pasa_altos = estimulos_table(estimulos_table.tipo == 'up' , :);
% pasa_bajos = estimulos_table(estimulos_table.tipo == 'down' , :);
% 
% 
% % Plotear INT PASA-ALTOS
% figure();
% plot(pasa_altos.frec_corte, pasa_altos.int_norm, '-o')
% title({strcat('INT_PASA-ALTOS_', datestr(now, 'yyyy-mm-dd')); ...
% string(directorio) ; ...
% strcat(string(puerto_canal), " = ",string(puerto_canal_custom),"  |  ",...
% string(thr), "uV", "  |  ", "ntrials:", string(ntrials), "  |  ", ...
% "t_inter_estimulo:", string(tiempo_file)) }, 'Interpreter','None')
% legend
% set(gca,'FontSize',20)
% 
% 
% % Plotear INT PASA-BAJOS
% figure();
% plot(pasa_bajos.frec_corte, pasa_bajos.int_norm, '-o')
% title({strcat('INT_PASA-BAJOS_', datestr(now, 'yyyy-mm-dd')); ...
% string(directorio) ; ...
% strcat(string(puerto_canal), " = ",string(puerto_canal_custom),"  |  ", ...
% string(thr), "uV", "  |  ", "ntrials:", string(ntrials), "  |  ", ...
% "t_inter_estimulo:", string(tiempo_file)) }, 'Interpreter','None')
% legend
% set(gca,'FontSize',20)
% 
% 
% % Plotear CORR PASA-ALTOS
% figure();
% plot(pasa_altos.frec_corte, pasa_altos.corr, '-o')
% title({strcat('CORR_PASA-ALTOS_', datestr(now, 'yyyy-mm-dd')); ...
% string(directorio) ; ...
% strcat(string(puerto_canal), " = ",string(puerto_canal_custom),"  |  ", ...
% string(thr), "uV", "  |  ", "ntrials:", string(ntrials), "  |  ", ...
% "t_inter_estimulo:", string(tiempo_file)) }, 'Interpreter','None')
% legend
% set(gca,'FontSize',20)
% 
% 
% % Plotear CORR PASA-BAJOS
% figure();
% plot(pasa_bajos.frec_corte, pasa_bajos.corr, '-o')
% title({strcat('CORR_PASA-BAJOS_', datestr(now, 'yyyy-mm-dd')); ...
% string(directorio) ; ...
% strcat(string(puerto_canal), " = ",string(puerto_canal_custom),"  |  ",...
% string(thr), "uV", "  |  ", "ntrials:", string(ntrials), "  |  ", ...
% "t_inter_estimulo:", string(tiempo_file)) }, 'Interpreter','None')
% legend
% set(gca,'FontSize',20)


% Guardo
if guardar == 1
%     print_png(1, directorio, strcat('_',string(puerto_canal_custom),...
%         '_spike-shape_', string(round(thr)), 'uV'))

    print_pdf(1, directorio, strcat('_',string(puerto_canal_custom),...
        '_grilla_PSTH-LFP-tetrode', string(round(thr)), 'uV.pdf'))
    
%     print_pdf(2, directorio, strcat('_',string(puerto_canal_custom),...
%         '_INT_pasa-ALTOS', '.pdf'))
%     print_pdf(3, directorio, strcat('_',string(puerto_canal_custom),...
%         '_INT_pasa-BAJOS', '.pdf'))
%     print_pdf(4, directorio, strcat('_',string(puerto_canal_custom),...
%         '_CORR_pasa-ALTOS', '.pdf'))
%     print_pdf(5, directorio, strcat('_',string(puerto_canal_custom),...
%         '_CORR_pasa-BAJOS', '.pdf'))
end

clear estimulos_aux j i  