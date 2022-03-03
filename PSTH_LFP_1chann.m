% Calcula PSTH y LFP de cada estimulo para UN SOLO CANAL (TUNGSTENO o NNx)
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

% Carga vector con parametros del analisis de datos
%params_info = dir(horzcat(directorio, '*parametros_protocolo*.txt'));
params = readtable(horzcat(directorio,'/parametros_protocolo.txt'), ...
    'Delimiter','\t');
clear params_info

% Carga vector con parametros del analisis de datos
%params_info = dir(horzcat(directorio, '*parametros_analisis*.txt'));
params_analisis = readtable( ...
    horzcat(directorio,'/parametros_analisis.txt'), ...
    'Delimiter','\t');
clear params_info

% Pregunto si determino canal de manera manual o automatica (esta
% explicitado en parametros que se cargan)
canal_automatico = ...
    input('\n�Determino canal automaticamente? (1 = SI / 0 = NO) : ');

if canal_automatico == 1
    % Cargo valores de puerto-canal del archivo parametros dentro del dir
    puerto = char(params.Puerto);
    canal = params.Canal;
    puerto_canal = horzcat(puerto, '-0', num2str(canal,'%.2d'))
    puerto_canal_custom = puerto_canal;
    clear puerto canal 

else
    % Defino canal a levantar con nombre custom name
    peine = ...
        input('\nDefino canal a levantar (custom name) \n \nPeine (X): ');
    tetrodo = input('\nTetrodo (X): ');
    canal = input('\nCanal (X): ');

    puerto_canal_custom = horzcat('P',num2str(peine),'-', ...
        'T',num2str(tetrodo),'-',num2str(canal));

    % Traduzco custom_channel_name a native_channel_name
    traduccion = strcmp(puerto_canal_custom, ...
        {amplifier_channels(:).custom_channel_name});
    
    puerto_canal = amplifier_channels(traduccion).native_channel_name;
    clear traduccion peine tetrodo canal
end

% Pregunto si ploteo toda la grilla o solo algunos estimulos
plot_grilla = input('\nPloteo toda la grilla? (1 = SI / 0 = NO): ');
if plot_grilla == 0
    grilla_psth = input('\nMatris lineal con numero ID estimulos : ');
end

% Pregunto si el umbral se determina manualmente o automaticamente
thr_automatico = ...
    input('\n�Busqueda de thr automatica? (1 = SI / 0 = NO) : ');

% Definimos manualmente un umbral para deteccion de spikes (en uV)
if thr_automatico == 0 
    thr = input('\nThreshold para el threshold cutting (en uV):  ');
end

% Guardo figuras?
guardar = input('\n�Guardo? (1 = SI / 0 = NO) : ');

% Cargamos cantidad de trials y tiempo que dura cada uno
ntrials = params.Ntrials
tiempo_file = params.tiempo_entre_estimulos

% Especifico numero de id del BOS y REV
id_BOS = params_analisis.id_bos(1)
id_REV = params_analisis.id_rev(1)

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

% Levanto se�al neuronal y analizo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Levanto el canal de interes
raw = read_INTAN_channel(directorio, puerto_canal, amplifier_channels);

% Define el filtro para SPIKES y LFP
filt_spikes = designfilt('highpassiir','DesignMethod','butter',...
    'FilterOrder', 4,'HalfPowerFrequency',500,...
    'SampleRate',frequency_parameters.amplifier_sample_rate);

filt_LFP = designfilt('lowpassiir','DesignMethod','butter',...
    'FilterOrder', 4,'HalfPowerFrequency', 300, ...
    'SampleRate', frequency_parameters.amplifier_sample_rate);

% Aplica filtros
raw_filtered = filtfilt(filt_spikes, raw);
LFP = filtfilt(filt_LFP, raw);

clear filt_spikes

% Genero struct con nombre de los estimulos y el momento de presentacion
estimulos = find_t0s(estimulos, ntrials, tiempo_file, ...
    board_adc_channels, frequency_parameters, directorio, false);

% Definimos umbral de deteccion de spikes
if thr_automatico == 1
    thr =find_thr(raw_filtered,estimulos,tiempo_file,frequency_parameters);
end
clear thr_automatico

% Buscamos spike times (en samples, NO unidades de tiempo) por threshold 
spike_times = find_spike_times(raw_filtered, thr, frequency_parameters);

% Genero objeto con raster de todos los estimulos
estimulos = generate_raster(spike_times, estimulos , tiempo_file, ...
    ntrials, frequency_parameters);

% Calculo LFP promediado por estimulo todos los trials
estimulos = trialAverage_LFP(LFP, estimulos, tiempo_file, ntrials, ...
    frequency_parameters);

% Calculo scores
estimulos = score_calculator(id_BOS, estimulos, frequency_parameters, ...
    spike_times, ntrials);

% Calculo sliding window para cada estimulo
t_window = 0.015; % 15 ms
step = 0.001; % 1 ms
for i = (1:length(estimulos))
    [sw_data, sw_times] = sliding_window(estimulos(i).spikes_norm, ...
        frequency_parameters.amplifier_sample_rate, ...
        t_window, step);
    psth_sw = [sw_data, sw_times];
    estimulos(i).psth_sw = psth_sw;
end
clear i psth_sw


% Guardo txt?
guardar_txt = input('\n�Guardo PSTH_sw y LFP_promedio BOS? (1 = SI / 0 = NO) : ');

if guardar_txt == 1
    
    writematrix(estimulos(id_BOS).psth_sw, ...
        [directorio '/PSTHsw_1chann_BOS_' puerto_canal_custom '.txt'],...
        'Delimiter','tab')
    
    writematrix(estimulos(id_BOS).LFP_promedio,...
        [directorio '/LFP_1chann_BOS_' puerto_canal_custom '.txt'],...
        'Delimiter','tab')
end


% PLOTEO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Ploteo spike shapes
plot_spikes_shapes(raw_filtered, spike_times, thr, frequency_parameters,...
    directorio)
sgtitle({datestr(now, 'yyyy-mm-dd'); ...
string(directorio) ; ...
strcat(string(puerto_canal), " = ",string(puerto_canal_custom),"  |  ", ...
string(thr), "uV", "  |  ", "ntrials:", string(ntrials), "  |  ",...
"t_inter_estimulo:", string(tiempo_file)) }, 'Interpreter','None',...
'FontSize',8)


% Ploteo Grilla PSTH
plot_some_raster_LFP_1channel(grilla_psth, id_BOS, estimulos, estimulos,...
    frequency_parameters, tiempo_file, ntrials, puerto_canal, thr, ...
    directorio, spike_times);
sgtitle({datestr(now, 'yyyy-mm-dd'); ...
string(directorio) ; ...
strcat(string(puerto_canal), " = ",string(puerto_canal_custom),"  |  ", ...
string(thr), "uV", "  |  ", "ntrials:", string(ntrials), "  |  ", ...
"t_inter_estimulo:", string(tiempo_file)) }, 'Interpreter','None',...
'FontSize',20)


% Selecciono datos de ese protocolo 
estimulos_table = struct2table(estimulos);
pasa_altos = estimulos_table(estimulos_table.tipo == 'up' , :);
pasa_bajos = estimulos_table(estimulos_table.tipo == 'down' , :);


% Plotear INT PASA-ALTOS
figure();
plot(pasa_altos.frec_corte, pasa_altos.int_norm, '-o')
title({strcat('INT_PASA-ALTOS_', datestr(now, 'yyyy-mm-dd')); ...
string(directorio) ; ...
strcat(string(puerto_canal), " = ",string(puerto_canal_custom),"  |  ",...
string(thr), "uV", "  |  ", "ntrials:", string(ntrials), "  |  ", ...
"t_inter_estimulo:", string(tiempo_file)) }, 'Interpreter','None')
legend
set(gca,'FontSize',20)


% Plotear INT PASA-BAJOS
figure();
plot(pasa_bajos.frec_corte, pasa_bajos.int_norm, '-o')
title({strcat('INT_PASA-BAJOS_', datestr(now, 'yyyy-mm-dd')); ...
string(directorio) ; ...
strcat(string(puerto_canal), " = ",string(puerto_canal_custom),"  |  ", ...
string(thr), "uV", "  |  ", "ntrials:", string(ntrials), "  |  ", ...
"t_inter_estimulo:", string(tiempo_file)) }, 'Interpreter','None')
legend
set(gca,'FontSize',20)


% Plotear CORR PASA-ALTOS
figure();
plot(pasa_altos.frec_corte, pasa_altos.corr, '-o')
title({strcat('CORR_PASA-ALTOS_', datestr(now, 'yyyy-mm-dd')); ...
string(directorio) ; ...
strcat(string(puerto_canal), " = ",string(puerto_canal_custom),"  |  ", ...
string(thr), "uV", "  |  ", "ntrials:", string(ntrials), "  |  ", ...
"t_inter_estimulo:", string(tiempo_file)) }, 'Interpreter','None')
legend
set(gca,'FontSize',20)


% Plotear CORR PASA-BAJOS
figure();
plot(pasa_bajos.frec_corte, pasa_bajos.corr, '-o')
title({strcat('CORR_PASA-BAJOS_', datestr(now, 'yyyy-mm-dd')); ...
string(directorio) ; ...
strcat(string(puerto_canal), " = ",string(puerto_canal_custom),"  |  ",...
string(thr), "uV", "  |  ", "ntrials:", string(ntrials), "  |  ", ...
"t_inter_estimulo:", string(tiempo_file)) }, 'Interpreter','None')
legend
set(gca,'FontSize',20)


% Guardo
if guardar == 1
    print_png(1, directorio, strcat('_',string(puerto_canal_custom),...
        '_spike-shape_', string(round(thr)), 'uV'))
    print_pdf(2, directorio, strcat('_',string(puerto_canal_custom),...
        '_grilla_PSTH-LFP', string(round(thr)), 'uV.pdf'))
    print_pdf(3, directorio, strcat('_',string(puerto_canal_custom),...
        '_INT_pasa-ALTOS', '.pdf'))
    print_pdf(4, directorio, strcat('_',string(puerto_canal_custom),...
        '_INT_pasa-BAJOS', '.pdf'))
    print_pdf(5, directorio, strcat('_',string(puerto_canal_custom),...
        '_CORR_pasa-ALTOS', '.pdf'))
    print_pdf(6, directorio, strcat('_',string(puerto_canal_custom),...
        '_CORR_pasa-BAJOS', '.pdf'))
end

clear estimulos_aux j i  