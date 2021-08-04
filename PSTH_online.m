% Script que hace todo de una

close all
clear all

% Defino directorio
directorio = input('Directorio: ','s');
directorio = horzcat(directorio , '/');

% Pregunto si el umbral se determina manualmente o automaticamente
thr_automatico = input('\n¿Busqueda de thr automatica? (1 = SI / 0 = NO) : ');

% Definimos manualmente un umbral para deteccion de spikes (en uV)
if thr_automatico == 0 
    thr = input('\nThreshold para el threshold cutting (en uV):  ');
end

% Guardo figuras?
guardar = input('\n¿Guardo? (1 = SI / 0 = NO) : ');

% Carga vector con parametros del analisis de datos
params_info = dir(horzcat(directorio, '*parametros_protocolo*.txt'));
params = readtable(horzcat(directorio,params_info.name),'Delimiter','\t');
clear params_info

% Carga vector con parametros del analisis de datos
params_info = dir(horzcat(directorio, '*parametros_analisis*.txt'));
params_analisis = readtable(horzcat(directorio,params_info.name),'Delimiter','\t');
clear params_info

% Cargo valores de puerto-canal
puerto = char(params.Puerto);
canal = params.Canal;
puerto_canal = horzcat(puerto, '-0', num2str(canal,'%.2d'))
clear puerto canal

% Cargamos cantidad de trials y tiempo que dura cada uno
ntrials = params.Ntrials
tiempo_file = params.tiempo_entre_estimulos

% Especifico numero de id del BOS y REV
id_BOS = params_analisis.id_bos(1)
id_REV = params_analisis.id_rev(1)

% Cargo orden de la grilla
grilla_psth = str2num(string(params_analisis.grilla_psth(1)))

% Genero songs.mat a partir de las canciones
estimulos = carga_songs(directorio);    

% cargo id_estimulos 
for i = (1:1:length(estimulos))
    estimulos(i).id = params_analisis.orden(i);
    estimulos(i).frec_corte = params_analisis.freq_corte(i);
    estimulos(i).tipo = categorical(params_analisis.tipo_estimulo(i));
%     estimulos(i).protocolo_id = categorical({directorio_nombre_corto});
end
clear i 

% Leer info INTAN
read_Intan_RHD2000_file(horzcat(directorio, 'info.rhd'));
clear notes spike_triggers supply_voltage_channels aux_input_channels 

% Levanto el canal de interes
raw = read_INTAN_channel(directorio, puerto_canal, amplifier_channels);

% Define el filtro
filt_spikes = designfilt('highpassiir','DesignMethod','butter','FilterOrder',...
    4,'HalfPowerFrequency',500,'SampleRate',frequency_parameters.amplifier_sample_rate);

% Aplica filtro
raw_filtered = filtfilt(filt_spikes, raw);
clear filt_spikes

% Genero diccionario con nombre de los estimulos y el momento de presentacion
estimulos = find_t0s(estimulos, ntrials, tiempo_file, board_adc_channels, frequency_parameters, directorio, false);

% Definimos umbral de deteccion de spikes
if thr_automatico == 1
    thr = find_thr(raw_filtered, estimulos, tiempo_file, frequency_parameters);
end

% Buscamos spike times (en samples, NO unidades de tiempo) por threshold cutting 
spike_times = find_spike_times(raw_filtered, thr, frequency_parameters);

% Genero objeto con raster de todos los estimulos
estimulos = generate_raster(spike_times, estimulos , tiempo_file, ntrials, frequency_parameters);

% Calculo scores
estimulos = score_calculator(id_BOS, estimulos, frequency_parameters, spike_times, ntrials);

% Ploteo spike shapes
plot_spikes_shapes(raw_filtered, spike_times, thr, frequency_parameters, directorio)

% Ploteo Grilla PSTH
plot_some_raster(grilla_psth, id_BOS, estimulos, estimulos, frequency_parameters, tiempo_file, ntrials, puerto_canal, thr, directorio, spike_times);

estimulos = struct2table(estimulos);



% Selecciono datos de ese protocolo 
pasa_altos = estimulos(estimulos.tipo == 'up' , :);
pasa_bajos = estimulos(estimulos.tipo == 'down' , :);

% Plotear INT PASA-ALTOS
figure();
plot(pasa_altos.frec_corte, pasa_altos.int_norm, '-o')
title({strcat('INT_PASA-ALTOS_', datestr(now, 'yyyy-mm-dd')); ...
string(directorio) ; ...
strcat(string(puerto_canal), "  " , string(thr), "uV", "  ntrials:", string(ntrials), "  t_inter_estimulo:", string(tiempo_file)) }, 'Interpreter','None')
legend

% Plotear INT PASA-BAJOS
figure();
plot(pasa_bajos.frec_corte, pasa_bajos.int_norm, '-o')
title({strcat('INT_PASA-BAJOS_', datestr(now, 'yyyy-mm-dd')); ...
string(directorio) ; ...
strcat(string(puerto_canal), "  " , string(thr), "uV", "  ntrials:", string(ntrials), "  t_inter_estimulo:", string(tiempo_file)) }, 'Interpreter','None')
legend


% Plotear CORR PASA-ALTOS
figure();
plot(pasa_altos.frec_corte, pasa_altos.corr, '-o')
title({strcat('CORR_PASA-ALTOS_', datestr(now, 'yyyy-mm-dd')); ...
string(directorio) ; ...
strcat(string(puerto_canal), "  " , string(thr), "uV", "  ntrials:", string(ntrials), "  t_inter_estimulo:", string(tiempo_file)) }, 'Interpreter','None')
legend

% Plotear CORR PASA-BAJOS
figure();
plot(pasa_bajos.frec_corte, pasa_bajos.corr, '-o')
title({strcat('CORR_PASA-BAJOS_', datestr(now, 'yyyy-mm-dd')); ...
string(directorio) ; ...
strcat(string(puerto_canal), "  " , string(thr), "uV", "  ntrials:", string(ntrials), "  t_inter_estimulo:", string(tiempo_file)) }, 'Interpreter','None')
legend


% Guardo
if guardar == 1
    print_png(1, directorio, strcat('_spike-shape_', string(round(thr)), 'uV'))
    print_pdf(2, directorio, strcat('_grilla_', string(round(thr)), 'uV.pdf'))
    print_pdf(3, directorio, strcat('_INT_pasa-ALTOS', '.pdf'))
    print_pdf(4, directorio, strcat('_INT_pasa-BAJOS', '.pdf'))
    print_pdf(5, directorio, strcat('_CORR_pasa-ALTOS', '.pdf'))
    print_pdf(6, directorio, strcat('_CORR_pasa-BAJOS', '.pdf'))
end

clear estimulos_aux j i spike_times  