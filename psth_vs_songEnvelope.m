% Calcula PSTH y LFP de todos los canales de un tetrodo especificado y la
% se�al promedio (solo funciona con NNx)
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
thr_automatico = input('\Busqueda thr automatica? (1 = SI / 0 = NO) : ');

% Definimos manualmente un umbral para deteccion de spikes (en uV)
if thr_automatico == 0 
    thr = input('\nThreshold para el threshold cutting (en uV):  ');
end

% Guardo figuras?
guardar = input('\Guardo? (1 = SI / 0 = NO) : ');

% Cargamos cantidad de trials y tiempo que dura cada uno
ntrials = params.Ntrials
tiempo_file = params.tiempo_entre_estimulos

% Especifico numero de id del BOS y REV
id_BOS = params_analisis.id_bos(1)

% Tamano sliding window
t_window = 0.015; % 15 ms tamaño ventana
step = 0.001; % 1 ms step de ventana

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

% Levanta senal neuronal y la filtra para obtener: LFP de cada canal del
% tetrodo , LFP promediando todos los canales y SPIKES de cada canal
[LFP_tetrodo, LFP_canales, spikes_canales, sr_lfp]= LFP_1tetrode(directorio,...
    amplifier_channels, frequency_parameters, puerto_canal_custom, 1000);


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
          frequency_parameters, sr_lfp);

    % Calculo scores
    estimulos = score_calculator(id_BOS, estimulos, ...
        frequency_parameters, spike_times, ntrials, tiempo_file);
    
    % Calculo sliding window para cada estimulo
    for i = (1:length(estimulos))
        [sw_data, sw_times] = sliding_window(estimulos(i).spikes_norm, ...
            frequency_parameters.amplifier_sample_rate, ...
            t_window, step, tiempo_file);
        psth_sw = [sw_data, sw_times];
        estimulos(i).psth_sw = psth_sw;
    end
    clear i psth_sw

    % Guardo resultados de este canal en una struct con todos los datos
    estimulos_tetrodos(c).canal = estimulos;
end
clear c


% Calculo envolvente BOS
song = estimulos(10).song;
[t, song_env] = song_envelope(song);


song_aux = zeros(uint64(tiempo_file*estimulos(10).freq),1);
song_aux((1:length(song_env)),1) = song_env;

% Calculo MUA (Sliding Window de PSTH)
    % Ventana y step del sliding window
t_window = 0.015; % 15 ms
step = 1/estimulos(10).freq;
% step = 0.001; % 1 ms

i = 10; % id estimulo BOS
% Para cada canal recorro para calcular promedio psth
    for c = (1:1:length(estimulos_tetrodos))
        
        estimulos = estimulos_tetrodos(c).canal;

        % Calculo y ploteo sliding window para cada estimulo
        [sw_data, sw_times] = sliding_window(estimulos(i).spikes_norm, ...
            frequency_parameters.amplifier_sample_rate, ...
            t_window, step, tiempo_file);

        % Guardo PSTH de cada canal
        PSTH_avgTetrodo(1:length(sw_data),c) = sw_data;
    end     
% Promedio PSTH de todos los canales de este estimulo
PSTH_avgTetrodo_mean = mean(PSTH_avgTetrodo, 2);

PSTH_avg_aux = zeros(length(song_aux),1);
PSTH_avg_aux((1:length(PSTH_avgTetrodo_mean)), 1) = PSTH_avgTetrodo_mean;

figure()
plot(PSTH_avg_aux/max(PSTH_avg_aux));
hold on 
plot(song_aux)
xlim([0 20000])
title('checkear que tengan MISMA escala en ejeX')


% Ploteo
figure()
font_sz = 20;

% Song 
h(1) = subplot(3,1,1);
plot((1:length(song))/estimulos(10).freq, song, '-k');
legend('sound', 'Interpreter', 'none', 'FontSize',font_sz)
hold on
max_song = max(abs(song));
ylim([-max_song max_song]);
set(gca,'FontSize',font_sz)

title({'ENVOLVENTES SONIDO vs MUA para el BOS'; directorio ; ...
    strcat('tetrodo:', puerto_canal_custom)}, 'Interpreter', 'none', ...
    'FontSize',font_sz)

% Song envelope
h(2) = subplot(3,1,2);
plot(t/estimulos(10).freq, song_env, 'LineWidth', 2);
legend('envolvente automatic 2015', 'FontSize',font_sz)
hold on
ylim([0 1]);
set(gca,'FontSize',font_sz)

% MUA (PSTH sliding window)
h(3) = subplot(3,1,3);
plot(sw_times, PSTH_avgTetrodo_mean, '-r', 'LineWidth', 2');
legend('MUA promediada de 4 canales de un tetrodo', 'FontSize',font_sz)
max_mua = max(PSTH_avgTetrodo_mean);
ylim([0 max_mua]);
set(gca,'FontSize',font_sz)

linkaxes(h,'x')
xlim([0 6])
set(gca,'FontSize',font_sz)


% Guardo
if guardar == 1
    
    print_pdf(1, directorio, strcat('_',string(puerto_canal_custom),...
        '_MUA_vs_ENVOLVENTE_BOS', string(round(thr)), 'uV.pdf'))
   
end

   
    

    
    






