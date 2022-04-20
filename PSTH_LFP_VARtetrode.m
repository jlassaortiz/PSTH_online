% Calcula PSTH y LFP de todos los canales de un tetrodo especificado y la
% senal promedio (solo funciona con NNx)
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

tetrodos_list = struct();

full_peine = input('\nGrafico TODO el peine?(1 = SI / 0 = NO): ');

if full_peine == 1
    count = 1;
    for tetrodo = [4, 3, 2, 1]
        for peine = [1,2,3,4]
            puerto_canal_custom = horzcat('P',num2str(peine),'-','T',num2str(tetrodo));
            tetrodos_list(count).puerto_canal_custom = puerto_canal_custom;
            count = count + 1;
        end 
    end
    
else
    masTetrodos = 1;
    count = 1;
    while masTetrodos == 1
        % Defino tetrodos a levantar con nombre custom name
        peine = input('\nDefino tetrodo a levantar (custom name) \n\nPeine (X): ');
        tetrodo = input('\nTetrodo (X): ');
        puerto_canal_custom = horzcat('P',num2str(peine),'-','T',num2str(tetrodo));

        tetrodos_list(count).puerto_canal_custom = puerto_canal_custom;

        masTetrodos = input('\nAgrego mas tetrodos? (1 = SI / 0 = NO): ');

        count = count + 1;
    end
end

% Indico estimulo a graficar (solo uno)
estimulo_ID = input('\nNumero ID del estimulo a graficar: ');

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

% Tamano del vector sw una vez que recorre todo el canto
size_sw = 0;
tf = t_window;
while tf <= tiempo_file
   size_sw = size_sw + 1;
   tf = tf + step;
end

% Vector de tiempos del sw (en segundos!)
times_sw = zeros(size_sw,1);
times_sw(1,1) = t_window /2;
for i = (2:1:size_sw)
    times_sw(i,1) = times_sw(i -1,1) + step;
end
clear i

% Genero songs.mat a partir de las canciones
estimulos = carga_songs(directorio);

% Cargo id_estimulos
for i = (1:1:length(estimulos))
    estimulos(i).id = params_analisis.orden(i);
    estimulos(i).frec_corte = params_analisis.freq_corte(i);
    estimulos(i).tipo = categorical(params_analisis.tipo_estimulo(i));
end
clear i


% Levanto senal neuronal y analizo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Genero struct con nombre de los estimulos y el momento de presentacion
estimulos = find_t0s(estimulos, ntrials, tiempo_file, ...
    board_adc_channels, frequency_parameters, directorio, false);

% Conservo solo el estimulo de interes
estimulos = estimulos(estimulo_ID);

% Inicializo struct donde guardo los struct de cada tetrodo a analizar
estimulos_VARchann = estimulos;
estimulos_VARchann.VARchann = struct();

for t = (1:1:length(tetrodos_list))

    estimulos_1chan = estimulos;

    % Indico el nombre del tetrodo a analizar
    puerto_canal_custom = tetrodos_list(t).puerto_canal_custom;
    estimulos_VARchann.VARchann(t).puerto_canal_custom = puerto_canal_custom;

    % Levanta senal neuronal y la filtra para obtener: LFP de cada canal del
    % tetrodo , LFP promediando todos los canales y SPIKES de cada canal
    % PASO MUY LENTO, MEJORAR !
    [LFP_tetrodo, LFP_canales, spikes_canales]= LFP_1tetrode(directorio,...
        amplifier_channels, frequency_parameters, puerto_canal_custom);

    % Genero struct donde guardo datos de todos los canales
    estimulos_tetrodos = struct();

    % Inicializo variables donde guardar psth y lfp de cada canal
    PSTH_avgTetrodo = zeros(size_sw, 4);
    LFP_avgTetrodo = [];

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
        estimulos_1chan = generate_raster(spike_times, estimulos_1chan ,...
            tiempo_file, ntrials, frequency_parameters);

        % Calculo LFP promediado por estimulo todos los trials de cada canal
        estimulos_1chan = trialAverage_LFP(LFP, estimulos_1chan, ...
            tiempo_file, ntrials, frequency_parameters);

       % Guardo LFP de cada canal
       LFP_avgTetrodo(:,c) = estimulos_1chan.LFP_promedio;

%         % Calculo scores
%         estimulos = score_calculator(1, estimulos, ...
%             frequency_parameters, spike_times, ntrials);

        % Calculo sliding window de cada canal
        [sw_data, sw_times] = sliding_window(estimulos_1chan.spikes_norm, ...
            frequency_parameters.amplifier_sample_rate, ...
            t_window, step);
        psth_sw = [sw_data, sw_times];

        % Guardo PSTH de cada canal
        PSTH_avgTetrodo(1:length(sw_data),c) = sw_data;
        estimulos_1chan.PSTH_1chann = [sw_data, sw_times];

        clear i psth_sw

        % Guardo resultados de este canal en una struct con todos los datos
        estimulos_tetrodos(c).canal = estimulos_1chan;
    end
    clear c

    % Promedio LFP de todos los canales del tetrodo de este estimulo
    LFP_avgTetrodo = mean(LFP_avgTetrodo, 2);
    LFP_avgTetrodo(:,2) = (0:1:length(LFP_avgTetrodo)-1) / ...
        frequency_parameters.amplifier_sample_rate;
    % Agrego vector de tiempos (segundos!)

    % Promedio PSTH de todos los canales de este estimulo
    PSTH_avgTetrodo = mean(PSTH_avgTetrodo, 2);
    PSTH_avgTetrodo(:,2) = times_sw;

    % Guardo curvas PSTH y LFP promedio tetrodos
    estimulos_VARchann.VARchann(t).LFP_tetrodo = LFP_avgTetrodo;
    estimulos_VARchann.VARchann(t).PSTH_tetrodo = PSTH_avgTetrodo;
    estimulos_VARchann.VARchann(t).canales_tet = estimulos_tetrodos;
end

clear estimulos estimulos_1chan estimulos_tetrodos

% Inicializo struct para plotear mas facil
plotear = struct();

% Song la defino una vez porque es la misma siempre en este caso
song = estimulos_VARchann.song;
song_times = (1:1:length(song)) / estimulos_VARchann.freq;
song = horzcat(song, song_times');

% Completo struct con datos
for t = (1:1:length(tetrodos_list))

    % Song
    plotear(t).song = song;

    % Subtitulo
    plotear(t).subTitle = estimulos_VARchann.VARchann(t).puerto_canal_custom;

    % PSTH
    plotear(t).psth = estimulos_VARchann.VARchann(t).PSTH_tetrodo;

    % LFP
    plotear(t).lfp = estimulos_VARchann.VARchann(t).LFP_tetrodo;
end

clear song song_times


% PLOTEO
plotSimple_song_psth_lfp(plotear, false, true)

suptitle2({datestr(now, 'yyyy-mm-dd'); ...
string(directorio) ; ...
estimulos_VARchann.name;
strcat("ntrials:",string(ntrials), ...
"  |  t_inter_estimulo:",string(tiempo_file))})


% Guardo
if guardar == 1
    print_pdf(1, directorio,strcat('_PSTH-LFP-VARtetrode_', ...
        estimulos_VARchann.name,'.pdf'))
end

clear estimulos_aux j i
