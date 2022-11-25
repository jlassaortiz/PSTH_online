% Script para hacer graf de psth para BOS rapido (posta)
% Necesito que el el directorio donde estan los archivos INTAN esten
% tambien los archivos de estimulos.txt y los .wav de los estimulos

close all
clear all

% Defino directorio
directorio = input('\nDirectorio: ','s');
directorio = horzcat(directorio , '/');

% Cargamos cantidad de trials y tiempo que dura cada uno
ntrials =  input('\nnúmero de trials : ');
tiempo_file = input('\ntiempo file (en seg.): ');

% Pregunto si el umbral se determina manualmente o automaticamente
thr_automatico = input('\n¿Busqueda de thr automatica? (1 = SI / 0 = NO): ');

% Definimos manualmente un umbral para deteccion de spikes (en uV)
if thr_automatico == 0 
    thr = input('\nThreshold para el threshold cutting (en uV):  ');
end

% Leer info INTAN
read_Intan_RHD2000_file(horzcat(directorio, 'info.rhd'));
clear notes spike_triggers supply_voltage_channels aux_input_channels 

% Genero songs.mat a partir de las canciones
estimulos = carga_songs(directorio);   

% Genero diccionario con nombre de los estimulos y el momento de presentacion
estimulos = find_t0s(estimulos, ntrials, tiempo_file, board_adc_channels,...
    frequency_parameters, directorio, false);

% Carga vector con parametros del analisis de datos
params_info = dir(horzcat(directorio, 'parametros_light.txt'));
params_light = readtable(horzcat(directorio,params_info.name),'Delimiter','\t');
clear params_info

% Load id BOS (index of BOS in the list of stimulus)
id_BOS = str2num(string(params_light.id_BOS));

% Load the geometry (N x M) of the probe
% The list of channels will be assign in the following order:
% (order is like matlab asigned order to subplots)
% first, by column from left to right
% then, by row from up to down
% i.e. after the first row is full from left to right, then we
% continue with the next column

geom_probe = str2num(string(params_light.geometria));

lista_canales_custom = split(params_light.canales_custom);

% Inicializo struct donde guardo la info de los bos a graficar
estimulos_arreglo = struct();

% Para cada canal de la lista de canales
for c = (1:length(lista_canales_custom))

    puerto_canal_custom = lista_canales_custom{c};

    % Traduzco custom_channel_name a native_channel_name
    traduccion = strcmp(puerto_canal_custom, ...
        {amplifier_channels(:).custom_channel_name});
    puerto_canal = amplifier_channels(traduccion).native_channel_name;
    
    % Levanto el canal de interes PASO MUY MUY LENTO!
    raw = read_INTAN_channel(directorio, puerto_canal, amplifier_channels);

    % Define el filtro
    filt_spikes = designfilt('highpassiir','DesignMethod','butter','FilterOrder',...
        4,'HalfPowerFrequency',500,...
        'SampleRate',frequency_parameters.amplifier_sample_rate);

    % Aplica filtro
    raw_filtered = filtfilt(filt_spikes, raw);
    clear filt_spikes

    % Definimos umbral de deteccion de spikes
    if thr_automatico == 1
        thr = find_thr(raw_filtered, estimulos, tiempo_file, frequency_parameters);
    end

    % Buscamos spike times (en samples, NO unidades de tiempo) por threshold cutting 
    spike_times = find_spike_times(raw_filtered, thr, frequency_parameters);

    % Genero objeto con raster de todos los estimulos
    estimulos = generate_raster_BOS(spike_times, estimulos , tiempo_file, ...
        ntrials, frequency_parameters, id_BOS);
    
    estimulos_arreglo(c).canal_nativo = puerto_canal;
    estimulos_arreglo(c).canal_custom = puerto_canal_custom;
    
    estimulos_arreglo(c).raster = estimulos(id_BOS);
    
    clear traduccion peine tetrodo canal
end


for columna = (1: geom_probe(1,2))
    
    for fila = (1:geom_probe(1,1))
        
        
        
    end 
    
end 

% Guardo figuras?
guardar = input('\n¿Guardo? (1 = SI / 0 = NO) : ');

% Cargo orden de la grilla
grilla_psth = str2num(string(params_analisis.grilla_psth(1)))

% cargo id_estimulos 
for i = (1:1:length(estimulos))
    estimulos(i).id = params_analisis.orden(i);
    estimulos(i).frec_corte = params_analisis.freq_corte(i);
    estimulos(i).tipo = categorical(params_analisis.tipo_estimulo(i));
%     estimulos(i).protocolo_id = categorical({directorio_nombre_corto});
end
clear i 


