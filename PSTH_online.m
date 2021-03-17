% Script que hace todo de una

close all

% Defino directorio
directorio = input('Directorio: ','s');
directorio = horzcat(directorio , '/');

% Pregunto si el umbral se determina manualmente o automaticamente
thr_automatico = input('¿Busqueda de thr automatica? (1 = SI / 0 = NO) : ');

% Definimos manualmente un umbral para deteccion de spikes (en uV)
if thr_automatico == 0 
    thr = input('Threshold para el threshold cutting (en uV):  ');
end

% Carga vector con parametros del analisis de datos
params_info = dir(horzcat(directorio, '*parametros*.txt'));
params = readtable(horzcat(directorio,params_info.name),'Delimiter','\t','ReadVariableNames',false);
clear params_info

% Cargo valores de puerto-canal
puerto = char(params.Var2(1));
canal = char(params.Var2(2));
puerto_canal = horzcat(puerto, '-0', num2str(canal,'%.2d'))
clear puerto canal

% Cargamos cantidad de trials y tiempo que dura cada uno
ntrials = str2num(char(params.Var2(3)))
tiempo_file = str2num(char(params.Var2(4)))

% Especifico numero de id del BOS
id_BOS = str2num(char(params.Var2(5)))

% Cargo orden de la grilla
% grilla = str2num(string(params.Var2(6)))

% Genero songs.mat a partir de las canciones
estimulos = carga_songs(directorio);

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
clear puerto canal filt_spikes raw

% Genero diccionario con nombre de los estimulos y el momento de presentacion
t0s_dictionary = find_t0s(estimulos, ntrials, tiempo_file, board_adc_channels, frequency_parameters, directorio, false);

    
% Definimos un umbral para threshold cutting de manera automatica (en uV)
if thr_automatico == 1
    thr = find_thr(raw_filtered, t0s_dictionary, tiempo_file, frequency_parameters);
end
clear thr_automatico

% Buscamos spike por threshold cutting
spike_times = find_spike_times(raw_filtered, thr, frequency_parameters);

% Genero objeto con raster de todos los estimulos
rasters = generate_raster(spike_times, t0s_dictionary, tiempo_file, ntrials, frequency_parameters);

% Carga datos filtrados y hace un threshold cutting
% plot_spikes_shapes(raw_filtered, spike_times, thr, frequency_parameters, directorio)

% Grafica raster tidi
% plot_all_raster(estimulos, id_BOS, rasters, frequency_parameters, tiempo_file, ntrials, puerto_canal, thr, directorio)
