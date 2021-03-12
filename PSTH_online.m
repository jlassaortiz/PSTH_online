% Script que hace todo de una

close all

% Defino directorio
directorio = input('Directorio: ','s');
directorio = horzcat(directorio , '/');

% Eligo un canal de un puerto en especifico
puerto = input('Puerto a levantar: ','s');
canal = input('Canal de INTAN a filtrar (X-0XX):  ');
puerto_canal = [puerto '-0' num2str(canal,'%.2d')];

% Cargamos cantidad de trials y tiempo que dura cada uno
ntrials = input('Numero de trials: ');
tiempo_file = input('Tiempo entre estimulos (en s): ');

% Definimos un umbral para threshold cutting (en uV)
thr = input('Threshold para el threshold cutting (en uV):  ');

% Especifico numero de id del BOS
id_BOS = input('id BOS: ');

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
clear puerto canal filt_spikes

% Buscamos spike por threshold cutting
spike_times = find_spike_times(raw_filtered, thr, frequency_parameters);

% Carga datos filtrados y hace un threshold cutting
% plot_spikes_shapes(raw_filtered, spike_times, thr, frequency_parameters, directorio)

% Genero diccionario con nombre de los estimulos y el momento de presentacion
t0s_dictionary = find_t0s(estimulos, ntrials, tiempo_file, board_adc_channels, frequency_parameters, directorio);

% Genero objeto con raster de todos los estimulos
rasters = generate_raster(spike_times, t0s_dictionary, tiempo_file, ntrials, frequency_parameters);

% Grafica raster
% plot_all_raster(estimulos, id_BOS, rasters, frequency_parameters, tiempo_file, ntrials, puerto_canal, thr, directorio)

% Grafica grilla de rasters ordenada por valores de C y L-traquea
plot_some_raster([1, 2, 3, 7, 10, 4, 8, 11, 5, 9, 12, 6], id_BOS,  estimulos, rasters, frequency_parameters, tiempo_file, ntrials, puerto_canal, thr, directorio)
xlim([0, 5000]);

% Guarda la grilla
print_pdf(1, directorio, strcat('_grilla_umbral_', string(thr)))
