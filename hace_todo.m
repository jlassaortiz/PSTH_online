% Script que hace todo de una
close all
clear all

directorio = input('Directorio: ','s');
directorio = horzcat(directorio , '/');

% Genero songs.mat a partir de las canciones
estimulos = carga_songs(directorio);

% Leer info INTAN
read_Intan_RHD2000_file(horzcat(directorio, 'info.rhd'));
clear notes spike_triggers supply_voltage_channels aux_input_channels 

% Eligo un canal de un puerto en especifico
puerto = input('Puerto a levantar: ','s');
canal = input('Canal de INTAN a filtrar (X-0XX):  ');
puerto_canal = [puerto '-0' num2str(canal,'%.2d')];

% Levanto el canal de interes
raw = read_INTAN_channel(directorio, puerto_canal, amplifier_channels);

% Define el filtro
filt_spikes = designfilt('highpassiir','DesignMethod','butter','FilterOrder',...
    4,'HalfPowerFrequency',500,'SampleRate',frequency_parameters.amplifier_sample_rate);

% Aplica filtro
raw_filtered = filtfilt(filt_spikes, raw);
clear puerto canal filt_spikes

% Carga datos filtrados y hace un threshold cutting
plot_view_spikes_ht

% Grafica raster
loader_raster_raw_ht
