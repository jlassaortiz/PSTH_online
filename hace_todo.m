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

% Filtra un canal de INTAN
Filtrar_raw_data_un_canal_ht;

% Carga datos filtrados y hace un threshold cutting
plot_view_spikes_ht

% Grafica raster
loader_raster_raw_ht
