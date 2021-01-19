% Script que hace todo de una

% directorio = 'D:\Experimentos\zf-JL007_VeNa\Registros\Penetracion 3\4_170915_155938\TEST psth online\data';
directorio = '/Users/javi_lassaortiz/Documents/LSD/Finch dormido/PSTH_online/data';
directorio = horzcat(directorio , '/');

% Genero songs.mat a partir de las canciones
carga_songs;

% Leer info INTAN
read_Intan_RHD2000_file(horzcat(directorio, 'info.rhd'));
clear notes spike_triggers supply_voltage_channels aux_input_channels% no tienen info de interes

% Filtra un canal de INTAN
Filtrar_raw_data_un_canal_ht;

% Carga datos filtrados y hace un threshold cutting
plot_view_spikes_ht

% Grafica raster
loader_raster_raw_ht
