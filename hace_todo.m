% Script que hace todo de una

directorio = 'E:\Javi - finch durmiendo\zf-JL016_NaVe\Datos\12-01-21\3 penetracion\zf-JL016-NaVe_p3_id_1_210112_170231';
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
