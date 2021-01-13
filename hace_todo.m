% Script que hace todo de una

directorio = 'D:\Experimentos\zf-JL007_VeNa\Registros\Penetracion 3\4_170915_155938\TEST psth online\data';
directorio = horzcat(directorio , '\');

% Genero songs.mat a partir de las canciones
carga_songs;

% % figure();
% % plot(estimulo(1).song);
% % title(estimulo(1).name, 'Interpreter', 'none');
% % 
% % figure();
% % plot(estimulo(2).song);
% % title(estimulo(2).name, 'Interpreter', 'none');
% % 
% % figure();
% % plot(estimulo(3).song);
% % title(estimulo(3).name, 'Interpreter', 'none');

% Leer info INTAN
read_Intan_RHD2000_file(horzcat(directorio, 'info.rhd'));
num_channels = length(amplifier_channels);
clear notes spike_triggers supply_voltage_channels % no tienen info de interes

% Filtra un canal de INTAN
Filtrar_raw_data_un_canal_ht;

% Carga datos filtrados y hace un threshold cutting
plot_view_spikes_ht

% Grafica raster
loader_raster_raw_ht

set(figure(1),'Visible','on')
movegui(figure(1),'east')
movegui(figure(2),'center')