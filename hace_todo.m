% Script que hace todo de una
clear all
close all

% Genero songs.mat a partir de las canciones
carga_songs

% Filtra un canal de INTAN
Filtrar_raw_data_un_canal_ht

% Carga datos filtrados y hace un threshold cutting
plot_view_spikes_ht

% Pregunta el tipo de protocolo
protocolo = input('Tipo de protocolo (1): solo BOS, (2): BOS, CON, REV:  ');

if protocolo==1;
    
    loader_raster_BOS_raw_ht

else
    
    loader_raster_raw_ht
    
end

set(figure(1),'Visible','on')
movegui(figure(1),'east')
movegui(figure(2),'center')