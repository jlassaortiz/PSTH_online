function [LFP_tetrodo, LFP_canales, spikes_canales]= LFP_1tetrode(...
    directorio,amplifier_channels,frequency_parameters,puerto_canal_custom)

% Calcula LFP de un tetrodo especificado (promedia 4 canales)
%
%   ! NECESITA QUE ANTES SE CORRA read_Intan_RHD2000_file
%
%   INPUTS:
%   directorio: (str) directorio donde esta el archivo INTAN
%
%   amplifier_channels: salida read_Intan_RHD2000_file
%
%   frequency_parameters: salida read_Intan_RHD2000_file
%
%   puerto_canal_custom: (str) nombre custom del tetrodo a levantar, el
%   nombre del canal debe ser de la forma puerto_canal_custom + _x
%   donde x es el nombre del canal y varia de 1 a 4
%
%   OUTPUT:
%   LFP_tetrodo: (matris Nx1) vector fila con señal filtrada y promediada 
%   de los 4 canales del tetrodo
%   
%   LFP_canales: (matris Nx4) matris donde c/columna es el LFP de 1 canal
%   
%   spikes_canales: (matris Nx4) matris donde c/columna es la señal
%   filtrada para conservar spikes de 1 canal


% Define el filtro para LFP y SPIKES
filt_LFP = designfilt('lowpassiir','DesignMethod','butter',...
    'HalfPowerFrequency',35,'FilterOrder', 4, ...
    'SampleRate', frequency_parameters.amplifier_sample_rate);

% Define el filtro para LFP y SPIKES
filt_LFP_h = designfilt('highpassiir','DesignMethod','butter',...
    'HalfPowerFrequency',25,'FilterOrder', 4, ...
    'SampleRate', frequency_parameters.amplifier_sample_rate);


filt_spikes = designfilt('highpassiir','DesignMethod','butter',...
    'HalfPowerFrequency',500,'FilterOrder', 4, ...
    'SampleRate',frequency_parameters.amplifier_sample_rate);
    
% Creo objeto donde guardo trazas de LFP y SPIKES
LFP_canales = [];
spikes_canales = []; 

% Para cada canal del tetrodo
for i = (1:1:4)

    % Genero nombre de canal a levantar
    puerto_canal_custom_aux = horzcat(puerto_canal_custom, ...
        '-',num2str(i));

    % Traduzco custom_channel_name a native_channel_name
    traduccion = strcmp(puerto_canal_custom_aux, ...
        {amplifier_channels(:).custom_channel_name});
    puerto_canal = amplifier_channels(traduccion).native_channel_name;

    % Levanto el canal de interes (TARDA MUCHO)
    raw = read_INTAN_channel(directorio, puerto_canal, ...
        amplifier_channels);

    % Aplica el filtro para LFP y SPIKES
    LFP = filtfilt(filt_LFP, raw);
    LFP = filtfilt(filt_LFP_h, raw);
    
    spikes = filtfilt(filt_spikes, raw);

    % Adjunto traza LFP y SPIKES de este canal a la "lista" de LFP y SPIKES
    LFP_canales = [LFP_canales,LFP];
    spikes_canales = [spikes_canales, spikes];
    
end

% Promedio los 4 canales de los tetrodos para obtener un LFP promediado
LFP_tetrodo = mean(LFP_canales, 2);

end

