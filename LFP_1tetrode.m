function [LFP_tetrodo, LFP_canales] = LFP_1tetrode(directorio,amplifier_channels, ...
    frequency_parameters, puerto_canal_custom)

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
%   LFP: (matris) vector fila con señal filtrada y promediada de los 4
%   canales del tetrodo



% Define el filtro para LFP
filt_LFP = designfilt('lowpassiir','DesignMethod','butter',...
    'HalfPowerFrequency',300,'FilterOrder', 4, ...
    'SampleRate', frequency_parameters.amplifier_sample_rate);
    
% Creo objeto donde guardo trazas de LFP de manera auxilar
LFP_canales = [];

% Para cada canal del tetrodo
for i = (1:1:4)

    % Genero nombre de canal a levantar
    puerto_canal_custom_aux = horzcat(puerto_canal_custom, ...
        '-',num2str(i));

    % Traduzco custom_channel_name a native_channel_name
    traduccion = strcmp(puerto_canal_custom_aux, ...
        {amplifier_channels(:).custom_channel_name});
    puerto_canal = amplifier_channels(traduccion).native_channel_name;

    % Levanto el canal de interes
    raw = read_INTAN_channel(directorio, puerto_canal, ...
        amplifier_channels);

    % Aplica el filtro para LFP
    LFP = filtfilt(filt_LFP, raw);

    % Adjunto LFP de este canal a la "lista" de LFP
    LFP_canales = [LFP_canales,LFP];
    
end

% Promedio los 4 canales de los tetrodos para obtener un LFP promediado
LFP_tetrodo = mean(LFP_canales, 2);

end

