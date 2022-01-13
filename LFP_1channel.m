function LFP = LFP_1channel(directorio,amplifier_channels, ...
    frequency_parameters, puerto_canal_custom)

% Calcula LFP de un canal especificado
%
%   ! NECESITA QUE ANTES SE CORRA read_Intan_RHD2000_file
%
%   INPUTS
%   directorio: (str) directorio donde esta el archivo INTAN
%   amplifier_channels: salida read_Intan_RHD2000_file
%   frequency_parameters: salida read_Intan_RHD2000_file
%   puerto_canal_custom: (str) nombre custom del canal a levantar
%
%   OUTPUT
%   LFP: (matris) vector fila con señal filtrada


% Define el filtro para LFP
filt_LFP = designfilt('lowpassiir','DesignMethod','butter',...
    'FilterOrder', 4,'HalfPowerFrequency', 300, ...
    'SampleRate', frequency_parameters.amplifier_sample_rate);

% Traduzco custom_channel_name a native_channel_name
traduccion = strcmp(puerto_canal_custom, ...
    {amplifier_channels(:).custom_channel_name});
puerto_canal = amplifier_channels(traduccion).native_channel_name;

% Levanto el canal de interes
raw = read_INTAN_channel(directorio, puerto_canal, ...
    amplifier_channels);

% Aplica el filtro para LFP
LFP = filtfilt(filt_LFP, raw);

end

