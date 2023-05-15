function LFP = LFP_1channel(raw, frequency_parameters, ...
    downsample_sr, banda, b_inf, b_sup)

% Calcula LFP de un canal especificado (FILTRA)
%
%   ! NECESITA QUE ANTES SE CORRA read_Intan_RHD2000_file
%
%   INPUTS
%   frequency_parameters: salida read_Intan_RHD2000_file
%
%   OUTPUT
%   LFP: (matris) vector fila con señal filtrada

% Define el filtro para LFP
filt_LFP = designfilt('lowpassiir','DesignMethod','butter',...
    'FilterOrder', 4,'HalfPowerFrequency', 300, ...
    'SampleRate', frequency_parameters.amplifier_sample_rate);

% Aplica el filtro para LFP
LFP = filtfilt(filt_LFP, raw);

% Downsampleo LFP pq al pedo tenerlo a 30kHz, lo llevo a downsample_sr
factor_conv = int8(frequency_parameters.amplifier_sample_rate / downsample_sr);
LFP = downsample(LFP, factor_conv, uint64(factor_conv-1));

% Filtro banda en particular
if banda
    LFP = filt_and_normalize(LFP, b_inf, b_sup, downsample_sr);
end

end
