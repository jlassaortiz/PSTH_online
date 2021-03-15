function dict_score = score_calculator(id_BOS, estimulos, rasters, frequency_parameters)
% Calcula la int_normalizada y la correlacion del protocolo experimental
%   
%   Entradas:
% 
%   Salidas:
%   dict_score = (struct)
%   dict_score.name = (string) nombre del estimulo
%   dict_score.int_norm = (double) integral de spikes durante la
%   presentacion del estimulo restandole la actividad espontanea y
%   estandarizando por la inetegral de spikes del BOS.
%   dict_score.corr = (double) correlacion de pearson de la slinding window
%   del estimulo con la sliding window del BOS

t_window = 0.015; % 15 ms
step = 0.001; % 1 ms

% Calculo la sw del BOS para poder hacer correlaciones con el resto
[sw_data_BOS, sw_times_BOS] = sliding_window(rasters(id_BOS).spikes_norm, frequency_parameters.amplifier_sample_rate, ...
        t_window, step);
    
% Conservo solo la seccion donde se presenta el estimulo auditivo
duracion_BOS = length(estimulos(id_BOS).song) / estimulos(id_BOS).freq; % en seg
sw_data_BOS = sw_data_BOS(sw_times_BOS < duracion_BOS);
sw_data_BOS_norm = sw_data_BOS / max(sw_data_BOS);

% Calculo integracion de spikes normalizada del BOS para poder estandarizar
% el resto
integral = sum(rasters(id_BOS).spikes_norm < duracion_BOS * frequency_parameters.amplifier_sample_rate);
ruido    = sum(rasters(id_BOS).spikes_norm > duracion_BOS * frequency_parameters.amplifier_sample_rate & ...
    rasters(id_BOS).spikes_norm < duracion_BOS * 2 * frequency_parameters.amplifier_sample_rate); 
integral_norm_BOS = integral - ruido;

dict_score = struct;

for i = (1:1:length(estimulos)) % para cada estímulo
    
    % Guardo el nombre del estimulo
    id_estimulo = estimulos(i).name;
    
    % Guardo la frecuencia de sampleo y el largo (en seg) de este estimulo
    song_freq = estimulos(i).freq;
    song_len = length(estimulos(i).song) / song_freq; % unidades: seg
    
    % Integracion de spikes normalizada
    integral = sum(rasters(i).spikes_norm < song_len * frequency_parameters.amplifier_sample_rate);
    ruido    = sum(rasters(i).spikes_norm > song_len * frequency_parameters.amplifier_sample_rate & ...
        rasters(i).spikes_norm < song_len * 2 * frequency_parameters.amplifier_sample_rate); 
    integral_norm = (integral - ruido)/integral_norm_BOS;
    
    % Calculo sliding window
    [sw_data, sw_times] = sliding_window(rasters(i).spikes_norm, frequency_parameters.amplifier_sample_rate, ...
        t_window, step);
    
    % Calculo correlación de sw normalizada con la sw normalizada del BOS
    sw_data_norm = sw_data / max(sw_data);
    sw_data_norm = sw_data_norm(sw_times < duracion_BOS);
    correlacion_pearson = corrcoef(sw_data_norm, sw_data_BOS_norm);
    
    % Guardo resultados
    dict_score(i).name = id_estimulo;
    dict_score(i).int_norm = integral_norm;
    dict_score(i).corr = correlacion_pearson(1,2);
end 

end

