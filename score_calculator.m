function rasters = score_calculator(id_BOS, rasters, frequency_parameters, ...
    spike_times, ntrials, tiempo_file)
% Calcula la int_normalizada y la correlacion del protocolo experimental
%   
%   Entradas:
%   id_BOS = (double) numero de orden que le corresponde al BOS en el
%   objeto "estimulos" generado por la funcion carga_songs.m
%   rasters = (struct) diccionarios con los rasters separados por estimulo
%   frecuency_parameters = (struct) diccionario generado por
%   read_Intan_RHD2000_file.m con info del amplificador y los datos grabados
% 
%   Salidas:
%   dict_score = (struct)
%   dict_score.name = (string) nombre del estimulo
%   dict_score.int_norm = (double) integral de spikes durante la
%   presentacion del estimulo restandole la actividad espontanea y
%   estandarizando por la inetegral de spikes del BOS.
%   dict_score.corr = (double) correlacion de pearson de la slinding window
%   del estimulo con la sliding window del BOS


% Tamano de ventana y step de la sliding window
t_window = 0.015; % 15 ms
step = 0.001; % 1 ms

% Calculo la sw del BOS para poder hacer correlaciones con el resto
[sw_data_BOS, sw_times_BOS] = sliding_window(rasters(id_BOS).spikes_norm, ...
    frequency_parameters.amplifier_sample_rate, ...
        t_window, step, tiempo_file);
    
% Conservo solo la seccion donde se presenta el estimulo auditivo
duracion_BOS = length(rasters(id_BOS).song) / rasters(id_BOS).freq; % en seg
sw_data_BOS = sw_data_BOS(sw_times_BOS < duracion_BOS);
sw_data_BOS_norm = sw_data_BOS / max(sw_data_BOS);

% Calculo integracion de spikes normalizada del BOS para poder estandarizar
% el resto
integral = sum(rasters(id_BOS).spikes_norm < duracion_BOS * ...
    frequency_parameters.amplifier_sample_rate);

% RUIDO = determinado desde la ACTIVIDAD ESPONTANEA (primemos 60 seg.)
ruido = sum(spike_times < 60*frequency_parameters.amplifier_sample_rate);
ruido = (ntrials*ruido) * duracion_BOS/60;

round(ruido)
round(integral)
round(integral / ruido, 2)

integral = integral - ruido;

integral_norm_BOS = integral;

for i = (1:1:length(rasters)) % para cada estímulo
    
%     % Guardo la frecuencia de sampleo y el largo (en seg) de este estimulo
%     song_freq = rasters(i).freq;
%     song_len = length(rasters(i).song) / song_freq; % unidades: seg
    
    % Integracion de spikes normalizada
    integral = sum(rasters(i).spikes_norm < duracion_BOS * ...
        frequency_parameters.amplifier_sample_rate);
    integral = integral - ruido;
    
    integral_norm = (integral)/integral_norm_BOS;
    
    % Calculo sliding window
    [sw_data, sw_times] = sliding_window(rasters(i).spikes_norm, ...
        frequency_parameters.amplifier_sample_rate, ...
        t_window, step, tiempo_file);
    
    % Calculo correlación de sw normalizada con la sw normalizada del BOS
    sw_data_norm = sw_data / max(sw_data);
    sw_data_norm = sw_data_norm(sw_times < duracion_BOS);
    correlacion_pearson = corrcoef(sw_data_norm, sw_data_BOS_norm);
    
    % Guardo resultados
    rasters(i).int_norm = integral_norm;
    rasters(i).corr = correlacion_pearson(1,2);
    rasters(i).act_esp = ruido;
    rasters(i).sw = horzcat(sw_times, sw_data);
end 

end

