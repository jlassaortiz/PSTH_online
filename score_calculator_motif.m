function [int_norm, corr, sw_motif] = score_calculator_motif(...
    spikes_norm, ...
    spikes_norm_BOS, ...
    duracion_BOS, ...
    sr_spikes, ...
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
%   int_norm = (double) integral de spikes durante la
%   presentacion del estimulo restandole la actividad espontanea y
%   estandarizando por la inetegral de spikes del BOS.
%   corr = (double) correlacion de pearson de la slinding window
%   del estimulo con la sliding window del BOS
%   sw = (nx2 double) tiempos (:,1) y valores del sw (:,2)


% Tamano de ventana y step de la sliding window
t_window = 0.015; % 15 ms
step = 0.001; % 1 ms

% Calculo la sw del BOS para poder hacer correlaciones con el resto
[sw_data_BOS, sw_times_BOS] = sliding_window(spikes_norm_BOS, sr_spikes, ...
    t_window, step, tiempo_file);
    
% Conservo solo la seccion donde se presenta el estimulo auditivo
% duracion_BOS = length(rasters(id_BOS).song) / rasters(id_BOS).freq; % en seg
sw_data_BOS = sw_data_BOS(sw_times_BOS < duracion_BOS);
sw_data_BOS_norm = sw_data_BOS / max(sw_data_BOS);

% Calculo integracion de spikes normalizada del BOS para poder estandarizar
% el resto
integral = sum(spikes_norm_BOS < duracion_BOS * sr_spikes);

% RUIDO = determinado desde la ACTIVIDAD ESPONTANEA (primemos 60 seg.)
ruido = sum(spike_times < 60 * sr_spikes);
ruido = (ntrials*ruido*3) * duracion_BOS/60; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% HARCODEO !

% Resto ruido a la integral del BOS
integral_norm_BOS = integral - ruido;
    
% Integracion de spikes normalizada
integral = sum(spikes_norm < duracion_BOS * sr_spikes);
integral = integral - ruido;
integral_norm = (integral)/integral_norm_BOS;

% Calculo sliding window
[sw_data, sw_times] = sliding_window(spikes_norm, sr_spikes, ...
    t_window, step, tiempo_file);

% Calculo correlación de sw normalizada con la sw normalizada del BOS
sw_data_norm = sw_data / max(sw_data);
sw_data_norm = sw_data_norm(sw_times < duracion_BOS);
correlacion_pearson = corrcoef(sw_data_norm, sw_data_BOS_norm);

% Guardo resultados
int_norm = integral_norm;
corr = correlacion_pearson(1,2);
sw_motif = horzcat(sw_times, sw_data);
sw_motif = sw_motif(sw_times < duracion_BOS,:);

% act_esp = ruido;
 

end


