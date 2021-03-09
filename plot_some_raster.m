function plot_some_raster(id_estimulos, id_BOS, estimulos, rasters, frequency_parameters, tiempo_file, ntrials, puerto_canal, thr, directorio)

% plot_raster plotea el raster y psth del numero de estimulos indicados
%   Detailed explanation goes here

t_window = 0.015; % 15 ms
step = 0.001; % 1 ms

% Busco el maximo de los psth
hist_aux = histcounts(rasters(id_BOS).spikes_norm * 1000/frequency_parameters.amplifier_sample_rate , ...
        (1000/frequency_parameters.amplifier_sample_rate) * (-1000:(0.015* ... 
        frequency_parameters.amplifier_sample_rate):(tiempo_file*frequency_parameters.amplifier_sample_rate)) );
psth_max = max(hist_aux) * 1.2; % ylim de los psth es un 20% mas que el maximo del BOS

% Inicializo figura
figure()

% Formula para armar grilla segun la cantidad de estimulos a analizar
if mod(length(id_estimulos), 3) == 0
    n = 5 * length(id_estimulos)/3;
    
elseif mod(length(id_estimulos), 3) == 2
    n = 5 * round(length(id_estimulos)/3);
    
else
    n = 5 * (round(length(id_estimulos)/3) + 1);
end
    
m = 3;    
j = 0;
k = 0;

% Calculo la sw del BOS para poder hacer correlaciones con el resto
[sw_data_BOS, sw_times_BOS] = sliding_window(rasters(id_BOS).spikes_norm, frequency_parameters.amplifier_sample_rate, ...
        t_window, step);
    
% Conservo solo la seccion donde se presenta el estimulo auditivo
duracion_BOS = length(estimulos(id_BOS).song) / estimulos(id_BOS).freq; 
sw_data_BOS = sw_data_BOS(sw_times_BOS < duracion_BOS);
sw_times_BOS = sw_times_BOS(sw_times_BOS < duracion_BOS);
sw_data_BOS_norm = sw_data_BOS / max(sw_data_BOS);

for i = id_estimulos % para cada estímulo
    
    k = k + 1;
    
    if mod(k, 3) == 1
        p = (k - 1)/3 * 15 + 1;
        
    elseif mod(k, 3) == 2
        p = (k-2)/3 * 15 + 2;
        
    else
        p = ((k / 3) - 1) * 15 + 3;
    end
    
    % sonido
    j = j + 1;
    h(j) = subplot(n, m , p);
    plot(1000/estimulos(i).freq * (0:1:(length(estimulos(i).song) -1)), estimulos(i).song,'black')
    hold on;
    line([0 tiempo_file*1000],[0 0],'color',[0 0 0]);
    xlim([0 tiempo_file * 1000])
    xticks([])
    
    title(strcat(string(i), " - ",estimulos(i).name), 'Interpreter','None', 'FontSize',6)

    % psth
    j = j + 1;
    h(j) = subplot(n, m, [p + 3, p + 6]);
    hist = histogram(rasters(i).spikes_norm * 1000/frequency_parameters.amplifier_sample_rate , ...
        (1000/frequency_parameters.amplifier_sample_rate) * (-1000:(0.015* ... 
        frequency_parameters.amplifier_sample_rate):(tiempo_file*frequency_parameters.amplifier_sample_rate)) );
    ylim([0 psth_max]);
    xlim([0 tiempo_file * 1000]);
    hold on;
   
    % Busco maximos en cada psth y los grafico
    [maximo, id_max] = max(hist.Values);
    pos_max = hist.BinEdges(id_max);
    pos_max = pos_max + hist.BinWidth /2;
    plot(pos_max, maximo, 'or');
    
    % Guardo la frecuencia de sampleo y el largo (en seg) de este estimulo
    song_freq = estimulos(i).freq;
    song_len = length(estimulos(i).song) / song_freq; % unidades: seg
    
    % Integracion de spikes normalizada
    integral = sum(rasters(i).spikes_norm < song_len * frequency_parameters.amplifier_sample_rate);
    ruido    = sum(rasters(i).spikes_norm > song_len * frequency_parameters.amplifier_sample_rate & ...
        rasters(i).spikes_norm < song_len * 2 * frequency_parameters.amplifier_sample_rate); 
    integral_norm = integral - ruido;
    integral_text = strcat('integral_norm: ', string(integral_norm));
    
    % Calculo sliding window para cada estimulo
    [sw_data, sw_times] = sliding_window(rasters(i).spikes_norm, frequency_parameters.amplifier_sample_rate, ...
        t_window, step);
    plot(sw_times * 1000, sw_data, '-b');
    plot(sw_times_BOS * 1000 , sw_data_BOS, '-r');
    
    % Calculo correlación de sw normalizada con la sw normalizada del BOS
    sw_data_norm = sw_data / max(sw_data);
    sw_data_norm = sw_data_norm(sw_times < duracion_BOS);
    R2 = corrcoef(sw_data_norm, sw_data_BOS_norm);
    R2_text = strcat('R2 sw BOS norm: ' , string(round(R2(1,2), 2)));

    
    % Escribo en el titulo los valores de integral y correlacion de sw
    % normalizadas
    title(strcat(integral_text, ' / ' , R2_text) , 'Interpreter','None')

    % raster
    j = j + 1;
    h(j) = subplot(n, m, [p + 9, p + 12]);
    plot((1000/frequency_parameters.amplifier_sample_rate) * rasters(i).spikes_norm, rasters(i).trials_id, '.')  
    xlim([0 tiempo_file * 1000])
    ylim([0 ntrials + 1])
    xticks([])
end 

% Linkeo eje x (no se pueden hacer varios links independientes juntos)
linkaxes(h, 'x');

% Titulo general
sgtitle({datestr(now, 'yyyy-mm-dd'); ...
    string(directorio) ; ...
    strcat(string(puerto_canal), "  " , string(thr), "uV", "  ntrials:", string(ntrials),  ...
    "  t_inter_estimulo:", string(tiempo_file)) }, 'Interpreter','None')

end

