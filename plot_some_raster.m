function plot_some_raster(id_estimulos, estimulos, rasters, frequency_parameters, tiempo_file, ntrials, puerto_canal, thr, directorio)
% 
% plot_raster plotea el raster y psth del número de estímulo indicado
%   Detailed explanation goes here

psth_max = ntrials;

figure()

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

for i = id_estimulos % para cada est�mulo
    
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
    title(strcat(string(i), " - ",estimulos(i).name), 'Interpreter','None', 'FontSize',6)

    % psth
    j = j + 1;
    h(j) = subplot(n, m, [p + 3, p + 6]);
    histogram(rasters(i).spikes_norm * 1000/frequency_parameters.amplifier_sample_rate , ...
        (1000/frequency_parameters.amplifier_sample_rate) * (-1000:(0.015*frequency_parameters.amplifier_sample_rate):(tiempo_file*frequency_parameters.amplifier_sample_rate)) );
    ylim([0 psth_max]);
    xlim([0 tiempo_file * 1000]);
    
    % Guardo la frecuencia de sampleo y el largo (en seg) de este estimulo
    song_freq = estimulos(i).freq;
    song_len = length(estimulos(i).song) / song_freq; % unidades: seg
    
    % Integracion de spikes normalizada
    integral = sum(rasters(i).spikes_norm < song_len * frequency_parameters.amplifier_sample_rate);
    ruido    = sum(rasters(i).spikes_norm > song_len * frequency_parameters.amplifier_sample_rate & ...
        rasters(i).spikes_norm < song_len * 2 * frequency_parameters.amplifier_sample_rate); 
    integral_norm = integral - ruido;
    integral_text = strcat('integral_norm: ', string(integral_norm));
    
    title(integral_text , 'Interpreter','None')

    % raster
    j = j + 1;
    h(j) = subplot(n, m, [p + 9, p + 12]);
    plot((1000/frequency_parameters.amplifier_sample_rate) * rasters(i).spikes_norm, rasters(i).trials_id, '.')  
    xlim([0 tiempo_file * 1000])
    ylim([0 ntrials + 1])
end 

% Linkeo eje x (no se pueden hacer varios links independientes juntos)
linkaxes(h, 'x');

% Titulo general
sgtitle({datestr(now, 'yyyy-mm-dd'); ...
    string(directorio) ; ...
    strcat(string(puerto_canal), "  " , string(thr), "uV", "  ntrials:", string(ntrials), "  t_inter_estimulo:", string(tiempo_file)) }, 'Interpreter','None')

end

