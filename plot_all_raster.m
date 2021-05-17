function plot_all_raster(estimulos, id_BOS, rasters, frequency_parameters, tiempo_file, ntrials, puerto_canal, thr, directorio)

% Plotea todos los estimulos
%   No importa la cantidad de estimulos


% Busco el maximo de los psth
hist_aux = histcounts(rasters(id_BOS).spikes_norm * 1000/frequency_parameters.amplifier_sample_rate , ...
        (1000/frequency_parameters.amplifier_sample_rate) * (-1000:(0.015* ... 
        frequency_parameters.amplifier_sample_rate):(tiempo_file*frequency_parameters.amplifier_sample_rate)) );
psth_max = max(hist_aux) * 1.2; % ylim de los psth es un 20% mas que el maximo del BOS

% Inicializo figura
figure()

% Para no ver tanta actividad espontanea grafico 75% mas de lo que dura el BOS
limite_eje_x = (1000 * length(estimulos(id_BOS).song) / estimulos(id_BOS).freq) * 1.75;

n = 5 * round(length(estimulos)/2);
m = 2;
j = 0;

% Para cada estimulo
for i = (1:1: length(rasters))
    
    if mod(i, 2) == 1
        p = (i - 1)/2 * 10 + 1;
    else
        p = ((i / 2) - 1) * 10 + 2;
    end
        
    % sonido
    j = j + 1;
    h(j) = subplot(n, m , p);
    plot(1000/estimulos(i).freq * (0:1:(length(estimulos(i).song) -1)), estimulos(i).song,'black')
    hold on;
    line([0 tiempo_file*1000],[0 0],'color',[0 0 0]);
    xlim([0 limite_eje_x])
    title(strcat(string(i), " - ",estimulos(i).name), 'Interpreter','None')
    xticks([]);
    
    % psth
    j = j + 1;
    h(j) = subplot(n, m, [p + 2, p + 4]);
    histogram(rasters(i).spikes_norm * 1000/frequency_parameters.amplifier_sample_rate , ...
        (1000/frequency_parameters.amplifier_sample_rate) * (-1000:(0.015*frequency_parameters.amplifier_sample_rate):(tiempo_file*frequency_parameters.amplifier_sample_rate)) );
    ylim([0 psth_max]);
    xlim([0 limite_eje_x]);
    xticks([]);
    
    % raster
    j = j + 1;
    h(j) = subplot(n, m, [p + 6, p + 8]);
    plot((1000/frequency_parameters.amplifier_sample_rate) * rasters(i).spikes_norm, rasters(i).trials_id, '.')  
    xlim([0 limite_eje_x])
    ylim([0 ntrials + 1])
end

% Linkeo eje x (no se pueden hacer varios links independientes juntos)
linkaxes(h, 'x');

% Titulo general
sgtitle({datestr(now, 'yyyy-mm-dd'); ...
    string(directorio) ; ...
    strcat(string(puerto_canal), "  " , string(thr), "uV", "  ntrials:", string(ntrials), "  t_inter_estimulo:", string(tiempo_file)) }, 'Interpreter','None')

end

