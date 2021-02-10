function plot_raster(i, estimulos, rasters, frequency_parameters, tiempo_file, ntrials, puerto_canal, thr, directorio)
%
% plot_raster plotea el raster y psth del número de estímulo indicado
%   Detailed explanation goes here

figure()

% sonido
h(1) = subplot(3, 1 , 1);
plot(1000/estimulos(i).freq * (0:1:(length(estimulos(i).song) -1)), estimulos(i).song,'black')
hold on;
line([0 tiempo_file*1000],[0 0],'color',[0 0 0]);
xlim([0 tiempo_file * 1000])
title(strcat(string(i), " - ",estimulos(i).name), 'Interpreter','None')

% psth
h(2) = subplot(3, 1, 2);
histogram(rasters(i).spikes_norm * 1000/frequency_parameters.amplifier_sample_rate , ...
    (1000/frequency_parameters.amplifier_sample_rate) * (-1000:(0.015*frequency_parameters.amplifier_sample_rate):(tiempo_file*frequency_parameters.amplifier_sample_rate)) );
xlim([0 tiempo_file * 1000]);

% raster
h(3) = subplot(3, 1, 3);
plot((1000/frequency_parameters.amplifier_sample_rate) * rasters(i).spikes_norm, rasters(i).trials_id, '.')
xlim([0 tiempo_file * 1000])
ylim([0 ntrials + 1])


% Linkeo eje x (no se pueden hacer varios links independientes juntos)
linkaxes(h, 'x');

% Titulo general
sgtitle({datestr(now, 'yyyy-mm-dd'); ...
    string(directorio) ; ...
    strcat(string(puerto_canal), "  " , string(thr), "uV", "  ntrials:", string(ntrials), "  t_inter_estimulo:", string(tiempo_file)) }, 'Interpreter','None')

end
