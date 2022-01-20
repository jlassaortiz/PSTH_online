function plot_some_raster_LFP_1channel(id_estimulos, id_BOS, estimulos, rasters, ...
    frequency_parameters, tiempo_file, ntrials, puerto_canal, thr, ...
    directorio, spike_times)

% Plotea el psth y LFP del numero de estimulos indicados
%   
%   INPUTS:
%   id_estimulos: (matris) lista de numeros id de los estimulos a graficar
%   id_BOS: (num) numero de id del estimulos BOS
%   estimulos: (struct) objeto donde estan los audios de los estimulos
%   rasters: (struct) objeto donde estan los rasters de cada estimulo
%   frequency_parameters: (struct) objeto donde esta el sampling rate INTAN
%   tiempo_file: (num) tiempo en segundos que duran los estimulos
%   ntrials: (num) numero de trials que se pasa cada estimulo
%   puerto_canal: (string) nombre custom del canal que se analiza
%   thr: (num) umbral de deteccion de spikes
%   directorio: (string) directorio donde estan los datos
%   spike_times: (matris) lista de tiempos de todos los spikes detectados
%
%   OUTPUTS:
%

% Ventana y step del sliding window
t_window = 0.015; % 15 ms
step = 0.001; % 1 ms

% Busco el maximo de los psth para saber que valor poner de ylim
hist_aux = histcounts(rasters(id_BOS).spikes_norm * 1000/...
    frequency_parameters.amplifier_sample_rate , ...
        (1000/frequency_parameters.amplifier_sample_rate) * ...
        (-1000:(0.015 * frequency_parameters.amplifier_sample_rate): ...
        (tiempo_file*frequency_parameters.amplifier_sample_rate)) );

% ylim de los psth es un 20% mas que el maximo del BOS
psth_max = max(hist_aux) * 1.2; 

% Inicializo figura
figure()

% Para no ver tanta actividad espontanea grafico 75% mas de lo que dura BOS
limite_eje_x = (1000 * length(estimulos(id_BOS).song) / ...
    estimulos(id_BOS).freq) * 1.75;

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

% Calculo scores de cada estimulo
dict_score = score_calculator(id_BOS, estimulos, frequency_parameters, ...
    spike_times, ntrials);

% Calculo la sw del BOS para poder hacer correlaciones con el resto
[sw_data_BOS, sw_times_BOS] = sliding_window(rasters(id_BOS).spikes_norm, ...
    frequency_parameters.amplifier_sample_rate, ...
        t_window, step);
    
% Conservo solo la seccion donde se presenta el estimulo auditivo
duracion_BOS = length(estimulos(id_BOS).song) / estimulos(id_BOS).freq; 
sw_data_BOS = sw_data_BOS(sw_times_BOS < duracion_BOS);
sw_times_BOS = sw_times_BOS(sw_times_BOS < duracion_BOS);

for i = id_estimulos % para cada estímulo
    
    k = k + 1;
    
    if mod(k, 3) == 1
        p = (k - 1)/3 * 15 + 1;
        
    elseif mod(k, 3) == 2
        p = (k-2)/3 * 15 + 2;
        
    else
        p = ((k / 3) - 1) * 15 + 3;
    end
    
    % SONIDO
    j = j + 1;
    if length(id_estimulos) == 1
        h(1) = subplot(5,1,1);
    else
        h(j) = subplot(n, m , p);
    end 
    plot(1000/estimulos(i).freq * (0:1:(length(estimulos(i).song) -1)), ... 
        estimulos(i).song,'black')
    hold on;
    line([0 tiempo_file*1000],[0 0],'color',[0 0 0]);
    xlim([0 limite_eje_x])
    xticks([])
    
    title(strcat(string(i), " - ",estimulos(i).name), 'Interpreter','None', ...
        'FontSize', 6)

    % PSTH
    j = j + 1;
    if length(id_estimulos) == 1
        h(2) = subplot(5,1, [2,3]);
    else
        h(j) = subplot(n, m, [p + 3, p + 6]);
    end
    
%     % Ploteo histograma del PSTH
%     histogram(rasters(i).spikes_norm * 1000/ ...
%         frequency_parameters.amplifier_sample_rate ,...
%         (1000/frequency_parameters.amplifier_sample_rate) * ...
%         (-1000:(0.015 * frequency_parameters.amplifier_sample_rate): ...
%         (tiempo_file*frequency_parameters.amplifier_sample_rate)) );
%     hold on;
    
    % Calculo y ploteo sliding window para cada estimulo
    [sw_data, sw_times] = sliding_window(rasters(i).spikes_norm, ...
        frequency_parameters.amplifier_sample_rate, ...
        t_window, step);
    plot(sw_times * 1000, sw_data, '-b');
    hold on;
    plot(sw_times_BOS * 1000 , sw_data_BOS, '-r');
    ylim([0 psth_max]);
    xlim([0 limite_eje_x]);
    
    % Integracion de spikes normalizada
    integral_text = strcat('Integral_norm : ', ...
        string(dict_score(i).int_norm));
    
    % Calculo correlación de sw normalizada con la sw normalizada del BOS
    R2_text = strcat(' Coef Pearson sw_BOS_norm : ' , ...
        string(round(dict_score(i).corr, 2)));

    % Escribo en el titulo los valores de integral y correlacion de sw
    % normalizadas
    title(strcat(integral_text, ' / ' , R2_text) , 'Interpreter','None')

    
    % LFP promediado
    j = j + 1;
    if length(id_estimulos) == 1
        h(3) = subplot(5,1,[4,5]);
    else 
        h(j) = subplot(n, m, [p + 9, p + 12]);
    end
    
    plot((1:1: tiempo_file * frequency_parameters.amplifier_sample_rate) *...
        1000 / frequency_parameters.amplifier_sample_rate, ...
        rasters(i).LFP_promedio)
    xlim([0 limite_eje_x])
%     ylim([-2000 2000])
    xticks([])
end 

% Linkeo eje x (no se pueden hacer varios links independientes juntos)
linkaxes(h, 'x');

% Titulo general
sgtitle({datestr(now, 'yyyy-mm-dd'); ...
    string(directorio) ; ...
    strcat(string(puerto_canal), "  " , string(thr), "uV", "  ntrials:", ...
    string(ntrials), "  t_inter_estimulo:", string(tiempo_file)) }, ...
    'Interpreter','None')

end

