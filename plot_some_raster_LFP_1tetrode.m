function plot_some_raster_LFP_1tetrode(id_estimulos, id_BOS, ...
    estimulos_tetrodos, frequency_parameters, tiempo_file, ntrials, thr,...
    directorio, spike_times)

% Plotea el psth y LFP del numero de estimulos indicados
%   
%   INPUTS:
%   id_estimulos: (matris) lista de numeros id de los estimulos a graficar
%   id_BOS: (num) numero de id del estimulos BOS
%   estimulos_tetrodos: (struct) donde audios y rasters de cada estimulo
%   frequency_parameters: (struct) objeto donde esta el sampling rate INTAN
%   tiempo_file: (num) tiempo en segundos que duran los estimulos
%   ntrials: (num) numero de trials que se pasa cada estimulo
%   puerto_canal: (string) nombre custom del canal que se analiza
%   thr: (num) umbral de deteccion de spikes
%   directorio: (string) directorio donde estan los datos
%   spike_times: (matris) lista de tiempos de todos los spikes detectados
%
%   OUTPUTS:


% Ventana y step del sliding window
t_window = 0.015; % 15 ms
step = 0.001; % 1 ms

% Tamaï¿½o del sw una vez que recorre todo el canto
size_sw = 0;
tf = t_window;
while tf <= tiempo_file
   size_sw = size_sw + 1;
   tf = tf + step;
end

% Inicializo y genero vector de tiempos del PSTH
t_PSTH = zeros(size_sw,1); % va a estar en SEGUNDOS
ti = 0;
tf = t_window;
for i = (1:1:size_sw)
    
    t_PSTH(i) = (ti + tf)/2;
    
    ti = ti + step;
    tf = tf + step;
end

% Busco el maximo de los psth para saber que valor poner de ylim
hist_aux = histcounts(estimulos_tetrodos(1).canal(id_BOS).spikes_norm * ...
    1000/frequency_parameters.amplifier_sample_rate , ...
        (1000/frequency_parameters.amplifier_sample_rate) * ...
        (-1000:(0.015 * frequency_parameters.amplifier_sample_rate): ...
        (tiempo_file*frequency_parameters.amplifier_sample_rate)) );

% ylim de los psth es un 20% mas que el maximo del BOS
psth_max = max(hist_aux) * 1.2;

% Inicializo figura
figure()

% Para no ver tanta actividad espontanea grafico 75% mas de lo que dura BOS
limite_eje_x = (1000* length(estimulos_tetrodos(1).canal(id_BOS).song) /...
    estimulos_tetrodos(1).canal(id_BOS).freq) * 1.75;

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

% % Calculo scores de cada estimulo
% dict_score = score_calculator(id_BOS, estimulos, frequency_parameters,...
%     spike_times, ntrials);

% % Calculo la sw del BOS para graficarla y comparar con el resto
% [sw_data_BOS, sw_times_BOS] = sliding_window(...
%     estimulos_tetrodos(1).canal(id_BOS).spikes_norm, ...
%     frequency_parameters.amplifier_sample_rate, t_window, step);
%     
% % Conservo solo la seccion donde se presenta el estimulo auditivo
% duracion_BOS = length(estimulos_tetrodos(1).canal(id_BOS).song) / ...
%     estimulos_tetrodos(1).canal(id_BOS).freq; 
% sw_data_BOS = sw_data_BOS(sw_times_BOS < duracion_BOS);
% sw_times_BOS = sw_times_BOS(sw_times_BOS < duracion_BOS);

for i = id_estimulos % para cada estï¿½mulo
    
    k = k + 1;
    
    if mod(k, 3) == 1
        p = (k - 1)/3 * 15 + 1;
        
    elseif mod(k, 3) == 2
        p = (k-2)/3 * 15 + 2;
        
    else
        p = ((k / 3) - 1) * 15 + 3;
    end
    
    % Inicializo vector donde guardo señales y correlaciones
    PSTH_avgTetrodo = zeros(size_sw, length(estimulos_tetrodos));
    corr_error_PSTH = zeros(2, length(estimulos_tetrodos));
    
    LFP_avgTetrodo = ones(...
        length(estimulos_tetrodos(1).canal(1).LFP_promedio), ...
        length(estimulos_tetrodos));
    corr_error_LFP = ones(2, length(estimulos_tetrodos)); 
    
    
    % SONIDO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    j = j + 1;
    if length(id_estimulos) == 1
        h(1) = subplot(5,1,1);
    else
        h(j) = subplot(n, m , p);
    end
    plot(1000/estimulos_tetrodos(1).canal(i).freq * ...
        (0:1:(length(estimulos_tetrodos(1).canal(i).song) -1)), ...
        estimulos_tetrodos(1).canal(i).song,'black')
    hold on;
    
    line([0 tiempo_file*1000],[0 0],'color',[0 0 0]);
    xlim([0 limite_eje_x])
    xticks([])
    
    title(strcat(string(i), " - ",estimulos_tetrodos(1).canal(i).name), ...
        'Interpreter','None','FontSize', 6)

    % PSTH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    j = j + 1;
    if length(id_estimulos) == 1
        h(2) = subplot(5,1, [2,3]);
    else
        h(j) = subplot(n, m, [p + 3, p + 6]);
    end
    
    % Para cada canal
    for c = (1:1:length(estimulos_tetrodos))
        
        estimulos = estimulos_tetrodos(c).canal;

        % Calculo y ploteo sliding window para cada estimulo
        [sw_data, sw_times] = sliding_window(estimulos(i).spikes_norm, ...
            frequency_parameters.amplifier_sample_rate, ...
            t_window, step);
        
        % Ploteo cada PSTH de cada canal
        plot(sw_times * 1000, sw_data);
        hold on;

        % Guardo PSTH de cada canal
        PSTH_avgTetrodo(1:length(sw_data),c) = sw_data;
    end
        
    % Promedio PSTH de todos los canales de este estimulo
    PSTH_avgTetrodo_mean = mean(PSTH_avgTetrodo, 2);
    
    % Calculo correlaciones y error cuadratico medio del PSTH de cada canal
    for c = (1:1:length(estimulos_tetrodos))
        corr_error_PSTH(1,c) = corr(PSTH_avgTetrodo(:,c), PSTH_avgTetrodo_mean);
        corr_error_PSTH(2,c) = sum(abs(PSTH_avgTetrodo(:,c) - ...
            PSTH_avgTetrodo_mean)) / length(PSTH_avgTetrodo_mean);
    end 
    corr_error_PSTH = round(corr_error_PSTH, 2);
    
    % Ploteo PSTH promedio
    plot(t_PSTH * 1000, PSTH_avgTetrodo_mean, '-r', 'LineWidth', 2);
    ylim([0 psth_max]);
    xlim([0 limite_eje_x]);
    title(strcat('corr chann: ', num2str(corr_error_PSTH(1,:)),...
        ' MAE chann: ', num2str(corr_error_PSTH(2, :))), ...
        'FontSize', 8)
    
    % LFP promediado %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    j = j + 1;
    if length(id_estimulos) == 1
        h(3) = subplot(5,1,[4,5]);
    else 
        h(j) = subplot(n, m, [p + 9, p + 12]);
    end
    
    % Para cada canal
    for c = (1:1:length(estimulos_tetrodos))
        
        estimulos = estimulos_tetrodos(c).canal;      
        
        % Ploteo LFP de cada canal
        plot((1:1: ...
            tiempo_file * frequency_parameters.amplifier_sample_rate) * ...
            1000 / frequency_parameters.amplifier_sample_rate, ...
            estimulos(i).LFP_promedio)
        hold on;
        
        % Guardo LFP de cada canal
        LFP_avgTetrodo(:,c) = estimulos(i).LFP_promedio;
        
    end
  
    % Promedio LFP de todos los canales del tetrodo de este estimulo
    LFP_avgTetrodo_mean = mean(LFP_avgTetrodo, 2);
    lfp_max = max(LFP_avgTetrodo_mean);
    
    % Calculo correlaciones y error cuadratico medio del PSTH de cada canal
    for c = (1:1:length(estimulos_tetrodos))
        corr_error_LFP(1,c) = corr(LFP_avgTetrodo(:,c), LFP_avgTetrodo_mean);
        corr_error_LFP(2,c) = sum(abs(LFP_avgTetrodo(:,c) - ...
            LFP_avgTetrodo_mean)) / length(LFP_avgTetrodo_mean);
    end 
    corr_error_LFP = round(corr_error_LFP, 2);
    
    % Calculo tiempos de LFP promedio
    t_LFP = ( 1:1:length(LFP_avgTetrodo_mean) ) * ...
        1000/frequency_parameters.amplifier_sample_rate;
    
    % Ploteo LFP promedio del tetrodo
    plot(t_LFP, LFP_avgTetrodo_mean, '-b', 'LineWidth', 2)
    title(strcat('corr chann: ', num2str(corr_error_LFP(1,:)), ...
    ' MAE chann: ', num2str(corr_error_LFP(2, :))),...
    'FontSize', 8, 'Position', [limite_eje_x/2, lfp_max*(0.95), 0])
    
    xlim([0 limite_eje_x])
%     ylim([-2000 2000])
    xticks([])
end 

% Linkeo eje x (no se pueden hacer varios links independientes juntos)
linkaxes(h, 'x');

% % Titulo general
% sgtitle({datestr(now, 'yyyy-mm-dd'); ...
%     string(directorio) ; ...
%     strcat(string(thr), "uV", ...
%     "  ntrials:", string(ntrials), ...
%     "  t_inter_estimulo:", string(tiempo_file)) }, ...
%     'Interpreter','None')

end

