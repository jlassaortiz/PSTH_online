function plot_some_raster_v2(id_estimulos, id_BOS, estimulos, dict_score, ...
    frequency_parameters, tiempo_file, ntrials, puerto_canal, thr, directorio,...
    spike_times, acumulo_motivo, motif_t1, motif_dur)

% plot_raster plotea el raster y psth del numero de estimulos indicados
%   Detailed explanation goes here

t_window = 0.015; % 15 ms
step = 0.001; % 1 ms
fz = 16;

sr_spikes = frequency_parameters.amplifier_sample_rate;

% Busco el maximo de los psth para saber que valor poner de ylim
if acumulo_motivo == 0
    hist_aux = histcounts(dict_score(id_BOS).spikes_norm * 1000/sr_spikes , ...
            (1000/sr_spikes) * (-1000:(0.015*sr_spikes): ...
            (tiempo_file*sr_spikes)) );
    psth_max = max(hist_aux) * 1.2; % ylim de los psth es un 20% mas que el maximo del BOS
else
    hist_aux = histcounts(dict_score(id_BOS).spikes_norm_motif * 1000/sr_spikes , ...
            (1000/sr_spikes) * (-1000:(0.015*sr_spikes): ...
            (tiempo_file*sr_spikes)) );
    psth_max = max(hist_aux) * 1.2; % ylim de los psth es un 20% mas que el maximo del BOS
end

% Inicializo figura
figure()

% Para no ver tanta actividad espontanea grafico 75% mas de lo que dura el BOS
limite_eje_x = (1000 * length(estimulos(id_BOS).song) / estimulos(id_BOS).freq) * 1.75;

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
% dict_score = score_calculator(id_BOS, rasters, frequency_parameters, ...
%     spike_times, ntrials, tiempo_file);
% 
% % Calculo la sw del BOS para poder hacer correlaciones con el resto
% [sw_data_BOS, sw_times_BOS] = sliding_window(rasters(id_BOS).spikes_norm, ...
%     sr_spikes, ...
%         t_window, step, tiempo_file);
    
% Conservo solo la seccion donde se presenta el estimulo auditivo
if acumulo_motivo == 0
    duracion_BOS = length(estimulos(id_BOS).song) / estimulos(id_BOS).freq; 
    sw_data_BOS = dict_score(id_BOS).sw(dict_score(id_BOS).sw(:,1) < duracion_BOS,2);
    sw_times_BOS = dict_score(id_BOS).sw(dict_score(id_BOS).sw(:,1) < duracion_BOS,1);
else
    sw_data_BOS = dict_score(id_BOS).sw_motif(:,2);
    sw_times_BOS = dict_score(id_BOS).sw_motif(:,1);     
end

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
    h(j) = subplot(n, m , p);
    plot(1000/estimulos(i).freq * (0:1:(length(estimulos(i).song) -1)), ...
        estimulos(i).song,'black')
    hold on;
    line([0 tiempo_file*1000],[0 0],'color',[0 0 0]);
    if acumulo_motivo == 0
        xlim([0 limite_eje_x])
    else
        xlim([motif_t1*1000 (motif_t1 + motif_dur)*1000]);
    end
    xticks([])
    set(gca, 'FontSize', fz)
    
    title(strcat(string(i), " - ",estimulos(i).name), ...
        'Interpreter','None', 'FontSize', 8)

    % PSTH
    j = j + 1;
    h(j) = subplot(n, m, [p + 3, p + 6]);
    if acumulo_motivo == 0
        % ploteo histograma
        histogram( dict_score(i).spikes_norm * 1000/sr_spikes , ...
        (1000/sr_spikes) * (-1000:(0.015* sr_spikes): ...
        (tiempo_file * sr_spikes)) ); 
        ylim([0 psth_max]);
        xlim([0 limite_eje_x]);
        hold on
        
        % ploteo sw del bos y del estimulo
        sw_times = dict_score(i).sw(:,1);
        sw_data = dict_score(i).sw(:,2);
        plot(sw_times * 1000, sw_data, '-b');
        plot(sw_times_BOS * 1000 , sw_data_BOS, '-r');
        set(gca, 'FontSize', fz)
        
        % Integracion de spikes normalizada
        integral_text = strcat('Integral_norm : ', ...
            string(dict_score(i).int_norm));
    
        % Calculo correlación de sw normalizada con la sw normalizada del BOS
        R2_text = strcat(' Coef Pearson sw_BOS_norm : ' , ...
            string(round(dict_score(i).corr, 2)));
    else
        % ploteo histograma del motivo
        histogram( (dict_score(i).spikes_norm_motif + motif_t1*sr_spikes)...
            * 1000/sr_spikes , ...
        (1000/sr_spikes) * (-1000:(0.015* sr_spikes): ...
        (tiempo_file * sr_spikes)) );
        hold on
        
        % ploteo sw (motivo acumulado) del bos y del estimulo
        sw_times = dict_score(i).sw_motif(:,1) + motif_t1;
        sw_data = dict_score(i).sw_motif(:,2);
        plot(sw_times * 1000, sw_data, '-b');
        plot((sw_times_BOS + motif_t1) * 1000 , sw_data_BOS, '-r');
        set(gca, 'FontSize', fz)
        
        ylim([0 psth_max]);
        xlim([motif_t1 (motif_t1 + motif_dur)]);
       
        % Integracion de spikes normalizada
        integral_text = strcat('Integral_norm : ', ...
            string(dict_score(i).int_norm_motif));
    
        % Calculo correlación de sw normalizada con la sw normalizada del BOS
        R2_text = strcat(' Coef Pearson sw_BOS_norm : ' , ...
            string(round(dict_score(i).corr_motif, 2)));
    end
    
    % Escribo en el titulo los valores de integral y correlacion de sw
    % normalizadas
    title(strcat(integral_text, ' / ' , R2_text) , 'Interpreter','None',...
        'Fontsize', 14)

    % RASTER
    j = j + 1;
    h(j) = subplot(n, m, [p + 9, p + 12]);
    plot((1000/sr_spikes) * ...
        dict_score(i).spikes_norm, dict_score(i).trials_id, '.')
    if acumulo_motivo == 0
        xlim([0 limite_eje_x])
    else
        xlim([motif_t1*1000 (motif_t1 + motif_dur)*1000]);
    end
    ylim([0 ntrials + 1])
    xticks([])
    set(gca, 'FontSize', fz)
end 

% Linkeo eje x (no se pueden hacer varios links independientes juntos)
linkaxes(h, 'x');
if acumulo_motivo == 0
    xlim([0 limite_eje_x])
else
    xlim([motif_t1*1000 (motif_t1 + motif_dur)*1000]);
end


% sgtitle({datestr(now, 'yyyy-mm-dd'); ...
%      string(directorio) ; ...
%      strcat(string(puerto_canal), "  " , string(thr), "uV", "  ntrials:", string(ntrials),  ...
%      "  t_inter_estimulo:", string(tiempo_file)) }, 'Interpreter','None')

suptitle({datestr(now, 'yyyy-mm-dd'); ...
    string(directorio) ; ...
    strcat(string(puerto_canal), "  " , string(thr), "uV", "  ntrials:", string(ntrials),  ...
    "  t_inter_estimulo:", string(tiempo_file)) })

end

