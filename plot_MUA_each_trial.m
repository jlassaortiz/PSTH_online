
close all
clear all

% Cargo y defino parametros %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Defino directorio
directorio = input('Directorio: ','s');
directorio = horzcat(directorio , '/');

% Leer info INTAN
read_Intan_RHD2000_file(horzcat(directorio, 'info.rhd'));
clear notes spike_triggers supply_voltage_channels aux_input_channels 

% Carga vector con parametros del analisis de datos
%params_info = dir(horzcat(directorio, '*parametros_protocolo*.txt'));
params = readtable(horzcat(directorio,'/parametros_protocolo.txt'), ...
    'Delimiter','\t');
clear params_info

% Carga vector con parametros del analisis de datos
%params_info = dir(horzcat(directorio, '*parametros_analisis*.txt'));
params_analisis = readtable( ...
    horzcat(directorio,'/parametros_analisis.txt'), ...
    'Delimiter','\t');
clear params_info

% Pregunto si determino canal de manera manual o automatica (esta
% explicitado en parametros que se cargan)
canal_automatico = ...
    input('\nDetermino canal automaticamente? (1 = SI / 0 = NO) : ');

if canal_automatico == 1
    % Cargo valores de puerto-canal del archivo parametros dentro del dir
    puerto = char(params.Puerto);
    canal = params.Canal;
    puerto_canal = horzcat(puerto, '-0', num2str(canal,'%.2d'))
    puerto_canal_custom = puerto_canal;
    clear puerto canal 

else
    % Defino canal a levantar con nombre custom name
    peine = ...
        input('\nDefino canal a levantar (custom name) \n \nPeine (X): ');
    tetrodo = input('\nTetrodo (X): ');
    canal = input('\nCanal (X): ');

    puerto_canal_custom = horzcat('P',num2str(peine),'-', ...
        'T',num2str(tetrodo),'-',num2str(canal));

    % Traduzco custom_channel_name a native_channel_name
    traduccion = strcmp(puerto_canal_custom, ...
        {amplifier_channels(:).custom_channel_name});
    
    puerto_canal = amplifier_channels(traduccion).native_channel_name;
    clear traduccion peine tetrodo canal
end

% Pregunto si ploteo toda la grilla o solo algunos estimulos
plot_grilla = input('\nPloteo toda la grilla? (1 = SI / 0 = NO): ');
if plot_grilla == 0
    grilla_psth = input('\nMatris lineal con numero ID estimulos : ');
end

% Pregunto si el umbral se determina manualmente o automaticamente
thr_automatico = ...
    input('\nBusqueda de thr automatica? (1 = SI / 0 = NO) : ');

% Definimos manualmente un umbral para deteccion de spikes (en uV)
if thr_automatico == 0 
    thr = input('\nThreshold para el threshold cutting (en uV):  ');
end

% Analizo por bandas?
bandas = input('\nFiltro banda particular? (1 = SI / 0 = NO) : ');
if bandas == 1
    b_inf = input('\nLimite inferior (Hz): ');
    b_sup = input('\nLimite superior (Hz): ');
end 

% Cargamos cantidad de trials y tiempo que dura cada uno
ntrials = params.Ntrials
tiempo_file = params.tiempo_entre_estimulos

% Especifico numero de id del BOS y REV
id_BOS = params_analisis.id_bos(1)
id_REV = params_analisis.id_rev(1)

% Cargo orden de la grilla
if plot_grilla == 1
    grilla_psth = str2num(string(params_analisis.grilla_psth(1)))
end

% Genero songs.mat a partir de las canciones
estimulos = carga_songs(directorio);    

% Cargo id_estimulos 
for i = (1:1:length(estimulos))
    estimulos(i).id = params_analisis.orden(i);
    estimulos(i).frec_corte = params_analisis.freq_corte(i);
    estimulos(i).tipo = categorical(params_analisis.tipo_estimulo(i));
end
clear i 

% Levanto seï¿½al neuronal y analizo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Levanto el canal de interes
raw = read_INTAN_channel(directorio, puerto_canal, amplifier_channels);

% Define el filtro para SPIKES y LFP
filt_spikes = designfilt('highpassiir','DesignMethod','butter',...
    'FilterOrder', 4,'HalfPowerFrequency',500,...
    'SampleRate',frequency_parameters.amplifier_sample_rate);

filt_LFP = designfilt('lowpassiir','DesignMethod','butter',...
    'FilterOrder', 4,'HalfPowerFrequency', 300, ...
    'SampleRate', frequency_parameters.amplifier_sample_rate);

% Aplica filtros
raw_filtered = filtfilt(filt_spikes, raw);
LFP = filtfilt(filt_LFP, raw);
sr_lfp = frequency_parameters.amplifier_sample_rate;
downsample_sr = 1000;

% Downsampleo LFP pq al pedo tenerlo a 30kHz, lo llevo a downsample_sr
factor_conv = int8(sr_lfp / downsample_sr);
LFP = downsample(LFP, factor_conv, uint64(factor_conv-1));

% Analizo por bandas?
if bandas == 1
    LFP = filt_and_normalize(LFP, b_inf, b_sup, downsample_sr);
end 

clear filt_spikes

% Genero struct con nombre de los estimulos y el momento de presentacion
estimulos = find_t0s(estimulos, ntrials, tiempo_file, ...
    board_adc_channels, frequency_parameters, directorio, false);

% Definimos umbral de deteccion de spikes
if thr_automatico == 1
    thr = find_thr(raw_filtered,estimulos,tiempo_file,frequency_parameters);
end
clear thr_automatico

% Buscamos spike times (en samples, NO unidades de tiempo) por threshold 
spike_times = find_spike_times(raw_filtered, thr, frequency_parameters);

% Genero objeto con raster de todos los estimulos
estimulos = generate_raster(spike_times, estimulos , tiempo_file, ...
    ntrials, frequency_parameters);

% ntrials = 5;

estimulos = score_calculator(id_BOS, estimulos, ...
        frequency_parameters, spike_times, ntrials, tiempo_file);
    


% Calculo sliding window para cada estimulo
t_window = 0.015; % 15 ms
step = 0.001; % 1 ms
for i = (1:length(estimulos))
    [sw_data, sw_times] = sliding_window(estimulos(i).spikes_norm, ...
        frequency_parameters.amplifier_sample_rate, ...
        t_window, step, tiempo_file);
    psth_sw = [sw_data, sw_times];
    estimulos(i).psth_sw = psth_sw;
    
    % Calculo sw promediando el sw de cada trial
    sw_trial_aux = zeros(length(estimulos(1).spikes_trials(1).sw(:,2)), ntrials);
    for t = (1:ntrials)
        sw_trial_aux(:,t) = estimulos(i).spikes_trials(t).sw(:,1);
    end
    sw_trial_aux = mean(sw_trial_aux, 2);
    psth_sw_trial_avg = horzcat(sw_trial_aux, estimulos(1).spikes_trials(1).sw(:,2));
    estimulos(i).psth_sw_trial_avg = psth_sw_trial_avg;
end
clear i psth_sw


% Agrego a cada trial un numero de id para identificar si el numero de
% spikes de cada trial es mayor o menor que la media
for e = (1:length(estimulos))
    
    % Inicializo vector con numero de spikes por trial
    spikes_por_t = zeros(ntrials,1);
    for i = (1:ntrials)
        n = length(estimulos(e).spikes_trials(i).spikes);
        spikes_por_t(i,1) = n;
        spikes_por_t(i,2) = i;
    end 
    % separo trials si son mayor o menor que el punto medio de la
    % diferencia entre el trial con mas spikes y el de menos spikes
    limite = uint32(max(spikes_por_t(:,1)) -  min(spikes_por_t(:,1)));
    limite = min(spikes_por_t(:,1)) + limite/2; 
    trials_few_spikes = spikes_por_t(:,1) < limite ;
    
    for t = (1:ntrials)
        estimulos(e).spikes_trials(t).trial_id_aux = trials_few_spikes(t,1);
        estimulos(e).spikes_trials(t).limite = limite;
    end 
end 




% Calculo sw para cada sub-grupo de trials (pocos spikes vs muchos spikes)
for e = (1:length(estimulos))
    
    indx_aux = [estimulos(e).spikes_trials.trial_id_aux]' == 1;
    
    few_spikes = sort(...
        vertcat(estimulos(e).spikes_trials(indx_aux).spikes));
    
    many_spikes = sort(...
        vertcat(estimulos(e).spikes_trials(not(indx_aux)).spikes));
    
    % Calculo sliding window para cada estimulo
    t_window = 0.015; % 15 ms
    step = 0.001; % 1 ms

    [sw_few_spikes, sw_t_few_spikes] = sliding_window(few_spikes, ...
        frequency_parameters.amplifier_sample_rate, ...
        t_window, step, tiempo_file);
    
    [sw_many_spikes, sw_t_many_spikes] = sliding_window(many_spikes, ...
        frequency_parameters.amplifier_sample_rate, ...
        t_window, step, tiempo_file);
    
    estimulos(e).sw_few_spikes = [sw_few_spikes, sw_t_few_spikes];
    estimulos(e).sw_many_spikes = [sw_many_spikes, sw_t_many_spikes];
end 
 





% PLOTEO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Ploteo sw de cada trial para cada estímulo
for i = (1:length(grilla_psth))
    
    e = grilla_psth(i);
    leyendas = {};
    
    figure(i)
    t_sound = (1:length(estimulos(e).song))/estimulos(e).freq;
    song = estimulos(e).song;

    % Ploteo sonograma estimulo para tenerlo de referencia
    ax(1) = subplot(8,1,1);
    plot(t_sound', song);
    est_name = estimulos(e).name;
    title({directorio; [est_name ' - ' puerto_canal_custom]}, ...
        'Interpreter', 'None', 'Fontsize', 20)

    % Ploteo cada trial
    ax(2) = subplot(8,1,[2 3 4 5 6 7 8]);
    dist = 13; % distancia entre trazas.
    for t = (1:ntrials)
        if estimulos(e).spikes_trials(t).trial_id_aux == 1
            plot(estimulos(e).spikes_trials(t).sw(:,2), ...
                estimulos(e).spikes_trials(t).sw(:,1) + dist*t, 'b',...
                'Linewidth', 2)
            hold on
            leyendas = [leyendas, ['trial : ', num2str(t)]];
        else
            plot(estimulos(e).spikes_trials(t).sw(:,2), ...
                estimulos(e).spikes_trials(t).sw(:,1) + dist*t, 'r', ...
                'Linewidth', 2)
            hold on
            leyendas = [leyendas, ['trial : ', num2str(t)]];
        end 
    end
    
    % Ploteo la sw total (sumo todos los spikes de todos los trials y
    % calculo la sw)
    k = 15/max(estimulos(e).psth_sw(:,1));
    media_sw = plot(estimulos(e).psth_sw(:,2), estimulos(e).psth_sw(:,1)*k ...
        + dist*(ntrials + 1), '-g', ...
        'Linewidth', 4);
    leyendas = [leyendas, ['SW total']];
    
    % Ploteo la sw proveneniente de promediar cada trial average y luego
    % escaleo (da la misma forma que la sw total)
    media_sw_trials = plot(estimulos(e).psth_sw_trial_avg(:,2), ...
        estimulos(e).psth_sw_trial_avg(:,1)*5 + dist*(ntrials + 2), '-k', ...
        'Linewidth', 4);
    leyendas = [leyendas, ['SW mean all trial']];
    linkaxes(ax, 'x');
    xlim([0 tiempo_file]);
    media_sw.Color(4) = 1;
    
    % Ploteo la sw proveneniente de promediar trial average con pocos spikes
    sw_trials_few_spikes = plot(estimulos(e).sw_few_spikes(:,2), ...
        estimulos(e).sw_few_spikes(:,1) + dist*(ntrials + 4), '-b', ...
        'Linewidth', 5);
    leyendas = [leyendas, ['SW trials few spikes']];
    linkaxes(ax, 'x');
    xlim([0 tiempo_file]);
    media_sw.Color(4) = 0.5;
    
    % Ploteo la sw proveneniente de promediar trial average con pocos spikes
    sw_trials_many_spikes = plot(estimulos(e).sw_many_spikes(:,2), ...
        estimulos(e).sw_many_spikes(:,1) + dist*(ntrials + 6), '-r', ...
        'Linewidth', 5);
    leyendas = [leyendas, ['SW trials many spikes']];
    linkaxes(ax, 'x');
    xlim([0 tiempo_file]);
    media_sw.Color(4) = 0.5;
    
    
    legend(leyendas)

    print_pdf(i, directorio, strcat('_',string(puerto_canal_custom),...
        '_spike_each_trial_', string(round(thr)), 'uV'))
end 




% Ploteo raw data y t0s de estimulos elegidos
figure()
leyendas = {'raw filtered'};

k = 15/300;
plot(raw_filtered*k);
hold on
plot(((1:length(LFP))/downsample_sr)*frequency_parameters.amplifier_sample_rate,...
    15*abs(hilbert(LFP))/max(abs(hilbert(LFP))))
leyendas = [leyendas 'LFP'];
t_window = 5*60; % ventana de 10 minutos (en segundos)
step = 30; % paso de 1 minuto (en segundos)
for e = grilla_psth
    plot(estimulos(e).t0s, zeros(1,20), 'o','Markersize', 20, 'LineWidth',3); 
    
    [sw_t0, sw_t0_t] = sliding_window( ...
        estimulos(e).t0s, frequency_parameters.amplifier_sample_rate,...
        t_window, step, tiempo_file*20*15); 
    
    plot(sw_t0_t * frequency_parameters.amplifier_sample_rate, sw_t0, 'LineWidth',6);

    leyendas = [leyendas estimulos(e).name estimulos(e).name];
end
yline(thr*k, '-r')
leyendas = [leyendas 'thr normalizado'];
title({directorio; [puerto_canal_custom]}, 'Interpreter', 'None', 'Fontsize', 8)
ylim([-15 15]);
legend(leyendas, 'Interpreter', 'None','Location','Southwest');

print_png(i+1, directorio, strcat('_',string(puerto_canal_custom),...
    '_t0s_raw_filtered_', string(round(thr)), 'uV'))





% PLOTEO INT PASA ALTOS/BAJOS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

estimulos_table = struct2table(estimulos);

% Selecciono datos de ese protocolo 
pasa_altos = estimulos_table(estimulos_table.tipo == 'up' , :);
pasa_bajos = estimulos_table(estimulos_table.tipo == 'down' , :);

% Plotear INT PASA-ALTOS
figure();
plot(pasa_altos.frec_corte, pasa_altos.int_norm, '-o')
title({strcat('INT_PASA-ALTOS_', datestr(now, 'yyyy-mm-dd')); ...
string(directorio) ; ...
strcat(string(puerto_canal_custom), "  " , string(thr), "uV",...
"  ntrials:", string(ntrials), "  t_inter_estimulo:", string(tiempo_file)) }, ...
'Interpreter','None')
legend

% Plotear INT PASA-BAJOS
figure();
plot(pasa_bajos.frec_corte, pasa_bajos.int_norm, '-o')
title({strcat('INT_PASA-BAJOS_', datestr(now, 'yyyy-mm-dd')); ...
string(directorio) ; ...
strcat(string(puerto_canal_custom), "  " , string(thr), "uV", "  ntrials:",...
string(ntrials), "  t_inter_estimulo:", string(tiempo_file)) }, 'Interpreter','None')
legend


% Plot auxiliar
figure()
plot(estimulos(10).psth_sw_trial_avg(:,1), 'Linewidth', 3)
hold on
plot(estimulos(12).psth_sw_trial_avg(:,1), 'Linewidth', 3)
title({directorio;[puerto_canal_custom, ' | BOS vs BOS_5-5_DOWN']}, ...
    'Interpreter', 'None')
legend({'BOS','BOS 5.5 UP'})



% Ploteo trial vs spikes por trial
% Para cada estimulo
for i = (1:length(grilla_psth))
    
    e = grilla_psth(i);
    
    % inicializo vector donde guardo numero de spikes por trial y t0s
    spikes_por_t = zeros(ntrials, 3);
    for t = (1:ntrials)
        spikes_por_t(t,1) = length(estimulos(e).spikes_trials(t).spikes);
        spikes_por_t(t,2) = (estimulos(e).t0s(t,1)/60) / ...
        frequency_parameters.amplifier_sample_rate; % t0 en min
        spikes_por_t(t,3) = estimulos(e).spikes_trials(t).trial_id_aux;
    end 
    
    trials_few_spikes = spikes_por_t(:,3) == 1;
    trias_high_spikes = spikes_por_t(:,3) == 0;
    
    figure()
    scatter(spikes_por_t(trials_few_spikes,2), spikes_por_t(trials_few_spikes,1))
    hold on 
    scatter(spikes_por_t(trias_high_spikes,2), spikes_por_t(trias_high_spikes,1))
    yline(estimulos(e).spikes_trials(t).limite);
    xlabel('t0 (min)')
    ylabel('spikes / trial')
    xlim([0 45]);
    title({estimulos(e).name;...
        'spikes/trial a medida que avanzan los trials'},...
        'Interpreter', 'None');
end

