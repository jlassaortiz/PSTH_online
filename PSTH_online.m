% Script que hace todo de una

close all
clear all

% Defino directorio
directorio = input('Directorio: ','s');
directorio = horzcat(directorio , '/');

% Pregunto si el umbral se determina manualmente o automaticamente
thr_automatico = input('\n¿Busqueda de thr automatica? (1 = SI / 0 = NO) : ');

% Definimos manualmente un umbral para deteccion de spikes (en uV)
if thr_automatico == 0 
    thr = input('\nThreshold para el threshold cutting (en uV):  ');
end

% Plot sabana?
sabana = input('\n¿Ploteo sabanas? (1 = SI / 0 = NO) : ');

if sabana == 1
    % NO hay sabana entonces no hay diag?
    no_diag = input('\n¿Ploteo diag ext? (1 = SI / 0 = NO) : ');
end

% Carga vector con parametros del analisis de datos
params_info = dir(horzcat(directorio, 'parametros.txt'));
params = readtable(horzcat(directorio,params_info.name),'Delimiter','\t',...
    'ReadVariableNames',false);
clear params_info

% Cargo valores de puerto-canal
puerto = char(params.Var2(1));
canal = char(params.Var2(2));
puerto_canal = horzcat(puerto, '-0', num2str(canal,'%.2d'))
clear puerto canal

% Cargamos cantidad de trials y tiempo que dura cada uno
ntrials = str2num(char(params.Var2(3)))
tiempo_file = str2num(char(params.Var2(4)))

% Debo especificar que trial va a conservar para analisis (harcodeo que
% analice todos)
trials = (1:ntrials);
% trials = (1:10)
% trials = (11:20)


% Especifico numero de id del BOS
id_BOS = str2num(char(params.Var2(5)))

% Cargo orden de la grilla
grilla_sabana = str2num(string(params.Var2(6)))
grilla_psth = str2num(string(params.Var2(7)))

% Cargo el nombre de los parametros que varian por fila y columna de la grilla
char(params.Var1(8))
ejeX_fila = char(params.Var2(8))

char(params.Var1(9))
ejeY_col  = char(params.Var2(9))

% Genero songs.mat a partir de las canciones
estimulos = carga_songs(directorio);
duracion_BOS = length(estimulos(id_BOS).song) / estimulos(id_BOS).freq;
BOS_times = (0:1:length(estimulos(1).song)-1) / estimulos(1).freq; 

% Leer info INTAN
read_Intan_RHD2000_file(horzcat(directorio, 'info.rhd'));
clear notes spike_triggers supply_voltage_channels aux_input_channels 

% Levanto el canal de interes
raw = read_INTAN_channel(directorio, puerto_canal, amplifier_channels);

% Define el filtro
filt_spikes = designfilt('highpassiir','DesignMethod','butter','FilterOrder',...
    4,'HalfPowerFrequency',500,'SampleRate',frequency_parameters.amplifier_sample_rate);

% Aplica filtro
raw_filtered = filtfilt(filt_spikes, raw);
clear puerto canal filt_spikes raw

% Genero diccionario con nombre de los estimulos y el momento de presentacion
t0s_dictionary = find_t0s(estimulos, ntrials, tiempo_file, board_adc_channels, ...
    frequency_parameters, directorio, true, trials);

print_png(1, directorio, '_t0s')

close all

ntrials = length(trials)

% Definimos un umbral para threshold cutting de manera automatica (en uV)
if thr_automatico == 1
    thr = find_thr(raw_filtered, t0s_dictionary, tiempo_file, frequency_parameters);
end
clear thr_automatico

% Buscamos spike por threshold cutting
spike_times = find_spike_times(raw_filtered, thr, frequency_parameters);

% Genero objeto con raster de todos los estimulos
rasters = generate_raster(spike_times, t0s_dictionary, tiempo_file, ntrials, ...
    frequency_parameters);

% para cada estimulo acumulo los spikes por motivo
motif_ti = [1.4509, 2.2344, 3.0151];  
motif_dur = 0.61;
for e = (1:1:length(rasters))
    [spikes_norm_motif, motif] = group_by_motif(rasters(e).spikes_norm,...
        motif_ti, motif_dur, frequency_parameters.amplifier_sample_rate, ...
        rasters(e).song, rasters(e).freq);
    
    rasters(e).spikes_norm_motif = spikes_norm_motif;
    rasters(e).motif = motif;
end

% Evaluo desempleño de los distintos estimulos
dict_score = score_calculator(id_BOS, rasters, frequency_parameters, ...
    spike_times, ntrials, tiempo_file);

% Evaluo desempeño de los distintos estímulos luego de acumular los spikes
% por motivo
for e = (1:1:length(rasters))
    [int_norm_motif, corr_motif, sw_motif] = score_calculator_motif(rasters(e).spikes_norm_motif, ...
        rasters(1).spikes_norm_motif, ...
        duracion_BOS, ...
        frequency_parameters.amplifier_sample_rate, ...
        spike_times, ntrials, tiempo_file);
    
    dict_score(e).int_norm_motif = int_norm_motif;
    dict_score(e).corr_motif = corr_motif;
    dict_score(e).sw_motif = sw_motif;
    
    clear int_norm_motif corr_motif sw_motif
end




% Transformo alguno de los resultados en grillas (si quiero graficar sabanas)
if sabana == 1
    [mat_scores, cell_estimulos] = scores_struct2mat(grilla_sabana,dict_score);
end 

% Machete algunos ploteos

% Carga datos filtrados y hace un threshold cutting
plot_spikes_shapes(raw_filtered, spike_times, thr, frequency_parameters, directorio)

% Grafica raster de todos los estimulos
plot_some_raster_v2(grilla_psth, id_BOS, estimulos, dict_score, ...
    frequency_parameters, tiempo_file, ntrials, puerto_canal, thr, ...
    directorio, spike_times)

% Ploteo sabana (x4 plots)
if sabana == 1
    plot_sabana(mat_scores, directorio, ejeY_col, ejeX_fila)
    
    if no_diag == 1
    figure()
    plot(mat_scores(:,1), mat_scores(:,3))
    xticks([1 2 3 4 5])
    xticklabels({'0.5', '1.0', '1.5', '2.0', '2.5'})
    xlabel('Lambda')
    ylabel('INT')
    title('INT vs id lambda')
    ylim([-0.2, 1.2])
    end
end

% Guardo
print_png(1, directorio, strcat('_spike-shape_', string(round(thr)), 'uV'))
print_pdf(2, directorio, strcat('_grilla_', string(round(thr)), 'uV.pdf'))

if sabana == 1
    print_pdf(3, directorio, strcat('_sabana_INT_', string(round(thr)), 'uV.pdf'))
    print_pdf(4, directorio, strcat('_sabana_CORR_', string(round(thr)), 'uV.pdf'))
    print_pdf(5, directorio, strcat('_CORTE_sabana_INT_', string(round(thr)), 'uV.pdf'))
    print_pdf(6, directorio, strcat('_CORTE_sabana_CORR_', string(round(thr)), 'uV.pdf'))
    
    if no_diag == 1
    print_pdf(7, directorio, strcat('_INT_vs_lambda', string(round(thr)), 'uV.pdf'))
    end

end

