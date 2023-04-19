% Calculo scores de varios directorios y ploteo los graficos tipo sabana

close all % Como ploteo y guardo, mejor asegurarme que no haya plots viejos
clear all

% Defino directorio donde esta archivo de parametros
directorio_params = input('Directorio parametros: ','s');
directorio_params = horzcat(directorio_params , '/');

% Guardo graficos sabanas de cada protocolo individual?
guardar_graf_protocolos_ind = input('Guardo graf sabanas de cada protocolo individual? (1 = SI / 0 = NO) : ');

% Carga vector con parametros del analisis de datos
params_info = dir(horzcat(directorio_params, 'parametros.txt'));
params = readtable(horzcat(directorio_params,params_info.name),'Delimiter','\t','ReadVariableNames',false);

% Posicion del primer directorio en el archivo de parametros
d = 10;

% Cargo "estimulos" usando el primer directorio de protocolos de la lista 
% de parametros
% TODOS LOS PROTOCOLOS DEBEN TENER LOS MIMOS ESTIMULOS
directorio_aux = horzcat(char(params.Var2(d)), '/');
estimulos = carga_songs(directorio_aux);

directorios = params(d:end, :);

% Genero diccionario donde se va a guardar el score de todos
score_total = struct;

% Para cada directorio (protocolo)
for j = (1:1:height(directorios))

    % Defino el directorio del protocolo
    directorio = horzcat(char(directorios.Var2(j)), '/') % directorio protocolo
    
    % Carga vector con parametros del analisis de datos
    params_info = dir(horzcat(directorio, 'parametros.txt'));
    params = readtable(horzcat(directorio,params_info.name),'Delimiter','\t','ReadVariableNames',false);
    clear params_info

    % Cargo valores de puerto-canal
    puerto = char(params.Var2(1));
    canal = char(params.Var2(2));
    puerto_canal = horzcat(puerto, '-0', num2str(canal,'%.2d'))
    clear puerto canal

    % Cargamos cantidad de trials y tiempo que dura cada uno
    ntrials = str2num(char(params.Var2(3)))
    tiempo_file = str2num(char(params.Var2(4)))
    
    trials = (1:ntrials);

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
    clear filt_spikes

    % Genero diccionario con nombre de los estimulos y el momento de presentacion
    t0s_dictionary = find_t0s(estimulos, ntrials, tiempo_file, board_adc_channels, ...
        frequency_parameters, directorio, false, trials);
    
    % Definimos umbral de deteccion de spikes
    thr = find_thr(raw_filtered, t0s_dictionary, tiempo_file, frequency_parameters);

    % Buscamos spike por threshold cutting
    spike_times = find_spike_times(raw_filtered, thr, frequency_parameters);

    % Genero objeto con raster de todos los estimulos
    rasters = generate_raster(spike_times, t0s_dictionary, tiempo_file, ntrials, frequency_parameters);

    % Evaluo desempleño de los distintos estimulos
    dict_score = score_calculator(id_BOS, rasters, frequency_parameters, ...
    spike_times, ntrials, tiempo_file);
    
    % Selecciono scores con los que me quedo y hago matriz para graficar
    [mat_scores, cell_estimulos] = scores_struct2mat(grilla_sabana,dict_score);
    
    % Agrego estos valores a la struct score_total que recopila todo
    score_total(j).id  = char(directorios.Var1(j)); % nombre corto protocolo
    score_total(j).dir = char(directorios.Var2(j)); % directorio protocolo
    score_total(j).grilla_scores = mat_scores; % array con valores XYZ para graficar
    score_total(j).grilla_nombre_estimulos = cell_estimulos; % cell con nombre de estimulos para no perderles el restro
   
    % Sabana x4: integral y correlacion
    plot_sabana(mat_scores, directorio, ejeY_col, ejeX_fila); % Hace 4 plots

    % Grafica raster de todos los estimulos
    plot_some_raster(grilla_psth, id_BOS, estimulos, rasters, ...
    frequency_parameters, tiempo_file, ntrials, puerto_canal, thr, ...
    directorio, spike_times)
    
    % Guardo todo en el directorio del protocolo
    if guardar_graf_protocolos_ind == 1
        print_pdf(1, directorio, strcat('_sabana_INT_', string(round(thr)),'uV', '.pdf'))
        print_pdf(2, directorio, strcat('_sabana_CORR_', string(round(thr)),'uV', '.pdf'))
        print_pdf(3, directorio, strcat('_CORTE_sabana_INT_', string(round(thr)),'uV', '.pdf'))
        print_pdf(4, directorio, strcat('_CORTE_sabana_CORR_', string(round(thr)),'uV', '.pdf'))
        print_pdf(5, directorio, strcat('_grilla_', string(round(thr)),'uV', '.pdf')) 
    end
    
    close all
    clear amplifier_channels board_adc_channels frequency_parameters 
end

% Calculo el promedio de todas las grillas
mat_avg = zeros(numel(grilla_sabana), 4);

for i = (1:1:length(score_total))
    mat_avg = score_total(i).grilla_scores + mat_avg;   
end

mat_avg = mat_avg / length(score_total);


% Ploteo
plot_some_sabana(score_total, mat_avg, ejeX_fila, ejeY_col);

% figure(1)
% xticks([1 2 3 4 5])
% xticklabels({'0.5', '1.0', '1.5', '2.0', '2.5'})
% xlabel('Lambda')
% ylabel('INT')
% title('INT vs id lambda')
% view(0,0)

plot_sabana(mat_avg, directorio_params, ejeY_col, ejeX_fila);

% Guardo
print_pdf(1, directorio_params, strcat('_sabana_INT_', 'uV.pdf'))
print_pdf(2, directorio_params, strcat('_sabana_CORR_', 'uV.pdf'))
print_pdf(5, directorio_params, strcat('_CORTE_sabana_INT_', 'uV.pdf'))
print_pdf(6, directorio_params, strcat('_CORTE_sabana_CORR_', 'uV.pdf'))

clear ans params_info d directorio directorio_aux  puerto canal i j 
clear mat_scores cell_estimulos rasters raw raw_filtered spike_times
clear t0s_dictionary thr dict_score 
