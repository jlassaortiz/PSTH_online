% Calculo scores de varios directorios y ploteo los graficos tipo sabana

% Defino directorio donde esta archivo de parametros
directorio_params = input('Directorio parametros: ','s');
directorio_params = horzcat(directorio_params , '/');

% Carga vector con parametros del analisis de datos
params_info = dir(horzcat(directorio_params, '*parametros*.txt'));
params = readtable(horzcat(directorio_params,params_info.name),'Delimiter','\t','ReadVariableNames',false);

% Cargo valores de puerto-canal
puerto = char(params.Var2(1));
canal = char(params.Var2(2));
puerto_canal = horzcat(puerto, '-0', num2str(canal,'%.2d'))

% Cargamos cantidad de trials y tiempo que dura cada uno
ntrials = str2num(char(params.Var2(3)))
tiempo_file = str2num(char(params.Var2(4)))

% Especifico numero de id del BOS
id_BOS = str2num(char(params.Var2(5)))

% Cargo orden de la grilla
grilla = str2num(string(params.Var2(6)))

% Cargo el nombre de los parametros que varian por fila y columna de la grilla
params(7,:)
ejeX_fila = char(params.Var2(7))

params(8,:)
ejeY_col  = char(params.Var2(8))

% Posicion del primer directorio en el archivo de parametros
d = 9;

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
    directorio = horzcat(char(directorios.Var2(j)), '/');% directorio protocolo

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
    t0s_dictionary = find_t0s(estimulos, ntrials, tiempo_file, board_adc_channels, frequency_parameters, directorio, false);

    % Definimos umbral de deteccion de spikes
    thr = find_thr(raw_filtered, t0s_dictionary, tiempo_file, frequency_parameters);

    % Buscamos spike por threshold cutting
    spike_times = find_spike_times(raw_filtered, thr, frequency_parameters);

    % Genero objeto con raster de todos los estimulos
    rasters = generate_raster(spike_times, t0s_dictionary, tiempo_file, ntrials, frequency_parameters);

    % Calculo scores
    dict_score = score_calculator(id_BOS, estimulos, rasters, frequency_parameters);
    
    % Selecciono scores con los que me quedo y hago matriz para graficar
    [mat_scores, cell_estimulos] = scores_struct2mat(grilla,dict_score);
    
    % Agrego estos valores a la struct score_total que recopila todo
    score_total(j).id  = char(directorios.Var1(j)); % nombre corto protocolo
    score_total(j).dir = char(directorios.Var2(j)); % directorio protocolo
    score_total(j).grilla_scores = mat_scores; % array con valores XYZ para graficar
    score_total(j).grilla_nombre_estimulos = cell_estimulos; % cell con nombre de estimulos para no perderles el restro
    
    % Ploteo y guard
    plot_sabana(mat_scores, directorio, ejeY_col, ejeX_fila);
    
    plot_some_raster([1, 2, 3, 7, 10, 4, 8, 11, 5, 9, 12, 6], id_BOS,  estimulos, ...
        rasters, frequency_parameters, tiempo_file, ntrials, puerto_canal, thr, directorio)
    
    print_pdf(1, directorio, strcat('_sabana_INT_', string(round(thr)),'uV', '.pdf'))
    print_pdf(2, directorio, strcat('_sabana_CORR_', string(round(thr)),'uV', '.pdf'))
    print_pdf(3, directorio, strcat('_grilla_', string(round(thr)),'uV', '.pdf')) 
    
    close all
    
    clear amplifier_channels board_adc_channels frequency_parameters
    
end


% Calculo el promedio de todas las grillas
mat_avg = zeros(numel(grilla), 4);

for i = (1:1:length(score_total))
    mat_avg = score_total(i).grilla_scores + mat_avg;   
end

mat_avg = mat_avg / length(score_total);


% Ploteo
plot_some_sabana(score_total, mat_avg, ejeX_fila, ejeY_col);

clear ans params_info d directorio directorio_aux  puerto canal i j 
clear mat_scores cell_estimulos rasters raw raw_filtered spike_times
clear t0s_dictionary thr dict_score 
