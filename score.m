% Calcula scores

% Defino directorio donde esta archivo de parametros
directorio_params = input('Directorio parametros: ','s');
directorio_params = horzcat(directorio_params , '/');

% Carga vector con parametros del analisis de datos
params_info = dir(horzcat(directorio_params, '*parametros*.txt'));
params = readtable(horzcat(directorio_params,params_info.name),'Delimiter','\t','ReadVariableNames',false);
clear params_info

% Cargo valores de puerto-canal
puerto = char(params.Var2(1));
canal = char(params.Var2(2));
puerto_canal = horzcat(puerto, '-0', num2str(canal,'%.2d'))
clear puerto canal

% Cargamos cantidad de trials y tiempo que dura cada uno
ntrials = str2num(char(params.Var2(3)))
tiempo_file = str2num(char(params.Var2(4)))

% Especifico numero de id del BOS
id_BOS = str2num(char(params.Var2(5)))

% Cargo orden de la grilla
grilla = str2num(string(params.Var2(6)))

% Cargo "estimulos" usando el primer directorio de protocolos de la lista 
% de parametros
% TODOS LOS PROTOCOLOS DEBEN TENER LOS MIMOS ESTIMULOS
directorio_aux = horzcat(char(params.Var2(7)), '/');
estimulos = carga_songs(directorio_aux);

% Genero diccionario donde se va a guardar el score de todos
score_total = struct;

% Para cada directorio (protocolo)
for j = (1:1:length(params.Var2(7:end)))
    
    % Primero guardo los directorios y el nombre corto de los protocolos
    score_total(j).id = char(params.Var1(j+6));
    score_total(j).dir = char(params.Var2(j+6));

    % Defino el directorio del protocolo
    directorio = horzcat(score_total(j).dir, '/');

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
    clear puerto canal filt_spikes

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
    
    % Inicializo la matriz donde van a ir los valores
    % Esta pensado para analizar "SABANAS" es decir variaciones de DOS
    % parametros del SYN y se le calcular DOS scores 
    mat_scores = zeros(numel(grilla), 4);
    
    % Incializo la fila
    fila = 1;
    
    % Recorro cada elemento de la grilla
    for x = (1:1:length(grilla))
        for y = (1:1:length(grilla))
            mat_scores(fila, 1) = x;
            mat_scores(fila, 2) = y;
            
            id_estimulo = grilla(x,y);
            mat_scores(fila, 3) = dict_score(id_estimulo).int_norm;
            mat_scores(fila, 4) = dict_score(id_estimulo).corr;
            
            fila = fila + 1;
        end
    end
    
    score_total(j).grilla_scores = mat_scores;
end

% Calculo el promedio de todas las grillas
mat_avg = zeros(numel(grilla), 4);
for i = (1:1:length(score_total))
    mat_avg = score_total(i).grilla_scores + mat_avg;   
end

mat_avg = mat_avg / length(score_total);


clear x y j id_estimulo id_aux dict_score directorio fila mat_scores directorio_aux
clear rasters raw raw_filtered t0s_dictionary spike_times thr 

% Vectores auxiliares
X  = [];
Y  = [];
Z1 = [];
Z2 = [];

for i = (1:1:length(score_total))
    
X = vertcat(X, score_total(i).grilla_scores(:,1));
Y = vertcat(Y, score_total(i).grilla_scores(:,2));
Z1 = vertcat(Z1, score_total(i).grilla_scores(:,3));
Z2 = vertcat(Z2, score_total(i).grilla_scores(:,4));

end 
clear i

[xq,yq] = meshgrid(1:0.1:3);
z = griddata(X,Y,Z1,xq,yq,'natural');
figure()
plot3(X,Y,Z1,'mo')
hold on
scatter3(mat_avg(:,1), mat_avg(:,2), mat_avg(:,3), 100, 'ro', 'filled')
mesh(xq,yq,z)
ylabel('C')
xlabel('L-traquea')
title('integral')


[xq,yq] = meshgrid(1:0.1:3);
z = griddata(X,Y,Z2,xq,yq,'natural');
figure()
plot3(X,Y,Z2,'mo')
hold on
scatter3(mat_avg(:,1), mat_avg(:,2), mat_avg(:,4), 100, 'ro', 'filled')
mesh(xq,yq,z)
ylabel('C')
xlabel('L-traquea')
title('correlacion')

clear X Y Z1 Z2 z xq yq z2
