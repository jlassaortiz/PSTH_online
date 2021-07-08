function estimulos = find_t0s(estimulos, ntrials, tiempo_file, board_adc_channels, frequency_parameters, directorio, plotear) 

% Genero diccionario con nombre de archivos de auido y tiempos en que se presentaron 
%
%   Entrada:
%   estimulos = (struct) es la salida de carga_songs
%   ntrials = (double) cantidad de trials del protocolo
%   tiempo_file = (double) tiempo entre presentacion de estimulos
%   board_adc_channels = (struct) datos del canal analogico
%   frequency_parameters = (struct) archivo generado por read_INTAN_channel.m con
%   info de frecuencias de sampleo de los datos adquiridos
%   directorio = (str) directorio donde estan los datos
%   plotear = (bool) indica si se plotea canal analogico junto con la posicion 
%   de los t0s
%
%   Salida:
%   t0s_dictionary (struct)
%   t0s_dictionary.id_estimulo = (str) nombre del archivo de audio
%   t0s_dictionary.t0s = (vector columna) tiempos en los que se presenta el
%   estimulo 
%   Los tiempos estan en samples, NO en unidades de tiempo


% Levanto la data de los canales analogicos (sacado del manual Intan y modificado levemente)
num_channels_analog = length(board_adc_channels); % ADC input info from header file
fileinfo_analog = dir(horzcat(directorio,'analogin.dat'));
num_samples_analog = fileinfo_analog.bytes/(num_channels_analog * 2); % uint16 = 2 bytes
fid = fopen(horzcat(directorio,'analogin.dat'), 'r');
analog = fread(fid, [num_channels_analog, num_samples_analog], 'uint16');
fclose(fid);
analog = analog * 0.000050354; % convert to volts
analog = analog';

% IDENTIFICA TIEMPOS EN QUE COMENZO CADA ESTIMULO
% Armo el vector de los t_0 en los que se presento cada estimulo.
% Pre-acondiciono la senal analogica
bpf = designfilt('bandpassiir','designmethod','butter','halfpowerfrequency1'...
    ,1000,'halfpowerfrequency2',3500,'filterorder',4,'samplerate',frequency_parameters.board_adc_sample_rate);

analog_filt = filtfilt(bpf,analog);

clear analog

% Identifico con un umbral los timestamps de analog_filt donde ocurren los estimulos
umbral = 0.05;
locs_estimulos = find(analog_filt > umbral);

% Inicializo vector donde voy a guardar t0s
t0s = zeros(ntrials * length(estimulos), 1);

% Para cada presentacion de estimulo
for i = (1:1: ntrials * length(estimulos))
    
    % Identifico el primer evento que supera el umbral
    if i == 1
        
        t0s(1,1) = locs_estimulos(1,1);
        
    else   
        % Uso el t0s anterior como referencia para buscar el siguiente
        condicion = locs_estimulos > (t0s(i-1,1) + tiempo_file * frequency_parameters.board_adc_sample_rate) - 1.0 * frequency_parameters.board_adc_sample_rate ;
        t_umbral = locs_estimulos(logical(condicion));

        % Conservo solo el primer elemento que supera el umbral
        t0s(i,1) = t_umbral(1,1);  
    end  
end

clear condicion t_umbral locs_estimulos


% Plotea senal analogica filtrada junto a posicion de los t0s
if plotear
    
    figure()
    t_analog = (0:1:length(analog_filt)-1) / frequency_parameters.board_adc_sample_rate;
    plot(t_analog, analog_filt)
    hold on
    plot(t0s / frequency_parameters.board_adc_sample_rate, ones(length(t0s),1)* umbral, 'or')

    title('analog filtrada');
end

clear pks lcs test found a ans bpf umbral;


% CARGA EL VECTOR CON EL ORDEN DE LOS ESTIMULOS
estimulos_log_info = dir(horzcat(directorio, '*estimulos*.txt'));
estimulos_log = readtable(horzcat(directorio, estimulos_log_info.name),'Delimiter','\t','ReadVariableNames',false);

clear estimulos_log_info

% struct donde guardo el nombre del estimulo y los t0s
t0s_dictionary = struct;

% Recorro cada uno de los nombres de los estimulos presentados
for i = (1:1:length(estimulos))
    % Tomo el nombre de uno de los estimulos
    estimulo = string(estimulos(i).name);
    
    % Inicializo el vector orden
    orden = zeros(height(estimulos_log), 1);
    
    % Recorro cada nombre de la lista de estimulos_log
    % y chequeo que elementos coinciden con 'estimulo'
    for j = (1:1: height(estimulos_log))
        orden(j,1) = string(estimulos_log.Var2(j)) == estimulo;
    end
    
  % Guardo el nombre del estimulo y los t0s en que se presento
  estimulos(i).t0s = t0s(logical(orden));
    
end 

end

