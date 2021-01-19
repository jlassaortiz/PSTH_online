% PARAMETROS GLOBALES
ntrials = input('Numero de trials: ');
tiempo_file = input('Tiempo entre estimulos (en s): ');

% Levanto la data de los canales analogicos (sacado del manual Intan y modificado levemente)
num_channels_analog = length(board_adc_channels); % ADC input info from header file
fileinfo_analog = dir(horzcat(directorio,'analogin.dat'));
num_samples_analog = fileinfo_analog.bytes/(num_channels_analog * 2); % uint16 = 2 bytes
fid = fopen(horzcat(directorio,'analogin.dat'), 'r');
analog = fread(fid, [num_channels_analog, num_samples_analog], 'uint16');
fclose(fid);
analog = analog * 0.000050354; % convert to volts
analog = analog';
clear fileinfo_analog num_samples_analog num_channels_analog fid

% IDENTIFICA TIEMPOS EN QUE COMENZO CADA ESTIMULO
% Armo el vector de los t_0 en los que se presento cada estimulo.
% Pre-acondiciono la senal analogica
bpf = designfilt('bandpassiir','designmethod','butter','halfpowerfrequency1'...
    ,1000,'halfpowerfrequency2',3500,'filterorder',4,'samplerate',frequency_parameters.board_adc_sample_rate);

analog_filt = filtfilt(bpf,analog);

% No tengo idea como pero funciona bien
[pks,lcs] = findpeaks(analog_filt);
test = diff(pks);
found = find(test > 0.08);
t0s = lcs(found); % ESTO ES LO QUE ME IMPORTA (me quedo con # de datos)
a = find(diff(t0s) < 6*frequency_parameters.board_adc_sample_rate) + 1;
t0s(a) = [];

% Identifico con un umbral los timestamps de analog_filt donde ocurren los
% estimulos
t0s_beta = zeros(ntrials * length(estimulos), 1);
t0s_picos = find(analog_filt > 0.05);

for i = (0:1: ntrials * length(estimulos) -1)
    
    % Identifico el primer evento que supera el umbral y lo uso de guía
    % Luego me quedo con el primer timestamp que cumple con la condicion:
    % es mayor a = (numero_estimulo * tiempo_file) * 0.9
    condicion = t0s_picos > (t0s_beta(1,1) + i * tiempo_file * frequency_parameters.amplifier_sample_rate) - 0.5 * frequency_parameters.amplifier_sample_rate ;
    t_umbral = t0s_picos(logical(condicion));
   
    % Conservo solo el primer elemento que supera el umbral
    t0s_beta(i+1,1) = t_umbral(1,1);
    
end

% t0s_beta = findchangepts(analog_filt, ...
%     'MaxNumChanges', length(estimulos) * ntrials, ...
%     'MinDistance', 0.9 * tiempo_file * frequency_parameters.board_adc_sample_rate);

% Verifica que cantidad de t0s sea la correcta
if (length(t0s) ~= (ntrials * length(estimulos)))
    
    figure()
    t_analog = (0:1:length(analog)-1) / frequency_parameters.amplifier_sample_rate;
    plot(t_analog, analog_filt)
    hold on
    plot(t0s / frequency_parameters.amplifier_sample_rate, ones(length(t0s),1)*0.1, 'or')
    hold on
    plot(t0s_beta / frequency_parameters.amplifier_sample_rate, ones(length(t0s_beta),1)*0.2, 'or')

    title('analog filtrada');
    
    error('ERROR EN CANTIDAD DE T0s.')
    
end

clear pks lcs test found a ans bpf;

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
  t0s_dictionary(i).id_estimulo = estimulo;
  t0s_dictionary(i).t0s = t0s(logical(orden));
    
end 

clear i j orden estimulo t0s

% Guardo los spikes separados por estimulo y por trial en un struct
raster = struct();

% Para cada estimulo
for i = (1:1:length(t0s_dictionary))
    
    % Guardo el nombre de cada estimulo
    estimulo = string(t0s_dictionary(i).id_estimulo);
    raster(i).estimulo = estimulo;
    
    % Inicializo la lista donde guardo los spikestimes normalizados de cada
    % estimulo y el id del trial
    spikes_norm = ones(1,1);
    trial_id = ones(1,1);
    
    % Para cada trial
    for j = (1:1: ntrials)
        
        % Defino tiempo inicial y final del trial
        t_inicial = t0s_dictionary(i).t0s(j);
        t_final = t_inicial + tiempo_file * frequency_parameters.amplifier_sample_rate;
        
        % Busco los spikes que ocurrieron durante el trial
        spikes_trial = spike_times(spike_times > t_inicial & spike_times < t_final);
        
        % Normalizo los spikes times de cada trial
        spikes_trial = spikes_trial - t_inicial;
        
        % Genero vector con la informacion del numero de trial
        trial = ones(length(spikes_trial), 1) * j;
        
        % Apendeo los spiketimes de este trial a la lista de spike times
        % del estimuolo
        spikes_norm = vertcat(spikes_norm, spikes_trial);
        
        % Apendeo el numero de trial a la lista de trial_id del estimulo
        trial_id = vertcat(trial_id, trial);
        
        % Guardo los spikes y los trial_id de este estimulo
        % raster(i).trials(j).spikes = spikes_trial;
        
    end
    
    % Guardo los spikes y los trial_id de este estimulo
    raster(i).spikes_norm = spikes_norm(2:end); % elimino el primer cero
    raster(i).trials_id = trial_id(2:end); % elimino el primer cero
    
end

clear estimulo t_inicial t_final spikes_trial trial trial_id spikes_norm i j


% PLOTEO

figure()
n = 5 * round(length(estimulos)/2);
m = 2;
j = 0;

for i = (1:1: length(raster))
    
    if mod(i, 2) == 1
        p = (i - 1)/2 * 10 + 1;
    else
        p = ((i / 2) - 1) * 10 + 2;
    end
        
    % sonido
    j = j + 1;
    h(j) = subplot(n, m , p);
    plot(1000/estimulos(i).freq * (0:1:(length(estimulos(i).song) -1)), estimulos(i).song,'black')
    hold on;
    line([0 tiempo_file*1000],[0 0],'color',[0 0 0]);
    xlim([0 tiempo_file * 1000])
    title(estimulos(i).name, 'Interpreter','None')
    
    % psth
    j = j + 1;
    h(j) = subplot(n, m, [p + 2, p + 4]);
    histogram(raster(i).spikes_norm * 1000/frequency_parameters.amplifier_sample_rate , ...
        (1000/frequency_parameters.amplifier_sample_rate) * (-1000:(0.015*frequency_parameters.amplifier_sample_rate):(tiempo_file*frequency_parameters.amplifier_sample_rate)) );
    ylim([0 ntrials + 10]);
    xlim([0 tiempo_file * 1000]);
    
    % raster
    j = j + 1;
    h(j) = subplot(n, m, [p + 6, p + 8]);
    plot((1000/frequency_parameters.amplifier_sample_rate) * raster(i).spikes_norm, raster(i).trials_id, '.')  
    xlim([0 tiempo_file * 1000])
    ylim([0 ntrials + 1])
end

% Linkeo eje x (no se pueden hacer varios links independientes juntos)
linkaxes(h, 'x');

% Titulo general
sgtitle(strcat(string(puerto_canal), " " , string(thr), "uV"))

clear i j m n p h