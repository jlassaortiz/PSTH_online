function t0s_dict = generate_raster(spike_times, t0s_dict, tiempo_file, ...
    ntrials, frequency_parameters)

% Genera un struct con el raster de cada estimulo
% Necesita como input:
% 1) los time stamp de los spikes (en samples, NO en unidades de tiempo)
% 2) el diccionario con los nombres de todos los estimuolos y el tiempo 
% en que fue presentado
% 3) tiempo que dura la presentacion del estimulo incluyendo el silencio
% 4) objeto frequency_parameters generado por read_Intan_RHD2000_file.m
%
% Devuelve un struct donde las entradas son los nombres de los estimulos, 
% el tiempo (en samples) en que ocurren los spikes "normalizados" al comienzo del
% estimulo y para cada spike se guarda el numero de trial en el que ocurrio

% Guardo los spikes separados por estimulo y por trial en un struct

% Para cada estimulo
for i = (1:1:length(t0s_dict))
    
    % Inicializo la lista donde guardo los spikestimes normalizados de cada
    % estimulo y el id del trial
    spikes_norm = ones(1,1);
    trial_id = ones(1,1);
    spike_trial_aux = struct();
    
    % Para cada trial
    for j = (1:1: ntrials)
        
        % Defino tiempo inicial y final del trial
        t_inicial = t0s_dict(i).t0s(j);
        t_final = t_inicial + ...
            tiempo_file * frequency_parameters.amplifier_sample_rate;
        
        % Busco los spikes que ocurrieron durante el trial
        spikes_trial = spike_times(spike_times > t_inicial & ...
            spike_times < t_final);
        
        % Normalizo los spikes times de cada trial
        spikes_trial = spikes_trial - t_inicial;
        
        % Genero vector con la informacion del numero de trial
        trial = ones(length(spikes_trial), 1) * j;
        
        % Apendeo los spiketimes de este trial a la lista de spike times
        % del estimuolo
        spikes_norm = vertcat(spikes_norm, spikes_trial);
        
        % Apendeo el numero de trial a la lista de trial_id del estimulo
        trial_id = vertcat(trial_id, trial);
        
        % Apendeo los spikes normalizados de esa presentacion aparete
        spike_trial_aux(j).spikes = spikes_trial;
        
        t_window = 0.015; % 15 ms
        step = 0.001; % 1 ms
        limite = tiempo_file;
        [sliding_window_data, sliding_window_tiempo] = sliding_window( ...
            spikes_trial, frequency_parameters.amplifier_sample_rate,...
            t_window, step, limite);

        spike_trial_aux(j).sw =[sliding_window_data, sliding_window_tiempo];
        
    end
    
    % Guardo los spikes y los trial_id de este estimulo
    t0s_dict(i).spikes_norm = spikes_norm(2:end); % elimino el primer cero
    t0s_dict(i).trials_id = trial_id(2:end); % elimino el primer cero
    t0s_dict(i).spikes_trials = spike_trial_aux;
    
end

end

