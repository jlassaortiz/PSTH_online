function t0s_dict = trialAverage_LFP(LFP, t0s_dict, tiempo_file, ntrials, ...
    frequency_parameters)
% "Apila" (promedia) LFP de todos los trials para cada estimulo
%   
%   INPUTS
%   LFP = (matris Nx1) señal de un canal filtrada para conservar el LFP
%
%   t0s_dict = (struct) diccionario con los nombres de todos los estimulos ...
%   y el tiempo en que fue presentado
%
%   tiempo_file = (num) tiempo que dura el estimulo incluyendo el silencio
%
%   ntrials = (num) cantidad de trials que se pasan por estimulos
%
%   frequency_parameters = (num) frecuencia de sampleo de INTAN
% 
%   OUTPUT
%   Agrego LFP promediado de cada estimulo a la struct t0s_dict


% Guardo los spikes separados por estimulo y por trial en un struct

% Para cada estimulo
for i = (1:1:length(t0s_dict))
    
    % Inicializo lista donde guardo los LFP de los trials de este estimulo
    allTrials_LFP = zeros(tiempo_file * ...
        frequency_parameters.amplifier_sample_rate, ntrials);
    
    % Para cada trial
    for j = (1:1: ntrials)
        
        % Defino tiempo inicial y final del trial
        t_inicial = t0s_dict(i).t0s(j);
        t_final = t_inicial + ...
            tiempo_file * frequency_parameters.amplifier_sample_rate -1 ;
        
        % Busco señal LFP correspondiente a ese trial
        trial_LFP = LFP(t_inicial : t_final, 1);
        
        % Agrego LFP de este trial a la lista de LFPs
        allTrials_LFP(:,j) = trial_LFP;
        
    end
    
    % Promedio todos los LFP del trial
    LFP_promedio = mean(allTrials_LFP, 2);
    
    % Guardo los spikes y los trial_id de este estimulo
    t0s_dict(i).LFP_promedio = LFP_promedio; % elimino el primer cero    
end

end

