% Script que filtra la raw data de un canal (ideal para protocolos largos)

% Eligo un canal de un puerto en especifico
puerto = input('Puerto a levantar: ','s');
canal = input('Canal de INTAN a filtrar (X-0XX):  ');
puerto_canal = [puerto '-0' num2str(canal,'%.2d')];

% Comparo "puerto_canal" vs la lista de canales en amplifier_channels
aux = strcmp(puerto_canal,{amplifier_channels(:).native_channel_name});
filtrar = find(aux); % indice del canal a levantar (no estan ordenados)
    if(isempty(filtrar))
    disp('NO HUBO COINCIDENCIA!')  
        return
    end  
    
% Levanto la raw data
fid = fopen(horzcat(directorio, 'amplifier.dat'),'r');
fseek(fid, 2*(filtrar -1), 'bof'); 
% 2*(filtrar -1) = off set: Number of bytes to move from origin, specified 
% as an integer. The value of offset can be positive, negative, or zero.
raw = fread(fid, 'int16', 2*(length(amplifier_channels) -1)); 
%2*(length(amplifier_channels) -1) = skip: cantidad de bytes que se saltean antes de llegar al siguiente 
fclose(fid);
raw = raw*0.195; % convierte en microvots

% Define el filtro
filt_spikes = designfilt('highpassiir','DesignMethod','butter','FilterOrder',...
    4,'HalfPowerFrequency',500,'SampleRate',frequency_parameters.amplifier_sample_rate);

% Aplica filtro
raw_filtered = filtfilt(filt_spikes, raw);

% BORRO variables inservibles
clear puerto canal aux fid raw filt_spikes ans filtrar