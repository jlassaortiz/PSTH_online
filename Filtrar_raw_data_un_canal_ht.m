% Script que filtra la raw data de un canal (ideal para protocolos largos)

% Eligo un canal de un puerto en especifico
puerto = input('Puerto a levantar: ','s');
puerto_canal = input('Canal de INTAN a filtrar (X-0XX):  ');
canal = [puerto '-0' num2str(puerto_canal,'%.2d')];

% Leer info INTAN
read_Intan_RHD2000_file('info.rhd');
num_channels=length(amplifier_channels);
sampling_freq=frequency_parameters.amplifier_sample_rate;

% comparo "canal" vs la lista de canales en amplifier_channels

aux=strcmp(canal,{amplifier_channels(:).native_channel_name});
filtrar=find(aux);
    if(isempty(filtrar))
    disp('NO HUBO COINCIDENCIA!')    
        return
    end
    
% me quedo donde hubo una coincidencia    


% Cuando pense que el amplifier.dat estaba estructurado de otra manera
% propuse esta forma de levantar la data

% if(puerto == 'A')
%     filtrar = 2*(puerto_canal + 1) -1 ;
% else
%     filtrar = 2*(puerto_canal +1);
% end

% Levanto la raw data
fileinfo_amplifier = dir('amplifier.dat');
num_samples_raw = fileinfo_amplifier.bytes/(num_channels * 2); 
fid=fopen('amplifier.dat','r'); %fileID = fopen(filename,permission)
fseek(fid,2*(filtrar -1),'bof'); 
% 2*(filtrar -1): Off set (Number of bytes to move from origin, specified 
% as an integer. The value of offset can be positive, negative, or zero.
raw = fread(fid, 'int16', 2*(num_channels -1) ); 
fclose(fid);
raw=raw*0.195;


% Define el filtro
filt_spikes=designfilt('highpassiir','DesignMethod','butter','FilterOrder',...
    4,'HalfPowerFrequency',500,'SampleRate',sampling_freq);

% Aplica filtro
raw_filtered=filtfilt(filt_spikes,raw);

% Guarda
%save(['raw_filtered_' canal '.mat'],'raw_filtered');

% 
% % Plotea pre-view
% t_ms = [1:1:length(raw_filtered)];
% t_ms = t_ms / (sampling_freq*1E-3); % creo vector de tiempo (en ms)
% plot(t_ms, raw_filtered);
% ylabel('mV');
% xlabel('ms');
% ylim([-300 300]);

% BORRO variables inservibles
chann_number=filtrar;
clear raw fid filt_spikes fileinfo_amplifier canal ans aux filtrar spike_triggers supply_voltage_channels notes
clear num_samples_raw t_ms 
 