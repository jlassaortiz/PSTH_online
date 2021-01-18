% PARAMETROS GLOBALES
ntrials = input('Numero de trials: ');

tiempo_file = input('Tiempo entre estimulos (en s): ');

% CARGA DATOS
% % Levanto el timestamp de la raw data (sacado de manual Intan y modificado levemente)
% fileinfo_time = dir(horzcat(directorio,'time.dat'));
% num_samples_time = fileinfo_time.bytes/4; % int32 = 4 bytes
% fid = fopen(horzcat(directorio,'time.dat'), 'r');
% t= fread(fid, num_samples_time, 'int32');
% fclose(fid);
% t = t/frequency_parameters.amplifier_sample_rate; % sample rate from header file
% clear fileinfo_time num_samples_time

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

% No tengo idea como pero funciona bien
[pks,lcs] = findpeaks(filtfilt(bpf,analog));
test = diff(pks);
found = find(test > 0.5);
t0s = lcs(found); % ESTO ES LO QUE ME IMPORTA (me quedo con # de datos)
a = find(diff(t0s) < 6*frequency_parameters.board_adc_sample_rate) + 1;
t0s(a) = [];

if (length(t0s) ~= (ntrials * n_estimulos))
    disp('ERROR EN CANTIDAD DE T0s.')
    return
end

clear pks lcs test found pksf a ans bpf;

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
% raster.trials = struct();

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
        
%         % Guardo el numero de cada trial
%         raster(i).trials(j).ntrial = j;
        
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

figure()
title(raster(1).estimulo , 'Interpreter','None')

% PLOTEO

figure()
n = 5 * round(length(estimulos)/2);
m = 2;

for i = (1:1: length(raster))
    
    % sonido
    subplot(n, m , i);
    plot(1000/estimulos(i).freq * (0:1:(length(estimulos(i).song) -1)), estimulos(i).song,'black')
    hold on;
    line([0 tiempo_file*1000],[0 0],'color',[0 0 0]);
    title(estimulos(i).name, 'Interpreter','None')
    
    % psth
    subplot(n, m, [n_estimulos + 2, n_estimulos + 4])
    histogram(raster(i).spikes_norm, 0.015 * frequency_parameters.amplifier_sample_rate);
    
    % raster
    subplot(n, m, [n_estimulo + 6, n_estimulos + 8])
    plot(raster(i).spikes_norm, raster(i).trials_id, '.')
   
    % titulo general
    sgtitle(strcat(string(puerto_canal), " " , string(thr), "uV"))
    
end


% % RASTER + HISTOGRAMAS (alinea con t0s obtenidos). Un poco nesteado, pero sale!
% % BOS,CON y REV (trials)
% for i=1:ntrials
%   trialBOS{i}=cluster(cluster(cluster<=t0BOS(i)+tiempo_file*sampling_freq)>(t0BOS(i)))-t0BOS(i);
%     trialCON{i}=cluster(cluster(cluster<=t0CON(i)+tiempo_file*sampling_freq)>(t0CON(i)))-t0CON(i);
%       trialREV{i}=cluster(cluster(cluster<=t0REV(i)+tiempo_file*sampling_freq)>(t0REV(i)))-t0REV(i);
% end
% 
% % BOS,CON y REV (histogramas);
% % Esto calcula los histogramas 'suavizados' a partir del raster
% psth.t=0:1:(tiempo_file*sampling_freq-1);
% BOS_spikes=cell2mat(trialBOS(:));
% CON_spikes=cell2mat(trialCON(:));
% REV_spikes=cell2mat(trialREV(:));
% % for i=1:length(psth.t)
% %     if i>timer && i<=(length(psth.t)-timer)
% %         psth.BOS(i)=length(find(BOS_spikes(BOS_spikes>psth.t(i-timer))<=psth.t(i+timer)));
% %         psth.CON(i)=length(find(CON_spikes(CON_spikes>psth.t(i-timer))<=psth.t(i+timer)));
% %         psth.REV(i)=length(find(REV_spikes(REV_spikes>psth.t(i-timer))<=psth.t(i+timer)));
% %     end
% %     if i<=timer
% %         psth.BOS(i)=length(find(BOS_spikes(BOS_spikes>psth.t(1))<=psth.t(i+timer)));
% %         psth.CON(i)=length(find(CON_spikes(CON_spikes>psth.t(1))<=psth.t(i+timer)));
% %         psth.REV(i)=length(find(REV_spikes(REV_spikes>psth.t(1))<=psth.t(i+timer)));
% %     end
% %     if i>length(psth.t)-timer
% %         psth.BOS(i)=length(find(BOS_spikes(BOS_spikes>psth.t(i-timer))<=psth.t(end)));
% %         psth.CON(i)=length(find(CON_spikes(CON_spikes>psth.t(i-timer))<=psth.t(end)));
% %         psth.REV(i)=length(find(REV_spikes(REV_spikes>psth.t(i-timer))<=psth.t(end)));
% %     end
% % end
% 
% % PSTH discretos (subsampling del PSTH 'suave')
% % j=1;
% % for i=timer:(2*timer):length(psth.BOS);
% % psth_BOS(j)=psth.BOS(i);
% % psth_CON(j)=psth.CON(i);
% % psth_REV(j)=psth.REV(i);
% % j=j+1;
% % end
% % psth_t=timer:(2*timer):length(psth.BOS);
% % clear j
% 
% % % % % % % % % FIGURA % % % % % % % % %
% figure(2)
% %load songs.mat
% %%%%%%%
% % BOS %
% %%%%%%%
% h(1)=subplot(10,2,1); % Senal BOS
% plot((1000/(sound_fs))*(0:1:(length(BOS)-1)),BOS,'black')
% hold on
% line([0 tiempo_file*1000],[0 0],'color',[0 0 0])
% 
% h(2)=subplot(10,2,[3,5]); % PSTH BOS
% histogram((1000/sampling_freq)*cell2mat(trialBOS(:)),(1000/sampling_freq)*(-1000:(0.015*sampling_freq):(tiempo_file*sampling_freq)),'FaceColor','blue')
% 
% h(3)=subplot(10,2,[7,9]); % RASTER
% for i=1:ntrials
%     if(isempty(trialBOS{:,i})==1)
%         continue
%     end
% line([trialBOS{:,i} trialBOS{:,i}]/(sampling_freq/1000),[i-0.45 i+0.45],'color',[1 0 0])
% hold on
% end
% 
% %%%%%%%
% % CON %
% %%%%%%%
% h(4)=subplot(10,2,2); % Senal CON
% plot((1000/44100)*(0:1:length(CON)-1),CON,'black')
% hold on
% line([0 tiempo_file*1000],[0 0],'color',[0 0 0])
% 
% h(5)=subplot(10,2,[4,6]); % PSTH CON
% histogram((1000/sampling_freq)*cell2mat(trialCON(:)),(1000/sampling_freq)*(-1000:(0.015*sampling_freq):(tiempo_file*sampling_freq)),'FaceColor','blue')
% 
% h(6)=subplot(10,2,[8,10]); % RASTER CON
% for i=1:ntrials
%     if(isempty(trialCON{:,i})==1)
%         continue
%     end
% line([trialCON{:,i} trialCON{:,i}]/(sampling_freq/1000),[i-0.45 i+0.45],'color',[1 0 0])
% end
% 
% %%%%%%%
% % REV %
% %%%%%%%
% h(7)=subplot(10,2,12); % Senal REV
% plot((1000/sound_fs)*(0:1:length(REV)-1),REV,'black')
% hold on
% line([0 tiempo_file*1000],[0 0],'color',[0 0 0])
% 
% h(8)=subplot(10,2,[14,16]); % PSTH REV
% histogram((1000/sampling_freq)*cell2mat(trialREV(:)),(1000/sampling_freq)*(-1000:(0.015*sampling_freq):(tiempo_file*sampling_freq)),'FaceColor','blue')
% 
% h(9)=subplot(10,2,[18,20]); % RASTER REV
% for i=1:ntrials
%     if(isempty(trialREV{:,i})==1)
%         continue
%     end
% line([trialREV{:,i} trialREV{:,i}]/(sampling_freq/1000),[i-0.45 i+0.45],'color',[1 0 0])
% hold on
% end
% 
% % PROPIEDADES DEL GRAFICO
% linkaxes(h,'x')
% linkaxes([h(2),h(5),h(8)],'y')
% xlim([0 tiempo_file*1000])
% ylim([0 ntrials+120])
% set([h(1),h(4),h(7)],'XTick',[],'YTick',[])
% %set([h(2),h(5),h(8)],'YLim',[0 1.5*max(psth_BOS)])
% set([h(3),h(6),h(9)],'YLim',[0 ntrials+1])
% %set([h(2),h(5),h(8)],'YLim',[0 30])
% set(h,'FontSize',8)
% %saveas(figure(2),'Protocolo_raster+histo.eps','epsc')
% 
% clear a ans i timer bpf st_vec t0s
