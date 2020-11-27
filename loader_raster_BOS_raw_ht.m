% ESTE SCRIPT HACE LO SIGUIENTE:
% i)   Carga datos post-proceso con klusta
% ii)  Se le da como input UN 'buen' cluster.
% iii) Arma rasters y calcula PSTHs.
% iv) Necesita que en el directorio haya un .mat con BOS, CON y REV.
% ('songs.mat').
% Todavia esta medio desprolijo, pero esta bastante comentado paso a paso.
% También está documentado en el cuaderno, a partir de pag. 40.
clear all
%close all

% PARAMETROS GLOBALES
ntrials = input('Numero de trials: ');
tiempo_file = input('Tiempo entre estimulos (en s): ');


%tiempo_file=25; %(en segundos, duracion del trial)
read_Intan_RHD2000_file('info.rhd'); % Carga la info Intan
timer=100; %(+/- timer = intervalo de 200 samples == 10 ms)
           % Este es el timer para los PSTH

% CARGA DATOS
% Levanto el timestamp de la raw data (sacado de manual Intan)
fileinfo_time = dir('time.dat');
num_samples_time = fileinfo_time.bytes/4; % int32 = 4 bytes
fid = fopen('time.dat', 'r');
t= fread(fid, num_samples_time, 'int32');
fclose(fid);
t = t/frequency_parameters.amplifier_sample_rate; % sample rate from header file
clear fileinfo_time num_samples_time

% Levanto la data de los canales analogicos (sacado del manual Intan)
num_channels_analog = length(board_adc_channels); % ADC input info from header file
fileinfo_analog = dir('analogin.dat');
num_samples_analog = fileinfo_analog.bytes/(num_channels_analog * 2); % uint16 = 2 bytes
fid = fopen('analogin.dat', 'r');
analog = fread(fid, [num_channels_analog, num_samples_analog], 'uint16');
fclose(fid);
analog = analog * 0.000050354; % convert to volts
analog = analog';
clear fileinfo_analog num_samples_analog num_channels_analog

% Levanto los spike times y cluster numbers
%spike_times=hdf5read([nombre_experimento '.kwik'], '/channel_groups/0/spikes/time_samples');
%spikes=hdf5read([nombre_experimento '.kwik'], '/channel_groups/0/spikes/clusters/main');
%load spikes.mat
clear fid

% IDENTIFICA TIEMPOS EN QUE COMENZO CADA ESTIMULO
% Armo el vector de los t_0 en los que se presento cada estimulo.
[pks,lcs]=findpeaks(analog);
test=pks-circshift(pks,[1 -1]);
found=find(test>0.5);
pksf=pks(found);
t0s=lcs(found); % ESTO ES LO QUE ME IMPORTA (me quedo con # de datos)
clear pks lcs test found pksf;

% CARGA EL VECTOR CON EL ORDEN DE LOS ESTIMULOS
% (1=BOS ; 2=CON ; 3=REV)
estimulos=readtable('estimulos.txt','Delimiter','\t','ReadVariableNames',false);
st_vec=table2array(estimulos(:,1));
%idBOS=find(st_vec==1);
%idCON=find(st_vec==2);
%idREV=find(st_vec==3);
t0BOS=t0s(st_vec==1);
%t0CON=t0s(st_vec==2);
%t0REV=t0s(st_vec==3);

% SEPARA EL CLUSTER A ANALIZAR
% Falta buscar una manera de levantarlos automaticamente a todos los 'good'
% Si no va a ser un embole cada post-procesamiento.
%cluster=spike_times(find(spikes==clu));
for k=1;
cluster=spikes;

% RASTER + HISTOGRAMAS (alinea con t0s obtenidos). Un poco nesteado, pero sale!
% BOS,CON y REV (trials)
for i=1:ntrials
  trialBOS{i}=cluster(cluster(cluster<=t0BOS(i)+tiempo_file*sampling_freq)>(t0BOS(i)))-t0BOS(i);
%    trialCON{i}=cluster(cluster(cluster<=t0CON(i)+tiempo_file*sampling_freq)>(t0CON(i)))-t0CON(i);
%      trialREV{i}=cluster(cluster(cluster<=t0REV(i)+tiempo_file*sampling_freq)>(t0REV(i)))-t0REV(i);
end

% BOS,CON y REV (histogramas);
% Esto calcula los histogramas 'suavizados' a partir del raster
psth.t=0:1:(tiempo_file*sampling_freq-1);
BOS_spikes=cell2mat(trialBOS(:));
%CON_spikes=cell2mat(trialCON(:));
%REV_spikes=cell2mat(trialREV(:));
for i=1:length(psth.t)    
    if i>timer && i<=(length(psth.t)-timer)
        psth.BOS(i)=length(find(BOS_spikes(BOS_spikes>psth.t(i-timer))<=psth.t(i+timer)));
 %       psth.CON(i)=length(find(CON_spikes(CON_spikes>psth.t(i-timer))<=psth.t(i+timer)));
  %      psth.REV(i)=length(find(CON_spikes(CON_spikes>psth.t(i-timer))<=psth.t(i+timer)));
    end
    if i<=timer
        psth.BOS(i)=length(find(BOS_spikes(BOS_spikes>psth.t(1))<=psth.t(i+timer)));
   %     psth.CON(i)=length(find(CON_spikes(CON_spikes>psth.t(1))<=psth.t(i+timer)));
    %    psth.REV(i)=length(find(REV_spikes(REV_spikes>psth.t(1))<=psth.t(i+timer)));
    end
    if i>length(psth.t)-timer
        psth.BOS(i)=length(find(BOS_spikes(BOS_spikes>psth.t(i-timer))<=psth.t(end)));
     %   psth.CON(i)=length(find(CON_spikes(CON_spikes>psth.t(i-timer))<=psth.t(end)));
      %  psth.REV(i)=length(find(REV_spikes(REV_spikes>psth.t(i-timer))<=psth.t(end)));
    end
end

% PSTH discretos (subsampling del PSTH 'suave')
j=1;
for i=timer:(2*timer):length(psth.BOS);
psth_BOS(j)=psth.BOS(i);
%psth_CON(j)=psth.CON(i);
%psth_REV(j)=psth.REV(i);
j=j+1;
end
psth_t=timer:(2*timer):length(psth.BOS);
clear j

% % % % % % % % FIGURA % % % % % % % % %

figure(2)
%%%%%%%
% BOS %
%%%%%%%
h(1)=subplot(10,1,[1 2]); % Senal BOS
plot((1000/(sound_fs))*(0:1:(length(BOS)-1)),BOS,'black')
hold on
line([0 tiempo_file*1000],[0 0],'color',[0 0 0])

h(2)=subplot(10,1,[3,5]); % PSTH BOS
%plot(psth.t/20,psth.BOS,'red') % smoothed psth
bar((1000/sound_fs)*psth_t,psth_BOS,1,'r')

h(3)=subplot(10,1,[7,9]); % RASTER
for i=1:ntrials
    if(isempty(trialBOS{:,i})==1)
        continue
    end
line([trialBOS{:,i} trialBOS{:,i}]/(sound_fs/1000),[i-0.45 i+0.45],'color',[1 0 0])
%line([trialBOS{:,i} trialBOS{:,i}],[i-0.45 i+0.45],'color',[1 0 0])
hold on
end

% %%%%%%%
% % CON %
% %%%%%%%
% h(4)=subplot(10,2,2); % Senal CON
% plot((1000/sound_fs)*(0:1:length(CON)-1),CON,'black')
% hold on
% line([0 tiempo_file*1000],[0 0],'color',[0 0 0])
% 
% h(5)=subplot(10,2,[4,6]); % PSTH CON
% %plot(psth.t/20,psth.CON,'red')
% bar((1000/sound_fs)*psth_t,psth_CON,1,'r')
% 
% h(6)=subplot(10,2,[8,10]); % RASTER CON
% for i=1:ntrials
%     if(isempty(trialCON{:,i})==1)
%         continue
%     end
% line([trialCON{:,i} trialCON{:,i}]/(sound_fs/1000),[i-0.45 i+0.45],'color',[1 0 0])
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
% %plot(psth.t/20,psth.REV,'red')
% bar((1000/sound_fs)*psth_t,psth_REV,1,'r')
% 
% h(9)=subplot(10,2,[18,20]); % RASTER REV
% for i=1:ntrials
%     if(isempty(trialREV{:,i})==1)
%         continue
%     end
% line([trialREV{:,i} trialREV{:,i}]/(sound_fs/1000),[i-0.45 i+0.45],'color',[1 0 0])
% hold on
% end

% PROPIEDADES DEL GRAFICO
linkaxes(h,'x')
xlim([0 tiempo_file*1000])
ylim([0 ntrials+1])
set(h(1),'XTick',[],'YTick',[])
set(h(2),'YLim',[0 max(psth_BOS)])
set(h(3),'YLim',[0 ntrials+1])
set(h,'FontSize',8)
end
% saveas(figure(2),'Protocolo_raster+histo.eps','epsc')