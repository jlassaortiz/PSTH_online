% Script para detectar spikes en la raw  data y plotearlos
% Sirve para analizar la calidad de los datos obtenidos con threshold
% cutting
%clear all
%close all
thr =input('Threshold para el threshold cutting (en muV):  ');
deadtime=(1E-3)*sampling_freq; % 1 ms deadtime

% Busca elementos que crucen el umbral
spike_times=find(raw_filtered<thr);
if isempty(spike_times);
    disp('NO HAY EVENTOS QUE SUPEREN ESE UMBRAL')
    return
end
% Se fija cuales estan separados mas de 1 ms
prueba=diff(spike_times);
spike_times=spike_times(prueba>deadtime);

% Calcula los ISI entre spikes
ISI=prueba;

% Extrae los spikes de la raw data
for i=1:length(spike_times)
    if spike_times(i)<(5E-3)*sampling_freq;
        continue
    end
    spike_samples(:,i)=raw_filtered((spike_times(i)-(0.8E-3)*sampling_freq):(spike_times(i)+(0.8E-3)*sampling_freq));
end

%%%%%%%%%%
% FIGURA %
%%%%%%%%%%

% Raw data
figure(1);
time_scale=(1/sampling_freq)*(0:1:(length(raw_filtered)-1));
h(1)=subplot(2,2,[1 2]);
plot(time_scale,raw_filtered);
hold on
plot([0,time_scale(end)],thr*[1 1],'r')
xlim([0 time_scale(end)])
ylim([1.5*min(raw_filtered) 1.5*max(raw_filtered)])
ylabel('V ($\mu$V)','Interpreter','Latex')
xlabel('s')

% Spikes
h(2)=subplot(2,2,3);
for i=1:(length(spike_times))
    plot(spike_samples(:,i));
    hold on
end
plot(mean(spike_samples,2),'LineWidth',3,'color','black')
plot(mean(spike_samples,2)+std(spike_samples,0,2),'color',1/100*[1 1 1])
plot(mean(spike_samples,2)-std(spike_samples,0,2),'color',1/100*[0 1 1])
ylabel('V ($\mu$V)','Interpreter','Latex')
xlabel('Sample #')
xlim([0 2*(0.8E-3)*sampling_freq])


%Distribucion de ISIs
h(3)=subplot(2,2,4);
histogram(ISI/(sampling_freq/1000),0:1:100)
ylabel('bin count')
xlabel('Time (ms)')
xlim([-5 105])

%save(['spikes_' amplifier_channels(chann_number).native_channel_name '.mat'],'spike_times');


print('spikes_ISIS.pdf','-fillpage','-r300','-dpdf')

% hace un clear
clear spike_samples ISI time_scale prueba i deadtime