% Script para detectar spikes en la raw  data y plotearlos
% Sirve para analizar la calidad de los datos obtenidos con threshold cutting

% Definimos un umbral y un deadtime de 1 ms
thr = input('Threshold para el threshold cutting (en muV):  ');
deadtime = (1E-3)*frequency_parameters.amplifier_sample_rate; % 1 ms deadtime

% Busca elementos que crucen el umbral
if thr < 0
    spike_times = find(raw_filtered < thr);
else
    spike_times = find(raw_filtered > thr);
end

if isempty(spike_times)
    disp('NO HAY EVENTOS QUE SUPEREN ESE UMBRAL')
    return
end

% Se fija cuales estan separados mas de 1 ms
prueba = diff(spike_times);
spike_times = spike_times(prueba > deadtime);

% Calcula los ISI entre spikes
ISI = diff(spike_times);

% Extrae los spikes shapes de la raw data
for i = 1:length(spike_times)
    if spike_times(i)<(5E-3)*frequency_parameters.amplifier_sample_rate
        continue
    end
    spike_samples(:,i) = raw_filtered((spike_times(i)-(0.8E-3)*frequency_parameters.amplifier_sample_rate):(spike_times(i)+(0.8E-3)*frequency_parameters.amplifier_sample_rate));
end

%%%%%%%%%%
% FIGURA %
%%%%%%%%%%

figure(1);
time_scale = (1/frequency_parameters.amplifier_sample_rate)*(0:1:(length(raw_filtered)-1));

% Raw data y umbral
h(1)=subplot(2,2,[1 2]);
plot(time_scale,raw_filtered);
hold on
plot([0,time_scale(end)],thr*[1 1],'r')
xlim([0 time_scale(end)])
ylim([1.5*min(raw_filtered) 1.5*max(raw_filtered)])
ylabel('V ($\mu$V)','Interpreter','Latex')
xlabel('seg')

title(datestr(now));

% Spikes
h(2)=subplot(2,2,3);
t = (0:1:length(spike_samples(:,1)) -1)/frequency_parameters.amplifier_sample_rate * 1000;
for i=1:(length(spike_times))
    plot(t, spike_samples(:,i));
    hold on
end
plot(t, mean(spike_samples,2),'LineWidth',3,'color','black')
plot(t, mean(spike_samples,2)+std(spike_samples,0,2),'color',1/100*[1 1 1])
plot(t, mean(spike_samples,2)-std(spike_samples,0,2),'color',1/100*[0 1 1])
ylabel('V ($\mu$V)','Interpreter','Latex')
xlabel('mseg')
xlim([0 1.6])

% Distribucion de ISIs
h(3)=subplot(2,2,4);
histogram(ISI/(frequency_parameters.amplifier_sample_rate /1000),0:1:100)
ylabel('bin count')
xlabel('Time (ms)')
xlim([-5 105])

sgtitle({datestr(now, 'yyyy-mm-dd'); ...
    string(directorio) }, 'Interpreter','None')

clear ISI time_scale prueba i deadtime t spike_samples