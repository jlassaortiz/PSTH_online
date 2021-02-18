function plot_spikes_shapes(raw_filtered, spike_times, thr, frequency_parameters, directorio)

% Plotea spike shapes, ISI y raw data filtrada
% Necesita:
% 1) la señal neuronal filtrada
% 2) los time stamps (en samples NO en unidades de tiempo) donde estan los spikes 
% 3) el umbral en uV
% 4) el objeto frequency_parameters generado por read_Intan_RHD2000_file.m
% 5) directorio donde estan los archivos de datos analizados

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
title(thr)

% Distribucion de ISIs
h(3)=subplot(2,2,4);
histogram(ISI/(frequency_parameters.amplifier_sample_rate /1000),0:1:100)
ylabel('bin count')
xlabel('Time (ms)')
xlim([-5 105])

sgtitle({datestr(now, 'yyyy-mm-dd'); ...
    string(directorio) }, 'Interpreter','None')

end