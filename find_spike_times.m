function spike_times = find_spike_times(raw_filtered, thr, frequency_parameters)

% Devuelve spikes times (en samples, NO en unidades de tiempo) que superen el umbral thr y el deadtime
% Necesita como input:
% 1) la señan neuronal filtrada con pasa altos de 500 Hz
% 2) umbral en uV (microvolts)
% 3) objeto frequency_parameters generado por read_Intan_RHD2000_file.m

% Definimos deadtime de 1ms
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
spike_times = spike_times(prueba > deadtime); %(en samples, NO en unidades de tiempo)

end