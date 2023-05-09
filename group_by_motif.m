function [spikes_norm_motif, motif] = group_by_motif(spikes_norm, ...
    motif_ti, motif_dur, sr_spikes, song, sr_song)
% Acumula los spikes por motivo para un mismo estimulo
%   Necesita los spikes normalizados (en samples)
%   motif_ti es una lista en segundos de los lugares donde empiezan los motif
%   motif_dur es un numero en segundos que indica cuanto duran los motif
%   sr es el sampling rate de los spikes

% inicializo el vector donde guardo el resultado final
spikes_norm_motif = ones(1,1);

% para cada motivo
for i = (1:1:length(motif_ti))
    
    % Separo spikes del motivo i
    spikes_aux = spikes_norm(spikes_norm >  motif_ti(i)*sr_spikes & ...
        spikes_norm < (motif_ti(i) + motif_dur)*sr_spikes);
    
%     Normalizo los spikes de este motivo
    for j = (1:1:length(spikes_aux))
        spikes_aux(j) = uint64(spikes_aux(j) - motif_ti(i)*sr_spikes);
    end
    
    % Agrego spikes normalizados del motivo i a la lista final 
    spikes_norm_motif = vertcat(spikes_norm_motif, spikes_aux);
end

% Guardo los spikes apilados y normalizados del motivo
spikes_norm_motif = spikes_norm_motif(2:end); 

% extraigo del audio el motivo
ti_aux = uint64( motif_ti(1)*sr_song);
tf_aux = uint64((motif_ti(1) + motif_dur) * sr_song);

motif = song( ti_aux : tf_aux );

end

