function [t, song_env] = song_envelope(song)
% song_envelope() Busca envolvente con estrategia de Boari 2015 (automatic)

% Calculamos la transformada de Hilbert y nos quedamos con la parte
% absoluta (envolvente) y normalizamos
hilb = abs(hilbert(song));
hilb = hilb / max(hilb);

% Integramos sistema dinámico lineal unidimensional para suavizar.
    % tiempos de integracion: uno por cada punto en hilb
t_spam = linspace(1,length(hilb),length(hilb));
    % integramos y normalizamos
[t, y] = ode45(@(t,y) takens(t,y, hilb), t_spam, 0);
y = y/max(y);

% Suavizamos pasando un filto Savitzky - Golay y normalizamos
song_env = sgolayfilt(y, 4, 513);
song_env = song_env/max(song_env);
end

