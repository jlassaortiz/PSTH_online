function plot_espectro(estimulo_name, estimulo_sound, estimulo_freq)

% Plotea sonido y espectrograma del sonido
% 
% Recibe:
% estimulo_name: nombre del estimulo (string)
% estimulo_sound: señal de sonido (vector numerico unidimensional)
% estimulo_freq: frecuencia de sampleo de sonido (int o float)
% 
% Devuelve: 
% figura con señal de sonido y espectrograma

sound_max = max(abs(estimulo_sound));

figure()


s(1) = subplot(2,1,1)
plot((0:1:length(estimulo_sound)-1)/estimulo_freq, estimulo_sound)
ylim([-sound_max sound_max])
ylabel('Sound (a.u.)')
xlim([0 length(estimulo_sound)/estimulo_freq])

set(gca,'FontSize',12*3)
title(estimulo_name, 'Interpreter', 'None', 'FontSize', 32)

s(2) = subplot(2,1,2)
spectrogram(estimulo_sound, gausswin((10E-3)*estimulo_freq, 5), ...
    ceil(0.97*(10E-3)*estimulo_freq), [], estimulo_freq, 'yaxis', 'psd')
colormap(flipud(hot))
colorbar off
ylim([0 16])
xlim([0 length(estimulo_sound)/estimulo_freq])
xlabel('seconds')

linkaxes(s, 'x')

set(gca,'FontSize',12*3)





end