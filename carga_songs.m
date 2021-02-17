function estimulos = carga_songs(directorio)

% Crea struct 'estimulos' con los .wav que encuentra en 'directorio'
% Los atributos de la struct son:
% name: nombre del archivo de audio
% dir: directorio completo del archivo de audio
% song: señal de sonido en formato matriz de matlab
% freq: frecuencia de sampleo del archivo de audio

files = dir(horzcat(directorio,'*.wav'));
n = length(files);

% Cargo todos los estimulos y el nombre de los archivos en un struct
estimulos = struct;

for i = 1:n
   estimulos(i).name = files(i).name;
   estimulos(i).dir = horzcat(directorio, files(i).name);
   [y, Fs] = audioread(estimulos(i).dir);
   estimulos(i).song = y;
   estimulos(i).freq = Fs;
end

end