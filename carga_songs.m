function estimulos = carga_songs(directorio)

% Crea struct 'estimulos' con los .wav que encuentra en 'directorio'
%     
%   Entrada:
%   directorio = (string) directorio donde estan los estimulos
%
%   Salida:
%   estimulos = (struct) tiene toda la info de los estimulos
%   estimulos.name = (string) nombre del archivo de audio
%   estimulos.dir = (string) directorio completo del archivo de audio
%   estimulos.song = (vector columna) señal de sonido en formato matriz de matlab
%   estimulos.freq = (double) frecuencia de sampleo del archivo de audio

files = dir(strcat(directorio,'*.wav'));
n = length(files);

% Cargo todos los estimulos y el nombre de los archivos en un struct
estimulos = struct;

for i = 1:n
   estimulos(i).name = files(i).name;
   estimulos(i).dir = strcat(directorio, files(i).name);
   [y, Fs] = audioread(estimulos(i).dir);
   estimulos(i).song = y;
   estimulos(i).freq = Fs;
end

end