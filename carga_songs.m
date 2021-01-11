
% Crea el archivo 'songs.mat' con los .wav que encuentra en el directorio
files = dir(fullfile('*.wav'));
n = length(files);
% Esto carga los archivos BOS,CON y REV (orden alfabetico)
% Carga BOS
BOS = audioread(files(1).name);
% Carga CON
CON = audioread(files(2).name);
% Carga REV y frecuencia de sampleo de los estímulos
[REV, sound_fs] = audioread(files(3).name);

% Cargo todos los SYN y el nombre de los archivos en un struct
SYN = struct;

for i = 4:n
   SYN(i-3).name = files(i).name;
   SYN(i-3).song = audioread(files(i).name);
end

clear files n
%save songs.mat