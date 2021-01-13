
% Crea el archivo 'songs.mat' con los .wav que encuentra en el directorio
files = dir(horzcat(directorio,'*.wav'));
n = length(files);

% Cargo todos los estimulos y el nombre de los archivos en un struct
estimulo = struct;

for i = 1:n
   estimulo(i).name = files(i).name;
   estimulo(i).dir = horzcat(directorio, files(i).name);
   estimulo(i).song = audioread(estimulo(i).dir);
end

clear files n i
