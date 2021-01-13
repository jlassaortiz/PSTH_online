
% Crea el archivo 'songs.mat' con los .wav que encuentra en el directorio
files = dir(horzcat(directorio,'*.wav'));
n = length(files);

% Cargo todos los estimulos y el nombre de los archivos en un struct
estimulos = struct;

for i = 1:n
   estimulos(i).name = files(i).name;
   estimulos(i).dir = horzcat(directorio, files(i).name);
   estimulos(i).song = audioread(estimulos(i).dir);
end

clear files n i