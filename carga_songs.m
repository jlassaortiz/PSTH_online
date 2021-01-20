function estimulos = carga_songs(directorio)
% Crea el archivo 'songs.mat' con los .wav que encuentra en el directorio
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
%clear files n i y Fs