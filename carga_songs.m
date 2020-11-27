
% Crea el archivo 'songs.mat' con los .wav que encuentra en el directorio
files=dir(fullfile('*.wav'));
% Esto carga los archivos BOS,CON y REV (orden alfabetico)
% Carga BOS
BOS=audioread(files(1).name);
% Carga CON
CON=audioread(files(2).name);
% Carga REV
[REV, sound_fs]=audioread(files(3).name);

clear files
%save songs.mat