function raw = read_INTAN_channel(directorio, puerto_canal, amplifier_channels)

% Levanta UN canal del amplificador INTAN
%
%   Se debe especificar: 
%   1) el directorio donde esta el archivo amplifier.dat
%   2) el canal a levantar en formato "PUERTO-NUMERO_CANAL" por ej: A-019
%   ATENCION: el nombre del canal se asume que es el NATIVO
%   3) objeto amplifier_channerls generado por read_INTAN_RHD2000_file.m
%   
%   Devuelve:
%   raw = (vector columna double) es la raw data del canal especificado

% Comparo "puerto_canal" vs la lista de canales en amplifier_channels
aux = strcmp(puerto_canal,{amplifier_channels(:).native_channel_name});
filtrar = find(aux); % indice del canal a levantar (no estan ordenados)
    if(isempty(filtrar))
    disp('NO HUBO COINCIDENCIA!')  
        return
    end  
    
% Levanto la raw data
fid = fopen(horzcat(directorio, 'amplifier.dat'),'r');
fseek(fid, 2*(filtrar -1), 'bof'); 
% 2*(filtrar -1) = off set: Number of bytes to move from origin, specified 
% as an integer. The value of offset can be positive, negative, or zero.
raw = fread(fid, 'int16', 2*(length(amplifier_channels) -1)); 
%2*(length(amplifier_channels) -1) = skip: cantidad de bytes que se saltean antes de llegar al siguiente 
fclose(fid);
raw = raw*0.195; % convierte en microvots

end
