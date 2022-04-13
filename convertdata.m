% Convert data 
% Guarda los registros neuronales del Intan en un archivo binario,
% separando por tetrodos

clear all
close all

% Cargo y defino parametros %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Defino directorio
directorio = input('Directorio: ','s');
directorio = horzcat(directorio , '/');

% Leer info INTAN
read_Intan_RHD2000_file(horzcat(directorio, 'info.rhd'));
clear notes spike_triggers supply_voltage_channels aux_input_channels 
clear board_adc_channels

fileinfo = dir([directorio 'amplifier.dat']);
num_channels = length(amplifier_channels); 
largo_raw = fileinfo.bytes /(num_channels * 2);



% Genero struct nombres de tetrodos (P+T) y nombres de canales ind.
NNx = struct();
t = 0; % Cuenta numero de tetrodos totales (16 en total)
for peine = (1:1:4)
    for tetrodo = (1:1:4)
        t = t + 1;
        
        % Guardo nombre del tetrodo (P + T)
        NNx(t).tetrodo_name = horzcat('P',num2str(peine),'-',...
                'T',num2str(tetrodo));
        
        for canal = (1:1:4)
            
            % Genero nombre custom del canal 
            puerto_canal_custom = horzcat(...
                'P',num2str(peine),'-','T',num2str(tetrodo), '-', ...
                num2str(canal)); 
            
            % Traduzco custom_channel_name a native_channel_name
            trad = strcmp(puerto_canal_custom, ...
                {amplifier_channels(:).custom_channel_name});
            
            % Nombre nativo del canal
            puerto_canal = amplifier_channels(trad).native_channel_name;
            
            % Guardo nombre del canal custom y nativo
            NNx(t).chann_list(canal).custom_name = puerto_canal_custom;
            NNx(t).chann_list(canal).native_name = puerto_canal;
        end 
    end
end
clear peine tetrodo canal t trad puerto_canal puerto_canal_custom


% Extraigo y guardo datos agrupando por tetrodos (4 canales por archivo)
% Para cada tetrodo
for t = (1:1:length(NNx))
    
    % Para el primer canal del tetrodo
    c = 1;
    puerto_canal = NNx(t).chann_list(c).native_name;
    
    % Levanto el canal de interes
    raw = read_INTAN_channel(directorio, puerto_canal,...
        amplifier_channels);
    raw = raw'; % convierto en vector fila
    largo_raw = length(raw);
    
    % Inicializo array donde guardo data de todos los canales del tetrodo
    amplifier_1tetrodo = zeros(4,largo_raw);
    
    % Guardo la data del primer canal del tetrodo
    amplifier_1tetrodo(1,:) = raw;
    
    % Para el resto de los canales del tetrodo
    for c = (2:1:4)
        
        puerto_canal = NNx(t).chann_list(c).native_name;
        
        % Levanto el canal de interes y convierto a uV
        raw = read_INTAN_channel(directorio, puerto_canal,...
            amplifier_channels);
        raw = raw';
        
        % Guardo data del canal 
        amplifier_1tetrodo(c,:) = raw;
    end 
    
    % Defino nombre del archivo binario con datos de 1 tetrodo
    filename = ['amplifier_' NNx(t).tetrodo_name '.bin'];
    % Guardo datos de 1 tetrodo
    fid = fopen([directorio  filename],'w');
    fwrite(fid, amplifier_1tetrodo ,'int16'); 
    fclose(fid);
    
    % Chequeo que se guardo bien
    fid = fopen([directorio  filename],'r');
    test = fread(fid,[4 largo_raw], 'int16');
    fclose(fid);

    % PLOTEO CADA TETRODO
    figure()
    
    % Para cada canal
    c = 1;
    for c_aux = [1, 3, 5, 7]
        
        % Ploteo datos que extraje
        a(c_aux) = subplot(8,1,c_aux);
        plot(amplifier_1tetrodo(c,1:300000));
        title([filename '_' num2str(c)], 'Interpreter', 'none')

        % Ploteo datos que guarde
        a(c_aux + 1) = subplot(8,1,c_aux + 1);
        plot(test(c,1:300000));
        title([filename '_' num2str(c) ' - GUARDADA'],'Interpreter','none')
        
        c = c + 1;
    end
    linkaxes(a, 'x');
    ylim(a, [-500 500]);
end
clear c c_aux raw filename t test puerto_canal fid


for p = 1:4
    
    tetrodos_1peine = [];
    
    for t = 1:4
        
        filename = ['amplifierP_' num2str(p) '-T' num2str(t) '.bin'];
     
        % Levanto canales guaradados
        fid = fopen([directorio  filename],'r');
        test = fread(fid,[4 largo_raw], 'int16');
        fclose(fid);
        
        tetrodos_1peine = [tetrodos_1peine; test];
    end
    
    filename2 = ['amplifier_P' num2str(p) '.bin'];
    % Guardo datos de 1 tetrodo
    fid = fopen([directorio  filename2],'w');
    fwrite(fid, tetrodos_1peine ,'int16'); 
    fclose(fid);
end 
        
   
% Me fijo si se guardo bien
for p = 1:4
    
    filename = ['amplifier_P' num2str(p) '.bin'];

    % Levanto canales guaradados
    fid = fopen([directorio  filename],'r');
    test = fread(fid,[16 largo_raw], 'int16');
    fclose(fid);
    
    for t = (0:4:15)
        figure()
        for c = 1:4
            
            c_aux = c + t;
            a(c) = subplot(4,1,c);
            plot(test(c_aux,1:300000));
            title([filename '_T' num2str(round(t/4) + 1) '-' ...
                num2str(c) ' - GUARDADA'], 'Interpreter','none')
        
        end
    end
end