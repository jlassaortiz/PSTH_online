% Calcula PSTH y LFP de todos los canales de un tetrodo especificado y la
% se�al promedio (solo funciona con NNx)
% Hace varios graficos

close all
clear all


% Cargo y defino parametros %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Defino directorio
directorio = input('Directorio: ','s');
directorio = horzcat(directorio , '/');

% Leer info INTAN
read_Intan_RHD2000_file(horzcat(directorio, 'info.rhd'));
clear notes spike_triggers supply_voltage_channels aux_input_channels 

% Carga vector con parametros del protocolo experimental
params = readtable(horzcat(directorio,'/parametros_protocolo.txt'), ...
    'Delimiter','\t');

% Carga vector con parametros para el analisis de datos
params_analisis = readtable(...
    horzcat(directorio,'/parametros_analisis.txt'), 'Delimiter','\t');

tetrodos_list = struct();
masTetrodos = 1;
count = 1;
while masTetrodos == 1
    % Defino tetrodos a levantar con nombre custom name
    peine = input('\nDefino tetrodo a levantar (custom name) \n\nPeine (X): ');
    tetrodo = input('\nTetrodo (X): ');
    puerto_canal_custom = horzcat('P',num2str(peine),'-','T',num2str(tetrodo));
    
    tetrodos_list(count).puerto_canal_custom = puerto_canal_custom;
    
    masTetrodos = input('\nAgrego mas tetrodos? (1 = SI / 0 = NO): ');
    
    count = count + 1;
end 

% % Pregunto si ploteo toda la grilla o solo algunos estimulos
% plot_grilla = input('\nPloteo toda la grilla? (1 = SI / 0 = NO): ');
% if plot_grilla == 0
%     grilla_psth = input('\nMatris lineal con numero ID estimulos : ');
% end

% Indico estimulo a graficar (solo uno)
estimulo_ID = input('\nNumero ID del estimulo a graficar: ');

% Pregunto si el umbral se determina manualmente o automaticamente
thr_automatico = input('\Busqueda thr automatica? (1 = SI / 0 = NO) : ');

% Definimos manualmente un umbral para deteccion de spikes (en uV)
if thr_automatico == 0 
    thr = input('\nThreshold para el threshold cutting (en uV):  ');
end

% Guardo figuras?
guardar = input('\Guardo? (1 = SI / 0 = NO) : ');

% % Guardo txt?
% guardar_txt = input('\Guardo PSTHsw_1tet y LFP_1tet BOS? (1 = SI / 0 = NO) : ');

% Cargamos cantidad de trials y tiempo que dura cada uno
ntrials = params.Ntrials
tiempo_file = params.tiempo_entre_estimulos

% Especifico numero de id del BOS y REV
id_BOS = params_analisis.id_bos(1)

% Tamano sliding window
t_window = 0.015; % 15 ms tamaño ventana
step = 0.001; % 1 ms step de ventana

% Tamano del vector sw una vez que recorre todo el canto
size_sw = 0;
tf = t_window;
while tf <= tiempo_file
   size_sw = size_sw + 1;
   tf = tf + step;
end

% Vector de tiempos del sw (en segundos!)
times_sw = zeros(size_sw,1);
times_sw(1,1) = t_window /2;
for i = (2:1:size_sw)
    times_sw(i,1) = times_sw(i -1,1) + step;
end 
clear i

% % Cargo orden de la grilla
% if plot_grilla == 1
%     grilla_psth = str2num(string(params_analisis.grilla_psth(1)))
% end

% Genero songs.mat a partir de las canciones
estimulos = carga_songs(directorio);    

% Cargo id_estimulos 
for i = (1:1:length(estimulos))
    estimulos(i).id = params_analisis.orden(i);
    estimulos(i).frec_corte = params_analisis.freq_corte(i);
    estimulos(i).tipo = categorical(params_analisis.tipo_estimulo(i));
end
clear i 


% Levanto senal neuronal y analizo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Genero struct con nombre de los estimulos y el momento de presentacion
estimulos = find_t0s(estimulos, ntrials, tiempo_file, ...
    board_adc_channels, frequency_parameters, directorio, false);

% Conservo solo el estimulo de interes
estimulos = estimulos(estimulo_ID);

% Inicializo struct donde guardo los struct de cada tetrodo a analizar
estimulos_VARchann = estimulos;
estimulos_VARchann.VARchann = struct();

for t = (1:1:length(tetrodos_list))
    
    estimulos_1chan = estimulos;
    
    % Indico el nombre del tetrodo a analizar
    puerto_canal_custom = tetrodos_list(t).puerto_canal_custom;
    estimulos_VARchann.VARchann(t).puerto_canal_custom = puerto_canal_custom;
        
    % Levanta senal neuronal y la filtra para obtener: LFP de cada canal del
    % tetrodo , LFP promediando todos los canales y SPIKES de cada canal
    % PASO MUY LENTO, MEJORAR !
    [LFP_tetrodo, LFP_canales, spikes_canales]= LFP_1tetrode(directorio,...
        amplifier_channels, frequency_parameters, puerto_canal_custom);

    % Genero struct donde guardo datos de todos los canales
    estimulos_tetrodos = struct();
    
    % Inicializo variables donde guardar psth y lfp de cada canal
    PSTH_avgTetrodo = zeros(size_sw, 4);
    LFP_avgTetrodo = [];

    % Para cada canal del tetrodo
    for c = ( 1:1: size(spikes_canales,2) )

        raw_filtered = spikes_canales(:,c);
        LFP = LFP_canales(:,c);

        % Calculamos el umbral de deteccion de spikes
        if thr_automatico == 1
            thr = find_thr(raw_filtered, estimulos, tiempo_file, ...
            frequency_parameters);
        end

        % Buscamos spike times (en samples, NO unidades tiempo) por umbral
        spike_times =find_spike_times(raw_filtered, thr, frequency_parameters);

        % Genero objeto con RASTERS de todos los estimulos
        estimulos_1chan = generate_raster(spike_times, estimulos_1chan ,...
            tiempo_file, ntrials, frequency_parameters);

        % Calculo LFP promediado por estimulo todos los trials de cada canal
        estimulos_1chan = trialAverage_LFP(LFP, estimulos_1chan, ...
            tiempo_file, ntrials, frequency_parameters);
  
       % Guardo LFP de cada canal 
       LFP_avgTetrodo(:,c) = estimulos_1chan.LFP_promedio;

%         % Calculo scores
%         estimulos = score_calculator(1, estimulos, ...
%             frequency_parameters, spike_times, ntrials);

        % Calculo sliding window de cada canal
        [sw_data, sw_times] = sliding_window(estimulos_1chan.spikes_norm, ...
            frequency_parameters.amplifier_sample_rate, ...
            t_window, step);
        psth_sw = [sw_data, sw_times];
        
        % Guardo PSTH de cada canal
        PSTH_avgTetrodo(1:length(sw_data),c) = sw_data;
        estimulos_1chan.PSTH_1chann = [sw_data, sw_times];
     
        clear i psth_sw

        % Guardo resultados de este canal en una struct con todos los datos
        estimulos_tetrodos(c).canal = estimulos_1chan;
    end
    clear c
    
    % Promedio LFP de todos los canales del tetrodo de este estimulo
    LFP_avgTetrodo = mean(LFP_avgTetrodo, 2);
    LFP_avgTetrodo(:,2) = (0:1:length(LFP_avgTetrodo)-1) / ...
        frequency_parameters.amplifier_sample_rate; % Agrego vector de tiempos (segundos!)
    
    % Promedio PSTH de todos los canales de este estimulo
    PSTH_avgTetrodo = mean(PSTH_avgTetrodo, 2);  
    PSTH_avgTetrodo(:,2) = times_sw;
    
    % Guardo curvas PSTH y LFP promedio tetrodos
    estimulos_VARchann.VARchann(t).LFP_tetrodo = LFP_avgTetrodo;
    estimulos_VARchann.VARchann(t).PSTH_tetrodo = PSTH_avgTetrodo;
    estimulos_VARchann.VARchann(t).canales_tet = estimulos_tetrodos; 

%     % Inicializo y genero vector de tiempos del PSTH
%     t_PSTH = zeros(size_sw,1); % va a estar en SEGUNDOS
%     ti = 0;
%     tf = t_window;
%     for i = (1:1:size_sw)
% 
%         t_PSTH(i) = (ti + tf)/2;
% 
%         ti = ti + step;
%         tf = tf + step;
%     end
%     clear ti tf
% 
%     % Inicializo vectores de PSTH y LFP del BOS promediado por tetrodo
%     PSTHsw_1tet_BOS_aux = zeros(size_sw, length(estimulos_tetrodos));
%     LFP_1tet_BOS_aux = ones(...
%             length(estimulos_tetrodos(1).canal(1).LFP_promedio), ...
%             length(estimulos_tetrodos));
% 
%     for c = 1:4 
%         l_aux = length(estimulos_tetrodos(c).canal(id_BOS).psth_sw(:,1));
%         PSTHsw_1tet_BOS_aux(1:l_aux,c) = estimulos_tetrodos(c).canal(id_BOS).psth_sw(:,1);
%         LFP_1tet_BOS_aux(:,c) = estimulos_tetrodos(c).canal(id_BOS).LFP_promedio;
%     end 
% 
%     PSTHsw_1tet_BOS = mean(PSTHsw_1tet_BOS_aux, 2);
%     PSTHsw_1tet_BOS(:,2) = t_PSTH;
%     LFP_1tet_BOS = mean(LFP_1tet_BOS_aux, 2);
% 
% 
%     if guardar_txt == 1
% 
%         csvwrite([directorio '/PSTHsw_1tet_BOS_' puerto_canal_custom '.txt'], ...
%             PSTHsw_1tet_BOS)
% 
%         csvwrite([directorio '/LFP_1tet_BOS_' puerto_canal_custom '.txt'], ...
%             LFP_1tet_BOS)
%     end
end


%%%%%%%%%%%%%%%%%%%%%% hasta aca codeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee
%%%%%%%%%%%%%%%%%%%%%% hasta aca codeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee

% PLOTEO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Ploteo Grilla PSTH
plot_some_raster_LFP_1tetrode(grilla_psth, id_BOS, estimulos_tetrodos, ...
    frequency_parameters, tiempo_file, ntrials,thr,directorio,spike_times);

suptitle2({datestr(now, 'yyyy-mm-dd'); ...
string(directorio) ; ...
strcat('tetrodo = ',string(puerto_canal_custom),"  |  ", ...
string(thr), "uV", "  |  ", "ntrials:", string(ntrials), "  |  ", ...
"t_inter_estimulo:", string(tiempo_file)) })


% sgtitle({datestr(now, 'yyyy-mm-dd'); ...
% string(directorio) ; ...
% strcat('tetrodo = ',string(puerto_canal_custom),"  |  ", ...
% string(thr), "uV", "  |  ", "ntrials:", string(ntrials), "  |  ", ...
% "t_inter_estimulo:", string(tiempo_file)) }, 'Interpreter','None',...
% 'FontSize',20)

%%%%%%%%%%%%%%%%% desarrolle hasta aca.
% Falta ver si puedo calcular corr de PSTH_sw y LFP de estimulos vs BOS

% % Selecciono datos de ese protocolo 
% estimulos_table = struct2table(estimulos);
% pasa_altos = estimulos_table(estimulos_table.tipo == 'up' , :);
% pasa_bajos = estimulos_table(estimulos_table.tipo == 'down' , :);
% 
% 
% % Plotear INT PASA-ALTOS
% figure();
% plot(pasa_altos.frec_corte, pasa_altos.int_norm, '-o')
% title({strcat('INT_PASA-ALTOS_', datestr(now, 'yyyy-mm-dd')); ...
% string(directorio) ; ...
% strcat(string(puerto_canal), " = ",string(puerto_canal_custom),"  |  ",...
% string(thr), "uV", "  |  ", "ntrials:", string(ntrials), "  |  ", ...
% "t_inter_estimulo:", string(tiempo_file)) }, 'Interpreter','None')
% legend
% set(gca,'FontSize',20)
% 
% 
% % Plotear INT PASA-BAJOS
% figure();
% plot(pasa_bajos.frec_corte, pasa_bajos.int_norm, '-o')
% title({strcat('INT_PASA-BAJOS_', datestr(now, 'yyyy-mm-dd')); ...
% string(directorio) ; ...
% strcat(string(puerto_canal), " = ",string(puerto_canal_custom),"  |  ", ...
% string(thr), "uV", "  |  ", "ntrials:", string(ntrials), "  |  ", ...
% "t_inter_estimulo:", string(tiempo_file)) }, 'Interpreter','None')
% legend
% set(gca,'FontSize',20)
% 
% 
% % Plotear CORR PASA-ALTOS
% figure();
% plot(pasa_altos.frec_corte, pasa_altos.corr, '-o')
% title({strcat('CORR_PASA-ALTOS_', datestr(now, 'yyyy-mm-dd')); ...
% string(directorio) ; ...
% strcat(string(puerto_canal), " = ",string(puerto_canal_custom),"  |  ", ...
% string(thr), "uV", "  |  ", "ntrials:", string(ntrials), "  |  ", ...
% "t_inter_estimulo:", string(tiempo_file)) }, 'Interpreter','None')
% legend
% set(gca,'FontSize',20)
% 
% 
% % Plotear CORR PASA-BAJOS
% figure();
% plot(pasa_bajos.frec_corte, pasa_bajos.corr, '-o')
% title({strcat('CORR_PASA-BAJOS_', datestr(now, 'yyyy-mm-dd')); ...
% string(directorio) ; ...
% strcat(string(puerto_canal), " = ",string(puerto_canal_custom),"  |  ",...
% string(thr), "uV", "  |  ", "ntrials:", string(ntrials), "  |  ", ...
% "t_inter_estimulo:", string(tiempo_file)) }, 'Interpreter','None')
% legend
% set(gca,'FontSize',20)


% Guardo
if guardar == 1
%     print_png(1, directorio, strcat('_',string(puerto_canal_custom),...
%         '_spike-shape_', string(round(thr)), 'uV'))

    print_pdf(1, directorio, strcat('_',string(puerto_canal_custom),...
        '_grilla_PSTH-LFP-tetrode', string(round(thr)), 'uV.pdf'))
    
%     print_pdf(2, directorio, strcat('_',string(puerto_canal_custom),...
%         '_INT_pasa-ALTOS', '.pdf'))
%     print_pdf(3, directorio, strcat('_',string(puerto_canal_custom),...
%         '_INT_pasa-BAJOS', '.pdf'))
%     print_pdf(4, directorio, strcat('_',string(puerto_canal_custom),...
%         '_CORR_pasa-ALTOS', '.pdf'))
%     print_pdf(5, directorio, strcat('_',string(puerto_canal_custom),...
%         '_CORR_pasa-BAJOS', '.pdf'))
end

clear estimulos_aux j i  