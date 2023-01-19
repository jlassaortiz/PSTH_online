% Calcula PSTH y LFP de todos los canales de un tetrodo especificado y la
% senal promedio (solo funciona con NNx)
% Hace varios graficos

close all
clear all

% Cargo todos los directorios

dir_aux = input('Directorio: ','s');
dir_aux = horzcat(dir_aux , '/');

% Pregunto si ploteo toda la grilla o solo algunos estimulos
plot_grilla = input('\nPloteo toda la grilla? (1 = SI / 0 = NO): ');
if plot_grilla == 0
    grilla_psth = input('\nMatris lineal con numero ID estimulos : ');
end

% Pregunto si el umbral se determina manualmente o automaticamente
thr_automatico = input('\nBusqueda thr automatica? (1 = SI / 0 = NO) : ');

% Definimos manualmente un umbral para deteccion de spikes (en uV)
if thr_automatico == 0 
    thr = input('\nThreshold para el threshold cutting (en uV):  ');
end

% Guardo txt?
% guardar_txt = input('\nGuardo PSTHsw_1tet y LFP_1tet BOS? (1 = SI / 0 = NO) : ');

% Guardo figuras?
guardar = input('\nGuardo? (1 = SI / 0 = NO) : ');

% Analizo por bandas?
bandas = input('\nFiltro banda particular? (1 = SI / 0 = NO) : ');
    if bandas == 1
        b_inf = input('\nLimite inferior (Hz): ');
        b_sup = input('\nLimite superior (Hz): ');
    else 
        b_inf = 0;
        b_sup = 400;
    end 

% apt : para cada Ave, Protocolo, Tetrodo (con buena señal)
dir_list_apt = dir(horzcat(dir_aux, 'lista_protocolos_tetrodos*.txt'));
list_apt = readtable(horzcat(dir_aux,dir_list_apt.name),'Delimiter',',',...
    'ReadVariableNames',true);

estimulos_apt = struct()

% apt : para cada Ave, Protocolo, Tetrodo (con buena señal)
for apt = (1:height(list_apt))
    
    apt_id = char(list_apt.leyenda(apt))
    
    puerto_canal_custom = char(list_apt.tetrodo(apt))
    
    directorio = char(list_apt.dir_apt(apt))
    
    
    estimulos_apt(apt).id = apt_id;
    estimulos_apt(apt).puerto_canal_custom = puerto_canal_custom;
    estimulos_apt(apt).directorio = directorio;
    
    % Cargo y defino parametros %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % % Defino directorio
    % directorio = input('Directorio: ','s');
    % directorio = horzcat(directorio , '/');

    % Leer info INTAN
    read_Intan_RHD2000_file(horzcat(directorio, 'info.rhd'));
    clear notes spike_triggers supply_voltage_channels aux_input_channels 

    % Carga vector con parametros del protocolo experimental
    params = readtable(horzcat(directorio,'/parametros_protocolo.txt'), ...
        'Delimiter','\t');

    % Carga vector con parametros para el analisis de datos
    params_analisis = readtable(...
        horzcat(directorio,'/parametros_analisis.txt'), 'Delimiter','\t');

    % Cargamos cantidad de trials y tiempo que dura cada uno
    ntrials = params.Ntrials
    tiempo_file = params.tiempo_entre_estimulos
    trials = (1:ntrials) % subset trials de interes, hardcodeado guardo todos

    % Especifico numero de id del BOS y REV
    id_BOS = params_analisis.id_bos(1)

    % Tamano sliding window
    t_window = 0.015; % 15 ms tamaÃ±o ventana
    step = 0.001; % 1 ms step de ventana

    % Cargo orden de la grilla
    if plot_grilla == 1
        grilla_psth = str2num(string(params_analisis.grilla_psth(1)))
    end

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

    % Levanta senal neuronal y la filtra para obtener: LFP de cada canal del
    % tetrodo , LFP promediando todos los canales y SPIKES de cada canal
    [LFP_tetrodo, LFP_canales, spikes_canales, sr_lfp]= LFP_1tetrode(directorio,...
        amplifier_channels, frequency_parameters, puerto_canal_custom, 1000, ...
        true, b_inf, b_sup);

    % Genero struct con nombre de los estimulos y el momento de presentacion
    estimulos = find_t0s(estimulos, ntrials, tiempo_file, ...
        board_adc_channels, frequency_parameters, directorio, false, trials);

    % Genero struct donde guardo datos de todos los canales
    estimulos_tetrodos = struct();

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
        estimulos = generate_raster(spike_times, estimulos , tiempo_file, ... 
            ntrials, frequency_parameters);

        % Calculo LFP promediado por estimulo todos los trials
        estimulos = trialAverage_LFP(LFP, estimulos, tiempo_file, ntrials, ...
              frequency_parameters, sr_lfp);

%         ntrials = 5;
        % Calculo scores
        estimulos = score_calculator(id_BOS, estimulos, ...
            frequency_parameters, spike_times, ntrials, tiempo_file);
%         ntrials = 20;

        % Calculo sliding window para cada estimulo
        for i = (1:length(estimulos))
            [sw_data, sw_times] = sliding_window(estimulos(i).spikes_norm, ...
                frequency_parameters.amplifier_sample_rate, ...
                t_window, step, tiempo_file);
            psth_sw = [sw_data, sw_times];
            estimulos(i).psth_sw = psth_sw;
        end
        clear i psth_sw

        % Guardo resultados de este canal en una struct con todos los datos
        estimulos_tetrodos(c).canal = estimulos;
    end
    clear c


    % Tamano del vector sw una vez que recorre todo el canto
    size_sw = 0;
    tf = t_window;
    while tf <= tiempo_file
       size_sw = size_sw + 1;
       tf = tf + step;
    end

    % Inicializo y genero vector de tiempos del PSTH
    t_PSTH = zeros(size_sw,1); % va a estar en SEGUNDOS
    ti = 0;
    tf = t_window;
    for i = (1:1:size_sw)

        t_PSTH(i) = (ti + tf)/2;

        ti = ti + step;
        tf = tf + step;
    end
    clear ti tf

    % Inicializo vectores de PSTH y LFP promediado por tetrodo
    PSTHsw_1tet_BOS_aux = zeros(size_sw, length(estimulos_tetrodos));
    LFP_1tet_BOS_aux = ones(...
            length(estimulos_tetrodos(1).canal(1).LFP_promedio), ...
            length(estimulos_tetrodos));

    for c = 1:4 
        l_aux = length(estimulos_tetrodos(c).canal(id_BOS).psth_sw(:,1));
        PSTHsw_1tet_BOS_aux(1:l_aux,c) = ...
            estimulos_tetrodos(c).canal(id_BOS).psth_sw(:,1);
        LFP_1tet_BOS_aux(:,c) = ...
            estimulos_tetrodos(c).canal(id_BOS).LFP_promedio;
    end 

    PSTHsw_1tet_BOS = mean(PSTHsw_1tet_BOS_aux, 2);
    PSTHsw_1tet_BOS(:,2) = t_PSTH;
    LFP_1tet_BOS = mean(LFP_1tet_BOS_aux, 2);

    % Cuantifiaciones LFP_30Hz y estímulos aud %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    LFP_mean = struct();
    MUA_mean = struct();

    % Para cada estimulo seleccionado
    for i = (1:length(estimulos_tetrodos(1).canal))

        % Guardo datos de interes de cada estimulo
        LFP_mean(i).name = estimulos_tetrodos(1).canal(i).name;
        LFP_mean(i).dir = estimulos_tetrodos(1).canal(i).dir;
        LFP_mean(i).song = estimulos_tetrodos(1).canal(i).song;
        LFP_mean(i).song_freq = estimulos_tetrodos(1).canal(i).freq;
        MUA_mean(i).name = estimulos_tetrodos(1).canal(i).name;
        MUA_mean(i).dir = estimulos_tetrodos(1).canal(i).dir;
        MUA_mean(i).song = estimulos_tetrodos(1).canal(i).song;  
        MUA_mean(i).song_freq = estimulos_tetrodos(1).canal(i).freq;


        % Inicializo vectores donde guardo señales de cada canal
        LFP_aux = zeros(length(estimulos_tetrodos(1).canal(i).LFP_promedio),4);
        MUA_aux = zeros(length(estimulos_tetrodos(1).canal(i).psth_sw(:,2)),4);

        % Para cada canal de un tetrodo
        for c = (1:4)
            LFP_aux(:,c) = estimulos_tetrodos(c).canal(i).LFP_promedio;
            MUA_aux = estimulos_tetrodos(c).canal(i).psth_sw(:,2);
        end 

        % Promedio señales de los 4 canales del tetrodo
        LFP_aux = mean(LFP_aux, 2);
        MUA_aux = mean(MUA_aux, 2);
        MUA_aux = horzcat(estimulos_tetrodos(1).canal(i).psth_sw(:,1), MUA_aux);

        % Agrego señales promediadas por tetrodo
        LFP_mean(i).LFP_tet = LFP_aux;
        MUA_mean(i).MUA_tet = MUA_aux;

        % Calculo score de LFP con estimulo y post-estimulo
        t_sil = 4.5*sr_lfp;
        h = abs(hilbert(LFP_aux));
        LFP_score_aud = mean(h(1:t_sil,1));
        LFP_score_sil = mean(h(t_sil:t_sil*2,1));

        LFP_mean(i).LFP_score_aud = LFP_score_aud;
        LFP_mean(i).LFP_score_sil = LFP_score_sil;
        LFP_mean(i).LFP_score_dif = LFP_score_aud - LFP_score_sil;
        LFP_mean(i).LFP_env = h;

        clear LFP_aux MUA_aux
    end 
% 
%     for i = (1:length(LFP_mean))
%         figure()
%         subplot(3, 1, 1)
%         t_song = (1:length(LFP_mean(i).song))/LFP_mean(i).song_freq;
%         plot(t_song, LFP_mean(i).song, 'black');
%         xlim([0, 10])
%         diferencia = LFP_mean(i).LFP_score_dif;
%         texto = ['diferencia: ' , num2str(diferencia)];
%         title({[LFP_mean(i).name,' | ',puerto_canal_custom]; texto}, ...
%         'Interpreter', 'none');
% 
%         subplot(3, 1, [2,3])
%         plot(LFP_mean(i).LFP_env, 'blue');
%         hold on
%         yline(LFP_mean(i).LFP_score_aud ,'r', {'FONACION'})
%         yline(LFP_mean(i).LFP_score_sil, 'r:', {'NO FONACION'})
%     end
% 
% 
%     print_pdf(1, directorio, strcat('_',string(puerto_canal_custom),...
%         '_CON_power_LFP-', string(b_inf),'-',string(b_sup),...
%         'Hz_',string(round(thr)), 'uV.pdf'))
% 
%     print_pdf(10, directorio, strcat('_',string(puerto_canal_custom),...
%         '_BOS_power_LFP-', string(b_inf),'-',string(b_sup),...
%         'Hz_',string(round(thr)), 'uV.pdf'))
% 
%     print_pdf(13, directorio, strcat('_',string(puerto_canal_custom),...
%         '_REV_power_LFP-', string(b_inf),'-',string(b_sup),...
%         'Hz_',string(round(thr)), 'uV.pdf'))
% 
%     close all

    estimulos_tetrodos_avg = estimulos;
    estimulos_tetrodos_avg = rmfield(estimulos_tetrodos_avg, ...
        {'t0s', 'spikes_norm','trials_id', 'LFP_promedio', 'int_norm', ...
        'corr', 'psth_sw'});

    % para cada estimulo
    for e = (1:length(estimulos_tetrodos(1).canal()))

        lfp_aux = zeros(length(estimulos_tetrodos(1).canal(1).LFP_promedio),4);
        int_aux = zeros(1,4);
        corr_aux = zeros(1,4);
        psth_aux = zeros(length(estimulos_tetrodos(1).canal(1).psth_sw(:,1)),4);

        % para cada canal
        for c = (1:length(estimulos_tetrodos))
            lfp_aux(:,c) = estimulos_tetrodos(c).canal(e).LFP_promedio;
            int_aux(1,c) = estimulos_tetrodos(c).canal(e).int_norm;
            corr_aux(1,c) = estimulos_tetrodos(c).canal(e).corr;
            psth_aux(:,c) = estimulos_tetrodos(c).canal(e).psth_sw(:,1);
        end 

        % Promedio los 4 canales
        estimulos_tetrodos_avg(e).LFP_promedio_tet = mean(lfp_aux, 2);
        estimulos_tetrodos_avg(e).int_norm_chan = int_aux;
        estimulos_tetrodos_avg(e).int_norm_tet = mean(int_aux, 2);
        estimulos_tetrodos_avg(e).corr_tet = mean(corr_aux, 2);
        psth_aux = mean(psth_aux, 2);


        % Agrego vector de tiempos a PSTH_SW
        t_aux = estimulos_tetrodos(1).canal(1).psth_sw(:,2);
        psth_aux2 = horzcat(psth_aux, t_aux);
        estimulos_tetrodos_avg(e).psth_sw_tet = psth_aux2;

        % Para cada estimulo agrego cuantificaciones de power LFP
        estimulos_tetrodos_avg(e).LFP_score_aud_tet =  LFP_mean(e).LFP_score_aud;
        estimulos_tetrodos_avg(e).LFP_score_sil_tet = LFP_mean(e).LFP_score_sil; 
        estimulos_tetrodos_avg(e).LFP_score_dif_tet = LFP_mean(e).LFP_score_dif;
        estimulos_tetrodos_avg(e).LFP_env_tet = LFP_mean(e).LFP_env;
    end  
    estimulos_apt(apt).apt = estimulos_tetrodos_avg;
end

beep

% Inicializo vector de 5 filas (total de frecuencias de corte) por n (aves,
% protocolos, tetrodos buenos)
int_altos_avg = zeros(5, length(estimulos_apt));
int_bajos_avg = zeros(5, length(estimulos_apt));
diff_altos_avg = zeros(5, length(estimulos_apt));
diff_bajos_avg = zeros(5, length(estimulos_apt));

f_corte = [1.5, 2.5, 3.5, 4.5, 5.5];

leyendas = {};

% para cada Ave Protocolo Tetrodo bueno
for apt = (1:length(estimulos_apt))
    
    estimulos_table = struct2table(estimulos_apt(apt).apt);
    leyendas =  [leyendas , estimulos_apt(apt).id];

    % Selecciono datos de ese protocolo 
    pasa_altos = estimulos_table(estimulos_table.tipo == 'up' , :);
    pasa_bajos = estimulos_table(estimulos_table.tipo == 'down' , :);

    % Plotear INT PASA-ALTOS
    figure(1);
    plot(pasa_altos.frec_corte, pasa_altos.int_norm_tet, '-o')
    hold on
    title({strcat('INT_PASA-ALTOS_', datestr(now, 'yyyy-mm-dd')); ...
    string(dir_aux)}, ...
    'Interpreter','None')
    int_altos_avg(:,apt) = pasa_altos.int_norm_tet; 
    
    % Plotear INT PASA-BAJOS
    figure(2);
    plot(pasa_bajos.frec_corte, pasa_bajos.int_norm_tet, '-o')
    hold on
    title({strcat('INT_PASA-BAJOS_', datestr(now, 'yyyy-mm-dd')); ...
    string(dir_aux)}, ...
    'Interpreter','None')
    int_bajos_avg(:,apt) = pasa_bajos.int_norm_tet;

    % Plotear DIF PASA-ALTOS
    figure(3);
    plot(pasa_altos.frec_corte, pasa_altos.LFP_score_dif_tet, '-o')
    hold on
    title({strcat('DIF_PASA-ALTOS_', datestr(now, 'yyyy-mm-dd')); ...
    string(dir_aux)}, ...
    'Interpreter','None')
    diff_altos_avg(:,apt) = pasa_altos.LFP_score_dif_tet;

    % Plotear DIF PASA-BAJOS
    figure(4);
    plot(pasa_bajos.frec_corte, pasa_bajos.LFP_score_dif_tet, '-o')
    hold on
    title({strcat('DIF_PASA-BAJOS_', datestr(now, 'yyyy-mm-dd')); ...
    string(dir_aux)}, ...
    'Interpreter','None')
    diff_bajos_avg(:,apt) = pasa_bajos.LFP_score_dif_tet;    
    
end

leyendas = [leyendas, 'Mean'];

int_altos_avg = mean(int_altos_avg, 2);
int_bajos_avg = mean(int_bajos_avg, 2);
diff_altos_avg = mean(diff_altos_avg, 2);
diff_bajos_avg = mean(diff_bajos_avg, 2); 

figure(1);
plot(pasa_altos.frec_corte, int_altos_avg, '-o', 'Linewidth', 6)
legend(leyendas, 'Interpreter', 'None')

figure(2);
plot(pasa_altos.frec_corte, int_bajos_avg, '-o', 'Linewidth', 6)
legend(leyendas, 'Interpreter', 'None', 'Location','northwest')

figure(3);
plot(pasa_altos.frec_corte, diff_altos_avg, '-o', 'Linewidth', 6)
legend(leyendas, 'Interpreter', 'None')

figure(4);
plot(pasa_altos.frec_corte, diff_bajos_avg, '-o', 'Linewidth', 6)
legend(leyendas, 'Interpreter', 'None', 'Location','northwest')

% Guardo
print_pdf(1, dir_aux, '_MEAN_INT_pasa_ALTOS.pdf')
print_pdf(2, dir_aux, '_MEAN_INT_pasa_BAJOS.pdf')
print_pdf(3, dir_aux, '_MEAN_DIF_pasa_ALTOS.pdf')
print_pdf(4, dir_aux, '_MEAN_DIF_pasa_BAJOS.pdf')

clear i j k apt apt_id bandas c corr_aux dif_avg diff_altos_avg diff_bajos_avg
clear e estimulos estimulos_table estimulos_tetrodos estimulos_tetrodos_avg
clear guardar i int_altos_avg int_aux int_bajos_avg leyendas LFP LFP_1tet_BOS
clear LFP_1tet_BOS_aux lfp_aux LFP_canales LFP_mean LFP_score_aud LFP_score_sil
clear MUA_mean params params_analisis pasa_altos pasa_bajos plot_grilla psth_aux
clear psth_aux2 PSTHsw_1tet_BOS PSTHsw_1tet_BOS_aux raw_filtered size_sw 
clear spike_times sr_lfp step sw_data sw_times t_aux t_psth t_sil t_window
clear thr thr_automatico amplifier_channels board_adc_channels
clear frequency_parameters l_aux LFP_tetrodo list_apt puerto_canal_custom t_PSTH
clear ans dir_list_apt spikes_canales


save([dir_aux datestr(now, 'yyyy-mm-dd_HH_MM_SS') '_estimulos_apt.mat'], 'estimulos_apt')

beep