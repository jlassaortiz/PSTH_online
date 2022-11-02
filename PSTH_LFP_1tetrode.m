% Calcula PSTH y LFP de todos los canales de un tetrodo especificado y la
% seï¿½al promedio (solo funciona con NNx)
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

% Defino tetrodo a levantar con nombre custom name
peine = input('\nDefino tetrodo a levantar (custom name) \n\nPeine (X): ');
tetrodo = input('\nTetrodo (X): ');
puerto_canal_custom = horzcat('P',num2str(peine),'-','T',num2str(tetrodo));

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

% Guardo figuras?
guardar = input('\nGuardo? (1 = SI / 0 = NO) : ');

% Guardo txt?
guardar_txt = input('\nGuardo PSTHsw_1tet y LFP_1tet BOS? (1 = SI / 0 = NO) : ');

% Analizo por bandas?
bandas = input('\nFiltro banda particular? (1 = SI / 0 = NO) : ');
    if bandas == 1
        b_inf = input('\nLimite inferior (Hz): ');
        b_sup = input('\nLimite superior (Hz): ');
    else 
        b_inf = 0;
        b_sup = 400;
    end 

% Cargamos cantidad de trials y tiempo que dura cada uno
ntrials = params.Ntrials
tiempo_file = params.tiempo_entre_estimulos

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
    board_adc_channels, frequency_parameters, directorio, false);

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
    
    % Calculo scores
    estimulos = score_calculator(id_BOS, estimulos, ...
        frequency_parameters, spike_times, ntrials, tiempo_file);

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
    PSTHsw_1tet_BOS_aux(1:l_aux,c) = estimulos_tetrodos(c).canal(id_BOS).psth_sw(:,1);
    LFP_1tet_BOS_aux(:,c) = estimulos_tetrodos(c).canal(id_BOS).LFP_promedio;
end 

PSTHsw_1tet_BOS = mean(PSTHsw_1tet_BOS_aux, 2);
PSTHsw_1tet_BOS(:,2) = t_PSTH;
LFP_1tet_BOS = mean(LFP_1tet_BOS_aux, 2);


if guardar_txt == 1
    
    csvwrite(strcat(directorio, 'PSTHsw_1tet_BOS','_',string(puerto_canal_custom),...
        '_BANDA-', string(b_inf),'-',string(b_sup) ,'Hz_', ...
        string(round(thr)), 'uV', '.txt'), ...
        PSTHsw_1tet_BOS);
    
    csvwrite(strcat(directorio, 'LFP_1tet_BOS','_',string(puerto_canal_custom),...
        '_BANDA-', string(b_inf),'-',string(b_sup) ,'Hz_', ...
        string(round(thr)), 'uV', '.txt'), ...
        LFP_1tet_BOS)
end


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

for i = (1:length(LFP_mean))
    figure()
    subplot(3, 1, 1)
    t_song = (1:length(LFP_mean(i).song))/LFP_mean(i).song_freq;
    plot(t_song, LFP_mean(i).song, 'black');
    xlim([0, 10])
    diferencia = LFP_mean(i).LFP_score_dif;
    texto = ['diferencia: ' , num2str(diferencia)];
    title({[LFP_mean(i).name,' | ',puerto_canal_custom]; texto}, ...
    'Interpreter', 'none');
    
    subplot(3, 1, [2,3])
    plot(LFP_mean(i).LFP_env, 'blue');
    hold on
    yline(LFP_mean(i).LFP_score_aud ,'r', {'FONACION'})
    yline(LFP_mean(i).LFP_score_sil, 'r:', {'NO FONACION'})
end


print_pdf(1, directorio, strcat('_',string(puerto_canal_custom),...
    '_CON_power_LFP-', string(b_inf),'-',string(b_sup),...
    'Hz_',string(round(thr)), 'uV.pdf'))

print_pdf(10, directorio, strcat('_',string(puerto_canal_custom),...
    '_BOS_power_LFP-', string(b_inf),'-',string(b_sup),...
    'Hz_',string(round(thr)), 'uV.pdf'))

print_pdf(13, directorio, strcat('_',string(puerto_canal_custom),...
    '_REV_power_LFP-', string(b_inf),'-',string(b_sup),...
    'Hz_',string(round(thr)), 'uV.pdf'))

close all

% PLOTEO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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


%%%%%%%%%%%%%%%%%%%%

close all
% Ploteo Grilla PSTH
plot_some_raster_LFP_1tetrode(grilla_psth, id_BOS, estimulos_tetrodos, ...
    frequency_parameters, sr_lfp, tiempo_file, ntrials,thr,directorio,spike_times);

suptitle2({datestr(now, 'yyyy-mm-dd'); ...
string(directorio) ; ...
strcat('tetrodo = ',string(puerto_canal_custom),"  |  ", ...
string(thr), "uV", "  |  ", "ntrials:", string(ntrials), "  |  ", ...
"t_inter_estimulo:", string(tiempo_file), "  |  ", ...
'BANDA: ', string(b_inf),'-',string(b_sup),'Hz') })

% Otra estrategia espectrograma
figure()
nwin = 63;
wind = kaiser(nwin,17);
nlap = nwin-10;
nfft = 256;

spectrogram(LFP_1tet_BOS,wind,nlap,nfft,sr_lfp,'yaxis')
ylim([0 50])

figure('DefaultAxesFontSize',18)
fft_plot(LFP_1tet_BOS,sr_lfp);
xlim([10 50])



% PLOTEO INT PASA ALTOS/BAJOS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

estimulos_table = struct2table(estimulos_tetrodos_avg);

% Selecciono datos de ese protocolo 
pasa_altos = estimulos_table(estimulos_table.tipo == 'up' , :);
pasa_bajos = estimulos_table(estimulos_table.tipo == 'down' , :);

% Plotear INT PASA-ALTOS
figure();
plot(pasa_altos.frec_corte, pasa_altos.int_norm_tet, '-o')
title({strcat('INT_PASA-ALTOS_', datestr(now, 'yyyy-mm-dd')); ...
string(directorio) ; ...
strcat(string(puerto_canal_custom), "  " , string(thr), "uV",...
"  ntrials:", string(ntrials), "  t_inter_estimulo:", string(tiempo_file)) }, ...
'Interpreter','None')
legend

% Plotear INT PASA-BAJOS
figure();
plot(pasa_bajos.frec_corte, pasa_bajos.int_norm_tet, '-o')
title({strcat('INT_PASA-BAJOS_', datestr(now, 'yyyy-mm-dd')); ...
string(directorio) ; ...
strcat(string(puerto_canal_custom), "  " , string(thr), "uV", "  ntrials:",...
string(ntrials), "  t_inter_estimulo:", string(tiempo_file)) }, 'Interpreter','None')
legend

% Plotear DIF PASA-ALTOS
figure();
plot(pasa_altos.frec_corte, pasa_altos.LFP_score_dif_tet, '-o')
title({strcat('DIF_PASA-ALTOS_', datestr(now, 'yyyy-mm-dd')); ...
string(directorio) ; ...
strcat(string(puerto_canal_custom), "  " , string(thr), "uV",...
"  ntrials:", string(ntrials), "  t_inter_estimulo:", string(tiempo_file)) }, ...
'Interpreter','None')
legend

% Plotear DIF PASA-BAJOS
figure();
plot(pasa_bajos.frec_corte, pasa_bajos.LFP_score_dif_tet, '-o')
title({strcat('DIF_PASA-BAJOS_', datestr(now, 'yyyy-mm-dd')); ...
string(directorio) ; ...
strcat(string(puerto_canal_custom), "  " , string(thr), "uV", "  ntrials:",...
string(ntrials), "  t_inter_estimulo:", string(tiempo_file)) }, 'Interpreter','None')
legend



% Ploteo diff power en funcion de freq corte %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% Guardo
if guardar == 1

    print_pdf(1, directorio, strcat('_',string(puerto_canal_custom),...
        '_grilla_PSTH-LFP-tetrode_BANDA-', string(b_inf),'-',string(b_sup),...
        'Hz_', string(round(thr)), 'uV.pdf'))
    
    print_pdf(3, directorio, strcat('_',string(puerto_canal_custom),...
        '_FFT-LFP-tetrode_BANDA-', string(b_inf),'-',string(b_sup) ,'Hz_', ...
        string(round(thr)), 'uV.pdf'))

    print_pdf(4, directorio, strcat('_',string(puerto_canal_custom),...
        '_INT_pasa-ALTOS', '.pdf'))
    
    print_pdf(5, directorio, strcat('_',string(puerto_canal_custom),...
        '_INT_pasa-BAJOS', '.pdf'))
    
    print_pdf(6, directorio, strcat('_',string(puerto_canal_custom),...
        '_DIF_pasa-ALTOS', '.pdf'))
    
    print_pdf(7, directorio, strcat('_',string(puerto_canal_custom),...
        '_DIF_pasa-BAJOS', '.pdf'))
   
end

clear estimulos_aux j i  

beep