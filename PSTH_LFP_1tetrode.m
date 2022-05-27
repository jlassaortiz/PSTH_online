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
thr_automatico = input('\Busqueda thr automatica? (1 = SI / 0 = NO) : ');

% Definimos manualmente un umbral para deteccion de spikes (en uV)
if thr_automatico == 0 
    thr = input('\nThreshold para el threshold cutting (en uV):  ');
end

% Guardo figuras?
guardar = input('\Guardo? (1 = SI / 0 = NO) : ');

% Guardo txt?
guardar_txt = input('\Guardo PSTHsw_1tet y LFP_1tet BOS? (1 = SI / 0 = NO) : ');

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
    amplifier_channels, frequency_parameters, puerto_canal_custom, 1000);

% % Genero n bandas de frecuencia equiespacidadas en log() superpuestas p
% n = 100; % cantidad de bandas superpuestas equiespaciadas
% p = 0.2; % superposición de las bandas (0 a 1)
% lista_bandas = list_multiple_bands(500, p, n, false);
% 
% % Inicializo vectores donde voy a guardar
% LFP_aux = zeros(length(LFP_tetrodo), 1); % vector momentaneo para calculos
% LFP_suma = zeros(length(LFP_tetrodo), 1); % suma de todas las bandas normalizadas
% LFP_list = zeros(length(LFP_tetrodo), n); % las bandas por separado
% 
% % Filtro y normalizo por bandas y luego vuelvo a sumar
% % para cada señal de cada canal
% for c = (1: size(spikes_canales,2))
%     figure()
%     for i = (1:size(lista_bandas,1))
%         % Filtro
%         LFP_aux = filt_and_normalize(LFP_canales(:,c), lista_bandas(i,1), lista_bandas(i,2), sr_lfp);
%         % Sumo banda al resto
%         LFP_suma = LFP_suma + LFP_aux;
%         fft_plot(LFP_aux, sr_lfp)
%         hold on
%     end
%     LFP_canales(:,c) = LFP_suma;
%     LFP_suma = zeros(length(LFP_tetrodo), 1);
% end


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
    
    csvwrite([directorio '/PSTHsw_1tet_BOS_' puerto_canal_custom '.txt'], ...
        PSTHsw_1tet_BOS)
    
    csvwrite([directorio '/LFP_1tet_BOS_' puerto_canal_custom '.txt'], ...
        LFP_1tet_BOS)
    
    filename = [directorio '/LFP_1tet_BOS_' puerto_canal_custom '.wav'];
    max_abs = max(abs(LFP_1tet_BOS));
    LFP_1tet_BOS_norm = (LFP_1tet_BOS / max_abs) * 0.9;
    audiowrite(filename,LFP_1tet_BOS_norm, sr_lfp)
    
end


% PLOTEO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Ploteo Grilla PSTH
plot_some_raster_LFP_1tetrode(grilla_psth, id_BOS, estimulos_tetrodos, ...
    frequency_parameters, sr_lfp, tiempo_file, ntrials,thr,directorio,spike_times);

suptitle2({datestr(now, 'yyyy-mm-dd'); ...
string(directorio) ; ...
strcat('tetrodo = ',string(puerto_canal_custom),"  |  ", ...
string(thr), "uV", "  |  ", "ntrials:", string(ntrials), "  |  ", ...
"t_inter_estimulo:", string(tiempo_file)) })

% % Espectrograma. no lo probé aun
% sample_rate=sr_lfp;
% window_width=300;
% [~,f,t,p]=spectrogram(saveme,gausswin(window_width,2.5),round(0.85*window_width),linspace(0,round(sample_rate/2),round(sample_rate/window_width)),sample_rate,'yaxis');
% %MinThreshold',-110. El minthreshold es opcional, lo uso para limitar la cantidad de gris que aparece.
% 
% % Y en general ploteo con esto:
% figure()
% imagesc('XData',t,'YData',f,'CData',10*log10(p(1:100,:)));
% colormap(flipud(gray));
% ylim([0 10000]);

% Otra estrategia espectrograma
figure()
nwin = 63;
wind = kaiser(nwin,17);
nlap = nwin-10;
nfft = 256;

spectrogram(LFP_1tet_BOS,wind,nlap,nfft,sr_lfp,'yaxis')
ylim([0 100])

figure()
fft_plot(LFP_1tet_BOS,sr_lfp);


% Guardo
if guardar == 1

    print_pdf(1, directorio, strcat('_',string(puerto_canal_custom),...
        '_grilla_PSTH-LFP-tetrode_BANDA_rara-30-50Hz_', string(round(thr)), 'uV.pdf'))
    
    print_pdf(3, directorio, strcat('_',string(puerto_canal_custom),...
    '_FFT-LFP-tetrode_BANDA_rara-30-50Hz_', string(round(thr)), 'uV.pdf'))
   
end

clear estimulos_aux j i  