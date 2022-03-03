% Calcula ISI violation de un cluster especifico

close all
clear all

% Cargo y defino parametros %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
refDur = 0.002; % duracion de tiempo refractario

% Defino directorios
directorio = input('Directorio INTAN: ','s');
directorio = horzcat(directorio , '/');

% Leer info INTAN
read_Intan_RHD2000_file(horzcat(directorio, 'info.rhd'));
clear notes spike_triggers supply_voltage_channels aux_input_channels 

directorio_clus = input('\nDirectorio clustering: ','s');
directorio_clus = horzcat(directorio_clus , '/');

% Pregunto que cluster grafico
clus = input('\nCluster ID: ');

% Defino canal a levantar con nombre custom name
peine = ...
    input('\nDefino canal a levantar (custom name) \n \nPeine (X): ');
tetrodo = input('\nTetrodo (X): ');
canal = input('\nCanal (X): ');

puerto_canal_custom = horzcat('P',num2str(peine),'-', ...
    'T',num2str(tetrodo),'-',num2str(canal));

% Traduzco custom_channel_name a native_channel_name
traduccion = strcmp(puerto_canal_custom, ...
    {amplifier_channels(:).custom_channel_name});

puerto_canal = amplifier_channels(traduccion).native_channel_name;
clear traduccion peine tetrodo canal

% Carga vector con parametros del analisis de datos
%params_info = dir(horzcat(directorio, '*parametros_protocolo*.txt'));
params = readtable(horzcat(directorio,'/parametros_protocolo.txt'), ...
    'Delimiter','\t');
clear params_info

% Cargo sampling rate de INTAN
sr = frequency_parameters.amplifier_sample_rate;

% Cargamos cantidad de trials y tiempo que dura cada uno
ntrials = params.Ntrials;
tiempo_file = params.tiempo_entre_estimulos;



% Proceso datos %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Levanto el canal de interes
raw = read_INTAN_channel(directorio, puerto_canal, amplifier_channels);

% Define el filtro para SPIKES
filt_spikes = designfilt('highpassiir','DesignMethod','butter',...
    'FilterOrder', 4,'HalfPowerFrequency',500,...
    'SampleRate',sr);

% Aplica filtros para quedarse con SPIKES y saca LFP
raw_filtered = filtfilt(filt_spikes, raw);

% Definimos umbral de deteccion de spikes
thr = 0;

% Levantamos spike times (en samples, NO unidades de tiempo) 
spike_times = readNPY([directorio_clus 'spike_times.npy']);
spike_clusters = readNPY([directorio_clus 'spike_clusters.npy']);

spike_clusters = spike_clusters == clus;
spike_times = spike_times(spike_clusters);
spike_times = double(spike_times);
spikeTrain = spike_times / sr;

[violationFraction, numViolations] = isiViolations(spikeTrain, refDur);
% 
% violationFraction
% 
% numViolations

plot_spikes_shapes(raw_filtered, spike_times, thr, frequency_parameters,...
    directorio);