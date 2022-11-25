
close all
clear all

% Cargo y defino parametros %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Defino directorio
directorio = input('Directorio: ','s');
directorio = horzcat(directorio , '/');

% Leer info INTAN
read_Intan_RHD2000_file(horzcat(directorio, 'info.rhd'));
clear notes spike_triggers supply_voltage_channels aux_input_channels 

% Pregunto si determino canal de manera manual o automatica (esta
% explicitado en parametros que se cargan)

% Defino canal a levantar con nombre custom name
peine = ...
    input('\nDefino canal a levantar (custom name) \n \nPeine (X): ');
tetrodo = input('\nTetrodo (X): ');

puerto_canal_custom = horzcat('P',num2str(peine),'-', ...
    'T',num2str(tetrodo),'-');

% Levanto el canal de interes
for c = (1:4)
    
    % Traduzco custom_channel_name a native_channel_name
    traduccion = strcmp([puerto_canal_custom num2str(c)], ...
        {amplifier_channels(:).custom_channel_name});

    puerto_canal = amplifier_channels(traduccion).native_channel_name;
    clear traduccion peine tetrodo canal

    raw = read_INTAN_channel(directorio, puerto_canal, amplifier_channels);

    % Define el filtro para SPIKES y LFP
    filt_spikes = designfilt('highpassiir','DesignMethod','butter',...
        'FilterOrder', 4,'HalfPowerFrequency',500,...
        'SampleRate',frequency_parameters.amplifier_sample_rate);

    % Aplica filtros
    raw_filtered = filtfilt(filt_spikes, raw);
    
    t = (1:length(raw_filtered))*frequency_parameters.amplifier_sample_rate;
    
    ax(c) = subplot(4,1,c);
    plot(t', raw_filtered);
    title([puerto_canal_custom num2str(c)], 'Interpreter', 'None', 'Fontsize', 20)
    ylim([-250 200]);
end
linkaxes(ax, 'x');

suptitle2(directorio);

