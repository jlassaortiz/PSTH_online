% Calculo scores de varios directorios y ploteo los graficos tipo sabana

close all % Como ploteo y guardo, mejor asegurarme que no haya plots viejos
clear all

% Defino directorio donde esta archivo de parametros
directorio_params = input('Directorio parametros: ','s');
directorio_params = horzcat(directorio_params , '/');

% Guardo graficos sabanas de cada protocolo individual?
guardar_graf_protocolos_ind = ...
    input('\nGuardo graf sabanas de cada protocolo individual? (1 = SI / 0 = NO) : ');

% Grafico sabana o diag extendida?
no_diag = input('\n¿Ploteo diag ext? (1 = SI / 0 = NO) : ');

% Analizo por bandas?
bandas = input('\nFiltro banda particular? (1 = SI / 0 = NO) : ');
    if bandas == 1
        b_inf = input('\nLimite inferior (Hz): ');
        b_sup = input('\nLimite superior (Hz): ');
    else 
        b_inf = 0;
        b_sup = 400;
    end    

% Carga vector con parametros del analisis de datos
params_info = dir(horzcat(directorio_params, 'parametros.txt'));
params = readtable(horzcat(directorio_params,params_info.name),'Delimiter','\t', ...
    'ReadVariableNames',false);

% Posicion del primer directorio en el archivo de parametros
d = 10;

% Cargo "estimulos" usando el primer directorio de protocolos de la lista 
% de parametros
% TODOS LOS PROTOCOLOS DEBEN TENER LOS MIMOS ESTIMULOS
directorio_aux = horzcat(char(params.Var2(d)), '/');
estimulos = carga_songs(directorio_aux);

directorios = params(d:end, :);

% Genero diccionario donde se va a guardar el score de todos
score_total = struct;

% Para cada directorio (protocolo)
for j = (1:1:height(directorios))

    % Defino el directorio del protocolo
    directorio = horzcat(char(directorios.Var2(j)), '/') % directorio protocolo
    
    % Carga vector con parametros del analisis de datos
    params_info = dir(horzcat(directorio, 'parametros.txt'));
    params = readtable(horzcat(directorio,params_info.name),'Delimiter','\t', ...
        'ReadVariableNames',false);
    clear params_info

    % Cargo valores de puerto-canal
    puerto = char(params.Var2(1));
    canal = char(params.Var2(2));
    puerto_canal = horzcat(puerto, '-0', num2str(canal,'%.2d'))
    clear puerto canal

    % Cargamos cantidad de trials y tiempo que dura cada uno
    ntrials = str2num(char(params.Var2(3)))
    tiempo_file = str2num(char(params.Var2(4)))
    
    trials = (1:ntrials);

    % Especifico numero de id del BOS
    id_BOS = str2num(char(params.Var2(5)))

    % Cargo orden de la grilla
    grilla_sabana = str2num(string(params.Var2(6)))
    grilla_psth = str2num(string(params.Var2(7)))

    % Cargo el nombre de los parametros que varian por fila y columna de la grilla
    char(params.Var1(8))
    ejeX_fila = char(params.Var2(8))

    char(params.Var1(9))
    ejeY_col  = char(params.Var2(9))
    
    % Genero songs.mat a partir de las canciones
    estimulos = carga_songs(directorio);

    % Leer info INTAN
    read_Intan_RHD2000_file(horzcat(directorio, 'info.rhd'));
    clear notes spike_triggers supply_voltage_channels aux_input_channels 

    % Levanto el canal de interes
    raw = read_INTAN_channel(directorio, puerto_canal, amplifier_channels);

    % Define el filtro
    filt_spikes = designfilt('highpassiir','DesignMethod','butter','FilterOrder',...
        4,'HalfPowerFrequency',500,'SampleRate',frequency_parameters.amplifier_sample_rate);

    % Aplica filtro
    downsample_sr = 10000; % Hz
    raw_filtered = filtfilt(filt_spikes, raw);
    LFP = LFP_1channel(raw, frequency_parameters, ...
        downsample_sr, bandas, b_inf, b_sup);
    
    sr_lfp = downsample_sr;
    clear puerto canal filt_spikes raw

    % Genero diccionario con nombre de los estimulos y el momento de presentacion
    t0s_dictionary = find_t0s(estimulos, ntrials, tiempo_file, board_adc_channels, ...
        frequency_parameters, directorio, false, trials);
    
    % Definimos umbral de deteccion de spikes
    thr = find_thr(raw_filtered, t0s_dictionary, tiempo_file, frequency_parameters);

    % Buscamos spike por threshold cutting
    spike_times = find_spike_times(raw_filtered, thr, frequency_parameters);

    % Genero objeto con raster de todos los estimulos
    rasters = generate_raster(spike_times, t0s_dictionary, tiempo_file, ntrials, ...
        frequency_parameters);
    
    % Calculo LFP promediado por estimulo todos los trials 
    rasters = trialAverage_LFP(LFP, rasters, tiempo_file, ntrials, ...
    frequency_parameters, sr_lfp);

    % Evaluo desempleño de los distintos estimulos
    dict_score = score_calculator(id_BOS, rasters, frequency_parameters, ...
    spike_times, ntrials, tiempo_file);
    
    % Selecciono scores con los que me quedo y hago matriz para graficar
    [mat_scores, cell_estimulos] = scores_struct2mat_motif(grilla_sabana,dict_score, 1);
    
    % Agrego estos valores a la struct score_total que recopila todo
    score_total(j).id  = char(directorios.Var1(j)); % nombre corto protocolo
    score_total(j).dir = char(directorios.Var2(j)); % directorio protocolo
    score_total(j).grilla_scores = mat_scores; % array con valores XYZ para graficar
    score_total(j).grilla_nombre_estimulos = cell_estimulos; 
    % nombre estimulos para no perderles rastro
   
    % Sabana x4: integral y correlacion
    plot_sabana(mat_scores, directorio, ejeY_col, ejeX_fila); % Hace 4 plots

    % Grafica raster de todos los estimulos
    plot_some_raster(grilla_psth, id_BOS, estimulos, rasters, ...
    frequency_parameters, tiempo_file, ntrials, puerto_canal, thr, ...
    directorio, spike_times)
    
    % Guardo todo en el directorio del protocolo
    if guardar_graf_protocolos_ind == 1
        print_pdf(1, directorio, strcat('_sabana_INT_', string(round(thr)),'uV', '.pdf'))
        print_pdf(2, directorio, strcat('_sabana_CORR_', string(round(thr)),'uV', '.pdf'))
        print_pdf(3, directorio, strcat('_CORTE_sabana_INT_', string(round(thr)),'uV', '.pdf'))
        print_pdf(4, directorio, strcat('_CORTE_sabana_CORR_', string(round(thr)),'uV', '.pdf'))
        print_pdf(5, directorio, strcat('_grilla_', string(round(thr)),'uV', '.pdf')) 
    end
    
    close all
    clear amplifier_channels board_adc_channels frequency_parameters 
end

% Calculo el promedio de todas las grillas
mat_avg = zeros(numel(grilla_sabana), 4);

for i = (1:1:length(score_total))
    mat_avg = score_total(i).grilla_scores + mat_avg;   
end

mat_avg = mat_avg / length(score_total);





% Guardo todos los datos ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if no_diag == 1  % guardo datos de Diagonal Extendida
    score_total_toTable(score_total, directorio_params);

else % guardo datos de Grilla/Sabana

    sabana_INT_table = struct();

    % Inicializo contador que va a determinar el numero de fila de mi tabla
    % final 'sabana_INT_table'
    fila = 0;

    % para cada protocolo
    for p = (1:1:length({score_total.id}))

        % para cada estimulo
        for e = (1:length(score_total(p).grilla_scores(:,1)))

            x = score_total(p).grilla_scores(e,1);
            y = score_total(p).grilla_scores(e,2);

            % genero estimulo_id generico basado en tamaño cabeza y cuello
            if x == 1
                ca = 'caCh';
            elseif x == 2
                ca = 'caNo';
            else
                ca = 'caGr';
            end
            if y == 1
                cu = 'cuCh';
            elseif y == 2
                cu = 'cuNo';
            else 
                cu = 'cuGr';
            end

            estimulo_id = strcat(ca, '_', cu);

            % levanto valor del INT para estimulo e
            INT = score_total(p).grilla_scores(e,3);

            % actualizo numero fila y guardo la data
            fila = fila + 1;

            % guardo nombre protocolo
            sabana_INT_table(fila).protocolo = score_total(p).id;

            % guardo nombre del estimulo
            sabana_INT_table(fila).estimulo_id = estimulo_id;

            % guardo valor INT
            sabana_INT_table(fila).INT = INT;

            clear x y ca cu estimulo_id INT
        end 
    end 
    clear p e fila

    writetable(struct2table(sabana_INT_table), strcat(directorio_params, ...
        datestr(now, 'yyyy-mm-dd_HH_MM_SS'),'_tabla_all_datos.csv'));
end







% Ploteo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_some_sabana(score_total, mat_avg, ejeX_fila, ejeY_col);

% PLOTEO INT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if no_diag == 1 % agrego rayas para que se vea la diagonal extendida
    figure(1)
    hold on
    for i = (1:1:length(score_total))  
        X  = score_total(i).grilla_scores(:,1);
        Y  = score_total(i).grilla_scores(:,2);
        Z1 = score_total(i).grilla_scores(:,3);
        Z_all(:,i) = Z1;

        plot3(X,Y,Z1,'-')
        
        % if i < 11
        %     plot3(X,Y,Z1,'-r', 'MarkerSize',18,'LineWidth',2)
        %     hold on
        % elseif (i > 10 && i < 18)
        %     plot3(X,Y,Z1,'-g', 'MarkerSize',18,'LineWidth',2)
        %     hold on
        % else
        %     plot3(X,Y,Z1,'-b', 'MarkerSize',25,'LineWidth',2)
        %     hold on
        % end
        hold on
    end
    
    % Calculo error
    Z_std = zeros(size(Z_all, 1), 1);
    Z_mean = zeros(size(Z_all, 1), 1);
    for fila = (1:1:size(Z_all, 1))
        Z_std(fila,1) = std(Z_all(fila, :))/sqrt(length(Z_all));
        Z_mean(fila,1) = mean(Z_all(fila, :));
    end 
    errl = Z_mean - Z_std;
    errh = Z_mean + Z_std;

    plot3(X,Y,Z_mean, 'k.', 'MarkerSize',20 );
    plot3(X,Y,Z_mean, 'k-', 'LineWidth', 5 );
    plot3([X(:),X(:)]', [Y(:),Y(:)]', [errl(:),errh(:)]', '-r','LineWidth',5) 
    xticks([1 2 3 4 5])
    xticklabels({'0.5', '1.0', '1.5', '2.0', '2.5'})
    xlabel('Lambda')
    ylabel('INT')
    title('INT vs id lambda')
    view(0,0)
 
end


% PLOTEO LFP dif-corr %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if no_diag == 1 % agrego rayas para que se vea la diagonal extendida
    figure(2)
    hold on
    for i = (1:1:length(score_total))  
        X  = score_total(i).grilla_scores(:,1);
        Y  = score_total(i).grilla_scores(:,2);
        Z1 = score_total(i).grilla_scores(:,4);
        Z_all(:,i) = Z1;

        plot3(X,Y,Z1,'-')
        
        % if i < 11
        %     plot3(X,Y,Z1,'-r', 'MarkerSize',18,'LineWidth',2)
        %     hold on
        % elseif (i > 10 && i < 18)
        %     plot3(X,Y,Z1,'-g', 'MarkerSize',18,'LineWidth',2)
        %     hold on
        % else
        %     plot3(X,Y,Z1,'-b', 'MarkerSize',25,'LineWidth',2)
        %     hold on
        % end
        hold on
    end
    
    % Calculo error
    Z_std = zeros(size(Z_all, 1), 1);
    Z_mean = zeros(size(Z_all, 1), 1);
    for fila = (1:1:size(Z_all, 1))
        Z_std(fila,1) = std(Z_all(fila, :))/sqrt(length(Z_all));
        Z_mean(fila,1) = mean(Z_all(fila, :));
    end 
    errl = Z_mean - Z_std;
    errh = Z_mean + Z_std;

    plot3(X,Y,Z_mean, 'k.', 'MarkerSize',20 );
    plot3(X,Y,Z_mean, 'k-', 'LineWidth', 5 );
    plot3([X(:),X(:)]', [Y(:),Y(:)]', [errl(:),errh(:)]', '-r','LineWidth',5)    
    xticks([1 2 3 4 5])
    xticklabels({'0.5', '1.0', '1.5', '2.0', '2.5'})
    xlabel('Lambda')
    ylabel('LFP dif-corr')
    title(strcat('LFP dif-corr vs lambda | banda: ', ...
        num2str(b_inf), ' - ', num2str(b_sup), ' Hz'))
    view(0,0)
end



plot_sabana(mat_avg, directorio_params, ejeY_col, ejeX_fila);
set(gca, 'FontSize', 36)


% Guardo
print_pdf(1, directorio_params, strcat('_sabana_INT_', '.pdf'))
print_pdf(2, directorio_params, strcat('_sabana_CORR_', '.pdf'))
print_pdf(5, directorio_params, strcat('_CORTE_sabana_INT_', '.pdf'))
print_pdf(6, directorio_params, strcat('_CORTE_sabana_CORR_', '.pdf'))

clear ans params_info d directorio directorio_aux  puerto canal i j 
clear mat_scores cell_estimulos rasters raw raw_filtered spike_times
clear t0s_dictionary thr dict_score 
