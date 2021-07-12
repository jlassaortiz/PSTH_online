% Calculo scores de varios directorios y ploteo los graficos tipo sabana

close all % Como ploteo y guardo, mejor asegurarme que no haya plots viejos
clear all

% Defino directorio donde esta archivo de parametros
directorio_params = input('Directorio parametros: ','s');
directorio_params = horzcat(directorio_params , '/');

% Guardo graficos sabanas de cada protocolo individual?
guardar_graf_protocolos_ind = input('\nGuardo graf sabanas de cada protocolo individual? (1 = SI / 0 = NO) : ');

% Carga vector con parametros del analisis de datos
params_info = dir(horzcat(directorio_params, '*parametros*.txt'));
params = readtable(horzcat(directorio_params,params_info.name),'Delimiter','\t','ReadVariableNames',false);

% Posicion del primer directorio en el archivo de parametros
d = 1;

directorios = params(d:end, :);

% Genero diccionario donde se va a guardar el score de todos
score_total = struct;

clear d directorio_aux

% Para cada directorio (protocolo)
for j = (1:1:height(directorios))

    % Defino el directorio del protocolo
    directorio = horzcat(char(directorios.Var2(j)), '/') % directorio protocolo
    
    % Estraigo nombre corto directorio
    directorio_nombre_corto = char(directorios.Var1(j));
    
    % Carga vector con parametros del analisis de datos
    params_info = dir(horzcat(directorio, '*parametros_protocolo*.txt'));
%     params = readtable(horzcat(directorio,params_info.name),'Delimiter','\t','ReadVariableNames',false);
    params = readtable(horzcat(directorio,params_info.name),'Delimiter','\t');
    clear params_info
    
    % Carga vector con parametros del analisis de datos
    params_info = dir(horzcat(directorio, '*parametros_analisis*.txt'));
    params_analisis = readtable(horzcat(directorio,params_info.name),'Delimiter','\t');
    clear params_info

    % Cargo valores de puerto-canal
    puerto = char(params.Puerto);
    canal = params.Canal;
    puerto_canal = horzcat(puerto, '-0', num2str(canal,'%.2d'))
    clear puerto canal

    % Cargamos cantidad de trials y tiempo que dura cada uno
    ntrials = params.Ntrials
    tiempo_file = params.tiempo_entre_estimulos
    
    % Especifico numero de id del BOS y REV
    id_BOS = params_analisis.id_bos(1)
    id_REV = params_analisis.id_rev(1)
    
    % Cargo orden de la grilla
    grilla_psth = str2num(string(params_analisis.grilla_psth(1)))
    
    % Genero songs.mat a partir de las canciones
    estimulos = carga_songs(directorio);    
    
    % cargo id_estimulos 
    for i = (1:1:length(estimulos))
        estimulos(i).id = params_analisis.orden(i);
        estimulos(i).frec_corte = params_analisis.freq_corte(i);
        estimulos(i).tipo = categorical(params_analisis.tipo_estimulo(i));
        estimulos(i).protocolo_id = categorical({directorio_nombre_corto});
    end
    clear i 
   
    % Leer info INTAN
    read_Intan_RHD2000_file(horzcat(directorio, 'info.rhd'));
    clear notes spike_triggers supply_voltage_channels aux_input_channels 

    % Levanto el canal de interes
    raw = read_INTAN_channel(directorio, puerto_canal, amplifier_channels);

    % Define el filtro
    filt_spikes = designfilt('highpassiir','DesignMethod','butter','FilterOrder',...
        4,'HalfPowerFrequency',500,'SampleRate',frequency_parameters.amplifier_sample_rate);

    % Aplica filtro
    raw_filtered = filtfilt(filt_spikes, raw);
    clear filt_spikes

    % Genero diccionario con nombre de los estimulos y el momento de presentacion
    estimulos = find_t0s(estimulos, ntrials, tiempo_file, board_adc_channels, frequency_parameters, directorio, false);

    % Definimos umbral de deteccion de spikes
    thr = find_thr(raw_filtered, estimulos, tiempo_file, frequency_parameters);

    % Buscamos spike por threshold cutting
    spike_times = find_spike_times(raw_filtered, thr, frequency_parameters);

    % Genero objeto con raster de todos los estimulos
    estimulos = generate_raster(spike_times, estimulos , tiempo_file, ntrials, frequency_parameters);

    % Calculo scores
    estimulos = score_calculator(id_BOS, id_REV, estimulos, frequency_parameters);
    
    % Selecciono sub-set de estimulos para guardar
    estimulos_resumen = rmfield(estimulos, {'song', 't0s','spikes_norm', 'trials_id', 'id'});
    estimulos_resumen = struct2table(estimulos_resumen);
    
    if j == 1
        score_total_tabla = estimulos_resumen;
    else
        score_total_tabla = [score_total_tabla ; estimulos_resumen];
    end 
    
    % Agrego estos valores a la struct score_total que recopila todo
    score_total(j).id  = char(directorios.Var1(j)); % nombre corto protocolo
    score_total(j).dir = char(directorios.Var2(j)); % directorio protocolo
    score_total(j).scores = estimulos_resumen; % nombre estimulos y respectivos scores
    
    
    % Guardo todo en el directorio del protocolo
    if guardar_graf_protocolos_ind == 1

        % Selecciono datos de ese protocolo 
        pasa_altos = estimulos_resumen(estimulos_resumen.tipo == 'up' , :);
        pasa_bajos = estimulos_resumen(estimulos_resumen.tipo == 'down' , :);

        % Plotear INT PASA-ALTOS
        figure();
        plot(pasa_altos.frec_corte, pasa_altos.int_norm, '-o', 'DisplayName', score_total(j).id)
        title({strcat('INT_PASA-ALTOS_', datestr(now, 'yyyy-mm-dd')); ...
        string(directorio) ; ...
        strcat(string(puerto_canal), "  " , string(thr), "uV", "  ntrials:", string(ntrials), "  t_inter_estimulo:", string(tiempo_file)) }, 'Interpreter','None')
        legend
        
        % Plotear INT PASA-BAJOS
        figure();
        plot(pasa_bajos.frec_corte, pasa_bajos.int_norm, '-o','DisplayName', score_total(j).id)
        title({strcat('INT_PASA-BAJOS_', datestr(now, 'yyyy-mm-dd')); ...
        string(directorio) ; ...
        strcat(string(puerto_canal), "  " , string(thr), "uV", "  ntrials:", string(ntrials), "  t_inter_estimulo:", string(tiempo_file)) }, 'Interpreter','None')
        legend
        
        % Plotear CORR PASA-ALTOS
        figure();
        plot(pasa_altos.frec_corte, pasa_altos.corr, '-o', 'DisplayName', score_total(j).id)
        title({strcat('CORR_PASA-ALTOS_', datestr(now, 'yyyy-mm-dd')); ...
        string(directorio) ; ...
        strcat(string(puerto_canal), "  " , string(thr), "uV", "  ntrials:", string(ntrials), "  t_inter_estimulo:", string(tiempo_file)) }, 'Interpreter','None')
        legend
        
        % Plotear CORR PASA-BAJOS
        figure();
        plot(pasa_bajos.frec_corte, pasa_bajos.corr, '-o','DisplayName', score_total(j).id)
        title({strcat('CORR_PASA-BAJOS_', datestr(now, 'yyyy-mm-dd')); ...
        string(directorio) ; ...
        strcat(string(puerto_canal), "  " , string(thr), "uV", "  ntrials:", string(ntrials), "  t_inter_estimulo:", string(tiempo_file)) }, 'Interpreter','None')
        legend
        
        % Guardo
        print_pdf(1, directorio, strcat('_INT_pasa-ALTOS', '.pdf'))
        print_pdf(2, directorio, strcat('_INT_pasa-BAJOS', '.pdf'))
        print_pdf(3, directorio, strcat('_CORR_pasa-ALTOS', '.pdf'))
        print_pdf(4, directorio, strcat('_CORR_pasa-BAJOS', '.pdf'))
    end
    
    close all
    clear amplifier_channels board_adc_channels frequency_parameters estimulos_aux estimulos estimulos_resumen directorio_nombre_corto
    clear grilla_psth j i id_BOS ntrials params params_analisis pasa_altos pasa_bajos puerto_canal raw raw_filtered spike_times thr tiempo_file directorio

end



% CONVIERTO FREC CORTES CONTROLES Y DEMAS EN CATEGORIAS DISTINTAS DE CERO 
score_total_tabla.frec_corte = categorical(score_total_tabla.frec_corte);
for i = (1:1:height(score_total_tabla))
    if score_total_tabla.frec_corte(i) == categorical(0)
        score_total_tabla.frec_corte(i) = score_total_tabla.tipo(i);
    end
end



% Calculo media y mediana
score_statistics = struct();

frec_cortes = categories(score_total_tabla.frec_corte);
tipo = categories(score_total_tabla.tipo);
k = 1;

for i = (1:1:length(frec_cortes))
    
    for j = (1:1:length(tipo))
        
        f = categorical(frec_cortes(i,1));
        t = categorical(tipo(j,1));
        chosen = score_total_tabla.frec_corte == f & score_total_tabla.tipo == t;
        
        if sum(chosen) > 0
            
            score_statistics(k).frec_cortes = f;
            score_statistics(k).tipo = t;
        
            score_statistics(k).mean_int = mean(score_total_tabla.int_norm(score_total_tabla.frec_corte == f & score_total_tabla.tipo == t));
            score_statistics(k).mean_corr = mean(score_total_tabla.corr(score_total_tabla.frec_corte == f & score_total_tabla.tipo == t));

            score_statistics(k).median_int = median(score_total_tabla.int_norm(score_total_tabla.frec_corte == f & score_total_tabla.tipo == t));
            score_statistics(k).median_corr = median(score_total_tabla.corr(score_total_tabla.frec_corte == f & score_total_tabla.tipo == t));

            score_statistics(k).sd_int = std(score_total_tabla.int_norm(score_total_tabla.frec_corte == f & score_total_tabla.tipo == t));
            score_statistics(k).sd_corr = std(score_total_tabla.corr(score_total_tabla.frec_corte == f & score_total_tabla.tipo == t));
            
            k = k+1;
        end
    end 
end

score_statistics = struct2table(score_statistics);

pasa_altos_stat = score_statistics(score_statistics.tipo == categorical({'up'}),:);
pasa_altos_stat.frec_cortes = double(string(pasa_altos_stat.frec_cortes));
pasa_bajos_stat = score_statistics(score_statistics.tipo == categorical({'down'}),:);
pasa_bajos_stat.frec_cortes = double(string(pasa_bajos_stat.frec_cortes));

error_width = 6;







% PLOTEO !

tamano_letra = 12*3;

% Grafico INT PASA-ALTOS
figure()
for i = (1:1:height(directorios))
    
    % Levanto nombre de protocolo
    dir_id = score_total(i).id;
    
    % Levanto tabla de scores de ese protocolo
    tabla_scores = score_total(i).scores;
    
    % Selecciono datos de ese protocolo 
    pasa_altos = tabla_scores(tabla_scores.tipo == 'up' , :);

    % Plotear por separado (en dos figuras) pasa altos y pasa bajo con
    % leyenda con PROTOCOLO_ID y titulo PASA ALTOS / BAJOS
    plot(pasa_altos.frec_corte, pasa_altos.int_norm, '-o','DisplayName', char(dir_id))
    
    hold on;
end
errorbar(pasa_altos_stat.frec_cortes, pasa_altos_stat.mean_int, pasa_altos_stat.sd_int,'LineWidth', error_width, 'DisplayName', 'MEAN INT pasa ALTOS')
xlim([0.5, 6.5])
xlabel('frecuencia de corte (kHz)');
ylabel('integral normalizada');
set(gca,'FontSize',tamano_letra)
title({strcat('INT_PASA-ALTOS_', datestr(now, 'yyyy-mm-dd')); ...
string(directorio_params)}, 'Interpreter','None', 'FontSize', 24)
legend('FontSize', 10)


% Grafico INT PASA-BAJOS
figure()
for i = (1:1:height(directorios))
    
    % Levanto nombre de protocolo
    dir_id = score_total(i).id;
    
    % Levanto tabla de scores de ese protocolo
    tabla_scores = score_total(i).scores;
    
    % Selecciono datos de ese protocolo 
    pasa_bajos = tabla_scores(tabla_scores.tipo == 'down' , :);

    % Plotear por separado (en dos figuras) pasa altos y pasa bajo con
    % leyenda con PROTOCOLO_ID y titulo PASA ALTOS / BAJOS
    plot(pasa_bajos.frec_corte, pasa_bajos.int_norm, '-o','DisplayName', char(dir_id))
    
    hold on;
end
errorbar(pasa_bajos_stat.frec_cortes, pasa_bajos_stat.mean_int, pasa_bajos_stat.sd_int,'LineWidth', error_width, 'DisplayName', 'MEAN INT pasa BAJOS')
xlim([0.5, 6.5])
xlabel('frecuencia de corte (kHz)');
ylabel('integral normalizada');
set(gca,'FontSize',tamano_letra)
title({strcat('INT_PASA-BAJOS_', datestr(now, 'yyyy-mm-dd')); ...
string(directorio_params)}, 'Interpreter','None', 'FontSize', 24)
legend('Location','southeast','FontSize', 10)


clear tabla_scores dir_id


% Grafico CORR PASA-ALTOS
figure()
for i = (1:1:height(directorios))
    
    % Levanto nombre de protocolo
    dir_id = score_total(i).id;
    
    % Levanto tabla de scores de ese protocolo
    tabla_scores = score_total(i).scores;
    
    % Selecciono datos de ese protocolo 
    pasa_altos = tabla_scores(tabla_scores.tipo == 'up' , :);

    % Plotear por separado (en dos figuras) pasa altos y pasa bajo con
    % leyenda con PROTOCOLO_ID y titulo PASA ALTOS / BAJOS
    plot(pasa_altos.frec_corte, pasa_altos.corr, '-o','DisplayName', char(dir_id))
    
    hold on;
end
errorbar(pasa_altos_stat.frec_cortes, pasa_altos_stat.mean_corr, pasa_altos_stat.sd_corr,'LineWidth', error_width, 'DisplayName', 'MEAN CORR pasa ALTOS')
xlim([0.5, 6.5])
xlabel('frecuencia de corte (kHz)');
ylabel('Correlacion');
set(gca,'FontSize',tamano_letra)
title({strcat('CORR_PASA-ALTOS_', datestr(now, 'yyyy-mm-dd')); ...
string(directorio_params)}, 'Interpreter','None', 'FontSize', 24)
legend('FontSize', 10)



% Grafico CORR PASA-BAJOS
figure()
for i = (1:1:height(directorios))
    
    % Levanto nombre de protocolo
    dir_id = score_total(i).id;
    
    % Levanto tabla de scores de ese protocolo
    tabla_scores = score_total(i).scores;
    
    % Selecciono datos de ese protocolo 
    pasa_bajos = tabla_scores(tabla_scores.tipo == 'down' , :);

    % Plotear por separado (en dos figuras) pasa altos y pasa bajo con
    % leyenda con PROTOCOLO_ID y titulo PASA ALTOS / BAJOS
    plot(pasa_bajos.frec_corte, pasa_bajos.corr, '-o','DisplayName', char(dir_id))
    
    hold on;
end
errorbar(pasa_bajos_stat.frec_cortes, pasa_bajos_stat.mean_corr, pasa_bajos_stat.sd_corr,'LineWidth', error_width, 'DisplayName', 'MEAN CORR pasa BAJOS')
xlim([0.5, 6.5])
xlabel('frecuencia de corte (kHz)');
ylabel('Correlacion');
set(gca,'FontSize',tamano_letra)
title({strcat('CORR_PASA-BAJOS_', datestr(now, 'yyyy-mm-dd')); ...
string(directorio_params)}, 'Interpreter','None', 'FontSize', 24)
legend('Location','southeast','FontSize', 10)


clear tabla_scores dir_id



% MEAN ALL INT
figure()
errorbar(pasa_altos_stat.frec_cortes, pasa_altos_stat.mean_int, pasa_altos_stat.sd_int,'LineWidth', 2, 'DisplayName', 'MEAN INT pasa ALTOS')
hold on
errorbar(pasa_bajos_stat.frec_cortes, pasa_bajos_stat.mean_int, pasa_bajos_stat.sd_int,'LineWidth', 2, 'DisplayName', 'MEAN INT pasa BAJOS')
xlim([0.5, 6.5])
xlabel('frecuencia de corte (kHz)');
ylabel('integral normalizada');
set(gca,'FontSize',tamano_letra)
title({strcat('MEAN_ALL_INT_', datestr(now, 'yyyy-mm-dd')); ...
string(directorio_params)}, 'Interpreter','None', 'FontSize', 24)
legend('Location','southeast','FontSize', 20)


% MEAN ALL CORR
figure()
errorbar(pasa_altos_stat.frec_cortes, pasa_altos_stat.mean_corr, pasa_altos_stat.sd_corr,'LineWidth', 2, 'DisplayName', 'MEAN CORR pasa ALTOS')
hold on
errorbar(pasa_bajos_stat.frec_cortes, pasa_bajos_stat.mean_corr, pasa_bajos_stat.sd_corr,'LineWidth', 2, 'DisplayName', 'MEAN CORR pasa BAJOS')
xlim([0.5, 6.5])
ylabel('Correlacion');
xlabel('frecuencia de corte (kHz)');
set(gca,'FontSize',tamano_letra)
title({strcat('MEAN_ALL_CORR_', datestr(now, 'yyyy-mm-dd')); ...
string(directorio_params)}, 'Interpreter','None', 'FontSize', 24)
legend('Location','southeast','FontSize', 20)



% BOX PLOT
pasa_altos_box = score_total_tabla(score_total_tabla.tipo == 'up' | score_total_tabla.tipo == 'CON' |score_total_tabla.tipo == 'REV',:);
pasa_bajos_box = score_total_tabla(score_total_tabla.tipo == 'down'| score_total_tabla.tipo == 'CON' |score_total_tabla.tipo == 'REV',:);

% BOXPLOT INT pasa ALTOS
figure()
boxplot(pasa_altos_box.int_norm,pasa_altos_box.frec_corte)
xlabel('frecuencia de corte (kHz)');
ylabel('integral normalizada');
set(gca,'FontSize',tamano_letra)
title({strcat('BOX_INT_pasa_ALTOS_', datestr(now, 'yyyy-mm-dd')); ...
string(directorio_params)}, 'Interpreter','None', 'FontSize', 24)
legend('Location','southeast','FontSize', 10)


% Boxplot INT pasa BAJOS
figure()
boxplot(pasa_bajos_box.int_norm,pasa_bajos_box.frec_corte)
xlabel('frecuencia de corte (kHz)');
ylabel('integral normalizada');
set(gca,'FontSize',tamano_letra)
title({strcat('BOX_INT_pasa_BAJOS_', datestr(now, 'yyyy-mm-dd')); ...
string(directorio_params)}, 'Interpreter','None', 'FontSize', 24)
legend('Location','southeast','FontSize', 10)


% BOXPLOT CORR pasa ALTOS
figure()
boxplot(pasa_altos_box.corr,pasa_altos_box.frec_corte)
xlabel('frecuencia de corte (kHz)');
ylabel('Correlacion');
set(gca,'FontSize',tamano_letra)
title({strcat('BOX_CORR_pasa_ALTOS_', datestr(now, 'yyyy-mm-dd')); ...
string(directorio_params)}, 'Interpreter','None', 'FontSize', 24)
legend('Location','southeast','FontSize', 10)


% Boxplot CORR pasa BAJOS
figure()
boxplot(pasa_bajos_box.corr,pasa_bajos_box.frec_corte)
xlabel('frecuencia de corte (kHz)');
ylabel('Correlacion');
set(gca,'FontSize',tamano_letra)
title({strcat('BOX_CORR_pasa_BAJOS_', datestr(now, 'yyyy-mm-dd')); ...
string(directorio_params)}, 'Interpreter','None', 'FontSize', 24)
legend('Location','southeast','FontSize', 10)



% Guardo
print_pdf(1, directorio_params, '_INT_pasa_ALTOS.pdf')
print_pdf(2, directorio_params, '_INT_pasa_BAJOS.pdf')
print_pdf(3, directorio_params, '_CORR_pasa_ALTOS.pdf')
print_pdf(4, directorio_params, '_CORR_pasa_BAJOS.pdf')
print_pdf(5, directorio_params, '_MEAN_ALL_INT.pdf')
print_pdf(6, directorio_params, '_MEAN_ALL_CORR.pdf')
print_pdf(7, directorio_params, '_BOX_INT_pasa_ALTOS.pdf')
print_pdf(8, directorio_params, '_BOX_INT_pasa_BAJOS.pdf')
print_pdf(9, directorio_params, '_BOX_CORR_pasa_ALTOS.pdf')
print_pdf(10, directorio_params, '_BOX_CORR_pasa_BAJOS.pdf')

clear i j k 
