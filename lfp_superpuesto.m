close all
clear all

%% Cargo los datos 

% Cargo directorios y nombre custom de cada protocolos a manopla
protocolos = {{'/Volumes/AUDIOS TEAM/Javi- finch durmierdo (NNx)/zf-JL037-VB/Datos/2021-12-20/zf-JL037-VB_p1_id1_211220_132913/',...
    '037-VB_p1_id1_211220'}, ...
    {'/Volumes/AUDIOS TEAM/Javi- finch durmierdo (NNx)/zf-JL037-VB/Datos/2021-12-21/zf-JL037-VB_p2_id1_211221_140619/',...
    '037-VB_p2_id1_211221'}, ...
    {'/Volumes/AUDIOS TEAM/Javi- finch durmierdo (NNx)/zf-JL037-VB/Datos/2021-12-21/zf-JL037-VB_p4_id1_211221_162306/',...
    '037-VB_p4_id1_211221'}};

datos = struct();
songs = struct();
i = 1;

% CUIDADO HARCODEO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dur_song = 4.5; % durancion del estímulo auditivo en 
sr_LFP = 1000; % sampling rate lfp en Hz


% Cargo nombre_id protocolo, directorios, nombres archivos con lfp y mua
for j = (1:length(protocolos))
    
    dir_aux = protocolos{j}{1};
    id_aux = protocolos{j}{2};
    
    % guardo cancion
    songs(j).protocolo = id_aux;
    songs(j).sound = carga_songs(dir_aux);
    
    % inicializo max_lfp y max_mua de cada protocolo
    max_env_lfp = 0;
    max_mua_protocolo = 0;
    
    env_lfp_int_aux = [];
    mua_int_aux = [];
    
    % Guardo archivos con LFP y MUA del BOS
    for p = (1:4)
        for t = (1:4)
            
            tetrodo = strcat('P', num2str(p),'-', 'T', num2str(t));

            id = strcat(id_aux,'_', tetrodo);
            file_name_lfp = dir(strcat(dir_aux, 'LFP_1tet_BOS_', tetrodo, '_BANDA-25-35Hz_*.txt')).name;
            file_name_mua = dir(strcat(dir_aux, 'PSTHsw_1tet_BOS_', tetrodo, '_BANDA-25-35Hz_*.txt')).name;

            datos(i).id = id;
            datos(i).protocolo = id_aux;
            datos(i).peine = p;
            datos(i).tetrodo = t;
            
            datos(i).dir_lfp = dir_aux;
            datos(i).file_lfp = strcat(dir_aux, file_name_lfp);
            
            datos(i).dir_mua = dir_aux;
            datos(i).file_mua = strcat(dir_aux, file_name_mua);
            
            % Levanto y guardo LFP
            lfp_aux = readtable(datos(i).file_lfp);
            lfp_aux = lfp_aux{:,1}; % paso de table a vector de double
            datos(i).lfp = lfp_aux;   
    
            % Calculo y guardo envolvente LFP 
            env_aux = abs(hilbert(lfp_aux));
            datos(i).env = env_aux;
            
            % Calculo integral de envolvente - ruido de la env_LFP
            env_LFP_int = 0;
            noise_aux = 0;
            lim_sound = uint64(dur_song*sr_LFP);
            for k = (1:lim_sound)
                env_LFP_int = env_LFP_int + env_aux(k,1);
            end
            for k = (lim_sound:lim_sound*2)
                noise_aux = noise_aux + env_aux(k,1);
            end
            env_LFP_int = env_LFP_int - noise_aux;
            datos(i).env_LFP_int = env_LFP_int;
            env_lfp_int_aux = [env_lfp_int_aux; env_LFP_int];
            
            
            % Calculo y guardo envolvente normalizada LFP 
            datos(i).env_norm = env_aux / max(env_aux);
            
            % Guardo amplitud max del tetrodo y de todos los tetrodos del protocolo
            datos(i).max_env_lfp = max(env_aux);
            if max(env_aux) > max_env_lfp
                max_env_lfp = max(env_aux);
                max_lfp_id = strcat('P', num2str(p), '-T', num2str(t));
            end

            % Levanto y guardo MUA
            mua_aux = readtable(datos(i).file_mua);
            mua_aux = table2array(mua_aux);
            datos(i).mua = mua_aux;
            
            % Calculo integral de envolvente - ruido de la MUA
            t_with_sound = mua_aux(:,2) < dur_song;
            t_without_sound = mua_aux(:,2) > dur_song & mua_aux(:,2) < dur_song*2;
            mua_int = sum(mua_aux(t_with_sound,1));
            noise_aux = sum(mua_aux(t_without_sound,1));
            mua_int = mua_int - noise_aux;
            datos(i).mua_int = mua_int; 
            mua_int_aux = [mua_int_aux; mua_int];

            % MUA SUAVIZADA
            mua_smooth = sgolayfilt(mua_aux(:,1), 4, 129);
            datos(i).mua_smooth = [mua_smooth, mua_aux(:,2)];

            mua_smooth2 = sgolayfilt(mua_aux(:,1), 4, 65);
            datos(i).mua_smooth2 = [mua_smooth2, mua_aux(:,2)];

            % Guardo MUA normalizada 
            datos(i).mua_norm = mua_aux;
            datos(i).mua_norm(:,1) = mua_aux(:,1)/max(mua_aux(:,1));
            
            % Guardo amplitud max del tetrodo y de todos los tetrodos del protocolo
            datos(i).max_mua = max(mua_aux(:,1));
            if max(mua_aux(:,1)) > max_mua_protocolo
                max_mua_protocolo = max(mua_aux(:,1));
                max_mua_id = strcat('P', num2str(p), '-T', num2str(t));
            end 
            
            i = i + 1;    
        end 
    end
    
    % Habiendo calculado el max LFP y MUA de todos los tetrodos, normalizo
    for k = (i-16 : i -1)
        datos(k).mua_max_protocolo = max_mua_protocolo;
        datos(k).mua_max_protocolo_id = max_mua_id;
        datos(k).mua_norm_protocolo = datos(k).mua;
        datos(k).mua_norm_protocolo(:,1) = (datos(k).mua(:,1) ./max_mua_protocolo) .*100;
        
        datos(k).mua_int_max_protocolo = max(mua_int_aux);
        datos(k).env_lfp_int_max_protocolo = max(env_lfp_int_aux);
        
        datos(k).env_lfp_max_protocolo = max_env_lfp;
        datos(k).env_lfp_max_protocolo_id = max_lfp_id;
        datos(k).env_lfp_norm_protocolo = (datos(k).env ./max_env_lfp) .*100;
    end   
end

clear id_aux p t i file_name_lfp file_name_mua tetrodo dir_aux j id max_lfp_id max_mua_id
clear mua_aux env_aux i j mua_smooth mua_smooth2 lfp_aux k max_lfp max_mua j noise_aux mua_int
clear mua_int_aux env_lfp_int_aux env_LFP_int


%% Guardo todos los datos

datos_all = datos;


%% Genero sub-set
close all

protocolo_analizar = 3; % indicar numero id del protocolo a analizar

% Separa automaticamente los datos del protocolo indicado
inicio_aux = (protocolo_analizar -1)*16 + 1;
fin_aux = inicio_aux + 15;
subset = (inicio_aux:1:fin_aux);
datos = datos_all(subset);

clear subset fin_aux inicio_aux

%% Corr Boari_2021

norm = input('\nCorrelacionamos señales normalizadas por max del protocolo? (1 = SI, 0 = NO): ');
norm = norm == 1;

corr_all_matrix_LFP = zeros(length(datos));
corr_all_matrix_MUA = zeros(length(datos));

% Matriz donde guardo valor max de cada señal de cada tetrodo de todo el peine
lfp_max = zeros(round(sqrt(length(datos))));
mua_max = zeros(round(sqrt(length(datos))));

for i = (1:1:length(datos))
    
    lfp_max(i) = datos(i).env_LFP_int;
    
    mua_max(i) = datos(i).mua_int;
    
end 

clear i 


% OJO HARCODEO FEO 
% Conservo seccion env lfp durante la presentacion del estimulo auditivo
fona_lfp = (1:1:length(datos(1).env))'; % vector de indices lfp
fona_lfp = fona_lfp < 4.5 * 1000; % transformo indices a t y me quedo con los menores a 4.5
% Conservo seccion mua durante la presentacion del estimulo
fona_mua = datos(1).mua(:,2) < 4.5;

% Calculo corr
labels = cell(1,length(datos));

peso_mua = zeros(length(datos)^2, 1);
peso_lfp = zeros(length(datos)^2, 1);

for i = (1:1:length(datos))
    
    % Hago lista de id de estimulos para usarlos como label en futuro grafico
    labels(1,i) = {datos(i).id};
    
    for j = (1:1:length(datos))
        % Calculo correlacion entre señales LFP
        if i == j
            corr_all_matrix_LFP(i,j) = 1;
        else
            if norm
                corr_all_matrix_LFP(i,j) = weighted_corr(datos(i).env(fona_lfp), ...
                    datos(j).env(fona_lfp), 100, ...
                    100*datos(i).env_LFP_int/datos(i).env_lfp_int_max_protocolo,  ...
                    100*datos(j).env_LFP_int/datos(j).env_lfp_int_max_protocolo);
            else
                corr_all_matrix_LFP(i,j) = weighted_corr(datos(i).env(fona_lfp), ...
                    datos(j).env(fona_lfp), ...
                    max(datos(i).env_lfp_int_max_protocolo, datos(j).env_lfp_int_max_protocolo),...
                    datos(i).env_LFP_int,  datos(j).env_LFP_int);
            end 
        end
        
        % Calculo correlacion entre señales MUA
        if i == j 
            corr_all_matrix_MUA(i,j) = 1; 
        else 
            if norm 
                corr_all_matrix_MUA(i,j) = weighted_corr(datos(i).mua(fona_mua,1), ...
                    datos(j).mua(fona_mua,1), 100, ...
                    100*datos(i).mua_int/datos(i).mua_int_max_protocolo, ...
                    100*datos(j).mua_int/datos(j).mua_int_max_protocolo);
            else
                corr_all_matrix_MUA(i,j) = weighted_corr(datos(i).mua(fona_mua,1), ...
                    datos(j).mua(fona_mua,1), ...
                    max(datos(i).mua_int_max_protocolo, datos(j).mua_int_max_protocolo), ...
                    datos(i).mua_int, datos(j).mua_int);
                
            end 
        end
    end
end




clear i j 

% plot LFP 
figure()
imagesc(corr_all_matrix_LFP)
colormap(gca,'parula');
colorbar();
caxis([0,1]); % or [-1,1]
ticks = 1:1:length(datos); 
set(gca,'TickLabelInterpreter','none')
set(gca, 'YTick', ticks, 'YTickLabel', labels);
set(gca, 'XTick', ticks, 'XTickLabel', labels, 'XTickLabelRotation',45);
axis equal
axis tight
title('correlacion pesada LFP')

% plot MUA 
figure()
imagesc(corr_all_matrix_MUA)
colormap(gca,'parula');
colorbar();
caxis([0,1]); % or [-1,1]
ticks = 1:1:length(datos); 
set(gca,'TickLabelInterpreter','none')
set(gca, 'YTick', ticks, 'YTickLabel', labels);
set(gca, 'XTick', ticks, 'XTickLabel', labels, 'XTickLabelRotation',45);
axis equal
axis tight
title('correlacion pesada MUA')


% plot max MUA 
figure()
imagesc(mua_max)
colormap(gca,'parula');
colorbar();
% caxis([0,1]); % or [-1,1]
labels_x = {'P1', 'P2', 'P3', 'P4'}; 
labels_y = {'T1', 'T2', 'T3', 'T4'};
ticks = (1:1:4);
set(gca,'TickLabelInterpreter','none')
set(gca, 'YTick', ticks, 'YTickLabel', labels_y);
set(gca, 'XTick', ticks, 'XTickLabel', labels_x);
axis equal
axis tight
title('INT de tetrodos de MUA (mientras fona - silencio)')

% plot max LFP 
figure()
imagesc(lfp_max)
colormap(gca,'parula');
colorbar();
% caxis([0,1]); % or [-1,1]
labels_x = {'P1', 'P2', 'P3', 'P4'}; 
labels_y = {'T1', 'T2', 'T3', 'T4'};
ticks = (1:1:4);
set(gca,'TickLabelInterpreter','none')
set(gca, 'YTick', ticks, 'YTickLabel', labels_y);
set(gca, 'XTick', ticks, 'XTickLabel', labels_x);
axis equal
axis tight
title('INT de tetrodos de LFP (mientras fona - silencio)')


%%

% Saco 1 de la diagonal
corr_all_matrix_LFP(corr_all_matrix_LFP > 0.99) = NaN;
corr_all_matrix_LFP
corr_all_matrix_MUA(corr_all_matrix_MUA > 0.99) = NaN;
corr_all_matrix_MUA

% Busco cual es la max corr y entre quienes se da
max_corr_LFP = max(corr_all_matrix_LFP ,[],'all')
[r, c] = find( corr_all_matrix_LFP == max_corr_LFP, 1)
max_corr_pair_LFP = {r, c ; datos(r).id, datos(c).id}

% Busco cual es la max corr y entre quienes se da
max_corr_MUA = max(corr_all_matrix_MUA, [], 'all')
[r, c] = find( corr_all_matrix_MUA == max_corr_MUA, 1)
max_corr_pair_MUA = {r, c; datos(r).id, datos(c).id}


% Busco cual es la min corr y entre quienes se da
min_corr_LFP = min(corr_all_matrix_LFP,[],'all')
[r, c] = find( corr_all_matrix_LFP == min_corr_LFP, 1)
min_corr_pair_LFP = {r, c ; datos(r).id, datos(c).id}

% Busco cual es la min corr y entre quienes se da
min_corr_MUA = min(corr_all_matrix_MUA,[],'all')
[r, c] = find( corr_all_matrix_MUA == min_corr_MUA, 1)
min_corr_pair_MUA = {r, c ; datos(r).id, datos(c).id}


%% Conservo solo las corr del tetrodo con mayor amplitud 
% Además calculo distancia del tetrodo con mayor amplitud al resto

% Dejo este proyecto en stand by

% Busco indice del tetrodo con mayor amplitud de lfp
index_id_max_LFP = 0;
for i = (1:length(datos))
    
    id_aux = datos(i).id;
    
    if id_aux == lfp_max_all_id
        index_id_max_LFP = i;
    end
end 

% De datos conservo la columna de correlaciones del tetrodo con mayor
% amplitud contra el resto
corr_tet_max_lfp = corr_all_matrix_LFP(:,index_id_max_LFP);





%% SONG vs MUA y LFP (CALCULOS)

indiceS_protocolo = [1,2,3]; % Protocolo a analizar
indiceS_BOS = [10, 10, 10]; % indice BOS
indiceS_tetrodo = [7, 22, 42]; % tetrodo a analizar

for i = (1:length(indiceS_protocolo))
    
    indice_protocolo = indiceS_protocolo(i);
    indice_BOS = indiceS_BOS(i);
    indice_tetrodo = indiceS_tetrodo(i);
    
    BOS_protocolo = songs(indice_protocolo).protocolo
    BOS_name = songs(indice_protocolo).sound(indice_BOS).name
    BOS_sr = songs(indice_protocolo).sound(indice_BOS).freq
    BOS_sound = songs(indice_protocolo).sound(indice_BOS).song;
    [BOS_times, BOS_env] = song_envelope(BOS_sound);
    BOS_times = BOS_times/BOS_sr;

    BOS_dur = (length(BOS_sound)/BOS_sr) + 0.5; % seconds

    LFP_times = (1:length(datos(indice_tetrodo).lfp))'/1000;
    lw = 2;

    
    
    
    % SONG vs MUA y LFP (PLOTEOS)

    % PLOT SONG - MUA
    figure()
    plot(BOS_times, BOS_env, 'LineWidth', lw)
    hold on
    plot(datos(indice_tetrodo).mua(:,2), datos(indice_tetrodo).mua(:,1)/max(datos(indice_tetrodo).mua(:,1)), ...
       ':','LineWidth', lw + 1)
    legend(BOS_name, strcat('MUA : ', datos(indice_tetrodo).id), 'Interpreter' , 'none');
    title(strcat(BOS_protocolo, ' | ' , datos(indice_tetrodo).id), 'Interpreter' , 'none')
    ylim([0, 1]);
    xlim([0, BOS_dur]);

    % PLOT SONG - LFP
    figure()
    plot(BOS_times, BOS_env, 'LineWidth', lw)
    hold on
    plot(LFP_times, datos(indice_tetrodo).env/max(datos(indice_tetrodo).env), ...
       ':','LineWidth', lw + 1 )
    legend(BOS_name, strcat('LFP : ', datos(indice_tetrodo).id), 'Interpreter' , 'none');
    title(strcat(BOS_protocolo, ' | ' , datos(indice_tetrodo).id), 'Interpreter' , 'none')
    ylim([0, 1]);
    xlim([0, BOS_dur]);


%     % PLOT SONG - MUA - LFP
%     figure()
%     plot(BOS_times, BOS_env, 'LineWidth', lw)
%     hold on
%     plot(datos(indice_tetrodo).mua(:,2), datos(indice_tetrodo).mua(:,1)/max(datos(indice_tetrodo).mua(:,1)),'LineWidth', lw + 1 )
%     plot(LFP_times, datos(indice_tetrodo).env/max(datos(indice_tetrodo).env),'LineWidth', lw + 1 )
%     legend(BOS_name, strcat('MUA : ', datos(indice_tetrodo).id),strcat('LFP : ', datos(indice_tetrodo).id), 'Interpreter' , 'none');
%     title(strcat(BOS_protocolo, ' | ' , datos(indice_tetrodo).id), 'Interpreter' , 'none')
%     ylim([0, 1]);
%     xlim([0, BOS_dur]);

%     %Guardo
%     guardar = input('\nGuardo? (1 = SI / 0 = NO) : ');
%     if guardar
%         dir = directorios{indice_protocolo}{1};
%         id = datos(indice_protocolo).id;
% 
%         print_pdf(1, dir, strcat('_env-BOS_vs_MUA_', id, '.pdf'))
%         print_pdf(2, dir, strcat('_env-BOS_vs_LFP_', id, '.pdf'))
%         print_pdf(3, dir, strcat('_env-BOS_vs_MUA_LFP_', id, '.pdf'))
%     end 
% 
%     clear guardar
    
end 


%% PLOTEO LFP, LFP_env, LFP_env_norm_protocolo
close all

% Elijo quienes plotear
ploteo = [7, 22, 45];
dur_sound = 5*1000; % para determinar xlim()


% Guardo figuras?
guardar = input('\nGuardo? (1 = SI / 0 = NO) : ');

% Defino directorio donde guardo (como son correlaciones de dos archivos
% distintos, guardo en primer directorio para evitar duplicados
dir_aux = datos(ploteo(1)).dir_lfp

% Defino alpha (transparencia) y LineWidth
alpha = 0.50;
lw = 2;


% LFP 
figure()
leyendas = cell(length(ploteo),1);
for k = (1:1:length(ploteo))
    
    p= plot(datos(ploteo(k)).lfp, 'LineWidth', lw);
    hold on

    p.Color(4) = alpha;

    leyendas{k,1} = datos(ploteo(k)).id;
end 
legend(leyendas, 'Interpreter' , 'none');
title('LFP banda 25-35Hz promediando tetrodo para el BOS')
set(gca,'FontSize',25)
xlim([0, dur_sound])


% LFP env
figure()
for k = (1:1:length(ploteo))
    
    % Ploteo envolvente
    p= plot(datos(ploteo(k)).env, 'LineWidth', lw);
    hold on
    
    p.Color(4) = alpha;
end 
legend(leyendas, 'Interpreter' , 'none');
title('Envolvente LFP banda 25-35Hz promediando tetrodo para el BOS')
set(gca,'FontSize',25)
xlim([0, dur_sound])


% LFP env normalizada por protocolo
figure()
for k = (1:1:length(ploteo))
    
    % Ploteo envolvente
    p= plot(datos(ploteo(k)).env_lfp_norm_protocolo, 'LineWidth', lw);
    hold on
    
    p.Color(4) = alpha;
end 
legend(leyendas, 'Interpreter' , 'none');
title('Envolvente LFP NORMALIZADA POR PROTOCOLO banda 25-35Hz promediando tetrodo para el BOS')
set(gca,'FontSize',25)
xlim([0, dur_sound])

% Guardo
if guardar
    titulo = '_LFP_promedioTetrodo_BOS';

    print_pdf(1, dir_aux, strcat(titulo, '.pdf'))
    print_pdf(2, dir_aux, strcat(titulo,'_env.pdf'))
    print_pdf(3, dir_aux, strcat(titulo,'_env_NORM-PROTOCOLO.pdf'))
end

clear directorio 


%% PLOTEO MUA, MUA_norm_protocolo

close all

% Elijo quienes plotear
ploteo = [7, 28, 42];

% Guardo figuras?
guardar = input('\nGuardo? (1 = SI / 0 = NO) : ');
dur_sound = 5;

% Defino directorio donde guardo (como son correlaciones de dos archivos
% distintos, guardo en primer directorio para evitar duplicados
dir_aux = datos(ploteo(1)).dir_lfp

% Defino alpha (transparencia) y LineWidth
alpha = 0.50;
lw = 2;


% MUA 
figure()
for k = (1:1:length(ploteo))
    
    % Ploteo envolvente
    p= plot(datos(ploteo(k)).mua(:,2), datos(ploteo(k)).mua(:,1), 'LineWidth', lw);
    hold on
    
    leyendas{k,1} = datos(ploteo(k)).id;
    
    p.Color(4) = alpha;
end 
legend(leyendas, 'Interpreter' , 'none');
title('MUA promediando tetrodo para el BOS')
set(gca,'FontSize',25)
xlim([0, dur_sound])


% MUA normalizada por protocolo
figure()
for k = (1:1:length(ploteo))
    
    % Ploteo envolvente
    p= plot(datos(ploteo(k)).mua_norm_protocolo(:,2), datos(ploteo(k)).mua_norm_protocolo(:,1), 'LineWidth', lw);
    hold on
    
    p.Color(4) = alpha;
end 
legend(leyendas, 'Interpreter' , 'none');
title('MUA NORMALIZADA POR PROTOCOLO promediando tetrodo para el BOS')
set(gca,'FontSize',25)
xlim([0, dur_sound])


% Guardo
if guardar
    titulo = '_MUA_promedioTetrodo_BOS';
    print_pdf(1, dir_aux, strcat(titulo, '.pdf'))
    print_pdf(2, dir_aux, strcat(titulo,'_NORM-PROTOCOLO.pdf'))
end 



