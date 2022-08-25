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
dur_song = 4.5; % durancion del est�mulo auditivo en
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
    max_mua_protocolo_smooth = 0;

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
            datos(i).coordinates = [(p-1)*200 , (t-1)*150];

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
            for k = (lim_sound:length(env_aux))
                noise_aux = noise_aux + env_aux(k,1);
            end
            env_LFP_int = env_LFP_int - noise_aux;
            datos(i).env_LFP_int = env_LFP_int;

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
            t_without_sound = mua_aux(:,2) > dur_song;
            mua_int = sum(mua_aux(t_with_sound,1));
            noise_aux = sum(mua_aux(t_without_sound,1));
            mua_int = mua_int - noise_aux;
            datos(i).mua_int = mua_int;

            % MUA SUAVIZADA
            mua_smooth = sgolayfilt(mua_aux(:,1), 4, 65);
            datos(i).mua_smooth = [mua_smooth, mua_aux(:,2)];

            % Guardo MUA normalizada
            datos(i).mua_norm = mua_aux;
            datos(i).mua_norm(:,1) = mua_aux(:,1) ./ max(mua_aux(:,1));

            % Guardo amplitud max del tetrodo y de todos los tetrodos del protocolo
            datos(i).max_mua = max(mua_aux(:,1));
            if max(mua_aux(:,1)) > max_mua_protocolo
                max_mua_protocolo = max(mua_aux(:,1));
                max_mua_id = strcat('P', num2str(p), '-T', num2str(t));
            end

            % Guardo amplitud max del tetrodo y de todos los tetrodos del
            % protocolo SMOOTH
            datos(i).max_mua_smooth = max(mua_smooth(:,1));
            if max(mua_smooth(:,1)) > max_mua_protocolo_smooth
                max_mua_protocolo_smooth = max(mua_smooth(:,1));
                max_mua_id_smooth = strcat('P', num2str(p), '-T', num2str(t));
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

        datos(k).mua_max_protocolo_smooth = max_mua_protocolo_smooth;
        datos(k).mua_max_protocolo_id_smooth = max_mua_id_smooth;
        datos(k).mua_norm_protocolo_smooth = datos(k).mua_smooth;
        datos(k).mua_norm_protocolo_smooth(:,1) = (datos(k).mua_smooth(:,1) ./max_mua_protocolo_smooth) .*100;

        datos(k).env_lfp_max_protocolo = max_env_lfp;
        datos(k).env_lfp_max_protocolo_id = max_lfp_id;
        datos(k).env_lfp_norm_protocolo = (datos(k).env ./max_env_lfp) .*100;
    end
end

clear id_aux p t i file_name_lfp file_name_mua tetrodo dir_aux j id max_lfp_id max_mua_id
clear mua_aux env_aux i j mua_smooth mua_smooth2 lfp_aux k max_lfp max_mua
clear noise_aux mua_int max_mua_protocolo max_mua_protocolo_smooth max_env_lfp max_mua_id_smooth


%% Guardo todos los datos

datos_all = datos;


%% Genero sub-set

protocolo_analizar = 1; % indicar numero id del protocolo a analizar

% Separa automaticamente los datos del protocolo indicado
inicio_aux = (protocolo_analizar -1)*16 + 1;
fin_aux = inicio_aux + 15;
subset = (inicio_aux:1:fin_aux);
datos = datos_all(subset);

clear subset fin_aux inicio_aux


%% Corr de un canal de un protocolo contra resto de los canales

% Para un solo protocolo, busco el canal de mayor amplitud y calcula la
% corr de este canal contra el resto. Grafico las corr en una tabla que
% respeta la forma del NNx.

tet_max_ampl_mua = datos(1).mua_max_protocolo_id;
tet_max_ampl_lfp = datos(1).env_lfp_max_protocolo_id;

i_max_mua = find(contains({datos.id}, strcat(datos(1).protocolo, '_', tet_max_ampl_mua )));
i_max_lfp =  find(contains({datos.id}, strcat(datos(1).protocolo, '_', tet_max_ampl_lfp )));

% norm = input('\nCorrelacionamos seniales normalizadas por max del protocolo? (1 = SI, 0 = NO): ');
% norm = norm == 1;

corr_all_matrix_LFP = zeros(sqrt(length(datos)));
corr_all_matrix_MUA = zeros(sqrt(length(datos)));


% % Matriz donde guardo valor max de cada senial de cada tetrodo de todo el peine
% lfp_max = zeros(round(sqrt(length(datos))));
% mua_max = zeros(round(sqrt(length(datos))));
% 
% for i = (1:1:length(datos))
% 
%     lfp_max(i) = datos(i).max_env_lfp;
% 
%     mua_max(i) = datos(i).max_mua;
% 
% end
% 
% clear i


% OJO HARCODEO FEO
% Conservo seccion env lfp durante la presentacion del estimulo auditivo
fona_lfp = (1:1:length(datos(1).env))'; % vector de indices lfp
fona_lfp = fona_lfp < 4.5 * 1000; % transformo indices a t y me quedo con los menores a 4.5
% Conservo seccion mua durante la presentacion del estimulo
fona_mua = datos(1).mua(:,2) < 4.5;

% Calculo corr
count = 1;
for col = (1:1:4)

    for fila = (4:-1:1)
        
        % Calculo correlacion entre se�ales LFP
        corr_all_matrix_LFP(fila,col) = weighted_corr(datos(count).env_lfp_norm_protocolo(fona_lfp) , ...
            datos(i_max_lfp).env_lfp_norm_protocolo(fona_lfp), ...
            100);

        % Calculo correlacion entre se�ales MUA
       corr_all_matrix_MUA(fila,col) = weighted_corr(datos(count).mua_norm_protocolo(fona_mua,1), ...
            datos(i_max_mua).mua_norm_protocolo(fona_mua,1), ...
            100);
        
        count = count + 1;
    end
end


clear i j


labels_x = {'T4 | 600um', 'T3 | 400um', 'T2 | 200um', 'T1 | 0um'};
labels_y = {'P1 | 0um', 'P2 | 200um', 'P3 | 400um', 'P4 | 600um'};

% plot LFP
figure()
x = repmat(1:4,4,1); % generate x-coordinates
y = x'; % generate y-coordinates
% Generate Labels
t = num2cell(round(corr_all_matrix_LFP, 2)); % extact values into cells
t = cellfun(@num2str, t, 'UniformOutput', false); % convert to string

imagesc(corr_all_matrix_LFP)
colormap(gca,'parula');
colorbar();
text(x(:), y(:), t, 'HorizontalAlignment', 'Center') % Draw Image and Label Pixels
caxis([0,1]); % or [-1,1]
ticks = 1:1:length(datos);
set(gca,'TickLabelInterpreter','none')
set(gca, 'YTick', ticks, 'YTickLabel', labels_x);
set(gca, 'XTick', ticks, 'XTickLabel', labels_y, 'XTickLabelRotation',45);
axis equal
axis tight
title({datos(1).protocolo;'correlacion pesada LFP'}, 'Interpreter', 'None')


% plot MUA
figure()
x = repmat(1:4,4,1); % generate x-coordinates
y = x'; % generate y-coordinates
% Generate Labels
t = num2cell(round(corr_all_matrix_MUA,2)); % extact values into cells
t = cellfun(@num2str, t, 'UniformOutput', false); % convert to string

imagesc(corr_all_matrix_MUA)
colormap(gca,'parula');
colorbar();
text(x(:), y(:), t, 'HorizontalAlignment', 'Center') % Draw Image and Label Pixels
caxis([0,1]); % or [-1,1]
ticks = 1:1:length(datos);
set(gca,'TickLabelInterpreter','none')
set(gca, 'YTick', ticks, 'YTickLabel', labels_x);
set(gca, 'XTick', ticks, 'XTickLabel', labels_y, 'XTickLabelRotation',45);
axis equal
axis tight
title({datos(1).protocolo;'correlacion pesada MUA'}, 'Interpreter', 'None')


print_png(1, datos(1).dir_lfp, strcat('_corr_lfp_tet_max'))
print_png(2, datos(1).dir_mua, strcat('_corr_mua_tet_max'))



%% Corr Boari_2021

norm = input('\nCorrelacionamos se�ales normalizadas por max del protocolo? (1 = SI, 0 = NO): ');
norm = norm == 1;

corr_all_matrix_LFP = zeros(length(datos));
corr_all_matrix_MUA = zeros(length(datos));

% Matriz donde guardo valor max de cada se�al de cada tetrodo de todo el peine
lfp_max = zeros(round(sqrt(length(datos))));
mua_max = zeros(round(sqrt(length(datos))));

for i = (1:1:length(datos))

    lfp_max(i) = datos(i).max_env_lfp;

    mua_max(i) = datos(i).max_mua;

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

for i = (1:1:length(datos))

    % Hago lista de id de estimulos para usarlos como label en futuro grafico
    labels(1,i) = {datos(i).id};

    for j = (1:1:length(datos))
        % Calculo correlacion entre se�ales LFP
        if i == j
            corr_all_matrix_LFP(i,j) = 1;
        else
            if norm
                corr_all_matrix_LFP(i,j) = weighted_corr(datos(i).env_lfp_norm_protocolo(fona_lfp) , ...
                    datos(j).env_lfp_norm_protocolo(fona_lfp), ...
                    100);
            else
                corr_all_matrix_LFP(i,j) = weighted_corr(datos(i).env(fona_lfp), ...
                    datos(j).env(fona_lfp), ...
                    max(datos(i).env_lfp_max_protocolo, datos(j).env_lfp_max_protocolo));
            end
        end

        % Calculo correlacion entre se�ales MUA
        if i == j
            corr_all_matrix_MUA(i,j) = 1;
        else
            if norm
               corr_all_matrix_MUA(i,j) = weighted_corr(datos(i).mua_norm_protocolo(fona_mua,1), ...
                    datos(j).mua_norm_protocolo(fona_mua,1), ...
                    100);
            else
                corr_all_matrix_MUA(i,j) = weighted_corr(datos(i).mua(fona_mua,1), ...
                    datos(j).mua(fona_mua,1), ...
                    max(datos(i).mua_max_protocolo, datos(j).mua_max_protocolo));
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
caxis([0.5,1]); % or [-1,1]
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
caxis([0.5,1]); % or [-1,1]
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
title('MAX de tetrodos de MUA')

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
title('MAX de tetrodos de LFP')

clear ticks


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


%% CORR en funcion de la distancia entre tetrodos

close all

% Busco indice del tetrodo con mayor amplitud de MUA
index_max_mua = find([datos.max_mua] == datos(1).mua_max_protocolo);

col_corr_max_mua = corr_all_matrix_MUA(:,index_max_mua);

col_posiciones = zeros(length(datos),2);
for i = (1:length(datos))
    x = datos(i).coordinates(1);
    y = datos(i).coordinates(2);
    v = [x,y];
    col_posiciones(i,:) = v;
end 

col_distancias_max_mua = zeros(length(datos),1);
v1 = col_posiciones(index_max_mua,:);
for i = (1:length(datos))
    v2 = col_posiciones(i,:);
    d = distancia(v1,v2);
    col_distancias_max_mua(i) = d;
end 
% Ajustamos a modelo lineal con MCO
lm = fitlm(col_distancias_max_mua, col_corr_max_mua);
lm.Coefficients

p = plot(lm);
set(p,'LineWidth',5,'MarkerSize', 20);
xlabel('distancia al tetrodo max amplitud MUA (um)')
ylabel('corr al tetrodo max amplitud MUA')
title('dist vs corr al tetrodo de max amplitud MUA')
set(gca,'FontSize',30)

titulo = ['_dist_vs_corr_max-', datos(1).mua_max_protocolo_id, '_'];
dir_aux = datos(1).dir_mua;
print_pdf(1, dir_aux, strcat(titulo, '.pdf'))


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

% close all
% 
% % Elijo quienes plotear
% ploteo = [7, 28, 42];

% ploteo = [7, 6, 16];
% ploteo = [12, 11, 6];
ploteo = [7, 28];


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
count = 1;
for k = (1:1:length(ploteo))

    % Ploteo envolvente
    p = plot(datos(ploteo(k)).mua(:,2), datos(ploteo(k)).mua(:,1), 'LineWidth', lw);
    hold on

    leyendas{count,1} = datos(ploteo(k)).id;
    count = count +1;

    p.Color(4) = alpha;
end
legend(leyendas, 'Interpreter' , 'none');
title('MUA promediando tetrodo para el BOS')
set(gca,'FontSize',25)
xlim([0, dur_sound])

clear count


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
