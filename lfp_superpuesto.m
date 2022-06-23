close all
clear all

%% Cargo los datos 

% Cargo directorios y nombre custom de protocolos a manopla
directorios = {{'/Volumes/AUDIOS TEAM/Javi- finch durmierdo (NNx)/zf-JL037-VB/Datos/2021-12-21/zf-JL037-VB_p4_id1_211221_162306/',...
    '037-VB_p4_id1_211221'},...
    {'/Volumes/AUDIOS TEAM/Javi- finch durmierdo (NNx)/zf-JL037-VB/Datos/2021-12-21/zf-JL037-VB_p2_id1_211221_140619/',...
    '037-VB_p2_id1_211221'}, ...
    {'/Volumes/AUDIOS TEAM/Javi- finch durmierdo (NNx)/zf-JL037-VB/Datos/2021-12-20/zf-JL037-VB_p1_id1_211220_132913/',...
    '037-VB_p1_id1_211220'}};

datos = struct();
i = 1;

% Cargo nombre_id protocolo, directorios, nombres archivos con lfp y mua
for j = (1:length(directorios))
    
    dir_aux = directorios{j}{1};
    id_aux = directorios{j}{2};
    for p = (1:4)
        for t = (1:4)
            tetrodo = strcat('P', num2str(p),'-', 'T', num2str(t));

            id = strcat(id_aux,'_', tetrodo);
            file_name_lfp = dir(strcat(dir_aux, 'LFP_1tet_BOS_', tetrodo, '_BANDA-25-35Hz_*.txt')).name;
            file_name_mua = dir(strcat(dir_aux, 'PSTHsw_1tet_BOS_', tetrodo, '_BANDA-25-35Hz_*.txt')).name;

            datos(i).id = id;
            datos(i).file_lfp = strcat(dir_aux, file_name_lfp);
            datos(i).file_mua = strcat(dir_aux, file_name_mua);
            i = i + 1;    
        end 
    end 
    
end

clear id_aux p t i file_name_lfp file_name_mua tetrodo dir_aux j id

for j = (1:1:length(datos))
    
    % Levanto y guardo LFP
    lfp_aux = readtable(datos(j).file_lfp);
    lfp_aux = lfp_aux{:,1}; % paso de table a vector de double
    datos(j).lfp = lfp_aux;
    
    % Calculo y guardo envolvente LFP 
    env_aux = abs(hilbert(lfp_aux));
    datos(j).env = env_aux;
    
    % Calculo y guardo envolvente normalizada LFP 
    datos(j).env_norm = env_aux / max(env_aux);
    
    % Levanto y guardo MUA
    mua_aux = readtable(datos(j).file_mua);
    mua_aux = table2array(mua_aux);
    datos(j).mua = mua_aux;
    
    % MUA SUAVIZADA
    mua_smooth = sgolayfilt(mua_aux(:,1), 4, 129);
    datos(j).mua_smooth = [mua_smooth, mua_aux(:,2)];

    mua_smooth2 = sgolayfilt(mua_aux(:,1), 4, 65);
    datos(j).mua_smooth2 = [mua_smooth2, mua_aux(:,2)];
        
    % Guardo MUA normalizada 
    datos(j).mua_norm = mua_aux;
    datos(j).mua_norm(:,1) = mua_aux(:,1)/max(mua_aux(:,1));
    
%     % Busco picos en envolvente LFP
%     thr = 0.5;
%     
%     % Encuentro picos con altura mayor a thr
%     [pks, locs] = findpeaks(env_aux);
%     pks_subset = pks > thr;
%     plot(locs(pks_subset,:), pks(pks_subset,:), 'or')
    
end 



clear mua_aux env_aux i j mua_smooth mua_smooth2 lfp_aux 

%% Guardo todos los datos

datos_all = datos;

%% Genero sub-set
subset = (33:1:48);
datos = datos_all(subset);

% lw = 1;
% close all
% 
% plot(datos(10).mua(:,2), datos(10).mua(:,1), 'LineWidth', lw)
% hold on 
% lw = 2;
% plot(datos(10).mua_smooth(:,2), datos(10).mua_smooth(:,1), ':', 'LineWidth', lw)
% plot(datos(10).mua_smooth2(:,2), datos(10).mua_smooth2(:,1), ':','LineWidth', lw)
% 
% xlim([0, 4.5])
% title(datos(10).id, 'Interpreter' , 'none')
% legend('normal', ' |  g:4, w:129', ' |  g:, w: 65')
% set(gca,'FontSize',25)

%% Corr Boari_2021

corr_all_matrix_LFP = zeros(length(datos));
corr_all_matrix_MUA = zeros(length(datos));

% Busco maximo de todas las señales (lo usare para normalizar)
lfp_max_all = 0;
mua_max_all = 0;

% Matriz donde guardo valor max de cada señal de cada tetrodo de todo el peine
lfp_max = zeros(4);
mua_max = zeros(4);

for i = (1:1:length(datos))
    
    % Conservo seccion env lfp durante la presentacion del estimulo auditivo
    fona_lfp = (1:1:length(datos(i).env))'; % vector de indices lfp
    fona_lfp = fona_lfp < 4.5 * 1000; % transformo indices a t y me quedo con los menores a 4.5
    
    max_aux_LFP = max(datos(i).env(fona_lfp));
    if max_aux_LFP > lfp_max_all
        lfp_max_all = max_aux_LFP;
    end
    
    lfp_max(i) = max_aux_LFP;
    
    
    % Conservo seccion mua durante la presentacion del estimulo auditivo
    fona_mua = datos(i).mua(:,2) < 4.5; % 4.5seg. es a ojo por lo que veo en graficos.
    
    max_aux_MUA = max(datos(i).mua(fona_mua,1));
    if max_aux_MUA > mua_max_all
        mua_max_all = max_aux_MUA;
    end
    
    mua_max(i) = max_aux_MUA;
    
end 

clear i rms_aux_LFP rms_aux_MUA 


for i = (1:1:length(datos))
    for j = (1:1:length(datos))
       
        corr_all_matrix_LFP(i,j) = weighted_corr(datos(i).env(fona_lfp), ...
            datos(j).env(fona_lfp), lfp_max_all);
        
        % Conservo seccion mua durante la presentacion del estimulo
        fona_mua = datos(i).mua(:,2) < 4.5;
    
        corr_all_matrix_MUA(i,j) = weighted_corr(datos(i).mua(fona_mua,1), ...
            datos(j).mua(fona_mua,1), mua_max_all);
    end
end


% Hago lista de id de estimulos para usarlos como label en futuro grafico
labels = cell(1,length(datos));
for i = (1:1:length(datos))
    labels(1,i) = {datos(i).id};
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


%%


% saco 1 de la diagonal
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
min_corr_pair_LFP = {r, c ;datos(r).id, datos(c).id}

% Busco cual es la min corr y entre quienes se da
min_corr_MUA = min(corr_all_matrix_MUA,[],'all')
[r, c] = find( corr_all_matrix_MUA == min_corr_MUA, 1)
min_corr_pair_MUA = {r, c ;datos(r).id, datos(c).id}


%% Diferencia absoluta media (MAD) todos contra todos 
MAD_all_matrix_LFP = zeros(length(datos));
MAD_all_matrix_MUA = zeros(length(datos));
n = length(datos(1).env_norm);

for  i = (1:1:length(datos))
    for j = (1:1:length(datos))
        env1_norm_LFP = datos(i).env_norm;
        env2_norm_LFP = datos(j).env_norm;
        MAD_all_matrix_LFP(i,j) = sum(abs(double(env2_norm_LFP) - double(env1_norm_LFP)))/n;
        
        norm1_MUA = datos(i).mua_norm(:,1);
        norm2_MUA = datos(j).mua_norm(:,1);
        MAD_all_matrix_MUA(i,j) = sum(abs(double(norm1_MUA) - double(norm2_MUA)))/n;
    end
end

% hago lista de id de estimulos para usarlos como label en futuro grafico
labels = cell(1,length(datos));
for i = (1:1:length(datos))
    labels(1,i) = {datos(i).id};
end 

% plot LFP 
figure()
imagesc(MAD_all_matrix_LFP)
colormap(gca,'parula');
colorbar();
%caxis([0,1]); % or [-1,1]
ticks = 1:1:length(datos); 
set(gca,'TickLabelInterpreter','none')
set(gca, 'YTick', ticks, 'YTickLabel', labels);
set(gca, 'XTick', ticks, 'XTickLabel', labels, 'XTickLabelRotation',45);
axis equal
axis tight
title('Mean Absolute Error LFP')

% plot MUA 
figure()
imagesc(MAD_all_matrix_MUA)
colormap(gca,'parula');
colorbar();
%caxis([0,1]); % or [-1,1]
ticks = 1:1:length(datos); 
set(gca,'TickLabelInterpreter','none')
set(gca, 'YTick', ticks, 'YTickLabel', labels);
set(gca, 'XTick', ticks, 'XTickLabel', labels, 'XTickLabelRotation',45);
axis equal
axis tight
title('Mean Abosulute Error MUA')


%%



% saco diagonal de 0
MAD_all_matrix_LFP(MAD_all_matrix_LFP == 0) = NaN;
MAD_all_matrix_LFP

MAD_all_matrix_MUA(MAD_all_matrix_MUA == 0) = NaN;
MAD_all_matrix_MUA

% Busco cual es la max diff y veo entre quienes se da
max_diff_LFP = max(MAD_all_matrix_LFP,[],'all')
[r, c] = find(MAD_all_matrix_LFP == max_diff_LFP, 1)
max_dif_pair_LFP = {r, c ;datos(r).id, datos(c).id}

max_diff_MUA = max(MAD_all_matrix_MUA,[],'all')
[r, c] = find(MAD_all_matrix_MUA == max_diff_MUA, 1)
max_dif_pair_MUA = {r, c ;datos(r).id, datos(c).id}

% Busco cual es el min diff y veo entre quienes se da
min_dif_LFP = min(MAD_all_matrix_LFP,[],'all')
[r, c] = find( MAD_all_matrix_LFP == min_dif_LFP, 1)
min_dif_pair_LFP = {r, c ;datos(r).id, datos(c).id}

min_dif_MUA = min(MAD_all_matrix_MUA,[],'all')
[r, c] = find( MAD_all_matrix_MUA == min_dif_MUA, 1)
min_dif_pair_MUA = {r, c ;datos(r).id, datos(c).id}


%% Picos 



%% PLOTEO LFP 30Hz NORMALIZDO

% Determino umbral (altura min) de picos en la envolvente normalizada (0-1)
thr = 0.5;

close all

% Elijo quienes plotear
ploteo = [8, 10, 11];

% Guardo figuras?
guardar = input('\nGuardo? (1 = SI / 0 = NO) : ');

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


% LFP env
figure()
for k = (1:1:length(ploteo))
    
    % Ploteo envolvente
    env_aux = datos(ploteo(k)).env;
    p= plot(datos(ploteo(k)).env, 'LineWidth', lw);
    hold on
    
    p.Color(4) = alpha;
end 
legend(leyendas, 'Interpreter' , 'none');
title('Envolvente LFP banda 25-35Hz promediando tetrodo para el BOS')
set(gca,'FontSize',25)


% LFP env normalizada y picos!
figure()
leyendas_env = cell(2*length(ploteo),1);
count = 1;
for k = (1:1:length(ploteo))
    
    % Ploteo env normalizada
    env_aux = datos(ploteo(k)).env_norm;
    p= plot(env_aux, 'LineWidth', lw);
    hold on
    
    % Encuentro picos con altura mayor a thr
    [pks, locs] = findpeaks(env_aux);
    pks_subset = pks > thr;
    plot(locs(pks_subset,:), pks(pks_subset,:), 'or')

    p.Color(4) = alpha;
    
    % Escribo leyendas
    leyendas_env{count,1} = datos(ploteo(k)).id;
    leyendas_env{count + 1,1} = strcat('pico > thr -',datos(ploteo(k)).id);
    count = count + 2; 
end 
legend(leyendas_env, 'Interpreter' , 'none');
title('Envolvente LFP NORMALIZADA banda 25-35Hz promediando tetrodo para el BOS')
set(gca,'FontSize',25)



%% PLOTEO MUA

% Determino umbral (altura min) de picos en la envolvente normalizada (0-1)
thr = 0.5;

% Elijo quienes plotear
ploteo = [9, 10, 11];

% Guardo figuras?
guardar = input('\nGuardo? (1 = SI / 0 = NO) : ');

% Defino alpha (transparencia) y LineWidth
alpha = 0.50;
lw = 2;

% MUA
figure()
leyendas = cell(length(ploteo),1);
for k = (1:1:length(ploteo))
    
    t = datos(ploteo(k)).mua_smooth(:,2);
    mua = datos(ploteo(k)).mua_smooth(:,1);
    p = plot(t , mua , 'LineWidth', lw);
    hold on

    p.Color(4) = alpha;

    leyendas{k,1} = datos(ploteo(k)).id;
end 
legend(leyendas, 'Interpreter' , 'none');
title('MUA banda 25-35Hz promediando tetrodo para el BOS')
set(gca,'FontSize',25)

% MUA norm
figure()
leyendas = cell(length(ploteo),1);
for k = (1:1:length(ploteo))
    
    t = datos(ploteo(k)).mua_norm(:,2);
    mua_norm = datos(ploteo(k)).mua_norm(:,1);
    p = plot(t , mua_norm , 'LineWidth', lw);
    hold on

    p.Color(4) = alpha;

    leyendas{k,1} = datos(ploteo(k)).id;
end 
legend(leyendas, 'Interpreter' , 'none');
title('MUA banda 25-35Hz promediando tetrodo para el BOS')
set(gca,'FontSize',25)


%% Guardo
if guardar
    titulo = '_LFP_promedioTetrodo_BOS';

    print_pdf(1, dir1, strcat(titulo, '.pdf'))
    print_pdf(2, dir1, strcat(titulo,'_env.pdf'))
    print_pdf(3, dir1, strcat(titulo,'_env-norm.pdf'))
    print_pdf(4, dir1, strcat(titulo, '_mua.pdf'))
end 



