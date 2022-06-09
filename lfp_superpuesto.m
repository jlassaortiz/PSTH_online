close all
clear all


%% Cargo los datos 

datos = struct();
i = 1;

% 1
dir1 = '/Volumes/AUDIOS TEAM/Javi- finch durmierdo (NNx)/zf-JL037-VB/Datos/2021-12-21/zf-JL037-VB_p4_id1_211221_162306/';
datos(i).id = '037-VB_p4_id1_211221_P2_T2';
datos(i).file_lfp = strcat(dir1, 'LFP_1tet_BOS_P2-T2_BANDA-25-35Hz_-71uV.txt');
datos(i).file_mua = strcat(dir1, 'PSTHsw_1tet_BOS_P2-T2_BANDA-25-35Hz_-71uV.txt');
i = i + 1; 

% 2
datos(i).id = '037-VB_p4_id1_211221_P3_T2';
datos(i).file_lfp = strcat(dir1, 'LFP_1tet_BOS_P3-T2_BANDA-25-35Hz_-72uV.txt');
datos(i).file_mua = strcat(dir1, 'PSTHsw_1tet_BOS_P3-T2_BANDA-25-35Hz_-72uV.txt');
i = i + 1; 

% 3
datos(i).id = '037-VB_p4_id1_211221_P4_T1';
datos(i).file_lfp = strcat(dir1, 'LFP_1tet_BOS_P4-T1_BANDA-25-35Hz_-68uV.txt');
datos(i).file_mua = strcat(dir1, 'PSTHsw_1tet_BOS_P4-T1_BANDA-25-35Hz_-68uV.txt');
i = i + 1; 

% 4
datos(i).id = '037-VB_p4_id1_211221_P1_T2';
datos(i).file_lfp = strcat(dir1, 'LFP_1tet_BOS_P1-T2_BANDA-25-35Hz_-67uV.txt');
datos(i).file_mua = strcat(dir1, 'PSTHsw_1tet_BOS_P1-T2_BANDA-25-35Hz_-67uV.txt');
i = i + 1; 

% 5
datos(i).id = '037-VB_p4_id1_211221_P1_T3';
datos(i).file_lfp = strcat(dir1, 'LFP_1tet_BOS_P1-T3_BANDA-25-35Hz_-71uV.txt');
datos(i).file_mua = strcat(dir1, 'PSTHsw_1tet_BOS_P1-T3_BANDA-25-35Hz_-71uV.txt');
i = i + 1; 

% 6
datos(i).id = '037-VB_p4_id1_211221_P1_T1';
datos(i).file_lfp = strcat(dir1, 'LFP_1tet_BOS_P1-T1_BANDA-25-35Hz_-61uV.txt');
datos(i).file_mua = strcat(dir1, 'PSTHsw_1tet_BOS_P1-T1_BANDA-25-35Hz_-61uV.txt');
i = i + 1; 

% 7
dir2 = '/Volumes/AUDIOS TEAM/Javi- finch durmierdo (NNx)/zf-JL037-VB/Datos/2021-12-20/zf-JL037-VB_p1_id1_211220_132913/';
datos(i).id = '037-VB_p1_id1_211220_P2_T3';
datos(i).file_lfp = strcat(dir2, 'LFP_1tet_BOS_P2-T3_BANDA-25-35Hz_-67uV.txt');
datos(i).file_mua = strcat(dir2, 'PSTHsw_1tet_BOS_P2-T3_BANDA-25-35Hz_-67uV.txt');
i = i + 1; 

% 8
dir3 = '/Volumes/AUDIOS TEAM/Javi- finch durmierdo (NNx)/zf-JL037-VB/Datos/2021-12-21/zf-JL037-VB_p2_id1_211221_140619/';
datos(i).id = '037-VB_p2_id1_211221_P3_T3';
datos(i).file_lfp = strcat(dir3, 'LFP_1tet_BOS_P3-T3_BANDA-25-35Hz_-79uV.txt');
datos(i).file_mua = strcat(dir3, 'PSTHsw_1tet_BOS_P3-T3_BANDA-25-35Hz_-79uV.txt');
i = i + 1; 

for j = (1:1:length(datos))
    
    % Levanto y guardo LFP
    lfp_aux = readtable(datos(j).file_lfp);
    lfp_aux = lfp_aux{:,1}; % paso de table a vector de double
    datos(j).lfp = lfp_aux;
    
    % Calculo y guardo envolvente LFP 30Hz
    env_aux = abs(hilbert(lfp_aux));
    datos(j).env = env_aux;
    
    % Calculo y guardo envolvente normalizada LFP 30Hz
    datos(j).env_norm = env_aux / max(env_aux);
    
    % Levanto y guardo MUA
    mua_aux = readtable(datos(j).file_mua);
    mua_aux = table2array(mua_aux);
    datos(j).mua = mua_aux;
    
    % Guardo MUA normalizada
    datos(j).mua_norm = mua_aux/max(mua_aux);
end 

clear mua_aux env_aux 

% SUBSET DE DATOS 
subset = [1,2,3,4,5,6];
datos = datos(subset);


%% Estrategias comparacion 


% Corr todos contra todos %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
corr_all_matrix = zeros(length(datos));

for i = (1:1:length(datos))
    for j = (1:1:length(datos))
        corr_all_matrix(i,j) = corr(datos(i).env_norm, datos(j).env_norm);  
    end  
end
% saco 1 de la diagonal
corr_all_matrix(corr_all_matrix > 0.99) = NaN;
corr_all_matrix

% Busco cual es la max corr y entre quienes se da
max_corr = max(corr_all_matrix ,[],'all')
[r, c] = find( corr_all_matrix == max_corr, 1)
max_corr_pair = {r, c ; datos(r).id, datos(c).id}

% Busco cual es la min corr y entre quienes se da
min_corr = min(corr_all_matrix,[],'all')
[r, c] = find( corr_all_matrix == min_corr, 1)
min_corr_pair = {r, c ;datos(r).id, datos(c).id}



% Diferencia absoluta media todos contra todos %%%%%%%%%%%%%%%%%%%%%%%%%%%%
diff_abs_media_all_matrix = zeros(length(datos));
n = length(datos(1).env_norm);

for  i = (1:1:length(datos))
    for j = (1:1:length(datos))
        env1_norm = datos(i).env_norm;
        env2_norm = datos(j).env_norm;
        diff_abs_media_all_matrix(i,j) = sum(abs(double(env2_norm) - double(env1_norm)))/n;
    end
end
% saco diagonal de 0
diff_abs_media_all_matrix(diff_abs_media_all_matrix == 0) = NaN;
diff_abs_media_all_matrix

% Busco cual es la max diff y veo entre quienes se da
max_diff = max(diff_abs_media_all_matrix,[],'all')
[r, c] = find( diff_abs_media_all_matrix == max_diff, 1)
max_dif_pair = {r, c ;datos(r).id, datos(c).id}


% Busco cual es el min diff y veo entre quienes se da
min_dif = min(diff_abs_media_all_matrix,[],'all')
[r, c] = find( diff_abs_media_all_matrix == min_dif, 1)
min_dif_pair = {r, c ;datos(r).id, datos(c).id}


% Picos %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%% PLOTEO LFP 30Hz NORMALIZDO

% Determino umbral (altura min) de picos en la envolvente normalizada (0-1)
thr = 0.5;

close all

% Elijo quienes plotear
ploteo = [8, 2];

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

figure()
leyendas = cell(length(ploteo),1);
for k = (1:1:length(ploteo))
    
    t = datos(ploteo(k)).mua(:,2);
    mua = datos(ploteo(k)).mua(:,1);
    p = plot(t , mua , 'LineWidth', lw);
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



