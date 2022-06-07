close all
clear all


%% Cargo los datos 

datos = struct();
i = 1;

% 1
dir1 = '/Volumes/AUDIOS TEAM/Javi- finch durmierdo (NNx)/zf-JL037-VB/Datos/2021-12-21/zf-JL037-VB_p4_id1_211221_162306/';
datos(i).id = '037-VB_p4_id1_211221_P2_T2';
datos(i).dir = strcat(dir1, 'LFP_1tet_BOS_P2-T2_BANDA-25-35Hz_-71uV.txt');
i = i + 1; 

% 2
datos(i).id = '037-VB_p4_id1_211221_P3_T2';
datos(i).dir = strcat(dir1, 'LFP_1tet_BOS_P3-T2_BANDA-25-35Hz_-72uV.txt');
i = i + 1; 

% 3
datos(i).id = '037-VB_p4_id1_211221_P4_T1';
datos(i).dir = strcat(dir1, 'LFP_1tet_BOS_P4-T1_BANDA-25-35Hz_-68uV.txt');
i = i + 1; 

% 4
datos(i).id = '037-VB_p4_id1_211221_P1_T2';
datos(i).dir = strcat(dir1, 'LFP_1tet_BOS_P1-T2_BANDA-25-35Hz_-67uV.txt');
i = i + 1; 

% 5
datos(i).id = '037-VB_p4_id1_211221_P1_T3';
datos(i).dir = strcat(dir1, 'LFP_1tet_BOS_P1-T3_BANDA-25-35Hz_-71uV.txt');
i = i + 1; 

% 6
datos(i).id = '037-VB_p4_id1_211221_P1_T1';
datos(i).dir = strcat(dir1, 'LFP_1tet_BOS_P1-T1_BANDA-25-35Hz_-61uV.txt');
i = i + 1; 

% 7
dir2 = '/Volumes/AUDIOS TEAM/Javi- finch durmierdo (NNx)/zf-JL037-VB/Datos/2021-12-20/zf-JL037-VB_p1_id1_211220_132913/';
datos(i).id = '037-VB_p1_id1_211220_P2_T3';
datos(i).dir = strcat(dir2, 'LFP_1tet_BOS_P2-T3_BANDA-25-35Hz_-67uV.txt');
i = i + 1; 

% 8
dir3 = '/Volumes/AUDIOS TEAM/Javi- finch durmierdo (NNx)/zf-JL037-VB/Datos/2021-12-21/zf-JL037-VB_p2_id1_211221_140619/';
datos(i).id = '037-VB_p2_id1_211221_P3_T3';
datos(i).dir = strcat(dir3, 'LFP_1tet_BOS_P3-T3_BANDA-25-35Hz_-79uV.txt');
i = i + 1; 

ploteo = [2, 7 , 8];

for j = (1:1:length(datos))
    
    % Levanto y guardo LFP
    lfp_aux = readtable(datos(j).dir);
    lfp_aux = lfp_aux{:,1}; % paso de table a vector de double
    datos(j).lfp = lfp_aux;
    
    % Calculo y guardo envolvente
    env_aux = abs(hilbert(lfp_aux));
    datos(j).env = env_aux;
    
end 


%% PLOTEO

% Defino alpha (transparencia) y LineWidth
alpha = 0.50;
lw = 2;

% Determino umbral (altura min) de picos en la envolvente normalizada (0-1)
thr = 0.5;


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


% LFP env normalizada
figure()
leyendas_env = cell(2*length(ploteo),1);
count = 1;
for k = (1:1:length(ploteo))
    
    % Calculo env normalizada y ploteo
    env_aux = datos(ploteo(k)).env;
    env_aux = env_aux/max(env_aux);
    p= plot(env_aux, 'LineWidth', lw);
    hold on
    
    % Encuentro picos con altura mayor a thr
    [pks, locs] = findpeaks(env_aux);
    pks_subset = pks > thr;
    plot(locs(pks_subset,:), pks(pks_subset,:), 'or')

    p.Color(4) = alpha;
    
    % Escribo leyendas
    leyendas_env{count,1} = datos(ploteo(k)).id;
    leyendas_env{count + 1,1} = strcat('pico > thr-',datos(ploteo(k)).id);
    count = count + 2; 
    end 
legend(leyendas_env, 'Interpreter' , 'none');
title('Envolvente LFP NORMALIZADA banda 25-35Hz promediando tetrodo para el BOS')
set(gca,'FontSize',25)


% Guardo
titulo = '_LFP_promedioTetrodo_BOS';

print_pdf(1, dir1, strcat(titulo, '.pdf'))
print_pdf(2, dir1, strcat(titulo,'_env', '.pdf'))
print_pdf(3, dir1, strcat(titulo,'_env-norm', '.pdf'))




%% Estrategias comparacion 
env2_n = env2/max(env2);
env1_n = env1/max(env1);
% env1_n = env7/max(env7);

% corr
corr(env2_n, env1_n)

% cross corr
figure()
plot(xcorr(env2_n, env1_n))

% suma diferencia cuadratica media
sum(abs(env2_n - env1_n))

% find peaks
[pks2, locs2] = findpeaks(env2_n);
[pks1, locs1] = findpeaks(env1_n);

thr2 = pks2 > 0.5;
thr1 = pks1 > 0.5;

figure()
plot(env1_n)
hold on
plot(env2_n)
plot(locs1(thr1,:), pks1(thr1,:), 'ob')
plot(locs2(thr2,:), pks2(thr2,:), 'or')






