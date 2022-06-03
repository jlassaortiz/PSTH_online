close all
clear all

directorio = '/Volumes/AUDIOS TEAM/Javi- finch durmierdo (NNx)/zf-JL037-VB/Datos/2021-12-21/zf-JL037-VB_p4_id1_211221_162306/';

psth1 = readtable(strcat(directorio, 'PSTHsw_1tet_BOS_P2-T2_BANDA-25-35Hz_-71uV.txt'));
psth2 = readtable(strcat(directorio, 'PSTHsw_1tet_BOS_P3-T2_BANDA-25-35Hz_-72uV.txt'));
psth3 = readtable(strcat(directorio, 'PSTHsw_1tet_BOS_P4-T1_BANDA-25-35Hz_-68uV.txt'));
psth4 = readtable(strcat(directorio, 'PSTHsw_1tet_BOS_P1-T2_BANDA-25-35Hz_-67uV.txt'));
psth5 = readtable(strcat(directorio, 'PSTHsw_1tet_BOS_P1-T3_BANDA-25-35Hz_-71uV.txt'));
psth6 = readtable(strcat(directorio, 'PSTHsw_1tet_BOS_P1-T1_BANDA-25-35Hz_-61uV.txt'));

psth1 = psth1{:,1};
psth2 = psth2{:,1};
psth3 = psth3{:,1};
psth4 = psth4{:,1};
psth5 = psth5{:,1};
psth6 = psth6{:,1};

alpha = 0.5;
lw = 2;


% Ploteo PSTH superpuestos
figure()
% p1 = plot(psth1, 'LineWidth', lw);
% hold on
p2 = plot(psth2, 'LineWidth', lw);
hold on
% p3 = plot(psth3, 'LineWidth', lw);
% p4 = plot(psth4, 'LineWidth', lw);
% p5 = plot(psth5, 'LineWidth', lw);
p6 = plot(psth6, 'LineWidth', lw);

% p1.Color(4) = alpha;
p2.Color(4) = alpha;
% p3.Color(4) = alpha;
% p4.Color(4) = alpha;
% p5.Color(4) = alpha;
p6.Color(4) = alpha; 

% legend('P2-T2', 'P3-T2', 'P4-T1', 'P1-T2', 'P1-T3');
legend('P3-T2', 'P1-T1');
title('PSTH promediando tetrodo para el BOS')
set(gca,'FontSize',25)


% Ploteo PSTH superpuestos normalizados
figure()
% p11 = plot(psth1/max(psth1), 'LineWidth', lw);
% hold on
p22 = plot(psth2/max(psth2), 'LineWidth', lw);
hold on 
% p33 = plot(psth3/max(psth3), 'LineWidth', lw);
% p44 = plot(psth4/max(psth4), 'LineWidth', lw);
% p55 = plot(psth5/max(psth5), 'LineWidth', lw);
p66 = plot(psth6/max(psth6), 'LineWidth', lw);

% p11.Color(4) = alpha;
p22.Color(4) = alpha;
% p33.Color(4) = alpha;
% p44.Color(4) = alpha;
% p55.Color(4) = alpha;
p66.Color(4) = alpha; 

% legend('P2-T2', 'P3-T2', 'P4-T1', 'P1-T2', 'P1-T3');
legend('P3-T2', 'P1-T1');
title('PSTH normalizado promediando tetrodo para el BOS')
set(gca,'FontSize',25)

% Guardo
titulo = '_PSTH_promedioTetrodo_BOS';

print_pdf(1, directorio, strcat(titulo, '.pdf'))
print_pdf(2, directorio, strcat(titulo,'_normalizado', '.pdf'))
