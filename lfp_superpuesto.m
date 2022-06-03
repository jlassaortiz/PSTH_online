close all
clear all

directorio = '/Volumes/AUDIOS TEAM/Javi- finch durmierdo (NNx)/zf-JL037-VB/Datos/2021-12-21/zf-JL037-VB_p4_id1_211221_162306/';

lfp1 = readtable(strcat(directorio, 'LFP_1tet_BOS_P2-T2_BANDA-25-35Hz_-71uV.txt'));
lfp2 = readtable(strcat(directorio, 'LFP_1tet_BOS_P3-T2_BANDA-25-35Hz_-72uV.txt'));
lfp3 = readtable(strcat(directorio, 'LFP_1tet_BOS_P4-T1_BANDA-25-35Hz_-68uV.txt'));
lfp4 = readtable(strcat(directorio, 'LFP_1tet_BOS_P1-T2_BANDA-25-35Hz_-67uV.txt'));
lfp5 = readtable(strcat(directorio, 'LFP_1tet_BOS_P1-T3_BANDA-25-35Hz_-71uV.txt'));
lfp6 = readtable(strcat(directorio, 'LFP_1tet_BOS_P1-T1_BANDA-25-35Hz_-61uV.txt'));

lfp1 = lfp1{:,1};
lfp2 = lfp2{:,1};
lfp3 = lfp3{:,1};
lfp4 = lfp4{:,1};
lfp5 = lfp5{:,1};
lfp6 = lfp6{:,1};

env1 = abs(hilbert(lfp1));
env2 = abs(hilbert(lfp2));
env3 = abs(hilbert(lfp3));
env4 = abs(hilbert(lfp4));
env5 = abs(hilbert(lfp5));
env6 = abs(hilbert(lfp6));

alpha = 0.50;
lw = 2;


% LFP
figure()
p1 = plot(lfp1, 'LineWidth', lw);
hold on
p2 = plot(lfp2, 'LineWidth', lw);
hold on
p3 = plot(lfp3, 'LineWidth', lw);
% p4 = plot(lfp4, 'LineWidth', lw);
% p5 = plot(lfp5, 'LineWidth', lw);
% p6 = plot(lfp6, 'LineWidth', lw);

p1.Color(4) = alpha;
p2.Color(4) = alpha;
p3.Color(4) = alpha;
% p4.Color(4) = alpha;
% p5.Color(4) = alpha;

% legend('P2-T2', 'P3-T2', 'P4-T1', 'P1-T2', 'P1-T3');
legend('P2-T2', 'P3-T2', 'P4-T1');
title('LFP banda 25-35Hz promediando tetrodo para el BOS')
set(gca,'FontSize',25)


% LFP env
figure()
p11 = plot(env1, 'LineWidth', lw);
hold on
p22 = plot(env2, 'LineWidth', lw);
hold on 
p33 = plot(env3, 'LineWidth', lw);
% p44 = plot(env4, 'LineWidth', lw);
% p55 = plot(env5, 'LineWidth', lw);
% p66 = plot(env6, 'LineWidth', lw);

p11.Color(4) = alpha;
p22.Color(4) = alpha;
p33.Color(4) = alpha;
% p44.Color(4) = alpha;
% p55.Color(4) = alpha;

% legend('P3-T2','P1-T2','P1-T3','P1-T1');
legend('P2-T2', 'P3-T2', 'P4-T1');
title('LFP banda 25-35Hz promediando tetrodo para el BOS')
set(gca,'FontSize',25)


% LFP env normalizacion
figure()
p11 = plot(env1/max(env1), 'LineWidth', lw);
hold on
p22 = plot(env2/max(env2), 'LineWidth', lw);
hold on 
p33 = plot(env3/max(env3), 'LineWidth', lw);
% p44 = plot(env4/max(env4), 'LineWidth', lw);
% p55 = plot(env5/max(env5), 'LineWidth', lw);
% p66 = plot(env6/max(env6), 'LineWidth', lw);

p11.Color(4) = alpha;
p22.Color(4) = alpha;
p33.Color(4) = alpha;
% p44.Color(4) = alpha;
% p55.Color(4) = alpha;

% legend('P2-T2', 'P3-T2', 'P4-T1', 'P1-T2', 'P1-T3');
legend('P2-T2', 'P3-T2', 'P4-T1');
title('LFP banda 25-35Hz promediando tetrodo para el BOS')
set(gca,'FontSize',25)


% Guardo
titulo = '_LFP_promedioTetrodo_BOS_ZOOM';

print_pdf(1, directorio, strcat(titulo, '.pdf'))
print_pdf(2, directorio, strcat(titulo,'_env', '.pdf'))
print_pdf(3, directorio, strcat(titulo,'_env-norm', '.pdf'))

