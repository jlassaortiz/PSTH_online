

function plotSimple_song_psth_lfp(plotear)
% Plotea varios (Canto + PSTH_sw + LFP) en grilla 3x:
%   Necesita de entrada un solo struct de 4xN (plotear) estructurado asi:
%   plotear(i).subTitle = (str) subtitulo de plot i
%   plotear(i).song = (mat :x2) col1 = song, col2 = tiempos (s) de plot i
%   plotear(i).psth = (mat :x2) col1 = sliding window psth, col2 = t (s)
%   plotear(i).lfp = (mat :x2) col1 = lfp, col2 = t (s) de plot i

% Numero de plots a graficar en la grilla
n_plots = length(plotear);

% Inicializo figura
figure()

% Formula para armar grilla segun la cantidad de estimulos a analizar
if mod(n_plots, 3) == 0
    n = 5 * n_plots/3;
    
elseif mod(n_plots, 3) == 2
    n = 5 * round(n_plots/3);
    
else
    n = 5 * (round(n_plots/3) + 1);
end
    
m = 3;    
j = 0;
k = 0;

% Para caso 2 graficos
pos2 = 1;
graf2 = 1; 

% Para caso 4 graficos
pos4 = 1;
graf4 = 1;

% Busco maximo del PSTH para determinar ylim
psth_max = 0;
for count = (1:n_plots)
    a = max(plotear(count).psth(:,1));
    if a > psth_max
        psth_max = a;
    end 
end

% Busco maximo y min del LFP para determinar ylim
lfp_max = 0;
lfp_min = 0;
for count = (1:n_plots)
    a = max(plotear(count).lfp(:,1));
    b = min(plotear(count).lfp(:,1));
    if a > lfp_max
        lfp_max = a;
    elseif b < lfp_min
        lfp_min = b; 
    end
end
    

% Para cada estimulo/plot
for i = (1:n_plots)
    
    k = k + 1;
    
    if mod(k, 3) == 1
        p = (k - 1)/3 * 15 + 1;
        
    elseif mod(k, 3) == 2
        p = (k-2)/3 * 15 + 2;
        
    else
        p = ((k / 3) - 1) * 15 + 3;
    end
     
    
    % SONIDO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    j = j + 1;
    if n_plots == 1
        h(1) = subplot(5,1,1);
    elseif n_plots == 2
        h(graf2) = subplot(10,1,pos2);
    elseif n_plots == 4
        h(graf4) = subplot(10,2, pos4);
    else
        h(j) = subplot(n, m , p);
    end
    
    plot(plotear(i).song(:,2), plotear(i).song(:,1),'black')
    xticks([])
    title(plotear(i).subTitle, 'Interpreter','None','FontSize', 6)

    
    % PSTH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    j = j + 1;
    if n_plots == 1
        h(2) = subplot(5,1, [2,3]);
    elseif n_plots == 2
        h(graf2 + 1) = subplot(10,1, [pos2 + 1,pos2 + 2]);
    elseif n_plots == 4
        h(graf4 + 1) = subplot(10, 2, [pos4 + 2, pos4 + 4]);
    else
        h(j) = subplot(n, m, [p + 3, p + 6]);
    end
    
    plot(plotear(i).psth(:,2), plotear(i).psth(:,1), '-r', 'LineWidth', 2);
    ylim([0 psth_max]);
  
    
    % LFP promediado %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    j = j + 1;
    if n_plots == 1
        h(3) = subplot(5,1,[4,5]);
    elseif n_plots == 2
        h(graf2 + 2) = subplot(10,1, [pos2 + 3, pos2 + 4]);
    elseif n_plots == 4
        h(graf4 + 2) = subplot(10, 2, [pos4 + 6, pos4 + 8]);
    else 
        h(j) = subplot(n, m, [p + 9, p + 12]);
    end
    
    % Avance para caso 2 graf
    pos2 = pos2 + 5;
    graf2 = graf2 + 3;
    
    % Avance para el caso 4 graf
    if i == 1
        graf4 = 4;
        pos4 = 2;
    elseif i == 2
        graf4 = 7;
        pos4 = 11;
    elseif i == 3
        graf4 = 10;
        pos4 = 12;
    end
    
    plot(plotear(i).lfp(:,2), plotear(i).lfp(:,1), '-b', 'LineWidth', 2)
    ylim([lfp_min lfp_max]);
    xticks([])
end 

% Linkeo eje x
linkaxes(h, 'x');

end