function [mat_scores, cell_estimulos] = scores_struct2mat(grilla,dict_score)

% Convierte salida de score_calculator en matriz para graficar facilmente.
%   
%   Entradas:
%
%   Salidas:
%   mat_scores = (array) en la columna 1 y 2 tiene los valores de X e Y de
%   la sabana. En las columnas 3 y 4 tiene los valores Z1 (integral) y Z2
%   (correlacion) de la sabana.
%   cell_estimulos = (cell) en la posicion cel_estimulos(x,y) tiene el
%   nombre del estimulo que esta en la posicion grilla(x,y)


% Inicializo la matriz donde van a ir los valores
% Esta pensado para analizar "SABANAS" es decir variaciones de DOS
% parametros del SYN y se le calcular DOS scores 
mat_scores = zeros(numel(grilla), 4);
cell_estimulos = cell(length(grilla));

% Incializo la fila
fila = 1;

% Recorro cada elemento de la grilla
[n_x, n_y] = size(grilla);
for x = (1:1:n_x)
    for y = (1:1:n_y)
        
        % Agrego valores X e Y de la sabana
        mat_scores(fila, 1) = x;
        mat_scores(fila, 2) = y;

        % Agrego valores Z de la sabana
        id_estimulo = grilla(x,y);
        mat_scores(fila, 3) = dict_score(id_estimulo).int_norm;
        mat_scores(fila, 4) = dict_score(id_estimulo).corr;

        % Guardo el nombre del estimulo (responsable de los scores de
        % arriba) en la misma posicion X  e Y que los scores
        cell_estimulos{x,y} = dict_score(id_estimulo).name;

        % Paso a la fila siguiente en la matriz
        fila = fila + 1;
    end
end
end

