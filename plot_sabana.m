function plot_sabana(mat_scores, directorio, ejeY_col, ejeX_fila)

% Plotea Sabana x2 (integral y correlacio) y perfiles de sabana x2
%
%   En total hace 4 ploteos:
%     - sabana usando integral como variable respuesta
%     - sabana usando r-pearson como variable respuesta
%     - perfiles de las sabanas (x2)
   

X  = mat_scores(:,1);
Y  = mat_scores(:,2);
Z1 = mat_scores(:,3);
Z2 = mat_scores(:,4);

% Plot sabana INTEGRAL
[xq,yq] = meshgrid(1:0.1:3);
z = griddata(X,Y,Z1,xq,yq,'natural');
figure()
plot3(X,Y,Z1,'mo')
hold on
mesh(xq,yq,z)
ylabel(ejeY_col)
xlabel(ejeX_fila)
title({'integral', directorio}, 'Interpreter','None')

% Plot sabana CORRELACION
[xq,yq] = meshgrid(1:0.1:3);
z = griddata(X,Y,Z2,xq,yq,'natural');
figure()
plot3(X,Y,Z2,'mo')
hold on
mesh(xq,yq,z)
ylabel(ejeY_col)
xlabel(ejeX_fila)
title({'correlacion', directorio}, 'Interpreter','None')



% Hago diagonales para plotear perfiles
diag_escala = mat_scores(:,1) == 1 & mat_scores(:,2) == 1 | ... 
    mat_scores(:,1) == 2 & mat_scores(:,2) == 2 | ...
    mat_scores(:,1) == 3 & mat_scores(:,2) == 3 ;

diag_mostruo = mat_scores(:,1) == 1 & mat_scores(:,2) == 3 | ... 
    mat_scores(:,1) == 2 & mat_scores(:,2) == 2 | ...
    mat_scores(:,1) == 3 & mat_scores(:,2) == 1 ;

% Ploteo perfiles INTEGRAL
figure()
plot([0, 1, 2], mat_scores(diag_escala, 3))
hold on
plot([0, 1, 2], mat_scores(diag_mostruo, 3))
legend('escala', 'mostruo')
ylabel('integral')
title({'integral', directorio}, 'Interpreter','None')

% Ploteo perfiles CORRELACION
figure()
plot([0, 1, 2], mat_scores(diag_escala, 4))
hold on
plot([0, 1, 2], mat_scores(diag_mostruo, 4))
legend('escala', 'mostruo')
ylabel('correlacion')
title({'correlacion', directorio}, 'Interpreter','None')

end