function plot_sabana(mat_scores, directorio, ejeY_col, ejeX_fila)

%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

X  = mat_scores(:,1);
Y  = mat_scores(:,2);
Z1 = mat_scores(:,3);
Z2 = mat_scores(:,4);

[xq,yq] = meshgrid(1:0.1:3);
z = griddata(X,Y,Z1,xq,yq,'natural');
figure()
plot3(X,Y,Z1,'mo')
hold on
mesh(xq,yq,z)
ylabel(ejeY_col)
xlabel(ejeX_fila)
title({'integral', directorio})

[xq,yq] = meshgrid(1:0.1:3);
z = griddata(X,Y,Z2,xq,yq,'natural');
figure()
plot3(X,Y,Z2,'mo')
hold on
mesh(xq,yq,z)
ylabel(ejeY_col)
xlabel(ejeX_fila)
title({'correlacion', directorio})

end

