function plot_some_sabana(score_total, mat_avg, ejeX_fila, ejeY_col)

%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here


% Ploteo integral
figure()

for i = (1:1:length(score_total))  
X  = score_total(i).grilla_scores(:,1);
Y  = score_total(i).grilla_scores(:,2);
Z1 = score_total(i).grilla_scores(:,3);
plot3(X,Y,Z1,'o')
legends{i} = score_total(i).id;
hold on
end

legends{length(legends) +1} = 'avg';

[xq,yq] = meshgrid(1:0.1:max(score_total(1).grilla_scores(:,1)));
z = griddata(mat_avg(:,1),mat_avg(:,2),mat_avg(:,3),xq,yq,'natural');
scatter3(mat_avg(:,1), mat_avg(:,2), mat_avg(:,3), 100, 'ro', 'filled')
mesh(xq,yq,z)
legend(legends, 'Interpreter','None')
ylabel(ejeY_col)
xlabel(ejeX_fila)
title('integral', 'Interpreter','None')

clear legends


% Ploteo correlacion
figure()

for i = (1:1:length(score_total))  
X  = score_total(i).grilla_scores(:,1);
Y  = score_total(i).grilla_scores(:,2);
Z2 = score_total(i).grilla_scores(:,4);

plot3(X,Y,Z2,'o')
legends{i} = score_total(i).id;
hold on
end

legends{length(legends) +1} = 'avg';

[xq,yq] = meshgrid(1:0.1:max(score_total(1).grilla_scores(:,1)));
z = griddata(mat_avg(:,1),mat_avg(:,2),mat_avg(:,4),xq,yq,'natural');
scatter3(mat_avg(:,1), mat_avg(:,2), mat_avg(:,4), 100, 'ro', 'filled')
mesh(xq,yq,z)
legend(legends, 'Interpreter','None')
ylabel(ejeY_col)
xlabel(ejeX_fila)
title('correlacion', 'Interpreter','None')

end

