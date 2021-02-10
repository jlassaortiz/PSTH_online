function print_png(fig_number, directorio, nombre_extra)

% Guarda la figura especificada en formato PDF

f = figure(fig_number);
nombre = strcat(directorio, datestr(now, 'yyyy-mm-dd_HH_MM_SS'), nombre_extra);
nombre = strcat(nombre, '.png');

% set(f,...
%     'Units','centimeters', ...
%     'PaperUnits', 'centimeters' , ...
%     'PaperSize',[60 46])

% print(nombre,f,'-dpng') % ,'-r300'  hace referencia a los dpi
exportgraphics(f,nombre, 'Resolution', 100)

end

