function print_png(fig_number, directorio, nombre_extra)

% Guarda la figura especificada en formato PDF

f = figure(fig_number);
nombre = strcat(directorio, datestr(now, 'yyyy-mm-dd_HH_MM_SS'), nombre_extra);
nombre = strcat(nombre, '.png');

print(f, nombre,'-dpng','-r100')

end

