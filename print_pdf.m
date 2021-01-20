function print_pdf(fig_number, nombre, directorio)

% Guarda la figura especificada en formato PDF

f = figure(fig_number);
nombre = strcat(directorio, datestr(now, 'yyyy-mm-dd_HH_MM_SS')," ", nombre);

set(f,...
    'Units','centimeters', ...
    'PaperUnits', 'centimeters' , ...
    'PaperSize',[60 46])

print(nombre,f,'-dpdf','-bestfit','-painters') % ,'-r300'  hace referencia a los dpi

end

