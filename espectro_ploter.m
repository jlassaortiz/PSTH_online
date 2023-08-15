
directorio = 'C:\Users\jlassa\Desktop\029-RoN\ESTIMULOS-2\';
estimulos = carga_songs(directorio);

for i = (1:length(estimulos))
    plot_espectro(estimulos(i).name, estimulos(i).song, estimulos(i).freq)
end


for i = (1:length(estimulos))
    print_pdf(i, directorio, strcat('_', erase(estimulos(i).name, '.wav'),'.pdf'))
end

