function lim_bandas = list_multiple_bands(banda_sup, p, n, ploteo)
% Genera lista con n bandas superpuestas
%
%  Determinamos ancho de banda a dividir y cantidad de divisiones, el
%  programa devuelve los limites de las bandas superpuestas de igual lago
%  en el espacio de logaritmico de frecuencias
%
%  Input: banda_inf y banda_sup en Hz
%  p: proporcion de superposicion (0 a 1)
%  n: cantidad de bandas a obtener (numero)
%  ploteo: si ploteo o no (True or False)

% Defino liminte de frecuencia a analizar en multiples bandas
k = log10(banda_sup);

% Relacion entre a y b para determinar bordes (ver cuaderno labo)
a = k/(n + (n*p/(1-2*p)));
b = p*a/(1-2*p);

% Inicializo lista de bandas y determino la primera
% primer columna es comienzo y segunda es fin de la banda
% tercer columna es para graficar
lim_bandas = zeros(n,3);
lim_bandas(1,:) = [0, b+a+b, 1]; 

% Calculo el resto de los limites de las bandas
for i = (2:n)
    lim_bandas(i,:) = [lim_bandas(i-1,2)-b, lim_bandas(i-1,2)+a+b, i];
end

% Ploteo para ver las bandas
if ploteo
    for i = (1:n)
        plot(lim_bandas(i,1:2),[i,i], '-')
        hold on
    end
    xline(k, 'r');
    xlabel('log10(Hz)')
    ylabel('# banda')
end

% Saco la columna auxiliar 3 que solo era para graficar.
lim_bandas = lim_bandas(:,1:2);

% Devuelvo los limites a Hz (estan en log)
lim_bandas = 10.^lim_bandas;

end

