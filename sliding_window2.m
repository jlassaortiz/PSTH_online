function [sliding_window_data, sliding_window_tiempo] = sliding_window2( ...
    data, fs, t_window, step)

% Toma un vector unidimensional y le calcula la slinding window (SW)
% con los parametros de tamaño de ventana (en seg.) y tamaño del paso  
% especificados(en seg.). Además necesita que se especifique la frecuencia 
% de sampleo del vector unidimensional.
%
% data = vector unidimensional al que se le va a calcular la SW
% fs = frecuencia de sampleo de 'data'
% t_window = tamaño de la ventana de la SW (en seg.)
% step = tamaño del paso de corrimiento de la SW (en seg.)

% Calculo el tamaño final de: sliding_window_data y sliding_window_tiempo
tf = t_window;
n = 0;
limite = length(data)/fs;

while tf <= limite
   n = n + 1;
   tf = tf + step;
end

% Inicializo vector donde se guardan resultados y tiempos
sliding_window_data = zeros(n, 1);
sliding_window_tiempo = zeros(n,1);
tf = t_window;
ti = 0;

% Calculo el valor de cada punto de la sw y el tiempo que le corresponde
for i = (1:1:n)

    punto = mean(data( uint32(ti*fs)+1 : uint32(tf*fs) , 1));
    sliding_window_data(i) = punto;

    sliding_window_tiempo(i) = (ti + tf)/2;

    ti = ti + step;
    tf = tf + step;
end
    
end