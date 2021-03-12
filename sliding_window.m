function [sliding_window_data, sliding_window_tiempo] = sliding_window(data, fs, t_window, step)

% Toma un vector unidimensional y le calcula la slinding window (SW)
% con los parametros de tama�o de ventana (en seg.) y tama�o del paso  
% especificados(en seg.). Adem�s necesita que se especifique la frecuencia 
% de sampleo del vector unidimensional.
%
% data = vector unidimensional al que se le va a calcular la SW
% fs = frecuencia de sampleo de 'data'
% t_window = tama�o de la ventana de la SW (en seg.)
% step = tama�o del paso de corrimiento de la SW (en seg.)


% Calculo el tama�o final de: sliding_window_data y sliding_window_tiempo
tf = t_window;
n = 0;
limite = max(data)/fs;

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

    punto = data( data > ti*fs & data <= tf*fs );
    sliding_window_data(i) = numel(punto);

    sliding_window_tiempo(i) = (ti + tf)/2;

    ti = ti + step;
    tf = tf + step;
end
    
end
