function thr = find_thr(raw_filtered, t0s_dictionary, tiempo_file, frequency_parameters)

% Determina umbral para detectar actividad neuronal "multi-unit" (MUA)
%
%   Utilizo la forma de determinar umbral de quiam quiroga en momenotos de
%   silencio (no se presentan estimulos auditivos). Esta decisión se basa
%   en que durante la presentación de estímulos la tasa de disparo de las
%   neuronas aumenta respecto a los momentos de silencio. Para determinar
%   menor el desvío estandar del ruido es mejor hacerlo en momentos donde
%   la actividad neuronal es baja. De esta manera se evita sobre estimar el
%   desvío estandar del ruido. 
%   (Quiroga, R. Q., Nadasdy, Z., & Ben-Shaul, Y. (2004). Unsupervised Spike Detection and Sorting with Wavelets and Superparamagnetic Clustering. Neural Computation, 16(8), 1661?1687. https://doi.org/10.1162/089976604774201631

% Busco t0s pseudo-aleatorios (en sampling rate no en unidad de tiempo) y 
% lo transformo de tal manera que este al mismo sr que el raw data
t1 = t0s_dictionary(1).t0s(1,1) / frequency_parameters.board_adc_sample_rate * frequency_parameters.amplifier_sample_rate;
t2 = t0s_dictionary(1).t0s(2,1) / frequency_parameters.board_adc_sample_rate * frequency_parameters.amplifier_sample_rate;
t3 = t0s_dictionary(1).t0s(3,1) / frequency_parameters.board_adc_sample_rate * frequency_parameters.amplifier_sample_rate;

% Transformo el tiempo entre estimulos (en seg.) a sampling rate
tiempo_file_sample = tiempo_file * frequency_parameters.amplifier_sample_rate;

% Calculo duracion de usual de BOS (2 seg.) en sampling rate
duracion_BOS = 2 * frequency_parameters.amplifier_sample_rate;

% Extraigo 3 segmentos de raw data filtrada para calcular el std del ruido
% Extraigo segmentos donde hay silencio (no se presentan estimulos)
noise1 = raw_filtered(int64(t1 + tiempo_file_sample - duracion_BOS) : int64(t1 + tiempo_file_sample));
noise2 = raw_filtered(int64(t2 + tiempo_file_sample - duracion_BOS) : int64(t2 + tiempo_file_sample));
noise3 = raw_filtered(int64(t3 + tiempo_file_sample - duracion_BOS) : int64(t3 + tiempo_file_sample));

% Calculo el desvio estandar del ruido (Quiroga 2004) para tres segmentos 
std1 = median( abs(noise1)/0.6745 );
std2 = median( abs(noise2)/0.6745 );
std3 = median( abs(noise3)/0.6745 );

% Promedio los tres desvios para sacar efectos de outliers
stdmean = (std1 + std2 + std3) / 3;

% Determino el umbral a partir desvio estandar del ruido (CASI Quiroga 2004)
thr = - 5 * stdmean;

% figure()
% histogram(abs(noise1))
% hold on
% xline(median( abs(noise1)/0.6745 ))
% xline( 4 * median( abs(noise1)/0.6745 ), ':r')
% 
% figure()
% histogram(noise1)
% hold on
% xline(std(noise1))
% xline(- 4 * std(noise1), ':r')

% histogram(noise2)
% histogram(noise3)

end

