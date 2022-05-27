function LFP_filt_norm = filt_and_normalize(LFP,lower_filter_bound, upper_filter_bound, sr)

% filt_and_normalize() filtra por banda y normaliza por amplitud media
%
% Algoritmo filtro pasabandas tomado de "Analyzing Neural Time Series Data" 
% de Mike X Cohen

% Armo filtro con pasabandas firls
nyquist          = sr/2;
transition_width = 0.2; 
filter_order     = round(3*(sr/lower_filter_bound));

% create the filter shape (this is explained more in the text around figure 14.4)
ffrequencies  = [ 0 (1-transition_width)*lower_filter_bound lower_filter_bound upper_filter_bound (1+transition_width)*upper_filter_bound nyquist ]/nyquist;
idealresponse = [ 0 0 1 1 0 0 ];
filterweights = firls(filter_order,ffrequencies,idealresponse);

% aplico filtro
LFP_filt = filtfilt(filterweights,1,LFP);

% calculo transformada de Hilbert
h = hilbert(LFP_filt);

% para evitar efecto borde, calculo media hilbert recortando inicio y fin señal
inicio = uint32( length(LFP_filt) * 0.1 );
fin = uint32( length(LFP_filt) * 0.9 );
m = median( abs( h(inicio:fin,1) ) );

% normalizo por la media de la amplitud de la banda
LFP_filt_norm = LFP_filt/m; 
 end

