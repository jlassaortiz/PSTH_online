function tabla_all_datos = score_total_toTable(score_total, directorio_params)
% Toma el objeto score_total (sale de score_all_and_plot)
%   
%   IN:
%   score_total (tabla)
%
%   OUT:
%   tabla_all_datos (tabla "cuadrada")

protocolo = (1:length(score_total));
all_lambdas = [0.5, 1.0, 1.5, 2.0, 2.5];

id = [];
lambda = [];
integral = [];
lfp_dif = [];

for p = protocolo
    count = 1;
    for l = all_lambdas
        id = vertcat(id, convertCharsToStrings(score_total(p).id));
        lambda = vertcat(lambda, l);
        
        index_int = score_total(p).grilla_scores(:,1) == count;
        int = unique(score_total(p).grilla_scores(index_int, 3));
        
        index_lfp = score_total(p).grilla_scores(:,1) == count;
        lfp_d = unique(score_total(p).grilla_scores(index_lfp, 4));
        
        integral = vertcat(integral, int);
        lfp_dif = vertcat(lfp_dif, lfp_d);
        count = count + 1;
    end   
end 

tabla_all_datos = table(id, lambda, integral, lfp_dif);
tabla_all_datos.integral = round(tabla_all_datos.integral, 3);
tabla_all_datos.lfp_dif = round(tabla_all_datos.lfp_dif, 3);

time_stamp = convertCharsToStrings(datestr(now, "yyyy-mm-dd_HH_MM_SS"));

writetable(tabla_all_datos, directorio_params + ...
     time_stamp + "_tabla_all_datos.csv")

end

