function wc = weighted_corr(signal_1,signal_2, A)
% Modified version of weighted correlation coefficient from Boari 2021
    % signal_1 and signal_2: (matrix nx1)
    
    r = corr(signal_1, signal_2);
    max_1 = max(signal_1);
    max_2 = max(signal_2);
    
    wc = r*(max_1 * max_2)/(A*A);
end