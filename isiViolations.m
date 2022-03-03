function [violationFraction, numViolations]  = isiViolations(spikeTrain, refDur)
% Calcula fraccion de neuronas que no respentan el tiempo refractario
%
% INPUTS
%   spikeTrain: (mat 1xn) [s] spike times a evaluar 
%   refDur: (num)   [s] tiempo refractario 
%
% OUTPUTS
%   violationFRaction: (num) [%] fraccion spikes q violan el t refractario
%   numViolations: (num) [un.] cantidad spikes q violan el t refractario

numViolations = sum(diff(spikeTrain) <= refDur)

totalSpikes = length(spikeTrain)

violationFraction = numViolations / totalSpikes

if violationFraction < 0.01
    disp('GOOD spike train!')  
else
    disp('BAD spike train!')
end 
    
end