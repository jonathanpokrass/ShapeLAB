function [F] = calcHKS(evecs1, evals1, numSamples)
        tn = linspace(65,90, numSamples);% 1783 2353 3104]; % full
        F =  (evecs1.^2)*exp(-bsxfun(@times, abs(evals1(:)), tn(:)'));    
        ;%F = F';
%        F  = normalize(F, 'l1', 2);
end