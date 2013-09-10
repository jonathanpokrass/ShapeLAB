function [hkm] = calcHKM(evecs1, evals1, ts, idx)
 hkm = evecs1 * diag(exp(-(ts(idx)) * abs(evals1))) * evecs1';
end