function [F] = calcWKS(evecs1, evals1, numSamples)
  evals1 = abs(evals1);
  emin = log(abs(evals1(2)));
  emax = log(abs(evals1(end)));
  s = 7*(emax - emin)/numSamples;
  emin = emin + 2*s;
  emax = emax - 2*s;
  tn = linspace(emin, emax, numSamples);
  ce = sum( exp( -((bsxfun(@minus, log(evals1(:)) , tn(:)')).^2)./(2*s*s)));
  ce(ce == 0) = 1;
  ce = 1./ce;
  tmp = bsxfun(@minus, log(evals1(:)), tn(:)');
  lognorm = (evecs1.^2)*exp( -(tmp.^2)./(2*s*s));
   F = lognorm * diag(ce);
  % F = F(:,1:numSamples);
  % F  = normalize(F, 'l1', 2);
end