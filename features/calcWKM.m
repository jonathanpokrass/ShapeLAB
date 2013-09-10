function [wkm] = calcWKM(evecs1, evals1, numSamples, idx, segs)
if nargin < 4
    doSegs = 0;
else
    doSegs = 1;
end
if doSegs, wkm = cell(numSamples,1); end;
  
  evals1 = abs(evals1);
  emin = log(abs(evals1(2)));
  emax = log(abs(evals1(end)));
  s = 7*(emax - emin)/numSamples;
  emin = emin + 2*s;
  emax = emax - 2*s;
  tn = linspace(emin, emax, numSamples);
  ce = sum( exp( -((bsxfun(@minus, log(evals1(:)) , tn(:)')).^2)./(2*s*s)));
  ce(ce == 0) = 1;
  
  tmp = bsxfun(@minus, log(evals1(:)), tn(:)');
  
  calcWkm = @(k) evecs1*diag(exp( -(tmp(:,k).^2)./(2*s*s)) ./ ce(k)) * evecs1';
  if doSegs
      parfor i = 1:numSamples
        wkm{i} = feval(calcWkm, i) * segs;
        disp(i)
      end
  else
     wkm = calcWkm(idx);
  end
end