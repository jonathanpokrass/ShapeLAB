function [indF, wksF, wkmF]  = generateFuncGrpsFromSeg(evecs1, evals1, areas1, seg)
 %1 - indicators
 indF = seg;
 if nargout > 1
   numSamples = 100;
  %2 - WKS Field
  wksF = cell(numSamples,1);
  wks = calcWKS(evecs1, evals1, numSamples);
  for i=1:size(wks,2)
    wksF{i} = repmat(wks(:,i),1,size(seg,2)) .* seg; 
  end
  
  if nargout > 2
    %3 - WKM
    wkmF = calcWKM(evecs1, evals1, numSamples, i, seg);
  end
 end
end
