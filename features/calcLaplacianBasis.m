function [evecs evals areas W1] = calcLaplacianBasis(shape, numVecsInBasis)
  %calculate Laplacian
  [W1, A1] = mshlp_matrix(shape,struct('dtype','cotangent'));
  %[W1, A1] = mshlp_matrix(shape,struct('dtype','euclidean', 'hs', 1.5, 'rho', 10, 'htype','psp'));
  W1 = (W1 + W1') * 0.5;
  %geodesic, %'euclidian'
  
  W1 = W1; %%%!!!
  %eigendecomposition
  Am1 = sparse([1:length(A1)], [1:length(A1)], A1);

  warning off;
  [evecs evals] = eigs(W1, Am1, numVecsInBasis, -1e-5, struct('disp',0));
  evals = diag(evals);
  warning on;
  areas = Am1;
end
