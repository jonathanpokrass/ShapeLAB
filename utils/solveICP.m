% X = argmin || A * X - B || s.t X' * X = I
function [X, U, S, V] = solveICP(A, B)
  [U,S,V] = svd(A'*B);
  U(:,end)=U(:,end)*det(U*V');
  X=U*V';
end