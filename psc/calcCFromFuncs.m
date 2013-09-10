function [C] = calcCFromFuncs(fg1, fg2, basis1, basis2, W1, W2)  

numFuncGroups1 = 1; 
if iscell(fg1)
  numFuncGroups1 = numel(fg1);
else
  tmp = fg1; fg1 = {}; fg1{1} = tmp;
end

numFuncGroups2 = 1; 
if iscell(fg2)
  numFuncGroups2 = numel(fg2);
else
  tmp = fg2; fg2 = {}; fg2{1} = tmp;
end


if nargin > 4
    %For commutativity
    D1 = basis1'*W1*basis1;
    D2 = basis2'*W2*basis2;
end
a1 = [];
b1 = [];
%b*D2 = a*D1*C 
for i=1:numel(fg1)
 a = fg1{i}' * basis1;
 b = fg2{i}' * basis2;
  
 if nargin > 4
    a1 = [a1; a; a*D1];
    b1 =  [b1; b; b*D2];
 else
   a1 = [a1; a];
   b1 = [b1; b];
 end
end
%Y = kron(eye(size(a1,1) / size(fg1{1},1)),P) * b1;
Y = b1;
C = solveICP(a1, Y);

end
