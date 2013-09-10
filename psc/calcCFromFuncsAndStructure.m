function [C, O, basis1, basis2] = calcCFromFuncsAndStructure(shape1, shape2, f1, f2, varargin)
% Finds C between the functions f1 on shape1 and f2 on shape2
% solves the the problem
%  argmin_{C} 
%      |basis2 * f2 - C * basis1 * f1 + O|^2_2 
%          + lam1 * |W.*C|_1 + lam2 * |O|_{21} 
% W is chosen so that C will be diagonal
% f1(:,1) is the first function on shape1.
% f2(:,2) is the first function on shape2.
%
% C is the resulting ... 
% Author Jonathan Pokrass
defaultOpts.lam2 = 2; %0.5 ; 2
defaultOpts.lam1 = 0.1;%0.05;% 0.05; %0.05; 1
defaultOpts.numVecsInBasis = 100;
defaultOpts.basis1 = [];
defaultOpts.basis2 = [];
defaultOpts.areas1 = [];
defaultOpts.areas2 = [];
defaultOpts.W1 = [];
defaultOpts.W2 = [];
defaultOpts.C0 = [];
defaultOpts.debug = 0;
defaultOpts.useCommutativity = 0;
opt = parseOpt(defaultOpts, varargin{:});

numFuncGroups1 = 1; 
if iscell(f1)
  numFuncGroups1 = numel(f1);
else
  tmp = f1; f1 = {}; f1{1} = tmp;
end

numFuncGroups2 = 1; 
if iscell(f2)
  numFuncGroups2 = numel(f2);
else
  tmp = f2; f2 = {}; f2{1} = tmp;
end


timeit = [];
if isempty(opt.basis1) || isempty(opt.basis2)
  tic;
  [basis1 evals1 areas1 W1] = calcLaplacianBasis(shape1, opt.numVecsInBasis);
  basis1 = sqrt(areas1) * basis1;
  
  [basis2 evals2 areas2 W2] = calcLaplacianBasis(shape2, opt.numVecsInBasis);
  basis2 = sqrt(areas2) * basis2;
  timeit.basis = toc;
else
  basis1 = opt.basis1; basis2 = opt.basis2;
  areas1 = opt.areas1; areas2 = opt.areas2;
  if opt.useCommutativity
    W1 = opt.W1;W2 = opt.W2;
  end
  opt.numVecsInBasis = size(opt.basis1,2);
end

if opt.useCommutativity
    D1 =(sqrt(areas1) * basis1)' * W1 * (sqrt(areas1) * basis1) ;%eye(opt.numVecsInBasis);% ; %diag(evals1)
    D2 =(sqrt(areas2) * basis2)' * W2 * (sqrt(areas2) * basis2);%eye(opt.numVecsInBasis);%basis2' * W2 * basis2; %diag(evals2)
end
numFuncs1 = size(f1{1}, 2);
numFuncs2 = size(f2{1}, 2);
switchedShapes = 0;
if (numFuncs1 ~= numFuncs2)
  error('Number of functions must be the same, use findPermutation)')
end
a1 = [];
b1 = [];
Y = [];
for grp=1:numFuncGroups1  
  a = f1{grp}' * basis1;
  a1 = [a1;a];
  if opt.useCommutativity
    a1 = [a1; a * D1];
  end
end
for grp=1:numFuncGroups2
  b = f2{grp}' * basis2;
  b1 = [b1; b];
  if opt.useCommutativity
    b1 = [b1; b * D2];
  end
end
if opt.useCommutativity
  numFuncGroups2 = numFuncGroups2 * 2;
  numFuncGroups1 = numFuncGroups1 * 2;
end

if isempty(opt.C0)
    C = zeros(opt.numVecsInBasis, opt.numVecsInBasis);%repmat(1/opt.numVecsInBasis, opt.numVecsInBasis, opt.numVecsInBasis);
else
    C = opt.C0;
end


printlog('Strting minimization');
tic;

    %% Find best C for best permutation (eye in this case)
   Y = b1;
   lam = calcOptimalLam() * opt.lam1;
   
   [C, funVals, tmp2, O] = solveLSRobustL1(a1, Y, lam, opt.lam2, 'x0', C, 'rsL2', 0);%, 'maxIter', 1e5);
  
timeit.minimization = toc;
printlog('Finished minimization');

%%  
if opt.debug
  figure;
  subplot(1, 2, 1);
  imagesc(abs(C));
  title('C');
  subplot(1,2,2);
  imagesc(abs(O));
  title('O');
end


function [optimalLam] = calcOptimalLam()
    optimalLambda = 0;
    for iii=1:1%size(Y,2)
      optimalLambda = max(max(abs(a1' * Y(:, iii))), optimalLambda);
    end
    %mask = (1 - Cmask(1:opt.numVecsInBasis,1:opt.numVecsInBasis));
    %mask(mask < 0.5) = 0;
    nvib = opt.numVecsInBasis;
    optimalLam = (ones(nvib) - eye(nvib));%optimalLambda .* mask;%(ones(size(C)) - eye(size(C)));
     
     for k = [1:nvib-1]
       if k < nvib/2
         coeff = 1/k;
        else
         coeff = sqrt(1/(k));
        end
        optimalLam([k+1:nvib], k) = optimalLam([k+1:nvib], k) .* optimalLambda * coeff;
        optimalLam(k,[k+1:nvib]) = optimalLam(k,[k+1:nvib]) .* optimalLambda * coeff;
     end
  end
end


