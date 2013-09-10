%
% Sparse Modeling of Intrinsic Correspondences. Computer Graphics Forum '13
% 
addpath('.\utils');
addpath('.\shape');
addpath('.\fmm');
addpath('.\features');
addpath('.\psc');

%% Load shapes and parts
s = {};
s{1} = loadShape('examples\gorilla1.off');
s{2} = loadShape('examples\gorilla14.off');

s{1}.parts = loadSegments('examples\gorilla1_seg.txt');
s{2}.parts = loadSegments('examples\gorilla14_seg.txt');
for k = 1:2
    figure(k);
    drawPartsOnShape(s{k}, s{k}.parts, 'strip', 0);light;lighting flat;    
end

%% Generate functions
useIndicators = 1;
useWKM = 0;
useWKS = 0;

nV = 20; %num eigenfunctions
numFeatureSamples = 100;

for k = 1:2
    s{k}.funcs = {};
    [s{k}.evecs, s{k}.evals, s{k}.areas, s{k}.W] = calcLaplacianBasis(s{k}, nV);
    s{k}.basis = sqrt(s{k}.areas) * s{k}.evecs;

    if useIndicators
        indF = s{k}.parts; 
        s{k}.funcs = {s{k}.funcs{:}, indF};
    end
    
    if useWKM
        wkmF = calcWKM(s{k}.evecs, s{k}.evals, numFeatureSamples, 1, s{k}.parts);
        s{k}.funcs = {s{k}.funcs{:}, wkmF{:}};
    end
    
    if useWKS
        wksF = cell(numFeatureSamples,1);
        wks = calcWKS(s{k}.evecs, s{k}.evals, numFeatureSamples);
        for i=1:size(wks,2)
            wksF{i} = repmat(wks(:,i),1,size(s{k}.parts,2)) .* s{k}.parts; 
        end
        s{k}.funcs = {s{k}.funcs{:}, wksF{:}};
    end
end  


%% Calculate correspondence

%not using C structure.
C = calcCFromFuncs(s{1}.funcs, s{2}.funcs, s{1}.basis, s{2}.basis); 
[shape1ToShape2 shape2ToShape1 C] = calcP2PFromC(s{1}, s{2}, C, s{1}.evecs, s{2}.evecs, 'debug', 1,'numRefinements', 30);
fprintf('Results without structre constraint\n Press any key to continue\n');
pause;
%using C structure.
[C2 O] = calcCFromFuncsAndStructure(s{1}, s{2}, s{1}.funcs, s{2}.funcs, 'basis1', s{1}.basis, 'basis2', s{2}.basis, 'debug',1);
[shape1ToShape2 shape2ToShape1 C2] = calcP2PFromC(s{1}, s{2}, C2, s{1}.evecs, s{2}.evecs, 'debug', 1,'numRefinements', 30);
fprintf('Results with structre constraint\n Press any key to continue\n');
