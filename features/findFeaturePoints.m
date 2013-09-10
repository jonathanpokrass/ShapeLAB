function [pts] = findFeaturePoints(s1, varargin)
defaultOpt.evecs = [];
defaultOpt.evals = [];

opt = parseOpt(defaultOpt, varargin{:});
if isempty(opt.evecs) || isempty(opt.evals)
  [opt.evecs opt.evals] = calcLaplacianBasis(s1, 100);
end

wks1 = calcWKS(opt.evecs, opt.evals, 100);

nv = length(s1.X(:));

% 1-ring neighbours
r1neigh = tri2adj(s1.TRIV);
%[i1, j1, k] = find(r1neigh);
% 2-ring neighbours
r2neigh = r1neigh;
parfor v = 1:nv
  r2neigh(v,:) = r2neigh(v,:) + sum(r1neigh(r1neigh(v,:)>0,:));
end

% find local maximas
t = 1;
neigh = r2neigh;
nextCandi = 1:nv;
for t = 80:5:80
  candi = nextCandi;
  nextCandi = [];
  for i = 1:length(candi) 
    x = candi(i);
    if all(wks1((neigh(x,:) > 0), t) <= wks1(x,t))
  %    disp('found')
      nextCandi = [nextCandi, x];
    end
  end
end
pts = nextCandi;
desc = wks1(pts,:);
end
