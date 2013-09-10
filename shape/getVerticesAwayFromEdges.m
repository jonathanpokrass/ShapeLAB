
function [idxs,d] = getVerticesAwayFromEdges(shape,dist,w)

if nargin < 3
	w = [];
end
N = length(shape.X(:));

if dist == 0 && nargout == 1
    idxs = 1:N;
    return
end

idxs = getEdges(shape);

%filter out idxs by w
if(length(w) > 0)
	idxs = idxs(find(w(idxs) == 0)); 
end

  f = fastmarchmex('init', int32(shape.TRIV-1), double(shape.X(:)), double(shape.Y(:)), double(shape.Z(:)));
  d = zeros(N,size(idxs,1));
  for k = 1:size(idxs,1)
    source = repmat(Inf, [N 1]);
    source(idxs(k)) = 0;
    d(:,k) = fastmarchmex('march', f, double(source));
    %d(:,k) = sqrt(sum([(shape.X - shape.X(idxs(k))).^2, (shape.Y - shape.Y(idxs(k))).^2, (shape.Z - shape.Z(idxs(k))).^2],2));
  end
   d(d>=9999999) = Inf; %!!

   fastmarchmex('deinit', f);

d = min(d');

if isempty(d)
  idxs = 1:length(shape.X(:));
else
  [tmp idxs] = find(d > dist);
end
    
end
