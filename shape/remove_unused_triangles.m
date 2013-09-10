function [newshape] = remove_unused_triangles(shape, partIdxs)
  % keep only triangles that have all the vertices in partIdxs
  loc = find(sum(ismember(shape.TRIV', partIdxs)) == 3);
%  loc = unique(floor((loc -1) ./ 3) + 1);
  shape.TRIV = shape.TRIV(loc,:);
  newshape = shape;
  %![newshape lookuptbl] = remove_unused_vertices(shape);
end

  