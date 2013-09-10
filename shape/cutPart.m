function [newshape, lookuptbl] = cutPart(shape, partIdxs)
  newshape = remove_unused_triangles(shape, partIdxs);
  [newshape lookuptbl] = remove_unused_vertices(newshape);
end

  