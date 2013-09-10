function [F] = loadSegments(filename)
f = fopen(filename);
if f < 0
  F = [];
  return
end
labels = fscanf(f, '%d\n');
nsegs = max(labels);
n = length(labels);
segs = find(labels > 0);
F = sparse(labels(segs), segs, ones(length(segs),1), nsegs, n);
F = F';
fclose(f);
end