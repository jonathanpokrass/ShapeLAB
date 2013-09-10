function [] = saveoff(shape, filename)

f = fopen(filename, 'wt');

fprintf(f, 'OFF\n');
nv = length(shape.X);
nt = length(shape.TRIV(:,1));
fprintf(f, '%d %d %d\n', nv, nt, nv + nt - 1);

%write points
fprintf(f, '%.6f %.6f %.6f\n', [shape.X(:) shape.Y(:) shape.Z(:)]');

%write triangles
fprintf(f, '3 %d %d %d\n', shape.TRIV' - 1);

fclose(f);

