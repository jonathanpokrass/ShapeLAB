function R = rotation_matrix(tx,ty,tz)

R = [ 1 0 0; 0 cos(tx) sin(tx); 0 -sin(tx) cos(tx)] * [ cos(ty) 0 sin(ty); 0 1 0; -sin(ty) 0 cos(ty)] * [ cos(tz) sin(tz) 0; -sin(tz) cos(tz) 0; 0 0 1];

