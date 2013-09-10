function shape1 = rotate_shape(shape, R,t)

if nargin < 3, t = [0 0 0]'; end

shape1 = shape;

shape1.X = R(1,:)*[shape.X(:) shape.Y(:) shape.Z(:)]' + t(1);
shape1.Y = R(2,:)*[shape.X(:) shape.Y(:) shape.Z(:)]' + t(2);
shape1.Z = R(3,:)*[shape.X(:) shape.Y(:) shape.Z(:)]' + t(3);


shape1.X = shape1.X(:);
shape1.Y = shape1.Y(:);
shape1.Z = shape1.Z(:);


