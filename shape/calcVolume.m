function [vol] = calcVolume(shape)
V = [shape.X(:) shape.Y(:) shape.Z(:)];
V = bsxfun(@minus, V, mean(V));

nt = size(shape.TRIV,1);
vol = 0;
for t = 1:nt
		x1 = shape.TRIV(t,1);
		x2 = shape.TRIV(t,2);
		x3 = shape.TRIV(t,3);
        
        vol = vol + abs(dot(V(x1,:), cross(V(x2,:), V(x3,:))));
end
vol = vol / 6;
end