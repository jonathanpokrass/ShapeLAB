function [areas] =  calcAreas(shape) 
ntri = size(shape.TRIV,1);
nv = size(shape.X,1);
areas = zeros(nv,1);
V = [shape.X(:) shape.Y(:) shape.Z(:)]';
nt = size(shape.TRIV,1);
for t = 1:nt
		xi = shape.TRIV(t,1);
		xj = shape.TRIV(t,2);
		xk = shape.TRIV(t,3);
		
		Xj = [V(:,xj)-V(:,xi), V(:,xk)-V(:,xi)];
		areas([xi xj xk]) = areas([xi xj xk]) +  0.5*sqrt(det(Xj'*Xj))/3;
end
