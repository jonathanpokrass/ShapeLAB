function [shape] = loadShape(filename, varargin)

defaultOpt.scale = 1;
defaultOpt.numVertices = -1;
defaultOpt.switchXZ = 0;
options = parseOpt(defaultOpt, varargin{:});

[path, name, ext] = fileparts(filename);

switch ext
case '.mat'
	shapest = load(filename);
	if(isfield(shapest, 'surface'))
		shape = shapest.surface;
	else
		if(isfield(shapest, 'shape'))
			shape = shapest.shape;
		end
    end
    if(isfield(shapest, 't'))
        shape.t = shapest.t;
    end
    if(isfield(shapest, 'u'))
        shape.u = shapest.u;
    end
case '.off'
	shape = loadoff(filename);
	shape.X = shape.X;% * 500;
	shape.Y = shape.Y;% * 500;
	shape.Z = shape.Z;% * 500;
  if options.switchXZ
    tmpshape = shape;
    shape.X = tmpshape.Z;
    shape.Z = tmpshape.X;
  end
case '.ply'
	[TRIV, V] = ply_read(filename,'tri');
	shape.TRIV = TRIV';
	shape.X = permute(V(3,:),[2 1])*100;
	shape.Y = permute(V(2,:),[2 1])*100;
	shape.Z = permute(V(1,:),[2 1])*100;
end

if options.numVertices > 0
  shape = remesh(shape ,struct('vertices', options.numVertices, 'placement', 0));
	shape = remove_unused_vertices(shape);
end
shape.X = shape.X * options.scale;
shape.Y = shape.Y * options.scale;
shape.Z = shape.Z * options.scale;
end
