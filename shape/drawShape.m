function [] = drawShape(shape, varargin)
  
	Nv = length(shape.X);

  defaultOpt.useTriv = 1;
  if isempty(varargin) || isstruct(varargin{1}) || ischar(varargin{1})
    defaultOpt.ptColor = shape.X(:);%1:Nv;
  else
    defaultOpt.ptColor = varargin{1};
    varargin = {};
  end
	opt = parseOpt(defaultOpt, varargin{:});

  %opt.ptColor = (opt.ptColor - min(opt.ptColor)) ./ (max(opt.ptColor) - min(opt.ptColor));

	if opt.useTriv && isfield(shape, 'TRIV')
		% Visualize current sample
    if length(opt.ptColor) == 1
      opt.ptColor = repmat(opt.ptColor, Nv, 1);
    end
		h = trisurf(shape.TRIV, shape.X, shape.Y, shape.Z, opt.ptColor);
	else
		plot3k([shape.X(:) shape.Y(:) shape.Z(:)]);
	end
	axis equal;
	shading interp;
	lighting phong;
	%camlight head;
	axis off;
	%hold on;
end


