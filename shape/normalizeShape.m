function [shapeOut] = normalizeShape(shape, varargin)
  defaultOpt.areas = [];
  opt = parseOpt(defaultOpt, varargin{:});
  if isempty(opt.areas)
    opt.areas = calcAreas(shape);
  end
  sa = sqrt(sum(opt.areas(:)));
  Xc = mean(shape.X);
  Yc = mean(shape.Y);
  Zc = mean(shape.Z);
  shapeOut = shape;
  shapeOut.X = (shape.X - Xc) / sa;% + Xc;
  shapeOut.Y = (shape.Y - Yc) / sa;% + Yc;
  shapeOut.Z = (shape.Z - Zc) / sa;% + Zc;
  shapeOut.Z = shapeOut.Z - min(shapeOut.Z);
end