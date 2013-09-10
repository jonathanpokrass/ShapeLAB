function [] = plotMatchingPoints(shape1, shape2, points1, points2, varargin)

  Nv1 = length(shape1.X);
  Nv2 = length(shape2.X);

  defaultOpt.ptColor1 = 1:Nv1;
  defaultOpt.ptColor2 = 1:Nv2;
  defaultOpt.t2 = [0];
  defaultOpt.drawSphere = 1;
	opt = parseOpt(defaultOpt, varargin{:});

  bbox1 = [];
  bbox1.tl = [min(shape1.X) min(shape1.Y) min(shape1.Z)];
  bbox1.br = [max(shape1.X) max(shape1.Y) max(shape1.Z)];

  shape1width = bbox1.br(1) - bbox1.tl(1);
  shape1centerX = 0.5 * (bbox1.br(1) + bbox1.tl(1));

  bbox2 = [];
  bbox2.tl = [min(shape2.X) min(shape2.Y) min(shape2.Z)];
  bbox2.br = [max(shape2.X) max(shape2.Y) max(shape2.Z)];

  shape2width = bbox2.br(1) - bbox2.tl(1);
  shape2centerX = 0.5 * (bbox2.br(1) + bbox2.tl(1));


  shape2offX = 0.5 * ( shape1width  + 1.5 * shape2width) + shape1centerX - shape2centerX;

  shape2.X = shape2.X + shape2offX + opt.t2;
%  points2.X = points2.X + shape2offX;
  cmap = colormap;
  idxToColMap = size(cmap,1) / (max(opt.ptColor1) - min(opt.ptColor1));

  drawShape(shape1, 'ptColor', opt.ptColor1);
  lighting none;
  hold on;
  drawShape(shape2, 'ptColor', opt.ptColor2);
  lighting none;

  if opt.drawSphere
    [sphX sphY sphZ] = sphere;
  end
  for i = 1:length(points1)
    if ~(points1(i)==0 || points2(i) == 0)
      
    color = cmap(min([size(cmap,1), ceil(idxToColMap * opt.ptColor1(abs(points1(i))))]),:);
    
    plot3([shape1.X(abs(points1(i))), shape2.X(abs(points2(i)))], [shape1.Y(abs(points1(i))) shape2.Y(abs(points2(i)))], [shape1.Z(abs(points1(i))) shape2.Z(abs(points2(i)))], ...
       '-o', ...
       'Color', color, ...
       'LineWidth', 2,...
       'MarkerEdgeColor','k',...
       'MarkerFaceColor',color,...
       'MarkerSize',12);
     end
    
    if opt.drawSphere
      radius = 0.03;
      cax = caxis;
      
      if points2(i) <= 0 || points1(i) <= 0
        cax = cax(2)*0.9;
      else
        cax = cax(1)*1.1;
      end
        
      surf(sphX*radius + shape1.X(abs(points1(i))), sphY*radius + shape1.Y(abs(points1(i))), sphZ*radius + shape1.Z(abs(points1(i))), repmat(cax,length(sphZ), length(sphZ)));
      shading flat
      if points2(i) ~= 0
        surf(sphX*radius + shape2.X(abs(points2(i))), sphY*radius + shape2.Y(abs(points2(i))), sphZ*radius + shape2.Z(abs(points2(i))), repmat(cax,length(sphZ), length(sphZ)));
        shading flat
      end
    end
  end

  lighting none;
  hold off;
end


