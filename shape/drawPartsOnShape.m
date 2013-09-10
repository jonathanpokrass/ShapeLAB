% f is the membership function: each column says whether the vertex belongs or not
function [] = drawPartsOnShape(shape, f, varargin)
  defaultOpts = [];
  defaultOpts.strip = 1;
  defaultOpts.distFromEdge = 0;
  defaultOpts.sortBySize = 0;
  defaultOpts.xoff = 0;
  defaultOpts.drawEdge = 0;
  defaultOpts.drawRest = 1;
  
  opt = parseOpt(defaultOpts, varargin{:});
  
  numParts = length(f(1,:));
  nvertices = length(f(:,1));
  
  if opt.sortBySize || ~opt.strip
    % sort draw order by number of vertices
    areas = calcAreas(shape);
    A = repmat(areas', size(f, 1), 1);
    
    [tmp, drawOrder] = sort(-sum(A * f));
  else
    drawOrder = 1:numParts;
  end
  

  xoff = opt.xoff;
  if opt.strip
  
    if opt.distFromEdge > 0
      %find existing edges
      existingEdges = getEdges(shape);
    end
    
    hold on;
    bbox = [];
    bbox.minx = Inf; bbox.maxx = 0;
    bbox.miny = Inf; bbox.maxy = -Inf;
    bbox.minz = Inf; bbox.maxz = -Inf;
    bboxList(1) = bbox;
    prevXoff = xoff;
  for i = 1:numParts
    j = drawOrder(i);

    
    [part, shapeToPart] = cutPart(shape, find(f(:,j)));
    
    if opt.distFromEdge > 0
      keepedEdgesOnPart = shapeToPart(existingEdges);
      keepedEdgesOnPart = keepedEdgesOnPart(keepedEdgesOnPart > 0);
    else
      keepedEdgesOnPart = [];
    end
    
    if opt.strip
      tmpPart = part;
      
      tmpPart.X = tmpPart.X(:) + prevXoff;
      partBbox = getBbox(tmpPart);
      
      overlapping = 0;
      for bboxnum = 1:length(bboxList)
        if isPointInsideBox(tmpPart, bboxList(bboxnum), 'minPercent', 0.3)
          overlapping = 1;
          break;
        end
      end
      % Move only if overlapping 
      if (overlapping || opt.drawRest) || i == 1
        %if i > 1
          xoff = (bbox.maxx - min(part.X(:))); %5%       
        %end
          part.X(:) = part.X(:) + xoff;
          prevXoff = xoff;
          bbox = getBbox(part);
          bboxList = bbox;
      else
        bboxList(end + 1) = partBbox;
        
        part = tmpPart;
        bbox = uniteBbox(partBbox, bbox);
        
      end
    else
      part.X(:) = part.X(:) + opt.xoff;
    end

    drawPart(part, [], 'color', j, 'distFromEdge', opt.distFromEdge, 'drawEdge', opt.drawEdge, 'ignoredEdges', keepedEdgesOnPart);
  end
  else
    % generate new shape that is split nicely between parts
    partsColoring = repmat(-1, length(shape.X), 1);
    new_shape = shape;
    for i = 1:numParts
      j = drawOrder(i);
     
      %find existing edges
      existingEdges = getEdges(new_shape);

      partIdxs = find(f(:,j));
      part = remove_unused_triangles(new_shape, partIdxs);
      
      %figure;
      %drawShape(part);
      % need the following because if we have unused vertices they will 
      % produce edges!
      %[part2 newShapeToPart partToNewShape] = remove_unused_vertices(part, 'keepDublicates', 1);
      
      ei = zeros(length(part.X(:)), 1);
      if ~isempty(existingEdges)
        ei(existingEdges) = 1;
      end
      partDistEdge = [];
      if opt.distFromEdge > 0
        [tmp, partDistEdge] = getVerticesAwayFromEdges(part, opt.distFromEdge, ei);
      end
      if ~isempty(partDistEdge)
        distEdge = repmat(-1, length(new_shape.X(:)),1);
        %distEdge(partToNewShape) = partDistEdge;
        distEdge = partDistEdge;
        distEdge(distEdge == Inf) = -1;
      else
        distEdge = [];
      end
      
      if length(distEdge) > 0 && opt.distFromEdge > 0
        [new_shape, edges] = geodesic_contour(new_shape, distEdge, opt.distFromEdge);
      else
        edges.X = [];
        %vidxs = sort(unique(part.TRIV(:)));
        %new_shape = part;
        new_shape.tri_labels = ones(size(new_shape.TRIV(:,1)));
        new_shape.tri_labels(sum(ismember(shape.TRIV', partIdxs)) == 3) = 0;
      end

      select = new_shape;
    	rest = select;

      select.TRIV =new_shape.TRIV(find(new_shape.tri_labels==0),:);
      rest.TRIV = new_shape.TRIV(find(new_shape.tri_labels==1),:);
      
      if ~isempty(partsColoring)
        partsColoring(end : length(new_shape.X(:))) = -1;
      
        pc1 = partsColoring(rest.TRIV(:, 1));
        pc2 = partsColoring(rest.TRIV(:, 2));
        pc3 = partsColoring(rest.TRIV(:, 3));
        pcRest = max(max(pc1, pc2), pc3);

        partsColoring(rest.TRIV(:, 1)) = pcRest;
        partsColoring(rest.TRIV(:, 2)) = pcRest;
        partsColoring(rest.TRIV(:, 3)) = pcRest;
      else
       % partsColoring = repmat(-1, size(rest.X(:),1));
      end
      
      %new_shape = select;
      vidxs = sort(unique(select.TRIV(:)));
      %vidxs2 = sort(unique(rest.TRIV(:)));
      partsColoring(vidxs) = j;
      %partsColoring(vidxs2) = 200;
    end
    %partsColoring(end:length(new_shape.X(:))) = 0;

    new_shape.X(:) = new_shape.X(:) + xoff;
    drawShape(new_shape, 'ptColor', partsColoring);
    
    %if (length(edges.X(:)) > 0) && opt.drawEdge
    %  plot3(edges.X,edges.Y,edges.Z,'k');
    %
    % end
  end
  
  %hold off;
%  h = trisurf(shape.TRIV, shape.X, shape.Y, shape.Z);
%  set(h, 'FaceVertexCData',  partsColoring)%, 'EdgeColor', 'none');
%
%  hold on;
%   for i = 1:numParts
%     %  drawPart(shape1, 0, f1(:, i), 'none', 0.99, 0.7, struct('drawRest', 0));
%     hold on;
%   end
%  axis image; axis off;
%  shading flat;
%  lighting phong;
	axis image;axis off;
	view([-20 30]);	
  shading flat;
end

function [bbox] = getBbox(shape)
 bbox.minx = min(shape.X(:));
 bbox.maxx = max(shape.X(:));
 bbox.miny = min(shape.Y(:));
 bbox.maxy = max(shape.Y(:));
 bbox.minz = min(shape.Z(:));
 bbox.maxz = max(shape.Z(:));
end

function [isinside] = isPointInsideBox(point, bbox, varargin)
  defaultOpt.minPercent = 0.0;
  opt = parseOpt(defaultOpt, varargin{:});
  
  if ~isstruct(point)
    x = point(:,1); y = point(:,2); z = point(:, 3);
  else
    x = point.X(:); y = point.Y(:); z = point.Z(:);
  end
  ox = and( x  > bbox.minx, x < bbox.maxx);
  oy = and( y  > bbox.miny, y < bbox.maxy);
  oz = and( z  > bbox.minz, z < bbox.maxz);
  numPointsInside = sum(and(and(ox, oy), oz));
  isinside = (numPointsInside / length(x))  > opt.minPercent;
     
end

function [newBbox] = uniteBbox(bbox1, bbox2)
  newBbox.minx = min(bbox1.minx, bbox2.minx);
  newBbox.maxx = max(bbox1.maxx, bbox2.maxx);
  newBbox.miny = min(bbox1.miny, bbox2.miny);
  newBbox.maxy = max(bbox1.maxy, bbox2.maxy);
  newBbox.minz = min(bbox1.minz, bbox2.minz);
  newBbox.maxz = max(bbox1.maxz, bbox2.maxz);
end

function [intersects] = isIntersecting(bbox1, bbox2)
if isPointInsideBox([bbox1.minx bbox1.miny bbox1.minz], bbox2) || ...
   isPointInsideBox([bbox1.minx bbox1.miny bbox1.maxz], bbox2) || ...
   isPointInsideBox([bbox1.minx bbox1.maxy bbox1.minz], bbox2) || ...        
   isPointInsideBox([bbox1.minx bbox1.maxy bbox1.maxz], bbox2) || ...
   isPointInsideBox([bbox1.maxx bbox1.miny bbox1.minz], bbox2) || ...
   isPointInsideBox([bbox1.maxx bbox1.miny bbox1.maxz], bbox2) || ...
   isPointInsideBox([bbox1.maxx bbox1.maxy bbox1.minz], bbox2) || ...        
   isPointInsideBox([bbox1.maxx bbox1.maxy bbox1.maxz], bbox2)

   intersects = 1;
 else
   intersects = 0;
 end 
end
