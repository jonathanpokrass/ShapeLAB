function [] = drawPart(part, w, varargin) 
		defaultOpt = [];
    defaultOpt.drawRest = 0;
    defaultOpt.opacity = 0.7;
    defaultOpt.th = 0.99;
    defaultOpt.distFromEdge = 0;
	  defaultOpt.color = 1;
    defaultOpt.drawEdge = 0;
    defaultOpt.ignoredEdges = [];
    opt = parseOpt(defaultOpt, varargin{:});

	if(opt.distFromEdge ~= 0)
    ei = zeros(length(part.X(:)), 1);
    ei(opt.ignoredEdges) = 1;
	  [idxs,distEdge] = getVerticesAwayFromEdges(part, opt.distFromEdge, ei);
    
   if length(distEdge) ~= 0
    [new_part,edges] = geodesic_contour(part, distEdge, opt.distFromEdge);
   else
      edges.X = [];
			new_part = part;
			new_part.tri_labels = zeros(size(new_part.TRIV(:,1)));
   end
  else
		if(opt.th >= 0 && length(w) > 0)
	 		[new_part,edges] = geodesic_contour(part,w,opt.th);
 		else
      edges.X = [];
			new_part = part;
			new_part.tri_labels = zeros(size(new_part.TRIV(:,1)));
 		end
  end

	select = new_part;
	rest = select;

	select.TRIV =new_part.TRIV(find(new_part.tri_labels==0),:);
	
	rest.TRIV = new_part.TRIV(find(new_part.tri_labels==1),:);

	if(length(select.X) > 0)
     if length(opt.color) == 1
      opt.color = repmat(opt.color, length(select.X(:)), 1);
    end
		h1 = trisurf(select.TRIV, select.X, select.Y, select.Z, opt.color); 
		%shading flat; 
	end
	hold on
	if (length(edges.X(:)) > 0) && opt.drawEdge
		plot3(edges.X,edges.Y,edges.Z,'k');
	end
	if(opt.drawRest > 0)
    if(length(rest.X) > 0)
      h = trisurf(rest.TRIV, rest.X, rest.Y, rest.Z); 
		%shading flat;

      set(h,'FaceAlpha',opacity);
      set(h,'FaceColor', [1 1 0.95], 'EdgeColor', 'none');
    end
  end
	if(length(select.X) > 0) && 0
			if(opt.th >= 0)
        cmap = colormap;
			 %set(h1, 'FaceColor', cmap(opt.color,:), 'EdgeColor', 'none');
		  else
				map = real2rgb(w, thermal);
			 	map = reshape(map,size(map,1),3);
				set(h1,'FaceVertexCData', map);
				colormap(thermal);
				shading interp;
				colorbar;
			end
	end	

	lighting phong;
  shading flat;
	%camlight head;		
	axis image;axis off;
	view([-15 25]);	
end
	
