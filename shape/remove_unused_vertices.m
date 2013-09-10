% @author Jonathan Pokrass
function [newsurf, shapeToPartTbl partToShapeTbl] = remove_unused_vertices(surface, varargin)

defaultOpt.keepDublicates = 0;
defaultOpt.use_ann = 0;
opt = parseOpt(defaultOpt, varargin{:});
newsurf = surface;
Ntri = size(surface.TRIV,1);

vidxs = sort(unique(surface.TRIV(:)));

Nv = size(vidxs,1);
newsurf.X = surface.X(vidxs);
newsurf.Y = surface.Y(vidxs);
newsurf.Z = surface.Z(vidxs);

lookuptbl = sparse(size(surface.X,1),1);
lookuptbl(vidxs) = 1:Nv;

%remove duplicated vertices
verts = [newsurf.X(:) newsurf.Y(:) newsurf.Z(:)];
if ~opt.keepDublicates
if opt.use_ann
	tree = ann('init', [newsurf.X(:) newsurf.Y(:) newsurf.Z(:)]');
	[treeIdx] = ann('search', tree, verts',10,'radius',0.001, 'search_sch','fr' );
	ann('deinit', tree);
	ann('close');
	treeIdx = sort(treeIdx); 

	%uniq_lookup = treeIdx(10,:)
	[uniq_verts,tmp, uniq_lookup] = unique(treeIdx(10,:));
	newsurf.X = newsurf.X(uniq_verts);
	newsurf.Y = newsurf.Y(uniq_verts);
	newsurf.Z = newsurf.Z(uniq_verts);
  newsurf.TRIV = zeros(Ntri*3,1);

  newsurf.TRIV(:) = uniq_lookup(lookuptbl(surface.TRIV(:)));
else
	[uniq_verts,tmp,uniq_lookup] = unique(verts,'rows','first');
	newsurf.X = uniq_verts(:,1);
	newsurf.Y = uniq_verts(:,2);
	newsurf.Z = uniq_verts(:,3);
  newsurf.TRIV = zeros(Ntri*3,1);

  newsurf.TRIV(:) = uniq_lookup(lookuptbl(surface.TRIV(:)));
end
  uniq_lookup(end + 1) = 0;
  dummy = length(uniq_lookup);
  lookuptbl(lookuptbl == 0) = dummy;
  shapeToPartTbl = uniq_lookup(lookuptbl);
else
  newsurf.TRIV(:) = lookuptbl(surface.TRIV(:));
  shapeToPartTbl = lookuptbl;
end

[tmp partToShapeTbl] = sort(unique(shapeToPartTbl));
partToShapeTbl = partToShapeTbl(find(tmp > 0):end);
newsurf.TRIV = reshape(newsurf.TRIV,Ntri,3);
end
