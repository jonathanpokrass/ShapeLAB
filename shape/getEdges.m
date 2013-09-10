function [idxs] = getEdges(shape)

    %find edges
    %edge is edge when it belongs only to one triangle
    %first remove duplicate vertices from triangles
    %% remove duplicated vertices
%    verts = [shape.X(:) shape.Y(:) shape.Z(:)];
%    [uniq_verts,tmp,uniq_lookup] = unique(verts,'rows','first');
%    shape.X = uniq_verts(:,1);
%    shape.Y = uniq_verts(:,2);
%    shape.Z = uniq_verts(:,3);
%%    shape.TRIV = zeros(Ntri*3,1);
    Ntri = size(shape.TRIV, 1);
%    shape.TRIV(:) = uniq_lookup(shape.TRIV(:));
%    shape.TRIV = reshape(shape.TRIV,Ntri,3);
    %%
    
    E = [shape.TRIV(:,1) shape.TRIV(:,2); 
     shape.TRIV(:,1) shape.TRIV(:,3);
     shape.TRIV(:,2) shape.TRIV(:,3)];
    E = sort(E,2);

    [E m n] = unique(E,'rows');
    c = accumarray(n,1);

    idxs = E( c == 1,:);

    idxs = unique(idxs); %% this is the vertex idxs of the edges
end
