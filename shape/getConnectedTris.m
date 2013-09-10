%% 
% v - list of triangles that share vertex n
%     around n
%%
function [tris] = getConnectedTris(tris,n)

tris = tris';
tris = tris(:);
tris = floor((find(tris == n) - 1)/3  + 1);

end

