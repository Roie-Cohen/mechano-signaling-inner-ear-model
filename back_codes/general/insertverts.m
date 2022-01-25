function gout=insertverts(v,gin)
gout=gin;
gout.verts(:,1:2) = reshape(v,[2 length(gout.verts)])';
end