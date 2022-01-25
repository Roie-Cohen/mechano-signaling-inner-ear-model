function v=extractverts(g)
l = length(g.verts);
v=zeros(2*l,1);
v(1:2:2*l-1)=g.verts(:,1);
v(2:2:2*l)=g.verts(:,2);
end