function L=cellPerimeter(g,c)

L=0;
vidx=getVerts(g, c);
vert = getRelativePosition(g,vidx);
nv=length(vidx);
for j=1:nv
 n = norm(vert(j,:)-vert(mod(j,nv)+1,:));
 L=L+n;
end
end