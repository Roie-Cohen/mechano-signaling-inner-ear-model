function area=cellarea(g,c)
% return the area of cell 'c'
vidx = getVerts(g, c);
vert = getRelativePosition(g,vidx);
vert = vert*g.scale; 
area = polyarea(vert(:,1),vert(:,2));
end