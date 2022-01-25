function dE=dHarea(g,c)
% area term: 0.5*alpha*(A-A0)^2
% derivative with respect to 'x' coordinate:
% alpha*A0*(A/A0 - 1)*(dA/dx)
linked_cell = 0;
if isfield(g, "linkedCells")
    linked_cell = g.linkedCells(c);
end

A0 = g.areas(c);
dE = zeros(2*length(g.verts),1);
vidx = getVerts(g, c);
vlist = getRelativePosition(g,vidx);
l = length(vidx);
for i=1:l
    dE(2*vidx(i)-1) = 0.5*(vlist(mod(i-2,l)+1,2)-vlist(mod(i,l)+1,2)); % dA/dx
    dE(2*vidx(i)) = 0.5*(vlist(mod(i,l)+1,1)-vlist(mod(i-2,l)+1,1)); % dA/dy
end
A = cellarea(g,c);
if linked_cell==0
    dE = g.cmpr(c)*(A - A0)*dE;
else
    c1 = min(c, linked_cell); % the original cell
    c2 = max(c, linked_cell); % the new cell
    A1 = cellarea(g, c1);
    A2 = cellarea(g, c2);
    dE = g.cmpr(c)*(A1 + (1+3*(c==c2))*(A2 - A0))*dE;
end

dE = g.paras(1)*dE/(g.A0^2);

end