function dE=dHlength(g,c)
% calculates the energy deferential derived from the bonds length and
% perimeter elements in the energy function.
l_fac = g.paras(2)/sqrt(g.A0);
L_fac = g.paras(3)*g.ctrc(c)/(g.A0);
dE = zeros(2*length(g.verts),1);
L = cellPerimeter(g,c);
vidx=getVerts(g, c); % an array of the vertices indices of the cell
vert = getRelativePosition(g,vidx); % the position of the vertices
nv=length(vidx);
bs = g.cells{c+1}; % bonds of cell 'c' (real bonds)
nb = length(bs);
real_vidx = g.bonds(bs, 1);
tensions = g.bonds(bs,5);
b_prev = nb;
b = 0;
for j=1:nv
    prev = mod(j-2,nv)+1; % the previous vertex index (in vert)
    next = mod(j,nv)+1;   % the next vertex
    n1 = norm(vert(j,:)-vert(prev,:));
    n2 = norm(vert(j,:)-vert(next,:));
    
    if ismember(vidx(j), real_vidx), b = mod(b,nb)+1; end
    if ismember(vidx(prev), real_vidx), b_prev = mod(b_prev,nb)+1; end
    
    % changed caused due to the bond connected to the previous vertex
    if(n1>0.0001)
        gamma = tensions(b_prev);
        dE(2*vidx(j)-1) = (L_fac*L + l_fac*gamma)*(vert(j,1)-vert(prev,1))/n1;
        dE(2*vidx(j)) =   (L_fac*L + l_fac*gamma)*(vert(j,2)-vert(prev,2))/n1;
    end
    
    % changed caused due to the bond connected to the next vertex
    if(n2>0.0001)
        gamma = tensions(b);
        dE(2*vidx(j)-1) =  dE(2*vidx(j)-1) + (L_fac*L + l_fac*gamma)*(vert(j,1)-vert(next,1))/n2;
        dE(2*vidx(j)) =   dE(2*vidx(j)) + (L_fac*L + l_fac*gamma)*(vert(j,2)-vert(next,2))/n2;
    end
    
end

end
