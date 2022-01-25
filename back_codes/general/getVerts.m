function vidx = getVerts(g, c)
% returns all the vertices of cell 'c' in a clockwise direction.

bs = g.cells{c+1}; % bonds of cell 'c'
if ~g.is_multiVM
    vidx = g.bonds(bs,1);
else
    vidx = [];
    for i=1:length(bs)
        b = bs(i); % current bond
        vidx = [vidx; g.bonds(b,1); g.bverts{b}']; 
    end
end

end