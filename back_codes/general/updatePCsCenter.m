function g = updatePCsCenter(g)
% updates the centroid position of the HCs
if ~isfield(g,'top_boundary_cells') , return; end

PCs = g.top_boundary_cells;
if size(PCs, 1) == 1, PCs = PCs'; end
for c=PCs
    if ~g.dead(c)
        g.centroid(c,:) = cellCOM(g,c);
    end
end

end