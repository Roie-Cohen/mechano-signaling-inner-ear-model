function g = updateHCsCenter(g)
% updates the centroid position of the HCs
if ~g.is_LImodel, return; end

HCs = g.LImodel.high_delta_cells;
if size(HCs, 1) == 1, HCs = HCs'; end
for c=HCs
    if ~g.dead(c)
        g.centroid(c,:) = cellCOM(g,c);
    end
end

end