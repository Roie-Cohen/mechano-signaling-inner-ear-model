function g = updateCellsCenter(g)
% updates the centroid position of all cells

for c=1:length(g.cells)-1
    if ~g.dead(c)
        g.centroid(c,:) = cellCOM(g,c);
    end
end

end