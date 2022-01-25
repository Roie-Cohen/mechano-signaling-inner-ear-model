function g = setCellDeltaLevel(g, cs, delta_level)
% sets delta level for cells
% 'cs': cells indices
% 'delta_level': the level of delta to set for each cell

for i=1:length(cs)
    c = cs(i);
    g.LImodel.cell_delta(c) = delta_level;
    % distribute over the cell bonds
    bs = g.cells{c+1}; % the bonds of the cell
    L = cellPerimeter(g, c);
    g.LImodel.bond_delta(bs) = g.LImodel.cell_delta(c)/L;
end

end