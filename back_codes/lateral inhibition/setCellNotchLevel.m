function g = setCellNotchLevel(g, cs, notch_level)
% sets notch level for cells
% 'cs': cells indices
% 'notch_level': the level of notch to set for each cell

for i=1:length(cs)
    c = cs(i);
    g.LImodel.cell_notch(c) = notch_level;
    % distribute over the cell bonds
    bs = g.cells{c+1}; % the bonds of the cell
    L = cellPerimeter(g, c);
    g.LImodel.bond_notch(bs) = g.LImodel.cell_notch(c)/L;
end

end