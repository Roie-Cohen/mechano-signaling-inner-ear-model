function g = setTopPCs(g)
% sets any cell that borders the lateral inhibition region from the top
% (type=10 cells that border type=1 cells) as PCs (type = 2).

above_cells = find(g.type == 10);
for i=1:length(above_cells)
    c = above_cells(i);
    neighs = g.bonds(g.cells{c+1}, 4);
    neigh_type = g.type(neighs);
    if ismember(1, neigh_type)
        g.type(c) = 2;
    end
end


end