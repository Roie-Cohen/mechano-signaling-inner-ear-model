function g = findNeighboringLinkedCells(g)
% find linked cells that contact each other
% Merges these cells into a single cell

linkedCells = find(g.linkedCells ~= 0);
while ~isempty(linkedCells)
    c1 = linkedCells(1);
    c2 = g.linkedCells(c1);
    c1_neighs = g.bonds(g.cells{c1+1},4);
    if ismember(c2, c1_neighs)
        g = mergeNeighboringCells(g, c1, c2);
        g.linkedCells(c1)=0;
        g.linkedCells(c2)=0;
    end
    linkedCells(linkedCells == c1 | linkedCells == c2) = [];
end

end