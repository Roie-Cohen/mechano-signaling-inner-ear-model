function  g = findHoppingIntercalations(g, area4hopping)
% initiates Hopping intercalation for
% area4hopping: if the cell area is less than this with respect to A0, it hops.
bond_length_thresh = 0.6*sqrt(g.A0);
HCs = g.LImodel.high_delta_cells;
for i=1:length(HCs)
    c = HCs(i);
    if g.linkedCells(c)~=0, continue; end
    AoA0 = cellarea(g,c)/g.areas(c);
    if AoA0 < area4hopping
        c_bonds = g.cells{c+1};
        if sum(ismember(g.bonds(c_bonds, 4), g.top_boundary_cells))~=0, continue; end
        for bi=1:length(c_bonds)
            b = c_bonds(bi);
            binv = find(g.bonds(:,1) == g.bonds(b,2) & g.bonds(:,2) == g.bonds(b,1), 1);
            neigh_bonds = g.cells{g.bonds(b,4)+1};
            binv_n_i = find(neigh_bonds==binv,1);
            b_out = neigh_bonds(mod(binv_n_i,length(neigh_bonds))+1);
            edge_cells = find_bond_edge_cells(g, b_out);
            edge_cells(edge_cells==c) = [];
            if ismember(edge_cells, g.top_boundary_cells)
                bl = bondLength(g, b_out);
                if bl > bond_length_thresh, continue; end
                vert = g.bonds(b_out, 2);
                g = HoppingIntercalation(g, c, vert);
                disp(['Hopping: cell ', num2str(c)]);
                break;
            end
        end
    end
end

end