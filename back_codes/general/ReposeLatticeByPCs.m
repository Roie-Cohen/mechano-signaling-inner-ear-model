function g = ReposeLatticeByPCs(g, g_prev)
% repositions the lattice according to previous position of the PCs
if isfield(g, 'top_boundary_cells')
    
    y_only = 1;
    
    PCs = g.top_boundary_cells;
    if isempty(PCs), return; end
    PCs_pos = cellCOM(g, PCs);
    current_pos = mean(PCs_pos,1);
    
    % previous PCs position
    prev_PCs = g_prev.top_boundary_cells;
    if isempty(prev_PCs), return; end
    prev_PCs_pos = cellCOM(g_prev, prev_PCs);
    prev_pos = mean(prev_PCs_pos,1);
    
    dpos = prev_pos - current_pos;
    
    g.verts(:,2) = g.verts(:,2) + dpos(2);
    if ~y_only
        g.verts(:,1) = g.verts(:,1) + dpos(1);
    end

end
end