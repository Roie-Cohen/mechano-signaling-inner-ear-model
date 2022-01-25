function g = updateParameters1(g)

if g.is_LImodel == 1
    
    %% sets the contractility and compressibility proportional to delta
    is_rigid_HCs = 1;
    if is_rigid_HCs
        cs = g.LImodel.high_delta_cells;
        for ci=1:length(cs)
            c = cs(ci);
            D = g.LImodel.cell_delta(c);
            lc = g.linkedCells(c);
            if lc~=0
                D = D + g.LImodel.cell_delta(lc);
            end
            g.ctrc(c) = min(max(1*D, 1), 2);
            g.cmpr(c) = min(max(2*D, 1), 5);
        end
    end
    %% sets adhesion rules
    bottomPCTension = g.tension(2,1); 
    topPCTension = g.tension(2,1); 
    PCPCTension = g.tension(2,2); 
    HCPCTension = g.tension(2,3); 
    HCSC_tension = g.tension(3,1); 
    HCHC_tension = g.tension(3,3); 
    
    g.tension(:, 5) = 1; 
    
    boundary_cells = [];
    totalTopNeighs = [];
    % sets high tension in the bonds on the top boundary
    cs = find(g.LImodel.abolished == 0);
    for ci=1:length(cs)
        c = cs(ci); % cell index
        bs = g.cells{c+1}; % the bonds of the cell
        neighs = g.bonds(bs, 4); % the neighbors of the cell
        topNeighs = neighs ~= 0; % logical array marking the neighbors of type "top boundary"
        for i=1:length(topNeighs)
            if topNeighs(i) == 1, topNeighs(i) = g.type(neighs(i))==10; end
        end
        for bi = 1:length(bs)
            b = bs(bi);
            neigh = neighs(bi); % the neighboring cell on the other side of b
            binv = find(g.bonds(:,4)==c & g.bonds(:,3)==neigh , 1);
            isTop = topNeighs(bi);
            if isTop
                if ismember(c, g.LImodel.high_delta_cells)
                    g.bonds([b, binv], 5) = HCPCTension;
                else
                    g.bonds([b, binv], 5) = bottomPCTension;
                end
            end
        end
        
        % updates the boundary cells vector
        if sum(topNeighs)>0
            boundary_cells = [boundary_cells, c];
            totalTopNeighs = [totalTopNeighs; neighs(topNeighs)];
        end
    end
    g.top_boundary_cells = unique(totalTopNeighs);
    
    % set top tension for PCS and PC:PC tension
    cs = g.top_boundary_cells;
    for ci=1:length(cs)
        c = cs(ci); % cell index
        bs = g.cells{c+1}; % the bonds of the cell
        neighs = g.bonds(bs, 4); % the neighbors of the cell
        topNeighs = neighs ~= 0; % logical array marking the neighbors of type "top boundary"
        for i=1:length(topNeighs)
            if topNeighs(i) == 1, topNeighs(i) = g.type(neighs(i)) == 10 & ~ismember(neighs(i), cs); end
        end
        PCneighs = ismember(neighs, cs);
        for bi = 1:length(bs)
            b = bs(bi);
            neigh = neighs(bi); % the neighboring cell on the other side of b
            binv = find(g.bonds(:,4)==c & g.bonds(:,3)==neigh , 1);
            isTop = topNeighs(bi);
            if isTop
                g.bonds([b, binv], 5) = topPCTension;
            end
            isPC = PCneighs(bi);
            if isPC
                g.bonds([b, binv], 5) = PCPCTension;
            end
        end
    end
    
    % set adhesion between HCs and SCs
    
    SCs = find(g.type == 1 & ~g.dead);
    PCs = g.top_boundary_cells;
    HCs = g.LImodel.high_delta_cells;
    SCs(ismember(SCs,PCs)) = [];
    SCs(ismember(SCs,HCs)) = [];
    
    for ci=1:length(SCs)
        c = SCs(ci);
        bs = g.cells{c+1}; % the bonds of the cell
        neighs = g.bonds(bs, 4); % the neighbors of the cell
        for bi=1:length(bs)
            neigh = neighs(bi);
            if ismember(neigh, HCs)
                b = bs(bi);
                binv = find(g.bonds(:,4)==c & g.bonds(:,3)==neigh , 1);
                g.bonds([b, binv], 5) = HCSC_tension;
            end
        end
    end

    % set adhesion between HCs and HCs
    for ci=1:length(HCs)
        c = HCs(ci);
        bs = g.cells{c+1}; % the bonds of the cell
        neighs = g.bonds(bs, 4); % the neighbors of the cell
        for bi=1:length(bs)
            neigh = neighs(bi);
            if ismember(neigh, HCs)
                b = bs(bi);
                binv = find(g.bonds(:,4)==c & g.bonds(:,3)==neigh , 1);
                g.bonds([b, binv], 5) = HCHC_tension;
            end
        end
    end
    
end

end