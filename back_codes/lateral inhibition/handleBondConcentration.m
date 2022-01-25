function g = handleBondConcentration(g, g_prev, is_diffusion_model)
% updates notch and delta bond concentrations after a geometry change
% is_diffusion_model = 1 for a model where concentration is instantly
% diffused over all cell circumference.

if isfield(g,'LImodel')
    
    if is_diffusion_model
        for c=1:length(g.cells)-1
            if g.dead(c) || g.linkedCells(c)~=0 , continue; end
            bs = g.cells{c+1}; % the bonds of the cell
            L = cellPerimeter(g, c);
            g.LImodel.bond_notch(bs) = g.LImodel.cell_notch(c)/L;
            g.LImodel.bond_delta(bs) = g.LImodel.cell_delta(c)/L;
            g.LImodel.bond_jag(bs) = g.LImodel.cell_jag(c)/L;
        end
        lc = find(g.linkedCells~=0);
        while ~isempty(lc)
            c1 = lc(1); 
            c2 = g.linkedCells(c1);
            L = cellPerimeter(g, c1) + cellPerimeter(g, c2);
            N = g.LImodel.cell_notch(c1) + g.LImodel.cell_notch(c2);
            D = g.LImodel.cell_delta(c1) + g.LImodel.cell_delta(c2);
            J = g.LImodel.cell_jag(c1) + g.LImodel.cell_jag(c2);
            bs1 = g.cells{c1+1};
            g.LImodel.bond_notch(bs1) = N/L;
            g.LImodel.bond_delta(bs1) = D/L;
            g.LImodel.bond_jag(bs1) = J/L;
            bs2 = g.cells{c2+1};
            g.LImodel.bond_notch(bs2) = N/L;
            g.LImodel.bond_delta(bs2) = D/L;
            g.LImodel.bond_jag(bs2) = J/L;
            lc(lc==c1 | lc==c2) = [];
        end
    else
        nb = length(g.bonds);
        bs = bondLength(g, 1:nb);
        bs_prev = bondLength(g_prev, 1:nb);
        cond = bs~=-1 & bs_prev~=-1;
        g.LImodel.bond_notch(cond) = g.LImodel.bond_notch(cond).*bs_prev(cond)./bs(cond);
        g.LImodel.bond_delta(cond) = g.LImodel.bond_delta(cond).*bs_prev(cond)./bs(cond);
        g.LImodel.bond_jag(cond) = g.LImodel.bond_jag(cond).*bs_prev(cond)./bs(cond);
    end
end


end