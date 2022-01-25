function g = killLosingHCs(g)
% if a hair cell does not touch the PC row for a long time it delaminates

period_thresh = 250; % the time after which a lonely HC will die
HCs = g.LImodel.high_delta_cells;
PCs = g.top_boundary_cells;
DH = g.LImodel.delta_history;
for i=1:length(HCs)
    c = HCs(i);
    neighs = g.bonds(g.cells{c+1}, 4);
    if sum(ismember(neighs, PCs))==0
        c_history = DH(:, c);
        risk_period = g.globs.timer - find(c_history ~= 1, 1, 'last');
        if risk_period > period_thresh
            g.ctrc(c) = 6;
            g.cmpr(c) = 0.5;
            g.LImodel.beta_d(c) = 0;
            g.LImodel.cell_delta(c) = 0;
            g.LImodel.bond_delta(g.cells{c+1}) = 0;
            g.LImodel.high_delta_cells(g.LImodel.high_delta_cells == c) = [];
        end
    end
end

end