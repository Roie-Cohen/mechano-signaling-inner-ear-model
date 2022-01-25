function g = relaxSignaling(g, epsilon)
if isfield(g,'LImodel')
    g = LIinitiateNewBonds(g);
    g = setCellDeltaLevel(g, find(g.type==20), 0.2); % sets basal delta level for GER
%     g = setCellNotchLevel(g, find(g.LImodel.abolished), mean(g.LImodel.cell_notch(~g.LImodel.abolished))); 
    y = [g.LImodel.bond_notch; g.LImodel.bond_delta; g.LImodel.bond_jag; g.LImodel.cell_repressor];
    
    % advancing the signaling with a small step. can use a ode solver instead
    dy = LIdiff([],y,g);
    y = y + epsilon*dy;
    
    nb = length(g.bonds);
    nc = length(g.cells)-1;
    g.LImodel.bond_notch = y(1:nb);
    g.LImodel.bond_delta = y(nb+1:2*nb);
    g.LImodel.bond_jag = y(2*nb+1:3*nb);
    g.LImodel.cell_repressor = y(3*nb+1:3*nb+nc);
    
    % fix false negatives due to overshoot
    g.LImodel.bond_notch(g.LImodel.bond_notch<0) = 0;
    g.LImodel.bond_delta(g.LImodel.bond_delta<0) = 0;
    g.LImodel.bond_jag(g.LImodel.bond_jag<0) = 0;
    g.LImodel.cell_repressor(g.LImodel.cell_repressor<0) = 0;
    
    Ns = zeros(nc,1);
    Ds = zeros(nc,1);
    Js = zeros(nc,1);
    for c = 1:nc
        Ns(c) = sum(g.LImodel.bond_notch(g.cells{c+1}).*bondLength(g, g.cells{c+1}));
        Ds(c) = sum(g.LImodel.bond_delta(g.cells{c+1}).*bondLength(g, g.cells{c+1}));
        Js(c) = sum(g.LImodel.bond_jag(g.cells{c+1}).*bondLength(g, g.cells{c+1}));
    end
    
    g.LImodel.cell_notch = Ns;
    g.LImodel.cell_delta = Ds;
    g.LImodel.cell_jag = Js;
    
    g.LImodel.high_delta_cells = find(Ds > g.LImodel.high_delta_thresh);
    g.LImodel.prev_delta_cells = unique([g.LImodel.prev_delta_cells; g.LImodel.high_delta_cells]);
    g.LImodel.delta_history(g.globs.timer, Ds >= g.LImodel.high_delta_thresh) = 1;
    g.LImodel.delta_history(g.globs.timer, Ds > 0.5*g.LImodel.high_delta_thresh & Ds < g.LImodel.high_delta_thresh) = 0.5;
end
end