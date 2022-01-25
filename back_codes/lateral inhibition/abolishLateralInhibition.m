function g = abolishLateralInhibition(g, cs)
% abolishes notch-delta signaling for cell 'cs'

g.LImodel.beta_n(cs) = 0;
g.LImodel.beta_d(cs) = 0;
g.LImodel.beta_j(cs) = 0;

for ci=1:length(cs)
    c = cs(ci);
    bs = g.cells{c+1};
    g.LImodel.bond_notch(bs) = 0;
    g.LImodel.bond_delta(bs) = 0;
    g.LImodel.bond_jag(bs) = 0;
    g.LImodel.cell_notch(c) = 0;
    g.LImodel.cell_delta(c) = 0;
    g.LImodel.cell_jag(c) = 0;
end

end