function g = HoppingIntercalation(g, c, vert)
% creates a linked cell for cell 'c' in vertex 'vert'

gi = g;
% create a new cell
new_c = length(g.cells);
g = insertNewCell(g, new_c, vert);

% link new cell to 'c'
if ~isfield(g,'linkedCells')
    g.linkedCells = zeros(new_c, 1);
end
g.linkedCells(c) = new_c;
g.linkedCells(new_c) = c;

% initiate parameters for the linked new cell
g.dead(new_c) = 0;
g.type(new_c) = g.type(c);
g.areas(new_c) = g.areas(c);
g.cmpr(new_c) = g.cmpr(c);
g.ctrc(new_c) = g.ctrc(c);

if isfield(g,'LImodel')
    g.LImodel.abolished(new_c) = g.LImodel.abolished(c);
    g.LImodel.beta_n(new_c) = g.LImodel.beta_n(c);
    g.LImodel.beta_d(new_c) = g.LImodel.beta_d(c);
    g.LImodel.beta_j(new_c) = g.LImodel.beta_j(c);
    nb = length(g.bonds);
    g.LImodel.bond_notch(nb-5:nb) = 0;
    g.LImodel.bond_delta(nb-5:nb) = 0;
    g.LImodel.bond_jag(nb-5:nb) = 0;
    g.LImodel.cell_notch(new_c) = 0;
    g.LImodel.cell_delta(new_c) = 0;
    g.LImodel.cell_jag(new_c) = 0;
    g.LImodel.cell_repressor(new_c) = 0;
    g = handleBondConcentration(g, gi, 1);
end

g.relaxedCells = [g.relaxedCells; new_c, g.globs.timer];

end