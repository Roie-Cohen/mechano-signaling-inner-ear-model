function g = initLateralInhibitionModel(g)
% Expanded model from Udi's paper. 
% Includes Jagged1 and Lfng effects.

nb = length(g.bonds); % number of bonds
nc = length(g.cells)-1; % number of cells

% lateral inhibition model parameters: (taken from Udi's paper)
g.LImodel.beta_n = 3.9*ones(nc, 1);
g.LImodel.beta_d = 3.9*ones(nc, 1);
g.LImodel.beta_j = 1*ones(nc, 1);
g.LImodel.beta_r = 194; %194.1
g.LImodel.gamma_r = 60; % 1
g.LImodel.kt = 0.2;
g.LImodel.l = 3; % 3
g.LImodel.m = 3;
g.LImodel.n = 3;
g.LImodel.tau = 0; % time delay for repressor production

g.LImodel.high_delta_thresh = 2;
g.LImodel.high_delta_cells = [];
g.LImodel.prev_delta_cells = []; % cells that were high delta at some point
g.LImodel.delta_history = zeros(1, nc);

% initiating components
bond_notch = zeros(nb, 1);
bond_delta = zeros(nb, 1);
bond_jag = zeros(nb, 1);

% initial total notch and ligands
N0 = 2;
D0 = 0.1;
J0 = 1;

% signal to noise ratio
N_snr = 0.1;
D_snr = 0.1;
J_snr = 0.1;

% total notch and ligands
N = N0*awgn(ones(nc,1), N_snr);
D = D0*awgn(ones(nc,1), D_snr);
J = J0*awgn(ones(nc,1), J_snr);
N(N<0)=0; D(D<0)=0; J(J<0)=0;

for c = 1:nc
    bs = g.cells{c+1}; % bonds of cell c
    L = cellPerimeter(g,c); % cell perimeter
    bond_notch(bs) = N(c)/L;
    bond_delta(bs) = D(c)/L;
    bond_jag(bs) = J(c)/L;
end

g.LImodel.bond_notch = bond_notch;
g.LImodel.bond_delta = bond_delta;
g.LImodel.bond_jag = bond_jag;
g.LImodel.cell_notch = N;
g.LImodel.cell_delta = D;
g.LImodel.cell_jag = J;

% initialize repressor as zero
g.LImodel.cell_repressor = zeros(nc, 1);

% initiates Fringe level to one
g.LImodel.cell_fringe = ones(nc, 1);
end