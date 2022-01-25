function dy = LIdiff(t,y,g)
% returns the differentials of the lateral inhibition model:
% the size of y is 3*nb+nc where nb, nc are the number of bonds and cells
% y(1:nb) = notch concentation on the bond
% y(nb+1:2*nb) = delta concentration on the bond
% y(2*nb+1:3*nb) = jagged concentration on the bond
% y(3*nb+1:3*nb+nc) = repressor amount in the cell

nb = length(g.bonds); % number of bonds
nc = length(g.cells)-1; % number of cells

% lateral inhibition model parameters:
Bn = g.LImodel.beta_n; % array of size nc
Bd = g.LImodel.beta_d; % array of size nc
Bj = g.LImodel.beta_j; % array of size nc
Br = g.LImodel.beta_r;
Gr = g.LImodel.gamma_r;
K = g.LImodel.kt;
l = g.LImodel.l;
m = g.LImodel.m;
n = g.LImodel.n;
tau = g.LImodel.tau; % time delay for repressor production

% initiating differential vector
dy = zeros(3*nb+nc,1);
bond_notch = y(1:nb); % g.LImodel.bond_notch
bond_delta = y(nb+1:2*nb); % g.LImodel.bond_delta
bond_jag = y(2*nb+1:3*nb); % g.LImodel.bond_jag
cell_repressor = y(3*nb+1:3*nb+nc); % g.LImodel.cell_repressor

% cell_fringe = g.LImodel.cell_fringe;

for c=1:nc
    if g.dead(c)==1, continue; end
    L = cellPerimeter(g, c); % the cell perimeter
    R = cell_repressor(c); % repressor amount in the cell
%     Fc = cell_fringe(c); % fringe level in the cell
    
    S = 0; % the total signal in the cell
    bonds = g.cells{c+1}; % the list of bonds in cell c
    for i=1:length(bonds)
        b = bonds(i);
        binv = find(g.bonds(:,1)==g.bonds(b,2) & g.bonds(:,2)==g.bonds(b,1),1); % inverse bond
%         fijn = (Fc)^n;
        
        % notch concentration differential
        nij = bond_notch(b); % notch concentration on the bond
        if isempty(binv)
            dji = 0; % delta concentration on the inverse bond
            jji = 0; % jagged concentration on the inverse bond
        else
            dji = bond_delta(binv); 
            jji = bond_jag(binv); 
        end
        dy(b) = Bn(c)/L - nij - K*nij*dji;
%         dy(b) = Bn(c)/L - nij - K*nij*dji*fijn/(1+fijn) - K*nij*jji/(1+fijn);
        
        
        % delta and jagged concentration differential
        if isempty(binv)
            nji = 0; % notch concentration on the inverse bond
        else
            nji = bond_notch(binv); 
        end
        dij = bond_delta(b); % delta concentration on the bond
        dy(b+nb) = (Bd(c)/L)/(1+R^l) - dij - K*nji*dij;
%         dy(b+nb) = (Bd(c)/L)/(1+R^l) - dij - K*nji*dij*fijn/(1+fijn);

%         jij = bond_jag(b); % jagged concentration on the bond
%         dy(2*b+nb) = (Bj(c)/L) - jij - K*nji*jij/(1+fijn);
        
        lij = bondLength(g, b); % the length of the bond
        S = S + lij*nij*dji;
%         S = S + lij*nij*(dji*fijn + jji)/(1+fijn);
    end
    
    % repressor differential
    hill = (g.globs.timer > tau)*(S^m)/(1+S^m);
    dy(3*nb+c) = Br*hill - Gr*R;
end

end