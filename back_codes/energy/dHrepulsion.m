function dE=dHrepulsion(g,c)
dE = 0;
if g.is_LImodel
    iHCs = g.LImodel.high_delta_cells;
    if ~ismember(c, iHCs), return; end % if this is not a HC
    
    dE = zeros(2*length(g.verts),1);
    % parameters for repulsion potential
    kappa = g.kappa; 
    sigma = 2*sqrt(g.A0/pi);
    dErep = @(r) max(-kappa*sigma^kappa./(r.^(kappa+1)), -100); % limiting the energy for the digit accuracy
    
    vidx = getVerts(g,c);
    vert = getRelativePosition(g,vidx);
    vn = length(vidx); % number of vertices in the cell
    Cn = g.centroid(c,:);
    
    iHCs(iHCs == c) = [];  % removing the current cell
    iHCs(iHCs == g.linkedCells(c)) = []; % removing the linked cell (if exists)
    nHCs = length(iHCs);
    A = cellarea(g, c);
    for j = 1:nHCs
        m = iHCs(j);    % HC index
        Cm = g.centroid(m, :);
        r = Cn - Cm;
        
        % if we're with a periodic boundary conditions we need to check
        % if the distance between 'n' and 'm' is actually closer
        if g.bc == 1
            % the scaling of the lattice in periodic BC is 2pi, meaning
            % the width of the lattice is 2pi.
            % if there's a closer matching vertex in the periodic
            % lattice, then r is replaced with the vector pointing from
            % the closer point.
            r(1) = r(1)-(2*pi-abs(r(1)) < abs(r(1)))*sign(r(1))*2*pi;
        end
        if norm(r) > 1 % cutoff distance
            continue;
        end
        % dEnm = dE/dr * dr/dx
        dr = r/norm(r);
        dEnmdr = dErep(norm(r));
        for k=1:vn
            kn = mod(k,vn)+1; % next vertex
            kp = mod(k-2,vn)+1; % previous vertex
            r0 = vert(kp,:); r1 = vert(k,:); r2 = vert(kn,:); % positions of vertices k-1, k, k+1 respectively
            xp = r0(1); x = r1(1); xn = r2(1); yp = r0(2); y = r1(2); yn = r2(2);
            dAdx = 0.5*( yn-yp );
            dAdy = 0.5*( xp-xn );
            dXcmdx = ( -dAdx*Cn(1) + (1/6)*(2*x*(yn-yp)+y*(xp-xn)+xn*yn-xp*yp) )/A;
            dYcmdx = ( -dAdx*Cn(2) + (1/6)*(yn*(y+yn)-yp*(y+yp)) )/A;
            dXcmdy = ( -dAdy*Cn(1) + (1/6)*(xp*(x+xp)-xn*(x+xn)) )/A;
            dYcmdy = ( -dAdy*Cn(2) + (1/6)*(2*y*(xp-xn)+x*(yn-yp)+yp*xp-yn*xn) )/A;
            dEnm = dEnmdr*[dr(1)*dXcmdx+dr(2)*dYcmdx, dr(1)*dXcmdy+dr(2)*dYcmdy];
            dE(2*vidx(k)-1) = dE(2*vidx(k)-1) - dEnm(1);
            dE(2*vidx(k)) = dE(2*vidx(k)) - dEnm(2);
        end
    end
    dE = g.paras(4)*dE;
    
end

end