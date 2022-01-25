function g = convert2multi(g)
% gets a 'g' structure of a lattice in the format of a vertex model
% returns a generalized structure with virtual vertices

% initialize parameters:
% Lmax = maximal bond length above which the bond breaks into two bonds
% Lmin = minimal bond length below which the bond the bond vanishes
Lmax = 0.5*sqrt(g.A0/pi);
Lmin = 0.1*sqrt(g.A0/pi);
g.virtual_parameters.Lmax = Lmax;
g.virtual_parameters.Lmin = Lmin;

% a cell array that containd the virtual vertices indices of each bond
% the vertices are ordered  in the direction from the first vertex of the
% bond (g.bonds(i,1)) to the second vertex (g.bonds(i,2))
g.bverts = cell(length(g.bonds), 1);
treated_bonds = zeros(length(g.bonds), 1, 'logical');
for i=1:length(g.cells)-1
    if g.dead(i), continue; end
    bs = g.cells{i+1}; % bonds of the current cell
    vidx = g.bonds(bs,1); % an array of the vertices indices of the cell
    vert = getRelativePosition(g,vidx); % the position of the vertices
    nv = length(vidx);
    for j=1:nv
        b = bs(j); % the current bond
        next = mod(j,nv)+1;   % the next vertex index in "vert"
        
        % checks if the reveresed bond has already been treated
        b_rev = find(g.bonds(:,1)==vidx(next) & g.bonds(:,2)==vidx(j));
        if treated_bonds(b_rev)
            g.bverts{b} = flip(g.bverts{b_rev});
            treated_bonds(b) = 1;
            continue;
        end
        
        % if the reversed bond has not been treated yet, continue with the
        % current bond
        L = norm(vert(j,:)-vert(next,:)); % bond length
        if L>2*Lmin 
            nvv = ceil(mean([ceil(L/Lmax),floor(L/Lmin)]))-1; % number of virtual vertices
            new_verts = [1:nvv] + length(g.verts); % new vertices indeced
            g.bverts{b} = new_verts;
            
            % update the position of the new verts in g.verts
            % divide the line that passes between 'v1' to 'v2' into 'nvv' break points.
            % the break points are ordered in the direction from 'v1' to 'v2'.
            v1 = vert(j,:);
            v2 = vert(next,:);
            vv_pos = [linspace(v1(1), v2(1), nvv+2)',linspace(v1(2), v2(2), nvv+2)']; 
            vv_pos(1, :) = []; vv_pos(end, :) = [];  % removes the two original vertices
            for k=1:nvv
                g.verts(new_verts(k), 1:2) = vv_pos(k, :);
            end
        end
        treated_bonds(b) = 1;
    end
end

end