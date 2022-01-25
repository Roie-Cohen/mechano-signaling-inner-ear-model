function g = handleVirtualVertices(g)
% adds or removes virtual vertices
% 1) a virtual vertex close to a real vertex (<Lmin) is removed.
% 2) a small bond (<Lmin) created by two virtual vertices is removed by
% averaging the two vertices into one.
% 3) a large bond (>Lmax) is splited by adding a virtual vertex in the middle.
if ~g.is_multiVM, return; end

Lmax = g.virtual_parameters.Lmax;
Lmin = g.virtual_parameters.Lmin;

nb = length(g.bverts);
treated_bonds = zeros(nb, 1, 'logical');
for b=1:nb
    rverts = g.bonds(b, 1:2); % the two real vertices
    if rverts(1)==0, continue; end
    
    % check if inverse bond was treated
    b_rev = find(g.bonds(:,1)==rverts(2) & g.bonds(:,2)==rverts(1));
    if treated_bonds(b_rev)
        g.bverts{b} = flip(g.bverts{b_rev});
        treated_bonds(b) = 1;
        continue;
    end
        
    vverts = g.bverts{b}; % virtual vertices
    verts = [rverts(1), vverts, rverts(2)];
    nv = length(verts);
    vpos = getRelativePosition(g, verts);
    offset = 0;
    for i=1:nv-1
        dr = vpos(i+1,:) - vpos(i,:);
        L = norm(dr);
        
        if L>Lmax
            % add new virtual vertex in the middle
            new_pos = vpos(i,:) + dr/2;
            new_vind = length(g.verts)+1;
            g.verts(new_vind, 1:2) = new_pos;
            vverts = [vverts(1:i+offset-1), new_vind, vverts(i+offset:end)];
            offset = offset + 1;
        elseif L<Lmin
            if isempty(vverts)
                treated_bonds(b)=1;
                continue;
            end
            % remove a virtual vertex
            switch i
                case 1 % the first vertex is real
                    % delete the virtual vertex 
                    rem_vind = 1;
                case nv-1 % the last vertex is real
                    % delete the virtual vertex
                    rem_vind = nv-2 + offset;
                otherwise
                    % average the two virtual vertices into one
                    rem_vind = i + offset;
                    g.verts(verts(i), 1:2) = vpos(i,:) + dr/2;
            end
            g.verts(vverts(rem_vind), :) = 0;
            vverts(rem_vind) = [];
            offset = offset - 1;
        end
    end
    g.bverts{b} = vverts;
    treated_bonds(b) = 1;
end

end