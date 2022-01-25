function g = insertNewCell(g, new_c, vert)
% insets a small triangular cell "c" into the vertex "vert"

if new_c<=length(g.cells)-1
    if ~isdead(new_c)
        disp("Error in insertNewCell");
        return;
    end
end

% find crossed bonds and cells
bs = find(g.bonds(:,2) == vert);
cs = g.bonds(bs, 3);
bsp = zeros(size(bs));
for i=1:length(bs)
    b = bs(i);
    c = cs(i);
    c_bonds = g.cells{c+1};
    bi_c = find(c_bonds == b,1);
    bsp(i) = c_bonds( mod(bi_c, length(c_bonds)) + 1 );
end

% reorder so that bs are counter-clockwise
if g.bonds(bsp(1), 4) ~= cs(2)
    % switch 2 and 3
    bs([2 3]) = bs([3 2]);
    bsp([2 3]) = bsp([3 2]);
    cs([2 3]) = cs([3 2]);
end

% create new vertices in the position of vert
nverts = length(g.verts);
vs = [1:3] + nverts;
g.verts(vs, 1) = g.verts(vert,1);
g.verts(vs, 2) = g.verts(vert,2);

% slightly move each new vertex in the direction of one bond
eps = 0.3; % ratio of bond length
rel_pos = getRelativePosition(g, [g.bonds(bs, 1); vert]);
for i=1:length(vs)
    v = vs(i);
    p_far = rel_pos(i,:);
    p_close = rel_pos(4,:);
    dr = eps*(p_far - p_close);
    g.verts(v,1:2) = g.verts(v,1:2) + dr;
end

% update bs and bsp and add new bonds
bspp = [1:3] + length(g.bonds);
b4s = [4:6] + length(g.bonds);
g.cells{new_c+1} = flip(b4s);
for i=1:length(vs)
    v = vs(i);
    v_next = vs(mod(i,length(vs))+1);
    b = bs(i);
    bp = bsp(i);
    bpp = bspp(i);
    c = cs(i);
    b4 = b4s(i);
    
    % update existing bonds
    g.bonds(b, 2) = v;
    g.bonds(bp, 1) = v_next;
    
    % add new bonds
    g.bonds(bpp, 1:4) = [v, v_next, c, new_c];
    g.bonds(b4, 1:4) = [v_next, v, new_c, c];
    
    % add bpp to g.cells
    c_bonds = g.cells{c+1};
    bi_c = find(c_bonds == b,1);
    g.cells{c+1} = [c_bonds(1:bi_c), bpp, c_bonds(bi_c+1:end)];
    
    % add new virtual bonds
    if isfield(g,"bverts")
        g.bverts{bpp} = [];
        g.bverts{b4} = [];
    end
end

% override "vert" in g.verts
g.verts(vert, :) = 0;

end