function g = mergeNeighboringCells(g, c1, c2)
% Merges the neighboring cells c1 and c2
% The merged cell index is c1, and c2 is considered as dead

c1_bonds = g.cells{c1+1};
c2_bonds = g.cells{c2+1};

nb1 = length(c1_bonds);
nb2 = length(c2_bonds);

b1_i = find(g.bonds(c1_bonds, 4)==c2,1);
b2_i = find(g.bonds(c2_bonds, 4)==c1,1);

b1 = c1_bonds(b1_i);
b2 = c2_bonds(b2_i);

b1p = c1_bonds(mod(b1_i-2,nb1)+1);
b2p = c2_bonds(mod(b2_i-2,nb2)+1);
b1n = c1_bonds(mod(b1_i,nb1)+1);
b2n = c2_bonds(mod(b2_i,nb2)+1);

c3 = g.bonds(b1p, 4);
c4 = g.bonds(b2p, 4);

v1 = g.bonds(b1, 1);
v2 = g.bonds(b2, 1);
v13 = g.bonds(b1p, 1);
v14 = g.bonds(b1n, 2);
v24 = g.bonds(b2p, 1);
v23 = g.bonds(b2n, 2);

b3p = find(g.bonds(:,1)==v23 & g.bonds(:,2)==v1,1);
b3n = find(g.bonds(:,1)==v1 & g.bonds(:,2)==v13,1);
b4p = find(g.bonds(:,1)==v14 & g.bonds(:,2)==v2,1);
b4n = find(g.bonds(:,1)==v2 & g.bonds(:,2)==v24,1);

% update merged cell
g.cells{c1+1} = [c1_bonds(b1_i+1:nb1), c1_bonds(1:b1_i-1), c2_bonds(b2_i+1:nb2), c2_bonds(1:b2_i-1)];
g.cells{c2+1} = [];
g.cells{c1+1}(g.cells{c1+1}==b1n | g.cells{c1+1}==b2n) = [];
g.bonds(g.cells{c1+1}, 3) = c1;

% update bonds
g.bonds(b1, :) = 0;
g.bonds(b2, :) = 0;
g.bverts{b1} = [];
g.bverts{b2} = [];

g.bonds(b1n, :) = 0;
g.bonds(b2n, :) = 0;
g.bonds(b1p, 2) = v23;
g.bonds(b2p, 2) = v14;
g.bonds(b2p, 3) = c1;
g.bverts{b1p} = [g.bverts{b1p}, v1, g.bverts{b2n}];
g.bverts{b2n} = [];
g.bverts{b2p} = [g.bverts{b2p}, v2, g.bverts{b1n}];
g.bverts{b1n} = [];

g.bonds(b3n, :) = 0;
g.cells{c3+1}(g.cells{c3+1}==b3n) = [];
g.bonds(b4n, :) = 0;
g.cells{c4+1}(g.cells{c4+1}==b4n) = [];
g.bonds(b3p, 2) = v13;
g.bonds(b4p, 2) = v24;
g.bonds(b3p, 4) = c1;
g.bverts{b3p} = [g.bverts{b3p}, v1, g.bverts{b3n}];
g.bverts{b3n} = [];
g.bverts{b4p} = [g.bverts{b4p}, v2, g.bverts{b4n}];
g.bverts{b4n} = [];


g.bonds(g.bonds(:,4)==c2, 4) = c1;

g.dead(c2) = 1;

end