function g = abolishLIregion(g, y_top, y_bottom, x_right, x_left)
% marks the abolished cells as type 10 (for above aboilished region) and type 20 (for below) 

nc = length(g.cells)-1;
pos = cellCOM(g, 1:nc);
removeLI = zeros(nc, 1);

if nargin>=2
    removeLI(pos(:,2) > y_top) = 1;
    g.type(pos(:,2) > y_top) = 10; 
end
if nargin>=3
    removeLI(pos(:,2) < y_bottom) = 1;
    g.type(pos(:,2) < y_bottom) = 20;
end
if nargin>=4
    removeLI(pos(:,1) > x_right) = 1;
end
if nargin==5
    removeLI(pos(:,1) < x_left) = 1;
end

cs = find(removeLI);
g = abolishLateralInhibition(g, cs);
g.LImodel.abolished = removeLI;

end