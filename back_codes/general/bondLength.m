function l = bondLength(g, bs)
% returns the bond length of bs
% if there is an error with a bond the function returns -1

l = zeros(length(bs),1); % initializing bond length
for bi=1:length(bs)
    b = bs(bi);
    vidx = [g.bonds(b,1), g.bverts{b}, g.bonds(b,2)];
    if vidx(1)==0
        l(b) = -1;
        continue;
    end
    vpos = getRelativePosition(g,vidx);
    nv = length(vidx);
    for i=1:nv-1
        n = norm(vpos(i,:)-vpos(i+1,:));
        l(bi) = l(bi)+n;
    end
end

end