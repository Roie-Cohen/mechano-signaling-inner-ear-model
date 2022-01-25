function [ng, torem] = T2transition(g,ce)
ng = g;
if(length(g.cells{ce+1})<=3)
    rc = find(g.relaxedCells(:,1)==ce);
    if ~isempty(rc)
        rc = rc(end);
        if g.globs.timer < g.relaxedCells(rc, 2) + 100
            return;
        end
    end
    
    disp(['removing cell ' num2str(ce)]);
    
    bo = ng.cells{ce+1};
    verts = ng.bonds(bo,1);
    pos = getRelativePosition(ng,verts);
    v = mean(pos);
    ng.verts(verts(1),1:2) = v;
    % replace vert(2) and vert(3) with vert(1)
    ng.bonds(ng.bonds(:,1) == verts(2),1)=verts(1);
    ng.bonds(ng.bonds(:,2) == verts(2),2)=verts(1);
    if(length(bo)==3)
        ng.bonds(ng.bonds(:,1) == verts(3),1)=verts(1);
        ng.bonds(ng.bonds(:,2) == verts(3),2)=verts(1);
    end
    binv = zeros(1, length(bo));
    for i=1:length(bo)
        bi = find(g.bonds(:,1) == g.bonds(bo(i),2) & g.bonds(:,2) == g.bonds(bo(i),1));
        if ~isempty(bi)
            bidx  =1;
            if(length(bi)>1) %% take the right one if many exist
                bidx = find(g.bonds(bi,4)==ce);
            end
            binv(i) =bi(bidx);
        end
    end
    torem = [bo binv(binv~=0)];
    ng = remove_bond(ng,torem);
    ng.cells{ce+1} = [];
    ng.dead(ce)=1;
    if isfield(ng,'linkedCells')
        ng.linkedCells(ce) = 0;
        ng.linkedCells(ng.linkedCells==ce) = 0;
    end
end