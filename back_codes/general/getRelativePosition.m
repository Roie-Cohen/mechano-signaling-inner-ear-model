function pos = getRelativePosition(g,v)
% gets a list of vertices (indices), can be both real and virtual.
% returns the reltive position of the vertices with respect to each other
% (as one cluster and not seperated by boundary periodicity)
pos = g.verts(v,1:2);
if(g.bc>=1) % relavent only for periodic boundary conditions
    for d=1:2
        ap = g.verts(v,d); % x or y position of the verts
        [~, pid] = min(abs(ap));
        p = ap(pid); % the x\y position of the vertex who's farthest from the boundary
        idx = mod(pid,length(v))+1; % next vertex index in 'v'
        d_period = g.latticeDims(d);
        while (idx~=pid)
            if(abs(pos(idx,d)-p + d_period)< abs(pos(idx,d)-p))
                pos(idx,d) =  pos(idx,d) + d_period;
            end
            if(abs(pos(idx,d)-p - d_period)< abs(pos(idx,d)-p))
                pos(idx,d) =  pos(idx,d) - d_period;
            end
            p = pos(idx,d);
            idx = mod(idx,length(v))+1; % gets the next vertex index
            %[idx pid]
        end
        
    end
end