function g=findTransitions(g,eps,probT1,probT2)
%% generates T1 transition with probability  probT1 if edge is  shorter than eps
%% if cell has less than 4  edges,  generates a T2 transition is probability probT2
cand = randperm(size(g.bonds,1)); %% random  selection order
while ~isempty(cand)
    bo = cand(1);
    if (g.bonds(bo,1)==0 || g.bonds(bo,3)==0)
        cand = cand(2:end);
        continue;
    end
    
    v1 = g.verts(g.bonds(bo,1),:);
    v2 = g.verts(g.bonds(bo,2),:);
    vec = v1-v2;
    
    c1 = g.bonds(bo,3);
    c2 = g.bonds(bo,4);
    if(length(g.cells{c1+1})>3 && length(g.cells{c2+1})>3)
        len = norm(vec);
        if len<eps
            if rand() < probT1
                g = T1transition(g,bo);
            end
        end
    else
        c1 = g.bonds(bo,3);
        c2 = g.bonds(bo,4);
%         torem = [];
        if(length(g.cells{c1+1})<=3) % 3
            if(rand()<probT2)
                g=T2transition(g,c1);
%                 g = redistributeAreas(g);
            end
        end
        if c2~=0 && length(g.cells{c2+1})<=3 % 3
            if(rand()<probT2)
                g=T2transition(g,c2);
%                 g = redistributeAreas(g);
            end
        end
    end
    cand = cand(2:end);
end

g = findNeighboringLinkedCells(g);

