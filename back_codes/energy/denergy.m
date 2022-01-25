 function dE = denergy(g)

g = updateCellsCenter(g);
dE = zeros(2*length(g.verts),1);
for c=1:length(g.cells)-1
    if(g.dead(c) == 0)
        
        Earea = dHarea(g,c); % area term (compressibility)
        
        Elength = dHlength(g,c); % tension and contractility terms
        
        if g.paras(4)~=0 % HC repulsion term
            Erep = dHrepulsion(g,c);
        else
            Erep = 0;
        end
        
        if g.paras(5) ~=0 % lateral compression term
            Epcs = dHpcs(g,c);
        else
            Epcs = 0;
        end
        
        dEi = Earea + Elength + Erep + Epcs;
        dE = dE+dEi;

    end
end

dE=dE';
end