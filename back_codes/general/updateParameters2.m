function g = updateParameters2(g)

if g.is_LImodel == 1
    
    for b=1:length(g.bonds)
        if g.bonds(b,1) ==0, continue; end
        neighs = g.bonds(b, 3:4);
        typs = g.type(neighs);
        typs(typs==10 | typs==20) = 4;
        g.bonds(b, 5) = g.tension(typs(1), typs(2));
    end
    
end

end