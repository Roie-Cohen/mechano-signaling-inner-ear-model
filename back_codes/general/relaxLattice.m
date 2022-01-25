function g = relaxLattice(g,n,vid, eps)
if nargin < 4, eps = 0.02; end
if nargin < 3, vid = []; end
for i=1:n
    g_prev = g;
    g.globs.timer = g.globs.timer + 1;
    g.A0 = average_area(g, 1:length(g.cells)-1);
    
    % lateral inhibition model
    eps_LI = 0.02; % 0.02
    g = relaxSignaling(g, eps_LI);
    
    % gradient decent
    ve = extractverts(g); 
    dE=denergy(g); 
    noE = norm(dE);

    if g.with_noise
        dE(dE~=0) = dE(dE~=0) + g.globs.noise(dE~=0)*0.2*noE;
        if mod(g.globs.timer, 30)
            g.globs.noise = rand(1, 2*length(g.verts)) - 0.5;
        end
    end

    if noE>1E-10
        
        dE = eps*dE/noE; 
        ve = ve-dE';
        g_prev = g;
        g = insertverts(ve,g);
        
        g = ReposeLatticeByPCs(g, g_prev);
         
        % update notch and delta concentration after geometry change
        diffusion_model = 1;
        g = handleBondConcentration(g, g_prev, diffusion_model);
    end
    

        if ~isempty(vid)
            % capturing a frame every nv steps
            nv = 3;
            if (mod(g.globs.timer,nv)==0 || g.globs.timer==1)
                figure(6);
                LatticePresentation(g,0, 6);
                frame = getframe;
                writeVideo(vid,frame);
            end
        end
    
end

end