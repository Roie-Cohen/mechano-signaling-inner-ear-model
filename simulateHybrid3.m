function g = simulateHybrid3(ini_latpath, latnum, with_video)
% runs hybrid3 model starting from the final point in a hybrid2 lattice
% latnum = the ID of the lattice
% set with_video=1 if you want a movie of the simulation
%% general
warning off
addpath(genpath(fileparts(which(mfilename))));

if nargin <= 2, with_video = 1; end
if nargin <= 1, latnum = 1; end
if nargin == 0, ini_latpath = ['lattices\hybrid2\lat_', num2str(latnum),'_final.mat']; end

load(ini_latpath, 'g');

% define cell types
PCs = g.top_boundary_cells;
PCs(logical(g.dead(PCs))) = [];
g.type(PCs) = 2;

HCs = g.LImodel.high_delta_cells;
HCs(logical(g.dead(HCs))) = [];
g.type(HCs) = 3;

% define SCs as non-HCs that touch PCs
SCs = find(~g.LImodel.abolished);
SCs(ismember(SCs, HCs)) = [];
g.type(SCs) = 1;
for j=1:length(SCs)
    c = SCs(j);
    neighs = g.bonds(g.cells{c+1}, 4);
    if ~ismember(2, g.type(neighs))
        g = abolishLateralInhibition(g, c);
        g.LImodel.abolished(c) = 1;
        g.type(c) = 20;
    end
end


% types: 1: SC, 2: PC, 3: HC, 4: other
g.tension = [1  10  1  1 ; ... % tension of bond seperating populations n and m is g.tension(n,m)
             10  1  7  10; ...
             1   7  10 1; ...
             1  10  1  1 ];
         

% simulation parameters and global parameters
timer = 0;
g.with_noise = 0;
noise = rand(1, 2*length(g.verts)) - 0.5;
g.globs = struct('timer',timer,'noise',noise);

T1prob = 1;
T1eps = 0.3*sqrt(g.A0/pi);
g.transitionedBonds = [0, 0];
g.relaxedCells = [0, 0];

% video
if with_video
    sim_folder = 'simulations\hybrid3\';
    vid = VideoWriter([sim_folder, 'sim_hybrid3_', num2str(latnum),'.avi']);
%     vid.Quality = 100; % 100 for best
    vid.FrameRate = 12;
    open(vid);
else
    vid = [];
end

% weights: incompessibility, tension, contractility, HC-HC repulsion,
% HC pull to PC
g.paras = [1 ; 0.05; 0.1; 1; 1; 0];

%% run simulation
g = updateParameters2(g);
ts = 30; t=1;
while t<= ts
    savename = ['lattices\hybrid3\lat_', num2str(latnum),'_t=', num2str(t-1),'.mat'];
    save(savename, 'g');
    
    g = relaxLattice(g,3,vid);
    g = findTransitions(g,T1eps,T1prob,0.02);
    g = updateParameters2(g);
    g = relaxLattice(g,3, vid);
    g = killSmallCells(g,0.1);
    g = killLosingHCs(g);
    g = updateParameters2(g);
    g = relaxLattice(g,3, vid);
    g = handleVirtualVertices(g);
    g = findHoppingIntercalations(g, 0.9);
    
    disp(t)
    t = t+1;
    
end
close(vid);

savename = ['lattices\hybrid3\lat_', num2str(latnum),'_final.mat'];
save(savename, 'g');

end