function g = simulateHybrid2(ini_latpath, latnum, with_video)
% runs hybrid2 model starting from an initial random lattice
% latnum = the ID of the lattice
% set with_video=1 if you want a movie of the simulation
%% general
warning off
addpath(genpath(fileparts(which(mfilename))));

if nargin <= 2, with_video = 1; end
if nargin <= 1, latnum = 1; end
if nargin == 0, ini_latpath = 'lattices\random\rand_lat_10x25_1.mat'; end

load(ini_latpath, 'g');
nc = length(g.cells)-1;
% setting parameters
g.A0 = average_area(g, 1:nc);
g.areas(:) = g.A0;
g.type(:) = 1;
g.tension = [1  1  1  1 ; ... % tension of bond seperating populations n and m is g.tension(n,m)
             1  1  1  1; ...
             1  1  1  1; ...
             1  1  1  1 ];

if ~g.is_multiVM
    g.is_multiVM = 1;
    g = convert2multi(g);
end

g.is_LImodel = 1;
if g.is_LImodel
    g = initLateralInhibitionModel(g);
    y_bigger = sqrt(g.A0); 
    y_smaller = -sqrt(g.A0); 
    g = abolishLIregion(g, y_bigger, y_smaller);
%     g = setTopPCs(g, y_bigger, 1);
end
g.linkedCells = zeros(length(g.cells)-1,1);

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
    sim_folder = 'simulations\hybrid2\';
    vid = VideoWriter([sim_folder, 'sim_hybrid2_', num2str(latnum),'.avi']);
%     vid.Quality = 100; % 100 for best
    vid.FrameRate = 12;
    open(vid);
else
    vid = [];
end

% weights: incompessibility, tension, contractility, HC-HC repulsion,
% HC pull to PC
g.paras = [1 ; 0.05; 0.1; 1; 1];

%% run simulation
g = updateParameters1(g);
ts = 100; t=1;
while t<= ts
    savename = ['lattices\hybrid2\lat_', num2str(latnum),'_t=', num2str(t-1),'.mat'];
    save(savename, 'g');
    
    g = relaxLattice(g,3,vid);
    g = findTransitions(g,T1eps,T1prob,0.02);
    g = updateParameters1(g);
    g = relaxLattice(g,3, vid);
    g = killSmallCells(g,0.1);
    g = killLosingHCs(g);
    g = updateParameters1(g);
    g = relaxLattice(g,3, vid);
    g = handleVirtualVertices(g);
    g = findHoppingIntercalations(g, 0.9);
    
    disp(t)
    t = t+1;
    
end
close(vid);

savename = ['lattices\hybrid2\lat_', num2str(latnum),'_final.mat'];
save(savename, 'g');

end