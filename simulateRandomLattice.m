function simulateRandomLattice(nrow, ncol, lat_num, with_video)
% initially creates a hexagonal lattice of size nrowXncol and then
% shuffles the cells by randomizing the mechanical parameters every few
% timesteps. 
% lat_id = the ID of the lattice
% set with_video=1 if you want a movie of the simulation
%% default code
addpath(genpath(fileparts(which(mfilename))));
warning off

if nargin <= 3, with_video = 0; end
if nargin <= 2, lat_num = 1; end
if nargin <= 1, ncol = 12; end
if nargin == 0, nrow = 12; end
%% create hexagonal lattice and set general lattice parameters
[cc,ve] = mHexLattice(nrow+2,ncol+2,0);
ve(~isfinite(ve))=nan;

% remove unused vertices
vidx = zeros(length(ve),1);
for i = 1:length(cc)
    vidx(cc{i}) = 1;
end
vindex = find(vidx == 1);
vertices = ve(vindex,1:2);
vshift = zeros(length(ve),1);
vshift(vindex) = 1:length(vindex);
for i = 1:length(cc)
    cc{i} = vshift(cc{i})';
    cc{i}(end+1) = cc{i}(1);
end

disp(['Number of cells: ', num2str(length(cc))]);
clear g
g(1) = GLattConversion(cc,vertices);
nc = length(g.cells)-1;
g.xboundary=zeros(nc,2);
g.yboundary=zeros(nc,2);
g.bc=1;
g.scale =eye(2);
g.dead =zeros(nc,1);
g.is_multiVM = 0;
g.latticeDims = 2*pi*[1, nrow/ncol];
if(g.bc==1)
    g = periodicBC(g,nrow,ncol);
end

g = rescale(g);
g.dead =zeros(nc,1);
hex_area = cellarea(g,1);
g.A0 = average_area(g, 1:nc);
g.areas = ones(nc,1)*hex_area;

g.type = 2*ones(nc,1);
g.cmpr = 0.8 + 0.4*rand(nc, 1); %ones(nc, 1); % compressibility of cells
g.ctrc = 0.8 + 0.4*rand(nc, 1); %ones(nc, 1); % contractility of cells
g.bonds(:,5) = 1;
g.Gamma = [1 1 2 1 1]; % contractility of each cell type
g.alpha = [1 1 2 1 1]; % compressibility of each cell type
g.kappa = 12;


g.is_multiVM = 1;
if g.is_multiVM
    g = convert2multi(g);
end
g.is_LImodel = 0;

g.linkedCells = zeros(length(g.cells)-1,1);

% simulation parameters and global parameters
timer = 0;
g.with_noise = 0;
noise = rand(1, 2*length(g.verts)) - 0.5;
g.globs = struct('timer',timer,'noise',noise);

T1prob = 1;
T1eps = 0.2*sqrt(g.A0);
g.transitionedBonds = [0, 0];
g.relaxedCells = [0, 0];

% video
if with_video
    sim_folder = 'simulations/random';
    vid = VideoWriter([sim_folder, '/sim_rand_lat', num2str(lat_num),'.avi']);
    % vid.Quality = 100; % 100 for best
    vid.FrameRate = 12;
    open(vid);
else
    vid = [];
end

% weights: compressibility, tension, contractility, repulsion, lateral
g.paras = [1 ; 0.05; 0.1; 0; 0];

%% run simulation
delam_cells = randperm(nc, round(0.2*nc, 0));
ts = 100; t=1;
while t<= ts
 
    if mod(t-1,10)==0
        g.bonds(:,5) = 4*rand(length(g.bonds),1);
        for b=1:length(g.bonds)
            if g.bonds(b,1) ~=0
                invb = find(g.bonds(:,1) == g.bonds(b,2) & g.bonds(:,2) == g.bonds(b,1));
                g.bonds(invb, 5) = g.bonds(b,5);
            end
        end
        nc = length(g.cells)-1;
        g.areas = (0.5+2*rand(nc,1))*average_area(g, 1:nc);
        g.cmpr = 0.8 + 0.4*rand(nc, 1);
        g.ctrc = 0.8 + 0.4*rand(nc, 1);
    end
    if t==ts-9
        g.bonds(:,5) = 1;
        g.cmpr(:) = 1;
        g.ctrc(:) = 1;
        g.areas(:) = average_area(g, 1:nc);
    end
    g.cmpr(delam_cells) = 0;
    g.ctrc(delam_cells) = 5;
    
    
    g = relaxLattice(g,3,vid, 0.05);
    g = findTransitions(g,T1eps,T1prob,0.02);
    g = relaxLattice(g,3, vid, 0.05);
    g = killSmallCells(g,0.1);
    g = relaxLattice(g,3, vid, 0.05);
    g = handleVirtualVertices(g);

    disp(t)
    t = t+1;
    
end
close(vid);
latname = ['lattices\random\rand_lat_',num2str(nrow),'x',num2str(ncol),'_', num2str(lat_num),'.mat'];
save(latname, 'g');
end