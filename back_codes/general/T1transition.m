function g1 = T1transition(g,b,pert)
timer = g.globs.timer;

if (nargin<4)
    pert = 0;
end
g1=g;
binv = find(g.bonds(:,1) == g.bonds(b,2) & g.bonds(:,2) == g.bonds(b,1));
if isempty(binv) || g.bonds(b,4)==0
    g1 = border_T1transition(g, b);
    return;
end

% check if the minimal time for transition has passed
T1delay = 40; % arbitraty period
bt = find(g1.transitionedBonds(:,1) == b);
if ~isempty(bt) && max(g1.transitionedBonds(bt, 2)) + T1delay > timer
    return;
end

edges = [b binv]; 
c1 = g.bonds(b,3); % the cell that contains the bond b
c3 = g.bonds(b,4); % the cell that contains the bond binv

% finds the bond that is counterclock-wise to b
f1 = find(g.cells{c1+1}==edges(1)); % the index of b in the cell's bonds array
if(f1==1)
    b1 = g.cells{c1+1}(length(g.cells{c1+1}));
else
    b1 = g.cells{c1+1}(f1-1);
end

% finds the bond that is counterclock-wise to binv
f2 = find(g.cells{c3+1}==edges(2)); % the index of binv in the cell's bonds array
if(f2==1)
    b2 = g.cells{c3+1}(length(g.cells{c3+1}));
else
    b2 = g.cells{c3+1}(f2-1);
end

c2 = g.bonds(b1,4); % the cell that is counterclock-wise to c1
c4 = g.bonds(b2,4); % the cell that is counterclock-wise to c3

% disp(strcat('--c1: ',num2str(c1),'--c2: ',num2str(c2),'--c3: ',num2str(c3),'c4: ',num2str(c4),'--'))
if(c2+c4==0 || c1==c4)
    return;
end
if(c1*c3==0 && c2*c4==0)
    return;
end

%% removing edge from cell c1 and c3
g1.cells{c1+1}=g.cells{c1+1}(~ismember(g.cells{c1+1},edges(1)));
g1.cells{c3+1}=g.cells{c3+1}(~ismember(g.cells{c3+1},edges(2)));

%% add edge to cell c2 and c4 between edges bidx
if c2~=0
    b1v = find(g.bonds(:,1)== g.bonds(b1,2) & g.bonds(:,2)== g.bonds(b1,1));
    bic = find(g.cells{c2+1} == b1v);
    g1.cells{c2+1} = [g.cells{c2+1}(1:bic-1) edges(1) g.cells{c2+1}(bic:length(g.cells{c2+1}))];
end
if c4~=0
    b2v = find(g.bonds(:,1)== g.bonds(b2,2) & g.bonds(:,2)== g.bonds(b2,1));
    bic = find(g.cells{c4+1} == b2v);
    g1.cells{c4+1} = [g.cells{c4+1}(1:bic-1) edges(2) g.cells{c4+1}(bic:length(g.cells{c4+1}))];
end
%% move edge extremities
g1.bonds(b1,2) = g.bonds(edges(1),2);
g1.bonds(b2,2) = g.bonds(edges(1),1);
if c2~=0, g1.bonds(b1v,1) = g.bonds(edges(2),1); end
if c4~=0, g1.bonds(b2v,1) = g.bonds(edges(2),2); end

%% update edge neighbors
if c2~=0
    g1.bonds(edges(1),3) = c2;
    g1.bonds(edges(1),4) = c4;
end
if c4~=0
    g1.bonds(edges(2),3) = c4;
    g1.bonds(edges(2),4) = c2;
end

%% remove residual virtual vertices
if g1.is_multiVM
    for i=1:2
        g1.bverts{edges(i)} = [];
    end
end

%% update vertex position
vert = getRelativePosition(g,g.bonds(b,1:2));
[blen, mid] = getBoundaryLength(g,b);
pvec1 =vert(2,:)-vert(1,:);
pvec1 = [-pvec1(2) pvec1(1)];
pvec2= -pvec1;


no1 = norm(pvec1);
no2 =norm(pvec2);
if(no1>0)
    pvec1 = pvec1/no1;
end
if(no2>0)
    pvec2 = pvec2/no2;
end
g1.verts(g.bonds(b,1),1:2)= mid+pert*blen*pvec1;
g1.verts(g.bonds(b,2),1:2)= mid+pert*blen*pvec2;

% if(g.bc == 1)
if(g.bc == 11)
    for i = 1:2
        if(g.xboundary(c2,i) && g.xboundary(c1,3-i) && g.xboundary(c3,3-i)),
            g1.xboundary(c4,3-i)=1;
        end
        if(g.yboundary(c2,i) && g.yboundary(c1,3-i) && g.yboundary(c3,3-i)),
            g1.yboundary(c4,3-i)=1;
        end
        
        if(g.xboundary(c4,i) && g.xboundary(c1,3-i) && g.xboundary(c3,3-i)),
            g1.xboundary(c2,3-i)=1;
        end
        if(g.yboundary(c4,i) && g.yboundary(c1,3-i) && g.yboundary(c3,3-i)),
            g1.yboundary(c2,3-i)=1;
        end
        if(g.xboundary(c1,i) && g.xboundary(c2,i)&& g.xboundary(c4,i) && g.xboundary(c3,3-i)),
            g1.xboundary(c1,i)=0;
        end
        if(g.yboundary(c1,i) && g.yboundary(c2,i)&& g.yboundary(c4,i) && g.yboundary(c3,3-i)),
            g1.yboundary(c1,i)=0;
        end
        if(g.xboundary(c3,i) && g.xboundary(c2,i)&& g.xboundary(c4,i) && g.xboundary(c1,3-i)),
            g1.xboundary(c3,i)=0;
        end
        if(g.yboundary(c3,i) && g.yboundary(c2,i)&& g.yboundary(c4,i) && g.yboundary(c1,3-i)),
            g1.yboundary(c3,i)=0;
        end
    end
    
    
end
% adding the bond to the delayed bonds list
g1.transitionedBonds = [g1.transitionedBonds; b timer; binv timer];

if c2==0, g1=remove_bond(g1, edges(1)); end
if c4==0, g1=remove_bond(g1, edges(2)); end
end