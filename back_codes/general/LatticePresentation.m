function LatticePresentation(g,num, fig_num)
%% Draws the lattice g. If num is true, display the cell numbers
if nargin==1,
    num=0;
end

if nargin==3
    clf(fig_num);
end

cells = [1:length(g.cells)-1];
if(isfield(g,'LImodel'))
    isLI = 1;
    Ns = g.LImodel.cell_notch;
    Ds = g.LImodel.cell_delta;
    Js = g.LImodel.cell_jag;
    Nmax = 10; Nmin = 0;
    Dmax = 4; Dmin = 0;
    Jmax = 4; Jmin = 0;
    Ns = (Ns-Nmin)/(Nmax-Nmin);
    Ds = (Ds-Dmin)/(Dmax-Dmin);
    Js = (Js-Jmin)/(Jmax-Jmin);
    green = max(min(2*(Ns-Ds), 1),0);
    red = max(min(2*(Ds-Ns), 1),0);
    blue = max(min(Js, 1), 0);
else
    isLI = 0;
end
edge_color = 'k';
for i = cells
    if(length(g.cells{i+1})>2)
        if(g.dead(i)==0)
            verts = getVerts(g, i);
            if ismember(0,verts), disp(i); end
            if(length(verts)>2)
                v = getRelativePosition(g,verts);
                if(isfield(g,'scale'))
                    v = v*g.scale;
                end
                if(isfield(g,'type'))
                    if (isLI)
                        if ~g.LImodel.abolished(i)
                            col = ([red(i), green(i), 0]+1)/2;
                        else
                            switch g.type(i)
                                case 10
                                    col = [0.9 0.9 0.9];
                                case 20
                                    col = [0.8 0.8 0.8];
                                otherwise
                                    col = [1, 1, 1];
                            end
                            if(isfield(g,'top_boundary_cells'))
                                if ismember(i,g.top_boundary_cells)
                                    col = [0.5 0.72 0.87];
                                end
                            end
                        end
                        patch(v(:,1),v(:,2),col,'EdgeColor','k', 'LineWidth',1.5);
                    else
                        switch g.type(i)
                            case 0
                                face_col = [0.95 0.95 0.95];
                            case 1
                                face_col = [0.9 0.9 0.9];
                            case 2
                                face_col = [0.8 0.8 0.8];
                            case 3
                                face_col = [0 0.75 0];
                            case 4
                                face_col = [0.5 0.5 0.5];
                            case 5
                                face_col = [0.7 0.7 0.7];
                        end
%                         if g.linkedCells(i) ~= 0, face_col = [0.7 0.7 0.7]; end
                        patch('XData',v(:,1), 'YData', v(:,2),'FaceColor',face_col,'EdgeColor',edge_color,'LineWidth',1.5);
                    end
                else
                    patch(v(:,1),v(:,2),'w','EdgeColor','k');
                end
                hold on;
                if num == 1
                    text(mean(v(:,1)),mean(v(:,2)),num2str(i),'HorizontalAlignment','center','color','k');
                end
            end
        end
    end
end

axis equal
xl = 5; yl = 3.5; 
xlim([-xl, xl])
ylim([-yl, yl])
set(gca,'visible','off');
imsize = 800;
set(gcf,'Position',[50 300 imsize imsize*yl/xl]);
set(gca,'position',[.01 .01 .98 .98])
hold off;
end
